#ifndef GCRE_METHODS_H
#define GCRE_METHODS_H

class JoinMethod_Base : public JoinMethod {

public:

  JoinMethod_Base(const JoinExec* exec, const UidRelSet& uids, float* p_perm_scores) :

  exec(exec),
  uids(uids),
  value_table(exec->value_table) {

    scores.push(Score());

    perm_scores = p_perm_scores;
    for(int k = 0; k < exec->iterations; k++)
      perm_scores[k] = 0;

  }

  virtual void score_permute(int idx, int loc, const uint64_t* path0, const uint64_t* path1, uint64_t* path_res, bool keep_paths) = 0;

  // run at end of worker thread to collect results
  virtual void merge_scores() {

    while(!scores.empty()) {
      auto score = scores.top();
      if(score.score > exec->scores.top().score)
        exec->scores.push(score);
      scores.pop();
    }

    for(int r = 0; r < exec->iterations; r++) {
      if(perm_scores[r] > exec->perm_scores[r])
        exec->perm_scores[r] = perm_scores[r];
    }

  }

protected:

  const JoinExec* exec;
  const UidRelSet& uids;
  const vec2d_d& value_table;

  priority_queue<Score> scores;
  float* perm_scores;

};

class JoinMethod1 : public JoinMethod_Base {

public:

  JoinMethod1(const JoinExec* exec, const UidRelSet& uids, float* p_perm_scores) : JoinMethod_Base(exec, uids, p_perm_scores) {}

  void score_permute(int idx, int loc, const uint64_t* path0, const uint64_t* path1, uint64_t* path_res, bool keep_paths) {

    const int top_k = exec->top_k;
    const int iters = exec->iterations;
    const int width_ul = exec->width_ul;
    const auto case_mask = exec->case_mask;
    const auto perm_case_mask = exec->perm_case_mask;

    uint64_t joined = 0;
    int cases = 0;
    int ctrls = 0;

    uint32_t perm_cases[iters];
    memset(perm_cases, 0, iters * 4);

    for(int k = 0; k < width_ul; k++) {

      if((joined = path0[k] | path1[k]) != 0) {

        cases += __builtin_popcountl(joined &  case_mask[k]);
        ctrls += __builtin_popcountl(joined & ~case_mask[k]);

        const uint64_t* p_mask = perm_case_mask + k * iters;
        for(int r = 0; r < iters; r++)
          perm_cases[r] += __builtin_popcountl(joined & p_mask[r]);

        if(keep_paths)
          path_res[k] = joined;
      }

    }

    double score = value_table[cases][ctrls];
    if(score > scores.top().score)
      scores.push(Score(score, idx, loc, cases, ctrls));
    while(scores.size() > top_k)
      scores.pop();

    int total = cases + ctrls;
    auto perm_scores = this->perm_scores;
    for(int r = 0; r < iters; r++){
      auto p_cases = perm_cases[r];
      auto p_score = value_table[p_cases][total - p_cases];
      if(p_score > perm_scores[r])
        perm_scores[r] = p_score;
    }

  }

};

// precompute max, including flipped indices
static vec2d_d compute_value_table_max(const vec2d_d& value_table) {
  auto value_table_max = value_table;
  for(int r = 0; r < value_table.size(); r++) {
    for(int c = 0; c < value_table[r].size(); c++) {
      value_table_max[r][c] = max(value_table[r][c], value_table[c][r]);
    }
  }
  return value_table_max;
}

class JoinMethod2 : public JoinMethod_Base {

public:

  // JoinMethod2(const JoinExec* exec, const UidRelSet& uids, const int flip_pivot_len, float* p_perm_scores) :
  JoinMethod2(const JoinExec* exec, const UidRelSet& uids, float* p_perm_scores) :
  JoinMethod_Base(exec, uids, p_perm_scores),
  // flip_pivot_len(flip_pivot_len),
  value_table_max(compute_value_table_max(exec->value_table)) {}

  void score_permute(int idx, int loc, const uint64_t* path0, const uint64_t* path1, uint64_t* path_res, bool keep_paths) {

    const int iters = exec->iterations;
    const int width_ul = exec->width_ul;
    const auto case_mask = exec->case_mask;
    const auto perm_case_mask = exec->perm_case_mask;

    const auto path_pos0 = path0;
    const auto path_neg0 = path0 + width_ul;

    const auto do_flip = uids.need_flip(idx, loc);
    const auto path_pos1 = (do_flip ? path1 : path1 + width_ul);
    const auto path_neg1 = (do_flip ? path1 + width_ul : path1);

    auto joined_pos = path_res;
    auto joined_neg = path_res + width_ul;

    uint32_t perm_count_pos[iters * 2] ALIGNED;
    memset(perm_count_pos, 0, iters * 8);
    auto perm_case_pos = perm_count_pos;
    auto perm_ctrl_pos = perm_count_pos + iters;

    uint32_t case_pos = 0;
    uint32_t case_neg = 0;
    uint32_t ctrl_pos = 0;
    uint32_t ctrl_neg = 0;
    uint32_t total_pos = 0;
    uint32_t total_neg = 0;

    // uint64_t bit_pos, bit_neg, bit_con, true_pos, true_neg, mask;
    uint64_t bit_pos, bit_neg, mask;

    for(int k = 0; k < width_ul; k++) {

      bit_pos = path_pos0[k] | path_pos1[k];
      bit_neg = path_neg0[k] | path_neg1[k];

      if(bit_pos != 0 || bit_neg != 0) {

        // bit_con = bit_pos & bit_neg;
        // true_pos = bit_pos & ~bit_con;
        // true_neg = bit_neg & ~bit_con;
        mask = case_mask[k];

        // total_pos += __builtin_popcountl(true_pos);
        // total_neg += __builtin_popcountl(true_neg);
        // case_pos += __builtin_popcountl(true_pos &  mask);
        // case_neg += __builtin_popcountl(true_neg & ~mask);
        // ctrl_pos += __builtin_popcountl(true_neg &  mask);
        // ctrl_neg += __builtin_popcountl(true_pos & ~mask);
        total_pos += __builtin_popcountl(bit_pos);
        total_neg += __builtin_popcountl(bit_neg);
        case_pos += __builtin_popcountl(bit_pos &  mask);
        case_neg += __builtin_popcountl(bit_neg & ~mask);
        ctrl_pos += __builtin_popcountl(bit_neg &  mask);
        ctrl_neg += __builtin_popcountl(bit_pos & ~mask);

        // if(true_pos != 0) {
        if(bit_pos != 0) {
          const uint64_t* p_mask = perm_case_mask + k * iters;
          for(int r = 0; r < iters; r++) {
            // perm_case_pos[r] += __builtin_popcountl(true_pos & p_mask[r]);
            perm_case_pos[r] += __builtin_popcountl(bit_pos & p_mask[r]);
          }
        }

        // if(true_neg != 0) {
        if(bit_neg != 0){
          const uint64_t* p_mask = perm_case_mask + k * iters;
          for(int r = 0; r < iters; r++){
            // perm_ctrl_pos[r] += __builtin_popcountl(true_neg & p_mask[r]);
            perm_ctrl_pos[r] += __builtin_popcountl(bit_neg & p_mask[r]);
          }
        }

        if(keep_paths) {
          joined_pos[k] = bit_pos;
          joined_neg[k] = bit_neg;
        }

      }

    }

    // int cases = case_pos + case_neg;
    // int ctrls = ctrl_pos + ctrl_neg;

    // keep_score(idx, loc, cases, ctrls);
    keep_score(idx, loc, case_pos, ctrl_neg, case_neg, ctrl_pos);

    for(int r = 0; r < iters; r++){
      // int perm_cases = perm_case_pos[r] + (total_neg - perm_ctrl_pos[r]);
      // int perm_ctrls = perm_ctrl_pos[r] + (total_pos - perm_case_pos[r]);
      int perm_case_neg = total_neg - perm_ctrl_pos[r];
      int perm_ctrl_neg = total_pos - perm_case_pos[r];

      // double p_score = value_table_max[perm_cases][perm_ctrls];
      double p_score = value_table_max[perm_case_pos[r]][perm_ctrl_neg] + value_table_max[perm_case_neg][perm_ctrl_pos[r]];
      if(p_score > perm_scores[r])
        perm_scores[r] = p_score;
    }

  }

protected:

  // const int flip_pivot_len;
  const vec2d_d value_table_max;

  void keep_score(int idx, int loc, int cases, int ctrls) {

    double score = value_table[cases][ctrls];
    if(score > scores.top().score)
      scores.push(Score(score, idx, loc, cases, ctrls));

    // double flips = value_table[ctrls][cases];
    // if(flips > scores.top().score)
    //   scores.push(Score(flips, idx, loc + flip_pivot_len, cases, ctrls));
    while(scores.size() > exec->top_k)
      scores.pop();

  }

  void keep_score(int idx, int loc, int case_pos, int ctrl_neg, int case_neg, int ctrl_pos) {

    double score = value_table[case_pos][ctrl_neg] + value_table[case_neg][ctrl_pos];
    int cases = case_pos + case_neg;
    int ctrls = ctrl_pos + ctrl_neg;
    if(score > scores.top().score)
      scores.push(Score(score, idx, loc, cases, ctrls));

    while(scores.size() > exec->top_k)
      scores.pop();

  }

};

#endif
