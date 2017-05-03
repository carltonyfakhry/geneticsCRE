
size_t JoinMethod2::workspace_size_b() {
  size_t size = 0;
  // joined path
  size += exec->width_ul * 2 * sizeof(uint64_t);
  // per-iteration pos/neg case-counts
  size += exec->iterations * 2 * sizeof(uint32_t);
  return size;
}

void JoinMethod1::score_permute_cpu(int idx, int loc, const uint64_t* path0, const uint64_t* path1, uint64_t* path_res, bool keep_paths) {

  // avoid de-refs like the plague
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

  keep_score(idx, loc, cases, ctrls);

  int total = cases + ctrls;
  auto perm_scores = this->perm_scores;
  for(int r = 0; r < iters; r++){
    auto p_cases = perm_cases[r];
    auto p_score = value_table_max[p_cases][total - p_cases];
    if(p_score > perm_scores[r])
      perm_scores[r] = p_score;
  }

}

const uint64_t* JoinMethod2::score_permute_cpu(int idx, int loc, const uint64_t* path_pos0, const uint64_t* path_neg0, const uint64_t* path_pos1, const uint64_t* path_neg1, const size_t work_size, void* work) {

  // avoid de-refs like the plague
  const int iters = exec->iterations;
  const int width_ul = exec->width_ul;
  const auto case_mask = exec->case_mask;
  const auto perm_case_mask = exec->perm_case_mask;

  // clear workspace
  memset(work, 0, work_size);

  auto joined_full = (uint64_t*) work;
  auto joined_pos = joined_full;
  auto joined_neg = joined_full + width_ul;
  auto perm_cases = (uint32_t*) (joined_full + width_ul * 2);
  auto perm_case_pos = perm_cases;
  auto perm_ctrl_pos = perm_cases + iters;

  int case_pos = 0;
  int case_neg = 0;
  int ctrl_pos = 0;
  int ctrl_neg = 0;
  int total_pos = 0;
  int total_neg = 0;

  uint64_t bit_pos, bit_neg, bit_con, true_pos, true_neg, mask;

  // zero out join path holder
  // memset(joined_pos, 0, width_ul * sizeof(uint64_t));
  // memset(joined_neg, 0, width_ul * sizeof(uint64_t));

  // for(int r = 0; r < iters; r++){
  //   perm_case_pos[r] = 0;
  //   perm_ctrl_pos[r] = 0;
  // }

  for(int k = 0; k < width_ul; k++) {

    bit_pos = path_pos0[k] | path_pos1[k];
    bit_neg = path_neg0[k] | path_neg1[k];

    if(bit_pos != 0 || bit_neg != 0) {

      bit_con = bit_pos & bit_neg;
      true_pos = bit_pos & ~bit_con;
      true_neg = bit_neg & ~bit_con;
      mask = case_mask[k];

      total_pos += __builtin_popcountl(true_pos);
      total_neg += __builtin_popcountl(true_neg);
      case_pos += __builtin_popcountl(true_pos &  mask);
      case_neg += __builtin_popcountl(true_neg & ~mask);
      ctrl_pos += __builtin_popcountl(true_neg &  mask);
      ctrl_neg += __builtin_popcountl(true_pos & ~mask);

      if(true_pos != 0) {
        const uint64_t* p_mask = perm_case_mask + k * iters;
        for(int r = 0; r < iters; r++) {
          perm_case_pos[r] += __builtin_popcountl(true_pos & p_mask[r]);
        }
      }

      if(true_neg != 0) {
        const uint64_t* p_mask = perm_case_mask + k * iters;
        for(int r = 0; r < iters; r++){
          perm_ctrl_pos[r] += __builtin_popcountl(true_neg & p_mask[r]);
        }
      }

      joined_pos[k] = bit_pos;
      joined_neg[k] = bit_neg;

    }

  }

  int cases = case_pos + case_neg;
  int ctrls = ctrl_pos + ctrl_neg;

  keep_score(idx, loc, cases, ctrls);

  for(int r = 0; r < iters; r++){
    int perm_cases = perm_case_pos[r] + (total_neg - perm_ctrl_pos[r]);
    int perm_ctrls = perm_ctrl_pos[r] + (total_pos - perm_case_pos[r]);

    double p_score = value_table_max[perm_cases][perm_ctrls];
    if(p_score > perm_scores[r])
      perm_scores[r] = p_score;
  }

  // pointer to work segment we got from the caller
  return joined_full;
}
