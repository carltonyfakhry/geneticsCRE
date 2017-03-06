#include <cstdio>
#include <cstring>
#include "gcre.h"

// 16 byte aligned
static int vector_width_ul(int count) {
  return 2 * (int) ceil(count / 128.0);
}

struct paths_block : paths_base {
  uint64_t* pos = NULL;
  uint64_t* neg = NULL;
  void resize(int size, int width_ul, int num_cases);
};

void paths_block::resize(int size, int width_ul, int num_cases) {
  paths_block::size = size;
  paths_block::width_ul = width_ul;
  paths_block::num_cases = num_cases;

  if(pos != NULL)
    delete[] pos;
  if(neg != NULL)
    delete[] neg;

  pos = new uint64_t[size * width_ul]();
  neg = new uint64_t[size * width_ul]();

  printf("  ** resized and zeroed array: %d x %d\n", size, width_ul);
}

paths_type* JoinMethod2Native::createPathSet() const {
  paths_block* paths = new paths_block;
  paths->size = 0;
  paths->width_ul = 0;
  paths->num_cases = 0;
  return paths;
}

paths_type* JoinMethod2Native::createPathSet(vec2d_i& data, int num_cases, int num_controls) const {

  int vlen = vector_width_ul(num_cases + num_controls);

  // allocate on heap
  paths_block* paths = new paths_block;
  paths->resize(data.size(), vlen, num_cases);

  for(int r = 0; r < data.size(); r++) {
    for(int c = 0; c < data[r].size(); c++) {
      if(data[r][c] != 0) {
        paths->pos[r * vlen + c/64] |= ONE_UL << c % 64;
      }
    }
  }

  return paths;
}

// create directly from bitset vectors
paths_type* JoinMethod2Native::createPathSet(vec2d_u64& pos, vec2d_u64& neg, int num_cases, int num_controls) const {

  int vlen = vector_width_ul(num_cases + num_controls);

  // allocate on heap
  paths_block* paths = new paths_block;
  paths->resize(pos.size(), pos.front().size(), num_cases);
  for(int r = 0; r < paths->size; r++){
    for(int k = 0; k < paths->width_ul; k++){
      paths->pos[r * vlen + k] = pos[r][k];
      paths->neg[r * vlen + k] = neg[r][k];
    }
  }
  return paths;
}

joined_res JoinMethod2Native::join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_i& cases,
  paths_type* p_paths0, paths_type* p_paths1, paths_type* p_pathsr, uint64_t total_paths) const {

  int vlen = vector_width_ul(conf.num_cases + conf.num_controls);
  // not safe, but don't see a better way with R
  paths_block* paths0 = (paths_block*) p_paths0;
  paths_block* paths1 = (paths_block*) p_paths1;
  paths_block* pathsr = (paths_block*) p_pathsr;

  // will be deallocated automatically
  paths_block tpaths;
  if(paths0 == NULL){
    tpaths.resize(uids.size(), vlen, conf.num_cases);
    paths0 = &tpaths;
  }

  bool keep_paths = pathsr != NULL;

  if(keep_paths)
    pathsr->resize(total_paths, paths1->width_ul, conf.num_cases);

  // vec2d_u64 case_mask = create_case_mask(cases, conf.num_cases, conf.num_controls);

  // TODO temp
  vec_u64 case_mask(vlen, 0);
  for(int k = 0; k < conf.num_cases; k++)
    case_mask[k/64] |= ONE_UL << k % 64;

  // priority queue for the indices of the top K paths in the data
  // add dummy score to avoid empty check
  priority_queue<Score> scores;
  scores.push(Score());

  // TODO conf.nthreads
  int nthreads = 1;

  int flipped_pivot_length = paths1->size;
  int prev_src = -1;
  int total_srcs = 0;
  int path_idx = 0;

  uint64_t joined_pos[vlen];
  uint64_t joined_neg[vlen];
  uint64_t joined_tp[vlen];
  uint64_t joined_tn[vlen];

  for(int i = 0; i < uids.size(); i++){

    // Check if source node has changed
    uid_ref& uid = uids[i];
    if(uid.src != prev_src){
      prev_src = uid.src;
      total_srcs++;
    }

    if(uid.count == 0)
      continue;

    uint64_t* path_pos0 = paths0->pos + i * vlen;
    uint64_t* path_neg0 = paths0->neg + i * vlen;

    for(int j = uid.location; j < (uid.location + uid.count); j++){

      int sign = 0;
      if(conf.path_length > 3)
        sign = join_gene_signs[i];
      else if(conf.path_length < 3)
        sign = join_gene_signs[j];
      else if(conf.path_length == 3)
        sign = (join_gene_signs[i] + join_gene_signs[j] == 0) ? -1 : 1;

      uint64_t* path_pos1 = (sign == 1 ? paths1->pos : paths1->neg) + j * vlen;
      uint64_t* path_neg1 = (sign == 1 ? paths1->neg : paths1->pos) + j * vlen;

      int cpos = 0;
      int cneg = 0;
      int tpos = 0;
      int tneg = 0;

      uint64_t p, n, tp, tn, c, m;

      for(int k = 0; k < vlen; k++){

        p = path_pos0[k] | path_pos1[k];
        n = path_neg0[k] | path_neg1[k];
        c = p & n;
        tp = p & ~c;
        tn = n & ~c;
        m = case_mask[k];

        cpos += __builtin_popcountl(tp & m);
        cneg += __builtin_popcountl(tn & ~m);
        tpos += __builtin_popcountl(tn & m);
        tneg += __builtin_popcountl(tp & ~m);

        joined_pos[k] = p;
        joined_neg[k] = n;

      }

      if(keep_paths){
        memcpy(pathsr->pos + path_idx * vlen, joined_pos, vlen * sizeof(uint64_t));
        memcpy(pathsr->neg + path_idx * vlen, joined_neg, vlen * sizeof(uint64_t));
      }

      path_idx += 1;

      int cases = cpos + cneg;
      int controls = tpos + tneg;

      double score = value_table[cases][controls];
      double flipped_score = value_table[controls][cases];

      if(score > scores.top().score){
        scores.push(Score(score, i, j));
        if(scores.size() > conf.top_k)
          scores.pop();
      }
      if(flipped_score > scores.top().score){
        scores.push(Score(flipped_score, i, j + flipped_pivot_length));
        if(scores.size() > conf.top_k)
          scores.pop();
      }

    //   for(int m = 0; m < conf.iterations; m++){

    //     double perm_score = 0;
    //     double perm_flipped_score = 0;

    //     vec_u64& caseorcontrol = case_mask[m];
    //     int cases_pos_m = 0;
    //     int controls_pos_m = 0;
    //     int cases_neg_m = 0;
    //     int controls_neg_m = 0;

    //     for(int k = 0; k < vlenc; k++){
    //       uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
    //       cases_pos_m += __builtin_popcountll(permuted_path_k);
    //       permuted_path_k = joined_neg[k] & caseorcontrol[k];
    //       controls_neg_m += __builtin_popcountll(permuted_path_k);
    //     }
    //     for(int k = vlenc; k < vlenc + vlent; k++){
    //       uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
    //       controls_pos_m += __builtin_popcountll(permuted_path_k);
    //       permuted_path_k = joined_neg[k] & caseorcontrol[k];
    //       cases_neg_m += __builtin_popcountll(permuted_path_k);
    //     }

    //     int new_cases = cases_pos_m + cases_neg_m + (controls_pos - controls_pos_m) + (controls_neg - controls_neg_m);
    //     int new_controls = controls_pos_m + controls_neg_m + (cases_pos - cases_pos_m) + (cases_neg - cases_neg_m);

    //     perm_score = value_table[new_cases][new_controls];
    //     perm_flipped_score = value_table[new_controls][new_cases];
    //   }
    }

  }

  joined_res res;
  res.permuted_scores.resize(conf.iterations, 0);
  res.scores.clear();
  while(!scores.empty()){
    res.scores.push_back(scores.top());
    scores.pop();
  }

  return res;
}
