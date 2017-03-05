#include "gcre.h"

// 16 byte aligned
static int vector_width_ul(int count) {
  return 2 * (int) ceil(count / 128.0);
}

struct paths_vect : paths_base {
  vec2d_u64 pos;
  vec2d_u64 neg;
  void resize(int size, int width_ul, int num_cases);
};

struct paths_block : paths_base {
  struct path {
    uint64_t hash = 0;
    short lengh = 0;
    short* index = NULL;
    uint64_t* pos;
    uint64_t* neg;
    uint64_t* con;
  };
};

void paths_vect::resize(int size, int width_ul, int num_cases) {
  paths_vect::size = size;
  paths_vect::width_ul = width_ul;
  paths_vect::num_cases = num_cases;
  pos.clear();
  neg.clear();
  pos.resize(size, vec_u64(width_ul, 0));
  neg.resize(size, vec_u64(width_ul, 0));
  printf("  ** resized and zeroed matrix: %d x %d\n", size, width_ul);
}

// static vec2d_u64 create_case_mask(vec2d_i& cases, int num_cases, int num_controls){

//   int vlen = vector_width_ul(num_cases + num_controls);

//   vec2d_u64 mask(cases.size(), vec_u64(vlen, 0));
//   for(int i = 0; i < cases.size(); i++){
//     for(int j = 0; j < num_cases; j++){
//       int index = j / 64;
//       if(cases[i][j] == 1)
//         mask[i][index] |= ONE_UL << j % 64;
//     }
//     for(int j = num_cases; j < num_cases + num_controls; j++){
//       int index = vlenc + (j-num_cases) / 64;
//       if(cases[i][j] == 1)
//         mask[i][index] |= ONE_UL << (j-num_cases) % 64;
//     }
//   }

//   return mask;
// }

paths_type* JoinMethod2Native::createPathSet() const {
  paths_vect* paths = new paths_vect;
  paths->size = 0;
  paths->width_ul = 0;
  return paths;
}

paths_type* JoinMethod2Native::createPathSet(vec2d_i& data, int num_cases, int num_controls) const {

  int vlen = vector_width_ul(num_cases + num_controls);

  // allocate on heap
  paths_vect* paths = new paths_vect;
  paths->resize(data.size(), vlen, num_cases);

  for(int r = 0; r < data.size(); r++) {
    for(int c = 0; c < data[r].size(); c++) {
      if(data[r][c] != 0) {
        paths->pos[r][c/64] |= ONE_UL << c % 64;
      }
    }
  }

  return paths;
}

joined_res JoinMethod2Native::join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_i& cases, paths_type* p_paths0, paths_type* p_paths1, paths_type* p_pathsr, uint64_t total_paths) const {

  int vlen = vector_width_ul(conf.num_cases + conf.num_controls);

  // not safe, but don't see a better way with R
  paths_vect* paths0 = (paths_vect*) p_paths0;
  paths_vect* paths1 = (paths_vect*) p_paths1;
  paths_vect* pathsr = (paths_vect*) p_pathsr;

  // will be deallocated automatically
  paths_vect tpaths;
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

  int flipped_pivot_length = paths1->pos.size();
  int prev_src = -1;
  int total_srcs = 0;
  int path_idx = 0;

  for(int i = 0; i < uids.size(); i++){
    // Check if source node has changed
    uid_ref& uid = uids[i];
    if(uid.src != prev_src){
      prev_src = uid.src;
      total_srcs++;
    }

    if(uid.count == 0) continue;
    int sign;
    if(conf.path_length > 3) sign = join_gene_signs[i];

    vec_u64 &path_pos0 = paths0->pos[i];
    vec_u64 &path_neg0 = paths0->neg[i];

    for(int j = uid.location; j < (uid.location + uid.count); j++){

      if(conf.path_length < 3) sign = join_gene_signs[j];

      if(conf.path_length == 3){
        int sign2 = join_gene_signs[i];
        int sign3 = join_gene_signs[j];
        sign = ((sign2 == -1 && sign3 == 1) || (sign2 == 1 && sign3 == -1)) ? -1 : 1;
      }

      vec_u64 joined_pos(path_pos0.size());
      vec_u64 joined_neg(path_neg0.size());

      // TODO ??
      vec_u64 &path_pos1 = (sign == 1) ? paths1->pos[j] : paths1->neg[j];
      vec_u64 &path_neg1 = (sign == 1) ? paths1->neg[j] : paths1->pos[j];

      for(int k = 0; k < vlen; k++){
        joined_pos[k] = path_pos0[k] | path_pos1[k];
        joined_neg[k] = path_neg0[k] | path_neg1[k];
      }

      if(keep_paths){
        pathsr->pos[path_idx] = joined_pos;
        pathsr->neg[path_idx] = joined_neg;
      }

      path_idx += 1;

      vec_u64 joined_tp(path_pos0.size(), 0);
      vec_u64 joined_tn(path_neg0.size(), 0);

      for(int k = 0; k < vlen; k++) {
        joined_tp[k] = joined_pos[k] & ~(joined_pos[k] & joined_neg[k]);
        joined_tn[k] = joined_neg[k] & ~(joined_pos[k] & joined_neg[k]);
      }

      int cpos = 0;
      int cneg = 0;
      int tpos = 0;
      int tneg = 0;

      for(int k = 0; k < vlen; k++) {
        cpos += __builtin_popcountl(joined_tp[k] & case_mask[k]);
        cneg += __builtin_popcountl(joined_tn[k] & ~case_mask[k]);
        tpos += __builtin_popcountl(joined_tn[k] & case_mask[k]);
        tneg += __builtin_popcountl(joined_tp[k] & ~case_mask[k]);
      }
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
