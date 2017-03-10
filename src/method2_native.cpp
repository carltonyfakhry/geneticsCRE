#include <cstdio>
#include <cstring>
#include "gcre.h"

#define ALIGN __attribute__ (( aligned(16) ))

// 16 byte aligned
static int vector_width_ul(int count) {
  return 2 * (int) ceil(count / 128.0);
}

static inline void set_bit(uint64_t* vec, int idx){
  vec[idx/64] |= ONE_UL << idx % 64;      
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

joined_res JoinMethod2Native::join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u16& permute_cases,
  paths_type* p_paths0, paths_type* p_paths1, paths_type* p_pathsr, uint64_t total_paths) const {

  int vlen = vector_width_ul(conf.num_cases + conf.num_controls);
  int iters = 8 * ceil(max(1, conf.iterations) / 8.0);

  printf("adjusting width: %d -> %d\n", conf.num_cases + conf.num_controls, vlen * 64);
  printf("adjusting iterations: %d -> %d\n\n", conf.iterations, iters);

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

  const int vlen_b = vlen * sizeof(uint64_t);
  uint64_t ZERO_L_VLEN[vlen];
  for(int k = 0; k < vlen; k++)
    ZERO_L_VLEN[k] = ZERO_UL;

  // TODO temp
  // no clean way to init dynamic allocation
  uint64_t case_mask[vlen];
  for(int k = 0; k < vlen; k++)
    case_mask[k] = ZERO_UL;

  for(int k = 0; k < conf.num_cases; k++)
    case_mask[k/64] |= ONE_UL << k % 64;

  // TODO too big for the stack?
  uint64_t permute_mask[iters * vlen] ALIGN;
  for(int m = 0; m < vlen * iters; m++)
    permute_mask[m] = ZERO_UL;

  // ignoring input from R for now and just generating random permutations
  printf("generating simplified case masks for permutation (it = %d)\n", conf.iterations);
  srand((unsigned) time(0));

  permute_cases.clear();
  for(int k = 0; k < conf.iterations; k++){
    set<int> cases;
    permute_cases.push_back(vector<uint16_t>());
    while(cases.size() < conf.num_cases)
      cases.insert(rand() % (conf.num_cases + conf.num_controls));
    for(auto c : cases)
      permute_cases.back().push_back(c);
  }
  printf("permuted_cases: %lu (%lu)\n", permute_cases.size(), permute_cases.back().size());

  for(int k = 0; k < permute_cases.size(); k++){
    for(auto i : permute_cases[k])
      set_bit(permute_mask + k * vlen, (int) i);
  }

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

  double perm_score[conf.iterations];
  double perm_flips[conf.iterations];
  for(int k = 0; k < conf.iterations; k++){
    perm_score[k] = 0;
    perm_flips[k] = 0;
  }

  uint64_t joined_pos[vlen];
  uint64_t joined_neg[vlen];
  uint64_t joined_tp[vlen];
  uint64_t joined_tn[vlen];

  int p_case_pos[conf.iterations];
  int p_ctrl_pos[conf.iterations];

  for(int i = 0; i < uids.size(); i++){

    uid_ref& uid = uids[i];
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

      int case_pos = 0;
      int case_neg = 0;
      int ctrl_pos = 0;
      int ctrl_neg = 0;
      int total_pos = 0;
      int total_neg = 0;

      uint64_t bit_pos, bit_neg, bit_con, true_pos, true_neg, mask;

      // zero out join path holder
      memcpy(joined_pos, ZERO_L_VLEN, vlen_b);
      memcpy(joined_neg, ZERO_L_VLEN, vlen_b);

      for(int r = 0; r < conf.iterations; r++){
        p_case_pos[r] = 0;
        p_ctrl_pos[r] = 0;
      }

      for(int k = 0; k < vlen; k++){

        bit_pos = path_pos0[k] | path_pos1[k];
        bit_neg = path_neg0[k] | path_neg1[k];

        if(bit_pos != 0 || bit_neg != 0){

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

          if(true_pos != 0){
            for(int r = 0; r < conf.iterations; r++){
              uint64_t* p_mask = permute_mask + r * vlen;
              p_case_pos[r] += __builtin_popcountl(true_pos & p_mask[k]);
            }
          }

          if(true_neg != 0)
            for(int r = 0; r < conf.iterations; r++){
              uint64_t* p_mask = permute_mask + r * vlen;
              p_ctrl_pos[r] += __builtin_popcountl(true_neg & p_mask[k]);
            }
          }

          joined_pos[k] = bit_pos;
          joined_neg[k] = bit_neg;

        }

        if(keep_paths){
          memcpy(pathsr->pos + path_idx * vlen, joined_pos, vlen * sizeof(uint64_t));
          memcpy(pathsr->neg + path_idx * vlen, joined_neg, vlen * sizeof(uint64_t));
        }

        path_idx += 1;

        int cases = case_pos + case_neg;
        int ctrls = ctrl_pos + ctrl_neg;

        double score = value_table[cases][ctrls];
        double flips = value_table[ctrls][cases];

        if(score > scores.top().score)
          scores.push(Score(score, i, j));
        if(flips > scores.top().score)
          scores.push(Score(flips, i, j + flipped_pivot_length));
        while(scores.size() > conf.top_k)
          scores.pop();

        for(int r = 0; r < conf.iterations; r++){
          int p_cases = p_case_pos[r] + (total_neg - p_ctrl_pos[r]);
          int p_ctrls = p_ctrl_pos[r] + (total_pos - p_case_pos[r]);

          double p_score = value_table[p_cases][p_ctrls];
          if(p_score > perm_score[r])
            perm_score[r] = p_score;
          double p_flips = value_table[p_ctrls][p_cases];
          if(p_flips > perm_flips[r])
            perm_flips[r] = p_flips;
        }

      }

    }

    joined_res res;
    res.permuted_scores.resize(conf.iterations, 0);
    for(int k = 0; k < conf.iterations; k++)
      res.permuted_scores[k] = max(perm_score[k], perm_flips[k]);
    res.scores.clear();
    while(!scores.empty()){
      res.scores.push_back(scores.top());
      scores.pop();
    }

    return res;
  }
