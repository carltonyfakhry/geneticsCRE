
#warning "using sse2"

#include <emmintrin.h>

static inline void zero_vector_dq(void* p_vec, int size_dq) {
  auto vec = (__m128i*) p_vec;
  auto end = vec + size_dq;
  while(vec < end)
    _mm_store_si128(vec++, _mm_setzero_si128());
}

size_t JoinMethod1::workspace_size_b() {
  size_t size = 0;
  // joined path
  size += exec->width_ul * sizeof(uint64_t);
  // joined path holder for permutations
  size += exec->width_ul * sizeof(uint64_t);
  // per-iteration case-counts
  size += exec->iterations * sizeof(uint32_t);
  return size;
}

// local vec_joined may have been faster
const uint64_t* JoinMethod1::score_permute_sse2(int idx, int loc, const uint64_t* path0, const uint64_t* path1, const size_t work_size, void* work) {

  // avoid de-refs like the plague
  const int iters = exec->iterations;
  const int width_qw = exec->width_ul;
  const auto case_mask = exec->case_mask;
  const auto perm_case_mask = exec->perm_case_mask;

  const int width_dq = (width_qw + 1) / 2;

  // clear workspace
  // memset(work, 0, work_size);
  zero_vector_dq(work, work_size / 16);

  auto joined_full = (uint64_t*) work;
  auto vec_joined = (__m128i*) (joined_full + width_qw);
  auto perm_cases = (uint32_t*) (vec_joined + width_dq);
  auto joined = (uint64_t*) vec_joined;

  int cases = 0;
  int ctrls = 0;

  int nkept = 0;
  uint32_t kept[width_qw];
  memset(kept, 0, width_qw * 4);

  for(int k = 0; k < width_qw; k++) {

    if((joined[nkept] = path0[k] | path1[k]) != 0) {

      joined_full[k] = joined[nkept];

      cases += __builtin_popcountl(joined[nkept] &  case_mask[k]);
      ctrls += __builtin_popcountl(joined[nkept] & ~case_mask[k]);

      kept[nkept++] = k;
    }

  }

  // auto vec_2qw =  _mm_setzero_si128();
  // auto t2qw = (uint64_t*) &vec_2qw;

  uint64_t t2qw[2] ALIGNED;
  auto vec_2qw = (__m128i*) t2qw;

  for(int k = 0; k < nkept; k += 2){
    const int k0 = kept[k+0];
    const int k1 = kept[k+1];
    const auto p_mask_k0 = perm_case_mask + k0 * iters;
    const auto p_mask_k1 = perm_case_mask + k1 * iters;
    for(int r = 0; r < iters; r++) {
      _mm_store_si128(vec_2qw, _mm_and_si128(_mm_set_epi64x(p_mask_k1[r], p_mask_k0[r]), vec_joined[k/2]));
      perm_cases[r] += __builtin_popcountl(t2qw[0]) + __builtin_popcountl(t2qw[1]);
    }
  }

  keep_score(idx, loc, cases, ctrls);

  // doesn't win much, really
  int total = cases + ctrls;
  const auto& vt_max = value_table_max;
  auto perm_scores = (__m128*) this->perm_scores;
  for(int r = 0; r < iters; r += 4) {

    auto pc0 = perm_cases[r+0];
    auto pc1 = perm_cases[r+1];
    auto pc2 = perm_cases[r+2];
    auto pc3 = perm_cases[r+3];
    
    auto packed = _mm_setr_ps(vt_max[pc0][total - pc0], vt_max[pc1][total - pc1], vt_max[pc2][total - pc2], vt_max[pc3][total - pc3]);
    _mm_store_ps((float*) perm_scores, _mm_max_ps(*perm_scores, packed));
    perm_scores++;    
  }

  return joined_full;
}

/*
void JoinMethod2::score_permute_sse2(int idx, int loc, const uint64_t* path_pos0, const uint64_t* path_neg0, const uint64_t* path_pos1, const uint64_t* path_neg1) {

  // avoid de-refs like the plague
  const int top_k = exec->top_k;
  const int iters = exec->iterations;
  const int width_ul = exec->width_ul;
  const int flip_pivot_len = this->flip_pivot_len;

  int case_pos = 0;
  int case_neg = 0;
  int ctrl_pos = 0;
  int ctrl_neg = 0;
  int total_pos = 0;
  int total_neg = 0;

  uint64_t bit_pos, bit_neg, bit_con, true_pos, true_neg, mask;

  // zero out join path holder
  memset(joined_pos, 0, width_ul * sizeof(uint64_t));
  memset(joined_neg, 0, width_ul * sizeof(uint64_t));

  for(int r = 0; r < iters; r++){
    perm_case_pos[r] = 0;
    perm_ctrl_pos[r] = 0;
  }

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
        for(int r = 0; r < iters; r++) {
          const uint64_t* p_mask = perm_case_mask + r * width_ul;
          perm_case_pos[r] += __builtin_popcountl(true_pos & p_mask[k]);
        }
      }

      if(true_neg != 0) {
        for(int r = 0; r < iters; r++){
          const uint64_t* p_mask = perm_case_mask + r * width_ul;
          perm_ctrl_pos[r] += __builtin_popcountl(true_neg & p_mask[k]);
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

}
*/