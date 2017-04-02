
// TODO check flags sse2
// sse2 and sse4.1
#include <emmintrin.h>
// width: 128bit (double quad)

static inline void zero_vector(void* p_vec, int size) {
  auto vec = (__m128i*) p_vec;
  auto end = vec + size;
  while(vec < end)
    _mm_store_si128(vec++, _mm_setzero_si128());
}

void JoinMethod1::score_permute_sse2(int idx, int loc, const uint64_t* path0, const uint64_t* path1) {

  // avoid de-refs like the plague
  const int top_k = exec->top_k;
  const int iters = exec->iterations;
  // TODO should never need ul
  const int width_ul = exec->width_ul;
  const int width_dq = exec->width_dq;
  const int flip_pivot_len = this->flip_pivot_len;

  // TODO set consistently
  const int iters_dw = iters / 4;

  int cases = 0;
  int ctrls = 0;

  zero_vector(joined_block, width_dq);
  zero_vector(perm_count_block, iters_dw);

  auto vec_joined_block = (__m128i*) joined_block;

  auto joined = _mm_setzero_si128();

  const auto vec_path0 = (__m128i*) path0;
  const auto vec_path1 = (__m128i*) path1;
  const auto vec_case_mask = (__m128i*) case_mask;
  const auto vec_perm_case_mask = (__m128i*) perm_case_mask;

  for(int k = 0; k < width_dq; k++) {

    joined = _mm_or_si128(vec_path0[k], vec_path1[k]);

    if(!_mm_test_all_zeros(joined, joined)) {

      auto vec_case_mask_k = _mm_load_si128(vec_case_mask + k);

      auto vec_cases = _mm_and_si128(vec_case_mask_k, joined);
      cases += __builtin_popcountll(_mm_extract_epi64(vec_cases, 0));
      cases += __builtin_popcountll(_mm_extract_epi64(vec_cases, 1));

      auto vec_ctrls = _mm_andnot_si128(vec_case_mask_k, joined);
      ctrls += __builtin_popcountll(_mm_extract_epi64(vec_ctrls, 0));
      ctrls += __builtin_popcountll(_mm_extract_epi64(vec_ctrls, 1));

      const auto p_mask = perm_case_mask + k;
      for(int r = 0; r < iters; r++) {
        auto perm_cases = _mm_and_si128(_mm_load_si128(vec_perm_case_mask + r * width_dq), joined);
        perm_count_block[r] += __builtin_popcountll(_mm_extract_epi64(perm_cases, 0)) + __builtin_popcountll(_mm_extract_epi64(perm_cases, 1));
      }

      _mm_storeu_si128(vec_joined_block + k, joined);
    }

  }

  double score = value_table[cases][ctrls];
  double flips = value_table[ctrls][cases];

  if(score > scores.top().score)
    scores.push(Score(score, idx, loc));
  if(flips > scores.top().score)
    scores.push(Score(flips, idx, loc + flip_pivot_len));
  while(scores.size() > top_k)
    scores.pop();


  int total = cases + ctrls;
  for(int r = 0; r < iters; r++){
    int p_cases = perm_count_block[r];
    double p_score = value_table_max[p_cases][total - p_cases];
    if(p_score > perm_scores[r])
      perm_scores[r] = p_score;
  }

}

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

  double score = value_table[cases][ctrls];
  double flips = value_table[ctrls][cases];

  if(score > scores.top().score)
    scores.push(Score(score, idx, loc));
  if(flips > scores.top().score)
    scores.push(Score(flips, idx, loc + flip_pivot_len));
  while(scores.size() > top_k)
    scores.pop();

  for(int r = 0; r < iters; r++){
    int perm_cases = perm_case_pos[r] + (total_neg - perm_ctrl_pos[r]);
    int perm_ctrls = perm_ctrl_pos[r] + (total_pos - perm_case_pos[r]);

    double p_score = value_table_max[perm_cases][perm_ctrls];
    if(p_score > perm_scores[r])
      perm_scores[r] = p_score;
  }

}
