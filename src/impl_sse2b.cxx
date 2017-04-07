
// sse2
#include <emmintrin.h>

void JoinMethod1::score_permute_sse2(int idx, int loc, const uint64_t* path0, const uint64_t* path1) {

  // avoid de-refs like the plague
  const int top_k = exec->top_k;
  const int iters = exec->iterations;
  const int width_ul = exec->width_ul;
  const int flip_pivot_len = this->flip_pivot_len;
  float* perm_scores = this->perm_scores;

  int cases = 0;
  int ctrls = 0;

  int width_dq = (width_ul + 1) / 2;

  // zero out join path holder
  memset(joined_block, 0, width_ul * sizeof(uint64_t));

  __m128i vec_joined[width_dq];
  memset(vec_joined, 0, width_dq * 16);
  uint64_t* joined = (uint64_t*) vec_joined;

  int nkept = 0;
  uint16_t kept[width_ul];

  for(int k = 0; k < width_ul; k++) {

    joined[nkept] = path0[k] | path1[k];

    if(joined[nkept] != 0) {

      joined_block[k] = joined[nkept];

      cases += __builtin_popcountl(joined[nkept] &  case_mask[k]);
      ctrls += __builtin_popcountl(joined[nkept] & ~case_mask[k]);

      kept[nkept++] = k;
    }


  }

  // uint64_t j = 0;
  uint64_t t2qw[2];

  int perm_counts[iters];
  memset(perm_counts, 0, iters * 4);
  auto vec_2qw = (__m128i*) t2qw;

  // for(int r = 0; r < iters; r++) {
  //   perm_counts[r] = 0;
  //   const uint64_t* p_mask = perm_case_mask + r * width_ul;
  //   for(int k = 0; k < nkept; k += 2){
  //     _mm_storeu_si128(vec_2qw, _mm_and_si128(_mm_set_epi64x(p_mask[kept[k+1]], p_mask[kept[k]]), vec_joined[k/2]));
  //     perm_counts[r] += __builtin_popcountl(t2qw[0]) + __builtin_popcountl(t2qw[1]);
  //   }
  // }

  // for(int p = 0; p < iters * width_ul; p += 16)
    // _mm_prefetch((char*) (perm_case_mask + p), _MM_HINT_T1);

// marginally faster with reordered mask elements
  for(int k = 0; k < nkept; k += 2){
    int k0 = kept[k];
    int k1 = kept[k+1];
    const uint64_t* p_mask_k0 = perm_case_mask + k0;
    const uint64_t* p_mask_k1 = perm_case_mask + k1;
    for(int r = 0; r < iters; r++) {
      // const int off = r * width_ul;
      _mm_storeu_si128(vec_2qw, _mm_and_si128(_mm_set_epi64x(p_mask_k1[r], p_mask_k0[r]), vec_joined[k/2]));
      // _mm_storeu_si128(vec_2qw, _mm_and_si128(_mm_set_epi64x(k, k+1), vec_joined[k/2]));
      perm_counts[r] += __builtin_popcountl(t2qw[0]) + __builtin_popcountl(t2qw[1]);
    }
  }

/* not faster than cpu
  for(int k = 0; k < nkept; k++){
    int kidx = kept[k];
    int k1 = kept[k+1];

    __m128i j = _mm_set1_epi64x(joined[k]);

    const __m128i* p_mask = (__m128i*) perm_case_mask + kept[k];
    // const uint64_t* p_mask_k0 = perm_case_mask + k0;
    // const uint64_t* p_mask_k1 = perm_case_mask + k1;
    for(int r = 0; r < iters / 2; r++) {
      // const int off = r * width_ul;
      _mm_storeu_si128(vec_2qw, _mm_and_si128(p_mask[r], j));
      // _mm_storeu_si128(vec_2qw, _mm_and_si128(_mm_set_epi64x(k, k+1), vec_joined[k/2]));
      perm_counts[r*2+0] += __builtin_popcountl(t2qw[0]);
      perm_counts[r*2+1] += __builtin_popcountl(t2qw[1]);
    }
  }
*/

  // printf("nkept: %d\n", nkept);

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
    int p_cases = perm_counts[r];
    float p_score = (float) value_table_max[p_cases][total - p_cases];
    // if(p_score > perm_scores[r])
    //   perm_scores[r] = p_score;
  }

  // float p_score = (float) value_table_max[p_cases][total - p_cases];
  // if(p_score > perm_scores[r])
  //   perm_scores[r] = p_score;

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
