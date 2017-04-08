
// configure instruction set based on cpu flags

#ifndef GCRE_INSTR_H
#define GCRE_INSTR_H

#if defined __AVX512__

#define COMPILE_EVEX
#define SCORE_METHOD_NAME score_permute_evex
constexpr int gs_vec_width = 512;
const std::string gs_impl_label = "AVX-512";

#elif defined __AVX2__

#define COMPILE_AVX2
#define SCORE_METHOD_NAME score_permute_avx2
constexpr int gs_vec_width = 256;
const std::string gs_impl_label = "AVX2";

// SSE4 also requires AVX (don't think that affects many CPUs)
#elif defined __AVX__

#define COMPILE_SSE4
#define SCORE_METHOD_NAME score_permute_sse4
constexpr int gs_vec_width = 256;
const std::string gs_impl_label = "SSE4";

#else

#define COMPILE_SSE2
#define SCORE_METHOD_NAME score_permute_sse2
constexpr int gs_vec_width = 128;
const std::string gs_impl_label = "SSE2";

#endif

#endif