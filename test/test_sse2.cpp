
#include <emmintrin.h>

#include <cstdint>
#include <cstdio>

int main() {

  double d = 0.00000000000000000000000000000000000000000000034355123;
  float f = (float) d;
  printf("%g --> %g\n", d, f);

  exit(0);

  uint64_t local1[] = {1,46};
  
  __m128i vec1 = _mm_loadu_si128((__m128i*)local1);
  __m128i vec2 = _mm_set_epi64x((unsigned long)46, (unsigned long)1);
  
  _mm_storeu_si128((__m128i*)local1, vec2);
  
  printf("\n");
  for(int k = 0; k < 2; k++)
    printf("k: %lu\n", local1[k]);
  printf("\n");


  printf("vec1\n");
  for(int k = 0; k < 2; k++)
    printf("k: %lu\n", ((uint64_t*)&vec1)[k]);
  printf("\n");

  printf("vec2\n");
  for(int k = 0; k < 2; k++)
    printf("k: %lu\n", ((uint64_t*)&vec2)[k]);
  printf("\n");


  printf("done.\n");
}