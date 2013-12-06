#ifdef __AVX__
#define NSSE 8 /* Number of floats in an AVX register */
#define Arch(a) a##_avx8
typedef float v8sf __attribute__ ((vector_size (32)));
typedef float vsf __attribute__ ((vector_size (32)));
#define vsf_hsum(a) (a[0]+a[1]+a[2]+a[3]+a[4]+a[5]+a[6]+a[7])
#define vsf_scalar(a) {a, a, a, a, a, a, a, a}
#define vsf_rsqrt(_r2) __builtin_ia32_rsqrtps256(_r2)
#else
#define NSSE 4 /* Number of floats in an SSE register */
#define Arch(a) a##_sse4
typedef float v4sf __attribute__ ((vector_size (16)));
typedef float vsf __attribute__ ((vector_size (16)));
#define vsf_hsum(a) (a[0]+a[1]+a[2]+a[3])
#define vsf_scalar(a) {a, a, a, a}
#define vsf_rsqrt(_r2) __builtin_ia32_rsqrtps(_r2)
#endif

