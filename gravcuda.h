#ifdef __CUDACC__
#define EXTERN extern "C"
#define CUDA
#else
#define EXTERN
#endif

#ifdef CUDA
EXTERN void CUDA_Init(void *device_prop);
EXTERN void CUDA_Finalize(void);
EXTERN int vecwidthCUDA(void);
EXTERN int qallocCUDA(int *inuse);
EXTERN void WalkInitSrcCUDA(float *qtab, int stride, int64_t nobj);
EXTERN void WalkTerminateSrcCUDA(void);
EXTERN void WalkInitSinkCUDA(float *btab, int stride, int64_t nobj);
EXTERN void WalkTerminateSinkCUDA(float *btab, int stride, int64_t nobj);
EXTERN void grav_mn_CUDA(const char *routine, const float *p, float *accp, const int n, const int stride, 
			 const float *f, const int source_n, const int sz, 
			 const float e, int *ncut, int q);
EXTERN void grav_mns_CUDA(const char *routine, const int base, const int m,
			  const segment *seg, const int seg_n, const int source_n, 
			  const float mmass, float e, int *ncut, int q);
EXTERN void grav_mnss_CUDA(const char *routine, const int sink_base, const int m,
			   const uint16_t *ss_index, const segment *ss_seg, const int ss_len,
			   const segment *seg, const int seg_len, 
			   const int source_base, const int *source_n, 
			   const float mmass, float e, int *ncut, int q);
EXTERN void grav_qns_CUDA(const char *routine, const int base, const int m,
			  const int *source_list, const int source_n, int q);
EXTERN void grav_qnss_CUDA(const char *routine, const int base, const int m,
			   const segment *seg, const float *source, const int source_n, int q);

#else
/* stubs */
inline static void CUDA_Init(void *device_prop) {}
inline static void CUDA_Finalize(void) {}
inline static int vecwidthCUDA(void) {return 0;}
inline static int qallocCUDA(int *inuse) {return 0;}
inline static void WalkInitSrcCUDA(float *qtab, int stride, int64_t nobj) {}
inline static void WalkTerminateSrcCUDA(void) {}
inline static void WalkInitSinkCUDA(float *btab, int stride, int64_t nobj) {}
inline static void WalkTerminateSinkCUDA(float *btab, int stride, int64_t nobj) {};
inline static void grav_mn_CUDA(const char *routine, const float *p, float *accp, const int n, const int stride, 
				const float *f, const int source_n, const int sz, 
				const float e, int *ncut, int q) {}
inline static void grav_mns_CUDA(const char *routine, const int base, const int m,
				 const segment *seg, const int seg_n, const int source_n, 
				 const float mmass, float e, int *ncut, int q) {}
inline static void grav_mnss_CUDA(const char *routine, const int sink_base, const int m,
				  const uint16_t *ss_index, const segment *ss_seg, const int ss_len,
				  const segment *seg, const int seg_len, 
				  const int source_base, const int *source_n, 
				  const float mmass, float e, int *ncut, int q) {}
inline static void grav_qns_CUDA(const char *routine, const int base, const int m,
				 const int *source_list, const int source_n) {};
inline static void grav_qnss_CUDA(const char *routine, const int base, const int m,
				  const segment *seg, const float *source, const int source_n, int q) {}
#endif
