/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* Prototypes for all user accessible RANLIB routines */

extern void advnst(long k);
extern float genbet(float aa,float bb);
extern float genchi(float df);
extern float genexp(float av);
extern float genf(float dfn, float dfd);
extern float gengam(float a,float r);
extern void genmn(float *parm,float *x,float *work);
extern void genmul(long n,float *p,long ncat,long *ix);
extern float gennch(float df,float xnonc);
extern float gennf(float dfn, float dfd, float xnonc);
extern float gennor(float av,float sd);
extern void genprm(long *iarray,int larray);
extern float genunf(float low,float high);
extern void getsd(long *iseed1,long *iseed2);
extern void gscgn(long getset,long *g);
extern long ignbin(long n,float pp);
extern long ignnbn(long n,float p);
extern long ignlgi(void);
extern long ignpoi(float mu);
extern long ignuin(long low,long high);
extern void initgn(long isdtyp);
extern long mltmod(long a,long s,long m);
extern void phrtsd(char* phrase,long* seed1,long* seed2);
extern float ranf(void);
extern void setall(long iseed1,long iseed2);
extern void setant(long qvalue);
extern void setgmn(float *meanv,float *covm,long p,float *parm);
extern void setsd(long iseed1,long iseed2);
extern float sexpo(void);
extern float sgamma(float a);
extern float snorm(void);
