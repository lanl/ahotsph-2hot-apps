/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

void UpdateX(float *xptr, int xstride, float *yptr, int ystride, int n, 
	     float dt, float h);
void UpdateSX(float *xptr, int xstride, float *yptr, int ystride, int n, 
	float dt, float h);
void UpdateSXd(double *xptr, int xstride, double *yptr, int ystride, int n, 
	float dt, float h);
void ABUpdateX(float *x, int xstride, float *xdot, int xdot_stride, 
	       float *xdot_last, int xdot_last_stride, int n, 
	       float dt, float h);
void PUpdateV(float *v, int vstride, float *x, int xstride, 
	      float *xlast, int xlast_stride, float *xddot,
	      int xddot_stride, int n, float dt, float h);
void PUpdateVd(float *v, int vstride, double *x, int xstride, 
	      double *xlast, int xlast_stride, float *xddot,
	      int xddot_stride, int n, float dt, float h);
void PUpdateX(float *x, int xstride, float *xlast, int xlast_stride, 
	      float *xddot, int xddot_stride, int n, float dt, float h);
void PUpdateXd(double *x, int xstride, double *xlast, int xlast_stride, 
	      float *xddot, int xddot_stride, int n, float dt, float h);
void UpdateXs(float *xptr, int xstride, float *yptr, int ystride, 
	unsigned int *select, int select_stride, int n, float dt, float h);
void UpdateSXs(float *xptr, int xstride, float *yptr, int ystride, 
	unsigned int *select, int select_stride, int n, float dt, float h);
void ABUpdateXs(float *x, int xstride, float *xdot, int xdot_stride, 
	   float *xdot_last, int xdot_last_stride,
	   unsigned int *select, int select_stride, int n, float dt, float h);
void PUpdateVs(float *v, int vstride, float *x, int xstride, 
	 float *xlast, int xlast_stride, float *xddot, int xddot_stride,
	  unsigned int *select, int select_stride, int n, float dt, float h);
void PUpdateXs(float *x, int xstride, float *xlast, int xlast_stride, 
	  float *xddot, int xddot_stride, 
	  unsigned int *select, int select_stride, int n, float dt, float h);

void SetBoundary(int bound);
