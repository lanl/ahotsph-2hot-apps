void UpdateX(float *xptr, int xstride, float *yptr, int ystride, int n, 
	     float dt, float h);
void UpdateSX(float *xptr, int xstride, float *yptr, int ystride, int n, 
	float dt, float h);
void ABUpdateX(float *x, int xstride, float *xdot, int xdot_stride, 
	       float *xdot_last, int xdot_last_stride, int n, 
	       float dt, float h);
void PUpdateV(float *v, int vstride, float *x, int xstride, 
	      float *xlast, int xlast_stride, float *xddot,
	      int xddot_stride, int n, float dt, float h);
void PUpdateX(float *x, int xstride, float *xlast, int xlast_stride, 
	      float *xddot, int xddot_stride, int n, float dt, float h);
void UpdateXs(float *xptr, int xstride, float *yptr, int ystride, 
	unsigned int *select, int select_stride, int n, float dt, float h);
void UpdateSXs(float *xptr, int xstride, float *yptr, int ystride, 
	unsigned int *select, int select_stride, int n, float dt, float h);
void ABUpdateXs(float *x, int xstride, float *xdot, int xdot_stride, 
	   float *xdot_last, int xdot_last_stride,
	   unsigned int *select, int select_stride, int n, float dt, float h);
void ABUpdateVs(float *x, int xstride, float *xdot, int xdot_stride, 
	   float *xdot_last, int xdot_last_stride,
	   unsigned int *select, int select_stride, int n, float dt, float h);
void PUpdateVs(float *v, int vstride, float *x, int xstride, 
	 float *xlast, int xlast_stride, float *xddot, int xddot_stride,
	  unsigned int *select, int select_stride, int n, float dt, float h);
void PUpdateXs(float *x, int xstride, float *xlast, int xlast_stride, 
	  float *xddot, int xddot_stride, 
	  unsigned int *select, int select_stride, int n, float dt, float h);


