#include "ndim.h"
#include "integrate.h"
#include "vop.h"

static int boundary_bit;

void
SetBoundary(int bound)
{
    boundary_bit = bound;
}

/* Leapfrog, 2nd order for 1st order ode, vector variable */

void
UpdateX(float *xptr, int xstride, float *yptr, int ystride, int n, 
	float dt, float h)
{
    float *end = xptr + n * xstride;
    float dt_h = (float)0.5*(dt+h);

    while (xptr < end) {
	VV(xptr, += dt_h * yptr);
	xptr += xstride;
	yptr += ystride;
    }
}

void
UpdateSX(float *xptr, int xstride, float *yptr, int ystride, int n, 
	float dt, float h)
{
    float *end = xptr + n * xstride;
    float dt_h = (float)0.5*(dt+h);

    while (xptr < end) {
	*xptr += dt_h * *yptr;
	xptr += xstride;
	yptr += ystride;
    }
}

void
UpdateSXd(double *xptr, int xstride, double *yptr, int ystride, int n, 
	float dt, float h)
{
    double *end = xptr + n * xstride;
    float dt_h = (float)0.5*(dt+h);

    while (xptr < end) {
	*xptr += dt_h * *yptr;
	xptr += xstride;
	yptr += ystride;
    }
}


/* Press method, 2nd order for 2nd order ode, vector variable */
/* Equivalent to leapfrog */

void
PUpdateX(float *x, int xstride, float *xlast, int xlast_stride, 
	 float *xddot, int xddot_stride, int n, float dt, float h)
{
    float *end = x + n * xstride;
    float tmp[NDIM];
    float half_dt2_plus_hdt = (float)0.5*dt*(dt+h);
    float dt_h = dt/h;

    while (x < end) {
	VVV(tmp, = x, - xlast);
	VV(xlast, = x);		/* Don't swap this line */
	VVV(x, += dt_h * tmp, + half_dt2_plus_hdt * xddot);
	x += xstride;
	xlast += xlast_stride;
	xddot += xddot_stride;
    }
}

void
PUpdateXd(double *x, int xstride, double *xlast, int xlast_stride, 
	 float *xddot, int xddot_stride, int n, float dt, float h)
{
    double *end = x + n * xstride;
    float tmp[NDIM];
    float half_dt2_plus_hdt = (float)0.5*dt*(dt+h);
    float dt_h = dt/h;

    while (x < end) {
	VVV(tmp, = x, - xlast);
	VV(xlast, = x);		/* Don't swap this line */
	VVV(x, += dt_h * tmp, + half_dt2_plus_hdt * xddot);
	x += xstride;
	xlast += xlast_stride;
	xddot += xddot_stride;
    }
}

/* Press method, 2nd order for xdot of 2nd order ode, vector variable */
/* Equivalent to leapfrog */

void
PUpdateV(float *v, int vstride, float *x, int xstride, 
	 float *xlast, int xlast_stride, float *xddot,
	 int xddot_stride, int n, float dt, float h)
{
    float *end = v + n * xstride;
    float tmp[NDIM];
    float h_2dt = (float)0.5*h+dt;
    float oneoh = 1.0 / h;

    while (v < end) {
	VVV(tmp, = x, - xlast);
	VVV(v, = oneoh * tmp, + h_2dt * xddot);
	v += vstride;
	x += xstride;
	xlast += xlast_stride;
	xddot += xddot_stride;
    }
}

void
PUpdateVd(float *v, int vstride, double *x, int xstride, 
	  double *xlast, int xlast_stride, float *xddot,
	  int xddot_stride, int n, float dt, float h)
{
    float *end = v + n * vstride;
    float tmp[NDIM];
    float h_2dt = (float)0.5*h+dt;
    float oneoh = 1.0 / h;

    while (v < end) {
	VVV(tmp, = x, - xlast);
	VVV(v, = oneoh * tmp, + h_2dt * xddot);
	v += vstride;
	x += xstride;
	xlast += xlast_stride;
	xddot += xddot_stride;
    }
}

/* Adams-Bashforth 2nd order for 1st order ode, scalar variable */
void
ABUpdateX(float *x, int xstride, float *xdot, int xdot_stride, 
	  float *xdot_last, int xdot_last_stride, int n, float dt, float h)
{
    float *end = x + n * xstride;
    float c1 = dt*((float)1.0 + (float)0.5*dt/h);
    float c2 = (float)0.5*dt*dt/h;
    while (x < end) {
	*x += c1 * *xdot - c2 * *xdot_last;
	*xdot_last = *xdot;
	x += xstride;
	xdot += xdot_stride;
	xdot_last += xdot_last_stride;
    }
}

/* Leapfrog, 2nd order for 1st order ode, vector variable */
/* With selection vector */

void
UpdateXs(float *xptr, int xstride, float *yptr, int ystride, 
	unsigned int *select, int select_stride, int n, float dt, float h)
{
    float *end = xptr + n * xstride;
    float dt_h = (float)0.5*(dt+h);

    while (xptr < end) {
	if (!(*select & (1 << boundary_bit))) {
	    VV(xptr, += dt_h * yptr);
	}
	xptr += xstride;
	yptr += ystride;
	select += select_stride;
    }
}

/* Leapfrog, 2nd order for 1st order ode, scalar variable */
/* With selection vector */

void
UpdateSXs(float *xptr, int xstride, float *yptr, int ystride, 
	unsigned int *select, int select_stride, int n, float dt, float h)
{
    float *end = xptr + n * xstride;
    float dt_h = (float)0.5*(dt+h);

    while (xptr < end) {
	if (!(*select & (1 << boundary_bit))) {
	    *xptr += dt_h * *yptr;
	}
	xptr += xstride;
	yptr += ystride;
	select += select_stride;
    }
}

/* Press method, 2nd order for 2nd order ode, vector variable */
/* Equivalent to leapfrog */

void
PUpdateXs(float *x, int xstride, float *xlast, int xlast_stride, 
	  float *xddot, int xddot_stride, 
	  unsigned int *select, int select_stride, int n, float dt, float h)
{
    float *end = x + n * xstride;
    float tmp[NDIM];
    float half_dt2_plus_hdt = (float)0.5*dt*(dt+h);
    float dt_h = dt/h;

    while (x < end) {
	if (!(*select & (1 << boundary_bit))) {
	    VVV(tmp, = x, - xlast);
	    VV(xlast, = x);
	    VVV(x, += dt_h * tmp, + half_dt2_plus_hdt * xddot);
	}
	x += xstride;
	xlast += xlast_stride;
	xddot += xddot_stride;
	select += select_stride;
    }
}

/* Press method, 2nd order for xdot of 2nd order ode, vector variable */
/* Equivalent to leapfrog */

void
PUpdateVs(float *v, int vstride, float *x, int xstride, 
	  float *xlast, int xlast_stride, float *xddot, int xddot_stride,
	  unsigned int *select, int select_stride, int n, float dt, float h)
{
    float *end = v + n * xstride;
    float tmp[NDIM];
    float h_2dt = (float)0.5*h+dt;
    float oneoh = 1.0 / h;

    while (v < end) {
	if (!(*select & (1 << boundary_bit))) {
	    VVV(tmp, = x, - xlast);
	    VVV(v, = oneoh * tmp, + h_2dt * xddot);
	}
	v += vstride;
	x += xstride;
	xlast += xlast_stride;
	xddot += xddot_stride;
	select += select_stride;
    }
}

/* Adams-Bashforth 2nd order for 1st order ode, scalar variable */
void
ABUpdateXs(float *x, int xstride, float *xdot, int xdot_stride, 
	   float *xdot_last, int xdot_last_stride,
	   unsigned int *select, int select_stride, int n, float dt, float h)
{
    float *end = x + n * xstride;
    float c1 = dt*((float)1.0 + (float)0.5*dt/h);
    float c2 = (float)0.5*dt*dt/h;

    while (x < end) {
	if (!(*select & (1 << boundary_bit))) {
	    *x += c1 * *xdot - c2 * *xdot_last;
	    *xdot_last = *xdot;
	}
	x += xstride;
	xdot += xdot_stride;
	xdot_last += xdot_last_stride;
	select += select_stride;
    }
}
