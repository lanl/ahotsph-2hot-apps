#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "physics.h"
#include "vop.h"
#include "moments.h"
#include "Malloc.h"
#include "mpmy.h"
#include "Msgs.h"

static double ewx4, ewx2y2, ewx6, ewx4y2, ewx2y2z2, ewx8, ewx6y2, ewx4y4;
static cartesian_moments CQ;
static int subtract_background;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void
ewald_background(const float xx[NDIM], float mass, double phicorr, double acc[NDIM], double *phi)
{
    const cartesian_moments *q = &CQ;
    float r[NDIM];

    r[0] = xx[0]-q->x; 
    r[1] = xx[1]-q->y; 
    r[2] = xx[2]-q->z;

    acc[0] = 4.188790204786*r[0];
    acc[1] = 4.188790204786*r[1];
    acc[2] = 4.188790204786*r[2];
    *phi = phicorr + 2.0*1.4186487533 - (2.0*M_PI/3.0)*(Dot(xx, xx) - 2.0*xx[0]*q->x - 2.0*xx[1]*q->y - 2.0*xx[2]*q->z + q->x2 + q->y2 + q->z2);
    /* mass terms are for convention without self-interactions in images */
    *phi -= 2.0*1.4186487533*mass/q->m;
#if 0
    *phi += (2.0*M_PI/3.0)*Dot(xx, xx)*mass/q->m;
#endif
 }

void o4(double *a0, double *a1, double *a2, double *pp, float x, float y, float z, float *Q)
{
    float x2, y2, z2;
    float x3, y3, z3;
    float x4, y4, z4;
    float t;

    x2 = x*x; y2 = y*y; z2 = z*z;
    x3 = x2*x; y3 = y2*y; z3 = z2*z;
    x4 = x3*x; y4 = y3*y; z4 = z3*z;

    t = 5.25*ewx2y2 - 1.75*ewx4;

    *pp += t * (Q[1] + Q[4]*x3 + x4*Q[0] - 0.6666666666666666*Q[7]*y3 + y4*Q[0] + Q[9]*z + y2*(Q[8] + Q[11]*z - 3.*Q[0]*z2) + x2*(Q[3] + Q[7]*y - 3.*Q[0]*y2 + Q[11]*z - 3.*Q[0]*z2) + Q[13]*z2 + x*(Q[2] + Q[6]*y - 1.5*Q[4]*y2 + Q[10]*z - 1.5*Q[4]*z2) + y*(Q[5] + Q[12]*z + Q[7]*z2) - 0.6666666666666666*Q[11]*z3 + z4*Q[0]);
    *a0 -= t * (Q[2] + 3.*Q[4]*x2 + 4.*Q[0]*x3 + Q[6]*y - 1.5*Q[4]*y2 + Q[10]*z + 2.*x*(Q[3] + Q[7]*y - 3.*Q[0]*y2 + Q[11]*z - 3.*Q[0]*z2) - 1.5*Q[4]*z2);
    *a1 -= t * (Q[5] + x2*(Q[7] - 6.*Q[0]*y) + x*(Q[6] - 3.*Q[4]*y) - 2.*Q[7]*y2 + 4.*Q[0]*y3 + Q[12]*z + 2.*y*(Q[8] + Q[11]*z - 3.*Q[0]*z2) + Q[7]*z2);
    *a2 -= t * (Q[9] + x2*(Q[11] - 6.*Q[0]*z) + y2*(Q[11] - 6.*Q[0]*z) + 2.*Q[13]*z + x*(Q[10] - 3.*Q[4]*z) + y*(Q[12] + 2.*Q[7]*z) - 2.*Q[11]*z2 + 4.*Q[0]*z3);
}

void o6(double *a0, double *a1, double *a2, double *pp, float x, float y, float z, float *Q)
{
    float x2, y2, z2;
    float x3, y3, z3;
    float x4, y4, z4;
    float x5, y5, z5;
    float x6, y6, z6;
    float t;

    x2 = x*x; y2 = y*y; z2 = z*z;
    x3 = x2*x; y3 = y2*y; z3 = z2*z;
    x4 = x3*x; y4 = y3*y; z4 = z3*z;
    x5 = x4*x; y5 = y4*y; z5 = z4*z;
    x6 = x5*x; y6 = y5*y; z6 = z5*z;

    t = -11.25*ewx2y2z2 + 5.625*ewx4y2 - 0.375*ewx6;

    *pp += t * (Q[14] + 1.5*Q[4]*x5 + x6*Q[0] - 1.*Q[7]*y5 + y6*Q[0] + Q[22]*z + Q[28]*z2 - 2.5*y4*(-1.*Q[8] - 1.*Q[11]*z + 3.*Q[0]*z2) - 2.5*x4*(-1.*Q[3] - 1.*Q[7]*y + 3.*Q[0]*y2 - 1.*Q[11]*z + 3.*Q[0]*z2) - 2.5*x3*(-2.*Q[2] - 2.*Q[6]*y + 3.*Q[4]*y2 - 2.*Q[10]*z + 3.*Q[4]*z2) + 5.*y3*(Q[5] + Q[12]*z + Q[7]*z2) + 5.*Q[9]*z3 + y2*(Q[20] + Q[27]*z + 15.*Q[3]*z2 + 5.*Q[11]*z3 - 7.5*z4) + x2*(Q[16] + 5.*Q[7]*y3 - 7.5*y4 + Q[24]*z + 15.*Q[8]*z2 + 15.*y2*(Q[13] - 2.*Q[11]*z + 6.*z2) + y*(Q[19] - 30.*Q[12]*z - 30.*Q[7]*z2) + 5.*Q[11]*z3 - 7.5*z4) + 2.5*Q[13]*z4 + x*(Q[15] + 5.*Q[6]*y3 - 3.75*Q[4]*y4 + Q[23]*z + y*(Q[18] + z*(Q[26] - 30.*Q[6]*z)) + Q[29]*z2 + y2*(Q[21] - 30.*Q[10]*z + 45.*Q[4]*z2) + 5.*Q[10]*z3 - 3.75*Q[4]*z4) + y*(Q[17] + Q[25]*z + Q[30]*z2 + 5.*Q[12]*z3 + 2.5*Q[7]*z4) - 1.*Q[11]*z5 + z6*Q[0]);
    *a0 -= t * (Q[15] + 7.5*Q[4]*x4 + 6.*Q[0]*x5 + 5.*Q[6]*y3 - 3.75*Q[4]*y4 + Q[23]*z + y*(Q[18] + z*(Q[26] - 30.*Q[6]*z)) + Q[29]*z2 - 10.*x3*(-1.*Q[3] - 1.*Q[7]*y + 3.*Q[0]*y2 - 1.*Q[11]*z + 3.*Q[0]*z2) - 7.5*x2*(-2.*Q[2] - 2.*Q[6]*y + 3.*Q[4]*y2 - 2.*Q[10]*z + 3.*Q[4]*z2) + y2*(Q[21] - 30.*Q[10]*z + 45.*Q[4]*z2) + 5.*Q[10]*z3 + x*(2.*Q[16] + 10.*Q[7]*y3 - 15.*Q[0]*y4 + 2.*Q[24]*z + 30.*Q[8]*z2 + 30.*y2*(Q[13] - 2.*Q[11]*z + 6.*Q[0]*z2) + y*(2.*Q[19] - 60.*Q[12]*z - 60.*Q[7]*z2) + 10.*Q[11]*z3 - 15.*Q[0]*z4) - 3.75*Q[4]*z4);
    *a1 -= t * (Q[17] + 2.5*x4*(Q[7] - 6.*Q[0]*y) + 5.*x3*(Q[6] - 3.*Q[4]*y) - 5.*Q[7]*y4 + 6.*Q[0]*y5 + Q[25]*z + 10.*y3*(Q[8] + Q[11]*z - 3.*Q[0]*z2) + Q[30]*z2 + 15.*y2*(Q[5] + Q[12]*z + Q[7]*z2) + x2*(Q[19] + 15.*Q[7]*y2 - 30.*Q[0]*y3 - 30.*Q[12]*z - 30.*Q[7]*z2 + 30.*y*(Q[13] - 2.*Q[11]*z + 6.*Q[0]*z2)) + x*(Q[18] + 15.*Q[6]*y2 - 15.*Q[4]*y3 + z*(Q[26] - 30.*Q[6]*z) + 2.*y*(Q[21] - 30.*Q[10]*z + 45.*Q[4]*z2)) + 5.*Q[12]*z3 + y*(2.*Q[20] + 2.*Q[27]*z + 30.*Q[3]*z2 + 10.*Q[11]*z3 - 15.*Q[0]*z4) + 2.5*Q[7]*z4);
    *a2 -= t * (Q[22] + 2.5*x4*(Q[11] - 6.*Q[0]*z) + 2.5*y4*(Q[11] - 6.*Q[0]*z) + 2.*Q[28]*z + 5.*x3*(Q[10] - 3.*Q[4]*z) + 5.*y3*(Q[12] + 2.*Q[7]*z) + 15.*Q[9]*z2 + y2*(Q[27] + 30.*Q[3]*z + 15.*Q[11]*z2 - 30.*Q[0]*z3) + x2*(Q[24] + 30.*Q[8]*z + 30.*y2*(-1.*Q[11] + 6.*Q[0]*z) - 30.*y*(Q[12] + 2.*Q[7]*z) + 15.*Q[11]*z2 - 30.*Q[0]*z3) + 10.*Q[13]*z3 + x*(Q[23] + 2.*Q[29]*z - 30.*y2*(Q[10] - 3.*Q[4]*z) + y*(Q[26] - 60.*Q[6]*z) + 15.*Q[10]*z2 - 15.*Q[4]*z3) + y*(Q[25] + 2.*Q[30]*z + 15.*Q[12]*z2 + 10.*Q[7]*z3) - 5.*Q[11]*z4 + 6.*Q[0]*z5);
}

void o8(double *a0, double *a1, double *a2, double *pp, float x, float y, float z, float *Q)
{
    float x2, y2, z2;
    float x3, y3, z3;
    float x4, y4, z4;
    float x5, y5, z5;
    float x6, y6, z6;
    float x7, y7, z7;
    float x8, y8, z8;
    float t;

    x2 = x*x; y2 = y*y; z2 = z*z;
    x3 = x2*x; y3 = y2*y; z3 = z2*z;
    x4 = x3*x; y4 = y3*y; z4 = z3*z;
    x5 = x4*x; y5 = y4*y; z5 = z4*z;
    x6 = x5*x; y6 = y5*y; z6 = z5*z;
    x7 = x6*x; y7 = y6*y; z7 = z6*z;
    x8 = x7*x; y8 = y7*y; z8 = z7*z;

    t = -54.140625*ewx4y4 + 43.3125*ewx6y2 - 1.546875*ewx8;

    *pp += t * (Q[31] + 2.*Q[4]*x7 + x8*Q[0] - 1.3333333333333333*Q[7]*y7 + y8*Q[0] + Q[48]*z + Q[57]*z2 - 4.666666666666667*y6*(-1.*Q[8] - 1.*Q[11]*z + 3.*Q[0]*z2) - 4.666666666666667*x6*(-1.*Q[3] - 1.*Q[7]*y + 3.*Q[0]*y2 - 1.*Q[11]*z + 3.*Q[0]*z2) - 7.*x5*(-2.*Q[2] - 2.*Q[6]*y + 3.*Q[4]*y2 - 2.*Q[10]*z + 3.*Q[4]*z2) + 14.*y5*(Q[5] + Q[12]*z + Q[7]*z2) + Q[66]*z3 + Q[67]*z4 + x4*(Q[35] + Q[40]*y + Q[45]*y2 - 23.333333333333332*Q[7]*y3 + 35.*Q[0]*y4 + Q[52]*z + Q[61]*z2 - 23.333333333333332*Q[11]*z3 + 35.*Q[0]*z4) + y4*(Q[47] + Q[56]*z + Q[65]*z2 - 23.333333333333332*Q[11]*z3 + 35.*Q[0]*z4) + x3*(Q[34] + Q[39]*y + Q[44]*y2 - 46.666666666666664*Q[6]*y3 + 35.*Q[4]*y4 + Q[51]*z + Q[60]*z2 - 46.666666666666664*Q[10]*z3 + 35.*Q[4]*z4) + y3*(Q[46] + Q[55]*z + Q[64]*z2 - 46.666666666666664*Q[12]*z3 - 23.333333333333332*Q[7]*z4) + 14.*Q[9]*z5 + x2*(Q[33] + Q[38]*y + Q[43]*y2 - 2.*Q[40]*y3 - 1.*Q[45]*y4 + 14.*Q[7]*y5 - 14.*Q[0]*y6 + Q[50]*z + Q[59]*z2 - 2.*Q[52]*z3 - 1.*Q[61]*z4 + 14.*Q[11]*z5 - 14.*Q[0]*z6) + y2*(Q[41] + Q[54]*z + Q[63]*z2 - 2.*Q[56]*z3 - 1.*Q[65]*z4 + 14.*Q[11]*z5 - 14.*Q[0]*z6) + 4.666666666666667*Q[13]*z6 + x*(Q[32] + Q[37]*y + Q[42]*y2 - 1.*Q[39]*y3 - 0.5*Q[44]*y4 + 14.*Q[6]*y5 - 7.*Q[4]*y6 + Q[49]*z + Q[58]*z2 - 1.*Q[51]*z3 - 0.5*Q[60]*z4 + 14.*Q[10]*z5 - 7.*Q[4]*z6) + y*(Q[36] + Q[53]*z + Q[62]*z2 - 1.*Q[55]*z3 - 0.5*Q[64]*z4 + 14.*Q[12]*z5 + 4.666666666666667*Q[7]*z6) - 1.3333333333333333*Q[11]*z7 + z8*Q[0]);
    *a0 -= t * (Q[32] + 14.*Q[4]*x6 + 8.*Q[0]*x7 + Q[37]*y + Q[42]*y2 - 1.*Q[39]*y3 - 0.5*Q[44]*y4 + 14.*Q[6]*y5 - 7.*Q[4]*y6 + Q[49]*z + Q[58]*z2 - 28.*x5*(-1.*Q[3] - 1.*Q[7]*y + 3.*Q[0]*y2 - 1.*Q[11]*z + 3.*Q[0]*z2) - 35.*x4*(-2.*Q[2] - 2.*Q[6]*y + 3.*Q[4]*y2 - 2.*Q[10]*z + 3.*Q[4]*z2) - 1.*Q[51]*z3 - 0.5*Q[60]*z4 + 1.3333333333333333*x3*(3.*Q[35] + 3.*Q[40]*y + 3.*Q[45]*y2 - 70.*Q[7]*y3 + 105.*Q[0]*y4 + 3.*Q[52]*z + 3.*Q[61]*z2 - 70.*Q[11]*z3 + 105.*Q[0]*z4) + x2*(3.*Q[34] + 3.*Q[39]*y + 3.*Q[44]*y2 - 140.*Q[6]*y3 + 105.*Q[4]*y4 + 3.*Q[51]*z + 3.*Q[60]*z2 - 140.*Q[10]*z3 + 105.*Q[4]*z4) + 14.*Q[10]*z5 + 2.*x*(Q[33] + Q[38]*y + Q[43]*y2 - 2.*Q[40]*y3 - 1.*Q[45]*y4 + 14.*Q[7]*y5 - 14.*Q[0]*y6 + Q[50]*z + Q[59]*z2 - 2.*Q[52]*z3 - 1.*Q[61]*z4 + 14.*Q[11]*z5 - 14.*Q[0]*z6) - 7.*Q[4]*z6);
    *a1 -= t * (Q[36] + 4.666666666666667*x6*(Q[7] - 6.*Q[0]*y) + 14.*x5*(Q[6] - 3.*Q[4]*y) + x3*(Q[39] + 2.*y*(Q[44] + 70.*y*(-1.*Q[6] + Q[4]*y))) + x4*(Q[40] + 2.*Q[45]*y - 70.*Q[7]*y2 + 140.*Q[0]*y3) + x*(Q[37] + y*(2.*Q[42] - 3.*Q[39]*y - 2.*Q[44]*y2 + 70.*Q[6]*y3 - 42.*Q[4]*y4)) + x2*(Q[38] + 2.*Q[43]*y - 6.*Q[40]*y2 - 4.*Q[45]*y3 + 70.*Q[7]*y4 - 84.*Q[0]*y5) - 9.333333333333334*Q[7]*y6 + 8.*Q[0]*y7 + Q[53]*z + 28.*y5*(Q[8] + Q[11]*z - 3.*Q[0]*z2) + Q[62]*z2 + 70.*y4*(Q[5] + Q[12]*z + Q[7]*z2) - 1.*Q[55]*z3 - 0.5*Q[64]*z4 + 1.3333333333333333*y3*(3.*Q[47] + 3.*Q[56]*z + 3.*Q[65]*z2 - 70.*Q[11]*z3 + 105.*Q[0]*z4) + y2*(3.*Q[46] + 3.*Q[55]*z + 3.*Q[64]*z2 - 140.*Q[12]*z3 - 70.*Q[7]*z4) + 14.*Q[12]*z5 + 2.*y*(Q[41] + Q[54]*z + Q[63]*z2 - 2.*Q[56]*z3 - 1.*Q[65]*z4 + 14.*Q[11]*z5 - 14.*Q[0]*z6) + 4.666666666666667*Q[7]*z6);
    *a2 -= t * (Q[48] + 4.666666666666667*x6*(Q[11] - 6.*Q[0]*z) + 4.666666666666667*y6*(Q[11] - 6.*Q[0]*z) + 2.*Q[57]*z + 14.*x5*(Q[10] - 3.*Q[4]*z) + 14.*y5*(Q[12] + 2.*Q[7]*z) + 3.*Q[66]*z2 + 4.*Q[67]*z3 + x4*(Q[52] + 2.*Q[61]*z - 70.*Q[11]*z2 + 140.*Q[0]*z3) + y4*(Q[56] + 2.*Q[65]*z - 70.*Q[11]*z2 + 140.*Q[0]*z3) + x3*(Q[51] + 2.*Q[60]*z - 140.*Q[10]*z2 + 140.*Q[4]*z3) + y3*(Q[55] + 2.*Q[64]*z - 140.*Q[12]*z2 - 93.33333333333333*Q[7]*z3) + 70.*Q[9]*z4 + x2*(Q[50] + 2.*Q[59]*z - 6.*Q[52]*z2 - 4.*Q[61]*z3 + 70.*Q[11]*z4 - 84.*Q[0]*z5) + y2*(Q[54] + 2.*Q[63]*z - 6.*Q[56]*z2 - 4.*Q[65]*z3 + 70.*Q[11]*z4 - 84.*Q[0]*z5) + 28.*Q[13]*z5 + x*(Q[49] + 2.*Q[58]*z - 3.*Q[51]*z2 - 2.*Q[60]*z3 + 70.*Q[10]*z4 - 42.*Q[4]*z5) + y*(Q[53] + 2.*Q[62]*z - 3.*Q[55]*z2 - 2.*Q[64]*z3 + 70.*Q[12]*z4 + 28.*Q[7]*z5) - 9.333333333333334*Q[11]*z6 + 8.*Q[0]*z7);
}


void 
ewald_le(const float xx[NDIM], double acc[NDIM], double *phi, float *Q, int nimage)
{
    float x, y, z;
    double a0 = 0.0, a1 = 0.0, a2 = 0.0, pp = 0.0;

    x = xx[0];
    y = xx[1];
    z = xx[2];

    if (nimage == 2) {
	/* dps NIMAGE NN L 64 2 2000 1.000000 */
	ewx4     = 1.4396314829106346e-01;
	ewx2y2   = 4.0171409156351564e-02;
	ewx6     = 6.5766562442844188e-03;
	ewx4y2   = 8.6523239706610977e-04;
	ewx2y2z2 = 2.2149598086679247e-04;
	ewx8     = 4.5164568468964107e-04;
	ewx6y2   = 3.6866246504650772e-05;
	ewx4y4   = 1.6274461485715282e-05;
    } else {
	ewx4     = 3.753808257962605e-01;
	ewx2y2   = 1.028100333036129e-01;
	ewx6     = 4.309234766008302e-02;
	ewx4y2   = 5.441256840176889e-03;
	ewx2y2z2 = 1.351866101467558e-03;
	ewx8     = 7.289731443807441e-03;
	ewx6y2   = 5.254647680142501e-04;
	ewx4y4   = 2.677806065828672e-04;
    }

    o8(&a0, &a1, &a2, &pp, x, y, z, Q);
    o6(&a0, &a1, &a2, &pp, x, y, z, Q);
    o4(&a0, &a1, &a2, &pp, x, y, z, Q);
    if (subtract_background) {
	a0 -= -(4.0*M_PI/3.0)*Q[4]/4;
	a1 -= (4.0*M_PI/3.0)*Q[7]/6;
	a2 -= (4.0*M_PI/3.0)*Q[11]/6;
    }

    acc[0] = a0;
    acc[1] = a1;
    acc[2] = a2;
    *phi = pp;
}

void
calculate_cartesian_moments(body *btab, int nobj, double L, float *Q, int msb)
{
    body *p;
    double m;
    double x, y, z;
    double x2, y2, z2;
    double x3, y3, z3;
    double x4, y4, z4;
    double x5, y5, z5;
    double x6, y6, z6;
    double x7, y7, z7;
    double x8, y8, z8;
    cartesian_moments_d *q;
    cartesian_moments *c = &CQ;
    double Linv = 1.0/L;

    subtract_background = msb;
    q = Calloc(1, sizeof(cartesian_moments_d));

    for (p = btab; p < btab+nobj; p++) {
	m = p->mass;
	q->m += m;

	x = Linv*p->pos[0];
	y = Linv*p->pos[1];
	z = Linv*p->pos[2];
	
	q->x += m*x;
	q->y += m*y;
	q->z += m*z;

	x2 = x*x; y2 = y*y; z2 = z*z;
	q->x2 += m*x2;
	q->xy += m*x*y;
	q->y2 += m*y2;
	q->xz += m*x*z;
	q->yz += m*y*z;
	q->z2 += m*z2;
	
	x3 = x2*x; y3 = y2*y; z3 = z2*z;
	q->x3 += m*x3;
	q->x2y += m*x2*y;
	q->xy2 += m*x*y2;
	q->y3 += m*y3;
	q->x2z += m*x2*z;
	q->xyz += m*x*y*z;
	q->y2z += m*y2*z;
	q->xz2 += m*x*z2;
	q->yz2 += m*y*z2;
	q->z3 += m*z3;
	
	x4 = x3*x; y4 = y3*y; z4 = z3*z;
	q->x4 += m*x4;
	q->x3y += m*x3*y;
	q->x2y2 += m*x2*y2;
	q->xy3 += m*x*y3;
	q->y4 += m*y4;
	q->x3z += m*x3*z;
	q->x2yz += m*x2*y*z;
	q->xy2z += m*x*y2*z;
	q->y3z += m*y3*z;
	q->x2z2 += m*x2*z2;
	q->xyz2 += m*x*y*z2;
	q->y2z2 += m*y2*z2;
	q->xz3 += m*x*z3;
	q->yz3 += m*y*z3;
	q->z4 += m*z4;
	
	x5 = x4*x; y5 = y4*y; z5 = z4*z;
	q->x5 += m*x5;
	q->x4y += m*x4*y;
	q->x3y2 += m*x3*y2;
	q->x2y3 += m*x2*y3;
	q->xy4 += m*x*y4;
	q->y5 += m*y5;
	q->x4z += m*x4*z;
	q->x3yz += m*x3*y*z;
	q->x2y2z += m*x2*y2*z;
	q->xy3z += m*x*y3*z;
	q->y4z += m*y4*z;
	q->x3z2 += m*x3*z2;
	q->x2yz2 += m*x2*y*z2;
	q->xy2z2 += m*x*y2*z2;
	q->y3z2 += m*y3*z2;
	q->x2z3 += m*x2*z3;
	q->xyz3 += m*x*y*z3;
	q->y2z3 += m*y2*z3;
	q->xz4 += m*x*z4;
	q->yz4 += m*y*z4;
	q->z5 += m*z5;
	
	x6 = x5*x; y6 = y5*y; z6 = z5*z;
	q->x6 += m*x6;
	q->x5y += m*x5*y;
	q->x4y2 += m*x4*y2;
	q->x3y3 += m*x3*y3;
	q->x2y4 += m*x2*y4;
	q->xy5 += m*x*y5;
	q->y6 += m*y6;
	q->x5z += m*x5*z;
	q->x4yz += m*x4*y*z;
	q->x3y2z += m*x3*y2*z;
	q->x2y3z += m*x2*y3*z;
	q->xy4z += m*x*y4*z;
	q->y5z += m*y5*z;
	q->x4z2 += m*x4*z2;
	q->x3yz2 += m*x3*y*z2;
	q->x2y2z2 += m*x2*y2*z2;
	q->xy3z2 += m*x*y3*z2;
	q->y4z2 += m*y4*z2;
	q->x3z3 += m*x3*z3;
	q->x2yz3 += m*x2*y*z3;
	q->xy2z3 += m*x*y2*z3;
	q->y3z3 += m*y3*z3;
	q->x2z4 += m*x2*z4;
	q->xyz4 += m*x*y*z4;
	q->y2z4 += m*y2*z4;
	q->xz5 += m*x*z5;
	q->yz5 += m*y*z5;
	q->z6 += m*z6;
	
	x7 = x6*x; y7 = y6*y; z7 = z6*z;
	q->x7 += m*x7;
	q->x6y += m*x6*y;
	q->x5y2 += m*x5*y2;
	q->x4y3 += m*x4*y3;
	q->x3y4 += m*x3*y4;
	q->x2y5 += m*x2*y5;
	q->xy6 += m*x*y6;
	q->y7 += m*y7;
	q->x6z += m*x6*z;
	q->x5yz += m*x5*y*z;
	q->x4y2z += m*x4*y2*z;
	q->x3y3z += m*x3*y3*z;
	q->x2y4z += m*x2*y4*z;
	q->xy5z += m*x*y5*z;
	q->y6z += m*y6*z;
	q->x5z2 += m*x5*z2;
	q->x4yz2 += m*x4*y*z2;
	q->x3y2z2 += m*x3*y2*z2;
	q->x2y3z2 += m*x2*y3*z2;
	q->xy4z2 += m*x*y4*z2;
	q->y5z2 += m*y5*z2;
	q->x4z3 += m*x4*z3;
	q->x3yz3 += m*x3*y*z3;
	q->x2y2z3 += m*x2*y2*z3;
	q->xy3z3 += m*x*y3*z3;
	q->y4z3 += m*y4*z3;
	q->x3z4 += m*x3*z4;
	q->x2yz4 += m*x2*y*z4;
	q->xy2z4 += m*x*y2*z4;
	q->y3z4 += m*y3*z4;
	q->x2z5 += m*x2*z5;
	q->xyz5 += m*x*y*z5;
	q->y2z5 += m*y2*z5;
	q->xz6 += m*x*z6;
	q->yz6 += m*y*z6;
	q->z7 += m*z7;
	
	x8 = x7*x; y8 = y7*y; z8 = z7*z;
	q->x8 += m*x8;
	q->x7y += m*x7*y;
	q->x6y2 += m*x6*y2;
	q->x5y3 += m*x5*y3;
	q->x4y4 += m*x4*y4;
	q->x3y5 += m*x3*y5;
	q->x2y6 += m*x2*y6;
	q->xy7 += m*x*y7;
	q->y8 += m*y8;
	q->x7z += m*x7*z;
	q->x6yz += m*x6*y*z;
	q->x5y2z += m*x5*y2*z;
	q->x4y3z += m*x4*y3*z;
	q->x3y4z += m*x3*y4*z;
	q->x2y5z += m*x2*y5*z;
	q->xy6z += m*x*y6*z;
	q->y7z += m*y7*z;
	q->x6z2 += m*x6*z2;
	q->x5yz2 += m*x5*y*z2;
	q->x4y2z2 += m*x4*y2*z2;
	q->x3y3z2 += m*x3*y3*z2;
	q->x2y4z2 += m*x2*y4*z2;
	q->xy5z2 += m*x*y5*z2;
	q->y6z2 += m*y6*z2;
	q->x5z3 += m*x5*z3;
	q->x4yz3 += m*x4*y*z3;
	q->x3y2z3 += m*x3*y2*z3;
	q->x2y3z3 += m*x2*y3*z3;
	q->xy4z3 += m*x*y4*z3;
	q->y5z3 += m*y5*z3;
	q->x4z4 += m*x4*z4;
	q->x3yz4 += m*x3*y*z4;
	q->x2y2z4 += m*x2*y2*z4;
	q->xy3z4 += m*x*y3*z4;
	q->y4z4 += m*y4*z4;
	q->x3z5 += m*x3*z5;
	q->x2yz5 += m*x2*y*z5;
	q->xy2z5 += m*x*y2*z5;
	q->y3z5 += m*y3*z5;
	q->x2z6 += m*x2*z6;
	q->xyz6 += m*x*y*z6;
	q->y2z6 += m*y2*z6;
	q->xz7 += m*x*z7;
	q->yz7 += m*y*z7;
	q->z8 += m*z8;
    }
    VV(c->pos, = q->pos);
    MPMY_Combine(q, q, sizeof(cartesian_moments_d)/sizeof(double), MPMY_DOUBLE, MPMY_SUM);
    if (subtract_background) {
	double a;

	a = -q->m;

	q->x2 += a/12.;
	q->y2 += a/12.;
	q->z2 += a/12.;

	q->x4 += a/80.;
	q->x2y2 += a/144.;
	q->y4 += a/80.;
	q->x2z2 += a/144.;
	q->y2z2 += a/144.;
	q->z4 += a/80.;

	q->x6 += a/448.;
	q->x4y2 += a/960.;
	q->x2y4 += a/960.;
	q->y6 += a/448.;
	q->x4z2 += a/960.;
	q->x2y2z2 += a/1728.;
	q->y4z2 += a/960.;
	q->x2z4 += a/960.;
	q->y2z4 += a/960.;
	q->z6 += a/448.;

	q->x8 += a/2304.;
	q->x6y2 += a/5376.;
	q->x4y4 += a/6400.;
	q->x2y6 += a/5376.;
	q->y8 += a/2304.;
	q->x6z2 += a/5376.;
	q->x4y2z2 += a/11520.;
	q->x2y4z2 += a/11520.;
	q->y6z2 += a/5376.;
	q->x4z4 += a/6400.;
	q->x2y2z4 += a/11520.;
	q->y4z4 += a/6400.;
	q->x2z6 += a/5376.;
	q->y2z6 += a/5376.;
	q->z8 += a/2304.;
    }
    {
	float *f;
	double *d;
	VV(c->pos, = q->pos);
	c->m = q->m;
	for (f = &c->x, d = &q->x; f <= &c->z8; f++, d++) {
	    *f = *d/q->m;
	}
    }
    Free(q);
    Msgf(("mass %g\n", c->m));
    Msgf(("dipole %g %g %g\n", c->x, c->y, c->z));
    Msgf(("quadrupole %g %g %g %g %g %g\n", c->x2, c->y2, c->z2, c->xy, c->xz, c->yz));
    Msgf(("octupole %g %g %g ...\n", c->x3, c->x2y, c->xy2));
    Msgf(("hexadecapole %g %g %g %g ...\n", 
	  c->x4, c->x3y, c->x2y2, c->xy3));
    Msgf(("5-pole %g %g %g %g %g ...\n", 
	  c->x5, c->x4y, c->x3y2, c->x2y3, c->xy4));
    Msgf(("6-pole %g %g %g %g %g %g %g...\n", 
	  c->x6, c->x5y, c->x4y2, c->x3y3, c->x2y4, c->xy5, c->x2y2z2));
    Msgf(("7-pole %g %g %g %g %g %g %g ...\n", 
	  c->x7, c->x6y, c->x5y2, c->x4y3, c->x3y4, c->x2y5, c->xy6));
    Msgf(("8-pole %g %g %g %g %g %g %g %g %g ...\n", 
	  c->x8, c->x7y, c->x6y2, c->x5y3, c->x4y4, c->x3y5, c->x2y6, c->xy7, c->x4y2z2));
    
    if (subtract_background) Q[0] = 0;
    else Q[0] = 1;
    Q[1] = -3*c->x2y2 - 3*c->x2z2 + c->x4 - 3*c->y2z2 + c->y4 + c->z4;
    Q[2] = -4*c->x3 + 6*c->xy2 + 6*c->xz2;
    Q[3] = 6*c->x2 - 3*c->y2 - 3*c->z2;
    Q[4] = -4*c->x;
    Q[5] = 6*c->x2y - 4*c->y3 + 6*c->yz2;
    Q[6] = -12*c->xy;
    Q[7] = 6*c->y;
    Q[8] = -3*c->x2 + 6*c->y2 - 3*c->z2;
    Q[9] = 6*c->x2z + 6*c->y2z - 4*c->z3;
    Q[10] = -12*c->xz;
    Q[11] = 6*c->z;
    Q[12] = -12*c->yz;
    Q[13] = -3*c->x2 - 3*c->y2 + 6*c->z2;
    Q[14] = 90*c->x2y2z2 - (15*c->x2y4)/2. - (15*c->x2z4)/2. - (15*c->x4y2)/2. - (15*c->x4z2)/2. + c->x6 - (15*c->y2z4)/2. - (15*c->y4z2)/2. + c->y6 + c->z6;
    Q[15] = 30*c->x3y2 + 30*c->x3z2 - 6*c->x5 - 180*c->xy2z2 + 15*c->xy4 + 15*c->xz4;
    Q[16] = -45*c->x2y2 - 45*c->x2z2 + 15*c->x4 + 90*c->y2z2 - (15*c->y4)/2. - (15*c->z4)/2.;
    Q[17] = 30*c->x2y3 - 180*c->x2yz2 + 15*c->x4y + 30*c->y3z2 - 6*c->y5 + 15*c->yz4;
    Q[18] = -60*c->x3y - 60*c->xy3 + 360*c->xyz2;
    Q[19] = 90*c->x2y + 30*c->y3 - 180*c->yz2;
    Q[20] = -45*c->x2y2 + 90*c->x2z2 - (15*c->x4)/2. - 45*c->y2z2 + 15*c->y4 - (15*c->z4)/2.;
    Q[21] = 30*c->x3 + 90*c->xy2 - 180*c->xz2;
    Q[22] = -180*c->x2y2z + 30*c->x2z3 + 15*c->x4z + 30*c->y2z3 + 15*c->y4z - 6*c->z5;
    Q[23] = -60*c->x3z + 360*c->xy2z - 60*c->xz3;
    Q[24] = 90*c->x2z - 180*c->y2z + 30*c->z3;
    Q[25] = 360*c->x2yz - 60*c->y3z - 60*c->yz3;
    Q[26] = -720*c->xyz;
    Q[27] = -180*c->x2z + 90*c->y2z + 30*c->z3;
    Q[28] = 90*c->x2y2 - 45*c->x2z2 - (15*c->x4)/2. - 45*c->y2z2 - (15*c->y4)/2. + 15*c->z4;
    Q[29] = 30*c->x3 - 180*c->xy2 + 90*c->xz2;
    Q[30] = -180*c->x2y + 30*c->y3 + 90*c->yz2;
    Q[31] = -14*c->x2y6 - 14*c->x2z6 + 35*c->x4y4 + 35*c->x4z4 - 14*c->x6y2 - 14*c->x6z2 + c->x8 - 14*c->y2z6 + 35*c->y4z4 - 14*c->y6z2 + c->y8 + c->z8;
    Q[32] = -140*c->x3y4 - 140*c->x3z4 + 84*c->x5y2 + 84*c->x5z2 - 8*c->x7 + 28*c->xy6 + 28*c->xz6;
    Q[33] = 210*c->x2y4 + 210*c->x2z4 - 210*c->x4y2 - 210*c->x4z2 + 28*c->x6 - 14*c->y6 - 14*c->z6;
    Q[34] = 280*c->x3y2 + 280*c->x3z2 - 56*c->x5 - 140*c->xy4 - 140*c->xz4;
    Q[35] = -210*c->x2y2 - 210*c->x2z2 + 70*c->x4 + 35*c->y4 + 35*c->z4;
    Q[36] = 84*c->x2y5 - 140*c->x4y3 + 28*c->x6y - 140*c->y3z4 + 84*c->y5z2 - 8*c->y7 + 28*c->yz6;
    Q[37] = 560*c->x3y3 - 168*c->x5y - 168*c->xy5;
    Q[38] = -840*c->x2y3 + 420*c->x4y + 84*c->y5;
    Q[39] = -560*c->x3y + 560*c->xy3;
    Q[40] = 420*c->x2y - 140*c->y3;
    Q[41] = -210*c->x2y4 + 210*c->x4y2 - 14*c->x6 + 210*c->y2z4 - 210*c->y4z2 + 28*c->y6 - 14*c->z6;
    Q[42] = -840*c->x3y2 + 84*c->x5 + 420*c->xy4;
    Q[43] = 1260*c->x2y2 - 210*c->x4 - 210*c->y4;
    Q[44] = 280*c->x3 - 840*c->xy2;
    Q[45] = -210*c->x2 + 210*c->y2;
    Q[46] = 280*c->x2y3 - 140*c->x4y + 280*c->y3z2 - 56*c->y5 - 140*c->yz4;
    Q[47] = -210*c->x2y2 + 35*c->x4 - 210*c->y2z2 + 70*c->y4 + 35*c->z4;
    Q[48] = 84*c->x2z5 - 140*c->x4z3 + 28*c->x6z + 84*c->y2z5 - 140*c->y4z3 + 28*c->y6z - 8*c->z7;
    Q[49] = 560*c->x3z3 - 168*c->x5z - 168*c->xz5;
    Q[50] = -840*c->x2z3 + 420*c->x4z + 84*c->z5;
    Q[51] = -560*c->x3z + 560*c->xz3;
    Q[52] = 420*c->x2z - 140*c->z3;
    Q[53] = 560*c->y3z3 - 168*c->y5z - 168*c->yz5;
    Q[54] = -840*c->y2z3 + 420*c->y4z + 84*c->z5;
    Q[55] = -560*c->y3z + 560*c->yz3;
    Q[56] = 420*c->y2z - 140*c->z3;
    Q[57] = -210*c->x2z4 + 210*c->x4z2 - 14*c->x6 - 210*c->y2z4 + 210*c->y4z2 - 14*c->y6 + 28*c->z6;
    Q[58] = -840*c->x3z2 + 84*c->x5 + 420*c->xz4;
    Q[59] = 1260*c->x2z2 - 210*c->x4 - 210*c->z4;
    Q[60] = 280*c->x3 - 840*c->xz2;
    Q[61] = -210*c->x2 + 210*c->z2;
    Q[62] = -840*c->y3z2 + 84*c->y5 + 420*c->yz4;
    Q[63] = 1260*c->y2z2 - 210*c->y4 - 210*c->z4;
    Q[64] = 280*c->y3 - 840*c->yz2;
    Q[65] = -210*c->y2 + 210*c->z2;
    Q[66] = 280*c->x2z3 - 140*c->x4z + 280*c->y2z3 - 140*c->y4z - 56*c->z5;
    Q[67] = -210*c->x2z2 + 35*c->x4 - 210*c->y2z2 + 35*c->y4 + 70*c->z4;
}

void
cubic_acc(const float *r, float a, float *acc)
{
    float mx = -a-r[0];
    float my = -a-r[1];
    float mz = -a-r[2];
    float px =  a-r[0];
    float py =  a-r[1];
    float pz =  a-r[2];

    float mx2 = mx*mx;
    float my2 = my*my;
    float mz2 = mz*mz;
    float px2 = px*px;
    float py2 = py*py;
    float pz2 = pz*pz;

    float mmm = mx2+my2+mz2;
    float mmp = mx2+my2+pz2;
    float mpm = mx2+py2+mz2;
    float pmm = px2+my2+mz2;
    float mpp = mx2+py2+pz2;
    float ppm = px2+py2+mz2;
    float pmp = px2+my2+pz2;
    float ppp = px2+py2+pz2;

    float smmm = sqrtf(mx2+my2+mz2);
    float smmp = sqrtf(mx2+my2+pz2);
    float smpm = sqrtf(mx2+py2+mz2);
    float spmm = sqrtf(px2+my2+mz2);
    float smpp = sqrtf(mx2+py2+pz2);
    float sppm = sqrtf(px2+py2+mz2);
    float spmp = sqrtf(px2+my2+pz2);
    float sppp = sqrtf(px2+py2+pz2);

    float smmm3 = smmm*smmm*smmm;
    float smmp3 = smmp*smmp*smmp;
    float smpm3 = smpm*smpm*smpm;
    float spmm3 = spmm*spmm*spmm;
    float smpp3 = smpp*smpp*smpp;
    float sppm3 = sppm*sppm*sppm;
    float spmp3 = spmp*spmp*spmp;
    float sppp3 = sppp*sppp*sppp;

    acc[0] = -(((-1 - mx/smmm)*my*mz)/(smmm + mx)) + (mx2*mz)/(smmm*(smmm + my)) + (mx2*my)/(smmm*(smmm + mz)) - (mz*px2)/((my + spmm)*spmm) - (my*px2)/((mz + spmm)*spmm) + (my*mz*(-1 - px/spmm))/(spmm + px) + ((-1 - mx/smpm)*mz*py)/(smpm + mx) - (mx2*py)/(smpm*(smpm + mz)) + (px2*py)/((mz + sppm)*sppm) - (mz*(-1 - px/sppm)*py)/(sppm + px) - (mx2*mz)/(smpm*(smpm + py)) + (mz*px2)/(sppm*(sppm + py)) + ((-1 - mx/smmp)*my*pz)/(smmp + mx) - (mx2*pz)/(smmp*(smmp + my)) + (px2*pz)/((my + spmp)*spmp) - (my*(-1 - px/spmp)*pz)/(spmp + px) - ((-1 - mx/smpp)*py*pz)/(smpp + mx) + ((-1 - px/sppp)*py*pz)/(sppp + px) + (mx2*pz)/(smpp*(smpp + py)) - (px2*pz)/(sppp*(sppp + py)) - (mx2*my)/(smmp*(smmp + pz)) + (mx2*py)/(smpp*(smpp + pz)) + (my*px2)/(spmp*(spmp + pz)) - (px2*py)/(sppp*(sppp + pz)) + (((-(my/(smmm*mz)) + (mx2*my)/(smmm3*mz))*mz2)/(1 + (mx2*my2)/(mmm*mz2)) + (my2*(-(mz/(smmm*my)) + (mx2*mz)/(smmm3*my)))/(1 + (mx2*mz2)/(mmm*my2)) + (mx2*((my*mz)/smmm3 + (my*mz)/(smmm*mx2)))/(1 + (my2*mz2)/(mmm*mx2)) - (((my*mz)/spmm3 + (my*mz)/(spmm*px2))*px2)/(1 + (my2*mz2)/(pmm*px2)) - (mz2*(-(my/(mz*spmm)) + (my*px2)/(mz*spmm3)))/(1 + (my2*px2)/(mz2*pmm)) - (my2*(-(mz/(my*spmm)) + (mz*px2)/(my*spmm3)))/(1 + (mz2*px2)/(my2*pmm)) - ((-(mz/(smpm*py)) + (mx2*mz)/(smpm3*py))*py2)/(1 + (mx2*mz2)/(mpm*py2)) + ((-(mz/(sppm*py)) + (mz*px2)/(sppm3*py))*py2)/(1 + (mz2*px2)/(ppm*py2)) - (mz2*(-(py/(smpm*mz)) + (mx2*py)/(smpm3*mz)))/(1 + (mx2*py2)/(mpm*mz2)) - (mx2*((mz*py)/smpm3 + (mz*py)/(smpm*mx2)))/(1 + (mz2*py2)/(mpm*mx2)) + (px2*((mz*py)/sppm3 + (mz*py)/(sppm*px2)))/(1 + (mz2*py2)/(ppm*px2)) + (mz2*(-(py/(mz*sppm)) + (px2*py)/(mz*sppm3)))/(1 + (px2*py2)/(mz2*ppm)) - ((-(my/(smmp*pz)) + (mx2*my)/(smmp3*pz))*pz2)/(1 + (mx2*my2)/(mmp*pz2)) + ((-(my/(spmp*pz)) + (my*px2)/(spmp3*pz))*pz2)/(1 + (my2*px2)/(pmp*pz2)) + ((-(py/(smpp*pz)) + (mx2*py)/(smpp3*pz))*pz2)/(1 + (mx2*py2)/(mpp*pz2)) - ((-(py/(sppp*pz)) + (px2*py)/(sppp3*pz))*pz2)/(1 + (px2*py2)/(ppp*pz2)) - (my2*(-(pz/(smmp*my)) + (mx2*pz)/(smmp3*my)))/(1 + (mx2*pz2)/(mmp*my2)) - (mx2*((my*pz)/smmp3 + (my*pz)/(smmp*mx2)))/(1 + (my2*pz2)/(mmp*mx2)) + (px2*((my*pz)/spmp3 + (my*pz)/(spmp*px2)))/(1 + (my2*pz2)/(pmp*px2)) + (my2*(-(pz/(my*spmp)) + (px2*pz)/(my*spmp3)))/(1 + (px2*pz2)/(my2*pmp)) + (py2*(-(pz/(smpp*py)) + (mx2*pz)/(smpp3*py)))/(1 + (mx2*pz2)/(mpp*py2)) - (py2*(-(pz/(sppp*py)) + (px2*pz)/(sppp3*py)))/(1 + (px2*pz2)/(ppp*py2)) + (mx2*((py*pz)/smpp3 + (py*pz)/(smpp*mx2)))/(1 + (py2*pz2)/(mpp*mx2)) - (px2*((py*pz)/sppp3 + (py*pz)/(sppp*px2)))/(1 + (py2*pz2)/(ppp*px2)) - 2*mx*atanf((my*mz)/(smmm*mx)) + 2*px*atanf((my*mz)/(spmm*px)) + 2*mx*atanf((mz*py)/(smpm*mx)) - 2*px*atanf((mz*py)/(sppm*px)) + 2*mx*atanf((my*pz)/(smmp*mx)) - 2*px*atanf((my*pz)/(spmp*px)) - 2*mx*atanf((py*pz)/(smpp*mx)) + 2*px*atanf((py*pz)/(sppp*px)))/2. + mz*logf(smmm + my) - pz*logf(smmp + my) + my*logf(smmm + mz) - py*logf(smpm + mz) - mz*logf(my + spmm) - my*logf(mz + spmm) + pz*logf(my + spmp) + py*logf(mz + sppm) - mz*logf(smpm + py) + pz*logf(smpp + py) + mz*logf(sppm + py) - pz*logf(sppp + py) - my*logf(smmp + pz) + py*logf(smpp + pz) + my*logf(spmp + pz) - py*logf(sppp + pz);

  acc[1] = (my2*mz)/(smmm*(smmm + mx)) - (mx*(-1 - my/smmm)*mz)/(smmm + my) + (mx*my2)/(smmm*(smmm + mz)) + (mz*(-1 - my/spmm)*px)/(my + spmm) - (my2*px)/((mz + spmm)*spmm) - (my2*mz)/(spmm*(spmm + px)) - (mz*py2)/(smpm*(smpm + mx)) - (mx*py2)/(smpm*(smpm + mz)) + (px*py2)/((mz + sppm)*sppm) + (mz*py2)/(sppm*(sppm + px)) + (mx*mz*(-1 - py/smpm))/(smpm + py) - (mz*px*(-1 - py/sppm))/(sppm + py) - (my2*pz)/(smmp*(smmp + mx)) + (mx*(-1 - my/smmp)*pz)/(smmp + my) - ((-1 - my/spmp)*px*pz)/(my + spmp) + (my2*pz)/(spmp*(spmp + px)) + (py2*pz)/(smpp*(smpp + mx)) - (py2*pz)/(sppp*(sppp + px)) - (mx*(-1 - py/smpp)*pz)/(smpp + py) + (px*(-1 - py/sppp)*pz)/(sppp + py) - (mx*my2)/(smmp*(smmp + pz)) + (mx*py2)/(smpp*(smpp + pz)) + (my2*px)/(spmp*(spmp + pz)) - (px*py2)/(sppp*(sppp + pz)) + (((-(mx/(smmm*mz)) + (mx*my2)/(smmm3*mz))*mz2)/(1 + (mx2*my2)/(mmm*mz2)) + (my2*((mx*mz)/smmm3 + (mx*mz)/(smmm*my2)))/(1 + (mx2*mz2)/(mmm*my2)) + (mx2*(-(mz/(smmm*mx)) + (my2*mz)/(smmm3*mx)))/(1 + (my2*mz2)/(mmm*mx2)) - (((my2*mz)/(spmm3*px) - mz/(spmm*px))*px2)/(1 + (my2*mz2)/(pmm*px2)) - (mz2*((my2*px)/(mz*spmm3) - px/(mz*spmm)))/(1 + (my2*px2)/(mz2*pmm)) - (my2*((mz*px)/spmm3 + (mz*px)/(my2*spmm)))/(1 + (mz2*px2)/(my2*pmm)) - (((mx*mz)/smpm3 + (mx*mz)/(smpm*py2))*py2)/(1 + (mx2*mz2)/(mpm*py2)) + (((mz*px)/sppm3 + (mz*px)/(sppm*py2))*py2)/(1 + (mz2*px2)/(ppm*py2)) - (mz2*(-(mx/(smpm*mz)) + (mx*py2)/(smpm3*mz)))/(1 + (mx2*py2)/(mpm*mz2)) - (mx2*(-(mz/(smpm*mx)) + (mz*py2)/(smpm3*mx)))/(1 + (mz2*py2)/(mpm*mx2)) + (px2*(-(mz/(sppm*px)) + (mz*py2)/(sppm3*px)))/(1 + (mz2*py2)/(ppm*px2)) + (mz2*(-(px/(mz*sppm)) + (px*py2)/(mz*sppm3)))/(1 + (px2*py2)/(mz2*ppm)) - ((-(mx/(smmp*pz)) + (mx*my2)/(smmp3*pz))*pz2)/(1 + (mx2*my2)/(mmp*pz2)) + (((my2*px)/(spmp3*pz) - px/(spmp*pz))*pz2)/(1 + (my2*px2)/(pmp*pz2)) + ((-(mx/(smpp*pz)) + (mx*py2)/(smpp3*pz))*pz2)/(1 + (mx2*py2)/(mpp*pz2)) - ((-(px/(sppp*pz)) + (px*py2)/(sppp3*pz))*pz2)/(1 + (px2*py2)/(ppp*pz2)) - (my2*((mx*pz)/smmp3 + (mx*pz)/(smmp*my2)))/(1 + (mx2*pz2)/(mmp*my2)) - (mx2*(-(pz/(smmp*mx)) + (my2*pz)/(smmp3*mx)))/(1 + (my2*pz2)/(mmp*mx2)) + (px2*((my2*pz)/(spmp3*px) - pz/(spmp*px)))/(1 + (my2*pz2)/(pmp*px2)) + (my2*((px*pz)/spmp3 + (px*pz)/(my2*spmp)))/(1 + (px2*pz2)/(my2*pmp)) + (py2*((mx*pz)/smpp3 + (mx*pz)/(smpp*py2)))/(1 + (mx2*pz2)/(mpp*py2)) - (py2*((px*pz)/sppp3 + (px*pz)/(sppp*py2)))/(1 + (px2*pz2)/(ppp*py2)) + (mx2*(-(pz/(smpp*mx)) + (py2*pz)/(smpp3*mx)))/(1 + (py2*pz2)/(mpp*mx2)) - (px2*(-(pz/(sppp*px)) + (py2*pz)/(sppp3*px)))/(1 + (py2*pz2)/(ppp*px2)) - 2*my*atanf((mx*mz)/(smmm*my)) + 2*my*atanf((mz*px)/(my*spmm)) + 2*py*atanf((mx*mz)/(smpm*py)) - 2*py*atanf((mz*px)/(sppm*py)) + 2*my*atanf((mx*pz)/(smmp*my)) - 2*my*atanf((px*pz)/(my*spmp)) - 2*py*atanf((mx*pz)/(smpp*py)) + 2*py*atanf((px*pz)/(sppp*py)))/2. + mz*logf(smmm + mx) - pz*logf(smmp + mx) - mz*logf(smpm + mx) + pz*logf(smpp + mx) + mx*logf(smmm + mz) - mx*logf(smpm + mz) - px*logf(mz + spmm) + px*logf(mz + sppm) - mz*logf(spmm + px) + pz*logf(spmp + px) + mz*logf(sppm + px) - pz*logf(sppp + px) - mx*logf(smmp + pz) + mx*logf(smpp + pz) + px*logf(spmp + pz) - px*logf(sppp + pz);

  acc[2] = (my*mz2)/(smmm*(smmm + mx)) + (mx*mz2)/(smmm*(smmm + my)) - (mx*my*(-1 - mz/smmm))/(smmm + mz) + (my*(-1 - mz/spmm)*px)/(mz + spmm) - (mz2*px)/((my + spmm)*spmm) - (my*mz2)/(spmm*(spmm + px)) - (mz2*py)/(smpm*(smpm + mx)) + (mx*(-1 - mz/smpm)*py)/(smpm + mz) - ((-1 - mz/sppm)*px*py)/(mz + sppm) + (mz2*py)/(sppm*(sppm + px)) - (mx*mz2)/(smpm*(smpm + py)) + (mz2*px)/(sppm*(sppm + py)) - (my*pz2)/(smmp*(smmp + mx)) - (mx*pz2)/(smmp*(smmp + my)) + (px*pz2)/((my + spmp)*spmp) + (my*pz2)/(spmp*(spmp + px)) + (py*pz2)/(smpp*(smpp + mx)) - (py*pz2)/(sppp*(sppp + px)) + (mx*pz2)/(smpp*(smpp + py)) - (px*pz2)/(sppp*(sppp + py)) + (mx*my*(-1 - pz/smmp))/(smmp + pz) - (mx*py*(-1 - pz/smpp))/(smpp + pz) - (my*px*(-1 - pz/spmp))/(spmp + pz) + (px*py*(-1 - pz/sppp))/(sppp + pz) + ((((mx*my)/smmm3 + (mx*my)/(smmm*mz2))*mz2)/(1 + (mx2*my2)/(mmm*mz2)) + (my2*(-(mx/(smmm*my)) + (mx*mz2)/(smmm3*my)))/(1 + (mx2*mz2)/(mmm*my2)) + (mx2*(-(my/(smmm*mx)) + (my*mz2)/(smmm3*mx)))/(1 + (my2*mz2)/(mmm*mx2)) - (((my*mz2)/(spmm3*px) - my/(spmm*px))*px2)/(1 + (my2*mz2)/(pmm*px2)) - (mz2*((my*px)/spmm3 + (my*px)/(mz2*spmm)))/(1 + (my2*px2)/(mz2*pmm)) - (my2*((mz2*px)/(my*spmm3) - px/(my*spmm)))/(1 + (mz2*px2)/(my2*pmm)) - ((-(mx/(smpm*py)) + (mx*mz2)/(smpm3*py))*py2)/(1 + (mx2*mz2)/(mpm*py2)) + (((mz2*px)/(sppm3*py) - px/(sppm*py))*py2)/(1 + (mz2*px2)/(ppm*py2)) - (mz2*((mx*py)/smpm3 + (mx*py)/(smpm*mz2)))/(1 + (mx2*py2)/(mpm*mz2)) - (mx2*(-(py/(smpm*mx)) + (mz2*py)/(smpm3*mx)))/(1 + (mz2*py2)/(mpm*mx2)) + (px2*((mz2*py)/(sppm3*px) - py/(sppm*px)))/(1 + (mz2*py2)/(ppm*px2)) + (mz2*((px*py)/sppm3 + (px*py)/(mz2*sppm)))/(1 + (px2*py2)/(mz2*ppm)) - (((mx*my)/smmp3 + (mx*my)/(smmp*pz2))*pz2)/(1 + (mx2*my2)/(mmp*pz2)) + (((my*px)/spmp3 + (my*px)/(spmp*pz2))*pz2)/(1 + (my2*px2)/(pmp*pz2)) + (((mx*py)/smpp3 + (mx*py)/(smpp*pz2))*pz2)/(1 + (mx2*py2)/(mpp*pz2)) - (((px*py)/sppp3 + (px*py)/(sppp*pz2))*pz2)/(1 + (px2*py2)/(ppp*pz2)) - (my2*(-(mx/(smmp*my)) + (mx*pz2)/(smmp3*my)))/(1 + (mx2*pz2)/(mmp*my2)) - (mx2*(-(my/(smmp*mx)) + (my*pz2)/(smmp3*mx)))/(1 + (my2*pz2)/(mmp*mx2)) + (px2*(-(my/(spmp*px)) + (my*pz2)/(spmp3*px)))/(1 + (my2*pz2)/(pmp*px2)) + (my2*(-(px/(my*spmp)) + (px*pz2)/(my*spmp3)))/(1 + (px2*pz2)/(my2*pmp)) + (py2*(-(mx/(smpp*py)) + (mx*pz2)/(smpp3*py)))/(1 + (mx2*pz2)/(mpp*py2)) - (py2*(-(px/(sppp*py)) + (px*pz2)/(sppp3*py)))/(1 + (px2*pz2)/(ppp*py2)) + (mx2*(-(py/(smpp*mx)) + (py*pz2)/(smpp3*mx)))/(1 + (py2*pz2)/(mpp*mx2)) - (px2*(-(py/(sppp*px)) + (py*pz2)/(sppp3*px)))/(1 + (py2*pz2)/(ppp*px2)) - 2*mz*atanf((mx*my)/(smmm*mz)) + 2*mz*atanf((my*px)/(mz*spmm)) + 2*mz*atanf((mx*py)/(smpm*mz)) - 2*mz*atanf((px*py)/(mz*sppm)) + 2*pz*atanf((mx*my)/(smmp*pz)) - 2*pz*atanf((my*px)/(spmp*pz)) - 2*pz*atanf((mx*py)/(smpp*pz)) + 2*pz*atanf((px*py)/(sppp*pz)))/2. + my*logf(smmm + mx) - my*logf(smmp + mx) - py*logf(smpm + mx) + py*logf(smpp + mx) + mx*logf(smmm + my) - mx*logf(smmp + my) - px*logf(my + spmm) + px*logf(my + spmp) - my*logf(spmm + px) + my*logf(spmp + px) + py*logf(sppm + px) - py*logf(sppp + px) - mx*logf(smpm + py) + mx*logf(smpp + py) + px*logf(sppm + py) - px*logf(sppp + py);
}
