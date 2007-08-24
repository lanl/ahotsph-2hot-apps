#define NO_MSGS
#define NOTIMERS  /* Timers are a major performance hit on the delta */
#include "physics.h"
#include "vop.h"
#include "tensop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"
#include "stk.h"


void 
do_grav(const float *p, const float *end, const float *pos0, float *mass0, 
	float *acc0, float *phi0, const float *eps2p, int *ncut)
{
    float dr2;
    Vxd(float r);
    float phii, mor3, mass;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;

    VxV(a, = acc0);

    while (p < end) {
	mass = *p++;
	r0 = *p++;
	r1 = *p++;
#if NDIM > 2
	r2 = *p++;
#endif
	VxVx(r, -= ppos);	/* 3 flops */

	dr2 = Dotx(r, r);	/* 5 flops */
	dr2 += eps2;

	phii = recipsqrtf(dr2);	/* 8 flops */
	
	mor3 = phii * phii;	/* 5 flops */
	phii *= mass;
	total_mass += mass;
	mor3 *= phii;
	phi -= phii;

	VxVx(a, += mor3 * r);	/* 6 flops */
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}

void
update_point_mass(body *btab, int nobj, 
		  body *p, float smooth2, float newt)
{
    body *r;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxV(ppos, = p->pos);

    for (r = btab; r < btab+nobj; r++) {
	VxVVx(r, = r->pos, - ppos); /* 3 flops */
	
	dr2 = Dotx(r, r);	/* 5 flops */
	
	if (dr2 != (float)0.0) {
	
	  dr2 += smooth2;
	
	  oneor = recipsqrtf(dr2);	/* 8 flops */
	
	  oneor2 = oneor * oneor;	/* 17 flops */
	  phii = newt * oneor * r->mass;
	  p->phi -= phii;
	  VVx(p->acc, += oneor2 * phii * r);
	}
    }
}

void 
do_grav_u2(const float *p, const float *end, const float *pos0, float *mass0,
	   float *acc0,	float *phi0, const float *eps2p, int *ncut)
{
    float dr2a, dr2b;
    Vxd(float ra);
    Vxd(float rb);
    float phiia, phiib;
    float mor3a, mor3b;
    float massa, massb;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;

    VxV(a, = acc0);

    while (p < end) {
	massa = *p++;
	ra0 = *p++;
	ra1 = *p++;
#if NDIM > 2
	ra2 = *p++;
#endif

	massb = *p++;
	rb0 = *p++;
	rb1 = *p++;
#if NDIM > 2
	rb2 = *p++;
#endif

	VxVx(ra, -= ppos);
	VxVx(rb, -= ppos);

	dr2a = Dotx(ra, ra);
	dr2b = Dotx(rb, rb);

	dr2a += eps2;
	dr2b += eps2;

	phiia = recipsqrtf(dr2a);
	phiib = recipsqrtf(dr2b);
	
	mor3a = phiia * phiia;
	mor3b = phiib * phiib;
	phiia *= massa;
	phiib *= massb;
	total_mass += massa + massb;
	mor3a *= phiia;
	mor3b *= phiib;
	phi -= phiia;
	phi -= phiib;

	VxVx(a, += mor3a * ra);
	VxVx(a, += mor3b * rb);
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}

#if defined(__T3D__) || defined (_IBMR2)
#define USE_CHEB_RSQRT
#endif

#ifdef USE_CHEB_RSQRT
#include "karp.h"

/* This does 38 actual flops per interaction */
/* The approximate rsqrt uses 13 multiplies and 6 adds */

void 
do_cheb_u2(const float *p, const float *end, const float *pos0, float *mass0,
	   float *acc0,	float *phi0, const float *eps2p, int *ncut)
{
    float dr2a, dr2b;
    Vxd(float ra);
    Vxd(float rb);
    float mor3a, mor3b;
    float massa, massb;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;
    float xa, xb;
    unsigned int ita, itb;
    const int mbits = 23;

    VxV(a, = acc0);

    while (p < end) {
	massa = *p++;
	ra0 = *p++;
	ra1 = *p++;
	ra2 = *p++;

	massb = *p++;
	rb0 = *p++;
	rb1 = *p++;
	rb2 = *p++;

	VxVx(ra, -= ppos);
	VxVx(rb, -= ppos);

	dr2a = Dotx(ra, ra);
	dr2b = Dotx(rb, rb);
	dr2a += eps2;
	dr2b += eps2;

	ita = (*(FLOAT_PUN *)&dr2a) >> mbits;
	itb = (*(FLOAT_PUN *)&dr2b) >> mbits; 
	xa = dr2a*u[ita];
	xb = dr2b*u[itb];
	mor3a = g0 + xa*(g1 + xa*(g2 + xa*(g3 + xa*(g4 + xa*g5))));
	mor3b = g0 + xb*(g1 + xb*(g2 + xb*(g3 + xb*(g4 + xb*g5))));
	mor3a *= (float)1.5 - (float)0.5*xa*xa*xa*mor3a*mor3a;
	mor3b *= (float)1.5 - (float)0.5*xb*xb*xb*mor3b*mor3b;
	mor3a *= t[ita];
	mor3b *= t[itb];

	mor3a *= massa;
	mor3b *= massb;
	total_mass += massa + massb;
	phi -= dr2a*mor3a;
	phi -= dr2b*mor3b;

	VxVx(a, += mor3a * ra);
	VxVx(a, += mor3b * rb);
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}
#endif

#ifdef __alpha

#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38F
#endif

/* Table to reduce domain of input for Chebychev polynomial approx.
   The first column is u = (1/2)^k.  The second column is t=u^3/2 */

static const float ut[2*256] = {
  1.70141e+38,         FLT_MAX,  
8.50706e+37,	     FLT_MAX,    
4.25353e+37,	     FLT_MAX,    
2.12676e+37,	     FLT_MAX,    
1.06338e+37,	     FLT_MAX,    
5.31691e+36,	     FLT_MAX,    
  2.65846e+36,	      FLT_MAX,   
1.32923e+36,	     FLT_MAX,    
6.64614e+35,	     FLT_MAX,    
3.32307e+35,	     FLT_MAX,    
1.66153e+35,	     FLT_MAX,    
8.30767e+34,	     FLT_MAX,    
  4.15384e+34,	       FLT_MAX,  
2.07692e+34,	     FLT_MAX,    
1.03846e+34,	     FLT_MAX,    
5.19230e+33,	     FLT_MAX,    
2.59615e+33,	     FLT_MAX,    
1.29807e+33,	     FLT_MAX,    
  6.49037e+32,	       FLT_MAX,  
3.24519e+32,	     FLT_MAX,    
1.62259e+32,	     FLT_MAX,    
8.11296e+31,	     FLT_MAX,    
4.05648e+31,	     FLT_MAX,    
2.02824e+31,	     FLT_MAX,    
  1.01412e+31,	       FLT_MAX,  
5.07060e+30,	     FLT_MAX,    
2.53530e+30,	     FLT_MAX,    
1.26765e+30,	     FLT_MAX,    
6.33825e+29,	     FLT_MAX,    
3.16913e+29,	     FLT_MAX,    
  1.58456e+29,	       FLT_MAX,  
7.92282e+28,	     FLT_MAX,    
3.96141e+28,	     FLT_MAX,    
1.98070e+28,	     FLT_MAX,    
9.90352e+27,	     FLT_MAX,    
4.95176e+27,	     FLT_MAX,    
  2.47588e+27,	       FLT_MAX,  
1.23794e+27,	     FLT_MAX,    
6.18970e+26,	     FLT_MAX,    
3.09485e+26,	     FLT_MAX,    
1.54743e+26,	     FLT_MAX,    
7.73713e+25,	     FLT_MAX,    
  3.86856e+25,	   2.40616e+38,  
1.93428e+25,	 8.50706e+37,    
9.67141e+24,	 3.00770e+37,    
4.83570e+24,	 1.06338e+37,    
2.41785e+24,	 3.75962e+36,    
1.20893e+24,	 1.32923e+36,    
  6.04463e+23,	   4.69953e+35,  
3.02231e+23,	 1.66153e+35,    
1.51116e+23,	 5.87441e+34,    
7.55579e+22,	 2.07692e+34,    
3.77789e+22,	 7.34302e+33,    
1.88895e+22,	 2.59615e+33,    
  9.44473e+21,	   9.17877e+32,  
4.72237e+21,	 3.24519e+32,    
2.36118e+21,	 1.14735e+32,    
1.18059e+21,	 4.05648e+31,    
5.90296e+20,	 1.43418e+31,    
2.95148e+20,	 5.07060e+30,    
  1.47574e+20,	   1.79273e+30,  
7.37870e+19,	 6.33825e+29,    
3.68935e+19,	 2.24091e+29,    
1.84467e+19,	 7.92282e+28,    
9.22337e+18,	 2.80114e+28,    
4.61169e+18,	 9.90352e+27,    
  2.30584e+18,	   3.50142e+27,  
1.15292e+18,	 1.23794e+27,    
5.76461e+17,	 4.37678e+26,    
2.88230e+17,	 1.54743e+26,    
1.44115e+17,	 5.47097e+25,    
7.20576e+16,	 1.93428e+25,    
  3.60288e+16,	   6.83872e+24,  
1.80144e+16,	 2.41785e+24,    
9.00720e+15,	 8.54840e+23,    
4.50360e+15,	 3.02231e+23,    
2.25180e+15,	 1.06855e+23,    
1.12590e+15,	 3.77789e+22,    
  5.62950e+14,	   1.33569e+22,  
2.81475e+14,	 4.72237e+21,    
1.40737e+14,	 1.66961e+21,    
7.03687e+13,	 5.90296e+20,    
3.51844e+13,	 2.08701e+20,    
1.75922e+13,	 7.37870e+19,    
  8.79609e+12,	   2.60876e+19,  
4.39805e+12,	 9.22337e+18,    
2.19902e+12,	 3.26095e+18,    
1.09951e+12,	 1.15292e+18,    
5.49756e+11,	 4.07619e+17,    
2.74878e+11,	 1.44115e+17,    
  1.37439e+11,	   5.09524e+16,  
6.87195e+10,	 1.80144e+16,    
3.43597e+10,	 6.36905e+15,    
1.71799e+10,	 2.25180e+15,    
8.58993e+09,	 7.96131e+14,    
4.29497e+09,	 2.81475e+14,    
  2.14748e+09,	   9.95164e+13,  
1.07374e+09,	 3.51844e+13,    
5.36871e+08,	 1.24396e+13,    
2.68435e+08,	 4.39805e+12,    
1.34218e+08,	 1.55494e+12,    
6.71089e+07,	 5.49756e+11,    
  3.35544e+07,	   1.94368e+11,  
1.67772e+07,	 6.87195e+10,    
8.38861e+06,	 2.42960e+10,    
4.19430e+06,	 8.58993e+09,    
2.09715e+06,	 3.03700e+09,    
1.04858e+06,	 1.07374e+09,    
      524288.,	   3.79625e+08,  
    262144.,	 1.34218e+08,    
    131072.,	 4.74531e+07,    
    65536.0,	 1.67772e+07,    
    32768.0,	 5.93164e+06,    
    16384.0,	 2.09715e+06,    
      8192.00,	       741455.,  
    4096.00,	     262144.,    
    2048.00,	     92681.9,    
    1024.00,	     32768.0,    
    512.000,	     11585.2,    
    256.000,	     4096.00,    
      128.000,	       1448.15,  
    64.0000,	     512.000,    
    32.0000,	     181.019,    
    16.0000,	     64.0000,    
    8.00000,	     22.6274,    
    4.00000,	     8.00000,    
      2.00000,	       2.82843,  
    1.00000,	     1.00000,    
   0.500000,	    0.353553,    
   0.250000,	    0.125000,    
   0.125000,	   0.0441942,    
  0.0625000,	   0.0156250,    
    0.0312500,	    0.00552427,  
  0.0156250,	  0.00195312,    
 0.00781250,	 0.000690534,    
 0.00390625,	 0.000244141,    
 0.00195312,	 8.63167e-05,    
0.000976562,	 3.05176e-05,    
  0.000488281,	   1.07896e-05,  
0.000244141,	 3.81470e-06,    
0.000122070,	 1.34870e-06,    
6.10352e-05,	 4.76837e-07,    
3.05176e-05,	 1.68587e-07,    
1.52588e-05,	 5.96046e-08,    
  7.62939e-06,	   2.10734e-08,  
3.81470e-06,	 7.45058e-09,    
1.90735e-06,	 2.63418e-09,    
9.53674e-07,	 9.31323e-10,    
4.76837e-07,	 3.29272e-10,    
2.38419e-07,	 1.16415e-10,    
  1.19209e-07,	   4.11590e-11,  
5.96046e-08,	 1.45519e-11,    
2.98023e-08,	 5.14488e-12,    
1.49012e-08,	 1.81899e-12,    
7.45058e-09,	 6.43110e-13,    
3.72529e-09,	 2.27374e-13,    
  1.86265e-09,	   8.03887e-14,  
9.31323e-10,	 2.84217e-14,    
4.65661e-10,	 1.00486e-14,    
2.32831e-10,	 3.55271e-15,    
1.16415e-10,	 1.25607e-15,    
5.82077e-11,	 4.44089e-16,    
  2.91038e-11,	   1.57009e-16,  
1.45519e-11,	 5.55112e-17,    
7.27596e-12,	 1.96262e-17,    
3.63798e-12,	 6.93889e-18,    
1.81899e-12,	 2.45327e-18,    
9.09495e-13,	 8.67362e-19,    
  4.54747e-13,	   3.06659e-19,  
2.27374e-13,	 1.08420e-19,    
1.13687e-13,	 3.83323e-20,    
5.68434e-14,	 1.35525e-20,    
2.84217e-14,	 4.79154e-21,    
1.42109e-14,	 1.69407e-21,    
  7.10543e-15,	   5.98943e-22,  
3.55271e-15,	 2.11758e-22,    
1.77636e-15,	 7.48678e-23,    
8.88178e-16,	 2.64698e-23,    
4.44089e-16,	 9.35848e-24,    
2.22045e-16,	 3.30872e-24,    
  1.11022e-16,	   1.16981e-24,  
5.55112e-17,	 4.13590e-25,    
2.77556e-17,	 1.46226e-25,    
1.38778e-17,	 5.16988e-26,    
6.93889e-18,	 1.82783e-26,    
3.46945e-18,	 6.46235e-27,    
  1.73472e-18,	   2.28479e-27,  
8.67362e-19,	 8.07794e-28,    
4.33681e-19,	 2.85598e-28,    
2.16840e-19,	 1.00974e-28,    
1.08420e-19,	 3.56998e-29,    
5.42101e-20,	 1.26218e-29,    
  2.71051e-20,	   4.46247e-30,  
1.35525e-20,	 1.57772e-30,    
6.77626e-21,	 5.57809e-31,    
3.38813e-21,	 1.97215e-31,    
1.69407e-21,	 6.97261e-32,    
8.47033e-22,	 2.46519e-32,    
  4.23516e-22,	   8.71576e-33,  
2.11758e-22,	 3.08149e-33,    
1.05879e-22,	 1.08947e-33,    
5.29396e-23,	 3.85186e-34,    
2.64698e-23,	 1.36184e-34,    
1.32349e-23,	 4.81482e-35,    
  6.61744e-24,	   1.70230e-35,  
3.30872e-24,	 6.01853e-36,    
1.65436e-24,	 2.12787e-36,    
8.27181e-25,	 7.52316e-37,    
4.13590e-25,	 2.65984e-37,    
2.06795e-25,	 9.40395e-38,    
  0.0000000,     0.00000000};	/* let the compiler extend with 0. */

/* Chebychev coefficients of the monomials on [1,2] for minimax fit.  */

/* These should be in registers for the duration of do_grav_cheb.  If not,
   they should be copied into automatic variables */
static const float g0 =   7.05452470f,
		   g1 = -14.85088557f,
		   g2 =  14.70832310f,
		   g3 =  -7.83703555f,
		   g4 =   2.17094576f,
		   g5 =  -0.24598174f;
/* end of karp.h */

#define NELEM (NDIM+1)		/* mass and 3 positions */
#define LAST (NVEC*NELEM)	/* index of last float in btab */
#define LOOPCNT(nfloat, unroll)  (((nfloat)/NELEM + (unroll-1))/unroll)

void
do_grav_cheb_unroll3(const float *p, const float *end, const float *pos0, 
	       float *mass0, float *acc0, float *phi0, const float 
	       *eps2p, int *ncut)
{
    float dr2a;
    float dr2b;
    float dr2c;
    Vxd(float ra);
    Vxd(float rb);
    Vxd(float rc);
    float massa;
    float massb;
    float massc;
    float mor3a;
    float mor3b;
    float mor3c;
    Vxd(float ppos);
    Vxd(float aa);
    Vxd(float ab);
    Vxd(float ac);
    float phia = 0.F;
    float phib = 0.F;
    float phic = 0.F;
    float total_massa = 0.F;
    float total_massb = 0.F;
    float total_massc = 0.F;
    const float eps2 = *eps2p;
    const int n = LOOPCNT(end-p, 3);
    int i;
    float xa, xb, xc;
    float xa3, xb3, xc3;
    unsigned int ita, itb, itc;
    float ua = 0.0, ub = 0.0, uc = 0.0;
    float ta = 0.0, tb = 0.0, tc = 0.0;
    unsigned int italast = UINT_MAX;
    unsigned int itblast = UINT_MAX;
    unsigned int itclast = UINT_MAX;

    VxS(aa, = 0.F);
    VxS(ab, = 0.F);
    VxS(ac, = 0.F);
    VxV(ppos, = pos0);

    for (i = 0; i < n; i++) {
	massa = p[0];
	ra0 = p[1];
	ra1 = p[2];
	ra2 = p[3];

	massb = p[4];
	rb0 = p[5];
	rb1 = p[6];
	rb2 = p[7];

	massc = p[8];
	rc0 = p[9];
	rc1 = p[10];
	rc2 = p[11];
	p += 12;

	VxVx(ra, -= ppos);
	VxVx(rb, -= ppos);
	VxVx(rc, -= ppos);

	dr2a = eps2 + Dotx(ra, ra);
	dr2b = eps2 + Dotx(rb, rb);
	dr2c = eps2 + Dotx(rc, rc);

	ita = (*(int *)&dr2a) >> 23;  /* Pun float to int */
	if( ita != italast ){
	  ua = ut[2*ita];
	  ta = ut[2*ita+1];
	  italast = ita;
	}
	itb = (*(int *)&dr2b) >> 23;  /* Pun float to int */
	if( itb != itblast ){
	  ub = ut[2*itb];
	  tb = ut[2*itb+1];
	  itblast = itb;
	}
	itc = (*(int *)&dr2c) >> 23;  /* Pun float to int */
	if( itc != itclast ){
	  uc = ut[2*itc];
	  tc = ut[2*itc+1];
	  itclast = itc;
	}

	xa = dr2a*ua;  /* Scale input domain to [1,2] */
	xb = dr2b*ub;  /* Scale input domain to [1,2] */
	xc = dr2c*uc;  /* Scale input domain to [1,2] */

	xa3 = xa*xa*xa;
	xb3 = xb*xb*xb;
	xc3 = xc*xc*xc;

	mor3a = g0 + xa*(g1 + xa*(g2)) + xa3*(g3 + xa*(g4 + xa*g5));
	mor3b = g0 + xb*(g1 + xb*(g2)) + xb3*(g3 + xb*(g4 + xb*g5));
	mor3c = g0 + xc*(g1 + xc*(g2)) + xc3*(g3 + xc*(g4 + xc*g5));

	mor3a *= 1.5f - 0.5f*xa3*mor3a*mor3a;  /* N-R iteration */
	mor3a *= ta*massa;  /* Restore result scaling */

	mor3b *= 1.5f - 0.5f*xb3*mor3b*mor3b;  /* N-R iteration */
	mor3b *= tb*massb;  /* Restore result scaling */

	mor3c *= 1.5f - 0.5f*xc3*mor3c*mor3c;  /* N-R iteration */
	mor3c *= tc*massc;  /* Restore result scaling */
	
	total_massa += massa;
	total_massb += massb;
	total_massc += massc;

	phia -= dr2a*mor3a;
	phib -= dr2b*mor3b;
	phic -= dr2c*mor3c;

	VxVx(aa, += mor3a * ra);
	VxVx(ab, += mor3b * rb);
	VxVx(ac, += mor3c * rc);
    }
    VxVx(aa, += ab);
    VxVx(aa, += ac);
    VVx(acc0, += aa);
    *mass0 += total_massa + total_massb + total_massc;
    *phi0 += (phia + phib + phic);
}

#endif
