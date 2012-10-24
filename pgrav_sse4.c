/*
 * Copyright 2012 Michael S. Warren. All Rights Reserved.
 */
#ifndef __AVX__
typedef float v4sf __attribute__ ((vector_size (16)));

#include <math.h>
#include "error.h"

#define massp f[j+0]
#define xp f[j+1]
#define yp f[j+2]
#define zp f[j+3]
#define R  f[j+4]
#define qxx f[j+5]
#define qxy f[j+6]
#define qyy f[j+7]
#define qxz f[j+8]
#define qyz f[j+9]
#define qxxx f[j+10]
#define qxxy f[j+11]
#define qxyy f[j+12]
#define qyyy f[j+13]
#define qxxz f[j+14]
#define qxyz f[j+15]
#define qyyz f[j+16]
#define qxxxx f[j+17]
#define qxxxy f[j+18]
#define qxxyy f[j+19]
#define qxyyy f[j+20]
#define qyyyy f[j+21]
#define qxxxz f[j+22]
#define qxxyz f[j+23]
#define qxyyz f[j+24]
#define qyyyz f[j+25]

void
pHinteract(const float *p, float *accp, const int n, const int stride, 
	   const v4sf *f, const int source_n)
{
    int i, j;
    v4sf ppos0, ppos1, ppos2;
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0, a1, a2, phi, mass;
    const v4sf zero = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    /* 217 flops, 104 bytes per interaction*/
    for (i = 0; i < n*stride; i += stride) {
	a0 = a1 = a2 = phi = mass = zero;
	ppos0[0] = ppos0[1] = ppos0[2] = ppos0[3] = p[i+1];
	ppos1[0] = ppos1[1] = ppos1[2] = ppos1[3] = p[i+2];
	ppos2[0] = ppos2[1] = ppos2[2] = ppos2[3] = p[i+3];

	for (j = 0; j < source_n*26; j += 26) {
	    mass += massp;
	    x = ppos0 - xp;
	    y = ppos1 - yp;
	    z = ppos2 - zp;
	    
	    r2 = x*x + y*y + z*z;
	    
	    /* Prefetching improves uncached performance by 10-20% */
	    asm("prefetcht0 512(%rdi)");
	    rinv = __builtin_ia32_rsqrtps(r2);
	    
	    asm("prefetcht0 576(%rdi)");
	    /* Newton-Raphson */
	    t = rinv;
	    r2 *= rinv;
	    rinv *= r2;
	    rinv -= three;
	    rinv *= t;
	    rinv *= half;		/* flips sign to avoid storing -0.5 */
	    /* end Newton-Raphson */
	    
	    asm("prefetcht0 640(%rdi)");
	    t = rinv;
	    eqe = t*massp;
	    rinv2 = t*t;
	    t *= rinv2;
	    phi += eqe;
	    eqe *= rinv2;
	    a0 += x*eqe;
	    a1 += y*eqe;
	    a2 += z*eqe;
	    
	    asm("prefetcht0 704(%rdi)");
	    t *= three*rinv2*R*R;
	    eq0 = t*(qxx*x + qxy*y + qxz*z);
	    eq1 = t*(qyy*y + qxy*x + qyz*z);
	    eq2 = t*(-(qxx + qyy)*z + qxz*x + qyz*y);
	    eqe = half*(eq0*x + eq1*y + eq2*z);
	    phi += eqe;
	    eqe *= five*rinv2;
	    a0 += x*eqe - eq0;
	    a1 += y*eqe - eq1;
	    a2 += z*eqe - eq2;
	    
	    asm("prefetcht0 768(%rdi)");
	    t *= five*rinv2*R;
	    xx = half*x*x;
	    xy = x*y;
	    xz = x*z;
	    yy = half*y*y;
	    yz = y*z;
	    zz = half*z*z;
	    xxx = x*(third*xx - zz);
	    xxz = z*(xx - third*zz);
	    yyy = y*(third*yy - zz);
	    yyz = z*(yy - third*zz);
	    xx -= zz;
	    yy -= zz;
	    
	    asm("prefetcht0 832(%rdi)");
	    eq0 = t*(qxxx*xx + qxyy*yy + qxxy*xy + qxxz*xz + qxyz*yz);
	    eq1 = t*(qxyy*xy + qxxy*xx + qyyy*yy + qyyz*yz + qxyz*xz);
	    eq2 = t*(-(qxxx + qxyy)*xz - (qxxy + qyyy)*yz + qxxz*xx + qyyz*yy + qxyz*xy);
	    eqe = third*(eq0*x + eq1*y + eq2*z);
	    phi += eqe;
	    eqe *= seven*rinv2;
	    a0 += x*eqe - eq0;
	    a1 += y*eqe - eq1;
	    a2 += z*eqe - eq2;
	    
	    asm("prefetcht0 896(%rdi)");
	    t *= seven*rinv2*R;
	    xxy = y*xx;
	    xyy = x*yy;
	    xyz = xy*z;
	    eq0 = t*(qxxxx*xxx + qxyyy*yyy + qxxxy*xxy + qxxxz*xxz + qxxyy*xyy + qxxyz*xyz + qxyyz*yyz);
	    eq1 = t*(qxyyy*xyy + qxxxy*xxx + qyyyy*yyy + qyyyz*yyz + qxxyy*xxy + qxxyz*xxz + qxyyz*xyz);
	    eq2 = t*(-qxxxx*xxz - (qxyyy + qxxxy)*xyz - qyyyy*yyz + qxxxz*xxx + qyyyz*yyy  - qxxyy*(xxz + yyz) + qxxyz*xxy + qxyyz*xyy);
	    eqe = quarter*(eq0*x + eq1*y + eq2*z);
	    phi += eqe;
	    eqe *= nine*rinv2;
	    a0 += x*eqe - eq0;
	    a1 += y*eqe - eq1;
	    a2 += z*eqe - eq2;
	}
	a0 = __builtin_ia32_haddps(a0, a0);
	a1 = __builtin_ia32_haddps(a1, a1);
	a2 = __builtin_ia32_haddps(a2, a2);
	phi = __builtin_ia32_haddps(phi, phi);
	accp[i+0] += __builtin_ia32_haddps(a0, a0)[0];
	accp[i+1] += __builtin_ia32_haddps(a1, a1)[0];
	accp[i+2] += __builtin_ia32_haddps(a2, a2)[0];
	accp[i+3] += __builtin_ia32_haddps(phi, phi)[0];
    }
}

void
pQinteract(const float *p, float *accp, const int n, const int stride, 
	   const v4sf *f, const int source_n)
{
    int i, j;
    v4sf ppos0, ppos1, ppos2;
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0, a1, a2, phi, mass;
    const v4sf zero = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};

    /* 66 flops , 40 bytes per interactions */
    for (i = 0; i < n*stride; i += stride) {
	a0 = a1 = a2 = phi = mass = zero;
	ppos0[0] = ppos0[1] = ppos0[2] = ppos0[3] = p[i+1];
	ppos1[0] = ppos1[1] = ppos1[2] = ppos1[3] = p[i+2];
	ppos2[0] = ppos2[1] = ppos2[2] = ppos2[3] = p[i+3];
	
	for (j = 0; j < source_n*10; j += 10) {
	    mass += massp;
	    x = ppos0 - xp;
	    y = ppos1 - yp;
	    z = ppos2 - zp;

	    r2 = x*x + y*y + z*z;

	    asm("prefetcht0 512(%rdi)");
	    rinv = __builtin_ia32_rsqrtps(r2);
	
	    asm("prefetcht0 576(%rdi)");
	    /* Newton-Raphson */
	    t = rinv;
	    r2 *= rinv;
	    rinv *= r2;
	    rinv -= three;
	    rinv *= t;
	    rinv *= half;		/* flips sign to avoid storing -0.5 */
	    /* end Newton-Raphson */
	    
	    asm("prefetcht0 640(%rdi)");
	    t = rinv;
	    eqe = t*massp;
	    rinv2 = t*t;
	    t *= rinv2;
	    phi += eqe;
	    eqe *= rinv2;
	    a0 += x*eqe;
	    a1 += y*eqe;
	    a2 += z*eqe;
	    
	    t *= three*rinv2*R*R;
	    eq0 = t*(qxx*x + qxy*y + qxz*z);
	    eq1 = t*(qyy*y + qxy*x + qyz*z);
	    eq2 = t*(-(qxx + qyy)*z + qxz*x + qyz*y);
	    eqe = half*(eq0*x + eq1*y + eq2*z);
	    phi -= eqe;
	    eqe *= five*rinv2;
	    a0 += x*eqe - eq0;
	    a1 += y*eqe - eq1;
	    a2 += z*eqe - eq2;
	}
	a0 = __builtin_ia32_haddps(a0, a0);
	a1 = __builtin_ia32_haddps(a1, a1);
	a2 = __builtin_ia32_haddps(a2, a2);
	phi = __builtin_ia32_haddps(phi, phi);
	accp[i+0] += __builtin_ia32_haddps(a0, a0)[0];
	accp[i+1] += __builtin_ia32_haddps(a1, a1)[0];
	accp[i+2] += __builtin_ia32_haddps(a2, a2)[0];
	accp[i+3] += __builtin_ia32_haddps(phi, phi)[0];
    }
}
#endif
