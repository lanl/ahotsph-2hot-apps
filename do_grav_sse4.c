/*
 * Copyright 2012-2014 Michael S. Warren. All Rights Reserved.
 */
#ifndef __AVX__
#include <assert.h>
#include "order.h"
#include "segment.h"
#include "vec.h"

#define mass f[0]
#define xp f[1]
#define yp f[2]
#define zp f[3]
#define R  f[4]
#define qx f[5]
#define qy f[6]
#define qz f[7]
#define qxx f[8]
#define qxy f[9]
#define qyy f[10]
#define qxz f[11]
#define qyz f[12]
#define QSZ 13
#define qxxx f[13]
#define qxxy f[14]
#define qxyy f[15]
#define qyyy f[16]
#define qxxz f[17]
#define qxyz f[18]
#define qyyz f[19]
#define qxxxx f[20]
#define qxxxy f[21]
#define qxxyy f[22]
#define qxyyy f[23]
#define qyyyy f[24]
#define qxxxz f[25]
#define qxxyz f[26]
#define qxyyz f[27]
#define qyyyz f[28]
#define HSZ 29

/* 217 flops, 100 bytes */
void
do_gravdh_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

        t *= R*rinv2;
        eq0 = t*qx;
        eq1 = t*qy;
        eq2 = t*qz;
        eqe = x*eq0 + y*eq1 + z*eq2;
        phi += eqe;
        eqe *= three*rinv2;
        a0 += x*eqe - eq0;
        a1 += y*eqe - eq1;
        a2 += z*eqe - eq2;

	__asm__("prefetcht0 704(%rdi)");
        t *= three*rinv2*R;
        eq0 = t*(qxx*x + qxy*y + qxz*z);
        eq1 = t*(qyy*y + qxy*x + qyz*z);
        eq2 = t*(-(qxx + qyy)*z + qxz*x + qyz*y);
        eqe = half*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= five*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 768(%rdi)");
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

	__asm__("prefetcht0 832(%rdi)");
        eq0 = t*(qxxx*xx + qxyy*yy + qxxy*xy + qxxz*xz + qxyz*yz);
        eq1 = t*(qxyy*xy + qxxy*xx + qyyy*yy + qyyz*yz + qxyz*xz);
        eq2 = t*(-(qxxx + qxyy)*xz - (qxxy + qyyy)*yz + qxxz*xx + qyyz*yy + qxyz*xy);
        eqe = third*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= seven*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 896(%rdi)");
        t *= seven*rinv2*R;
        xxy = y*xx;
        xyy = x*yy;
        xyz = xy*z;
        eq0 = t*(qxxxx*xxx + qxyyy*yyy + qxxxy*xxy + qxxxz*xxz + qxxyy*xyy + qxxyz*xyz + qxyyz*yyz);
        eq1 = t*(qxyyy*xyy + qxxxy*xxx + qyyyy*yyy + qyyyz*yyz + qxxyy*xxy + qxxyz*xxz + qxyyz*xyz);
        eq2 = t*(-qxxxx*xxz - (qxyyy + qxxxy)*xyz - qyyyy*yyz + qxxxz*xxx + qyyyz*yyy  - qxxyy*(xxz + yyz) + qxxyz*xxy + qxyyz*xyy);
	f += HSZ;
        eqe = quarter*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= nine*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

void
do_gravdh_amd6100_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

	t *= R*rinv2;
        eq0 = t*qx;
        eq1 = t*qy;
        eq2 = t*qz;
        eqe = x*eq0;
	eqe += y*eq1;
	eqe += z*eq2;
        phi += eqe;
        eqe *= three;
	eqe *= rinv2;
        a0 += x*eqe - eq0;
        a1 += y*eqe - eq1;
        a2 += z*eqe - eq2;

	__asm__("prefetcht0 704(%rdi)");
        t *= three;
        t *= rinv2;
        t *= R;
        eq0 = qxx*x;
        eq1 = qyy*y;
        eq2 = -qxx;
        eq2 -= qyy;
        eq2 *= z;
        eq0 += qxy*y;
        eq1 += qxy*x;
        eq2 += qxz*x;
        eq0 += qxz*z;
        eq1 += qyz*z;
        eq2 += qyz*y;
	eq0 *= t;
        eq1 *= t;
        eq2 *= t;
        eqe = eq0*x;
        eqe += eq1*y;
        eqe += eq2*z;
        eqe *= half;
	phi += eqe;
	eqe *= five;
	eqe *= rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 768(%rdi)");
        t *= five;
        t *= rinv2;
        t *= R;
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

	__asm__("prefetcht0 832(%rdi)");
        eq0 = qxxx*xx;
        eq1 = qxyy*xy;
        eq2 = -(qxxx + qxyy);
	eq2 *= xz;
        eq0 += qxyy*yy;
        eq1 += qxxy*xx;
        eq2 -= (qxxy + qyyy)*yz;
        eq0 += qxxy*xy;
        eq1 += qyyy*yy;
        eq2 += qxxz*xx;
        eq0 += qxxz*xz;
        eq1 += qyyz*yz;
        eq2 += qyyz*yy;
        eq0 += qxyz*yz;
        eq1 += qxyz*xz;
        eq2 += qxyz*xy;
        eq0 *= t;
        eq1 *= t;
        eq2 *= t;
        eqe = eq0*x;
        eqe += eq1*y;
	eqe += eq2*z;
	eqe *= third;
	phi += eqe;
	eqe *= seven*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 896(%rdi)");
        t *= seven*rinv2*R;
        xxy = y*xx;
        xyy = x*yy;
        xyz = xy*z;
        eq2 = -qxxxx*xxz;
        eq0 = qxxxx*xxx;
        eq1 = qxyyy*xyy;
        eq2 -= (qxyyy + qxxxy)*xyz;
        eq0 += qxyyy*yyy;
        eq1 += qxxxy*xxx;
        eq2 -= qyyyy*yyz;
        eq0 += qxxxy*xxy;
        eq1 += qyyyy*yyy;
        eq2 += qxxxz*xxx;
        eq0 += qxxxz*xxz;
        eq1 += qyyyz*yyz;
        eq2 += qyyyz*yyy;
        eq0 += qxxyy*xyy;
        eq1 += qxxyy*xxy;
	eq2 -= qxxyy*(xxz + yyz);
        eq0 += qxxyz*xyz;
        eq1 += qxxyz*xxz;
	eq2 += qxxyz*xxy;
        eq0 += qxyyz*yyz;
        eq1 += qxyyz*xyz;
	eq2 += qxyyz*xyy;
	eq0 *= t;
	eq1 *= t;
        eq2 *= t;
	f += HSZ;
        eqe = eq0*x;
        eqe += eq1*y;
        eqe += eq2*z;
        eqe *= quarter;
	phi += eqe;
	eqe *= nine*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

/* 66 flops, 40 bytes */
void
do_gravdq_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	t *= rinv2;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

        t *= R;
        eq0 = t*qx;
        eq1 = t*qy;
        eq2 = t*qz;
        eqe = x*eq0 + y*eq1 + z*eq2;
        phi += eqe;
        eqe *= three*rinv2;
        a0 += x*eqe - eq0;
        a1 += y*eqe - eq1;
        a2 += z*eqe - eq2;

        t *= three*rinv2*R;
        eq0 = t*(qxx*x + qxy*y + qxz*z);
        eq1 = t*(qyy*y + qxy*x + qyz*z);
        eq2 = t*(-(qxx + qyy)*z + qxz*x + qyz*y);
	f += QSZ;
        eqe = half*(eq0*x + eq1*y + eq2*z);
	phi -= eqe;
	eqe *= five*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

#ifndef DIPOLE
#undef mass
#undef xp
#undef yp
#undef zp
#undef R
#undef qx
#undef qy
#undef qz
#undef qxx
#undef qxy
#undef qyy
#undef qxz
#undef qyz
#undef QSZ
#undef qxxx
#undef qxxy
#undef qxyy
#undef qyyy
#undef qxxz
#undef qxyz
#undef qyyz
#undef qxxxx
#undef qxxxy
#undef qxxyy
#undef qxyyy
#undef qyyyy
#undef qxxxz
#undef qxxyz
#undef qxyyz
#undef qyyyz
#undef HSZ

#define mass f[0]
#define xp f[1]
#define yp f[2]
#define zp f[3]
#define R  f[4]
#define qxx f[5]
#define qxy f[6]
#define qyy f[7]
#define qxz f[8]
#define qyz f[9]
#define QSZ 10
#define qxxx f[10]
#define qxxy f[11]
#define qxyy f[12]
#define qyyy f[13]
#define qxxz f[14]
#define qxyz f[15]
#define qyyz f[16]
#define qxxxx f[17]
#define qxxxy f[18]
#define qxxyy f[19]
#define qxyyy f[20]
#define qyyyy f[21]
#define qxxxz f[22]
#define qxxyz f[23]
#define qxyyz f[24]
#define qyyyz f[25]
#define HSZ 26
#endif

/* 217 flops, 100 bytes */
void
do_gravh_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	t *= rinv2;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

	__asm__("prefetcht0 704(%rdi)");
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

	__asm__("prefetcht0 768(%rdi)");
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

	__asm__("prefetcht0 832(%rdi)");
        eq0 = t*(qxxx*xx + qxyy*yy + qxxy*xy + qxxz*xz + qxyz*yz);
        eq1 = t*(qxyy*xy + qxxy*xx + qyyy*yy + qyyz*yz + qxyz*xz);
        eq2 = t*(-(qxxx + qxyy)*xz - (qxxy + qyyy)*yz + qxxz*xx + qyyz*yy + qxyz*xy);
        eqe = third*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= seven*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 896(%rdi)");
        t *= seven*rinv2*R;
        xxy = y*xx;
        xyy = x*yy;
        xyz = xy*z;
        eq0 = t*(qxxxx*xxx + qxyyy*yyy + qxxxy*xxy + qxxxz*xxz + qxxyy*xyy + qxxyz*xyz + qxyyz*yyz);
        eq1 = t*(qxyyy*xyy + qxxxy*xxx + qyyyy*yyy + qyyyz*yyz + qxxyy*xxy + qxxyz*xxz + qxyyz*xyz);
        eq2 = t*(-qxxxx*xxz - (qxyyy + qxxxy)*xyz - qyyyy*yyz + qxxxz*xxx + qyyyz*yyy  - qxxyy*(xxz + yyz) + qxxyz*xxy + qxyyz*xyy);
	f += HSZ;
        eqe = quarter*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= nine*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

void
do_gravh_amd6100_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	t *= rinv2;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

	__asm__("prefetcht0 704(%rdi)");
        t *= three;
        t *= rinv2;
        t *= R;
        t *= R;
        eq0 = qxx*x;
        eq1 = qyy*y;
        eq2 = -qxx;
        eq2 -= qyy;
        eq2 *= z;
        eq0 += qxy*y;
        eq1 += qxy*x;
        eq2 += qxz*x;
        eq0 += qxz*z;
        eq1 += qyz*z;
        eq2 += qyz*y;
	eq0 *= t;
        eq1 *= t;
        eq2 *= t;
        eqe = eq0*x;
        eqe += eq1*y;
        eqe += eq2*z;
        eqe *= half;
	phi += eqe;
	eqe *= five;
	eqe *= rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 768(%rdi)");
        t *= five;
        t *= rinv2;
        t *= R;
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

	__asm__("prefetcht0 832(%rdi)");
        eq0 = qxxx*xx;
        eq1 = qxyy*xy;
        eq2 = -(qxxx + qxyy);
	eq2 *= xz;
        eq0 += qxyy*yy;
        eq1 += qxxy*xx;
        eq2 -= (qxxy + qyyy)*yz;
        eq0 += qxxy*xy;
        eq1 += qyyy*yy;
        eq2 += qxxz*xx;
        eq0 += qxxz*xz;
        eq1 += qyyz*yz;
        eq2 += qyyz*yy;
        eq0 += qxyz*yz;
        eq1 += qxyz*xz;
        eq2 += qxyz*xy;
        eq0 *= t;
        eq1 *= t;
        eq2 *= t;
        eqe = eq0*x;
        eqe += eq1*y;
	eqe += eq2*z;
	eqe *= third;
	phi += eqe;
	eqe *= seven*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 896(%rdi)");
        t *= seven*rinv2*R;
        xxy = y*xx;
        xyy = x*yy;
        xyz = xy*z;
        eq2 = -qxxxx*xxz;
        eq0 = qxxxx*xxx;
        eq1 = qxyyy*xyy;
        eq2 -= (qxyyy + qxxxy)*xyz;
        eq0 += qxyyy*yyy;
        eq1 += qxxxy*xxx;
        eq2 -= qyyyy*yyz;
        eq0 += qxxxy*xxy;
        eq1 += qyyyy*yyy;
        eq2 += qxxxz*xxx;
        eq0 += qxxxz*xxz;
        eq1 += qyyyz*yyz;
        eq2 += qyyyz*yyy;
        eq0 += qxxyy*xyy;
        eq1 += qxxyy*xxy;
	eq2 -= qxxyy*(xxz + yyz);
        eq0 += qxxyz*xyz;
        eq1 += qxxyz*xxz;
	eq2 += qxxyz*xxy;
        eq0 += qxyyz*yyz;
        eq1 += qxyyz*xyz;
	eq2 += qxyyz*xyy;
	eq0 *= t;
	eq1 *= t;
        eq2 *= t;
	f += HSZ;
        eqe = eq0*x;
        eqe += eq1*y;
        eqe += eq2*z;
        eqe *= quarter;
	phi += eqe;
	eqe *= nine*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

/* 66 flops, 40 bytes */
void
do_gravq_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
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
	f += QSZ;
        eqe = half*(eq0*x + eq1*y + eq2*z);
	phi -= eqe;
	eqe *= five*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

/* 26 flops, 16 bytes */
void
do_grav_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf eps2 = {*e**e, *e**e, *e**e, *e**e};

    while (f < fend) {
	v4sf x = ppos0 - xp;
	v4sf y = ppos1 - yp;
	v4sf z = ppos2 - zp;

	v4sf r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	v4sf rinv = __builtin_ia32_rsqrtps(r2);
	v4sf mask = __builtin_ia32_cmpleps(eps2, r2);

	/* Newton-Raphson */
	v4sf t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	rinv = __builtin_ia32_andps(mask, rinv);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi += rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

typedef struct {
    v4sf phi; 
    v4sf a;
} sret_s;

/* Dehnen K1 Compensating Kernel */
static sret_s
sK1(const v4sf r2, const float e)
{
    sret_s r;
    const v4sf one = {1.0f, 1.0f, 1.0f, 1.0f};
    const v4sf f1 = {-3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f};
    const v4sf f2 = {135.0f/16.0f, 135.0f/16.0f, 135.0f/16.0f, 135.0f/16.0f};
    const v4sf p1 = {-1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f};
    const v4sf p2 = {3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f};
    const v4sf p3 = {-45.0f/32.0f, -45.0f/32.0f, -45.0f/32.0f, -45.0f/32.0f};
    const v4sf eps_inv = {e, e, e, e};

    v4sf u2 = r2 * eps_inv * eps_inv;
    v4sf mask = __builtin_ia32_cmpleps(u2, one);
    u2 -= one;

    v4sf t = p3*u2;
    t += p2;
    t *= u2;
    t += p1;
    t *= u2;
    t += one;
    t *= eps_inv;
    r.phi = -__builtin_ia32_andps(mask, t);
	
    t = f2*u2;
    t += f1;
    t *= u2;
    t += one;
    t *= eps_inv * eps_inv * eps_inv;
    r.a = -__builtin_ia32_andps(mask, t);

    return r;
}

void
do_grav_sK1_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf eps2 = {1.0f/(*e**e), 1.0f/(*e**e), 1.0f/(*e**e), 1.0f/(*e**e)};

    while (f < fend) {
	v4sf x = ppos0 - xp;
	v4sf y = ppos1 - yp;
	v4sf z = ppos2 - zp;

	v4sf r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	v4sf rinv = __builtin_ia32_rsqrtps(r2);
	v4sf mask = __builtin_ia32_cmpleps(eps2, r2);
	unsigned int nmask = vsf_count(mask);
	if (nmask != NSSE) {
	    *ncut += NSSE-nmask;
	    sret_s r = sK1(r2, *e);
	    phi += mass * r.phi;
	    r.a *= mass;
	    a0 += x * r.a;
	    a1 += y * r.a;
	    a2 += z * r.a;
	}

	/* Newton-Raphson */
	v4sf t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	rinv = __builtin_ia32_andps(mask, rinv);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi += rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

void
do_gravmm_sK1_sse4(const float *xyz, const int stride, const float pmass, const segment *mm, const int mm_n, 
		   const int source_n, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    vsf mmass = vsf_scalar(pmass);
    vsf xx, yy, zz;
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf eps2 = vsf_scalar(1.0f/(*e**e));
    int source_interactions = 0;

    int i = 0;
    int j = 0;
    const float *ff = xyz + mm[0].base * stride;
    while (i < mm_n) {
	/* Load vsf vectors from pairs of start, count indices into array of x,y,z positions */
	for (int k = 0; k < NSSE; k++) {
	    xx[k] = ff[j*stride+0];
	    yy[k] = ff[j*stride+1];
	    zz[k] = ff[j*stride+2];
	    source_interactions++;
	    if (++j == mm[i].length) {
		ff = xyz + mm[++i].base * stride;
		j = 0;
		if (i >= mm_n) {
		    while (++k < NSSE) {
			mmass[k] = 0.0f;
		    }
		    break;
		}
	    } 
	}

	vsf x = ppos0 - xx;
	vsf y = ppos1 - yy;
	vsf z = ppos2 - zz;

	vsf r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	vsf rinv = vsf_rsqrt(r2);
	vsf mask = vsf_cmple(eps2, r2);
	unsigned int nmask = vsf_count(mask);
	if (nmask != NSSE) {
	    *ncut += NSSE-nmask;
	    sret_s r = sK1(r2, *e);
	    phi += mmass * r.phi;
	    r.a *= mmass;
	    a0 += x * r.a;
	    a1 += y * r.a;
	    a2 += z * r.a;
	}

	/* Newton-Raphson */
	vsf t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	rinv = vsf_and(mask, rinv);

	t = rinv*rinv;
	rinv *= mmass;
	phi += rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
    assert(source_interactions == source_n);
}


/* Spline Kernel */

void
do_gravsS_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf u2, t, rinv, r;
    v4sf mask, tm, tp;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps_inv = {*e, *e, *e, *e};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf minus_half = {-0.5f, -0.5f, -0.5f, -0.5f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf c48o5 = {48.0/5.0, 48.0/5.0, 48.0/5.0, 48.0/5.0};
    const v4sf c32o5 = {32.0/5.0, 32.0/5.0, 32.0/5.0, 32.0/5.0};
    const v4sf c16o3 = {16.0/3.0, 16.0/3.0, 16.0/3.0, 16.0/3.0};
    const v4sf c14o5 = {14.0/5.0, 14.0/5.0, 14.0/5.0, 14.0/5.0};
    const v4sf c32o15 = {32.0/15.0, 32.0/15.0, 32.0/15.0, 32.0/15.0};
    const v4sf c16 = {16.0, 16.0, 16.0, 16.0};
    const v4sf c32o3 = {32.0/3.0, 32.0/3.0, 32.0/3.0, 32.0/3.0};
    const v4sf c1o15 = {1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0};
    const v4sf c16o5 = {16.0/5.0, 16.0/5.0, 16.0/5.0, 16.0/5.0};
    const v4sf c192o5 = {192.0/5.0, 192.0/5.0, 192.0/5.0, 192.0/5.0};
    const v4sf c32 = {32.0, 32.0, 32.0, 32.0};
    const v4sf c48 = {48.0, 48.0, 48.0, 48.0};
    const v4sf c64o3 = {64.0/3.0, 64.0/3.0, 64.0/3.0, 64.0/3.0};
	

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;

	__asm__("prefetcht0 512(%rdi)");

	rinv = __builtin_ia32_rsqrtps(u2);
	/* Newton-Raphson */
	t = rinv;
	rinv *= rinv;
	rinv *= u2;
	rinv -= three;
	rinv *= t;
	rinv *= minus_half;
	/* end Newton-Raphson */

	r = u2*rinv;
	mask = __builtin_ia32_cmpltps(r, half);

	tm = u2 * (u2 * (c48o5 - c32o5*r) - c16o3) + c14o5;
	tp = u2 * ((u2 * (c32o15*r - c48o5) + c16*r) - c32o3) - c1o15*rinv + c16o5;
#if 0				/* SSE4.1 */
	t = __builtin_ia32_blendvps(tm, tp, mask);
#else
	tm = __builtin_ia32_andps(mask, tm);
	tp = __builtin_ia32_andnps(mask, tp);
	t = __builtin_ia32_orps(tm, tp);
#endif

	phi -= mass * eps_inv * t;

	tm = u2 * (c32*r - c192o5) + c32o3;
	tp = u2 * (c192o5 - c32o3*r) - c48*r + c64o3 - c1o15*rinv*rinv*rinv;
#if 0
	t = __builtin_ia32_blendvps(tm, tp, mask);
#else
	tm = __builtin_ia32_andps(mask, tm);
	tp = __builtin_ia32_andnps(mask, tp);
	t = __builtin_ia32_orps(tm, tp);
#endif
	
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

void
do_gravsF1_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf u2, t;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps_inv = {*e, *e, *e, *e};
    const v4sf one = {1.0f, 1.0f, 1.0f, 1.0f};
    const v4sf f1 = {-3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f};
    const v4sf p1 = {-1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f};
    const v4sf p2 = {3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;
	u2 -= one;

	__asm__("prefetcht0 512(%rdi)");

	t = p2*u2;
	t += p1;
	t *= u2;
	t += one;

	phi -= mass * eps_inv * t;
	
	t = f1*u2;
	t += one;
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}


/* F2 "Epanechnikov" Kernel */

void
do_gravsF2_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf u2, t;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps_inv = {*e, *e, *e, *e};
    const v4sf one = {1.0f, 1.0f, 1.0f, 1.0f};
    const v4sf f1 = {-3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f};
    const v4sf f2 = {15.0f/8.0f, 15.0f/8.0f, 15.0f/8.0f, 15.0f/8.0f};
    const v4sf p1 = {-1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f};
    const v4sf p2 = {3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f};
    const v4sf p3 = {-5.0f/16.0f, -5.0f/16.0f, -5.0f/16.0f, -5.0f/16.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;
	u2 -= one;

	__asm__("prefetcht0 512(%rdi)");

	t = p3*u2;
	t += p2;
	t *= u2;
	t += p1;
	t *= u2;
	t += one;

	phi -= mass * eps_inv * t;
	
	t = f2*u2;
	t += f1;
	t *= u2;
	t += one;
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

/* Dehnen K1 Compensating Kernel */
void
do_gravsK1_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf u2, t;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps_inv = {*e, *e, *e, *e};
    const v4sf one = {1.0f, 1.0f, 1.0f, 1.0f};
    const v4sf f1 = {-3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f, -3.0f/2.0f};
    const v4sf f2 = {135.0f/16.0f, 135.0f/16.0f, 135.0f/16.0f, 135.0f/16.0f};
    const v4sf p1 = {-1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f, -1.0f/2.0f};
    const v4sf p2 = {3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f, 3.0f/8.0f};
    const v4sf p3 = {-45.0f/32.0f, -45.0f/32.0f, -45.0f/32.0f, -45.0f/32.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;
	u2 -= one;

	__asm__("prefetcht0 512(%rdi)");

	t = p3*u2;
	t += p2;
	t *= u2;
	t += p1;
	t *= u2;
	t += one;

	phi -= mass * eps_inv * t;
	
	t = f2*u2;
	t += f1;
	t *= u2;
	t += one;
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

/* Uniform density F0 kernel */

void
do_gravsU_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf u2, t;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps_inv = {*e, *e, *e, *e};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	u2 = x*x + y*y + z*z;
	t = eps_inv * eps_inv;
	u2 *= t;

	__asm__("prefetcht0 512(%rdi)");
	
	phi -= mass * eps_inv * half * (three - u2);

	t *= mass * eps_inv;
	f += 4;
	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

/* Compact Plummer kernel */

#define fac98pow32 1.1932426932523f
#define fac98pow12 1.0606601717798f

void
do_gravsCP_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    v4sf t, r2, rinv;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps2 = {*eps2p*0.125f, *eps2p*0.125f, *eps2p*0.125f, *eps2p*0.125f};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi += rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0]*fac98pow12;
    acc[0] += __builtin_ia32_haddps(a0, a0)[0]*fac98pow32;
    acc[1] += __builtin_ia32_haddps(a1, a1)[0]*fac98pow32;
    acc[2] += __builtin_ia32_haddps(a2, a2)[0]*fac98pow32;
}


void
do_grav11bits_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    v4sf t, r2, rinv;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi -= rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] -= __builtin_ia32_haddps(a0, a0)[0];
    acc[1] -= __builtin_ia32_haddps(a1, a1)[0];
    acc[2] -= __builtin_ia32_haddps(a2, a2)[0];
}


void
do_gravph_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps2 = {*eps2p, *eps2p, *eps2p, *eps2p};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	t *= rinv2;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

	__asm__("prefetcht0 704(%rdi)");
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

	__asm__("prefetcht0 768(%rdi)");
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

	__asm__("prefetcht0 832(%rdi)");
        eq0 = t*(qxxx*xx + qxyy*yy + qxxy*xy + qxxz*xz + qxyz*yz);
        eq1 = t*(qxyy*xy + qxxy*xx + qyyy*yy + qyyz*yz + qxyz*xz);
        eq2 = t*(-(qxxx + qxyy)*xz - (qxxy + qyyy)*yz + qxxz*xx + qyyz*yy + qxyz*xy);
        eqe = third*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= seven*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 896(%rdi)");
        t *= seven*rinv2*R;
        xxy = y*xx;
        xyy = x*yy;
        xyz = xy*z;
        eq0 = t*(qxxxx*xxx + qxyyy*yyy + qxxxy*xxy + qxxxz*xxz + qxxyy*xyy + qxxyz*xyz + qxyyz*yyz);
        eq1 = t*(qxyyy*xyy + qxxxy*xxx + qyyyy*yyy + qyyyz*yyz + qxxyy*xxy + qxxyz*xxz + qxyyz*xyz);
        eq2 = t*(-qxxxx*xxz - (qxyyy + qxxxy)*xyz - qyyyy*yyz + qxxxz*xxx + qyyyz*yyy  - qxxyy*(xxz + yyz) + qxxyz*xxy + qxyyz*xyy);
	f += HSZ;
        eqe = quarter*(eq0*x + eq1*y + eq2*z);
	phi += eqe;
	eqe *= nine*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

void
do_gravph_amd6100_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf xx, xy, yy, xz, yz, zz;
    v4sf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps2 = {*eps2p, *eps2p, *eps2p, *eps2p};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf third = {(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0),(float)(1.0/3.0)};
    const v4sf quarter = {0.25f, 0.25f, 0.25f, 0.25f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};
    const v4sf seven = {7.0f, 7.0f, 7.0f, 7.0f};
    const v4sf nine = {9.0f, 9.0f, 9.0f, 9.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
	rinv2 = t*t;
	t *= rinv2;
	phi += eqe;
	eqe *= rinv2;
	a0 += x*eqe;
	a1 += y*eqe;
	a2 += z*eqe;

	__asm__("prefetcht0 704(%rdi)");
        t *= three;
        t *= rinv2;
        t *= R;
        t *= R;
        eq0 = qxx*x;
        eq1 = qyy*y;
        eq2 = -qxx;
        eq2 -= qyy;
        eq2 *= z;
        eq0 += qxy*y;
        eq1 += qxy*x;
        eq2 += qxz*x;
        eq0 += qxz*z;
        eq1 += qyz*z;
        eq2 += qyz*y;
	eq0 *= t;
        eq1 *= t;
        eq2 *= t;
        eqe = eq0*x;
        eqe += eq1*y;
        eqe += eq2*z;
        eqe *= half;
	phi += eqe;
	eqe *= five;
	eqe *= rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 768(%rdi)");
        t *= five;
        t *= rinv2;
        t *= R;
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

	__asm__("prefetcht0 832(%rdi)");
        eq0 = qxxx*xx;
        eq1 = qxyy*xy;
        eq2 = -(qxxx + qxyy);
	eq2 *= xz;
        eq0 += qxyy*yy;
        eq1 += qxxy*xx;
        eq2 -= (qxxy + qyyy)*yz;
        eq0 += qxxy*xy;
        eq1 += qyyy*yy;
        eq2 += qxxz*xx;
        eq0 += qxxz*xz;
        eq1 += qyyz*yz;
        eq2 += qyyz*yy;
        eq0 += qxyz*yz;
        eq1 += qxyz*xz;
        eq2 += qxyz*xy;
        eq0 *= t;
        eq1 *= t;
        eq2 *= t;
        eqe = eq0*x;
        eqe += eq1*y;
	eqe += eq2*z;
	eqe *= third;
	phi += eqe;
	eqe *= seven*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;

	__asm__("prefetcht0 896(%rdi)");
        t *= seven*rinv2*R;
        xxy = y*xx;
        xyy = x*yy;
        xyz = xy*z;
        eq2 = -qxxxx*xxz;
        eq0 = qxxxx*xxx;
        eq1 = qxyyy*xyy;
        eq2 -= (qxyyy + qxxxy)*xyz;
        eq0 += qxyyy*yyy;
        eq1 += qxxxy*xxx;
        eq2 -= qyyyy*yyz;
        eq0 += qxxxy*xxy;
        eq1 += qyyyy*yyy;
        eq2 += qxxxz*xxx;
        eq0 += qxxxz*xxz;
        eq1 += qyyyz*yyz;
        eq2 += qyyyz*yyy;
        eq0 += qxxyy*xyy;
        eq1 += qxxyy*xxy;
	eq2 -= qxxyy*(xxz + yyz);
        eq0 += qxxyz*xyz;
        eq1 += qxxyz*xxz;
	eq2 += qxxyz*xxy;
        eq0 += qxyyz*yyz;
        eq1 += qxyyz*xyz;
	eq2 += qxyyz*xyy;
	eq0 *= t;
	eq1 *= t;
        eq2 *= t;
	f += HSZ;
        eqe = eq0*x;
        eqe += eq1*y;
        eqe += eq2*z;
        eqe *= quarter;
	phi += eqe;
	eqe *= nine*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}

void
do_gravpq_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    v4sf t, r2, rinv, rinv2;
    v4sf x, y, z;
    v4sf eqe, eq0, eq1, eq2;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps2 = {*eps2p, *eps2p, *eps2p, *eps2p};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    const v4sf five = {5.0f, 5.0f, 5.0f, 5.0f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	__asm__("prefetcht0 576(%rdi)");
	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	
	__asm__("prefetcht0 640(%rdi)");
	t = rinv;
	eqe = t*mass;
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
	f += QSZ;
        eqe = half*(eq0*x + eq1*y + eq2*z);
	phi -= eqe;
	eqe *= five*rinv2;
	a0 += x*eqe - eq0;
	a1 += y*eqe - eq1;
	a2 += z*eqe - eq2;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}


void
do_gravp_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    v4sf t, r2, rinv;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps2 = {*eps2p, *eps2p, *eps2p, *eps2p};
    const v4sf three = {3.0f, 3.0f, 3.0f, 3.0f};
    const v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	/* Newton-Raphson */
	t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi += rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] += __builtin_ia32_haddps(a0, a0)[0];
    acc[1] += __builtin_ia32_haddps(a1, a1)[0];
    acc[2] += __builtin_ia32_haddps(a2, a2)[0];
}


void
do_gravp11bits_sse4(const v4sf *f, const v4sf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    v4sf t, r2, rinv;
    v4sf x, y, z;
    v4sf a0 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a1 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf a2 = {0.0f, 0.0f, 0.0f, 0.0f};
    v4sf phi = {0.0f, 0.0f, 0.0f, 0.0f};
    const v4sf ppos0 = {pos0[0], pos0[0], pos0[0], pos0[0]};
    const v4sf ppos1 = {pos0[1], pos0[1], pos0[1], pos0[1]};
    const v4sf ppos2 = {pos0[2], pos0[2], pos0[2], pos0[2]};
    const v4sf eps2 = {*eps2p, *eps2p, *eps2p, *eps2p};

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = __builtin_ia32_rsqrtps(r2);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi -= rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    phi = __builtin_ia32_haddps(phi, phi);
    a0 = __builtin_ia32_haddps(a0, a0);
    a1 = __builtin_ia32_haddps(a1, a1);
    a2 = __builtin_ia32_haddps(a2, a2);
    *phi0 += __builtin_ia32_haddps(phi, phi)[0];
    acc[0] -= __builtin_ia32_haddps(a0, a0)[0];
    acc[1] -= __builtin_ia32_haddps(a1, a1)[0];
    acc[2] -= __builtin_ia32_haddps(a2, a2)[0];
}
#endif
