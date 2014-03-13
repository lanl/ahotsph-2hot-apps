/*
 * Copyright 2012-2014 Michael S. Warren. All Rights Reserved.
 */
#ifdef __AVX__
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
do_gravdh_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0/3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravdh_amd6100_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0/3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

/* 66 flops, 40 bytes */
void
do_gravdq_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf five = vsf_scalar(5.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravh_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0/3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravh_amd6100_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0/3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
        eq0 = qxxxx*xxx;
        eq1 = qxyyy*xyy;
        eq2 = -qxxxx*xxz;
        eq0 += qxyyy*yyy;
        eq1 += qxxxy*xxx;
        eq2 -= (qxyyy + qxxxy)*xyz;
        eq0 += qxxxy*xxy;
        eq1 += qyyyy*yyy;
        eq2 -= qyyyy*yyz;
        eq0 += qxxxz*xxz;
        eq1 += qyyyz*yyz;
        eq2 += qxxxz*xxx;
        eq0 += qxxyy*xyy;
        eq1 += qxxyy*xxy;
        eq2 += qyyyz*yyy;
        eq0 += qxxyz*xyz;
        eq1 += qxxyz*xxz;
	eq2 -= qxxyy*(xxz + yyz);
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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravq_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf five = vsf_scalar(5.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


void
do_grav_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf eps2 = vsf_scalar(*e**e);

    while (f < fend) {
	vsf x = ppos0 - xp;
	vsf y = ppos1 - yp;
	vsf z = ppos2 - zp;

	vsf r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	vsf rinv = vsf_rsqrt(r2);
	vsf mask = vsf_cmple(eps2, r2);

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
	rinv *= mass;
	f += 4;
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
}

typedef struct {
    vsf phi; 
    vsf a;
} sret_s;

/* Dehnen K1 Compensating Kernel */
static sret_s
sK1(const vsf r2, const float e)
{
    sret_s r;
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f/2.0f);
    const vsf f2 = vsf_scalar(135.0f/16.0f);
    const vsf p1 = vsf_scalar(-1.0f/2.0f);
    const vsf p2 = vsf_scalar(3.0f/8.0f);
    const vsf p3 = vsf_scalar(-45.0f/32.0f);
    const vsf eps_inv = vsf_scalar(e);

    vsf u2 = r2 * eps_inv * eps_inv;
    vsf mask = vsf_cmple(u2, one);
    u2 -= one;

    vsf t = p3*u2;
    t += p2;
    t *= u2;
    t += p1;
    t *= u2;
    t += one;
    t *= eps_inv;
    r.phi = -vsf_and(mask, t);
	
    t = f2*u2;
    t += f1;
    t *= u2;
    t += one;
    t *= eps_inv * eps_inv * eps_inv;
    r.a = -vsf_and(mask, t);

    return r;
}

void
do_grav_sK1_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf eps2 = vsf_scalar(1.0f/(*e**e));

    while (f < fend) {
	vsf x = ppos0 - xp;
	vsf y = ppos1 - yp;
	vsf z = ppos2 - zp;

	vsf r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	vsf rinv = vsf_rsqrt(r2);
	vsf mask = vsf_cmple(eps2, r2);
	unsigned int nmask = __builtin_popcount(__builtin_ia32_movmskps256(mask));
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
	vsf t = rinv;
	r2 *= rinv;
	rinv *= r2;
	rinv -= three;
	rinv *= t;
	rinv *= half;		/* flips sign to avoid storing -0.5 */
	/* end Newton-Raphson */
	rinv = vsf_and(mask, rinv);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
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
}

void
do_gravmm_sK1_avx8(const float *xyz, const int stride, const float pmass, const segment *mm, const int mm_n, 
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
do_gravsS_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t, rinv, r;
    vsf mask, tm, tp;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf three = vsf_scalar(3.0f);
    const vsf minus_half = vsf_scalar(-0.5f);
    const vsf half = vsf_scalar(0.5f);
    const vsf c48o5 = vsf_scalar(48.0/5.0);
    const vsf c32o5 = vsf_scalar(32.0/5.0);
    const vsf c16o3 = vsf_scalar(16.0/3.0);
    const vsf c14o5 = vsf_scalar(14.0/5.0);
    const vsf c32o15 = vsf_scalar(32.0/15.0);
    const vsf c16 = vsf_scalar(16.0);
    const vsf c32o3 = vsf_scalar(32.0/3.0);
    const vsf c1o15 = vsf_scalar(1.0/15.0);
    const vsf c16o5 = vsf_scalar(16.0/5.0);
    const vsf c192o5 = vsf_scalar(192.0/5.0);
    const vsf c32 = vsf_scalar(32.0);
    const vsf c48 = vsf_scalar(48.0);
    const vsf c64o3 = vsf_scalar(64.0/3.0);
	
    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;

	__asm__("prefetcht0 512(%rdi)");

	rinv = vsf_rsqrt(u2);
	/* Newton-Raphson */
	t = rinv;
	rinv *= rinv;
	rinv *= u2;
	rinv -= three;
	rinv *= t;
	rinv *= minus_half;
	/* end Newton-Raphson */

	r = u2*rinv;
	mask = __builtin_ia32_cmpps256(r, half, 0x11); /* _CMP_LT_OQ */

	tm = u2 * (u2 * (c48o5 - c32o5*r) - c16o3) + c14o5;
	tp = u2 * ((u2 * (c32o15*r - c48o5) + c16*r) - c32o3) - c1o15*rinv + c16o5;
	t = __builtin_ia32_blendvps256(tm, tp, mask);

	phi -= mass * eps_inv * t;

	tm = u2 * (c32*r - c192o5) + c32o3;
	tp = u2 * (c192o5 - c32o3*r) - c48*r + c64o3 - c1o15*rinv*rinv*rinv;
	t = __builtin_ia32_blendvps256(tm, tp, mask);
	
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravsF1_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f/2.0f);
    const vsf p1 = vsf_scalar(-1.0f/2.0f);
    const vsf p2 = vsf_scalar(3.0f/8.0f);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


/* F2 "Epanechnikov" Kernel */

void
do_gravsF2_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f/2.0f);
    const vsf f2 = vsf_scalar(15.0f/8.0f);
    const vsf p1 = vsf_scalar(-1.0f/2.0f);
    const vsf p2 = vsf_scalar(3.0f/8.0f);
    const vsf p3 = vsf_scalar(-5.0f/16.0f);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

/* Dehnen K1 Compensating Kernel */
void
do_gravsK1_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f/2.0f);
    const vsf f2 = vsf_scalar(135.0f/16.0f);
    const vsf p1 = vsf_scalar(-1.0f/2.0f);
    const vsf p2 = vsf_scalar(3.0f/8.0f);
    const vsf p3 = vsf_scalar(-45.0f/32.0f);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

/* Dehnen K2 Compensating Kernel */
void
do_gravsK2_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t, mask;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf zero = vsf_scalar(0.0f);
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f/2.0f);
    const vsf f2 = vsf_scalar(15.0f/8.0f);
    const vsf f3 = vsf_scalar(-385.0f/32.0f);
    const vsf p1 = vsf_scalar(-1.0f/2.0f);
    const vsf p2 = vsf_scalar(3.0f/8.0f);
    const vsf p3 = vsf_scalar(-5.0f/16.0f);
    const vsf p4 = vsf_scalar(385.0f/256.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	/* t will overflow for m,x,y,z=0 padded entries */
	mask = __builtin_ia32_cmpps256(mass, zero, 0x04); /* _CMP_NEQ_UQ */
	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;
	u2 -= one;
	u2 = (vsf)((vsi)mask & (vsi)u2);

	__asm__("prefetcht0 512(%rdi)");

	t = p4*u2;
	t += p3;
	t *= u2;
	t += p2;
	t *= u2;
	t += p1;
	t *= u2;
	t += one;

	phi -= mass * eps_inv * t;
	
	t = f3*u2;
	t += f2;
	t *= u2;
	t += f1;
	t *= u2;
	t += one;
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

/* Dehnen K3 Compensating Kernel */
void
do_gravsK3_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t, mask;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf zero = vsf_scalar(0.0f);
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f/2.0f);
    const vsf f2 = vsf_scalar(15.0f/8.0f);
    const vsf f3 = vsf_scalar(-35.0f/16.0f);
    const vsf f4 = vsf_scalar(4095.0f/256.0f);
    const vsf p1 = vsf_scalar(-1.0f/2.0f);
    const vsf p2 = vsf_scalar(3.0f/8.0f);
    const vsf p3 = vsf_scalar(-5.0f/16.0f);
    const vsf p4 = vsf_scalar(35.0f/128.0f);
    const vsf p5 = vsf_scalar(-819.0f/512.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	/* t will overflow for m,x,y,z=0 padded entries */
	mask = __builtin_ia32_cmpps256(mass, zero, 0x04); /* _CMP_NEQ_UQ */
	u2 = x*x + y*y + z*z;
	u2 *= eps_inv * eps_inv;
	u2 -= one;
	u2 = (vsf)((vsi)mask & (vsi)u2);

	__asm__("prefetcht0 512(%rdi)");

	t = p5*u2;
	t += p4;
	t *= u2;
	t += p3;
	t *= u2;
	t += p2;
	t *= u2;
	t += p1;
	t *= u2;
	t += one;

	phi -= mass * eps_inv * t;
	
	t = f4*u2;
	t += f3;
	t *= u2;
	t += f2;
	t *= u2;
	t += f1;
	t *= u2;
	t += one;
	t *= mass * eps_inv * eps_inv * eps_inv;
	f += 4;

	a0 -= x*t;
	a1 -= y*t;
	a2 -= z*t;
    }
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


/* Uniform density F0 kernel */

void
do_gravsU_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf u2, t;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps_inv = vsf_scalar(*e);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}



/* Compact Plummer kernel */

#define fac98pow32 1.1932426932523f
#define fac98pow12 1.0606601717798f

void
do_gravsCP_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    vsf t, r2, rinv;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps2 = vsf_scalar(*eps2p*0.125f);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


void
do_grav11bits_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)
{
    vsf t, r2, rinv;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi -= rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


void
do_gravph_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps2 = vsf_scalar(*eps2p);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0/3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravph_amd6100_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps2 = vsf_scalar(*eps2p);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0/3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	/* Prefetching improves uncached performance by 10-20% */
	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
        eq0 = qxxxx*xxx;
        eq1 = qxyyy*xyy;
        eq2 = -qxxxx*xxz;
        eq0 += qxyyy*yyy;
        eq1 += qxxxy*xxx;
        eq2 -= (qxyyy + qxxxy)*xyz;
        eq0 += qxxxy*xxy;
        eq1 += qyyyy*yyy;
        eq2 -= qyyyy*yyz;
        eq0 += qxxxz*xxz;
        eq1 += qyyyz*yyz;
        eq2 += qxxxz*xxx;
        eq0 += qxxyy*xyy;
        eq1 += qxxyy*xxy;
        eq2 += qyyyz*yyy;
        eq0 += qxxyz*xyz;
        eq1 += qxxyz*xxz;
	eq2 -= qxxyy*(xxz + yyz);
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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

void
do_gravpq_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf eqe, eq0, eq1, eq2;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps2 = vsf_scalar(*eps2p);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf five = vsf_scalar(5.0f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


void
do_gravp_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    vsf t, r2, rinv;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps2 = vsf_scalar(*eps2p);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

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
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}


void
do_gravp11bits_avx8(const vsf *f, const vsf *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *eps2p, int *ncut)
{
    vsf t, r2, rinv;
    vsf x, y, z;
    vsf a0 = vsf_scalar(0.0f);
    vsf a1 = vsf_scalar(0.0f);
    vsf a2 = vsf_scalar(0.0f);
    vsf phi = vsf_scalar(0.0f);
    const vsf ppos0 = vsf_scalar(pos0[0]);
    const vsf ppos1 = vsf_scalar(pos0[1]);
    const vsf ppos2 = vsf_scalar(pos0[2]);
    const vsf eps2 = vsf_scalar(*eps2p);

    while (f < fend) {
	x = ppos0 - xp;
	y = ppos1 - yp;
	z = ppos2 - zp;

	r2 = x*x + y*y + z*z + eps2;

	__asm__("prefetcht0 512(%rdi)");
	rinv = vsf_rsqrt(r2);

	t = rinv*rinv;
	rinv *= mass;
	f += 4;
	phi -= rinv;
	t *= rinv;
	a0 += x*t;
	a1 += y*t;
	a2 += z*t;
    }
    *phi0 += vsf_hsum(phi);
    acc[0] += vsf_hsum(a0);
    acc[1] += vsf_hsum(a1);
    acc[2] += vsf_hsum(a2);
}

#endif
