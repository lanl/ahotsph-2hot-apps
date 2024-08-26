/*
 * Copyright 2012-2013 Michael S. Warren. All Rights Reserved.
 */
#include "vec.h"

/* Assumes -DDIPPOLE */

#define mass f[j + 0]
#define xp f[j + 1]
#define yp f[j + 2]
#define zp f[j + 3]
#define MSZ 4
#define R f[j + 4]
#define qx f[j + 5]
#define qy f[j + 6]
#define qz f[j + 7]
#define qxx f[j + 8]
#define qxy f[j + 9]
#define qyy f[j + 10]
#define qxz f[j + 11]
#define qyz f[j + 12]
#define QSZ 13
#define qxxx f[j + 13]
#define qxxy f[j + 14]
#define qxyy f[j + 15]
#define qyyy f[j + 16]
#define qxxz f[j + 17]
#define qxyz f[j + 18]
#define qyyz f[j + 19]
#define qxxxx f[j + 20]
#define qxxxy f[j + 21]
#define qxxyy f[j + 22]
#define qxyyy f[j + 23]
#define qyyyy f[j + 24]
#define qxxxz f[j + 25]
#define qxxyz f[j + 26]
#define qxyyz f[j + 27]
#define qyyyz f[j + 28]
#define HSZ 29

extern void WalkPoll(void);

void pHinteract(
    const float *p, float *accp, const int n, const int stride, const vsf *f, const int source_n) {
    int i, j;
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf xx, xy, yy, xz, yz, zz;
    vsf xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    vsf eqe, eq0, eq1, eq2;
    vsf a0, a1, a2, phi;
    const vsf zero = vsf_scalar(0.0f);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf third = vsf_scalar((float)(1.0 / 3.0));
    const vsf quarter = vsf_scalar(0.25f);
    const vsf five = vsf_scalar(5.0f);
    const vsf seven = vsf_scalar(7.0f);
    const vsf nine = vsf_scalar(9.0f);

    for (i = 0; i < n * stride; i += stride) {
        a0 = a1 = a2 = phi = zero;
        vsf ppos0 = vsf_scalar(p[i + 1]);
        vsf ppos1 = vsf_scalar(p[i + 2]);
        vsf ppos2 = vsf_scalar(p[i + 3]);

        for (j = 0; j < source_n * HSZ; j += HSZ) {
            x = ppos0 - xp;
            y = ppos1 - yp;
            z = ppos2 - zp;

            r2 = x * x + y * y + z * z;

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
            rinv *= half; /* flips sign to avoid storing -0.5 */
            /* end Newton-Raphson */

            __asm__("prefetcht0 640(%rdi)");
            t = rinv;
            eqe = t * mass;
            rinv2 = t * t;
            phi += eqe;
            eqe *= rinv2;
            a0 += x * eqe;
            a1 += y * eqe;
            a2 += z * eqe;

            t *= R * rinv2;
            eq0 = t * qx;
            eq1 = t * qy;
            eq2 = t * qz;
            eqe = x * eq0 + y * eq1 + z * eq2;
            phi += eqe;
            eqe *= three * rinv2;
            a0 += x * eqe - eq0;
            a1 += y * eqe - eq1;
            a2 += z * eqe - eq2;

            __asm__("prefetcht0 704(%rdi)");
            t *= three * rinv2 * R;
            eq0 = t * (qxx * x + qxy * y + qxz * z);
            eq1 = t * (qyy * y + qxy * x + qyz * z);
            eq2 = t * (-(qxx + qyy) * z + qxz * x + qyz * y);
            eqe = half * (eq0 * x + eq1 * y + eq2 * z);
            phi += eqe;
            eqe *= five * rinv2;
            a0 += x * eqe - eq0;
            a1 += y * eqe - eq1;
            a2 += z * eqe - eq2;

            __asm__("prefetcht0 768(%rdi)");
            t *= five * rinv2 * R;
            xx = half * x * x;
            xy = x * y;
            xz = x * z;
            yy = half * y * y;
            yz = y * z;
            zz = half * z * z;
            xxx = x * (third * xx - zz);
            xxz = z * (xx - third * zz);
            yyy = y * (third * yy - zz);
            yyz = z * (yy - third * zz);
            xx -= zz;
            yy -= zz;

            __asm__("prefetcht0 832(%rdi)");
            eq0 = t * (qxxx * xx + qxyy * yy + qxxy * xy + qxxz * xz + qxyz * yz);
            eq1 = t * (qxyy * xy + qxxy * xx + qyyy * yy + qyyz * yz + qxyz * xz);
            eq2 = t
                  * (-(qxxx + qxyy) * xz - (qxxy + qyyy) * yz + qxxz * xx + qyyz * yy + qxyz * xy);
            eqe = third * (eq0 * x + eq1 * y + eq2 * z);
            phi += eqe;
            eqe *= seven * rinv2;
            a0 += x * eqe - eq0;
            a1 += y * eqe - eq1;
            a2 += z * eqe - eq2;

            __asm__("prefetcht0 896(%rdi)");
            t *= seven * rinv2 * R;
            xxy = y * xx;
            xyy = x * yy;
            xyz = xy * z;
            eq0 = t
                  * (qxxxx * xxx + qxyyy * yyy + qxxxy * xxy + qxxxz * xxz + qxxyy * xyy
                     + qxxyz * xyz + qxyyz * yyz);
            eq1 = t
                  * (qxyyy * xyy + qxxxy * xxx + qyyyy * yyy + qyyyz * yyz + qxxyy * xxy
                     + qxxyz * xxz + qxyyz * xyz);
            eq2 = t
                  * (-qxxxx * xxz - (qxyyy + qxxxy) * xyz - qyyyy * yyz + qxxxz * xxx + qyyyz * yyy
                     - qxxyy * (xxz + yyz) + qxxyz * xxy + qxyyz * xyy);
            eqe = quarter * (eq0 * x + eq1 * y + eq2 * z);
            phi += eqe;
            eqe *= nine * rinv2;
            a0 += x * eqe - eq0;
            a1 += y * eqe - eq1;
            a2 += z * eqe - eq2;
        }
        accp[i + 0] += vsf_hsum(a0);
        accp[i + 1] += vsf_hsum(a1);
        accp[i + 2] += vsf_hsum(a2);
        accp[i + 3] += vsf_hsum(phi);
    }
}

void pQinteract(
    const float *p, float *accp, const int n, const int stride, const vsf *f, const int source_n) {
    int i, j;
    vsf t, r2, rinv, rinv2;
    vsf x, y, z;
    vsf eqe, eq0, eq1, eq2;
    vsf a0, a1, a2, phi;
    const vsf zero = vsf_scalar(0.0f);
    const vsf three = vsf_scalar(3.0f);
    const vsf half = vsf_scalar(0.5f);
    const vsf five = vsf_scalar(5.0f);

    for (i = 0; i < n * stride; i += stride) {
        a0 = a1 = a2 = phi = zero;
        vsf ppos0 = vsf_scalar(p[i + 1]);
        vsf ppos1 = vsf_scalar(p[i + 2]);
        vsf ppos2 = vsf_scalar(p[i + 3]);

        for (j = 0; j < source_n * QSZ; j += QSZ) {
            x = ppos0 - xp;
            y = ppos1 - yp;
            z = ppos2 - zp;

            r2 = x * x + y * y + z * z;

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
            rinv *= half; /* flips sign to avoid storing -0.5 */
            /* end Newton-Raphson */

            __asm__("prefetcht0 640(%rdi)");
            t = rinv;
            eqe = t * mass;
            rinv2 = t * t;
            phi += eqe;
            eqe *= rinv2;
            a0 += x * eqe;
            a1 += y * eqe;
            a2 += z * eqe;

            t *= R * rinv2;
            eq0 = t * qx;
            eq1 = t * qy;
            eq2 = t * qz;
            eqe = x * eq0 + y * eq1 + z * eq2;
            phi += eqe;
            eqe *= three * rinv2;
            a0 += x * eqe - eq0;
            a1 += y * eqe - eq1;
            a2 += z * eqe - eq2;

            __asm__("prefetcht0 704(%rdi)");
            t *= three * rinv2 * R;
            eq0 = t * (qxx * x + qxy * y + qxz * z);
            eq1 = t * (qyy * y + qxy * x + qyz * z);
            eq2 = t * (-(qxx + qyy) * z + qxz * x + qyz * y);
            eqe = half * (eq0 * x + eq1 * y + eq2 * z);
            phi += eqe;
            eqe *= five * rinv2;
            a0 += x * eqe - eq0;
            a1 += y * eqe - eq1;
            a2 += z * eqe - eq2;
        }
        accp[i + 0] += vsf_hsum(a0);
        accp[i + 1] += vsf_hsum(a1);
        accp[i + 2] += vsf_hsum(a2);
        accp[i + 3] += vsf_hsum(phi);
    }
}

#define PREFETCH_AHEAD 128

void pMinteract(
    const float *p, float *accp, const int m, const int stride, const vsf *f, const int n) {
    int i, j;
    vsf a0, a1, a2, phi;

    for (i = 0; i < m * stride; i += stride) {
        vsf ppos0 = (vsf)vsf_scalar(p[i + 1]);
        vsf ppos1 = (vsf)vsf_scalar(p[i + 2]);
        vsf ppos2 = (vsf)vsf_scalar(p[i + 3]);
        a0 = a1 = a2 = phi = (vsf)vsf_scalar(0.0f);

        for (j = 0; j < n * MSZ; j += MSZ) {
            vsf x = ppos0 - xp;
            vsf y = ppos1 - yp;
            vsf z = ppos2 - zp;
            vsf r2 = x * x + y * y + z * z;
            vsf rinv = vsf_rsqrt(r2); /* 11 bit precision */

            __builtin_prefetch(&f[j + PREFETCH_AHEAD / MSZ], 0, 3);
            /* Newton-Raphson */
            vsf t = rinv;
            r2 *= rinv;
            rinv *= r2;
            rinv -= (vsf)vsf_scalar(3.0f);
            rinv *= t;
            rinv *= (vsf)vsf_scalar(0.5f); /* flips sign to avoid storing -0.5 */
            /* end Newton-Raphson */

            vsf eqe = rinv * mass;
            vsf rinv2 = rinv * rinv;
            phi += eqe;
            eqe *= rinv2;
            a0 += x * eqe;
            a1 += y * eqe;
            a2 += z * eqe;
        }
        accp[i + 0] += vsf_hsum(a0);
        accp[i + 1] += vsf_hsum(a1);
        accp[i + 2] += vsf_hsum(a2);
        accp[i + 3] += vsf_hsum(phi);
    }
}

typedef struct {
    vsf phi;
    vsf a;
} sret_s;

/* Dehnen K1 Compensating Kernel */
static sret_s sK1(const vsf r2, const float e) {
    sret_s ret;
    const vsf one = vsf_scalar(1.0f);
    const vsf f1 = vsf_scalar(-3.0f / 2.0f);
    const vsf f2 = vsf_scalar(135.0f / 16.0f);
    const vsf p1 = vsf_scalar(-1.0f / 2.0f);
    const vsf p2 = vsf_scalar(3.0f / 8.0f);
    const vsf p3 = vsf_scalar(-45.0f / 32.0f);
    const vsf eps_inv = vsf_scalar(e);

    vsf u2 = r2 * eps_inv * eps_inv;
    vsf mask = vsf_cmple(u2, one);
    u2 -= one;

    vsf t = p3 * u2;
    t += p2;
    t *= u2;
    t += p1;
    t *= u2;
    t += one;
    t *= eps_inv;
    ret.phi = -vsf_and(mask, t);

    t = f2 * u2;
    t += f1;
    t *= u2;
    t += one;
    t *= eps_inv * eps_inv * eps_inv;
    ret.a = -vsf_and(mask, t);

    return ret;
}

void pMinteract_sK1(const float *p,
                    float *accp,
                    const int m,
                    const int stride,
                    const vsf *f,
                    const int n,
                    const float e,
                    int *ncut) {
    int i, j;
    vsf a0, a1, a2, phi;
    const vsf eps2 = vsf_scalar(1.0f / (e * e));

    for (i = 0; i < m * stride; i += stride) {
        vsf ppos0 = (vsf)vsf_scalar(p[i + 1]);
        vsf ppos1 = (vsf)vsf_scalar(p[i + 2]);
        vsf ppos2 = (vsf)vsf_scalar(p[i + 3]);
        a0 = a1 = a2 = phi = (vsf)vsf_scalar(0.0f);

        for (j = 0; j < n * MSZ; j += MSZ) {
            vsf x = ppos0 - xp;
            vsf y = ppos1 - yp;
            vsf z = ppos2 - zp;
            vsf r2 = x * x + y * y + z * z;
            vsf rinv = vsf_rsqrt(r2); /* 11 bit precision */
            vsf mask = vsf_cmple(eps2, r2);
            unsigned int nmask = vsf_count(mask);
            if (nmask != NSSE) {
                *ncut += NSSE - nmask;
                sret_s ret = sK1(r2, e);
                phi += mass * ret.phi;
                ret.a *= mass;
                a0 += x * ret.a;
                a1 += y * ret.a;
                a2 += z * ret.a;
            }

            __builtin_prefetch(&f[j + PREFETCH_AHEAD / MSZ], 0, 3);
            /* Newton-Raphson */
            vsf t = rinv;
            r2 *= rinv;
            rinv *= r2;
            rinv -= (vsf)vsf_scalar(3.0f);
            rinv *= t;
            rinv *= (vsf)vsf_scalar(0.5f); /* flips sign to avoid storing -0.5 */
            /* end Newton-Raphson */
            rinv = vsf_and(mask, rinv);

            vsf eqe = rinv * mass;
            vsf rinv2 = rinv * rinv;
            phi += eqe;
            eqe *= rinv2;
            a0 += x * eqe;
            a1 += y * eqe;
            a2 += z * eqe;
        }
        accp[i + 0] += vsf_hsum(a0);
        accp[i + 1] += vsf_hsum(a1);
        accp[i + 2] += vsf_hsum(a2);
        accp[i + 3] += vsf_hsum(phi);
    }
}
