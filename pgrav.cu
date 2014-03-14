/*
 * Copyright 2013-2014 Michael S. Warren.
 * All Rights Reserved.
 */

#include <stdio.h>
#include <cuda.h>
#include <assert.h>
#include "cudavec.h"
#include "segment.h"
#include "gravcuda.h"

#define mass f[ii+0*NSSE]
#define xp f[ii+1*NSSE]
#define yp f[ii+2*NSSE]
#define zp f[ii+3*NSSE]
#define MSZ 4
#define R  f[ii+4*NSSE]
#define qx f[ii+5*NSSE]
#define qy f[ii+6*NSSE]
#define qz f[ii+7*NSSE]
#define qxx f[ii+8*NSSE]
#define qxy f[ii+9*NSSE]
#define qyy f[ii+10*NSSE]
#define qxz f[ii+11*NSSE]
#define qyz f[ii+12*NSSE]
#define QSZ 13
#define qxxx f[ii+13*NSSE]
#define qxxy f[ii+14*NSSE]
#define qxyy f[ii+15*NSSE]
#define qyyy f[ii+16*NSSE]
#define qxxz f[ii+17*NSSE]
#define qxyz f[ii+18*NSSE]
#define qyyz f[ii+19*NSSE]
#define qxxxx f[ii+20*NSSE]
#define qxxxy f[ii+21*NSSE]
#define qxxyy f[ii+22*NSSE]
#define qxyyy f[ii+23*NSSE]
#define qyyyy f[ii+24*NSSE]
#define qxxxz f[ii+25*NSSE]
#define qxxyz f[ii+26*NSSE]
#define qxyyz f[ii+27*NSSE]
#define qyyyz f[ii+28*NSSE]
#define HSZ 29

__global__ void
pH(const float *p, float *ret, const int n, const int stride, const float *f, const int source_n)
{
    int i, ii;
    float t[VECWIDTH], r2[VECWIDTH], rinv[VECWIDTH], rinv2[VECWIDTH];
    float x[VECWIDTH], y[VECWIDTH], z[VECWIDTH];
    float xx[VECWIDTH], xy[VECWIDTH], yy[VECWIDTH], xz[VECWIDTH], yz[VECWIDTH], zz[VECWIDTH];
    float xxx[VECWIDTH], xxy[VECWIDTH], xyy[VECWIDTH], yyy[VECWIDTH], xxz[VECWIDTH], xyz[VECWIDTH], yyz[VECWIDTH];
    float eqe[VECWIDTH], eq0[VECWIDTH], eq1[VECWIDTH], eq2[VECWIDTH];
    float4 accp[VECWIDTH] = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos[VECWIDTH] = Vdecl(p, 3);

    if (VECWIDTH*index >= n) return;

    for (i = 0; i < source_n; i++) {
	ii = (i/NSSE)*NSSE*HSZ + i%NSSE;
    	VVS(x, = ppos, .x - xp);	
    	VVS(y, = ppos, .y - yp);	
    	VVS(z, = ppos, .z - zp);	
	VVV(r2, = x, * x);
	VVV(r2, += y, * y);
	VVV(r2, += z, * z);
	VVS(rinv, = -rsqrtf LPAREN r2, RPAREN);
	
	VV(t, = rinv);
	VV(eqe, = mass * t);
	VVV(rinv2, = t, * t);
	VV(accp,.Phi += eqe);
	VV(eqe, *= rinv2);
	VVV(accp,.Ax += x, * eqe);
	VVV(accp,.Ay += y, * eqe);
	VVV(accp,.Az += z, * eqe);

        VV(t, *= R * rinv2);
        VV(eq0, = qx * t);
        VV(eq1, = qy * t);
        VV(eq2, = qz * t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VV(accp,.Phi += eqe);
        VV(eqe, *= 3.0f * rinv2);
	VVVV(accp,.Ax += x, * eqe, - eq0);
	VVVV(accp,.Ay += y, * eqe, - eq1);
	VVVV(accp,.Az += z, * eqe, - eq2);

        VV(t, *= 3.0f * R * rinv2);
        VV(eq0,  = qxx * x);
        VV(eq1,  = qyy * y);
	VV(eq2, = -LPAREN qxx + qyy RPAREN * z);
        VV(eq0, += qxy * y);
	VV(eq1, += qxy * x);
	VV(eq2, += qxz * x);
	VV(eq0, += qxz * z);
	VV(eq1, += qyz * z);
	VV(eq2, += qyz * y);
	VV(eq0, *= t);
	VV(eq1, *= t);
	VV(eq2, *= t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VS(eqe, *= 0.5f);
	VV(accp,.Phi -= eqe);
	VV(eqe, *= 5.0f * rinv2);
	VVVV(accp,.Ax += x, * eqe, - eq0);
	VVVV(accp,.Ay += y, * eqe, - eq1);
	VVVV(accp,.Az += z, * eqe, - eq2);

	VV(t, *= 5.0f * R * rinv2);
        VVV(xx, = 0.5f * x, * x);
        VVV(xy, = x, * y);
        VVV(xz, = x, * z);
        VVV(yy, = 0.5f * y, * y);
        VVV(yz, = y, * z);
        VVV(zz, = 0.5f * z, * z);
	VVV(xxx, = (1.0f/3.0f) * xx, - zz);
	VV(xxx, *= x);
	VVV(xxz, = xx, - (1.0f/3.0f) * zz);
	VV(xxz, *= z);
        VVV(yyy, = (1.0f/3.0f) * yy, - zz);
	VV(yyy, *= y);
	VVV(yyz, = yy, - (1.0f/3.0f) * zz);
	VV(yyz, *= z);
        VV(xx, -= zz);
        VV(yy, -= zz);

        VV(eq0, = qxxx * xx);
        VV(eq1, = qxyy * xy);
        VV(eq2, = -LPAREN qxxx + qxyy RPAREN * xz);
        VV(eq0, += qxyy * yy);
        VV(eq1, += qxxy * xx);
        VV(eq2, -= LPAREN qxxy + qyyy RPAREN * yz);
        VV(eq0, += qxxy * xy);
        VV(eq1, += qyyy * yy);
        VV(eq2, += qxxz * xx);
        VV(eq0, += qxxz * xz);
        VV(eq1, += qyyz * yz);
        VV(eq2, += qyyz * yy);
        VV(eq0, += qxyz * yz);
        VV(eq1, += qxyz * xz);
        VV(eq2, += qxyz * xy);
        VV(eq0, *= t);
        VV(eq1, *= t);
        VV(eq2, *= t);
        VVV(eqe, = eq0, * x);
        VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VS(eqe, *= (1.0f/3.0f));
	VV(accp,.Phi += eqe);
	VV(eqe, *= 7.0f * rinv2);
	VVVV(accp,.Ax += x, * eqe, - eq0);
	VVVV(accp,.Ay += y, * eqe, - eq1);
	VVVV(accp,.Az += z, * eqe, - eq2);

        VV(t, *= 7.0f * R * rinv2);
        VVV(xxy, = y, * xx);
        VVV(xyy, = x, * yy);
        VVV(xyz, = xy, * z);
        VV(eq2, = -qxxxx * xxz);
        VV(eq0, = qxxxx * xxx);
        VV(eq1, = qxyyy * xyy);
        VV(eq2, -= LPAREN qxyyy + qxxxy RPAREN * xyz);
        VV(eq0, += qxyyy * yyy);
        VV(eq1, += qxxxy * xxx);
        VV(eq2, -= qyyyy * yyz);
        VV(eq0, += qxxxy * xxy);
        VV(eq1, += qyyyy * yyy);
        VV(eq2, += qxxxz * xxx);
        VV(eq0, += qxxxz * xxz);
        VV(eq1, += qyyyz * yyz);
        VV(eq2, += qyyyz * yyy);
        VV(eq0, += qxxyy * xyy);
        VV(eq1, += qxxyy * xxy);
	VVV(eq2, -= qxxyy * xxz, + qxxyy * yyz);
        VV(eq0, += qxxyz * xyz);
        VV(eq1, += qxxyz * xxz);
	VV(eq2, += qxxyz * xxy);
        VV(eq0, += qxyyz * yyz);
        VV(eq1, += qxyyz * xyz);
	VV(eq2, += qxyyz * xyy);
	VV(eq0, *= t);
	VV(eq1, *= t);
        VV(eq2, *= t);
        VVV(eqe,  = eq0, * x);
        VVV(eqe, += eq1, * y);
        VVV(eqe, += eq2, * z);
        VS(eqe, *= 0.25f);
	VV(accp,.Phi += eqe);
	VV(eqe, *= 9.0f * rinv2);
	VVVV(accp,.Ax += x, * eqe, - eq0);
	VVVV(accp,.Ay += y, * eqe, - eq1);
	VVVV(accp,.Az += z, * eqe, - eq2);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (VECWIDTH*index+i < n) {
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+0], accp[i].Ax);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+1], accp[i].Ay);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+2], accp[i].Az);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+3], accp[i].Phi);
	}
    }
}

__global__ void
pQ(const float *p, float *ret, const int n, const int stride, 
   const float *f, const int source_n)
{
    int i, ii;
    float t[VECWIDTH], r2[VECWIDTH], rinv[VECWIDTH], rinv2[VECWIDTH];
    float x[VECWIDTH], y[VECWIDTH], z[VECWIDTH];
    float eqe[VECWIDTH], eq0[VECWIDTH], eq1[VECWIDTH], eq2[VECWIDTH];
    float4 accp[VECWIDTH] = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos[VECWIDTH] = Vdecl(p, 3);

    if (VECWIDTH*index >= n) return;

    for (i = 0; i < source_n; i++) {
	ii = (i/NSSE)*NSSE*QSZ + i%NSSE;
    	VVS(x, = ppos, .x - xp);	
    	VVS(y, = ppos, .y - yp);	
    	VVS(z, = ppos, .z - zp);	
	VVV(r2, = x, * x);
	VVV(r2, += y, * y);
	VVV(r2, += z, * z);
	VVS(rinv, = -rsqrtf LPAREN r2, RPAREN);
	
	VV(t, = rinv);
	VV(eqe, = mass * t);
	VVV(rinv2, = t, * t);
	VV(accp,.Phi += eqe);
	VV(eqe, *= rinv2);
	VVV(accp,.Ax += x, * eqe);
	VVV(accp,.Ay += y, * eqe);
	VVV(accp,.Az += z, * eqe);

        VV(t, *= R * rinv2);
        VV(eq0, = qx * t);
        VV(eq1, = qy * t);
        VV(eq2, = qz * t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VV(accp,.Phi += eqe);
        VV(eqe, *= 3.0f * rinv2);
	VVVV(accp,.Ax += x, * eqe, - eq0);
	VVVV(accp,.Ay += y, * eqe, - eq1);
	VVVV(accp,.Az += z, * eqe, - eq2);

        VV(t, *= 3.0f * R * rinv2);
        VV(eq0,  = qxx * x);
        VV(eq1,  = qyy * y);
	VV(eq2, = -LPAREN qxx + qyy RPAREN * z);
        VV(eq0, += qxy * y);
	VV(eq1, += qxy * x);
	VV(eq2, += qxz * x);
	VV(eq0, += qxz * z);
	VV(eq1, += qyz * z);
	VV(eq2, += qyz * y);
	VV(eq0, *= t);
	VV(eq1, *= t);
	VV(eq2, *= t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VS(eqe, *= 0.5f);
	VV(accp,.Phi -= eqe);
	VV(eqe, *= 5.0f * rinv2);
	VVVV(accp,.Ax += x, * eqe, - eq0);
	VVVV(accp,.Ay += y, * eqe, - eq1);
	VVVV(accp,.Az += z, * eqe, - eq2);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (VECWIDTH*index+i < n) {
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+0], accp[i].Ax);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+1], accp[i].Ay);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+2], accp[i].Az);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+3], accp[i].Phi);
	}
    }
}

__global__ void
pQ1_mnss(const float * __restrict__ sink, const int sink_m,
	 const segment * __restrict__ seg,  const float * __restrict__ ff, float *accp_ret)
{
    float t, r2, rinv, rinv2;
    float x, y, z;
    float eqe, eq0, eq1, eq2;
    float4 accp = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos = {sink[3*index], sink[3*index+1], sink[3*index+2]};
    const float * __restrict__ f = ff + seg[index].base * QSZ;

    if (index >= sink_m) return;

    for (int i = 0; i < seg[index].length; i++) {
	const int ii = (i/NSSE)*NSSE*QSZ + i%NSSE;
    	x = ppos.x - xp;	
 	y = ppos.y - yp;	
    	z = ppos.z - zp;	
	r2 = x * x;
	r2 += y * y;
	r2 += z * z;
	rinv = -rsqrtf(r2);
	
	t = rinv;
	eqe = mass * t;
	rinv2 = t * t;
	accp.Phi += eqe;
	eqe *= rinv2;
	accp.Ax += x * eqe;
	accp.Ay += y * eqe;
	accp.Az += z * eqe;

        t *= R * rinv2;
        eq0 = qx * t;
        eq1 = qy * t;
        eq2 = qz * t;
	eqe  = eq0 * x;
	eqe += eq1 * y;
	eqe += eq2 * z;
	accp.Phi += eqe;
        eqe *= 3.0f * rinv2;
	accp.Ax += x * eqe - eq0;
	accp.Ay += y * eqe - eq1;
	accp.Az += z * eqe - eq2;

        t *= 3.0f * R * rinv2;
        eq0  = qxx * x;
        eq1  = qyy * y;
	eq2 = -(qxx + qyy) * z;
        eq0 += qxy * y;
	eq1 += qxy * x;
	eq2 += qxz * x;
	eq0 += qxz * z;
	eq1 += qyz * z;
	eq2 += qyz * y;
	eq0 *= t;
	eq1 *= t;
	eq2 *= t;
	eqe  = eq0 * x;
	eqe += eq1 * y;
	eqe += eq2 * z;
	eqe *= 0.5f;
	accp.Phi -= eqe;
	eqe *= 5.0f * rinv2;
	accp.Ax += x * eqe - eq0;
	accp.Ay += y * eqe - eq1;
	accp.Az += z * eqe - eq2;
    }
    atomicAdd(&accp_ret[4*index+0], accp.Ax);
    atomicAdd(&accp_ret[4*index+1], accp.Ay);
    atomicAdd(&accp_ret[4*index+2], accp.Az);
    atomicAdd(&accp_ret[4*index+3], accp.Phi);
}

__global__ void
pM(const float *p, float *ret, const int n, const int stride, 
   const float *f, const int source_n)
{
    int i, ii;
    float t[VECWIDTH], r2[VECWIDTH], rinv[VECWIDTH], rinv2[VECWIDTH];
    float x[VECWIDTH], y[VECWIDTH], z[VECWIDTH];
    float eqe[VECWIDTH];
    float4 accp[VECWIDTH] = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos[VECWIDTH] = Vdecl(p, 3);

    if (VECWIDTH*index >= n) return;

    for (i = 0; i < source_n; i++) {
	ii = (i/NSSE)*NSSE*MSZ + i%NSSE;
    	VVS(x, = ppos, .x - xp);	
    	VVS(y, = ppos, .y - yp);	
    	VVS(z, = ppos, .z - zp);	
	VVV(r2, = x, * x);
	VVV(r2, += y, * y);
	VVV(r2, += z, * z);
	VVS(rinv, = -rsqrtf LPAREN r2, RPAREN);
	
	VV(t, = rinv);
	VV(eqe, = mass * t);
	VVV(rinv2, = t, * t);
	VV(accp,.Phi += eqe);
	VV(eqe, *= rinv2);
	VVV(accp,.Ax += x, * eqe);
	VVV(accp,.Ay += y, * eqe);
	VVV(accp,.Az += z, * eqe);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (VECWIDTH*index+i < n) {
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+0], accp[i].Ax);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+1], accp[i].Ay);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+2], accp[i].Az);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+3], accp[i].Phi);
	}
    }
}

__global__ void
pM_sK1(const float *p, float *ret, const int n, const int stride, 
       const float *f, const int source_n, const float eps_inv, int *ncut)
{
    int i, ii;
    float t[VECWIDTH], r2[VECWIDTH], rinv[VECWIDTH], u2[VECWIDTH];
    float x[VECWIDTH], y[VECWIDTH], z[VECWIDTH];
    float eqe[VECWIDTH];
    float4 accp[VECWIDTH] = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos[VECWIDTH] = Vdecl(p, 3);
    const float eps_inv2 = eps_inv*eps_inv;
    const float eps2 = 1.0f/eps_inv2;

    if (VECWIDTH*index >= n) return;

    for (i = 0; i < source_n; i++) {
	ii = (i/NSSE)*NSSE*MSZ + i%NSSE;
    	VVS(x, = ppos, .x - xp);	
    	VVS(y, = ppos, .y - yp);	
    	VVS(z, = ppos, .z - zp);	
	VVV(r2, = x, * x);
	VVV(r2, += y, * y);
	VVV(r2, += z, * z);
	VVS(rinv, = -rsqrtf LPAREN r2, RPAREN);
	for (int j = 0; j < VECWIDTH; j++) {
	    if (r2[j] <= eps2) {
		u2[j] = r2[j] * eps_inv2 - 1.0f;
		t[j] = (((45.0f/32.0f)*u2[j] + (-3.0f/8.0f)) * u2[j] + (1.0f/2.0f)) * u2[j] - 1.0f;
		t[j] *= eps_inv;
		eqe[j] = ((-135.0f/16.0f) * u2[j] + (3.0f/2.0f)) * u2[j] - 1.0f;
		eqe[j] *= mass * eps_inv * eps_inv2;
	    } else {
		t[j] = rinv[j];
		eqe[j] = mass * rinv[j] * rinv[j] * rinv[j];
	    }
	}
	VV(accp,.Phi += mass * t);
	VVV(accp,.Ax += x, * eqe);
	VVV(accp,.Ay += y, * eqe);
	VVV(accp,.Az += z, * eqe);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (VECWIDTH*index+i < n) {
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+0], accp[i].Ax);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+1], accp[i].Ay);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+2], accp[i].Az);
	    atomicAdd(&ret[stride*(index*VECWIDTH+i)+3], accp[i].Phi);
	}
    }
}

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

#define mass f[i*QSZ+0]
#define xp f[i*QSZ+1]
#define yp f[i*QSZ+2]
#define zp f[i*QSZ+3]
#define R  f[i*QSZ+4]
#define qx f[i*QSZ+5]
#define qy f[i*QSZ+6]
#define qz f[i*QSZ+7]
#define qxx f[i*QSZ+8]
#define qxy f[i*QSZ+9]
#define qyy f[i*QSZ+10]
#define qxz f[i*QSZ+11]
#define qyz f[i*QSZ+12]

__global__ void
pQQ1(const float * __restrict__ sink, const int sink_m,
     const int * __restrict__ ii, const int source_n,
     const float * __restrict__ f, float *accp_ret)
{
    float t, r2, rinv, rinv2;
    float x, y, z;
    float eqe, eq0, eq1, eq2;
    float4 accp = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos = {sink[3*index], sink[3*index+1], sink[3*index+2]};

    if (index >= source_n) return;

    for (int i = 0; i < source_n; i++) {
    	x = ppos.x - xp;	
 	y = ppos.y - yp;	
    	z = ppos.z - zp;	
	r2 = x * x;
	r2 += y * y;
	r2 += z * z;
	rinv = -rsqrtf(r2);
	
	t = rinv;
	eqe = mass * t;
	rinv2 = t * t;
	accp.Phi += eqe;
	eqe *= rinv2;
	accp.Ax += x * eqe;
	accp.Ay += y * eqe;
	accp.Az += z * eqe;

        t *= R * rinv2;
        eq0 = qx * t;
        eq1 = qy * t;
        eq2 = qz * t;
	eqe  = eq0 * x;
	eqe += eq1 * y;
	eqe += eq2 * z;
	accp.Phi += eqe;
        eqe *= 3.0f * rinv2;
	accp.Ax += x * eqe - eq0;
	accp.Ay += y * eqe - eq1;
	accp.Az += z * eqe - eq2;

        t *= 3.0f * R * rinv2;
        eq0  = qxx * x;
        eq1  = qyy * y;
	eq2 = -(qxx + qyy) * z;
        eq0 += qxy * y;
	eq1 += qxy * x;
	eq2 += qxz * x;
	eq0 += qxz * z;
	eq1 += qyz * z;
	eq2 += qyz * y;
	eq0 *= t;
	eq1 *= t;
	eq2 *= t;
	eqe  = eq0 * x;
	eqe += eq1 * y;
	eqe += eq2 * z;
	eqe *= 0.5f;
	accp.Phi -= eqe;
	eqe *= 5.0f * rinv2;
	accp.Ax += x * eqe - eq0;
	accp.Ay += y * eqe - eq1;
	accp.Az += z * eqe - eq2;
    }
    atomicAdd(&accp_ret[4*index+0], accp.Ax);
    atomicAdd(&accp_ret[4*index+1], accp.Ay);
    atomicAdd(&accp_ret[4*index+2], accp.Az);
    atomicAdd(&accp_ret[4*index+3], accp.Phi);
}


__global__ void
pMM_sK1(const float * __restrict__ sink, const int sink_m,
	const segment * __restrict__ seg, const int seg_len, 
	const float * __restrict__ source, const int source_n, 
	const float mmass, const float eps_inv, float *accp, int *ncut)
{
    float t[VECWIDTH], r2[VECWIDTH], rinv[VECWIDTH], u2[VECWIDTH];
    float x[VECWIDTH], y[VECWIDTH], z[VECWIDTH];
    float eqe[VECWIDTH];
    float4 acc[VECWIDTH] = {};
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    const float3 ppos[VECWIDTH] = Vdecl(sink, 3);
    const float eps_inv2 = eps_inv*eps_inv;
    const float eps2 = 1.0f/eps_inv2;
    int sum_n = 0;

    if (VECWIDTH*index >= sink_m) return;

    for (int seg_i = 0; seg_i < seg_len; seg_i++) {
	sum_n += seg[seg_i].length;
	for (int seg_j = 0; seg_j < seg[seg_i].length; seg_j++) {
	    const float3 qpos = ((float3 *)source)[seg[seg_i].base + seg_j];
	    for (int j = 0; j < VECWIDTH; j++) {
		x[j] = ppos[j].x - qpos.x;
		y[j] = ppos[j].y - qpos.y;
		z[j] = ppos[j].z - qpos.z;
		r2[j] = x[j]*x[j] + y[j]*y[j] + z[j]*z[j];
		if (r2[j] > eps2) {
		    rinv[j] = -rsqrtf(r2[j]);
		    eqe[j] = mmass * rinv[j];
		    acc[j].Phi += eqe[j];
		    eqe[j] *= rinv[j] * rinv[j];
		} else {
		    u2[j] = r2[j] * eps_inv2 - 1.0f;
		    t[j] = (((45.0f/32.0f)*u2[j] + (-3.0f/8.0f)) * u2[j] + (1.0f/2.0f)) * u2[j] - 1.0f;
		    t[j] *= mmass * eps_inv;
		    acc[j].Phi += t[j];
		    eqe[j] = ((-135.0f/16.0f) * u2[j] + (3.0f/2.0f)) * u2[j] - 1.0f;
		    eqe[j] *= mmass * eps_inv * eps_inv2;
		}
		acc[j].Ax += x[j] * eqe[j];
		acc[j].Ay += y[j] * eqe[j];
		acc[j].Az += z[j] * eqe[j];
	    }
	}
    }
    for (int j = 0; j < VECWIDTH; j++) {
	if (VECWIDTH*index+j < sink_m) {
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+0], acc[j].Ax);
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+1], acc[j].Ay);
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+2], acc[j].Az);
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+3], acc[j].Phi);
	}
    }
    assert(sum_n == source_n);
}


__global__ void
pMM_mnss_sK1(const float * __restrict__ sink, const int sink_m,
	     const uint16_t * __restrict__ ss_index, const segment * __restrict__ ss_seg, const int ss_len,
	     const segment * __restrict__ s_seg, const int s_seg_len, 
	     const float * __restrict__ source, const int *s_source_n, 
	     const float mmass, const float eps_inv, float *accp, int *ncut)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (VECWIDTH*index >= sink_m) return;
    const int sindex = ss_index[VECWIDTH*index];
    assert(sindex < ss_len);
    assert(ss_seg[sindex].base + ss_seg[sindex].length <= s_seg_len);
    const segment *seg = s_seg + ss_seg[sindex].base;
    const int seg_len = ss_seg[sindex].length;
    const int source_n = s_source_n[sindex];

    float t[VECWIDTH], r2[VECWIDTH], rinv[VECWIDTH], u2[VECWIDTH];
    float x[VECWIDTH], y[VECWIDTH], z[VECWIDTH];
    float eqe[VECWIDTH];
    float4 acc[VECWIDTH] = {};
    const float3 ppos[VECWIDTH] = Vdecl(sink, 3);
    const float eps_inv2 = eps_inv*eps_inv;
    const float eps2 = 1.0f/eps_inv2;
    int source_interactions = 0;

    for (int seg_i = 0; seg_i < seg_len; seg_i++) {
	source_interactions += seg[seg_i].length;
	for (int seg_j = 0; seg_j < seg[seg_i].length; seg_j++) {
	    const float3 qpos = ((float3 *)source)[seg[seg_i].base + seg_j];
	    for (int j = 0; j < VECWIDTH; j++) {
		x[j] = ppos[j].x - qpos.x;
		y[j] = ppos[j].y - qpos.y;
		z[j] = ppos[j].z - qpos.z;
		r2[j] = x[j]*x[j] + y[j]*y[j] + z[j]*z[j];
		if (r2[j] > eps2) {
		    rinv[j] = -rsqrtf(r2[j]);
		    eqe[j] = mmass * rinv[j];
		    acc[j].Phi += eqe[j];
		    eqe[j] *= rinv[j] * rinv[j];
		} else {
		    u2[j] = r2[j] * eps_inv2 - 1.0f;
		    t[j] = (((45.0f/32.0f)*u2[j] + (-3.0f/8.0f)) * u2[j] + (1.0f/2.0f)) * u2[j] - 1.0f;
		    t[j] *= mmass * eps_inv;
		    acc[j].Phi += t[j];
		    eqe[j] = ((-135.0f/16.0f) * u2[j] + (3.0f/2.0f)) * u2[j] - 1.0f;
		    eqe[j] *= mmass * eps_inv * eps_inv2;
		}
		acc[j].Ax += x[j] * eqe[j];
		acc[j].Ay += y[j] * eqe[j];
		acc[j].Az += z[j] * eqe[j];
	    }
	}
    }
    for (int j = 0; j < VECWIDTH; j++) {
	if (VECWIDTH*index+j < sink_m) {
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+0], acc[j].Ax);
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+1], acc[j].Ay);
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+2], acc[j].Az);
	    atomicAdd(&accp[4*(index*VECWIDTH+j)+3], acc[j].Phi);
	}
    }
    assert(source_interactions == source_n);
}

__global__ void
pMM1_mnss_sK1(const float * __restrict__ sink, const int sink_m,
	      const uint16_t * __restrict__ ss_index, const segment * __restrict__ ss_seg, const int ss_len,
	      const segment * __restrict__ s_seg, const int s_seg_len, 
	      const float * __restrict__ source, const int *s_source_n, 
	      const float mmass, const float eps_inv, float *accp, int *ncut, int proc)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index >= sink_m) return;
    const int sindex = ss_index[index];
    assert(sindex < ss_len);
    assert(ss_seg[sindex].base + ss_seg[sindex].length <= s_seg_len);
    const segment *seg = s_seg + ss_seg[sindex].base;
    const int seg_len = ss_seg[sindex].length;
    const int source_n = s_source_n[sindex];
    // if (debug) printf("%d %d %d %d %d\n", index, sindex, ss_len, ss_seg[sindex].base, ss_seg[sindex].length);

    const float3 ppos = {sink[3*index], sink[3*index+1], sink[3*index+2]};
    const float eps_inv2 = eps_inv*eps_inv;
    const float eps2 = 1.0f/eps_inv2;
    float4 acc = {};
    int source_interactions = 0;

    for (int seg_i = 0; seg_i < seg_len; seg_i++) {
	source_interactions += seg[seg_i].length;
	for (int seg_j = 0; seg_j < seg[seg_i].length; seg_j++) {
	    const float3 qpos = ((float3 *)source)[seg[seg_i].base + seg_j];
	    float x = ppos.x - qpos.x;
	    float y = ppos.y - qpos.y;
	    float z = ppos.z - qpos.z;
	    float r2 = x*x + y*y + z*z;
	    float eqe;
	    if (r2 > eps2) {
		float rinv = -rsqrtf(r2);
		eqe = mmass * rinv;
		acc.Phi += eqe;
		eqe *= rinv * rinv;
	    } else {
		float u2 = r2 * eps_inv2 - 1.0f;
		float t = (((45.0f/32.0f)*u2 + (-3.0f/8.0f)) * u2 + (1.0f/2.0f)) * u2 - 1.0f;
		t *= mmass * eps_inv;
		acc.Phi += t;
		eqe = ((-135.0f/16.0f) * u2 + (3.0f/2.0f)) * u2 - 1.0f;
		eqe *= mmass * eps_inv * eps_inv2;
	    }
	    acc.Ax += x * eqe;
	    acc.Ay += y * eqe;
	    acc.Az += z * eqe;
	}
    }
    if (index < sink_m) {
	atomicAdd(&accp[4*index+0], acc.Ax);
	atomicAdd(&accp[4*index+1], acc.Ay);
	atomicAdd(&accp[4*index+2], acc.Az);
	atomicAdd(&accp[4*index+3], acc.Phi);
    }
    if (source_interactions != source_n) 
	printf("CUDA ERROR proc %d: %d %d %d %d %d %d %d %d\n", 
	       proc, index, sindex, ss_len, ss_seg[sindex].base, seg_len, s_seg_len, source_n, source_interactions);
}


#include <stdio.h>
#include <unistd.h>
#include <cuda_profiler_api.h>
#include "error.h"
#include "Msgs.h"

typedef struct cudaq_t {
    int inuse;
    cudaStream_t stream;
    size_t dev_size;
    void *host;
    void *dev;
} cudaq_t;

#define NQ 32
static cudaq_t cudaq[NQ];

static int
check_cudaq(int *inuse)
{
    int i;
    cudaError_t err;
    int ninuse = 0;

    for (i = 0; i < NQ; i++) {
	if (cudaq[i].inuse) {
	    err = cudaStreamQuery(cudaq[i].stream);
	    if (err == cudaSuccess) {
		err = cudaStreamDestroy(cudaq[i].stream);
		if (err != cudaSuccess) 
		    Error("cudaStreamDestroy failed, %d %s\n", err, cudaGetErrorString(err));
		cudaq[i].inuse = 0;
	    } else if (err != cudaErrorNotReady) {
		Error("cudaStreamQuery failed, %d %s\n", err, cudaGetErrorString(err));
	    } else ninuse++;
	}
    }
    if (inuse) *inuse = ninuse;
    if (ninuse == NQ) return -1; // All slots in use
    
    for (i = 0; i < NQ; i++) {
	if (!cudaq[i].inuse) break;
    }

    if (cudaq[i].dev_size == 0) {
	cudaq[i].dev_size = 1024*1024;
	err = cudaMallocHost((void **)&cudaq[i].host, cudaq[i].dev_size, cudaHostAllocDefault);
	if (err != cudaSuccess) 
	    Error("cudaMallocHost failed, %d %s\n", err, cudaGetErrorString(err));
	err = cudaMalloc((void **)&cudaq[i].dev, cudaq[i].dev_size);
	if (err != cudaSuccess) 
	    Error("cudaMalloc failed, %d %s\n", err, cudaGetErrorString(err));
    }
    return i;
}

static void
wait_cudaq(void)
{
    int i;
    cudaError_t err;

    for (i = 0; i < NQ; i++) {
	if (cudaq[i].inuse) {
	    err = cudaStreamSynchronize(cudaq[i].stream);
	    if (err == cudaSuccess) {
		err = cudaStreamDestroy(cudaq[i].stream);
		if (err != cudaSuccess) 
		    Error("cudaStreamDestroy failed, %d %s\n", err, cudaGetErrorString(err));
		cudaq[i].inuse = 0;
	    } else {
		Error("cudaStreamSynchronize failed, %d %s\n", err, cudaGetErrorString(err));
	    }
	}
    }
}

static float *devpos, *devaccp;
static float *Btab;
static float *qdev;

extern "C" void
WalkInitSrcCUDA(float *qtab, int stride, int64_t ncells)
{
    Msgf(("WalkInitSrcCUDA stride %d nobj %ld\n", stride, ncells));

    cudaError_t err = cudaMalloc((void **)&qdev, ncells*QSZ*sizeof(float));
    if (err != cudaSuccess) 
	Error("cudaMalloc failed, %d %s\n", err, cudaGetErrorString(err));

    float *host;		// pinned
    err = cudaMallocHost((void **)&host, ncells*QSZ*sizeof(float), cudaHostAllocDefault);
    if (err != cudaSuccess) 
	Error("cudaMallocHost failed, %d %s\n", err, cudaGetErrorString(err));

    for (int64_t i = 0; i < ncells; i++) {
	host[i*QSZ+0] = qtab[i*stride+0];
	host[i*QSZ+1] = qtab[i*stride+1];
	host[i*QSZ+2] = qtab[i*stride+2];
	host[i*QSZ+3] = qtab[i*stride+3];
	host[i*QSZ+5] = qtab[i*stride+4];
	host[i*QSZ+12] = qtab[i*stride+5];
	host[i*QSZ+13] = qtab[i*stride+6];
	host[i*QSZ+14] = qtab[i*stride+7];
	host[i*QSZ+15] = qtab[i*stride+8];
	host[i*QSZ+16] = qtab[i*stride+9];
	host[i*QSZ+17] = qtab[i*stride+10];
	host[i*QSZ+18] = qtab[i*stride+11];
	host[i*QSZ+19] = qtab[i*stride+12];
    }

    err = cudaMemcpy(qdev, host, ncells*QSZ*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaFreeHost(host);
    if (err != cudaSuccess) 
	Error("cudaFreeHost failed, %d %s\n", err, cudaGetErrorString(err));
}

extern "C" void
WalkTerminateSrcCUDA(void)
{
    cudaError_t err = cudaFree(qdev);
    if (err != cudaSuccess) 
	Error("cudaFree failed, %d %s\n", err, cudaGetErrorString(err));
}

extern "C" void
WalkInitSinkCUDA(float *btab, int stride, int64_t nobj)
{
    float *hostpos;		// pinned
    size_t pos_bytes = 3 * sizeof(float) * nobj;
    size_t accp_bytes = 4 * sizeof(float) * nobj;
    cudaError_t err;

    Msgf(("WalkInitSinkCUDA stride %d nobj %ld sizeof(segment) %ld\n", 
	  stride, nobj, sizeof(segment)));
    // cudaProfilerInitialize("profile.in", "profile.out", cudaKeyValuePair);
    // cudaProfilerStart();

    Btab = btab;

    /* Copy x, y, z to device through host buffer */
    err = cudaMalloc((void **)&devpos, pos_bytes);
    if (err != cudaSuccess) 
	Error("cudaMalloc failed, %d %s\n", err, cudaGetErrorString(err));

    /* Allocate ax, ay, az, phi and initialize to zero */
    err = cudaMalloc((void **)&devaccp, accp_bytes);
    if (err != cudaSuccess) 
	Error("cudaMalloc failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaMallocHost((void **)&hostpos, pos_bytes, cudaHostAllocDefault);
    if (err != cudaSuccess) 
	Error("cudaMallocHost failed, %d %s\n", err, cudaGetErrorString(err));

    for (int64_t i = 0; i < nobj; i++) {
	hostpos[i*3+0] = btab[i*stride+1];
	hostpos[i*3+1] = btab[i*stride+2];
	hostpos[i*3+2] = btab[i*stride+3];
    }

    err = cudaMemcpy(devpos, hostpos, pos_bytes, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaFreeHost(hostpos);
    if (err != cudaSuccess) 
	Error("cudaFreeHost failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaMemset(devaccp, 0, accp_bytes);
    if (err != cudaSuccess) 
	Error("cudaMemset failed, %d %s\n", err, cudaGetErrorString(err));
}

extern "C" void
WalkTerminateSinkCUDA(float *btab, int stride, int64_t nobj)
{
    float *hostaccp;		// pinned
    size_t accp_bytes = 4 * sizeof(float) * nobj;
    cudaError_t err;
    int i;

    wait_cudaq();

    for (i = 0; i < NQ; i++) {
	if (cudaq[i].dev_size) {
	    cudaq[i].dev_size = 0;
	    err = cudaFreeHost(cudaq[i].host);
	    if (err != cudaSuccess) 
		Error("cudaFreeHost failed, %d %s\n", err, cudaGetErrorString(err));
	    err = cudaFree(cudaq[i].dev);
	    if (err != cudaSuccess) 
		Error("cudaFree failed, %d %s\n", err, cudaGetErrorString(err));
	}
    }

    err = cudaFree(devpos);
    if (err != cudaSuccess) 
	Error("cudaFree failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaMallocHost((void **)&hostaccp, accp_bytes, cudaHostAllocDefault);
    if (err != cudaSuccess) 
	Error("cudaMallocHost failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaMemcpy(hostaccp, devaccp, accp_bytes, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaFree(devaccp);
    if (err != cudaSuccess) 
	Error("cudaFree failed, %d %s\n", err, cudaGetErrorString(err));

    /* Add accp to btab */
    for (int64_t i = 0; i < nobj; i++) {
	btab[i*stride+7] += hostaccp[i*4+0];
	btab[i*stride+8] += hostaccp[i*4+1];
	btab[i*stride+9] += hostaccp[i*4+2];
	btab[i*stride+10] += hostaccp[i*4+3];
	// Msgf(("b %5ld %12g %12g %12g\n", i, btab[i*stride+7], btab[i*stride+8], btab[i*stride+9]));
    }

    err = cudaFreeHost(hostaccp);
    if (err != cudaSuccess) 
	Error("cudaFreeHost failed, %d %s\n", err, cudaGetErrorString(err));

    // cudaProfilerStop();
}

extern "C" int
qallocCUDA(int *inuse)
{
    return check_cudaq(inuse);
}

extern "C" void
grav_mn_CUDA(const char *routine, const float *p, float *accp, const int m, const int stride, 
	     const float *f, const int source_n, const int sz, float e, int *ncut, int q)
{
    cudaError_t err;

    cudaq[q].inuse = 1;
    err = cudaStreamCreate(&cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaStreamCreate failed, %d %s\n", err, cudaGetErrorString(err));

    float *dev = (float *)cudaq[q].dev;
    size_t total_size = source_n * sz * sizeof(float);

    if (cudaq[q].dev_size < total_size) {
	Error("dev_size (%ld) too small, source_n is %d, total_size is %ld\n", 
	      cudaq[q].dev_size, source_n, total_size);
    }
    memcpy(cudaq[q].host, f, total_size);
    err = cudaMemcpyAsync(dev, cudaq[q].host, total_size, 
			  cudaMemcpyHostToDevice, cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    float *pos = devpos + 3*(p-Btab)/stride;
    float *accp_ret = devaccp + 4*(p-Btab)/stride;;

    int threads = 64;
    int blocks = 1 + (m-1) / (VECWIDTH * threads);
    cudaStream_t stream = cudaq[q].stream;

    Msgf(("%s %d sinks, %d sources, %d blocks, %d threads. Offset %ld\n",
	  routine, m, source_n, blocks, threads, (p-Btab)/stride));
    if (strcmp(routine, "pM") == 0)
	pM<<<blocks,threads,0,stream>>>(pos, accp_ret, m, 4, dev, source_n);
    else if (strcmp(routine, "pQ") == 0)
	pQ<<<blocks,threads,0,stream>>>(pos, accp_ret, m, 4, dev, source_n);
    else if (strcmp(routine, "pH") == 0)
	pH<<<blocks,threads,0,stream>>>(pos, accp_ret, m, 4, dev, source_n);
    else if (strcmp(routine, "pM_sK1") == 0)
	pM_sK1<<<blocks,threads,0,stream>>>(pos, accp_ret, m, 4, dev, source_n, e, ncut);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
	Error("CUDA error, %d %s\n", err, cudaGetErrorString(err));
    }
}

extern "C" void
grav_mns_CUDA(const char *routine, const int base, const int m,
	      const segment *seg, const int seg_n, const int source_n, 
	      const float mmass, float e, int *ncut, int q)
{
    cudaError_t err;

    cudaq[q].inuse = 1;
    err = cudaStreamCreate(&cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaStreamCreate failed, %d %s\n", err, cudaGetErrorString(err));
    
    segment *devseg = (segment *)cudaq[q].dev;
    size_t total_size = seg_n * sizeof(segment);

    if (cudaq[q].dev_size < total_size) {
	Error("dev_size (%ld) too small, seg_n is %d, total_size is %ld\n", 
	      cudaq[q].dev_size, seg_n, total_size);
    }
    memcpy(cudaq[q].host, seg, total_size);
    err = cudaMemcpyAsync(devseg, cudaq[q].host, total_size, 
			  cudaMemcpyHostToDevice, cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    float *sink = devpos + base * 3;
    float *source = devpos;
    float *accp = devaccp + base * 4;

    int threads = 64;
    int blocks = 1 + (m-1) / (VECWIDTH * threads);
    cudaStream_t stream = cudaq[q].stream;

    Msgf(("M %d sinks, %d segments %d sources, %d blocks, %d threads. Base %d\n",
	  m, seg_n, source_n, blocks, threads, base));
    if (strcmp(routine, "pMM_sK1") == 0)
	pMM_sK1<<<blocks,threads,0,stream>>>
	    (sink, m, devseg, seg_n, source, source_n, mmass, e, accp, ncut);
    else Error("routine %s not found\n", routine);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
	Error("CUDA error, %d %s\n", err, cudaGetErrorString(err));
    }
}

extern "C" void
grav_qns_CUDA(const char *routine, const int base, const int m,
	      const int *source_list, const int source_n, int q)
{
    cudaError_t err;

    cudaq[q].inuse = 1;
    err = cudaStreamCreate(&cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaStreamCreate failed, %d %s\n", err, cudaGetErrorString(err));
    
    int *dev = (int *)cudaq[q].dev;
    size_t total_size = source_n * sizeof(int);

    if (cudaq[q].dev_size < total_size) {
	Error("dev_size (%ld) too small, source_n is %d, total_size is %ld\n", 
	      cudaq[q].dev_size, source_n, total_size);
    }
    memcpy(cudaq[q].host, source_list, total_size);
    err = cudaMemcpyAsync(dev, cudaq[q].host, total_size, 
			  cudaMemcpyHostToDevice, cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    float *sink = devpos + base * 3;
    int *source = dev;
    float *accp = devaccp + base * 4;

    int threads = 256;
    int blocks = 1 + (m-1) / threads;
    cudaStream_t stream = cudaq[q].stream;

    Msgf(("QQ %d sinks, %d sources, %d blocks, %d threads. Base %d\n",
	  m, source_n, blocks, threads, base));
    if (strcmp(routine, "pQQ1") == 0)
	pQQ1<<<blocks,threads,0,stream>>>
	    (sink, m, source, source_n, qdev, accp);
    else Error("routine %s not found\n", routine);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
	Error("CUDA error, %d %s\n", err, cudaGetErrorString(err));
    }
}

extern "C" void
grav_mnss_CUDA(const char *routine, const int sink_base, const int m,
	       const uint16_t *ss_index, const segment *ss_seg, const int ss_len,
	       const segment *seg, const int seg_len, 
	       const int source_base, const int *source_n, 
	       const float mmass, float e, int *ncut, int q)
{
    int align_mask = 127;	// align device pointers to 128 byte boundary
    cudaError_t err;

    cudaq[q].inuse = 1;
    err = cudaStreamCreate(&cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaStreamCreate failed, %d %s\n", err, cudaGetErrorString(err));

    char *dev = (char *)cudaq[q].dev;

    size_t ss_index_offset = 0;
    size_t ss_index_size = m * sizeof(uint16_t);
    uint16_t *ss_index_dev = (uint16_t *)(dev + ss_index_offset);

    size_t ss_seg_offset = ss_index_size + align_mask & ~align_mask;
    size_t ss_seg_size = ss_len * sizeof(segment);
    segment *ss_seg_dev = (segment *)(dev + ss_seg_offset);

    size_t source_n_offset = ss_seg_offset + ss_seg_size + align_mask & ~align_mask;
    size_t source_n_size = ss_len * sizeof(int);
    int *source_n_dev = (int *)(dev + source_n_offset);

    size_t seg_offset = source_n_offset + source_n_size + align_mask & ~align_mask;
    size_t seg_size = seg_len * sizeof(segment);
    segment *seg_dev = (segment *)(dev + seg_offset);

    size_t total_size = seg_offset + seg_size;

    if (cudaq[q].dev_size < total_size) {
	Error("dev_size (%ld) too small, total size is %ld\n", 
	      cudaq[q].dev_size, total_size);
    }
    memcpy(cudaq[q].host, ss_index, ss_index_size);
    memcpy((char *)cudaq[q].host + ss_seg_offset, ss_seg, ss_seg_size);
    memcpy((char *)cudaq[q].host + source_n_offset, source_n, source_n_size);
    memcpy((char *)cudaq[q].host + seg_offset, seg, seg_size);
    err = cudaMemcpyAsync(dev, cudaq[q].host, total_size, 
			  cudaMemcpyHostToDevice, cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    float *sink = devpos + sink_base * 3;
    float *source = devpos + source_base * 3;
    float *accp = devaccp + sink_base * 4;
    cudaStream_t stream = cudaq[q].stream;

    int threads = 256;
    int blocks = 1 + (m-1) / threads;

    Msgf(("%s %d sinks, %d ss_segments %d segments, %d blocks, %d threads. Base %d\n",
	  routine, m, ss_len, seg_len, blocks, threads, sink_base));
    extern int _MPMY_procnum_;
    if (strcmp(routine, "pMM_mnss_sK1") == 0)
	pMM1_mnss_sK1<<<blocks,threads,0,stream>>>
	    (sink, m, ss_index_dev, ss_seg_dev, ss_len,
	     seg_dev, seg_len, source, source_n_dev,
	     mmass, e, accp, ncut, _MPMY_procnum_);
    else Error("routine %s not found\n", routine);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
	Error("CUDA error, %d %s\n", err, cudaGetErrorString(err));
    }
}

extern "C" void
grav_qnss_CUDA(const char *routine, const int sink_base, const int m,
	       const segment *ss_seg, const float *source_base, const int source_n, 
	       int q)
{
    int align_mask = 127;	// align device pointers to 128 byte boundary
    cudaError_t err;

    cudaq[q].inuse = 1;
    err = cudaStreamCreate(&cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaStreamCreate failed, %d %s\n", err, cudaGetErrorString(err));

    char *dev = (char *)cudaq[q].dev;

    size_t ss_seg_offset = 0;
    size_t ss_seg_size = m * sizeof(segment);
    segment *ss_seg_dev = (segment *)(dev + ss_seg_offset);

    size_t source_offset = ss_seg_offset + ss_seg_size + align_mask & ~align_mask;
    size_t source_size = source_n * QSZ * sizeof(float);
    float *source_dev = (float *)(dev + source_offset);

    size_t total_size = source_offset + source_size;

    if (cudaq[q].dev_size < total_size) {
	Error("dev_size (%ld) too small, total size is %ld\n", 
	      cudaq[q].dev_size, total_size);
    }
    memcpy(cudaq[q].host, ss_seg, ss_seg_size);
    memcpy((char *)cudaq[q].host + source_offset, source_base, source_size);
    err = cudaMemcpyAsync(dev, cudaq[q].host, total_size, 
			  cudaMemcpyHostToDevice, cudaq[q].stream);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    float *sink = devpos + sink_base * 3;
    float *accp = devaccp + sink_base * 4;
    cudaStream_t stream = cudaq[q].stream;

    int threads = 64;
    int blocks = 1 + (m-1) / threads;

    Msgf(("%s %d sinks, %d sources, %d blocks, %d threads. Base %d\n",
	  routine, m, source_n, blocks, threads, sink_base));

    if (strcmp(routine, "pQ1_qnss") == 0)
	pQ1_mnss<<<blocks,threads,0,stream>>>
	    (sink, m, ss_seg_dev, source_dev, accp);
    else Error("routine %s not found\n", routine);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
	Error("CUDA error, %d %s\n", err, cudaGetErrorString(err));
    }
}
