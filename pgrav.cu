/*
 * Copyright 2013 Michael S. Warren.
 * All Rights Reserved.
 */

#include <cuda.h>
#include "cudavec.h"
#include "stdio.h"

typedef unsigned long long ticks;

static __inline__ ticks getticks(void)
{
     unsigned a, d; 
     __asm__ volatile("rdtsc" : "=a" (a), "=d" (d)); 
     return ((ticks)a) | (((ticks)d) << 32); 
}

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
pH(const float *p, float *ret, const int n, const int stride, const float *f, const int source_n, int offset)
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
	VV(accp, Phi += eqe);
	VV(eqe, *= rinv2);
	VVV(accp, Ax += x, * eqe);
	VVV(accp, Ay += y, * eqe);
	VVV(accp, Az += z, * eqe);

        VV(t, *= R * rinv2);
        VV(eq0, = qx * t);
        VV(eq1, = qy * t);
        VV(eq2, = qz * t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VV(accp, Phi += eqe);
        VV(eqe, *= 3.0f * rinv2);
	VVVV(accp, Ax += x, * eqe, - eq0);
	VVVV(accp, Ay += y, * eqe, - eq1);
	VVVV(accp, Az += z, * eqe, - eq2);

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
	VV(accp, Phi -= eqe);
	VV(eqe, *= 5.0f * rinv2);
	VVVV(accp, Ax += x, * eqe, - eq0);
	VVVV(accp, Ay += y, * eqe, - eq1);
	VVVV(accp, Az += z, * eqe, - eq2);

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
	VV(accp, Phi += eqe);
	VV(eqe, *= 7.0f * rinv2);
	VVVV(accp, Ax += x, * eqe, - eq0);
	VVVV(accp, Ay += y, * eqe, - eq1);
	VVVV(accp, Az += z, * eqe, - eq2);

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
	VV(accp, Phi += eqe);
	VV(eqe, *= 9.0f * rinv2);
	VVVV(accp, Ax += x, * eqe, - eq0);
	VVVV(accp, Ay += y, * eqe, - eq1);
	VVVV(accp, Az += z, * eqe, - eq2);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (index+i < n) {
	    ret[stride*(index*VECWIDTH+i)+0] += accp[i] Ax;
	    ret[stride*(index*VECWIDTH+i)+1] += accp[i] Ay;
	    ret[stride*(index*VECWIDTH+i)+2] += accp[i] Az;
	    ret[stride*(index*VECWIDTH+i)+3] += accp[i] Phi;
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
	VV(accp, Phi += eqe);
	VV(eqe, *= rinv2);
	VVV(accp, Ax += x, * eqe);
	VVV(accp, Ay += y, * eqe);
	VVV(accp, Az += z, * eqe);

        VV(t, *= R * rinv2);
        VV(eq0, = qx * t);
        VV(eq1, = qy * t);
        VV(eq2, = qz * t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VV(accp, Phi += eqe);
        VV(eqe, *= 3.0f * rinv2);
	VVVV(accp, Ax += x, * eqe, - eq0);
	VVVV(accp, Ay += y, * eqe, - eq1);
	VVVV(accp, Az += z, * eqe, - eq2);

        VV(t, *= 3.0f * R * rinv2);
        VV(eq0,  = qxx * x);
        VV(eq1,  = qyy * y);
	VV(eq2, = -LPAREN qxx + qyy RPAREN * z);
        VV(eq0, += qxy * y);
	VV(eq1, += qxy * x);
	VV(eq2, += qxz * x);
	VV(eq0, += qxz * z);
	VV(eq1, += qyz * z);
	VV(eq2, += qyz * z);
	VV(eq0, *= t);
	VV(eq1, *= t);
	VV(eq2, *= t);
	VVV(eqe,  = eq0, * x);
	VVV(eqe, += eq1, * y);
	VVV(eqe, += eq2, * z);
	VS(eqe, *= 0.5f);
	VV(accp, Phi -= eqe);
	VV(eqe, *= 5.0f * rinv2);
	VVVV(accp, Ax += x, * eqe, - eq0);
	VVVV(accp, Ay += y, * eqe, - eq1);
	VVVV(accp, Az += z, * eqe, - eq2);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (index+i < n) {
	    ret[stride*(index*VECWIDTH+i)+0] += accp[i] Ax;
	    ret[stride*(index*VECWIDTH+i)+1] += accp[i] Ay;
	    ret[stride*(index*VECWIDTH+i)+2] += accp[i] Az;
	    ret[stride*(index*VECWIDTH+i)+3] += accp[i] Phi;
	}
    }
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
	VV(accp, Phi += eqe);
	VV(eqe, *= rinv2);
	VVV(accp, Ax += x, * eqe);
	VVV(accp, Ay += y, * eqe);
	VVV(accp, Az += z, * eqe);
    }
    for (i = 0; i < VECWIDTH; i++) {
	if (index+i < n) {
	    ret[stride*(index*VECWIDTH+i)+0] += accp[i] Ax;
	    ret[stride*(index*VECWIDTH+i)+1] += accp[i] Ay;
	    ret[stride*(index*VECWIDTH+i)+2] += accp[i] Az;
	    ret[stride*(index*VECWIDTH+i)+3] += accp[i] Phi;
	}
    }
}

#include <stdio.h>
#include <unistd.h>
#include "error.h"
#include "Msgs.h"

typedef struct cudaq_t {
    int inuse;
    cudaStream_t stream;
    int devf_size;
    float *hostf;
    float *devf;
    float *pos;
    float *accp;
} cudaq_t;

#define NQ 32
static cudaq_t cudaq[NQ];

static int
check_cudaq(void)
{
    int i;
    cudaError_t err;

    do {
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
	/* if (ninuse == NQ) printf("All slots in use\n"); */

	for (i = 0; i < NQ; i++) {
	    if (!cudaq[i].inuse) break;
	}
    } while (i == NQ);

    if (cudaq[i].devf_size == 0) {
	cudaq[i].devf_size = 16*1024*1024;
	err = cudaMallocHost((void **)&cudaq[i].hostf, cudaq[i].devf_size);
	if (err != cudaSuccess) 
	    Error("cudaMallocHost failed, %d %s\n", err, cudaGetErrorString(err));
	err = cudaMalloc((void **)&cudaq[i].devf, cudaq[i].devf_size);
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

extern "C" void
WalkInitSinkCUDA(float *btab, int stride, int64_t nobj)
{
    float *hostpos;		// pinned
    size_t pos_bytes = 3 * sizeof(float) * nobj;
    size_t accp_bytes = 4 * sizeof(float) * nobj;
    cudaError_t err;

    Btab = btab;

    /* Copy x, y, z to device through host buffer */
    err = cudaMalloc((void **)&devpos, pos_bytes);
    if (err != cudaSuccess) 
	Error("cudaMalloc failed, %d %s\n", err, cudaGetErrorString(err));

    /* Allocate ax, ay, az, phi and initialize to zero */
    err = cudaMalloc((void **)&devaccp, accp_bytes);
    if (err != cudaSuccess) 
	Error("cudaMalloc failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaMallocHost((void **)&hostpos, pos_bytes);
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

    wait_cudaq();

    err = cudaFree(devpos);
    if (err != cudaSuccess) 
	Error("cudaFree failed, %d %s\n", err, cudaGetErrorString(err));

    err = cudaMallocHost((void **)&hostaccp, accp_bytes);
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
    }

    err = cudaFreeHost(hostaccp);
    if (err != cudaSuccess) 
	Error("cudaFreeHost failed, %d %s\n", err, cudaGetErrorString(err));
}


extern "C" void
pinteractCUDA(const float *p, float *accp, const int m, const int stride, 
	      const float *f, const int source_n, const int sz)
{
    cudaError_t err;
    int blocks, threads;
    int i = check_cudaq();

    cudaq[i].inuse = 1;
    err = cudaStreamCreate(&cudaq[i].stream);
    if (err != cudaSuccess) 
	Error("cudaStreamCreate failed, %d %s\n", err, cudaGetErrorString(err));

    if (cudaq[i].devf_size < source_n*sz*sizeof(float)) {
	Error("devf_size too small\n");
    }
    memcpy(cudaq[i].hostf, f, source_n*sz*sizeof(float));
    err = cudaMemcpyAsync(cudaq[i].devf, cudaq[i].hostf, source_n*sz*sizeof(float), 
			  cudaMemcpyHostToDevice, cudaq[i].stream);
    if (err != cudaSuccess) 
	Error("cudaMemcpy failed, %d %s\n", err, cudaGetErrorString(err));

    cudaq[i].pos = devpos + 3*(p-Btab)/stride;
    cudaq[i].accp = devaccp + 4*(p-Btab)/stride;;

    if (m <= 512) {
	blocks = 1;
	threads = (m+VECWIDTH-1)/(VECWIDTH);
    } else if (m <= 1024) {
	blocks = 2;
	threads = (m+2*VECWIDTH-1)/(2*VECWIDTH);
    } else if (m <= 2048) {
	blocks = 4;
	threads = (m+4*VECWIDTH-1)/(4*VECWIDTH);
    } else if (m <= 4096) {
	blocks = 8;
	threads = (m+8*VECWIDTH-1)/(8*VECWIDTH);
    } else if (m <= 8192) {
	blocks = 16;
	threads = (m+16*VECWIDTH-1)/(16*VECWIDTH);
    } else {
	blocks = m/32;
	threads = (m+blocks*VECWIDTH-1)/(blocks*VECWIDTH);
    }

    if (sz == MSZ) {
	Msg_do("M %d sinks, %d sources, %d blocks, %d threads. Offset %ld\n",
	       m, source_n, blocks, threads, (p-Btab)/stride);
	pM<<<blocks,threads,0,cudaq[i].stream>>>(cudaq[i].pos, cudaq[i].accp, m, 4, cudaq[i].devf, source_n);
    } else if (sz == QSZ) {
	Msg_do("Q %d sinks, %d sources, %d blocks, %d thread.  Offset %ld\n",
	       m, source_n, blocks, threads, (p-Btab)/stride);
	pQ<<<blocks,threads,0,cudaq[i].stream>>>(cudaq[i].pos, cudaq[i].accp, m, 4, cudaq[i].devf, source_n);
    } else if (sz == HSZ) {
	Msg_do("H %d sinks, %d sources, %d blocks, %d threads.  Offset %ld\n",
	       m, source_n, blocks, threads, (p-Btab)/stride);
	pH<<<blocks,threads,0,cudaq[i].stream>>>(cudaq[i].pos, cudaq[i].accp, m, 4, cudaq[i].devf, source_n, (p-Btab)/stride);
    } else Error("Unknown size %d\n", sz);

    // wait_cudaq();
}
