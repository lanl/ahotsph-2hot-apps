#include <cm/cdpeac.h>
#include <ac.h>
#define NDIM 3
#include "vop.h"
#include "error.h"
#include "Msgs.h"

static int aux_alloc_setup = 0;
extern int last_icnt;
static void *Heap;

#define VLEN 8
#define N_DP 4

#define Vsum(x) \
  x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7]

#define DPsum(x) \
  AC_read_float(x, DP_0) + \
  AC_read_float(x, DP_1) + \
  AC_read_float(x, DP_2) + \
  AC_read_float(x, DP_3)

static void 
dp_copy(const double *pp, double *const end, aux double *space, int dp)
{
    double a0;
    space = AC_change_dp(space, dp);

    while (pp < end) {
	a0 = *pp++;
	dfdpst(a0, space++);
    }
}

static void 
do_grav_vu(aux const float *pp, aux float *const end, const float *pos0,
	   float *mass0, float *acc0, float *phi0, const float *eps2p)
{
    aux float dr2[VLEN];
    aux float r0[VLEN];
    aux float r1[VLEN];
    aux float r2[VLEN];
    aux float mass[VLEN];
    aux float mor3[VLEN];
    aux float phii[VLEN];
    aux float phi[VLEN];
    aux float a0[VLEN];
    aux float a1[VLEN];
    aux float a2[VLEN];
    aux unsigned int v[VLEN];
    aux float mtot[VLEN];
    VxdV(aux const float ppos, = pos0);
    aux const float eps2 = *eps2p;
    aux float sum;
    unsigned int i;

    phi = (float)0.0;
    mtot = (float)0.0;
    a0 = (float)0.0;
    a1 = (float)0.0;
    a2 = (float)0.0;
    for (i = 0; i < VLEN; i++) {
	v[i] = 4*i;
    }
    while (pp < end) {
	mass = (pp+0)[v];
	r0 = (pp+1)[v];
	r1 = (pp+2)[v];
	r2 = (pp+3)[v];
	pp += VLEN*4;
	VxVx(r, -= ppos);	/* 3 flops */
	dr2 = Dotx(r, r);	/* 5 flops */
	dr2 += eps2;		/* 1 flop */
	asm ("fisqtv %1,%0" : "=v" (phii) : "v" (dr2));	/* 8 flops */
	mor3 = phii * phii;	/* 4 flops */
	phii *= mass;
	mtot += mass;
	phi -= phii;
	mor3 *= phii;
	VxS(r, *= mor3);	/* 6 flops */
	VxVx(a, += r);
    }
    {

	sum = Vsum(a0);
	acc0[0] = DPsum(sum);
	sum = Vsum(a1);
	acc0[1] = DPsum(sum);
	sum = Vsum(a2);
	acc0[2] = DPsum(sum);
	sum = Vsum(phi);
	*phi0 = DPsum(sum);
	sum = Vsum(mtot);
	*mass0 = DPsum(sum);
    }
}

void
cm5_alloc_heap(int n)
{
    aux_alloc_setup = n*sizeof(int);

    /* aux_alloc_heap() does global communication */
    /* If it is called during an async phase, nodes will block */

    Heap = aux_alloc_heap(n);
}

void 
do_grav_cm5(const float *p, const float *end, const float *pos0, float *mass0,
	    float *acc0, float *phi0, const float *eps2p, int *ncut)
{
    int i, j;
    int cnt;
    int n = end-p;
    aux float *bodies;
    const int dp[N_DP] = {DP_0, DP_1, DP_2, DP_3};
    float acc[NDIM];
    float mass;
    float phi;

    if (n*sizeof(float)/N_DP > aux_alloc_setup)
      Error("Not enough aux space allocated. Have %d, need %d\n",
	    aux_alloc_setup, n*sizeof(float));

    /* n is # of floats, not bodies */
    cnt = n - last_icnt*4;
    p += last_icnt*4;

    /* Msg_do("n is %d, cnt is %d, last_icnt is %d\n", n, cnt, last_icnt); */
    
    bodies = (aux float *)Heap;
    for (j = 0; j < cnt/(4*VLEN*N_DP); j++) {
	for (i = 0; i < N_DP; i++) {
	    dp_copy((double *)(p+(i+j*N_DP)*4*VLEN),
		    (double *)(p+(i+1+j*N_DP)*4*VLEN),
		    (aux double *)(bodies+(last_icnt*4/N_DP)+j*4*VLEN), dp[i]);
	}
    }
    bodies = AC_change_dp(bodies, ALL_DPS);
    do_grav_vu(bodies, bodies+n/N_DP, pos0, &mass, acc, &phi, eps2p);
    VV(acc0, += acc);
    *mass0 = mass;
    *phi0 += phi;
    last_icnt = n/4;
}
