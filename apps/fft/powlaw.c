#include <math.h>
#include "singlio.h"
#include "randoms.h"
#include "error.h"

extern float Wsq(int i, int j, int k);

#ifndef M_PI
#define	M_PI	3.14159265358979323846
#endif
#define MAXNDIM 3

#define c_km_s 2.99792458e5	/* speed of light in km/s */
#define H0_km_s_Kpc 0.1		/* Hubble constant in km/s/kpc */

static setup_done;
static ran_state *rs;
static float sqrtA;
static float powlaw_n;

void
setup_powlaw(float L0, float Omega0, float h, float normalization, float Tquad,
	     float n,  ran_state *ranstate)
{
    float norm;
    float recipsqrt_vol, recipsqrt_two;
    float kfac;

    rs = ranstate;
    powlaw_n = n;

    norm = normalization;

    /* 1/sqrt(volume) accounts for volume factor */
    recipsqrt_vol = 1.0/sqrt(L0*L0*L0);
    /* 1/sqrt(2) accounts for real and imaginary unit variance */
    recipsqrt_two = 1.0/sqrt(2.0);
    kfac = (float)2.0*M_PI/L0;

    sqrtA = recipsqrt_vol*norm*recipsqrt_two*pow(kfac, powlaw_n*0.5);
	
    singlPrintf("Building n=%f spectrum sqrtA=%g\n", powlaw_n, sqrtA);
    setup_done = 1;
}

void
powlaw(int i, int j, int k, float *real, float *imag)
{
    float a, b;
    float ksqf;
    float power;

    if (setup_done == 0) Error("You must call setup before pspec!\n");

    /* You must call the rng, even if you don't need values for a */
    /* given mode.  Otherwise, the parallel random numbers don't sync */
    a = normal_rand(rs);
    b = normal_rand(rs);

    if ((ksqf = i*i+j*j+k*k) == (float)0.0) {
	*real = *imag = (float)0.0;
	return;
    }
    power = sqrtA * pow(ksqf, powlaw_n*0.25);

    *real = a*power;
    *imag = b*power;
}
