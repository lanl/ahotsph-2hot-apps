#include <math.h>
#include "singlio.h"
#include "randoms.h"
#include "error.h"

#ifndef M_PI
#define	M_PI	3.14159265358979323846
#endif

#define c_km_s 2.99792458e5	/* speed of light in km/s */
#define H0_km_s_Kpc 0.1		/* Hubble constant (h_100) in km/s/kpc */

static setup_done;
static ran_state *rs;
static float sqrtA;
static float nu_inverse, nuprime, aprime, bprime, cprime;

void setup_EBW_cdm(float L0, float Omega0, float h, float T0, float Tquad,
		float nu, float a, float b, float c, ran_state *ranstate)
{
    float linv;
    float q;
    float two_c_H0;
    float norm;
    float h2Omega0;
    float recipsqrt_vol, recipsqrt_two;
    float kfac;

    rs = ranstate;

    /* See Efstathiou, Bond and White MNRAS 258 (1992) page 2p, eqn. (6) */
    two_c_H0 = 2. * c_km_s / (h*H0_km_s_Kpc);
    q = Tquad / T0;
    norm = q * sqrt(6./5.) * M_PI * two_c_H0*two_c_H0 * pow(Omega0, -0.77);

    /* 1/sqrt(volume) accounts for volume factor */
    recipsqrt_vol = 1.0/sqrt(L0*L0*L0);
    /* 1/sqrt(2) accounts for real and imaginary unit variance */
    recipsqrt_two = 1.0/sqrt(2.0);
    kfac = (float)2.0*M_PI/L0;

    sqrtA = recipsqrt_vol*norm*recipsqrt_two*sqrt(kfac);

    linv = 1000.0*kfac;		/* convert from 1/kpc to 1/Mpc */
    nuprime = nu;
    nu_inverse = 1.0/nu;
    h2Omega0 = h * h * Omega0;
    aprime = a * linv / h2Omega0;
    bprime = b / h2Omega0;
    bprime *= sqrt(bprime) * linv * sqrt(linv);
    cprime = c / h2Omega0;
    cprime *= cprime * linv * linv;
	
    singlPrintf("Building EBW cdm spec q=%g, h=%g, sqrtA=%g\n"
		"nu=%g, a'=%g, b'=%g, c'=%g\n" , 
		q, h, sqrtA, nu, aprime, bprime, cprime);
    setup_done = 1;
}

void
EBW_cdm(int i, int j, int k, float *real, float *imag)
{
    float a, b;
    float ksqf, kmag, ksqrt;
    float power;

    if (setup_done == 0) Error("You must call cdm_setup before pspec!\n");

    /* You must call the rng, even if you don't need values for a */
    /* given mode.  Otherwise, the parallel random numbers don't sync */
    a = normal_rand(rs);
    b = normal_rand(rs);

    if ((ksqf = i*i+j*j+k*k) == (float)0.0) {
	*real = *imag = (float)0.0;
	return;
    }
    kmag = sqrt(ksqf);
    ksqrt = sqrt(kmag);
    power = sqrtA * ksqrt * 
      pow((float)1.0+pow(aprime*kmag+bprime*kmag*ksqrt+cprime*ksqf,nuprime),
	  -nu_inverse);

    *real = a*power;
    *imag = b*power;
}
