/* Stolen and modified from johns ic/cdm2.c by msw on Fri Jan 27, 1995 */
/* Note that the normalization from that code is incorrect by a factor */
/* of sqrt(5/4pi) = 0.630783, since Qrms = (5*C_2/4pi)T_0  */

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
static float t2prime, t3prime, t4prime, t5prime;

void setup_Holtzman_cdm(float L0, float Omega0, float h, float T0, float Tquad,
		float t2, float t3, float t4, float t5, ran_state *ranstate)
{
    float linv;
    float q;
    float two_c_H0;
    float norm;
    float recipsqrt_vol, recipsqrt_two;
    float kfac;

    rs = ranstate;

    if (h != 0.5) 
      Error("h must be 0.5 for Holtzman, unless you get the new numbers\n");

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
    
    linv = 1000.0*kfac;		/* convert from kpc to inverse Mpc */
    t2prime = t2 * sqrt(linv);
    t3prime = t3 * linv;
    t4prime = t4 * pow(linv, 1.5);
    t5prime = t5 * linv*linv;
	
    singlPrintf("Building Holtzman cdm spec q=%g, h=%g, sqrtA=%g\n"
		"t2'=%g, t3'=%g, t4'=%g, t5'=%g\n" , 
		q, h, sqrtA, t2prime, t3prime, t4prime, t5prime);
    setup_done = 1;
}

void
Holtzman_cdm(int i, int j, int k, float *real, float *imag)
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
    power = sqrtA*ksqrt / 
      ((float)1.0 + t2prime*ksqrt + t3prime*kmag + 
       t4prime*kmag*ksqrt + t5prime*ksqf);

    *real = a*power;
    *imag = b*power;
}
