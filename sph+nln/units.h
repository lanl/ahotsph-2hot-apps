/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* all physical constants go in here */
/* all physical constants IN CGS UNITS PLEASE!!! */
/* the conversion factors assume that all internal physics is done 
 * in  user-units, and convert FROM CGS TO USER-UNITS. They are 
 * set in the .ctl file. ADD */

extern double massCF;
extern double lenCF;
extern double timeCF;

extern double ivmassCF, ivtimeCF, ivlenCF;
extern double timeCF2, ivtimeCF2;
extern double lenCF2, ivlenCF2, ivlenCF3;
extern double ldivtCF, tdivlCF;

extern double grav_c, c_light;

#ifndef CONSTS
#define CONSTS
static const double GRAV_C =6.67428e-8;
static const double A_RAD =7.565700e-15;
static const double C_LIGHT =2.99792458e+10;
static const double KES_COEFF =6.6524586e-25;
static const double KFF_COEFF =3.680000e+22;
static const double KBF_COEFF =4.340000e+25;
static const double A_NOUGHT =5.291772108e-9;
static const double K_BOLTZ =1.3806503e-16;
static const double MH =1.67262158e-24;
static const double N_AVOG = 6.02214179e23;
#endif
