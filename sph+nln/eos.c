/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <stdio.h>
#include <math.h>
#include <Msgs.h>
#include "singlio.h"
#include "units.h"

#define MAX_ITER 200

/* Need to change printf's to Msgfs or whatever */

float newtraph(double xl, double xr, double prec, double (*f)(double x), 
	     double (*df)(double x)) {
    double xguess = 2.0*xl*xr/(xl+xr);
    int i;
   
    if ((*f)(xl)*(*f)(xr) > 0.0) {
	singlPrintf("Bisect: %1.1e and %1.1e do not bracket a root\n", xl, xr);
        xguess = -99.0;
        return (float)xguess;
    }

    for(i = 0; i <= MAX_ITER; i++) {
	xguess = xguess - (*f)(xguess)/(*df)(xguess); /* Try Newton-Raphson */

	if (xguess < xl || xguess > xr)
	    xguess = (xl + xr) / 2.0; /* Fall back to bisection */
    
	if ( fabs((*f)(xguess)) < prec )
	    break;
	else if ((*f)(xl)*(*f)(xguess) < 0.0)
	    xr = xguess;
	else
	    xl = xguess;
    }

    if (i == MAX_ITER) {
	singlPrintf("Bisect: max iterations exceeded\n");
        xguess = -99.0;
    }

    return (float)xguess;
}


/* changed: these now need to be in cgs */
extern double eos_n;
extern double eos_u;

#include "physics_sph.h"

double uvst(double t) {
    return 1.5*eos_n*(K_BOLTZ)*t + t*t*t*t*(A_RAD) - eos_u;
}


double duvst(double t) {
    return 1.5*eos_n*(K_BOLTZ) + 4.0*t*t*t*(A_RAD);
}
