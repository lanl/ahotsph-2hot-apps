/*
 * Copyright 1992 Michael S. Warren, John K. Salmon, and
 * Gregoire S. Winckelmans. All Rights Reserved.
 */

/* Panel code */



#include <math.h>
#include "physics_panel.h"
#include "vop.h"
#include "mpmy.h"


void
Update(bodyptr btab, int n, float relax, float *residual)
{
    int i;
    bodyptr bp;
    double sumsq, sumarea;
    float Veltot[3], uzl;
    MPMY_Comm_request req;

    sumsq = 0.;
    sumarea = 0.;
    for(i=0; i<n; i++)
    {
	bp = btab+i;


/* GSW Apr16 93: So far, it uses Over-Relaxed Jacobi */

/* velocity induced by all panels with their present sources (this was
   obtained from the ""Tree Code"!) + freestream velocity.
   ... in absolute coordinate system:  */

	VVV(Veltot, = Vel(bp), + Uext(bp) );

/* ... in local coordinate system: */

	uzl = Dot (Ez(bp), Veltot);

/* Change in source

   Also, Recall that dzlsphi=2*pi... 
*/
/* GSW Apr16 93: So far, it uses Over-Relaxed Jacobi ! */

	/* sign?  There's only two possibilities... */
	bp->dsigma = relax*uzl/(6.283185308);
	bp->sigma += bp->dsigma;
	sumsq += uzl*uzl*Ip(bp);
	sumarea += Ip(bp);
	/* This is as good a place as any to re-zero the 'nterms' counter */
	/* and errsum accumulators in each body */
	bp->nterms = 1;
	bp->errsum = 0.;
	bp->errsum2 = 0.;
    }
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&sumsq, &sumsq, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&sumarea, &sumarea, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    *residual = sqrt(sumsq/sumarea);
}

