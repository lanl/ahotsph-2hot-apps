/* relaxomega.c */

/*
Copyright 1992, 1993, 1994, 1995. All Rights Reserved.
Michael S. Warren, John K. Salmon, Gregoire S. Winckelmans
*/

/* 
May 1995: particle strength update in Winckelmans' relaxation  scheme
(see thesis and JCP paper).
*/


#include "physics_vrtx.h"
#include "vop.h"

void RelaxOmega(bodyptr btab, int n, float relaxw)
{
    int i;
    bodyptr bp;
    float omega[3];

    for(i=0; i<n; i++)
    {
	bp = btab+i;

        /* divergence free vorticity field that the new particle strengths
           will ``best" approximate after all iterations of the relaxation
           scheme.  */

        omega[0] = Gradvel(bp)[2][1]-Gradvel(bp)[1][2];
        omega[1] = Gradvel(bp)[0][2]-Gradvel(bp)[2][0];
        omega[2] = Gradvel(bp)[1][0]-Gradvel(bp)[0][1];

        /* incremental change of particle strength for this iteration */

        Strength(bp)[0] += relaxw*Vol(bp)*(omega[0]-Omegat(bp)[0]);
        Strength(bp)[1] += relaxw*Vol(bp)*(omega[1]-Omegat(bp)[1]);
        Strength(bp)[2] += relaxw*Vol(bp)*(omega[2]-Omegat(bp)[2]);
    }

}

