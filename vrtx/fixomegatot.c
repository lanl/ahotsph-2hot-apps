/* fixomegatot.c */

/*
Copyright 1992, 1993, 1994, 1995. All Rights Reserved.
Michael S. Warren, John K. Salmon, Gregoire S. Winckelmans
*/

/* 
October 1995: Fix the particle strengths so as to make sure that the sum of all
particle strengths is zero.
*/

#include "physics_vrtx.h"
#include "vop.h"

extern double omega_tot[3];

void FixOmegaTot(bodyptr btab, int nobj, int gnobj)
{
    int i;
    bodyptr bp;
    double term[3];


    term[0]=omega_tot[0]/gnobj;
    term[1]=omega_tot[1]/gnobj;
    term[2]=omega_tot[2]/gnobj;

    for(i=0; i<nobj; i++)
    {
	bp = btab+i;

        VV(Strength(bp), -= term);
    }

}

