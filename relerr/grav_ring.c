#define NO_MSGS
/* Msgs here really slow things down, even if they aren't activated */
/*
 * Copyright 1994 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#include "physics_n.h"
#include "vop.h"
#include "tensop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"

Counter_t GravCnt;

static float Eps2;

void
set_eps(float eps)
{
    Eps2 = eps*eps;
}

void
set_body(void *o, void *p)
{
    memcpy(o, p, TBODYSZ);
}

void 
do_grav2(void *p0, void *list, int bsize, int n)
{
    body *p = p0;
    int ncut;
    float total_mass;
    char *clist = (char *)list;
#ifdef __PARAGON__
    int nn;
#endif

    AddCounter(&GravCnt, n);

#ifdef __PARAGON__
    nn = n - (n % 3);
    /* Use the interface to the fast assembly code */
    if ((int)list & 07 || bsize & 07)
      Error("btab not aligned for asm code\n");
    do_grav_fast(list, (float *)(clist + nn * bsize), p->pos, 
		 &total_mass, p->acc, &p->phi, &Eps2, &ncut);
    if (n % 3)
      do_grav((char *)list+nn*bsize, (float *)(clist + n * bsize), p->pos, 
	      &total_mass, p->acc, &p->phi, &Eps2, &ncut);
#else
    do_grav(list, (float*)(clist + n * bsize), 
	    p->pos, &total_mass, p->acc,
	    &p->phi, &Eps2, &ncut);
#endif
}

