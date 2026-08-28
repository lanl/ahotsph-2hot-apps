/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <math.h>
#include "physics_sph.h"
#include "verify.h"
#include "mpmy.h"
#include "gc.h"
#include "key.h"
#include "vop.h"

#ifndef FLT_MAX
#define FLT_MAX 1.e38
#endif

/* Call FindBbox to learn what the extent of the system is. */
void
SPHFindBbox(SPHbody *bp, int n, float *rmin, float *rmax)
{
    SPHbody *b;
    MPMY_Comm_request req;

    VS(rmax, = -FLT_MAX);
    VS(rmin, = FLT_MAX);
    for(b=bp; b<&bp[n]; b++){
	VVVV(if LPAREN rmin, > b->pos, RPAREN rmin, = b->pos);
	VVVV(if LPAREN rmax, < b->pos, RPAREN rmax, = b->pos);
	VS(if LPAREN !isfinite LPAREN b->pos, RPAREN RPAREN 
	     	  Error("Bad value for particle %d of %d\n",
		b-bp, n));
    }
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(rmin, rmin, NDIM, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(rmax, rmax, NDIM, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine_Wait(req);
}    

void SPHFixKeys(SPHbody *btab, int nobj, Key_t (*func)(const void *)){
    SPHbody *btabend= btab+nobj;

    while(btab<btabend){
	btab->key = func(btab);
	btab++;
    }
}

float SPHGetCost(const SPHbody *ptr)
{
    return (float) ptr->nterms;
}

void SPHFixId(SPHbody *btab, int nobj, int gnobj){
    int start;
    int mynobj;
    int i;

    NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &mynobj, &start);
/*     VerifyX(mynobj == nobj, Shout("mynobj=%d, nobj=%d, start=%d, gnobj=%d, nproc=%d, procnum=%d\n", mynobj, nobj, start, gnobj, MPMY_Nproc(), MPMY_Procnum())); */
    for(i=0; i<nobj; i++){
	if (btab[i].ident & (1 << 30))
	    btab[i].ident = (2*start+i) | (1 << 30);
	else
	    btab[i].ident = 2*start+i;
    }
/*     for(i=0; i<nobj; i++) { */
/* 	btab[i].ident = start + i; */
/*     } */
}

void SPHFixNterms(SPHbody *btab, int nobj){
    int i;
    for(i=0; i<nobj; i++){
	btab[i].nterms = 1;
	btab[i].grav_nterms = 1;
    }
}

/* Presumably ptr->key has been previously filled with either */
/* GetKey() or GetKeyPH */
Key_t SPHGetKeyFromStruct(const SPHbody *ptr){
    return ptr->key;
}

Key_t accbodyGetKey(const void *ptr){
    return ((accbody *)ptr)->key;
}

Key_t SPHOutIdentKey(const SPHoutbody *bp)
{
    return KeyLshift(KeyInt(bp->ident), KEYBITS/2);
}

Key_t SPHShortOutIdentKey(const SPHshortoutbody *bp)
{
    return KeyLshift(KeyInt(bp->ident), KEYBITS/2);
}
