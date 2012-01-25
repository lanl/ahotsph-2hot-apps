/* Add and subtract SPH particles from Nbody*/
/*-SD 
Modification History, starting 08/23/2005
-----------------------------------------
  V 2.0: Original version
  V 2.1: Added function SPHPlusSPH
  V 2.2: Added function NbodyPlusNbody
  
*/


#include "stk.h"
#include "physics.h"
#include "physics_sph.h"
#include "vop.h"
#include "bigmalloc.h"
#include "error.h"
#include "singlio.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "mpmy_abnormal.h"

#define SPH_FLAG (1<<31)

void
GravPlusSPH(void **btabp, int *nobj, SPHbody *SPHbtab, int SPHnobj)
{
    int i;
    body *btab, *p;
    SPHbody *q;
    int grav_nobj = *nobj;

    *nobj += SPHnobj;
    btab = Realloc(*btabp, *nobj * sizeof(body));
    for (i = 0; i < SPHnobj; i++) {
	q = SPHbtab+i;
	p = btab + grav_nobj + i;
	p->mass = q->mass;
	p->h = q->h;
	VV(p->pos, = q->pos);
	if (q->ident & SPH_FLAG) Error("SPH flag already set\n");
	p->nterms = q->grav_nterms;
	p->ident = q->ident | SPH_FLAG;
    }
    *btabp = btab;
}



void
GravMinusSPH(void **btabp, int *nobj, accbody **atab, int *anobj)
{
    body *btab = *btabp;
    body *p, *next;
    Stk s;
    accbody *q;

    StkInitEz(&s);

    /* Shrink btab, taking out SPH particles and copying them to atab */
    for (p = next = btab; p < btab+*nobj; p++) {
	if (p->ident & SPH_FLAG) {
	    q = StkPush(&s, sizeof(accbody));
	    VV(q->grav_acc, = p->acc);
	    q->phi = p->phi;
	    q->grav_nterms = p->nterms;
	    q->ident = p->ident & ~SPH_FLAG;
	    q->key = p->key;
	} else *next++ = *p;
    }
    StkCrunch(&s);
    *anobj = StkSz(&s)/sizeof(accbody);
    *atab = StkBase(&s);
    *nobj = next-btab;
    *btabp = Realloc(btab, *nobj * sizeof(body));
}


void
SPHPlusSPH(void **btabp, int *nobj, SPHbody *SPHbtab, int SPHnobj)
{
    int i,j;
    SPHbody *btab;
    int sumnobj, grav_nobj = *nobj;
    MPMY_Comm_request req;
        /* Note: I didn't bother to copy every single tag of the structure */

    sumnobj=grav_nobj;
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&sumnobj, &sumnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    *nobj += SPHnobj;
    btab = Realloc(*btabp, *nobj * sizeof(SPHbody));
    for (i = 0; i < SPHnobj; i++) {
	j=i+grav_nobj;
	btab[j].mass = SPHbtab[i].mass;
	VV(btab[j].pos, = SPHbtab[i].pos);
	btab[j].grav_nterms = SPHbtab[i].grav_nterms;
	btab[j].ident = SPHbtab[i].ident+sumnobj;
	btab[j].type = SPHbtab[i].type;
	btab[j].h=SPHbtab[i].h;
	VV(btab[j].vel,= SPHbtab[i].vel);
	VV(btab[j].lvel,= SPHbtab[i].lvel);
	btab[j].u=SPHbtab[i].u;
	btab[j].pr=SPHbtab[i].pr;
	btab[j].rho=SPHbtab[i].rho;
	btab[j].rho_est=SPHbtab[i].rho_est;
	btab[j].vsound=SPHbtab[i].vsound;
	btab[j].udot=SPHbtab[i].udot;
	btab[j].udot_last=SPHbtab[i].udot_last;
/* 	btab[j].dt=SPHbtab[i].dt; */
/* 	btab[j].tacc=SPHbtab[i].tacc; */
	btab[j].min_nbr_dt=SPHbtab[i].min_nbr_dt;
/* 	VV(btab[j].acc,= SPHbtab[i].acc); */
/* 	VV(btab[j].grav_acc,= SPHbtab[i].grav_acc); */
/* 	VV(btab[j].acc_last,= SPHbtab[i].acc_last); */
    }
    *btabp = btab;
}



void
NbodyPlusNbody(void **btabp, int *nobj, body *btab2, int nobj2)
{
    int i,j;
    body *btab;
    int sumnobj, grav_nobj = *nobj;
    MPMY_Comm_request req;

    sumnobj=grav_nobj;
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&sumnobj, &sumnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    *nobj += nobj2;
    btab = Realloc(*btabp, *nobj * sizeof(body));
    for (i = 0; i < nobj2; i++) {
	j=i+grav_nobj;
	btab[j].mass = btab2[i].mass;
	btab[j].h = btab2[j].h;
	VV(btab[j].pos, = btab2[i].pos);
	btab[j].nterms = btab2[i].nterms;
	btab[j].ident = btab2[i].ident+sumnobj;
	btab[j].type = btab2[i].type;
	VV(btab[j].vel,= btab2[i].vel);
	VV(btab[j].acc,= btab2[i].acc);
	btab[j].phi=btab2[i].phi;
	btab[j].key=btab2[i].key;
	VV(btab[j].pos_last,= btab2[i].pos_last);
    }
    *btabp = btab;
}



void
NbodyPlusNbodyold(void **btabp, int *nobj, body *btab2, int nobj2)
{
    int i;
    body *btab, *p;
    body *q;
    int grav_nobj = *nobj;
    /* Note: I didn't bother to copy every single tag of the structure */


    *nobj += nobj2;
    btab = Realloc(*btabp, *nobj * sizeof(body));
    for (i = 0; i < nobj2; i++) {
	q = btab2+i;
	p = btab + grav_nobj + i;
	p->mass = q->mass;
	VV(p->pos, = q->pos);
	p->nterms = q->nterms;
	p->ident = q->ident+grav_nobj;
	VV(p->vel,= q->vel);
	VV(p->acc,= q->acc);
        p->phi=q->phi;
        p->key=q->key;
	VV(p->pos_last,= q->pos_last);
	
    }
    *btabp = btab;
}
