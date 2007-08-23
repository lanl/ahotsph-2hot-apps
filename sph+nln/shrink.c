#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <strings.h>
#include <math.h>
#include "fastflpt.h"
#include "Msgs.h"
#include "physics.h"
#include "physics_sph.h"
#include "stk.h"
#include "vop.h"
#include "singlio.h"

void
ShrinkBtab (SPHbody **SPHbtabp, body *btabp, int *nobj, float r_limit)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float r2;

    StkInitEz(&s);
    r2 = r_limit*r_limit;

    for (p = btab; p < btab+*nobj; p++) {
      if (Dot(p->pos, p->pos) >= r2) { /* acceptable */
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;
      } else {
	btabp->accmass += p->mass;
	btabp->l[0] += p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1];
	btabp->l[1] += p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2];
	btabp->l[2] += p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0];
	Msgf(("Point mass gobbled m = %e; total = %e\nAccreted ang momentum = (%e, %e, %e)\n", 
	      p->mass, btabp->accmass, 
	      btabp->l[0], btabp->l[1], btabp->l[2]));
      }
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}

void
ShrinkBtab2 (SPHbody **SPHbtabp, int *nobj, float r_limit)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float r2;

    StkInitEz(&s);
    r2 = r_limit*r_limit;

    for (p = btab; p < btab+*nobj; p++) {
      if (Dot(p->pos, p->pos) >= r2) { /* acceptable */
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;
      } 
/*        else { */
/*  	btabp->accmass += p->mass; */
/*  	btabp->l[0] += p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]; */
/*  	btabp->l[1] += p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]; */
/*  	btabp->l[2] += p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]; */
/*  	Msgf(("Point mass gobbled m = %e; total = %e\nAccreted ang momentum = (%e, %e, %e)\n",  */
/*  	      p->mass, btabp->accmass,  */
/*  	      btabp->l[0], btabp->l[1], btabp->l[2])); */
/*        } */
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}

void
AdjustBtab (SPHbody **SPHbtabp, int *nobj, SPHbody *w, int nwind, 
	    float dt)
{

    /* Edit this to remove particles that go outside the boundaries, 
       or inside the inner radius - both are hardcoded here right now.
       Also adjust to add particles at location of wind source; 
       also hardcoded. */

    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float rinner = 5.0;  /* Don't hardcode this here...  */
    float r2 = rinner*rinner;

    StkInitEz(&s);

    for (p = btab; p < btab+*nobj; p++) {
	/* Keep all particles outside of BH at origin, and
	   keep all particles inside reasonable volume of solution
	   Needs to be more easily adjustable... */

	if ( (Dot(p->pos, p->pos) >= r2)  /* acceptable */
	     && (fabs(p->pos[0]) <= 150.0) 
	     && (fabs(p->pos[1]) <= 150.0) 
	     && (fabs(p->pos[2]) <= 750.0) ) 
	    { 

	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	} 

    }

    /* Add particles around wind source; need to take current dt into
       account to calculate total mass to be added */
    /* Put in loop over particles in SPHwind */

    for (p = w; p < w+nwind; p++) {
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;

	/* Adjust x,y,z positions */
	/* To do: Add random rotation to wind source sphere */

	VS(q->pos, *= rinner);
	q->pos[2] -= 500.0;

	/* Set pos_last to reflect desired velocity */
	VVV(q->pos_last, = q->pos, - dt*q->vel);

	q->mass = 1.0e-4*dt/(nwind);  /* 1e-5 Msun/yr, total */
	q->rho = q->mass/(4.0/3.0*M_PI*((rinner+q->h/2.0)-
					(rinner-q->h/2.0)));
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));

/*      for (p = *SPHbtabp; p < *SPHbtabp+*nobj; p++) { */
/*  	Msgf(("p->pos: %f %f %f; id: %d\n", p->pos[0], p->pos[1], p->pos[2], */
/*  	      p->ident)); */
/*      } */
}
