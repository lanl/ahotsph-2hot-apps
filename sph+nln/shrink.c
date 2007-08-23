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
AdjustBtab (SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
	    int windnobj, float r_limit, float dt)
{

    /* Edit this to remove particles that go outside the boundaries,
       hardcoded here right now.  Also adjust to add particles at
       location of wind source; also hardcoded. */
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    int id;
    float r2 = r_limit*r_limit;
    double x, y, z, r;

    StkInitEz(&s);

    for (p = btab; p < btab+*nobj; p++) {
	/* Keep all particles outside of BH at origin, and
	   keep all particles inside reasonable volume of solution
	   Needs to be more easily adjustable... */

	if ( (Dot(p->pos, p->pos) >= r2)
	     && (fabs(p->pos[0]) <= 200.0) 
	     && (fabs(p->pos[1]) <= 200.0) 
	     && (fabs(p->pos[2]) <= 800.0) ) 
	    { 

	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;

	    if ( p->ident <= windnobj ) {  /* Add extra particle? */
		x = p->pos[0];
		y = p->pos[1];
		z = p->pos[2]+500.0;  /* HARDCODED WIND SOURCE */
		r = sqrt(x*x+y*y+z*z);

		if (r > r_limit + 1.0) {  /* Change this when h changes */
		    id = q->ident;
		    q->ident += windnobj;  /* Turn off particle addition for
					      recently pushed particle */
		    q = StkPush(&s, sizeof(SPHbody));

		    /* Be aware that some quantities not set here are set
		       only when exact_rho = 1 */

		    q->mass = p->mass;

		    q->pos[0] = x * r_limit / r;
		    q->pos[1] = y * r_limit / r;
		    q->pos[2] = z * r_limit / r;

		    VV(q->vel, = 0.5711/r_limit*q->pos);  /* HARDCODED */
		    q->pos[2] -= 500.0;  /* HARDCODED WIND SOURCE */
		    VVV(q->pos_last, = q->pos, - dt*q->vel);

		    q->h = 1.324;  /* Needs to be adjusted to match original
				      particles */

		    q->u = 3.230928e-04;  /* HARDCODED WIND CONDITIONS */
		    q->pr = 0.0;  /* Fixed in update_intermediate */

		    VS(q->acc, = 0.0);
		    VS(q->acc_last, = 0.0);
		    VS(q->grav_acc, = 0.0);

		    q->nterms = 1;  /* Equivalent to SPHFixNterms */

		    /* Lots of possibly-unnecessary initializations */
		    /* Without diffusion, these should all stay 0 */
		    q->dt = p->dt;
		    q->du = 0.0;
		    q->du_r = 0.0;
		    q->u_r = 0.0;
		    q->phi = 0.0;  /* Set this correctly? */

		    q->ident = id;

		    Msgf(("p->pos: %f %f %f; ident: %d\n", q->pos[0], 
			  q->pos[1], q->pos[2], q->ident));
		}
	    }
	} 
	/* Else track accreted/ejected material */
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}
