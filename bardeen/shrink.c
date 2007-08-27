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
AdjustBtab (SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
	    int windnobj, int windpartpershell, float dt, 
	    int iter, float tpos, int *added_particles)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    unsigned int id;
    float r_wind = 5.0; /* Decoupled from inner boundary size */
    float wr, wr2;
    double wpos[NDIM];

    float d = r_wind * 
	sqrt( 4.0-1.0 / (pow( sin( M_PI*(windpartpershell)/
				   (6.0*((windpartpershell)-2)) ), 2.0 )) );

    StkInitEz(&s);

    for (p = btab; p < btab+*nobj; p++) {
	/* Keep all particles outside of BH at origin, and
	   keep all particles inside reasonable volume of solution */

	if ( (fabs(p->pos[0]) <= 3400.0) 
	     && (fabs(p->pos[1]) <= 3400.0) 
	     && (fabs(p->pos[2]) <= 3400.0) ) { 

	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;

	    if ( p->windid < windnobj ) {  /* Particle on inner shell? */
		VVV(wpos, = p->pos, - windbtab[p->windid].pos);
		wr = sqrt(Dot(wpos, wpos));

		if (wr > r_wind + 0.8*d){ /* Particle far from source? */
		    *added_particles = 1;  /* Indicate particle addition */

		    id = q->windid;
		    q->windid += windnobj;  /* Turn off addition for
					       recently pushed particle */
		    q = StkPush(&s, sizeof(SPHbody));

		    /* Be aware that some quantities not set here are set
		       only when exact_rho = 1 */

		    q->mass = p->mass;

		    VVS(q->pos, = wpos, * r_wind / wr);
		    VV(q->vel, = windbtab[id].vwind/r_wind*q->pos);
		    VV(q->pos, += windbtab[id].pos);

		    VVV(q->pos_last, = q->pos, - dt*q->vel);

		    q->h = 1.8*d;  /* Match calculation in writewind.c */

		    q->u = (p->u + windbtab[id].uwind) / 2.0;
		    q->udot = 0.0;
		    q->udot_last = 0.0;  /* Just in case */
		    q->pr = 0.0;  /* Fixed in update_intermediate */

		    VS(q->acc, = 0.0);
		    VS(q->acc_last, = 0.0);
		    VS(q->grav_acc, = 0.0);

		    q->nterms = 1;  /* Equivalent to SPHFixNterms */
		    q->grav_nterms = 1;

		    q->tacc = -1e30;

		    /* Lots of possibly-unnecessary initializations */
		    /* Without diffusion, these should all stay 0 */
		    q->dt = q->dt_next = dt;  /* CORRECT?? */
		    q->min_nbr_dt = 1e30;  /* Just testing */
		    q->du = 0.0;
		    q->du_r = 0.0;
		    q->u_r = 0.0;
		    q->phi = 0.0;  /* Set this correctly? */

		    q->windid = id;
		    q->ident = 0;  /* Fix in subsequent call to SPHFixId */

/* 		    Msgf(("p->pos: %f %f %f; windid: %d; u: %e\n",  */
/* 			  q->pos[0], q->pos[1], q->pos[2], q->windid,  */
/* 			  q->u)); */
		}
	    }
	}
	/* Else track accreted/ejected material; do this right sometime */
/* 	else { */
/* 	    Msgf(("%d: %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %u\n", iter, tpos, p->pos[0], p->pos[1], p->pos[2], p->vel[0], p->vel[1], p->vel[2], p->mass, p->rho, p->u, p->h, p->windid)); */
/* 	} */

    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}


void
AdjustBtab2 (SPHbody **SPHbtabp, int *nobj, float r_inner)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p, *q;
    Stk s;
    float r2 = r_inner*r_inner;

    StkInitEz(&s);

    for (p = btab; p < btab+*nobj; p++) {

	if (Dot(p->pos, p->pos) > r2) {
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}
