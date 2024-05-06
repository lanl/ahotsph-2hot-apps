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
AdjustBtab(SPHbody **SPHbtabp, int *nobj, SPHbody **accbtabp, int *accnobj,
	   bndry_t b, float *newmass, float *newj, float G, float tpos, 
	   int iter, float *newr)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *atab = *accbtabp;
    SPHbody *p, *q;
    Stk s, a;
    float r2, b2, minb2 = 1e30;
    float j[NDIM], jhat[NDIM];
    float jm, jmax;

    StkInitEz(&s);
    StkInitEz(&a);

    if (b.force_r)
	r2 = b.force_r*b.force_r;
    else
	r2 = b.r*b.r;

    for(p = atab; p < atab+*accnobj; p++) {
	q = StkPush(&a, sizeof(SPHbody));
	*q = *p;
    }

    for (*newmass = 0.0, newj[0] = 0.0, newj[1] = 0.0, newj[2] = 0.0, p = btab;
	 p < btab+*nobj; p++) {

	b2 = (p->pos[0] - b.pos[0])*(p->pos[0] - b.pos[0]) + 
	    (p->pos[1] - b.pos[1])*(p->pos[1] - b.pos[1]) + 
	    (p->pos[2] - b.pos[2])*(p->pos[2] - b.pos[2]);

	if ( b2 >= ( (r2 > 1.0e-6) ? r2 : 1.0e-6 ) ) {  /* b2 less than
							   max(r2, 1.0e-6)? */
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	    if (b2 < minb2) minb2 = b2;
	} else {
	    q = StkPush(&a, sizeof(SPHbody));
	    *q = *p;

	    q->taccreted = tpos;
	    q->iteraccreted = iter;

	    *newmass += p->mass;

	    j[0] = p->mass * (p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	    j[1] = p->mass * (p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
            j[2] = p->mass * (p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);

	    jm = sqrt(j[0]*j[0] + j[1]*j[1] + j[2]*j[2]);
	    jhat[0] = j[0]/jm;
	    jhat[1] = j[1]/jm;
	    jhat[2] = j[2]/jm;

	    jmax = sqrt(G*b.mass*b.r) * p->mass;  /* jmax^2/m^2 =
						     G*M_bh*r_isco */

	    jm = ( jm < jmax ? jm : jmax );  /* jm = min(jm, jmax) */

	    j[0] = jm*jhat[0];
	    j[1] = jm*jhat[1];
	    j[2] = jm*jhat[2];

	    newj[0] += j[0];
	    newj[1] += j[1];
	    newj[2] += j[2];

	    Msgf(("t: %g: #%d: m: %g; x: %g; y: %g; z: %g; vx: %g; vy: %g; vz: %g\n", tpos, p->ident, p->mass, p->pos[0], p->pos[1], p->pos[2], p->vel[0], p->vel[1], p->vel[2]));
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));

    Free(atab);
    StkCrunch(&a);
    *accnobj = StkSz(&a)/sizeof(SPHbody);
    atab = StkBase(&a);
    *accbtabp = Realloc(atab, *accnobj * sizeof(SPHbody));

    *newr = 0.5*sqrt(minb2);
    if (*newr < b.r) *newr = b.r;  /* Never shrink boundary */
}
