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
AdjustBtab(SPHbody **SPHbtabp, int *nobj, bndry_t b, float *newmass, 
	   float *newr, float newt, float tpos)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p, *q;
    Stk s;
    float r2, v2, b2, minb2 = 1e30;

    StkInitEz(&s);

    for (*newmass = 0.0, p = btab; p < btab+*nobj; p++) {

	v2 = (p->vel[0] - b.vel[0])*(p->vel[0] - b.vel[0]) + 
	    (p->vel[1] - b.vel[1])*(p->vel[1] - b.vel[1]) + 
	    (p->vel[2] - b.vel[2])*(p->vel[2] - b.vel[2]);
	
	/* One option: adjust r2 based on particle velocities to
	   simulate capture-radius behavior */
	/* r2 = 4.0*newt*newt*b.mass*b.mass / (v2 * v2); */

	/* Another option: start small and move r2 out after eating
	   all particles to 10% of the radius of the next-nearest
	   particle */

	r2 = b.r*b.r;

	b2 = (p->pos[0] - b.pos[0])*(p->pos[0] - b.pos[0]) + 
	    (p->pos[1] - b.pos[1])*(p->pos[1] - b.pos[1]) + 
	    (p->pos[2] - b.pos[2])*(p->pos[2] - b.pos[2]);

	if ( b2 >= r2 ) {  /* If distance to bndry > capture radius */
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	    if (b2 < minb2) minb2 = b2;
	} else {
	    *newmass += p->mass;

	    Msgf(("t: %g: #%d: m: %g; x: %g; y: %g; z: %g; vx: %g; vy: %g; vz: %g\n", tpos, p->ident, p->mass, p->pos[0], p->pos[1], p->pos[2], p->vel[0], p->vel[1], p->vel[2]));
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));

    *newr = 0.5*sqrt(minb2);  /* Candidate new boundary radius =
				 innermost particle's
				 distance-to-boundary * 25% */
    if (*newr < b.r) *newr = b.r;  /* Never shrink boundary */
}
