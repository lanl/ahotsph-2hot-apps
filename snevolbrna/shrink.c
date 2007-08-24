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
	   float newt, float tpos)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p, *q;
    Stk s;
    float r2, v2, b2;

    StkInitEz(&s);

    for (*newmass = 0.0, p = btab; p < btab+*nobj; p++) {

	v2 = (p->vel[0] - b.vel[0])*(p->vel[0] - b.vel[0]) + 
	    (p->vel[1] - b.vel[1])*(p->vel[1] - b.vel[1]) + 
	    (p->vel[2] - b.vel[2])*(p->vel[2] - b.vel[2]);

	r2 = 4.0*newt*newt*b.mass*b.mass / (v2 * v2);

	b2 = (p->pos[0] - b.pos[0])*(p->pos[0] - b.pos[0]) + 
	    (p->pos[1] - b.pos[1])*(p->pos[1] - b.pos[1]) + 
	    (p->pos[2] - b.pos[2])*(p->pos[2] - b.pos[2]);

	if ( b2 >= r2 ) {  /* If distance to bndry > BH capture radius */
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
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
}
