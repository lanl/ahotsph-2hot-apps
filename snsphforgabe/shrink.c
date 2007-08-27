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
AdjustBtab(SPHbody **SPHbtabp, int *nobj, bndry_t b, float *newmass, float *newj, 
	   float G, float tpos)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p, *q;
    Stk s;
    float r2, b2;

    StkInitEz(&s);

    r2 = b.r*b.r;

    for (*newmass = 0.0, newj[0] = 0.0, newj[1] = 0.0, newj[2] = 0.0, p = btab; 
	 p < btab+*nobj; p++) {

	b2 = (p->pos[0] - b.pos[0])*(p->pos[0] - b.pos[0]) + 
	    (p->pos[1] - b.pos[1])*(p->pos[1] - b.pos[1]) + 
	    (p->pos[2] - b.pos[2])*(p->pos[2] - b.pos[2]);

	if ( b2 >= r2 ) {  /* If distance to bndry > capture radius */
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	} else {
	    *newmass += p->mass;
	    newj[0] += p->mass * (p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	    newj[1] += p->mass * (p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	    newj[2] += p->mass * (p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);

	    Msgf(("t: %g: #%d: m: %g; x: %g; y: %g; z: %g; vx: %g; vy: %g; vz: %g\n", tpos, p->ident, p->mass, p->pos[0], p->pos[1], p->pos[2], p->vel[0], p->vel[1], p->vel[2]));
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}
