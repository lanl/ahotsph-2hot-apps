#include "Msgs.h"
#include "physics.h"
#include "physics_sph.h"
#include "stk.h"
#include "vop.h"
#include "singlio.h"

void
ShrinkBtab
(SPHbody **SPHbtabp, body *btabp, int *nobj, float r_limit)
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
	Msgf(("Point mass gobbled m = %e; total = %e\n", p->mass,
	      btabp->accmass));
	/* Still need to update specangmom */
      }
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}
