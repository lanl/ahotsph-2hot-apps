#include "physics_sph.h"
#include "stk.h"
#include "vop.h"

void
ShrinkBtab
(body **btabp, int *nobj, float r_limit)
{
    body *btab = *btabp;
    body *p;
    Stk s;
    body *q;
    float r2;

    StkInitEz(&s);
    r2 = r_limit*r_limit;

    for (p = btab; p < btab+*nobj; p++) {
      if (Dot(p->pos, p->pos) >= r2) { /* acceptable */
	q = StkPush(&s, sizeof(body));
	*q = *p;
      }
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(body);
    btab = StkBase(&s);
    *btabp = Realloc(btab, *nobj * sizeof(body));
}
