#include "stk.h"
#include "physics.h"
#include "physics_sph.h"
#include "vop.h"
#include "bigmalloc.h"
#include "error.h"

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
#ifdef SPH_GRAV
	p->h = q->h;
#endif
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
