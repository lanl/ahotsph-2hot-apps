#include <stdio.h>
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
	p->key = q->key;
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

#include "mpmy.h"
#include "SDF.h"
#include "SDFwrite.h"
#include "singlio.h"

static Stk outstk;
static float ReduceBmax;
static float ReduceRmax2;

/* Minimum amount of information required */
typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float u;
    float h;
    unsigned int ident;		/* unique? identifier */
    int ifleos;
    float abar;
    float temp;
    float ye;
    float xp;
    float xn;
    float u2;
    float ynue;
    float ynueb;
    float ynux;
    float unue;
    float unueb;
    float unux;
    float ufreez;
} SPHminoutbody;

#define SPHMINOUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;		/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    unsigned int ident;		/* unique identifier */\n\
    int ifleos;			\n\
    float abar;			\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float ynue;			\n\
    float ynueb;		\n\
    float ynux;			\n\
    float unue;			\n\
    float unueb;		\n\
    float unux;			\n\
    float ufreez;		\n\
}"

static int 
reduce_output(tree_t *tp, hcellptr p)
{
    SPHcell *c;
    SPHbody *b;
    SPHminoutbody *q;
  
    if (p == NULL) return 0;

    if (Sub_Flags(p)) {
      c = p->ptr;
      if (c->bmax < ReduceBmax && Dot(c->pos, c->pos) < ReduceRmax2) {
	q = StkPush(&outstk, sizeof(SPHminoutbody));
	q->mass = c->mass;
	VV(q->pos, = c->pos);
	VV(q->vel, = c->vel);
	q->u = c->u;
	q->h = (c->lap + c->bmax) * 1.4;
	q->ident = c->ident;
	if (c->ifleos > 2.5) q->ifleos = 3;
	else if (c->ifleos > 1.5) q->ifleos = 2;
	else q->ifleos = 1;
	q->abar = c->abar;
	q->temp = c->temp;
	q->ye = c->ye;
	q->xp = c->xp;
	q->xn = c->xn;
	q->u2 = c->u2;
	q->ynue = c->ynue;
	q->ynueb = c->ynueb;
	q->ynux = c->ynux;
	q->unue = c->unue;
	q->unueb = c->unueb;
	q->unux = c->unux;
	q->ufreez = c->ufreez;
	return 0;
      } else {
	return 1;
      }
    } else {
      b = p->ptr;
      q = StkPush(&outstk, sizeof(SPHminoutbody));

      q->mass = b->mass;
      VV(q->pos, = b->pos);
      VV(q->vel, = b->vel);
      q->u = b->u;
      q->h = b->h;
      q->ident = b->ident;
      q->ifleos = b->ifleos;
      q->abar = b->abar;
      q->temp = b->temp;
      q->ye = b->ye;
      q->xp = b->xp;
      q->xn = b->xn;
      q->u2 = b->u2;
      q->ynue = b->ynue;
      q->ynueb = b->ynueb;
      q->ynux = b->ynux;
      q->unue = b->unue;
      q->unueb = b->unueb;
      q->unux = b->unux;
      q->ufreez = b->ufreez;
      return 0;
    }
}

void
SPHreduce(tree_t *sphtree, float bmax, float rmax,
	  char *outnamebase, int iter, 
	  float gnewt, float dt, float tpos, float tvel,
	  float rmaxnue, float rmaxnueb, float rmaxnux, 
	  float enue, float enueb, float enux, 
	  float e2nue, float e2nueb, float e2nux,
	  float ftrape, float ftrapb, float ftrapx)
{
    SPHminoutbody *output_btab;
    int output_nobj, output_gnobj;
    char outname[256];

    if (MPMY_Nproc() != 1) {
      Error("This function does not work in parallel\n");
    }
    ReduceBmax = bmax;
    ReduceRmax2 = rmax*rmax;

    StkInitEz(&outstk);

    Traverse(sphtree, Find(sphtree, KeyInt(1)), reduce_output, NULL);

    StkCrunch(&outstk);

    output_nobj = StkSz(&outstk)/sizeof(SPHminoutbody);
    output_btab = StkBase(&outstk);

    MPMY_Combine(&output_nobj, &output_gnobj, 1, MPMY_INT, MPMY_SUM);

    sprintf(outname, "%s_sphred.%04d", outnamebase, iter);

    SDFwrite(outname, output_gnobj, 
	     output_nobj, output_btab, sizeof(SPHminoutbody),
	     SPHMINOUTBODYDESC,
	     "npart", SDF_INT, output_gnobj,
	     "reduction_bmax", SDF_FLOAT, bmax,
	     "reduction_rmax", SDF_FLOAT, rmax,
	     "dt", SDF_FLOAT, dt,
	     "Gnewt", SDF_FLOAT, gnewt,
	     "iter", SDF_INT, iter,
	     "ndim", SDF_INT, NDIM,
	     "tpos", SDF_FLOAT, tpos,
	     "tvel", SDF_FLOAT, tvel,
	     "rmaxnue", SDF_FLOAT, rmaxnue,
	     "rmaxnueb", SDF_FLOAT, rmaxnueb,
	     "rmaxnux", SDF_FLOAT, rmaxnux,
	     "enue", SDF_FLOAT, enue,
	     "enueb", SDF_FLOAT, enueb,
	     "enux", SDF_FLOAT, enux,
	     "e2nue", SDF_FLOAT, e2nue,
	     "e2nueb", SDF_FLOAT, e2nueb,
	     "e2nux", SDF_FLOAT, e2nux,
	     "ftrape", SDF_FLOAT, ftrape,
	     "ftrapb", SDF_FLOAT, ftrapb,
	     "ftrapx", SDF_FLOAT, ftrapx,
	     NULL);
    Free(output_btab);
    singlPrintf("\nReduced output done, now gnob = %d\n", output_gnobj);
}
