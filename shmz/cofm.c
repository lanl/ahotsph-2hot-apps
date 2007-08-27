#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tree.h"
#include "physics.h"
#include "vop.h"
#include "bigmalloc.h"
#include "Msgs.h"
#include "key.h"
#include "fmm.h"

static double K;
static int Order[KEYBITS];
static int Lowest_level;
static int Nobj;
static int Nd = 1;

#define N_u0 72
#define N_v0 144
static double u0[N_u0], v0[N_v0];
static double w[32];
static int old_nu = -1;

#define MUMX 128
#define MVMX 256

static complex do2o[MVMX*MUMX];
static complex do2i[MVMX*MUMX];
static complex di2i[MVMX*MUMX];

static int stenciloo[2][2][2], isoo, isstoo[3][4];
static int *pdo2oo, *pdo2ii;
static complex *stdo2oo, *stdo2ii;

static void st_out_out(void)
{
    int i;
    int order[KEYBITS];
    double rotate[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double dimen = 0.13;
    int level = 3;
    
    for (i = 0; i < 5; i++)
	order[i] = Order[4-i];	/* reverse to Fortran concept */
    
    pdo2oo = Malloc(level*4 * sizeof(int));
    pdo2ii = Malloc(level*4 * sizeof(int));
    stdo2oo = Malloc(14720 * sizeof(complex));
    stdo2ii = Malloc(30496 * sizeof(complex));
    Fortran2(st_out_out)(&level, &K, order, rotate, &dimen, stdo2oo, pdo2oo,
		 stdo2ii, pdo2ii, stenciloo, &isoo, isstoo);
}

void SetupCofm(int order, double lambda)
{
    int i1;
    K = 2.0*M_PI/lambda;
	
    /* Tree Levels are opposite Fortran code.  Level 0 is at the top */
    Lowest_level = 4;
    /* plan 6 - 12 - 20 - 36 - 64 */
    Order[5] = 6;
    Order[4] = 6;
    Order[3] = 12;
    Order[2] = 20;
    Order[1] = 36;
    Order[0] = 64;
    /* This generates some shared data which is used deep inside cfix2y */
    Fortran(genabm)();
    st_out_out();
    i1 = 2*Order[4]-1;
    Fortran(shfqwt)(&i1, w);
}

/* should be in key.c */
static void 
CellICorner(Key_t key, int *icorner)
{
    unsigned int iscale = 1;
    int i;

    VS(icorner, = 0);
    while(KeyGT(key, KeyInt(1))) {
	for(i=0; i<NDIM; i++){
	    if( KeyAndInt(key, (1<<i)) )
	      icorner[i] |= iscale;
	}
	key = KeyRshift(key, NDIM);
	iscale <<= 1;
    }
}

static void cfix2y (int mkpure, int ier, int n3, 
		    int ndegl, int nul, int nvl, complex *ffsg, 
		    int ndeg, int nu, int nv,  complex *ffsf)
{
    complex *tmp;
    
    /* cfix2y changes ffsg, which is not really a good idea */
    tmp = Malloc(nul * nvl * sizeof(complex));
    memcpy(tmp, ffsg, nul * nvl * sizeof(complex));
    ier = 0;
    Fortran(cfix2y)(&mkpure, &ier, &n3, &ndegl, &nul, &nvl, tmp,
	    &ndeg, &nu, &nv, ffsf);
    if (ier != 0) Error("cfix2y reports an error\n");
    Free(tmp);
}



static void ex_out_out(Key_t key, Key_t dkey, int *pdo2xx, 
		       complex *stdo2xx, complex *dx2x)
{
    int ll;
    int level = 3;
    int cubei[3], cubej[3];

    /* Convert Key coords to Fortran code equivalent */

    ll = 1+level-TreeLevel(key, NDIM);
    CellICorner(KeyLshift(dkey, (ll-1)*NDIM), cubei);
    CellICorner(KeyLshift(key, ll*NDIM), cubej);
    VS(cubei, += 1 << (ll-1));
    VS(cubej, += 1 << ll);		/* I'm not sure why this is right */

    Fortran2(ex_out_out)(cubei, cubej, &level, &ll, stenciloo, &isoo, isstoo, 
		 pdo2xx, stdo2xx, dx2x);
}


/* This is where the outer to outer translations happen */
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[])
{
    int i;
    float center[NDIM], cellsz;
    double r[NDIM], strength;
    cofmdata *cmp, *dp;
    body *bp;
    complex *ffsf;
    int nu, nv;
    int nul, nvl, ndeg;
    int level;
    Key_t daughter;
    
    level = TreeLevel(hptr->key, NDIM);
    if (level > Lowest_level) {
	/* SeriousWarning("TreeLevel deeper than desired\n"); */
	nu = Order[Lowest_level];
	nv = 2*nu;
	CellCorner(KeyRshift(hptr->key, (level-Lowest_level)*NDIM), 
		   center, &cellsz);
    } else {
	nu = Order[level];
	nv = 2*nu;
	CellCorner(hptr->key, center, &cellsz);
    }
    cellsz *= (float)0.5;
    VS(center, += cellsz);
    cmp = hptr->ptr;
    cmp->strength = 0.0;
    VV(cmp->pos, = center);
    cmp->sz = cellsz;
    cmp->nu = cmp->nv = cmp->ndeg = 0;
    cmp->ffsf = NULL;
    cmp->ndaughters = 0;

    if (level <= 1) {
	for(i=0; i<(1<<NDIM); i++){
	    if (daughters[i] == NULL)
		continue;
	    if (Sub_Flags(daughters[i]) == 0) {
		bp = daughters[i]->ptr;
		cmp->strength += bp->strength;
		cmp->ndaughters++;
	    } else {
		dp = daughters[i]->ptr;
		cmp->strength += dp->strength;
		cmp->ndaughters += dp->ndaughters;
	    }
	}
	return;
    }

    if (nu != old_nu) {
	old_nu = nu;
	Fortran(gnthph)(&nu, &nv, u0, v0);
    }

    ffsf = Malloc(nu * nv * sizeof(complex));
    cmp->nu = nu;
    cmp->nv = nv;
    cmp->ndeg = nu-1;
    cmp->ffsf = Malloc(nu * nv * sizeof(complex));
    Fortran(sfzero)(&nu, &nv, cmp->ffsf, &Nd);

    for(i=0; i<(1<<NDIM); i++){
	if (daughters[i] == NULL)
	  continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    VVV(r, = center, - bp->pos);	/* seems backwards */
	    cmp->strength += bp->strength;
	    cmp->ndaughters++;
	    Fortran(gendto)(&nu, &nv, do2o, u0, v0, r, &K);
	    strength = bp->strength;
	    Fortran(genffsf)(&strength, &nu, &nv, ffsf, do2o, &Nd);
	    Fortran(sfadd)(&nu, &nv, ffsf, cmp->ffsf, cmp->ffsf, &Nd);
	} else {
	    dp = daughters[i]->ptr;
	    cmp->strength += dp->strength;
	    cmp->ndaughters += dp->ndaughters;
	    if (dp->nu == nu) {
		/* This happens when there are more levels in the tree */
		/* than we want. */
		if (   dp->pos[0] != center[0] 
		    || dp->pos[1] != center[1] 
		    || dp->pos[2] != center[2]) Error("Offset mismatch\n");
		Fortran(sfadd)(&nu, &nv, dp->ffsf, cmp->ffsf, cmp->ffsf, &Nd);
	    } else {
		if (level >= 6) Error("level too large\n");
		nul = Order[level+1];
		nvl = 2*nul;
		ndeg = nul-1;

		if (dp->nu != nul || dp->ndeg != ndeg)
		    Error("nu or ndeg mismatch\n");
		cfix2y(0, 0, 1, ndeg, nul, nvl, dp->ffsf, ndeg, nu, nv, ffsf);
		daughter = KeyOrInt(KeyLshift(hptr->key, NDIM), i);
		ex_out_out(hptr->key, daughter, pdo2oo, stdo2oo, do2o);
		Fortran(sfmult)(&nu, &nv, ffsf, do2o, ffsf, &Nd);
		Fortran(sfadd)(&nu, &nv, ffsf, cmp->ffsf, cmp->ffsf, &Nd);
	    }
	} 
    }
    Free(ffsf);
}

/* Turn the ptr from a cofmdata to a cell. */
void CellFromCofm(cell *cp, cofmdata *cmp)
{
    cp->strength = cmp->strength;
    VV(cp->pos, = cmp->pos);
    cp->sz = cmp->sz;
    cp->daughters = cmp->ndaughters;
    cp->nu = cmp->nu;
    cp->nv = cmp->nv;
    cp->ndeg = cmp->ndeg;
    cp->ffsf = cmp->ffsf;
}

void SetTol(int gnobj)
{
    Nobj = gnobj;
}

void InheritSink(const Sink *from, Sink *to, hcell *pp)
{
    int level;

    if (to == NULL) {
	int one = 1;
	int nu, nv, ndeg;
	complex *ffsf;
	complex phi_f;
	/* Copy data from terminal sink to body */
	body *bp = pp->ptr;
	/* Make sure these are initialized to zero externally */
	bp->phi_r += from->phi.r;
	bp->phi_i += from->phi.i;
	bp->nterms += from->nterms;
	nu = 2*Order[4];
	nv = 2 * nu;
	ndeg = nu -1;
	ffsf = Malloc(nu * nv * sizeof(complex));
	Msgf(("%f %f %d\n", from->ffsf[0].r, from->ffsf[0].i, from->nu));
	cfix2y(0, 0, 1, from->nu-1, from->nu, from->nv, from->ffsf, 
	       ndeg, nu, nv, ffsf);
	Msgf(("%f %f %d\n",  ffsf[0].r, ffsf[0].i, nu));
	phi_f.r = phi_f.i = 0.0;
	Fortran(eval)(&nu, &nv, w, ffsf, &phi_f, &one);
	Msgf(("far %f %f\n", phi_f.r, phi_f.i));
	Msgf(("near %f %f\n", bp->phi_r, bp->phi_i));
	bp->phi_r += phi_f.r;
	bp->phi_i += phi_f.i;
	Free(ffsf);
	Free(from->ffsf);
	if (from->interactions != Nobj)
	    Error("Ninteract is %d, should be %d\n", 
		  from->interactions, Nobj);
	Msgf(("done %f %f %f, %f %f\n",
	       bp->pos[0], bp->pos[1], bp->pos[2], bp->phi_r, bp->phi_i));
	/* printf("done %f %f %f, %f %f\n",
	       bp->pos[0], bp->pos[1], bp->pos[2], bp->phi_r, bp->phi_i); */
	return;
    }

    to->key = pp->key;
    level = TreeLevel(to->key, NDIM);
    if (Sub_Flags(pp)) {
	cell *cp = pp->ptr;
	VV(to->pos, = cp->pos);
	to->sz = cp->sz;
	to->isbody = 0;
	to->daughters = cp->daughters;
	to->nu = Order[level];
	to->nv = 2*to->nu;
	to->ndeg = to->nu-1;
	if (level >= 2) {
	    /* Allocate space for inner signature function */
	    to->ffsf = Malloc(to->nu * to->nv * sizeof(complex));
	    Fortran(sfzero)(&to->nu, &to->nv, to->ffsf, &Nd);
	} else {
	    to->ffsf = NULL;
	    to->nu = to->nv = to->ndeg = 0;
	}
    } else {
	body *bp = pp->ptr;
	VV(to->pos, = bp->pos);
	to->phi.r = to->phi.i = 0.0;
	to->sz = 0.0;
	to->isbody = 1;
	to->daughters = 1.0;
	to->nu = Order[(level < 4) ? level : 4];
	to->nv = 2*to->nu;
	to->ndeg = to->nu-1;
	to->ffsf = Malloc(to->nu * to->nv * sizeof(complex));
	Fortran(sfzero)(&to->nu, &to->nv, to->ffsf, &Nd);
    }

    if (from && from->ffsf) {
	/* This is where the inner to inner translation happens */
	int level;
	int nu, nv, nul, nvl, nun, nvn;
	int ndeg, ndegi, ndegt;
	complex *ffsf, *ffsg;
	double r[NDIM];

	to->interactions = from->interactions;
	to->nterms = from->nterms;

	level = TreeLevel(from->key, NDIM);

	if (level < 2 && !to->isbody) return;

	Msgf(("Inherit %d %s %d\n", level, PrintKey(to->key), to->isbody));

	if (level >= 3) {
	    if (!to->isbody) {
		Msgf(("Deferring level %d\n", level));
		to->ffsf = from->ffsf;
		to->nu = from->nu;
		to->nv = from->nv;
		to->ndeg = from->ndeg;
		VV(to->pos, = from->pos);
		return;
	    } else level = 3;
	}

	nu = Order[level+1];
	nv = 2 * nu;
	nul = Order[level];
	nvl = 2 * nul;
	nun = nul + (nul-nu);
	nvn = 2 * nun;
	ndeg = nul-1;
	ndegi = ndeg + (nul-nu);
	ndegt = nu - 1;

	ffsf = Malloc(nun * nvn * sizeof(complex));
	ffsg = Malloc(nun * nvn * sizeof(complex));
	if (from->nu != nul) Error("nu mismatch\n");
	if (from->ndeg != ndeg) Error("ndeg mismatch\n");
	Msgf(("from %f %f %f, %f %f %d\n", 
	      from->pos[0], from->pos[1], from->pos[2],
	      from->ffsf[0].r, from->ffsf[0].i, from->nu));
	cfix2y(0, 0, 1, ndeg, nul, nvl, from->ffsf, ndeg, nun, nvn, ffsg);
	if (to->isbody) {
	    if (nun != old_nu) {
		old_nu = nun;
		Fortran(gnthph)(&nun, &nvn, u0, v0);
	    }
	    VVV(r, = to->pos, - from->pos);
	    Fortran(gendto)(&nun, &nvn, di2i, u0, v0, r, &K);
	} else {
	    ex_out_out(from->key, to->key, pdo2ii, stdo2ii, di2i);
	}
	Msgf(("%f %f\n", di2i[0].r, di2i[0].i));
	Fortran(sfmult)(&nun, &nvn, ffsg, di2i, ffsf, &Nd);
	Msgf(("%f %f\n", ffsf[0].r, ffsf[0].i));
	if (to->nu != nu) Error("nu mismatch\n");
	if (to->ndeg != ndegt) Error("ndeg mismatch\n");
	cfix2y(0, 0, 1, ndegi, nun, nvn, ffsf, ndegt, nu, nv, to->ffsf);
	Msgf(("  to %f %f %f, %f %f %d\n", 
	      to->pos[0], to->pos[1], to->pos[2],
	      to->ffsf[0].r, to->ffsf[0].i, to->nu));
	Free(ffsg);
	Free(ffsf);
    } else {
	to->interactions = 0;
	to->nterms = 0;
    }
}

/* This is where the Outer to Inner translations happen */
void OutToIn(Sink *sink, const hcell **source_vec, int *result, int n)
{
    double r[NDIM];
    int i;
    int level, slevel;
    int nu, nv, nul, nvl, ndeg, ndegi;
    complex *ffsf, *ffsg, *ffsh;
    int l1, l2, l3;

    level = TreeLevel(sink->key, NDIM);

    if (level < 2) {
	/* Msgf(("splitting sink for %d sources\n", n)); */
	for (i = 0; i < n; i++) result[i] = MAC_SPLIT_SINK;
	return;
    }

    Msgf(("%d sink %s %d\n", level, PrintKey(sink->key), n));

    if (level > 3 && !sink->isbody) {
	/* Msgf(("splitting sink for %d sources\n", n)); */
	for (i = 0; i < n; i++) result[i] = MAC_SPLIT_SINK;
	return;
    }
    
    if (sink->isbody && level > 4) level = 4;

    nu = Order[level-1];
    nv = 2*nu;
    nul = Order[level];
    nvl = 2*nul;
    ndeg = nul-1;
    ndegi = nu-1;
    l1 = Order[level-1];
    l2 = 2 * l1;
    l3 = Order[level-1]-Order[level];
    if (sink->nu != nul) Error("nu mismatch\n");

    if (old_nu != l1) {
	old_nu = l1;
	Fortran(gnthph)(&l1, &l2, u0, v0);
    }

    ffsf = Malloc(nu * nv * sizeof(complex));
    ffsg = Malloc(nu * nv * sizeof(complex));
    ffsh = Malloc(nul * nvl * sizeof(complex));

    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];
	
	if (Sub_Flags(source_vec[i])) {
	    const cell *cp = source_vec[i]->ptr;
	    float sz;
	    /* Msgf(("%s\n", PrintKey(source->key))); */
	    slevel = TreeLevel(source->key, NDIM);
	    if (slevel < level && !sink->isbody) {
		result[i] = MAC_SPLIT_SRC;
		/* Msgf(("splitting source %d\n", cp->daughters)); */
		continue;
	    }
	    if (!sink->isbody) sz = sink->sz;
	    else sz = cp->sz;
	    VVV(r, = sink->pos, - cp->pos);
	    if (Dot(r, r) < 15*sz*sz) {
		if (sink->isbody)
		    result[i] = MAC_SPLIT_SRC;
		else
		    result[i] = MAC_SPLIT_SINK;
		/* Msgf(("%d reject, %f %d\n", i, sqrt(Dot(r, r))/sz, 
		      cp->daughters)); */
		continue;
	    }
	    /* Msgf(("%d accept, %f %d\n", i, sqrt(Dot(r, r))/sz, 
		  cp->daughters)); */
	    /* Msgf(("%s\n", PrintCellContents(cp))); */
	    sink->interactions += cp->daughters;
	    result[i] = MAC_ACCEPT;
	    Msgf(("%s %d %d %f %f\n", PrintKey(source->key), 
		  cp->nu, cp->nv, cp->ffsf[0].r, cp->ffsf[0].i));
	    cfix2y(0, 0, 1, cp->ndeg, cp->nu, cp->nv, cp->ffsf, 
		   ndeg, nu, nv, ffsf);
	    if (l1 != nu) Error("nu mismatch\n");
	    Fortran(out2ind)(&l1, &l2, &l3, do2i, u0, v0, r, &K);
	    Fortran(sfmult)(&nu, &nv, ffsf, do2i, ffsg, &Nd);
	    cfix2y(0, 0, 1, ndegi, nu, nv, ffsg, ndeg, nul, nvl, ffsh);
	    if (sink->ndeg != ndeg) Error("ndeg mismatch\n");
	    Fortran(sfadd)(&nul, &nvl, ffsh, sink->ffsf, sink->ffsf, &Nd);
	    Msgf(("%.3f %.3f %.3f -- %.3f %.3f %.3f, %f %f %d\n", 
		  sink->pos[0], sink->pos[1], sink->pos[2],
		  cp->pos[0], cp->pos[1], cp->pos[2], 
		  ffsh[0].r, ffsh[0].i, nul));
	} else if (sink->isbody) {
	    /* This is where direct interactions are calculated */
	    const body *bp = source_vec[i]->ptr;
	    float stren;
	    double dr2;
	    double kx, kx_inv;
	    double s, c;
	    
	    sink->interactions++;
	    /* Msgf(("%s\n", PrintKey(bp->key))); */
	    /* Msgf(("%d interact\n", i)); */
	    /* Msgf(("%s\n", PrintBodyContents(bp))); */
	    result[i] = MAC_ACCEPT;
	    stren = bp->strength;
	    VVV(r, = sink->pos, - bp->pos);
	    dr2 = Dot(r, r)*K*K;
	    if (dr2 == 0.0) continue;
	    kx_inv = 1.0/sqrt(dr2);
	    kx = kx_inv * dr2;
#ifdef HAS_SINCOS
	    sincos(kx, &s, &c);
#else
	    s = sin(kx); c = cos(kx); 
#endif
	    sink->phi.r += s*kx_inv*stren;
	    sink->phi.i += c*kx_inv*stren;
	} else {
	    const body *bp = source_vec[i]->ptr;
	    VVV(r, = sink->pos, - bp->pos);
	    if (Dot(r, r) < 15*sink->sz*sink->sz) {
		/* Msgf(("%s\n", PrintKey(bp->key))); */
		/* Msgf(("%d splitting sink 1\n", i)); */
		result[i] = MAC_SPLIT_SINK;
		continue;
	    } else {
		double r0[3] = {0,0,0};
		double strength;
		Msgf(("Doing cell-body interaction\n"));
		sink->interactions++;
		result[i] = MAC_ACCEPT;

		Fortran(gendto)(&nu, &nv, do2o, u0, v0, r0, &K);
		strength = bp->strength;
		Fortran(genffsf)(&strength, &nu, &nv, ffsf, do2o, &Nd);
		if (l1 != nu) Error("nu mismatch\n");
		Fortran(out2ind)(&l1, &l2, &l3, do2i, u0, v0, r, &K);
		Fortran(sfmult)(&nu, &nv, ffsf, do2i, ffsg, &Nd);
		cfix2y(0, 0, 1, ndegi, nu, nv, ffsg, ndeg, nul, nvl, ffsh);
		Fortran(sfadd)(&nul, &nvl, ffsh, sink->ffsf, sink->ffsf, &Nd);
		Msgf(("%.3f %.3f %.3f -- %.3f %.3f %.3f, %f %f %f %f\n", 
		      sink->pos[0], sink->pos[1], sink->pos[2],
		      bp->pos[0], bp->pos[1], bp->pos[2], 
		      ffsh[0].r, ffsh[0].i, ffsh[1].r, ffsh[1].i));
	    }
	}
    }
    Free(ffsh);
    Free(ffsg);
    Free(ffsf);
}

