#include <math.h>
#include <stdio.h>
#include "tree.h"
#include "physics.h"
#include "vop.h"
#include "chn.h"
#include "bigmalloc.h"
#include "Msgs.h"
#include "error.h"
#include "Assert.h"
#include "fastflpt.h"
#include "protos.h"

int MACtype = BMAX_MAC;		/* default */
float Tol;
float invTol;
float RelTol;
float invRelTol;

void SetupCofm(int type, float tol, float rel_tol)
{
    MACtype = type;
    Tol = tol;
    invTol = 1.0/tol;
    RelTol = rel_tol;
    invRelTol = 1.0/rel_tol;
}

void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]){
    int i;
    cofmdata *dp;
    cofmdata *cmp;
    body *bp;
    float dmass;
    float newbmax;
    float center[NDIM], cellsz;
    Vxd(float dx);

    assert(Sub_Flags(hptr));

    cmp = hptr->ptr;
    assert(cmp);
    cmp->mass = 0.;
    VS(cmp->pos, = 0.);
    cmp->B2 = 0.;
    cmp->bmax = 0.;
    cmp->ndaughters = 0;

    /* First get the cm of the new cell. */
    for(i=0; i<(1<<NDIM); i++){
	if (daughters[i] == NULL)
	  continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    dmass = bp->mass;
	    cmp->mass += dmass;
#ifdef SPH_GRAV
	    /* Gabe+Steven Dirty trick: Fix bmax to be 2h for the case
	       that you are on the lowest level to make sure that
	       rcrit is at least as big as the maximum of all h within
	       the cell. */
	    cmp->bmax = (cmp->bmax > 2.0*fabs(bp->h)) ? 
		cmp->bmax : 2.0*fabs(bp->h);
#endif
	    VV(cmp->pos, += dmass * bp->pos);
	    cmp->ndaughters++;
	} else {
	    dp = daughters[i]->ptr;
	    dmass = dp->mass;
	    cmp->mass += dmass;
	    VV(cmp->pos, += dmass * dp->pos);
	    cmp->ndaughters += dp->ndaughters;
	} 
    }
    /* Divide out the total mass */
    if (cmp->mass != (float)0.) {
	cmp->massinv = recipf(cmp->mass);
	VS(cmp->pos, *= cmp->massinv);
    } else {
	Error("Zero mass in BranchFromDaughters!\n");
    }
    /* Now loop again to pick up B2, etc.  */
    for (i=0; i<(1<<NDIM); i++) {
	float tmp[NDIM];
	float tmpsq;

	if(daughters[i] == NULL)
	    continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    VVV(tmp, = cmp->pos, - bp->pos);
	    dmass = bp->mass;
	    newbmax = (float)0.;
	} else {
	    dp = daughters[i]->ptr;
	    VVV(tmp, = cmp->pos, - dp->pos);
	    dmass = dp->mass;
	    newbmax = dp->bmax;
	    cmp->B2 += dp->B2;
	}
	tmpsq = Dot(tmp, tmp);
	if (tmpsq != 0.F) {
	    /* avoid doing a sqrtf_fast(0).  and don't bother multiplying
	       by and adding zero either */
	    cmp->B2 += dmass * tmpsq;
	    newbmax += sqrtf_fast(tmpsq);
	}
	if (newbmax > cmp->bmax)
	  cmp->bmax = newbmax;
    }
    /* This is an alternative bound on bmax, which is sometimes tighter */
    /* than the cumulative bound computed above. */
#ifndef SPH_GRAV
    CELLCORNER(hptr->key, center, &cellsz);
    cmp->sz = cellsz;		/* for pure Barnes-But */
    cellsz *= (float)0.5;
    VS(center, += cellsz);
    VxVVS(dx, = cellsz+ fabs LPAREN cmp->pos, - center,  RPAREN);
    newbmax = sqrtf_fast(Dotx(dx, dx));
    cmp->bmax = (newbmax < cmp->bmax) ? newbmax : cmp->bmax;
    /* Gabe&Steven Dirty Fix: Let's leave this out for SPH particles,
       as the new bmax has no clue about particle SPHs, and we don't
       want to overwrite that information. */
#endif

    hptr->ptr = cmp;
}

static double a[6];		/* coef of error poly */
static void rcrit_poly(double r, double *value, double *deriv);
static double rtnewt(void (*funcd)(double, double *, double *),  
		    double x1, double xacc);

/* Turn the ptr from a cofmdata to a cell. */
void CellFromCofm(cell *cp, cofmdata *cmp)
{
    cp->mass = cmp->mass;
    VV(cp->pos, = cmp->pos);
    cp->bmax = cmp->bmax;
    if (MACtype == AREL_MAC) {
	float abs_rcrit;
	float rel_rcrit;
	float B3, bmaxhalf, rcritmax;
	float B2 = cmp->B2;

	bmaxhalf = cmp->bmax * (float)0.5;
	rcritmax = bmaxhalf + sqrtf_fast(bmaxhalf*bmaxhalf
					 + sqrtf_fast((float)3.*invTol*B2));
	if (!isfinite(rcritmax))
	  Error("Bad rcritmax, q->bmax = %g, B2 = %g\n", cmp->bmax, B2);
	if (B2 == (float)0.0) Error("B2 is zero\n");
	B3 = B2 * sqrtf_fast(B2*cmp->massinv); 
	if (!isfinite(B2) || !isfinite(B3) || !isfinite(cmp->bmax))
	  Error("Bad value B2 = %g, B3 = %g, bmax = %g\n", B2, B3, cmp->bmax);
	a[0] = 2.*B3;
	a[1] = -3. * B2;
	a[2] = 0.;
	a[3] = Tol*cmp->bmax*cmp->bmax;
	a[4] = -2. * Tol * cmp->bmax;
	a[5] = Tol;
	abs_rcrit =  rtnewt(rcrit_poly, rcritmax, .01*rcritmax);
	abs_rcrit += 0.01*rcritmax;

	rcritmax = cmp->bmax + sqrtf_fast((float)3.*invRelTol*B2*cmp->massinv);
	a[0] = 2.*B3;
	a[1] = -3. * B2 + RelTol * cmp->bmax*cmp->bmax * cmp->mass;
	a[2] = -2. * RelTol * cmp->bmax * cmp->mass;
	a[3] = RelTol * cmp->mass;
	a[4] = a[5] = 0.0;	/* should write rcrit_poly3 if it matters */
	rel_rcrit =  rtnewt(rcrit_poly, rcritmax, .01*rcritmax);
	rel_rcrit += 0.01*rcritmax;

	cp->rcrit = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
    } else if(MACtype == BH_MAC) cp->rcrit = cmp->sz*invTol;
    else if (MACtype == BMAX_MAC) cp->rcrit = cmp->bmax*invTol;
    else Error("Bad MAC type (%d)\n", MACtype);

    cp->daughters = cmp->ndaughters;
    Msgf(("Cell: %s\n", PrintCellContents(cp)));
}

/* Use doubles here to avoid catastrophe from roundoff. */
static void rcrit_poly(double r, double *value, double *deriv)
{
    /* Do we care that a[2] and da[1] are zero? */
    double p = a[5];
    double dp = 0.;
    int n = 5;

    /* See pg. 149 of Numerical Rec. */
    /* We could unroll it... */
    while(n>0) { dp = dp*r + p; p = p*r + a[--n]; }
    if (!isfinite(p) || !isfinite(dp))
      Error("Bad p or dp, p = %g, dp = %g, n = %d\n", p, dp, n);
    *value = p;
    *deriv = dp;
}

/* From Numerical Recipes, rtnewt.c (modified) */
/* It's even more dangerous than the version that NR says is too */
/* dangerous to use...No checking of bounds.  We might just run off */
/* to infinity... */
#define JMAX 40

static double
rtnewt(void (*funcd)(double, double *, double *), double x1, double xacc)
{
    int j;
    double df,dx,f,rtn;
    
    rtn=x1;
    for (j=1;j<=JMAX;j++) {
	(*funcd)(rtn,&f,&df);
	dx=f/df;
	Msg("rcrit", ("f(%g)=%g, dx=%g\n", rtn, f, -dx));
	rtn -= dx;
	if (fabs(dx) < xacc) return rtn;
    }
    Error("Maximum number of iterations exceeded in RTNEWT\n");
}
