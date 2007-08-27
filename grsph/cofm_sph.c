#include <math.h>
#include <stdio.h>
#include "tree.h"
#include "physics_sph.h"
#include "vop.h"
#include "chn.h"
#include "bigmalloc.h"
#include "Msgs.h"
#include "error.h"
#include "Assert.h"
#include "fastflpt.h"
#include "protos.h"

static float tol_stat;
static float three_tol;

void cofm_setup(float tol){
    tol_stat = tol;
    three_tol = 3.F/tol;
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
    cmp->mass = 0.F;
    VS(cmp->pos, = 0.F);
    cmp->B2 = 0.F;
    cmp->bmax = 0.F;
    cmp->lap = 0.F;
    cmp->ndaughters = 0;
    /* First get the cm of the new cell. */
    for(i=0; i<(1<<NDIM); i++){
	if( daughters[i] == NULL )
	    continue;
	if( Sub_Flags(daughters[i]) == 0 ){
	    bp = daughters[i]->ptr;
	    assert(bp);
	    cmp->ndaughters++;
	    dmass = bp->mass;
	    cmp->mass += dmass;
	    if (bp->h > cmp->lap)
	      cmp->lap = bp->h;
	    VV(cmp->pos, += dmass * bp->pos);
	}else{
	    dp = daughters[i]->ptr;
	    /* If this daughter doesn't exist, bail out now. */
	    assert(dp);
	    cmp->ndaughters += dp->ndaughters;
	    dmass = dp->mass;
	    cmp->mass += dmass;
	    if (dp->lap > cmp->lap)
	      cmp->lap = dp->lap;
	    VV(cmp->pos, += dmass * dp->pos);
	}
    }
    /* Divide out the total mass */
    if( cmp->mass != 0.F ){
	cmp->massinv = recipf(cmp->mass);
	VS(cmp->pos, *= cmp->massinv);
    }else{
	Error("Zero mass in BranchFromDaughters!\n");
    }
    /* Now loop again to pick up B2, etc.  */
    for(i=0; i<(1<<NDIM); i++){
	float tmp[NDIM];
	float tmpsq;

	if(daughters[i] == NULL)
	    continue;
	if( Sub_Flags(daughters[i]) == 0 ){
	    bp = daughters[i]->ptr;
	    VVV(tmp, = cmp->pos, - bp->pos);
	    dmass = bp->mass;
	    newbmax = 0.F;
	}else{
	    dp = daughters[i]->ptr;
	    VVV(tmp, = cmp->pos, - dp->pos);
	    dmass = dp->mass;
	    newbmax = dp->bmax;
	    cmp->B2 += dp->B2;
	}
	tmpsq = Dot(tmp, tmp);
	cmp->B2 += dmass * tmpsq;
	newbmax += sqrtf_fast(tmpsq);
	if( newbmax > cmp->bmax )
	    cmp->bmax = newbmax;
    }
    /* This is an alternative bound on bmax, which is sometimes tighter */
    /* than the cumulative bound computed above. */
    CELLCORNER(hptr->key, center, &cellsz);
    cellsz *= 0.5F;
    VS(center, += cellsz);
    VxVVS(dx, = cellsz+ fabs LPAREN cmp->pos, - center,  RPAREN);
    newbmax = sqrtf_fast(Dotx(dx, dx));
    cmp->bmax = ( newbmax < cmp->bmax ) ? newbmax : cmp->bmax;
    hptr->ptr = cmp;
}

static double a[6];		/* coef of error poly */
static void rcrit_poly(double r, double *value, double *deriv);
static double rtnewt(void (*funcd)(double, double *, double *),  
		    double x1, double xacc);

/* Turn the ptr from a cofmdata to a cell. */
void CellFromCofm(cell *cp, cofmdata *cmp){
    float B2, B3;
    float bmaxhalf, rcritmax;

    B2 = cmp->B2;

    if( B2 == 0.F ){
	/* Error("Does this happen?\n"); */
        /* q->rcrit = q->bmax; */
        cp->rcrit = cp->bmax;
    }else{
	bmaxhalf = cmp->bmax * 0.5F;
	rcritmax = bmaxhalf + sqrtf_fast(bmaxhalf*bmaxhalf + sqrtf_fast(three_tol*B2));
	if (!finite(rcritmax))
	  Error("Bad rcritmax, q->bmax = %g, B2 = %g\n", cmp->bmax, B2);
	B3 = B2 * sqrtf_fast(B2*cmp->massinv); 
	cp->B3 = 2.F*B3;
	if (!finite(B2) || !finite(B3) || !finite(cmp->bmax))
	  Error("Bad value B2 = %g, B3 = %g, bmax = %g\n", B2, B3, cmp->bmax);
	a[0] = 2.*B3;
	a[1] = -3. * B2;
	a[2] = 0.;
	a[3] = tol_stat*cmp->bmax*cmp->bmax;
	a[4] = -2. * tol_stat *cmp->bmax;
	a[5] = tol_stat;
	cp->rcrit =  rtnewt(rcrit_poly, rcritmax, .01*rcritmax);
    }
    cp->B2 = 3.F*B2;

    cp->mass = cmp->mass;
    VV(cp->pos, = cmp->pos);
    cp->daughters = cmp->ndaughters;
    cp->bmax = cmp->bmax;
    cp->lap = cmp->lap;
    Msgf(("Cell: %s\n", PrintCellContents(cp)));
}


/* Use doubles here to avoid catastrophe from roundoff. */
static void rcrit_poly(double r, double *value, double *deriv){
    /* Do we care that a[2] and da[1] are zero? */
    double p = a[5];
    double dp = 0.;
    int n = 5;

    /* See pg. 149 of Numerical Rec. */
    /* We could unroll it... */
    while(n>0) { dp = dp*r + p; p = p*r + a[--n]; }
    if (!finite(p) || !finite(dp))
      Error("Bad p or dp, p = %g, dp = %g, n = %d\n", p, dp, n);
    *value = p;
    *deriv = dp;
}

/* From Numerical Recipes, rtnewt.c (modified) */
/* It's even more dangerous than the version that NR says is too */
/* dangerous to use...No checking of bounds.  We might just run off */
/* to infinity... */
#define JMAX 20

static double
    rtnewt(void (*funcd)(double, double *, double *),
	   double x1, double xacc)
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

#undef JMAX
