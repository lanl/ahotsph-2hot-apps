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

void SPHSetupCofm(int type, double tol, double rel_tol)
{
}

void SPHCofmFromDaugh(hcellptr hptr, hcellptr daughters[]){
    int i;
    SPHcofmdata *dp;
    SPHcofmdata *cmp;
    SPHbody *bp;
    double dmass;
    double newbmax;
    double center[NDIM], cellsz;
    Vxd(double dx);

    assert(Sub_Flags(hptr));

    cmp = hptr->ptr;
    assert(cmp);
    cmp->mass = 0.;
    VS(cmp->pos, = 0.);
    cmp->bmax = 0.;
    cmp->lap = 0.;
    cmp->ndaughters = 0;

    /* First get the cm of the new cell. */
    for(i=0; i<(1<<NDIM); i++){
	if (daughters[i] == NULL)
	  continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    dmass = bp->mass;
	    cmp->mass += dmass;
	    if (bp->h > cmp->lap)
	      cmp->lap = bp->h;
	    VV(cmp->pos, += dmass * bp->pos);
	    cmp->ndaughters++;
	} else {
	    dp = daughters[i]->ptr;
	    dmass = dp->mass;
	    cmp->mass += dmass;
	    if (dp->lap > cmp->lap)
	      cmp->lap = dp->lap;
	    VV(cmp->pos, += dmass * dp->pos);
	    cmp->ndaughters += dp->ndaughters;
	} 
    }
    /* Divide out the total mass */
    if (cmp->mass != (double)0.) {
	cmp->massinv = recipf(cmp->mass);
	VS(cmp->pos, *= cmp->massinv);
    } else {
	Error("Zero mass in BranchFromDaughters!\n");
    }
    /* Now loop again to pick up B2, etc.  */
    for (i=0; i<(1<<NDIM); i++) {
	double tmp[NDIM];
	double tmpsq;

	if(daughters[i] == NULL)
	    continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    VVV(tmp, = cmp->pos, - bp->pos);
	    newbmax = (double)0.;
	} else {
	    dp = daughters[i]->ptr;
	    VVV(tmp, = cmp->pos, - dp->pos);
	    newbmax = dp->bmax;
	}
	tmpsq = Dot(tmp, tmp);
	if (tmpsq != 0.F) {
	    /* avoid doing a sqrtf_fast(0).  and don't bother multiplying
	       by and adding zero either */
	    newbmax += sqrtf_fast(tmpsq);
	}
	if (newbmax > cmp->bmax)
	  cmp->bmax = newbmax;
    }
    /* This is an alternative bound on bmax, which is sometimes tighter */
    /* than the cumulative bound computed above. */
    CellCorner(hptr->key, center, &cellsz);
    cmp->sz = cellsz;		/* for pure Barnes-But */
    cellsz *= (double)0.5;
    VS(center, += cellsz);
    VxVVS(dx, = cellsz+ fabs LPAREN cmp->pos, - center,  RPAREN);
    newbmax = sqrtf_fast(Dotx(dx, dx));
    cmp->bmax = (newbmax < cmp->bmax) ? newbmax : cmp->bmax;
    hptr->ptr = cmp;
}

/* Turn the ptr from a cofmdata to a cell. */
void SPHCellFromCofm(SPHcell *cp, SPHcofmdata *cmp)
{
    cp->mass = cmp->mass;
    VV(cp->pos, = cmp->pos);
    cp->bmax = cmp->bmax;
    cp->daughters = cmp->ndaughters;
    cp->lap = cmp->lap;
    Msgf(("Cell: %s\n", PrintSPHCellContents(cp)));
}

