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

void SPHSetupCofm(int type, float tol, float rel_tol)
{
}

void SPHCofmFromDaugh(hcellptr hptr, hcellptr daughters[]){
    int i;
    SPHcofmdata *dp;
    SPHcofmdata *cmp;
    SPHbody *bp;
    float dmass;
    float newbmax;
    float center[NDIM], cellsz;
    Vxd(float dx);

    assert(Sub_Flags(hptr));

    cmp = hptr->ptr;
    assert(cmp);
    cmp->mass = 0.;
    VS(cmp->pos, = 0.);
    VS(cmp->vel, = 0.);
    cmp->u = 0.;
    cmp->ident = 0;
    cmp->ifleos = 0.;
    cmp->abar = 0;
    cmp->temp = 0;
    cmp->ye = 0;
    cmp->xp = 0;
    cmp->xn = 0;
    cmp->u2 = 0;
    cmp->ynue = 0;
    cmp->ynueb = 0;
    cmp->ynux = 0;
    cmp->unue = 0;
    cmp->unueb = 0;
    cmp->unux = 0;
    cmp->ufreez = 0;
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
	    VV(cmp->vel, += dmass * bp->vel);
	    cmp->ifleos += bp->ifleos;
	    if (bp->ident > cmp->ident) cmp->ident = bp->ident;
	    cmp->ndaughters++;
	    cmp->u += dmass * bp->u;
	    cmp->abar += dmass * bp->abar;
	    cmp->temp += dmass * bp->temp;
	    cmp->ye += dmass * bp->ye;
	    cmp->xp += dmass * bp->xp;
	    cmp->xn += dmass * bp->xn;
	    cmp->u2 += dmass * bp->u2;
	    cmp->ynue += dmass * bp->ynue;
	    cmp->ynueb += dmass * bp->ynueb;
	    cmp->ynux += dmass * bp->ynux;
	    cmp->unue += dmass * bp->unue;
	    cmp->unueb += dmass * bp->unueb;
	    cmp->unux += dmass * bp->unux;
	    cmp->ufreez += dmass * bp->ufreez;
	} else {
	    dp = daughters[i]->ptr;
	    dmass = dp->mass;
	    cmp->mass += dmass;
	    if (dp->lap > cmp->lap)
	      cmp->lap = dp->lap;
	    VV(cmp->pos, += dmass * dp->pos);
	    VV(cmp->vel, += dmass * dp->vel);
	    cmp->ifleos += dp->ifleos;
	    if (dp->ident > cmp->ident) cmp->ident = dp->ident;
	    cmp->ndaughters += dp->ndaughters;
	    cmp->u += dmass * dp->u;
	    cmp->abar += dmass * dp->abar;
	    cmp->temp += dmass * dp->temp;
	    cmp->ye += dmass * dp->ye;
	    cmp->xp += dmass * dp->xp;
	    cmp->xn += dmass * dp->xn;
	    cmp->u2 += dmass * dp->u2;
	    cmp->ynue += dmass * dp->ynue;
	    cmp->ynueb += dmass * dp->ynueb;
	    cmp->ynux += dmass * dp->ynux;
	    cmp->unue += dmass * dp->unue;
	    cmp->unueb += dmass * dp->unueb;
	    cmp->unux += dmass * dp->unux;
	    cmp->ufreez += dmass * dp->ufreez;
	} 
    }
    /* Divide out the total mass */
    if (cmp->mass != (float)0.) {
	cmp->massinv = recipf(cmp->mass);
	VS(cmp->pos, *= cmp->massinv);
	VS(cmp->vel, *= cmp->massinv);
	cmp->u *= cmp->massinv;
	cmp->abar *= cmp->massinv;
	cmp->temp *= cmp->massinv;
	cmp->ye *= cmp->massinv;
	cmp->xp *= cmp->massinv;
	cmp->xn *= cmp->massinv;
	cmp->u2 *= cmp->massinv;
	cmp->ynue *= cmp->massinv;
	cmp->ynueb *= cmp->massinv;
	cmp->ynux *= cmp->massinv;
	cmp->unue *= cmp->massinv;
	cmp->unueb *= cmp->massinv;
	cmp->unux *= cmp->massinv;
	cmp->ufreez *= cmp->massinv;
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
	    newbmax = (float)0.;
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
    cellsz *= (float)0.5;
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
    VV(cp->vel, = cmp->vel);
    cp->ident = cmp->ident;
    cp->ifleos = cmp->ifleos/cmp->ndaughters;
    cp->bmax = cmp->bmax;
    cp->daughters = cmp->ndaughters;
    cp->lap = cmp->lap;
    cp->u = cmp->u;
    cp->abar = cmp->abar;
    cp->temp = cmp->temp;
    cp->ye = cmp->ye;
    cp->xp = cmp->xp;
    cp->xn = cmp->xn;
    cp->u2 = cmp->u2;
    cp->ynue = cmp->ynue;
    cp->ynueb = cmp->ynueb;
    cp->ynux = cmp->ynux;
    cp->unue = cmp->unue;
    cp->unueb = cmp->unueb;
    cp->unux = cmp->unux;
    cp->ufreez = cmp->ufreez;
    Msgf(("Cell: %s\n", PrintSPHCellContents(cp)));
}

