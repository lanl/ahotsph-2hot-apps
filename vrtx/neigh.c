/* Auxiliary routines to do tree-based neighbor finding.  */

#include <math.h>
#include "physics_vrtx.h"
#include "tree.h"
#include "neigh.h"
#include "vop.h"
#include "Msgs.h"
#include "timers.h"

Counter_t NfindTestsCnt;
Counter_t NfindAcceptsCnt;

#ifndef SQRT3
#define SQRT3 1.7320508F
#endif

void NeighCofmFromDaugh(hcell *hptr, hcell **dlist){
    ncofm *cofmp = hptr->ptr;
    /* I could be compulsive and do a good job of bmax-finding, or I
       can just be lazy.  We could avoid the cofm data structure
       altogether if CellFromCofm knew the hcell it was dealing
       with... */
    /* DANGER! DANGER!  CellCorner uses static Rsize and Rmin values
       that are left over from previous calls to FixRsize, which may well
       have been dealing with a completely different tree!  Not to
       mention the compiled-in value of NDIM, which might be wrong.
       physics_generic.[ch] really needs to be cleaned up! */
    CellCorner(hptr->key, cofmp->center, &cofmp->size);
    cofmp->size *= 0.5F;
    VS(cofmp->center, += cofmp->size);
}

void NeighCellFromCofm(ncell *cellp, ncofm *cofmp){
    float sz;
    VV(cellp->pos, = cofmp->center);
    sz = SQRT3 * cofmp->size + kc; /* kernel cutoff */
    cellp->r2 = sz*sz;
}

/* These look a lot like the code in mac.c. */
void
NeighInherit(const nsink *from, nsink *to, hcell *pp){
    if( to == NULL ){
	/* this is the last sink */
	/* we would normally call 'sink-to-body', but that does nothing */
	/* in this case anyway. */
	return;
    }

    /* If we're terminal, then copy bp. */
    if( Sub_Flags(pp) ){
	to->bp = NULL;
    }else{
	to->bp = pp->ptr;
	/* Zero all accumulators. */
	VS(Omegat(to->bp), = (float)0.0);
    }
}

static int 
NeighMAC(nsink *sink, const hcell *source){
    Vxd(float dr);
    ncell *cp;
    body *bp;

    if( sink->bp == NULL ){
	return MAC_SPLIT_SINK;
    }
    if( Sub_Flags(source) == 0 ){
	bp = (body *)(source->ptr);
	OmegatBody(bp, sink->bp, epsinv);
	return MAC_ACCEPT;
    }
    cp = (ncell *)(source->ptr);

    VxVV(dr, = sink->bp->pos, - cp->pos);
    dr2 = Dotx(dr, dr);
    IncrCounter(&NfindTestsCnt);
    if( dr2 < cp->r2 ){
	/* At least part of the cell may be close enough to contain
           neighbors. */
	/* In theory, it's possible to streamline this when we are
	   completely inside the kernel-cutoff. */
	Msgf(("Mac: dr=%g, cp->dr = %g\n", sqrt(dr2), sqrt(cp->r2)));
	return MAC_SPLIT_SRC;
    }

    /* There's nothing to do if we're far enough away... */
    return MAC_ACCEPT;
}

void 
NeighMACv(nsink *sink, const hcell **srcs, int *results, int nsrc){
    while(nsrc--){
	*results++ = NeighMAC(sink, *srcs++);
    }
}

