#include "physics_vrtx.h"
#include "tree.h"
#include "vop.h"

Timer_t VrtxTm;

void 
NlgNInherit(const Sink *from, Sink *to, hcell *pp){
    body *bp;

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
	bp = to->bp = pp->ptr;
	VS(Psi(bp), = (float)0.0);
	VS(Dstr(bp), = (float)0.0);
	VS(Vel(bp), = (float)0.0);
	MS(Gradvel(bp), = (float)0.0);
	Errsum(bp) = (float)0.;
	Errsum2(bp) = (float)0.;
	bp->nterms = 0;
    }
}

void NlgNMACv(Sink *sink, const hcell **srcs, int *results, int nsrc)
{
    int i;
    body *bp = sink->bp;
    Vxd(float pos);
    Vxd(float dr);
    float r2;

    StartTimer(&VrtxTm);
    if (bp == NULL) {
	for (i = 0; i < nsrc; i++) results[i] = MAC_SPLIT_SINK;
	return;
    }
    VxV(pos, = bp->pos);
    for (i = 0; i < nsrc; i++) {
	const hcell *source = srcs[i];
	if (Sub_Flags(source) == 0) {
	    InteractBody(source->ptr, bp, epsinv, nu);
	    results[i] = MAC_ACCEPT;
	} else {
	    cell *cp = (cell *)(source->ptr);
	    VxVxV(dr, = pos, - cp->pos);
	    r2 = Dotx(dr, dr);
	    if (r2 < cp->rcrit2) {
		results[i] = MAC_SPLIT_SRC;
	    } else {
		InteractCell(cp, bp);
		results[i] = MAC_ACCEPT;
	    }
	}
    }
    StopTimer(&VrtxTm);
}

