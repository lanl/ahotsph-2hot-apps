#include "physics_panel.h"
#include "tree.h"
#include "vop.h"

void 
NlgNInherit(const Sink *from, Sink *to, hcell *pp){
    body *bp;

    if (to == NULL) return;

    /* If we're terminal, then copy bp. */
    if( Sub_Flags(pp) ){
	to->bp = NULL;
    }else{
	bp = to->bp = pp->ptr;
	VS(Vel(bp), = (float)0.0);
	Phi(bp) = (float)0.0;
    }
}

void 
NlgNSinkToBody(const Sink *sink, body *to){
    return;
}

int 
NlgNMAC(Sink *sink, hcell *source){
    Vxd(float dr);
    cell *cp;

    if( sink->bp == NULL ){
	return MAC_SPLIT_SINK;
    }
    if( Sub_Flags(source) == 0 ){
	Binter(sink->bp, (body *)(source->ptr));
	return MAC_ACCEPT;
    }
    cp = (cell *)(source->ptr);

    VxVV(dr, = sink->bp->pos, - cp->pos);
    dr2 = Dotx(dr, dr);
    if( dr2 < cp->rcrit2 )
	return MAC_SPLIT_SRC;

    Cinter(sink->bp, cp);
    return MAC_ACCEPT;
}

void
NlgNMACv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    Vxd(float dr);
    cell *cp;
    int i;

    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];

	if (sink->bp == NULL) {
	    result[i] = MAC_SPLIT_SINK;
	} else if (Sub_Flags(source) == 0) {
	    Binter(sink->bp, (body *)(source->ptr));
	    result[i] = MAC_ACCEPT;
	} else {
	    cp = (cell *)(source->ptr);

	    VxVV(dr, = sink->bp->pos, - cp->pos);
	    dr2 = Dotx(dr, dr);
	    if (dr2 < cp->rcrit2) {
		result[i] = MAC_SPLIT_SRC;
	    } else {
		Cinter(sink->bp, cp);
		result[i] = MAC_ACCEPT;
	    }
	}
    }
}
