#include "physics_vrtx.h"
#include "tree.h"
#include "vop.h"

Timer_t VrtxTm;

/* These should be dynamically extensible */
#define BVECSZ 1024
#define CVECSZ 2048

body **Bvec[BVECSZ];
cell **Cvec[CVECSZ];

void 
DLNlgNInherit(const Sink *from, Sink *to, hcell *pp)
{
    body *bp;

    if( to == NULL ){
	/* this is the last sink */
	int i;
	bp = pp->ptr;
	VS(Psi(bp), = (float)0.0);
	VS(Dstr(bp), = (float)0.0);
	VS(Vel(bp), = (float)0.0);
	MS(Gradvel(bp), = (float)0.0);
	Errsum(bp) = (float)0.;
	Errsum2(bp) = (float)0.;
	bp->nterms = 0;
	for (i = 0; i < from->bcnt; i++)
	  InteractBody(Bvec[i], bp, epsinv, nu);
	for (i = 0; i < from->ccnt; i++)
	  InteractCell(Cvec[i], bp);
	return;
    }

    if (Sub_Flags(pp)) {
	cell *cp = pp->ptr;
 	VV(to->pos, = cp->pos);
	to->bmax = cp->bmax;
    } else {
	bp = pp->ptr;
	VV(to->pos = bp->pos);
	to->bmax = 0.F;
    }
    if (from) {
	to->bcnt = from->bcnt;
	if (from->bcnt >= BVECSZ) Error("Bvec overflow\n");
	to->ccnt = from->ccnt;
	if (from->ccnt >= CVECSZ) Error("Cvec overflow\n");
    } else {
	to->bcnt = to->ccnt = 0;
    }
}

#define RcritFac ((float)4.0)	/* should be >= 2.0 */

/* RcritMAC with Don't Laugh-like traversal */
void DLNlgNMACv(Sink *sink, const hcell **srcs, int *results, int nsrc)
{
    int i;
    VxdV(float pos_sink, = sink->pos);
    float bmax = sink->bmax;
    Vxd(float dr);
    float r2, rcrit_bmax;
    cell *cp;

    StartTimer(&VrtxTm);
    for (i = 0; i < nsrc; i++) {
	const hcell *source = srcs[i];
	if (Sub_Flags(source) == 0) {
	    Bvec[bcnt++] = source->ptr;
	    results[i] = MAC_ACCEPT;
	} else {
	    cp = source->ptr;
	    rcrit_bmax = cp->rcrit + bmax;     /* can't use cp->rcrit2 */
	    VxVxV(dr, = pos_sink, - cp->pos);
	    r2 = Dotx(dr, dr);
	    if (r2 > rcrit_bmax*rcrit_bmax) {
		Cvec[ccnt++] = source->ptr;
		results[i] = MAC_ACCEPT;
	    } else if (RcritFac * bmax > rcrit_bmax) {
		results[i] = MAC_SPLIT_SINK;
		if (sink->bmax == 0.F) Error("Trying to split body\n");
	    } else {
		results[i] = MAC_SPLIT_SRC;
	    }
	}
    }
    StopTimer(&VrtxTm);
}

