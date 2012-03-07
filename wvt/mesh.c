#include "ndim.h"
#include "stk.h"
#include "vop.h"
#include "tree.h"
#include "mesh.h"
#include "physics_sph.h"
#include "fastflpt.h"

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 40000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

Counter_t SPHCnt, SPHrej, nbrMACCnt;

static double dvtable; /* == 0.0001 ... */
static double invdvtable; /* == 10000.0 ... */
static double cnormk;
static double wij[MAX_INDEX];
static double grwij[MAX_INDEX];
static double fmass[MAX_INDEX];
static double fpoten[MAX_INDEX];
static int ndim;


void MeshInit(Meshbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
	      double max[NDIM], double delta[NDIM])
{
  /* Here you define the mesh structure! In principle, this could also be 
     an adaptive mesh (use larger squares on outskirts?).*/
    Stk s;
    Meshbody *q;
    double pos[NDIM];

    StkInitEz(&s);

    for(pos[0] = min[0]; pos[0] <= max[0]; pos[0] += delta[0]) {

	for(pos[1] = min[1]; pos[1] <= max[1]; pos[1] += delta[1]) {

	    for(pos[2] = min[2]; pos[2] <= max[2]; pos[2] += delta[2]) {

		q = StkPush(&s, sizeof(Meshbody));
		VV(q->pos, = pos);
/* 		q->nterms = 1; */
	    }

	}

    }

    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(Meshbody);
    q = StkBase(&s);
    *btabp = Realloc(q, *nobj * sizeof(Meshbody));
}


void MeshFixKeys(Meshbody *btab, int nobj, Key_t (*func)(const void *))
{
    Meshbody *btabend = btab+nobj;

    while(btab<btabend) {
	btab->key = func(btab);
	btab++;
    }

    return;
}


double MeshGetCost(const Meshbody *ptr)
{
    return 1.0;
}


Key_t MeshGetKeyFromStruct(const Meshbody *ptr)
{
    return ptr->key;
}

p
void MeshCofmFromDaugh(hcellptr hptr, hcellptr daughters[])
     /* This may not even be necessary */
{
/*     Meshbody *bp; */
/*     Meshcofmdata *dp; */
/*     Meshcofmdata *cmp; */

/*     int i; */

/*     cmp = hptr->ptr; */
/*     assert(cmp); */
/*     VS(cmp->pos, = 0.0); */
/*     cmp->ndaughters = 0; */

/*     for(i = 0; i < (1<<NDIM); ++i) { */
/* 	if (daughters[i] == NULL)  */
/* 	    continue; */

/* 	if (Sub_Flags(daughters[i]) == 0) { */
/* 	    bp = daughters[i]->ptr; */
/* 	    VV(cmp->pos, += bp->pos); */
/* 	    cmp->ndaughters++; */
/* 	} */
/* 	else { */
/* 	    dp = daughters[i]->ptr; */
/* 	    VV(cmp->pos, += dp->pos); */
/* 	    cmp->ndaughters += dp->ndaughters; */
/* 	} */
/*     } */

    return;
}


void MeshCellFromCofm(Meshcell *cp, Meshcofmdata *cmp) 
{
/*     VV(cp->pos, = cmp->pos); */
/*     cp->ndaughters = cmp->ndaughters; */
}


void InheritMesh(const Meshbody *from, Meshbody *to, hcell *pp)
{
    Meshbody *bp = pp->ptr;

    if ( to == NULL ) {
	bp->rho += from->rho;
    }
    else {
	VV(to->pos, = bp->pos);
	to->rho = (double)0.0;
    }

    return;
}


void Meshgate(Meshbody *sink, hcell **source_vec, int *result, int n)
{
    VxdV(const double pos_sink, = sink->pos);
    Vxd(double r);
    /* Vxd(double dv); */
    double extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    /* double projv; */
    double v2;
    double rij;
    double wtij, grwtij;
    double hmean11, hmean21;
    double dxx, dwdx, dgrwdx;
    int index;
    int nbrs = 0;
    double rhoi = (double)0.0;
    /* double divvi = (double)0.0; */
    double dr2;
    int interactions = 0;
      
    for (i = 0; i < n; i++) {
	const hcellptr source = source_vec[i];

	if (Sub_Flags(source)) {
	    const SPHcell *cp = source->ptr;
	    VxV(r, = cp->pos);	
	    extent_src = cp->bmax + cp->lap;
	    daughters = cp->daughters;
	} else {
	    bp = source->ptr;
	    VxV(r, = bp->pos);
	    extent_src = bp->h;
	    daughters = 1;
	}

	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r); /* == 400.0000000032623 */

	/* What about the h vs. 2*h thing */
	if ((rij = sqrtf_fast(dr2)) > extent_src /* + extent_sink */
	    || dr2 == (double)0.0) {
	    goto accept;
	} else if (daughters != 1) {
	    goto failed;
	}

	/* hmean11 = (double)2.0 / (h + bp->h); */
	hmean11 = (double)1.0/(bp->h);
	hmean21 = hmean11 * hmean11; /* == 0.01 */
	    
	v2 = dr2 * hmean21;
	index = v2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	dxx = v2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	wtij = (wij[index] + dwdx * dxx ) * hmean21 * hmean11;
	if (wtij < (double)0.0) Error("Negative wtij (macRho) = %g %g\n", wtij, hmean11); 
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;

	rhoi += bp->mass * wtij;

	/* velocity divergence times density */
	/* VxVV(dv, = bp->vel, - sink->vel); */
/* 	projv = grwtij * Dotx(dv, r) * recipsqrtf(dr2); */
/* 	divvi -= bp->mass * projv; */

	nbrs++;
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
      failed:
	result[i] = MAC_SPLIT_SRC;
    }

/*     sink->interactions += interactions; */
    sink->rho += rhoi;
/*     sink->nbrs += nbrs; */
/*     sink->drho_dt -= divvi; */
}
