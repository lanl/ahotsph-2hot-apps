#include <math.h>
#include "tree.h"
#include "Msgs.h"
#include "key.h"
#include "chn.h"
#include "physics_vrtx.h"
#include "physics_generic.h"
#include "vop.h"
#include "pqsort.h"
#include "Assert.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BTAB_EXTEND 1000
#define MSG_TAG 1325

typedef struct{
    float pos[NDIM];
    float strength[NDIM];
    int ident;
    void *next;
} rmshbody;

int ilog2(unsigned int num);
void exch_bounds(int n, Key_t hi, Key_t lo, Key_t *hiboundp, Key_t *loboundp);
static Key_t KeyOffset(Key_t key, int *ii);

/* Note: We reset all btab information here */

void
remesh(body **btabp, int *nobjp, float remesh_h, float remesh_min_str)
{
    float remesh_min_str2 = remesh_min_str*remesh_min_str;
    int nobj = *nobjp;
    body *btab = *btabp;
    body *end;
    Chn remeshed_stk, *s;
    int rmsh_level;
    body *p, *q;
    rmshbody *b;
    int nmesh;
    tree_t thetree;
    float center[NDIM];
    float corner[NDIM], corner0[NDIM];
    float cellsz;
    float rmin[NDIM], rmax[NDIM];
    float r2, dx[NDIM];
    hcell *pp;
    Key_t key, tkey;
    int ii[NDIM];
    float weight, wgt[NDIM];
    float u, d0[NDIM];
    float vol, inverse_h;
    rmshbody bi, *bprev;	/* bi.next will point to first rmshbody */
    int btab_n;
    sortresult_t srt;
    Key_t last_k, hibound, lobound;
    
    FindBbox(btab, nobj, rmin, rmax);
    nmesh = 2 + FixRsize(rmin, rmax)/remesh_h; /* add 2 for overlap */
    rmsh_level = 1+ilog2(nmesh);
    nmesh = 1 << rmsh_level;

    VVVS(center, = LPAREN rmax, + rmin, RPAREN*0.5);
    VV(rmin, = -0.5*nmesh*remesh_h + center);
    VV(rmax, =  0.5*nmesh*remesh_h + center);
    /* It would be better to use FixRsizeExact here */
    FixRsizeExact(rmin, rmax);
    FixKeys(btab, nobj, GetKey);

    thetree.htab = Calloc(HASH_TABLE_SIZE, sizeof(hcellptr));
    thetree.ndim = NDIM;
    thetree.hash_mask = HASH_MASK;
    ChnInit(&thetree.hcellchn, sizeof(hcell), 4096, Realloc_f);

    CellCorner(KeyLshift(KeyInt(1),rmsh_level*thetree.ndim), corner, &cellsz);
    vol = cellsz*cellsz*cellsz/(4.0*M_PI);
    inverse_h = 1.0/cellsz;

    s = &remeshed_stk;
    ChnInit(s, sizeof(rmshbody), nobj/8, Realloc_f);

    /* This uses a temporary array of about nobj * sizeof(rmshbody) */
    /* It could be reduced somewhat by going through btab in reverse */
    /* and realloc'ing to smaller sizes as we use the data */
    
    bprev = &bi;		/* linked list of remeshed bodies */
    for (p = btab; p < btab+nobj; p++) {
	key = KeyRshift(GetKey(p), KEYBITS-1-rmsh_level*thetree.ndim);
	CellCorner(key, corner0, &cellsz);
	VVV(d0, = p->pos, - (float)0.5 * cellsz - corner0);
	for (ii[0] = -1; ii[0] <= 1; ii[0]++) {
	    u = inverse_h * d0[0];
	    if (ii[0] == -1) wgt[0] = 0.5*u*(u-1.0);
	    else if (ii[0] == 1) wgt[0] = 0.5*u*(u+1.0);
	    else wgt[0] = 1.0-u*u;
	    for (ii[1] = -1; ii[1] <= 1; ii[1]++) {
		u = inverse_h * d0[1];
		if (ii[1] == -1) wgt[1] = 0.5*u*(u-1.0);
		else if (ii[1] == 1) wgt[1] = 0.5*u*(u+1.0);
		else wgt[1] = 1.0-u*u;
		for (ii[2] = -1; ii[2] <= 1; ii[2]++) {
		    u = inverse_h * d0[2];
		    if (ii[2] == -1) wgt[2] = 0.5*u*(u-1.0);
		    else if (ii[2] == 1) wgt[2] = 0.5*u*(u+1.0);
		    else wgt[2] = 1.0-u*u;
		    weight = wgt[0] * wgt[1] * wgt[2];
		    assert(weight <= 1.0 && weight >= -0.125);
		    tkey = KeyOffset(key, ii);
		    CellCorner(tkey, corner, &cellsz);
		    if ((pp = Find(&thetree, tkey))) {
			b = pp->ptr;
			b->ident = KeyAndInt(tkey, ~0);
			VV(b->strength, += weight * p->strength);
		    } else {
			bprev->next = b = ChnAlloc(s);
			Enter(&thetree, tkey, b, 0);
			VV(b->strength, = weight * p->strength);
			VV(b->pos, = (float)0.5 * cellsz + corner);
			b->next = NULL;
			bprev = b;
		    }
		    /* sanity check, not required */
		    VVV(dx, = corner, - corner0);
		    r2 = Dot(dx, dx);
		    assert(r2*inverse_h*inverse_h <= 3.0001);
		}
	    }
	}
    }
    Free(thetree.htab);
    ChnTerminate(&thetree.hcellchn);

    /* re-create btab from linked list of remeshed lattice points */

    btab_n = nobj;
    nobj = 0;
    p = btab;
    for (b = bi.next; b;  b = b->next) {
	VV(p->pos, = b->pos);
	VV(p->strength, = b->strength);
	/* All particles get same vol */
	p->vol = vol;
	p->ident = b->ident;
	p++;
	nobj++;
	if (nobj >= btab_n) {
	    btab_n += BTAB_EXTEND;
	    btab = Realloc(btab, btab_n*sizeof(body));
	    p = btab + nobj;
	}
    }
    ChnTerminate(s);

    /* re-sort the list of bodies */

    FindBbox(btab, nobj, rmin, rmax);
    FixRsize(rmin, rmax);
    FixKeys(btab, nobj, GetKey);

    pqsortsetup_order(&srt, btab, nobj, sizeof(body), .01, 1, Realloc_f);
    btab = pqsort(&srt, (pq_wgtproto)UnityCost, (pq_keyproto)GetKey);
    nobj = srt.nobj;

    /* Sum pieces that were assigned across a domain boundary */
    
    end = btab + nobj;
    nobj = 0;
    for (q = p = btab; p < end; p++) {
	if (q != p) {
	    VV(q->pos, = p->pos);
	    VV(q->strength, = p->strength);
	    q->ident = p->ident;
	    q->key = p->key;
	}
	while (KeyEQ(q->key, (p+1)->key)) {
	    assert(q->pos[0] == (p+1)->pos[0]);
	    assert(q->pos[1] == (p+1)->pos[1]);
	    assert(q->pos[2] == (p+1)->pos[2]);
	    VV(q->strength, += (p+1)->strength);
	    p++;
	}
	if (Dot(q->strength, q->strength) > remesh_min_str2) {
	    nobj++;
	    q++;
	}
    }

    /* Check if there is a point that straddles a processor boundary */
    
    last_k = (btab+nobj-1)->key;
    exch_bounds(NDIM, last_k, btab->key, &hibound, &lobound);
    if (KeyEQ(last_k, hibound)) {
	/* I'm not sure if this presents a problem or not */
	Shout("Key %s straddles procs\n", PrintKey(last_k));
    }

    *btabp = Realloc(btab, nobj*sizeof(body));
    *nobjp = nobj;
}


/* This will fail if the most significant bit of the key is more than 31 */
/* This was written somewhere else in a better fashion, but I can't find it */
static Key_t
KeyOffset(Key_t key, int *ii)
{
    Key_t outkey;
    unsigned int dim, bit;
    
    outkey = key;
    for (dim = 0; dim < NDIM; dim++) {
	bit = 1 << dim;
	if (ii[dim] == 1) {
	    while (KeyAndInt(outkey, bit)) {
		outkey = KeyXOR(outkey, KeyInt(bit));
		bit <<= NDIM;
	    }
	    outkey = KeyOrInt(outkey, bit);
	} else if (ii[dim] == -1) {
	    while (KeyAndInt(KeyNot(outkey), bit)) {
		outkey = KeyOrInt(outkey, bit);
		bit <<= NDIM;
	    }
	    outkey = KeyXOR(outkey, KeyInt(bit));
	} else {
	    assert(ii[dim] == 0);
	}
    }
    return outkey;
}
	


#if 0 /* This is the loop for Nearest Grid Point Interpolation */

    for (p = btab; p < btab+nobj; p++) {
	key = KeyRshift(GetKey(p), KEYBITS-1-rmsh_level*thetree.ndim);
	if ((pp = Find(&thetree, key))) {
	    b = pp->ptr;
	    b->ident = KeyAndInt(key, ~0);
	    VV(b->strength, += p->strength);
	} else {
	    b = StkPush(s, sizeof(rmshbody));
	    Enter(&thetree, key, b, 0);
	    VV(b->strength, = p->strength);
	    CellCorner(key, corner, &cellsz);
	    VV(b->pos, = (float)0.5 * cellsz + corner);
	}
    }
#endif
