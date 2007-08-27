/* #define SHOW_WGTS */
/* SetupDecomp() figures out a way to assign every item to a processor */
/* The assignments are available by using DestDecomp() */

#include <stdlib.h>
#include "Assert.h"
#include "key.h"
#include "decomp.h"
#include "bigmalloc.h"
#include "mpmy.h"
#include "timers.h"
#include "stk.h"
#include "Msgs.h"
#include "gc.h"
#include "protos.h"

Timer_t DecompTm;
Timer_t DecompWaitTm;
Timer_t DecompCommTm;

static Key_t (*getkey_s)(const void *);
static float (*weight_s)(const void *);
static Key_t *decomptab;
static int save_decomp;
#ifdef SHOW_WGTS
static float *decomp_wgt;
#endif

#define CLEAR 0
#define SAVE 1
#define SET 2

int
MPMY_Doc(void)
{
    int doc = ilog2(MPMY_Nproc());
    if (MPMY_Nproc() != 1 << doc)
      doc++;			/* for non power-of-two sizes */
    return doc;
}

int
MPMY_PowOf2(void)
{
    int doc = ilog2(MPMY_Nproc());
    return (MPMY_Nproc() == 1 << doc);
}

static int
cmpkey2(const void *a, const void *b){
    Key_t ka;
    Key_t kb;
    ka = *(Key_t *)a;
    kb = *(Key_t *)b;
    if( KeyGT(ka, kb) )
	return 1;
    else if( KeyLT(ka, kb) )
	return -1;
    else
	return 0;
}

/* This technology needs reworking, with the decomp stored in a sortresult_t */
/* As it is, you can only keep one decomposition, which works for now */
void *
SaveDecomp(void)
{
    save_decomp = SAVE;
    decomptab = Malloc(MPMY_Nproc() * sizeof(Key_t));
    return decomptab;
}

void
SetDecomp(void *ptr)
{
    if (ptr == 0) return;
    save_decomp = SET;
    decomptab = ptr;
}

void
ClearDecomp(void *ptr)
{
    save_decomp = CLEAR;
    Free(ptr);
}


void
SetupDecomp(sortresult_t *decompp, 
	    float (*weight)(const void *), Key_t (*getkey)(const void *))
{
    char *b;
    int size = decompp->size;
    int nobj = decompp->nobj;
    char *data = decompp->data;
    Key_t *key_n, tmp;
    float wtfac;
    int nin;
    Stk ostk;
    float total_wgt;
    struct {
	Key_t key;
	float n;
    } *keydata, *kn;
    int i;

    getkey_s = getkey;
    weight_s = weight;

    Msgf(("SetupDecomp: starting in mode %d\n", save_decomp));    
    if (save_decomp == SET) {
	Msgf(("SetupDecomp: save decomp is set, returning\n"));
	return;
    }
    StartTimer(&DecompTm);
    total_wgt = 0.0;
    keydata = Malloc(nobj * sizeof(*keydata));
    i = 0;
    for(b = data; b < data + nobj * size; b += size) {
	float wt = weight(b);
	total_wgt += wt;
	keydata[i].key = getkey(b);
	keydata[i++].n = wt;
    }
    StopTimer(&DecompTm);
    MPMY_Combine(&total_wgt, &total_wgt, 1, MPMY_FLOAT, MPMY_SUM);
    StartTimer(&DecompTm);
    wtfac = total_wgt/(MPMY_Nproc()*400);
    Msgf(("total weight %f, wtfac %f\n", total_wgt, wtfac));
    total_wgt = 0.0;
    StkInitEz(&ostk);
    Msgf(("sorting %d keys in SetupDecomp\n", nobj));
    if (nobj) qsort(keydata, nobj, sizeof(*keydata), cmpkey2);
    for(kn = keydata; kn < keydata + nobj; kn++) {
	total_wgt += kn->n;
	if (total_wgt > wtfac) {
	    StkPushType(&ostk, kn->key, Key_t);
	    total_wgt = 0.0;
	}
    }
    Free(keydata);
    Msgf(("Sending %d bytes\n", StkSz(&ostk)));
    StopTimer(&DecompTm);
    nin = MPMY_NGather(StkBase(&ostk), StkSz(&ostk), MPMY_CHAR, 
		       (void **)&key_n, 0);
    StartTimer(&DecompTm);
    nin /= sizeof(Key_t);
    StkTerminate(&ostk);

    if (!save_decomp) decomptab = Malloc(MPMY_Nproc() * sizeof(Key_t));
#ifdef SHOW_WGTS
    decomp_wgt = Calloc(MPMY_Nproc(), sizeof(float));
#endif
    if (MPMY_Procnum() == 0) {
	float i;
	float n = (float)nin/MPMY_Nproc();
	int j = 0;
	Msgf(("nin is %d\n", nin));
	qsort(key_n, nin, sizeof(Key_t), cmpkey2);
	for (i = n; i < nin; i += n) {
	    assert((int)i < nin);
	    decomptab[j++] = KeyAnd(key_n[(int)i], KeyNot(KeyInt((1<<12)-1)));
	}
	tmp = KeyLshift(KeyInt(1), KEYBITS-1);
	tmp = KeyOr(tmp, KeySub(tmp, KeyInt(1)));
	decomptab[MPMY_Nproc()-1] = tmp;
	Free(key_n);
    }
    Msgf(("Doing decomptab Bcast\n"));
    StopTimer(&DecompTm);
    StartTimer(&DecompWaitTm);
    MPMY_Bcast(decomptab, MPMY_Nproc()*sizeof(Key_t)/sizeof(int), MPMY_INT, 0);
    StopTimer(&DecompWaitTm);
    Msgf(("SetupDecomp done\n"));
}

int
DestDecomp(void *b)
{
    const Key_t key = getkey_s(b);
    int nproc = MPMY_Nproc();
    int procnum = MPMY_Procnum();
    int i = (nproc > 3) ? nproc/2-1 : 1;
    int inc = (i > 2) ? i/2 : 1;

    /* Be fast if nothing is going to move */
    if (KeyLE(key, decomptab[procnum]) && (!procnum || KeyGT(key, decomptab[procnum-1])))
      i = procnum;
    else while (1) {
	if (KeyGT(key, decomptab[i])) i += inc;
	else if (i && KeyLE(key, decomptab[i-1])) i -= inc;
	else break;
	if (inc > 1) inc >>= 1;
    }

    assert (i >= 0 && i < MPMY_Nproc());
    assert (KeyLE(key, decomptab[i]) && (!i || KeyGT(key, decomptab[i-1])));
#ifdef SHOW_WGTS
    decomp_wgt[i] += weight_s(b);
#endif
    return i;
}

void
FinishDecomp(void)
{
#ifdef SHOW_WGTS
    int i;

    MPMY_Combine(decomp_wgt, decomp_wgt, MPMY_Nproc(), MPMY_FLOAT, MPMY_SUM);
    for (i = 0; i < MPMY_Nproc(); i++) {
	singlPrintf("%3d %s %f\n", i, PrintKey(decomptab[i]), decomp_wgt[i]);
    }
    Free(decomp_wgt);
#endif
    if (!save_decomp) {
	Free(decomptab);
	save_decomp = 0;
    }
}
