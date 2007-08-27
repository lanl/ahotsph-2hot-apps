/* SetupDecomp() figures out a way to assign every item to a processor */
/* The assignments are available by using DestDecomp() */

#include "Assert.h"
#include "key.h"
#include "decomp.h"
#include "bigmalloc.h"
#include "mpmy.h"
#include "timers.h"
#include "stk.h"
#include "chn.h"
#include "gc.h"
#include "Msgs.h"

#define HASH_BITS 12		/* override default in tree.h */
#include "tree.h"

Timer_t DecompTm;

static Stk ostk;
static tree_t dtree, *t;
static double total_wgt, to_left;	/* prone to roundoff problems */
static Key_t (*getkey_s)(const void *);

/* Make an hcell ptr into a float variable */

#define Wgt(p) (*(float *)(&(p->ptr)))

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
put_ostk(hcellptr p)
{
    if (p->type) {
	StkPushType(&ostk, p->key, Key_t);
	StkPushType(&ostk, Wgt(p), float);
    }
    return 1;
}

void 
set_dest(hcellptr pp, float weight_above)
{
    hcell *p0, *p1;
    Key_t key;

    weight_above /= 2.0;
    pp->type = -1;
    key = KeyLshift(pp->key, 1);
    p0 = Find(t, key);
    p1 = Find(t, KeyOrInt(key,1));

    if (!p0 || !p1) {
	pp->type = ((to_left+weight_above)*MPMY_Nproc())/total_wgt;
	assert (pp->type >= 0 && pp->type < MPMY_Nproc());
    }
    if (p0) set_dest(p0, weight_above+Wgt(pp));
    to_left += Wgt(pp)/2.0;
    if (p1) set_dest(p1, weight_above+Wgt(pp));
    to_left += Wgt(pp)/2.0;
}

void
SetupDecomp(sortresult_t *decompp, 
	    float (*weight)(const void *), Key_t (*getkey)(const void *))
{
    char *b;
    int size = decompp->size;
    int nobj = decompp->nobj;
    char *data = decompp->data;
    hcell *parent = 0, *p;
    Key_t key;
    int level;
    int nin;
    struct {
	Key_t key;
	float n;
    } *key_n, *kn;
    int chan;
    int doc = MPMY_Doc();
    float w;
    
    StartTimer(&DecompTm);
    getkey_s = getkey;

    /* set up 1d tree, with 0 size bodies and cells */
    /* Overload the ptr to accumulate weights */
    t = &dtree;
    SetupTree(t, 1, 0, 0, 0, 0, getkey, 0, 0, 0);
    t->hash_mask = HASH_MASK;
    t->htab = Calloc(HASH_TABLE_SIZE, sizeof(hcellptr));
    ChnInit(&t->hcellchn, sizeof(hcell), sizeof(hcell), Realloc_f);

    Enter(t, KeyInt(1), 0, 0);

    total_wgt = 0.0;
    for(b = data; b < data + nobj * size; b += 10*size) {
	const Key_t key1 = getkey(b);
	for (level = KEYBITS-1; level >=0 ; parent = p, level--) {
	    key = KeyRshift(key1, level);
	    if ((p = Find(t, key)) == NULL) break;
	}
	w = weight(b);
	total_wgt += w;
	assert(level > 0);
	p = Enter(t, key, 0, 1);
	Wgt(p) = w+Wgt(parent)/2.0;
	p = Enter(t, KeyXOR(key,KeyInt(1)), 0, 1);
	Wgt(p) = Wgt(parent)/2.0;
	Wgt(parent) = 0.0;
	parent->type = 0;
    }

    for (chan = 0; chan < doc; chan++) {
	int sendproc = MPMY_Procnum()^(1<<chan);
	int nout;
	if (sendproc < 0 || sendproc >= MPMY_Nproc()) continue;

	StkInitEz(&ostk);
	Tr0(t, KeyInt(1), put_ostk, 0);
	nout = StkSz(&ostk);
	MPMY_Shift(sendproc, &nin, sizeof(int), &nout, sizeof(int), 0);
	key_n = Malloc(nin);
	MPMY_Shift(sendproc, key_n, nin, StkBase(&ostk), nout, 0);
	StkTerminate(&ostk);
	assert (nin % sizeof(*key_n) == 0);
	nin /= sizeof(*key_n);

	for (kn = key_n; kn < key_n + nin; kn++) {
	    key = kn->key;
	    total_wgt += kn->n;
	    if ((p = Find(t, key)) == NULL) {
		p = Enter(t, key, 0, 1);
		Wgt(p) = kn->n;
		while (Find(t, key = KeyRshift(key, 1)) == NULL)
		  Enter(t, key, 0, 0);
	    } else {
		p->type = 1;
		Wgt(p) += kn->n;
	    }
	}
	Free(key_n);
    }
    if (!MPMY_PowOf2()) {
	StkInitEz(&ostk);
	Tr0(t, KeyInt(1), put_ostk, 0);
	nin = StkSz(&ostk);
	MPMY_Bcast(&nin, 1, MPMY_INT, 0);
	if (MPMY_Procnum()) key_n = Malloc(nin);
	else key_n = StkBase(&ostk);
	MPMY_Bcast(key_n, nin, MPMY_CHAR, 0);
	assert (nin % sizeof(*key_n) == 0);
	nin /= sizeof(*key_n);

	total_wgt = 0.0;
	for (kn = key_n; kn < key_n + nin; kn++) {
	    key = kn->key;
	    total_wgt += kn->n;
	    if ((p = Find(t, key)) == NULL) {
		p = Enter(t, key, 0, 1);
		Wgt(p) = kn->n;
		while (Find(t, key = KeyRshift(key, 1)) == NULL)
		  Enter(t, key, 0, 0);
	    } else {
		p->type = 1;
		Wgt(p) = kn->n;
	    }
	}
	StkTerminate(&ostk);
	if (MPMY_Procnum()) Free(key_n);
    }
    to_left = 0.0;
    set_dest(Find(t, KeyInt(1)), 0);
    Msgf(("total_wgt = %f, to_left = %f\n", total_wgt, to_left));

    StopTimer(&DecompTm);
}

int
DestDecomp(void *b)
{
    Key_t key;
    int level;
    hcell *p, *parent;
    const Key_t key1 = getkey_s(b);

    for (level = KEYBITS-1; level >=0 ; parent = p, level--) {
	key = KeyRshift(key1, level);
	p = Find(t, key);
	assert(p);
	if (p->type != -1) {
	    assert (p->type >= 0 && p->type < MPMY_Nproc());
	    return p->type;
	}
    }
    assert(0);
}

void
FinishDecomp(void)
{
    Free(t->htab);
    ChnTerminate(&t->hcellchn);
}
