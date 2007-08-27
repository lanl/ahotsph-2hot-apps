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
#define HASH_BITS 12		/* override default in tree.h */
#include "tree.h"

#define LEAF_MAX 200

Timer_t DecompTm;

static Stk ostk;
static tree_t dtree, *t;
static int total_wgt;
static int to_left;
static int leaf_min;
static Key_t (*getkey_s)(const void *);

static int 
put_ostk(hcellptr p)
{
    StkPushType(&ostk, p->key, Key_t);
    StkPushType(&ostk, p->type, int);
    return (p->type >= leaf_min);
}

static int 
fixup_n(hcellptr p)
{
    Key_t key = KeyLshift(p->key, 1);
    hcellptr p0, p1;
    int n = p->type;

    p0 = Find(t, key);
    p1 = Find(t, KeyOrInt(key,1));
    if (p0 && p1) {
	int diff = p0->type + p1->type - n;
	if (diff) {
	    p0->type += diff/2;
	    p1->type += diff - diff/2;
	}
    } else if (p0) {
	if (p0->type != n) p0->type = n;
    } else if (p1) {
	if (p1->type != n) p1->type = n;
    } else {
	return 0;
    }
    return 1;
}

static int 
set_dest(hcellptr p)
{
    Key_t key = KeyLshift(p->key, 1);
    hcellptr p0, p1;
    int n = p->type;

    p0 = Find(t, key);
    p1 = Find(t, KeyOrInt(key,1));
    if (n < LEAF_MAX || !p0 || !p1) {
	to_left += n;
	p->type = (to_left*MPMY_Nproc())/total_wgt;
	assert (p->type >= 0 && p->type < MPMY_Nproc());
	return 0;
    } else {
	p->type = -1;
	return 1;
    }
}

static int 
print_dest(hcellptr p)
{
    printf("%d %s %d\n", MPMY_Procnum(), PrintKey(p->key), p->type);
    if (p->type != -1) return 0;
    else return 1;
}

void
SetupDecomp(sortresult_t *decompp, 
	    float (*weight)(const void *), Key_t (*getkey)(const void *))
{
    char *b;
    int size = decompp->size;
    int nobj = decompp->nobj;
    char *data = decompp->data;
    hcell *parent, *p;
    Key_t key;
    int level;
    int nin;
    int sum;
    struct {
	Key_t key;
	int n;
    } *key_n, *kn;
    
    StartTimer(&DecompTm);
    getkey_s = getkey;
    leaf_min = nobj/128;
    if (leaf_min < 2) leaf_min = 2;
    leaf_min &= ~1;
    singlPrintf("leaf min is %d\n", leaf_min);

    /* set up 1d tree, with 0 size bodies and cells */
    /* Use type to count number of daughters */
    t = &dtree;
    SetupTree(t, 1, 0, 0, 0, 0, getkey, 0, 0, 0);
    t->hash_mask = HASH_MASK;
    t->htab = Calloc(HASH_TABLE_SIZE, sizeof(hcellptr));
    ChnInit(&t->hcellchn, sizeof(hcell), sizeof(hcell), Realloc_f);

    Enter(t, KeyInt(1), 0, 0);

    for(b = data; b < data + nobj * size; b += size) {
	for (level = KEYBITS-1; level >=0 ; parent = p, level--) {
	    key = getkey(b);
	    key = KeyRshift(key, level);
	    if ((p = Find(t, key)) == NULL) 
	      break;
	    else
	      p->type++;
	}
	assert(level > 0);
	if (parent->type == leaf_min) { /* give half to each daughter */
	    Enter(t, key, 0, leaf_min/2);
	    Enter(t, KeyXOR(key,KeyInt(1)), 0, leaf_min/2);
	} 
    }
    StkInitEz(&ostk);
    Tr0(t, KeyInt(1), put_ostk, 0);
    printf("Stk sz is %d\n", StkSz(&ostk));
    nin = MPMY_NGather(StkBase(&ostk), StkSz(&ostk), MPMY_CHAR, &key_n, 0);

    MPMY_Bcast(&nin, 1, MPMY_INT, 0);
    if (MPMY_Procnum()) key_n = Malloc(nin);
    MPMY_Bcast(key_n, nin, MPMY_CHAR, 0);

    assert (nin % sizeof(*key_n) == 0);
    nin /= sizeof(*key_n);
    StkTerminate(&ostk);
    Free(t->htab);
    ChnTerminate(&t->hcellchn);

    SetupTree(t, 1, 0, 0, 0, 0, getkey, 0, 0, 0);
    t->hash_mask = HASH_MASK;
    t->htab = Calloc(HASH_TABLE_SIZE, sizeof(hcellptr));
    ChnInit(&t->hcellchn, sizeof(hcell), sizeof(hcell), Realloc_f);

    for (kn = key_n; kn < key_n + nin; kn++) {
	key = kn->key;
	if ((p = Find(t, key))) {
	    p->type += kn->n;
	} else {
	    Enter(t, key, 0, kn->n);
	}
    }
    singlPrintf("total weight is %d\n", Find(t, KeyInt(1))->type);
    Free(key_n);
    Tr0(t, KeyInt(1), fixup_n, 0);
    total_wgt = Find(t, KeyInt(1))->type + 1;
    to_left = 0;
    Tr0(t, KeyInt(1), set_dest, 0);
    StopTimer(&DecompTm);
}

int
DestDecomp(void *b)
{
    Key_t key;
    int level;
    hcell *p, *parent;

    for (level = KEYBITS-1; level >=0 ; parent = p, level--) {
	key = getkey_s(b);
	key = KeyRshift(key, level);
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
