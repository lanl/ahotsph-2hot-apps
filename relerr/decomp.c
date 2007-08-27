/* SetupDecomp() figures out a way to assign every item to a processor */
/* The assignments are available by using DestDecomp() */

#include "Assert.h"
#include "key.h"
#include "decomp.h"
#include "bigmalloc.h"
#include "mpmy.h"
#include "timers.h"

#define DECOMP_PBITS 10

Timer_t DecompTm;
Timer_t DecompWaitTm;
Timer_t DecompCommTm;

static float *decomptab;
static int decomp_bits;
static int decomp_mask;
static int decomp_size;
static Key_t (*getkey_d)(const void *);

void
SetupDecomp(sortresult_t *decompp, 
	    float (*weight)(const void *), Key_t (*getkey)(const void *))
{
    char *p;
    int i;
    float total_weight = (float)0.0;
    float fac;
    float w;
    int bin;
    int nproc = MPMY_Nproc();
    int size = decompp->size;
    int nobj = decompp->nobj;
    char *data = decompp->data;

    StartTimer(&DecompTm);

    getkey_d = getkey;
    decomp_bits = DECOMP_PBITS+ilog2(nproc);
    decomp_size = 1 << decomp_bits;
    decomp_mask = decomp_size-1;
    decomptab = Calloc(decomp_size, sizeof(float));

    for(p = data; p < data + nobj * size; p += size){
	bin = KeyAndInt(KeyRshift(getkey(p),(KEYBITS-1)-decomp_bits),
			decomp_mask);
	assert(bin >= 0 && bin < decomp_size);
	w = weight(p);
	assert(w > (float)0.0);
	decomptab[bin] += w;
    }

    MPMY_Combine(decomptab, decomptab, decomp_size, MPMY_FLOAT, MPMY_SUM);
    
    for (i = 0; i < decomp_size; i++) {
	assert(decomptab[i] >= (float)0.0);
	total_weight += decomptab[i];
	decomptab[i] = total_weight;
    }
    fac = total_weight/nproc;
    for (i = 0; i < decomp_size; i++) {
	decomptab[i] = (int)(decomptab[i]/fac);
	if (decomptab[i] > nproc-1)
	  decomptab[i] = nproc-1;
    }
    StopTimer(&DecompTm);
}

int
DestDecomp(void *p)
{
    int bin = KeyAndInt(KeyRshift(getkey_d(p),(KEYBITS-1)
				  -decomp_bits),decomp_mask);
    assert(bin >= 0 && bin < decomp_size);
    return decomptab[bin];
}

void
FinishDecomp(void)
{
    Free(decomptab);
}
