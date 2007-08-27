#include <stdlib.h>
#include <stdio.h>
#include "mpmy.h"
#include "singlio.h"
#include "bigmalloc.h"
#include "timers.h"
#include "Assert.h"
#include "key.h"
#include "gc.h"
#include "pdrmsort.h"

#ifndef STANDALONE				/* poor man's template */
#include "physics.h"
#define Tag_t body
#define weight nterms
#else
typedef struct {
    Key_t key;
    float weight;
    void *ptr;
} Tag_t;
#endif

Timer_t SortTm, Sort1Tm, MergeTm, CommTm, ProbeSortTm, ProbeTm;
Counter_t Nsort;

void
check_sort(char *place, Tag_t *tag, int n)
{
    int i;
    char msg[256];

    for (i = 1; i < n; i++) {
	if (KeyGT(tag[i-1].key, tag[i].key)) {
	    sprintf(msg, "%s\n", PrintKey(tag[i-1].key));
	    sprintf(msg, "%s%s\n", msg, PrintKey(tag[i].key));
	    sprintf(msg, "%s%s\n", msg, PrintKey(tag[i+1].key));
	    Error("%s, Failed at %d!\n%s", place, i, msg);
	}
    }
}

static void
merge_t(Tag_t *a, Tag_t *b, int na, int nb, Tag_t *t)
{
    Tag_t *alast, *blast;
    alast = a + na;
    blast = b + nb;

#if !defined(__PARAGON__)
    while (a < alast && b < blast)
      *t++ = (KeyLT(a->key,b->key)) ? *a++ : *b++;
#else
    /* icc has an apparent compiler optimization bug */
    while (a < alast && b < blast) {
	if (KeyLT(a->key, b->key)) {
	    *t = *a;	    a++;
	} else {
	    *t = *b;
	    b++;
	}
	t++;
    }
#endif
    while (a < alast)
      *t++ = *a++;
    while (b < blast)
      *t++ = *b++;
}

static Tag_t *
rsort_t(Tag_t *tag, int n, int radixBits)
{
    int i, shift, nk;
    unsigned int *keyden;
    _KTYPE *k;
    Tag_t *tag2, *ktmp;
    unsigned int radix = 1 << radixBits;
    unsigned int mask = radix - 1;
    unsigned int stride = sizeof(Tag_t)/sizeof(_KTYPE);

    if (sizeof(Tag_t) % sizeof(_KTYPE)) Error("sizes incompatible\n");
    keyden = Malloc(radix * sizeof(int));
    tag2 = Malloc(n * sizeof(Tag_t));

    for (nk = 0; nk < NK; nk++) {
	for (shift = 0; shift < KEYBITS/NK; shift += radixBits) {
	    for (i = 0; i < radix; i++) 
	      keyden[i] = 0;
	    for (k = &tag[0].key.k[nk]; k < &tag[n].key.k[nk]; k += stride) 
	      ++keyden[*k >> shift & mask];
	    for (i = 1; i < radix; i++) 
	      keyden[i] += keyden[i-1];
	    for (ktmp = tag+n-1; ktmp >= tag; ktmp--)  {
		_KTYPE t = ktmp->key.k[nk];
		tag2[--keyden[t >> shift & mask]] = *ktmp;
	    }
	    ktmp = tag; tag = tag2; tag2 = ktmp;
	}
    }
    Free(tag2);
    Free(keyden);
    return tag;
}


#define NPROBES 128

#ifdef STANDALONE
int
main(int argc, char *argv[])
{
    int i;
    Tag_t *tagtbl;
    int radixBits = 8;
    int n = 1024*1024;
    int nproc, procnum;
    int nk;
    Tag_t check[2], *checktbl = 0;

    MPMY_Init(&argc, &argv);
    nproc = MPMY_Nproc();
    procnum = MPMY_Procnum();
    if (argc >= 2) n = atoi(argv[1])/nproc;
    if (argc == 3) radixBits = atoi(argv[2]);
    EnableTimer(&SortTm, "Sort Total");
    EnableTimer(&Sort1Tm, "Initial Sort");
    EnableTimer(&ProbeSortTm, "Probe Sort");
    EnableTimer(&ProbeTm, "Probe");
    EnableTimer(&CommTm, "Shuffle");
    EnableTimer(&MergeTm, "Merge");
    EnableCounter(&Nsort, "Nsort");
    singlPrintf("Welcome to the Key sort running on %d procs\n", nproc);
    singlPrintf("%d keys (%d per proc)\n", n*nproc, n);
    singlPrintf("radixBits = %d\n", radixBits);

    tagtbl = Malloc(n * sizeof(Tag_t));

    srandom (1+procnum);
    for (i = 0; i < n; i++) {
	for (nk = 0; nk < NK; nk++) {
	    tagtbl[i].key.k[nk] = random();
	    tagtbl[i].weight = 1.0;
	    tagtbl[i].ptr = 0;
	}
    }
    StartTimer(&SortTm);
    pdrmsort(tagtbl, n, radixBits);
    StopTimer(&SortTm);

    /* printf("%d %8ld %8ld %s\n", procnum, tagtbl.key[0].k[1], 
       tagtbl.key[n-1].k[1], PrintKey(key[0])); */
    /* We still have the edge problem for proc 0 and nproc-1 */
    /* printf("%d %d\n", procnum, n); */

    check[0] = tagtbl[0];
    check[1] = tagtbl[n-1];
    if (procnum == 0) checktbl = Malloc(2*nproc*sizeof(Tag_t));
    MPMY_Gather(check, sizeof(Tag_t)*2, MPMY_CHAR, checktbl, 0);

    check_sort("final", tagtbl, n);
    if (procnum == 0 && nproc != 1) check_sort("global", checktbl, nproc*2);
    singlPrintf("Passed.\n");
    OutputTimers(singlPrintf);
    OutputCounters(singlPrintf);
    exit(0);
}
#endif STANDALONE

void
pdrmsort(void *p, int n, int radixBits)
{
    int i, j;
    Tag_t *tagtbl2;
    Tag_t *k, *l;
    int nproc, procnum;
    Tag_t *probes, *splits;
    int *idx, *nintbl;
    int dest;
    int nin, chunk;
    int doc;
    int relative;
    MPMY_Status stat;
    float weight, frac_weight;
    Tag_t *tagtbl = p;

    StartTimer(&SortTm);
    nproc = MPMY_Nproc();
    procnum = MPMY_Procnum();

    StartTimer(&Sort1Tm);
    tagtbl = rsort_t(tagtbl, n, radixBits);
    StopTimer(&Sort1Tm);

    if (nproc == 1) goto done;

    StartTimer(&ProbeTm);
    probes = Malloc(NPROBES*sizeof(Tag_t));
    splits = Malloc(nproc * sizeof(Tag_t));
    nintbl = Malloc(nproc * sizeof(int));

    if (n < NPROBES) Error("n less than NPROBES\n");
    for (weight = 0.0, i = 0; i < n; i++)
      weight += tagtbl[i].weight;
    frac_weight = weight/NPROBES;
    for (weight = tagtbl[0].weight, j = 0, i = 1; i < n; i++) {
	if (weight >= frac_weight) {
	    probes[j] = tagtbl[i];
	    probes[j++].weight = weight;
	    weight = 0.0;
	}
	weight += tagtbl[i].weight;
    }
    if (j == NPROBES-1) {
	probes[j] = tagtbl[n-1];
	probes[j].weight = weight;
    } else Error("Only got %d of %d probes\n", j+1, NPROBES);

    if (procnum == 0) {
	Tag_t *probetbl;
	int nsamp = nproc*NPROBES;

	probetbl = Malloc(nsamp*sizeof(Tag_t));
	MPMY_Gather(probes, sizeof(Tag_t)*NPROBES, MPMY_CHAR, probetbl, 0);
	StartTimer(&ProbeSortTm);
	probetbl = rsort_t(probetbl, nsamp, radixBits);
	StopTimer(&ProbeSortTm);

	for (weight = 0.0, i = 0; i < nsamp; i++)
	  weight += probetbl[i].weight;
	frac_weight = weight/nproc;
	for (weight = probetbl[0].weight, j = 0, i = 1; i < nsamp; i++) {
	    if (weight >= frac_weight) {
		splits[j] = probetbl[i];
		splits[j++].weight = weight;
		weight = 0.0;
	    }
	    weight += probetbl[i].weight;
	}
	if (j == nproc-1) {
	    splits[j] = probetbl[n-1];
	    splits[j].weight = weight;
	} else Error("Only got %d of %d splits\n", j+1, nproc);
	
	Free(probetbl);
    } else MPMY_Gather(probes, sizeof(Tag_t)*NPROBES, MPMY_CHAR, 0, 0);

    MPMY_Bcast(splits, sizeof(Tag_t)*nproc, MPMY_CHAR, 0);

#if 0
    if (procnum == 0) {
	for (i = 0; i < nproc; i++) 
	  printf("%5.2f ", splits[i].weight);
	printf("\n%f\n", frac_weight*nproc);
    }
#endif

    idx = Malloc((nproc+1) * sizeof(int));
    i = idx[0] = 0;
    j = 1;
    while (j < nproc) {
	while (KeyLT(tagtbl[i].key,splits[j-1].key)) i++;
	idx[j] = i;
	j++;
    }
    idx[nproc] = n;
    StopTimer(&ProbeTm);

    StartTimer(&CommTm);
    doc = ilog2(nproc-1) + 1;
    nin = idx[procnum+1]-idx[procnum];
    tagtbl2 = Malloc(nin * sizeof(Tag_t));
    memcpy(tagtbl2, tagtbl+idx[procnum], nin*sizeof(Tag_t));
    nintbl[0] = nin;
    chunk = 1;
    for (relative = 1; relative < 1<<doc; relative++) {
	int nout, nrecv;
	dest = relative ^ procnum;
	if (dest >= nproc) continue;
	nout = idx[dest+1]-idx[dest];
	MPMY_Shift(dest, &nrecv, sizeof(int), &nout, sizeof(int), &stat);
	/* printf("%d-%d send %d recv %d\n", procnum, dest, nout, nrecv); */
	nintbl[chunk++] = nrecv;
	tagtbl2 = Realloc(tagtbl2, (nin+nrecv) * sizeof(Tag_t));
	MPMY_Shift(dest, tagtbl2+nin, nrecv*sizeof(Tag_t), 
		   tagtbl+idx[dest], nout*sizeof(Tag_t), &stat);
	if (MPMY_Count(&stat)/sizeof(Tag_t) != nrecv) 
	  Error("Wrong recv size\n");;
	nin += nrecv;
    }
    Free(tagtbl);
    n = nin;
    tagtbl = Realloc(tagtbl2, n * sizeof(Tag_t));
    AddCounter(&Nsort, n);
    StopTimer(&CommTm);

    StartTimer(&MergeTm);
    tagtbl2 = Malloc(n * sizeof(Tag_t));
    for (j = 1; j < nproc; j <<= 1) {
	k = tagtbl2;
	l = tagtbl;
	for (i = 0; i < nproc; i += 2*j) {
	    merge_t(l, l+nintbl[i], nintbl[i], nintbl[i+j], k);
	    nintbl[i] += nintbl[i+j];
	    k += nintbl[i];
	    l += nintbl[i];
	}
	k = tagtbl; tagtbl = tagtbl2; tagtbl2 = k;
    }

    Free(tagtbl2);
    StopTimer(&MergeTm);
  done:
    StopTimer(&SortTm);
}
