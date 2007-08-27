#include "extra.h"
#include "physics.h"
#include "error.h"
#include "bigmalloc.h"
#include "chn.h"

#define N_PER_BLOCK 16
#define NMAX 16
static Chn ids[NMAX];
static int sizes[NMAX];

static void *
chn_extra_alloc(int n)
{
    int i;
    if (n <= 0) Error("Bad size\n");
    if (n % 8) Error("Alignment problems possible\n");

    for (i = 0; i < NMAX; i++) {
	if (sizes[i] == n) return ChnAlloc(ids+i);
	else if (sizes[i] == 0) {
	    sizes[i] = n;
	    ChnInit(ids+i, n, N_PER_BLOCK, Realloc_f);    
	    return ChnAlloc(ids+i);
	}
    }
    Error("Too many different sizes, increase NMAX\n");
}

void
CellExtraFree(void)
{
    int i;

    for (i = 0; i < NMAX; i++) {
	if (sizes[i] != 0) ChnTerminate(ids+i);
    }
}

int 
CellExtraSz(void *pp)
{
  int sz;
  cell *p = pp;
  sz = p->nu * p->nv * sizeof(complex);
  if (sz < 0 || sz > 1024*1024) Error("Bad sz (%d)\n", sz);
  return (sz);
}

void *
CellExtraPtr(void *pp)
{
  cell *p = pp;
  return (p->ffsf);
}

void *
CellExtraAlloc(void *pp)
{
  cell *p = pp;
  void *ptr;
  if (p->nu > 256 || p->nu < 0) Error("Bad nu (%d)\n", p->nu);
  if (p->nv > 256 || p->nv < 0) Error("Bad nv (%d)\n", p->nv);
  ptr = chn_extra_alloc(p->nu * p->nv * sizeof(complex));
  p->ffsf = ptr;
  return(ptr);
}


