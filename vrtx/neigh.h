#ifndef NEIGH_dot_H
#define NEIGH_dot_H

/* structures and prototypes for the neighbor-finding code */

#include "physics_vrtx.h"
#include "timers.h"

#ifdef NDIM
# if NDIM != 3 
# error this code assumes NDIM is 3
# endif
#else
#define NDIM 3
#endif

typedef struct {
    float pos[NDIM];
    float r2;
} ncell;

typedef struct {
    float center[NDIM];
    float size;
} ncofm;

typedef struct {
    body *bp;
} nsink;


#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

void NeighCofmFromDaugh(hcell *hptr, hcell **dlist);
void NeighCellFromCofm(ncell *cellp, ncofm *cofmp);
void NeighInherit(const nsink *from, nsink *to, hcell *pp);
void NeighMACv(nsink *sink, const hcell **srcs, int *results, int nsrc);

extern Counter_t NfindTestsCnt;
extern Counter_t NfindAcceptsCnt;

#ifdef __cplusplus
}
#endif

#endif
