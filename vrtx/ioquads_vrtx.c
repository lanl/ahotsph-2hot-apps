/* Ensure that physics_vrtx.h won't be included again by quads_vrtx.h */
#include "physics_vrtx.h"

#define body iobody
#define bodyptr iobodyptr

typedef struct {
    float pos[3];
    int nterms;
    float psi[3];
    float vel[3];
    float gradvel[3][3];
    float dstr[3];
    float errsum, errsum2;
} body, *bodyptr;

/* And tweak all the external symbols in quads_vrtx.c */
#define InteractCell InteractIOCell
#define CellCnt IOCellCnt

#include "quads_vrtx.c"
