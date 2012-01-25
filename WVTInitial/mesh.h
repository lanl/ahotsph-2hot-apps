#include "key.h"
#include "ndim.h"

typedef struct {
    double pos[NDIM];
    double rho;
    double u;
/*     double nterms; */
    double T;
    double s;
    Key_t key;
} Meshbody;

#define MESHTBODYSZ offsetof(Meshbody, key)

typedef struct {
    double pos[NDIM];
    int ndaughters;
} Meshcell;

typedef struct {
    double pos[NDIM];
    int ndaughters;
} Meshcofmdata;

typedef struct {
    double pos[NDIM];
    double rho;
    double u;
    double T;
    double s;
} Meshoutbody;

typedef struct {
    double pos[NDIM];
    double rho;
    double u;
    int isbody;
} Meshsink;


#define MESHOUTBODYDESC \
"struct {\n\
    double x, y, z;		/* position of body */\n\
    double rho;			/* density */\n\
    double u;                    /* specific internal energy */\n\
    double T;                    /* temperature */\n\
    double s;                    /* specific entropy */\n\
}"


void MeshInit(Meshbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
	      double max[NDIM], int num[NDIM]);
void MeshFixKeys(Meshbody *btab, int nobj, Key_t (*func)(const void *));
double MeshGetCost(const Meshbody *ptr);
Key_t MeshGetKeyFromStruct(const Meshbody *ptr);
void MeshCofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void MeshCellFromCofm(Meshcell *cp, Meshcofmdata *cmp);
void InheritMesh(const Meshsink *from, Meshsink *to, hcell *pp);
void Meshgate(Meshsink *sink, hcell **src_vec, int *result, int n);
char *PrintMeshBody(const Meshbody *p);
char *PrintMeshCell(const Meshcell *p);
