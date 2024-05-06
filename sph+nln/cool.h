/*#include"physics_sph.h"*/

#ifndef NISO
#define NISO 20 	/* number of isotopes tracked */
#endif

#ifndef COOLING
#define COOLING
extern float **tablep; /*array to hold cooling curve table values*/
extern float **ionfracp; /*array to hold ionfraction table values*/
#else
extern float **tablep;
extern float **ionfracp;
extern int do_cooling;
#endif

#ifndef BURNING
#define BURNING
extern int NNW;	/* number of isotopes in network */
extern int **inNW;
extern int nparr[NISO], nnarr[NISO];
#else
extern int NNW;
extern int **inNW;
extern int nparr[NISO];
extern int nnarr[NISO];
extern int do_burning;
#endif

double find_ne(float abundarr[],int nparr[],int nnarr[],double rho,double temp,int Gridpts,int Nel);
double calc_lcool1(float abundarr[],int nparr[],int nnarr[],double rho,double temp,int Gridpts,int Nel,int extrapolate);
double calc_lcool2(double temp,double uint,double dens,int numu,int numD);

double analytic_cool(double temp);

void locate(float xx[], long Nel, float x, long *j);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

void polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy);

void init_CoolTable(int *Gridpts, int *Nel);

int prep_cool_burn(SPHbody *p, float tlo, float tup, int Gridpts, int Nel, int rho_or_rhoest);

float burning(SPHbody *p, float dt, int rank);

float cooling(SPHbody *p, float dt, float frac, int Gridpts, int Nel, int *notprinted);

