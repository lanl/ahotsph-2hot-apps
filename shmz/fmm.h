/* External Fortran linkage */

#define Fortran(x) x##_

/* GNU Fortran adds two underscores if there is an underscore in the name */
#ifdef __GNUC__
#define Fortran2(x) x##__
#else
#define Fortran2(x) x##_
#endif

void Fortran(genabm)(void);
void Fortran(gnthph)(int *nth, int *nph, double *csth, double *csph);
void Fortran(gendto)(int *nx, int *mx, complex *d, double *csu, double *csv, 
		     double *r, double *k);
void Fortran(genffsf)(double *wt, int *nu, int *nv, complex *ffsf, complex *c, 
		      int *nd);
void Fortran(cfix2y)(int *mkpure, int *ier, int *n3, 
		     int *ndegl, int *nul, int *nvl, complex *ffsg, 
		     int *ndeg, int *nu, int *nv,  complex *ffsf);
void Fortran(sfzero)(int *nu, int *nv, complex *f, int *nd);
void Fortran(sfadd)(int *nu, int *nv, complex *f, complex *g, complex *h, 
		    int *nd);
void Fortran(sfmult)(int *nu, int *nv, complex *f, complex *g, complex *h, 
		     int *nd);
void Fortran(out2ind)(int *nx, int *mx, int *n, complex *d, double *csu,
		      double *csv, double *r, double *k);
void Fortran(shfqwt)(int *N, double *w);
void Fortran(eval)(int *nu, int *nv, double *w, complex *f, 
		   complex *quad, int *nd);
void Fortran2(st_out_out)(int *L, double *k, int *order, double rotate[3][3],
			  double *dimen, complex *st_do2o, int *pdo2o,
			  complex *st_do2i, int *pdo2i, int stencil[2][2][2],
			  int *is, int isst[3][4]);
void Fortran2(ex_out_out)(int *cubei, int *cubej, int *L, int *ll, 
			  int stencil[2][2][2], int *is, int isst[3][4], 
			  int *pdo2o, complex *st_do2o, complex *do2o);





