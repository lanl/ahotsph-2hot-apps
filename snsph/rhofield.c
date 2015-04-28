#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "macr.h"
#include "SDF.h"
#include "mpmy.h"
#define NDIM 3
#include "vop.h"
#include "fastflpt.h"
#include "singlio.h"
#include "SDFreadf.h"
#include "SDFwrite.h"
#include "error.h"
#include "ring.h"
#include "randoms.h"
#include "bigmalloc.h"
#include "pqsort.h"
#include "Assert.h"
#include "Msgs.h"

#define SMOOTH 2

#define Index(i,j,k) ((((i)*nn[1]+(j))*nn[2]+(k)))

void SDFwrite(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);
void memfile_init(int sz);
void memfile_vfprintf(void *junk, const char *fmt, va_list args);

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float temp;
    float u2;
    float ye;
} body;


static double rmin[NDIM], Rmin_pm[NDIM], Rsize_pm, keyfactor_pm, physfactor_pm;

/* return 1 and set xp if inside box */
static int
MeshAssign(const body *p, unsigned int *ip, float *xp){
    double pos[NDIM];

    VV(pos, = p->pos);
    if (pos[0] < Rmin_pm[0]) return 0;
    if (pos[1] < Rmin_pm[1]) return 0;
    if (pos[2] < Rmin_pm[2]) return 0;
    if (pos[0] > Rsize_pm+Rmin_pm[0]) return 0;
    if (pos[1] > Rsize_pm+Rmin_pm[1]) return 0;
    if (pos[2] > Rsize_pm+Rmin_pm[2]) return 0;

    ip[0] = keyfactor_pm*(pos[0]-Rmin_pm[0]);
    ip[1] = keyfactor_pm*(pos[1]-Rmin_pm[1]);
    ip[2] = keyfactor_pm*(pos[2]-Rmin_pm[2]);

    xp[0] = physfactor_pm*ip[0]+Rmin_pm[0];
    xp[1] = physfactor_pm*ip[1]+Rmin_pm[1];
    xp[2] = physfactor_pm*ip[2]+Rmin_pm[2];

    return 1;
}

static void
FixKeyfactor(double *rmin, double size, int Nm)
{
    VV(Rmin_pm, = rmin);
    keyfactor_pm = Nm/size;
    physfactor_pm = size/Nm;
    Rsize_pm = size;
}    


int
main(int argc, char *argv[])
{
    int massconf;
    int xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int tempconf, u2conf, yeconf;
    int i;
    char name[256], basename[256], outbase[256], outname[256], msgfile[256];
    SDF *csdfp, *sdfp;
    body *btab;
    int nobj, gnobj;
    float R0;
    int seed;
    ran_state ranstate;
    float shrink_fac;
    int Nmesh;
    int npts;
    float sysradius;
    float *rhomesh, *vrmesh, *tempmesh, *u2mesh, *yemesh;
    body *p;
    int nn[NDIM];
    float rho_max, rho_min;
    int empty = 0;
    double mtot = 0.0;
    double sum = 0.0;
    FILE *outfp;
    int ii, jj, kk, index; 
    int start_iter, end_iter, stride;
    int iter;

    MPMY_Init(&argc, &argv);
    if (argc != 2) {
	singlPrintf("usage: %s ctlfile\n", argv[0]);
	exit(1);
    }
    singlPrintf("Welcome to the machine\n");
    sprintf(msgfile, "msgs/msg.%d", MPMY_Procnum());
    MsgdirInit(msgfile);
    
    csdfp = SDFopen(0, argv[1]);
    if (csdfp == 0)
      SinglError("Sorry, couldn't SDFopen %s\n%s\n", argv[2], SDFerrstring);

    SDFgetstring(csdfp, "datafile", basename, sizeof(basename));
    SDFgetstring(csdfp, "outfile", outbase, sizeof(outbase));
    SDFgetintOrDie(csdfp, "Nmesh", &Nmesh);
    SDFgetfloatOrDie(csdfp, "R0", &R0);
    SDFgetintOrDefault(csdfp, "seed", &seed, 123);
    SDFgetintOrDefault(csdfp, "start_iter", &start_iter, 0);
    SDFgetintOrDefault(csdfp, "end_iter", &end_iter, 100);
    SDFgetintOrDefault(csdfp, "stride", &stride, 1);
    SDFgetfloatOrDefault(csdfp, "shrink_fac", &shrink_fac, 0);

    SDFclose(csdfp);

    for (iter = start_iter; iter <= end_iter; iter += stride) {

    sprintf(name, "%s.%04d", basename, iter);
    singlPrintf("Reading \"%s\"\n", name);

    sdfp = SDFreadf(name, (void **)&btab, &gnobj, &nobj, sizeof(body),
		   "mass", offsetof(body, mass), &massconf,
		   "x", offsetof(body, pos[0]), &xconf,
		   "y", offsetof(body, pos[1]), &yconf,
		   "z", offsetof(body, pos[2]), &zconf,
		   "vx", offsetof(body, vel[0]), &vxconf,
		   "vy", offsetof(body, vel[1]), &vyconf,
		   "vz", offsetof(body, vel[2]), &vzconf,
		   "temp", offsetof(body, temp), &tempconf,
		   "u2", offsetof(body, u2), &u2conf,
		   "ye", offsetof(body, ye), &yeconf,
		   NULL);
    if( massconf==0 || xconf==0 || yconf==0 || zconf==0 ){
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }

    if (shrink_fac != (float)0.0) {
	int k, nout = 0;
	float p;

	ran_init(seed*(MPMY_Procnum()+1), &ranstate);

	p = 2.0 * 1.0/shrink_fac;
	for(k=0; k<nobj; ) {
	    btab[nout] = btab[k];
	    k++;
	    k += (int)(p * uniform_rand(&ranstate));
	    nout++;
	}
	nobj = nout;
	btab = Realloc(btab, nout*sizeof(body));
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM);
	singlPrintf("sample reduced to %d particles\n", gnobj);
    }


    SDFclose(sdfp);

    sysradius = R0;
    VS(rmin, = -sysradius);
    singlPrintf("R0 is %12.2f\n", R0);

    nn[0] = nn[1] = nn[2] = Nmesh;
    npts = Nmesh*Nmesh*Nmesh;
    FixKeyfactor(rmin, 2.0*sysradius, Nmesh);
    rhomesh = Calloc(npts, sizeof(float));
    vrmesh = Calloc(npts, sizeof(float));
    tempmesh = Calloc(npts, sizeof(float));
    u2mesh = Calloc(npts, sizeof(float));
    yemesh = Calloc(npts, sizeof(float));

    for (p = btab; p < btab+nobj; p++) {
      unsigned int ip[NDIM];
      float xp[NDIM];
      float dr[NDIM], v;
      float h3 = keyfactor_pm*keyfactor_pm*keyfactor_pm;
      float vfac;
      float rr;

      if (MeshAssign(p, ip, xp) 
	  && ip[0] >= SMOOTH
	  && ip[1] >= SMOOTH
	  && ip[2] >= SMOOTH) {
	mtot += p->mass;
	for (ii = -SMOOTH; ii <= SMOOTH; ii++) {
	  for (jj = -SMOOTH; jj <= SMOOTH; jj++) {
	    for (kk = -SMOOTH; kk <= SMOOTH; kk++) {
	      rr = Dot(p->pos, p->pos);
	      if (rr > R0*R0) continue;
	      index = Index(ip[0]+ii, ip[1]+jj, ip[2]+kk);
	      if (index >= npts) {
		continue;
	      }
	      VVV(dr, = p->pos, - xp);
	      v = sqrt(Dot(dr, dr))*keyfactor_pm;
	      
	      if (v < 1.0) {
		vfac = 0.25*(1.0f-1.5f*v*v+0.75f*v*v*v);
		rhomesh[index] += p->mass*vfac;
		tempmesh[index] += p->mass*p->temp*vfac;
		u2mesh[index] += p->mass*p->u2*vfac;
		yemesh[index] += p->mass*p->ye*vfac;
		vrmesh[index] += p->mass * vfac * Dot(p->pos, p->vel)/sqrt(rr);
		sum += p->mass*vfac;
	      } else if (v < 2.0) {
		float tmp = 2.0f-v;
		vfac = 0.25*(0.25f*tmp*tmp*tmp);
		rhomesh[index] += p->mass*vfac;
		tempmesh[index] += p->mass*p->temp*vfac;
		u2mesh[index] += p->mass*p->u2*vfac;
		yemesh[index] += p->mass*p->ye*vfac;
		vrmesh[index] += p->mass * vfac * Dot(p->pos, p->vel)/sqrt(rr);
		sum += p->mass*vfac;
	      }
	    }
	  }
	}
      }
    }
    rho_max = 0.0;
    rho_min = 1e30;
    for (i = 0; i < npts; i++) {
      if (rhomesh[i] == 0.0) {
	empty++;
      } else {
	vrmesh[i] /= rhomesh[i];
	tempmesh[i] /= rhomesh[i];
	u2mesh[i] /= rhomesh[i];
	yemesh[i] /= rhomesh[i];
	rhomesh[i] *= Nmesh*Nmesh*Nmesh/(R0*R0*R0);
	if (rhomesh[i] > rho_max) rho_max = rhomesh[i];
	if (rhomesh[i] != 0.0 && rhomesh[i] < rho_min) rho_min = rhomesh[i];
      }
    }
    singlPrintf("rho_max is %g, rho_min is %g\n", rho_max, rho_min);
    singlPrintf("mtot is %f, sum is %f, full fraction is %.2f\n", 
		mtot, sum, 1.0-(float)empty/npts);


    sprintf(outname, "%s.%04d", outbase, iter);
    Fopen(outfp, outname, "w");

    for (i = 0; i < npts; i++) {
#if 0
      fprintf(outfp, "%g %g %g %g %g\n", rhomesh[i], vrmesh[i], 
	      tempmesh[i], u2mesh[i], yemesh[i]);
#else
      fprintf(outfp, "%g\n", vrmesh[i]);
#endif
    }
    Fclose(outfp);

#if 0
    sprintf(outname, "%s_ye.%04d", outbase, iter);
    Fopen(outfp, outname, "w");
    for (i = 0; i < npts; i++) {
      fprintf(outfp, "%g\n", yemesh[i]);
    }
    Fclose(outfp);
#endif


    Free(btab);
    Free(rhomesh);
    Free(vrmesh);
    Free(tempmesh);
    Free(u2mesh);
    }

    exit(0);
}
