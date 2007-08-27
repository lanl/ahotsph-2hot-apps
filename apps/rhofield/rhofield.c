#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include "SDF.h"
#include "mpmy.h"
#define NDIM 3
#include "vop.h"
#include "fastflpt.h"
#include "singlio.h"
#include "SDFread.h"
#include "SDFwrite.h"
#include "error.h"
#include "ring.h"
#include "randoms.h"
#include "bigmalloc.h"
#include "pqsort.h"
#include "Assert.h"
#include "Msgs.h"

#define Index(i,j,k) ((((i)*nn[1]+(j))*nn[2]+(k)))

void SDFwrite(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);
void memfile_init(int sz);
void memfile_vfprintf(void *junk, const char *fmt, va_list args);

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
} body;

#define OUTBODYDESC \
"struct {\n\
    unsigned char log_rho;\n\
}"

static double rmin[NDIM], Rmin_pm[NDIM], Rsize_pm, keyfactor_pm;

/* return 1 and set xp if inside box */
static int
MeshAssign(const body *p, unsigned int *xp){
    double pos[NDIM];

    VV(pos, = p->pos);
    if (pos[0] < Rmin_pm[0]) return 0;
    if (pos[1] < Rmin_pm[1]) return 0;
    if (pos[2] < Rmin_pm[2]) return 0;
    if (pos[0] > Rsize_pm+Rmin_pm[0]) return 0;
    if (pos[1] > Rsize_pm+Rmin_pm[1]) return 0;
    if (pos[2] > Rsize_pm+Rmin_pm[2]) return 0;
    xp[0] = keyfactor_pm*(pos[0]-Rmin_pm[0]);
    xp[1] = keyfactor_pm*(pos[1]-Rmin_pm[1]);
    xp[2] = keyfactor_pm*(pos[2]-Rmin_pm[2]);
    return 1;
}

static void
FixKeyfactor(double *rmin, double size, int Nm)
{
    VV(Rmin_pm, = rmin);
    keyfactor_pm = Nm/size;
    Rsize_pm = size;
}    


void
main(int argc, char *argv[])
{
    int massconf;
    int xconf, yconf, zconf;
    int i;
    char name[256], outname[256], hdrname[256], msgfile[256];
    SDF *csdfp, *sdfp;
    body *btab;
    int nobj, gnobj;
    float R0;
    int seed;
    ran_state ranstate;
    float redshift;
    float shrink_fac;
    int Nmesh;
    int npts;
    float sysradius;
    float *mesh;
    unsigned char *cmesh;
    body *p;
    int nn[NDIM];
    float rho_fac, rho_max, rho_min;
    int empty = 0;
    double mtot = 0.0;
    float offset, scale;
    unsigned char pix_min, pix_max;
 
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

    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    SDFgetstring(csdfp, "outfile", outname, sizeof(outname));
    SDFgetintOrDie(csdfp, "Nmesh", &Nmesh);
    SDFgetintOrDie(csdfp, "seed", &seed);
    SDFgetfloatOrDefault(csdfp, "shrink_fac", &shrink_fac, 0);
    SDFgetfloatOrDefault(csdfp, "offset", &offset, 0.0);
    SDFgetfloatOrDefault(csdfp, "scale", &scale, 60.0);

    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		   "mass", offsetof(body, mass), &massconf,
		   "x", offsetof(body, pos[0]), &xconf,
		   "y", offsetof(body, pos[1]), &yconf,
		   "z", offsetof(body, pos[2]), &zconf,
		   NULL);
    if( massconf==0 || xconf==0 || yconf==0 || zconf==0 ){
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }

    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    if (shrink_fac != (float)0.0) {
	int k, nout = 0;
	float p;

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


    SDFgetfloatOrDefault(sdfp, "redshift",  &redshift, 0.0);
    SDFgetfloatOrDefault(sdfp, "R0",  &R0, 0.0);
    if (R0 == 0.0) 
        SDFgetfloatOrDefault(csdfp, "R0",  &R0, 0.0);
    SDFclose(csdfp);
    SDFclose(sdfp);

    sysradius = R0 * (1 + 1e-6) / (1.0 + redshift);
    VS(rmin, = -sysradius);
    singlPrintf("R0 is %12.2f, redshift is %4.2f\n", R0, redshift);

    nn[0] = nn[1] = nn[2] = Nmesh;
    npts = Nmesh*Nmesh*Nmesh;
    FixKeyfactor(rmin, 2.0*sysradius, Nmesh);
    mesh = Calloc(npts, sizeof(float));
    cmesh = Calloc(npts, sizeof(unsigned char));

    for (p = btab; p < btab+nobj; p++) {
      unsigned int ii;
      unsigned int xp[NDIM];
      if (MeshAssign(p, xp)) {
	ii = Index(xp[0], xp[1], xp[2]);
	if (ii >= npts) {
	    SeriousWarning("index too large\n");
	    continue;
	}
	mesh[ii] += p->mass;
	mtot += p->mass;
      }
    }
    sysradius /= 1000.0;  /* 1e10 Msol per Mpc */
    rho_fac = (Nmesh*Nmesh*Nmesh)/(sysradius*sysradius*sysradius);
    rho_max = 0.0;
    rho_min = 1e30;
    for (i = 0; i < npts; i++) {
      mesh[i] *= rho_fac;
      if (mesh[i] > rho_max) rho_max = mesh[i];
      if (mesh[i] != 0.0 && mesh[i] < rho_min) rho_min = mesh[i];
      if (mesh[i] == 0.0) empty++;
    }
    singlPrintf("rho_max is %g, rho_min is %g\n", rho_max, rho_min);
    singlPrintf("mtot is %f, full fraction is %.2f\n", 
		mtot, 1.0-(float)empty/npts);
    pix_min = 255;
    pix_max = 0;
    for (i = 0; i < npts; i++) {
	if (mesh[i] == 0.0) cmesh[i] = 0;
	else {
	    float val = 1.0+(log10(mesh[i])+offset)*scale;
	    if (val > pix_max) pix_max = val;
	    if (val < pix_min) pix_min = val;
	    if (val > 255.0) val = 255.0;
	    cmesh[i] = val;
	}
    }
    singlPrintf("pix_max is %d, pix_min is %d\n", (int)pix_max, (int)pix_min);

    sprintf(hdrname, "%s.hdr", outname);
    SDFwritehdr(hdrname, OUTBODYDESC, 
	"rho_min", SDF_FLOAT, rho_min, 
	"rho_max", SDF_FLOAT, rho_max, 
	NULL);

    SDFwrite(outname, npts, 
	     npts, cmesh, sizeof(unsigned char),
	     OUTBODYDESC, NULL);
    singlPrintf("\nOutput to %s done.\n", outname);
    
    exit(0);
}
