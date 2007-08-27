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

static float Radius, GaussRadius;
static float Physical_L0;

void SDFwrite(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);
void memfile_init(int sz);
void memfile_vfprintf(void *junk, const char *fmt, va_list args);

typedef struct {
    float mass;			/* mass of body */
    float gauss_mass;		/* mass of body */
    float pos[NDIM];		/* position of body */
    int n;
} body;

#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of halo within top-hat Radius */\n\
    float gauss_mass;		/* gauss. window of gauss_radius */\n\
    float x, y, z;		/* position of halo */\n\
    int n;			/* n particles within Radius */\n\
}"

/* Use this to sort by "ident" for output */
float UnityCost(const void *ptr){
    return 1.0;
}

/* Sort by x coord. We just want loadbalance anyway */
Key_t
OutFloatKey(const body *bp)
{
    int i;
    float x = bp->pos[0]/(Physical_L0+1e-6);
    assert (x >= 0.0 && x < 1.0);
    i = (1<<31)*x;
    return KeyLshift(KeyInt(i), KEYBITS/2);

}

void
init(void *o, void *p)
{
    memcpy(o, p, sizeof(body));
}


void 
interact(void *p0, void *list, int bsize, int n)
{
    body *p = p0;
    body *q;
    void *last = (char *)list + bsize*n;
    Vxd(float dr);
    float r2;
    VxdV(const float ppos, = p->pos);
    const float r2cut = Radius*Radius;
    const float r2gcut = GaussRadius*GaussRadius;
    const float compact_r2gcut = (float)16.0*r2gcut;
    const float fac = (float)1.0/((float)2.0*r2gcut);
    assert(compact_r2gcut > r2cut);

    while (list < last) {
	q = list;
	VxVxV(dr, = ppos, - q->pos);
	r2 = Dotx(dr, dr);
	if (r2 < compact_r2gcut) {
	  p->gauss_mass += q->mass * exp(-r2*fac);
	}
	if (r2 < r2cut) {
	  p->mass += q->mass;
	  p->n++;
	}
	list = (char *)list + bsize;
    }
}

void
main(int argc, char *argv[])
{
    int massconf;
    int xconf, yconf, zconf;
    char name[256], outname[256];
    SDF *csdfp, *sdfp;
    body *btab, *ctab;
    int i;
    int nobj, gnobj;
    int cnobj, cgnobj;
    float tpos, redshift, R0, rphys;
    int iter;
    int read_now, level, min_particles;
    float rho_now, cut, h;
    int Nsamples;
    int seed;
    ran_state ranstate;
    float H0;
    float Radius_h100;
    float h_100 = 0.102274;	/* 100 km/s/Mpc in units of 1/Gyr */
    float shrink_fac;
    float total_mass, total_gmass, mean_mass, mean_gmass;
    float dm, sumdm2, sumgdm2, dmm, gdmm;
    float sample_radius;
    float buffer_fac;
 
    MPMY_Init(&argc, &argv);
    if (argc != 2) {
	singlPrintf("usage: %s ctlfile\n", argv[0]);
	exit(1);
    }
    singlPrintf("Welcome to the machine\n");
    
    csdfp = SDFopen(0, argv[1]);
    if (csdfp == 0)
      SinglError("Sorry, couldn't SDFopen %s\n%s\n", argv[2], SDFerrstring);

    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    SDFgetstring(csdfp, "outfile", outname, sizeof(outname));
    SDFgetfloatOrDie(csdfp, "radius_h100", &Radius_h100);
    SDFgetintOrDie(csdfp, "Nsamples", &Nsamples);
    SDFgetintOrDie(csdfp, "seed", &seed);
    SDFgetfloatOrDefault(csdfp, "shrink_fac", &shrink_fac, 0);
    SDFgetfloatOrDefault(csdfp, "buffer_fac", &buffer_fac, 0.05);

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

    SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
    SDFgetfloatOrDefault(sdfp, "redshift",  &redshift, -1.0);
    SDFgetfloatOrDefault(sdfp, "R0",  &R0, 0.0);
    if (R0 == 0.0) 
        SDFgetfloatOrDefault(csdfp, "R0",  &R0, 0.0);
    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
    SDFgetintOrDefault  (sdfp, "read_now",  &read_now, 0);
    SDFgetintOrDefault  (sdfp, "level",  &level, 0);
    SDFgetintOrDefault  (sdfp, "min_particles",  &min_particles, 0);
    SDFgetfloatOrDefault  (sdfp, "rho_now",  &rho_now, (float)0.0);
    SDFgetfloatOrDefault  (sdfp, "cut",  &cut, (float)0.0);
    SDFgetfloatOrDefault  (sdfp, "h",  &h, (float)0.0);
    SDFgetfloatOrDefault(sdfp, "H0", &H0, 0.5*h_100);

    SDFclose(csdfp);
    SDFclose(sdfp);

    Radius = Radius_h100*h_100/H0;
    GaussRadius = Radius*0.47375;

    /* Keep samples a bit away from the edge */
    sample_radius = R0*(1.0-buffer_fac)-Radius;
    singlPrintf("Leaving %.2f buffer from R0 = %.2f\n", R0*buffer_fac, R0);

    rphys = R0;
    if (redshift != -1.0) {
	Radius /= (1.0+redshift);
	GaussRadius /= (1.0+redshift);
	rphys /= (1.0 + redshift);
	sample_radius /= (1.0 + redshift);
    }

    Nsamples /= MPMY_Nproc();
    cnobj = Nsamples;
    MPMY_Combine(&cnobj, &cgnobj, 1, MPMY_INT, MPMY_SUM);
    singlPrintf("Nsamples is %d\n", cgnobj);
    ctab = Malloc(cnobj * sizeof(body));
    for (i = 0; i < cnobj; i++) {
	ctab[i].mass = (float)0.0;
	ctab[i].gauss_mass = (float)0.0;
	ctab[i].n = 0;
	sphere_rand(&ranstate, 3, ctab[i].pos);
	VS(ctab[i].pos, *= sample_radius);
    }

    Physical_L0 = rphys * 2.0;

    singlPrintf("Doing ring\n");
    Ring(ctab, sizeof(body), cnobj, btab, sizeof(body), nobj, sizeof(body),
	 init, interact);

    {
	total_mass = total_gmass = 0.0;
	for (i = 0; i < cnobj; i++) {
	    total_mass += ctab[i].mass;
	    total_gmass += ctab[i].gauss_mass;
	}
	MPMY_Combine(&total_mass, &total_mass, 1, MPMY_FLOAT, MPMY_SUM);
	MPMY_Combine(&total_gmass, &total_gmass, 1, MPMY_FLOAT, MPMY_SUM);
	mean_mass = total_mass / cgnobj;
	mean_gmass = total_gmass / cgnobj;

	sumdm2 = sumgdm2 = 0.0;
	for (i = 0; i < cnobj; i++) {
	    dm = 1.0 - ctab[i].mass/mean_mass;
	    sumdm2 += dm*dm;
	    dm = 1.0 - ctab[i].gauss_mass/mean_gmass;
	    sumgdm2 += dm*dm;
	}
	MPMY_Combine(&sumdm2, &sumdm2, 1, MPMY_FLOAT, MPMY_SUM);
	MPMY_Combine(&sumgdm2, &sumgdm2, 1, MPMY_FLOAT, MPMY_SUM);
	dmm = sqrt(sumdm2/cgnobj);
	gdmm = sqrt(sumgdm2/cgnobj);

	singlPrintf("%8f %8g %8g %8g %8g %8g %8g\n",
		    redshift, mean_mass, mean_gmass, dmm, gdmm, 
		    dmm*(1.0+redshift), gdmm*(1.0+redshift));
    }

#if 0
    SDFwrite(outname, cgnobj, 
	     cnobj, ctab, sizeof(body),
	     OUTBODYDESC,
	     "npart", SDF_INT, cgnobj,
	     "tpos", SDF_FLOAT, tpos,
	     "radius", SDF_FLOAT, Radius,
	     "gauss_radius", SDF_FLOAT, GaussRadius,
	     "sample_radius", SDF_FLOAT, sample_radius,
	     "shrink_fac", SDF_FLOAT, shrink_fac,
	     "nsource", SDF_INT, gnobj,
	     "redshift", SDF_FLOAT, redshift,
	     "R0", SDF_FLOAT, R0,
	     "iter", SDF_INT, iter,
	     "seed", SDF_INT, seed,
	     "nproc", SDF_INT, MPMY_Nproc(),
	     "mean_mass", SDF_FLOAT, mean_mass,
	     "mean_gmass", SDF_FLOAT, mean_gmass,
	     "dmm", SDF_FLOAT, dmm,
	     "gdmm", SDF_FLOAT, gdmm,
	      NULL);
    singlPrintf("\nOutput to %s done.\n", outname);
#endif
    
    exit(0);
}
