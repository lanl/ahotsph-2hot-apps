#include <math.h>
#include "SDF.h"
#include "mpmy.h"
#define NDIM 3
#include "vop.h"
#include "stk.h"
#include "fastflpt.h"
#include "singlio.h"
#include "SDFread.h"
#include "SDFwrite.h"
#include "pqsort.h"
#include "timers.h"

void SDFwrite2(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);

Timer_t WaitTm;

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    int ident;
} body, *bodyptr;

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    int ident;
} outbody, *outbodyptr;

#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;			/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    int ident;			/* idensity */\n\
}"

Key_t OutIdentKey(const outbody *bp)
{
    return KeyLshift(KeyInt(bp->ident), KEYBITS/2);
}

float UnityCost(const void *ptr){
    return 1.0;
}

void
main(int argc, char *argv[])
{
    int massconf;
    int xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int identconf;
    Vxd(float center);
    Vxd(float dr);
    float r2;
    float pick_radius;
    outbodyptr output_btab, p;
    sortresult_t outputsort;
    Stk outstk;
    char name[256], outname[256];
    SDF *csdfp, *sdfp;
    body *btab;
    int nobj, gnobj, i;
    float tpos, redshift, R0;
    int iter;
    int zero_center;
    int do_cube;
    float newt;

    MPMY_Init(&argc, &argv);
    if (argc != 2) {
	singlPrintf("usage: %s ctlfile\n", argv[0]);
	exit(1);
    }
    singlPrintf("Welcome to the machine\n");

    csdfp = SDFopen(0, argv[1]);
    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		   "mass", offsetof(body, mass), &massconf,
		   "x", offsetof(body, pos[0]), &xconf,
		   "y", offsetof(body, pos[1]), &yconf,
		   "z", offsetof(body, pos[2]), &zconf,
		   "vx", offsetof(body, vel[0]), &vxconf,
		   "vy", offsetof(body, vel[1]), &vyconf,
		   "vz", offsetof(body, vel[2]), &vzconf,
		   "ident", offsetof(body, ident), &identconf,
		   NULL);
    if( massconf==0 || xconf==0 || yconf==0 || zconf==0 ){
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    if( identconf==0 || vxconf==0 || vyconf==0 || vzconf==0 ){
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (identconf==0)? "ident" : "",
		   (vxconf==0)? "vx" : "",
		   (vyconf==0)? "vy" : "",
		   (vzconf==0)? "vz" : "");
    }
    SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
    SDFgetfloatOrDefault(sdfp, "redshift",  &redshift, -1.0);
    SDFgetfloatOrDefault(sdfp, "R0",  &R0, 0.0);
    SDFgetfloatOrDefault(sdfp, "Gnewt",  &newt, 1.0);
    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
    
    if(sdfp) SDFclose(sdfp);

    SDFgetfloatOrDie(csdfp, "center_x", &center0);
    SDFgetfloatOrDie(csdfp, "center_y", &center1);
    SDFgetfloatOrDie(csdfp, "center_z", &center2);
    SDFgetfloatOrDie(csdfp, "pick_radius", &pick_radius);
    SDFgetintOrDefault(csdfp, "do_cube", &do_cube, 0);
    SDFgetintOrDefault(csdfp, "zero_center",  &zero_center, 0);

    if (redshift != -1.0)
      pick_radius = pick_radius/(1.0+redshift);

    StkInitEz(&outstk);

    for(i=0; i< nobj; i++) {
	VxVVx(dr, = btab[i].pos, - center);
	if (do_cube) {
	    if (fabs(dr0) < pick_radius && fabs(dr1) < pick_radius 
		&& fabs(dr2) < pick_radius) {
		StkPushData(&outstk, &(btab[i].mass), sizeof(float));
		StkPushData(&outstk, btab[i].pos , NDIM*sizeof(float));
		StkPushData(&outstk, btab[i].vel , NDIM*sizeof(float));
		StkPushData(&outstk, &(btab[i].ident), sizeof(int));
	    }
	} else {
	    r2 = Dotx(dr,dr);
	    if (r2 < pick_radius*pick_radius) {
		StkPushData(&outstk, &(btab[i].mass), sizeof(float));
		StkPushData(&outstk, btab[i].pos , NDIM*sizeof(float));
		StkPushData(&outstk, btab[i].vel , NDIM*sizeof(float));
		StkPushData(&outstk, &(btab[i].ident), sizeof(int));
	    }
	}
    }
    Free(btab);

    SDFgetstring(csdfp, "outfile", outname, sizeof(outname));

    output_btab = StkBase(&outstk);
    nobj = StkSz(&outstk)/sizeof(outbody);
    MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM);

    singlPrintf("Found %d bodies\n", gnobj);

    if (zero_center) {
	singlPrintf("zeroing center about (%f, %f, %f)\n", 
		    center0, center1, center2);
	for (p = output_btab; p < output_btab+nobj; p++) {
	    p->pos[0] -= center0;
	    p->pos[1] -= center1;
	    p->pos[2] -= center2;
	}
    }

#if 0
    singlPrintf("Trying to id sort output\n");
    pqsortsetup_order(&outputsort, output_btab, nobj,
		      sizeof(outbody), 0.1F, 1, Realloc_f);
    output_btab = pqsort(&outputsort,
			 (pq_wgtproto)UnityCost, 
			 (pq_keyproto)OutIdentKey);
    nobj = outputsort.nobj;
#endif

    SDFwrite(outname, gnobj, 
	     nobj, output_btab, sizeof(outbody),
	     OUTBODYDESC,
	     "npart", SDF_INT, gnobj,
	     "center_x", SDF_FLOAT, center0,
	     "center_y", SDF_FLOAT, center1,
	     "center_z", SDF_FLOAT, center2,
	     "pick_radius", SDF_FLOAT, pick_radius,
	     "tpos", SDF_FLOAT, tpos,
	     "redshift", SDF_FLOAT, redshift,
	     "Gnewt", SDF_FLOAT, newt,
	     "R0", SDF_FLOAT, R0,
	     "iter", SDF_INT, iter,
	     "zero_center", SDF_INT, zero_center,
	      NULL);

    singlPrintf("\nOutput done.\n");
    exit(0);
}



