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
#include "Msgs.h"

void memfile_init(int sz);
void memfile_vfprintf(void *junk, const char *fmt, va_list args);

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
} halo_s;

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
    Stk outstk;
    char name[256], outname[256], thisname[256];
    SDF *csdfp, *sdfp;
    body *btab;
    halo_s *halos;
    int nhalo;
    int nobj, gnobj, i;
    float tpos, redshift, R0;
    int iter;
    int zero_center;
    float newt;
    int halo;
    int npart, gnpart;
    char msg_turn_on[512];
    char msgdir[256];
    int Msg_memfile;

    MPMY_Init(&argc, &argv);
    if (argc != 2) {
	singlPrintf("usage: %s ctlfile\n", argv[0]);
	exit(1);
    }
    singlPrintf("Welcome to the machine\n");

    csdfp = SDFopen(0, argv[1]);
    if (csdfp == 0)
      SinglError("Sorry, couldn't SDFopen %s\n%s\n", argv[2], SDFerrstring);
    SDFgetintOrDefault(csdfp, "Msg_memfile", &Msg_memfile, 0);
    if (Msg_memfile) {
#ifdef __PARAGON__
	sigio_setup();
#endif
	memfile_init(Msg_memfile);
	Msg_addfile(0, memfile_vfprintf, 0);
	singlPrintf("Putting all Msgs in memfile\n");
    } else {
	sprintf(msgdir, "msgs/msg.%d", MPMY_Procnum());
	MsgdirInit(msgdir);
    }
    SDFgetstringOrDefault(csdfp, "Msg_turn_on", msg_turn_on, 
			  sizeof(msg_turn_on), "");
    Msg_turnon(msg_turn_on);

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

    SDFgetfloatOrDie(csdfp, "pick_radius", &pick_radius);
    SDFgetintOrDefault(csdfp, "zero_center",  &zero_center, 0);
    SDFgetstring(csdfp, "outfile", outname, sizeof(outname));
    SDFgetstringOrDie(csdfp, "halofile", name, sizeof(name));

    {
	int ret;
	int nnames = 4;
	char *names[4];
	int nns[4];
	void *addrs[4];
	int sizes[4];

	SDF_Setiomode(SDF_ASYNC);	/* same data on each node */
	singlPrintf("Reading \"%s\"\n", name);
	sdfp = SDFopen(0, name);
	nhalo = SDFnrecs("x", sdfp);
	halos = Malloc(nhalo * sizeof(halo_s));
	names[0] = "mass";
	names[1] = "x";
	names[2] = "y";
	names[3] = "z";
	addrs[0] = (char *)halos + offsetof(halo_s, mass);
	addrs[1] = (char *)halos + offsetof(halo_s, pos[0]);
	addrs[2] = (char *)halos + offsetof(halo_s, pos[1]);
	addrs[3] = (char *)halos + offsetof(halo_s, pos[2]);
	for(i=0; i<nnames; i++){
	    nns[i] = nhalo;
	    sizes[i] = sizeof(halo_s);
	}
	ret = SDFrdvecsarr(sdfp, nnames, names, nns, addrs, sizes);
	if (ret != 0 ) Error("SDFrdvecs failed: %s\n", SDFerrstring);
	SDFclose(sdfp);
	SDF_Setiomode(SDF_SYNC);
    }

    if (redshift != -1.0)
      pick_radius = pick_radius/(1.0+redshift);
    StkInitEz(&outstk);
    singlPrintf("Finding %d halos\n", nhalo);

    for (halo = 0; halo < nhalo; halo++) {

	VxV(center, = halos[halo].pos);

	for(i=0; i< nobj; i++) {
	    VxVVx(dr, = btab[i].pos, - center);
	    r2 = Dotx(dr,dr);
	    if (r2 < pick_radius*pick_radius) {
		StkPushData(&outstk, &(btab[i].mass), sizeof(float));
		StkPushData(&outstk, btab[i].pos , NDIM*sizeof(float));
		StkPushData(&outstk, btab[i].vel , NDIM*sizeof(float));
		StkPushData(&outstk, &(btab[i].ident), sizeof(int));
	    }
	}

	output_btab = StkBase(&outstk);
	npart = StkSz(&outstk)/sizeof(outbody);
	MPMY_Combine(&npart, &gnpart, 1, MPMY_INT, MPMY_SUM);

	singlPrintf("Found %d bodies\n", gnpart);
	
	if (zero_center) {
	    singlPrintf("zeroing center about (%f, %f, %f)\n", 
			center0, center1, center2);
	    for (p = output_btab; p < output_btab+npart; p++) {
		p->pos[0] -= center0;
		p->pos[1] -= center1;
		p->pos[2] -= center2;
	    }
	}
	
	sprintf(thisname, "%s.%03d", outname, halo);

	SDFwrite(thisname, gnpart, 
		 npart, output_btab, sizeof(outbody),
		 OUTBODYDESC,
		 "npart", SDF_INT, gnpart,
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

	StkClear(&outstk);
    }
    singlPrintf("\nDone.\n");
    exit(0);
}



