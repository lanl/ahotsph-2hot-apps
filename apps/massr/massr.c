#include <stddef.h>
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

static float Radius;

void SDFwrite(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);
void memfile_init(int sz);
void memfile_vfprintf(void *junk, const char *fmt, va_list args);

typedef struct {
    float mass;			/* mass of body */
    float mass2;
    float mass4;
    float mass8;
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* position of body */
    int ident;
    int n;
} body;

#define OUTBODYDESC\
"struct {\n\
    float mass;			/* mass of halo within Radius */\n\
    float mass2;		/* mass of halo within 2*Radius */\n\
    float mass4;		/* mass of halo within 4*Radius */\n\
    float mass8;		/* mass of halo within 8*Radius */\n\
    float x, y, z;		/* position of halo */\n\
    float vx, vy, vz;		/* velocity of halo */\n\
    int ident;			/* identity */\n\
    int n;			/* n particles within Radius */\n\
}"

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
    float r2cut = Radius*Radius;
    float r2cut2 = 4.0*Radius*Radius;
    float r2cut4 = 16.0*Radius*Radius;
    float r2cut8 = 64.0*Radius*Radius;

    while (list < last) {
	q = list;
	VxVV(dr, = p->pos, - q->pos);
	r2 = Dotx(dr, dr);
	if (r2 < r2cut8) {
	    p->mass8 += q->mass;
	    if (r2 < r2cut4) {
		p->mass4 += q->mass;
		if (r2 < r2cut2) {
		    p->mass2 += q->mass;
		    if (r2 < r2cut) {
			p->mass += q->mass;
			p->n++;
		    }
		}
	    }
	}
	list = (char *)list + bsize;
    }
}

void
main(int argc, char *argv[])
{
    int massconf;
    int xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int identconf;
    char name[256], outname[256];
    SDF *csdfp, *sdfp;
    body *btab, *ctab;
    int i;
    int nobj, gnobj;
    int cnobj, cgnobj;
    char msg_turn_on[512];
    char msgdir[256];
    char cfile[256];
    int Msg_memfile;
    float tpos, redshift, R0;
    int iter;
    int read_now, level, min_particles;
    float rho_now, cut, h;
 
    MPMY_Init(&argc, &argv);
    if (argc != 3) {
	singlPrintf("usage: %s ctlfile halodata_ctlfile\n", argv[0]);
	exit(1);
    }
    singlPrintf("Welcome to the machine\n");
    
    /* This is necessary since we can't pass "datafile" without a new .ctl */
    csdfp = SDFopen(0, argv[2]);
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
    sdfp = SDFread(csdfp, (void **)&ctab, &cgnobj, &cnobj, sizeof(body),
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
    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
    SDFgetintOrDefault  (sdfp, "read_now",  &read_now, 0);
    SDFgetintOrDefault  (sdfp, "level",  &level, 0);
    SDFgetintOrDefault  (sdfp, "min_particles",  &min_particles, 0);
    SDFgetfloatOrDefault  (sdfp, "rho_now",  &rho_now, (float)0.0);
    SDFgetfloatOrDefault  (sdfp, "cut",  &cut, (float)0.0);
    SDFgetfloatOrDefault  (sdfp, "h",  &h, (float)0.0);

    SDFclose(csdfp);
    SDFclose(sdfp);

    csdfp = SDFopen(0, argv[1]);
    if (csdfp == 0)
      SinglError("Sorry, couldn't SDFopen %s\n%s\n", argv[1], SDFerrstring);
    SDFgetfloatOrDie(csdfp, "radius", &Radius);
    SDFgetstringOrDie(csdfp, "datafile", name, sizeof(name));
    SDFgetstringOrDie(csdfp, "outfile", outname, sizeof(outname));
    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		   "mass", offsetof(body, mass), &massconf,
		   "x", offsetof(body, pos[0]), &xconf,
		   "y", offsetof(body, pos[1]), &yconf,
		   "z", offsetof(body, pos[2]), &zconf,
		   0);
    if( massconf==0 || xconf==0 || yconf==0 || zconf==0 ){
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    SDFclose(sdfp);

    for (i = 0; i < cnobj; i++) {
	ctab[i].mass = (float)0.0;
	ctab[i].n = 0;
    }

    if (redshift != -1.0)
      Radius = Radius/(1.0+redshift);

    Ring(ctab, sizeof(body), cnobj, btab, sizeof(body), nobj, sizeof(body),
	 init, interact);

    SDFwrite(outname, cgnobj, 
	     cnobj, ctab, sizeof(body),
	     OUTBODYDESC,
	     "npart", SDF_INT, cgnobj,
	     "tpos", SDF_FLOAT, tpos,
	     "radius", SDF_FLOAT, Radius,
	     "redshift", SDF_FLOAT, redshift,
	     "R0", SDF_FLOAT, R0,
	     "iter", SDF_INT, iter,
	     "read_now", SDF_INT, read_now,
	     "level", SDF_INT, level,
	     "min_particles", SDF_INT, min_particles,
	     "rho_now", SDF_FLOAT, rho_now,
	     "cut", SDF_FLOAT, cut,
	     "h", SDF_FLOAT, h,
	      NULL);
    singlPrintf("\nOutput to %s done.\n", outname);
    
    exit(0);
}
