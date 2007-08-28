#include <math.h>
#include "randoms.h"
#include "physics.h"
#include "physics_sph.h"
#include "bigmalloc.h"
#include "vop.h"
#include "singlio.h"
#include "fastflpt.h"
#include "mpmy.h"
#include "gc.h"
#include "Msgs.h"
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"

void *
InitRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
	 SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
	 int set_id, int setpvel, float new_h, float new_u)
{
    int i;
    SDF *sdfp;
    int massconf, xconf, yconf, zconf = 1;
    int vxconf, vyconf, vzconf = 1;
    int identconf;
    body *btab, *p; 
    int nobj, gnobj;
    SPHbody *SPHbtab, *q; 
    float hubble;
    
    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFread(csdfp, (void **)btabp, gnobjp, nobjp, sizeof(body),
		      "mass", offsetof(body, mass), &massconf,
		      "x", offsetof(body, pos[0]), &xconf,
		      "y", offsetof(body, pos[1]), &yconf,
#if NDIM==3
		      "z", offsetof(body, pos[2]), &zconf,
#endif
		      "vx", offsetof(body, vel[0]), &vxconf,
		      "vy", offsetof(body, vel[1]), &vyconf,
#if NDIM==3
		      "vz", offsetof(body, vel[2]), &vzconf,
#endif
		      "ident", offsetof(body, ident), &identconf,
		      NULL);
    nobj = *nobjp;
    gnobj = *gnobjp;
    btab = *(body **)btabp;
    Msgf(("Data read, nobj=%d, gnobj=%d\n", *nobjp, *gnobjp));
    Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
	  MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
    if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    if (vxconf != vyconf || vxconf != vzconf){
	if (setpvel) SinglError("Missing velocity components!\n");
    }
    if (identconf == 0 || set_id){
	SinglWarning("No \"ident\" in file, numbering sequentially\n");
	FixId(btab, nobj, gnobj);
    }

    SPHbtab = Malloc(nobj * sizeof(SPHbody));
    singlPrintf("Setting h to %f\n", new_h);
    singlPrintf("Setting u to %f\n", new_u);
    SDFgetfloatOrDie(sdfp, "hubble",  &hubble);
    for (i = 0; i < nobj; i++) {
	p = btab+i;
	q = SPHbtab+i;
	q->mass = p->mass * 0.1; /* 10 percent baryons */
	p->mass *= 0.9;
	VV(q->pos, = p->pos);
	/* Offset a little so tree build doesn't fail due to identical pos */
	q->pos[0] += 3.9;
	VV(q->vel, = hubble * p->pos);
	q->ident = p->ident + gnobj;
	q->h = new_h;
	q->u = new_u;
    }
    *SPHgnobjp = gnobj;
    *SPHnobjp = nobj;
    *SPHbtabp = SPHbtab;
    return sdfp;
}

void
DarkSPHTestData(void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
		SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, int periodic)
{
    int i;
    ran_state ranstate;
    int seed, cencon;
    int start;
    int gnobj, nobj;
    body *btab, *p;
    int SPHgnobj, SPHnobj;
    SPHbody *SPHbtab, *q;
    float new_u;
    float h, rsq;

    singlPrintf("Generating random dataset\n");
    if (SDFgetint(csdfp, "nobj", &gnobj))
      SinglError("Sorry, you've got to have an \"nobj\"\n");
    SDFgetintOrDefault(csdfp, "seed", &seed, 123);
    SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
    SDFgetfloatOrDefault(csdfp, "new_u", &new_u, 0.0);
    singlPrintf("int seed = %d;\n", seed);
    singlPrintf("int cencon = %d;\n", cencon);

    NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
    btab = (body *) Malloc(nobj*sizeof(SPHbody));
    SPHgnobj = gnobj;
    SPHnobj = nobj;
    SPHbtab = (SPHbody *) Malloc(SPHnobj*sizeof(SPHbody));
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    for (p = &btab[0]; p < &btab[nobj]; p++) {
#ifdef __PARAGON__
	    clear_tregs();	/* avoid system bug */
#endif
	p->mass = 1.0 / gnobj;		 /*   set masses equal */
	if (periodic)
	  rsq = cube_rand(&ranstate, NDIM, p->pos);
	else
	  rsq = sphere_rand(&ranstate, NDIM, p->pos);
	VS(p->vel, = 0.0);
    }
    h = pow((float)8.5/SPHgnobj, .333333);
    for (i = 0; i < nobj; i++) {
	p = btab+i;
	q = SPHbtab+i;
	q->mass = p->mass * 0.1; /* 10 percent baryons */
	p->mass = p->mass * 0.9;
	/* Offset a little so tree build doesn't fail due to identical pos */
	VVS(q->pos, = p->pos, + .001);
	VS(q->vel, = (float)0.0);
	q->ident = p->ident + gnobj;
	q->h = h;
	q->u = new_u;
    }
    singlPrintf("Extracted 10%% baryons from dark matter input\n");
    FixId(btab, nobj, gnobj);
    FixNterms(btab, nobj);
    SPHFixId(SPHbtab, SPHnobj, SPHgnobj);
    SPHFixNterms(SPHbtab, SPHnobj);
    *gnobjp = gnobj;
    *nobjp = nobj;
    *btabp = btab;
    *SPHgnobjp = SPHgnobj;
    *SPHnobjp = SPHnobj;
    *SPHbtabp = SPHbtab;
}


void *
SPHRead(char *name, void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp,
	int set_id, int setpvel, float new_h, float new_u)
{
    SDF *sdfp;
    int i;
    int massconf, xconf, yconf, zconf = 1;
    int vxconf, vyconf, vzconf = 1;
    int hconf, uconf, prconf;
    int identconf;
    int rhoconf, abarconf, tempconf, yeconf;
    int xpconf, xnconf, u2conf, ifleosconf;
    int ynueconf, ynuebconf, ynuxconf;
    int unueconf, unuebconf, unuxconf;
    int ufreezconf;
    int etanueconf, xpfconf, p2conf, p3conf, p4conf;
    int taccconf, iteraccconf;
    SPHbody *btab, *p; 
    SPHoutbody *obtab;
    int nobj, gnobj;
    
    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFreadf(name, (void **)&obtab, gnobjp, nobjp, sizeof(SPHoutbody),
		      "mass", offsetof(SPHoutbody, mass), &massconf,
		      "x", offsetof(SPHoutbody, pos[0]), &xconf,
		      "y", offsetof(SPHoutbody, pos[1]), &yconf,
#if NDIM==3
		      "z", offsetof(SPHoutbody, pos[2]), &zconf,
#endif
		      "vx", offsetof(SPHoutbody, vel[0]), &vxconf,
		      "vy", offsetof(SPHoutbody, vel[1]), &vyconf,
#if NDIM==3
		      "vz", offsetof(SPHoutbody, vel[2]), &vzconf,
#endif
		      "u", offsetof(SPHoutbody, u), &uconf,
		      "h", offsetof(SPHoutbody, h), &hconf,
		      "rho", offsetof(SPHoutbody, rho), &rhoconf,
		      "pr", offsetof(SPHoutbody, pr), &prconf,
		      "ident", offsetof(SPHoutbody, ident), &identconf,
		      "abar", offsetof(SPHoutbody, abar), &abarconf,
		      "temp", offsetof(SPHoutbody, temp), &tempconf,
		      "ye", offsetof(SPHoutbody, ye), &yeconf,
		      "xp", offsetof(SPHoutbody, xp), &xpconf,
		      "xn", offsetof(SPHoutbody, xn), &xnconf,
		      "ifleos", offsetof(SPHoutbody, ifleos), &ifleosconf,
		      "u2", offsetof(SPHoutbody, u2), &u2conf,
		      "ynue", offsetof(SPHoutbody, ynue), &ynueconf,
		      "ynueb", offsetof(SPHoutbody, ynueb), &ynuebconf,
		      "ynux", offsetof(SPHoutbody, ynux), &ynuxconf, 
		      "unue", offsetof(SPHoutbody, unue), &unueconf,
		      "unueb", offsetof(SPHoutbody, unueb), &unuebconf,
		      "unux", offsetof(SPHoutbody, unux), &unuxconf,
		      "ufreez", offsetof(SPHoutbody, ufreez), &ufreezconf,
		      "etanue", offsetof(SPHoutbody, etanue), &etanueconf,
		      "xpf", offsetof(SPHoutbody, xpf), &xpfconf,
		      "p2", offsetof(SPHoutbody, p2), &p2conf,
		      "p3", offsetof(SPHoutbody, p3), &p3conf,
		      "p4", offsetof(SPHoutbody, p4), &p4conf,
		    "taccreted", offsetof(SPHoutbody, taccreted), &taccconf,
		    "iteraccreted", offsetof(SPHoutbody, iteraccreted), &iteraccconf,
		    NULL);
    nobj = *nobjp;
    gnobj = *gnobjp;
    btab = Calloc(nobj, sizeof(SPHbody));
    *btabp = btab;
    for(i=0; i<nobj; i++){
        btab[i].mass = obtab[i].mass;
	VV(btab[i].pos, = obtab[i].pos);
	VV(btab[i].vel, = obtab[i].vel);
	btab[i].u = obtab[i].u;
	btab[i].h = obtab[i].h;
	btab[i].rho = obtab[i].rho;
	btab[i].pr = obtab[i].pr;
	btab[i].ident = obtab[i].ident;
	btab[i].abar = obtab[i].abar;
	btab[i].temp = obtab[i].temp;
	btab[i].ye = obtab[i].ye;
	btab[i].xp = obtab[i].xp;
	btab[i].xn = obtab[i].xn;
	btab[i].ifleos = obtab[i].ifleos;
	btab[i].u2 = obtab[i].u2;
	btab[i].ynue = obtab[i].ynue;
	btab[i].ynueb = obtab[i].ynueb;
	btab[i].ynux = obtab[i].ynux;
	btab[i].unue = obtab[i].unue;
	btab[i].unueb = obtab[i].unueb;
	btab[i].unux = obtab[i].unux;
	btab[i].ufreez = obtab[i].ufreez;
	btab[i].etanue = obtab[i].etanue;
	btab[i].xpf = obtab[i].xpf;
	btab[i].p2 = obtab[i].p2;
	btab[i].p3 = obtab[i].p3;
	btab[i].p4 = obtab[i].p4;
    }
    Free(obtab);

    Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", *nobjp, *gnobjp));
    Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
	  MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
    if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    if (vxconf != vyconf || vxconf != vzconf){
	if (setpvel) SinglError("Missing velocity components!\n");
    }
    /* Initialize missing members to 0.0 (11/8/2004 --gmr) */
    if (ynueconf == 0 || ynuebconf == 0 || ynuxconf == 0)
/*       SinglError("Missing ynus in data file\n"); */
	for (p = btab; p < btab+nobj; ++p) { p->ynue = 0.0;p->ynueb = 0.0;p->ynux = 0.0;}
    if (unueconf == 0 || unuebconf == 0 || unuxconf == 0)
/*       SinglError("Missing unus in data file\n"); */
	for (p = btab; p < btab+nobj; ++p) { p->unue = 0.0;p->unueb = 0.0;p->unux = 0.0;}
    if (identconf == 0 || set_id){
	SinglWarning("No \"ident\" in file, numbering sequentially\n");
	SPHFixId(btab, nobj, gnobj);
    }
    if (new_h != (float)0.0) {
	singlPrintf("Setting h to %f\n", new_h);
	for (p = btab; p < btab+nobj; p++) p->h = new_h;
    } else if (hconf == 0) {
	SinglError("No h in data file\n");
    }
    if (prconf == 0) {
      /* pr needed for ghosts */
      SinglError("No pr in data file\n");
    }
    if (new_u != (float)0.0) {
	singlPrintf("Setting u to %f\n", new_u);
	for (p = btab; p < btab+nobj; p++)  p->u = new_u;
    } else if (uconf == 0) {
	SinglError("No u in data file\n");
    }
    if (!(abarconf && tempconf && yeconf && xpconf && xnconf 
	  && u2conf &&ifleosconf && ufreezconf))
      SinglError("Missing required field in initial data\n");
    if (etanueconf == 0) {
      for (p = btab; p < btab+nobj; p++) p->etanue = 0.0;
    } 
    if (xpfconf == 0) {
      for (p = btab; p < btab+nobj; p++) p->xpf = p->rho * 1.204e-9 * p->ye;
    }
    if (p2conf == 0) {
      for (p = btab; p < btab+nobj; p++) p->p2 = 0.155;
    }
    if (p3conf == 0) {
      for (p = btab; p < btab+nobj; p++) p->p3 = -15.0;
    }
    if (p4conf == 0) {
      for (p = btab; p < btab+nobj; p++) p->p4 = -10.0;
    }
    if (taccconf == 0) {
	Msgf(("No \"taccreted\" in file, setting to zero\n"));
	for (p = btab; p < btab+nobj; p++) p->taccreted = 0.0;
    }
    if (iteraccconf == 0) {
	Msgf(("No \"iteraccreted\" in file, setting to zero\n"));
	for (p = btab; p < btab+nobj; p++) p->iteraccreted = 0;
    }
    return sdfp;
}

void *
DarkRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp,
	int set_id, int setpvel)
{
    SDF *sdfp;
    int massconf, xconf, yconf, zconf = 1;
    int vxconf, vyconf, vzconf = 1;
    int identconf;
    body *btab;
    int nobj, gnobj;
    
    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFreadf(name, (void **)btabp, gnobjp, nobjp, sizeof(body),
		      "mass", offsetof(body, mass), &massconf,
		      "x", offsetof(body, pos[0]), &xconf,
		      "y", offsetof(body, pos[1]), &yconf,
#if NDIM==3
		      "z", offsetof(body, pos[2]), &zconf,
#endif
		      "vx", offsetof(body, vel[0]), &vxconf,
		      "vy", offsetof(body, vel[1]), &vyconf,
#if NDIM==3
		      "vz", offsetof(body, vel[2]), &vzconf,
#endif
		      "ident", offsetof(body, ident), &identconf,
		      NULL);
    nobj = *nobjp;
    gnobj = *gnobjp;
    btab = *btabp;
    Msgf(("Data read, nobj=%d, gnobj=%d\n", *nobjp, *gnobjp));
    Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
	  MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
    if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    if (vxconf != vyconf || vxconf != vzconf){
	if (setpvel) SinglError("Missing velocity components!\n");
    }
    if (identconf == 0 || set_id){
	SinglWarning("No \"ident\" in file, numbering sequentially\n");
	FixId(btab, nobj, gnobj);
    }
    return sdfp;
}

void
SPHTestData(void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, int periodic)
{
    ran_state ranstate;
    int seed, cencon;
    int start;
    int gnobj, nobj;
    SPHbody *btab, *p;
    float new_u;
    float h, rsq;

    singlPrintf("Generating random dataset\n");
    if (SDFgetint(csdfp, "nobj", gnobjp))
      SinglError("Sorry, you've got to have an \"nobj\"\n");
    gnobj = *gnobjp;
    SDFgetintOrDefault(csdfp, "seed", &seed, 123);
    SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
    SDFgetfloatOrDefault(csdfp, "new_u", &new_u, 0.0);
    singlPrintf("int seed = %d;\n", seed);
    singlPrintf("int cencon = %d;\n", cencon);

    NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
    btab = (SPHbody *) Malloc(nobj*sizeof(SPHbody));
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    h = pow((float)8.5/gnobj, .333333);
    for (p = &btab[0]; p < &btab[nobj]; p++) {
#ifdef __PARAGON__
	clear_tregs();	/* avoid system bug */
#endif
	p->mass = 1.0 / gnobj;		 /*   set masses equal */
	if (periodic) {
	    rsq = cube_rand(&ranstate, NDIM, p->pos);
	} else {
	    rsq = sphere_rand(&ranstate, NDIM, p->pos);
	}
	if (cencon == 1) {
	    rsq = -1.0/sqrt(rsq);
	    VV(p->vel, = rsq * p->pos);
	} else {
	    VS(p->vel, = 0.0);
	}
	p->h = h;
	p->u = new_u;
    }
    SPHFixId(btab, nobj, gnobj);
    SPHFixNterms(btab, nobj);
    *nobjp = nobj;
    *btabp = btab;
}
