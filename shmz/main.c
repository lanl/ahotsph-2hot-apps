#include <stdlib.h>
#include <stdio.h>
#include "physics.h"
#include "singlio.h"
#include "getparam.h"
#include "mpmy.h"
#include "SDF.h"
#include "SDFread.h"
#include "SDFwrite.h"
#include "bigmalloc.h"
#include "ring.h"
#include "vop.h"
#include "timers.h"
#include "Msgs.h"
#include "bigmalloc.h"
#include "malloc.h"

Timer_t StepTot, StepTotWC, BuildTot;
Timer_t FindForcesTm;
Counter_t MemCnt, HeapCnt_;
static int maxmem(void);
static int maxheap(void);

int
main(int argc, char *argv[])
{
    int gnobj, nobj;
    bodyptr btab;
    char cfile[256];
    char msg_turn_on[256];
    char name[256], outname[256];
    SDF *csdfp;			/* SDF pointer to control file */
    SDF *sdfp;			/* SDF pointer to data file */
    int sconf, xconf, yconf, zconf, idconf;
    int output_nobj;
    int i;
    float lambda;
    outbody *output_btab;
    int do_nsquared;
    float rmin[NDIM], rmax[NDIM];
    sortresult_t sortedbtab, outputsort;
    tree_t thetree;
    inherit_t inherit;
    macv_t mac;
    
    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the scalar Helmholtz solver running on %d procs\n",
		MPMY_Nproc());
    if (argc > 1)
      strncpy(cfile, argv[1], sizeof(cfile));
    else
      Getsparam("control file", cfile);
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
      SinglError("Sorry, couldn't SDFopen %s\n%s\n",
		 cfile, SDFerrstring);
    }

    MsgdirInit("msgs/msg.0");
    SDFgetstringOrDefault(csdfp, "Msg_turn_on", msg_turn_on, 
			  sizeof(msg_turn_on), "");
    Msg_turnon(msg_turn_on);

    singlPrintf("cfile \"%s\" opened\n", cfile);
    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    SDFgetstring(csdfp, "outfile", outname, sizeof(outname));
    SDFgetfloatOrDefault(csdfp, "lambda", &lambda, (float)0.4);
    SDFgetintOrDefault(csdfp, "do_nsquared", &do_nsquared, 0);
    singlPrintf("Reading \"%s\"\n", name);

    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		   "strength", offsetof(body, strength), &sconf,
		   "x", offsetof(body, pos[0]), &xconf,
		   "y", offsetof(body, pos[1]), &yconf,
		   "z", offsetof(body, pos[2]), &zconf,
		   "ident", offsetof(body, ident), &idconf,
		   NULL);

    if (sconf==0 || xconf==0 || yconf==0 || zconf==0) {
      SinglError("Could not find %s %s %s %s in data file!\n",
		 (sconf==0)? "strength" : "",
		 (xconf==0)? "x" : "",
		 (yconf==0)? "y" : "",
		 (zconf==0)? "z" : "");
    }
    FixNterms(btab, nobj);

    EnableTimer(&StepTot, "Step Total");
    EnableWCTimer(&StepTotWC, "Step Tot(WC)");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&FindForcesTm, "Force Eval");
    EnableCounter(&MemCnt, "Mem Used (K)");
    EnableCounter(&HeapCnt_, "Heap Sz (K)");

    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), 0.01, Realloc_f);
    SetupTree(&thetree, NDIM, 
	      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofmdata),
	      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
	      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);

    StartTimer(&StepTot);

    for (i=0; i < nobj; i++) {
	btab[i].phi_r = 0.0;
	btab[i].phi_i = 0.0;
    }

    inherit = (inherit_t)InheritSink;
    mac = (macv_t)OutToIn;

    if (do_nsquared) {
	set_k(lambda);
	Ring(btab, sizeof(body), nobj, btab, sizeof(body), nobj, 
	     TBODYSZ, set_body, do_shmz);
    } else {
	/* FindBbox(btab, nobj, rmin, rmax); */
	VS(rmin, = -1.0);
	VS(rmax, = 1.08);
	FixRsizeExact(rmin, rmax);
	Msgf(("rmin=(%g, %g, %g) rmax=(%g, %g, %g)\n",
	      rmin[0], rmin[1], rmin[2], 
	      rmax[0], rmax[1], rmax[2]));

	set_k(lambda);
	SetupCofm(6, lambda);

	FixKeys(btab, nobj, GetKey);
	singlPrintf("BuildTree\n");

	StartTimer(&BuildTot);
	BuildTree(&thetree, &sortedbtab);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	StopTimer(&BuildTot);
	singlPrintf("BuildTree done %d (%d)\n", maxmem(), maxheap());
	PrintTree(&thetree, PrintBodyContents, PrintCellContents);

	SetTol(gnobj);
	StartTimer(&FindForcesTm);
	Walk(&thetree, &thetree, sizeof(Sink), mac, inherit);
	StopTimer(&FindForcesTm);

	/* This should be the high-water mark for memory use */
	AddCounter(&MemCnt, malloc_used()/1024);
	AddCounter(&HeapCnt_, malloc_heapsz()/1024);

	FreeTree(&thetree);
	singlPrintf("FreeTree done %d (%d)\n", maxmem(), maxheap());
    }

    StopTimer(&StepTot);
    OutputTimer(&StepTot, singlPrintf);

    output_nobj = nobj;
    output_btab = Malloc(output_nobj * sizeof(outbody));
    for (i=0; i < output_nobj; i++){
      output_btab[i].strength = btab[i].strength;
      VV(output_btab[i].pos, = btab[i].pos);
      output_btab[i].phi_r = btab[i].phi_r;
      output_btab[i].phi_i = btab[i].phi_i;
      output_btab[i].ident = btab[i].ident;
    }
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
		      sizeof(outbody), 0.1, 1, Realloc_f);
    output_btab = pqsort(&outputsort,
			 (pq_wgtproto)UnityCost, 
			 (pq_keyproto)OutIdentKey);
    output_nobj = outputsort.nobj;
    SDFwrite(outname, gnobj, 
	     output_nobj, output_btab, sizeof(outbody), OUTBODYDESC,
	     "npart", SDF_INT, gnobj,
	     "lambda", SDF_FLOAT, lambda,
	     NULL);
    Free(output_btab);
    singlPrintf("Bye!\n");
    MPMY_Finalize();
}

static int 
maxheap(void)
{
    int memused = malloc_heapsz()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

static int 
maxmem(void)
{
    int memused = malloc_used()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}
