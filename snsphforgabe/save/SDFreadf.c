#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdarg.h>
#include "mpmy.h"
#include "SDF.h"
#include "Assert.h"
#include "bigmalloc.h"
#include "Msgs.h"
#include "error.h"
#include "verify.h"
#include "SDFreadf.h"
#include "gc.h"
#include "singlio.h"
#include "timers.h"

#define MAXNAMES 64

extern Timer_t SDFreadTm;

SDF *SDFreadf(char *name, void **btabp, int *gnobjp, int *nobjp, 
	    int stride,
	    /* char *name, offset_t offset, int *confirm */...)
{
    va_list ap;
    int start;
    SDF *sdfp;
    int gnobj, nobj;
    void *btab;
    void *addrs[MAXNAMES];
    char *names[MAXNAMES];
    int strides[MAXNAMES];
    int nobjs[MAXNAMES];
    int starts[MAXNAMES];
    int *confirm;
    int nnames;

    EnableTimer(&SDFreadTm, "SDFread");
    StartTimer(&SDFreadTm);

    VerifySX(sdfp = SDFopen(0, name),SinglShout("%s", SDFerrstring));
    
    if( SDFgetint(sdfp, "npart", &gnobj) ){
	/* Hopefully calling va_start and va_end in here won't disturb */
	/* the real loop over arguments below... */
	va_start(ap, stride);
	names[0] = va_arg(ap, char *);
	gnobj = SDFnrecs(names[0], sdfp);
	va_end(ap);
	if( MPMY_Procnum() == 0 ){
	    SinglShout("%s does not have an \"npart\".\n", name);
	    SinglShout("Guessing npart=%d from SDFnrecs(., %s)\n", 
		       gnobj, names[0]);
	}
    }
    
    NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
    btab = Calloc(nobj, stride);
    Msgf(("Proc %d starting at %d in file, reading %d of %d\n",
	  MPMY_Procnum(), start, nobj, gnobj));

    nnames = 0;
    va_start(ap, stride);
    while(( names[nnames] = va_arg(ap, char *)) != NULL ){
	assert(nnames < MAXNAMES);
	addrs[nnames] = va_arg(ap, int) + (char *)btab;
	confirm = va_arg(ap, int *);
	if( !SDFhasname(names[nnames], sdfp) ){
	    *confirm = 0;
	    Msgf(("SDF file does not have %s\n", names[nnames]));
	    continue;
	}else{
	    *confirm = 1;
	}
	starts[nnames] = start;
	nobjs[nnames] = nobj;
	strides[nnames] = stride;
	nnames++;
    }
    va_end(ap);
    
    VerifyX(0==SDFseekrdvecsarr(sdfp, nnames,
			   names, starts, nobjs, addrs, strides),
	    Shout("%s", SDFerrstring));
    
    *nobjp = nobj;
    *gnobjp = nobj;
    MPMY_Combine(nobjp, gnobjp, 1, MPMY_INT, MPMY_SUM);
    Msgf(("nobj=%d, gnobj=%d\n", nobj, gnobj));

    *btabp = btab;
    StopTimer(&SDFreadTm);
    OutputTimer(&SDFreadTm, singlPrintf); /* global sync and sets timer->max */
    singlPrintf("read speed %.0f kb/s\n", gnobj*nnames*sizeof(float)/(1000.0*SDFreadTm.max));
    DisableTimer(&SDFreadTm);
    return sdfp;
}

