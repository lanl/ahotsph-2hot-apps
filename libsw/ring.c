#include "ring.h"
#include "bigmalloc.h"
#include "mpmy.h"
#include "gc.h"
#include "singlio.h"
#include <string.h>
#include "Msgs.h"

#define MSGTYPE 142

void
Ring(void *bptr, int bsize, int bnobj,
     void *optr, int osize, int onobj, int oused, 
     void initf(void *, void *), void interactf(void *, void *, int, int))
{
    char *p;
    void *travel_btab;
    void *tmpbuf;
    int travel_size;
    int n, i;
    int from_proc, to_proc;
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int max_nobj = onobj;

    MPMY_Combine(&onobj, &max_nobj, 1, MPMY_INT, MPMY_MAX);

    travel_size = max_nobj * oused;
    travel_btab = Malloc(travel_size);
    tmpbuf = Malloc(travel_size);
    for (i = 0; i < onobj; i++) {
	p = (char *)optr+i*osize;
	initf((char *)travel_btab+i*oused, p);
    }

#if GRAYDECOMP
    /* There should be functions like Gcup() which are periodic */
    to_proc = Gcup(procnum, nproc);
    from_proc = Gcdown(procnum, nproc);
    if (to_proc == -1) to_proc = bin2gray(0);
    if (from_proc == -1) from_proc = bin2gray(nproc-1);
#else
    to_proc = (procnum+1)%nproc;
    from_proc = (procnum+nproc-1)%nproc;
#endif

    /* local part */
    for (p = bptr; p < (char *)bptr + bnobj * bsize; p += bsize) {
	interactf(p, travel_btab, oused, onobj);
    }

    for(n = 1; n < nproc; n++) {
	MPMY_Comm_request req, req2;
	MPMY_Status stat;

	singlPrintf("cycle %d starting\n", n);
	Msgf(("communicate, cycle %d\n", n));
	/* This uses a lot more memory than packets would */
	memcpy(tmpbuf, travel_btab, travel_size);
	MPMY_Irecv(travel_btab, travel_size, from_proc, MSGTYPE, &req);
	MPMY_Isend(tmpbuf, onobj * oused, to_proc, MSGTYPE, &req2);
	MPMY_Wait2(req, &stat, req2, 0);
	onobj = MPMY_Count(&stat)/oused;

	Msgf(("compute, cycle %d\n", n));
	for (p = bptr; p < (char *)bptr + bnobj * bsize; p += bsize) {
	    interactf(p, travel_btab, oused, onobj);
	}
    }
    Free(tmpbuf);
    Free(travel_btab);
}



