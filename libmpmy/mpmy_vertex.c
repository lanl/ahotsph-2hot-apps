/* 
   Implementation of the mpmy interface in terms of nx primitives.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bigmalloc.h"
#include "chn.h"
#include "protos.h"
#include "Msgs.h"
#include "mpmy.h"
#include "mpmy_abnormal.h"
#include "mpmy_io.h"
#include "verify.h"

/* vertex doesn't have asynchronous routines and it doesn't have "handles" */
/* It does have, ntest(), however, and that should be enough */

#define NCOMM 100

struct comm_s{
    void *buf;
    int cnt;
    int src;
    int tag;
};
static Chn commchn;
static struct comm_s send_comm;

int MPMY_Isend(const void *buf, int cnt, int dest, int tag,
	       MPMY_Comm_request *req)
{
    Msgf(("nwrite %d bytes to %d, tag=%d\n", cnt, dest, tag));
    nwrite((void *)buf, cnt, dest, tag, NULL);	/* blocking! */
    *req = &send_comm;
    return MPMY_SUCCESS;
}

int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp)
{
    struct comm_s *comm = ChnAlloc(&commchn);

    if( comm == NULL ){
	return MPMY_FAILED;
    }
    *reqp = comm;
    comm->buf = buf;
    comm->cnt = cnt;
    comm->tag = (tag==MPMY_TAG_ANY)? -1 : tag;
    comm->src = (src==MPMY_SOURCE_ANY)? -1 : src;
    Msgf(("Irecv %d bytes from %d, tag=%d\n", cnt, src, tag));
    return MPMY_SUCCESS;
}

int
MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat)
{
    int inbytes;
    int intid, intag, bufid;
    struct comm_s *comm = req;
    int ret;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	*flag = 1;
	return MPMY_SUCCESS;
    }

    inbytes = ntest(&comm->src, &comm->tag);
    if( inbytes > comm->cnt){
	Error("Message too long, req=%d, src=%d, tag=%d, inbytes=%d, expecting %d\n",
	      req, comm->src, comm->tag, inbytes, comm->cnt);
	return MPMY_FAILED;
    }
    if( inbytes >= 0 ){
	nread(comm->buf, inbytes, &comm->src, &comm->tag, NULL);
	ChnFree(&commchn, comm);
	if( stat ){
	    stat->src = comm->src;
	    stat->tag = comm->tag;
	    stat->count = inbytes;
	}
	*flag = 1;
    }else{
	Msgf(("Msg_Test(src=%d, tag=%d) not ready\n", comm->src, comm->tag));
	*flag = 0;
    }
    return MPMY_SUCCESS;
}

/* nnodes() and mynode() stopped working in vertex 3.2 ???!!! */
int MPMY_Init(int *argcp, char ***argvp){
    int host, proc, dim;

    whoami(&_MPMY_procnum_, &proc, &host, &dim);
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(MPMY_Abannounce);
    _MPMY_nproc_ = 1<<dim;
    _MPMY_initialized_ = 1;
    return MPMY_SUCCESS;
}

#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
