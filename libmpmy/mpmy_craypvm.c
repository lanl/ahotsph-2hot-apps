/* 
   Implementation of the mpmy interface in terms of pvm (T3D) primitives.
   This code is based on the vertex implementation.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <pvm3.h>
#include "bigmalloc.h"
#include "chn.h"
#include "protos.h"
#include "Msgs.h"
#include "mpmy.h"
#include "mpmy_abnormal.h"
#include "mpmy_io.h"
#include "verify.h"

struct comm_s{
    void *buf;
    int cnt;
    int src;
    int tag;
};

#define NCOMM 100
static Chn commchn;
static struct comm_s send_comm;

int MPMY_Isend(const void *buf, int cnt, int dest, int tag, MPMY_Comm_request *req){
    int ret;
    ret = pvm_psend(dest, tag, buf, cnt, PVM_BYTE);
    if (ret < 0) Error("pvm_psend returned %d\n", ret);
    *req = &send_comm;
    return MPMY_SUCCESS;
}

int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp){
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

int MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat){
    int Src, inbytes;
    int intid, intag, bufid;
    struct comm_s *comm = req;
    int ret;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	*flag = 1;
	return MPMY_SUCCESS;
    }

    bufid = pvm_nrecv(comm->src, comm->tag);
    if( bufid > 0 ){
	Verify0(pvm_bufinfo(bufid, &inbytes, &intag, &intid));
	Msgf(("pvm_bufinfo(bufid=%d, inbytes=%d, intag=%d, intid=%d)\n", 
	      bufid, inbytes, intag, intid));
	if( inbytes <= comm->cnt ) {
	    ret = pvm_upkbyte(comm->buf, inbytes, 1);
	    if (ret < 0) Error("pvm_upkbyte returned %d\n", ret);
	}		   
	else
	    Error("Incoming message too long, tag=%d, inbytes=%d\n",
		  intag, inbytes);
	ChnFree(&commchn, comm);
	if( stat ){
	    stat->src = pvm_get_PE(intid); /* CRAY-ism */
	    stat->tag = intag;
	    stat->count = inbytes;
	}
	*flag = 1;
    }else{
	Msgf(("Msg_Test(src=%d, tag=%d) not ready\n", comm->src, comm->tag));
	*flag = 0;
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int
MPMY_Wait(MPMY_Comm_request req, MPMY_Status *stat)
{
    int inbytes;
    int intid, intag, bufid;
    struct comm_s *comm = req;
    int ret;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	return MPMY_SUCCESS;
    }

    if (comm->src < -1 || comm->src > _MPMY_nproc_)
      Error("Bad source in MPMY_Wait (%d)\n", comm->src);

    ret = pvm_precv(comm->src, comm->tag, comm->buf, comm->cnt, PVM_BYTE, 
		    &intid, &intag, &inbytes);

    if (ret < 0) Error("pvm_precv(%d, %d, %d) returned %d\n", 
		       comm->src, comm->tag, comm->cnt, ret);
    Msgf(("pvm_precv(inbytes=%d, intag=%d, intid=%d)\n", 
	  inbytes, intag, intid));
    if (inbytes > comm->cnt)
      Error("Incoming message too long, tag=%d, inbytes=%d\n",intag, inbytes);
    ChnFree(&commchn, comm);
    if( stat ){
	stat->src = pvm_get_PE(intid);
	stat->tag = intag;
	stat->count = inbytes;
    }
    return MPMY_SUCCESS;
}

extern void PrintMemfile(void);

int MPMY_Init(int *argcp, char ***argvp){
    int host, proc, dim;

    _MPMY_procnum_ = pvm_get_PE(pvm_mytid());	/* CRAY-ism */
    _MPMY_nproc_ = pvm_gsize(0);
    _MPMY_initialized_ = 1;
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    singlAutoflush(1);
    MPMY_OnAbnormal(PrintMemfile);
    MPMY_OnAbnormal((Abhndlr)Msg_flush);
    return MPMY_SUCCESS;
}

#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
