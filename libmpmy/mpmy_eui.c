#include "chn.h"
#include "bigmalloc.h"
#include "mpmy.h"
#include "Assert.h"
#include "timers.h"
#include "mpmy_abnormal.h"
#include "Msgs.h"
#include "error.h"

#define EUI_SUCCEED 0
#define EUI_ACTIVE -1
#define EUI_ERROR  -2

/* This is some crap that mpcc creates */
int mpmondata=0; int mp_linked_euilib=0;

/* #include "mpceui.h" */
/* #include "mpctof.c" */
#include "mpproto.h"

struct comm_s {
    int hndl;
    int inout;
    int src;
    int tag;
};

#define IN 1
#define OUT 2

#define NCOMM 200
static Chn commchn;

int MPMY_Isend(const void *buf, int cnt, int dest, int tag,
	       MPMY_Comm_request *reqp) {
    struct comm_s *comm = ChnAlloc(&commchn);

    Msgf(("Isend: buf=%p, dest=%d, tag=%d\n",
	  buf, dest, tag));
    if (mpc_send(buf, cnt, dest, tag, &comm->hndl) < EUI_SUCCEED)
	Error("MPMY_Isend mpc_send failed\n");
    Msgf(("Isend: hndl=%d\n", comm->hndl));
    comm->inout = OUT;
    *reqp = (MPMY_Comm_request)comm;
    return MPMY_SUCCESS;
}

int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp) {
    struct comm_s *comm = ChnAlloc(&commchn);

    if (tag == MPMY_TAG_ANY)
	comm->tag = DONTCARE;
    else
	comm->tag = tag;

    if (src == MPMY_SOURCE_ANY)
	comm->src = DONTCARE;
    else
	comm->src = src;

    Msgf(("Irecv: buf=%p, src=%d, tag=%d\n",
	  buf, src, tag));
    if (mpc_recv(buf, cnt, &comm->src, &comm->tag, &comm->hndl) < EUI_SUCCEED)
	Error("MPMY_Irecv mpc_recv failed\n");
    comm->inout = IN;
    Msgf(("Irecv: hndl=%d\n", comm->hndl));
    *reqp = (MPMY_Comm_request)comm;
    return MPMY_SUCCESS;
}

int MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat) {
    struct comm_s *comm = req;

    Msgf(("MPMY_Test hndl=%d\n", comm->hndl));
    if ((*flag = mpc_status(comm->hndl)) == EUI_ERROR)
	Error("MPMY_Test mpc_status failed\n");
    Msgf(("Tested (%s), %d\n", 
	  (comm->inout==IN)?"in":"out", *flag));
    if (*flag >= 0) {
	if(comm->inout == IN && stat) {
	    stat->count = *flag;
	    stat->src = comm->src;
	    stat->tag = comm->tag;
	}
	/* Convert EUI zero-length count to MPMY "true" flag.  */
	if (*flag == 0)
	    *flag = 1;
	ChnFree(&commchn, comm);
    }
    /* Convert EUI active/undefined (-1) to MPMY "false" flag.  */
    if (*flag == EUI_ACTIVE)
	*flag = 0;
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int MPMY_Wait(MPMY_Comm_request req, MPMY_Status *stat) {
    struct comm_s *comm = req;
    int inbytes;

    Msgf(("Wait for %d\n", comm->hndl));
    if (mpc_wait(&comm->hndl, &inbytes) < EUI_SUCCEED)
	Error("MPMY_Wait mpc_wait failed\n");
    Msgf(("Waited for (%s), deallocated\n", 
	  (comm->inout==IN)?"in":"out"));
    if(comm->inout == IN && stat) {
	stat->count = inbytes;
	stat->src = comm->src;
	stat->tag = comm->tag;
    }
    ChnFree(&commchn, comm);
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SHIFT
#define SHIFT_TAG 0x1492
int MPMY_Shift(int proc, void *recvbuf, int recvcnt, 
	       const void *sendbuf, int sendcnt, MPMY_Status *stat) {
    int count;

    Msgf(("Starting MPMY_Shift(proc=%d, recvcnt=%d, sendcnt=%d, recvbuf=%p, sendbuf=%p\n",
	  proc, recvcnt, sendcnt, recvbuf, sendbuf));

    if (mpc_bsendrecv(sendbuf, sendcnt, proc, SHIFT_TAG,
		      recvbuf, recvcnt, &proc, &count) < EUI_SUCCEED)
	Error("MPMY_Shift mpc_bsendrecv failed\n");
    Msgf(("MPMY_Shift done, received=%d\n", count));
    if (stat) {
	stat->count = count;
	stat->src = proc;
	stat->tag = SHIFT_TAG;
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SYNC
int MPMY_Sync(void) {
    mpc_sync(ALLGRP);
    return MPMY_SUCCESS;
}

int MPMY_Init(int *argcp, char ***argvp) {
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    mpc_environ(&_MPMY_nproc_, &_MPMY_procnum_);
    _MPMY_initialized_ = 1;
    system("uptime");
    return MPMY_SUCCESS;
}

#include "timers_readrtc.c"
#define CANT_USE_ALARM
#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
