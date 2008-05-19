#ifdef _AIX
#define Malloc MPI_Malloc
#define Realloc MPI_Realloc
#define Calloc MPI_Calloc
#define Free MPI_Free
#endif

#ifdef _SWAMPI
#include <swampi.h>
#else
#include <mpi.h>
#endif

#ifdef _AIX
#undef Malloc
#undef Realloc
#undef Calloc
#undef Free
#endif
#include "bigmalloc.h"

#include "chn.h"
#include "mpmy.h"
#include "Assert.h"
#include "timers.h"
#include "Msgs.h"
#include "error.h"

#ifdef _AIX
/* This is some crap that mpcc creates */
int mpmondata=0; int mp_linked_euilib=0;
#endif

struct comm_s {
    MPI_Request hndl;
    int inout;
};

static Chn commchn;
#define NCOMM 128
#define IN 1
#define OUT 2

int MPMY_Isend(const void *buf, int cnt, int dest, int tag,
	       MPMY_Comm_request *reqp) {
    struct comm_s *comm = ChnAlloc(&commchn);

    Msgf(("Isend: buf=%p, dest=%d, tag=%d\n",
	  buf, dest, tag));
    if (MPI_Isend(buf, cnt, MPI_BYTE, dest, tag, MPI_COMM_WORLD,
		  &comm->hndl) != MPI_SUCCESS)
	Error("MPMY_Isend MPI_Isend failed\n");
    comm->inout = OUT;
    Msgf(("Isend: hndl=%d\n", (int) comm->hndl));
    *reqp = comm;
    return MPMY_SUCCESS;
}

#if 0
#define HAVE_MPMY_IRSEND
int MPMY_Irsend(const void *buf, int cnt, int dest, int tag,
	       MPMY_Comm_request *reqp) {
    struct comm_s *comm = ChnAlloc(&commchn);

    Msgf(("Irsend: buf=%p, dest=%d, tag=%d\n",
	  buf, dest, tag));
    if (MPI_Irsend(buf, cnt, MPI_BYTE, dest, tag, MPI_COMM_WORLD,
		   &comm->hndl) != MPI_SUCCESS)
	Error("MPMY_Isend MPI_Irsend failed\n");
    comm->inout = OUT;
    Msgf(("Irsend: hndl=%d\n", (int) comm->hndl));
    *reqp = comm;
    return MPMY_SUCCESS;
}
#endif /* 0 */

int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp) {
    struct comm_s *comm = ChnAlloc(&commchn);

    if (tag == MPMY_TAG_ANY) tag = MPI_ANY_TAG;
    if (src == MPMY_SOURCE_ANY) src = MPI_ANY_SOURCE;
    Msgf(("Irecv: buf=%p, src=%d, tag=%d\n",
	  buf, src, tag));
    if (MPI_Irecv(buf, cnt, MPI_BYTE, src, tag, MPI_COMM_WORLD,
		  &comm->hndl) != MPI_SUCCESS)
	Error("MPMY_Irecv MPI_Irecv failed\n");
    comm->inout = IN;
    Msgf(("Irecv: hndl=%d\n", (int) comm->hndl));
    *reqp = comm;
    return MPMY_SUCCESS;
}

int MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat) {
    struct comm_s *comm = req;
    MPI_Status status;
    int cnt;

    Msgf(("MPMY_Test hndl=%d\n",  (int) comm->hndl));
    if (MPI_Test(&comm->hndl, flag, &status) != MPI_SUCCESS)
	Error("MPMY_Test MPI_Test failed\n");
    Msgf(("Tested (%s), %d\n", 
	  (comm->inout==IN)?"in":"out", *flag));
    if (*flag) {
	if(comm->inout == IN) {
	    MPI_Get_count(&status, MPI_BYTE, &cnt);
	    Msgf(("Recvd(T) from %d, tag %d, count: %d\n",
		  status.MPI_SOURCE, status.MPI_TAG, cnt));
	    if (stat) {
		stat->src = status.MPI_SOURCE;
		stat->tag = status.MPI_TAG;
		stat->count = cnt;
	    }
	}
	ChnFree(&commchn, comm);
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int MPMY_Wait(MPMY_Comm_request req, MPMY_Status *stat) {
    struct comm_s *comm = req;
    MPI_Status status;
    int cnt;

    Msgf(("Wait for %d\n", (int) comm->hndl));
    if (MPI_Wait(&comm->hndl, &status) != MPI_SUCCESS)
	Error("MPMY_Wait MPI_Wait failed\n");
    Msgf(("Waited for (%s), deallocated\n", 
	  (comm->inout==IN)?"in":"out"));
    if(comm->inout == IN) {
	MPI_Get_count(&status, MPI_BYTE, &cnt);
	Msgf(("Recvd(W) from %d, tag %d, count: %d\n", 
	      status.MPI_SOURCE, status.MPI_TAG, cnt));
	if (stat) {
	    stat->src = status.MPI_SOURCE;
	    stat->tag = status.MPI_TAG;
	    stat->count = cnt;
	}
    }
    ChnFree(&commchn, comm);
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SHIFT
#define SHIFT_TAG 0x1492
int MPMY_Shift(int proc, void *recvbuf, int recvcnt, 
	       const void *sendbuf, int sendcnt, MPMY_Status *stat) {
    MPI_Status status;
    int count;

    Msgf(("Starting MPMY_Shift(proc=%d, recvcnt=%d, sendcnt=%d, recvbuf=%p, sendbuf=%p\n",
	  proc, recvcnt, sendcnt, recvbuf, sendbuf));

    if (MPI_Sendrecv(sendbuf, sendcnt, MPI_BYTE, proc, SHIFT_TAG,
		     recvbuf, recvcnt, MPI_BYTE, proc, SHIFT_TAG,
		     MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Error("MPMY_Shift MPI_Sendrecv failed\n");
    MPI_Get_count(&status, MPI_BYTE, &count);
    Msgf(("MPMY_Shift done, received=%d\n", count));
    if (stat) {
	stat->count = count;
	stat->src = status.MPI_SOURCE;
	stat->tag = status.MPI_TAG;
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SYNC
int MPMY_Sync(void) {
    MPI_Barrier(MPI_COMM_WORLD);
    return MPMY_SUCCESS;
}

int MPMY_Init(int *argcp, char ***argvp) {
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    MPI_Initialized(&_MPMY_initialized_);
    if (_MPMY_initialized_ == 0)
	if (MPI_Init(argcp, argvp) != MPI_SUCCESS)
	    Error("MPMY_Init MPI_Init failed\n");
    if (MPI_Comm_size(MPI_COMM_WORLD, &_MPMY_nproc_) != MPI_SUCCESS)
	Error("MPMY_Init MPI_Comm_size failed\n");
    if (MPI_Comm_rank(MPI_COMM_WORLD, &_MPMY_procnum_) != MPI_SUCCESS)
	Error("MPMY_Init MPI_Comm_rank failed\n");
    _MPMY_initialized_ = 1;
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_FINALIZE
int MPMY_Finalize(void){
  return (MPI_Finalize() == MPI_SUCCESS) ? MPMY_SUCCESS : MPMY_FAILED ;
}

#ifndef _AIX	/* We need a CPU timer, and MPI doesn't have one */
#define HAVE_MPMY_TIMERS
#include "timers_mpi.c"
#endif

#if defined(__CM5__) || defined(_AIX) || defined(__AP1000__)
#define CANT_USE_ALARM
#endif
#if defined(__CM5__) || defined(__INTEL_SSD__)
#include "mpmy_pario.c"
#else
#if defined(USE_MPIIO)
#include "mpmy_mpiio.c"
#else
#include "mpmy_io.c"
#endif
#endif
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
