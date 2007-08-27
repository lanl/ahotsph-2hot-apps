/* batch.c:  Collect a series of small sends into larger ones */

#include "batch.h"
#include "mpmy.h"
#include "stk.h"
#include "bigmalloc.h"
#include "Msgs.h"

void PollWait(MPMY_Comm_request req, int tag);

static Stk **stks;
static int tag;
static int batch_size;

void
SetupBatch(int ttag, int size)
{
    int dest;

    tag = ttag;
    batch_size = size;
    stks = Calloc(MPMY_Nproc(), sizeof(Stk *));
    /* allocate all memory beforehand */
    /* Otherwise the incoming poll buffer will fight with the batch stks */
    /* for heap space, and we end up with a bunch of holes in the heap */
    for (dest = 0; dest < MPMY_Nproc(); dest++) {
	stks[dest] = Malloc(sizeof(Stk));
	StkInit(stks[dest], batch_size, Realloc_f, 0);
    }
}

void
FinishBatch(void)
{
    int i;
    MPMY_Comm_request req;
    int nproc = MPMY_Nproc();
    int procnum = MPMY_Procnum();

    for (i = 0; i < nproc; i++) {
	int dest = (procnum + i) % nproc;
	Stk *s = stks[dest];
	if (StkSz(s) > 0) {
	    MPMY_Isend(StkBase(s), StkSz(s), dest, tag, &req);
	    PollWait(req, tag);
	}
	StkTerminate(s);
	Free(s);
    }
    Free(stks);
}

void
SendBatch(void *outbuf, int size, int dest)
{
    Stk *s = stks[dest];

    StkPushData(s, outbuf, size);
    if (StkSz(s) > batch_size - size) {
	MPMY_Comm_request req;

	Msgf(("SendBatch: %d to %d\n", StkSz(s), dest));
	MPMY_Isend(StkBase(s), StkSz(s), dest, tag, &req);
	PollWait(req, tag);
	Msgf(("SendBatch: waited\n"));
	StkClear(s);
    }
}
