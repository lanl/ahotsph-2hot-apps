#include "mpmy.h"
#include "Msgs.h"
#include "assert.h"
#include "gc.h"
#include "timers.h"
#include "error.h"
#include "dll.h"
#include "chn.h"

#define INBUFSZ (256*sizeof(int))

Timer_t PollWaitTm;

static int done;
static Dll dlldone;
static Chn dlldonechn;
static void DllWait(Dll *dll);
static void NLDone(int who, int tag);
static void IBcast(const void *buf, int count, int tag);

static void (*func)();
static int size;
static MPMY_Comm_request inreq;
static int allbitsdone;
static int inbuf[INBUFSZ/sizeof(int)]; /* avoid using malloc */
static const int junk[2];

void
PollSetup(void put(void *buf, int size), int max_size, int tag)
{
    int doc;

    doc = ilog2(MPMY_Nproc());
    if (MPMY_Nproc() != 1 << doc)
      doc++;			/* for non power-of-two sizes */
    allbitsdone = (1 << (doc+1))-1;
    done = 0;
    func = put;
    if (max_size > INBUFSZ) SinglError("INBUFSZ too small\n");

    /* In fact, this test is insufficient if somebody decides to send
       a short message anyway we'll still be confused! */
    if( max_size == sizeof(int) || max_size == 2*sizeof(int) )
	SinglError("Poll uses size for message sorting.  You can't use size=%ld or %ld without some new coding\n", (long)sizeof(int), (long)2*sizeof(int));

    size = max_size;
#if 0
    /* This generates a lot of Shouts from MPMY */
    {int flag;
    if (MPMY_Test(inreq, &flag, 0) != MPMY_FAILED)
      Error("PollSetup found leftover Irecv\n");}
#endif
    /* This will call malloc.  Is that a problem?? */
    DllCreateChn(&dlldonechn, sizeof(MPMY_Comm_request), doc);
    DllCreate(&dlldone, &dlldonechn);
    MPMY_Irecv(&inbuf, size, MPMY_SOURCE_ANY, tag, &inreq);
}

void
Poll(int tag)
{
    int flag;
    MPMY_Status stat;
    int reset=0;
    
    Msgf(("P(tag=%d)\n", tag));
    if( Msg_test(__FILE__) && !Msg_test("mpmy_cm5.c") ){
	Msg_on("mpmy_cm5.c");
	reset = 1;
    }
    MPMY_Flick();
    while (MPMY_Test(inreq, &flag, &stat), flag) {
	if (stat.count == sizeof(int)) {
	    NLDone(stat.src, tag);
	} else {
	    func(&inbuf, stat.count);
	} 
	MPMY_Irecv(&inbuf, size, MPMY_SOURCE_ANY, tag, &inreq);
    }
    if( reset ){
	Msg_off("mpmy_cm5.c");
    }
}

/* If we use a plain MPMY_Wait() during a poll session, deadlock may */
/* result from isends blocking */
void
PollWait(MPMY_Comm_request req, int tag)
{
    int flag;
    MPMY_Status stat;
    int reset=0;

    if( Msg_test(__FILE__) && !Msg_test("mpmy_cm5.c") ){
	Msg_on("mpmy_cm5.c");
	reset = 1;
    }
    Msgf(("PW(tag=%d)\n", tag));
    while (1) {
	while (MPMY_Test(inreq, &flag, &stat), flag) {
	    if (stat.count == sizeof(int)) {
		NLDone(stat.src, tag);
	    } else {
		func(&inbuf, stat.count);
	    } 
	    MPMY_Irecv(&inbuf, size, MPMY_SOURCE_ANY, tag, &inreq);
	}
	if (MPMY_Test(req, &flag, 0), flag) return;
	MPMY_Flick();
    }
    if( reset ){
	Msg_off("mpmy_cm5.c");
    }
}

void
PollUntilDone(int tag)
{
    static MPMY_Comm_request req;
    MPMY_Status stat;

    Msgf(("PUD(tag=%d)\n", tag));
    NLDone(MPMY_Procnum(), tag);

    StartTimer(&PollWaitTm);
    while (done != allbitsdone || MPMY_Procnum() != 0) {
	MPMY_Wait(inreq, &stat);
	if (stat.count == 2*sizeof(int)) {
	    /* This must have come from IBcast */
	    break;
	} else if (stat.count == sizeof(int)) {
	    NLDone(stat.src, tag);
	} else {
	    func(&inbuf, stat.count);
	}
	MPMY_Irecv(&inbuf, size, MPMY_SOURCE_ANY, tag, &inreq);
    }
    if (MPMY_Procnum() == 0) {
	/* cancel last Irecv */
	MPMY_Isend(junk, sizeof(int), 0, tag, &req); /* self-send */
	/* Could this deadlock?  The send and recv have both been posted. */
	/* It must complete eventually, right? */
	MPMY_Wait(req, 0);
	MPMY_Wait(inreq, 0);
    }
    Msgf(("done\n"));

    IBcast(junk, 2, tag);
    DllWait(&dlldone);
    StopTimer(&PollWaitTm);
}

static void NLDone(int who, int tag) {
    int procnum= MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int mask, mask2;
    int reset=0;

    if( Msg_test(__FILE__) && !Msg_test("mpmy_cm5.c") ){
	Msg_on("mpmy_cm5.c");
	reset = 1;
    }
    Msgf(("poll.c:NLDone(%d, %d)\n", who, tag));
    mask = who ^ procnum;
    if (mask == 0) mask = 1;
    else mask <<= 1;
    done |= mask;

    /* Make sure that all the bits in 'done' below mask are set */
    mask2 = mask-1;	/* all bits below mask */
    if( (mask2 & done) != mask2 )
	goto done;

    while ( (mask & done) && (mask<nproc) ) {
	int dest = procnum^mask;
	if (dest < nproc) {	/* for non power-of-two sizes */
	    Msgf(("NLDone:  informing %d\n", dest));
	    MPMY_Isend(junk, sizeof(int), dest, tag, 
		       DllData(DllInsertAtBottom(&dlldone)));
	} else {
	    Msgf(("NLDone:  ghost informed %d\n", dest));
	    done |= (dest ^ procnum) << 1;
	}
	mask <<= 1;
    }
 done:
    if( reset ){
	Msg_off("mpmy_cm5.c");
    }
    return;
}

static void DllWait(Dll *dll){
    Dll_elmt *p;
    int flag;

    while( DllLength(dll) > 0 ){
	MPMY_Flick();
	for(p=DllTop(dll); p!=DllInf(dll); p = DllDown(p)){
	    if( MPMY_Test(*(MPMY_Comm_request *)DllData(p), &flag, 0) != MPMY_SUCCESS ){
		Error("MPMY_Test failed on %lx\n", (unsigned long)DllData(p));
	    }
	    if( flag ){
		p = DllDeleteUp(dll, p);
	    }
	}
    }
}

static void
IBcast(const void *buf, int count, int tag)
{
    int chan, doc;
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();

    Msgf(("IBcast(tag=%d)\n", tag));
    doc = ilog2(MPMY_Nproc());
    if (MPMY_Nproc() != 1 << doc)
      doc++;			/* for non power-of-two sizes */

    for (chan = hibit(MPMY_Procnum())+1; chan < doc; chan++) {
	int sendproc = procnum ^ (1 << chan);
	if (!(procnum & (1 << chan)) && sendproc < nproc) {
	    Msgf(("IBcast chan=%d, sendproc=%d\n", chan, sendproc));
	    MPMY_Isend(buf, count*sizeof(int), sendproc, tag,
		       DllData(DllInsertAtBottom(&dlldone)));
	    Msgf(("IBcast Send done!\n"));
	}
    }
}

