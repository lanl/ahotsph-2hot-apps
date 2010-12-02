#include <string.h>

#include "mpmy.h"
#include "Msgs.h"
#include "bigmalloc.h"
#include "gc.h"			/* for ilog2 */
#include "verify.h"

/* Could we implement gather using MPMY_Combine with MPMY_Op == MPMY_GATHER? */


#define BCAST_DEFAULT_TAG 0x47

#define GATHER_BCAST_TAG 0x1145
#define GATHER_TAG 0x2145
#define NGATHER_TAG 0x3145

/* Should this be public? */
static unsigned int MPMY_Datasize[] = 
{ sizeof(float), sizeof(double), sizeof(int), sizeof(char), sizeof(short),
  sizeof(long), sizeof(unsigned int), sizeof(unsigned char), 
  sizeof(unsigned short), sizeof(unsigned long), 1/*user_data*/
};

void
MPMY_send(const void *buf, int cnt, int dest, int tag)
{
    MPMY_Comm_request req;

    MPMY_Isend(buf, cnt, dest, tag, &req);
    MPMY_Wait(req, 0);
    if( Msg_test(__FILE__)){
	int i;
	int sum = 0;
	const char *cbuf = buf;
	for(i=0; i<cnt; i++){
	    sum ^= cbuf[i];
	}
	Msg_do("mpmy_gather: send(cnt=%d, dest=%d, tag=%d), sum=%d\n", 
	       cnt, dest, tag, sum);
    }

}

void
MPMY_recvn(void *buf, int cnt, int src, int tag)
{
    MPMY_Status stat;
    MPMY_Comm_request req;

    Verify(MPMY_Irecv(buf, cnt, src, tag, &req)==MPMY_SUCCESS);
    Verify(MPMY_Wait(req, &stat)==MPMY_SUCCESS);
    if (MPMY_Count(&stat) != cnt) 
      Error("Recv failed, expected %d got %d\n", cnt, MPMY_Count(&stat));
    if( Msg_test(__FILE__)){
	int i;
	int sum = 0;
	char *cbuf = buf;
	for(i=0; i<cnt; i++){
	    sum ^= cbuf[i];
	}
	Msg_do("mpmy_gather: recvn(cnt=%d, dest=%d, tag=%d), sum=%d\n", 
	       cnt, src, tag, sum);
    }
}

int
MPMY_AllGather(const void *sendbuf, int count, MPMY_Datatype type, 
	       void *recvbuf)
{
    int chan;
    int doc;
    MPMY_Status stat;
    int sendproc;
    int bufsz;
    unsigned int mask;
    int ret;
    void *inptr;
    const void *outptr;
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int nbytes = MPMY_Datasize[type] * count;

    doc = ilog2(nproc);
    if (nproc != 1 << doc) {	/* for non power-of-two sizes */
	MPMY_Gather(sendbuf, count, type, recvbuf, 0);
	MPMY_BcastTag(recvbuf, nproc*count, type, 0, GATHER_BCAST_TAG);
	return MPMY_SUCCESS;
    }

    inptr = (char *)recvbuf + procnum * nbytes;
    outptr = sendbuf;

    /* Self contribution */
    if (inptr != outptr)
      memcpy(inptr, outptr, nbytes);

    mask = ~0;
    for (chan = 0; chan < doc; chan++) {
	sendproc = procnum ^ (1 << chan);
	bufsz = (1 << chan) * nbytes;
	inptr = (char *)recvbuf + (sendproc & mask) * nbytes;
	outptr = (char *)recvbuf + (procnum & mask) * nbytes;
	MPMY_Shift(sendproc, inptr, bufsz, outptr, bufsz, &stat);
	ret = MPMY_Count(&stat);
	if (ret != bufsz) 
	  Error("Shift failed, expected %d got %d\n", bufsz, ret);
	mask <<= 1;
    }
    return MPMY_SUCCESS;
}

int
MPMY_Gather(const void *sendbuf, int count, MPMY_Datatype type, void *recvbuf,
	    int recvproc)
{
    int chan;
    int doc;
    int sendproc;
    int inbufsz, outbufsz;
    int nin, nout;
    unsigned int mask;
    void *inptr;
    const void *outptr;
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int nbytes = MPMY_Datasize[type] * count;

    if (recvproc != 0) Error("Gather to procnum != 0 not supported yet.\n");

    doc = ilog2(nproc);
    if (nproc != 1 << doc)
      doc++;			/* for non power-of-two sizes */

    if (procnum != recvproc)
      recvbuf = Malloc(nproc*nbytes); /* overkill */

    inptr = (char *)recvbuf + procnum * nbytes;
    outptr = sendbuf;
    if (inptr != outptr) memcpy(inptr, outptr, nbytes);

    mask = ~0;
    for (chan = 0; chan < doc; chan++) {
	sendproc = procnum ^ (1 << chan);
	if (sendproc >= 0 && sendproc < nproc) {
	    nin = nout = (1 << chan);
	    if (procnum & (1 << chan)) {
		if (nproc - procnum < nout) nout = nproc - procnum;
		outbufsz = nout * nbytes;
		outptr = (char *)recvbuf + (procnum & mask) * nbytes;
		Msgf(("Gather: to %d, outidx %d, outsz %d\n", 
		       sendproc, (procnum & mask), nout));
		MPMY_send(outptr, outbufsz, sendproc, GATHER_TAG+chan);
		break;
	    } else {
		if (nproc - sendproc < nin) nin = nproc - sendproc;
		inbufsz = nin * nbytes;
		inptr = (char *)recvbuf + (sendproc & mask) * nbytes;
		Msgf(("Gather: from %d, inidx %d, insz %d\n", 
		       sendproc, (sendproc & mask), nin));
		MPMY_recvn(inptr, inbufsz, sendproc, GATHER_TAG+chan);
	    }
	}
	mask <<= 1;
    }
    if (procnum != recvproc)
      Free(recvbuf);
    return MPMY_SUCCESS;
}

/* "count" can vary for each processor in NGather */

int
MPMY_NGather(const void *sendbuf, int count, MPMY_Datatype type, 
	     void **recvhndl, int recvproc)
{
    int chan;
    int doc;
    int sendproc;
    long bufsz;
    long inbytes;
    unsigned int mask;
    void *buf;
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int nbytes = MPMY_Datasize[type] * count;

    if (recvproc != 0) Error("NGather to procnum != 0 not supported yet.\n");

    doc = ilog2(nproc);
    if (nproc != 1 << doc)
      doc++;			/* for non power-of-two sizes */

    buf = Malloc(nbytes); 
    memcpy(buf, sendbuf, nbytes);
    bufsz = nbytes;

    mask = ~0;
    for (chan = 0; chan < doc; chan++) {
	sendproc = procnum ^ (1 << chan);
	if (sendproc >= 0 && sendproc < nproc) {
	    if (procnum & (1 << chan)) {
		Msgf(("NGather: to %d, bufsz %d\n", sendproc, bufsz));
		MPMY_send(&bufsz, sizeof(long), sendproc, NGATHER_TAG+chan);
		MPMY_send(buf, bufsz, sendproc, NGATHER_TAG+chan);
		break;
	    } else {
		MPMY_recvn(&inbytes, sizeof(long), sendproc, NGATHER_TAG+chan);
		Msgf(("NGather: from %d, inbytes %d\n", sendproc, inbytes));
		buf = Realloc(buf, bufsz+inbytes);
		MPMY_recvn(bufsz+(char*)buf, inbytes, sendproc, NGATHER_TAG+chan);
		bufsz += inbytes;
	    }
	}
	mask <<= 1;
    }
    if (procnum != recvproc) {
	Free(buf);
	return 0;
    } else {
	Msgf(("NGather: Final bufsz %d\n", bufsz));
	*recvhndl = buf;
	return bufsz/MPMY_Datasize[type];
    }
}

int MPMY_Bcast(void *buf, int count, MPMY_Datatype type, int srcproc)
{
    return MPMY_BcastTag(buf, count, type, srcproc, BCAST_DEFAULT_TAG);
}

int
MPMY_BcastTag(void *buf, int count, MPMY_Datatype type, int srcproc, int tag)
{
    int chan;
    int doc;
    int sendproc;
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int nbytes = MPMY_Datasize[type] * count;

    if (srcproc != 0) Error("Bcast from procnum != 0 not supported yet.\n");

    Msgf(("MPMYBcast(buf=%p, count=%d, type=%d, srcproc=%d\n",
	  buf, count, type, srcproc)); Msg_flush();
    doc = ilog2(nproc);
    if (nproc != 1 << doc)
      doc++;			/* for non power-of-two sizes */

    if( procnum != 0 ){		/* DEBUG ONLY! */
	unsigned char *cbuf = buf;
	cbuf[nbytes-1] = procnum + 128;
    }

    for (chan = 0; chan < doc; chan++) {
	sendproc = procnum ^ (1 << chan);
	if (sendproc >= 0 && sendproc < nproc) {
	    if (procnum & (1 << chan)) {
		Msgf(("Bcast: recv from %d\n", sendproc));
		MPMY_recvn(buf, nbytes, sendproc, tag+chan);
	    } else {
		Msgf(("Bcast: send to %d\n", sendproc));
		MPMY_send(buf, nbytes, sendproc, tag+chan);
	    }
	}
    }
    return MPMY_SUCCESS;
}

