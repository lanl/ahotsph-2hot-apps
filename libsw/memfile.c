/* This version of memfile has special code for the paragon to
   avoid NORMA crashes resulting from too much simultaneous I/O activity.

   In retrospect, what I ended up doing was implement a distributed
   semaporhe on node 0 by using an hrecv handler to enforce the critical
   section.  Crude but effective.

   Unfortunately, it's still not enough.  Even though we guarantee that
   writes do not overlap, and that node P does fflush() before node P'
   attempts to do fwrite, we still crash due to NORMA.  My guess is that
   fflush() returns to node P before the i/o proc is has really had a
   chance to "catch up" (whatever that means). 

   Next attempt:  add a call to fsync() to the output code.

   Alternative:  actually send the data through to node 0 for output.
   Then we KNOW we can't be overwhelming the i/o proc because we can't
   send data any faster than we can move it through node0. */


#include <stdio.h>
#include <string.h>
#include "bigmalloc.h"
#include "Msgs.h"
#include "error.h"
#include "protos.h"
#include "mpmy.h"

/* Here we try to implement a cicular memory buffer which we can use as */
/* the vfprintf-like arg to Msg_init */

#if 1				/* IPD-hack? */
static char *memfile;
static int memfile_offset;
static int memfile_bufsz;

#ifdef __PARAGON__
static void HrecvSetup(int);
#endif

void
memfile_init(int sz)
{
    memfile = Malloc(sz);
    memfile_offset = 0;
    memfile_bufsz = sz-1;
    memfile[sz-1] = 0;   /* final null to make wrapped output cleaner */
#ifdef __PARAGON__ 
    HrecvSetup(sz);
#endif
}

void
memfile_delete(void)
{
    Free(memfile);
    memfile_offset = 0;
    memfile_bufsz = 0;
}
#else  /* IPD-hack? */
/* IPD will not output the buffer in a reasonable manner using code above */
char memfile[65536];
int memfile_offset;
int memfile_bufsz = 65535;
void memfile_init(){
}
void memfile_delete(){
}
#endif /* IPD-hack? */

#define BUFSZ 1024

void
memfile_vfprintf(void *junk, const char *fmt, va_list args)
{
    char tbuf[BUFSZ];	/* This might overflow, but msgs should be small */
    int i, len;

    vsprintf(tbuf, fmt, args);
    len = strlen(tbuf);
    if (len >= BUFSZ) Error("Buffer overflowed in mem_vfprintf\n");
    for (i = 0; i <= len; i++) {
	memfile[(memfile_offset+i)%memfile_bufsz] = tbuf[i];
    }
    memfile_offset += len;
}

void
PrintMemfile(void)
{
    if (memfile_offset == 0) return;
    printf("----- Messages from procnum %d -----\n", MPMY_Procnum());
    if (memfile_offset < memfile_bufsz) {
	printf("%s\n", memfile);
    } else {
	printf("Buffer has wrapped\n");
	printf("%s\n%s\n", memfile+(memfile_offset+1)%memfile_bufsz,
	       memfile);
    }
    fflush(stdout);
}

#ifndef __PARAGON__
/* I used lots of paragon-specific code in the implementation below, so
   non-paragon callers won't actually go 'via' node 0. */

void PrintMemfileVia0(void)
{
    PrintMemfile();
}

#else

/* The paragon will crash if all nodes do PrintMemfile simultaneously.
   The dreaded NORMA bug.  So what to do?  */
#include <nx.h>
#include <errno.h>

/* It's >>WAY<< too late to be thinking about Msgs.  But this code is
   twisted enough that the only way to debug it is to print messages.
   So we use the following HACK to run printf.  This causes a lot of
   useless clutter, and the clock is ticking (60sec), so use with
   discretion. */
#define Xprintf if( Msg_test(__FILE__) ) printf

/* 'random' numbers.  Hoping they don't conflict with anything... */
#define PRTMEMTYPE 48923
#define PRTMEMTYPE1 (PRTMEMTYPE + 1)
#define PRTMEMTYPE2 (PRTMEMTYPE + 2)
#define PRTMEMTYPE3 (PRTMEMTYPE + 3)

static int hrecv_initized = 0;
static int node_0_inbuf;

/* This is the hrecv handler set up only in node 0 for PRTMEMTYPE1 messages.
   It handshakes with the requesting node, allowing the requestor to dump
   to stdout without interfering with other processors. */
static void
PrintMemfileHndlr(long type, long count, long node, long ptype){
    int ok, done;

    if( count != sizeof(int) || type != PRTMEMTYPE1 )
	/* Returning will cause the sender to hang.  But what else
	   can we do.  His alarm clock will wake him... */
	return;
    Xprintf("PMH: for node=%d\n", node);	/* meta-debug */
    ok = 1;
    csend(PRTMEMTYPE2, &ok, sizeof(int), node, 0);
    Xprintf("PMH: node=%d, handshake done!\n", node); /* meta-debug */
    crecv(PRTMEMTYPE3, &done, sizeof(int));
    Xprintf("PMH: node=%d, output done!\n", node); /* meta-debug */
    /* Arm the handler again! */
    hrecv(PRTMEMTYPE1, &node_0_inbuf, sizeof(int), PrintMemfileHndlr); 
}

/* This is the alternative version of PrintMemfile.  Perhaps it should
   simply replace the other one?  */
void
PrintMemfileVia0(void){
    int begin, reply, done;
    long smid, rmid;
    int rdone, sdone;
    double give_up_time;
    int start;

    /* This probably means the caller never set up a handler on node 0.
       The whole exercise is pointless */
    if( !hrecv_initized )
	return;

    /* Bail out in 60 sec. */
    give_up_time = dclock() + 60.0;

    Xprintf("Printmemfilevia0 on proc %d, give up at %g\n", mynode(), give_up_time); /* meta-debug */
    begin = 1;
    smid = isend(PRTMEMTYPE1, &begin, sizeof(int), 0, 0);
    rmid = irecv(PRTMEMTYPE2, &reply, sizeof(int));
    sdone = rdone = 0;
    do{
	if(!sdone) sdone = msgdone(smid);
	if(!rdone) rdone = msgdone(rmid);
    }while(dclock() < give_up_time && (!sdone || !rdone));
    Xprintf("Printmemfilevia0 on proc %d, sdone=%d, rdone=%d, time=%g\n", 
	   mynode(), sdone, rdone, dclock()); /* meta-debug */

    if( sdone && rdone ){
	printf("--- Begin memfile from node %d ---\n", mynode());
	if (memfile_offset < memfile_bufsz) {
	    fwrite(memfile, 1, memfile_offset, stdout);
	} else {
	    printf("Buffer has wrapped\n");
	    start = (memfile_offset+1)%memfile_bufsz;
	    fwrite(memfile+start, 1, memfile_bufsz-start, stdout); 
	    printf("-- WRAP --");
	    start = memfile_offset%memfile_bufsz;
	    fwrite(memfile, 1, start, stdout);
	}
	printf("\n--- End memfile from node %d ---\n", mynode());
	fflush(stdout);
	if( fsync(fileno(stdout)) < 0 )
	    Xprintf("fsync fails on node %d, errno=%d\n", mynode(), errno); /* meta-debug */
	done = 1;
	csend(PRTMEMTYPE3, &done, sizeof(int), 0, 0);
    }else{
	printf("Node %d:  No confirmation from node 0.  Abandoning PrintMemfile\n", mynode());
    }
}

static void
HrecvSetup(int sz){
    if( mynode() == 0 ){
	hrecv(PRTMEMTYPE1, &node_0_inbuf, sizeof(int), PrintMemfileHndlr);
    }
    hrecv_initized = 1;
    gsync();
}

#endif /* __PARAGON__ */
