#define NO_MSGS
#include "abm.h"
#include "tree.h"
#include "mpmy.h"
#include "key.h"
#include "Msgs.h"
#include "bigmalloc.h"
#include "error.h"
#include "chn.h"
#include "stk.h"
#include "timers.h"
#include "Assert.h"

/* Emit a warning if we Poll this many times without making any progress */
#define NOPROGRESS 100000

/* Make the packet size a couple of times bigger than the largest possible
   reply.  Don't sweat the details.  Just make it nice and big. */
#define REPLIES_PER_PKT 20 /* used to be 5. */
Counter_t DeferCnt;
Counter_t RequestCnt;
Timer_t WalkDeferTm;
		
/* Should there be a data structure with all this in it??? */
/* Should that data structure be tree_t?? */
static int level;
static int common_level;
static int done_first;
static Stk *unacc;
static void **sink_tbl;
static Stk def1, def2;
static Stk arrived;
static Stk walkstk1, walkstk2;
static hcell **hc_tbl;
static int max_pp_vec;
static int *result_vec;
static hcell **pp_vec;
static tree_t *Srctp, *Sinktp;
static macv_t MACv;
static inherit_t Inherit;
static float Sinksz;

#if NK==1
#define KEYNEQ(key1, key2) ((key1).k[0] != (key2).k[0])
#else
#define KEYNEQ(key1, key2) \
   ((key1).k[0] != (key2).k[0] || (key1).k[1] != (key2).k[1])
#endif

static ABM Abm;
static ABMhndlr_t reqhndlr;
static ABMhndlr_t replyhnlr;

/* These correspond to entries in walkhndlrs.  
   A REQUESTTYPE will invoke walkhndlr[REQUESTTYPE] (i.e., reqhndlr).
   A REPLYTYPE will invoke walkhndlr[REPLYTYPE] (i.e., replyhnlr)
   */
#define REQUESTTYPE 0
#define REPLYTYPE 1
/* Some compilers can't handle this initiaizer.  I'll put it in 
   setup for their benefit.  Who knows, maybe it's really non-ANSI... */
static ABMhndlr_t *walkhndlrs[2]; /* = {&reqhndlr, &replyhnlr}; */

/* The tag to use for all asynchronous walk messages. */
#define WALKTAG 3210
static ABMpktz_t copyContents;
static int szContents(hcell *pp);

#if 0 /* This didn't have any beneficial effect.  Perhaps more
	 experimentation with the nskip values would help */
/* Do it as a macro so every call has a 'private' defer_cnt */
#define ABMFlushMaybe(abmp, nskip) do{\
  static int defer_cnt = 0;\
  if( defer_cnt++ >= (nskip) ){\
    defer_cnt = 0;\
    ABMFlush(abmp);\
  }\
}while(0)
#else
#define ABMFlushMaybe(abmp, nskip) ABMFlush(abmp)
#endif

#define ABMPollMaybe(abmp, nskip) do{\
  static int defer_cnt = 0;\
  if( defer_cnt-- <= 0 ){\
    defer_cnt = (nskip);\
    if (ABMPoll(abmp) < 0) Error("ABMPoll failed\n");\
  }\
}while(0)

static void 
setupWalk(tree_t *srctp, tree_t *sinktp, int sinksz,
	  inherit_t inherit,
	  macv_t macV)
{
    int i, chubits2;

    level = 0;
    common_level = 0;
    done_first = 0;
    StkInitEz(&def1);
    StkInitEz(&def2);
    StkInitEz(&arrived);
    StkInitEz(&walkstk1);
    StkInitEz(&walkstk2);
    chubits2 = ((KEYBITS-1)/(sinktp->ndim)) + 2;
    sink_tbl = Malloc(chubits2*sizeof(void*));
    hc_tbl = Malloc(chubits2*sizeof(hcell *));
    unacc = Malloc(chubits2*sizeof(Stk));

    sink_tbl[0] = Malloc(chubits2*sinksz);
    inherit(NULL, sink_tbl[0], Find(sinktp, KeyInt(1)));
    for(i=0; i<chubits2; i++){
	if( i>0 ){
	    sink_tbl[i] = sinksz + (char *)(sink_tbl[i-1]);
	}
	StkInitEz(&unacc[i]);
	hc_tbl[i] = NULL;
    }
    StkPushType(&unacc[0], KeyInt(1), Key_t);
    max_pp_vec = 0;
    result_vec = NULL;
    pp_vec = NULL;
}

static void terminateWalk(void){
    int chubits2, i;

    StkTerminate(&def1);
    StkTerminate(&def2);
    StkTerminate(&arrived);
    StkTerminate(&walkstk1);
    StkTerminate(&walkstk2);
    chubits2 = ((KEYBITS-1)/(Sinktp->ndim)) + 2;
    for(i=0; i<chubits2; i++){
	StkTerminate(&unacc[i]);
    }
    Free(pp_vec);
    Free(result_vec);
    Free(sink_tbl[0]);
    Free(sink_tbl);
    Free(hc_tbl);
    Free(unacc);
    pp_vec = NULL;
    result_vec = NULL;
    sink_tbl = NULL;
    hc_tbl = NULL;
    unacc = NULL;
    max_pp_vec = 0;
}

static void
findv(const tree_t *tp, const Key_t *key_vec, hcellptr *out, const int n)
{
    hcellptr *firstp, *prevnext, np;
    int i;

    for (i = 0; i < n; i++) {
	firstp = tp->htab + (key_vec[i].k[0] & HASH_MASK);

	for (np = *(prevnext = firstp);
	     np && KEYNEQ(np->key, key_vec[i]);
	     np = *(prevnext = &np->next)){}

	out[i] = np;

	if (np != *firstp && np) { 
	    *prevnext = np->next; /* this shouldn't happen if np was at top */
	    np->next = *firstp;
	    *firstp = np;
	}
    }
}

static void 
walkv(void *sink, const Stk *parent_unacc, Stk *unacceptable, 
     Stk *deferred)
{
    Key_t key;
    hcell_type type;
    int n;
    const Stk *in;
    Stk *nextin;
    Key_t *key_vec;
    int i, nvec;
    int toggle;
    int ndim = Srctp->ndim;

    in = parent_unacc;
    nextin = &walkstk1;
    toggle = 0;
    while ((nvec = StkSz(in)/sizeof(Key_t)) > 0) {
        ABMPollMaybe(&Abm, 20);
	if (nvec > max_pp_vec) {
	    max_pp_vec = nvec + 512;
	    pp_vec = Realloc(pp_vec, max_pp_vec * sizeof(hcellptr));
	    result_vec = Realloc(result_vec, max_pp_vec * sizeof(int));
	}
	StkClear(nextin);
	key_vec = StkBase(in);
	findv(Srctp, key_vec, pp_vec, nvec);
	MACv(sink, (const hcell **)pp_vec, result_vec, nvec); /* why cast? */

	for (i = 0; i < nvec; i++) {
	    key = *(Key_t *) ((char *)key_vec + i * StkAlign(in, sizeof(Key_t)));
	    switch(result_vec[i]){
	    case MAC_SPLIT_SINK:
		StkPushType(unacceptable, key, Key_t);
		break;
	    case MAC_SPLIT_SRC:
		type = Type(pp_vec[i]);
	    
		if (TreeKidsOK(type)) {
		    key = KeyLshift(key, ndim);
		    for (n = Sub_Flags_Type(type); n; n >>= 1,key.k[0]++) {
			if ((n & 1) == 0) continue;
			StkPushType(nextin, key, Key_t);
		    }
		} else {
		    StkPushType(deferred, key, Key_t);
		    IncrCounter(&DeferCnt);
		    if ((type & REQUESTED) == 0) {
			Type(pp_vec[i]) |= REQUESTED;
			IncrCounter(&RequestCnt);
			ABMPost(&Abm, GetSource(type), sizeof(Key_t), 
				REQUESTTYPE, (ABMpktz_t *)&memcpy, 
				&(pp_vec[i]->key));
		    }
		}
		break;
	    }
	}
	/* Toggle in and nextin */
	if( toggle ){
	    in = &walkstk2;
	    nextin = &walkstk1;
	    toggle = 0;
	}else{
	    in = &walkstk1;
	    nextin = &walkstk2;
	    toggle = 1;
	}
    }
    ABMFlushMaybe(&Abm, 30);	/* every time we descend a level */
}

/* This version does deferrals only for the duration of the MACv. */
/* It might be worthwhile to defer for longer by keeping separate defer */
/* lists at each level and putting another loop around the whole thing. */
static void Walkbody(int commonlev, int bodylev){
    int l;
    Stk *deffrom, *defto, *deftmp;
    void *sink;
    int ndef, i;
    Key_t *key_vec;
    int foundone;

    for(l=commonlev; l<bodylev; l++){
	Inherit(sink_tbl[l], sink_tbl[l+1], hc_tbl[l+1]);
	sink = sink_tbl[l+1];
	StkClear(&def1);
	StkClear(&unacc[l+1]);
	ABMPollMaybe(&Abm, 200);
	walkv(sink, &unacc[l], &unacc[l+1], &def1);
	deffrom = &def1;
	defto = &def2;
	foundone = 1;		/* ???? could be 0 the first time ???  */
	while( StkSz(deffrom) ){
	    StkClear(defto);
	    ABMFlushMaybe(&Abm, 100); /* every time we are deferred */
	    /* We can't get out of here until something arrives, so
	       we might be better off calling PollWait, which >>might<<
	       be more sociable about letting somebody else have the CPU */
	    /* But we shouldn't wait if some new stuff has arrived since 
	       last time because there might be oodles of good stuff
	       available in memory to chew on... */
	    if( !foundone ){
		StartTimer(&WalkDeferTm);
		if( ABMPollWait(&Abm) < 0 )
		    Error("ABMPollWait failed\n");
		StopTimer(&WalkDeferTm);
	    }else{
	        ABMPollMaybe(&Abm, 200);
	    }
	    /* This code is almost identical to the code in walkv...*/
	    ndef = StkSz(deffrom)/sizeof(Key_t);
	    key_vec = StkBase(deffrom);
	    if (ndef > max_pp_vec) {
		max_pp_vec = ndef + 512;
		pp_vec = Realloc(pp_vec, max_pp_vec * sizeof(hcellptr));
		result_vec = Realloc(result_vec, max_pp_vec * sizeof(int));
	    }
	    findv(Srctp, key_vec, pp_vec, ndef);
	    StkClear(&arrived);
	    foundone = 0;
	    for(i=0; i<ndef; i++){
		int type;
		Key_t key;
		int ndim = Srctp->ndim;
		type = Type(pp_vec[i]);
		/* key = key_vec[i] */
		key = *(Key_t *)( (char *)key_vec + i*StkAlign(deffrom, sizeof(Key_t)));
		if(TreeKidsOK(type)){
		    int n;
		    foundone = 1;
		    key = KeyLshift(key, ndim);
		    for (n = Sub_Flags_Type(type); n; n >>= 1,key.k[0]++) {
			if ((n & 1) == 0) continue;
			StkPushType(&arrived, key, Key_t);
		    }
		} else {
		    StkPushType(defto, key, Key_t);
		    IncrCounter(&DeferCnt);
		    if (DeferCnt.counter > 10000000) {
			SeriousWarning("Defer count is very large\n");
			DeferCnt.counter = 0;
		    }
		    assert(type & REQUESTED);
		}
	    }
	    walkv(sink, &arrived, &unacc[l+1], defto);
	    deftmp = deffrom;
	    deffrom = defto;
	    defto = deftmp;
	}
    }
    if( TreeLocal(hc_tbl[bodylev]->type) ) 
	Inherit(sink_tbl[bodylev], NULL, hc_tbl[bodylev]);
    ABMFlushMaybe(&Abm, 2);		/* every time we finish a body */
}

static int preWalk(tree_t *tp, hcell *hp){
    level++;
    hc_tbl[level] = hp;
    if( Sub_Flags(hp) == 0 && TreeLocal(hp->type) ){
	if( Msg_test(__FILE__) ){
	    Msg_do("Body at %s.\n", hcellPrint(hp));
	    Msg_do("level=%d, common level=%d\n", level, common_level);
	}
	if( done_first ){
	    Walkbody(common_level, level);
	}else{
	    Walkbody(0, level);
	    done_first = 1;
	}
	common_level = --level;
	return 0;		/* don't look for kids */
    }else{
	if( TreeLocal(hp->type) ){
	    return 1;
	}else{
	    common_level = --level;
	    return 0;
	}
    }
}

static void postWalk(tree_t *tp, hcell *hp, hcell **daughters){
    common_level = --level;
}


void Walk(tree_t *srctp, tree_t *sinktp, int sinksz,
	  macv_t MAC, inherit_t InheritSink){
    WalkInit(srctp, sinktp, sinksz, MAC, InheritSink);
    WalkNT(sinktp);
    WalkTerminate();
}

void
WalkInit(tree_t *srctp, tree_t *sinktp, int sinksz, macv_t MAC, 
	 inherit_t InheritSink)
{
    walkhndlrs[REQUESTTYPE] = reqhndlr;
    walkhndlrs[REPLYTYPE] = replyhnlr;
    ABMSetup(&Abm, REPLIES_PER_PKT*(1<<srctp->ndim)*(srctp->cell_sz+srctp->tbody_sz),
	     WALKTAG, 2, walkhndlrs);
    Srctp = srctp;
    MACv = MAC;
    Inherit = InheritSink;
    Sinktp = sinktp;
    Sinksz = sinksz;
}

void 
WalkNT(tree_t *sinktp)
{
    setupWalk(Srctp, Sinktp, Sinksz, Inherit, MACv);
    Traverse(sinktp, Find(sinktp, KeyInt(1)), preWalk, postWalk);
    terminateWalk();
}

void 
WalkTerminate(void)
{
    ABMIamDone(&Abm);
    while( !ABMAllDone(&Abm) ){
	if( ABMPoll(&Abm) < 0 )
	    Error("ABMPoll fails while waiting for Alldone\n");
	ABMFlush(&Abm);
    }
    ABMShutdown(&Abm);
    Srctp = NULL;
    MACv = NULL;
}

/* How big will the message be that returns the "Contents" of pp? */
static int 
szContents(hcell *pp){
    unsigned int sub_flags, n;
    Key_t key, key0;
    int ndim = Srctp->ndim;
    int answer;

    answer = sizeof(Key_t);
    sub_flags = Sub_Flags(pp);
    key0 = KeyLshift(pp->key, ndim);
    answer += sizeof(int);
    for(n=0; sub_flags; n++, sub_flags>>=1){
	if( (sub_flags&1) == 0 )
	    continue;
	key = KeyOrInt(key0, n);
	pp = Find(Srctp, key);
	assert(pp);
	answer += sizeof(int);	/* 'type' */
	if( Sub_Flags(pp) )
	    answer += Srctp->cell_sz;
	else
	    answer += Srctp->tbody_sz;
    }
    Msgf(("szContents(%s)=%d\n", PrintKey(pp->key), answer));
    return answer;
}

static void 
copyContents(void *to, void *vp, int sz){
    hcell *pp = vp;
    unsigned int sub_flags, n;
    Key_t key, key0;
    int ndim = Srctp->ndim;
    int type;
    char *p = to;

    Msgf(("copyContents %s\n", PrintKey(pp->key)));

#if FORCE_KEY_ALIGNMENT
    memcpy(p, &pp->key, sizeof(Key_t));
#else
    *(Key_t *)p = pp->key;
#endif
    p += sizeof(Key_t);
    sub_flags = Sub_Flags(pp);
    key0 = KeyLshift(pp->key, ndim);
    *(int *)p = sub_flags;
    p += sizeof(int);
    for(n=0; sub_flags; n++, sub_flags>>=1){
	if( (sub_flags&1) == 0 )
	    continue;
	key = KeyOrInt(key0, n);
	pp = Find(Srctp, key);
	assert(pp);
	type = pp->type | DATAHERE;
	if( GetSource(pp->type) == -1 )
	    type |= NONLOCAL|PutSource(MPMY_Procnum());
	*(int *)p = type;
	p += sizeof(type);
	if( Sub_Flags(pp) ){
	    memcpy(p, pp->ptr, Srctp->cell_sz);
	    p += Srctp->cell_sz;
	}else{
	    memcpy(p, pp->ptr, Srctp->tbody_sz);
	    p += Srctp->tbody_sz;
	}
    }
    assert ((p - (char *)to) == sz );
}

static void 
reqhndlr(int src, int sz, void *p){
    Key_t key;
    hcell *pp;
    int repsz;

    Msgf(("WalkReply @%p %x %x \n", p, ((int *)p)[0], ((int *)p)[1]));
#if FORCE_KEY_ALIGNMENT
    memcpy(&key, p, sizeof(Key_t));
#else
    key = *(Key_t *)p;
#endif
    /* We have to Find it in order to get the size.  That's too bad because
       otherwise we could postpone the whole thing until we get called-back
       via copyContents */
    pp = Find(Srctp, key);
    assert(pp);
    repsz = szContents(pp);
    ABMPost(&Abm, src, repsz, REPLYTYPE, copyContents, pp);
}

static void 
replyhnlr(int src, int sz, void *vp){
    char *p = vp;		/* so we can do arithmetic */
    Key_t key0;
    Key_t key;
    unsigned int n, sub_flags;
    hcell *parent;
    int type;
    void *c;
    
#if FORCE_KEY_ALIGNMENT
    memcpy(&key0, p, sizeof(Key_t));
#else
    key0 = *(Key_t *)p;
#endif
    Msgf(("WalkReact to %s from %d\n", PrintKey(key0), src));
    p += sizeof(Key_t);
    parent = Find(Srctp, key0);

    sub_flags = *(int *)p;
    p += sizeof(int);
    Set_Sub_Flag(parent, sub_flags);
    key0 = KeyLshift(key0, Srctp->ndim);

    for(n=0; sub_flags; n++, sub_flags>>=1){
	if( (sub_flags&1) == 0 )
	    continue;
	key = KeyOrInt(key0, n);
	Msgf(("%d ", n));
	type = *(int *)p;
	p += sizeof(int);
	if( Sub_Flags_Type(type) ){
	    c = ChnAlloc(&Srctp->cellchn);
	    memcpy(c, p, Srctp->cell_sz);
	    p += Srctp->cell_sz;
	}else{
	    c = ChnAlloc(&Srctp->tbodychn);
	    memcpy(c, p, Srctp->tbody_sz);
	    p += Srctp->tbody_sz;
	}
	Enter(Srctp, key, c, type);
    }
    Msgf(("\n"));
    parent->type |= KIDSHERE;
    parent->type &= ~REQUESTED;
}    

