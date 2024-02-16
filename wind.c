#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "Assert.h"
#include "Malloc.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "Msgs.h"
#include "SDF.h"
#include "singlio.h"
#include "timers.h"
#include "verify.h"

extern Timer_t SDFreadTm;
/* I'm guessing this should be 'extern'? It's also defined in libSDF/SDFwrite.c. CIE */
extern Timer_t SDFwriteTm;
static char *header_buf;
static int header_size;
static int header_len;

/* How much to increase header buffer size when realloced */
#define BUF_INC 4096

/* How big is our line buffer */
#define LINE_LEN 512

/* Align the data segment on this size boundary. */
/* It MUST be less than LINE_LEN */
#define DATAALIGN 32

static void
outstr(const char *str)
     /* This is an obstack, but who's counting? */
{
    int len;

    len = strlen(str);
    if (header_size - header_len <= len) { /* <= deals with terminal null */
        header_size += BUF_INC + len;
        header_buf = Realloc(header_buf, header_size);
    }
    strcpy(header_buf + header_len, str);
    header_len += len;
}

void
SDFwritewind(const char *filename, int gnobj, int nobj, 
	     const void *btab, int windnobj, const void *windbtab, 
	     int bsize, int wsize, const char *winddesc, 
	     const char *bodydesc, ...){
    va_list alist;
    MPMYFile *myfd;
    int mode;
    int i, pad;
    char line[LINE_LEN];
    int ival;
    double dval;
    char *sval;
    char *name;
    const char *buf;
    int len;
    int ok, allok, retried;

    Msgf(("In Wtdata\n"));
    header_len = 0;

    header_size = BUF_INC;
    header_buf = Malloc(header_size);

    EnableTimer(&SDFwriteTm, "SDFwrite");
    StartTimer(&SDFwriteTm);
    va_start(alist, bodydesc);

    if (MPMY_Procnum() == 0) {
	outstr ("# SDF\n");
	sprintf(line, "parameter byteorder = 0x%x;\n", 
		SDFcpubyteorder()); outstr(line); 
	while( (name = va_arg(alist, char *)) ){
	    Msgf(("name(%lx)=%s\n", (unsigned long int)name, name));
	    switch( va_arg(alist, enum SDF_type_enum) ){
	      case SDF_INT:
		ival = va_arg(alist, int);
		sprintf(line, "int %s = %d;\n", name, ival); outstr(line);
		break;
	      case SDF_FLOAT:
		dval = va_arg(alist, double);
		sprintf(line, "float %s = %g;\n", name, dval); outstr(line);
		break;
	      case SDF_DOUBLE:
		dval = va_arg(alist, double);
		sprintf(line, "double %s = %.16g;\n", name, dval); 
		outstr(line);
		break;
	      case SDF_STRING:
		sval = va_arg(alist, char *);
		sprintf(line, "char %s[] = \"%s\";\n", name, sval);
		outstr(line);
		break;
	      default:
		Shout("Unexpected type in wtdata\n");
		break;
	    }
	}
	if( winddesc ){
	    outstr(winddesc);
	    if( windnobj > 0 )
	      sprintf(line, "[%d];\n", windnobj);
	    else
	      sprintf(line, "[];\n");
	    outstr(line);
	}
	if( bodydesc ){
	    outstr(bodydesc);
	    if( gnobj > 0 )
	      sprintf(line, "[%d];\n", gnobj);
	    else
	      sprintf(line, "[];\n");
	    outstr(line);
	}
	outstr("#\f\n");
	outstr ("# SDF-EOH ");
	/* This little bit of magic will cause the first word of data */
	/* to be aligned.  This isn't required by anything, but it makes */
	/* it a lot easier to use really primitive tools like od. */
	pad = (header_len+1)%DATAALIGN;	/* the +1 is to account for the '\n' */
	if( pad )
	  pad = DATAALIGN - pad;
	for(i=0; i<pad; i++){
	    line[i] = ' ';
	}
	line[pad] = '\n';
	line[pad+1] = '\0';
	outstr(line);
	/* Avoid separate write for header due to paragon limitations */
	/* It's only memory, after all */
	len = header_len+(wsize*windnobj)+bsize*nobj;
	buf = header_buf = Realloc(header_buf, len);
	memcpy(header_buf+header_len, windbtab, wsize*windnobj);
	memcpy(header_buf+header_len+(wsize*windnobj), btab, bsize*nobj);
    } else {
	len = bsize*nobj;
	buf = btab;
    }

    va_end(alist);

    /* The best number here is somewhat arbitrary, depends on nproc, */
    /* and should probably be externally setable */
    mode = MPMY_WRONLY|MPMY_CREAT|MPMY_TRUNC|MPMY_MULTI;
#if 0
    /* This is based on a faulty understanding of what MPMY_IOZERO does */
    if (gnobj > 0 && gnobj * bsize < 2000000) {
	mode |= MPMY_IOZERO;
    }else{
	mode |= MPMY_MULTI;
    }
#endif
    retried = 0;
 retry:
    myfd = MPMY_Fopen(filename, mode);
    if( myfd == NULL ){
	SeriousWarning("MPMY_Fopen(%s, 0x%x) returns NULL, errno=%d\n",
		       filename, mode, errno);
	goto outahere;
    }

    i = MPMY_Fwrite(buf, 1, len, myfd);
    if (i != len){
	SeriousWarning("MPMY_Fwrite(btab, len=%d) only wrote %d, errno=%d\n", 
	      len, i, errno);
	SeriousWarning("\"%s\" is probably corrupt!\n", filename);
    }
    ok = (i==len);
    /* Should there be an MPMY_LAND and MPMY_LOR ? */
    allok = 0;
    MPMY_Combine(&ok, &allok, 1, MPMY_INT, MPMY_BAND);
    Msgf(("ok=%d, allok=%d\n", ok, allok));
    if( !allok ){
	int retryable;
#ifdef ETIMEDOUT  /* This seems to be a solaris thing */
	retryable = !retried && (ok || errno==ETIMEDOUT);
#else
	retryable = 0;
#endif
	Msgf(("retryable (local) = %d\n", retryable));
	MPMY_Combine(&retryable, &retryable, 1, MPMY_INT, MPMY_BAND);
	Msgf(("retryable (global) = %d\n", retryable));
	if( retryable ){
	    SinglShout("Fingers crossed, we are going to retry!\n");
	    retried = 1;		/* only retry once! */
	    MPMY_Fclose(myfd);
#ifdef ETIMEDOUT
	    /* The failure mode we are trying to recover from is an
	       NFS timeout which might go away in a few seconds.  We
	       trust that if ETIMEDOUT exists, then sleep does too...  */
	    SinglShout("sleep(30), maybe the timeout will go away!\n");
	    sleep(30);
#endif
	    goto retry;
	}
    }
    MPMY_Fclose(myfd);
 outahere:
    Free(header_buf);
    header_size = header_len = 0;
    header_buf = NULL;
    Msgf(("SDFwrite2 done\n"));

    StopTimer(&SDFwriteTm);
    OutputTimer(&SDFwriteTm, singlPrintf); /* global sync and set timer->max */
    if (SDFwriteTm.max != 0.0) 
      singlPrintf("write speed %.0f kb/s\n", 
		  gnobj*bsize/(1000.0*SDFwriteTm.max));
    DisableTimer(&SDFwriteTm);	/* suppress printing again in OutputTimers */
}

#define MAXNAMES 64

SDF *SDFreadwind(char *name, void **btabp, int *gnobjp, int *nobjp, 
		 int stride,
		 /* char *name, offset_t offset, int *confirm */...)
{
    va_list ap;
    int start;
    SDF *sdfp;
    int gnobj, nobj;
    void *btab;
    void *addrs[MAXNAMES];
    char *names[MAXNAMES];
    int strides[MAXNAMES];
    int nobjs[MAXNAMES];
    int64_t starts[MAXNAMES];
    int *confirm;
    int nnames;

    EnableTimer(&SDFreadTm, "SDFread");
    StartTimer(&SDFreadTm);

    VerifySX(sdfp = SDFopen(0, name),SinglShout("%s", SDFerrstring));
    
    if( SDFgetint(sdfp, "npart", &gnobj) ){
	/* Hopefully calling va_start and va_end in here won't disturb */
	/* the real loop over arguments below... */
	va_start(ap, stride);
	names[0] = va_arg(ap, char *);
	gnobj = SDFnrecs(names[0], sdfp);
	va_end(ap);
	if( MPMY_Procnum() == 0 ){
	    SinglShout("%s does not have an \"npart\".\n", name);
	    SinglShout("Guessing npart=%d from SDFnrecs(., %s)\n", 
		       gnobj, names[0]);
	}
    }
    
/*      NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start); */
    nobj = gnobj;
    start = 0;

    btab = Calloc(nobj, stride);
    Msgf(("Proc %d starting at %d in file, reading %d of %d\n",
	  MPMY_Procnum(), start, nobj, gnobj));

    nnames = 0;
    va_start(ap, stride);
    while(( names[nnames] = va_arg(ap, char *)) != NULL ){
	assert(nnames < MAXNAMES);
	addrs[nnames] = va_arg(ap, int) + (char *)btab;
	confirm = va_arg(ap, int *);
	if( !SDFhasname(names[nnames], sdfp) ){
	    *confirm = 0;
	    Msgf(("SDF file does not have %s\n", names[nnames]));
	    continue;
	}else{
	    *confirm = 1;
	}
	starts[nnames] = start;
	nobjs[nnames] = nobj;
	strides[nnames] = stride;
	nnames++;
    }
    va_end(ap);
    
    VerifyX(0==SDFseekrdvecsarr(sdfp, nnames,
			   names, starts, nobjs, addrs, strides),
	    Shout("%s", SDFerrstring));

    /* Don't worry about this too much for now, but be aware that gnobj 
       doesn't actually indicate how many wind sources there are in the 
       simulation.  I probably don't need to do an MPMY_Combine here.  */
    *nobjp = nobj;
    *gnobjp = nobj;
    MPMY_Combine(nobjp, gnobjp, 1, MPMY_INT, MPMY_SUM);
    Msgf(("nobj=%d, gnobj=%d\n", nobj, gnobj));

    *btabp = btab;
    StopTimer(&SDFreadTm);
    OutputTimer(&SDFreadTm, singlPrintf); /* global sync and sets timer->max */
    singlPrintf("read speed %.0f kb/s\n", gnobj*nnames*sizeof(float)/(1000.0*SDFreadTm.max));
    DisableTimer(&SDFreadTm);
    return sdfp;
}
