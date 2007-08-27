
/* The default I/O model here is that only proc 0 can read/write etc. */
/* There are two I/O modes for read and write, MULTI and SINGL */
/* All other functions do not have an I/O mode associated with them */
/* All calls must be loosely synchronous. There is a single file pointer */
/* This model works on all machines, but it may not be the most efficient */

/* MPMY_NFILE extends the mpmy_io model to include multiple file segments */

#include <stdarg.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include "Assert.h"
#include "protos.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "Msgs.h"
#include "bigmalloc.h"
#include "iozero.h"

#ifndef EINVAL
/* just in case... */
#define EINVAL 0
#endif

#define NFILES 256
static struct _File{
    int fd;
    int iomode;
    int iotype;
} _files[NFILES];

#define DEFAULT_PERMS 0644
static int do_nfileio;

int 
MPMY_Nfileio(int val){
    int oldval = do_nfileio;
    do_nfileio = val;
    return oldval;
}

MPMYFile *
MPMY_Fopen(const char *path, int mpmy_flags)
{
    int fd;
    int flags = 0;
    int iomode = MPMY_SINGL;
    int iotype = 0;

    Msgf(("Fopen %s\n", path));
    if (mpmy_flags & MPMY_RDONLY) flags |= O_RDONLY;
    if (mpmy_flags & MPMY_WRONLY) flags |= O_WRONLY;
    if (mpmy_flags & MPMY_RDWR)   flags |= O_RDWR;
    if (mpmy_flags & MPMY_APPEND) flags |= O_APPEND;
    if (mpmy_flags & MPMY_TRUNC)  flags |= O_TRUNC;
    if (mpmy_flags & MPMY_CREAT)  flags |= O_CREAT;

    if (mpmy_flags & MPMY_MULTI) iomode = MPMY_MULTI;
    if (mpmy_flags & MPMY_SINGL) iomode = MPMY_SINGL;

    /* This is a sub-mode which can work with either of the above */
    if (mpmy_flags & MPMY_NFILE) iotype = MPMY_NFILE;

    if (mpmy_flags & MPMY_INDEPENDENT) {
	iotype = MPMY_INDEPENDENT;
	iomode = MPMY_MULTI;
    }

    /* We need an external control for nfile mode */
    if (do_nfileio) iotype = MPMY_NFILE;

    if (iotype == MPMY_NFILE) {
	char real_path[256];
	sprintf(real_path, "%s.p%03d", path, MPMY_Procnum());
	fd = open(real_path, flags, DEFAULT_PERMS);
    } else if (iotype == MPMY_INDEPENDENT) {
	fd = open(path, flags, DEFAULT_PERMS);
    } else {
	if( strcmp(path, "-") != 0 ){
	    fd = open0(path, flags, DEFAULT_PERMS);
	}else{
	    /* Be very afraid. */
	    /* We made MPMY_RDONLY zero because O_RDONLY==0.  But why
	       did 'they' do that?! */
	    if ( (mpmy_flags & MPMY_RDWR) )
		Error("Can't open '-' RDWR\n");
#if MPMY_RDONLY != 0 
 # error Aargh.
#endif
	    if(mpmy_flags & MPMY_WRONLY )
		fd = 1;		/* stdout */
	    else		/* we can't test for MPMY_RDONLY! */
		fd = 0;		/* stdout */
	}
    }

    if (fd >= NFILES) Error("fd too large\n");
    if(fd >= 0){
	_files[fd].iomode = iomode;
	_files[fd].iotype = iotype;
	_files[fd].fd = fd;
	Msgf(("Fopen returns %d, iomode=%d, flags=0x%x\n", fd, iomode, flags));
	return &(_files[fd]);
    }else{
	Msgf(("Fopen fails, errno=%d\n", errno));
	return NULL;
    }
}

int
MPMY_Fclose(MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (fp->iotype & MPMY_NFILE || fp->iotype & MPMY_INDEPENDENT) {
	ret = close(fp->fd);
    } else {
	ret = close0(fp->fd);
    }
    fp->fd = -1;
    Msgf(("Close returns %d\n", ret));
    return ret;
}

int 
MPMY_Mkdir(const char *path, int mode){
    return mkdir0(path, mode);
}

int
MPMY_Fread(void *ptr, int size, int nitems, MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret = -1;

    Msgf(("Fread(ptr=%p, size=%d, nitems=%d, FILE=%p)\n",
	  ptr, size, nitems, Fp));
    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }

    switch (fp->iomode) {
      case MPMY_SINGL:
	ret = read0(fp->fd, ptr, size*nitems);
	break;
      case MPMY_MULTI:
	if (fp->iotype & MPMY_NFILE || fp->iotype & MPMY_INDEPENDENT)
	  ret = read(fp->fd, ptr, size*nitems);
	else
	  ret = read0_multi(fp->fd, ptr, size*nitems);
	break;
      default:
	ret =  -1;
	break;
    }

    Msgf(("Fread returns %d.\n", ret));
    if( ret < 0 ){
	Warning("MPMY_Fread:  read returns %d, errno=%d\n", ret, errno);
	return -1;
    }
    if (ret % size) 
      Shout("MPMY_Fread has a problem, ret=%d, size=%d\n", ret, size);
    return ret/size;
}

int
MPMY_Fwrite(const void *ptr, int size, int nitems, MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret = -1;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    
    switch (fp->iomode) {
      case MPMY_SINGL:
	ret = write0(fp->fd, ptr, size*nitems);
	break;
      case MPMY_MULTI:
	if (fp->iotype & MPMY_NFILE || fp->iotype & MPMY_INDEPENDENT)
	  ret = write(fp->fd, ptr, size*nitems);
	else
	  ret = write0_multi(fp->fd, ptr, size*nitems);
	break;
      default:
	ret = -1;
	Shout("Bad iomode in Fwrite\n");
    }

    Msgf(("Fwrite returns %d.\n", ret));
    if( ret < 0 ){
	Warning("MPMY_Fwrite: write returns %d, errno=%d\n", ret, errno);
	return -1;
    }
    if (ret % size) 
	Shout("MPMY_Fwrite has a problem: ret=%d, size=%d\n", ret, size);
    return ret/size;
}

int
MPMY_Fseek(MPMYFile *Fp, long offset, int whence)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;
    int real_whence = 0;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }

    if (whence == MPMY_SEEK_SET) real_whence = SEEK_SET;
    if (whence == MPMY_SEEK_CUR) real_whence = SEEK_CUR;
    if (whence == MPMY_SEEK_END) real_whence = SEEK_END;

    if (fp->iotype == MPMY_INDEPENDENT) {
	ret = lseek(fp->fd, offset, real_whence);
    } else {
	ret = lseek0(fp->fd, offset, real_whence);
    }
    if (ret != -1) ret = 0;
    Msgf(("Fseek returns %d\n", ret));
    return ret;
}

int
MPMY_Ftell(MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (fp->iotype == MPMY_INDEPENDENT) {
	ret = lseek(fp->fd, 0L, SEEK_CUR);
    } else {
	ret = tell0(fp->fd);
    }
    Msgf(("Ftell returns %d\n", ret));
    return ret;
}

int
MPMY_Flen(MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (fp->iotype & MPMY_NFILE) {
	ret = flen(fp->fd);
	Msgf(("flen segment %d\n", ret));
	MPMY_Combine(&ret, &ret, 1, MPMY_INT, MPMY_SUM);
    } else if (fp->iotype == MPMY_INDEPENDENT) {
	ret = flen(fp->fd);
    } else {
	ret = flen0(fp->fd);
    }
    Msgf(("Flen returns %d\n", ret));
    return ret;
}

#define MAXREQ 4

static int 
fseekrd_nfile(int fd, long offset, int whence, void *buf, int reclen,
	 int nrecs)
{
    int offend;
    void *tmpbuf = 0;		/* realloc null first time */
    int i;
    int nread;
    int size;
    struct {
	long offset;
	long reclen;
	long nrecs;
	long start;
	long flen;
    } parbuf, *allbuf;
    MPMY_Status stat;
    MPMY_Comm_request req[MAXREQ];	/* should be dynamic? */
    int nreq = 0;
    int procnum = MPMY_Procnum();
    int my_start, my_size, my_end;
    int nproc = MPMY_Nproc();

    if (whence != MPMY_SEEK_SET)
      Error("Fseekrd Nfiles must use SEEK_SET\n");

    parbuf.offset = offset;
    parbuf.reclen = reclen;
    parbuf.nrecs = nrecs;
    parbuf.flen = flen(fd);

    allbuf = Malloc(sizeof(parbuf)*nproc);
    
    MPMY_AllGather(&parbuf, sizeof(parbuf)/sizeof(long), MPMY_LONG, allbuf);

    my_start = 0;
    for (i = 0; i < nproc; i++) {
	allbuf[i].start = my_start;
	my_start += allbuf[i].flen;
    }
    offend = offset+nrecs*reclen;
    my_start = allbuf[procnum].start;
    my_size = allbuf[procnum].flen;
    my_end = my_start+my_size;

    Msgf(("offset %ld, nrecs*reclen %d\n", offset, nrecs*reclen));
    Msgf(("my_start %d, my_end %d\n", my_start, my_end));

    /* post Irecvs for my data */
    for (i = 0; i < nproc; i++) {
	if (allbuf[i].start<=offset && offset<allbuf[i].start+allbuf[i].flen) {
	    size = allbuf[i].start+allbuf[i].flen-offset;
	    if (size > reclen*nrecs) size = reclen*nrecs;
	    Msgf(("irecv %d from %d\n", size, i));
	    if (nreq >= MAXREQ) Error("Too man Irecvs in fseekrd\n");
	    MPMY_Irecv(buf, size, i, MPMY_IOTAG+nreq, &req[nreq]);
	    ++nreq;
	} else if (allbuf[i].start<offend && offend<=allbuf[i].start+allbuf[i].flen) {
	    size = offend-allbuf[i].start;
	    Msgf(("irecv %d from %d\n", size, i));
	    if (nreq >= MAXREQ) Error("Too man Irecvs in fseekrd\n");
	    MPMY_Irecv((char *)buf+reclen*nrecs-size, size, i, 
		       MPMY_IOTAG+nreq, &req[nreq]);
	    ++nreq;
	}
    }

    for (i = 0; i < nproc; i++) {
	offset = allbuf[i].offset;
	reclen = allbuf[i].reclen;
	nrecs = allbuf[i].nrecs;
	offend = offset+reclen*nrecs;

	if (my_start <= offset && offset < my_end) {
	    size = my_end-offset;
	    if (size > reclen*nrecs) size = reclen*nrecs;
	    lseek(fd, offset-my_start, SEEK_SET);
	    tmpbuf = Realloc(tmpbuf, size);
	    nread = read(fd, tmpbuf, size);
	    Msgf(("send %d to %d\n", size, i));
	    MPMY_send(tmpbuf, size, i, MPMY_IOTAG);
	} else if (my_start < offend && offend <= my_end) {
	    size = offend-my_start;
	    lseek(fd, 0, SEEK_SET);
	    tmpbuf = Realloc(tmpbuf, size);
	    nread = read(fd, tmpbuf, size);
	    Msgf(("send %d to %d\n", size, i));
	    MPMY_send(tmpbuf, size, i, MPMY_IOTAG+1);
	}
    }
    Free(tmpbuf);
    Free(allbuf);
    nread = 0;
#if 0
    for (i = 0; i < nreq; i++) {
	MPMY_Wait(req[i], &stat);
	Msgf(("Got %d from %d\n", stat.count, stat.src));
	nread += stat.count;
    }
#else
    i = 0;
    while( nreq ){
	int done;
	Msgf(("Testing req[%d] = %p\n", i, req[i]));
	MPMY_Test(req[i], &done, &stat);
	if( done ){
	    Msgf(("Got %d from %d\n", stat.count, stat.src));
	    nread += stat.count;
	    Msgf(("Moving req[%d] = req[%d] = ", i, nreq-1));
	    req[i] = req[--nreq];
	    Msgf(("%p\n", req[i]));
	}else{
	    Msgf(("Req[%d] not done yet\n", i));
	    i++;
	}
	assert(i <= nreq );
	if( i == nreq )
	    i = 0;
    }
#endif
    return nread/reclen;
}

/* This should be called from fseekrd0 instead of replicating code */

static int 
fseekrd(int fd, long offset, int whence, void *buf, int reclen,
	 int nrecs)
{
    int doseek;
    int real_whence;
    int nread = 0;
    int len;
    
    if( whence == SEEK_CUR ){
        doseek = (offset != 0);
    }else if( whence == SEEK_SET ){
        /* don't worry about errors.  If ftell returns -1, */
        /* doseek will be turned on, and the fseek below will */
        /* (probably) fail */
        doseek = (lseek(fd, 0L, SEEK_CUR) != offset);
    }else{
        doseek = 1;
    }

    if( doseek ){
	switch(whence){
	  case MPMY_SEEK_SET:
	    real_whence = SEEK_SET;
	    break;
	  case MPMY_SEEK_CUR:
	    real_whence = SEEK_CUR;
	    break;
	  case MPMY_SEEK_END:
	    real_whence = SEEK_END;
	    break;
	  default:
	    Shout("Illegal value of whence (%d) in fseekrd\n", whence);
	    return -1;
	}
	if (lseek(fd, offset, real_whence) == -1) {
	    Error("fseekrd: lseek(%d, %ld, %d) failed, errno=%d\n",
		  fd, offset, whence, errno);
	    return -1;
	}
    }
    len = 0;
    while (len < reclen*nrecs) {
	nread = read(fd, (char *)buf+len, reclen*nrecs-len);
	if (nread == -1) {
	    Error("fseekrd: read(%d, %d) failed, errno=%d\n", 
		  fd, reclen*nrecs-len, errno);
	    return -1;
	} else if (nread == 0) {
	    Error("fseekrd: read(%d, %d) got EOF\n", 
		  fd, reclen*nrecs-len);
	    return -1;
	} else {
	    printf("%d fseekrd(%d): got %d\n", fd, MPMY_Procnum(), nread);
	    len += nread;
	}
    }
    if (len != reclen*nrecs) Error("fseekrd: Wrong amount of data\n");
    return nread/reclen;
}

int 
MPMY_Fseekrd(MPMYFile *Fp, long offset, int whence, void *buf, int reclen,
	     int nrecs)
{
    struct _File *fp = (struct _File *)Fp;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (fp->iotype & MPMY_NFILE) {
	nrecs = fseekrd_nfile(fp->fd, offset, whence, buf, reclen, nrecs);
    } else if (fp->iotype == MPMY_INDEPENDENT) {
	nrecs = fseekrd(fp->fd, offset, whence, buf, reclen, nrecs);
    } else {
	nrecs = fseekrd0(fp->fd, offset, whence, buf, reclen, nrecs);
    }
    return nrecs;
}

#include "iozero.c"
#include "io_generic.c"
