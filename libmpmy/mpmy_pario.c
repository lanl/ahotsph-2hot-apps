/* This file contains the parallel I/O suitable for CMMD and NX.
   It wouldn't be hard to add cubix-syntax as well.  Are there any
   other options?  Would they fit in this structure?   The only
   difference between CMMD and NX is in Fopen, where one
   system calls gopen() and the other calls CMMD_set_io_mode().
   We use a pre-processor symbol (__INTEL_SSD__) or (__CM5__) to decide
   which one to use.  __INTEL_SSD__ is set by the ARCH-specific Makefiles,
   while __CM5__ is set by mpmy_cm5.c, which #includes this file.
*/
#include <stdarg.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include "protos.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "Msgs.h"
#include "iozero.h"

#ifndef EINVAL
/* just in case... */
#define EINVAL 0
#endif

#define NFILES 64
static struct _File{
    int fd;
    int iomode;
    int flags;
    char ungot;
    int ungot_active;
} _files[NFILES];

#ifdef __INTEL_SSD__
#if !defined( __PARAGON__ )
/* On delta and gamma, we have to provide our own 
   symbolic names for the setiomode modes */
#define M_UNIX 0
#define M_LOG 1
#define M_SYNC 2
#define M_RECORD 3
#define M_GLOBAL 4
extern void setiomode(int, int);
#endif /* !__PARAGON__ */

#if !defined( __PARAGON__ ) || defined( __SUNMOS__ )
/* We also have to define gopen.  Under sunmos, too. */
int gopen(const char *path, int flags, int mode, /* mode_t */ unsigned long int perms){
    int fd;

    if(flags & O_CREAT){
	if( MPMY_Procnum() == 0 ){
	    fd = open(path, flags, perms);
	    gsync();
	}else{
	    gsync();
	    fd = open(path, flags&~O_CREAT, perms);
	}
    }else{
	fd = open(path, flags, perms);
	gsync();
    }
    if( fd >= 0 ){
	setiomode(fd, mode);
    }else{
	Shout("open(\"%s\", %#x, 0%04o) failed, gopen returning -1, errno=%d\n", 
	      path, flags, perms, errno);
    }
    return fd;
}
#endif /* !__PARAGON__ || __SUNMOS__*/

MPMYFile *
MPMY_Fopen(const char *path, int flags)
{
    int fd;
    int mode = 0;
    int iomode = M_UNIX;	/* default */
    int real_flags = 0;

    Msgf(("Fopen %s\n", path));
    if (flags & MPMY_RDONLY) real_flags |= O_RDONLY;
    if (flags & MPMY_WRONLY) real_flags |= O_WRONLY;
    if (flags & MPMY_RDWR)   real_flags |= O_RDWR;
    if (flags & MPMY_APPEND) real_flags |= O_APPEND;
    if (flags & MPMY_TRUNC)  real_flags |= O_TRUNC;
    if (flags & MPMY_CREAT){
	real_flags |= O_CREAT;
	mode = 0644;
    }

    /* Should we make sure that only one of them is on?? */
    if (flags & MPMY_MULTI) iomode = M_SYNC;
    if (flags & MPMY_SINGL) iomode = M_GLOBAL;
    if (flags & MPMY_UNIX) iomode = M_UNIX;
    if (flags & MPMY_INDEPENDENT) iomode = M_UNIX;

    if (flags & MPMY_INDEPENDENT) fd = open(path, real_flags, mode);
    else if (flags & MPMY_IOZERO) fd = open0(path, real_flags, mode);
    else fd = gopen(path, real_flags, iomode, mode);

    if( fd >= 0 ){
	_files[fd].flags = flags;
	_files[fd].iomode = iomode;
	_files[fd].ungot_active = 0;
	Msgf(("Fopen returns %d\n", fd));
	_files[fd].fd = fd;
	MPMY_Sync();
	return &(_files[fd]);
    }else{
	Msgf(("Fopen returns NULL, errno=%d\n", errno));
	return NULL;
    }
}

#endif /* __INTEL_SSD__ */

#ifdef __CM5__
#include <cm/cmmd.h>
MPMYFile *
MPMY_Fopen(const char *path, int flags)
{
    int fd;
    int mode = 0644;
    CMMD_file_mode_t iomode = CMMD_independent;	/* default */
    int real_flags = 0;

    Msgf(("Fopen %s\n", path));
    if (flags & MPMY_RDONLY) real_flags |= O_RDONLY;
    if (flags & MPMY_WRONLY) real_flags |= O_WRONLY;
    if (flags & MPMY_RDWR)   real_flags |= O_RDWR;
    if (flags & MPMY_APPEND) real_flags |= O_APPEND;
    if (flags & MPMY_TRUNC)  real_flags |= O_TRUNC;
    if (flags & MPMY_CREAT)  real_flags |= O_CREAT;

    /* Should we make sure that only one of them is on?? */
    if (flags & MPMY_MULTI) iomode = CMMD_sync_seq;
    if (flags & MPMY_SINGL) iomode = CMMD_sync_bc;
    if (flags & MPMY_UNIX) iomode = CMMD_independent;
    if (flags & MPMY_IOZERO) iomode = CMMD_local;
    if (flags & MPMY_INDEPENDENT) iomode = CMMD_local;

    if (flags & MPMY_INDEPENDENT)
	fd = open(path, real_flags, mode);
    else if (flags & MPMY_IOZERO)
	fd = open0(path, real_flags, mode);
    else 
	fd = CMMD_global_open(path, real_flags, mode);

    if( fd >= 0 ){
	CMMD_set_io_mode(fd, iomode);
	_files[fd].flags = flags;
	_files[fd].iomode = iomode;
	_files[fd].ungot_active = 0;
	Msgf(("Fopen returns %d\n", fd));
	_files[fd].fd = fd;
	MPMY_Sync();
	return &(_files[fd]);
    }else{
	Msgf(("Fopen returns NULL, errno=%d\n", errno));
	return NULL;
    }
}

#endif /* __CM5__ */

int
MPMY_Fclose(MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (IoZero(fp)) ret = close0(fp->fd);
    else ret = close(fp->fd);
    Msgf(("Fclose of %d returns %d\n", fp->fd, ret));
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
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (IoZero(fp)) ret = read0(fp->fd, ptr, size*nitems);
    else ret = read(fp->fd, ptr, size*nitems);

    /* Msgf(("Fread from %d returns %d\n", fp->fd, ret)); */
    if (ret % size) Error("MPMY_Fread has a problem\n");
    return ret/size;
}

int
MPMY_Fwrite(void *ptr, int size, int nitems, MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (IoZero(fp)) ret= write0(fp->fd, ptr, size*nitems);
    else ret = write(fp->fd, ptr, size*nitems);
    /* Msgf(("Fwrite to %d returns %d.\n", fp->fd, ret)); */
    if (ret % size) Error("MPMY_Fwrite has a problem\n");
    return ret/size;
}

int
MPMY_Fseek(MPMYFile *Fp, long offset, int whence)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;
    int real_whence = 0;
    int iomode;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    iomode = fp->iomode;
    if (iomode == MPMY_MULTI)
      Error("Illegal iomode for Fseek\n");
    
    fp->ungot_active = 0;	/* read the manual for ungetc */
    if (whence == MPMY_SEEK_SET) real_whence = SEEK_SET;
    if (whence == MPMY_SEEK_CUR) real_whence = SEEK_CUR;
    if (whence == MPMY_SEEK_END) real_whence = SEEK_END;

    if (IoZero(fp)) ret = lseek0(fp->fd, offset, real_whence);
    else ret = lseek(fp->fd, offset, real_whence);
    if (ret != -1) ret = 0;
    Msgf(("Fseek of %ld on %d returns %d\n", offset, fp->fd, ret));
    return ret;
}

int
MPMY_Ftell(MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    int ret;
    int iomode;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    iomode = fp->iomode;
    if (iomode == MPMY_MULTI)
      Error("Illegal iomode for Ftell\n");

    if (IoZero(fp)) ret = tell0(fp->fd);
    else ret = lseek(fp->fd, 0L, SEEK_CUR);

    Msgf(("Ftell returns %d\n", ret));
    return ret;
}

int
MPMY_Flen(MPMYFile *Fp)
{
    struct _File *fp = (struct _File *)Fp;
    struct stat buf;
    int ret;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (IoZero(fp)) {
	ret = flen0(fp->fd);
    } else {
	/* the warning is because protos.h is included before sys/stat.h */
	/* it't too hard to fix. */
	fstat(fp->fd, &buf);
	ret = buf.st_size;
    }
    Msgf(("Flen returns %d\n", ret));
    return ret;
}

int 
MPMY_Fseekrd(MPMYFile *Fp, long offset, int whence, void *buf, int reclen,
	     int nrecs)
{
    struct _File *fp = (struct _File *)Fp;
    int doseek;

    if( fp == NULL ){
	errno = EINVAL;
	return -1;
    }
    if (IoZero(fp)) {
	nrecs = fseekrd0(fp->fd, offset, whence, buf, reclen, nrecs);
	return nrecs;
    }

    if( whence == MPMY_SEEK_CUR ){
	doseek = (offset != 0);
    }else if( whence == MPMY_SEEK_SET ){
	/* don't worry about errors.  If ftell returns -1, */
	/* doseek will be turned on, and the fseek below will */
	/* (probably) fail */
	doseek = (MPMY_Ftell(Fp) != offset);
    }else{
	doseek = 1;
    }

    if( doseek ){
	if( MPMY_Fseek(Fp, offset, whence) ){
	    if( whence == MPMY_SEEK_CUR && offset > 0 ){
		/* Make a final heroic effort to seek by reading forward! */
		char junk[BUFSIZ];
		int nleft = offset;
		while( nleft ){
		    int ntry = ( nleft > sizeof(junk) ) ? sizeof(junk) : nleft;
		    if( MPMY_Fread(junk, ntry, 1, Fp) != 1 ){
			Error("fseekrd: incremental fread(%#lx, %d, 1, %#lx) failed, errno=%d\n",
			      (unsigned long)junk, ntry, 
			      (unsigned long)fp, errno);
			return -1;
		    }
		    nleft -= ntry;
		}
	    }else{
		Error("fseekrd: fseek(%#lx, %ld, %d) failed, errno=%d\n",
		      (unsigned long)fp, offset, whence, errno);
		return -1;
	    }
	}
    }
    if( MPMY_Fread(buf, reclen, nrecs, Fp) != nrecs ){
	Error("fseekrd: fread(%#lx, %d, %d, %#lx) failed, errno=%d\n", 
	      (unsigned long)buf, reclen, nrecs, (unsigned long)fp, errno);
	return -1;
    }
    return nrecs;
}

#include "iozero.c"
#include "io_generic.c"
