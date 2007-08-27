#include <stdarg.h>
#include <fcntl.h>
#include <unistd.h>
/* Disregard warning in stdio.h.  It's because unistd.h #defines SEEK_CUR */
/* etc., and it's not the same as in stdio...  "Intel inside (tm)" */
#include <stdio.h>
#include <errno.h>
#include "protos.h"
#include "mpmy_io.h"
#include "Msgs.h"
#include "mpmy.h"
#include "bigmalloc.h"
#include "gc.h"

static unsigned long seekptr;		/* Is this a good idea? */

/* I'm assuming nobody would want an MPMY_Creat() */

int
MPMY_Global_Open(const char *path, int flags, ...)
{
    int fd;

    /* One must make sure that only one process creates or truncates a file */

    Msgf(("Global_Open %s\n", path));
    if (flags & O_CREAT) {
	va_list alist;
	int mode;

	va_start(alist, flags);
	mode = va_arg(alist, int);
	va_end(alist);
	if (MPMY_Procnum() == 0) {
	    fd = open(path, flags, mode);
	    MPMY_Sync();
	} else {
	    flags &= ~(O_CREAT|O_TRUNC);
	    MPMY_Sync();
	    fd = open(path, flags);
	}
    } else if (flags & O_TRUNC) {
	if (MPMY_Procnum() == 0) {
	    fd = open(path, flags);
	    MPMY_Sync();
	} else {
	    flags &= ~(O_TRUNC);
	    MPMY_Sync();
	    fd = open(path, flags);
	}
    } else {
	fd = open(path, flags);
    }
    seekptr = 0;
    Msgf(("Global_Open returns %d\n", fd));
    return fd;
}

int
MPMY_Global_Write(int fd, const void *buf, unsigned int nbytes)
{
    int ret;
    unsigned long seek_cur;
    int nproc = MPMY_Nproc();
    int procnum = MPMY_Procnum();
    int gcup, gcdown;

    gcup = Gcup(procnum, nproc);
    gcdown = Gcdown(procnum, nproc);

    Msgf(("Global_Write %d...", nbytes));
    
    if (procnum == 0) {
	Msgf(("seeklen=%lu\n", seekptr));
        if ((seek_cur = lseek(fd, seekptr, SEEK_CUR)) != seekptr)
            Error("lseek(seek_cur=%lu) failed, errno=%d\n", seek_cur, errno);
	ret = write(fd, buf, nbytes);
	(void) fsync(fd);
    } else {
	MPMY_Shift(gcdown, &seek_cur, sizeof(long), NULL, 0, NULL);
	Msgf(("seeklen=%lu\n", seek_cur));
	if (lseek(fd, seek_cur, SEEK_SET) != seek_cur)
	    Error("lseek(seek_cur=%lu) failed, errno=%d\n", seek_cur, errno);
	ret = write(fd, buf, nbytes);
	(void) fsync(fd);
    }

    seek_cur += ret;
    if (gcup != -1) {
	MPMY_Shift(gcup, NULL, 0, &seek_cur, sizeof(long), NULL);
#if 0
        MPMY_Recv(&seekptr, sizeof(long), MPMY_SOURCE_ANY, MPMY_TAG_ANY);
    } else {
	seekptr = seek_cur;
	MPMY_Bcast(&seekptr, 1, MPMY_UNSIGNED_LONG, procnum);
#endif
    }

    /* Should use MPMY_Bcast above but srcproc=0 enforced for now.  */
    MPMY_Combine(&seek_cur, &seekptr, 1, MPMY_UNSIGNED_LONG, MPMY_MAX);

    Msgf(("Global_Write returns %d...", ret));
    return ret;
}

int
MPMY_Global_Write0(int fd, const void *buf, unsigned int nbytes)
{
    int ret = 0;
    int allbytes;
    char *allbuf;

    Msgf(("Global_Write0 %d...", nbytes));
    
    allbytes = MPMY_NGather(buf, nbytes, MPMY_CHAR, (void **)&allbuf, 0);
    if (MPMY_Procnum() == 0) {
	ret = write(fd, allbuf, allbytes);
	Free(allbuf);
    } 
    Msgf(("Global_Write0 returns %d...", ret));
    return ret;
}

int
MPMY_Global_Read(int fd, void *buf, unsigned int nbytes)
{
    int ret;
    unsigned long seeklen = seekptr;
    int *sizes;
    int i;
    int nproc = MPMY_Nproc();
    int procnum = MPMY_Procnum();

    Msgf(("Global_Read %d...", nbytes));

    sizes = Malloc(sizeof(int)*nproc);
    MPMY_AllGather(&nbytes, 1, MPMY_INT, sizes);
    for (i = 0; i < procnum; i++)
      seeklen += sizes[i];

    Msgf(("seeklen=%lu\n", seeklen));
    if (lseek(fd, seeklen, SEEK_SET) != seeklen)
      Error("lseek(seeklen=%lu) failed, errno=%d\n", seeklen, errno);

    ret = read(fd, buf, nbytes);

    for (i = 0; i < nproc; i++)
      seekptr += sizes[i];

    Free(sizes);
    Msgf(("Global_Read returns %d...", ret));
    return ret;
}

int
MPMY_Global_Close(int fd)
{
    int ret;
    ret = close(fd);
    Msgf(("Global_Close returns %d\n", ret));
    return ret;
}
