#ifndef IOzeroDOTh
#define IOzeroDOTh

#define MPMY_IOTAG (0x183)
#define IoZero(fp) (fp->flags & MPMY_IOZERO)

static int open0(const char *path, int flags, int mode);
static int close0(int fd);
static int mkdir0(const char *path, int mode);
static int read0(int fd, void *buf, unsigned int nbytes);
static int write0(int fd, const void *buf, unsigned int nbytes);
static int lseek0(int fd, long offset, int whence);
static int tell0(int fd);
static int flen(int fd);
static int flen0(int fd);
static int  fseekrd0(int fd, long offset, int whence, void *buf, 
		     int reclen, int nrecs);
static int write0_multi(int fd, const void *buf, unsigned int nbytes);
static int read0_multi(int fd, void *buf, unsigned int nbytes);
#endif
