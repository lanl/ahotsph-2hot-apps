# Make.$(ARCH) sets many of the variables used below including:
# CC, CFLAGS, AS, RANLIB, objdir, objsuf, asmdir

treedir=..
treedir_sed=\.\.
appexcludes=
libname=libsw

# The "Salmon Warren Utility library"
src:= \
bigmalloc.c 	key.c \
Msgs.c		error.c 	files.c \
chn.c		dll.c \
stk.c		heap.c		randoms.c \
malloc.c	gc.c		byteswap.c \
timers.c	SDFwrite.c	SDFread.c \
ring.c		obstack.c \
mpmy_combine.c	counters.c \
mpmy_gather.c	memfile.c \
msgdirinit.c	singlio.c	abm.c \
batch.c		poll.c	\
mpi_bcast.c	mpi_reduce.c	keycvt.c \
peano.c		SDFreadf.c
 
include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/Malloc$(objsuf):
$(objdir)/Malloc$(objsuf):
$(objdir)/Malloc$(objsuf):
$(objdir)/Malloc$(objsuf):
$(objdir)/Malloc$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/Malloc$(objsuf):
$(objdir)/Malloc$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/bigmalloc.h
$(objdir)/Malloc$(objsuf): $(treedir)/include/error.h $(treedir)/include/malloc.h
$(objdir)/key$(objsuf):
$(objdir)/key$(objsuf):
$(objdir)/key$(objsuf): $(treedir)/include/key.h
$(objdir)/key$(objsuf):
$(objdir)/key$(objsuf): $(treedir)/include/protos.h
$(objdir)/Msgs$(objsuf):
$(objdir)/Msgs$(objsuf):
$(objdir)/Msgs$(objsuf):
$(objdir)/Msgs$(objsuf):
$(objdir)/Msgs$(objsuf):
$(objdir)/Msgs$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/Msgs$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/Msgs$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/Msgs.h
$(objdir)/Msgs$(objsuf): $(treedir)/include/protos.h
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf):
$(objdir)/error$(objsuf): $(treedir)/include/error.h
$(objdir)/error$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/Msgs.h
$(objdir)/error$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/error$(objsuf): $(treedir)/include/mpmy_abnormal.h $(treedir)/include/protos.h
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf):
$(objdir)/files$(objsuf): $(treedir)/include/protos.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/files$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/mpmy.h
$(objdir)/files$(objsuf): $(treedir)/include/timers.h
$(objdir)/chn$(objsuf):
$(objdir)/chn$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/chn$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h $(treedir)/include/Msgs.h
$(objdir)/chn$(objsuf):
$(objdir)/dll$(objsuf): $(treedir)/include/chn.h
$(objdir)/dll$(objsuf):
$(objdir)/dll$(objsuf): $(treedir)/include/dll.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/dll$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/stk$(objsuf): $(treedir)/include/stk.h
$(objdir)/stk$(objsuf):
$(objdir)/stk$(objsuf):
$(objdir)/stk$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/stk$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/heap$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/heap$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/heap$(objsuf):
$(objdir)/heap$(objsuf):
$(objdir)/heap$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/heap.h
$(objdir)/randoms$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/randoms$(objsuf): $(treedir)/include/randoms.h
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf): $(treedir)/include/malloc.h $(treedir)/include/protos.h $(treedir)/include/Msgs.h
$(objdir)/malloc$(objsuf):
$(objdir)/malloc$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/error.h
$(objdir)/gc$(objsuf): $(treedir)/include/gc.h
$(objdir)/byteswap$(objsuf):
$(objdir)/byteswap$(objsuf):
$(objdir)/byteswap$(objsuf):
$(objdir)/byteswap$(objsuf):
$(objdir)/byteswap$(objsuf):
$(objdir)/byteswap$(objsuf): $(treedir)/include/byteswap.h
$(objdir)/timers$(objsuf):
$(objdir)/timers$(objsuf):
$(objdir)/timers$(objsuf): $(treedir)/include/timers.h $(treedir)/include/mpmy.h $(treedir)/include/mpmy_time.h
$(objdir)/timers$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/timers$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf):
$(objdir)/SDFwrite$(objsuf): $(treedir)/include/SDF.h
$(objdir)/SDFwrite$(objsuf): $(treedir)/include/SDFwrite.h $(treedir)/include/bigmalloc.h
$(objdir)/SDFwrite$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/SDFwrite$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/mpmy_io.h
$(objdir)/SDFwrite$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/protos.h $(treedir)/include/singlio.h
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf):
$(objdir)/SDFread$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/SDF.h
$(objdir)/SDFread$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/SDFread$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/bigmalloc.h
$(objdir)/SDFread$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/verify.h $(treedir)/include/SDFread.h
$(objdir)/SDFread$(objsuf): $(treedir)/include/gc.h $(treedir)/include/singlio.h
$(objdir)/ring$(objsuf): $(treedir)/include/ring.h $(treedir)/include/bigmalloc.h
$(objdir)/ring$(objsuf):
$(objdir)/ring$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/ring$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/gc.h
$(objdir)/ring$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/ring$(objsuf):
$(objdir)/obstack$(objsuf): $(treedir)/include/obstack.h
$(objdir)/obstack$(objsuf):
$(objdir)/obstack$(objsuf):
$(objdir)/obstack$(objsuf):
$(objdir)/obstack$(objsuf):
$(objdir)/mpmy_combine$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/mpmy_combine$(objsuf):
$(objdir)/mpmy_combine$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/mpmy_combine$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/stk.h
$(objdir)/mpmy_combine$(objsuf):
$(objdir)/mpmy_combine$(objsuf): $(treedir)/include/gc.h op_template.c
$(objdir)/counters$(objsuf):
$(objdir)/counters$(objsuf):
$(objdir)/counters$(objsuf): $(treedir)/include/timers.h $(treedir)/include/mpmy.h $(treedir)/include/Assert.h
$(objdir)/counters$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/mpmy_gather$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/Msgs.h
$(objdir)/mpmy_gather$(objsuf):
$(objdir)/mpmy_gather$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/bigmalloc.h
$(objdir)/mpmy_gather$(objsuf):
$(objdir)/mpmy_gather$(objsuf): $(treedir)/include/error.h $(treedir)/include/gc.h $(treedir)/include/verify.h
$(objdir)/memfile$(objsuf):
$(objdir)/memfile$(objsuf):
$(objdir)/memfile$(objsuf):
$(objdir)/memfile$(objsuf):
$(objdir)/memfile$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/memfile$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/Msgs.h
$(objdir)/memfile$(objsuf):
$(objdir)/memfile$(objsuf): $(treedir)/include/protos.h $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf):
$(objdir)/msgdirinit$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/msgdirinit$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/msgdirinit$(objsuf): $(treedir)/include/mpmy_abnormal.h $(treedir)/include/protos.h
$(objdir)/msgdirinit$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/msgdirinit$(objsuf): $(treedir)/include/files.h $(treedir)/include/mpmy_io.h
$(objdir)/singlio$(objsuf):
$(objdir)/singlio$(objsuf):
$(objdir)/singlio$(objsuf):
$(objdir)/singlio$(objsuf): $(treedir)/include/protos.h
$(objdir)/singlio$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/Msgs.h
$(objdir)/singlio$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/singlio.h
$(objdir)/abm$(objsuf):
$(objdir)/abm$(objsuf):
$(objdir)/abm$(objsuf):
$(objdir)/abm$(objsuf):
$(objdir)/abm$(objsuf):
$(objdir)/abm$(objsuf): $(treedir)/include/abm.h $(treedir)/include/chn.h
$(objdir)/abm$(objsuf): $(treedir)/include/dll.h $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/abm$(objsuf): $(treedir)/include/protos.h $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/abm$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/bigmalloc.h
$(objdir)/abm$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/abm$(objsuf):
$(objdir)/abm$(objsuf): $(treedir)/include/verify.h $(treedir)/include/gc.h
$(objdir)/batch$(objsuf): $(treedir)/include/batch.h $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/batch$(objsuf): $(treedir)/include/stk.h
$(objdir)/batch$(objsuf):
$(objdir)/batch$(objsuf):
$(objdir)/batch$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/batch$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/batch$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/batch$(objsuf):
$(objdir)/poll$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/Msgs.h
$(objdir)/poll$(objsuf):
$(objdir)/poll$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/poll$(objsuf):
$(objdir)/poll$(objsuf): $(treedir)/include/gc.h $(treedir)/include/error.h $(treedir)/include/dll.h
$(objdir)/poll$(objsuf): $(treedir)/include/chn.h
$(objdir)/poll$(objsuf):
$(objdir)/mpi_bcast$(objsuf): $(treedir)/include/swampi.h $(treedir)/include/Msgs.h
$(objdir)/mpi_bcast$(objsuf):
$(objdir)/mpi_bcast$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/gc.h
$(objdir)/mpi_bcast$(objsuf): $(treedir)/include/error.h
$(objdir)/mpi_reduce$(objsuf):
$(objdir)/mpi_reduce$(objsuf):
$(objdir)/mpi_reduce$(objsuf):
$(objdir)/mpi_reduce$(objsuf):
$(objdir)/mpi_reduce$(objsuf):
$(objdir)/mpi_reduce$(objsuf): $(treedir)/include/swampi.h $(treedir)/include/stk.h
$(objdir)/mpi_reduce$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/mpi_reduce$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/mpi_reduce$(objsuf): $(treedir)/include/gc.h $(treedir)/include/Msgs.h
$(objdir)/mpi_reduce$(objsuf):
$(objdir)/mpi_reduce$(objsuf): mpi_template.c
$(objdir)/keycvt$(objsuf):
$(objdir)/keycvt$(objsuf):
$(objdir)/keycvt$(objsuf): $(treedir)/include/protos.h $(treedir)/include/mpmy.h $(treedir)/include/timers.h
$(objdir)/keycvt$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/keycvt$(objsuf):
$(objdir)/keycvt$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/Assert.h
$(objdir)/keycvt$(objsuf): $(treedir)/include/error.h $(treedir)/include/key.h
$(objdir)/keycvt$(objsuf):
$(objdir)/keycvt$(objsuf):
$(objdir)/keycvt$(objsuf): $(treedir)/include/peano.h
$(objdir)/keycvt$(objsuf): $(treedir)/include/keycvt.h
$(objdir)/peano$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/peano$(objsuf):
$(objdir)/peano$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/Assert.h
$(objdir)/peano$(objsuf): $(treedir)/include/error.h $(treedir)/include/key.h
$(objdir)/peano$(objsuf):
$(objdir)/peano$(objsuf):
$(objdir)/peano$(objsuf): $(treedir)/include/vop.h
