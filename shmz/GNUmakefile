TREEHOME=..
treedir_sed=\.\.
##### Application-specific stuff goes here

ifndef prog
src = main.c shmz_ring.c cofm.c physics.c print.c f77malloc.c \
	walk.c tree.c extra.c
override PRELIBS:=$(PRELIBS) libzft.a libutil.a liblocalbcslib.a
programname=shmz
override LD:=$(F77)
else
src=$(prog).c
programname=$(prog)
endif

#$(programname).$(ARCH)-p: libzft-p.a libutil-p.a liblocalbcslib-p.a
#$(programname).$(ARCH): libzft.a libutil.a liblocalbcslib.a

#serial: $(programname).$(ARCH)

treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): physics.h ndim.h $(treedir)/include/tree.h
$(objdir)/main$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/main$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/main$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h complex.h
$(objdir)/main$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/singlio.h
$(objdir)/main$(objsuf): $(treedir)/include/getparam.h $(treedir)/include/mpmy.h $(treedir)/include/SDF.h
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/SDFread.h $(treedir)/include/SDFwrite.h $(treedir)/include/ring.h
$(objdir)/main$(objsuf): $(treedir)/include/vop.h $(treedir)/include/Msgs.h $(treedir)/include/malloc.h
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf): physics.h ndim.h $(treedir)/include/tree.h
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf):
$(objdir)/shmz_ring$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/shmz_ring$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/shmz_ring$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h complex.h
$(objdir)/shmz_ring$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/shmz_ring$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/cofm$(objsuf): $(treedir)/include/key.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h
$(objdir)/cofm$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics.h ndim.h
$(objdir)/cofm$(objsuf): complex.h $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/cofm$(objsuf): fmm.h
$(objdir)/physics$(objsuf): physics.h ndim.h $(treedir)/include/tree.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/stk.h
$(objdir)/physics$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/physics$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/physics$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h complex.h
$(objdir)/physics$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/physics_generic.c
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/protos.h
$(objdir)/physics$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/vop.h $(treedir)/include/Msgs.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/verify.h $(treedir)/include/files.h $(treedir)/include/gc.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): physics.h ndim.h $(treedir)/include/tree.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/stk.h
$(objdir)/print$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/print$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/print$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h complex.h
$(objdir)/print$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/protos.h
$(objdir)/print$(objsuf): $(treedir)/include/vop.h
$(objdir)/f77malloc$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/f77malloc$(objsuf):
$(objdir)/f77malloc$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
