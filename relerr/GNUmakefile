TREEHOME=..
treedir_sed=\.\.
# Application-code might want to use  something like instead:
#ifeq (X$(TREEHOME), X)
#TREEHOME=$(HOME)/nbody/tree16
#endif
##### Application-specific stuff goes here

src=	cofm_n.c	grav_nv.c	main_nv.c	physics_n.c \
	print.c 	grav_ring.c	# decomp.c

programname=relerr
treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/cofm_n$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/cofm_n$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/cofm_n$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/cofm_n$(objsuf): $(treedir)/include/pqsort.h physics_n.h $(treedir)/include/physics_generic.h
$(objdir)/cofm_n$(objsuf): $(treedir)/include/vop.h $(treedir)/include/Msgs.h $(treedir)/include/Assert.h
$(objdir)/cofm_n$(objsuf): $(treedir)/include/fastflpt.h $(treedir)/include/protos.h
$(objdir)/grav_nv$(objsuf): physics_n.h $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/grav_nv$(objsuf): $(treedir)/include/key.h $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h
$(objdir)/grav_nv$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/grav_nv$(objsuf): $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/grav_nv$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/grav_nv$(objsuf): $(treedir)/include/tensop.h $(treedir)/include/fastflpt.h $(treedir)/include/Msgs.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/fastflpt.h $(treedir)/include/Assert.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/SDF.h $(treedir)/include/protos.h $(treedir)/include/macr.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/malloc.h $(treedir)/include/bigmalloc.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/SDFwrite.h $(treedir)/include/SDFread.h physics_n.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/stk.h $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/singlio.h $(treedir)/include/mpmy.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/mpmy_io.h $(treedir)/include/mpmy_abnormal.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/gc.h $(treedir)/include/files.h $(treedir)/include/getparam.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/verify.h $(treedir)/include/randoms.h $(treedir)/include/ring.h
$(objdir)/main_nv$(objsuf): $(treedir)/include/decomp.h $(treedir)/include/image.h $(treedir)/include/memfile.h
$(objdir)/physics_n$(objsuf): physics_n.h $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/physics_n$(objsuf): $(treedir)/include/key.h $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h
$(objdir)/physics_n$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/physics_n$(objsuf): $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/physics_n$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/physics_generic.c
$(objdir)/physics_n$(objsuf): $(treedir)/include/protos.h $(treedir)/include/mpmy.h $(treedir)/include/vop.h
$(objdir)/physics_n$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/Assert.h $(treedir)/include/verify.h
$(objdir)/physics_n$(objsuf): $(treedir)/include/files.h $(treedir)/include/gc.h
$(objdir)/print$(objsuf): physics_n.h $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/print$(objsuf): $(treedir)/include/key.h $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h
$(objdir)/print$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/print$(objsuf): $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/print$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/protos.h
$(objdir)/print$(objsuf): $(treedir)/include/vop.h
$(objdir)/grav_ring$(objsuf): physics_n.h $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/grav_ring$(objsuf): $(treedir)/include/key.h $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h
$(objdir)/grav_ring$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/grav_ring$(objsuf): $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/grav_ring$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/grav_ring$(objsuf): $(treedir)/include/tensop.h $(treedir)/include/fastflpt.h
$(objdir)/grav_ring$(objsuf): $(treedir)/include/Msgs.h
