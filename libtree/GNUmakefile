# Make.$(ARCH) sets many of the variables used below including:
# CC, CFLAGS, AS, RANLIB, objdir, objsuf, asmdir

treedir=..
treedir_sed=\.\.
appexcludes=
libname=libtree

src:= \
tree.c pqsort.c traverse.c walk.c decomp.c image.c

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/tree$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/tree$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/verify.h
$(objdir)/tree$(objsuf): $(treedir)/include/key.h $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/tree$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/chn.h
$(objdir)/tree$(objsuf): $(treedir)/include/pqsort.h $(treedir)/include/gc.h $(treedir)/include/mpmy.h
$(objdir)/tree$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/protos.h
$(objdir)/pqsort$(objsuf): $(treedir)/include/pqsort.h $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/pqsort$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/pqsort$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/bigmalloc.h
$(objdir)/pqsort$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/gc.h $(treedir)/include/Msgs.h
$(objdir)/pqsort$(objsuf): $(treedir)/include/stk.h $(treedir)/include/decomp.h $(treedir)/include/abm.h
$(objdir)/pqsort$(objsuf): $(treedir)/include/chn.h $(treedir)/include/dll.h
$(objdir)/traverse$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/traverse$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/bigmalloc.h
$(objdir)/traverse$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/traverse$(objsuf): $(treedir)/include/stk.h $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/traverse$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/walk$(objsuf): $(treedir)/include/abm.h $(treedir)/include/chn.h $(treedir)/include/dll.h
$(objdir)/walk$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/timers.h $(treedir)/include/tree.h
$(objdir)/walk$(objsuf): $(treedir)/include/key.h $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h
$(objdir)/walk$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/walk$(objsuf): $(treedir)/include/pqsort.h $(treedir)/include/Msgs.h $(treedir)/include/Assert.h
$(objdir)/decomp$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/decomp$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/key.h
$(objdir)/decomp$(objsuf): $(treedir)/include/decomp.h $(treedir)/include/pqsort.h $(treedir)/include/timers.h
$(objdir)/decomp$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/mpmy.h $(treedir)/include/stk.h
$(objdir)/decomp$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/gc.h $(treedir)/include/protos.h
$(objdir)/image$(objsuf): $(treedir)/include/SDF.h $(treedir)/include/protos.h $(treedir)/include/malloc.h
$(objdir)/image$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/image$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/SDFwrite.h
$(objdir)/image$(objsuf): $(treedir)/include/SDFread.h $(treedir)/include/pixel.h $(treedir)/include/tree.h
$(objdir)/image$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h $(treedir)/include/stk.h
$(objdir)/image$(objsuf): $(treedir)/include/chn.h $(treedir)/include/pqsort.h $(treedir)/include/vop.h
$(objdir)/image$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/files.h $(treedir)/include/getparam.h
$(objdir)/image$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/randoms.h $(treedir)/include/fastflpt.h
$(objdir)/image$(objsuf): $(treedir)/include/singlio.h
