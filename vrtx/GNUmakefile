##### Application-specific stuff goes here

src = main_vrtx.c cofm_vrtx.c mac.c physics_vrtx.c relaxomega.c omegat.c \
  mono_vrtx.c quads_vrtx.c print.c diags.c remesh.c neigh.c fixomegatot.c \
  $(integrate)

## Note: $(integrate) does not get included in "make depends"
## This will cause the code to fail if you change a header file

treedir=$(TREEDHOME)

ifndef programname
vrtx vrtx_wing vrtx_helic vrtx_euler vrtx_rk : 
	$(MAKE) programname=$@
endif


# The generic version uses the rk integrator.  In order for this all to
# work properly, I think we will need separate 'main's for the different
# programs, so they can use slightly different calling conventions for
# Update...
ifeq "$(programname)" "vrtx"
       integrate=integrate_vrtx_rk.c
endif

ifeq "$(programname)" "vrtx_wing"
	integrate=integrate_vrtx_wing.c
endif

ifeq "$(programname)" "vrtx_helic"
	integrate=integrate_vrtx_helic.c
endif

ifeq "$(programname)" "vrtx_euler"
	integrate=integrate_vrtx_euler.c
endif

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

$(objdir)/quads_vrtx$(objsuf) : quads_vrtx.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c quads_vrtx.c
	mv quads_vrtx$(objsuf) $(objdir)

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/main_vrtx$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/timers.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/key.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/stk.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/error.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/chn.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/SDF.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/SDFwrite.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/Assert.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/protos.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/vop.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/singlio.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/files.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/SDFread.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/getparam.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/mpmy.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/malloc.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/memfile.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/abm.h
$(objdir)/main_vrtx$(objsuf): $(treedir)/include/dll.h neigh.h
$(objdir)/cofm_vrtx$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/timers.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/key.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/stk.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/error.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/chn.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/Assert.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/vop.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/cofm_vrtx$(objsuf): $(treedir)/include/protos.h
$(objdir)/mac$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/mac$(objsuf): $(treedir)/include/timers.h
$(objdir)/mac$(objsuf): $(treedir)/include/key.h
$(objdir)/mac$(objsuf): $(treedir)/include/stk.h
$(objdir)/mac$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/mac$(objsuf): $(treedir)/include/error.h
$(objdir)/mac$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/mac$(objsuf): $(treedir)/include/chn.h
$(objdir)/mac$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/mac$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/mac$(objsuf): $(treedir)/include/vop.h
$(objdir)/physics_vrtx$(objsuf): physics_vrtx.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/tree.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/timers.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/key.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/stk.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/error.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/chn.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/pqsort.h
$(objdir)/physics_vrtx$(objsuf): structures.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/physics_generic.c
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/protos.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/mpmy.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/vop.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/Assert.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/verify.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/files.h
$(objdir)/physics_vrtx$(objsuf): $(treedir)/include/gc.h
$(objdir)/relaxomega$(objsuf): physics_vrtx.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/tree.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/timers.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/key.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/stk.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/error.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/chn.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/pqsort.h
$(objdir)/relaxomega$(objsuf): structures.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/relaxomega$(objsuf): $(treedir)/include/vop.h
$(objdir)/omegat$(objsuf): $(treedir)/include/protos.h
$(objdir)/omegat$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/omegat$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/omegat$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/omegat$(objsuf): $(treedir)/include/timers.h
$(objdir)/omegat$(objsuf): $(treedir)/include/key.h
$(objdir)/omegat$(objsuf): $(treedir)/include/stk.h
$(objdir)/omegat$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/omegat$(objsuf): $(treedir)/include/error.h
$(objdir)/omegat$(objsuf): $(treedir)/include/chn.h
$(objdir)/omegat$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/omegat$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/omegat$(objsuf): $(treedir)/include/vop.h
$(objdir)/omegat$(objsuf): $(treedir)/include/fastflpt.h neigh.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/mono_vrtx$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/timers.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/key.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/stk.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/error.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/chn.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/vop.h
$(objdir)/mono_vrtx$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/quads_vrtx$(objsuf): physics_vrtx.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/tree.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/timers.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/key.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/stk.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/error.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/chn.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/pqsort.h
$(objdir)/quads_vrtx$(objsuf): structures.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/protos.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/vop.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/Assert.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/quads_vrtx$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/print$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/print$(objsuf): $(treedir)/include/timers.h
$(objdir)/print$(objsuf): $(treedir)/include/key.h
$(objdir)/print$(objsuf): $(treedir)/include/stk.h
$(objdir)/print$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/print$(objsuf): $(treedir)/include/error.h
$(objdir)/print$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/print$(objsuf): $(treedir)/include/chn.h
$(objdir)/print$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/print$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/print$(objsuf): $(treedir)/include/protos.h
$(objdir)/print$(objsuf): $(treedir)/include/vop.h
$(objdir)/diags$(objsuf): $(treedir)/include/mpmy.h
$(objdir)/diags$(objsuf): $(treedir)/include/timers.h physics_vrtx.h
$(objdir)/diags$(objsuf): $(treedir)/include/tree.h
$(objdir)/diags$(objsuf): $(treedir)/include/key.h
$(objdir)/diags$(objsuf): $(treedir)/include/stk.h
$(objdir)/diags$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/diags$(objsuf): $(treedir)/include/error.h
$(objdir)/diags$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/diags$(objsuf): $(treedir)/include/chn.h
$(objdir)/diags$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/diags$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/diags$(objsuf): $(treedir)/include/vop.h
$(objdir)/diags$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/diags$(objsuf): $(treedir)/include/singlio.h
$(objdir)/diags$(objsuf): $(treedir)/include/protos.h
$(objdir)/remesh$(objsuf): $(treedir)/include/tree.h
$(objdir)/remesh$(objsuf): $(treedir)/include/timers.h
$(objdir)/remesh$(objsuf): $(treedir)/include/key.h
$(objdir)/remesh$(objsuf): $(treedir)/include/stk.h
$(objdir)/remesh$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/remesh$(objsuf): $(treedir)/include/error.h
$(objdir)/remesh$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/remesh$(objsuf): $(treedir)/include/chn.h
$(objdir)/remesh$(objsuf): $(treedir)/include/pqsort.h
$(objdir)/remesh$(objsuf): $(treedir)/include/Msgs.h physics_vrtx.h
$(objdir)/remesh$(objsuf): structures.h
$(objdir)/remesh$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/remesh$(objsuf): $(treedir)/include/vop.h
$(objdir)/remesh$(objsuf): $(treedir)/include/Assert.h
$(objdir)/neigh$(objsuf): physics_vrtx.h $(treedir)/include/tree.h
$(objdir)/neigh$(objsuf): $(treedir)/include/timers.h
$(objdir)/neigh$(objsuf): $(treedir)/include/key.h
$(objdir)/neigh$(objsuf): $(treedir)/include/stk.h
$(objdir)/neigh$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/neigh$(objsuf): $(treedir)/include/error.h
$(objdir)/neigh$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/neigh$(objsuf): $(treedir)/include/chn.h
$(objdir)/neigh$(objsuf): $(treedir)/include/pqsort.h structures.h
$(objdir)/neigh$(objsuf): $(treedir)/include/physics_generic.h neigh.h
$(objdir)/neigh$(objsuf): $(treedir)/include/vop.h
$(objdir)/neigh$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/fixomegatot$(objsuf): physics_vrtx.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/tree.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/timers.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/key.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/stk.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/bigmalloc.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/error.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/gccextensions.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/chn.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/pqsort.h
$(objdir)/fixomegatot$(objsuf): structures.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/physics_generic.h
$(objdir)/fixomegatot$(objsuf): $(treedir)/include/vop.h
