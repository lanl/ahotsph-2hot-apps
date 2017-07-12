TREEHOME=..
treedir_sed=\.\.
##### Application-specific stuff goes here

ifndef prog
src=cofm.c eos.c grav.c mac.c main.c physics.c print.c integrate.c \
	sphcofm.c sph.c sphinit.c shrink.c physics_sph.c sphprint.c sphplus.c \
	SDFreadf.c nrutil.c coolandburn.c params.c

#fsrc=build.f solven.f
fsrc=find.f getfkt.f abinit.f masses.f sparse_ma28.f naray.f weakread.f build.f locate.f partfun.f bilenshort.f weakintp.f rate.f rhside.f jacob.f dtnuc.f nse.f leqs.f rated.f solven.f

programname=sph
else
src=$(prog).c
fsrc=$(prog).f
programname=$(prog)
endif

treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

$(objdir)/grav$(objsuf) : grav.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c grav.c
	-@mv grav$(objsuf) $(objdir)

$(objdir)/sph$(objsuf) : sph.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c sph.c
	-@mv sph$(objsuf) $(objdir)

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/find$(objsuf) : find.f
	$(FC) $(FFLAGS) $(LDFLAGS) -c find.f
	-@mv find$(objsuf) $(objdir)

$(objdir)/abinit$(objsuf) : abinit.f
	$(FC) $(FFLAGS) $(LDFLAGS) -c abinit.f
	-@mv abinit$(objsuf) $(objdir)

$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/tree.h
$(objdir)/cofm$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/cofm$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics.h ndim.h
$(objdir)/cofm$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/fastflpt.h $(treedir)/include/protos.h
$(objdir)/grav$(objsuf): physics.h $(treedir)/include/tree.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/grav$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h ndim.h
$(objdir)/grav$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/vop.h
$(objdir)/grav$(objsuf): $(treedir)/include/tensop.h $(treedir)/include/fastflpt.h $(treedir)/include/Msgs.h
$(objdir)/mac$(objsuf): physics.h $(treedir)/include/tree.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h $(treedir)/include/Assert.h
$(objdir)/mac$(objsuf): $(treedir)/include/pqsort.h ndim.h $(treedir)/include/physics_generic.h
$(objdir)/mac$(objsuf): $(treedir)/include/vop.h $(treedir)/include/fastflpt.h $(treedir)/include/Msgs.h
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/fastflpt.h $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/main$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/SDF.h
$(objdir)/main$(objsuf): $(treedir)/include/protos.h $(treedir)/include/macr.h $(treedir)/include/malloc.h
$(objdir)/main$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/SDFwrite.h $(treedir)/include/SDFread.h
$(objdir)/main$(objsuf): $(treedir)/include/timers.h physics.h $(treedir)/include/tree.h 
$(objdir)/main$(objsuf): $(treedir)/include/key.h
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/stk.h $(treedir)/include/chn.h $(treedir)/include/pqsort.h ndim.h
$(objdir)/main$(objsuf): $(treedir)/include/physics_generic.h physics_sph.h $(treedir)/include/vop.h
$(objdir)/main$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/singlio.h $(treedir)/include/mpmy.h
$(objdir)/main$(objsuf): $(treedir)/include/mpmy_io.h $(treedir)/include/mpmy_abnormal.h
$(objdir)/main$(objsuf): $(treedir)/include/gc.h $(treedir)/include/files.h $(treedir)/include/getparam.h
$(objdir)/main$(objsuf): $(treedir)/include/verify.h $(treedir)/include/randoms.h $(treedir)/include/ring.h
$(objdir)/main$(objsuf): $(treedir)/include/decomp.h $(treedir)/include/image.h $(treedir)/include/memfile.h
$(objdir)/main$(objsuf): integrate.h
$(objdir)/main$(objsuf): $(treedir)/sph+nln/nrutil.h $(treedir)/sph+nln/cool.h
$(objdir)/physics$(objsuf): physics.h $(treedir)/include/tree.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/physics$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h ndim.h
$(objdir)/physics$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/physics_generic.c
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/protos.h $(treedir)/include/mpmy.h $(treedir)/include/vop.h
$(objdir)/physics$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/verify.h $(treedir)/include/files.h
$(objdir)/physics$(objsuf): $(treedir)/include/gc.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): physics.h $(treedir)/include/tree.h
$(objdir)/print$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/print$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/print$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h ndim.h
$(objdir)/print$(objsuf): $(treedir)/include/physics_generic.h $(treedir)/include/protos.h
$(objdir)/print$(objsuf): $(treedir)/include/vop.h
$(objdir)/integrate$(objsuf): ndim.h integrate.h $(treedir)/include/vop.h
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf): $(treedir)/include/tree.h
$(objdir)/sphcofm$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf):
$(objdir)/sphcofm$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/sphcofm$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/sphcofm$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_sph.h
$(objdir)/sphcofm$(objsuf): $(treedir)/include/vop.h $(treedir)/include/Msgs.h $(treedir)/include/fastflpt.h
$(objdir)/sphcofm$(objsuf): $(treedir)/include/protos.h
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf): physics_sph.h $(treedir)/include/tree.h
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/sph$(objsuf):
$(objdir)/sph$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h $(treedir)/include/Assert.h
$(objdir)/sph$(objsuf): $(treedir)/include/pqsort.h $(treedir)/include/vop.h $(treedir)/include/fastflpt.h
$(objdir)/sph$(objsuf): physics.h ndim.h $(treedir)/include/physics_generic.h 
$(objdir)/sph$(objsuf): $(treedir)/sph+nln/nrutil.h $(treedir)/sph+nln/cool.h
$(objdir)/sph$(objsuf): $(treedir)/sph+nln/units.h
$(objdir)/sphinit$(objsuf):
$(objdir)/sphinit$(objsuf):
$(objdir)/sphinit$(objsuf):
$(objdir)/sphinit$(objsuf): $(treedir)/include/randoms.h
$(objdir)/sphinit$(objsuf): physics.h $(treedir)/include/tree.h 
$(objdir)/sphinit$(objsuf):
$(objdir)/sphinit$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/sphinit$(objsuf):
$(objdir)/sphinit$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/sphinit$(objsuf):
$(objdir)/sphinit$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/sphinit$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h ndim.h
$(objdir)/sphinit$(objsuf): $(treedir)/include/physics_generic.h physics_sph.h $(treedir)/include/vop.h
$(objdir)/sphinit$(objsuf): $(treedir)/include/singlio.h $(treedir)/include/fastflpt.h
$(objdir)/sphinit$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/gc.h $(treedir)/include/Msgs.h
$(objdir)/sphinit$(objsuf): $(treedir)/include/SDF.h $(treedir)/include/SDFread.h SDFreadf.h
$(objdir)/sphinit$(objsuf): $(treedir)/sph+nln/nrutil.h $(treedir)/sph+nln/cool.h
$(objdir)/physics_sph$(objsuf): physics_sph.h $(treedir)/include/tree.h
$(objdir)/physics_sph$(objsuf):
$(objdir)/physics_sph$(objsuf):
$(objdir)/physics_sph$(objsuf):
$(objdir)/physics_sph$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/physics_sph$(objsuf):
$(objdir)/physics_sph$(objsuf):
$(objdir)/physics_sph$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/physics_sph$(objsuf):
$(objdir)/physics_sph$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/physics_sph$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h
$(objdir)/physics_sph$(objsuf): $(treedir)/include/verify.h $(treedir)/include/mpmy.h $(treedir)/include/gc.h
$(objdir)/physics_sph$(objsuf): $(treedir)/include/vop.h
$(objdir)/physics_sph$(objsuf): $(treedir)/sph+nln/units.h
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf): physics_sph.h
$(objdir)/sphprint$(objsuf): $(treedir)/include/tree.h
$(objdir)/sphprint$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/sphprint$(objsuf):
$(objdir)/sphprint$(objsuf): $(treedir)/include/stk.h $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/sphprint$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/sphprint$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h
$(objdir)/sphprint$(objsuf): $(treedir)/include/protos.h $(treedir)/include/vop.h
$(objdir)/sphplus$(objsuf): $(treedir)/include/stk.h
$(objdir)/sphplus$(objsuf):
$(objdir)/sphplus$(objsuf):
$(objdir)/sphplus$(objsuf):
$(objdir)/sphplus$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/error.h
$(objdir)/sphplus$(objsuf):
$(objdir)/sphplus$(objsuf): $(treedir)/include/gccextensions.h physics.h $(treedir)/include/tree.h
$(objdir)/sphplus$(objsuf): $(treedir)/include/timers.h $(treedir)/include/key.h
$(objdir)/sphplus$(objsuf):
$(objdir)/sphplus$(objsuf): $(treedir)/include/chn.h $(treedir)/include/Assert.h $(treedir)/include/pqsort.h
$(objdir)/sphplus$(objsuf): ndim.h $(treedir)/include/physics_generic.h physics_sph.h
$(objdir)/sphplus$(objsuf): $(treedir)/include/vop.h
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf):
$(objdir)/SDFreadf$(objsuf): $(treedir)/include/mpmy.h
$(objdir)/SDFreadf$(objsuf): $(treedir)/include/timers.h $(treedir)/include/SDF.h $(treedir)/include/Assert.h
$(objdir)/SDFreadf$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/SDFreadf$(objsuf): $(treedir)/include/bigmalloc.h $(treedir)/include/Msgs.h $(treedir)/include/verify.h
$(objdir)/SDFreadf$(objsuf): SDFreadf.h $(treedir)/include/gc.h $(treedir)/include/singlio.h
$(objdir)/nrutils$(objsuf):
$(objdir)/coolandburn$(objsuf): $(treedir)/include/error.h $(treedir)/include/singlio.h $(treedir)/include/fastflpt.h $(treedir)/sph+nln/nrutil.h $(treedir)/sph+nln/cool.h
$(objdir)/params$(objsuf):
