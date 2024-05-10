TREEHOME=..
treedir_sed=\.\.
##### Application-specific stuff goes here

ifdef PROFILE
EXTRACFLAGS = -I$(PREFIX)/include -DGPERF
PRELIBS = -lrt -Wl,-rpath -Wl,$(PREFIX)/lib -L$(PREFIX)/lib -lprofiler
endif

ifeq ($(ARCH),$(filter $(ARCH), cudaxk7 cudaxk7-g cudai7 cudaxk7-4.6 cudaxk7-4.8))
cusrc = cuda.cu pgrav.cu trace.cu
PRELIBS +=-L/opt/cudatoolkit-5.5/lib64 -lcudart
endif

programname=nlna

src = cofm.c grav.c mac.c main.c physics.c print.c output.c integrate.c ewald_le2.c do_grav_sse4.c do_grav_avx8.c pgrav_vec.c ewald.c grav_n2.c version.c

treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

$(objdir)/pgrav_vec$(objsuf) : pgrav_vec.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c pgrav_vec.c
	mv pgrav_vec$(objsuf) $(objdir)

$(objdir)/do_grav_sse4$(objsuf) : do_grav_sse4.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c do_grav_sse4.c
	mv do_grav_sse4$(objsuf) $(objdir)

$(objdir)/do_grav_avx8$(objsuf) : do_grav_avx8.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c do_grav_avx8.c
	mv do_grav_avx8$(objsuf) $(objdir)

$(objdir)/ewald_le2$(objsuf) : ewald_le2.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c ewald_le2.c
	mv ewald_le2$(objsuf) $(objdir)

.PHONY: version.proto
version.proto:
	@echo char Version\[\] = \"`git describe --tags --dirty`\"\; > version.proto
	@echo char Arch\[\] = \"`printenv ARCH`\"\; >> version.proto
	@echo char Compiled_date\[\] = __DATE__\; >> version.proto
	@echo char Compiled_time\[\] = __TIME__\; >> version.proto
	@echo char Compiler\[\] = __VERSION__\; >> version.proto

version.c: version.proto
	@cmp -s $< $@ || cp -p $< $@

# DO NOT DELETE THIS LINE -- make depend depends on it.

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
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/key.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/stk.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/cofm$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h order.h physics.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): vec.h segment.h
$(objdir)/cofm$(objsuf): physics_generic.h $(treedir)/include/vop.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/fastflpt.h $(treedir)/include/protos.h
$(objdir)/grav$(objsuf): physics.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/key.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/stk.h $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/grav$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/grav$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/grav$(objsuf): physics_generic.h $(treedir)/include/vop.h
$(objdir)/grav$(objsuf): $(treedir)/include/tensop.h $(treedir)/include/fastflpt.h $(treedir)/include/Msgs.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): order.h physics.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/tree.h
$(objdir)/mac$(objsuf): $(treedir)/include/timers.h
$(objdir)/mac$(objsuf): $(treedir)/include/key.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/stk.h
$(objdir)/mac$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/mac$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h $(treedir)/include/Assert.h
$(objdir)/mac$(objsuf): $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/mac$(objsuf): physics_generic.h $(treedir)/include/vop.h $(treedir)/include/fastflpt.h
$(objdir)/mac$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/mpmy.h $(treedir)/include/ewald_le.h
$(objdir)/mac$(objsuf): gravcuda.h
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
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/fastflpt.h
$(objdir)/main$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/error.h
$(objdir)/main$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/SDF.h
$(objdir)/main$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/SDFwrite.h $(treedir)/include/SDFread.h
$(objdir)/main$(objsuf): $(treedir)/include/timers.h
$(objdir)/main$(objsuf): order.h physics.h $(treedir)/include/tree.h
$(objdir)/main$(objsuf): $(treedir)/include/key.h
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/stk.h
$(objdir)/main$(objsuf): $(treedir)/include/chn.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/main$(objsuf): physics_generic.h $(treedir)/include/vop.h
$(objdir)/main$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/singlio.h $(treedir)/include/mpmy.h
$(objdir)/main$(objsuf): $(treedir)/include/mpmy_io.h $(treedir)/include/mpmy_abnormal.h
$(objdir)/main$(objsuf): $(treedir)/include/gc.h $(treedir)/include/files.h $(treedir)/include/getparam.h
$(objdir)/main$(objsuf): $(treedir)/include/verify.h $(treedir)/include/randoms.h $(treedir)/include/decomp.h
$(objdir)/main$(objsuf): $(treedir)/include/image.h $(treedir)/include/memfile.h $(treedir)/include/ewald_le.h
$(objdir)/main$(objsuf): $(treedir)/include/cosmo.h integrate.h output.h $(treedir)/include/cpu.h
$(objdir)/main$(objsuf): gravcuda.h version.h
$(objdir)/physics$(objsuf): physics.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/key.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/stk.h $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/physics$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/physics$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/physics$(objsuf): physics_generic.h physics_generic.c
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
$(objdir)/print$(objsuf): physics.h
$(objdir)/print$(objsuf): $(treedir)/include/tree.h
$(objdir)/print$(objsuf): $(treedir)/include/timers.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/key.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/stk.h
$(objdir)/print$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/print$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/print$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/print$(objsuf): physics_generic.h $(treedir)/include/protos.h
$(objdir)/print$(objsuf): $(treedir)/include/vop.h
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf): $(treedir)/include/SDF.h
$(objdir)/output$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/SDFwrite.h physics.h
$(objdir)/output$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/output$(objsuf): $(treedir)/include/key.h
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf): $(treedir)/include/stk.h $(treedir)/include/chn.h $(treedir)/include/Assert.h
$(objdir)/output$(objsuf): $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/output$(objsuf): physics_generic.h $(treedir)/include/vop.h $(treedir)/include/Msgs.h
$(objdir)/output$(objsuf): $(treedir)/include/singlio.h $(treedir)/include/mpmy.h $(treedir)/include/mpmy_io.h
$(objdir)/output$(objsuf): $(treedir)/include/cosmo.h integrate.h output.h version.h
$(objdir)/output$(objsuf): $(treedir)/include/Random123/threefry.h
$(objdir)/output$(objsuf): $(treedir)/include/Random123/features/compilerfeatures.h
$(objdir)/output$(objsuf): $(treedir)/include/Random123/features/gccfeatures.h
$(objdir)/output$(objsuf): $(treedir)/include/Random123/array.h
$(objdir)/output$(objsuf): $(treedir)/include/Random123/features/sse.h
$(objdir)/output$(objsuf): $(treedir)/include/Random123/u01fixedpt.h
$(objdir)/integrate$(objsuf): physics.h
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf): $(treedir)/include/tree.h
$(objdir)/integrate$(objsuf): $(treedir)/include/timers.h
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf): $(treedir)/include/key.h
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf): $(treedir)/include/stk.h
$(objdir)/integrate$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/integrate$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/integrate$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/integrate$(objsuf): physics_generic.h
$(objdir)/integrate$(objsuf): $(treedir)/include/vop.h $(treedir)/include/cosmo.h integrate.h
$(objdir)/integrate$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf): physics.h $(treedir)/include/tree.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/timers.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/key.h
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf):
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/stk.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/ewald_le2$(objsuf): physics_generic.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/vop.h moments.h $(treedir)/include/mpmy.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/do_grav_sse4$(objsuf): order.h vec.h
$(objdir)/pgrav_vec$(objsuf): vec.h
$(objdir)/ewald$(objsuf):
$(objdir)/ewald$(objsuf):
$(objdir)/ewald$(objsuf):
$(objdir)/ewald$(objsuf):
$(objdir)/ewald$(objsuf): $(treedir)/include/vop.h
$(objdir)/ewald$(objsuf): $(treedir)/include/error.h $(treedir)/include/gccextensions.h
$(objdir)/grav_n2$(objsuf): physics.h
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf): $(treedir)/include/key.h
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf): $(treedir)/include/stk.h $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h vec.h segment.h
$(objdir)/grav_n2$(objsuf): physics_generic.h $(treedir)/include/vop.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/randoms.h $(treedir)/include/ring.h $(treedir)/include/mpmy.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/Msgs.h $(treedir)/include/singlio.h
