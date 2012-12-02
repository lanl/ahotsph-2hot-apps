TREEHOME=..
treedir_sed=\.\.
##### Application-specific stuff goes here

ARCH=amd6100

EXTRACFLAGS=-I/usr/projects/cosmo/amd6100/include
PRELIBS=-lrt -Wl,-rpath -Wl,/usr/projects/cosmo/amd6100/lib -L/usr/projects/cosmo/amd6100/lib -lprofiler
programname=../devel/nln

src = cofm.c grav.c mac.c main.c physics.c print.c output.c integrate.c ewald_le2.c do_grav_sse4.c pgrav_sse4.c ewald.c grav_n2.c version.c
treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

$(objdir)/grav$(objsuf) : grav.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c grav.c
	mv grav$(objsuf) $(objdir)

$(objdir)/do_grav_sse4$(objsuf) : do_grav_sse4.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c do_grav_sse4.c
	mv do_grav_sse4$(objsuf) $(objdir)

$(objdir)/ewald_le2$(objsuf) : ewald_le2.c
	$(defaultCC) $(CFLAGS) $(AGGRESSIVE_OPT) -c ewald_le2.c
	mv ewald_le2$(objsuf) $(objdir)

.PHONY: version.proto
version.proto:
	@echo char Version\[\] = \"`git describe --dirty`\"\; > version.proto
	@echo char Compiled_date\[\] = __DATE__\; >> version.proto
	@echo char Compiled_time\[\] = __TIME__\; >> version.proto

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
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/tree.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/timers.h
$(objdir)/cofm$(objsuf): $(treedir)/include/key.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): $(treedir)/include/stk.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/cofm$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/cofm$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h order.h physics.h
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf):
$(objdir)/cofm$(objsuf): physics_generic.h
$(objdir)/cofm$(objsuf): $(treedir)/include/vop.h $(treedir)/include/Msgs.h $(treedir)/include/fastflpt.h
$(objdir)/cofm$(objsuf): $(treedir)/include/protos.h
$(objdir)/grav$(objsuf): physics.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/tree.h
$(objdir)/grav$(objsuf): $(treedir)/include/timers.h
$(objdir)/grav$(objsuf): $(treedir)/include/key.h
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf):
$(objdir)/grav$(objsuf): $(treedir)/include/stk.h
$(objdir)/grav$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/grav$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/grav$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_generic.h
$(objdir)/grav$(objsuf): $(treedir)/include/vop.h $(treedir)/include/tensop.h $(treedir)/include/fastflpt.h
$(objdir)/grav$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): order.h physics.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/tree.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/timers.h
$(objdir)/mac$(objsuf): $(treedir)/include/key.h
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf):
$(objdir)/mac$(objsuf): $(treedir)/include/stk.h
$(objdir)/mac$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/mac$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h $(treedir)/include/Assert.h
$(objdir)/mac$(objsuf): $(treedir)/include/pqsort.h physics_generic.h $(treedir)/include/vop.h
$(objdir)/mac$(objsuf): $(treedir)/include/fastflpt.h $(treedir)/include/Msgs.h $(treedir)/include/mpmy.h
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
$(objdir)/main$(objsuf): $(treedir)/include/protos.h $(treedir)/include/macr.h $(treedir)/include/malloc.h
$(objdir)/main$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/SDFwrite.h $(treedir)/include/SDFread.h
$(objdir)/main$(objsuf): $(treedir)/include/timers.h
$(objdir)/main$(objsuf): order.h physics.h $(treedir)/include/tree.h
$(objdir)/main$(objsuf): $(treedir)/include/key.h
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf):
$(objdir)/main$(objsuf): $(treedir)/include/stk.h $(treedir)/include/chn.h $(treedir)/include/pqsort.h
$(objdir)/main$(objsuf): physics_generic.h $(treedir)/include/vop.h $(treedir)/include/Msgs.h
$(objdir)/main$(objsuf): $(treedir)/include/singlio.h $(treedir)/include/mpmy.h $(treedir)/include/mpmy_io.h
$(objdir)/main$(objsuf): $(treedir)/include/mpmy_abnormal.h $(treedir)/include/gc.h $(treedir)/include/files.h
$(objdir)/main$(objsuf): $(treedir)/include/getparam.h $(treedir)/include/verify.h $(treedir)/include/randoms.h
$(objdir)/main$(objsuf): $(treedir)/include/decomp.h $(treedir)/include/image.h $(treedir)/include/memfile.h
$(objdir)/main$(objsuf): $(treedir)/include/ewald_le.h $(treedir)/include/cosmo.h integrate.h output.h
$(objdir)/physics$(objsuf): physics.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/tree.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/timers.h
$(objdir)/physics$(objsuf): $(treedir)/include/key.h
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/stk.h
$(objdir)/physics$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/physics$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/physics$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_generic.h
$(objdir)/physics$(objsuf): physics_generic.c
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf):
$(objdir)/physics$(objsuf): $(treedir)/include/protos.h
$(objdir)/physics$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/vop.h $(treedir)/include/Msgs.h
$(objdir)/physics$(objsuf): $(treedir)/include/verify.h $(treedir)/include/files.h $(treedir)/include/gc.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): physics.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/tree.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/timers.h
$(objdir)/print$(objsuf): $(treedir)/include/key.h
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf):
$(objdir)/print$(objsuf): $(treedir)/include/stk.h
$(objdir)/print$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/print$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/print$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_generic.h
$(objdir)/print$(objsuf): $(treedir)/include/protos.h $(treedir)/include/vop.h
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
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf): $(treedir)/include/SDF.h
$(objdir)/output$(objsuf): $(treedir)/include/malloc.h $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/output$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/SDFwrite.h physics.h
$(objdir)/output$(objsuf): $(treedir)/include/tree.h $(treedir)/include/timers.h
$(objdir)/output$(objsuf): $(treedir)/include/key.h
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf):
$(objdir)/output$(objsuf): $(treedir)/include/stk.h
$(objdir)/output$(objsuf): $(treedir)/include/chn.h $(treedir)/include/Assert.h $(treedir)/include/pqsort.h
$(objdir)/output$(objsuf): physics_generic.h $(treedir)/include/vop.h $(treedir)/include/Msgs.h
$(objdir)/output$(objsuf): $(treedir)/include/singlio.h $(treedir)/include/mpmy.h $(treedir)/include/mpmy_io.h
$(objdir)/output$(objsuf): $(treedir)/include/cosmo.h integrate.h output.h
$(objdir)/integrate$(objsuf): physics.h
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf): $(treedir)/include/tree.h
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf): $(treedir)/include/timers.h
$(objdir)/integrate$(objsuf): $(treedir)/include/key.h
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf):
$(objdir)/integrate$(objsuf): $(treedir)/include/stk.h
$(objdir)/integrate$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/integrate$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/integrate$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_generic.h
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
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/stk.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_generic.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/vop.h moments.h $(treedir)/include/mpmy.h
$(objdir)/ewald_le2$(objsuf): $(treedir)/include/Msgs.h
$(objdir)/ewald$(objsuf):
$(objdir)/ewald$(objsuf):
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
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf): $(treedir)/include/tree.h
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf): $(treedir)/include/timers.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/key.h
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf):
$(objdir)/grav_n2$(objsuf): $(treedir)/include/stk.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/Malloc.h $(treedir)/include/error.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/gccextensions.h $(treedir)/include/chn.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/Assert.h $(treedir)/include/pqsort.h physics_generic.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/vop.h $(treedir)/include/randoms.h $(treedir)/include/ring.h
$(objdir)/grav_n2$(objsuf): $(treedir)/include/mpmy.h $(treedir)/include/Msgs.h $(treedir)/include/singlio.h
