#ARCH = p4icc-coyote
#PAROS = mpi
#ARCH=p4icc
#PAROS=swampi
treedir:=../2hot

ifndef prog
src= grav.c interp_cylgrid.c mac.c physics.c physics_sph.c \
	sphcofm.c sphinit.c sph.c sphprint.c sphplus.c ghosts.c com.c \
	ranlib.c wvt.c polint.c initial.c

programname=interp_cylgrid
else
src=$(prog).c
fsrc=$(prog).f
programname=$(prog)
endif

appexcludes:=-name data

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

$(objdir)/grav$(objsuf) : grav.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c grav.c
	-@mv grav$(objsuf) $(objdir)

$(objdir)/sph$(objsuf) : sph.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c sph.c
	-@mv sph$(objsuf) $(objdir)

# DO NOT DELETE THIS LINE -- make depend depends on it.
#
#.PHONY: all clean
#
#all: ${PROGS}
#
#clean: 
#	-${RM} *.o ${PROGS}
#
#
#interp: grav.o interp.o mac.o physics.o physics_sph.o sphcofm.o sphinit.o \
#	sph.o sphprint.o sphplus.o ghosts.o com.o ranlib.o wvt.o polint.o \
#	initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_seq.o -o $@ \
#	$(LDLIBS)
#
#interp.swampi: grav.o interp.o mac.o physics.o physics_sph.o sphcofm.o \
#	sphinit.o sph.o sphprint.o sphplus.o ghosts.o com.o ranlib.o wvt.o \
#	polint.o initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_swampi.o -o $@ \
#	$(LDLIBS)
#
#interp.mpi: grav.o interp.o mac.o physics.o physics_sph.o sphcofm.o sphinit.o \
#	sph.o sphprint.o sphplus.o ghosts.o com.o ranlib.o wvt.o polint.o \
#	initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_mpi.o -o $@ \
#	$(LDLIBS)
#
#gridtosdf: grav.o gridtosdf.o mac.o physics.o physics_sph.o sphcofm.o \
#	sphinit.o sph.o sphprint.o sphplus.o ghosts.o com.o ranlib.o wvt.o \
#	polint.o initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_seq.o -o $@ \
#	$(LDLIBS)
#
#wvt_regrid: grav.o wvt_regrid.o mac.o physics.o physics_sph.o sphcofm.o \
#	sphinit.o sph.o sphprint.o sphplus.o ghosts.o com.o ranlib.o wvt.o \
#	polint.o initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_seq.o -o $@ \
#	$(LDLIBS)
#
#interp_cylgrid: grav.o interp_cylgrid.o mac.o physics.o physics_sph.o \
#	sphcofm.o sphinit.o sph.o sphprint.o sphplus.o ghosts.o com.o \
#	ranlib.o wvt.o polint.o initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_seq.o -o $@ \
#	$(LDLIBS)
#
#interp_cylgrid.mpi: grav.o interp_cylgrid.o mac.o physics.o physics_sph.o \
#	sphcofm.o sphinit.o sph.o sphprint.o sphplus.o ghosts.o com.o \
#	ranlib.o wvt.o polint.o initial.o
#	$(CC) $(LDFLAGS) $(PAROSCFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_mpi.o -o $@ \
#	$(LDLIBS) ${LOADLIBES} ${PAROSCFLAGS}
#
#interp_cylgrid.swampi: grav.o interp_cylgrid.o mac.o physics.o physics_sph.o \
#	sphcofm.o sphinit.o sph.o sphprint.o sphplus.o ghosts.o com.o \
#	ranlib.o wvt.o polint.o initial.o
#	$(CC) $(LDFLAGS) $^ ${TREEHOME}/Objfiles/${ARCH}/mpmy_swampi.o -o $@ \
#	$(LDLIBS)
#
