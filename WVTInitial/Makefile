TREEHOME=..

ifndef prog
src=grav.c interp_cylgrid.c mac.c physics.c physics_sph.c sphcofm.c	\
	sphinit.c sph.c sphprint.c sphplus.c ghosts.c com.c ranlib.c	\
	wvt.c polint.c initial.c

programname=interp_cylgrid
else
src=$(prog).c
programname=$(prog)
endif

treedir=$(TREEHOME)

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic
