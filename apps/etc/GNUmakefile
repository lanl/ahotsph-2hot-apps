ifeq (X$(TREEHOME), X)
TREEHOME=../../
endif
##### Application-specific stuff goes here
treedir=$(TREEHOME)

ifeq ($(prog),link)
src=link.c schecter.c
else
ifeq ($(prog),pview)
src=pview.c
EXTRACFLAGS=-I$(HOME)/lib
LOADLIBES=$(HOME)/lib/libvogl.a -lX11 -lm
LDFLAGS=-static
endif
ifeq ($(prog),rotate_movie)
src=rotate_movie.c decomp.c
else
ifeq ($(prog),rho_image)
src=rho_image.c rimage.c
else
src=$(prog).c
endif
endif
endif

programname=$(prog)

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.
