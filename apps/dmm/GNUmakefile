ifeq (X$(TREEHOME), X)
TREEHOME=../../
endif
##### Application-specific stuff goes here
treedir=$(TREEHOME)

ifeq ($(prog),sdmm)
src=sdmm.c
programname=sdmm
else
src=dmm.c
programname=dmm
endif

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.
