ifeq (X$(TREEHOME), X)
TREEHOME=../../
endif
##### Application-specific stuff goes here
treedir=$(TREEHOME)

ifeq ($(prog),mrpick)
src=mrpick.c
programname=mrpick
else
src=rpick.c
programname=rpick
endif

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.
