TREEHOME=..
treedir_sed=\.\.
# Application-code might want to use  something like instead:
#ifeq (X$(TREEHOME), X)
#TREEHOME=$(HOME)/nbody/tree16
#endif
##### Application-specific stuff goes here

src=	cofm_n.c	grav_nv.c	main_x.c	physics_n.c \
	print.c 	grav_ring.c	# decomp.c

programname=relerr_x
treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE THIS LINE -- make depend depends on it.

