##### Application-specific stuff goes here

src=	swampiIH.c
programname=swampi
treedir=..

appexcludes:=-name data
ifeq ($(ARCH), sun5)
LOADLIBES:=-lsocket -lnsl
endif

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# DO NOT DELETE

