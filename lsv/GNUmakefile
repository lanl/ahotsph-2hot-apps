ifeq (X$(TREEHOME), X)
TREEHOME=$(HOME)/tree16.oz
endif
##### Application-specific stuff goes here

src=	lsv.c
programname=lsv
treedir=..
treedir_sed=\.\.

appexcludes:=-name data
ifeq ($(ARCH), sun5)
LOADLIBES:=-lsocket -lnsl
endif

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

