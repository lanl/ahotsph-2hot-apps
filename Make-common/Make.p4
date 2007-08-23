defaultCC:=gcc

CC_SPECIFIC:=-g -Wall
ARCH_SPECIFIC:=-D__iX86__=686 -ffloat-store -ffast-math -malign-double -march=pentium4
OPTIMIZE=-O2
override AGGRESSIVE_OPT=-O3
LEX:=flex
YACC:=bison -y

include $(treedir)/Make-common/Make.default

swsrc:=lsv.c swampi.c

ifeq ($(PAROS),mpi)
LOADLIBES:=-L/ss1-raid/local/mpich-1.2.4-ssh/lib -lmpich -lpmpich
PAROSCFLAGS:=-I/ss1-raid/local/mpich-1.2.4-ssh/include
endif

ifeq ($(PAROS),lam)
LOADLIBES:=-L/usr/local/lam/lib -lmpi -ltstdio -ltrillium -largs -lt
PAROSCFLAGS:=-I/usr/local/lam/h
endif

ifeq ($(PAROS),pvm)
LOADLIBES:=-L/usr/local/pvm3/lib/LINUX -lpvm3
PAROSCFLAGS:=-I/usr/local/pvm3/include
endif
