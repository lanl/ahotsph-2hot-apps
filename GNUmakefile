xctarname:=tree19

# Make.$(ARCH) sets many of the variables that are then used in
# Make.generic.  Leaving it out can cause problems, for example,
# Make.$(ARCH) wants to redefine 'TAR'.  On the other hand, leaving
# it in can cause problems (apparently) with older versions of GNUmake
# because, e.g., assignements like LOADLIBES:=-foo $(LOADLIBES) get
# 'executed' multiple times.  Once in this Makefile, and once in the
# daughter makefiles.

###include Make-common/Make.$(ARCH)

# we should also go into SDFcvt, SDF2fld, commtst, lsv, lsvtst, (anything
# else?) and build them.  Possibly under a different target?

all: All

# Make.generic has targets for clean, all, etc., so we need to
# spell them a little differently in this file...
include 2hot/Make-common/Make.generic

subdirs:= 2hot/libsw 2hot/libSDF 2hot/libtree 2hot/libmpmy 2hot/nln sph+nln snsph wvt
#nln # currently has compile errors in sse stuff

All:
	for dir in $(subdirs); do (cd $$dir; $(MAKE) all); done

Depends :
	for dir in $(subdirs); do (cd $$dir; $(MAKE) depends); done

Clean : 
	$(RM) $(objdir)/*
	for dir in $(subdirs); do (cd $$dir; $(MAKE) clean); done

Clean-all:
	$(RM) $(objdir)/*
	for dir in $(subdirs); do (cd $$dir; $(MAKE) clean-all); done

.PHONY: version.proto
version.proto:
	@echo char Version\[\] = \"`git describe --tags --dirty`\"\; > version.proto
	@echo char Arch\[\] = \"`printenv ARCH`\"\; >> version.proto
	@echo char Compiled_date\[\] = __DATE__\; >> version.proto
	@echo char Compiled_time\[\] = __TIME__\; >> version.proto
	@echo char Compiler\[\] = __VERSION__\; >> version.proto

version.c: version.proto
	@cmp -s $< $@ || cp -p $< $@

