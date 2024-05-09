# Contents

- 2hot
  The new version (2.0) of the Hashed Oct Tree.

- sph+nln
  An SPH implementation on top of 2HOT. This started out as being implemented on top of HOT, as is currently being ported to work with 2HOT.

- snsph
  An SPH implementation on top of 2HOT. This started out as being implemented on top of HOT, as is currently being ported to work with 2HOT.

- wvt
  A setup code for creating particle setups for running in sph+nln or snsph. This started out as being implemented on top of HOT, as is currently being ported to work with 2HOT.

# How to get

```
git clone --recursive git@gitlab.lanl.gov:SPH/2hot-apps
```

or

```
git clone git@gitlab.lanl.gov:SPH/2hot-apps
cd 2hot-apps
git submodule update --init --recursive
```

# How to build

## General

```
cd 2hot-apps
module load <supported compiler> <supported mpi>
make ARCH=<my-arch> PAROS=[mpi|seq] [CC=<c-compiler-cmd> FC=<fortran-compiler-cmd>]
```

## Rocinante

## Darwin - Haswell node

### Prerequisites
Build libunwind: 
- clone from the official libunwind repo or from the fork found in SPH/3rd-party (coming soon)
- load a GCC module
- follow the build instructions in the Readme, set the prefix to somewhere you can write to (e.g. $HOME, $HOME/my_installs, ...)
- note the path to the `lib` and `include` folders, let's call them LIBUNWIND_LIB AND LIBUNWIND_INCLUDE

```
module load gcc/12.2.0 openmpi/4.1.5-gcc_12.2.0
cd 2hot-apps
vim 2hot/Make.common/Make.x86_64
# edit so that LDFLAGS points to the lib folder in the path where you built libunwind: "LDFLAGS:=-L$(LIBUNWIND_LIB) ..."
# edit PAROSFLAGS so it point to the include folder from your libunwind install: "PAROSFLAGS:=-I$(LIBUNWIND_INCLUDE) ..."
make ARCH=x86_64 PAROS=mpi
```
