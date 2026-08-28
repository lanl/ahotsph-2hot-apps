# Astrophysical Hashed Oct Tree-based Smoothed Particle Hydrodynamics code submodule: '2hot-apps'

## Description

This is the original SPH + tree code that I received many moons ago in grad school. It is based on the old version of the hashed oct tree code. It will eventually be replaced by what's in the '2hot-apps' repo (once that's working).

# Repository Contents

- nln  
  The O(n*log(n)) n-body code written by Mike Warren, adapted to work with 2HOT.

- sph+nln  
  An SPH implementation on top of 2HOT. This started out as being implemented on top of HOT, as is currently being ported to work with 2HOT.

- snsph  
  An SPH implementation on top of 2HOT geared towards supernova calculations. Less maintained than sph+nln. This started out as being implemented on top of HOT, as is currently being ported to work with 2HOT.

- wvt  
  A setup code for creating particle setups for running in sph+nln or snsph. This started out as being implemented on top of HOT, as is currently being ported to work with 2HOT.

# How to get

Either clone the repo directly:
```
git clone git@git@github.com:lanl/ahotsph-2hot-apps.git
```

or

clone the super-repo recursively:
```
git clone --recursive git@github.com:lanl/ahotsph.git

```

# How to build

## General

```
cd 2hot-apps
module load <supported compiler> <supported mpi>
<modify a make file in Make-Common/>
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

## Contributing
tbd

## License and copyright


BSD 3-Clause License

Copyright (c) 2026, Los Alamos National Laboratory. O5196.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
