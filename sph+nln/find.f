c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine find(name,id,inam,blank)
      implicit none

      integer*4 niso, n, id
      parameter(niso = 5410)
      character*5 name
      character*5 inam(niso),blank
c-----------------------------------------------------

      if( name .ne. blank )then
        do n = 1, niso
          if( name .eq. inam(n) )then
             id = n
             goto 10
           endif
         enddo
       else
          id = 0
       endif

10     continue

       return
       end




