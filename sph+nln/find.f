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




