c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


      subroutine partfun(t9)

      implicit none

      include 'dimen'

c     Sets up the partition functions for the target and product
c     nuclei for each reaction by using Friedel's NETWINV file.  The
c     partition functions are interpolated to the temperature T9
c     using the values for partition function at 24 temp specified in
c     NETWINV.

      include 'comcsolve'
      include 'cbug'

      Real*8 mxcess(nnuc), p(nnuc,24), t(24), w(nnuc), anuc(nnuc), t9 
      real*8 p1, p2, p3

      integer*4 pfread
      integer*4  nt(24)
      integer*4 z(nnuc), n(nnuc), mm, nn, i, j

      character*5 nuc
      character*5 cdummy
      character*72 line

      data pfread/1/
      save
c------------------------------------------------------------------
      if( pfread .eq. 1 )then
c..first call to pfread
         open(12, file='netwinv', status='old')
c     Read in the temperature array on the second line of NETWINV:
c     These are the T's at which Friedel's partition functions
c     are defined.
         read (12,'(a72)')line
c..   deconstruct them into a useable form
         read (12,'(24i3)') (nt(i), i=1, 24)
         do i = 1, 23
            t(i) = dble( nt(i) )/100.0d0
         enddo
         t(24) = dble( nt(24) )/10.0d0

c     Read in the info from NETWINV
c     Skip the nuclear symbols (look at NETWINV)
         
         do i=1,itot
            read (12,'(a5)')cdummy
         enddo
         do i = 1, nnuc
            if( pfread .eq. 1 )then
               read (12,'(a5,4x,f9.3,i3,1x,i3,2x,f5.1,2x,f7.3)') 
     1              nuc, anuc(i), z(i), n(i), w(i), mxcess(i)
               read (12,'(8f9.2)') (p(i,j), j = 1, 24)
            endif
         enddo

         close(12)

c..reset to avoid reading partition function data again
         pfread = 0

      endif

c     Identify indices m and n such that T9 lies between T(m)
c     and T(n).  This information is needed for interpolation later

      mm = 1
      nn = mm
      do i = 2, 24
         if (t9 .ge. t(i-1) .and. t9 .le. t(i) ) then
            mm = i-1
            nn = i
         else
            mm = mm
            nn = nn
         endif
      enddo
      if(t9 .gt. t(24))then
         mm = 24
         nn = mm
      endif


c     Interpolate to obtain the value of partition function 
c     at the given temperature T9 
      do i = 1, nnuc
         if ( p(i,mm) .eq. p(i,nn) .or. t9 .eq. t(mm) ) then
            pf(i) = p(i,mm)
         else
            pf(i) = (p(i,mm)-p(i,nn)) * (t9-t(mm)) / (t(mm)-t(nn))
            pf(i) = pf(i)+p(i,mm)
         endif
      enddo
      inuc = i-1


c     Rearrange with n, p, and alpha at the end
      p1 = pf(1)
      p2 = pf(2)
      p3 = pf(3)
      
      do i = 4, nnuc
         pf(i-3) = pf(i)
      enddo
      
      pf(nnuc-2) = p1
      pf(nnuc-1) = p2
      pf(nnuc)   = p3

      return
      end
