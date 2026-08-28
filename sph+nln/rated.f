c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine rated(t90,rho)

c..   fkt rate formula with T9, V derivatives
c..   deck 1, decay only

      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'

      real*8 t9, rho, t913, t953, t9ln, expon
      real*8 t92, t923, t943, factor, third
      real*8 ye, t90, t9lim

      parameter( third = 1.0d0/3.0d0, t9lim = 0.01d0 )

      integer*4 k, k1, k2, j
c----------------------------------------------------
c..   avoid bad fkt extrapolations to low temperature
      if( t90 .gt. t9lim )then
         t9 = t90
      else
         t9 = t9lim
      endif

c..   electron mole number Ye = Ne/rho*avagadro
      ye = 0.0d0
      do j = 1, nnuc
         ye = ye + y(j) * dble( nz(j) )
      enddo

      t913 = t9**third
      t953 = t9 * t913 *t913
      t9ln = dlog(t9)

      t92  = t9*t9
      t923 = t913*t913
      t943 = t913*t9

      k1 = k1deck(1)
      k2 = k2deck(2)
      do k = k1,k2
         expon = rcoef(1,k) + rcoef(2,k)/t9 + rcoef(3,k)/t913
     1        + rcoef(4,k)*t913 + rcoef(5,k)*t9
     2        + rcoef(6,k)*t953 + rcoef(7,k)*t9ln
         if( expon .gt. -300.0d0 )then
            sig(k) = exp(expon)
            factor = -rcoef(2,k)/t92 -third*rcoef(3,k)/t943
     1           +third*rcoef(4,k)/t923 + rcoef(5,k)
     2           +5.0d0*third*rcoef(6,k)*t923 + rcoef(7,k)/t9

            sigt(k) = sig(k) * factor
         else
            sig(k)  = 0.0d0
            sigt(k) = 0.0d0
         endif

c         if( k .eq. 26 )then
c            write(*,'(a20,1p8e12.3)')'RATED 1',sig(k)
c            stop'rated'
c         endif
ccccccccccccc

      enddo

c..   electron capture density*Ye factor
c..   applied to both ec and positron decay
      do k = k1, k2
         if( qval(k) .gt. 0.0d0 )then
c..   only use fkt rates, not ffnu here
c..   previous version had only first coef here (wda 3-13-08)
               expon = rcoef(1,k) + rcoef(2,k)/t9 + rcoef(3,k)/t913
     1              + rcoef(4,k)*t913 + rcoef(5,k)*t9
     2              + rcoef(6,k)*t953 + rcoef(7,k)*t9ln
            if( expon .gt. -300.0d0 )then
               sig(k)  = exp(expon)
               sigt(k) = 0.0d0
c               write(*,'(a10,i5,7a5,a5,1p9e12.3)')'RATED 3',k,
c     1              (rname(j,k),j=1,7),rlkh(k),sig(k),qval(k),
c     2              (rcoef(j,k),j=1,7)
c               stop'ggg'
cccccccc
            endif
         else
c..these are ffn reverse rates
            sig(k) = 0.0d0
            sigt(k)= 0.0d0
         endif
         sigv(k) = -sig(k) * rho
      enddo
      do k = k1, k2
         if( rlkh(k) .eq. '  ec' )then
c..density factors for electron captures
            sig(k)  = ye * rho * sig(k)
            sigt(k) = ye * rho * sig(k)
            sigv(k) = -sig(k) * rho
         endif
      enddo

c      write(*,'(a20,a5,1p8e12.3)')'RATED 2',rlkh(26),sig(26),ye,rho
cccccccc

c..   decks 2 to 8
      k1 = k1deck(2)
      k2 = k2deck(8)
      do k = k1, k2
         sigv(k)  =  0.0d0
         sig(k)   =  0.0d0
         sigt(k)  =  0.0d0
      enddo

c      k = 26
c      t9 = 0.01
c      do j = 1, 11
c         expon = rcoef(1,k) + rcoef(2,k)/t9 + rcoef(3,k)/t913
c     1        + rcoef(4,k)*t913 + rcoef(5,k)*t9
c     2        + rcoef(6,k)*t953 + rcoef(7,k)*t9ln
c      write(*,'(a20,i5,1p12e12.3)')'RATED',k,t9,sig(26),sig(135),
c     1        ye*rho*exp(expon),
c     1        expon, rcoef(1,k) , rcoef(2,k)/t9 , rcoef(3,k)/t913,
c     1        + rcoef(4,k)*t913 , rcoef(5,k)*t9,
c     2        + rcoef(6,k)*t953 , rcoef(7,k)*t9ln
c         t9 = t9 - 0.1
c      enddo
c      write(*,'(1p8e12.3)')(rcoef(j,k),j=1,7)
c      k = 135
c      t9 = 0.01
c      do j = 1, 11
c         expon = rcoef(1,k) + rcoef(2,k)/t9 + rcoef(3,k)/t913
c     1        + rcoef(4,k)*t913 + rcoef(5,k)*t9
c     2        + rcoef(6,k)*t953 + rcoef(7,k)*t9ln
c      write(*,'(a20,i5,1p12e12.3)')'RATED',k,t9,sig(26),sig(135),
c     1        ye*rho*exp(expon),
c     1        expon, rcoef(1,k) , rcoef(2,k)/t9 , rcoef(3,k)/t913,
c     1        + rcoef(4,k)*t913 , rcoef(5,k)*t9,
c     2        + rcoef(6,k)*t953 , rcoef(7,k)*t9ln
c         t9 = t9 - 0.1
c      enddo
c      write(*,'(1p8e12.3)')(rcoef(j,k),j=1,7)
cccccccccccc

      return
      end








