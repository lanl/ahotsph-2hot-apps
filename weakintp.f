      subroutine weakintp(t9,rho,rloc,expon,exponu,tdex,rhodex,k)

c..   read netweak and perform interpolation in t9, rho

      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'

      integer*4 rloc, i, j, k, tdex, rhodex
c..   defined in dimenfile      parameter(tsize = 13, rhosize = 11)
      real*8 t9array(tsize), rhoarray(rhosize),  
     1       t9, rho, expon, exponu
      real*8 rcoeff, rholim
      

      data t9array/0.01d0,0.10d0,0.20d0,0.40d0,0.70d0,
     1     1.0d0,1.5d0,2.0d0,3.0d0,5.0d0,10.0d0,30.0d0,100.0d0/
      data rhoarray/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,
     1     7.0d0,8.0d0,9.0d0,10.0d0,11.0d0/
      data rholim/1.0d0/
c-----------------------------------------------------------------
      if( rho .lt. rholim)rho = rholim

c..use indices from rate.f

      call bilenshort(t9array,rhoarray,ratarray,t9,rho,expon,
     1     tsize,rhosize,wkreac,tdex,rhodex,rloc)

      call bilenshort(t9array,rhoarray,enuarray,t9,rho,exponu,
     1     tsize,rhosize,wkreac,tdex,rhodex,rloc)

      return
      end

      
