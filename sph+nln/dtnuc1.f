c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine dtnuc(dth0,tin,dth1,nflag)

c..   time step choice for subcycling
c..   dth1 is the updated value
c..   dth0 is the previously used value
c..   nflag = 0, before leqs, b(n) = dY(n)/dt
c..   else        after leqs, b(n) = dY(n)

      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'
      include 'cdtnuc'
      include 'cinput'

      real*8    dth0,fak,dth1,ydsum,aaa,taug,tin

      integer*4 n, nflag

      data delchi/1.0d-2/, chimin/1.0d-5/, fdtn/1.2d0/
      data fdysum/2.0d-3/
c-------------------------------------------------------------
c..DELCHI is maximum fractional change allowed in number density
c..CHIMIN is smallest y for delchi restriction
c..FDTN is the largest increase factor allowed for timestep dth
c..NUCLEU is the fastest changing nucleus in dth evaluation
c..TAUS is fdysum/sum(abs(Ai*dYi/dt))
c..FDYSUM is of order 0.01 or less, in practice
c..this keeps the matrix solution happy (dth < taus)
c-------------------------------------------------------------

c..limit fractional increase in time step: tau(1)
      tau(1) = fdtn * dth0
      nucleu = 0
c      itot = nnuc
     
c..limit fractional change in abundance:   tau(2)
      tau(2) = tau(1)
      taug   = tau(1)

      if( nflag .eq. 0 )then
c..   before leqs, so b(n) is dY(n)/dt
c..nuclei
         do  n = 1, itot-3
            if( y(n) .gt. chimin  )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
                  taug  =  dabs( y(n)/fak )*delchi
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo
         do  n = itot-2, itot
c            if( y(n) .gt. chimin*chimin  )then
            if( y(n) .gt. chimin )then
c..ignore neutrons on explicit guess (it has large errors)
               fak = b(n)
               if( fak .ne. 0.0d0 )then
                  taug  =  dabs( y(n)/fak )*delchi
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo

      else
c..   after leqs, so b(n) is dY(n)
         do  n = 1, itot-3
c..   guess at best time step
            if( y(n) .gt. chimin )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
c..delchi
                  taug  =  dabs( ( 1.0d-3 + y(n) )/fak )*delchi*dth0

               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            elseif( y(n) .gt. chimin**2 )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
c..delchi => 1
                  taug  =  dabs( ( 1.0d-3 + y(n) )/fak )*dth0
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo
         do  n =  itot-2,itot
c..   guess at best time step
            if( y(n) .gt. chimin )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
                  taug  =  dabs( ( 1.0d-3 + y(n) )/fak )*delchi*dth0
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
c            elseif( y(n) .gt. 1.0d-15 )then
            elseif( y(n) .gt. 1.0d-6 )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
c..delchi = 1 for very low abundances (try to avoid negative y)
                  taug  =  dabs( ( y(n) )/fak )*dth0
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo

      endif

c..limit fractional change in identity of nucleons:   tau(3)
c..limit total fractional number change to fdysum per step
      ydsum = 0
      do n = 1,itot
        aaa   = dble( nz(n) + nn(n) )
        ydsum = ydsum + abs( b(n)     )*aaa
      enddo

      if( abs(ydsum) .ne. 0.0d0 )then
        taus = fdysum/ydsum*dth0
      else
        taus = fdtn*dth0
      endif
      tau(3) = taus

c..decide on new value
      dth1 = tau(1)
      ndth = 1
      if( tau(2) .lt. dth1 )then
         dth1 =  tau(2)
         ndth = 2
      endif

c      if(tin .gt. 9.0d9)then 
c         if(dth .lt. 1.0d-10)dth=1.0d-10
c         if(dth .gt. 1.0d-5)dth=1.0d-5
c      elseif(tin .gt. 5.0d9)then
c         if(dth .lt. 1.0d-8)dth=1.0d-10
c      endif

c      if(tin.gt.1.0d9)write(*,*)'dtnuc: found new dth =',dth1
      return
      end

