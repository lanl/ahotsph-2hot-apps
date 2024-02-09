      subroutine solven(dtstar,temp, rho, xin, deltah, irank, partid)

c..   evaluates rates, sets up matrix equations,
c..   and solves reaction network over time interval dth,
c..   updating mole fractions y and nucleon fractions x
c..   Takes timestep, temperature, density, and mole fractions
c..   xin in a 1d array as input. Outputs new y and energy generated
c..   (deltah)

      implicit none

      include 'dimen'

      include 'crate'
      include 'comcsolve'
      include 'caeps'
      include 'cdtnuc'
      include 'cdeuter'
      include 'cburn'
      include 'cconst'

      integer*4 irank
      integer*4 iDbug
      integer*4 partid

      real*8  yold(ndim), xin(ndim)
c..   instantaneous rate of energy generation
      real*8  ssi(kdm)

      real*8 temp
      real*8 tiny, t9, rho, dtstar, dth
      real*8 time, dtleft, dth0, rdt, sumd, suma, eb, dn, en, deltah
      real*8 fak, dne, yneu, ypro, xalf,  xnucleu
      real*8 dxnucl
      real*8 yestim,worsty
      real*8 enc,yoldsum,ynewsum,yesum, t9old(kdm)

      integer*4 inegs
      integer*4 idebug, ncycle, k, j, negflag, leq, n
      integer*4 it, kk, kreac

      integer*4 indx(ndim), i
      real*8 foo,alud(ndim,ndim)
cccccccccccccccccc


      integer*4 k1,k2,i1,i2
      integer*4 ncytest,ncymax
      integer*4 icall

c..   sparse solver variables
      integer*4 nelem,iflag
      parameter(nelem = nreac)
      real*8 berr
      
      parameter(ncytest = 3000, ncymax = 10000)

      save
c---------------------------------------------------------------
c..   floor on abundance variable used in network
      parameter( tiny=1.0d-100 )

c..   input: idebug,t9,rho,dtstar,y
c..   output: dth,yold,ncycle,y, deltah
c     aeps,y               (via commons caeps, comcsolve)
c     aepst,aepsv

c-------------------------------------------------------------
      iDbug = 1
c      if(temp.gt.1.d0) iDbug = 1
c-------------------------------------------------------------

      dth0   = 0.0d0
      inegs  = 0

c---------------------------------------------------------------
c      if(irank.eq.0 .and. iDbug) write(*,'(a,$)')
c     1   'solven, received T:'
c-------------------------------------------------------------
c..   full network solution with abundance update (implicit), or
c..   energy generation rate from right hand side evaluation (explicit)

      jnb = 0
      k = 1
      ncycle = 0
      ncytot = 0
c..   T,rho passed through call
      t9     = temp*1.0d-9
c      if(irank.eq.0 .and. iDbug) write(*,*)'T=',temp

c..   xin passed through call

      do n =1, nnuc
         y(n) = xin(n) / dble( nz(n) + nn(n))
c         if(irank.eq.0 .and. iDbug) write(*,*)y(n),nz(n),nn(n)
      enddo

c     if1: t9 sufficient for reactions
      if( (t9 .ge. 1.0d-2) .and. (t9 .lt. tnse) ) then
c     1    .and. (rho .gt. 1.0d3))then
c..   full network
         do j = 1, itot
            b(j)    = 0.0d0
            bt(j)   = 0.0d0
            bv(j)   = 0.0d0
            yold(j) = y(j)
            ydcon(j) = axdcon(j,k)
         enddo
         
c..   evaluate rates and put in sig array

         call rate(t9,rho,dtstar)

c..   al26 rate: fix at low T
c..   reaclib fit is bad below 2.9d7 K
c         if( t9 .lt. 0.029d0 )then
c            sig(lkal26) = 3.051d-14 
c         endif
c.. commented out by CIE - lkal26 not defined anywhere
cccccccccccccc


c..   weak, intermediate, and strong screening
         call screen(t9,rho)

         time   = 0.0d0
         dtleft = dtstar
c..   most optimistic choice for time step
         dth    = dtleft
         deltah = 0.0d0
c..   ncycle: subcycle over shorter timesteps for network if needed
 100     continue

         ncycle = ncycle + 1

         if( ncycle .gt. ncymax )goto 101

c..   insure time interval does not overshoot
         if( dth .gt. dtleft )then
            dth = dtleft
         endif

         if( dth .le. 0.0d0 )then
            write(*,'(7a12)')
     1           'dtstar','dth','dth0','dtleft','T9','rho','ncycle','k'
            write(*,'(1p6e12.5,2i12)')dtstar,dth,dth0,dtleft,t9,rho,ncycle,k
            write(*,*)'particle',partid
            do n =1, nnuc
               if(iDbug.eq.0) write(*,*)y(n),nz(n),nn(n)
            enddo
            stop'solven dth'
         endif

         dth0   = dth

c..   setup matrices for linear equation solver
c..   need rates, abundances Y, partition functions (2J+1=w factors)
         call rhside(idebug)

c..   continue with network solution
         call jacob

c..   update diagonal terms
         rdt = 1.0d0/dth
         do j = 1, itot
            a(j,j) = a(j,j) + rdt
         enddo
         
c..   SPARSE SOLVER SECTION
c..   put nonzero entries of A matrix into sparse_dfdy vector
         do j=1,nlinks
            sparse_dfdy(j) = a(iloc(j),jloc(j))
         enddo

c..   sparseu controls bias towards sparsity or numerical pivoting
c..   see sparse_ma28.f
         sparseu = 0.1d0

         call ma28bd(nnuc,nlinks,sparse_dfdy,nelem,iloc,jloc,jvect,
     1               ikeep,iw,sparsew,iflag)
         if (iflag .lt. 0) then
            write(*,*) 'error in ma28bd flag',iflag
            stop 'error in ma28bd in solven'
        endif

        call ma28cd(nnuc,sparse_dfdy,nelem,jvect,ikeep,b,
     1               sparsew,1)
         

c..   evaluate next time step
         call dtnuc(dth0,t9*1.d9,dth,1)
c         if(t9.ge.1.d-1) write(*,*)'new dth=',dth

c..   check result for excessive changes
         sumd = 0
         suma = 0
         sumx = 0
         sumz = 0
         do  j = 1, itot
            sumd = sumd +         b(j)  * dble( nn(j) + nz(j) )
            suma = suma +     abs(b(j)) * dble( nn(j) + nz(j) )
            sumx = sumx + (y(j) + b(j)) * dble( nz(j) + nn(j) )
         enddo

c..   catches worst negative abundance
         negflag = 0
         worsty  = 0.0d0
         do j = 1, itot
            yestim = y(j) + b(j)
            if( yestim .lt. -3.0d-16 )then
               negflag = j
               worsty = yestim
            endif
         enddo
         
c..   adjust time step
         if( negflag .gt. 0 )then
            if( b(negflag) .lt. 0.0d0 )then

c..   would give 0.5 of initial value if linear decline
               if( y(negflag) .lt. 0.0d0 )then
c..   avoid zero abundance
                  dth = -0.5d0*dth0*y(negflag)/b(negflag)
               endif

            else
               write(*,*)'solven error: negflag',
     1              negflag,b(negflag)
               stop' solven negflag'
            endif
         endif

         if(  abs(sumd) .lt. 1.0d-5
     1        .and. negflag .eq. 0 
     2        .and. dth .ge. dth0*0.5d0 )then

c..   tests were successful
c..   so update particle densities, time elapsed
            sumx = 0
            sumz = 0
            eta  = 0
            do  j = 1, itot
               y(j) = y(j) + b(j)
c..   positive or zero
               y(j) = dmax1( y(j), 0.0d0 )
               x(j) = y(j) * dble( nz(j)+nn(j) )
               sumx = sumx + x(j)
               sumz = sumz + y(j)* dble( nz(j) )
               eta  = eta  + y(j)* dble(nn(j) - nz(j))

            enddo
c..   Brand new thesis paranoia
c..   renormalize abundances
         sumx = -1.0d0
         do j = 1, nnuc
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
         enddo

         do j = 1, ndim-1
            x(j) = x(j)/( 1.0d0 + sumx )
            y(j) = x(j)/dble( nz(j) + nn(j) )
         enddo
c..   reset Ye
         yesum = 0.0d0
         sumz   = -1.0d0
         do j = 1, ndim-1
            yesum = yesum + y(j) * dble( nz(j) )
            sumz = sumz + x(j)
         enddo
         x(ndim) = yesum
         y(ndim) = yesum


c..   elapsed time
            time     = time + dth0
            dtleft   = dtstar - time
            checksum = sumx - 1.0d0

c..   output if excessive subcycles
            if( ncycle .ge. ncytest )then
               if( nucleu .gt. 0 .and. nucleu .le. itot )then

                  if( ncycle .eq. ncytest )
     1                 write(*,'(/2a4,a5,a11,10a9)')
     1                 'k','ncy','nuc','y','tau1','tau2','tau3',
     1                 'time','dth0','dth','dtleft'
     2                 ,'suma','sumd','check'
                  write(*,'(2i4,a5,1pe11.4,1p10e9.2)')
     1                 k,ncycle,xid(nucleu),y(nucleu),tau,time,
     1                 dth0,dth,dtleft
     2                 ,suma,sumd,checksum
                  write(3,'(2i4,a5,1pe11.4,1p10e9.2)')
     1                 k,ncycle,xid(nucleu),y(nucleu),tau,time,
     1                 dth0,dth,dtleft
     2                 ,suma,sumd,checksum
               else
                  write(*,'(2i6,12x,1p10e11.3)')
     1                 k,nucleu,tau,time,dth0,dth,dtleft
     2                 ,suma,sumd
                  write(3,'(2i6,12x,1p10e11.3)')
     1                 k,nucleu,tau,time,dth0,dth,dtleft
     2                 ,suma,sumd
               endif
            endif

c..   reset Ye (copied from NSE-loop. ~CIE)
         yesum = 0.0d0
         sumz   = -1.0d0
         do j = 1, ndim-1
            yesum = yesum + y(j) * dble( nz(j) )
            sumz = sumz + x(j)
         enddo
         x(ndim) = yesum
         y(ndim) = yesum

c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
            eb   = 0.0d0
            dn   = 0.0d0
            en   = 0.0d0
            do j = 1, itot
               fak =  y(j) - yold(j)
               eb  = eb - qq(j)*fak
               dn  = dn +       fak
               k1 = k1deck(1)
               k2 = k2deck(2)
               do kreac = k1, k2
                  if(nrr(1,kreac) .eq. j)then
                     en = en - avagadro*1.6021d-6*signue(kreac)*fak
                  endif
               enddo
            enddo
            eb      = eb *   9.65d+17    /(dtstar-dtleft)
            dne     = dn * 1.2476d+17 *t9/(dtstar-dtleft)
            en      = en / (dtstar-dtleft)
            aeps = eb - dne + en
            deltah = deltah + aeps*(dtstar-dtleft)

         else
c..   tests were unsuccessful
            dth = dmin1( dth, 0.5d0*dth0 )
            if( ncycle .ge. ncytest )then
               if( ncycle .eq. ncytest )
     1              write(*,'(a7,a10,a4,a10,
     1              3a5,7a10,a5)')' ','dth0',
     1              ' ','dth', 'k','ncy','xid','X','dX','dtlef',
     1              'yp','xa','suma','sumd','negf'
               yneu = y(itot-2) + b(itot-2)
               ypro = y(itot-1) + b(itot-1)
               xalf = 4.0d0*( y(itot) + b(itot) )
               if( nucleu .le. 0 .or. nucleu .gt. itot )then
                  write(*,'(2(a7,1pe9.1),3i5,1p7e10.2,i5)')
     1                 'reject',dth0,' try',dth, k,ncycle,nucleu, 
     2                 t9,dtleft,yneu,ypro,xalf,suma,sumd,negflag 
               else
                  xnucleu = y(nucleu)*
     1                 dble( nz(nucleu) + nn(nucleu))
                  dxnucl  = b(nucleu)*
     1                 dble( nz(nucleu) + nn(nucleu))
                  write(*,'(a7,1pe10.3,a4,1pe10.3,
     1                 2i5,a5,1p8e10.2,i5)')
     1                 'reject',dth0,' try',dth, k,ncycle,
     1                 xid(nucleu),xnucleu,dxnucl,dtleft,
     2                 yneu,ypro,xalf,suma,sumd,negflag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  if( negflag .gt. 0 )
     1                 write(*,*)xid(negflag),y(negflag),
     1                 b(negflag),y(negflag)+b(negflag)

                  do j = 1, nreac
                     if( nrr(1,j) .eq. negflag .and.
     1                    sig(j) .gt. 0.0d0 )then
                        write(*,'(4i5,2a5,1p8e12.3)')negflag,nrr(1,j),
     1                       j,ideck(j),xid(nrr(1,j)),
     1                       xid(nrr(2,j)),
     1                       sig(j),sig(j)*dth,y(nrr(1,j))
                     endif
                  enddo
                  write(3,'(a7,1pe10.3,a4,1pe10.3,2i5,a5,
     1                 1p7e10.2,i5)')
     1                 'reject',dth0,' try',dth, k,ncycle,
     1                 xid(nucleu),xnucleu,dxnucl,dtleft,
     2                 ypro,xalf,suma,sumd,negflag
                  write(*,*)'particle',partid
                  do n =1, nnuc
                     if(iDbug.eq.0) write(*,*)y(n),nz(n),nn(n)
                  enddo
                  stop'solven: no success'
               endif
            endif
         endif
         
         if( dtleft .gt. 1.0d-14*dtstar )goto 100

         if( ncycle .gt. ncytest )then
            yneu = y(itot-2) + b(itot-2)
            ypro = y(itot-1) + b(itot-1)
            xalf = 4.0d0*( y(itot) + b(itot) )
            if( nucleu .le. 0 .or. nucleu .gt. itot )then
               write(*,'(2(a7,1pe9.1),3i5,1p7e10.2,i5)')
     1              'last',dth0,' try',dth, k,ncycle,nucleu,
     2              t9,dtleft,yneu,ypro,xalf,suma,sumd,negflag
            else
               xnucleu = y(nucleu)*dble( nz(nucleu) + nn(nucleu)) 
               dxnucl  = b(nucleu)*dble( nz(nucleu) + nn(nucleu)) 
               write(*,'(2(a7,1pe11.3),2i5,a5,1p7e10.2,i5)')
     1              'last',dth0,' try',dth, k,ncycle,xid(nucleu), 
     2              xnucleu,dxnucl,dtleft,ypro,xalf,suma,sumd,
     3              negflag
               write(3,'(a5,1pe10.2,a4,1pe10.2,2i5,a5,1p7e10.2,
     1              i5)')
     1              'last',dth0,' try',dth, k,ncycle,xid(nucleu), 
     2              xnucleu,dxnucl,dtleft,ypro,xalf,suma,sumd,
     3              negflag
            endif
         endif


 101  continue
         if(ncycle .gt. ncymax)then
           write(*,*)'more than ',ncymax,' subcycles'
           write(*,*)ncycle,' subcycles in zone ',k
           write(*,'(a5,2a20,6a12)')'k','t9','rho','dtstar','dtleft',
     1       'dth','checksum','eta','sumd'
           write(*,'(i5,1p2e20.13,1p6e12.3)') k,t9,rho,dtstar,dtleft,
     1       dth,checksum,eta,sumd
           if( nucleu .gt. 0 .and. nucleu .le. itot )then
              write(*,*)'fastest nucleus is ',xid(nucleu)
           endif

           write(*,*)'particle',partid
           do n =1, nnuc
              if(iDbug.eq.0) write(*,*)y(n),nz(n),nn(n)
           enddo
           stop'solven too many subcycles'
         endif
c     end if1?     

c..   total matrix solutions
         ncytot = ncytot + ncycle


c     if2: t9 sufficient for NSE
      elseif( (t9 .ge. tnse) ) then
c     1       .and. (rho .gt. 1.0d2))then
c..   itbrn(i) = 3 ...............................................
c..   nse solver plus weak interactions......................
         
         ncycle = 0

         do j = 1, itot
            b(j)    = 0.0d0
            yold(j) = y(j)
         enddo

c..   evaluate rates and put in sig array
         call rate(t9,rho,dtstar)

         call screen(t9,rho)

cccccccccccccccccccccc

c..   make sure Ye is set
      if( y(ndim) .ne. y(ndim) ) then
         yesum = 0.0d0
         do j = 1, ndim-1
            yesum = yesum + y(j) * dble( nz(j) )
         enddo
         x(ndim) = yesum
         y(ndim) = yesum
      endif

         call nse(t9,rho,y,enc)

c..   most optimistic choice for time step
         dth    = dtstar

c..   deck 1: i ---> j
c..   beta decays, positron decays,
c..   i + e ---> j for electron captures (see rate.f)
         weaksum = 0.0d0
         k1 = k1deck(1)
         k2 = k2deck(1)
         if( k1 .ne. 0 .and. k2 .ne. 0 )then
            do j = k1,k2
               i1 = nrr(1,j)
               i2 = nrr(2,j)
               b(i1)    = b(i1)    - y(i1)*sig(j)
               b(i2)    = b(i2)    + y(i1)*sig(j)
               weaksum = weaksum + ( nz(i2)-nz(i1) )*y(i1)*sig(j)
c..   dYe/dt is weaksum
               bt(i1)    = bt(i1)    - y(i1)*sigt(j)
               bv(i1)    = bv(i1)    - y(i1)*sigv(j)
               bt(i2)    = bt(i2)    + y(i1)*sigt(j)
               bv(i2)    = bv(i2)    + y(i1)*sigv(j)
            enddo
         endif

c..   setup matrices for linear equation solver
c..   need rates, abundances Y, partition functions (2J+1=w factors)
         call rhside(idebug)

c..   continue with network solution
         call jacob

c..   update diagonal terms
         rdt = 1.0d0/dth
         do j = 1, itot
            a(j,j) = a(j,j) + rdt
         enddo

c..   leqs solves linear equations
         call leqs(a,b,itot,nnuc)

c..   total matrix solutions done
         ncytot = ncytot + 1

c..   update particle densities, time elapsed
         sumx = 0
         sumz = 0
         eta  = 0
         do  j = 1, itot
            y(j) = y(j) + b(j)
c..   positive or zero
            y(j) = dmax1( y(j), 0.0d0 )
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
            sumz = sumz + y(j)* dble( nz(j) )
            eta  = eta  + y(j)* dble(nn(j) - nz(j))
         enddo
         checksum = sumx - 1.0d0
c         write(*,*)'checksum in nse ',checksum
c..   Brand new thesis paranoia
c..   renormalize abundances
         sumx = -1.0d0
         do j = 1, nnuc
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
         enddo

         do j = 1, ndim-1
            x(j) = x(j)/( 1.0d0 + sumx )
            y(j) = x(j)/dble( nz(j) + nn(j) )
         enddo
c..   reset Ye
         yesum = 0.0d0
         sumz   = -1.0d0
         do j = 1, ndim-1
            yesum = yesum + y(j) * dble( nz(j) )
            sumz = sumz + x(j)
         enddo
         x(ndim) = yesum
         y(ndim) = yesum

c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
         eb   = 0
         dn   = 0
         yoldsum = 0.0d0
         ynewsum = 0.0d0
         do j = 1, itot
            yoldsum = yoldsum+yold(j)
            ynewsum = ynewsum+y(j)
            fak =  y(j) - yold(j)
            eb  = eb - qq(j)*fak
            dn  = dn +       fak
         enddo

         eb      = eb *   9.65d+17   /dtstar
         dne      = dn * 1.2476d+17 *t9/dtstar
         aeps = (eb - dne)/rho
         deltah = aeps*dtstar


c     end if2?


c     if3: t9 not sufficient for any burning, just decays
      elseif( (t9 .lt. 1.0d-2) ) then
c             .or. (rho .le. 1.0d1))then
c..   itbrn(i) = 0 ...............................................
c..   network solution for decays only......................
         ncycle = 0
c         rho = den
         do j = 1, itot
c            y(j)    = aex(j)
            b(j)    = 0.0d0
            yold(j) = y(j)
         enddo

c..   evaluate rates and put in sig array
         call rated(t9,rho)

c..   most optimistic choice for time step
         dth    = dtstar

c..   deck 1: i ---> j
c..   beta decays, positron decays,
c..   i + e ---> j for electron captures (see rate.f)
         weaksum = 0.0d0
         k1 = k1deck(1)
         k2 = k2deck(1)
         if( k1 .ne. 0 .and. k2 .ne. 0 )then
            do j = k1,k2
               i1 = nrr(1,j)
               i2 = nrr(2,j)
               b(i1)    = b(i1)    - y(i1)*sig(j)
               b(i2)    = b(i2)    + y(i1)*sig(j)
               weaksum = weaksum + ( nz(i2)-nz(i1) )*y(i1)*sig(j)
c..   dYe/dt is weaksum
               bt(i1)    = bt(i1)    - y(i1)*sigt(j)
               bv(i1)    = bv(i1)    - y(i1)*sigv(j)
               bt(i2)    = bt(i2)    + y(i1)*sigt(j)
               bv(i2)    = bv(i2)    + y(i1)*sigv(j)
            enddo
         endif

c..   setup matrices for linear equation solver
c..   need rates, abundances Y, partition functions (2J+1=w factors)
         call rhside(idebug)

c..   continue with network solution
         call jacob

c..   update diagonal terms
         rdt = 1.0d0/dth
         do j = 1, itot
            a(j,j) = a(j,j) + rdt
         enddo

c..   leqs solves linear equations
         call leqs(a,b,itot,nnuc)

c..   total matrix solutions done
         ncytot = ncytot + 1

c..   update particle densities, time elapsed
         sumx = 0
         sumz = 0
         eta  = 0
         do  j = 1, itot
            y(j) = y(j) + b(j)
c..   positive or zero
            y(j) = dmax1( y(j), 0.0d0 )
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
            sumz = sumz + y(j)* dble( nz(j) )
            eta  = eta  + y(j)* dble(nn(j) - nz(j))
         enddo
         checksum = sumx - 1.0d0

c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
         eb   = 0
         dn   = 0
         do j = 1, itot
            fak =  y(j) - yold(j)
            eb  = eb - qq(j)*fak
            dn  = dn +       fak
         enddo
         eb      = eb *   9.65d+17    /dtstar
         dne      = dn * 1.2476d+17 *t9/dtstar
         aeps = eb - dne
         deltah = aeps*dtstar

c     end if3
 
c     just in case
      else
         deltah = 0.0

      endif



      icall = 1
c      if(irank.eq.0) write(*,*)'deltah=',deltah

c     update xin to be passed back to calling routine. ~CIE
c     I'm guessing in original that was unnecessary since global array
      do n =1, ndim
         xin(n) = x(n) 
      enddo

      return

      end


