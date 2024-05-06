      subroutine nse(t9,rho,yeq,enc)

c..   NUCLEAR STATISTICAL EQUILIBRIUM (modified: wda 4/2/00)

c..   input
c     itot = total number of n,p,alpha,nuclei
c     inuc = itot-3 = number of "nuclei"
c     idebug
c     t9,rho,eta
c..   output
c     niter = number of iterations used
c     uaaa  = 2*( mu(n) + mu(p) ) = mu(alpha) in nse
c     uhat  =     mu(n) - mu(p)
c     ierr  = flag for nonconvergent return

      implicit none

      include 'dimen'
      include 'crate'
      include 'cburn'
      include 'comcsolve'
      include 'caeps'
      include 'ceoset'
      include 'cconst'

      real*8      u(nnuc), yeq(ndim), xeq(nnuc), qeq(nnuc)
      character*5 plus, goes, blank

      integer*4 nloop,k,ierr,niter
      real*8    delmx, tiny
      parameter( nloop = 500, delmx = 0.2d0, tiny = 1.0d-13 )

c..   niso = number of isotopes in FKT list
c      integer*4 niso
c      parameter(niso = 806) 
c      real*8      p(7)
c      character*4 lkh                                     
c      character*5 nam(6),inam(niso),rnam(6)                      
c      character*1 vw,nr

c..   cxid = character id for isotope for summary i/o
c..   xxeq = nucleon fraction in nse
      character*5 cxid(nnuc)
      real*8      xxeq(nnuc)

      real*8 bpern,fact,rho9,t932,t9,tin,rho,theta,tk,uhlow,ya,xal
      real*8 uahi,omxa,uhhi,uhat,uaaa,uaaa0,aa,zz,arat,ueff,xeqm,yense
      real*8 dxdmu,dedmu,etaq,dxdmh,dedmh,dudmu,dydmu,etai
      real*8 dudmh,dydmh,det,errx,yeeq,errz,anum,hnum,ddu,delua,ddh
      real*8 deluh,dfracu,dfrach,dfrac,df0,reduct,enc,pnc,gam4
      real*8 xpp,xnn,xaa,ynn,ypp,xxl,uplow,xp,yp,uphi,adelua
      real*8 ualow,uhat0
      real*8 up,un,yheold,ww(nnuc),yeqold(ndim)
      real*8 chemfak, chemcon
      real*8 dlog_ww_p, dlog_ww_n, dlog_ww_a
      integer*4 ib, in,ip,ia,ni56,nc12,nfe54,nfe56,nucn,nucpr
      integer*4 nuc0, i, n, iitot

      data blank/'     '/
      data plus /'  +  '/,goes/' <-> '/

      integer*4 icall, iachain
      data icall/0/,iachain/1/

c..   icall flags first call for initiation
c..   iachain flags lack of N .ne. Z nuclei
c..   ib is index of most bound nucleus
c..   bpern is binding energy of most bound nucleus
c..   ni56, nc12 are indices for special nuclei
      save icall,iachain,ib,bpern
      save in,ip,ia,ni56,nc12
c----------------------------------------------------

c      write(*,*)'NSE, got t9', t9

      ierr = 0
      ye = yeq(ndim)
      if( icall .eq. 0 )then
c..   only on initial call
         icall = 1
c..   determine if network is alpha-chain
         nucn = 0
         nucpr = 0
         nuc0 = 0
         do i = 1, inuc
            if( nn(i) .gt. nz(i) )then
               nucn = nucn + 1
            elseif( nn(i) .lt. nz(i) )then
               nucpr = nucpr + 1
            else
               nuc0 = nuc0 + 1
            endif
         enddo
c         write(*,'(4(a15,i5))')"n-rich nuclei",nucn,
c     1        "p-rich nuclei",nucp,
c     2        "Z=N nuclei",nuc0,
c     3        "total nuclei",inuc

         if( inuc .eq. nuc0 )then
c..   true, it is alpha chain
            iachain = 0
            open(33,file='data.a')
         else
c..   false, it is general network
            iachain = 1
         endif

c..   define "particles"
         ia   = itot
         ip   = itot - 1
         in   = itot - 2

c..   find most bound nucleus in network
c..   start with He4
         bpern = (-qq(ia) + qq(ip)*2.0 + qq(in)*2.0)/4.0d0
         do i = 1, itot
            fact = (-qq(i) + qq(ip)*float( nz(i) ) 
     1           + qq(in)*float( nn(i) ) )
     2           /float( nz(i) + nn(i) ) 
c..   remember index of special nuclei
            if( nz(i) .eq. 28 .and. nn(i) .eq. 28 )ni56 = i
            if( nz(i) .eq.  6 .and. nn(i) .eq.  6 )nc12 = i
            if( fact .gt. bpern )then
               bpern = fact
               ib = i
            endif
         enddo
c         write(*,'(i5,a5,f9.4,a60/)')ib,xid(ib),bpern,
c     1        " excess Mev/nucleon; most bound nucleus in network"
      endif

c..   executed once on each call
      do i=1,itot
         ww(i) = w(i) * pf(i)
      enddo
      rho9  = rho * 1.0d-9
      t932  = t9 * sqrt( t9 )
      theta = 0.1013d0 * rho9 /t932
      tk    = t9/11.605d0
      tin   = t9*1.0d9
      yense    = ye(1)
      eta = 1.0d0 - 2.0d0*ye(1)
c----------------------------------
      if( iachain .ne. 0 )then
c----------------------------------
         if( ni56 .eq. 0 )then
            write(*,*)'error in qse, no Ni56 in network'
            stop
         endif
c         nfe54  = 0
c         do i = 1, itot
c            yeqold(i) = yeq(i)
c            if( nz(i) .eq. 26 .and.  nn(i) .eq. 28 )then
c               nfe54 = i
c            endif
c         enddo
c         if( nfe54 .eq. 0 )then
c            write(*,*)'error in qse, no Fe54 in network'
c            stop
c         endif
         nfe56  = 0
         do i = 1, itot
            yeqold(i) = yeq(i)
            if( nz(i) .eq. 26 .and.  nn(i) .eq. 30 )then
               nfe56 = i
            endif
         enddo
         if( nfe56 .eq. 0 )then
            write(*,*)'error in qse, no Fe54 in network'
            stop
         endif
c..   NUCLEAR STATISTICAL EQUILIBRIUM
c..   try dissociated value for Yn/Yp for first uhat guess
c..   initial guess for mu(alpha) = u(ia)
c..   low T value = all ni56 approximation
         ualow = ( qq(ni56) + tk*dlog( theta/56.0d0**2.5 ) )
     1        *   4.0d0 / ( 56.0d0 )
         ya    = 8.0d0/theta * exp( (ualow - qq(itot) )/tk  )
c..   add saturation at high T
         ya    = ya/(1.0 + 4.0*ya)
         xal    = 4.0d0*ya
c..   high T value = all He4, n, p approximation
         uahi = qq(itot) + tk*dlog( theta/32.0d0 )
c..   compromise value
         chemfak =  dlog( rho/tin**1.5d0 )
         chemcon = dlog(
     1        avagadro*(planck**2*avagadro/(2.0d0*pi*boltz))**1.5d0 )
         u(ia) = ualow*(1.0d0 - xal) + uahi*xal
c         u(ia) = qq(nnuc) + tk*(etanuc(nnuc,k)-dlog(ww(nnuc)))
c cie: due to a bug in masses.f (somthing stomping over the memory of 
c      nz/n when or after those are read in), the last three elements
c      in 'ww' can be zero, which gives an 'inf' when taking the log
         if (ww (nnuc-1) .gt. 0.0) then
             dlog_ww_a = dlog (ww (nnuc))
         else 
             dlog_ww_a = 0.0
         endif
         u(ia) = qq(nnuc) + tk*(dlog(yeq(nnuc))+chemfak+chemcon
     1   - 1.5d0*dlog(xa(nnuc)) - dlog_ww_a)
c..   initial guess for muhat = mu(n) - mu(p)
c..   high T value = free n have all neutron excess
         omxa  = 1.0d0
         uhhi  = ( qq(in) - qq(ip) )
     1        + tk*dlog( (omxa + eta)/(omxa - eta) )
c..   low T value = Fe54 has all neutron excess
c         uhlow = qq(nfe54) + tk*dlog( eta*theta*0.5d0/54.0d0**1.5 )
c     1        - 54.0d0/4.0d0*ualow
c..   Small SNSPH network has no Fe54. Use Fe56 instead 
         uhlow = qq(nfe56) + tk*dlog( eta*theta*0.5d0/56.0d0**1.5 )
     1        - 56.0d0/4.0d0*ualow
c..   compromise value
         uhat  = uhlow*(1.0d0 - xal) + uhhi*xal
c..   use chemical potential values from state.f
c         up = qq(nnuc-1) + tk*(etanuc(nnuc-1,k)-dlog(ww(nnuc-1)))
c         un = qq(nnuc-2) + tk*(etanuc(nnuc-2,k)-dlog(ww(nnuc-2)))
c cie: due to a bug in masses.f (somthing stomping over the memory of 
c      nz/n when or after those are read in), the last three elements
c      in 'ww' can be zero, which gives an 'inf' when taking the log
         if (ww (nnuc-1) .gt. 0.0) then
             dlog_ww_p = dlog (ww (nnuc-1))
         else 
             dlog_ww_p = 0.0
         endif
         if (ww (nnuc-2) .gt. 0.0) then
             dlog_ww_n = dlog (ww (nnuc-2))
         else
             dlog_ww_n = 0.0
         endif
         up = qq(nnuc-1) + tk*(dlog(yeq(nnuc-1))+chemfak+chemcon
     1   - dlog(2.0d0) - dlog_ww_p)
         un = qq(nnuc-2) + tk*(dlog(yeq(nnuc-2))+chemfak+chemcon
     1   - dlog(2.0d0) - dlog_ww_n)
         uhat = (un - up)/1.0d0
         uaaa = 2.0*(un + up)/1.0d0
c..   uses u(ia) = u + m here, so it is "mu" not "u"
         uaaa = u(ia)
c..   save first guesses for convergence diagnostics
         uhat0 = uhat
         uaaa0 = uaaa

c..   interate for nucleon number (sum Xi = 1).....................
         do n = 1, nloop
c..   generate trial nse abundance values
            do i = 1, itot  
               aa     =  nz(i) + nn(i) 
               zz     =  nz(i) 
               arat   = aa**1.5d0 * ww(i) / theta
               ueff   = (aa*0.25d0     ) * uaaa
     1              + (aa*0.50d0 - zz) * uhat
     2              - qq(i)
               yeq(i) = arat * exp( ueff / tk )
               xeq(i) = yeq(i)*aa
            enddo
c..   generate Ye and sum Xi, and their derivatives for Newton-Raphson
            xeqm    = 0
            dxdmu   = 0
            dedmu   = 0
            etaq    = 0
            dxdmh   = 0
            dedmh   = 0
            do i = 1, itot  
c..   complete sum for all nucleons
               aa    = nz(i) + nn(i) 
               xeqm  = xeqm + aa * yeq(i)
               dudmu = aa * 0.25d0
               dydmu = yeq(i) / tk * dudmu
               dxdmu = dxdmu + aa*dydmu
               zz    =  nz(i) 
               if( nz(i) .ne. nn(i) )then
c..   Z .ne. N sum for neutron excesses
                  etai  = aa - 2.0d0*zz
                  dedmu = dedmu + etai*dydmu
                  etaq  = etaq + etai * yeq(i)
                  dudmh = aa*0.50d0 - zz
                  dydmh = yeq(i) / tk * dudmh
                  dxdmh = dxdmh + aa  *dydmh
                  dedmh = dedmh + etai*dydmh
               endif
            enddo 
            det    = dxdmu*dedmh - dxdmh*dedmu
c..   safety check for redundant equations:
c..   attempts to avoid unrealistically large corrections
            if( n .eq. 1 )then
               if( abs(det) .le. 1.0d-1 )then
                  det = dsign( 1.0d-1, det ) 
               endif
            endif
            errx   = xeqm - 1.0d0
            yeeq   = (1.0d0 - etaq)*0.5d0
            errz   = etaq - eta
            anum   = -errx*dedmh + errz*dxdmh
            hnum   =  errx*dedmu - errz*dxdmu
            ddu    = anum/det
            delua  = -errx/dxdmu
            ddh    = hnum/det
            deluh  = -errz/dedmh
c..   safety check for large changes:
c..   reduces increments simultaneously
c..   uses n*kT as measure of reasonable step size
            dfracu = abs( ddu )
            dfrach = abs( ddh )
            dfrac  = dmax1( dfracu, dfrach )
            df0    = 3.0d0*tk
            if( dfrac .gt. df0 )then
               reduct = df0/dfrac
               ddu = ddu * reduct
               ddh = ddh * reduct
            else
               reduct = 1.0d0
            endif
            uaaa  = uaaa + ddu
            uhat  = uhat + ddh

            if( abs(xeqm - 1.0d0) .lt. tiny )goto 998
         enddo

         write(*,'(a33,i4,4(a8,1pe11.3))')
     1        "NONCONVERGENCE in NSE: iteration ",n-1,
     2        " d(X)", xeqm-1,  " d(eta)", errz,
     3        " mu(a)",   uaaa,  " mu(hat)",  uhat
         ierr  = 1
         niter = n
         return

 998     continue
c...................................end of iteration logic......
c..   NSE values are now determined
         niter = n
c..   uses mass excesses for energy calculation....................
c..   set zero energy to cold, pure composition of
c..   the most bound nucleus in network
c     enc =  ( - qq(ib)/float(nz(ib)+nn(ib)) )
c     1    * 9.65d+17
c..   energy relative to C12 nuclei
         enc  = 0
         pnc  = 0
         yeeq = 0
         do i = 1, itot
            enc  = enc  + (yeq(i)-yeqold(i))*(1.2476d+17*t9 -
     1          qq(i)*9.65d+17 )
            pnc  = pnc  + yeq(i)* 1.2476d+17*t9 * 1.5d0
            yeeq = yeeq + float( nz(i) )*yeq(i)
         enddo
         pnc  = pnc * rho
         gam4 = pnc/(enc * rho) + 1.0d0
c         write(*,*)'gam4',gam4,k
         yeq(ndim) = yeeq

c..   sum n-rich,alpha,p-rich isotopes
         xpp = 0
         xnn = 0
         xaa = 0
         ynn = 0
         ypp = 0
         do i = 1, itot
            aa = nz(i) + nn(i)
            zz = nz(i)
            b(i) = yeqold(i) - yeq(i)
            if( nz(i) .eq. nn(i) )then
               xaa = xaa + xeq(i)
            elseif( nz(i) .lt. nn(i) )then
               xnn = xnn + xeq(i)
               ynn = ynn + yeq(i)* float( nn(i) - nz(i) )
            else
               xpp = xpp + xeq(i)
               ypp = ypp + yeq(i)* float( nz(i) - nn(i) )
            endif
         enddo
c..   most abundant nuclei (above xxl)
         iitot = 0
         xxl   = 1.0d-1
         do i = 1, itot-3
            if( xeq(i) .gt. xxl )then
               iitot        = iitot + 1
               cxid(iitot)  = xid(i)
               xxeq(iitot)  = xeq(i)
            endif
         enddo

c         if( iitot .gt. 0 )
c     1        write(*,'(8(a5,1pe10.2))')(cxid(i),xxeq(i),i=1,iitot)

c..   n,p,alphas
c         write(*,'(3(a5,1pe10.2))')(xid(i),xeq(i),i=itot-2,itot)
c         write(*,*)' '
c         write(33,'(2i5,1p13e11.3)')ierr,niter,t9,rho,eta,
c     1        uaaa,uaaa0,uhat,uhat0,xeq(itot-1),xeq(itot-2),xeq(itot)
c     2        ,tk,enc*1.0d-18,pnc/rho*1.0d-18



c         write(*,'(i5,0pf8.3,1p2e10.2,1pe12.4,1p6e10.2,a5,1pe10.2)')
c     1        n,t9,rho,eta,enc,xeq(in),xeq(ip),xeq(ia),uaaa,uhat,det,
c     2        cxid(1),xxeq(1)

         return

c-------------------------------------------------
      else
c-------------------------------------------------
         do i = 1, itot
            qeq(i) = qq(i)
         enddo
c..   represent n,p as "average nucleons" so Yp = Yn
         qeq(in) = 0.5*( qq(in) + qq(ip) )
         qeq(ip) = qeq(in)
c..   alpha chain equilibrium
c..   initial guess for mu(alpha) = u(ia)  
c..   low T value = all ni56 approximation
         ualow = ( qeq(ni56) + tk*dlog( theta/56.0d0**2.5 ) )
     1        *   4.0d0 / ( 56.0d0 )
         ya    = 8.0d0/theta * exp( (ualow - qeq(itot) )/tk  )
c..   add saturation at high T
c     ya    = ya/(1.0 + 4.0*ya)
         ya    = dmin1(ya, 0.25d0)
         xal    = 4.0d0*ya
c..   higher T value = all He4, n, p approximation
         uahi = qeq(itot) + tk*dlog( theta/32.0d0 )
c..   compromise value
         u(ia) = ualow*(1.0d0 - xal) + uahi*xal
c..   correct for very high temperature were n and p dominate
c..   low T value = all He4,n,p approximation
         uplow = u(ia)
         yp    = 2.0d0/theta * exp( (uplow*0.25d0 - qeq(in) )/tk  )
c..   add saturation at high T
c     yp    = yp/(1.0 + 2.0*yp)
         yp    = dmin1(yp, 0.5d0)
         xp    = 2.0d0*yp
c..   highest T value = all n, p approximation (Yp=yn=1/2)
         uphi = 4.0d0*(qeq(in) + tk*dlog( theta/4.0d0 ) )
c..   compromise value
         u(ia) = uplow*(1.0d0 - xp) + uphi*xp
c..   uses u(ia) = u + m here
         uaaa  = u(ia)
c..   save guess for convergence diagnostics
         uaaa0 = uaaa
c..   interate for nucleon number (sum Xi = 1).....................
         do n = 1, nloop
c..   generate trial yeq values
            do i = 1, itot  
               aa     =  float( nz(i) + nn(i) )
               arat   = aa**1.5d0 * w(i) / theta
               ueff   = aa*0.25d0 * uaaa - qeq(i)
               yeq(i) = arat * exp( ueff / tk )
               xeq(i) = yeq(i)*aa
            enddo
            xeqm  = -1.0d0
            dxdmu = 0
            do i = 1, itot  
               aa    =  float( nz(i) + nn(i) )
               xeqm  = xeqm + aa * yeq(i)
               dudmu = aa*0.25
               dydmu = yeq(i) / tk * dudmu
               dxdmu = dxdmu + aa*dydmu
            enddo 
            if( dxdmu .ne. 0.0d0 )then
               delua  = - xeqm / dxdmu
c..   limit size of change in mu(alpha)
               adelua  = abs( delua )
               if( adelua .gt. delmx )then
                  delua = delua/adelua* delmx
               endif
               uaaa = uaaa + delua
            else
               write(*,*)' error: dxdmu = ',dxdmu
               stop
            endif
            if( abs(xeqm) .lt. tiny )goto 999
         enddo
         write(*,'(a25,i5,2(a10,1pe12.3))')
     1        "nonconvergence (equil), n =",n," mu(a)=",uaaa,
     1        "xeq - 1",xeqm
         ierr  = 1
         niter = n
         return
 999     continue
c...................................end of iteration logic......
c..   NSE values are now determined
         niter = n
c..   uses mass excesses for energy calculation....
c..   set zero energy to cold, pure most bound nucleus in network
c     enc =  bpern * 9.65d+17
c..   relative to C12 nuclei
         enc = 0
         pnc = 0
         do i = 1, itot
            enc = enc + yeq(i)*( 1.2476d+17*t9 + qeq(i)*9.65d+17 )
            pnc = pnc + yeq(i)*  1.2476d+17*t9 *1.5d0 
         enddo
         pnc  = pnc * rho
         gam4 = pnc/(enc*rho) + 1.0d0
         yeq(ndim) = yeeq
         write(*,'(i5,1p13e11.3)')n,t9,rho,enc*1.0d-18,xeqm,
     1        xeq(ia),xeq(ni56),xeq(ip),uaaa,pnc*1.0d-18/rho,uaaa0
c         write(33,'(i5,1p10e11.3)')n,t9,rho,enc*1.0d-18,xeqm,
c     1        xeq(ia),xeq(ni56),xeq(ip),uaaa,pnc*1.0d-18/rho,uaaa0

c----------------------------
      endif
c----------------------------

      return
      end





