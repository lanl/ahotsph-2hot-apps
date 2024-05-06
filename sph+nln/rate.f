      subroutine rate(t90,rho,dtstar)

c..   fkt rate formula with T9, V derivatives (wda 5/9/99)
c..   weak rates added (pay 11/14/2004)
c..   weak screening (wda 11/14/04)
c..   statistical inverse rates (wda 11/14/04)

      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'
      include 'ceoset'
      include 'cburn'
      include 'cconst'

      integer*4 tdex, rhodex, rhodexw
      real*8 t9array(tsize), rhoarray(rhosize)

      real*8 dtstar, linklim
      real*8 t9, rho, t913, t953, t9ln, expon, exponu
      real*8 t92, t923, t943, factor, third
      real*8 yee, t90, t9lim, logrho, logrhoye

      integer*4 i1, i2, i3, i4, i5, i6

c..   irev  is the index of the reverse rate
c..   iline is the ordering number for subrates of a given reaction
c..   (this includes nonresonant and resonant rates)
c..   isum  is a dummy variable for counting
c      integer*4 isum
      real*8 fak

c..t9lim is the same value as tburnlo parameter in cdeuter
      parameter( third = 1.0d0/3.0d0, t9lim = 0.01d0 )

      integer*4 k, k1, k2, j, kzone

      integer*4 i

      data t9array/0.01d0,0.10d0,0.20d0,0.40d0,0.70d0,
     1     1.0d0,1.5d0,2.0d0,3.0d0,5.0d0,10.0d0,30.0d0,100.0d0/
      data rhoarray/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,
     1     7.0d0,8.0d0,9.0d0,10.0d0,11.0d0/
      data linklim/1.0d14/
c----------------------------------------------------
c..   avoid bad fkt extrapolations to very low temperature by setting
c..   a basement value 
      if( t90 .gt. t9lim )then
         t9 = t90
      else
         t9 = t9lim
      endif
      
c..   electron mole number Ye = Ne/rho*avagadro
      yee = 0.0d0
      do j = 1, nnuc
         yee = yee + y(j) * dble( nz(j) )
      enddo

      t913 = t9**third
      t953 = t9 * t913 *t913
      t9ln = dlog(t9)

      t92  = t9*t9
      t923 = t913*t913
      t943 = t913*t9

c..   calculate nuclear partition functions
      logrho = dlog10(rho)
      if(logrho .le. 1.0d0)logrho=1.0d0
      call locate(t9array,tsize,t9,tdex)
      call locate(rhoarray,rhosize,logrho,rhodex)
      logrho = dlog10(rho)
      call partfun(t9)
      
c..   define sig, sigt (sigv is calculated from sig below)
      do k = 1,ireac
         if( irev(k) .ge. 0 )then
c..   avoid reverse ffnu rates
            expon = rcoef(1,k) + rcoef(2,k)/t9 + rcoef(3,k)/t913
     1           + rcoef(4,k)*t913 + rcoef(5,k)*t9
     2           + rcoef(6,k)*t953 + rcoef(7,k)*t9ln
            if( expon .gt. -300.0d0 )then
               sig(k) = dexp(expon)
               factor = -rcoef(2,k)/t92 -third*rcoef(3,k)/t943
     1              +third*rcoef(4,k)/t923 + rcoef(5,k)
     2              +5.0d0*third*rcoef(6,k)*t923 + rcoef(7,k)/t9
               sigt(k) = sig(k) * factor
            else
               sig(k)  = 0.0d0
               sigt(k) = 0.0d0
            endif
         endif
        
      enddo

c..   initialize
      do k = 1, ireac
         signue(k) = 0.0d0
      enddo
      
c..   calculate density factors, deck by deck
c..   deck 1: electron capture density*Ye factor
      k1 = k1deck(1)
      k2 = k2deck(2)
      do k = k1, k2
c         if( rlkh(k) .eq. ' ffn' )then
         if( iffn(k) .ne. 0 )then
            logrhoye = dlog10(rho*x(ndim))
            call locate(rhoarray,rhosize,logrho,rhodexw)
            if(t9 .gt. t9lim .and. logrhoye .gt. 1.0d0)then

               call weakintp(t9,logrhoye,iffn(k),expon,exponu,
     1              tdex,rhodexw,k)

               if( expon .gt. -99.999d0 )then
                  sig(k)  = 10.0d0**(expon)
                  sigt(k) = 0.0d0
               else
                  sig(k)  = 0.0d0
                  sigt(k) = 0.0d0
               endif
               if( exponu .gt. -99.999d0 )then
                  signue(k) = 10.0d0**(exponu)
               else
                  signue(k) = 0.0d0
               endif
            else
               if( irev(k) .lt. 0 )then
c..   ffnu reverse rates set to zero at low temperature
                  sigt(k)   = 0.0d0
                  sig(k)    = 0.0d0
                  signue(k) = 0.0d0
               else
c..   only lab forward rates
c..   pass through
               endif
            endif
         else
c..   no fuller,fowler,newman
            signue(k) = 0.0d0
         endif
      enddo
      
      do k = k1, k2
         if( rlkh(k) .eq. '  ec' )then
            sig(k)  = yee * rho * sig(k)
            sigt(k) = yee * rho * sigt(k)
            sigv(k) = -sig(k) * rho
         endif
      enddo

c..   decks 2 and 3 have no density dependence
      k1 = k1deck(2)
      k2 = k2deck(3)
      do k = k1, k2
         sigv(k) = 0.0d0
      enddo
c..   decks 4 to 7 have binary interaction in entrance channel
      k1 = k1deck(4)
      k2 = k2deck(7)
      do k = k1, k2
         if( nrr(1,k) .eq. nrr(2,k) )then
c..   pairs reduced for identical particles to avoid double counting
c..   sigv is constructed from sig
            sig(k)  = 0.5d0 * sig(k)
            sigt(k) = 0.5d0 * sigt(k)
         endif

c..   this is another power of rho
         sigv(k)  = -rho * rho * sig(k)
         sig(k)   =  rho * sig(k)
         sigt(k)  =  rho * sigt(k)
      enddo
c..   deck 8 has triple (tertiary collision) in entrance channel
      k1 = k1deck(8)
      k2 = k2deck(8)
      do k = k1, k2
         if( nrr(1,k) .eq. nrr(2,k) .and. nrr(2,k) .eq. nrr(3,k))then
c..   waf, grc, baz (1967) ann rev a & a 5, 558 definition
c     has an extra factor of 1/2 for 3-alpha
c..   from 3 He4 destroyed per reaction, n/3! triples --> 1/2
c..   d/dV gives a compensating factor of 2
c..   factor of 1/3 moved from rhside and jacob to here 8/20/00
            sigv(k)  = -rho * rho * rho * sig(k)/3.0d0
            sig(k)   =        rho * rho * sig(k)  /6.0d0
            sigt(k)  =        rho * rho * sigt(k) /6.0d0
         else
            sigv(k) = -rho * rho * rho * sig(k) * 2.0d0
            sig(k)  =        rho * rho * sig(k)
            sigt(k) =        rho * rho * sigt(k)
         endif
      enddo
      

c.. deck 4: captures, i + j ---> k
      
      do k= k1deck(4),k2deck(4)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
c..   forward rates k
c..   irev = 0 for no reverse rate (eg. weak decays)
c..   irev < 0 for reverse rate k 
c..   (to be dealt with from forward rate elsewhere)
            if( irev(k) .gt. 0 )then
c..   k is forward rate, irev(k) is reverse rate
               fak =9.8678d9 * pf(i1)*pf(i2)/pf(i3) / rho
     1              *( anuc(i1)*anuc(i2)/anuc(i3) * t9 )**1.5d0
     2              *exp( -11.605d0*qval(k)/t9 )
c               fak =pf(i1)*pf(i2)/pf(i3) 
c     1              *( anuc(i1)*anuc(i2)/anuc(i3) )**1.5d0
c     2              *exp( -11.605d0*qval(k)/t9 )
c..   overwrite fkt value
               sig(irev(k)) = fak * sig(k)
               sigt(irev(k)) = fak * sigt(k)
               sigv(irev(k)) = fak * sigv(k)
               if( nrr(1,k) .eq. nrr(2,k) )then
                  sig(irev(k)) = 2.0d0*sig(irev(k))
                  sigt(irev(k)) = 2.0d0*sigt(irev(k))
                  sigv(irev(k)) = 2.0d0*sigv(irev(k))
               endif
            elseif( irev(k) .lt. 0 )then
               write(*,*)k,irev(k),' irev ERROR',
     1              ideck(k),ideck(-irev(k))
               write(*,'(2i5,3a6,5x,a5,a1,1p3e12.3,3i5,1p6e11.3)')
     1              k,ideck(k),
     1              (rname(i,k),i=1,3),rlkh(k),rvw(k),sig(k),
     2              sig(irev(k))*pf(i1)*pf(i2)/pf(i3) ,
     2              fak,k,irev(k),iline(k)
               stop'rate 1'
            else

            endif

      enddo

c deck 5: exchange, i + j ---> k + l
      do k= k1deck(5),k2deck(5)
         if ( irev(k) .gt. 0 )then
c..forward rate, overwrite fkt value with detailed balance result
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)

            fak = pf(i1)*pf(i2)/(pf(i3)*pf(i4))
     1           *( anuc(i1)*anuc(i2)/(anuc(i3)*anuc(i4)))**1.5d0
     2           *exp( -11.605d0*qval(k)/t9 ) 

c..   overwrite fkt value
            sig(irev(k)) = fak * sig(k)
            sigt(irev(k)) = fak * sigt(k)
            sigv(irev(k)) = fak * sigv(k)
            if( nrr(1,k) .eq. nrr(2,k) )then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 2.0d0*sigv(irev(k))
            endif

         endif
      enddo

c deck 6: exchange, i + j ---> k + l + m
c..use fkt for now ccccccccccccccccccccccccccccccccccccccccccccc
      do k= k1deck(6),k2deck(6)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            if( irev(k) .gt. 0 )then
c..add partition function ratio for inverse reaction a la friedel
               fak=(pf(i3)*pf(i4)*pf(i5)/pf(i1)/pf(i2))*rho/9.8678d9
     1          *(anuc(i1)*anuc(i2)/(anuc(i3)*anuc(i4)*anuc(i5)))**1.5d0
     2          *exp( -11.605d0*qval(k)/t9 ) / t9**1.5d0
               sig(irev(k)) = fak * sig(k)
               sigt(irev(k)) = fak * sigt(k)
               sigv(irev(k)) = fak * sigv(k)
            endif


      enddo

c deck 7: exchange, i + j ---> k + l + m + n
c..use fkt for now ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k= k1deck(7),k2deck(7)
c..no reverse rates in 176 nuc network
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            if( irev(k) .gt. 0 )then
              fak=(pf(i3)*pf(i4)*pf(i5)*pf(i6)*(rho/9.8678d9)**2.0d0
     1        /pf(i1)/pf(i2))*(anuc(i1)*anuc(i2)/
     2        (anuc(i3)*anuc(i4)*anuc(i5)*anuc(i6)))**1.5d0
     2        *exp( -11.605d0*qval(k)/t9 ) /t9**3.0d0
              sig(irev(k)) = fak * sig(k)
              sigt(irev(k)) = fak * sigt(k)
              sigv(irev(k)) = fak * sigv(k)
            endif
c            write(*,'(2i5,6a6,5x,a5,a1,1p3e12.3,2i5)')
c     1           k,ideck(k),(rname(i,k),i=1,6),rlkh(k),rvw(k),
c     2           sig(k),sig(irev(k)),qval(k),irev(k),iline(k)

      enddo

c deck 8:  3 ---> 1  and  3 ---> 2
c..use fkt for now, it looks ok for triple alpha
      do k= k1deck(8),k2deck(8)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            if( irev(k) .gt. 0 )then
               if (i5 .eq. 0) then
                 fak=(pf(i4)/pf(i1)/pf(i2)/pf(i3))*(9.8678d9/rho)**2.0d0
     1           *(anuc(i1)*anuc(i2)*anuc(i3)/(anuc(i4)))**1.5d0
     2           *exp( -11.605d0*qval(k)/t9 ) *t9**3.0d0
                 sig(irev(k)) = fak * sig(k)                 
                 sigt(irev(k)) = fak * sigt(k)
                 sigv(irev(k)) = fak * sigv(k)
               else
                fak=(pf(i4)*pf(i5)/pf(i1)/pf(i2)/pf(i3))*9.8678d9/rho
     1          *(anuc(i1)*anuc(i2)*anuc(i3)/(anuc(i4)*anuc(i5)))**1.5d0
     2          *exp( -11.605d0*qval(k)/t9 )*t9**1.5d0
                sig(irev(k)) = fak * sig(k)
                sigt(irev(k)) = fak * sigt(k)
                sigv(irev(k)) = fak * sigv(k)
               endif
               
           endif
           if( nrr(1,k) .eq. nrr(2,k) .and. 
     1          nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           endif
c..triple alpha for this deck

      enddo
      
c..limit rate links to mitigate roundoff
      do k= k1deck(4),k2deck(8)
         if( irev(k) .eq. 0 )then
c..no reverse rate so only forward rate needs scaling
            fak = sig(k) * dtstar
            if( fak .gt. linklim )then
               sig(k) = linklim/dtstar
            endif
c            write(*,'(2i6,7a6,1p8e12.3)')k,irev(k),
c     1           (rname(i,k),i=1,7),sig(k),sig(k)*dtstar,
c     2           fak,linklim/fak,linklim/dtstar
         elseif( irev(k) .gt. 0 )then
c..reverse rate has been defined, so it will need scaling too
            if( sig(k)*dtstar .gt. linklim .or. 
     1           sig(irev(k))*dtstar .gt. linklim )then
c               write(*,'(2i6,7a6,1p8e12.3)')k,irev(k),
c     1              (rname(i,k),i=1,7),sig(k),sig(irev(k)),
c     2              linklim/fak,linklim/dtstar,
c     3              sig(k)/( sig(k) + sig(irev(k)))*linklim/dtstar,
c     3              (1.0d0-sig(k)/( sig(k) + sig(irev(k))))
c     4              *linklim/dtstar
               fak = sig(k)/( sig(k) + sig(irev(k)))
c..   scale both rates, preserving ratio
               sig(k)       = linklim/dtstar* fak
               sig(irev(k)) = linklim/dtstar*(1.0d0 - fak)
            endif
         endif
      enddo
      

c............................................................
      
      return
      
      end


      subroutine screen(t9,rho)

c..   weak, intermediate and strong screening wda 12/22/04
c..   after Graboske, et al, ApJ 181, page 465 (1973), table 4
c..   see also DeWitt, Graboske, Cooper, 1973, Apj companion paper
c..   and John Bahcall's exportenergy.f subroutine

c..   kz is zone index

      implicit none

      include 'dimen'
      include 'crate'
      include 'cburn'
      include 'comcsolve'

      real*8 t9,rho,opbexp,z1,z2,z12,z3,amu,emu,xtr,zet,pfmc2,efmkt
      real*8 fprf,degd,t9m32
      real*8 weakscrn
      parameter( opbexp= 1.86d0 , weakscrn = 0.03d0 )
      real*8 xxl,xxl6,xxl8,zcurl,zbar,z58,z28,z33,tm1,uwk
      real*8 uint,ustr

      real*8 f0(nreac),zprd(nreac),z86(nreac)
      real*8 zprd3(nreac),z863(nreac)
      real*8 z53(nreac),z43(nreac),z23(nreac)
      real*8 z533(nreac),z433(nreac),z233(nreac)
      real*8 utot(nreac),dscr(nreac),dsct(nreac)

      integer*4 kz,i,k,calls
      integer*4 j

      save

      data calls/0/
c---------------------------------------------------------------

      if( calls .eq. 0 )then
c..only for first time through
         calls = calls + 1
c         write(*,*) 'first time in screen ',calls,kz
c..   set up constants
         do k = 1, nreac
            zprd(k) = 0.0d0
            z86(k)  = 0.0d0
            z53(k)  = 0.0d0
            z43(k)  = 0.0d0
            z23(k)  = 0.0d0
         enddo
c..deck 4 a + b ---> c
c..deck 5 a + b ---> c + d
c..deck 6 a + b ---> c + d + e
c..deck 7 a + b ---> c + d + e + f
         do k = k1deck(4),k2deck(7)
            z1      = nz(nrr(1,k))
            z2      = nz(nrr(2,k)) 
            zprd(k) = z1*z2 
            z86(k)  = (z1+z2)**opbexp - z1**opbexp - z2**opbexp
            z53(k)  = (z1+z2)**(5.0d0/3.0d0) - z1**(5.0d0/3.0d0) 
     1           - z2**(5.0d0/3.0d0)
            z43(k)  = (z1+z2)**(4.0d0/3.0d0) - z1**(4.0d0/3.0d0) 
     1           - z2**(4.0d0/3.0d0)
            z23(k)  = (z1+z2)**(2.0d0/3.0d0) - z1**(2.0d0/3.0d0) 
     1           - z2**(2.0d0/3.0d0)
         enddo
c..deck 8 a + b + c = d or d + e
         do k = k1deck(8),k2deck(8)
            z1      = nz(nrr(1,k))
            z2      = nz(nrr(2,k))
            z3      = nz(nrr(3,k))
            z12     = nz(nrr(1,k))+nz(nrr(2,k))
            zprd(k) = z1*z2 
            zprd3(k)= z12*z3
            z86(k)  = (z1+z2)**opbexp - z1**opbexp - z2**opbexp
            z863(k) = (z12+z3)**opbexp-z12**opbexp - z3**opbexp
            z53(k)  = (z1+z2)**(5.0d0/3.0d0) - z1**(5.0d0/3.0d0) 
     1           - z2**(5.0d0/3.0d0)
            z533(k) = (z12+z3)**(5.0d0/3.0d0) - z12**(5.0d0/3.0d0) 
     1           - z3**(5.0d0/3.0d0)
            z43(k)  = (z1+z2)**(4.0d0/3.0d0) - z1**(4.0d0/3.0d0) 
     1           - z2**(4.0d0/3.0d0)
            z433(k) = (z12+z3)**(4.0d0/3.0d0) - z12**(4.0d0/3.0d0) 
     1           - z3**(4.0d0/3.0d0)
            z23(k)  = (z1+z2)**(2.0d0/3.0d0) - z1**(2.0d0/3.0d0) 
     1           - z2**(2.0d0/3.0d0)
            z233(k)  = (z12+z3)**(2.0d0/3.0d0) - z12**(2.0d0/3.0d0) 
     1           - z3**(2.0d0/3.0d0)
         enddo
      endif

c..deck 8 a + b + c ---> d or
c..       a + b + c ---> d + e
c..treat 3 body screening like 3 body rate after Salpeter & Grebokse
c..if rate(tot) = rate(rxn1)*rate(rxn2) multiply screening f1*f2
c..now that constants are defined, calculate variables
      do i=1, nreac
         f0(i)   = 1.0d0
         utot(i) = 0.0d0
      enddo
c..   amu=Ytot,emu=Ye,xtr is an interpolation variable,zet=ave Z**2
      amu = 0.0d0
      emu = 0.0d0
      xtr = 0.0d0
      zet = 0.0d0
      do i = 1, inuc
         amu = amu + y(i)
         emu = emu + y(i)*dble( nz(i) )
         xtr = xtr + y(i)*dble( nz(i) )**1.58d0
         zet = zet + y(i)*dble( nz(i) )**2
      enddo

      pfmc2 = 1.017677d-4*( rho*emu )**0.666666667d0
      efmkt = 5.92986d0/t9*( sqrt( 1.0d0 + pfmc2 ) - 1.0d0 )
      if( efmkt .le. 1.0d-2 )then
         fprf = 1.0d0
      else
         degd = log10( efmkt )
         if( degd .ge. 1.5d0 )then
            fprf = 0.0d0
         else
            fprf =0.75793-0.54621*DEGD-0.30964*DEGD**2+0.12535*DEGD**3+
     1      0.1203*DEGD**4-0.012857*DEGD**5-0.014768*DEGD**6
         endif
      endif

      t9m32 = 1.0d0/t9**1.5d0
      XXL   = 5.9426d-6*T9M32*DSQRT(RHO*AMU)
      XXL6  = XXL**0.666667d0
      XXL8  = XXL**0.86d0
      ZCURL = DSQRT((ZET+FPRF*EMU)/AMU)
      ZBAR  = EMU/AMU
      Z58   = ZCURL**0.58d0
      Z28   = ZBAR**0.28d0
      Z33   = ZBAR**0.333333d0
      TM1   = XXL*ZCURL

c..for each reaction
c..deck 4 a + b ---> c
c..deck 5 a + b ---> c + d
      do i = k1deck(4),k2deck(7)
         uwk  = tm1*zprd(i)
         uint = 0
         ustr = 0
         if( uwk .le. weakscrn )then
            utot(i) = uwk
            dscr(i) =  0.5d0*uwk
            dsct(i) = -1.5d0*uwk
         else
            uint = 0.38d0*xxl8*xtr*z86(i)/(amu*z58*z28)
            if( uint .le. 2.0d0 )then
               utot(i) = uint
               dscr(i) =  0.43d0*uint
               dsct(i) = -1.29d0*uint
             else
                ustr = 0.624d0*z33*xxl6*(z53(i)+0.316d0*z33*z43(i)
     1               +0.737d0*z23(i)/(zbar*xxl6))
                if( ustr .lt. uint .or. uwk .ge. 5.0d0 )then
                   utot(i) = ustr
                   dscr(i) = 0.208d0*z33*(z53(i) + 0.316d0*z33*z43(i))
     1                  *xxl6
                   dsct(i) = -3.0d0*dscr(i)
                else
                   utot(i) = uint
                   dscr(i) =  0.43d0*uint
                   dsct(i) = -1.29d0*uint
                endif
             endif
         endif
      enddo

      do i = k1deck(8),k2deck(8)
         uwk  = tm1*zprd(i)+tm1*zprd3(i)
         uint = 0
         ustr = 0
         if( uwk .le. weakscrn )then
            utot(i) = uwk
            dscr(i) =  0.5d0*uwk
            dsct(i) = -1.5d0*uwk
         else
            uint = 0.38d0*xxl8*xtr*z86(i)/(amu*z58*z28)
     1             +0.38d0*xxl8*xtr*z863(i)/(amu*z58*z28)
            if( uint .le. 2.0d0 )then
               utot(i) = uint
               dscr(i) =  0.43d0*uint
               dsct(i) = -1.29d0*uint
             else
                ustr = 0.624d0*z33*xxl6*(z53(i)+0.316d0*z33*z43(i)
     1               +0.737d0*z23(i)/(zbar*xxl6))
     2               +0.624d0*z33*xxl6*(z533(i)+0.316d0*z33*z433(i)
     1               +0.737d0*z233(i)/(zbar*xxl6))    
                if( ustr .lt. uint .or. uwk .ge. 5.0d0 )then
                   utot(i) = ustr
                   dscr(i) = 0.208d0*z33*(z53(i) + 0.316d0*z33*z43(i))
     1                  *xxl6
     2                  +0.208d0*z33*(z533(i) + 0.316d0*z33*z433(i))
     1                  *xxl6
                   dsct(i) = -3.0d0*dscr(i)
                else
                   utot(i) = uint
                   dscr(i) =  0.43d0*uint
                   dsct(i) = -1.29d0*uint
                endif
             endif
         endif
      enddo

      do i = k1deck(4),k2deck(8)         
            f0(i) = exp (utot(i))

      enddo


c..   adjust for screening
c..   can add T and V derivatives from f0 here
      do i = 1,nreac
           sig(i)  = sig(i)  * f0(i)
           sigt(i) = sigt(i) * f0(i)
           sigv(i) = sigv(i) * f0(i)
      enddo



      return
      end



