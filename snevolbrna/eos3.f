      subroutine eossetup
c
      common /typef/ ieos
      ieos=4
c
      call unit
c
c--get nuclear data
c
      call nucdata
c
c--initialize the swesty eos
c--note that this assumes that the density are known for the dump
c
      call loadmx
c
      return
      end
c
      subroutine eos3(dtburn,massi,rhoi,ui,dui,u2i,
     $     yei,tempi,ifleosi,abari,xpi,xni,
     $     xpfi, p2i, p3i, p4i,temprev,rhoprev,
     $     xpprev,xnprev,yeprev,ufreezi,xhe4ni,xc12ni,
     $     xo16ni,xne20ni,xmg24ni,xsi28ni,xni56ni,
     $     xco56ni,xfe56ni, 
     $     dtbrni)
c*************************************************************
c
c     compute pressures and temperatures with the    
c     Ocean eos assuming NSE
c
c************************************************************
c
      double precision tgo
      double precision umass, uinput
      double precision totdecay, totco, totfe
      double precision rhoi, ui, u2i, tempi, yei, ptot, pri, cs, etai, 
     $     abari, xpi, xni, xai, xhi, yehi, rhoold, yeold,
     $     xpfi, p2i, p3i, p4i, ufreezi, xmuhi, stot,vsoundi,
     $     temprev,rhoprev,xpprev,xnprev,yeprev,xmuei,xmuhati,
     $     xalphai,xheavyi,pmassi
      real dtburn, massi
      double precision dtbrni, uburni,dui,udecayi,dt,
     $     xhe4ni,xc12ni,xo16ni,xne20ni,xmg24ni,xsi28ni,xni56ni,
     $     xco56ni,xfe56ni
c
      common /output/ vsoundi,pri,etai,yehi,xmuei,xmuhati,xalphai,
     $     xheavyi
c
      common /konst/ gg, clight, arad, bigr, xsecnn, xsecne
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /typef/ ieos
      common /tmpc/ tgo 
c     
      pmassi=dble(massi)
      if (xhe4ni.lt.1.e-25) then
         xhe4ni=0.
      end if
      if (xc12ni.lt.1.e-25) then
         xc12ni=0.
      end if
      if (xo16ni.lt.1.e-25) then
         xo16ni=0.
      end if
      if (xne20ni.lt.1.e-25) then
         xne20ni=0.
      end if
      if (xmg24ni.lt.1.e-25) then
         xmg24ni=0.
      end if
      if (xsi28ni.lt.1.e-25) then
         xsi28ni=0.
      end if
      if (xni56ni.lt.1.e-25) then
         xni56ni=0.
      end if
      if (xco56ni.lt.1.e-25) then
         xco56ni=0.
      end if
      if (xfe56ni.lt.1.e-25) then
         xfe56ni=0.
      end if
      if (yei.lt.0.) then
         print *,'yei',yei
         stop
      endif
      tgo=tgo+1.

c     Commented out 8/18/2004 (gmr, clf)
c      abari=1d0/(xhe4ni/4.d0+xc12ni/12.d0+
c     $     xo16ni/16.d0+
c     $     xne20ni/20.d0+xmg24ni/24.d0+
c     $     xsi28ni/28.d0+xni56ni/56.d0+
c     $     xco56ni/56.d0+xfe56ni/56.d0+
c     $     xpi+xni)
c
c--assume chemical freeze-out
c
c      goto 20
      ui=max(1.d-5,ui)
c      print *, 'calling eosfl',rhoi,ui
      call eosfl(rhoi,pri,
     $     ui,u2i,yei,tempi,ifleosi,abari,
     $     xpi,xni,xpfi,p2i,p3i,p4i,ufreezi,
     $     xmuei,xmuhi,etai,temprev,
     $     yeprev,xpprev,xnprev)

c     Commented out 8/18/2004 (gmr, clf)
c      abari=1d0/(xhe4ni/4.d0+xc12ni/12.d0+
c     $     xo16ni/16.d0+
c     $     xne20ni/20.d0+xmg24ni/24.d0+
c     $     xsi28ni/28.d0+xni56ni/56.d0+
c     $     xco56ni/56.d0+xfe56ni/56.d0+
c     $     xpi+xni)
      if (ifleosi.eq.0) then
         call rootemp2u(rhoi,ui,tempi,yei,abari,
     1        ptot,cs,etai,stot)
c         xpi=0.
         xni=0.
         xmuei=etai*tempi
         xmuhati=0.0
         xalphai=0.0
         xheavyi=1.
         yehi=yei
c
c--ocean eos + burning
c
      elseif (ifleosi.eq.1) then
         call rootemp2u(rhoi,ui,tempi,yei,abari,
     $        ptot,cs,etai,stot)
         xmuei=etai*tempi
         xmuhati=0.0
         xalphai=0.0
         xheavyi=1.
         yehi=yei
c     
c--Swesty and Lattimer eos
c
      else
         if (ieos.eq.4) then
            uinput=u2i
         end if
         call rootemp4(rhoi,uinput,yei,tempi,xpfi,p2i,p3i,p4i,
     $        ptot,cs,etai,xpi,xni,xai,xhi,yehi,abari,xmuhi,stot)
c-- assume all dissociated, this could be changed in the future
         xalphai=xai
         xheavyi=xhi
         xmuei=etai*tempi
         xmuhati=xmuhi
      endif
c
c--store values
c
      if (xpi.le.1d-10) xpi=0.
      if (xni.le.1d-10) xni=0.
      if (xhe4ni.le.1d-10) xhe4ni=0.
      if (xc12ni.le.1d-10) xc12ni=0.
      if (xo16ni.le.1d-10) xo16ni=0.
      if (xne20ni.le.1d-10) xne20ni=0.
      if (xmg24ni.le.1d-10) xmg24ni=0.
      if (xsi28ni.le.1d-10) xsi28ni=0.
      if (xni56ni.le.1d-10) xni56ni=0.
      xtot=xpi+xni+xhe4ni+xc12ni+xo16ni+xne20ni
     $     +xmg24ni+xsi28ni+xni56ni
      xhe4ni=xhe4ni/xtot
      xc12ni=xc12ni/xtot
      xo16ni=xo16ni/xtot
      xne20ni=xne20ni/xtot
      xmg24ni=xmg24ni/xtot
      xsi28ni=xsi28ni/xtot
      xni56ni=xni56ni/xtot
      vsoundi=min(cs,0.3333*dble(clight))
      pri=ptot
      temprev=tempi
      rhoprev=rhoi
      xnprev=xni
      xpprev=xpi
      yeprev=yei
c
c-- call burning network
c
c      if (mod(tgo,195000).eq.43197) then
c         print *,rhoi,ui
c         ui=135.678787
c         print *, xhe4ni,xc12ni,xo16ni,xne20ni,xmg24ni,xsi28ni,
c     $        xni56ni,xpi,xni
c      end if
c      if (mod(tgo,195000).gt.40000.and.
c     $     mod(tgo,195000).lt.60000) then
c         print *, mod(tgo,195000), tempi, rhoi, ui
c      end if
c      if (tempi.gt.13.0) read *, dgoingon
c      if (tgo.gt.1953) read *, dgoingon
c      if (tgo.gt. 120000) then
c         print *, tgo,dtburn,dtbrni,tempi,rhoi
c         print *, xhe4ni, xc12ni,xo16ni,xne20ni,xmg24ni,
c     $        xsi28ni,xni56ni,xpi,xni
c      end if

c     Turned off burning - 8/18/2004 (gmr, clf)
c      call burner(dtburn,dtbrni,tempi,rhoi,ui,
c     $     uburni,dui,abari,pmassi,etai,ifleosi,
c     $     xhe4ni,xc12ni,xo16ni,xne20ni,xmg24ni,
c     $     xsi28ni,xni56ni,xpi,xni,iexpl)
c
      dt=dble(dtburn)
c      call nidecay(dt,pmassi,xni56ni,xco56ni,xfe56ni,udecayi)

      udecayi=0.
      xco56ni=0.
      xfe56ni=0.
      totdecay=totdecay+pmassi*udecayi
      totco=totco+pmassi*xco56ni
      totfe=totfe+pmassi*xfe56ni
c
c      if (xni56ni.gt.0.5) then
c         print *, dt, ui, uburni, udecayi, xni56ni, xco56ni
c      end if
c      if ((abs(0.05-rhoi)).lt.0.0002.and.(abs(2.5-tempi)).lt.0.002) then
c         write(68,102) 2.d6*rhoi,1.d9*tempi, yei, dt, uburni, ui,
c     $        udecay, dtbrni,xhe4ni,xc12ni,xo16ni,xne20ni,
c     $        xmg24ni,
c     $        xsi28ni,xni56ni,xpi,xni
c         print *, xhe4ni,xc12ni,xo16ni,xne20ni,xmg24ni,
c     $        xsi28ni,xni56ni,xpi,xni,dtbrni
c      end if
 102  format(18(1pe12.4))

c     uburni is only modified in burner - commented out 8/18/2004 (gmr, clf)
c      ui=ui+uburni+udecayi
c      if (uburni.gt.0) print *, ui, uburni
      if (xhe4ni.lt.1.d-25) then
         xhe4ni=0.
      end if
      if (xc12ni.lt.1.d-25) then
         xc12ni=0.
      end if
      if (xo16ni.lt.1.d-25) then
         xo16ni=0.
      end if
      if (xne20ni.lt.1.d-25) then
         xne20ni=0.
      end if
      if (xmg24ni.lt.1.d-25) then
         xmg24ni=0.
      end if
      if (xsi28ni.lt.1.d-25) then
         xsi28ni=0.
      end if
      if (xni56ni.lt.1.e-25) then
         xni56ni=0.
      end if
      if (xco56ni.lt.1.e-25) then
         xco56ni=0.
      end if
      if (xfe56ni.lt.1.e-25) then
         xfe56ni=0.
      end if
c      if (mod(tgo,195000).eq.43197) then
c         rhoi=rhoiold
c      end if
c      if (mod(tgo,195000).eq.43198) then
c         rhoi=rhoiold
c      end if
c
c      print *, 'out of eos3',dt
 20   continue
      return
      end
c
      subroutine rootemp2s(rhoi,ui,tempi,yei,
     1                    abar,ptot,cs,eta,si)
c*****************************************************************
c                                                                *
c  Given rho, u, ye and an initial T, this                       *
c  subroutine iterates over temperatures with a Newton-Raphson   *
c  scheme coupled with bissection to determine thermodynamical   *
c  variables.                                                    *
c  The ocean eos is called to determine perfect gas, radiation   *
c  and electron/positron contributions, assuming complete        *
c  ionization. The coulomb subroutine determines Coulomb         *
c  corrections                                                   *
c                                                                *
c*****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (itmax=80)
      parameter (dtol=1d-2)
c
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/iarg/lst,kentr,kpar,jurs,jkk
c
c--change to the appropriate units
c
      t9=tempi*uotemp1
      rho=rhoi*uorho1
      s=si*uou1/uotemp1
c
c--set brackets for T iterations
c
      t9l=3d-9
      t9h=3d1
c
c--compute zbar from abar and ye
c
      zbar=yei*abar
c
c--compute coulomb correction (since coulomb corr. not dependant 
c  on T, and freeze-out is assumed call only once)
c
      call coulomb(rhoi,zbar,yei,ucoul,pcoul)
c
c--use Newton-Raphson to find T
c
      call nados(t9,rho,zbar,abar,pel,eel,sel,
     1           ptot,etot,stot,dpt,det,dst,dpd,ded,gamm,eta)
      sres=stot-s
c
      dt9=dabs(t9h-t9l)
      do k=1,itmax
         dt9old=dt9
         if (((t9-t9h)*dst-sres)*
     1      ((t9-t9l)*dst-sres).ge.0.d0.or.
     2      dabs(2.d0*sres).gt.dabs(dt9old*dst)) then
            dt9=0.5d0*(t9h-t9l)
            t9=t9l+dt9
         else
            dt9=sres/dst
            t9=t9-dt9
         endif
         if (dabs(dt9/t9).lt.dtol) goto 20
         call nados(t9,rho,zbar,abar,pel,eel,sel,
     1              ptot,etot,stot,dpt,det,dst,dpd,ded,gamm,eta)

         sres=stot-s
         if (sres.lt.0.d0) then
            t9l=t9
         else
            t9h=t9
         endif
      enddo
c
c--did not converge, print out error message and stop
c
      print *,'rootemp2: no convergence'
      stop
c
c--iteration sucessful, transform back in code units 
c
   20 continue
      tempi=t9*uotemp
      ptot=ptot*uopr+pcoul
      cs=dsqrt(gamm*ptot/rhoi)
      ui=etot/uou1+ucoul
c
      return
      end
c
      subroutine rootemp2u(rhoi,ui,tempi,yei,
     1                    abar,ptot,cs,eta,stot)
c*****************************************************************
c                                                                *
c  Given rho, u, ye and an initial T, this                       *
c  subroutine iterates over temperatures with a Newton-Raphson   *
c  scheme coupled with bissection to determine thermodynamical   *
c  variables.                                                    *
c  The ocean eos is called to determine perfect gas, radiation   *
c  and electron/positron contributions, assuming complete        *
c  ionization. The coulomb subroutine determines Coulomb         *
c  corrections                                                   *
c                                                                *
c*****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (itmax=80)
      parameter (dtol=1d-2)
c
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/iarg/lst,kentr,kpar,jurs,jkk
c
c--change to the appropriate units
c
      t9=tempi*uotemp1
      rho=rhoi*uorho1
      u=ui*uou1
c
c--set brackets for T iterations
c
      t9l=3d-9
      t9h=3d1
c
c--compute zbar from abar and ye
c
      zbar=yei*abar
c
c--compute coulomb correction (since coulomb corr. not dependant 
c  on T, and freeze-out is assumed call only once)
c
      call coulomb(rhoi,zbar,yei,ucoul,pcoul)
      ucoul=ucoul*uou1
c
c--use Newton-Raphson to find T
c
      call nados(t9,rho,zbar,abar,pel,eel,sel,
     1           ptot,etot,stot,dpt,det,dst,dpd,ded,gamm,eta)
      ures=ucoul+etot-u
      dut=det
c
      dt9=dabs(t9h-t9l)
      do k=1,itmax
         dt9old=dt9
         if (((t9-t9h)*dut-ures)*
     1      ((t9-t9l)*dut-ures).ge.0.d0.or.
     2      dabs(2.d0*ures).gt.dabs(dt9old*dut)) then
            dt9=0.5d0*(t9h-t9l)
            t9=t9l+dt9
         else
            dt9=ures/dut
            t9=t9-dt9
         endif
         if (dabs(dt9/t9).lt.dtol) goto 20
         call nados(t9,rho,zbar,abar,pel,eel,sel,
     1              ptot,etot,stot,dpt,det,dst,dpd,ded,gamm,eta)
         ures=etot+ucoul-u
         dut=det
         if (ures.lt.0.d0) then
            t9l=t9
         else
            t9h=t9
         endif
      enddo
c
c--did not converge, print out error message and stop
c
      print *,'rootemp2: no convergence'
      stop
c
c--iteration sucessful, transform back in code units 
c
   20 continue
      tempi=t9*uotemp
      ptot=ptot*uopr+pcoul
      cs=dsqrt(gamm*ptot/rhoi)
      stot=stot*uotemp1/uou1
c
      return
      end
c
      subroutine rootemp3(rhoi,ui,tempi,yei,rhoold,yeold,
     1                    ptot,cs,eta,yp,yn,xa,xh,yeh,abar,stot)
c****************************************************************
c
c  Given rho, u, an old T, rho, yp and yn this                   *
c  subroutine iterates over temperatures with a Newton-Raphson   *
c  scheme coupled with bissection to determine thermodynamical   *
c  variables.                                                    *
c  The ocean eos is called to determine perfect gas, radiation   *
c  and electron/positron contributions, assuming complete        *
c  ionization. The coulomb subroutine determines Coulomb         *
c  corrections. The nserho and nsetemp subroutines compute
c  abar, zbar, and the amount of energy gone into dissociating
c  nuclei.                                                       *
c
c*****************************************************************
c
      implicit double precision (a-h,o-z)
      parameter (itmax=80)
      parameter (dtol=1d-2)
c
      real  udist, udens, utime, uergg, uergcc
c
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/iarg/lst,kentr,kpar,jurs,jkk
c
c--change to the appropriate units
c
      t9=tempi*uotemp1
      t9old=t9
      rho=rhoi*uorho1
      u=ui*uou1
      t9l=1d-1
      t9h=2d2
c
c--transform in cgs for the NSE table
c
      rhocgs=rhoi*udens
      ucgs=ui*uergg
      call nserho(rhoold,yeold,rhocgs,yei,t9,yp,yn,
     1            xa,xh,yeh,zbar,abar,ubind,dubind)
      xfn=yp+yn
c
c--transform back in Ocean's units
c
      ediss=ubind/uergg*uou1
      dediss=dubind/uergg*uou1
c
c--add coulomb correction
c
      call coulomb(rhoi,zbar,yei,ucoul,pcoul)
      ucoul=ucoul*uou1
c
c--use Newton-Raphson to iterate for T
c
c--note abar is averaged over massive nuclei, zbar overall
      abar2=zbar/yei
      call nados(t9,rho,zbar,abar2,pel,eel,sel,
     1           ptot,etot,stot,dpt,det,dst,dpd,ded,gamm,eta)
      ures=(etot+ediss+ucoul)-u
      dut=det+dediss
      dt9old=dabs(t9h-t9l)
      dt9=dt9old
c
      do iter=1,itmax
         if (((t9-t9h)*dut-ures)*
     1      ((t9-t9l)*dut-ures).ge.0.d0.or.
     2      dabs(2.d0*ures).gt.dabs(dt9old*dut)) then
            dt9old=dt9
            t9old=t9
            dt9=0.5d0*(t9h-t9l)
            t9=t9l+dt9
         else
            dt9old=dt9
            t9old=t9
            dt9=ures/dut
            t9=t9-dt9
c-- this is because coming down from high temperatures,
c-- dediss can be too small
            if (t9.lt.0.5d0*t9old) t9=0.5d0*t9old
         endif
         if (t9.lt.0.5d0) then
            print *,'t9,t9old',t9,t9old
            print *,'u,ures,dut',u,ures,dut
            print *,'rhocgs,rhoold',rhocgs,rhoold
            print *,'yei,yeold',yei,yeold
            print *,'ediss,etot,ucoul',ediss,etot,ucoul
            stop
         endif
c
c--solve nse 
c
         call nsetemp(t9old,rhocgs,yei,t9,yp,yn,
     1                xa,xh,yeh,zbar,abar,ubind,dubind)
c
c--if tolerance is satisfied, exit
c
         if (dabs(dt9/t9).lt.dtol) goto 20
c
c--get derivatives for new iteration
c
         call coulomb(rhoi,zbar,yei,ucoul,pcoul)
         ucoul=ucoul*uou1
         abar2=zbar/yei
         call nados(t9,rho,zbar,abar2,pel,eel,sel,
     1              ptot,etot,stot,dpt,det,dst,dpd,ded,gamm,eta)
c
         ediss=ubind/uergg*uou1
         dediss=dubind/uergg*uou1
         ures=(etot+ediss+ucoul)-u
         dut=det+dediss
         if (ures.lt.0.d0) then
            t9l=t9
         else
            t9h=t9
         endif
      enddo
c
c--did not converge, print out error message and stop
c
      print *,'rootemp3: no convergence for part. i,rho(cgs)',
     1         rhoi
      stop
c
c--iteration sucessful, transform back in code units 
c
   20 continue
      tempi=t9*uotemp
      ptot=ptot*uopr+pcoul
      cs=dsqrt(gamm*ptot/rhoi)
      stot=stot*uotemp1/uou1
c
      return
      end
c
      subroutine rootemp4(rhoi,ui,yei,tempi,xpfi,p2i,p3i,p4i,
     1                   press,cs,eta,yp,yn,xa,xh,yeh,abar,xmuhi,stot)
c**************************************************************
c
c     This subroutine computes the temperature
c     with Doug Swesty's eos using
c     a Newton-Raphson procedure coupled with
c     bissection to prevent convergence problems.
c
c******************************************************
      implicit double precision (a-h,o-z)
c
      parameter(itmax=80)
      parameter(dtol=1d-2)
c
      real utemp, utmev, ufoe, umevnuc, umeverg
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
c
      double precision inpvar(4)
c
      templ=.1d0
      temph=1d3
      dtemp=dabs(temph-templ)
c-- all the u's in this routine will be MeV/baryons
c-- all the u's in this routine will be kB/baryons
      u=ui*uslu
      temp=tempi*usltemp
      inpvar(1)=temp
      inpvar(2)=p2i
      inpvar(3)=p3i
      inpvar(4)=p4i
      pprev=xpfi
      brydns=rhoi*uslrho
      call slwrap(inpvar,yei,brydns,pprev,
     1      psl,usl,dusl,gamsl,eta,yp,yn,xa,xh,yeh,abar,xmuh,u2sl)
      ures=usl-u
      do 10 k=1,itmax
         dtempold=dtemp
         if (((temp-temph)*dusl-ures)*
     1      ((temp-templ)*dusl-ures).ge.0.d0.or.
     2      dabs(2.d0*ures).gt.dabs(dtempold*dusl)) then
            dtemp=0.5d0*(temph-templ)
            temp=templ+dtemp
         else
            dtemp=ures/dusl
            temp=temp-dtemp
         endif
         if (dabs(dtemp/temp).lt.dtol) goto 20
         inpvar(1)=temp
         call slwrap(inpvar,yei,brydns,pprev,
     1            psl,usl,dusl,gamsl,eta,yp,yn,xa,xh,yeh,abar,xmuh,u2sl)
         ures=usl-u
         if (ures.lt.0.0d0) then
            templ=temp
         else
            temph=temp
         endif
   10 continue
      print *,'rootemp4: no convergence for particle i, rho',
     1         rhoi
   20 continue
      p2i=inpvar(2)
      p3i=inpvar(3)
      p4i=inpvar(4)
      xpfi=pprev
c-- convert back to code units
      press=psl/uslp
      tempi=temp/usltemp
      if (gamsl.gt.0.0) then
         cs=dsqrt(gamsl*press/rhoi)
      else
         write (*,*) 'gamsl is ', gamsl
         cs = 1e30
      endif
c--xmuhat is eta*T with T in code units
      xmuhi=xmuh/utmev
      stot=u2sl/u2slu
      return
      end
c
      subroutine coulomb(rhoi,zbar,ye,ucoul,pcoul)
c***********************************************************
c
c  compute Coulomb corrections as given in Shapiro 
c  and Teukolsky. p. 31 (2.4.9) and (2.4.11)
c  in cgs:
c        ucoul=-1.45079*e**2*avo**4/3*ye**4/3*rho**1/3*Z**2/3
c             =-1.70e13 Ye**4/3 * rho**1/3 * Z**2/3
c  code units: mulitply by udens**1/3 / uergg 
c
c        pcoul=-0.4836*e**2*avo**4/3*Ye**4/3*rho**4/3*Z**2/3
c             =-5.67e12 Ye**4/3 * rho**4/3 * Z**2/3
c  code units: mulitply by udens**4/3 / uergcc 
c
c***********************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter(ufac=-0.214d0)
      parameter(pfac=-0.0714d0)
c    
      rho13=rhoi**0.333333333333d0
      rho43=rho13*rhoi
      ye2=ye*ye
      y43z23=(ye2*zbar)**0.66666666666d0 
      ucoul=ufac*rho13*y43z23
      pcoul=pfac*rho43*y43z23
c
      return
      end
c
      subroutine nserho(rhoold,yeold,rho,ye,t9,yp,yn,
     1                  xa,xh,yeh,zbar,abar,ubind,dubind)
c*************************************************************
c
c this subroutine figures out the NSE eq. assuming that yp
c and yn were previously known at different density and ye,
c but SAME temperature
c
c**************************************************************
c
      implicit double precision(a-h,o-z)
      parameter (tolnse=1d-5,kmax=10)
c
c
      common /testnse/ testk,testzy,testay,testyp,testyn
c
c--finding zero point of binding energy
c
      ider=2
      call nsesolv(ider,t9,rhoold,yeold,yp,yn,kit,kmax,ubind0,
     &             xa,xh,yeh,zbar,abar)
c
      If(kit.ge.kmax) Then
          write(*,*) 'NSE mis-stored entering nserho'
          write(*,*) 'T9, rho, ye',t9,rhoold,yeold
          write(*,*) 'inconsistent with yp, yn',yp,yn
      Endif
      ypold=yp
      ynold=yn
      delye=ye-yeold
      rhovar=rho-rhoold
      ider=0
c
c--If delye small, skip ye variation.
c
      If(delye.eq.0.) GOTO 50
      yelast=yeold
      Do 40 i=1,100
          delye=dsign(min(abs(ye-yelast),abs(delye)),delye)
          yetmp=yelast+delye
          call nsesolv(ider,t9,rhoold,yetmp,yp,yn,kit,kmax,ubind,
     &                 xa,xh,yeh,zbar,abar)
          If (dabs(yetmp-ye).le.tolnse.and.kit.lt.kmax) goto 50
          If (kit.ge.kmax) then
             delye=0.5d0*delye
             yp=ypold
             yn=ynold
          Elseif(kit.lt.4) Then
             yelast=yetmp
             delye=2.d0*delye
             ypold=yp
             ynold=yn
          Else
             yelast=yetmp
             ypold=yp
             ynold=yn
          Endif
   40 Continue
      write(*,*)'Ye loop failure'
      write(*,*)t9,rhoold,rho
      write(*,*)yetmp,yelast,ye
      write(*,*)yeold,delye,yp,yn
      write(*,*)kit,zbar,abar
      write(*,*)kmax,ubind
      write(*,*)testk,testzy,testay,testyp,testyn
   50 Continue
c
c--Begin rho iteration
c
      ypold=yp
      ynold=yn
      rholast=rhoold
      If(dabs(rhovar).gt.1d7) Then
          delrho=dsign(max(1d7,0.125d0*rhovar),rhovar)
      Else
          delrho=rhovar
      Endif
c
      Do 60 i=1,500
          delrho=dsign(min(abs(rho-rholast),abs(delrho)),delrho)
          rhotmp=rholast + delrho
          call nsesolv(ider,t9,rhotmp,ye,yp,yn,kit,kmax,ubind,
     &                 xa,xh,yeh,zbar,abar)
c         write(*,90)'nserho: rhotmp,yp,yn,kit',rhotmp,yp,yn,kit
          If(dabs((rhotmp-rho)/rho).lt.tolnse.and.kit.lt.kmax) goto 70
          If (kit.ge.kmax) Then
              delrho=.5d0*delrho
              yp=ypold
              yn=ynold
          Elseif(kit.lt.4) Then
             rholast=rhotmp
             delrho=2.d0*delrho
             ypold=yp
             ynold=yn
          Else
              rholast=rhotmp
              ypold=yp
              ynold=yn
          Endif
   60 Continue
      write(*,*)'rho loop failure'
      write(*,*)t9,rhoold,rho,rhotmp
      write(*,*)yetmp,yelast,ye
      write(*,*)yeold,delye,yp,yn
      write(*,*)kit,zbar,abar
      write(*,*)kmax,ubind
      write(*,*)testk,testzy,testay,testyp,testyn
   70 Continue
   90 format(A25,3(1pe12.4),I3)
c
c--solve for actual rho and Ye, calcuate binding energy, average A and Z 
c
      ider=2
      call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,
     &             xa,xh,yeh,zbar,abar)
      If(kit.ge.kmax) Then
          write(*,*) 'Final step failure in nserho'
          write(*,*) kit,rhotmp,rho
      Endif
c
c--Overstep in T9 to give initial estimate of dUb/dT9
c
      If(t9.lt.12.d0) Then
          delt9=1d-4
      Elseif(t9.gt.25.d0) Then
          delt9=1.d-1
      Else
          delt9=1d-2
      Endif
   80 t9d=t9+delt9
      ypd=yp
      ynd=yn
      call nsesolv(ider,t9d,rho,ye,ypd,ynd,kit,kmax,ubindd,
     &             xa,xh,yeh,zbar,abar)
      If(kit.ge.kmax) Then
          write(*,*) 'dUb/dT step failure in nserho'
          write(*,*) kit,t9d,t9
          delt9=.5d0*delt9
          goto 80
      Endif
      dubind=(ubindd-ubind)/(t9d-t9)
c
      return
      end
c 
      subroutine nsetemp(t9old,rho,ye,t9,yp,yn,
     1                   xa,xh,yeh,zbar,abar,ubind,dubind)
c*************************************************************
c
c this subroutine figures out the NSE eq. assuming that yp
c and yn were previously know at the SAME density and ye,
c but different temperature
c
c**************************************************************
c
c
      implicit double precision(a-h,o-z)
      parameter (tolnse=1d-5,kmax=10)
c
      common /testnse/ testk,testzy,testay,testyp,testyn
c
c--finding zero point of binding energy
c
      ider=2
      call nsesolv(ider,t9old,rho,ye,yp,yn,kit,kmax,ubind0,
     &             xa,xh,yeh,zbar,abar)
c
      If(kit.ge.kmax) Then
          write(*,*) 'NSE mis-stored entering nsetemp'
          write(*,*) 'T9, rho, ye',t9old,rho,ye
          write(*,*) 'inconsistent with yp, yn',yp,yn
      Endif
      ypold=yp
      ynold=yn
      t9min=min(t9,t9old)
      ider=0
c
c--pick initial delt9
c
      delt9=0.5d0
      if(t9min.lt.12.d0)then
         delt9=0.05d0
      elseif(t9min.gt.20.d0) then
         delt9=2.d0
      end if
      delt9=dsign(delt9,t9-t9old)
      ypold=yp
      ynold=yn
      t9last=t9old
c
c--Begin temp iteration
c
      do i=1,1000
          delt9=dsign(min(dabs(t9-t9last),dabs(delt9)),delt9)
          t9tmp=t9last+delt9
          call nsesolv(ider,t9tmp,rho,ye,yp,yn,kit,kmax,ubind,
     &                   xa,xh,yeh,zbar,abar)
          If(dabs((t9tmp-t9)/t9).lt.tolnse.and.kit.lt.kmax) goto 70
          If (kit.ge.kmax) Then
              delt9=.5d0*delt9
              yp=ypold
              yn=ynold
          Elseif(kit.lt.4) Then
              delt9=2.d0*delt9
              t9last=t9tmp
              ypold=yp
              ynold=yn
          Else
              t9last=t9tmp
              ypold=yp
              ynold=yn
          Endif
      enddo
c
c--did not converge, print out error
c
      write(*,*)'nsetemp(1) did not convege!'
      write(*,*)t9old,t9tmp,t9
      write(*,*)i,delt9,rho
      write(*,*)ye,yp,yn
      write(*,*)kit,kmax,ubind
      write(*,*)abar,zbar
      write(*,*)'before nsestart',yp,yn
      call nsestart(t9,rho,ye,yp,yn)
      write(*,*)'after nsestart',yp,yn
   70 Continue
c
c--solve for actual t9, calculate binding energy, average A, Z
c
      ider=2
      call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,
     &                   xa,xh,yeh,zbar,abar)
      If(kit.ge.kmax) Then
          write(*,*) 'NSEtemp failed for final T9' 
          write(*,*) kit,t9,t9tmp
          STOP
      Endif
c
c--Overstep in T9 to calculate dUb/dT9
c
      If(t9.lt.12.d0) Then
          delt9=1d-4
      Elseif(t9.gt.25.d0) Then
          delt9=1.d0
      Else
          delt9=1d-2
      Endif
   80 t9d=t9+delt9
      ypd=yp
      ynd=yn
      call nsesolv(ider,t9d,rho,ye,ypd,ynd,kit,kmax,ubindd,
     &             xa,xh,yeh,zbar,abar)
      If(kit.ge.kmax) Then
          write(*,*) 'dUb/dT step failure in nsetemp, particle '
          write(*,*) kit,t9d,t9
          delt9=.5d0*delt9
          goto 80
      Endif
      dubind=(ubindd-ubind)/(t9d-t9)
c
c
c--step back in T9, in order to calculate dUbind/dT9
c
c     call nsesolv(ider,t9last,rho,ye,ypold,ynold,kit,kmax,ubindlast,
c    &                   xa,xh,yeh,zbar,abar)
c     If(kit.ge.kmax) Then
c         write(*,*) 'NSEtemp failed for deriv T9last' 
c         write(*,*) kit,t9,t9tmp
c         STOP
c     Endif
c     dubind=(ubindlast-ubind)/(t9last-t9)
      Return
      End
c     
      subroutine slwrap(inpvar,yesl,brydns,pprev,
     1                  psl,usl,dusl,gamsl,etasl,
     2                  ypsl,ynsl,xasl,xhsl,yehsl,abar,xmuh,u2sl)
c******************************************************************
c
c  This is the wrapper routine for the swesty-lattimer eos.
c
c******************************************************************
c
c-- 0.46*(mn-mp-me) + bindfe56 MeV/nucleon
      double precision ushift
c     parameter(ushift=0.46d0*0.783d0+8.7904d0)
      parameter(ushift=0.0d0)
c
      common /typef/ ieos
c
      include 'eos_m4a.inc'
      include 'el_eos.inc'
c
c--these variables are needed but not declared in the include files
      integer sf
      double precision told, pprev, u2sl
      double precision yesl, psl, usl, dusl, gamsl, etasl, ypsl, ynsl
      double precision xasl, xhsl, yehsl, abar, xmuh
c
      ye=dmax1(yesl,0.031d0)
      call inveos(inpvar,told,ye,brydns,1,eosflg,0,sf,
     1            xprev,pprev)
c
      if (sf.ne.1) print *,'inveos fails for particle'
      psl=ptot
      if (ieos.eq.3) then
         usl=utot+ushift
         dusl=dudt
         u2sl=stot
      else
c--entropy as variable of state
         usl=stot
         dusl=dsdt
         u2sl=utot+ushift
      endif
c
      gamsl=gam_s
      etasl=musube/inpvar(1)
      xmuh=muhat
c
c-- free (exterior) nucleon fractions
      ypsl=xprot
      ynsl=xnut
c--mass fraction and proton fraction of heavy nuclei
      xhsl=xh
      yehsl=x
      xasl=xalfa
      abar=a
c
      return
      end
c
      subroutine unit
c************************************************************
c                                                           *
c  this routine computes the transformation between the     *
c  physical units (cgs) and the units used in the code.     *
c  And the value of physical constants in code units        *
c                                                           *
c************************************************************
c
      double precision umass
      double precision usltemp, uslrho, uslu, uslp,u2slu
      double precision ud1,ut1,ue1,ueg1,uec1
      double precision gg1, arad1, bigr1
      double precision uopr, uotemp, uorho1, uotemp1, uou1
      double precision ufoe1, sigma1, sigma2, xsecnn1, xsecne1,fermi
c
      common /typef/ ieos
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
      common /epcap/ betafac, c2cu, c3cu
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /konst/ gg, clight, arad, bigr, xsecnn, xsecne
c
      data pi/3.14159/, ggcgs/6.67e-8/, avo/6.02e23/
      data aradcgs /7.565e-15/, boltzk/1.381e-16/
      data hbar/1.055e-27/, ccgs/3e10/
      data emssmev/0.511e0/, boltzmev/8.617e-11/
      data ergmev /6.2422e5/, sigma1/9d-44/, sigma2/5.6d-45/
      data c2cgs /6.15e-4/, c3cgs /5.04e-10/, fermi/1d-13/
c
c--1) work out code units:
c
c--specifie mass unit (g)
c
      umass=2d33
c
c--specifie distance unit (cm)
c
      udist=1e9
c
c--transformation factor for :
c
c 1a) density
c
      ud1=dble(umass)/dble(udist)**3
      udens=ud1
c
c 1b) time
c
c     ut1=dsqrt(dble(udist)**3/(dble(gg)*umass))
c     utime=ut1
      ut1=1.d1
      utime=10.
c
c 1c) ergs
c
      ue1=dble(umass)*dble(udist)**2/dble(utime)**2
c--uerg will overflow
c      uerg=ue1
c
c 1d) ergs per gram
c
      ueg1=dble(udist)**2/dble(utime)**2
      uergg=ueg1
c
c 1e) ergs per cc
c
      uec1=dble(umass)/(dble(udist)*dble(utime)**2)
      uergcc=uec1
c
c--2) constants
c
c 2a) gravitational
c
      gg1=dble(ggcgs)*umass*dble(utime)**2/(dble(udist)**3)
      gg=gg1
c
c 2aa) velocity of light
c
      clight=dble(ccgs*utime/udist)
c
c 2b) Stefan Boltzmann (note that T unit is 1e9 K here)
c
      utemp=1e9
      arad1=dble(aradcgs)/ue1*dble(udist)**3*dble(utemp)**4
      arad=arad1
c
c 2c) Perfect Gas "R" constant
c
      bigr1=dble(avo)*dble(boltzk)/ue1*umass*dble(utemp)
      bigr=bigr1
c
c 2d) nucleon+nu x-section/4pi, in udist**2/umass/Mev**2
c
      xsecnn1=sigma1*umass*dble(avo)/dble(udist*udist)
      xsecnn=xsecnn1/(4d0*3.14159d0)
c
c 2e) e+nu x-section/4pi, in Enu(Mev)*udist**2/umass/Mev**2/utemp**4
c
      xsecne1=sigma2*umass*dble(ergmev)*(dble(utemp)*dble(boltzk))**4/
     1     (dble(udist)**2*ud1*3.14159d0*(dble(hbar)*dble(ccgs))**3)
      xsecne=xsecne1/(4d0*3.14159d0)
c
c--3a) Conversion to the Ocean eos units from code units:
c     in the ocean eos; density unit=1e7g/cc
c                       temperature unit= 1e9K
c                       energy/mass=1.0e17cgs
c                       energy/vol =1.0e24cgs
c
      uotemp=1d9/dble(utemp)
      uotemp1=1.d0/uotemp
      uopr=1.d24/dble(uergcc)
      uorho1=ud1/1.d7
      uou1=dble(uergg)/1.d17
c
c--3b) Conversion to the Swesty-Lattimer units from code units:
c     
      usltemp=utemp*boltzmev
      uslrho=ud1*avo*fermi**3
      uslp=uec1*fermi**3/1.602d-6
      if (ieos.eq.3) then
c--if u is internal energy
         uslu=dble(ergmev)*dble(uergg)/dble(avo)
         u2slu=dble(uergg)/dble(utemp)/(dble(avo)*dble(boltzk))
      elseif (ieos.eq.4) then
c--if u is specific entropy
         uslu=dble(uergg)/dble(utemp)/(dble(avo)*dble(boltzk))
         u2slu=dble(ergmev)*dble(uergg)/dble(avo)
      endif
c
c--4) common unit2 stuff
c
c 4a) temp code unit in Mev
c
      utmev=utemp*boltzmev
c
c 4b) energy code unit in foes
c
      ufoe1=ue1/1d51
      ufoe=ufoe1
c
c 4c) Mev/nucleon in code units
c
      umevnuc=avo/(ergmev*uergg)
c
c 4e) Mev to errgs times avogadro's number
c
      umeverg=avo*1.602e-6
c
c--5) epcapture betafac, c2, and c3. 
c
      betafac=emssmev/utmev
      c2cu=c2cgs*utime
      c3cu=c3cgs*utime*ergmev
c
   99 return
      end

      subroutine eosfl(rhoi,pri,
     $     ui,u2i,yei,tempi,ifleosi,abari,
     $     xpi,xni,xpfi,p2i,p3i,p4i,ufreezi,
     $     xmuei,xmuhi,etai,temprev,
     $     yeprev,xpprev,xnprev)
c****************************************************************
c
c this subroutine determines what kind of eos to use:
c       eosflg = 1: freeze-out, just Ocean's eos + Coul corr.
c       eosflg = 2: NSE with Raph's routines, + Ocean eos + Coul
c       eosflg = 3: Swesty's eos
c
c**************************************************************
c
      implicit double precision (a-h,o-z)
      parameter (avokb=6.02e23*1.381e-16)
c
      real  udist, udens, utime, uergg, uergcc
      common /units/ umass, udist, udens, utime, uergg, uergcc
      real  utemp,utmev,ufoe,umevnuc,umeverg
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      double precision uopr, uotemp, uorho1, uotemp1, uou1
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      double precision usltemp, uslrho, uslu, uslp, u2slu
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
      common /typef/ ieos
      double precision inpvar(4)
      common /burnop/ t9nse,t9burn,tfreeze
c
      rhoswe=1.0d5
c      rhoswe=1.05d3
c      rhoswe=2.05d3
      t9nse=1.0d1
      t9burn=21.d1
      tfreeze=t9nse-1.0d0
      rhoburn=1.0d5
c
c--entropy conversion factor
      sfac=avokb*utemp/uergg      
c
      rhocgs=rhoi*udens
      if(ifleosi.eq.0.and.tempi.gt.t9burn.and.
     $     rhoi.gt.rhoburn) then
         ifleosi=1
c         print *, ifleosi
      endif
      if(ifleosi.eq.1.and.(tempi.lt.t9burn.or.
     $     rhoi.lt.rhoburn)) then
         ifleosi=0
c         print *, ifleosi
      endif
      if(ifleosi.eq.1.and.tempi.gt.t9nse) then
         call nsestart(tempi,rhocgs,yei,xpi,xni)
c
c--for NSE, add in the nuclear component to the thermal
c--energy to get the total available internal energy
c
         ui=ui+ufreezi
         ufreezi=0.0
         ifleosi=2
         xnprev=xni
         xpprev=xpi
         yeprev=yei
         temprev=tempi
      endif
      if(ifleosi.eq.2.and.tempi.le.tfreeze) then
         ifleosi=1
         call nsetemp(tempi,rhocgs,yei,tempi,xpi,xni,
     1        xai,xhi,yehi,zbari,abari,ubind,dubind)
c
c--for freeze-out, remove the nuclear component from the total
c--energy to get the thermal energy
c
         ui=ui-ubind/uergg
         ufreezi=ubind/uergg
         xnprev=xni
         xpprev=xpi
         yeprev=yei
         temprev=tempi
      elseif(ifleosi.eq.2.and.rhoi.gt.rhoswe) then
         ifleosi=3
c--make call to sleos to get the entropy or intenal energy 
c--at present rho, ye, T
         brydns=rhoi*uslrho
         pprev=xpfi
         inpvar(1)=tempi*usltemp
         inpvar(2)=p2i
         inpvar(3)=p3i
         inpvar(4)=p4i
         call slwrap(inpvar,yei,brydns,pprev,psl,usl,
     $        dusl,gamsl,etai,xpi,xni,xai,xhi,yehi,abari,xmuh,stot)
         xpfi=pprev
         p2i=inpvar(2)
         p3i=inpvar(3)
         p4i=inpvar(4)
         u2i=usl/uslu
         xmuei=etai*tempi
         xmuhi=xmuh/utmev
         xnprev=xni
         xpprev=xpi
         yeprev=yei
         temprev=tempi
      elseif(ifleosi.eq.3.and.rhoi.lt.rhoswe) then
c-- switch back to internal energy variable of state
         call nsestart(tempi,rhocgs,yei,xpi,xni)
         call nsetemp(tempi,rhocgs,yei,tempi,xpi,xni,
     $        xai,xhi,yehi,zbari,abari,ubind,dubind)
         dens=rhoi*uorho1
         abar2=zbari/yei
         call nados(tempi,dens,zbari,abar2,pel,eel,sel,
     $        ptot,etot,stot,dpt,det,dst,dpd,ded,gamsl,etai)
         dens=rhoi
         call coulomb(dens,zbari,yei,ucoul,pcoul)
         ui=ucoul+ubind/uergg+etot/uou1
         ifleosi=2
         xnprev=xni
         xpprev=xpi
         yeprev=yei
         temprev=tempi
      endif
c
      return
      end
c
      subroutine nidecay(dtburn,pmassi,xni56ni,xco56ni,
     $     xfe56ni,udecayi)
c************************************************************
c
c  This subroutine estimates the decay of nickel into 
c  cobalt and changes the nickel/cobalt abundances as well 
c  as calculating the energy (udecayi) injected into the 
c  system due to this decay.
c
c*************************************************************
c
      implicit double precision(a-h,o-z)
      real  udist, udens, utime, uergg, uergcc
      common /units/ umass, udist, udens, utime, uergg, uergcc
c
c--dufacni/co is the energy per gram of nickel/cobalt decayed
c--from Colgate, Petschek, & Kriese 1980 ApJ 237, L81
c--Ni decay is 32.8% of total energy
      dufacni=7.6272d15/dble(uergg)
c--Co decay is 67.2% of total energy
      dufacco=1.5623d16/dble(uergg)
c
c--taudni is the exponential decay time of nickel
c--6.1d half life converted to 10s for an exponential decay
c
      taudni=3.6532d4
c--taudco is the exponential decay time of cobalt
c--77.25d half life converted to 10s for an exponential decay
      taudco=4.6263d5
c
      decexp=dexp(-dt/taudni)
      xni56ni=decexp*xni56ni
      xco56ni=xco56ni+(1.d0-decexp)*xni56ni
      udecayi=(1.d0-decexp)*xni56ni*dufacni
c      if (xni56ni.gt.0.5) then
c         print *, decexp, 1.0-decexp, dt, taudni, 
c     $        xni56ni, xco56ni,udecayi
c      end if
c
      decexp=dexp(-dt/taudco)
      xco56ni=decexp*xco56ni
      xfe56ni=xfe56ni+(1.d0-decexp)*xco56ni
      udecayi=udecayi+(1.d0-decexp)*xco56ni*dufacco
c
      return
      end
c
c
c
      subroutine burner(dtburn,dtbrni,tempi,rhoi,ui,
     $     ubi,dui,abari,pmassi,etai,ifleosi,
     $     xhe4ni,xc12ni,xo16ni,xne20ni,xmg24ni,
     $     xsi28ni,xni56ni,xpi,xni,iexpl)
c**************************************************************
c
c This is the subcycling burning routine
c
c**************************************************************
c
      implicit double precision(a-h,o-z)
c
c--mass conservation tolerance
c
      parameter (xoff=1.d-3)
c
      real  udist, udens, utime, uergg, uergcc
      common /units/ umass, udist, udens, utime, uergg, uergcc
      real  utemp,utmev,ufoe,umevnuc,umeverg
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /burnop/ t9nse,t9burn,tfreeze
      common /tmpc/ tgo 
      real dtburn

      dtimecgs=dble(dtburn*utime)
      ubi=0.0
c      print *, ifleosi, tempi,dtbrni
      if (ifleosi.eq.1) then
         uburn=0.0d0
         tburn=0.0d0
         dtb=min(dtbrni,dble(dtburn))
         rhocgs=dble(udens)*rhoi
         tempold=tempi
         ui0=ui
         yei=0.5d0
         duhydro=dui
         xhe4_old=xhe4ni
         xc12_old=xc12ni
         xo16_old=xo16ni
         xne20_old=xne20ni
         xmg24_old=xmg24ni
         xsi28_old=xsi28ni
         xni56_old=xni56ni
         xp1_old=max(xpi,1.d-30)
         xn1_old=max(xni,1.d-30)
         abar_old=1d0/(xhe4_old/4.d0+xc12_old/12.d0+
     $        xo16_old/16.d0+
     $        xne20_old/20.d0+xmg24_old/24.d0+
     $        xsi28_old/28.d0+xni56_old/56.d0+
     $        xp1_old/1.d0+xn1_old/1.d0)
c         print *, 'calling rootemp2u',tempold,rhoi,ui0
         call rootemp2u(rhoi,ui0,tempold,yei,abar_old,
     $        ptot,cs,etai,stot)
c--return here if failed
 10      continue
         tempin=tempold
         xhe4=xhe4_old
         xc12=xc12_old
         xo16=xo16_old
         xne20=xne20_old
         xmg24=xmg24_old
         xsi28=xsi28_old
         xni56=xni56_old
         xp1=xp1_old
         xn1=xn1_old
         abari=abar_old
         if (tempin.gt.dble(1.5)) goto 200
c--return here if success
 20      continue
         tempk=dble(utemp)*tempin
         dtbcgs=dble(utime)*dtb
c         print *, 'calling burnwrap',tempold,rhocgs,ui0,tburn
         call burnwrap(tburn,dtimecgs,dtbcgs,rhocgs,tempk,abari,
     $         xhe4,xc12,xo16,xne20,xmg24,xsi28,xni56,
     $         xp1,xn1,dubdtcgs,isucc,xsum)
c         print *, 'out of burnwrap',tempold,rhoi,ui0,tburn
c         if (tgo.gt.1953) read *, dttx
         dtime=dtimecgs/utime
         dubdt=dubdtcgs*dble(utime/uergg)
c--check for abundance variations
         dxhe4=dabs(xhe4-xhe4_old)/dmax1(xhe4_old,xoff)
c         print *,dxhe4,xhe4,xhe4_old,xoff
         dxc12=dabs(xc12-xc12_old)/dmax1(xc12_old,xoff)
         dxo16=dabs(xo16-xo16_old)/dmax1(xo16_old,xoff)
         dxne20=dabs(xne20-xne20_old)/dmax1(xne20_old,xoff)
         dxmg24=dabs(xmg24-xmg24_old)/dmax1(xmg24_old,xoff)
         dxsi28=dabs(xsi28-xsi28_old)/dmax1(xsi28_old,xoff)
         dxni56=dabs(xni56-xni56_old)/dmax1(xni56_old,xoff)
         dxp1=dabs(xp1-xp1_old)/dmax1(xp1_old,xoff)
         dxn1=dabs(xn1-xn1_old)/dmax1(xn1_old,xoff)
         dxmx=max(dxhe4,dxc12,dxo16,dxne20,dxmg24,dxsi28,dxni56,
     $        xp1,xn1)
c         print *, 'dxmx',dxmx, dxhe4,dxc12,dxo16,dxne20,dxmg24,dxsi28,
c     $        dxni56,xp1,xn1
         uin=ui0+(tburn+dtb)*duhydro+uburn+dtb*dubdt
         abari=1.d0/(xhe4/4.d0+xc12/12.d0+xo16/16.d0
     $         +xne20/20.d0+
     $         xmg24/24.d0+xsi28/28.d0+xni56/56.d0+
     $         xp1/1.d0+xn1/1.d0)


c         print *, 'calling rootemp2u',isucc,dxmx,dtemp,dtb
         call rootemp2u(rhoi,uin,tempin,yei,abari,
     $         ptot,cs,etai,stot)
         dtemp=dabs(tempin-tempold)/tempold
c--bad matter conservation or 
c--if abundance varies > 5% or 
c--dtemp varies more than 3%
         if (isucc.eq.-1.or.dxmx.gt.0.5d0.or.dtemp.gt.0.3d0) then
c--chop time step by 2, fail timestep
            dtb=0.5d0*dtb
            goto 10
         endif
c--check for T--->NSE
 200     continue
         if (tempin.gt.dble(1.5)) then
c--set explosive burning flag
            iexpl=1
c--burn to nickel
            utoni56=-umevnuc*(xhe4*28.297d0/4.d0+
     $            xc12*92.16d0/12.d0+xo16*127.62d0/16.d0+
     $            xne20*160.65d0/20.d0+xmg24*198.26/24.d0+
     $            xsi28*236.54d0/28.d0+(xni56-1.0d0)*484.0/56.d0)
            ubi=utoni56+uburn
c            ubi=0.
            xhe4ni=1.0d-30
            xc12ni=1.0d-30
            xo16ni=1.0d-30
            xne20ni=1.0d-30
            xmg24ni=1.0d-30
            xsi28ni=1.0d-30
            xni56ni=1.0d0
            xpi=1.0d-30
            xni=1.0d-30
            abari=56.0
            dtbrni=1.d-9
            abari=56.
            uiold=uin
            uin=ui0+ubi+dtime*duhydro
            call rootemp2u(rhoi,uin,tempin,yei,abari,
     $           ptot,cs,etai,stot)
         else
c--time step accepted update tburn
c            print *, 'tburn',tburn, dtb
            tburn=tburn+dtb
c            print *, 'tburn',tburn, dtb
            uburn=uburn+dtb*dubdt
            tempold=tempin
            xhe4_old=xhe4
            xc12_old=xc12
            xo16_old=xo16
            xne20_old=xne20
            xmg24_old=xmg24
            xsi28_old=xsi28
            xni56_old=xni56
            xp1_old=xp1
            xn1_old=xn1
            abar_old=abari
c--if good matter conservation and abundances vary < 2% and 
c--temperature < 1%
            if((isucc.eq.1).and.(dxmx.lt.0.02d0).and.
     $            (dtemp.lt.0.01d0)) then
               dtb=2.0*dtb
            endif
c
            tleft=dtime-tburn
c            print *, 'tleft=',tleft,dtime,tburn,dtb
            if (tleft.gt.(1d-5*dtime)) then
               if (dtb.lt.tleft) then
                  goto 20
               else
                  dtbrni=real(dtb)
                  dtb=tleft
c                  print *, 'dtb gt tleft',tburn,tleft
                  goto 20
               endif
            endif
c            print *, 'tleft=',tleft,tgo
c--burn energy 
            ubi=uburn
            xhe4ni=xhe4
            xc12ni=xc12
            xo16ni=xo16
            xne20ni=xne20
            xmg24ni=xmg24
            xsi28ni=xsi28
            xni56ni=xni56
            xpi=xp1
            xni=xn1
         endif
      endif
c      print *, 'done'
c
      return
      end
c
      subroutine burnwrap(tburn,dtimecgs,dtbcgs,rhocgs,tempk,abar,
     1                    xhe4,xc12,xo16,xne20,xmg24,xsi28,xni56,
     2                    xp1,xn1,dubdtcgs,isucc,xsum)
c************************************************************
c
c this is the wrapper routine for Stan & Rob's 9 isotope net
c
c************************************************************
c
      implicit double precision(a-h,o-z)
c
      parameter (xoff=1d-2)
c
      common /iso9d/dens,temp,dt,
     1       rnep,time
      common/xnuc/ xn(9),energy
      common/scrn/sfdpg,sf3a,sfcag,sfoag,sfneag,sfmgag,sfcaag,sf1212,
     1 sf1216,sf1616
c
      dens=rhocgs
      temp=tempk
      dt=dtbcgs
      rnep=0.0d0
      time=dtimecgs
      xn(1)=xc12
      xn(2)=xo16
      xn(3)=xne20
      xn(4)=xmg24
      xn(5)=xsi28
      xn(6)=xni56
      xn(7)=xhe4
      xn(8)=xn1
      xn(9)=xp1
      zbar=0.5d0*abar
      z2bar=xc12*(6.0d0*6.0d0/12.0d0)+xo16*(8.0d0*8.0d0/16.0d0)+
     1      xne20*(10.0d0*10.0d0/20.0d0)+
     2      xmg24*(12.0d0*12.0d0/24.0d0)+
     3      xsi28*(14.0d0*14.0d0/28.0d0)+
     4      xni56*(28.0d0*28.0d0/56.0d0)+
     5      xhe4*(2.0d0*2.0d0/4.0d0)+xp1
      z2bar=z2bar*abar
      ytot1=1./abar

c
      call screen(temp,dens,abar,zbar,z2bar,ytot1)
c
      call rate(rdpg,rhegp,r3a,rg3a,rcag,roga,roag,rnega,
     1 rneag,rmgga,rmgag,rsiga,rcaag,rtiga,r1212,r1216,r1616)
c      print *, crdpg,rhegp,r3a,rg3a,rcag,roga,roag,rnega, 
c     $     rneag,rmgga,rmgag,rsiga,rcaag,rtiga,r1212,r1216,r1616
      call snuc(rdpg,rhegp,r3a,rg3a,rcag,roga,roag,rnega,
     1 rneag,rmgga,rmgag,rsiga,rcaag,rtiga,r1212,r1216,r1616)
c      print *, xn(1),xn(2),xn(3),xn(4),xn(5),xn(6),xn(7),xn(8),
c     $     xn(9)
c
      xc12=xn(1)
      xo16=xn(2)
      xne20=xn(3)
      xmg24=xn(4)
      xsi28=xn(5)
      xni56=xn(6)
      xhe4=xn(7)
      xn1=xn(8)
      xp1=xn(9)
c      print *, xhe4,xn(7)
c
      xsum=1.d0-(xp1+xn1+xhe4+xc12+xo16+xne20+xmg24+xsi28+xni56)
c      print *, 'xsum', dabs(xsum)
      if (dabs(xsum).lt.xoff) then
         isucc=1
      elseif (dabs(xsum).lt.(10.0d0*xoff)) then
         isucc=0
c--reset largest abundance
         xmx=dmax1(xn1,xp1,xhe4,xc12,xo16,xne20,xmg24,xsi28,xni56)
         if (xo16.ge.xmx) then
            xo16=xo16+xsum
         elseif (xsi28.ge.xmx) then
            xsi28=xsi28+xsum
         elseif (xni56.ge.xmx) then
            xni56=xni56+xsum
         elseif (xmg24.ge.xmx) then
            xmg24=xmg24+xsum
         elseif (xne20.ge.xmx) then
            xne20=xne20+xsum
         elseif (xhe4.ge.xmx) then
            xhe4=xhe4+xsum
         elseif (xc12.ge.xmx) then
            xc12=xc12+xsum
         elseif (xn1.ge.xmx) then
            xn1=xn1+xsum
         else
            xp1=xp1+xsum
         endif
      else
         isucc=-1
      endif
      dubdtcgs=energy
c
      return
      end
c

      subroutine rate(rdpg,rhegp,r3a,rg3a,rcag,roga,roag,rnega,
     1 rneag,rmgga,rmgag,rsiga,rcaag,rtiga,r1212,r1216,r1616)

c- change 1
      implicit double precision(a-h,o-z)
      common/iso9d/dens,temp,dt,
     1       rnep,time
      common/scrn/sfdpg,sf3a,sfcag,sfoag,sfneag,sfmgag,sfcaag,sf1212,
     1 sf1216,sf1616

c.... rates to be evaluated
c.... * rdpg  = 2h(p,g)3he       rhegp = 3he(g,p)2h
c.... * r3a = 3-alpha reaction   rg3a  = 12c(g,a)2alpha
c.... * rcag  = 12c(a,g)16o      roga  = 16o(g,a)12c
c.... * roag  = 16o(a,g)ne20     rnega = 20ne(g,a)16o 
c.... * rneag = 20ne(a,g)24mg    rmgga = 24mg(g,a)20ne
c.... * rmgag = 24mg(a,g)28si    rsiga = 28si(g,a)24mg  
c.... * rcaag = 40ca(a,g)44ti    rtiga = 44ti(g,a)40ca
c.... * r1212 = c12+c12                   
c.... * r1216 = c12+o16          omit inverses for these 3        
c.... * r1616 = o16+o16                   

c.. in terms of screening, the rates with a * above are to be screened
c.. screen factor come in from common /scrn/ 

      t9r= temp/1.0e+09
      t9= min(10.d00,t9r)
      t92= t9*t9
      t912= sqrt(t9)
      t913= t9**(1./3.)
      t923= t913*t913
      t943= t9*t913
      t953= t9*t923
      t932= t9*t912
      t93= t9*t9*t9
      t9i= 1./t9
      t9i13= 1./t913
      t9i23= 1./t923
      t9i32= 1./t932
      t9ri= 1./t9r
      t9r32= t9r*sqrt(t9r)
      theta= 0.1
c
c.... reaction rates
c
c      triple alpha reaction rate
c      rate is for formation of c12 compound nucleus
c      (has been divided by 6)
c
c.... rate revised according to caughlan et al. oap400 june 1984.
c
      r2abe= (7.40d+05*t9i32)*exp(-1.0663*t9i)
     1      +3.714d+09*t9i23*exp(-13.49*t9i13-(t9/0.098)**2)
     2      *(1.+0.031*t913+8.431*t923+1.824*t9+53.313
     3      *t943+29.311*t953)
      rbeac= (130./t932)*exp(-3.3364*t9i)
     1      +3.27d+07*t9i23*exp(-23.57*t9i13-(t9/0.235)**2)
     2      *(1.+0.018*t913+4.708*t923+0.583*t9+16.621
     3      *t943+5.23*t953)
      r3a= 2.90d-16*(r2abe*rbeac)
     1    +0.1*1.35d-07*t9i32*exp(-24.811*t9i)
      r3a= sf3a*r3a*(dens*dens)/6.
      rev= 2.0d+20*exp(-84.424*t9ri)
      rg3a= rev*(t9r**3)*6.*r3a/(dens*dens)
c
c.... r1212= c12+c12                                     fczii(75)
c
      t9a= t9/(1.+0.067*t9)
      t9a13= t9a**(1./3.)
      t9a23= t9a13*t9a13
      t9a56= t9a/sqrt(t9a13)
      denom= exp(-0.010*t9a*t9a*t9a*t9a)+5.56e-03*exp(1.685*t9a23)
      r1212= sf1212*0.5*1.26e+27*(t9a56/(t932*denom))
     1      *exp(-84.165/t9a13)*dens
c
c.... r1216 = c12+o16                                    fczii(75)
c
      t9a= t9/(1.+0.055*t9)
      t9a13= t9a**(1./3.)
      t9a23= t9a13*t9a13
      t9a56= t9a/sqrt(t9a13)
      denom= exp(-0.180*t9a*t9a)+1.06e-03*exp(2.562*t9a23)
      r1216= sf1216*1.72e+31*(t9a56/(t932*denom))
     1      *exp(-106.594/t9a13)*dens
c
c.... r1616= o16+o16                                    fczii(75)
c
      t9a= t9/(1.+0.067*t9)
      t9a13= t9a**(1./3.)
      t9a23= t9a13*t9a13
      t9a56= t9a/sqrt(t9a13)
      denom= exp(-0.032*t9a*t9a*t9a*t9a)+3.89e-04*exp(2.659*t9a23)
      r1616= sf1616*0.5*3.61e+37*(t9a56/(denom*t932))
     1      *exp(-135.930/t9a13)*dens
c
c        now do rates for alpha chain
c             
c... 12c(a,g)16o cf88
c.... updated as specified by caughlan & fowler tnrr5 1988
c    q = 7.162
c     
      rcag= sfcag*dens*(1.04e+08/t9**2/(1.+0.0489/t923)**2
     1        *exp(-32.120/t913-(t9/3.496)**2)
     2        +1.76e+08/t9**2/(1.+0.2654/t923)**2*exp(-32.120/t913)
     3        +1.25e+03/t932*exp(-27.499/t9)+1.43e-02*t9**5
     4        *exp(-15.541/t9))
c  see comment in cfhz85 
      roga= rcag/dens*5.13e+10*t9r32*exp(-83.111*t9ri)
c
c.... 16o(a,g)ne20
c.... updated as specified by caughlan et al oap400 june 1984. 
c

c... 16o(a,g)20ne cf88
c.... updated as specified by caughlan & fowler tnrr5 1988
c..  q = 4.734
      roag= sfoag*dens*(9.37e+09/t923*exp(-39.757/t913-(t9/1.586)**2)
     1+6.21e+01/t932*exp(-10.297/t9)+5.38e+02/t932*exp(-12.226/t9)
     2+1.30e+01*t9**2*exp(-20.093/t9))
      rnega= roag*5.65e+10*t9r32*exp(-54.937*t9ri)/dens
c
c.... 20ne(ag)24mg + inverse                            fczii(75)
c
      term1= 4.11e+11*t9i23*exp(-46.766*t9i13-(t9/2.219)**2)
     1      *(1.+0.009*t913+0.882*t923+0.055*t9+0.749*t943
     2      +0.119*t953)
      term2= 5.27e+03*t9i32*exp(-15.869*t9i)+6.51e+03*t912
     1      *exp(-16.223*t9i)
      term3= 4.21*t9i32*exp(-9.115*t9i)+3.2*t9i23*exp(-9.383*t9i)
      rneag= sfneag*dens*(term1+term2+term3)/(1.+5.*exp(-18.960*t9i))
      rev=   6.01e+10*t9r32*exp(-108.094*t9ri)
      rmgga= (rneag/dens)*rev
c
c.... 24mg(a,g)28si                                    fcz3(83)
c
      gt9cd= (1.+5.*exp(-15.882*t9i))
      rmgag= sfmgag*dens*(47.8*t9i32*exp(-13.506*t9i)+2.38d+03*t9i32
     1      *exp(-15.218*t9i)+2.47d+02*t932*exp(-15.147*t9i)
     2      +theta*(1.72d-09*t9i32*exp(-5.028*t9i)+1.25d-03*t9i32
     3      *exp(-7.929*t9i)+2.43d+01*t9i*exp(-11.523*t9i)))/gt9cd
      rsiga=rmgag*6.27d+10*t9r32*exp(-1.15882d+02*t9ri)/dens
c
c     include rate necessary for alpha photodisintegration and coupling
c     si to ni.
c 
c.... 2d(pg)3he                                        fczii(75)
c
      term= 2.65e+03*t9i23*exp(-3.720*t9i13)
     1     *(1.+0.122*t913+1.99*t923+1.56*t9+0.162*t943
     2     +0.324*t953)
      rev= 1.63e+10*t9r32*exp(-63.755*t9ri)
      rdpg= sfdpg*dens*term
      rhegp= rev*term
c
c
c.... 40ca(ag)44ti + inverse                           woosley(78)
c
      term=4.66e+24*t9i23*exp(-76.435*t9i13*(1.+1.650e-02*t9
     1     +5.973e-03*t92-3.889e-04*t93))
      rev= 6.843e+10*t9r32*exp(-59.510*t9ri)
      rcaag= sfcaag*dens*term
      rtiga= rev*term
c
      return
      end

      subroutine screen(t,d,abar,zbar,z2bar,ytot1)
c      implicit real*8 (a-h,o-z)
      implicit double precision (a-h,o-z)
      save
c.... revised 10/94 to screen only d(p,g), various (a,g), and heavy ion rates 
c.... for 9iso application - screen factors applied in subr rate
c.... and passed via common /scrn/
c.... revised august 10, 1992 to speed up computation
c....
c.... last earlier revision 11 nov 1982 .
c....
c..   this subroutine calculates screening factors for nuclear reaction
c..   rates in the weak, intermediate , and strong regimes
c..   given the temperature (t--degk), the density (d--g/cc), the
c..   atomic numbers and weights of the elements in the reaction
c..   channel with the largest coulomb barrier (z1,z2,a1,a2),
c..   and the mean plasma parameters
c..   calculated in main and passed in arg list are:
c..   (mean atomic number--zbar, mean square of the atomic number--z2bar
c..   mean atomic weight--abar, and total number of moles of
c..   nuclei per gram--ytot1).
c..   the unscreened rate is to be multiplied by the dimensionless
c..   output parameter scfac to get the corrected rate.
c..   the treatment is based on graboske, dewit, grossman, and cooper
c..   ap j. 181,457 (1973) for weak screening and on
c..   alastuey and jancovici, ap.j. 226, 1034, 1978, with plasma
c..   parameters from itoh, totsuji, setsuo, and dewitt, ap.j. 234,
c..   1079,1979,  for strong screening (rkw modification).

      parameter (nscf=10)
      parameter (nscf1=11)
      parameter (x13=1.d+0/3.d+0,x14=1.d+0/4.d+0)

      dimension scfac(nscf),rz1(nscf1),rz2(nscf1),ra1(nscf1),
     1          ra2(nscf1),scfac1(nscf1),zhat(nscf1),zhat2(nscf1),
     2          gamefac(nscf1),tau12fac(nscf1),cfac(nscf1)

      common/scrn/sfdpg,sf3a,sfcag,sfoag,sfneag,sfmgag,sfcaag,sf1212,
     1 sf1216,sf1616

      data  rz1 /1. ,2. ,2.  ,2.  ,2.  ,2.  ,2.  ,6.  ,6.  ,8.  ,2./
      data  ra1 /1. ,4. ,4.  ,4.  ,4.  ,4.  ,4.  ,12. ,12. ,16. ,4./
      data  rz2 /1. ,4. ,6.  ,8.  ,10. ,12. ,20. ,6.  ,8.  ,8.  ,2./
      data  ra2 /2. ,8. ,12. ,16. ,20. ,24. ,40. ,12. ,16. ,16. ,4./

      data ipass / 0 /

c.... initialize constants

      if (ipass.gt.0) go to 20
      theta=1.
      x53=5./3.
      x512=5./12.

c....
c.... calculate individual screening factors
c.... approx. for strong screening only good for alpha .lt. 1.6

c.. loop0 (label 10)
      do 10 i1=1,nscf1
      z1=rz1(i1)
      z2=rz2(i1)
      a1=ra1(i1)
      a2=ra2(i1)
      zhat(i1)=(z1+z2)**x53-z1**x53-z2**x53
      zhat2(i1)=(z1+z2)**x512-z1**x512-z2**x512
      gamefac(i1)=2.**x13*z1*z2/(z1+z2)**x13
      tau12fac(i1)=(z1**2*z2**2*a1*a2/(a1+a2))**x13
      cfac(i1)=x53*log(z1*z2/(z1+z2))
10    continue
      ipass=1

c.... calculate average plasma parameters

c..strt (label 20)  
20    qlam0=1.88d+8*sqrt(d/(abar*t*t*t))
      ztilda=sqrt(z2bar+zbar*theta)
      qlam0z=qlam0*ztilda
      gamp=2.27493d+5*(d*zbar*ytot1)**x13/t
      taufac=4.248719d+3/t**x13

c.... calculate individual screening factors
c.... approx. for strong screening only good for alpha .lt. 1.6

c.. loop1 (label 30)
      do 30 i1=1,nscf1
      gamef=gamp*gamefac(i1)
      tau12=taufac*tau12fac(i1)
      alph12=3.*gamef/tau12
c....
c.... limit alph12 to 1.6 to prevent unphysical behavior
c.... (h dec. as rho inc.) at high rho.  this should really
c.... be replaced by a pycnonuclear reaction rate formula.
c....
c.. ncor (label 40)
      if(alph12.le.1.6) go to 40
      alph12=1.6
      gamef=1.6*tau12*x13
      gamp=gamef/gamefac(i1)

40    h12w=rz1(i1)*rz2(i1)*qlam0z
c.. weak (lable 41)
41    h12=h12w
      if(gamef.gt.0.3) go to 50
      go to 60

c.. strng (lable 50)
50    c=0.896434*gamp*zhat(i1)-3.44740*gamp**x14*zhat2(i1)-
     1  0.5551*(log(gamp)+cfac(i1))-2.996

      alph122=alph12*alph12
      alph123=alph12*alph122
      alph124=alph12*alph123
      alph125=alph12*alph124
      alph126=alph12*alph125
      h12=c-(tau12/3.)*(5.*alph123/32.-0.014*alph124
     1    -0.0128*alph125)-gamef*(0.0055*alph124
     2    -0.0098*alph125+0.0048*alph126)

      xlgfac=1.-0.0562*alph123
      if (xlgfac.lt.0.77) xlgfac=0.77
      h12=log(xlgfac)+h12
      if(gamef.gt.0.8) go to 60
      h12=h12w*((0.8-gamef)*2.0)+h12*((gamef-0.3)*2.0)

c.. scalc (lable 60)
60    if(h12.gt.300.) h12=300.
      if(h12.lt.0.) h12=0.
      scfac1(i1)=exp(h12)

30    continue

c.... complete triple alpha screening factor

      scfac1(2)=scfac1(2)*scfac1(nscf1)
      scfac1(2)=min(scfac1(2),1.0d+130)

c.. loop2 (lable 70)
      do 70 i1=1,nscf
70    scfac(i1)=scfac1(i1)
c.. assign to specific screen factor names.
c.. passed in common /scrn/ to subr rate
      sfdpg=scfac(1)
      sf3a=scfac(2)
      sfcag=scfac(3)
      sfoag=scfac(4)
      sfneag=scfac(5)
      sfmgag=scfac(6)
      sfcaag=scfac(7)
      sf1212=scfac(8)
      sf1216=scfac(9)
      sf1616=scfac(10)
      return
      end

      subroutine leqs(n)
c      implicit real*8 (a-h,o-z)
      implicit double precision (a-h,o-z)
c
c.... leqs performs matrix inversion
c
c     n-dimension of a
      common /acl/ amat(9,9),bvec(9)
c
      n1=n-1
c     find maximum element in each row, and divide
      do 1000 i=1,n
         r=abs(amat(i,1))
         do 100 l=2,n
            r=max(r,abs(amat(i,l)))
  100    continue
         ri=1.d0/r
         do 200 l=1,n
            amat(i,l)=amat(i,l)*ri
  200    continue
         bvec(i)=bvec(i)*ri
 1000 continue
c
      do 2001 k=1,n1
         l=k+1
         damatkk=1.0d0/amat(k,k)
         do 2000 i=l,n
            r=-amat(i,k)*damatkk
            do 300 m=l,n
               amat(i,m)=amat(i,m) + r*amat(k,m)
  300       continue
            bvec(i)=bvec(i) + r*bvec(k)
 2000    continue
 2001 continue
c
c	the matrix is now in upper triangular form
c	start the back substitution and find the solution
c
      bvec(n)=bvec(n)/amat(n,n)
c
      do 3000 i1=1,n1
         i=n-i1
         r=0.
         do 500 k1=1,i1
            k= n + 1 - k1
            r=r + amat(i,k)*bvec(k)
  500    continue
         bvec(i)=(bvec(i)-r)/amat(i,i)
 3000 continue
c
      return
      end


      subroutine snuc(rdpg,rhegp,r3a,rg3a,rcag,roga,roag,rnega,
     1 rneag,rmgga,rmgag,rsiga,rcaag,rtiga,r1212,r1216,r1616)

c.... should run double precision on vax or ibm equipment
c.... compile g-float if available on vaxes

c      implicit real*8(a-h,o-z)
      implicit double precision(a-h,o-z)

c.... input via common: mass fractions of 12c, 16o, 20ne, 
c....       24mg, 28si, 56ni, alpha, n, and p; elapsed 
c....       time; time step; temperature and density. 
c....       new mass fractions and energy generation rate 
c....       in erg/(g sec) calculated.

      common /acl/ ab(9,9),bb(9)
      common/iso9d/dens,temp,dt,
     1       rnep,time
      common/xnuc/ xn(9),enuc

      dimension yn(9),be(9),zn(9),an(9)
      data be/92.16294,127.62093,160.64788,198.2579,
     1        236.5379,483.98,28.29603,0.,0./
      data an/12.,16.,20.,24.,28.,56.,4.,1.,1./
      data zn/6.,8.,10.,12.,14.,28.,2.,0.,1./
c
c.... user should provide subroutine to generate reaction 
c.... rates as a function of temperature, t9 (t/1.e+09 k),
c.... and density, dens (gm/cc)

c.... rates to be evaluated
c.... * rdpg  = 2h(p,g)3he       rhegp = 3he(g,p)2h
c.... * r3a = 3-alpha reaction   rg3a  = 12c(g,a)2alpha
c.... * rcag  = 12c(a,g)16o      roga  = 16o(g,a)12c
c.... * roag  = 16o(a,g)ne20     rnega = 20ne(g,a)16o 
c.... * rneag = 20ne(a,g)24mg    rmgga = 24mg(g,a)20ne
c.... * rmgag = 24mg(a,g)28si    rsiga = 28si(g,a)24mg  
c.... * rcaag = 40ca(a,g)44ti    rtiga = 44ti(g,a)40ca
c.... * r1212 = c12+c12                   
c.... * r1216 = c12+o16          omit inverses for these 3        
c.... * r1616 = o16+o16                   

c.. in terms of screening, the rates with a * above are to be screened
c.. screen factors calculated in subr screen and applied in subr rate

c.... need to define t9, rates passed in arg list above
      t9=temp/1.0e+09
c.... start stans snuc from here down
      t9i=1./t9
      t932=t9*sqrt(t9)
      t9i32=1./t932
      ft9= 0.
      gt9= 0.
      ht9= 0.
      xjt9= 0.
      tauq= 4.0e-2*((t9/3.6)**(-33.3))
      if (time.le.tauq) go to 10
      ft9= (t9i32**3)*exp(239.42*t9i-74.741)
      gt9= (t932**3)*exp(-274.12*t9i+74.914)
      if (t9.lt. 3.0) go to 10
      ht9= t9i32*exp(25.762*t9i-22.260)
      xjt9= t932*exp(-238.79*t9i+23.967)
 10   continue      

      do 20 i=1,9
  20  yn(i)=xn(i)/an(i)
c
      dti=1./dt
c
      w1= yn(7)*rcag
      w2= yn(1)*r1212
      w3= yn(2)*r1216
      w4= yn(1)*r1216
      w5= yn(2)*r1616
      w6= yn(7)*roag
      w7= 0.
      w8= yn(7)*rneag
      w9= 0.
      w10= yn(7)*rmgag
      w11= 0.
      w12= 0.
      w13= yn(7)*yn(7)*r3a
      w14= yn(1)*rcag
      w15= 0.
      w16= 0.
      w17= yn(2)*roag
      w18= 0.
      w19= yn(3)*rneag
      w20= yn(4)*rmgag
      w21= 0.
      w22= 0.
      w23= 0.
      w24= 0.
      w25= 0. 
      w26= 0.
      w27= 0.        
      w28= 0.
      w29= 0.
      if (yn(7).eq.0.) go to 30
      w24= ft9*(dens*yn(7))**3*rcaag
      w25= w24*yn(5)
      w26= w24*yn(7)
      w27= gt9*rtiga/((dens*yn(7))**3)
      w28= w27/yn(7)
      w29= w28*yn(6)
 30   continue
      if (yn(9).le.0.) yn(9)=1.0e-30
      w30= dens*yn(9)*yn(8)*ht9*rdpg
      w31= yn(8)*w30
      w32= w31/yn(9)
      w33= xjt9*rhegp/(dens*yn(9))
      w34= w33*yn(7)
      w35= w34/yn(9)
c  
c.... construct the matrix elements
c
      do 40 i=1,9
      do 41 j=1,9
      ab(i,j)=0.
 41   continue
 40   continue
c
c.... carbon
c
      ab(1,1)= dti+rg3a+w1+4.*w2+w3
      ab(1,2)= -roga+w4
      ab(1,7)= -3.*w13+w14
c
c.... oxygen
c
      ab(2,1)= -w1+w3
      ab(2,2)= dti+roga+4.*w5+w4+w6
      ab(2,3)= -rnega
      ab(2,7)= w17-w14
c
c.... neon
c
      ab(3,1)= -2.*w2
      ab(3,2)= -w6
      ab(3,3)= dti+rnega+w8
      ab(3,4)= -rmgga
      ab(3,7)= -w17+w19
c
c.... magnesium
c
      ab(4,1)= -w3
      ab(4,2)= -w4
      ab(4,3)= -w8
      ab(4,4)= dti+rmgga+w10
      ab(4,5)= -rsiga
      ab(4,7)= -w19+w20
c
c.... silicon
c
      ab(5,2)= -2.*w5
      ab(5,4)= -w10
      ab(5,5)= dti+w26+rsiga
      ab(5,6)= -w27
      ab(5,7)= -w20+4.*w25+3.*w29
c
c.... nickel
c
      ab(6,5)= -w26
      ab(6,6)= dti+w27
      ab(6,7)= -4.*w25-3.*w29
c
c.... helium
c
      ab(7,1)= -3.*rg3a+w1-2.*w2-w3
      ab(7,2)= -roga+w6-2.*w5-w4
      ab(7,3)= -rnega+w8
      ab(7,4)= -rmgga+w10
      ab(7,5)= -rsiga+7.*w26
      ab(7,6)= -7.*w27
      ab(7,7)= dti+9.*w13+w14+w17+w19
     1       +w20+28.*w25+21.*w29+w33
      ab(7,8)= -2.*w30
      ab(7,9)= -w32-w35
c
c.... neutrons
c
      ab(8,7)= -2.*w33
      ab(8,8)= dti+rnep+4.*w30
      ab(8,9)= 2.*w32+2.*w35
c
c.... protons
c
      ab(9,7)= -2.*w33
      ab(9,8)= -rnep+4.*w30
      ab(9,9)= dti+2.*w32+2.*w35
c
      z1= yn(7)*w13
      z2= yn(1)*rg3a
      z3= yn(2)*roga
      z4= yn(1)*w1
      z5= yn(1)*w2
      z6= 0.
      z7= yn(1)*w3
      z8= 0.
      z9= yn(2)*w5
      z10= 0.
      z11= yn(2)*w6
      z12= yn(3)*rnega
      z13= yn(4)*rmgga
      z14= yn(3)*w8
      z15= yn(5)*rsiga
      z16= yn(4)*w10
      z17= w25*yn(7)-w27*yn(6)
      z18= yn(8)*w30-w34
c
      bb(1)= z1-z2+z3-z4+2.*(-z5)-z7
      bb(2)= -z3+z4-z7+z8+2.*(-z9)-z11+z12
      bb(3)= z5+z11-z12+z13-z14
      bb(4)= z7-z13+z14+z15-z16
      bb(5)= z9-z15+z16-z17
      bb(6)= z17
      bb(7)= 3.*(z2-z1)+z3-z4+z5+z7+z9-z11+z12
     1       +z13-z14+z15-z16-7.*z17+z18
      bb(8)= -2.*z18
      bb(9)= -2.*z18
c      print *, 'bb8',bb(8),z18,w30,w34
c
c.... perform matrix inversion to obtain solution.
c....      for low temperatures (t9 .lt. 3) invert
c....      only a subset of the matrix for speed.
c.... user supply standard (gaussian elimination) 
c.... matrix inverter. will contain common /acl/. 
c.... returns new vector in bb.
c
      mdim= 9
c      if (t9.le.3.0) mdim= 7
      call leqs(mdim)
c      print *, 'bb8',bb(8)
c
c.... calculate energy generation rate
c
      enuc=0.
      do 50 i=1,9
 50   enuc=enuc+be(i)*bb(i)
      enuc=enuc*9.647e+17/dt
c      print *, 'bb8',bb(8)
c
c.... calculate mass fractions
c
c      print *, 'bb8',bb(8)
      do 60 i=1,9
c         print *, 'bb,an',i,bb(i),an(i)
 60   xn(i)=xn(i)+an(i)*bb(i)
c
      return
      end






