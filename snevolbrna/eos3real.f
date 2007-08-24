      subroutine eossetup
c
      common /typef/ ieos
c
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
      subroutine eos3(rhoi,ui,u2i,yei,tempi,ifleosi,abari,xpi,xni,
     $     xpfi, p2i, p3i, p4i,temprev,rhoprev,
     $     xpprev,xnprev,yeprev,ufreezi)
c*************************************************************
c
c     compute pressures and temperatures with the    
c     Ocean eos assuming NSE
c
c************************************************************
c
      double precision umass, uinput
      double precision rhoi, ui, u2i, tempi, yei, ptot, pri, cs, etai, 
     $     abari, xpi, xni, xai, xhi, yehi, rhoold, yeold,
     $     xpfi, p2i, p3i, p4i, ufreezi, xmuhi, stot,vsoundi,
     $     temprev,rhoprev,xpprev,xnprev,yeprev,xmuei,xmuhati,
     $     xalphai,xheavyi
c
      common /output/ vsoundi,pri,etai,yehi,xmuei,xmuhati,xalphai,
     $     xheavyi
c
      common /konst/ gg, clight, arad, bigr, xsecnn, xsecne
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /typef/ ieos
c
      if (yei.lt.0.) then
         print *,'yei',yei
         stop
      endif
c
c--assume chemical freeze-out
c
      ui=max(1.e-5,ui)
      call eosfl(rhoi,pri,
     $     ui,u2i,yei,tempi,ifleosi,abari,
     $     xpi,xni,xpfi,p2i,p3i,p4i,ufreezi,
     $     xmuei,xmuhi,etai,temprev,
     $     yeprev,xpprev,xnprev)
      if (ifleosi.eq.1) then
         call rootemp2(rhoi,ui,tempi,yei,abari,
     1        ptot,cs,etai,stot)
         xpi=0.
         xni=0.
         xmuei=etai*tempi
         xmuhati=0.0
         xalphai=0.0
         xheavyi=1.
         yehi=yei
c
c--ocean eos + nse
c
      elseif (ifleosi.eq.2) then
         rhoold=rhoprev*udens
         yeold=yeprev
         call rootemp3(rhoi,ui,temprev,yei,rhoold,yeold,
     $        ptot,cs,etai,xpprev,xnprev,xai,xhi,yehi,
     $        abari,stot)
         xalphai=xai
         xheavyi=xhi
         xmuei=etai*tempi
         xmuhati=0.0
         pri=ptot
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
      if (xpi.le.1d-15) xpi=0.
      if (xni.le.1d-15) xni=0.
      vsoundi=min(cs,0.3333*dble(clight))
      pri=ptot
      temprev=tempi
      rhoprev=rhoi
      xnprev=xni
      xpprev=xpi
      yeprev=yei
c         
      return
      end
c
      subroutine rootemp2(rhoi,ui,tempi,yei,
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
     1           ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)
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
     1              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)
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
     1           ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)
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
     1              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)
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
c
      rhoswe=1.05d3
c      rhoswe=2.05d3
      t9nse=21.d0
      tfreeze=t9nse-7.0d0
c
c--entropy conversion factor
      sfac=avokb*utemp/uergg      
c
      rhocgs=rhoi*udens
      if(ifleosi.eq.1.and.(tempi.gt.t9nse.or.rhoi.gt.4.0*rhoswe)) then
c	 print *,tempi,rhocgs,yei,xpi,xni,t9nse
         call nsestart(tempi,rhocgs,yei,xpi,xni)
c
c--for NSE, add in the nuclear component to the thermal
c--energy to get the total available internal energy
c
c	 print *,tempi,rhocgs,yei,xpi,xni,ufreezi
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
c	 print *, tempi,rhocgs,yei,xpi,xni,tfreeze
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
     $        ptot,etot,stot,dpt,det,dpd,ded,gamsl,etai)
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
