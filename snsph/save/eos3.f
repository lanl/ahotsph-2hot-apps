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
     $     xpprev,xnprev,yeprev,ufreezi,iident,iprocnum)
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
      call eosfl(rhoi,pri,
     $     ui,u2i,yei,tempi,ifleosi,abari,
     $     xpi,xni,xpfi,p2i,p3i,p4i,ufreezi,
     $     xmuei,xmuhi,etai,rhoprev,temprev,
     $     yeprev,xpprev,xnprev,iident,iprocnum)
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
         u2i=stot
c
c--ocean eos + nse
c
      elseif (ifleosi.eq.2) then
         rhoold=rhoprev*udens
         yeold=yeprev
         call rootemp3(rhoi,ui,tempi,yei,rhoold,yeold,
     $        ptot,cs,etai,xpi,xni,xai,xhi,yehi,
     $        abari,stot,iident,iprocnum)
         xalphai=xai
         xheavyi=xhi
         xmuei=etai*tempi
         xmuhati=0.0
         u2i=stot
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
         ui=stot
      endif
c
c--store values
c
      if (xpi.le.1d-15) xpi=0.
      if (xni.le.1d-15) xni=0.
      pri=ptot
      vsoundi=min(cs,0.3333*dble(clight))
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
      t9l=3d-4
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
     1     ptot,cs,eta,yp,yn,xa,xh,yeh,abar,stot,
     1     iident,iprocnum)
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
     1            xa,xh,yeh,zbar,abar,ubind,dubind,iident,iprocnum)
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
            print *,'t9 too small'
            print *,'t9,t9old',t9,t9old
            print *,'u,ures,dut',u,ures,dut
            print *,'rhocgs,rhoold',rhocgs,rhoold
            print *,'yei,yeold',yei,yeold
            print *,'ediss,etot,ucoul',ediss,etot,ucoul
            t9=tempi*uotemp1
            t9old=t9
            dt9old=dabs(5.-t9l)
            dt9=dt9old
         endif
c
c--solve nse 
c
         call nsetemp(t9old,rhocgs,yei,t9,yp,yn,
     1                xa,xh,yeh,zbar,abar,ubind,dubind,iident,iprocnum)
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
      print *,'rootemp3: no convergence for part. rho(cgs),t9',
     1         rhoi,t9
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
     1     xa,xh,yeh,zbar,abar,ubind,dubind,
     1     iident,iprocnum)
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
          write(*,*) 'NSE mis-stored entering nserho',iprocnum,iident
          write(*,*) 'T9, rho, ye',t9,rhoold,yeold
          write(*,*) 'inconsistent with yp, yn',yp,yn
          ider=0
          call nsestart(t9,rhoold,yeold,yp,yn,iident,iprocnum)
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
      write(*,*)'Ye loop failure',iprocnum,iident
      write(*,*)t9,rhoold,rho
      write(*,*)yetmp,yelast,ye
      write(*,*)yeold,delye,yp,yn
      write(*,*)kit,zbar,abar
      write(*,*)kmax,ubind
      write(*,*)testk,testzy,testay,testyp,testyn
      call swerror
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
      write(*,*)'rho loop failure',iprocnum,iident
      write(*,*)t9,rhoold,rho,rhotmp
c     yetmp, yelast can SIGFPE in write
c      write(*,*)yetmp,yelast,ye
      write(*,*)yeold,delye,yp,yn
      write(*,*)kit,zbar,abar
      write(*,*)kmax,ubind
      write(*,*)testk,testzy,testay,testyp,testyn
      call swerror
   70 Continue
   90 format(A25,3(1pe12.4),I3)
c
c--solve for actual rho and Ye, calcuate binding energy, average A and Z 
c
      ider=2
      call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,
     &             xa,xh,yeh,zbar,abar)
      If(kit.ge.kmax) Then
          write(*,*) 'Final step failure in nserho',iprocnum,iident
          write(*,*) kit,rhotmp,rho
          call swerror
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
	  if (delt9.gt.1d-14) then
              goto 80
	  else
	      print *,'delt9 too small',iprocnum,iident
              call swerror
c	      stop
	  endif
      Endif
      dubind=(ubindd-ubind)/(t9d-t9)
c
      return
      end
c 
      subroutine nsetemp(t9old,rho,ye,t9,yp,yn,
     1     xa,xh,yeh,zbar,abar,ubind,dubind,
     1     iident,iprocnum)
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
          write(*,*) 'NSE mis-stored entering nsetemp',iprocnum,iident
          write(*,*) 'T9, rho, ye',t9old,rho,ye
          write(*,*) 'inconsistent with yp, yn',yp,yn
          call swerror
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
      write(*,*)'nsetemp(1) did not convege!',iprocnum,iident
      write(*,*)t9old,t9tmp,t9
      write(*,*)i,delt9,rho
      write(*,*)ye,yp,yn
c     write below will fail with FPE if ubind has underflowed
c      write(*,*)kit,kmax,ubind
      write(*,*)abar,zbar
      write(*,*)'before nsestart',yp,yn
      call nsestart(t9,rho,ye,yp,yn,iident,iprocnum)
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
	  if (delt9.gt.1d-14) then
              goto 80
	  else
	      print *,'delt9 too small',iprocnum,iident
              call swerror
c	      stop
	  endif
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
      subroutine neutrino(steps,rhoi,yei,xpi,xni,hi,
     $     xheavyk,xalphak,yehk,etai,
     $     tempk,abari,gshifti,ri,
     $     pmassi, vsoundi, xmuhati,
     $     ynuei,ynuebi,ynuxi,unuei,unuebi,unuxi,
     $     rmaxnue,rmaxnueb,rmaxnux,
     $     ftrape,ftrapb,ftrapx,jtrape,jtrapb,jtrapx,iident)
c
      double precision xheavyk,xalphak,yehk
      real ftrape,ftrapb,ftrapx,jtrape,jtrapb,jtrapx
      real rmaxnue,rmaxnueb,rmaxnux
      integer ebetaeqi, pbetaeqi
      common /beta/ ebetaeqi, pbetaeqi
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /nulums/ rlumnue, rlumnueb, rlumnux,
     $     enue, enueb, enux, e2nue, e2nueb, e2nux,
     $     enues, enuebs, enuxs, dee, deeb, dex
      common /nuout/ rmxnue,rmxnueb,rmxnux
c
c     tempi is changed in nsesolv
      tempi=tempk
      xheavyi=real(xheavyk)
      xalphai=real(xalphak)
      yehi=real(yehk)
      call nuinit(rhoi,yei,abari,xpi,xni,
     $     tempi,etai,xalphai,xheavyi,yehi,
     $     ynuei,ynuebi,ynuxi,unuei,unuebi,unuxi,iident)
c
      call nucheck(ri,hi,dnuei,dnuebi,dnuxi,
     $     ynuei,ynuebi,ynuxi,unuei,unuebi,unuxi,
     $     rmaxnue, rmaxnueb, rmaxnux,
     $     ftrape,ftrapb,ftrapx,jtrape,jtrapb,jtrapx)
c
      call nupress(ri,rhoi,unuei,unuebi,unuxi,
     $     prnui,rmaxnue,rmaxnueb,rmaxnux)
c
c--emission from electron/proton capture
c
      call nuecap(steps,tempi,rhoi,etai,xpi,xni,ri,gshifti,pmassi,yei,
     $     rmaxnue,rmaxnueb,rmaxnux,iident)
c
c--emission from pair/plasma effects
c
      call nupp(tempi,rhoi,etai,yei,ri,gshifti,pmassi,
     $     rmaxnue,rmaxnueb,rmaxnux)
c
      call nusphere(steps,hi,ri,ynuei,ynuebi,ynuxi,
     $     unuei,unuebi,unuxi,gshifti,pmassi,
     $     rmaxnue,rmaxnueb,rmaxnux)
c
      return
      end
      subroutine eosfl(rhoi,pri,
     $     ui,u2i,yei,tempi,ifleosi,abari,
     $     xpi,xni,xpfi,p2i,p3i,p4i,ufreezi,
     $     xmuei,xmuhi,etai,rhoprev,temprev,
     $     yeprev,xpprev,xnprev,iident,iprocnum)

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
c      rhoswe=0.55d3
c      rhoswe=1.05d3
c      rhoswe=2.05d3
c      t9nse=21.d0
c      tfreeze=t9nse-7.0d0
c  These are the standard numbers
      rhoswe=2.5d4
c remember to change this in eosgen as well
      t9nse=8.d0
      tfreeze=t9nse-1.5
c
c--entropy conversion factor
      sfac=avokb*utemp/uergg      
c
      rhocgs=rhoi*udens
      if(ifleosi.eq.1.and.(tempi.gt.t9nse.or.rhoi.gt.4.0*rhoswe)) then
c	 print *,tempi,rhocgs,yei,xpi,xni,t9nse
         call nsestart(tempi,rhocgs,yei,xpi,xni,iident,iprocnum)
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
         rhoprev=rhoi
         temprev=tempi
      endif
      if(ifleosi.eq.2.and.tempi.le.tfreeze) then
         ifleosi=1
c	 print *, tempi,rhocgs,yei,xpi,xni,tfreeze
         call nsetemp(tempi,rhocgs,yei,tempi,xpi,xni,
     1        xai,xhi,yehi,zbari,abari,ubind,dubind,iident,iprocnum)
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
         call nsestart(tempi,rhocgs,yei,xpi,xni,iident,iprocnum)
         call nsetemp(tempi,rhocgs,yei,tempi,xpi,xni,
     $        xai,xhi,yehi,zbari,abari,ubind,dubind,iident,iprocnum)
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
      subroutine nuinit(rhoi,yei,abari,xpi,xni,
     $     tempi,etai,xalphai,xheavyi,yehi,
     $     ynuei,ynuebi,ynuxi,unuei,unuebi,unuxi,iident)
c*************************************************************
c
c  This subroutine zeroes out whatever is necessary for the
c  neutrino physics and computes the MeV mean energies
c  Note that ynux is the abundance of any single specie of 
c  tau, mu, antineutrino or neutrino, but unux is the total
c  energy in all 4 fields
c  We also compute opacities, diffusion coefficients and
c  degeneracy for each species
c
c*************************************************************
c
      parameter(rcrit=1.0)
c
c--so that we don't have to go double precision:
c--avo*1e-44=6.02e-20
      parameter(avosig=6.02e-21)
c--same thing for 1e13 (fm/cm) avo**-1/3
      parameter(fachn=1.2e5)
c--hbar*c/kb in fermi*Kelvin: why a cap for Kelvin and not for fermi?
      parameter(facdeb=2.3e12)
c--e**2/kb/(1fm)
      parameter(gamfac=1.67e10)
c--nucleon elastic xsection from Bowers and Wilson 1982, ApJSup 50
      parameter(sigpel=1.79)
      parameter(signel=1.64)
c--Mayle's thesis for alpha/nu elastic scattering
      parameter(sigalel=0.048)
c--prefactor to BBAL formula (Bethe 1990) for coherent scattering
      parameter(sigheel=1.7)
c--neutrino absorptions by free nucleons:
      parameter(sigabs=9.0)
c--xneutrino-el scatterings: neutral currents
      parameter(sigelx=0.14)
c--eneutrino-el scatterings: charged currents + neutral currents
      parameter(sigele=0.92)
      parameter(sigeleb=0.39)
c--1/(avo*pi**2*(hbar*c)**3 in Mev-3 cm-3 nucleon g-1)
      parameter(prefac=2.19e7)
c--factors for S3(eta)=F3(eta)+F3(-eta)
      parameter(pi=3.141592)
      parameter(s3a=7.*pi*pi*pi*pi/60.)
      parameter(s3b=0.5*pi*pi)
      parameter(pi43=4.*pi/3.)
c
      double precision umass
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /nulums/ rlumnue, rlumnueb, rlumnux,
     $     enue, enueb, enux, e2nue, e2nueb, e2nux,
     $     enues, enuebs, enuxs, dee, deeb, dex
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
c
      fac=avosig*udens*udist
c
c--zero-out dunu, dye, dynus, denus
c
      uconv=1./umevnuc
      yfac=prefac/udens
c
      dunui=0.0
      dyei=0.0
      dynuei=0.0
      dynuebi=0.0
      dynuxi=0.0
      dunuei=0.0
      dunuebi=0.0
      dunuxi=0.0

c

      if (ynuei.gt.0.) then
         enueti=uconv*unuei/ynuei
      else
         enueti=3.0
      endif
c
      if (ynuebi.gt.0.) then
         enuebti=uconv*unuebi/ynuebi
      else
         enuebti=3.0
      endif
      if (ynuxi.gt.0.) then
         enuxti=uconv*unuxi/(4.*ynuxi)
      else
         enuxti=3.0
      endif
c
      facrho=fac*rhoi
c--Debye radius
      if (etai*tempi.ge.15.) then
         rd=10.0*facdeb/(etai*utemp*tempi)
      else
         rd=89.0*sqrt(utemp*tempi/(rhoi*udens))
      endif
      if (xheavyi.ge.1e-3) then
c--mean distance between heavy nuclei (fm)
         rhn=fachn/(pi43*xheavyi*udens*rhoi/abari)**0.3333333
c--ratio of Coulomb inter-nuclei energy to thermal energy:
         rgamma=gamfac*(abari*yehi)**2*exp(-rhn/rd)/
     $        (rhn*utemp*tempi)
      else
c--a big number...
         rhn=1.e9
         rgamma=0.0
      endif
      cohfac=abari*(1.-1.08*yehi)/6.
c--neutral current/heavies scattering without corrections
      alphanh=facrho*xheavyi*sigheel*cohfac
c--neutral current/nucleon and alpha  scattering
      alphann=facrho*(xpi*sigpel+xni*signel+
     $     xalphai*sigalel)
c--charged current/neutron scattering
      alphacn=facrho*xni*sigabs
c--charged current/proton scattering
      alphacp=facrho*xpi*sigabs
c--xnu neutral current/e+e- scattering
      alphaxe=facrho*sigelx
c--enu charged+neutral current/e+e- scattering
      alphaee=facrho*sigele
      alphaebe=facrho*sigeleb
      tmev=utmev*tempi
      tmev2=tmev*tmev
c--compute local electron energies
      if (etai*tmev.ge.1.) then
         call integrals(etai,f0,f1,f2,f3,f4,f5,0,df2,df3)
         elmean=tmev*f3/f2
         elmean2=tmev2*f4/f2
         etai2=etai*etai
         etai4=etai2*etai2
         yemean=yfac*tmev2*tmev2/rhoi*(s3a+0.25*etai4+s3b*etai2)
      else
         elmean=3.*tmev
         elmean2=9.*tmev2
         yemean=yei*elmean
      endif
      if (enueti.gt.5.) then
c--figure out neutrino degeneracy 
         call rooteta(tmev,rhoi,ynuei,unuei,etanuei,tempnuei,iident)
         ediff=enueti
         ediff2=enueti*enueti
      else
         etanuei=0.0
         tempnuei=tmev
         ediff=amax1(elmean,enue)
         ediff2=amax1(elmean2,e2nue)
      endif
c
c--de Broglie wavelength fm
      debro=1240./ediff
      rrd=rd/debro
      shield=(rrd+.5/rrd)/(rrd+1./rrd+3.)
      rrhn=rhn/debro
      struct=(rrhn**2+(1./rrhn)**3)/(rrhn**2+rgamma)
      opacnh=struct*shield*alphanh*ediff2
      opac=(alphann+alphacn)*ediff2+alphaee*ediff*yemean
      dnuei=1./(opac+opacnh)
      if (enuebti.gt.5.) then
         call rooteta(tmev,rhoi,ynuebi,unuebi,etanuebi,tempnuebi,0)
         ediff=enuebti
         ediff2=enuebti*enuebti
      else
         etanuebi=0.0
         tempnuebi=tmev
         ediff=amax1(elmean,enueb)
         ediff2=amax1(elmean2,e2nueb)
      endif
c--de Broglie wavelength fm
      debro=1240./ediff
      rrd=rd/debro
      shield=(rrd+0.5/rrd)/(rrd+1./rrd+3.)
      rrhn=rhn/debro
      struct=(rrhn**2+(1./rrhn)**3)/(rrhn**2+rgamma)
      opacnh=struct*shield*alphanh*ediff2
      opac=(alphann+alphacp)*ediff2+alphaebe*ediff*yemean
      dnuebi=1./(opac+opacnh)
      if (enuxti.gt.5.) then
         call rooteta(tmev,rhoi,ynuxi,.25*unuxi,etanuxi,tempnuxi,0)
         ediff=enuxti
         ediff2=enuxti*enuxti
      else
         etanuxi=0.0
         tempnuxi=tmev
         ediff=amax1(elmean,enux)
         ediff2=amax1(elmean2,e2nux)
      endif
c--de Broglie wavelength fm
      debro=1240./ediff
      rrhn=rhn/debro
      struct=(rrhn**2+(1./rrhn)**3)/(rrhn**2+rgamma)
      opacnh=struct*alphanh*ediff2
      opac=alphann*ediff2+alphaxe*ediff*yemean
      dnuxi=1./(opac+opacnh)
c
c--set neutrino luminosities to zero.
c
      rlumnue=0.
      rlumnueb=0.
      rlumnux=0.
c
c--initialize energy sums
c
      enue=0.
      enueb=0.
      enux=0.
      e2nue=0.
      e2nueb=0.
      e2nux=0.
c
      return
      end
c
      subroutine nucheck(ri,hi,dnuei,dnuebi,dnuxi,
     $     ynuei,ynuebi,ynuxi,unuei,unuebi,unuxi,
     $     rmaxnue, rmaxnueb, rmaxnux, 
     $     ftrape,ftrapb,ftrapx,jtrape,jtrapb,jtrapx)
c*************************************************************
c
c  
c  This routine determines trapping status by computing opacities 
c  and diffusion coefficients. Also mops up residual neutrinos
c  from untrapped particles if they are present in tiny 
c  quantities
c
c*************************************************************
c
      parameter(rcrit=1.0)
      parameter(tiny=1e-8)
c
      common /nuout/ rmxnue,rmxnueb,rmxnux
      real ftrape,ftrapb,ftrapx,jtrape,jtrapb,jtrapx
      double precision umass
c
c-- check for nu trapping
c
      rmxnue=0.0
      rmxnueb=0.0
      rmxnux=0.0
      che=ftrape*2.*hi
      chb=ftrapb*2.*hi
      chx=ftrapx*2.*hi
      rnue=che/dnuei
      rnueb=chb/dnuebi
      rnux=chx/dnuxi
c--note if nux trapped, then nue trapped and contraposition
      if (rnue.gt.rcrit) rmxnue=ri
      if (rnueb.gt.rcrit) rmxnueb=ri
      if (rnux.gt.rcrit) rmxnux=ri
c
c--mop up neutrino garbage
c
      if (ynuei.lt.tiny) then
         ynuei=0.
         unuei=0.
      endif
      if (ynuebi.lt.tiny) then
         ynuebi=0.
         unuebi=0.
      endif
      if (ynuxi.lt.tiny) then
         ynuxi=0.
         unuxi=0.
      endif
c
      return
      end
c

      subroutine nupress(ri,rhoi,unuei,unuebi,unuxi,
     $     prnui,rmaxnue,rmaxnueb,rmaxnux)
c*******************************************************
c
c This subroutine coputes the energy trapped in the form
c of neutrinos and derives a pressure from that
c
c******************************************************
c
c--gamma=4/3 since nus are always relativistic
      rhoi3=0.3333333333333*rhoi
      prnui=0.0
      if (ri.lt.rmaxnue) prnui=prnui+rhoi3*unuei
      if (ri.lt.rmaxnueb) prnui=prnui+rhoi3*unuebi
      if (ri.lt.rmaxnux) prnui=prnui+rhoi3*unuxi
c
      return
      end
c
      subroutine nuecap(steps,tempi,rhoi,etai,xpi,xni,ri,
     $     gshifti,pmassi,yei,rmaxnue,rmaxnueb,rmaxnux,iident)
c****************************************************
c
c this subroutine computes the neutrino production
c by e+/e- capture on nucleons
c Note: all neutrino energies are in MeV          
c
c****************************************************
c
c--1/(avo*pi**2*(hbar*c)**3 in Mev-3 cm-3 nucleon g-1)
      parameter(prefac=2.19e7)
c-- tffac=(6pi^2/2)^2/3 hbar*2/(2 mp kb)*avo^2/3
      parameter (tfermi=164.)
c
      double precision umass
c
      integer ebetaeqi, pbetaeqi
      common /beta/ ebetaeqi, pbetaeqi
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /nulums/ rlumnue, rlumnueb, rlumnux,
     $     enue, enueb, enux, e2nue, e2nueb, e2nux,
     $     enues, enuebs, enuxs, dee, deeb, dex
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
c
      yfac=prefac/udens
      tf=tfermi/utemp
c
      dt9=9.*steps
      if(tempi.ge.6..or.etai*tempi.ge.10.)then
         xfn=xpi+xni
c
c--compute electron and positron capture
c  -------------------------------------
c
         if(xfn.ne.0)then
            call epcapture(erate,prate,due,dup,tempi,etai,2)
            call integrals(etai,f0,f1,f2,f3,f4,f5,0,df2,df3)
            tmev=tempi*utmev
c
c--supress rates by end state degeneracy
c
            expon=exp(amin1(f3/f2*tmev/tempnuei-etanuei,
     $           50.))
c--bnue: neutrino end-state blocking
            bnue=expon/(1.+expon)
c
            expon=exp(amin1(3.*tmev/tempnuebi-etanuebi,50.))
c--   bnueb: anti-neutrino end-state blocking
            bnueb=expon/(1.+expon)
c
            tfne=tf*(xni*rhoi*udens)**(0.666666666)
            tfpr=tf*(xpi*rhoi*udens)**(0.666666666)
            tprot=amax1(tempi,tfpr)
            tneut=amax1(tempi,tfne)
            expon=exp(amin1((tprot-tfne)/tempi, 50.))
c--bne: neutron end-state blocking
            bne=expon/(1.+expon)
            expon=exp(amin1((tneut-tfpr)/tempi, 50.))
c--bpr: proton end-state blocking
            bpr=expon/(1.+expon)
            blocke=bnue*bne
            blockp=bnueb*bpr
            erate=erate*blocke
            prate=prate*blockp
            due=due*blocke
            dup=dup*blockp
            if (prate.lt.0.) then
               prate=1e-10
               dup=0.
            endif
            facp=xpi*erate
            facn=xni*prate
c--   check for beta eq.
            tmev2=tmev*tmev
            tmev3=tmev2*tmev
c
c-- if rate*(9 time steps) is larger than the e- abundance
c-- beta eq. is declared, nuecap zeroed-out, 
c-- will be taken care of in nubeta
c
            if (facp*dt9.gt.yei.and.ri.lt.rmaxnue) then
               ebetaeqi=1
               facp=0.0
               due=0.0
            else
               ebetaeqi=0
            endif
            call integrals(-etai,f0,f1,f2,f3,f4,f5,0,df2,df3)
            ypo=yfac*tmev3*f2/rhoi
            if (facn*dt9.gt.yei.and.ri.lt.rmaxnueb) then
               pbetaeqi=1
               facn=0.0
               dup=0.0
            else
               pbetaeqi=0
            endif
c
c--change in Ye
c
            dyei=dyei+(facn-facp)
c
c--energy loss and mean neutrino energies in MeV
c
            enumean=due/(erate+1e-30)
            enubmean=dup/(prate+1e-30)
            due=due*umevnuc*xpi
            dup=dup*umevnuc*xni
            dunui=dunui+due+dup
            if (ri.lt.rmaxnue) then
               dynuei=dynuei+facp
               dunuei=dunuei+due
            else
c--weigh the mean energy sums by luminosity
               rlnue=due*pmassi*gshifti
               rlumnue=rlumnue+rlnue
               enue=enue+enumean*gshifti*rlnue
               e2nue=e2nue+enumean*enumean*gshifti*gshifti*rlnue
            end if
            if (ri.lt.rmaxnueb) then
               dynuebi=dynuebi+facn
               dunuebi=dunuebi+dup
            else
c--weigh the mean energy sums by luminosity
               rlnueb=dup*pmassi*gshifti
               rlumnueb=rlumnueb+rlnueb
               enueb=enueb+enubmean*gshifti*rlnueb
               e2nueb=e2nueb+enubmean*enubmean*gshifti*gshifti*rlnueb
            end if
         endif
      endif
c
      return 
      end
c
      subroutine nupp(tempi,rhoi,etai,yei,ri,gshifti,pmassi,
     $     rmaxnue,rmaxnueb,rmaxnux)
c****************************************************
c
c this subroutine computes the neutrino production
c by e+/e- capture on nucleons
c Note: all neutrino energies are in MeV          
c
c****************************************************
c
      double precision umass,ugserg
      double precision rhocgs, tempk, yek, deta
      double precision dunuel,dunuxl,enuel,enuxl,enuel2,enuxl2
c
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /nulums/ rlumnue, rlumnueb, rlumnux,
     $     enue, enueb, enux, e2nue, e2nueb, e2nux,
     $     enues, enuebs, enuxs, dee, deeb, dex
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
c
      ugserg=dble(utime/uergg)
      uconv=1./umevnuc
c
c--compute pair plasma processes
c  -----------------------------
      if(tempi.ge.6.)then
         rhocgs=rhoi*udens
         tempk=utemp*tempi
         deta=etai
         yek=yei
         call pppb(rhocgs,tempk,yek,deta,dunuel,dunuxl,
     $        enuel,enuxl,enuel2,enuxl2)
c
c--supress rates by end state degeneracy
c
         tmev=tempi*utmev
         expon=exp(amin1(real(enuel)/tempnuei-etanuei,50.))
c--bnue: e-neutrino end-state blocking
         bnue=expon/(1.+expon)
         expon=exp(amin1(real(enuel)/tempnuebi-etanuebi,50.))
c--bnueb: anti e-neutrino end-state blocking
         bnueb=expon/(1.+expon)
         expon=exp(amin1(real(enuxl)/tempnuxi-etanuxi,50.))
c--bnux: x neutrino end-state blocking
         bnux=expon/(1.+expon)
         blocke=bnue*bnueb
         blockx=bnux*bnux
         dunuel=dunuel*ugserg*blocke
         dunuxl=dunuxl*ugserg*blockx
         dunui=dunui+dunuxl+dunuel
         if (ri.lt.rmaxnue) then
            facnue=0.5*uconv*dunuel/enuel
            dynuei=dynuei+facnue
            dunuei=dunuei+0.5*dunuel
         else
            rlnue=0.5*dunuel*pmassi*gshifti
            rlumnue=rlumnue+rlnue
            enue=enue+enuel*gshifti*rlnue
            e2nue=e2nue+enuel2*gshifti*gshifti*rlnue
         endif
         if (ri.lt.rmaxnueb) then
            facnue=0.5*uconv*dunuel/enuel
            dynuebi=dynuebi+facnue
            dunuebi=dunuebi+0.5*dunuel
         else
            rlnueb=0.5*dunuel*pmassi*gshifti
            rlumnueb=rlumnueb+rlnueb
            enueb=enueb+enuel*gshifti*rlnueb
            e2nueb=e2nueb+enuel2*gshifti*gshifti*rlnueb
         endif
         if (ri.lt.rmaxnux) then
            facnux=uconv*dunuxl/enuxl
c-- 0.25 to divide the production among 4 nu species
            dynuxi=dynuxi+0.25*facnux
            dunuxi=dunuxi+dunuxl
         else
            rlnux=dunuxl*pmassi*gshifti
            rlumnux=rlumnux+rlnux
            enux=enux+enuxl*gshifti*rlnux
            e2nux=e2nux+enuxl2*gshifti*gshifti*rlnux
         endif
      endif
c
      return 
      end
c
      subroutine nusphere(steps,hi,ri,ynuei,ynuebi,ynuxi,
     $     unuei,unuebi,unuxi,gshifti,
     $     pmassi,rmaxnue,rmaxnueb,rmaxnux)
c************************************************************
c
c  This subroutine takes care of neutrino depletion from
c  particles which are optically thin and have non-zero
c  ynus either because:
c     1) They have gone from optically think to optically thin
c  or
c     2) They have optically thick neighbours diffusing neutrinos
c        into them
c
c  The losses are setup so that the e-folding length is
c  3 time steps or the diffusion time accross the particle,
c  whichever is larger.
c
c************************************************************
c
      parameter(tiny=1e-8)
c
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /nulums/ rlumnue, rlumnueb, rlumnux,
     $     enue, enueb, enux, e2nue, e2nueb, e2nux,
     $     enues, enuebs, enuxs, dee, deeb, dex
      common /konst/ gg, clight, arad, bigr, xsecnn, xsecne
c
      dee=0.0
      deeb=0.0
      dex=0.0
c-- minimum time scale of emission =3 time steps
      t3step=3.*steps
c
      enues=0.
      enuebs=0.
      enuxs=0.
      e2nues=0.
      e2nuebs=0.
      e2nuxs=0.
c     cthird=clight/3.
c    3d only!!!
      cthird=clight
      twave=4.*hi/clight
      if (ri.gt.rmaxnue.and.ynuei.ge.tiny) then
         tdif=(4.*hi)**2/(cthird*dnuei)
         ttran=amax1(tdif,twave)
         if (ttran.lt.t3step) then
            ttran=t3step
         endif
         fac=1./ttran
         facye=fac*ynuei
         dynuei=dynuei-facye
         facue=fac*unuei
         dunuei=dunuei-facue
         facee=facue*pmassi*gshifti
         rlumnue=rlumnue+facee
         dee=dee+facee
         enues=enues+facee*enueti*gshifti
         e2nues=e2nues+facee*enueti*enueti*gshifti*gshifti
      endif
      if (ri.gt.rmaxnueb.and.ynuebi.ge.tiny) then
         tdif=(4.*hi)**2/(cthird*dnuebi)
         ttran=amax1(tdif,twave)
         if (ttran.lt.t3step) ttran=t3step
         fac=1./ttran
         facyeb=fac*ynuebi
         dynuebi=dynuebi-facyeb
         facueb=fac*unuebi
         dunuebi=dunuebi-facueb
         faceeb=facueb*pmassi*gshifti
         rlumnueb=rlumnueb+faceeb
         deeb=deeb+faceeb
         enuebs=enuebs+faceeb*enuebti*gshifti
         e2nuebs=e2nuebs+faceeb*enuebti*enuebti*gshifti*gshifti
      endif
      if (ri.gt.rmaxnux.and.ynuxi.ge.tiny) then
         tdif=(4.*hi)**2/(cthird*dnuxi)
         ttran=amax1(tdif,twave)
         if (ttran.lt.t3step) ttran=t3step
         fac=1./ttran
         facyx=fac*ynuxi
         dynuxi=dynuxi-facyx
         facux=fac*unuxi
         dunuxi=dunuxi-facux
         facex=facux*pmassi*gshifti
         rlumnux=rlumnux+facex
         dex=dex+facex
         enuxs=enuxs+facex*enuxti*gshifti
         e2nuxs=e2nuxs+facex*enuxti*enuxti*gshifti*gshifti
      endif

      enue=enue+enues
      e2nue=e2nue+e2nues
      enueb=enueb+enuebs
      e2nueb=e2nueb+e2nuebs
      enux=enux+enuxs
      e2nux=e2nux+e2nuxs
c
      return
      end
      subroutine rooteta(tmev,rho,ynu,unu,eta,temp,iident)
c***********************************************************
c
c  subroutine computes the degree of degeneracy in a nu-field
c  assuming a thermalized Fermi distribution, i.e.:
c     n = ynu*rho ~ F2(eta)
c     E = unu*rho ~ F3(eta)
c
c***********************************************************
c
c const= 2*pi**2*hbar**3*c**3
c--2.*avo*pi**2*(hbar*c)**3 in Mev+3 cm+3 nucleon g-1)
      parameter(prefac=9.2e-8)
      parameter(etamin=1e-2)
      parameter(etamax=50.0)
      parameter(tol=1e-3)
c
      double precision umass
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
c
c      if(iident.eq.6369)then
c         write(*,*) 'tmev rho ynu unu eta',tmev,rho,ynu,unu,eta
c      endif
      etain=eta
      yfac=prefac*udens
      ufac=yfac/umevnuc
      a=(rho*ynu*yfac)**(0.333333333333)
      b=sqrt(sqrt(rho*unu*ufac))
      if (eta.lt.etamin) eta=etamin
c
      do i=1,80
         call integrals(eta,f0,f1,f2,f3,f4,f5,1,df2,df3)
         f3_14=sqrt(sqrt(f3))
         f2_13=f2**(0.33333333333)
         t3=a*f3_14
         t2=b*f2_13
         f=a*f3_14 - b*f2_13
         df=a*0.25*f3_14/f3*df3 - b*0.33333333*f2_13/f2*df2
         deta=f/(df+1e-30)
         eta=eta-deta
c--"effective" neutrino temperature in MeVs
         temp=0.5*(a/f2_13+b/f3_14)
c
c         if(iident.eq.6369)then
c            write(*,*) 'eta deta temp f df',eta,deta,temp,f,df
c            call swflush
c         endif
         if (eta.ge.etamax) then
            eta=etain
            call integrals(eta,f0,f1,f2,f3,f4,f5,1,df2,df3)
            f3_14=sqrt(sqrt(f3))
            f2_13=f2**(0.33333333333)
            temp=0.5*(a/f2_13+b/f3_14)
            return
         elseif (eta.lt.etamin.and.etain.gt.10.) then
c           print *,'eta convergence problem'
            eta=etain
            call integrals(eta,f0,f1,f2,f3,f4,f5,1,df2,df3)
            f3_14=sqrt(sqrt(f3))
            f2_13=f2**(0.33333333333)
            temp=0.5*(a/f2_13+b/f3_14)
            return
         elseif (eta.lt.etamin) then
            eta=0.
            temp=tmev
            return
         elseif (abs(deta/eta).lt.tol) then
            return
         endif
      enddo
c
      eta=etain
      call integrals(eta,f0,f1,f2,f3,f4,f5,1,df2,df3)
      f3_14=sqrt(sqrt(f3))
      f2_13=f2**(0.33333333333)
      temp=0.5*(a/f2_13+b/f3_14)
      write(*,*) 'eta convergence failed!'
      write(*,*) 'id tmev rho ynu unu eta',iident,tmev,rho,ynu,unu,eta
c
      end
c         
      subroutine pppb(rho,temp,ye,eta,dunuel,dunuxl,
     1                enuel,enuxl,enuel2,enuxl2)
c******************************************************
c
c Compute energy loss (ergs/g/s) via pair, photo, plasma and
c bremsstahlung processes with formulae from:
c Itoh et al. (1989), ApJ 339, p.354
c Average energies (MeV) eyeballed from:
c Schinder et al. (1987) , ApJ 313, p.531
c rho and temp should be in cgs
c Note: I only implemented plasma and pair processes,
c       but anyone should feel free to add the rest in...
c Note: These processes do not change Ye, neutrinos
c       of all families are created in pairs
c       denuel=denue+denueb=2*denue
c       denuxl=denumu+denumub+denutau+denutaub=4*denumu
c
c******************************************************
      implicit double precision(a-h,o-z)
c
      parameter(sinw2=0.23d0)
      parameter(cv=0.5d0+2.d0*sinw2)
      parameter(cv2=cv*cv)
      parameter(ca=0.5d0)
      parameter(ca2=ca*ca)
      parameter(cvp=1.d0-cv)
      parameter(cvp2=cvp*cvp)
      parameter(cap=1.d0-ca)
      parameter(cap2=cap*cap)
      parameter(dn=2.d0)
      parameter(facq=((cv2-ca2)+dn*(cvp2-cap2))/
     1               ((cv2+ca2)+dn*(cvp2+cap2)))
      parameter(emassk1=1.d0/5.9302d9)
c
      rhoye=rho*ye
      tme=temp*emassk1
      tme2=tme*tme
      tme3=tme2*tme
      tme4=tme3*tme
      tme6=tme4*tme2
      tme8=tme6*tme2
      tmep5=dsqrt(tme)
      tmem1=1.d0/tme
      tmem2=tmem1*tmem1
      tmem3=tmem2*tmem1
      xsi=(rhoye*(1.d-9))**0.3333333d0*tmem1
      xsi2=xsi*xsi
      xsi3=xsi2*xsi
c--temp in Mev
      tmev=temp*8.62d-11
c--electron chemical potential in MeV
      emumev=eta*tmev
c
c-- pairs 
c
      if (temp.ge.1d10) then
         fpair=dexp(-4.9924d0*xsi)*
     1         (6.002d19 + 2.084d20*xsi + 1.872d21*xsi2)/
     2         (xsi3 + 1.238d0*tmem1 - 0.8141d0*tmem2)
      else
         fpair=dexp(-4.9924d0*xsi)*
     1         (6.002d19 + 2.084d20*xsi + 1.872d21*xsi2)/
     2         (xsi3 + 9.383d-1*tmem1 - 4.141d-1*tmem2 + 5.829d-2*tmem3)
      endif
      gpair=1.d0 - 13.04d0*tme2 + 133.5d0*tme4 + 1534.d0*tme6 +
     1      918.6d0*tme8
      qpair=((1.d0 +
     1        rhoye/(7.692d7*tme3 + 9.715d6*tmep5))**(-0.3d0))/
     2      (10.748d0*tme2 + 0.3967d0*tmep5 + 1.005d0)
      pairfac=0.5d0*(1.d0+facq*qpair)*gpair*fpair*dexp(-2.d0*tmem1)
      deepair=(cv2+ca2)*pairfac
      dexpair=dn*(cvp2+cap2)*pairfac
      epair=dmax1(4.d0*tmev,.75d0*emumev)
c--I don't know where 4 comes from (process heavily weighted towards
c--high energies). 0.75 because <Kinetic E of el>=3/4efermi
      epair2=epair*epair
c
c-- plasma 
c
      fplas=dexp(-0.56457d0*xsi)*
     1      (2.32d-7 + 8.449d-8*xsi + 1.787d-8*xsi2)/
     2      (xsi3 + 2.581d-2*tmem1 + 1.734d-2*tmem2 + 6.99d-4*tmem3)
      plasfac=rhoye*rhoye*rhoye*fplas
      deeplas=cv2*plasfac
      dexplas=dn*cvp2*plasfac
c--why 0.05d0, don't know have to trust Schindler
c--tmev because that's the energy of the photons
      eplas=dmax1(0.5d0*tmev,0.05d0*emumev)
      eplas2=eplas*eplas
c
c--total emission ergs/cc/s
c
      denuel=deepair+deeplas
      denuxl=dexpair+dexplas
c
c--emission ergs/g/s
c
      dunuel=denuel/rho
      dunuxl=denuxl/rho
c
c--luminosity weighted energies:
c
      enuel=(deepair*epair+deeplas*eplas)/denuel
      enuel2=(deepair*epair2+deeplas*eplas2)/denuel
      enuxl=(dexpair*epair+dexplas*eplas)/denuxl
      enuxl2=(dexpair*epair2+dexplas*eplas2)/denuxl
c
      return
      end
c

      subroutine integrals(eta,f0,f1,f2,f3,f4,f5,ider,df2,df3)
c************************************************************
c                                                           *
c  This subroutine calculates numerical approximations to   *
c  the fermi integrals (Takahashi, El Eid, Hillebrandt, 1978*
c  AA, 67, 185.                                             *
c                                                           *
c************************************************************
c
c--case where eta > 1e-3
c
      if(eta.ge.1.e-3) then
         eta2=eta*eta
         eta3=eta*eta2
         eta4=eta*eta3
         eta5=eta*eta4
         eta6=eta*eta5
         if (eta.le.30.) then
            f1=(0.5*eta2+1.6449)/
     1         (1.+exp(-1.6855*eta))
            f2num=eta3/3.+3.2899*eta
            f2den=1.-exp(-1.8246*eta)
            f2=f2num/f2den
            f3num=0.25*eta4+4.9348*eta2+11.3644
            f3den=1.+exp(-1.9039*eta)
            f3=f3num/f3den
            f4=(0.2*eta5+6.5797*eta3+45.4576*eta)/
     1         (1.-exp(-1.9484*eta))
            f5=(eta6/6.+8.2247*eta4+113.6439*eta2+236.5323)/
     1         (1.+exp(-1.9727*eta))
            f0=eta+alog(1.+exp(-eta))
            if (ider.eq.1) then
               df2=((eta2+3.2899)*f2den -
     1              f2num*1.8246*exp(-1.8246*eta))/
     1             (f2den*f2den)
               df3=((eta3+9.8696*eta)*f3den +
     1              f3num*1.9039*exp(-1.9039*eta))/
     2             (f3den*f3den)
            endif
         else
            f1=0.5*eta2+1.6449
            f2=eta3/3.+3.2899*eta
            f3=0.25*eta4+4.9348*eta2+11.3644
            f4=0.2*eta5+6.5797*eta3+45.4576*eta
            f5=eta6/6.+8.2247*eta4+113.6439*eta2+236.5323
            f0=eta
            if (ider.eq.1) then
               df2=eta2+3.2899
               df3=eta3+9.8696*eta
            endif
         endif
      elseif (eta.lt.-40.) then
         f0=eta
         f1=0.0
         f2=0.0
         f3=0.0
         f4=0.0
         f5=0.0
      else
         expeta=exp(eta)
         f1=expeta/(1.+0.2159*exp(0.8857*eta))
         f2=2.*expeta/(1.+0.1092*exp(0.8908*eta))
         f3=6.*expeta/(1.+0.0559*exp(0.9069*eta))
         f4=24.*expeta/(1.+0.0287*exp(0.9257*eta))
         f5=120.*expeta/(1.+0.0147*exp(0.9431*eta))
         f0=eta+alog(1.+exp(-eta))
      end if
c
      return
      end
c
      subroutine epcapture(erate,prate,due,dup,tempi,eta,iflag)
c************************************************************
c                                                           *
c  subroutine to compute the energy loss due to electron    *
c  capture on protons. Formalism taken from Takahashi, El   *
c  Eid, Hillebrandt, 1978, AA, 67, 185. (The energy release *
c  is in Mev/utime/nucleon and the rate in /utime/nucleon   *
c                                                           *
c************************************************************
c
      parameter (delta=1.531)
      parameter (q1=2.-2.*delta)
      parameter (q2=(1.-8.*delta+2.*delta**2)/2.)
      parameter (q3=2.*delta**2-delta)
      parameter (q4=(4.*delta**2-1.)/8.)
      parameter (p1=(2.-3.*delta))
      parameter (p2=(1.-12.*delta+6.*delta**2)/2.)
      parameter (p3=(-3.*delta+12.*delta**2-2.*delta**3)/2.)
      parameter (p4=(-1.+12.*delta**2-16.*delta**3)/8.)
      parameter (p5=(2.+3*delta-4.*delta**3)/8.)
      common /epcap/ betafac, c2cu, c3cu
c
      beta=betafac/tempi
      beta2=beta*beta
      beta3=beta*beta2
      beta4=beta*beta3
      beta5=beta*beta4
      beta6=beta*beta5
c
c--compute fermi integrals for electrons
c
      call integrals(eta,fe0,fe1,fe2,fe3,fe4,fe5,0,df2,df3)
c
c--compute fermi integrals for positrons
c
      call integrals(-eta,fp0,fp1,fp2,fp3,fp4,fp5,0,df2,df3)
c
c--e- + p -> n + nu rate / nucleons
c
      erate=fe4 + q1*fe3*beta + q2*fe2*beta2 + q3*fe1*beta3 +
     1      q4*fe0*beta4
      erate=c2cu*erate/beta5
c 
c--e+ + n -> p + nubar rate / nucleon
c
      prate=fp4 - q1*fp3*beta + q2*fp2*beta2 - q3*fp1*beta3 +
     1      q4*fp0*beta4
      prate=c2cu*prate/beta5
c
c--if nu trapped, bail out without computing energy loss
      if(iflag.eq.1) then
         due=0.
         dup=0.
         return
      endif
c
c--compute cooling due to electron capture (per nucleon)
c
      due=fe5 + p1*fe4*beta + p2*fe3*beta2 + p3*fe2*beta3 +
     1    p4*fe1*beta4 + p5*fe0*beta5
      due=c3cu*due/beta6
c
c--compute cooling due to positron capture (per nucleon)
c
      dup=fp5 - p1*fp4*beta + p2*fp3*beta2 - p3*fp2*beta3 +
     1    p4*fp1*beta4 - p5*fp0*beta5
      dup=c3cu*dup/beta6
c
      return
      end
c
      subroutine neutrino2(steps,rhoi,xpi,xni,etai,tempi,
     $     ri,pmassi,vsoundi,xmuhati,ynuei,ynuebi,ynuxi,
     $     unuei,unuebi,rlumnue,rlumnueb,rlumnux,enue,enueb,enux,
     $     e2nue,e2nueb,e2nux,gshifti,
     $     rmaxnue, rmaxnueb, rmaxnux, hi)
c
      integer ebetaeqi, pbetaeqi
      real rmaxnue,rmaxnueb,rmaxnux
      real hi
      common /beta/ ebetaeqi, pbetaeqi
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /nutrap/ dtrapnuei, dtrapnuebi
c
      call nuabs(steps,ri,etai,tempi,rhoi,xni,xpi,pmassi,
     $     ynuei,ynuebi,ynuxi,gshifti,rlumnue,rlumnueb,
     $     rlumnux,enue,enueb,enux,e2nue,e2nueb,e2nux,
     $     rmaxnue, rmaxnueb, rmaxnux)
c
      call nubeta(steps,vsoundi,hi,tempi,xmuhati,
     $     etai,rhoi,ynuei,ynuebi,unuei,unuebi)
c
      call nuann(hi,rhoi,tempi,etai,ri,
     $     ynuei,ynuebi,ynuxi,rmaxnue,rmaxnueb,rmaxnux)

      call nuconv(ri,hi,rhoi,tempi,gshifti,pmassi,
     $     ynuei,ynuebi,ynuxi,rmaxnue,rmaxnueb,rmaxnux,
     $     enue,enueb,enux,e2nue,e2nueb,e2nux,rlumnue,rlumnueb,
     $     rlumnux)
c
      return
      end
      subroutine nuabs(steps,ri,etai,tempi,rhoi,xni,xpi,pmassi,
     $     ynuei,ynuebi,ynuxi,gshifti,rlumnue,rlumnueb,
     $     rlumnux,enue,enueb,enux,e2nue,e2nueb,e2nux,
     $     rmaxnue,rmaxnueb,rmaxnux)
c****************************************************
c
c this subroutine computes the neutrino absorption
c by nucleons.
c Note: all neutrino energies are in MeV          
c
c****************************************************
      parameter (delta=0.783)
      parameter (deltab=1.805)
c-- tffac=(6pi^2/2)^2/3 hbar*2/(2 mp kb)*avo^2/3
      parameter (tfermi=164.)
c
      double precision umass
      integer ebetaeqi, pbetaeqi
      real rmaxnue,rmaxnueb,rmaxnux
      common /beta/ ebetaeqi, pbetaeqi
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /konst/ gg, clight, arad, bigr, xsecnn, xsecne
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /nutrap/ dtrapnuei, dtrapnuebi
c
c--the energy input is modified by delta=mn-mp-me=0.783 Mev
c
      if (enue.gt.0.) then
         facnue=(enue+delta)/enue
         xnorme=uergg/(facnue*enue*umeverg)
      else
         facnue=0.
         xnorme=0.
      end if
      if (enueb.gt.0.) then
         facnueb=(enueb+delta)/enueb
         xnormb=uergg/(facnueb*enueb*umeverg)
      else
         facnueb=0.
         xnormb=0.
      end if
c
      tf=tfermi/utemp
      prefac=4.*3.14159*xsecnn*clight
      dt9=9.*steps
      if (ri.lt.rmaxnue) then
         dtrapnuei = 0.
c
c--e neutrino absorption by neutrons
c  ---------------------------------
c
c-- if beta equilibrium declared, set all changes to zero
c-- particle will be taken care of in nubeta 
         if (ebetaeqi.eq.1) then
            facn=0.
         else
c--compute degeneracy blocking
            call integrals(etai,f0,f1,f2,f3,f4,f5,0,df2,df3)
            call integrals(etanuei,g0,g1,g2,g3,g4,g5,0,dg2,dg3)
            expon=exp(amax1(amin1(enueti/(tempi*utmev)-etai,50.),-50.))
c--   bel: electron end-state blocking
            bel=expon/(1.+expon)
            tneut=amax1(tempi,tf*(xni*rhoi*udens)**(0.66666))
            tfpr=tf*(xpi*rhoi*udens)**(0.66666)
            expon=exp(amin1((tneut-tfpr)/tempi,50.))
c--bpr: proton end-state blocking
            bpr=expon/(1.+expon)
            block=bel*bpr
            fac=prefac*rhoi
            facn=fac*enueti*enueti*ynuei*xni*block
            dy9=ynuei/dt9
            if (facn.gt.dy9) facn=dy9
         endif
         dynuei=dynuei-facn
         enuei=enueti+delta
         facunue=facn*umevnuc*enuei
         dunuei=dunuei-facunue
         dunui=dunui-facunue
         dyei=dyei+facn
      else
c--correct enue and e2nue for gravitational redshift
         fshift=1./gshifti
         cenue=enue*fshift
         ce2nue=e2nue*fshift*fshift
         expon=exp(amin1(amax1(cenue/(tempi*utmev)-etai,-50.),50.))
c--bel: electron end-state blocking
         bel=expon/(1.+expon)
         block=bel
         fac=xsecnn/(ri*ri)
c
c  a) energy absorption
c
         facunue=xni*fac*ce2nue*rlumnue*fshift*facnue*block
         dtrapnuei=facunue*pmassi/fshift
         dunui=dunui - facunue
c
c  b) change in Ye
c
         dyei=dyei + facunue*xnorme
      endif
c
c--e anti-neutrino absorption by protons
c  -------------------------------------
c
      if (ri.lt.rmaxnueb) then
         dtrapnuebi = 0.
         if (pbetaeqi.eq.1) then
            facp=0.
         else
c--since all the energy from the nueb goes into e+ which is never
c--degenerate, the only term that counts is xmuhat
            tprot=amax1(tempi,tf*(xpi*rhoi*udens)**(0.66666))
            tfne=tf*(xni*rhoi*udens)**(0.666666666)
            expon=exp(amin1((tprot-tfne)/tempi, 50.))
c--bne: neutron end-state blocking
            bne=expon/(1.+expon)
            block=bne
            fac=prefac*rhoi
            facp=fac*enuebti*enuebti*ynuebi*xpi*block
            dy9=ynuebi/dt9
            if (facp.gt.dy9) facp=dy9
         endif
         dynuebi=dynuebi-facp
         enuebi=enuebti-deltab
         facunueb=facp*umevnuc*enuebi
         dunuebi=dunuebi-facunueb
         dunui=dunui-facunueb
         dyei=dyei-facp
      else
         fshift=1./gshifti
         ce2nueb=e2nueb*fshift*fshift
         fac=xsecnn/(ri*ri)
c
c  a) energy absorption
c
         facunueb=xpi*fac*ce2nueb*rlumnueb*fshift*facnueb
         dtrapnuebi=facunueb*pmassi/fshift
         dunui=dunui - facunueb
c
c  b) change in Ye
c    
         dyei=dyei - facunueb*xnormb
c
      endif
c
      return
      end
c
      subroutine nubeta(steps,vsoundi,hi,tempi,xmuhati,
     $     etai,rhoi,ynuei,ynuebi,unuei,unuebi)
c*************************************************************
c
c This subroutine treats cases where beta eq. has occurred.
c In beta eq.: munue(beta)=mue-muhat
c so we compute Ynue(munue(beta)) and unue(munue(beta))
c assuming thermal distribution at matter temperature,
c compare with actual Ynue and unue, and move
c things in the right direction
c
c************************************************************
c
c--1/(2.*avo*pi**2*(hbar*c)**3 in Mev-3 cm-3 nucleon g-1)
      parameter(prefac=1.09e7)
c--mn-mp-me in MeV
      parameter (delta=0.783)
      parameter (delta2=delta*delta)
      parameter (deltab=1.805)
      parameter (deltab2=deltab*deltab)
c
      integer ebetaeqi, pbetaeqi
      common /beta/ ebetaeqi, pbetaeqi
      double precision umass
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
      common /konst/ gg, clight, arad, bigr, xsecnn, xsecne
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
c
      yfac=prefac/udens
      ufac=umevnuc*yfac
c-- time scale to eq. = 10 sound crossing time
      deltinv=0.1*steps*(vsoundi/hi)**2
      if (ebetaeqi.eq.1) then
         etabeta=etai-xmuhati/tempi
         call integrals(etabeta,f0,f1,f2,f3,f4,f5,0,df2,df3)
         tmev=tempi*utmev
         tmev2=tmev*tmev
         tmev3=tmev2*tmev
         tmev4=tmev3*tmev
         ynuebeta=yfac*f2*tmev3/rhoi
         unuebeta=ufac*f3*tmev4/rhoi
         dybeta=(ynuebeta-ynuei)*deltinv
         dynuei=dynuei+dybeta
         dyei=dyei-dybeta
         dubeta=(unuebeta-unuei)*deltinv
         dunuei=dunuei+dubeta
         dunui=dunui+dubeta
      endif        
      if (pbetaeqi.eq.1) then
         etabeta=xmuhati/tempi-etai
         call integrals(etabeta,f0,f1,f2,f3,f4,f5,0,df2,df3)
         tmev=tempi*utmev
         tmev2=tmev*tmev
         tmev3=tmev2*tmev
         tmev4=tmev3*tmev
         ynuebeta=yfac*f2*tmev3/rhoi
         unuebeta=ufac*f3*tmev4/rhoi
         dybeta=(ynuebeta-ynuebi)*deltinv
         dynuebi=dynuebi+dybeta
         dyei=dyei+dybeta
         dubeta=(unuebeta-unuebi)*deltinv
         dunuebi=dunuebi+dubeta
         dunui=dunui+dubeta
      endif 
c
      return
      end
c          
      subroutine nuconv(ri,hi,rhoi,tempi,gshifti,pmassi,
     $     ynuei,ynuebi,ynuxi,rmaxnue,rmaxnueb,rmaxnux,
     $     enue,enueb,enux,e2nue,e2nueb,e2nux,rlumnue,rlumnueb,
     $     rlumnux)
c**********************************************************
c
c     This subroutine computes the annhilation of
c     neutrino antineutrino pairs into neutrino antineutrino
c     pairs of other species
c     Calculation follows e+/e- scattering discussion in
c     Mandl and Shaw, p. 315
c     total cross-section is Gf^2*s/(12pi)
c
c*****************************************************
c
      parameter(fs=1./(12.*3.14159))
c-- cross section Gf / gram is 6.02e23*5.29e-44=3.2e-20
      parameter(sigma=fs*3.2e-20)
c
      double precision umass
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
c
      double precision dsigcu
c
      dsigcu=dble(sigma)*umass/dble(udist*udist)
      sigcu=real(dsigcu)
      uconv=1./umevnuc
      tmev=tempi*utmev
      h2i=2.*hi
      if (ri.lt.rmaxnux.and.dnuxi.lt.h2i) then
         expon=exp(min(enuxti/tempnuei-etanuei,50.))
c--   bnue: neutrino end-state blocking
         bnue=expon/(1.+expon)
         expon=exp(min(enuxti/tempnuebi-etanuebi,50.))
c--bnueb: anti-neutrino end-state blocking
         bnueb=expon/(1.+expon)
         blockx=bnue*bnueb
         sigrho=sigcu*rhoi
         facx=sigrho*ynuxi*ynuxi*blockx
c--facnux: number of nux/nuxb annihilation into nue/nueb (per channel)
         facnux=facx*enuxti*enuxti
         outnue=-2.*facnux
         outnux=-.5*outnue
         dynuxi=dynuxi-outnux
         facunux=facnux*umevnuc*2.*enuxti
         outunue=-facunux
         outunueb=-facunux
         outunux=-outunue-outunueb
         dunuxi=dunuxi-outunux
         if (facnux.ne.0.0) then
            if (ri.lt.rmaxnue) then
c--contribution to the nue field if trapped:
               dynuei=dynuei-outnue
               dunuei=dunuei-outunue
            else
c--contribution to the nue luminosities if nue not trapped:
               dunui=dunui-outunue
               shift=gshifti
               tle=0.5*facunux*pmassi*shift
               tenue=enuxti*shift
               te2nue=tenue*tenue
               rlumnue=rlumnue+tle
               enue=enue+tle*tenue
               e2nue=e2nue+tle*te2nue
            endif
            if (ri.lt.rmaxnueb) then
c--contribution to the nueb field if trapped:
               dynuebi=dynuebi-outnue
               dunuebi=dunuebi-outunueb
            else
c--contribution to the nueb luminosities if nueb not trapped:
               dunui=dunui-outunueb
               shift=gshifti
               tleb=0.5*facunux*pmassi*shift
               tenueb=enuxti*shift
               te2nueb=tenueb*tenueb
               rlumnueb=rlumnueb+tleb
               enueb=enueb+tleb*tenueb
               e2nueb=e2nueb+tleb*te2nueb
            endif
         endif
      endif
c--nue-nueb conversion only if both trapped
      if (ri.lt.rmaxnue.and.ri.lt.rmaxnueb.and.
     1     dnuei.lt.h2i.and.dnuebi.lt.h2i) then
         expon=exp(min(.5*(enueti+enuebti)/tempnuxi-
     1        etanuxi,50.))
c--bnux: xneutrino end-state blocking
         bnux=expon/(1.+expon)
         blocke=bnux*bnux
         sigrho=sigcu*rhoi
         face=sigrho*ynuei*ynuebi*blocke
c--facnue: number of nue/nueb annihilation into nux/nuxb (per channel)
         facnue=face*enueti*enuebti
         outnue=2.*facnue
         dynuei=dynuei-outnue
         dynuebi=dynuebi-outnue
         facunue=facnue*umevnuc*enueti
         facunueb=facnue*umevnuc*enuebti
         outunue=facunue
         outunueb=facunueb
         dunuei=dunuei-outunue
         dunuebi=dunuebi-outunueb
c
         if (facnue.ne.0.0) then
            if (ri.lt.rmaxnux) then
c--contribution to the nux field if trapped:
               dynuxi=dynuxi+0.5*outnue
               dunuxi=dunuxi+outunue
            else
c--contribution to the nux luminosities if nux not trapped:
               dunui=dunui-outunue-outunueb
               shift=gshifti
               tlx=(facunue+facunueb)*pmassi*shift
               tenux=0.5*(enueti+enuebti)*shift
               te2nux=tenux*tenux
               rlumnux=rlumnux+tlx
               enux=enux+tlx*tenux
               e2nux=e2nux+tlx*te2nux
            endif
         endif
      endif

c
      return
      end
c
      subroutine nuann(hi,rhoi,tempi,etai,ri,
     $     ynuei,ynuebi,ynuxi,rmaxnue,rmaxnueb,rmaxnux)
c************************************************************
c
c  This subroutine computes the rate of neutrino anti neutrino
c  annihilation into e+/e- pairs
c  see Goodman, Dar, Nussinov, ApJ 314 L7
c
c*********************************************************
c
      parameter(sinw2=0.23)
      parameter(fe=(1.+4.*sinw2+8.*sinw2*sinw2)/(6.*3.14159))
      parameter(fx=(1.-4.*sinw2+8.*sinw2*sinw2)/(6.*3.14159))
c-- cross section Gf / gram is 6.02e23*5.29e-44=3.2e-20
      parameter(sigmae=fe*3.2e-20)
      parameter(sigmax=fx*3.2e-20)
c
      double precision umass
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /neutout/ dyei,dynuei,dynuebi,dynuxi,
     $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
     $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
     $     etanuei,etanuebi,etanuxi,prnui
c
c
      double precision dsigcue, dsigcux
c
      dsigcue=dble(sigmae)*umass/dble(udist*udist)
      dsigcux=dble(sigmax)*umass/dble(udist*udist)
      sigcue=dsigcue
      sigcux=dsigcux
      uconv=1./umevnuc
      tmev=utmev*tempi
      h2i=2.*hi
      if (ri.lt.rmaxnux.and.dnuxi.lt.h2i) then
         expon=exp(amin1(enuxti/tmev-etai,50.))
c--bel: electron end-state blocking
         bel=expon/(1.+expon)
         facx=sigcux*rhoi*ynuxi*ynuxi*bel
         outnux=facx*enuxti*enuxti
         dynuxi=dynuxi-outnux
c--4 because 4 x neutrinos
         facunux=umevnuc*4.*outnux*enuxti
         dunuxi=dunuxi-facunux
         dunui=dunui-facunux
      endif
      if (ipart.eq.1) print *, ri, rmaxnue,rmaxnueb,h2i,dnuei,dnuebi
c--we require BOTH nue and nueb to be trapped
      if (ri.lt.rmaxnue.and.ri.lt.rmaxnueb.and.
     $     dnuei.lt.h2i.and.dnuebi.lt.h2i) then
         expon=exp(amin1(0.5*(enueti+enuebti)/tmev-
     $        etai,50.))
c--bel: electron end-state blocking
         bel=expon/(1.+expon)
         face=sigcue*rhoi*ynuei*ynuebi*bel
         if (ipart.eq.1) print *, 'inside', expon, tmev, etai, 
     $        enueti,enuebti,bel,face,sigcue,rhoi,ynuei,ynuebi
         outnue=face*enueti*enuebti
         dynuei=dynuei-outnue
         dynuebi=dynuebi-outnue
         facunue=umevnuc*outnue*enueti
         facunueb=umevnuc*outnue*enuebti
         dunuei=dunuei-facunue
         dunuebi=dunuebi-facunueb
         dunui=dunui-facunue-facunueb
      endif
c
      return
      end
c
