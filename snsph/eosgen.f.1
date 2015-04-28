      subroutine eosgen(rhoi,tempi,yei,abari,ui,u2i,pri,xpi,
     $     xni,ufreezi,ifleosi,iident,iprocnum)
c
c--rhoi,tempi,yei,abari are inputs
c--ui,u2i,pri,xpi,xni,ufreezi are outputs
c
      implicit double precision (a-h,o-z)

      integer iflag,ifleosi,iident
      double precision inpvar(4), xa, xh, yeh

      tkelv=tempi*1.d9
      rhocgs=rhoi*2.d6

c The number here must match the one in eosfl
      if (tkelv.lt.8.d9) then
         iflag=1
      else
         if (rhocgs.lt.5.d10) then
            iflag=2
         else
            iflag=3
         endif
      endif
      if (iflag.eq.1) then
         zbar=yei*abari
         call coulomb2(rhocgs,zbar,yei,ucoul,pcoul)
         tno=tkelv/1d9
         rhono=rhocgs/1d7
         call nados(tno,rhono,zbar,abari,pel,eel,sel,
     1              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)
         ucgs=etot*1d17+ucoul
         pcgs=ptot*1d24+pcoul
         scgs=stot*1d8
         xpi=0.0
         xni=0.0
      elseif (iflag.eq.2) then
         t9=tkelv/1d9
         call nsestart(t9,rhocgs,yei,xpi,xni,iident,iprocnum)
         call nsetemp(t9,rhocgs,yei,t9,xpi,xni,xa,xh,yeh,
     1                zbar,abari,ubind,dubind,iident,iprocnum)
         zbar=yei*abari
         call coulomb2(rhocgs,zbar,yei,ucoul,pcoul)
         tno=tkelv/1d9
         rhono=rhocgs/1d7
         call nados(tno,rhono,zbar,abari,pel,eel,sel,
     1              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)
         ucgs=etot*1d17+ucoul+ubind
         pcgs=ptot*1d24+pcoul
         scgs=stot*1d8
      elseif (iflag.eq.3) then
         tswe=tkelv/1.16d10
         inpvar(1)=tswe
         inpvar(2)=0.155d0
         inpvar(3)=-15.0d0
         inpvar(4)=-10.d0
         brydns=rhocgs*6.02d-16
         pprev=yei*brydns
         call slwrap(inpvar,yei,brydns,pprev,
     $        psl,usl,dusl,gamsl,eta,xpi,xni,
     $        xa,xh,yeh,abar,xmuh,ssl)
         abari=1.d0
         ucgs=usl*9.644d18
         pcgs=psl*1.602d33
         scgs=ssl*8.25d7
      endif
      ui=ucgs/1.d16
      pri=pcgs/2.d22
      u2i=scgs/1.d7
      ifleosi=iflag
      if (iflag.eq.1.) then
         if (abari.lt.2.5) then
            ufr=-3.3d17/1.d16
         elseif (abari.lt.4.5) then
            ufr=-6.1d17/1.d16
         elseif (abari.lt.20.) then
            ufr=-7.7d18/1.d16
         else
            ufr=-8.4d18/1.d16
c            ufr=-8.2d18/1.d16
         endif
      else
         ufr=0.0
      endif
      ufreezi=ufr
c
      return
      end

      subroutine coulomb2(rhoi,zbar,ye,ucoul,pcoul)
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
c     parameter(ufac=-0.214d0)
c     parameter(pfac=-0.0714d0)
      parameter(ufac=-1.70d13)
      parameter(pfac=-5.67d12)
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


