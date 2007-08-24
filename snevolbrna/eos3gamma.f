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
      double precision umass, uinput,gamma,gamma1
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
      
      gamma=4.d0/3.d0
      gamma1=gamma-1.d0
      ptot=ui*gama1*rhoi
      cs=dsqrt(gamma*ptot/rhoi)
c
      vsoundi=min(cs,0.3333*dble(clight))
      pri=ptot
      temprev=tempi
      rhoprev=rhoi
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


