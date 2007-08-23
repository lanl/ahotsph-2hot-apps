      program rhomap
c*********************************************************************
c*    a program to produce close packed polytropes:
c*    the particle positions are first set up then their
c*    masses are adjusted to fit rho
c*    based on earlier work of M. Herant and M. Davies
c*    12th Feb, 1992 
c*    in this version, updated 26th April, 1993, say only pmass
c*    propto rho 
c*********************************************************************

      parameter (idim=1000000)

      character *10 reply

      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc

c*****set up sphere of particles
10    call cpack

      write(*,*) 'total number of particles: ',npart
      write(*,*) 'is this OK?'
      read(*,100) reply
      if ((reply.ne.'y').and.(reply.ne.'Y')) goto 10

100   format(a10)

c*****now adjust particle masses to fit rho
      call init
      print *,'calling readrho'
      call readrho
      print *,'calling rhotheo'
      call rhotheo
c*****set initial mass
      call initmass
      
      call intern
c**** call wascii
      call wdump

      close(20)

      stop
      end

      subroutine init
c************************************************************
c                                                           *
c  this routine initialises variables for rho mapping       *
c                                                           *
c************************************************************
c
      parameter (idim=1000000)
      common /smooth/ h0
      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc
      common /ener2/ trot, tkin, tgrav, tterm, xunuc
      common /grav/ dgrav(idim),npm
      common /bodys/ n1, fmas1, n2, fmas2
      common /time / t, dt
      common /cgas / gamma
      common /carac/ pmass(idim), matter(idim)
c
      pmass0=1./npart
      write(*,*) 'in init, h0 :',h0
c*****set up arrays
      do i=1,npart
         pmass(i)=pmass0
         matter(i)=0
         dgrav(i)=0.0
         u(i)=0.0
         h(i)=h0
         vx(i)=0.0
         vy(i)=0.0
         vz(i)=0.0
      end do
      npm=0
      n1=npart
      n2=0
      escap=0.0
      gamma=1.+1./1.5
      t=0.0
      tgrav=0.0
      tkin=0.0
      tterm=0.0
      return
      end


c***********************************************************

      subroutine readrho

      character*80 s1,s2
      common /map/ rhotab(10000),grhotab(10000),ptab(10000),step

      open(19,file='poly1.5')
      print *,'in readrho'
c      read(19,*) s1
      read(19,*) imax,step,pindex,tmass,radius
      print *,imax,step,pindex,tmass,radius
c      read(19,*) s2
      print *,'stepsize=',step
      do 10 i=1,imax
         read(19,*) r, rhotab(i), ptab(i)
   10 continue
      close(19)
      print *,'got ',imax,' densities and pressures'
      do 20 i=imax+1,10000
         rhotab(i)=1.0e-10
         ptab(i)=1.0e-25
   20 continue

      grhotab(1)=0.0
      grhotab(2)=(rhotab(2)-8.0*rhotab(1)+
     1        8.0*rhotab(3)-rhotab(4))/(12.0*step)
      do 30 i=3,9998
         grhotab(i)=(rhotab(i-2)-8.0*rhotab(i-1)+
     1           8.0*rhotab(i+1)-rhotab(i+2))/(12.0*step)
   30 continue
      grhotab(9999)=0.0
      grhotab(10000)=0.0

      return
      end

c*************************************************************
      subroutine rhotheo

      parameter(idim=1000000)

      common /map/ rhotab(10000),grhotab(10000),ptab(10000),step
      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc
      common /rho/ rhot(idim)

      do 20 i=1,npart
         r=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         rstep=r/step
         index1=int(rstep)+1
         index2=int(rstep)+2
         if (index2 .lt. 10000) then
            t=rstep-index1+1.0
            rhot(i)=(1-t)*rhotab(index1)+t*rhotab(index2)
         else
            rhot(i)=1e-10
            write(*,*) index1,index2,i,r,step
            stop
         endif
  20  continue

      return
      end

      subroutine initmass

c*****set initial masses such that pmass is proportional to
c*****rhotheo

      parameter(idim=1000000)
      parameter(neighmx=60)


      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc
      common /rho/ rhot(idim)

      common /carac/ pmass(idim), matter(idim)

      sum=0.

      do 10 i=1,npart
         pmass(i)=rhot(i)*h(i)**3
         sum=sum+pmass(i)
         if (pmass(i).eq.0.) write(69,*) i,pmass(i),h(i),rhot(i)
10    continue

      write(*,*) 'sum after init is....',sum
      do 20 i=1,npart
         pmass(i)=pmass(i)/sum
20    continue

      return
      end

      subroutine intern

      parameter (idim=1000000)

      common /map/ rhotab(10000),grhotab(10000),ptab(10000),step
      common /cgas / gamma
      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc
      common /rho/ rhot(idim)
      common /carac/ pmass(idim), matter(idim)

      tterm=0.0
      do 20 i=1,npart
         r=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         rstep=r/step
         index1=int(rstep)+1
         index2=int(rstep)+2
         if (index2 .lt. 10000) then
            t=rstep-index1+1.0
            press=(1-t)*ptab(index1)+t*ptab(index2)
         else
            press=1e-25
         endif
         u(i)=press/rhot(i)/(gamma-1.)
  20  continue
      write(*,*) 'gamma is ....',gamma
c*****calculate thermal energy
      do i=1,npart
         tterm=tterm+u(i)*pmass(i)
      end do
      write(*,*) 'tterm is ...', tterm

      return
      end


      subroutine wascii

      parameter(idim=1000000)
      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc
      common /carac/ pmass(idim), matter(idim)
      common /densi/ rho(idim), temp(idim)
      common /rho/ rhot(idim)
      character*30 string
      save kount
      data kount/0/

      k10=kount/10
      k1=kount-k10*10
      string='wasciiE'//char(48+k10)//char(48+k1)
      open(30,file=string)
      do 10 i=1,npart
         r=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         write(30,*) r,pmass(i),rho(i),rhot(i)
   10 continue
      close(30)
      kount=kount+1
      return
      end



      subroutine wdump
c************************************************************
c                                                           *
c  this routine writes a dump on disk                       *
c                                                           *
c************************************************************
c
      parameter (idim=1000000)
      parameter (nel=1)
c
      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc
      common /rho/   rhot(idim)
      common /carac/ pmass(idim), matter(idim)
      common /typef/ igrp, igphi, ifsvi, ifcor, ichoc, iener,
     1               ibound, iexf, iburn, idiff, ipenet, iexpan
      common /compo/ cc(nel,idim), ccave(nel), eccm
      common /cgas / gamma
      common /fracg/ fgas, escap
      common /kerne/ cnormk, ftomas, fnbtot
      common /time / t, dt
      common /bodys/ n1, fmas1, n2, fmas2

      common /ener2/ trot, tkin, tgrav, tterm, xunuc
      common /recor/ irec, ipos
      common /ptpos/ px,py,pz
      common /point/ iptms,pvx,pvy,pvz,ptms
      common /grav/ dgrav(idim),npm
c
      character*7 where
      character*12 filen
c
      data where/'wdump'/
c
c--write
c

      if (t .lt. 0) then
         t=0.0
         do 5 i=1,npart
            vx(i)=0.0
            vy(i)=0.0
            vz(i)=0.0
   5     continue
      endif

      irec=irec + 1
      iptms=0
      ichkl=0
c
      write(*,*) 'name of output file'
      read(*,101) filen
  101 format(a12)
      open(20,file=filen,form='unformatted')

      write(20,iostat=io)npart,
     1   t,gamma,(h(i),i=1,npart),tkin,
     2   tgrav,tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,
     3   npart),(vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,
     4   npart),(u(i),i=1,npart),(pmass(i),i=1,npart),
     5   (rhot(i),i=1,npart)
c
      return
      end


      subroutine cpack
c*****construct a close-packed sphere
c*****do it the easier way, rather than using rotations

      parameter (idim=1000000)

      common /smooth/ h0
      common /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), h(idim), unuc

      write(*,*) 'please give radius, particle spacing,and h'

      read(*,*) r,d,h0

      r2=r*r

      dx=0.5*d
      dy=d/sqrt(12.)
      dz=d*sqrt(2./3.)

      xstep=0.5*d
      ystep=xstep*sqrt(3.)
      zstep=2.*dz 

      rbig=r+d

      ipart=0

      iy=rbig/ystep
      ix=rbig/xstep
      iz=rbig/zstep

      xlim=ix*xstep
      ylim=iy*ystep
      zlim=iz*zstep

c*****if iy is even, then have to offset x by d/2
      if (mod(iy,2).eq.0) xlim=xlim+0.5*d

      factor=0.5

      do y0=-ylim,ylim,ystep
         xlim=xlim+factor*d
         factor=-factor
         do x0=-xlim,xlim,d
            do z0=-zlim,zlim,zstep
               ri2=x0*x0+y0*y0+z0*z0
               if (ri2.lt.r2) then
                  ipart=ipart+1
                  x(ipart)=x0
                  y(ipart)=y0
                  z(ipart)=z0
               end if
               x1=x0+dx
               y1=y0+dy
               z1=z0+dz
               ri2=x1*x1+y1*y1+z1*z1
               if (ri2.lt.r2) then
                  ipart=ipart+1
                  x(ipart)=x1
                  y(ipart)=y1
                  z(ipart)=z1
               end if
            end do
         end do
      end do

      npart=ipart

      return
      end
