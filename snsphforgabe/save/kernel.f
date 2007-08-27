      subroutine kernel 
c***********************************************************
c                                                          *
c  This subroutine generates a smooth polynomial spline    *
c  kernel W(u/h) with the following optional properties:   *
c                                                          *
c  * u = sqrt (sum ((i=1,d) x_i^2)), d = dimension of f    *
c                                                          *
c  * The spline is given on two domains [0,h] and [h,2h].  *
c                                                          *
c  * DW ---> 0 at u=0 and u=2h     (Flatness)              *
c                                                          *
c  * W  ---> 0 at u=2h             (Compactness)           *
c                                                          *
c  * W is in C^n (0,2h)            (Smoothness)            *
c                                                          *
c  * Min DW is at least stiff*V_d < 0     (Stiffness)      *
c                                                          *
c  * W yields an Nth order (N even) kernel estimate of any *
c  analytic function of d variables.                       *
c                                                          *
c  These requirements for W yield L linearly independent   *
c  linear constraints on W.  Each constraint               *
c  equation gives a linear relation between the            *
c  coefficients of the polynomial spline.                  *
c                                                          *
c  On the primary interval [0,h]                           *
c                                                          *
c      W(u) = 1/Vd (sum (i=0,i_alpha) alpha(i) (u/h)^i )   *
c                                                          *
c  On the second interval [h,2h]                           *
c                                                          *
c      W(u) = 1/Vd (sum (i=0,i_beta) beta(i) (u/h)^i )     *
c                                                          *
c  where Vd = 'volume of an d-Ball'  and                   *
c  i_alpha + i_beta = L.                                   *
c                                                          *
c  After the options are selected (flatness, stiffness,    *
c  smoothness, compactness), then the coefficients alpha(i)*
c  and beta(i) are computed by inverting an L x L          *
c  matrix (for the Nth order kernel estimate.)             *
c  The code determines the number of coefficients which    *
c  are necessary to satisfy all of the constraints on the  *
c  kernel W.  The options allow for arbitrary stiffness    *
c  (or none), arbitrary smoothness of W, and a choice of   *
c  compact, noncompact, flat or nonflat kernels, and ones  *
c  which yield estimate for functions of R^d.              *
c                                                          *
c***********************************************************
c                                                          *
c  created: 8.01.88  Greg Willette     KERN (F77)          *
c  revised: 9.20.89  GTW               KERNEL (F77)        *
c                                                          * 
c  KERNEL generater written in CRAY fortran.               *
c                                                          *
c  Load lib=(none)                                         *
c                                                          *
c***********************************************************
c
      parameter (iker=50)
      parameter (idum=20)
      parameter (ncoefs=6)
c
      implicit double precision (a-h,o-z)
c
      common /kernn/ ncoef1, ncoef2
      common /kernc/ wcoef1(ncoefs), wcoef2(ncoefs)
      common /kern1/ stiff,wrstf,ismw,inor 
      common /kern2/ iflat, icompact, istiff, ismooth, iorder
c
      dimension gamma(iker), trix(iker,iker), trix1(iker,iker)
c
      character*7 iflat,icompact,istiff,ismooth,iorder
c
c--set options for kernel
c
c  the standard 2D SPH kernel is given by:
c  iflat=flat; icompact=compact; istiff=stiff; ismooth=smfin;
c  iorder==norder
c  -1.333333333 0.666666666666 3 0 3
c
c--peaked kernel
c     iflat='noflat'
c     icompact='compact'
c     istiff='stiff'
c     ismooth='sminf'
c     iorder='norder'
c     stiff=0.
c     wrstf=2.
c     ismw=3
c     inor=0
c     idimf=2
c
c--4rth order  (sl)
c     iflat='noflat'
c     icompact='compact'
c     istiff='stiff'
c     ismooth='sminf'
c     iorder='norder'
c     stiff=0.
c     wrstf=1.999999
c     ismw=3
c     inor=2
c     idimf=2
c
c-- yet another kernel (si01)
c     iflat='flat'
c     icompact='compact'
c     istiff='stiff'
c     ismooth='smfin'
c     iorder='norder'
c     stiff=0.0
c     wrstf=2.
c     ismw=3
c     inor=2
c     idimf=2
c
c--another kernel (sk)
      iflat='flat'
      icompact='compact'
      istiff='stiff'
      ismooth='smfin'
      iorder='norder'
      stiff=-1.75
      wrstf=0.3
      ismw=4
      inor=0
      idimf=2
c
c--classic kernel
c     iflat='flat'
c     icompact='compact'
c     istiff='stiff'
c     ismooth='smfin'
c     iorder='norder'
c     stiff=-10./7.
c     wrstf=0.66666666666666666666
c     ismw=3
c     inor=0
c     idimf=2
c
c--Compute Size of Matrix
c
      imat=0
      if(iflat(1:4).eq.'flat') then
         imat=imat+2
      end if
      if(icompact(1:7).eq.'compact') then
         imat=imat+1
      end if
      if(istiff(1:5).eq.'stiff') then
         imat=imat+1
      end if
      if(iorder(1:6).eq.'norder') then
         imat=imat+1+inor/2
      end if
      if(ismooth(1:5).eq.'sminf') then
         imat=imat*2
      end if
      if(ismooth(1:5).eq.'smfin') then
         imat=imat+ismw
      end if
c 
c--Compute number of alphas and betas
c
      ibeta=int(imat/2.d0)
      ialpha=imat-ibeta
      ialph1=ialpha+1
c
c--Compute Matrix
c
      irow=1
c
c--Flatness
c
      if(iflat(1:4).eq.'flat') then
         do icol=1,ialpha
            if(icol.eq.2) then
               trix(irow,icol)=1.d0
            else
               trix(irow,icol)=0.d0
            end if
         enddo
         do icol=ialph1,imat
            trix(irow,icol)=0.d0
         enddo
         gamma(irow)=0.d0
         irow=irow+1
         do icol=1,ialpha
            trix(irow,icol)=0.d0
         enddo
         do icol=ialph1,imat
            ibetr=ibeta+icol-imat-1
            trix(irow,icol)=(ibetr)*(2**(ibetr-1))
         enddo
         gamma(irow)=0.d0
         irow=irow+1
      end if
c
c--Compactness
c
      if(icompact(1:7).eq.'compact') then
         do icol=1,ialpha
            trix(irow,icol)=0.d0
         enddo
         do icol=ialph1,imat
            ibetr=icol+ibeta-imat-1
            trix(irow,icol)=2**(ibetr)
         enddo
         gamma(irow)=0.d0
         irow=irow+1
      end if
c
c--Stiffness
c
      if(istiff(1:5).eq.'stiff') then
         do icol=1,ialpha
             trix(irow,icol)=(icol-1)*((wrstf)**(icol-2))
         enddo
         do icol=ialph1,imat
             trix(irow,icol)=0.d0
         enddo
         gamma(irow)=stiff
         irow=irow+1
      end if
c
c--Smoothness
c
      if(ismooth(1:5).eq.'sminf') then
         ismrow=irow+ialpha-1
         do irowc=irow,ismrow
            intemp=irowc-irow
            do icol=1,ialpha
               if(icol.le.intemp) then
                  trix(irowc,icol)=0.d0
               else
                  call factor(icol-1,fac1)
                  call factor(icol-intemp-1,fac2)
                  trix(irowc,icol)=fac1/fac2
               end if
            enddo
            do icol=ialph1,imat
               icltmp=icol-ialpha
               ifac=icltmp-intemp-1
               if(icltmp.le.intemp) then
                  trix(irowc,icol)=0.d0
               else
                  call factor(icltmp-1,fac1)
                  call factor(ifac,fac2)
                  trix(irowc,icol)=-1.d0*fac1/fac2
               end if
            enddo
            gamma(irowc)=0.d0
         enddo
         irow=irow+ialpha
      end if
c
      if(ismooth(1:5).eq.'smfin') then
         ismrow=irow+ismw-1
         do irowc=irow,ismrow
            intemp=irowc-irow
            do icol=1,ialpha
               if(icol.le.intemp) then
                  trix(irowc,icol)=0.d0
               else
                  ifac=icol-intemp-1
                  call factor(icol-1,fac1)
                  call factor(ifac,fac2)
                  trix(irowc,icol)=fac1/fac2
               end if
            enddo
            do icol=ialph1,imat
               icltmp=icol-ialpha
               ifac=icltmp-intemp-1
               if(icltmp.le.intemp) then
                  trix(irowc,icol)=0.d0
               else
                  call factor(icltmp-1,fac1)
                  call factor(ifac,fac2)
                  trix(irowc,icol)=-1.d0*fac1/fac2
               end if
            enddo
            gamma(irowc)=0.d0
         enddo
         irow=irow+ismw
      end if
c
c--Accuracy
c
      if(iorder(1:6).eq.'norder') then
         inorow=irow+inor/2
         if(idimf.eq.1)then
            geotemp=2.0d0
         else
            geotemp=1.0d0
         end if
         do irowc=irow,inorow
            intemp=2*(irowc-irow)
            do icol=1,ialpha
               trix(irowc,icol)=(geotemp/(icol+intemp+idimf-1))
            enddo
            do icol=ialph1,imat
               icltmp=icol-ialpha
               twofix=((2**(icltmp+intemp+idimf-1))-1.d0)
               trix(irowc,icol)=(geotemp/(icltmp+intemp+idimf-1))*twofix
            enddo
            if(irowc.eq.irow) then
               gamma(irowc)=1.0d0/idimf
            else
               gamma(irowc)=0.0d0
            end if
         enddo
         irow=irow+inor/2
      end if
c
c--Invert Matrix by Gauss-Jordan Reduction
c
c
      do i=1,irow
         do j=1,imat
            if(i.eq.j) then
               trix1(i,j)=1.d0
            else
               trix1(i,j)=0.d0
            end if
         enddo
      enddo
      do i=1,irow
         if(trix(i,i).ne.0.d0) then
            go to 220
         else
            j=i+1
  210       continue
            if(trix(j,i).ne.0.d0) then
               do itemp=1,imat
                  temp=trix(i,itemp)
                  trix(i,itemp)=trix(j,itemp)
                  trix(j,itemp)=temp
                  temp=trix1(i,itemp)
                  trix1(i,itemp)=trix1(j,itemp)
                  trix1(j,itemp)=temp
               enddo
               goto 220
            else
               j=j+1
               if(j.gt.imat) then
                  write(6,*)'singular matrix'
                  stop
               end if
               go to 210
            end if
         end if
  220    continue
         sub=trix(i,i)
         do idiv=1,imat
            trix1(i,idiv)=trix1(i,idiv)/sub
            trix(i,idiv)=trix(i,idiv)/sub
         enddo
         do isub=1,irow
            if(isub.ne.i) then
               sub0=trix(isub,i)
               do irun=1,imat
                  sub1=sub0*trix1(i,irun)
                  trix1(isub,irun)=trix1(isub,irun)-sub1
                  sub=sub0*trix(i,irun)
                  trix(isub,irun)=trix(isub,irun)-sub
               enddo
            end if
         enddo
      enddo
c
c--Compute coefficients
c
      ncoef1=ialpha
      do i=1,ncoef1
         wcoef1(i)=0.0d0
      enddo
      ncoef2=ibeta
      do i=1,ncoef2
         wcoef2(i)=0.0d0
      enddo
c
      do i=1,ncoef1
         do j=1,imat
            wcoef1(i)=wcoef1(i)+trix1(i,j)*gamma(j)
         enddo
      enddo
c
      do i=1,ncoef2
         ishift=i+ialpha
         do j=1,imat
            wcoef2(i)=wcoef2(i)+trix1(ishift,j)*gamma(j)
         enddo
      enddo
c
   99 return
      end
c
      subroutine factor(i,fac) 
      implicit double precision (a-h,o-z)
      fac=1.d0
      do j=1,i
         fac=j*fac
      end do
      return
      end
c
      subroutine ktable
c***********************************************************
c                                                          *
c  This subroutine builds a table for the various values   *
c  of the kernel, and the gradient of the kernel.          *
c                                                          *
c***********************************************************
c
      parameter (ncoefs=6)
      parameter (itable=80000)
c
      double precision wcoef1, wcoef2, v2, v, w, dw
c
      common /kernn/ ncoef1, ncoef2
      common /kernc/ wcoef1(ncoefs), wcoef2(ncoefs)
      common /kerne/ cnormk
      common /table/ wij(0:itable), grwij(0:itable), dvtable
      common /typef/ iextf, ieos
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
      common /debug/ idebug, itrace
c
      character*7 idebug
      character*3 itrace
c
      data pi/3.141592654/
c
c--allow for tracing flow
c
      if(itrace.eq.'all')write(iprint,250)
  250 format(' entry subroutine ktable')
c
c--generate kernel
c
      call kernel
      write(*,*)ncoef1
      do i=1,ncoef1
         write(*,255)wcoef1(i)*7.0/10.0
      enddo
      write(*,*)ncoef2
      do i=1,ncoef2
         write(*,255)wcoef2(i)*7.0/10.0
      enddo
 255  format('f :', f18.14)
c
c--maximum interaction length and step size
c
      v2max=4.
      dvtable=v2max/itable
      i1=1./dvtable
c
c--normalisation constant
c
      cnormk=1./pi
c
c--build tables
c
c  a) v less than 1
c
      wij(0)=cnormk*wcoef1(1)
      grwij(0)=0.
      do i=1,i1
         v2=i*dvtable
         v=dsqrt(v2)
         w=wcoef1(ncoef1)
         do j=ncoef1-1,1,-1
            w=w*v + wcoef1(j)
         enddo
         dw=(ncoef1-1)*wcoef1(ncoef1)
         do j=ncoef1-2,1,-1
            dw=dw*v + j*wcoef1(j+1)
         enddo
         wij(i)=cnormk*w
         grwij(i)=cnormk*dw/v
      enddo
c
c  b) v greater than 1
c
      do i=i1+1,itable
         v2=i*dvtable
         v=dsqrt(v2)
         w=wcoef2(ncoef2)
         do j=ncoef2-1,1,-1
            w=w*v + wcoef2(j)
         enddo
         dw=(ncoef2-1)*wcoef2(ncoef2)
         do j=ncoef2-2,1,-1
            dw=dw*v + j*wcoef2(j+1)
         enddo
         wij(i)=cnormk*w
         grwij(i)=cnormk*dw/v
      enddo
c
      return
      end
