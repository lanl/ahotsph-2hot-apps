      subroutine abinit(irank,netrcfile)

      implicit none

      include 'dimen'
      include 'comod'

      integer*4 irank
      character*20 netrcfile

      real*8 y(ndim),xx(ndim)

c..   solar system abundance data
      integer*4 iagdim
      parameter( iagdim = 286 )
      integer*4 nzag(iagdim),naag(iagdim)
      real*8 xxag(iagdim)
      character*2 chag(iagdim)
      integer*4 k, i,j,n

      include 'crate'
      include 'caeps'
      include 'cgen'
      include 'cburn'

      real*8 summ, zhe

      integer*4 izbu(ndim),inbu(ndim),ibu,idummy
      character*5 cbu

      integer*4 nsp
      parameter( nsp=20 )
c..Z and N for special nuclei
      integer*4 nspz(nsp), nspn(nsp)
c..   nscr=scratch array for index reordering
c cie: nscr should have N=nucpg elements?
      integer*4 nscr(nucpg),iscr,itno,jz,ja
 
      integer*4 lun

      logical tobe

c..default values for new network; used to reset abundances
      data zpop/0.01886d0/,zhyd/0.70683d0/

c..special nuclei to be used in diagnostics and i/o
      data nspz/ 1, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 26, 
     1     26, 28, 28, 28, 28/
      data nspn/ 0, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28, 28, 30,
     1     32, 30, 32, 34, 36/

c-------------------------------------------------------------------

c      write(*,*) netrcfile

      if(irank.eq.0) write(*,*)'ENTERING ABINIT, netsize ',netsize
ccccccccccccccccccccccccccccc

c..   resource file for analysis programs
      inquire(file=trim(netrcfile),exist=tobe)
      if( .not. tobe )then
         if(irank.eq.0) write(*,*)'abinit: no ',trim(netrcfile)
         newnet = 1
      endif

      if( newnet .ne. 0 )then

         if(irank.eq.0) write(*,*)
     1      "abinit: newnet =",newnet,', redefining abundances'

         do j = 1, netsize
            cnuc(j) = xid(j)
            lz(j) = nz(j)
            ln(j) = nn(j)
         enddo
       
         do j = netsize +2, ndim
            cnuc(j) = '     '
            lz(j) = 0
            ln(j) = 0
         enddo
         cnuc(netsize+1) = '   Ye'

c..   solar  "mass" (nucleon) fractions
         open(10,file='ssystem.dat')
c..naag is the "mass number" A=Z+N, not the mass in amu
         do i = 1, iagdim
            read(10,'(i5,i3,a3,i4,1pe11.3)')idummy,nzag(i),chag(i),
     1           naag(i),xxag(i)
         enddo
         do n = 1, netsize+1
            xx(n) = 0.0d0
         enddo
         do n = 1, netsize
            do i = 1, iagdim
               if( lz(n) .eq. nzag(i) .and.
     1              ln(n)+lz(n) .eq. naag(i))then
                  xx(n) = xxag(i)
               endif
            enddo
         enddo
         do n = 1, netsize
            solarx(n) = xx(n)
         enddo
c..   determine which nuclei are stable
c..   chose them for pgplot.f and for dout3.f
         do j = 1, nucpg
            nucp(j) = 0
         enddo
c..   always choose H1 and He4
         nucp(1) = netsize-1
         nucp(2) = netsize
         j = 2
         do n = 1, netsize
            if( xx(n) .gt. 0.0d0 .and. j+1 .le. nucpg )then
               j = j+1
               nucp(j) = n
c               if(irank.eq.0) write(*,'(5i5,a5,1p8e12.3)')
c     1              j,n,nucp(j),lz(n),ln(n),
c     1              cnuc(n),xx(n)
            endif
         enddo

c..   number of nonradioactive nuclei
         mnucpg = j


c..force values if new network (ag88m)
         zhyd = 0.0d0
         zpop = 0.0d0
         zhe  = 0.0d0
         do n = 1, nnuc
            if( lz(n) .eq. 1 )then
               zhyd = zhyd + solarx(n)
            elseif( lz(n) .eq. 2 )then
               zhe  = zhe  + solarx(n)
            else
               zpop = zpop + solarx(n)
            endif
         enddo

c         write(*,'(a30,2(a5,1pe12.3))')'forcing values:','zpop',zpop,
c     1        'zhyd',zhyd

c..   adjust elemental helium to fit
         do n = 1, netsize
            if( lz(n) .eq. 2 )then
               xx(n) = xx(n)*(1.0d0-zhyd-zpop)/zhe
            endif
         enddo

         summ = 0.0d0
         do n = 1, netsize
            summ = summ + xx(n)
         enddo
c         write(*,'(a30,1pe12.3)')'error after He adjustment:',
c     1        summ-1.0d0

c..   mole fractions
c.. cie: I think this block is unnecessary?!
c         do n = 1, netsize
c            xa(n) = nz(n) + nn(n) + qex(n)/931.487d0
c            y(n) = xx(n)/xa(n)
c         enddo
c         y(netsize+1) = 0.0d0
c         do n = 1, netsize
c            y(netsize+1) =  y(netsize+1) + y(n)*dble(lz(n))
c         enddo
c         if(irank.eq.0) write(*,*)'Ye =',y(netsize+1)

c..   spead over spatial grid 
c         do k = 2, kk+1
c            do n = 1, netsize+1
c               x(n,k) = y(n)
c            enddo
c         enddo

c..   redefine net.rc values
         lun = 30+irank
         open(lun,file=trim(netrcfile))
         do j = 1, netsize
            write(lun,'(3i5,a5,0pf10.4,1pe12.4)')
     1           j,nz(j),nn(j),xid(j),qex(j),solarx(j)
         enddo
c..   write whole array, including zeros, for ease in reading
c..   by many routines
         do j = 1, nucpg, 10
            write(lun,'(10i5)')(nucp(i),i=j,j+9)
         enddo


         close(lun)

         if(irank.eq.0) write(*,*)
     1      'abinit: new network and abundances defined'
         newnet = 0

      else

         if(irank.eq.0) write(*,*)'abinit: checking netrc'

         open(lun,file=trim(netrcfile),status='old')
         ibu = 1
c..initializing qex and solarx
 100     read(lun,'(3i5,a5,0pf10.4,1pe12.4)',end=101)
     1        idummy,izbu(ibu),inbu(ibu),cbu,qex(ibu),solarx(ibu)
         if( izbu(ibu) .ne. nz(ibu) .or.
     1        inbu(ibu) .ne. nn(ibu) )then
c..   different net.rc
            if(irank.eq.0) write(*,*)trim(netrcfile)
            if(irank.eq.0) write(*,'(3i5,a5)')
     1         ibu,izbu(ibu),inbu(ibu),cbu
            if(irank.eq.0) write(*,*)'abinit'
            if(irank.eq.0) write(*,'(3i5,a5)')
     1         ibu,nz(ibu),nn(ibu),xid(ibu)
            if(irank.eq.0) write(*,*)
     1         'abinit: different ',trim(netrcfile)
            if(irank.eq.0) write(*,*)
     1         'rm ',trim(netrcfile),' or change its name, and rerun'
            stop'abinit: reading net.rc'
         endif
         if( izbu(ibu) .eq. 2 .and. inbu(ibu) .eq. 2 )goto 101
         ibu = ibu + 1
         goto 100
 101     continue
c..   net.rc is consistent
         do j = 1, netsize
            cnuc(j) = xid(j)
            lz(j) = nz(j)
            ln(j) = nn(j)
         enddo
       
         do j = netsize +2, ndim
            cnuc(j) = '     '
            lz(j) = 0
            ln(j) = 0
         enddo
         cnuc(netsize+1) = '   Ye'

         read(lun,'(10i5)')nucp

c..   determine nonzero entries
         mnucpg = 0
         do j = 1, nucpg
            if( nucp(j) .ne. 0 )then
               mnucpg = mnucpg + 1
            endif
         enddo

c..define solar metallicity from solar tables for consistency
         zsol = 0.0d0
         do j = 1, ibu
            if( lz(j) .ne. 1 .and. lz(j) .ne. 2 )then
               zsol = zsol + solarx(j)
            endif
         enddo

         if(irank.eq.0) write(*,*)
     1      'abinit: ',trim(netrcfile),' exists and is left unchanged'
         close(lun)
      endif

c..   define Ye
      do k = 2, kk
         x(netsize+1,k) = 0.0d0
         do n = 1, netsize
            x(netsize+1,k) = x(netsize+1,k) + x(n,k)*dble(lz(n))
         enddo
         if( x(netsize+1,k) .lt. 0.0d0 
     1        .or. x(netsize+1,k) .gt. 1.0d0)then
            if(irank.eq.0) write(*,*)
     1         ' abinit: Ye error, k ', x(netsize+1,k), k
            stop'abinit: Ye error'
         endif
      enddo

c..adjust nucp array to control which nuclei are monitored.............

c..find special nuclei
      j = 0
      do i = 1, nsp
         do n = 1, netsize
            if( nspz(i) .eq. lz(n) .and. nspn(i) .eq. ln(n) )then
               nucp(i) = n
               j = j + 1
               go to 110
            endif
         enddo
         if(irank.eq.0) write(*,'(a8,i5,a3,i5,a3,i5,a15)')'nucleus',i,
     1        'Z',nspz(i),'N',nspn(i),'not found'
 110     continue
      enddo
      if(irank.eq.0) write(*,*)j
      if( j .eq. nsp )then
         do i = 1, nsp
               if(irank.eq.0) write(*,'(3i5,a5,1p8e12.3)')
     1            i,nspz(i),nspn(i),
     1              cnuc(nucp(i))
         enddo
         if(irank.eq.0) write(*,*)'ABINIT: All special nuclei found'
      else
         if(irank.eq.0) write(*,*)
     1      'ABINIT: ',nsp-mnucpg,'  special nuclei NOT found'
      endif

c      j = mnucpg
c..   j is number of special nuclei which were found
      do n = 1, netsize
         if( solarx(n) .gt. 0.0d0 .and. j+1 .le. nucpg )then
            do i = 1, j
               if( nucp(i) .eq. n )then
c..   avoid duplication
                  go to 120
               endif
            enddo
            j = j+1
            nucp(j) = n
         endif
 120     continue
      enddo
c..j has been increased to include some stable isotopes, up to nucpg=60
      if(irank.eq.0) write(*,*)j,' nuclei chosen'
      if(irank.eq.0) write(*,'(20a6)')(cnuc(nucp(i)),i=1, j)

c..reorder by charge and atomic number
      if(irank.eq.0) write(*,*)'begin reordering'
      do n = 1, j
         nscr(n) = nucp(n)
      enddo
      do n = j+1, ndim
         nscr(n) = -1
      enddo

c      do j = 1, nucpg !cie: why is this nucpg?? lz/ln are declared in
cburn to have ndim elements
      do j = 1, ndim
         itno = 0
c         do n = 1, nucpg-1
         do n = 1, ndim-1
            if( lz(nscr(n)) .gt. lz(nscr(n+1)) )then
c..switch n and n+1 in nscr (index array)
               iscr = nscr(n)
               nscr(n) = nscr(n+1)
               nscr(n+1)   = iscr
               itno = itno + 1
            endif
         enddo
         if( itno .le. 0 )then
            if(irank.eq.0) write(*,'(20a6)')(cnuc(nscr(i)),i=1, nucpg)
            goto 130
         endif
      enddo
 130  continue
      jz = j-1
      if(irank.eq.0) write(*,*)'reordered in Z in ',jz,' steps'

c cie: ditto, see ~20l above      
      do j = 1, ndim
         itno = 0
         do n = 1, ndim-1
            if(  lz(nscr(n))   + ln(nscr(n)) .gt. 
     1           lz(nscr(n+1)) + ln(nscr(n+1)) .and.
     2           lz(nscr(n)) .eq. lz(nscr(n+1)))then
c..switch n and n+1 in nscr (index array)
               iscr = nscr(n)
               nscr(n) = nscr(n+1)
               nscr(n+1)   = iscr
               itno = itno + 1
            endif
         enddo

         if( itno .le. 0 )then
            if(irank.eq.0) write(*,'(20a6)')(cnuc(nscr(i)),i=1, nucpg)
            goto 140
         endif
      enddo
 140  continue
      ja = j-1
      if(irank.eq.0) write(*,*)'reordered in A in ',ja,' steps'

      do n = 1, nucpg
         nucp(n) = nscr(n)
      enddo
      if(irank.eq.0) write(*,*)'replaced original by reordered index'
c..   reordering finished

      if( jz .ne. 0 .or. ja .ne. 0 )then
c..   redefine net.rc values
         open(lun,file=trim(netrcfile))
         do j = 1, netsize
            read(lun,'(3i5,a5,0pf10.4,1pe12.4)')
     1           i,nz(j),nn(j),xid(j),qex(j),solarx(j)
         enddo
c..   overwrite array, including zeros, for ease in reading
c..   by many routines
         do j = 1, nucpg, 10
            write(lun,'(10i5)')(nucp(i),i=j,j+9)
         enddo
         close(lun)
         if(irank.eq.0) write(*,*) trim(netrcfile),
     1      ' adjusted for new nucp index array'
      endif

      if(irank.eq.0) write(*,*)'LEAVING ABINIT'
      return
      end



