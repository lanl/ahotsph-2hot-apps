      subroutine build(irank,idbug,filenm_orig)

c..reads Thielemann file format
c..identifies symbol in Z,N
c..and echo to standard i/o

      implicit none

      include 'dimen'

      include 'crate'
      include 'comcsolve'

      include 'cburn'
      include 'ctmstp'

      integer*4 irank,idbug
      integer*4 idebug,n,i,j
      character*20 filenm_orig
      character*20 filenm

c---------------------------------------------------------------
c..get friedel's reaction rate parameters, build network inferred
c..from reaction rates available
      filenm = trim(filenm_orig)   
      idebug = 0
      call getfkt(idebug,irank)

      if(irank.eq.0) write(*,*)'return from getfkt'
      
ccccccccccccccccccccccc
      itot = nnuc
      
c..actual size of network (number of nuclei)
      netsize = itot
c      filenm = 'net.rc.1'
      call abinit(irank,filenm_orig)
      
c..setup binding energies and 2J+1 factors
      call masses(qq,w,itot,irank)
c..   save mass excesses for cburn
      do n = 1, ndim
         qex(n) = qq(n)
         gspin(n) = w(n)
      enddo
      
c..   initialize reaction rate array
      do i = 1, ireac
         sig(i) = 0.0d0
      enddo

c..   setup identifier integers for reactants
      call naray(idebug,irank)

c..   find inverse reactions
      call inverses(irank)

c     .weak interaction data
      call weakread(irank)

c..set up identification for nuclei to use in tycho timestep control
      do j = 1, itot
c..define nucleon number  for each isotope 
c        xa(j) = nz(j) + nn(j)
c..   define atomic weight for each isotope from mass excesses
        xa(j) = nz(j) + nn(j) + qex(j)/931.487d0

        if( nz(j) .eq. 1 .and. nn(j) .eq. 0 )then
c..h1
          nxdt(1) = j
        elseif( nz(j) .eq. 2 .and. nn(j) .eq. 1 )then
c..he3
          nxdt(2) = j
        elseif( nz(j) .eq. 2 .and. nn(j) .eq. 2 )then
c..he4
          nxdt(3) = j
        elseif( nz(j) .eq. 6 .and. nn(j) .eq. 6 )then
c..c12
          nxdt(4) = j
        elseif( nz(j) .eq. 7 .and. nn(j) .eq. 7 )then
c..n14
          nxdt(5) = j
        elseif( nz(j) .eq. 8 .and. nn(j) .eq. 8 )then
c..o16
          nxdt(6) = j
        elseif( nz(j) .eq. 10 .and. nn(j) .eq. 10 )then
c..ne20
          nxdt(7) = j
        elseif( nz(j) .eq. 12 .and. nn(j) .eq. 12 )then
c..mg24
          nxdt(8) = j
        elseif( nz(j) .eq. 14 .and. nn(j) .eq. 14 )then
c..si28
          nxdt(9) = j
        endif
      enddo

      nnxdt = 9

      if(irank.eq.0) write(*,*)'BUILD: ending, called ',
     1     'GETFKT,MASSES,NARAY,INVERSES'

      return
      end


      subroutine inverses(irank)

c..   written 11/17/04 wda
c..   find reaction identity of inverse rates
c..   find multiple reaction lines per reaction link
 
      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'
      include 'ceoset'
      include 'cburn'

      integer*4 irank

      integer*4 i, j, k
      real*8 i1, i2, i3, i4, i5, i6
c..   irev  is the index of the reverse rate
c..   iline is the ordering number for subrates of a given reaction
c..   (this includes nonresonant and resonant rates)
c..   isum  is a dummy variable for counting
      integer*4 isum

c---------------------------------------------------------------
c..set up A=Z+N in double precision
      do i = 1, itot
         anuc(i) = dble( nz(i) + nn(i) )
      enddo
c..find multiple lines
      do k = 1, nreac
         iline(k) = 0
      enddo
      if(irank.eq.0) write(*,*)
     1     'INVERSES: finding rates with more than one reaction line'

c..deck 1: 1-->1
      iline(k1deck(1)) = 1
      do k = k1deck(1)+1,k2deck(1)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(1),k2deck(1)
         isum = isum + iline(k)-1
c         write(*,'(i5,2a6,2i5,a5,f8.3)')k,rname(1,k),rname(2,k),
c     1           iline(k),iffn(k),rlkh(k),qval(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 1'

c..deck 2: 1-->2
      iline(k1deck(2)) = 1
      do k = k1deck(2)+1,k2deck(2)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(2),k2deck(2)
         isum = isum + iline(k)-1
c         write(*,'(i5,3a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 2'
c..deck 3: 1-->3
      iline(k1deck(3)) = 1
      do k = k1deck(3)+1,k2deck(3)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) .and.
     3        rname(4,k) .eq. rname(4,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
c..new rate
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(3),k2deck(3)
         isum = isum + iline(k)-1
c         write(*,'(i5,4a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        rname(4,k),iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 3'

c..deck 4: 2-->1
      iline(k1deck(4)) = 1
      do k = k1deck(4)+1,k2deck(4)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
c..new rate
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(4),k2deck(4)
         isum = isum + iline(k)-1
c         write(*,'(i5,3a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 4'

c..deck 5: 2-->2
      iline(k1deck(5)) = 1
      do k = k1deck(5)+1,k2deck(5)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) .and.
     3        rname(4,k) .eq. rname(4,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
c..new rate
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(5),k2deck(5)
         isum = isum + iline(k)-1
c         write(*,'(i5,4a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        rname(4,k),iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 5'


c..deck 6: 2-->3
      iline(k1deck(6)) = 1
      do k = k1deck(6)+1,k2deck(6)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) .and.
     3        rname(4,k) .eq. rname(4,k-1) .and.
     1        rname(5,k) .eq. rname(5,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
c..new rate
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(6),k2deck(6)
         isum = isum + iline(k)-1
c         write(*,'(i5,5a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        rname(4,k),rname(5,k),iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 6'

c..deck 7: 2-->4
      iline(k1deck(7)) = 1
      do k = k1deck(7)+1,k2deck(7)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) .and.
     3        rname(4,k) .eq. rname(4,k-1) .and.
     4        rname(5,k) .eq. rname(5,k-1) .and.
     5        rname(6,k) .eq. rname(6,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
c..new rate
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(7),k2deck(7)
         isum = isum + iline(k)-1
c         write(*,'(i5,6a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        rname(4,k),rname(5,k),rname(6,k),iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 7'

c..deck 8: 3-->1 and 2
      iline(k1deck(8)) = 1
      do k = k1deck(8)+1,k2deck(8)
         if( rname(1,k) .eq. rname(1,k-1) .and.
     1        rname(2,k) .eq. rname(2,k-1) .and.
     2        rname(3,k) .eq. rname(3,k-1) .and.
     3        rname(4,k) .eq. rname(4,k-1) .and.
     4        rname(5,k) .eq. rname(5,k-1) )then
c..repeat
            iline(k) = iline(k-1)+1
         else
c..new rate
            iline(k) = 1
         endif
      enddo
      isum = 0
      do k = k1deck(8),k2deck(8)
         isum = isum + iline(k)-1
c         write(*,'(i5,5a6,i5)')k,rname(1,k),rname(2,k),rname(3,k),
c     1        rname(4,k),rname(5,k),iline(k)
      enddo
      if(irank.eq.0) write(*,*)isum,' multiline rates in deck 8'

c.................................................................
c..find inverse rates
      do k = 1, nreac
         irev(k) = 0
      enddo


c deck 1: decays, i ---> j
      do k= k1deck(1),k2deck(1)
         i1 = nrr(1,k)
         i2 = nrr(2,k)
c         write(*,'(2i5,2a6,5x,a6,2a2,f8.3)')k,ideck(k),
c     1        (rname(i,k),i=1,2),rlkh(k),rvw(k),rnr(k),qval(k)
c..   look for reverse rate
c..   ffnu convention is higher Z is first for a given A
c..   "forward" rate has a positive Q-value
         do j = k1deck(1), k2deck(1)
            if(  rname(1,k) .eq. rname(2,j) .and.
     1           rname(2,k) .eq. rname(1,j) .and.
     2           iline(k)   .eq. iline(j))then
              if(  qval(k) .gt. 0.0d0 )then
c..count only once 
                  irev(k) = j
                  irev(j) = -k
              else
                  irev(k) = -j
                  irev(j) = k
              endif
c                  write(*,'(2(2i5,2a6,3x,a6,2a2,3i5,a5),i5)')
c     1                 k,ideck(k),(rname(i,k),i=1,2),rlkh(k),
c     2                 rvw(k),rnr(k),irev(k),iline(k),iffn(k),'*',
c     3                 j,ideck(j),(rname(i,j),i=1,2),rlkh(j),
c     4                 rvw(j),rnr(j),irev(j),iline(j),iffn(j),'*',
c     5                 nz(i1)-nz(i2)

            endif
         enddo
      enddo

c deck 2: dissociations, i ---> j + k; deck 4: j + k ---> i
      do k= k1deck(2),k2deck(2)
         i1 = nrr(1,k)
         i2 = nrr(2,k)
         i3 = nrr(3,k)
c            write(*,'(2i5,7a6,5x,a6,2a2)')k,ideck(k),(rname(i,k),i=1,3),
c     1           rlkh(k),rvw(k),rnr(k)
c..look for its forward rate
         do j = k1deck(4),k2deck(4)
            if( rname(1,k) .eq. rname(3,j) .and. 
     1           rname(2,k) .eq. rname(1,j) .and.
     2           rname(3,k) .eq. rname(2,j) )then
               if( iline(k) .eq. iline(j) )then
                  irev(j) = k
                  irev(k) = -j
c                  write(*,'(10x,2i5,3a6,5x,a6,2a2,4i5)')j,ideck(j),
c     1              (rname(i,j),i=1,3),rlkh(j),rvw(j),rnr(j),
c     2              j,irev(j),iline(j),iline(k)
               endif
               
            endif
         enddo
      enddo

c deck 5: exchanges, i + j ---> k + l
      do k= k1deck(5),k2deck(5)
         i1 = nrr(1,k)
         i2 = nrr(2,k)
         i3 = nrr(3,k)
         i4 = nrr(4,k)
c         write(*,'(2i5,4a6,5x,a6,2a2)')k,ideck(k),(rname(i,k),i=1,4),
c     1        rlkh(k),rvw(k),rnr(k)

c..look for its forward rate
         do j = k1deck(5),k2deck(5)
            if(  rname(1,k) .eq. rname(3,j) .and. 
     1           rname(2,k) .eq. rname(4,j) .and.
     2           rname(3,k) .eq. rname(1,j) .and.
     3           rname(4,k) .eq. rname(2,j) )then
               if( iline(k) .eq. iline(j) )then
                  if( qval(k) .gt. 0.0d0 )then
                     irev(j) = -k
                     irev(k) = j
                  else
                     irev(j) = k
                     irev(k) = -j
                  endif
c                  write(*,'(10x,2i5,4a6,5x,a6,2a2,4i5)')j,ideck(j),
c     1                 (rname(i,j),i=1,4),rlkh(j),rvw(j),rnr(j),
c     2                 j,irev(j),iline(j),iline(k)
               endif
            endif
         enddo
      enddo

c deck 3 and 8: fission and fusion 1--->3 and 3-->1 
      do k= k1deck(3),k2deck(3)
         i1 = nrr(1,k)
         i2 = nrr(2,k)
         i3 = nrr(3,k)
         i4 = nrr(4,k)
c           write(*,'(2i5,4a6,5x,a6,2a2)')k,ideck(k),(rname(i,k),i=1,4),
c     1           rlkh(k),rvw(k),rnr(k)

c..look for its forward rate
         do j = k1deck(8),k2deck(8)
            if(  rname(1,k) .eq. rname(4,j) .and. 
     1           rname(2,k) .eq. rname(1,j) .and.
     2           rname(3,k) .eq. rname(2,j) .and.
     3           rname(4,k) .eq. rname(3,j) )then
               if( iline(k) .eq. iline(j) )then
                  if( qval(k) .gt. 0.0d0 )then
                     irev(j) = -k
                     irev(k) = j
                  else
                     irev(j) = k
                     irev(k) = -j
                  endif
c                     write(*,'(10x,2i5,4a6,5x,a6,2a2,4i5)')j,ideck(j),
c     1                    (rname(i,j),i=1,4),rlkh(j),rvw(j),rnr(j),
c     2                    j,irev(j),iline(j),iline(k)
               endif
            endif
         enddo
      enddo

c deck 6 and 8: fission and fusion 2--->3 and 3-->2
      do k= k1deck(6),k2deck(6)
         i1 = nrr(1,k)
         i2 = nrr(2,k)
         i3 = nrr(3,k)
         i4 = nrr(4,k)
         i5 = nrr(5,k)
c            write(*,'(2i5,5a6,5x,a6,2a2)')k,ideck(k),(rname(i,k),i=1,5),
c     1           rlkh(k),rvw(k),rnr(k)

c..look for its forward rate
         do j = k1deck(8),k2deck(8)
            if(  rname(1,k) .eq. rname(4,j) .and. 
     1           rname(2,k) .eq. rname(5,j) .and.
     2           rname(3,k) .eq. rname(1,j) .and.
     3           rname(4,k) .eq. rname(2,j) .and.
     4           rname(5,k) .eq. rname(3,j)  )then
               if( iline(k) .eq. iline(j) )then
                  if( qval(k) .gt. 0.0d0 )then
                     irev(j) = -k
                     irev(k) = j
                  else
                     irev(j) = k
                     irev(k) = -j
                  endif
c                     write(*,'(5x,2i5,5a6,5x,a6,2a2,4i5)')j,
c     1                    ideck(j),
c     1                    (rname(i,j),i=1,5),rlkh(j),rvw(j),rnr(j),
c     2                    j,irev(j),iline(j),iline(k)
               endif
            endif
         enddo
      enddo

      if(irank.eq.0) write(*,*)
     1   'These rates do not have inverses in netsu or netweak'
      if(irank.eq.0) write(*,'(10x,2a5,6a6,5x,a6,2a2,2a5)')
     1   'j','deck','i1','i2',
     1     'i3','i4','i5','i6','rlkh ','vw','rn','irev',
     2     'iffn'
c..   seek unmatched rates
      do j  = k1deck(1),k2deck(8)
         if( irev(j) .eq. 0 )then
         if(irank.eq.0) write(*,'(10x,2i5,6a6,5x,a6,2a2,2i5)')
     1      j,ideck(j),
     1        (rname(i,j),i=1,6),rlkh(j),rvw(j),rnr(j),
     2        irev(j),iffn(j)
         endif
      enddo
      if(irank.eq.0) write(*,'(10x,2a5,6a6,5x,a6,2a2,2a5)')
     1   'j','deck','i1','i2',
     1     'i3','i4','i5','i6','rlkh ','vw','rn','irev',
     2     'iffn'
      i = 0
      do j = 1, nreac
         if( irev(j) .ne. 0 )i=i+1
      enddo
      if(irank.eq.0) write(*,*)i,' reactions have inverses'
      if(irank.eq.0) write(*,*)'leaving ',k2deck(8)-i,
     1     ' extra lone rates listed above, which do not'
      

      return
      end
