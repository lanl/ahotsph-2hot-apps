      subroutine naray(idebug,irank)

      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'

      integer*4 irank

      integer*4 idebug, i, j, k, l, k1, k2, inv, iflag
      character*5 blank

      real*8 sparse_tmp(nnuc*nnuc*10)
      integer*4 iloc_tmp(nnuc*nnuc*10),jloc_tmp(nnuc*nnuc*10)
c..set up interger identification of reactants

c..1-->1 reaction (deck=1)
c..1-->2 reaction (deck=2)
c..1-->3 reaction (deck=3)
c..2-->1 reaction (deck=4)
c..2-->2 reaction (deck=5)
c..2-->3 reaction (deck=6)
c..2-->4 reaction (deck=7)
c..3-->1 (and 2) reaction (deck=8)

      data blank/'     '/
c--------------------------------------------------------------
c..zero nrr array and write nonzero indexes for real reactants
      do k = 1, ireac
        do j = 1, 6
          nrr(j,k) = 0
        enddo
      enddo

c      idebug = 1
      
      do j = 1, nreac
          iloc(j) = 0
          iloc_tmp(j) = 0
          jloc(j) = 0
          jloc_tmp(j) = 0
          sparse_dfdy(j) = 0.0d0
          sparse_tmp(j) = 0.0d0
      enddo
      nvecs = 0
   
c..identify reactants by index in nz,nn,xid vectors

c..beta decays, positron decays, electron captures
c      write(*,*)'deck 1'
c..1-->1 reaction (deck=1)
      k1        = 1
      k1deck(1) = k1
      k2        = ndeck(1)
      k2deck(1) =  k2
      do k = k1,k2
        do j = 1,2
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,2
          do i = 1,2
             nvecs = nvecs+1
             iloc_tmp(nvecs) = nrr(j,k)
             jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo

      if( idebug .ne. 0 )then
      do k = k1,k2
       if(irank.eq.0)  write(*,'(i5,a5,2x,a5,2x,a5,5x,a4,2a1)') k,
     1   xid(nrr(1,k))," --> ",
     1   xid(nrr(2,k)),rlkh(k),rnr(k),rvw(k)
      enddo
      endif

c      write(*,*)'deck 2'
c..1-->2 reaction (deck=2)
      k1 =  ndeck(1)+1
      k2 =  ndeck(1) + ndeck(2)
      k1deck(2) = k1
      k2deck(2) = k2
      do k = k1, k2
        do j = 1,3
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            nvecs = nvecs+1
            iloc_tmp(nvecs) = nrr(j,k)
            jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo

      if( idebug .ne. 0 )then
      do k = k1,k2
        if(irank.eq.0) write(*,'(i5,4(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k))," --> ",
     1   xid(nrr(2,k)),xid(nrr(3,k)),rlkh(k),rnr(k),rvw(k)
      enddo
      endif


c      write(*,*)'deck 3'
c..1-->3 reaction (deck=3)
      k1 = ndeck(1) + ndeck(2) + 1
      k2 = k1 - 1 + ndeck(3)
      k1deck(3) = k1
      k2deck(3) = k2
      do k = k1, k2
        do j = 1,4
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,4
          do i = 1,4
            nvecs = nvecs+1
            iloc_tmp(nvecs) = nrr(j,k)
            jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo

      if( idebug .ne. 0 )then
      do k = k1,k2
        if(irank.eq.0) write(*,'(i5,5(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k))," --> ",
     1   xid(nrr(2,k)),xid(nrr(3,k)),xid(nrr(4,k)),
     1   rlkh(k),rnr(k),rvw(k)
      enddo
      endif


c      write(*,*)'deck 4'
c..2-->1 reaction (deck=4)
      k1 = ndeck(1) + ndeck(2) + ndeck(3) + 1
      k2 = k1 - 1 + ndeck(4)
      k1deck(4) = k1
      k2deck(4) = k2
      do k = k1, k2
        do j = 1,3
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            nvecs = nvecs+1
            iloc_tmp(nvecs) = nrr(j,k)
            jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo

      if( idebug .ne. 0 )then
      do k = k1,k2
        if(irank.eq.0) write(*,'(i5,4(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k)),xid(nrr(2,k))," --> ",
     1   xid(nrr(3,k)),
     1   rlkh(k),rnr(k),rvw(k)
      enddo
      endif
c
c
c      write(*,*)'deck 5'
c..2-->2 reaction (deck=5)
      k1 = ndeck(1) + ndeck(2) + ndeck(3) + ndeck(4) + 1
      k2 = k1 - 1 + ndeck(5)
      k1deck(5) = k1
      k2deck(5) = k2
      do k = k1, k2
        do j = 1,4
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,4
          do i = 1,4
            nvecs = nvecs+1
            iloc_tmp(nvecs) = nrr(j,k)
            jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo
c
      if( idebug .ne. 0 )then
      do k = k1,k2
        if(irank.eq.0) write(*,'(i5,5(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k)),xid(nrr(2,k))," --> ",
     1   xid(nrr(3,k)),xid(nrr(4,k)),
     1   rlkh(k),rnr(k),rvw(k)
      enddo
      endif

      if( ndeck(6) .ne. 0 )then
c..2-->3 reaction (deck=6)
      k1 = ndeck(1) + ndeck(2) + ndeck(3) + ndeck(4) + ndeck(5) + 1
      k2 = k1 - 1 + ndeck(6)
      k1deck(6) = k1
      k2deck(6) = k2
      do k = k1, k2
        do j = 1,5
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,5
          do i = 1,5
            nvecs = nvecs+1
            iloc_tmp(nvecs) = nrr(j,k)
            jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo

        if( idebug .ne. 0 )then
        do k = k1,k2
          if(irank.eq.0) write(*,'(i5,6(a5,2x),3x,a4,2a1)') k,
     1     xid(nrr(1,k)),xid(nrr(2,k))," --> ",
     1     xid(nrr(3,k)),xid(nrr(4,k)),xid(nrr(5,k)),
     1     rlkh(k),rnr(k),rvw(k)
        enddo
        endif

      endif

      if( ndeck(7) .ne. 0 )then
c..2-->4 reaction (deck=7)
      k1 = ndeck(1) + ndeck(2) + ndeck(3) + ndeck(4) + ndeck(5)
     1   + ndeck(6) + 1
      k2 = k1 - 1 + ndeck(7)
      k1deck(7) = k1
      k2deck(7) = k2
      do k = k1, k2
        do j = 1,6
          do i = 1, itot
            if( rname(j,k) .eq. xid(i) )then
              nrr(j,k) = i
            endif
          enddo
        enddo
        do j = 1,6
          do i = 1,6
            nvecs = nvecs+1
            iloc_tmp(nvecs) = nrr(j,k)
            jloc_tmp(nvecs) = nrr(i,k)
          enddo
        enddo
      enddo

      if( idebug .ne. 0 )then
      do k = k1,k2
        if(irank.eq.0) write(*,'(i5,7(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k)),xid(nrr(2,k))," --> ",
     1   xid(nrr(3,k)),xid(nrr(4,k)),xid(nrr(5,k)),xid(nrr(6,k)),
     1   rlkh(k),rnr(k),rvw(k)
      enddo
      endif

      endif

c      write(*,*)'deck 8'
c..3-->1 (and 2) reaction (deck=8)
      k1 = ndeck(1) + ndeck(2) + ndeck(3) + ndeck(4) + ndeck(5)
     1   + ndeck(6) + ndeck(7)           + 1
      k2 = k1 - 1 + ndeck(8)
      k1deck(8) = k1
      k2deck(8) = k2
      do k = k1, k2
        if( rname(5,k) .eq. blank )then
          do j = 1,4
            do i = 1, itot
              if( rname(j,k) .eq. xid(i) )then
                nrr(j,k) = i
              endif
            enddo
          enddo
          do j = 1,4
            do i = 1,4
              nvecs = nvecs+1
              iloc_tmp(nvecs) = nrr(j,k)
              jloc_tmp(nvecs) = nrr(i,k)
            enddo
          enddo
        else
          do j = 1,5
            do i = 1, itot
              if( rname(j,k) .eq. xid(i) )then
                nrr(j,k) = i
              endif
            enddo
          enddo
          do j = 1,5
            do i = 1,5
              nvecs = nvecs+1
              iloc_tmp(nvecs) = nrr(j,k)
              jloc_tmp(nvecs) = nrr(i,k)
            enddo
          enddo
        endif
      enddo
      
      nlinks = 0
      do j = 1,nvecs
         sparse_tmp(j) = 1.0d0
      enddo
      do j=1,nvecs
         if(sparse_tmp(j) .ne. 0.0d0)then
            do i=j+1,nvecs
               if((iloc_tmp(j) .eq. iloc_tmp(i)) .and.
     1             (jloc_tmp(j) .eq. jloc_tmp(i)))then
                  sparse_tmp(i) = 0.0d0
               endif
            enddo
         endif
      enddo
      

      do j=1,nvecs
         if(sparse_tmp(j) .ne. 0.0d0)then
            nlinks = nlinks+1
            sparse_dfdy(nlinks) = 1.0d0
            iloc(nlinks) = iloc_tmp(j)
            jloc(nlinks) = jloc_tmp(j)
            ivect(nlinks) = iloc(nlinks)
            jvect(nlinks) = jloc(nlinks)
            if(iloc(nlinks) .ne. jloc(nlinks))then
               sparse_dfdy(nlinks) = 1.0d-10
            endif
         endif
      enddo

      sparseu = 0.1d0

      call ma28ad(nnuc,nlinks,sparse_dfdy,nreac,ivect,nreac,jvect,
     1               sparseu,ikeep,iw,sparsew,iflag)
      
      if( idebug .ne. 0 )then
      do k = k1,k2
        if( nrr(5,k) .gt. 0 )then
        if(irank.eq.0) write(*,'(i5,6(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k)),xid(nrr(2,k)),xid(nrr(3,k))," --> ",
     1   xid(nrr(4,k)),xid(nrr(5,k)),
     1   rlkh(k),rnr(k),rvw(k)
        else
        if(irank.eq.0) write(*,'(i5,5(a5,2x),3x,a4,2a1)') k,
     1   xid(nrr(1,k)),xid(nrr(2,k)),xid(nrr(3,k))," --> ",
     1   xid(nrr(4,k)),
     1   rlkh(k),rnr(k),rvw(k)
        endif
      enddo
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      inv = 0
      if( inv .eq. 0 )return
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..identify inverse rates
      if(irank.eq.0) write(*,*)'deck = 1'
      inv = 0
      do k = k1deck(1), k2deck(1)
         do i = k1deck(1), k2deck(1)
            if( nrr(1,k) .eq. nrr(2,i) .and. 
     1           nrr(2,k) .eq. nrr(1,i) )then
               if(irank.eq.0) write(*,'(i5,5x,6i5,5x,i5,5x,6i5)') 
     1            k,(nrr(j,k),j=1,6),
     1             i, (nrr(j,i), j=1,6 )
               inv = inv + 1
            endif
         enddo
      enddo
      if(irank.eq.0) write(*,*)inv,' inverses found in deck 1'

      if(irank.eq.0) write(*,*)'deck = 2 and deck = 4'

      inv = 0
      do k = k1deck(2), k2deck(2)
c..multiple resonances for F18 <--> O17 + p
c       do k = 139, 143
         do i = k1deck(4), k2deck(4)

            if(  nrr(1,k) .eq. nrr(3,i) .and.  nrr(2,k) .eq. nrr(1,i) 
     1           .and. nrr(3,k) .eq. nrr(2,i) )then
               
               if( rnr(k) .eq. rnr(i) )then
                  if( rcoef(5,k) .eq. rcoef(5,i) )then
                     if(irank.eq.0) write(*,'(1p7e12.4)')
     1                  (rcoef(j,k),j=1,7)
                     if(irank.eq.0) write(*,'(1p7e12.4)')
     1                  (rcoef(j,i),j=1,7)
                     inv = inv + 1
                     if(irank.eq.0) write(*,
     1                  '(i5,2x,2(i5,4(a5,2x),3x,a4,2a1,10x))') 
     1                    inv,k,
     1                    xid(nrr(1,k))," --> ",xid(nrr(2,k)),
     1                    xid(nrr(3,k)),
     1                    rlkh(k),rnr(k),rvw(k),
     2                    i,
     2                    xid(nrr(1,i)), xid(nrr(2,i))," --> ",
     2                    xid(nrr(3,i)),rlkh(i),rnr(i),rvw(i)
                  endif
               endif
            endif

         enddo
      enddo
      if(irank.eq.0) write(*,*)inv,' inverses found in deck 2'
      if(irank.eq.0) write(*,*)'deck 2 ',k2deck(2)-k1deck(1)+1
      if(irank.eq.0) write(*,*)'deck 4 ',k2deck(4)-k1deck(4)+1


      if(irank.eq.0) write(*,*)ndeck
      if(irank.eq.0) write(*,*)k1deck
      if(irank.eq.0) write(*,*)k2deck




      stop'naray'
ccccccccccccccccccccccccc

      return
      end

