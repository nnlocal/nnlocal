      subroutine masscuts(p,*)
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'process.f'
      include 'cutoff.f'
      logical first
      double precision p(mxpart,4),s34,s56,s36,s45
      integer nqcdjets,nqcdstart,nproc
      common/nqcdjets/nqcdjets,nqcdstart 
      common/nproc/nproc
      data first/.true./
      save first

      if (first) then
      first=.false.

     
      write(6,*)
      write(6,*) '****************** Basic mass cuts *****************'
      write(6,*) '*                                                  *'
      write(6,99) dsqrt(wsqmin),'m34',dsqrt(wsqmax)
      if (nqcdjets .lt. 2) then
      write(6,99) dsqrt(bbsqmin),'m56',dsqrt(bbsqmax)
      else
      write(6,98) dsqrt(bbsqmin),'m(jet1,jet2)',dsqrt(bbsqmax)
      endif
      write(6,*) '****************************************************'
      endif
            
c--- only apply cuts on s34 if vectors 3 and 4 are defined
      if ((abs(p(3,4)) .gt. 1d-8) .and. (abs(p(4,4)) .gt. 1d-8)) then
        s34=+(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .      -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
c--- do not accept s34<cutoff either
        if ((s34 .lt. max(wsqmin,cutoff)).or.(s34 .gt. wsqmax)) return 1
      endif
         
c--- only apply cuts on s56 if vectors 5 and 6 are defined
      if ((abs(p(5,4)) .gt. 1d-8) .and. (abs(p(6,4)) .gt. 1d-8)) then
        s56=+(p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     .      -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2
c--- do not accept s56<cutoff either
        if ((s56 .lt. max(bbsqmin,cutoff)).or.(s56.gt.bbsqmax)) return 1
      endif
     

   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')

      return
      end

