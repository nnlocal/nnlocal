      subroutine setnotag()
      implicit none
      include 'removebr.f'
      integer nproc,notag
      common/nproc/nproc
      common/notag/notag
c--- this routine sets the value of "notag", the number of jets
c--- that may be safely ignored without affecting finiteness of result;
c--- the minimum number of jets allowed by the code is equal to
c---    nqcdjets - notag
c--- where nqcdjets is itself process-dependent

c--- Modifying this routine to allowe *larger* values of notag
c--- than the defaults below (or >0 for processes not listed here)
c--- should be done with care

c--- Modifying this routine to allow *smaller* values of notag,
c--- i.e. stricter constraints on the number of jets observed,
c--- should not cause problems

      if (nproc .eq. 709) then
c---- H+J production
        notag=0
      endif
      
      return
      end
      
