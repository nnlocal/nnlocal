      subroutine bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot) 
      implicit none
      include 'nplot.f'
      include 'part.f'
      include 'outputflags.f'
      include 'vegas_common.f'
      integer n
      character*(*) titlex
      character*3 llplot
      character*4 tag
      double precision var,wt,wt2,xmin,xmax,dx

      if     (tag .eq. 'book') then
         call mbook(n,titlex,dx,xmin,xmax)
c---  Traditional MCFM histograms
c---  also book the errors now (in maxhisto+n,2*maxhisto+n)
         call mbook(maxhisto+n,titlex,dx,xmin,xmax)
         call mbook(2*maxhisto+n,titlex,dx,xmin,xmax)
         if ( (part .eq. 'virt') .or. (part .eq. 'real')
     &        .or.(part .eq. 'tota')) then
            call mbook(3*maxhisto+n,titlex,dx,xmin,xmax)
         endif
      elseif (tag .eq. 'plot') then
c---  Traditional MCFM histograms
c---  also book the errors now (in maxhisto+n); fill temp histos for real
         if (part .eq. 'born') then
            call mfill(n,var,wt)
            call mfill(maxhisto+n,var,wt2)
         else
            call mfill(3*maxhisto+n,var,wt)
         endif
         linlog(n)=llplot
         titlearray(n)=titlex
      endif
      
      return
      end

      subroutine ebookplot(n,tag,var,wt) 
      implicit none
      include 'PDFerrors.f'
      include 'outputflags.f'
      integer n
      double precision var,wt
      character tag*4

      if (PDFerrors .eqv. .false.) return

      if (tag.eq.'book') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call ebook(n)
        endif
      elseif (tag .eq. 'plot') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call efill(n,var,wt)
        endif
      endif

      return
      end

