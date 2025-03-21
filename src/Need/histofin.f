      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/nplot*'lin'/
      end

      subroutine histofin(xsec,xsec_err,itno,itmx)
c--- This outputs the final histograms for itno=0
c--- For itno>0, this is an intermediate result only
      implicit none
      include 'verbose.f'
      include 'PDFerrors.f'
      include 'histo.f'
      include 'ehisto.f'
      include 'outputoptions.f'
      include 'irregbins_incl.f'
      integer j,nlength,itno,itmx,nplotmax,nempty
      character*255 outfiledat,outfiletop,outfileerr
      character*9 runname
c--F  Add gnuplot output and root output
      character*255, outfilegnuplot, outfileps
      character*255, outfileroot, outfilerootC
c--F
c      character*255 outfilepwg
      character*3 oldbook
      character mop
      double precision xsec,xsec_err,scalefac,itscale
      logical scaleplots                  
      integer itmx1,ncall1,itmx2,ncall2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/runname/runname
      common/nlength/nlength
      common/nplotmax/nplotmax
      common/scaleplots/scalefac,scaleplots

      if (itno .eq. 0) then
      write(6,*)
      write(6,*) '****************************************************'
      write(6,*) 'output files  ',runname(1:nlength)
      write(6,*) '****************************************************'
      call flush(6)
      scaleplots=.false.
      else
      scaleplots=.true.
      scalefac=1d0/dfloat(itno)
      endif

      outfiledat=runname
      outfiletop=runname
      outfilegnuplot=runname
      outfileps=runname
      outfileroot=runname
      outfilerootC=runname
      outfileerr=runname
c      outfilepwg=runname
      outfiledat(nlength+1:nlength+4)='.dat'
      outfiletop(nlength+1:nlength+4)='.top'
      outfilegnuplot(nlength+1:nlength+4)='.gnu'
      outfileps(nlength+1:nlength+3)='.ps'
      outfileroot(nlength+1:nlength+5)='.root'
      outfilerootC(nlength+1:nlength+2)='.C'
      outfileerr(nlength+1:nlength+10)='_error.top'
c      outfilepwg(nlength+1:nlength+7)='pwg.top'
      

      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        open(unit=95,file=outfileerr,status='unknown')
      endif

      if (writedat) then
        open(unit=98,file=outfiledat,status='unknown')
      endif
      if (writetop) then
      open(unit=99,file=outfiletop,status='unknown')
      endif
      if (writegnu) then
      open(unit=97, file=outfilegnuplot,status='unknown')
      endif
      if (writeroot) then
      open(unit=96, file=outfilerootC, status='unknown')
      endif
c      if (writepwg) then
c      open(unit=100,file=outfilepwg,status='unknown')
c      endif
      
c--- write out run info to top of files
      if (writedat) then
      call writeinfo(98,' (',xsec,xsec_err,itno)      
      endif
      if (writetop) then
      call writeinfo(99,' (',xsec,xsec_err,itno)      
      endif
c      if (writegnu) then
c      call writeinfo(97,'# ',xsec,xsec_err,itno) 
c      write(97,120) outfileps(1:nlength+3)
c      endif

  120 FORMAT (/1x,
     & ' set terminal postscript col enhanced', /1x,
     & ' set output "', A, '"', /1X,
     & ' set style data points', /1X,
     & ' set key off')


c--- make sure to scale results by the maximum number of iterations
      if (itno .eq. 0) then
        itscale=1d0/dfloat(itmx)
        mop='V'
      else
        itscale=1d0/dfloat(itno)
        mop='U'
      endif
      
c--- calculate the errors in each plot (and store in 2*maxhisto+j)    
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Calculating errors for plot ',j
        call flush(6)
      endif
      call mopera(j,mop,maxhisto+j,2*maxhisto+j,itscale,1d0)
      enddo

      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Finalizing plot ',j
        call flush(6)
      endif
c--- ensure that MFINAL doesn't turn off booking for intermediate results
      oldbook=book(j)
      call mfinal(j)
      if (itno .gt. 0) then
      book(j)=oldbook
      endif
      enddo

c--- Perform integrals on plots if required (at end of run only)
      if (itno .eq. 0) then
        do j=1,nplotmax
          if (index(title(j),'+INTEGRAL+') .gt. 0) then
          call integratehisto(j)
          endif
        enddo
      endif
      
      nempty=0
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Writing .dat for plot ',j
        call flush(6)
      endif
c      if (writedat) then
c        call mprint(j,2*maxhisto+j)
c      endif
c      if (book(j) .ne. 'YES') nempty=nempty+1
      enddo
      
      if (writedat) close(unit=98)

c---generate topdrawer file - only for non-empty plots
c--F also gnuplot files
      do j=1,nplotmax-nempty
      if (verbose) then
c        write(6,*) 'Writing .top for plot ',j
        call flush(6)
      endif
c      if (writetop) then
c        call mtop(j,2*maxhisto+j,'x','y',linlog(j))
c        if (irregbin(j)) call getirregbins(j)
c      endif
      if (writegnu) then
        call mgnuplot(j,2*maxhisto+j,'x','y',linlog(j))
      endif
c      if (writeroot) then
c        call mrootplot(j,2*maxhisto+j,'x','y')
c      endif
      if ((PDFerrors) .and. (IHISTOMATCH(j) .ne. 0)) then
        call emtop(j,2*maxhisto,'x','y',linlog(j))
      endif
      enddo
      if (writetop) close(unit=99)
      if (writegnu) close(unit=97)


c---generate error file
      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        do j=1,nplotmax
          if (IHISTOMATCH(j) .ne. 0) then
            if (verbose) then
c              write(6,*) 'Writing .top for plot ',j
              call flush(6)
            endif
            call etop(j,2*maxhisto,'x','y',linlog(j))
          endif
        enddo
        close (unit=95)
      endif
      
      return
      end

