      subroutine writeinfo(unitno,commchars,xsec,xsec_err,itno)
************************************************************************
*   Routine to write out run information to a desired unit             *
************************************************************************
      implicit none
      include 'PDFerrors.f'
      include 'process.f'
      include 'outputflags.f'
      include 'runstring.f'
      include 'finitemcorr.f'
      integer unitno,j,k,itno
      double precision xsec,xsec_err
      double precision lord_bypart(-1:1,-1:1),lordnorm,rescale
      double precision ggpart,gqpart,qgpart,qqpart,qqbpart,
     . gqbpart,qbgpart,qbqbpart,qbqpart
      
      character*2 commchars
      logical dryrun,makecuts
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,origij
      integer NPTYPE,NGROUP,NSET
      double precision sqrts
      double precision Rcut 

      common/nproc/nproc
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/origij/origij

      common/bypart/lord_bypart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart

      if (itno .gt. 0) then
c--- write warning that result is only intermediate; populate the
c--- variables in finalpart (normally done in nnlocal_exit)
      write(unitno,*) commchars//
     & ' Intermediate result for iteration',itno,')'
      lordnorm=0d0
      do j=-1,1
      do k=-1,1
        lordnorm=lordnorm+lord_bypart(j,k)
      enddo
      enddo
      ggpart=lord_bypart( 0, 0)/lordnorm
      gqpart=lord_bypart( 0,+1)/lordnorm
      gqbpart=lord_bypart( 0,-1)/lordnorm
      qgpart=lord_bypart(+1, 0)/lordnorm
      qbgpart=lord_bypart(-1, 0)/lordnorm
      qqpart=lord_bypart(+1,+1)/lordnorm
      qbqbpart=lord_bypart(-1,-1)/lordnorm
      qqbpart=lord_bypart(+1,-1)/lordnorm
      qbqpart=lord_bypart(-1,+1)/lordnorm
      endif
      write(unitno,55) commchars//
     & ' Cross-section is: ',xsec,' +/-',xsec_err,')'
      write(unitno,*)

c--- for gg->H+X processes, also write out the cross section
c---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
      if ((case(1:2) .eq. 'H_') .and. (fmcorr)) then
        call finitemtcorr(rescale)
        write(unitno,55) commchars//'Rescaled x-sec is:',
     &     xsec*rescale,' +/-',xsec_err*rescale,')'
        write(unitno,*)
      endif
     
      write(unitno,*) commchars,
     &                ' Contribution from parton sub-processes:'
      write(unitno,95)commchars,'   GG    ',ggpart*xsec,ggpart*100d0
      write(unitno,95)commchars,'   GQ    ',gqpart*xsec,gqpart*100d0
      write(unitno,95)commchars,'   GQB   ',gqbpart*xsec,gqbpart*100d0
      write(unitno,95)commchars,'   QG    ',qgpart*xsec,qgpart*100d0
      write(unitno,95)commchars,'   QBG   ',qbgpart*xsec,qbgpart*100d0
      write(unitno,95)commchars,'   QQ    ',qqpart*xsec,qqpart*100d0
      write(unitno,95)commchars,'   QBQB  ',qbqbpart*xsec,qbqbpart*100d0
      write(unitno,95)commchars,'   QQB   ',qqbpart*xsec,qqbpart*100d0
      write(unitno,95)commchars,'   QBQ   ',qbqpart*xsec,qbqpart*100d0
      write(unitno,*)

      if (PDFerrors) then
        do j=0,maxPDFsets
          write(unitno,56) j,PDFxsec(j)
        enddo
        write(unitno,*)
      endif

      if (commchars .eq. ' (') then
c--- new routine for writing out contents of input file
        call writeinput(unitno,' (',' )','WRITEALL')
      else
        call writeinput(unitno,commchars,'  ','WRITEALL')
      endif

      return

c--- 55 format
   55 format(a20,f24.10,a4,f24.10,a1)
c--- 56 character format
   56 format('( PDF error set ',i3,'  --->',f13.3,' fb  )')
c--- 95 character format
   95 format(a2,5x,a9,' |',f18.5,f8.2,'%')
c--- 96 character format      
   96 format(' (',a20,12x,'[',a,']',' )')  
c--- 97 integer format      
   97 format(' (',i20,12x,'[',a,']',' )')  
c--- 98 logical format      
   98 format(' (',L20,12x,'[',a,']',' )')  
c--- 99 floating point format
   99 format(' (',f20.4,12x,'[',a,']',' )')  
      
      end
      
      
