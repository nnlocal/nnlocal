      subroutine coupling2
c--- this routine calculates alpha-s using the now-determined
c--- value of nflav and makes CKM matrix diagonal if necessary
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'verbose.f'
      include 'nlooprun.f'
      include 'nflav.f'
      include 'b0.f'
      include 'dynamicscale.f'
      include 'fourthgen.f'
      include 'couple.f'
      include 'part.f'
      include 'order.f'
      integer nproc
      double precision alphas,cmass,bmass
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
      common/nproc/nproc

c--- set up the beta-function
      b0=(xn*11d0-2d0*nflav)/6d0
      write(6,*)'*******************************************'
      write(6,*)'*     No light quarks: nflav set to 0     *'
      write(6,*)'*******************************************'
      b0=(xn*11d0)/6d0

c--- initialize the pdf set
      nlooprun=0
      call pdfwrap      

      musq=scale**2
c--- set up masses used in running of alphas (alfamz.f)
      cmass=dsqrt(mcsq)
      if (fourthgen) then
        continue ! normal b mass already set in chooser.f
      else
        bmass=dsqrt(mbsq) ! use the mass specified in the input file
      endif
      
      if (nflav .lt. 5) then
        bmass=1001d0
      endif
      if (nflav .lt. 4) then
        cmass=1000d0
      endif
 
c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfwrap already
      if (order .eq. 0) then
         nlooprun=1
         if (nproc.eq.709) nlooprun=3
      elseif (abs(order) .eq. 1) then
         nlooprun=2
         if (nproc.eq.709) nlooprun=3
      elseif (abs(order) .eq. 2) then
         nlooprun=3
      endif

c--- initialize alpha_s
      as=alphas(abs(scale),amz,nlooprun)
      
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

      if (verbose) then
      write(6,*)
      write(6,*) '***************** CKM mixing matrix ****************'
      write(6,*) '*                                                  *'
      write(6,47) Vud,Vus,Vub
      write(6,48) Vcd,Vcs,Vcb
      write(6,*) '****************************************************'
 47   format(' *      Vud=',g10.5,'Vus=',g10.5,'Vub=',g10.5,'  *')
 48   format(' *      Vcd=',g10.5,'Vcs=',g10.5,'Vcb=',g10.5,'  *')
      endif      

      if ((verbose) .and. (scale .gt. 0d0)) then      
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      if (dynamicscale .eqv. .false.) then
      write(6,49) 'alpha_s (scale)',gsq/fourpi
      write(6,49) 'alpha_s (zmass)',amz
      else
      write(6,*) '*  Dynamic scale - alpha_s changed event-by-event  *'
      write(6,49) 'alpha_s (zmass)',amz
      endif
      write(6,50) ' (using ',nlooprun,'-loop running of alpha_s)'  
      write(6,*) '****************************************************'
 49   format(' *  ',a20,f12.8,16x,'*')
 50   format(' *  ',6x,a8,i1,a25,8x,'*')
      endif
      
      return
      end
