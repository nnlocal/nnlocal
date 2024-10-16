      double precision function realint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'phasemin.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'realwt.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'process.f'
      include 'PDFerrors.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'dipolescale.f'
      include 'decay1q2a.f'
      include 'outputoptions.f'
      include 'breit.f'
      include 'outputflags.f'
      include 'runstring.f'
      include 'qcdcouple.f'
      include 'reweight.f'

      integer nproc
      common/nproc/nproc
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec,ii,jj
      double precision A1sub, A2sub, A12sub
      double precision vector(mxdim),W,val,val2,valsum,xint,ptmp
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf),
     . dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf),
     . fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
!---  DEBUG
      double precision phrans(mxpart,4),dot,sprod
!---  END DEBUG
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqa(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision msqd(maxd,-nf:nf,-nf:nf)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4)
      double precision m3,m4,m5,R,Rbbmin
      double precision xmsq_bypart(0:maxd,-1:1,-1:1),xmsqjk,
     . lord_bypart(-1:1,-1:1),plo(mxpart,4),pswtdip
      integer sgnj,sgnk
      common/xreal/xreal,xreal2
      common/Rbbmin/Rbbmin
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      double precision ph(mxpart,4),sh(mxpart,mxpart)
      integer kk,ll
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/bypart/lord_bypart
      common/incldip/incldip
      double precision realeventp(mxpart,4)
      common/realeventp/realeventp
      data p/56*0d0/
      data first/.true./
      save first,rscalestart,fscalestart

      include 'singular.f'
      
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale

c         open(1278,file='bad_points.txt',status='NEW',action='write')
         
      endif
      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      

      W=sqrts**2
      
      if (first) then
         write(6,*)
         write(6,*) 'nmin=',nmin,',nmax=',nmax
         write(6,*)
         first=.false.
      endif
      
      call gen_real_ps(vector,p,pswt,*999)

      nvec=npart+2
      call dotem(nvec,p,s)
      
c----calculate the x's for the incoming partons from generated momenta
      xx1=-two*p(1,4)/sqrts
      xx2=-two*p(2,4)/sqrts
      
      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if (check.eqv..false.) then
      
      if ((xx1 .gt.  1d0) .or. (xx2 .gt.  1d0)
     &.or.(xx1 .lt. xmin) .or. (xx2 .lt. xmin)) then
         goto 999
      endif

c---- reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
c---- check also for mapped momenta
      call transformC2_colorful(p,ph,5,6)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)
      
      call transformC2_colorful(p,ph,1,5)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)
      
      call transformC2_colorful(p,ph,2,5)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)
      
      call transformC2_colorful(p,ph,1,6)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)
      
      call transformC2_colorful(p,ph,2,6)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)
      
      call transformS1_colorful(p,ph,5)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)
      
      call transformS1_colorful(p,ph,6)
      call dotem(npart-1+2,ph,sh)
      call smalls(sh,npart-1,*999)

      
      
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal = includedipole(0,p)
      incldip(0)  = includereal 
      
      else
         includereal = .true.
      endif

      
      if (includereal .eqv. .false.) then
        do j=-nf,nf
           do k=-nf,nf
              msq(j,k)=0d0
           enddo
        enddo
      endif
      
      do j=1,mxpart
         do k=1,4
            realeventp(j,k)=p(j,k)
         enddo
      enddo
      
      
      if (dynamicscale) then
        call scaleset(rscalestart,fscalestart,p)
      dipscale(0)=facscale
      endif
      
      
c--- Calculate the required matrix elements      
      if (case .eq. 'H_gaga') then

         if (includereal) call gg_hgagagg(p,msq)
         
            call gg_hgagag_gs_colorful(p,msqc)

            if (check) then
         msqd=msqc

         write(6,*) 'CF    dipole 1',msqd(1,0,0)
         write(6,*) 'CF    dipole 2',msqd(2,0,0)
         write(6,*) 'CF    dipole 3',msqd(3,0,0)
         write(6,*) 'CF    dipole 4',msqd(4,0,0)
         write(6,*) 'CF    dipole 5',msqd(5,0,0)
         write(6,*) 'CF    dipole 6',msqd(6,0,0)
         write(6,*) 'CF    dipole 7',msqd(7,0,0)
         write(6,*) 'CF    dipole 8',msqd(8,0,0)
         write(6,*) 'CF    dipole 9',msqd(9,0,0)
         write(6,*) 'CF    dipole 10',msqd(10,0,0)
         write(6,*) 'CF    dipole 11',msqd(11,0,0)
         write(6,*) 'CF    dipole 12',msqd(12,0,0)
         write(6,*) 'CF    dipole 13',msqd(13,0,0)
         write(6,*) 'CF    dipole 14',msqd(14,0,0)
         write(6,*) 'CF    dipole 15',msqd(15,0,0)
         write(6,*) 'CF    dipole 16',msqd(16,0,0)
         write(6,*) 'CF    dipole 17',msqd(17,0,0)
         write(6,*) 'CF    dipole 18',msqd(18,0,0)
         write(6,*) 'CF    dipole 19',msqd(19,0,0)
         write(6,*) 'CF    dipole 20',msqd(20,0,0)
         write(6,*) 'CF    dipole 21',msqd(21,0,0)
         write(6,*) 'CF    dipole 22',msqd(22,0,0)
         write(6,*) 'CF    dipole 23',msqd(23,0,0)
         write(6,*) 'CF    dipole 24',msqd(24,0,0)
         write(6,*) 'CF    dipole 25',msqd(25,0,0)
         write(6,*) 'CF    dipole 26',msqd(26,0,0)
         write(6,*) 'CF    dipole 27',msqd(27,0,0)
         write(6,*) 'CF    dipole 28',msqd(28,0,0)         
         write(6,*) 'CF    dipole 29',msqd(29,0,0)
         write(6,*) 'CF    dipole 30',msqd(30,0,0)
         write(6,*) 'CF    dipole 31',msqd(31,0,0)
         write(6,*) 'CF    dipole 32',msqd(32,0,0)
         write(6,*) 'CF    dipole 33',msqd(33,0,0)
         write(6,*) 'CF    dipole 34',msqd(34,0,0)
         write(6,*) 'CF    dipole 35',msqd(35,0,0)
         write(6,*) 'CF    dipole 36',msqd(36,0,0)
         write(6,*) 'CF    dipole 37',msqd(37,0,0)
         write(6,*) 'CF    dipole 38',msqd(38,0,0)
         write(6,*) 'CF    dipole 39',msqd(39,0,0)
         write(6,*) 'CF    dipole 40',msqd(40,0,0)
         write(6,*) 'CF    dipole 41',msqd(41,0,0)
         
c     write(6,*) 'CF      ratio',(msqd(1,0,0)+msqd(2,0,0)
c     .        +msqd(3,0,0)+msqd(4,0,0)+msqd(5,0,0)+msqd(6,0,0)
c     .        +msqd(7,0,0))/msq(0,0)
c         write(6,*) 'CF   s4 ratio',(msqd(1,0,0)
c     .        +msqd(3,0,0)+msqd(5,0,0)+msqd(6,0,0))/msq(0,0)
c         write(6,*) 'CF   s5 ratio',(msqd(2,0,0)
c     .        +msqd(4,0,0)+msqd(5,0,0)+msqd(7,0,0))/msq(0,0)
c         write(6,*) 'CF      ratio',(msqd(10,0,0)
c     .        )/msq(0,0)
c         write(6,*) 'CF      ratio',(msqd(15,0,0)
c     .        )/msq(0,0)

c         write(6,*) 'ratio A2 dipoles ',(msqd(8,0,0)+msqd(9,0,0)
c     .        +msqd(10,0,0)+msqd(11,0,0)+msqd(12,0,0)+msqd(13,0,0)
c     .        +msqd(14,0,0)+msqd(15,0,0)+msqd(16,0,0))/msq(0,0)

c         write(6,*)'ratio 1 = ',msqd(17,0,0)/msqd(8,0,0)
c         write(6,*)'ratio 2 = ',msqd(18,0,0)/msqd(9,0,0)
c         write(6,*)'ratio 3 = ',msqd(17,0,0)/msqd(5,0,0)
c         write(6,*)'ratio 4 = ',msqd(18,0,0)/msqd(5,0,0)

c         write(6,*)'ratio c_15c_156/c_156  = ',msqd(19,0,0)/msqd(8,0,0)
c         write(6,*)'ratio c_16c_165/c_165  = ',msqd(20,0,0)/msqd(8,0,0)
c         write(6,*)'ratio c_25c_256/c_256  = ',msqd(21,0,0)/msqd(9,0,0)
c         write(6,*)'ratio c_26c_265/c_265  = ',msqd(22,0,0)/msqd(9,0,0)

c$$$         write(6,*)'ratio c_15c_1526/c_1526  = '
c$$$     .        ,msqd(28,0,0)/msqd(10,0,0)
c$$$         write(6,*)'ratio c_16c_1625/c_1625  = '
c$$$     .        ,msqd(29,0,0)/msqd(11,0,0)
c$$$         write(6,*)'ratio c_25c_2516/c_2516  = '
c$$$     .        ,msqd(30,0,0)/msqd(11,0,0)
c$$$         write(6,*)'ratio c_26c_2615/c_2615  = '
c$$$     .        ,msqd(31,0,0)/msqd(10,0,0)
c$$$
c$$$         write(6,*) 'sum A2 dipoles   ',(msqd(8,0,0)+msqd(9,0,0)
c$$$     .        +msqd(10,0,0)+msqd(11,0,0)+msqd(12,0,0)+msqd(13,0,0)
c$$$     .        +msqd(14,0,0)+msqd(15,0,0)+msqd(16,0,0))
c$$$         
c$$$         
c$$$         write(6,*)'sum of c15A2 in A12 = ',(msqd(19,0,0)+msqd(24,0,0)
c$$$     .        +msqd(28,0,0))
c$$$         write(6,*)'sum of c16A2 in A12 = ',(msqd(20,0,0)+msqd(25,0,0)
c$$$     .        +msqd(29,0,0))
c$$$         write(6,*)'sum of c25A2 in A12 = ',(msqd(21,0,0)+msqd(26,0,0)
c$$$     .        +msqd(30,0,0))
c$$$         write(6,*)'sum of c26A2 in A12 = ',(msqd(22,0,0)+msqd(27,0,0)
c$$$  .        +msqd(31,0,0))

c         write(6,*)'ratio s6c156/c156 = ',msqd(32,0,0)/msqd(8,0,0)
c         write(6,*)'ratio s5c165/c156 = ',msqd(33,0,0)/msqd(8,0,0)
c         write(6,*)'ratio s6c256/c256 = ',msqd(34,0,0)/msqd(9,0,0)
c         write(6,*)'ratio s5c265/c256 = ',msqd(35,0,0)/msqd(9,0,0)

c         write(6,*)'rat s6cs156/cs156 = ',msqd(36,0,0)
c     .        /(msqd(12,0,0))
c         write(6,*)'rat s5cs165/cs165 = ',msqd(37,0,0)
c     .        /(msqd(13,0,0))
c         write(6,*)'rat s6s56/s56 = ',msqd(40,0,0)
c     .        /(msqd(16,0,0))
c         write(6,*)'rat s5s65/s56 = ',msqd(41,0,0)
c     .        /(msqd(16,0,0))
c         write(6,*)'rat c26s6cs156/s6cs156 = '
c     .        ,msqd(32,0,0)/msqd(36,0,0)
c         write(6,*)'rat c25s5cs165/s5cs165 = '
c     .        ,msqd(33,0,0)/msqd(37,0,0)
c         write(6,*)'rat c16s6cs256/s6cs256 = '
c     .        ,msqd(34,0,0)/msqd(38,0,0)
c         write(6,*)'rat c15s5cs265/s5cs265 = '
c     .        ,msqd(35,0,0)/msqd(39,0,0)

         A1sub = 0d0
         do jj=1,7
            A1sub = A1sub + msqc(jj,0,0)
         enddo
         A2sub = 0d0
         do jj=8,16
            A2sub = A2sub + msqc(jj,0,0)
         enddo
         A12sub = 0d0
         do jj=17,41
            A12sub = A12sub + msqc(jj,0,0)
         enddo

         write(6,*)'p1 = ',p(1,:)
         write(6,*)'p2 = ',p(2,:)
         write(6,*)'p3 = ',p(3,:)+p(4,:)
         write(6,*)'p4 = ',p(5,:)
         write(6,*)'p5 = ',p(6,:)
         
         write(6,*) 'matrix element ',msq(0,0)
         write(6,*) 'sum A1         ',A1sub
         write(6,*) 'sum A2         ',A2sub
         write(6,*) 'sum A12        ',A12sub
         write(6,*) 'sum A1+A2+A12  ',A1sub+A2sub+A12sub
         write(6,*) '   ratio       ',(A1sub+A2sub+A12sub)/msq(0,0)
         write(6,*) '    diff       ',msq(0,0)-(A1sub+A2sub+A12sub)
         
         pause
         endif
         
      endif

      
      do nd=0,ndmax
      xmsq(nd)=0d0
      do j=-1,1
      do k=-1,1
      xmsq_bypart(nd,j,k)=0d0
      enddo
      enddo

      enddo
      
      currentPDF=0
            
      flux=fbGeV2/(two*xx1*xx2*W)

c--- initialize a PDF set here, if calculating errors
  777 continue    
      do nd=0,ndmax
      xmsq(nd)=0d0
      enddo
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
         
c--- calculate PDF's  
      if (dynamicscale) then
        do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
          if (dipscale(nd) .lt. 1d-8) then        
c--- in case dipole is not used, set up dummy value of scale for safety
c--- and set all PDF entries to zero
          dipscale(nd)=dipscale(0)
          do j=-nf,nf
            fx1(j)=0d0
            fx2(j)=0d0
          enddo
        else
            call fdist(ih1,xx1,dipscale(nd),fx1)
            call fdist(ih2,xx2,dipscale(nd),fx2)
            do j=-nf,nf
            dipfx1(nd,j)=fx1(j)
            dipfx2(nd,j)=fx2(j)
          enddo
        endif
      enddo

      else
            call fdist(ih1,xx1,facscale,fx1)
            call fdist(ih2,xx2,facscale,fx2)

      endif      
            
      do j=-nflav,nflav
         do k=-nflav,nflav

            if ((j.ne.0).or.(k.ne.0)) cycle
            
      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif      
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (omitgg) then 
      if ((j.eq.0) .and. (k.eq.0)) goto 20
      endif


      if (realonly) then 
        xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        do nd=1,ndmax
        xmsq(nd)=0d0
        enddo
      elseif (virtonly) then
         xmsq(0)=0d0
         do nd=1,ndmax
         if (dynamicscale) then         
             xmsq(nd)=xmsq(nd)+dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
         else
             xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
         endif
         enddo
      else

         if     (j .gt. 0) then
           sgnj=+1
         elseif (j .lt. 0) then
           sgnj=-1
         else
           sgnj=0
         endif
         if     (k .gt. 0) then
           sgnk=+1
         elseif (k .lt. 0) then
           sgnk=-1
         else
           sgnk=0
         endif


         xmsqjk=fx1(j)*fx2(k)*msq(j,k)
       
         xmsq(0)=xmsq(0)+xmsqjk

         
         if (currentPDF .eq. 0) then
           xmsq_bypart(0,sgnj,sgnk)=xmsq_bypart(0,sgnj,sgnk)+xmsqjk
         endif
         do nd=1,ndmax
         if (dynamicscale) then         
             xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
          else
             xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
         endif
           xmsq(nd)=xmsq(nd)+xmsqjk
           if (currentPDF .eq. 0) then
             xmsq_bypart(nd,sgnj,sgnk)=xmsq_bypart(nd,sgnj,sgnk)+xmsqjk
           endif
        enddo
         
      endif
 20   continue

      enddo
      enddo

c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        do nd=0,ndmax 
          PDFxsec_nd(currentPDF,nd)=xmsq(nd)
        enddo
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
c--- reset xmsq to the central PDF values
        do nd=0,ndmax 
          xmsq(nd)=PDFxsec_nd(0,nd)
        enddo
      endif    

      realint=0d0
      xint=0d0
      valsum=0d0 ! running total of weights at this point

c--- zero out temporary histograms
      if (bin) call zerorealhistos

C---Trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
         xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat         
         failed=.false.
        
        if (nd .eq. 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) .eq. 0d0) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) .eq. 0d0) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo

c---  the following logical should be activated inside dipoles
c---  if alpha parameters are activated the computation could
c---  be skipped
c---  incldip(nd)=.true.
          if (incldip(nd)) incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
        endif

 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997         
        endif

c---if it does, add to total
        xint=xint+xmsq(nd)!*reweight

        do j=-1,1
        do k=-1,1
          lord_bypart(j,k)=lord_bypart(j,k)+
     .         wgt*flux*pswt*xmsq_bypart(nd,j,k)/BrnRat
        enddo
        enddo

        val=xmsq(nd)*wgt
        val2=val**2
      
        valsum=valsum+val
      
        
c--- update PDF errors
        if (PDFerrors) then
          do currentPDF=0,maxPDFsets        
          PDFwgt(currentPDF)=
     .       flux*pswt*PDFxsec_nd(currentPDF,nd)/BrnRat*wgt/itmx
          PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .       +PDFwgt(currentPDF)
          enddo           
       endif



c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)          
          call nplotter(pjet,val,val2,nd)
        endif
c---otherwise, skip contribution
 997    continue
      enddo


c--- add temporary histograms to cumulative totals
      if (bin) call addrealhistos(wgt)

c--- update the maximum weight so far, if necessary
      if (dabs(valsum) .gt. wtmax) then
         wtmax=dabs(valsum)

c$$$!---  DEBUG
c$$$         if (wtmax.gt.60d0) then
c$$$            call dotem(6,p,s)
c$$$            write(1278,*) s
c$$$            call mywriteout(p,1278)
c$$$         write(1278,*) 'weight is ',wtmax
c$$$         write(1278,*) 
c$$$         write(1278,*) 
c$$$      
c$$$         write(1278,*) 'wtmax          ',wtmax,wgt
c$$$         write(1278,*) 'val            ',val
c$$$         do jj=0,100
c$$$            write(1278,*) 'xmsq           ',jj,xmsq(jj)
c$$$         enddo
c$$$         write(1278,*) 'matrix element ',msq(0,0)
c$$$         A1sub = 0d0
c$$$         do jj=1,7
c$$$            A1sub = A1sub + msqc(jj,0,0)
c$$$         enddo
c$$$         A2sub = 0d0
c$$$         do jj=8,16
c$$$            A2sub = A2sub + msqc(jj,0,0)
c$$$         enddo
c$$$         A12sub = 0d0
c$$$         do jj=17,41
c$$$            A12sub = A12sub + msqc(jj,0,0)
c$$$         enddo
c$$$         write(1278,*) 'sum A1         ',A1sub
c$$$         write(1278,*) 'sum A2         ',A2sub
c$$$         write(1278,*) 'sum A12        ',A12sub
c$$$         write(1278,*) 'sum A1+A2+A12  ',A1sub+A2sub+A12sub
c$$$         write(1278,*) '   ratio       ',(A1sub+A2sub+A12sub)/msq(0,0)
c$$$         write(1278,*) '    diff       ',msq(0,0)-(A1sub+A2sub+A12sub)
c$$$         call transformC2_colorful(p,phrans,1,6)
c$$$         write(1278,*) 'y16   = ',-dot(p,1,6)/dot(p,1,2)
c$$$         write(1278,*) 'y2h5h = ',-dot(phrans,2,5)/dot(phrans,1,2)
c$$$         write(1278,*) 'x6    = ',-(dot(p,1,6)+dot(p,2,6))/dot(p,1,2)
c$$$         write(1278,*) 'x5h   = ',-(sprod(phrans(5,:),p(1,:)+p(2,:)))
c$$$     .        /(sprod(phrans(2,:),p(1,:)+p(2,:)))
c$$$         call mywriteout(phrans,1278)
c$$$!--- END DEBUG
c$$$         stop
c$$$      endif
         
      endif
      
      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
            
      return

 999  realint=0d0
      ntotzero=ntotzero+1
 
      return
      
      end
















