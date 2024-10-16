      double precision function virtint(r,wgt)
      implicit none
      include 'constants.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'msq_cs.f'
      include 'msq_struc.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'efficiency.f'
      include 'process.f'
      include 'maxwt.f'
      include 'limits.f'
      include 'nflav.f'
      include 'b0.f'
      include 'PDFerrors.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'nores.f'
      include 'ewcouple.f'
      include 'phasemin.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'runstring.f' 
      include 'ptilde.f'
      include 'dipolescale.f'
      include 'order.f'
      include 'part.f'
      include 'reweight.f'
      include 'singular.f'

      integer nproc
      common/nproc/nproc
      double precision mqq(0:2,fn:nf,fn:nf),mx
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision AP(-1:1,-1:1,3),APz(-1:1,-1:1,3),APt(-1:1,-1:1,3)
     .     ,APqg_mass
      integer ih1,ih2,i,j,k,m,n,cs,ics,csmax,
     .     nvec,is,iq,ia,ib,ic,ii,jj,nd
      double precision A1sub
      double precision p(mxpart,4),qv(mxpart,4),pjet(mxpart,4),r(mxdim),
     . val,val2,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf),
     . fx1t(-nf:nf),fx2t(-nf:nf),xmsq(0:maxd),xint,valsum,W,sqrts,
     . dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf)
      double precision pswt,xjac,rscalestart,fscalestart,m3,m4,m5,
     . wgt,msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . msqs(maxd,-nf:nf,-nf:nf),msqx(maxd,-nf:nf,-nf:nf,0:2),
     . msqvs(maxd,-nf:nf,-nf:nf),msqd(maxd,-nf:nf,-nf:nf),
     . msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr
      double precision xx(2),z,t,x1onz,x2onz,x1ont,x2ont,flux,omz,omt,
     . BrnRat,xmsq_old,zmsq_old,tmp,ptmp,pttwo,xmsqjk,ymsq,zmsq,xx1,xx2
      double precision xmsq_bypart(0:maxd,-1:1,-1:1),
     . lord_bypart(-1:1,-1:1)
      integer nshot,sgnj,sgnk
      logical bin,first,includedipole,checkpiDpjk,failed
      logical incldip(0:maxd),includevirt
      common/incldip/incldip
      double precision QandGint
      character*4 mypart
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/bypart/lord_bypart
      common/mypart/mypart
      data p/56*0d0/
      data nshot/1/
      data first/.true./
      save first,rscalestart,fscalestart
      logical purevirt_cf_new
      common/useropt_cf_new/purevirt_cf_new
      data purevirt_cf_new/.false./

      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      pswt=0d0
      virtint=0d0

      W=sqrts**2

      call gen_virt_ps(r,p,pswt,*999)
            
      nvec=npart+2
      call dotem(nvec,p,s)

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)


      if (check.eqv..false.) then
      
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
         
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      includevirt = includedipole(0,p)
      incldip(0)  = includevirt

      else
         includevirt = .true.
      endif

      if (includevirt .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
          msqv(j,k)=0d0
        enddo
        enddo
        if ((order.eq.1).and.(nproc.eq.709)) goto 999
      endif
      
      

      if (dynamicscale) then
        call scaleset(rscalestart,fscalestart,p)
        dipscale(0)=facscale
      endif


      z=0.5d0
      t=0.5d0
      xjac=1d0
      if ((nproc.eq.709.and.order.eq.1).or.order.gt.1) then
         z=r(ndim-1)**2
         t=r(ndim)**2
         if (nshot .eq. 1) then
            z=0.95d0
            t=0.95d0
         endif
         xjac=four*dsqrt(z*t)
      endif
      
      omz=1d0-z
      omt=1d0-t

      flux=fbGeV2/(two*xx(1)*xx(2)*W)
            
   12 continue
c--- point to restart from when checking epsilon poles

c--- correction to epinv from AP subtraction when mu_FAC != mu_REN,
c--- corresponding to subtracting -1/epinv*Pab*log(musq_REN/musq_FAC)
      epcorr=epinv+2d0*dlog(scale/facscale)
      
c--- debug: note that we still have to cancel the sing. of integrated A1
      if (order.le.1) then
         APz(g,g,:)=0d0
         APt(g,g,:)=0d0
      endif

      if ((nproc.eq.709.and.order.eq.1).or.order.gt.1) then
         APz(g,g,1)=+ason2pi*b0*epcorr
         APz(g,g,2)=+ason2pi*xn*2d0*(1d0/z+z*omz-2d0)*epcorr
         APz(g,g,3)=+ason2pi*xn*2d0/omz*epcorr
         
         APt(g,g,1)=+ason2pi*b0*epcorr
         APt(g,g,2)=+ason2pi*xn*2d0*(1d0/t+t*omt-2d0)*epcorr
         APt(g,g,3)=+ason2pi*xn*2d0/omt*epcorr

c---  DEBUG: return ep^(-2) part
c         APz(g,g,:)=0d0
c         APt(g,g,:)=0d0   
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         epcorr=1d0
c               
c         APz(g,g,1)=+ason2pi*b0*epcorr
c         APz(g,g,2)=+ason2pi*xn*2d0*(1d0/z+z*omz-2d0)*epcorr
c         APz(g,g,3)=+ason2pi*xn*2d0/omz*epcorr
c         
c         APt(g,g,1)=+ason2pi*b0*epcorr
c         APt(g,g,2)=+ason2pi*xn*2d0*(1d0/t+t*omt-2d0)*epcorr
c         APt(g,g,3)=+ason2pi*xn*2d0/omt*epcorr
c---  END DEBUG
         
      endif
      
      do ia=-1,+1
      do ib=-1,+1
      do ic=-1,+1
      do is=1,4
        T1(ia,ib,ic,is)=0d0
        T2(ia,ib,ic,is)=0d0
        Tp1(:,ia,ib,ic,is,:)=0d0
        Tp2(:,ia,ib,ic,is,:)=0d0
      enddo
      enddo
      enddo
      enddo      
      
      
c 44   continue
      
      
c--- Calculate the required matrix elements      
      if (case .eq. 'H_gaga') then
         
c--- debug: to check only the correction
         if (includevirt) call gg_hgagag(p,msq)
         if (nproc.eq.710) then
            call gg_hgamgam_gs_cf(p,msqs,msqvs,msqx)
            if (order.eq.2) then
               if (includevirt) call gg_hgagag_v(p,msqv)
               call gg_hg_z_cf(p,z,t)
               do nd=1,ndmax
                  call gg_hg_z_wrap(nd,p,z,t)
               enddo
            endif
         else if (nproc.eq.709) then
            if (includevirt.and.(order.eq.1)) call gg_hgagag_v(p,msqv)
            if (order.eq.1) call gg_hg_z_cf(p,z,t)
         endif

         if (check) then
            xx(1) = 0.5d0
            xx(2) = 0.25d0
         if (order.eq.1) msqd=msqs
         if (order.eq.2) msqd=msqvs

            
         write(6,*) 'CF    dipole 1',msqd(1,0,0)
         write(6,*) 'CF    dipole 2',msqd(2,0,0)
         write(6,*) 'CF    dipole 3',msqd(3,0,0)

         A1sub = 0d0
         do jj=1,3
            A1sub = A1sub + msqd(jj,0,0)
         enddo

         write(6,*)'p1 = ',p(1,:)
         write(6,*)'p2 = ',p(2,:)
         write(6,*)'p3 = ',p(3,:)+p(4,:)
         write(6,*)'p4 = ',p(5,:)

         if (order.eq.1) mx=msq(0,0)
         if (order.eq.2) mx=msqv(0,0)
         write(6,*) 'matrix element ',mx
         write(6,*) 'sum A1         ',A1sub
         write(6,*) '   ratio       ',A1sub/mx
         write(6,*) '    diff       ',mx-A1sub
         

         
c         pause
         endif

      endif

      
C---initialize to zero
      do j=-1,1
         do k=-1,1
            xmsq_bypart(:,j,k)=0d0
         enddo
      enddo

      currentPDF=0
            
c--- initialize a PDF set here, if calculating errors
  777 continue
      do nd=0,ndmax
         xmsq(nd)=0d0
      enddo

      
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif


c---  calculate PDF's  
      if (dynamicscale) then
         do nd=ndmax,0,-1       ! so that fx1,fx2 correct for real kinematics
            if (dipscale(nd) .lt. 1d-8) then        
c---  in case dipole is not used, set up dummy value of scale for safety
c---  and set all PDF entries to zero
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
         
c---  calculate PDF's         
         call fdist(ih1,xx(1),facscale,fx1)
         call fdist(ih2,xx(2),facscale,fx2)
         
         
      endif
      
      
      do j=-nf,nf
         fx1z(j)=0d0
         fx2z(j)=0d0
         fx1t(j)=0d0
         fx2t(j)=0d0
      enddo

      if ((nproc.eq.709.and.order.eq.1).or.order.gt.1) then
         
         if (z .gt. xx(1)) then
            x1onz=xx(1)/z
            call fdist(ih1,x1onz,facscale,fx1z)
         endif
         if (z .gt. xx(2)) then
            x2onz=xx(2)/z
            call fdist(ih2,x2onz,facscale,fx2z)
         endif         
         if (t .gt. xx(1)) then
            x1ont=xx(1)/t
            call fdist(ih1,x1ont,facscale,fx1t)
         endif         
         if (t .gt. xx(2)) then
            x2ont=xx(2)/t
            call fdist(ih2,x2ont,facscale,fx2t)
         endif         
         
      endif

      ymsq=0d0
      zmsq=0d0
      
      do j=-nflav,nflav
         do k=-nflav,nflav
            
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
            
            
c---  SUM BY TOTAL MATRIX ELEMENTS: everything else
            if     ((j .gt. 0) .and. (k.gt.0)) then
               
            elseif ((j .eq. g) .and. (k.eq.g)) then
               msq_qg=msq(+5,g)+msq(+4,g)+msq(+3,g)+msq(+2,g)+msq(+1,g)
     1              +msq(-5,g)+msq(-4,g)+msq(-3,g)+msq(-2,g)+msq(-1,g)
               msq_gq=msq(g,+5)+msq(g,+4)+msq(g,+3)+msq(g,+2)+msq(g,+1)
     1              +msq(g,-5)+msq(g,-4)+msq(g,-3)+msq(g,-2)+msq(g,-1)
               ymsq=ymsq+(msqv(g,g)
     1              +msq(g,g)*(one+APz(g,g,1)-APz(g,g,3)+T1(g,g,g,1)
     2              +APt(g,g,1)-APt(g,g,3)+T2(g,g,g,1)))*fx1(g)*fx2(g)
     3              +(msq(g,g)*(APz(g,g,2)+APz(g,g,3)+T1(g,g,g,2))
     4              )*fx1z(g)/z*fx2(g)
     5              +(msq(g,g)*(APt(g,g,2)+APt(g,g,3)+T2(g,g,g,2))
     6              )*fx1(g)*fx2t(g)/t
     7              +(msq(g,g)*T1(g,g,g,3)*fx1z(g)/z*fx2z(g)/z)
     8              +(msq(g,g)*T2(g,g,g,3)*fx1t(g)/t*fx2t(g)/t)    
     9              +(msq(g,g)*T2(g,g,g,4)*fx1z(g)/z*fx2t(g)/t)
     1              +(msq(g,g)*T1(g,g,g,4)*fx1z(g)/z*fx2t(g)/t)
               
            endif
            
            if (purevirt_cf_new) then
               ymsq=ymsq-msq(j,k)*fx1(j)*fx2(k)
            endif  
            
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
            
            
            xmsq(0) = xmsq(0) + ymsq
            
            
            if (currentPDF .eq. 0) then
               xmsq_bypart(0,sgnj,sgnk)=xmsq_bypart(0,sgnj,sgnk)+ymsq
            endif
            
            if ((nproc.eq.710.and.order.eq.1).or.order.gt.1) then
               
               do nd=1,ndmax
                  zmsq=0d0
                  if (dynamicscale) then         
                     xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqvs(nd,j,k))
                  else
                     
                     if     ((j .gt. 0) .and. (k.gt.0)) then
                        
                     elseif ((j .eq. g) .and. (k.eq.g)) then
C--   gg
                        if (nd.lt.3) then

                           zmsq=zmsq+(-msqvs(nd,g,g)
     1                          -msqs(nd,g,g)*(one
     2                          +APz(g,g,1)-APz(g,g,3)
     3                          +Tp1(nd,g,g,g,1,0)
     4                          +APt(g,g,1)-APt(g,g,3)
     5                          +Tp2(nd,g,g,g,1,0)
     6                          ))
     7                          *fx1(g)*fx2(g)
     8                          -(msqs(nd,g,g)*(APz(g,g,2)+APz(g,g,3)
     9                          +Tp1(nd,g,g,g,2,0))
     1                          )*fx1z(g)/z*fx2(g)
     2                          -(msqs(nd,g,g)*(APt(g,g,2)+APt(g,g,3)
     3                          +Tp2(nd,g,g,g,2,0))
     4                          )*fx1(g)*fx2t(g)/t
     5                          -msqs(nd,g,g)*Tp1(nd,g,g,g,3,0)
     6                          *fx1z(g)/z*fx2z(g)/z
     7                          -msqs(nd,g,g)*Tp2(nd,g,g,g,3,0)
     8                          *fx1t(g)/t*fx2t(g)/t    
     9                          -msqs(nd,g,g)*Tp2(nd,g,g,g,4,0)
     1                          *fx1z(g)/z*fx2t(g)/t
     2                          -msqs(nd,g,g)*Tp1(nd,g,g,g,4,0)
     3                          *fx1z(g)/z*fx2t(g)/t
                           
                        elseif (nd.eq.3) then 
                           zmsq=zmsq+(-msqvs(nd,g,g)
     1                          -msqs(nd,g,g)*(one+APz(g,g,1)-APz(g,g,3)
     2                          +APt(g,g,1)-APt(g,g,3)
     3                          ))
     4                          *fx1(g)*fx2(g)
     5                          -(msqs(nd,g,g)*(APz(g,g,2)+APz(g,g,3))
     6                          )*fx1z(g)/z*fx2(g)
     7                          -(msqs(nd,g,g)*(APt(g,g,2)+APt(g,g,3))
     8                          )*fx1(g)*fx2t(g)/t
                           
                           
                           do i=0,2
                              zmsq=zmsq+(
     1                             -msqx(nd,g,g,i)*(
     2                             Tp1(nd,g,g,g,1,i)
     3                             +Tp2(nd,g,g,g,1,i)
     4                             ))
     5                             *fx1(g)*fx2(g)
     6                             -(msqx(nd,g,g,i)*(
     7                             Tp1(nd,g,g,g,2,i))
     8                             )*fx1z(g)/z*fx2(g)
     9                             -(msqx(nd,g,g,i)*(
     1                             Tp2(nd,g,g,g,2,i))
     2                             )*fx1(g)*fx2t(g)/t
     3                             -msqx(nd,g,g,i)*Tp1(nd,g,g,g,3,i)
     4                             *fx1z(g)/z*fx2z(g)/z
     5                             -msqx(nd,g,g,i)*Tp2(nd,g,g,g,3,i)
     6                             *fx1t(g)/t*fx2t(g)/t    
     7                             -msqx(nd,g,g,i)*Tp2(nd,g,g,g,4,i)
     8                             *fx1z(g)/z*fx2t(g)/t
     9                             -msqx(nd,g,g,i)*Tp1(nd,g,g,g,4,i)
     1                             *fx1z(g)/z*fx2t(g)/t
                           enddo
                           
                        endif
                        
                     endif
            
                  endif

                  xmsq(nd)=xmsq(nd)+zmsq
                  if (currentPDF .eq. 0) then
                     xmsq_bypart(nd,sgnj,sgnk)=xmsq_bypart(nd,sgnj,sgnk)+zmsq
                  endif

               enddo    ! do nd=1,ndmax
               
            endif    ! if (order.gt.1) then
            
 20         continue
            
         enddo
      enddo
      
      
      if (check) then
         write(6,*)'ymsq = ',xmsq(0)
         write(6,*)'zmsq = ',xmsq(1)+xmsq(2)+xmsq(3)
         write(6,*)'sum  = ',xmsq(0)+xmsq(1)+xmsq(2)+xmsq(3)
         write(6,*)'rat  = ',(xmsq(1)+xmsq(2)+xmsq(3))/xmsq(0)
         pause
      endif
      
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

      virtint=0d0
      xint=0d0
      
      
c      if (currentPDF .eq. 0) then
c        virtint=flux*xjac*pswt*xmsq/BrnRat
c      endif
            
      valsum=0d0 ! running total of weights at this point


c--- code to check that epsilon poles cancel      
      if (nshot .eq. 1) then
        if (xmsq(0) .eq. 0d0) goto 999
        xmsq_old=xmsq(0)
        zmsq_old=xmsq(1)+xmsq(2)+xmsq(3)
        nshot=nshot+1
        epinv=0d0
        epinv2=0d0
        goto 12
      elseif (nshot .eq. 2) then
        nshot=nshot+1        
        if (abs(xmsq_old/xmsq(0)-1d0) .gt. 1d-6) then
          write(6,*) 'epsilon fails to cancel'
          write(6,*) 'xmsq (epinv=large) = ',xmsq_old
          write(6,*) 'xmsq (epinv=zero ) = ',xmsq(0)
C          stop
        else
          write(6,*) 'Poles cancelled!'
c          write(6,*) 'xmsq (epinv=large) = ',xmsq_old
c          write(6,*) 'xmsq (epinv=zero ) = ',xmsq
c          pause
       endif

       if (abs(zmsq_old/(xmsq(1)+xmsq(2)+xmsq(3))-1d0) .gt. 1d-6) then
          write(6,*) 'epsilon fails to cancel in z'
          write(6,*) 'zmsq (epinv=large) = ',zmsq_old
          write(6,*) 'zmsq (epinv=zero ) = ',zmsq
C          stop
        else
          write(6,*) 'Poles cancelled in z!'
c          write(6,*) 'zmsq (epinv=large) = ',zmsq_old
c          write(6,*) 'zmsq (epinv=zero ) = ',zmsq
c          pause
        endif
      endif



      
c--- zero out temporary histograms
      if (bin) call zerorealhistos
      
c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
         xmsq(nd)=xmsq(nd)*flux*xjac*pswt/BrnRat               
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
                 qv(j,k)=ptilde(nd,j,k)
              enddo
           enddo

c---  the following logical should be activated inside dipoles
c---  if alpha parameters are activated the computation could
c---  be skipped
c      incldip(nd)=.true.
           if (incldip(nd)) incldip(nd)=includedipole(nd,qv)
           if (incldip(nd) .eqv. .false.) failed=.true.
c          write(6,*) nd,'incldip(',nd,')=',incldip(nd)
          
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


        
        xint=xint+xmsq(nd)*reweight        

      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=lord_bypart(j,k)+
     .       wgt*flux*xjac*pswt*xmsq_bypart(nd,j,k)/BrnRat
      enddo
      enddo

      val=xmsq(nd)*wgt 
      val2=val**2 

      valsum = valsum + val


c--- update PDF errors
        if (PDFerrors) then
          do currentPDF=0,maxPDFsets        
          PDFwgt(currentPDF)=
     .       flux*xjac*pswt*PDFxsec_nd(currentPDF,nd)/BrnRat*wgt/itmx
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
      endif


      virtint=xint
      virtint = virtint*reweight

c---  update the maximum weight so far, if necessary
      if (abs(val) .gt. wtmax) then
        wtmax=abs(val)
      endif
      
      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


