      double precision function bornint(r,wgt)
      implicit none
      include 'constants.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'agq.f'
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
      include 'Ioperators.f'
      include 'reweight.f'
      include 'singular.f'
      include 'factiny.f'

      double precision mqq(0:2,fn:nf,fn:nf),mx
      double precision msqtot(-nf:nf,-nf:nf,4),tau
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
      double precision QandGint
      double precision tmp1(-nf:nf,-nf:nf,-2:2),tmp2(-nf:nf,-nf:nf,-4:0)
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
      double precision 
     1     res(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0),
     1   resab(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0),
     1   resba(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0)
      character *37 FMT1,FMT2
      character *60 FMT3
      integer myexp,myshot,origord
      save myshot
      logical keepterm1, keepterm2
      double precision tiny,tin2
      tiny=factiny
      tin2=1d-2
            
      if (first) myshot=0
      
      ndmax=0

      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      pswt=0d0
      bornint=0d0

      W=sqrts**2

      if (case.eq.'H_tota') then
         tau=(hmass/sqrts)**2
         xx(1)=(1d0-tau)*r(1)+tau
         xx(2)=tau/xx(1)
         npart=0
         pswt=1d0-tau
         reweight = 1.0d0        
      else
         call gen_born_ps(r,p,pswt,*999)
      endif

      nvec=npart+2
      call dotem(nvec,p,s)

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)

      if (case.ne.'H_tota') then
         
c----reject event if any s(i,j) is too small
         call smalls(s,npart,*999)
      
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
         if (includedipole(0,p) .eqv. .false.) then
            goto 999
         endif

      endif
      

      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)


      z=0.75123d0
      t=0.51231d0
      xjac=1d0
      if (order.gt.0) then
c         if (case.eq.'H_tota') then
c            z=r(2)
c            xjac=1d0
c         else
            z=r(ndim-1)**2
            t=r(ndim)**2
            if (nshot .eq. 1) then
               z=0.75123d0
               t=0.51231d0
            endif
            xjac=four*dsqrt(z*t)
c         endif
      endif

      if (check) myexp=1
      
 101  continue
      
      if (check) then
c         t=1-0.6d0
         t=1-0.12d0**myexp
         z=1-0.1d0**myexp
c         t=1-0.11d0**myexp
         
      endif

      omz=1d0-z
      omt=1d0-t
      
      if (case.eq.'H_tota') then
         flux=fbGeV2*(hmass/sqrts)**2/xx(1)
      else   
         flux=fbGeV2/(two*xx(1)*xx(2)*W)
      endif
      
   12 continue
c--- point to restart from when checking epsilon poles

c--- correction to epinv from AP subtraction when mu_FAC != mu_REN,
c--- corresponding to subtracting -1/epinv*Pab*log(musq_REN/musq_FAC)
      epcorr=epinv+2d0*dlog(scale/facscale)
      
c--- debug: note that we still have to cancel the sing. of integrated A1
      if (order.eq.0) then
         APz(g,g,:)=0d0
         APt(g,g,:)=0d0

      elseif (order.ge.1) then
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
            
      I10op1(:,:,:,:,:)=0d0
      I10op2(:,:,:,:,:)=0d0

      
c--- Calculate the required matrix elements      
      if (case .eq. 'H_gaga') then
        call gg_hgamgam(p,msq)
        if (order.ge.1) then
           call gg_hgamgam_v(p,msqv)
           call iopmat1(p,z,t)
        endif
        if (order.gt.1) then

           keepterm1=.true.
           keepterm2=.true.

           if (keepterm2.eqv..false.) goto 555
           
           if (
     1          z.lt.tiny.or.
     2          t.lt.tiny.or.
     3          z.gt.(1d0-tiny/100d0).or.
     4          t.gt.(1d0-tiny/100d0).or.
     5          abs(9d0*z-t).lt.tiny.or.
     6          abs(9d0*t-z).lt.tiny.or.
     7          abs(z-t).lt.tiny
     8          ) then
              keepterm2=.false.
              goto 555
           endif
           
           call gg_hgaga_v_cf(p,tmp1)
           call gg_hgaga_vv_cf(p,tmp2)
           
c           if (
c     1          abs(res(0,0,0,0,2,2,-1)).gt.1d-7.or.
c     1          z.lt.tin2.or.
c     2          t.lt.tin2.or.
c     3          z.gt.(1d0-tin2).or.
c     4          t.gt.(1d0-tin2).or.
c     5          abs(9d0*z-t).lt.tin2.or.
c     6          abs(9d0*t-z).lt.tin2.or.
c     7          abs(z-t).lt.tin2
c     8          ) then
              call iopmat2(p,z,t)
              call BuildICT(p,tmp1,tmp2,resab,1d0,0d0)
              call BuildICT(p,tmp1,tmp2,resba,0d0,1d0)
c              write(6,*) resab(0,0,0,0,2,2,0)
c     1                  +resba(0,0,0,0,2,2,0)!+tmp2(g,g,-1)
c              write(6,*) 'hp'
c           else
c              call hp_iopmat2(p,z,t)
c              write(6,*) resab(0,0,0,0,1,1,-1)
c     1             +resba(0,0,0,0,1,1,-1) +tmp2(g,g,-1)              
c              call hp_BuildICT(p,tmp1,tmp2,resab,1d0,0d0)
c              call hp_BuildICT(p,tmp1,tmp2,resba,0d0,1d0)
c              write(6,*) 'lp'
c              write(6,*) resab(0,0,0,0,2,2,-1)
c     1                  +resba(0,0,0,0,2,2,-1)
c              write(6,*) resab(0,0,0,0,2,2,0)
c     1                  +resba(0,0,0,0,2,2,0)!+tmp2(g,g,-1)
c              write(6,*) 
c           endif


c           write(6,*)
c           FMT1="('msqv1={',4(F23.17,','),F23.17,'};')"
c           FMT2="('msqv2={',4(F23.17,','),F23.17,'};')"
c           write(6,FMT1) tmp1(0,0,:)
c           write(6,FMT2) tmp2(0,0,:)

 555       continue

              
           if (check) then
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
              write(6,*) '1,1,ab', resab(g,g,g,g,1,1,:) !*fx1(g)*fx2(g)
              write(6,*) '1,1,ba', resba(g,g,g,g,1,1,:) !*fx1(g)*fx2(g)
              write(6,*) '1,1tot', resab(g,g,g,g,1,1,:)+resba(g,g,g,g,1,1,:)+tmp2(0,0,:)
              write(6,*)
              write(6,*) '1,2,ab', resab(g,g,g,g,1,2,:) !*fx1(g)*fx2t(g)/t
              write(6,*) '1,2,ba', resba(g,g,g,g,1,2,:) !*fx1(g)*fx2t(g)/t
              write(6,*) '1,2tot', resab(g,g,g,g,1,2,:)+resba(g,g,g,g,1,2,:)
              write(6,*)              
              write(6,*) '2,1,ab', resab(g,g,g,g,2,1,:) !*fx1z(g)/z*fx2(g)
              write(6,*) '2,1,ba', resba(g,g,g,g,2,1,:) !*fx1z(g)/z*fx2(g)
              write(6,*) '2,1tot', resab(g,g,g,g,2,1,:)+resba(g,g,g,g,2,1,:)
              write(6,*)              
              write(6,*) '2,2,ab', resab(g,g,g,g,2,2,:) !*fx1z(g)/z*fx2t(g)/t
              write(6,*) '2,2,ba', resba(g,g,g,g,2,2,:) !*fx1z(g)/z*fx2t(g)/t
              write(6,*) '2,2tot', resab(g,g,g,g,2,2,:)+resba(g,g,g,g,2,2,:)
              write(6,*)              
              write(6,*) '0,0,ab', resab(g,g,g,g,0,0,:) !*fx1z(g)/z*fx2z(g)/z
              write(6,*) '0,0,ba', resba(g,g,g,g,0,0,:) !*fx1t(g)/t*fx2t(g)/t
              write(6,*)
              write(6,*) 'v     ',tmp1(g,g,:)/msq(g,g)/ason2pi
              write(6,*) 'vv    ',tmp2(g,g,:)/msq(g,g)/ason2pi**2
              write(6,*)
              write(6,*) 'z=',z
              write(6,*) 't=',t
              pause
              myexp = 1+ myexp
              goto 101              
           endif

        endif

      elseif (case .eq. 'H_tota') then
         origord=order
c---  cut to simulate loss of cross section
c         if (order.gt.1) then
c            if (
c     1           z.lt.tiny.or.
c     2           t.lt.tiny.or.
c     3           z.gt.(1d0-tiny/100d0).or.
c     4           t.gt.(1d0-tiny/100d0).or.
c     5           abs(9d0*z-t).lt.4*tiny.or.
c     6           abs(9d0*t-z).lt.4*tiny.or.
c     7           abs(z-t).lt.4*tiny
c     8           ) then
c               order=1
c            endif
c         endif
         call gg_htot(z,msqtot)
c         order=origord
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

      if (order.gt.0) then
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
            


      
c--- SUM BY TOTAL MATRIX ELEMENTS: everything else
            if     ((j .gt. 0) .and. (k.gt.0)) then
               
            elseif ((j .eq. g) .and. (k.eq.g)) then

               if (case.eq.'H_tota') then

                  ymsq=ymsq+msqtot(g,g,1)*fx1(g)*fx2(g)
                  if (z.gt.tau) then
                     ymsq=ymsq-msqtot(g,g,2)*fx1(g)*fx2(g)
                     if (xx(1).gt.tau/z) then
                        ymsq=ymsq+(msqtot(g,g,2)+msqtot(g,g,3))*fx1(g)*fx2z(g)/z**2
                     endif
                  endif

               else

                  if (keepterm1) then
                  ymsq=ymsq+(msqv(g,g)
     1                 +msq(g,g)*(one
     2                 +APz(g,g,1)-APz(g,g,3)+I10op1(g,g,g,1,1)
     3                 +APt(g,g,1)-APt(g,g,3)+I10op2(g,g,g,1,1)))*fx1(g)*fx2(g)
     4                 +(msq(g,g)*(APz(g,g,2)+APz(g,g,3)+I10op1(g,g,g,2,1))
     5                 )*fx1z(g)/z*fx2(g)
     6                 +(msq(g,g)*(APt(g,g,2)+APt(g,g,3)+I10op2(g,g,g,2,1))
     7                 )*fx1(g)*fx2t(g)/t
     8                 +(msq(g,g)*I10op2(g,g,g,2,2)*fx1z(g)/z*fx2t(g)/t)
     9                 +(msq(g,g)*I10op1(g,g,g,2,2)*fx1z(g)/z*fx2t(g)/t)
                  endif
                  
                  if (keepterm2) then
                     ymsq=ymsq+tmp2(g,g,0)*fx1(g)*fx2(g)
     1                    +resab(g,g,g,g,1,1,0)*fx1(g)*fx2(g)
     2                    +resba(g,g,g,g,1,1,0)*fx1(g)*fx2(g)
     3                    +resab(g,g,g,g,1,2,0)*fx1(g)*fx2t(g)/t
     4                    +resba(g,g,g,g,1,2,0)*fx1(g)*fx2t(g)/t
     5                    +resab(g,g,g,g,2,1,0)*fx1z(g)/z*fx2(g)
     6                    +resba(g,g,g,g,2,1,0)*fx1z(g)/z*fx2(g)
     7                    +resab(g,g,g,g,2,2,0)*fx1z(g)/z*fx2t(g)/t
     8                    +resba(g,g,g,g,2,2,0)*fx1z(g)/z*fx2t(g)/t
     9                    +resab(g,g,g,g,0,0,0)*fx1z(g)/z*fx2z(g)/z
     1                    +resba(g,g,g,g,0,0,0)*fx1t(g)/t*fx2t(g)/t
                  endif
               endif
            endif
      
c     if (purevirt_cf_new) then
c     ymsq=ymsq-msq(j,k)*fx1(j)*fx2(k)
c     endif  
      
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
            
      
 20         continue

         enddo
      enddo


      
c---  loop over all PDF error sets, if necessary
      if (PDFerrors) then
         do nd=0,ndmax 
            PDFxsec_nd(currentPDF,nd)=xmsq(nd)
         enddo
         currentPDF=currentPDF+1
         if (currentPDF .le. maxPDFsets) goto 777
c---  reset xmsq to the central PDF values
         do nd=0,ndmax 
            xmsq(nd)=PDFxsec_nd(0,nd)
         enddo
      endif    
      
      bornint=0d0
      xint=0d0
      
      
c      if (currentPDF .eq. 0) then
c        bornint=flux*xjac*pswt*xmsq/BrnRat
c      endif
            
      valsum=0d0 ! running total of weights at this point


c--- code to check that epsilon poles cancel      
      if (order.gt.0) then
         if (nshot .eq. 1) then
            if (xmsq(0) .eq. 0d0) goto 999
            xmsq_old=xmsq(0)
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
c     stop
            else
               write(6,*) 'Poles cancelled!'
c     write(6,*) 'xmsq (epinv=large) = ',xmsq_old
c     write(6,*) 'xmsq (epinv=zero ) = ',xmsq
c     pause
            endif
            
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
     .             wgt*flux*xjac*pswt*xmsq_bypart(nd,j,k)/BrnRat
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

      bornint=xint

c--- reweight if implemented      
      bornint = bornint*reweight

c--- update the maximum weight so far, if necessary
      if (abs(val) .gt. wtmax) then
        wtmax=abs(val)
      endif
      
      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end




      subroutine BuildICT(p,tmp1,tmp2,res,tag1,tag2)
      implicit none
      include 'types.h'
      include 'Ioperators2.f'
      include 'agq.f'
      include 'constants.f'
      include 'qcdcouple.f'
      double precision tmp1(-nf:nf,-nf:nf,-2:2),tmp2(-nf:nf,-nf:nf,-4:0)
      double precision res(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0)
      integer i,j
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),tag1,tag2,as2
      call gg_hgamgam(p,msq)
      as2=ason2pi**2
      do i=0,2
         do j=0,2
            res(g,g,g,g,i,j,-4)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-4)+tag2*I20op2eps(g,g,g,g,j,i,-4))*msq(g,g)*as2
            res(g,g,g,g,i,j,-3)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-3)+tag2*I20op2eps(g,g,g,g,j,i,-3))*msq(g,g)*as2
            res(g,g,g,g,i,j,-2)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-2)+tag2*I20op2eps(g,g,g,g,j,i,-2))*msq(g,g)*as2
            res(g,g,g,g,i,j,-1)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-1)+tag2*I20op2eps(g,g,g,g,j,i,-1))*msq(g,g)*as2
            res(g,g,g,g,i,j, 0)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j, 0)+tag2*I20op2eps(g,g,g,g,j,i, 0))*msq(g,g)*as2
         enddo           
      enddo           
      do i=1,2
         do j=1,2
            res(g,g,g,g,i,j,-4)=res(g,g,g,g,i,j,-4)+
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
            res(g,g,g,g,i,j,-3)=res(g,g,g,g,i,j,-3)+                                  
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
            res(g,g,g,g,i,j,-2)=res(g,g,g,g,i,j,-2)+                                  
     .           +tmp1(g,g, 0)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j, 0)+tag2*I10op2eps(g,g,g,j,i, 0))*ason2pi
            res(g,g,g,g,i,j,-1)=res(g,g,g,g,i,j,-1)+                                  
     .           +tmp1(g,g, 1)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g, 0)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j, 0)+tag2*I10op2eps(g,g,g,j,i, 0))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j, 1)+tag2*I10op2eps(g,g,g,j,i, 1))*ason2pi
            res(g,g,g,g,i,j, 0)=res(g,g,g,g,i,j, 0)+                                  
     .           +tmp1(g,g, 2)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g, 1)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
     .           +tmp1(g,g, 0)*(tag1*I10op1eps(g,g,g,i,j, 0)+tag2*I10op2eps(g,g,g,j,i, 0))*ason2pi
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j, 1)+tag2*I10op2eps(g,g,g,j,i, 1))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j, 2)+tag2*I10op2eps(g,g,g,j,i, 2))*ason2pi
         enddo           
      enddo           
      return
      end
      
      
      subroutine hp_BuildICT(p,tmp1,tmp2,res,tag1,tag2)
      implicit none
      include 'hp_types.h'
      include 'hp_Ioperators2.f'
      include 'agq.f'
      include 'constants.f'
      include 'qcdcouple.f'
      double precision tmp1(-nf:nf,-nf:nf,-2:2),tmp2(-nf:nf,-nf:nf,-4:0)
      double precision res(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0)
      integer i,j
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),tag1,tag2,as2
      call gg_hgamgam(p,msq)
      as2=ason2pi**2
      do i=0,2
         do j=0,2
            res(g,g,g,g,i,j,-4)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-4)+tag2*I20op2eps(g,g,g,g,j,i,-4))*msq(g,g)*as2
            res(g,g,g,g,i,j,-3)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-3)+tag2*I20op2eps(g,g,g,g,j,i,-3))*msq(g,g)*as2
            res(g,g,g,g,i,j,-2)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-2)+tag2*I20op2eps(g,g,g,g,j,i,-2))*msq(g,g)*as2
            res(g,g,g,g,i,j,-1)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j,-1)+tag2*I20op2eps(g,g,g,g,j,i,-1))*msq(g,g)*as2
            res(g,g,g,g,i,j, 0)=
     .           +(tag1*I20op1eps(g,g,g,g,i,j, 0)+tag2*I20op2eps(g,g,g,g,j,i, 0))*msq(g,g)*as2
         enddo           
      enddo           
      do i=1,2
         do j=1,2
            res(g,g,g,g,i,j,-4)=res(g,g,g,g,i,j,-4)+
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
            res(g,g,g,g,i,j,-3)=res(g,g,g,g,i,j,-3)+                                  
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
            res(g,g,g,g,i,j,-2)=res(g,g,g,g,i,j,-2)+                                  
     .           +tmp1(g,g, 0)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j, 0)+tag2*I10op2eps(g,g,g,j,i, 0))*ason2pi
            res(g,g,g,g,i,j,-1)=res(g,g,g,g,i,j,-1)+                                  
     .           +tmp1(g,g, 1)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g, 0)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j, 0)+tag2*I10op2eps(g,g,g,j,i, 0))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j, 1)+tag2*I10op2eps(g,g,g,j,i, 1))*ason2pi
            res(g,g,g,g,i,j, 0)=res(g,g,g,g,i,j, 0)+                                  
     .           +tmp1(g,g, 2)*(tag1*I10op1eps(g,g,g,i,j,-2)+tag2*I10op2eps(g,g,g,j,i,-2))*ason2pi
     .           +tmp1(g,g, 1)*(tag1*I10op1eps(g,g,g,i,j,-1)+tag2*I10op2eps(g,g,g,j,i,-1))*ason2pi
     .           +tmp1(g,g, 0)*(tag1*I10op1eps(g,g,g,i,j, 0)+tag2*I10op2eps(g,g,g,j,i, 0))*ason2pi
     .           +tmp1(g,g,-1)*(tag1*I10op1eps(g,g,g,i,j, 1)+tag2*I10op2eps(g,g,g,j,i, 1))*ason2pi
     .           +tmp1(g,g,-2)*(tag1*I10op1eps(g,g,g,i,j, 2)+tag2*I10op2eps(g,g,g,j,i, 2))*ason2pi
         enddo           
      enddo           
      return
      end
      
      
