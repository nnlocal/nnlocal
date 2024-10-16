      subroutine dipsC3_colorful(nd,p,ip,rp,sp,sub,subv,msq,msqv,
     .     subr_born,subr_corr)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the tripe collinear subtraction terms in              *
*     colorful nnlo with momentum p                                    *
*     Automatically chooses subtraction kind (IFF or FFF)              *
*     Currently only IFF subtraction C_ars implemented                 *    
*     Returns the subs in sub,subv and matrix elements in msq,msqv     *
*     nd labels the collinear configurations                           *
*     ip labels one daughter parton                                    *
*     rp labels the second daughter parton                             *
*     sp labels the third daughter parton                              *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4),subv(6)
      double precision msq(-nf:nf,-nf:nf),msqv(6,-nf:nf,-nf:nf)

      double precision sirs,fac,dot
      double precision vec0(4),vec1(4),vec2(4),vec3(4),vec4(4),
     .     vec5(4),vec6(4),pggg0,pggg1,pggg2,pggg3,pggg4,pggg5,pggg6
      integer nd,ip,rp,sp,nu,i,j,k
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubC3_colorful: ",ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.sp
      if (ip.eq.sp) then
         write(6,*)"Indices ip and sp in dipolesubC3_colorful: ",ip,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubC3_colorful: ",rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubC3_colorful: ",ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check ip, sp not 1 and 2
      if (((ip.eq.1).and.(sp.eq.2)).or.((ip.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices ip and sp in dipolesubC3_colorful: ",ip,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, sp not 1 and 2
      if (((rp.eq.1).and.(sp.eq.2)).or.((rp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices rp and sp in dipolesubC3_colorful: ",rp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif            

c---order, so only ip can be initial state
      if ((ip.gt.rp).or.(ip.gt.sp).or.(rp.gt.sp)) then
         write(6,*)"Wrong order of indices in dipolesubC3_colorful."
         write(6,*)"Correct is ip < rp < sp."
         write(6,*)"Currently ip, rp, sp = ",ip, rp, sp
         write(6,*)"Stopping."
         stop
      endif
      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      subv=0d0
      vec0=0d0
      vec1=0d0
      vec2=0d0
      vec3=0d0
      vec4=0d0
      vec5=0d0
      vec6=0d0
c      call zeromsq(msq,msqv)
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0d0
            do i=1,6
               msqv(i,j,k)=0d0
            enddo
         enddo      
      enddo      
      incldip(nd)=.true.

      sirs=two*(dot(p,ip,rp)+dot(p,ip,sp)+dot(p,rp,sp))

***********************************************************************
************************ INITIAL-FINAL-FINAL **************************
***********************************************************************
      if ((ip .le. 2) .and. (rp .gt. 2) .and. (sp .gt. 2)) then
                  
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC3_colorful(p,ptrans,ip,rp,sp)
         call storeptilde(nd,ptrans)
         
c--   Check to see if this dipole will be included
c     incldip(nd)=includedipole(nd,ptrans)
C--   if not return
c     if (incldip(nd) .eqv. .false.) return
         
c---  if using a dynamic scale, set that scale with dipole kinematics      
         if (dynamicscale) then
            call scaleset(initscale,initfacscale,ptrans)
            dipscale(nd)=facscale
         endif
      
      call pggg_eval(0,p,ip,rp,sp,pggg0,vec0)
      sub(gg) = pggg0
      call subr_born(ptrans,msq)

      call pggg_eval(1,p,ip,rp,sp,pggg1,vec1)
      subv(1) = pggg1
      call subr_corr(ptrans,vec1,ip,msqv(1,:,:))

      call pggg_eval(2,p,ip,rp,sp,pggg2,vec2)
      subv(2) = pggg2
      call subr_corr(ptrans,vec2,ip,msqv(2,:,:))
      
      call pggg_eval(3,p,ip,rp,sp,pggg3,vec3)
      subv(3) = pggg3
      call subr_corr(ptrans,vec3,ip,msqv(3,:,:))

      call pggg_eval(12,p,ip,rp,sp,pggg4,vec4)
      subv(4) = pggg4
      call subr_corr(ptrans,vec4,ip,msqv(4,:,:))

      call pggg_eval(13,p,ip,rp,sp,pggg5,vec5)
      subv(5) = pggg5
      call subr_corr(ptrans,vec5,ip,msqv(5,:,:))

      call pggg_eval(23,p,ip,rp,sp,pggg6,vec6)
      subv(6) = pggg6
      call subr_corr(ptrans,vec6,ip,msqv(6,:,:))

      fac = (2d0*gsq/sirs)**2
      sub = fac*sub
c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
      subv = -fac*subv

      
***********************************************************************
************************* FINAL-FINAL-FINAL ***************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (rp .gt. 2) .and. (sp .gt. 2)) then
         write(6,*)"C_irs^FFF not implemented yet!"
      endif
      
      return
      end


      subroutine dipsC22_colorful(nd,p,ip,rp,jp,sp,sub,subv,subvv,
     .     msq,msqv,msqvv,subr_born,subr_corr,subr_corr_2)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the double unresolved collinear subtraction terms in  *
*     colorful nnlo with momentum p                                    *
*     Automatically chooses subtraction kind (IF,IF, etc.)             *
*     Currently only IF,IF subtracgtion C_ar,bs implemented            *      
*     Returns the subs in sub,subv,subvv                               *
*     and matrix elements in msq,msqv,msqvv                            *
*     nd labels the collinear configurations                           *
*     ip labels one daughter parton of the first coll. pair            *
*     rp labels the other daughter parton of the first coll. pair      *
*     jp labels one daughter parton of the second coll. pair           *
*     sp labels the other daughter parton of the second coll. pair     *      
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
*     subr_corr_2 is the subroutine which call the born process dotted *
*     with vec1 and vec2 for two emitted gluons                        *      
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4),subv(2),
     .     subvv(4)
      double precision msq(-nf:nf,-nf:nf),msqv(2,-nf:nf,-nf:nf),
     .     msqvv(-nf:nf,-nf:nf)

      double precision sir,sjs,fac,dot
      double precision vec0a(4),vec0b(4),vec1a(4),vec1b(4),
     .     vec2a(4),vec2b(4),vec3a(4),vec3b(4),
     .     pggpgg0,pggpgg1,pggpgg2,pggpgg3
      integer nd,ip,rp,jp,sp,nu,i,j,k
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr,subr_corr_2

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubC22_colorful: ",ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.sp
      if (jp.eq.sp) then
         write(6,*)"Indices jp and sp in dipolesubC22_colorful: ",jp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.jp
      if (ip.eq.jp) then
         write(6,*)"Indices ip and jp in dipolesubC22_colorful: ",ip,jp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubC22_colorful: ",rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubC22_colorful: ",ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check jp, sp not 1 and 2
      if (((jp.eq.1).and.(sp.eq.2)).or.((jp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices jp and sp in dipolesubC22_colorful: ",jp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, sp not 1 and 2
      if (((rp.eq.1).and.(sp.eq.2)).or.((rp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices rp and sp in dipolesubC22_colorful: ",rp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif            

c---order, so only ip can be initial state
      if ((ip.gt.rp).or.(jp.gt.sp)) then
         write(6,*)"Wrong order of indices in dipolesubC22_colorful."
         write(6,*)"Correct is ip < rp and jp < sp."
         write(6,*)"Currently ip, rp, jp, sp = ",ip, rp, jp, sp
         write(6,*)"Stopping."
         stop
      endif
      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      subv=0d0
      vec0a=0d0
      vec0b=0d0
      vec1a=0d0
      vec1b=0d0
      vec2a=0d0
      vec2b=0d0
      vec3a=0d0
      vec3b=0d0
      call zeromsq(msq,msqvv)
      do j=-nf,nf
         do k=-nf,nf
            do i=1,2
               msqv(i,j,k)=0d0
            enddo
         enddo      
      enddo      
      incldip(nd)=.true.

      sir=two*dot(p,ip,rp)
      sjs=two*dot(p,jp,sp)

***********************************************************************
******************** INITIAL-FINAL, INITIAL-FINAL *********************
***********************************************************************
      if ((ip.le.2).and.(rp.gt.2).and.(jp.le.2).and.(sp.gt.2)) then
                  
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC22_colorful(p,ptrans,ip,rp,jp,sp)
         call storeptilde(nd,ptrans)
         
c--   Check to see if this dipole will be included
c     incldip(nd)=includedipole(nd,ptrans)
C--   if not return
c     if (incldip(nd) .eqv. .false.) return
         
c---  if using a dynamic scale, set that scale with dipole kinematics      
         if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif
      
      call pgg_pgg_eval(0,p,ip,rp,jp,sp,pggpgg0,vec0a,vec0b)
      sub(gg) = pggpgg0
      call subr_born(ptrans,msq)

      call pgg_pgg_eval(1,p,ip,rp,jp,sp,pggpgg1,vec1a,vec1b)
      subv(1) = pggpgg1
      call subr_corr(ptrans,vec1a,ip,msqv(1,:,:))

      call pgg_pgg_eval(2,p,ip,rp,jp,sp,pggpgg2,vec2a,vec2b)
      subv(2) = pggpgg2
      call subr_corr(ptrans,vec2b,jp,msqv(2,:,:))
      
      call pgg_pgg_eval(3,p,ip,rp,jp,sp,pggpgg3,vec3a,vec3b)
      subvv(gg) = pggpgg3
      call subr_corr_2(ptrans,vec3a,vec3b,ip,jp,msqvv(:,:))

      fac = (2d0*gsq)**2/sir/sjs
      sub = fac*sub
c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
      subv = -fac*subv
      subvv = fac*subvv
      
***********************************************************************
******************** INITIAL-FINAL, FINAL-FINAL ***********************
***********************************************************************
      elseif ((ip.lt.2).and.(rp.gt.2)
     .     .and.(jp.gt.2).and.(sp .gt. 2)) then
         write(6,*)"C_ir,bs^IF,FF not implemented yet!"

***********************************************************************
******************** FINAL-FINAL, INITIAL-FINAL ***********************
***********************************************************************
      elseif ((ip.gt.2).and.(rp.gt.2)
     .     .and.(jp.lt.2).and.(sp .gt. 2)) then
         write(6,*)"C_ir,bs^FF,IF not implemented yet!"
      
***********************************************************************
********************* FINAL-FINAL, FINAL-FINAL ************************
***********************************************************************
      elseif ((ip.gt.2).and.(rp.gt.2)
     .     .and.(jp.gt.2).and.(sp .gt. 2)) then
         write(6,*)"C_ir,bs^FF,FF not implemented yet!"
      endif
      
      return
      end
      

      subroutine dipsCS2_colorful(nd,p,ip,rp,sp,sub,subv,
     .     msq,msqv,subr_born,subr_corr)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the double unresolved collinear-soft subtraction term *
*     including collinear overlaps in colorful nnlo with               *
*     momentum p                                                       *
*     Returns the subs in sub, subv and matrix elements in msq, msqv   *
*     nd labels the collinear configurations                           *
*     ip labels one of the daughter partons of the coll. splitting     *
*     rp labels the other daughter parton                              *
*     sp labels the soft parton                                        *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *      
*     Note: treatment of color correlations is non-generic  for now    *
************************************************************************
      implicit none
      include 'constants.f'
c      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      double precision p(mxpart,4),phrans(mxpart,4),ptrans(mxpart,4),
     .     sub(4),eik,subv
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)

      double precision sir,xa,vec(4),vecsq,sjhkh,sjhs,skhs,dot,
     .     stot,sis,sjs,srs,ssQ,xs
      integer nd,ip,rp,sp,idec1,idec2,nu,j,k,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubCS2_colorful: ",ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.sp
      if (ip.eq.sp) then
         write(6,*)"Indices ip and sp in dipolesubCS2_colorful: ",ip,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubCS2_colorful: ",rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubCS2_colorful: ",ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check ip, sp not 1 and 2
      if (((ip.eq.1).and.(sp.eq.2)).or.((ip.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices ip and sp in dipolesubCS2_colorful: ",ip,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, sp not 1 and 2
      if (((rp.eq.1).and.(sp.eq.2)).or.((rp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices rp and sp in dipolesubCS2_colorful: ",rp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif            

c---order, so only ip can be initial state
      if ((ip.gt.rp).or.(ip.gt.sp)) then
         write(6,*)"Wrong order of indices in dipolesubCS2_colorful."
         write(6,*)"Correct is ip < rp and ip < sp."
         write(6,*)"Currently ip, rp, sp = ",ip, rp, sp
         write(6,*)"Stopping."
         stop
      endif
      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      subv=0d0
      vec=0d0
      eik=0d0
      call zeromsq(msq,msqv)
      incldip(nd)=.true.

c---  note convention: p(1) and p(2) have negative energy
      if (ip.le.2) then
         sir=-two*dot(p,ip,rp)
      else
         sir=two*dot(p,ip,rp)
      endif

***********************************************************************
************************ INITIAL-FINAL-FINAL **************************
***********************************************************************
      if ((ip .le. 2) .and. (rp .gt. 2) .and. (sp .gt. 2)) then

c---  note convention: p(1) and p(2) have negative energy
         xa = 1d0 + (dot(p,rp,1)+dot(p,rp,2))/dot(p,1,2)
         vec(:) = 0d0
         vec(1) = p(rp,1)
         vec(2) = p(rp,2)
         vecsq  = -vec(1)**2 - vec(2)**2

         vec(:) = vec(:)/dsqrt(-vecsq)   
      
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC2_colorful(p,phrans,ip,rp)
         if (rp.lt.sp) then
            call transformS1_colorful(phrans,ptrans,sp-1)
         elseif (sp.lt.rp) then
            call transformS1_colorful(phrans,ptrans,sp)
         endif
         call storeptilde(nd,ptrans)
         
c--   Check to see if this dipole will be included
c     incldip(nd)=includedipole(nd,ptrans)
C--   if not return
c     if (incldip(nd) .eqv. .false.) return
      
c---  if using a dynamic scale, set that scale with dipole kinematics      
         if (dynamicscale) then
            call scaleset(initscale,initfacscale,ptrans)
            dipscale(nd)=facscale
         endif
         
         call subr_born(ptrans,msq)
         call subr_corr(ptrans,vec,ip,msqv)
         
c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
         sub(gg)=2d0*gsq/xa/sir*(2d0*(xa*(1d0-xa) + xa/(1d0-xa)))
         subv   =2d0*gsq/xa/sir*(4d0*(1d0-xa)/xa)      

      endif
      
c---  main loop over eikonal
      nmax = 2
      idec1 = 3
      idec2 = 4
      
      do j = 1,nmax
         if ((j.eq.idec1).or.(j.eq.idec2)) cycle
         sjhs = 2d0*(p(sp,4)*phrans(j,4)
     .        -p(sp,1)*phrans(j,1)
     .        -p(sp,2)*phrans(j,2)
     .        -p(sp,3)*phrans(j,3))
         do k = j+1,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            
            sjhkh = 2d0*dot(phrans,j,k)
            skhs = 2d0*(p(sp,4)*phrans(k,4)
     .           -p(sp,1)*phrans(k,1)
     .           -p(sp,2)*phrans(k,2)
     .           -p(sp,3)*phrans(k,3))
c---  note sum over j>k so extra factor of two
c---  soft
            eik = eik + 2d0*sjhkh/sjhs/skhs
            
         enddo         
      enddo

c---  initial-final-final triple collinear limit of soft-collinear
c      stot=s(1,2)
c      sir=s(ip,rp)
c      sis=s(ip,sp)
c      srs=s(rp,sp)
c      ssQ=s(1,sp)+s(2,sp)

      stot=two*dot(p,1,2)
      sir=two*dot(p,ip,rp)
      sis=two*dot(p,ip,sp)
      srs=two*dot(p,rp,sp)
      ssQ=two*(dot(p,1,sp)+dot(p,2,sp))
      
c---  note convention: p(1) and p(2) have negative energy
      xs = -ssQ/stot - srs/stot*sis/(sir+sis)

c---  sign exchanged so overall color factor is ca**2
      eik = eik + 2d0/sis/xs

c---  initial-final, initial-final double collinear limit of soft-collinear
c      if (ip.eq.1) sjs=s(2,sp)
c      if (ip.eq.2) sjs=s(1,sp)
      if (ip.eq.1) sjs=two*dot(p,2,sp)
      if (ip.eq.2) sjs=two*dot(p,1,sp)
c---  note convention: p(1) and p(2) have negative energy      
      xs = -ssQ/stot
c---  sign exchanged so overall color factor is ca**2
      eik = eik + 2d0/sjs/xs

!---  NEW DEFINITION
      eik=0d0
      sub = -2d0*gsq*eik*sub
      subv = -2d0*gsq*eik*subv
      
      return
      end         
      
      
      subroutine dipsS2_colorful(nd,p,rp,sp,sub,msq,subr_born)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the double unresolved soft subtraction terms          *
*     including soft-collinear overlaps in colorful nnlo with          *
*     momentum p                                                       *
*     Returns the subs in sub and matrix elements in msq               *
*     nd labels the collinear configurations                           *
*     rp labels one of the soft partons                                *
*     sp labels the other soft parton                                  *
*     subr_born is the subroutine which call the born process          *
*     Note: treatment of color correlations is non-generic for now     *
************************************************************************
      implicit none
      include 'constants.f'
c      include 'sprods_com.f'      
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4)
      double precision msq(-nf:nf,-nf:nf)

      double precision sitkt,sitr,sits,sktr,skts,
     .     sjtlt,sjtr,sjts,sltr,slts,srs,symfac,dot,sprod,
     .     xr,xs,xri,xsk,xsi,xrk,xrj,xsj
      double precision tmpSrs(4),tmpCarsSrs(4),tmpCarbsSrs(4),
     .     tmpCSarsSrs(4),tmpCarsCSarsSrs(4)

      
      integer nd,rp,sp,idec1,idec2,nu,i,k,j,l,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born

      tmpSrs=0d0
      tmpCarsSrs=0d0
      tmpCarbsSrs=0d0
      tmpCSarsSrs=0d0
      tmpCarsCSarsSrs=0d0
      
      
c---check rp.gt.2 and sp.gt.2
      if ((rp.le.2).or.(sp.le.2)) then
         write(6,*)"Index rp, sp in dipolesubS2_colorful: ",rp,sp
         write(6,*)"At least one not in the final state, stopping."
         stop
      endif
      
C---  Initialize the dipoles to zero
      sub(:) = 0d0
      msq(:,:) = 0d0
      incldip(nd)=.true.

c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
      call transformS2_colorful(p,ptrans,rp,sp)
      call storeptilde(nd,ptrans)
         
c--   Check to see if this dipole will be included
c     incldip(nd)=includedipole(nd,ptrans)
C--   if not return
c     if (incldip(nd) .eqv. .false.) return
      
c---  if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
         call scaleset(initscale,initfacscale,ptrans)
         dipscale(nd)=facscale
      endif
            
      call subr_born(ptrans,msq)
            
c---  main loop over eikonal
      nmax = 2
      idec1 = 3
      idec2 = 4
      srs = 2d0*dot(p,rp,sp)

      do i = 1,nmax
         if ((i.eq.idec1).or.(i.eq.idec2)) cycle
         sitr = 2d0*(p(rp,4)*ptrans(i,4)
     .              -p(rp,1)*ptrans(i,1)
     .              -p(rp,2)*ptrans(i,2)
     .              -p(rp,3)*ptrans(i,3))
         sits = 2d0*(p(sp,4)*ptrans(i,4)
     .              -p(sp,1)*ptrans(i,1)
     .              -p(sp,2)*ptrans(i,2)
     .              -p(sp,3)*ptrans(i,3))
         do k = i,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            sitkt = 2d0*dot(ptrans,i,k)
            sktr = 2d0*(p(rp,4)*ptrans(k,4)
     .                 -p(rp,1)*ptrans(k,1)
     .                 -p(rp,2)*ptrans(k,2)
     .                 -p(rp,3)*ptrans(k,3))
            skts = 2d0*(p(sp,4)*ptrans(k,4)
     .                 -p(sp,1)*ptrans(k,1)
     .                 -p(sp,2)*ptrans(k,2)
     .                 -p(sp,3)*ptrans(k,3))
            do j = 1,nmax
               if ((j.eq.idec1).or.(j.eq.idec2)) cycle
               sjtr = 2d0*(p(rp,4)*ptrans(j,4)
     .                    -p(rp,1)*ptrans(j,1)
     .                    -p(rp,2)*ptrans(j,2)
     .                    -p(rp,3)*ptrans(j,3))
               sjts = 2d0*(p(sp,4)*ptrans(j,4)
     .                    -p(sp,1)*ptrans(j,1)
     .                    -p(sp,2)*ptrans(j,2)
     .                    -p(sp,3)*ptrans(j,3))
               do l = j,nmax
                  if ((l.eq.idec1).or.(l.eq.idec2)) cycle
                  sjtlt = 2d0*dot(ptrans,j,l)
                  sltr = 2d0*(p(rp,4)*ptrans(l,4)
     .                       -p(rp,1)*ptrans(l,1)
     .                       -p(rp,2)*ptrans(l,2)
     .                       -p(rp,3)*ptrans(l,3))
                  slts = 2d0*(p(sp,4)*ptrans(l,4)
     .                       -p(sp,1)*ptrans(l,1)
     .                       -p(sp,2)*ptrans(l,2)
     .                       -p(sp,3)*ptrans(l,3))                  
c---  abelian soft
c---  color factor is {(Ti.Tk),(Tj.Tl)}=2CA^2, only i.ne.k and j.ne.l conribute
c---  add 2 so full soft gets CA^2
c---  note sum over i <= k and j <= l so extra factor of two for each
                  tmpSrs(gg) = tmpSrs(gg)
     .                 + 4d0*sitkt/sitr/sktr*sjtlt/sjts/slts
               enddo

c---  initial-final soft-collinear limit of double soft
c---  overall color factor ca**2 correct only for f(i) = g
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$               if ((i.le.2).and.(k.le.2).and.(j.le.2)) then
c$$$                  xrj = -sprod(p(rp,:),ptrans(1,:)+ptrans(2,:))
c$$$     .                 /sprod(ptrans(j,:),ptrans(1,:)+ptrans(2,:))
c$$$                  xsj = -sprod(p(sp,:),ptrans(1,:)+ptrans(2,:))
c$$$     .                 /sprod(ptrans(j,:),ptrans(1,:)+ptrans(2,:))
c$$$c---  symmetrized in i and k so i <= k summation gets both contributions
c$$$                  tmpCSarsSrs(gg) = tmpCSarsSrs(gg)
c$$$     .                 + 4d0*sitkt/sits/skts/sjtr/xrj
c$$$     .                 + 4d0*sitkt/sitr/sktr/sjts/xsj
c$$$               endif
               
            enddo
c---  non-abelian soft
c---  color factor is ca*Ti.Tk = ca**2 if i.eq.k, and -ca**2 if i.ne.k
c---  add appropriate (-1) so full soft gets CA^2
c---  note sum over i <= k so extra factor of two for i.ne.k
            symfac = 1d0
            if (i.ne.k) symfac = -2d0
            tmpSrs(gg) = tmpSrs(gg) - symfac*
     .           ((sitkt/sitr/skts/srs
     .           + sitkt/sktr/sits/srs
     .           - sitkt**2/sitr/sktr/sits/skts)
     .           *(1d0 - 1d0/2d0*(sitr*skts+sits*sktr)
     .           /((sitr+sits)*(sktr+skts)))                 
     .           + (sitr*skts+sits*sktr)
     .           /((sitr+sits)*(sktr+skts))/srs**2
     .           - 2d0/srs*sitkt/((sitr+sits)*(sktr+skts)))

c---  initial-final, initial-final double collinear limit of double soft
c---  overall color factor ca**2 correct only for f(i) = g
            if ((i.le.2).and.(k.le.2).and.(i.ne.k)) then
               xri = -sprod(p(rp,:),ptrans(1,:)+ptrans(2,:))
     .              /sprod(ptrans(i,:),ptrans(1,:)+ptrans(2,:))
               xsk = -sprod(p(sp,:),ptrans(1,:)+ptrans(2,:))
     .              /sprod(ptrans(k,:),ptrans(1,:)+ptrans(2,:))
               xsi = -sprod(p(sp,:),ptrans(1,:)+ptrans(2,:))
     .              /sprod(ptrans(i,:),ptrans(1,:)+ptrans(2,:))
               xrk = -sprod(p(rp,:),ptrans(1,:)+ptrans(2,:))
     .              /sprod(ptrans(k,:),ptrans(1,:)+ptrans(2,:))

c---  symmetrized in i and k so i <= k summation gets both contributions
               tmpCarbsSrs(gg) = tmpCarbsSrs(gg)
     .              + 4d0/sitr/skts/xri/xsk
     .              + 4d0/sktr/sits/xrk/xsi
            endif            
         enddo
      
c---  initial-final-final triple collinear limit of double soft
c---  overall color factor ca**2 correct only for f(i) = g
         if (i.le.2) then
            xr = -sprod(p(rp,:),ptrans(1,:)+ptrans(2,:))
     .           /sprod(ptrans(i,:),ptrans(1,:)+ptrans(2,:))
            xs = -sprod(p(sp,:),ptrans(1,:)+ptrans(2,:))
     .           /sprod(ptrans(i,:),ptrans(1,:)+ptrans(2,:))

            tmpCarsSrs(gg) = tmpCarsSrs(gg) - (sits**2*srs*(xr + xs)*
     .           (-(sits*xr*(xr + 2*xs)) + srs*(3*xr + 2*xs)) - 
     .           sitr**3*(-2*sits*xr*xs**3 + 
     .           srs*xs*(2*xr**2 + 3*xr*xs + xs**2)) + 
     .           sitr*sits*(2*sits**2*xr**3*xs + 
     .           5*srs**2*(xr + xs)**2 - 
     .           sits*srs*(3*xr**3 + 3*xr**2*xs + 2*xr*xs**2 + 
     .           2*xs**3)) + 
     .           sitr**2*(-4*sits**2*xr**2*xs**2 + 
     .           srs**2*(2*xr**2 + 5*xr*xs + 3*xs**2) - 
     .           sits*srs*(2*xr**3 + 2*xr**2*xs + 3*xr*xs**2 + 
     .           3*xs**3)))/
     .           (sitr*sits*(sitr + sits)**2*srs**2*xr*xs*
     .           (xr + xs)**2)            
            
c---  triple coll. of soft-coll. of double soft
c---  factor of 2 from symmetrization in r and s
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$            tmpCarsCSarsSrs(gg) = tmpCarsCSarsSrs(gg)
c$$$     .           + 8d0/sitr/sits/xr/xs
         endif 
      enddo
      
!---  NEW CANCELLATIONS IMPLEMENTED: NOTE SIGN OF CarbsSrs
c$$$      sub =
c$$$     .     tmpSrs
c$$$     .     + tmpCarsSrs
c$$$     .     + tmpCSarsSrs
c$$$     .     + tmpCarbsSrs
c$$$     .     + tmpCarsCSarsSrs
      sub =
     .     tmpSrs
     .     + tmpCarsSrs
     .     - tmpCarbsSrs
      
      sub = (2d0*gsq)**2*sub
      
      return
      end         


***********************************************************************
***********************************************************************
      subroutine pggg_eval(tens,p,ip,rp,sp,pggg,vec)

      implicit none
      include 'constants.f'
c      include 'sprods_com.f'
      double precision p(mxpart,4)
      
      integer tens,ip,rp,sp
      double precision stot,s12,s13,s23,srQ,ssQ,xr,xs,xa,xia
      double precision z1,z2,z3,t12,s123,pggg,vec(4),vecsq,dot
      double precision t123, t132, t213, t231, t312, t321

      if ((tens.ne.0).and.(tens.ne.1).and.(tens.ne.2).and.(tens.ne.3)
     .     .and.(tens.ne.12).and.(tens.ne.13).and.(tens.ne.23)) then
         write(6,*)"Non-existent tensor structure ",tens
         write(6,*)"Allowed values are 0, 1, 2, 3, 12, 13, 23"
         write(6,*)"Stopping."
         stop
      endif
      
c      stot=s(1,2)
c      s12=s(ip,rp)
c      s13=s(ip,sp)
c      s23=s(rp,sp)

      stot=two*dot(p,1,2)
      s12=two*dot(p,ip,rp)
      s13=two*dot(p,ip,sp)
      s23=two*dot(p,rp,sp)
     
      s123=s12+s13+s23

      if (ip.le.2) then
c---  INITIAL-FINAL-FINAL
c         srQ=s(1,rp)+s(2,rp)
c         ssQ=s(1,sp)+s(2,sp)

         srQ=two*(dot(p,1,rp)+dot(p,2,rp))
         ssQ=two*(dot(p,1,sp)+dot(p,2,sp))

         
!         xr = -srQ/stot - s23/stot*s12/(s12+s13)
!         xs = -ssQ/stot - s23/stot*s13/(s12+s13)
!---  NEW DEFINITION
         xr = -srQ/stot
         xs = -ssQ/stot
         xa = 1d0-xr-xs
         xia = 1d0+srQ/stot+ssQ/stot+s23/stot

         z1 = 1d0/xa
         z2 = -xr/xa
         z3 = -xs/xa
      endif

      t123 = 2d0*(z1*s23 - z2*s13)/(z1 + z2) + (z1 - z2)/(z1 + z2)*s12 
      t132 = 2d0*(z1*s23 - z3*s12)/(z1 + z3) + (z1 - z3)/(z1 + z3)*s13
      t213 = 2d0*(z2*s13 - z1*s23)/(z2 + z1) + (z2 - z1)/(z2 + z1)*s12
      t231 = 2d0*(z2*s13 - z3*s12)/(z2 + z3) + (z2 - z3)/(z2 + z3)*s23
      t321 = 2d0*(z3*s12 - z2*s13)/(z3 + z2) + (z3 - z2)/(z3 + z2)*s23
      t312 = 2d0*(z3*s12 - z1*s23)/(z3 + z1) + (z3 - z1)/(z3 + z1)*s13
      
      if (tens.eq.0) then
c---  tensor structure 0: -g^{mu nu}
         pggg = 4.5 + t123**2/(4.*s12**2) + t132**2/(4.*s13**2) + 
     -        t213**2/(4.*s12**2) + t231**2/(4.*s23**2) + 
     -        t312**2/(4.*s13**2) + t321**2/(4.*s23**2) + 
     -        (4*s123)/(s23*(1 - z1)) - 
     -        (2*s123)/(s12*(1 - z1)*z1) - 
     -        (2*s123)/(s13*(1 - z1)*z1) - 
     -        (4*s123)/(s23*(1 - z1)*z1) - 
     -        (8*s123*z1)/(s23*(1 - z1)) + 
     -        (4*s123)/(s13*(1 - z2)) + 
     -        s123**2/(s13*s23*(1 - z1)*(1 - z2)) + 
     -        s123**2/(s13*s23*z1*z2) + 
     -        s123/(s13*(1 - z1)*z1*z2) - 
     -        (2*s123)/(s12*(1 - z2)*z2) - 
     -        (4*s123)/(s13*(1 - z2)*z2) - 
     -        (2*s123)/(s23*(1 - z2)*z2) + 
     -        s123/(s23*z1*(1 - z2)*z2) + 
     -        (2*s123*z1)/(s23*(1 - z2)*z2) + 
     -        (2*s123*z2)/(s13*(1 - z1)*z1) - 
     -        (8*s123*z2)/(s13*(1 - z2)) - 
     -        (4*s123**2*z1*z2)/(s13*s23*(1 - z1)*(1 - z2)) + 
     -        (4*s123)/(s12*(1 - z3)) + 
     -        s123**2/(s12*s23*(1 - z1)*(1 - z3)) + 
     -        s123**2/(s12*s13*(1 - z2)*(1 - z3)) - 
     -        (2*s123**2*z1)/(s12*s13*(1 - z2)*(1 - z3)) + 
     -        (2*s123**2*z1**2)/(s12*s13*(1 - z2)*(1 - z3)) - 
     -        (2*s123**2*z2)/(s12*s23*(1 - z1)*(1 - z3)) + 
     -        (2*s123**2*z2**2)/(s12*s23*(1 - z1)*(1 - z3)) + 
     -        s123**2/(s12*s23*z1*z3) + 
     -        s123/(s12*(1 - z1)*z1*z3) + 
     -        s123**2/(s12*s13*z2*z3) - 
     -        (2*s123**2*z1)/(s12*s13*z2*z3) + 
     -        (2*s123**2*z1**2)/(s12*s13*z2*z3) + 
     -        s123/(s12*(1 - z2)*z2*z3) - 
     -        (2*s123**2*z2)/(s12*s23*z1*z3) + 
     -        (2*s123**2*z2**2)/(s12*s23*z1*z3) - 
     -        (4*s123)/(s12*(1 - z3)*z3) - 
     -        (2*s123)/(s13*(1 - z3)*z3) - 
     -        (2*s123)/(s23*(1 - z3)*z3) + 
     -        s123/(s23*z1*(1 - z3)*z3) + 
     -        (2*s123*z1)/(s23*(1 - z3)*z3) + 
     -        s123/(s13*z2*(1 - z3)*z3) + 
     -        (2*s123*z2)/(s13*(1 - z3)*z3) + 
     -        (2*s123*z3)/(s12*(1 - z1)*z1) - 
     -        (2*s123**2*z3)/(s13*s23*(1 - z1)*(1 - z2)) - 
     -        (2*s123**2*z3)/(s13*s23*z1*z2) + 
     -        (2*s123*z3)/(s12*(1 - z2)*z2) - 
     -        (8*s123*z3)/(s12*(1 - z3)) - 
     -        (4*s123**2*z1*z3)/(s12*s23*(1 - z1)*(1 - z3)) - 
     -        (4*s123**2*z2*z3)/(s12*s13*(1 - z2)*(1 - z3)) + 
     -        (2*s123**2*z3**2)/(s13*s23*(1 - z1)*(1 - z2)) + 
     -        (2*s123**2*z3**2)/(s13*s23*z1*z2)

         vec(:)=0d0
         
      elseif (tens.eq.1) then
c---  tensor structure 1: kt1^mu kt1^nu
         pggg = (2*s123*z1)/(s13*(1 - z2)) - 
     -        (2*s123**2*z1)/(s13*s23*(1 - z2)) - 
     -        (2*s123*z1**2)/(s13*(1 - z2)) + 
     -        (4*s123**2*z1**2)/(s13*s23*(1 - z2)) - 
     -        (2*s123**2*z1**3)/(s13*s23*(1 - z2)) - 
     -        (6*s123*z2)/s13 + (6*s123**2*z2)/(s13*s23) - 
     -        (6*s123**2*z1*z2)/(s13*s23) + 
     -        (2*s123**2*z2**2)/(s13*s23*(1 - z1)) + 
     -        (2*s123*z2**2)/(s13*(1 - z1)*z1) - 
     -        (2*s123**2*z2**2)/(s13*s23*(1 - z1)*z1) - 
     -        (2*s123**2*z2**3)/(s13*s23*(1 - z1)) - 
     -        (2*s123*z2**3)/(s13*(1 - z1)*z1) + 
     -        (2*s123**2*z2**3)/(s13*s23*(1 - z1)*z1) + 
     -        (2*s123*z1)/(s12*(1 - z3)) - 
     -        (2*s123**2*z1)/(s12*s23*(1 - z3)) - 
     -        (2*s123*z1**2)/(s12*(1 - z3)) + 
     -        (4*s123**2*z1**2)/(s12*s23*(1 - z3)) - 
     -        (2*s123**2*z1**3)/(s12*s23*(1 - z3)) - 
     -        (8*s123*z1*z2)/(s12*(1 - z3)) + 
     -        (8*s123**2*z1*z2)/(s12*s23*(1 - z3)) - 
     -        (8*s123**2*z1**2*z2)/(s12*s23*(1 - z3)) + 
     -        (4*s123*z1*z2)/(s12*(1 - z3)*z3) - 
     -        (4*s123**2*z1*z2)/(s12*s23*(1 - z3)*z3) + 
     -        (4*s123**2*z1**2*z2)/(s12*s23*(1 - z3)*z3) - 
     -        (6*s123*z3)/s12 + (6*s123**2*z3)/(s12*s23) - 
     -        (6*s123**2*z1*z3)/(s12*s23) - 
     -        (8*s123*z1*z3)/(s13*(1 - z2)) + 
     -        (8*s123**2*z1*z3)/(s13*s23*(1 - z2)) - 
     -        (8*s123**2*z1**2*z3)/(s13*s23*(1 - z2)) + 
     -        (4*s123*z1*z3)/(s13*(1 - z2)*z2) - 
     -        (4*s123**2*z1*z3)/(s13*s23*(1 - z2)*z2) + 
     -        (4*s123**2*z1**2*z3)/(s13*s23*(1 - z2)*z2) + 
     -        (2*s123**2*z3**2)/(s12*s23*(1 - z1)) + 
     -        (2*s123*z3**2)/(s12*(1 - z1)*z1) - 
     -        (2*s123**2*z3**2)/(s12*s23*(1 - z1)*z1) - 
     -        (2*s123**2*z3**3)/(s12*s23*(1 - z1)) - 
     -        (2*s123*z3**3)/(s12*(1 - z1)*z1) + 
     -        (2*s123**2*z3**3)/(s12*s23*(1 - z1)*z1)

         if (ip.le.2) then
c---  iff
            vec(:) = 0d0
            vec(1) = -p(rp,1)-p(sp,1)
            vec(2) = -p(rp,2)-p(sp,2)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
      
      elseif (tens.eq.2) then
c---  tensor structure 2: kt2^mu kt2^nu
         pggg = (-6*s123*z1)/s23 + (6*s123**2*z1)/(s13*s23) + 
     -        (2*s123**2*z1**2)/(s13*s23*(1 - z2)) - 
     -        (2*s123**2*z1**3)/(s13*s23*(1 - z2)) + 
     -        (2*s123*z1**2)/(s23*(1 - z2)*z2) - 
     -        (2*s123**2*z1**2)/(s13*s23*(1 - z2)*z2) - 
     -        (2*s123*z1**3)/(s23*(1 - z2)*z2) + 
     -        (2*s123**2*z1**3)/(s13*s23*(1 - z2)*z2) + 
     -        (2*s123*z2)/(s23*(1 - z1)) - 
     -        (2*s123**2*z2)/(s13*s23*(1 - z1)) - 
     -        (6*s123**2*z1*z2)/(s13*s23) - 
     -        (2*s123*z2**2)/(s23*(1 - z1)) + 
     -        (4*s123**2*z2**2)/(s13*s23*(1 - z1)) - 
     -        (2*s123**2*z2**3)/(s13*s23*(1 - z1)) + 
     -        (2*s123*z2)/(s12*(1 - z3)) - 
     -        (2*s123**2*z2)/(s12*s13*(1 - z3)) - 
     -        (8*s123*z1*z2)/(s12*(1 - z3)) + 
     -        (8*s123**2*z1*z2)/(s12*s13*(1 - z3)) - 
     -        (2*s123*z2**2)/(s12*(1 - z3)) + 
     -        (4*s123**2*z2**2)/(s12*s13*(1 - z3)) - 
     -        (8*s123**2*z1*z2**2)/(s12*s13*(1 - z3)) - 
     -        (2*s123**2*z2**3)/(s12*s13*(1 - z3)) + 
     -        (4*s123*z1*z2)/(s12*(1 - z3)*z3) - 
     -        (4*s123**2*z1*z2)/(s12*s13*(1 - z3)*z3) + 
     -        (4*s123**2*z1*z2**2)/(s12*s13*(1 - z3)*z3) - 
     -        (6*s123*z3)/s12 + (6*s123**2*z3)/(s12*s13) - 
     -        (6*s123**2*z2*z3)/(s12*s13) - 
     -        (8*s123*z2*z3)/(s23*(1 - z1)) + 
     -        (8*s123**2*z2*z3)/(s13*s23*(1 - z1)) + 
     -        (4*s123*z2*z3)/(s23*(1 - z1)*z1) - 
     -        (4*s123**2*z2*z3)/(s13*s23*(1 - z1)*z1) - 
     -        (8*s123**2*z2**2*z3)/(s13*s23*(1 - z1)) + 
     -        (4*s123**2*z2**2*z3)/(s13*s23*(1 - z1)*z1) + 
     -        (2*s123**2*z3**2)/(s12*s13*(1 - z2)) + 
     -        (2*s123*z3**2)/(s12*(1 - z2)*z2) - 
     -        (2*s123**2*z3**2)/(s12*s13*(1 - z2)*z2) - 
     -        (2*s123**2*z3**3)/(s12*s13*(1 - z2)) - 
     -        (2*s123*z3**3)/(s12*(1 - z2)*z2) + 
     -        (2*s123**2*z3**3)/(s12*s13*(1 - z2)*z2)

         if (ip.le.2) then
c---  iff
            vec(:) = 0d0
            vec(1) = (1d0-xs)*p(rp,1)+xr*p(sp,1)
            vec(2) = (1d0-xs)*p(rp,2)+xr*p(sp,2)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
         
      elseif (tens.eq.3) then
c---  tensor structure 3: kt3^mu kt3^nu
         pggg = (-6*s123*z1)/s23 + (6*s123**2*z1)/(s12*s23) - 
     -        (6*s123*z2)/s13 + (6*s123**2*z2)/(s12*s13) + 
     -        (2*s123**2*z1**2)/(s12*s23*(1 - z3)) - 
     -        (2*s123**2*z1**3)/(s12*s23*(1 - z3)) + 
     -        (2*s123**2*z2**2)/(s12*s13*(1 - z3)) - 
     -        (2*s123**2*z2**3)/(s12*s13*(1 - z3)) + 
     -        (2*s123*z1**2)/(s23*(1 - z3)*z3) - 
     -        (2*s123**2*z1**2)/(s12*s23*(1 - z3)*z3) - 
     -        (2*s123*z1**3)/(s23*(1 - z3)*z3) + 
     -        (2*s123**2*z1**3)/(s12*s23*(1 - z3)*z3) + 
     -        (2*s123*z2**2)/(s13*(1 - z3)*z3) - 
     -        (2*s123**2*z2**2)/(s12*s13*(1 - z3)*z3) - 
     -        (2*s123*z2**3)/(s13*(1 - z3)*z3) + 
     -        (2*s123**2*z2**3)/(s12*s13*(1 - z3)*z3) + 
     -        (2*s123*z3)/(s23*(1 - z1)) - 
     -        (2*s123**2*z3)/(s12*s23*(1 - z1)) - 
     -        (6*s123**2*z1*z3)/(s12*s23) + 
     -        (2*s123*z3)/(s13*(1 - z2)) - 
     -        (2*s123**2*z3)/(s12*s13*(1 - z2)) - 
     -        (8*s123*z1*z3)/(s13*(1 - z2)) + 
     -        (8*s123**2*z1*z3)/(s12*s13*(1 - z2)) + 
     -        (4*s123*z1*z3)/(s13*(1 - z2)*z2) - 
     -        (4*s123**2*z1*z3)/(s12*s13*(1 - z2)*z2) - 
     -        (6*s123**2*z2*z3)/(s12*s13) - 
     -        (8*s123*z2*z3)/(s23*(1 - z1)) + 
     -        (8*s123**2*z2*z3)/(s12*s23*(1 - z1)) + 
     -        (4*s123*z2*z3)/(s23*(1 - z1)*z1) - 
     -        (4*s123**2*z2*z3)/(s12*s23*(1 - z1)*z1) - 
     -        (2*s123*z3**2)/(s23*(1 - z1)) + 
     -        (4*s123**2*z3**2)/(s12*s23*(1 - z1)) - 
     -        (2*s123*z3**2)/(s13*(1 - z2)) + 
     -        (4*s123**2*z3**2)/(s12*s13*(1 - z2)) - 
     -        (8*s123**2*z1*z3**2)/(s12*s13*(1 - z2)) + 
     -        (4*s123**2*z1*z3**2)/(s12*s13*(1 - z2)*z2) - 
     -        (8*s123**2*z2*z3**2)/(s12*s23*(1 - z1)) + 
     -        (4*s123**2*z2*z3**2)/(s12*s23*(1 - z1)*z1) - 
     -        (2*s123**2*z3**3)/(s12*s23*(1 - z1)) - 
     -        (2*s123**2*z3**3)/(s12*s13*(1 - z2))

         if (ip.le.2) then
c---  iff
            vec(:) = 0d0
            vec(1) = (1d0-xr)*p(sp,1)+xs*p(rp,1)
            vec(2) = (1d0-xr)*p(sp,2)+xs*p(rp,2)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
         
      elseif (tens.eq.12) then
c---  tensor structure 4: kt12^mu kt12^nu
         pggg = (-6*s12*s123)/(s13*s23) + 
     -        (2*s12*s123*z1)/(s13*s23*(1 - z2)*z2) - 
     -        (2*s12*s123*z1**2)/(s13*s23*(1 - z2)*z2) + 
     -        (2*s12*s123*z2)/(s13*s23*(1 - z1)*z1) - 
     -        (2*s12*s123*z2**2)/(s13*s23*(1 - z1)*z1) - 
     -        (8*s123*z1*z2)/(s12*(1 - z3)*z3)
         
         if (ip.le.2) then
c---  iff
            vec(:) = 0d0
            vec(1) = xa**2/xr*p(rp,1)
            vec(2) = xa**2/xr*p(rp,1)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
         
      elseif (tens.eq.13) then
c---  tensor structure 5: kt13^mu kt13^nu
         pggg = (-6*s123*s13)/(s12*s23) + 
     -        (2*s123*s13*z1)/(s12*s23*(1 - z3)*z3) - 
     -        (2*s123*s13*z1**2)/(s12*s23*(1 - z3)*z3) + 
     -        (2*s123*s13*z3)/(s12*s23*(1 - z1)*z1) - 
     -        (8*s123*z1*z3)/(s13*(1 - z2)*z2) - 
     -        (2*s123*s13*z3**2)/(s12*s23*(1 - z1)*z1)
         
         if (ip.le.2) then
c---  iff
            vec(:) = 0d0
            vec(1) = xa**2/xs*p(sp,1)
            vec(2) = xa**2/xs*p(sp,1)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
         
      elseif (tens.eq.23) then
c---  tensor structure 6: kt23^mu kt23^nu
         pggg = (-6*s123*s23)/(s12*s13) + 
     -        (2*s123*s23*z2)/(s12*s13*(1 - z3)*z3) - 
     -        (2*s123*s23*z2**2)/(s12*s13*(1 - z3)*z3) + 
     -        (2*s123*s23*z3)/(s12*s13*(1 - z2)*z2) - 
     -        (8*s123*z2*z3)/(s23*(1 - z1)*z1) - 
     -        (2*s123*s23*z3**2)/(s12*s13*(1 - z2)*z2)

         if (ip.le.2) then
c---  iff
            vec(:) = 0d0
            vec(1) = xa**2*(p(sp,1)/xs - p(rp,1)/xr)
            vec(2) = xa**2*(p(sp,2)/xs - p(rp,2)/xr)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
      endif

c---  iff
!--   NEW
      if (ip.le.2) pggg = xa**2/xia**2*pggg
      
      return
      end

***********************************************************************
***********************************************************************
      subroutine pgg_pgg_eval(tens,p,ip,rp,jp,sp,pggpgg,vecr,vecs)

      implicit none
      include 'constants.f'
c      include 'sprods_com.f'
      double precision p(mxpart,4)
      
      integer tens,ip,rp,jp,sp
      double precision stot,srQ,ssQ,xa,xb,xr,xs,dot
      double precision yrQ,ysQ,yar,ybs,yrs,AA
      double precision pggpgg,vecr(4),vecs(4),vecsq

      if ((tens.ne.0).and.(tens.ne.1)
     .     .and.(tens.ne.2).and.(tens.ne.3)) then
         write(6,*)"Non-existent tensor structure ",tens
         write(6,*)"Allowed values are 0, 1, 2, 3"
         write(6,*)"Stopping."
         stop
      endif
      
      if ((ip.le.2).and.(jp.le.2)) then
c---  INITIAL-FINAL, INITIAL-FINAL
c         stot=s(1,2)
c         srQ=s(1,rp)+s(2,rp)
c         ssQ=s(1,sp)+s(2,sp)
c         yrs=s(rp,sp)/stot

         stot=two*dot(p,1,2)
         srQ=two*(dot(p,1,rp)+dot(p,2,rp))
         ssQ=two*(dot(p,1,sp)+dot(p,2,sp))
         yrs=two*dot(p,rp,sp)/stot
         
         xr = -srQ/stot
         xs = -ssQ/stot
         xa = 1d0 - xr
         xb = 1d0 - xs
                  
         if (tens.eq.0) then
c---  tensor structure 0: (-g^{mu nu})(-g^{al be})
            pggpgg = (4*(1 + xr**2)*(1 + xs**2))/(xr*xs)
            
            vecr(:)=0d0
            vecs(:)=0d0
            
         elseif (tens.eq.1) then
c---  tensor structure 1: ktr^mu ktr^nu (-g^{al be})
            pggpgg = (-8*xr*(1 + xs**2))/(xa**2*xs)
            
            vecr(:) = 0d0
            vecr(1) = -p(rp,1)
            vecr(2) = -p(rp,2)
            vecsq = vecr(1)**2 + vecr(2)**2
            vecr(:) = vecr(:)/dsqrt(vecsq)
            
            vecs = vecr
            
         elseif (tens.eq.2) then
c---  tensor structure 2: (-g^{mu nu}) kts^mu kts^nu
            pggpgg = (-8*(1 + xr**2)*xs)/(xb**2*xr)
            
            vecs(:) = 0d0
            vecs(1) = -p(sp,1)
            vecs(2) = -p(sp,2)
            vecsq = vecs(1)**2 + vecs(2)**2
            vecs(:) = vecs(:)/dsqrt(vecsq)
            
            vecr = vecs        
            
         elseif (tens.eq.3) then
c---  tensor structure 3: kt3^mu kt3^nu
            pggpgg = (16*xr*xs)/(xa**2*xb**2)
            
            vecr(:) = 0d0
            vecr(1) = -p(rp,1)
            vecr(2) = -p(rp,2)
            vecsq = vecr(1)**2 + vecr(2)**2
            vecr(:) = vecr(:)/dsqrt(vecsq)
            
            vecs(:) = 0d0
            vecs(1) = -p(sp,1)
            vecs(2) = -p(sp,2)
            vecsq = vecs(1)**2 + vecs(2)**2
            vecs(:) = vecs(:)/dsqrt(vecsq)
            
         endif

!--   NEW
         pggpgg = (xa*xb)**2/(1-xr-xs+yrs)**2*pggpgg
      
      endif
      return
      end
      
************************************************************************
************************************************************************


      
      
