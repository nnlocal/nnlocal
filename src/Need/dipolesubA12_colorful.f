      subroutine dipsS1S2_colorful(nd,p,rp,sp,sub,msq,subr_born)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single soft limit of the double soft subtraction  *
*     terms including soft-collinear overlaps in colorful nnlo with    *
*     momentum p                                                       *
*     Returns the subs in sub and matrix elements in msq               *
*     nd labels the collinear configurations                           *
*     rp labels one of the soft partons                                *
*     sp labels the other soft parton - this is softer                 *
*     subr_born is the subroutine which call the born process          *
*     Note: treatment of color correlations is non-generic for now     *
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'      
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      double precision p(mxpart,4),ptrans(mxpart,4),phrans(mxpart,4),
     .     sub(4)
      double precision msq(-nf:nf,-nf:nf)

      double precision sitkt,sitrh,sktrh,sihkh,sihrh,skhrh,sihs,skhs,
     .     sjhlh,sjhs,slhs,srhs,sjhrh,xrh,xs,symfac,dot,sprod
      double precision tmpSsSrs(4),tmpSsCarsSrs(4),tmpCbsSsCSarsSrs(4),
     .     tmpSsCSarsSrs(4),tmpSsCarsCSarsSrs(4),tmpCrsSsSrs(4),
     .     tmpCrsSsCarsSrs(4),tmpCasSsSrs(4),tmpCasSsCarsSrs(4)

      
      integer nd,rp,rind,sp,idec1,idec2,nu,i,k,j,l,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born

      tmpSsSrs=0d0
      tmpSsCarsSrs=0d0
      tmpCbsSsCSarsSrs=0d0
      tmpSsCSarsSrs=0d0
      tmpSsCarsCSarsSrs=0d0
      tmpCrsSsSrs=0d0
      tmpCrsSsCarsSrs=0d0
      tmpCasSsSrs=0d0
      tmpCasSsCarsSrs=0d0
      
c---check rp.gt.2 and sp.gt.2
      if ((rp.le.2).or.(sp.le.2)) then
         write(6,*)"Index rp, sp in dipolesubS1S2_colorful: ",rp,sp
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

      call transformS1_colorful(p,phrans,sp)
      if (sp.lt.rp) rind = rp-1
      if (rp.lt.sp) rind = rp
      call transformS1_colorful(phrans,ptrans,rind)
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
      srhs = 2d0*sprod(phrans(rind,:),p(sp,:))
    
      do i = 1,nmax
         if ((i.eq.idec1).or.(i.eq.idec2)) cycle
         sihrh = 2d0*dot(phrans,i,rind)
         sihs = 2d0*sprod(phrans(i,:),p(sp,:))
         sitrh = 2d0*sprod(ptrans(i,:),phrans(rind,:))
         do k = i,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            sihkh = 2d0*dot(phrans,i,k)
            sitkt = 2d0*dot(ptrans,i,k)
            skhrh = 2d0*dot(phrans,k,rind)
            skhs = 2d0*sprod(phrans(k,:),p(sp,:))
            sktrh = 2d0*sprod(ptrans(k,:),phrans(rind,:))
            do j = 1,nmax
               if ((j.eq.idec1).or.(j.eq.idec2)) cycle
               sjhs = 2d0*sprod(phrans(j,:),p(sp,:))
               do l = j,nmax
                  if ((l.eq.idec1).or.(l.eq.idec2)) cycle
                  sjhlh = 2d0*dot(phrans,j,l)
                  slhs = 2d0*sprod(phrans(l,:),p(sp,:))
c---  abelian soft
c---  color factor is {(Ti.Tk),(Tj.Tl)}=2CA^2, only i.ne.k and j.ne.l conribute
c---  add 2 so full soft gets CA^2
c---  note sum over i <= k and j <= l so extra factor of two for each
                  tmpSsSrs(gg) = tmpSsSrs(gg)
     .                 + 4d0*sitkt/sitrh/sktrh*sjhlh/sjhs/slhs
               enddo

c---  initial-final soft-collinear limit of double soft
c---  overall color factor ca**2 correct only for f(i) = g
!---  NEW CANCELLATIONS IMPLEMENTED
               if ((i.le.2).and.(k.le.2).and.(j.le.2)) then
c$$$                  sjhrh = 2d0*dot(phrans,j,rind)
c$$$                  xrh = -(dot(phrans,rind,1)+dot(phrans,rind,2))
c$$$     .                 /dot(phrans,1,2)
c$$$c---  note sum over i <= k so extra factor of two for i.ne.k
c$$$                  tmpSsCSarsSrs(gg) = tmpSsCSarsSrs(gg)
c$$$     .                 + 4d0*sihkh/sihs/skhs/sjhrh/xrh
                  
c---  initial-final collinear-soft of double soft
                  xs = -sprod(p(sp,:),phrans(1,:)+phrans(2,:))
     .                 /dot(phrans,1,2)
c---  note sum over i <= k so extra factor of two for i.ne.k
                  tmpCasSsSrs(gg) = tmpCasSsSrs(gg)
     .                 + 4d0*sihkh/sihrh/skhrh/sjhs/xs
               endif               
            enddo
c---  non-abelian soft
c---  color factor is ca*Ti.Tk = ca**2 if i.eq.k, and -ca**2 if i.ne.k
c---  add appropriate (-1) so full soft gets CA^2
c---  note sum over i <= k so extra factor of two for i.ne.k
            symfac = 1d0
            if (i.ne.k) symfac = -2d0
            tmpSsSrs(gg) = tmpSsSrs(gg) - symfac*
     .           (sitkt/sitrh/sktrh*sihrh/sihs/srhs
     .           + sitkt/sitrh/sktrh*skhrh/skhs/srhs
     .           - sitkt/sitrh/sktrh*sihkh/sihs/skhs)

c---  initial-final soft-collinear limit of double soft
c---  overall color factor ca**2 correct only for f(i) = g
            if ((i.le.2).and.(k.le.2).and.(i.ne.k)) then
               xrh = -(dot(phrans,rind,1)+dot(phrans,rind,2))
     .              /dot(phrans,1,2)  
               xs = -sprod(p(sp,:),phrans(1,:)+phrans(2,:))
     .              /dot(phrans,1,2)
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$c---  symmetrized in i and k so i <= k summation gets both contributions
c$$$               tmpCbsSsCSarsSrs(gg) = tmpCbsSsCSarsSrs(gg)
c$$$     .              + 4d0/sihrh/skhs/xrh/xs
c$$$     .              + 4d0/skhrh/sihs/xrh/xs
               
c---  final-final collinear-soft of double soft
c---  note sum over i <= k so extra factor of two for i.ne.k
c---  note: (1-zs)/zs = xrh/xs in this case
c$$$               tmpCrsSsSrs(gg) = tmpCrsSsSrs(gg)
c$$$     .              - 4d0*sihkh/sihrh/skhrh*xrh/xs/srhs               
            endif            
         enddo
      
c---  initial-final-final triple collinear limit of double soft
c---  overall color factor ca**2 correct only for f(i) = g
         if (i.le.2) then
            xrh = -(dot(phrans,rind,1)+dot(phrans,rind,2))
     .           /dot(phrans,1,2)
            xs = -sprod(p(sp,:),phrans(1,:)+phrans(2,:))
     .           /dot(phrans,1,2)
c---  note convention: p(1) and p(2) have negative energy
            tmpSsCarsSrs(gg) = tmpSsCarsSrs(gg) + 2d0/(sihs*srhs*xrh)
     .           + 2d0/(sihrh*srhs*xs) - 2d0/(sihrh*sihs*xrh*xs)

c---  triple coll. of soft-coll. of double soft
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$            tmpSsCarsCSarsSrs(gg) = tmpSsCarsCSarsSrs(gg)
c$$$     .           + 4d0/sihs/sihrh/xrh/xs

c---  final-final collinear-soft of triple coll. of double soft
c---  note: zr/zs = xrh/xs in this case
c$$$            tmpCrsSsCarsSrs(gg) = tmpCrsSsCarsSrs(gg)
c$$$     .           - 4d0/sihrh/srhs/xrh*xrh/xs
            
c---  initial-final collinear-soft of triple coll. of double soft
c---  note convention: p(1) and p(2) have negative energy
            tmpCasSsCarsSrs(gg) = tmpCasSsCarsSrs(gg)
     .           + 4d0/sihrh/xrh/sihs/xs
         endif 
      enddo      

!---  NEW CANCELLATIONS IMPLEMENTED
c$$$      sub =
c$$$     .     tmpSsSrs
c$$$     .     + tmpSsCarsSrs
c$$$     .     + tmpSsCSarsSrs
c$$$     .     + tmpCbsSsCSarsSrs
c$$$     .     + tmpSsCarsCSarsSrs
c$$$     .     + tmpCrsSsSrs
c$$$     .     + tmpCrsSsCarsSrs
c$$$     .     + tmpCasSsSrs
c$$$     .     + tmpCasSsCarsSrs

      sub =
     .     tmpSsSrs
     .     + tmpSsCarsSrs
     .     + tmpCasSsSrs
     .     + tmpCasSsCarsSrs
            
      sub = -(2d0*gsq)**2*sub
      
      return
      end         


      subroutine dipsS1CS2_colorful(nd,p,ip,rp,sp,sub,subv,
     .     msq,msqv,subr_born,subr_corr)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single soft limit of the double unresolved        *
*     collinear-soft subtraction term in colorful nnlo with momentum p *
*     Returns the subs in sub, subv and matrix elements in msq, msqv   *
*     nd labels the collinear configurations                           *
*     ip labels one of the daughter partons of the coll. splitting     *
*     rp labels the other daughter parton                              *
*     sp labels the soft parton                                        *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *      
*     Note: treatment of color correlations is non-generic for now     *
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
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

      double precision sahrh,xah,vec(4),vecsq,sjhkh,sjhs,skhs,dot,sprod
      integer nd,ip,rp,sp,rind,idec1,idec2,nu,j,k,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubS1CS2_colorful: "
     .        ,ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.sp
      if (ip.eq.sp) then
         write(6,*)"Indices ip and sp in dipolesubS1CS2_colorful: "
     .        ,ip,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubS1CS2_colorful: "
     .        ,rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubS1CS2_colorful: "
     .        ,ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check ip, sp not 1 and 2
      if (((ip.eq.1).and.(sp.eq.2)).or.((ip.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices ip and sp in dipolesubS1CS2_colorful: "
     .        ,ip,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, sp not 1 and 2
      if (((rp.eq.1).and.(sp.eq.2)).or.((rp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices rp and sp in dipolesubS1CS2_colorful: "
     .        ,rp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif            

c---order, so only ip can be initial state
      if ((ip.gt.rp).or.(ip.gt.sp)) then
         write(6,*)"Wrong order of indices in dipolesubS1CS2_colorful."
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

***********************************************************************
************************ INITIAL-FINAL-FINAL **************************
***********************************************************************
      if ((ip .le. 2) .and. (rp .gt. 2) .and. (sp .gt. 2)) then

      
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif

         call transformS1_colorful(p,phrans,sp)
         if (sp.lt.rp) rind = rp-1
         if (rp.lt.sp) rind = rp
         call transformC2_colorful(phrans,ptrans,ip,rind)
         call storeptilde(nd,ptrans)         
         
c---  note convention: p(1) and p(2) have negative energy
         sahrh=two*dot(phrans,ip,rind)

c---  note convention: p(1) and p(2) have negative energy
         xah = 1d0 + (dot(phrans,rind,1)+dot(phrans,rind,2))
     .        /dot(phrans,1,2)
         vec(:) = 0d0
         vec(1) = phrans(rind,1)
         vec(2) = phrans(rind,2)
         vecsq  = vec(1)**2 + vec(2)**2
         
         vec(:) = vec(:)/dsqrt(vecsq)          
         
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
         sub(gg)=2d0*gsq/xah/sahrh*(2d0*(xah*(1d0-xah) + xah/(1d0-xah)))
         subv   =2d0*gsq/xah/sahrh*(4d0*(1d0-xah)/xah)      

      endif
      
c---  main loop over eikonal
      nmax = 2
      idec1 = 3
      idec2 = 4
      
      do j = 1,nmax
         if ((j.eq.idec1).or.(j.eq.idec2)) cycle
         sjhs = 2d0*sprod(phrans(j,:),p(sp,:))
         do k = j+1,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            
            sjhkh = 2d0*dot(phrans,j,k)
            skhs = 2d0*sprod(phrans(k,:),p(sp,:))
c---  note sum over j>k so extra factor of two
c---  soft
            eik = eik + 2d0*sjhkh/sjhs/skhs
            
         enddo         
      enddo

!---  NEW DEFINITION
      eik=0d0
      sub = -2d0*gsq*eik*sub
      subv = -2d0*gsq*eik*subv
      
      return
      end         
      

      subroutine dipsC2C3_colorful(nd,p,ip,rp,jp,sub,subv,msq,msqv,
     .     subr_born,subr_corr)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single collinear limit of the tripe               *
*     collinear subtraction terms in colorful nnlo with momentum p     *
*     Automatically chooses subtraction kind (IF,IFF, FF,IFF, FF,FFF)  *
*     Currently only FF,IFF and IF,IFF subtraction C_rs C_ars and      *
*     C_ar C_ars implemented                                           *      
*     Returns the subs in sub,subv and matrix elements in msq,msqv     *
*     nd labels the collinear configurations                           *
*     ip labels one of the single collinear partons                    *
*     rp labels the other single collinear parton                      *
*     jp labels the third parton in the triple collinear configuration *
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
      double precision p(mxpart,4),phrans(mxpart,4),ptrans(mxpart,4),
     .     sub(4),subv(2)
      double precision msq(-nf:nf,-nf:nf),msqv(2,-nf:nf,-nf:nf)

      double precision sir,sirhjh,fac,dot
      double precision vec0(4),vec1(4),vec2(4),
     .     pggpggg0,pggpggg1,pggpggg2
      integer nd,ip,rp,jp,jind,nu,i,j,k
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubC2C3_colorful: ",ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.jp
      if (ip.eq.jp) then
         write(6,*)"Indices ip and jp in dipolesubC2C3_colorful: ",ip,jp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.jp
      if (rp.eq.jp) then
         write(6,*)"Indices rp and jp in dipolesubC2C3_colorful: ",rp,jp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubC2C3_colorful: ",ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check ip, jp not 1 and 2
      if (((ip.eq.1).and.(jp.eq.2)).or.((ip.eq.2).and.(jp.eq.1))) then
         write(6,*)"Indices ip and jp in dipolesubC2C3_colorful: ",ip,jp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, jp not 1 and 2
      if (((rp.eq.1).and.(jp.eq.2)).or.((rp.eq.2).and.(jp.eq.1))) then
         write(6,*)"Indices rp and jp in dipolesubC2C3_colorful: ",rp,jp
         write(6,*)"Both initial state, stopping."
         stop
      endif

c---order, so in the single collinear pair ip < rp
      if (ip.gt.rp) then
         write(6,*)"Wrong order of indices in dipolesubC2C3_colorful."
         write(6,*)"Correct is ip < rp."
         write(6,*)"Currently ip, rp = ",ip, rp
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
c      call zeromsq(msq,msqv)
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0d0
            do i=1,2
               msqv(i,j,k)=0d0
            enddo
         enddo      
      enddo      

      incldip(nd)=.true.

      sir=two*(dot(p,ip,rp))
      
***********************************************************************
*************** INITIAL-FINAL of INITIAL-FINAL-FINAL ******************
***********************************************************************
      if ((ip .le. 2) .and. (rp .gt. 2) .and. (jp .gt. 2)) then

c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC2_colorful(p,phrans,ip,rp)
         if (rp.lt.jp) jind = jp-1
         if (jp.lt.rp) jind = jp
c---  note: here ip < jp, rp
         sirhjh=two*dot(phrans,ip,jind)
         call transformC2_colorful(phrans,ptrans,ip,jind)
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
            
      call pgg_pggg_eval(0,p,phrans,ip,rp,jp,pggpggg0,vec0)
      sub(gg) = pggpggg0
      call subr_born(ptrans,msq)
    
      call pgg_pggg_eval(1,p,phrans,ip,rp,jp,pggpggg1,vec1)
      subv(1) = pggpggg1
      call subr_corr(ptrans,vec1,ip,msqv(1,:,:))

      call pgg_pggg_eval(2,p,phrans,ip,rp,jp,pggpggg2,vec2)
      subv(2) = pggpggg2
      call subr_corr(ptrans,vec2,ip,msqv(2,:,:))

      fac = -(2d0*gsq)**2/sir/sirhjh      
      sub = fac*sub
c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
      subv = -fac*subv

***********************************************************************
**************** FINAL-FINAL of INITIAL-FINAL-FINAL *******************
***********************************************************************
      elseif ((ip .gt. 2) .and. (rp .gt. 2) .and. (jp .le. 2)) then
         
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC2_colorful(p,phrans,ip,rp)
         sirhjh=two*dot(phrans,jp,ip)
c--- note: here jp < ip < rp
         call transformC2_colorful(phrans,ptrans,jp,ip)
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
            
      call pgg_pggg_eval(0,p,phrans,ip,rp,jp,pggpggg0,vec0)
      sub(gg) = pggpggg0
      call subr_born(ptrans,msq)
    
      call pgg_pggg_eval(1,p,phrans,ip,rp,jp,pggpggg1,vec1)
      subv(1) = pggpggg1
      call subr_corr(ptrans,vec1,ip,msqv(1,:,:))

      call pgg_pggg_eval(2,p,phrans,ip,rp,jp,pggpggg2,vec2)
      subv(2) = pggpggg2
      call subr_corr(ptrans,vec2,ip,msqv(2,:,:))

      fac = -(2d0*gsq)**2/sir/sirhjh      
      sub = fac*sub
c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
      subv = -fac*subv
      
***********************************************************************
************************* FINAL-FINAL-FINAL ***************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (rp .gt. 2) .and. (jp .gt. 2)) then
         write(6,*)"C^{FF}C^{FFF} not yet implemented in dipA12_col."
         write(6,*)"Stopping."
         stop
      endif
      
      return
      end
      

      subroutine dipsC2C22_colorful(nd,p,ip,rp,jp,sp,sub,subv,subvv,
     .     msq,msqv,msqvv,subr_born,subr_corr,subr_corr_2)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single collinear limit of the double collinear    *
*     subtraction terms in colorful nnlo with momentum p               *
*     Automatically chooses subtraction kind (IF,IF, etc.)             *
*     Currently only IF,IF subtracgtion C_ar C_ar,bs implemented       *      
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
      double precision p(mxpart,4),ptrans(mxpart,4),phrans(mxpart,4),
     .     sub(4),subv(2),subvv(4)
      double precision msq(-nf:nf,-nf:nf),msqv(2,-nf:nf,-nf:nf),
     .     msqvv(-nf:nf,-nf:nf)

      double precision sir,sjhsh,fac,dot,sprod
      double precision vec0a(4),vec0b(4),vec1a(4),vec1b(4),
     .     vec2a(4),vec2b(4),vec3a(4),vec3b(4),
     .     pggpggpgg0,pggpggpgg1,pggpggpgg2,pggpggpgg3
      integer nd,ip,rp,jp,sp,sind,nu,i,j,k
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

***********************************************************************
******************** INITIAL-FINAL, INITIAL-FINAL *********************
***********************************************************************
      if ((ip.le.2).and.(rp.gt.2).and.(jp.le.2).and.(sp.gt.2)) then
                  
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif

         call transformC2_colorful(p,phrans,ip,rp)
         if (rp.lt.sp) sind = sp-1
         if (sp.lt.rp) sind = sp
         sjhsh=two*dot(phrans,jp,sind)
         call transformC2_colorful(phrans,ptrans,jp,sind)
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
            
      call pgg_pgg_pgg_eval(0,p,phrans,ip,rp,jp,sp,
     .     pggpggpgg0,vec0a,vec0b)
      sub(gg) = pggpggpgg0
      call subr_born(ptrans,msq)

      call pgg_pgg_pgg_eval(1,p,phrans,ip,rp,jp,sp,
     .     pggpggpgg1,vec1a,vec1b)
      subv(1) = pggpggpgg1
      call subr_corr(ptrans,vec1a,ip,msqv(1,:,:))

      call pgg_pgg_pgg_eval(2,p,phrans,ip,rp,jp,sp,
     .     pggpggpgg2,vec2a,vec2b)
      subv(2) = pggpggpgg2
      call subr_corr(ptrans,vec2b,jp,msqv(2,:,:))
      
      call pgg_pgg_pgg_eval(3,p,phrans,ip,rp,jp,sp,
     .     pggpggpgg3,vec3a,vec3b)
      subvv(gg) = pggpggpgg3
      call subr_corr_2(ptrans,vec3a,vec3b,ip,jp,msqvv(:,:))

      fac = -(2d0*gsq)**2/sir/sjhsh
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


      subroutine dipsC2CS2_colorful(nd,p,ip,rp,sp,sub,subv,
     .     msq,msqv,subr_born,subr_corr)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single collinear limit of the double unresolved   *
*     collinear-soft subtraction term including collinear overlaps in  *
*     colorful nnlo with momentum p                                    *
*     Returns the subs in sub, subv and matrix elements in msq, msqv   *
*     nd labels the collinear configurations                           *
*     ip labels one of the daughter partons of the coll. splitting     *
*     rp labels the other daughter parton                              *
*     sp labels the soft parton                                        *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *      
*     Note: treatment of color correlations is non-generic for now     *
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
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

      double precision sir,xa,vec(4),vecsq,sjhkh,sjhsh,skhsh,
     .     sihsh,xsh,dot,sprod
      integer nd,ip,rp,sp,sind,idec1,idec2,nu,j,k,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubC2CS2_colorful: "
     .        ,ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.sp
      if (ip.eq.sp) then
         write(6,*)"Indices ip and sp in dipolesubC2CS2_colorful: "
     .        ,ip,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubC2CS2_colorful: "
     .        ,rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubC2CS2_colorful: "
     .        ,ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check ip, sp not 1 and 2
      if (((ip.eq.1).and.(sp.eq.2)).or.((ip.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices ip and sp in dipolesubC2CS2_colorful: "
     .        ,ip,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, sp not 1 and 2
      if (((rp.eq.1).and.(sp.eq.2)).or.((rp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices rp and sp in dipolesubC2CS2_colorful: "
     .        ,rp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif            

c---order, so only ip can be initial state
      if ((ip.gt.rp).or.(ip.gt.sp)) then
         write(6,*)"Wrong order of indices in dipolesubC2CS2_colorful."
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
         vecsq  = vec(1)**2 + vec(2)**2

         vec(:) = vec(:)/dsqrt(vecsq)   
      
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC2_colorful(p,phrans,ip,rp)
         if (rp.lt.sp) sind = sp-1
         if (sp.lt.rp) sind = sp
         call transformS1_colorful(phrans,ptrans,sind)
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
         sjhsh = 2d0*dot(phrans,sind,j)
         do k = j+1,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            
            sjhkh = 2d0*dot(phrans,j,k)
            skhsh = 2d0*dot(phrans,sind,k)
c---  note sum over j>k so extra factor of two
c---  soft
            eik = eik + 2d0*sjhkh/sjhsh/skhsh
            
         enddo         
      enddo

c---  initial-final-final triple collinear limit of soft-collinear
      sihsh = 2d0*dot(phrans,ip,sind)
c---  note convention: p(1) and p(2) have negative energy
      xsh = -(dot(phrans,sind,1)+dot(phrans,sind,2))/dot(phrans,1,2)

c---  sign exchanged so overall color factor is ca**2
      eik = eik + 2d0/sihsh/xsh

c---  initial-final, initial-final double collinear limit of soft-collinear
      if (ip.eq.1) then
         sjhsh=2d0*dot(phrans,2,sind)
c---  note convention: p(1) and p(2) have negative energy
         xsh = -(sprod(phrans(sind,:),p(1,:)+p(2,:)))
     .        /(sprod(phrans(2,:),p(1,:)+p(2,:)))
      elseif (ip.eq.2) then
         sjhsh=2d0*dot(phrans,1,sind)
c---  note convention: p(1) and p(2) have negative energy
         xsh = -(sprod(phrans(sind,:),p(1,:)+p(2,:)))
     .        /(sprod(phrans(1,:),p(1,:)+p(2,:)))
      endif      
c---  sign exchanged so overall color factor is ca**2
      eik = eik + 2d0/sjhsh/xsh

!---  NEW DEFINITION
      eik=0d0
      sub = 2d0*gsq*eik*sub
      subv = 2d0*gsq*eik*subv
      
      return
      end         

      
      subroutine dipsC2S2_colorful(nd,p,rp,sp,sub,msq,subr_born)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single collinear limit of the double soft         *
*     subtraction terms in colorful nnlo with momentum p               *
*     Returns the subs in sub and matrix elements in msq               *
*     nd labels the collinear configurations                           *
*     rp labels one soft parton                                        *
*     sp labels the other soft parton in the double soft configuration *
*     subr_born is the subroutine which call the born process          *
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
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
     .     sub(4)
      double precision msq(-nf:nf,-nf:nf)

      double precision srs,sitkt,sitrsh,sktrsh,sitktrs,sktktrs,
     .     zr,zs,xrs,yrsab,yrshab,yrs,al,vec(4),vecsq,dot,sprod,
     .     fac,symfac
      integer nd,rp,sp,nu,i,j,k,nmax,idec1,idec2
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born

c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubC2S2_colorful: ",rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.jp
      if ((rp.le.2).or.(sp.le.2)) then
         write(6,*)"Indices rp and sp in dipolesubC2S2_colorful: ",rp,sp
         write(6,*)"Some not in the final state, stopping."
         stop
      endif
c---order, so in the single collinear pair rp < sp
      if (rp.gt.sp) then
         write(6,*)"Wrong order of indices in dipolesubC2S2_colorful."
         write(6,*)"Correct is rp < sp."
         write(6,*)"Currently rp, sp = ",rp, sp
         write(6,*)"Stopping."
         stop
      endif

      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      msq = 0d0
      vec=0d0
      incldip(nd)=.true.

c      srs=s(rp,sp)
      srs=two*dot(p,rp,sp)
c---  note convention: p(1) and p(2) have negative energy
c      zr = (s(rp,1)+s(rp,2))/(s(rp,1)+s(rp,2)+s(sp,1)+s(sp,2))
      zr = (dot(p,rp,1)+dot(p,rp,2))/
     .     (dot(p,rp,1)+dot(p,rp,2)+dot(p,sp,1)+dot(p,sp,2))
      zs = 1d0-zr
c      yrsab = -(s(rp,1)+s(rp,2)+s(sp,1)+s(sp,2))/s(1,2)
      yrsab = -(dot(p,rp,1)+dot(p,rp,2)+dot(p,sp,1)+dot(p,sp,2))/
     .         dot(p,1,2)
c      yrs = s(rp,sp)/s(1,2) 
      yrs = dot(p,rp,sp)/dot(p,1,2) 
      al = 0.5d0*(yrsab - dsqrt(yrsab**2 - 4d0*yrs))
      yrshab = yrsab+2d0*al
      
      vec(:) = (zr-yrs/(al*yrsab))*p(sp,:)
     .     -(1d0-zr-yrs/(al*yrsab))*p(rp,:)
     .     +yrs/(al*yrshab)*(1d0-2d0*zr)
     .     *(p(rp,:)+p(sp,:)+al*(p(1,:)+p(2,:)))
      vecsq  = vec(4)**2 - vec(1)**2 - vec(2)**2 - vec(3)**2
      vec(:) = vec(:)/dsqrt(-vecsq)

c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
      call transformC2_colorful(p,phrans,rp,sp)
c---  note: here rp < sp
      call transformS1_colorful(phrans,ptrans,rp)
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

      do i = 1,nmax
         if ((i.eq.idec1).or.(i.eq.idec2)) cycle         
         sitrsh = 2d0*sprod(ptrans(i,:),phrans(rp,:))
         sitktrs = 2d0*sprod(ptrans(i,:),vec(:))
         
         do k = i,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            sitkt = 2d0*dot(ptrans,i,k)
            sktrsh = 2d0*sprod(ptrans(k,:),phrans(rp,:))
            sktktrs = 2d0*sprod(ptrans(k,:),vec(:))

c---  color factor is ca*Ti.Tk = ca**2 if i.eq.k, and -ca**2 if i.ne.k
c---  add appropriate (-1) so full soft gets CA^2
c---  note sum over i <= k so extra factor of two for i.ne.k
c---  note vec already includes 1/(-kt^2) so add (-1)
            symfac = 1d0
            if (i.ne.k) symfac = -2d0            
            sub(gg) = sub(gg) + symfac * (2d0*(zr/zs + zs/zr)
     .           * sitkt/sitrsh/sktrsh
     .           -2d0*zr*zs*sitktrs*sktktrs/sitrsh/sktrsh)
         enddo
c---  note vec already includes 1/(-kt^2) so add (-1)
c---  note factor of srs added at the end
         xrs = sprod(phrans(rp,:),ptrans(1,:)+ptrans(2,:))
     .        /sprod(ptrans(i,:),ptrans(1,:)+ptrans(2,:))
         
         sub(gg) = sub(gg) + 2d0/sitrsh*
     .        (2d0/xrs*(zr/zs + zs/zr)
     .        +zr*zs*sitktrs**2/sitrsh) 
         
      enddo

      fac = (2d0*gsq)**2/srs      
      sub = fac*sub
      
      return
      end


      subroutine dipsS1C3_colorful(nd,p,ip,rp,sp,sub,subv,
     .     msq,msqv,subr_born,subr_corr)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single soft limit of the triple collinear         *
*     subtraction term including collinear overlaps in colorful nnlo   *
*     with momentum p                                                  *
*     Returns the subs in sub, subv and matrix elements in msq, msqv   *
*     nd labels the collinear configurations                           *
*     ip labels one of the daughter partons of the coll. splitting     *
*     rp labels the other daughter parton                              *
*     sp labels the third daughter and this is the soft parton         *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *      
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
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
     .     sub(4),subv
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)

      double precision sahrh,xah,vec(4),vecsq,
     .     satrh,sats,srhs,xat,xrt,xst,sbs,xs,psggg,dot,sprod
      integer nd,ip,rp,rind,sp,jp,sind,idec1,idec2,nu,j,k,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and rp in dipolesubS1C3_colorful: "
     .        ,ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip.ne.sp
      if (ip.eq.sp) then
         write(6,*)"Indices ip and sp in dipolesubS1C3_colorful: "
     .        ,ip,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check rp.ne.sp
      if (rp.eq.sp) then
         write(6,*)"Indices rp and sp in dipolesubS1C3_colorful: "
     .        ,rp,sp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubS1C3_colorful: "
     .        ,ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check ip, sp not 1 and 2
      if (((ip.eq.1).and.(sp.eq.2)).or.((ip.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices ip and sp in dipolesubS1C3_colorful: "
     .        ,ip,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif
c---check rp, sp not 1 and 2
      if (((rp.eq.1).and.(sp.eq.2)).or.((rp.eq.2).and.(sp.eq.1))) then
         write(6,*)"Indices rp and sp in dipolesubS1C3_colorful: "
     .        ,rp,sp
         write(6,*)"Both initial state, stopping."
         stop
      endif            

c---order, so only ip can be initial state
      if ((ip.gt.rp).or.(ip.gt.sp)) then
         write(6,*)"Wrong order of indices in dipolesubS1C3_colorful."
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
      psggg=0d0
      call zeromsq(msq,msqv)
      incldip(nd)=.true.

***********************************************************************
************************ INITIAL-FINAL-FINAL **************************
***********************************************************************
      if ((ip .le. 2) .and. (rp .gt. 2) .and. (sp .gt. 2)) then
      
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformS1_colorful(p,phrans,sp)
         if (sp.lt.rp) rind = rp-1
         if (rp.lt.sp) rind = rp
         call transformC2_colorful(phrans,ptrans,ip,rind)
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
         
c---  note convention: p(1) and p(2) have negative energy
         sahrh = 2d0*dot(phrans,ip,rind)
         xah = 1d0 + (dot(phrans,rind,1)+dot(phrans,rind,2))
     .        /dot(phrans,1,2)
         
         vec(:) = 0d0
         vec(1) = phrans(rind,1)
         vec(2) = phrans(rind,2)
         vecsq  = vec(1)**2 + vec(2)**2
         vec(:) = vec(:)/dsqrt(vecsq)

         call subr_born(ptrans,msq)
         call subr_corr(ptrans,vec,ip,msqv)

c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
         sub(gg)=2d0*gsq/xah/sahrh*(2d0*(xah*(1d0-xah) + xah/(1d0-xah)))
         subv   =2d0*gsq/xah/sahrh*(4d0*(1d0-xah)/xah)

c---  note convention: p(1) and p(2) have negative energy
         satrh = 2d0*sprod(ptrans(ip,:),phrans(rind,:))
         sats = 2d0*sprod(ptrans(ip,:),p(sp,:))
         srhs = 2d0*sprod(phrans(rind,:),p(sp,:))
         xrt = sprod(phrans(rind,:),p(1,:)+p(2,:))
     .        /sprod(ptrans(ip,:),p(1,:)+p(2,:))
         xst = sprod(p(sp,:),p(1,:)+p(2,:))
     .        /sprod(ptrans(ip,:),p(1,:)+p(2,:))
         xat = 1d0 - xrt - xst

         psggg = (satrh/srhs/sats + 1d0/sats/xst + xrt/srhs/xst)
         
c---  final-final collinear-soft limit of the triple collinear
c---  sign exchanged so overall color factor is ca**2
c---  note: zr/zs = xr/xs in this case
         psggg = psggg - 2d0/srhs*xrt/xst

c---  soft limit of the triple collinear of double soft-collinear
c---  sign exchanged so overall color factor is ca**2
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$         psggg = psggg - 2d0/sats/xst
         
c---  initial-final collinear-soft limit of the triple collinear
c---  sign exchanged so overall color factor is ca**2
         psggg = psggg - 2d0/sats/xst

c---  initial-final collinear-soft limit of the double soft-collinear
c---  sign exchanged so overall color factor is ca**2
         if (ip.eq.1) jp=2
         if (ip.eq.2) jp=1
c         sbs = s(jp,sp)
c         xs = (s(sp,1)+s(sp,2))/(s(jp,1)+s(jp,2))
         sbs = two*dot(p,jp,sp)
         xs = (dot(p,sp,1)+dot(p,sp,2))/(dot(p,jp,1)+dot(p,jp,2))

!---  NEW CANCELLATIONS IMPLEMENTED
c$$$      psggg = psggg - 2d0/sbs/xs
      endif
      
      sub(gg) = 2d0*gsq*psggg*sub(gg)
      subv = 2d0*gsq*psggg*subv
      
      return
      end         
      

***********************************************************************
***********************************************************************
      subroutine pgg_pggg_eval(tens,p,ph,ip,rp,jp,pggpggg,vec)

      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      double precision p(mxpart,4),ph(mxpart,4)
      
      integer tens,ip,rp,jp,jind
      double precision sirhjh,siQ,srQ,xa,zi,zr,xah,xiah,zirh,zjh,
     .     yirab,yir,al,yirhab,sjhktir,dot,sprod
      double precision pggpggg,vec(4),vecsq

      if ((tens.ne.0).and.(tens.ne.1).and.(tens.ne.2)) then
         write(6,*)"Non-existent tensor structure ",tens
         write(6,*)"Allowed values are 0, 1, 2"
         write(6,*)"Stopping."
         stop
      endif
            
      if ((ip.gt.2).and.(jp.le.2)) then
c---  FINAL-FINAL OF INITIAL-FINAL-FINAL
         sirhjh = two*dot(ph,ip,jp)
         
         siQ=two*(dot(p,1,ip)+dot(p,2,ip))
         srQ=two*(dot(p,1,rp)+dot(p,2,rp))
c         siQ=s(1,ip)+s(2,ip)
c         srQ=s(1,rp)+s(2,rp)
         zi = siQ/(siQ+srQ)
         zr = 1d0-zi
c---  note convention: p(1) and p(2) have negative energy
         xah = 1d0 + (dot(ph,ip,1)+dot(ph,ip,2))/dot(ph,1,2)
         zjh  = 1d0/xah
         zirh = 1d0 - zjh

c         yirab = -(siQ+srQ)/s(1,2)
c         yir = s(ip,rp)/s(1,2)
         yirab = -(siQ+srQ)/two/dot(p,1,2)
         yir = dot(p,ip,rp)/dot(p,1,2)
         al = 0.5d0*(yirab - dsqrt(yirab**2 - 4d0*yir))
         yirhab = yirab+2d0*al
         
         vec(:) = (zi-yir/(al*yirab))*p(rp,:)
     .        -(1d0-zi-yir/(al*yirab))*p(ip,:)
     .        +yir/(al*yirhab)*(1d0-2d0*zi)
     .        *(p(ip,:)+p(rp,:)+al*(p(1,:)+p(2,:)))
         vecsq  = vec(4)**2 - vec(1)**2 - vec(2)**2 - vec(3)**2
         vec(:) = vec(:)/dsqrt(-vecsq)
         
         sjhktir = 2d0*sprod(ph(jp,:),vec(:))
         
      elseif ((ip.le.2).and.(jp.gt.2)) then
c---  INITIAL-FINAL OF INITIAL-FINAL-FINAL
         if (rp.lt.jp) jind = jp-1
         if (jp.lt.rp) jind = jp
         sirhjh = two*dot(ph,ip,jind)      

c         srQ=s(1,rp)+s(2,rp)
         srQ=two*(dot(p,1,rp)+dot(p,2,rp))
c---  note convention: p(1) and p(2) have negative energy         
c         xa = 1d0 + srQ/s(1,2)
         xa = 1d0 + srQ/two/dot(p,1,2)
         zi = 1d0/xa
         zr = 1d0-zi
c---  note convention: p(1) and p(2) have negative energy         
!---  NEW
         xiah = 1d0 + (dot(ph,jind,1)+dot(ph,jind,2))/dot(ph,1,2)
         xah = 1d0 + sprod(ph(jind,:),p(1,:)+p(2,:))
     .        /sprod(ph(ip,:),p(1,:)+p(2,:))
         zirh  = 1d0/xah
         zjh = 1d0 - zirh
         
         vec(:) = 0d0
         vec(1) = p(rp,1)
         vec(2) = p(rp,2)
         vecsq  = vec(1)**2 + vec(2)**2
         vec(:) = vec(:)/dsqrt(vecsq)
         
         sjhktir = 2d0*sprod(ph(jind,:),vec(:))         

      endif
      
      if (tens.eq.0) then
c---  tensor structure 0: -g^{mu nu}
c---  note: (-ktir^2) already included with sjhktir, so correct for (-1)
         pggpggg = 4d0*(zi/zr+zr/zi)*(zjh/zirh+zirh/zjh)
     .        + 2d0*zi*zr*sjhktir**2/sirhjh
         
         vec(:)=0d0         
      elseif (tens.eq.1) then         
c---  tensor structure 1: kt_ir^mu kt_ir^nu
         pggpggg = -8d0*zi*zr*zirh/zjh

         if ((ip.gt.2).and.(jp.le.2)) then
c--- ff-iff
c            vec(:) = vec(:) - two*sprod(vec(:),p(jp,:))/s(1,2)
            vec(:) = vec(:) - sprod(vec(:),p(jp,:))/dot(p,1,2)
     .           *(p(1,:)+p(2,:))
            vecsq = vec(4)**2 - vec(1)**2 - vec(2)**2 - vec(3)**2
            vec(:) = vec(:)/dsqrt(-vecsq)
         elseif (ip.le.2) then
c---  if-iff
            vec(:) = vec(:)
         endif
      elseif (tens.eq.2) then
c---  tensor structure 2: kt2^mu kt2^nu
         pggpggg = -8d0*zirh*zjh*(zi/zr+zr/zi+zi*zr)
         
         if ((ip.gt.2).and.(jp.le.2)) then
c---  ff-iff
            vec(:) = 0d0
            vec(1) = ph(ip,1)
            vec(2) = ph(ip,2)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         elseif (ip.le.2) then
c---  if-iff
            vec(:) = 0d0
            vec(1) = ph(jind,1)
            vec(2) = ph(jind,2)
            vecsq = vec(1)**2 + vec(2)**2
            vec(:) = vec(:)/dsqrt(vecsq)
         endif
      endif

c---  if-iff
!---  NEW
      if (ip.le.2) pggpggg = (xah/xiah)**2*pggpggg
      
      return
      end


***********************************************************************
***********************************************************************
      subroutine pgg_pgg_pgg_eval(tens,p,ph,ip,rp,jp,sp,
     .     pggpggpgg,vecr,vecs)

      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      double precision p(mxpart,4),ph(mxpart,4)
      
      integer tens,ip,rp,jp,sp,sind
      double precision stot,srQ,ssQ,xa,xbh,xr,xsh,xish,sprod
      double precision pggpggpgg,vecr(4),vecs(4),vecsq,dot

      if ((tens.ne.0).and.(tens.ne.1)
     .     .and.(tens.ne.2).and.(tens.ne.3)) then
         write(6,*)"Non-existent tensor structure ",tens
         write(6,*)"Allowed values are 0, 1, 2, 3"
         write(6,*)"Stopping."
         stop
      endif
      
      if ((ip.le.2).and.(jp.le.2)) then
c---  INITIAL-FINAL, INITIAL-FINAL
         if (rp.lt.sp) sind = sp-1
         if (sp.lt.rp) sind = sp
         
c         stot=s(1,2)
c         srQ=s(1,rp)+s(2,rp)
         stot=two*dot(p,1,2)
         srQ=two*(dot(p,1,rp)+dot(p,2,rp))
c---  note convention: p(1) and p(2) have negative energy
         xr = -srQ/stot
         xsh = -(sprod(ph(sind,:),p(1,:)+p(2,:)))
     .        /(sprod(ph(jp,:),p(1,:)+p(2,:)))
!---  NEW
         xish = -(sprod(ph(sind,:),ph(1,:)+ph(2,:)))
     .        /(sprod(ph(jp,:),ph(1,:)+ph(2,:)))
         
         xa = 1d0-xr
         xbh = 1d0-xsh
         
         if (tens.eq.0) then
c---  tensor structure 0: (-g^{mu nu})(-g^{al be})
            pggpggpgg = (4*(1 + xr**2)*(1 + xsh**2))/(xr*xsh)
            
            vecr(:)=0d0
            vecs(:)=0d0
            
         elseif (tens.eq.1) then
c---  tensor structure 1: ktr^mu ktr^nu (-g^{al be})
            pggpggpgg = (-8*xr*(1 + xsh**2))/(xa**2*xsh)
            
            vecr(:) = 0d0
            vecr(1) = -p(rp,1)
            vecr(2) = -p(rp,2)
            vecsq = vecr(1)**2 + vecr(2)**2
            vecr(:) = vecr(:)/dsqrt(vecsq)
            
            vecs = vecr
            
         elseif (tens.eq.2) then
c---  tensor structure 2: (-g^{mu nu}) kts^mu kts^nu
            pggpggpgg = (-8*(1 + xr**2)*xsh)/(xbh**2*xr)
            
            vecs(:) = 0d0
            vecs(1) = -ph(sind,1)
            vecs(2) = -ph(sind,2)
            vecsq = vecs(1)**2 + vecs(2)**2
            vecs(:) = vecs(:)/dsqrt(vecsq)
            
            vecr = vecs        
            
         elseif (tens.eq.3) then
c---  tensor structure 3: kt3^mu kt3^nu
            pggpggpgg = (16*xr*xsh)/(xa**2*xbh**2)
            
            vecr(:) = 0d0
            vecr(1) = -p(rp,1)
            vecr(2) = -p(rp,2)
            vecsq = vecr(1)**2 + vecr(2)**2
            vecr(:) = vecr(:)/dsqrt(vecsq)
            
            vecs(:) = 0d0
            vecs(1) = -ph(sind,1)
            vecs(2) = -ph(sind,2)
            vecsq = vecs(1)**2 + vecs(2)**2
            vecs(:) = vecs(:)/dsqrt(vecsq)
            
         endif

!---  NEW
         pggpggpgg = pggpggpgg*(1d0-xsh)**2/(1d0-xish)**2
         
      endif
      return
      end
      
************************************************************************
************************************************************************


      
      
      
