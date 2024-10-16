      subroutine dipsRVC2_colorful(nd,p,ip,rp,sub,subv,msq,msqv,
     .     subrv,subrvv,msqrv,msqrvv,
     .     subr_born,subr_corr,subr_virt,subr_vcor)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single unresolved collinear subtraction terms in  *
*     colorful nnlo with momentum p                                    *
*     Automatically chooses subtraction kind (IF or FF)                *
*     Returns the subs in sub,subv and matrix elements in msq,msqv     *
*     nd labels the collinear configurations                           *
*     ip labels one daughter parton                                    *
*     rp labels the other daughter parton                              *
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
      include 'scale.f'
      include 'order.f'
      include 'part.f'
      include 'b0.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4),subv
      double precision subrv(4,-2:0),subrvv(-2:0)
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)
      double precision msqrv(-nf:nf,-nf:nf,-2:0),
     .                 msqrvv(-nf:nf,-nf:nf,-2:0)
        
      double precision sir,xa,zi,yir,yirab,yirhab,al,vec(4),vecsq,dot
      integer nd,ip,rp,indtmp,nu,j,k,ipt
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr

c---check ip.ne.rp
      if (ip.eq.rp) then
         write(6,*)"Indices ip and pr in dipolesubC2_colorful: ",ip,rp
         write(6,*)"Not distinct, stopping."
         stop
      endif
c---check ip, rp not 1 and 2
      if (((ip.eq.1).and.(rp.eq.2)).or.((ip.eq.2).and.(rp.eq.1))) then
         write(6,*)"Indices ip and pr in dipolesubC2_colorful: ",ip,rp
         write(6,*)"Both initial state, stopping."
         stop
      endif      

c---order, so only ip can be initial state
      if (ip.gt.rp) then
         indtmp = ip
         ip = rp
         rp = indtmp
      endif
      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      subrv(j,:)=0d0
      enddo
      subv=0d0
      subrvv(:)=0d0
      vec=0d0
      call zeromsq(msq,msqv)
      msqrv(:,:,:)=0d0
      msqrvv(:,:,:)=0d0
      incldip(nd)=.true.

c---  note convention: p(1) and p(2) have negative energy
      if (ip.le.2) then
         sir=-two*dot(p,ip,rp)
      else
         sir=two*dot(p,ip,rp)
      endif

***********************************************************************
*************************** INITIAL-FINAL *****************************
***********************************************************************
      if ((ip .le. 2) .and. (rp .gt. 2)) then
         
c---  note convention: p(1) and p(2) have negative energy
         xa = 1d0 + (dot(p,rp,1)+dot(p,rp,2))/dot(p,1,2)
         vec(:) = 0d0
         vec(1) = p(rp,1)
         vec(2) = p(rp,2)
         vecsq  = -vec(1)**2 - vec(2)**2

         vec(:) = vec(:)/dsqrt(-vecsq)
         
c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aii) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
         call transformC2_colorful(p,ptrans,ip,rp)
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
      
c$$$  c--- for "gamgam" process, gvec contribution is not written in
c$$$  c--- a canonical way; instead, pass emitted vector directly as "vec"
c$$$  if (case .eq. 'gamgam') then
c$$$  do nu=1,4
c$$$  vec(nu)=p(jp,nu)
c$$$  enddo
c$$$  endif
      
      call subr_born(ptrans,msq)
      call subr_corr(ptrans,vec,ip,msqv)
      
      sub(gg)=2d0*gsq/xa/sir*(2d0*(xa*(1d0-xa) + xa/(1d0-xa)))
      subv   =2d0*gsq/xa/sir*(4d0*(1d0-xa)/xa)

      if ((part .eq. 'virt').and.(order .eq. 2)) then
         call subr_virt(ptrans,msqrv)
         call subr_vcor(ptrans,vec,ip,msqrvv)
         
         subrv(gg, -2)=ca
         subrv(gg, -1)=ca*dlog(musq*xa**2/(sir*(1d0-xa))) + b0
         subrv(gg,  0)=ca*(half*dlog(musq*xa**2/(sir*(1d0-xa)))**2
     .        -pisq/3d0+2d0*dlog(xa)*dlog((1d0-xa)/xa))

         subrvv(-2)=subrv(gg, -2)
         subrvv(-1)=subrv(gg, -1)
         subrvv( 0)=subrv(gg, 0)
     .        + (ca-0d0*2d0*dfloat(nf)*tr)*xa**2/(1d0-xa)/6d0

         subrv(gg,  :)= -gsq/(8d0*pisq)*subrv(gg,:)*sub(gg)
         subrvv( :)= -gsq/(8d0*pisq)*subrvv(:)*subv
      endif
      
***********************************************************************
**************************** FINAL-FINAL ******************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (rp .gt. 2)) then

c---  note convention: p(1) and p(2) have negative energy
         zi = (dot(p,ip,1)+dot(p,ip,2))
     .        /(dot(p,ip,1)+dot(p,ip,2)+dot(p,rp,1)+dot(p,rp,2))

         yirab = -(dot(p,ip,1)+dot(p,ip,2)+dot(p,rp,1)+dot(p,rp,2))
     .        /dot(p,1,2)
         yir = dot(p,ip,rp)/dot(p,1,2)
         al = 0.5d0*(yirab - dsqrt(yirab**2 - 4d0*yir))
         yirhab = yirab+2d0*al
         
         vec(:) = (zi-yir/(al*yirab))*p(rp,:)
     .        -(1d0-zi-yir/(al*yirab))*p(ip,:)
     .        +yir/(al*yirhab)*(1d0-2d0*zi)
     .        *(p(ip,:)+p(rp,:)+al*(p(1,:)+p(2,:)))
         vecsq  = vec(4)**2 - vec(1)**2 - vec(2)**2 - vec(3)**2
         vec(:) = vec(:)/dsqrt(-vecsq)

c$$$C---Modification so that only close to singular subtracted
c$$$  if (variable .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
C---calculate the ptrans-momenta 

        call transformC2_colorful(p,ptrans,ip,rp)
        call storeptilde(nd,ptrans)

c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
c        if (incldip(nd) .eqv. .false.) return
        
c--- if using a dynamic scale, set that scale with dipole kinematics      
        if (dynamicscale) then
         call scaleset(initscale,initfacscale,ptrans)
         dipscale(nd)=facscale
        endif
      
        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ip,msqv)

c---  note: msqv includes a division by (-kt^2) already
c---  compenstate (-1) here
        sub(gg)=2d0*gsq/sir*(2d0*(zi/(1d0-zi)+(1d0-zi)/zi))
        subv   =2d0*gsq/sir*(4d0*zi*(1-zi))

      endif
      
      return
      end
      

      subroutine dipsRVS1_colorful(nd,p,rp,sub,msq,subrv,msqrv,
     .     subr_born,subr_virt)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020                                                        *
*     Calculates the single unresolved soft subtraction terms          *
*     including soft-collinear overlaps in colorful nnlo with          *
*     momentum p                                                       *
*     Returns the subs in sub and matrix elements in msq               *
*     nd labels the collinear configurations                           *
*     rp labels the soft parton                                        *
*     subr_born is the subroutine which call the born process          *
*     Note: treatment of color correlations is non-generich for now    *
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
      include 'scale.f'
      include 'order.f'
      include 'part.f'
      include 'b0.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(0:2),subrv(-2:0)
      double precision msq(-nf:nf,-nf:nf),msqrv(-nf:nf,-nf:nf,-2:0)

      double precision sjhkh,sjhr,skhr,dot,eik,colleik
      integer nd,rp,idec1,idec2,nu,j,k,nmax
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born

c---check rp.gt.2
      if (rp.le.2) then
         write(6,*)"Index pr in dipolesubS1_colorful: ",rp
         write(6,*)"Not in the final state, stopping."
         stop
      endif
      
C---  Initialize the dipoles to zero
      sub(:) = 0d0
      subrv(:) = 0d0
      msq(:,:) = 0d0
      msqrv(:,:,:) = 0d0
      incldip(nd)=.true.

c$$$  C---Modification so that only close to singular subtracted
c$$$  if (-varialbe .gt. aff) then
c$$$  incldip(nd)=.false.
c$$$  return
c$$$  endif
         
      call transformS1_colorful(p,ptrans,rp)
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

c--- hard coded for double real
      nmax = 5

      if (part .eq. 'virt') then
         nmax = 2
         if (order .eq. 2) call subr_virt(ptrans,msqrv)
      endif
c---  main loop over eikonal
      idec1 = 3
      idec2 = 4
      
      do j = 1,nmax
         if ((j.eq.idec1).or.(j.eq.idec2)) cycle
         sjhr = 2d0*(p(rp,4)*ptrans(j,4)
     .        -p(rp,1)*ptrans(j,1)
     .        -p(rp,2)*ptrans(j,2)
     .        -p(rp,3)*ptrans(j,3))
         do k = j+1,nmax
            if ((k.eq.idec1).or.(k.eq.idec2)) cycle
            
            sjhkh = 2d0*dot(ptrans,j,k)
            skhr = 2d0*(p(rp,4)*ptrans(k,4)
     .           -p(rp,1)*ptrans(k,1)
     .           -p(rp,2)*ptrans(k,2)
     .           -p(rp,3)*ptrans(k,3))
c---  note sum over j>k so extra factor of two
c---  soft
            eik = sjhkh/sjhr/skhr
            sub(0) = sub(0) + 2d0*eik
            
            if ((part .eq. 'virt').and.(order .eq. 2)) then
               subrv(-2) = subrv(-2) + two*ca*eik
               subrv(-1) = subrv(-1) + two*eik*(ca*dlog(musq*eik) + b0)
               subrv( 0) = subrv( 0) + ca*eik*(dlog(musq*eik)**2
     .              - 2d0*pisq/3d0)
            endif
         enddo
         
c---  initial-final soft-collinear and final-final collinear
c---  have the same expression when explicitly written out in
c---  terms of momenta p and ptrans

         colleik = 1d0/sjhr*
     .        ((p(1,4)+p(2,4))*ptrans(j,4)
     .        -(p(1,1)+p(2,1))*ptrans(j,1)
     .        -(p(1,2)+p(2,2))*ptrans(j,2)
     .        -(p(1,3)+p(2,3))*ptrans(j,3))
     .        /(dot(p,rp,1)+dot(p,rp,2))

c---  Ti^2/Ti.Tk = 1 here, so no extra factors required
c---  Ti.Tk added in _gs

c         sub = sub - 2d0*colleik
         sub(j) = -2d0*colleik
         
         if ((part .eq. 'virt').and.(order .eq. 2)) then
            subrv(-2) = subrv(-2) - two*ca*colleik
            subrv(-1) = subrv(-1)
     .           - two*colleik*(ca*dlog(musq*colleik) + b0)
            subrv( 0) = subrv( 0)
     .           -ca*colleik*(dlog(musq*colleik)**2 - 2d0*pisq/3d0)
            
         endif
      enddo      
      
      sub(:) = -2d0*gsq*sub(:)
      subrv(:) = gsq**2/(four*pisq)*subrv(:)
      
      return
      end         
