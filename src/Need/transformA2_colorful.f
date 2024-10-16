      subroutine transformC3_colorful(p,q,ip,rp,sp)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020.                                                       *
*     Triple collinear momentum mappings c_ars, c_irs in colorful nnlo *
*     Correct branch chosen automatically                              *
*     Currently only IFF mapping c_ars implemented                     *
*     p is original, q is transformed momenta                          *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'npart.f'
      double precision p(mxpart,4),q(mxpart,4),xa,xb,yirab,yir,al,dot,
     .     k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      double precision yar,ybr,yrQ,yas,ybs,ysQ,yrs

      integer ip,rp,sp,jp,j,nu,ipart

      do j=1,npart+2
         do nu=1,4
            q(j,nu)=0d0
         enddo
      enddo
      
      if ((ip .le. 2) .and. (rp .gt. 2) .and. (sp .gt. 2)) then
c---  initial-final-final
         if (ip.eq.1) then
            jp = 2
         elseif (ip.eq.2) then
            jp = 1
         endif
         
c---  colorful prescription
c---  note convention: p(1) and p(2) have negative energy
!---  NEW DEFINITION
c$$$         xa = 1d0 + (dot(p,rp,1)+dot(p,rp,2) + dot(p,sp,1)+dot(p,sp,2)
c$$$     .        + dot(p,rp,sp))/dot(p,1,2)
c$$$         xb = 1d0
         yar = -dot(p,ip,rp)/dot(p,1,2)
         ybr = -dot(p,jp,rp)/dot(p,1,2)
         yrQ = yar+ybr
         yas = -dot(p,ip,sp)/dot(p,1,2)
         ybs = -dot(p,jp,sp)/dot(p,1,2)
         ysQ = yas+ybs
         yrs = dot(p,rp,sp)/dot(p,1,2)
         xa = sqrt((1d0-ybr-ybs)/(1d0-yar-yas)*(1d0-yrQ-ysQ+yrs))
         xb = sqrt((1d0-yar-yas)/(1d0-ybr-ybs)*(1d0-yrQ-ysQ+yrs))
         
         do nu=1,4
            q(ip,nu) = xa*p(ip,nu)
            q(jp,nu) = xb*p(jp,nu)
            k(nu)  = -p(ip,nu)-p(jp,nu)-p(rp,nu)-p(sp,nu)
            kt(nu) = -xa*p(ip,nu)-xb*p(jp,nu)
            ks(nu) = k(nu)+kt(nu)
         enddo
         
         kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
         ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2
         
         ipart=3
         do j=3,npart+2
            if ((j .eq. rp) .or. (j .eq. sp)) then
               go to 19
            else
             kDp(j)=k(4)*p(j,4)-k(1)*p(j,1)-k(2)*p(j,2)-k(3)*p(j,3)
             ksDp(j)=ks(4)*p(j,4)-ks(1)*p(j,1)-ks(2)*p(j,2)-ks(3)*p(j,3)
               do nu=1,4
                  q(ipart,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks
     .                 +two*kDp(j)*kt(nu)/kDk
               enddo
               ipart=ipart+1
            endif
 19         continue
         enddo
         return
      endif
      end


      subroutine transformC22_colorful(p,q,ip,rp,jp,sp)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020.                                                       *
*     Double collinear momentum mappings c_ar,bs, etc in colorful nnlo *
*     Correct branch chosen automatically                              *
*     Currently only IF,IF mapping c_ar,bs implemented                 *
*     p is original, q is transformed momenta                          *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'npart.f'
      double precision p(mxpart,4),q(mxpart,4),xa,xb,dot,
     .     yar,yrab,ybs,ysab,yrs,A,
     .     k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      double precision ybr,yrQ,yas,ysQ
      integer ip,rp,sp,jp,j,nu,ipart

      do j=1,npart+2
         do nu=1,4
            q(j,nu)=0d0
         enddo
      enddo
      
      if ((ip .le. 2).and.(jp.le.2).and.(rp .gt. 2).and.(sp.gt.2)) then
c---  initial-final, initial-final         
c---  colorful prescription
c---  note convention: p(1) and p(2) have negative energy
!---  NEW DEFINITION
c$$$         yar = -dot(p,ip,rp)/dot(p,1,2)
c$$$         yrab = -(dot(p,ip,rp)+dot(p,jp,rp))/dot(p,1,2)
c$$$         ybs = -dot(p,jp,sp)/dot(p,1,2)
c$$$         ysab = -(dot(p,ip,sp)+dot(p,jp,sp))/dot(p,1,2)
c$$$         yrs = dot(p,rp,sp)/dot(p,1,2)
c$$$         
c$$$         A = 1d0/(2d0*yar*ybs)*((1d0-yrab)*ybs+(1d0-ysab)*yar
c$$$     .        - dsqrt(((1d0-yrab)*ybs+(1d0-ysab)*yar)**2
c$$$     .        - 4d0*yar*ybs*(yrab*ysab-yrs)))
c$$$         
c$$$         xa = 1d0 - yrab - yar * A
c$$$         xb = 1d0 - ysab - ybs * A
         yar = -dot(p,ip,rp)/dot(p,1,2)
         ybr = -dot(p,jp,rp)/dot(p,1,2)
         yrQ = yar+ybr
         yas = -dot(p,ip,sp)/dot(p,1,2)
         ybs = -dot(p,jp,sp)/dot(p,1,2)
         ysQ = yas+ybs
         yrs = dot(p,rp,sp)/dot(p,1,2)
         xa = sqrt((1d0-ybr-ybs)/(1d0-yar-yas)*(1d0-yrQ-ysQ+yrs))
         xb = sqrt((1d0-yar-yas)/(1d0-ybr-ybs)*(1d0-yrQ-ysQ+yrs))
         
         do nu=1,4
            q(ip,nu) = xa*p(ip,nu)
            q(jp,nu) = xb*p(jp,nu)
            k(nu)  = -p(ip,nu)-p(jp,nu)-p(rp,nu)-p(sp,nu)
            kt(nu) = -xa*p(ip,nu)-xb*p(jp,nu)
            ks(nu) = k(nu)+kt(nu)
         enddo
         
         kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
         ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2
         
         ipart=3
         do j=3,npart+2
            if ((j .eq. rp) .or. (j .eq. sp)) then
               go to 19
            else
             kDp(j)=k(4)*p(j,4)-k(1)*p(j,1)-k(2)*p(j,2)-k(3)*p(j,3)
             ksDp(j)=ks(4)*p(j,4)-ks(1)*p(j,1)-ks(2)*p(j,2)-ks(3)*p(j,3)
               do nu=1,4
                  q(ipart,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks
     .                 +two*kDp(j)*kt(nu)/kDk
               enddo
               ipart=ipart+1
            endif
 19         continue
         enddo
         return
      endif
      end
      
      
      subroutine transformS2_colorful(p,q,rp,sp)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020.                                                       *
*     Double soft momentum mapping s_rs in colorful nnlo               *
*     p is original, q is transformed momenta                          *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'npart.f'
      double precision p(mxpart,4),q(mxpart,4),lam,dot,
     . k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      integer rp,sp,j,nu,ipart
      
      do j=1,npart+2
         do nu=1,4
            q(j,nu)=0d0
         enddo
      enddo
               
c---  colorful prescription
c---  note convention: p(1) and p(2) have negative energy
      lam = dsqrt(1d0 + (dot(p,rp,1)+dot(p,rp,2)
     .     + dot(p,sp,1)+dot(p,sp,2) + dot(p,rp,sp))/dot(p,1,2))
      
      do nu=1,4
         q(1,nu) = lam*p(1,nu)
         q(2,nu) = lam*p(2,nu)
         k(nu)  = -p(1,nu)-p(2,nu)-p(rp,nu)-p(sp,nu)
         kt(nu) = -lam*p(1,nu)-lam*p(2,nu)
         ks(nu) = k(nu)+kt(nu)
      enddo
      
      kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
      ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2
      
      ipart=3
      do j=3,npart+2
         if ((j .eq. rp) .or. (j .eq. sp)) then
            go to 23
         else
            kDp(j)=k(4)*p(j,4)-k(1)*p(j,1)-k(2)*p(j,2)-k(3)*p(j,3)
            ksDp(j)=ks(4)*p(j,4)-ks(1)*p(j,1)-ks(2)*p(j,2)-ks(3)*p(j,3)
            do nu=1,4
               q(ipart,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks
     .              +two*kDp(j)*kt(nu)/kDk
            enddo
            ipart=ipart+1
         endif
 23      continue
      enddo
      return
      end      
      
