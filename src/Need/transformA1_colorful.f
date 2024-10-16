      subroutine transformC2_colorful(p,q,ip,rp)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020.                                                       *
*     Single collinear momentum mappings c_ar, c_ir in colorful nnlo   *
*     Correct branch chosen automatically                              *
*     p is original, q is transformed momenta                          *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'npart.f'
      double precision p(mxpart,4),q(mxpart,4),xa,xb,yirab,yir,al,dot,
     .     k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      double precision yar,ybr,yrQ
      integer ip,rp,jp,j,nu,ipart

      do j=1,npart+2
         do nu=1,4
            q(j,nu)=0d0
         enddo
      enddo
      
      if ((ip .le. 2) .and. (rp .gt. 2)) then
c---  initial-final
         if (ip.eq.1) then
            jp = 2
         elseif (ip.eq.2) then
            jp = 1
         endif
         
c---  colorful prescription
c---  note convention: p(1) and p(2) have negative energy
!---  NEW DEFINITION
c$$$         xa = 1d0 + (dot(p,rp,1)+dot(p,rp,2))/dot(p,1,2)
c$$$         xb = 1d0
         yar = -dot(p,ip,rp)/dot(p,1,2)
         ybr = -dot(p,jp,rp)/dot(p,1,2)
         yrQ = yar+ybr
         xa = sqrt((1d0-ybr)/(1d0-yar)*(1d0-yrQ))
         xb = sqrt((1d0-yar)/(1d0-ybr)*(1d0-yrQ))
         
         do nu=1,4
            q(ip,nu) = xa*p(ip,nu)
            q(jp,nu) = xb*p(jp,nu)
            k(nu)  = -p(ip,nu)-p(jp,nu)-p(rp,nu)
            kt(nu) = -xa*p(ip,nu)-xb*p(jp,nu)
            ks(nu) = k(nu)+kt(nu)
         enddo
         
         kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
         ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2
         
         ipart=3
         do j=3,npart+2
            if (j .eq. rp) then
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
      elseif ((ip .gt. 2) .and. (rp .gt. 2)) then
c---  final-final

c---  note convention: p(1) and p(2) have negative energy
         yirab = -(dot(p,ip,1)+dot(p,ip,2)+dot(p,rp,1)+dot(p,rp,2))
     .        /dot(p,1,2)
         yir = dot(p,ip,rp)/dot(p,1,2)
         al = 0.5d0*(yirab - dsqrt(yirab**2 - 4d0*yir))

         do nu=1,4
            q(1,nu) = (1d0-al)*p(1,nu)
            q(2,nu) = (1d0-al)*p(2,nu)
         enddo

         ipart=3
         do j=3,npart+2
            do nu=1,4
               if (j.eq.ip) then
                  q(ipart,nu)=p(ip,nu)+p(rp,nu)+al*(p(1,nu)+p(2,nu))
               elseif (j.eq.rp) then
                  goto 21
               else
                  q(ipart,nu)=p(j,nu)
               endif
            enddo
            ipart=ipart+1
 21         continue
         enddo
         return
      endif
      end

      
      subroutine transformS1_colorful(p,q,rp)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020.                                                       *
*     Single soft momentum mapping s_r in colorful nnlo                *
*     p is original, q is transformed momenta                          *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'npart.f'
      double precision p(mxpart,4),q(mxpart,4),lam,dot,
     . k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      integer rp,j,nu,ipart
      
      do j=1,npart+2
         do nu=1,4
            q(j,nu)=0d0
         enddo
      enddo
               
c---  colorful prescription
c---  note convention: p(1) and p(2) have negative energy
      lam = dsqrt(1d0 + (dot(p,rp,1)+dot(p,rp,2))/dot(p,1,2))
      
      do nu=1,4
         q(1,nu) = lam*p(1,nu)
         q(2,nu) = lam*p(2,nu)
         k(nu)  = -p(1,nu)-p(2,nu)-p(rp,nu)
         kt(nu) = -lam*p(1,nu)-lam*p(2,nu)
         ks(nu) = k(nu)+kt(nu)
      enddo
      
      kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
      ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2
      
      ipart=3
      do j=3,npart+2
         if (j .eq. rp) then
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
      
