      subroutine gg_hgaga_v_cf(p,msqv)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf,-2:2),
     . p(mxpart,4),dot,xl12,xlf
      integer j,k
      do j=-nf,nf
         do k=-nf,nf
            msqv(j,k,:)=0d0
         enddo
      enddo
      call gg_hgamgam(p,msq)
      xl12=log(two*dot(p,1,2)/musq)
      xlf=dble(nf)
      xlf=zip

      msqv(0,0,-2)=-6d0

      msqv(0,0,-1)=(-33+18*xl12+2*xlf)/3d0

      msqv(0,0, 0)=11-3*xl12**2+3*pisq

      msqv(0,0, 1)=xl12**3-3*pisq*xl12+xlf*pisq/18d0
     1     +3*(4*zeta3-11*pisq/36d0-2)

      msqv(0,0, 2)=(-9*xl12**4+33*pisq+54*pisq*xl12**2
     1     -3*(44*zeta3+3*pisq**2/5d0+288-72)
     2     +8*xlf*zeta3-216*xl12*(2*zeta3-1))/36d0

      do j=-2,2
         msqv(0,0,j)=msqv(0,0,j)*msq(0,0)*ason2pi
      enddo
      
      return
      end


      subroutine gg_hgaga_vv_cf(p,msqvv)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'masses.f'
      double precision msq(-nf:nf,-nf:nf),msqvv(-nf:nf,-nf:nf,-4:0),
     . p(mxpart,4),dot,xl12,xlf,lt
      integer j,k
      do j=-nf,nf
         do k=-nf,nf
            msqvv(j,k,:)=0d0
         enddo
      enddo
      call gg_hgamgam(p,msq)
      xl12=log(two*dot(p,1,2)/musq)
c      xlf=dble(nf)
      xlf=zip
      lt=log(musq/mt**2)
      
      msqvv(0,0,-4)=18d0
      
      msqvv(0,0,-3)=(363-144*xl12-22*xlf)/4d0
      
      msqvv(0,0,-2)=8+36*xl12**2-69*pisq/4d0
     1     +5*xl12*(2*xlf-33)/2d0-61*xlf/6d0
     2     +xlf**2/3d0

      msqvv(0,0,-1)=-428/3d0-24*xl12**3+xl12**2*(33-2*xlf)
     1     - 22*pisq+4*xlf*(26+3*pisq)/9d0-135*zeta3/2d0
     2     + xl12*(597-10*xlf+207*pisq)/6d0

      msqvv(0,0, 0)=lt*(57+16*xlf)/6d0+(-399+10*xlf-207*pisq)
     1     *xl12**2/6d0+(-33+2*xlf)*xl12**3/6d0 + 12*xl12**4
     2     +(476490+88245*pisq+60*xlf**2*pisq+3699*pisq**2
     3     -10*xlf*(3671+408*pisq-360*zeta3)-178200*zeta3)/1080d0
     5     +xl12*(-2112+xlf*(152-6*pisq)+99*pisq+4860*zeta3)/36d0
     6     -18
      
      do j=-4,0
         msqvv(0,0,j)=msqvv(0,0,j)*msq(0,0)*ason2pi**2
      enddo

      return
      end
      
