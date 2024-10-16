      subroutine gg_hgamgam_v_cf(p,msqv)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf,-2:0),
     . p(mxpart,4),dot,xl12
      integer j,k
      
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k,:)=0d0
      enddo
      enddo

      call gg_hgamgam(p,msq)
      xl12=log(two*dot(p,1,2)/musq)

      msqv(0,0,-2)=ason2pi*xn*2d0*(-1d0)*msq(0,0)

c---  Note nf=0 in beta0 for now!
      msqv(0,0,-1)=ason2pi*xn*2d0*(xl12
     .  -11d0/6d0)*msq(0,0)
c     .  -((11d0-two*dble(nf)/xn))/6d0)*msq(0,0)

      msqv(0,0,0)=ason2pi*xn*2d0*(-0.5d0*xl12**2+11d0/6d0+0.5d0*pisq)
     .     *msq(0,0)

      return
      end
     
