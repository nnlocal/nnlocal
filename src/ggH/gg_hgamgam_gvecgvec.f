      subroutine gg_hgamgam_gvecgvec(p,n,m,in,jm,msq)
      implicit none
      include 'constants.f'
C  in is the label of the momentum contracted with n
C  jm is the label of the momentum contracted with m
      integer j,k,in,jm
      double precision msq(-nf:nf,-nf:nf),msqt(-nf:nf,-nf:nf)
      double precision n(4),m(4),nDm,p(mxpart,4)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      nDm=n(4)*m(4)-n(3)*m(3)-n(2)*m(2)-n(1)*m(1)
      call gg_hgamgam(p,msqt)

      msq(0,0)=0.5d0*nDm**2*msqt(0,0)      

      return
      end


      
