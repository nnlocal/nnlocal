      subroutine gg_hgamgam_gs_cf(p,msq,msqv,msqx)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c----for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H --> (b(p3)+b~(p4))+g(p5)
c---
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'order.f'
      include 'epinv.f'
      include 'epinv2.f'

      integer j,k,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
     .     ,msqx(maxd,-nf:nf,-nf:nf,0:2)
     .     ,msqv(maxd,-nf:nf,-nf:nf)
      double precision msq15(-nf:nf,-nf:nf),msq25(-nf:nf,-nf:nf),
     .     msq15_v(-nf:nf,-nf:nf),msq25_v(-nf:nf,-nf:nf),
     .     msq5(-nf:nf,-nf:nf),
     .     msqrv15(-nf:nf,-nf:nf,-2:0),msqrv25(-nf:nf,-nf:nf,-2:0),
     .     msqrv15_v(-nf:nf,-nf:nf,-2:0),msqrv25_v(-nf:nf,-nf:nf,-2:0),
     .     msqrv5(-nf:nf,-nf:nf,-2:0),
     .     sub15(4),sub25(4),sub15_v,sub25_v,sub5(0:2),
     .     subrv15(4,-2:0),subrv25(4,-2:0),
     .     subrv15_v(-2:0),subrv25_v(-2:0),subrv5(-2:0)
      external gg_hgamgam,gg_hgamgam_gvec,
     .     gg_hgamgam_v_cf,gg_hgamgam_vgvec_cf
      integer iglue
      parameter(iglue=5)
      ndmax=3
      ndmax1=ndmax

      
c---  A1 terms
c---  calculate single collinear subtractions
      call dipsRVC2_colorful(1,p,1,iglue,sub15,sub15_v,msq15,msq15_v,
     .     subrv15,subrv15_v,msqrv15,msqrv15_v,
     .     gg_hgamgam,gg_hgamgam_gvec,
     .     gg_hgamgam_v_cf,gg_hgamgam_vgvec_cf)

      
      call dipsRVC2_colorful(2,p,2,iglue,sub25,sub25_v,msq25,msq25_v,
     .     subrv25,subrv25_v,msqrv25,msqrv25_v,
     .     gg_hgamgam,gg_hgamgam_gvec,
     .     gg_hgamgam_v_cf,gg_hgamgam_vgvec_cf)
      

      
c---  calculate single soft subtractions      
      call dipsRVS1_colorful(3,p,iglue,sub5,msq5,subrv5,msqrv5,
     .     gg_hgamgam,gg_hgamgam_v_cf)

      
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0d0
      msqv(nd,j,k)=0d0
      msqx(nd,j,k,:)=0d0
      enddo
      enddo
      enddo

      
      do j=-nf,nf
      do k=-nf,nf
         
      if     ((j .ne. 0) .and. (k .eq. 0)) then
c         msq(1,j,k)=2d0*cf
c     .   *(msq17_2(0,0)*sub17_2(gq)+msq17_2v(0,0)*sub17_2v)
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
c         msq(2,j,k)=2d0*cf
c     .   *(msq27_1(0,0)*sub27_1(gq)+msq27_1v(0,0)*sub27_1v)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
         msq(1,j,k)=ca*(sub15(gg)*msq15(j,k)+sub15_v*msq15_v(j,k))
         msq(2,j,k)=ca*(sub25(gg)*msq25(j,k)+sub25_v*msq25_v(j,k))
         msq(3,j,k)=-ca*(sub5(0)+sub5(1)+sub5(2))*msq5(j,k)
         
c--- debug: we are adding only the subtraction to the correction
         if (abs(order) .eq. 2) then

c---  DEBUG: return ep^(-2) part
c---  DEBUG: return ep^(-1) part

c---  ep^(-2)
            msqv(1,j,k)=ca*(subrv15(gg,-2)*msq15(j,k)
     .           + subrv15_v(-2)*msq15_v(j,k)
     .           + sub15(gg)*msqrv15(j,k,-2)
     .           + sub15_v*msqrv15_v(j,k,-2))
            msqv(2,j,k)=ca*(subrv25(gg,-2)*msq25(j,k)
     .           + subrv25_v(-2)*msq25_v(j,k)
     .           + sub25(gg)*msqrv25(j,k,-2)
     .           + sub25_v*msqrv25_v(j,k,-2))
            msqv(3,j,k)=-ca*(subrv5(-2)*msq5(j,k)
     .           + (sub5(0)+sub5(1)+sub5(2))*msqrv5(j,k,-2))
            
c---  ep^(-1)
            msqv(1,j,k)= epinv2*msqv(1,j,k) +
     .           ca*(subrv15(gg,-1)*msq15(j,k)
     .           + subrv15_v(-1)*msq15_v(j,k)
     .           + sub15(gg)*msqrv15(j,k,-1)
     .           + sub15_v*msqrv15_v(j,k,-1))
            msqv(2,j,k)= epinv2*msqv(2,j,k) +
     .           ca*(subrv25(gg,-1)*msq25(j,k)
     .           + subrv25_v(-1)*msq25_v(j,k)
     .           + sub25(gg)*msqrv25(j,k,-1)
     .           + sub25_v*msqrv25_v(j,k,-1))
            msqv(3,j,k)= epinv2*msqv(3,j,k) +
     .           -ca*(subrv5(-1)*msq5(j,k)
     .           + (sub5(0)+sub5(1)+sub5(2))*msqrv5(j,k,-1))            

c---  ep^(0)
            msqv(1,j,k)= epinv*msqv(1,j,k) +
     .           ca*(subrv15(gg,0)*msq15(j,k)
     .           + subrv15_v(0)*msq15_v(j,k)
     .           + sub15(gg)*msqrv15(j,k,0)
     .           + sub15_v*msqrv15_v(j,k,0))
            msqv(2,j,k)= epinv*msqv(2,j,k) +
     .           ca*(subrv25(gg,0)*msq25(j,k)
     .           + subrv25_v(0)*msq25_v(j,k)
     .           + sub25(gg)*msqrv25(j,k,0)
     .           + sub25_v*msqrv25_v(j,k,0))
            msqv(3,j,k)= epinv*msqv(3,j,k) +
     .           -ca*(subrv5(0)*msq5(j,k)
     .           + (sub5(0)+sub5(1)+sub5(2))*msqrv5(j,k,0))

            msqx(3,j,k,:)=-ca*sub5(:)*msq5(j,k)
            
         endif

      endif
c---  END DEBUG


      enddo
      enddo

      return      
      end
