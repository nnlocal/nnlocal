      subroutine gg_hgagag_gs_colorful(p,msq)
************************************************************************
*     Author: G.Somogyi, F.Tramontano                                  *
*     Nov, 2020.                                                       *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     g(-p1)+g(-p2) -->  H+ parton(p5) + parton(p6)
c                        |
c                        -->gamma(p3)+gamma(p4)

      implicit none 
      include 'constants.f'
      include 'order.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd,iglue1,iglue2
c --- remember: nd will count the dipoles
      
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision 
     &     msq15(-nf:nf,-nf:nf),msq25(-nf:nf,-nf:nf),
     &     msq16(-nf:nf,-nf:nf),msq26(-nf:nf,-nf:nf),
     &     msq56(-nf:nf,-nf:nf),
     &     msq15_v(-nf:nf,-nf:nf),msq25_v(-nf:nf,-nf:nf),
     &     msq16_v(-nf:nf,-nf:nf),msq26_v(-nf:nf,-nf:nf),
     &     msq56_v(-nf:nf,-nf:nf),
     &     msq5(-nf:nf,-nf:nf),msq6(-nf:nf,-nf:nf),
     &     sub15(4),sub25(4),sub16(4),sub26(4),sub56(4),
     &     sub15_v,sub25_v,sub16_v,sub26_v,sub56_v,
     &     sub5,sub6
      double precision
     &     msq156(-nf:nf,-nf:nf),msq156_v(6,-nf:nf,-nf:nf),
     &     msq256(-nf:nf,-nf:nf),msq256_v(6,-nf:nf,-nf:nf),
     &     msq1526(-nf:nf,-nf:nf),msq1526_v(2,-nf:nf,-nf:nf),
     &     msq1526_vv(-nf:nf,-nf:nf),
     &     msq1625(-nf:nf,-nf:nf),msq1625_v(2,-nf:nf,-nf:nf),
     &     msq1625_vv(-nf:nf,-nf:nf),
     &     msq_s56(-nf:nf,-nf:nf),
     &     msq_cs156(-nf:nf,-nf:nf),msq_cs156_v(-nf:nf,-nf:nf),
     &     msq_cs165(-nf:nf,-nf:nf),msq_cs165_v(-nf:nf,-nf:nf),
     &     msq_cs256(-nf:nf,-nf:nf),msq_cs256_v(-nf:nf,-nf:nf),
     &     msq_cs265(-nf:nf,-nf:nf),msq_cs265_v(-nf:nf,-nf:nf),
     &     sub156(4),sub156_v(6),
     &     sub256(4),sub256_v(6),
     &     sub1526(4),sub1526_v(2),sub1526_vv(4),
     &     sub1625(4),sub1625_v(2),sub1625_vv(4),
     &     sub_cs156(4),sub_cs156_v,
     &     sub_cs165(4),sub_cs165_v,
     &     sub_cs256(4),sub_cs256_v,
     &     sub_cs265(4),sub_cs265_v,
     &     sub_s56(4)
      double precision
     &     msq_c56c156(-nf:nf,-nf:nf),msq_c56c156_v(2,-nf:nf,-nf:nf),
     &     msq_c56c256(-nf:nf,-nf:nf),msq_c56c256_v(2,-nf:nf,-nf:nf),
     &     msq_c15c156(-nf:nf,-nf:nf),msq_c15c156_v(2,-nf:nf,-nf:nf),
     &     msq_c16c165(-nf:nf,-nf:nf),msq_c16c165_v(2,-nf:nf,-nf:nf),
     &     msq_c25c256(-nf:nf,-nf:nf),msq_c25c256_v(2,-nf:nf,-nf:nf),
     &     msq_c26c265(-nf:nf,-nf:nf),msq_c26c265_v(2,-nf:nf,-nf:nf),
     &     msq_c56s56(-nf:nf,-nf:nf),
     &     msq_c15cs156(-nf:nf,-nf:nf),msq_c15cs156_v(-nf:nf,-nf:nf),
     &     msq_c16cs165(-nf:nf,-nf:nf),msq_c16cs165_v(-nf:nf,-nf:nf),
     &     msq_c25cs256(-nf:nf,-nf:nf),msq_c25cs256_v(-nf:nf,-nf:nf),
     &     msq_c26cs265(-nf:nf,-nf:nf),msq_c26cs265_v(-nf:nf,-nf:nf),
     &     msq_c15c1526(-nf:nf,-nf:nf),msq_c15c1526_v(2,-nf:nf,-nf:nf),
     &     msq_c15c1526_vv(-nf:nf,-nf:nf),
     &     msq_c16c1625(-nf:nf,-nf:nf),msq_c16c1625_v(2,-nf:nf,-nf:nf),
     &     msq_c16c1625_vv(-nf:nf,-nf:nf),
     &     msq_c25c2516(-nf:nf,-nf:nf),msq_c25c2516_v(2,-nf:nf,-nf:nf),
     &     msq_c25c2516_vv(-nf:nf,-nf:nf),
     &     msq_c26c2615(-nf:nf,-nf:nf),msq_c26c2615_v(2,-nf:nf,-nf:nf),
     &     msq_c26c2615_vv(-nf:nf,-nf:nf),
     &     msq_s6c156(-nf:nf,-nf:nf),msq_s6c156_v(-nf:nf,-nf:nf),
     &     msq_s5c165(-nf:nf,-nf:nf),msq_s5c165_v(-nf:nf,-nf:nf),
     &     msq_s6c256(-nf:nf,-nf:nf),msq_s6c256_v(-nf:nf,-nf:nf),
     &     msq_s5c265(-nf:nf,-nf:nf),msq_s5c265_v(-nf:nf,-nf:nf),
     &     msq_s6cs156(-nf:nf,-nf:nf),msq_s6cs156_v(-nf:nf,-nf:nf),
     &     msq_s5cs165(-nf:nf,-nf:nf),msq_s5cs165_v(-nf:nf,-nf:nf),
     &     msq_s6cs256(-nf:nf,-nf:nf),msq_s6cs256_v(-nf:nf,-nf:nf),
     &     msq_s5cs265(-nf:nf,-nf:nf),msq_s5cs265_v(-nf:nf,-nf:nf),
     &     msq_s6s56(-nf:nf,-nf:nf),msq_s5s56(-nf:nf,-nf:nf),
     &     sub_c56c156(4),sub_c56c156_v(2),
     &     sub_c56c256(4),sub_c56c256_v(2),
     &     sub_c15c156(4),sub_c15c156_v(2),
     &     sub_c16c165(4),sub_c16c165_v(2),
     &     sub_c25c256(4),sub_c25c256_v(2),
     &     sub_c26c265(4),sub_c26c265_v(2),
     &     sub_c56s56(4),
     &     sub_c15cs156(4),sub_c15cs156_v,
     &     sub_c16cs165(4),sub_c16cs165_v,
     &     sub_c25cs256(4),sub_c25cs256_v,
     &     sub_c26cs265(4),sub_c26cs265_v,
     &     sub_c15c1526(4),sub_c15c1526_v(2),sub_c15c1526_vv(4), 
     &     sub_c16c1625(4),sub_c16c1625_v(2),sub_c16c1625_vv(4), 
     &     sub_c25c2516(4),sub_c25c2516_v(2),sub_c25c2516_vv(4),
     &     sub_c26c2615(4),sub_c26c2615_v(2),sub_c26c2615_vv(4), 
     &     sub_s6c156(4),sub_s6c156_v,
     &     sub_s5c165(4),sub_s5c165_v,
     &     sub_s6c256(4),sub_s6c256_v,
     &     sub_s5c265(4),sub_s5c265_v,
     &     sub_s6cs156(4),sub_s6cs156_v,
     &     sub_s5cs165(4),sub_s5cs165_v,
     &     sub_s6cs256(4),sub_s6cs256_v,
     &     sub_s5cs265(4),sub_s5cs265_v,
     &     sub_s6s56(4),sub_s5s56(4)


      parameter(iglue1=5,iglue2=6)
      external gg_hgagag,gg_hgagag_gvec
      external gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec
      ndmax=7

c--- A1 terms
c--- calculate single collinear subtractions      
      call dipsC2_colorful(1,p,1,iglue1,sub15,sub15_v,msq15,msq15_v,
     .     gg_hgagag,gg_hgagag_gvec)
      call dipsC2_colorful(2,p,1,iglue2,sub16,sub16_v,msq16,msq16_v,
     .     gg_hgagag,gg_hgagag_gvec)
      call dipsC2_colorful(3,p,2,iglue1,sub25,sub25_v,msq25,msq25_v,
     .     gg_hgagag,gg_hgagag_gvec)
      call dipsC2_colorful(4,p,2,iglue2,sub26,sub26_v,msq26,msq26_v,
     .     gg_hgagag,gg_hgagag_gvec)
      call dipsC2_colorful(5,p,iglue1,iglue2,sub56,sub56_v,
     .     msq56,msq56_v,gg_hgagag,gg_hgagag_gvec)
      
c---  calculate single soft subtractions      
      call dipsS1_colorful(6,p,iglue1,sub5,msq5,gg_hgagag)
      call dipsS1_colorful(7,p,iglue2,sub6,msq6,gg_hgagag)

      if (order.eq.1) goto 10

      
c--- A2 terms
c--- calculate triple collinear subtractions      
      call dipsC3_colorful(8,p,1,iglue1,iglue2,sub156,sub156_v,
     .     msq156,msq156_v,gg_hgamgam,gg_hgamgam_gvec)      
      call dipsC3_colorful(9,p,2,iglue1,iglue2,sub256,sub256_v,
     .     msq256,msq256_v,gg_hgamgam,gg_hgamgam_gvec)
c--- calculate double collinear subtractions      
      call dipsC22_colorful(10,p,1,iglue1,2,iglue2,sub1526,sub1526_v,
     .     sub1526_vv,msq1526,msq1526_v,msq1526_vv,
     .     gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec)
      call dipsC22_colorful(11,p,1,iglue2,2,iglue1,sub1625,sub1625_v,
     .     sub1625_vv,msq1625,msq1625_v,msq1625_vv,
     .     gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec)
c--- calculate double soft-collinear subtractions      
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$      call dipsCS2_colorful(12,p,1,iglue1,iglue2,sub_cs156,sub_cs156_v,
c$$$     .     msq_cs156,msq_cs156_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$      call dipsCS2_colorful(13,p,1,iglue2,iglue1,sub_cs165,sub_cs165_v,
c$$$     .     msq_cs165,msq_cs165_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$      call dipsCS2_colorful(14,p,2,iglue1,iglue2,sub_cs256,sub_cs256_v,
c$$$     .     msq_cs256,msq_cs256_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$      call dipsCS2_colorful(15,p,2,iglue2,iglue1,sub_cs265,sub_cs265_v,
c$$$     .     msq_cs265,msq_cs265_v,gg_hgamgam,gg_hgamgam_gvec)      
c---  calculate double soft subtractions      
      call dipsS2_colorful(16,p,iglue1,iglue2,sub_s56,msq_s56,
     .     gg_hgamgam)
c---  A12 terms
c---  calculate single collinear of the triple collinear subtractions
      call dipsC2C3_colorful(17,p,5,6,1,sub_c56c156,sub_c56c156_v,
     .     msq_c56c156,msq_c56c156_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsC2C3_colorful(18,p,5,6,2,sub_c56c256,sub_c56c256_v,
     .     msq_c56c256,msq_c56c256_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsC2C3_colorful(19,p,1,5,6,sub_c15c156,sub_c15c156_v,
     .     msq_c15c156,msq_c15c156_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsC2C3_colorful(20,p,1,6,5,sub_c16c165,sub_c16c165_v,
     .     msq_c16c165,msq_c16c165_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsC2C3_colorful(21,p,2,5,6,sub_c25c256,sub_c25c256_v,
     .     msq_c25c256,msq_c25c256_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsC2C3_colorful(22,p,2,6,5,sub_c26c265,sub_c26c265_v,
     .     msq_c26c265,msq_c26c265_v,gg_hgamgam,gg_hgamgam_gvec)
c---  calculate single collinear of the double soft subtractions
      call dipsC2S2_colorful(23,p,5,6,sub_c56s56,msq_c56s56,gg_hgamgam)
c---  calculate single collinear of the double soft-collinear subtractions
c---  note: dipole 24 uses same mapping as dipole 12 above 
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$      call dipsC2CS2_colorful(24,p,1,5,6,sub_c15cs156,sub_c15cs156_v,
c$$$     .     msq_c15cs156,msq_c15cs156_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$c---  note: dipole 25 uses same mapping as dipole 13 above 
c$$$      call dipsC2CS2_colorful(25,p,1,6,5,sub_c16cs165,sub_c16cs165_v,
c$$$     .     msq_c16cs165,msq_c16cs165_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$c---  note: dipole 26 uses same mapping as dipole 14 above 
c$$$      call dipsC2CS2_colorful(26,p,2,5,6,sub_c25cs256,sub_c25cs256_v,
c$$$     .     msq_c25cs256,msq_c25cs256_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$c---  note: dipole 27 uses same mapping as dipole 15 above 
c$$$      call dipsC2CS2_colorful(27,p,2,6,5,sub_c26cs265,sub_c26cs265_v,
c$$$     .     msq_c26cs265,msq_c26cs265_v,gg_hgamgam,gg_hgamgam_gvec)
c---  calculate single collinear of the double collinear subtractions
      call dipsC2C22_colorful(28,p,1,iglue1,2,iglue2,
     .     sub_c15c1526,sub_c15c1526_v,sub_c15c1526_vv,
     .     msq_c15c1526,msq_c15c1526_v,msq_c15c1526_vv,
     .     gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec)
      call dipsC2C22_colorful(29,p,1,iglue2,2,iglue1,
     .     sub_c16c1625,sub_c16c1625_v,sub_c16c1625_vv,
     .     msq_c16c1625,msq_c16c1625_v,msq_c16c1625_vv,
     .     gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec)
      call dipsC2C22_colorful(30,p,2,iglue1,1,iglue2,
     .     sub_c25c2516,sub_c25c2516_v,sub_c25c2516_vv,
     .     msq_c25c2516,msq_c25c2516_v,msq_c25c2516_vv,
     .     gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec)
      call dipsC2C22_colorful(31,p,2,iglue2,1,iglue1,
     .     sub_c26c2615,sub_c26c2615_v,sub_c26c2615_vv,
     .     msq_c26c2615,msq_c26c2615_v,msq_c26c2615_vv,
     .     gg_hgamgam,gg_hgamgam_gvec,gg_hgamgam_gvecgvec)
c---  calculate single soft of the triple collinear subtractions
      call dipsS1C3_colorful(32,p,1,5,6,sub_s6c156,sub_s6c156_v,
     .     msq_s6c156,msq_s6c156_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsS1C3_colorful(33,p,1,6,5,sub_s5c165,sub_s5c165_v,
     .     msq_s5c165,msq_s5c165_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsS1C3_colorful(34,p,2,5,6,sub_s6c256,sub_s6c256_v,
     .     msq_s6c256,msq_s6c256_v,gg_hgamgam,gg_hgamgam_gvec)
      call dipsS1C3_colorful(35,p,2,6,5,sub_s5c265,sub_s5c265_v,
     .     msq_s5c265,msq_s5c265_v,gg_hgamgam,gg_hgamgam_gvec)
c---  calculate single soft of the doube soft-collinear subtractions
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$c---  note: dipole 36 uses same mapping as dipole 32 above 
c$$$      call dipsS1CS2_colorful(36,p,1,5,6,sub_s6cs156,sub_s6cs156_v,
c$$$     .     msq_s6cs156,msq_s6cs156_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$c---  note: dipole 37 uses same mapping as dipole 33 above 
c$$$      call dipsS1CS2_colorful(37,p,1,6,5,sub_s5cs165,sub_s5cs165_v,
c$$$     .     msq_s5cs165,msq_s5cs165_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$c---  note: dipole 38 uses same mapping as dipole 34 above 
c$$$      call dipsS1CS2_colorful(38,p,2,5,6,sub_s6cs256,sub_s6cs256_v,
c$$$     .     msq_s6cs256,msq_s6cs256_v,gg_hgamgam,gg_hgamgam_gvec)
c$$$c---  note: dipole 39 uses same mapping as dipole 35 above 
c$$$      call dipsS1CS2_colorful(39,p,2,6,5,sub_s5cs265,sub_s5cs265_v,
c$$$     .     msq_s5cs265,msq_s5cs265_v,gg_hgamgam,gg_hgamgam_gvec)
c---  calculate single soft of double soft subtractions      
      call dipsS1S2_colorful(40,p,iglue1,iglue2,sub_s6s56,msq_s6s56,
     .     gg_hgamgam)
      call dipsS1S2_colorful(41,p,iglue2,iglue1,sub_s5s56,msq_s5s56,
     .     gg_hgamgam)


 10   continue
      
      do j=-nf,nf
      do k=-nf,nf      
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo

c---  fully gluonic subprocess only for now      
      do j=-nf,nf
      do k=-nf,nf
         
         if ((j .eq. 0).and.(k .eq. 0)) then
c---  A1 single unresolved
            msq(1,j,k) = msq(1,j,k)
     .           +half*ca*(msq15(j,k)*sub15(gg)+msq15_v(j,k)*sub15_v)
            msq(2,j,k) = msq(2,j,k)
     .           +half*ca*(msq16(j,k)*sub16(gg)+msq16_v(j,k)*sub16_v)
            msq(3,j,k) = msq(3,j,k)
     .           +half*ca*(msq25(j,k)*sub25(gg)+msq25_v(j,k)*sub25_v)
            msq(4,j,k) = msq(4,j,k)
     .           +half*ca*(msq26(j,k)*sub26(gg)+msq26_v(j,k)*sub26_v)
            msq(5,j,k) = msq(5,j,k)
     .           +half*ca*(msq56(j,k)*sub56(gg)+msq56_v(j,k)*sub56_v)
            msq(6,j,k) = msq(6,j,k) + half*(-ca/2d0)*(msq5(j,k)*sub5)
            msq(7,j,k) = msq(7,j,k) + half*(-ca/2d0)*(msq6(j,k)*sub6)

            endif
         enddo
      enddo

      
      ndmax1=ndmax
      if (order.eq.1) return
      ndmax=41

      do j=-nf,nf
      do k=-nf,nf      
      do nd=ndmax1+1,ndmax
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo

      
      
      do j=-nf,nf
      do k=-nf,nf
         
         if ((j .eq. 0).and.(k .eq. 0)) then
            
c---  A2 double unresolved
c---  C_ars
            msq(8,j,k) = msq(8,j,k) + half*ca**2*(
     .           msq156(j,k)*sub156(gg)
     .           + msq156_v(1,j,k)*sub156_v(1)
     .           + msq156_v(2,j,k)*sub156_v(2)
     .           + msq156_v(3,j,k)*sub156_v(3)
     .           + msq156_v(4,j,k)*sub156_v(4)
     .           + msq156_v(5,j,k)*sub156_v(5)
     .           + msq156_v(6,j,k)*sub156_v(6))
            msq(9,j,k) = msq(9,j,k) + half*ca**2*(
     .           msq256(j,k)*sub256(gg)
     .           + msq256_v(1,j,k)*sub256_v(1)
     .           + msq256_v(2,j,k)*sub256_v(2)
     .           + msq256_v(3,j,k)*sub256_v(3)
     .           + msq256_v(4,j,k)*sub256_v(4)
     .           + msq256_v(5,j,k)*sub256_v(5)
     .           + msq256_v(6,j,k)*sub256_v(6))
c---  C_ar,bs
            msq(10,j,k) = msq(10,j,k) + half*ca**2*(
     .           msq1526(j,k)*sub1526(gg)
     .           + msq1526_v(1,j,k)*sub1526_v(1)
     .           + msq1526_v(2,j,k)*sub1526_v(2)
     .           + msq1526_vv(j,k)*sub1526_vv(gg))
            msq(11,j,k) = msq(11,j,k) + half*ca**2*(
     .           msq1625(j,k)*sub1625(gg)
     .           + msq1625_v(1,j,k)*sub1625_v(1)
     .           + msq1625_v(2,j,k)*sub1625_v(2)
     .           + msq1625_vv(j,k)*sub1625_vv(gg))
c---  CS_ar,s
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$            msq(12,j,k) = msq(12,j,k) + half*(-ca**2)*(
c$$$     .           msq_cs156(j,k)*sub_cs156(gg)
c$$$     .           + msq_cs156_v(j,k)*sub_cs156_v)
c$$$            msq(13,j,k) = msq(13,j,k) + half*(-ca**2)*(
c$$$     .           msq_cs165(j,k)*sub_cs165(gg)
c$$$     .           + msq_cs165_v(j,k)*sub_cs165_v)
c$$$            msq(14,j,k) = msq(14,j,k) + half*(-ca**2)*(
c$$$     .           msq_cs256(j,k)*sub_cs256(gg)
c$$$     .           + msq_cs256_v(j,k)*sub_cs256_v)
c$$$            msq(15,j,k) = msq(15,j,k) + half*(-ca**2)*(
c$$$     .           msq_cs265(j,k)*sub_cs265(gg)
c$$$     .           + msq_cs265_v(j,k)*sub_cs265_v)
            msq(12,j,k) = 0d0
            msq(13,j,k) = 0d0
            msq(14,j,k) = 0d0
            msq(15,j,k) = 0d0
c---  S_rs
            msq(16,j,k) = msq(16,j,k)
     .           + half*(ca**2)*(msq_s56(j,k)*sub_s56(gg))


c---  A12 iterated single unresolved
c---  C_rs C_ars
            msq(17,j,k) = msq(17,j,k) + half*ca**2*(
     .           msq_c56c156(j,k)*sub_c56c156(gg)
     .           + msq_c56c156_v(1,j,k)*sub_c56c156_v(1)
     .           + msq_c56c156_v(2,j,k)*sub_c56c156_v(2))
            msq(18,j,k) = msq(18,j,k) + half*ca**2*(
     .           msq_c56c256(j,k)*sub_c56c256(gg)
     .           + msq_c56c256_v(1,j,k)*sub_c56c256_v(1)
     .           + msq_c56c256_v(2,j,k)*sub_c56c256_v(2))
c---  C_ar C_ars
            msq(19,j,k) = msq(19,j,k) + half*ca**2*(
     .           msq_c15c156(j,k)*sub_c15c156(gg)
     .           + msq_c15c156_v(1,j,k)*sub_c15c156_v(1)
     .           + msq_c15c156_v(2,j,k)*sub_c15c156_v(2))            
            msq(20,j,k) = msq(20,j,k) + half*ca**2*(
     .           msq_c16c165(j,k)*sub_c16c165(gg)
     .           + msq_c16c165_v(1,j,k)*sub_c16c165_v(1)
     .           + msq_c16c165_v(2,j,k)*sub_c16c165_v(2))
            msq(21,j,k) = msq(21,j,k) + half*ca**2*(
     .           msq_c25c256(j,k)*sub_c25c256(gg)
     .           + msq_c25c256_v(1,j,k)*sub_c25c256_v(1)
     .           + msq_c25c256_v(2,j,k)*sub_c25c256_v(2))
            msq(22,j,k) = msq(22,j,k) + half*ca**2*(
     .           msq_c26c265(j,k)*sub_c26c265(gg)
     .           + msq_c26c265_v(1,j,k)*sub_c26c265_v(1)
     .           + msq_c26c265_v(2,j,k)*sub_c26c265_v(2))
c---  C_rs S_rs
            msq(23,j,k) = msq(23,j,k) + half*ca**2*(
     .           msq_c56s56(j,k)*sub_c56s56(gg))
c---  C_ar CS_ar,s
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$            msq(24,j,k) = msq(24,j,k) + half*(-ca**2)*(
c$$$     .           msq_c15cs156(j,k)*sub_c15cs156(gg)
c$$$     .           + msq_c15cs156_v(j,k)*sub_c15cs156_v)
c$$$            msq(25,j,k) = msq(25,j,k) + half*(-ca**2)*(
c$$$     .           msq_c16cs165(j,k)*sub_c16cs165(gg)
c$$$     .           + msq_c16cs165_v(j,k)*sub_c16cs165_v)
c$$$            msq(26,j,k) = msq(26,j,k) + half*(-ca**2)*(
c$$$     .           msq_c25cs256(j,k)*sub_c25cs256(gg)
c$$$     .           + msq_c25cs256_v(j,k)*sub_c25cs256_v)
c$$$            msq(27,j,k) = msq(27,j,k) + half*(-ca**2)*(
c$$$     .           msq_c26cs265(j,k)*sub_c26cs265(gg)
c$$$     .           + msq_c26cs265_v(j,k)*sub_c26cs265_v)
            msq(24,j,k) = 0d0
            msq(25,j,k) = 0d0
            msq(26,j,k) = 0d0
            msq(27,j,k) = 0d0
c---  C_ar C_ar,bs
            msq(28,j,k) = msq(28,j,k) + half*ca**2*(
     .           msq_c15c1526(j,k)*sub_c15c1526(gg)
     .           + msq_c15c1526_v(1,j,k)*sub_c15c1526_v(1)
     .           + msq_c15c1526_v(2,j,k)*sub_c15c1526_v(2)
     .           + msq_c15c1526_vv(j,k)*sub_c15c1526_vv(gg))
            msq(29,j,k) = msq(29,j,k) + half*ca**2*(
     .           msq_c16c1625(j,k)*sub_c16c1625(gg)
     .           + msq_c16c1625_v(1,j,k)*sub_c16c1625_v(1)
     .           + msq_c16c1625_v(2,j,k)*sub_c16c1625_v(2)
     .           + msq_c16c1625_vv(j,k)*sub_c16c1625_vv(gg))
            msq(30,j,k) = msq(30,j,k) + half*ca**2*(
     .           msq_c25c2516(j,k)*sub_c25c2516(gg)
     .           + msq_c25c2516_v(1,j,k)*sub_c25c2516_v(1)
     .           + msq_c25c2516_v(2,j,k)*sub_c25c2516_v(2)
     .           + msq_c25c2516_vv(j,k)*sub_c25c2516_vv(gg))
            msq(31,j,k) = msq(31,j,k) + half*ca**2*(
     .           msq_c26c2615(j,k)*sub_c26c2615(gg)
     .           + msq_c26c2615_v(1,j,k)*sub_c26c2615_v(1)
     .           + msq_c26c2615_v(2,j,k)*sub_c26c2615_v(2)
     .           + msq_c26c2615_vv(j,k)*sub_c26c2615_vv(gg))
c---  S_s C_ars
            msq(32,j,k) = msq(32,j,k) + half*(ca**2)*(
     .           msq_s6c156(j,k)*sub_s6c156(gg)
     .           + msq_s6c156_v(j,k)*sub_s6c156_v)            
            msq(33,j,k) = msq(33,j,k) + half*(ca**2)*(
     .           msq_s5c165(j,k)*sub_s5c165(gg)
     .           + msq_s5c165_v(j,k)*sub_s5c165_v)
            msq(34,j,k) = msq(34,j,k) + half*(ca**2)*(
     .           msq_s6c256(j,k)*sub_s6c256(gg)
     .           + msq_s6c256_v(j,k)*sub_s6c256_v)            
            msq(35,j,k) = msq(35,j,k) + half*(ca**2)*(
     .           msq_s5c265(j,k)*sub_s5c265(gg)
     .           + msq_s5c265_v(j,k)*sub_s5c265_v)

c---  S_s CS_ars
!---  NEW CANCELLATIONS IMPLEMENTED
c$$$            msq(36,j,k) = msq(36,j,k) + half*(-ca**2)*(
c$$$     .           msq_s6cs156(j,k)*sub_s6cs156(gg)
c$$$     .           + msq_s6cs156_v(j,k)*sub_s6cs156_v)            
c$$$            msq(37,j,k) = msq(37,j,k) + half*(-ca**2)*(
c$$$     .           msq_s5cs165(j,k)*sub_s5cs165(gg)
c$$$     .           + msq_s5cs165_v(j,k)*sub_s5cs165_v)
c$$$            msq(38,j,k) = msq(38,j,k) + half*(-ca**2)*(
c$$$     .           msq_s6cs256(j,k)*sub_s6cs256(gg)
c$$$     .           + msq_s6cs256_v(j,k)*sub_s6cs256_v)            
c$$$            msq(39,j,k) = msq(39,j,k) + half*(-ca**2)*(
c$$$     .           msq_s5cs265(j,k)*sub_s5cs265(gg)
c$$$     .           + msq_s5cs265_v(j,k)*sub_s5cs265_v)
            msq(36,j,k) = 0d0
            msq(37,j,k) = 0d0
            msq(38,j,k) = 0d0
            msq(39,j,k) = 0d0
c---  S_s S_rs
            msq(40,j,k) = msq(40,j,k)
     .           + half*(ca**2)*(msq_s6s56(j,k)*sub_s6s56(gg))
            msq(41,j,k) = msq(41,j,k)
     .           + half*(ca**2)*(msq_s5s56(j,k)*sub_s5s56(gg))            
         endif
      
      enddo
      enddo

      return
      end

