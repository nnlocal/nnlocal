      subroutine hp_a1a1crssscarggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      real(ki) xa,xb,res(2,2,-4:0)
      if (xa.gt.xb) call hp_AGTBa1a1crssscarggg(xa,xb,res)
      if (xb.gt.xa) call hp_BGTAa1a1crssscarggg(xa,xb,res)
      return
      end

c$$$      subroutine hp_a1a1crssscargqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa1a1crssscargqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa1a1crssscargqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a1a1crssscarqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa1a1crssscarqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa1a1crssscarqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a1a1crssscarqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa1a1crssscarqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa1a1crssscarqqg(xa,xb,res)
c$$$      return
c$$$      end
