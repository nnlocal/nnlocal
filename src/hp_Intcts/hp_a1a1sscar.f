      subroutine hp_a1a1sscarggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      if (xa.gt.xb) call hp_AGTBa1a1sscarggg(xa,xb,res)
      if (xb.gt.xa) call hp_BGTAa1a1sscarggg(xa,xb,res)
      return
      end

c$$$      subroutine hp_a1a1sscargqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa1a1sscargqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa1a1sscargqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a1a1sscarqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa1a1sscarqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa1a1sscarqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a1a1sscarqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa1a1sscarqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa1a1sscarqqg(xa,xb,res)
c$$$      return
c$$$      end
