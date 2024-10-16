      subroutine hp_a12carsssggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      if (xa.gt.xb) call hp_AGTBa12carsssggg(xa,xb,res)
      if (xb.gt.xa) call hp_BGTAa12carsssggg(xa,xb,res)
      return
      end

c$$$      subroutine hp_a12carsssgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carsssgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carsssgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carsssqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carsssqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carsssqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carsssqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carsssqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carsssqqg(xa,xb,res)
c$$$      return
c$$$      end
