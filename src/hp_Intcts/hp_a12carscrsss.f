      subroutine hp_a12carscrsssggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      real(ki) xa,xb,res(2,2,-4:0)
      if (xa.gt.xb) call hp_AGTBa12carscrsssggg(xa,xb,res)
      if (xb.gt.xa) call hp_BGTAa12carscrsssggg(xa,xb,res)
      return
      end

c$$$      subroutine hp_a12carscrsssgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsssgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsssgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscrsssqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsssqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsssqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscrsssqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsssqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsssqqg(xa,xb,res)
c$$$      return
c$$$      end
