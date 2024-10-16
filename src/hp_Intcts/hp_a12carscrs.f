      subroutine hp_a12carscrsggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      real(ki) xa,xb,res(2,2,-4:0)
      if (xa.gt.xb) call hp_AGTBa12carscrsggg(xa,xb,res)
      if (xb.gt.xa) call hp_BGTAa12carscrsggg(xa,xb,res)
      return
      end

c$$$      subroutine hp_a12carscrsgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscrsgqq(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsgqq(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsgqq(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscrsqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscrsqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsqqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscrsqqq(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscrsqqq(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscrsqqq(xa,xb,res)
c$$$      return
c$$$      end
      
