      subroutine a12carscrsggg(xa,xb,res)
      implicit none
      include 'types.h'
      real(ki) xa,xb,res(2,2,-4:0)
      if (xa.gt.xb) call AGTBa12carscrsggg(xa,xb,res)
      if (xb.gt.xa) call BGTAa12carscrsggg(xa,xb,res)
      return
      end

c$$$      subroutine a12carscrsgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscrsgqq(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsgqq(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsgqq(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscrsqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscrsqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsqqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscrsqqq(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsqqq(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsqqq(xa,xb,res)
c$$$      return
c$$$      end
      
