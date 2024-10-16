      subroutine a12carscrsssggg(xa,xb,res)
      implicit none
      include 'types.h'
      real(ki) xa,xb,res(2,2,-4:0)
      if (xa.gt.xb) call AGTBa12carscrsssggg(xa,xb,res)
      if (xb.gt.xa) call BGTAa12carscrsssggg(xa,xb,res)
      return
      end

c$$$      subroutine a12carscrsssgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsssgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsssgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscrsssqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsssqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsssqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscrsssqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscrsssqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscrsssqqg(xa,xb,res)
c$$$      return
c$$$      end
