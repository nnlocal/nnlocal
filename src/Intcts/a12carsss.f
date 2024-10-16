      subroutine a12carsssggg(xa,xb,res)
      implicit none
      include 'types.h'
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      if (xa.gt.xb) call AGTBa12carsssggg(xa,xb,res)
      if (xb.gt.xa) call BGTAa12carsssggg(xa,xb,res)
      return
      end

c$$$      subroutine a12carsssgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carsssgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carsssgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carsssqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carsssqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carsssqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carsssqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carsssqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carsssqqg(xa,xb,res)
c$$$      return
c$$$      end
