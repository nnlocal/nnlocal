      subroutine a1a1crscasggg(xa,xb,res)
      implicit none
      include 'types.h'
      real(ki) xa,xb,res(2,2,-4:0)
      if (xa.gt.xb) call AGTBa1a1crscasggg(xa,xb,res)
      if (xb.gt.xa) call BGTAa1a1crscasggg(xa,xb,res)
      return
      end

c$$$      subroutine a1a1crscasgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa1a1crscasgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa1a1crscasgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a1a1crscasgqq(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa1a1crscasgqq(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa1a1crscasgqq(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a1a1crscasqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa1a1crscasqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa1a1crscasqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a1a1crscasqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa1a1crscasqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa1a1crscasqqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a1a1crscasqqq(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(2,2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa1a1crscasqqq(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa1a1crscasqqq(xa,xb,res)
c$$$      return
c$$$      end

      
      
