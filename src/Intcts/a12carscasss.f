      subroutine a12carscasssggg(xa,xb,res)
      implicit none
      include 'types.h'
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      if (xa.gt.xb) call AGTBa12carscasssggg(xa,xb,res)
      if (xb.gt.xa) call BGTAa12carscasssggg(xa,xb,res)
      return
      end

c$$$      subroutine a12carscasssgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscasssgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscasssgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscasssqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscasssqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscasssqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine a12carscasssqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call AGTBa12carscasssqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call BGTAa12carscasssqqg(xa,xb,res)
c$$$      return
c$$$      end
