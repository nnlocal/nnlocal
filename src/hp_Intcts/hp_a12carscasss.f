      subroutine hp_a12carscasssggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      if (xa.gt.xb) call hp_AGTBa12carscasssggg(xa,xb,res)
      if (xb.gt.xa) call hp_BGTAa12carscasssggg(xa,xb,res)
      return
      end

c$$$      subroutine hp_a12carscasssgqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscasssgqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscasssgqg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscasssqgg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscasssqgg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscasssqgg(xa,xb,res)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine hp_a12carscasssqqg(xa,xb,res)
c$$$      implicit none
c$$$      include 'hp_types.h'
c$$$      real(ki) xa,xb,res(0:2,0:2,-4:0)
c$$$      if (xa.gt.xb) call hp_AGTBa12carscasssqqg(xa,xb,res)
c$$$      if (xb.gt.xa) call hp_BGTAa12carscasssqqg(xa,xb,res)
c$$$      return
c$$$      end
