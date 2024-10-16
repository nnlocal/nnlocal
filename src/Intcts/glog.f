      function glog(a, b, z)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) glog
      complex(ki) a,b,z
      real(ki) absa,absb,absz,tol
      tol=1e-28_ki
      absa=abs(a)
      absb=abs(b)
      absz=abs(z)
      if (absa.lt.tol.and.absb.lt.tol) then
         call glog00z(z,glog)
         return
      else if (absz.lt.tol.and.absb.lt.tol) then
         call glogz00(a,glog)
         return
      else if (absz.lt.tol) then
         call glogaz0(a,z,glog)
         return
      else if (absa.lt.tol) then
         call glog0az(b,z,glog)
         return
      else if (absb.lt.tol) then
         call gloga0z(a,z,glog)
         return
      else if (abs((b-a)/b).lt.tol) then
         call glogaaz(a,z,glog)
         return
      else if (abs((z-b)/z).lt.tol) then
         call glogazz(a,z,glog)
         return
      else
         call glogabz(a,b,z,glog)
      endif
      return
      end

ccccccccccccccccccccccc

      subroutine glog00z(z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) z, out, mylog
      out = half*mylog(z)**2
      return
      end

ccccccccccccccccccccccc

      subroutine glogz00(z, out)
      implicit none
      include 'types.h'
      complex(ki) z, out
      out = 0._ki
      return
      end

ccccccccccccccccccccccccc      

      subroutine glogaz0(a, z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) a, z, out
      out = 0._ki
      return
      end

ccccccccccccccccccccccccc

      subroutine glog0az(a, z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) a, z, out, cli2
      out = -cli2(z/a)
      return
      end

cccccccccccccccccccccccccc
      
      subroutine gloga0z(a, z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) a, z, out, cli2, mylog
      real(ki) ar, zr, ai, zi, r
      real(ki) tiny

      if (ki==dp) then
         tiny = 1.0e-14_ki
      else if (ki==qp) then
         tiny = 1.0e-28_ki
      endif
      
      ar = real(a)
      zr = real(z)
      ai = aimag(a)
      zi = aimag(z)
      r = zr/ar
      
      if(zi/zr<tiny .and. ai/ar<tiny) then
         if (r<1) then
            out = mylog(z)*mylog(1-z/a) + cli2(z/a)
         else if (r>1) then
            out = mylog(a)*mylog(1-z/a) + pisqo6 - cli2(1-z/a)
         endif
      else
         out = mylog(z)*mylog(1-z/a) + cli2(z/a)
      endif
      return
      end

cccccccccccccccccccccccccc

      subroutine glogaaz(a, z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) a, z, out, cli2, mylog
      
      out = half*mylog(1-z/a)**2
      return
      end

cccccccccccccccccccccccccc

      subroutine glogazz(a, z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) a, z, out, cli2, mylog, aonz
      real(ki) r, aonzr
      aonz = a/z
      aonzr = real(aonz)
      r = abs(aonzr)
      if(r>1) then
         out=-cli2(z/(-a+z))
         return
      else
         out=-pisqo6+mylog((-a+z)/z)*mylog(1-z/a)+cli2(a/(a-z))
         return
      endif
      end

ccccccccccccccccccccccccccc
      
      subroutine glogabz(a, b, z, out)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) a, b, z, out, cli2, mylog
      complex(ki) aonz, bonz
      real(ki) aonzr, aonzi, bonzr, bonzi, absa, absb, absai, absbi, mod
     .     a, modb, tiny, ar, br, zr, ai, bi, zi

      if (ki==dp) then
         tiny = 1.0e-14_ki
      else if (ki==qp) then
         tiny = 1.0e-28_ki
      endif

      ar = real(a)
      br = real(b)
      zr = real(z)
      ai = aimag(a)
      bi = aimag(b)
      zi = aimag(z)
      
      aonz = a/z
      bonz = b/z
      aonzr = real(aonz)
      aonzi = aimag(aonz)
      bonzr = real(bonz)
      bonzi = aimag(bonz)
      
      absa = abs(aonzr)
      absb = abs(bonzr)
      absai = abs(aonzi)
      absbi = abs(bonzi)

      if (abs(ai/ar)<tiny .and. abs(bi/br)<tiny .and. abs(zi/zr)<tiny) then
         if (0<br .and. br<ar .and. ar<zr) then
            out = mylog((b-a)/b)*mylog(1-z/a) + cli2(a/(a-b))-cli2((a-z)/(a-b))
         else if((absa>absb .and. abs(aonzi)<tiny .and. abs(bonzi)<tiny
     .           ) .or.
     .           (absai>absbi)) then
            out = mylog((-a+z)/(-a+b))*mylog(1-z/b) - cli2(b/(-a+b))
     .           + cli2((b-z)/(-a+b))
         else if (zr<br .and. br<ar .and. ar<0) then
            out = mylog((b-a)/b)*mylog(1-z/a) + cli2(a/(a-b)) -
     .           (-cli2(1-(a-z)/(a-b))+pisqo6-mylog((a-z)/(a-b))
     .           *mylog(1-(a-z)/(a-b)))     
         else
            out = mylog((b-a)/b)*mylog(1-z/a) + cli2(a/(a-b))-cli2((a-z)/(a-b))
         endif
      else
         if (absai>absbi) then
            out = cli2((b-z)/(b-a)) - cli2(b/(b-a))
     .            + mylog(1-z/b)*mylog((z-a)/(b-a))
         else
            out = cli2(a/(a-b)) - cli2((a-z)/(a-b))
     .           + mylog(1-z/a)*mylog((b-a)/b)
         endif
      endif
      return
      end
