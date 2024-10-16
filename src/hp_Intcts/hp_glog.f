      function hp_glog(a, b, z)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) hp_glog
      complex(ki) a,b,z
      real(ki) absa,absb,absz,tol
      tol=1e-14_ki
      absa=abs(a)
      absb=abs(b)
      absz=abs(z)
      if (absa.lt.tol.and.absb.lt.tol) then
         call hp_glog00z(z,hp_glog)
         return
      else if (absz.lt.tol.and.absb.lt.tol) then
         call hp_glogz00(a,hp_glog)
         return
      else if (absz.lt.tol) then
         call hp_glogaz0(a,z,hp_glog)
         return
      else if (absa.lt.tol) then
         call hp_glog0az(b,z,hp_glog)
         return
      else if (absb.lt.tol) then
         call hp_gloga0z(a,z,hp_glog)
         return
      else if (abs((b-a)/b).lt.tol) then
         call hp_glogaaz(a,z,hp_glog)
         return
      else if (abs((z-b)/z).lt.tol) then
         call hp_glogazz(a,z,hp_glog)
         return
      else
         call hp_glogabz(a,b,z,hp_glog)
      endif
      return
      end

ccccccccccccccccccccccc

      subroutine hp_glog00z(z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) z, out
      out = half*log(z)**2
      return
      end

ccccccccccccccccccccccc

      subroutine hp_glogz00(z, out)
      implicit none
      include 'hp_types.h'
      complex(ki) z, out
      out = 0._ki
      return
      end

ccccccccccccccccccccccccc      

      subroutine hp_glogaz0(a, z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) a, z, out
      out = 0._ki
      return
      end

ccccccccccccccccccccccccc

      subroutine hp_glog0az(a, z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) a, z, out, hp_cli2
      out = -hp_cli2(z/a)
      return
      end

cccccccccccccccccccccccccc
      
      subroutine hp_gloga0z(a, z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) a, z, out, hp_cli2
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
            out = log(z)*log(1-z/a) + hp_cli2(z/a)
         else if (r>1) then
            out = log(a)*log(1-z/a) + pisqo6 - hp_cli2(1-z/a)
         endif
      else
         out = log(z)*log(1-z/a) + hp_cli2(z/a)
      endif
      return
      end

cccccccccccccccccccccccccc

      subroutine hp_glogaaz(a, z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) a, z, out, hp_cli2
      
      out = half*log(1-z/a)**2
      return
      end

cccccccccccccccccccccccccc

      subroutine hp_glogazz(a, z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) a, z, out, hp_cli2, aonz
      real(ki) r, aonzr
      aonz = a/z
      aonzr = real(aonz)
      r = abs(aonzr)
      if(r>1) then
         out=-hp_cli2(z/(-a+z))
         return
      else
         out=-pisqo6+log((-a+z)/z)*log(1-z/a)+hp_cli2(a/(a-z))
         return
      endif
      end

ccccccccccccccccccccccccccc
      
      subroutine hp_glogabz(a, b, z, out)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) a, b, z, out, hp_cli2
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
            out = log((b-a)/b)*log(1-z/a) + hp_cli2(a/(a-b))-hp_cli2((a-z)/(a-b))
         else if((absa>absb .and. abs(aonzi)<tiny .and. abs(bonzi)<tiny
     .           ) .or.
     .           (absai>absbi)) then
            out = log((-a+z)/(-a+b))*log(1-z/b) - hp_cli2(b/(-a+b))
     .           + hp_cli2((b-z)/(-a+b))
         else if (zr<br .and. br<ar .and. ar<0) then
            out = log((b-a)/b)*log(1-z/a) + hp_cli2(a/(a-b)) -
     .           (-hp_cli2(1-(a-z)/(a-b))+pisqo6-log((a-z)/(a-b))
     .           *log(1-(a-z)/(a-b)))     
         else
            out = log((b-a)/b)*log(1-z/a) + hp_cli2(a/(a-b))-hp_cli2((a-z)/(a-b))
         endif
      else
         if (absai>absbi) then
            out = hp_cli2((b-z)/(b-a)) - hp_cli2(b/(b-a))
     .            + log(1-z/b)*log((z-a)/(b-a))
         else
            out = hp_cli2(a/(a-b)) - hp_cli2((a-z)/(a-b))
     .           + log(1-z/a)*log((b-a)/b)
         endif
      endif
      return
      end
