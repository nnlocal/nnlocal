      function mylog(z)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) mylog
      complex(ki) z
      real(ki) zr,zi,tol
      tol=1e-28_ki
      zr = real(z)
      zi = aimag(z)
      if ((abs(zi/zr).lt.tol).and.(zr.lt.zip)) then
         mylog=log(-zr)+im*pi
      else
         mylog=log(z)
      endif
      return
      end
      
            
