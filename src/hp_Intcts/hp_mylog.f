      function hp_mylog(z)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) hp_mylog
      complex(ki) z
      real(ki) zr,zi,tol
      tol=1e-14_ki
      zr = real(z)
      zi = aimag(z)
      if ((abs(zi/zr).lt.tol).and.(zr.lt.zip)) then
         hp_mylog=log(-zr)+im*pi
      else
         hp_mylog=log(z)
      endif
      return
      end
      
            
