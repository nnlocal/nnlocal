      function hp_Li3(x)
      implicit none
      include 'hp_types.h'
C     returns Li_3(x) for real x, minf < x < 1
      complex(ki) hp_Li3
      real(ki) x
      complex(ki) hp_wgplg

      if (x .gt. 1d0) then
        write(6,*) 'x>1 in hp_Li3 function'
        stop
      endif
      
      hp_Li3 = hp_wgplg(2,1,x)
c      hp_Li3 = real(hp_wgplg(2,1,x),ki)

      end
