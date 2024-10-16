      function Li3(x)
      implicit none
      include 'types.h'
C     returns Li_3(x) for real x, minf < x < 1
      complex(ki) Li3
      real(ki) x
      complex(ki) wgplg

      if (x .gt. 1d0) then
        write(6,*) 'x>1 in Li3 function, src/Lib/Li3.f'
        stop
      endif
      
c      Li3 = wgplg(2,1,x)
      Li3 = real(wgplg(2,1,x),ki)

      end
