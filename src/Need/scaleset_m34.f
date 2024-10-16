      subroutine scaleset_m34(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3 and 4
      implicit none
      include 'constants.f'
      include 'process.f'
      double precision p(mxpart,4),mu0

      if (case .eq. 'H_gaga') then
        mu0=(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &     -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2       
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*) 'dynamicscale m(34) not supported for this process.'
        stop
      endif
      
      return
      end
      
