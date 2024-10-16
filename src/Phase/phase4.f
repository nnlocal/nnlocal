      subroutine phase4(r,p1,p2,p3,p4,p5,p6,wt,*)
      implicit none
      include 'constants.f'
      include 'heavyflav.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'process.f'
      include 'breit.f'
      include 'zerowidth.f'
      include 'limits.f'
c---- generate phase space for 2-->4 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p12(4),p34(4),p56(4),p345(4)
      double precision wt,wt3456,wt34,wt56,wt0,jac
      double precision p35(4),p46(4),wt35,wt46,s3min,mres,mres2
      integer j
      parameter(wt0=1d0/twopi**2)
      
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

      call breitw(r(6),wsqmin,wsqmax,hmass,hwidth,mres2,jac)
      mres=dsqrt(mres2)
      
      call phi1_2m_nobw(0d0,r(1),r(2),r(3),mres2,p12,p6,p345,wt3456,*99)
      call phi3m(r(4),r(5),p345,p5,p34,0d0,mres,wt56,*99)
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)
      wt=wt0*wt3456*wt34*wt56*jac
      
      if (debug) write(6,*) 'wt in phase4',wt
      return

 99   wt=0d0
      return 1
      end

