      subroutine phase3(r1,r2,r3,r6,r7,p1,p2,p3,p4,p5,p6,p7,wt)
c----generate phase space for 2-->3 process
c----r(mxdim),p1(4),p2(4) are inputs 
c----incoming p1 and p2 reversed in sign from physical values 
c----i.e. phase space for -p1-p2 --> p3+p4+p5
c----with all 2 pi's (ie 1/(2*pi)^5)
c----(p4,p5) are dummies
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'masses.f'
      include 'process.f'
      include 'hdecaymode.f'

      integer j
      double precision r1,r2,r3,r6,r7
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision p12(4),p34(4),smin
      double precision wt,wt125,wt34,wt0,m5
      parameter(wt0=1d0/twopi)

      double precision s12,m34,s34,w3
      
      m5=0d0

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      p6(j)=0d0
      p7(j)=0d0
      enddo

      call breitw(r1,wsqmin,wsqmax,hmass,hwidth,s34,w3)
      m34=dsqrt(s34)
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      if (s12.lt.s34) goto 99
      call phi3m(r2,r3,p12,p5,p34,m5,m34,wt125,*99)
      call phi3m0(r6,r7,p34,p3,p4,wt34,*99)
      wt=wt0*wt125*wt34*w3

      return

c--- decay 34-system
      if (hdecaymode == 'tlta') then
      call phi3m(r6,r7,p34,p3,p4,mtau,mtau,wt34,*99)
      elseif (hdecaymode == 'bqba') then
      call phi3m(r6,r7,p34,p3,p4,mb,mb,wt34,*99)
      else
      call phi3m0(r6,r7,p34,p3,p4,wt34,*99)
      endif
      wt=wt0*wt125*wt34
      return
      
 99   continue
      wt=0d0
      return
      end

