      subroutine iopmat1(p,xa,xb)
      implicit none
      include 'constants.f'
      include 'order.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'qcdcouple.f'
      include 'Ioperators.f'
      include 'scale.f'
      include 'facscale.f'
      include 'agq.f'
      include 'b0.f'
      integer i,j,k
      double precision xl12,p(mxpart,4),dot
      double precision xa,xb
      double precision resab(2,2,-2:2),resba(2,2,-2:2)

      character *60 FMT3
      integer eps
      
      FMT3="(A20,4(F30.17,','),F30.16,'};')"

      
      xl12=dlog(two*dot(p,1,2)/musq)
      call  a1ncar0gg_loprec(xa,xb,resab)
      call  a1ncar0gg_loprec(xb,xa,resba)
      ! (a,b,c) stand for incoming, hard, final
      do i=1,2
         do j=1,2
        I10op1(g,g,g,i,j)=ason2pi*(
     1    resab(i,j,-2)*(epinv*(epinv2-xl12)+xl12**2/2d0)
     2   +resab(i,j,-1)*(epinv-xl12)
     3   +resab(i,j, 0))
        I10op2(g,g,g,i,j)=ason2pi*(
     1    resba(i,j,-2)*(epinv*(epinv2-xl12)+xl12**2/2d0)
     2   +resba(i,j,-1)*(epinv-xl12)
     3   +resba(i,j, 0))
         enddo
      enddo

c      write(6,*) ason2pi
c      do i=1,2
c         do j=1,2
c            write(6,*) 'i=',i, 'j=',j
c            write(6,*) resab(i,j,-2),resab(i,j,-1),resab(i,j,0)
c            write(6,*) resba(i,j,-2),resba(i,j,-1),resba(i,j,0)
c            write(6,*)
c         enddo
c      enddo
      
      
      return
      end

      
      
