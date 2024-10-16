      subroutine gencheck(r,p,wt4,*)
************************************************************************
*     Generate PS points tending to any given limit                    *
*     Gabor Somogyi March 2014                                         *
*     Adapted by F. Tramontano                                         *
************************************************************************
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'seed.f'
      include 'npart.f'
      include 'singular.f'
      include 'chkseed.f'

      double precision r(mxdim),p(mxpart,4)
      double precision sqrts,wt0,wt4,wt
      parameter(wt0=1d0/twopi**2)
      common/energy/sqrts
      double precision q(4,mxpart),rn,dummy
      integer i,j,nu
      integer mynpart


      sqrts=200d0
      mynpart=npart
      npart=npart-1
      
      do j=1,mxpart     
      do nu=1,4     
      p(j,nu)=0d0
      enddo     
      enddo     

c      p(1,4)=-0.5d0*sqrts
c      p(1,1)=0d0
c      p(1,2)=0d0
c      p(1,3)=-0.5d0*sqrts
     
c      p(2,4)=-0.5d0*sqrts
c      p(2,1)=0d0
c      p(2,2)=0d0
c      p(2,3)=+0.5d0*sqrts

c--- 3 is problematic
      do i=1,cseed
      dummy=rn(1)
      enddo


      
C---Which limit to approach?
C---Types of limits implemented
c     coll2:      FF    ip || rp
c     icoll2:     IF    ip || rp
c     coll3:      FFF   ip || rp || sp
c     icoll3:     IFF   ip || rp || sp      
c     coll22:     FF,FF ip || rp and jp || sp
c     icoll22:    IF,FF ip || rp and jp || sp
c     iicoll22:   IF,IF ip || rp and jp || sp            
c     soft1:      F     rp -> 0
c     soft2:      FF    rp ->0 and sp -> 0
c     softcoll3:  FF,F ip || rp and sp -> 0
c     isoftcoll3: IF,F ip || rp and sp -> 0
      
c      limit='coll2'
c      ip = 4
c      rp = 5

c      limit='icoll2'
c      ip = 2
c      rp = 4
    
c      limit='coll3'
c      ip = 4
c      rp = 5
c      sp = 6

c      limit='icoll3'
c      ip = 1
c      rp = 4
c      sp = 5

c      limit='coll22'
c      ip = 3
c      rp = 4
c      jp = 5
c      sp = 6

c      limit='icoll22'
c      ip = 1
c      rp = 4
c      jp = 5
c      sp = 6

c      limit='iicoll22'
c      ip = 1
c      rp = 4
c      jp = 2
c      sp = 5            

c      limit='soft1'
c      rp = 4

c      limit='soft2'
c      rp = 4
c      sp = 5

c      limit='softcoll3'
c      ip = 4
c      rp = 5
c      sp = 6

c      limit='isoftcoll3'
c      ip = 1
c      rp = 4
c      sp = 5      



      call genlimit(limit,q,ip,rp,jp,sp)

      npart=mynpart

c      write(6,*) q(:,1),sqrt(abs(q(4,1)**2-q(1,1)**2
c     .                          -q(2,1)**2-q(3,1)**2))
c      write(6,*) q(:,2),sqrt(abs(q(4,2)**2-q(1,2)**2
c     .                          -q(2,2)**2-q(3,2)**2))
c      write(6,*) q(:,3),sqrt(abs(q(4,3)**2-q(1,3)**2
c     .                          -q(2,3)**2-q(3,3)**2))
c      write(6,*) q(:,4),sqrt(abs(q(4,4)**2-q(1,4)**2
c     .                          -q(2,4)**2-q(3,4)**2))
c      write(6,*) q(:,5),sqrt(abs(q(4,5)**2-q(1,5)**2
c     .                          -q(2,5)**2-q(3,5)**2))


      
      call phi3(rn(1),rn(1),q(:,3),p(3,:),p(4,:),wt)
      p(1,:)=q(:,1)
      p(2,:)=q(:,2)
      p(5,:)=q(:,4)
      p(6,:)=q(:,5)

      wt4=wt0*wt

      return
 999  continue
      return 1
      end


