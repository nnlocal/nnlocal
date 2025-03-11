************************************************************************
*     Generate PS points tending to any given limit                    *
*     Gabor Somogyi March 2014                                         *
************************************************************************
      subroutine genlimit(limit,p1,ip,rp,jp,sp)
      implicit none
      include 'seed.f'
      include 'npart.f'
      include 'breit.f'
      include 'constants.f'
      
      double precision p(4,mxpart),wgt,masses(mxpart),p1(4,mxpart)
     .     ,ptmp(4,mxpart),sqrts
      common/masses/masses

      integer nrambo,ip,rp,jp,sp,step,nit,indices(4),verbosity,j,k
      parameter(nit=5)

      double precision yir,zr,yjs,zs,phi,yirs,zrs,phi2
      double precision yir0,zr0,yjs0,zs0,yirs0,zrs0

      double precision p_sg_in(100,1:5), p_sg_out(100,1:5)
     .                ,p_sg_tmp(100,1:5)

      character*1 ccont
      character*20 limit

      double precision rn, rnmin, rndel, scale
      parameter(rnmin=2d-1, rndel=6d-1, scale=8d0)
      common/energy/sqrts
      logical first
      data first/.true./
      save first,step
      external rn


      if (first) then
      step=0
      first=.false.
      else
      step=step+1
      endif


C---  Basic data for process set here by hand
c---  Number of particles in the final state is npart
c---  Momenta 1 and 2 are incoming (here with negative energy)
c---  Momenta 3 thorught npart+2 are outgoing
c---  Momenta 1 and 2 are massless,
c---  Momenta 3 thorugh npart+2 have masses as in masses(1:npart)
c---  Note: when there is only one massive particle in the LO final
c---        state, sqrts is automatically reset to its mass
      masses(:)=0d0
      masses(1)=mass3
      iseed = 12345


      
C---Setup for the limit
      if ((limit.eq.'coll2').or.(limit.eq.'soft1').or.
     &     (limit.eq.'icoll2')) then
         nrambo = npart-1
      else
         nrambo = npart-2
      endif

      
C---Generate starting phase space with Rambo
 10   call psbanner(limit,ip,rp,jp,sp)
      if (nrambo.gt.1) then
         call rambo(nrambo,sqrts,masses,ptmp,wgt)
      elseif (nrambo.eq.1) then
         sqrts=masses(1)
         ptmp(4,1) = sqrts
         ptmp(1,1) = 0d0
         ptmp(2,1) = 0d0
         ptmp(3,1) = 0d0
      endif
      do j=1,mxpart-2
         do k=1,4
            p(k,j+2) = ptmp(k,j)
         enddo
      enddo
      p(4,1)=-0.5d0*sqrts
      p(1,1)=0d0
      p(2,1)=0d0
      p(3,1)=-0.5d0*sqrts
      
      p(4,2)=-0.5d0*sqrts
      p(1,2)=0d0
      p(2,2)=0d0
      p(3,2)=+0.5d0*sqrts
c      call printmom(p,nrambo+2,1)
C---Limits
c      do step=0,nit
         write(6,100)step
         write(6,*)'=================================================='
         write(6,*)
C---coll2
         if (limit.eq.'coll2') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
            else
               yir = yir0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-1,(0d0,0d0),p_sg_in,
     &           min(ip,rp),yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=ip
            indices(2)=rp
            call finalizemom(p1,npart+2,2,indices)
C---icoll2
         elseif (limit.eq.'icoll2') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
            else
               yir = yir0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-1,(0d0,0d0),p_sg_in,
     &           ip,yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=rp
            call finalizemom(p1,npart+2,1,indices)            
C---coll3
         elseif (limit.eq.'coll3') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yirs0 = rndel*rn(4)+rnmin
               yirs = yirs0
               zrs = rndel*rn(5)+rnmin
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               yirs = yirs0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           npart+2-2,yir,zr,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           npart+2-1,yirs,zrs,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=ip
            indices(2)=rp
            indices(3)=sp
            call finalizemom(p1,npart+2,3,indices)
C---icoll3
         elseif (limit.eq.'icoll3') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yirs0 = rndel*rn(4)+rnmin
               yirs = yirs0
               zrs = rndel*rn(5)+rnmin
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               yirs = yirs0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           ip,yir,zr,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           npart+2-1,yirs,zrs,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=rp
            indices(2)=sp
            call finalizemom(p1,npart+2,2,indices)            
C---coll22
         elseif (limit.eq.'coll22') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yjs0 = rndel*rn(4)+rnmin
               yjs = yjs0
               zs = rndel*rn(5)+rnmin
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               yjs = yjs0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           npart+2-3,yir,zr,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           npart+2-3,yjs,zs,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=ip
            indices(2)=rp
            indices(3)=jp
            indices(4)=sp
            call finalizemom(p1,npart+2,4,indices)
C---icoll22
         elseif (limit.eq.'icoll22') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yjs0 = rndel*rn(4)+rnmin
               yjs = yjs0
               zs = rndel*rn(5)+rnmin
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               yjs = yjs0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           npart+2-3,yjs,zs,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           ip,yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=jp
            indices(2)=sp
            indices(3)=rp
            call finalizemom(p1,npart+2,3,indices)
C---iicoll22
         elseif (limit.eq.'iicoll22') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yjs0 = rndel*rn(4)+rnmin
               yjs = yjs0
               zs = rndel*rn(5)+rnmin
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               yjs = yjs0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           jp,yjs,zs,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           ip,yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=sp
            indices(2)=rp
            call finalizemom(p1,npart+2,2,indices)            
C---soft1
         elseif (limit.eq.'soft1') then
            if (step.eq.0) then
               yir = rndel*rn(1)+rnmin
               zr0 = rndel*rn(2)+rnmin
               zr = zr0
               phi = TWOPI*rn(3)
            else
               zr = zr0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-1,(0d0,0d0),p_sg_in,
     &           1,yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=rp
            call finalizemom(p1,npart+2,1,indices)
C---soft2
         elseif (limit.eq.'soft2') then
            if (step.eq.0) then
               yir = rndel*rn(1)+rnmin
               zr0 = rndel*rn(2)+rnmin
               zr = zr0
               phi = TWOPI*rn(3)
               yjs = rndel*rn(4)+rnmin
               zs = rndel*rn(5)+rnmin
               phi2 = TWOPI*rn(6)
            else
               zr = zr0*scale**(-step)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           1,yir,zr,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           npart+2-1,yjs,zs,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=rp
            indices(2)=sp
            call finalizemom(p1,npart+2,2,indices)
C---softcoll3
         elseif (limit.eq.'softcoll3') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yjs = rndel*rn(4)+rnmin
               zs0 = rndel*rn(5)+rnmin
               zs = zs0
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               zs = zs0*scale**(-step/2.)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           1,yjs,zs,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           npart+2-2,yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=sp
            indices(2)=ip
            indices(3)=rp
            call finalizemom(p1,npart+2,3,indices)
C---isoftcoll3
         elseif (limit.eq.'isoftcoll3') then
            if (step.eq.0) then
               yir0 = rndel*rn(1)+rnmin
               yir = yir0
               zr = rndel*rn(2)+rnmin
               phi = TWOPI*rn(3)
               yjs = rndel*rn(4)+rnmin
               zs0 = rndel*rn(5)+rnmin
               zs = zs0
               phi2 = TWOPI*rn(6)
            else
               yir = yir0*scale**(-step)
               zs = zs0*scale**(-step/2.)
            endif
            call pTOsingen(p,p_sg_in)
            call singen(npart+2-2,(0d0,0d0),p_sg_in,
     &           1,yjs,zs,p_sg_tmp)
            call singen(npart+2-1,(0d0,0d0),p_sg_tmp,
     &           ip,yir,zr,p_sg_out)
            call pFROMsingen(p_sg_out,p1)
            indices(:)=0
            indices(1)=sp
            indices(2)=rp
            call finalizemom(p1,npart+2,2,indices)            
         endif
C---Phase space is ready, print momenta
c   verbosity.ge.0 prints the set of momenta
c   verbosity.ge.1 prints also scaled energies, scaled dotproducts and 
c   cosines of angles
         verbosity=1
         call printmom(p1,npart+2,verbosity)
c      enddo
C---Do another point?
c      write(6,*)'Continue y/n?'
c      read(5,*)ccont
c      if (ccont.eq.'y') then
c         goto 10
c      endif
      
 100  format(' Iteration no.: ',I2)
      end

      subroutine psbanner(limit,i,r,j,s)
      implicit none
      character*20 limit
      integer i,r,j,s

      write(6,*)
      write(6,*)'**********************************************'
      write(6,*)'**********************************************'
      write(6,*)'*                                            *'
      if (limit.eq.'coll2') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   FF SINGLE COLLINEAR limit:               *'
         write(6,*)'*                                            *'
         write(6,100)i,r
      elseif (limit.eq.'icoll2') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   IF SINGLE COLLINEAR limit:               *'
         write(6,*)'*                                            *'
         write(6,100)i,r
      elseif (limit.eq.'coll3') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   FFF TRIPLE COLLINEAR limit:              *'
         write(6,*)'*                                            *'
         write(6,200)i,r,s
      elseif (limit.eq.'icoll3') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   IFF TRIPLE COLLINEAR limit:              *'
         write(6,*)'*                                            *'
         write(6,200)i,r,s
      elseif (limit.eq.'coll22') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   FF,FF DOUBLE COLLINEAR limit:            *'
         write(6,*)'*                                            *'
         write(6,300)i,r,j,s
      elseif (limit.eq.'icoll22') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   IF,FF DOUBLE COLLINEAR limit:            *'
         write(6,*)'*                                            *'
         write(6,300)i,r,j,s
      elseif (limit.eq.'iicoll22') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   IF,IF DOUBLE COLLINEAR limit:            *'
         write(6,*)'*                                            *'
         write(6,300)i,r,j,s
      elseif (limit.eq.'soft1') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   F SINGLE SOFT limit:                     *'
         write(6,*)'*                                            *'
         write(6,400)r
      elseif (limit.eq.'soft2') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   FF DOUBLE SOFT limit:                    *'
         write(6,*)'*                                            *'
         write(6,500)r,s
      elseif (limit.eq.'softcoll3') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   FF,F DOUBLE SOFT COLLINEAR limit:        *'
         write(6,*)'*                                            *'
         write(6,600)i,r,s
      elseif (limit.eq.'isoftcoll3') then 
         write(6,*)'*   Generating sequence of PS points for     *'
         write(6,*)'*   IF,F DOUBLE SOFT COLLINEAR limit:        *'
         write(6,*)'*                                            *'
         write(6,600)i,r,s
      endif
      write(6,*)'*                                            *'
      write(6,*)'**********************************************'
      write(6,*)'**********************************************'
      write(6,*)
      write(6,*)

 100  format(' *   P(',I2,') || P(',I2,')                           *')
 200  format(' *   P(',I2,') || P(',I2,') || P(',I2,
     .       ')                  *')
 300  format(' *   P(',I2,') || P(',I2,')  and  P(',I2,
     .       ') || P(',I2,')      *')
 400  format(' *   P(',I2,') -> 0                               *')
 500  format(' *   P(',I2,') -> 0  and  P(',I2,') -> 0              *')      
 600  format(' *   P(',I2,') || P(',I2,')  and  P(',I2,
     .       ') -> 0          *')
       end


      subroutine finalizemom(p,npart,n,ind)
      implicit none
      include 'constants.f'
      double precision p(4,mxpart),ptmp(4,mxpart)
      integer n,npart,ind(4),i,j,slot
      logical isfilled(mxpart)
      
      ptmp(:,:) = 0d0
      isfilled(:)=.false.
      do j=1,n
         do i=1,4
            ptmp(i,ind(j)) = p(i,npart-n+j)
         enddo
         isfilled(ind(j)) = .true.
      enddo
       
      j=1
      do slot=1,npart
         if (.not.isfilled(slot)) then
            do i=1,4
               ptmp(i,slot) = p(i,j)
            enddo
            j = j+1
         endif
      enddo

      do i=1,4
         do j=1,npart
            p(i,j) = ptmp(i,j)
         enddo
      enddo
      
      end
