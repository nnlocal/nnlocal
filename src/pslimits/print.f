************************************************************************
*     Print momenta and other info on PS point                         *
*                                                                      *
************************************************************************
      subroutine printmom(p,n,verbosity)
      implicit none
      include 'constants.f'
      double precision sqrts
      common/energy/sqrts
      double precision p(4,mxpart),xi(mxpart),yij(mxpart,mxpart),
     .     cosij(mxpart,mxpart),sum(4)
      integer n,i,j,verbosity
      
C---Compute scaled energies, scaled dotproducts, cosines, summed momentum
      sum(:)=0d0
      do i=1,n
         xi(i)=p(4,i)/sqrts
         do j = 1,4
            sum(j)=sum(j)+p(j,i)
         enddo
         do j=1,n
            yij(i,j)=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)
     .           -p(3,i)*p(3,j)
            yij(i,j)=2d0*yij(i,j)/sqrts/sqrts
            cosij(i,j)=p(1,i)*p(1,j)+p(2,i)*p(2,j)+p(3,i)*p(3,j)
            if (((p(1,i)*p(1,i)+p(2,i)*p(2,i)
     .           +p(3,i)*p(3,i)).gt.0d0)
     .           .and.((p(1,j)*p(1,j)+p(2,j)*p(2,j)
     .           +p(3,j)*p(3,j)).gt.0d0)) then
               cosij(i,j)=cosij(i,j)
     .              /sqrt(p(1,i)*p(1,i)+p(2,i)*p(2,i)+p(3,i)*p(3,i))
     .              /sqrt(p(1,j)*p(1,j)+p(2,j)*p(2,j)+p(3,j)*p(3,j))
            endif
         enddo
      enddo
      
      if (verbosity.ge.0) then
C---Write momenta
         write(6,*)'Momenta: p(x),   p(y),   p(z),   p(E),   p^2'
         write(6,*)
         do i=1,n
            write(6,100) p(1,i),p(2,i),p(3,i),p(4,i),
     .           p(4,i)*p(4,i)-p(1,i)*p(1,i)-p(2,i)*p(2,i)
     .           -p(3,i)*p(3,i)
         enddo
         write(6,100) sum(1),sum(2),sum(3),sum(4),
     .        sum(4)*sum(4)-sum(1)*sum(1)-sum(2)*sum(2)-sum(3)*sum(3)
         write(6,*)
      endif
      
      if (verbosity.ge.1) then
C---Write scaled energies
         write(6,*)'Scaled energies xi = pi(E)/E_tot'
         write(6,*)
         do i=1,n
            write(6,200) i,xi(i)
         enddo
         write(6,*)
C---Write scaled invariants
         write(6,*)'Scaled dotproducts yij = 2p_i.p_j/E_tot^2'
         write(6,*)
         do i=1,n
            do j=i+1,n
               write(6,300) i,j,yij(i,j)
            enddo
         enddo
         write(6,*)
C---Write cosines
         write(6,*)'Cosines of three-momenta'
         write(6,*)
         do i=1,n
            do j=i+1,n
               write(6,400) i,j,cosij(i,j)
            enddo
         enddo
         write(6,*)
      endif

 100  format(G16.10,' , ',G16.10,' , ',G16.10,' , ',G16.10,' , ',G16.10)
 200  format('P('I2,') = ',G16.10)
 300  format('yij('I2,' , ',I2') = ',G16.10)
 400  format('cosij('I2,' , ',I2') = ',G16.10)
      end
