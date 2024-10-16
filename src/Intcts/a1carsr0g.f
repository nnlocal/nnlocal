      subroutine a1carsr0g(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(2,2,-2:2)
      complex(ki) ab(4)
      include 'a1carsr0_functions.h'

!##### DOUBLE POLE #####

      
      ab(1)= - 1 + lam

      tmp = z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - lam + 1
      ab(1)=l1*ab(1)
      ab(1)= - 1 + ab(1)

      tmp = 2*z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)= - 1.0_ki/6.0_ki*lam + 1.0_ki/6.0_ki
      ab(1)=ab(1)*Pi**2
      ab(2)= - 1 + lam
      ab(2)=l1*ab(2)
      ab(2)=2 + ab(2)
      ab(2)=l1*ab(2)
      ab(1)=2*ab(2) + 4*l2 + ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1, 0) = real(tmp,ki)

!##### ORDER EPS #####

      
      ab(1)=l2**2
      ab(2)=lam*zeta3
      ab(1)= - ab(2) - 2*ab(1) + zeta3
      ab(2)=lam - 1
      ab(3)=Pi**2
      ab(3)=1.0_ki/3.0_ki*ab(3)
      ab(4)=ab(3)*ab(2)
      ab(2)=l1*ab(2)
      ab(2)= - 1 - 1.0_ki/3.0_ki*ab(2)
      ab(2)=l1*ab(2)
      ab(2)=4*ab(2) - 8*l2 + ab(4)
      ab(2)=l1*ab(2)
      ab(1)=ab(2) + ab(3) + 2*ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1, 1) = real(tmp,ki)

!##### ORDER EPS^2 #####

      
      ab(1)=lam - 1
      ab(2)=zeta3*ab(1)
      ab(3)=l2**2
      ab(2)=2*ab(3) + ab(2)
      ab(3)=Pi**2
      ab(2)=2*ab(2) - 1.0_ki/3.0_ki*ab(3)
      ab(4)=ab(1)*ab(3)
      ab(1)=l1*ab(1)
      ab(1)=4 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)=2.0_ki/3.0_ki*ab(1) + 8*l2 - 1.0_ki/3.0_ki*ab(4)
      ab(1)=l1*ab(1)
      ab(1)=2*ab(2) + ab(1)
      ab(1)=l1*ab(1)
      ab(2)= - 2.0_ki/3.0_ki*l2 - 1.0_ki/40.0_ki*ab(4)
      ab(2)=ab(2)*ab(3)
      ab(3)=l2**3
      ab(3)=zeta3 + 2.0_ki/3.0_ki*ab(3)
      ab(1)=ab(1) + 4*ab(3) + ab(2)

      tmp = z1*Ta2*ab(1)
      res(1,1, 2) = real(tmp,ki)

!##### DOUBLE POLE #####

      

      tmp =  0
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      

      tmp = 4*z1*z2*Ta2*lam**3
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)= - l4 + l3

      tmp = 8*z1*z2*Ta2*lam**3*ab(1)
      res(2,2, 0) = real(tmp,ki)

!##### ORDER EPS #####

      
      ab(1)=l4**2
      ab(2)= - 2*l4 + l3
      ab(2)=l3*ab(2)
      ab(1)=ab(1) + ab(2)
      ab(2)=Pi**2
      ab(1)=4*ab(1) - 1.0_ki/3.0_ki*ab(2)

      tmp = 2*z1*z2*Ta2*lam**3*ab(1)
      res(2,2, 1) = real(tmp,ki)

!##### ORDER EPS^2 #####

      
      ab(1)=l3**2
      ab(1)=4*ab(1)
      ab(2)=Pi**2
      ab(3)=l3 - 1.0_ki/3.0_ki*l4
      ab(3)=l4*ab(3)
      ab(3)=4*ab(3) + 1.0_ki/3.0_ki*ab(2) - ab(1)
      ab(3)=l4*ab(3)
      ab(1)= - ab(2) + ab(1)
      ab(1)=l3*ab(1)
      ab(1)=ab(3) - 2*zeta3 + 1.0_ki/3.0_ki*ab(1)

      tmp = 4*z1*z2*Ta2*lam**3*ab(1)
      res(2,2, 2) = real(tmp,ki)

      return
      end

