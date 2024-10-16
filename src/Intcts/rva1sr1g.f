      subroutine rva1sr1g(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(6)
      include 'rva1sr1_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)=1 - lam

      tmp = 1.0_ki/2.0_ki*CA*z1*Ta2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=lam - 1
      ab(1)=ab(1)*l1
      ab(1)=1 + ab(1)

      tmp = 2*CA*z1*Ta2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=1.0_ki/2.0_ki*lam - 1.0_ki/2.0_ki
      ab(1)=ab(1)*Pi**2
      ab(2)=1 - lam
      ab(2)=l1*ab(2)
      ab(2)= - 2 + ab(2)
      ab(2)=l1*ab(2)
      ab(1)=4*ab(2) - 8*l2 + ab(1)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=Pi**2
      ab(3)= - ab(2)*ab(1)
      ab(4)=l1*ab(1)
      ab(4)=1 + 1.0_ki/3.0_ki*ab(4)
      ab(4)=l1*ab(4)
      ab(3)=8*ab(4) + 16*l2 + ab(3)
      ab(3)=l1*ab(3)
      ab(2)=ab(2) - ab(3)
      ab(3)=9*zeta3
      ab(1)=ab(3)*ab(1)
      ab(3)=l2**2
      ab(1)=16*ab(3) + ab(1) - 2*ab(2)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=lam - 1
      ab(2)=9*zeta3
      ab(3)= - ab(2)*ab(1)
      ab(4)=Pi**2
      ab(5)=ab(1)*ab(4)
      ab(1)= - l1*ab(1)
      ab(1)= - 4 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)=4.0_ki/3.0_ki*ab(1) - 16*l2 + ab(5)
      ab(1)=l1*ab(1)
      ab(6)=l2**2
      ab(1)=ab(1) + 2*ab(4) - 16*ab(6) + ab(3)
      ab(1)=l1*ab(1)
      ab(3)=l2**3
      ab(1)=ab(1) - ab(2) - 16.0_ki/3.0_ki*ab(3)
      ab(2)=8*l2 + 7.0_ki/60.0_ki*ab(5)
      ab(2)=ab(2)*ab(4)
      ab(1)=ab(2) + 4*ab(1)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  - 4*CA*z1*z2*Ta2*lam**3
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l4 - l3

      tmp = 16*CA*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - l4**2
      ab(2)=2*l4 - l3
      ab(2)=l3*ab(2)
      ab(1)=ab(1) + ab(2)
      ab(2)=Pi**2
      ab(1)=8*ab(1) + ab(2)

      tmp = 4*CA*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=Pi**2
      ab(2)=l3**2
      ab(3)=ab(1) - 8.0_ki/3.0_ki*ab(2)
      ab(3)=l3*ab(3)
      ab(4)= - l3 + 1.0_ki/3.0_ki*l4
      ab(4)=l4*ab(4)
      ab(2)=ab(2) + ab(4)
      ab(1)= - ab(1) + 8*ab(2)
      ab(1)=l4*ab(1)
      ab(1)=ab(3) + ab(1)
      ab(1)=9*zeta3 + 2*ab(1)

      tmp = 8*CA*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

