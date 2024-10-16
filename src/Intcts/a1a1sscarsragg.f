      subroutine a1a1sscarsragg(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(8)
      include 'a1a1sscarsr_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)=1.0_ki/2.0_ki*CA + 2*Ta2
      ab(2)=lam - 1

      tmp = z1*Ta2*ab(2)*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - CA - 4*Ta2
      ab(2)=lam - 1
      ab(2)=ab(2)*l1
      ab(2)=ab(2) + 1

      tmp = 2*z1*Ta2*ab(2)*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - 1.0_ki/2.0_ki*lam + 1.0_ki/2.0_ki
      ab(2)=Pi**2
      ab(1)=ab(2)*ab(1)
      ab(3)=lam - 1
      ab(4)=ab(3)*l1
      ab(4)=ab(4) + 2
      ab(4)=ab(4)*l1
      ab(1)=4*ab(4) + 8*l2 + ab(1)
      ab(1)=CA*ab(1)
      ab(2)= - ab(2)*ab(3)
      ab(2)=8*ab(4) + 16*l2 + ab(2)
      ab(2)=Ta2*ab(2)
      ab(1)=2*ab(2) + ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=5*zeta3
      ab(2)=ab(1)*ab(2)
      ab(3)=l2**2
      ab(4)=Pi**2
      ab(2)=ab(2) - ab(4) + 8*ab(3)
      ab(3)=1.0_ki/3.0_ki*l1
      ab(3)=ab(1)*ab(3)
      ab(3)=ab(3) + 1
      ab(5)=8*l1
      ab(3)=ab(3)*ab(5)
      ab(1)=ab(1)*ab(4)
      ab(1)=ab(3) - ab(1) + 16*l2
      ab(3)= - l1*ab(1)
      ab(3)=ab(3) - ab(2)
      ab(3)=Ta2*ab(3)
      ab(2)= - CA*ab(2)
      ab(1)= - l1*CA*ab(1)
      ab(1)=4*ab(3) + ab(2) + ab(1)

      tmp = 2*z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=lam - 1
      ab(2)=5*zeta3
      ab(3)=ab(1)*ab(2)
      ab(4)=l2**2
      ab(5)=Pi**2
      ab(3)=ab(3) - ab(5) + 8*ab(4)
      ab(3)=2*ab(3)
      ab(4)=CA*ab(3)
      ab(6)=ab(1)*l1
      ab(6)=ab(6) + 4
      ab(7)=4.0_ki/3.0_ki*l1
      ab(6)=ab(6)*ab(7)
      ab(1)=ab(1)*ab(5)
      ab(7)=16*l2
      ab(6)=ab(6) - ab(1) + ab(7)
      ab(8)=l1*CA*ab(6)
      ab(4)=ab(4) + ab(8)
      ab(4)=l1*ab(4)
      ab(6)=l1*ab(6)
      ab(3)=ab(3) + ab(6)
      ab(3)=l1*ab(3)
      ab(6)=l2**3
      ab(2)=ab(2) + 8.0_ki/3.0_ki*ab(6)
      ab(6)= - ab(7) - 1.0_ki/5.0_ki*ab(1)
      ab(6)=ab(6)*ab(5)
      ab(3)=8*ab(3) + 16*ab(2) + ab(6)
      ab(3)=Ta2*ab(3)
      ab(1)= - 8*l2 - 1.0_ki/10.0_ki*ab(1)
      ab(1)=ab(1)*ab(5)
      ab(1)=8*ab(2) + ab(1)
      ab(1)=CA*ab(1)
      ab(1)=2*ab(3) + ab(1) + 4*ab(4)

      tmp = z1*Ta2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=4*Ta2
      ab(2)=CA + ab(1)

      tmp = z1*z2*lam**3*ab(2)*ab(1)
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=4*Ta2
      ab(2)= - CA - ab(1)
      ab(3)=l4 + l2
      ab(3)= - 5*l3 + 4*ab(3)

      tmp = z1*z2*lam**3*ab(3)*ab(2)*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=2*l4 - 5*l3 + 4*l2
      ab(2)=4*l4
      ab(1)=ab(1)*ab(2)
      ab(2)= - 11*l3 + 20*l2
      ab(2)=ab(2)*l3
      ab(3)=Pi**2
      ab(4)=l2**2
      ab(1)=ab(1) - ab(2) - ab(3) + 8*ab(4)
      ab(2)=4*Ta2
      ab(3)=CA + ab(2)

      tmp = z1*z2*lam**3*ab(3)*ab(2)*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=2*l2
      ab(2)=ab(1) + l4
      ab(3)=10*l4
      ab(2)=ab(2)*ab(3)
      ab(3)=d2 + d1
      ab(4)=Pi**2
      ab(5)=l2**2
      ab(2)=ab(2) + 3*ab(3) - ab(4) + 10*ab(5)
      ab(3)=l4 + l2
      ab(3)= - 23.0_ki/3.0_ki*l3 + 22*ab(3)
      ab(3)=ab(3)*l3
      ab(2)=ab(3) - 2*ab(2)
      ab(3)= - l3*CA*ab(2)
      ab(2)= - l3*ab(2)
      ab(6)=l2 + 1.0_ki/3.0_ki*l4
      ab(7)=8*l4
      ab(6)=ab(6)*ab(7)
      ab(6)=ab(6) - ab(4) + 8*ab(5)
      ab(7)=2*l4
      ab(6)=ab(6)*ab(7)
      ab(4)= - ab(4) + 8.0_ki/3.0_ki*ab(5)
      ab(1)=ab(4)*ab(1)
      ab(4)=t2 + t1
      ab(1)=ab(6) + ab(1) + 7*zeta3 + 12*ab(4)
      ab(2)=ab(2) - ab(1)
      ab(2)=Ta2*ab(2)
      ab(1)= - CA*ab(1)
      ab(1)=4*ab(2) + ab(1) + ab(3)

      tmp = 8*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

