      subroutine a12srsssgg(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(12)
      include 'a12srsss_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)=CA + 4*Ta2
      ab(2)=lam - 1

      tmp = z1*Ta2*ab(2)*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=4*Ta2
      ab(2)= - CA - ab(1)
      ab(3)=lam - 1
      ab(3)=ab(3)*l1
      ab(3)=ab(3) + 1

      tmp = z1*ab(3)*ab(2)*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - 1.0_ki/3.0_ki*lam + 1.0_ki/3.0_ki
      ab(2)=Pi**2
      ab(1)=ab(2)*ab(1)
      ab(3)=lam - 1
      ab(4)=ab(3)*l1
      ab(4)=ab(4) + 2
      ab(4)=ab(4)*l1
      ab(1)=2*ab(4) + 4*l2 + ab(1)
      ab(1)=CA*ab(1)
      ab(2)= - ab(2)*ab(3)
      ab(2)=8*ab(4) + 16*l2 + ab(2)
      ab(2)=Ta2*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=Pi**2
      ab(3)=ab(2)*ab(1)
      ab(4)=1.0_ki/3.0_ki*l1
      ab(4)=ab(1)*ab(4)
      ab(4)=ab(4) + 1
      ab(5)=8*l1
      ab(6)= - ab(4)*ab(5)
      ab(3)=ab(6) - 16*l2 + ab(3)
      ab(3)=l1*ab(3)
      ab(1)=zeta3*ab(1)
      ab(6)=l2**2
      ab(3)=ab(3) + ab(2) - 8*ab(6) - 5*ab(1)
      ab(3)=Ta2*ab(3)
      ab(1)=8.0_ki/3.0_ki*ab(2) - 16*ab(6) - 15*ab(1)
      ab(1)=CA*ab(1)
      ab(6)=1.0_ki/3.0_ki*lam - 1.0_ki/3.0_ki
      ab(2)=ab(2)*ab(6)
      ab(4)= - l1*ab(4)
      ab(2)=2*ab(4) - 4*l2 + ab(2)
      ab(2)=ab(5)*CA*ab(2)
      ab(1)=8*ab(3) + ab(1) + ab(2)

      tmp = 2*z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=lam - 1
      ab(2)=5*zeta3
      ab(3)=ab(2)*ab(1)
      ab(4)=l2**2
      ab(5)=Pi**2
      ab(3)= - ab(5) + 8*ab(4) + ab(3)
      ab(6)=ab(1)*ab(5)
      ab(7)=ab(1)*l1
      ab(7)=ab(7) + 4
      ab(7)=ab(7)*l1
      ab(8)=16*l2
      ab(9)=4.0_ki/3.0_ki*ab(7) + ab(8) - ab(6)
      ab(9)=l1*ab(9)
      ab(3)=2*ab(3) + ab(9)
      ab(9)=8*l1
      ab(3)=ab(3)*ab(9)
      ab(10)=l2**3
      ab(2)=ab(2) + 8.0_ki/3.0_ki*ab(10)
      ab(8)= - ab(8) - 1.0_ki/5.0_ki*ab(6)
      ab(8)=ab(8)*ab(5)
      ab(2)=ab(3) + 16*ab(2) + ab(8)
      ab(2)=Ta2*ab(2)
      ab(3)=15*zeta3
      ab(8)=ab(3) + 16.0_ki/3.0_ki*ab(10)
      ab(10)= - 64*l2 - 19.0_ki/30.0_ki*ab(6)
      ab(11)=1.0_ki/3.0_ki*ab(5)
      ab(10)=ab(10)*ab(11)
      ab(8)=8*ab(8) + ab(10)
      ab(8)=CA*ab(8)
      ab(1)=ab(3)*ab(1)
      ab(1)= - 8.0_ki/3.0_ki*ab(5) + 16*ab(4) + ab(1)
      ab(1)=CA*ab(1)
      ab(3)=1.0_ki/3.0_ki*ab(7) + 4*l2 - 1.0_ki/3.0_ki*ab(6)
      ab(3)=l1*CA*ab(3)
      ab(1)=ab(1) + 4*ab(3)
      ab(1)=ab(1)*ab(9)
      ab(1)=4*ab(2) + ab(8) + ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=CA + 4*Ta2

      tmp = 8*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - CA - 4*Ta2
      ab(2)=l4 + l2
      ab(2)= - 5*l3 + 4*ab(2)

      tmp = 8*z1*z2*Ta2*lam**3*ab(2)*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=2*l4 - 5*l3 + 4*l1
      ab(2)=4*l4
      ab(1)=ab(1)*ab(2)
      ab(2)= - 11*l3 + 20*l2
      ab(2)=ab(2)*l3
      ab(1)=ab(1) - ab(2)
      ab(2)=l1**2
      ab(2)=d1 + ab(2) + d5
      ab(3)=l2 + 2*l1
      ab(3)=ab(3)*l2
      ab(4)=Pi**2
      ab(5)=8*ab(3) + 5.0_ki/3.0_ki*ab(4) + ab(1) - 16*ab(2)
      ab(5)=Ta2*ab(5)
      ab(2)=2*ab(3) + 1.0_ki/3.0_ki*ab(4) - 4*ab(2)
      ab(1)=4*ab(2) + ab(1)
      ab(1)=CA*ab(1)
      ab(1)=ab(1) + 4*ab(5)

      tmp = 8*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=CA*Ta2
      ab(2)=Ta2**2
      ab(3)=ab(1) + 4*ab(2)
      ab(4)=lam**3*z1*z2
      ab(3)=ab(3)*ab(4)
      ab(5)=ab(3)*l2
      ab(6)=ab(3)*l3
      ab(7)=4*l4
      ab(8)= - ab(3)*ab(7)
      ab(8)= - 4.0_ki/3.0_ki*ab(5) + 7*ab(6) + ab(8)
      ab(9)=2*l2
      ab(8)=ab(8)*ab(9)
      ab(10)=2*l4
      ab(11)=ab(3)*ab(10)
      ab(11)=ab(6) + ab(11)
      ab(7)=ab(11)*ab(7)
      ab(11)=l3**2
      ab(11)=11*ab(11)
      ab(11)=ab(3)*ab(11)
      ab(7)=ab(8) - ab(11) + ab(7)
      ab(7)=ab(7)*ab(9)
      ab(8)=ab(3)*l4
      ab(9)=5*l3
      ab(12)=ab(3)*ab(9)
      ab(12)=ab(12) - 4.0_ki/3.0_ki*ab(8)
      ab(12)=ab(12)*ab(10)
      ab(11)= - ab(11) + ab(12)
      ab(10)=ab(11)*ab(10)
      ab(8)=ab(8) - ab(6)
      ab(11)=l1 - l4 - l2
      ab(11)=l1*ab(8)*ab(11)
      ab(12)=ab(3)*l3**3
      ab(7)=32*ab(11) + ab(7) + 23.0_ki/3.0_ki*ab(12) + ab(10)
      ab(10)=16.0_ki/3.0_ki*ab(2) + ab(1)
      ab(9)=ab(9)*ab(10)
      ab(10)= - 5*ab(2) - ab(1)
      ab(10)=l4*ab(10)
      ab(11)=2*ab(2) + ab(1)
      ab(11)=l2*ab(11)
      ab(9)=8.0_ki/3.0_ki*ab(11) + ab(9) + 16.0_ki/3.0_ki*ab(10)
      ab(9)=Pi**2*ab(4)*ab(9)
      ab(8)=64*ab(8)
      ab(10)=d1 + d5
      ab(8)=ab(8)*ab(10)
      ab(10)=3*ab(6) + 4*ab(5)
      ab(10)=d4*ab(10)
      ab(11)=40*ab(3)
      ab(12)= - t1 - t2
      ab(11)=ab(11)*ab(12)
      ab(1)= - 7*ab(2) - 3*ab(1)
      ab(1)=zeta3*ab(4)*ab(1)
      ab(2)=4*ab(3)
      ab(2)=t3*ab(2)
      ab(3)=d3*ab(6)
      ab(4)=d2*ab(5)
      ab(1)=16*ab(4) + 12*ab(3) + ab(2) + 8*ab(1) + 4*ab(10) + 2*ab(7)
     &  + ab(9) + ab(11) + ab(8)

      tmp = 8*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

