      subroutine AGTBa12carsssggg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      complex(ki) ab(70)
      include 'AGTBa12carsss_functions.h'

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(1,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xb - xa

      tmp = 8*CA**2*l1*z6*z7*z9*xa**3*ab(1)
      res(1,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l2 + l4 - l3
      ab(2)= - xb*ab(1)
      ab(2)= - 1 + ab(2)
      ab(3)=5*l1
      ab(4)=xb*ab(3)
      ab(2)=8*ab(2) + ab(4)
      ab(2)=l1*ab(2)
      ab(1)=8*ab(1) - ab(3)
      ab(1)=l1*ab(1)
      ab(3)=Pi**2
      ab(3)=1.0_ki/3.0_ki*ab(3)
      ab(1)=ab(3) + ab(1)
      ab(1)=xa*ab(1)
      ab(3)= - xb*ab(3)
      ab(1)=ab(1) + ab(3) + ab(2)

      tmp = 2*CA**2*z6*z7*z9*xa**3*ab(1)
      res(1,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)= - l3 + 2*l2
      ab(2)=16*l3
      ab(1)=ab(1)*ab(2)
      ab(2)=l3 - l2
      ab(3)= - l4 + 2*ab(2)
      ab(4)=16*l4
      ab(3)=ab(3)*ab(4)
      ab(4)=l2**2
      ab(5)=Pi**2
      ab(1)= - 16*ab(4) + 1.0_ki/3.0_ki*ab(5) + ab(1) + ab(3)
      ab(3)= - xb*ab(1)
      ab(2)=ab(2) - l4
      ab(4)=xb*ab(2)
      ab(4)= - 1 + ab(4)
      ab(6)=2*l1
      ab(7)=xb*ab(6)
      ab(4)=5*ab(4) + ab(7)
      ab(7)=4*l1
      ab(4)=ab(4)*ab(7)
      ab(8)=l4 + l2
      ab(9)= - 2*l3 + l5 + ab(8)
      ab(3)=ab(4) + 32*ab(9) + ab(3)
      ab(3)=l1*ab(3)
      ab(2)= - 5*ab(2) - ab(6)
      ab(2)=ab(2)*ab(7)
      ab(1)=ab(2) + ab(1)
      ab(1)=l1*ab(1)
      ab(2)=ab(8) - l3
      ab(2)=4.0_ki/3.0_ki*ab(2)
      ab(2)=ab(2)*ab(5)
      ab(2)=ab(2) + 17.0_ki/2.0_ki*zeta3
      ab(1)=ab(1) - ab(2)
      ab(1)=xa*ab(1)
      ab(2)=xb*ab(2)
      ab(1)=ab(1) + ab(3) + 4.0_ki/3.0_ki*ab(5) + ab(2)

      tmp = CA**2*z6*z7*z9*xa**3*ab(1)
      res(1,0,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 6*CA**2*z1*z2*xb**3
      res(0,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - l6 - l8 + l7
      ab(1)=3*ab(1) + 2*l1

      tmp = 4*CA**2*z1*z2*xb**3*ab(1)
      res(0,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l8**2
      ab(2)= - 2*l8 + l7
      ab(2)=l7*ab(2)
      ab(1)=ab(1) + ab(2)
      ab(2)=CA**2*xb**3*z1*z2
      ab(1)=ab(2)*ab(1)
      ab(3)=l7 - l8
      ab(4)= - 2*ab(3) + l6
      ab(4)=l6*ab(2)*ab(4)
      ab(1)=ab(4) + ab(1)
      ab(3)= - l6 + ab(3)
      ab(3)=8*ab(3) + 5*l1
      ab(3)=l1*ab(2)*ab(3)
      ab(2)=Pi**2*ab(2)
      ab(1)=2.0_ki/3.0_ki*ab(2) + 6*ab(1) + ab(3)

      tmp = 2*ab(1)
      res(0,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=3*l6
      ab(2)= - ab(1) + 4*l1
      ab(3)=2*l6
      ab(2)=ab(2)*ab(3)
      ab(4)=Pi**2
      ab(5)=l1**2
      ab(6)=2.0_ki/3.0_ki*ab(4) + 5*ab(5)
      ab(2)=ab(2) - ab(6)
      ab(7)=2*l1
      ab(1)=ab(7) - ab(1)
      ab(8)=3*l7
      ab(9)= - 2*ab(1) - ab(8)
      ab(10)=2*l7
      ab(9)=ab(9)*ab(10)
      ab(8)= - l8 + ab(8) + ab(1)
      ab(8)=l8*ab(8)
      ab(8)=2*ab(8) + ab(9) + ab(2)
      ab(8)=l8*ab(8)
      ab(7)=ab(7) - l6
      ab(3)=ab(7)*ab(3)
      ab(3)=ab(3) - ab(6)
      ab(3)=l6*ab(3)
      ab(1)=l7 + ab(1)
      ab(1)=ab(1)*ab(10)
      ab(1)=ab(1) - ab(2)
      ab(1)=l7*ab(1)
      ab(1)=ab(8) + ab(3) + ab(1)
      ab(2)=5.0_ki/3.0_ki*ab(4) + 8*ab(5)
      ab(2)=l1*ab(2)
      ab(1)= - 17.0_ki/2.0_ki*zeta3 + ab(2) + 4*ab(1)

      tmp = CA**2*z1*z2*xb**3*ab(1)
      res(0,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      
      ab(1)=3 - xb
      ab(2)=xa - 1
      ab(1)=xb*ab(2)*ab(1)
      ab(2)=xa - 3
      ab(2)=1.0_ki/2.0_ki*ab(2)
      ab(3)=xa**2*ab(2)
      ab(1)=1.0_ki/2.0_ki*ab(1) + 1 + ab(3)
      ab(1)=xb*ab(1)
      ab(2)= - xa*ab(2)
      ab(2)= - 1 + ab(2)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 3*CA**2*z1*z6*z7*z8*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=l1 - 1
      ab(2)=xa*l1
      ab(3)=ab(2) - ab(1)
      ab(3)=xb*ab(3)
      ab(4)=3*l1
      ab(5)=ab(4) - 4
      ab(4)=ab(4) - 1
      ab(6)= - xa*ab(4)
      ab(3)=ab(3) + ab(6) + ab(5)
      ab(3)=xb*ab(3)
      ab(2)= - ab(2) + ab(4)
      ab(2)=ab(2)*xa**2
      ab(4)=l1 - 2
      ab(4)=2*ab(4)
      ab(2)=ab(3) - ab(4) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) - ab(5)
      ab(1)=xa*ab(1)
      ab(1)=ab(4) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + ab(2)

      tmp = 3*CA**2*z1*z6*z7*z8*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=Pi**2
      ab(2)=7.0_ki/2.0_ki*ab(1)
      ab(3)=12*l2
      ab(4)=9*l1
      ab(5)=10 - ab(4)
      ab(5)=l1*ab(5)
      ab(5)=ab(5) + ab(3) + ab(2)
      ab(5)=xb*ab(5)
      ab(6)=7.0_ki/6.0_ki*ab(1)
      ab(7)=l1**2
      ab(7)= - ab(6) + 3*ab(7)
      ab(7)=ab(7)*xb
      ab(8)=l7 + l6
      ab(7)=ab(7) + 6*ab(8)
      ab(9)=3*l1
      ab(10)=2 - ab(9)
      ab(10)=l1*ab(10)
      ab(10)=ab(10) + ab(6) + ab(7)
      ab(10)=xa*ab(10)
      ab(11)= - 1 - 3*l2
      ab(11)=2*ab(11) - 9*ab(8)
      ab(12)= - 16 + ab(4)
      ab(12)=l1*ab(12)
      ab(5)=ab(10) + ab(5) + ab(12) + 2*ab(11) - ab(2)
      ab(5)=xa*ab(5)
      ab(10)=1 + 3*ab(8)
      ab(11)=6*l2
      ab(12)=ab(10) + ab(11)
      ab(13)=l1 - 4
      ab(14)=6*l1
      ab(13)=ab(13)*ab(14)
      ab(1)= - ab(13) + 7.0_ki/3.0_ki*ab(1) + 4*ab(12)
      ab(12)= - 2 + ab(4)
      ab(12)=l1*ab(12)
      ab(7)=ab(12) - ab(2) - ab(7)
      ab(7)=xb*ab(7)
      ab(10)= - 4*l1 - ab(11) + ab(10)
      ab(7)=4*ab(10) + ab(7)
      ab(7)=xb*ab(7)
      ab(5)=ab(5) + ab(7) + ab(1)
      ab(5)=xa*ab(5)
      ab(7)= - 10 + ab(9)
      ab(7)=l1*ab(7)
      ab(3)=ab(7) - ab(3) - ab(6)
      ab(3)=xb*ab(3)
      ab(6)=ab(11) + ab(8)
      ab(4)=32 - ab(4)
      ab(4)=l1*ab(4)
      ab(2)=ab(3) + ab(4) + 6*ab(6) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) - ab(1)
      ab(1)=xb*ab(1)
      ab(1)=ab(1) + ab(5)

      tmp = CA**2*z1*z6*z7*z8*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=3*ab(3)
      ab(5)=xa*ab(3)
      ab(6)=ab(4) - ab(5)
      ab(7)=ab(6)*xa
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(7) - ab(1)
      ab(2)=ab(2)*xa
      ab(8)=ab(5) - ab(3)
      ab(9)=ab(8)*xb
      ab(10)=xb - 2
      ab(10)=ab(10)*ab(9)
      ab(10)=ab(2) + ab(10)
      ab(11)=3*l6
      ab(11)=ab(10)*ab(11)
      ab(12)=2*l1
      ab(13)=ab(10)*ab(12)
      ab(14)= - 15*ab(3) + 7*ab(5)
      ab(14)=xa*ab(14)
      ab(15)=8*ab(3)
      ab(13)=ab(11) + ab(13) + ab(9) + ab(15) + ab(14)
      ab(13)=l6*ab(13)
      ab(14)=l1*ab(10)
      ab(11)=ab(14) + ab(11)
      ab(10)=l7*ab(10)
      ab(10)=2*ab(11) + 3*ab(10)
      ab(10)=l7*ab(10)
      ab(11)=9*ab(3) - ab(5)
      ab(11)=xa*ab(11)
      ab(11)= - 7*ab(9) - ab(15) + ab(11)
      ab(11)=l9*ab(11)
      ab(10)=ab(11) + ab(13) + ab(10)
      ab(11)=xb*ab(4)
      ab(11)=ab(11) - 10*ab(3) + ab(5)
      ab(11)=xb*ab(11)
      ab(13)=4*ab(3) - 3*ab(5)
      ab(13)=xa*ab(13)
      ab(11)=ab(11) + ab(15) + ab(13)
      ab(11)=xb*ab(11)
      ab(13)=6*ab(3) - ab(5)
      ab(13)=xa*ab(13)
      ab(13)= - ab(15) + ab(13)
      ab(13)=xa*ab(13)
      ab(11)=ab(13) + ab(11)
      ab(13)=xb - 3
      ab(13)=ab(13)*ab(9)
      ab(14)=xa**2
      ab(6)=ab(6)*ab(14)
      ab(6)=ab(13) + ab(6)
      ab(13)= - ab(1) + ab(6)
      ab(13)=xb*ab(13)
      ab(2)= - ab(2) + ab(13)
      ab(2)=ab(2)*ab(12)
      ab(2)=3*ab(11) + ab(2)
      ab(2)=l1*ab(2)
      ab(11)=ab(1) - ab(5)
      ab(11)=ab(11)*xa
      ab(12)=ab(3)*xb
      ab(13)= - ab(1) + ab(12)
      ab(13)=xb*ab(13)
      ab(13)=ab(13) + ab(1) - ab(11)
      ab(4)=ab(12) - ab(4)
      ab(4)=ab(4)*xb
      ab(4)=ab(4) + ab(11) + ab(1)
      ab(4)=ab(4)*xb
      ab(4)=ab(4) - ab(11)
      ab(11)=8*ab(4)
      ab(12)=l3*ab(11)
      ab(2)=ab(12) + 8*ab(13) + ab(2)
      ab(2)=l1*ab(2)
      ab(12)= - xb*ab(1)
      ab(12)=ab(12) + 5*ab(3) + ab(5)
      ab(13)=1.0_ki/3.0_ki*xb
      ab(12)=ab(12)*ab(13)
      ab(5)= - ab(3) + 1.0_ki/3.0_ki*ab(5)
      ab(15)=ab(5)*xa
      ab(16)= - 1.0_ki/3.0_ki*ab(3) + ab(15)
      ab(12)=2*ab(16) + ab(12)
      ab(12)=xb*ab(12)
      ab(16)= - xa*ab(8)
      ab(1)=ab(1) + ab(16)
      ab(1)=xa*ab(1)
      ab(1)=1.0_ki/3.0_ki*ab(1) + ab(12)
      ab(5)=ab(5)*ab(14)
      ab(12)= - ab(13) + 1
      ab(9)=ab(12)*ab(9)
      ab(12)=2.0_ki/3.0_ki*ab(3)
      ab(5)=ab(9) + ab(12) + ab(5)
      ab(5)=xb*ab(5)
      ab(9)= - ab(12) - ab(15)
      ab(9)=xa*ab(9)
      ab(5)=ab(9) + ab(5)
      ab(5)=l1*ab(5)
      ab(1)=2*ab(1) + 7*ab(5)
      ab(1)=ab(1)*Pi**2
      ab(5)=l1*ab(11)
      ab(4)=l2*ab(4)
      ab(9)=xb - xa
      ab(8)=ab(9)*ab(8)
      ab(4)=6*ab(4) - 3*ab(8) + ab(5)
      ab(4)=l2*ab(4)
      ab(5)= - ab(3) + 1.0_ki/2.0_ki*ab(6)
      ab(5)=xb*ab(5)
      ab(3)=ab(3) - 1.0_ki/2.0_ki*ab(7)
      ab(3)=xa*ab(3)
      ab(3)=ab(3) + ab(5)
      ab(3)=zeta3*ab(3)
      ab(5)=l5*ab(8)

      tmp = ab(1) + ab(2) + 97.0_ki/2.0_ki*ab(3) + 4*ab(4) + 8*ab(5) + 2*
     & ab(10)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=xa*ab(3)
      ab(5)= - ab(4) + 9*ab(3)
      ab(5)=ab(5)*xa
      ab(6)=ab(4) - ab(3)
      ab(7)=ab(6)*xb
      ab(8)=7*ab(7)
      ab(9)=8*ab(3)
      ab(5)=ab(5) - ab(8) - ab(9)
      ab(10)=2*l1
      ab(11)= - 3*l2 - ab(10)
      ab(11)=ab(5)*ab(11)
      ab(12)=11*ab(7)
      ab(13)=13*ab(3)
      ab(14)= - ab(13) + ab(4)
      ab(14)=xa*ab(14)
      ab(15)=12*ab(3)
      ab(14)=ab(12) + ab(15) + ab(14)
      ab(14)=l6*ab(14)
      ab(16)=7*ab(3)
      ab(17)= - ab(16) + ab(4)
      ab(17)=xa*ab(17)
      ab(18)=6*ab(3)
      ab(17)=5*ab(7) + ab(18) + ab(17)
      ab(17)=l9*ab(17)
      ab(11)=ab(17) + ab(14) + ab(11)
      ab(11)=l9*ab(11)
      ab(14)=3*ab(3)
      ab(17)= - xb*ab(14)
      ab(17)=ab(17) + 10*ab(3) - ab(4)
      ab(17)=xb*ab(17)
      ab(19)=3*ab(4)
      ab(20)=4*ab(3)
      ab(21)= - ab(20) + ab(19)
      ab(21)=xa*ab(21)
      ab(17)=ab(17) - ab(9) + ab(21)
      ab(17)=xb*ab(17)
      ab(21)=1.0_ki/3.0_ki*xb
      ab(22)=ab(21) - 1
      ab(22)=ab(22)*ab(7)
      ab(23)= - ab(3) + 1.0_ki/3.0_ki*ab(4)
      ab(24)=xa**2
      ab(25)=ab(23)*ab(24)
      ab(22)=ab(22) - ab(25)
      ab(25)=2.0_ki/3.0_ki*ab(3)
      ab(26)= - ab(25) + ab(22)
      ab(26)=ab(26)*xb
      ab(23)=ab(23)*xa
      ab(25)=ab(23) + ab(25)
      ab(25)=ab(25)*xa
      ab(26)=ab(26) + ab(25)
      ab(27)= - ab(26)*ab(10)
      ab(18)= - ab(18) + ab(4)
      ab(18)=xa*ab(18)
      ab(18)=ab(9) + ab(18)
      ab(18)=xa*ab(18)
      ab(17)=ab(27) + ab(18) + ab(17)
      ab(17)=l1*ab(17)
      ab(18)=5*ab(4)
      ab(27)=ab(3)*xb
      ab(28)= - ab(27) + ab(16) - ab(18)
      ab(28)=xb*ab(28)
      ab(29)=23*ab(3) - 4*ab(4)
      ab(29)=xa*ab(29)
      ab(17)=ab(17) + 3*ab(28) - 22*ab(3) + ab(29)
      ab(28)=l1**2
      ab(17)=ab(17)*ab(28)
      ab(29)=ab(27) - ab(14)
      ab(29)=ab(29)*xb
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(1) - ab(4)
      ab(2)=ab(2)*xa
      ab(30)=ab(2) + ab(1)
      ab(29)=ab(29) + ab(30)
      ab(29)=ab(29)*xb
      ab(29)=ab(29) - ab(2)
      ab(28)=3*ab(28)
      ab(31)=8*l1
      ab(32)= - l2*ab(31)
      ab(33)=4*l1
      ab(34)= - l3*ab(33)
      ab(32)=ab(34) + ab(32) - ab(28)
      ab(32)=ab(29)*ab(32)
      ab(34)=ab(3) - 1.0_ki/3.0_ki*ab(27)
      ab(34)=xb*ab(34)
      ab(30)= - 1.0_ki/3.0_ki*ab(30) + ab(34)
      ab(30)=xb*ab(30)
      ab(2)=1.0_ki/3.0_ki*ab(2) + ab(30)
      ab(30)=Pi**2
      ab(2)=ab(2)*ab(30)
      ab(2)=ab(2) + ab(32)
      ab(2)=l3*ab(2)
      ab(32)=5*ab(3)
      ab(34)= - ab(32) - ab(19)
      ab(34)=xa*ab(34)
      ab(12)=ab(12) + ab(9) + ab(34)
      ab(12)=d1*ab(12)
      ab(2)=ab(17) + ab(2) + ab(11) + ab(12)
      ab(11)=1.0_ki/2.0_ki*ab(4)
      ab(12)=17.0_ki/6.0_ki*ab(27) - ab(9) - ab(11)
      ab(12)=xb*ab(12)
      ab(17)=20*ab(3) - 17.0_ki/2.0_ki*ab(4)
      ab(17)=xa*ab(17)
      ab(17)=14*ab(3) + ab(17)
      ab(12)=1.0_ki/3.0_ki*ab(17) + ab(12)
      ab(12)=xb*ab(12)
      ab(17)=ab(26)*ab(31)
      ab(11)=4.0_ki/3.0_ki*ab(3) + ab(11)
      ab(11)=xa*ab(11)
      ab(11)= - 14.0_ki/3.0_ki*ab(3) + ab(11)
      ab(11)=xa*ab(11)
      ab(11)=ab(17) + ab(11) + ab(12)
      ab(11)=l1*ab(11)
      ab(12)=ab(1)*xb
      ab(17)= - ab(12) + 17*ab(3) - 13*ab(4)
      ab(17)=xb*ab(17)
      ab(26)=25*ab(3)
      ab(31)= - ab(26) + 19*ab(4)
      ab(31)=xa*ab(31)
      ab(17)=ab(17) + ab(20) + ab(31)
      ab(31)= - ab(7) + 2*ab(6)
      ab(21)=ab(31)*ab(21)
      ab(21)=ab(21) + ab(25)
      ab(25)=4*l6
      ab(34)=ab(21)*ab(25)
      ab(35)=1.0_ki/3.0_ki*ab(3)
      ab(22)= - ab(35) + 1.0_ki/2.0_ki*ab(22)
      ab(22)=xb*ab(22)
      ab(23)=ab(35) + 1.0_ki/2.0_ki*ab(23)
      ab(23)=xa*ab(23)
      ab(22)=ab(23) + ab(22)
      ab(22)=ab(22)*ab(30)
      ab(11)=77.0_ki/24.0_ki*ab(22) + ab(34) + 1.0_ki/3.0_ki*ab(17) + ab(11)
      ab(11)=ab(11)*ab(30)
      ab(17)=ab(6)*xa
      ab(17)=ab(17) - ab(7)
      ab(22)=l5 - l2
      ab(22)=2*ab(22)
      ab(17)=ab(17)*ab(22)
      ab(22)= - ab(20) + ab(4)
      ab(22)=xa*ab(22)
      ab(23)=2*ab(4) - ab(27)
      ab(23)=xb*ab(23)
      ab(22)=ab(23) + ab(1) + ab(22)
      ab(22)=ab(22)*ab(33)
      ab(23)=7*ab(4)
      ab(34)= - ab(23) + 15*ab(3)
      ab(34)=ab(34)*xa
      ab(34)=ab(34) - ab(7) - ab(9)
      ab(35)= - l6*ab(34)
      ab(5)=l9*ab(5)
      ab(5)=ab(5) + ab(22) + ab(35) + ab(17)
      ab(5)=l5*ab(5)
      ab(17)=ab(31)*xb
      ab(22)=ab(14) - ab(4)
      ab(31)=ab(22)*xa
      ab(35)=ab(31) - ab(1)
      ab(36)=ab(35)*xa
      ab(17)=ab(36) - ab(17)
      ab(37)=l1*ab(17)
      ab(34)=2*ab(34) - 3*ab(37)
      ab(34)=ab(34)*ab(10)
      ab(25)= - ab(25) - ab(33)
      ab(25)=ab(17)*ab(25)
      ab(33)=53*ab(3) - 25*ab(4)
      ab(33)=xa*ab(33)
      ab(25)= - 3*ab(7) - 28*ab(3) + ab(33) + ab(25)
      ab(25)=l6*ab(25)
      ab(25)=ab(34) + ab(25)
      ab(25)=l6*ab(25)
      ab(12)= - ab(12) + 11*ab(3) - ab(23)
      ab(12)=xb*ab(12)
      ab(13)=ab(13) - ab(19)
      ab(13)=xa*ab(13)
      ab(12)=ab(12) - ab(15) + ab(13)
      ab(13)=l1*ab(29)
      ab(12)=4*ab(12) - 15*ab(13)
      ab(12)=l1*ab(12)
      ab(13)=ab(26) - 9*ab(4)
      ab(13)=xa*ab(13)
      ab(8)= - ab(8) - 16*ab(3) + ab(13)
      ab(13)=2*l6
      ab(8)=ab(8)*ab(13)
      ab(15)=ab(29)*ab(30)
      ab(8)=3*ab(15) + ab(12) + ab(8)
      ab(12)= - 32*l2 - 56*l1
      ab(12)=ab(29)*ab(12)
      ab(15)= - ab(14) - ab(18)
      ab(15)=xa*ab(15)
      ab(9)=13*ab(7) + ab(9) + ab(15) + ab(12)
      ab(9)=l2*ab(9)
      ab(8)=2*ab(8) + ab(9)
      ab(8)=l2*ab(8)
      ab(9)=3*l6
      ab(10)= - ab(10) - ab(9)
      ab(10)=ab(13)*ab(10)
      ab(12)=2*l7
      ab(9)= - l7 - l1 - ab(9)
      ab(9)=ab(12)*ab(9)
      ab(9)=ab(9) - ab(28) + ab(10)
      ab(9)=ab(17)*ab(9)
      ab(10)=ab(21)*ab(30)
      ab(9)=2*ab(10) + ab(9)
      ab(9)=ab(9)*ab(12)
      ab(10)=1.0_ki/4.0_ki*ab(4)
      ab(12)= - 9.0_ki/4.0_ki*ab(27) + ab(16) - ab(10)
      ab(12)=xb*ab(12)
      ab(4)= - ab(20) + 9.0_ki/4.0_ki*ab(4)
      ab(4)=xa*ab(4)
      ab(4)=ab(12) - ab(32) + ab(4)
      ab(4)=xb*ab(4)
      ab(10)= - ab(14) + ab(10)
      ab(10)=xa*ab(10)
      ab(10)=ab(32) + ab(10)
      ab(10)=xa*ab(10)
      ab(4)=ab(10) + ab(4)
      ab(10)= - ab(7) + 3*ab(6)
      ab(10)=ab(10)*xb
      ab(12)=ab(22)*ab(24)
      ab(10)=ab(10) - ab(12)
      ab(12)=ab(3) + 1.0_ki/2.0_ki*ab(10)
      ab(12)=xb*ab(12)
      ab(3)= - ab(3) + 1.0_ki/2.0_ki*ab(31)
      ab(3)=xa*ab(3)
      ab(3)=ab(3) + ab(12)
      ab(3)=l1*ab(3)
      ab(3)=17*ab(4) + 111*ab(3)
      ab(3)=zeta3*ab(3)
      ab(1)=ab(1) + ab(10)
      ab(1)=xb*ab(1)
      ab(1)=ab(36) + ab(1)
      ab(1)=plg4half*ab(1)
      ab(4)=ab(7) - ab(35)
      ab(4)=d3*ab(4)
      ab(6)= - ab(7) + ab(6)
      ab(6)=d2*ab(6)

      tmp = 8*ab(1) + 2*ab(2) + ab(3) + 32*ab(4) + 4*ab(5) + 64*ab(6)
     &  + ab(8) + ab(9) + ab(11) + ab(25)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - 2 + xa
      ab(1)=xa*ab(1)
      ab(1)=3 + ab(1)
      ab(1)=xa*ab(1)
      ab(1)= - 2 + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=1 + ab(1)

      tmp = 3*CA**2*z6*z10*ab(1)
      res(2,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=CA**2*z6*z10
      ab(2)=ab(1)*xa
      ab(3)=2*ab(1)
      ab(2)=ab(2) - ab(3)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) + 3*ab(1)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - ab(3)
      ab(2)=ab(2)*xa
      ab(1)=ab(2) + ab(1)
      ab(2)= - 5*l1 - 6*l2 + 2*l4 + l3

      tmp = 2*ab(2)*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=8*l2
      ab(2)= - l3 + ab(1) + 7*l1
      ab(3)=3*l4
      ab(2)= - ab(3) + 2*ab(2)
      ab(4)=ab(2)*l4
      ab(5)=3*l2 + 5*l1
      ab(1)=ab(5)*ab(1)
      ab(6)=4*l2 + 3*l1
      ab(6)= - l3 + 2*ab(6)
      ab(6)=ab(6)*l3
      ab(7)=Pi**2
      ab(8)=l1**2
      ab(1)=ab(4) + ab(6) - ab(1) + 5.0_ki/3.0_ki*ab(7) - 17*ab(8)
      ab(4)=2*ab(1)
      ab(9)= - xa*ab(1)
      ab(9)=ab(4) + ab(9)
      ab(9)=xa*ab(9)
      ab(2)= - ab(2)*ab(3)
      ab(3)=l2*ab(5)
      ab(2)=ab(9) + ab(2) - 3*ab(6) + 24*ab(3) - 5*ab(7) + 51*ab(8)
      ab(2)=xa*ab(2)
      ab(2)=ab(4) + ab(2)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) - ab(1)

      tmp = CA**2*z6*z10*ab(1)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=CA**2*z6*z10
      ab(2)=ab(1)*xa
      ab(3)=2*ab(1)
      ab(2)=ab(2) - ab(3)
      ab(4)=ab(2)*xa
      ab(4)=ab(4) + 3*ab(1)
      ab(4)=ab(4)*xa
      ab(3)=ab(4) - ab(3)
      ab(3)=ab(3)*xa
      ab(3)=ab(3) + ab(1)
      ab(4)=ab(3)*l4
      ab(5)=1.0_ki/3.0_ki*xa
      ab(2)=ab(2)*ab(5)
      ab(2)=ab(2) + ab(1)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - 2.0_ki/3.0_ki*ab(1)
      ab(2)=ab(2)*xa
      ab(1)=ab(2) + 1.0_ki/3.0_ki*ab(1)
      ab(2)=ab(1)*l1
      ab(5)=ab(3)*l2
      ab(6)=ab(3)*l3
      ab(7)= - 59*ab(2) - 68*ab(5) + 25*ab(4) + 9*ab(6)
      ab(7)=l1*ab(7)
      ab(8)=2*l4
      ab(8)=ab(3)*ab(8)
      ab(8)=ab(8) + ab(6)
      ab(9)=3*l3
      ab(10)= - ab(8)*ab(9)
      ab(11)=7*l4
      ab(9)=ab(11) + ab(9)
      ab(9)=ab(3)*ab(9)
      ab(9)= - 10*ab(5) + ab(9)
      ab(9)=l2*ab(9)
      ab(12)=ab(3)*l4**2
      ab(7)=ab(7) + 8*ab(9) - 11*ab(12) + ab(10)
      ab(7)=l1*ab(7)
      ab(9)= - l3*ab(8)
      ab(5)= - 2*ab(5) + ab(8)
      ab(8)=4*l2
      ab(5)=ab(5)*ab(8)
      ab(5)=ab(5) - 3*ab(12) + ab(9)
      ab(5)=ab(5)*ab(8)
      ab(8)=20*l2 - ab(11)
      ab(8)=ab(1)*ab(8)
      ab(2)=17*ab(2) - ab(6) + ab(8)
      ab(2)=ab(2)*Pi**2
      ab(6)=l3*ab(1)
      ab(4)=ab(4) + ab(6)
      ab(4)=l3*ab(4)
      ab(4)=ab(12) + ab(4)
      ab(4)=l3*ab(4)
      ab(1)=ab(1)*l4**3
      ab(3)=zeta3*ab(3)

      tmp = 5*ab(1) + ab(2) - 34*ab(3) + ab(4) + ab(5) + ab(7)
      res(2,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,2,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(2,2,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xa**2
      ab(2)= - 1 + ab(1)
      ab(2)=ab(2)*ab(1)
      ab(1)=ab(1) + 1
      ab(3)=ab(1)*xa**4
      ab(4)= - xb*xa**5
      ab(3)=ab(3) + ab(4)
      ab(3)=xb*ab(3)
      ab(4)= - ab(1)*xa**3
      ab(3)=ab(4) + ab(3)
      ab(3)=xb*ab(3)
      ab(2)=ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(3)=xa*ab(1)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) - ab(1)
      ab(1)=xb*ab(1)
      ab(1)=xa + ab(1)

      tmp = 8*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l8 + l6
      ab(2)=l4 + l2
      ab(2)= - l7 + 3*ab(2)
      ab(3)= - ab(2) - 2 + 2*l5
      ab(3)=2*ab(3)
      ab(4)=9*l3
      ab(5)= - l10 - ab(1) - ab(3) - ab(4) + 7*l11
      ab(6)=xb**2
      ab(7)=4*ab(6)
      ab(8)=ab(7) + ab(5)
      ab(8)=ab(8)*ab(6)
      ab(5)=ab(8) - ab(5)
      ab(5)=ab(5)*ab(6)
      ab(3)=l11 + ab(3) + l3 + 9*l10 - 7*ab(1)
      ab(8)=l11 + 5*l10 - 3*ab(1)
      ab(4)=ab(4) + ab(8)
      ab(9)=2*ab(2)
      ab(10)=ab(9) - ab(4)
      ab(11)= - ab(10)*ab(6)
      ab(11)=ab(11) + ab(3)
      ab(11)=ab(11)*ab(6)
      ab(1)= - ab(1) + l10 + l5
      ab(12)= - ab(1) + 2*l3
      ab(13)=1 - ab(12)
      ab(11)=4*ab(13) + ab(11)
      ab(11)=xb*ab(11)
      ab(2)=ab(2) - 2
      ab(2)=2*ab(2)
      ab(4)=ab(2) - ab(4)
      ab(4)=ab(4)*ab(6)
      ab(10)=ab(4) + ab(10)
      ab(10)=ab(10)*ab(6)
      ab(13)=1 + ab(12)
      ab(10)=4*ab(13) + ab(10)
      ab(10)=ab(10)*ab(6)
      ab(4)= - 4 - ab(4)
      ab(4)=ab(4)*ab(6)
      ab(12)=4*ab(12)
      ab(4)= - ab(12) + ab(4)
      ab(4)=xa*xb*ab(4)
      ab(4)=ab(4) + ab(12) + ab(10)
      ab(4)=xa*ab(4)
      ab(4)=ab(11) + ab(4)
      ab(4)=xa*ab(4)
      ab(4)=ab(4) - 4 + ab(5)
      ab(4)=xa*ab(4)
      ab(5)= - 1 - ab(1)
      ab(5)=ab(5)*ab(7)
      ab(3)=ab(5) - ab(3)
      ab(3)=ab(3)*ab(6)
      ab(5)=ab(8) + l3
      ab(8)=ab(9) - ab(5)
      ab(3)=ab(3) + ab(8)
      ab(3)=xb*ab(3)
      ab(3)=ab(3) + ab(4)
      ab(3)=xa*ab(3)
      ab(4)=ab(1)*ab(6)
      ab(4)=ab(4) - 1
      ab(1)=ab(1) + ab(4)
      ab(1)=ab(1)*ab(7)
      ab(1)=ab(1) - ab(8)
      ab(1)=ab(1)*ab(6)
      ab(2)=ab(2) - ab(5)
      ab(1)=ab(3) + ab(1) - ab(2)
      ab(1)=xa*ab(1)
      ab(3)= - ab(4)*ab(7)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l4 + l2
      ab(2)=l11 - 3
      ab(3)=2*ab(2)
      ab(4)=ab(3) - l5
      ab(5)=3*l3
      ab(6)=ab(5) + 7*l7
      ab(4)= - ab(6) - 2*ab(4) - 32*l10 + 25*l8 + 21*ab(1)
      ab(4)=9*l6 + 2*ab(4)
      ab(4)=ab(4)*l6
      ab(7)=2*l5
      ab(8)=ab(7) - l11
      ab(9)=2*ab(1)
      ab(6)= - ab(9) - 12 - 22*l8 + ab(6) - ab(8)
      ab(10)=2*l8
      ab(6)=ab(6)*ab(10)
      ab(11)=4*ab(1)
      ab(12)=ab(11) - l11
      ab(13)= - l5 + ab(12) + 8
      ab(14)=9*l7
      ab(13)=32*l8 - 15*l3 - ab(14) - 11*l10 + 2*ab(13)
      ab(15)=2*l10
      ab(13)=ab(13)*ab(15)
      ab(16)=3*l2
      ab(17)=ab(16) - l11
      ab(18)=4*l5
      ab(19)=ab(18) - 4
      ab(20)=ab(19) - ab(17)
      ab(21)=6*l2
      ab(20)=ab(20)*ab(21)
      ab(22)=3*l4
      ab(23)=ab(22) + ab(21)
      ab(24)=ab(23) - l11
      ab(25)=ab(24) - ab(19)
      ab(26)=6*l4
      ab(25)=ab(25)*ab(26)
      ab(27)=ab(26) + ab(21)
      ab(28)=ab(27) - l7
      ab(29)=ab(28) - l11
      ab(30)=ab(29) - ab(19)
      ab(31)=2*l7
      ab(30)=ab(30)*ab(31)
      ab(32)=ab(2) - ab(7)
      ab(32)=ab(32)*ab(18)
      ab(33)=d10 - d9
      ab(34)= - 6*ab(33)
      ab(35)=Pi**2
      ab(36)= - ab(34) + 7*ab(35)
      ab(37)=d16 + d13
      ab(38)=d15 - d14
      ab(39)=d8 + d7
      ab(40)=l11 - 4
      ab(40)=ab(40)*l11
      ab(41)=d12 - d6
      ab(42)= - 32*ab(41)
      ab(4)= - ab(4) + ab(13) - ab(40) + ab(20) - ab(25) + ab(30) + 
     & ab(42) + ab(36) + 11*d11 + ab(6) - ab(32) - 22*ab(39) - 38*
     & ab(38) - 6*ab(37)
      ab(6)=l6 + l8
      ab(13)=3*l5
      ab(6)= - 8*l10 - ab(13) + ab(27) - ab(31) + 6*ab(6)
      ab(20)= - ab(5) - 2 - l11 + ab(6)
      ab(25)=xa**2
      ab(20)=ab(20)*ab(25)
      ab(30)=9*l3
      ab(32)=11*l2
      ab(40)=ab(30) + 14*l4 + ab(32) - l7
      ab(43)=5*l5
      ab(44)= - 7 + ab(43)
      ab(44)=2*ab(44) - ab(40)
      ab(45)=2*l3
      ab(44)=ab(44)*ab(45)
      ab(46)=d5 + d4
      ab(47)=16*ab(46)
      ab(48)=ab(47) + 24
      ab(20)=4*ab(20) + ab(44) - ab(48) - ab(4)
      ab(20)=ab(20)*ab(25)
      ab(44)= - 3 - ab(43)
      ab(40)=2*ab(44) + ab(40)
      ab(40)=ab(40)*ab(45)
      ab(4)=ab(20) + ab(40) + 8 + ab(47) + ab(4)
      ab(4)=ab(4)*ab(25)
      ab(20)= - l7 - l3 + 14*l5
      ab(40)=11*l11
      ab(44)=4*l8
      ab(49)=ab(20) - ab(40) + ab(44) - ab(9)
      ab(49)=ab(49)*ab(10)
      ab(50)=9*l11
      ab(51)=ab(11) - 9*l5 + ab(50) + 2
      ab(52)=l10 + l7
      ab(44)= - ab(44) - ab(52) - 11*l3 + 2*ab(51)
      ab(44)=ab(44)*ab(15)
      ab(22)=ab(16) + ab(22)
      ab(20)= - 4*l10 + 3*l8 + ab(22) + ab(20)
      ab(20)= - 13*l6 + 2*ab(20)
      ab(20)=ab(20)*l6
      ab(50)=ab(50) - ab(7)
      ab(50)=ab(50)*ab(18)
      ab(19)= - ab(19) + 7*l11
      ab(51)=ab(19) + ab(28)
      ab(51)=ab(51)*ab(31)
      ab(53)=ab(19) + ab(23)
      ab(53)=ab(53)*ab(26)
      ab(19)=ab(19) + ab(16)
      ab(19)=ab(19)*ab(21)
      ab(40)=ab(40) - 8
      ab(40)=ab(40)*l11
      ab(19)= - ab(34) - ab(42) + ab(40) - ab(50) + ab(51) - ab(53) + 
     & ab(44) + ab(20) + 64*d17 + 3*d11 + ab(49) - ab(19) - 6*ab(39) + 
     & 10*ab(38) - 22*ab(37)
      ab(20)=3*l7
      ab(34)=9*l8
      ab(40)=l3 + l11
      ab(9)= - ab(40) + ab(34) + ab(9) - ab(20)
      ab(42)=ab(9) + 6
      ab(42)=ab(42)*ab(10)
      ab(37)=2*ab(37)
      ab(41)=4*d17 + ab(37) + 14*ab(38) + 12*ab(41)
      ab(44)=l11 + 4
      ab(49)=ab(44)*l11
      ab(49)=ab(49) + ab(41)
      ab(1)= - l3 - ab(20) - 14*l10 + 11*l8 + 9*ab(1)
      ab(20)=ab(3) - ab(1)
      ab(50)=5*l6
      ab(20)= - ab(50) + 2*ab(20)
      ab(20)=ab(20)*l6
      ab(20)=ab(20) - ab(42) - ab(49) + 16*l5
      ab(29)=ab(29)*ab(31)
      ab(24)=ab(24)*ab(26)
      ab(17)=ab(17)*ab(21)
      ab(17)= - ab(17) + ab(29) - ab(24)
      ab(24)= - ab(36) - 13*d11 + 26*ab(39) + 32*ab(46)
      ab(29)=ab(17) - ab(24)
      ab(36)=12*l4
      ab(32)=ab(32) + ab(36)
      ab(42)=4*l3
      ab(51)=ab(42) + ab(32) - ab(14)
      ab(53)= - 14 + ab(51)
      ab(53)=ab(53)*ab(45)
      ab(54)=ab(2) - ab(11)
      ab(54)=2*ab(54)
      ab(55)=23*l3
      ab(56)=6*l10 - 14*l8 + 5*l7
      ab(57)=ab(55) + ab(56)
      ab(58)= - ab(54) - ab(57)
      ab(58)=ab(58)*ab(15)
      ab(53)=ab(58) + ab(53) + ab(29) + ab(20)
      ab(53)=ab(53)*ab(25)
      ab(58)=5*l3
      ab(14)= - 46*l4 - 35*l2 + ab(58) + ab(14)
      ab(3)= - ab(3) + ab(13)
      ab(3)=2*ab(3) - ab(14)
      ab(3)=ab(3)*ab(45)
      ab(3)=ab(53) + ab(3) - ab(35) + ab(48) + ab(19)
      ab(3)=ab(3)*ab(25)
      ab(48)=l3 - l5
      ab(53)=ab(48) + ab(31)
      ab(27)=l6 + ab(27) - ab(53) - 9*l10 + 7*l8
      ab(2)=ab(27) - ab(2)
      ab(59)=2*l6
      ab(2)=ab(2)*ab(59)
      ab(60)=ab(53) - 3
      ab(61)=13*l8
      ab(60)= - ab(61) + 2*ab(60)
      ab(60)=ab(60)*l8
      ab(62)=ab(8) - 1
      ab(62)=ab(62)*ab(7)
      ab(36)=ab(36) + 12*l2
      ab(63)=4*l7
      ab(64)= - ab(63) + ab(36)
      ab(65)=l5 - 1
      ab(65)=ab(64)*ab(65)
      ab(37)=10*d12 + ab(37) + 12*ab(38)
      ab(2)=ab(60) - ab(2) + ab(62) + 4*l11 + ab(65) - ab(37)
      ab(38)=d17 + 5*d6
      ab(60)= - ab(35) - ab(38) + 5*ab(39) - ab(33)
      ab(62)= - 6 + ab(60)
      ab(43)=ab(43) - l4
      ab(65)= - ab(30) + 2*ab(43)
      ab(65)=ab(65)*l3
      ab(66)=5*d11
      ab(65)=ab(65) - ab(66)
      ab(34)= - l5 + ab(34) - ab(31)
      ab(67)=ab(34) + 5
      ab(68)=6*l3
      ab(69)=ab(68) - ab(67)
      ab(70)=5*l10
      ab(69)=2*ab(69) + ab(70)
      ab(69)=l10*ab(69)
      ab(62)=ab(69) + 2*ab(62) + ab(65) - ab(2)
      ab(3)=2*ab(62) + ab(3)
      ab(3)=xa*ab(3)
      ab(9)=ab(9) - 6
      ab(9)=ab(9)*ab(10)
      ab(10)=l11 + 3
      ab(1)= - ab(1) + 2*ab(10)
      ab(1)= - ab(50) + 2*ab(1)
      ab(1)=ab(1)*l6
      ab(1)=ab(1) - ab(9)
      ab(9)=ab(28) - ab(44)
      ab(9)=ab(9)*ab(31)
      ab(23)=ab(23) - ab(44)
      ab(23)=ab(23)*ab(26)
      ab(16)=ab(16) - ab(44)
      ab(16)=ab(16)*ab(21)
      ab(9)= - ab(9) + ab(23) + ab(16) - ab(1) + ab(49)
      ab(16)=ab(51) - 18
      ab(16)=ab(16)*ab(45)
      ab(12)=ab(12) - 5
      ab(12)= - ab(56) + 2*ab(12)
      ab(21)=ab(12) - ab(55)
      ab(21)=ab(21)*ab(15)
      ab(16)= - ab(21) - ab(16) + ab(24) + ab(9) + 8
      ab(16)=ab(16)*ab(25)
      ab(21)=l11 + 8
      ab(21)=ab(21)*l11
      ab(1)= - ab(21) - ab(18) + ab(1) - ab(41)
      ab(21)=16 - ab(51)
      ab(21)=ab(21)*ab(45)
      ab(11)=ab(11) - ab(10)
      ab(11)=2*ab(11)
      ab(23)= - ab(11) + ab(57)
      ab(23)=ab(23)*ab(15)
      ab(21)=ab(16) + ab(23) + ab(21) - ab(29) - ab(1)
      ab(21)=ab(21)*ab(25)
      ab(23)=ab(27) - ab(10)
      ab(23)=ab(23)*ab(59)
      ab(24)=ab(53) + 3
      ab(24)= - ab(61) + 2*ab(24)
      ab(24)=ab(24)*l8
      ab(23)=ab(24) - ab(23) - ab(37)
      ab(8)=ab(8) - 2
      ab(8)=ab(8)*ab(7)
      ab(24)=l5 + 1
      ab(24)=ab(64)*ab(24)
      ab(8)=ab(8) + 2*l11 + ab(24) + ab(23)
      ab(24)= - 2 - ab(60)
      ab(26)= - 3 - ab(43)
      ab(26)=2*ab(26) + ab(30)
      ab(26)=l3*ab(26)
      ab(27)=ab(34) - 5
      ab(28)= - ab(68) + ab(27)
      ab(28)=2*ab(28) - ab(70)
      ab(28)=l10*ab(28)
      ab(24)=ab(28) + ab(26) + 2*ab(24) + ab(66) + ab(8)
      ab(21)=2*ab(24) + ab(21)
      ab(21)=ab(21)*ab(25)
      ab(22)=ab(52) - ab(22)
      ab(24)=1 + ab(48) + ab(22)
      ab(16)=8*ab(24) - ab(16)
      ab(16)=ab(16)*ab(25)
      ab(10)=ab(10) - ab(7)
      ab(7)=ab(10)*ab(7)
      ab(10)= - l5*ab(36)
      ab(18)=ab(18)*l7
      ab(7)=ab(7) - ab(23) + ab(18) + ab(10)
      ab(10)=ab(34) - 3
      ab(18)=ab(10) - ab(68)
      ab(18)= - ab(70) + 2*ab(18)
      ab(18)=ab(18)*l10
      ab(18)=ab(18) - ab(7) - ab(65) - 2*ab(60)
      ab(18)=2*ab(18)
      ab(16)= - ab(18) + ab(16)
      ab(16)=xb*xa*ab(16)
      ab(16)=ab(16) + ab(18) + ab(21)
      ab(16)=xb*ab(16)
      ab(3)=ab(3) + ab(16)
      ab(3)=xb*ab(3)
      ab(6)= - ab(6) + ab(40) + 6
      ab(3)=ab(3) + 4*ab(6) + ab(4)
      ab(3)=xb*ab(3)
      ab(4)=ab(38) + ab(39)
      ab(6)=1 + 6*ab(46)
      ab(6)=2*ab(6) + ab(4)
      ab(16)=2 - l5
      ab(16)=ab(63) + 5*ab(16) + l4
      ab(16)=2*ab(16) + ab(58)
      ab(16)=l3*ab(16)
      ab(18)=ab(42) + ab(67)
      ab(18)=2*ab(18) - ab(70)
      ab(18)=l10*ab(18)
      ab(2)=ab(18) + ab(16) + 2*ab(6) - d11 + ab(2)
      ab(6)=2*ab(25)
      ab(2)=ab(2)*ab(6)
      ab(16)= - 1 + l11
      ab(13)=2*ab(16) - ab(13)
      ab(13)=2*ab(13) + ab(14)
      ab(13)=ab(13)*ab(45)
      ab(14)=ab(35) - 8
      ab(2)=ab(2) + ab(13) + ab(14) - ab(47) - ab(19)
      ab(2)=ab(2)*ab(25)
      ab(13)=ab(33) - ab(39)
      ab(13)= - d11 - ab(47) - 2*ab(13)
      ab(16)= - ab(13) + ab(17) + 3*ab(35)
      ab(17)=ab(56) + ab(5)
      ab(18)=ab(54) + ab(17)
      ab(18)=ab(18)*ab(15)
      ab(19)=ab(32) - l7
      ab(21)= - 2 - ab(19)
      ab(21)=ab(21)*ab(45)
      ab(2)=ab(2) + ab(18) + ab(21) - ab(16) - ab(20)
      ab(2)=xa*ab(2)
      ab(2)=ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(3)=ab(19) - 2
      ab(3)=ab(3)*ab(45)
      ab(12)=ab(12) - ab(5)
      ab(12)=ab(12)*ab(15)
      ab(3)=ab(12) + ab(3) - ab(13) - ab(9) + 3*ab(14)
      ab(9)=ab(43) - ab(63)
      ab(12)=ab(9) - 6
      ab(12)= - ab(58) + 2*ab(12)
      ab(12)=ab(12)*l3
      ab(10)=ab(10) + ab(42)
      ab(10)= - ab(70) + 2*ab(10)
      ab(10)=ab(10)*l10
      ab(13)=ab(4) + 12*ab(46)
      ab(7)= - ab(12) + 2*ab(13) + ab(10) - ab(7) - d11
      ab(7)=ab(7)*ab(25)
      ab(10)=1 - 2*ab(46)
      ab(4)=6*ab(10) - ab(4)
      ab(9)= - 1 + ab(9)
      ab(9)=2*ab(9) - ab(58)
      ab(9)=l3*ab(9)
      ab(10)= - ab(42) - ab(27)
      ab(10)=2*ab(10) + ab(70)
      ab(10)=l10*ab(10)
      ab(4)= - ab(7) + ab(10) + ab(9) + 2*ab(4) + d11 - ab(8)
      ab(4)=ab(4)*ab(6)
      ab(8)=ab(11) - ab(17)
      ab(8)=ab(8)*ab(15)
      ab(9)=12 + ab(19)
      ab(9)=ab(9)*ab(45)
      ab(1)=ab(4) + ab(8) + ab(9) + ab(16) + ab(1)
      ab(1)=ab(1)*ab(25)
      ab(1)=ab(2) + ab(1) + ab(3)
      ab(1)=xb*ab(1)
      ab(2)= - ab(5) - 3 + l5 - ab(22)
      ab(2)=4*ab(2) + ab(7)
      ab(2)=ab(2)*ab(6)
      ab(2)=ab(2) - ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

