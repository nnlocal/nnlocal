      subroutine BGTAa12carsssggg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      complex(ki) ab(76)
      include 'BGTAa12carsss_functions.h'

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(1,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      

      tmp =  - 8*CA**2*l1*z6*z9*xa**3
      res(1,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=Pi**2
      ab(2)=l4 + l2 - l3
      ab(2)=8*ab(2) - 5*l1
      ab(2)=l1*ab(2)
      ab(1)=1.0_ki/3.0_ki*ab(1) + ab(2)

      tmp = 2*CA**2*z6*z9*xa**3*ab(1)
      res(1,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l2**2
      ab(2)=2*l2 - l3
      ab(2)=l3*ab(2)
      ab(3)=l3 - l2
      ab(4)=2*ab(3) - l4
      ab(4)=l4*ab(4)
      ab(1)= - ab(4) + ab(1) - ab(2)
      ab(2)=Pi**2
      ab(3)=ab(3) - l4
      ab(4)= - 5*ab(3) - 2*l1
      ab(4)=l1*ab(4)
      ab(1)=4*ab(4) + 1.0_ki/3.0_ki*ab(2) - 16*ab(1)
      ab(1)=l1*ab(1)
      ab(3)=4.0_ki/3.0_ki*ab(3)
      ab(2)=ab(2)*ab(3)
      ab(1)=ab(1) - 17.0_ki/2.0_ki*zeta3 + ab(2)

      tmp = CA**2*z6*z9*xa**3*ab(1)
      res(1,0,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - xa + xb

      tmp = 6*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xb - xa
      ab(2)=l5 + l7 - l6
      ab(2)=2*l1 - 3*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(1)=3 + ab(1)

      tmp = 4*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=3*l5
      ab(2)= - ab(1) + 2*l1
      ab(3)=3*l6
      ab(4)=ab(2) + ab(3)
      ab(5)=3*l7
      ab(4)= - ab(5) + 2*ab(4)
      ab(6)=2*l7
      ab(4)=ab(4)*ab(6)
      ab(3)=ab(3) + 2*ab(2)
      ab(6)=2*l6
      ab(3)=ab(3)*ab(6)
      ab(1)= - ab(1) + 4*l1
      ab(6)=2*l5
      ab(1)=ab(1)*ab(6)
      ab(6)=Pi**2
      ab(7)=l1**2
      ab(1)=ab(4) - ab(3) + ab(1) - 2.0_ki/3.0_ki*ab(6) - 5*ab(7)
      ab(3)= - xb + xa
      ab(1)=ab(1)*ab(3)
      ab(2)= - ab(5) + 6*l6 - 3*l8 + ab(2)
      ab(1)=4*ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=3*xa
      ab(1)=ab(1)*l6
      ab(2)=l1*xa
      ab(3)=2*ab(2)
      ab(4)=ab(3) - 3
      ab(5)=3*l5
      ab(6)=ab(5)*xa
      ab(7)=l7*xa
      ab(7)=ab(7) - ab(1) + ab(6) - ab(4)
      ab(8)=2*l7
      ab(7)=ab(7)*ab(8)
      ab(9)=Pi**2
      ab(10)=ab(9)*xa
      ab(11)=1.0_ki/3.0_ki*ab(10)
      ab(12)=ab(11) + 6*l8
      ab(13)=5*ab(2)
      ab(14)=ab(13) - 8
      ab(14)=ab(14)*l1
      ab(12)=ab(14) + 2*ab(12)
      ab(14)= - 2*ab(4) + ab(6)
      ab(15)=2*l5
      ab(14)=ab(14)*ab(15)
      ab(16)=ab(2) - 3
      ab(17)=2*ab(16) - ab(6)
      ab(1)=2*ab(17) + ab(1)
      ab(17)=2*l6
      ab(1)=ab(1)*ab(17)
      ab(1)=ab(7) + ab(1) + ab(14) + ab(12)
      ab(7)=4*l7
      ab(1)=ab(1)*ab(7)
      ab(14)=4*l1
      ab(18)=ab(14) - ab(5)
      ab(18)=ab(18)*ab(15)
      ab(19)=l1**2
      ab(20)=5*ab(19) + 2.0_ki/3.0_ki*ab(9)
      ab(18)=ab(18) - ab(20)
      ab(21)=2*l1
      ab(5)=ab(21) - ab(5)
      ab(22)=3*l6
      ab(23)= - 2*ab(5) - ab(22)
      ab(23)=ab(23)*ab(17)
      ab(22)= - l7 + ab(22) + ab(5)
      ab(8)=ab(22)*ab(8)
      ab(8)=ab(8) + ab(23) + ab(18)
      ab(7)=ab(8)*ab(7)
      ab(8)=ab(21) - l5
      ab(8)=ab(8)*ab(15)
      ab(8)=ab(8) - ab(20)
      ab(20)=4*l5
      ab(8)=ab(8)*ab(20)
      ab(5)=l6 + ab(5)
      ab(5)=ab(5)*ab(17)
      ab(5)=ab(5) - ab(18)
      ab(18)=4*l6
      ab(5)=ab(5)*ab(18)
      ab(19)=5.0_ki/3.0_ki*ab(9) + 8*ab(19)
      ab(19)=l1*ab(19)
      ab(21)=17.0_ki/2.0_ki*zeta3
      ab(5)=ab(7) + ab(5) + ab(8) - ab(21) + ab(19)
      ab(5)=xb*ab(5)
      ab(2)=6 - ab(2)
      ab(7)= - l6*xa
      ab(2)=ab(7) + 2*ab(2) + ab(6)
      ab(2)=ab(2)*ab(17)
      ab(6)=4*ab(16) - ab(6)
      ab(6)=ab(6)*ab(15)
      ab(7)= - 12*l8 - ab(11)
      ab(8)=16 - ab(13)
      ab(8)=l1*ab(8)
      ab(2)=ab(2) + ab(6) + 2*ab(7) + ab(8)
      ab(2)=ab(2)*ab(18)
      ab(6)=l5*xa
      ab(4)=ab(6) - ab(4)
      ab(4)=ab(4)*ab(15)
      ab(4)=ab(4) + ab(12)
      ab(4)=ab(4)*ab(20)
      ab(6)=l8**2
      ab(6)=3*ab(6) + 1.0_ki/3.0_ki*ab(9)
      ab(3)=5 - ab(3)
      ab(3)=ab(3)*ab(14)
      ab(3)=ab(3) - 32*l8 - 5.0_ki/3.0_ki*ab(10)
      ab(3)=l1*ab(3)
      ab(7)=xa*ab(21)
      ab(1)=ab(5) + ab(1) + ab(2) + ab(4) + ab(3) + 8*ab(6) + ab(7)

      tmp = CA**2*z1*z2*z4*xb**3*ab(1)
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
      ab(8)=l6 + l5
      ab(7)=ab(7) + 6*ab(8)
      ab(9)=3*l1
      ab(10)=2 - ab(9)
      ab(10)=l1*ab(10)
      ab(10)=ab(10) + ab(6) + ab(7)
      ab(10)=xa*ab(10)
      ab(8)=3*ab(8)
      ab(11)=1 - 2*l2 - ab(8)
      ab(12)= - 16 + ab(4)
      ab(12)=l1*ab(12)
      ab(5)=ab(10) + ab(5) + ab(12) + 6*ab(11) - ab(2)
      ab(5)=xa*ab(5)
      ab(10)= - 2 + 3*l2
      ab(10)=ab(8) + 2*ab(10)
      ab(11)=l1 - 4
      ab(12)=6*l1
      ab(11)=ab(11)*ab(12)
      ab(1)= - ab(11) + 7.0_ki/3.0_ki*ab(1) + 4*ab(10)
      ab(10)= - 2 + ab(4)
      ab(10)=l1*ab(10)
      ab(7)=ab(10) - ab(2) - ab(7)
      ab(7)=xb*ab(7)
      ab(10)= - 4*l1 + 1 - 6*l2 + ab(8)
      ab(7)=4*ab(10) + ab(7)
      ab(7)=xb*ab(7)
      ab(5)=ab(5) + ab(7) + ab(1)
      ab(5)=xa*ab(5)
      ab(7)= - 10 + ab(9)
      ab(7)=l1*ab(7)
      ab(3)=ab(7) - ab(3) - ab(6)
      ab(3)=xb*ab(3)
      ab(6)= - 5 + 18*l2 + ab(8)
      ab(4)=32 - ab(4)
      ab(4)=l1*ab(4)
      ab(2)=ab(3) + ab(4) + 2*ab(6) + ab(2)
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
      ab(11)=3*l5
      ab(11)=ab(10)*ab(11)
      ab(12)=3*ab(5)
      ab(13)=ab(3)*xb
      ab(14)=17*ab(13) - 31*ab(3) - ab(12)
      ab(14)=xb*ab(14)
      ab(15)=2*l1
      ab(16)=ab(10)*ab(15)
      ab(17)=5*ab(3)
      ab(18)=ab(17) - 2*ab(5)
      ab(18)=xa*ab(18)
      ab(19)=8*ab(3)
      ab(14)=ab(11) + ab(16) + ab(14) + ab(19) + 3*ab(18)
      ab(14)=l5*ab(14)
      ab(16)=l1*ab(10)
      ab(11)=ab(16) + ab(11)
      ab(10)=l6*ab(10)
      ab(10)=2*ab(11) + 3*ab(10)
      ab(10)=l6*ab(10)
      ab(11)=ab(13) - ab(4)
      ab(11)=ab(11)*xb
      ab(16)=ab(1) - ab(5)
      ab(16)=ab(16)*xa
      ab(11)=ab(11) + ab(16) + ab(1)
      ab(11)=ab(11)*xb
      ab(11)=ab(11) - ab(16)
      ab(18)=12*l2 + 16*l1
      ab(18)=ab(11)*ab(18)
      ab(20)=4*ab(3)
      ab(21)=ab(20) - ab(5)
      ab(21)=xa*ab(21)
      ab(22)= - ab(1) - ab(5)
      ab(4)=ab(4)*xb
      ab(22)=2*ab(22) + ab(4)
      ab(22)=xb*ab(22)
      ab(18)=ab(21) + ab(22) + ab(18)
      ab(18)=l2*ab(18)
      ab(21)=9*ab(3) - ab(5)
      ab(21)=xa*ab(21)
      ab(21)= - 7*ab(9) - ab(19) + ab(21)
      ab(21)=l9*ab(21)
      ab(10)=ab(14) + ab(10) + ab(18) + ab(21)
      ab(4)=ab(4) - 10*ab(3) + ab(5)
      ab(4)=xb*ab(4)
      ab(12)=ab(20) - ab(12)
      ab(12)=xa*ab(12)
      ab(4)=ab(4) + ab(19) + ab(12)
      ab(4)=xb*ab(4)
      ab(12)=6*ab(3) - ab(5)
      ab(12)=xa*ab(12)
      ab(12)= - ab(19) + ab(12)
      ab(12)=xa*ab(12)
      ab(4)=ab(12) + ab(4)
      ab(12)=xb - 3
      ab(12)=ab(12)*ab(9)
      ab(14)=xa**2
      ab(6)=ab(6)*ab(14)
      ab(6)=ab(12) + ab(6)
      ab(12)= - ab(1) + ab(6)
      ab(12)=xb*ab(12)
      ab(2)= - ab(2) + ab(12)
      ab(2)=ab(2)*ab(15)
      ab(2)=3*ab(4) + ab(2)
      ab(2)=l1*ab(2)
      ab(4)= - ab(1) + ab(13)
      ab(4)=xb*ab(4)
      ab(4)=5*ab(4) + ab(20) + ab(16)
      ab(11)=l3*ab(11)
      ab(2)=8*ab(11) + 4*ab(4) + ab(2)
      ab(2)=l1*ab(2)
      ab(4)= - xb*ab(1)
      ab(4)=ab(4) + ab(17) + ab(5)
      ab(11)=1.0_ki/3.0_ki*xb
      ab(4)=ab(4)*ab(11)
      ab(5)= - ab(3) + 1.0_ki/3.0_ki*ab(5)
      ab(12)=ab(5)*xa
      ab(13)= - 1.0_ki/3.0_ki*ab(3) + ab(12)
      ab(4)=2*ab(13) + ab(4)
      ab(4)=xb*ab(4)
      ab(13)= - xa*ab(8)
      ab(1)=ab(1) + ab(13)
      ab(1)=xa*ab(1)
      ab(1)=1.0_ki/3.0_ki*ab(1) + ab(4)
      ab(4)=ab(5)*ab(14)
      ab(5)= - ab(11) + 1
      ab(5)=ab(5)*ab(9)
      ab(9)=2.0_ki/3.0_ki*ab(3)
      ab(4)=ab(5) + ab(9) + ab(4)
      ab(4)=xb*ab(4)
      ab(5)= - ab(9) - ab(12)
      ab(5)=xa*ab(5)
      ab(4)=ab(5) + ab(4)
      ab(4)=l1*ab(4)
      ab(1)=2*ab(1) + 7*ab(4)
      ab(1)=ab(1)*Pi**2
      ab(4)= - ab(3) + 1.0_ki/2.0_ki*ab(6)
      ab(4)=xb*ab(4)
      ab(3)=ab(3) - 1.0_ki/2.0_ki*ab(7)
      ab(3)=xa*ab(3)
      ab(3)=ab(3) + ab(4)
      ab(3)=zeta3*ab(3)
      ab(4)= - xa + xb
      ab(4)=l8*ab(8)*ab(4)

      tmp = ab(1) + ab(2) + 97.0_ki/2.0_ki*ab(3) + 8*ab(4) + 2*ab(10)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=xa*ab(3)
      ab(5)=ab(4) - ab(3)
      ab(6)=ab(5)*xb
      ab(7)=9*ab(3)
      ab(8)=ab(7) - ab(4)
      ab(8)=ab(8)*xa
      ab(9)=8*ab(3)
      ab(8)=ab(8) - ab(9) - 7*ab(6)
      ab(10)=2*l1
      ab(11)= - l2 - ab(10)
      ab(11)=ab(8)*ab(11)
      ab(12)=25*ab(6)
      ab(13)=3*ab(4)
      ab(14)= - 31*ab(3) + ab(13)
      ab(14)=xa*ab(14)
      ab(14)=ab(12) + 28*ab(3) + ab(14)
      ab(14)=l5*ab(14)
      ab(15)=5*ab(6)
      ab(16)=7*ab(3)
      ab(17)= - ab(16) + ab(4)
      ab(17)=xa*ab(17)
      ab(18)=6*ab(3)
      ab(17)=ab(15) + ab(18) + ab(17)
      ab(17)=l9*ab(17)
      ab(11)=ab(17) + ab(14) + ab(11)
      ab(11)=l9*ab(11)
      ab(14)=10*ab(3)
      ab(17)=3*ab(3)
      ab(19)=ab(17)*xb
      ab(20)= - ab(19) + ab(14) - ab(4)
      ab(20)=xb*ab(20)
      ab(21)=4*ab(3)
      ab(22)= - ab(21) + ab(13)
      ab(22)=xa*ab(22)
      ab(20)=ab(20) - ab(9) + ab(22)
      ab(20)=xb*ab(20)
      ab(22)=1.0_ki/3.0_ki*xb
      ab(23)=ab(22) - 1
      ab(23)=ab(23)*ab(6)
      ab(24)= - ab(3) + 1.0_ki/3.0_ki*ab(4)
      ab(25)=xa**2
      ab(26)=ab(25)*ab(24)
      ab(23)=ab(23) - ab(26)
      ab(26)=2.0_ki/3.0_ki*ab(3)
      ab(27)= - ab(26) + ab(23)
      ab(27)=ab(27)*xb
      ab(24)=ab(24)*xa
      ab(26)=ab(24) + ab(26)
      ab(26)=ab(26)*xa
      ab(27)=ab(27) + ab(26)
      ab(28)= - ab(27)*ab(10)
      ab(18)= - ab(18) + ab(4)
      ab(18)=xa*ab(18)
      ab(18)=ab(9) + ab(18)
      ab(18)=xa*ab(18)
      ab(18)=ab(28) + ab(18) + ab(20)
      ab(18)=l1*ab(18)
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)= - xb*ab(1)
      ab(2)=ab(2) + ab(16) - ab(13)
      ab(2)=xb*ab(2)
      ab(20)=ab(17) + ab(4)
      ab(28)=3*xa
      ab(20)=ab(20)*ab(28)
      ab(29)=22*ab(3)
      ab(2)=ab(18) + 5*ab(2) - ab(29) + ab(20)
      ab(18)=l1**2
      ab(2)=ab(2)*ab(18)
      ab(20)=ab(3)*xb
      ab(30)=ab(20) - ab(17)
      ab(30)=ab(30)*xb
      ab(31)=ab(1) - ab(4)
      ab(31)=ab(31)*xa
      ab(32)=ab(31) + ab(1)
      ab(30)=ab(30) + ab(32)
      ab(30)=ab(30)*xb
      ab(30)=ab(30) - ab(31)
      ab(18)=3*ab(18)
      ab(33)=8*l1
      ab(34)= - l2*ab(33)
      ab(35)=4*l1
      ab(36)= - l3*ab(35)
      ab(34)=ab(36) + ab(34) - ab(18)
      ab(34)=ab(30)*ab(34)
      ab(36)=ab(3) - 1.0_ki/3.0_ki*ab(20)
      ab(36)=xb*ab(36)
      ab(32)= - 1.0_ki/3.0_ki*ab(32) + ab(36)
      ab(32)=xb*ab(32)
      ab(31)=1.0_ki/3.0_ki*ab(31) + ab(32)
      ab(32)=Pi**2
      ab(31)=ab(31)*ab(32)
      ab(31)=ab(31) + ab(34)
      ab(31)=l3*ab(31)
      ab(34)=23*ab(3)
      ab(36)= - ab(34) - ab(4)
      ab(36)=xa*ab(36)
      ab(37)=24*ab(3)
      ab(12)=ab(12) + ab(37) + ab(36)
      ab(12)=d1*ab(12)
      ab(2)=ab(12) + ab(11) + ab(2) + ab(31)
      ab(11)=1.0_ki/2.0_ki*ab(4)
      ab(12)=17.0_ki/6.0_ki*ab(20) - ab(9) - ab(11)
      ab(12)=xb*ab(12)
      ab(31)=20*ab(3)
      ab(36)=ab(31) - 17.0_ki/2.0_ki*ab(4)
      ab(36)=xa*ab(36)
      ab(36)=14*ab(3) + ab(36)
      ab(12)=1.0_ki/3.0_ki*ab(36) + ab(12)
      ab(12)=xb*ab(12)
      ab(27)=ab(27)*ab(33)
      ab(11)=4.0_ki/3.0_ki*ab(3) + ab(11)
      ab(11)=xa*ab(11)
      ab(11)= - 14.0_ki/3.0_ki*ab(3) + ab(11)
      ab(11)=xa*ab(11)
      ab(11)=ab(27) + ab(11) + ab(12)
      ab(11)=l1*ab(11)
      ab(12)=17.0_ki/3.0_ki*ab(20) - 25.0_ki/3.0_ki*ab(3) - ab(13)
      ab(12)=xb*ab(12)
      ab(27)= - ab(6) + 2*ab(5)
      ab(22)=ab(27)*ab(22)
      ab(22)=ab(22) + ab(26)
      ab(26)=4*l5
      ab(36)=ab(22)*ab(26)
      ab(38)=1.0_ki/3.0_ki*ab(3)
      ab(23)= - ab(38) + 1.0_ki/2.0_ki*ab(23)
      ab(23)=xb*ab(23)
      ab(24)=ab(38) + 1.0_ki/2.0_ki*ab(24)
      ab(24)=xa*ab(24)
      ab(23)=ab(24) + ab(23)
      ab(23)=ab(23)*ab(32)
      ab(24)=4*ab(4)
      ab(38)=ab(3) + ab(24)
      ab(38)=xa*ab(38)
      ab(11)=77.0_ki/24.0_ki*ab(23) + ab(36) + ab(11) + ab(12) + ab(21) 
     & + 1.0_ki/3.0_ki*ab(38)
      ab(11)=ab(11)*ab(32)
      ab(12)=ab(5)*xa
      ab(12)=ab(12) - ab(6)
      ab(23)=l8 + l2
      ab(23)=2*ab(23)
      ab(12)=ab(12)*ab(23)
      ab(23)=ab(6) - ab(5)
      ab(33)=ab(23)*ab(33)
      ab(36)= - 11*ab(3) + ab(13)
      ab(36)=xa*ab(36)
      ab(9)=ab(15) + ab(9) + ab(36)
      ab(9)=l5*ab(9)
      ab(8)=l9*ab(8)
      ab(8)=ab(8) + ab(33) + ab(9) + ab(12)
      ab(8)=l8*ab(8)
      ab(9)= - 17*ab(20) + 35*ab(3) - ab(4)
      ab(9)=xb*ab(9)
      ab(12)= - ab(17) + 2*ab(4)
      ab(12)=xa*ab(12)
      ab(15)=16*ab(3)
      ab(9)=ab(9) - ab(15) + ab(12)
      ab(12)=ab(27)*xb
      ab(27)=ab(17) - ab(4)
      ab(33)=ab(27)*xa
      ab(36)=ab(33) - ab(1)
      ab(36)=ab(36)*xa
      ab(12)=ab(36) - ab(12)
      ab(38)=l1*ab(12)
      ab(9)=2*ab(9) - 3*ab(38)
      ab(9)=ab(9)*ab(10)
      ab(26)= - ab(26) - ab(35)
      ab(26)=ab(12)*ab(26)
      ab(35)= - 59*ab(20) + 115*ab(3) + ab(13)
      ab(35)=xb*ab(35)
      ab(7)= - ab(7) + ab(24)
      ab(7)=ab(7)*ab(28)
      ab(7)=ab(35) - 44*ab(3) + ab(7) + ab(26)
      ab(7)=l5*ab(7)
      ab(7)=ab(9) + ab(7)
      ab(7)=l5*ab(7)
      ab(9)=12*ab(3) + ab(4)
      ab(9)=xa*ab(9)
      ab(14)=ab(14) - 7*ab(4)
      ab(14)=2*ab(14) - ab(19)
      ab(14)=xb*ab(14)
      ab(9)=ab(14) - ab(15) + ab(9)
      ab(14)=l1*ab(30)
      ab(9)=2*ab(9) - 15*ab(14)
      ab(9)=l1*ab(9)
      ab(13)=ab(31) + ab(13)
      ab(13)=xa*ab(13)
      ab(14)=ab(29) - 13*ab(4)
      ab(14)=2*ab(14) - 9*ab(20)
      ab(14)=xb*ab(14)
      ab(13)=ab(14) - 32*ab(3) + ab(13)
      ab(13)=l5*ab(13)
      ab(14)=ab(30)*ab(32)
      ab(9)=3*ab(14) + ab(9) + ab(13)
      ab(13)= - 32*l2 - 56*l1
      ab(13)=ab(30)*ab(13)
      ab(14)= - ab(19) - 25*ab(3) + 31*ab(4)
      ab(14)=xb*ab(14)
      ab(15)= - ab(34) - ab(24)
      ab(15)=xa*ab(15)
      ab(13)=ab(14) + ab(37) + ab(15) + ab(13)
      ab(13)=l2*ab(13)
      ab(9)=2*ab(9) + ab(13)
      ab(9)=l2*ab(9)
      ab(13)=3*l5
      ab(10)= - ab(10) - ab(13)
      ab(10)=l5*ab(12)*ab(10)
      ab(14)=ab(22)*ab(32)
      ab(10)=ab(10) + ab(14)
      ab(14)=2*l6
      ab(13)= - l6 - l1 - ab(13)
      ab(13)=ab(14)*ab(13)
      ab(13)=ab(13) - ab(18)
      ab(12)=ab(12)*ab(13)
      ab(10)=2*ab(10) + ab(12)
      ab(10)=ab(10)*ab(14)
      ab(12)=1.0_ki/4.0_ki*ab(4)
      ab(13)= - 9.0_ki/4.0_ki*ab(20) + ab(16) - ab(12)
      ab(13)=xb*ab(13)
      ab(4)= - ab(21) + 9.0_ki/4.0_ki*ab(4)
      ab(4)=xa*ab(4)
      ab(14)=5*ab(3)
      ab(4)=ab(13) - ab(14) + ab(4)
      ab(4)=xb*ab(4)
      ab(12)= - ab(17) + ab(12)
      ab(12)=xa*ab(12)
      ab(12)=ab(14) + ab(12)
      ab(12)=xa*ab(12)
      ab(4)=ab(12) + ab(4)
      ab(5)= - ab(6) + 3*ab(5)
      ab(5)=ab(5)*xb
      ab(6)=ab(25)*ab(27)
      ab(5)=ab(5) - ab(6)
      ab(6)=ab(3) + 1.0_ki/2.0_ki*ab(5)
      ab(6)=xb*ab(6)
      ab(3)= - ab(3) + 1.0_ki/2.0_ki*ab(33)
      ab(3)=xa*ab(3)
      ab(3)=ab(3) + ab(6)
      ab(3)=l1*ab(3)
      ab(3)=17*ab(4) + 111*ab(3)
      ab(3)=zeta3*ab(3)
      ab(4)=d2 - d3
      ab(4)= - 64*ab(4)
      ab(4)=ab(23)*ab(4)
      ab(1)=ab(1) + ab(5)
      ab(1)=xb*ab(1)
      ab(1)=ab(36) + ab(1)
      ab(1)=plg4half*ab(1)

      tmp = 8*ab(1) + 2*ab(2) + ab(3) + ab(4) + ab(7) + 4*ab(8) + ab(9)
     &  + ab(10) + ab(11)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=xb - 1
      ab(2)=2*ab(1)
      ab(3)=xa*ab(1)
      ab(3)= - ab(2) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)=3*ab(1) + ab(3)
      ab(3)=xa*ab(3)
      ab(2)= - ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 3*CA**2*z1*z6*z10*ab(1)
      res(2,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=CA**2*z1*z6
      ab(2)=2*z10
      ab(2)=ab(1)*ab(2)
      ab(3)=ab(1)*z10
      ab(4)=xa*ab(3)
      ab(4)=ab(4) - ab(2)
      ab(4)=ab(4)*xa
      ab(5)=3*z10
      ab(1)=ab(5)*ab(1)
      ab(1)=ab(4) + ab(1)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) - ab(2)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) + ab(3)
      ab(2)=xb - 1
      ab(3)= - 5*l1 - 6*l2 + l3 + 2*l4
      ab(2)=ab(3)*ab(2)
      ab(2)= - 5 + ab(2)

      tmp = 2*ab(2)*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=CA**2*z1*z6
      ab(2)=2*z10
      ab(2)=ab(1)*ab(2)
      ab(3)=ab(1)*z10
      ab(4)=xa*ab(3)
      ab(4)=ab(4) - ab(2)
      ab(5)=ab(4)*xa
      ab(6)=3*z10
      ab(1)=ab(6)*ab(1)
      ab(1)=ab(5) + ab(1)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) - ab(2)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) + ab(3)
      ab(2)=xb - 1
      ab(5)=ab(2)*ab(1)
      ab(6)=ab(5)*l2
      ab(7)=ab(5)*l3
      ab(8)=ab(5)*l4
      ab(9)=3*ab(1)
      ab(10)=12*ab(6) - 4*ab(7) + ab(9) - 8*ab(8)
      ab(10)=l2*ab(10)
      ab(6)=20*ab(6) - 3*ab(7) + 17*ab(1) - 7*ab(8)
      ab(5)=l1*ab(5)
      ab(5)=2*ab(6) + 17*ab(5)
      ab(5)=l1*ab(5)
      ab(6)= - 14*ab(1) + 3*ab(8)
      ab(6)=l4*ab(6)
      ab(8)= - ab(9) + ab(8)
      ab(7)=2*ab(8) + ab(7)
      ab(7)=l3*ab(7)
      ab(1)=l5*ab(1)
      ab(8)=1.0_ki/3.0_ki*xa
      ab(4)=ab(4)*ab(8)
      ab(4)=ab(4) + ab(3)
      ab(4)=ab(4)*xa
      ab(4)=ab(4) - 2.0_ki/3.0_ki*ab(3)
      ab(4)=ab(4)*xa
      ab(3)=ab(4) + 1.0_ki/3.0_ki*ab(3)
      ab(2)= - Pi**2*ab(3)*ab(2)

      tmp = 34*ab(1) + 5*ab(2) + ab(5) + ab(6) + ab(7) + 2*ab(10)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=12*l2
      ab(2)=xb - 1
      ab(1)=ab(1)*ab(2)
      ab(3)=ab(2)*l1
      ab(4)=ab(3) + 1
      ab(5)=ab(2)*l3
      ab(1)= - ab(5) + ab(1) + 11*ab(4)
      ab(6)=ab(2)*l4
      ab(7)= - ab(1) + 5.0_ki/3.0_ki*ab(6)
      ab(7)=ab(7)*l4
      ab(8)=Pi**2
      ab(9)=ab(8)*xb
      ab(10)=ab(9) - ab(8)
      ab(11)=4*l2
      ab(11)=ab(11)*ab(2)
      ab(4)=ab(11) + 3*ab(4)
      ab(11)= - ab(5) + 2*ab(4)
      ab(11)=ab(11)*l3
      ab(12)=16*l2
      ab(12)=ab(12)*ab(2)
      ab(12)=ab(12) + 3 + 28*ab(3)
      ab(13)=2*l2
      ab(14)=ab(12)*ab(13)
      ab(15)=ab(3) + 2
      ab(15)=ab(15)*l1
      ab(7)=ab(7) - ab(11) + 50*l5 + 25*ab(15) + ab(14) - 7.0_ki/3.0_ki*
     & ab(10)
      ab(7)=ab(7)*l4
      ab(14)=9 + 34*ab(3)
      ab(14)=ab(14)*l1
      ab(16)= - 9*l5 - ab(14) + 10.0_ki/3.0_ki*ab(10)
      ab(17)=ab(2)*l2
      ab(17)=32*ab(17) + 3 + 80*ab(3)
      ab(17)=ab(17)*l2
      ab(16)= - ab(17) + 2*ab(16)
      ab(16)=ab(16)*l2
      ab(18)=8*l2
      ab(2)=ab(18)*ab(2)
      ab(18)=1 + 4*ab(3)
      ab(2)=ab(2) + 3*ab(18)
      ab(2)=ab(2)*ab(13)
      ab(2)=ab(2) + 9*ab(15) + 18*l5 - ab(10)
      ab(13)= - ab(4) + 1.0_ki/3.0_ki*ab(5)
      ab(13)=ab(13)*l3
      ab(13)=ab(13) + ab(2)
      ab(13)=ab(13)*l3
      ab(18)=1 + 1.0_ki/3.0_ki*ab(3)
      ab(19)=59*l1
      ab(18)=ab(18)*ab(19)
      ab(20)=17.0_ki/3.0_ki*ab(8)
      ab(18)=ab(18) + 118*l5 + ab(20) - 17.0_ki/3.0_ki*ab(9)
      ab(18)=ab(18)*l1
      ab(21)=l5**2
      ab(22)=34*zeta3
      ab(21)= - ab(22) + 59*ab(21)
      ab(22)=ab(22)*xb
      ab(7)=ab(7) - ab(18) + ab(20) - ab(22) - ab(21) + ab(16) + ab(13)
      ab(13)=2*ab(7)
      ab(16)=xa*ab(7)
      ab(16)= - ab(13) + ab(16)
      ab(16)=xa*ab(16)
      ab(1)= - 3*ab(1) + 5*ab(6)
      ab(1)=l4*ab(1)
      ab(6)=l2*ab(12)
      ab(1)=ab(1) - 3*ab(11) + 6*ab(6) + 75*ab(15) + 150*l5 - 7*ab(10)
      ab(1)=l4*ab(1)
      ab(3)= - 3 - ab(3)
      ab(3)=ab(3)*ab(19)
      ab(6)=17*ab(8)
      ab(3)=ab(3) + 17*ab(9) - 354*l5 - ab(6)
      ab(3)=l1*ab(3)
      ab(8)= - 3*ab(14) - 27*l5 + 10*ab(10)
      ab(8)=2*ab(8) - 3*ab(17)
      ab(8)=l2*ab(8)
      ab(4)= - 3*ab(4) + ab(5)
      ab(4)=l3*ab(4)
      ab(2)=3*ab(2) + ab(4)
      ab(2)=l3*ab(2)
      ab(4)=xb*zeta3
      ab(1)=ab(16) + ab(1) + ab(2) + ab(8) + ab(3) - 102*ab(4) - 3*
     & ab(21) + ab(6)
      ab(1)=xa*ab(1)
      ab(1)= - ab(13) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + ab(7)

      tmp = CA**2*z1*z6*z10*ab(1)
      res(2,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,2,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(2,2,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xa**2
      ab(2)=2*ab(1)
      ab(3)=1 - ab(2)
      ab(3)=ab(3)*ab(1)
      ab(3)=3 + ab(3)
      ab(3)=xa*ab(3)
      ab(4)=ab(1) + 1
      ab(2)=ab(4)*ab(2)
      ab(2)= - 3 + ab(2)
      ab(2)=ab(2)*ab(1)
      ab(5)=xa**4
      ab(6)=3 - 2*ab(5)
      ab(6)=xb*xa*ab(6)
      ab(2)=ab(6) - 3 + ab(2)
      ab(2)=xb*ab(2)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)= - 1 + ab(1)
      ab(3)=ab(3)*ab(1)
      ab(2)=5*ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=3*ab(1)
      ab(6)= - 1 - ab(3)
      ab(6)=ab(6)*ab(1)
      ab(6)=2 + ab(6)
      ab(6)=xa*ab(6)
      ab(2)=ab(6) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=ab(4)*ab(3)
      ab(3)= - 2 + ab(3)
      ab(1)=ab(3)*ab(1)
      ab(1)=ab(2) - 2 + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=2 - 3*ab(5)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l7 + l5
      ab(2)=5*l10 + l3 + l11 - 7*ab(1)
      ab(3)=l4 + l2
      ab(4)=ab(3) - 2
      ab(4)= - ab(2) + 2*ab(4)
      ab(5)=10*l6
      ab(6)=ab(4) - ab(5)
      ab(7)=xa**2
      ab(8)=ab(7) + 1
      ab(9)= - ab(7)*ab(6)*ab(8)
      ab(10)=l3 - l8 - 5*ab(1)
      ab(11)=2*l10
      ab(12)=ab(10) + ab(11) + 2
      ab(13)=5*l6
      ab(14)=ab(13) + ab(12)
      ab(14)=2*ab(14)
      ab(9)= - ab(14) + ab(9)
      ab(9)=ab(9)*ab(7)
      ab(15)=xa**4
      ab(6)=ab(6)*ab(15)
      ab(6)=ab(14) + ab(6)
      ab(6)=xb*xa*ab(6)
      ab(6)=ab(6) - ab(14) + ab(9)
      ab(6)=xb*ab(6)
      ab(9)= - ab(3) + 5*l8
      ab(9)= - l3 + 2*ab(9) - 7*l11 + l10 + 3*ab(1)
      ab(14)=ab(3) + 2
      ab(2)= - ab(2) + 2*ab(14)
      ab(5)= - ab(5) + ab(2)
      ab(5)=ab(5)*ab(7)
      ab(5)=ab(5) - ab(9)
      ab(5)=ab(5)*ab(7)
      ab(10)=ab(10) + ab(11) - 2
      ab(11)=ab(13) + ab(10)
      ab(5)=2*ab(11) + ab(5)
      ab(5)=xa*ab(5)
      ab(5)=ab(5) + ab(6)
      ab(5)=xb*ab(5)
      ab(3)=l8 + ab(3) + 4
      ab(1)= - 2*ab(3) + 3*l3 + 12*l6 + l11 + 9*l10 - 17*ab(1)
      ab(3)= - 8*ab(7) + ab(1)
      ab(3)=ab(3)*ab(7)
      ab(1)=ab(3) - ab(1)
      ab(1)=ab(1)*ab(7)
      ab(1)=ab(5) + 8 + ab(1)
      ab(1)=xb*ab(1)
      ab(3)= - l6 - ab(10)
      ab(5)=2*ab(7)
      ab(3)=ab(3)*ab(5)
      ab(3)=ab(3) + ab(9)
      ab(3)=ab(3)*ab(7)
      ab(6)=2*l6
      ab(2)=ab(3) + ab(6) - ab(2)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=ab(12) + l6
      ab(3)=ab(5)*ab(2)*ab(8)
      ab(4)=ab(4) - ab(6)
      ab(3)=ab(3) + ab(4)
      ab(3)=ab(3)*ab(7)
      ab(1)=ab(1) + ab(3) + ab(4)
      ab(1)=xb*ab(1)
      ab(2)=ab(2)*ab(15)
      ab(2)= - 2*ab(2) - ab(4)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=25*l6
      ab(2)=3*l3
      ab(3)= - l4 + ab(2) - ab(1) + 4*l10
      ab(4)=3*l2
      ab(5)= - ab(4) + 17*l8
      ab(6)=9*l5
      ab(5)= - 2*ab(5) - ab(6) + ab(3) + 25*l11
      ab(5)= - 9*l7 + 2*ab(5)
      ab(5)=ab(5)*l7
      ab(7)=2*l8
      ab(8)=ab(7) - ab(4)
      ab(9)=4*l11
      ab(10)=2*l4
      ab(8)=ab(9) - ab(10) + 3*ab(8)
      ab(8)= - l3 + 2*ab(8)
      ab(8)=ab(8)*l3
      ab(11)=4*l4 + 4*l2
      ab(12)=9*l8
      ab(13)= - 2 - 9*l11 + ab(12) - ab(11)
      ab(14)=11*l6
      ab(13)=l3 + ab(14) + l10 + 2*ab(13)
      ab(15)=2*l10
      ab(13)=ab(13)*ab(15)
      ab(12)=ab(12) + 2
      ab(16)=7*l2
      ab(12)=ab(16) + 2*ab(12)
      ab(12)= - 29*l11 + 2*ab(12)
      ab(12)=ab(12)*l11
      ab(3)=32*l11 - 34*l8 + ab(3) - l2
      ab(3)= - ab(6) + 2*ab(3)
      ab(3)=ab(3)*l5
      ab(6)=l8 - 1
      ab(17)=2*ab(6)
      ab(18)=ab(17) - ab(4)
      ab(18)=2*ab(18)
      ab(19)= - ab(18) + 7*l11
      ab(19)=13*l4 + 2*ab(19)
      ab(19)=ab(19)*l4
      ab(20)=d16 - d13
      ab(21)=7*d4
      ab(22)= - 8*d11 - ab(21) + 7*ab(20)
      ab(23)=4 + 5*l8
      ab(24)=4*l8
      ab(23)=ab(23)*ab(24)
      ab(25)=4*ab(6)
      ab(26)=ab(25) - ab(4)
      ab(27)=2*l2
      ab(26)=ab(26)*ab(27)
      ab(28)=d15 - d14
      ab(29)=Pi**2
      ab(30)=6*ab(29)
      ab(31)=7*d7
      ab(32)=2*d8
      ab(3)=ab(5) - 2*ab(22) - 64*d17 - 32*d9 - 46*d10 - ab(23) - 
     & ab(26) + ab(30) - ab(31) + ab(19) - ab(32) + ab(13) + ab(12) + 
     & ab(8) + ab(3) - 14*ab(28)
      ab(5)=13*l6
      ab(8)=3*l8
      ab(12)= - ab(10) + 9*l10 + ab(5) - ab(8)
      ab(13)=ab(2) + ab(12) - 7
      ab(19)=17*l5
      ab(22)=ab(13) - ab(19)
      ab(23)=17*l7
      ab(22)= - ab(23) + 2*ab(22)
      ab(22)=ab(22)*l7
      ab(13)= - ab(19) + 2*ab(13)
      ab(13)=ab(13)*l5
      ab(33)=ab(17) + l4
      ab(33)=ab(33)*ab(10)
      ab(34)=l8 - 2
      ab(34)= - l11 + 2*ab(34)
      ab(34)=ab(34)*l11
      ab(35)=l8 + 1
      ab(36)=2*ab(35)
      ab(37)=ab(36) - l3
      ab(37)=ab(37)*l3
      ab(13)=ab(22) - ab(34) + ab(33) + ab(13) + ab(37)
      ab(22)=d12 - d6
      ab(33)=ab(22) - d7
      ab(34)=d10 + d9
      ab(37)= - 2*ab(28)
      ab(33)=d8 + d17 + ab(37) - 1.0_ki/3.0_ki*ab(29) + 4*d4 + 10*ab(34) + 
     & 2*ab(33)
      ab(34)=2*d5
      ab(38)=ab(34) + 1
      ab(38)= - d11 + 2*ab(38)
      ab(39)= - l8*ab(35)
      ab(40)=ab(17) - l2
      ab(41)=l2*ab(40)
      ab(39)=ab(41) + ab(39) + ab(38) + ab(33)
      ab(41)=l11 + 4*l6
      ab(42)=ab(41) - ab(24)
      ab(43)=ab(42) + l3
      ab(44)=11 - ab(43)
      ab(45)=2*l6
      ab(44)=ab(44)*ab(45)
      ab(46)=l3 + l8
      ab(47)=ab(46) - 5
      ab(48)=3*l6
      ab(49)=ab(48) - ab(47)
      ab(50)=5*l10
      ab(49)=2*ab(49) - ab(50)
      ab(49)=l10*ab(49)
      ab(39)=ab(49) + ab(44) + 2*ab(39) + ab(13)
      ab(44)=xa**2
      ab(49)=2*ab(44)
      ab(39)=ab(39)*ab(49)
      ab(51)=2*l3
      ab(52)= - ab(10) - ab(51) + 11*l11
      ab(53)=1 + ab(7)
      ab(53)=4*ab(53) - l2
      ab(53)=2*ab(53) + ab(52)
      ab(53)=2*ab(53) - ab(5)
      ab(53)=l6*ab(53)
      ab(39)=ab(39) + ab(53) + ab(3)
      ab(39)=ab(39)*ab(44)
      ab(53)=3*l11
      ab(54)=15*l6 - 3*l4 + ab(2) + 14*l10
      ab(55)=25*l5
      ab(56)= - ab(55) + ab(53) + ab(54)
      ab(57)=ab(4) + 7
      ab(57)= - ab(56) + 2*ab(57)
      ab(58)=25*l7
      ab(57)=ab(58) + 2*ab(57)
      ab(57)=ab(57)*l7
      ab(54)=ab(9) + ab(54)
      ab(59)=l2 + 2
      ab(60)= - ab(54) + 7*ab(59)
      ab(60)=ab(55) + 2*ab(60)
      ab(60)=ab(60)*l5
      ab(59)=ab(59) + ab(10)
      ab(59)= - l3 + 2*ab(59)
      ab(59)=ab(59)*l3
      ab(61)= - l11 + 6*l2
      ab(62)=9*l4
      ab(61)=ab(62) + 2*ab(61)
      ab(61)=ab(61)*l4
      ab(63)=l2**2
      ab(61)=ab(61) + 2*ab(63)
      ab(63)=l2 - 2
      ab(64)=2*ab(63)
      ab(65)=ab(64) - l11
      ab(65)=ab(65)*l11
      ab(66)=4*d17 + 8*d12 - 6*ab(28)
      ab(65)=ab(65) - ab(66)
      ab(57)=ab(57) + ab(61) - ab(65) + ab(60) - ab(59)
      ab(59)=ab(32) - d7
      ab(60)=ab(59) + 8.0_ki/3.0_ki*ab(29) + 14*d10 + 16*d9
      ab(67)=ab(10) - l11
      ab(68)= - 8 - l2
      ab(68)=2*ab(68) - ab(67)
      ab(68)=2*ab(68) + ab(48)
      ab(68)=l6*ab(68)
      ab(69)=5*l6
      ab(70)=ab(2) + 6*l10
      ab(71)=ab(69) + ab(70)
      ab(11)= - l11 + ab(11)
      ab(72)=ab(11) + 3
      ab(72)=2*ab(72)
      ab(73)= - ab(72) + ab(71)
      ab(73)=ab(73)*ab(15)
      ab(74)= - ab(20) + 4*d6
      ab(75)=ab(74) + d4
      ab(76)= - 12 - ab(75)
      ab(39)=ab(39) + ab(73) + ab(68) + 2*ab(76) - ab(60) + ab(57)
      ab(39)=xa*ab(39)
      ab(31)=ab(32) - ab(31) - 22.0_ki/3.0_ki*ab(29) + 42*d10 + 40*d9
      ab(21)=ab(21) - ab(74)
      ab(32)=ab(34) - 1
      ab(32)= - d11 + 2*ab(32)
      ab(34)= - 2*ab(32) - ab(21)
      ab(68)= - ab(67) + 4*l3
      ab(73)=ab(27) - ab(68)
      ab(73)=2*ab(73) - ab(48)
      ab(73)=l6*ab(73)
      ab(74)=ab(70) + ab(1)
      ab(72)=ab(72) - ab(74)
      ab(72)=ab(72)*ab(15)
      ab(34)=ab(72) + ab(73) + 2*ab(34) - ab(31) - ab(57)
      ab(34)=ab(34)*ab(44)
      ab(57)= - 8*l8 + l2
      ab(52)=2*ab(57) - ab(52)
      ab(5)=2*ab(52) + ab(5)
      ab(5)=l6*ab(5)
      ab(3)=ab(34) + ab(5) - ab(3)
      ab(3)=ab(3)*ab(44)
      ab(5)= - 5.0_ki/3.0_ki*ab(29) + 8*d10 + 4*ab(28)
      ab(28)=ab(35)*ab(7)
      ab(29)= - ab(40)*ab(27)
      ab(34)=ab(2) - 1 + ab(42)
      ab(34)=ab(34)*ab(45)
      ab(52)=7*l6
      ab(47)=ab(52) + ab(47)
      ab(47)=2*ab(47) + ab(50)
      ab(47)=l10*ab(47)
      ab(22)=2*d9 - ab(22)
      ab(57)= - 3 + ab(22)
      ab(57)=2*ab(57) - d17
      ab(13)=ab(47) + ab(34) + ab(29) + ab(28) + 2*ab(57) + ab(5) - 
     & ab(13)
      ab(3)=2*ab(13) + ab(3)
      ab(3)=xa*ab(3)
      ab(13)=ab(2) + 7
      ab(12)=ab(13) + ab(12)
      ab(28)=ab(12) - ab(19)
      ab(23)= - ab(23) + 2*ab(28)
      ab(23)=ab(23)*l7
      ab(12)= - ab(19) + 2*ab(12)
      ab(12)=ab(12)*l5
      ab(19)=ab(17) - l3
      ab(19)=ab(19)*l3
      ab(12)=ab(23) + ab(12) + ab(19)
      ab(19)=ab(22) + 3
      ab(19)= - d17 + 2*ab(19)
      ab(5)=ab(5) - ab(12) + 2*ab(19)
      ab(19)=ab(46) + 3
      ab(22)=ab(19) + ab(52)
      ab(22)=ab(50) + 2*ab(22)
      ab(22)=ab(22)*l10
      ab(13)=ab(42) + ab(13)
      ab(13)=ab(13)*ab(45)
      ab(23)=ab(7)*ab(6)
      ab(28)=ab(7) + l4
      ab(28)=ab(28)*ab(10)
      ab(29)=ab(7) - l11
      ab(29)=ab(29)*l11
      ab(28)=ab(28) - ab(29)
      ab(29)=ab(7) - l2
      ab(34)=ab(29)*ab(27)
      ab(13)=ab(22) + ab(13) + ab(5) + ab(23) - ab(28) - ab(34)
      ab(13)=2*ab(13)
      ab(22)= - ab(2) + ab(25) - ab(41)
      ab(22)=ab(22)*ab(45)
      ab(23)=ab(36) + l4
      ab(23)=ab(23)*ab(10)
      ab(17)=ab(17) - l11
      ab(17)=ab(17)*l11
      ab(17)=ab(23) - ab(17)
      ab(23)=l8 + 2
      ab(7)= - ab(23)*ab(7)
      ab(25)=ab(36) - l2
      ab(34)=ab(25)*ab(27)
      ab(36)=ab(46) + 5
      ab(41)= - ab(52) - ab(36)
      ab(41)=2*ab(41) - ab(50)
      ab(41)=l10*ab(41)
      ab(5)=ab(41) + ab(22) + ab(34) + ab(7) - ab(5) + ab(17)
      ab(7)=ab(4) - 7
      ab(7)= - ab(56) + 2*ab(7)
      ab(7)=ab(58) + 2*ab(7)
      ab(7)=ab(7)*l7
      ab(22)= - ab(54) + 7*ab(63)
      ab(22)=ab(55) + 2*ab(22)
      ab(22)=ab(22)*l5
      ab(34)=ab(63) + ab(10)
      ab(34)= - l3 + 2*ab(34)
      ab(34)=ab(34)*l3
      ab(7)=ab(7) + ab(22) - ab(34)
      ab(22)=ab(4) - 2
      ab(22)= - l11 + 2*ab(22)
      ab(22)=ab(62) + 2*ab(22)
      ab(22)=ab(22)*l4
      ab(34)=l2 - 4
      ab(41)=ab(34)*ab(27)
      ab(22)=ab(41) - ab(65) + ab(22) + ab(7)
      ab(21)=ab(21) + 2*ab(38)
      ab(21)=ab(31) + 2*ab(21)
      ab(31)=l2 - 10
      ab(31)= - ab(68) + 2*ab(31)
      ab(31)= - ab(48) + 2*ab(31)
      ab(31)=ab(31)*l6
      ab(38)=ab(11) - 5
      ab(38)= - ab(70) + 2*ab(38)
      ab(1)=ab(38) - ab(1)
      ab(1)=ab(1)*ab(15)
      ab(1)=ab(31) + ab(1) - ab(21) - ab(22)
      ab(1)=ab(1)*ab(44)
      ab(31)= - l11 + 2*ab(34)
      ab(31)=ab(31)*l11
      ab(7)= - ab(31) + ab(61) - 12*l8 + ab(7) + ab(66)
      ab(31)=15 - l2
      ab(31)=2*ab(31) + ab(68)
      ab(31)=2*ab(31) + ab(48)
      ab(31)=l6*ab(31)
      ab(34)=ab(11) - 3
      ab(34)=2*ab(34)
      ab(41)= - ab(34) + ab(74)
      ab(41)=ab(41)*ab(15)
      ab(21)= - ab(1) + ab(41) + ab(31) + ab(21) + ab(7)
      ab(21)=ab(21)*ab(44)
      ab(5)=2*ab(5) + ab(21)
      ab(5)=ab(5)*ab(44)
      ab(21)= - l8 - l10 + l4 + l2
      ab(31)= - ab(45) - ab(21)
      ab(1)=8*ab(31) + ab(1)
      ab(1)=ab(1)*ab(44)
      ab(1)=ab(13) + ab(1)
      ab(1)=xb*xa*ab(1)
      ab(1)=ab(1) - ab(13) + ab(5)
      ab(1)=xb*ab(1)
      ab(1)=ab(3) + ab(1)
      ab(1)=xb*ab(1)
      ab(3)=ab(8) + 14
      ab(4)=ab(3) + ab(4)
      ab(5)= - 32*l10 - 41*l6 - 9*l3 + 7*l4
      ab(8)=59*l5
      ab(4)=2*ab(4) - ab(53) + ab(5) + ab(8)
      ab(4)=59*l7 + 2*ab(4)
      ab(4)=ab(4)*l7
      ab(3)=ab(5) + 2*ab(3) + ab(16) - ab(9)
      ab(3)=ab(8) + 2*ab(3)
      ab(3)=ab(3)*l5
      ab(5)=l2 + ab(10) + 2*ab(23)
      ab(2)= - ab(2) + 2*ab(5)
      ab(2)=ab(2)*l3
      ab(5)= - l8 + ab(11) + 8
      ab(8)=19*l6
      ab(5)= - ab(8) - 11*l10 - 5*l3 + 2*ab(5)
      ab(5)=ab(5)*ab(15)
      ab(9)=ab(18) + l11
      ab(9)= - 5*l4 + 2*ab(9)
      ab(9)=ab(9)*l4
      ab(11)=ab(35)*ab(24)
      ab(13)= - l11 + 2*ab(40)
      ab(13)=ab(13)*l11
      ab(2)= - ab(4) - ab(3) + ab(9) - ab(13) - ab(11) + ab(26) + 
     & ab(30) + ab(59) + ab(37) - 2*d10 + ab(2) + ab(5)
      ab(3)=l7 + l5
      ab(3)=ab(51) - ab(67) - ab(27) - l8 + 8*l10 - 14*ab(3)
      ab(4)= - ab(14) - 4 - ab(3)
      ab(4)=ab(4)*ab(44)
      ab(5)= - ab(10) + ab(53) + 6*l3
      ab(9)=ab(24) + l2
      ab(10)= - 15 - ab(9)
      ab(10)=2*ab(10) + ab(5)
      ab(10)=2*ab(10) + ab(8)
      ab(10)=l6*ab(10)
      ab(11)=ab(20) - d4
      ab(13)= - 8 + ab(11)
      ab(4)=4*ab(4) + ab(10) + 2*ab(13) - ab(2)
      ab(4)=ab(4)*ab(44)
      ab(9)=5 + ab(9)
      ab(5)=2*ab(9) - ab(5)
      ab(5)=2*ab(5) - ab(8)
      ab(5)=l6*ab(5)
      ab(8)=24 - ab(11)
      ab(2)=ab(4) + ab(5) + 2*ab(8) + ab(2)
      ab(2)=ab(2)*ab(44)
      ab(3)=9*l6 + 12 + ab(3)
      ab(1)=ab(1) + 4*ab(3) + ab(2)
      ab(1)=xb*ab(1)
      ab(1)=ab(39) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=ab(33) + ab(32)
      ab(3)=ab(6)*l8
      ab(4)=ab(29)*l2
      ab(3)= - ab(2) + ab(3) - ab(4)
      ab(4)=ab(19) - ab(48)
      ab(4)=ab(50) + 2*ab(4)
      ab(4)=ab(4)*l10
      ab(5)=ab(43) + 1
      ab(5)=ab(5)*ab(45)
      ab(3)=ab(4) + 2*ab(3) + ab(5) - ab(28) - ab(12)
      ab(3)=ab(3)*ab(44)
      ab(4)=l8*ab(23)
      ab(5)= - l2*ab(25)
      ab(2)=ab(5) + ab(4) - ab(2)
      ab(4)=ab(43)*ab(45)
      ab(5)= - ab(48) + ab(36)
      ab(5)=2*ab(5) + ab(50)
      ab(5)=l10*ab(5)
      ab(2)=ab(3) + ab(5) + ab(4) + 2*ab(2) - ab(12) - ab(17)
      ab(2)=ab(2)*ab(49)
      ab(4)=ab(75) - 12
      ab(4)=ab(60) + 2*ab(4)
      ab(5)= - 1 + l2
      ab(5)=2*ab(5) + ab(67)
      ab(5)=2*ab(5) - ab(48)
      ab(5)=l6*ab(5)
      ab(6)=ab(34) - ab(71)
      ab(6)=ab(6)*ab(15)
      ab(2)=ab(2) + ab(6) + ab(5) + ab(4) - ab(7)
      ab(2)=ab(2)*ab(44)
      ab(5)=ab(64) + ab(67)
      ab(5)= - ab(48) + 2*ab(5)
      ab(5)=ab(5)*l6
      ab(6)=ab(38) - ab(69)
      ab(6)=ab(6)*ab(15)
      ab(4)=ab(5) + ab(6) + ab(4) - ab(22)
      ab(1)=ab(1) + ab(2) + ab(4)
      ab(1)=xb*ab(1)
      ab(2)=4*ab(21) - ab(3)
      ab(2)=ab(2)*ab(49)
      ab(2)=ab(2) - ab(4)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

