      subroutine AGTBa1a1sscarggg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      complex(ki) ab(77)
      include 'AGTBa1a1sscar_functions.h'

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(1,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xb - xa

      tmp = 12*CA**2*l1*z6*z7*z9*xa**3*ab(1)
      res(1,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l2 + l4 - l3
      ab(2)= - xb*ab(1)
      ab(2)= - 1 + ab(2)
      ab(3)=7*l1
      ab(4)=xb*ab(3)
      ab(2)=12*ab(2) + ab(4)
      ab(4)=2*l1
      ab(2)=ab(2)*ab(4)
      ab(1)=12*ab(1) - ab(3)
      ab(1)=ab(1)*ab(4)
      ab(3)=Pi**2
      ab(1)=ab(3) + ab(1)
      ab(1)=xa*ab(1)
      ab(3)= - xb*ab(3)
      ab(1)=ab(1) + ab(3) + ab(2)

      tmp = CA**2*z6*z7*z9*xa**3*ab(1)
      res(1,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l3 - l2
      ab(2)= - l4 + 2*ab(1)
      ab(3)=24*l4
      ab(2)=ab(2)*ab(3)
      ab(3)=2*l2
      ab(4)=ab(3) - l3
      ab(5)=24*l3
      ab(4)=ab(4)*ab(5)
      ab(5)=l2**2
      ab(6)=Pi**2
      ab(2)=ab(2) + ab(4) - 24*ab(5) + 1.0_ki/3.0_ki*ab(6)
      ab(4)= - xb*ab(2)
      ab(1)=ab(1) - l4
      ab(5)=xb*ab(1)
      ab(5)= - 1 + ab(5)
      ab(7)=8.0_ki/3.0_ki*l1
      ab(8)=xb*ab(7)
      ab(5)=7*ab(5) + ab(8)
      ab(8)=4*l1
      ab(5)=ab(5)*ab(8)
      ab(9)=2*l3
      ab(10)= - ab(9) + l5 + l4 + l2
      ab(4)=ab(5) + 48*ab(10) + ab(4)
      ab(4)=l1*ab(4)
      ab(1)= - 7*ab(1) - ab(7)
      ab(1)=ab(1)*ab(8)
      ab(1)=ab(1) + ab(2)
      ab(1)=l1*ab(1)
      ab(2)= - 2*l4 + ab(9) - ab(3)
      ab(2)=ab(2)*ab(6)
      ab(2)=ab(2) - 25.0_ki/2.0_ki*zeta3
      ab(1)=ab(1) + ab(2)
      ab(1)=xa*ab(1)
      ab(2)= - xb*ab(2)
      ab(1)=ab(1) + ab(4) + 2*ab(6) + ab(2)

      tmp = CA**2*z6*z7*z9*xa**3*ab(1)
      res(1,0,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 10*CA**2*z1*z2*xb**3
      res(0,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - l6 - l8 + l7
      ab(1)=5*ab(1) + 3*l1

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
      ab(3)=12*ab(3) + 7*l1
      ab(3)=l1*ab(2)*ab(3)
      ab(1)=10*ab(1) + ab(3)
      ab(2)=Pi**2*ab(2)

      tmp = 2*ab(1) + 5.0_ki/3.0_ki*ab(2)
      res(0,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=5*l6
      ab(2)= - ab(1) + 6*l1
      ab(3)=4*l6
      ab(2)=ab(2)*ab(3)
      ab(4)=Pi**2
      ab(5)=l1**2
      ab(6)=5.0_ki/3.0_ki*ab(4) + 14*ab(5)
      ab(2)=ab(2) - ab(6)
      ab(7)=3*l1
      ab(1)=ab(7) - ab(1)
      ab(8)=5*l7
      ab(9)= - 2*ab(1) - ab(8)
      ab(10)=4*l7
      ab(9)=ab(9)*ab(10)
      ab(8)= - 5.0_ki/3.0_ki*l8 + ab(8) + ab(1)
      ab(8)=l8*ab(8)
      ab(8)=4*ab(8) + ab(9) + ab(2)
      ab(8)=l8*ab(8)
      ab(7)=ab(7) - 5.0_ki/3.0_ki*l6
      ab(3)=ab(7)*ab(3)
      ab(3)=ab(3) - ab(6)
      ab(3)=l6*ab(3)
      ab(1)=5.0_ki/3.0_ki*l7 + ab(1)
      ab(1)=ab(1)*ab(10)
      ab(1)=ab(1) - ab(2)
      ab(1)=l7*ab(1)
      ab(1)=ab(8) + ab(3) + ab(1)
      ab(2)=5*ab(4) + 32*ab(5)
      ab(2)=l1*ab(2)
      ab(1)= - 41.0_ki/2.0_ki*zeta3 + 1.0_ki/3.0_ki*ab(2) + 2*ab(1)

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

      tmp = 5*CA**2*z1*z6*z7*z8*ab(1)
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

      tmp = 5*CA**2*z1*z6*z7*z8*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=3*ab(3)
      ab(5)=xa*ab(3)
      ab(6)= - ab(4) + 2*ab(5)
      ab(6)=xa*ab(6)
      ab(6)= - 5*ab(3) + ab(6)
      ab(7)= - 4*xb + 13
      ab(7)=ab(3)*ab(7)
      ab(7)= - ab(5) + ab(7)
      ab(7)=xb*ab(7)
      ab(6)=2*ab(6) + ab(7)
      ab(6)=xb*ab(6)
      ab(7)= - 7*ab(3) + ab(5)
      ab(7)=xa*ab(7)
      ab(7)=10*ab(3) + ab(7)
      ab(7)=xa*ab(7)
      ab(6)=ab(7) + ab(6)
      ab(7)=ab(5) - ab(3)
      ab(8)=xb - 3
      ab(8)=ab(8)*xb
      ab(9)=ab(7)*ab(8)
      ab(4)=ab(4) - ab(5)
      ab(10)=ab(4)*xa**2
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(9) + ab(10) - ab(1)
      ab(2)=ab(2)*xb
      ab(4)=ab(4)*xa
      ab(4)=ab(4) - ab(1)
      ab(4)=ab(4)*xa
      ab(2)=ab(2) - ab(4)
      ab(9)=l1*ab(2)
      ab(6)=4*ab(6) - 5*ab(9)
      ab(6)=l1*ab(6)
      ab(9)=xb - 2
      ab(9)=ab(9)*ab(7)*xb
      ab(4)=ab(4) + ab(9)
      ab(9)=l6 + l7
      ab(9)= - 10*ab(9)
      ab(4)=ab(4)*ab(9)
      ab(2)=ab(2)*Pi**2
      ab(5)=ab(1) - ab(5)
      ab(5)=ab(5)*xa
      ab(3)= - ab(3)*ab(8)
      ab(1)=ab(3) - ab(1) - ab(5)
      ab(1)=xb*ab(1)
      ab(1)=ab(5) + ab(1)
      ab(1)=l2*ab(1)
      ab(3)= - xa + xb
      ab(3)=ab(7)*ab(3)

      tmp = 20*ab(1) + 2*ab(2) + 6*ab(3) + ab(4) + ab(6)
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
      ab(11)=2*l1
      ab(12)=ab(10)*ab(11)
      ab(13)= - 11*ab(3) + 5*ab(5)
      ab(13)=xa*ab(13)
      ab(14)=6*ab(3)
      ab(13)=ab(12) + ab(9) + ab(14) + ab(13)
      ab(15)=5*l6
      ab(15)=ab(10)*ab(15)
      ab(13)=2*ab(13) + ab(15)
      ab(13)=l6*ab(13)
      ab(12)=ab(12) + ab(15)
      ab(10)=l7*ab(10)
      ab(10)=2*ab(12) + 5*ab(10)
      ab(10)=l7*ab(10)
      ab(12)=ab(3)*xb
      ab(15)=ab(12) - ab(4)
      ab(15)=ab(15)*xb
      ab(16)=ab(1) - ab(5)
      ab(16)=ab(16)*xa
      ab(15)=ab(15) + ab(16) + ab(1)
      ab(15)=ab(15)*xb
      ab(15)=ab(15) - ab(16)
      ab(17)=20*l2 + 26*l1
      ab(17)=ab(15)*ab(17)
      ab(18)=xb - xa
      ab(8)=ab(18)*ab(8)
      ab(17)= - 9*ab(8) + ab(17)
      ab(17)=l2*ab(17)
      ab(18)=13*ab(3) - ab(5)
      ab(18)=xa*ab(18)
      ab(18)= - 11*ab(9) - 12*ab(3) + ab(18)
      ab(18)=l9*ab(18)
      ab(10)=ab(13) + ab(10) + ab(17) + ab(18)
      ab(13)=l3*ab(15)
      ab(15)= - ab(1) + ab(12)
      ab(15)=xb*ab(15)
      ab(13)=ab(13) + ab(15) + ab(1) - ab(16)
      ab(4)=xb*ab(4)
      ab(4)=ab(4) - 10*ab(3) + ab(5)
      ab(4)=xb*ab(4)
      ab(15)=4*ab(3) - 3*ab(5)
      ab(15)=xa*ab(15)
      ab(16)=8*ab(3)
      ab(4)=ab(4) + ab(16) + ab(15)
      ab(4)=xb*ab(4)
      ab(15)= - ab(3) + 1.0_ki/3.0_ki*ab(5)
      ab(17)=xa**2
      ab(18)= - ab(15)*ab(17)
      ab(19)=1.0_ki/3.0_ki*xb
      ab(20)=ab(19) - 1
      ab(20)=ab(20)*ab(9)
      ab(21)=2.0_ki/3.0_ki*ab(3)
      ab(18)=ab(20) - ab(21) + ab(18)
      ab(18)=xb*ab(18)
      ab(15)=xa*ab(15)
      ab(15)=ab(21) + ab(15)
      ab(15)=xa*ab(15)
      ab(15)=ab(15) + ab(18)
      ab(11)=ab(15)*ab(11)
      ab(15)=ab(14) - ab(5)
      ab(15)=xa*ab(15)
      ab(15)= - ab(16) + ab(15)
      ab(15)=xa*ab(15)
      ab(4)=ab(11) + ab(15) + ab(4)
      ab(4)=l1*ab(4)
      ab(4)=5*ab(4) + 12*ab(13)
      ab(4)=l1*ab(4)
      ab(11)=5.0_ki/2.0_ki*ab(5)
      ab(12)= - 13.0_ki/2.0_ki*ab(12) + 17*ab(3) + ab(11)
      ab(12)=ab(12)*ab(19)
      ab(5)= - ab(14) + 13.0_ki/6.0_ki*ab(5)
      ab(5)=xa*ab(5)
      ab(5)=ab(12) - 8.0_ki/3.0_ki*ab(3) + ab(5)
      ab(5)=xb*ab(5)
      ab(12)=xb - 3
      ab(9)=ab(12)*ab(9)
      ab(6)=ab(6)*ab(17)
      ab(6)=ab(9) + ab(6)
      ab(1)=ab(1) - ab(6)
      ab(1)=xb*ab(1)
      ab(1)=ab(2) + ab(1)
      ab(1)=l1*ab(1)
      ab(2)=ab(3) - ab(11)
      ab(2)=xa*ab(2)
      ab(2)=ab(16) + ab(2)
      ab(2)=xa*ab(2)
      ab(1)=4*ab(1) + 1.0_ki/3.0_ki*ab(2) + ab(5)
      ab(1)=ab(1)*Pi**2
      ab(2)= - ab(3) + 1.0_ki/2.0_ki*ab(6)
      ab(2)=xb*ab(2)
      ab(3)=ab(3) - 1.0_ki/2.0_ki*ab(7)
      ab(3)=xa*ab(3)
      ab(2)=ab(3) + ab(2)
      ab(2)=zeta3*ab(2)
      ab(3)=l5*ab(8)

      tmp = ab(1) + 169.0_ki/2.0_ki*ab(2) + 12*ab(3) + ab(4) + 2*ab(10)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=xa*ab(3)
      ab(5)=ab(4) - ab(3)
      ab(6)=ab(5)*xb
      ab(7)= - ab(6) + 2*ab(5)
      ab(8)=ab(7)*xb
      ab(9)=3*ab(3)
      ab(10)=ab(9) - ab(4)
      ab(11)=ab(10)*xa
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(11) - ab(1)
      ab(12)=ab(2)*xa
      ab(8)=ab(12) - ab(8)
      ab(13)=ab(8)*l1
      ab(14)=1.0_ki/3.0_ki*xb
      ab(7)=ab(7)*ab(14)
      ab(15)=1.0_ki/3.0_ki*ab(4)
      ab(16)=ab(15) - ab(3)
      ab(17)=ab(16)*xa
      ab(18)=2.0_ki/3.0_ki*ab(3)
      ab(19)=ab(17) + ab(18)
      ab(19)=ab(19)*xa
      ab(7)=ab(19) + ab(7)
      ab(20)=l6*ab(7)
      ab(21)=3*ab(6)
      ab(22)=37*ab(3) - 17*ab(4)
      ab(22)=xa*ab(22)
      ab(20)=10*ab(20) - 4*ab(13) - ab(21) - 20*ab(3) + ab(22)
      ab(20)=l6*ab(20)
      ab(22)=11*ab(3)
      ab(23)= - ab(22) + 5*ab(4)
      ab(23)=ab(23)*xa
      ab(23)=ab(23) + ab(6) + 6*ab(3)
      ab(13)= - 4*ab(23) - 5*ab(13)
      ab(13)=l1*ab(13)
      ab(13)=ab(13) + ab(20)
      ab(13)=l6*ab(13)
      ab(20)=13*ab(3)
      ab(24)=ab(20) - ab(4)
      ab(24)=ab(24)*xa
      ab(25)=12*ab(3)
      ab(24)= - ab(24) + ab(25) + 11*ab(6)
      ab(26)=2*l1
      ab(27)=3*l2
      ab(28)=ab(27) + ab(26)
      ab(28)=ab(24)*ab(28)
      ab(29)= - 17*ab(3) + ab(4)
      ab(29)=xa*ab(29)
      ab(30)=16*ab(3)
      ab(29)=15*ab(6) + ab(30) + ab(29)
      ab(29)=l6*ab(29)
      ab(31)= - ab(22) + ab(4)
      ab(31)=xa*ab(31)
      ab(32)=10*ab(3)
      ab(31)=9*ab(6) + ab(32) + ab(31)
      ab(31)=l9*ab(31)
      ab(28)=ab(31) + ab(29) + ab(28)
      ab(28)=l9*ab(28)
      ab(29)= - 4.0_ki/3.0_ki*ab(3) + ab(4)
      ab(29)=xa*ab(29)
      ab(31)=ab(32) - ab(4)
      ab(32)=ab(3)*xb
      ab(31)=1.0_ki/3.0_ki*ab(31) - ab(32)
      ab(31)=xb*ab(31)
      ab(33)=8.0_ki/3.0_ki*ab(3)
      ab(29)=ab(31) - ab(33) + ab(29)
      ab(29)=xb*ab(29)
      ab(15)= - ab(1) + ab(15)
      ab(15)=xa*ab(15)
      ab(15)=ab(33) + ab(15)
      ab(15)=xa*ab(15)
      ab(15)=ab(15) + ab(29)
      ab(29)= - ab(6) + 3*ab(5)
      ab(29)=ab(29)*xb
      ab(31)=xa**2
      ab(10)=ab(10)*ab(31)
      ab(10)=ab(29) - ab(10)
      ab(29)=ab(1) + ab(10)
      ab(29)=ab(29)*xb
      ab(12)=ab(29) + ab(12)
      ab(29)=l1*ab(12)
      ab(15)=5*ab(15) + ab(29)
      ab(15)=l1*ab(15)
      ab(29)= - 5*ab(32) + 33*ab(3) - 23*ab(4)
      ab(29)=xb*ab(29)
      ab(33)=35*ab(3) - 6*ab(4)
      ab(33)=xa*ab(33)
      ab(15)=ab(15) + ab(29) - 34*ab(3) + ab(33)
      ab(29)=l1**2
      ab(15)=ab(15)*ab(29)
      ab(13)=ab(28) + ab(15) + ab(13)
      ab(15)=32.0_ki/3.0_ki*ab(3) - 9.0_ki/2.0_ki*ab(4)
      ab(15)=xa*ab(15)
      ab(28)= - 38*ab(3) - 5.0_ki/2.0_ki*ab(4)
      ab(28)=1.0_ki/3.0_ki*ab(28) + 9.0_ki/2.0_ki*ab(32)
      ab(28)=xb*ab(28)
      ab(33)=22.0_ki/3.0_ki*ab(3)
      ab(15)=ab(28) + ab(33) + ab(15)
      ab(15)=xb*ab(15)
      ab(14)=ab(14) - 1
      ab(14)=ab(14)*xb
      ab(28)=ab(5)*ab(14)
      ab(16)=ab(16)*ab(31)
      ab(16)=ab(28) - ab(16)
      ab(18)= - ab(18) + ab(16)
      ab(18)=xb*ab(18)
      ab(18)=ab(19) + ab(18)
      ab(18)=l1*ab(18)
      ab(19)=ab(1) + 5.0_ki/6.0_ki*ab(4)
      ab(19)=xa*ab(19)
      ab(19)= - ab(33) + ab(19)
      ab(19)=xa*ab(19)
      ab(15)=13*ab(18) + ab(19) + ab(15)
      ab(15)=l1*ab(15)
      ab(18)=7*ab(4)
      ab(19)= - ab(32) + 9*ab(3) - ab(18)
      ab(19)=xb*ab(19)
      ab(28)=5*l6
      ab(31)=ab(7)*ab(28)
      ab(33)=1.0_ki/3.0_ki*ab(3)
      ab(16)= - ab(33) + 1.0_ki/2.0_ki*ab(16)
      ab(16)=xb*ab(16)
      ab(17)=ab(33) + 1.0_ki/2.0_ki*ab(17)
      ab(17)=xa*ab(17)
      ab(16)=ab(17) + ab(16)
      ab(17)=Pi**2
      ab(16)=ab(16)*ab(17)
      ab(33)= - ab(20) + 10*ab(4)
      ab(33)=xa*ab(33)
      ab(15)=611.0_ki/120.0_ki*ab(16) + ab(31) + ab(15) + ab(19) + ab(1) + 
     & ab(33)
      ab(15)=ab(15)*ab(17)
      ab(16)=ab(32) - ab(9)
      ab(16)=ab(16)*xb
      ab(19)=ab(1) - ab(4)
      ab(19)=ab(19)*xa
      ab(31)=ab(19) + ab(1)
      ab(16)=ab(16) + ab(31)
      ab(16)=ab(16)*xb
      ab(16)=ab(16) - ab(19)
      ab(33)=ab(16)*l1
      ab(34)= - ab(3) - ab(4)
      ab(34)=xa*ab(34)
      ab(21)=ab(21) + ab(1) + ab(34)
      ab(14)=ab(3)*ab(14)
      ab(14)=ab(14) + 1.0_ki/3.0_ki*ab(31)
      ab(14)=ab(14)*xb
      ab(14)=ab(14) - 1.0_ki/3.0_ki*ab(19)
      ab(19)=l2*ab(14)
      ab(19)= - 80*ab(19) + 3*ab(21) - 46*ab(33)
      ab(19)=l2*ab(19)
      ab(21)= - xb*ab(1)
      ab(18)=ab(21) + ab(22) - ab(18)
      ab(18)=xb*ab(18)
      ab(20)=ab(20) - 3*ab(4)
      ab(20)=xa*ab(20)
      ab(18)=ab(18) - ab(25) + ab(20)
      ab(18)=6*ab(18) - 25*ab(33)
      ab(18)=l1*ab(18)
      ab(20)=25*ab(3) - 9*ab(4)
      ab(20)=xa*ab(20)
      ab(20)= - 7*ab(6) - ab(30) + ab(20)
      ab(20)=l6*ab(20)
      ab(18)=ab(19) + ab(18) + 3*ab(20)
      ab(14)=ab(14)*ab(17)
      ab(14)=29*ab(14) + 2*ab(18)
      ab(14)=l2*ab(14)
      ab(18)=5*ab(8)
      ab(18)= - ab(29)*ab(18)
      ab(19)= - ab(28) - 4*l1
      ab(19)=l6*ab(8)*ab(19)
      ab(18)=ab(18) + 2*ab(19)
      ab(19)= - ab(26) - ab(28)
      ab(8)=ab(8)*ab(19)
      ab(7)=5*ab(7)
      ab(19)=l7*ab(7)
      ab(8)=ab(19) + ab(8)
      ab(8)=l7*ab(8)
      ab(7)=ab(17)*ab(7)
      ab(7)=4*ab(8) + 2*ab(18) + ab(7)
      ab(7)=l7*ab(7)
      ab(8)=4*ab(3)
      ab(18)= - ab(8) + ab(4)
      ab(18)=xa*ab(18)
      ab(19)=2*ab(4) - ab(32)
      ab(19)=xb*ab(19)
      ab(1)=ab(19) + ab(1) + ab(18)
      ab(1)=l1*ab(1)
      ab(18)=l6*ab(23)
      ab(1)=3*ab(1) + ab(18)
      ab(18)=ab(5)*xa
      ab(18)=ab(18) - ab(6)
      ab(19)=3*l5 - ab(27)
      ab(18)=ab(18)*ab(19)
      ab(19)= - l9*ab(24)
      ab(1)=ab(19) + 2*ab(1) + ab(18)
      ab(1)=l5*ab(1)
      ab(17)= - 10*ab(29) - ab(17)
      ab(16)=ab(16)*ab(17)
      ab(17)= - 12*l3 - 24*l2
      ab(17)=ab(17)*ab(33)
      ab(16)=ab(17) + ab(16)
      ab(16)=l3*ab(16)
      ab(17)=41.0_ki/4.0_ki*ab(4)
      ab(18)= - 225.0_ki/4.0_ki*ab(32) + 179*ab(3) - ab(17)
      ab(18)=xb*ab(18)
      ab(19)= - 92*ab(3) + 225.0_ki/4.0_ki*ab(4)
      ab(19)=xa*ab(19)
      ab(20)=133*ab(3)
      ab(18)=ab(18) - ab(20) + ab(19)
      ab(18)=xb*ab(18)
      ab(10)=ab(3) + 1.0_ki/2.0_ki*ab(10)
      ab(10)=xb*ab(10)
      ab(11)= - ab(3) + 1.0_ki/2.0_ki*ab(11)
      ab(11)=xa*ab(11)
      ab(10)=ab(11) + ab(10)
      ab(10)=l1*ab(10)
      ab(3)= - 87*ab(3) + ab(17)
      ab(3)=xa*ab(3)
      ab(3)=ab(20) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)=183*ab(10) + ab(3) + ab(18)
      ab(3)=zeta3*ab(3)
      ab(10)=plg4half*ab(12)
      ab(2)=ab(6) - ab(2)
      ab(2)=d3*ab(2)
      ab(4)= - ab(9) - ab(4)
      ab(4)=xa*ab(4)
      ab(4)=5*ab(6) + ab(8) + ab(4)
      ab(4)=d1*ab(4)
      ab(5)= - ab(6) + ab(5)
      ab(5)=d2*ab(5)

      tmp = 4*ab(1) + 48*ab(2) + ab(3) + 6*ab(4) + 96*ab(5) + ab(7) + 8
     & *ab(10) + 2*ab(13) + ab(14) + ab(15) + ab(16)
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

      tmp = 5*CA**2*z6*z10*ab(1)
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
      ab(2)= - 8*l1 - 10*l2 + 3*l4 + 2*l3

      tmp = 2*ab(2)*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - 6*l2 - 5*l1 + l4 + l3
      ab(2)=2*l4
      ab(2)=ab(1)*ab(2)
      ab(3)=5*l2 + 8*l1
      ab(4)=4*l2
      ab(5)=ab(3)*ab(4)
      ab(4)=ab(4) + 3*l1
      ab(4)= - l3 + 2*ab(4)
      ab(4)=ab(4)*l3
      ab(6)=Pi**2
      ab(7)=l1**2
      ab(2)= - ab(4) + ab(2) + ab(5) - 4.0_ki/3.0_ki*ab(6) + 13*ab(7)
      ab(5)=2*ab(2)
      ab(8)=xa*ab(2)
      ab(8)= - ab(5) + ab(8)
      ab(8)=xa*ab(8)
      ab(3)=l2*ab(3)
      ab(1)=l4*ab(1)
      ab(1)=ab(8) + 6*ab(1) - 3*ab(4) + 12*ab(3) - 4*ab(6) + 39*ab(7)
      ab(1)=xa*ab(1)
      ab(1)= - ab(5) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + ab(2)

      tmp = 2*CA**2*z6*z10*ab(1)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=CA**2*z6*z10
      ab(2)=ab(1)*xa
      ab(3)=2*ab(1)
      ab(2)=ab(2) - ab(3)
      ab(4)=1.0_ki/3.0_ki*xa
      ab(4)=ab(2)*ab(4)
      ab(4)=ab(4) + ab(1)
      ab(4)=ab(4)*xa
      ab(4)=ab(4) - 2.0_ki/3.0_ki*ab(1)
      ab(4)=ab(4)*xa
      ab(4)=ab(4) + 1.0_ki/3.0_ki*ab(1)
      ab(5)=ab(4)*l1
      ab(2)=ab(2)*xa
      ab(2)=ab(2) + 3*ab(1)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - ab(3)
      ab(2)=ab(2)*xa
      ab(1)=ab(2) + ab(1)
      ab(2)=ab(1)*l4
      ab(3)=ab(1)*l3
      ab(6)=l2*ab(1)
      ab(6)= - 43*ab(5) - 52*ab(6) + 17*ab(2) + 9*ab(3)
      ab(6)=l1*ab(6)
      ab(7)=2*l4
      ab(7)=ab(1)*ab(7)
      ab(7)=ab(7) + ab(3)
      ab(8)=3*l3
      ab(9)= - ab(7)*ab(8)
      ab(10)=5*l4
      ab(11)=8*l2
      ab(8)= - ab(11) + ab(10) + ab(8)
      ab(8)=ab(11)*ab(1)*ab(8)
      ab(11)=ab(1)*l4**2
      ab(6)=ab(6) + ab(8) - 7*ab(11) + ab(9)
      ab(6)=l1*ab(6)
      ab(8)=ab(4)*l2
      ab(9)= - ab(4)*ab(10)
      ab(5)=13*ab(5) + 16*ab(8) + ab(9) - ab(3)
      ab(5)=ab(5)*Pi**2
      ab(3)= - 10*ab(8) + 3*ab(2) + 2*ab(3)
      ab(3)=l2*ab(3)
      ab(3)=ab(11) - ab(3)
      ab(7)= - l3*ab(7)
      ab(3)=ab(7) - 2*ab(3)
      ab(3)=l2*ab(3)
      ab(4)=l3*ab(4)
      ab(2)=ab(2) + ab(4)
      ab(2)=l3*ab(2)
      ab(2)=ab(11) + ab(2)
      ab(2)=l3*ab(2)
      ab(4)=l4**3
      ab(4)= - 25*zeta3 + ab(4)
      ab(1)=ab(1)*ab(4)
      ab(1)=ab(5) + ab(6) + 4*ab(3) + ab(2) + ab(1)

      tmp = 2*ab(1)
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

      tmp = 12*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l8 + l6
      ab(2)=l4 + l2
      ab(2)= - l7 + 3*ab(2)
      ab(3)= - ab(2) - 2 + 2*l5
      ab(3)=3*ab(3)
      ab(4)=14*l3
      ab(5)=ab(3) - 11*l11 + ab(4) + l10 + 2*ab(1)
      ab(6)=xb**2
      ab(7)=6*ab(6)
      ab(8)=ab(7) - ab(5)
      ab(8)=ab(8)*ab(6)
      ab(5)=ab(8) + ab(5)
      ab(5)=ab(5)*ab(6)
      ab(8)= - 7*l10 + 4*ab(1)
      ab(4)=l11 + ab(4) - ab(8)
      ab(9)=ab(2) - 2
      ab(9)=3*ab(9)
      ab(10)=ab(4) - ab(9)
      ab(10)=ab(10)*ab(6)
      ab(2)=3*ab(2)
      ab(4)=ab(2) - ab(4)
      ab(11)= - ab(10) + ab(4)
      ab(11)=ab(11)*ab(6)
      ab(12)=l5 + l10 - ab(1)
      ab(13)=2*l3
      ab(14)=ab(12) - ab(13)
      ab(15)=1 - ab(14)
      ab(11)=6*ab(15) + ab(11)
      ab(11)=ab(11)*ab(6)
      ab(15)=6*ab(14)
      ab(10)= - 6 + ab(10)
      ab(10)=ab(10)*ab(6)
      ab(10)=ab(15) + ab(10)
      ab(10)=xa*xb*ab(10)
      ab(10)=ab(10) - ab(15) + ab(11)
      ab(10)=xa*ab(10)
      ab(11)=ab(13) + l11
      ab(1)= - ab(11) - ab(3) - 13*l10 + 10*ab(1)
      ab(3)= - ab(4)*ab(6)
      ab(3)=ab(3) - ab(1)
      ab(3)=ab(3)*ab(6)
      ab(4)=1 + ab(14)
      ab(3)=6*ab(4) + ab(3)
      ab(3)=xb*ab(3)
      ab(3)=ab(3) + ab(10)
      ab(3)=xa*ab(3)
      ab(3)=ab(3) - 6 + ab(5)
      ab(3)=xa*ab(3)
      ab(4)= - 1 - ab(12)
      ab(4)=ab(4)*ab(7)
      ab(1)=ab(4) + ab(1)
      ab(1)=ab(1)*ab(6)
      ab(4)=ab(11) - ab(8)
      ab(2)=ab(2) - ab(4)
      ab(1)=ab(1) + ab(2)
      ab(1)=xb*ab(1)
      ab(1)=ab(1) + ab(3)
      ab(1)=xa*ab(1)
      ab(3)=ab(12)*ab(6)
      ab(3)=ab(3) - 1
      ab(5)=ab(12) + ab(3)
      ab(5)=ab(5)*ab(7)
      ab(2)=ab(5) - ab(2)
      ab(2)=ab(2)*ab(6)
      ab(4)=ab(9) - ab(4)
      ab(1)=ab(1) + ab(2) - ab(4)
      ab(1)=xa*ab(1)
      ab(2)= - ab(3)*ab(7)
      ab(2)=ab(2) + ab(4)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l4 + l2
      ab(2)=l11 - 4
      ab(3)=ab(2) - l5
      ab(4)=3*l3
      ab(3)=17*l8 - 5*l7 - ab(4) - 22*l10 - 2*ab(3) + 15*ab(1)
      ab(3)=5*l6 + 2*ab(3)
      ab(5)=2*l6
      ab(3)=ab(3)*ab(5)
      ab(6)=3*l4
      ab(7)=3*l2
      ab(8)=ab(6) + ab(7)
      ab(9)=ab(8) - l11
      ab(10)=6*l3
      ab(11)=4*l5
      ab(12)= - 10*l7 - ab(10) + 16 + ab(9) + ab(11)
      ab(12)=61*l8 + 2*ab(12)
      ab(12)=ab(12)*l8
      ab(1)=6*ab(1)
      ab(13)=ab(1) - l11
      ab(14)= - l5 + ab(13) + 11
      ab(14)=2*ab(14) - 23*l3 - 13*l7 - 15*l10 + 44*l8
      ab(15)=2*l10
      ab(14)=ab(14)*ab(15)
      ab(16)=9*l2
      ab(17)=ab(16) - l11
      ab(18)=9*l4
      ab(19)=ab(17) + ab(18)
      ab(20)=6*l5
      ab(21)=ab(20) - 6
      ab(22)=ab(21) - ab(19)
      ab(23)=3*l7
      ab(22)=ab(23) + 2*ab(22)
      ab(22)=ab(22)*l7
      ab(24)=ab(21) - ab(17)
      ab(24)= - ab(18) + 2*ab(24)
      ab(24)=ab(24)*ab(6)
      ab(25)=ab(21) + l11
      ab(25)= - ab(16) + 2*ab(25)
      ab(25)=ab(25)*ab(7)
      ab(26)=3*l5
      ab(27)=ab(2) - ab(26)
      ab(27)=ab(27)*ab(11)
      ab(28)=d15 - d14
      ab(29)=d16 + d13
      ab(30)=ab(2)*l11
      ab(31)=48*d12
      ab(3)=ab(14) - ab(3) - ab(12) - ab(22) - ab(27) + 15*d11 - ab(30)
     &  - ab(31) + ab(24) + ab(25) - 6*ab(29) - 54*ab(28)
      ab(12)=8*l8
      ab(14)=8*l6 - 11*l10 + ab(12) + ab(19)
      ab(22)=ab(11) + 3
      ab(10)= - ab(10) - ab(23) - ab(22) + ab(14)
      ab(24)=xa**2
      ab(10)=ab(10)*ab(24)
      ab(25)=d10 - d9
      ab(27)= - 3*ab(25)
      ab(30)=Pi**2
      ab(32)= - ab(27) + 5*ab(30)
      ab(33)=d8 + d7
      ab(34)=24*d6
      ab(35)=ab(34) + ab(32) - 15*ab(33)
      ab(36)=d5 + d4
      ab(37)=12*ab(36)
      ab(38)=ab(37) + 18
      ab(39)= - ab(38) - ab(35)
      ab(40)=2*l5
      ab(41)= - 3 + ab(40)
      ab(41)=4*ab(41) - ab(16)
      ab(42)=2*l7
      ab(43)= - ab(42) + 21*l4
      ab(41)=2*ab(41) - ab(43)
      ab(44)=25*l3
      ab(41)=2*ab(41) - ab(44)
      ab(41)=l3*ab(41)
      ab(10)=4*ab(10) + ab(41) + 2*ab(39) - ab(3)
      ab(10)=ab(10)*ab(24)
      ab(37)=ab(37) + 6
      ab(35)=ab(37) + ab(35)
      ab(39)=ab(40) + 1
      ab(39)=4*ab(39)
      ab(41)= - ab(39) + ab(16)
      ab(41)=2*ab(41) + ab(43)
      ab(41)=2*ab(41) + ab(44)
      ab(41)=l3*ab(41)
      ab(3)=ab(10) + ab(41) + 2*ab(35) + ab(3)
      ab(3)=ab(3)*ab(24)
      ab(10)= - ab(11) + 3*l11
      ab(35)=2*l3
      ab(10)=ab(35) + ab(8) + ab(42) + 5*ab(10)
      ab(10)= - 15*l8 + 2*ab(10)
      ab(10)=ab(10)*l8
      ab(8)=ab(8) - l7
      ab(41)=ab(8) - l3
      ab(43)= - 3*l8 - 10*l5 - ab(41) + 4*l10
      ab(43)=9*l6 + 2*ab(43)
      ab(43)=ab(43)*ab(5)
      ab(44)=13*l11
      ab(45)=ab(1) - 13*l5 + ab(44) + 3
      ab(12)= - 15*l3 - ab(12) - l10 - l7 + 2*ab(45)
      ab(12)=ab(12)*ab(15)
      ab(45)=ab(18) + ab(16)
      ab(21)= - ab(21) + 11*l11
      ab(46)=ab(21) + ab(45)
      ab(46)= - ab(23) + 2*ab(46)
      ab(46)=ab(46)*l7
      ab(47)=ab(21) + ab(16)
      ab(47)=ab(18) + 2*ab(47)
      ab(47)=ab(47)*ab(6)
      ab(21)=ab(16) + 2*ab(21)
      ab(21)=ab(21)*ab(7)
      ab(48)=ab(26) - 1
      ab(44)=ab(48) - ab(44)
      ab(44)=ab(44)*ab(11)
      ab(49)= - 8 + 15*l11
      ab(49)=ab(49)*l11
      ab(28)= - 18*ab(28)
      ab(10)= - ab(12) + ab(10) + ab(43) - ab(46) - 3*d11 - ab(31) + 
     & ab(28) + ab(21) - ab(44) + ab(47) - ab(49) + 30*ab(29)
      ab(12)=2*l11
      ab(21)=ab(12) - ab(16)
      ab(21)=ab(21)*ab(7)
      ab(31)=ab(18) + 2*ab(17)
      ab(31)=ab(31)*ab(6)
      ab(43)= - ab(23) + 2*ab(19)
      ab(43)=ab(43)*l7
      ab(21)=ab(21) - ab(31) + ab(43)
      ab(31)=2*d17
      ab(43)=8*d6
      ab(44)=ab(31) - ab(43)
      ab(32)=ab(44) - ab(32) + 17*ab(33) + 24*ab(36)
      ab(46)=17*d11
      ab(47)= - ab(46) - ab(21) + 2*ab(32)
      ab(9)=ab(35) - ab(9) + 4*l7
      ab(49)=ab(9) - 8
      ab(50)=23*l8
      ab(49)= - ab(50) + 2*ab(49)
      ab(49)=ab(49)*l8
      ab(2)=ab(2) - ab(1)
      ab(51)= - l3 - ab(42) - 9*l10 + 7*l8
      ab(52)=ab(2) - ab(51)
      ab(53)=3*l6
      ab(52)= - ab(53) + 2*ab(52)
      ab(52)=ab(52)*ab(5)
      ab(28)=16*d12 - ab(28) + 2*ab(29)
      ab(29)=l11 + 4
      ab(54)=ab(29)*l11
      ab(54)=ab(54) + ab(28)
      ab(49)= - ab(54) + 20*l5 + ab(49) + ab(52)
      ab(52)=14*l7
      ab(55)=19*l4
      ab(56)=ab(52) - ab(55)
      ab(57)= - 1 + l2
      ab(57)=18*ab(57) - ab(56)
      ab(58)=7*l3
      ab(57)=2*ab(57) + ab(58)
      ab(57)=l3*ab(57)
      ab(59)=8*l10 - 18*l8 + 7*l7
      ab(2)=ab(59) + 2*ab(2)
      ab(60)=33*l3
      ab(61)= - ab(60) - ab(2)
      ab(61)=ab(61)*ab(15)
      ab(57)=ab(61) + ab(57) - ab(47) + ab(49)
      ab(57)=ab(57)*ab(24)
      ab(27)=ab(30) - 48*d17 + ab(34) + ab(27) + 3*ab(33)
      ab(34)=ab(38) - ab(27)
      ab(38)= - ab(12) + 27*l2
      ab(11)=ab(11) + 9
      ab(61)=ab(11) + ab(38)
      ab(52)= - ab(52) + 69*l4
      ab(61)=2*ab(61) + ab(52)
      ab(62)=19*l3
      ab(61)=2*ab(61) - ab(62)
      ab(61)=l3*ab(61)
      ab(34)=ab(57) + ab(61) + 2*ab(34) - ab(10)
      ab(34)=ab(34)*ab(24)
      ab(35)=ab(35) + ab(23)
      ab(57)=l6 + ab(40) - ab(35) - 13*l10 + 10*l8
      ab(61)=ab(19) + ab(57) + 4
      ab(61)=ab(61)*ab(5)
      ab(63)=l5 + 2
      ab(63)= - ab(35) + 2*ab(63)
      ab(64)=19*l8
      ab(63)=ab(64) + 2*ab(63)
      ab(63)=ab(63)*l8
      ab(65)=18*l2
      ab(66)=18*l4
      ab(67)=ab(65) + ab(66)
      ab(68)=6*l7
      ab(69)= - ab(68) + ab(67)
      ab(70)=l5 - 1
      ab(69)=ab(69)*ab(70)
      ab(48)=ab(48) - l11
      ab(48)=ab(48)*ab(40)
      ab(48)=ab(63) + ab(61) - ab(69) - ab(48) + ab(28) - 4*l11
      ab(31)= - ab(31) - 7*d11 - 16*d6 + 14*ab(33) - 2*ab(25)
      ab(61)=13*l3
      ab(63)=8*l5
      ab(69)=ab(63) - l4
      ab(71)= - 1 + ab(69)
      ab(71)=2*ab(71) - ab(61)
      ab(71)=l3*ab(71)
      ab(72)= - l5 - ab(23) + 13*l8
      ab(73)=ab(72) + 7
      ab(74)=9*l3
      ab(75)=ab(74) - ab(73)
      ab(76)=7*l10
      ab(75)=2*ab(75) + ab(76)
      ab(75)=l10*ab(75)
      ab(77)= - 6 - ab(30)
      ab(71)=ab(75) + ab(71) + 3*ab(77) + ab(31) + ab(48)
      ab(34)=2*ab(71) + ab(34)
      ab(34)=xa*ab(34)
      ab(1)=ab(1) - ab(29)
      ab(51)=ab(51) + ab(1)
      ab(51)=ab(53) + 2*ab(51)
      ab(51)=ab(51)*ab(5)
      ab(9)=ab(9) + 8
      ab(9)= - ab(50) + 2*ab(9)
      ab(9)=ab(9)*l8
      ab(9)=ab(51) - ab(9)
      ab(19)=ab(19) - 6
      ab(19)= - ab(23) + 2*ab(19)
      ab(19)=ab(19)*l7
      ab(17)=ab(17) - 6
      ab(17)=ab(18) + 2*ab(17)
      ab(6)=ab(17)*ab(6)
      ab(17)=l11 + 6
      ab(17)= - ab(16) + 2*ab(17)
      ab(17)=ab(17)*ab(7)
      ab(6)= - ab(54) - ab(9) + ab(17) + ab(19) - ab(6)
      ab(17)=ab(32) + 6
      ab(18)=ab(16) - 14
      ab(18)= - ab(56) + 2*ab(18)
      ab(18)=ab(58) + 2*ab(18)
      ab(18)=ab(18)*l3
      ab(13)=ab(13) - 7
      ab(13)= - ab(59) + 2*ab(13)
      ab(19)=ab(13) - ab(60)
      ab(19)=ab(19)*ab(15)
      ab(17)=ab(18) + ab(6) + ab(46) + ab(19) - 2*ab(17)
      ab(17)=ab(17)*ab(24)
      ab(18)=l11 + 8
      ab(18)=ab(18)*l11
      ab(9)=ab(28) + ab(63) + ab(18) + ab(9)
      ab(7)=4 - ab(7)
      ab(7)=6*ab(7) + ab(56)
      ab(7)=2*ab(7) - ab(58)
      ab(7)=l3*ab(7)
      ab(1)= - ab(59) + 2*ab(1)
      ab(18)=ab(60) - ab(1)
      ab(18)=ab(18)*ab(15)
      ab(7)= - ab(17) + ab(18) + ab(7) + ab(47) + ab(9)
      ab(7)=ab(7)*ab(24)
      ab(18)= - ab(29) + ab(57) + ab(45)
      ab(5)=ab(18)*ab(5)
      ab(18)=l5 - 2
      ab(18)= - ab(35) + 2*ab(18)
      ab(18)=ab(64) + 2*ab(18)
      ab(18)=ab(18)*l8
      ab(5)=ab(18) + ab(5) + ab(28)
      ab(18)=ab(65) + ab(66) - ab(68)
      ab(19)=l5 + 1
      ab(18)=ab(19)*ab(18)
      ab(19)= - l11 + ab(26) - 2
      ab(19)=ab(19)*ab(40)
      ab(12)=ab(19) + ab(12) - ab(5) + ab(18)
      ab(18)= - ab(39) + l4
      ab(18)=2*ab(18) + ab(61)
      ab(18)=l3*ab(18)
      ab(19)=ab(72) - 7
      ab(23)= - ab(74) + ab(19)
      ab(23)=2*ab(23) - ab(76)
      ab(23)=l10*ab(23)
      ab(28)= - 2 + ab(30)
      ab(18)=ab(23) + ab(18) + 3*ab(28) - ab(31) + ab(12)
      ab(7)=2*ab(18) + ab(7)
      ab(7)=ab(7)*ab(24)
      ab(18)=l10 - ab(70) - ab(41)
      ab(17)=12*ab(18) + ab(17)
      ab(17)=ab(17)*ab(24)
      ab(18)=l5*ab(67)
      ab(23)=ab(26) - ab(29)
      ab(23)=ab(23)*ab(40)
      ab(20)=ab(20)*l7
      ab(5)=ab(23) - ab(5) - ab(20) + ab(18)
      ab(18)=ab(72) - 4
      ab(20)=ab(18) - ab(74)
      ab(20)= - ab(76) + 2*ab(20)
      ab(20)=ab(20)*l10
      ab(23)= - ab(61) + 2*ab(69)
      ab(23)=ab(23)*l3
      ab(20)=ab(20) - ab(23) + 3*ab(30) + ab(5) - ab(31)
      ab(20)=2*ab(20)
      ab(17)= - ab(20) + ab(17)
      ab(17)=xb*xa*ab(17)
      ab(7)=ab(17) + ab(20) + ab(7)
      ab(7)=xb*ab(7)
      ab(7)=ab(34) + ab(7)
      ab(7)=xb*ab(7)
      ab(11)=ab(35) + ab(11) - ab(14)
      ab(3)=ab(7) + 4*ab(11) + ab(3)
      ab(3)=xb*ab(3)
      ab(7)=d17 + ab(33) + ab(43)
      ab(11)=6*ab(36)
      ab(14)=1 + ab(11)
      ab(14)=3*ab(14) + ab(7)
      ab(17)=ab(69) - ab(68)
      ab(20)=15 - ab(17)
      ab(20)=2*ab(20) + ab(74)
      ab(20)=l3*ab(20)
      ab(23)=5*l3
      ab(26)=ab(23) + ab(73)
      ab(26)=2*ab(26) - ab(76)
      ab(26)=l10*ab(26)
      ab(14)=ab(26) + ab(20) + 2*ab(14) - d11 - ab(48)
      ab(20)=2*ab(24)
      ab(14)=ab(14)*ab(20)
      ab(26)= - ab(37) + ab(27)
      ab(22)= - ab(22) - ab(38)
      ab(22)=2*ab(22) - ab(52)
      ab(22)=2*ab(22) + ab(62)
      ab(22)=l3*ab(22)
      ab(10)=ab(14) + ab(22) + 2*ab(26) + ab(10)
      ab(10)=ab(10)*ab(24)
      ab(14)=ab(33) + ab(44) - ab(25)
      ab(11)=ab(11) + ab(30)
      ab(22)= - ab(14) + 2*ab(11)
      ab(21)=ab(21) + d11 + 2*ab(22)
      ab(22)=ab(55) - ab(42)
      ab(25)= - 1 - ab(16)
      ab(25)=2*ab(25) - ab(22)
      ab(25)=2*ab(25) + l3
      ab(25)=l3*ab(25)
      ab(2)=ab(23) + ab(2)
      ab(2)=ab(2)*ab(15)
      ab(2)=ab(10) + ab(2) + ab(25) - ab(49) - ab(21)
      ab(2)=xa*ab(2)
      ab(2)=ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(3)=ab(11) - 9
      ab(3)= - ab(14) + 2*ab(3)
      ab(10)=ab(16) - 2
      ab(10)=ab(22) + 2*ab(10)
      ab(10)= - l3 + 2*ab(10)
      ab(10)=ab(10)*l3
      ab(11)=ab(13) - ab(23)
      ab(11)=ab(11)*ab(15)
      ab(3)=ab(10) + ab(11) + ab(6) + d11 + 2*ab(3)
      ab(6)= - l4 - ab(68) + 8*ab(70)
      ab(6)= - ab(74) + 2*ab(6)
      ab(6)=ab(6)*l3
      ab(10)=ab(18) + ab(23)
      ab(10)= - ab(76) + 2*ab(10)
      ab(10)=ab(10)*l10
      ab(11)=ab(7) + 18*ab(36)
      ab(5)=ab(5) - d11 - ab(6) + ab(10) + 2*ab(11)
      ab(5)=ab(5)*ab(24)
      ab(6)=1 - 2*ab(36)
      ab(6)=9*ab(6) - ab(7)
      ab(7)=2*ab(17) - ab(74)
      ab(7)=l3*ab(7)
      ab(10)= - ab(23) - ab(19)
      ab(10)=2*ab(10) + ab(76)
      ab(10)=l10*ab(10)
      ab(6)= - ab(5) + ab(10) + ab(7) + 2*ab(6) + d11 - ab(12)
      ab(6)=ab(6)*ab(20)
      ab(7)=8 + ab(16)
      ab(7)=2*ab(7) + ab(22)
      ab(7)=2*ab(7) - l3
      ab(7)=l3*ab(7)
      ab(1)= - ab(23) + ab(1)
      ab(1)=ab(1)*ab(15)
      ab(1)=ab(6) + ab(1) + ab(7) + ab(21) - ab(9)
      ab(1)=ab(1)*ab(24)
      ab(1)=ab(2) + ab(1) + ab(3)
      ab(1)=xb*ab(1)
      ab(2)= - l10 - ab(4) - 3 + l5 + ab(8)
      ab(2)=6*ab(2) + ab(5)
      ab(2)=ab(2)*ab(20)
      ab(2)=ab(2) - ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

