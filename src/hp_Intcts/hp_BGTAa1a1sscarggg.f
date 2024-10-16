      subroutine hp_BGTAa1a1sscarggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) tmp,hp_cli2,hp_Li3
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      complex(ki) ab(80)
      include 'hp_BGTAa1a1sscar_functions.h'

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(1,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      

      tmp =  - 12*CA**2*l1*z6*z9*xa**3
      res(1,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=Pi**2
      ab(2)=l4 + l2 - l3
      ab(2)=12*ab(2) - 7*l1
      ab(2)=l1*ab(2)
      ab(1)=ab(1) + 2*ab(2)

      tmp = CA**2*z6*z9*xa**3*ab(1)
      res(1,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l2**2
      ab(2)=2*l2
      ab(3)=ab(2) - l3
      ab(3)=l3*ab(3)
      ab(4)=l3 - l2
      ab(5)=2*ab(4) - l4
      ab(5)=l4*ab(5)
      ab(1)= - ab(5) + ab(1) - ab(3)
      ab(3)=Pi**2
      ab(4)=l4 - ab(4)
      ab(4)=7*ab(4) - 8.0_ki/3.0_ki*l1
      ab(4)=l1*ab(4)
      ab(1)=4*ab(4) + 1.0_ki/3.0_ki*ab(3) - 24*ab(1)
      ab(1)=l1*ab(1)
      ab(4)=l4 - l3
      ab(2)= - ab(2) - 2*ab(4)
      ab(2)=ab(3)*ab(2)
      ab(1)=ab(1) - 25.0_ki/2.0_ki*zeta3 + ab(2)

      tmp = CA**2*z6*z9*xa**3*ab(1)
      res(1,0,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - xa + xb

      tmp = 10*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xb - xa
      ab(2)=l5 + l7 - l6
      ab(2)=3*l1 - 5*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(1)=5 + ab(1)

      tmp = 4*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=5*l5
      ab(2)= - ab(1) + 3*l1
      ab(3)=5*l6
      ab(4)=ab(2) + ab(3)
      ab(5)=5*l7
      ab(4)= - ab(5) + 2*ab(4)
      ab(6)=4*l7
      ab(4)=ab(4)*ab(6)
      ab(3)=ab(3) + 2*ab(2)
      ab(6)=4*l6
      ab(3)=ab(3)*ab(6)
      ab(1)= - ab(1) + 6*l1
      ab(6)=4*l5
      ab(1)=ab(1)*ab(6)
      ab(6)=Pi**2
      ab(7)=l1**2
      ab(1)=ab(4) - ab(3) + ab(1) - 5.0_ki/3.0_ki*ab(6) - 14*ab(7)
      ab(3)= - xb + xa
      ab(1)=ab(1)*ab(3)
      ab(2)= - ab(5) + 10*l6 - 5*l8 + ab(2)
      ab(1)=8*ab(2) + ab(1)

      tmp = CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=CA**2*z1*z2
      ab(2)=ab(1)*xb**4
      ab(3)=l7*z4
      ab(4)=ab(2)*ab(3)
      ab(1)=ab(1)*xb**3
      ab(5)=z4*ab(1)
      ab(6)= - ab(4) + 2*ab(5)
      ab(2)=ab(2)*z4
      ab(7)=ab(2)*l6
      ab(8)=ab(7) + 2*ab(6)
      ab(8)=ab(8)*l6
      ab(6)=ab(6)*l7
      ab(6)=ab(8) - ab(6)
      ab(8)= - ab(7) + ab(4) - ab(5)
      ab(9)=ab(2)*l5
      ab(10)= - 1.0_ki/3.0_ki*ab(9) - ab(8)
      ab(10)=l5*ab(10)
      ab(10)=ab(10) - ab(6)
      ab(10)=l5*ab(10)
      ab(11)= - ab(4) + 4*ab(5)
      ab(12)= - l7*ab(11)
      ab(7)=1.0_ki/3.0_ki*ab(7) + ab(11)
      ab(7)=l6*ab(7)
      ab(7)=ab(12) + ab(7)
      ab(7)=l6*ab(7)
      ab(4)=ab(5) - 1.0_ki/3.0_ki*ab(4)
      ab(11)=l7**2
      ab(4)=ab(4)*ab(11)
      ab(4)=ab(10) + ab(4) + ab(7)
      ab(7)=2*ab(8) + ab(9)
      ab(7)=l5*ab(7)
      ab(6)=ab(7) + ab(6)
      ab(7)=ab(8) + ab(9)
      ab(8)=ab(2)*l1
      ab(9)= - 7*ab(7) + 8.0_ki/3.0_ki*ab(8)
      ab(9)=l1*ab(9)
      ab(6)=6*ab(6) + ab(9)
      ab(6)=l1*ab(6)
      ab(9)=ab(11)*ab(5)
      ab(1)=ab(3)*ab(1)
      ab(3)=ab(5)*l6
      ab(10)=ab(1) - 1.0_ki/3.0_ki*ab(3)
      ab(10)=l6*ab(10)
      ab(10)= - ab(9) + ab(10)
      ab(10)=l6*ab(10)
      ab(11)= - ab(3) + 2*ab(1)
      ab(11)=ab(11)*l6
      ab(9)=ab(11) - ab(9)
      ab(11)=ab(3) - ab(1)
      ab(12)=ab(5)*l5
      ab(13)=1.0_ki/3.0_ki*ab(12) - ab(11)
      ab(13)=l5*ab(13)
      ab(13)=ab(13) - ab(9)
      ab(13)=l5*ab(13)
      ab(14)=l7**3*ab(5)
      ab(10)=ab(13) + 1.0_ki/3.0_ki*ab(14) + ab(10)
      ab(13)=2*ab(11) - ab(12)
      ab(13)=l5*ab(13)
      ab(9)=ab(13) + ab(9)
      ab(11)=ab(11) - ab(12)
      ab(13)=ab(5)*l1
      ab(14)= - 7*ab(11) - 8.0_ki/3.0_ki*ab(13)
      ab(14)=l1*ab(14)
      ab(9)=6*ab(9) + ab(14)
      ab(9)=l1*ab(9)
      ab(9)=10*ab(10) + ab(9)
      ab(9)=xa*ab(9)
      ab(4)=ab(9) + 10*ab(4) + ab(6)
      ab(6)= - 2*ab(11) - ab(13)
      ab(6)=xa*ab(6)
      ab(6)=ab(6) - 2*ab(7) + ab(8)
      ab(6)=ab(6)*Pi**2
      ab(1)= - 2*ab(3) + ab(12) + ab(1)
      ab(1)=5*ab(1) - 3*ab(13)
      ab(3)=l8*ab(5)
      ab(1)=2*ab(1) + 5*ab(3)
      ab(1)=l8*ab(1)
      ab(3)=xa*ab(5)
      ab(2)= - ab(2) + ab(3)
      ab(2)=zeta3*ab(2)

      tmp = 8*ab(1) + 41.0_ki/2.0_ki*ab(2) + 4*ab(4) + 5.0_ki/3.0_ki*ab(6)
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
      ab(4)=xa*ab(3)
      ab(5)=ab(4) - ab(3)
      ab(6)=xb - 3
      ab(6)=ab(6)*xb
      ab(7)=ab(5)*ab(6)
      ab(8)=3*ab(3)
      ab(9)=ab(8) - ab(4)
      ab(10)=ab(9)*xa**2
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(7) + ab(10) - ab(1)
      ab(2)=ab(2)*xb
      ab(7)=ab(9)*xa
      ab(7)=ab(7) - ab(1)
      ab(7)=ab(7)*xa
      ab(2)=ab(2) - ab(7)
      ab(9)=ab(2)*Pi**2
      ab(10)= - 13*ab(3) + 5*ab(4)
      ab(10)=xa*ab(10)
      ab(11)= - 8*xb + 13
      ab(11)=ab(3)*ab(11)
      ab(11)=3*ab(4) + ab(11)
      ab(11)=xb*ab(11)
      ab(9)=ab(9) + ab(10) + ab(11)
      ab(8)= - ab(8) + 2*ab(4)
      ab(8)=xa*ab(8)
      ab(8)= - 5*ab(3) + ab(8)
      ab(10)= - 4*xb + 13
      ab(10)=ab(3)*ab(10)
      ab(10)= - ab(4) + ab(10)
      ab(10)=xb*ab(10)
      ab(8)=2*ab(8) + ab(10)
      ab(8)=xb*ab(8)
      ab(10)= - 7*ab(3) + ab(4)
      ab(10)=xa*ab(10)
      ab(10)=10*ab(3) + ab(10)
      ab(10)=xa*ab(10)
      ab(8)=ab(10) + ab(8)
      ab(2)=l1*ab(2)
      ab(2)=4*ab(8) - 5*ab(2)
      ab(2)=l1*ab(2)
      ab(8)=xb - 2
      ab(5)=ab(8)*ab(5)*xb
      ab(5)=ab(7) + ab(5)
      ab(7)=l5 + l6
      ab(7)= - 10*ab(7)
      ab(5)=ab(5)*ab(7)
      ab(4)=ab(1) - ab(4)
      ab(4)=ab(4)*xa
      ab(3)= - ab(3)*ab(6)
      ab(1)=ab(3) - ab(1) - ab(4)
      ab(1)=xb*ab(1)
      ab(1)=ab(4) + ab(1)
      ab(1)=l2*ab(1)

      tmp = 20*ab(1) + ab(2) + ab(5) + 2*ab(9)
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
      ab(13)=12*ab(3)
      ab(14)=ab(13) - 5*ab(5)
      ab(14)=xa*ab(14)
      ab(15)= - ab(13) - ab(5)
      ab(16)=ab(3)*xb
      ab(15)=2*ab(15) + 13*ab(16)
      ab(15)=xb*ab(15)
      ab(17)=6*ab(3)
      ab(14)=ab(12) + ab(15) + ab(17) + ab(14)
      ab(15)=5*l5
      ab(15)=ab(10)*ab(15)
      ab(14)=2*ab(14) + ab(15)
      ab(14)=l5*ab(14)
      ab(12)=ab(12) + ab(15)
      ab(10)=l6*ab(10)
      ab(10)=2*ab(12) + 5*ab(10)
      ab(10)=l6*ab(10)
      ab(12)=ab(16) - ab(4)
      ab(12)=ab(12)*xb
      ab(15)=ab(1) - ab(5)
      ab(15)=ab(15)*xa
      ab(12)=ab(12) + ab(15) + ab(1)
      ab(12)=ab(12)*xb
      ab(12)=ab(12) - ab(15)
      ab(18)=20*l2 + 26*l1
      ab(18)=ab(12)*ab(18)
      ab(19)=xb*ab(1)
      ab(19)=ab(19) - ab(4) - ab(5)
      ab(19)=xb*ab(19)
      ab(19)=ab(7) + ab(19)
      ab(18)=3*ab(19) + ab(18)
      ab(18)=l2*ab(18)
      ab(19)=13*ab(3) - ab(5)
      ab(19)=xa*ab(19)
      ab(13)= - 11*ab(9) - ab(13) + ab(19)
      ab(13)=l9*ab(13)
      ab(10)=ab(14) + ab(10) + ab(18) + ab(13)
      ab(13)=xb*ab(4)
      ab(13)=ab(13) - 10*ab(3) + ab(5)
      ab(13)=xb*ab(13)
      ab(14)=4*ab(3) - 3*ab(5)
      ab(14)=xa*ab(14)
      ab(18)=8*ab(3)
      ab(13)=ab(13) + ab(18) + ab(14)
      ab(13)=xb*ab(13)
      ab(14)= - ab(3) + 1.0_ki/3.0_ki*ab(5)
      ab(19)=xa**2
      ab(20)= - ab(14)*ab(19)
      ab(21)=1.0_ki/3.0_ki*xb
      ab(22)=ab(21) - 1
      ab(22)=ab(22)*ab(9)
      ab(23)=2.0_ki/3.0_ki*ab(3)
      ab(20)=ab(22) - ab(23) + ab(20)
      ab(20)=xb*ab(20)
      ab(14)=xa*ab(14)
      ab(14)=ab(23) + ab(14)
      ab(14)=xa*ab(14)
      ab(14)=ab(14) + ab(20)
      ab(11)=ab(14)*ab(11)
      ab(14)=ab(17) - ab(5)
      ab(14)=xa*ab(14)
      ab(14)= - ab(18) + ab(14)
      ab(14)=xa*ab(14)
      ab(11)=ab(11) + ab(14) + ab(13)
      ab(11)=l1*ab(11)
      ab(13)= - ab(1) + ab(16)
      ab(13)=xb*ab(13)
      ab(4)=4*ab(13) + ab(4) + ab(15)
      ab(12)=l3*ab(12)
      ab(4)=12*ab(12) + 8*ab(4) + 5*ab(11)
      ab(4)=l1*ab(4)
      ab(11)=5.0_ki/2.0_ki*ab(5)
      ab(12)= - 13.0_ki/2.0_ki*ab(16) + 17*ab(3) + ab(11)
      ab(12)=ab(12)*ab(21)
      ab(5)= - ab(17) + 13.0_ki/6.0_ki*ab(5)
      ab(5)=xa*ab(5)
      ab(5)=ab(12) - 8.0_ki/3.0_ki*ab(3) + ab(5)
      ab(5)=xb*ab(5)
      ab(12)=xb - 3
      ab(9)=ab(12)*ab(9)
      ab(6)=ab(6)*ab(19)
      ab(6)=ab(9) + ab(6)
      ab(1)=ab(1) - ab(6)
      ab(1)=xb*ab(1)
      ab(1)=ab(2) + ab(1)
      ab(1)=l1*ab(1)
      ab(2)=ab(3) - ab(11)
      ab(2)=xa*ab(2)
      ab(2)=ab(18) + ab(2)
      ab(2)=xa*ab(2)
      ab(1)=4*ab(1) + 1.0_ki/3.0_ki*ab(2) + ab(5)
      ab(1)=ab(1)*Pi**2
      ab(2)= - ab(3) + 1.0_ki/2.0_ki*ab(6)
      ab(2)=xb*ab(2)
      ab(3)=ab(3) - 1.0_ki/2.0_ki*ab(7)
      ab(3)=xa*ab(3)
      ab(2)=ab(3) + ab(2)
      ab(2)=zeta3*ab(2)
      ab(3)= - xa + xb
      ab(3)=l8*ab(8)*ab(3)

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
      ab(2)=ab(2)*xa
      ab(8)=ab(2) - ab(8)
      ab(12)=ab(8)*l1
      ab(13)=ab(3)*xb
      ab(14)= - 43*ab(13) + 85*ab(3) + ab(4)
      ab(14)=xb*ab(14)
      ab(15)=1.0_ki/3.0_ki*xb
      ab(7)=ab(7)*ab(15)
      ab(16)=1.0_ki/3.0_ki*ab(4)
      ab(17)=ab(16) - ab(3)
      ab(18)=ab(17)*xa
      ab(19)=2.0_ki/3.0_ki*ab(3)
      ab(20)=ab(18) + ab(19)
      ab(20)=ab(20)*xa
      ab(7)=ab(20) + ab(7)
      ab(21)=l5*ab(7)
      ab(22)= - 21*ab(3) + 10*ab(4)
      ab(22)=xa*ab(22)
      ab(14)=10*ab(21) - 4*ab(12) + ab(14) - 32*ab(3) + ab(22)
      ab(14)=l5*ab(14)
      ab(21)= - 13*ab(13) + 27*ab(3) - ab(4)
      ab(21)=xb*ab(21)
      ab(22)=2*ab(4)
      ab(23)= - ab(9) + ab(22)
      ab(23)=xa*ab(23)
      ab(24)=12*ab(3)
      ab(21)=ab(21) - ab(24) + ab(23)
      ab(12)=4*ab(21) - 5*ab(12)
      ab(12)=l1*ab(12)
      ab(12)=ab(12) + ab(14)
      ab(12)=l5*ab(12)
      ab(14)=13*ab(3)
      ab(21)=ab(14) - ab(4)
      ab(21)=ab(21)*xa
      ab(21)= - ab(21) + ab(24) + 11*ab(6)
      ab(23)=2*l1
      ab(24)=l2 + ab(23)
      ab(24)=ab(21)*ab(24)
      ab(25)=37*ab(6)
      ab(26)=3*ab(4)
      ab(27)= - 43*ab(3) + ab(26)
      ab(27)=xa*ab(27)
      ab(27)=ab(25) + 40*ab(3) + ab(27)
      ab(27)=l5*ab(27)
      ab(28)=11*ab(3)
      ab(29)= - ab(28) + ab(4)
      ab(29)=xa*ab(29)
      ab(30)=10*ab(3)
      ab(29)=9*ab(6) + ab(30) + ab(29)
      ab(29)=l9*ab(29)
      ab(24)=ab(29) + ab(27) + ab(24)
      ab(24)=l9*ab(24)
      ab(27)= - 4.0_ki/3.0_ki*ab(3) + ab(4)
      ab(27)=xa*ab(27)
      ab(29)=ab(30) - ab(4)
      ab(29)=1.0_ki/3.0_ki*ab(29) - ab(13)
      ab(29)=xb*ab(29)
      ab(30)=8.0_ki/3.0_ki*ab(3)
      ab(27)=ab(29) - ab(30) + ab(27)
      ab(27)=xb*ab(27)
      ab(29)= - ab(1) + ab(16)
      ab(29)=xa*ab(29)
      ab(29)=ab(30) + ab(29)
      ab(29)=xa*ab(29)
      ab(27)=ab(29) + ab(27)
      ab(29)= - ab(6) + 3*ab(5)
      ab(29)=ab(29)*xb
      ab(30)=xa**2
      ab(10)=ab(10)*ab(30)
      ab(10)=ab(29) - ab(10)
      ab(29)=ab(1) + ab(10)
      ab(29)=ab(29)*xb
      ab(2)=ab(29) + ab(2)
      ab(29)=l1*ab(2)
      ab(27)=5*ab(27) + ab(29)
      ab(27)=l1*ab(27)
      ab(29)=23*ab(4)
      ab(31)= - 16*ab(13) + 55*ab(3) - ab(29)
      ab(31)=xb*ab(31)
      ab(14)=ab(14) + 5*ab(4)
      ab(14)=xa*ab(14)
      ab(14)=ab(27) + ab(31) - 34*ab(3) + ab(14)
      ab(27)=l1**2
      ab(14)=ab(14)*ab(27)
      ab(31)= - 35*ab(3) - ab(4)
      ab(31)=xa*ab(31)
      ab(25)=ab(25) + 36*ab(3) + ab(31)
      ab(25)=d1*ab(25)
      ab(12)=ab(25) + ab(24) + ab(14) + ab(12)
      ab(14)=32.0_ki/3.0_ki*ab(3) - 9.0_ki/2.0_ki*ab(4)
      ab(14)=xa*ab(14)
      ab(24)= - 38*ab(3) - 5.0_ki/2.0_ki*ab(4)
      ab(24)=1.0_ki/3.0_ki*ab(24) + 9.0_ki/2.0_ki*ab(13)
      ab(24)=xb*ab(24)
      ab(25)=22.0_ki/3.0_ki*ab(3)
      ab(14)=ab(24) + ab(25) + ab(14)
      ab(14)=xb*ab(14)
      ab(15)=ab(15) - 1
      ab(15)=ab(15)*xb
      ab(24)=ab(5)*ab(15)
      ab(17)=ab(17)*ab(30)
      ab(17)=ab(24) - ab(17)
      ab(19)= - ab(19) + ab(17)
      ab(19)=xb*ab(19)
      ab(19)=ab(20) + ab(19)
      ab(19)=l1*ab(19)
      ab(20)=ab(1) + 5.0_ki/6.0_ki*ab(4)
      ab(20)=xa*ab(20)
      ab(20)= - ab(25) + ab(20)
      ab(20)=xa*ab(20)
      ab(14)=13*ab(19) + ab(20) + ab(14)
      ab(14)=l1*ab(14)
      ab(16)=2.0_ki/3.0_ki*ab(13) - ab(3) - ab(16)
      ab(16)=xb*ab(16)
      ab(19)=5*l5
      ab(20)=ab(7)*ab(19)
      ab(24)=1.0_ki/3.0_ki*ab(3)
      ab(17)= - ab(24) + 1.0_ki/2.0_ki*ab(17)
      ab(17)=xb*ab(17)
      ab(18)=ab(24) + 1.0_ki/2.0_ki*ab(18)
      ab(18)=xa*ab(18)
      ab(17)=ab(18) + ab(17)
      ab(18)=Pi**2
      ab(17)=ab(17)*ab(18)
      ab(24)=ab(3) + 5.0_ki/3.0_ki*ab(4)
      ab(24)=xa*ab(24)
      ab(14)=611.0_ki/120.0_ki*ab(17) + ab(20) + ab(14) + 13*ab(16) + 6*
     & ab(3) + ab(24)
      ab(14)=ab(14)*ab(18)
      ab(16)=ab(13) - ab(9)
      ab(16)=ab(16)*xb
      ab(17)=ab(1) - ab(4)
      ab(17)=ab(17)*xa
      ab(20)=ab(17) + ab(1)
      ab(16)=ab(16) + ab(20)
      ab(16)=ab(16)*xb
      ab(16)=ab(16) - ab(17)
      ab(24)=ab(16)*l1
      ab(25)= - xb*ab(9)
      ab(25)=ab(25) - 17*ab(3) + ab(29)
      ab(25)=xb*ab(25)
      ab(15)=ab(3)*ab(15)
      ab(15)=ab(15) + 1.0_ki/3.0_ki*ab(20)
      ab(15)=ab(15)*xb
      ab(15)=ab(15) - 1.0_ki/3.0_ki*ab(17)
      ab(17)=l2*ab(15)
      ab(20)= - 19*ab(3) - ab(22)
      ab(20)=xa*ab(20)
      ab(17)= - 80*ab(17) - 46*ab(24) + ab(25) + 18*ab(3) + ab(20)
      ab(17)=l2*ab(17)
      ab(1)= - xb*ab(1)
      ab(1)=ab(1) + ab(28) - 7*ab(4)
      ab(1)=xb*ab(1)
      ab(20)=5*ab(3) + ab(4)
      ab(20)=xa*ab(20)
      ab(1)=ab(1) - 8*ab(3) + ab(20)
      ab(1)=6*ab(1) - 25*ab(24)
      ab(1)=l1*ab(1)
      ab(20)= - 6*ab(13) + 25*ab(3) - 13*ab(4)
      ab(20)=xb*ab(20)
      ab(22)=7*ab(3) + ab(26)
      ab(22)=xa*ab(22)
      ab(20)=ab(20) - 16*ab(3) + ab(22)
      ab(20)=l5*ab(20)
      ab(1)=ab(17) + ab(1) + 3*ab(20)
      ab(15)=ab(15)*ab(18)
      ab(1)=29*ab(15) + 2*ab(1)
      ab(1)=l2*ab(1)
      ab(15)=5*ab(8)
      ab(15)= - ab(27)*ab(15)
      ab(17)= - ab(19) - 4*l1
      ab(17)=l5*ab(8)*ab(17)
      ab(15)=ab(15) + 2*ab(17)
      ab(17)= - ab(23) - ab(19)
      ab(8)=ab(8)*ab(17)
      ab(7)=5*ab(7)
      ab(17)=l6*ab(7)
      ab(8)=ab(17) + ab(8)
      ab(8)=l6*ab(8)
      ab(7)=ab(18)*ab(7)
      ab(7)=4*ab(8) + 2*ab(15) + ab(7)
      ab(7)=l6*ab(7)
      ab(8)=ab(6) - ab(5)
      ab(15)=l1*ab(8)
      ab(17)= - 4*ab(3) + ab(4)
      ab(17)=xa*ab(17)
      ab(9)=2*ab(6) + ab(9) + ab(17)
      ab(9)=l5*ab(9)
      ab(9)=3*ab(15) + ab(9)
      ab(5)=ab(5)*xa
      ab(5)=ab(5) - ab(6)
      ab(6)=l8 + l2
      ab(6)=3*ab(6)
      ab(5)=ab(5)*ab(6)
      ab(6)= - l9*ab(21)
      ab(5)=ab(6) + 4*ab(9) + ab(5)
      ab(5)=l8*ab(5)
      ab(6)= - 10*ab(27) - ab(18)
      ab(6)=ab(16)*ab(6)
      ab(9)= - 12*l3 - 24*l2
      ab(9)=ab(9)*ab(24)
      ab(6)=ab(9) + ab(6)
      ab(6)=l3*ab(6)
      ab(9)=41.0_ki/4.0_ki*ab(4)
      ab(13)= - 225.0_ki/4.0_ki*ab(13) + 179*ab(3) - ab(9)
      ab(13)=xb*ab(13)
      ab(4)= - 92*ab(3) + 225.0_ki/4.0_ki*ab(4)
      ab(4)=xa*ab(4)
      ab(15)=133*ab(3)
      ab(4)=ab(13) - ab(15) + ab(4)
      ab(4)=xb*ab(4)
      ab(10)=ab(3) + 1.0_ki/2.0_ki*ab(10)
      ab(10)=xb*ab(10)
      ab(11)= - ab(3) + 1.0_ki/2.0_ki*ab(11)
      ab(11)=xa*ab(11)
      ab(10)=ab(11) + ab(10)
      ab(10)=l1*ab(10)
      ab(3)= - 87*ab(3) + ab(9)
      ab(3)=xa*ab(3)
      ab(3)=ab(15) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)=183*ab(10) + ab(3) + ab(4)
      ab(3)=zeta3*ab(3)
      ab(4)=d2 - d3
      ab(4)= - 96*ab(4)
      ab(4)=ab(8)*ab(4)
      ab(2)=plg4half*ab(2)

      tmp = ab(1) + 8*ab(2) + ab(3) + ab(4) + 4*ab(5) + ab(6) + ab(7)
     &  + 2*ab(12) + ab(14)
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

      tmp = 5*CA**2*z1*z6*z10*ab(1)
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
      ab(3)= - 8*l1 - 10*l2 + 2*l3 + 3*l4
      ab(2)=ab(3)*ab(2)
      ab(2)= - 8 + ab(2)

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
      ab(7)=ab(5)*l4
      ab(8)=ab(5)*l3
      ab(9)=16*ab(6) - 3*ab(8) + 13*ab(1) - 5*ab(7)
      ab(10)=l1*ab(5)
      ab(9)=2*ab(9) + 13*ab(10)
      ab(9)=l1*ab(9)
      ab(10)=2*l4
      ab(5)= - ab(5)*ab(10)
      ab(5)=ab(5) + ab(1)
      ab(5)=10*ab(6) + 3*ab(5) - 4*ab(8)
      ab(5)=l2*ab(5)
      ab(6)= - 5*ab(1) + ab(7)
      ab(6)=ab(6)*ab(10)
      ab(7)= - 3*ab(1) + ab(7)
      ab(7)=2*ab(7) + ab(8)
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
      ab(1)=26*ab(1) + 4*ab(2) + ab(9) + 2*ab(5) + ab(6) + ab(7)

      tmp = 2*ab(1)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=8*l2
      ab(2)=xb - 1
      ab(1)=ab(1)*ab(2)
      ab(3)=ab(2)*l4
      ab(4)=ab(2)*l1
      ab(5)=ab(4) + 1
      ab(6)=ab(2)*l3
      ab(3)= - ab(6) + ab(1) - ab(3) + 7*ab(5)
      ab(3)=ab(3)*l4
      ab(7)=4*l2
      ab(7)=ab(7)*ab(2)
      ab(5)=ab(7) + 3*ab(5)
      ab(7)= - ab(6) + 2*ab(5)
      ab(7)=ab(7)*l3
      ab(3)=ab(3) + ab(7)
      ab(7)=Pi**2
      ab(8)=ab(7)*xb
      ab(9)=ab(8) - ab(7)
      ab(10)=12*l2
      ab(10)=ab(10)*ab(2)
      ab(10)=ab(10) + 3 + 20*ab(4)
      ab(11)=2*l2
      ab(12)=ab(10)*ab(11)
      ab(13)=ab(4) + 2
      ab(14)=ab(13)*l1
      ab(12)=ab(12) + 17*ab(14) + 34*l5 - 5.0_ki/3.0_ki*ab(9) - ab(3)
      ab(12)=ab(12)*l4
      ab(15)=4*l1
      ab(15)=ab(15)*ab(2)
      ab(15)=ab(15) + 1
      ab(1)=ab(1) + 3*ab(15)
      ab(1)=ab(1)*ab(11)
      ab(11)=9*l1
      ab(11)=ab(11)*ab(13)
      ab(1)=ab(1) + ab(11) + 18*l5 - ab(9)
      ab(11)= - ab(5) + 1.0_ki/3.0_ki*ab(6)
      ab(11)=ab(11)*l3
      ab(11)=ab(11) + ab(1)
      ab(11)=ab(11)*l3
      ab(13)=3 + 64*ab(4)
      ab(2)=ab(2)*l2
      ab(15)=ab(13) + 80.0_ki/3.0_ki*ab(2)
      ab(15)=ab(15)*l2
      ab(16)=9 + 26*ab(4)
      ab(16)=ab(16)*l1
      ab(17)= - 9*l5 - ab(16) + 8.0_ki/3.0_ki*ab(9)
      ab(15)=ab(15) - 2*ab(17)
      ab(15)=ab(15)*l2
      ab(17)=1 + 1.0_ki/3.0_ki*ab(4)
      ab(18)=43*l1
      ab(17)=ab(17)*ab(18)
      ab(19)=13.0_ki/3.0_ki*ab(7)
      ab(17)=ab(17) + ab(19) - 13.0_ki/3.0_ki*ab(8) + 86*l5
      ab(17)=ab(17)*l1
      ab(20)=l5**2
      ab(21)=25*zeta3
      ab(20)= - ab(21) + 43*ab(20)
      ab(21)=ab(21)*xb
      ab(11)=ab(15) - ab(12) - ab(11) - ab(19) + ab(17) + ab(21) + 
     & ab(20)
      ab(12)=2*ab(11)
      ab(15)= - xa*ab(11)
      ab(15)=ab(12) + ab(15)
      ab(15)=xa*ab(15)
      ab(5)= - 3*ab(5) + ab(6)
      ab(5)=l3*ab(5)
      ab(1)=3*ab(1) + ab(5)
      ab(1)=l3*ab(1)
      ab(5)=l2*ab(10)
      ab(3)=6*ab(5) + 51*ab(14) + 102*l5 - 5*ab(9) - 3*ab(3)
      ab(3)=l4*ab(3)
      ab(4)= - 3 - ab(4)
      ab(4)=ab(4)*ab(18)
      ab(5)=13*ab(7)
      ab(4)=ab(4) + 13*ab(8) - 258*l5 - ab(5)
      ab(4)=l1*ab(4)
      ab(6)= - 3*ab(16) - 27*l5 + 8*ab(9)
      ab(2)= - 3*ab(13) - 80*ab(2)
      ab(2)=l2*ab(2)
      ab(2)=2*ab(6) + ab(2)
      ab(2)=l2*ab(2)
      ab(6)=xb*zeta3
      ab(1)=ab(15) + ab(3) + ab(1) + ab(2) + ab(4) - 75*ab(6) - 3*
     & ab(20) + ab(5)
      ab(1)=xa*ab(1)
      ab(1)=ab(12) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) - ab(11)

      tmp = 2*CA**2*z1*z6*z10*ab(1)
      res(2,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,2,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(2,2,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xa**2
      ab(2)=3*ab(1)
      ab(3)=2 - ab(2)
      ab(3)=ab(3)*ab(1)
      ab(3)=5 + ab(3)
      ab(3)=xa*ab(3)
      ab(4)=ab(1) + 1
      ab(2)=ab(4)*ab(2)
      ab(2)= - 5 + ab(2)
      ab(2)=ab(2)*ab(1)
      ab(5)=xa**4
      ab(6)=5 - 3*ab(5)
      ab(6)=xb*xa*ab(6)
      ab(2)=ab(6) - 5 + ab(2)
      ab(2)=xb*ab(2)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)= - 1 + ab(1)
      ab(3)=ab(3)*ab(1)
      ab(2)=8*ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=5*ab(1)
      ab(6)= - 2 - ab(3)
      ab(6)=ab(6)*ab(1)
      ab(6)=3 + ab(6)
      ab(6)=xa*ab(6)
      ab(2)=ab(6) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=ab(4)*ab(3)
      ab(3)= - 3 + ab(3)
      ab(1)=ab(3)*ab(1)
      ab(1)=ab(2) - 3 + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=3 - 5*ab(5)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l7 + l5
      ab(2)=2*l3
      ab(3)= - l11 - ab(2) - 7*l10 + 10*ab(1)
      ab(4)=l4 + l2
      ab(5)=ab(4) - 2
      ab(5)=ab(3) + 3*ab(5)
      ab(6)=15*l6
      ab(7)=ab(5) - ab(6)
      ab(8)=xa**2
      ab(9)=ab(8) + 1
      ab(10)= - ab(8)*ab(7)*ab(9)
      ab(11)= - ab(2) - 3*l10 + 2*l8 + 8*ab(1)
      ab(12)=3 - ab(11)
      ab(13)=8*l6
      ab(14)=ab(13) + ab(12)
      ab(14)=2*ab(14)
      ab(10)= - ab(14) + ab(10)
      ab(10)=ab(10)*ab(8)
      ab(15)=xa**4
      ab(7)=ab(7)*ab(15)
      ab(7)=ab(14) + ab(7)
      ab(7)=xb*xa*ab(7)
      ab(7)=ab(7) - ab(14) + ab(10)
      ab(7)=xb*ab(7)
      ab(2)= - ab(2) - 3*ab(4) + l10 - l6 - 11*l11 + 16*l8 + 6*ab(1)
      ab(10)=ab(4) + 2
      ab(3)=ab(3) + 3*ab(10)
      ab(6)= - ab(6) + ab(3)
      ab(6)=ab(6)*ab(8)
      ab(6)=ab(6) - ab(2)
      ab(6)=ab(6)*ab(8)
      ab(10)=3 + ab(11)
      ab(11)=ab(13) - ab(10)
      ab(6)=2*ab(11) + ab(6)
      ab(6)=xa*ab(6)
      ab(6)=ab(6) + ab(7)
      ab(6)=xb*ab(6)
      ab(4)=ab(4) + 4
      ab(1)=3*ab(4) - 19*l6 - 13*l10 - l11 - 6*l3 + 4*l8 + 26*ab(1)
      ab(4)= - 12*ab(8) - ab(1)
      ab(4)=ab(4)*ab(8)
      ab(1)=ab(4) + ab(1)
      ab(1)=ab(1)*ab(8)
      ab(1)=ab(6) + 12 + ab(1)
      ab(1)=xb*ab(1)
      ab(4)=2*l6
      ab(6)= - ab(4) + ab(10)
      ab(7)=2*ab(8)
      ab(6)=ab(6)*ab(7)
      ab(2)=ab(6) + ab(2)
      ab(2)=ab(2)*ab(8)
      ab(6)=3*l6
      ab(2)=ab(2) + ab(6) - ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=ab(4) + ab(12)
      ab(3)=ab(7)*ab(2)*ab(9)
      ab(4)=ab(5) - ab(6)
      ab(3)=ab(3) + ab(4)
      ab(3)=ab(3)*ab(8)
      ab(1)=ab(1) + ab(3) + ab(4)
      ab(1)=xb*ab(1)
      ab(2)=ab(2)*ab(15)
      ab(2)= - 2*ab(2) - ab(4)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=8*l10
      ab(2)=9*l2
      ab(3)=2*l4
      ab(4)=6*l3
      ab(5)= - ab(1) + 9*l7 - ab(2) + ab(3) - 37*l11 + 52*l8 - ab(4) + 
     & 18*l5
      ab(6)=2*l7
      ab(5)=ab(5)*ab(6)
      ab(7)=l4 + l2
      ab(8)=3*l3
      ab(9)= - ab(8) + 26*l8
      ab(10)=4*l10
      ab(11)= - ab(10) - 24*l11 + ab(9) + ab(7)
      ab(11)=9*l5 + 2*ab(11)
      ab(12)=2*l5
      ab(11)=ab(11)*ab(12)
      ab(13)=13*l8
      ab(14)=6*l2
      ab(15)=6*l4
      ab(16)=13*l11 - ab(13) + 3 + ab(14) + ab(15)
      ab(16)= - l10 - l3 + 2*ab(16)
      ab(17)=2*l10
      ab(16)=ab(16)*ab(17)
      ab(18)=2*l11
      ab(19)=ab(18) - 7*l2 + 6*l8
      ab(20)=3*l4
      ab(19)= - ab(20) - l3 + 2*ab(19)
      ab(21)=2*l3
      ab(19)=ab(19)*ab(21)
      ab(22)=d16 - d13
      ab(23)=11*d4
      ab(24)= - 12*d11 - ab(23) + 11*ab(22)
      ab(25)=l8 - 1
      ab(26)=2*ab(25)
      ab(27)=3*l2
      ab(28)=ab(26) - ab(27)
      ab(28)=3*ab(28)
      ab(29)=10*l4
      ab(30)= - 11*l11 + ab(28) - ab(29)
      ab(30)=ab(30)*ab(3)
      ab(31)=8*l8
      ab(32)=ab(31) + 5
      ab(33)=4*l8
      ab(32)=ab(32)*ab(33)
      ab(25)= - ab(27) + 4*ab(25)
      ab(25)=ab(25)*ab(27)
      ab(34)=d15 - d14
      ab(35)=ab(13) + 2
      ab(35)=11*l2 + 2*ab(35)
      ab(35)= - 41*l11 + 2*ab(35)
      ab(35)=ab(35)*l11
      ab(36)=11*d7
      ab(37)=2*d8
      ab(38)=Pi**2
      ab(5)=ab(19) - ab(5) - 96*d17 - 48*d9 - 70*d10 - ab(32) - ab(25)
     &  + 29.0_ki/3.0_ki*ab(38) - ab(36) - ab(37) - ab(30) - ab(11) - 
     & ab(16) + ab(35) - 2*ab(24) - 22*ab(34)
      ab(11)=ab(21) + l8
      ab(16)=ab(11) - 7
      ab(19)=7*l10
      ab(16)=ab(19) + 2*ab(16)
      ab(16)=ab(16)*l10
      ab(24)=13*l10 + ab(4) - ab(20)
      ab(30)=3*l8
      ab(32)=ab(30) + 5
      ab(32)= - ab(24) + 2*ab(32)
      ab(35)=13*l7 + 26*l5
      ab(39)=ab(32) + ab(35)
      ab(39)=ab(39)*ab(6)
      ab(40)=13*l5
      ab(32)=ab(32) + ab(40)
      ab(32)=ab(32)*ab(12)
      ab(41)=l8 + 1
      ab(42)=2*ab(41)
      ab(43)=ab(42) - l3
      ab(43)=ab(43)*ab(21)
      ab(44)=ab(26) + l4
      ab(44)=ab(44)*ab(20)
      ab(45)=ab(26) - l2
      ab(46)=ab(45)*ab(27)
      ab(47)=l8 - 2
      ab(47)= - l11 + 2*ab(47)
      ab(47)=ab(47)*l11
      ab(16)=ab(16) + ab(39) + ab(32) - ab(43) - ab(44) + ab(47) - 
     & ab(46)
      ab(32)=d12 - d6
      ab(39)=ab(32) - d7
      ab(43)=d10 + d9
      ab(44)=d17 - 3*ab(34)
      ab(39)=d8 + ab(44) - 1.0_ki/3.0_ki*ab(38) + 6*d4 + 14*ab(43) + 3*
     & ab(39)
      ab(43)=2*l8
      ab(46)=ab(43) + 1
      ab(46)=ab(46)*l8
      ab(47)=2*d5
      ab(48)=ab(47) + 1
      ab(48)= - d11 + 3*ab(48)
      ab(49)= - ab(46) + ab(48) + ab(39)
      ab(50)=l7 + l5
      ab(51)= - l11 + 20*ab(50)
      ab(52)=7*l8
      ab(53)=ab(51) + ab(52)
      ab(10)= - ab(21) + ab(10) + ab(53)
      ab(54)=16 + ab(10)
      ab(55)=13*l6
      ab(54)=2*ab(54) - ab(55)
      ab(54)=l6*ab(54)
      ab(49)=ab(54) + 2*ab(49) - ab(16)
      ab(54)=xa**2
      ab(56)=2*ab(54)
      ab(49)=ab(49)*ab(56)
      ab(57)=ab(27) + ab(20)
      ab(58)=15*l10
      ab(59)=34*l5
      ab(60)= - 11*l6 - 34*l7 + ab(58) - ab(57) - ab(59) + 15*l11
      ab(13)=6 + ab(13)
      ab(13)= - ab(8) + 2*ab(13) + ab(60)
      ab(61)=2*l6
      ab(13)=ab(13)*ab(61)
      ab(13)=ab(49) + ab(13) + ab(5)
      ab(13)=ab(13)*ab(54)
      ab(49)=3*l11
      ab(62)=ab(49) - ab(2)
      ab(63)=4*l4
      ab(59)=ab(4) - ab(59) - ab(63) + ab(62) - 17*l7 + 18*l10
      ab(64)=ab(59) - 20
      ab(64)=ab(64)*ab(6)
      ab(65)= - ab(3) + ab(8) + ab(18) + 9*l10
      ab(66)=l2 + 2
      ab(67)= - ab(65) + 5*ab(66)
      ab(68)=17*l5
      ab(67)=ab(68) + 2*ab(67)
      ab(67)=ab(67)*ab(12)
      ab(69)=ab(27) + 2
      ab(15)=ab(15) - l11
      ab(69)=ab(15) + 2*ab(69)
      ab(1)=ab(1) + 5*l3
      ab(69)= - ab(1) + 2*ab(69)
      ab(69)=ab(69)*ab(17)
      ab(70)=ab(20) - l3
      ab(66)=ab(70) + 2*ab(66)
      ab(66)=ab(66)*ab(21)
      ab(71)= - l11 + 7*l4
      ab(2)=ab(71) + ab(2)
      ab(2)=ab(2)*ab(3)
      ab(72)=l2**2
      ab(2)=ab(2) + 3*ab(72)
      ab(72)=l2 - 2
      ab(73)=2*ab(72)
      ab(74)=ab(73) - l11
      ab(74)=ab(74)*l11
      ab(75)=4*d17 + 12*d12 - 10*ab(34)
      ab(74)=ab(74) - ab(75)
      ab(33)=ab(69) - ab(67) + ab(66) + ab(64) + ab(74) - ab(2) - 
     & ab(33)
      ab(64)= - l11 - ab(61) + 22*ab(50)
      ab(66)= - ab(19) + ab(64) + ab(70)
      ab(67)=ab(66) + ab(27)
      ab(69)= - 22 - ab(67)
      ab(69)=ab(69)*ab(61)
      ab(76)=ab(37) - d7
      ab(77)=ab(76) + 11.0_ki/3.0_ki*ab(38) + 22*d10 + 24*d9
      ab(78)= - ab(22) + 6*d6
      ab(79)=ab(78) + d4
      ab(80)= - 18 - ab(79)
      ab(13)=ab(13) + ab(69) + 2*ab(80) - ab(77) - ab(33)
      ab(13)=xa*ab(13)
      ab(36)=ab(37) - ab(36) - 31.0_ki/3.0_ki*ab(38) + 58*d10 + 56*d9
      ab(37)=9*l3
      ab(64)=ab(37) - ab(64) + 35*l10
      ab(69)=ab(64) - ab(57)
      ab(80)=2 - ab(69)
      ab(80)=ab(80)*ab(61)
      ab(23)=ab(23) - ab(78)
      ab(47)=ab(47) - 1
      ab(47)= - d11 + 3*ab(47)
      ab(78)= - 2*ab(47) - ab(23)
      ab(33)=ab(80) + 2*ab(78) - ab(36) + ab(33)
      ab(33)=ab(33)*ab(54)
      ab(9)= - ab(9) - ab(60)
      ab(9)=ab(9)*ab(61)
      ab(5)=ab(33) + ab(9) - ab(5)
      ab(5)=ab(5)*ab(54)
      ab(9)= - ab(44) - 4.0_ki/3.0_ki*ab(38) + 6*d10
      ab(32)=2*d9 - ab(32)
      ab(33)= - 3 + ab(32)
      ab(33)=ab(46) + 3*ab(33) + ab(9)
      ab(4)=ab(4) - ab(53) + 10*l10
      ab(44)= - 2 + ab(4)
      ab(44)=2*ab(44) + ab(55)
      ab(44)=l6*ab(44)
      ab(16)=ab(44) + 2*ab(33) + ab(16)
      ab(5)=2*ab(16) + ab(5)
      ab(5)=xa*ab(5)
      ab(16)=ab(59) + 20
      ab(16)=ab(16)*ab(6)
      ab(33)= - ab(65) + 5*ab(72)
      ab(33)=ab(68) + 2*ab(33)
      ab(33)=ab(33)*ab(12)
      ab(44)=ab(70) + ab(73)
      ab(44)=ab(44)*ab(21)
      ab(16)=ab(33) - ab(16) - ab(44)
      ab(33)=ab(27) - 2
      ab(44)=ab(71) + 3*ab(33)
      ab(44)=ab(44)*ab(3)
      ab(14)=ab(14) + ab(15)
      ab(46)=ab(14) - 7
      ab(46)= - ab(1) + 2*ab(46)
      ab(46)=ab(46)*ab(17)
      ab(53)=l2 - 4
      ab(59)=ab(53)*ab(27)
      ab(44)=ab(44) + ab(59) - ab(74) - ab(46) + ab(16)
      ab(46)=l2 - 10
      ab(46)= - ab(64) + ab(20) + 3*ab(46)
      ab(46)=ab(46)*ab(61)
      ab(23)=ab(23) + 2*ab(48)
      ab(23)=ab(36) + 2*ab(23)
      ab(36)=ab(46) - ab(23) - ab(44)
      ab(36)=ab(36)*ab(54)
      ab(15)=ab(15) + 2*ab(33)
      ab(1)= - ab(1) + 2*ab(15)
      ab(1)=ab(1)*ab(17)
      ab(15)= - l11 + 2*ab(53)
      ab(15)=ab(15)*l11
      ab(1)= - ab(1) + ab(75) - ab(15) - 16*l8 + ab(2) + ab(16)
      ab(2)=44 + ab(69)
      ab(2)=ab(2)*ab(61)
      ab(2)= - ab(36) + ab(2) + ab(23) + ab(1)
      ab(2)=ab(2)*ab(54)
      ab(15)=ab(30) - 5
      ab(15)= - ab(24) + 2*ab(15)
      ab(16)=ab(35) + ab(15)
      ab(16)=ab(16)*ab(6)
      ab(15)=ab(40) + ab(15)
      ab(15)=ab(15)*ab(12)
      ab(23)=ab(26) - l3
      ab(23)=ab(23)*ab(21)
      ab(15)= - ab(23) + ab(16) + ab(15)
      ab(16)=ab(11) + 7
      ab(16)=ab(19) + 2*ab(16)
      ab(16)=ab(16)*l10
      ab(23)=ab(42) + l4
      ab(23)=ab(23)*ab(20)
      ab(24)=ab(42) - l2
      ab(24)=ab(24)*ab(27)
      ab(33)=ab(26) - l11
      ab(33)=ab(33)*l11
      ab(16)=ab(16) - ab(23) - ab(24) + ab(33) + ab(15)
      ab(23)=ab(32) + 3
      ab(9)=ab(9) + 3*ab(23)
      ab(23)=ab(41)*ab(43)
      ab(24)= - ab(23) - ab(9)
      ab(32)= - 5 - ab(4)
      ab(32)=2*ab(32) - ab(55)
      ab(32)=l6*ab(32)
      ab(24)=ab(32) + 2*ab(24) - ab(16)
      ab(2)=2*ab(24) + ab(2)
      ab(2)=ab(2)*ab(54)
      ab(24)= - l8 + ab(7) - l10
      ab(32)= - ab(61) - ab(24)
      ab(32)=12*ab(32) + ab(36)
      ab(32)=ab(32)*ab(54)
      ab(33)=ab(43) + l4
      ab(33)=ab(33)*ab(20)
      ab(35)=ab(43) - l2
      ab(27)=ab(35)*ab(27)
      ab(11)=ab(11) + 4
      ab(11)=ab(19) + 2*ab(11)
      ab(11)=ab(11)*l10
      ab(19)=ab(43) - l11
      ab(19)=ab(19)*l11
      ab(11)=ab(33) - ab(19) + ab(27) - ab(11) - ab(15)
      ab(4)=ab(4) + 10
      ab(4)=ab(55) + 2*ab(4)
      ab(4)=ab(4)*l6
      ab(15)=ab(26)*l8
      ab(9)=ab(9) + ab(15)
      ab(4)=ab(4) - ab(11) + 2*ab(9)
      ab(4)=2*ab(4)
      ab(9)=ab(4) + ab(32)
      ab(9)=xb*xa*ab(9)
      ab(2)=ab(9) - ab(4) + ab(2)
      ab(2)=xb*ab(2)
      ab(2)=ab(5) + ab(2)
      ab(2)=xb*ab(2)
      ab(4)=ab(30) + 10
      ab(5)= - ab(18) - ab(37) + 2*ab(4) - 22*l10 + 5*ab(7)
      ab(5)=43*l5 + 2*ab(5)
      ab(5)=ab(5)*ab(12)
      ab(7)=l8 + 2
      ab(7)=l2 + 2*ab(7)
      ab(7)= - ab(8) + ab(20) + 2*ab(7)
      ab(7)=ab(7)*ab(21)
      ab(4)= - 86*l5 - 4*ab(4) - 43*l7 - ab(29) + ab(62) + 44*l10 + 18*
     & l3
      ab(4)=ab(4)*ab(6)
      ab(6)= - l8 + ab(14) + 11
      ab(6)= - ab(58) - ab(37) + 2*ab(6)
      ab(6)=ab(6)*ab(17)
      ab(8)=l11 + ab(28) - ab(63)
      ab(3)=ab(8)*ab(3)
      ab(8)=ab(31)*ab(41)
      ab(9)=ab(34) + d10
      ab(12)= - l11 + 2*ab(45)
      ab(12)=ab(12)*l11
      ab(3)= - ab(6) + ab(5) - ab(7) - ab(4) - ab(3) + ab(8) - 9*ab(38)
     &  - ab(25) - ab(76) + ab(12) + 2*ab(9)
      ab(4)=15*l6 - ab(57) + ab(49) + 27*l10 + 13*l3 - 62*ab(50)
      ab(5)= - 23 - ab(52)
      ab(5)=2*ab(5) + ab(4)
      ab(5)=ab(5)*ab(61)
      ab(6)= - 11*l10 - 4*l3 + ab(57) + ab(51)
      ab(7)= - 3 + l8
      ab(7)= - 17*l6 + 2*ab(7) + ab(6)
      ab(7)=ab(7)*ab(54)
      ab(8)=ab(22) - d4
      ab(9)= - 12 + ab(8)
      ab(5)=4*ab(7) + ab(5) + 2*ab(9) + ab(3)
      ab(5)=ab(5)*ab(54)
      ab(4)=14*ab(41) - ab(4)
      ab(4)=ab(4)*ab(61)
      ab(7)=36 - ab(8)
      ab(3)=ab(5) + ab(4) + 2*ab(7) - ab(3)
      ab(3)=ab(3)*ab(54)
      ab(4)=9 - l8
      ab(4)=ab(55) + 2*ab(4) - ab(6)
      ab(2)=ab(2) + 4*ab(4) + ab(3)
      ab(2)=xb*ab(2)
      ab(2)=ab(13) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=ab(66) + 3*ab(72)
      ab(3)=ab(3)*ab(61)
      ab(4)=ab(79) - 18
      ab(4)=ab(77) + 2*ab(4)
      ab(3)=ab(3) + ab(4) - ab(44)
      ab(5)=ab(39) + ab(47)
      ab(6)=ab(5) - ab(15)
      ab(7)=ab(10) - 2
      ab(7)= - ab(55) + 2*ab(7)
      ab(7)=ab(7)*l6
      ab(6)=ab(7) + ab(11) + 2*ab(6)
      ab(6)=ab(6)*ab(54)
      ab(5)=ab(23) - ab(5)
      ab(7)=1 - ab(10)
      ab(7)=2*ab(7) + ab(55)
      ab(7)=l6*ab(7)
      ab(5)= - ab(6) + ab(7) + 2*ab(5) + ab(16)
      ab(5)=ab(5)*ab(56)
      ab(7)= - 4 + ab(67)
      ab(7)=ab(7)*ab(61)
      ab(1)=ab(5) + ab(7) + ab(4) - ab(1)
      ab(1)=ab(1)*ab(54)
      ab(1)=ab(2) + ab(1) + ab(3)
      ab(1)=xb*ab(1)
      ab(2)=6*ab(24) + ab(6)
      ab(2)=ab(2)*ab(56)
      ab(2)=ab(2) - ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

