      subroutine hp_BGTAa1a1crssscarggg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) tmp,hp_cli2,hp_Li3
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      complex(ki) ab(49)
      include 'hp_BGTAa1a1crssscar_functions.h'

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(1,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      

      tmp =  - 4*CA**2*l1*z6*z9*xa**3
      res(1,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=Pi**2
      ab(2)=l4 + l2 - l3
      ab(2)=2*ab(2) - l1
      ab(2)=l1*ab(2)
      ab(1)=1.0_ki/3.0_ki*ab(1) + 4*ab(2)

      tmp = CA**2*z6*z9*xa**3*ab(1)
      res(1,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l2**2
      ab(2)=2*l2 - l3
      ab(2)=l3*ab(2)
      ab(3)=l3 - l2
      ab(4)=2*ab(3) - l4
      ab(4)=l4*ab(4)
      ab(3)=ab(3) - l4
      ab(5)= - 1.0_ki/3.0_ki*l1 - ab(3)
      ab(5)=l1*ab(5)
      ab(1)=ab(5) + ab(4) - ab(1) + ab(2)
      ab(1)=l1*ab(1)
      ab(2)=1.0_ki/3.0_ki*ab(3)
      ab(2)=ab(2)*Pi**2
      ab(1)=4*ab(1) - 2*zeta3 + ab(2)

      tmp = 2*CA**2*z6*z9*xa**3*ab(1)
      res(1,0,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - xa + xb

      tmp = 4*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - xb + xa
      ab(2)=l5 + l7 - l6
      ab(2)= - l1 + 2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(1)=2 + ab(1)

      tmp = 4*CA**2*z1*z2*z4*xb**3*ab(1)
      res(0,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - l1 + 2*l5
      ab(2)= - ab(1) - l7 + 2*l6
      ab(3)=8*l7
      ab(2)=ab(2)*ab(3)
      ab(3)=ab(1) - l6
      ab(4)=8*l6
      ab(3)=ab(3)*ab(4)
      ab(4)=l5 - l1
      ab(5)=8*l5
      ab(4)=ab(4)*ab(5)
      ab(5)=Pi**2
      ab(6)=l1**2
      ab(2)= - ab(2) - ab(3) + ab(4) + 1.0_ki/3.0_ki*ab(5) + 4*ab(6)
      ab(3)=xb - xa
      ab(2)=ab(2)*ab(3)
      ab(3)=l7 + l8
      ab(1)=4*l6 - ab(1) - 2*ab(3)
      ab(1)=8*ab(1) + ab(2)

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
      ab(10)=2*ab(8) + ab(9)
      ab(10)=l5*ab(10)
      ab(11)=ab(8) + ab(9)
      ab(12)=1.0_ki/3.0_ki*ab(2)
      ab(12)=l1*ab(12)
      ab(12)=ab(12) - ab(11)
      ab(12)=l1*ab(12)
      ab(10)=ab(12) + ab(10) + ab(6)
      ab(10)=l1*ab(10)
      ab(12)=l7**2
      ab(13)=ab(12)*ab(5)
      ab(1)=ab(3)*ab(1)
      ab(3)=ab(5)*l6
      ab(14)=ab(1) - 1.0_ki/3.0_ki*ab(3)
      ab(14)=l6*ab(14)
      ab(14)= - ab(13) + ab(14)
      ab(14)=l6*ab(14)
      ab(15)= - ab(3) + 2*ab(1)
      ab(15)=ab(15)*l6
      ab(13)=ab(15) - ab(13)
      ab(15)=ab(3) - ab(1)
      ab(16)=ab(5)*l5
      ab(17)=1.0_ki/3.0_ki*ab(16) - ab(15)
      ab(17)=l5*ab(17)
      ab(17)=ab(17) - ab(13)
      ab(17)=l5*ab(17)
      ab(18)=l7**3*ab(5)
      ab(14)=ab(17) + 1.0_ki/3.0_ki*ab(18) + ab(14)
      ab(17)=2*ab(15) - ab(16)
      ab(17)=l5*ab(17)
      ab(15)=ab(15) - ab(16)
      ab(18)=ab(5)*l1
      ab(19)= - 1.0_ki/3.0_ki*ab(18) - ab(15)
      ab(19)=l1*ab(19)
      ab(13)=ab(19) + ab(17) + ab(13)
      ab(13)=l1*ab(13)
      ab(13)=2*ab(14) + ab(13)
      ab(13)=xa*ab(13)
      ab(8)= - 1.0_ki/3.0_ki*ab(9) - ab(8)
      ab(8)=l5*ab(8)
      ab(6)=ab(8) - ab(6)
      ab(6)=l5*ab(6)
      ab(8)= - ab(4) + 4*ab(5)
      ab(9)= - l7*ab(8)
      ab(7)=1.0_ki/3.0_ki*ab(7) + ab(8)
      ab(7)=l6*ab(7)
      ab(7)=ab(9) + ab(7)
      ab(7)=l6*ab(7)
      ab(4)=ab(5) - 1.0_ki/3.0_ki*ab(4)
      ab(4)=ab(4)*ab(12)
      ab(4)=ab(6) + ab(4) + ab(7)
      ab(4)=ab(13) + 2*ab(4) + ab(10)
      ab(6)= - xa*ab(15)
      ab(6)=ab(6) - ab(11)
      ab(6)=ab(6)*Pi**2
      ab(1)= - 2*ab(3) + ab(16) + ab(1)
      ab(3)=l8*ab(5)
      ab(1)=ab(3) + 2*ab(1) - ab(18)
      ab(1)=l8*ab(1)
      ab(3)=xa*ab(5)
      ab(2)= - ab(2) + ab(3)
      ab(2)=zeta3*ab(2)
      ab(1)=6*ab(2) + 8*ab(1) + 4*ab(4) + 1.0_ki/3.0_ki*ab(6)

      tmp = 2*ab(1)
      res(0,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      
      ab(1)=3 - xb
      ab(2)=xa - 1
      ab(1)=xb*ab(2)*ab(1)
      ab(2)=xa - 3
      ab(3)=ab(2)*xa**2
      ab(1)=ab(1) + 2 + ab(3)
      ab(1)=xb*ab(1)
      ab(2)= - xa*ab(2)
      ab(2)= - 2 + ab(2)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = CA**2*z1*z6*z7*z8*ab(1)
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

      tmp = 2*CA**2*z1*z6*z7*z8*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l1 - 1
      ab(2)=6*l1
      ab(3)= - ab(1)*ab(2)
      ab(4)=Pi**2
      ab(5)=5.0_ki/2.0_ki*ab(4)
      ab(6)=8*l2
      ab(3)=ab(3) + ab(6) + ab(5)
      ab(3)=xb*ab(3)
      ab(7)=5.0_ki/6.0_ki*ab(4)
      ab(8)=l1**2
      ab(8)= - ab(7) + 2*ab(8)
      ab(8)=ab(8)*xb
      ab(9)=l6 + l5
      ab(10)=4*ab(9)
      ab(8)=ab(8) + ab(10)
      ab(11)=2*l1
      ab(1)= - ab(1)*ab(11)
      ab(1)=ab(1) + ab(7) + ab(8)
      ab(1)=xa*ab(1)
      ab(12)= - 2 + l1
      ab(2)=ab(12)*ab(2)
      ab(12)=1 - 2*l2 - 3*ab(9)
      ab(1)=ab(1) + ab(3) + ab(2) + 4*ab(12) - ab(5)
      ab(1)=xa*ab(1)
      ab(2)=4*l1
      ab(3)= - ab(2) + 1 - ab(6) + ab(10)
      ab(12)=3*l1
      ab(13)= - 1 + ab(12)
      ab(13)=ab(13)*ab(11)
      ab(8)=ab(13) - ab(5) - ab(8)
      ab(8)=xb*ab(8)
      ab(3)=2*ab(3) + ab(8)
      ab(3)=xb*ab(3)
      ab(8)=ab(10) + ab(6) - 5
      ab(10)=l1 - 4
      ab(2)=ab(10)*ab(2)
      ab(2)= - 2*ab(8) + ab(2) - 5.0_ki/3.0_ki*ab(4)
      ab(1)=ab(1) + ab(3) - ab(2)
      ab(1)=xa*ab(1)
      ab(3)= - 3 + l1
      ab(3)=ab(3)*ab(11)
      ab(3)=ab(3) - ab(6) - ab(7)
      ab(3)=xb*ab(3)
      ab(4)= - 1 + 4*l2
      ab(4)=3*ab(4) + 2*ab(9)
      ab(6)=10 - ab(12)
      ab(6)=ab(6)*ab(11)
      ab(3)=ab(3) + ab(6) + 2*ab(4) + ab(5)
      ab(3)=xb*ab(3)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = CA**2*z1*z6*z7*z8*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=3*ab(3)
      ab(5)=ab(4)*xb
      ab(6)=xa*ab(3)
      ab(7)=ab(5) - 10*ab(3) + ab(6)
      ab(7)=xb*ab(7)
      ab(8)=4*ab(3)
      ab(9)=ab(8) - 3*ab(6)
      ab(9)=xa*ab(9)
      ab(10)=8*ab(3)
      ab(7)=ab(7) + ab(10) + ab(9)
      ab(7)=xb*ab(7)
      ab(9)=ab(6) - ab(3)
      ab(11)=ab(9)*xb
      ab(12)=1.0_ki/3.0_ki*xb
      ab(13)=ab(12) - 1
      ab(13)=ab(13)*ab(11)
      ab(14)= - ab(3) + 1.0_ki/3.0_ki*ab(6)
      ab(15)=xa**2
      ab(16)=ab(14)*ab(15)
      ab(17)=2.0_ki/3.0_ki*ab(3)
      ab(13)= - ab(13) + ab(16) + ab(17)
      ab(13)=ab(13)*xb
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + ab(17)
      ab(14)=ab(14)*xa
      ab(13)=ab(13) - ab(14)
      ab(14)=2*l1
      ab(16)= - ab(13)*ab(14)
      ab(17)=6*ab(3) - ab(6)
      ab(17)=xa*ab(17)
      ab(10)= - ab(10) + ab(17)
      ab(10)=xa*ab(10)
      ab(7)=ab(16) + ab(10) + ab(7)
      ab(7)=l1*ab(7)
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(1) - ab(6)
      ab(2)=ab(2)*xa
      ab(10)=ab(2) + ab(1)
      ab(16)=ab(3)*xb
      ab(17)= - ab(1) + ab(16)
      ab(18)=3*xb
      ab(17)=ab(17)*ab(18)
      ab(17)=ab(17) + ab(10)
      ab(7)=2*ab(17) + ab(7)
      ab(7)=l1*ab(7)
      ab(17)= - ab(11) + 2*ab(9)
      ab(17)=ab(17)*xb
      ab(18)=ab(4) - ab(6)
      ab(19)=ab(18)*xa
      ab(19)=ab(19) - ab(1)
      ab(19)=ab(19)*xa
      ab(17)=ab(19) - ab(17)
      ab(20)=2*l5
      ab(14)=ab(20) + ab(14)
      ab(14)=ab(17)*ab(14)
      ab(21)=9*ab(16) - 17*ab(3) - ab(6)
      ab(21)=xb*ab(21)
      ab(22)=9*ab(3) - 4*ab(6)
      ab(22)=xa*ab(22)
      ab(14)=ab(21) + ab(8) + ab(22) + ab(14)
      ab(14)=l5*ab(14)
      ab(4)=ab(16) - ab(4)
      ab(4)=ab(4)*xb
      ab(4)=ab(4) + ab(10)
      ab(4)=ab(4)*xb
      ab(2)=ab(4) - ab(2)
      ab(4)=8*l2 + 10*l1
      ab(4)=ab(2)*ab(4)
      ab(10)=5*ab(3)
      ab(5)=ab(5) - ab(10) - ab(6)
      ab(5)=xb*ab(5)
      ab(10)=ab(10) - 2*ab(6)
      ab(10)=xa*ab(10)
      ab(4)=ab(10) + ab(5) + ab(4)
      ab(4)=l2*ab(4)
      ab(5)=l6 + l1 + ab(20)
      ab(5)=l6*ab(17)*ab(5)
      ab(4)=ab(4) + 2*ab(5) + ab(7) + ab(14)
      ab(2)=l3*l1*ab(2)
      ab(5)= - xa*ab(9)
      ab(5)=ab(5) + ab(11)
      ab(5)=l8*ab(5)
      ab(2)=ab(2) + ab(5)
      ab(5)=1.0_ki/2.0_ki*ab(6)
      ab(7)= - 5.0_ki/2.0_ki*ab(16) + 7*ab(3) + ab(5)
      ab(7)=ab(7)*ab(12)
      ab(6)= - ab(1) + 5.0_ki/6.0_ki*ab(6)
      ab(6)=xa*ab(6)
      ab(6)=ab(7) - 4.0_ki/3.0_ki*ab(3) + ab(6)
      ab(6)=xb*ab(6)
      ab(7)=l1*ab(13)
      ab(3)= - ab(3) - ab(5)
      ab(3)=xa*ab(3)
      ab(3)=ab(8) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)=5*ab(7) + 1.0_ki/3.0_ki*ab(3) + ab(6)
      ab(3)=ab(3)*Pi**2
      ab(5)= - 3*ab(9) + ab(11)
      ab(5)=xb*ab(5)
      ab(6)=ab(18)*ab(15)
      ab(1)=ab(5) - ab(1) + ab(6)
      ab(1)=xb*ab(1)
      ab(1)= - ab(19) + ab(1)
      ab(1)=zeta3*ab(1)
      ab(5)= - ab(11) + ab(9)
      ab(5)=l9*ab(5)

      tmp = 18*ab(1) + 4*ab(2) + ab(3) + 2*ab(4) + 8*ab(5)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=ab(3)*xb
      ab(5)=xa*ab(3)
      ab(6)=5*ab(4) - 14*ab(3) - ab(5)
      ab(7)=1.0_ki/3.0_ki*xb
      ab(6)=ab(6)*ab(7)
      ab(8)=4*ab(3)
      ab(9)=ab(8) - 5.0_ki/3.0_ki*ab(5)
      ab(9)=xa*ab(9)
      ab(10)=8.0_ki/3.0_ki*ab(3)
      ab(6)=ab(6) + ab(10) + ab(9)
      ab(6)=xb*ab(6)
      ab(9)=ab(5) - ab(3)
      ab(11)=ab(7) - 1
      ab(11)=ab(11)*xb
      ab(12)=ab(9)*ab(11)
      ab(13)=1.0_ki/3.0_ki*ab(5)
      ab(14)=ab(13) - ab(3)
      ab(15)=xa**2
      ab(16)=ab(15)*ab(14)
      ab(12)=ab(12) - ab(16)
      ab(16)=2.0_ki/3.0_ki*ab(3)
      ab(17)= - ab(16) + ab(12)
      ab(17)=ab(17)*xb
      ab(14)=ab(14)*xa
      ab(16)=ab(14) + ab(16)
      ab(16)=ab(16)*xa
      ab(17)=ab(17) + ab(16)
      ab(18)=5*l1
      ab(18)=ab(17)*ab(18)
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(1) + ab(5)
      ab(19)=ab(2)*xa
      ab(20)=8*ab(3)
      ab(21)= - ab(20) + ab(19)
      ab(22)=1.0_ki/3.0_ki*xa
      ab(21)=ab(21)*ab(22)
      ab(6)=ab(18) + ab(21) + ab(6)
      ab(6)=l1*ab(6)
      ab(18)=ab(9)*xb
      ab(21)= - ab(18) + 2*ab(9)
      ab(7)=ab(21)*ab(7)
      ab(7)=ab(16) + ab(7)
      ab(16)=ab(7)*l5
      ab(23)=1.0_ki/3.0_ki*ab(3)
      ab(12)= - ab(23) + 1.0_ki/2.0_ki*ab(12)
      ab(12)=xb*ab(12)
      ab(14)=ab(23) + 1.0_ki/2.0_ki*ab(14)
      ab(14)=xa*ab(14)
      ab(12)=ab(14) + ab(12)
      ab(14)=Pi**2
      ab(12)=ab(12)*ab(14)
      ab(2)=ab(2)*ab(22)
      ab(23)=2*ab(5)
      ab(24)= - 7*ab(3) - ab(23)
      ab(25)=3*ab(3)
      ab(26)=ab(25)*xb
      ab(24)=2.0_ki/3.0_ki*ab(24) + ab(26)
      ab(24)=xb*ab(24)
      ab(2)=113.0_ki/60.0_ki*ab(12) + ab(16) + ab(6) + ab(24) + ab(1) + 
     & ab(2)
      ab(2)=ab(2)*ab(14)
      ab(6)= - 27*ab(4) + 55*ab(3) - ab(5)
      ab(6)=xb*ab(6)
      ab(12)=ab(21)*xb
      ab(21)=ab(25) - ab(5)
      ab(24)=ab(21)*xa
      ab(24)=ab(24) - ab(1)
      ab(24)=ab(24)*xa
      ab(12)=ab(24) - ab(12)
      ab(27)=4*l1
      ab(28)= - ab(12)*ab(27)
      ab(29)= - 15*ab(3) + 8*ab(5)
      ab(29)=xa*ab(29)
      ab(6)=8*ab(16) + ab(28) + ab(6) - 20*ab(3) + ab(29)
      ab(6)=l5*ab(6)
      ab(16)=9*ab(4)
      ab(28)= - ab(16) + 19*ab(3) - ab(5)
      ab(28)=xb*ab(28)
      ab(29)= - l1*ab(12)
      ab(30)= - ab(25) + ab(23)
      ab(30)=xa*ab(30)
      ab(28)=ab(29) + ab(28) - ab(20) + ab(30)
      ab(28)=ab(28)*ab(27)
      ab(6)=ab(28) + ab(6)
      ab(6)=l5*ab(6)
      ab(28)=ab(4) - ab(25)
      ab(28)=ab(28)*xb
      ab(29)=ab(1) - ab(5)
      ab(30)=ab(29)*xa
      ab(31)=ab(30) + ab(1)
      ab(28)=ab(28) + ab(31)
      ab(28)=ab(28)*xb
      ab(28)=ab(28) - ab(30)
      ab(30)=ab(28)*l1
      ab(32)= - ab(26) + 13*ab(3) - 7*ab(5)
      ab(32)=xb*ab(32)
      ab(33)=ab(25) + ab(23)
      ab(33)=xa*ab(33)
      ab(32)= - 5*ab(30) + ab(32) - ab(20) + ab(33)
      ab(33)=2*l1
      ab(32)=ab(32)*ab(33)
      ab(16)= - ab(16) + 31*ab(3) - 13*ab(5)
      ab(16)=xb*ab(16)
      ab(34)=ab(3) + 6*ab(5)
      ab(34)=xa*ab(34)
      ab(16)=ab(16) - 16*ab(3) + ab(34)
      ab(16)=l5*ab(16)
      ab(16)=ab(32) + ab(16)
      ab(32)=5*ab(5)
      ab(25)= - ab(4) - ab(25) + ab(32)
      ab(25)=xb*ab(25)
      ab(25)= - 12*ab(30) + ab(25) + ab(8) - ab(32)
      ab(11)=ab(3)*ab(11)
      ab(11)=ab(11) + 1.0_ki/3.0_ki*ab(31)
      ab(11)=ab(11)*xb
      ab(22)=ab(22)*ab(29)
      ab(11)=ab(11) - ab(22)
      ab(22)=l2*ab(11)
      ab(22)=3*ab(25) - 64*ab(22)
      ab(22)=l2*ab(22)
      ab(11)=ab(11)*ab(14)
      ab(16)=ab(22) + 2*ab(16) + 11*ab(11)
      ab(16)=l2*ab(16)
      ab(22)=ab(9)*xa
      ab(22)=ab(22) - ab(18)
      ab(25)=l8 + l2
      ab(22)=ab(22)*ab(25)
      ab(25)=ab(18) - ab(9)
      ab(29)=ab(25)*ab(27)
      ab(30)=5*ab(3)
      ab(31)= - ab(30) + ab(5)
      ab(31)=xa*ab(31)
      ab(32)=3*xb
      ab(32)=ab(9)*ab(32)
      ab(31)=ab(32) + ab(8) + ab(31)
      ab(31)=l5*ab(31)
      ab(22)=ab(29) + ab(31) + ab(22)
      ab(22)=l8*ab(22)
      ab(29)= - 4.0_ki/3.0_ki*ab(3) + ab(5)
      ab(29)=xa*ab(29)
      ab(31)= - ab(5) + 10*ab(3)
      ab(4)=1.0_ki/3.0_ki*ab(31) - ab(4)
      ab(4)=xb*ab(4)
      ab(4)=ab(4) - ab(10) + ab(29)
      ab(4)=xb*ab(4)
      ab(13)= - ab(1) + ab(13)
      ab(13)=xa*ab(13)
      ab(10)=ab(10) + ab(13)
      ab(10)=xa*ab(10)
      ab(4)=ab(10) + ab(4)
      ab(10)= - l1*ab(17)
      ab(4)=2*ab(4) + ab(10)
      ab(4)=l1*ab(4)
      ab(10)=ab(30) - ab(23)
      ab(10)=2*ab(10) - ab(26)
      ab(10)=xb*ab(10)
      ab(3)=6*ab(3)
      ab(10)=ab(10) - ab(3) + ab(19)
      ab(4)=2*ab(10) + ab(4)
      ab(10)=l1**2
      ab(4)=ab(4)*ab(10)
      ab(13)=2*l5
      ab(17)= - l1 - ab(13)
      ab(17)=ab(12)*ab(17)
      ab(19)=l6*ab(7)
      ab(17)=2*ab(19) + ab(17)
      ab(17)=l6*ab(17)
      ab(19)= - l1 - l5
      ab(13)=ab(13)*ab(19)
      ab(13)= - ab(10) + ab(13)
      ab(12)=ab(12)*ab(13)
      ab(12)=ab(12) + ab(17)
      ab(7)=ab(7)*ab(14)
      ab(7)=ab(7) + 4*ab(12)
      ab(7)=l6*ab(7)
      ab(12)=l2*l1
      ab(13)= - l3*ab(27)
      ab(10)=ab(13) - 8*ab(12) - 4*ab(10)
      ab(10)=ab(28)*ab(10)
      ab(10)= - ab(11) + ab(10)
      ab(10)=l3*ab(10)
      ab(8)= - ab(8) + 3*ab(5)
      ab(8)=xa*ab(8)
      ab(11)= - ab(26) + ab(31)
      ab(11)=xb*ab(11)
      ab(8)=ab(11) - ab(20) + ab(8)
      ab(8)=xb*ab(8)
      ab(9)=3*ab(9) - ab(18)
      ab(9)=xb*ab(9)
      ab(11)= - ab(21)*ab(15)
      ab(1)=ab(9) + ab(1) + ab(11)
      ab(1)=xb*ab(1)
      ab(1)=ab(24) + ab(1)
      ab(1)=l1*ab(1)
      ab(3)= - ab(3) + ab(5)
      ab(3)=xa*ab(3)
      ab(3)=ab(20) + ab(3)
      ab(3)=xa*ab(3)
      ab(1)=6*ab(1) + ab(3) + ab(8)
      ab(1)=zeta3*ab(1)
      ab(3)=d2 - d3
      ab(5)=l9 - 2*l8 + l2 + ab(33) + 3*l5
      ab(5)=l9*ab(5)
      ab(3)=8*ab(5) + 24*d1 - 32*ab(3)
      ab(3)=ab(25)*ab(3)

      tmp = 6*ab(1) + ab(2) + ab(3) + 2*ab(4) + ab(6) + ab(7) + ab(10)
     &  + ab(16) + 4*ab(22)
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

      tmp = 2*CA**2*z1*z6*z10*ab(1)
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
      ab(3)= - 3*l1 - 4*l2 + l3 + l4
      ab(2)=ab(3)*ab(2)
      ab(2)= - 3 + ab(2)

      tmp = 2*ab(2)*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
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
      ab(2)=ab(1)*ab(2)
      ab(3)=ab(2)*l3
      ab(4)=ab(2)*l4
      ab(5)=ab(3) + ab(4)
      ab(6)=ab(2)*l2
      ab(7)=3*ab(1)
      ab(5)=8*ab(6) + ab(7) - 4*ab(5)
      ab(5)=l2*ab(5)
      ab(6)=ab(4) - ab(7)
      ab(7)=4*ab(2)
      ab(7)=l2*ab(7)
      ab(7)=ab(7) - ab(3) - ab(6)
      ab(8)=3*l1
      ab(9)=ab(2)*ab(8)
      ab(7)=2*ab(7) + ab(9)
      ab(7)=ab(7)*ab(8)
      ab(4)= - 6*ab(1) + ab(4)
      ab(4)=l4*ab(4)
      ab(3)=2*ab(6) + ab(3)
      ab(3)=l3*ab(3)
      ab(2)= - ab(2)*Pi**2
      ab(1)=l5*ab(1)

      tmp = 18*ab(1) + ab(2) + ab(3) + ab(4) + 2*ab(5) + ab(7)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=xb - 1
      ab(2)=4*l1
      ab(2)=ab(2)*ab(1)
      ab(2)=ab(2) + 1
      ab(3)=8*l2
      ab(3)=ab(3)*ab(1)
      ab(2)=ab(3) + 3*ab(2)
      ab(3)=2*l2
      ab(2)=ab(2)*ab(3)
      ab(3)=ab(1)*l1
      ab(4)=ab(3) + 2
      ab(5)=9*l1
      ab(4)=ab(4)*ab(5)
      ab(6)=Pi**2
      ab(7)=ab(6)*ab(1)
      ab(8)=18*l5 - ab(7)
      ab(2)=ab(2) + ab(4) + ab(8)
      ab(4)=4*l2
      ab(4)=ab(4)*ab(1)
      ab(9)=ab(3) + 1
      ab(4)=ab(4) + 3*ab(9)
      ab(9)=ab(1)*l3
      ab(10)= - ab(9) + 2*ab(4)
      ab(10)=ab(10)*l3
      ab(10)=ab(10) - ab(2)
      ab(11)=ab(4) - ab(9)
      ab(12)=1.0_ki/3.0_ki*ab(1)
      ab(12)=ab(12)*l4
      ab(12)=ab(12) - ab(11)
      ab(12)=ab(12)*l4
      ab(12)=ab(10) - ab(12)
      ab(12)=ab(12)*l4
      ab(13)= - ab(4) + 1.0_ki/3.0_ki*ab(9)
      ab(13)=ab(13)*l3
      ab(13)=ab(13) + ab(2)
      ab(13)=ab(13)*l3
      ab(14)=ab(3) + 3
      ab(15)=3*l1
      ab(14)=ab(14)*ab(15)
      ab(8)=ab(14) + ab(8)
      ab(8)=ab(8)*ab(15)
      ab(14)=16*zeta3
      ab(14)=ab(14)*ab(1)
      ab(15)=l5**2
      ab(6)=ab(8) - 3*ab(6) + 27*ab(15) + ab(14)
      ab(8)=1 + 2*ab(3)
      ab(5)=ab(8)*ab(5)
      ab(5)= - ab(5) - 9*l5 + 2*ab(7)
      ab(3)=1 + 16*ab(3)
      ab(7)=ab(1)*l2
      ab(8)=3*ab(3) + 64.0_ki/3.0_ki*ab(7)
      ab(8)=ab(8)*l2
      ab(8)=ab(8) - 2*ab(5)
      ab(8)=ab(8)*l2
      ab(8)= - ab(6) - ab(12) + ab(13) - ab(8)
      ab(12)=2*ab(8)
      ab(13)=xa*ab(8)
      ab(13)= - ab(12) + ab(13)
      ab(13)=xa*ab(13)
      ab(4)= - 3*ab(4) + ab(9)
      ab(4)=l3*ab(4)
      ab(2)=3*ab(2) + ab(4)
      ab(2)=l3*ab(2)
      ab(1)=l4*ab(1)
      ab(1)= - 3*ab(11) + ab(1)
      ab(1)=l4*ab(1)
      ab(1)= - 3*ab(10) + ab(1)
      ab(1)=l4*ab(1)
      ab(3)= - 9*ab(3) - 64*ab(7)
      ab(3)=l2*ab(3)
      ab(3)=6*ab(5) + ab(3)
      ab(3)=l2*ab(3)
      ab(1)=ab(13) + ab(1) + ab(2) - 3*ab(6) + ab(3)
      ab(1)=xa*ab(1)
      ab(1)= - ab(12) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + ab(8)

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
      ab(2)=ab(1) - 1
      ab(2)=ab(2)*ab(1)
      ab(3)=2 - ab(2)
      ab(3)=xa*ab(3)
      ab(4)=ab(1) + 1
      ab(5)=ab(4)*ab(1)
      ab(5)= - 2 + ab(5)
      ab(5)=ab(5)*ab(1)
      ab(6)=xa**4
      ab(7)=2 - ab(6)
      ab(7)=xb*xa*ab(7)
      ab(5)=ab(7) - 2 + ab(5)
      ab(5)=xb*ab(5)
      ab(3)=ab(3) + ab(5)
      ab(3)=xb*ab(3)
      ab(2)=3*ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(3)=2*ab(1)
      ab(5)= - 1 - ab(3)
      ab(5)=ab(5)*ab(1)
      ab(5)=1 + ab(5)
      ab(5)=xa*ab(5)
      ab(2)=ab(5) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=ab(4)*ab(3)
      ab(3)= - 1 + ab(3)
      ab(1)=ab(3)*ab(1)
      ab(1)=ab(2) - 1 + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=1 - 2*ab(6)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l7 + l5
      ab(2)= - l3 + 3*ab(1)
      ab(3)= - l6 + ab(2)
      ab(4)=l4 + l2
      ab(5)=ab(4) - ab(3) - 6*l8 + 4*l11
      ab(6)= - ab(4) + 2*l10
      ab(7)=ab(6) + 5*l6 - ab(2)
      ab(8)=2 - ab(7)
      ab(9)=xa**2
      ab(8)=ab(8)*ab(9)
      ab(8)=ab(8) + ab(5)
      ab(8)=ab(8)*ab(9)
      ab(2)=l10 - l8 + 3*l6 - ab(2)
      ab(10)= - 1 + ab(2)
      ab(8)=2*ab(10) + ab(8)
      ab(8)=xa*ab(8)
      ab(7)=ab(7) + 2
      ab(10)=ab(9) + 1
      ab(11)=ab(9)*ab(7)*ab(10)
      ab(2)=ab(2) + 1
      ab(2)=2*ab(2)
      ab(11)= - ab(2) + ab(11)
      ab(11)=ab(11)*ab(9)
      ab(12)=xa**4
      ab(7)= - ab(7)*ab(12)
      ab(7)=ab(2) + ab(7)
      ab(7)=xb*xa*ab(7)
      ab(2)=ab(7) - ab(2) + ab(11)
      ab(2)=xb*ab(2)
      ab(2)=ab(8) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)= - 3*l3 - 7*l6 + ab(4) - 4*l10 + 4 + 2*l8 + 9*ab(1)
      ab(4)= - 4*ab(9) - ab(1)
      ab(4)=ab(4)*ab(9)
      ab(1)=ab(4) + ab(1)
      ab(1)=ab(1)*ab(9)
      ab(1)=ab(2) + 4 + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=2*ab(9)
      ab(4)=l8 + ab(3) - l10
      ab(7)=1 + ab(4)
      ab(7)=ab(7)*ab(2)
      ab(5)=ab(7) - ab(5)
      ab(5)=ab(5)*ab(9)
      ab(3)=ab(3) - ab(6)
      ab(5)=ab(5) - 2 - ab(3)
      ab(5)=xa*ab(5)
      ab(1)=ab(5) + ab(1)
      ab(1)=xb*ab(1)
      ab(4)= - 1 + ab(4)
      ab(2)= - ab(2)*ab(4)*ab(10)
      ab(3)=2 - ab(3)
      ab(2)=ab(2) - ab(3)
      ab(2)=ab(2)*ab(9)
      ab(1)=ab(1) + ab(2) - ab(3)
      ab(1)=xb*ab(1)
      ab(2)=ab(4)*ab(12)
      ab(2)=2*ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=2*l11
      ab(2)=3*l8
      ab(3)=ab(1) - ab(2)
      ab(3)=2*ab(3)
      ab(4)=ab(3) + l2
      ab(5)=3*l3
      ab(6)=4*l10
      ab(7)=ab(5) + ab(6)
      ab(4)= - l4 + ab(7) + 3*ab(4)
      ab(8)=5*l6
      ab(9)= - ab(8) + ab(4)
      ab(10)=9*l5
      ab(9)= - ab(10) + 2*ab(9)
      ab(9)=ab(9)*l5
      ab(11)=9*l6
      ab(4)= - ab(10) - ab(11) + ab(4)
      ab(12)=9*l7
      ab(4)= - ab(12) + 2*ab(4)
      ab(4)=ab(4)*l7
      ab(13)=4*l11
      ab(14)=l4 + l2
      ab(15)=2*ab(14)
      ab(16)=4*l8
      ab(17)= - ab(16) + ab(15) + ab(13) + 1
      ab(17)=ab(17)*ab(6)
      ab(18)=l8 - 1
      ab(1)=ab(1) - ab(18)
      ab(19)=3*l2
      ab(1)=ab(19) + 2*ab(1)
      ab(1)=7*l4 + 2*ab(1)
      ab(1)=ab(1)*l4
      ab(20)=d5 + d4
      ab(21)= - d11 - d14 + 4*d15 - ab(20)
      ab(22)=l11**2
      ab(21)=ab(22) + d7 + 2*ab(21)
      ab(3)=ab(3) + ab(14)
      ab(3)=l3 + 2*ab(3)
      ab(3)=ab(3)*l3
      ab(22)=ab(2) + 1
      ab(13)=ab(22) - ab(13)
      ab(13)=ab(13)*ab(16)
      ab(23)= - ab(19) + 4*ab(18)
      ab(23)=ab(23)*l2
      ab(24)=d13 - d12
      ab(25)=d9 + d8
      ab(26)=16*ab(25)
      ab(27)=Pi**2
      ab(1)=ab(1) - 4*ab(21) - ab(3) - ab(13) - ab(23) - ab(26) - 
     & ab(17) + ab(27) + ab(9) + ab(4) - 8*ab(24)
      ab(3)=7*l6
      ab(4)=ab(7) + ab(3) - l4
      ab(9)=ab(4) - ab(10)
      ab(13)=l8 + 1
      ab(17)=3*ab(13)
      ab(21)=ab(17) - ab(9)
      ab(21)=ab(12) + 2*ab(21)
      ab(21)=ab(21)*l7
      ab(17)=ab(17) - ab(4)
      ab(17)=ab(10) + 2*ab(17)
      ab(17)=ab(17)*l5
      ab(28)=l10 - 2
      ab(29)=2*l10
      ab(28)=ab(28)*ab(29)
      ab(30)=ab(13) - l10
      ab(30)= - l3 + 2*ab(30)
      ab(30)=ab(30)*l3
      ab(31)=2*ab(18)
      ab(32)=ab(31) + l4
      ab(32)=ab(32)*l4
      ab(33)=ab(31) - l2
      ab(33)=ab(33)*l2
      ab(17)=ab(21) + ab(28) - ab(32) - ab(33) + ab(17) - ab(30)
      ab(21)= - ab(24) + d10 - d6
      ab(24)=4*ab(25)
      ab(28)=ab(24) + ab(21)
      ab(20)= - d7 + 2*ab(20)
      ab(30)=ab(20) + 1
      ab(32)=l8**2
      ab(33)= - ab(32) + ab(30) + ab(28)
      ab(34)=ab(2) - l3
      ab(35)=l10 + 5 + ab(34)
      ab(35)=2*ab(35) - ab(8)
      ab(35)=l6*ab(35)
      ab(33)=ab(35) + 2*ab(33) - ab(17)
      ab(35)=xa**2
      ab(36)=2*ab(35)
      ab(33)=ab(33)*ab(36)
      ab(37)=l3 + ab(14) - ab(6)
      ab(38)=2 + 5*l8
      ab(38)=2*ab(38) - ab(37)
      ab(38)=2*ab(38) - ab(11)
      ab(38)=l6*ab(38)
      ab(33)=ab(33) + ab(38) + ab(1)
      ab(33)=ab(33)*ab(35)
      ab(38)=l2 + 2
      ab(38)=3*ab(38)
      ab(39)=ab(38) - ab(9)
      ab(39)=ab(12) + 2*ab(39)
      ab(39)=ab(39)*l7
      ab(38)=ab(38) - ab(4)
      ab(38)=ab(10) + 2*ab(38)
      ab(38)=ab(38)*l5
      ab(40)=ab(29) - ab(14)
      ab(41)=ab(40) - 2
      ab(41)=l3 + 2*ab(41)
      ab(41)=ab(41)*l3
      ab(42)=l10 - 1
      ab(43)=ab(42) - ab(15)
      ab(43)=ab(43)*ab(6)
      ab(44)=5*l4
      ab(45)=ab(44) + 6*l2
      ab(45)=ab(45)*l4
      ab(46)=l2**2
      ab(45)=ab(45) + ab(46)
      ab(46)=4*ab(21)
      ab(38)=ab(38) + ab(39) + ab(43) + ab(41) + ab(16) + ab(45) + 
     & ab(46)
      ab(25)=ab(27) + 8*ab(25)
      ab(39)=l3 - 6 + ab(40)
      ab(39)=2*ab(39) + l6
      ab(39)=l6*ab(39)
      ab(33)=ab(33) + ab(39) - 12 - ab(25) + ab(38)
      ab(33)=xa*ab(33)
      ab(39)=3*ab(27)
      ab(26)=ab(39) - ab(26)
      ab(20)=ab(20) - 1
      ab(41)= - ab(14) + 5*l3 + 10*l10
      ab(43)=2 - ab(41)
      ab(43)=2*ab(43) - l6
      ab(43)=l6*ab(43)
      ab(38)=ab(43) - 4*ab(20) + ab(26) - ab(38)
      ab(38)=ab(38)*ab(35)
      ab(37)= - 10*l8 + ab(37)
      ab(11)=2*ab(37) + ab(11)
      ab(11)=l6*ab(11)
      ab(1)=ab(38) + ab(11) - ab(1)
      ab(1)=ab(1)*ab(35)
      ab(11)= - ab(27) + ab(24) - 2*ab(21)
      ab(21)=3*l10
      ab(24)=ab(21) + ab(5)
      ab(22)= - ab(22) + ab(24)
      ab(22)=2*ab(22) + ab(8)
      ab(22)=l6*ab(22)
      ab(37)=2*ab(32)
      ab(17)=ab(22) + ab(37) - 6 + ab(11) + ab(17)
      ab(1)=2*ab(17) + ab(1)
      ab(1)=xa*ab(1)
      ab(17)=3*ab(18)
      ab(22)=ab(9) - ab(17)
      ab(22)= - ab(12) + 2*ab(22)
      ab(22)=ab(22)*l7
      ab(17)=ab(17) - ab(4)
      ab(17)=ab(10) + 2*ab(17)
      ab(17)=ab(17)*l5
      ab(38)=ab(18) - l10
      ab(43)= - l3 + 2*ab(38)
      ab(43)=ab(43)*l3
      ab(17)=ab(22) + ab(43) - ab(17)
      ab(22)=l10 + 1
      ab(43)=ab(22)*ab(29)
      ab(47)=2*l8
      ab(48)=ab(47) + l4
      ab(48)=ab(48)*l4
      ab(49)=ab(47) - l2
      ab(49)=ab(49)*l2
      ab(43)=ab(17) - ab(43) + ab(48) + ab(49)
      ab(11)=ab(11) + 6
      ab(38)=ab(38) - l3
      ab(38)= - ab(8) + 6*ab(38)
      ab(38)=ab(38)*l6
      ab(47)=ab(47)*ab(18)
      ab(38)=ab(38) - ab(11) - ab(47) + ab(43)
      ab(38)=2*ab(38)
      ab(47)=l2 - 2
      ab(47)=3*ab(47)
      ab(9)=ab(9) - ab(47)
      ab(9)= - ab(12) + 2*ab(9)
      ab(9)=ab(9)*l7
      ab(4)=ab(47) - ab(4)
      ab(4)=ab(10) + 2*ab(4)
      ab(4)=ab(4)*l5
      ab(40)=ab(40) + 2
      ab(47)=l3 + 2*ab(40)
      ab(47)=ab(47)*l3
      ab(4)= - ab(9) + ab(46) + ab(4) + ab(47)
      ab(9)=ab(22) - ab(15)
      ab(9)=ab(9)*ab(6)
      ab(9)=ab(9) - ab(16) + ab(45) + ab(4)
      ab(22)=ab(19) - 2
      ab(22)=ab(44) + 2*ab(22)
      ab(22)=ab(22)*l4
      ab(44)=ab(14) - 1
      ab(44)= - l10 + 2*ab(44)
      ab(44)=ab(44)*ab(6)
      ab(45)=l2 - 4
      ab(45)=ab(45)*l2
      ab(4)=ab(45) + ab(4) + ab(22) - ab(44)
      ab(22)= - ab(26) + 4*ab(30)
      ab(26)=ab(41) + 10
      ab(26)=l6 + 2*ab(26)
      ab(26)=ab(26)*l6
      ab(26)=ab(26) + ab(22) + ab(4)
      ab(26)=ab(26)*ab(35)
      ab(30)=14 + ab(41)
      ab(30)=2*ab(30) + l6
      ab(30)=l6*ab(30)
      ab(22)=ab(26) + ab(30) + ab(22) + ab(9)
      ab(22)=ab(22)*ab(35)
      ab(30)=l10 + 2
      ab(30)=ab(30)*ab(29)
      ab(41)=2*ab(13)
      ab(44)=ab(41) + l4
      ab(44)=ab(44)*l4
      ab(41)=ab(41) - l2
      ab(41)=ab(41)*l2
      ab(17)= - ab(30) + ab(17) + ab(44) + ab(41)
      ab(24)= - 1 + ab(2) - ab(24)
      ab(24)=2*ab(24) - ab(8)
      ab(24)=l6*ab(24)
      ab(11)=ab(24) - ab(37) - ab(11) + ab(17)
      ab(11)=2*ab(11) + ab(22)
      ab(11)=ab(11)*ab(35)
      ab(22)= - l8 + ab(14) - l10
      ab(24)= - 2*l6 - ab(22)
      ab(24)=4*ab(24) - ab(26)
      ab(24)=ab(24)*ab(35)
      ab(24)= - ab(38) + ab(24)
      ab(24)=xb*xa*ab(24)
      ab(11)=ab(24) + ab(38) + ab(11)
      ab(11)=xb*ab(11)
      ab(1)=ab(1) + ab(11)
      ab(1)=xb*ab(1)
      ab(11)=l8 + 2
      ab(11)=ab(14) + 2*ab(11)
      ab(3)= - ab(3) + ab(11) - ab(7)
      ab(7)=ab(10) + ab(3)
      ab(7)=ab(12) + 2*ab(7)
      ab(12)=3*l7
      ab(7)=ab(7)*ab(12)
      ab(3)=ab(10) + 2*ab(3)
      ab(10)=3*l5
      ab(3)=ab(3)*ab(10)
      ab(10)=ab(11) - ab(6)
      ab(5)= - ab(5) + 2*ab(10)
      ab(5)=ab(5)*l3
      ab(10)=ab(31) - ab(19)
      ab(10)= - 3*l4 + 2*ab(10)
      ab(10)=ab(10)*l4
      ab(11)= - ab(29) + ab(15) + 3
      ab(6)=ab(11)*ab(6)
      ab(11)=ab(16)*ab(13)
      ab(3)=ab(23) + ab(6) - ab(11) - ab(7) - ab(3) + ab(5) + ab(10)
      ab(5)=l7 + l5
      ab(5)= - 2*l3 - ab(21) + ab(14) + l8 + 6*ab(5)
      ab(6)= - 6*l6 - 2 + ab(5)
      ab(6)=ab(6)*ab(35)
      ab(7)= - ab(14) + 7*l3 + 8*l10
      ab(10)= - 8 - ab(2)
      ab(10)=2*ab(10) + ab(7)
      ab(11)=11*l6
      ab(10)=2*ab(10) + ab(11)
      ab(10)=l6*ab(10)
      ab(6)=4*ab(6) + ab(10) - 8 - ab(39) - ab(3)
      ab(6)=ab(6)*ab(35)
      ab(2)=2 + ab(2)
      ab(2)=2*ab(2) - ab(7)
      ab(2)=2*ab(2) - ab(11)
      ab(2)=l6*ab(2)
      ab(7)=8 + ab(27)
      ab(2)=ab(6) + ab(2) + 3*ab(7) + ab(3)
      ab(2)=ab(2)*ab(35)
      ab(3)=4*l6 + 6 - ab(5)
      ab(1)=ab(1) + 4*ab(3) + ab(2)
      ab(1)=xb*ab(1)
      ab(1)=ab(33) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=ab(28) + ab(20)
      ab(3)=ab(18)*l8
      ab(3)=ab(3) - ab(2)
      ab(5)=ab(34) + ab(42)
      ab(5)= - ab(8) + 2*ab(5)
      ab(5)=ab(5)*l6
      ab(3)= - ab(5) - ab(43) + 2*ab(3)
      ab(3)=ab(3)*ab(35)
      ab(2)=ab(32) - ab(2)
      ab(2)=ab(3) - ab(5) + 2*ab(2) - ab(17)
      ab(2)=ab(2)*ab(36)
      ab(5)=ab(40) + l3
      ab(5)=l6 + 2*ab(5)
      ab(5)=ab(5)*l6
      ab(5)=ab(5) - ab(25) + 12
      ab(2)=ab(2) - ab(5) - ab(9)
      ab(2)=ab(2)*ab(35)
      ab(4)=ab(5) + ab(4)
      ab(1)=ab(1) + ab(2) - ab(4)
      ab(1)=xb*ab(1)
      ab(2)=2*ab(22) - ab(3)
      ab(2)=ab(2)*ab(36)
      ab(2)=ab(2) + ab(4)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

