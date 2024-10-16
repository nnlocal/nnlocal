      subroutine AGTBa1a1crssscarggg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(0:2,0:2,-4:0)
      complex(ki) ab(58)
      include 'AGTBa1a1crssscar_functions.h'

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(1,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xb - xa

      tmp = 4*CA**2*l1*z6*z7*z9*xa**3*ab(1)
      res(1,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l2 + l4 - l3
      ab(2)= - xb*ab(1)
      ab(2)= - 1 + ab(2)
      ab(3)=l1*xb
      ab(2)=2*ab(2) + ab(3)
      ab(3)=4*l1
      ab(2)=ab(2)*ab(3)
      ab(1)=2*ab(1) - l1
      ab(1)=ab(1)*ab(3)
      ab(3)=Pi**2
      ab(3)=1.0_ki/3.0_ki*ab(3)
      ab(1)=ab(3) + ab(1)
      ab(1)=xa*ab(1)
      ab(3)= - xb*ab(3)
      ab(1)=ab(1) + ab(3) + ab(2)

      tmp = CA**2*z6*z7*z9*xa**3*ab(1)
      res(1,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)= - l3 + 2*l2
      ab(1)=ab(1)*l3
      ab(2)=l3 - l2
      ab(3)= - l4 + 2*ab(2)
      ab(3)=ab(3)*l4
      ab(4)=l2**2
      ab(1)= - ab(4) + ab(1) + ab(3)
      ab(3)= - xb*ab(1)
      ab(4)=l4 + l2
      ab(5)= - 2*l3 + l5 + ab(4)
      ab(2)=1.0_ki/3.0_ki*l1 + ab(2) - l4
      ab(6)=xb*ab(2)
      ab(6)= - 1 + ab(6)
      ab(6)=l1*ab(6)
      ab(3)=ab(6) + 2*ab(5) + ab(3)
      ab(5)=4*l1
      ab(3)=ab(3)*ab(5)
      ab(2)= - l1*ab(2)
      ab(1)=ab(2) + ab(1)
      ab(1)=ab(1)*ab(5)
      ab(2)=ab(4) - l3
      ab(2)=1.0_ki/3.0_ki*ab(2)
      ab(4)=Pi**2
      ab(2)=ab(2)*ab(4)
      ab(2)=ab(2) + 2*zeta3
      ab(1)=ab(1) - ab(2)
      ab(1)=xa*ab(1)
      ab(2)=xb*ab(2)
      ab(1)=ab(1) + ab(3) + 1.0_ki/3.0_ki*ab(4) + ab(2)

      tmp = 2*CA**2*z6*z7*z9*xa**3*ab(1)
      res(1,0,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 4*CA**2*z1*z2*xb**3
      res(0,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - l6 - l8 + l7
      ab(1)=2*ab(1) + l1

      tmp = 4*CA**2*z1*z2*xb**3*ab(1)
      res(0,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - l1 + l6
      ab(1)=l6*ab(1)
      ab(2)= - l1 + 2*l6
      ab(3)=l7 - ab(2)
      ab(3)=l7*ab(3)
      ab(2)=l8 - 2*l7 + ab(2)
      ab(2)=l8*ab(2)
      ab(1)=ab(2) + ab(1) + ab(3)
      ab(2)=Pi**2
      ab(3)=l1**2
      ab(1)=1.0_ki/3.0_ki*ab(2) + 4*ab(3) + 8*ab(1)

      tmp = CA**2*z1*z2*xb**3*ab(1)
      res(0,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=Pi**2
      ab(2)=l1**2
      ab(1)=1.0_ki/3.0_ki*ab(1) + 4*ab(2)
      ab(2)=l6 - l1
      ab(3)=8*l6
      ab(2)=ab(2)*ab(3)
      ab(2)=ab(1) + ab(2)
      ab(3)= - l1 + 2*l6
      ab(4)=2.0_ki/3.0_ki*l7 - ab(3)
      ab(4)=l7*ab(4)
      ab(4)=4*ab(4) + ab(2)
      ab(4)=l7*ab(4)
      ab(5)= - l7 + ab(3)
      ab(5)=l7*ab(5)
      ab(3)= - 2.0_ki/3.0_ki*l8 + 2*l7 - ab(3)
      ab(3)=l8*ab(3)
      ab(2)=4*ab(3) + 8*ab(5) - ab(2)
      ab(2)=l8*ab(2)
      ab(3)=l1 - 2.0_ki/3.0_ki*l6
      ab(3)=l6*ab(3)
      ab(1)=4*ab(3) - ab(1)
      ab(1)=l6*ab(1)
      ab(3)=l1**3
      ab(3)= - 3*zeta3 + 2.0_ki/3.0_ki*ab(3)
      ab(1)=ab(2) + ab(4) + 2*ab(3) + ab(1)

      tmp = 2*CA**2*z1*z2*xb**3*ab(1)
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
      ab(9)=l7 + l6
      ab(10)=4*ab(9)
      ab(8)=ab(8) + ab(10)
      ab(11)=2*l1
      ab(1)= - ab(1)*ab(11)
      ab(1)=ab(1) + ab(7) + ab(8)
      ab(1)=xa*ab(1)
      ab(12)= - 2 + l1
      ab(2)=ab(12)*ab(2)
      ab(12)= - 1 - 4*l2 - 6*ab(9)
      ab(1)=ab(1) + ab(3) + ab(2) + 2*ab(12) - ab(5)
      ab(1)=xa*ab(1)
      ab(2)=1 + ab(10)
      ab(3)=4*l1
      ab(10)= - ab(3) - ab(6) + ab(2)
      ab(12)=3*l1
      ab(13)= - 1 + ab(12)
      ab(13)=ab(13)*ab(11)
      ab(8)=ab(13) - ab(5) - ab(8)
      ab(8)=xb*ab(8)
      ab(8)=2*ab(10) + ab(8)
      ab(8)=xb*ab(8)
      ab(2)=ab(2) + ab(6)
      ab(10)=l1 - 4
      ab(3)=ab(10)*ab(3)
      ab(2)= - 2*ab(2) + ab(3) - 5.0_ki/3.0_ki*ab(4)
      ab(1)=ab(1) + ab(8) - ab(2)
      ab(1)=xa*ab(1)
      ab(3)= - 3 + l1
      ab(3)=ab(3)*ab(11)
      ab(3)=ab(3) - ab(6) - ab(7)
      ab(3)=xb*ab(3)
      ab(4)=10 - ab(12)
      ab(4)=ab(4)*ab(11)
      ab(6)=6*l2 + ab(9)
      ab(3)=ab(3) + ab(4) + 4*ab(6) + ab(5)
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
      ab(4)=xa*ab(3)
      ab(5)=3*ab(3)
      ab(6)=xb*ab(5)
      ab(6)=ab(6) - 10*ab(3) + ab(4)
      ab(6)=xb*ab(6)
      ab(7)=3*ab(4)
      ab(8)=4*ab(3)
      ab(9)=ab(8) - ab(7)
      ab(9)=xa*ab(9)
      ab(10)=8*ab(3)
      ab(6)=ab(6) + ab(10) + ab(9)
      ab(6)=xb*ab(6)
      ab(9)=ab(4) - ab(3)
      ab(11)=ab(9)*xb
      ab(12)=1.0_ki/3.0_ki*xb
      ab(13)=ab(12) - 1
      ab(13)=ab(13)*ab(11)
      ab(14)= - ab(3) + 1.0_ki/3.0_ki*ab(4)
      ab(15)=xa**2
      ab(16)=ab(14)*ab(15)
      ab(17)=2.0_ki/3.0_ki*ab(3)
      ab(13)=ab(13) - ab(16) - ab(17)
      ab(13)=ab(13)*xb
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + ab(17)
      ab(14)=ab(14)*xa
      ab(13)=ab(13) + ab(14)
      ab(14)=2*l1
      ab(16)=ab(13)*ab(14)
      ab(17)=6*ab(3) - ab(4)
      ab(17)=xa*ab(17)
      ab(10)= - ab(10) + ab(17)
      ab(10)=xa*ab(10)
      ab(6)=ab(16) + ab(10) + ab(6)
      ab(6)=l1*ab(6)
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(1) - ab(4)
      ab(2)=ab(2)*xa
      ab(10)=ab(3)*xb
      ab(16)= - ab(1) + ab(10)
      ab(16)=xb*ab(16)
      ab(16)=ab(16) + ab(1) - ab(2)
      ab(6)=2*ab(16) + ab(6)
      ab(6)=l1*ab(6)
      ab(16)=ab(5) - ab(4)
      ab(17)=ab(16)*xa
      ab(17)=ab(17) - ab(1)
      ab(17)=ab(17)*xa
      ab(18)=xb - 2
      ab(18)=ab(18)*ab(11)
      ab(18)=ab(17) + ab(18)
      ab(19)=2*l6
      ab(14)=ab(19) + ab(14)
      ab(14)=ab(18)*ab(14)
      ab(20)=7*ab(3)
      ab(7)= - ab(20) + ab(7)
      ab(7)=xa*ab(7)
      ab(7)=ab(11) + ab(8) + ab(7) + ab(14)
      ab(7)=l6*ab(7)
      ab(14)=l7 + l1 + ab(19)
      ab(14)=l7*ab(18)*ab(14)
      ab(5)=ab(10) - ab(5)
      ab(5)=ab(5)*xb
      ab(5)=ab(5) + ab(2) + ab(1)
      ab(5)=ab(5)*xb
      ab(2)=ab(5) - ab(2)
      ab(5)=8*l2 + 10*l1
      ab(5)=ab(2)*ab(5)
      ab(18)=xb - xa
      ab(18)=ab(18)*ab(9)
      ab(5)= - 3*ab(18) + ab(5)
      ab(5)=l2*ab(5)
      ab(5)=ab(5) + 2*ab(14) + ab(6) + ab(7)
      ab(6)=1.0_ki/2.0_ki*ab(4)
      ab(7)= - 5.0_ki/2.0_ki*ab(10) + ab(20) + ab(6)
      ab(7)=ab(7)*ab(12)
      ab(4)= - ab(1) + 5.0_ki/6.0_ki*ab(4)
      ab(4)=xa*ab(4)
      ab(4)=ab(7) - 4.0_ki/3.0_ki*ab(3) + ab(4)
      ab(4)=xb*ab(4)
      ab(7)=l1*ab(13)
      ab(3)= - ab(3) - ab(6)
      ab(3)=xa*ab(3)
      ab(3)=ab(8) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)= - 5*ab(7) + 1.0_ki/3.0_ki*ab(3) + ab(4)
      ab(3)=ab(3)*Pi**2
      ab(2)=l3*l1*ab(2)
      ab(4)=l5*ab(18)
      ab(2)=ab(2) + ab(4)
      ab(4)=ab(16)*ab(15)
      ab(6)= - 3 + xb
      ab(6)=ab(6)*ab(11)
      ab(1)=ab(6) - ab(1) + ab(4)
      ab(1)=xb*ab(1)
      ab(1)= - ab(17) + ab(1)
      ab(1)=zeta3*ab(1)
      ab(4)= - xb + 1
      ab(4)=l9*ab(9)*ab(4)

      tmp = 18*ab(1) + 4*ab(2) + ab(3) + 8*ab(4) + 2*ab(5)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z6*z7*z8
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=ab(3)*xb
      ab(5)=14*ab(3)
      ab(6)=xa*ab(3)
      ab(7)=5*ab(4) - ab(5) - ab(6)
      ab(8)=1.0_ki/3.0_ki*xb
      ab(7)=ab(7)*ab(8)
      ab(9)=4*ab(3)
      ab(10)=ab(9) - 5.0_ki/3.0_ki*ab(6)
      ab(10)=xa*ab(10)
      ab(11)=8.0_ki/3.0_ki*ab(3)
      ab(7)=ab(7) + ab(11) + ab(10)
      ab(7)=xb*ab(7)
      ab(10)=ab(6) - ab(3)
      ab(12)=ab(8) - 1
      ab(12)=ab(12)*xb
      ab(13)=ab(10)*ab(12)
      ab(14)=1.0_ki/3.0_ki*ab(6)
      ab(15)=ab(14) - ab(3)
      ab(16)=xa**2
      ab(17)=ab(15)*ab(16)
      ab(13)=ab(13) - ab(17)
      ab(17)=2.0_ki/3.0_ki*ab(3)
      ab(18)= - ab(17) + ab(13)
      ab(18)=ab(18)*xb
      ab(15)=ab(15)*xa
      ab(17)=ab(15) + ab(17)
      ab(17)=ab(17)*xa
      ab(18)=ab(18) + ab(17)
      ab(19)=5*l1
      ab(20)=ab(18)*ab(19)
      ab(2)=2*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(1) + ab(6)
      ab(2)=xa*ab(2)
      ab(21)=8*ab(3)
      ab(2)= - ab(21) + ab(2)
      ab(22)=1.0_ki/3.0_ki*xa
      ab(2)=ab(2)*ab(22)
      ab(2)=ab(20) + ab(2) + ab(7)
      ab(2)=l1*ab(2)
      ab(7)=ab(10)*xb
      ab(20)= - ab(7) + 2*ab(10)
      ab(8)=ab(20)*ab(8)
      ab(8)=ab(17) + ab(8)
      ab(17)=ab(8)*l6
      ab(5)= - ab(5) + 11*ab(6)
      ab(5)=xa*ab(5)
      ab(22)=5*ab(3) - 4*ab(6)
      ab(22)=2*ab(22) - ab(4)
      ab(22)=xb*ab(22)
      ab(5)=ab(22) + ab(1) + ab(5)
      ab(22)=1.0_ki/3.0_ki*ab(3)
      ab(13)= - ab(22) + 1.0_ki/2.0_ki*ab(13)
      ab(13)=xb*ab(13)
      ab(15)=ab(22) + 1.0_ki/2.0_ki*ab(15)
      ab(15)=xa*ab(15)
      ab(13)=ab(15) + ab(13)
      ab(15)=Pi**2
      ab(13)=ab(13)*ab(15)
      ab(2)=113.0_ki/60.0_ki*ab(13) + ab(17) + 1.0_ki/3.0_ki*ab(5) + ab(2)
      ab(2)=ab(2)*ab(15)
      ab(5)=3*ab(6)
      ab(13)= - ab(5) + 7*ab(3)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) - ab(7) - ab(9)
      ab(20)=ab(20)*xb
      ab(22)=3*ab(3)
      ab(23)=ab(22) - ab(6)
      ab(24)=ab(23)*xa
      ab(24)=ab(24) - ab(1)
      ab(25)=ab(24)*xa
      ab(20)=ab(25) - ab(20)
      ab(26)=4*l1
      ab(27)= - ab(20)*ab(26)
      ab(17)=8*ab(17) + 3*ab(13) + ab(27)
      ab(17)=l6*ab(17)
      ab(27)= - l1*ab(20)
      ab(27)=ab(27) + ab(13)
      ab(27)=ab(27)*ab(26)
      ab(17)=ab(27) + ab(17)
      ab(17)=l6*ab(17)
      ab(27)=ab(4) - ab(22)
      ab(27)=ab(27)*xb
      ab(28)=ab(1) - ab(6)
      ab(28)=ab(28)*xa
      ab(29)=ab(28) + ab(1)
      ab(27)=ab(27) + ab(29)
      ab(27)=ab(27)*xb
      ab(27)=ab(27) - ab(28)
      ab(30)=l1*ab(27)
      ab(12)=ab(3)*ab(12)
      ab(12)=ab(12) + 1.0_ki/3.0_ki*ab(29)
      ab(12)=ab(12)*xb
      ab(12)=ab(12) - 1.0_ki/3.0_ki*ab(28)
      ab(28)=l2*ab(12)
      ab(29)= - ab(22) - ab(6)
      ab(29)=xa*ab(29)
      ab(28)= - 64*ab(28) - 36*ab(30) + 5*ab(7) + ab(9) + ab(29)
      ab(28)=l2*ab(28)
      ab(29)= - xb*ab(1)
      ab(29)=ab(29) + 11*ab(3) - 7*ab(6)
      ab(29)=xb*ab(29)
      ab(19)= - ab(27)*ab(19)
      ab(30)=13*ab(3) - ab(5)
      ab(30)=xa*ab(30)
      ab(19)=ab(19) + ab(29) - 12*ab(3) + ab(30)
      ab(29)=2*l1
      ab(19)=ab(19)*ab(29)
      ab(30)=25*ab(3) - 9*ab(6)
      ab(30)=xa*ab(30)
      ab(30)= - 7*ab(7) - 16*ab(3) + ab(30)
      ab(30)=l6*ab(30)
      ab(19)=ab(19) + ab(30)
      ab(12)=ab(12)*ab(15)
      ab(19)=ab(28) + 2*ab(19) + 11*ab(12)
      ab(19)=l2*ab(19)
      ab(28)=ab(10)*xa
      ab(28)=ab(28) - ab(7)
      ab(30)=l5 - l2
      ab(28)=ab(28)*ab(30)
      ab(30)= - ab(9) + ab(6)
      ab(30)=xa*ab(30)
      ab(31)=2*ab(6)
      ab(32)=ab(31) - ab(4)
      ab(32)=xb*ab(32)
      ab(30)=ab(32) + ab(1) + ab(30)
      ab(30)=ab(30)*ab(29)
      ab(13)= - l6*ab(13)
      ab(13)=ab(30) + ab(13) + ab(28)
      ab(13)=l5*ab(13)
      ab(28)= - 4.0_ki/3.0_ki*ab(3) + ab(6)
      ab(28)=xa*ab(28)
      ab(30)= - ab(6) + 10*ab(3)
      ab(32)=1.0_ki/3.0_ki*ab(30) - ab(4)
      ab(32)=xb*ab(32)
      ab(28)=ab(32) - ab(11) + ab(28)
      ab(28)=xb*ab(28)
      ab(14)= - ab(1) + ab(14)
      ab(14)=xa*ab(14)
      ab(11)=ab(11) + ab(14)
      ab(11)=xa*ab(11)
      ab(11)=ab(11) + ab(28)
      ab(14)= - l1*ab(18)
      ab(11)=2*ab(11) + ab(14)
      ab(11)=l1*ab(11)
      ab(3)=6*ab(3)
      ab(6)=ab(3) - ab(6)
      ab(6)=ab(6)*xa
      ab(14)=ab(22) - ab(31)
      ab(4)=2*ab(14) - ab(4)
      ab(4)=xb*ab(4)
      ab(3)=ab(4) - ab(3) + ab(6)
      ab(3)=2*ab(3) + ab(11)
      ab(4)=l1**2
      ab(3)=ab(3)*ab(4)
      ab(11)=2*l6
      ab(14)= - l1 - ab(11)
      ab(14)=ab(20)*ab(14)
      ab(18)=l7*ab(8)
      ab(14)=2*ab(18) + ab(14)
      ab(14)=l7*ab(14)
      ab(18)= - l1 - l6
      ab(11)=ab(11)*ab(18)
      ab(11)= - ab(4) + ab(11)
      ab(11)=ab(20)*ab(11)
      ab(11)=ab(11) + ab(14)
      ab(8)=ab(8)*ab(15)
      ab(8)=ab(8) + 4*ab(11)
      ab(8)=l7*ab(8)
      ab(11)=l2*l1
      ab(14)= - l3*ab(26)
      ab(4)=ab(14) - 8*ab(11) - 4*ab(4)
      ab(4)=ab(27)*ab(4)
      ab(4)= - ab(12) + ab(4)
      ab(4)=l3*ab(4)
      ab(5)= - ab(9) + ab(5)
      ab(5)=xa*ab(5)
      ab(9)= - xb*ab(22)
      ab(9)=ab(9) + ab(30)
      ab(9)=xb*ab(9)
      ab(5)=ab(9) - ab(21) + ab(5)
      ab(5)=xb*ab(5)
      ab(9)=3*ab(10) - ab(7)
      ab(9)=xb*ab(9)
      ab(11)= - ab(23)*ab(16)
      ab(1)=ab(9) + ab(1) + ab(11)
      ab(1)=xb*ab(1)
      ab(1)=ab(25) + ab(1)
      ab(1)=l1*ab(1)
      ab(6)=ab(21) - ab(6)
      ab(6)=xa*ab(6)
      ab(1)=6*ab(1) + ab(6) + ab(5)
      ab(1)=zeta3*ab(1)
      ab(5)=ab(7) - ab(10)
      ab(6)=l9 - 2*l5 + 3*l2 + ab(29) + l6
      ab(6)=l9*ab(6)
      ab(6)=8*ab(6) + 8*d1 - 32*d2
      ab(5)=ab(5)*ab(6)
      ab(6)=ab(7) - ab(24)
      ab(6)=d3*ab(6)

      tmp = 6*ab(1) + ab(2) + 2*ab(3) + ab(4) + ab(5) + 16*ab(6) + 
     & ab(8) + 4*ab(13) + ab(17) + ab(19)
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

      tmp = 2*CA**2*z6*z10*ab(1)
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
      ab(2)= - 3*l1 - 4*l2 + l4 + l3

      tmp = 2*ab(2)*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=3*l1
      ab(2)=ab(1) + 2*l2
      ab(3)=8*l2
      ab(2)=ab(2)*ab(3)
      ab(1)=ab(1) + 4*l2
      ab(3)=ab(1) - l3
      ab(3)= - l4 + 2*ab(3)
      ab(3)=ab(3)*l4
      ab(1)= - l3 + 2*ab(1)
      ab(1)=ab(1)*l3
      ab(4)=Pi**2
      ab(5)=l1**2
      ab(1)=ab(3) + ab(1) - ab(2) + ab(4) - 9*ab(5)
      ab(2)=2*ab(1)
      ab(3)= - xa*ab(1)
      ab(3)=ab(2) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)= - 3*ab(1) + ab(3)
      ab(3)=xa*ab(3)
      ab(2)=ab(2) + ab(3)
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
      ab(5)=ab(3)*l3
      ab(6)=1.0_ki/3.0_ki*xa
      ab(2)=ab(2)*ab(6)
      ab(2)=ab(2) + ab(1)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - 2.0_ki/3.0_ki*ab(1)
      ab(2)=ab(2)*xa
      ab(1)=ab(2) + 1.0_ki/3.0_ki*ab(1)
      ab(2)=4*l2
      ab(6)= - ab(1)*ab(2)
      ab(6)=ab(6) + ab(4) + ab(5)
      ab(6)=ab(6)*ab(2)
      ab(7)=2*l4
      ab(7)=ab(3)*ab(7)
      ab(5)=ab(7) + ab(5)
      ab(5)=ab(5)*l3
      ab(7)=ab(3)*l4**2
      ab(5)=ab(5) + ab(7)
      ab(6)=ab(6) - ab(5)
      ab(6)=ab(6)*ab(2)
      ab(8)=l3 + l4
      ab(2)=ab(8) - ab(2)
      ab(9)= - l1 + ab(2)
      ab(10)=3*l1
      ab(9)=ab(10)*ab(9)
      ab(8)= - 2*l2 + ab(8)
      ab(8)=l2*ab(8)
      ab(8)=ab(9) + 8*ab(8)
      ab(8)=ab(3)*ab(8)
      ab(5)= - ab(5) + ab(8)
      ab(5)=ab(5)*ab(10)
      ab(2)=ab(10) - ab(2)
      ab(2)=Pi**2*ab(2)
      ab(2)= - 16*zeta3 + ab(2)
      ab(2)=ab(3)*ab(2)
      ab(3)=l3*ab(1)
      ab(3)=ab(4) + ab(3)
      ab(3)=l3*ab(3)
      ab(3)=ab(7) + ab(3)
      ab(3)=l3*ab(3)
      ab(1)=ab(1)*l4**3

      tmp = ab(1) + ab(2) + ab(3) + ab(5) + ab(6)
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

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l8 + l6
      ab(2)=1 + 2*l11
      ab(3)=l4 + l2
      ab(3)=3*ab(3)
      ab(4)= - ab(3) + l7 + 2*l5
      ab(5)=5*l3
      ab(2)= - 2*ab(2) + ab(1) + ab(4) + ab(5)
      ab(6)=xb**2
      ab(7)=2*ab(6)
      ab(8)=ab(7) - ab(2)
      ab(8)=ab(8)*ab(6)
      ab(2)=ab(8) + ab(2)
      ab(2)=ab(2)*ab(6)
      ab(4)= - l3 + 2 - 4*l10 - ab(4) + 3*ab(1)
      ab(3)=ab(3) - l7 + ab(1) - 2*l10
      ab(8)= - ab(5) + ab(3)
      ab(9)= - ab(8)*ab(6)
      ab(9)=ab(9) - ab(4)
      ab(9)=ab(9)*ab(6)
      ab(1)= - l5 + ab(1) - l10
      ab(10)=ab(1) + 2*l3
      ab(11)=1 - ab(10)
      ab(9)=2*ab(11) + ab(9)
      ab(9)=xb*ab(9)
      ab(11)= - 2 + ab(3)
      ab(5)=ab(11) - ab(5)
      ab(5)=ab(5)*ab(6)
      ab(8)=ab(5) + ab(8)
      ab(8)=ab(8)*ab(6)
      ab(12)=1 + ab(10)
      ab(8)=2*ab(12) + ab(8)
      ab(8)=ab(8)*ab(6)
      ab(5)= - 2 - ab(5)
      ab(5)=ab(5)*ab(6)
      ab(10)=2*ab(10)
      ab(5)= - ab(10) + ab(5)
      ab(5)=xa*xb*ab(5)
      ab(5)=ab(5) + ab(10) + ab(8)
      ab(5)=xa*ab(5)
      ab(5)=ab(9) + ab(5)
      ab(5)=xa*ab(5)
      ab(2)=ab(5) - 2 + ab(2)
      ab(2)=xa*ab(2)
      ab(5)= - 1 + ab(1)
      ab(5)=ab(5)*ab(7)
      ab(4)=ab(5) + ab(4)
      ab(4)=ab(4)*ab(6)
      ab(3)= - l3 + ab(3)
      ab(4)=ab(4) + ab(3)
      ab(4)=xb*ab(4)
      ab(2)=ab(4) + ab(2)
      ab(2)=xa*ab(2)
      ab(4)=ab(1)*ab(6)
      ab(4)=ab(4) + 1
      ab(1)= - ab(1) - ab(4)
      ab(1)=ab(1)*ab(7)
      ab(1)=ab(1) - ab(3)
      ab(1)=ab(1)*ab(6)
      ab(3)=ab(11) - l3
      ab(1)=ab(2) + ab(1) - ab(3)
      ab(1)=xa*ab(1)
      ab(2)=ab(4)*ab(7)
      ab(2)=ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l5 + 2
      ab(2)=3*l3
      ab(1)= - ab(2) - 12*l10 - 3*l7 + 2*ab(1)
      ab(3)=l8 + l2
      ab(4)=9*l4
      ab(3)=ab(4) + ab(1) + 9*ab(3)
      ab(3)=l6 + 2*ab(3)
      ab(3)=ab(3)*l6
      ab(5)=2*l2
      ab(6)=ab(5) + 2*l4
      ab(7)=2*l7
      ab(8)=2*l10
      ab(9)= - ab(6) - 3 + ab(7) + ab(8)
      ab(10)=4*l10
      ab(9)=ab(9)*ab(10)
      ab(11)=3*l4
      ab(12)=3*l2
      ab(13)=ab(11) + ab(12)
      ab(14)=l5 - 1
      ab(15)=2*ab(14)
      ab(16)=ab(13) - ab(15)
      ab(16)= - l7 + 2*ab(16)
      ab(16)=ab(16)*l7
      ab(17)=l4 + l2
      ab(1)=ab(1) + ab(17)
      ab(1)=17*l8 + 2*ab(1)
      ab(1)=ab(1)*l8
      ab(18)=ab(15) - ab(12)
      ab(18)= - ab(11) + 2*ab(18)
      ab(18)=ab(18)*ab(11)
      ab(19)= - ab(12) + 4*ab(14)
      ab(19)=ab(19)*ab(12)
      ab(20)=4*l5
      ab(21)=l5 + 1
      ab(22)=ab(20)*ab(21)
      ab(23)=d13 - d12
      ab(24)=d10 - d6
      ab(25)= - 16*ab(24)
      ab(26)=Pi**2
      ab(27)=3*ab(26)
      ab(28)=d5 + d4
      ab(29)=8*ab(28)
      ab(1)= - ab(19) + ab(3) - ab(16) + ab(9) - ab(22) - ab(25) - 
     & ab(27) - ab(29) + ab(1) - ab(18) + 16*ab(23)
      ab(3)=2*l8
      ab(9)=3*l10
      ab(16)=2*l6
      ab(18)= - l5 + ab(16) + ab(3) - ab(9)
      ab(19)=l7 + 1
      ab(22)= - ab(2) - ab(19) + ab(13) + ab(18)
      ab(30)=xa**2
      ab(31)=4*ab(30)
      ab(22)=ab(22)*ab(31)
      ab(32)=7*ab(17)
      ab(33)=ab(32) - l7
      ab(34)= - ab(33) + 8*l10
      ab(35)=3*l5
      ab(36)= - 5 + ab(35)
      ab(36)=2*ab(36) + ab(34)
      ab(37)=7*l3
      ab(36)=2*ab(36) - ab(37)
      ab(36)=l3*ab(36)
      ab(38)=d8 + d7
      ab(38)= - d9 + 2*ab(38)
      ab(39)=ab(38) - 3
      ab(22)=ab(22) + ab(36) + 4*ab(39) + ab(1)
      ab(22)=ab(22)*ab(30)
      ab(36)= - 1 - ab(35)
      ab(34)=2*ab(36) - ab(34)
      ab(34)=2*ab(34) + ab(37)
      ab(34)=l3*ab(34)
      ab(36)=1 - ab(38)
      ab(1)=ab(22) + ab(34) + 4*ab(36) - ab(1)
      ab(1)=ab(1)*ab(30)
      ab(22)=2*l11
      ab(34)=ab(22) - ab(35)
      ab(36)=l3 + l7
      ab(37)=ab(36) + ab(10)
      ab(34)=ab(37) + ab(17) + 2*ab(34)
      ab(34)= - 7*l8 + 2*ab(34)
      ab(34)=ab(34)*l8
      ab(22)=ab(22) - ab(14)
      ab(40)=2*ab(22)
      ab(41)=ab(13) + ab(40)
      ab(41)= - l7 + 2*ab(41)
      ab(41)=ab(41)*l7
      ab(42)=ab(6) + 1
      ab(43)=4*l11
      ab(44)= - ab(20) + ab(42) + ab(43)
      ab(44)=ab(44)*ab(10)
      ab(45)=ab(37) - 3*l8
      ab(46)=ab(17) + 2*l5
      ab(46)= - ab(45) + 3*ab(46)
      ab(46)= - 5*l6 + 2*ab(46)
      ab(46)=ab(46)*l6
      ab(40)=ab(40) + ab(12)
      ab(40)=ab(11) + 2*ab(40)
      ab(40)=ab(40)*ab(11)
      ab(22)=ab(12) + 4*ab(22)
      ab(22)=ab(22)*ab(12)
      ab(43)=ab(43) - ab(14)
      ab(43)=ab(43)*ab(20)
      ab(22)= - ab(46) + ab(34) - ab(41) + ab(25) - ab(29) + ab(26) + 
     & ab(22) + ab(43) + ab(40) - ab(44) - 8*ab(23)
      ab(25)=ab(45) - ab(13)
      ab(34)=ab(25) - 2
      ab(34)= - l6 + 2*ab(34)
      ab(34)=ab(34)*l6
      ab(37)=ab(17) - ab(37)
      ab(40)= - 2 - ab(37)
      ab(41)=5*l8
      ab(40)= - ab(41) + 2*ab(40)
      ab(40)=ab(40)*l8
      ab(43)=l10 + l7
      ab(42)=ab(42) - ab(43)
      ab(42)=ab(42)*ab(10)
      ab(23)=ab(24) + ab(23)
      ab(24)= - 4*ab(23)
      ab(34)=ab(40) + ab(34) + ab(42) + ab(24) + ab(20)
      ab(5)=ab(5) + l4
      ab(4)=ab(5)*ab(4)
      ab(5)=6*ab(17)
      ab(40)= - l7 + ab(5)
      ab(40)=ab(40)*l7
      ab(42)=l2**2
      ab(4)= - ab(4) + ab(40) - 9*ab(42)
      ab(27)= - ab(27) + 16*ab(28)
      ab(40)=ab(27) - ab(4) + 4*ab(38)
      ab(42)=5*l7
      ab(32)= - ab(32) + ab(42) + 10*l10
      ab(44)= - 4 - ab(32)
      ab(44)=2*ab(44) - l3
      ab(44)=l3*ab(44)
      ab(44)=ab(44) - ab(40) + ab(34)
      ab(44)=ab(44)*ab(30)
      ab(46)=d14 + d11
      ab(47)=l11**2
      ab(46)=ab(47) + 8*d15 - 2*ab(46)
      ab(47)=3 + ab(46)
      ab(42)= - 23*l4 - 19*l2 + ab(42) + ab(10)
      ab(48)=l5 + 3
      ab(48)=2*ab(48) - ab(42)
      ab(49)=9*l3
      ab(48)=2*ab(48) - ab(49)
      ab(48)=l3*ab(48)
      ab(44)=ab(44) + ab(48) + 4*ab(47) - ab(22)
      ab(44)=ab(44)*ab(30)
      ab(47)=ab(25) - ab(21)
      ab(48)= - ab(47)*ab(16)
      ab(5)= - ab(7) + ab(5)
      ab(50)= - ab(14)*ab(5)
      ab(51)=ab(26) - 6*ab(23)
      ab(52)=ab(43) - 2
      ab(53)=ab(52)*ab(8)
      ab(54)=2*l3
      ab(35)=ab(35) - ab(54)
      ab(9)=ab(9) + ab(35)
      ab(55)= - 1 + ab(9)
      ab(55)=ab(55)*ab(54)
      ab(56)=ab(21) - ab(45)
      ab(3)=ab(56)*ab(3)
      ab(57)=l5**2
      ab(58)=2*ab(57)
      ab(3)=ab(48) + ab(3) + ab(55) + ab(53) - ab(58) + 2*ab(39) - 
     & ab(51) + ab(50)
      ab(3)=2*ab(3) + ab(44)
      ab(3)=xa*ab(3)
      ab(39)=ab(25) - ab(14)
      ab(16)=ab(39)*ab(16)
      ab(44)=ab(45) - ab(14)
      ab(44)=ab(44)*l8
      ab(16)=ab(51) + ab(16) + 2*ab(44)
      ab(45)=ab(21)*ab(5)
      ab(48)=ab(43) + 2
      ab(50)= - ab(48)*ab(8)
      ab(9)= - 1 - ab(9)
      ab(9)=ab(9)*ab(54)
      ab(51)=ab(38) + 1
      ab(9)=ab(9) + ab(50) + ab(58) - 2*ab(51) + ab(45) + ab(16)
      ab(25)=ab(25) + 2
      ab(25)= - l6 + 2*ab(25)
      ab(25)=ab(25)*l6
      ab(37)=2 - ab(37)
      ab(37)= - ab(41) + 2*ab(37)
      ab(37)=ab(37)*l8
      ab(24)=ab(37) + ab(25) + ab(24)
      ab(25)=ab(13) - 2
      ab(25)= - l7 + 2*ab(25)
      ab(25)=ab(25)*l7
      ab(17)=ab(17) - 1
      ab(17)= - ab(43) + 2*ab(17)
      ab(17)=ab(17)*ab(10)
      ab(37)=ab(12) - 2
      ab(37)=ab(11) + 2*ab(37)
      ab(11)=ab(37)*ab(11)
      ab(37)=ab(12) - 4
      ab(12)=ab(37)*ab(12)
      ab(11)=ab(25) - ab(12) + ab(24) + ab(17) - ab(11)
      ab(12)=ab(32) + 10
      ab(12)=l3 + 2*ab(12)
      ab(12)=ab(12)*l3
      ab(12)=ab(12) + ab(27) - ab(11) + 4*ab(51)
      ab(12)=ab(12)*ab(30)
      ab(17)=ab(19) + l10
      ab(6)=ab(17) - ab(6)
      ab(6)=ab(6)*ab(10)
      ab(6)=ab(6) - ab(24) + ab(20)
      ab(10)=8 + ab(32)
      ab(10)=2*ab(10) + l3
      ab(10)=l3*ab(10)
      ab(10)=ab(12) + ab(10) + ab(40) + ab(6)
      ab(10)=ab(10)*ab(30)
      ab(9)=2*ab(9) + ab(10)
      ab(9)=ab(9)*ab(30)
      ab(10)=l10 + ab(36) - ab(13) - ab(14)
      ab(10)=4*ab(10) - ab(12)
      ab(10)=ab(10)*ab(30)
      ab(12)=l10 + l5
      ab(12)= - ab(54) + 3*ab(12)
      ab(12)=ab(12)*ab(54)
      ab(5)=ab(5) + ab(15)
      ab(5)=ab(5)*l5
      ab(15)=ab(17)*ab(8)
      ab(5)=ab(12) - ab(5) - ab(16) + ab(15) + 2*ab(38)
      ab(5)=2*ab(5)
      ab(10)=ab(5) + ab(10)
      ab(10)=xb*xa*ab(10)
      ab(5)=ab(10) - ab(5) + ab(9)
      ab(5)=xb*ab(5)
      ab(3)=ab(3) + ab(5)
      ab(3)=xb*ab(3)
      ab(5)=ab(13) - 3
      ab(9)=ab(36) - ab(5) - ab(18)
      ab(1)=ab(3) + 4*ab(9) + ab(1)
      ab(1)=xb*ab(1)
      ab(3)=l6*ab(47)
      ab(9)=ab(13) - l7
      ab(10)=ab(14)*ab(9)
      ab(7)= - l10 + ab(35) - ab(7)
      ab(12)=5 - ab(7)
      ab(12)=l3*ab(12)
      ab(13)= - l8*ab(56)
      ab(15)= - l10*ab(52)
      ab(3)=ab(3) + ab(13) + ab(12) + ab(15) + ab(57) + 1 + ab(10) + 6*
     & ab(28) - 3*ab(23)
      ab(3)=ab(3)*ab(31)
      ab(10)= - 1 - ab(46)
      ab(12)= - 2*ab(21) + ab(42)
      ab(12)=2*ab(12) + ab(49)
      ab(12)=l3*ab(12)
      ab(3)=ab(3) + ab(12) + 4*ab(10) + ab(22)
      ab(3)=ab(3)*ab(30)
      ab(10)=ab(29) + ab(26)
      ab(4)=ab(4) + ab(10)
      ab(8)=ab(33) - ab(8)
      ab(12)= - 2*ab(8) + l3
      ab(12)=l3*ab(12)
      ab(3)=ab(3) + ab(12) - ab(4) - ab(34)
      ab(3)=xa*ab(3)
      ab(1)=ab(3) + ab(1)
      ab(1)=xb*ab(1)
      ab(3)=ab(39)*l6
      ab(3)=ab(3) + ab(44)
      ab(12)= - ab(23) + 2*ab(28)
      ab(13)=ab(7) - 2
      ab(13)=ab(13)*l3
      ab(14)=ab(9) + ab(14)
      ab(14)=ab(14)*l5
      ab(15)=ab(17)*l10
      ab(13)= - 3*ab(12) + ab(15) - ab(3) + ab(13) - ab(14)
      ab(13)=ab(13)*ab(30)
      ab(9)= - ab(21)*ab(9)
      ab(12)=1 - ab(12)
      ab(7)=1 + ab(7)
      ab(7)=l3*ab(7)
      ab(14)=l10*ab(48)
      ab(3)=ab(13) + ab(7) + ab(14) + 3*ab(12) - ab(57) + ab(9) - ab(3)
      ab(3)=ab(3)*ab(31)
      ab(7)=4 + ab(8)
      ab(7)=2*ab(7) - l3
      ab(7)=l3*ab(7)
      ab(3)=ab(3) + ab(7) + ab(4) - ab(6)
      ab(3)=ab(3)*ab(30)
      ab(4)=ab(8) - 2
      ab(4)= - l3 + 2*ab(4)
      ab(4)=ab(4)*l3
      ab(4)=ab(4) + ab(10) + ab(11) - 12
      ab(1)=ab(1) + ab(3) + ab(4)
      ab(1)=xb*ab(1)
      ab(2)= - ab(13) - ab(2) + l5 - ab(43) + ab(5)
      ab(2)=ab(2)*ab(31)
      ab(2)=ab(2) - ab(4)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*z5*z7*z10*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

