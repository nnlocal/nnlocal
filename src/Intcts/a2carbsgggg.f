      subroutine a2carbsgggg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(2,2,-4:0)
      complex(ki) ab(18)
      include 'a2carbs_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)= - 4 + xa
      ab(1)=xa*ab(1)
      ab(2)= - 1 + xa
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + 3 + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=3 - xa
      ab(2)=xa*ab(2)
      ab(1)=ab(1) - 2 + ab(2)

      tmp = CA**2*z1*z3*z5*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=4 - xa
      ab(1)=xa*ab(1)
      ab(2)=2 - xa
      ab(2)=2*ab(2) - xb
      ab(2)=xb*ab(2)
      ab(1)=ab(2) - 4 + ab(1)

      tmp = 2*CA**2*z1*z3*z5*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=1 - l2
      ab(2)=l3 + l2
      ab(2)=xa*ab(2)
      ab(3)=xb*l2
      ab(1)=ab(3) + ab(2) + 3*ab(1) - l3
      ab(1)=xb*ab(1)
      ab(2)=xa - 3
      ab(2)=l3*ab(2)
      ab(3)=l2 - 3
      ab(2)= - ab(3) + ab(2)
      ab(2)=xa*ab(2)
      ab(3)=l3 + ab(3)
      ab(1)=ab(1) + 2*ab(3) + ab(2)

      tmp = 2*CA**2*z1*z3*z5*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l2**2
      ab(2)=l3**2
      ab(3)=Pi**2
      ab(4)=13*zeta3
      ab(5)=xa*ab(4)
      ab(5)=ab(5) - 52*zeta3 - ab(2) + 2.0_ki/3.0_ki*ab(3) - ab(1)
      ab(5)=xa*ab(5)
      ab(6)=xa - 1
      ab(6)=ab(4)*ab(6)
      ab(3)=1.0_ki/3.0_ki*ab(3)
      ab(1)=ab(3) - ab(1) + ab(6)
      ab(1)=xb*ab(1)
      ab(6)=ab(3) + l4
      ab(6)= - 39*zeta3 + 4*ab(6)
      ab(7)= - 2 + l2
      ab(7)=l2*ab(7)
      ab(8)= - 14 + l3
      ab(8)=l3*ab(8)
      ab(1)=ab(1) + ab(5) + ab(8) + 3*ab(7) - ab(6)
      ab(1)=xb*ab(1)
      ab(5)=2*l4 + ab(3)
      ab(7)=10 - l2
      ab(7)=l2*ab(7)
      ab(8)=10 - l3
      ab(8)=l3*ab(8)
      ab(5)= - ab(4) + ab(8) + 2*ab(5) + ab(7)
      ab(2)= - ab(4) + ab(3) - ab(2)
      ab(2)=xa*ab(2)
      ab(3)= - 14 + l2
      ab(3)=l2*ab(3)
      ab(4)= - 2 + l3
      ab(4)=l3*ab(4)
      ab(2)=ab(2) + 3*ab(4) + ab(3) - ab(6)
      ab(2)=xa*ab(2)
      ab(1)=ab(1) + 2*ab(5) + ab(2)

      tmp = CA**2*z1*z3*z5*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z5*CA**2*z1*z3
      ab(2)=xb*ab(1)
      ab(3)= - ab(1) + 1.0_ki/3.0_ki*ab(2)
      ab(4)=ab(3)*xb
      ab(5)=2.0_ki/3.0_ki*ab(1)
      ab(4)=ab(4) + ab(5)
      ab(6)=ab(2) - ab(1)
      ab(7)=1.0_ki/3.0_ki*xa
      ab(7)=ab(7)*ab(6)
      ab(8)=ab(7) + ab(4)
      ab(9)=l2*ab(8)
      ab(10)=26*ab(1)
      ab(11)=xa*ab(1)
      ab(9)=ab(9) + 23*ab(11) - ab(10) + 3*ab(2)
      ab(9)=l2*ab(9)
      ab(12)=Pi**2
      ab(8)= - ab(8)*ab(12)
      ab(13)=3*ab(11) - 10*ab(1) + 7*ab(2)
      ab(13)=l3*ab(13)
      ab(8)=ab(9) + ab(8) + 2*ab(13)
      ab(8)=l2*ab(8)
      ab(6)=ab(6)*xa
      ab(9)=ab(1) - 1.0_ki/4.0_ki*ab(2)
      ab(9)=xb*ab(9)
      ab(9)= - 1.0_ki/4.0_ki*ab(6) - 3.0_ki/4.0_ki*ab(1) + ab(9)
      ab(9)=xa*ab(9)
      ab(13)=3*ab(1)
      ab(14)=ab(13) - ab(2)
      ab(14)=ab(14)*xb
      ab(15)=ab(1) - 1.0_ki/2.0_ki*ab(14)
      ab(9)=1.0_ki/2.0_ki*ab(15) + ab(9)
      ab(9)=ab(9)*ab(12)
      ab(15)=1.0_ki/3.0_ki*ab(11)
      ab(5)=3.0_ki/5.0_ki*ab(9) + ab(15) + ab(5) - ab(2)
      ab(5)=ab(5)*ab(12)
      ab(9)=2*ab(1)
      ab(16)=ab(9) - ab(2)
      ab(3)=ab(3) + ab(15)
      ab(3)=ab(3)*xa
      ab(3)=ab(3) + 1.0_ki/3.0_ki*ab(16)
      ab(15)=l3*ab(3)
      ab(10)=ab(15) + 7*ab(11) - ab(10) + 19*ab(2)
      ab(10)=l3*ab(10)
      ab(3)= - ab(3)*ab(12)
      ab(3)=ab(3) + ab(10)
      ab(3)=l3*ab(3)
      ab(10)=4*ab(1)
      ab(15)=ab(10) - ab(2)
      ab(17)=ab(15)*xb
      ab(18)= - 2*ab(16) + ab(11)
      ab(18)=xa*ab(18)
      ab(10)=ab(18) + ab(10) - ab(17)
      ab(6)= - ab(17) + ab(6) + ab(13)
      ab(6)=ab(6)*xa
      ab(6)=ab(6) + ab(14) - ab(9)
      ab(9)=l1*ab(6)
      ab(9)=6*ab(10) + 35*ab(9)
      ab(9)=zeta3*ab(9)
      ab(10)=1.0_ki/3.0_ki*xb
      ab(10)=ab(10)*ab(15)
      ab(1)= - ab(10) + ab(7) + ab(1)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) - ab(4)
      ab(4)=l1**2
      ab(7)= - ab(12) + ab(4)
      ab(1)=ab(4)*ab(1)*ab(7)
      ab(4)=plg4half*ab(6)
      ab(6)=ab(11) - ab(16)
      ab(7)=l4 + l3 + l2
      ab(6)=l4*ab(6)*ab(7)
      ab(2)= - ab(2) + ab(11)
      ab(2)=d1*ab(2)

      tmp = 10*ab(1) + 8*ab(2) + ab(3) + 80*ab(4) + ab(5) + 4*ab(6) + 
     & ab(8) + 2*ab(9)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(1,2,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=xa - 1
      ab(2)=2*ab(1)
      ab(3)=xb*ab(1)
      ab(3)= - ab(2) + ab(3)
      ab(3)=xb*ab(3)
      ab(3)=3*ab(1) + ab(3)
      ab(3)=xb*ab(3)
      ab(2)= - ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z3*z4*ab(1)
      res(1,2,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - xa + 1
      ab(2)= - l6 - l5 + l3 + l1
      ab(1)=ab(2)*ab(1)
      ab(1)= - 3 + ab(1)
      ab(2)=2*ab(1)
      ab(3)=xb*ab(1)
      ab(3)= - ab(2) + ab(3)
      ab(3)=xb*ab(3)
      ab(3)=3*ab(1) + ab(3)
      ab(3)=xb*ab(3)
      ab(2)= - ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z3*z4*ab(1)
      res(1,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=CA**2*z1*z3
      ab(2)=2*z4
      ab(2)=ab(1)*ab(2)
      ab(3)=ab(1)*z4
      ab(4)=xb*ab(3)
      ab(4)=ab(4) - ab(2)
      ab(5)=ab(4)*xb
      ab(6)=3*z4
      ab(1)=ab(6)*ab(1)
      ab(1)=ab(5) + ab(1)
      ab(1)=ab(1)*xb
      ab(1)=ab(1) - ab(2)
      ab(1)=ab(1)*xb
      ab(1)=ab(1) + ab(3)
      ab(2)=xa - 1
      ab(5)=ab(2)*ab(1)
      ab(6)=ab(5)*l6
      ab(7)= - 6*ab(1) + ab(6)
      ab(7)=l6*ab(7)
      ab(6)=ab(6) - 3*ab(1)
      ab(8)=ab(5)*l5
      ab(9)=2*ab(6) + ab(8)
      ab(9)=l5*ab(9)
      ab(6)=ab(6) + ab(8)
      ab(8)=ab(5)*l3
      ab(10)= - 2*ab(6) + ab(8)
      ab(10)=l3*ab(10)
      ab(6)=ab(8) - ab(6)
      ab(5)=l1*ab(5)
      ab(5)=2*ab(6) + ab(5)
      ab(5)=l1*ab(5)
      ab(1)=l2*ab(1)
      ab(6)=1.0_ki/3.0_ki*xb
      ab(4)=ab(4)*ab(6)
      ab(4)=ab(4) + ab(3)
      ab(4)=ab(4)*xb
      ab(4)=ab(4) - 2.0_ki/3.0_ki*ab(3)
      ab(4)=ab(4)*xb
      ab(3)=ab(4) + 1.0_ki/3.0_ki*ab(3)
      ab(2)= - Pi**2*ab(3)*ab(2)

      tmp = 18*ab(1) + ab(2) + ab(5) + ab(7) + ab(9) + ab(10)
      res(1,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=xa - 1
      ab(2)=ab(1)*l1
      ab(3)=ab(2) + 3
      ab(4)=ab(1)*l3
      ab(5)=ab(4) + 2*ab(3)
      ab(5)=ab(5)*l3
      ab(6)=1.0_ki/3.0_ki - 1.0_ki/3.0_ki*xa
      ab(7)=Pi**2
      ab(6)=ab(7)*ab(6)
      ab(6)=18*l2 + ab(6)
      ab(8)=ab(2) + 6
      ab(8)=ab(8)*l1
      ab(9)=ab(8) + ab(6)
      ab(10)=ab(5) + ab(9)
      ab(11)= - l1 - l3 + l5 + 1.0_ki/3.0_ki*l6
      ab(11)=ab(11)*ab(1)
      ab(11)=ab(11) - 3
      ab(11)=ab(11)*l6
      ab(12)=ab(3) + ab(4)
      ab(13)=ab(1)*l5
      ab(14)= - ab(13) + 2*ab(12)
      ab(14)=ab(14)*l5
      ab(11)=ab(11) - ab(14) + ab(10)
      ab(11)=ab(11)*l6
      ab(15)=1.0_ki/3.0_ki*ab(1)
      ab(15)=ab(15)*l5
      ab(15)=ab(12) - ab(15)
      ab(15)=ab(15)*l5
      ab(10)=ab(10) - ab(15)
      ab(10)=ab(10)*l5
      ab(15)=ab(3) + 1.0_ki/3.0_ki*ab(4)
      ab(15)=ab(15)*l3
      ab(9)=ab(15) + ab(9)
      ab(9)=ab(9)*l3
      ab(15)=l2**2
      ab(15)= - 4*zeta3 + 9*ab(15)
      ab(16)=xa*zeta3
      ab(15)= - ab(7) + 12*ab(16) + 3*ab(15)
      ab(16)=3 + 1.0_ki/3.0_ki*ab(2)
      ab(16)=ab(16)*l1
      ab(6)=ab(16) + ab(6)
      ab(6)=ab(6)*l1
      ab(6)=ab(11) - ab(15) + ab(10) - ab(6) - ab(9)
      ab(9)=2*ab(6)
      ab(10)=xb*ab(6)
      ab(10)= - ab(9) + ab(10)
      ab(10)=xb*ab(10)
      ab(7)= - ab(7)*ab(1)
      ab(7)=54*l2 + ab(7)
      ab(8)=ab(7) + 3*ab(8)
      ab(5)=ab(8) + 3*ab(5)
      ab(11)=ab(13) - ab(12)
      ab(1)=l6*ab(1)
      ab(1)=3*ab(11) + ab(1)
      ab(1)=l6*ab(1)
      ab(1)=ab(1) - 3*ab(14) + ab(5)
      ab(1)=l6*ab(1)
      ab(3)= - 3*ab(3) - ab(4)
      ab(3)=l3*ab(3)
      ab(3)=ab(3) - ab(8)
      ab(3)=l3*ab(3)
      ab(4)= - 3*ab(12) + ab(13)
      ab(4)=l5*ab(4)
      ab(4)=ab(4) + ab(5)
      ab(4)=l5*ab(4)
      ab(2)= - 9 - ab(2)
      ab(2)=l1*ab(2)
      ab(2)=ab(2) - ab(7)
      ab(2)=l1*ab(2)
      ab(1)=ab(10) + ab(1) + ab(4) + ab(3) - 3*ab(15) + ab(2)
      ab(1)=xb*ab(1)
      ab(1)= - ab(9) + ab(1)
      ab(1)=xb*ab(1)
      ab(1)=ab(1) + ab(6)

      tmp = CA**2*z1*z3*z4*ab(1)
      res(1,2,0) = real(tmp,ki)

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

      tmp = 2*CA**2*z1*z2*z3*ab(1)
      res(2,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - xb + 1
      ab(2)= - l8 - l7 + l2 + l1
      ab(1)=ab(2)*ab(1)
      ab(1)= - 3 + ab(1)
      ab(2)=2*ab(1)
      ab(3)=xa*ab(1)
      ab(3)= - ab(2) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)=3*ab(1) + ab(3)
      ab(3)=xa*ab(3)
      ab(2)= - ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA**2*z1*z2*z3*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=CA**2*z1*z2
      ab(2)=2*z3
      ab(2)=ab(1)*ab(2)
      ab(3)=ab(1)*z3
      ab(4)=xa*ab(3)
      ab(4)=ab(4) - ab(2)
      ab(5)=ab(4)*xa
      ab(6)=3*z3
      ab(1)=ab(6)*ab(1)
      ab(1)=ab(5) + ab(1)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) - ab(2)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) + ab(3)
      ab(2)=xb - 1
      ab(5)=ab(2)*ab(1)
      ab(6)=ab(5)*l8
      ab(7)= - 6*ab(1) + ab(6)
      ab(7)=l8*ab(7)
      ab(6)=ab(6) - 3*ab(1)
      ab(8)=ab(5)*l7
      ab(9)=2*ab(6) + ab(8)
      ab(9)=l7*ab(9)
      ab(6)=ab(6) + ab(8)
      ab(8)=ab(5)*l2
      ab(10)= - 2*ab(6) + ab(8)
      ab(10)=l2*ab(10)
      ab(6)=ab(8) - ab(6)
      ab(5)=l1*ab(5)
      ab(5)=2*ab(6) + ab(5)
      ab(5)=l1*ab(5)
      ab(1)=l3*ab(1)
      ab(6)=1.0_ki/3.0_ki*xa
      ab(4)=ab(4)*ab(6)
      ab(4)=ab(4) + ab(3)
      ab(4)=ab(4)*xa
      ab(4)=ab(4) - 2.0_ki/3.0_ki*ab(3)
      ab(4)=ab(4)*xa
      ab(3)=ab(4) + 1.0_ki/3.0_ki*ab(3)
      ab(2)= - Pi**2*ab(3)*ab(2)

      tmp = 18*ab(1) + ab(2) + ab(5) + ab(7) + ab(9) + ab(10)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=xb - 1
      ab(2)=ab(1)*l1
      ab(3)=ab(2) + 3
      ab(4)=ab(1)*l2
      ab(5)=ab(4) + 2*ab(3)
      ab(5)=ab(5)*l2
      ab(6)=1.0_ki/3.0_ki - 1.0_ki/3.0_ki*xb
      ab(7)=Pi**2
      ab(6)=ab(7)*ab(6)
      ab(6)=18*l3 + ab(6)
      ab(8)=ab(2) + 6
      ab(8)=ab(8)*l1
      ab(9)=ab(8) + ab(6)
      ab(10)=ab(5) + ab(9)
      ab(11)= - l1 - l2 + l7 + 1.0_ki/3.0_ki*l8
      ab(11)=ab(11)*ab(1)
      ab(11)=ab(11) - 3
      ab(11)=ab(11)*l8
      ab(12)=ab(3) + ab(4)
      ab(13)=ab(1)*l7
      ab(14)= - ab(13) + 2*ab(12)
      ab(14)=ab(14)*l7
      ab(11)=ab(11) - ab(14) + ab(10)
      ab(11)=ab(11)*l8
      ab(15)=1.0_ki/3.0_ki*ab(1)
      ab(15)=ab(15)*l7
      ab(15)=ab(12) - ab(15)
      ab(15)=ab(15)*l7
      ab(10)=ab(10) - ab(15)
      ab(10)=ab(10)*l7
      ab(15)=ab(3) + 1.0_ki/3.0_ki*ab(4)
      ab(15)=ab(15)*l2
      ab(9)=ab(15) + ab(9)
      ab(9)=ab(9)*l2
      ab(15)=l3**2
      ab(15)= - 4*zeta3 + 9*ab(15)
      ab(16)=xb*zeta3
      ab(15)= - ab(7) + 12*ab(16) + 3*ab(15)
      ab(16)=3 + 1.0_ki/3.0_ki*ab(2)
      ab(16)=ab(16)*l1
      ab(6)=ab(16) + ab(6)
      ab(6)=ab(6)*l1
      ab(6)=ab(11) - ab(15) + ab(10) - ab(6) - ab(9)
      ab(9)=2*ab(6)
      ab(10)=xa*ab(6)
      ab(10)= - ab(9) + ab(10)
      ab(10)=xa*ab(10)
      ab(7)= - ab(7)*ab(1)
      ab(7)=54*l3 + ab(7)
      ab(8)=ab(7) + 3*ab(8)
      ab(5)=ab(8) + 3*ab(5)
      ab(11)=ab(13) - ab(12)
      ab(1)=l8*ab(1)
      ab(1)=3*ab(11) + ab(1)
      ab(1)=l8*ab(1)
      ab(1)=ab(1) - 3*ab(14) + ab(5)
      ab(1)=l8*ab(1)
      ab(3)= - 3*ab(3) - ab(4)
      ab(3)=l2*ab(3)
      ab(3)=ab(3) - ab(8)
      ab(3)=l2*ab(3)
      ab(4)= - 3*ab(12) + ab(13)
      ab(4)=l7*ab(4)
      ab(4)=ab(4) + ab(5)
      ab(4)=l7*ab(4)
      ab(2)= - 9 - ab(2)
      ab(2)=l1*ab(2)
      ab(2)=ab(2) - ab(7)
      ab(2)=l1*ab(2)
      ab(1)=ab(10) + ab(1) + ab(4) + ab(3) - 3*ab(15) + ab(2)
      ab(1)=xa*ab(1)
      ab(1)= - ab(9) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + ab(6)

      tmp = CA**2*z1*z2*z3*ab(1)
      res(2,1,0) = real(tmp,ki)

      return
      end

