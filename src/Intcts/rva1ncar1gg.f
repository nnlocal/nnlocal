      subroutine rva1ncar1gg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(2,2,-4:0)
      real(ki) xlf
      complex(ki) ab(30)
      include 'rva1ncar1_functions.h'
      xlf=dble(nf)
      xlf=0d0
      
!##### QUARTIC POLE #####

      
      ab(1)=1/( - 8 + 4*z3**(-1))
      ab(2)=1/( - 4 + 2*z3**(-1))
      ab(3)=1/( - 2 + z3**(-1))
      ab(4)=xa**2
      ab(5)=1 - xa
      ab(5)=xb*ab(5)
      ab(4)=ab(5) - 3 - ab(4)
      ab(4)=xb*ab(4)
      ab(5)= - 3 + xa
      ab(5)=xa*ab(5)
      ab(4)=ab(5) + ab(4)
      ab(4)=ab(1)*ab(4)
      ab(5)=xb*xa*ab(3)
      ab(4)=ab(4) + ab(2) + ab(5)

      tmp = CA**2*z1*z6*ab(4)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=1/( - 2 + z3**(-1))
      ab(2)=xb - 3 + xa
      ab(2)=xb*ab(2)
      ab(2)=ab(2) + 2 - xa

      tmp = CA**2*z1*z6*ab(2)*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=1/( - 2 + z3**(-1))
      ab(2)= - xb + 3 - xa
      ab(2)=xb*ab(2)
      ab(2)=ab(2) - 2 + xa
      ab(2)=l2*ab(2)
      ab(2)=ab(2) + 1 - xb

      tmp = 2*CA**2*z1*z6*ab(2)*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=1/( - 8 + 4*z3**(-1))
      ab(2)=1/( - 6 + 3*z3**(-1))
      ab(3)=1/( - 4 + 2*z3**(-1))
      ab(4)=1/( - 2 + z3**(-1))
      ab(5)=2*l2
      ab(6)=2 - 3*l2
      ab(6)=ab(6)*ab(5)
      ab(7)=l2**2
      ab(7)=2*ab(7)
      ab(8)=5*zeta3
      ab(9)=ab(7) + ab(8)
      ab(9)=xa*ab(9)
      ab(10)=Pi**2
      ab(7)=xb*ab(7)
      ab(6)=ab(7) + ab(9) + ab(6) + 4*l3 + ab(10)
      ab(6)=xb*ab(6)
      ab(5)= - 3 + ab(5)
      ab(5)=l2*ab(5)
      ab(7)=1 - l2
      ab(7)=l2*ab(7)
      ab(7)= - l4 + ab(7)
      ab(7)=xa*ab(7)
      ab(5)=ab(7) + ab(5) - 2*l3 + l4
      ab(5)=2*ab(5) + ab(6)
      ab(5)=ab(4)*ab(5)
      ab(6)=5*xa
      ab(7)=zeta3*ab(1)
      ab(9)= - ab(7)*ab(6)
      ab(10)=ab(10)*ab(2)
      ab(9)= - ab(10) + ab(9)
      ab(9)=xa*ab(9)
      ab(11)= - ab(6) + 5
      ab(11)=ab(7)*ab(11)
      ab(11)= - ab(10) + ab(11)
      ab(11)=xb*ab(11)
      ab(9)=ab(11) - 15*ab(7) + ab(9)
      ab(9)=xb*ab(9)
      ab(6)=ab(6) - 15
      ab(6)=ab(7)*ab(6)
      ab(6)=ab(10) + ab(6)
      ab(6)=xa*ab(6)
      ab(7)=ab(3)*ab(8)
      ab(5)=ab(5) + ab(9) + ab(6) - 2*ab(10) + ab(7)

      tmp = CA**2*z1*z6*ab(5)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=1/( - 80 + 40*z3**(-1))
      ab(2)=1/( - 40 + 20*z3**(-1))
      ab(3)=1/( - 20 + 10*z3**(-1))
      ab(4)=1/( - 8 + 4*z3**(-1))
      ab(5)=1/( - 6 + 3*z3**(-1))
      ab(6)=1/( - 4 + 2*z3**(-1))
      ab(7)=1/( - 2 + z3**(-1))
      ab(8)=CA**2*z1*z6
      ab(9)=ab(8)*xb
      ab(10)=ab(9) - ab(8)
      ab(11)=ab(10)*xa
      ab(12)= - 4*ab(8) + ab(9)
      ab(12)=xb*ab(12)
      ab(13)=3*ab(8)
      ab(12)=ab(12) + ab(11) + ab(13)
      ab(12)=xa*ab(12)
      ab(14)=ab(13) - ab(9)
      ab(14)=ab(14)*xb
      ab(15)=2*ab(8)
      ab(16)=ab(14) - ab(15)
      ab(12)=ab(12) + ab(16)
      ab(12)=plg4half*ab(12)
      ab(17)=ab(9)*xa
      ab(18)= - 1.0_ki/2.0_ki*ab(8) - ab(17)
      ab(19)=l1**4
      ab(18)=ab(18)*ab(19)
      ab(20)=Pi**2
      ab(21)=ab(20)*l1**2
      ab(22)=ab(17)*ab(21)
      ab(12)=6*ab(12) + ab(18) + ab(22)
      ab(12)=ab(7)*ab(12)
      ab(18)=ab(8)*xb**2
      ab(18)=ab(11) + ab(18)
      ab(22)=ab(13) + ab(18)
      ab(22)=ab(22)*xa
      ab(14)=ab(22) + ab(14)
      ab(22)=ab(14)*ab(4)
      ab(23)= - ab(7)*ab(17)
      ab(23)=ab(22) + ab(23)
      ab(24)=21*l1
      ab(23)=ab(24)*ab(23)
      ab(11)= - ab(11) + ab(16)
      ab(16)=2*ab(7)
      ab(11)=ab(11)*ab(16)
      ab(11)=ab(11) + ab(23)
      ab(11)=zeta3*ab(11)
      ab(14)=ab(1)*ab(14)
      ab(17)=ab(3)*ab(17)
      ab(14)=ab(14) - ab(17)
      ab(17)=ab(2)*ab(13)
      ab(14)=ab(17) - 3*ab(14)
      ab(14)=ab(14)*Pi**4
      ab(17)=ab(19) - ab(21)
      ab(17)=ab(17)*ab(22)
      ab(19)= - xa*ab(13)
      ab(19)=ab(19) + 7*ab(8) - 4*ab(9)
      ab(19)=ab(7)*ab(19)
      ab(22)=ab(9)*ab(7)
      ab(23)=l2*ab(22)
      ab(19)=ab(19) + 4*ab(23)
      ab(19)=l2*ab(19)
      ab(22)=ab(20)*ab(22)
      ab(19)= - 2*ab(22) + ab(19)
      ab(19)=l2*ab(19)
      ab(22)=l2**2
      ab(22)=ab(20) - 2*ab(22)
      ab(15)=ab(15) + ab(18)
      ab(15)=l2*ab(15)*ab(22)
      ab(18)=ab(10)*ab(20)
      ab(15)=ab(18) + ab(15)
      ab(15)=ab(5)*ab(15)
      ab(18)=ab(8)*xa
      ab(9)= - ab(18) + ab(13) - 2*ab(9)
      ab(13)=l2*ab(7)
      ab(9)=ab(9)*ab(13)
      ab(10)= - l3*ab(7)*ab(10)
      ab(9)=ab(9) + ab(10)
      ab(9)=l3*ab(9)
      ab(10)=l3*ab(16)
      ab(10)=ab(13) + ab(10)
      ab(13)=l4*ab(7)
      ab(10)=2*ab(10) + ab(13)
      ab(13)=ab(18) - ab(8)
      ab(10)=l4*ab(13)*ab(10)
      ab(13)= - zeta3*ab(24)
      ab(13)=ab(21) + ab(13)
      ab(8)=ab(6)*ab(8)*ab(13)

      tmp = ab(8) + 4*ab(9) + ab(10) + ab(11) + ab(12) + ab(14) + 2*
     & ab(15) + ab(17) + ab(19)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=xb - 1
      ab(2)=2*ab(1)
      ab(3)= - xa*ab(1)
      ab(3)=ab(2) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)= - 3*ab(1) + ab(3)
      ab(3)=xa*ab(3)
      ab(2)=ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) - ab(1)

      tmp = CA**2*z1*z4*z6*ab(1)
      res(2,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - xb + 1
      ab(2)= - l6 - l5 + l2 + l1
      ab(1)=ab(2)*ab(1)
      ab(1)= - 1 + ab(1)
      ab(2)=2*ab(1)
      ab(3)= - xa*ab(1)
      ab(3)=ab(2) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)= - 3*ab(1) + ab(3)
      ab(3)=xa*ab(3)
      ab(2)=ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) - ab(1)

      tmp = 2*CA**2*z1*z4*z6*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=z1*z4
      ab(2)=ab(1)*CA**2
      ab(3)=xa*z6
      ab(4)=ab(2)*ab(3)
      ab(5)=2*z6
      ab(5)=ab(2)*ab(5)
      ab(6)=ab(4) - ab(5)
      ab(7)=ab(6)*xa
      ab(8)=3*z6
      ab(8)=ab(8)*ab(2)
      ab(7)=ab(7) + ab(8)
      ab(7)=ab(7)*xa
      ab(5)=ab(7) - ab(5)
      ab(5)=ab(5)*xa
      ab(2)=ab(2)*z6
      ab(5)=ab(5) + ab(2)
      ab(7)=xb - 1
      ab(8)=ab(7)*ab(5)
      ab(9)=ab(8)*l6
      ab(10)=2*ab(5) - ab(9)
      ab(10)=l6*ab(10)
      ab(9)=ab(9) - ab(5)
      ab(11)=ab(8)*l2
      ab(12)= - ab(11) + 2*ab(9)
      ab(13)=l2*ab(12)
      ab(14)=2*l1
      ab(15)=ab(8)*ab(14)
      ab(12)=ab(15) - ab(12)
      ab(12)=l5*ab(12)
      ab(10)=ab(12) + ab(10) + ab(13)
      ab(9)= - ab(11) + ab(9)
      ab(8)= - l1*ab(8)
      ab(8)=2*ab(9) + ab(8)
      ab(8)=ab(8)*ab(14)
      ab(4)=ab(4) - ab(2)
      ab(3)=z6 - ab(3)
      ab(1)=TR*xlf*CA*ab(1)*ab(3)
      ab(1)=1.0_ki/3.0_ki*ab(1) + 1.0_ki/6.0_ki*ab(4)
      ab(1)=ab(1)*ab(7)*xa**2
      ab(3)=1.0_ki/3.0_ki*xa
      ab(3)=ab(3)*ab(6)
      ab(3)=ab(3) + ab(2)
      ab(3)=ab(3)*xa
      ab(3)=ab(3) - 2.0_ki/3.0_ki*ab(2)
      ab(3)=ab(3)*xa
      ab(2)=ab(3) + 1.0_ki/3.0_ki*ab(2)
      ab(2)=Pi**2*ab(2)*ab(7)
      ab(3)=l3*ab(5)

      tmp = ab(1) + ab(2) - 4*ab(3) + ab(8) + 2*ab(10)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z4
      ab(2)=ab(1)*CA**2
      ab(3)=xa*z6
      ab(4)=ab(2)*ab(3)
      ab(5)=2*z6
      ab(6)=ab(2)*ab(5)
      ab(7)=ab(4) - ab(6)
      ab(8)=ab(7)*xa
      ab(9)=3*z6
      ab(9)=ab(9)*ab(2)
      ab(8)=ab(8) + ab(9)
      ab(8)=ab(8)*xa
      ab(8)=ab(8) - ab(6)
      ab(8)=ab(8)*xa
      ab(2)=ab(2)*z6
      ab(8)=ab(8) + ab(2)
      ab(9)=ab(8)*l2
      ab(10)=ab(8)*xb
      ab(11)=l2*ab(10)
      ab(11)=ab(11) - ab(9)
      ab(12)=ab(11) + ab(8)
      ab(13)=ab(10) - ab(8)
      ab(14)=ab(13)*l6
      ab(15)= - 2*ab(12) + ab(14)
      ab(16)=4*l6
      ab(15)=ab(15)*ab(16)
      ab(17)=1.0_ki/3.0_ki*xa
      ab(7)=ab(7)*ab(17)
      ab(7)=ab(7) + ab(2)
      ab(7)=ab(7)*xa
      ab(7)=ab(7) - 2.0_ki/3.0_ki*ab(2)
      ab(7)=ab(7)*xa
      ab(7)=ab(7) + 1.0_ki/3.0_ki*ab(2)
      ab(18)=xb - 1
      ab(19)=ab(7)*ab(18)
      ab(20)=l1*ab(19)
      ab(20)=ab(20) - ab(14) + ab(12)
      ab(21)=4*l1
      ab(20)=ab(20)*ab(21)
      ab(22)=2*ab(8)
      ab(23)=ab(22) - ab(9)
      ab(24)=4*l2
      ab(23)=ab(23)*ab(24)
      ab(25)=l2**2
      ab(26)=4*ab(25)
      ab(27)=ab(8)*ab(26)
      ab(28)=ab(4) - ab(2)
      ab(29)=xa**2
      ab(30)=1.0_ki/3.0_ki*ab(29)
      ab(30)=ab(30)*ab(28)
      ab(27)=ab(27) - ab(30)
      ab(27)=ab(27)*xb
      ab(23)=ab(30) + ab(23) + ab(27)
      ab(15)=ab(20) + ab(15) + ab(23)
      ab(15)=l1*ab(15)
      ab(20)=ab(11) + ab(22)
      ab(14)= - ab(14) + ab(20)
      ab(14)=ab(14)*ab(16)
      ab(22)=2*l6
      ab(27)= - l1 + ab(22)
      ab(27)=ab(13)*ab(27)
      ab(20)=ab(27) - ab(20)
      ab(20)=ab(20)*ab(21)
      ab(21)= - 4*ab(8) - ab(9)
      ab(21)=l2*ab(21)
      ab(10)=ab(25)*ab(10)
      ab(10)=ab(20) + ab(14) + ab(21) + ab(10)
      ab(10)=l5*ab(10)
      ab(1)=ab(1)*CA
      ab(3)=ab(3)*ab(1)
      ab(14)=z6*ab(1)
      ab(20)=ab(3) - ab(14)
      ab(20)=ab(20)*ab(29)
      ab(21)= - ab(22) + 2*l1
      ab(18)=ab(21)*ab(18)*ab(20)
      ab(1)=ab(5)*ab(1)
      ab(5)= - 13*ab(14) + 10*ab(3)
      ab(5)=ab(5)*ab(17)
      ab(5)=ab(1) + ab(5)
      ab(5)=xa*ab(5)
      ab(5)= - ab(14) + ab(5)
      ab(3)=10*ab(14) - 7*ab(3)
      ab(3)=ab(3)*ab(17)
      ab(1)= - ab(1) + ab(3)
      ab(1)=xa*ab(1)
      ab(1)=ab(14) + ab(1)
      ab(3)=l2*ab(20)
      ab(1)=2*ab(1) + ab(3)
      ab(1)=xb*ab(1)
      ab(1)=ab(1) + 2*ab(5) - ab(3) + ab(18)
      ab(1)=TR*xlf*ab(1)
      ab(3)=23*ab(2) - 17*ab(4)
      ab(5)=1.0_ki/6.0_ki*xa
      ab(3)=ab(3)*ab(5)
      ab(3)= - ab(6) + ab(3)
      ab(3)=xa*ab(3)
      ab(1)=ab(1) + ab(2) + ab(3)
      ab(3)= - l6*ab(19)
      ab(3)=ab(3) + ab(12)
      ab(3)=ab(3)*ab(16)
      ab(3)=ab(3) - ab(23)
      ab(3)=l6*ab(3)
      ab(12)= - ab(19)*ab(21)
      ab(12)= - 2*ab(7) - ab(11) + ab(12)
      ab(12)=ab(12)*Pi**2
      ab(14)=d2 + d3
      ab(14)=2*ab(14)
      ab(11)=ab(11)*ab(14)
      ab(14)=t2 + t3
      ab(16)=t1 + zeta3
      ab(14)=t4 + 2*ab(16) - 4*ab(14)
      ab(13)=ab(13)*ab(14)
      ab(14)= - l2*ab(7)
      ab(14)=ab(14) + ab(8)
      ab(14)=ab(14)*ab(24)
      ab(16)=ab(28)*ab(29)
      ab(16)=1.0_ki/6.0_ki*ab(16)
      ab(14)=ab(16) + ab(14)
      ab(14)=l2*ab(14)
      ab(4)= - 17*ab(2) + 11*ab(4)
      ab(4)=ab(4)*ab(5)
      ab(4)=ab(6) + ab(4)
      ab(4)=xa*ab(4)
      ab(2)= - ab(2) + ab(4)
      ab(4)=ab(7)*ab(26)
      ab(4)= - ab(16) + ab(4)
      ab(4)=l2*ab(4)
      ab(2)=1.0_ki/3.0_ki*ab(2) + ab(4)
      ab(2)=xb*ab(2)
      ab(4)= - l5 + l1 - l6
      ab(4)=ab(8)*ab(4)
      ab(4)=ab(9) + ab(4)
      ab(5)=l3*ab(8)
      ab(4)=2*ab(4) + ab(5)
      ab(4)=l3*ab(4)

      tmp = 1.0_ki/3.0_ki*ab(1) + ab(2) + ab(3) + 4*ab(4) + ab(10) + ab(11)
     &  + ab(12) + ab(13) + ab(14) + ab(15)
      res(2,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,2,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp =  0
      res(2,2,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=xa**4
      ab(2)= - xb*xa**5
      ab(1)=ab(1) + ab(2)
      ab(1)=xb*ab(1)
      ab(2)=xa**3
      ab(1)= - ab(2) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=xa**2
      ab(1)= - ab(2) + ab(1)
      ab(1)=xb*ab(1)
      ab(1)=xa + ab(1)
      ab(1)=xb*ab(1)
      ab(1)= - 1 + ab(1)

      tmp = 4*CA**2*z1*z2*z3*z4*z5*ab(1)
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l8 + l3 - l5
      ab(2)=l7 + l9
      ab(1)= - l2 - l6 - l10 + 3*ab(2) - 2*ab(1)
      ab(2)=xa*xb**5
      ab(3)=xb**4
      ab(2)=ab(3) - ab(2)
      ab(2)=xa*ab(2)
      ab(3)=xb**3
      ab(2)=ab(2) - ab(3)
      ab(3)=z1*z2*z3*z4*z5*CA**2
      ab(2)=xa*ab(3)*ab(2)
      ab(4)=ab(3)*xb**2
      ab(2)=ab(2) - ab(4)
      ab(2)=ab(2)*xa
      ab(4)=ab(3)*xb
      ab(2)=ab(2) + ab(4)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - ab(3)

      tmp = 4*ab(2)*ab(1)
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l6 + l2
      ab(2)= - l10 + 2*l5
      ab(1)=6*l3 - ab(2) + 3*ab(1)
      ab(1)= - 5*l7 + 2*ab(1)
      ab(1)=ab(1)*l7
      ab(3)=l2 - ab(2)
      ab(4)=ab(3) + 2*l3
      ab(5)= - l6 + 3*l7 - ab(4)
      ab(6)=ab(5) - 2*l8
      ab(7)=3*l9
      ab(6)=ab(7) + 2*ab(6)
      ab(6)=ab(6)*ab(7)
      ab(5)=ab(5) - l8
      ab(7)=4*l8
      ab(5)=ab(5)*ab(7)
      ab(3)=ab(3) + l3
      ab(7)=4*l3
      ab(3)=ab(3)*ab(7)
      ab(4)=l6 + 2*ab(4)
      ab(4)=ab(4)*l6
      ab(2)= - l2 + 2*ab(2)
      ab(2)=ab(2)*l2
      ab(7)=l10**2
      ab(1)= - ab(7) - ab(4) + ab(1) - ab(6) - ab(3) + ab(2) + ab(5)
      ab(2)=Pi**2
      ab(3)= - 1 + 2*ab(2)
      ab(3)= - 1.0_ki/3.0_ki*ab(3) - ab(1)
      ab(4)=xa**4
      ab(5)=ab(3)*ab(4)
      ab(1)=2.0_ki/3.0_ki*ab(2) + ab(1)
      ab(2)=xb*ab(1)*xa**5
      ab(2)=ab(5) + ab(2)
      ab(2)=xb*ab(2)
      ab(5)=ab(1)*xa**3
      ab(2)=ab(5) + ab(2)
      ab(2)=xb*ab(2)
      ab(5)=xa**2
      ab(3)= - ab(3)*ab(5)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(3)= - xa*ab(1)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)
      ab(1)=CA*ab(1)
      ab(2)=xb**2
      ab(3)= - ab(4)*ab(2)
      ab(3)=ab(5) + ab(3)
      ab(2)=ab(2)*xlf*TR*ab(3)
      ab(1)=2.0_ki/3.0_ki*ab(2) + ab(1)

      tmp = 2*CA*z1*z2*z3*z4*z5*ab(1)
      res(2,2,0) = real(tmp,ki)

      return
      end

