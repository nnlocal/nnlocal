      subroutine a2carsggg(xa,xb,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) xa,xb,res(2,2,-4:0)
      complex(ki) ab(51)
      include 'a2cars_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)=1.0_ki/2.0_ki*xa
      ab(2)= - 7 + xa
      ab(2)=ab(2)*ab(1)
      ab(2)=7 + ab(2)
      ab(2)=xa*ab(2)
      ab(3)= - 7.0_ki/2.0_ki + xa
      ab(3)=xa*ab(3)
      ab(4)= - 1 + xa
      ab(4)=xb*ab(4)
      ab(3)=1.0_ki/2.0_ki*ab(4) + 5.0_ki/2.0_ki + ab(3)
      ab(3)=xb*ab(3)
      ab(2)=ab(3) - 4 + ab(2)
      ab(2)=xb*ab(2)
      ab(3)=5 - xa
      ab(1)=ab(3)*ab(1)
      ab(1)= - 4 + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(2) + 2 + ab(1)
      ab(2)=z3*CA

      tmp = 3*z1*z2*ab(2)**2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=1.0_ki/2.0_ki*xa
      ab(2)=11*xa
      ab(3)= - 149 + ab(2)
      ab(3)=ab(3)*ab(1)
      ab(3)=293 + ab(3)
      ab(1)=ab(3)*ab(1)
      ab(3)=1.0_ki/2.0_ki*xb
      ab(4)= - 83 + ab(2)
      ab(4)=ab(4)*ab(3)
      ab(5)= - 221.0_ki/2.0_ki + ab(2)
      ab(5)=xa*ab(5)
      ab(4)=ab(4) + 415.0_ki/2.0_ki + ab(5)
      ab(3)=ab(4)*ab(3)
      ab(1)=ab(3) - 166 + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=127 - ab(2)
      ab(2)=xa*ab(2)
      ab(2)= - 94 + 1.0_ki/4.0_ki*ab(2)
      ab(2)=xa*ab(2)
      ab(1)=ab(1) + 83 + ab(2)
      ab(2)=z3*CA

      tmp = 1.0_ki/3.0_ki*z1*z2*ab(2)**2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=Pi**2
      ab(2)=3*ab(1)
      ab(3)=32*l2
      ab(4)=ab(2) - 163.0_ki/9.0_ki + ab(3)
      ab(5)=ab(2) - 61.0_ki/9.0_ki
      ab(6)=1.0_ki/4.0_ki*xa
      ab(6)=ab(5)*ab(6)
      ab(6)=ab(6) - 16*l2
      ab(7)=15.0_ki/4.0_ki*ab(1)
      ab(8)= - ab(7) + 437.0_ki/36.0_ki + ab(6)
      ab(8)=xa*ab(8)
      ab(4)=2*ab(4) + ab(8)
      ab(4)=xa*ab(4)
      ab(8)=21.0_ki/4.0_ki*ab(1)
      ab(5)=xa*ab(5)
      ab(3)= - 1.0_ki/2.0_ki*ab(5) + ab(8) - 691.0_ki/36.0_ki + ab(3)
      ab(3)=xa*ab(3)
      ab(5)=3.0_ki/4.0_ki*ab(1) - 193.0_ki/36.0_ki - ab(6)
      ab(5)=xb*ab(5)
      ab(3)=ab(5) + ab(3) - ab(7) + 1253.0_ki/36.0_ki - 80*l2
      ab(3)=xb*ab(3)
      ab(5)=ab(8) - 559.0_ki/36.0_ki - ab(6)
      ab(5)=xa*ab(5)
      ab(1)=ab(5) - 21.0_ki/2.0_ki*ab(1) + 967.0_ki/18.0_ki - 96*l2
      ab(1)=xa*ab(1)
      ab(2)=ab(2) + 64*l2
      ab(5)= - 301.0_ki/9.0_ki + ab(2)
      ab(1)=ab(3) + 2*ab(5) + ab(1)
      ab(1)=xb*ab(1)
      ab(1)=ab(1) + ab(4) + 337.0_ki/9.0_ki - ab(2)
      ab(2)=z3*CA

      tmp = z1*z2*ab(2)**2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=14*l1
      ab(2)=Pi**2
      ab(3)= - 2113.0_ki/3.0_ki + 293.0_ki/2.0_ki*ab(2)
      ab(3)=1.0_ki/6.0_ki*ab(3) + ab(1)
      ab(4)=9*l2
      ab(5)=ab(4) - 11.0_ki/3.0_ki
      ab(6)=4*l2
      ab(7)= - ab(5)*ab(6)
      ab(8)=147*zeta3
      ab(9)= - 17.0_ki/3.0_ki + 1.0_ki/2.0_ki*ab(2)
      ab(10)= - 11.0_ki/3.0_ki*ab(9) - 4*l1
      ab(10)=1.0_ki/3.0_ki*ab(10) - 42*zeta3
      ab(10)=xa*ab(10)
      ab(3)=ab(10) + ab(7) + 1.0_ki/3.0_ki*ab(3) + ab(8)
      ab(3)=xa*ab(3)
      ab(7)=2*l1
      ab(9)=ab(7) + 11.0_ki/6.0_ki*ab(9)
      ab(10)=21*zeta3
      ab(9)=ab(10) + 1.0_ki/3.0_ki*ab(9)
      ab(9)=ab(9)*xa
      ab(11)=2*l2
      ab(5)=ab(5)*ab(11)
      ab(5)=ab(9) + ab(5)
      ab(12)= - 589.0_ki/3.0_ki + 119.0_ki/2.0_ki*ab(2)
      ab(7)=1.0_ki/6.0_ki*ab(12) + ab(7)
      ab(7)=1.0_ki/3.0_ki*ab(7) + ab(10) - ab(5)
      ab(7)=xb*ab(7)
      ab(10)= - 97.0_ki/3.0_ki + 45*l2
      ab(10)=ab(10)*ab(11)
      ab(12)=105*zeta3 + 10.0_ki/3.0_ki*l1
      ab(13)=4*l3
      ab(3)=ab(7) + ab(3) + ab(10) - 595.0_ki/36.0_ki*ab(2) + 3341.0_ki/54.0_ki
     &  - ab(13) - ab(12)
      ab(3)=xb*ab(3)
      ab(7)= - 1711.0_ki/3.0_ki + 185.0_ki/2.0_ki*ab(2)
      ab(1)=1.0_ki/6.0_ki*ab(7) + ab(1)
      ab(1)=1.0_ki/3.0_ki*ab(1) + ab(8) - ab(5)
      ab(1)=xa*ab(1)
      ab(5)=12*l4
      ab(7)= - 7 + ab(4)
      ab(7)=l2*ab(7)
      ab(1)=ab(1) + 12*ab(7) - 294*zeta3 - 28.0_ki/3.0_ki*l1 - 401.0_ki/18.0_ki
     & *ab(2) + ab(5) + 2695.0_ki/27.0_ki - ab(13)
      ab(1)=xa*ab(1)
      ab(7)=84*zeta3 + 8.0_ki/3.0_ki*l1
      ab(8)=ab(7) + 119.0_ki/9.0_ki*ab(2)
      ab(10)=ab(4) - 29.0_ki/3.0_ki
      ab(13)=8*l2
      ab(14)= - ab(10)*ab(13)
      ab(15)=l4 - l3
      ab(15)= - 733.0_ki/27.0_ki - 3*ab(15)
      ab(14)=ab(14) + 2*ab(15) + ab(8)
      ab(1)=ab(3) + 2*ab(14) + ab(1)
      ab(1)=xb*ab(1)
      ab(3)=ab(10)*ab(11)
      ab(3)=ab(9) + ab(3) - 163.0_ki/36.0_ki*ab(2) + 1337.0_ki/54.0_ki + ab(5)
     &  - ab(12)
      ab(3)=xa*ab(3)
      ab(5)=35.0_ki/3.0_ki - ab(4)
      ab(5)=ab(5)*ab(6)
      ab(6)= - 9*l4 - 433.0_ki/27.0_ki + l3
      ab(2)=ab(5) + 2*ab(6) + 65.0_ki/9.0_ki*ab(2) + ab(7)
      ab(2)=2*ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(3)= - 41.0_ki/3.0_ki + ab(4)
      ab(3)=ab(3)*ab(13)
      ab(4)=6*l4 + 389.0_ki/27.0_ki - 2*l3
      ab(1)=ab(1) + ab(2) + ab(3) + 4*ab(4) - ab(8)
      ab(2)=z3*CA

      tmp = z1*z2*ab(2)**2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l1**2
      ab(2)=Pi**2
      ab(1)=ab(1) - ab(2)
      ab(3)=19.0_ki/3.0_ki*l1
      ab(3)=ab(1)*ab(3)
      ab(3)=ab(3) + 20.0_ki/9.0_ki + 133*zeta3
      ab(4)=2*l1
      ab(4)=ab(3)*ab(4)
      ab(5)=7*l3
      ab(6)=10.0_ki/3.0_ki + 9*l4
      ab(6)=2*ab(6) + ab(5)
      ab(6)=l3*ab(6)
      ab(7)=1.0_ki/3.0_ki*l2
      ab(8)=19*l2
      ab(9)= - 53 + ab(8)
      ab(9)=ab(9)*ab(7)
      ab(10)=3*l4
      ab(11)= - 5*l3 + 41.0_ki/9.0_ki + ab(10)
      ab(9)=2*ab(11) + ab(9)
      ab(11)=4*l2
      ab(9)=ab(9)*ab(11)
      ab(12)=1.0_ki/3.0_ki*ab(2)
      ab(13)=7.0_ki/15.0_ki*ab(2)
      ab(14)= - ab(13) + 211.0_ki/9.0_ki - 40*l2
      ab(14)=ab(14)*ab(12)
      ab(15)=l4**2
      ab(16)=152*plg4half
      ab(17)=ab(16) - 2621.0_ki/81.0_ki - 12*d1
      ab(6)=ab(4) + ab(14) + ab(9) + 200.0_ki/3.0_ki*zeta3 + ab(6) + 2*
     & ab(17) + 15*ab(15)
      ab(9)=19.0_ki/6.0_ki*l1
      ab(1)=ab(1)*ab(9)
      ab(1)=ab(1) + 10.0_ki/9.0_ki + 133.0_ki/2.0_ki*zeta3
      ab(9)=ab(1)*l1
      ab(14)= - 73.0_ki/3.0_ki + 7.0_ki/5.0_ki*ab(2)
      ab(17)=1.0_ki/36.0_ki*ab(2)
      ab(17)=ab(14)*ab(17)
      ab(18)=76*plg4half
      ab(19)=11.0_ki/3.0_ki*zeta3 + ab(18) - 565.0_ki/81.0_ki
      ab(17)= - ab(9) + ab(17) - ab(19)
      ab(17)=ab(17)*xa
      ab(20)=5*ab(9) + 380*plg4half
      ab(21)=50 - ab(8)
      ab(21)=ab(21)*ab(7)
      ab(22)=8*l4
      ab(21)=ab(21) + 6*l3 - 67.0_ki/9.0_ki - ab(22)
      ab(23)=2*l2
      ab(21)=ab(21)*ab(23)
      ab(24)=20*l2
      ab(25)=7.0_ki/12.0_ki*ab(2)
      ab(26)=ab(25) - 497.0_ki/36.0_ki + ab(24)
      ab(26)=ab(26)*ab(12)
      ab(27)= - 10*ab(15) + 16*d1
      ab(28)=l3*l4
      ab(21)= - ab(17) + ab(26) + ab(21) - 133.0_ki/3.0_ki*zeta3 - 12*
     & ab(28) + 4037.0_ki/81.0_ki + ab(27) - ab(20)
      ab(21)=xa*ab(21)
      ab(6)=2*ab(6) + ab(21)
      ab(6)=xa*ab(6)
      ab(14)=ab(14)*ab(2)
      ab(21)= - l1*ab(3)
      ab(14)=ab(21) - 2*ab(19) + 1.0_ki/18.0_ki*ab(14)
      ab(14)=xa*ab(14)
      ab(19)=7*l1
      ab(1)=ab(1)*ab(19)
      ab(21)=532*plg4half
      ab(1)=ab(1) + ab(21)
      ab(26)=49.0_ki/60.0_ki*ab(2)
      ab(28)=8*l2
      ab(29)=155.0_ki/36.0_ki - ab(28)
      ab(29)=5*ab(29) - ab(26)
      ab(29)=ab(29)*ab(12)
      ab(30)=ab(8) - 11
      ab(30)=ab(30)*l2
      ab(30)=ab(30) + 67.0_ki/3.0_ki
      ab(31)=l2*ab(30)
      ab(14)=ab(14) + ab(29) + 4.0_ki/3.0_ki*ab(31) + 233.0_ki/3.0_ki*zeta3 - 
     & 6379.0_ki/81.0_ki + ab(1)
      ab(14)=xa*ab(14)
      ab(29)=2.0_ki/3.0_ki*l2
      ab(30)=ab(30)*ab(29)
      ab(17)=ab(30) + ab(17)
      ab(11)=41.0_ki/36.0_ki - ab(11)
      ab(11)=5*ab(11) - 7.0_ki/60.0_ki*ab(2)
      ab(11)=ab(11)*ab(12)
      ab(9)=ab(9) + ab(11) + 89.0_ki/3.0_ki*zeta3 - 1777.0_ki/81.0_ki + ab(18)
     &  + ab(17)
      ab(9)=xb*ab(9)
      ab(11)=2*l3
      ab(18)= - 22.0_ki/3.0_ki - ab(5)
      ab(18)=ab(18)*ab(11)
      ab(30)=106 - 95*l2
      ab(30)=ab(30)*ab(7)
      ab(30)=ab(30) - 401.0_ki/9.0_ki + 22*l3
      ab(30)=ab(30)*ab(23)
      ab(25)=ab(25) - 1601.0_ki/36.0_ki + 100*l2
      ab(25)=ab(25)*ab(12)
      ab(9)=ab(9) + ab(14) + ab(25) + ab(30) - 445.0_ki/3.0_ki*zeta3 + 
     & ab(18) + 10091.0_ki/81.0_ki - ab(20)
      ab(9)=xb*ab(9)
      ab(14)= - 35 + ab(8)
      ab(7)=ab(14)*ab(7)
      ab(7)=ab(7) - 9*l3 + 91.0_ki/9.0_ki + l4
      ab(7)=ab(7)*ab(28)
      ab(13)=ab(13) + 80*l2
      ab(14)=421.0_ki/9.0_ki - ab(13)
      ab(14)=ab(14)*ab(12)
      ab(15)=5*ab(15)
      ab(18)=38*plg4half - 1109.0_ki/81.0_ki - d1
      ab(20)=32.0_ki/3.0_ki + ab(10)
      ab(20)=2*ab(20) + 21*l3
      ab(20)=l3*ab(20)
      ab(7)=ab(4) + ab(14) + ab(7) + 356.0_ki/3.0_ki*zeta3 + ab(20) + 8*
     & ab(18) + ab(15)
      ab(10)= - 10.0_ki/3.0_ki - ab(10)
      ab(10)=2*ab(10) - ab(5)
      ab(10)=l3*ab(10)
      ab(14)=26 - ab(8)
      ab(14)=l2*ab(14)
      ab(14)=ab(14) + 14*l3 - 77.0_ki/3.0_ki - 4*l4
      ab(14)=ab(14)*ab(23)
      ab(10)=ab(14) - 311.0_ki/3.0_ki*zeta3 + ab(10) - ab(15) - ab(21) + 
     & 8149.0_ki/81.0_ki + 8*d1
      ab(14)= - ab(26) + 643.0_ki/36.0_ki - ab(24)
      ab(14)=ab(14)*ab(12)
      ab(1)=ab(14) + 155.0_ki/3.0_ki*zeta3 - 5167.0_ki/81.0_ki + ab(1) + ab(17)
      ab(1)=xa*ab(1)
      ab(14)= - 239.0_ki/54.0_ki + ab(28)
      ab(14)=5*ab(14) + 49.0_ki/90.0_ki*ab(2)
      ab(2)=ab(14)*ab(2)
      ab(3)= - ab(3)*ab(19)
      ab(1)=ab(1) + ab(3) + 2*ab(10) + ab(2)
      ab(1)=xa*ab(1)
      ab(1)=ab(9) + 2*ab(7) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)= - ab(5) - 7 - 6*l4
      ab(2)=ab(2)*ab(11)
      ab(3)=56 - ab(8)
      ab(3)=ab(3)*ab(29)
      ab(3)=ab(3) + 28*l3 - 197.0_ki/9.0_ki - ab(22)
      ab(3)=ab(3)*ab(23)
      ab(2)=ab(3) - 178.0_ki/3.0_ki*zeta3 + ab(2) - ab(16) + 4715.0_ki/81.0_ki
     &  + ab(27)
      ab(3)= - 493.0_ki/9.0_ki + ab(13)
      ab(3)=ab(3)*ab(12)
      ab(1)=ab(1) + ab(6) - ab(4) + 2*ab(2) + ab(3)
      ab(2)=z3*CA

      tmp = z1*z2*ab(2)**2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(2,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=xa - 1
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + 1
      ab(2)=CA**2*z1*z2*z4*z5
      ab(1)=xa*ab(2)*ab(1)
      ab(1)=ab(1) + ab(2)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) - ab(2)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) + ab(2)
      ab(2)=xb - 1

      tmp = 6*ab(2)*ab(1)
      res(2,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)= - 8*l6 + 16*l2 + 8*l1
      ab(2)=4*l5
      ab(3)=ab(1) - ab(2)
      ab(4)=ab(3) - 55.0_ki/3.0_ki
      ab(4)=ab(4)*xb
      ab(3)=ab(4) - ab(3) + 79.0_ki/3.0_ki
      ab(4)= - xa*ab(3)
      ab(5)=12*l5 - ab(1)
      ab(6)= - 47.0_ki/3.0_ki - ab(5)
      ab(6)=xb*ab(6)
      ab(4)=ab(4) + ab(6) + 71.0_ki/3.0_ki + ab(5)
      ab(4)=xa*ab(4)
      ab(1)=ab(1) + ab(2)
      ab(2)=1 - ab(1)
      ab(2)=xb*ab(2)
      ab(2)=ab(4) + ab(2) - 9 + ab(1)
      ab(2)=xa*ab(2)
      ab(4)=1 + ab(5)
      ab(4)=xb*ab(4)
      ab(2)=ab(2) + ab(4) - 9 - ab(5)
      ab(2)=xa*ab(2)
      ab(4)= - 47.0_ki/3.0_ki + ab(1)
      ab(4)=xb*ab(4)
      ab(1)=ab(2) + ab(4) + 71.0_ki/3.0_ki - ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) - ab(3)

      tmp = CA**2*z1*z2*z4*z5*ab(1)
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - 12*d5 + 16*d6
      ab(2)=Pi**2
      ab(3)= - ab(1) + 20.0_ki/3.0_ki*ab(2)
      ab(4)=2*l5
      ab(5)=ab(4) - 37.0_ki/3.0_ki
      ab(5)=ab(5)*l5
      ab(5)=ab(5) - ab(3)
      ab(6)=9*l2 + 10*l1
      ab(7)=6*l5
      ab(8)=ab(6) - ab(7)
      ab(9)=ab(8) - 47.0_ki/3.0_ki
      ab(10)=2*l2
      ab(9)=ab(9)*ab(10)
      ab(11)=2*l1
      ab(12)= - l1 + 14*l5
      ab(13)=23.0_ki/3.0_ki + ab(12)
      ab(13)=ab(13)*ab(11)
      ab(14)=2*l6
      ab(15)= - l6 + ab(11) + 10*l2
      ab(16)=ab(15) - 8*l5
      ab(17)= - 23.0_ki/3.0_ki + ab(16)
      ab(17)=ab(17)*ab(14)
      ab(13)=ab(17) - ab(9) + ab(13) + 7.0_ki/18.0_ki - ab(5)
      ab(13)=xb*ab(13)
      ab(17)=ab(8) - 55.0_ki/3.0_ki
      ab(17)=ab(17)*l2
      ab(7)=ab(7) - l1
      ab(18)=ab(7) - 1.0_ki/3.0_ki
      ab(18)=ab(18)*l1
      ab(17)=ab(18) - ab(17)
      ab(18)=4*l5
      ab(15)=ab(15) - ab(18)
      ab(19)= - 1.0_ki/3.0_ki - ab(15)
      ab(19)=l6*ab(19)
      ab(20)=2*d5
      ab(21)= - 8 - 5*ab(2)
      ab(19)=ab(19) + 26*l5 + 1.0_ki/3.0_ki*ab(21) - ab(20) - ab(17)
      ab(19)=xb*ab(19)
      ab(21)=ab(8) - 97.0_ki/3.0_ki
      ab(21)=ab(21)*l2
      ab(7)=ab(7) + 17.0_ki/3.0_ki
      ab(7)=ab(7)*l1
      ab(7)= - 47.0_ki/3.0_ki + 2*l3 + ab(7) - ab(21)
      ab(21)= - 17.0_ki/3.0_ki + ab(15)
      ab(21)=l6*ab(21)
      ab(19)=ab(19) + ab(21) - 28*l5 + ab(20) + 5.0_ki/3.0_ki*ab(2) + ab(7)
      ab(19)=xa*ab(19)
      ab(21)=ab(4) + 23.0_ki/3.0_ki
      ab(21)=ab(21)*l5
      ab(22)=4*l3
      ab(3)=ab(21) - ab(3) - ab(22)
      ab(21)=ab(8) - 89.0_ki/3.0_ki
      ab(21)=ab(21)*ab(10)
      ab(23)= - 41.0_ki/3.0_ki - ab(12)
      ab(23)=ab(23)*ab(11)
      ab(24)=41.0_ki/3.0_ki - ab(16)
      ab(24)=ab(24)*ab(14)
      ab(13)=2*ab(19) + ab(13) + ab(24) + ab(21) + ab(23) + 557.0_ki/18.0_ki
     &  + ab(3)
      ab(13)=xa*ab(13)
      ab(8)=ab(8) - 1
      ab(8)=ab(8)*ab(10)
      ab(1)= - ab(1) + 2.0_ki/3.0_ki*ab(2)
      ab(19)=ab(4) + l1
      ab(23)= - 35.0_ki/3.0_ki + ab(19)
      ab(23)=ab(23)*ab(11)
      ab(24)=35.0_ki/3.0_ki - ab(15)
      ab(24)=ab(24)*ab(14)
      ab(25)= - 35.0_ki/3.0_ki - ab(18)
      ab(25)=l5*ab(25)
      ab(23)=ab(24) + ab(8) + ab(23) + ab(25) + 217.0_ki/18.0_ki + ab(1)
      ab(23)=xb*ab(23)
      ab(24)=ab(4) + 5
      ab(6)= - ab(6) + 3*ab(24)
      ab(6)=ab(6)*ab(10)
      ab(10)=ab(1) - ab(22)
      ab(22)=53.0_ki/3.0_ki - ab(19)
      ab(22)=ab(22)*ab(11)
      ab(24)= - 53.0_ki/3.0_ki + ab(15)
      ab(24)=ab(24)*ab(14)
      ab(25)=71.0_ki/3.0_ki + ab(18)
      ab(25)=l5*ab(25)
      ab(13)=ab(13) + ab(23) + ab(24) + ab(6) + ab(22) + ab(25) - 253.0_ki
     & /18.0_ki - ab(10)
      ab(13)=xa*ab(13)
      ab(4)= - 5.0_ki/3.0_ki - ab(4)
      ab(4)=7*ab(4) + l1
      ab(4)=ab(4)*ab(11)
      ab(22)=35.0_ki/3.0_ki - ab(16)
      ab(22)=ab(22)*ab(14)
      ab(4)=ab(22) + ab(8) + ab(4) + 223.0_ki/18.0_ki + ab(5)
      ab(4)=xb*ab(4)
      ab(5)=53.0_ki/3.0_ki + ab(12)
      ab(5)=ab(5)*ab(11)
      ab(8)= - 53.0_ki/3.0_ki + ab(16)
      ab(8)=ab(8)*ab(14)
      ab(3)=ab(13) + ab(4) + ab(8) + ab(6) + ab(5) - 259.0_ki/18.0_ki - 
     & ab(3)
      ab(3)=xa*ab(3)
      ab(4)=23.0_ki/3.0_ki - ab(19)
      ab(4)=ab(4)*ab(11)
      ab(5)= - 23.0_ki/3.0_ki + ab(15)
      ab(5)=ab(5)*ab(14)
      ab(6)= - 121.0_ki/3.0_ki + ab(18)
      ab(6)=l5*ab(6)
      ab(1)=ab(5) - ab(9) + ab(4) + ab(6) + 13.0_ki/18.0_ki - ab(1)
      ab(1)=xb*ab(1)
      ab(4)= - 41.0_ki/3.0_ki + ab(19)
      ab(4)=ab(4)*ab(11)
      ab(5)=41.0_ki/3.0_ki - ab(15)
      ab(5)=ab(5)*ab(14)
      ab(6)=85.0_ki/3.0_ki - ab(18)
      ab(6)=l5*ab(6)
      ab(1)=ab(3) + ab(1) + ab(5) + ab(21) + ab(4) + ab(6) 
     & + 551.0_ki/18.0_ki + ab(10)
      ab(1)=xa*ab(1)
      ab(3)= - 2 - ab(2)
      ab(3)=2.0_ki/3.0_ki*ab(3) + d5
      ab(4)=l5**2
      ab(5)= - 1.0_ki/3.0_ki - ab(16)
      ab(5)=l6*ab(5)
      ab(3)=ab(5) + 2*ab(3) - ab(4) - ab(17)
      ab(3)=xb*ab(3)
      ab(4)= - 2 + l5
      ab(4)=l5*ab(4)
      ab(5)= - 17.0_ki/3.0_ki + ab(16)
      ab(5)=l6*ab(5)
      ab(2)=ab(3) + ab(5) + ab(4) - ab(20) + 4.0_ki/3.0_ki*ab(2) + ab(7)
      ab(1)=2*ab(2) + ab(1)

      tmp = CA**2*z1*z2*z4*z5*ab(1)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=z1*z2*z4*z5
      ab(2)=CA**2
      ab(3)=ab(1)*ab(2)
      ab(4)=xa*ab(3)
      ab(5)=ab(4) - ab(3)
      ab(5)=ab(5)*xa
      ab(6)=5*ab(3)
      ab(7)=ab(5) + ab(6)
      ab(7)=ab(7)*xa
      ab(7)=ab(7) + ab(3)
      ab(7)=ab(7)*xa
      ab(7)=ab(7) - ab(6)
      ab(7)=ab(7)*xa
      ab(7)=ab(7) + ab(6)
      ab(8)=xb - 1
      ab(7)=ab(7)*ab(8)
      ab(9)=2*l6
      ab(10)=2*l1
      ab(11)=ab(9) - ab(10)
      ab(12)=ab(7)*ab(11)
      ab(13)= - 8*ab(3) + 91.0_ki/3.0_ki*ab(4)
      ab(13)=xa*ab(13)
      ab(13)= - 100.0_ki/3.0_ki*ab(3) + ab(13)
      ab(13)=xa*ab(13)
      ab(13)= - 22.0_ki/3.0_ki*ab(3) + ab(13)
      ab(13)=xa*ab(13)
      ab(13)=10.0_ki/3.0_ki*ab(3) + ab(13)
      ab(13)=xa*ab(13)
      ab(14)= - 4*ab(3) - 79.0_ki/3.0_ki*ab(4)
      ab(14)=xa*ab(14)
      ab(14)=64.0_ki/3.0_ki*ab(3) + ab(14)
      ab(14)=xa*ab(14)
      ab(14)=58.0_ki/3.0_ki*ab(3) + ab(14)
      ab(14)=xa*ab(14)
      ab(15)=26.0_ki/3.0_ki*ab(3)
      ab(14)=ab(15) + ab(14)
      ab(14)=xa*ab(14)
      ab(16)=35.0_ki/3.0_ki*ab(3)
      ab(14)= - ab(16) + ab(14)
      ab(14)=xb*ab(14)
      ab(17)=9*ab(3)
      ab(18)=ab(17) - ab(4)
      ab(18)=ab(18)*xa
      ab(2)=3*ab(2)
      ab(1)=ab(2)*ab(1)
      ab(2)=ab(18) + ab(1)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - ab(17)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - ab(1)
      ab(2)=ab(2)*xa
      ab(2)=ab(2) - ab(6)
      ab(18)=2*l5
      ab(19)=ab(8)*ab(18)
      ab(2)=ab(2)*ab(19)
      ab(20)=ab(1) - ab(4)
      ab(20)=ab(20)*xa
      ab(20)=ab(20) + ab(1)
      ab(20)=ab(20)*xa
      ab(20)=ab(20) - ab(1)
      ab(20)=ab(20)*xa
      ab(20)=ab(20) - ab(1)
      ab(20)=ab(20)*xa
      ab(20)=ab(20) + ab(3)
      ab(20)=ab(8)*ab(20)
      ab(21)=l2*ab(20)
      ab(2)= - 4*ab(21) + ab(2) + ab(14) + 23.0_ki/3.0_ki*ab(3) + ab(13) + 
     & ab(12)
      ab(2)=d5*ab(2)
      ab(12)=3*ab(4)
      ab(13)=ab(12) + ab(3)
      ab(13)=ab(13)*xa
      ab(14)=7*ab(3)
      ab(13)=ab(13) + ab(14)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) - ab(3)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) - ab(14)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) + ab(1)
      ab(13)= - ab(11)*ab(8)*ab(13)
      ab(21)=101.0_ki/3.0_ki*ab(4)
      ab(22)= - 12*ab(3) - ab(21)
      ab(22)=xa*ab(22)
      ab(22)=86.0_ki/3.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(22)= - 8.0_ki/3.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(15)=ab(15) + ab(22)
      ab(15)=xa*ab(15)
      ab(21)=28*ab(3) + ab(21)
      ab(21)=xa*ab(21)
      ab(21)= - 38.0_ki/3.0_ki*ab(3) + ab(21)
      ab(21)=xa*ab(21)
      ab(21)= - 40.0_ki/3.0_ki*ab(3) + ab(21)
      ab(21)=xa*ab(21)
      ab(21)= - 74.0_ki/3.0_ki*ab(3) + ab(21)
      ab(21)=xa*ab(21)
      ab(22)=55.0_ki/3.0_ki*ab(3)
      ab(21)= - ab(22) + ab(21)
      ab(21)=xb*ab(21)
      ab(23)=ab(6) - ab(4)
      ab(23)=ab(23)*xa
      ab(24)=ab(23) + ab(1)
      ab(24)=ab(24)*xa
      ab(24)=ab(24) - ab(6)
      ab(24)=ab(24)*xa
      ab(24)=ab(24) - ab(1)
      ab(24)=ab(24)*xa
      ab(24)=ab(24) - ab(3)
      ab(24)= - ab(24)*ab(19)
      ab(25)=ab(4) + ab(3)
      ab(25)=ab(25)*xa
      ab(26)=ab(25) - ab(3)
      ab(26)=ab(26)*xa
      ab(26)=ab(26) - ab(3)
      ab(27)=xb*xa
      ab(28)= - xa + ab(27)
      ab(26)=l2*ab(26)*ab(28)
      ab(13)=16*ab(26) + ab(24) + ab(21) + ab(22) + ab(15) + ab(13)
      ab(13)=d6*ab(13)
      ab(15)= - 55*ab(4) + 47*ab(3)
      ab(21)=1.0_ki/3.0_ki*xa
      ab(15)=ab(15)*ab(21)
      ab(15)=ab(15) - ab(3)
      ab(15)=ab(15)*xa
      ab(15)=ab(15) - ab(3)
      ab(15)=ab(15)*xa
      ab(15)=ab(15) + 47.0_ki/3.0_ki*ab(3)
      ab(15)=ab(15)*xa
      ab(15)=ab(15) - ab(22)
      ab(22)=ab(5) + ab(3)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) + ab(3)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) - ab(3)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) + ab(3)
      ab(24)=l6 - l1
      ab(26)= - 11*l2 + 3*ab(24)
      ab(26)=ab(22)*ab(26)
      ab(28)= - ab(14) - ab(4)
      ab(28)=xa*ab(28)
      ab(28)= - ab(17) + ab(28)
      ab(28)=xa*ab(28)
      ab(28)=ab(14) + ab(28)
      ab(28)=xa*ab(28)
      ab(28)=ab(17) + ab(28)
      ab(28)=xa*ab(28)
      ab(28)= - ab(3) + ab(28)
      ab(28)=l5*ab(28)
      ab(26)=ab(28) + ab(26) - ab(15)
      ab(28)=l3*ab(22)
      ab(26)=2*ab(26) + 7*ab(28)
      ab(26)=l3*ab(26)
      ab(28)=1.0_ki/3.0_ki*ab(3)
      ab(29)= - ab(28) + 19*ab(4)
      ab(29)=ab(29)*xa
      ab(30)=67.0_ki/3.0_ki*ab(3)
      ab(29)=ab(29) - ab(30)
      ab(29)=ab(29)*xa
      ab(29)=ab(29) - ab(30)
      ab(29)=ab(29)*xa
      ab(29)=ab(29) - ab(28)
      ab(29)=ab(29)*xa
      ab(29)=ab(29) + 19*ab(3)
      ab(29)=ab(29)*ab(8)
      ab(22)=ab(22)*ab(8)
      ab(11)=4*l5 - ab(11)
      ab(11)=ab(22)*ab(11)
      ab(11)=ab(11) + ab(29)
      ab(11)=d3*ab(11)
      ab(2)=ab(26) + ab(11) + ab(2) + ab(13)
      ab(11)=xa**2
      ab(13)=ab(11)*ab(3)
      ab(26)=2*ab(3)
      ab(29)=ab(13) + ab(26)
      ab(29)=ab(29)*ab(11)
      ab(29)=ab(29) - ab(26)
      ab(29)=ab(29)*xa
      ab(29)=ab(29) + ab(3)
      ab(31)=4*l1
      ab(29)= - ab(31)*ab(29)
      ab(32)= - ab(6) + 13.0_ki/3.0_ki*ab(4)
      ab(32)=ab(32)*xa
      ab(33)=11.0_ki/3.0_ki*ab(3)
      ab(32)=ab(32) + ab(33)
      ab(32)=ab(32)*xa
      ab(32)=ab(32) + ab(6)
      ab(32)=ab(32)*xa
      ab(32)=ab(32) - ab(33)
      ab(32)=ab(32)*xa
      ab(34)=13.0_ki/3.0_ki*ab(3)
      ab(32)=ab(32) + ab(34)
      ab(35)=2*l2
      ab(32)= - ab(35)*ab(32)
      ab(29)=ab(32) + ab(29)
      ab(29)=ab(8)*ab(29)
      ab(15)=ab(15)*xb
      ab(32)=49*ab(3) - 53*ab(4)
      ab(32)=ab(32)*ab(21)
      ab(32)= - ab(17) + ab(32)
      ab(32)=xa*ab(32)
      ab(32)= - ab(17) + ab(32)
      ab(32)=xa*ab(32)
      ab(36)=49.0_ki/3.0_ki*ab(3)
      ab(32)=ab(36) + ab(32)
      ab(32)=xa*ab(32)
      ab(32)= - 53.0_ki/3.0_ki*ab(3) + ab(32)
      ab(37)=l6*ab(22)
      ab(38)=7*ab(4)
      ab(39)=ab(38) - ab(6)
      ab(39)=ab(39)*xa
      ab(40)=ab(39) + ab(17)
      ab(40)=ab(40)*xa
      ab(40)=ab(40) + ab(6)
      ab(40)=ab(40)*xa
      ab(40)=ab(40) - ab(17)
      ab(40)=ab(40)*xa
      ab(40)=ab(40) + ab(14)
      ab(19)=ab(40)*ab(19)
      ab(19)=ab(19) + 11*ab(37) + 2*ab(32) - ab(15) + ab(29)
      ab(19)=ab(19)*ab(35)
      ab(29)= - ab(18) + 6*l1
      ab(22)=ab(22)*ab(29)
      ab(29)= - 77*ab(3) + 85*ab(4)
      ab(29)=ab(29)*ab(21)
      ab(32)=11*ab(3)
      ab(29)=ab(32) + ab(29)
      ab(29)=xa*ab(29)
      ab(29)=ab(32) + ab(29)
      ab(29)=xa*ab(29)
      ab(35)=77.0_ki/3.0_ki*ab(3)
      ab(29)= - ab(35) + ab(29)
      ab(29)=xa*ab(29)
      ab(37)=23*ab(3)
      ab(40)=ab(37) + 15*ab(5)
      ab(40)=ab(40)*xa
      ab(41)=15*ab(3)
      ab(40)=ab(40) + ab(41)
      ab(40)=ab(40)*xa
      ab(40)=ab(40) - ab(37)
      ab(40)=ab(40)*xa
      ab(40)=ab(40) + ab(37)
      ab(42)=ab(8)*l6
      ab(40)= - ab(40)*ab(42)
      ab(43)=85.0_ki/3.0_ki*ab(3)
      ab(15)=ab(40) + ab(15) + ab(43) + ab(29) + ab(22)
      ab(15)=ab(15)*ab(18)
      ab(22)= - 557.0_ki/3.0_ki*ab(3) + 188*ab(4)
      ab(22)=xa*ab(22)
      ab(22)=253.0_ki/3.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(22)=259.0_ki/3.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(22)= - 551.0_ki/3.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(29)=7.0_ki/3.0_ki*ab(3)
      ab(40)= - ab(29) + 32*ab(4)
      ab(40)=xa*ab(40)
      ab(40)= - 217.0_ki/3.0_ki*ab(3) + ab(40)
      ab(40)=xa*ab(40)
      ab(40)= - 223.0_ki/3.0_ki*ab(3) + ab(40)
      ab(40)=xa*ab(40)
      ab(40)= - ab(34) + ab(40)
      ab(40)=xa*ab(40)
      ab(40)=32*ab(3) + ab(40)
      ab(40)=xb*ab(40)
      ab(22)=ab(40) + 188*ab(3) + ab(22)
      ab(40)=13*ab(4)
      ab(44)=25*ab(3)
      ab(45)=ab(40) - ab(44)
      ab(45)=ab(45)*xa
      ab(46)=31*ab(3)
      ab(45)=ab(45) + ab(46)
      ab(45)=ab(45)*xa
      ab(45)=ab(45) + ab(46)
      ab(45)=ab(45)*xa
      ab(44)=ab(45) - ab(44)
      ab(44)=ab(44)*xa
      ab(45)=13*ab(3)
      ab(44)=ab(44) + ab(45)
      ab(47)=ab(37) + ab(4)
      ab(47)=ab(47)*xa
      ab(48)=35*ab(3)
      ab(47)=ab(47) - ab(48)
      ab(47)=ab(47)*xa
      ab(47)=ab(47) - ab(48)
      ab(47)=ab(47)*xa
      ab(37)=ab(47) + ab(37)
      ab(37)=ab(37)*xa
      ab(37)=ab(37) + ab(3)
      ab(37)=ab(37)*xb
      ab(44)=ab(37) + 2*ab(44)
      ab(47)=9*ab(4)
      ab(48)=ab(47) - ab(45)
      ab(48)=ab(48)*xa
      ab(48)=ab(48) + ab(3)
      ab(48)=ab(48)*xa
      ab(48)=ab(48) + ab(45)
      ab(48)=ab(48)*xa
      ab(48)=ab(48) - ab(3)
      ab(48)=ab(48)*xa
      ab(48)=ab(48) + ab(6)
      ab(49)=ab(8)*l1
      ab(48)= - ab(48)*ab(49)
      ab(48)= - 2.0_ki/3.0_ki*ab(44) + ab(48)
      ab(48)=ab(48)*ab(10)
      ab(50)=ab(5) - ab(3)
      ab(50)=ab(50)*xa
      ab(50)=ab(50) + ab(3)
      ab(50)=ab(50)*xa
      ab(50)=ab(50) + ab(3)
      ab(50)=ab(50)*xa
      ab(50)=ab(50) - ab(3)
      ab(51)= - l6 + ab(10)
      ab(50)=ab(51)*ab(8)*ab(50)
      ab(44)=1.0_ki/3.0_ki*ab(44) + ab(50)
      ab(50)=4*l6
      ab(44)=ab(44)*ab(50)
      ab(15)=ab(19) + ab(15) + ab(44) + 1.0_ki/3.0_ki*ab(22) + ab(48)
      ab(15)=l2*ab(15)
      ab(19)=25.0_ki/6.0_ki*ab(3)
      ab(22)= - ab(19) - 24*ab(4)
      ab(22)=xa*ab(22)
      ab(22)=95.0_ki/6.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(19)=ab(19) + ab(22)
      ab(19)=xa*ab(19)
      ab(19)=61.0_ki/6.0_ki*ab(3) + ab(19)
      ab(19)=xa*ab(19)
      ab(22)=37.0_ki/6.0_ki*ab(3)
      ab(44)=ab(22) + 26*ab(4)
      ab(44)=xa*ab(44)
      ab(48)=35.0_ki/6.0_ki*ab(3)
      ab(44)= - ab(48) + ab(44)
      ab(44)=xa*ab(44)
      ab(22)= - ab(22) + ab(44)
      ab(22)=xa*ab(22)
      ab(22)= - 121.0_ki/6.0_ki*ab(3) + ab(22)
      ab(22)=ab(22)*ab(27)
      ab(44)=ab(20)*ab(10)
      ab(13)=ab(13) - ab(3)
      ab(11)=ab(13)*ab(11)
      ab(11)=ab(11) - ab(3)
      ab(11)= - ab(50)*ab(11)*ab(8)
      ab(13)=ab(26) + ab(4)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) - ab(3)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) - ab(26)
      ab(13)=ab(13)*xa
      ab(13)=ab(13) + ab(3)
      ab(50)=2.0_ki/3.0_ki*l5
      ab(50)=ab(8)*ab(50)
      ab(13)= - ab(13)*ab(50)
      ab(11)=ab(13) + ab(11) + ab(44) + ab(22) + 6*ab(3) + ab(19)
      ab(11)=l5*ab(11)
      ab(13)=41.0_ki/3.0_ki*ab(3) - ab(40)
      ab(13)=xa*ab(13)
      ab(13)= - ab(30) + ab(13)
      ab(13)=xa*ab(13)
      ab(19)=29*ab(3)
      ab(13)= - ab(19) + ab(13)
      ab(13)=xa*ab(13)
      ab(13)=ab(35) + ab(13)
      ab(13)=xa*ab(13)
      ab(22)=ab(28) + ab(38)
      ab(22)=xa*ab(22)
      ab(22)=ab(36) + ab(22)
      ab(22)=xa*ab(22)
      ab(22)=ab(41) + ab(22)
      ab(22)=xa*ab(22)
      ab(22)= - 59.0_ki/3.0_ki*ab(3) + ab(22)
      ab(22)=xa*ab(22)
      ab(16)= - ab(16) + ab(22)
      ab(16)=xb*ab(16)
      ab(22)=ab(12) - ab(14)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) + ab(1)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) + ab(14)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) - ab(1)
      ab(22)=ab(22)*xa
      ab(22)=ab(22) + ab(14)
      ab(35)=ab(8)*ab(10)
      ab(22)= - ab(22)*ab(35)
      ab(36)=ab(39) + ab(32)
      ab(36)=ab(36)*xa
      ab(36)=ab(36) + ab(6)
      ab(36)=ab(36)*xa
      ab(36)=ab(36) - ab(32)
      ab(36)=ab(36)*xa
      ab(36)=ab(36) + ab(17)
      ab(36)=ab(36)*ab(42)
      ab(13)=ab(36) + ab(22) + ab(16) - ab(29) + ab(13)
      ab(13)=ab(13)*ab(9)
      ab(16)=2*ab(4)
      ab(22)=ab(16) - ab(3)
      ab(22)=ab(22)*xa
      ab(36)= - ab(32) + 65.0_ki/3.0_ki*ab(22)
      ab(36)=xa*ab(36)
      ab(36)=65.0_ki/3.0_ki*ab(3) + ab(36)
      ab(36)=xa*ab(36)
      ab(30)= - ab(30) + ab(36)
      ab(30)=xa*ab(30)
      ab(36)= - ab(45) - 100*ab(4)
      ab(36)=ab(36)*ab(21)
      ab(36)=ab(6) + ab(36)
      ab(36)=xa*ab(36)
      ab(34)=ab(34) + ab(36)
      ab(34)=xa*ab(34)
      ab(34)=ab(43) + ab(34)
      ab(34)=ab(34)*ab(27)
      ab(36)=17*ab(3)
      ab(39)=ab(36) - ab(47)
      ab(39)=ab(39)*xa
      ab(39)=ab(39) - ab(3)
      ab(39)=ab(39)*xa
      ab(36)=ab(39) - ab(36)
      ab(36)=ab(36)*xa
      ab(36)=ab(36) + ab(3)
      ab(36)=ab(36)*xa
      ab(36)=ab(36) - ab(17)
      ab(36)= - ab(36)*ab(49)
      ab(30)=ab(36) + ab(34) + 10*ab(3) + ab(30)
      ab(30)=ab(30)*ab(10)
      ab(34)= - 463.0_ki/6.0_ki*ab(3) - 308*ab(4)
      ab(34)=xa*ab(34)
      ab(34)=737.0_ki/6.0_ki*ab(3) + ab(34)
      ab(34)=xa*ab(34)
      ab(34)=439.0_ki/6.0_ki*ab(3) + ab(34)
      ab(34)=xa*ab(34)
      ab(34)=1087.0_ki/6.0_ki*ab(3) + ab(34)
      ab(34)=xa*ab(34)
      ab(36)=4*ab(4)
      ab(39)=19.0_ki/6.0_ki*ab(3) - ab(36)
      ab(39)=xa*ab(39)
      ab(39)= - 317.0_ki/6.0_ki*ab(3) + ab(39)
      ab(39)=xa*ab(39)
      ab(39)=5.0_ki/6.0_ki*ab(3) + ab(39)
      ab(39)=xa*ab(39)
      ab(39)=365.0_ki/6.0_ki*ab(3) + ab(39)
      ab(27)=ab(39)*ab(27)
      ab(27)=ab(34) + ab(27)
      ab(11)=ab(11) + ab(13) + 1.0_ki/3.0_ki*ab(27) + ab(30)
      ab(11)=l5*ab(11)
      ab(13)= - 71*ab(3) + 499.0_ki/2.0_ki*ab(4)
      ab(13)=xa*ab(13)
      ab(13)= - 85*ab(3) + ab(13)
      ab(13)=ab(13)*ab(21)
      ab(13)=16*ab(3) + ab(13)
      ab(13)=xa*ab(13)
      ab(13)= - 116.0_ki/3.0_ki*ab(3) + ab(13)
      ab(13)=xa*ab(13)
      ab(27)= - 43*ab(3) - 391.0_ki/2.0_ki*ab(4)
      ab(27)=xa*ab(27)
      ab(27)=67*ab(3) + ab(27)
      ab(21)=ab(27)*ab(21)
      ab(21)=22*ab(3) + ab(21)
      ab(21)=xa*ab(21)
      ab(21)=134.0_ki/3.0_ki*ab(3) + ab(21)
      ab(21)=xa*ab(21)
      ab(21)= - ab(48) + ab(21)
      ab(21)=xb*ab(21)
      ab(13)=ab(21) + 119.0_ki/6.0_ki*ab(3) + ab(13)
      ab(21)= - ab(45) + 19.0_ki/3.0_ki*ab(4)
      ab(21)=ab(21)*xa
      ab(21)=ab(21) + ab(29)
      ab(21)=ab(21)*xa
      ab(21)=ab(21) + ab(45)
      ab(21)=ab(21)*xa
      ab(21)=ab(21) - ab(29)
      ab(21)=ab(21)*xa
      ab(21)=ab(21) + ab(17)
      ab(21)=ab(21)*ab(49)
      ab(27)=ab(36) - ab(14)
      ab(27)=ab(27)*xa
      ab(27)=ab(27) + ab(26)
      ab(27)=ab(27)*xa
      ab(27)=ab(27) + ab(14)
      ab(27)=ab(27)*xa
      ab(27)=ab(27) - ab(26)
      ab(27)=ab(27)*xa
      ab(27)=ab(27) + ab(6)
      ab(27)= - ab(27)*ab(42)
      ab(30)=ab(36) - ab(6)
      ab(30)=ab(30)*xa
      ab(30)=ab(30) - ab(26)
      ab(30)=ab(30)*xa
      ab(30)=ab(30) + ab(6)
      ab(30)=ab(30)*xa
      ab(30)=ab(30) + ab(26)
      ab(30)=ab(30)*xa
      ab(30)=ab(30) - ab(3)
      ab(30)= - ab(30)*ab(50)
      ab(16)=ab(16) - ab(6)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) - ab(26)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) + ab(6)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) + ab(26)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) + ab(3)
      ab(16)=l2*ab(16)*ab(8)
      ab(13)=4.0_ki/3.0_ki*ab(16) + ab(30) + 2.0_ki/3.0_ki*ab(27) + 1.0_ki/3.0_ki*
     & ab(13) + ab(21)
      ab(13)=ab(13)*Pi**2
      ab(5)=ab(17) + 5*ab(5)
      ab(5)=ab(5)*xa
      ab(5)=ab(5) + ab(6)
      ab(5)=ab(5)*xa
      ab(5)=ab(5) - ab(17)
      ab(5)=ab(5)*xa
      ab(5)=ab(5) + ab(17)
      ab(5)=t9*ab(5)
      ab(16)=ab(38) - ab(32)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) + ab(1)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) + ab(32)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) - ab(1)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) + ab(14)
      ab(16)=t6*ab(16)
      ab(17)=ab(1) + ab(4)
      ab(17)=ab(17)*xa
      ab(17)=ab(17) + ab(3)
      ab(17)=ab(17)*xa
      ab(17)=ab(17) - ab(1)
      ab(17)=ab(17)*xa
      ab(17)=ab(17) - ab(3)
      ab(17)=ab(17)*xa
      ab(17)=ab(17) - ab(1)
      ab(17)=t2*ab(17)
      ab(21)=ab(23) + ab(14)
      ab(21)=ab(21)*xa
      ab(21)=ab(21) - ab(6)
      ab(21)=ab(21)*xa
      ab(14)=ab(21) - ab(14)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + ab(1)
      ab(14)=t1*ab(14)
      ab(5)=ab(5) + ab(16) + ab(17) + ab(14)
      ab(14)= - ab(46) + 25*ab(4)
      ab(16)=3*xa
      ab(14)=ab(14)*ab(16)
      ab(16)=107*ab(3)
      ab(14)=ab(14) + ab(16)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + 93*ab(3)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) - ab(16)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + 125*ab(3)
      ab(14)= - zeta3*ab(14)
      ab(5)=1.0_ki/2.0_ki*ab(14) + 4*ab(5)
      ab(5)=ab(8)*ab(5)
      ab(14)= - ab(19) + 5*ab(4)
      ab(14)=ab(14)*xa
      ab(17)=41*ab(3)
      ab(14)=ab(14) + ab(17)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + ab(17)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) - ab(19)
      ab(14)=ab(14)*xa
      ab(6)=ab(14) + ab(37) + ab(6)
      ab(6)=1.0_ki/3.0_ki*ab(6)
      ab(14)=ab(22) + ab(26)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + ab(3)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) - ab(26)
      ab(14)=ab(14)*xa
      ab(14)=ab(14) + ab(3)
      ab(14)=ab(8)*ab(14)
      ab(17)=ab(14)*ab(31)
      ab(19)= - ab(3) + 5.0_ki/3.0_ki*ab(4)
      ab(19)=ab(19)*xa
      ab(19)=ab(19) + ab(1)
      ab(19)=ab(19)*xa
      ab(19)=ab(19) + ab(3)
      ab(19)=ab(19)*xa
      ab(1)=ab(19) - ab(1)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) + ab(29)
      ab(1)= - ab(9)*ab(1)*ab(8)
      ab(1)=ab(1) - ab(6) + ab(17)
      ab(1)=l6*ab(1)
      ab(9)=352*ab(4) + 701*ab(3)
      ab(9)=ab(9)*xa
      ab(9)=ab(9) - ab(16)
      ab(9)=ab(9)*xa
      ab(9)=ab(9) - 263*ab(3)
      ab(9)=ab(9)*xa
      ab(9)=ab(9) - 367*ab(3)
      ab(9)=ab(9)*xa
      ab(9)=ab(9) - 560*ab(3)
      ab(9)=ab(9)*xb
      ab(16)=364*ab(4) + 977*ab(3)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) - 527*ab(3)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) - 683*ab(3)
      ab(16)=ab(16)*xa
      ab(16)=ab(16) - 91*ab(3)
      ab(16)=ab(16)*xa
      ab(9)=ab(9) - ab(16) + 548*ab(3)
      ab(9)=1.0_ki/9.0_ki*ab(9)
      ab(16)=l1*ab(20)
      ab(16)=ab(6) + ab(16)
      ab(10)=ab(16)*ab(10)
      ab(1)=ab(1) - ab(9) + ab(10)
      ab(1)=l6*ab(1)
      ab(10)=ab(33) - ab(12)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) + ab(3)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) - ab(33)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) - ab(3)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) + ab(28)
      ab(10)= - ab(10)*ab(35)
      ab(6)= - ab(6) + ab(10)
      ab(6)=l1*ab(6)
      ab(6)=ab(9) + ab(6)
      ab(6)=l1*ab(6)
      ab(9)= - l1 - l6
      ab(9)=ab(18)*ab(9)
      ab(10)= - l1**2
      ab(12)=l6**2
      ab(9)=ab(9) + ab(10) + ab(12)
      ab(10)=ab(25) + ab(3)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) - ab(3)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) - ab(3)
      ab(10)=ab(10)*xa
      ab(10)=ab(10) - ab(3)
      ab(8)=ab(10)*ab(8)
      ab(9)=l7*ab(8)*ab(9)
      ab(7)=t11*ab(7)
      ab(7)=ab(9) + ab(7)
      ab(9)=ab(18) - ab(24)
      ab(10)=d2 + d4
      ab(10)= - 4*ab(10)
      ab(9)=ab(10)*ab(9)
      ab(10)= - t8 + t3 - t7
      ab(12)=t4 - t10
      ab(9)=ab(9) - 12*ab(12) + 4*ab(10)
      ab(8)=ab(8)*ab(9)
      ab(9)=628*ab(3) - 505*ab(4)
      ab(9)=xa*ab(9)
      ab(9)= - 1523*ab(3) + 2*ab(9)
      ab(9)=xa*ab(9)
      ab(9)= - 512*ab(3) + ab(9)
      ab(9)=xa*ab(9)
      ab(9)=1727*ab(3) + ab(9)
      ab(9)=xa*ab(9)
      ab(4)= - 1277*ab(3) + 1298*ab(4)
      ab(4)=xa*ab(4)
      ab(4)=872*ab(3) + ab(4)
      ab(4)=xa*ab(4)
      ab(4)= - 157*ab(3) + ab(4)
      ab(4)=xa*ab(4)
      ab(4)= - 1766*ab(3) + ab(4)
      ab(4)=xa*ab(4)
      ab(4)=1838*ab(3) + ab(4)
      ab(4)=xb*ab(4)
      ab(3)=ab(4) - 1550*ab(3) + ab(9)
      ab(4)=t5*ab(14)

      tmp = ab(1) + 2*ab(2) + 1.0_ki/27.0_ki*ab(3) + 8*ab(4) + ab(5) + 
     & ab(6) + 4*ab(7) + ab(8) + ab(11) + ab(13) + ab(15)
      res(2,1,0) = real(tmp,ki)

      return
      end

