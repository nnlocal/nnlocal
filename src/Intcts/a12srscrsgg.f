      subroutine a12srscrsgg(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(18)
      include 'a12srscrs_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)= - 1 + lam

      tmp = 2*CA*z1*Ta2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=8*l1
      ab(2)=11.0_ki/3.0_ki - ab(1)
      ab(2)=lam*ab(2)
      ab(1)=ab(2) - 35.0_ki/3.0_ki + ab(1)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=35 - 11*lam
      ab(2)= - 1 + lam
      ab(2)=l1*ab(2)
      ab(1)=1.0_ki/3.0_ki*ab(1) + 4*ab(2)
      ab(1)=l1*ab(1)
      ab(1)=2*ab(1) + 35.0_ki/9.0_ki*lam - 101.0_ki/9.0_ki + 16*l2

      tmp = 2*CA*z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=8*l1
      ab(2)=1 - lam
      ab(2)=ab(2)*ab(1)
      ab(2)=ab(2) - 35 + 11*lam
      ab(2)=l1*ab(2)
      ab(3)=Pi**2
      ab(4)= - 35.0_ki/9.0_ki + ab(3)
      ab(4)=lam*ab(4)
      ab(2)=1.0_ki/3.0_ki*ab(2) + ab(4) - ab(3) + 101.0_ki/9.0_ki - 16*l2
      ab(1)=ab(2)*ab(1)
      ab(2)=11.0_ki/3.0_ki*ab(3)
      ab(3)=29.0_ki/27.0_ki + zeta3
      ab(3)=16*ab(3) - ab(2)
      ab(3)=lam*ab(3)
      ab(4)=11.0_ki/3.0_ki - 4*l2
      ab(4)=l2*ab(4)
      ab(4)=ab(4) - zeta3
      ab(4)= - 163.0_ki/27.0_ki + 2*ab(4)
      ab(1)=ab(1) + ab(3) + 8*ab(4) + ab(2)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=Pi**2
      ab(2)=35.0_ki/3.0_ki - ab(1)
      ab(2)=lam*ab(2)
      ab(2)=ab(1) + ab(2)
      ab(3)=35 - 11*lam
      ab(4)= - 1 + lam
      ab(4)=l1*ab(4)
      ab(3)=1.0_ki/3.0_ki*ab(3) + ab(4)
      ab(3)=l1*ab(3)
      ab(2)=2.0_ki/3.0_ki*ab(3) - 101.0_ki/9.0_ki + 16*l2 + 1.0_ki/3.0_ki*ab(2)
      ab(3)=4*l1
      ab(2)=ab(2)*ab(3)
      ab(4)= - 29.0_ki/27.0_ki - zeta3
      ab(4)=16*ab(4) + 11.0_ki/3.0_ki*ab(1)
      ab(4)=lam*ab(4)
      ab(5)= - 11.0_ki/3.0_ki + 4*l2
      ab(5)=l2*ab(5)
      ab(5)=ab(5) + zeta3
      ab(5)=163.0_ki/27.0_ki + 2*ab(5)
      ab(2)=ab(2) + ab(4) + 8*ab(5) - 35.0_ki/3.0_ki*ab(1)
      ab(2)=ab(2)*ab(3)
      ab(3)=32*plg4half
      ab(4)= - 1109.0_ki/81.0_ki + ab(3)
      ab(5)= - 11 + 8*l2
      ab(5)=l2*ab(5)
      ab(5)=35.0_ki/3.0_ki + ab(5)
      ab(5)=l2*ab(5)
      ab(4)=8.0_ki/3.0_ki*ab(5) + 2*ab(4) + 7.0_ki/3.0_ki*zeta3
      ab(5)=1.0_ki/9.0_ki*ab(1)
      ab(1)=163.0_ki/5.0_ki*ab(1)
      ab(6)=202 - ab(1)
      ab(6)=ab(6)*ab(5)
      ab(1)= - 70 + ab(1)
      ab(1)=ab(1)*ab(5)
      ab(3)=413.0_ki/81.0_ki - ab(3)
      ab(3)=2*ab(3) - 55.0_ki/3.0_ki*zeta3
      ab(1)=4*ab(3) + ab(1)
      ab(1)=lam*ab(1)
      ab(1)=ab(2) + ab(1) + 4*ab(4) + ab(6)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=z1*z2
      ab(2)=lam**2
      ab(2)= - 1 + ab(2)

      tmp = 16*CA*Ta2*lam**3*ab(2)*ab(1)**2
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l4 + l2 - l3
      ab(1)=11.0_ki/3.0_ki - 8*ab(1)
      ab(2)=lam**2
      ab(2)=ab(2) - 1
      ab(3)=z1*z2

      tmp = 8*CA*Ta2*lam**3*ab(3)**2*ab(2)*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l4 + l2
      ab(1)=55.0_ki/3.0_ki + 16*l3 - 36*ab(1)
      ab(2)=2*l3
      ab(1)=ab(1)*ab(2)
      ab(2)=8*l2
      ab(3)=ab(2) - 11.0_ki/3.0_ki + 4*l4
      ab(4)=8*l4
      ab(3)=ab(3)*ab(4)
      ab(4)= - 11.0_ki/3.0_ki + 4*l2
      ab(2)=ab(4)*ab(2)
      ab(4)=d3 + d2
      ab(1)=ab(3) + ab(1) + ab(2) - 8*ab(4)
      ab(2)=Pi**2
      ab(2)= - 137.0_ki/3.0_ki + 4*ab(2)
      ab(3)= - 1.0_ki/3.0_ki*ab(2) + ab(1)
      ab(3)=lam*ab(3)
      ab(4)=l4 - l2
      ab(4)=ab(4)*l3
      ab(4)=ab(4) - d3 + d2
      ab(4)=8*ab(4)
      ab(3)= - ab(4) + ab(3)
      ab(3)=lam*ab(3)
      ab(5)=l5*r1
      ab(2)=ab(5) + ab(2)
      ab(1)=ab(3) + 1.0_ki/3.0_ki*ab(2) - ab(1)
      ab(2)=z2*z1
      ab(2)=Ta2*CA*ab(2)**2
      ab(1)=lam*ab(2)*ab(1)
      ab(2)=ab(2)*ab(4)
      ab(1)=ab(2) + ab(1)

      tmp = 4*lam**2*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=Pi**2
      ab(2)=4*ab(1)
      ab(3)=d3 + d2
      ab(3)=8*ab(3)
      ab(4)=ab(3) + ab(2) - 137.0_ki/9.0_ki
      ab(5)= - 11.0_ki/3.0_ki + 4*l2
      ab(6)=8*l2
      ab(7)= - ab(5)*ab(6)
      ab(8)=ab(6) - 11.0_ki/3.0_ki
      ab(9)=ab(8) + 3*l4
      ab(10)=4*l4
      ab(11)= - ab(9)*ab(10)
      ab(7)=ab(11) + ab(7) + ab(4)
      ab(11)=2*l4
      ab(7)=ab(7)*ab(11)
      ab(12)=ab(6) - 11
      ab(13)=4.0_ki/3.0_ki*l2
      ab(12)=ab(12)*ab(13)
      ab(2)=ab(2) - 137.0_ki/3.0_ki
      ab(2)=ab(3) - ab(12) + 1.0_ki/3.0_ki*ab(2)
      ab(12)=2*l2
      ab(13)=ab(2)*ab(12)
      ab(8)=ab(10) + ab(8)
      ab(8)=l4*ab(8)
      ab(5)=ab(5)*l2
      ab(8)=ab(5) + ab(8)
      ab(10)= - 347.0_ki/3.0_ki + 8*ab(1)
      ab(3)=ab(3) + 1.0_ki/3.0_ki*ab(10)
      ab(10)=l4 + l2
      ab(14)=64.0_ki/3.0_ki*l3 + 121.0_ki/3.0_ki - 76*ab(10)
      ab(14)=l3*ab(14)
      ab(8)=ab(14) - ab(3) + 20*ab(8)
      ab(8)=l3*ab(8)
      ab(14)= - t2 + 2*t4
      ab(15)= - t3 + 2*t1
      ab(16)=ab(15) - 3*zeta3 + ab(14) + 52.0_ki/27.0_ki
      ab(16)= - 11.0_ki/3.0_ki*ab(1) + 8*ab(16)
      ab(7)=ab(8) + ab(7) + ab(13) + ab(16)
      ab(7)=lam*ab(7)
      ab(8)=d3 - d2
      ab(13)= - ab(8) + 1.0_ki/3.0_ki*ab(1)
      ab(17)=l4**2
      ab(13)= - 1.0_ki/3.0_ki*ab(17) + 2*ab(13)
      ab(13)=ab(13)*l4
      ab(18)=ab(12)*ab(8)
      ab(13)=ab(13) - ab(15) - ab(18) + ab(14)
      ab(14)=l2**2
      ab(14)=ab(14) - ab(17)
      ab(8)=ab(8) - 2*ab(14)
      ab(14)=l4 - l2
      ab(15)=3*l3
      ab(14)=ab(14)*ab(15)
      ab(8)= - ab(14) + 2*ab(8)
      ab(8)=ab(8)*l3
      ab(8)=ab(8) + 2*ab(13)
      ab(7)=4*ab(8) + ab(7)
      ab(7)=lam*ab(7)
      ab(13)=1.0_ki/3.0_ki*r1
      ab(14)= - 2*l1 - 5*l5
      ab(14)=ab(14)*ab(13)
      ab(9)=l4*ab(9)
      ab(4)=16*ab(9) + 32*ab(5) - 4*ab(4) + ab(14)
      ab(4)=ab(4)*ab(11)
      ab(5)=220 - r1
      ab(5)=1.0_ki/3.0_ki*ab(5) - 80*l2
      ab(5)=ab(5)*ab(12)
      ab(9)=1.0_ki/3.0_ki*l1 + l5
      ab(9)=r1*ab(9)
      ab(11)= - 160*l4 - 320*l2 + 440.0_ki/3.0_ki - r1
      ab(11)=l4*ab(11)
      ab(3)=ab(11) + ab(5) + 2*ab(3) + ab(9)
      ab(5)= - 484 - r1
      ab(5)= - 256.0_ki/3.0_ki*l3 + 1.0_ki/3.0_ki*ab(5) + 304*ab(10)
      ab(5)=l3*ab(5)
      ab(3)=2*ab(3) + ab(5)
      ab(3)=l3*ab(3)
      ab(5)=1.0_ki/3.0_ki*l5
      ab(9)=20.0_ki/3.0_ki - l1
      ab(9)=2*ab(9) + l5
      ab(5)=ab(9)*ab(5)
      ab(9)=d1 + 4*d4
      ab(9)= - 4.0_ki/3.0_ki*d6 + 1.0_ki/3.0_ki*ab(9) - 2*d5
      ab(1)= - 4.0_ki/3.0_ki*d3 - 8.0_ki/3.0_ki*d2 + 2.0_ki/9.0_ki*ab(1) + 2*ab(9)
     &  + ab(5)
      ab(1)=r1*ab(1)
      ab(5)= - l5*ab(13)
      ab(2)=ab(5) - ab(2)
      ab(2)=ab(2)*ab(6)
      ab(1)=4*ab(7) + ab(3) + ab(4) + ab(2) - 4*ab(16) + ab(1)
      ab(2)=z2*z1
      ab(2)=Ta2*CA*ab(2)**2
      ab(1)=lam*ab(2)*ab(1)
      ab(2)=ab(8)*ab(2)
      ab(1)= - 16*ab(2) + ab(1)

      tmp = 2*lam**2*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

