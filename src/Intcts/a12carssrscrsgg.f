      subroutine a12carssrscrsgg(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(16)
      include 'a12carssrscrs_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)= - 1 + lam

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=4*l1
      ab(2)=11.0_ki/6.0_ki - ab(1)
      ab(2)=lam*ab(2)
      ab(1)=ab(2) - 35.0_ki/6.0_ki + ab(1)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=35 - 11*lam
      ab(2)= - 1 + lam
      ab(2)=l1*ab(2)
      ab(1)=1.0_ki/3.0_ki*ab(1) + 4*ab(2)
      ab(1)=l1*ab(1)
      ab(1)=2*ab(1) + 67.0_ki/18.0_ki*lam - 199.0_ki/18.0_ki + 16*l2

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=Pi**2
      ab(2)=2*ab(1)
      ab(3)= - 67.0_ki/9.0_ki + ab(2)
      ab(3)=lam*ab(3)
      ab(4)=1 - lam
      ab(4)=l1*ab(4)
      ab(4)=8*ab(4) - 35 + 11*lam
      ab(4)=l1*ab(4)
      ab(2)=2.0_ki/3.0_ki*ab(4) + ab(3) - ab(2) + 199.0_ki/9.0_ki - 32*l2
      ab(2)=l1*ab(2)
      ab(3)=4*l2
      ab(4)=11.0_ki/3.0_ki - ab(3)
      ab(3)=ab(4)*ab(3)
      ab(2)=ab(2) + ab(3) - 311.0_ki/27.0_ki - 4*zeta3
      ab(1)=11.0_ki/6.0_ki*ab(1)
      ab(3)=55.0_ki/27.0_ki + 2*zeta3
      ab(3)=4*ab(3) - ab(1)
      ab(3)=lam*ab(3)
      ab(1)=ab(3) + ab(1) + 2*ab(2)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=35 - 11*lam
      ab(2)= - 1 + lam
      ab(2)=l1*ab(2)
      ab(1)=1.0_ki/3.0_ki*ab(1) + ab(2)
      ab(1)=l1*ab(1)
      ab(2)=Pi**2
      ab(3)=67.0_ki/3.0_ki - 2*ab(2)
      ab(3)=lam*ab(3)
      ab(1)=4.0_ki/3.0_ki*ab(1) + 1.0_ki/3.0_ki*ab(3) + 2.0_ki/3.0_ki*ab(2) - 199.0_ki
     & /9.0_ki + 32*l2
      ab(3)=2*l1
      ab(1)=ab(1)*ab(3)
      ab(4)=4*l2
      ab(5)= - 11.0_ki/3.0_ki + ab(4)
      ab(4)=ab(5)*ab(4)
      ab(4)=ab(4) + 311.0_ki/27.0_ki + 4*zeta3
      ab(5)= - 55.0_ki/27.0_ki - 2*zeta3
      ab(5)=8*ab(5) + 11.0_ki/3.0_ki*ab(2)
      ab(5)=lam*ab(5)
      ab(1)=ab(1) + ab(5) + 4*ab(4) - 35.0_ki/3.0_ki*ab(2)
      ab(1)=ab(1)*ab(3)
      ab(3)=32*plg4half
      ab(4)= - 1067.0_ki/81.0_ki + ab(3)
      ab(5)= - 11 + 8*l2
      ab(5)=l2*ab(5)
      ab(5)=67.0_ki/3.0_ki + 2*ab(5)
      ab(5)=l2*ab(5)
      ab(4)=4.0_ki/3.0_ki*ab(5) + 2*ab(4) + 7.0_ki/3.0_ki*zeta3
      ab(5)=1.0_ki/18.0_ki*ab(2)
      ab(2)=163.0_ki/5.0_ki*ab(2)
      ab(6)=199 - ab(2)
      ab(6)=ab(6)*ab(5)
      ab(2)= - 67 + ab(2)
      ab(2)=ab(2)*ab(5)
      ab(3)=407.0_ki/81.0_ki - ab(3)
      ab(3)=2*ab(3) - 55.0_ki/3.0_ki*zeta3
      ab(2)=2*ab(3) + ab(2)
      ab(2)=lam*ab(2)
      ab(1)=ab(1) + ab(2) + 2*ab(4) + ab(6)

      tmp = CA*z1*Ta2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 8*CA*z1*z2*Ta2*lam**3
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=1/(1 - lam**2)
      ab(2)=lam**2
      ab(2)=1 - ab(2)
      ab(2)=l3*ab(1)*ab(2)
      ab(2)= - l4 + ab(2) - l2
      ab(2)=11.0_ki/3.0_ki + 8*ab(2)

      tmp = 4*CA*z1*z2*Ta2*lam**3*ab(2)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=1/(1 - lam**2)
      ab(2)=1/(3 - 3*lam**2)
      ab(3)=1/((1 - lam**2)**(3.0_ki/2.0_ki))
      ab(4)=1.0_ki/3.0_ki*ab(3)
      ab(5)=ab(4)*l5
      ab(6)=ab(2)*Pi**2
      ab(5)=ab(5) + 4*ab(6)
      ab(6)=32*ab(1)
      ab(7)=ab(6)*l3
      ab(4)=ab(4) + ab(7) + 110*ab(2)
      ab(7)= - l3*ab(4)
      ab(6)= - l4**2*ab(6)
      ab(6)=ab(7) + ab(6) + ab(5)
      ab(7)=lam**2
      ab(6)=ab(6)*ab(7)
      ab(8)=8*l2
      ab(9)=l4*ab(1)
      ab(9)=4*ab(9) - 11.0_ki/3.0_ki + ab(8)
      ab(9)=l4*ab(9)
      ab(9)= - d1 + ab(9) - d5
      ab(10)= - 11.0_ki/3.0_ki + 4*l2
      ab(8)=ab(10)*ab(8)
      ab(10)=l4 + l2
      ab(4)=ab(4) - 72*ab(10)
      ab(4)=l3*ab(4)
      ab(4)=ab(6) + ab(4) + ab(8) + 131.0_ki/9.0_ki - ab(5) + 8*ab(9)
      ab(4)=lam*ab(4)
      ab(5)=l2 - l4
      ab(5)=l3*ab(5)
      ab(5)=ab(5) + d1 - d5
      ab(4)=8*ab(5) + ab(4)

      tmp = 2*CA*z1*z2*Ta2*ab(7)*ab(4)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=(1 - lam**2)**( - 3.0_ki/2.0_ki)
      ab(2)=1/(3 - 3*lam**2)
      ab(3)=1/(9 - 9*lam**2)
      ab(4)=1/((1 - lam**2)**(3.0_ki/2.0_ki))
      ab(5)=d1 + d5
      ab(5)= - 131.0_ki/9.0_ki + 8*ab(5)
      ab(6)=Pi**2
      ab(7)=4*ab(6) + ab(5)
      ab(8)=l5*ab(1)
      ab(9)=4*l2
      ab(10)=11.0_ki/3.0_ki - ab(9)
      ab(10)=l2*ab(10)
      ab(11)=8*l2
      ab(12)= - 3*l4 + 11*ab(2) - ab(11)
      ab(12)=l4*ab(12)
      ab(7)=16*ab(12) + 32*ab(10) + 4*ab(7) + ab(8)
      ab(7)=l4*ab(7)
      ab(10)=l1**2
      ab(12)=1.0_ki/3.0_ki - ab(2)
      ab(12)=ab(12)*ab(10)
      ab(13)=11 - ab(11)
      ab(13)=l2*ab(13)
      ab(13)=ab(6) + ab(13)
      ab(5)=ab(5) + 4.0_ki/3.0_ki*ab(13)
      ab(5)=ab(5)*ab(9)
      ab(9)=t3 + t2
      ab(13)=t4 + t1
      ab(14)=ab(6)*ab(2)
      ab(5)=ab(7) + ab(5) + 8*ab(12) - 22*ab(14) + 775.0_ki/27.0_ki - 48*
     & zeta3 + 32*ab(13) - 16*ab(9)
      ab(7)=l1 - 23.0_ki/3.0_ki
      ab(7)=ab(11) + l5 + 2*ab(7)
      ab(7)=l5*ab(7)
      ab(9)=2*d1
      ab(12)=1.0_ki/3.0_ki*ab(6)
      ab(13)=2*d3 + ab(9) - ab(12) + d2 + 4*d4
      ab(15)= - l4 + 2*l1
      ab(15)=ab(15)*l4
      ab(7)= - 2*ab(13) + ab(15) - ab(10) + ab(7)
      ab(13)=1.0_ki/3.0_ki*ab(4)
      ab(7)=ab(7)*ab(13)
      ab(10)=ab(2)*ab(10)
      ab(10)=11*ab(14) + 4*ab(10)
      ab(14)=l4*ab(2)
      ab(8)= - ab(8) - 176*ab(14)
      ab(8)=l4*ab(8)
      ab(8)=2*ab(10) + ab(8)
      ab(10)=4*l4
      ab(14)=4*l5 + ab(10) - 23.0_ki/3.0_ki
      ab(13)=ab(14)*ab(13)
      ab(14)=l2*ab(1)
      ab(14)= - 332*ab(3) + ab(14)
      ab(15)=ab(1) + 242*ab(2)
      ab(16)= - l3*ab(15)
      ab(14)=ab(16) + 2*ab(14) + ab(13)
      ab(16)=2*l3
      ab(14)=ab(14)*ab(16)
      ab(8)=ab(14) + 2*ab(8) - ab(7)
      ab(14)=lam**2
      ab(8)=ab(8)*ab(14)
      ab(6)= - ab(9) - 2.0_ki/3.0_ki*ab(6) + 83*ab(3) - 2*d5
      ab(9)=ab(10) - 11.0_ki/3.0_ki + ab(11)
      ab(9)=l4*ab(9)
      ab(10)=80*l2 - 220.0_ki/3.0_ki - ab(1)
      ab(10)=l2*ab(10)
      ab(6)=20*ab(9) + 4*ab(6) + ab(10)
      ab(9)=l4 + l2
      ab(9)=128.0_ki/3.0_ki*l3 + ab(15) - 152*ab(9)
      ab(9)=l3*ab(9)
      ab(6)=ab(9) + 2*ab(6) - ab(13)
      ab(6)=ab(6)*ab(16)
      ab(5)=ab(8) + ab(6) + 2*ab(5) + ab(7)
      ab(5)=lam*ab(5)
      ab(6)=d1 - d5
      ab(7)=ab(12) - ab(6)
      ab(8)=l4**2
      ab(7)=2*ab(7) - 1.0_ki/3.0_ki*ab(8)
      ab(7)=l4*ab(7)
      ab(9)=l2*ab(6)
      ab(9)=t1 + ab(9) - t4
      ab(7)=ab(7) + t3 - t2 - 2*ab(9)
      ab(9)=l2**2
      ab(8)=ab(8) - ab(9)
      ab(6)=ab(6) + 2*ab(8)
      ab(8)=l2 - l4
      ab(8)=l3*ab(8)
      ab(6)=2*ab(6) + 3*ab(8)
      ab(6)=l3*ab(6)
      ab(6)=2*ab(7) + ab(6)
      ab(5)=16*ab(6) + ab(5)

      tmp = CA*z1*z2*Ta2*ab(14)*ab(5)
      res(0,0,0) = real(tmp,ki)

      return
      end

