      subroutine a2srsgg(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(13)
      include 'a2srs_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)=CA + 4*Ta2
      ab(2)=lam - 1

      tmp = z1*Ta2*ab(2)*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=lam - 1
      ab(1)=ab(1)*l1
      ab(2)= - 35 + 11*lam
      ab(2)=1.0_ki/6.0_ki*ab(2) - 4*ab(1)
      ab(2)=CA*ab(2)
      ab(1)= - 1 - ab(1)
      ab(1)=Ta2*ab(1)
      ab(1)=ab(2) + 16*ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=ab(1)*l1
      ab(3)=35 - 11*lam
      ab(3)=1.0_ki/3.0_ki*ab(3) + 4*ab(2)
      ab(3)=l1*ab(3)
      ab(4)=16*l2
      ab(5)=Pi**2
      ab(6)=67.0_ki/3.0_ki - 7*ab(5)
      ab(6)=lam*ab(6)
      ab(3)=2*ab(3) + 1.0_ki/6.0_ki*ab(6) + 7.0_ki/6.0_ki*ab(5) - 199.0_ki/18.0_ki
     &  + ab(4)
      ab(3)=CA*ab(3)
      ab(1)= - ab(5)*ab(1)
      ab(2)=2 + ab(2)
      ab(2)=l1*ab(2)
      ab(1)=8*ab(2) + ab(4) + ab(1)
      ab(1)=Ta2*ab(1)
      ab(1)=ab(3) + 4*ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=8*l1
      ab(3)= - ab(1)*ab(2)
      ab(3)=ab(3) - 35 + 11*lam
      ab(3)=l1*ab(3)
      ab(4)=Pi**2
      ab(5)= - 67.0_ki/3.0_ki + 7*ab(4)
      ab(5)=lam*ab(5)
      ab(3)=2.0_ki/3.0_ki*ab(3) + 1.0_ki/3.0_ki*ab(5) - 7.0_ki/3.0_ki*ab(4) + 199.0_ki
     & /9.0_ki - 32*l2
      ab(3)=l1*ab(3)
      ab(5)=25*zeta3
      ab(6)= - 11.0_ki/6.0_ki*ab(4) + 202.0_ki/27.0_ki - ab(5)
      ab(6)=lam*ab(6)
      ab(7)=11.0_ki/3.0_ki - 4*l2
      ab(7)=l2*ab(7)
      ab(3)=2*ab(3) + ab(6) + 13.0_ki/2.0_ki*ab(4) + 8*ab(7) - 604.0_ki/27.0_ki
     &  + ab(5)
      ab(3)=CA*ab(3)
      ab(5)=ab(4)*ab(1)
      ab(6)=l1*ab(1)
      ab(6)= - 1 - 1.0_ki/3.0_ki*ab(6)
      ab(2)=ab(6)*ab(2)
      ab(2)=ab(2) - 16*l2 + ab(5)
      ab(2)=l1*ab(2)
      ab(5)=5*zeta3
      ab(1)= - ab(5)*ab(1)
      ab(5)=l2**2
      ab(1)=ab(2) + ab(4) - 8*ab(5) + ab(1)
      ab(1)=Ta2*ab(1)
      ab(1)=ab(3) + 16*ab(1)

      tmp = z1*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=Ta2**2
      ab(2)=ab(1)*z1
      ab(3)=CA*z1
      ab(4)=ab(3)*Ta2
      ab(5)=8*ab(2) + 7.0_ki/3.0_ki*ab(4)
      ab(6)=2*l1
      ab(7)=ab(5)*ab(6)
      ab(8)=32*ab(2)
      ab(9)=ab(7) - ab(8) - 13*ab(4)
      ab(9)=ab(9)*ab(6)
      ab(10)=11.0_ki/3.0_ki*Ta2
      ab(3)=ab(10)*ab(3)
      ab(7)=ab(3) - ab(7)
      ab(7)=ab(7)*ab(6)
      ab(7)= - 67.0_ki/18.0_ki*ab(4) + ab(7)
      ab(7)=lam*ab(7)
      ab(10)=4*z1
      ab(1)=ab(10)*ab(1)
      ab(10)=ab(1) + 37.0_ki/36.0_ki*ab(4)
      ab(11)= - lam + 1
      ab(12)=Pi**2
      ab(10)=ab(12)*ab(10)*ab(11)
      ab(7)=1.0_ki/5.0_ki*ab(10) + ab(7) + 199.0_ki/18.0_ki*ab(4) + ab(9)
      ab(7)=ab(7)*ab(12)
      ab(1)=ab(1) + ab(4)
      ab(9)=ab(1)*ab(6)
      ab(10)= - ab(9) + ab(8) + 35.0_ki/3.0_ki*ab(4)
      ab(11)=4*l1
      ab(10)=ab(10)*ab(11)
      ab(10)= - 199.0_ki/3.0_ki*ab(4) + ab(10)
      ab(10)=l1*ab(10)
      ab(10)=604.0_ki/9.0_ki*ab(4) + ab(10)
      ab(10)=ab(10)*ab(6)
      ab(9)= - ab(3) + ab(9)
      ab(9)=ab(9)*ab(11)
      ab(9)=67.0_ki/3.0_ki*ab(4) + ab(9)
      ab(9)=l1*ab(9)
      ab(9)= - 202.0_ki/9.0_ki*ab(4) + ab(9)
      ab(9)=ab(9)*ab(6)
      ab(9)=607.0_ki/27.0_ki*ab(4) + ab(9)
      ab(9)=lam*ab(9)
      ab(9)=ab(9) - 1819.0_ki/27.0_ki*ab(4) + ab(10)
      ab(10)=ab(1)*ab(11)
      ab(10)= - ab(3) + ab(10)
      ab(10)=ab(10)*ab(11)
      ab(11)=8.0_ki/3.0_ki*l2 + 8*l1
      ab(1)=ab(1)*ab(11)
      ab(1)= - ab(3) + ab(1)
      ab(1)=l2*ab(1)
      ab(5)= - ab(5)*ab(12)
      ab(1)=2*ab(1) + ab(5) + 67.0_ki/9.0_ki*ab(4) + ab(10)
      ab(1)=l2*ab(1)
      ab(2)=16*ab(2) + 5*ab(4)
      ab(2)=ab(2)*ab(6)
      ab(3)= - ab(3) + ab(2)
      ab(3)=lam*ab(3)
      ab(2)=ab(3) - ab(2) + ab(8) + 41.0_ki/3.0_ki*ab(4)
      ab(2)=zeta3*ab(2)

      tmp = 8*ab(1) + 10*ab(2) + ab(7) + 2.0_ki/3.0_ki*ab(9)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=CA + 4*Ta2

      tmp = 8*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l4 + l2
      ab(2)=ab(1) - l3
      ab(2)=11.0_ki/3.0_ki - 8*ab(2)
      ab(2)=CA*ab(2)
      ab(1)=5*l3 - 4*ab(1)
      ab(1)=Ta2*ab(1)
      ab(1)=ab(2) + 8*ab(1)

      tmp = 4*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=Ta2**2
      ab(2)=4*ab(1)
      ab(3)=lam**3*z1*z2
      ab(2)=ab(2)*ab(3)
      ab(4)=CA*Ta2
      ab(5)=ab(3)*ab(4)
      ab(2)=ab(2) + ab(5)
      ab(6)=8*l4
      ab(7)=4*l2
      ab(8)=ab(7) + ab(6)
      ab(8)=ab(2)*ab(8)
      ab(9)=5*ab(1)
      ab(9)=ab(9)*ab(3)
      ab(9)=ab(9) + ab(5)
      ab(10)=l3*ab(9)
      ab(4)=11.0_ki/3.0_ki*ab(4)
      ab(4)=ab(4)*ab(3)
      ab(8)= - 8*ab(10) - ab(4) + ab(8)
      ab(7)=ab(8)*ab(7)
      ab(8)=4*l4
      ab(2)=ab(2)*ab(8)
      ab(2)= - ab(4) + ab(2)
      ab(2)=ab(2)*ab(8)
      ab(8)=4*l3
      ab(6)=ab(8) - ab(6)
      ab(6)=ab(9)*ab(6)
      ab(4)=ab(4) + ab(6)
      ab(4)=ab(4)*ab(8)
      ab(1)=ab(1)*ab(3)
      ab(1)= - 8*ab(1) - 7.0_ki/3.0_ki*ab(5)
      ab(1)=ab(1)*Pi**2
      ab(1)=ab(1) + ab(7) + ab(4) + 67.0_ki/9.0_ki*ab(5) + ab(2)

      tmp = 4*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=2*l2
      ab(2)= - 16*l2 + 44.0_ki/3.0_ki - l1
      ab(2)=ab(2)*ab(1)
      ab(3)=l1**2
      ab(3)=ab(3) - 134.0_ki/9.0_ki
      ab(4)=Pi**2
      ab(5)= - 31.0_ki/3.0_ki*l4 + 44.0_ki/3.0_ki - 31*l2
      ab(5)=l4*ab(5)
      ab(2)=ab(5) + ab(2) + 4*ab(4) + ab(3)
      ab(5)=4*l4
      ab(2)=ab(2)*ab(5)
      ab(6)=d2 + d1
      ab(7)= - 7.0_ki/3.0_ki*ab(4) + 67.0_ki/9.0_ki + 2*ab(6)
      ab(8)=2*l3
      ab(9)=l4 + l2
      ab(10)=16.0_ki/3.0_ki*l3 + 22.0_ki/3.0_ki - 15*ab(9)
      ab(10)=ab(10)*ab(8)
      ab(11)= - 33 + 32*l2
      ab(11)=l2*ab(11)
      ab(12)=32*l4 - 33 + 64*l2
      ab(12)=l4*ab(12)
      ab(7)=ab(10) + ab(12) + 2*ab(7) + ab(11)
      ab(7)=l3*ab(7)
      ab(10)=t4 + t3
      ab(11)=ab(10) + t5 + t2
      ab(11)=t7 + t6 + 101.0_ki/27.0_ki + t1 - 2*ab(11)
      ab(12)=l1**3
      ab(11)= - 4.0_ki/3.0_ki*ab(12) + 4*ab(11) - 45*zeta3 - 22.0_ki/3.0_ki*
     & ab(6)
      ab(12)=11 - 8*l2
      ab(12)=l2*ab(12)
      ab(3)=4.0_ki/3.0_ki*ab(12) + 14.0_ki/3.0_ki*ab(4) + ab(3)
      ab(12)=4*l2
      ab(3)=ab(3)*ab(12)
      ab(13)= - 55.0_ki/3.0_ki + 4*l1
      ab(13)=ab(13)*ab(4)
      ab(2)=4*ab(7) + ab(2) + ab(3) + 2*ab(11) + 1.0_ki/3.0_ki*ab(13)
      ab(2)=CA*ab(2)
      ab(3)=l2**2
      ab(1)=ab(1) + l4
      ab(1)=l4*ab(1)
      ab(1)=ab(3) + ab(1)
      ab(7)= - 21*ab(9) + 20.0_ki/3.0_ki*l3
      ab(7)=l3*ab(7)
      ab(1)=ab(7) + 8*ab(6) - 5.0_ki/3.0_ki*ab(4) + 20*ab(1)
      ab(1)=ab(1)*ab(8)
      ab(6)=ab(4) - 8.0_ki/3.0_ki*ab(3)
      ab(6)=ab(6)*ab(12)
      ab(7)= - l2 - 1.0_ki/3.0_ki*l4
      ab(7)=l4*ab(7)
      ab(3)=ab(3) - ab(7)
      ab(3)=ab(4) - 8*ab(3)
      ab(3)=ab(3)*ab(5)
      ab(1)=ab(1) + ab(3) + ab(6) - 36*ab(10) - 11*zeta3
      ab(1)=Ta2*ab(1)
      ab(1)=ab(2) + 16*ab(1)

      tmp = 2*z1*z2*Ta2*lam**3*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

