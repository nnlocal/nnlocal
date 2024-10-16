      subroutine hp_a2carssrsgg(Ta2,lam,res)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) tmp,hp_cli2,hp_Li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(16)
      include 'hp_a2carssrs_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)=1.0_ki/2.0_ki*CA + Ta2
      ab(2)=lam - 1

      tmp = z2*Ta2*ab(2)*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)=lam - 1
      ab(1)=ab(1)*l1
      ab(2)= - 35 + 11*lam
      ab(2)=1.0_ki/12.0_ki*ab(2) - 2*ab(1)
      ab(2)=CA*ab(2)
      ab(1)= - 1 - ab(1)
      ab(1)=Ta2*ab(1)
      ab(1)=ab(2) + 4*ab(1)

      tmp = z2*Ta2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=ab(1)*l1
      ab(3)=35 - 11*lam
      ab(3)=1.0_ki/3.0_ki*ab(3) + 4*ab(2)
      ab(3)=l1*ab(3)
      ab(4)=Pi**2
      ab(5)=61.0_ki/3.0_ki - 5*ab(4)
      ab(5)=lam*ab(5)
      ab(3)=ab(3) + 1.0_ki/12.0_ki*ab(5) + 5.0_ki/12.0_ki*ab(4) - 193.0_ki/36.0_ki
     &  + 8*l2
      ab(3)=CA*ab(3)
      ab(1)= - ab(4)*ab(1)
      ab(2)=2 + ab(2)
      ab(2)=l1*ab(2)
      ab(1)=8*ab(2) + 16*l2 + ab(1)
      ab(1)=Ta2*ab(1)
      ab(1)=ab(3) + ab(1)

      tmp = z2*Ta2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=8*l1
      ab(3)= - ab(1)*ab(2)
      ab(3)=ab(3) - 35 + 11*lam
      ab(3)=l1*ab(3)
      ab(4)=Pi**2
      ab(5)= - 61.0_ki/3.0_ki + 5*ab(4)
      ab(5)=lam*ab(5)
      ab(3)=2.0_ki/3.0_ki*ab(3) + 1.0_ki/3.0_ki*ab(5) - 5.0_ki/3.0_ki*ab(4) + 193.0_ki
     & /9.0_ki - 32*l2
      ab(3)=l1*ab(3)
      ab(5)=13*zeta3
      ab(6)= - 11.0_ki/6.0_ki*ab(4) + 151.0_ki/27.0_ki - ab(5)
      ab(6)=lam*ab(6)
      ab(5)=ab(6) - 517.0_ki/27.0_ki + ab(5)
      ab(6)=4*l2
      ab(7)=11.0_ki/3.0_ki - ab(6)
      ab(6)=ab(7)*ab(6)
      ab(3)=ab(3) + 31.0_ki/12.0_ki*ab(4) + ab(6) + 1.0_ki/2.0_ki*ab(5)
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
      ab(1)=ab(3) + 4*ab(1)

      tmp = z2*Ta2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=Ta2**2
      ab(2)=ab(1)*z2
      ab(3)=CA*z2
      ab(4)=ab(3)*Ta2
      ab(5)=4*ab(2) + 5.0_ki/3.0_ki*ab(4)
      ab(6)=2*l1
      ab(7)=ab(5)*ab(6)
      ab(8)=16*ab(2)
      ab(9)=ab(7) - ab(8) - 31.0_ki/3.0_ki*ab(4)
      ab(9)=l1*ab(9)
      ab(10)=11.0_ki/3.0_ki*Ta2
      ab(3)=ab(10)*ab(3)
      ab(7)=ab(3) - ab(7)
      ab(7)=l1*ab(7)
      ab(7)= - 61.0_ki/36.0_ki*ab(4) + ab(7)
      ab(7)=lam*ab(7)
      ab(10)=ab(2) + 1.0_ki/12.0_ki*ab(4)
      ab(11)= - lam + 1
      ab(12)=Pi**2
      ab(10)=ab(12)*ab(10)*ab(11)
      ab(7)=1.0_ki/5.0_ki*ab(10) + ab(7) + 193.0_ki/36.0_ki*ab(4) + ab(9)
      ab(7)=ab(7)*ab(12)
      ab(9)=2*z2
      ab(1)=ab(9)*ab(1)
      ab(1)=ab(1) + ab(4)
      ab(9)=ab(1)*ab(6)
      ab(8)= - ab(9) + ab(8) + 35.0_ki/3.0_ki*ab(4)
      ab(10)=4*l1
      ab(8)=ab(8)*ab(10)
      ab(8)= - 193.0_ki/3.0_ki*ab(4) + ab(8)
      ab(8)=l1*ab(8)
      ab(8)=517.0_ki/9.0_ki*ab(4) + ab(8)
      ab(8)=ab(8)*ab(6)
      ab(9)= - ab(3) + ab(9)
      ab(9)=ab(9)*ab(10)
      ab(9)=61.0_ki/3.0_ki*ab(4) + ab(9)
      ab(9)=l1*ab(9)
      ab(9)= - 151.0_ki/9.0_ki*ab(4) + ab(9)
      ab(9)=ab(9)*ab(6)
      ab(9)=259.0_ki/27.0_ki*ab(4) + ab(9)
      ab(9)=lam*ab(9)
      ab(8)=ab(9) - 1165.0_ki/27.0_ki*ab(4) + ab(8)
      ab(9)=ab(1)*ab(10)
      ab(9)= - ab(3) + ab(9)
      ab(9)=ab(9)*ab(10)
      ab(10)=8.0_ki/3.0_ki*l2 + 8*l1
      ab(1)=ab(1)*ab(10)
      ab(1)= - ab(3) + ab(1)
      ab(1)=l2*ab(1)
      ab(3)= - ab(5)*ab(12)
      ab(1)=2*ab(1) + ab(3) + 61.0_ki/9.0_ki*ab(4) + ab(9)
      ab(1)=l2*ab(1)
      ab(3)=40*ab(2) + 13*ab(4)
      ab(3)=ab(3)*ab(6)
      ab(5)= - 55.0_ki/3.0_ki*ab(4) + ab(3)
      ab(5)=lam*ab(5)
      ab(2)=ab(5) - ab(3) + 80*ab(2) + 133.0_ki/3.0_ki*ab(4)
      ab(2)=zeta3*ab(2)

      tmp = 4*ab(1) + ab(2) + ab(7) + 1.0_ki/3.0_ki*ab(8)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - lam**2
      ab(1)=1 + ab(1)
      ab(2)=CA + 2*Ta2

      tmp = 4*z1**2*z2*z3*Ta2*lam**3*ab(2)*ab(1)
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l2 + l4 - l3
      ab(2)=16*Ta2
      ab(2)=ab(1)*ab(2)
      ab(1)= - 11.0_ki/3.0_ki + 8*ab(1)
      ab(1)=ab(1)*CA
      ab(1)=ab(1) + ab(2)
      ab(2)=lam**2
      ab(2)= - 1 + ab(2)

      tmp = 2*z1**2*z2*z3*Ta2*lam**3*ab(2)*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam**2
      ab(1)=ab(1) - 1
      ab(1)=z1**2*ab(1)*lam**3*z2*z3
      ab(2)=Ta2**2*ab(1)
      ab(1)=CA*Ta2*ab(1)
      ab(3)=5*ab(2) + 2*ab(1)
      ab(4)=4*l4
      ab(5)=ab(3)*ab(4)
      ab(6)=ab(1) + 2*ab(2)
      ab(7)=4*l3
      ab(8)= - ab(6)*ab(7)
      ab(9)=11.0_ki/3.0_ki*ab(1)
      ab(5)=ab(8) - ab(9) + ab(5)
      ab(5)=ab(5)*ab(7)
      ab(8)=4*l2
      ab(10)= - ab(8) - 8*l4
      ab(10)=ab(6)*ab(10)
      ab(3)=ab(3)*ab(7)
      ab(3)=ab(3) + ab(9) + ab(10)
      ab(3)=ab(3)*ab(8)
      ab(7)=d1 + d2
      ab(7)=16*ab(7)
      ab(7)=ab(2)*ab(7)
      ab(6)= - ab(6)*ab(4)
      ab(6)=ab(9) + ab(6)
      ab(4)=ab(6)*ab(4)
      ab(2)=8*ab(2) + 5*ab(1)
      ab(2)=ab(2)*Pi**2
      ab(1)=1.0_ki/3.0_ki*ab(2) + ab(3) + ab(5) - 61.0_ki/9.0_ki*ab(1) + ab(4)
     &  + ab(7)

      tmp = 2*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=32*l3
      ab(2)= - l1 - 16*l2 + ab(1) + 44.0_ki/3.0_ki
      ab(3)=2*l2
      ab(2)=ab(2)*ab(3)
      ab(1)=ab(1) + 33
      ab(1)=ab(1)*l3
      ab(4)=l1**2
      ab(1)=ab(1) - ab(4)
      ab(5)=8*l3
      ab(6)=ab(5) + 11.0_ki/3.0_ki
      ab(7)= - 4*ab(6) + 31.0_ki/3.0_ki*l4 + 31*l2
      ab(7)=ab(7)*l4
      ab(8)=Pi**2
      ab(9)= - 61.0_ki/3.0_ki + 4*ab(8)
      ab(2)=ab(2) - ab(7) - ab(1) + 2.0_ki/3.0_ki*ab(9)
      ab(7)=4*l4
      ab(2)=ab(2)*ab(7)
      ab(7)= - t6 - t7 - t1 + 2*t2
      ab(9)=t4 + t3
      ab(10)= - 2*t5 + 80.0_ki/27.0_ki + 4*ab(9) - ab(7)
      ab(11)=d2 + d1
      ab(12)=27*zeta3
      ab(10)= - ab(12) + 4*ab(10) - 22.0_ki/3.0_ki*ab(11)
      ab(6)= - ab(6) + 8.0_ki/3.0_ki*l2
      ab(13)=4*l2
      ab(6)=ab(6)*ab(13)
      ab(14)= - 61.0_ki/3.0_ki + 5*ab(8)
      ab(1)= - ab(1) - ab(6) + 2.0_ki/3.0_ki*ab(14)
      ab(1)=ab(1)*ab(13)
      ab(6)=ab(5) + 11
      ab(13)=4.0_ki/3.0_ki*l3
      ab(6)=ab(6)*ab(13)
      ab(13)=4*ab(11)
      ab(6)= - ab(6) - 119.0_ki/9.0_ki + 10.0_ki/3.0_ki*ab(8) + ab(13)
      ab(14)=4*l3
      ab(6)=ab(6)*ab(14)
      ab(14)= - ab(8) + 2*ab(4)
      ab(14)=ab(14)*l1
      ab(15)=l5*z1
      ab(1)= - 2.0_ki/3.0_ki*ab(15) - ab(2) - 2*ab(10) - ab(1) + ab(6) + 4.D
     & 0/3.0_ki*ab(14) + 55.0_ki/9.0_ki*ab(8)
      ab(2)=z1**2*z3*z2
      ab(1)=ab(2)*CA*ab(1)
      ab(6)=ab(5) + l1
      ab(10)=8*l2
      ab(6)= - ab(10) + 3*ab(6)
      ab(6)=ab(6)*ab(3)
      ab(13)=5.0_ki/3.0_ki*ab(8) + ab(13)
      ab(15)=l3**2
      ab(4)=22*ab(15) + 3*ab(4)
      ab(16)= - 19.0_ki/3.0_ki*l4 - 19*l2 + 24*l3
      ab(16)=ab(16)*l4
      ab(6)=ab(6) + ab(16) - ab(4) + 2*ab(13)
      ab(13)=2*l4
      ab(6)=ab(6)*ab(13)
      ab(13)= - 2.0_ki/3.0_ki*l2 + 3*l3
      ab(10)=ab(13)*ab(10)
      ab(8)=1.0_ki/3.0_ki*ab(8)
      ab(13)=ab(8) + 2*ab(11)
      ab(4)=ab(10) - ab(4) + 4*ab(13)
      ab(3)=ab(4)*ab(3)
      ab(4)= - 4.0_ki/3.0_ki*ab(15) + ab(8) + ab(11)
      ab(4)=ab(4)*ab(5)
      ab(5)= - 6*t5 + 2*ab(9) - 3*ab(7)
      ab(3)= - ab(6) + ab(4) + ab(12) - 2*ab(14) - ab(3) + 4*ab(5)
      ab(2)=ab(2)*Ta2
      ab(2)=8*ab(2)
      ab(2)=ab(3)*ab(2)
      ab(1)=ab(1) + ab(2)
      ab(2)=lam**2
      ab(2)= - 1 + ab(2)

      tmp = Ta2*lam**3*ab(2)*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

