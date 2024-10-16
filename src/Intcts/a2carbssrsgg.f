      subroutine a2carbssrsgg(Ta2,Tb2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,Tb2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(10)
      include 'a2carbssrs_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)= - 1 + lam

      tmp = z1*Ta2*Tb2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - lam + 1
      ab(1)=ab(1)*l1
      ab(1)= - 1 + ab(1)

      tmp = 4*z1*Ta2*Tb2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=lam - 1
      ab(2)= - ab(1)*Pi**2
      ab(1)=l1*ab(1)
      ab(1)=2 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)=8*ab(1) + 16*l2 + ab(2)

      tmp = z1*Ta2*Tb2*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=lam - 1
      ab(2)=5*zeta3
      ab(2)= - ab(2)*ab(1)
      ab(3)=Pi**2
      ab(4)=ab(3)*ab(1)
      ab(1)=l1*ab(1)
      ab(1)= - 1 - 1.0_ki/3.0_ki*ab(1)
      ab(1)=l1*ab(1)
      ab(1)=8*ab(1) - 16*l2 + ab(4)
      ab(1)=l1*ab(1)
      ab(4)=l2**2
      ab(1)=ab(1) + ab(3) - 8*ab(4) + ab(2)

      tmp = 4*z1*Ta2*Tb2*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=lam - 1
      ab(2)=5*zeta3
      ab(3)=ab(2)*ab(1)
      ab(4)=Pi**2
      ab(5)=l2**2
      ab(3)= - ab(4) + 8*ab(5) + ab(3)
      ab(5)=ab(1)*ab(4)
      ab(6)=16*l2
      ab(1)=l1*ab(1)
      ab(1)=4 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)=4.0_ki/3.0_ki*ab(1) + ab(6) - ab(5)
      ab(1)=l1*ab(1)
      ab(1)=2*ab(3) + ab(1)
      ab(1)=l1*ab(1)
      ab(3)=l2**3
      ab(2)=ab(2) + 8.0_ki/3.0_ki*ab(3)
      ab(3)= - ab(6) - 1.0_ki/5.0_ki*ab(5)
      ab(3)=ab(3)*ab(4)
      ab(1)=8*ab(1) + 16*ab(2) + ab(3)

      tmp = z1*Ta2*Tb2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 8*z1*z2*Ta2*Tb2*lam**3
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l2 + l4
      ab(1)=3*l3 - 2*ab(1)

      tmp = 16*z1*z2*Ta2*Tb2*lam**3*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l4**2
      ab(2)= - 5*l4 + 3*l3
      ab(2)=l3*ab(2)
      ab(3)=2*l2 + 4*l4 - 5*l3
      ab(3)=l2*ab(3)
      ab(4)=Pi**2
      ab(1)=ab(3) + ab(2) - 1.0_ki/3.0_ki*ab(4) + d1 + d2 + 2*ab(1)

      tmp = 32*z1*z2*Ta2*Tb2*lam**3*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l1**2
      ab(2)=l3**2
      ab(3)= - 3*ab(1) + 20*ab(2)
      ab(4)=Pi**2
      ab(5)=1.0_ki/3.0_ki*ab(4)
      ab(6)=d2 + d1
      ab(7)= - 4*ab(6) + ab(5)
      ab(8)=2*l2
      ab(9)=16*l3
      ab(10)= - 8*l2 - 3*l1 + ab(9)
      ab(10)=ab(10)*ab(8)
      ab(9)= - 13.0_ki/3.0_ki*l4 + ab(9) - 13*l2
      ab(9)=l4*ab(9)
      ab(7)=ab(9) + ab(10) + 2*ab(7) - ab(3)
      ab(7)=l4*ab(7)
      ab(1)=ab(4) - 2*ab(1)
      ab(1)=l1*ab(1)
      ab(1)=ab(1) + ab(7)
      ab(2)=4*ab(2) + 10*ab(6) - ab(4)
      ab(2)=l3*ab(2)
      ab(4)=t6 + t7 + t1 - 2*t2
      ab(7)=t4 + t3
      ab(2)=ab(2) - 6*t5 - 16*ab(7) + 3*ab(4)
      ab(4)=ab(5) - ab(6)
      ab(5)=l3 - 1.0_ki/3.0_ki*l2
      ab(5)=l2*ab(5)
      ab(3)=16*ab(5) + 8*ab(4) - ab(3)
      ab(3)=ab(3)*ab(8)
      ab(1)=ab(3) + 5*zeta3 + 4*ab(2) + 2*ab(1)

      tmp = 8*z1*z2*Ta2*Tb2*lam**3*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

