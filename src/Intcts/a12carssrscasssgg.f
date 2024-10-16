      subroutine a12carssrscasssgg(Ta2,lam,res)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki) tmp,cli2,li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(7)
      include 'a12carssrscasss_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)= - 1 + lam

      tmp = z1*Ta2**2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - lam + 1
      ab(1)=ab(1)*l1
      ab(1)= - 1 + ab(1)

      tmp = 4*z1*Ta2**2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=lam - 1
      ab(2)= - ab(1)*Pi**2
      ab(1)=l1*ab(1)
      ab(1)=2 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)=8*ab(1) + 16*l2 + ab(2)

      tmp = z1*Ta2**2*ab(1)
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

      tmp = 4*z1*Ta2**2*ab(1)
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

      tmp = z1*Ta2**2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 8*z1*z2*Ta2**2*lam**3
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l2 + l4
      ab(1)=5*l3 - 4*ab(1)

      tmp = 8*z1*z2*Ta2**2*lam**3*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=d1 + d5
      ab(2)=l4**2
      ab(3)= - 20*l4 + 11*l3
      ab(3)=l3*ab(3)
      ab(4)= - 5*l3 + 2*l2
      ab(4)=l2*ab(4)
      ab(5)= - l1 + l4 + l2
      ab(5)=l1*ab(5)
      ab(6)=Pi**2
      ab(1)=16*ab(5) + 4*ab(4) + ab(3) + 5.0_ki/3.0_ki*ab(6) + 8*ab(2) - 16
     & *ab(1)

      tmp = 8*z1*z2*Ta2**2*lam**3*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=d5 + d1
      ab(2)=l1**2
      ab(3)=Pi**2
      ab(1)=16*ab(2) - 5.0_ki/3.0_ki*ab(3) + 16*ab(1)
      ab(2)=2*l2
      ab(4)=8*l1 + 7*l2
      ab(4)=ab(4)*ab(2)
      ab(5)=2*l4
      ab(6)=4*l1 + l2
      ab(6)=2*ab(6) + 5*l4
      ab(6)=ab(6)*ab(5)
      ab(7)=d4 + d3
      ab(4)=ab(6) + ab(4) + 3*ab(7) - ab(1)
      ab(6)= - l2 - l4
      ab(6)=22*ab(6) + 23.0_ki/3.0_ki*l3
      ab(6)=l3*ab(6)
      ab(4)=2*ab(4) + ab(6)
      ab(4)=l3*ab(4)
      ab(6)=d2 + d4
      ab(7)=l2**2
      ab(3)= - 8.0_ki/3.0_ki*ab(7) + 4*ab(6) + 1.0_ki/3.0_ki*ab(3)
      ab(2)=ab(3)*ab(2)
      ab(3)=2*l1
      ab(6)= - ab(3) - l2
      ab(6)=l2*ab(6)
      ab(3)= - 1.0_ki/3.0_ki*l4 - ab(3) + l2
      ab(3)=l4*ab(3)
      ab(3)=ab(6) + ab(3)
      ab(1)=8*ab(3) + ab(1)
      ab(1)=ab(1)*ab(5)
      ab(3)= - t1 - t2
      ab(3)=10*ab(3) + t3
      ab(1)=ab(4) + ab(1) + ab(2) + 2*ab(3) - 7*zeta3

      tmp = 16*z1*z2*Ta2**2*lam**3*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

