      subroutine hp_a1a1cassssragg(Ta2,lam,res)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) tmp,hp_cli2,hp_Li3
      real(ki) Ta2,lam,res(0:2,0:2,-4:0)
      complex(ki) ab(6)
      include 'hp_a1a1cassssr_functions.h'

!##### QUARTIC POLE #####

      
      ab(1)= - 1 + lam

      tmp = 2*z1*Ta2**2*ab(1)
      res(1,1,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      
      ab(1)= - lam + 1
      ab(1)=ab(1)*l1
      ab(1)= - 1 + ab(1)

      tmp = 8*z1*Ta2**2*ab(1)
      res(1,1,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=lam - 1
      ab(2)= - ab(1)*Pi**2
      ab(1)=l1*ab(1)
      ab(1)=2 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)=8*ab(1) + 16*l2 + ab(2)

      tmp = 2*z1*Ta2**2*ab(1)
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

      tmp = 8*z1*Ta2**2*ab(1)
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

      tmp = 2*z1*Ta2**2*ab(1)
      res(1,1,0) = real(tmp,ki)

!##### QUARTIC POLE #####

      

      tmp =  0
      res(0,0,-4) = real(tmp,ki)

!##### TRIPLE POLE #####

      

      tmp = 16*z1*z2*Ta2**2*lam**3
      res(0,0,-3) = real(tmp,ki)

!##### DOUBLE POLE #####

      
      ab(1)=l2 + l4
      ab(1)=5*l3 - 4*ab(1)

      tmp = 16*z1*z2*Ta2**2*lam**3*ab(1)
      res(0,0,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=l4**2
      ab(2)= - 20*l4 + 11*l3
      ab(2)=l3*ab(2)
      ab(3)=2*l2 + 4*l4 - 5*l3
      ab(3)=l2*ab(3)
      ab(4)= - Pi**2
      ab(1)=4*ab(3) + ab(2) + ab(4) + 8*ab(1)

      tmp = 16*z1*z2*Ta2**2*lam**3*ab(1)
      res(0,0,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=l2**2
      ab(2)=2*l2
      ab(3)=ab(2) + l4
      ab(3)=l4*ab(3)
      ab(3)=ab(1) + ab(3)
      ab(4)=Pi**2
      ab(5)=d1 + d2
      ab(3)=3*ab(5) - ab(4) + 10*ab(3)
      ab(5)= - l2 - l4
      ab(5)=22*ab(5) + 23.0_ki/3.0_ki*l3
      ab(5)=l3*ab(5)
      ab(3)=2*ab(3) + ab(5)
      ab(3)=l3*ab(3)
      ab(5)=ab(4) - 8.0_ki/3.0_ki*ab(1)
      ab(2)=ab(5)*ab(2)
      ab(5)= - l2 - 1.0_ki/3.0_ki*l4
      ab(5)=l4*ab(5)
      ab(1)=ab(1) - ab(5)
      ab(1)=ab(4) - 8*ab(1)
      ab(1)=l4*ab(1)
      ab(4)= - t1 - t2
      ab(1)=ab(3) + 2*ab(1) + ab(2) + 12*ab(4) - 7*zeta3

      tmp = 32*z1*z2*Ta2**2*lam**3*ab(1)
      res(0,0,0) = real(tmp,ki)

      return
      end

