      subroutine hp_a1ncar0gg(xa,xb,res)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki) tmp,hp_cli2,hp_Li3
      real(ki) xa,xb,res(2,2,-2:2)
      complex(ki) ab(11)
      include 'hp_a1ncar0_functions.h'

!##### DOUBLE POLE #####

      
      ab(1)=xa - 3 + xb
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + 2 - xb

      tmp = CA*z6*z7*ab(1)
      res(1,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)= - xa + 2 - xb

      tmp = 2*CA*z6*z7*ab(1)
      res(1,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=xa + xb - 2
      ab(1)=ab(1)*l1
      ab(1)=1 + ab(1)

      tmp = 2*CA*z6*z7*ab(1)
      res(1,1, 0) = real(tmp,ki)

!##### ORDER EPS #####

      
      ab(1)= - xb + 2 - xa
      ab(1)=l1*ab(1)
      ab(1)= - 2 + ab(1)
      ab(1)=l1*ab(1)
      ab(1)= - 2*l2 + ab(1)

      tmp = CA*z6*z7*ab(1)
      res(1,1, 1) = real(tmp,ki)

!##### ORDER EPS^2 #####

      
      ab(1)=l2**2
      ab(2)=xb - 2 + xa
      ab(2)=l1*ab(2)
      ab(2)=1 + 1.0_ki/3.0_ki*ab(2)
      ab(2)=l1*ab(2)
      ab(2)=2*l2 + ab(2)
      ab(2)=l1*ab(2)
      ab(1)=ab(1) + ab(2)

      tmp = CA*z6*z7*ab(1)
      res(1,1, 2) = real(tmp,ki)

!##### DOUBLE POLE #####

      

      tmp =  0
      res(2,1,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      
      ab(1)=xa - 2
      ab(1)=xa*ab(1)
      ab(1)=ab(1) + 3
      ab(2)=CA*z1*z4*z6
      ab(1)=xa*ab(2)*ab(1)
      ab(1)=ab(1) - 2*ab(2)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) + ab(2)
      ab(2)=xb - 1

      tmp = 2*ab(2)*ab(1)
      res(2,1,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=xb - 1
      ab(2)= - l1 - l3 + l4
      ab(1)=ab(2)*ab(1)
      ab(1)= - 1 + ab(1)
      ab(2)=CA*z1*z4*z6
      ab(3)=ab(2)*xa
      ab(4)=2*ab(2)
      ab(3)=ab(3) - ab(4)
      ab(3)=ab(3)*xa
      ab(3)=ab(3) + 3*ab(2)
      ab(3)=ab(3)*xa
      ab(3)=ab(3) - ab(4)
      ab(3)=ab(3)*xa
      ab(2)=ab(3) + ab(2)

      tmp = 2*ab(2)*ab(1)
      res(2,1, 0) = real(tmp,ki)

!##### ORDER EPS #####

      
      ab(1)=xb - 1
      ab(2)=ab(1)*l1
      ab(3)=ab(2) + 1
      ab(4)=ab(1)*l3
      ab(5)=ab(3) + ab(4)
      ab(1)=ab(1)*l4
      ab(1)= - ab(1) + 2*ab(5)
      ab(1)=ab(1)*l4
      ab(3)=ab(4) + 2*ab(3)
      ab(3)=ab(3)*l3
      ab(2)=ab(2) + 2
      ab(2)=ab(2)*l1
      ab(1)= - ab(1) + ab(3) + ab(2) + 2*l2
      ab(2)=2*ab(1)
      ab(3)=xa*ab(1)
      ab(3)= - ab(2) + ab(3)
      ab(3)=xa*ab(3)
      ab(3)=3*ab(1) + ab(3)
      ab(3)=xa*ab(3)
      ab(2)= - ab(2) + ab(3)
      ab(2)=xa*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = CA*z1*z4*z6*ab(1)
      res(2,1, 1) = real(tmp,ki)

!##### ORDER EPS^2 #####

      
      ab(1)=xb - 1
      ab(2)=ab(1)*l1
      ab(3)=ab(2) + 1
      ab(4)=ab(1)*l3
      ab(5)=ab(4) + 2*ab(3)
      ab(5)=ab(5)*l3
      ab(6)=ab(2) + 2
      ab(6)=ab(6)*l1
      ab(7)=2*l2
      ab(6)=ab(6) + ab(7)
      ab(5)=ab(5) + ab(6)
      ab(8)=ab(3) + ab(4)
      ab(9)=1.0_ki/3.0_ki*l4
      ab(9)=ab(9)*ab(1)
      ab(9)=ab(9) - ab(8)
      ab(9)=ab(9)*l4
      ab(9)=ab(9) + ab(5)
      ab(9)=ab(9)*l4
      ab(10)=1.0_ki/3.0_ki*l3
      ab(10)=ab(10)*ab(1)
      ab(10)=ab(3) + ab(10)
      ab(10)=ab(10)*l3
      ab(10)=ab(10) + ab(6)
      ab(10)=ab(10)*l3
      ab(11)=1.0_ki/3.0_ki*l1
      ab(11)=ab(11)*ab(1)
      ab(11)=ab(11) + 1
      ab(11)=ab(11)*l1
      ab(7)=ab(11) + ab(7)
      ab(7)=ab(7)*l1
      ab(11)=l2**2
      ab(7)= - ab(9) + ab(10) + ab(7) + ab(11)
      ab(9)=2*ab(7)
      ab(10)= - xa*ab(7)
      ab(10)=ab(9) + ab(10)
      ab(10)=xa*ab(10)
      ab(3)= - 3*ab(3) - ab(4)
      ab(3)=l3*ab(3)
      ab(3)= - 3*ab(6) + ab(3)
      ab(3)=l3*ab(3)
      ab(1)=l4*ab(1)
      ab(1)= - 3*ab(8) + ab(1)
      ab(1)=l4*ab(1)
      ab(1)=3*ab(5) + ab(1)
      ab(1)=l4*ab(1)
      ab(2)= - 3 - ab(2)
      ab(2)=l1*ab(2)
      ab(2)= - 6*l2 + ab(2)
      ab(2)=l1*ab(2)
      ab(1)=ab(10) + ab(1) + ab(3) - 3*ab(11) + ab(2)
      ab(1)=xa*ab(1)
      ab(1)=ab(9) + ab(1)
      ab(1)=xa*ab(1)
      ab(1)=ab(1) - ab(7)

      tmp = CA*z1*z4*z6*ab(1)
      res(2,1, 2) = real(tmp,ki)

!##### DOUBLE POLE #####

      

      tmp =  0
      res(2,2,-2) = real(tmp,ki)

!##### SINGLE POLE #####

      

      tmp =  0
      res(2,2,-1) = real(tmp,ki)

!##### FINITE PART #####

      
      ab(1)=xa**4
      ab(2)=xb*xa**5
      ab(1)= - ab(1) + ab(2)
      ab(1)=xb*ab(1)
      ab(2)=xa**3
      ab(1)=ab(2) + ab(1)
      ab(1)=xb*ab(1)
      ab(2)=xa**2
      ab(1)=ab(2) + ab(1)
      ab(1)=xb*ab(1)
      ab(1)= - xa + ab(1)
      ab(1)=xb*ab(1)
      ab(1)=1 + ab(1)

      tmp = 4*CA*z1*z2*z3*z4*z5*ab(1)
      res(2,2, 0) = real(tmp,ki)

!##### ORDER EPS #####

      
      ab(1)=xa*xb**5
      ab(2)=xb**4
      ab(1)=ab(2) - ab(1)
      ab(1)=xa*ab(1)
      ab(2)=xb**3
      ab(1)=ab(1) - ab(2)
      ab(2)=CA*z1*z2*z3*z4*z5
      ab(1)=xa*ab(2)*ab(1)
      ab(3)=ab(2)*xb**2
      ab(1)=ab(1) - ab(3)
      ab(1)=ab(1)*xa
      ab(3)=ab(2)*xb
      ab(1)=ab(1) + ab(3)
      ab(1)=ab(1)*xa
      ab(1)=ab(1) - ab(2)
      ab(2)=l1 + l2 + l4 - 2*l6 + l5

      tmp = 4*ab(2)*ab(1)
      res(2,2, 1) = real(tmp,ki)

!##### ORDER EPS^2 #####

      
      ab(1)=l2 + l1
      ab(2)=ab(1) + l4
      ab(3)=ab(2) - l6 + l5
      ab(4)=4*l6
      ab(3)=ab(3)*ab(4)
      ab(2)=l5 + 2*ab(2)
      ab(2)=ab(2)*l5
      ab(1)=l4 + 2*ab(1)
      ab(1)=ab(1)*l4
      ab(4)=l2 + 2*l1
      ab(4)=ab(4)*l2
      ab(5)=l1**2
      ab(1)= - ab(3) + ab(2) + ab(1) + ab(4) + ab(5)
      ab(2)=ab(1)*xa**2
      ab(3)= - xa**4
      ab(4)=xb*xa**5
      ab(3)=ab(3) + ab(4)
      ab(3)=xb*ab(3)
      ab(4)=xa**3
      ab(3)=ab(4) + ab(3)
      ab(3)=xb*ab(1)*ab(3)
      ab(2)=ab(2) + ab(3)
      ab(2)=xb*ab(2)
      ab(3)= - xa*ab(1)
      ab(2)=ab(3) + ab(2)
      ab(2)=xb*ab(2)
      ab(1)=ab(2) + ab(1)

      tmp = 2*CA*z1*z2*z3*z4*z5*ab(1)
      res(2,2, 2) = real(tmp,ki)

      return
      end

