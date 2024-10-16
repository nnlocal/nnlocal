      subroutine AP0(xa)
      implicit none
      include 'types.h'
      include 'constants.h'
      include 'APinc.h'
      real(ki) xa
      real(ki) xlf
      xlf=real(nf,kind=ki)
      xlf=zip
      
      Pgg0(1,1) = (11*CA - 2*xlf)/6 + (2*CA)/(-1 + xa)

      Pgq0(1,1) = 0.0_ki

      Pqg0(1,1) = 0.0_ki

      Pqq0(1,1) = (3*CF)/2 + (2*CF)/(-1 + xa)

      Pgg0(2,1) = -4*CA - (2*CA)/(-1 + xa) + (2*CA)/xa + 2*CA*xa - 2*CA*xa**2

      Pgq0(2,1) = -2*CF + (2*CF)/xa + CF*xa

      Pqg0(2,1) = half - xa + xa**2

      Pqq0(2,1) = -CF - (2*CF)/(-1 + xa) - CF*xa

      Pgg0(1,2) = 0.0_ki

      Pgq0(1,2) = 0.0_ki

      Pqg0(1,2) = 0.0_ki

      Pqq0(1,2) = 0.0_ki

      Pgg0(2,2) = 0.0_ki

      Pgq0(2,2) = 0.0_ki

      Pqg0(2,2) = 0.0_ki

      Pqq0(2,2) = 0.0_ki

      return
      end
      

      subroutine AP1(xa)
      implicit none
      include 'types.h'
      include 'constants.h'
      include 'APinc.h'
      real(ki) xa
      real(ki) l1,l2,l3,l1sq,l2sq,d1
      complex(ki) cli2
      real(ki) xlf
      xlf=real(nf,kind=ki)
      xlf=zip

      l1=log(xa)
      l2=log(1-xa)
      l3=log(1+xa)
      l1sq=l1*l1
      l2sq=l2*l2
      d1=real(cli2(-xa*cone),kind=ki)

      Pgg1(1,1) = (67*CA**2 - 10*CA*xlf - 3*CA**2*pisq)/(9*(-1 + xa)) + 
     1    (16*CA**2 - 4*CA*xlf - 3*CF*xlf + 18*CA**2*zeta3)/6

      Pgq1(1,1) = 0.0_ki

      Pqg1(1,1) = 0.0_ki

      Pqq1(1,1) = (67*CA*CF - 10*CF*xlf - 3*CA*CF*pisq)/(9*(-1 + xa)) - 
     1    (CF*(-51*CA - 27*CF + 6*xlf - 44*CA*pisq + 36*CF*pisq +
     2    8*xlf*pisq + 216*CA*zeta3 - 432*CF*zeta3))/72

      Pgg1(2,1) = (-67*CA**2 + 10*CA*xlf + 3*CA**2*pisq + 36*CA**2*l2*l1 - 
     1    9*CA**2*l1sq)/(9*(-1 + xa)) + 
     2    (-(CA**2*pisq) + 3*CA**2*l1sq - 12*CA**2*l1*l3 - 
     3    12*CA**2*d1)/(3*(1 + xa)) + 
     4    (-23*CA*xlf + 6*CF*xlf - 36*CA**2*l2*l1 + 
     5    36*CA**2*l1*l3 + 36*CA**2*d1)/(9*xa) + 
     6    (xa**2*(23*CA*xlf + 30*CF*xlf + 6*CA**2*pisq - 132*CA**2*l1 + 
     7    36*CA**2*l2*l1 - 18*CA**2*l1sq + 
     8    36*CA**2*l1*l3 + 36*CA**2*d1))/9 + 
     9    (xa*(-109*CA**2 - 38*CA*xlf + 72*CF*xlf + 66*CA**2*l1 - 
     1    12*CA*xlf*l1 - 90*CF*xlf*l1 - 72*CA**2*l2*l1 + 
     2    72*CA**2*l1sq - 18*CF*xlf*l1sq + 72*CA**2*l1*l3 + 
     3    72*CA**2*d1))/18 + (-25*CA**2 + 58*CA*xlf - 144*CF*xlf + 
     4    24*CA**2*pisq - 150*CA**2*l1 - 12*CA*xlf*l1 - 54*CF*xlf*l1 + 
     5    144*CA**2*l2*l1 - 18*CF*xlf*l1sq + 
     6    144*CA**2*l1*l3 + 144*CA**2*d1)/18

      Pgq1(2,1) = (-4*CA*CF*xa**2*(-11 + 6*l1))/9 + 
     1    (CF*xa*(74*CA - 63*CF - 32*xlf + 102*CA*l2 - 90*CF*l2 - 
     2    12*xlf*l2 + 18*CA*l2sq - 18*CF*l2sq - 
     3    90*CA*l1 + 63*CF*l1 - 36*CA*l2*l1 + 
     4    18*CA*l1sq + 9*CF*l1sq + 36*CA*l1*l3 + 
     5    36*CA*d1))/18 + 
     6    (CF*(38*CA - 45*CF + 40*xlf + 12*CA*pisq - 132*CA*l2 + 
     7    108*CF*l2 + 24*xlf*l2 - 36*CA*l2sq + 
     8    36*CF*l2sq - 216*CA*l1 + 36*CF*l1 + 
     9    72*CA*l2*l1 + 36*CA*l1sq - 18*CF*l1sq + 
     1    72*CA*l1*l3 + 72*CA*d1))/18 + 
     2    (9*CA*CF - 20*CF*xlf + 66*CA*CF*l2 - 54*CF**2*l2 - 
     3    12*CF*xlf*l2 + 18*CA*CF*l2sq - 18*CF**2*l2sq - 
     4    36*CA*CF*l2*l1 + 36*CA*CF*l1*l3 + 
     5    36*CA*CF*d1)/(9*xa)

      Pqg1(2,1) = (10*CA)/(9*xa) + (-12*CA + 42*CF - 2*CF*pisq -
     1    6*CA*l2sq + 6*CF*l2sq + 6*CA*l1 + 9*CF*l1 - 
     2    12*CF*l2*l1 - 6*CA*l1sq + 3*CF*l1sq - 
     3    12*CA*l1*l3 - 12*CA*d1)/12 - 
     4    (xa**2*(109*CA - 45*CF + 3*CF*pisq - 18*CA*l2 + 
     5    18*CF*l2 + 9*CA*l2sq - 9*CF*l2sq - 
     6    66*CA*l1 - 18*CF*l1 + 18*CF*l2*l1 - 
     7    9*CF*l1sq + 18*CA*l1*l3 + 18*CA*d1))/
     8    9 - (xa*(-150*CA + 87*CF + 4*CA*pisq - 4*CF*pisq + 24*CA*l2 - 
     9    24*CF*l2 - 12*CA*l2sq + 12*CF*l2sq - 
     1    48*CA*l1 + 12*CF*l1 - 24*CF*l2*l1 + 
     2    12*CA*l1sq + 6*CF*l1sq + 24*CA*l1*l3 + 
     3    24*CA*d1))/12

      Pqq1(2,1) = (10*CF)/(9*xa) + (4*CF*xa**2*(-7 + 3*l1))/9 - 
     1    (CF*(18 - 53*CA + 90*CF + 2*xlf - 3*CA*pisq - 9*l1 + 15*CA*l1 - 
     2    6*xlf*l1 - 36*CF*l2*l1 + 9*l1sq + 
     3    9*CA*l1sq + 9*CF*l1sq))/18 - 
     4    (CF*xa*(-54 + 187*CA - 90*CF - 22*xlf - 3*CA*pisq - 45*l1 + 
     5    15*CA*l1 + 36*CF*l1 - 6*xlf*l1 - 
     6    36*CF*l2*l1 + 9*l1sq + 9*CA*l1sq + 
     7    9*CF*l1sq))/18 + (-67*CA*CF + 10*CF*xlf + 3*CA*CF*pisq - 
     8    33*CA*CF*l1 + 27*CF**2*l1 + 6*CF*xlf*l1 + 
     9    36*CF**2*l2*l1 - 9*CA*CF*l1sq)/(9*(-1 + xa))

      Pgg1(1,2) = 0.0_ki

      Pgq1(1,2) = 0.0_ki

      Pqg1(1,2) = 0.0_ki

      Pqq1(1,2) = 0.0_ki

      Pgg1(2,2) = 0.0_ki

      Pgq1(2,2) = 0.0_ki

      Pqg1(2,2) = 0.0_ki

      Pqq1(2,2) = 0.0_ki

      return
      end
      
      
      subroutine APxAP(xa)
      implicit none
      include 'types.h'
      include 'constants.h'
      include 'APinc.h'
      real(ki) xa
      real(ki) l1,l2
      real(ki) xlf
      xlf=real(nf,kind=ki)
      xlf=zip

      l1=log(xa)
      l2=log(1-xa)
      
      P0ggxP0gg(1,1) = (121*CA**2 - 44*CA*xlf + 4*xlf**2 - 24*xn**2*pisq)/36 + 
     1    (2*(11*CA*xn - 2*xn*xlf + 12*xn**2*l2))/(3*(-1 + xa))

      P0ggxP0gq(1,1) = 0.0_ki

      P0ggxP0qg(1,1) = 0.0_ki

      P0ggxP0qq(1,1) = -1/24.0_ki*((-1 + xn)*(1 + xn)*(-33*CA + 6*xlf + 8*xn*pisq))/xn + 
     1    (-11*CA - 9*xn + 11*CA*xn**2 + 9*xn**3 + 2*xlf - 2*xn**2*xlf - 
     2    24*xn*l2 + 24*xn**3*l2)/(6*xn*(-1 + xa))

      P0gqxP0gg(1,1) = 0.0_ki

      P0gqxP0gq(1,1) = 0.0_ki

      P0gqxP0qg(1,1) = 0.0_ki

      P0gqxP0qq(1,1) = 0.0_ki

      P0qgxP0gg(1,1) = 0.0_ki

      P0qgxP0gq(1,1) = 0.0_ki

      P0qgxP0qg(1,1) = 0.0_ki

      P0qgxP0qq(1,1) = 0.0_ki

      P0qqxP0gg(1,1) = -1/24.0_ki*((-1 + xn)*(1 + xn)*(-33*CA + 6*xlf + 8*xn*pisq))/xn + 
     1    (-11*CA - 9*xn + 11*CA*xn**2 + 9*xn**3 + 2*xlf - 2*xn**2*xlf - 
     2    24*xn*l2 + 24*xn**3*l2)/(6*xn*(-1 + xa))

      P0qqxP0gq(1,1) = 0.0_ki

      P0qqxP0qg(1,1) = 0.0_ki

      P0qqxP0qq(1,1) = -1/48.0_ki*((-1 + xn)**2*(1 + xn)**2*(-27 + 8*pisq))/xn**2 + 
     1    (3 - 6*xn**2 + 3*xn**4 + 4*l2 - 8*xn**2*l2 + 
     2    4*xn**4*l2)/(2*xn**2*(-1 + xa))

      P0ggxP0gg(2,1) = (-4*xn*(11*CA - 9*xn - 2*xlf + 12*xn*l2))/3 + 
     1    (2*xn*xa*(11*CA - 18*xn - 2*xlf + 12*xn*l2 - 18*xn*l1))/3 - 
     2    (2*xn*xa**2*(11*CA - 22*xn - 2*xlf + 12*xn*l2 - 6*xn*l1))/3 - 
     3    (2*(11*CA*xn - 2*xn*xlf + 12*xn**2*l2 - 6*xn**2*l1))/
     4    (3*(-1 + xa)) + (2*(11*CA*xn - 22*xn**2 - 2*xn*xlf + 12*xn**2*l2 - 
     5    6*xn**2*l1))/(3*xa)

      P0ggxP0gq(2,1) = (2*(-1 + xn**2)*xa**2)/3 + 
     1    ((-1 + xn**2)*xa*(11*CA + 6*xn - 2*xlf + 12*xn*l2 - 
     2    24*xn*l1))/(12*xn) - 
     3    ((-1 + xn**2)*(11*CA - 24*xn - 2*xlf + 12*xn*l2 + 12*xn*l1))/
     4    (6*xn) + (-11*CA + 31*xn + 11*CA*xn**2 - 31*xn**3 + 2*xlf - 2*xn**2*xlf - 
     5    12*xn*l2 + 12*xn**3*l2 + 12*xn*l1 - 
     6    12*xn**3*l1)/(6*xn*xa)

      P0ggxP0qg(2,1) = (2*xn)/(3*xa) + 
     1    (xa**2*(11*CA - 31*xn - 2*xlf + 12*xn*l2))/6 - 
     2    (xa*(11*CA - 24*xn - 2*xlf + 12*xn*l2 - 24*xn*l1))/6 + 
     3    (11*CA + 6*xn - 2*xlf + 12*xn*l2 + 12*xn*l1)/12

      P0ggxP0qq(2,1) = -1/12.0_ki*((-1 + xn**2)*xa*(11*CA + 12*xn - 2*xlf - 
     1    12*xn*l2))/xn + (2*(-l2 + xn**2*l2))/xa - 
     2    2*(-1 + xn**2)*xa**2*(l2 - l1) - 
     3    ((-1 + xn**2)*(11*CA + 6*xn - 2*xlf + 60*xn*l2 - 36*xn*l1))/
     4    (12*xn) + (11*CA + 9*xn - 11*CA*xn**2 - 9*xn**3 - 2*xlf + 2*xn**2*xlf + 
     5    24*xn*l2 - 24*xn**3*l2 - 12*xn*l1 + 
     6    12*xn**3*l1)/(6*xn*(-1 + xa))

      P0gqxP0gg(2,1) = (2*(-1 + xn**2)*xa**2)/3 + 
     1    ((-1 + xn**2)*xa*(11*CA + 6*xn - 2*xlf + 12*xn*l2 - 
     2    24*xn*l1))/(12*xn) - 
     3    ((-1 + xn**2)*(11*CA - 24*xn - 2*xlf + 12*xn*l2 + 12*xn*l1))/
     4    (6*xn) + (-11*CA + 31*xn + 11*CA*xn**2 - 31*xn**3 + 2*xlf - 2*xn**2*xlf - 
     5    12*xn*l2 + 12*xn**3*l2 + 12*xn*l1 - 
     6    12*xn**3*l1)/(6*xn*xa)

      P0gqxP0gq(2,1) = -1/4.0_ki*((-1 + xn**2)**2*xa*(-2 + l1))/xn**2 - 
     1    ((-1 + xn**2)**2*(-1 + l1))/xn**2 + 
     2    (-3 + 6*xn**2 - 3*xn**4 - 2*l1 + 4*xn**2*l1 - 2*xn**4*l1)/
     3    (2*xn**2*xa)

      P0gqxP0qg(2,1) = ((-1 + xn)*(1 + xn))/(3*xn*xa) - ((-1 + xn)*(1 + xn)*xa**2)/
     1    (3*xn) + ((-1 + xn)*(1 + xn)*xa*(-1 + 2*l1))/(4*xn) + 
     2    ((-1 + xn)*(1 + xn)*(1 + 2*l1))/(4*xn)
     
      P0gqxP0qq(2,1) = (l2 - 2*xn**2*l2 + xn**4*l2)/
     1    (xn**2*xa) + ((-1 + xn**2)**2*xa*(-1 + 4*l2 - 2*l1))/
     2    (8*xn**2) - ((-1 + xn**2)**2*(-1 + 2*l2 - l1))/(2*xn**2)
     
      P0qgxP0gg(2,1) = (2*xn)/(3*xa) + 
     1    (xa**2*(11*CA - 31*xn - 2*xlf + 12*xn*l2))/6 - 
     2    (xa*(11*CA - 24*xn - 2*xlf + 12*xn*l2 - 24*xn*l1))/6 + 
     3    (11*CA + 6*xn - 2*xlf + 12*xn*l2 + 12*xn*l1)/12
     
      P0qgxP0gq(2,1) = ((-1 + xn)*(1 + xn))/(3*xn*xa) - ((-1 + xn)*(1 + xn)*xa**2)/
     1    (3*xn) + ((-1 + xn)*(1 + xn)*xa*(-1 + 2*l1))/(4*xn) + 
     2    ((-1 + xn)*(1 + xn)*(1 + 2*l1))/(4*xn)
     
      P0qgxP0qg(2,1) = (-2 - l1)/4 - xa*(1 + l1) - (xa**2*(-3 + 2*l1))/2
     
      P0qgxP0qq(2,1) = ((-1 + xn**2)*(-1 + 4*l2 - 2*l1))/(8*xn) + 
     1    ((-1 + xn**2)*xa**2*(l2 - l1))/xn - 
     2    ((-1 + xn**2)*xa*(-1 + 2*l2 - l1))/(2*xn)
     
      P0qqxP0gg(2,1) = -1/12.0_ki*((-1 + xn**2)*xa*(11*CA + 12*xn - 2*xlf - 
     1    12*xn*l2))/xn + (2*(-l2 + xn**2*l2))/xa - 
     2    2*(-1 + xn**2)*xa**2*(l2 - l1) - 
     3    ((-1 + xn**2)*(11*CA + 6*xn - 2*xlf + 60*xn*l2 - 36*xn*l1))/
     4    (12*xn) + (11*CA + 9*xn - 11*CA*xn**2 - 9*xn**3 - 2*xlf + 2*xn**2*xlf + 
     5    24*xn*l2 - 24*xn**3*l2 - 12*xn*l1 + 
     6    12*xn**3*l1)/(6*xn*(-1 + xa))
       
      P0qqxP0gq(2,1) = (l2 - 2*xn**2*l2 + xn**4*l2)/
     1    (xn**2*xa) + ((-1 + xn**2)**2*xa*(-1 + 4*l2 - 2*l1))/
     2    (8*xn**2) - ((-1 + xn**2)**2*(-1 + 2*l2 - l1))/(2*xn**2)
       
      P0qqxP0qg(2,1) = ((-1 + xn**2)*(-1 + 4*l2 - 2*l1))/(8*xn) + 
     1    ((-1 + xn**2)*xa**2*(l2 - l1))/xn - 
     2    ((-1 + xn**2)*xa*(-1 + 2*l2 - l1))/(2*xn)
     
      P0qqxP0qq(2,1) = -1/4.0_ki*((-1 + xn**2)**2*xa*(1 + 4*l2 - 3*l1))/xn**2 - 
     1    ((-1 + xn**2)**2*(5 + 4*l2 - 3*l1))/(4*xn**2) + 
     2    (-3 + 6*xn**2 - 3*xn**4 - 4*l2 + 8*xn**2*l2 - 
     3    4*xn**4*l2 + 2*l1 - 4*xn**2*l1 + 2*xn**4*l1)/
     4    (2*xn**2*(-1 + xa))
      
      P0ggxP0gg(1,2) = 0.0_ki
      
      P0ggxP0gq(1,2) = 0.0_ki
      
      P0ggxP0qg(1,2) = 0.0_ki
      
      P0ggxP0qq(1,2) = 0.0_ki
      
      P0gqxP0gg(1,2) = 0.0_ki
      
      P0gqxP0gq(1,2) = 0.0_ki
      
      P0gqxP0qg(1,2) = 0.0_ki
      
      P0gqxP0qq(1,2) = 0.0_ki
      
      P0qgxP0gg(1,2) = 0.0_ki
      
      P0qgxP0gq(1,2) = 0.0_ki
      
      P0qgxP0qg(1,2) = 0.0_ki
      
      P0qgxP0qq(1,2) = 0.0_ki
      
      P0qqxP0gg(1,2) = 0.0_ki
      
      P0qqxP0gq(1,2) = 0.0_ki
      
      P0qqxP0qg(1,2) = 0.0_ki
      
      P0qqxP0qq(1,2) = 0.0_ki
      
      P0ggxP0gg(2,2) = 0.0_ki

      P0ggxP0gq(2,2) = 0.0_ki
      
      P0ggxP0qg(2,2) = 0.0_ki
      
      P0ggxP0qq(2,2) = 0.0_ki
      
      P0gqxP0gg(2,2) = 0.0_ki
      
      P0gqxP0gq(2,2) = 0.0_ki
      
      P0gqxP0qg(2,2) = 0.0_ki
      
      P0gqxP0qq(2,2) = 0.0_ki
      
      P0qgxP0gg(2,2) = 0.0_ki
      
      P0qgxP0gq(2,2) = 0.0_ki

      P0qgxP0qg(2,2) = 0.0_ki
      
      P0qgxP0qq(2,2) = 0.0_ki
      
      P0qqxP0gg(2,2) = 0.0_ki
      
      P0qqxP0gq(2,2) = 0.0_ki
      
      P0qqxP0qg(2,2) = 0.0_ki
      
      P0qqxP0qq(2,2) = 0.0_ki


      return
      end
      

      subroutine APAP(xa,xb)
      implicit none
      include 'types.h'
      include 'constants.h'
      include 'APinc.h'
      real(ki) xa,xb
      real(ki) xlf
      xlf=real(nf,kind=ki)
      xlf=zip
      
      P0ggP0gg(1,1) = (121*CA**2 - 44*CA*xlf + 4*xlf**2)/36 + 
     1    (132*CA**2 - 24*CA*xlf)/(36*(-1 + xa)) + (132*CA**2 - 24*CA*xlf)/
     2    (36*(-1 + xb)) + (4*CA**2)/((-1 + xa)*(-1 + xb))
     
      P0ggP0gq(1,1) = 0.0_ki
     
      P0ggP0qg(1,1) = 0.0_ki
     
      P0ggP0qq(1,1) = (33*CA*CF - 6*CF*xlf)/12 + (3*CA*CF)/(-1 + xa) + 
     1    (44*CA*CF - 8*CF*xlf)/(12*(-1 + xb)) + (4*CA*CF)/((-1 + xa)*(-1 + xb))
     
      P0gqP0gg(1,1) = 0.0_ki
     
      P0gqP0gq(1,1) = 0.0_ki
     
      P0gqP0qg(1,1) = 0.0_ki
     
      P0gqP0qq(1,1) = 0.0_ki
     
      P0qgP0gg(1,1) = 0.0_ki
     
      P0qgP0gq(1,1) = 0.0_ki
     
      P0qgP0qg(1,1) = 0.0_ki
     
      P0qgP0qq(1,1) = 0.0_ki
     
      P0qqP0gg(1,1) = (33*CA*CF - 6*CF*xlf)/12 + (44*CA*CF - 8*CF*xlf)/(12*(-1 + xa)) + 
     1    (3*CA*CF)/(-1 + xb) + (4*CA*CF)/((-1 + xa)*(-1 + xb))
     
      P0qqP0gq(1,1) = 0.0_ki
     
      P0qqP0qg(1,1) = 0.0_ki
     
      P0qqP0qq(1,1) = (9*CF**2)/4 + (3*CF**2)/(-1 + xa) + (3*CF**2)/(-1 + xb) + 
     1    (4*CF**2)/((-1 + xa)*(-1 + xb))
     
      P0ggP0gg(2,1) = (-11*CA**2 + 2*CA*xlf)/(3*(-1 + xa)) + 
     1    (11*CA**2 - 2*CA*xlf)/(3*xa) + (-22*CA**2 + 4*CA*xlf + 11*CA**2*xa - 
     2    2*CA*xlf*xa - 11*CA**2*xa**2 + 2*CA*xlf*xa**2)/3 - 
     3    (4*CA**2)/((-1 + xa)*(-1 + xb)) + (4*CA**2)/(xa*(-1 + xb)) + 
     4    (-24*CA**2 + 12*CA**2*xa - 12*CA**2*xa**2)/(3*(-1 + xb))
     
      P0ggP0gq(2,1) = 0.0_ki
     
      P0ggP0qg(2,1) = 0.0_ki
     
      P0ggP0qq(2,1) = -6*CA*CF - (3*CA*CF)/(-1 + xa) + (3*CA*CF)/xa + 3*CA*CF*xa - 
     1    3*CA*CF*xa**2 - (4*CA*CF)/((-1 + xa)*(-1 + xb)) + 
     2    (4*CA*CF)/(xa*(-1 + xb)) + (-8*CA*CF + 4*CA*CF*xa - 4*CA*CF*xa**2)/
     3    (-1 + xb)
     
      P0gqP0gg(2,1) = (22*CA*CF - 4*CF*xlf)/(6*xa) + 
     1    (-22*CA*CF + 4*CF*xlf + 11*CA*CF*xa - 2*CF*xlf*xa)/6 + 
     2    (4*CA*CF)/(xa*(-1 + xb)) + (-24*CA*CF + 12*CA*CF*xa)/(6*(-1 + xb))
     
      P0gqP0gq(2,1) = 0.0_ki
     
      P0gqP0qg(2,1) = 0.0_ki
     
      P0gqP0qq(2,1) = (3*CF**2)/xa + (-6*CF**2 + 3*CF**2*xa)/2 + 
     1    (4*CF**2)/(xa*(-1 + xb)) + (-8*CF**2 + 4*CF**2*xa)/(2*(-1 + xb))
     
      P0qgP0gg(2,1) = (11*CA - 2*xlf - 22*CA*xa + 4*xlf*xa + 22*CA*xa**2 - 4*xlf*xa**2)/
     1    12 + (12*CA - 24*CA*xa + 24*CA*xa**2)/(12*(-1 + xb))
     
      P0qgP0gq(2,1) = 0.0_ki
     
      P0qgP0qg(2,1) = 0.0_ki
     
      P0qgP0qq(2,1) = (3*CF - 6*CF*xa + 6*CF*xa**2)/4 + (4*CF - 8*CF*xa + 8*CF*xa**2)/
     1    (4*(-1 + xb))
     
      P0qqP0gg(2,1) = (-22*CA*CF + 4*CF*xlf)/(6*(-1 + xa)) + 
     1    (-11*CA*CF + 2*CF*xlf - 11*CA*CF*xa + 2*CF*xlf*xa)/6 - 
     2    (4*CA*CF)/((-1 + xa)*(-1 + xb)) + (-12*CA*CF - 12*CA*CF*xa)/(6*(-1 + xb))
     
      P0qqP0gq(2,1) = 0.0_ki
     
      P0qqP0qg(2,1) = 0.0_ki
     
      P0qqP0qq(2,1) = (-3*CF**2)/(-1 + xa) + (-3*CF**2 - 3*CF**2*xa)/2 - 
     1    (4*CF**2)/((-1 + xa)*(-1 + xb)) + (-4*CF**2 - 4*CF**2*xa)/(2*(-1 + xb))
     
      P0ggP0gg(1,2) = (-11*CA**2 + 2*CA*xlf)/(3*(-1 + xb)) - 
     1    (4*CA**2)/((-1 + xa)*(-1 + xb)) + (11*CA**2 - 2*CA*xlf)/(3*xb) + 
     2    (4*CA**2)/((-1 + xa)*xb) + (-24*CA**2 + 12*CA**2*xb - 12*CA**2*xb**2)/
     3    (3*(-1 + xa)) + (-22*CA**2 + 4*CA*xlf + 11*CA**2*xb - 2*CA*xlf*xb - 
     4    11*CA**2*xb**2 + 2*CA*xlf*xb**2)/3
       
      P0ggP0gq(1,2) = (22*CA*CF - 4*CF*xlf)/(6*xb) + (4*CA*CF)/((-1 + xa)*xb) + 
     1    (-24*CA*CF + 12*CA*CF*xb)/(6*(-1 + xa)) + 
     2    (-22*CA*CF + 4*CF*xlf + 11*CA*CF*xb - 2*CF*xlf*xb)/6
     
      P0ggP0qg(1,2) = (12*CA - 24*CA*xb + 24*CA*xb**2)/(12*(-1 + xa)) + 
     1    (11*CA - 2*xlf - 22*CA*xb + 4*xlf*xb + 22*CA*xb**2 - 4*xlf*xb**2)/12
     
      P0ggP0qq(1,2) = (-22*CA*CF + 4*CF*xlf)/(6*(-1 + xb)) - 
     1    (4*CA*CF)/((-1 + xa)*(-1 + xb)) + (-12*CA*CF - 12*CA*CF*xb)/
     2    (6*(-1 + xa)) + (-11*CA*CF + 2*CF*xlf - 11*CA*CF*xb + 2*CF*xlf*xb)/6
     
      P0gqP0gg(1,2) = 0.0_ki
     
      P0gqP0gq(1,2) = 0.0_ki
     
      P0gqP0qg(1,2) = 0.0_ki
     
      P0gqP0qq(1,2) = 0.0_ki
     
      P0qgP0gg(1,2) = 0.0_ki
     
      P0qgP0gq(1,2) = 0.0_ki
     
      P0qgP0qg(1,2) = 0.0_ki
     
      P0qgP0qq(1,2) = 0.0_ki
     
      P0qqP0gg(1,2) = -6*CA*CF - (3*CA*CF)/(-1 + xb) - 
     1    (4*CA*CF)/((-1 + xa)*(-1 + xb)) + (3*CA*CF)/xb + 
     2    (4*CA*CF)/((-1 + xa)*xb) + 3*CA*CF*xb - 3*CA*CF*xb**2 + 
     3    (-8*CA*CF + 4*CA*CF*xb - 4*CA*CF*xb**2)/(-1 + xa)
     
      P0qqP0gq(1,2) = (3*CF**2)/xb + (4*CF**2)/((-1 + xa)*xb) + 
     1    (-6*CF**2 + 3*CF**2*xb)/2 + (-8*CF**2 + 4*CF**2*xb)/(2*(-1 + xa))
     
      P0qqP0qg(1,2) = (3*CF - 6*CF*xb + 6*CF*xb**2)/4 + (4*CF - 8*CF*xb + 8*CF*xb**2)/
     1    (4*(-1 + xa))
     
      P0qqP0qq(1,2) = (-3*CF**2)/(-1 + xb) - (4*CF**2)/((-1 + xa)*(-1 + xb)) + 
     1    (-4*CF**2 - 4*CF**2*xb)/(2*(-1 + xa)) + (-3*CF**2 - 3*CF**2*xb)/2
     
      P0ggP0gg(2,2) = 16*CA**2 - 8*CA**2*xa + 8*CA**2*xa**2 + 
     1    (4*CA**2)/((-1 + xa)*(-1 + xb)) - (4*CA**2)/(xa*(-1 + xb)) + 
     2    (8*CA**2 - 4*CA**2*xa + 4*CA**2*xa**2)/(-1 + xb) - (4*CA**2)/((-1 + xa)*xb) + 
     3    (4*CA**2)/(xa*xb) + (-8*CA**2 + 4*CA**2*xa - 4*CA**2*xa**2)/xb - 8*CA**2*xb + 
     4    4*CA**2*xa*xb - 4*CA**2*xa**2*xb + 8*CA**2*xb**2 - 4*CA**2*xa*xb**2 + 
     5    4*CA**2*xa**2*xb**2 + (-8*CA**2 + 4*CA**2*xb - 4*CA**2*xb**2)/xa + 
     6    (8*CA**2 - 4*CA**2*xb + 4*CA**2*xb**2)/(-1 + xa)
     
      P0ggP0gq(2,2) = 8*CA*CF - 4*CA*CF*xa + 4*CA*CF*xa**2 - 
     1    (4*CA*CF)/((-1 + xa)*xb) + (4*CA*CF)/(xa*xb) + 
     2    (-8*CA*CF + 4*CA*CF*xa - 4*CA*CF*xa**2)/xb - 4*CA*CF*xb + 2*CA*CF*xa*xb - 
     3    2*CA*CF*xa**2*xb + (4*CA*CF - 2*CA*CF*xb)/(-1 + xa) + 
     4    (-4*CA*CF + 2*CA*CF*xb)/xa
     
      P0ggP0qg(2,2) = -2*CA + CA*xa - CA*xa**2 + 4*CA*xb - 2*CA*xa*xb + 2*CA*xa**2*xb - 
     1    4*CA*xb**2 + 2*CA*xa*xb**2 - 2*CA*xa**2*xb**2 + (-CA + 2*CA*xb - 2*CA*xb**2)/
     2    (-1 + xa) + (CA - 2*CA*xb + 2*CA*xb**2)/xa
     
      P0ggP0qq(2,2) = 4*CA*CF - 2*CA*CF*xa + 2*CA*CF*xa**2 + 
     1    (4*CA*CF)/((-1 + xa)*(-1 + xb)) - (4*CA*CF)/(xa*(-1 + xb)) + 
     2    (8*CA*CF - 4*CA*CF*xa + 4*CA*CF*xa**2)/(-1 + xb) + 4*CA*CF*xb - 
     3    2*CA*CF*xa*xb + 2*CA*CF*xa**2*xb + (-2*CA*CF - 2*CA*CF*xb)/xa + 
     4    (2*CA*CF + 2*CA*CF*xb)/(-1 + xa)
     
      P0gqP0gg(2,2) = 8*CA*CF - 4*CA*CF*xa - (4*CA*CF)/(xa*(-1 + xb)) + 
     1    (4*CA*CF - 2*CA*CF*xa)/(-1 + xb) + (4*CA*CF)/(xa*xb) + 
     2    (-4*CA*CF + 2*CA*CF*xa)/xb - 4*CA*CF*xb + 2*CA*CF*xa*xb + 4*CA*CF*xb**2 - 
     3    2*CA*CF*xa*xb**2 + (-8*CA*CF + 4*CA*CF*xb - 4*CA*CF*xb**2)/xa
     
      P0gqP0gq(2,2) = 4*CF**2 - 2*CF**2*xa + (4*CF**2)/(xa*xb) + 
     1    (-4*CF**2 + 2*CF**2*xa)/xb - 2*CF**2*xb + CF**2*xa*xb + 
     2    (-4*CF**2 + 2*CF**2*xb)/xa
     
      P0gqP0qg(2,2) = (2*CF - 4*CF*xb + 4*CF*xb**2)/(2*xa) + 
     1    (-2*CF + CF*xa + 4*CF*xb - 2*CF*xa*xb - 4*CF*xb**2 + 2*CF*xa*xb**2)/2
     
      P0gqP0qq(2,2) = 2*CF**2 - CF**2*xa - (4*CF**2)/(xa*(-1 + xb)) + 
     1    (4*CF**2 - 2*CF**2*xa)/(-1 + xb) + 2*CF**2*xb - CF**2*xa*xb + 
     2    (-2*CF**2 - 2*CF**2*xb)/xa
     
      P0qgP0gg(2,2) = -2*CA + 4*CA*xa - 4*CA*xa**2 + (-CA + 2*CA*xa - 2*CA*xa**2)/
     1    (-1 + xb) + (CA - 2*CA*xa + 2*CA*xa**2)/xb + CA*xb - 2*CA*xa*xb + 
     2    2*CA*xa**2*xb - CA*xb**2 + 2*CA*xa*xb**2 - 2*CA*xa**2*xb**2
     
      P0qgP0gq(2,2) = (2*CF - 4*CF*xa + 4*CF*xa**2)/(2*xb) + 
     1    (-2*CF + 4*CF*xa - 4*CF*xa**2 + CF*xb - 2*CF*xa*xb + 2*CF*xa**2*xb)/2
     
      P0qgP0qg(2,2) = 1/4.0_ki - xa/2 + xa**2/2 - xb/2 + xa*xb - xa**2*xb + xb**2/2 - 
     1    xa*xb**2 + xa**2*xb**2
     
      P0qgP0qq(2,2) = (-2*CF + 4*CF*xa - 4*CF*xa**2)/(2*(-1 + xb)) + 
     1    (-CF + 2*CF*xa - 2*CF*xa**2 - CF*xb + 2*CF*xa*xb - 2*CF*xa**2*xb)/2
     
      P0qqP0gg(2,2) = 4*CA*CF + 4*CA*CF*xa + (4*CA*CF)/((-1 + xa)*(-1 + xb)) + 
     1    (2*CA*CF + 2*CA*CF*xa)/(-1 + xb) - (4*CA*CF)/((-1 + xa)*xb) + 
     2    (-2*CA*CF - 2*CA*CF*xa)/xb - 2*CA*CF*xb - 2*CA*CF*xa*xb + 2*CA*CF*xb**2 + 
     3    2*CA*CF*xa*xb**2 + (8*CA*CF - 4*CA*CF*xb + 4*CA*CF*xb**2)/(-1 + xa)
     
      P0qqP0gq(2,2) = 2*CF**2 + 2*CF**2*xa - (4*CF**2)/((-1 + xa)*xb) + 
     1    (-2*CF**2 - 2*CF**2*xa)/xb - CF**2*xb - CF**2*xa*xb + 
     2    (4*CF**2 - 2*CF**2*xb)/(-1 + xa)
     
      P0qqP0qg(2,2) = (-2*CF + 4*CF*xb - 4*CF*xb**2)/(2*(-1 + xa)) + 
     1    (-CF - CF*xa + 2*CF*xb + 2*CF*xa*xb - 2*CF*xb**2 - 2*CF*xa*xb**2)/2
     
      P0qqP0qq(2,2) = CF**2 + CF**2*xa + (4*CF**2)/((-1 + xa)*(-1 + xb)) + 
     1    (2*CF**2 + 2*CF**2*xa)/(-1 + xb) + CF**2*xb + CF**2*xa*xb + 
     2    (2*CF**2 + 2*CF**2*xb)/(-1 + xa)

      return
      end

      
