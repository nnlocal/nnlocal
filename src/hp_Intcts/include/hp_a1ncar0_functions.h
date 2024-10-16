!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      real(ki) z3
      real(ki) z4
      real(ki) z5
      real(ki) z6
      real(ki) z7
      complex(ki) l1
      complex(ki) l2
      complex(ki) l3
      complex(ki) l4
      complex(ki) l5
      complex(ki) l6
      z1 = (1/(-1 + xb))
      z2 = (1/(1 + xb))
      z3 = (1/(xa + xb))
      z4 = (1/xa)
      z5 = (1/(-1 + xa*xb))
      z6 = (1/(-1 + xa))
      z7 = (1/(-2 + xa + xb))
      l1 = log(cone*(1 - xa))
      l2 = log(cone*(1 - xb))
      l3 = log(2*cone)
      l4 = log(cone*(1 + xa))
      l5 = log(cone*(1 + xb))
      l6 = log(cone*(xa + xb))
