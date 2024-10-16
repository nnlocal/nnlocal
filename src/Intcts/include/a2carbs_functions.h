!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      real(ki) z3
      real(ki) z4
      real(ki) z5
      complex(ki) l1
      complex(ki) l2
      complex(ki) l3
      complex(ki) l4
      complex(ki) l5
      complex(ki) l6
      complex(ki) l7
      complex(ki) l8
      complex(ki) d1
      z1 = (1/(-1 + xb))
      z2 = (1/xa)
      z3 = (1/(-1 + xa))
      z4 = (1/xb)
      z5 = (1/(-2 + xa + xb))
      l1 = log(2*cone)
      l2 = log(cone*(1 - xa))
      l3 = log(cone*(1 - xb))
      l4 = log(cone*(2 - xa - xb))
      l5 = log(cone*xb)
      l6 = log(cone*(1 + xb))
      l7 = log(cone*xa)
      l8 = log(cone*(1 + xa))
      d1 = cli2(-((cone*(-1 + xa))/(-1 + xb)))
