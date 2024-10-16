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
      complex(ki) l7
      complex(ki) l8
      complex(ki) l9
      complex(ki) l10
      complex(ki) d1
      complex(ki) d2
      complex(ki) d3
      complex(ki) d4
      complex(ki) d5
      complex(ki) t1
      complex(ki) t2
      complex(ki) t3
      complex(ki) t4
      z1 = (1/(-1 + xb))
      z2 = (1/(1 + xb))
      z3 = (1/(xa + xb))
      z4 = (1/xa)
      z5 = (1/(-1 + xa*xb))
      z6 = (1/(-1 + xa))
      z7 = (1/(-2 + xa + xb))
      l1 = log(2*cone)
      l2 = log(cone*(1 - xa))
      l3 = log(cone*(1 - xb))
      l4 = log(cone*(2 - xa - xb))
      l5 = log(cone*xa)
      l6 = log(cone*(1 + xa))
      l7 = log(cone*xb)
      l8 = log(cone*(1 + xb))
      l9 = log(cone*(xa + xb))
      l10 = log(cone*(1 - xa*xb))
      d1 = cli2((cone*(1 - xa))/2.0_ki)
      d2 = cli2(cone*(1 - xa))
      d3 = cli2(cone*xa)
      d4 = cli2((cone*(1 + xa))/2.0_ki)
      d5 = cli2(cone*(1 - 1/(xa*xb)))
      t1 = Li3(1 - xa)
      t2 = Li3(-xa)
      t3 = Li3(xa)
      t4 = Li3(xa**2)
