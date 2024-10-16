!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      real(ki) z3
      real(ki) z4
      real(ki) z5
      complex(ki) z6
      complex(ki) z7
      real(ki) z8
      real(ki) z9
      real(ki) z10
      real(ki) z11
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
      complex(ki) l11
      complex(ki) d1
      complex(ki) d2
      complex(ki) d3
      complex(ki) d4
      complex(ki) d5
      complex(ki) d6
      complex(ki) d7
      complex(ki) d8
      complex(ki) d9
      complex(ki) d10
      complex(ki) d11
      complex(ki) d12
      complex(ki) d13
      complex(ki) d14
      complex(ki) d15
      z1 = (1/xa)
      z2 = (1/(-1 + xb))
      z3 = (1/(1 + xb))
      z4 = (1/(xa + xb))
      z5 = (1/(-1 + xa))
      z6 = (1/(-im + xa))
      z7 = (1/(im + xa))
      z8 = (1/(1 + xa))
      z9 = (1/(-1 + xa*xb))
      z10 = (1/(-2 + xa + xb))
      z11 = (1/xb)
      l1 = log(cone*(1 - xa))
      l2 = log(2*cone)
      l3 = log(cone*xa)
      l4 = log(cone*(1 + xa))
      l5 = log(cone*(1 + xa**2))
      l6 = log(cone*(1 - xb))
      l7 = log(cone*xb)
      l8 = log(cone*(1 + xb))
      l9 = log(cone*(xa + xb))
      l10 = log(cone*(1 - xa*xb))
      l11 = log(cone*(2 - xa - xb))
      d1 = hp_cli2((cone*(1 - xa))/2.0_ki)
      d2 = hp_cli2(-(cone*xa))
      d3 = hp_cli2(cone*xa)
      d4 = hp_cli2(-(cone*xa**2))
      d5 = hp_cli2(-((cone*(-1 + xa))/(1 + xa**2)))
      d6 = hp_cli2((cone*(1 - xb))/2.0_ki)
      d7 = hp_cli2(-((cone*(-1 + xb))/(1 + xa)))
      d8 = hp_cli2(-((cone*xa)/xb))
      d9 = hp_cli2(-(cone*xb))
      d10 = hp_cli2(cone*xb)
      d11 = hp_cli2(-((cone*(-1 + xa))/(1 + xb)))
      d12 = hp_cli2((cone*(xa + xb))/((1 + xa**2)*xb))
      d13 = hp_cli2((cone*(-1 + xb))/(-1 + xa*xb))
      d14 = hp_cli2(-((cone*(-1 + xa*xb))/(1 + xb)))
      d15 = hp_cli2(-((cone*(-1 + xb))/(-1 + xa)))
