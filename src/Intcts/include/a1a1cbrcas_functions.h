!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      real(ki) z3
      real(ki) z4
      real(ki) z5
      real(ki) z6
      real(ki) z7
      real(ki) z8
      real(ki) z9
      real(ki) z10
      real(ki) z11
      real(ki) z12
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
      complex(ki) l12
      complex(ki) l13
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
      complex(ki) d16
      complex(ki) d17
      complex(ki) d18
      complex(ki) d19
      complex(ki) d20
      complex(ki) d21
      complex(ki) d22
      complex(ki) d23
      complex(ki) d24
      complex(ki) d25
      z1 = (1/xa)
      z2 = (1/(xa - xb))
      z3 = (1/xb)
      z4 = (1/(-1 + xa))
      z5 = (1/(1 + xa))
      z6 = (1/(1 + xa**2))
      z7 = (1/(xa + xb))
      z8 = (1/(-1 + xb))
      z9 = (1/(1 + xb))
      z10 = (1/(1 + xb**2))
      z11 = (1/(-1 + xa*xb))
      z12 = (1/(-2 + xa + xb))
      l1 = log(2*cone)
      l2 = log(cone*(1 - xa))
      l3 = log(cone*(1 - xb))
      l4 = log(cone*(2 - xa - xb))
      l5 = log(cone*(1 + xb))
      l6 = log(cone*xa)
      l7 = log(cone*(1 + xa))
      l8 = log(cone*(1 + xa**2))
      l9 = log(cone*xb)
      l10 = log(cone*(xa + xb))
      l11 = log(cone*(1 - xa*xb))
      l12 = log(cone*(1 + xa*xb))
      l13 = log(cone*(1 + xb**2))
      d1 = cli2(-((cone*(1 - xb))/(1 - xa)))
      d2 = cli2((cone*(1 - xa))/2.0_ki)
      d3 = cli2(-(cone*xa))
      d4 = cli2(cone*xa)
      d5 = cli2((-2*cone*xa)/(1 - xa)**2)
      d6 = cli2(cone*xa**4)
      d7 = cli2((2*cone*xa)/(1 + xa)**2)
      d8 = cli2((cone*(1 - xb))/2.0_ki)
      d9 = cli2(-((cone*xa)/xb))
      d10 = cli2(-(cone*xb))
      d11 = cli2(cone*xb)
      d12 = cli2(cone*xa*xb)
      d13 = cli2((-2*cone*xb)/(1 - xb)**2)
      d14 = cli2(-((cone*(1 + xa)*xb)/(1 - xb)))
      d15 = cli2(cone*xb**4)
      d16 = cli2((2*cone*xb)/(1 + xb)**2)
      d17 = cli2(-((cone*xa*(1 + xb))/(1 - xa)))
      d18 = cli2(-((cone*(-xa + xb))/(1 + xa)))
      d19 = cli2((cone*(-xa + xb))/xb)
      d20 = cli2(-((cone*(1 - xa))/(xa + xb)))
      d21 = cli2(-((cone*(1 - xb))/(xa + xb)))
      d22 = cli2((cone*(1 + xa**2)*xb)/(xa + xb))
      d23 = cli2((cone*(1 - xb))/(1 - xa*xb))
      d24 = cli2(-((cone*xb*(xa + xb))/(1 - xa*xb)))
      d25 = cli2((cone*(1 - xa*xb))/(1 + xb))
