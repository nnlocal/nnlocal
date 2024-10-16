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
      complex(ki) t1
      complex(ki) t2
      complex(ki) t3
      complex(ki) t4
      complex(ki) t5
      complex(ki) t6
      complex(ki) t7
      complex(ki) t8
      z1 = (1/xa)
      z2 = (1/(xa - xb))
      z3 = (1/(-1 + xb))
      z4 = (1/xb)
      z5 = (1/(1 + xb))
      z6 = (1/(xa + xb))
      z7 = (1/(1 + xb**2))
      z8 = (1/(-1 + xa*xb))
      z9 = (1/(-1 + xa))
      z10 = (1/(-2 + xa + xb))
      l1 = log(2*cone)
      l2 = log(cone*(1 - xa))
      l3 = log(cone*(1 - xb))
      l4 = log(cone*(2 - xa - xb))
      l5 = log(cone*xa)
      l6 = log(cone*(1 + xa))
      l7 = log(cone*(1 + xa**2))
      l8 = log(cone*xb)
      l9 = log(cone*(1 + xb))
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
      d9 = cli2((cone*(xa - xb))/(1 + xa))
      d10 = cli2(-((cone*xa)/xb))
      d11 = cli2(-((cone*(xa - xb))/xb))
      d12 = cli2(-(cone*xb))
      d13 = cli2(cone*xb)
      d14 = cli2(cone*xa*xb)
      d15 = cli2((-2*cone*xb)/(1 - xb)**2)
      d16 = cli2(-((cone*(1 + xa)*xb)/(1 - xb)))
      d17 = cli2(cone*xb**4)
      d18 = cli2((2*cone*xb)/(1 + xb)**2)
      d19 = cli2(-((cone*xa*(1 + xb))/(1 - xa)))
      d20 = cli2(-((cone*(1 - xa))/(xa + xb)))
      d21 = cli2(-((cone*(1 - xb))/(xa + xb)))
      d22 = cli2((cone*(1 + xa**2)*xb)/(xa + xb))
      d23 = cli2((cone*(1 - xb))/(1 - xa*xb))
      d24 = cli2(-((cone*xb*(xa + xb))/(1 - xa*xb)))
      d25 = cli2((cone*(1 - xa*xb))/(1 + xb))
      t1 = Li3((1 - xa)/2.0_ki)
      t2 = Li3(-xa)
      t3 = Li3(xa)
      t4 = Li3((2*xa)/(-1 + xa))
      t5 = Li3((1 - xa)/(1 + xa))
      t6 = Li3((1 + xa)/2.0_ki)
      t7 = Li3(1 - xa)
      t8 = Li3(1/(1 + xa))
