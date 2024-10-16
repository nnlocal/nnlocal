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
      complex(ki) r1
      complex(ki) r2
      complex(ki) r3
      complex(ki) r4
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
      complex(ki) l14
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
      z1 = (1/Sqrt(1 - xa*xb))
      z2 = (1/(-1 + xa))
      z3 = (1/(1 + xa))
      z4 = (1/(-1 + xb))
      z5 = (1/(1 + xb))
      z6 = (1/(-2 + xa + xb))
      z8 = (1/(2 + xa + xb))
      z7 = (1/(xa + xb))
      z9 = (1/(-1 + xa*xb))
      z11 = (1/(xa - xb))
      z10 = (1/xb)
      z12 = (1/xa)
      r1 = Sqrt(cone*xa*xb)
      r2 = Sqrt(cone*xa)
      r3 = Sqrt(cone*xb)
      r4 = Sqrt(cone*(1 - xa*xb))
      l1 = log(2*cone)
      l2 = log(cone*(1 - xa))
      l3 = log(cone*(1 - xb))
      l4 = log(cone*(2 - xa - xb))
      l5 = log(cone*(xa - xb))
      l6 = log(cone*xa)
      l7 = log(cone*(1 + xa))
      l8 = log((cone*(Sqrt(xa) - Sqrt(xb)))/(Sqrt(xa) + Sqrt(xb)))
      l9 = log((cone*(1 - Sqrt(xa)*Sqrt(xb)))/(1 + Sqrt(xa)*Sqrt(xb)))
      l10 = log(cone*xb)
      l11 = log(cone*(1 + xb))
      l12 = log(cone*(xa + xb))
      l13 = log(cone*(1 - xa*xb))
      l14 = log((cone*(-1 + xa + Sqrt(1 - xa*xb)))/(1 - xa + Sqrt(1 - xa*xb)))
      d1 = hp_cli2((cone*(1 - xb))/(-1 + xa))
      d2 = hp_cli2((cone*(xa - xb))/(-1 + xa))
      d3 = hp_cli2(-(cone*xa))
      d4 = hp_cli2(cone*xa)
      d5 = hp_cli2((cone*(1 - xa))/(1 + Sqrt(xa)*Sqrt(xb)))
      d6 = hp_cli2((cone*(1 - xa))/(1 - xb))
      d7 = hp_cli2((cone*Sqrt(xb))/Sqrt(xa))
      d8 = hp_cli2(-(cone*xb))
      d9 = hp_cli2(cone*xb)
      d10 = hp_cli2((cone*xb)/xa)
      d11 = hp_cli2(cone*xa*xb)
      d12 = hp_cli2((cone*(xa + xb))/(1 + xa))
      d13 = hp_cli2((cone*(xa + xb))/(1 + xb))
      d14 = hp_cli2((cone*(1 - xa))/(1 - xa*xb))
      d15 = hp_cli2((cone*(1 - xb))/(1 - xa*xb))
      d16 = hp_cli2((cone*(1 - xa*xb))/(1 + xb))
      d17 = hp_cli2((cone*(1 - xa)*xa)/(xa - xa*xb))
      d18 = hp_cli2((cone*(xa - xa*xb))/(xa + Sqrt(xa)*Sqrt(xb)))
