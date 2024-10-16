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
      complex(ki) d1
      complex(ki) d2
      complex(ki) d3
      complex(ki) d4
      complex(ki) d5
      complex(ki) d6
      complex(ki) t1
      complex(ki) t2
      complex(ki) t3
      complex(ki) t4
      complex(ki) t5
      complex(ki) t6
      complex(ki) t7
      complex(ki) t8
      complex(ki) t9
      complex(ki) t10
      complex(ki) t11
      z1 = (1/(-1 + xa))
      z2 = (1/(-1 + xb))
      z3 = (1/(-2 + xa + xb))
      z4 = (1/xa)
      z5 = (1/(1 + xa))
      l1 = log(2*cone)
      l2 = log(cone*(1 - xa))
      l3 = log(cone*(1 - xb))
      l4 = log(cone*(2 - xa - xb))
      l5 = log(cone*xa)
      l6 = log(cone*(1 + xa))
      l7 = log(cone*(1 + 2*xa))
      d1 = cli2(-((cone*(-1 + xb))/(-1 + xa)))
      d2 = cli2(cone*(-1 - 2*xa))
      d3 = cli2((cone*(1 - xa))/2.0_ki)
      d4 = cli2(-2*cone*xa)
      d5 = cli2(-(cone*xa))
      d6 = cli2(cone*xa)
      t1 = Li3((1 - xa)/2.0_ki)
      t2 = Li3(1 - xa)
      t3 = Li3(-2*xa)
      t4 = Li3(-xa)
      t5 = Li3(xa)
      t6 = Li3((2*xa)/(-1 + xa))
      t7 = Li3(1/(2.0_ki*(1 + xa)))
      t8 = Li3(1/(1 + xa))
      t9 = Li3(-((-1 + xa)/(1 + xa)))
      t10 = Li3(-(xa/(1 + xa)))
      t11 = Li3((1 + xa)/2.0_ki)
