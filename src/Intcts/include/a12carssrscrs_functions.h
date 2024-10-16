!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      complex(ki) r1
      complex(ki) l1
      complex(ki) l2
      complex(ki) l3
      complex(ki) l4
      complex(ki) l5
      complex(ki) d1
      complex(ki) d2
      complex(ki) d3
      complex(ki) d4
      complex(ki) d5
      complex(ki) t1
      complex(ki) t2
      complex(ki) t3
      complex(ki) t4
      z1 = (1/(-1 + lam))
      z2 = (1/(1 + lam))
      r1 = Sqrt(cone*(1 - lam**2))
      l1 = log(2*cone)
      l2 = log(cone*(1 - lam))
      l3 = log(cone*lam)
      l4 = log(cone*(1 + lam))
      l5 = log(cone*(1 - Sqrt(1 - lam**2)))
      d1 = cli2(cone*lam)
      d2 = cli2((cone*(-1 + lam))/(1 + lam))
      d3 = cli2(-(cone*Sqrt(1 - lam**2)))
      d4 = cli2((cone*(1 - lam - Sqrt(1 - lam**2)))/2.0_ki)
      d5 = cli2(-(cone*lam))
      t1 = Li3(1 - lam)
      t2 = Li3(-lam)
      t3 = Li3(lam)
      t4 = Li3(1/(1 + lam))
