!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      complex(ki) l1
      complex(ki) l2
      complex(ki) l3
      complex(ki) l4
      complex(ki) d1
      complex(ki) d2
      complex(ki) d3
      complex(ki) d4
      complex(ki) d5
      complex(ki) t1
      complex(ki) t2
      complex(ki) t3
      z1 = (1/(-1 + lam))
      z2 = (1/(1 + lam))
      l1 = log(2*cone)
      l2 = log(cone*(1 - lam))
      l3 = log(cone*lam)
      l4 = log(cone*(1 + lam))
      d1 = hp_cli2((cone*(1 - lam))/2.0_ki)
      d2 = hp_cli2(cone*(1 - lam))
      d3 = hp_cli2(-(cone*lam))
      d4 = hp_cli2(cone*lam)
      d5 = hp_cli2((cone*(1 + lam))/2.0_ki)
      t1 = hp_Li3(-lam)
      t2 = hp_Li3(lam)
      t3 = hp_Li3(lam**2)
