!--- include file for double soft ct's
      real(ki) z1
      real(ki) z2
      complex(ki) l1
      complex(ki) l2
      complex(ki) l3
      complex(ki) l4
      complex(ki) d1
      complex(ki) d2
      complex(ki) t1
      complex(ki) t2
      complex(ki) t3
      complex(ki) t4
      complex(ki) t5
      complex(ki) t6
      complex(ki) t7
      z1 = (1/(-1 + lam))
      z2 = (1/(1 + lam))
      l1 = log(2*cone)
      l2 = log(cone*(1 - lam))
      l3 = log(cone*lam)
      l4 = log(cone*(1 + lam))
      d1 = hp_cli2(-(cone*lam))
      d2 = hp_cli2(cone*lam)
      t1 = hp_Li3((1 - lam)/2.0_ki)
      t2 = hp_Li3(1 - lam)
      t3 = hp_Li3(-lam)
      t4 = hp_Li3(lam)
      t5 = hp_Li3(1/(1 + lam))
      t6 = hp_Li3((1 - lam)/(1 + lam))
      t7 = hp_Li3((1 + lam)/2.0_ki)
