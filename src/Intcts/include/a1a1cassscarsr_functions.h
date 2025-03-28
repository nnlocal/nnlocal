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
      z1 = (1/(-1 + lam))
      z2 = (1/(1 + lam))
      l1 = log(2*cone)
      l2 = log(cone*(1 - lam))
      l3 = log(cone*lam)
      l4 = log(cone*(1 + lam))
      d1 = cli2(-(cone*lam))
      d2 = cli2(cone*lam)
      t1 = Li3(-lam)
      t2 = Li3(lam)
