!--- include file for double soft ct's
      double precision z1
      double precision z2
      double precision z3
      double precision z4
      double precision z5
      double precision z6
      double precision z7
      double complex l1
      double complex l2
      double complex l3
      double complex l4
      double complex l5
      double complex l6
      z1 = (1/(-1 + xb))
      z2 = (1/(1 + xb))
      z3 = (1/(xa + xb))
      z4 = (1/xa)
      z5 = (1/(-1 + xa*xb))
      z6 = (1/(-1 + xa))
      z7 = (1/(-2 + xa + xb))
      l1 = log(cone*(1 - xa))
      l2 = log(cone*(1 - xb))
      l3 = log(2*cone)
      l4 = log(cone*(1 + xa))
      l5 = log(cone*(1 + xb))
      l6 = log(cone*(xa + xb))
