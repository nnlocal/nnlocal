!--- include file for double soft ct's
      complex(ki) r1
      complex(ki) r2
      complex(ki) r3
      complex(ki) r4
      complex(ki) r5
      complex(ki) z1
      complex(ki) z2
      complex(ki) z3
      complex(ki) z4
      complex(ki) z5
      complex(ki) z6
      complex(ki) z7
      complex(ki) z8
      complex(ki) z9
      complex(ki) z10
      complex(ki) z11
      complex(ki) z12
      complex(ki) z13
      complex(ki) z14
      complex(ki) z15
      complex(ki) z16
      complex(ki) z17
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
      complex(ki) l15
      complex(ki) l16
      complex(ki) l17
      complex(ki) l18
      complex(ki) l19
      complex(ki) l20
      complex(ki) l21
      complex(ki) l22
      complex(ki) l23
      complex(ki) l24
      complex(ki) l25
      complex(ki) l26
      complex(ki) l27
      complex(ki) l28
      complex(ki) l29
      complex(ki) l30
      complex(ki) l31
      complex(ki) g1
      complex(ki) g2
      complex(ki) g3
      complex(ki) g4
      complex(ki) g5
      complex(ki) g6
      complex(ki) g7
      complex(ki) g8
      complex(ki) g9
      complex(ki) g10
      complex(ki) g11
      complex(ki) g12
      complex(ki) g13
      complex(ki) g14
      complex(ki) g15
      complex(ki) g16
      complex(ki) g17
      complex(ki) g18
      complex(ki) g19
      complex(ki) g20
      complex(ki) g21
      complex(ki) g22
      complex(ki) g23
      complex(ki) g24
      complex(ki) g25
      complex(ki) g26
      complex(ki) g27
      complex(ki) g28
      complex(ki) g29
      complex(ki) g30
      complex(ki) g31
      complex(ki) g32
      complex(ki) g33
      complex(ki) g34
      complex(ki) g35
      complex(ki) g36
      complex(ki) g37
      complex(ki) g38
      complex(ki) g39
      complex(ki) g40
      complex(ki) g41
      complex(ki) g42
      complex(ki) g43
      complex(ki) g44
      complex(ki) g45
      complex(ki) g46
      complex(ki) g47
      complex(ki) g48
      complex(ki) g49
      complex(ki) g50
      complex(ki) g51
      complex(ki) g52
      complex(ki) g53
      complex(ki) g54
      complex(ki) g55
      complex(ki) g56
      complex(ki) g57
      complex(ki) g58
      complex(ki) g59
      complex(ki) g60
      complex(ki) g61
      complex(ki) g62
      complex(ki) g63
      complex(ki) g64
      complex(ki) g65
      complex(ki) g66
      complex(ki) g67
      complex(ki) g68
      complex(ki) g69
      complex(ki) g70
      complex(ki) g71
      complex(ki) g72
      complex(ki) g73
      complex(ki) g74
      complex(ki) g75
      complex(ki) g76
      complex(ki) g77
      complex(ki) g78
      complex(ki) g79
      complex(ki) g80
      complex(ki) g81
      complex(ki) g82
      complex(ki) g83
      complex(ki) g84
      complex(ki) g85
      complex(ki) g86
      complex(ki) g87
      complex(ki) g88
      complex(ki) g89
      complex(ki) g90
      complex(ki) g91
      complex(ki) g92
      complex(ki) g93
      complex(ki) g94
      complex(ki) g95
      complex(ki) g96
      complex(ki) g97
      complex(ki) g98
      complex(ki) g99
      complex(ki) g100
      complex(ki) g101
      complex(ki) g102
      complex(ki) g103
      complex(ki) g104
      complex(ki) g105
      complex(ki) g106
      complex(ki) g107
      complex(ki) g108
      complex(ki) g109
      complex(ki) g110
      complex(ki) g111
      complex(ki) g112
      complex(ki) g113
      complex(ki) g114
      complex(ki) g115
      complex(ki) g116
      complex(ki) g117
      complex(ki) g118
      complex(ki) g119
      complex(ki) g120
      complex(ki) g121
      complex(ki) g122
      complex(ki) g123
      complex(ki) g124
      complex(ki) g125
      complex(ki) g126
      complex(ki) g127
      complex(ki) g128
      complex(ki) g129
      complex(ki) g130
      complex(ki) g131
      complex(ki) g132
      complex(ki) g133
      complex(ki) g134
      complex(ki) g135
      complex(ki) g136
      complex(ki) g137
      complex(ki) g138
      complex(ki) g139
      complex(ki) g140
      complex(ki) g141
      complex(ki) g142
      complex(ki) g143
      complex(ki) g144
      complex(ki) g145
      complex(ki) g146
      complex(ki) g147
      complex(ki) g148
      complex(ki) g149
      complex(ki) hp_mylog
      r1 = Sqrt(cone*(1 + 4*xa**2 + 4*xa*xb))
      r2 = Sqrt(cone*xa*xb*(-1 + xa*xb + xb**2))
      r3 = Sqrt(cone*xa)
      r4 = Sqrt(cone*xb)
      r5 = Sqrt(cone*(xa**2*xb + xb**3 + 2*xa*(-2 + xb**2)))
      z1 = 1/r5
      z7 = 1/(-1 + xb*(xa + xb))
      z15 = 1/(-2 + xa + xb)
      z2 = 1/(1 + xb)
      z4 = 1/(-1 + xb)
      z5 = 1/(xa + xb)
      z6 = 1/(-1 + xa*xb)
      z8 = 1/(-1 + xa)
      z9 = 1/(xa - xb)
      z11 = 1/(-2 + xb*(xa + xb))
      z12 = 1/r1
      z13 = 1/(1 + xa)
      z14 = 1/(9*xa - xb)
      z3 = 1/xa
      z10 = 1/xb
      l1 = hp_mylog(2*cone)
      l2 = hp_mylog(cone*(1 - xa))
      l3 = hp_mylog(cone*xa)
      l4 = hp_mylog(cone*(1 + xa))
      l5 = hp_mylog(cone*(1 - xb))
      l6 = hp_mylog(cone*(2 - xb))
      l7 = hp_mylog(cone*(1 - xa/xb))
      l8 = hp_mylog(cone*(2 - xa - xb))
      l9 = hp_mylog(cone*xb)
      l10 = hp_mylog(cone*(1 + xb))
      l11 = hp_mylog(cone*(2 + xb))
      l12 = hp_mylog(cone*(xa + xb))
      l13 = hp_mylog(cone*(2 + xa + xb))
      l14 = hp_mylog((cone*(1 + r1 - 2*xa - 2*xb))/(-1 + r1 - 2*xa*xb))
      l15 = hp_mylog((cone*(-1 + r1 - 2*xa*xb))/(2.0_ki*(-1 + xa)*(-1 + xb)))
      l16 = hp_mylog(cone*(1 - xa*xb))
      l17 = hp_mylog(cone*(1 + xa*xb))
      l18 = hp_mylog((cone*(-1 + r1 + 2*xa + 2*xb))/(1 + r1 + 2*xa*xb))
      l19 = hp_mylog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb))
      l20 = hp_mylog((cone*(-2 - r4*r5 + xa*xb + xb**2))/(-(r4*r5) + xa*xb + xb**2))
      l21 = hp_mylog((cone*(-2 + r4*r5 + xa*xb + xb**2))/(r4*r5 + xa*xb + xb**2))
      l22 = hp_mylog((cone*(r2*(1 - xa) - (-1 + xa)**2*xb - (-1 + xa)*xb**2))/(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))
      l23 = hp_mylog(-((cone*(-1 + xa)**2*(-1 + xb)*xb)/(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2)))
      l24 = hp_mylog((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb))
      l25 = hp_mylog((cone*(r2*(1 - xa) + (-1 + xa)**2*xb + (-1 + xa)*xb**2))/(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))
      l26 = hp_mylog((cone*(-1 + xa)**2*(-1 + xb)*xb)/(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))
      l27 = hp_mylog(-((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)))
      l28 = hp_mylog(cone*(1 - xa - xb + 2*xa*xb - 2*xa**2*xb**2 + xa**3*xb**3))
      l29 = hp_mylog(cone*(2 - xb + xa**2*xb**2 - xa*(1 + xb)))
      l30 = hp_mylog((cone*(4*r4*xa + r5*xa - r4*xa**2 + r5*xb - 2*r4*xa*xb - xb**2.5_ki + (-r5 + r4*xa + r4*xb)*Abs(xa - xb)))/((r5 - r4*xa - r4*xb)*(xa + xb - Abs(xa - xb))))
      l31 = hp_mylog((cone*(-4*r4*xa + r5*xa + r4*xa**2 + r5*xb + 2*r4*xa*xb + xb**2.5_ki - (r5 + r4*xa + r4*xb)*Abs(xa - xb)))/((r5 + r4*xa + r4*xb)*(xa + xb - Abs(xa - xb))))
      g1 = hp_glog(-2*cone,-cone,cone*xb)
      g2 = hp_glog(-2*cone,czip,cone*xb)
      g3 = hp_glog(-2*cone,2*cone,cone*xb)
      g4 = hp_glog(-cone,-cone,cone*xa)
      g5 = hp_glog(-cone,-cone,cone*xb)
      g6 = hp_glog(-cone,czip,cone*xa)
      g7 = hp_glog(-cone,czip,cone*xb)
      g8 = hp_glog(-cone,cone,cone*xa)
      g9 = hp_glog(-cone,cone,cone*xb)
      g10 = hp_glog(-cone,2*cone,cone*xb)
      g11 = hp_glog(-cone,-(cone*xa),cone*xb)
      g12 = hp_glog(-cone,cone*xa,-(cone*xa*xb))
      g13 = hp_glog(-cone,cone*(2 - xb),cone*xa)
      g14 = hp_glog(-cone,-(cone/xb),cone*xa)
      g15 = hp_glog(-cone,cone/xb,cone*xa)
      g16 = hp_glog(-cone,-(cone*xb),cone*xa)
      g17 = hp_glog(czip,-cone,cone*xa)
      g18 = hp_glog(czip,-cone,cone*xb)
      g19 = hp_glog(czip,czip,2*cone)
      g20 = hp_glog(czip,czip,cone*xa)
      g21 = hp_glog(czip,czip,cone*xb)
      g22 = hp_glog(czip,czip,cone*(1 - xa*xb))
      g23 = hp_glog(czip,cone,cone*xa)
      g24 = hp_glog(czip,cone,cone*xb)
      g25 = hp_glog(czip,cone,cone*xa*xb)
      g26 = hp_glog(czip,2*cone,cone*xb)
      g27 = hp_glog(czip,-(cone*xa),cone*xb)
      g28 = hp_glog(czip,cone*(2 - xb),cone*xa)
      g29 = hp_glog(czip,-(cone/xb),cone*xa)
      g30 = hp_glog(czip,cone/xb,cone*xa)
      g31 = hp_glog(czip,-(cone*xb),cone*xa)
      g32 = hp_glog(czip,(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g33 = hp_glog(czip,(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g34 = hp_glog(czip,(cone*(-1 + xa + xb - xa*xb))/(-1 + xa*xb)**2,cone)
      g35 = hp_glog(czip,-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g36 = hp_glog(czip,(cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/(r2*(1 - xa) - (-1 + xa)**2*xb - (-1 + xa)*xb**2),cone)
      g37 = hp_glog(czip,(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g38 = hp_glog(czip,(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/(r2*(1 - xa) + (-1 + xa)**2*xb + (-1 + xa)*xb**2),cone)
      g39 = hp_glog(czip,(cone*(2 - xb + xa**2*xb**2 - xa*(1 + xb)))/(-1 + xa*xb)**2,cone)
      g40 = hp_glog(czip,(cone*(2 - xb + xa**2*xb**2 - xa*(1 + xb)))/(-1 + xa*xb)**2,cone*(1 - xa*xb))
      g41 = hp_glog(czip,(cone*(-1 + xa*xb)**3)/(1 - xb - 2*xa**2*xb**2 + xa**3*xb**3 + xa*(-1 + 2*xb)),cone)
      g42 = hp_glog(cone,-cone,cone*xa)
      g43 = hp_glog(cone,-cone,cone*xb)
      g44 = hp_glog(cone,czip,cone*xa)
      g45 = hp_glog(cone,czip,cone*xb)
      g46 = hp_glog(cone,cone,cone*xa)
      g47 = hp_glog(cone,cone,cone*xb)
      g48 = hp_glog(cone,cone,cone*xa*xb)
      g49 = hp_glog(cone,2*cone,cone*xb)
      g50 = hp_glog(cone,-(cone*xa),cone*xb)
      g51 = hp_glog(cone,-(cone*xa),cone*xa*xb)
      g52 = hp_glog(cone,cone*xa,cone*xa*xb)
      g53 = hp_glog(cone,cone*(2 - xb),cone*xa)
      g54 = hp_glog(cone,-(cone/xb),cone*xa)
      g55 = hp_glog(cone,cone/xb,cone*xa)
      g56 = hp_glog(cone,-(cone*xb),cone*xa)
      g57 = hp_glog(2*cone,-cone,cone*xb)
      g58 = hp_glog(2*cone,czip,cone*xb)
      g59 = hp_glog(2*cone,2*cone,cone*xb)
      g60 = hp_glog(-(cone*xa),cone,cone*xb)
      g61 = hp_glog(-(cone*xa),cone,cone*xa*xb)
      g62 = hp_glog(cone*xa,cone,cone*xa*xb)
      g63 = hp_glog(cone*(-2 - xb),-cone,cone*xa)
      g64 = hp_glog(cone*(-2 - xb),cone*(2 - xb),cone*xa)
      g65 = hp_glog(cone*(-2 - xb),cone/xb,cone*xa)
      g66 = hp_glog(cone*(-2 - xb),-(cone*xb),cone*xa)
      g67 = hp_glog(cone*(2 - xb),-cone,cone*xa)
      g68 = hp_glog(cone*(2 - xb),czip,cone*xa)
      g69 = hp_glog(cone*(2 - xb),cone*(2 - xb),cone*xa)
      g70 = hp_glog(cone*(2 - xb),cone/xb,cone*xa)
      g71 = hp_glog(cone*(2 - xb),-(cone*xb),cone*xa)
      g72 = hp_glog(-(cone/xb),-(cone/xb),cone*xa)
      g73 = hp_glog(-(cone/xb),-(cone*xb),cone*xa)
      g74 = hp_glog(cone/xb,czip,cone*xa)
      g75 = hp_glog(cone/xb,cone*(2 - xb),cone*xa)
      g76 = hp_glog(cone/xb,cone/xb,cone*xa)
      g77 = hp_glog(-(cone*xb),czip,cone*xa)
      g78 = hp_glog(-(cone*xb),cone,cone*xa)
      g79 = hp_glog(-(cone*xb),cone,cone*xa*xb)
      g80 = hp_glog(-(cone*xb),-(cone/xb),cone*xa)
      g81 = hp_glog(-(cone*xb),-(cone*xb),cone*xa)
      g82 = hp_glog(cone*xb,czip,cone*xa)
      g83 = hp_glog(cone*xb,cone,cone*xa*xb)
      g84 = hp_glog(cone*xb,cone*(2 - xb),cone*xa)
      g85 = hp_glog(cone*xb,cone/xb,cone*xa)
      g86 = hp_glog(cone*xb,-(cone*xb),cone*xa)
      g87 = hp_glog((cone*(1 + xb))/(1 - xa),czip,cone)
      g88 = hp_glog((cone*(1 + xb))/(1 - xa),cone,cone)
      g89 = hp_glog((cone*(1 + xb))/(1 - xa),(cone*(1 + xb))/(1 - xa),cone)
      g90 = hp_glog((cone*(1 + xb))/(1 - xa),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g91 = hp_glog((cone*(1 + xb))/(1 - xa),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g92 = hp_glog((cone*(1 + xb))/(1 - xa),(cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g93 = hp_glog((cone*(1 + xb))/(1 - xa),(cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g94 = hp_glog((cone*(1 + xb))/(1 - xa),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g95 = hp_glog((cone*(1 + xb))/(1 - xa),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g96 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),czip,cone)
      g97 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),cone,cone)
      g98 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g99 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g100 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g101 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g102 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g103 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g104 = hp_glog((cone*(-r5 + r4*xa + r4*xb))/(2.0_ki*r3),czip,cone/(r3*r4))
      g105 = hp_glog((cone*(-r5 + r4*xa + r4*xb))/(2.0_ki*r3),czip,(2*cone*r3*r4)/(xa + xb - Abs(xa - xb)))
      g106 = hp_glog((cone*(-r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone/(r3*r4),cone/(r3*r4))
      g107 = hp_glog((cone*(-r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone/(r3*r4),(2*cone*r3*r4)/(xa + xb - Abs(xa - xb)))
      g108 = hp_glog((cone*(-r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone*r3*r4,cone/(r3*r4))
      g109 = hp_glog((cone*(-r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone*r3*r4,(2*cone*r3*r4)/(xa + xb - Abs(xa - xb)))
      g110 = hp_glog((cone*(r5 + r4*xa + r4*xb))/(2.0_ki*r3),czip,cone/(r3*r4))
      g111 = hp_glog((cone*(r5 + r4*xa + r4*xb))/(2.0_ki*r3),czip,(2*cone*r3*r4)/(xa + xb - Abs(xa - xb)))
      g112 = hp_glog((cone*(r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone/(r3*r4),cone/(r3*r4))
      g113 = hp_glog((cone*(r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone/(r3*r4),(2*cone*r3*r4)/(xa + xb - Abs(xa - xb)))
      g114 = hp_glog((cone*(r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone*r3*r4,cone/(r3*r4))
      g115 = hp_glog((cone*(r5 + r4*xa + r4*xb))/(2.0_ki*r3),cone*r3*r4,(2*cone*r3*r4)/(xa + xb - Abs(xa - xb)))
      g116 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),czip,cone)
      g117 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),cone,cone)
      g118 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(1 + xb))/(1 - xa),cone)
      g119 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g120 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g121 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g122 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g123 = hp_glog((cone*(1 + xa*xb))/((-1 + xa)*(-1 + xb)),czip,cone)
      g124 = hp_glog((cone*(1 + xa*xb))/((-1 + xa)*(-1 + xb)),cone,cone)
      g125 = hp_glog((cone*(1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(1 + xb))/(1 - xa),cone)
      g126 = hp_glog((cone*(1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g127 = hp_glog((cone*(1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g128 = hp_glog((cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),czip,cone)
      g129 = hp_glog((cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone,cone)
      g130 = hp_glog((cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g131 = hp_glog((cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g132 = hp_glog((cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g133 = hp_glog((cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g134 = hp_glog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),czip,cone)
      g135 = hp_glog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone,cone)
      g136 = hp_glog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g137 = hp_glog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g138 = hp_glog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*(1 - r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g139 = hp_glog((cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),(cone*(1 + r1 + 2*xa*xb))/(2 - 2*xa - 2*xb + 2*xa*xb),cone)
      g140 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone,cone)
      g141 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g142 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g143 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g144 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g145 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone,cone)
      g146 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g147 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g148 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g149 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
