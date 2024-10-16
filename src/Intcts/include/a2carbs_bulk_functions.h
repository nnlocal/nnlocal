!--- include file for double soft ct's
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
      complex(ki) r1
      complex(ki) r2
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
      complex(ki) mylog
      r1 = Sqrt(cone*xa*xb*(-1 + xa**2 + xa*xb))
      r2 = Sqrt(cone*xa*xb*(-1 + xa*xb + xb**2))
      z2 = 1/(-1 + xb*(xa + xb))
      z3 = 1/(xa + xb)
      z5 = 1/(1 + xb)
      z6 = 1/(-1 + xb)
      z8 = 1/(xa - xb)
      z9 = 1/(-1 + xa*(xa + xb))
      z10 = 1/(1 + xa)
      z4 = 1/(-1 + xa*xb)
      z11 = 1/(-1 + xa)
      z1 = 1/xa
      z7 = 1/xb
      l1 = mylog(2*cone)
      l2 = mylog(cone*(1 - xa))
      l3 = mylog(cone*(2 - xa))
      l4 = mylog(cone*xa)
      l5 = mylog(cone*(1 + xa))
      l6 = mylog(cone*(2 + xa))
      l7 = mylog(cone*(1 - xb))
      l8 = mylog(cone*(2 - xb))
      l9 = mylog(cone*(1 - xa/xb))
      l10 = mylog(cone*(2 - xa - xb))
      l11 = mylog(cone*xb)
      l12 = mylog(cone*(1 + xb))
      l13 = mylog(cone*(2 + xb))
      l14 = mylog(cone*(xa + xb))
      l15 = mylog(cone*(2 + xa + xb))
      l16 = mylog(cone*(1 - xb/xa))
      l17 = mylog(cone*(1 - xa*xb))
      l18 = mylog((cone*(r1*(1 - xb) - xa**2*(-1 + xb) - xa*(-1 + xb)**2))/(r1*(1 - xb) - xa**2*(-1 + xb)*xb))
      l19 = mylog((cone*(r1*(1 - xb) - xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2))
      l20 = mylog((cone*(r1*(1 - xb) + xa**2*(-1 + xb) + xa*(-1 + xb)**2))/(r1*(1 - xb) + xa**2*(-1 + xb)*xb))
      l21 = mylog(-((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2)))
      l22 = mylog((cone*(r2*(1 - xa) - (-1 + xa)**2*xb - (-1 + xa)*xb**2))/(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))
      l23 = mylog((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb))
      l24 = mylog((cone*(r2*(1 - xa) + (-1 + xa)**2*xb + (-1 + xa)*xb**2))/(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))
      l25 = mylog(-((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)))
      g1 = glog(-2*cone,-cone,cone*xa)
      g2 = glog(-2*cone,-cone,cone*xb)
      g3 = glog(-2*cone,czip,cone*xa)
      g4 = glog(-2*cone,czip,cone*xb)
      g5 = glog(-2*cone,2*cone,cone*xa)
      g6 = glog(-2*cone,2*cone,cone*xb)
      g7 = glog(-cone,-cone,cone*xa)
      g8 = glog(-cone,-cone,cone*xb)
      g9 = glog(-cone,czip,cone*xa)
      g10 = glog(-cone,czip,cone*xb)
      g11 = glog(-cone,cone,cone*xa)
      g12 = glog(-cone,cone,cone*xb)
      g13 = glog(-cone,2*cone,cone*xa)
      g14 = glog(-cone,2*cone,cone*xb)
      g15 = glog(-cone,cone*(2 - xa),cone*xb)
      g16 = glog(-cone,cone/xa,cone*xb)
      g17 = glog(-cone,-(cone*xa),cone*xb)
      g18 = glog(-cone,cone*xa,-(cone*xa*xb))
      g19 = glog(-cone,cone*(2 - xb),cone*xa)
      g20 = glog(-cone,cone/xb,cone*xa)
      g21 = glog(-cone,-(cone*xb),cone*xa)
      g22 = glog(czip,-cone,cone*xa)
      g23 = glog(czip,-cone,cone*xb)
      g24 = glog(czip,czip,2*cone)
      g25 = glog(czip,czip,cone*xa)
      g26 = glog(czip,czip,cone*xb)
      g27 = glog(czip,cone,cone*xa)
      g28 = glog(czip,cone,cone*xb)
      g29 = glog(czip,cone,cone*xa*xb)
      g30 = glog(czip,2*cone,cone*xa)
      g31 = glog(czip,2*cone,cone*xb)
      g32 = glog(czip,cone*(2 - xa),cone*xb)
      g33 = glog(czip,cone/xa,cone*xb)
      g34 = glog(czip,-(cone*xa),cone*xb)
      g35 = glog(czip,cone*(2 - xb),cone*xa)
      g36 = glog(czip,cone/xb,cone*xa)
      g37 = glog(czip,-(cone*xb),cone*xa)
      g38 = glog(czip,(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g39 = glog(czip,(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g40 = glog(czip,(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g41 = glog(czip,(cone*(xb - xa**2*xb))/(xa*(-1 + xb**2)),cone)
      g42 = glog(czip,(cone*(xa - xa*xb**2))/((-1 + xa**2)*xb),cone)
      g43 = glog(cone,-cone,cone*xa)
      g44 = glog(cone,-cone,cone*xb)
      g45 = glog(cone,czip,cone*xa)
      g46 = glog(cone,czip,cone*xb)
      g47 = glog(cone,cone,cone*xa)
      g48 = glog(cone,cone,cone*xb)
      g49 = glog(cone,cone,cone*xa*xb)
      g50 = glog(cone,2*cone,cone*xa)
      g51 = glog(cone,2*cone,cone*xb)
      g52 = glog(cone,cone*(2 - xa),cone*xb)
      g53 = glog(cone,cone/xa,cone*xb)
      g54 = glog(cone,-(cone*xa),cone*xb)
      g55 = glog(cone,-(cone*xa),cone*xa*xb)
      g56 = glog(cone,cone*xa,cone*xa*xb)
      g57 = glog(cone,cone*(2 - xb),cone*xa)
      g58 = glog(cone,cone/xb,cone*xa)
      g59 = glog(cone,-(cone*xb),cone*xa)
      g60 = glog(2*cone,-cone,cone*xa)
      g61 = glog(2*cone,-cone,cone*xb)
      g62 = glog(2*cone,czip,cone*xa)
      g63 = glog(2*cone,czip,cone*xb)
      g64 = glog(2*cone,2*cone,cone*xa)
      g65 = glog(2*cone,2*cone,cone*xb)
      g66 = glog(cone*(-2 - xa),-cone,cone*xb)
      g67 = glog(cone*(-2 - xa),cone*(2 - xa),cone*xb)
      g68 = glog(cone*(-2 - xa),cone/xa,cone*xb)
      g69 = glog(cone*(-2 - xa),-(cone*xa),cone*xb)
      g70 = glog(cone*(2 - xa),-cone,cone*xb)
      g71 = glog(cone*(2 - xa),czip,cone*xb)
      g72 = glog(cone*(2 - xa),cone*(2 - xa),cone*xb)
      g73 = glog(cone*(2 - xa),cone/xa,cone*xb)
      g74 = glog(cone*(2 - xa),-(cone*xa),cone*xb)
      g75 = glog(cone/xa,cone*(2 - xa),cone*xb)
      g76 = glog(cone/xa,cone/xa,cone*xb)
      g77 = glog(-(cone*xa),-cone,cone*xb)
      g78 = glog(-(cone*xa),czip,cone*xb)
      g79 = glog(-(cone*xa),cone,cone*xb)
      g80 = glog(-(cone*xa),cone,cone*xa*xb)
      g81 = glog(-(cone*xa),-(cone*xa),cone*xb)
      g82 = glog(cone*xa,czip,cone*xb)
      g83 = glog(cone*xa,cone,cone*xa*xb)
      g84 = glog(cone*xa,cone*(2 - xa),cone*xb)
      g85 = glog(cone*xa,cone/xa,cone*xb)
      g86 = glog(cone*xa,-(cone*xa),cone*xb)
      g87 = glog(cone*(-2 - xb),-cone,cone*xa)
      g88 = glog(cone*(-2 - xb),cone*(2 - xb),cone*xa)
      g89 = glog(cone*(-2 - xb),cone/xb,cone*xa)
      g90 = glog(cone*(-2 - xb),-(cone*xb),cone*xa)
      g91 = glog(cone*(2 - xb),-cone,cone*xa)
      g92 = glog(cone*(2 - xb),czip,cone*xa)
      g93 = glog(cone*(2 - xb),cone*(2 - xb),cone*xa)
      g94 = glog(cone*(2 - xb),cone/xb,cone*xa)
      g95 = glog(cone*(2 - xb),-(cone*xb),cone*xa)
      g96 = glog(cone/xb,cone*(2 - xb),cone*xa)
      g97 = glog(cone/xb,cone/xb,cone*xa)
      g98 = glog(-(cone*xb),-cone,cone*xa)
      g99 = glog(-(cone*xb),czip,cone*xa)
      g100 = glog(-(cone*xb),cone,cone*xa)
      g101 = glog(-(cone*xb),cone,cone*xa*xb)
      g102 = glog(-(cone*xb),-(cone*xb),cone*xa)
      g103 = glog(cone*xb,czip,cone*xa)
      g104 = glog(cone*xb,cone,cone*xa*xb)
      g105 = glog(cone*xb,cone*(2 - xb),cone*xa)
      g106 = glog(cone*xb,cone/xb,cone*xa)
      g107 = glog(cone*xb,-(cone*xb),cone*xa)
      g108 = glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),czip,cone)
      g109 = glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone,cone)
      g110 = glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g111 = glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g112 = glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g113 = glog((cone*xa*(1 + xb))/((-1 + xa)*xb),czip,cone)
      g114 = glog((cone*xa*(1 + xb))/((-1 + xa)*xb),cone,cone)
      g115 = glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g116 = glog((cone*xa*(1 + xb))/((-1 + xa)*xb),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g117 = glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g118 = glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g119 = glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g120 = glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g121 = glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g122 = glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g123 = glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g124 = glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone,cone)
      g125 = glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g126 = glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g127 = glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g128 = glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g129 = glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone,cone)
      g130 = glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g131 = glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g132 = glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g133 = glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g134 = glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone,cone)
      g135 = glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g136 = glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g137 = glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g138 = glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g139 = glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone,cone)
      g140 = glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g141 = glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g142 = glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g143 = glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
