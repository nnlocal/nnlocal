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
      complex(ki) hp_mylog
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
      l1 = hp_mylog(2*cone)
      l2 = hp_mylog(cone*(1 - xa))
      l3 = hp_mylog(cone*(2 - xa))
      l4 = hp_mylog(cone*xa)
      l5 = hp_mylog(cone*(1 + xa))
      l6 = hp_mylog(cone*(2 + xa))
      l7 = hp_mylog(cone*(1 - xb))
      l8 = hp_mylog(cone*(2 - xb))
      l9 = hp_mylog(cone*(1 - xa/xb))
      l10 = hp_mylog(cone*(2 - xa - xb))
      l11 = hp_mylog(cone*xb)
      l12 = hp_mylog(cone*(1 + xb))
      l13 = hp_mylog(cone*(2 + xb))
      l14 = hp_mylog(cone*(xa + xb))
      l15 = hp_mylog(cone*(2 + xa + xb))
      l16 = hp_mylog(cone*(1 - xb/xa))
      l17 = hp_mylog(cone*(1 - xa*xb))
      l18 = hp_mylog((cone*(r1*(1 - xb) - xa**2*(-1 + xb) - xa*(-1 + xb)**2))/(r1*(1 - xb) - xa**2*(-1 + xb)*xb))
      l19 = hp_mylog((cone*(r1*(1 - xb) - xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2))
      l20 = hp_mylog((cone*(r1*(1 - xb) + xa**2*(-1 + xb) + xa*(-1 + xb)**2))/(r1*(1 - xb) + xa**2*(-1 + xb)*xb))
      l21 = hp_mylog(-((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2)))
      l22 = hp_mylog((cone*(r2*(1 - xa) - (-1 + xa)**2*xb - (-1 + xa)*xb**2))/(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))
      l23 = hp_mylog((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb))
      l24 = hp_mylog((cone*(r2*(1 - xa) + (-1 + xa)**2*xb + (-1 + xa)*xb**2))/(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))
      l25 = hp_mylog(-((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)))
      g1 = hp_glog(-2*cone,-cone,cone*xa)
      g2 = hp_glog(-2*cone,-cone,cone*xb)
      g3 = hp_glog(-2*cone,czip,cone*xa)
      g4 = hp_glog(-2*cone,czip,cone*xb)
      g5 = hp_glog(-2*cone,2*cone,cone*xa)
      g6 = hp_glog(-2*cone,2*cone,cone*xb)
      g7 = hp_glog(-cone,-cone,cone*xa)
      g8 = hp_glog(-cone,-cone,cone*xb)
      g9 = hp_glog(-cone,czip,cone*xa)
      g10 = hp_glog(-cone,czip,cone*xb)
      g11 = hp_glog(-cone,cone,cone*xa)
      g12 = hp_glog(-cone,cone,cone*xb)
      g13 = hp_glog(-cone,2*cone,cone*xa)
      g14 = hp_glog(-cone,2*cone,cone*xb)
      g15 = hp_glog(-cone,cone*(2 - xa),cone*xb)
      g16 = hp_glog(-cone,cone/xa,cone*xb)
      g17 = hp_glog(-cone,-(cone*xa),cone*xb)
      g18 = hp_glog(-cone,cone*xa,-(cone*xa*xb))
      g19 = hp_glog(-cone,cone*(2 - xb),cone*xa)
      g20 = hp_glog(-cone,cone/xb,cone*xa)
      g21 = hp_glog(-cone,-(cone*xb),cone*xa)
      g22 = hp_glog(czip,-cone,cone*xa)
      g23 = hp_glog(czip,-cone,cone*xb)
      g24 = hp_glog(czip,czip,2*cone)
      g25 = hp_glog(czip,czip,cone*xa)
      g26 = hp_glog(czip,czip,cone*xb)
      g27 = hp_glog(czip,cone,cone*xa)
      g28 = hp_glog(czip,cone,cone*xb)
      g29 = hp_glog(czip,cone,cone*xa*xb)
      g30 = hp_glog(czip,2*cone,cone*xa)
      g31 = hp_glog(czip,2*cone,cone*xb)
      g32 = hp_glog(czip,cone*(2 - xa),cone*xb)
      g33 = hp_glog(czip,cone/xa,cone*xb)
      g34 = hp_glog(czip,-(cone*xa),cone*xb)
      g35 = hp_glog(czip,cone*(2 - xb),cone*xa)
      g36 = hp_glog(czip,cone/xb,cone*xa)
      g37 = hp_glog(czip,-(cone*xb),cone*xa)
      g38 = hp_glog(czip,(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g39 = hp_glog(czip,(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g40 = hp_glog(czip,(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g41 = hp_glog(czip,(cone*(xb - xa**2*xb))/(xa*(-1 + xb**2)),cone)
      g42 = hp_glog(czip,(cone*(xa - xa*xb**2))/((-1 + xa**2)*xb),cone)
      g43 = hp_glog(cone,-cone,cone*xa)
      g44 = hp_glog(cone,-cone,cone*xb)
      g45 = hp_glog(cone,czip,cone*xa)
      g46 = hp_glog(cone,czip,cone*xb)
      g47 = hp_glog(cone,cone,cone*xa)
      g48 = hp_glog(cone,cone,cone*xb)
      g49 = hp_glog(cone,cone,cone*xa*xb)
      g50 = hp_glog(cone,2*cone,cone*xa)
      g51 = hp_glog(cone,2*cone,cone*xb)
      g52 = hp_glog(cone,cone*(2 - xa),cone*xb)
      g53 = hp_glog(cone,cone/xa,cone*xb)
      g54 = hp_glog(cone,-(cone*xa),cone*xb)
      g55 = hp_glog(cone,-(cone*xa),cone*xa*xb)
      g56 = hp_glog(cone,cone*xa,cone*xa*xb)
      g57 = hp_glog(cone,cone*(2 - xb),cone*xa)
      g58 = hp_glog(cone,cone/xb,cone*xa)
      g59 = hp_glog(cone,-(cone*xb),cone*xa)
      g60 = hp_glog(2*cone,-cone,cone*xa)
      g61 = hp_glog(2*cone,-cone,cone*xb)
      g62 = hp_glog(2*cone,czip,cone*xa)
      g63 = hp_glog(2*cone,czip,cone*xb)
      g64 = hp_glog(2*cone,2*cone,cone*xa)
      g65 = hp_glog(2*cone,2*cone,cone*xb)
      g66 = hp_glog(cone*(-2 - xa),-cone,cone*xb)
      g67 = hp_glog(cone*(-2 - xa),cone*(2 - xa),cone*xb)
      g68 = hp_glog(cone*(-2 - xa),cone/xa,cone*xb)
      g69 = hp_glog(cone*(-2 - xa),-(cone*xa),cone*xb)
      g70 = hp_glog(cone*(2 - xa),-cone,cone*xb)
      g71 = hp_glog(cone*(2 - xa),czip,cone*xb)
      g72 = hp_glog(cone*(2 - xa),cone*(2 - xa),cone*xb)
      g73 = hp_glog(cone*(2 - xa),cone/xa,cone*xb)
      g74 = hp_glog(cone*(2 - xa),-(cone*xa),cone*xb)
      g75 = hp_glog(cone/xa,cone*(2 - xa),cone*xb)
      g76 = hp_glog(cone/xa,cone/xa,cone*xb)
      g77 = hp_glog(-(cone*xa),-cone,cone*xb)
      g78 = hp_glog(-(cone*xa),czip,cone*xb)
      g79 = hp_glog(-(cone*xa),cone,cone*xb)
      g80 = hp_glog(-(cone*xa),cone,cone*xa*xb)
      g81 = hp_glog(-(cone*xa),-(cone*xa),cone*xb)
      g82 = hp_glog(cone*xa,czip,cone*xb)
      g83 = hp_glog(cone*xa,cone,cone*xa*xb)
      g84 = hp_glog(cone*xa,cone*(2 - xa),cone*xb)
      g85 = hp_glog(cone*xa,cone/xa,cone*xb)
      g86 = hp_glog(cone*xa,-(cone*xa),cone*xb)
      g87 = hp_glog(cone*(-2 - xb),-cone,cone*xa)
      g88 = hp_glog(cone*(-2 - xb),cone*(2 - xb),cone*xa)
      g89 = hp_glog(cone*(-2 - xb),cone/xb,cone*xa)
      g90 = hp_glog(cone*(-2 - xb),-(cone*xb),cone*xa)
      g91 = hp_glog(cone*(2 - xb),-cone,cone*xa)
      g92 = hp_glog(cone*(2 - xb),czip,cone*xa)
      g93 = hp_glog(cone*(2 - xb),cone*(2 - xb),cone*xa)
      g94 = hp_glog(cone*(2 - xb),cone/xb,cone*xa)
      g95 = hp_glog(cone*(2 - xb),-(cone*xb),cone*xa)
      g96 = hp_glog(cone/xb,cone*(2 - xb),cone*xa)
      g97 = hp_glog(cone/xb,cone/xb,cone*xa)
      g98 = hp_glog(-(cone*xb),-cone,cone*xa)
      g99 = hp_glog(-(cone*xb),czip,cone*xa)
      g100 = hp_glog(-(cone*xb),cone,cone*xa)
      g101 = hp_glog(-(cone*xb),cone,cone*xa*xb)
      g102 = hp_glog(-(cone*xb),-(cone*xb),cone*xa)
      g103 = hp_glog(cone*xb,czip,cone*xa)
      g104 = hp_glog(cone*xb,cone,cone*xa*xb)
      g105 = hp_glog(cone*xb,cone*(2 - xb),cone*xa)
      g106 = hp_glog(cone*xb,cone/xb,cone*xa)
      g107 = hp_glog(cone*xb,-(cone*xb),cone*xa)
      g108 = hp_glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),czip,cone)
      g109 = hp_glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone,cone)
      g110 = hp_glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g111 = hp_glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g112 = hp_glog((cone*(1 + xa)*xb)/(xa*(-1 + xb)),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g113 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),czip,cone)
      g114 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),cone,cone)
      g115 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g116 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g117 = hp_glog((cone*xa*(1 + xb))/((-1 + xa)*xb),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g118 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g119 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g120 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g121 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g122 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g123 = hp_glog((cone*(-1 + xa*xb))/((-1 + xa)*(-1 + xb)),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g124 = hp_glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone,cone)
      g125 = hp_glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g126 = hp_glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g127 = hp_glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g128 = hp_glog((cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g129 = hp_glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone,cone)
      g130 = hp_glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*xb)/(xa*(-1 + xb)),cone)
      g131 = hp_glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g132 = hp_glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(-(r1*(1 - xb)) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g133 = hp_glog((cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),(cone*(r1*(1 - xb) + xa**2*(-1 + xb)*xb))/((-1 + xa)*xa*(-1 + xb)**2),cone)
      g134 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone,cone)
      g135 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g136 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g137 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g138 = hp_glog(-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
      g139 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone,cone)
      g140 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*(1 + xa)*(1 + xb))/((-1 + xa)*(-1 + xb)),cone)
      g141 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*xa*(1 + xb))/((-1 + xa)*xb),cone)
      g142 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),-((cone*(r2*(1 - xa) + xa*xb**2 - xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb)),cone)
      g143 = hp_glog((cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),(cone*(r2*(1 - xa) - xa*xb**2 + xa**2*xb**2))/((-1 + xa)**2*(-1 + xb)*xb),cone)
