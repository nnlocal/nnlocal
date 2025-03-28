c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting 
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      double precision 
     . T1(-1:1,-1:1,-1:1,4),T2(-1:1,-1:1,-1:1,4),
     . Tp1(3,-1:1,-1:1,-1:1,4,0:2),Tp2(3,-1:1,-1:1,-1:1,4,0:2)
      common/RP_new/T1,T2,Tp1,Tp2
