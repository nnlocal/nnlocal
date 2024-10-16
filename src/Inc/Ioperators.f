c--- The variables I10op1 and I10op2 provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (I10op1(a,b,c,i1,i2) and leg 2 (I10op2(a,b,c,i1,i2)
c--- In each case the parton labelling is
c--- everything forward
c---       incoming parton =   a
c---       hard parton     =   b
c---       emitted parton  =   c
c--- i1 and i1 denote the mom.frac. for leg1 and leg2 respectively
c--- i1,i2 = 1 denote the original born level mom.frac. while
c--- i1,i2 = 2 denote the real (rescaled) momentum fraction

      double precision 
     . I10op1(-1:1,-1:1,-1:1,2,2), I10op2(-1:1,-1:1,-1:1,2,2)
      common/I10op/ I10op1, I10op2

      

      
