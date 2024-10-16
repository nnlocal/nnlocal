      FUNCTION HP_WGPLG(N,P,X)
c..
c..   Nielson's Generalized Polylogarithm
c..
      implicit none
      include 'hp_types.h'
      complex(ki) HP_WGPLG
      integer i,k,l,m,n,n1
      INTEGER P,P1,NC(10),INDEX(31)
      real(ki) FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4)
      real(ki) A(0:30,10)
      real(ki) X,X1,H,ALFA,R,Q,C1,C2,B0,B1,B2,ZERO,HALF

      complex(ki) V(0:5),SK,SM

      DATA FCT /1.0_ki,1.0_ki,2.0_ki,6.0_ki,24.0_ki/
      DATA SGN /1.0_ki,-1.0_ki,1.0_ki,-1.0_ki,1.0_ki/
      DATA ZERO /0.0_ki/, HALF /0.5_ki/
      DATA C1 /1.33333 33333 333_ki/, C2 /0.33333 33333 3333_ki/

      DATA S1(1,1) /1.64493 40668 482_ki/
      DATA S1(1,2) /1.20205 69031 596_ki/
      DATA S1(1,3) /1.08232 32337 111_ki/
      DATA S1(1,4) /1.03692 77551 434_ki/
      DATA S1(2,1) /1.20205 69031 596_ki/
      DATA S1(2,2) /2.70580 80842 778E-1_ki/
      DATA S1(2,3) /9.65511 59989 444E-2_ki/
      DATA S1(3,1) /1.08232 32337 111_ki/
      DATA S1(3,2) /9.65511 59989 444E-2_ki/
      DATA S1(4,1) /1.03692 77551 434_ki/

      DATA C(1,1) / 1.64493 40668 482_ki/
      DATA C(1,2) / 1.20205 69031 596_ki/
      DATA C(1,3) / 1.08232 32337 111_ki/
      DATA C(1,4) / 1.03692 77551 434_ki/
      DATA C(2,1) / 0.00000 00000 000_ki/
      DATA C(2,2) /-1.89406 56589 945_ki/
      DATA C(2,3) /-3.01423 21054 407_ki/
      DATA C(3,1) / 1.89406 56589 945_ki/
      DATA C(3,2) / 3.01423 21054 407_ki/
      DATA C(4,1) / 0.00000 00000 000_ki/

      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/

      DATA NC /24,26,28,30,22,24,26,19,22,17/

      DATA A( 0,1) / .96753 21504 3498_ki/
      DATA A( 1,1) / .16607 30329 2785_ki/
      DATA A( 2,1) / .02487 93229 2423_ki/
      DATA A( 3,1) / .00468 63619 5945_ki/
      DATA A( 4,1) / .00100 16274 9616_ki/
      DATA A( 5,1) / .00023 20021 9609_ki/
      DATA A( 6,1) / .00005 68178 2272_ki/
      DATA A( 7,1) / .00001 44963 0056_ki/
      DATA A( 8,1) / .00000 38163 2946_ki/
      DATA A( 9,1) / .00000 10299 0426_ki/
      DATA A(10,1) / .00000 02835 7538_ki/
      DATA A(11,1) / .00000 00793 8705_ki/
      DATA A(12,1) / .00000 00225 3670_ki/
      DATA A(13,1) / .00000 00064 7434_ki/
      DATA A(14,1) / .00000 00018 7912_ki/
      DATA A(15,1) / .00000 00005 5029_ki/
      DATA A(16,1) / .00000 00001 6242_ki/
      DATA A(17,1) / .00000 00000 4827_ki/
      DATA A(18,1) / .00000 00000 1444_ki/
      DATA A(19,1) / .00000 00000 0434_ki/
      DATA A(20,1) / .00000 00000 0131_ki/
      DATA A(21,1) / .00000 00000 0040_ki/
      DATA A(22,1) / .00000 00000 0012_ki/
      DATA A(23,1) / .00000 00000 0004_ki/
      DATA A(24,1) / .00000 00000 0001_ki/

      DATA A( 0,2) / .95180 88912 7832_ki/
      DATA A( 1,2) / .43131 13184 6532_ki/
      DATA A( 2,2) / .10002 25071 4905_ki/
      DATA A( 3,2) / .02442 41559 5220_ki/
      DATA A( 4,2) / .00622 51246 3724_ki/
      DATA A( 5,2) / .00164 07883 1235_ki/
      DATA A( 6,2) / .00044 40792 0265_ki/
      DATA A( 7,2) / .00012 27749 4168_ki/
      DATA A( 8,2) / .00003 45398 1284_ki/
      DATA A( 9,2) / .00000 98586 9565_ki/
      DATA A(10,2) / .00000 28485 6995_ki/
      DATA A(11,2) / .00000 08317 0847_ki/
      DATA A(12,2) / .00000 02450 3950_ki/
      DATA A(13,2) / .00000 00727 6496_ki/
      DATA A(14,2) / .00000 00217 5802_ki/
      DATA A(15,2) / .00000 00065 4616_ki/
      DATA A(16,2) / .00000 00019 8033_ki/
      DATA A(17,2) / .00000 00006 0204_ki/
      DATA A(18,2) / .00000 00001 8385_ki/
      DATA A(19,2) / .00000 00000 5637_ki/
      DATA A(20,2) / .00000 00000 1735_ki/
      DATA A(21,2) / .00000 00000 0536_ki/
      DATA A(22,2) / .00000 00000 0166_ki/
      DATA A(23,2) / .00000 00000 0052_ki/
      DATA A(24,2) / .00000 00000 0016_ki/
      DATA A(25,2) / .00000 00000 0005_ki/
      DATA A(26,2) / .00000 00000 0002_ki/

      DATA A( 0,3) / .98161 02799 1365_ki/
      DATA A( 1,3) / .72926 80632 0726_ki/
      DATA A( 2,3) / .22774 71490 9321_ki/
      DATA A( 3,3) / .06809 08329 6197_ki/
      DATA A( 4,3) / .02013 70118 3064_ki/
      DATA A( 5,3) / .00595 47848 0197_ki/
      DATA A( 6,3) / .00176 76901 3959_ki/
      DATA A( 7,3) / .00052 74821 8502_ki/
      DATA A( 8,3) / .00015 82746 1460_ki/
      DATA A( 9,3) / .00004 77492 2076_ki/
      DATA A(10,3) / .00001 44792 0408_ki/
      DATA A(11,3) / .00000 44115 4886_ki/
      DATA A(12,3) / .00000 13500 3870_ki/
      DATA A(13,3) / .00000 04148 1779_ki/
      DATA A(14,3) / .00000 01279 3307_ki/
      DATA A(15,3) / .00000 00395 9070_ki/
      DATA A(16,3) / .00000 00122 9055_ki/
      DATA A(17,3) / .00000 00038 2658_ki/
      DATA A(18,3) / .00000 00011 9459_ki/
      DATA A(19,3) / .00000 00003 7386_ki/
      DATA A(20,3) / .00000 00001 1727_ki/
      DATA A(21,3) / .00000 00000 3687_ki/
      DATA A(22,3) / .00000 00000 1161_ki/
      DATA A(23,3) / .00000 00000 0366_ki/
      DATA A(24,3) / .00000 00000 0116_ki/
      DATA A(25,3) / .00000 00000 0037_ki/
      DATA A(26,3) / .00000 00000 0012_ki/
      DATA A(27,3) / .00000 00000 0004_ki/
      DATA A(28,3) / .00000 00000 0001_ki/

      DATA A( 0,4) /1.06405 21184 614 _ki/
      DATA A( 1,4) /1.06917 20744 981 _ki/
      DATA A( 2,4) / .41527 19325 1768_ki/
      DATA A( 3,4) / .14610 33293 6222_ki/
      DATA A( 4,4) / .04904 73264 8784_ki/
      DATA A( 5,4) / .01606 34086 0396_ki/
      DATA A( 6,4) / .00518 88935 0790_ki/
      DATA A( 7,4) / .00166 29871 7324_ki/
      DATA A( 8,4) / .00053 05827 9969_ki/
      DATA A( 9,4) / .00016 88702 9251_ki/
      DATA A(10,4) / .00005 36832 8059_ki/
      DATA A(11,4) / .00001 70592 3313_ki/
      DATA A(12,4) / .00000 54217 4374_ki/
      DATA A(13,4) / .00000 17239 4082_ki/
      DATA A(14,4) / .00000 05485 3275_ki/
      DATA A(15,4) / .00000 01746 7795_ki/
      DATA A(16,4) / .00000 00556 7550_ki/
      DATA A(17,4) / .00000 00177 6234_ki/
      DATA A(18,4) / .00000 00056 7224_ki/
      DATA A(19,4) / .00000 00018 1313_ki/
      DATA A(20,4) / .00000 00005 8012_ki/
      DATA A(21,4) / .00000 00001 8579_ki/
      DATA A(22,4) / .00000 00000 5955_ki/
      DATA A(23,4) / .00000 00000 1911_ki/
      DATA A(24,4) / .00000 00000 0614_ki/
      DATA A(25,4) / .00000 00000 0197_ki/
      DATA A(26,4) / .00000 00000 0063_ki/
      DATA A(27,4) / .00000 00000 0020_ki/
      DATA A(28,4) / .00000 00000 0007_ki/
      DATA A(29,4) / .00000 00000 0002_ki/
      DATA A(30,4) / .00000 00000 0001_ki/

      DATA A( 0,5) / .97920 86066 9175_ki/
      DATA A( 1,5) / .08518 81314 8683_ki/
      DATA A( 2,5) / .00855 98522 2013_ki/
      DATA A( 3,5) / .00121 17721 4413_ki/
      DATA A( 4,5) / .00020 72276 8531_ki/
      DATA A( 5,5) / .00003 99695 8691_ki/
      DATA A( 6,5) / .00000 83806 4065_ki/
      DATA A( 7,5) / .00000 18684 8945_ki/
      DATA A( 8,5) / .00000 04366 6087_ki/
      DATA A( 9,5) / .00000 01059 1733_ki/
      DATA A(10,5) / .00000 00264 7892_ki/
      DATA A(11,5) / .00000 00067 8700_ki/
      DATA A(12,5) / .00000 00017 7654_ki/
      DATA A(13,5) / .00000 00004 7342_ki/
      DATA A(14,5) / .00000 00001 2812_ki/
      DATA A(15,5) / .00000 00000 3514_ki/
      DATA A(16,5) / .00000 00000 0975_ki/
      DATA A(17,5) / .00000 00000 0274_ki/
      DATA A(18,5) / .00000 00000 0077_ki/
      DATA A(19,5) / .00000 00000 0022_ki/
      DATA A(20,5) / .00000 00000 0006_ki/
      DATA A(21,5) / .00000 00000 0002_ki/
      DATA A(22,5) / .00000 00000 0001_ki/

      DATA A( 0,6) / .95021 85196 3952_ki/
      DATA A( 1,6) / .29052 52916 1433_ki/
      DATA A( 2,6) / .05081 77406 1716_ki/
      DATA A( 3,6) / .00995 54376 7280_ki/
      DATA A( 4,6) / .00211 73389 5031_ki/
      DATA A( 5,6) / .00047 85947 0550_ki/
      DATA A( 6,6) / .00011 33432 1308_ki/
      DATA A( 7,6) / .00002 78473 3104_ki/
      DATA A( 8,6) / .00000 70478 8108_ki/
      DATA A( 9,6) / .00000 18278 8740_ki/
      DATA A(10,6) / .00000 04838 7492_ki/
      DATA A(11,6) / .00000 01303 3842_ki/
      DATA A(12,6) / .00000 00356 3769_ki/
      DATA A(13,6) / .00000 00098 7174_ki/
      DATA A(14,6) / .00000 00027 6586_ki/
      DATA A(15,6) / .00000 00007 8279_ki/
      DATA A(16,6) / .00000 00002 2354_ki/
      DATA A(17,6) / .00000 00000 6435_ki/
      DATA A(18,6) / .00000 00000 1866_ki/
      DATA A(19,6) / .00000 00000 0545_ki/
      DATA A(20,6) / .00000 00000 0160_ki/
      DATA A(21,6) / .00000 00000 0047_ki/
      DATA A(22,6) / .00000 00000 0014_ki/
      DATA A(23,6) / .00000 00000 0004_ki/
      DATA A(24,6) / .00000 00000 0001_ki/

      DATA A( 0,7) / .95064 03218 6777_ki/
      DATA A( 1,7) / .54138 28546 5171_ki/
      DATA A( 2,7) / .13649 97959 0321_ki/
      DATA A( 3,7) / .03417 94232 8207_ki/
      DATA A( 4,7) / .00869 02788 3583_ki/
      DATA A( 5,7) / .00225 28408 4155_ki/
      DATA A( 6,7) / .00059 51608 9806_ki/
      DATA A( 7,7) / .00015 99561 7766_ki/
      DATA A( 8,7) / .00004 36521 3096_ki/
      DATA A( 9,7) / .00001 20747 4688_ki/
      DATA A(10,7) / .00000 33801 8176_ki/
      DATA A(11,7) / .00000 09563 2476_ki/
      DATA A(12,7) / .00000 02731 3129_ki/
      DATA A(13,7) / .00000 00786 6968_ki/
      DATA A(14,7) / .00000 00228 3195_ki/
      DATA A(15,7) / .00000 00066 7205_ki/
      DATA A(16,7) / .00000 00019 6191_ki/
      DATA A(17,7) / .00000 00005 8018_ki/
      DATA A(18,7) / .00000 00001 7246_ki/
      DATA A(19,7) / .00000 00000 5151_ki/
      DATA A(20,7) / .00000 00000 1545_ki/
      DATA A(21,7) / .00000 00000 0465_ki/
      DATA A(22,7) / .00000 00000 0141_ki/
      DATA A(23,7) / .00000 00000 0043_ki/
      DATA A(24,7) / .00000 00000 0013_ki/
      DATA A(25,7) / .00000 00000 0004_ki/
      DATA A(26,7) / .00000 00000 0001_ki/

      DATA A( 0,8) / .98800 01167 2229_ki/
      DATA A( 1,8) / .04364 06760 9601_ki/
      DATA A( 2,8) / .00295 09117 8278_ki/
      DATA A( 3,8) / .00031 47780 9720_ki/
      DATA A( 4,8) / .00004 31484 6029_ki/
      DATA A( 5,8) / .00000 69381 8230_ki/
      DATA A( 6,8) / .00000 12464 0350_ki/
      DATA A( 7,8) / .00000 02429 3628_ki/
      DATA A( 8,8) / .00000 00504 0827_ki/
      DATA A( 9,8) / .00000 00109 9075_ki/
      DATA A(10,8) / .00000 00024 9467_ki/
      DATA A(11,8) / .00000 00005 8540_ki/
      DATA A(12,8) / .00000 00001 4127_ki/
      DATA A(13,8) / .00000 00000 3492_ki/
      DATA A(14,8) / .00000 00000 0881_ki/
      DATA A(15,8) / .00000 00000 0226_ki/
      DATA A(16,8) / .00000 00000 0059_ki/
      DATA A(17,8) / .00000 00000 0016_ki/
      DATA A(18,8) / .00000 00000 0004_ki/
      DATA A(19,8) / .00000 00000 0001_ki/

      DATA A( 0,9) / .95768 50654 6350_ki/
      DATA A( 1,9) / .19725 24967 9534_ki/
      DATA A( 2,9) / .02603 37031 3918_ki/
      DATA A( 3,9) / .00409 38216 8261_ki/
      DATA A( 4,9) / .00072 68170 7110_ki/
      DATA A( 5,9) / .00014 09187 9261_ki/
      DATA A( 6,9) / .00002 92045 8914_ki/
      DATA A( 7,9) / .00000 63763 1144_ki/
      DATA A( 8,9) / .00000 14516 7850_ki/
      DATA A( 9,9) / .00000 03420 5281_ki/
      DATA A(10,9) / .00000 00829 4302_ki/
      DATA A(11,9) / .00000 00206 0784_ki/
      DATA A(12,9) / .00000 00052 2823_ki/
      DATA A(13,9) / .00000 00013 5066_ki/
      DATA A(14,9) / .00000 00003 5451_ki/
      DATA A(15,9) / .00000 00000 9436_ki/
      DATA A(16,9) / .00000 00000 2543_ki/
      DATA A(17,9) / .00000 00000 0693_ki/
      DATA A(18,9) / .00000 00000 0191_ki/
      DATA A(19,9) / .00000 00000 0053_ki/
      DATA A(20,9) / .00000 00000 0015_ki/
      DATA A(21,9) / .00000 00000 0004_ki/
      DATA A(22,9) / .00000 00000 0001_ki/

      DATA A( 0,10) / .99343 65167 1347_ki/
      DATA A( 1,10) / .02225 77012 6826_ki/
      DATA A( 2,10) / .00101 47557 4703_ki/
      DATA A( 3,10) / .00008 17515 6250_ki/
      DATA A( 4,10) / .00000 89997 3547_ki/
      DATA A( 5,10) / .00000 12082 3987_ki/
      DATA A( 6,10) / .00000 01861 6913_ki/
      DATA A( 7,10) / .00000 00317 4723_ki/
      DATA A( 8,10) / .00000 00058 5215_ki/
      DATA A( 9,10) / .00000 00011 4739_ki/
      DATA A(10,10) / .00000 00002 3652_ki/
      DATA A(11,10) / .00000 00000 5082_ki/
      DATA A(12,10) / .00000 00000 1131_ki/
      DATA A(13,10) / .00000 00000 0259_ki/
      DATA A(14,10) / .00000 00000 0061_ki/
      DATA A(15,10) / .00000 00000 0015_ki/
      DATA A(16,10) / .00000 00000 0004_ki/
      DATA A(17,10) / .00000 00000 0001_ki/

      IF(N .LT. 1 .OR. N .GT. 4 .OR. P .LT. 1 .OR. P .GT. 4 .OR.
     1   N+P .GT. 5) THEN
       HP_WGPLG=ZERO
       PRINT 1000, N,P
       RETURN
      END IF
      IF(X .EQ. SGN(0)) THEN
       HP_WGPLG=S1(N,P)
       RETURN
      END IF

      IF(X .GT. FCT(2) .OR. X .LT. SGN(1)) THEN
       X1=SGN(0)/X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       V(1)=LOG(cmplx(-X,ZERO,kind=ki))
       DO 33 L = 2,N+P
   33  V(L)=V(1)*V(L-1)/L
       SK=ZERO
       DO 34 K = 0,P-1
       P1=P-K
       R=X1**P1/(FCT(P1)*FCT(N-1))
       SM=ZERO
       DO 35 M = 0,K
       N1=N+K-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 31 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   31  B1=B0
       Q=(FCT(N1-1)/FCT(K-M))*(B0-H*B2)*R/P1**N1
   35  SM=SM+V(M)*Q
   34  SK=SK+SGN(K)*SM
       SM=ZERO
       DO 36 M = 0,N-1
   36  SM=SM+V(M)*C(N-M,P)
       HP_WGPLG=SGN(N)*SK+SGN(P)*(SM+V(N+P))
       RETURN
      END IF

      IF(X .GT. HALF) THEN
       X1=SGN(0)-X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       U(0)=SGN(0)
       V(1)=LOG(cmplx(X1,ZERO,kind=ki))
       U(1)=LOG(X)
       DO 23 L = 2,P
   23  V(L)=V(1)*V(L-1)/L
       DO 26 L = 2,N
   26  U(L)=U(1)*U(L-1)/L
       SK=ZERO
       DO 24 K = 0,N-1
       P1=N-K
       R=X1**P1/FCT(P1)
       SM=ZERO
       DO 25 M = 0,P-1
       N1=P-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 12 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   12  B1=B0
       Q=SGN(M)*(B0-H*B2)*R/P1**N1
   25  SM=SM+V(M)*Q
   24  SK=SK+U(K)*(S1(P1,P)-SM)
       HP_WGPLG=SK+SGN(P)*U(N)*V(P)
       RETURN
      END IF

      L=INDEX(10*N+P-10)
      H=C1*X+C2
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 11 I = NC(L),0,-1
      B0=A(I,L)+ALFA*B1-B2
      B2=B1
   11 B1=B0
      HP_WGPLG=(B0-H*B2)*X**P/(FCT(P)*P**N)
      RETURN
 1000 FORMAT(/' ***** CERN SUBROUTINE HP_WGPLG ... ILLEGAL VALUES',
     1        '   N = ',I3,'   P = ',I3)
      END

C-}}}


