************************************************************************

      SUBROUTINE RAMBO(N,ET,XM,P,WT)
C------------------------------------------------------
C
C              RAMBO
C
C    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
C
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
C    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
C
C    N    = NUMBER OF PARTICLES (>1, IN THIS VERSION <10)
C    ET = TOTAL CENTRE-OF-MASS ENERGY
C    XM = PARTICLE MASSES ( DIM=9 )
C    P    = PARTICLE MOMENTA ( DIM=(4,9) )
C    WT = WEIGHT OF THE EVENT
C
C------------------------------------------------------

C---Removed implicit declarations and N<10 limitation
C---Gabor Somogyi March 2014
      implicit none
      include 'constants.f'

      integer n
      double precision et,wt
      
      double precision xm(mxpart),p(4,mxpart),q(4,mxpart),z(mxpart)
     .     ,r(4), b(3),p2(mxpart),xm2(mxpart),e(mxpart),v(mxpart)      
      SAVE Z
      
      double precision acc
      integer itmax,ibegin,iwarn(5)
      DATA ACC/1.D-14/,ITMAX/6/,IBEGIN/0/,IWARN/5*0/

      integer i,iter,k,nm
      double precision a,accu,bq,c,f,f0,g,g0,rmas,s,wt2,wt3,wtm,x,x2,
     .     xmax,xmt,rn

C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      Z(2)=PO2LOG
      DO 101 K=3,9
  101 Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))
      DO 102 K=3,9
  102 Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))

C CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.10) GOTO 104
      PRINT 1001,N
      STOP

C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.
      NM=0
      DO 105 I=1,N
      IF(XM(I).NE.0.D0) NM=NM+1
  105 XMT=XMT+DABS(XM(I))
      IF(XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      STOP

C THE PARAMETER VALUES ARE NOW ACCEPTED

C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
  201 DO 202 I=1,N
      C=2.*RN(1)-1.
      S=DSQRT(1.-C*C)
      F=TWOPI*RN(2)
      Q(4,I)=-DLOG(RN(3)*RN(4))
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*DCOS(F)
  202 Q(1,I)=Q(4,I)*S*DSIN(F)

C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1./(1.+G)
      X=ET/RMAS

C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)

C CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.*N-4.)*DLOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
C      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
C      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1

C RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=DEXP(WT)
      RETURN

C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=DSQRT(1.-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.
      X2=X*X
      DO 303 I=1,N
      E(I)=DSQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(DABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
      PRINT 1006,ITMAX
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)

C CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.
      WT3=0.
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.*N-3.)*DLOG(X)+DLOG(WT2/WT3*ET)

C RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IF(IWARN(4).LE.5) PRINT 1005,WT
      IWARN(4)=IWARN(4)+1
  310 WT=DEXP(WT)
      RETURN

 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
     . ' DESIRED ACCURACY =',D15.6)

      END

************************************************************************

      function random(seed)
C---Removed implicit declarations
C---Gabor Somogyi March 2014
      implicit none
      integer m,a,q,r,hi,lo,seed
      double precision minv,random
      save
      parameter(m=2147483647,a=16807,q=127773,r=2836)
      parameter(minv=0.46566128752458d-09)
      hi = seed/q
      lo = mod(seed,q)
      seed = a*lo - r*hi
      if(seed.le.0) seed = seed + m
      random = seed*minv
      end

************************************************************************

c      double precision FUNCTION RN(IDUMMY)
cC---Removed implicit declarations
cC---Gabor Somogyi March 2014
c      implicit none
c      include 'seed.f'
c      INTEGER IDUMMY
c      double precision RANDOM
c      RN=RANDOM(ISEED)
c      RETURN
c      END
