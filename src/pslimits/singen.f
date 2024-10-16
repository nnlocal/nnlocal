************************************************************************************
* Phase-space generator of singular (soft and/or collinear) phase space points for *
* an arbitrary number of final state particles. The level of softness/collinearity *
* is tunable by the user.                                                          *
* Based on the RAMBO generator and the Nagy-Soper mapping from a n-parton to a     *
* (n+1)-parton phase space.                                                        *
*                                                                                  *
* INPUTS:                                                                          *
*  - n = number of external particles (in the n-parton phase space)                *
*  - rmass_split(1:2) = list of masses for the two splitting daughters.            *
*                       [NOTE: second element is assumed to be mj, i.e. the mass   *
*                       of the unconstrained parton]                               *
*  - pmomtilde(n,1:5) = array of momenta in the (m+1)-parton phase space           *
*  - ipos_emitter = position of the emitter in the n-particle list of partons      *
*  - P2tune = invariant-mass tuner. A number in the range [0,1].                   *
*             0/1 = minimum/maximum invariant mass allowed by kinematics           *
*  - softtune = soft tuner. A number in the range [0,1].                           *
*               0/1 = minimum/maximum energy of the unconstrained parton allowed   *
*               by kinematics                                                      *
*                                                                                  *
* OUTPUT:                                                                          *
*  - pmom(n+1,1:5) = array of momenta in the (n+1)-parton phase space              *
*                    Momenta generated in the center-of-mass frame.                *
*                    Last two momenta in the list belong to pi and pj              *
*                    (unconstrained parton) respectively                           *
*                                                                                  *
*                                                                                  *
* Author: G. Bevilacqua                                                            *
* Last update: 8 May 2014                                                          *
************************************************************************************


      subroutine singen(n,rmass_split,pmomtilde,ipos_emitter,
     .     P2tune,softtune,pmom)


      IMPLICIT NONE

      double precision aux
c... Input/output
      integer n                 !... number of external particles (n-kinematics)
      double precision rmass_split !... list of masses of the daughter particles
      dimension rmass_split(2)     !... --> 1=i ; 2=j (unconstrained)
      integer ipos_emitter      !... position of the emitter in rmass
      double precision P2tune   !... inv. mass tune (range: [0,1], max singular=0)
      double precision softtune !... softness tuner (range: [0,1], max softness=0)
      double precision pmomtilde     !... momenta in the (m)-kinematics
      double precision pmom     !... momenta in the (m+1)-kinematics
      dimension pmomtilde(100,1:5),pmom(100,1:5)
c... Internal variables
      integer i,mu,isign,ik
      double precision app,pmodul,xrandom
      double precision sqrts    !... CM energy for the process
      double precision rmass    !... list of masses (n-kinematics)
      dimension rmass(100)
      double precision Ktilde(1:5),K(1:5) !... collective spectator's momentum
      double precision P(1:4)   !... pi+pj
      double precision a,b,y,lambda
      double precision MK,P2min,P2max
      double precision Pmod,P2
      double precision pi3dotpj
      double precision Qsq,Ktsq
      double precision Q(1:4)
c...  Emitter
      double precision pitilde   !... 4-momentum
      dimension pitilde(1:4)
      double precision pitildemod !... modulus of 3-momentum
      double precision ctit,stit  !... cos and sin of the polar angle
      double precision cfit,sfit  !... cos and sin of the azimuthal angle
c...  Unconstrained parton (pj)
      double precision pj       !... 4-momentum
      dimension pj(1:5)
      double precision ctij,stij !... cos and sin of the angle between pi and pj
      double precision ctjit,stjit !... cos and sin of angle between pj and pitilde
      double precision cfjit,sfjit !... relative azimuth angle between pj and ptilde
      double precision ctjmin,ctjmax,ctj,stj !... cos and sin of the polar angle
      double precision cfj,sfj  !... cos and sin of the azimuthal angle
      double precision Ejmin,Ejmax,Ej !... Energy
      double precision pjmod    !... modulus of 3-momentum
      double precision betaj    !... pjmod/Ej
      double precision x        !... energy fraction for the initial-state splitting
c...  Other parton from the splitting (pi)
      double precision pi       !... 4-momentum
      dimension pi(1:5)
      double precision Ei       !... Energy
      double precision pimod    !... modulus of 3-momentum
      double precision betai    !... pimod/Ei
c...  Spectators
      double precision pktilde,pk
      dimension pktilde(100,1:5),pk(100,1:5)
c...  Initial state
      double precision pa,pb
      dimension pa(1:5),pb(1:5)
c...  "P2ness and softness"
      double precision P2ness,softness
c...  Boost parameters
      double precision z,ys,bcm,gcm,pmomtilde_boost
      dimension pmomtilde_boost(100,1:5)


      double precision RN
cx      external function RN
      external RN


      do i=1,n
         rmass(i) = pmomtilde(i,5)
      enddo
      sqrts = sqrt(
     .     (pmomtilde(1,4)+pmomtilde(2,4))**2-
     .     (pmomtilde(1,1)+pmomtilde(2,1))**2-
     .     (pmomtilde(1,2)+pmomtilde(2,2))**2-
     .     (pmomtilde(1,3)+pmomtilde(2,3))**2
     .     )


      if(ipos_emitter.ge.3)then

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                              FINAL-STATE SPLITTING                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c...  The initial state
      do mu=1,4
         pa(mu)=pmomtilde(1,mu)
         pb(mu)=pmomtilde(2,mu)
      enddo
      pa(5)=rmass(1)
      pb(5)=rmass(2)

c...  The emitter's 4-momentum
      do mu=1,4
         pitilde(mu)=pmomtilde(ipos_emitter,mu)
      enddo

c...  Compute momentum of collective spectator (m-parton kinematics)
      ik=0
      pktilde(1:n-3,1:5)=0.d0
      Ktilde(1:5)=0.d0
      do i=3,n
         if(i.ne.ipos_emitter)then
            ik=ik+1
            do mu=1,4
               pktilde(ik,mu)=pmomtilde(i,mu)
               Ktilde(mu)=Ktilde(mu)+pmomtilde(i,mu)
            enddo
            pktilde(ik,5)=rmass(i)
         endif
      enddo
      if(ik.eq.1)then
         MK=pktilde(1,5)
      else
         MK=sqrt(Ktilde(4)**2-Ktilde(1)**2-Ktilde(2)**2-Ktilde(3)**2)
      endif
      Ktilde(5)=MK

c...  Compute basic kinematical variables for the mapping together with their
c...  extrema allowed by kinematics

         P2min = (rmass_split(1)+rmass_split(2))**2
         P2max = (sqrts-MK)**2
         P2 = P2min + P2tune*(P2max-P2min)
         a = (sqrts**2)/(2.d0*sqrts*pitilde(4))
         b = rmass(ipos_emitter)**2/(2.d0*sqrts*pitilde(4))
         y = (P2-rmass(ipos_emitter)**2)/(2.d0*sqrts*pitilde(4))
         lambda = sqrt(((1.d0+y)**2-4.d0*a*(y+b))/(1.d0-4.d0*a*b))

c...  Total splitting momentum, P(mu) = pi(mu)+pj(mu), assuming total momentum Q 
c...  in the rest frame
         do mu=1,4
            if(mu.eq.4)then
               P(mu) = lambda*pitilde(mu) + 
     .              (1.d0-lambda+y)/(2.d0*a)*sqrts
            else
               P(mu) = lambda*pitilde(mu)
            endif
         enddo
         Pmod=sqrt(P(1)**2+P(2)**2+P(3)**2)

c...  Energy variable of the unconstrained parton (pj)
         Ejmin = (P2+rmass_split(2)**2-rmass_split(1)**2)/(2.d0*P2)*
     .        (P(4)-Pmod*sqrt(1.d0-(4.d0*P2*rmass_split(2)**2)/
     .        ((P2+rmass_split(2)**2-rmass_split(1)**2)**2)))
         Ejmax = (P2+rmass_split(2)**2-rmass_split(1)**2)/(2.d0*P2)*
     .        (P(4)+Pmod*sqrt(1.d0-(4.d0*P2*rmass_split(2)**2)/
     .        ((P2+rmass_split(2)**2-rmass_split(1)**2)**2)))
         Ej = Ejmin + softtune*(Ejmax-Ejmin)

        if(Ej.lt.0.d0)then
            print*,'ERROR: negative energy for pj !!!',Ej
            STOP
         endif

c...  Having fixed the energy, the modulus of pj is constrained by onshellness
         pjmod=sqrt(Ej**2-rmass_split(2)**2)

c...  Angular parameters of the emitter (pitilde)
         pitildemod = sqrt(pitilde(1)**2+pitilde(2)**2+pitilde(3)**2)
         ctit = pitilde(3)/pitildemod
         stit = sqrt(1.d0-ctit**2)
         call azimuth(pitilde,cfit,sfit)

c...  Determine the energy of the splitting parton pi.
c...  The condition P = pi+pj provides the constraint.
         Ei = P(4)-Ej

c...  Having fixed the energy, the modulus of pi is constrained by onshellness
         pimod = sqrt(Ei**2-rmass_split(1)**2)

c...  Determine the angle between the two splitting partons pi and pj.
c...  The condition (pi+pj)**2 = P2 provides the constraint.
         betai = pimod/Ei
         betaj = pjmod/Ej
         ctij = (2.d0*Ei*Ej+rmass_split(1)**2+rmass_split(2)**2-P2)/
     .        (2.d0*betai*betaj*Ei*Ej)
         stij=sqrt(1.d0-ctij**2)

c...  Determine the angle between pj and the emitter pitilde.
c...  The condition (P.pj) = lambda*(pitilde.pj) [here "." denotes the 
c...  3-dot product] provides the constraint.
         pi3dotpj = pimod*pjmod*ctij
         ctjit = (pi3dotpj + pjmod**2)/(lambda*pitildemod*pjmod)
         stjit=sqrt(1.d0-ctjit**2)

c...  Determine angular parameters of the unconstrained parton: polar and azimuthal
c...  angle.

c...  Using the trigonometric identity: ctjit = stit*stj*cfjit + ctit*ctj   (Eq.1).
c...  We adopt the following procedure to include contribution of a random 
c...  azimuthal angle:
c...   1. generate a random value of polar angle for pj in the allowed range
c...      [ (theta_itilde-theta_jitilde) , (theta_itilde+theta_jitilde) ]:
c...                   ctj = ctj_min + xrandom*(ctj_max-ctj_min)
c...                   stj = sqrt(1-ctj**2)
c...   2. derive cfjit using Eq.1
c...   3. get the azimuthal angle of the unconstrained parton as
c...                   phi_j = phi_itilde + random_sign*phi_jit
c...      or equivalently, in terms of the cosine:
c...                   cfj = cfi*cfjit - random_sign*sfi*sfjit

         ctjmin = min(ctit*ctjit + stit*stjit , ctit*ctjit - stit*stjit)
         ctjmax = max(ctit*ctjit + stit*stjit , ctit*ctjit - stit*stjit)
         xrandom = RN(100.d0)      !... random number in the range [0,1]
         ctj = ctjmin + xrandom*(ctjmax-ctjmin)
         stj = sqrt(1.d0-ctj**2)
         cfjit = (ctjit-ctit*ctj)/(stit*stj)
         sfjit = sqrt(1.d0-cfjit**2)
         xrandom = RN(101.d0)      !... random number in the range [0,1]
         if(xrandom.le.0.5d0)then
            isign=-1
         else
            isign=1
         endif
         cfj = cfit*cfjit - sfit*sfjit*isign
         sfj = sfit*cfjit + cfit*sfjit*isign

c...  Determine the 4-momentum of the daughter partons (pj,pi)
         pj(1)=pjmod*stj*cfj
         pj(2)=pjmod*stj*sfj
         pj(3)=pjmod*ctj
         pj(4)=Ej
         pj(5)=rmass_split(2)
         do mu=1,4
            pi(mu)=P(mu)-pj(mu)
         enddo
         pi(5)=rmass_split(1)

c...  Determine momentum of collective spectator ((m+1)-parton kinematics)
      K(1:5)=0.d0
      do mu=1,4
         if(mu.eq.4)then
            K(mu)=sqrts-P(mu)
         else
            K(mu)=-P(mu)
         endif
      enddo
      K(5)=MK

c...  Boost spectators according to Lorentz transformation
      call Lorentz(n-3,K,Ktilde,pktilde,pk)

c...  Fill out final momenta in the (m+1)-parton kinematics
c...  Last two particles are the splitting partons pi and pj
      do mu=1,5
         pmom(1,mu)=pa(mu)
         pmom(2,mu)=pb(mu)
      enddo
      do i=3,n-1
         do mu=1,5
            pmom(i,mu)=pk(i-2,mu)
         enddo
      enddo
      do mu=1,5
         pmom(n,mu)=pi(mu)
         pmom(n+1,mu)=pj(mu)
      enddo


      elseif(ipos_emitter.le.2)then

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                            INITIAL-STATE SPLITTING                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         x = 1.d0-softtune
         Ktsq = sqrts**2
         Qsq = Ktsq/x

c...  Energy of the unconstrained parton in the frame where Q=(pa+pb) is at rest
         Ej = sqrt(Qsq)/2.d0*(1.d0-x)
         if(Ej.lt.0.d0)then
            print*,'ERROR: negative energy for pj !!!',Ej
            STOP
         endif

c...  Boost momenta
         z = 1.d0/x
         ys = -0.5d0*log(x)     !...rapidity
c...  If emitter is parton b (coming from right), rapidity has opposite sign
         if(ipos_emitter.eq.2) ys=-ys
         bcm = (exp(ys)-exp(-ys))/(exp(ys)+exp(-ys)) !...beta
         gcm = (exp(ys)+exp(-ys))/2.d0 !...gamma
c... Boost ptilde momenta so that the (m+1) process is in the center-of-mass frame
         do i=1,n
            pmomtilde_boost(i,4) = gcm*(pmomtilde(i,4)-bcm*
     .           pmomtilde(i,3))
            pmomtilde_boost(i,1) = pmomtilde(i,1)
            pmomtilde_boost(i,2) = pmomtilde(i,2)
            pmomtilde_boost(i,3) = pmomtilde(i,3)/gcm-bcm*
     .           pmomtilde_boost(i,4)
         enddo

c...  The initial state
         do mu=1,4
            pa(mu)=pmomtilde_boost(1,mu)
            pb(mu)=pmomtilde_boost(2,mu)
         enddo
         pa(5)=rmass(1)
         pb(5)=rmass(2)

c...  Compute momentum of collective spectator (m-parton kinematics)
c...  NOTE: in case of initial-state splitting, the number of spectators is (n-2)
c...  and not (n-3) like in the final-state case!
         pktilde(1:n-2,1:5)=0.d0
         Ktilde(1:5)=0.d0
         do i=3,n
            do mu=1,4
               pktilde(i-2,mu)=pmomtilde_boost(i,mu)
               Ktilde(mu)=Ktilde(mu)+pmomtilde_boost(i,mu)
            enddo
            pktilde(i-2,5)=rmass(i)
         enddo
         if(n.eq.3)then
            MK=pktilde(1,5)
         else
            MK=sqrt(Ktilde(4)**2-Ktilde(1)**2-Ktilde(2)**2-Ktilde(3)**2)
         endif
         Ktilde(5)=MK

c...  Determine invariant mass Paj^2 and collinear variable
         if(ipos_emitter.eq.1)then !... splitting parton: pa
            do mu=1,4
               pa(mu) = pa(mu)/x
            enddo
            P2min = 0.d0
            P2max = -4.d0*pa(4)*Ej
            P2 = P2min + P2tune*(P2max-P2min)
            ctjit = 1.d0 + P2/(2.d0*pa(4)*Ej)
         elseif(ipos_emitter.eq.2)then !... splitting parton: pb
            do mu=1,4
               pb(mu) = pb(mu)/x
            enddo
            P2min = 0.d0
            P2max = -4.d0*pb(4)*Ej
            P2 = P2min + P2tune*(P2max-P2min)
            ctjit = 1.d0 + P2/(2.d0*pb(4)*Ej)
         endif

c...  Reconstruct the 4-momentum of the unconstrained parton (pj)
c...  Generate a random azimuthal angle
cx         xrandom = RN(102)      !... random number in the range [0,1]
         xrandom = RN(102.d0)      !... random number in the range [0,1]
         cfj = -1.d0 + 2.d0*xrandom
cx         xrandom = RN(203)      !... random number in the range [0,1]
         xrandom = RN(203.d0)      !... random number in the range [0,1]
         if(xrandom.le.0.5d0)then
            isign=-1
         else
            isign=1
         endif
         sfj = sqrt(1.d0-cfj**2)*isign

         pjmod = sqrt(Ej**2-rmass_split(2)**2)
         stjit = sqrt(1.d0-ctjit**2)
c         pj(1) = 0.d0
c         pj(2) = Ej*stjit
         pj(1) = Ej*stjit*sfj
         pj(2) = Ej*stjit*cfj
         pj(3) = Ej*ctjit
         if(ipos_emitter.eq.2) pj(3)=-pj(3)
         pj(4) = Ej
         pj(5) = rmass_split(2)

c...  Total momentum of the (m+1) process
         do mu=1,4
            Q(mu) = pa(mu)+pb(mu)
         enddo

c...  Determine the collective spectator momentum K
         do mu=1,4
            K(mu) = Q(mu)-pj(mu)
         enddo
         K(5) = Ktilde(5)

c...  Boost spectators according to Lorentz transformation
c...  NOTE: in case of initial-state splitting, the number of spectators is (n-2)
c...  and not (n-3) like in the final-state case!
         call Lorentz(n-2,K,Ktilde,pktilde,pk)

c...  Fill out final momenta in the (m+1)-parton kinematics
c...  Last particle is the unconstrained parton pj
         do mu=1,5
            pmom(1,mu)=pa(mu)
            pmom(2,mu)=pb(mu)
         enddo
      do i=3,n
         do mu=1,5
            pmom(i,mu)=pk(i-2,mu)
         enddo
      enddo
      do mu=1,5
         pmom(n+1,mu)=pj(mu)
      enddo

      endif




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...  Print check of momentum conservation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C$$$      print*,' '
C$$$      print*,'CHECK OF MOMENTUM CONSERVATION -- m-parton kinematics'
C$$$      aux=0.d0
C$$$      do i=3,n
C$$$         aux=aux+pmomtilde(i,1)
C$$$      enddo
C$$$      print*,'px IN: ',pmomtilde(1,1)+pmomtilde(2,1),'px FIN: ',aux
C$$$      aux=0.d0
C$$$      do i=3,n
C$$$         aux=aux+pmomtilde(i,2)
C$$$      enddo
C$$$      print*,'py IN: ',pmomtilde(1,2)+pmomtilde(2,2),'py FIN: ',aux
C$$$      aux=0.d0
C$$$      do i=3,n
C$$$         aux=aux+pmomtilde(i,3)
C$$$      enddo
C$$$      print*,'pz IN: ',pmomtilde(1,3)+pmomtilde(2,3),'pz FIN: ',aux
C$$$      aux=0.d0
C$$$      do i=3,n
C$$$         aux=aux+pmomtilde(i,4)
C$$$      enddo
C$$$      print*,'En IN: ',pmomtilde(1,4)+pmomtilde(2,4),'En FIN: ',aux
C$$$      print*,' '

C$$$      print*,' '
C$$$      print*,'CHECK OF MOMENTUM CONSERVATION -- (m+1)-parton kinematics'
C$$$      aux=0.d0
C$$$      do i=3,n+1
C$$$         aux=aux+pmom(i,1)
C$$$      enddo
C$$$      print*,'px IN: ',pmom(1,1)+pmom(2,1),'px FIN: ',aux
C$$$      aux=0.d0
C$$$      do i=3,n+1
C$$$         aux=aux+pmom(i,2)
C$$$      enddo
C$$$      print*,'py IN: ',pmom(1,2)+pmom(2,2),'py FIN: ',aux
C$$$      aux=0.d0
C$$$      do i=3,n+1
C$$$         aux=aux+pmom(i,3)
C$$$      enddo
C$$$      print*,'pz IN: ',pmom(1,3)+pmom(2,3),'pz FIN: ',aux
C$$$      aux=0.d0
C$$$      do i=3,n+1
C$$$         aux=aux+pmom(i,4)
C$$$      enddo
C$$$      print*,'En IN: ',pmom(1,4)+pmom(2,4),'En FIN: ',aux
C$$$      print*,' '

C$$$      if(ipos_emitter.ge.3)then
C$$$         P2ness = sqrt(
C$$$     .        (pmom(n,4)+pmom(n+1,4))**2-
C$$$     .        (pmom(n,1)+pmom(n+1,1))**2-
C$$$     .        (pmom(n,2)+pmom(n+1,2))**2-
C$$$     .        (pmom(n,3)+pmom(n+1,3))**2
C$$$     .        )
C$$$         softness = pmom(n+1,4)
C$$$      elseif(ipos_emitter.eq.1)then
C$$$         P2ness = sqrt(
C$$$     .        abs(
C$$$     .        (pmom(1,4)-pmom(n+1,4))**2-
C$$$     .        (pmom(1,1)-pmom(n+1,1))**2-
C$$$     .        (pmom(1,2)-pmom(n+1,2))**2-
C$$$     .        (pmom(1,3)-pmom(n+1,3))**2
C$$$     .        )
C$$$     .        )
C$$$         softness = pmom(n+1,4)
C$$$      elseif(ipos_emitter.eq.2)then
C$$$         P2ness = sqrt(
C$$$     .        abs(
C$$$     .        (pmom(2,4)-pmom(n+1,4))**2-
C$$$     .        (pmom(2,1)-pmom(n+1,1))**2-
C$$$     .        (pmom(2,2)-pmom(n+1,2))**2-
C$$$     .        (pmom(2,3)-pmom(n+1,3))**2
C$$$     .        )
C$$$     .        )
C$$$         softness = pmom(n+1,4)
C$$$      endif

C$$$      print*,' '
C$$$      print*,'sqrts',sqrts
C$$$      print*,'Inv. mass',P2ness
C$$$      print*,'Ej',softness
C$$$      print*,'ctj',ctjit
C$$$      print*,' '
c



      return
      end



c-----------------------------------------------------------------------------------


c-----------------------------------------------------------------------------------

      subroutine azimuth(p,cosphi,sinphi)

      IMPLICIT NONE
      double precision p,px,py,cosphi,sinphi
      dimension p(1:4)

      px=p(1)
      py=p(2)

      if(px.eq.0.d0 .AND. py.eq.0.d0)then
         cosphi=1.d0
         sinphi=0.d0
         return
      elseif(px.gt.0.d0 .AND. py.eq.0.d0)then
         cosphi=1.d0
         sinphi=0.d0
         return
      elseif(px.eq.0.d0 .AND. py.gt.0.d0)then
         cosphi=0.d0
         sinphi=1.d0
         return
      elseif(px.lt.0.d0 .AND. py.eq.0.d0)then
         cosphi=-1.d0
         sinphi=0.d0
         return
      elseif(px.eq.0.d0 .AND. py.lt.0.d0)then
         cosphi=0.d0
         sinphi=-1.d0
         return
      else
         if(px.gt.0.d0)then
            cosphi=cos(atan(py/px))
            if(py.gt.0.d0)then
               sinphi=sqrt(1.d0-cosphi**2)
            else
               sinphi=-sqrt(1.d0-cosphi**2)
            endif
         else                   !... px.lt.0
            cosphi=-cos(atan(py/px))
            if(py.gt.0.d0)then
               sinphi=sqrt(1.d0-cosphi**2)
            else
               sinphi=-sqrt(1.d0-cosphi**2)
            endif
         endif
      endif

      return

      end




c-----------------------------------------------------------------------

      subroutine Lorentz(nspec,K,Ktilde,ptk,pk)

      integer nspec,i,mu
      double precision K,Ktilde,ptk,pk
      dimension K(1:5),Ktilde(1:5),ptk(100,1:5),pk(100,1:5)
      double precision K2,Ktilde2,Kdotptk,Ktildedotptk,KdotKtilde
c      external dot4


      if(nspec.eq.1)then

         do mu=1,4
            pk(nspec,mu) = K(mu)
         enddo
         pk(nspec,5) = ptk(nspec,5)

      elseif(nspec.gt.1)then

         K2 = K(5)**2
         Ktilde2 = Ktilde(5)**2
         KdotKtilde = K(4)*Ktilde(4)-K(1)*Ktilde(1)-K(2)*Ktilde(2)-
     .        K(3)*Ktilde(3)

         do i=1,nspec
            Kdotptk = K(4)*ptk(i,4) - K(1)*ptk(i,1) - K(2)*ptk(i,2) - 
     .           K(3)*ptk(i,3)
            Ktildedotptk = Ktilde(4)*ptk(i,4) - Ktilde(1)*ptk(i,1) -
     .           Ktilde(2)*ptk(i,2) - Ktilde(3)*ptk(i,3)
            do mu=1,4
               pk(i,mu) = ptk(i,mu)
     .              -2.d0*(Kdotptk + Ktildedotptk)*
     .              (K(mu)+Ktilde(mu))/(K2 + Ktilde2 +
     .              2.d0*KdotKtilde) + 2.d0*Ktildedotptk*K(mu)/
     .              Ktilde2
            enddo
            pk(i,5) = ptk(i,5)
         enddo

      else

         print*,'ERROR in LORENTZ subroutine'
         STOP

      endif

      return
      end

c-----------------------------------------------------------------------
************************************************************************************
*     Some routines to translate between conventions for singen
*     by Gabor Somogyi
************************************************************************************

      subroutine pTOsingen(pin,pout)
      implicit none

      include 'constants.f'

      double precision pin(4,mxpart)
      double precision pout(100,1:5)

      double precision masses(mxpart)
      common/masses/masses

      integer j,k

      do j=1,2
         do k=1,4
            pout(j,k) = -pin(k,j)
         enddo
         pout(j,5) = 0d0
      enddo
      
      do j=3,mxpart
         do k=1,4
            pout(j,k) = pin(k,j)
         enddo
         pout(j,5) = masses(j-2)
      enddo

      end


      subroutine pFROMsingen(pin,pout)
      implicit none

      include 'constants.f'

      double precision pin(100,1:5)
      double precision pout(4,mxpart)

      integer j,k

      do j=1,2
         do k=1,4
            pout(k,j) = -pin(j,k)
         enddo
      enddo
      
      do j=3,mxpart
         do k=1,4
            pout(k,j) = pin(j,k)
         enddo
      enddo
      
      end
