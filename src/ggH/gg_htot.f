      subroutine gg_htot(z,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c    g(-p1)+g(-p2)-->H --> gamma(p3) + gamma(p4) 
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'order.f'
      include 'scale.f'
      include 'facscale.f'
      integer i,j,k
      double precision msqlo(-nf:nf,-nf:nf,4), msqnlo(-nf:nf,-nf:nf,4),
     1     msqnnlo(-nf:nf,-nf:nf,4), msq(-nf:nf,-nf:nf,4), z, s12, lt
      double precision decay,gg,Asq,msqgamgam,hp_wgplg
      double precision b0, b1, cw1, cw2, logr, logf
      integer switch

!---  switch = 0  : Anastasiou, Melnikov
!---  switch else : Mistlberger+...
      switch = 1
      

      b0 = 11d0/4d0             !-nf/6d0
      b1 = 51d0/8d0             !-19d0*nf/24d0
      lt=log(facscale**2/mt**2)
      cw1 = 11d0/4d0
      cw2 = (2777 + 6*lt*57)/288d0
      logr = -dlog(facscale**2/musq)
      logf = dlog(hmass**2/facscale**2)
      
c---  set msq=0 to initialize
      do i=1,4
         do j=-nf,nf
            do k=-nf,nf
               msqlo(j,k,i)   = 0d0
               msqnlo(j,k,i)  = 0d0
               msqnnlo(j,k,i) = 0d0
               msq(j,k,i)     = 0d0
            enddo
         enddo
      enddo

      Asq=(as/(pi))**2/vevsq*pi
      gg=Asq/576d0
      
c---  order-by-order
c---  lo
      msqlo(0,0,1)   = gg

      if (order.eq.0) then
         msq(0,0,1)=msqlo(0,0,1)
         msq(0,0,2)=msqlo(0,0,2)
         msq(0,0,3)=msqlo(0,0,3)
      endif

      if (order.eq.0) return

      
c---  nlo
      if (switch.eq.0) then
         msqnlo(0,0,1)  = gg*as/pi*(5.5d0+pisq)

         msqnlo(0,0,2)  = gg*as/pi*(12d0*log(1-z)/(1-z) + 6d0*logf/(1d0-z))
         
         msqnlo(0,0,3)  = gg*as/pi*(-5.5d0*(1-z)**3-12d0*z*(2d0-z+z**2)*log(1-z)
     1        -6d0*(1-z+z**2)**2*log(z)/(1-z)
     2        +logf*(0d0 - 6d0*z*(2d0 + (-1d0 + z)*z)))
      else
         msqnlo(0,0,1)  = gg*as/pi*(pisq)

         msqnlo(0,0,2)  = gg*as/pi*(12d0*log(1-z)/(1-z) + 6d0*logf/(1d0-z))
         
         msqnlo(0,0,3)  = gg*as/pi*((-11*(1 - z)**3)/2d0
     1        + 6*logf*(- 2*z + z**2 - z**3) + 12*(- 2*z + z**2 - z**3)*Log(1 - z)
     2        - (6*(1 - z + z**2)**2*Log(z))/(1 - z))
      endif


      if (order.eq.1) then
         if (switch.eq.0) then
            msq(0,0,1)=msqlo(0,0,1)*(1d0 + 2d0*b0*logr*(as/pi)) + msqnlo(0,0,1)
            msq(0,0,2)=msqlo(0,0,2)*(1d0 + 2d0*b0*logr*(as/pi)) + msqnlo(0,0,2)
            msq(0,0,3)=msqlo(0,0,3)*(1d0 + 2d0*b0*logr*(as/pi)) + msqnlo(0,0,3)
         else
            msq(0,0,1)=msqlo(0,0,1)*(1d0 + 2d0*(b0*logr+cw1)*(as/pi)) + msqnlo(0,0,1)
            msq(0,0,2)=msqlo(0,0,2)*(1d0 + 2d0*(b0*logr+cw1)*(as/pi)) + msqnlo(0,0,2)
            msq(0,0,3)=msqlo(0,0,3)*(1d0 + 2d0*(b0*logr+cw1)*(as/pi)) + msqnlo(0,0,3)
         endif
      endif

      if (order.eq.-1) then
         if (switch.eq.0) then
            msq(0,0,1)=msqlo(0,0,1)*(2d0*b0*logr*(as/pi)) + msqnlo(0,0,1)
            msq(0,0,2)=msqlo(0,0,2)*(2d0*b0*logr*(as/pi)) + msqnlo(0,0,2)
            msq(0,0,3)=msqlo(0,0,3)*(2d0*b0*logr*(as/pi)) + msqnlo(0,0,3)
         else
            msq(0,0,1)=msqlo(0,0,1)*(2d0*(b0*logr+cw1)*(as/pi)) + msqnlo(0,0,1)
            msq(0,0,2)=msqlo(0,0,2)*(2d0*(b0*logr+cw1)*(as/pi)) + msqnlo(0,0,2)
            msq(0,0,3)=msqlo(0,0,3)*(2d0*(b0*logr+cw1)*(as/pi)) + msqnlo(0,0,3)
         endif
      endif

      if (abs(order).eq.1) return

      
c---  nnlo
      if (switch.eq.0) then
      msqnnlo(0,0,1) = gg*(as/pi)**2*
     1   (11399/144+19*lt/8d0+133*pisq/12d0-pisq**2/80d0-165*zeta3/4d0)
      msqnnlo(0,0,2) = gg*(as/pi)**2*
     1  (((133 - 15*pisq)*log(1-z))/(1-z) 
     2   -(33*log(1-z)**2)/(1-z) + (72*log(1-z)**3)/(1-z)
     3   +(-101/3d0 + (11*pisq)/2d0 + (351*zeta3)/2d0)/(1-z))
      msqnnlo(0,0,3) = gg*(as/pi)**2*(
     1   -(18157/48) + (22879*z)/48d0 - (24107*z**2)/48d0 + (7539*z**3)/16d0 + 
     2  (Pi**2*(16 - 200*z + 564*z**2 - 886*z**3 + 495*z**4))/(8d0*(1 - z)) + 
     3  (506.75d0 - (2687*z)/4d0 + (1079*z**2)/2d0 - (2559*z**3)/4d0)*Log(1 - z) - 
     4  (3*Pi**2*(3 - 46*z + 69*z**2 - 48*z**3 + 23*z**4)*Log(1 - z))/(4d0*(1 - z)) + 
     5  (-297 + 381*z - 348*z**2 + 330*z**3)*Log(1 - z)**2 - 
     6  72*z*(2 - z + z**2)*Log(1 - z)**3 - 
     7  ((2333 - 5611*z + 6447*z**2 - 6955*z**3 + 4318*z**4)*Log(z))/(8d0*(1 - z)) + 
     8  (3*Pi**2*(2 + 6*z + 18*z**2 - 6*z**3 - 19*z**4 + 10*z**5)*Log(z))/
     9   (2d0*(1 - z**2)) + (3*(641 - 1238*z + 1227*z**2 - 1180*z**3 + 605*z**4)*
     1     Log(1 - z)*Log(z))/(4d0*(1 - z)) - 
     2  (9*(59 - 118*z + 177*z**2 - 116*z**3 + 59*z**4)*Log(1 - z)**2*Log(z))/
     3   (4d0*(1 - z)) - (3*(154 - 389*z + 699*z**2 - 827*z**3 + 374*z**4)*Log(z)**2)/
     4   (8d0*(1 - z)) + (9*(6 - 12*z + 18*z**2 - 11*z**3 + 6*z**4)*Log(1 - z)*
     5     Log(z)**2)/(1 - z) - (3*(4 + 7*z + 18*z**2 - 7*z**3 - 17*z**4 + 9*z**5)*
     6     Log(z)**3)/(2d0*(1 - z**2)) + 
     7  (3*Pi**2*(11 + 22*z + 33*z**2 + 16*z**3 + 8*z**4)*Log(1 + z))/(4d0*(1 + z)) + 
     8  (-21 - 18*z + (63*z**2)/2d0 + (33*z**3)/2d0)*Log(z)*Log(1 + z) - 
     9  (9*(9 + 18*z + 27*z**2 + 8*z**3 + 4*z**4)*Log(z)**2*Log(1 + z))/(4d0*(1 + z)) + 
     1  (27*(1 + 2*z + 3*z**2)*Log(z)*Log(1 + z)**2)/(2d0*(1 + z)) + 
     2  (-21 - 18*z + (63*z**2)/2d0 + (33*z**3)/2d0)*hp_wgplg(1,1,-z) + 
     3  (9*(-5 - 10*z - 15*z**2 + 2*z**4)*Log(z)*hp_wgplg(1,1,-z))/(2d0*(1 + z)) + 
     4  (27*(1 + 2*z + 3*z**2)*Log(1 + z)*hp_wgplg(1,1,-z))/(1 + z) - 
     5  (3*(-289 + 394*z + 21*z**2 - 280*z**3 + 143*z**4)*hp_wgplg(1,1,z))/(4d0*(1 - z)) + 
     6  (9*(3 - 6*z + 9*z**2 - 8*z**3 + 3*z**4)*Log(1 - z)*hp_wgplg(1,1,z))/(2d0*(1 - z)) - 
     7  (9*(16 + 36*z + 64*z**2 - 40*z**3 - 67*z**4 + 13*z**5)*Log(z)*hp_wgplg(1,1,z))/
     8   (2d0*(1 - z**2)) - (36*(1 + z + z**2)**2*Log(1 + z)*hp_wgplg(1,1,z))/(1 + z) - 
     9  (9*(-1 - 2*z - 3*z**2 + 8*z**3 + 8*z**4)*hp_wgplg(2,1,-z))/(2d0*(1 + z)) + 
     1  (9*(14 + 18*z + 38*z**2 - 20*z**3 - 39*z**4 + 7*z**5)*hp_wgplg(2,1,z))/(1 - z**2) + 
     2  (9*(7 + 14*z + 21*z**2 + 8*z**3 + 4*z**4)*hp_wgplg(1,2,-z))/(1 + z) - 
     3  (9*(-11 + 59*z + 53*z**2 - 57*z**3 - 51*z**4 + 5*z**5)*hp_wgplg(1,2,z))/
     4   (2d0*(1 - z**2)) - (18*(1 + z + z**2)**2*hp_wgplg(1,2,z**2))/(1 + z) + 
     5  (9*(-36 - 52*z + 19*z**2 + 13*z**3 - 15*z**4 + 33*z**5)*zeta3)/
     6   (2d0*(1 - z**2)))
         
      else
         
         msqnnlo(0,0,1) = gg*(as/pi)**2*
     1        (837/16d0 + 67*pisq/12d0 - pisq**2/80d0 - 165*zeta3/4d0
     2        -3d0*logf**2*pisq
     3        + logf*(-54 - 11*pisq + 342*zeta3)/4d0)
         
         msqnnlo(0,0,2) = gg*(as/pi)**2*
     1        ((67 - 15*pisq)*log(1-z)/(1-z) 
     2        -33*log(1-z)**2/(1-z) + 72*log(1-z)**3/(1-z)
     3        +(-101/3d0 + 11*pisq/2d0 + 351*zeta3/2d0)/(1-z)
     4        +(9d0*logf**2*(-33) - 6d0*logf*(-201 + 45*pisq))/36d0/(1-z)
     5        +(-33d0*logf + 36d0*logf**2)*log(1-z)/(1-z)
     6        +(108d0*logf)*log(1-z)**2/(1-z))


         msqnnlo(0,0,3) = gg*(as/pi)**2*(
     1        (logf*(2161 - 1781*z + 1513*z**2 - 2161*z**3))/8d0 + 
     2        (3*logf**2*(-99 + 94*z - 83*z**2 + 99*z**3))/4d0 + 
     3        (-18321 + 18523*z - 19751*z**2 + 21165*z**3)/48d0 + 
     4        (3*logf*Pi**2*(-4 - 9*z - 24*z**2 - 14*z**3 + 6*z**4))/(2d0*(1 + z)) + 
     5        (Pi**2*(60 - 244*z + 564*z**2 - 886*z**3 + 495*z**4))/(8 - 8*z) + 
     6        (12*(1 + z + z**2)**2*Log(2d0)**3)/(1 + z) - 
     7        (3*Pi**2*(1 + z + z**2)**2*Log(4d0))/(1 + z) + 
     8        ((2295 - 2159*z + 1894*z**2 - 2295*z**3)*Log(1 - z))/4d0 - 
     9        36*logf**2*(-1 + 2*z - z**2 + z**3)*Log(1 - z) + 
     1        logf*(-330 + 381*z - 348*z**2 + 330*z**3)*Log(1 - z) + 
     2        (3*Pi**2*(27 - 39*z + 7*z**2 + 17*z**3 - 29*z**4 + 19*z**5)*Log(1 - z))/
     3        (4d0*(-1 + z**2)) - (18*(1 + z + z**2)**2*Log(2d0)**2*Log(1 - z))/(1 + z) - 
     4        108*logf*(-1 + 2*z - z**2 + z**3)*Log(1 - z)**2 + 
     5        (-330 + 381*z - 348*z**2 + 330*z**3)*Log(1 - z)**2 - 
     6        (6*(1 + z + z**2)**2*Log((1 - z)/8d0)*Log(1 - z)**2)/(1 + z) - 
     7        72*(-1 + 2*z - z**2 + z**3)*Log(1 - z)**3 + 
     8        (18*logf**2*(1 + 3*z**2 - 4*z**3 + z**4)*Log(z))/(-1 + z) + 
     9        (logf*(231 - 537*z + 837*z**2 - 1059*z**3 + 561*z**4)*Log(z))/(2 - 2*z) + 
     1        ((2069 - 5083*z + 5655*z**2 - 6427*z**3 + 4054*z**4)*Log(z))/(8d0*(-1 + z)) + 
     2        (Pi**2*(6 + 18*z + 54*z**2 - 18*z**3 - 57*z**4 + 30*z**5)*Log(z))/
     3        (2 - 2*z**2) + (126*logf*(1 - z + z**2)**2*Log(1 - z)*Log(z))/(-1 + z) + 
     4        ((1923 - 3714*z + 3681*z**2 - 3540*z**3 + 1815*z**4)*Log(1 - z)*Log(z))/
     5        (4 - 4*z) + (9*(12 + 28*z**2 + z**3 - 27*z**4 + 16*z**5)*Log(1 - z)**2*Log(z))/
     6        (-1 + z**2) + (3*(154 - 389*z + 699*z**2 - 827*z**3 + 374*z**4)*Log(z)**2)/
     7        (8d0*(-1 + z)) - (9*logf*(2 + 3*z + 8*z**2 - 3*z**3 - 8*z**4 + 3*z**5)*
     8        Log(z)**2)/(-1 + z**2) - 
     9        (9*(6 - 12*z + 18*z**2 - 11*z**3 + 6*z**4)*Log(1 - z)*Log(z)**2)/(-1 + z) + 
     1        (3*(4 + 7*z + 18*z**2 - 7*z**3 - 17*z**4 + 9*z**5)*Log(z)**3)/
     2        (2d0*(-1 + z**2)) + (9*Pi**2*(3 + 6*z + 9*z**2 + 8*z**3 + 4*z**4)*Log(1 + z))/
     3        (4d0*(1 + z)) - (18*(1 + z + z**2)**2*Log(2d0)**2*Log(1 + z))/(1 + z) + 
     4        (18*logf*(1 + z + z**2)**2*Log(z)*Log(1 + z))/(1 + z) + 
     5        (3*(-14 - 12*z + 21*z**2 + 11*z**3)*Log(z)*Log(1 + z))/2d0 - 
     6        (9*(9 + 18*z + 27*z**2 + 8*z**3 + 4*z**4)*Log(z)**2*Log(1 + z))/(4d0*(1 + z)) + 
     7        (3*(1 + z + z**2)**2*Log(64d0)*Log(1 + z)**2)/(1 + z) + 
     8        (9*(1 + 2*z + 3*z**2)*Log(1 + z)**3)/(2d0*(1 + z)) - 
     9        (9*(3 - 70*z + 9*z**2 + 56*z**3 + 3*z**4)*Log(1 - z)*hp_wgplg(1,1,1 - z))/
     1        (2d0*(-1 + z)) + (18*logf*(1 + z + z**2)**2*hp_wgplg(1,1,-z))/(1 + z) + 
     2        (3*(-14 - 12*z + 21*z**2 + 11*z**3)*hp_wgplg(1,1,-z))/2d0 + 
     3        (36*(1 + z + z**2)**2*Log(1 - z)*hp_wgplg(1,1,-z))/(1 + z) + 
     4        (9*(-5 - 10*z - 15*z**2 + 2*z**4)*Log(z)*hp_wgplg(1,1,-z))/(2d0*(1 + z)) + 
     5        144*logf*z*(1 + z)*hp_wgplg(1,1,z) + 
     6        (3*(-289 + 394*z + 21*z**2 - 280*z**3 + 143*z**4)*hp_wgplg(1,1,z))/(4d0*(-1 + z)) - 
     7        (9*(3 - 6*z + 9*z**2 - 8*z**3 + 3*z**4)*Log(1 - z)*hp_wgplg(1,1,z))/
     8        (2d0*(-1 + z)) + (9*(16 + 36*z + 64*z**2 - 40*z**3 - 67*z**4 + 13*z**5)*Log(z)*
     9        hp_wgplg(1,1,z))/(2d0*(-1 + z**2)) - 
     1        (36*(1 + z + z**2)**2*hp_wgplg(2,1,(1 - z)/2d0))/(1 + z) + 
     2        (9*(-5 - 75*z - 69*z**2 + 73*z**3 + 67*z**4 + 11*z**5)*hp_wgplg(2,1,1 - z))/
     3        (2d0*(-1 + z**2)) - (9*(7 + 14*z + 21*z**2 + 24*z**3 + 16*z**4)*hp_wgplg(2,1,-z))/
     4        (2d0*(1 + z)) - (9*(14 + 18*z + 38*z**2 - 20*z**3 - 39*z**4 + 7*z**5)*
     5        hp_wgplg(2,1,z))/(-1 + z**2) + 
     6        (36*(1 + z + z**2)**2*hp_wgplg(2,1,(2*z)/(-1 + z)))/(1 + z) - 
     7        (27*(1 + 2*z + 3*z**2)*hp_wgplg(2,1,1/(1 + z)))/(1 + z) - 
     8        (36*(1 + z + z**2)**2*hp_wgplg(2,1,z/(1 + z)))/(1 + z) + 
     9        (36*(1 + z + z**2)**2*hp_wgplg(2,1,(2*z)/(1 + z)))/(1 + z) - 
     1        (36*(1 + z + z**2)**2*hp_wgplg(2,1,(1 + z)/2d0))/(1 + z) - 
     2        (9*(18 - 107*z - 69*z**2 + 54*z**3 + 38*z**4 + 30*z**5)*zeta3)/
     3        (2d0*(-1 + z**2))
     .        -((67 - 15*pisq)*log(1-z)
     2        -33*log(1-z)**2 + 72*log(1-z)**3
     3        +(-101/3d0 + 11*pisq/2d0 + 351*zeta3/2d0)
     4        +(9d0*logf**2*(-33) - 6d0*logf*(-201 + 45*pisq))/36d0
     5        +(-33d0*logf + 36d0*logf**2)*log(1-z)
     6        +(108d0*logf)*log(1-z)**2))
         
      endif

c---  full


      if (order.eq.2) then
         if (switch.eq.0) then
            msq(0,0,1)=msqlo(0,0,1)*(1d0 + 2d0*b0*logr*(as/pi)
     1           + (2d0*b1*logr+3d0*b0**2*logr**2)*(as/pi)**2)
     2           + msqnlo(0,0,1)*(1d0+3d0*b0*logr*(as/pi))
     3           + msqnnlo(0,0,1)
            msq(0,0,2)=msqlo(0,0,2)*(1d0 + 2d0*b0*logr*(as/pi)
     1           + (2d0*b1*logr+3d0*b0**2*logr**2)*(as/pi)**2)
     2           + msqnlo(0,0,2)*(1d0+3d0*b0*logr*(as/pi))
     3           + msqnnlo(0,0,2)
            msq(0,0,3)=msqlo(0,0,3)*(1d0 + 2d0*b0*logr*(as/pi)
     1           + (2d0*b1*logr+3d0*b0**2*logr**2)*(as/pi)**2)
     2           + msqnlo(0,0,3)*(1d0+3d0*b0*logr*(as/pi))
     3           + msqnnlo(0,0,3)
         else

            msq(0,0,1) = msqnnlo(0,0,1)
     1           + msqnlo(0,0,1)*(1d0 + (as/pi)*(2d0*cw1 + 3d0*b0*logr))
     2           + msqlo(0,0,1)*(1d0 + (as/pi)*(2d0*cw1 + 2d0*b0*logr)
     3           + (as/pi)**2*(cw1**2 + 2d0*b1*logr + 6d0*b0*cw1*logr
     4           + 3d0*b0**2*logr**2 + 2d0*cw2))
            msq(0,0,2) = msqnnlo(0,0,2)
     1           + msqnlo(0,0,2)*(1d0 + (as/pi)*(2d0*cw1 + 3d0*b0*logr))
     2           + msqlo(0,0,2)*(1d0 + (as/pi)*(2d0*cw1 + 2d0*b0*logr)
     3           + (as/pi)**2*(cw1**2 + 2d0*b1*logr + 6d0*b0*cw1*logr
     4           + 3d0*b0**2*logr**2 + 2d0*cw2))
            msq(0,0,3) = msqnnlo(0,0,3)
     1           + msqnlo(0,0,3)*(1d0 + (as/pi)*(2d0*cw1 + 3d0*b0*logr))
     2           + msqlo(0,0,3)*(1d0 + (as/pi)*(2d0*cw1 + 2d0*b0*logr)
     3           + (as/pi)**2*(cw1**2 + 2d0*b1*logr + 6d0*b0*cw1*logr
     4           + 3d0*b0**2*logr**2 + 2d0*cw2))
            
         endif
      endif

      if (order.eq.-2) then
         if (switch.eq.0) then
            msq(0,0,1)=msqlo(0,0,1)*(
     1           + (2d0*b1*logr+3d0*b0**2*logr**2)*(as/pi)**2)
     2           + msqnlo(0,0,1)*(3d0*b0*logr*(as/pi))
     3           + msqnnlo(0,0,1)
            msq(0,0,2)=msqlo(0,0,2)*(
     1           + (2d0*b1*logr+3d0*b0**2*logr**2)*(as/pi)**2)
     2           + msqnlo(0,0,2)*(3d0*b0*logr*(as/pi))
     3           + msqnnlo(0,0,2)
            msq(0,0,3)=msqlo(0,0,3)*(
     1           + (2d0*b1*logr+3d0*b0**2*logr**2)*(as/pi)**2)
     2           + msqnlo(0,0,3)*(3d0*b0*logr*(as/pi))
     3           + msqnnlo(0,0,3)
         else

            msq(0,0,1) = msqnnlo(0,0,1)
     1           + msqnlo(0,0,1)*((as/pi)*(2d0*cw1 + 3d0*b0*logr))
     2           + msqlo(0,0,1)*(
     3           + (as/pi)**2*(cw1**2 + 2d0*b1*logr + 6d0*b0*cw1*logr
     4           + 3d0*b0**2*logr**2 + 2d0*cw2))
            msq(0,0,2) = msqnnlo(0,0,2)
     1           + msqnlo(0,0,2)*((as/pi)*(2d0*cw1 + 3d0*b0*logr))
     2           + msqlo(0,0,2)*(
     3           + (as/pi)**2*(cw1**2 + 2d0*b1*logr + 6d0*b0*cw1*logr
     4           + 3d0*b0**2*logr**2 + 2d0*cw2))
            msq(0,0,3) = msqnnlo(0,0,3)
     1           + msqnlo(0,0,3)*((as/pi)*(2d0*cw1 + 3d0*b0*logr))
     2           + msqlo(0,0,3)*(
     3           + (as/pi)**2*(cw1**2 + 2d0*b1*logr + 6d0*b0*cw1*logr
     4           + 3d0*b0**2*logr**2 + 2d0*cw2))
            
         endif
      endif
      
      return
      end


