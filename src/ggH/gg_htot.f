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
      integer i,j,k
      double precision msq(-nf:nf,-nf:nf,4),z,s12,lt
      double precision decay,gg,Asq,msqgamgam,hp_wgplg
      
c---  set msq=0 to initialize
      do i=1,4
         do j=-nf,nf
            do k=-nf,nf
               msq(j,k,i)=0d0
            enddo
         enddo
      enddo

      Asq=(as/(pi))**2/vevsq*pi
      gg=Asq/576

      
c---calculate propagators
      msq(0,0,1)=gg

      if (order.ge.1) then
         
      msq(0,0,1)=msq(0,0,1)+gg*as/pi*(5.5d0+pisq)
      msq(0,0,2)=gg*as/pi*(12d0*log(1-z)/(1-z))
      msq(0,0,3)=gg*as/pi*(-5.5d0*(1-z)**3-12d0*z*(2d0-z+z**2)*log(1-z)
     1                     -6d0*(1-z+z**2)**2*log(z)/(1-z))

      endif

      if (order.eq.2) then
c         lt=log(musq/mt**2)
         lt=log(musq/173.2d0**2)
         

         msq(0,0,1)=msq(0,0,1)+gg*(as/pi)**2*
     1   (11399/144+19*lt/8d0+133*pisq/12d0-pisq**2/80d0-165*zeta3/4d0)

         msq(0,0,2)=msq(0,0,2)+gg*(as/pi)**2*
     1  (((133 - 15*pisq)*log(1-z))/(1-z) 
     2   -(33*log(1-z)**2)/(1-z) + (72*log(1-z)**3)/(1-z)
     3   +(-101/3d0 + (11*pisq)/2d0 + (351*zeta3)/2d0)/(1-z))

         msq(0,0,3)=msq(0,0,3)+gg*(as/pi)**2*(
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
         
      endif
         
      return
      end


