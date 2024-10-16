***********************************************************************
************************ FINAL-STATE COLLINEAR ************************
***********************************************************************

******************* Final state gluon splitting  **********************
      double precision function f_g(x,L,xir,vorz)
      implicit none
      integer vorz
      double precision x,L,xir,omx,lomx,lx,lxir
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for
c--- final state gluon splitting
c---  vorz=1: (Id - Ip[xir]) [double delta],
c---  vorz=3: (Ip[x xir] + Ir) [same z's]

      f_g = 0d0

      omx=one-x
      lomx=dlog(omx)
      lxir=dlog(xir)
    
      if (vorz .eq. 1) then
         f_g=two*epinv*(epinv2-L)+L**2
     .        +(two*11d0/6d0 - 4d0*lxir)*(epinv-L)
     .        +67d0/9d0-pisq-11d0/3d0*lxir+2d0*lxir**2
     .        +11d0/(3d0*omx) + 4d0*lomx/omx - 4d0*lxir/omx

c---  DEBUG: return ep^(-2) part
c         f_g=two         
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         f_g=two*(-L)
c     .        +(two*11d0/6d0 - 4d0*lxir)
c---  END DEBUG

         f_g=f_g/2d0
c---  Note: g -> qq not added yet
         return
      endif

      if (vorz .eq. 3) then
         lx=dlog(x)
         
         f_g=
     .        (x**2*(-(x*xir*(42d0 + 42d0*x*(-2d0 + xir) + x**2*(42d0
     .        + xir*(-42d0 + 11d0*xir)))) - 12d0*(2d0 + x*(-2d0
     .        + xir))**3*lomx + 12d0*(2d0 + x*(-2d0 + xir))**2*(1d0
     .        + omx + x*(-1d0 + xir))*dlog(omx + x*xir)))/
     .        (3d0*(2d0 + x*(-2d0 + xir))**2*(omx + omx*x*(-1d0 + xir)))        

c     .    ((-12d0*x**2*(2d0*omx + xir)*lomx)/(omx + xir) + 12d0*lx
c     .        +(x**2*(-(xir*(42d0*omx**2 + 42d0*omx*xir + 11d0*xir**2))
c     .        +12d0*(2d0*omx + xir)**3*dlog(omx + xir)))/((omx + xir)
c     .        *(2d0*omx + xir)**2))/(3d0*omx)
         f_g=f_g/2d0

c---  DEBUG: return ep^(-2) part
c         f_g=0d0    
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         f_g=0d0    
c---  END DEBUG

         return
      endif
c---  Note: g -> qq not added yet
      return

      end


***********************************************************************
*********************** INITIAL-STATE COLLINEAR ***********************
***********************************************************************

****************** Initial state gluon splitting  *********************
      double precision function i_g(xa,xb,L,vorz)
      implicit none
      integer vorz
      double precision xa,xb,L,omxa,omxb,lxa,lxb,lomxa,lomxb
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for
c--- final state gluon splitting
c---  vorz=1: (Idd - Ipd) [double delta] 
c---  vorz=2: (Ipd + Ird - Irp) [one delta, one z]
c---  vorz=4: (Irp + Irr) [different z's]

      i_g = 0d0

      omxa=1d0-xa
      lomxa=dlog(omxa)

      if (vorz .eq. 1) then
         i_g=epinv*(epinv2-L)+0.5d0*L**2
     .        + (epinv-L)*2d0/omxa
     .        + pisq/6d0-4d0*lomxa/omxa

c---  DEBUG: return ep^(-2) part
c         i_g=1d0    
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         i_g=(-L)
c     .        + 2d0/omxa
c---  END DEBUG
         
         return
      endif

      if (vorz .eq. 2) then
         i_g= (epinv-L)*(-2d0*(omxa + xa**2)**2)/(omxa*xa)
     .        + (-4d0 + 2d0*xb - 2d0*xa*(2d0 - omxa*xa)
     .        *(-2d0 + xa + xb))/(xa*(-1d0 + xb)*(-2d0 + xa + xb))
     .        + (2d0*(1d0 - omxa*xa)**2*dlog(2d0))/(omxa*xa)
     .        + 2d0*(1d0/xa + xa + (2d0*xa)/omxa - xa**2)*lomxa
     .        - (2d0*dlog(2d0 - xa))/omxa
     .        - (2d0*(1d0 - omxa*xa)**2*dlog(1d0 + xa))/(omxa*xa)

c---  DEBUG: return ep^(-2) part
c         i_g=0d0    
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         i_g= (-2d0*(omxa + xa**2)**2)/(omxa*xa)
c---  END DEBUG

         return
      endif

      if (vorz .eq. 4) then
         omxb=1d0-xb

         i_g=(-4d0*(1d0 + xa*xb)*(1d0 + xa*xb*(-1d0 + xa*xb))**2)
     .        /(omxb*xa*(1d0 + xb)*(xa + xb)*(-1d0 + xa*xb))

c---  DEBUG: return ep^(-2) part
c         i_g=0d0    
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         i_g=0d0    
c---  END DEBUG
         
         return
      endif

      end

***********************************************************************
**************************** SOFT RADIATION ***************************
***********************************************************************

************************* Soft type terms  ****************************
      double precision function f_s(x,L,Y,vorz)
      implicit none
      integer vorz
      double precision x,L,Y,omx,lY,ddilog
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for
c--- final state gluon splitting
c---  vorz=1: (Id - Ip) [double delta]
c---  vorz=3: (Ip + Ir) [same z's]

      f_s = 0d0
      
      lY=dlog(Y)
      omx=one-x
      
      if (vorz .eq. 1) then
         f_s=(epinv-L)*lY
     .        +(2d0*lY)/omx - 2d0*dlog(2d0)*lY - lY**2/2d0
     .        - ddilog(1d0 - Y)

c---  DEBUG: return ep^(-2) part
c         f_s=0d0    
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         f_s=lY
c---  END DEBUG
         
         return
      endif

      if (vorz .eq. 3) then
         f_s=(-2d0/omx + 4d0*x + 2d0/(1d0 + x))*lY

c---  DEBUG: return ep^(-2) part
c         f_s=0d0    
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c         f_s=0d0    
c---  END DEBUG
         
         return
      endif

      end
