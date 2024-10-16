      subroutine gg_hg_z_wrap(nd,p,z,t)
************************************************************************
*     Authors: G. Somogyi and F. Tramontano                            *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is, nd
      double precision z,t,xl12,p(mxpart,4),dot
      double precision x5,Y12,Y15,Y25
      double precision f_g,i_g,f_s


      xl12=dlog(two*dot(p,1,2)/musq)
      if (nd.eq.1) then
         x5=-(dot(p,1,5)+dot(p,2,5))/dot(p,1,2)
         Y12=1d0
         Y15=dot(p,1,5)/(dot(p,1,5)+dot(p,2,5))
         Y25=1d0
      elseif (nd.eq.2) then
         x5=-(dot(p,1,5)+dot(p,2,5))/dot(p,1,2)
         Y12=1d0
         Y15=1d0
         Y25=dot(p,2,5)/(dot(p,1,5)+dot(p,2,5))
      elseif (nd.eq.3) then
         x5=-(dot(p,1,5)+dot(p,2,5))/dot(p,1,2)
         Y12=1d0
         Y15=dot(p,1,5)/(dot(p,1,5)+dot(p,2,5))
         Y25=dot(p,2,5)/(dot(p,1,5)+dot(p,2,5))
      endif
         
c---  sum over all terms
      do is=1,4
            
c---  (g,g)
            
         Tp1(nd,g,g,g,is,0)=1d0*ason2pi*ca*(
     .        i_g(z,t,xl12,is)
     .        + f_g(z,xl12,x5,is)/2d0
     .        -half*(f_s(z,xl12,Y12,is)
     .        +f_s(z,xl12,Y15,is)+f_s(z,xl12,Y25,is))
     .        )
         
         Tp2(nd,g,g,g,is,0)=1d0*ason2pi*ca*(
     .        i_g(t,z,xl12,is)
     .        + f_g(t,xl12,x5,is)/2d0
     .        -half*(f_s(t,xl12,Y12,is)
     .        +f_s(t,xl12,Y15,is)+f_s(t,xl12,Y25,is))
     .        )

      enddo

      if (nd.eq.3) then

c---  sum over all terms
         do is=1,4
            
            x5=-(dot(p,1,5)+dot(p,2,5))/dot(p,1,2)
            Y12=1d0
            Y15=dot(p,1,5)/(dot(p,1,5)+dot(p,2,5))
            Y25=1d0

            Tp1(nd,g,g,g,is,1)=1d0*ason2pi*ca*(
     .           i_g(z,t,xl12,is)
     .           + f_g(z,xl12,x5,is)/2d0
     .           -half*(f_s(z,xl12,Y12,is)
     .           +f_s(z,xl12,Y15,is)+f_s(z,xl12,Y25,is))
     .           )
            
            Tp2(nd,g,g,g,is,1)=1d0*ason2pi*ca*(
     .           i_g(t,z,xl12,is)
     .           + f_g(t,xl12,x5,is)/2d0
     .           -half*(f_s(t,xl12,Y12,is)
     .           +f_s(t,xl12,Y15,is)+f_s(t,xl12,Y25,is))
     .           )

            x5=-(dot(p,1,5)+dot(p,2,5))/dot(p,1,2)
            Y12=1d0
            Y15=1d0
            Y25=dot(p,2,5)/(dot(p,1,5)+dot(p,2,5))

            Tp1(nd,g,g,g,is,2)=1d0*ason2pi*ca*(
     .           i_g(z,t,xl12,is)
     .           + f_g(z,xl12,x5,is)/2d0
     .           -half*(f_s(z,xl12,Y12,is)
     .           +f_s(z,xl12,Y15,is)+f_s(z,xl12,Y25,is))
     .           )
            
            Tp2(nd,g,g,g,is,2)=1d0*ason2pi*ca*(
     .           i_g(t,z,xl12,is)
     .           + f_g(t,xl12,x5,is)/2d0
     .           -half*(f_s(t,xl12,Y12,is)
     .           +f_s(t,xl12,Y15,is)+f_s(t,xl12,Y25,is))
     .           )

         enddo
      endif
              
      return
      end
