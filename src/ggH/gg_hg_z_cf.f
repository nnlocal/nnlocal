      subroutine gg_hg_z_cf(p,z,t)
************************************************************************
*     Authors: G. Somogyi and F. Tramontano                            *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is
      double precision z,t,xl12,p(mxpart,4),dot
      double precision x5,Y12,Y15,Y25
      double precision f_g,i_g,f_s

      xl12=dlog(two*dot(p,1,2)/musq)
      x5=-(dot(p,1,5)+dot(p,2,5))/dot(p,1,2)
      Y12=1d0
      Y15=dot(p,1,5)/(dot(p,1,5)+dot(p,2,5))
      Y25=dot(p,2,5)/(dot(p,1,5)+dot(p,2,5))
c      Y25 = 1d0-Y15
      
c--- sum over all terms
      do is=1,4

c---  (g,g)

         T1(g,g,g,is)=1d0*ason2pi*ca*(
     .        i_g(z,t,xl12,is)
     .        + f_g(z,xl12,x5,is)/2d0
     .        -half*(f_s(z,xl12,Y12,is)
     .        +f_s(z,xl12,Y15,is)+f_s(z,xl12,Y25,is))
     .        )

         T2(g,g,g,is)=1d0*ason2pi*ca*(
     .        i_g(t,z,xl12,is)
     .        + f_g(t,xl12,x5,is)/2d0
     .        -half*(f_s(t,xl12,Y12,is)
     .        +f_s(t,xl12,Y15,is)+f_s(t,xl12,Y25,is))
     .        )
      
      enddo
      
      return
      end
