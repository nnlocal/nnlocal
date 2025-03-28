      subroutine scaleset(rscalestart,fscalestart,p)
c--- wrapper routine to set a dynamic scale; please refer to individual
c--- routines for exact definitions of the scales.
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'nlooprun.f'
      include 'qcdcouple.f'
      include 'couple.f'
      double precision rscalestart,fscalestart,p(mxpart,4),mu0,
     & alphas
      logical first
      data first/.true./  
      save first
      
      if     (dynstring .eq. 'm(34)') then
        call scaleset_m34(p,mu0)
      elseif (dynstring .eq. 'sqrt(M^2+pt34^2)') then
        call scaleset_Msqpt34sq(p,mu0)
      elseif (dynstring .eq. 'HT') then
        call scaleset_HT(p,mu0)
      else
        write(6,*) 'Dynamic scale choice not recognized'
        write(6,*) '   dynamicscale = ',dynstring
        stop
      endif
      
      scale=rscalestart*mu0
      facscale=fscalestart*mu0
          
      if (first) then
        write(6,*)
        write(6,*)'************** Dynamic scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*                 RENORMALIZATION                  *'
        write(6,45) ' mu_ren  =',rscalestart,dynstring
        write(6,*)'*                                                  *'
        write(6,*)'*                  FACTORIZATION                   *'
        write(6,45) ' mu_fac  =',fscalestart,dynstring
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.      
      endif
  
c--- catch absurdly large and small scales      
      if  (scale .gt. 10000d0) scale=10000d0
      if  (facscale .gt. 10000d0) facscale=10000d0
      if  (scale .lt. 1d0) scale=1d0
      if  (facscale .lt. 1d0) facscale=1d0

c--- run alpha_s
      as=alphas(scale,amz,nlooprun)
      
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2

      return

 45   format(1x,'* ',a15,f6.2,' x ',a24,' *')

      end
      
