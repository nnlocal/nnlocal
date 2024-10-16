      subroutine gen_born_ps(r,p,pswt,*)      
c--- calls appropriate subroutine to choose a 
c--- leading-order phase space point for the process
c--- specified by the variable 'case'   
c---
c---    input: vector of random numbers, r
c---    output: momenta p and phase-space weight pswt
c---    note: common block 'npart' is also filled here
c---    
c---    alternate return if generated point should be discarded
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'process.f'
      include 'zerowidth.f'
      include 'singular.f'
      double precision r(mxdim),p(mxpart,4),pswt,
     & ptmp,m3,m4,m5,sqrts
      common/energy/sqrts
      
c--- processes that use "gen2"     
      if (case .eq. 'H_gaga') then
        npart=2
        call gen2(r,p,pswt,*999)
      endif

      return

c--- alternate return      
  999 continue   
      return 1
      
      end
      
