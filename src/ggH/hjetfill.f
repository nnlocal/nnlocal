      subroutine hjetfill(s,t,u,virtgg,virtqa,virtaq,virtqg,virtgq)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'b0.f'
      double precision virtgg,virtqa,virtaq,virtqg,virtgq,
     . logg,loqa,loaq,loqg,logq,ddilog,Li2s,Li2t,Li2u,
     . lnm,lns,lnt,lnu,ln2t,ln2u,mhsq,s,t,u,xlf,subuv,Delta
      
      mhsq=s+t+u

c--- we match results of other codes with nf=0
c     xlf=dfloat(nf)
      xlf=0d0
      Li2t=ddilog(t/mhsq)
      Li2u=ddilog(u/mhsq)
      Li2s=ddilog((s-mhsq)/s)
      lns=dlog(s/mhsq)
      lnt=dlog(-t/mhsq)
      lnu=dlog(-u/mhsq)
      lnm=dlog(musq/mhsq)
      ln2t=dlog((mhsq-t)/mhsq)
      ln2u=dlog((mhsq-u)/mhsq)

      logg=+V*xn*(mhsq**4+s**4+t**4+u**4)/(s*t*u)
      loqa=xn*cf/s*(t**2+u**2)
      loaq=loqa
      logq=-xn*cf/u*(s**2+t**2)
      loqg=-xn*cf/t*(s**2+u**2)

c--- UV counterterm in MS bar scheme. 
      subuv=-epinv*b0
      Delta=11d0
c--- See C.R.Schmidt, PLB (413) 391, eq. (16),(17)
c--- Factor of ason2pi included in gg_hg_v.f
C--- Three powers of as in Born --> 3      
      subuv=3d0*subuv+Delta
  
      virtgg=-3d0*(epinv*epinv2)*xn*logg 
     .     +epinv*xn*logg*(lns+lnt+lnu-3d0*lnm )
     .     +xn*logg
     .     *(2d0*(Li2t+Li2u+Li2s)
     .     +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
     .     +0.5d0*(lns**2-lnt**2-lnu**2)-1.5d0*lnm**2
     .     +2d0*(lnu*ln2u+lnt*ln2t)+4d0/3d0*pisq)
     .     +V*xn*(xn-xlf)/3d0*mhsq*(1d0+mhsq/s+mhsq/t+mhsq/u)
     .     +subuv*logg

      return
      end
