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
  
c      virtgg=-3d0*epinv**2*xn*logg 
      virtgg=-3d0*(epinv*epinv2)*xn*logg 
     .     +epinv*xn*logg*(lns+lnt+lnu-3d0*lnm )
     .     +xn*logg
     .     *(2d0*(Li2t+Li2u+Li2s)
     .     +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
     .     +0.5d0*(lns**2-lnt**2-lnu**2)-1.5d0*lnm**2
     .     +2d0*(lnu*ln2u+lnt*ln2t)+4d0/3d0*pisq)
     .     +V*xn*(xn-xlf)/3d0*mhsq*(1d0+mhsq/s+mhsq/t+mhsq/u)
     .     +subuv*logg

c---  DEBUG: return ep^(-2) part
c      virtgg = -3d0*xn*logg 
c---  END DEBUG
c---  DEBUG: return ep^(-1) part
c      subuv=-b0
c      Delta=0d0
c      subuv=3d0*subuv+Delta
c      virtgg=xn*logg*(lns+lnt+lnu-3d0*lnm )
c     .     +subuv*logg
c---  END DEBUG
      

c$$$      virtqa=+(-2d0*xn+1d0/xn)*loqa*epinv**2
c$$$     . -2d0/3d0*xlf*epinv*loqa
c$$$     . +epinv*xn*loqa*(13d0/6d0-2d0*lnm+lnt+lnu)
c$$$     . +epinv/xn*loqa*(1.5d0-lns+lnm)
c$$$     . +loqa*xlf*(-10d0/9d0+2d0/3d0*lns-2d0/3d0*lnm)
c$$$     . +xn*loqa* (40d0/9d0+Li2t+Li2u+2d0*Li2s-13d0/6d0*(lns-lnm)
c$$$     . +(lnm-lns)*(lnt+lnu)+lns**2-lnm**2-0.5d0*lnt**2-0.5d0*lnu**2
c$$$     . +lnt*ln2t+lnu*ln2u)
c$$$     . +loqa/xn
c$$$     . *(4d0-Li2t-Li2u-1.5d0*(lns-lnm)+0.5d0*(lns-lnm)**2
c$$$     . +lnt*lnu-lnt*ln2t-lnu*ln2u)
c$$$     . -4d0/3d0*pi**2/xn*loqa
c$$$     . -0.25d0*(xn**3-1d0/xn)*(t+u)
c$$$     . +subuv*loqa
c$$$
c$$$      virtaq=(-2d0*xn+1d0/xn)*loaq*epinv**2
c$$$     . -2d0/3d0*xlf*epinv*loaq
c$$$     . +epinv*xn*loaq*(13d0/6d0-2d0*lnm+lnu+lnt)
c$$$     . +epinv/xn*loaq*(1.5d0-lns+lnm)
c$$$     . +loaq*xlf*(-10d0/9d0+2d0/3d0*lns-2d0/3d0*lnm)
c$$$     . +xn*loaq* (40d0/9d0+Li2u+Li2t+2d0*Li2s-13d0/6d0*(lns-lnm)
c$$$     . +(lnm-lns)*(lnu+lnt)+lns**2-lnm**2-0.5d0*lnu**2-0.5d0*lnt**2
c$$$     . +lnu*ln2u+lnt*ln2t)
c$$$     . +loaq/xn
c$$$     . *(4d0-Li2u-Li2t-1.5d0*(lns-lnm)+0.5d0*(lns-lnm)**2
c$$$     . +lnu*lnt-lnu*ln2u-lnt*ln2t)
c$$$     . -4d0/3d0*pi**2/xn*loaq
c$$$     . -0.25d0*(xn**3-1d0/xn)*(u+t)
c$$$     . +subuv*loaq
c$$$ 
c$$$
c$$$      virtgq=(-2d0*xn+1d0/xn)*epinv**2*logq
c$$$     . -2d0/3d0*xlf*epinv*logq
c$$$     . +epinv*xn*logq*(13d0/6d0+lns-2d0*lnm+lnt)
c$$$     . +epinv/xn*logq*(3d0/2d0+lnm-lnu)
c$$$     . +logq*xlf*(-10d0/9d0-2d0/3d0*lnm+2d0/3d0*lnu)
c$$$     . +xn*logq*(40d0/9d0+Li2t+2d0*Li2u+Li2s
c$$$     . +lns*lnm-lns*lnu-13d0/6d0*(lnu-lnm)
c$$$     . +lnm*lnt-lnm**2-lnt*lnu-0.5d0*lnt**2
c$$$     . +2d0*lnu*ln2u+lnt*ln2t)
c$$$     . +logq/xn*(4d0-Li2t-Li2s+lns*lnt+0.5d0*lnu**2-0.5d0*lns**2
c$$$     . -lnm*lnu+0.5d0*lnm**2-lnt*ln2t-1.5d0*(lnu-lnm))
c$$$     . +4d0/3d0*pi**2*xn*logq
c$$$     . +0.25d0*(xn**3-1d0/xn)*(t+s)
c$$$     . +subuv*logq
c$$$
c$$$      virtqg=(-2d0*xn+1d0/xn)*epinv**2*loqg
c$$$     . -2d0/3d0*xlf*epinv*loqg
c$$$     . +epinv*xn*loqg*(13d0/6d0+lns-2d0*lnm+lnu)
c$$$     . +epinv/xn*loqg*(3d0/2d0+lnm-lnt)
c$$$     . +loqg*xlf*(-10d0/9d0-2d0/3d0*lnm+2d0/3d0*lnt)
c$$$     . +xn*loqg*(40d0/9d0+Li2u+2d0*Li2t+Li2s
c$$$     . +lns*lnm-lns*lnt-13d0/6d0*(lnt-lnm)
c$$$     . +lnm*lnu-lnm**2-lnu*lnt-0.5d0*lnu**2
c$$$     . +2d0*lnt*ln2t+lnu*ln2u)
c$$$     . +loqg/xn*(4d0-Li2u-Li2s+lns*lnu+0.5d0*lnt**2-0.5d0*lns**2
c$$$     . -lnm*lnt+0.5d0*lnm**2-lnu*ln2u-1.5d0*(lnt-lnm))
c$$$     . +4d0/3d0*pi**2*xn*loqg
c$$$     . +0.25d0*(xn**3-1d0/xn)*(u+s)
c$$$     . +subuv*loqg

      return
      end
