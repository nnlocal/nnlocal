      integer maxd,ndmax,ndmax1
!----maxd=The maximum possible number of dipoles
!----ndmax=The maximum number of dipoles for the problem at hand
      parameter (maxd=100)
      double precision ptilde(0:maxd,mxpart,4)
      double precision ptildejet(0:maxd,mxpart,4)
      common/ptildes/ptilde,ptildejet,ndmax,ndmax1
