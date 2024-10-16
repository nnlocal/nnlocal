      !----------------------------------------------------------------------
      !      This replaces the original includedipole
      !      It calls the original and then the user one
      logical function includedipole(nd,ptrans)
      implicit none
      include 'constants.f'
      double precision ptrans(mxpart,4)
      integer nd
      logical nnlocal_includedipole
      
      ! it looks like all (nd .ne. 0) automatically have their momenta
      ! stored in the ptilde common; do this also for nd=0
      if (nd .eq. 0) call storeptilde(nd,ptrans)

      ! first call the original NNLOCAL includedipole
      includedipole = nnlocal_includedipole(nd,ptrans)

      end 

      logical function nnlocal_includedipole(nd,ptrans) 
     &                 result(nnlocalincdipole)
c--- This function returns TRUE if the specified point ptrans,
c--- corresponding to dipole nd (nd=0 => real radiation),
c--- should be included 
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'process.f'
      include 'reweight.f'

      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut
      integer j,nd,nqcdjets,nqcdstart,notag,isub,nproc
      logical gencuts,failedgencuts,makecuts
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag
      common/nproc/nproc
c---- set default reweight to 1 (hence no reweighting)      
      reweight = 1.0d0

c--- default: include this contribution
      nnlocalincdipole=.true.
c---  isub=2 for double unresolved dipole subtractions
c---  isub=1 for single unresolved dipole subtractions
c---  isub=0 for real radiation
      if     (nd .gt. ndmax1) then
        isub=2
      elseif (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

      call genclust2(ptrans,rcut,pjet,isub)


c--- perform mass cuts
      call masscuts(pjet,*999)
          
c--- fill ptilde array as persistent storage for the jet momenta
      ptildejet(nd,1:npart+2,1:4)=pjet(1:npart+2,1:4)
          
c--- if the number of jets is not correct, then do not include dipole
      if ((clustering .and. (jets .ne. nqcdjets-notag)
     &       .and. (inclusive .eqv. .false.)) .or.
     &    (clustering .and. (jets .lt. nqcdjets-notag)
     &       .and. (inclusive .eqv. .true.))) then
          nnlocalincdipole=.false.
          return
      else
c--- otherwise check the lepton cuts, if necessary
         if (makecuts) then
            failedgencuts=gencuts(pjet,jets,nd)        
            if (failedgencuts) then
               nnlocalincdipole=.false.
               return
            endif
         endif
      endif
            
     
   99 continue
      
      return

  999 continue
      nnlocalincdipole=.false.
      return    
      
      end
            

