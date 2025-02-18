      subroutine chooser
c---- Note added 4/21/03
c---- plabel set to 'ig' (for 'ignore') means that this
c---- particle should not be subject to any cuts, so that the
c---- total cross-section comes out correctly when the BR is removed
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'vegas_common.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'nwz.f'
      include 'process.f'
      include 'flags.f'
      include 'nflav.f'
      include 'nodecay.f'
      include 'scale.f'
      include 'facscale.f'
      include 'nlooprun.f'
      include 'b0.f'
      include 'colstruc.f'
      include 'clustering.f'
      include 'plabel.f'
      include 'couple.f'
      include 'part.f'
      include 'hdecaymode.f'
      include 'breit.f'
      include 'verbose.f'
      include 'runstring.f'
      include 'vdecayid.f'
      include 'nuflav.f'
      include 'startstop.f'
      include 'order.f'
      double precision wwbr,zzbr,tautaubr,gamgambr,zgambr,Rcut,Rbbmin,
     . alphas,cmass,bmass
      double precision br,BrnRat,brwen,brzee,brznn,brtau,brtop,brcharm
      integer nproc,mproc,j,nqcdjets,nqcdstart,isub,notag,ilomomenta
      character*100 pname
      character*1 ord
      character*72 string
      double precision f0q,f2q,f4q
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/bitflags/f0q,f2q,f4q
      common/Rbbmin/Rbbmin
      common/Rcut/Rcut
      common/nproc/nproc
      common/BrnRat/BrnRat
      common/nqcdjets/nqcdjets,nqcdstart
      common/isub/isub
      common/notag/notag
      common/ilomomenta/ilomomenta
      common/qmass/cmass,bmass
      data hdecaymode/'xxxx'/
      do j=1,mxpart
      plabel(j)=''
      enddo

      string='process.DAT' 
      open(unit=21,file=string,status='old',err=43)
      call checkversion(21,string)
      
      if (verbose) write(6,*) 'Chooser:process chosen by nproc=',nproc

      do j=1,600
      read(21,*,err=44) mproc,pname,ord
      
      if (nproc .lt. 0) then 
      write(6,*) mproc,pname 
      endif

      if (mproc .eq. nproc) go to 42
      if (pname .eq. 'EOF') go to 44
      enddo
      goto 44

 42   continue
      if (verbose) then
      write(6,*)
      write(6,*) '*************************** f(p1)+f(p2) --> *****'//
     . '*************************************'
      write(6,*) '* ',pname(19:100),' *'
      write(6,*) '*************************************************'//
     . '*************************************'
      write(6,*)
      endif

      close(unit=21)

c--- check no. of momenta appearing in LO process, fill ilomomenta common block
      if (index(pname,'p3') .gt. 0)  ilomomenta=3
      if (index(pname,'p4') .gt. 0)  ilomomenta=4
      if (index(pname,'p5') .gt. 0)  ilomomenta=5
      if (index(pname,'p6') .gt. 0)  ilomomenta=6
      if (index(pname,'p7') .gt. 0)  ilomomenta=7
      if (index(pname,'p8') .gt. 0)  ilomomenta=8
      if (index(pname,'p9') .gt. 0)  ilomomenta=9
      if (index(pname,'p10') .gt. 0) ilomomenta=10
      if (index(pname,'p11') .gt. 0) ilomomenta=11
      if (index(pname,'p12') .gt. 0) ilomomenta=12
      
      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio
      BrnRat=1d0

c--- set up most parameters
      call coupling

      notag=0
      nqcdjets=0
      isub=0
      nodecay=.false.
      nfonly=.false.
      caonly=.false.

c-- Rbbmin is an additional variable, added so that the separation
c-- between two b jets can be controlled separately from the Delta_R
c-- cut between any other types of jet
c-- Default behaviour: the same value as for the other jets
      Rbbmin=Rcut
      
c-----------------------------------------------------------------------

      if ((nproc .ge. 700) .and. (nproc .le. 711)) then
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(5)='pp'
        plabel(6)='pp'        
        ndim=4
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth

        call checkminzmass(0)
        
        if     (nproc .eq. 711) then
c--  711 '  f(p1)+f(p2) -->H TOTAL'
           case='H_tota'
           ndim=2
        endif
        
        if     (nproc .eq. 710) then
c--  710 '  f(p1)+f(p2) -->H(-->gamma(p3)+gamma(p4))'
            case='H_gaga'
            plabel(3)='ga'
            plabel(4)='ga'  
            hdecaymode='gaga'

            runstart = 0
            if (abs(order).eq.1) runstop = 1
                        
            if (removebr) then
              BrnRat=gamgambr
              plabel(3)='ig'
              plabel(4)='ig'
            endif
        endif

        if     (nproc .eq. 709) then
c--  709 '  f(p1)+f(p2) -->H(-->gamma(p3)+gamma(p4))'
            case='H_gaga'
            plabel(3)='ga'
            plabel(4)='ga'  
            hdecaymode='gaga'
            nqcdjets=1

            runstart = 1
            
            if (removebr) then
              BrnRat=gamgambr
              plabel(3)='ig'
              plabel(4)='ig'
            endif
        endif
        
      else 
        call nprocinvalid()
      endif

c--- set notag (may be modified by user, with care!)
      call setnotag()

c--- set up alpha-s again (in case nflav was changed)
      call coupling2

c--- remove 2 dimensions from integration if decay is not included
      if (nodecay) ndim=ndim-2

c--- report on the removed BR, if necessary
      if (removebr) then
        write(6,*)'****************************************************'
        write(6,98) BrnRat
        write(6,*)'****************************************************'
      endif


c--- initialize arrays that are used in is_functions
      call init_is_functions()

c--- fill up CKM matrix
      call ckmfill(nwz)

      return

 43   write(6,*) 'problems opening process.DAT'
      stop

 44   write(6,*) 'Unimplemented process number, nproc = ',nproc, 
     . ' exiting...'
      stop
 
 98   format(' *             Brn.Rat. removed = ',  f11.7, '       *')
     
      end

      subroutine nprocinvalid()
      implicit none
      integer nproc
      common/nproc/nproc

      write(6,*) 'chooser: Unimplemented case'
      write(6,*) 'nproc=',nproc      
      stop
      
      return 
      end
      
      subroutine checkminzmass(i)
c--- Checks that the minimum invariant mass specified in the options
c--- file is not zero for boson 34 (i=1) or boson 56 (i=2)
      implicit none
      include 'limits.f'
      include 'zerowidth.f'
      include 'masses.f'
      integer i

c--- if generating exactly on-shell, there's nothing to worry about
      if (zerowidth) return
      
      if ((i .eq. 0) .and. (wsqmin .ge. hmass)) then
        write(6,*)
        write(6,*) 'Please set m34min lower then the Higgs mass'
        stop
      endif

      if ((i .eq. 0) .and. (wsqmax .le. hmass)) then
        write(6,*)
        write(6,*) 'Please set m34max greather then the Higgs mass'
        stop
      endif
      
      return
      end
      
      
