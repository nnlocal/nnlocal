      subroutine reader_input(inputfile,workdir)
************************************************************************
*     Routine to read in the file input.DAT                            *
************************************************************************
      implicit none
      include 'constants.f'
      include 'order.f'
      include 'debug.f'
      include 'couple.f'
      include 'part.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'cutoff.f'
      include 'factiny.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'process.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'flags.f'
      include 'clustering.f'
      include 'gridinfo.f'
      include 'verbose.f'
      include 'limits.f'
      include 'werkdir.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'lhapdf.f'
      include 'pdlabel.f'
      include 'qcdcouple.f'
      include 'nlooprun.f'
      include 'initialscales.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'breit.f'
      include 'vdecayid.f'
      include 'runstring.f'
      include 'singular.f'
      include 'chkseed.f'
      include 'parallel.f'
      include 'finitemcorr.f'
      integer seed
      
      character*72 workdir,inputfile
      character*90 line
      logical spira,dryrun,makecuts
      integer nmin,nmax,ii
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,origij
      integer NPTYPE,NGROUP,NSET
      integer vonbcalls,ronbcalls
      double precision rtsmin,sqrts,factor
      double precision mbbmin,mbbmax,Mrmin,Mrmax
      double precision Rcut
      logical technicalincluded
      double precision ran2,randummy
      double precision cmass,bmass
      double precision alphas
      
      common/spira/spira
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin 

      common/nproc/nproc
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/callrats/vonbcalls,ronbcalls
      common/ranno/idum
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/qmass/cmass,bmass

      common/origij/origij

      save /ranno/
      
      verbose=.true.

      werkdir=workdir
c--- work out the name of the input file and open it

      write(6,*) '* Using input file named ',inputfile

      open(unit=20,file=inputfile,status='old',err=999)
      call checkversion(20,inputfile)

c--- read-in the user inputs

      writegnu = .true.
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- general options

c--- read in whole line for nproc
      read(20,99) line
      ii=index(line,'.')
      if (ii .gt. 0) then
        vdecayid=.true.
c------ special string present to specify V decays
        read(line(1:ii-1),*) nproc
        v34id=line(ii+1:ii+2)
        v56id=line(ii+3:ii+4)
c        write(6,*) 'special nproc=',nproc
c        write(6,*) 'v34id',v34id
c        write(6,*) 'v56id',v56id
c        stop
c------ normal case
      else
        vdecayid=.false.
        read(line,*) nproc
      endif
      if (verbose) call writeinput(6,' * ',' ','nproc')
      read(20,*) order
      if (verbose) call writeinput(6,' * ',' ','order')
      read(20,*) part
      if (verbose) call writeinput(6,' * ',' ','part')
      read(20,*) runstring
      if (verbose) call writeinput(6,' * ',' ','runstring')
      read(20,*) sqrts
      if (verbose) call writeinput(6,' * ',' ','sqrts')
      read(20,*) ih1
      if (verbose) call writeinput(6,' * ',' ','ih1')
      read(20,*) ih2
      if (verbose) call writeinput(6,' * ',' ','ih2')
      read(20,*) hmass
      if (verbose) call writeinput(6,' * ',' ','hmass')
      read(20,*) scale
      initscale=scale
      if (verbose) call writeinput(6,' * ',' ','scale')
      read(20,*) facscale
      initfacscale=facscale
      if (verbose) call writeinput(6,' * ',' ','facscale')
      read(20,*) dynstring 
      if (verbose) call writeinput(6,' * ',' ','dynamicscale')
      read(20,*) zerowidth
      if (verbose) call writeinput(6,' * ',' ','zerowidth')
      read(20,*) removebr
      if (verbose) call writeinput(6,' * ',' ','removebr')
      read(20,*) itmx1
      if (verbose) call writeinput(6,' * ',' ','itmx1')
      read(20,*) ncall1
      if (verbose) call writeinput(6,' * ',' ','ncall1')
      read(20,*) itmx2
      if (verbose) call writeinput(6,' * ',' ','itmx2')
      read(20,*) ncall2
      if (verbose) call writeinput(6,' * ',' ','ncall2')
      read(20,*) vonbcalls
      read(20,*) ronbcalls
      read(20,*) origij
      read(20,*) dryrun
      read(20,*) parallel
      if (parallel.eq.1) then
         call get_command_argument(2,seedchar)
         read(seedchar,'(I4)') seed
         origij=seed
         if (verbose) call writeinput(6,' * ',' ','ij')
         call get_command_argument(3,seedchar)
         read(seedchar,'(I4)') seed
         write(seedchar,'(i4.4)') seed
         call get_command_argument(4,stage)
         itmx1=1
c         itmx2=1
         if (stage.ne.'st2') itmx2=0
      else
         if (verbose) call writeinput(6,' * ',' ','ij')
      endif
      if (verbose) call writeinput(6,' * ',' ','dryrun')
      read(20,*) fmcorr
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- heavy quark masses 
      read(20,*) mt
      if (verbose) call writeinput(6,' * ',' ','top mass')
      read(20,*) mb
      if (verbose) call writeinput(6,' * ',' ','bottom mass')
      read(20,*) mc
      if (verbose) call writeinput(6,' * ',' ','charm mass')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- pdf options 
      read(20,*) PDFname
      if (verbose) call writeinput(6,' * ',' ','LHAPDF group')
      read(20,*) PDFmember
      if (verbose) call writeinput(6,' * ',' ','LHAPDF set')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- jets and cuts options 
      read(20,*) Mrmin
      wsqmin=Mrmin**2
      if (verbose) call writeinput(6,' * ',' ','m34min')
      read(20,*) Mrmax 
      if (Mrmax .gt. sqrts) Mrmax=sqrts ! physical cap on m34max
      wsqmax=Mrmax**2
      if (verbose) call writeinput(6,' * ',' ','m34max')
      read(20,*) mbbmin
      bbsqmin=mbbmin**2
      if (verbose) call writeinput(6,' * ',' ','m56min')
      read(20,*) mbbmax 
      if (mbbmax .gt. sqrts) Mbbmax=sqrts ! physical cap on m56max
      bbsqmax=mbbmax**2
      if (verbose) call writeinput(6,' * ',' ','m56max')
      read(20,*) inclusive
      if (verbose) call writeinput(6,' * ',' ','inclusive')
      read(20,*) algorithm
      if (verbose) call writeinput(6,' * ',' ','algorithm')
      read(20,*) ptjetmin
      if (verbose) call writeinput(6,' * ',' ','ptjetmin')
      read(20,*) etajetmin 
      if (verbose) call writeinput(6,' * ',' ','etajetmin')
      read(20,*) etajetmax 
      if (verbose) call writeinput(6,' * ',' ','etajetmax')
      read(20,*) Rcut
      if (verbose) call writeinput(6,' * ',' ','Rcut')
      read(20,*) makecuts
      if (verbose) call writeinput(6,' * ',' ','makecuts')
      read(20,*) leptpt
      if (verbose) call writeinput(6,' * ',' ','leptpt')
      read(20,*) leptrap
      if (verbose) call writeinput(6,' * ',' ','leptrap')
      read(20,*) leptveto1min,leptveto1max
      if (verbose) call writeinput(6,' * ',' ','leptveto')
      read(20,*) misspt
      if (verbose) call writeinput(6,' * ',' ','misspt')
      read(20,*) leptpt2
      if (verbose) call writeinput(6,' * ',' ','leptpt2')
      read(20,*) leptrap2
      if (verbose) call writeinput(6,' * ',' ','leptrap2')
      read(20,*) leptveto2min,leptveto2max
      if (verbose) call writeinput(6,' * ',' ','leptveto2')
      read(20,*) mtrans34cut
      if (verbose) call writeinput(6,' * ',' ','mtrans34cut')
      read(20,*) Rjlmin
      if (verbose) call writeinput(6,' * ',' ','Rjlmin')
      read(20,*) Rllmin
      if (verbose) call writeinput(6,' * ',' ','Rllmin')
      read(20,*) delyjjmin
      if (verbose) call writeinput(6,' * ',' ','delyjjmin')
      read(20,*) jetsopphem 
      if (verbose) call writeinput(6,' * ',' ','jetsopphem')
      read(20,*) lbjscheme 
      if (verbose) call writeinput(6,' * ',' ','lbjscheme')
      read(20,*) ptbjetmin
      if (verbose) call writeinput(6,' * ',' ','ptbjetmin')
      read(20,*) etabjetmax
      if (verbose) call writeinput(6,' * ',' ','etabjetmax')
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- grid information 
      read(20,*) readin
      read(20,*) writeout
      if (parallel.eq.1) then
         if (stage.eq.'xg1') then
            readin = .false.
         else
            readin = .true.
         endif
         if (stage.eq.'st2') then
            writeout = .false.
         else
            writeout = .true.
         endif
      endif
      if (verbose) call writeinput(6,' * ',' ','readin')
      if (verbose) call writeinput(6,' * ',' ','writeout')

      if (verbose) write(6,*)

c--- check if the contents of technical.DAT are included here
c--- (the default behaviour going forward)
      technicalincluded=.false.

      read(20,99,end=88) line
      read(20,99,end=88) line
c      write(6,*) line
      if (line(2:10) .ne. 'Technical') goto 88
      technicalincluded=.true.

   88 continue
      
c--- if the contents of technical.DAT are not included, close the input
c--- file and open technical.DAT instead; otherwise continue on
      if (technicalincluded .eqv. .false.) then
        close(20)
        open(unit=20,file='technical.DAT',status='old',err=999)
        call checkversion(20,'technical.DAT')
      endif

      if (verbose) write(6,*) '* [Technical parameters that'//
     .                        ' should not normally be changed]'
      if (verbose) write(6,*)

c---- read-in the technical parameters

      read(20,*) debug
      if (verbose) call writeinput(6,' * ',' ','debug')
      read(20,*) verbose
      if (verbose) call writeinput(6,' * ',' ','verbose')
      read(20,*) virtonly
      if (verbose) call writeinput(6,' * ',' ','virtonly')
      read(20,*) realonly
      if (verbose) call writeinput(6,' * ',' ','realonly')
      read(20,*) spira
      if (verbose) call writeinput(6,' * ',' ','spira')
      read(20,*) noglue
      if (verbose) call writeinput(6,' * ',' ','noglue')
      read(20,*) ggonly
      if (verbose) call writeinput(6,' * ',' ','ggonly')
      read(20,*) gqonly
      if (verbose) call writeinput(6,' * ',' ','gqonly')
      read(20,*) omitgg
      if (verbose) call writeinput(6,' * ',' ','omitgg')
      read(20,*) nmin
      if (verbose) call writeinput(6,' * ',' ','nmin')
      read(20,*) nmax
      if (verbose) call writeinput(6,' * ',' ','nmax')
      read(20,*) clustering
      if (verbose) call writeinput(6,' * ',' ','clustering')
      read(20,*) realwt
      if (verbose) call writeinput(6,' * ',' ','realwt')
      read(20,*) rtsmin
      if (verbose) call writeinput(6,' * ',' ','rtsmin')
      read(20,*) cutoff
      if (verbose) call writeinput(6,' * ',' ','cutoff')
      if (verbose) write(6,*)
      read(20,*) factiny
      read(20,*) check
      read(20,*) cseed
      read(20,*) limit
      read(20,*) ip
      read(20,*) jp
      read(20,*) rp
      read(20,*) sp

      close(unit=20)

      if ((etajetmin .lt. 0d0) .or. (etajetmax .lt. 0d0)) then
        write(6,*) 'etajetmin and etajetmax are absolute values,'
      write(6,*) ' please reset to a positive value.'
      stop
      endif
      
c---  create logical variable dynamicscale for use in other routines
      if (  (dynstring .eq. 'no') .or. (dynstring .eq. '.false.')
     & .or. (dynstring .eq. 'none') ) then 
         dynamicscale=.false. 
      else
         dynamicscale=.true. 
      endif

c--- print warning messages if some parton fluxes are not included      
      if (noglue) then
        write(6,*) 'WARNING: no gluon contribution included in PDF'
      write(6,*)
      endif
      if (ggonly) then
        write(6,*) 'WARNING: only gluon-gluon flux included'
      write(6,*)
      endif
      if (gqonly) then
        write(6,*) 'WARNING: only gluon-quark flux included'
      write(6,*)
      endif
      if (omitgg) then
        write(6,*) 'WARNING: no gluon-gluon contribution included'
      write(6,*)
      endif
      
c--- assign squared masses for b- and c-quarks
      if (abs(mb) .gt. 1d-8) then
        mbsq=mb**2
      else
        mbsq=4.75d0**2
      endif
      if (abs(mc) .gt. 1d-8) then
        mcsq=mc**2
      else
        mcsq=1.5d0**2
      endif
      
c--- set-up the variables for the process we wish to consider
      call chooser

c--- set-up the random number generator with a negative seed
      idum=-abs(origij)
      randummy=ran2()

c--- initialize masses for alpha_s routine
      cmass=dsqrt(mcsq)
      bmass=dsqrt(mbsq)

c--- check that we have a valid value of 'part'
      if ( (part .ne. 'born') .and. (part .ne. 'real') .and.
     .     (part .ne. 'virt') .and. (part .ne. 'tota') ) then
         write(6,*) 'part=',part,' is not a valid option'
         write(6,*) 'for this process number.'
         stop     
      endif      


c--- set up the default choices of static scale, if required
      if (scale .lt. 0d0) then
         if     (scale .eq. -2d0) then
            factor=0.25d0
         elseif (scale .eq. -3d0) then
            factor=0.5d0
         elseif (scale .eq. -4d0) then
            factor=0.75d0
         elseif (scale .eq. -5d0) then
            factor=1d0
         elseif (scale .eq. -6d0) then
            factor=2d0
         elseif (scale .eq. -7d0) then
            factor=4d0
         else
            factor=1d0
         endif        
         if (n2+n3 .ne. 0) then
            if (case(1:1) .eq. 'H') then
               scale=factor*hmass
            else     
               scale=factor*(dfloat(n2)*mass2+dfloat(n3)*mass3)/dfloat(n2+n3)
            endif
            as=alphas(scale,amz,nlooprun)
            ason2pi=as/twopi
            ason4pi=as/fourpi
            gsq=fourpi*as
            musq=scale**2
            write(6,*)
            write(6,*)'************* Strong coupling, alpha_s  ************'
            write(6,*)'*                                                  *'
            write(6,49)'alpha_s (scale)',gsq/fourpi
            write(6,49)'alpha_s (zmass)',amz
            write(6,50)' (using ',nlooprun,'-loop running of alpha_s)'  
            write(6,*)'****************************************************'
            write(6,*)
            write(6,*)'****************************************************'
            write(6,76) scale
            write(6,*)'****************************************************'
         else
            write(6,*) 'Invalid choice of renormalization scale!'
            stop
         endif
      endif
      if (facscale .lt. 0d0) then
         if     (facscale .eq. -2d0) then
            factor=0.25d0
         elseif (facscale .eq. -3d0) then
            factor=0.5d0
         elseif (facscale .eq. -4d0) then
            factor=0.75d0
         elseif (facscale .eq. -5d0) then
            factor=1d0
         elseif (facscale .eq. -6d0) then
            factor=2d0
         elseif (facscale .eq. -7d0) then
            factor=4d0
         else
            factor=1d0
         endif        
         if (n2+n3 .ne. 0) then
            if (case(1:1) .eq. 'H') then
               facscale=factor*hmass
            else     
               facscale=factor*
     &              (dfloat(n2)*mass2+dfloat(n3)*mass3)/dfloat(n2+n3)
            endif
            write(6,*)
            write(6,*)'****************************************************'
            write(6,77) facscale
            write(6,*)'****************************************************'
         else
            write(6,*) 'Invalid choice of factorization scale!'
            stop
         endif
      endif
      
      return

   49 format(' *  ',a20,f12.8,16x,'*')
   50 format(' *  ',6x,a8,i1,a25,8x,'*')
   76 format(' *      Renormalization scale =',f7.2,'              *')
   77 format(' *        Factorization scale =',f7.2,'              *')
   99 format(a90)

  999 continue
      write(6,*) 'Problem reading ',inputfile
      write(6,*)
      write(6,*) 'Refer to documentation for the format of input.DAT'
      write(6,*)
      stop

      end
      
