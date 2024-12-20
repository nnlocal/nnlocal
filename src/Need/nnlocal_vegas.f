      subroutine nnlocal_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)
************************************************************************
*                                                                      *
*  This routine should perform the sweeps of vegasnr                   *
*                                                                      *
*    Input parameters:                                                 *
*       myinit  :  the vegasnr routine entry point                     *
*       myitmx  :  the number of vegasnr sweeps                        *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error                                   *
*                                                                      *
*    born -> Born kinematic + c.t. at the required order               *  
*    virt -> Real-Virt kinematic + c.t. at the required order          *
*    real -> Real kinematic + c.t. at the required order               *  
*                                                                      *
*    tota -> born + virt + real kin.                                   *
*                                                                      *
************************************************************************

      implicit none
      include 'gridinfo.f'
      include 'realwt.f'
      include 'order.f'
      include 'scale.f'
      include 'facscale.f'
      include 'vegas_common.f'
      include 'PDFerrors.f'
      include 'reset.f'
      include 'masses.f'
      include 'process.f'
      include 'part.f'
      include 'parallel.f'
      include 'startstop.f'
      logical dryrun
      common/dryrun/dryrun
      integer myitmx,myncall,myinit,i,j,k,nproc,mynproc,myncall_save
      integer vonbcalls,ronbcalls
      logical mybin,bin
      double precision sigb,sigv,sigr,sdb,sdv,sdr,chi,sigdk,sddk,chidk,
     & xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsigf,sumsd,sumsdr,sumsdf,
     & sigips(4),sdips(4),xcallwt
      character*4 mypart
      character*3 getstr,psgen
      common/callrats/vonbcalls,ronbcalls
      common/nproc/nproc
      common/mypart/mypart
      common/bin/bin
      double precision bornint,virtint,realint
      double precision region(2*mxdim),lord_bypart(-1:1,-1:1)
      logical first,myreadin
      common/bypart/lord_bypart
      external bornint,virtint,realint
      data first/.true./
      save first,sigips,sdips
           
c--- Initialize all integration results to zero, so that the
c--- total of virt and real may be combined at the end for 'tota'
      sigb=0d0
      sigv=0d0
      sigr=0d0
      sdb=0d0
      sdv=0d0
      sdr=0d0
      
      
      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=0d0
      enddo
      enddo
      if (PDFerrors) then
        do i=0,maxPDFsets
          PDFxsec(i)=0d0
        enddo
      endif


            
c--- Controls behaviour of gen_njets: need to reset phase-space
c--- boundaries when going from virt to real (using tota)
c--- need to reset scale also, for special scalestart values
      reset=.false.
      scalereset=.false.

c--- Put the vegasnr parameters in the common block
      itmx=myitmx
      ncall=myncall
      bin=mybin
      
c--- Store value of part in mypart, which will be retained;
c--- also store value of scale in myscale, which will be retained;
c--- part and scale can be changed to make sure that the tota option works.
      mypart=part
      myscale=scale
      myfacscale=facscale
      mynproc=nproc
      myncall_save=myncall

      if (runstart.gt.0) goto 10
      
c--- If we're doing the tota integration, then set up the grid info
      if ( (mypart .eq. 'tota') .or. (mypart .eq. 'born') ) then
         if (first .and. (myinit .eq. 1)) then
c-- special input name for born kinematic grid
            ingridfile='born.grid'
            myreadin=readin
         else
            if (first .eqv. .true.) then
               readin=.false.
               writeout=.true.
               outgridfile='born.grid'          
            else
               readin=.true.
               writeout=.false.
               ingridfile='born.grid'
            endif
         endif
      endif        
            
c--- Virtual and double virtual integrations have two extra dimensions
c--- (added and then taken away)
      if (  (mypart .eq. 'born') .or. (mypart .eq. 'tota') )  then
        part='born'        
        reset=.true.
        scalereset=.true.
        
        ndim=ndim+2
        if (parallel.eq.1) call setparallel(part,ndim)
        call boundregion(ndim,region)
        call vegasnr(region,ndim,bornint,myinit,myncall,myitmx,
     .       nprn,sigb,sdb,chi)
        ndim=ndim-2
        readin=myreadin
      endif

      if (order.eq.0) goto 20

 10    continue

      
c--- If we're doing the tota integration, then set up the grid info
      if ( (mypart .eq. 'tota') .or. (mypart .eq. 'virt') ) then
         if (first .and. (myinit .eq. 1)) then
c-- special input name for real-virt grid
            ingridfile='virt.grid'
            myreadin=readin
         else
            if (first .eqv. .true.) then
               readin=.false.
               writeout=.true.
               outgridfile='virt.grid'          
            else
               readin=.true.
               writeout=.false.
               ingridfile='virt.grid'
            endif
         endif
      endif        

c--- Real-Virtual integration have two extra dimensions
c--- (added and then taken away)
      if (mypart .eq. 'virt') then
        part='virt'        
        reset=.true.
        scalereset=.true.
        ndim=ndim+3
        
        ndim=ndim+2
        if (parallel.eq.1) call setparallel(part,ndim)
        call boundregion(ndim,region)
        call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
     .       nprn,sigv,sdv,chi)
        ndim=ndim-2
        readin=myreadin
        ndim=ndim-3
      endif

      if ( mypart .eq. 'tota' )  then
        part='virt'        
        reset=.true.
        scalereset=.true.

        if (vonbcalls.gt.0) then
           ncall=myncall*vonbcalls
        else           
           adjust=(dfloat(ndim+3))/(dfloat(ndim+1))
           ncall=int(dfloat(myncall)**adjust)/2
           if (ncall .gt. myncall*10) ncall=myncall*10
        endif
        write(6,*) 'Adjusting number of points for virt to',ncall
        ndim=ndim+3
        ndim=ndim+2
        if (parallel.eq.1) call setparallel(part,ndim)
        call boundregion(ndim,region)
        call vegasnr(region,ndim,virtint,myinit,ncall,myitmx,
     .       nprn,sigv,sdv,chi)
        ndim=ndim-2
        ndim=ndim-3        
        ncall=myncall
        readin=myreadin
      endif

      if ((order.eq.0).or.(runstop.eq.1)) goto 20
      
c--- If we're doing the tota integration, then set up the grid info
      if ( (mypart .eq. 'tota') .or. (mypart .eq. 'real') ) then
        if (first .and. (myinit .eq. 1)) then
c-- special input name for real grid
           ingridfile='real.grid'
           myreadin=readin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='real.grid'          
          else
            readin=.true.
            writeout=.false.
            ingridfile='real.grid'
          endif
        endif        
      endif 

c--- Real integration should have three extra dimensions
      if (mypart .eq. 'real') then
        part='real'
        scalereset=.true.

        ndim=ndim+6
        if (parallel.eq.1) call setparallel(part,ndim)
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,myncall,myitmx,
     .              nprn,sigr,sdr,chi)
        ndim=ndim-6
        readin=myreadin
      endif

      if (mypart .eq. 'tota') then
        scale=myscale
        facscale=myfacscale
        part='real'
        reset=.true.

        if (ronbcalls.gt.0) then
           ncall=myncall*ronbcalls
        else           
           adjust=(dfloat(ndim+6))/(dfloat(ndim+1))
           ncall=int(dfloat(myncall)**(adjust))/2
           if (ncall .gt. myncall*100) ncall=myncall*100
        endif
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+6
        if (parallel.eq.1) call setparallel(part,ndim)
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     .              nprn,sigr,sdr,chi)
        ndim=ndim-6
        ncall=myncall
        readin=myreadin
      endif      

 20    continue
      
      
c--- calculate integration variables to be returned
      xinteg=sigb+sigv+sigr
      xerr=dsqrt(sdb**2+sdv**2+sdr**2)      
      
c--- return part, scale and myncall to their real values
      part=mypart
      scale=myscale
      first=.false.
      
      return
      end
      
      
      subroutine boundregion(idim,region)
c--- Initializes integration region [0,1] for each variable
c--- in the idim-dimensional integration range
      implicit none
      include 'mxdim.f'
      integer i,idim
      double precision region(2*mxdim)
      
      do i=1,idim
      region(i)=0d0
      region(i+idim)=1d0
      enddo
      
      return
      end
      
      

      subroutine setparallel(part,ndim)
      implicit none
      include 'gridinfo.f'
      include 'parallel.f'
      character *9 prename
      character *4 part
      character *4 nfilech,igridch
      integer ndim
      call get_command_argument(4,nfilech)
      read(nfilech,'(I4)') nfile
      call get_command_argument(5,igridch)
      
      
      prename=stage//seedchar//'-'
      if (stage .eq. 'xg1-') then
         readin=.false.
         writeout=.true.
         outgridfile=prename//part//'.grid'
      else
         readin=.true.
         if (stage .eq. 'st2-') then
            write(6,*) 'building grid'
            call cmpnewgridpars(ndim,igridch,part,1)
            ingridfile='cfmc-'//seedchar//'-'//part//'.grid'
            writeout=.false.
         else
            write(6,*) 'building grid'
            call cmpnewgridpars(ndim,igridch,part,0)
            ingridfile='cfmc-'//seedchar//'-'//stage//part//'.grid'             
            writeout=.true.
            outgridfile=prename//part//'.grid'          
         endif
      endif
      return
      end

      
