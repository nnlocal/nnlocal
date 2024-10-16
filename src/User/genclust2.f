      subroutine genclust2(q,R,qfinal,isub)
c--- this is a wrapper routine for the jet clustering algorithm
c--- which re-routes according to the value of 'algorithm'  to:
c---  ('ktal') genclust_kt.f     for kt clustering
c---  ('ankt') genclust_kt.f     for "anti-kt" clustering
c---  ('cone') genclust_cone.f   for cone algorithm
c---  ('none') to perform no clustering at all
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'bbproc.f'
      include 'process.f'
      include 'part.f'
      integer nproc
      common/nproc/nproc
      character*4 mypart
      double precision q(mxpart,4),qfinal(mxpart,4),
     & qreorder(mxpart,4),R,Rbbmin
      integer nqcdjets,nqcdstart,notag,isub,i,nu,njetsmin,njetsmax
      logical first
      common/nqcdjets/nqcdjets,nqcdstart
      common/notag/notag
      common/Rbbmin/Rbbmin
      common/mypart/mypart
      data first/.true./
      save first
      
      if ((first) .and.
     &     ((nqcdjets .gt. 0).or.(part .eq. 'real').or.
     &    (part .eq. 'virt').or.(notag.gt.0))) then
        first=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
      write(6,*)
      write(6,*) '*********** Basic jet-defining parameters **********'
      if     (algorithm .eq. 'ktal') then
      write(6,*) '*          (Run II kT clustering algorithm)        *'
      elseif (algorithm .eq. 'ankt') then
      write(6,*) '*     (Anti-kt algorithm - see arXiv:0802.1189)    *'
      elseif (algorithm .eq. 'cone') then
      write(6,*) '*              (Run II cone algorithm)             *'
      elseif (algorithm .eq. 'none') then
      write(6,*) '*             (no clustering algorithm)            *'
      else
      write(6,*)
      write(6,*) 'Invalid selection of algorithm in input file.'
      write(6,*) 'Please select either ktal, cone, hqrk or none'
      stop
      endif
      write(6,*) '*                                                  *'
      write(6,79) ' *     pt(jet)         > ',ptjetmin
      write(6,79) ' *   |pseudo-rap(jet)| > ',etajetmin   
      write(6,79) ' *   |pseudo-rap(jet)| < ',etajetmax   
      if (bbproc) then
        ptbjetmin=max(ptjetmin,ptbjetmin)
        etabjetmax=min(etajetmax,etabjetmax)
      write(6,79) ' *   pt(b-jet)         > ',ptbjetmin
      write(6,79) ' * |pseudo-rap(b-jet)| < ',etabjetmax   
      endif
      write(6,79) ' * pseudo-cone size, R : ',R
      write(6,*) '*                                                  *'
      njetsmin=nqcdjets-notag

c---  TODO: Fix njetsmin and njetsmax as a function of
c---  process and order (notag)
      
      if (inclusive) then
        if( (part.eq.'real') .or. (mypart.eq.'tota')
     &  .or.(mypart.eq.'todk') )then
          njetsmax=nqcdjets+2
        else
          njetsmax=nqcdjets
        endif
      else
        njetsmax=nqcdjets-notag
      endif
      write(6,78) njetsmin,njetsmax
      write(6,*) '****************************************************'
      call flush(6)
      endif
   78 format(' *    Cross-section defined by:  ',i2,' <= jets <=',
     &        i2,'    *')
   79 format(a25,f8.4,'                   *')

      if     (algorithm .eq. 'ktal') then
        call genclust_kt(q,R,qfinal,isub,+1)
      elseif (algorithm .eq. 'ankt') then
        call genclust_kt(q,R,qfinal,isub,-1)
      elseif (algorithm .eq. 'cone') then
        call genclust_cone(q,R,qfinal,isub)
      elseif (algorithm .eq. 'none') then
        do i=1,mxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        jets=nqcdjets
        if ((part .eq. 'virt') .and. (isub .eq. 0)) jets=jets+1
        if ((part .eq. 'real') .and. (isub .eq. 0)) jets=jets+2
        if (nproc.eq.709) jets=jets-1
        return
      else
        write(6,*) 'Invalid choice of jet algorithm, must be'
        write(6,*) '   ktal, ankt, cone, none'
        stop
      endif


      return
      end
      
