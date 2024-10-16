      subroutine nplotter_generic(p,wt,wt2,switch)
      implicit none
      include 'vegas_common.f'
      include 'bbproc.f'
      include 'clustering.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'process.f'
      include 'removebr.f'
      include 'nodecay.f'
      include 'useet.f'
      include 'plabel.f'
      include 'outputflags.f'
cz Add single-top b fraction, Z. Sullivan 1/25/05
      double precision bwgt
      common/btagging/ bwgt

cz //

      integer n,switch,i5,i6,i7,nu,nplotmax
      character tag*4
      double precision DETABB,DPHIZJ
     & ,DPHIBB
     & ,ETAB1
     & ,ETAB2
     & ,ETANOB
     & ,ETARAP
     & ,ETARAPTWO
     & ,ETARAPTHREE
     & ,ETBIN
     & ,ETDOUBLEBIN
     & ,GETET
     & ,M56CLUST
     & ,MBB
     & ,PT
     & ,PTB1
     & ,PTB2
     & ,PTNOB
     & ,PTQ1
     & ,ETAQ1
     & ,PTOTHER
     & ,ETAOTHER
     & ,PTTWO
     & ,PTTHREE
     & ,R
     & ,R57
     & ,R67
     & ,RBB
     & ,SWAP
     & ,WT
     & ,WT2
     & ,YRAP
     & ,YRAPTWO
     & ,YRAPTHREE
     & ,DPHI_LL
     & ,M_LL
     & ,MTRANS
     & ,SCUT1
     & ,SCUT2
     & ,PT34_VETO
     & ,PTQ1_VETO
     & ,ETAQ1_VETO
     & ,PHI56
     & ,PTLEADINGB
     & ,PTLEADINGNONB

      integer jet(mxpart),jetstart,ibbar,iz,izj,iztmp,
     . inotb,ilight1,ilight2
cz //

      double precision wtbbar,wtnotb,wtlight1,wtlight2
      double precision m34,etmiss,misset,m345678,
     . p(mxpart,4),fphi,HT,etcharm,deltaeta,cosdeltaphi
      double precision eta3,eta4,eta5,eta6,eta7,eta8,eta9,eta10,
     . eta34,eta56,y34,ay34
      double precision r45,r56,m345,m348,m567,m678
      double precision pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,
     . pt34,pt56,oldpt(5:7)
      double precision pt345,eta345,y345,y3,y4
      double precision bclustmass,etvec(4),tmp5(4),tmp6(4),tmp7(4)
      integer ssi(4),tmpi,j,k
      double precision ssd(4),tmpd
      integer nproc,eventpart,ib1,ib2,nqcdjets,nqcdstart
      logical first,jetmerge
      logical jetevent
      character*2 ptet
      common/nplotmax/nplotmax
      common/nproc/nproc
      common/nqcdjets/nqcdjets,nqcdstart
      common/jetmerge/jetmerge
      common/hwwvars/dphi_ll,m_ll,mtrans,scut1,scut2
      data first/.true./
      save first
      
c--- Set up string for pt or Et
      if (useEt) then
        ptet='Et'
      else
        ptet='pt'
      endif
      
      if (first) then
        tag='book'
c--- ensure we initialize all possible histograms
        eventpart=npart+2
        eta3=1d3
        pt3=0d0
        eta4=1d3
        pt4=0d0
        eta5=1d3
        pt5=0d0
        eta6=1d3
        pt6=0d0
        eta7=1d3
        pt7=0d0
        eta8=1d3
        pt8=0d0
        eta9=1d3
        pt9=0d0
        eta10=1d3
        pt10=0d0
        eta34=1d3
        y34=1d3
        ay34=1d3
        pt34=0d0
        r45=0d0
        eta56=1d3
        pt56=0d0        
        m56clust=0d0
        r56=0d0
        r57=0d0
        r67=0d0
        misset=0d0
        etbin=0d0
        mbb=0d0
        etab1=1d3
        etab2=1d3
        etanob=1d3
        ptb1=0d0
        ptb2=0d0
        ptQ1=0d0
        etaQ1=0d0
        ptother=0d0
        etaother=1d3
        ptnob=0d0
        rbb=0d0
        detabb=1d3
        dphibb=1d3
        m345=0d0
        m348=0d0
        m567=0d0
        m678=0d0
        m345678=0d0
        HT=0d0
      deltaeta=99d0
      cosdeltaphi=99d0
      ssi(1)=5
      ssi(2)=6
        jetmerge=.true.
        jets=2
        goto 99
      else
        tag='plot'
      endif

c--- 'eventpart' will contain the number of actual particles that have
c--- a defined momentum. For most processes, this is calculated as follows:
c----  for lowest order and virtual terms switch=0 and eventpart=npart+2
c---   for real events switch=0 and eventpart=npart+2
c---   for real counter-events switch=1 and eventpart=npart+1
c--- There are some processes for which this is not correct and these
c---  are handled with reference to nproc  
c      eventpart=npart-switch+2

      eventpart=4
      if (jets .gt. 0) then
        eventpart=4+jets
      endif

      
c--- this variable should be set to .true. when the jets are reordered
c---  according to their pt (or Et)
      jetevent=.true.
      if (algorithm .eq. 'cone') then
         if (jets .gt. 0) pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))
         if (jets .gt. 1) pt6=getet(p(6,4),p(6,1),p(6,2),p(6,3))
c         if (jets .gt. 2) pt7=getet(p(7,4),p(7,1),p(7,2),p(7,3))
      else
         if (jets .gt. 0) pt5=pt(5,p)
         if (jets .gt. 1) pt6=pt(6,p)
c         if (jets .gt. 2) pt7=pt(7,p)        
      endif
      i5=5
      i6=6
      i7=7
      if (jets .gt. 0) oldpt(5)=pt5
      if (jets .gt. 1) oldpt(6)=pt6
c      if (jets .gt. 2) oldpt(7)=pt7
c--- sort for 2 jets 
      if (jets .eq. 2) then          
         if (pt6 .gt. pt5) then
            i5=6
            i6=5
         endif
      endif
c$$$c---  sort for 3 jets 
c$$$      if (jets .eq. 3) then
c$$$         if ((pt5 .gt. pt6) .and. (pt5 .gt. pt7)) then
c$$$            i5=5
c$$$            if (pt6 .gt. pt7) then
c$$$               i6=6
c$$$               i7=7
c$$$            else
c$$$               i6=7
c$$$               i7=6
c$$$            endif
c$$$         endif
c$$$         if ((pt6 .gt. pt5) .and. (pt6 .gt. pt7)) then
c$$$            i5=6
c$$$            if (pt5 .gt. pt7) then
c$$$               i6=5
c$$$               i7=7
c$$$            else
c$$$               i6=7
c$$$               i7=5
c$$$            endif
c$$$         endif
c$$$         if ((pt7 .gt. pt5) .and. (pt7 .gt. pt6)) then
c$$$            i5=7
c$$$            if (pt5 .gt. pt6) then
c$$$               i6=5
c$$$               i7=6
c$$$            else
c$$$               i6=6
c$$$               i7=5
c$$$            endif
c$$$         endif
c$$$      endif
c$$$c---  perform exchange
      do nu=1,4
         tmp5(nu)=p(i5,nu)
         tmp6(nu)=p(i6,nu)
c$$$         tmp7(nu)=p(i7,nu)
      enddo
      do nu=1,4
         p(5,nu)=tmp5(nu)
         p(6,nu)=tmp6(nu)
c$$$         p(7,nu)=tmp7(nu)
      enddo
      if (jets .gt. 0) pt5=oldpt(i5)
      if (jets .gt. 1) pt6=oldpt(i6)
c$$$      if (jets .gt. 2) pt7=oldpt(i7)


      m34=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)


      eta3=etarap(3,p)
      y3=yrap(3,p)
      pt3=pt(3,p)
      eta4=etarap(4,p)
      y4=yrap(4,p)
      pt4=pt(4,p)        
      eta34=etaraptwo(3,4,p)
      y34=yraptwo(3,4,p)
      ay34=dabs(y34)
      pt34=pttwo(3,4,p)
      HT=pt3+pt4
      
      if (eventpart .gt. 4) then        
      eta5=etarap(5,p)
      if (jetevent .eqv. .false.) pt5=pt(5,p)
      r45=R(p,4,5)
      m345=dsqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     .          -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
      pt345=ptthree(3,4,5,p)
      eta345=etarapthree(3,4,5,p)
      y345=yrapthree(3,4,5,p)
      HT=HT+pt5
      endif
      
      if (eventpart .gt. 5) then
         eta6=etarap(6,p)
         if (jetevent .eqv. .false.) pt6=pt(6,p)
         eta56=etaraptwo(5,6,p)
         eta56=etaraptwo(5,6,p)
         pt56=pttwo(5,6,p)
         r56=R(p,5,6)
         phi56=fphi(5,6,p)
         HT=HT+pt6
      endif

c      if (eventpart .gt. 6) then        
c      eta7=etarap(7,p)
c      if (jetevent .eqv. .false.) pt7=pt(7,p)
c      r57=R(p,5,7)
c      r67=R(p,6,7)
c      m567=dsqrt((p(5,4)+p(6,4)+p(7,4))**2-(p(5,1)+p(6,1)+p(7,1))**2
c     .          -(p(5,2)+p(6,2)+p(7,2))**2-(p(5,3)+p(6,3)+p(7,3))**2)
c      HT=HT+pt7
c      endif

c      misset=etmiss(p,etvec)


c--- find largest rapidity difference between the jets
      deltaeta=99d0
      cosdeltaphi=99d0
      if ((jets .eq. 2) .or. (jets .eq. 3)) then
c         if     (jets .eq. 2) then
c            i5=5
c            i6=6
c         elseif (jets .eq. 3) then
c             if (abs(eta5-eta6).gt.max(abs(eta5-eta7),abs(eta6-eta7))) then
c                 i5=5
c                 i6=6
c             else
c                 if (abs(eta5-eta7).gt.abs(eta6-eta7)) then
c                    i5=5
c                    i6=7
c                 else
c                    i5=6
c                    i6=7
c                 endif
c             endif
c         endif
      deltaeta=abs(etarap(i5,p)-etarap(i6,p))
      cosdeltaphi=(p(i5,1)*p(i6,1)+p(i5,2)*p(i6,2))
     .           /dsqrt((p(i5,1)**2+p(i5,2)**2)*(p(i6,1)**2+p(i6,2)**2))
        if (cosdeltaphi .lt. -0.999999999D0) cosdeltaphi=-1d0
c      cosdeltaphi=abs(cosdeltaphi)
      endif
        
   99 continue

      
c      write(6,*) jets,nqcdjets

c--- Otherwise, fill the histograms 
      n=nextnplot 
      
      call bookplot(n,tag,'HT',HT,wt,wt2,0d0,500d0,20d0,'lin')
      n=n+1
c --- Histograms to monitor exclusive/inclusive cross-sections:
      if ((jets .eq. nqcdjets)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'= #LO j ',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
      if ((jets .ge. nqcdjets)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'>= #LO j',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
c --- Histograms to monitor exclusive/inclusive cross-sections:
      if ((jets .eq. nqcdjets+1)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'= #NLO j ',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
      if ((jets .ge. nqcdjets+1)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'>= #NLO j',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1

c --- Histograms to monitor exclusive/inclusive cross-sections:
      if ((jets .eq. nqcdjets+2)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'= #NNLO j ',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1


c      write(6,*) 'nodecay=',nodecay
c      pause
c      if (nodecay .eqv. .false.) then      

c      call bookplot(n,tag,'eta3',eta3,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'y3',y3,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,ptet//'3',pt3,wt,wt2,0d0,150d0,5d0,'log')
c      n=n+1      
c      call bookplot(n,tag,'eta4',eta4,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'y4',y4,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,ptet//'4',pt4,wt,wt2,25d0,50d0,0.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,150d0,5d0,'log')
c      n=n+1
c      call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,600d0,20d0,'log')
c      n=n+1
      call bookplot(n,tag,'eta34',eta34,wt,wt2,-8d0,8d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5d0,5d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5d0,5d0,0.25d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5d0,5d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5d0,5d0,1d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-5.5d0,5.5d0,1d0,'lin')
      n=n+1
      call bookplot(n,tag,'ay34',ay34,wt,wt2,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'ay34',ay34,wt,wt2,-5d0,5d0,1d0,'lin')
      n=n+1
      call bookplot(n,tag,'ay34',ay34,wt,wt2,-5.5d0,5.5d0,1d0,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'34',pt34,wt,wt2,0d0,100d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'34',pt34,wt,wt2,0d0,40d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,102.5d0,152.5d0,5d0,'lin')
      n=n+1
c      call bookplot(n,tag,'m34',m34,wt,wt2,0d0,500d0,10d0,'log')
c      n=n+1

c      endif

c      call bookplot(n,tag,'misset',misset,wt,wt2,0d0,100d0,2d0,'lin')
c      n=n+1
      
      if (eventpart .gt. 4) then
c      call bookplot(n,tag,'eta5',eta5,wt,wt2,-6d0,6d0,0.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'eta5',eta5,wt,wt2,-10d0,10d0,0.5d0,'lin')
c      n=n+1
      call bookplot(n,tag,ptet//'5',pt5,wt,wt2,0d0,500d0,10d0,'lin')
      n=n+1
c      call bookplot(n,tag,ptet//'5',pt5,wt,wt2,0d0,100d0,2d0,'log')
c      n=n+1
c      call bookplot(n,tag,'pt345',pt345,wt,wt2,0d0,200d0,5d0,'log')
c      n=n+1
c      call bookplot(n,tag,'eta345',eta345,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'y345',y345,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'m345',m345,wt,wt2,0d0,200d0,5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'R45',R45,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
      endif
      
      if (eventpart .gt. 5) then
      call bookplot(n,tag,'eta6',eta6,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'6',pt6,wt,wt2,0d0,500d0,10d0,'log')
      n=n+1
c      call bookplot(n,tag,ptet//'6',pt6,wt,wt2,0d0,200d0,2d0,'log')
c      n=n+1
      call bookplot(n,tag,'eta56',eta56,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'56',pt56,wt,wt2,10d0,250d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'r56',r56,wt,wt2,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'phi56',phi56,wt,wt2,0d0,3.14d0,0.314d0,'lin')
      n=n+1
      endif


      if (nqcdjets .ge. 2) then
      call bookplot(n,tag,'delta(eta) >= 2 jets',deltaeta,
     . wt,wt2,0d0,10d0,0.5d0,'lin')
      n=n+1      
      endif

      if (eventpart .gt. 6) then
      call bookplot(n,tag,'eta7',eta7,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'7',pt7,wt,wt2,0d0,100d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'7',pt7,wt,wt2,0d0,600d0,20d0,'log')
      n=n+1
      call bookplot(n,tag,'r57',r57,wt,wt2,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'r67',r67,wt,wt2,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'m567',m567,wt,wt2,108d0,308d0,2d0,'log')
      n=n+1
      call bookplot(n,tag,'eta diff',eta7-(eta5+eta6)/2d0,
     . wt,wt2,-6d0,6d0,0.2d0,'lin')
      n=n+1
      endif      

c--- plots useful for processes with 2 or more jets only
      if (nqcdjets .ge. 2) then
      call bookplot(n,tag,'delta(eta) >= 2 jets',deltaeta,
     . wt,wt2,0d0,10d0,0.5d0,'lin')
      n=n+1      
      if ((jets .eq. 2) .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'delta(eta) - 2 jets',deltaeta,
     . wt,wt2,0d0,10d0,0.5d0,'lin')
      endif
      n=n+1      
      if ((jets .eq. 2) .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'for phi: delta(eta) - 2 jets',deltaeta,
     . wt*cosdeltaphi,wt2,0d0,10d0,0.5d0,'lin')
      endif
      n=n+1      
      if ((jets .eq. 2) .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'cos delta(phi) - 2 jets',cosdeltaphi,
     . wt,wt2,0d0,1d0,0.05d0,'lin')
      endif
      n=n+1      
      if ((jets .eq. 3) .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'delta(eta) - 3 jets',deltaeta,
     . wt,wt2,0d0,10d0,0.5d0,'lin')
      endif
      n=n+1      
      if ((jets .eq. 3) .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'for phi: delta(eta) - 3 jets',deltaeta,
     . wt*cosdeltaphi,wt2,0d0,10d0,0.5d0,'lin')
      endif
      n=n+1      
      if ((jets .eq. 3) .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'cos delta(phi) - 3 jets',cosdeltaphi,
     . wt,wt2,0d0,1d0,0.05d0,'lin')
      endif
      n=n+1
      endif
            
      n=n-1
      
c--- ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
cz      
      bwgt=0d0  ! for safety
cz //

      return 
      end
      
      subroutine cross(p,i,j,r)
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),r(3)
      
      r(1)=p(i,2)*p(j,3)-p(j,2)*p(i,3)
      r(2)=p(i,3)*p(j,1)-p(j,3)*p(i,1)
      r(3)=p(i,1)*p(j,2)-p(j,1)*p(i,2)
      
      return
      end
            
      double precision function fphi(n1,n2,p)
      implicit none
      include 'constants.f'
      integer n1,n2
      double precision p(mxpart,4)
    
      fphi=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)
      fphi=fphi/dsqrt(p(n1,1)**2+p(n1,2)**2)
      fphi=fphi/dsqrt(p(n2,1)**2+p(n2,2)**2)
      if     (fphi .gt. +0.9999999D0) then
        fphi=0d0
      elseif (fphi .lt. -0.9999999D0) then
        fphi=pi
      else
        fphi=dacos(fphi)
      endif

      return
      end
                     
      double precision function coslpairet(n1,n2,nm1,nm2,p)
      implicit none
      include 'constants.f'
      integer n1,n2,nm1,nm2,i
      double precision p(mxpart,4),misset(4),pp(4)
      
      if (nm2 .eq. 0) then
        do i=1,4
          misset(i)=p(nm1,i)
        enddo
      else
        do i=1,4
          misset(i)=p(nm1,i)+p(nm2,i)
        enddo
      endif
      
      do i=1,4
        pp(i)=p(n1,i)+p(n2,i)
      enddo
      
      coslpairet=pp(1)*misset(1)+pp(2)*misset(2)
      coslpairet=coslpairet/dsqrt(pp(1)**2+pp(2)**2)
      coslpairet=coslpairet/dsqrt(misset(1)**2+misset(2)**2)
      
      return
      end
            
      double precision function deltar(i,j,p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),phi1,phi2,etarap,dphi
      integer i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(p(j,1),p(j,2))
      dphi=phi1-phi2
      if (dphi .gt. pi) dphi=twopi-dphi
      if (dphi .lt. -pi) dphi=twopi+dphi
      deltar=(etarap(i,p)-etarap(j,p))**2+dphi**2
      deltar=dsqrt(deltar)
      
      return
      end
      
      subroutine checkmaxhisto(n)
c--- ensure the built-in maximum number of histograms is not exceeded    
      implicit none
      integer n
      include 'histo.f'
      
      if (n .gt. maxhisto) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > ',maxhisto,', which is the built-in maximum.'
      write(6,*) 'To use more histograms, change the value of the'
      write(6,*) 'constant MAXHISTO in src/Inc/nplot.f then do:'
      write(6,*)
      write(6,*) ' make clean; make        to recompile from scratch.'
        write(6,*)
      stop
      endif
      
      return
      end
      
