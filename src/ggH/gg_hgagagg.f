      subroutine gg_hgagagg(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2)-->H -->  gamma(p3)+gamma(p4)) +g(p_i5=5)+g(p_i6=6) 
c

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_struc.f'
      integer j,k,i5,i6
      double precision p(mxpart,4),Asq,fac
      double precision Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625
c     .                     ,Hgggg_1652,Hgggg_1562,Hgggg_1526
      double precision Hqagg,Haqgg,Hgqqg,Hgaag,Hqgqg,Hagag,Hggqa
      double precision Hggqa_ab,Hggqa_ba,Hggqa_sym
      double precision Hqgqg_ab,Hqgqg_ba,Hqgqg_sym
      double precision Hgqqg_ab,Hgqqg_ba,Hgqqg_sym
      double precision Hagag_ab,Hagag_ba,Hagag_sym
      double precision Hgaag_ab,Hgaag_ba,Hgaag_sym
      double precision Hqagg_ab,Hqagg_ba,Hqagg_sym
      double precision Haqgg_ab,Haqgg_ba,Haqgg_sym
      double precision Hqqqq_a,Hqqqq_b,Hqqqq_i
      double precision Hqaqa_a,Hqaqa_b,Hqaqa_i
      double precision Haqaq_a,Haqaq_b,Haqaq_i
      double precision Hqaaq_a,Hqaaq_b,Hqaaq_i
      double precision 
     . Hqrqr,Hqqqq,
     . Habab,Haaaa,
     . Hqarb,Hqaqa,Hqbqb,
     . Haqbr,Haqaq,Hbqbq,
     . Hqaaq,msqgamgam
      double precision msq(-nf:nf,-nf:nf),hdecay
      parameter(i5=5,i6=6)


C---fill spinor products up to maximum number
      call spinoru(i6,p,za,zb)  

C   Deal with Higgs decay to b-bbar
      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)
      Asq=(as/(3d0*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay

C--four gluon terms
      call h4g(1,2,i5,i6,Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)
      
cC--two quark two gluon terms
c      call hqqggdfm(1,2,i5,i6,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)
c      call hqqggdfm(2,1,i5,i6,Haqgg,Haqgg_ab,Haqgg_ba,Haqgg_sym)
cC====symmetric in first two arguments, but not the ab, ba terms
cc      Haqgg=Hqagg
c
c      call hqqggdfm(1,i5,2,i6,Hqgqg,Hqgqg_ab,Hqgqg_ba,Hqgqg_sym)
c      call hqqggdfm(i5,1,2,i6,Hagag,Hagag_ab,Hagag_ba,Hagag_sym)
cC====symmetric in first two arguments
cc      Hagag=Hqgqg
c
c      call hqqggdfm(2,i5,1,i6,Hgqqg,Hgqqg_ab,Hgqqg_ba,Hgqqg_sym)
c      call hqqggdfm(i5,2,1,i6,Hgaag,Hgaag_ab,Hgaag_ba,Hgaag_sym)
cC====symmetric in first two arguments
cc      Hgaag=Hgqqg
c
c      call hqqggdfm(i6,i5,1,2,Hggqa,Hggqa_ab,Hggqa_ba,Hggqa_sym)
c      
cC---four quark terms
c      call H4qn(1,2,i5,i6,Hqrqr)
c      call H4qi(1,2,i5,i6,Hqqqq,Hqqqq_a,Hqqqq_b,Hqqqq_i)
cC---four anti-quark terms
cc      call H4qn(i5,i6,1,2,Habab)
cc      call H4qi(i5,i6,1,2,Haaaa)
c      Habab=Hqrqr
c      Haaaa=Hqqqq
c
cC-qqb
c      call H4qn(1,i6,2,i5,Hqarb)
c      call H4qi(1,i6,i5,2,Hqaqa,Hqaqa_a,Hqaqa_b,Hqaqa_i)
c      call H4qn(1,i6,i5,2,Hqbqb)
cc      write(6,*) 'Hqaqa',Hqaqa_a,Hqaqa_b,Hqaqa_i
cc      write(6,*) 'Hqbqb',Hqbqb
cC-qbq
c      Haqbr=Hqarb
c      
c      Haqaq=Hqaqa
c      Haqaq_a=Hqaqa_a
c      Haqaq_b=Hqaqa_b
c      Haqaq_i=Hqaqa_i
c      Hbqbq=Hqbqb

      do j=fn,nf
      do k=fn,nf
      msq(j,k)=0d0
      msq_struc(iqr,j,k)=0d0

      if ((j.ne.0).or.(k.ne.0)) cycle

      
c$$$      if ((j.gt.0).and.(k.gt.0)) then 
c$$$        if (j.eq.k) then
c$$$          msq(j,k)=0.5d0*aveqq*fac*Hqqqq
c$$$          msq_struc(iqq_a,j,k)=0.5d0*aveqq*fac*Hqqqq_a
c$$$          msq_struc(iqq_b,j,k)=0.5d0*aveqq*fac*Hqqqq_b
c$$$          msq_struc(iqq_i,j,k)=0.5d0*aveqq*fac*Hqqqq_i
c$$$        else
c$$$          msq(j,k)=aveqq*fac*Hqrqr
c$$$          msq_struc(iqq_a,j,k)=msq(j,k)
c$$$          msq_struc(iqq_b,j,k)=0d0
c$$$          msq_struc(iqq_i,j,k)=0d0
c$$$        endif
c$$$      endif
c$$$      
c$$$      if ((j.lt.0).and.(k.lt.0)) then 
c$$$        if (j.eq.k) then
c$$$          msq(j,k)=0.5d0*aveqq*fac*Haaaa
c$$$        else
c$$$          msq(j,k)=aveqq*fac*Habab
c$$$          msq_struc(iqq_a,j,k)=msq(j,k)
c$$$          msq_struc(iqq_b,j,k)=0d0
c$$$          msq_struc(iqq_i,j,k)=0d0
c$$$        endif
c$$$      endif
c$$$
c$$$      if ((j.gt.0).and.(k.lt.0)) then
c$$$        if (j.eq.-k) then
c$$$          msq(j,k)=aveqq*fac*(0.5d0*Hqagg+Hqaqa+dfloat(nf-1)*Hqarb)
c$$$          msq_struc(iqr,j,k)=aveqq*fac*dfloat(nf-1)*Hqarb
c$$$          msq_struc(iqq_a,j,k)=aveqq*fac*Hqaqa_a
c$$$          msq_struc(iqq_b,j,k)=aveqq*fac*Hqaqa_b
c$$$          msq_struc(iqq_i,j,k)=aveqq*fac*Hqaqa_i
c$$$          msq_struc(igg_ab,j,k)=aveqq*fac*0.5d0*Hqagg_ab
c$$$          msq_struc(igg_ba,j,k)=aveqq*fac*0.5d0*Hqagg_ba
c$$$          msq_struc(igg_sym,j,k)=aveqq*fac*0.5d0*Hqagg_sym
c$$$        else
c$$$          msq(j,k)=aveqq*fac*Hqbqb
c$$$          msq_struc(iqq_a,j,k)=msq(j,k)
c$$$          msq_struc(iqq_b,j,k)=0d0
c$$$          msq_struc(iqq_i,j,k)=0d0
c$$$        endif
c$$$      endif
c$$$
c$$$      if ((j.lt.0).and.(k.gt.0)) then
c$$$        if (j.eq.-k) then
c$$$          msq(j,k)=aveqq*fac*(0.5d0*Haqgg+Haqaq+dfloat(nf-1)*Haqbr)
c$$$          msq_struc(iqr,j,k)=aveqq*fac*dfloat(nf-1)*Haqbr
c$$$          msq_struc(iqq_a,j,k)=aveqq*fac*Haqaq_a
c$$$          msq_struc(iqq_b,j,k)=aveqq*fac*Haqaq_b
c$$$          msq_struc(iqq_i,j,k)=aveqq*fac*Haqaq_i
c$$$          msq_struc(igg_ab,j,k)=aveqq*fac*0.5d0*Haqgg_ab
c$$$          msq_struc(igg_ba,j,k)=aveqq*fac*0.5d0*Haqgg_ba
c$$$          msq_struc(igg_sym,j,k)=aveqq*fac*0.5d0*Haqgg_sym
c$$$        else
c$$$          msq(j,k)=aveqq*fac*Hbqbq
c$$$          msq_struc(iqq_a,j,k)=msq(j,k)
c$$$          msq_struc(iqq_b,j,k)=0d0
c$$$          msq_struc(iqq_i,j,k)=0d0
c$$$        endif
c$$$      endif
c$$$
c$$$      if ((j.gt.0).and.(k.eq.0)) then
c$$$        msq(j,0)=aveqg*fac*Hqgqg
c$$$        msq_struc(igg_ab,j,0)=aveqg*fac*Hqgqg_ab
c$$$        msq_struc(igg_ba,j,0)=aveqg*fac*Hqgqg_ba
c$$$        msq_struc(igg_sym,j,0)=aveqg*fac*Hqgqg_sym
c$$$      endif
c$$$      
c$$$      if ((j.lt.0).and.(k.eq.0)) then
c$$$        msq(j,0)=aveqg*fac*Hagag
c$$$        msq_struc(igg_ab,j,0)=aveqg*fac*Hagag_ab
c$$$        msq_struc(igg_ba,j,0)=aveqg*fac*Hagag_ba
c$$$        msq_struc(igg_sym,j,0)=aveqg*fac*Hagag_sym
c$$$      endif
c$$$
c$$$      if ((j.eq.0).and.(k.gt.0)) then
c$$$        msq(0,k)=aveqg*fac*Hgqqg
c$$$        msq_struc(igg_ab,0,k)=aveqg*fac*Hgqqg_ab
c$$$        msq_struc(igg_ba,0,k)=aveqg*fac*Hgqqg_ba
c$$$        msq_struc(igg_sym,0,k)=aveqg*fac*Hgqqg_sym
c$$$      endif
c$$$
c$$$      if ((j.eq.0).and.(k.lt.0)) then
c$$$        msq(0,k)=aveqg*fac*Hgaag
c$$$        msq_struc(igg_ab,0,k)=aveqg*fac*Hgaag_ab
c$$$        msq_struc(igg_ba,0,k)=aveqg*fac*Hgaag_ba
c$$$        msq_struc(igg_sym,0,k)=aveqg*fac*Hgaag_sym
c$$$      endif

      if ((j.eq.0).and.(k.eq.0)) then
c        msq(0,0)=avegg*fac*(0.5d0*Hgggg+dfloat(nf)*Hggqa)
        msq(0,0)=avegg*fac*(0.5d0*Hgggg)
c        msq_struc(igg_ab,0,0)=avegg*fac*dfloat(nf)*Hggqa_ab
c        msq_struc(igg_ba,0,0)=avegg*fac*dfloat(nf)*Hggqa_ba
c        msq_struc(igg_sym,0,0)=avegg*fac*dfloat(nf)*Hggqa_sym
        msq_struc(igggg_a,0,0)=avegg*fac*0.5d0*Hgggg_1256
        msq_struc(igggg_b,0,0)=avegg*fac*0.5d0*Hgggg_1625
        msq_struc(igggg_c,0,0)=avegg*fac*0.5d0*Hgggg_1265
      endif
      
      enddo
      enddo

cc--- subtraction matrix elements use qa->aq; calculate this and
cc--- artificially store it in msq_struc(iqr,0,0), which is not
cc--- used for anything else
c      call H4qi(1,i5,i6,2,Hqaaq,Hqaaq_a,Hqaaq_b,Hqaaq_i)
c      msq_struc(iqr,0,0)=aveqq*fac*Hqaaq
      
      return
      end

 
