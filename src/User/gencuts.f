      logical function gencuts(pjet,njets,nd)
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'runstring.f'
      logical failed,gencuts_Zt,gencuts_WZjj,gencuts_input
      integer njets
      double precision pjet(mxpart,4)
      logical my_gencuts
      integer nd
      gencuts=.false.
      gencuts=my_gencuts(pjet,njets,nd)
      return
      end
  
      logical function my_gencuts(pjet,njets,nd)
      implicit none
      include 'constants.f'
      integer njets
      double precision pt,yrap
      double precision pt3,pt4,y3,y4
      double precision pjet(mxpart,4)
      integer nd
      my_gencuts=.false.

c---  symmetric cuts for the decay products of the resonance      
c      y3=yrap(3,pjet)
c      pt3=pt(3,pjet)
c      y4=yrap(4,pjet)
c      pt4=pt(4,pjet)
c      if (abs(y3).gt.2.5) my_gencuts = .true.
c      if (abs(pt3).lt.25) my_gencuts = .true.
c      if (abs(y4).gt.2.5) my_gencuts = .true.
c      if (abs(pt4).lt.25) my_gencuts = .true.

      if (njets.lt.1) my_gencuts = .true.

      return
      end
