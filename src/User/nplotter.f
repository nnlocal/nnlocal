      subroutine nplotter(p,wt,wt2,nd)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of leptons and jets in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c---     nd:  an integer specifying the dipole number of this contribution
c---          (if applicable), otherwise equal to zero
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'nplot.f'
      include 'ptilde.f'
      double precision p(mxpart,4),wt,wt2
      integer switch,nproc,nd
      common/nproc/nproc
      
      nextnplot = 1

c--- switch:  an integer equal to either 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      if     (nd .gt. ndmax1) then
        switch=2
      elseif (nd .gt. 0) then
        switch=1
      else
        switch=0
      endif

      call nplotter_generic(p,wt,wt2,switch)
      
      return
      end
      
