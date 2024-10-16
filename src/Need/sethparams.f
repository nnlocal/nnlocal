      subroutine sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
c--- set up the necessary parameters for a Standard Model Higgs boson
c---   hwidth : either the NLO value from Spira,
c---                or the LO value calculated here
c---   br,wwbr,zzbr,tautaubr : the LO calculated values
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'process.f' 
      include 'verbose.f'
      include 'cpscheme.f'
      double precision br,gamgambr,zgambr,wwbr,zzbr,tautaubr,x_w,x_z,
     & msqgamgam,hzgamwidth,msqhbb,msqhtautau,
     & pw_bb,pw_tautau,pw_gamgam,pw_ww,pw_zz,pw_zgam,
     & br_sp,tautaubr_sp,gamgambr_sp,wwbr_sp,zzbr_sp,zgambr_sp
      logical spira
      common/spira/spira

***************************** COMPUTE PARTIAL WIDTHS ***************************
*                                                                              *
*  Included decays are:                                                        *
*                                                                              *
*    H -> bb, H -> tau tau                                                     *
*    H -> WW, H -> ZZ, H -> gamma gamma                                        *
*                                                                              *
********************************************************************************

c--- Compute partial width of H -> bb : note that the mass is mostly retained
c--- in the coupling only, since mb=0 usually (the mass in the phase space)
      pw_bb=msqhbb(hmass**2)/(16d0*pi*hmass)
     &     *sqrt(1d0-4d0*mb**2/hmass**2)

c--- Compute partial width of H -> tau^- tau^+
      pw_tautau=msqhtautau(hmass**2)/(16d0*pi*hmass)
     &     *sqrt(1d0-4d0*mtau**2/hmass**2)

      x_w=4d0*wmass**2/hmass**2
      x_z=4d0*zmass**2/hmass**2
c--- Compute partial width of H -> WW
      if (x_w .lt. 1d0) then
        pw_ww=gwsq/64d0/pi*hmass**3/wmass**2
     &       *dsqrt(1d0-x_w)*(1d0-x_w+0.75d0*x_w**2)
      else
        pw_ww=0d0
      endif
c--- Compute partial width of H -> ZZ
      if (x_z .lt. 1d0) then
        pw_zz=gwsq/128d0/pi*hmass**3/wmass**2
     &       *dsqrt(1d0-x_z)*(1d0-x_z+0.75d0*x_z**2)
      else
        pw_zz=0d0
      endif

c--- Compute partial width of H -> gamma gamma
      pw_gamgam=msqgamgam(hmass)/(16d0*pi*hmass)
      
c--- Compute partial width of H -> Z gamma
c      pw_zgam=hzgamwidth(hmass)
      
****************************** COMPUTE TOTAL WIDTHS ****************************

c--- Calculate the value of hwidth, depending on "spira"
c--- note that the branching ratios "br_sp","gamgambr_sp","wwbr_sp","zzbr_sp"
c--- are not actually used in our calculations   
      if (spira) then
        call higgsp(br_sp,tautaubr_sp,gamgambr_sp,zgambr_sp,wwbr_sp,
     &              zzbr_sp)
      else
        hwidth=pw_bb+pw_tautau+pw_ww+pw_zz+pw_gamgam+pw_zgam
      endif 

c--- complex pole scheme, if desired
      CPscheme=.false.
      if (CPscheme) then
        call interpolate_hto(hmass,hwidth)
      endif


**************************** COMPUTE BRANCHING RATIOS **************************

c--- Branching ratio H -> bb
      br=pw_bb/hwidth
c--- Branching ratio H -> tau tau
      tautaubr=pw_tautau/hwidth
c--- Branching ratio H -> WW
      wwbr=pw_ww/hwidth
c--- Branching ratio H -> ZZ
      zzbr=pw_zz/hwidth
c--- Branching ratio H -> gamgam
      gamgambr=pw_gamgam/hwidth
c--- Branching ratio H -> Zgam
      zgambr=pw_zgam/hwidth

*************************** WRITE OUT BRANCHING RATIOS *************************


      if (verbose) then
       write(6,99) hmass,hwidth,br,tautaubr,wwbr,zzbr,gamgambr,zgambr
       if (spira) then
       write(6,*) '*                                                  *'
       write(6,*) '* Note: branching ratios reported here can be > 1  *'
       write(6,*) '*       since the total Higgs width is calculated  *'
       write(6,*) '*       at NLO but the BR calculation uses a       *'
       write(6,*) '*       partial width at LO only.                  *'
       write(6,*) '*                                                  *'
       write(6,*) '****************************************************'
       endif
      endif

      return

 99   format(/,
     .       ' ****************** Higgs parameters ****************'/, 
     .       ' *                                                  *'/, 
     .       ' *   mass(H) = ',f7.2,'      width(H) = ',e12.5,' *'/,
     .       ' *                                                  *'/, 
     .       ' *              Br( H -> b bbar)  = ',f9.5,'       *'/,
     .       ' *              Br( H -> tau tau) = ',f9.5,'       *'/,
     .       ' *              Br( H -> W W)     = ',f9.5,'       *'/,
     .       ' *              Br( H -> Z Z)     = ',f9.5,'       *'/,
     .       ' *              Br( H -> gam gam) = ',f9.5,'       *'/,
     .       ' *              Br( H -> Z gam)   = ',f9.5,'       *'/,
     .       ' ****************************************************')

      end
