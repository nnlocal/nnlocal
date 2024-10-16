      real(ki), parameter ::  zip=0._ki
      real(ki), parameter ::  zero=0._ki
      real(ki), parameter ::  one=1._ki
      real(ki), parameter ::  two=2._ki
      real(ki), parameter ::  three=3._ki
      real(ki), parameter ::  four=4._ki
      real(ki), parameter ::  five=5._ki
      real(ki), parameter ::  six=6._ki
      real(ki), parameter ::  seven=7._ki
      real(ki), parameter ::  eight=8._ki
      real(ki), parameter ::  nine=9._ki
      real(ki), parameter ::  ten=10._ki
      real(ki), parameter ::  eleven=11._ki
      real(ki), parameter ::  twelve=12._ki
      real(ki), parameter ::  sixteen=16._ki
      real(ki), parameter ::  half=0.5_ki
      real(ki), parameter ::  quarter=0.25_ki

      real(ki), parameter :: pi=3.1415926535897932384626433832795029_ki
      real(ki), parameter :: pisq=pi*pi
      real(ki), parameter :: pisqo6=pisq/six
      real(ki), parameter :: twopi=two*pi
      real(ki), parameter :: fourpi=four*pi
      real(ki), parameter :: pion4=pi/four
      real(ki), parameter :: pion10=pi/ten
      real(ki), parameter :: pisqm8=pisq-eight

      real(ki), parameter :: rt2=1.41421356237309504880168872420969798_ki
      real(ki), parameter :: twort2=two*rt2
      real(ki), parameter :: fourrt2=four*rt2
! sqrt(two/pi)
      real(ki), parameter ::  rt2onpi=0.79788456080286535587989211986876374_ki 
      real(ki), parameter ::    zeta3=1.20205690315959428539973816151144999_ki
      real(ki), parameter :: plg4half=0.5174790616738993863307581618988629_ki

      complex(ki), parameter :: im=(zip,one)
      complex(ki), parameter :: impi=(zip,pi)
      complex(ki), parameter :: czip=(zip,zip)
      complex(ki), parameter :: cone=(one,zip)
      complex(ki), parameter :: ctwo=(two,zip)
      complex(ki), parameter :: chalf=(half,zip)

      real(ki), parameter :: cf=four/three
      real(ki), parameter :: ca=three
      real(ki), parameter :: xn=three
      real(ki), parameter :: Nc=three
      real(ki), parameter :: Ncinv=one/three
      real(ki), parameter :: xnsq=nine
      real(ki), parameter :: v=eight
      real(ki), parameter :: tr=half
      real(ki), parameter :: Von4=two
      real(ki), parameter :: ninth=one/nine
      real(ki), parameter :: xn4=xnsq-four
      real(ki), parameter :: qu=two/three
      real(ki), parameter :: qd=-one/three
      real(ki), parameter :: qe=-one
      real(ki), parameter :: spinave=one/four
      real(ki), parameter :: aveqq=spinave/xnsq
      real(ki), parameter :: aveqg=spinave/xn/v
      real(ki), parameter :: avegg=spinave/v**2
      real(ki), parameter :: aem=one/137.035989_ki


      real(ki), parameter :: nbGeV2=0.389379e6_ki
      real(ki), parameter :: pbGeV2=0.389379e9_ki
      real(ki), parameter :: fbGeV2=0.389379e12_ki
!----decifemtobarns
      real(ki), parameter :: dfbGeV2=0.389379e13_ki
      real(ki), parameter :: overa=pbGeV2/xn/256._ki/pi
      integer, parameter :: nloop=2, fn=-5
      integer, parameter :: nf=5
      integer, parameter :: mxpart=14

