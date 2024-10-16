      function cli2(x)
      implicit none
      include 'types.h'
      complex(ki):: cli2
c--complex dilogarithm (spence-function)
      include 'constants.h'
      complex(ki):: x,y,li2taylor,xtmp
      real(ki):: rzeta2,rzeta3,tiny,rr,r2,xr,xi
      common/const/rzeta2,rzeta3
      logical,save:: first=.true.
      include 'cplx.h'
      integer nbes
      common/bersum/nbes
      nbes=18
      if ( ki == sp ) then
         tiny=1.e-4_ki
      elseif ( ki == dp ) then
         tiny=1.e-8_ki
      elseif ( ki == ex ) then
         tiny=1.e-10_ki
      elseif ( ki == qp ) then
         tiny=1.e-16_ki
         nbes=30
      endif
      xtmp=x
      if (first) then
      first=.false.
      call bernini
      endif
      xr=real(xtmp)
      xi=aimag(xtmp)
      r2=xr*xr+xi*xi
      cli2=cplx2(zip,zip)
      if (((xi/xr).eq.zip).and.(xr.gt.one)) then
         xi=-2e-28_ki
         r2=xr*xr+xi*xi
         xtmp=xr+im*xi
      endif
      if(r2<=tiny)then
        cli2=xtmp+xtmp**2/4._ki
        return
      endif
      rr=xr/r2
      if ((r2==one) .and. (xi==zip)) then
        if (xr==one) then
          cli2=cplx1(rzeta2)
        else
          cli2=-cplx1(rzeta2/two)
        endif
        return
      elseif ((r2>one) .and. (rr>half)) then
        y=(xtmp-one)/xtmp
        cli2=li2taylor(y)+rzeta2-log(xtmp)*log(one-xtmp)+half*log(xtmp)**2
        return
      elseif ((r2>one) .and. (rr<=half))then
        y=one/xtmp
        cli2=-li2taylor(y)-rzeta2-half*log(-xtmp)**2
        return
      elseif ((r2<=one) .and. (xr>half)) then
        y=one-xtmp
        cli2=-li2taylor(y)+rzeta2-log(xtmp)*log(one-xtmp)
       return
      elseif ((r2<=one) .and. (xr<=half)) then
        y=xtmp
        cli2=li2taylor(y)
        return
      endif
      end
 
      function li2taylor(x)
      implicit none
      include 'types.h'
      include 'constants.h'
      complex(ki):: li2taylor
c--taylor-expansion for complex dilogarithm (spence-function)
      integer:: nber,i,n
      parameter(nber=30)
      real(ki):: b2(nber)
      complex(ki):: x,z
      common/bernoulli/b2
      integer nbes
      common/bersum/nbes
      n=nbes-1
      z=-log(one-x)
      li2taylor=b2(nbes)
      do 111 i=n,1,-1
      li2taylor=z*li2taylor+b2(i)
111   continue
      li2taylor=z**2*li2taylor+z
      return
      end
 
      function facult(n)
c--real(ki):: version of faculty
      implicit none
      include 'types.h'
      real(ki):: facult
      integer:: i,n
      include 'constants.h'
      facult=one
      if(n==0)return
      do 999 i=1,n
        facult=facult*real(i,ki)
999   continue
      return
      end
 
      subroutine bernini
c--initialization of coefficients for polylogarithms
      implicit none
      include 'types.h'
      include 'constants.h'
      integer:: nber,i
      parameter(nber=30)
      real(ki):: b(nber),b2(nber),rzeta2,rzeta3,facult
      common/bernoulli/b2
      common/const/rzeta2,rzeta3
      b(1)=-1._ki/2._ki
      b(2)=1._ki/6._ki
      b(3)=0._ki
      b(4)=-1._ki/30._ki
      b(5)=0._ki
      b(6)=1._ki/42._ki
      b(7)=0._ki
      b(8)=-1._ki/30._ki
      b(9)=0._ki
      b(10)=5._ki/66._ki
      b(11)=0._ki
      b(12)=-691._ki/2730._ki
      b(13)=0._ki
      b(14)=7._ki/6._ki
      b(15)=0._ki
      b(16)=-3617._ki/510._ki
      b(17)=0._ki
      b(18)=43867._ki/798._ki
      b(19)=0._ki
      b(20)=-174611._ki/330._ki
      b(21)=0._ki
      b(22)=854513._ki/138._ki
      b(23)=0._ki
      b(24)=-236364091._ki/2730._ki
      b(25)=0._ki
      b(26)=8553103._ki/6._ki
      b(27)=0._ki
      b(28)=-23749461029._ki/870._ki
      b(29)=0._ki
      b(30)=8615841276005._ki/14322._ki
      rzeta2=pi*pi/6._ki
      rzeta3=1.2020569031595942853997381615114500_ki
      do 995 i=1,nber
        b2(i)=b(i)/facult(i+1)
995   continue
      return
      end
 
