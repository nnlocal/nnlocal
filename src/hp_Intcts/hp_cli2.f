      function hp_cli2(x)
      implicit none
      include 'hp_types.h'
      complex(ki):: hp_cli2
c--complex dilogarithm (spence-function)
      include 'hp_constants.h'
      complex(ki):: x,y,hp_li2taylor,xtmp
      real(ki):: rzeta2,rzeta3,tiny,rr,r2,xr,xi
      common/hp_const/rzeta2,rzeta3
      logical,save:: first=.true.
      include 'hp_cplx.h'
      integer nbes
      common/hp_bersum/nbes
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
      call hp_bernini
      endif
      xr=real(xtmp)
      xi=aimag(xtmp)
      r2=xr*xr+xi*xi
      hp_cli2=cplx2(zip,zip)
      if (((xi/xr).eq.zip).and.(xr.gt.one)) then
         xi=-2e-15_ki
         r2=xr*xr+xi*xi
         xtmp=xr+im*xi
      endif
      if(r2<=tiny)then
        hp_cli2=xtmp+xtmp**2/4._ki
        return
      endif
      rr=xr/r2
      if ((r2==one) .and. (xi==zip)) then
        if (xr==one) then
          hp_cli2=cplx1(rzeta2)
        else
          hp_cli2=-cplx1(rzeta2/two)
        endif
        return
      elseif ((r2>one) .and. (rr>half)) then
        y=(xtmp-one)/xtmp
        hp_cli2=hp_li2taylor(y)+rzeta2-log(xtmp)*log(one-xtmp)+half*log(xtmp)**2
        return
      elseif ((r2>one) .and. (rr<=half))then
        y=one/xtmp
        hp_cli2=-hp_li2taylor(y)-rzeta2-half*log(-xtmp)**2
        return
      elseif ((r2<=one) .and. (xr>half)) then
        y=one-xtmp
        hp_cli2=-hp_li2taylor(y)+rzeta2-log(xtmp)*log(one-xtmp)
       return
      elseif ((r2<=one) .and. (xr<=half)) then
        y=xtmp
        hp_cli2=hp_li2taylor(y)
        return
      endif
      end
 
      function hp_li2taylor(x)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      complex(ki):: hp_li2taylor
c--taylor-expansion for complex dilogarithm (spence-function)
      integer:: nber,i,n
      parameter(nber=30)
      real(ki):: b2(nber)
      complex(ki):: x,z
      common/hp_bernoulli/b2
      integer nbes
      common/hp_bersum/nbes
      n=nbes-1
      z=-log(one-x)
      hp_li2taylor=b2(nbes)
      do 111 i=n,1,-1
      hp_li2taylor=z*hp_li2taylor+b2(i)
111   continue
      hp_li2taylor=z**2*hp_li2taylor+z
      return
      end
 
      function hp_facult(n)
c--real(ki):: version of faculty
      implicit none
      include 'hp_types.h'
      real(ki):: hp_facult
      integer:: i,n
      include 'hp_constants.h'
      hp_facult=one
      if(n==0)return
      do 999 i=1,n
        hp_facult=hp_facult*real(i,ki)
999   continue
      return
      end
 
      subroutine hp_bernini
c--initialization of coefficients for polylogarithms
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      integer:: nber,i
      parameter(nber=30)
      real(ki):: b(nber),b2(nber),rzeta2,rzeta3,hp_facult
      common/hp_bernoulli/b2
      common/hp_const/rzeta2,rzeta3
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
        b2(i)=b(i)/hp_facult(i+1)
995   continue
      return
      end
 
