      subroutine setrunname()
      implicit none
      include 'part.f'
      include 'parallel.f'
      integer nlength,lenocc
      character*9 runname
      common/runname/runname
      common/nlength/nlength
      runname=part
      if (parallel.eq.1) runname=runname(1:4)//'-'//seedchar
      nlength=lenocc(runname)
      return
      end


      character*3 function getstr(no)
c returns a string of length 3 from an integer
      integer no,i1,i2,i3,zero
      
      zero=ichar('0')

      i1=abs(no)/100
      i2=(abs(no)-i1*100)/10
      i3=abs(no)-i1*100-i2*10

      if    (i1.eq.0.and.i2.eq.0) then
        if (no .lt. 0) then
        getstr='-'//char(i3+zero)//'_'
        else
        getstr=char(i3+zero)//'__'
        endif
      elseif(i1.eq.0) then
        getstr=char(i2+zero)//char(i3+zero)//'_'
      else
        getstr=char(i1+zero)//char(i2+zero)//char(i3+zero)
      endif
      
      return
      end

