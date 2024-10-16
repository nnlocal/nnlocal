      program gengridplot
      implicit none
      integer nd,mxdim,ndim,i,j,jj
      PARAMETER (nd=50,mxdim=9999)
      double precision xi(nd,mxdim)
      character *20 ingridfile
      character *39 var
      character *2 ll
      integer filen
      call get_command_argument(1,ingridfile)
      Filen=len(trim(ingridfile))
      write(6,*)
      write(6,*)'****************************************************'
      write(6,*)'      Reading in vegas grid from ',ingridfile
      write(6,*)'****************************************************'
      call flush(6)
c--- read-in grid
      open(unit=11,file=ingridfile,status='unknown')
      ndim=0
      do i=1,mxdim
         ndim=ndim+1
         read(11,*,end=10)
      enddo
 10   continue
      ndim=int(dfloat(ndim-1)/12d0)
      write(6,*) 'ndim=',ndim
      rewind(11)
      open(unit=12,file=ingridfile(1:filen)//'.gnu',status='unknown')
      write(12,'(A)') 'reset'
      write(12,'(A)') 'set terminal postscript col enhanced'
      write(12,'(A)') trim('set output "'//ingridfile(1:filen)//'.ps"')
      write(12,'(A)') 'set style data points'
      write(12,'(A)') 'set key off'
      do j = 1,ndim
         read(11,203) jj,(xi(i,j),i=1,nd)
         write(ll,'(I2)') j
         var='set title "var'//trim(ll)//'" font "Helvetica, 20"'
         write(12,'(A)') var
         write(12,'(A)') 'set xrange [0:1]'
         write(12,'(A)') 'set yrange [0:1]'
         write(12,'(A)') '$Data <<EOD'
         write(12,*) '   0 , 0 '
         do i=1,nd
            write(12,*) xi(i,j),',',dfloat(i)/nd
         enddo
         write(12,'(A)') 'EOD'
         write(12,'(A)') 'set datafile separator comma'
         write(12,*)
         write(12,'(A)') 'set table $Dummy'
         write(12,'(A)') 'plot myXtics=myYtics="" $Data using \'
         write(12,'(A)') '(myXtics = myXtics.($0==0?"":",").strcol(1), \'
         write(12,'(A)') 'myYtics = myYtics.($0==0?"":",").strcol(2)) w table'
         write(12,'(A)') 'unset table'
         write(12,*)
         write(12,'(A)') 'set offsets 1,1,1,1'
         write(12,'(A)') 'set x2tics (@myXtics)'
         write(12,'(A)') 'set grid x2tics ls 12 lc rgb "black"'
         write(12,'(A)') 'set xtics 0.2'
         write(12,'(A)') 'set mxtics 4'
         write(12,'(A)') 'set format x2 ""'
         write(12,'(A)') 'plot $Data u 1:2 with linespoint lw 2 lc rgb "black" ps 0 notitle ,\'
         write(12,'(A)') '     for [j=1:50] (j/50.) lc rgb "black"'
         write(12,*)
         write(12,*)
      enddo
      close(11)
      close(12)      
203   FORMAT(/(5z16))
      stop
      end
