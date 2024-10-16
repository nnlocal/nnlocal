      program cmptot
      implicit none
      double precision a,b,val,err
      integer i,j
      double precision res,rer
      character *24 filename,var
      character *60 line
      call get_command_argument(1,filename)
      open(20,file=filename,status='old',action='read')
      do i=1,10000
         read(20,'(A)',END=30) var
         if (var(1:1).eq.'#') then
            res=0d0
            rer=0d0
            do j=1,200
               read(20,'(A)') line
               if (line.eq."") goto 10
               read(line,'(4(1X,D14.8))') a,b,val,err
               res=res+val*(b-a)
               rer=rer+(err*(b-a))**2
            enddo
 10         write(6,*) var//'int.',res,'+/-',sqrt(rer)
         endif
      enddo      
 30   close(20)
      stop
      end
      
    
      
