module Dyn_grids
  double precision, DIMENSION(:, :, :), ALLOCATABLE :: grids  
  double precision, DIMENSION(:, :), ALLOCATABLE :: xi
endmodule Dyn_grids

subroutine cmpnewgridpars(ndim,igridch,part,switch)
  use Dyn_grids
  implicit none
  include 'parallel.f'
  integer ndim,i,j,k,nd,jj
  integer jdim
  common/csidim/jdim
  double precision cum
  common/csigri/cum
  character *32 ingridfile
  character *32 outgridfile
  character *4 ks,nfilchar,part
  character *1 igridch,im1ch
  integer xgn,ll,switch,igrid
  double precision err,suminterp,t,xmin,xmax
  external suminterp
  nd=50
  ALLOCATE (grids(1:nfile, 0:nd,1:ndim))
  ALLOCATE (xi(1:nd,1:ndim))
  ll=len_trim(igridch)
  read(igridch,'(I4)') igrid
  if (switch.eq.0) write(im1ch,'(i1)') igrid-1
  if (switch.eq.1) write(im1ch,'(i1)') igrid
  do k = 1,nfile
     write(ks,'(i4.4)') k
     ingridfile='xg'//im1ch(1:ll)//'-'//ks//'-'//part//'.grid'
     open(unit=100+k,file=ingridfile,status='unknown')
     do j=1,ndim
        grids(k,0,j)=0d0
        read(100+k,203,end=10) jj,(grids(k,i,j),i=1,nd)
     enddo
10   close(100+k)
  enddo
  if (jj.ne.ndim) then
     write(6,*) 'compnewgrid: grid not sound, exiting...'
     write(6,*) 'jj,ndim ',jj,ndim
     stop
  end if
  do jdim = 1,jj
     cum=0d0
     do i = 1,nd
        cum=dfloat(i*nfile)/dfloat(nd)
        call dzero(0d0,1d0,t,err,1d-10,10000000,suminterp)
        xi(i,jdim)=t
     enddo
  enddo
  if     (switch .eq. 0) then
     outgridfile='cfmc-'//seedchar//'-'//'xg'//igridch//'-'//''//part//'.grid'
  elseif (switch .eq. 1) then
     outgridfile='cfmc-'//seedchar//'-'//part//'.grid'
  endif
  open(unit=999,file=outgridfile,status='unknown')
  do j=1,jj
     write(999,203) j,(xi(i,j),i=1,nd)
  enddo
  close(999)      
DEALLOCATE (grids)
DEALLOCATE (xi)
203 FORMAT(/(5z16))
return
end subroutine cmpnewgridpars

!c---  Interpolation of the sum of all cumulative distributions
!c---  from all seed minus current cumulant
!c---  The integer specifying the variable is passed via common
double precision function suminterp(x,l)
  use Dyn_grids
  implicit none
  include 'parallel.f'
  integer l
  double precision x
  integer i,j,k,nd
  double precision amb,xmb,amx
  integer jdim
  common/csidim/jdim
  double precision cum
  common/csigri/cum  
  if(x.eq.1d0) then
     suminterp=nfile-cum
     return
  endif
  nd=50
  suminterp=0
  do i=1,nfile
     do j=1,nd
        if (grids(i,j,jdim).gt.x) then
           amb=grids(i,j-1,jdim)-grids(i,j,jdim)
           xmb=x                -grids(i,j,jdim)
           amx=grids(i,j-1,jdim)-x
           suminterp=suminterp &
 &         +xmb*(j-1)/amb/dfloat(nd) &
 &         +amx*j    /amb/dfloat(nd)
           exit
        endif
     enddo
  enddo
  suminterp = suminterp - cum
  return
end function suminterp








