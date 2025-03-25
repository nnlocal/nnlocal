program merge_gnuplot_data
  implicit none
  integer, parameter :: maxfiles = 500, maxlines = 25000
  character(len=100) :: files(maxfiles)
  character(len=100) :: line(maxlines, maxfiles)
  integer :: nlines(maxfiles)
  integer :: ifile, nfiles, ios, k, imethod
  character(len=1) :: cmethod
  real(kind=8) :: v1, v2, v3, v4, y, err
  integer :: ilength
  external :: ilength

  integer :: j, nfilest
  real(8), allocatable :: upval(:,:), uperr(:,:)
  real(8), allocatable :: dnval(:,:), dnerr(:,:)
  real(kind=8) :: tmpval, tmperr
  integer :: maxupt, maxdnt, maxuptrim, maxdntrim

  ! - kh modification started here >>>>>>>
  imethod = -1
  call getarg(1, cmethod)
  if (cmethod == '1') then
     imethod = 1
  elseif (cmethod == '2') then
     imethod = 2
  elseif (cmethod == '3') then
     imethod = 3
  elseif (cmethod == '4') then
     imethod = 4
  elseif (cmethod == '5') then
     imethod = 5
  elseif (cmethod == ' ') then
     imethod = 0 ! If nothing follows ./mergedata.exe on the command line this setting is acquired which then steers the code to read in the combination mode and input files from the read command prompts as it had been doing.
  else
     write(6,*) 'Combination mode must be "1-5" or " " : ', cmethod
     write(6,*) 'Quitting ...'
     stop
  endif

  if (imethod /= 0) then
     do ifile = 1, maxfiles
        call getarg(ifile + 1, files(ifile))
        if (trim(files(ifile)) == '') then
           nfiles = ifile - 1
           write(6,*) 'mergedata.exe found ', nfiles, ' files on the command line ...'
           goto 9
        endif
     enddo
  endif

9 continue
  if (imethod == 0) then
    write(*,*) ' enter 1 for combining sets with equal statistics'
    write(*,*) ' 2 to combine uneven sets'
    write(*,*) ' 3 to add sets (like born+virtual+real ... etc)'
    write(*,*) ' 4 to get maximum'
    write(*,*) ' 5 to get minimum'
    read(*,*) imethod, maxuptrim, maxdntrim
    write(*,*) ' enter files'
    do ifile = 1, maxfiles
       read(*,'(a)') files(ifile)
       if (files(ifile) == ' ') then
          nfiles = ifile - 1
          goto 10
       endif
    enddo
    write(*,*) ' too many files, increase maxfiles'
    call exit(-1)
  endif

10 continue
  ! set trim length
  maxupt = 0
  if (maxuptrim > 0 .and. maxuptrim < nfiles) maxupt = maxuptrim
  maxdnt = 0
  if (maxdntrim > 0 .and. maxdntrim < nfiles) maxdnt = maxdntrim
  if (maxupt + maxdnt > 0) write(*,*) maxupt + maxdnt, ' entries removed from each bin'

  allocate(upval(maxupt, maxlines))
  allocate(uperr(maxupt, maxlines))
  allocate(dnval(maxdnt, maxlines))
  allocate(dnerr(maxdnt, maxlines))
  
  ! load data
  do ifile = 1, nfiles
     open(unit=11, file=files(ifile), status='old')
     do k = 1, maxlines + 1
        read(unit=11, fmt='(a)', end=111) line(k, ifile)
        if (k == maxlines + 1) then
           write(*,*) ' too many lines in file, increase maxlines'
           call exit(-1)
        endif
        goto 12
111    nlines(ifile) = k - 1
        goto 11
12 continue
     enddo
11 continue
  enddo

  do ifile = 1, nfiles
     if (nlines(ifile) /= nlines(1)) then
        write(*,*) ' error: file', files(ifile), ' does not match in length'
        call exit(-1)
     endif
  enddo

  do k = 1, nlines(1)
     read(unit=line(k,1), fmt=*, iostat=ios) v1, v2, v3, v4
     if (ios /= 0) then
        write(12,'(a)') line(k, 1)(1:ilength(line(k,1)))
     else
        if (imethod == 1) then
           y = v3
           err = v4**2
           upval(1, k) = y
           uperr(1, k) = err
           do j = 2, maxdnt
              dnval(j, k) = 1d100
           enddo
           dnval(1, k) = y
           dnerr(1, k) = err
        elseif (imethod == 2) then
           if (v4 /= 0) then
              y = v3 / v4**2
              err = 1 / v4**2
           else
              y = 0
              err = 0
           endif
        elseif (imethod == 3) then
           y = v3
           err = v4**2
        elseif (imethod == 4) then
           y = v3
           err = v4**2
        elseif (imethod == 5) then
           y = v3
           err = v4**2
        endif

        do ifile = 2, nfiles
           read(unit=line(k, ifile), fmt=*, iostat=ios) v1, v2, v3, v4
           if (imethod == 1 .or. imethod == 3) then
              if (v3 > upval(maxupt, k)) then
                 upval(maxupt, k) = v3
                 uperr(maxupt, k) = v4**2
                 do j = maxupt-1, 1, -1
                    if (upval(j, k) < upval(j+1, k)) then
                       tmpval = upval(j, k)
                       tmperr = uperr(j, k)
                       upval(j, k) = upval(j+1, k)
                       uperr(j, k) = uperr(j+1, k)
                       upval(j+1, k) = tmpval
                       uperr(j+1, k) = tmperr
                    endif
                 enddo
              endif
              if (v3 < dnval(maxdnt, k)) then
                 dnval(maxdnt, k) = v3
                 dnerr(maxdnt, k) = v4**2
                 do j = maxdnt-1, 1, -1
                    if (dnval(j, k) > dnval(j+1, k)) then
                       tmpval = dnval(j, k)
                       tmperr = dnerr(j, k)
                       dnval(j, k) = dnval(j+1, k)
                       dnerr(j, k) = dnerr(j+1, k)
                       dnval(j+1, k) = tmpval
                       dnerr(j+1, k) = tmperr
                    endif
                 enddo
              endif
              y = y + v3
              err = err + v4**2
           elseif (imethod == 2) then
              if (v4 /= 0) then
                 y = y + v3 / v4**2
                 err = err + 1 / v4**2
              endif
           elseif (imethod == 3) then
              y = y + v3
              err = err + v4**2
           elseif (imethod == 4) then
              if (v3 > y) then
                 y = v3
                 err = v4**2
              endif
           elseif (imethod == 5) then
              if (v3 < y) then
                 y = v3
                 err = v4**2
              endif
           endif
        enddo

        if (imethod == 1) then
           y = y / nfiles
           err = sqrt(err / nfiles**2)
        elseif (imethod == 2) then
           if (err /= 0) then
              y = y / err
              err = 1 / sqrt(err)
           else
              y = 0
              err = 0
           endif
        elseif (imethod >= 3) then
           err = sqrt(err)
        endif
        write(12, '(4(1x, d14.8))') v1, v2, y, err
     endif
  enddo

  if (maxupt + maxdnt == 0) return
  ! --- trimming
  rewind 12
  do k = 1, nlines(1)
     if (trim(adjustl(line(k, 1))) == '') then
        write(13, *)
        cycle
     endif
     read(12, fmt=*, iostat=ios) v1, v2, v3, v4
     if (ios /= 0) then
        write(13, '(a)') line(k, 1)(1:ilength(line(k, 1)))
     else
        if (imethod == 1) then
           y = v3 * nfiles
           err = (v4 * nfiles)**2
           do j = 1, maxupt
              y = y - upval(j, k)
              err = err - uperr(j, k)
           enddo
           nfilest = nfiles - maxupt
           if (dnval(1, k) < 1d100) then
              do j = 1, maxdnt
                 y = y - dnval(j, k)
                 err = err - dnerr(j, k)
              enddo
              nfilest = nfilest - maxdnt
           endif
           y = y / nfilest
           if (err < 0) err = 0
           err = sqrt(err / nfilest**2)
           write(13, '(4(1x, d14.8))') v1, v2, y, err
        endif
     endif
  enddo

end program merge_gnuplot_data


function ilength(line)
  integer :: ilength
  character(len=*) :: line
  ilength = len(line)
  do j = ilength, 1, -1
     if (line(j:j) /= ' ') then
        ilength = j
        return
     endif
  enddo
  ilength = 0
end function ilength



