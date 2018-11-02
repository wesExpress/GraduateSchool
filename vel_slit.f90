program velocity
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! velocity(inname, outname1, outname2, len) !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                  !
  ! 'inname' will be name of file which has been     !
  ! formatted by hand and contains the wavelengths   !
  ! for the various lines per fiber.                 !
  !                                                  !
  ! 'outname' will be name of output file.           !
  !                                                  !
  ! The data file will contain a column to tell the  !
  ! code which spectral line each line is:           !
  !    0 : no line                                   !
  !    1 : H-alpha                                   !
  !    2 : N-II                                      !
  !    3 : S-II (first)                              !
  !    4 : S-II (second)                             !
  ! This will tell the code which sky-line to use    !
  ! as well.                                         !
  !                                                  !
  ! p.a. is the position angle of the observation.   !
  ! This will be in the log file of your observation.!
  !                                                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: pi = 3.1415926535

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The helio-centric velocity, Calculate beforehand !
  ! using 'rvcorrect' in IRAF. In km/s.              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: helio = -13.29

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Speed of light in km/s !
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: c = 299792

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Rest wavelengths in Angstroms !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: halpha = 6562.8
  real, parameter :: nitrogen = 6583.39
  real, parameter :: sulfur1 = 6716.47
  real, parameter :: sulfur2 = 6730.85

  real, parameter :: ox = 5006.843
  real, parameter :: hbeta = 4861.363

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sky-lines froms Osterbrok.                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real, parameter :: sky1 = 6300.304
  real, parameter :: sky2 = 0

  real, allocatable :: wavein(:), skyob(:), lines(:), skydif(:), waver(:), z(:), vel(:), veld(:), rwave(:), rsky(:), wavenum(:), velha(:), veln(:), vels1(:), vels2(:), line(:), avgvel(:), errha(:), errn(:), errs1(:), errs2(:), avgerr(:)

  real, allocatable :: velhbeta(:), velox(:), errhbeta(:), errox(:), errvec(:)

  integer, allocatable :: rwaved(:)

  integer len, lnum, mode

  real loop, tol

  real sum

  integer i, j, k, m, l, dum


  character inname*40
  character outname1*40
  character outname2*40
  character dimen*40
  character moded*40

  call getarg(1, inname)
  call getarg(2, outname1)
  call getarg(3, outname2)
  call getarg(4, dimen)
  call getarg(5, moded)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Converts the character input arguments into their !
  ! respective types. Remember: number before decimal !
  ! is total number of numbers, number after decimal  !
  ! is how many numbers after decimal.                !
  ! e.g.                                              !
  ! 0 '(f5.1)' = 0000.0                               !
  ! 62 '(f5.1)' = 0006.2                              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  read (dimen, '(i4)') len
  read (moded, '(i4)') mode

  if(mode .eq. 1) then
     allocate(errvec(4))
  elseif(mode .eq. 2) then
     allocate(errvec(2))
  endif

  write(6, *) len, mode

  allocate(wavein(len), skyob(len), skydif(len), waver(len), z(len), vel(len), veld(len), rwave(len), rsky(len), wavenum(len), rwaved(len), lines(len))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Determine what spectral line the code is dealing with  !
  ! per fiber.                                             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(42, file = inname, status = 'old')

  i = 1
  do i = 1, len
     read (42, *) wavein(i), skyob(i), lines(i), wavenum(i)
     !write(6, *) wavein(i), skyob(i), lines(i), wavenum(i)
  enddo

  j = 1
  k = 0
  do j = 1, len - 1
     if(lines(j) .ne. lines(j + 1)) then
        k = k + 1
     endif
  enddo
  lnum = k
  write(6, *) lnum

  do i = 1, len
     if(wavenum(i) .eq. 1) then
        rwave(i) = halpha
     elseif(wavenum(i) .eq. 2) then
        rwave(i) = nitrogen
     elseif(wavenum(i) .eq. 3) then
        rwave(i) = sulfur1
     elseif(wavenum(i) .eq. 4) then
        rwave(i) = sulfur2
     elseif(wavenum(i) .eq. 5) then
        rwave(i) = hbeta
     elseif(wavenum(i) .eq. 6) then
        rwave(i) = ox
     endif
  enddo

  allocate(velha(lnum), veln(lnum), vels1(lnum), vels2(lnum), line(lnum), avgvel(lnum), errha(lnum), errn(lnum), errs1(lnum), errs2(lnum), avgerr(lnum), velhbeta(lnum), velox(lnum), errhbeta(lnum), errox(lnum))
  
  k = 1
  m = 1
  do k = 1, lnum
10  do j = m, len - 1
        if(lines(m) .eq. lines(m + 1)) then
           m = m + 1
           goto 10
        elseif(lines(m) .ne. lines(m + 1)) then
           line(k) = lines(m)
        endif
        m = m + 1
        goto 400
     enddo
400  dum = 1
  enddo
 
  open(3, file = outname1, status = 'unknown')
  open(5, file = outname2, status = 'unknown')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Taking into account any differences between the sky !
  ! line measured and what it is supposed to be. Stored !
  ! as 'skydif'. Should be a positive value if line is  !
  ! redder than lab value, negative if bluer. 'skydif'  !
  ! then added to 'wavein' to give 'waver', the         !
  ! corrected wavelength.                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(mode .eq. 1) then
     if(skyob(i) .ne. 0) then 
        skydif(i) = sky1 - skyob(i)
     else
        skydif = 0
     endif
  elseif(mode .eq. 2) then
     skydif = 0
  endif

  waver = wavein + skydif

  z = (waver - rwave)/rwave

  veld = z*c + helio

  i = 1
  do i = 1, len
     if(veld(i) .lt. 0) then
        vel(i) = 0
     else
        vel(i) = veld(i)
     endif
     !write(6, *) lines(i), vel(i)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create velocities for each element !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  i = 1
  m = 1
  l = 1
  if(mode .eq. 1) then
30   do i = l, lnum
        do j = m, len
           !!!!!!!!!!!!!!!!!!!!!!
           ! If Halpha is first !
           !!!!!!!!!!!!!!!!!!!!!!
           if(wavenum(m) .eq. 1) then
              velha(l) = vel(m)
              m = m + 1

              if(lines(m) .eq. lines(m - 1)) then
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Check to see if [NII] is next !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 if(wavenum(m) .eq. 2) then
                    veln(l) = vel(m)
                    m = m + 1

                    if(lines(m) .eq. lines(m - 1)) then
                       if(wavenum(m) .eq. 3) then
                          vels1(l) = vel(m)
                          m = m + 1

                          if(lines(m) .eq. lines(m - 1)) then
                             if(wavenum(m) .eq. 4) then
                                vels2(l) = vel(m)
                                m = m + 1
                                goto 50
                             else
                                vels2(l) = 0
                                goto 50
                             endif
                          else
                             goto 50
                          endif

                       elseif(wavenum(m) .eq. 4) then
                          vels1(l) = 0
                          vels2(l) = vel(m)
                          m = m + 1
                          goto 50
                       else
                          vels1(l) = 0
                          vels2(l) = 0
                          goto 50
                       endif

                    else
                       goto 50
                    endif

                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Or if first [SII] is next !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 elseif(wavenum(m) .eq. 3) then
                    veln(l) = 0
                    vels1(l) = vel(m)
                    m = m + 1

                    if(lines(m) .eq. lines(m - 1)) then
                       if(wavenum(m) .eq. 4) then
                          vels2(l) = vel(m)
                          m = m + 1
                          goto 50
                       else
                          vels2(l) = 0
                          goto 50
                       endif

                    else
                       goto 50
                    endif

                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Or if second [SII] is next !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 elseif(wavenum(m) .eq. 4) then
                    veln(l) = 0
                    vels1(l) = 0
                    vels2(l) = vel(m)
                    m = m + 1
                    goto 50
                 endif

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! If there is nothing after this Halpha !
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              else
                 veln(l) = 0
                 vels1(l) = 0
                 vels2(l) = 0
                 goto 50
              endif

           !!!!!!!!!!!!!!!!!!!
           ! If NII is first !
           !!!!!!!!!!!!!!!!!!!
           elseif(wavenum(m) .eq. 2) then
              velha(l) = 0
              veln(l) = vel(m)
              m = m + 1

              if(lines(m) .eq. lines(m - 1)) then
                 if(wavenum(m) .eq. 3) then
                    vels1(l) = vel(m)
                    m = m + 1

                    if(lines(m) .eq. lines(m - 1)) then
                       if(wavenum(m) .eq. 4) then
                          vels2(l) = vel(m)
                          m = m + 1
                          goto 50
                       else
                          vels2(l) = 0
                          goto 50
                       endif
                    else
                       goto 50
                    endif

                 elseif(wavenum(m) .eq. 4) then
                    vels1(l) = 0
                    vels2(l) = vel(m)
                    m = m + 1
                    goto 50
                 else
                    vels1(l) = 0
                    vels2(l) = 0
                    goto 50
                 endif

                 if(lines(m) .eq. lines(m - 1)) then
                    if(wavenum(m) .eq. 4) then
                       vels2(l) = vel(m)
                       m = m + 1
                       goto 50
                    else
                       vels2(l) = 0
                       goto 50
                    endif
                 else
                    goto 50
                 endif

              else
                 vels1(l) = 0
                 vels2(l) = 0
                 goto 50
              endif

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! If the first S line is first !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           elseif(wavenum(m) .eq. 3) then
              velha(l) = 0
              veln(l) = 0
              vels1(l) = vel(m)
              m = m + 1

              if(lines(m) .eq. lines(m - 1)) then
                 if(wavenum(m) .eq. 4) then
                    vels2(l) = vel(m)
                    m = m + 1
                    goto 50
                 else
                    vels2(l) = 0
                    goto 50
                 endif
              else
                 vels2(l) = 0
                 goto 50
              endif

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! If the second S line is first !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           elseif(wavenum(m) .eq. 4) then
              velha(l) = 0
              veln(l) = 0
              vels1(l) = 0
              vels2(l) = vel(m)
              m = m + 1
              goto 50
           endif

50         l = l + 1
           goto 30
        enddo
     enddo

  elseif(mode .eq. 2) then
35   do i = l, lnum
     do j = m, len
           dum = 1

           !!!!!!!!!!!!!!!!!!!!!
           ! If HBeta is first !
           !!!!!!!!!!!!!!!!!!!!!
           if(wavenum(m) .eq. 5) then
              velox(l) = 0
              velhbeta(l) = vel(m)
              m = m + 1
              if(lines(m) .eq. lines(m-1)) then
                 if(wavenum(m) .eq. 6) then
                    velox(l) = vel(m)
                    m = m + 1
                    goto 500
                 else
                    velox(l) = 0
                    goto 500
                 endif
              else
                 goto 500
              endif

           !!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! If [OIII] 5007 is first !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!
           elseif(wavenum(m) .eq. 6) then
              velhbeta(l) = 0
              velox(l) = vel(m)
              m = m + 1
              goto 500
           endif
500        l = l + 1   
           goto 35
        enddo
     enddo
  endif

  write(6, *) l

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Determine the average velocity per fiber/line    ! 
  !                                                  !
  ! This is my crowning achievement as a programmer. !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  sum = 0
  m = 1
  do i = 1, lnum
     m = 1
     sum = 0
     do j = m, len
         k = 0
100      if(lines(m) .eq. line(i)) then
            sum = sum + vel(m)
            k = k + 1
            m = m + 1
            goto 100
         endif
         if(k .eq. 0) then
            k = 1
            m = m + 1
         endif
         sum = sum/k
         avgvel(i) = sum
         goto 200
200   enddo
      !write(6, *) line(i), avgvel(i)
   enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Computes errors in velocities by subtracting the   !
  ! element velocities from average. If only Ha then   !
  ! the error is 10 km/s                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   i = 1
   if(mode .eq. 1) then
      do i = 1, lnum
         if(velha(i) .ne. 0) then
            errha(i) = abs(avgvel(i) - velha(i))
         else
            errha(i) = 0
         endif
         if(veln(i) .ne. 0) then
            errn(i) = abs(avgvel(i) - veln(i))
         else
            errn(i) = 0
         endif
         if(vels1(i) .ne. 0) then
            errs1(i) = abs(avgvel(i) - vels1(i))
         else
            errs1(i) = 0
         endif
         if(vels2(i) .ne. 0) then
            errs2(i) = abs(avgvel(i) - vels2(i))
         else
            errs2(i) = 0
         endif
      enddo
   elseif(mode .eq. 2) then
      do j = 1, lnum
         if(velhbeta(j) .ne. 0) then
            errhbeta(j) = abs(avgvel(j) - velhbeta(j))
         else
            errhbeta(j) = 0
         endif
         if(velox(j) .ne. 0) then
            errox(j) = abs(avgvel(j) - velox(j))
         else
            errox(j) = 0
         endif
         !write(6, *) errhbeta(j), velox(j)
      enddo
   endif
   
   j = 1
   do j = 1, lnum
      if(mode .eq. 1) then
         errvec = (/ errha(j), errn(j), errs1(j), errs2(j) /)
         avgerr(j) = maxval(errvec)
         if(avgerr(j) .eq. 0 .and. avgvel(j) .ne. 0) then
            avgerr(j) = 10
         endif
      elseif(mode .eq. 2) then
         errvec = (/ errhbeta(j), errox(j) /)
         avgerr(j) = maxval(errvec)
         if(avgerr(j) .eq. 0 .and. avgvel(j) .ne. 0) then
            avgerr(j) = 10
         endif
      endif
   enddo
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write out the velocities                           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = 1, lnum
     write(3, *) line(i), avgvel(i), avgerr(i)
     if(mode .eq. 1) then
        write(5, *) line(i), velha(i), veln(i), vels1(i), vels2(i)
     elseif(mode .eq. 2) then
        write(5, *) line(i), velhbeta(i), velox(i)
     endif
  enddo

  
end program
