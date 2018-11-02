program velocity
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! velocity(inname, outname1, outname2, len,  p.a., xmove, ymove) !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  !                                                  !
  ! Purpose of this code: to give an output file of  !
  ! rotational velocities for a galaxy per fiber.    !
  ! This code requires an input file of wavelengths  !
  ! of a given line from reduced multi-fiber spectra.!
  ! Need to determine the helio-centric correction   !
  ! before-hand to put it into the code.             !
  !                                                  !
  ! After correcting for any differences between the !
  ! measured skyline and the value from Osterbrok    !
  ! '98, the 'redshift' per fiber is found in order  !
  ! to convert the wavelengths to velocities using   !
  ! z = v/c where c is the speed of light, z is the  !
  ! redshift and v is the measured velocity.         !
  !                                                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: pi = 3.1415926535

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Number of fibers on spectrograph. 82 for sparse- !
  ! pak.                                             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: fibnum = 82

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The helio-centric velocity, Calculate beforehand !
  ! using 'rvcorrect' in IRAF. In km/s.              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: helio = -13.14

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sky-lines froms Osterbrok. First one is for !
  ! Ha and NII. Second is for the S lines.      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real, parameter :: sky1 = 6577.284
  real, parameter :: sky2 = 6834.22

  real, allocatable :: wavein(:), skyob(:), fiber(:), skydif(:), waver(:), z(:), vel(:), veld(:), x3(:), y3(:), rwave(:), rsky(:), wavenum(:), avgfiber(:)

  integer, allocatable :: rwaved(:)

  integer len

  real loop, tol

  real sum

  real x(fibnum), y(fibnum)
  real x2(fibnum), y2(fibnum)

  real velha(fibnum), veln(fibnum), vels1(fibnum), vels2(fibnum)
  real avgvel(fibnum), avgerr(fibnum), errvec(4)
  real errha(fibnum), errn(fibnum), errs1(fibnum), errs2(fibnum)

  real xmove, ymove

  real posang, posangrad
  integer i, j, k, m, l, dum

  real spfiber(fibnum)

  character inname*40
  character outname1*40
  character outname2*40
  character pa*40
  character xdel*40
  character ydel*40
  character dimen*40

  call getarg(1, inname)
  call getarg(2, outname1)
  call getarg(3, outname2)
  call getarg(4, dimen)
  call getarg(5, pa)
  call getarg(6, xdel)
  call getarg(7, ydel)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Converts the character input arguments into their !
  ! respective types. Remember: number before decimal !
  ! is total number of numbers, number after decimal  !
  ! is how many numbers after decimal.                !
  ! e.g.                                              !
  ! 0 '(f5.1)' = 0000.0                               !
  ! 62 '(f5.1)' = 0006.2                              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  read (pa, '(f4.0)') posang
  read (xdel, '(f4.0)') xmove
  read (ydel, '(f4.0)') ymove
  read (dimen, '(i4)') len

  write(6, *) xmove, ymove

  write(6, *) len, posang

  allocate(wavein(len), skyob(len), fiber(len), skydif(len), waver(len), z(len), vel(len), veld(len), x3(len), y3(len), rwave(len), rsky(len), wavenum(len), avgfiber(len), rwaved(len))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Turn posang into radians !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  posangrad = -posang*pi/180

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Determine what spectral line the code is dealing with  !
  ! per fiber.                                             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(42, file = inname, status = 'old')

  i = 1
  do i = 1, len
     read (42, *) wavein(i), skyob(i), fiber(i), wavenum(i)
     if(wavenum(i) == 1) then
        rwave(i) = halpha
        rsky(i) = sky1
     else if(wavenum(i) == 2) then
        rwave(i) = nitrogen
        rsky(i) = sky1
     else if (wavenum(i) == 3) then
        rwave(i) = sulfur1
        rsky(i) = sky2
     else if(wavenum(i) == 4) then
        rwave(i) = sulfur2
        rsky(i) = sky2
     else
        rwave(i) = 1
     endif
  enddo
 
  open(3, file = outname1, status = 'unknown')
  open(4, file = 'sparsepakcoord.dat', status = 'old')
  open(5, file = outname2, status = 'unknown')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read in the sparsepak coordinates. x and y !
  ! in arcsecs.                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  j = 1
  do j = 1, fibnum
     read(4, *) spfiber(j), x(j), y(j)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Rotates the fiber positions with the given position !
  ! angle.                                              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  x2 = -x*cos(posangrad) + xmove + y*sin(posangrad)
  y2 = x*sin(posangrad) + y*cos(posangrad) + ymove 
  
  !y2 = (x + xmove)*sin(posangrad) + (y + ymove)*cos(posangrad)
  !x2 = (-x + xmove)*cos(posangrad) + (y + ymove)*sin(posangrad)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Taking into account any differences between the sky !
  ! line measured and what it is supposed to be. Stored !
  ! as 'skydif'. Should be a positive value if line is  !
  ! redder than lab value, negative if bluer. 'skydif'  !
  ! then added to 'wavein' to give 'waver', the         !
  ! corrected wavelength.                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  skydif = rsky - skyob

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
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create velocities for each element !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  i = 1
  m = 1
  l = 1
30 do i = l, fibnum
800  do j = m, len
        dum = 1
        if(wavenum(m) .eq. 0) then
           velha(l) = 0
           veln(l) = 0
           vels1(l) = 0
           vels2(l) = 0
           m = m + 1
           l = l + 1
           goto 800

        !!!!!!!!!!!!!!!!!!!!!!
        ! If Halpha is first !
        !!!!!!!!!!!!!!!!!!!!!!
        elseif(wavenum(m) .eq. 1) then
           velha(l) = vel(m)
           m = m + 1
           if(fiber(m) .eq. fiber(m - 1)) then
              if(wavenum(m) .eq. 2) then
                 veln(l) = vel(m)
                 m = m + 1
              elseif(wavenum(m) .eq. 3) then
                 veln(l) = 0
                 vels1(l) = vel(m)
                 m = m + 1
                 if(fiber(m) .eq. fiber(m - 1)) then
                    if(wavenum(m) .eq. 4) then
                       vels2(l) = vel(m)
                       m = m + 1
                       goto 500
                    else
                       vels2(l) = 0
                       goto 500
                    endif
                 else
                    goto 500
                 endif
              elseif(wavenum(m) .eq. 4) then
                 veln(l) = 0
                 vels1(l) = 0
                 vels2(l) = vel(m)
                 m = m + 1
                 goto 500
              else
                 veln(l) = 0
                 vels1(l) = 0
                 vels2(l) = 0 
                 goto 500
              endif
              if(fiber(m) .eq. fiber(m - 1)) then
                 if(wavenum(m) .eq. 3) then
                    vels1(l) = vel(m)
                    m = m + 1
                 elseif(wavenum(m) .eq. 4) then
                    vels1(l) = 0
                    vels2(l) = vel(m)
                    m = m + 1
                    goto 500
                 else
                    vels1(l) = 0
                    vels2(l) = 0
                    goto 500
                 endif
              else
                 goto 500
              endif
              if(fiber(m) .eq. fiber(m -1)) then
                 if(wavenum(m) .eq. 4) then
                    vels2(l) = vel(m)
                    m = m + 1
                    goto 500
                 else
                    vels2(l) = 0
                    goto 500
                 endif
              else
                 goto 500
              endif
           else
              goto 500
           endif

        !!!!!!!!!!!!!!!!!!!
        ! If NII is first !
        !!!!!!!!!!!!!!!!!!!
        elseif(wavenum(m) .eq. 2) then
           velha(l) = 0
           veln(l) = vel(m)
           m = m + 1
           if(fiber(m) .eq. fiber(m - 1)) then
              if(wavenum(m) .eq. 3) then
                 vels1(l) = vel(m)
                 m = m + 1
                 if(fiber(m) .eq. fiber(m - 1)) then
                    if(wavenum(m) .eq. 4) then
                       vels2(l) = vel(m)
                       m = m + 1
                       goto 500
                    else
                       vels2(l) = 0
                       goto 500
                    endif
                 else
                    goto 500
                 endif
              elseif(wavenum(m) .eq. 4) then
                 vels1(l) = 0
                 vels2(l) = vel(m)
                 m = m + 1
                 goto 500
              else
                 vels1(l) = 0
                 vels2(l) = 0
                 goto 500
              endif
              if(fiber(m) .eq. fiber(m - 1)) then
                 if(wavenum(m) .eq. 4) then
                    vels2(l) = vel(m)
                    m = m + 1
                    goto 500
                 else
                    vels2(l) = 0
                    goto 500
                 endif
              else
                 goto 500
              endif
           else
              goto 500
           endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! If the first S line is first !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif(wavenum(m) .eq. 3) then
           velha(l) = 0
           veln(l) = 0
           vels1(l) = vel(m)
           m = m + 1
           if(fiber(m) .eq. fiber(m - 1)) then
              if(wavenum(m) .eq. 4) then
                 vels2(l) = vel(m)
                 m = m + 1
                 goto 500
              else
                 vels2(l) = 0
                 goto 500
              endif
           else
              goto 500
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
        endif
500     l = l + 1
        goto 30
     enddo
  enddo

  do i = 1, fibnum
     write(5, *) spfiber(i), velha(i), veln(i), vels1(i), vels2(i)
     !write(6, *) spfiber(i), velha(i), veln(i), vels1(i), vels2(i)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Determine the average velocity per fiber         ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  sum = 0
  m = 1
  do i = 1, fibnum
     m = 1
     sum = 0
     do j = m, len
         k = 0
100      if(fiber(m) .eq. spfiber(i)) then
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
   enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Computes errors in velocities by subtracting the   !
  ! element velocities from average. If only Ha then   !
  ! the error is 10 km/s                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  i = 1
  do i = 1, fibnum
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
  j = 1
  do j = 1, fibnum
      errvec = (/ errha(j), errn(j), errs1(j), errs2(j) /)
      avgerr(j) = maxval(errvec)
      if(avgerr(j) .eq. 0 .and. avgvel(j) .ne. 0) then
         avgerr(j) = 10
      endif
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write out the velocities and fibers with corrected !
  ! x and y values and errors.                         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  i = 1
  do i = 1, fibnum
     write(3, *) x2(i), y2(i), avgvel(i), avgerr(i)
     !write(6, *) spfiber(i), avgvel(i), x2(i), y2(i)
  enddo

end program
