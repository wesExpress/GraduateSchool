program az_prof
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program takes in the output file from the  !
  ! elliptical panda region in ds9. This region     !
  ! creates azimuthal profiles and outputs a file   !
  ! that proceeds azimuthally moving out. In other  !
  ! words, if there are 10 azimuthal sections and   !
  ! 10 radial sections, then the data file will be  !
  ! 100 lines long, with every 10 lines composing   !
  ! one radial slice.                               !
  !                                                 !
  ! This code creates an azimuthal profile for each !
  ! radial slice. The azimuthal profile is output   !
  ! as columns of radii, with each row being an     !
  ! azimuthal slice.                                !
  !                                                 !
  ! The phase angle is also calculated by taking    !
  ! the Fourier Transform of the azimuthal profiles !
  ! following the method in Puerari & Dottori '97.  !
  !                                                 !
  ! The inputs for this code are:                   !
  !  1. inname  => input file                       !
  !  2. outname1 => output file for all profiles    !
  !  3. outname2 => averaged profile                !
  !  4. raddiv  => number of radial divisions       !
  !  5. angdiv  => number of azimuthal divisions    !
  !  6. rstart
  !  7. rend
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: pi = 3.1415926535

  real, allocatable :: reg(:), sumf(:), sumf_err(:), area(:), surf(:), surf_err(:)
  real, allocatable :: prof(:,:), phase(:), phi(:), r(:), amp(:,:), a_coeff(:,:), b_coeff(:,:)
  real, allocatable :: phi2(:), a_o(:), prof2(:,:), amp_rat(:), new_phase(:)
  real, allocatable :: i_bar(:), i_interbar(:)

  complex, allocatable :: fourier(:)
  
  complex img, fouriersum

  integer raddiv, angdiv
  integer n, dum1, dum2
  integer i, j, k
  integer xpt, r_start, r_end
  
  real fsum, top, dtheta, asum, bsum, a0sum

  character inname*32, outname1*32, outname2*32, rad*32, ang*32
  character zpd*32, xptd*32, pixd*32, incd*32
  character r_startd*32, r_endd*32

  !!!!!!!!!!!!!!!!!!!!!
  ! take in arguments !
  ! allocate arrays   !
  !!!!!!!!!!!!!!!!!!!!!

  call getarg(1, inname)
  call getarg(2, outname1)
  call getarg(3, outname2)
  call getarg(4, rad)
  call getarg(5, ang)
  call getarg(6, r_startd)
  call getarg(7, r_endd)

  read(rad, '(i4)') raddiv
  read(ang, '(i4)') angdiv
  read(r_startd, '(i4)') r_start
  read(r_endd, '(i4)') r_end

  n = raddiv*angdiv

  dtheta = 360./angdiv*pi/180.

  img = (0.0,1.0)

  allocate(reg(n), sumf(n), sumf_err(n), area(n), surf(n), surf_err(n))
  allocate(prof(raddiv-1,angdiv), fourier(raddiv-1), phase(raddiv-1), phi(angdiv), r(raddiv-1))
  allocate(a_coeff(raddiv-1,7), b_coeff(raddiv-1,7), amp(raddiv-1,7), a_o(raddiv-1), phi2(angdiv))
  allocate(prof2(raddiv-1,angdiv), amp_rat(raddiv-1), new_phase(raddiv-1))
  allocate(i_bar(raddiv-1),i_interbar(raddiv-1))

  write(6, *) ''
  write(6, *) 'inname = ', inname
  write(6, *) 'outname1 = ', outname1
  write(6, *) 'outname2 = ', outname2
  write(6, *) 'raddiv = ', raddiv
  write(6, *) 'angdiv = ', angdiv
  write(6, *) 'n = ', n
  write(6, *) 'dtheta = ', dtheta
  write(6, *) 'r_start = ', r_start
  write(6, *) 'r_end = ', r_end
  write(6, *) ''

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! read in the data file !
  ! create phi array      !
  !!!!!!!!!!!!!!!!!!!!!!!!!

  open(42, file = inname, status = 'old')
  open(43, file = outname1, status = 'unknown')
  open(44, file = outname2, status = 'unknown')

  do i = 1, raddiv-1
     read(42, *) prof(i,:)
  enddo

  do i = 1, raddiv-1
     r(i) = r_start + real(r_end - r_start)/real(raddiv)*(i)
  enddo

  do i = 1, angdiv
     phi(i) = -pi + 2*pi/real(angdiv)*(i-1)
     phi2(i) = 0 + 2*pi/real(angdiv)*(i-1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write out the profiles !
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = 1, raddiv-1
     write(43, *) prof(i,:)
     !write(6, *) prof(i,:)
  enddo

  !prof = prof + 20

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Caluclate Fourier Transform and Phase !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fouriersum = (0.0,0.0)
  write(6, *) ''
  write(6, *) 'Computing Fourier Transform...'
  do i = 1, raddiv-1
     do j = 1, angdiv
        fouriersum = prof(i,j)*exp(-2.*img*phi(j)) + fouriersum
        if(i .eq. 1) then
           !write(6, *) prof(i,j), phi(j)
        endif
     enddo
     fourier(i) = fouriersum*dtheta

     fouriersum = (0.0,0.0)
  enddo

  write(6, *) 'Computing Phase...'
  phase = atan2(realpart(fourier),imagpart(fourier))
  phase = phase*180./pi
  new_phase = phase
  write(6, *) 'Done.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate Fourier Coefficients for m 1,2,3,4,5,6 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(6, *) ''
  write(6, *) 'Computing Fourier Amplitudes...'

  ! rearranges the profile array to be from 0 to 2pi !
  do i = 1, raddiv-1
     do j = 1, angdiv
        dum1 = int(j-angdiv/2.)
        dum2 = int(j+angdiv/2.)

        if(j .gt. angdiv/2.) then
           prof2(i,dum1) = prof(i,j)
        elseif(j .le. angdiv/2.) then
           prof2(i,dum2) = prof(i,j)
        endif
        
     enddo
  enddo
  !do i = 1, raddiv-1
  !   do j = 1, angdiv
  !      if(i .eq. 1) then
  !         write(6, *) prof(i,j), phi(j), prof2(i,j), phi2(j)
  !      endif
  !   enddo
  !enddo
   
  asum = 0
  bsum = 0

  ! calculates a and b coefficients and amplitude !
  do i = 1, 7
     do j = 1, raddiv-1
        if(i .gt. 1) then
           do k = 1, angdiv
              asum = asum + prof2(j,k)*cos((i - 1)*phi2(k))
              bsum = bsum + prof2(j,k)*sin((i - 1)*phi2(k))
           enddo
           a_coeff(j,i) = asum*dtheta/pi
           b_coeff(j,i) = bsum*dtheta/pi

           amp(j,i) = sqrt(a_coeff(j,i)**2. + b_coeff(j,i)**2.)
        elseif(i .eq. 1) then
           do k = 1, angdiv
              asum = asum + prof2(j,k)*cos((i-1)*phi2(k))
           enddo
           a_coeff(j,i) = asum*dtheta/(2*pi)

           amp(j,i) = a_coeff(j,i)
        endif

        asum = 0
        bsum = 0
     enddo
  enddo

  ! Determine bar and interbar intensities !

  i_bar = amp(:,1) + amp(:,3) + amp(:,5) + amp(:,7)
  i_interbar = amp(:,1) - amp(:,3) + amp(:,5) - amp(:,7)
  
  write(6, *) 'Done.'
  write(6, *) ''

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write out the phase profile and fourier information !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = 1, raddiv-1
     write(44, *) r(i), new_phase(i), amp(i,:), i_bar(i), i_interbar(i)
  enddo

end program az_prof
