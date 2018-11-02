program az_prof
  implicit none

  !!!!!!!!!!!!!!!!!!!
  ! Inputs:         !
  !   1. inname     !
  !   2. outname1   !
  !   3. outname2   !
  !   4. angdiv     !
  !   5. raddiv     !
  !   6. rstart     !
  !   7. rmax       !
  !   8. xcen       !
  !   9. ycen       !
  !  10. nx         !
  !  11. ny         !
  !!!!!!!!!!!!!!!!!!!

  real, parameter :: pi = 3.1415926535

  character angdivd*32, raddivd*32, rstartd*32, rmaxd*32
  character xcend*32, ycend*32, nxd*32, nyd*32
  character inname*32, outname1*32, outname2*32

  real, allocatable :: phi(:), phi_deg(:), rad(:)
  real, allocatable :: image_data(:,:)
  real, allocatable :: az_profs(:,:), rad_profs(:,:)
  real, allocatable :: rad_vals(:,:), phi_vals(:,:)

  real runavg, numpix
  real radcond, angcond
  real area, areacheck
  real rad_grid, phi_grid

  integer angdiv, raddiv, rstart, rmax, xcen, ycen, nx, ny
  integer xcen_new, ycen_new

  integer i, j, x, y

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read in the input parameters !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1, inname)
  call getarg(2, outname1)
  call getarg(3, outname2)
  call getarg(4, angdivd)
  call getarg(5, raddivd)
  call getarg(6, rstartd)
  call getarg(7, rmaxd)
  call getarg(8, xcend)
  call getarg(9, ycend)
  call getarg(10, nxd)
  call getarg(11, nyd)

  read(angdivd, '(i4)') angdiv
  read(raddivd, '(i4)') raddiv
  read(rstartd, '(i4)') rstart
  read(rmaxd, '(i4)') rmax
  read(xcend, '(i4)') xcen
  read(ycend, '(i4)') ycen
  read(nxd, '(i4)') nx
  read(nyd, '(i4)') ny

  !!!!!!!!!!!!!!!!!!!!!!!
  ! allocate the arrays !
  !!!!!!!!!!!!!!!!!!!!!!!

  allocate(phi(angdiv), phi_deg(angdiv), rad(raddiv), az_profs(raddiv,angdiv)) 
  allocate(image_data(ny,nx), rad_profs(angdiv,raddiv))
  allocate(rad_vals(ny,nx), phi_vals(ny,nx))

  !!!!!!!!!!!!!!!!!!!!!
  ! Read in the image !
  !!!!!!!!!!!!!!!!!!!!!

  open(42, file = inname, status = 'old')
  write(6, *) ''
  write(6, *) 'Reading in file ', inname
  do y = 1, ny
     read(42, *) image_data(y,:)
  enddo
  close(42, status = 'delete', iostat = I)
  write(6, *) 'Done reading in file ', inname
  write(6, *) ''

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find the true xcen and ycen !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(6, *) 'Computing true center...'
  xcen_new = xcen
  ycen_new = ycen
  do y = 1, ny
     do x = 1, nx
        if(x .ge. (xcen - 20) .and. x .le. (xcen + 20) .and. y .ge. (ycen - 20) .and. y .le. (ycen + 20)) then
           if(image_data(y,x) .gt. image_data(ycen_new,xcen_new)) then
              xcen_new = x
              ycen_new = y
           endif
        endif
     enddo
  enddo
  write(6, *) xcen, ycen
  write(6, *) xcen_new, ycen_new
  xcen = xcen_new
  ycen = ycen_new
  write(6, *) 'Done.'
  write(6, *) ''

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create the radial and angular arrays !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = 1, angdiv
     phi(i) = -pi + 2*pi/real(angdiv)*(i-1)
  enddo
  do j = 1, raddiv
     rad(j) = rstart + real(rmax-rstart)/real(raddiv)*(j-1)
  enddo
  phi_deg = phi*180./pi

  write(6, *) 'Creating rad and phi arrays...'
  do x = 1, nx
     do y = 1, ny
        rad_vals(y,x) = sqrt((x-xcen)**2. + (y-ycen)**2.)
        phi_vals(y,x) = atan2(real(y-ycen),real(x-xcen))
     enddo
  enddo
  write(6, *) 'Done.'
  write(6, *) ''

  rad_grid = real(rmax-rstart)/real(raddiv)
  phi_grid = 2*pi/real(angdiv)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This set of loops runs through the entire image and determines the azimuthal profiles !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  runavg = 0
  numpix = 0
  areacheck = 0
  write(6, *) 'Computing azimuthal profiles...'
  do i = 1, raddiv
     do j = 1, angdiv

        ! Looping through (x,y) !
        do x = 1, nx      
           do y = 1, ny           
              if((abs(rad_vals(y,x) - rad(i)) .le. rad_grid) .and. &
                   abs(phi_vals(y,x) - phi(j)) .le. 1.5*phi_grid) then
                 runavg = runavg + image_data(y,x)
                 numpix = numpix + 1
                 areacheck = areacheck + 1
              endif
           enddo
        enddo

        if(runavg .ne. 0) then
           az_profs(i,j) = runavg/numpix
           rad_profs(j,i) = runavg/numpix
        else
           az_profs(i,j) = 0
           rad_profs(j,i) = 0
           write(6, *) rad(i), phi(j)
        endif

        runavg = 0 
        numpix = 0

        ! Done looping !
     enddo
     write(6, *) 'Completed radius = ', rad(i)
  enddo
  write(6, *) 'Done'
  
  write(6, *) ''
  area = pi*(rmax**2. - rstart**2.)
  write(6, *) 'Fraction of pixels used: ', areacheck/area
  write(6, *) ''

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write the profiles out to a file for python to plot !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(43, file = outname1, status = 'unknown')
  do i = 1, raddiv
     write(43, *) az_profs(i,:)
  enddo
  open(44, file = outname2, status = 'unknown')
  do j = 1, angdiv
     write(44, *) rad_profs(j,:)
  enddo

end program az_prof
