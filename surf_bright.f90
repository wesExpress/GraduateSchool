program surf_bright
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Turns DiskFit intensities into surface brightnesses   !
  ! Needs:                                                !
  !    1. length of input arrays                          !
  !    2. zero point magnitude for filter                 !
  !    3. input name of the file                          !
  !    4. desired outputname                              !
  !    5. switch for diskfit components or single profile !
  !    6. plate scale in arcsec per pixel                 !
  !    7. exp time in seconds                             !
  ! If not using in conjunction with the sbright SM macro !
  ! need to create a file from DiskFit output which only  !
  ! has radius and intensities of the 3 components        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer expt
  real pix 
  
  real, allocatable :: Idisk(:), Ibar(:), Ibulge(:), r(:)
  real, allocatable :: Itot(:)
  real, allocatable :: surf_tot(:)
  real, allocatable :: surf_disk(:), surf_bar(:), surf_bulge(:)
  real zp
  real dum1, dum2

  integer n
  integer i, j, k
  integer l, m, p
  integer switch

  character len*40
  character zpd*40
  character inname*40
  character outname*40
  character blah*40
  character pixd*40
  character exptd*40
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read in the arguments, allocate the input arrays !
  ! and then read them in from 'inname'              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1, len)
  call getarg(2, zpd)
  call getarg(3, inname)
  call getarg(4, outname)
  call getarg(5, blah)
  call getarg(6, pixd)
  call getarg(7, exptd)

  read (len, '(i4)') n
  read (zpd, '(f6.2)') zp
  read (blah, '(i4)') switch
  read (pixd, '(f6.2)') pix
  read (exptd, '(i4)') expt

  !write(6, *) n, zp, inname, outname, switch, pix, expt

  allocate(r(n), Idisk(n), Ibar(n), Ibulge(n), Itot(n))
  allocate(surf_disk(n), surf_bar(n), surf_bulge(n))
  allocate(surf_tot(n))

  open(42, file = inname, status = 'old')
  if(switch .eq. 1) then
     do i = 1, n
        read (42, *) r(i), Idisk(i), Ibar(i), Ibulge(i)
     enddo
  elseif(switch .eq. 2) then
     do i = 1, n
        read (42, *) r(i), Itot(i)
        !write(6, *) Itot(i)
     enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Change the intensities into surface brightnesses !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(switch .eq. 1) then
     do i = 1, n
        if(Idisk(i) .ne. 0) then
           surf_disk(i) = -2.5*log10(Idisk(i)/(expt*(pix**2.))) + zp
        else
           surf_disk(i) = 0
        endif
        if(Ibar(i) .ne. 0) then
           surf_bar(i) = -2.5*log10(Ibar(i)/(expt*(pix**2.))) + zp
        else
           surf_bar(i) = 0
        endif
        if(Ibulge(i) .ne. 0) then
           surf_bulge(i) = -2.5*log10(Ibulge(i)/(expt*(pix**2.))) + zp
        else
           surf_bulge(i) = 0
        endif
     enddo
  elseif(switch .eq. 2) then
     do i = 1, n
        if(Itot(i) .ne. 0) then
           surf_tot(i) = -2.5*log10(Itot(i)/(expt*(pix**2.))) + zp
           !write(6, *) surf_tot(i)
        else
           surf_tot(i) = 0
        endif
     enddo
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Combine the individual components into total profile !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(switch .eq. 1) then
     do j = 1, n
        if(surf_disk(j) .eq. 0) then
           if(surf_bulge(j) .eq. 0) then
              surf_tot(j) = surf_bar(j)
           elseif(surf_bar(j) .eq. 0) then
              surf_tot(j) = surf_bulge(j)
           else
              surf_tot(j) = -2.5*log10(10**(-surf_bulge(j)/2.5) + 10**(-surf_bar(j)/2.5))
           endif
        elseif(surf_bar(j) .eq. 0) then
           if(surf_disk(j) .eq. 0) then
              surf_tot(j) = surf_bulge(j)
           elseif(surf_bulge(j) .eq. 0) then
              surf_tot(j) = surf_disk(j)
           else
              surf_tot(j) = -2.5*log10(10**(-surf_disk(j)/2.5) + 10**(-surf_bulge(j)/2.5))
           endif
        elseif(surf_bulge(j) .eq. 0) then
           if(surf_disk(j) .eq. 0) then
              surf_tot(j) = surf_bar(j)
           elseif(surf_bar(j) .eq. 0) then
              surf_tot(j) = surf_disk(j)
           else
              surf_tot(j) = -2.5*log10(10**(-surf_disk(j)/2.5) + 10**(-surf_bar(j)/2.5))
           endif
        else
           surf_tot(j) = -2.5*log10(10**(-surf_disk(j)/2.5) + 10**(-surf_bar(j)/2.5) + 10**(-surf_bulge(j)/2.5))
        endif
        !write(6, *) surf_tot(j)
     enddo
  endif

  open(43, file = outname, status = 'unknown')
  if(switch .eq. 1) then
     do i = 1, n
        write(43, *) r(i), surf_disk(i), surf_bar(i), surf_bulge(i), surf_tot(i)
     enddo
  elseif(switch .eq. 2) then
     do i = 1, n
        write(43, *) r(i), surf_tot(i)
     enddo
  endif
  
end program surf_bright
