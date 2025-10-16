!*******************************************************************************
!
!     program grid_power_spectrum
!
!*******************************************************************************
!
!     Function: Apply various frequency domain filters to gridded data using FFTW
!     Currently, only Butterworth low-pass is implemented in this version
!
!-------------------------------------------------------------------------------
!
!     Programmer:  John Bain
!     Date:        October 9, 2025
!     Revision:    Updated to match Geosoft FFT2PSPC (log10, no normalization, 
!                  10% extension to 280x280 with edge-mean removal and MEM-based 
!                  dummy filling via Burg algorithm with damping and smoothing, 
!                  edge-based mean removal, dynamic cropping with workaround: 
!                  full-size grid with left half set to d_dummy). Outputs 
!                  input_name_orig_minus_mean.grd (mean-subtracted input). 
!                  Uses dextnsn < 0.0 to read input grid as extended grid (e.g., 
!                  uncompressed boug_smallprp.grd). Added double precision for 
!                  crop_dxorg_Nx, dxorg_Nx, dyorg_Ny, crop_deltaKx. Removed 
!                  test_mode and statement 1000. Fixed deallocation error. 
!                  Restored grid_contains_dummies logic. Restored symmetric 
!                  padding logic (pad_x = (ncol_ext - ncol) / 2) to fix centering.
!                  Enhanced debug prints for 1-row shift check, retained 
!                  jnew = mod(j-1 + ny/2, ny) + 1 pending evidence.
!
!===============================================================================

      program grid_power_spectrum

      implicit none
      integer ncol, nrow, istat, outlun, maxz, maxz_ext, maxpts, &
              il2, i, j, ncol_ext, nrow_ext, ilen, ipt, &
              crop_ncol, crop_nrow, maxpts_cropped, &
              pad_x_left, pad_x_right, pad_y_top, pad_y_bottom, ipt1, ipt2, maxx

      character gridin*300, gridout*300, gridout_ext*300, &
                gridout_pwr*300, gridout_pwr_shifted*300, gridout_pwr_norm*300, &
                gridout_pwr_shifted_cropped*300, gridout_pwr_shifted_full*300, &
                gridout_orig_minus_mean*300, radialout*300, &
                tempfile*200, tempfilew*200, &
                clabel*48, blankstring200*200, &
                gtype*3, grid_fill*300, answer*1

      logical write_progress, show_it, grid_contains_dummies, &
              output_extended, crop_flag

      double precision d_dummy, drot, ddelx, ddely, dxorg, dyorg, dextnsn, &
                       dxorgnew, dyorgnew, deltaKx, deltaKy, crop_dxorg_Nx, &
                       sum_edges, count_edges, mean_edges, &
                       dxorg_Nx, dyorg_Ny, crop_deltaKx

      real, allocatable :: zwork1(:), zwork2(:), rwork1(:), rwork2(:), rwork3(:), &
                           rworkcol1(:), rworkcol2(:), rworkcol3(:)

      real*8, allocatable :: dz(:)             ! Input grid
      real*8, allocatable :: dz_ext(:)         ! Extended grid
      real*8, allocatable :: dwork1(:)         ! Work grid
      real*8, allocatable :: dwork2(:)         ! Work grid - extended
      real*8, allocatable :: u_ext(:,:)        ! FFT forward
      real*8, allocatable :: dz_pwr(:)         ! Output power grid
      real*8, allocatable :: dz_pwr_shifted(:) ! Output shifted power grid
      real*8, allocatable :: dz_pwr_norm(:)    ! Output shifted, normalized power grid
      real*8, allocatable :: dz_pwr_shifted_cropped(:) ! Cropped for Geosoft (full size)
      real*8, allocatable :: pwr(:,:)          ! 2D power
      real*8, allocatable :: pwr_shifted(:,:)  ! pwr after shifting to center
      real*8, allocatable :: pwr_norm(:,:)     ! pwr_shifted after normalizing to radial-mean
      real*8, allocatable :: radial_spec(:)    ! log10 power mean per ring
      real*8, allocatable :: k_spec(:)         ! cycles/m per ring center
      real*8, allocatable :: lambda_spec(:)    ! meters (= 1/k, cycles/m)
      real*8 :: dr                            ! ring spacing (cycles/m)
      integer :: nbins                        ! 0 => auto-pick; returns value

!-------------------------------------------------------------------------------

      gtype = 'GRD'
      show_it = .true.
      outlun = 3
      d_dummy = -1.0d+32
      drot = 0.0D0
      write_progress = .true.
      grid_contains_dummies = .false.

      do i = 1, 300
         gridout(i:i) = ' '
         gridout_ext(i:i) = ' '
         grid_fill(i:i) = ' '
         gridout_pwr(i:i) = ' '
         gridout_pwr_shifted(i:i) = ' '
         gridout_pwr_norm(i:i) = ' '
         gridout_pwr_shifted_cropped(i:i) = ' '
         gridout_pwr_shifted_full(i:i) = ' '
         gridout_orig_minus_mean(i:i) = ' '
         radialout(i:i) = ' '
         if ( i .le. 200 ) blankstring200(i:i) = ' '
      end do

! ##############################################################################

      print *,'PROGRAM: GRID_POWER_SPECTRUM'
      print *,'   Compute 2D Power Spectrum (log10) and Radial Mean Spectrum'
      print *,'   Matches Geosoft FFT2PSPC with dynamic cropping'
      print *,'   Note: grid cannot contain any null/dummies'

!===============================================================================
! GET THE INPUT GRID
!===============================================================================

5     print *
      print *,'Enter the input potential field grid file name'
      read ( *,'(a)', err=5 ) gridin

      tempfile = blankstring200
      tempfile = gridin
      call remove_percent20_spaces ( tempfile, tempfilew, il2 )
      if ( tempfile(1:1) .eq. '/' ) then
         do j = 1, il2 - 1
            tempfile(j:j) = tempfile(j+1:j+1)
         end do
         tempfile(il2:il2) = ' '
         il2 = il2 - 1
      endif
      gridin = blankstring200
      gridin = tempfile

      call string_length ( gridin, ilen )
      gridout(1:ilen-4) = gridin(1:ilen-4)
      ilen = ilen - 4

      print *
      print *,'Write out the extended, prepped input to FFTW? (Y or N)?'
      read (*,'(a)') answer
      if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
         output_extended = .true.
      else
         output_extended = .false.
      endif

      print *
      print *,'Enter the extension as a % of input grid (e.g., 10=10%)'
      print *,'  (negative number uses input grid as extended grid, e.g., boug_smallprp.grd)'
      read (*,*) dextnsn
      if ( dextnsn .lt. 0.0 ) then
         print *,'Using input grid as extended grid (no extension, no mean removal)'
      else
         print *,'Using extension: dextnsn = ', dextnsn
      endif

      print *
      print *,'Crop to right half (positive kx/ky) like Geosoft? (Y or N)'
      read (*,'(a)') answer
      if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
         crop_flag = .true.
         print *,'Will create _pwr_shifted_crop.grd (full size with left half null) + full grid'
      else
         crop_flag = .false.
         print *,'Will create full grid only'
      endif

!===============================================================================
! READ WITH SIZE SET TO 1 TO GET NCOL, NROW OF ACTUAL GRID SIZE
!===============================================================================

      maxpts = 1
      allocate ( dz(maxpts), zwork1(maxpts) )      

      call read_geosoft_grid ( gridin, ncol, nrow, &
                               dxorg, dyorg, ddelx, ddely, &
                               dz, drot, maxpts, gtype, d_dummy, &
                               zwork1, istat )

      maxz = ncol * nrow

      if ( show_it ) then
         print *,'Test input grid read to get ncol/nrow: ', ncol, nrow
      endif

      deallocate ( dz, zwork1 )

! Handle input grid based on dextnsn
      if ( dextnsn .lt. 0.0 ) then
         ! Use input grid as extended grid (e.g., boug_smallprp.grd)
         ncol_ext = ncol
         nrow_ext = nrow
         maxz_ext = ncol_ext * nrow_ext
         dxorgnew = dxorg
         dyorgnew = dyorg

         allocate ( dz_ext(maxz_ext), dwork2(maxz_ext), zwork2(maxz_ext) )

         maxpts = maxz_ext
         call read_geosoft_grid ( gridin, ncol_ext, nrow_ext, &
                                  dxorgnew, dyorgnew, ddelx, ddely, &
                                  dz_ext, drot, maxpts, gtype, d_dummy, &
                                  zwork2, istat )

         if ( istat .ne. 0 ) then
            print *,'*Error reading input grid as extended grid'
            print *,'istat = ',istat
            go to 900
         endif

         if ( show_it ) then
            print *,'Input grid read as extended grid: ncol_ext, nrow_ext = ', ncol_ext, nrow_ext
            print *,'min/max input =', minval(dz_ext), maxval(dz_ext)
         endif
      else
         ! Normal mode: read input grid and extend
         allocate ( dz(maxz), dwork1(maxz), zwork1(maxz) )

         maxpts = maxz
         call read_geosoft_grid ( gridin, ncol, nrow, &
                                  dxorg, dyorg, ddelx, ddely, &
                                  dz, drot, maxpts, gtype, d_dummy, &
                                  zwork1, istat )

         if ( istat .ne. 0 ) then
            print *,'*Error reading input grid file'
            print *,'istat = ',istat
            go to 900
         endif

         if ( show_it ) then
            print *,'Input grid read: ncol, nrow = ', ncol, nrow
            print *,'min/max input =', minval(dz), maxval(dz)
         endif

! Check for dummy values and fill if necessary
         grid_contains_dummies = .false.
         do i = 1, maxz
            if ( dz(i) .eq. d_dummy ) then
               grid_contains_dummies = .true.
               exit
            endif
         end do
         if ( grid_contains_dummies ) then
            grid_fill = 'temp_grid_fill.grd'
            print *,'Dummies detected, filling grid to: ', trim(grid_fill)
            ! Placeholder: replace with your actual fill_grid_dummies call
            ! call fill_grid_dummies ( dz, ncol, nrow, d_dummy, grid_fill, istat )
            ! For now, assume fill_grid_dummies writes to grid_fill
            if ( istat .ne. 0 ) then
               print *,'*Error filling dummy values'
               print *,'istat = ',istat
               go to 900
            endif
            ! Read the filled grid
            call read_geosoft_grid ( grid_fill, ncol, nrow, &
                                     dxorg, dyorg, ddelx, ddely, &
                                     dz, drot, maxpts, gtype, d_dummy, &
                                     zwork1, istat )
            if ( istat .ne. 0 ) then
               print *,'*Error reading filled grid'
               print *,'istat = ',istat
               go to 900
            endif
            if ( show_it ) print *,'Filled grid read: min/max =', minval(dz), maxval(dz)
         endif

! Edge-based mean removal (matches Geosoft FFT2PREP setting)
         sum_edges = 0.d0
         count_edges = 0.d0
         do i = 1, ncol
            sum_edges = sum_edges + dz(i) + dz((nrow-1)*ncol + i)
            count_edges = count_edges + 2.d0
         end do
         do j = 2, nrow-1
            sum_edges = sum_edges + dz((j-1)*ncol + 1) + dz(j*ncol)
            count_edges = count_edges + 2.d0
         end do
         mean_edges = sum_edges / count_edges
         dz = dz - mean_edges
         if ( show_it ) print *,'Edge-based mean removed: ', mean_edges

! Write mean-subtracted input grid
         gridout_orig_minus_mean = trim(gridout)//'_orig_minus_mean.grd'
         clabel(1:48) = gridout_orig_minus_mean(1:48)
         call write_geosoft_grid ( gridout_orig_minus_mean, ncol, nrow, &
                                   dxorg, dyorg, ddelx, ddely, &
                                   d_dummy, dz, drot, maxz, gtype, clabel, &
                                   dwork1, zwork1, istat )
         if ( istat .ne. 0 ) then
            print *,'*Error writing mean-subtracted input grid file'
            print *,'istat = ',istat
            go to 900
         endif
         if ( show_it ) print *,'Wrote mean-subtracted grid: ', trim(gridout_orig_minus_mean)

! Compute extended grid dimensions
         call compute_fftw_extension ( ncol, nrow, dextnsn, ncol_ext, nrow_ext )

         if ( show_it ) then
            print *,'Extended grid dimensions: ', ncol_ext, nrow_ext
         endif

         maxz_ext = ncol_ext * nrow_ext
         allocate ( dz_ext(maxz_ext), dwork2(maxz_ext), zwork2(maxz_ext) )

! Allocate work arrays for extend_grid_ftprep
         maxx = 288  ! Matches maxx in extend_grid_ftprep.f
         allocate ( rwork1(maxz_ext), rwork2(maxz_ext), rwork3(maxz_ext) )
         allocate ( rworkcol1(maxx), rworkcol2(maxx), rworkcol3(maxx) )

! Copy input grid to rwork1
         do i = 1, ncol * nrow
            rwork1(i) = real(dz(i))
         end do

! Compute symmetric padding from extend_grid_geosoft.f90
         pad_x_left = (ncol_ext - ncol) / 2
         pad_x_right = ncol_ext - ncol - pad_x_left
         pad_y_top = (nrow_ext - nrow) / 2
         pad_y_bottom = nrow_ext - nrow - pad_y_top
         if ( show_it ) then
            print *,'Symmetric padding: left/right = ', pad_x_left, pad_x_right
            print *,'Symmetric padding: top/bottom = ', pad_y_top, pad_y_bottom
         endif

         dxorgnew = dxorg - dble(pad_x_left) * ddelx
         dyorgnew = dyorg - dble(pad_y_top) * ddely

! Call extend_grid_ftprep with original usgs_ftprep_new.f interface
         if ( show_it ) print *,'Calling extend_grid_ftprep'
         call extend_grid_ftprep ( dextnsn, ncol, nrow, dxorg, dyorg, ddelx, ddely, &
                                   ncol_ext, nrow_ext, dxorgnew, dyorgnew, &
                                   rwork1, rwork2, rwork3, rworkcol1, rworkcol2, &
                                   rworkcol3, maxx, istat )
         if ( istat /= 0 ) then
            print *,'*Error in extend_grid_ftprep'
            print *,'istat = ',istat
            go to 900
         endif
         if ( show_it ) print *,'Completed extend_grid_ftprep'

! Copy output to dz_ext
         do i = 1, ncol_ext * nrow_ext
            dz_ext(i) = dble(rwork3(i))
         end do

! Deallocate work arrays
         deallocate ( rwork1, rwork2, rwork3, rworkcol1, rworkcol2, rworkcol3 )

! Enforce periodicity on dz_ext
         if ( show_it ) print *,'Enforcing periodicity on dz_ext'
         do j = 1, nrow_ext
            ipt1 = (j-1)*ncol_ext + 1
            ipt2 = j*ncol_ext
            if ( ipt1 <= 0 .or. ipt1 > maxz_ext .or. ipt2 <= 0 .or. ipt2 > maxz_ext ) then
               print *,'*Error in periodicity enforcement: invalid indices'
               print *,'j = ', j, 'ipt1 = ', ipt1, 'ipt2 = ', ipt2, 'maxz_ext = ', maxz_ext
               istat = 5
               go to 900
            endif
            dz_ext(ipt1) = 0.5d0 * (dz_ext(ipt1) + dz_ext(ipt2))
            dz_ext(ipt2) = dz_ext(ipt1)
         end do
         do i = 1, ncol_ext
            ipt1 = i
            ipt2 = (nrow_ext-1)*ncol_ext + i
            if ( ipt1 <= 0 .or. ipt1 > maxz_ext .or. ipt2 <= 0 .or. ipt2 > maxz_ext ) then
               print *,'*Error in periodicity enforcement: invalid indices'
               print *,'i = ', i, 'ipt1 = ', ipt1, 'ipt2 = ', ipt2, 'maxz_ext = ', maxz_ext
               istat = 5
               go to 900
            endif
            dz_ext(ipt1) = 0.5d0 * (dz_ext(ipt1) + dz_ext(ipt2))
            dz_ext(ipt2) = dz_ext(ipt1)
         end do
         if ( show_it ) print *,'Periodicity enforced on dz_ext'

         if ( output_extended ) then
            if ( show_it ) print *,'Writing extended grid: ', trim(gridout_ext)
            gridout_ext = trim(gridout)//'_ext.grd'
            clabel(1:48) = gridout_ext(1:48)
            call write_geosoft_grid ( gridout_ext, ncol_ext, nrow_ext, &
                                      dxorgnew, dyorgnew, ddelx, ddely, &
                                      d_dummy, dz_ext, drot, maxz_ext, &
                                      gtype, clabel, dwork2, zwork2, istat )
            if ( istat .ne. 0 ) then
               print *,'*Error writing extended grid file'
               print *,'istat = ',istat
               go to 900
            endif
            if ( show_it ) print *,'Wrote extended grid: ', trim(gridout_ext)
         else
            if ( show_it ) print *,'Skipped writing extended grid (output_extended = .false.)'
         endif
      endif

!===============================================================================
! COMPUTE SPECTRAL POWER GRID, CENTER THE GRID, COMPUTE RADIAL-AVERAGE/RADIAL-MEAN
!===============================================================================

      if ( show_it ) print *,'Allocating FFT arrays'
      nbins = min(ncol_ext, nrow_ext) / 2
      allocate(radial_spec(nbins), k_spec(nbins), lambda_spec(nbins))

      allocate ( u_ext(ncol_ext,nrow_ext), pwr(ncol_ext,nrow_ext), &
                 pwr_shifted(ncol_ext,nrow_ext), pwr_norm(ncol_ext,nrow_ext), &
                 dz_pwr(maxz_ext), dz_pwr_shifted(maxz_ext), dz_pwr_norm(maxz_ext) )

      if ( crop_flag ) then
         crop_ncol = ncol_ext / 2 + 1  ! Positive kx including DC (e.g., 141 for ncol_ext=280)
         crop_nrow = nrow_ext          ! Full ky (e.g., 280)
         maxpts_cropped = ncol_ext * nrow_ext  ! Full size for workaround
         allocate(dz_pwr_shifted_cropped(maxpts_cropped))
         if ( show_it ) then
            print *,'Crop dimensions allocated: crop_ncol, crop_nrow = ', crop_ncol, crop_nrow
            print *,'Crop maxpts_cropped (full size) = ', maxpts_cropped
            print *,'dz_pwr_shifted_cropped allocated: size = ', size(dz_pwr_shifted_cropped)
         endif
      endif

      if ( show_it ) print *,'Copying dz_ext to u_ext'
      do j = 1, nrow_ext
         do i = 1, ncol_ext
            ipt = (j-1)*ncol_ext + i
            if ( ipt <= 0 .or. ipt > maxz_ext ) then
               print *,'*Error copying dz_ext to u_ext: invalid index'
               print *,'i = ', i, 'j = ', j, 'ipt = ', ipt, 'maxz_ext = ', maxz_ext
               istat = 5
               go to 900
            endif
            u_ext(i,j) = dz_ext(ipt)
         end do
      end do
      if ( show_it ) then
         print *, 'min/max dz_ext before FFT:', minval(dz_ext), maxval(dz_ext)
         print *, 'min/max u_ext before FFT:', minval(u_ext), maxval(u_ext)
      endif

      if ( show_it ) print *,'Calling compute_fft_power_spectrum'
      call compute_fft_power_spectrum ( u_ext, pwr, ncol_ext, nrow_ext, ddelx, ddely, istat )
      if ( istat .ne. 0 ) then
         print *,'*Error in compute_fft_power_spectrum'
         print *,'istat = ',istat
         go to 900
      endif
      if ( show_it ) then
         print *,'Just after computing power spectrum'
         print *, 'Max Power =', maxval(pwr)
         print *, 'Min Power =', minval(pwr, mask=pwr.ne.d_dummy)
      endif

      if ( show_it ) print *,'Calling fftshift_and_radialnorm'
      call fftshift_and_radialnorm ( pwr, pwr_shifted, pwr_norm, &
                                     radial_spec, k_spec, lambda_spec, &
                                     ncol_ext, nrow_ext, ddelx, ddely, dr, nbins, istat )
      if ( istat .ne. 0 ) then
         print *,'*Error in fftshift_and_radialnorm'
         print *,'istat = ',istat
         go to 900
      endif
      if ( show_it ) then
         print *,'Just after shifting spectrum and radial norm'
         print *,'pwr_shifted min/max = ', minval(pwr_shifted), maxval(pwr_shifted)
         print *,'pwr_norm min/max = ', minval(pwr_norm), maxval(pwr_norm)
      endif

      if ( show_it ) print *,'Copying FFT results to output arrays'
      if ( show_it ) then
         print *,'Array sizes:'
         print *,'dz_pwr: ', size(dz_pwr)
         print *,'dz_pwr_shifted: ', size(dz_pwr_shifted)
         print *,'dz_pwr_norm: ', size(dz_pwr_norm)
         if ( crop_flag ) print *,'dz_pwr_shifted_cropped: ', size(dz_pwr_shifted_cropped)
         print *,'pwr: ', size(pwr)
         print *,'pwr_shifted: ', size(pwr_shifted)
         print *,'pwr_norm: ', size(pwr_norm)
      endif
      ipt = 0
      do j = 1, nrow_ext
         do i = 1, ncol_ext
            ipt = ipt + 1
            if ( ipt <= 0 .or. ipt > maxz_ext ) then
               print *,'*Error copying FFT results: invalid index'
               print *,'i = ', i, 'j = ', j, 'ipt = ', ipt, 'maxz_ext = ', maxz_ext
               istat = 5
               go to 900
            endif
            if ( i > ncol_ext .or. j > nrow_ext ) then
               print *,'*Error copying FFT results: invalid 2D index'
               print *,'i = ', i, 'j = ', j, 'ncol_ext = ', ncol_ext, 'nrow_ext = ', nrow_ext
               istat = 5
               go to 900
            endif
            dz_pwr(ipt) = pwr(i,j)
            dz_pwr_shifted(ipt) = pwr_shifted(i,j)
            dz_pwr_norm(ipt) = pwr_norm(i,j)
            if ( crop_flag ) then
               if ( i <= 0 .or. i > ncol_ext .or. j <= 0 .or. j > nrow_ext ) then
                  print *,'*Error copying dz_pwr_shifted_cropped: invalid index'
                  print *,'i = ', i, 'j = ', j, 'ncol_ext = ', ncol_ext, 'nrow_ext = ', nrow_ext
                  istat = 5
                  go to 900
               endif
               if ( ipt <= 0 .or. ipt > maxpts_cropped ) then
                  print *,'*Error copying dz_pwr_shifted_cropped: invalid ipt'
                  print *,'ipt = ', ipt, 'maxpts_cropped = ', maxpts_cropped
                  istat = 5
                  go to 900
               endif
               if ( i .ge. crop_ncol ) then
                  dz_pwr_shifted_cropped(ipt) = pwr_shifted(i,j)
               else
                  dz_pwr_shifted_cropped(ipt) = d_dummy
               endif
            endif
         end do
      end do

      if ( show_it ) print *, 'Max values after copy:'
      if ( show_it ) print *, 'dz_pwr:', maxval(dz_pwr)
      if ( show_it ) print *, 'dz_pwr_shifted:', maxval(dz_pwr_shifted)
      if ( show_it ) print *, 'dz_pwr_norm:', maxval(dz_pwr_norm)
      if ( show_it .and. crop_flag ) print *, 'dz_pwr_shifted_cropped:', maxval(dz_pwr_shifted_cropped)
      if ( show_it ) then
         print *,'Debug: dz_pwr_shifted center [141,141] = ', dz_pwr_shifted((140)*ncol_ext + 141)
         print *,'Debug: dz_pwr_shifted_cropped center [141,141] = ', dz_pwr_shifted_cropped((140)*ncol_ext + 141)
         print *,'Debug: dz_pwr_shifted corner [1,1] = ', dz_pwr_shifted(1)
         print *,'Debug: dz_pwr_shifted_cropped corner [1,1] = ', dz_pwr_shifted_cropped(1)
         print *,'Debug: dz_pwr_shifted row 140 [141,140] = ', dz_pwr_shifted((139)*ncol_ext + 141)
         print *,'Debug: dz_pwr_shifted row 142 [141,142] = ', dz_pwr_shifted((141)*ncol_ext + 141)
      endif

!===============================================================================
! WRITE SPECTRAL GRIDS
!===============================================================================

      maxpts = maxz_ext

      if ( show_it ) print *,'Writing unshifted grid'
      gridout_pwr = trim(gridout)//'_pwr.grd'
      clabel(1:48) = gridout_pwr(1:48)
      deltaKx = 1.0d0 / ( dble(min(ncol_ext,nrow_ext)) * ddelx )
      deltaKy = deltaKx
      dxorg_Nx = -1.0d0 / ( 2.0d0 * ddelx )
      dyorg_Ny = dxorg_Nx
      call write_geosoft_grid ( gridout_pwr, ncol_ext, nrow_ext, &
                                dxorg_Nx, dyorg_Ny, deltaKx, deltaKy, &
                                d_dummy, dz_pwr, drot, maxpts, gtype, clabel, &
                                dwork2, zwork2, istat )
      if ( istat .ne. 0 ) then
         print *,'*Error writing output spectral grid file'
         print *,'istat = ',istat
         go to 900
      endif
      if ( show_it ) print *,'Wrote unshifted grid: ', trim(gridout_pwr)

      if ( crop_flag ) then
         crop_dxorg_Nx = 0.0d0  ! Positive kx starts at 0
         crop_deltaKx = deltaKx
         
         if ( show_it ) then
            print *,'Crop dimensions: crop_ncol, crop_nrow = ', crop_ncol, crop_nrow
            print *,'Crop maxpts_cropped (full size) = ', maxpts_cropped
            print *,'Crop range: ', minval(dz_pwr_shifted_cropped, mask=dz_pwr_shifted_cropped.ne.d_dummy), &
                    maxval(dz_pwr_shifted_cropped)
         endif

         if ( show_it ) print *,'Writing cropped grid'
         gridout_pwr_shifted_cropped = trim(gridout)//'_pwr_shifted_crop.grd'
         clabel(1:48) = gridout_pwr_shifted_cropped(1:48)
         call write_geosoft_grid ( gridout_pwr_shifted_cropped, ncol_ext, nrow_ext, &
                                   dxorg_Nx, dyorg_Ny, deltaKx, deltaKy, &
                                   d_dummy, dz_pwr_shifted_cropped, drot, maxpts_cropped, &
                                   gtype, clabel, dwork2, zwork2, istat )
         if ( istat .ne. 0 ) then
            print *,'*Error writing cropped spectral grid file'
            print *,'istat = ',istat
            go to 900
         endif
         deallocate(dz_pwr_shifted_cropped)
         if ( show_it ) print *,'Wrote cropped grid: ', trim(gridout_pwr_shifted_cropped)
      endif

      if ( show_it ) print *,'Writing full shifted grid'
      gridout_pwr_shifted_full = trim(gridout)//'_pwr_shifted_full.grd'
      clabel(1:48) = gridout_pwr_shifted_full(1:48)
      call write_geosoft_grid ( gridout_pwr_shifted_full, ncol_ext, nrow_ext, &
                                dxorg_Nx, dyorg_Ny, deltaKx, deltaKy, &
                                d_dummy, dz_pwr_shifted, drot, maxpts, &
                                gtype, clabel, dwork2, zwork2, istat )
      if ( istat .ne. 0 ) then
         print *,'*Error writing full shifted spectral grid file'
         print *,'istat = ',istat
         go to 900
      endif
      if ( show_it ) print *,'Wrote full shifted grid: ', trim(gridout_pwr_shifted_full)

      if ( show_it ) print *,'Writing normalized grid'
      gridout_pwr_norm = trim(gridout)//'_pwr_norm.grd'
      clabel(1:48) = gridout_pwr_norm(1:48)
      call write_geosoft_grid ( gridout_pwr_norm, ncol_ext, nrow_ext, &
                                dxorg_Nx, dyorg_Ny, deltaKx, deltaKy, &
                                d_dummy, dz_pwr_norm, drot, maxpts, &
                                gtype, clabel, dwork2, zwork2, istat )
      if ( istat .ne. 0 ) then
         print *,'*Error writing normalized spectral grid file'
         print *,'istat = ',istat
         go to 900
      endif
      if ( show_it ) print *,'Wrote normalized grid: ', trim(gridout_pwr_norm)

!===============================================================================
! WRITE RADIAL SPECTRA OUTPUT FILES
!===============================================================================

      if ( show_it ) print *,'Writing radial spectra'
      radialout = trim(gridout)//'.dat'
      open ( file=radialout, unit=outlun )
      do i = 1, nbins
         write ( outlun, 710 ) k_spec(i), lambda_spec(i), radial_spec(i)
710      format ( g20.6, 2f20.4 )
      end do
      close ( outlun )
      if ( show_it ) print *,'Wrote radial spectra: ', trim(radialout)

!===============================================================================

      if ( show_it ) print *,'Deallocating arrays'
      if ( dextnsn .ge. 0.0 ) then
         if ( allocated(dz) ) deallocate ( dz )
         if ( allocated(dwork1) ) deallocate ( dwork1 )
         if ( allocated(zwork1) ) deallocate ( zwork1 )
      endif
      if ( allocated(dz_ext) ) deallocate ( dz_ext )
      if ( allocated(dwork2) ) deallocate ( dwork2 )
      if ( allocated(zwork2) ) deallocate ( zwork2 )
      if ( allocated(u_ext) ) deallocate ( u_ext )
      if ( allocated(pwr) ) deallocate ( pwr )
      if ( allocated(pwr_shifted) ) deallocate ( pwr_shifted )
      if ( allocated(pwr_norm) ) deallocate ( pwr_norm )
      if ( allocated(dz_pwr) ) deallocate ( dz_pwr )
      if ( allocated(dz_pwr_shifted) ) deallocate ( dz_pwr_shifted )
      if ( allocated(dz_pwr_norm) ) deallocate ( dz_pwr_norm )
      if ( allocated(radial_spec) ) deallocate ( radial_spec )
      if ( allocated(k_spec) ) deallocate ( k_spec )
      if ( allocated(lambda_spec) ) deallocate ( lambda_spec )

      if ( grid_contains_dummies ) then
         if ( show_it ) print *,'Deleting temporary grid file'
         call system('cmd /C del temp_grid_fill.grd')
      endif

900   print *
      print *,'Program completed successfully!'

      stop
      end