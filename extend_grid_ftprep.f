************************************************************************
*
      subroutine extend_grid_ftprep ( pct, ncolin, nrowin, dxorg, dyorg,
     + ddx, ddy, ncnew, nrnew, dxorgnew, dyorgnew, rwork1, rwork2,
     + rwork3, rworkcol1, rworkcol2, rworkcol3, maxx, ierr )
*
************************************************************************
*
*     Function: From Lin Cordell / Jeff Phillips routines to prep a grid
*     Note: Bain comments are listed as "c#".  Original comments are "c"
*
*-----------------------------------------------------------------------
*
*     Programmer:  John Bain
*     Date:        May 16, 2019
*     Notes:       Passed array "rwork1" contains the grid in single precision
*                  rwork3 contains the (output) prepped grid
*
* ----------------------------------------------------------------------
*
*     Revised by:  John Bain 
*     Date:        August 7, 2019
*     Revision:    Number of arguments changed for subroutine getrow
*                  (note: I had already done this, evidently) 
*
* ----------------------------------------------------------------------
*
*     Revised by:  John Bain 
*     Date:        August 23, 2022
*     Revision:    Massive amount of testing and cleanup in the code.
*                  I discovered that the test (Trois) data grid was very
*                  goofy, with some noise, that was causing havoc in the
*                  "extend" subroutine.  I added the smoothing routine
*                  to the original grid before the first extend, and this
*                  seems to have fixed this issue, without causing any
*                  major changes to the data.
*
* ----------------------------------------------------------------------
*
*=======================================================================
*
*
c  ftprep.for
c
c  Extends a grid to reduce FFT wraparound.
c
c  prep:   Lin Cordell Jan 1990.
c          Revised June 1991.
c  prep3:  Jeff Phillips February 1996  prediction filter extension
c  prep4:  Jeff Phillips August 1999    added centering and min,max limits
c  ftprep: Jeff Phillips September 2001 modified for use with FTFWD
c          Jeff Phillips December 2002 modified for use as a Geosoft GX
c
*
*===============================================================================
*
* Bain (16May19) declarations
*
      integer ncolin, nrowin, ncolnew, nrownew, maxpts, itemp
     + ,ncnew, nrnew, maxx, ierr, i, istat, maxz
     + ,nc, nc1, ncol, ndim, nr, nr1, nrow, nx, ny
*
      real rwork1(*), rwork2(*), rwork3(*), pct
     +    ,rworkcol1(maxx), rworkcol2(*), rworkcol3(*), zmin, zmax
     +    ,rfactor
*
      double precision dxorg, dyorg, ddx, ddy, dxorgnew, dyorgnew
     +  ,dznull, drot
*
      character gridout*48, clabel*48, blank48*48, grid_type*3
*
*===============================================================================
*
        common /extreme/ gmin,gmax
*
*===============================================================================
*
* Local double precision and real arrays for temporary output for testing
*
      double precision, allocatable ::  dz(:), dwork(:)
      real, allocatable :: rwork(:)
*
*===============================================================================
*
      maxz = maxx * maxx
      allocate ( dz(maxz), dwork(maxz), rwork(maxz) )
*
*===============================================================================
*
* Bain (16May19): For above, grid is already open and data passed to ftprep
* Set some grid parameters for local variables
*
*===============================================================================
*
      nc = ncolin
      nr = nrowin
      xo = sngl ( dxorg )
      yo = sngl ( dyorg )
      dx = sngl ( ddx )
      dy = sngl ( ddy )
*
c      print *,'Sub: ftprep #1: nc,nr,xo,yo,dx,dy',nc,nr,xo,yo,dx,dy
*
* maxpts is maximum number of z values, so maxcol * maxrow
* But for now, I will set this to ndim, which should be maximum number
* of either rows or columns
*
      ndim = maxx
*
      do i = 1, 48
         blank48(i:i) = ' '
      enddo
*
*
*===============================================================================
*
*
      if(nc.gt.ndim.or.nr.gt.ndim) then
         print *,'Bain error - nc.gt.ndim.or.nr.gt.ndim'
         ierr=1
         goto 9999
      endif
      if(nc.lt.15.or.nr.lt.15) then
         print *,'Bain error - nc.lt.15.or.nr.lt.15'
         print *,'NC = ',nc, ' NR = ',nr
         ierr=2
         go to 9999
      endif
*
* Pass input grid array
*
* Bain (14Oct24): commenting out this check as conflict with grid utilities
*
c      call check(rwork1,nc,nr,rworkcol1,ierr)
c      if(ierr.eq.1) go to 9999
*
*
*###############################################################################
*
* Bain (14Oct25): use explicit ncol/nrow output sizes and jump below other logic
*
      ncol = ncnew
      nrow = nrnew
      go to 800
*
*###############################################################################
*
*
*=======================================================================
*
* Set new grid geometry
*
* Bain (27Feb21) - I put in a rounding factor after having a resolution
* issue with pct being single precision (? - seems strange to need this...)
*
      rfactor = 1.0+((pct/100.0)+(pct/100000.0))
      nc1=int(float(nc)*rfactor)
      nr1=int(float(nr)*rfactor)
*
* ----------------------------------------------------------------------
*
* Bain (27Feb21) - temporary - comment out and set ncol to nc1, nrow to nr1,
* so no extension beyond the pct specified
* The two lines surrounded by $$$ were only used for testing Curie_depth
* - Otherwise, leave prime5 as active
*
* Bain (21Apr21): for this local version I am setting this to be nc1/nr1
* if pct is zero
*
*
*=======================================================================
*
      if ( pct .lt. 0.00000001 ) then
         ncol = nc1
         nrow = nr1
      else
         call prime5(nc1,ncol)
         call prime5(nr1,nrow)
      endif
*
*=======================================================================
*
*
c        print *,'nc, nc1, nr, nr1, ncol, nrow',nc,nc1,nr,nr1,ncol,nrow
*
c        zmin = 1e20
c        zmax = -1e20
c        do i = 1, nc*nr
c           zmin = min ( zmin, rwork1(i) )
c           zmax = max ( zmax, rwork1(i) )
c        enddo
c        print *,'Just after prime5: zmin, zmax: ',zmin, zmax
*
*=======================================================================
*
      if(nc.gt.ncol.or.nr.gt.nrow) then
         ierr=3
         goto 9999          
      endif
      if(ncol.gt.ndim.or.nrow.gt.ndim) then
         print*,'Maximum rows or columns =',ndim
         ierr=4
         go to 9999          
      endif
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
      go to 800
*
* Temporary write of input grid for testing
*
      gridout = blank48
      clabel = blank48
      drot = 0.0
      grid_type = 'GRD'
      gridout(1:12) = 'ftprep_1.grd'
      clabel = gridout
      ncnew = nc
      nrnew = nr
      do i = 1, ncnew*nrnew
         dz(i) = dble ( rwork1(i) )
      enddo
*
      call write_geosoft_grid ( gridout, ncnew, nrnew, dxorg,
     +                          dyorg, ddx, ddy, dznull, dz,
     +                          drot, maxz, grid_type, clabel, dwork,
     +                          rwork, istat )
c      stop
*
800   continue
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*
*===============================================================================
*
* Bain (23Aug22): call smooth on original input grid, to eliminate issues
* caused by noise in the input grid
* Commented out - didn't really help with Trois test data set
*
*===============================================================================
*
c      call smooth ( rwork1, rwork2, nr, ndim, rworkcol1, rworkcol2,
c     +              nc, nr, ierr )
c      if(ierr.eq.1) go to 9999
* Bain (23Aug22): to keep things below consistent, copy rwork2 back in to rwork1
*
c      do i = 1, nc * nr
c         rwork1(i) = rwork2(i)
c      enddo
*
*
*===============================================================================
*
* Bain (16May19): call extend with new passed arrays, store output in arrays
*
*===============================================================================
*
c      print *,'Just before call extend #1'
*
      call extend (nc, nr, nr, ncol, rwork1, rwork2, rworkcol1,
     +             rworkcol2, rworkcol3, ierr)
*
c      print *,'Just after call extend #1, nc, nr, ncol',nc, nr, ncol
*
* We now have a grid that is ncol columns wide, and nr rows tall, with
* all extension applied to the right side of the grid (by ncol - nc nodes)
*
* This grid resides in rwork2, and the original grid resides in rwork1
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
      go to 810
*
* Temporary write of expanded grid for testing
*
      gridout = blank48
      clabel = blank48
      drot = 0.0
      grid_type = 'GRD'
      gridout(1:12) = 'ftprep_2.grd'
      clabel = gridout
      ncnew = ncol
      nrnew = nr
      do i = 1, ncnew*nrnew
         dz(i) = dble ( rwork2(i) )
      enddo
*
      call write_geosoft_grid ( gridout, ncnew, nrnew, dxorg,
     +                          dyorg, ddx, ddy, dznull, dz,
     +                          drot, maxz, grid_type, clabel, dwork,
     +                          rwork, istat )
c      stop
*
810   continue
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*
*===============================================================================
*
* Bain (16May19): Transpose grid with new calling arguments
* This now uses a transposed array with size: (nr, ncol) =
* the extended columns by original rows, including swapping the XY origins
*
*===============================================================================
*
* Read from extended grid (in rwork2), and fill column-oriented grid -
* placing in work array rwork3
*
c      print *,'Just before call transpos #1'
*
      call bain_transpose ( rwork2, rwork1, rworkcol1, ndim,
     +                      ncol, nr, ierr)
      if(ierr.eq.1) go to 9999
*
c      print *,'Just after call transpos #1, ndim, ncol',ndim, ncol, nr
*
* This has each column in the original grid turned on its side, and put
* out to rwork1 array.  So dimensions of this are now (nr, ncol),
* i.e., transposed version of the column-extended grid
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
      go to 820
*
* Temporary write of expanded grid for testing
*
      gridout = blank48
      clabel = blank48
      drot = 0.0
      grid_type = 'GRD'
      gridout(1:12) = 'ftprep_3.grd'
      clabel = gridout
      ncnew = nr
      nrnew = ncol
      do i = 1, ncnew*nrnew
         dz(i) = dble ( rwork1(i) )
      enddo
*
      call write_geosoft_grid ( gridout, ncnew, nrnew, dxorg,
     +                          dyorg, ddx, ddy, dznull, dz,
     +                          drot, maxz, grid_type, clabel, dwork,
     +                          rwork, istat )
c      stop
*
820   continue
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*
*===============================================================================
*
* Bain (16May19): extend grid on right side, which is now the top of the
* rows of the original grid, including smoothing before (as of 23Aug22)
*
* Didn't really help for Trois test data set
*
*===============================================================================
*
c      call smooth ( rwork1, rwork2, ncol, ndim, rworkcol1, rworkcol2,
c     +              nr, ncol, ierr )
c      if(ierr.eq.1) go to 9999
*
* Bain (23Aug22): to keep things below consistent, copy rwork2 back in to rwork1
*
c      do i = 1, ncol * nr
c         rwork1(i) = rwork2(i)
c      enddo
*
c      print *,'Just before call extend #2'
*
      call extend (nr, ncol, nc, nrow, rwork1, rwork2, rworkcol1,
     +             rworkcol2, rworkcol3, ierr)
      if(ierr.eq.1) go to 9999
*
c      print *,'Just after call extend #2'
*
* grid in rwork2 is now transposed, with width = nrow, and height = ncol
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
      go to 830
*
* Temporary write of expanded grid for testing
*
      gridout = blank48
      clabel = blank48
      drot = 0.0
      grid_type = 'GRD'
      gridout(1:12) = 'ftprep_4.grd'
      clabel = gridout
      ncnew = nrow
      nrnew = ncol
      do i = 1, ncnew*nrnew
         dz(i) = dble ( rwork2(i) )
      enddo
*
      call write_geosoft_grid ( gridout, ncnew, nrnew, dxorg,
     +                          dyorg, ddx, ddy, dznull, dz,
     +                          drot, maxz, grid_type, clabel, dwork,
     +                          rwork, istat )
c      stop
*
830   continue
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*
*===============================================================================
*
* Bain (16May19): smooth grid
*
*===============================================================================
*
*
c      print *,'Just before call smooth #1'
*
* Note: grid is transposed, so reversing the ncol and nrow values
*
      call smooth ( rwork2, rwork1, nc, ndim, rworkcol1, rworkcol2,
     +              nrow, ncol, ierr )
c      print *,'Just after call smooth #1'
*
      if(ierr.eq.1) go to 9999
*
*
*===============================================================================
*
* Bain (16May19): Transpose back
*
*===============================================================================
*
c      print *,'Just before call transpos #2'
*
* This grid is transposed - so ncol,nrow becomes nrow,ncol for this call
*
996   continue
*
      call bain_transpose ( rwork1, rwork2, rworkcol1, ndim,
     +                      nrow, ncol, ierr)
      if(ierr.eq.1) go to 9999
*
c      print *,'Just after call transpos #2'
*
* Bain (15Aug22): jump to 998 to test
*
c      go to 998
*
*
* rwork2 is now transposed back to have width = ncol and height = nrow
*
*
*===============================================================================
*
* Bain (16May19): smooth grid again
*
*===============================================================================
*
* Note: grid is no longer transposed, so going back to having ncol, nrow (proper order)
*
c      print *,'Just before call smooth #2'
      call smooth ( rwork2, rwork1, nr, ndim, rworkcol1, rworkcol2,
     +              ncol, nrow, ierr )
      if(ierr.eq.1) go to 9999
c      print *,'Just after call smooth #2'
*
*
* rwork1 now has smoothing in both directions
*
*
*===============================================================================
*
* Bain (16May19): center data
*
*===============================================================================
*
*
      nx=(ncol-nc)/2
      ny=(nrow-nr)/2
*
c      print *,'Sub: ftprep #2: nx,ncol,nc,xo,yo',nx, ncol, nc, xo, yo
*
      xo=xo-nx*dx
      yo=yo-ny*dy
*
c      print *,'Sub: ftprep #3: nx, ncol, nc, xo, yo, dx, dy'
c      print *,nx, ncol, nc, xo, yo, dx, dy
*
      go to 994
*
*
*###############################################################################
*
*
* Bain (15Aug22): test copying rwork1 to 3, instead of next call (which is now commented
*
998   do i = 1, ncnew * nrnew
         rwork3(i) = rwork1(i)
      enddo
*
997   continue
*
      go to 995
*
*-------------------------------------------------------------------------------
*
*
994   call center(rwork1,rwork3,nc,nr,ndim,rworkcol1,ncol,nrow,ierr)
*
      ncnew = ncol
      nrnew = nrow
*
*
*###############################################################################
*
*
      if(ierr.eq.1) go to 9999
*
* rwork3 contains the output grid, which is back in the proper orientation,
* with width = ncol and height = nrow
*
*
995   dxorgnew = dble ( xo )
      dyorgnew = dble ( yo )
*
c      print *,'Sub: ftprep #4: ncnew, nrnew, dxorgnew, dyorgnew'
c      print *,ncnew, nrnew, dxorgnew, dyorgnew
*
*
*=======================================================================
*
*
9999  return
      end
*
*
************************************************************************
*
      subroutine extend ( nc, nr, nrow, ncol, zgridin, zgridout, g, aa,
     +                    v, ierr )
*
************************************************************************
*
c     Extend grid to the right using prediction filters
c
      implicit none
      integer nc, nr, nrow, ncol, i, j, j2, iwrite, ierr
      real zgridin(nc*nr), zgridout(ncol*nr), g(ncol), aa(ncol), v(ncol)
      real we(16), wesav(16)
*
      ierr = 0
      we = 0.0
      wesav = 0.0
*
c Rows 1 to nrow
      do j2 = 1, nrow
         call getrow ( j2, nc, zgridin, g )
         call predfilt ( nc, ncol, g, 16, we, aa, v )
         do i = 1, 16
            wesav(i) = wesav(i) + we(i)
         end do
      end do
*
c Average prediction filter coefficients
      do i = 1, 16
         wesav(i) = wesav(i) / real(nrow)
      end do
*
c Apply filter to all rows
      do j = 1, nr
         call getrow ( j, nc, zgridin, g )
         iwrite = 0
         if ( j == 10 ) iwrite = 1
         call filter ( nc, ncol, g, 16, wesav, aa, iwrite )
         call putrow ( j, ncol, g, zgridout )
      end do
*
      return
      end
*
*
************************************************************************
*
      subroutine filter ( n, lc, a, nco, we, aa, iwrite )
*
************************************************************************
*
c     Extend profile using prediction filter coefficients with cosine weighting
c
      implicit none
      integer n, lc, nco, iwrite
      real a(lc), we(nco), aa(lc)
      integer i, j, lp
      real w, w1, pi
      parameter ( pi = 3.1415927 )
*
      lp = n + 1
      do i = 1, n
         aa(i) = a(n-i+1)
      end do
*
      do i = lp, lc
         a(i) = 0.0
         aa(i) = 0.0
         do j = 2, nco
            a(i) = a(i) - we(j) * a(i-j+1)
            aa(i) = aa(i) - we(j) * aa(i-j+1)
         end do
         w1 = pi * real(i-lp) / real(lc-lp)
         w = (cos(w1) + 1.0) / 2.0
         a(i) = w * a(i) + (1.0 - w) * aa(lc-i+lp)
      end do
*
      return
      end
*
*
************************************************************************
*
      subroutine predfilt ( n, lc, a, nco, we, aa, v )
*
************************************************************************
*
c     Compute prediction filter coefficients using Burg algorithm for MEM extrapolation
c
      implicit none
      integer n, lc, nco
      real a(lc), we(nco), aa(lc), v(lc)
      integer i, j, jh, k
      real ap, xind, rc, temp_val
*
      do i = 1, n
         aa(i) = a(i)
         v(i) = a(i)
      end do
      we(1) = 1.0
      do j = 2, nco
         ap = 0.0
         xind = 0.0
         do i = j, n
            ap = ap + aa(i) * aa(i) + v(i-j+1) * v(i-j+1)
            xind = xind + aa(i) * v(i-j+1)
         end do
         if ( ap .eq. 0.0 ) then
            rc = 0.0
         else
            rc = -2.0 * xind / ap
         endif
         do i = j, n
            temp_val = aa(i)
            aa(i) = aa(i) + rc * v(i-j+1)
            v(i-j+1) = v(i-j+1) + rc * temp_val
         end do
         we(j) = 0.0
         jh = (j + 1) / 2
         do i = 1, jh
            k = j - i + 1
            temp_val = we(k) + rc * we(i)
            we(i) = we(i) + rc * we(k)
            we(k) = temp_val
         end do
      end do
      return
      end
*
*
************************************************************************
*
      subroutine smooth ( zgridin, zgridout, nr, ndim, g, h, nc, nrow, 
     +                    ierr )
*
************************************************************************
*
c     Smooth grid across rows
c
      implicit none
      integer nr, ndim, nc, nrow, ierr
      real zgridin(nc*nrow), zgridout(nc*nrow)
      real g(ndim), h(ndim)
      integer i, j, k, ii, kk
      real g1, g2
      real DUMMY
      parameter ( DUMMY = -1.0e32 )
*
      ierr = 0
      do j = 1, nr
         call getrow ( j, nc, zgridin, g )
         call putrow ( j, nc, g, zgridout )
      end do
*
      do j = nr + 1, nrow
         call getrow ( j, nc, zgridin, g )
         k = min(j-nr, nrow-j+1)
         do i = 1, k
            h(i) = g(i)
            do ii = 1, k
               kk = mod(nc+i-ii, nc)
               if ( kk .eq. 0 ) kk = nc
               if ( g(i+ii) .eq. DUMMY ) then
                  g1 = 0.0
               else
                  g1 = g(i+ii)
               endif
               if ( g(kk) .eq. DUMMY ) then
                  g2 = 0.0
               else
                  g2 = g(kk)
               endif
               h(i) = h(i) + real(k-ii+1) * (g1 + g2) / real(k+1)
            end do
            h(i) = h(i) / real(k+1)
         end do
         do i = k + 1, nc - k
            if ( g(i) .eq. DUMMY ) then
               h(i) = 0.0
            else
               h(i) = g(i)
            endif
            do ii = 1, k
               if ( g(i+ii) .eq. DUMMY ) then
                  g1 = 0.0
               else
                  g1 = g(i+ii)
               endif
               if ( g(i-ii) .eq. DUMMY ) then
                  g2 = 0.0
               else
                  g2 = g(i-ii)
               endif
               h(i) = h(i) + real(k-ii+1) * (g1 + g2) / real(k+1)
            end do
            h(i) = h(i) / real(k+1)
         end do
         do i = nc - k + 1, nc
            if ( g(i) .eq. DUMMY ) then
               h(i) = 0.0
            else
               h(i) = g(i)
            endif
            do ii = 1, k
               kk = mod(i+ii, nc)
               if ( kk .eq. 0 ) kk = nc
               if ( g(i-ii) .eq. DUMMY ) then
                  g1 = 0.0
               else
                  g1 = g(i-ii)
               endif
               if ( g(kk) .eq. DUMMY ) then
                  g2 = 0.0
               else
                  g2 = g(kk)
               endif
               h(i) = h(i) + real(k-ii+1) * (g1 + g2) / real(k+1)
            end do
            h(i) = h(i) / real(k+1)
         end do
         call putrow ( j, nc, h, zgridout )
      end do
      return
      end
*
*
************************************************************************
*
      subroutine center ( zgridin, zgridout, nc, nr, ndim, g, ncol,
     +                    nrow, ierr )
*
************************************************************************
*
c     Center the original data in the output grid
c
      implicit none
      integer nc, nr, ndim, ncol, nrow, ierr
      real zgridin(ncol*nrow), zgridout(ncol*nrow)
      real g(ndim)
      integer i, j, k, jw, nx, ny
      real temp
*
      ierr = 0
      nx = (ncol - nc) / 2
      ny = (nrow - nr) / 2
*
      jw = 0
*
      do j = nrow - ny + 1, nrow
         jw = jw + 1
         call getrow ( j, ncol, zgridin, g )
         do i = 1, nx
            temp = g(ncol)
            do k = ncol, 2, -1
               g(k) = g(k-1)
            end do
            g(1) = temp
         end do
         call putrow ( jw, ncol, g, zgridout )
      end do
*
      do j = 1, nrow - ny
         jw = jw + 1
         call getrow ( j, ncol, zgridin, g )
         do i = 1, nx
            temp = g(ncol)
            do k = ncol, 2, -1
               g(k) = g(k-1)
            end do
            g(1) = temp
         end do
         call putrow ( jw, ncol, g, zgridout )
      end do
*
      return
      end
*
*
************************************************************************
*
      subroutine check ( rwork1, nc, nr, g, ierr )
*
************************************************************************
*
c     Check input grid for nulls and compute min/max
c
      implicit none
      integer nc, nr, ierr
      real rwork1(nc*nr)
      real g(nc)
      real gmin, gmax
      common /extreme/ gmin,gmax
      integer i, j
*
      ierr = 0
      gmin = 1.0e38
      gmax = -1.0e38
*
      do j = 1, nr
         call getrow ( j, nc, rwork1, g )
         do i = 1, nc
            if ( g(i) .lt. gmin ) gmin = g(i)
            if ( g(i) .gt. gmax ) gmax = g(i)
         end do
      end do
      return
      end
*
*
************************************************************************
*
      subroutine bain_transpose ( zgridin, zgridout, g, ndim, nc, nr, 
     +                            ierr )
*
************************************************************************
*
c     Transpose grid
c
      implicit none
      integer nc, nr, ndim, ierr
      real zgridin(nc*nr), zgridout(nr*nc)
      real g(ndim)
      integer i, j, ipt_out
*
      ierr = 0
      do j = 1, nr
         call getrow ( j, nc, zgridin, g )
         do i = 1, nc
            ipt_out = (i-1)*nr + j
            zgridout(ipt_out) = g(i)
         end do
      end do
      return
      end
*
*
************************************************************************
*
      subroutine getrow ( j, nc, zgrid, g )
*
************************************************************************
*
c     Get row from grid array
c
      implicit none
      integer j, nc
      real zgrid(*)
      real g(nc)
      integer i, ipt
*
      do i = 1, nc
         ipt = (j-1)*nc + i
         g(i) = zgrid(ipt)
      end do
      return
      end
*
*
************************************************************************
*
      subroutine putrow ( j, nc, g, zgrid )
*
************************************************************************
*
c     Put row into grid array
c
      implicit none
      integer j, nc
      real g(nc)
      real zgrid(*)
      integer i, ipt
*
      do i = 1, nc
         ipt = (j-1)*nc + i
         zgrid(ipt) = g(i)
      end do
      return
      end
*
*
************************************************************************
      subroutine prime5 ( n, nprime )
      integer n, nprime
      integer i
      logical found
*
      nprime = n
      found = .false.
      do while ( .not. found )
         nprime = nprime + 1
         found = .true.
         do i = 2, 4
            if ( mod(nprime,i) .eq. 0 ) then
               found = .false.
            endif
         enddo
      enddo
      return
      end
*
*