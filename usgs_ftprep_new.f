************************************************************************
*
      subroutine usgs_ftprep_new ( pct, ncolin, nrowin, dxorg, dyorg,
     + ddx, ddy,
     + ncnew, nrnew, dxorgnew, dyorgnew, rwork1, rwork2, rwork3,
     + rworkcol1, rworkcol2, rworkcol3, maxx, ierr )
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
*===============================================================================
*
* Set new grid geometry
*
* Bain (27Feb21) - I put in a rounding factor after having a resolution
* issue with pct being single precision (? - seems strange to need this...)
*
      rfactor = 1.0+((pct/100.0)+(pct/100000.0))
c        print *,'pct, nc, nr,rfactor: ',pct,nc,nr,rfactor
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
*===============================================================================
*
      if ( pct .lt. 0.00000001 ) then
         ncol = nc1
         nrow = nr1
      else
         call prime5(nc1,ncol)
         call prime5(nr1,nrow)
      endif
*
*===============================================================================
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
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
      call bain_transpose ( rwork2, rwork3, rworkcol1, ndim,
     +                      ncol, nr, ierr)
      if(ierr.eq.1) go to 9999
*
c      print *,'Just after call transpos #1, ndim, ncol',ndim, ncol, nr
*
* This has each column in the original grid turned on its side, and put
* out to rwork3 array.  So dimensions of this are now (nr, ncol),
* i.e., transposed version of the column-extended grid
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
         dz(i) = dble ( rwork3(i) )
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
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*
*=======================================================================
*
* Bain (16May19): extend grid on right side, which is now the top of the
* rows of the original grid, including smoothing before (as of 23Aug22)
*
* Didn't really help for Trois test data set
*
*=======================================================================
*
c      call smooth ( rwork3, rwork2, ncol, ndim, rworkcol1, rworkcol2,
c     +              nr, ncol, ierr )
c      if(ierr.eq.1) go to 9999
*
* Bain (23Aug22): to keep things below consistent, copy rwork2 back in to rwork1
*
c      do i = 1, ncol * nr
c         rwork3(i) = rwork2(i)
c      enddo
*
c      print *,'Just before call extend #2'
*
      call extend (nr, ncol, nc, nrow, rwork3, rwork1, rworkcol1,
     +             rworkcol2, rworkcol3, ierr)
      if(ierr.eq.1) go to 9999
*
c      print *,'Just after call extend #2'
*
* grid in rwork1 is now transposed, with width = nrow, and height = ncol
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
         dz(i) = dble ( rwork1(i) )
      enddo
*
      call write_geosoft_grid ( gridout, ncnew, nrnew, dxorg,
     +                          dyorg, ddx, ddy, dznull, dz,
     +                          drot, maxz, grid_type, clabel, dwork,
     +                          rwork, istat )
      stop
*
830   continue
*
*
*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*
*=======================================================================
*
* Bain (16May19): smooth grid
*
*=======================================================================
*
*
c      print *,'Just before call smooth #1'
*
* Note: grid is transposed, so reversing the ncol and nrow values
*
      call smooth ( rwork1, rwork2, nc, ndim, rworkcol1, rworkcol2,
     +              nrow, ncol, ierr )
c      print *,'Just after call smooth #1'
*
      if(ierr.eq.1) go to 9999
*
*
*=======================================================================
*
* Bain (16May19): Transpose back
*
*=======================================================================
*
c      print *,'Just before call transpos #2'
*
* This grid is transposed - so ncol,nrow becomes nrow,ncol for this call
*
996   continue
*
      call bain_transpose ( rwork2,rwork1,rworkcol1,ndim,
     +                      nrow,ncol,ierr)
      if(ierr.eq.1) go to 9999
*
c      print *,'Just after call transpos #2'
*
* Bain (15Aug22): jump to 998 to test
*
c      go to 998
*
*
* rwork1 is now transposed back to have width = ncol and height = nrow
*
*
*=======================================================================
*
* Bain (16May19): smooth grid again
*
*=======================================================================
*
* Note: grid is no longer transposed, so going back to having ncol, nrow (proper order)
*
c      print *,'Just before call smooth #2'
      call smooth ( rwork1, rwork2, nr, ndim, rworkcol1, rworkcol2,
     +              ncol, nrow, ierr )
      if(ierr.eq.1) go to 9999
c      print *,'Just after call smooth #2'
*
*
* rwork2 now has smoothing in both directions
*
*
*=======================================================================
*
* Bain (16May19): center data
*
*=======================================================================
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
* Bain (15Aug22): test copying rwork2 to 3, instead of next call (which is now commented
*
998   do i = 1, ncnew * nrnew
         rwork3(i) = rwork2(i)
      enddo
*
997   continue
*
      go to 995
*
*-------------------------------------------------------------------------------
*
*
994   call center(rwork2,rwork3,nc,nr,ndim,rworkcol1,ncol,nrow,ierr)
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
