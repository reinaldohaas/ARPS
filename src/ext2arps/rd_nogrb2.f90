!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDNMCGRB2                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,   &
           iproj_grb,nx_grb,ny_grb,dx_grb,dy_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, maxvar, nzsoilin_ext,                          &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl, var3dsoil, &
           var_grb2d, var_grb3d, lvldbg,                                &
           iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is a dummy subroutine for case when NCEP GRIB2 library
!  is not linked.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/01/2006
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE


  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext

  INTEGER, INTENT(IN) :: grbflen

  CHARACTER(LEN=*), INTENT(IN) :: gribfile
  CHARACTER(LEN=*), INTENT(IN) :: gribtime

  INTEGER  :: iproj_grb    ! Map projection indicator
                                       ! Already converted to ARPS map definitions
  INTEGER  :: nx_grb       ! Number of points along x-axis
  INTEGER  :: ny_grb       ! Number of points along y-axis
  REAL  :: dx_grb          ! x-direction increment or grid length
  REAL  :: dy_grb          ! y-direction increment or grid length
  REAL  :: latsw           ! Latitude  of South West corner point
  REAL  :: lonsw           ! Longitude of South West corner point
  REAL  :: lattru1         ! Latitude (1st) at which projection is true
  REAL  :: lattru2         ! Latitude (2nd) at which projection is true
  REAL  :: lontrue         ! Longitude      at which projection is true
  INTEGER  :: uvearth         ! = 0, Resolved u and v components of vector quantities relative to easterly and northerly directions
                                          ! = 1, Resolved u and v components of vector quantities relative to the defined grid in the direction of increasing x and y (or i and j) coordinates, respectively

  INTEGER, INTENT(IN)  :: n2dvs, n3dvs, maxvar
  INTEGER, INTENT(IN)  :: nzsoilin_ext, lvldbg
  INTEGER, INTENT(IN)  :: varids(4,maxvar)
  INTEGER, INTENT(IN)  :: var2dindx(maxvar), var3dindx(maxvar)
  REAL,    INTENT(IN)  :: var2dlvl(n2dvs), var3dlvl(nz_ext), var3dsoil(nzsoilin_ext)
  REAL  :: var_grb2d(nx_ext,ny_ext,n2dvs)
  REAL  :: var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs)

  INTEGER  :: iret         ! Return flag
!
!-----------------------------------------------------------------------
!
!  Temporary arrays to read GRIB file
!
!-----------------------------------------------------------------------
!

!  REAL    :: rcdata(nx_ext*ny_ext)           ! temporary data array
!  REAL    :: var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt)        ! GRIB 2-D variables
!  REAL    :: var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt) ! GRIB 3-D variables
!  INTEGER :: var_lev3d(nz_ext,n3dvs,n3dlvt)
                                       ! Levels (hybrid) for each 3-D variable

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  WRITE(6,'(/,2(a,/),/,a,/,/,3(a,/))')                                  &
  'WARNING: GRIB2 libary is not linked. To link NCEP GRIB2 library,',   &
  '         you should compile the program as:',                        &
  '         $> makearps -io grib2 [Other_options] <arps_program>',      &
  '   NOTE: You must first install libraries libjasper.a(JPEG2000),',   &
  '         libpng.a(PNG) and libz.a. See src/external/g2lib/README',   &
  '         for more details.'

  iret = -2

  RETURN
END SUBROUTINE rdnmcgrb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDGRB2DIMS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdgrb2dims(gribfile,grbflen,nx_ext,ny_ext,iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is a dummy subroutine for case when NCEP GRIB2 library
!  is not linked.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/18/2007
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE


  INTEGER  :: nx_ext,ny_ext

  INTEGER, INTENT(IN) :: grbflen

  CHARACTER(LEN=*), INTENT(IN) :: gribfile

  INTEGER  :: iret         ! Return flag

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  WRITE(6,'(/,2(a,/),/,a,/,/,3(a,/))')                                  &
  'WARNING: GRIB2 libary is not linked. To link NCEP GRIB2 library,',   &
  '         you should compile the program as:',                        &
  '         $> makearps -io grib2 [Other_options] <arps_program>',      &
  '   NOTE: You must first install libraries libjasper.a(JPEG2000),',   &
  '         libpng.a(PNG) and libz.a. See src/external/g2lib/README',   &
  '         for more details.'

  iret = -2

  RETURN
END SUBROUTINE rdgrb2dims
