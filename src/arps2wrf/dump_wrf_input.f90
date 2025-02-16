!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE define_wrf_variables          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE define_wrf_variables_V2(ncid,file_type,wrfversion,           &
                     nx,ny,nz,nzsoil,spec_bdy_width,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Define NetCDF variables for input file or boundary file of
!   WRF version 2.2
!
! NOTE:
!   This subroutine only defines the variables to be written in NetCDF
!   define mode. Please make sure to write the actual data in arps2wrf.f90.
!   because it set NF_SET_FILL to NF_NOFILL.
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: file_type             ! 1 for WRF input file
                                                ! 2 for WRF boundary file
  REAL,    INTENT(IN)  :: wrfversion
  INTEGER, INTENT(IN)  :: nx
  INTEGER, INTENT(IN)  :: ny
  INTEGER, INTENT(IN)  :: nz
  INTEGER, INTENT(IN)  :: nzsoil
  INTEGER, INTENT(IN)  :: spec_bdy_width
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER  :: dimunlim_id
  INTEGER  :: dimx_id, dimxs_id
  INTEGER  :: dimy_id, dimys_id
  INTEGER  :: dimz_id, dimzs_id
  INTEGER  :: dimzsoil_id
  INTEGER  :: dimstr_id
  INTEGER  :: dimext_id
  INTEGER  :: dimbdy_id
  INTEGER  :: dimsoil_id, dimland_id

  INTEGER  :: dimlenx, dimleny, dimlenz, dimlenzsoil
  INTEGER  :: dimlenstr, dimlenext
  INTEGER  :: dimlenbdy
  INTEGER  :: dimlensoil, dimlenland
  INTEGER  :: dimlandcat_id, dimsoilcat_id

  INTEGER  :: varid, oldfillmode

  TYPE(wrf_var_metadata) :: var_meta

!                                            West  East  South North
  CHARACTER(LEN=2), PARAMETER :: appd(4) = (/'XS', 'XE', 'YS', 'YE'/)
  CHARACTER(LEN=10)           :: varname
  INTEGER, DIMENSION(4)       :: dimids

  INTEGER :: n, nq

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code below ......
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! get dimension ids and dimension length
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'Time',dimunlim_id)

  istatus = NF_INQ_DIMID(ncid,'DateStrLen',dimstr_id)
  istatus = NF_INQ_DIMLEN(ncid,dimstr_id,dimlenstr)

  istatus = NF_INQ_DIMID(ncid,'west_east_stag',dimx_id)
  istatus = NF_INQ_DIMID(ncid,'west_east',     dimxs_id)

  istatus = NF_INQ_DIMLEN(ncid,dimx_id,dimlenx)
  IF(dimlenx /= nx) THEN
    WRITE(6,*) 'Mismatched dimension size in X direction.'
    STOP
  END IF

  istatus = NF_INQ_DIMID(ncid,'south_north_stag',dimy_id)
  istatus = NF_INQ_DIMID(ncid,'south_north',     dimys_id)

  istatus = NF_INQ_DIMLEN(ncid,dimy_id,dimleny)
  IF(dimleny /= ny) THEN
    WRITE(6,*) 'Mismatched dimension size in Y direction.'
    STOP
  END IF

  istatus = NF_INQ_DIMID(ncid,'bottom_top_stag',dimz_id)
  istatus = NF_INQ_DIMID(ncid,'bottom_top',     dimzs_id)

  istatus = NF_INQ_DIMLEN(ncid,dimz_id,dimlenz)
  IF(dimlenz /= nz) THEN
    WRITE(6,*) 'Mismatched dimension size in the 3rd dimension.'
    STOP
  END IF

  IF (file_type == 1) THEN        ! input file
    IF (wrfversion < 3.0) THEN
      istatus = NF_INQ_DIMID(ncid,'DIM0008',dimbdy_id)   !% version 2.2
    ELSE
      istatus = NF_INQ_DIMID(ncid,'DIM0009',dimbdy_id)   !% version 3.0
    END IF
    istatus = NF_INQ_DIMLEN(ncid,dimbdy_id,dimlenbdy)

    istatus = NF_INQ_DIMID(ncid,'soil_layers_stag',dimzsoil_id)
    istatus = NF_INQ_DIMLEN(ncid,dimzsoil_id,dimlenzsoil)
  ELSE                           ! boundary file
    istatus = NF_INQ_DIMID(ncid,'bdy_width',dimbdy_id)
    istatus = NF_INQ_DIMLEN(ncid,dimbdy_id,dimlenbdy)
  END IF

  IF (wrfversion >= 3.1) THEN
    istatus = NF_INQ_DIMID(ncid,'land_cat_stag',dimlandcat_id)
    istatus = NF_INQ_DIMID(ncid,'soil_cat_stag',dimsoilcat_id)
  END IF

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
!   Defines Times string
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_VAR(ncid,'Times',NF_CHAR,2,                     &
                         (/dimstr_id, dimunlim_id/),varid)

  IF (file_type == 1) THEN                ! WRF input file

      var_meta%name        = 'LU_INDEX'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LAND USE CATEGORY'
      var_meta%units       = ' '
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!-------------------------------------------------------------------
!
!  Defines U  and V
!
!-------------------------------------------------------------------
!
      var_meta%name        = 'U'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'x-wind component'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = 'X'
      var_meta%coordinates = 'XLONG_U XLAT_U'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimx_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'V'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'y-wind component'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = 'Y'
      var_meta%coordinates = 'XLONG_V XLAT_V'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimy_id,dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!-----------------------------------------------------------------------
!
!  Define vertical velocity
!
!-----------------------------------------------------------------------
      var_meta%name        = 'W'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'z-wind component'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = 'Z'
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimys_id,dimz_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define perturbation and base state geopotential
!
!----------------------------------------------------------------------
      var_meta%name        = 'PH'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'perturbation geopotential'
      var_meta%units       = 'm2 s-2'
      var_meta%stagger     = 'Z'
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimys_id,dimz_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'PHB'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'base-state geopotential'
      var_meta%units       = 'm2 s-2'
      var_meta%stagger     = 'Z'
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimys_id,dimz_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define perturbation potential temperature (theta - t0)
!
!----------------------------------------------------------------------
      var_meta%name        = 'T'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'perturbation potential temperature (theta-t0)'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'T_INIT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'initial potential temperature'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define air mass
!
!----------------------------------------------------------------------
      var_meta%name        = 'MU'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'perturbation dry air mass in column'
      var_meta%units       = 'Pa'
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'MUB'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'base state dry air mass in column'
      var_meta%units       = 'Pa'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      IF (wrfversion < 3.1) THEN
        var_meta%name        = 'MU0'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'initial dry mass in column'
        var_meta%units       = 'Pa'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      ELSE
        var_meta%name        = 'P'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'perturbation pressure'
        var_meta%units       = 'Pa'
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                         (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'PB'
        var_meta%description = 'BASE STATE PRESSURE'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                         (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END IF

      var_meta%name        = 'SR'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'fraction of frozen precipitation'
      var_meta%units       = '-'
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define  Vertical mass coordinate parameters
!
!----------------------------------------------------------------------

      var_meta%coordinates = ''

      var_meta%name        = 'FNM'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'upper weight for vertical stretching'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'FNP'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'lower weight for vertical stretching'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      !  Store the vertical coordinate information
      var_meta%name        = 'RDNW'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'inverse d(eta) values between full (w) levels'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'RDN'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'inverse d(eta) values between half (mass) levels'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'DNW'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'd(eta) values between full (w) levels'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'DN'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'd(eta) values between half (mass) levels'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'ZNU'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'eta values on half (mass) levels'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'ZNW'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'eta values on full (w) levels'
      var_meta%units       = ' '
      var_meta%stagger     = 'Z'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                        (/dimz_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'CFN'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'extrapolation constant'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'CFN1'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'extrapolation constant'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!-----------------------------------------------------------------------
!
! Misc. variables
!
!-----------------------------------------------------------------------

!      var_meta%name        = 'EPSTS'
!      var_meta%fieldType   = 104
!      var_meta%memoryOrder = '0  '
!      var_meta%description = 'leapfrog time filter constant'
!      var_meta%units       = ' '
!      var_meta%stagger     = ' '
!
!      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
!                                        (/dims_id,dimunlim_id/),varid)
!      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'STEP_NUMBER'
      var_meta%fieldType   = 106
      var_meta%memoryOrder = '0  '
      var_meta%description = 'step number'
      var_meta%units       = '-'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define  Q2, T2, TH2, PSFC, U10, V10 etc.
!
!----------------------------------------------------------------------

      var_meta%coordinates = 'XLONG XLAT'

      var_meta%name        = 'Q2'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'QV at 2 M'
      var_meta%units       = 'kg kg-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'T2'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'TEMP at 2 M'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'TH2'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'POT TEMP at 2 M'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'PSFC'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'SFC PRESSURE'
      var_meta%units       = 'Pa'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'U10'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'U at 10 M'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'V10'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'V at 10 M'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      IF (wrfversion < 3.1) THEN

        var_meta%name        = 'URATX'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Ratio of U over U10 on mass points'
        var_meta%units       = 'dimensionless'
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'VRATX'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Ratio of V over V10 on mass points'
        var_meta%units       = 'dimensionless'
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'TRATX'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Ratio of T over TH2 on mass points'
        var_meta%units       = 'dimensionless'
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      END IF

!-----------------------------------------------------------------------

      var_meta%coordinates = ''

      var_meta%name        = 'RDX'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'INVERSE X GRID LENGTH'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'RDY'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'INVERSE Y GRID LENGTH'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'DTS'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'SMALL TIMESTEP'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'DTSEPS'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'TIME WEIGHT CONSTANT FOR SMALL STEPS'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'RESM'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'TIME WEIGHT CONSTANT FOR SMALL STEPS'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'ZETATOP'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'ZETA AT MODEL TOP'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'CF1'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = '2nd order extrapolation constant'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'CF2'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = '2nd order extrapolation constant'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'CF3'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = '2nd order extrapolation constant'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define water vapor mixing ratio
!
!  From WRFV2.0 Registry.EM
!package   passiveqv     mp_physics==0       -       moist:qv
!package   kesslerscheme mp_physics==1       -       moist:qv,qc,qr
!package   linscheme     mp_physics==2       -       moist:qv,qc,qr,qi,qs,qg
!package   wsm3scheme    mp_physics==3       -       moist:qv,qc,qr
!package   wsm5scheme    mp_physics==4       -       moist:qv,qc,qr,qi,qs
!package   etampnew      mp_physics==5       -       moist:qv,qc
!package   wsm6scheme    mp_physics==6       -       moist:qv,qc,qr,qi,qs,qg
!package   thompson      mp_physics==8       -       moist:qv,qc,qr,qi,qs,qg;scalar:qni
!package   ncepcloud3    mp_physics==98      -       moist:qv,qc,qr
!package   ncepcloud5    mp_physics==99      -       moist:qv,qc,qr,qi,qs
!
!  At present, we only support, mp_physics = 0,1-6,8
!
!----------------------------------------------------------------------

      var_meta%coordinates = 'XLONG XLAT'

      IF (P_QV_wrf > 0) THEN
!     All microphysics packages need qv
      var_meta%name        = 'QVAPOR'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'Water vapor mixing ratio'
      var_meta%units       = 'kg kg-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)
      END IF

!----------------------------------------------------------------------
!
!  Define QCLOUD and QRAIN, QSNOW, QICE, QGRAUP
!
!----------------------------------------------------------------------


      DO nq = 1, nscalar_wrf

        var_meta%name        = TRIM(qnames_wrf(nq))
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = TRIM(qdescp_wrf(nq))
        var_meta%units       = 'kg kg-1'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      IF (P_QT_wrf > 0) THEN
        var_meta%name        = 'QT'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'Total condensate mixing ratio'
        var_meta%units       = 'kg kg-1'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END IF

      IF (P_QNI_wrf > 0) THEN
      ! Define QNICE for MP package 8

        var_meta%name        = 'QNICE'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'Ice Number concentration'
        var_meta%units       = '  kg(-1)'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END IF

!----------------------------------------------------------------------
!
!  Define FCX, GCX
!
!----------------------------------------------------------------------
      var_meta%coordinates = ''

      var_meta%name        = 'FCX'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'C  '
      var_meta%description = 'RELAXATION TERM FOR BOUNDARY ZONE'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                       (/dimbdy_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'GCX'
      var_meta%description = '2ND RELAXATION TERM FOR BOUNDARY ZONE'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                       (/dimbdy_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define DTBC (scalar)
!
!----------------------------------------------------------------------
      var_meta%name        = 'DTBC'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'TIME SINCE BOUNDARY READ'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                       (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%coordinates = 'XLONG XLAT'

      var_meta%name        = 'LANDMASK'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LAND MASK (1 FOR LAND, 0 FOR WATER)'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                            (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SHDMAX'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'ANNUAL MAX VEG FRACTION'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SHDMIN'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'ANNUAL MIN VEG FRACTION'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                            (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SNOALB'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'ANNUAL MAX SNOW ALBEDO IN FRACTION'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                            (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      IF (wrfversion >= 3.1) THEN
        var_meta%name        = 'LANDUSEF'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'LANDUSE FRACTION BY CATEGORY'
        var_meta%units       = ''
        var_meta%stagger     = 'Z'
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                         (/dimxs_id,dimys_id,dimlandcat_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'SOILCTOP'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'SOIL CAT FRACTION (TOP)'
        var_meta%units       = ' '
        var_meta%stagger     = 'Z'
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                         (/dimxs_id,dimys_id,dimsoilcat_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'SOILCBOT'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'SOIL CAT FRACTION (BOTTOM)'
        var_meta%units       = ' '
        var_meta%stagger     = 'Z'
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                         (/dimxs_id,dimys_id,dimsoilcat_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      END IF
!----------------------------------------------------------------------
!
!  Define soil layers and variables
!
!----------------------------------------------------------------------

      var_meta%name        = 'TSLB'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'SOIL TEMPERATURE'
      var_meta%units       = 'K'
      var_meta%stagger     = 'Z'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                    (/dimxs_id,dimys_id,dimzsoil_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'ZS'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z  '
      var_meta%description = 'DEPTHS OF CENTERS OF SOIL LAYERS'
      var_meta%units       = 'm'
      var_meta%stagger     = 'Z'
      var_meta%coordinates = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                       (/dimzsoil_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'DZS'
      var_meta%fieldType   = 104
      var_meta%description = 'THICKNESSES OF SOIL LAYERS'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                       (/dimzsoil_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SMOIS'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'SOIL MOISTURE'
      var_meta%units       = 'm3 m-3'
      var_meta%stagger     = 'Z'
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,         &
                    (/dimxs_id,dimys_id,dimzsoil_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SH2O'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'SOIL LIQUID WATER'
      var_meta%units       = 'm3 m-3'
      var_meta%stagger     = 'Z'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                  (/dimxs_id,dimys_id,dimzsoil_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      IF (wrfversion >= 3.1) THEN
        var_meta%name        = 'SEAICE'
      ELSE
        var_meta%name        = 'XICE'
      END IF
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'SEA ICE FLAG'
      var_meta%units       = ' '
      var_meta%stagger     = ' '
      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)


!----------------------------------------------------------------------
!
!  Define soil type and vegetation type
!
!----------------------------------------------------------------------

      var_meta%name        = 'IVGTYP'
      var_meta%fieldType   = 106
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'DOMINANT VEGETATION CATEGORY'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,3,           &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'ISLTYP'
      var_meta%fieldType   = 106
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'DOMINANT SOIL CATEGORY'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,3,           &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'VEGFRA'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'VEGETATION FRACTION'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define Water equivalent of accumulated snow
!
!----------------------------------------------------------------------

      var_meta%name        = 'SNOW'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'SNOW WATER EQUIVALENT'
      var_meta%units       = 'kg m-2'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SNOWH'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'PHYSICAL SNOW DEPTH'
      var_meta%units       = 'm'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                     (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'RHOSN'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = ' SNOW DENSITY'
      var_meta%units       = 'kg m-3'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                     (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'CANWAT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY'
      var_meta%description = 'CANOPY WATER'
      var_meta%units       = 'kg m-2'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!----------------------------------------------------------------------
!
!  Define SST
!
!----------------------------------------------------------------------

      var_meta%name        = 'SST'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'SEA SURFACE TEMPERATURE'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%coordinates = ''

      var_meta%name        = 'FNDSNOWH'
      var_meta%fieldType   = 106
      var_meta%memoryOrder = '0  '
      var_meta%description = 'SNOWH_LOGICAL'
      var_meta%units       = '-'
      var_meta%stagger     = ''

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,1,         &
                           (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'FNDSOILW'
      var_meta%fieldType   = 106
      var_meta%memoryOrder = '0  '
      var_meta%description = 'SOILW_LOGICAL'
      var_meta%units       = '-'
      var_meta%stagger     = ''

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,1,         &
                           (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%coordinates = 'XLONG XLAT'

!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------

      IF (wrfversion < 3.1) THEN
!       var_meta%name        = 'TOTSWDN'
!       var_meta%fieldType   = 104
!       var_meta%memoryOrder = 'XY '
!       var_meta%description = '-'
!       var_meta%units       = '-'
!       var_meta%stagger     = ' '
!
!       istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
!                             (/dimxs_id,dimys_id,dimunlim_id/),varid)
!       CALL write_var_meta(ncid,varid,var_meta)
!
!       var_meta%name        = 'TOTLWDN'
!       var_meta%fieldType   = 104
!       var_meta%memoryOrder = 'XY '
!       var_meta%description = '-'
!       var_meta%units       = '-'
!       var_meta%stagger     = ' '
!
!       istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
!                             (/dimxs_id,dimys_id,dimunlim_id/),varid)
!       CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'RSWTOA'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'RLWTOA'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CZMEAN'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CFRACL'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CFRACM'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CFRACH'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'ACFRST'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'NCFRST'
        var_meta%fieldType   = 106
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'ACFRCV'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'NCFRCV'
        var_meta%fieldType   = 106
        var_meta%memoryOrder = 'XY'
        var_meta%description = '-'
        var_meta%units       = '-'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,3,       &
                              (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      ELSE

        var_meta%name        = 'LAI'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Leaf area index'
        var_meta%units       = 'area/area'
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'VAR'
        var_meta%description = 'OROGRAPHIC VARIANCE'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CON'
        var_meta%description = 'OROGRAPHIC CONVEXITY'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OA1'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OA2'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OA3'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OA4'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OL1'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OL2'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OL3'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'OL4'
        var_meta%description = 'OROGRAPHIC DIRECTION ASYMMETRY FUNCTION'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      END IF

!-----------------------------------------------------------------------
!
!  Map scale factor
!
!-----------------------------------------------------------------------

!      IF (wrfversion < 3.0) THEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Version before WRFV3.0
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      var_meta%name        = 'MAPFAC_M'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Map scale factor on mass grid'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'MAPFAC_U'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Map scale factor on u-grid'
      var_meta%units       = ' '
      var_meta%stagger     = 'X'
      var_meta%coordinates = 'XLONG_U XLAT_U'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimx_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'MAPFAC_V'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Map scale factor on v-grid'
      var_meta%units       = ' '
      var_meta%stagger     = 'Y'
      var_meta%coordinates = 'XLONG_V XLAT_V'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimy_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! For WRFV3.0 to accomodate global WRF capbality. However, our code
! only work with conformal map projection at present
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (wrfversion >= 3.0) THEN

        var_meta%name        = 'MAPFAC_MX'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Map scale factor on mass grid, x direction'
        var_meta%units       = ' '
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'MAPFAC_MY'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Map scale factor on mass grid, y direction'
        var_meta%units       = ' '
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'MAPFAC_UX'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Map scale factor on u-grid, x direction'
        var_meta%units       = ' '
        var_meta%stagger     = 'X'
        var_meta%coordinates = 'XLONG_U XLAT_U'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimx_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'MAPFAC_UY'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Map scale factor on u-grid, y direction'
        var_meta%units       = ' '
        var_meta%stagger     = 'X'
        var_meta%coordinates = 'XLONG_U XLAT_U'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimx_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'MAPFAC_VX'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Map scale factor on v-grid'
        var_meta%units       = ' '
        var_meta%stagger     = 'Y'
        var_meta%coordinates = 'XLONG_V XLAT_V'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'MAPFAC_VY'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Map scale factor on v-grid, y direction'
        var_meta%units       = ' '
        var_meta%stagger     = 'Y'
        var_meta%coordinates = 'XLONG_V XLAT_V'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'MF_VX_INV'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'Inverse of map scale factor on v-grid, x direction'
        var_meta%units       = ' '
        var_meta%stagger     = 'Y'
        var_meta%coordinates = 'XLONG_V XLAT_V'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CLAT'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'COMPUTATIONAL GRID LATITUDE, SOUTH IS NEGATIVE'
        var_meta%units       = 'degree_north '
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'CLONG'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XY '
        var_meta%description = 'COMPUTATIONAL GRID LONGITUDE, WEST IS NEGATIVE'
        var_meta%units       = 'degree_east '
        var_meta%stagger     = ' '
        var_meta%coordinates = 'XLONG XLAT'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                         (/dimxs_id,dimys_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      END IF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! END of VERSION 3.0
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      var_meta%name        = 'F'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Coriolis sine latitude term'
      var_meta%units       = 's-1'
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'E'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Coriolis cosine latitude term'
      var_meta%units       = 's-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SINALPHA'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Local sine of map rotation'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'COSALPHA'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Local cosine of map rotation'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!---------------------------------------------------------------------
!
!  Terrain height
!
!---------------------------------------------------------------------

      var_meta%name        = 'HGT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'Terrain Height'
      var_meta%units       = 'm'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!-------------------------------------------------------------------
!
! Surface skin temperature TSK.
!
!-------------------------------------------------------------------

      var_meta%name        = 'TSK'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'SURFACE SKIN TEMPERATURE'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!---------------------------------------------------------------------
!
! Unused arrays --  for real case only
!
!---------------------------------------------------------------------

      var_meta%coordinates = ''

      var_meta%name        = 'U_BASE'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z '
      var_meta%description = 'BASE STATE X WIND IN IDEALIZED CASES'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                         (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'V_BASE'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z '
      var_meta%description = 'BASE STATE Y WIND IN IDEALIZED CASES'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                         (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'QV_BASE'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z '
      var_meta%description = 'BASE STATE QV IN IDEALIZED CASES'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                         (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'T_BASE'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z '
      var_meta%description = 'BASE STATE T IN IDEALIZED CASES'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                         (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'Z_BASE'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'Z '
      var_meta%description = 'BASE STATE HEIGHT IN IDEALIZED CASES'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,2,         &
                                         (/dimzs_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'U_FRAME'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0 '
      var_meta%description = 'FRAME X WIND'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                         (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'V_FRAME'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0 '
      var_meta%description = 'FRAME Y WIND'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                         (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!---------------------------------------------------------------------
!
!   Define variables again.
!
!---------------------------------------------------------------------

      var_meta%name        = 'P_TOP'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = '0  '
      var_meta%description = 'PRESSURE TOP OF THE MODEL'
      var_meta%units       = 'Pa'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                        (/dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!-----------------------------------------------------------------------
!
!  WRF V2.0.3 newly added variables, Latitude and Longitude at corners
!
!-----------------------------------------------------------------------
      IF (wrfversion < 3.1) THEN
        !
        ! Latitude, T point
        !
        var_meta%name        = 'LAT_LL_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower left, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UL_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up left, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UR_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up right, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_LR_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower right, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! Latitude, U point
        !
        var_meta%name        = 'LAT_LL_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower left, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UL_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up left, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UR_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up right, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_LR_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower right, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! Latitude, V point
        !
        var_meta%name        = 'LAT_LL_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower left, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UL_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up left, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UR_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up right, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_LR_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower right, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! Latitude, massless point
        !
        var_meta%name        = 'LAT_LL_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower left, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UL_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up left, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_UR_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude up right, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LAT_LR_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'latitude lower right, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! longitude, T point
        !
        var_meta%name        = 'LON_LL_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower left, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UL_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up left, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UR_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up right, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_LR_T'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower right, temp point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! longitude, U point
        !
        var_meta%name        = 'LON_LL_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower left, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UL_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up left, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UR_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up right, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_LR_U'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower right, u point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! longitude, V point
        !
        var_meta%name        = 'LON_LL_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower left, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UL_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up left, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UR_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up right, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_LR_V'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower right, v point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        !
        ! longitude, massless point
        !
        var_meta%name        = 'LON_LL_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower left, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UL_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up left, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_UR_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude up right, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'LON_LR_D'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'longitude lower right, massless point'
        var_meta%units       = 'degrees'
        var_meta%stagger     = ''

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,         &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      END IF

!-----------------------------------------------------------------------
!
!  WRF latitude and longitude
!
!-----------------------------------------------------------------------

      var_meta%name        = 'XLAT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LATITUDE, SOUTH IS NEGATIVE'
      var_meta%units       = 'degree_north'
      var_meta%stagger     = ' '
      var_meta%coordinates = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'XLONG'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LONGITUDE, WEST IS NEGATIVE'
      var_meta%units       = 'degree_east'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'XLAT_U'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LATITUDE, SOUTH IS NEGATIVE'
      var_meta%units       = 'degree_north'
      var_meta%stagger     = 'X'
      var_meta%coordinates = 'XLONG_U XLAT_U'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimx_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'XLONG_U'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LONGITUDE, WEST IS NEGATIVE'
      var_meta%units       = 'degree_east'
      var_meta%stagger     = 'X'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimx_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'XLAT_V'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LATITUDE, SOUTH IS NEGATIVE'
      var_meta%units       = 'degree_north'
      var_meta%stagger     = 'Y'
      var_meta%coordinates = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimy_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'XLONG_V'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LONGITUDE, WEST IS NEGATIVE'
      var_meta%units       = 'degree_east'
      var_meta%stagger     = 'Y'
      var_meta%coordinates = 'XLONG_V XLAT_V'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimy_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

!-------------------------------------------------------------------
!
!  Surface characteristics
!
!------------------------------------------------------------------

      var_meta%coordinates = 'XLONG XLAT'

      var_meta%name        = 'ALBBCK'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'BACKGROUND ALBEDO'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'TMN'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'SOIL TEMPERATURE AT LOWER BOUNDARY'
      var_meta%units       = 'K'
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'XLAND'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'LAND MASK (1 FOR LAND, 2 FOR WATER)'
      var_meta%units       = ' '
      var_meta%stagger     = ' '
      var_meta%coordinates = 'XLONG XLAT'

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,         &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      var_meta%name        = 'SNOWC'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER)'
      var_meta%units       = ' '
      var_meta%stagger     = ' '

      istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,        &
                               (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL write_var_meta(ncid,varid,var_meta)

      IF (wrfversion >= 3.1) THEN
        var_meta%name        = 'T00'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'BASE STATE TEMPERATURE'
        var_meta%units       = 'K'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,      &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'P00'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'BASE STATE PRESSURE'
        var_meta%units       = 'Pa'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,      &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'TLP'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'BASE STATE LAPSE REATE'
        var_meta%units       = 'K'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,      &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

        var_meta%name        = 'TISO'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = '0  '
        var_meta%description = 'TEMP AT WHICH THE BASE T TURNS CONST'
        var_meta%units       = 'K'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,1,      &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)


        var_meta%name        = 'SAVE_TOPO_FROM_REAL'
        var_meta%fieldType   = 106
        var_meta%memoryOrder = '0  '
        var_meta%description = '(1=original topo from real/0=topo modified by WRF) flag if input has original topography from real'
        var_meta%units       = 'flag'
        var_meta%stagger     = ' '

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_INT,1,        &
                                          (/dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)

      END IF

    ELSE IF(file_type == 2 ) THEN    ! boundary file

      istatus = NF_DEF_VAR(ncid,                                       &
          'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_',NF_CHAR,2, &
          (/dimstr_id, dimunlim_id/),varid)
      istatus = NF_DEF_VAR(ncid,                                       &
          'md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_',NF_CHAR,2, &
          (/dimstr_id, dimunlim_id/),varid)

      var_meta%coordinates = ' '

!-----------------------------------------------------------------------
!
!  Defind U related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimys_id    ! West first dimension
      dimids(2) = dimys_id    ! East
      dimids(3) = dimx_id     ! South
      dimids(4) = dimx_id     ! North

      varname              = 'U_B'       ! will change in do loop
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'       ! will change in do loop
      var_meta%description = 'bdy x-wind component'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = 'X'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'U_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy tend x-wind component'
      var_meta%units       = '(m s-1)/dt'
      var_meta%stagger     = 'X'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind V related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimy_id      ! West first dimension
      dimids(2) = dimy_id      ! East
      dimids(3) = dimxs_id     ! South
      dimids(4) = dimxs_id     ! North

      varname              = 'V_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy y-wind component'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = 'Y'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'V_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy tend y-wind component'
      var_meta%units       = '(m s-1)/dt'
      var_meta%stagger     = 'Y'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind W related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimys_id      ! West first dimension
      dimids(2) = dimys_id      ! East
      dimids(3) = dimxs_id      ! South
      dimids(4) = dimxs_id      ! North

      varname              = 'W_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy z-wind component'
      var_meta%units       = 'm s-1'
      var_meta%stagger     = 'Z'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                  (/dimids(n),dimz_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'W_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy tend z-wind component'
      var_meta%units       = '(m s-1)/dt'
      var_meta%stagger     = 'Z'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                  (/dimids(n),dimz_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind PH related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimys_id      ! West first dimension
      dimids(2) = dimys_id      ! East
      dimids(3) = dimxs_id      ! South
      dimids(4) = dimxs_id      ! North

      varname              = 'PH_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy perturbation geopotential'
      var_meta%units       = 'm2 s-2'
      var_meta%stagger     = 'Z'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimz_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'PH_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy tend perturbation geopotential'
      var_meta%units       = '(m2 s-2)/dt'
      var_meta%stagger     = 'Z'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimz_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind T related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimys_id      ! West first dimension
      dimids(2) = dimys_id      ! East
      dimids(3) = dimxs_id      ! South
      dimids(4) = dimxs_id      ! North

      varname              = 'T_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy perturbation potential temperature (theta-t0)'
      var_meta%units       = 'K'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'T_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'bdy tend perturbation potential temperature (theta-t0)'
      var_meta%units       = '(K)/dt'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind MU related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimys_id      ! West first dimension
      dimids(2) = dimys_id      ! East
      dimids(3) = dimxs_id      ! South
      dimids(4) = dimxs_id      ! North

      varname        = 'MU_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'bdy perturbation dry air mass in column'
      var_meta%units       = 'Pa'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                    (/dimids(n),dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'MU_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XY '
      var_meta%description = 'bdy tend perturbation dry air mass in column'
      var_meta%units       = '(Pa)/dt'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                    (/dimids(n),dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind QV related boundary variables
!
!-----------------------------------------------------------------------

      dimids(1) = dimys_id      ! West first dimension
      dimids(2) = dimys_id      ! East
      dimids(3) = dimxs_id      ! South
      dimids(4) = dimxs_id      ! North

      varname              = 'QVAPOR_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'Water vapor mixing ratio'
      var_meta%units       = 'kg kg-1'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'QVAPOR_BT'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XYZ'
      var_meta%description = 'Water vapor mixing ratio'
      var_meta%units       = '(kg kg-1)/dt'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)//'Z'

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                    (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

!-----------------------------------------------------------------------
!
!  Defind QC related boundary variables
!
!-----------------------------------------------------------------------

      DO nq = 1, nscalar_wrf

        dimids(1) = dimys_id      ! West first dimension
        dimids(2) = dimys_id      ! East
        dimids(3) = dimxs_id      ! South
        dimids(4) = dimxs_id      ! North

        varname              = TRIM(qnames_wrf(nq))//'_B'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = TRIM(qdescp_wrf(nq))
        var_meta%units       = 'kg kg-1'
        var_meta%stagger     = ''

        DO n = 1,4
          var_meta%name        = TRIM(varname)//appd(n)
          var_meta%memoryOrder = appd(n)//'Z'

          istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                      (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
          CALL write_var_meta(ncid,varid,var_meta)
        END DO

        varname              = TRIM(qnames_wrf(nq))//'_BT'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = TRIM(qdescp_wrf(nq))
        var_meta%units       = 'kg kg-1'
        var_meta%stagger     = ''

        DO n = 1,4
          var_meta%name        = TRIM(varname)//appd(n)
          var_meta%memoryOrder = appd(n)//'Z'

          istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                      (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
          CALL write_var_meta(ncid,varid,var_meta)
        END DO

      END DO

!-----------------------------------------------------------------------
!
!  Defind QNICE related boundary variables
!
!-----------------------------------------------------------------------

      IF (P_QNI_wrf > 0) THEN
        dimids(1) = dimys_id      ! West first dimension
        dimids(2) = dimys_id      ! East
        dimids(3) = dimxs_id      ! South
        dimids(4) = dimxs_id      ! North

        varname              = 'QNICE_B'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'Ice Number concentration'
        var_meta%units       = ''
        var_meta%stagger     = ''

        DO n = 1,4
          var_meta%name        = TRIM(varname)//appd(n)
          var_meta%memoryOrder = appd(n)//'Z'

          istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                      (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
          CALL write_var_meta(ncid,varid,var_meta)
        END DO

        varname              = 'QNICE_BT'
        var_meta%fieldType   = 104
        var_meta%memoryOrder = 'XYZ'
        var_meta%description = 'Ice Number concentration'
        var_meta%units       = ''
        var_meta%stagger     = ''

        DO n = 1,4
          var_meta%name        = TRIM(varname)//appd(n)
          var_meta%memoryOrder = appd(n)//'Z'

          istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,       &
                      (/dimids(n),dimzs_id,dimbdy_id,dimunlim_id/),varid)
          CALL write_var_meta(ncid,varid,var_meta)
        END DO
      END IF

!-----------------------------------------------------------------------
!
!  Defind bdy Height of orographic shadow
!
!-----------------------------------------------------------------------
    IF (wrfversion >= 3.1) THEN
      dimids(1) = dimys_id      ! West first dimension
      dimids(2) = dimys_id      ! East
      dimids(3) = dimxs_id      ! South
      dimids(4) = dimxs_id      ! North

      varname              = 'HT_SHAD_B'
      var_meta%fieldType   = 104
      var_meta%memoryOrder = 'XS'
      var_meta%description = 'bdy Height of orographic shadow'
      var_meta%units       = 'm'
      var_meta%stagger     = ''

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                             (/dimids(n),dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO

      varname              = 'HT_SHAD_BT'
      var_meta%description = 'bdy tend Height of orographic shadow'
      var_meta%units       = '(m)/dt'

      DO n = 1,4
        var_meta%name        = TRIM(varname)//appd(n)
        var_meta%memoryOrder = appd(n)

        istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,       &
                             (/dimids(n),dimbdy_id,dimunlim_id/),varid)
        CALL write_var_meta(ncid,varid,var_meta)
      END DO
    END IF

  ELSE
    WRITE(6,*) 'Define variables only for the first file.'
    WRITE(6,*) ' file # is ',file_type,' greater than 2.'
    STOP
  END IF

  istatus = NF_ENDDEF(ncid)

  RETURN
END SUBROUTINE define_wrf_variables_V2
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE write_wrf_input             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_wrf_input(FileHandler,io_form,DateStr,wrfversion,      &
                     sfcinitopt, sfcdtfl,                               &
                     nx_wrf,ny_wrf,nz_wrf,nxlg_wrf,nylg_wrf,            &
                     nzsoil_wrf,bdy_width,fzone_wrf,dx_wrf,dy_wrf,      &
                     zlevels_wrf, zlevels_half,ptop,                    &
                     mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,    &
                     ctrlat_wrf,ctrlon_wrf,lat_wrf,lon_wrf,             &
                     latu_wrf,lonu_wrf,latv_wrf,lonv_wrf,               &
                     lat_ll,lat_ul,lat_ur,lat_lr,                       &
                     lon_ll,lon_ul,lon_ur,lon_lr,                       &
                     msft_wrf,msfu_wrf,msfv_wrf,zs_wrf,dzs_wrf,         &
                     u_wrf,v_wrf,w_wrf,ph_wrf,phb_wrf,pt_wrf,           &
                     pt_init_wrf,p_wrf,pb_wrf,mup_wrf,mub_wrf,mu_wrf,   &
                     qv_wrf,qscalar_wrf,                                &
                     soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xland_wrf,      &
                     xice_wrf,tmn_wrf,shdmax_wrf,shdmin_wrf,snoalb_wrf, &
                     snowh_wrf,canwat_wrf,albbck_wrf,sst_wrf,           &
                     hterain_wrf,tsk_wrf,tslb_wrf,smois_wrf,            &
                     work1d,work2d,work2di,work3d,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Construct and write WRF input file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  01/12/2005
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  USE wrf_metadata             ! WRF constants and metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: FileHandler       ! file handler
  INTEGER, INTENT(IN)  :: io_form           ! file format
                                            ! IO_NET, IO_INT, IO_PHDF5
  CHARACTER(*), INTENT(IN) :: DateStr       ! Current date string

  REAL,    INTENT(IN)  :: wrfversion
  CHARACTER(LEN=*), INTENT(IN) :: sfcinitopt
  CHARACTER(LEN=*), INTENT(IN) :: sfcdtfl

  INTEGER, INTENT(IN)  :: nx_wrf
  INTEGER, INTENT(IN)  :: ny_wrf
  INTEGER, INTENT(IN)  :: nz_wrf
  INTEGER, INTENT(IN)  :: nzsoil_wrf
  INTEGER, INTENT(IN)  :: nxlg_wrf          ! WRF domain size
  INTEGER, INTENT(IN)  :: nylg_wrf
  INTEGER, INTENT(IN)  :: fzone_wrf         ! WRF fakezone, = 1
  REAL,    INTENT(IN)  :: dx_wrf,dy_wrf

  REAL,    INTENT(IN)  :: zlevels_wrf(nz_wrf)
  REAL,    INTENT(IN)  :: zlevels_half(nz_wrf-1)
  INTEGER, INTENT(IN)  :: bdy_width
  REAL,    INTENT(IN)  :: ptop

  INTEGER, INTENT(IN)  :: mapproj_wrf
  REAL,    INTENT(IN)  :: trulat1_wrf, trulat2_wrf
  REAL,    INTENT(IN)  :: trulon_wrf
  REAL,    INTENT(IN)  :: ctrlat_wrf,ctrlon_wrf

  REAL,    INTENT(IN)  :: lat_wrf (nx_wrf,ny_wrf),latu_wrf(nx_wrf,ny_wrf),latv_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: lon_wrf (nx_wrf,ny_wrf),lonu_wrf(nx_wrf,ny_wrf),lonv_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: msft_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: msfu_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: msfv_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: lat_ll(4),lat_ul(4),lat_ur(4),lat_lr(4)
  REAL,    INTENT(IN)  :: lon_ll(4),lon_ul(4),lon_ur(4),lon_lr(4)

  REAL,    INTENT(IN)  :: zs_wrf(nzsoil_wrf)
  REAL,    INTENT(IN)  :: dzs_wrf(nzsoil_wrf)

  REAL,    INTENT(IN)  :: u_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: v_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: w_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: ph_wrf (nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: phb_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: pt_wrf (nx_wrf,ny_wrf,nz_wrf)  ! Pot. Temp. - t0
  REAL,    INTENT(IN)  :: pt_init_wrf (nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: p_wrf (nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: pb_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: mup_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: mub_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: mu_wrf (nx_wrf,ny_wrf)

  REAL,    INTENT(IN)  :: qv_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: qscalar_wrf(nx_wrf,ny_wrf,nz_wrf,nscalar_wrf)

  INTEGER, INTENT(IN)  :: soiltyp_wrf(nx_wrf,ny_wrf)
  INTEGER, INTENT(IN)  :: vegtyp_wrf (nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: vegfrct_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: xland_wrf (nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: xice_wrf  (nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: tmn_wrf   (nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: snowh_wrf (nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: canwat_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: sst_wrf   (nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: shdmin_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: shdmax_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: snoalb_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: hterain_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: tsk_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: albbck_wrf(nx_wrf,ny_wrf)

  REAL,    INTENT(IN)  :: tslb_wrf(nx_wrf,ny_wrf,nzsoil_wrf)
  REAL,    INTENT(IN)  :: smois_wrf(nx_wrf,ny_wrf,nzsoil_wrf)

  REAL,    INTENT(OUT) :: work1d (nz_wrf)       ! temporary arrays
  REAL,    INTENT(OUT) :: work2d (nx_wrf,ny_wrf)
  INTEGER, INTENT(OUT) :: work2di(nx_wrf,ny_wrf)
  REAL,    INTENT(OUT) :: work3d (nx_wrf,ny_wrf,nz_wrf)
                                 ! assume nz_wrf > nzsoil_wrf

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Temporary local variables
!
!-----------------------------------------------------------------------

  TYPE(wrf_var_metadata)    :: var_meta

  REAL,         ALLOCATABLE :: snowc(:,:)

  REAL,         ALLOCATABLE :: dnw(:)
  REAL,         ALLOCATABLE :: rdnw(:)
  REAL,         ALLOCATABLE :: dn(:)
  REAL,         ALLOCATABLE :: rdn(:)
  REAL,         ALLOCATABLE :: fnp(:)
  REAL,         ALLOCATABLE :: fnm(:)

  REAL,         ALLOCATABLE :: tem3dlg1(:,:,:)   ! domain size array
  REAL,         ALLOCATABLE :: tem3dlg2(:,:,:)   ! Mem. Order 'XZY'

  REAL  :: projrot_latlon            ! external function

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER  :: i,j,k,nq
  REAL     :: cfn, cfn1
  REAL     :: projrot_deg
  REAL     :: p_surf
  REAL     :: t00, p00, a
  REAL     :: cof1, cof2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Allocate temporary arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(snowc(nx_wrf,ny_wrf),    STAT = istatus)

  ALLOCATE(dnw (nz_wrf-1),  STAT = istatus)
  ALLOCATE(rdnw(nz_wrf-1),  STAT = istatus)
  ALLOCATE(dn  (nz_wrf-1),  STAT = istatus)
  ALLOCATE(rdn (nz_wrf-1),  STAT = istatus)
  ALLOCATE(fnp (nz_wrf-1),  STAT = istatus)
  ALLOCATE(fnm (nz_wrf-1),  STAT = istatus)

  ALLOCATE(tem3dlg1(nxlg_wrf,nylg_wrf,nz_wrf),    STAT = istatus)
  ALLOCATE(tem3dlg2(nxlg_wrf,nz_wrf,nylg_wrf),    STAT = istatus)

!-----------------------------------------------------------------------
!
! Land use category
!
!-----------------------------------------------------------------------

  var_meta%name        = 'LU_INDEX'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LAND USE CATEGORY'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  work2d(:,:) = FLOAT(vegtyp_wrf(:,:))
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                        &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

!-------------------------------------------------------------------
!
!   Write out atmospheric variables U (in u_wrf), V (in v_wrf)
!   W, PH, PHB, T, MU, MUB, MU0.
!
!-------------------------------------------------------------------
!
  var_meta%name        = 'U'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'x-wind component'
  var_meta%units       = 'm s-1'
  var_meta%stagger     = 'X'
  var_meta%dimName1    = 'west_east_stag'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               u_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                    &
               tem3dlg1,tem3dlg2,nxlg_wrf,nylg_wrf-1,nz_wrf-1)

  var_meta%name        = 'V'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'y-wind component'
  var_meta%units       = 'm s-1'
  var_meta%stagger     = 'Y'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north_stag'
  var_meta%dimName3    = 'bottom_top'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               v_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                    &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf,nz_wrf-1)

  var_meta%name        = 'W'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'z-wind component'
  var_meta%units       = 'm s-1'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top_stag'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               w_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                    &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf)
  !
  !  Write out perturbation geopotential
  !
  var_meta%name        = 'PH'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'perturbation geopotential'
  var_meta%units       = 'm2 s-2'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top_stag'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               ph_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                   &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf)

  var_meta%name        = 'PHB'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'base-state geopotential'
  var_meta%units       = 'm2 s-2'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top_stag'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               phb_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                  &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf)
  !
  !  Write out perturbation potential temperature (theta - t0)
  !
  var_meta%name        = 'T'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'perturbation potential temperature (theta-t0)'
  var_meta%units       = 'K'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               pt_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                   &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)

  var_meta%name        = 'T_INIT'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'initial potential temperature'
  var_meta%units       = 'K'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               pt_init_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,              &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)
  !
  !
  !  Write out air mass
  !
  var_meta%name        = 'MU'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'perturbation dry air mass in column'
  var_meta%units       = 'Pa'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               mup_wrf,nx_wrf,ny_wrf,fzone_wrf,                         &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'MUB'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'base state dry air mass in column'
  var_meta%units       = 'Pa'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               mub_wrf,nx_wrf,ny_wrf,fzone_wrf,                         &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  IF (wrfversion < 3.1) THEN
    var_meta%name        = 'MU0'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'initial dry mass in column'
    var_meta%units       = 'Pa'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
                 mu_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
                 tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    !
    ! Added since WRFV2.2 but removed since WRFV3.1
    !
    work2d(:,:) = 0.0
    var_meta%name        = 'URATX'
    var_meta%description = 'Ratio of U over U10 on mass points'
    var_meta%units       = 'dimensionless'
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
                 work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
                 tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'VRATX'
    var_meta%description = 'Ratio of V over V10 on mass points'
    var_meta%units       = 'dimensionless'
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
                 work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
                 tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'TRATX'
    var_meta%description = 'Ratio of T over TH2 on mass points'
    var_meta%units       = 'dimensionless'
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
                 work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
                 tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  END IF

  var_meta%name        = 'SR'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'fraction of frozen precipitation'
  var_meta%units       = '-'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  work2d(:,:) = 0.0
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

!-----------------------------------------------------------------------
!
!   WRF mass coordinate varaiables (from real.exe)
!
!-----------------------------------------------------------------------
  !
  !  From the full level data, we can get the half levels, reciprocals,
  !  and layer thicknesses.  These are all defined at half level
  !  locations, so one less level.
  !

  DO k=1, nz_wrf-1
     dnw(k) = zlevels_wrf(k+1) - zlevels_wrf(k)
     rdnw(k) = 1./dnw(k)
  END DO

  !
  !  Now the same sort of computations with the half eta levels, even ANOTHER
  !  level less than the one above.
  !
  dn(1)  = 0.0
  rdn(1) = 0.0
  fnp(1) = 0.0
  fnm(1) = 0.0
  DO k=2, nz_wrf-1
    dn(k) = 0.5*(dnw(k)+dnw(k-1))
    rdn(k) = 1./dn(k)
    fnp(k) = .5* dnw(k  )/dn(k)
    fnm(k) = .5* dnw(k-1)/dn(k)
  END DO

  var_meta%name        = 'FNM'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'upper weight for vertical stretching'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               fnm,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'FNP'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'lower weight for vertical stretching'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               fnp,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'RDNW'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'inverse dn values on full (w) levels'
  var_meta%units       = 'Pa'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               rdnw,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'RDN'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'inverse dn values on half (mass) levels'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               rdn,nz_wrf-1,fzone_wrf)

  !  Scads of vertical coefficients.
  var_meta%name        = 'DNW'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'd(eta) values between full (w) levels'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               dnw,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'DN'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'd(eta) values between half (mass) levels'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               dn,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'ZNU'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'eta values on half (mass) levels'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               zlevels_half,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'ZNW'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'eta values on full (w) levels'
  var_meta%units       = ''
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'bottom_top_stag'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               zlevels_wrf,nz_wrf,fzone_wrf)

  work1d(:) = 0.0
  var_meta%name        = 'T_BASE'          ! set to 0.0
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'BASE STATET T IN IDEALIZED CASES'
  var_meta%units       = 'K'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work1d,nz_wrf-1,fzone_wrf)

  cfn  = (.5*dnw(nz_wrf-1)+dn(nz_wrf-1))/dn(nz_wrf-1)
  cfn1 = -.5*dnw(nz_wrf-1)/dn(nz_wrf-1)

  var_meta%name        = 'CFN'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = '0'
  var_meta%description = 'extrapolation constant'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),cfn,istatus)

  var_meta%name        = 'CFN1'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = '0'
  var_meta%description = 'extrapolation constant'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),cfn1,istatus)

  work1d(:) = 0.0
!  var_meta%name        = 'EPSTS'           ! set to 0.0
!  var_meta%fieldType   = WRF_REAL
!  var_meta%MemoryOrder = '0'
!  var_meta%description = 'leapfrog time filter constant'
!  var_meta%units       = ''
!  var_meta%stagger     = ''
!  var_meta%dimName1    = ''
!  var_meta%dimName2    = ''
!  var_meta%dimName3    = ''
!  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
!               work1d,1,fzone_wrf)

  var_meta%name        = 'STEP_NUMBER'     ! set to 0
  var_meta%fieldType   = WRF_integer
  var_meta%MemoryOrder = '0'
  var_meta%description = ''
  var_meta%units       = '-'
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write_int(FileHandler,io_form,var_meta,DateStr(1:19),0,istatus)

  !
  !  Write out Q2, T2, TH2, U10, V10 etc.
  !  They are only for plotting purposes according to "init_domain_rk" in REAL.EXE.
  !  We are not all that worried.
  !

  CALL const_module_initialize ( p00 , t00 , a )

  work2d(:,:) = 0.0

  var_meta%fieldType   = WRF_REAL   ! common to all 2D arrays below - 5
  var_meta%MemoryOrder = 'XY'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  work2d(:,:) = qv_wrf(:,:,1) ! see WRFV2.1 dyn_em/module_initialize_real.f
  var_meta%name        = 'Q2'
  var_meta%description = 'QV at 2 M'
  var_meta%units       = 'kg kg-1'
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  work2d(:,:) = pt_wrf(:,:,1)+300.  ! see WRFV2.1 dyn_em/module_initialize_real.f
  DO j = 1,ny_wrf
    DO i = 1,nx_wrf
      work2d(i,j) = work2d(i,j)*(((p_wrf(i,j,1)+pb_wrf(i,j,1))/base_pres)**(r_d/cp_wrf))
    END DO
  END DO
  var_meta%name        = 'T2'
  var_meta%description = 'TEMP at 2 M'
  var_meta%units       = 'K'
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  work2d(:,:) = pt_wrf(:,:,1)+300.  ! see WRFV2.1 dyn_em/module_initialize_real.f
  var_meta%name        = 'TH2'
  var_meta%description = 'POT TEMP at 2 M'
  var_meta%units       = 'K'
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  DO j = 1,ny_wrf                  ! see WRFV2.1 dyn_em/module_initialize_real.f
    DO i = 1,nx_wrf
      p_surf = p00* EXP (-t00/a + ( (t00/a)**2 -                        &
                                 2.*g_wrf*hterain_wrf(i,j)/a/r_d) **0.5)
      work2d(i,j) = p_surf + p_wrf(i,j,1)
    END DO
  END DO
  var_meta%name        = 'PSFC'
  var_meta%description = 'SFC PRESSURE'
  var_meta%units       = 'Pa'
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  DO j = 1,ny_wrf          ! should do message passing here, it did not
    DO i = 1,nx_wrf-1      ! because U10/V10 is not used in WRF.
      work2d(i,j) = (u_wrf(i,j,1)+u_wrf(i+1,j,1))*0.5
                           ! see WRFV2.1 dyn_em/module_initialize_real.f
    END DO
  END DO
  var_meta%name        = 'U10'
  var_meta%description = 'U at 10 M'
  var_meta%units       = 'm s-1'
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  DO j = 1,ny_wrf-1
    DO i = 1,nx_wrf
      work2d(i,j) = (v_wrf(i,j,1)+v_wrf(i,j+1,1))*0.5
                           ! see WRFV2.1 dyn_em/module_initialize_real.f
    END DO
  END DO
  var_meta%name        = 'V10'
  var_meta%description = 'V at 10 M'
  var_meta%units       = 'm s-1'
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%fieldType   = WRF_REAL   ! common for all 1D variable below - 9
  var_meta%MemoryOrder = '0'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''

  var_meta%name        = 'RDX'
  var_meta%description = 'INVERSE X GRID LENGTH'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  1/dx_wrf,istatus)

  var_meta%name        = 'RDY'
  var_meta%description = 'INVERSE Y GRID LENGTH'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  1/dy_wrf,istatus)

  var_meta%name        = 'DTS'
  var_meta%description = 'SMALL TIMESTEP'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  0.0,istatus)

  var_meta%name        = 'DTSEPS'
  var_meta%description = 'TIME WEIGHT CONSTANT FOR SMALL STEPS'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  0.0,istatus)

  var_meta%name        = 'RESM'
  var_meta%description = 'TIME WEIGHT CONSTANT FOR SMALL STEPS'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
               0.0,istatus)

  var_meta%name        = 'ZETATOP'
  var_meta%description = 'ZETA AT MODEL TOP'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  0.0,istatus)

  !  Scads of vertical coefficients.
  cof1 = (2.*dn(2)+dn(3))/(dn(2)+dn(3))*dnw(1)/dn(2)
  cof2 =     dn(2)       /(dn(2)+dn(3))*dnw(1)/dn(3)

  work1d(1)  = fnp(2) + cof1
  work1d(2)  = fnm(2) - cof1 - cof2
  work1d(3)  = cof2

  var_meta%name        = 'CF1'
  var_meta%description = '2nd order extrapolation constant'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  work1d(1),istatus)

  var_meta%name        = 'CF2'
  var_meta%description = '2nd order extrapolation constant'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  work1d(2),istatus)

  var_meta%name        = 'CF3'
  var_meta%description = '2nd order extrapolation constant'

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                  work1d(3),istatus)
!
!----------------------------------------------------------------------
!
!  Write out water vapor mixing ratio
!
!----------------------------------------------------------------------
!
  var_meta%fieldType   = WRF_REAL       ! common for all mixing ratios - 7
  var_meta%units       = 'kg kg-1'
  var_meta%MemoryOrder = 'XYZ'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'bottom_top'

  IF (P_QV_wrf > 0) THEN
    var_meta%name        = 'QVAPOR'
    var_meta%description = 'Water vapor mixing ratio'

    CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),            &
                 qv_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                 &
                 tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)
  END IF
  !
  !  Write out QCLOUD and QRAIN, QSNOW, QICE, QGRAUP
  !
  !  NOTE: They are not all zeros according to real.exe.
  !        We can initialize them here use those values from ARPS.
  !
  DO nq = 1,nscalar_wrf

     var_meta%name        = TRIM(qnames_wrf(nq))
     var_meta%description = TRIM(qdescp_wrf(nq))

     CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),           &
               qscalar_wrf(:,:,:,nq),nx_wrf,ny_wrf,nz_wrf,fzone_wrf,    &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)

  END DO

  work3d(:,:,:) = 0.0
  IF (P_QT_wrf > 0) THEN
    var_meta%name        = 'QT'
    var_meta%description = 'Total condensate mixing ratio'
    var_meta%units       = 'kg kg-1'

    CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),            &
              work3d,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                    &
              tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)
  END IF

  work3d(:,:,:) = 0.0
  IF (P_QNI_wrf > 0) THEN
    ! QNICE
    var_meta%name        = 'QNICE'
    var_meta%description = 'Ice Number concentration'
    var_meta%units       = '  kg(-1)'

    CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),            &
              work3d,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                    &
              tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)
  END IF

!----------------------------------------------------------------------
!
!  Write out Soil related variables
!
!-----------------------------------------------------------------------

  !
  !  Write out FCX, GCX, DTBC (0.0 as in real.exe)
  !
  work1d = 0.0

  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'C'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''

  var_meta%name        = 'FCX'
  var_meta%description = 'RELAXATION TERM FOR BOUNDARY ZONE'
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),               &
               work1d,bdy_width,fzone_wrf)

  var_meta%name        = 'GCX'
  var_meta%description = '2ND RELAXATION TERM FOR BOUNDARY ZONE'
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),               &
               work1d,bdy_width,fzone_wrf)

  var_meta%name        = 'DTBC'
  var_meta%MemoryOrder = '0'
  var_meta%description = 'TIME SINCE BOUNDARY READ'
  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),0.0,istatus)

  ! LANDMAK  0  water      XLAND     2  water
  !          1  land                 1  land
  !
  DO j = 1, ny_wrf-1
    DO i = 1, nx_wrf-1
      work2d(i,j) = MOD(NINT(xland_wrf(i,j)), 2)
    END DO
  END DO

  var_meta%name        = 'LANDMASK'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LAND MASK (1 FOR LAND, 0 FOR WATER)'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                          &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  !
  !  WRFSI provides these static variables, but ARPS does not.
  !
  var_meta%name        = 'SHDMAX'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'ANNUAL MAX VEG FRACTION'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               shdmax_wrf,nx_wrf,ny_wrf,fzone_wrf,                      &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'SHDMIN'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'ANNUAL MIN VEG FRACTION'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               shdmin_wrf,nx_wrf,ny_wrf,fzone_wrf,                      &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'SNOALB'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'ANNUAL MAX SNOW ALBEDO IN FRACTION'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''
  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               snoalb_wrf,nx_wrf,ny_wrf,fzone_wrf,                      &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  !
  !  Interpolate soil variables
  !
  var_meta%name        = 'TSLB'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'SOIL TEMPERATURE'
  var_meta%units       = 'K'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'soil_layers_stag'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               tslb_wrf,nx_wrf,ny_wrf,nzsoil_wrf,fzone_wrf,             &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nzsoil_wrf)

  !
  !  Soil layers and variables
  !
  var_meta%name        = 'ZS'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'DEPTHS OF CENTERS OF SOIL LAYERS'
  var_meta%units       = 'm'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'soil_layers_stag'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),               &
               zs_wrf,nzsoil_wrf,fzone_wrf)

  var_meta%name        = 'DZS'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'THICKNESSES OF SOIL LAYERS'
  var_meta%units       = 'm'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'soil_layers_stag'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),               &
               dzs_wrf,nzsoil_wrf,fzone_wrf)

  var_meta%name        = 'SMOIS'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'SOIL MOISTURE'
  var_meta%units       = 'm3 m-3'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'soil_layers_stag'

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               smois_wrf,nx_wrf,ny_wrf,nzsoil_wrf,fzone_wrf,            &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nzsoil_wrf)

  !   SH2O  soil liquid water, depends on swxxxxxx etc.
  !         we do not know how swxxxxxx  are defined.
  !
  work3d(:,:,:) = 0.0
  DO k = 1,nzsoil_wrf
    DO j = 1, ny_wrf-1
      DO i = 1, nx_wrf-1
        IF(xland_wrf(i,j) > 1.5) work3d(i,j,k) = 1.0   ! over water
      END DO
    END DO
  END DO

  var_meta%name        = 'SH2O'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XYZ'
  var_meta%description = 'SOIL LIQUID WATER'
  var_meta%units       = 'm3 m-3'
  var_meta%stagger     = 'Z'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = 'soil_layers_stag'
  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work3d,nx_wrf,ny_wrf,nzsoil_wrf,fzone_wrf,               &
               tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nzsoil_wrf)

  IF (wrfversion < 3.1) THEN
    var_meta%name        = 'XICE'
  ELSE
    var_meta%name        = 'SEAICE'
  END IF
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'SEA ICE FLAG'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               xice_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  !
  !  Soil type and vegetation type
  !
  var_meta%name        = 'IVGTYP'
  var_meta%fieldType   = WRF_INTEGER
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'DOMINANT VEGETATION CATEGORY'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''
  CALL write2di(FileHandler,io_form,var_meta,DateStr(1:19),             &
               vegtyp_wrf,nx_wrf,ny_wrf,fzone_wrf,                      &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'ISLTYP'
  var_meta%fieldType   = WRF_INTEGER
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'DOMINANT SOIL CATEGORY'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''
  CALL write2di(FileHandler,io_form,var_meta,DateStr(1:19),             &
               soiltyp_wrf,nx_wrf,ny_wrf,fzone_wrf,                     &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'VEGFRA'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'VEGETATION FRACTION'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               vegfrct_wrf,nx_wrf,ny_wrf,fzone_wrf,                     &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  !
  !  Process Water equivalent of accumulated snow
  !
    DO j = 1, ny_wrf-1
      DO i = 1, nx_wrf-1
        snowc(i,j) = 0
        IF(snowh_wrf(i,j) >= 10.0/100) snowc(i,j) = 1
      END DO
    END DO

!  Convert snow depth (in meter) to water equiv. of accum. snow depth(kg/m**2)
!  (where 1 meter liquid water is set equivqlent to 10 meters snow).
!          100 = ((m snow) / 10 ) * (1000 kg/m**3)

    work2d(:,:) = 100.*snowh_wrf(:,:)

    var_meta%name        = 'SNOW'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'SNOW WATER EQUIVALENT'
    var_meta%units       = 'kg m-2'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'SNOWH'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'PHYSICAL SNOW DEPTH'
    var_meta%units       = 'm'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    ! SNOWH       may be snowdpth of ARPS. Its unit is "m" as in ARPS
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             snowh_wrf,nx_wrf,ny_wrf,fzone_wrf,                         &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'RHOSN'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = ' SNOW DENSITY'
    var_meta%units       = 'kg m-3'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    work2d(:,:) = 0.0
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    !
    !  Canopy water amount (meter, in ARPS, kg/m**2 in WRF)
    !
    var_meta%name        = 'CANWAT'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'CANOPY WATER'
    var_meta%units       = 'kg m-2'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             canwat_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'SST'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'SEA SURFACE TEMPERATURE'
    var_meta%units       = 'K'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             sst_wrf,nx_wrf,ny_wrf,fzone_wrf,                           &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    ! ifndsnowh   depends on flag_snowh to indicate whether SNOWH
    !             is provided. set it to 0.

    var_meta%name        = 'FNDSNOWH'
    var_meta%fieldType   = WRF_INTEGER
    var_meta%MemoryOrder = '0'
    var_meta%description = 'SNOWH_LOGICAL'
    var_meta%units       = '-'
    var_meta%stagger     = ''
    var_meta%dimName1    = ''
    var_meta%dimName2    = ''
    var_meta%dimName3    = ''
    CALL write_int(FileHandler,io_form,var_meta,DateStr(1:19),0,istatus)

    ! ifndsoilw   depends on num_sw_levels_input to denotes the
    !             number of soilwxxx provided. No such variable
    !             provided in WRFSI 1.3.2
    !
    var_meta%name        = 'FNDSOILW'
    var_meta%fieldType   = WRF_INTEGER
    var_meta%MemoryOrder = '0'
    var_meta%description = 'SOILW_LOGICAL'
    var_meta%units       = '-'
    var_meta%stagger     = ''
    var_meta%dimName1    = ''
    var_meta%dimName2    = ''
    var_meta%dimName3    = ''
    CALL write_int(FileHandler,io_form,var_meta,DateStr(1:19),0,istatus)

    !
    !  These variable are intialized as 0.0 in real_em.f
    !  WRFSI also does not provide definitions for those variables
    !
    work2d = 0.0

!    var_meta%name        = 'TOTSWDN'
!    var_meta%fieldType   = WRF_REAL
!    var_meta%MemoryOrder = 'XY'
!    var_meta%description = '-'
!    var_meta%units       = '-'
!    var_meta%stagger     = ''
!    var_meta%dimName1    = 'west_east'
!    var_meta%dimName2    = 'south_north'
!    var_meta%dimName3    = ''
!
!    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
!             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
!             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)
!
!    var_meta%name        = 'TOTLWDN'
!    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
!             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
!             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    IF (wrfversion < 3.1) THEN
      var_meta%name        = 'RSWTOA'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'RLWTOA'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'CZMEAN'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'CFRACL'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'CFRACM'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'CFRACH'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'ACFRST'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      work2di(:,:) = 0
      var_meta%name        = 'NCFRST'
      var_meta%fieldType   = WRF_INTEGER
      CALL write2di(FileHandler,io_form,var_meta,DateStr(1:19),           &
                work2di,nx_wrf,ny_wrf,fzone_wrf,                          &
                tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      var_meta%name        = 'ACFRCV'
      CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
               work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

      work2di(:,:) = 0
      var_meta%name        = 'NCFRCV'
      var_meta%fieldType   = WRF_INTEGER
      CALL write2di(FileHandler,io_form,var_meta,DateStr(1:19),           &
               work2di,nx_wrf,ny_wrf,fzone_wrf,                           &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)
    END IF

!-----------------------------------------------------------------------
!
!  Map projection information
!
!-----------------------------------------------------------------------

!     IF (wrfversion < 3.0) THEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Version before WRFV3.0
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    var_meta%name        = 'MAPFAC_M'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Map scale factor on mass grid'
    var_meta%units       = ''
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msft_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'MAPFAC_U'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Map scale factor on u-grid'
    var_meta%units       = ''
    var_meta%stagger     = 'X'
    var_meta%dimName1    = 'west_east_stag'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msfu_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf,nylg_wrf-1)

    var_meta%name        = 'MAPFAC_V'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Map scale factor on v-grid'
    var_meta%units       = ''
    var_meta%stagger     = 'Y'
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north_stag'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msfv_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf-1,nylg_wrf)

    IF (wrfversion >= 3.0) THEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! For WRFV3.0 to accomodate global WRF capbality. However, our code
! only work with conformal map projection at present
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    var_meta%name        = 'MAPFAC_MX'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Map scale factor on mass grid'
    var_meta%units       = ''
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msft_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'MAPFAC_MY'
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msft_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'MAPFAC_UX'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Map scale factor on u-grid'
    var_meta%units       = ''
    var_meta%stagger     = 'X'
    var_meta%dimName1    = 'west_east_stag'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msfu_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf,nylg_wrf-1)

    var_meta%name        = 'MAPFAC_UY'
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msfu_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf,nylg_wrf-1)

    var_meta%name        = 'MAPFAC_VX'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Map scale factor on v-grid'
    var_meta%units       = ''
    var_meta%stagger     = 'Y'
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north_stag'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msfv_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf-1,nylg_wrf)

    var_meta%name        = 'MAPFAC_VY'
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             msfv_wrf,nx_wrf,ny_wrf,fzone_wrf,                          &
             tem3dlg1,nxlg_wrf-1,nylg_wrf)

    var_meta%name        = 'MF_VX_INV'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Inverse of map scale factor on v-grid'
    var_meta%units       = ''
    var_meta%stagger     = 'Y'
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north_stag'
    var_meta%dimName3    = ''

    DO j=1,ny_wrf
      DO i=1,nx_wrf-1
        IF (msfv_wrf(i,j) /= 0.0) THEN
          work2d(i,j)  = 1.0/ msfv_wrf(i,j)
        ELSE
          work2d(i,j)  = 0.0
        END IF
      END DO
    END DO
    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
             tem3dlg1,nxlg_wrf-1,nylg_wrf)

    var_meta%name        = 'CLAT'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'COMPUTATION GRID LATITUDE, SOUTH IS NEGATIVE'
    var_meta%units       = 'degree_north'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             lat_wrf,nx_wrf,ny_wrf,fzone_wrf,                           &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'CLONG'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'COMPUTATION GRID LONGITUDE, WEST IS NEGATIVE'
    var_meta%units       = 'degree_east'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             lon_wrf,nx_wrf,ny_wrf,fzone_wrf,                           &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    END IF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! END of VERSION 3.0
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DO j=1,ny_wrf-1
      DO i=1,nx_wrf-1
        work3d(i,j,1)=2*omega_ear*SIN(lat_wrf(i,j)*d2rfactor)
        work2d(i,j)  =2*omega_ear*COS(lat_wrf(i,j)*d2rfactor)
      END DO
    END DO

    var_meta%name        = 'F'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Coriolis sine latitude term'
    var_meta%units       = 's-1'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work3d(:,:,1),nx_wrf,ny_wrf,fzone_wrf,                     &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'E'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Coriolis cosine latitude term'
    var_meta%units       = 's-1'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    DO j=1,ny_wrf-1
      DO i=1,nx_wrf-1
        projrot_deg = projrot_latlon(mapproj_wrf,trulat1_wrf,           &
                           trulat2_wrf,trulon_wrf,ctrlat_wrf,ctrlon_wrf,&
                           lat_wrf(i,j),lon_wrf(i,j),istatus)
        IF(istatus == 1)THEN
          work3d(i,j,1)=SIN(d2rfactor*projrot_deg)
          work2d(i,j)  =COS(d2rfactor*projrot_deg)
        ELSE
          CALL wrf_debug(1,'WARNING: Problem in projrot_latlon.')
        END IF
      END DO
    END DO

    var_meta%name        = 'SINALPHA'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Local sine of map rotation'
    var_meta%units       = ''
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work3d(:,:,1),nx_wrf,ny_wrf,fzone_wrf,                     &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'COSALPHA'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Local cosine of map rotation'
    var_meta%units       = ''
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    !
    !  Terrain height
    !
    var_meta%name        = 'HGT'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'Terrain Height'
    var_meta%units       = 'm'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             hterain_wrf,nx_wrf,ny_wrf,fzone_wrf,                       &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

!
! Surface skin temperature TSK. ARPS does not provide such variable.
! use tsoil at the first soil layers to represent.
!
    DO j = 1,ny_wrf-1
      DO i = 1,nx_wrf-1
        IF (xland_wrf(i,j) < 1.5) THEN    ! 1 for land, 2 for water
          work2d(i,j) = tsk_wrf(i,j)
        ELSE
          work2d(i,j) = sst_wrf(i,j)
        END IF
      END DO
    END DO

  var_meta%name        = 'TSK'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'SURFACE SKIN TEMPERATURE'
  var_meta%units       = 'K'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
           work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
           tem3dlg1,nxlg_wrf-1,nylg_wrf-1)


  work1d(:) = 0.0

  var_meta%name        = 'U_BASE'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'Z'
  var_meta%description = 'BASE STATE X WIND IN IDEALIZED CASES'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'bottom_top'
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''

  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work1d,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'V_BASE'
  var_meta%description = 'BASE STATE Y WIND IN IDEALIZED CASES'
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work1d,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'QV_BASE'
  var_meta%description = 'BASE STATE QV IN IDEALIZED CASES'
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work1d,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'Z_BASE'
  var_meta%description = 'BASE STATE HEIGHT IN IDEALIZED CASES'
  CALL write1d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               work1d,nz_wrf-1,fzone_wrf)

  var_meta%name        = 'U_FRAME'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = '0'
  var_meta%description = 'FRAME X WIND'
  var_meta%units       = 'm s-1'
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''
  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),0.0,istatus)

  var_meta%name        = 'V_FRAME'
  var_meta%description = 'FRAME Y WIND'
  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),0.0,istatus)

  !
  ! Vertical mass coordinate parameters
  !
  var_meta%name        = 'P_TOP'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = '0'
  var_meta%description = 'PRESSURE TOP OF THE MODEL'
  var_meta%units       = 'Pa'
  var_meta%stagger     = ''
  var_meta%dimName1    = ''
  var_meta%dimName2    = ''
  var_meta%dimName3    = ''

  CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),ptop,istatus)

!-----------------------------------------------------------------------
!
!  Latitude and Longitude variables
!
!-----------------------------------------------------------------------

  IF (wrfversion < 3.1) THEN
    !
    ! WRFV2.0.3, Latitude & longitude at corners
    !

    var_meta%fieldType   = WRF_REAL
    var_meta%memoryOrder = '0'
    var_meta%units       = 'degrees'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'ext_scalar'
    var_meta%dimName2    = ''
    var_meta%dimName3    = ''

    !
    ! Latitude, T point
    !

    var_meta%name        = 'LAT_LL_T'
    var_meta%description = 'latitude lower left, temp point'

    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                    lat_ll(1),istatus)

    var_meta%name        = 'LAT_UL_T'
    var_meta%description = 'latitude up left, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ul(1),istatus)

    var_meta%name        = 'LAT_UR_T'
    var_meta%description = 'latitude up right, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ur(1),istatus)

    var_meta%name        = 'LAT_LR_T'
    var_meta%description = 'latitude lower right, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_lr(1),istatus)

    !
    ! Latitude, U point
    !
    var_meta%name        = 'LAT_LL_U'
    var_meta%description = 'latitude lower left, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ll(2),istatus)

    var_meta%name        = 'LAT_UL_U'
    var_meta%description = 'latitude up left, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ul(2),istatus)

    var_meta%name        = 'LAT_UR_U'
    var_meta%description = 'latitude up right, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ur(2),istatus)

    var_meta%name        = 'LAT_LR_U'
    var_meta%description = 'latitude lower right, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_lr(2),istatus)

    !
    ! latitude , V point
    !
    var_meta%name        = 'LAT_LL_V'
    var_meta%description = 'latitude lower left, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ll(3),istatus)

    var_meta%name        = 'LAT_UL_V'
    var_meta%description = 'latitude up left, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ul(3),istatus)

    var_meta%name        = 'LAT_UR_V'
    var_meta%description = 'latitude up right, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_ur(3),istatus)

    var_meta%name        = 'LAT_LR_V'
    var_meta%description = 'latitude lower right, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 lat_lr(3),istatus)

    !
    ! latitude, massless point
    !
    var_meta%name        = 'LAT_LL_D'
    var_meta%description = 'latitude lower left, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),lat_ll(4),istatus)

    var_meta%name        = 'LAT_UL_D'
    var_meta%description = 'latitude up left, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),lat_ul(4),istatus)

    var_meta%name        = 'LAT_UR_D'
    var_meta%description = 'latitude up right, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),lat_ur(4),istatus)

    var_meta%name        = 'LAT_LR_D'
    var_meta%description = 'latitude lower right, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),lat_lr(4),istatus)


    !
    ! longitude, T point
    !

    var_meta%name        = 'LON_LL_T'
    var_meta%description = 'longitude lower left, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ll(1),istatus)

    var_meta%name        = 'LON_UL_T'
    var_meta%description = 'longitude up left, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ul(1),istatus)

    var_meta%name        = 'LON_UR_T'
    var_meta%description = 'longitude up right, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ur(1),istatus)

    var_meta%name        = 'LON_LR_T'
    var_meta%description = 'longitude lower right, temp point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_lr(1),istatus)

    !
    ! longitude, U point
    !

    var_meta%name        = 'LON_LL_U'
    var_meta%description = 'longitude lower left, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ll(2),istatus)

    var_meta%name        = 'LON_UL_U'
    var_meta%description = 'longitude up left, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ul(2),istatus)

    var_meta%name        = 'LON_UR_U'
    var_meta%description = 'longitude up right, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ur(2),istatus)

    var_meta%name        = 'LON_LR_U'
    var_meta%description = 'longitude lower right, u point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_lr(2),istatus)

    !
    ! longitude , V point
    !
    var_meta%name        = 'LON_LL_V'
    var_meta%description = 'longitude lower left, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ll(3),istatus)

    var_meta%name        = 'LON_UL_V'
    var_meta%description = 'longitude up left, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ul(3),istatus)

    var_meta%name        = 'LON_UR_V'
    var_meta%description = 'longitude up right, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),           &
                 LON_ur(3),istatus)

    var_meta%name        = 'LON_LR_V'
    var_meta%description = 'longitude lower right, v point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),            &
                 LON_lr(3),istatus)

    !
    ! longitude, massless point
    !
    var_meta%name        = 'LON_LL_D'
    var_meta%description = 'longitude lower left, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),LON_ll(4),istatus)

    var_meta%name        = 'LON_UL_D'
    var_meta%description = 'longitude up left, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),LON_ul(4),istatus)

    var_meta%name        = 'LON_UR_D'
    var_meta%description = 'longitude up right, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),LON_ur(4),istatus)

    var_meta%name        = 'LON_LR_D'
    var_meta%description = 'longitude lower right, massless point'
    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),LON_lr(4),istatus)
  END IF

  !
  !  WRF latitude and longitude
  !
  var_meta%name        = 'XLAT'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LATITUDE, SOUTH IS NEGATIVE'
  var_meta%units       = 'degree_north'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               lat_wrf,nx_wrf,ny_wrf,fzone_wrf,                         &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'XLONG'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LONGITUDE, WEST IS NEGATIVE'
  var_meta%units       = 'degree_east'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               lon_wrf,nx_wrf,ny_wrf,fzone_wrf,                         &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

  var_meta%name        = 'XLAT_U'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LATITUDE, SOUTH IS NEGATIVE'
  var_meta%units       = 'degree_north'
  var_meta%stagger     = 'X'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               latu_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
               tem3dlg1,nxlg_wrf,nylg_wrf-1)

  var_meta%name        = 'XLONG_U'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LONGITUDE, WEST IS NEGATIVE'
  var_meta%units       = 'degree_east'
  var_meta%stagger     = 'X'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               lonu_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
               tem3dlg1,nxlg_wrf,nylg_wrf-1)

  var_meta%name        = 'XLAT_V'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LATITUDE, SOUTH IS NEGATIVE'
  var_meta%units       = 'degree_north'
  var_meta%stagger     = 'Y'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               latv_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
               tem3dlg1,nxlg_wrf-1,nylg_wrf)

  var_meta%name        = 'XLONG_V'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'LONGITUDE, WEST IS NEGATIVE'
  var_meta%units       = 'degree_east'
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               lonv_wrf,nx_wrf,ny_wrf,fzone_wrf,                        &
               tem3dlg1,nxlg_wrf-1,nylg_wrf)

!-----------------------------------------------------------------------
!
! Misc. surface variables
!
!-----------------------------------------------------------------------

  var_meta%name        = 'ALBBCK'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY'
  var_meta%description = 'BACKGROUND ALBEDO'
  var_meta%units       = ''
  var_meta%stagger     = ''
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ''

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),             &
               albbck_wrf,nx_wrf,ny_wrf,fzone_wrf,                     &
               tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

!  IF (tmn_wrf(1,1) >= 0) THEN
!    work2d(:,:) = tmn_wrf(:,:)
!  ELSE
    !
    ! Soil temperature at lower boundary, use tsoil at top layer
    !
    DO j = 1,ny_wrf-1
      DO i = 1,nx_wrf-1
        IF (xland_wrf(i,j) < 1.5) THEN    ! 1 for land, 2 for water
          work2d(i,j) = tsk_wrf(i,j)
        ELSE
          work2d(i,j) = sst_wrf(i,j)
        END IF
      END DO
    END DO
!  END IF

    var_meta%name        = 'TMN'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'SOIL TEMPERATURE AT LOWER BOUNDARY'
    var_meta%units       = 'K'
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             work2d,nx_wrf,ny_wrf,fzone_wrf,                            &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'XLAND'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'LAND MASK (1 FOR LAND, 2 FOR WATER)'
    var_meta%units       = ''
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             xland_wrf,nx_wrf,ny_wrf,fzone_wrf,                         &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

    var_meta%name        = 'SNOWC'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XY'
    var_meta%description = 'FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER)'
    var_meta%units       = ''
    var_meta%stagger     = ''
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = ''

    CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),            &
             snowc,nx_wrf,ny_wrf,fzone_wrf,                             &
             tem3dlg1,nxlg_wrf-1,nylg_wrf-1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! New fields for WRFV3.1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IF (wrfversion >= 3.1) THEN
    var_meta%name        = 'P'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = 'XYZ'
    var_meta%description = 'perturbation pressure'
    var_meta%units       = 'Pa'
    var_meta%stagger     = ' '
    var_meta%dimName1    = 'west_east'
    var_meta%dimName2    = 'south_north'
    var_meta%dimName3    = 'bottom_top'

    CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),            &
                 p_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                  &
                 tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)

    var_meta%name        = 'PB'
    var_meta%description = 'BASE STATE PRESSURE'

    CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),            &
                 pb_wrf,nx_wrf,ny_wrf,nz_wrf,fzone_wrf,                 &
                 tem3dlg1,tem3dlg2,nxlg_wrf-1,nylg_wrf-1,nz_wrf-1)

    var_meta%name        = 'T00'
    var_meta%fieldType   = WRF_REAL
    var_meta%MemoryOrder = '0  '
    var_meta%description = 'BASE STATE TEMPERATURE'
    var_meta%units       = 'K'
    var_meta%stagger     = ' '
    var_meta%dimName1    = ' '
    var_meta%dimName2    = ' '
    var_meta%dimName3    = ' '

    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),         &
                    base_temp,istatus)

    var_meta%name        = 'P00'
    var_meta%description = 'BASE STATE PRESSURE'
    var_meta%units       = 'Pa'

    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),         &
                    base_pres,istatus)

    var_meta%name        = 'TLP'
    var_meta%description = 'BASE STATE LAPSE RATE'
    var_meta%units       = ' '

    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),         &
                    base_lapse,istatus)

    var_meta%name        = 'TISO'
    var_meta%description = 'TEMP AT WHICH THE BASE T TURNS CONST'
    var_meta%units       = 'K'

    CALL write_real(FileHandler,io_form,var_meta,DateStr(1:19),         &
                    iso_temp,istatus)

    var_meta%name        = 'SAVE_TOPO_FROM_REAL'
    var_meta%fieldType   = WRF_INTEGER
    var_meta%description = '(1=original topo from real/0=topo modified by WRF) flag if input has original topography from real'
    var_meta%units       = 'flag'

    CALL write_int(FileHandler,io_form,var_meta,DateStr(1:19),1,istatus)

    IF (sfcinitopt == 'GEOGRID') THEN
      CALL write_new_gravity_drag_static(FileHandler,io_form,DateStr(1:19), &
               nxlg_wrf-1,nylg_wrf-1,nz_wrf,sfcinitopt,sfcdtfl,tem3dlg1,    &
               istatus)
    ELSE
      WRITE(6,'(/,1x,3a,/)') 'WARNING: Surface option "',               &
        TRIM(sfcinitopt),'" does not contain new data for the gravity drag option.'
    END IF

  END IF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! End of WRFV3.1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  DEALLOCATE(snowc,                             STAT = istatus)
  DEALLOCATE(dnw, rdnw, dn, rdn, fnp, fnm,      STAT = istatus)
  DEALLOCATE(tem3dlg1, tem3dlg2,                STAT = istatus)

  RETURN
END SUBROUTINE write_wrf_input
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE write_wrf_namelist          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_wrf_namelist(nchout,filename,io_form,wrfversion,       &
                   max_dom,input_from_file,nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,  &
                   nproc_x,nproc_y,fzone,dx_wrf,dy_wrf,                 &
                   dt,spec_bdy_width,parent_time_step_ratio,            &
                   i_parent_start,j_parent_start,parent_id,parent_grid_ratio,&
                   eyear,emonth,eday,ehour,eminute,esecond,tintv_bdyin, &
                   syear,smonth,sday,shour,sminute,ssecond,             &
                   mp_physics,ra_lw_physics,ra_sw_physics,sf_sfclay_physics, &
                   sf_surface_physics,bl_pbl_physics,cu_physics,        &
                   diff_opt,km_opt,khdif,kvdif,nprocx_wrf,nprocy_wrf,   &
                   frames_per_outfile,restart_interval,radt,cudt,       &
                   ifsnow,w_damping,io_form_history,io_form_restart,    &
                   history_interval,intval2d,dir2d,                     &
                   indir,outdir,staticdir,readyfl,                      &
                   moist_adv_opt,dampcoef,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!     Write out WRF namelist file
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER,      INTENT(OUT) :: nchout
  CHARACTER(*), INTENT(IN)  :: filename
  INTEGER,      INTENT(IN)  :: io_form
  REAL,         INTENT(IN)  :: wrfversion
  INTEGER,      INTENT(IN)  :: max_dom
  INTEGER, INTENT(IN)      :: nz_wrf
  INTEGER, INTENT(IN)      :: nzsoil_wrf
  INTEGER, INTENT(IN)      :: nproc_x, nproc_y, fzone
  LOGICAL, INTENT(IN), DIMENSION(max_dom) :: input_from_file
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: nx_wrf, ny_wrf
  REAL,    INTENT(IN), DIMENSION(max_dom) :: dx_wrf, dy_wrf
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: parent_time_step_ratio
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: i_parent_start, j_parent_start
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: parent_grid_ratio
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: parent_id
  REAL,    INTENT(IN)      :: dt
  INTEGER, INTENT(IN)      :: spec_bdy_width
  INTEGER, INTENT(IN)      :: tintv_bdyin
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: syear, smonth, sday, shour, sminute, ssecond
  INTEGER, INTENT(IN)                     :: eyear, emonth, eday, ehour, eminute, esecond
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: mp_physics
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: ra_lw_physics,ra_sw_physics
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: sf_sfclay_physics, sf_surface_physics
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: bl_pbl_physics
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: cu_physics
  REAL,    INTENT(IN), DIMENSION(max_dom) :: khdif,kvdif
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: frames_per_outfile
  REAL,    INTENT(IN), DIMENSION(max_dom) :: radt,cudt
  INTEGER, INTENT(IN)      :: diff_opt,km_opt
  INTEGER, INTENT(IN)      :: nprocx_wrf,nprocy_wrf
  INTEGER, INTENT(IN)      :: restart_interval
  INTEGER, INTENT(IN)      :: ifsnow, w_damping
  CHARACTER(*), INTENT(IN) :: staticdir,indir,outdir
  INTEGER, INTENT(IN)      :: readyfl
  INTEGER, INTENT(IN)      :: io_form_history, io_form_restart
  INTEGER, INTENT(IN), DIMENSION(max_dom) :: history_interval
  INTEGER,           INTENT(IN) :: intval2d
  CHARACTER(LEN=40), INTENT(IN) :: dir2d
  INTEGER, INTENT(IN)      :: moist_adv_opt
  REAL,    INTENT(IN)      :: dampcoef
  INTEGER, INTENT(OUT)     :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: abstime, abstime1
  INTEGER :: rday, rhour, rmin, rsec

  INTEGER :: io_form_wrf

  INTEGER :: nxlg_wrf, nylg_wrf
  INTEGER :: domid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE ( io_form )
  CASE ( IO_INT )
    io_form_wrf = IO_INT_WRF
  CASE ( IO_PHDF5 )
    io_form_wrf = IO_PHDF5_WRF
  CASE ( IO_NET )
    io_form_wrf = IO_NET_WRF
  CASE DEFAULT
    io_form_wrf = 0           ! should not be here
  END SELECT

  !
  ! get the end date and time
  !
  CALL ctim2abss( eyear,emonth,eday,ehour,eminute,esecond, abstime )
  abstime = abstime + tintv_bdyin
  CALL abss2ctim( abstime, eyear, emonth, eday, ehour, eminute, esecond )

  !
  ! get the initial date and time from saved variables
  !
  CALL ctim2abss( syear(1),smonth(1),sday(1),shour(1),sminute(1),ssecond(1), abstime1 )

  !
  ! Compute the run time in days, hours, minutes and seconds
  !
  ! NOTE: abstime and abstime1 are used as temporary variables below
  !
  abstime = abstime - abstime1      ! run time in seconds
  rday    = abstime/(24*3600)
  abstime1= MOD(abstime, 24*3600)   ! remain seconds < one day
  rhour   = abstime1 / 3600
  abstime1= MOD(abstime1,3600)      ! remain seconds < one hour
  rmin    = abstime1 / 60
  rsec    = MOD(abstime1, 60)

!-----------------------------------------------------------------------
!
! Generate a namelist file for WRF version 2.0.3
!
!-----------------------------------------------------------------------

  CALL getunit( nchout )

  OPEN(nchout,FILE=filename, STATUS='UNKNOWN')

  WRITE(6,'(/a,F4.2,a/)') ' Generating namelist for WRF version ',wrfversion,' ...'

  WRITE(nchout,'(a)') '!'
  WRITE(nchout,'(a)') '! Autogenerated namelist for WRF 2.2 from ARPS2WRF.'
  WRITE(nchout,'(a)') '! The variables are in three categories'
  WRITE(nchout,'(a)') '!'
  WRITE(nchout,'(a)') '!     1. WRF defaults, need to be checked/changed before run WRF'
  WRITE(nchout,'(a)') '!     2. ARPS2WRF determined, may be checked/changed before run WRF'
  WRITE(nchout,'(a)') '!     3. Must NOT changed because wrfinput was created based on them.'
  WRITE(nchout,'(a)') '!        -- Capitalized'
  WRITE(nchout,'(a)') '!'
  WRITE(nchout,'(a)') '! For detailed instructions, '
  WRITE(nchout,'(a)') '! see http://www.mmm.ucar.edu/wrf/users/wrfv2/wrf-namelist.html.'
  WRITE(nchout,'(a)') '!'

  WRITE(nchout,'(a   )')     '&time_control'
  WRITE(nchout,'(a,I8,a)')   '   run_days               = ',rday,','
  WRITE(nchout,'(a,I8,a)')   '   run_hours              = ',rhour,','
  WRITE(nchout,'(a,I8,a)')   '   run_minutes            = ',rmin,','
  WRITE(nchout,'(a,I8,a)')   '   run_seconds            = ',rsec,','
  WRITE(nchout,FMT='(a)',ADVANCE='NO')   '   START_YEAR             = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   syear(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   START_MONTH            = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   smonth(domid), ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   START_DAY              = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   sday(domid),   ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   START_HOUR             = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   shour(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   START_MINUTE           = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   sminute(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   START_SECOND           = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   ssecond(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   end_year               = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   eyear,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   end_month              = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   emonth,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   end_day                = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   eday,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   end_hour               = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   ehour,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   end_minute             = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   eminute,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   end_second             = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   esecond,  ','
  END DO
  WRITE(nchout,'(/,a,I8,a)')   '   interval_seconds       = ',tintv_bdyin,','
  WRITE(nchout,FMT='(a)',ADVANCE='NO')     '   input_from_file        = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(6x,L,a)',ADVANCE='NO')   input_from_file(domid),','
  END DO

  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   history_interval       = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   history_interval(domid),','
  END DO
  WRITE(nchout,'(a)')     '        ! minutes'

  IF (.FALSE.) THEN   ! True for HWT 2008 only
    WRITE(nchout,FMT='(a,I8,a)') '   output_interval_2d     = ',intval2d,','
    WRITE(nchout,FMT='(a,a,a)')  '   output_dir_2d          = ''',TRIM(dir2d),''','
  END IF

  WRITE(nchout,FMT='(a)',ADVANCE='NO')        '   frames_per_outfile     = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   frames_per_outfile(domid),  ','
  END DO
  WRITE(nchout,*)
  WRITE(nchout,'(a   )')     '   restart                =   .FALSE.,'
  WRITE(nchout,'(a,I8,a)')   '   restart_interval       = ',restart_interval,','
  WRITE(nchout,'(a,I8,a)')   '   io_form_history        = ',io_form_history,','
  WRITE(nchout,'(a,I8,a)')   '   io_form_restart        = ',io_form_restart,','
  WRITE(nchout,'(a,I8,a)')   '   IO_FORM_INPUT          = ',io_form_wrf,','

  WRITE(nchout,'(a,I8,a)')   '   IO_FORM_BOUNDARY       = ',io_form_wrf,','

  IF (.FALSE.) THEN   ! True for LEAD project only.
    WRITE(nchout,'(a,a,a)')  '   staticdir              = ''',TRIM(staticdir), ''','
    WRITE(nchout,'(a,a,a)')  '   indir                  = ''',TRIM(indir), ''','
    WRITE(nchout,'(a,a,a)')  '   outdir                 = ''',TRIM(outdir),''','
    WRITE(nchout,'(a,I8,a)') '   readyfl                = ',readyfl,','
  END IF

  WRITE(nchout,'(a   )')     '   debug_level            =        0,'
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&domains'
  WRITE(nchout,'(a,I8,a)')   '   time_step              = ',INT(dt),',        ! seconds'
  WRITE(nchout,'(a   )')     '   time_step_fract_num    =        0,'
  WRITE(nchout,'(a   )')     '   time_step_fract_den    =        1,'
  WRITE(nchout,'(a,I8,a)')   '   max_dom                = ',max_dom,','
  WRITE(nchout,FMT='(a)',ADVANCE='NO')        '   S_WE                   = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   1,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   E_WE                   = '
  DO domid = 1,max_dom
    nxlg_wrf = (nx_wrf(domid)-fzone)*nproc_x+fzone
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   nxlg_wrf,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   S_SN                   = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   1,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   E_SN                   = '
  DO domid = 1,max_dom
    nylg_wrf = (ny_wrf(domid)-fzone)*nproc_y+fzone
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   nylg_wrf,  ','
  END DO
  WRITE(nchout,FMT='(/a)',ADVANCE='NO')       '   S_VERT                 = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   1,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   E_VERT                 = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   nz_wrf,  ','
  END DO

  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   DX                     = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(F8.0,a)',ADVANCE='NO')   dx_wrf(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   DY                     = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(F8.0,a)',ADVANCE='NO')   dy_wrf(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   grid_id                = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   domid,  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')      '   parent_id              = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   parent_id(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   i_parent_start         = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   i_parent_start(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   j_parent_start         = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   j_parent_start(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   parent_grid_ratio      = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   parent_grid_ratio(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   parent_time_step_ratio = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   parent_time_step_ratio(domid),  ','
  END DO
  WRITE(nchout,'(/,a )')     '   feedback               =        1,'
  WRITE(nchout,'(a   )')     '   smooth_option          =        0,'
  WRITE(nchout,'(a,I8,a)')   '   nproc_x                = ',nprocx_wrf,','
  WRITE(nchout,'(a,I8,a)')   '   nproc_y                = ',nprocy_wrf,','
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&physics'
  WRITE(nchout,FMT='(a   )',ADVANCE='NO')     '   MP_PHYSICS             = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   mp_physics(domid),  ','
  END DO
  WRITE(nchout,'(/,a )')     '   mp_zero_out            =        0,'
  WRITE(nchout,'(a   )')     '   mp_zero_out_thresh     =    1.e-8,'
  WRITE(nchout,FMT='(a)',ADVANCE='NO')     '   ra_lw_physics          = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   ra_lw_physics(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   ra_sw_physics          = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   ra_sw_physics(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   radt                   = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(F8.2,a)',ADVANCE='NO')   radt(domid),  ','
  END DO
  WRITE(nchout,'(a)') '        ! minutes'
  WRITE(nchout,FMT='(a)',ADVANCE='NO')       '   sf_sfclay_physics      = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   sf_sfclay_physics(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   SF_SURFACE_PHYSICS     = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   sf_surface_physics(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   bl_pbl_physics         = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   bl_pbl_physics(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')   '   bldt                   = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   0,  ','
  END DO
  WRITE(nchout,'(a   )') '        ! minutes'
  WRITE(nchout,FMT='(a)',ADVANCE='NO')     '   cu_physics             = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   cu_physics(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a )',ADVANCE='NO')     '   cudt                   = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(F8.2,a)',ADVANCE='NO')   cudt(domid),  ','
  END DO
  WRITE(nchout,'(a   )') '        ! minutes'
  WRITE(nchout,'(a   )')     '   isfflx                 =        1,'
  WRITE(nchout,'(a,I8,a)')   '   ifsnow                 = ',ifsnow,','
  WRITE(nchout,'(a   )')     '   icloud                 =        1,'
  WRITE(nchout,'(a   )')     '   surface_input_source   =        1,'
  WRITE(nchout,'(a,I8,a)')   '   NUM_SOIL_LAYERS        = ',nzsoil_wrf,','
if(ra_lw_physics(1) == 3 .OR. ra_sw_physics(1) == 3) then
  WRITE(nchout,'(a     )')   '   LEVSIZ                 =       59,'
  WRITE(nchout,'(a     )')   '   PAERLEV                =       29,'
  WRITE(nchout,'(a     )')   '   CAM_ABS_DIM1           =        4,'
  WRITE(nchout,FMT='(a)',ADVANCE='NO')   '   CAM_ABS_DIM2           = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(I8,a)',ADVANCE='NO')   nz_wrf,  ','
  END DO
  WRITE(nchout,'(/,a   )')   '   CAM_ABS_FREQ_S         =    21600,'
endif
  WRITE(nchout,'(a     )')   '   maxiens                =        1,'
  WRITE(nchout,'(a     )')   '   maxens                 =        3,'
  WRITE(nchout,'(a     )')   '   maxens2                =        3,'
  WRITE(nchout,'(a     )')   '   maxens3                =       16,'
  WRITE(nchout,'(a     )')   '   ensdim                 =      144,'
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&dynamics'
!  WRITE(nchout,'(a   )')     '   DYN_OPT                =        2,'
!  WRITE(nchout,'(a   )')     '   rk_ord                 =        3,'
  WRITE(nchout,'(a,I8,a)')   '   w_damping              = ',w_damping,','
  WRITE(nchout,'(a,I8,a)')   '   diff_opt               = ',diff_opt,','
  WRITE(nchout,'(a,I8,a)')   '   km_opt                 = ',km_opt,  ','
!  WRITE(nchout,'(a     )')   '   damp_opt               =        0,'
!  WRITE(nchout,'(a     )')   '   zdamp                  =    5000.,   5000.,   5000.,'
  WRITE(nchout,'(a,F8.2,a)') '   dampcoef               = ',dampcoef,',     0.2,     0.2,'
  WRITE(nchout,'(a,F8.2,a)') '   base_temp              = ',base_temp,','
  WRITE(nchout,FMT='(a)',ADVANCE='NO')     '   khdif                  = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(F8.2,a)',ADVANCE='NO')   khdif(domid),  ','
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')     '   kvdif                  = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(F8.2,a)',ADVANCE='NO')   kvdif(domid),  ','
  END DO
  WRITE(nchout,*)

!  WRITE(nchout,'(a   )')     '   smdiv                  =      0.1,     0.1,     0.1,'
!  WRITE(nchout,'(a   )')     '   emdiv                  =     0.01,    0.01,    0.01,'
!  WRITE(nchout,'(a   )')     '   epssm                  =      0.1,     0.1,     0.1,'
  WRITE(nchout,FMT='(a)',ADVANCE='NO')     '   non_hydrostatic        = '
  DO domid = 1,max_dom
    WRITE(nchout,FMT='(a9)',ADVANCE='NO')   '.TRUE.,'
  END DO
  WRITE(nchout,*)
  IF (wrfversion < 3.1) THEN
    IF (moist_adv_opt == 1) THEN
      WRITE(nchout,FMT='(a,L,a)')  '   pd_moist               =  ',.TRUE.,','
    ELSE
      WRITE(nchout,FMT='(a,L,a)')  '   pd_moist               =  ',.FALSE.,','
    END IF
  ELSE
    WRITE(nchout,FMT='(a,I8,a)') '   moist_adv_opt          =  ',moist_adv_opt,','
  END IF
!  WRITE(nchout,'(a   )')     '   time_step_sound        =        4,       4,       4,'
!  WRITE(nchout,'(a   )')     '   h_mom_adv_order        =        5,       5,       5,'
!  WRITE(nchout,'(a   )')     '   v_mom_adv_order        =        3,       3,       3,'
!  WRITE(nchout,'(a   )')     '   h_sca_adv_order        =        5,       5,       5,'
!  WRITE(nchout,'(a   )')     '   v_sca_adv_order        =        3,       3,       3,'
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&bdy_control'
  WRITE(nchout,'(a,I8,a)')   '   SPEC_BDY_WIDTH         = ',spec_bdy_width,','
  WRITE(nchout,'(a   )')     '   spec_zone              =        1,'
  WRITE(nchout,'(a,I8,a)')   '   relax_zone             = ',spec_bdy_width-1,','
  WRITE(nchout,FMT='(a)',ADVANCE='NO')      '   specified              =   .TRUE.,'
  DO domid = 2,max_dom
    WRITE(nchout,FMT='(a9)',ADVANCE='NO')   '.FALSE.,'
  END DO
  WRITE(nchout,FMT='(/,a)',ADVANCE='NO')    '   nested                 =  .FALSE.,'
  DO domid = 2,max_dom
    WRITE(nchout,FMT='(a9)',ADVANCE='NO')   '.TRUE.,'
  END DO
  WRITE(nchout,*)
!  WRITE(nchout,'(a   )')     '   periodic_x             =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   symmetric_xs           =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   symmetric_xe           =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   open_xs                =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   open_xe                =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   periodic_y             =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   symmetric_ys           =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   symmetric_ye           =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   open_ys                =  .FALSE., .FALSE., .FALSE.,'
!  WRITE(nchout,'(a   )')     '   open_ye                =  .FALSE., .FALSE., .FALSE.,'
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&namelist_quilt'
  WRITE(nchout,'(a   )')     '   nio_tasks_per_group    =        0,'
  WRITE(nchout,'(a   )')     '   nio_groups             =        1,'
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&fdda'
  WRITE(nchout,'(a  /)')     '/'

  WRITE(nchout,'(a   )')     '&grib2'
  WRITE(nchout,'(a  /)')     '/'

  IF (wrfversion >= 3.0) THEN
    WRITE(nchout,'(a   )')     '&dfi_control'
    WRITE(nchout,'(a  /)')     '/'
  END IF

  CLOSE(nchout)

  CALL retunit (nchout )

  istatus = 0
  RETURN
END SUBROUTINE write_wrf_namelist
