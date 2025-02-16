!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE define_wrf_static             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE define_wrf_static(ncid,nx,ny,land_cat,soil_cat,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Define NetCDF variables for WRF static file.
!
! NOTE:
!   This subroutine only defines the variables to be written in NetCDF
!   define mode. Please make sure to write the actual data in wrfstatic.f90.
!   because it set NF_SET_FILL to NF_NOFILL.
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: nx
  INTEGER, INTENT(IN)  :: ny
  INTEGER, INTENT(IN)  :: land_cat
  INTEGER, INTENT(IN)  :: soil_cat
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER  :: dimunlim_id
  INTEGER  :: dimx_id, dimxs_id
  INTEGER  :: dimy_id, dimys_id
  INTEGER  :: dimm_id
  INTEGER  :: dimstr_id
  INTEGER  :: dimsoil_id, dimland_id

  INTEGER  :: dimlenx, dimleny
  INTEGER  :: dimlenm, dimlenstr
  INTEGER  :: dimlensoil, dimlenland

  INTEGER  :: varid, oldfillmode

  TYPE(wrf_var_metadata) :: var_meta

  CHARACTER(LEN=9)            :: varname

!
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

  istatus = NF_INQ_DIMID(ncid,'land_cat',dimland_id)
  istatus = NF_INQ_DIMLEN(ncid,dimland_id,dimlenland)

  istatus = NF_INQ_DIMID(ncid,'soil_cat',dimsoil_id)
  istatus = NF_INQ_DIMLEN(ncid,dimsoil_id,dimlensoil)

  istatus = NF_INQ_DIMID(ncid,'month',dimm_id)
  istatus = NF_INQ_DIMLEN(ncid,dimm_id,dimlenm)

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
!   Defines Times string
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_VAR(ncid,'Times',NF_CHAR,2,                     &
                         (/dimstr_id, dimunlim_id/),varid)

!-----------------------------------------------------------------------
!
!  WRF latitude and longitude
!
!-----------------------------------------------------------------------

  var_meta%name        = 'XLAT'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Latitude of each gridpoint of the mass grid'
  var_meta%units       = 'degree_north'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                       (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)

  var_meta%name        = 'XLONG'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Longitude of each gridpoint of the mass grid'
  var_meta%units       = 'degree_east'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                          (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)

!-------------------------------------------------------------------
!
!  Surface characteristics
!
!------------------------------------------------------------------

  var_meta%name        = 'LANDMASK'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Flag values indicating water(0) or land (1) points'
  var_meta%units       = 'flag: 1=land pt'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'HGT'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Topographic height on the mass grid points'
  var_meta%units       = 'm'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'TOPOSTDV'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Sub-gridscale standard deviation of topography'
  var_meta%units       = 'm m-1'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'TOPOSLPX'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Sub-gridscale mean topographic slope in x-direction'
  var_meta%units       = 'm m-1'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'TOPOSLPY'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Sub-gridscale mean topographic slope in y-direction'
  var_meta%units       = 'm m-1'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Map projection related variables
!
!-----------------------------------------------------------------------

  var_meta%name        = 'COSALPHA'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Cosine of angle between mass point longitudes and std longitude'
  var_meta%units       = ' '
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'SINALPHA'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Sin of angle between mass point longitudes and std longitude'
  var_meta%units       = ' '
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'F'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Horizontal component of the coriolis acceleration'
  var_meta%units       = 's-1'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'E'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Vertical component of the coriolis acceleration'
  var_meta%units       = 's-1'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'MAPFAC_M'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Map scale factor on mass grid'
  var_meta%units       = ' '
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'MAPFAC_U'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Map scale factor on u grid'
  var_meta%units       = ' '
  var_meta%stagger     = 'X'

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimx_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'MAPFAC_V'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Map scale factor on v grid'
  var_meta%units       = ' '
  var_meta%stagger     = 'Y'

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimy_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Greenness fraction
!
!-----------------------------------------------------------------------

  var_meta%name        = 'SHDMAX'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Maximum annual greenness fraction'
  var_meta%units       = 'Percentage'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'SHDMIN'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Minimum annual greenness fraction'
  var_meta%units       = 'Percentage'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Snow albedo
!
!-----------------------------------------------------------------------

  var_meta%name        = 'SNOALB'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Maximum albedo over snow'
  var_meta%units       = 'Fraction'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Average annual temperature
!
!-----------------------------------------------------------------------

  var_meta%name        = 'TMN'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Average annual temperature'
  var_meta%units       = 'K'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Topographical categorical slope
!
!-----------------------------------------------------------------------

  var_meta%name        = 'SLOPECAT'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XY '
  var_meta%description = 'Topographical categorical slope'
  var_meta%units       = 'Category 0-9'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,3,             &
                           (/dimxs_id,dimys_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Land use
!
!-----------------------------------------------------------------------

  var_meta%name        = 'LANDUSEF'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XYZ'
  var_meta%description = 'Land use categorical fraction on mass grid'
  var_meta%units       = 'Fraction'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,             &
                     (/dimxs_id,dimys_id,dimland_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Soil type
!
!-----------------------------------------------------------------------

  var_meta%name        = 'SOILCTOP'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XYZ'
  var_meta%description = 'Top layer soil type as a categorical fraction'
  var_meta%units       = 'Fraction'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,             &
                     (/dimxs_id,dimys_id,dimsoil_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
  var_meta%name        = 'SOILCBOT'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XYZ'
  var_meta%description = 'Bottom layer soil type as a categorical fraction'
  var_meta%units       = 'Fraction'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,             &
                     (/dimxs_id,dimys_id,dimsoil_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
!  Greenness fraction
!
!-----------------------------------------------------------------------

  var_meta%name        = 'GREEN12M'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XYZ'
  var_meta%description = 'Monthly climatological greenness fraction'
  var_meta%units       = 'Fraction'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,             &
                     (/dimxs_id,dimys_id,dimm_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)

!-----------------------------------------------------------------------
!
! Albedo
!
!-----------------------------------------------------------------------

  var_meta%name        = 'ALBDO12M'
  var_meta%fieldType   = WRF_REAL
  var_meta%memoryOrder = 'XYZ'
  var_meta%description = 'Monthly climatological albedo'
  var_meta%units       = 'Fraction'
  var_meta%stagger     = ' '

  istatus = NF_DEF_VAR(ncid,TRIM(var_meta%name),NF_FLOAT,4,             &
                     (/dimxs_id,dimys_id,dimm_id,dimunlim_id/),varid)
  CALL write_var_meta(ncid,varid,var_meta)
 
!-----------------------------------------------------------------------
!
! Before return
!
!-----------------------------------------------------------------------

  istatus = NF_ENDDEF(ncid)

  RETURN
END SUBROUTINE define_wrf_static
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_wrf_static              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_wrf_static(FileHandler,io_form,DateStr,                &
                     nx,ny,land_cat,soil_cat,                           &
                     lat,lon,f,e,msft,msfu,msfv,cosalpha,sinalpha,      &
                     hgt,topostdv,toposlpx,toposlpy,landmask,landusef,  &
                     soilctop,soilcbot,tmn,slopecat,snoalb,             &
                     shdmax,shdmin,green12m,albdo12m, istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!   Write static arrays
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: FileHandler
  INTEGER, INTENT(IN)  :: io_form
  CHARACTER(LEN=19), INTENT(IN) :: DateStr
  INTEGER, INTENT(IN)  :: nx,ny
  INTEGER, INTENT(IN)  :: land_cat, soil_cat
  REAL,    INTENT(IN)  :: lat(nx,ny), lon(nx,ny)
  REAL,    INTENT(IN)  :: msft(nx,ny), msfu(nx,ny),msfv(nx,ny)
  REAL,    INTENT(IN)  :: f(nx,ny),e(nx,ny)
  REAL,    INTENT(IN)  :: cosalpha(nx,ny),sinalpha(nx,ny)
  REAL,    INTENT(IN)  :: hgt(nx,ny),topostdv(nx,ny)
  REAL,    INTENT(IN)  :: toposlpx(nx,ny),toposlpy(nx,ny)
  REAL,    INTENT(IN)  :: landmask(nx,ny),landusef(nx,ny,land_cat)
  REAL,    INTENT(IN)  :: soilctop(nx,ny,soil_cat),soilcbot(nx,ny,soil_cat)
  REAL,    INTENT(IN)  :: tmn(nx,ny)
  REAL,    INTENT(IN)  :: slopecat(nx,ny)
  REAL,    INTENT(IN)  :: snoalb(nx,ny)
  REAL,    INTENT(IN)  :: shdmax(nx,ny)
  REAL,    INTENT(IN)  :: shdmin(nx,ny)
  REAL,    INTENT(IN)  :: green12m(nx,ny,12)
  REAL,    INTENT(IN)  :: albdo12m(nx,ny,12)
  
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER      :: fzone_wrf = 1
  TYPE(wrf_var_metadata)  :: var_meta

  INTEGER           :: nxlg, nylg
  REAL, ALLOCATABLE :: temlg(:,:,:)

!f@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = nx
  nylg = ny

  ALLOCATE(temlg(nxlg,nylg,MAX(land_cat,soil_cat)), STAT = istatus)


  var_meta%name        = 'XLAT'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Latitude of each gridpoint of the mass grid'
  var_meta%units       = 'deg north'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               lat,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)


  var_meta%name        = 'XLONG'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Longitude of each gridpoint of the mass grid'
  var_meta%units       = 'deg east'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               lon,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'LANDMASK'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Flag valuses indicationg water(0) or land(1) points'
  var_meta%units       = 'flag: 1=land pt'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               landmask,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'HGT'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Topographic height on the mass grid points'
  var_meta%units       = 'm'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               hgt,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'TOPOSTDV'
  var_meta%description = 'Sub-gridscale standard deviation of topography'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               topostdv,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'TOPOSLPX'
  var_meta%description = 'Sub-gridscale mean topographic slope in x-direction'
  var_meta%units       = 'm m-1'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               toposlpx,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'TOPOSLPY'
  var_meta%description = 'Sub-gridscale mean topographic slope in y-direction'
  var_meta%units       = 'm m-1'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               toposlpy,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'COSALPHA'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Cosine of angle between mass point longitudes and std lon'
  var_meta%units       = 'Dimensionless'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               cosalpha,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'SINALPHA'
  var_meta%description = 'Sin of angle between mass point longitude and std longitude'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               sinalpha,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'F'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Horizontal component of the coriolis acceleartion'
  var_meta%units       = 's-1'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               f,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'E'
  var_meta%description = 'Vertical component of the coriolis acceleartion'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               e,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'MAPFAC_M'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Map scale factor on mass grid'
  var_meta%units       = 'Dimensionless'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               msft,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'MAPFAC_U'
  var_meta%description = 'Map scale factor on u grid'
  var_meta%dimName1    = 'west_east_stag'
  var_meta%stagger     = 'X'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               msfu,nx,ny,fzone_wrf,temlg,nxlg,nylg-1)

  var_meta%name        = 'MAPFAC_V'
  var_meta%description = 'Map scale factor on v grid'
  var_meta%stagger     = 'Y'
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north_stag'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               msfv,nx,ny,fzone_wrf,temlg,nxlg-1,nylg)

  var_meta%name        = 'SHDMAX'
  var_meta%fieldType   = WRF_REAL
  var_meta%MemoryOrder = 'XY '
  var_meta%description = 'Maximum annual greenness fraction'
  var_meta%units       = 'Percentage'
  var_meta%stagger     = ' '
  var_meta%dimName1    = 'west_east'
  var_meta%dimName2    = 'south_north'
  var_meta%dimName3    = ' '

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               shdmax,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'SHDMIN'
  var_meta%description = 'Minimun annual greenness fraction'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               shdmin,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'SNOALB'
  var_meta%units       = 'Fraction'
  var_meta%description = 'Maximum albedo over snow'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               snoalb,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'TMN'
  var_meta%units       = 'K'
  var_meta%description = 'Average annual temperature'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               tmn,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'SLOPECAT'
  var_meta%units       = 'Category 0-9'
  var_meta%description = 'Topographical categorical slope'

  CALL write2d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               slopecat,nx,ny,fzone_wrf,temlg,nxlg-1,nylg-1)

  var_meta%name        = 'LANDUSEF'
  var_meta%MemoryOrder = 'XYZ'
  var_meta%units       = 'Fraction'
  var_meta%description = 'Land use categorical fraction on mass grid'
  var_meta%dimName3    = 'land_cat '
  var_meta%stagger     = 'Z'  ! actually no stagger

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               landusef,nx,ny,land_cat,fzone_wrf,                       &
               temlg,temlg,nxlg-1,nylg-1,land_cat)

  var_meta%name        = 'SOILCTOP'
  var_meta%MemoryOrder = 'XYZ'
  var_meta%units       = 'Fraction'
  var_meta%description = 'Top layer soil type as a categorical fraction'
  var_meta%dimName3    = 'soil_cat '
  var_meta%stagger     = 'Z'  ! actually no stagger

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               soilctop,nx,ny,soil_cat,fzone_wrf,                       &
               temlg,temlg,nxlg-1,nylg-1,soil_cat)

  var_meta%name        = 'SOILCBOT'
  var_meta%MemoryOrder = 'XYZ'
  var_meta%units       = 'Percentage'
  var_meta%description = 'Bottom layer soil type as a categorical fraction'
  var_meta%dimName3    = 'soil_cat '
  var_meta%stagger     = 'Z'  ! actually no stagger

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               soilcbot,nx,ny,soil_cat,fzone_wrf,                       &
               temlg,temlg,nxlg-1,nylg-1,soil_cat)

  var_meta%name        = 'GREEN12M'
  var_meta%MemoryOrder = 'XYZ'
  var_meta%units       = 'Fraction'
  var_meta%description = 'Monthly climatological greenness fraction'
  var_meta%dimName3    = 'month '
  var_meta%stagger     = 'Z'  ! actually no stagger

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               green12m,nx,ny,12,fzone_wrf,                             &
               temlg,temlg,nxlg-1,nylg-1,12)

  var_meta%name        = 'ALBDO12M'
  var_meta%MemoryOrder = 'XYZ'
  var_meta%units       = 'Fraction'
  var_meta%description = 'Monthly climatological albedo'
  var_meta%dimName3    = 'month '
  var_meta%stagger     = 'Z'  ! actually no stagger

  CALL write3d(FileHandler,io_form,var_meta,DateStr(1:19),              &
               albdo12m,nx,ny,12,fzone_wrf,                             &
               temlg,temlg,nxlg-1,nylg-1,12)

  DEALLOCATE(temlg)

  istatus = 0

  RETURN
END SUBROUTINE write_wrf_static
