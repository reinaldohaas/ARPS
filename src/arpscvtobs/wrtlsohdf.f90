!########################################################################
!########################################################################
!#########                                                      #########
!#########                 SUBROUTINE wrtlsohdf                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE wrtlsohdf(maxsta,n_obs_b,yrob,monthob,dayob,hrob,minob, &
                     stn2d,obstype2d,lat,lon, &
                     elev,elevunits, &
                     wmowx, &
                     t,tunits, &
                     td,tdunits, &
                     dd, ddunits, &
                     ff, ffunits, &
                     ddg, ddgunits, &
                     ffg, ffgunits, &
                     pmsl, pmslunits,  &
                     cover, &
                     ceil, ceilunits, &
                     lowcld, lowcldunits, &
                     vis, visunits, &
                     filename)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Writes lso data into an HDF file for use by NCL.
!
!----------------------------------------------------------------------- 
!
! AUTHOR:  Eric Kemp, January 2002.
!
!-----------------------------------------------------------------------
!
! Use modules.
!
!-----------------------------------------------------------------------

  USE hdf_constants

!----------------------------------------------------------------------- 
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Declare arguments.
! 
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: maxsta                 ! Size of array
  INTEGER, INTENT(IN) :: n_obs_b                ! Number of observations
  INTEGER, INTENT(IN) :: yrob(maxsta)           ! Year (UTC)
  INTEGER, INTENT(IN) :: monthob(maxsta)        ! Month (UTC)
  INTEGER, INTENT(IN) :: dayob(maxsta)          ! Day (UTC)
  INTEGER, INTENT(IN) :: hrob(maxsta)           ! Hour (UTC)
  INTEGER, INTENT(IN) :: minob(maxsta)          ! Minute (UTC)
  CHARACTER(LEN=1), INTENT(IN) :: stn2d(5,maxsta)    ! Station ID
  CHARACTER(LEN=*), INTENT(IN) :: obstype2d(maxsta)  ! Type of station.
  REAL, INTENT(IN) :: lat(maxsta)               ! Latitude
  REAL, INTENT(IN) :: lon(maxsta)               ! Longitude
  REAL, INTENT(IN) :: elev(maxsta)              ! Elevation.
  CHARACTER(LEN=*), INTENT(IN) :: elevunits     ! Elevation units.
  INTEGER, INTENT(IN) :: wmowx(maxsta)          ! WMO numerical weather
                                                ! code. 
  REAL, INTENT(IN) :: t(maxsta)                 ! Temperature.
  CHARACTER(LEN=*), INTENT(IN) :: tunits        ! Temperature units.
  REAL, INTENT(IN) :: td(maxsta)                ! Dew point.
  CHARACTER(LEN=*), INTENT(IN) :: tdunits       ! Dew point units.
  REAL, INTENT(IN) :: dd(maxsta)                ! Wind direction
  CHARACTER(LEN=*), INTENT(IN) :: ddunits       ! Wind direction units.
  REAL, INTENT(IN) :: ff(maxsta)                ! Wind speed.
  CHARACTER(LEN=*), INTENT(IN) :: ffunits       ! Wind speed units.
  REAL, INTENT(IN) :: ddg(maxsta)               ! Wind gust direction.
  CHARACTER(LEN=*), INTENT(IN) :: ddgunits      ! Wind gust direction 
                                                ! units.
  REAL, INTENT(IN) :: ffg(maxsta)               ! Wind gust speed.
  CHARACTER(LEN=*), INTENT(IN) :: ffgunits      ! Wind gust speed units.
  REAL, INTENT(IN) :: pmsl(maxsta)              ! MSL pressure.
  CHARACTER(LEN=*), INTENT(IN) :: pmslunits     ! MSL pressure units.
  REAL, INTENT(IN) :: cover(maxsta)             ! Cloud cover.
  REAL, INTENT(IN) :: ceil(maxsta)              ! Ceiling height.
  CHARACTER(LEN=*), INTENT(IN) :: ceilunits     ! Ceiling height units.
  REAL, INTENT(IN) :: lowcld(maxsta)            ! Lowest cloud height.
  CHARACTER(LEN=*), INTENT(IN) :: lowcldunits   ! Lowest cloud height
                                                ! units.
  REAL, INTENT(IN) :: vis(maxsta)               ! Visibility.
  CHARACTER(LEN=*), INTENT(IN) :: visunits      ! Visibility (miles).
  CHARACTER(LEN=*), INTENT(IN) :: filename      ! Filename of output
                                                ! HDF file.

!-----------------------------------------------------------------------
!
! HDF parameters and variables
!
!-----------------------------------------------------------------------
  
  INTEGER, PARAMETER :: maxrank = 2
  INTEGER :: dims(maxrank),start(maxrank),edges(maxrank),stride(maxrank)

  INTEGER :: sd_id, sds_id, dim_id

  REAL, PARAMETER :: missing = -9999.
  INTEGER :: status
  CHARACTER(LEN=132) :: string

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Open HDF file.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'-----------------------------------------------------------'
  PRINT *, 'wrtlsohdf: Creating HDF file ', TRIM(filename)
  sd_id = sfstart(filename, DFACC_CREATE)
  IF (sd_id == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
!
! Write n_obs_b as scalar attribute.
!
!-----------------------------------------------------------------------

  status = sfsnatt(sd_id, 'nobs', DFNT_INT32, 1, n_obs_b)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')

!-----------------------------------------------------------------------
!
! Write yrob array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing yrob array.'
  CALL write_sds (sds_id,sd_id, DFNT_INT32, 'yrob', yrob, 1, n_obs_b,   &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Year'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write monthob array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing monthob array.'
  CALL write_sds (sds_id,sd_id, DFNT_INT32, 'monthob', monthob, 1,      &
                  n_obs_b, missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Month'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write dayob array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing dayob array.'
  CALL write_sds (sds_id,sd_id, DFNT_INT32, 'dayob', dayob, 1,          &
                  n_obs_b, missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Day'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write hrob array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing hrob array.'
  CALL write_sds (sds_id,sd_id, DFNT_INT32, 'hrob', hrob, 1,            &
                  n_obs_b, missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Hour'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write minob array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing minob array.'
  CALL write_sds (sds_id,sd_id, DFNT_INT32, 'minob', minob, 1,          &
                  n_obs_b, missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Minute'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write stn2d array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing stn2d array.'

  dims(2) = n_obs_b
  dims(1) = 5
  CALL write_sds_char(sds_id,sd_id, 'stn', stn2d, 2, dims)

  dim_id = sfdimid(sds_id,1)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'maxnumchar_stn')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'StationID'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write obstype2d array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing obstype2d array.'

  dims(2) = n_obs_b
  dims(1) = 8
  CALL write_sds_char(sds_id,sd_id, 'obstype', obstype2d, 2, dims)

  dim_id = sfdimid(sds_id,1)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'maxnumchar_obstype')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'ObTypeID'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write lat array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing lat array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'lat', lat, 1, n_obs_b,   &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Latitude'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Deg N'
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write lon array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing lon array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'lon', lon, 1, n_obs_b,   &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Longitude'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Deg E'
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write elev array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing elev array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'elev', elev, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Elevation'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(elevunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write wmowx array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing wmowx array.'
  CALL write_sds (sds_id,sd_id, DFNT_INT32, 'wmowx', wmowx, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'WMO Weather Code'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write t array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing t array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 't', t, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Temperature'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(tunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write td array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing td array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'td', td, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Dew Point'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(tdunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write u array.
!
!-----------------------------------------------------------------------
!
!  WRITE(6,*)'wrtlsohdf:  Writing u array.'
!  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'u', u, 1, n_obs_b, &
!                  missing)
!
!  dim_id = sfdimid(sds_id,0)
!  status = sfsdmname(dim_id,'nobs')
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = 'U Wind'
!  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = TRIM(uunits)
!  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  status = sfendacc(sds_id)
!
!-----------------------------------------------------------------------
!
! Write v array.
!
!-----------------------------------------------------------------------
!
!  WRITE(6,*)'wrtlsohdf:  Writing v array.'
!  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'v', v, 1, n_obs_b, &
!                  missing)
!
!  dim_id = sfdimid(sds_id,0)
!  status = sfsdmname(dim_id,'nobs')
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = 'V Wind'
!  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = TRIM(vunits)
!  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  status = sfendacc(sds_id)
!
!-----------------------------------------------------------------------
!
! Write dd array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing dd array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'dd', dd, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Wind Direction'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(ddunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write ff array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing ff array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'ff', ff, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Wind Speed'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(ffunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write ug array.
!
!-----------------------------------------------------------------------
!
!  WRITE(6,*)'wrtlsohdf:  Writing ug array.'
!  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'ug', ug, 1, n_obs_b, &
!                  missing)
!
!  dim_id = sfdimid(sds_id,0)
!  status = sfsdmname(dim_id,'nobs')
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = 'U Wind Gust'
!  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = TRIM(uunits)
!  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  status = sfendacc(sds_id)
!
!-----------------------------------------------------------------------
!
! Write vg array.
!
!-----------------------------------------------------------------------
!
!  WRITE(6,*)'wrtlsohdf:  Writing vg array.'
!  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'vg', vg, 1, n_obs_b, &
!                  missing)
!
!  dim_id = sfdimid(sds_id,0)
!  status = sfsdmname(dim_id,'nobs')
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = 'V Wind Gust'
!  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  string = TRIM(vgunits)
!  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
!                   TRIM(string))
!  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
!
!  status = sfendacc(sds_id)
!
!-----------------------------------------------------------------------
!
! Write ddg array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing ddg array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'ddg', ddg, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Wind Gust Direction'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(ddgunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write ffg array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing ffg array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'ffg', ffg, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Wind Gust Speed'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(ffgunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write pmsl array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing pmsl array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'pmsl', pmsl, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'MSL Pressure'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(pmslunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write cover array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing cover array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32,'cover',cover,1,n_obs_b,   &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Cloud Cover'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'tenths'
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write ceil array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing ceil array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32, 'ceil', ceil, 1, n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Ceiling'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(ceilunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write lowcld array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing lowcld array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32,'lowcld',lowcld,1,n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Lowest Cloud Height'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(lowcldunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------
!
! Write vis array.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'wrtlsohdf:  Writing vis array.'
  CALL write_sds (sds_id,sd_id, DFNT_FLOAT32,'vis',vis,1,n_obs_b, &
                  missing)

  dim_id = sfdimid(sds_id,0)
  status = sfsdmname(dim_id,'nobs')
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'Visibility'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = TRIM(visunits)
  status = sfscatt(sds_id,'units',DFNT_CHAR8,LEN_TRIM(string),      &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)

!-----------------------------------------------------------------------  
!
! Close file and exit.
!   
!-----------------------------------------------------------------------
    
  status = sfend(sd_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

  PRINT *, 'wrtlsohdf: Done'
  WRITE(6,*)'-----------------------------------------------------------'

  RETURN
END SUBROUTINE wrtlsohdf

!########################################################################
!########################################################################
!######                                                            ######
!######                   SUBROUTINE HDF_FAIL                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################
  
SUBROUTINE hdf_fail (disposition, status, id, string)
  
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Gracefully exits program if there is an HDF library call error.
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
! Eric Kemp, January 2002
! Updated documentation.
!
!----------------------------------------------------------------------- 
! 
! Force explicit declarations.
! 
!-----------------------------------------------------------------------
 
  IMPLICIT NONE

!----------------------------------------------------------------------- 
! 
! Declare arguments.
! 
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: disposition, status, id
  CHARACTER(LEN=*), INTENT(IN) :: string
  CHARACTER(LEN=5) :: disposition_string(0:1) = (/ 'ERROR', 'FATAL' /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  PRINT *, 'hdf_fail: ', TRIM(disposition_string(disposition)), &
           ': HDF error. status=', status, ', id=', id
  PRINT *, 'hdf_fail: Message: ', TRIM(string)
  IF (disposition > 0) STOP

END SUBROUTINE hdf_fail

!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE WRITE_SDS                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE write_sds (sds_id, sd_id, type, name, var, rank, dims, missing)
  
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Writes numerical SDS data to a HDF file.
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
! Eric Kemp, 31 March 2000
! Corrected checks for file access errors.
!
! Eric Kemp, January 2002
! Updated documentation.
!
!-----------------------------------------------------------------------
!
! Use modules.
! 
!-----------------------------------------------------------------------

  USE hdf_constants

!-----------------------------------------------------------------------
!
! Force explicit declarations.
! 
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Declare arguments.
! 
!-----------------------------------------------------------------------

  INTEGER, INTENT(OUT) :: sds_id
  INTEGER, INTENT(IN) :: sd_id, type, rank
  INTEGER, INTENT(IN), DIMENSION(rank) :: dims
  REAL, INTENT(IN) :: var, missing
  CHARACTER(LEN=*) :: name

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxrank=8
  INTEGER :: status
  INTEGER, DIMENSION(maxrank) :: start, edges, stride
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  WRITE(6,*) 'write_sds: ', TRIM(name), dims(1:rank)

!-----------------------------------------------------------------------
!
! Create SDS
!
!-----------------------------------------------------------------------

  start(:) = 0
  edges(1:rank) = dims(1:rank)
  stride(:) = 1

  sds_id = sfcreate(sd_id, TRIM(name), type, rank, dims)
  IF (sds_id == FAIL) CALL hdf_fail (1, sds_id, sd_id, 'wrt1sds: sfcreate')
 
!-----------------------------------------------------------------------
! 
! Set compression
!
!-----------------------------------------------------------------------
 
  status = sfscompress(sds_id, comp_type, comp_prm(1))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1sds: sfcompr')

!-----------------------------------------------------------------------
!
! Write the data to the SDS.  
!
!-----------------------------------------------------------------------
  
  status = sfwdata(sds_id, start, stride, edges, var)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1sds: sfwdata')

!-----------------------------------------------------------------------
! 
! Set the missing value and exit.
!
!-----------------------------------------------------------------------
  
  status = sfsnatt(sds_id, '_FillValue', DFNT_FLOAT32, 1, missing)
  IF (status == FAIL) CALL hdf_fail (0, status, sds_id, 'sfsnatt')

END SUBROUTINE write_sds

!########################################################################
!########################################################################
!######                                                            ######
!######                SUBROUTINE WRITE_SDS_CHAR                   ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################
  
SUBROUTINE write_sds_char (sds_id, sd_id, name, var, rank, dims)  

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Writes character SDS data to a HDF file.
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
! Eric Kemp, 31 March 2000
! Corrected checks for file access errors.
!
! Eric Kemp, January 2002.
! Updated documentation.
!
!-----------------------------------------------------------------------
!
! Use modules
!
!-----------------------------------------------------------------------
 
  USE hdf_constants

!-----------------------------------------------------------------------
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------
 
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Declare arguments.  Note that the first dimension is the length of
! the strings, and the rank is increased by 1.
!
!-----------------------------------------------------------------------
 
  INTEGER, INTENT(OUT) :: sds_id
  INTEGER, INTENT(IN)  :: sd_id, rank
  INTEGER, INTENT(IN), DIMENSION(rank) :: dims
  CHARACTER(LEN=*)     :: var, name

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------
 
  INTEGER, PARAMETER :: maxrank=8
  INTEGER, DIMENSION(maxrank) :: start, edges, stride
  INTEGER :: status
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  WRITE(6,*) 'write_sds_char: ', TRIM(name), dims(1:rank)

!-----------------------------------------------------------------------
!
! Create SDS
!
!-----------------------------------------------------------------------
  
  start(:) = 0
  edges(1:rank) = dims(1:rank)
  stride(:) = 1

  sds_id = sfcreate(sd_id, TRIM(name), DFNT_CHAR8, rank, dims)
  IF (sds_id == FAIL) CALL hdf_fail (1, sds_id, sd_id, 'wrt1c: sfcreate')
 
!-----------------------------------------------------------------------
! 
! Set compression.
!
!-----------------------------------------------------------------------
 
  !status = sfscompress(sds_id, comp_type, comp_prm(1))
  !  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1c: sfcompr')

!-----------------------------------------------------------------------
!
! Write the data to the SDS.  
!
!-----------------------------------------------------------------------
  
  status = sfwcdata(sds_id, start, stride, edges, var)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1c: sfwcdata')

END SUBROUTINE write_sds_char

