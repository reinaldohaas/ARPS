!
!##################################################################
!##################################################################
!######                                                      ######
!######                     GRIBDIMS.INC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'gribdims.inc' for GRIB dataset decoder
!
!  This file defines the GRIB dimension parameters.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Craig A. Mattocks
!  6/1/1995.
!
!  MODIFICATION HISTORY:
!
!  6/1/1995, version 1.0 (C. Mattocks)
!    Wrote original program.
!
!  12/14/1995 (Yuhe Liu)
!  Modified to use ARPS GRIB definitions.
!
!  October 1999 (Eric Kemp)
!  Modified to process RUC-2 40km isobaric and
!  native coordinate GRIB files (Grid #236).
!
!  Nov. 4, 2000 (Ming Xue)
!  Added T, PT, U, V and RH at sigma=0.995 for AVN
!
!  2003-07-30 (Richard Carpenter)
!  Added GRIB Grid 2 for GFS and NCAR/NCEP Reanalysis
!
!-----------------------------------------------------------------------
!
!  Data dimension should be consistent with the definition in
!  extdims.inc (See extdims.inc).
!
!-----------------------------------------------------------------------
!
!  include 'extdims.inc'
!  integer nxgrb      ! Number of grid points in the x-direction
!  integer nygrb      ! Number of grid points in the y-direction
!  integer nzgrb      ! Number of grid points in the z-direction
!  integer nxygrb     ! Number of grid points in the x*y-directions
!  parameter (nxgrb = nx_ext, nygrb = ny_ext, nzgrb = nz_ext)
!  parameter (nxygrb = nx_ext*ny_ext)
!
!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid 87) RUC model:
!
!-----------------------------------------------------------------------
!
  INTEGER :: ruc87grid              ! Grid identification of RUC
  PARAMETER ( ruc87grid = 87 )   ! RUC Polar Stereo

  INTEGER :: ruc87proj             ! Map projection ID for RUC
  PARAMETER ( ruc87proj = 5 )   ! Polar stereographic projection

  INTEGER :: ruc87nlvt2d      ! Number of RUC level types to be processed
  INTEGER :: ruc87nvs2d       ! Number of RUC 2-D variables to be read
  PARAMETER ( ruc87nlvt2d = 1, ruc87nvs2d = 3 )

  INTEGER :: ruc87levs2d(ruc87nlvt2d) ! Level or Layer identification of RUC
  INTEGER :: ruc87var_id2d(ruc87nvs2d,ruc87nlvt2d)

  DATA ruc87levs2d  / 001/       ! Ground surface
  DATA ruc87var_id2d/ 11,  8, 223/
!  11 = Temperature (K)
!   8 = Geometric height (m)
! 223 = Plant canopy surface water (kg/m**2)

  INTEGER :: ruc87nvs3d    ! Number of RUC 3-D variables to be read
  INTEGER :: ruc87nlvt3d   ! Number of RUC 3-D level types to be read
  PARAMETER ( ruc87nlvt3d = 2, ruc87nvs3d = 6 )

  INTEGER :: ruc87levs3d(ruc87nlvt3d)  ! Number of RUC 3-D variables to be read
  INTEGER :: ruc87var_id3d(ruc87nvs3d, ruc87nlvt3d)

  DATA ruc87levs3d/ 109, 111 /
  DATA ruc87var_id3d/   1,  37,  13, 159,  33,  34,                     &
                     85,  86, 255, 255, 255, 255/
           !   1 = pressure (Pa)
!  37 = montgomery stream function (m**2/s**2)
!  13 = potential temperature (K)
! 159 = condensation pressure of lifted parcel (Pa)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!  85 = soil temperature (K)
!  86 = officially soil moisture content (kg/m**2),
!       but arrays contain volumetric (fractional) moisture,
!       i.e. should have used #144
! 255 = undefined
!
!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid 211) RUC model (awips grid):
!
!-----------------------------------------------------------------------
!
  INTEGER :: ruc211grid  ! Grid identification of ruc211
  PARAMETER ( ruc211grid = 211 )   ! RUC Lambert Conformal

  INTEGER :: ruc211proj  ! Map projection ID for ruc211
  PARAMETER ( ruc211proj = 3 )   ! for Lambert Conformal

  INTEGER :: ruc211nlvt2d      ! Number of ruc211 level types to be processed
  INTEGER :: ruc211nvs2d       ! Number of ruc211 2-D variables to be read
  PARAMETER ( ruc211nvs2d = 4, ruc211nlvt2d = 2 )

  INTEGER :: ruc211levs2d(ruc211nlvt2d) ! Level or Layer identification of RUC
  INTEGER :: ruc211var_id2d(ruc211nvs2d,ruc211nlvt2d)

  DATA ruc211levs2d  / 001, 105/       ! Ground surface
  DATA ruc211var_id2d/ 1,  7, 255, 255,                                 &
                      11, 52,  33,  34/
! 1 = Pressure (Pa)
! 7 = Geometric height (m)
! 11= T_2m (C)
! 52= RH_2m
! 33= U_10m
! 34= V_10m

  INTEGER :: ruc211nvs3d       ! Number of RUC 3-D variables to be read
  INTEGER :: ruc211nlvt3d      ! Number of RUC 3-D level types to be read
  PARAMETER ( ruc211nvs3d = 5, ruc211nlvt3d= 1 )

  INTEGER :: ruc211levs3d(ruc211nlvt3d)       ! Number of RUC 3-D
                                              ! variables to be read
  INTEGER :: ruc211var_id3d(ruc211nvs3d, ruc211nlvt3d)

  DATA ruc211levs3d/ 100 /
  DATA ruc211var_id3d/   7,  11,  52,  33,  34/
!   7 = Pressure surface height (m)
!  11 = Temperature (K)
!  52 = Relative humidity (%)
!  33 = u wind (m/s) (grid relative)
!  34 = v wind (m/s) (grid relative)
!
!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid #212) ETA model (AWIP3D):
!
!-----------------------------------------------------------------------
!
  INTEGER :: eta212grid              ! Grid identification of ETA
  PARAMETER ( eta212grid = 212 )     ! ETA 40 km Lambert Conformal projection

  INTEGER :: eta212proj             ! Map projection ID for RUC
  PARAMETER ( eta212proj = 3 )      ! Lambert Confornal projection

  INTEGER :: eta212nlvt2d  ! Number of ETA 2-D level types to be processed
  INTEGER :: eta212nvs2d   ! Number of ETA 2-D variables   to be processed
  PARAMETER ( eta212nlvt2d = 2, eta212nvs2d = 8 ) ! F.KONG

  INTEGER :: eta212levs2d(eta212nlvt2d) ! 2-D Level identification of ETA
  INTEGER :: eta212var_id2d(eta212nvs2d,eta212nlvt2d)

  DATA eta212levs2d  / 001, 105/   ! At ground surface and meters above
  DATA eta212var_id2d/   1,   7,  11, 223,  81,  65,  80,  61,      &
                        11,  51,  33,  34, 255, 255, 255, 255 /
 !   1 = Surface pressure (Pa)
 !   7 = Geopotential height (gpm)
 !  11 = Surface temperature (K)
 ! 223 = Plant canopy surface water (kg/m**2)
 !  81 = Land/sea mask
 !  65 = Water equiv. of accum snow depth (kg/m**2)
 !  80 = Water temperature (K)           ! EMK
 !  61 = Total precipitation (kg/m**2)
!! F.KONG add second row for:
 !  11 = T_2m (K)
 !  51 = Qv_2m (kg/kg)
 !  33 = U_10m (m/s)
 !  34 = V_10m (m/s)

  INTEGER :: eta212nlvt3d        ! Number of ETA 3-D level types to be processed
  INTEGER :: eta212nvs3d         ! Number of ETA 3-D variables to be processed
  PARAMETER ( eta212nlvt3d = 2, eta212nvs3d = 6 )

  INTEGER :: eta212levs3d(eta212nlvt3d)   ! 3-D Level identification of ETA
  INTEGER :: eta212var_id3d     (eta212nvs3d,eta212nlvt3d)
  INTEGER :: eta212var_id3d_sref(eta212nvs3d,eta212nlvt3d)

  DATA eta212levs3d / 100, 112 /          ! sigma vertical coordinates
  DATA eta212var_id3d/       11,  51,  33,  34,   7,  39,               &
                             85, 144, 255, 255, 255, 255/

  DATA eta212var_id3d_sref/  11,  52,  33,  34,   7,  39,               &
                             85, 144, 255, 255, 255, 255/
!  11 = Temperature (K)
!  51 = Specific humidity (kg/kg)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!   7 = Geopotential height (gpm)
!  39 = Pressure vertical velocity (Pa/s)
!  85 = Soil temperature (K)
! 144 = Volumetric soil moisture (m**3/m**3)
! 255 = undefined
!

!-----------------------------------------------------------------------
!
!  Parameters for NCEP GRIB (Grid #4) GFS model (AWIP3D):
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: gfs04grid = 4   ! Global 0.5 degree GFS data

  INTEGER, PARAMETER :: gfs04proj = 0   ! lat/lon

  ! Number of 2-D level types and Number of 2-D variables to be processed
  INTEGER, PARAMETER :: gfs04nlvt2d = 3, gfs04nvs2d = 6

  INTEGER :: gfs04levs2d  (gfs04nlvt2d) ! 2-D Level identification of ETA
  INTEGER :: gfs04var_id2d(gfs04nvs2d,gfs04nlvt2d)

  DATA gfs04levs2d  / 001, 105, 200 /   ! At ground surface and meters above
  DATA gfs04var_id2d/   1,   7,   11,  65,  81,  91,         &
                        11,  51,  33,  34, 255, 255,       &
                        54, 255, 255, 255, 255, 255 /
 !   1 = Surface pressure (Pa)
 !   7 = Geopotential height (gpm)
 !  11 = Surface temperature (K)
 !  65 = Water equiv. of accum snow depth (kg/m**2)
 !  81 = Land cover(land=1, sea=0)
 !  91 = Ice cover (ice=1,no ice=0)
! second row for:
 !  11 = T_2m (K)
 !  51 = Qv_2m (kg/kg)
 !  33 = U_10m (m/s)
 !  34 = V_10m (m/s)

  ! Number of 3-D level types and Number of 3-D variables to be processed
  INTEGER, PARAMETER :: gfs04nlvt3d = 2, gfs04nvs3d = 6

  INTEGER :: gfs04levs3d  (gfs04nlvt3d)   ! 3-D Level identification
  INTEGER :: gfs04var_id3d(gfs04nvs3d,gfs04nlvt3d)

  DATA gfs04levs3d / 100, 112 /          ! sigma vertical coordinates
  DATA gfs04var_id3d/  7, 11,   33,  34,  52, 153,                &
                      11, 144, 255, 255, 255, 255 /
!   7 = Geopotential height (gpm)
!  11 = Temperature (K)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!  39 = Pressure vertical velocity (Pa/s)
!  52 = Relative humidity (kg/kg)
!  11 = Soil temperature (K)
! 144 = Volumetric soil moisture (m**3/m**3)
! 153 = Cloud water (kg/kg)
! 255 = undefined
!

!-----------------------------------------------------------------------
!TINA
!  Parameters for NARR reanalysis (Grid #221) NAM model (AWIP3D):
!
!  for North American Regional Reanalysis Data
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: narr221grid = 221    ! Grid identification of NAM
                                             ! 32 km NARR reanalysis

  INTEGER, PARAMETER :: narr221proj = 3      ! Map projection ID
                                             ! Lambert Conformal projection

  INTEGER, PARAMETER :: narr221nlvt2d = 2 ! Number of 2-D level types to be processed
  INTEGER, PARAMETER :: narr221nvs2d  = 7 ! Number of 2-D variables   to be processed

  INTEGER :: narr221levs2d(narr221nlvt2d) ! 2-D Level identification
  INTEGER :: narr221var_id2d(narr221nvs2d,narr221nlvt2d)  ! variable IDs

  DATA narr221levs2d  / 001, 105/   ! At ground surface and meters above
  DATA narr221var_id2d/  1,   7,  11, 223,  81,  65,  61,         &
                        11,  51,  33,  34, 255, 255, 255/


 !TINA - these are defined in nceptab_131
 !   1 = Surface pressure (Pa)
 !   7 = Geopotential height (gpm)   ! defined in NARR fixed field file
 !  11 = Surface temperature (K)
 ! 223 = Plant canopy surface water (kg/m**2)
 !  80 = Water temperature (K)
 !  81 = Land/sea mask               ! defined in NARR fixed field file
 !  65 = Water equiv. of accum snow depth (kg/m**2)
 !  61 = Total precipitation (mm)
 !
 !!!!!!!!!!!!!!!!!!!!!! Second row for:
 !  11 = T_2m (K)
 !  51 = Qv_2m (kg/kg)
 !  33 = U_10m (m/s)
 !  34 = V_10m (m/s)

  INTEGER, PARAMETER :: narr221nlvt3d = 2  ! Number of 3-D level types to be processed
  INTEGER, PARAMETER :: narr221nvs3d  = 8  ! Number of 3-D variables to be processed

  INTEGER :: narr221levs3d(narr221nlvt3d)  ! 3-D Level identification of ETA
  INTEGER :: narr221var_id3d(narr221nvs3d,narr221nlvt3d)

  DATA narr221levs3d / 100, 112 /          ! Pressure vertical coordinates
  DATA narr221var_id3d/ 11,  51,  33,  34,   7,  39, 153, 178,          &
                        85, 144, 255, 255, 255, 255, 255, 255/


!!!!!!!!!!!!!!! 100 = isobaric level !!!!!!!!!!!!!!!!!!
!  11 = Temperature (K)
!  51 = Specific humidity (kg/kg)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!   7 = Geopotential height (gpm)
!  39 = Pressure vertical velocity (Pa/s)
! 153 = Cloud water moisture (kg/kg)      !NEW added 5/3/05
! 178 = Ice mixing ratio (kg/kg)          !NEW added 5/3/05
!!!!!!! 112 = layer between two depths below land surface !!!!!!!!!!
!  85 = Soil temperature (K)
! 144 = Volumetric soil moisture (m**3/m**3)
! 255 = undefined
!

!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid #218) ETA model (AWIP3D):
!
!-----------------------------------------------------------------------
!
  INTEGER :: eta218grid              ! Grid identification of ETA
  PARAMETER ( eta218grid = 218 )     ! ETA 12 km Lambert Conformal projection
                                     ! nostandard grid number

  INTEGER :: eta218proj             ! Map projection ID for RUC
  PARAMETER ( eta218proj = 3 )      ! Lambert Confornal projection

  INTEGER :: eta218nlvt2d  ! Number of ETA 2-D level types to be processed
  INTEGER :: eta218nvs2d   ! Number of ETA 2-D variables   to be processed
  PARAMETER ( eta218nlvt2d = 2, eta218nvs2d = 7 )

  INTEGER :: eta218levs2d(eta218nlvt2d) ! 2-D Level identification of ETA
  INTEGER :: eta218var_id2d(eta218nvs2d,eta218nlvt2d)

  DATA eta218levs2d  / 001, 105/   ! At ground surface and meters above
  DATA eta218var_id2d/   1,   7,  11, 223,  81,  65, 61,         &
                        11,  52,  33,  34, 255, 255, 255  /
 !   1 = Surface pressure (Pa)
 !   7 = Geopotential height (gpm)
 !  11 = Surface temperature (K)
 ! 223 = Plant canopy surface water (kg/m**2)
 !  81 = Land/sea mask
 !  65 = Water equiv. of accum snow depth (kg/m**2)
 !  61 = Total precipitation (mm)
 !!!!!!!!!!!!!!!!!!!!! Second row for:
 !  11 = T_2m (K)
 !  52 = Rh_2m (%)
 !  33 = U_10m (m/s)
 !  34 = V_10m (m/s)

  INTEGER :: eta218nlvt3d        ! Number of ETA 3-D level types to be processed
  INTEGER :: eta218nvs3d         ! Number of ETA 3-D variables to be processed
  PARAMETER ( eta218nlvt3d = 2, eta218nvs3d = 6 )

  INTEGER :: eta218levs3d(eta218nlvt3d)   ! 3-D Level identification of ETA
  INTEGER :: eta218var_id3d(eta218nvs3d,eta218nlvt3d)

  DATA eta218levs3d / 100, 112 /          ! Pressure vertical coordinates
  DATA eta218var_id3d/  11,  52,  33,  34,   7,  39,                    &
                        85, 144, 255, 255, 255, 255/
!  11 = Temperature (K)
!  52 = Relative humidity (%)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!   7 = Geopotential height (gpm)
!  39 = Pressure vertical velocity (Pa/s)
!  85 = Soil temperature (K)
! 144 = Volumetric soil moisture (m**3/m**3)
! 255 = undefined

!
!-----------------------------------------------------------------------
!
!  Parameters for Global Re-Analysis on T62 Gaussian grid:
!
!-----------------------------------------------------------------------
!
  INTEGER :: reanalt62grid  ! Grid identification of Global Reanalysis
  PARAMETER ( reanalt62grid = 255 )   ! Actually it uses grid 98
                                      ! but 255 shows up in PDS

  INTEGER :: reanalt62proj  ! Map projection ID for Re-Analaysis
  PARAMETER ( reanalt62proj = 4 )   ! for Gaussian latlon grid
  INTEGER :: reanalt62nlvt2d      ! Number of RUCawips level types to be processed
  INTEGER :: reanalt62nvs2d       ! Number of RUCawips 2-D variables to be read
  PARAMETER ( reanalt62nvs2d = 2, reanalt62nlvt2d = 1 )

  INTEGER :: reanalt62levs2d(reanalt62nlvt2d) ! Level or Layer identification of RUC
  INTEGER :: reanalt62var_id2d(reanalt62nvs2d,reanalt62nlvt2d)

  DATA reanalt62levs2d  / 1 /       ! Ground surface
  DATA reanalt62var_id2d/ 1,  7/
!   1 = Pressure (Pa)
!   7 = Geometric height (m)

  INTEGER :: reanalt62nvs3d       ! Number of 3-D variables to be read
  INTEGER :: reanalt62nlvt3d      ! Number of 3-D level types to be read
  PARAMETER ( reanalt62nvs3d = 4, reanalt62nlvt3d = 1 )

  INTEGER :: reanalt62levs3d(reanalt62nlvt3d)   ! Number of 3-D variables to be read
  INTEGER :: reanalt62var_id3d(reanalt62nvs3d, reanalt62nlvt3d)

  DATA reanalt62levs3d/ 107 /
  DATA reanalt62var_id3d/   11,  51,  33,  34/
 !  11 = Temperature (K)
 !  52 = Specific humidity (%)
 !  33 = u wind (m/s) (grid relative)
 !  34 = v wind (m/s) (grid relative)
 !

!-----------------------------------------------------------------------
!
!  Parameters for NCEP GRIB ( GRID #3) AVN global data
!
!-----------------------------------------------------------------------
!
  INTEGER :: avn3grid             ! Grid identification of AVN
  PARAMETER ( avn3grid = 3 )      ! AVN 1 x 1 degree global lat/lon

  INTEGER :: avn3proj             ! Map projection ID
  PARAMETER ( avn3proj = 0 )      ! lat/lon

  INTEGER :: avn3nlvt2d           ! Number of AVN 2-D level types to process
  INTEGER :: avn3nvs2d            ! Number of AVN 2-D variables   to process
  PARAMETER ( avn3nlvt2d = 4, avn3nvs2d = 6 )

  INTEGER :: avn3levs2d(avn3nlvt2d) ! 2-D Level identification of AVN
  INTEGER :: avn3var_id2d(avn3nvs2d,avn3nlvt2d)

  DATA avn3levs2d   / 001, 107, 102, 105/       ! At ground surface
  DATA avn3var_id2d /   1,   7,  11, 223,  81,   65,   &  ! values at surface
                       11,  33,  34,  52,  13,  255,   &  ! at sigma = 0.995
                        2,   2,   2,   2,   2,    2,   &  ! mean sea level
                       11,  51,  33,  34, 255,  255 /     ! meters above

  !   1 = Surface pressure (Pa)
  !   7 = Geopotential height (gpm)
  !  11 = Surface temperature (K)
  ! 223 = Plant canopy surface water (kg/m**2)
  !  81 = Land/sea mask
  !  65 = Water equiv. accum. snow depth (kg/m**2)

  !  11 = temperature at sigma = 0.995 (K)
  !  33 = u, east-west velocity at sigma = 0.995 (m/s)
  !  34 = v, north-south velocity at sigma = 0.995 (m/s)
  !  52 = relative humidity at sigma = 0.995
  !  13 = potential temperature at sigma = 0.995 (K)
  ! 255 = undefined, unused

  !  2  = mean sea level pressure (Pa)
  ! 255 = undefined, unused

  ! F.KONG add the fourth row for:
  !  11 = T_2m (K)
  !  51 = Qv_2m (kg/kg)
  !  33 = U_10m (m/s)
  !  34 = V_10m (m/s)

  INTEGER :: avn3nlvt3d        ! Number of AVN 3-D level types to be processed
  INTEGER :: avn3nvs3d         ! Number of AVN 3-D variables to be processed
  PARAMETER ( avn3nlvt3d = 2, avn3nvs3d = 7 )

  INTEGER :: avn3levs3d(avn3nlvt3d)   ! 3-D Level identification of AVN
  INTEGER :: avn3var_id3d(avn3nvs3d,avn3nlvt3d)

  DATA avn3levs3d   / 100, 112 /     ! pressure vertical coordinates
  DATA avn3var_id3d /  7,  33,  34,  11,  52,  39, 153,                  &
                      11, 144, 255, 255, 255, 255, 255/
!   7 = Geopotential height (gpm)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!  11 = Temperature (K)
!  52 = Relative humidity
!  39 = Pressure vertical velocity (Pa/s)
! 144 = Volumetric soil moisture (m**3/m**3)
! 153 = cloud water mixing ratio (kg/kg)
! 255 = undefined

!-----------------------------------------------------------------------
!
!  Parameters for NCEP GRIB ( GRID #2) AVN global data
!
!-----------------------------------------------------------------------

  INTEGER :: avn2grid          ! Grid identification of AVN
  PARAMETER ( avn2grid = 2 )   ! AVN 2.5 x 2.5 degree global lat/lon

  INTEGER :: avn2proj             ! Map projection ID
  PARAMETER ( avn2proj = 0 )   ! lat/lon

  INTEGER :: avn2nlvt2d  ! Number of AVN 2-D level types to be processed
  INTEGER :: avn2nvs2d   ! Number of AVN 2-D variables   to be processed
  PARAMETER ( avn2nlvt2d = 3, avn2nvs2d = 6 )

  INTEGER :: avn2levs2d(avn2nlvt2d) ! 2-D Level identification of AVN
  INTEGER :: avn2var_id2d(avn2nvs2d,avn2nlvt2d)

  DATA avn2levs2d  / 001, 107, 102/       ! At ground surface
  DATA avn2var_id2d/   1,   7,  11, 223,  81,   65,   &  ! values at surface
                       11, 33,  34,  52,  13,  255,   &  ! at sigma = 0.995
                       2,  2,   2,   2,   2,   2      /  ! mean sea level

  !   1 = Surface pressure (Pa)
  !   7 = Geopotential height (gpm)
  !  11 = Surface temperature (K)
  ! 223 = Plant canopy surface water (kg/m**2)
  !  81 = Land/sea mask
  !  65 = Water equiv. accum. snow depth (kg/m**2)

  !  11 = temperature at sigma = 0.995 (K)
  !  33 = u, east-west velocity at sigma = 0.995 (m/s)
  !  34 = v, north-south velocity at sigma = 0.995 (m/s)
  !  52 = relative humidity at sigma = 0.995
  !  13 = potential temperature at sigma = 0.995 (K)
  ! 255 = undefined, unused

  !  2  = mean sea level pressure (Pa)
  ! 255 = undefined, unused

  INTEGER :: avn2nlvt3d        ! Number of AVN 3-D level types to be processed
  INTEGER :: avn2nvs3d         ! Number of AVN 3-D variables to be processed
  PARAMETER ( avn2nlvt3d = 2, avn2nvs3d = 6 )

  INTEGER :: avn2levs3d(avn2nlvt3d)   ! 3-D Level identification of AVN
  INTEGER :: avn2var_id3d(avn2nvs3d,avn2nlvt3d)

  DATA avn2levs3d / 100, 112 /     ! pressure vertical coordinates
  DATA avn2var_id3d/  7,  33,  34,  11,  52,  39,                       &
                     11, 144, 255, 255, 255, 255/
!   7 = Geopotential height (gpm)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!  11 = Temperature (K)
!  52 = Relative humidity (kg/kg)
!  39 = Pressure vertical velocity (Pa/s)
! 144 = Volumetric soil moisture (m**3/m**3)
! 255 = undefined

!-----------------------------------------------------------------------
!
!  Parameters for NCEP RUC Native Coordinate GRIB (Grid 236 or Grid 252):
!
!-----------------------------------------------------------------------
!
  INTEGER :: rucn236grid              ! Grid identification of RUC
  PARAMETER ( rucn236grid = 236 )     ! RUC Lambert conformal
  INTEGER :: rucn252grid              ! Grid identification of RUC
  PARAMETER ( rucn252grid = 252 )     ! RUC Lambert conformal
  !  New parameters for 20-km NCEP RUC Native Coordinate GRIB (Grid 252):
  !  J.Case, ENSCO Inc. (8/01/2005)

  INTEGER :: rucn236proj             ! Map projection ID for RUC
  PARAMETER ( rucn236proj = 3 )      ! Lambert conformal projection

  INTEGER :: rucn236nlvt2d      ! Number of RUC level types to be
                                ! processed
  INTEGER :: rucn236nvs2d       ! Number of RUC 2-D variables to be read
  PARAMETER ( rucn236nlvt2d = 2, rucn236nvs2d = 4 )

  INTEGER :: rucn236levs2d(rucn236nlvt2d) ! Level or Layer identification
                                          ! of RUC
  INTEGER :: rucn236var_id2d(rucn236nvs2d,rucn236nlvt2d)

  DATA rucn236levs2d  / 001, 105/       ! Ground surface
  DATA rucn236var_id2d/ 223, 66, 85, 86,                &
                         11, 51, 33, 34 /

! 223 = Plant canopy surface water (kg/m**2)
!  66 = Water equivalent of snow depth (m)
!  85 = Soil temperature (K)
!  86 = officially soil moisture content (kg/m**2),
!       but arrays contain volumetric (fractional) moisture,
!       i.e. should have used #144

  INTEGER :: rucn236nvs3d    ! Number of RUC 3-D variables to be read
  INTEGER :: rucn236nlvt3d   ! Number of RUC 3-D level types to be read
  PARAMETER ( rucn236nlvt3d = 2, rucn236nvs3d = 11 )

  INTEGER :: rucn236levs3d(rucn236nlvt3d)  ! Number of RUC 3-D variables
                                           ! to be read
  INTEGER :: rucn236var_id3d(rucn236nvs3d, rucn236nlvt3d)

  DATA rucn236levs3d/ 109, 111 /
  DATA rucn236var_id3d/   1,   7, 189, 53,  33,  34,    &
                 153, 170, 178, 171, 179,               & ! Hybrid levels
                  85,  86, 255, 255, 255, 255,          &
                 255, 255, 255, 255, 255/                 ! Soil levels
  !   1 = Pressure (Pa)
  !   7 = height (m)
  ! 189 = virtual potential temperature (K)
  !  53 = water vapor mixing ratio (kg/kg)
  !  33 = u wind (m/s)
  !  34 = v wind (m/s)
  ! 153 = cloud water mixing ratio (kg/kg)
  ! 170 = rain water mixing ratio (kg/kg)
  ! 178 = ice mixing ratio (kg/kg)
  ! 171 = snow mixing ratio (kg/kg)
  ! 179 = graupel mixing ratio (kg/kg)
  !  85 = soil temperature (K)
  ! 144 = soil volumetric moisture content (fraction)
  !  86 = officially soil moisture content (kg/m**2),
  !       but arrays contain volumetric (fractional) moisture,
  !       i.e. should have used #144
  ! 255 = undefined
!
!-----------------------------------------------------------------------
!
!  Parameters for NCEP RUC2 Isobaric GRIB (Grid #236)
!
!-----------------------------------------------------------------------
!
!  INTEGER :: rucp236grid  ! Grid identification of ruc236
  INTEGER, PARAMETER :: rucp236grid = 236  ! RUC Lambert Conformal
  INTEGER, PARAMETER :: rucp252grid = 252  ! RUC Lambert Conformal

  INTEGER :: rucp236proj  ! Map projection ID for ruc236
  PARAMETER ( rucp236proj = 3 )   ! for Lambert Conformal

  INTEGER :: rucp236nlvt2d      ! Number of rucp236 level types to be
                                !processed
  INTEGER :: rucp236nvs2d       ! Number of rucp236 2-D variables to be
                                !read
  PARAMETER ( rucp236nvs2d = 4, rucp236nlvt2d = 2 )

  INTEGER :: rucp236levs2d(rucp236nlvt2d) ! Level or Layer identification
                                          ! of RUC
  INTEGER :: rucp236var_id2d(rucp236nvs2d,rucp236nlvt2d)

  DATA rucp236levs2d  / 001, 105/       ! Ground surface
  DATA rucp236var_id2d/ 1,  7, 66, 255,                                 &
                        11, 52, 33, 34/
  ! 1 = Pressure (Pa)
  !  7 = Geometric height (m)
  !  66= Water equivalent snow depth (m)
  ! 255= undefined
  !  11= T_2m (C)
  !  52= RH_2m
  !  33= U_10m
  !  34= V_10m

  INTEGER :: rucp236nvs3d       ! Number of RUC 3-D variables to be read
  INTEGER :: rucp236nlvt3d      ! Number of RUC 3-D level types to be
                                ! read
  PARAMETER ( rucp236nvs3d = 5, rucp236nlvt3d= 1 )

  INTEGER :: rucp236levs3d(rucp236nlvt3d)       ! Number of RUC 3-D
                                                ! variables to be read
  INTEGER :: rucp236var_id3d(rucp236nvs3d, rucp236nlvt3d)

  DATA rucp236levs3d   / 100 /
  DATA rucp236var_id3d /   7,  11,  52,  33,  34/

  !   7 = Pressure surface height (m)
  !  11 = Temperature (K)
  !  52 = Relative humidity (%)
  !  33 = u wind (m/s) (grid relative)
  !  34 = v wind (m/s) (grid relative)

!-----------------------------------------------------------------------
!
!  Parameters for the lastest Rapid Refresh and High-resolution rapid
!  refresh (HRRR) data sets (actually, NCEP RUC Native Coordinate GRIB in
!  Grid #130)
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: rucn130grid = 130     ! RUC Lambert conformal

  INTEGER, PARAMETER :: rucn130proj = 3       ! Lambert conformal projection

  INTEGER, PARAMETER :: rucn130nlvt2d = 2, rucn130nvs2d = 7
   ! Number of RUC level types to be processed
   ! Number of RUC 2-D variables to be read

  INTEGER :: rucn130levs2d(rucn130nlvt2d) ! Level or Layer identification
                                          ! of RUC
  INTEGER :: rucn130var_id2d(rucn130nvs2d,rucn130nlvt2d)

  DATA rucn130levs2d  / 001, 105/         ! Ground surface & specified height
  DATA rucn130var_id2d/ 62, 63, 66, 81,  85, 144, 223,   & ! Ground surface
                        11, 33, 34, 53, 255, 255, 255 /   ! 2m T & Q, 10m U & V

!  62  = Large scale precipitation (non-conv., kg/m2)
!  63  = Convective precipitation
!  66  = snow depth (m)
!  81  = Land cover (land=1, sea=0)
!  85  = Soil temperature (K)
!  144 = Volumetric soil moisture content
!  223 = Plant canopy surface water (kg/m**2)

  INTEGER :: rucn130nvs3d    ! Number of RUC 3-D variables to be read
  INTEGER :: rucn130nlvt3d   ! Number of RUC 3-D level types to be read
  PARAMETER ( rucn130nlvt3d = 2, rucn130nvs3d = 11 )

  INTEGER :: rucn130levs3d(rucn130nlvt3d)  ! Number of RUC 3-D variables
                                           ! to be read
  INTEGER :: rucn130var_id3d(rucn130nvs3d, rucn130nlvt3d)

  DATA rucn130levs3d/ 109, 111 /
  DATA rucn130var_id3d/    1,   7, 189,  53,  33,  34,    &
                         153, 170, 178, 171, 179,         & ! Hybrid levels
                          85, 144, 255, 255, 255, 255,    &
                         255, 255, 255, 255, 255/          ! Soil levels
  !   1 = Pressure (Pa)
  !   7 = height (m)
  ! 189 = virtual potential temperature (K)
  !  53 = water vapor mixing ratio (kg/kg)
  !  33 = u wind (m/s)
  !  34 = v wind (m/s)
  ! 153 = cloud water mixing ratio (kg/kg)
  ! 170 = rain water mixing ratio (kg/kg)
  ! 178 = ice mixing ratio (kg/kg)
  ! 171 = snow mixing ratio (kg/kg)
  ! 179 = graupel mixing ratio (kg/kg)
  !  85 = soil temperature (K)
  ! 144 = soil volumetric moisture content (fraction)
  ! 255 = undefined

!-----------------------------------------------------------------------
!
!  Parameters for ECMWF data:
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: ecmf128grid = 255  ! Grid identification of ETA
  INTEGER, PARAMETER :: ecmf128proj = 0    ! Map projection ID

  INTEGER, PARAMETER :: ecmf128nlvt2d = 2  ! Number of 2-D level types to be processed
  INTEGER, PARAMETER :: ecmf128nvs2d  = 8  ! Number of 2-D variables   to be processed

  INTEGER :: ecmf128levs2d(ecmf128nlvt2d) ! 2-D Level identification
  INTEGER :: ecmf128var_id2d(ecmf128nvs2d,ecmf128nlvt2d)

  DATA ecmf128levs2d  / 001, 112 /   ! At ground surface and meters above
  DATA ecmf128var_id2d/  34, 129, 134, 165, 166, 167, 172, 235,         &
                         39,  40,  41,  42, 139, 170, 183, 236  /

 !  34 = Sea Surface temperature (K)
 ! 129 = Geopotential height (m**2 s**-2)
 ! 134 = Surface pressure (Pa)
 ! 235 = Skin temperature (K)
 ! 172 = Land/sea mask
 ! 165 = U_10m (m/s)
 ! 166 = V_10m (m/s)
 ! 167 = T_2m (K)

! 112 layer between two depths below land surface
!
!  39 = Volumetric soil water (  0-7   cm, m**3/m**3)
!  40 = Volumetric soil water (  7-28  cm, m**3/m**3)
!  41 = Volumetric soil water ( 28-100 cm, m**3/m**3)
!  42 = Volumetric soil water (100-255 cm, m**3/m**3)
! 139 = Soil Temperature (  0-7   cm, K)
! 170 = Soil Temperature (  7-28  cm, K)
! 183 = Soil Temperature ( 28-100 cm, K)
! 236 = Soil Temperature (100-255 cm, K)
!

  INTEGER, PARAMETER :: ecmf128nlvt3d = 1  ! Number of 3-D level types to be processed
  INTEGER, PARAMETER :: ecmf128nvs3d  = 5  ! Number of 3-D variables to be processed

  INTEGER :: ecmf128levs3d(ecmf128nlvt3d)   ! 3-D Level identification of ETA
  INTEGER :: ecmf128var_id3d (ecmf128nvs3d,ecmf128nlvt3d)

  DATA ecmf128levs3d / 100 /                ! isobaric vertical coordinates
  DATA ecmf128var_id3d/     130, 131, 132, 156, 157 /

! 130 = Temperature (K)
! 131 = u wind (m/s)
! 132 = v wind (m/s)
! 156 = Geopotential height (gpm)
! 157 = Relative humidity (%)
!
!-----------------------------------------------------------------------
!
!  General definitions for GRIB parameters and variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nprods     ! Maximum number of GRIB products (2-D arrays)
  PARAMETER ( nprods = 2000 )

  INTEGER :: n2dvs      ! Maximum number of 2-D variables   to be read
  INTEGER :: n2dlvt     ! Maximum number of 2-D level types to be read
  PARAMETER ( n2dlvt = 4, n2dvs = 8 ) ! KFY

  INTEGER :: n3dvs      ! Maximum number of 3-D variables   to be read
  INTEGER :: n3dlvt     ! Maximum number of 3-D level types to be read
  PARAMETER ( n3dlvt = 2, n3dvs = 11 )

  INTEGER :: rcbytes(nprods)  ! record length in bytes
  INTEGER :: rcstart(nprods)  ! record starting byte in a GRIB fil

  INTEGER :: var_nr2d(n2dvs,n2dlvt)
                              ! number of record for all 2-D variables
  INTEGER :: var_nr3d(n3dvs,n3dlvt)  ! number of record for each 3-D variable

!  integer var_lev3d(nzgrb,n3dvs,n3dlvt)
                               ! Levels (hybrid) for each 3-D variable

!  real rcdata(nxygrb)      ! temporary data array
!  real var_grb2d(nxgrb,nygrb,n2dvs,n2dlvt) ! GRIB variables
!  real var_grb3d(nxgrb,nygrb,nzgrb,n3dvs,n3dlvt) ! GRIB 3-D variables
!
!-----------------------------------------------------------------------
!
!  The following variables and arrays must be set before calling
!  subroutine RDNMCGRB.
!
!  (See subroutine GETNMCRUC for RUC and GETNMCETA for ETA in file
!  getextd3d.f as examples.)
!
!-----------------------------------------------------------------------
!
  INTEGER :: gridtyp    ! Grid ID (working variable)
  INTEGER :: mproj_grb  ! Map projection ID

  INTEGER :: n2dvars    ! Number of 2-D variables (working variable)
  INTEGER :: n2dlvtps   ! Number of 2-D level types (working variable)
  INTEGER :: levtyp2d(n2dlvt)     ! 2-D Level IDs (working variable)
  INTEGER :: var_id2d(n2dvs,n2dlvt)
                         ! Working array for 2-D variable IDs

  INTEGER :: n3dvars    ! Number of 3-D variables (working variable)
  INTEGER :: n3dlvtps   ! Number of 3-D level types (working variable)
  INTEGER :: levtyp3d(n3dlvt)   ! 3-D Level IDs (working variable)
  INTEGER :: var_id3d(n3dvs,n3dlvt)
                         ! Working array for 3-D variable IDs

  COMMON /grib001/gridtyp, mproj_grb
  COMMON /grib002/n2dvars,n2dlvtps,levtyp2d,var_id2d
  COMMON /grib003/n3dvars,n3dlvtps,levtyp3d,var_id3d
  COMMON /grib004/rcbytes,rcstart,var_nr2d,var_nr3d

!  COMMON /grib004/rcbytes,rcstart,var_nr2d,var_nr3d,var_lev3d
!  COMMON /grib005/rcdata, var_grb2d,var_grb3d
!
!-----------------------------------------------------------------------
!
!  GRIB section parameters:
!
!-----------------------------------------------------------------------
!
  INTEGER :: mxnbufsz               ! Size of GRIB buffer
  INTEGER :: nbufsz                 ! Size of GRIB buffer
!  PARAMETER ( mxnbufsz = 800000, nbufsz = 4*(mxnbufsz/4) )
  PARAMETER ( mxnbufsz = 4800000, nbufsz = 4*(mxnbufsz/4) )

  INTEGER :: npdsz, ngdsz           ! sizes of PDS, GDS
  PARAMETER ( npdsz = 100, ngdsz = 200 )

  INTEGER :: ipdsz, igdsz           ! sizes of IPDS, IGDS
  PARAMETER ( ipdsz = 100, igdsz = 200 )

  INTEGER :: mptrs(15)  ! Integer array containing storage for
                        ! following parameters:
                        ! (1)  total length of grib message
                        ! (2)  length of indicator (SECTION  0)
                        ! (3)  length of pds       (SECTION  1)
                        ! (4)  length of gds       (SECTION  2)
                        ! (5)  length of bms       (SECTION  3)
                        ! (6)  length of bds       (SECTION  4)
                        ! (7)  value of current byte
                        ! (8)  bit pointer
                        ! (9)  GRIB start bit nr
                        ! (10) GRIB/GRID element count
                        ! (11) nr unused bits at end of SECTION 3
                        ! (12) bit map flag (copy of BMS octets 5,6)
                        ! (13) nr unused bits at end of SECTION 2
                        ! (14) BDS flags (right adj copyat end of SECTION 2
                        ! (14) BDS flags (right adj copy of octet 4)
                        ! (15) nr unused bits at end of SECTION 4

  INTEGER :: ipds(ipdsz)       ! PDS integer array
  INTEGER :: igds(igdsz)       ! GDS integer array

!  logical ibms(nxygrb)      ! BMS logical array for data bit map

  INTEGER :: ibdshd(4)         ! IBDS header

  CHARACTER (LEN=1) :: pds(npdsz)    ! PDS ( GRIB Section 1)
  CHARACTER (LEN=1) :: gds(npdsz)    ! GDS ( GRIB Section 2)
  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS ( GRIB Section 4)
  INTEGER :: ibds(nbufsz/4)    ! identical to BDS

  EQUIVALENCE ( bds,ibds )

  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Buffer to carry GRIB messages
!
!-----------------------------------------------------------------------
!
!  Parameters which define WMO GRIB Edition 2 format:
!
!  WARNING:  Don''t muck with these parameters unless the WMO changes
!            the definition of GRIB Edition 2 format!
!
!-----------------------------------------------------------------------
!
  INTEGER :: nctrs    ! Maximum number of meteorological centers
  INTEGER :: nsbctrs  ! Maximum number of meteorological sub-centers
  INTEGER :: nmods    ! Maximum number of numerical models/processes
  INTEGER :: ntbvars  ! Maximum number of variables in GRIB tables
  INTEGER :: ntimes   ! Maximum number of time units & range indicators
  INTEGER :: nproj    ! Maximum number of map projection names
  INTEGER :: nscan    ! Maximum number of directional scanning modes

  PARAMETER ( nctrs   = 99  )
  PARAMETER ( nsbctrs = 170 )
  PARAMETER ( nmods   = 153 )
  PARAMETER ( ntbvars = 255 )
  PARAMETER ( ntimes  = 254 )
  PARAMETER ( nproj   = 254 )
  PARAMETER ( nscan   = 8   )
!
!-----------------------------------------------------------------------
!
!  Table 0 of national/international originating centers, designated
!  by Parameter Description Section (PDS) in Section 1 of GRIB
!  message:
!     Center     = PDS octet # 5
!     Sub-center = PDS octet #26
!
!  NOTE:  (1) Center 1 is a sub-center value for Center 7, the US
!             National Meteorological Center.
!         (2) Centers 150 - 170 are sub-center values for Center 9,
!             the US National Weather Service Field Stations.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: center(0:nsbctrs)       ! Names of met. sub-centers
  CHARACTER (LEN=42) :: sbcenter(0:nsbctrs)     ! Names of met. sub-centers
!
!-----------------------------------------------------------------------
!
!  Table A of numerical model or generating process, designated by
!  Parameter Description Section (PDS) octet #6 in Section 1 of GRIB
!  message.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: model(nmods)         ! Names of numerical models
!
!-----------------------------------------------------------------------
!
!  Table 2 of GRIB variables, designated by Parameter Description
!  Section (PDS) octet #9 in Section 1 of GRIB message.
!
!  NOTE:  Variables 128 - 254 are RESERVED for use by the originating
!         center. NWS/NMC usage is assumed.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: varname(0:ntbvars)       ! Names  of GRIB variables
  CHARACTER (LEN=42) :: varunit(0:ntbvars)       ! Units  of GRIB variables
  CHARACTER (LEN=42) :: varabbr(0:ntbvars)       ! Abbrev of GRIB variables
!
!-----------------------------------------------------------------------
!
!  Table 3 of levels/layers where variables are defined, designated
!  by Parameter Description Section (PDS) octets #10-12 in Section 1
!  of GRIB message.
!
!  NOTE: The numbering allows for additions within this framework:
!
!        Octet #10   Numerical Precision
!        ---------   -------------------
!          100-119   normal
!          120-139   high
!          140-159   mixed
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: varlvl(0:ntbvars)       ! Levels of GRIB variables
!
!-----------------------------------------------------------------------
!
!  Tables 4, 5 of time units and range indicators, designated by
!  Parameter Description Section (PDS) in Section 1 of GRIB message:
!     Forecast Time Unit   = PDS octet #18
!     Time Range Indicator = PDS octet #21
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: timeunit(0:ntimes)      ! Forecast time units
  CHARACTER (LEN=42) :: timerang(0:ntimes)      ! Time range indicators
!
!-----------------------------------------------------------------------
!
!  Table 6 of map projection type, designated by Grid Description
!  Section (GDS) in Section 2 of GRIB message:
!     Projection = GDS octet #8
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: projs(0:nproj)          ! Names of map projections
!
!-----------------------------------------------------------------------
!
!  Table 8 of directional scanning modes, designated by Grid
!  Description Section (GDS) in Section 2 of GRIB message:
!     Scanning Mode = GDS octet #28
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: scans(0:1,nscan)    ! Directional scanning modes
!
!-----------------------------------------------------------------------
!
!  End of GRIBCST.INC
!
!-----------------------------------------------------------------------
!
