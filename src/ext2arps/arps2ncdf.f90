PROGRAM arps2ncdf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARSP2NCDF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Converts ARPS history files to a netCDF format file.
!
!  Reads in a history file produced by ARPS in any ARPS format.
!  Converts variables and interpolates to a fixed set of pressure
!  levels.  Sets up variables needed for LDADS/AWIPS D2d display
!  of data.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  02/19/1999
!
!  MODIFICATION HISTORY:
!  08/13/1999
!  Fixed MSLP calculation (J. Levit, K. Brewster)
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
!
!  11/22/2003 (K. Brewster)
!  Added logic to allow appending new time levels to existing file.
!  Includes additions to input file.  Corrected itim1970 constant.
!
!  12/02/2003 (T. Oram, NWS SFMG)
!  Modifications for AWIPS compatibility plus comments to help
!  myself understand what is going on
!  a.   This program works with the arps.cdl netCDF cdl file that is used
!       for the AWIPS localization process.  The arps.cdl file is new.
!  b.   Changed the tracking of the current and maximum number of files that
!       can be appended to a netCDF file using the record and n_valtimes
!       dimensions in the netCDF file.
!  c.   modified arps2ncdf.input for new variables and modified write of
!       1-d variables (valtime, reftime, valtimeMINUSreftime, etc) to only
!       be written once for netcdf append files.  This assumes that all forecast
!       times are known when the first history file is read and converted to
!       netCDF.  Plan to move the number of pressure levels and the pressure level
!       array into arps2ncdf.input for run-time determination rather than
!       netcdf.inc for compile-time determination.
!  d.   The left most dimension of the inventory variables also was
!       changed to n_valtimes rather than record.  AWIPS chokes
!       without a full inventory array for all valid times.  The new dimension of
!       the inventory array matches the AWIPS convention.
!  e.   added surface level to the 3-dimensional array to match the awips convention
!  f.   added cw, cice, and tke to 3-dimensional array (tdo, 02/02/2004)
!  g.   removed 2-d fields not included in the dataFieldTable.txt file.
!       also changed reftime to remove the record dimension to match awips (tdo, 02/05/2004)
!
!  02/22/2005 (Y. Wang)
!  Moved the pressure level specification to arps2ncdf.input for run-time
!  determination as planed by T. Oram.
!
!  Cleaned up documents and formats.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of soil model (m)
!
!    w        vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!    uprt     perturbation x velocity component (m/s)
!    vprt     perturbation y velocity component (m/s)
!    wprt     perturbation z velocity component (m/s)
!
!    qv       water vapor mixing ratio (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil Temperature (K)
!    qsoil    Soil moisture
!    wetcanp  Canopy water amount
!
!    rain     Total rain (raing + rainc)
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz          ! Grid dimensions.
  INTEGER :: nzsoil            ! Soil levels
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER, PARAMETER  :: nhisfile_max = 200
  INTEGER             :: hinfmt,nhisfile,lengbf,lenfil
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'netcdf.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x(:)       ! The x-coord. of the physical and
                                  ! computational grid.
                                  ! Defined at u-point.
  REAL, ALLOCATABLE :: y(:)       ! The y-coord. of the physical and
                                  ! computational grid.
                                  ! Defined at v-point.
  REAL, ALLOCATABLE :: z(:)       ! The z-coord. of the computational
                                  ! grid. Defined at w-point.

  REAL, ALLOCATABLE :: zp(:,:,:)  ! The height of the terrain.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height of the soil model.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar (:,:,:,:) ! Microphysics scalar array
  REAL, ALLOCATABLE :: qv     (:,:,:) ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qi     (:,:,:)

  REAL, ALLOCATABLE :: tke    (:,:,:) ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh    (:,:,:) ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv    (:,:,:) ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: ubar   (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:) ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:) ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:) ! Base state water vapor specific

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct(:,:,:) ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: lai    (:,:)   ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)   ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)   ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:) ! Soil Temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:) ! Soil Moisture
  REAL, ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (m)

  REAL, ALLOCATABLE :: rain (:,:)     ! Total rainfall
  REAL, ALLOCATABLE :: raing(:,:)     ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)     ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:) ! precipitation rate (kg/(m**2*s))
                                      ! prcrate(1,1,1) = total precip. rate
                                      ! prcrate(1,1,2) = grid scale precip. rate
                                      ! prcrate(1,1,3) = cumulus precip. rate
                                      ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)  ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)    ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)    ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)  ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)   ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: e_mb   (:,:,:)    ! vapor pressure in mb
  REAL, ALLOCATABLE :: mix    (:,:,:)    ! total mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: esat_mb(:,:,:)    ! saturation vapor pressure in mb
  REAL, ALLOCATABLE :: rh     (:,:,:)    ! Relative humidity in %
  REAL, ALLOCATABLE :: t_dew  (:,:,:)    ! dewpoint temp. in degrees K

!
!-----------------------------------------------------------------------
!
!  2-D stability index arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: lcl(:,:)   ! lifting condensation level
  REAL, ALLOCATABLE :: lfc(:,:)   ! level of free convection
  REAL, ALLOCATABLE :: el(:,:)    ! equilibrium level
  REAL, ALLOCATABLE :: twdf(:,:)  ! max. wet bulb pot. temp. difference
  REAL, ALLOCATABLE :: li(:,:)    ! lifted index
  REAL, ALLOCATABLE :: pbe(:,:)   ! CAPE
  REAL, ALLOCATABLE :: mbe(:,:)   ! Moist CAPE
  REAL, ALLOCATABLE :: nbe(:,:)   ! CIN
  REAL, ALLOCATABLE :: tcap(:,:)  ! Cap Strength

  REAL, ALLOCATABLE :: heli(:,:)  ! helicity
  REAL, ALLOCATABLE :: brn(:,:)   ! Bulk Richardson Number (Weisman and Klemp)
  REAL, ALLOCATABLE :: brnu(:,:)  ! Shear parameter of BRN, "U"
  REAL, ALLOCATABLE :: srlfl(:,:) ! storm-relative low-level flow (0-2km AGL)
  REAL, ALLOCATABLE :: srmfl(:,:) ! storm-relative mid-level flow (2-9km AGL)
  REAL, ALLOCATABLE :: shr37(:,:) ! 7km - 3km wind shear
  REAL, ALLOCATABLE :: ustrm(:,:) ! Estimated storm motion (modified Bob Johns)
  REAL, ALLOCATABLE :: vstrm(:,:) ! Estimated storm motion (modified Bob Johns)
  REAL, ALLOCATABLE :: blcon(:,:) ! Boundary layer convergence
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:,:)
  REAL, ALLOCATABLE :: tem2d1(:,:)
  REAL, ALLOCATABLE :: tem2d2(:,:)
  REAL, ALLOCATABLE :: tem2d3(:,:)
  REAL, ALLOCATABLE :: tem2d4(:,:)

  REAL, ALLOCATABLE :: wrk1(:),wrk2(:),wrk3(:),wrk4 (:),wrk5 (:),wrk6 (:)
  REAL, ALLOCATABLE :: wrk7(:),wrk8(:),wrk9(:),wrk10(:),wrk11(:),wrk12(:)

  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
  REAL, ALLOCATABLE :: zps_km(:,:,:)

! J.Case, ENSCO Inc. (9/29/2004)
! Additional variables.
! Note: radar_comp is maximum column reflecitivity and
!       echo_top is maximum height of 18 dBZ echo in a column.

  REAL, ALLOCATABLE :: radar_3d(:,:,:),radar_comp(:,:),echo_top(:,:)

  REAL, ALLOCATABLE :: mf2d(:,:)
  REAL, ALLOCATABLE :: coriolis(:,:)
  REAL, ALLOCATABLE :: spacing(:,:)
!
!-----------------------------------------------------------------------
!
!  Misc ARPS variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nch

  INTEGER :: first18                    ! J.Case, ENSCO Inc. (9/29/2004)

  REAL    :: time
  REAL    :: latnot(2)
  LOGICAL :: init
!
!-----------------------------------------------------------------------
!
!  netCDF dimensions
!
!  nprlvl       defined in arps2cdf.inc is the number of pressure levels
!  nxcdf        number of points in x direction in netcdf file
!  nycdf        number of points in y direction in netcdf file
!  nrecord      number of records in netcdf file
!  ntotaltime   number of valid times that are to be or can be written into file
!  ncdflvls     number of levels in 3-d grids in file; nprlvl + 1 (surface)
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mxtime = 40
  INTEGER :: nxcdf,nycdf,nrecord,nvaltimes,ntotaltimes

  INTEGER, PARAMETER :: nprlvl_max = 128
  INTEGER            :: iprlvl(nprlvl_max)
  INTEGER            :: nprlvl, ncdflvls                     ! = nprlvl+1

  INTEGER, PARAMETER :: n3dvar = 10, n2dvar = 16, nstvar = 3
  INTEGER, PARAMETER :: nvartot = (n3dvar+n2dvar)
                        ! J.Case, ENSCO Inc. (9/29/2004) -- Added 1 new 3D
                        ! variable (rr) and 2 new 2D variables (cxr, mret)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istatus,ncid
  INTEGER :: dimidrec,dimidtot,dimidnt,dimidnav
  INTEGER :: dimidnx,dimidny,dimid1,dimidnz,dimchrlvl,dimnamlen
  INTEGER :: vtmrefid,valtimid,reftimid,origid,modelid,ntwrtid

  INTEGER :: vecistr(2)
  INTEGER :: planeistr(4)
  INTEGER :: volistart(4)
  DATA vecistr /1,1/
  DATA planeistr /1,1,1,1/
  DATA volistart /1,1,1,1/

  INTEGER :: planedim(4)
  INTEGER :: voldim(4)
  INTEGER :: lvldim(2)
  INTEGER :: invdim(2)
  INTEGER :: planelen(4)
  INTEGER :: lvllen(2)
  INTEGER :: sfclen(2)
  INTEGER :: vollen(4)
  INTEGER :: invlen(2)
  INTEGER :: ntwrt

  REAL :: vrange(2)
!
!-----------------------------------------------------------------------
!
!  netCDF output variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: v3did(n3dvar)
  INTEGER :: v3dlvlid(n3dvar)
  INTEGER :: v2dinvid(n2dvar)
  INTEGER :: v3dinvid(n3dvar)
  INTEGER :: v2did(n2dvar)
  INTEGER :: vstid(nstvar)
  INTEGER :: v2dlvlid(n2dvar)
  INTEGER :: vtmref(mxtime)
  REAL*8  valtime(mxtime)
  REAL*8  reftime(mxtime)

  CHARACTER (LEN=4) :: v3dnam(n3dvar)
  CHARACTER (LEN=4) :: v2dnam(n2dvar)
  CHARACTER (LEN=14) :: vstnam(nstvar)
  CHARACTER (LEN=10) :: v3dlvl(nprlvl_max)
  CHARACTER (LEN=10) :: v2dlvl
  CHARACTER (LEN=30) :: v3dlngnam(n3dvar)
  CHARACTER (LEN=30) :: v2dlngnam(n2dvar)
  CHARACTER (LEN=30) :: vstlngnam(nstvar)
  CHARACTER (LEN=30) :: v3dunits(n3dvar)
  CHARACTER (LEN=30) :: v2dunits(n2dvar)
  CHARACTER (LEN=30) :: vstunits(nstvar)
  CHARACTER (LEN=nprlvl_max) :: inventory
  CHARACTER (LEN=1) :: inventory_sfc
  CHARACTER (LEN=1) :: inventory_blank
  REAL :: v3dmin(n3dvar)
  REAL :: v2dmin(n2dvar)
  REAL :: v3dmax(n3dvar)
  REAL :: v2dmax(n2dvar)

  CHARACTER (LEN=10) :: v2dlvlnam
  CHARACTER (LEN=10) :: v3dlvlnam
  CHARACTER (LEN=13) :: v2dinv
  CHARACTER (LEN=13) :: v3dinv

  DATA vtmref /mxtime*0/
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: idxgh,idxrh,idxt,idxuw,idxvw,idxww,idxcw,idxci,idxtk
                                         ! J.Case, ENSCO Inc. (9/29/2004)
  INTEGER :: idxrr
  PARAMETER (idxgh=1,idxrh=2,idxt=3,                                    &
             idxuw=4,idxvw=5,idxww=6,                                   &
             idxcw=7,idxci=8,idxtk=9,idxrr=10)
!
! J.Case, cont. (9/29/2004) -- Added radar reflectivity to 3D variables (rr).

  DATA v3dnam /'gh','rh','t','uw','vw','ww','cw','cice','tke','rr'/
  DATA v3dlngnam /'Geopotential Height',                                &
                  'Relative Humidity',                                  &
                  'Temperature',                                        &
                  'u Wind Component',                                   &
                  'v Wind Component',                                   &
                  'w Wind Component',                                   &
                  'Cloud Water',                                        &
                  'Cloud Ice',                                          &
                  'Turbulent Kinetic Energy',                           &
                  'Radar Reflectivity'/
  DATA v3dunits /'geopotential meters',                                 &
                 'percent',                                             &
                 'degrees K',                                           &
                 'meters/second',                                       &
                 'meters/second',                                       &
                 'meters/second',                                       &
                 'kg/kg',                                               &
                 'kg/kg',                                               &
                 '(meters/second)**2',                                  &
                 'dBZ'/
  DATA v3dmin /    0.,  0.,150.,-300.,-300.,-100.,0.,0.,0.,-200./
  DATA v3dmax /40000.,110.,400., 300., 300., 100.,1.,1.,1000.,200./

  INTEGER :: idxmslp,idxp,idxlgsp,idxcp,idxtp
  INTEGER :: idxdpt
  INTEGER :: idxsh,idxept,idxsli,idxcape,idxcin
  INTEGER :: idxustrm,idxvstrm
  INTEGER :: idxheli

         ! J.Case, ENSCO Inc. (9/29/2004) -- new 2D variables cxr and mret
  INTEGER :: idxcxr,idxmret

  PARAMETER (idxmslp=1,                                                 &
             idxp   =2,                                                 &
             idxlgsp=3,                                                 &
             idxcp  =4,                                                 &
             idxtp  =5,                                                 &
             idxdpt =6,                                                 &
             idxsh  =7,                                                 &
             idxept =8,                                                 &
             idxsli =9,                                                 &
             idxcape=10,                                                &
             idxcin =11,                                                &
             idxustrm=12,                                               &
             idxvstrm=13,                                               &
             idxheli =14,                                               &
             idxcxr  =15,                                               &
             idxmret =16)

!             idxsfu =6,                                                 &
!             idxsfv =7,                                                 &
!             idxsft =8,                                                 &
!             idxlcl =15,                                                &
!             idxlfc =16,                                                &
!             idxel  =17,                                                &
!             idxtwdf =18,                                               &
!             idxtcap =19,                                               &
!             idxshr37=20,                                               &
!             idxsrlfl=23,                                               &
!             idxsrmfl=24,                                               &
!             idxbrn  =26,                                               &
!             idxbrnu =27,                                               &
!             idxblcon=28)


! these 2d variable names must match cdl names in the AWIPS
! dataFieldTable.txt file
!  DATA v2dnam /'mslp','p','lgsp','cp','tp','sfu','sfv',                 &
!               'sft','dpt','sh','ept','sli','cape','cin','lcl',         &
!               'lfc','el','twdf','tcap','sh37','ustm','vstm',           &
!               'srlf','srmf','heli','brn','brnu','blcn'/

!JLC (8/23/2004) -- added cxr,mret to 2D variables.

  DATA v2dnam /'mslp','p','lgsp','cp','tp',                             &
               'dpt','sh','ept','sli','cape','cin',                     &
               'ustm','vstm',                                           &
               'heli','cxr','mret'/


  DATA v2dlngnam /'Mean Sea Level Pressure',                            &
                  'Surface Pressure',                                   &
                  'Grid Scale Precipitation',                           &
                  'Convective Precipitation',                           &
                  'Total Precipitation',                                &
                  'Surface Dew Point',                                  &
                  'Specific Humidity',                                  &
                  'Theta-e',                                            &
                  'Surface Lifted Index',                               &
                  'CAPE',                                               &
                  'CIN',                                                &
                  'Est Storm Motion u comp',                            &
                  'Est Storm Motion v comp',                            &
                  'Storm Rel Helicity',                                 &
                  'Column Max Reflectivity',                            &
                  'Maximum Radar Echo Tops'/

!                  'Surface u Wind Component',                           &
!                  'Surface v Wind Component',                           &
!                  'Surface Temperature',                                &
!                  'Lifting Condensation Lvl',                           &
!                  'Level of Free Convection',                           &
!                  'Equilibrium Level',                                  &
!                  'Max Wet-Bulb Temp Diff',                             &
!                  'Cap Strength',                                       &
!                  '3-7km Shear',                                        &
!                  'Storm Rel Low-Level Flow',                           &
!                  'Storm Rel Mid-Level Flow',                           &
!                  'Bulk Richardson Number',                             &
!                  'BRN Shear u',                                        &
!                  'Bound-Layer Conv x 1000'/

! J.Case, cont. (9/29/2004)

  DATA v2dunits /'Pascals','Pascals','mm','mm','mm',                    &
                 'degrees K','kg/kg',                                   &
                 'degrees K','degrees C','J/kg','J/kg',                 &
                 'meters/second','meters/second',                       &
                 'm2/s2','dBZ','meters'/

  DATA v2dmin / 80000.,50000.,0.,0.,0.,                                 &
                -100.,0.,                                               &
                200.,-100.,0.,0.,                                       &
                -100.,-100.,                                            &
                -900.,0.,0./

  DATA v2dmax / 110000.,120000.,1000.,1000.,1000.,                      &
                100.,1.,                                                &
                500.,100.,10000.,10000.,                                &
                100.,100.,                                              &
                900.,100.,20000./

  INTEGER :: idxtop,idxcor,idxspa
  PARAMETER (idxtop=1,idxcor=2,idxspa=3)
  DATA vstnam /'staticTopo',                                            &
               'staticCoriolis',                                        &
               'staticSpacing'/
  DATA vstlngnam /'Topography',                                         &
                  'Coriolis parameter',                                 &
                  'Grid spacing'/
  DATA vstunits /'meters',                                              &
                 '/second',                                             &
                 'meters'/

  CHARACTER (LEN=20) :: cproj
!
!  The following are the projections supported by LDAD/AWIPS
!  in their ordering/naming convention
!
  CHARACTER (LEN=25) :: chproj(17)
  DATA chproj /'STEREOGRAPHIC',                                         &
               'ORTHOGRAPHIC',                                          &
               'LAMBERT_CONFORMAL',                                     &
               'AZIMUTHAL_EQUAL_AREA',                                  &
               'GNOMONIC',                                              &
               'AZIMUTHAL_EQUIDISTANT',                                 &
               'SATELLITE_VIEW',                                        &
               'CYLINDRICAL_EQUIDISTANT',                               &
               'MERCATOR',                                              &
               'MOLLWEIDE',                                             &
               'NULL',                                                  &
               'NULL',                                                  &
               'NULL',                                                  &
               'NULL',                                                  &
               'NULL',                                                  &
               'AZIMUTH_RANGE',                                         &
               'RADAR_AZ_RAN'/
!
!  kproj is the index in the LDAD/AWIPS map projection list
!  for the mapproj number in the ARPS file
!
  INTEGER :: kproj(3)
  DATA kproj /1,3,9/

  REAL, ALLOCATABLE :: static2d(:,:,:)
  REAL, ALLOCATABLE :: cdf2d(:,:,:,:)
  REAL, ALLOCATABLE :: cdf3d(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Constants
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER  :: itim1970=315622800
  CHARACTER(LEN=8),   PARAMETER :: cdldate  = '20031202'
  CHARACTER(LEN=40),  PARAMETER :: depictor = 'ARPS'
  CHARACTER(LEN=132), PARAMETER :: origin   = 'CAPS'
  CHARACTER(LEN=132), PARAMETER :: model    = 'ARPS 5.3'
!
!-----------------------------------------------------------------------
!
!  Namelists
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: ipktyp,nbits,ncdapnd,ncdntime,ncdtinterval
  NAMELIST /ncdf_options/ nprlvl,iprlvl,cdfname,ipktyp,nbits,ncdapnd, &
                          ncdntime,ncdtinterval
!
!-----------------------------------------------------------------------
!
!  External functions
!
!-----------------------------------------------------------------------
!
  REAL  :: wmr2td,oe,dpt
  EXTERNAL wmr2td,oe,dpt
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=10) :: nzstr
  REAL    :: latll,lonll,latur,lonur,dxkm,dykm,latdxdy,londxdy,rotate
  INTEGER :: ivar,ls,itime,itimref,iproj,idxtime, ntime

  INTEGER :: i,j,k,klev,ifile,ireturn,imid,jmid,ncmode
  REAL    :: rlnpcdf,r2d
  REAL    :: ctrx,ctry,swx,swy
  REAL    :: gamma,ex1,ex2,rln700,p00,t00
  REAL    :: prmb,tdew
  LOGICAL :: adddefs,ncexist

  REAL, ALLOCATABLE :: p_mb(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
!
  WRITE(6,'(9(/a))')                                                    &
      '###############################################################',&
      '###############################################################',&
      '####                                                       ####',&
      '####                Welcome to ARPS2NCDF                   ####',&
      '####        A program that reads in history files          ####',&
      '####  generated by ARPS and produces a netCDF format file. ####',&
      '####                                                       ####',&
      '###############################################################',&
      '###############################################################'
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  get_input_filenames reads the history file input options
!  from the namelist then open the file and read some of
!  the grid dimension information
!
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lengbf = len_trim(grdbasfn)

  CALL get_dims_from_data(hinfmt,TRIM(hisfile(1)),                      &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(3(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz

!
!-----------------------------------------------------------------------
!
!  Read the netcdf output options from the netcdf namelist options
!       ncdapnd    option to append to netcdf file
!       ncdntime   the number of valid times in the file
!       the name of the grid/base data set.
!
!-----------------------------------------------------------------------
!
  ncdapnd=0
  ipktyp=1
  nbits=16
  READ(5,ncdf_options,END=900)

  lengbf=LEN_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

!-----------------------------------------------------------------------
! output grid is 3 gridpoints less than history file grid
! in both x and y directions
!-----------------------------------------------------------------------

  ncdflvls = nprlvl + 1       ! include surface layer

  nxcdf=(nx-3)
  nycdf=(ny-3)

  planelen(1) = nxcdf
  planelen(2) = nycdf
  planelen(3) = 1
  planelen(4) = 1
  lvllen(1)   = 10
  lvllen(2)   = ncdflvls
  sfclen(1)   = 10
  sfclen(2)   = 1
  vollen(1)   = nxcdf
  vollen(2)   = nycdf
  vollen(3)   = ncdflvls
  vollen(4)   = 1
  ncdntime    = max(nhisfile,ncdntime)

!-----------------------------------------------------------------------
!
! now that you've read all the dimensions at execute time
! allocate the memory for the arrays
!
!-----------------------------------------------------------------------
  ALLOCATE(x      (nx))
  ALLOCATE(y      (ny))
  ALLOCATE(z      (nz))
  ALLOCATE(zp     (nx,ny,nz))
  ALLOCATE(zpsoil (nx,ny,nzsoil))

  ALLOCATE(uprt   (nx,ny,nz))
  ALLOCATE(vprt   (nx,ny,nz))
  ALLOCATE(wprt   (nx,ny,nz))
  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))
  ALLOCATE(qv     (nx,ny,nz))
  ALLOCATE(qscalar (nx,ny,nz,nscalar))
  ALLOCATE(qi     (nx,ny,nz))
  ALLOCATE(tke    (nx,ny,nz))
  ALLOCATE(kmh    (nx,ny,nz))
  ALLOCATE(kmv    (nx,ny,nz))
  ALLOCATE(ubar   (nx,ny,nz))
  ALLOCATE(vbar   (nx,ny,nz))
  ALLOCATE(wbar   (nx,ny,nz))
  ALLOCATE(ptbar  (nx,ny,nz))
  ALLOCATE(pbar   (nx,ny,nz))
  ALLOCATE(rhobar (nx,ny,nz))
  ALLOCATE(qvbar  (nx,ny,nz))

  ALLOCATE(soiltyp (nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp (nx,ny))
  ALLOCATE(lai    (nx,ny))
  ALLOCATE(roufns (nx,ny))
  ALLOCATE(veg    (nx,ny))

  ALLOCATE(tsoil  (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil  (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(rain (nx,ny))
  ALLOCATE(raing(nx,ny))
  ALLOCATE(rainc(nx,ny))
  ALLOCATE(prcrate(nx,ny,4))

  ALLOCATE(radfrc(nx,ny,nz))
  ALLOCATE(radsw (nx,ny))
  ALLOCATE(rnflx (nx,ny))
  ALLOCATE(radswnet(nx,ny))
  ALLOCATE(radlwin(nx,ny))

  ALLOCATE(usflx (nx,ny))
  ALLOCATE(vsflx (nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))

  ALLOCATE(e_mb   (nx,ny,nz))
  ALLOCATE(mix    (nx,ny,nz))
  ALLOCATE(esat_mb(nx,ny,nz))
  ALLOCATE(rh     (nx,ny,nz))
  ALLOCATE(t_dew  (nx,ny,nz))

  ALLOCATE(lcl(nx,ny))
  ALLOCATE(lfc(nx,ny))
  ALLOCATE(el(nx,ny))
  ALLOCATE(twdf(nx,ny))
  ALLOCATE(li(nx,ny))
  ALLOCATE(pbe(nx,ny))
  ALLOCATE(mbe(nx,ny))
  ALLOCATE(nbe(nx,ny))
  ALLOCATE(tcap(nx,ny))

  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))
  ALLOCATE(tem4(nx,ny,nz))
  ALLOCATE(tem2d1(nx,ny))
  ALLOCATE(tem2d2(nx,ny))
  ALLOCATE(tem2d3(nx,ny))
  ALLOCATE(tem2d4(nx,ny))

  ALLOCATE(wrk1(nz),wrk2(nz),wrk3(nz),wrk4(nz),wrk5(nz),wrk6(nz))
  ALLOCATE(wrk7(nz),wrk8(nz),wrk9(nz),wrk10(nz),wrk11(nz),wrk12(nz))

  ALLOCATE(xs(nx))
  ALLOCATE(ys(ny))
  ALLOCATE(zps(nx,ny,nz))
  ALLOCATE(zps_km(nx,ny,nz))

  ! J.Case, ENSCO Inc. (9/29/2004)
  ALLOCATE(radar_3d(nx,ny,nz))
  ALLOCATE(radar_comp(nx,ny))
  ALLOCATE(echo_top(nx,ny))

  ALLOCATE(p_mb(nx,ny,nz))

  ALLOCATE(heli (nx,ny))
  ALLOCATE(brn  (nx,ny))
  ALLOCATE(brnu (nx,ny))
  ALLOCATE(srlfl(nx,ny))
  ALLOCATE(srmfl(nx,ny))
  ALLOCATE(shr37(nx,ny))
  ALLOCATE(ustrm(nx,ny))
  ALLOCATE(vstrm(nx,ny))
  ALLOCATE(blcon(nx,ny))

  ALLOCATE(mf2d(nx,ny))
  ALLOCATE(coriolis(nx,ny))
  ALLOCATE(spacing(nx,ny))


  ALLOCATE(static2d(nxcdf,nycdf,nstvar))
  ALLOCATE(cdf2d(nxcdf,nycdf,1,n2dvar))
  ALLOCATE(cdf3d(nxcdf,nycdf,ncdflvls,n3dvar))

!-----------------------------------------------------------------------
! initialize a bunch of stuff
!-----------------------------------------------------------------------
  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0

  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
  qscalar=0.0
  qi     = 0.0
  tke    =0.0
  kmh    =0.0
  kmv    =0.0
  ubar   =0.0
  vbar   =0.0
  wbar   =0.0
  ptbar  =0.0
  pbar   =0.0
  rhobar =0.0
  qvbar  =0.0
  qv     =0.0

  soiltyp =0.0
  stypfrct=0.0
  vegtyp =0.0
  lai    =0.0
  roufns =0.0
  veg    =0.0

  tsoil  =0.0
  qsoil  =0.0
  wetcanp=0.0
  snowdpth=0.0

  rain =0.0
  raing=0.0
  rainc=0.0
  prcrate=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet =0.0
  radlwin =0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  e_mb    =0.0
  mix     =0.0
  esat_mb =0.0
  rh      =0.0
  t_dew   =0.0

  lcl =0.0
  lfc =0.0
  el  =0.0
  twdf=0.0
  li  =0.0
  pbe =0.0
  mbe =0.0
  nbe =0.0
  tcap=0.0

  tem1 =0.0
  tem2 =0.0
  tem3 =0.0
  tem4 =0.0
  tem2d1=0.0
  tem2d2=0.0
  tem2d3=0.0
  tem2d4=0.0

  wrk1 = 0.0
  wrk2 = 0.0
  wrk3 = 0.0
  wrk4 = 0.0
  wrk5 = 0.0
  wrk6 = 0.0
  wrk7 = 0.0
  wrk8 = 0.0
  wrk9 = 0.0
  wrk10= 0.0
  wrk11= 0.0
  wrk12= 0.0

  xs=0.0
  ys=0.0
  zps=0.0
  zps_km=0.0

  ! J.Case, ENSCO Inc. (9/29/2004)
  radar_3d=0.0
  radar_comp=0.0
  echo_top=0.0

  p_mb =0.0

  heli =0.0
  brn  =0.0
  brnu =0.0
  srlfl=0.0
  srmfl=0.0
  shr37=0.0
  ustrm=0.0
  vstrm=0.0
  blcon=0.0

  mf2d=0.0
  coriolis=0.0
  spacing=0.0

  init=.false.

  nrecord=0
  nvaltimes=0

  DO ifile=1,nhisfile
    lenfil = LEN(hisfile(ifile))
    CALL strlnth( hisfile(ifile), lenfil)
    WRITE(6,'(/a,/1x,a)')' The data set name is ',hisfile(ifile)(1:lenfil)
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    CALL dtaread(nx,ny,nz,nzsoil,nstyps,hinfmt,nch,grdbasfn(1:lengbf),  &
               lengbf,                                                  &
               hisfile(ifile)(1:lenfil),lenfil,time,                    &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt, qscalar, tke, kmh, kmv,                           &
               ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,            &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
    IF( ireturn == 0 ) THEN   ! successful read
!
!===============================================================
!
!  J.Case, ENSCO, Inc. 9/29/2004 -- added subroutine call for
!                          derived radar reflectivity products.
!
!  D. Dawson, CAPS. 05/09/2008 -modified to choose reflectivity
!                          calculation based on microphysics
!                          scheme
!
!===============================================================
!

      ! Calculate temperature (store in tem1)

      !-----------------------------------------------------------------------
      !
      !  Calculate the temperature using Poisson's formula.
      !
      !-----------------------------------------------------------------------
      !
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1

             tem1(i,j,k) = (ptbar(i,j,k)+ptprt(i,j,k)) *    &
             (((pprt(i,j,k) + pbar(i,j,k)) / p0) ** rddcp)

          END DO
        END DO
      END DO

      WRITE(*,*) ' '
      WRITE (*,*) 'Deriving radar reflectivity.'
      IF(nscalar == 2) THEN ! Warm-rain only (Kessler) scheme
        CALL reflec_wr(nx,ny,nz,rhobar,qscalar(1,1,1,P_QR),radar_3d)
      ELSE IF(nscalar == 5) THEN ! all ice schemes except MY
        CALL reflec_ferrier (nx,ny,nz,rhobar,qscalar,tem1,radar_3d)
      ELSE
        CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tem1,radar_3d)
      END IF

      DO j=1,ny-1
        DO i=1,nx-1
          radar_comp(i,j) = 0.0
          echo_top(i,j)   = 0.0
          first18         = 1

          DO k=1,nz-1
            IF (radar_3d(i,j,k) > radar_comp(i,j)) THEN
              radar_comp(i,j) = radar_3d(i,j,k)
            END IF
          END DO

          DO k=nz-1,1,-1
            IF (radar_3d(i,j,k) >= 18.0 .AND. first18 == 1) THEN
              echo_top(i,j) = zp(i,j,k)
              first18       = 0
            END IF
          END DO

        END DO
      END DO

      IF( .NOT.init ) THEN
!
!-----------------------------------------------------------------------
!
!  Establish coordinate for scalar fields.
!
!-----------------------------------------------------------------------
!
        DO i=1,nx-1
          xs(i)=0.5*(x(i)+x(i+1))
        END DO
        xs(nx)=2.*xs(nx-1)-xs(nx-2)

        DO j=1,ny-1
          ys(j)=0.5*(y(j)+y(j+1))
        END DO
        ys(ny)=2.*ys(ny-1)-ys(ny-2)

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Global Attributes
!
!-----------------------------------------------------------------------
!
        latnot(1)=trulat1
        latnot(2)=trulat2
        CALL setmapr( mapproj, 1.0 , latnot , trulon)

        CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
        dx=x(3)-x(2)
        dy=y(3)-y(2)
        dxkm=0.001*dx
        dykm=0.001*dy
        swx = ctrx - (FLOAT(nx-3)/2.) * dx
        swy = ctry - (FLOAT(ny-3)/2.) * dy
        CALL setorig( 1, swx, swy)

        CALL xytoll(1,1,xs(2),ys(2),latll,lonll)
        CALL xytoll(1,1,xs(nx-2),ys(ny-2),latur,lonur)

        IF(mapproj > 0) THEN
          iproj=kproj(mapproj)
          cproj=chproj(iproj)
        ELSE
          iproj=0
          cproj='NONE'
        END IF

        IF(mapproj == 2) THEN
          latdxdy=0.5*(trulat1+trulat2)
        ELSE
          latdxdy=trulat1
        END IF
        londxdy=ctrlon

        CALL xytoll(nx,ny,xs,ys,coriolis,tem1)
        CALL xytomf(nx,ny,xs,ys,mf2d)

        r2d = ACOS(-1.0)/180.0
        DO j=1,ny-1
          DO i=1,nx-1
            coriolis(i,j) = 2.*omega*SIN( r2d * coriolis(i,j) )
            spacing(i,j) = dx*mf2d(i,j)
          END DO
        END DO
!-----------------------------------------------------------------------
!
!  Open netCDF file
!
!-----------------------------------------------------------------------
!
        INQUIRE(FILE=cdfname,EXIST=ncexist)

!-----------------------------------------------------------------------
!       if this is an append and the netCDF files exists
!       then read the record and n_valtimes dimensions from the
!       netcdf file
!-----------------------------------------------------------------------
        IF(ncdapnd > 0 .AND. ncexist) THEN
          ncmode=NF_WRITE
          adddefs=.false.
          istatus=NF_OPEN(cdfname,ncmode,ncid)
          IF (istatus /= NF_NOERR) CALL handle_err('nf_open',istatus)

          ! i think we want to replace with a read of the
          ! record dimension to determine the record to write to
          ! and the n_valtimes dimension to determine the maximum
          ! number of forecast times that can be included in the file (tdo)

          istatus=NF_INQ_DIMID(ncid,'record',dimidrec)
          IF (istatus /= NF_NOERR) CALL handle_err('nf_inq_dimid record',istatus)
          istatus=NF_INQ_DIMLEN(ncid,dimidrec,nrecord)

          istatus=NF_INQ_DIMID(ncid,'n_valtimes',dimidnt)
          IF (istatus /= NF_NOERR) CALL handle_err('nf_inq_dimid ntime',istatus)
          istatus=NF_INQ_DIMLEN(ncid,dimidnt,nvaltimes)

          WRITE(6,'(2(a,i5))') 'current record =',nrecord,              &
                ',number of valid times =',nvaltimes,                   &
                ',number of history files =',nhisfile
          ntotaltimes=nrecord+nhisfile
          IF (ntotaltimes > nvaltimes) THEN
            WRITE(6,'(a)') 'Trying to write too many hours into netcdf file'
            STOP
          END IF

          ntwrt=nrecord
          ! i think this is the simplest change to preserve Keith Brewsters
          ! append capability while living within the AWIPS netCDF conventions

!          istatus=nf_inq_varid(ncid,'ntwrite',ntwrtid)
!          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_var_id',istatus)
!          istatus=nf_get_var_int(ncid,ntwrtid,ntwrt)
!          IF (istatus /= nf_noerr) CALL handle_err('nf_get_var_int',istatus)

        ELSE
          ncmode=NF_CLOBBER
          adddefs=.true.
          istatus=NF_CREATE(cdfname,ncmode,ncid)
          ntwrt=0
          IF (istatus /= NF_NOERR) CALL handle_err('nf_create',istatus)
        END IF
!
!-----------------------------------------------------------------------
!
!  Define dimensions
!
!-----------------------------------------------------------------------
!
        IF(ncdflvls < 10) THEN
          WRITE(nzstr,'(a,i1)') 'levels_',ncdflvls
        ELSE IF(ncdflvls < 100) THEN
          WRITE(nzstr,'(a,i2)') 'levels_',ncdflvls
        ELSE
          WRITE(nzstr,'(a,i3)') 'levels_',ncdflvls
        END IF

        IF( adddefs ) THEN

          istatus=nf_def_dim(ncid,'record',nf_unlimited,dimidrec)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim record',istatus)

          istatus=nf_def_dim(ncid,'n_valtimes',ncdntime,dimidnt)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim ncdntime',istatus)

          istatus=nf_def_dim(ncid,'data_variables',nvartot,dimidtot)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim nvartot',istatus)

          istatus=nf_def_dim(ncid,'charsPerLevel',10,dimchrlvl)
          IF (istatus /= nf_noerr)                                        &
              CALL handle_err('nf_def_dim charsPerLevel',istatus)

          istatus=nf_def_dim(ncid,'namelen',132,dimnamlen)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim namelen',istatus)

          istatus=nf_def_dim(ncid,'x',nxcdf,dimidnx)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim nx',istatus)

          istatus=nf_def_dim(ncid,'y',nycdf,dimidny)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim ny',istatus)

          istatus=nf_def_dim(ncid,'levels_1',1,dimid1)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim 1',istatus)

          istatus=nf_def_dim(ncid,nzstr,ncdflvls,dimidnz)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim levels',istatus)

          istatus=nf_def_dim(ncid,'nav',1,dimidnav)
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_dim nav',istatus)
!
!       The following arrays are used to create the dimensions of the variables
!       that are defined next.  The arrays consist of dimension ids (from above)
!       When listed by the ncdump utility, the 1st dimension (i.e. arraydim(1)) is
!       the right most dimension listed while the last dimension is the left most
!       dimension.  The order of dimensions here is set to match the AWIPS convention
!
        planedim(1)=dimidnx
        planedim(2)=dimidny
        planedim(3)=dimid1
        planedim(4)=dimidrec
        voldim(1)=dimidnx
        voldim(2)=dimidny
        voldim(3)=dimidnz
        voldim(4)=dimidrec
        lvldim(1)=dimchrlvl
        lvldim(2)=dimidnz
        invdim(1)=dimidnz
        invdim(2)=dimidnt
!
!-----------------------------------------------------------------------
!
!  Define 3-d variables and attributes
!
!-----------------------------------------------------------------------
!
        DO ivar=1,n3dvar
          vrange(1)=v3dmin(ivar)
          vrange(2)=v3dmax(ivar)

          istatus=nf_def_var(ncid,v3dnam(ivar),nf_float,4,              &
                             voldim,v3did(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var 3d',istatus)

          ls=LEN(v3dlngnam(ivar))
          CALL strlnth(v3dlngnam(ivar),ls)
          istatus=nf_put_att_text(ncid,v3did(ivar),'long_name',         &
                          ls,v3dlngnam(ivar)(1:ls))
          ls=LEN(v3dunits(ivar))
          CALL strlnth(v3dunits(ivar),ls)
          istatus=nf_put_att_text(ncid,v3did(ivar),'units',             &
                          ls,v3dunits(ivar)(1:ls))
          istatus=nf_put_att_real(ncid,v3did(ivar),'valid_range',       &
                          nf_float,2,vrange)
          istatus=nf_put_att_real(ncid,v3did(ivar),'_FillValue',        &
                          nf_float,1,-99999.)
          istatus=nf_put_att_int(ncid,v3did(ivar),'_n3D',               &
                          nf_int,1,ncdflvls)
          istatus=nf_put_att_text(ncid,v3did(ivar),'levels',            &
                          10,'MB and SFC')

          ls=LEN(v3dnam(ivar))
          CALL strlnth(v3dnam(ivar),ls)
          WRITE(v3dlvlnam,'(a,a)') v3dnam(ivar)(1:ls),'Levels'
          istatus=nf_def_var(ncid,v3dlvlnam,nf_char,2,                  &
                             lvldim,v3dlvlid(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var lvl 3d',istatus)
!
!         create the inventory variable name by appending Inventory to variable
!         then define the inventory variable
!
          WRITE(v3dinv,'(a,a)') v3dnam(ivar)(1:ls),'Inventory'
          istatus=nf_def_var(ncid,v3dinv,nf_char,2,                     &
                             invdim,v3dinvid(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var inv 3d',istatus)

        END DO
!
!-----------------------------------------------------------------------
!
!  Define 2-d variables and attributes
!
!-----------------------------------------------------------------------
!
        invdim(1)=dimid1
        invdim(2)=dimidnt
        lvldim(2)=dimid1
        DO ivar=1,n2dvar
          vrange(1)=v2dmin(ivar)
          vrange(2)=v2dmax(ivar)
          istatus=nf_def_var(ncid,v2dnam(ivar),nf_float,4,              &
                             planedim,v2did(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var 2d',istatus)
          ls=LEN(v2dlngnam(ivar))
          CALL strlnth(v2dlngnam(ivar),ls)
          istatus=nf_put_att_text(ncid,v2did(ivar),'long_name',         &
                                  ls,v2dlngnam(ivar)(1:ls))
          ls=LEN(v2dunits(ivar))
          CALL strlnth(v2dunits(ivar),ls)
          istatus=nf_put_att_text(ncid,v2did(ivar),'units',             &
                                  ls,v2dunits(ivar)(1:ls))
          istatus=nf_put_att_real(ncid,v2did(ivar),'valid_range',       &
                                  nf_float,2,vrange)
          istatus=nf_put_att_real(ncid,v2did(ivar),'_FillValue',        &
                          nf_float,1,-99999.)
          istatus=nf_put_att_int(ncid,v2did(ivar),'_n3D',               &
                          nf_int,1,0)
          IF(ivar == idxmslp) THEN
            istatus=nf_put_att_text(ncid,v2did(ivar),'levels',          &
                          3,'MSL')
          ELSE
            istatus=nf_put_att_text(ncid,v2did(ivar),'levels',          &
                          3,'SFC')
          END IF

          ls=LEN(v2dnam(ivar))
          CALL strlnth(v2dnam(ivar),ls)
          WRITE(v2dlvlnam,'(a,a)') v2dnam(ivar)(1:ls),'Levels'
          istatus=nf_def_var(ncid,v2dlvlnam,nf_char,2,                  &
                             lvldim,v2dlvlid(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var lvl 2d',istatus)
!
!         create the inventory variable name by appending Inventory to variable
!         then define the inventory variable
!
          WRITE(v2dinv,'(a,a)') v2dnam(ivar)(1:ls),'Inventory'
          istatus=nf_def_var(ncid,v2dinv,nf_char,2,                     &
                             invdim,v2dinvid(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var inv 2d',istatus)
        END DO
!
!-----------------------------------------------------------------------
!
!  Define 1-d variables
!
!-----------------------------------------------------------------------
!
        istatus=nf_def_var(ncid,'valtimeMINUSreftime',nf_int,1,         &
                           dimidnt,vtmrefid)
        IF (istatus /= nf_noerr) CALL handle_err('nf_def_var vtmref',istatus)
        istatus=nf_put_att_text(ncid,vtmrefid,'units',                  &
                                7,'seconds')

        istatus=nf_def_var(ncid,'reftime',nf_double,0,                  &
                           0,reftimid)
        IF (istatus /= nf_noerr) CALL handle_err('nf_def_var reftim',istatus)
        istatus=nf_put_att_text(ncid,reftimid,'units',                  &
                   33,'seconds since 1970-1-1 00:00:00.0')

        istatus=nf_def_var(ncid,'valtime',nf_double,1,                  &
                           dimidnt,valtimid)
        IF (istatus /= nf_noerr) CALL handle_err('nf_def_var valtime',istatus)
        istatus=nf_put_att_text(ncid,vtmrefid,'units',                  &
                   33,'seconds since 1970-1-1 00:00:00.0')

        istatus=nf_def_var(ncid,'origin',nf_char,1,                     &
                           dimnamlen,origid)
        IF (istatus /= nf_noerr) CALL handle_err('nf_def_var orig',istatus)

        istatus=nf_def_var(ncid,'model',nf_char,1,                      &
                           dimnamlen,modelid)
        IF (istatus /= nf_noerr) CALL handle_err('nf_def_var model',istatus)

!       commented out the following definition of ntwrt variable
!       probably should be a global attribute (tdo)
!        istatus=nf_def_var(ncid,'ntwrite',nf_int,0,0,ntwrtid)
!        IF (istatus /= nf_noerr) CALL handle_err('nf_def_var ntwrite',istatus)

!
!  Define static variables
!
        DO ivar=1,nstvar
          istatus=nf_def_var(ncid,vstnam(ivar),nf_float,2,              &
                             planedim,vstid(ivar))
          IF (istatus /= nf_noerr) CALL handle_err('nf_def_var static',istatus)
          ls=LEN(vstlngnam(ivar))
          CALL strlnth(vstlngnam(ivar),ls)
          istatus=nf_put_att_text(ncid,vstid(ivar),'long_name',         &
                                  ls,vstlngnam(ivar)(1:ls))
          ls=LEN(vstunits(ivar))
          CALL strlnth(vstunits(ivar),ls)
          istatus=nf_put_att_text(ncid,vstid(ivar),'units',             &
                                  ls,vstunits(ivar)(1:ls))
        END DO

        WRITE(6,'(a,a)') '  cdl date ',cdldate

        istatus=nf_put_att_text(ncid,nf_global,'cdlDate',               &
                                8,cdldate)
        ls=LEN(depictor)
        CALL strlnth(depictor,ls)
        istatus=nf_put_att_text(ncid,nf_global,'depictorName',          &
                                ls,depictor)
        istatus=nf_put_att_int(ncid,nf_global,'projIndex',              &
                                nf_int,1,iproj)
        ls=LEN(cproj)
        CALL strlnth(cproj,ls)
        istatus=nf_put_att_text(ncid,nf_global,'projName',              &
                                ls,cproj)
        istatus=nf_put_att_real(ncid,nf_global,'centralLat',            &
                                nf_float,1,ctrlat)
        istatus=nf_put_att_real(ncid,nf_global,'centralLon',            &
                                nf_float,1,ctrlon)
        rotate=ctrlon-trulon
        istatus=nf_put_att_real(ncid,nf_global,'rotation',              &
                                nf_float,1,rotate)
        istatus=nf_put_att_real(ncid,nf_global,'lat00',                 &
                                nf_float,1,latll)
        istatus=nf_put_att_real(ncid,nf_global,'lon00',                 &
                                nf_float,1,lonll)
        istatus=nf_put_att_real(ncid,nf_global,'latNxNy',               &
                                nf_float,1,latur)
        istatus=nf_put_att_real(ncid,nf_global,'lonNxNy',               &
                                nf_float,1,lonur)
        istatus=nf_put_att_real(ncid,nf_global,'dxKm',                  &
                                nf_float,1,dxkm)
        istatus=nf_put_att_real(ncid,nf_global,'dyKm',                  &
                                nf_float,1,dykm)
        istatus=nf_put_att_real(ncid,nf_global,'latDxDy',               &
                                nf_float,1,latdxdy)
        istatus=nf_put_att_real(ncid,nf_global,'lonDxDy',               &
                                nf_float,1,londxdy)

!
!-----------------------------------------------------------------------
!
!  End define mode
!
!-----------------------------------------------------------------------
!
          istatus=nf_enddef(ncid)
          IF (istatus /= nf_noerr) CALL handle_err('nf_enddef',istatus)

        ELSE  ! get variable ids from existing dataset

!
!  Get dimension ids
!  The inquery for record and n_valtimes is redundant now, but ok so won't touch
!
          istatus=nf_inq_dimid(ncid,'record',dimidrec)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid record',istatus)

          istatus=nf_inq_dimid(ncid,'n_valtimes',dimidnt)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid ntime',istatus)

          istatus=nf_inq_dimid(ncid,'data_variables',dimidtot)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid nvartot',istatus)

          istatus=nf_inq_dimid(ncid,'charsPerLevel',dimchrlvl)
          IF (istatus /= nf_noerr)                                        &
            CALL handle_err('nf_inq_dimid charsPerLevel',istatus)

          istatus=nf_inq_dimid(ncid,'namelen',dimnamlen)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid namelen',istatus)

          istatus=nf_inq_dimid(ncid,'x',dimidnx)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid nx',istatus)

          istatus=nf_inq_dimid(ncid,'y',dimidny)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid ny',istatus)

          istatus=nf_inq_dimid(ncid,'levels_1',dimid1)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid 1',istatus)

          istatus=nf_inq_dimid(ncid,nzstr,dimidnz)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid levels',istatus)

          istatus=nf_inq_dimid(ncid,'nav',dimidnav)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_dimid nav',istatus)
!
!  Get variable ids
!
          DO ivar=1,n3dvar
            istatus=nf_inq_varid(ncid,v3dnam(ivar),v3did(ivar))
            IF (istatus /= nf_noerr) CALL handle_err('nf_inq_varid 3d',istatus)
            ls=LEN(v3dnam(ivar))
            CALL strlnth(v3dnam(ivar),ls)
            WRITE(v3dinv,'(a,a)') v3dnam(ivar)(1:ls),'Inventory'
            istatus=nf_inq_varid(ncid,v3dinv,v3dinvid(ivar))
            IF (istatus /= nf_noerr) CALL handle_err('nf_inq_varid 3dinv',istatus)
          END DO
          DO ivar=1,n2dvar
            istatus=nf_inq_varid(ncid,v2dnam(ivar),v2did(ivar))
            IF (istatus /= nf_noerr) CALL handle_err('nf_inq_varid 2d',istatus)
            ls=LEN(v2dnam(ivar))
            CALL strlnth(v2dnam(ivar),ls)
            WRITE(v2dinv,'(a,a)') v2dnam(ivar)(1:ls),'Inventory'
            istatus=nf_inq_varid(ncid,v2dinv,v2dinvid(ivar))
            IF (istatus /= nf_noerr) CALL handle_err('nf_inq_varid 2dinv',istatus)
          END DO

          istatus=nf_inq_varid(ncid,'valtimeMINUSreftime',vtmrefid)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_var_id',istatus)

! we're going to have trouble here need to replace with a read of the
! number of the number of records in the dimensions

          istatus=nf_get_vara_int(ncid,vtmrefid,1,ntwrt,vtmref)
          IF (istatus /= nf_noerr) CALL handle_err('nf_get_vara_int',istatus)

          istatus=nf_inq_varid(ncid,'reftime',reftimid)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_var_id',istatus)
          istatus=nf_get_vara_double(ncid,reftimid,1,1,reftime)
          IF (istatus /= nf_noerr) CALL handle_err('nf_get_vara_double',istatus)

          istatus=nf_inq_varid(ncid,'valtime',valtimid)
          IF (istatus /= nf_noerr) CALL handle_err('nf_inq_var_id',istatus)
          istatus=nf_get_vara_double(ncid,valtimid,1,ntwrt,valtime)
          IF (istatus /= nf_noerr) CALL handle_err('nf_get_vara_double',istatus)

        END IF
!
!-----------------------------------------------------------------------
!
!  Set inventory arrays.
!  inventory(k:k) appends the right number of 1's to the string
!
!-----------------------------------------------------------------------
!
        DO k=1,nprlvl
          inventory(k:k)='1'
          IF(iprlvl(k) > 999) THEN
            WRITE(v3dlvl(k),'(a,i4)') 'MB ',iprlvl(k)
          ELSE IF(iprlvl(k) > 99) THEN
            WRITE(v3dlvl(k),'(a,i3)') 'MB ',iprlvl(k)
          ELSE
            WRITE(v3dlvl(k),'(a,i2)') 'MB ',iprlvl(k)
          END IF
        END DO

        ! add surface to 3d level definition (tdo)
        k=ncdflvls
        inventory(k:k)='1'
        WRITE(v3dlvl(k),'(a)') 'SFC'

        !
        !       added a blank inventory (tdo)
        !
        inventory_sfc='1'
        inventory_blank=''

        CALL ctim2abss(year,month,day,hour,minute,second,itimref)
        itimref=itimref-itim1970
!
!-----------------------------------------------------------------------
!
!    End of first-read intializations
!
!-----------------------------------------------------------------------
!
        init=.true.
      END IF
!
!-----------------------------------------------------------------------
!
!    Begin ARPS data conversions
!
!-----------------------------------------------------------------------
!
      itime=itimref+nint(time)
      idxtime=ifile+ntwrt
      valtime(idxtime)=FLOAT(itime)
      vtmref(idxtime)=nint(time)
      reftime(idxtime)=FLOAT(itimref)

      DO i=1,nx
        DO j=1,ny
          rain(i,j)=raing(i,j)+rainc(i,j)
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k)=pprt(i,j,k)+pbar(i,j,k)
            ptprt(i,j,k)=ptprt(i,j,k)+ptbar(i,j,k)
            qvprt(i,j,k)=qvprt(i,j,k)+qvbar(i,j,k)
            tem1(i,j,k)=0.5*(uprt(i,j,k)+ubar(i,j,k)+                   &
                        uprt(i+1,j,k)+ubar(i+1,j,k))
            tem2(i,j,k)=0.5*(vprt(i,j,k)+vbar(i,j,k)+                   &
                        vprt(i,j+1,k)+vbar(i,j+1,k))
            tem3(i,j,k)=0.5*(wprt(i,j,k)+wbar(i,j,k)+                   &
                        wprt(i,j,k+1)+wbar(i,j,k+1))
!            qi(i,j,k)=qi(i,j,k)+qs(i,j,k)+qh(i,j,k)
            qi(i,j,k) = 0.0
            IF(P_QI > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QI)
            END IF
            IF(P_QS > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QS)
            END IF
            IF(P_QG > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QG)
            END IF
            IF(P_QH > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QH)
            END IF
          END DO
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!  Swap wind data back into wind arrays
!
!-----------------------------------------------------------------------
!
      uprt=0.
      vprt=0.
      wprt=0.
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            uprt(i,j,k)=tem1(i,j,k)
            vprt(i,j,k)=tem2(i,j,k)
            wprt(i,j,k)=tem3(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Put temperature into tem2 and -ln(p) into tem3
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=ptprt(i,j,k)*                                   &
                          (pprt(i,j,k)/p0)**rddcp
            tem3(i,j,k)=-ALOG(pprt(i,j,k))
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            zps_km(i,j,k)=zps(i,j,k)/1000.0
            p_mb(i,j,k)=0.01*pprt(i,j,k)
            mix(i,j,k)=1000.0*(qvprt(i,j,k)/(1.-qvprt(i,j,k)))
            e_mb(i,j,k)=qvprt(i,j,k)*p_mb(i,j,k)/(qvprt(i,j,k)+         &
                        287.0/461.0)
            esat_mb(i,j,k)=6.1078*EXP((2.500780E6/461.0)*((1.0/273.15)- &
                          (1.0/tem2(i,j,k))))
            t_dew(i,j,k)=wmr2td(p_mb(i,j,k),mix(i,j,k))                 &
                       + 273.15
            rh(i,j,k)=100.0*(e_mb(i,j,k)/esat_mb(i,j,k))
            rh(i,j,k)=MAX(0.,rh(i,j,k))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Calculate stability indices.
!  Use level k=2 as the "surface".
!
!-----------------------------------------------------------------------
!
      imid=(nx/2)+1
      jmid=(ny/2)+1
      CALL arps_be(nx,ny,nz,                                            &
           pprt,zps_km,tem2,qvprt,                                      &
           lcl,lfc,el,twdf,li,pbe,mbe,nbe,tcap,                         &
           wrk1,wrk2,wrk3,wrk4,wrk5,wrk6,                               &
           wrk7,wrk8,wrk9,wrk10,wrk11,wrk12,tem2d1)

      PRINT *, ' Sample stability values: '
      PRINT *, '   lcl,lfc : ',lcl(imid,jmid),lfc(imid,jmid)
      PRINT *, '   el, twdf: ',el(imid,jmid),twdf(imid,jmid)
      PRINT *, '   li, pbe : ',li(imid,jmid),pbe(imid,jmid)
      PRINT *, '   mbe, nbe: ',mbe(imid,jmid),nbe(imid,jmid)
      PRINT *, '   tcap    : ',tcap(imid,jmid)

      CALL calcshr(nx,ny,nz,x,y,zp,mf2d,                                &
          pprt,tem2,uprt,vprt,pbe,                                      &
          shr37,ustrm,vstrm,srlfl,srmfl,heli,brn,brnu,blcon,            &
          tem2d1,tem2d2,tem2d3,tem2d4,tem4)

      PRINT *, ' Sample shear values: '
      PRINT *, ' shr37,ustrm,vstrm: ',                                  &
          shr37(imid,jmid),ustrm(imid,jmid),vstrm(imid,jmid)
      PRINT *, ' srlfl,srmfl,heli: ',                                   &
          srlfl(imid,jmid),srmfl(imid,jmid),heli(imid,jmid)
      PRINT *, ' brn,brnu,blcon: ',                                     &
          brn(imid,jmid),brnu(imid,jmid),blcon(imid,jmid)
!
!  Store k=2 theta-e and dewpt in tem1,
!  level 1 and 2 respectively.
!
      DO j=1,ny-1
        DO i=1,nx-1
          prmb=0.01*pprt(i,j,2)
          tdew=wmr2td(prmb,(1000.*qvprt(i,j,2)))
          tem1(i,j,1)=oe((tem2(i,j,2)-273.15),tdew,prmb) + 273.15
          tem1(i,j,2)=tdew+273.15
          tem1(i,j,3)=prmb
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Output near-sfc data.
!  Simularity theory or something could be applied here
!  to make these data be valid at sfc instrument height,
!  for now we output level 2.  This is approximately 10m in
!  our configuration
!
!-----------------------------------------------------------------------
!

      ! addition of near surface gph (tdo)
      PRINT *, ' Storing near-sfc gh'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, zps(1,1,2),                      &
                    cdf3d(1,1,ncdflvls,idxgh) )

      ! addition of near surface rh (tdo)
      PRINT *, ' Storing near-sfc rh'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, rh(1,1,2),                       &
                    cdf3d(1,1,ncdflvls,idxrh) )

      PRINT *, ' Storing near-sfc pressure'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,3),                     &
                    cdf2d(1,1,1,idxp) )

      PRINT *,' Storing grid-scale rainfall'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, raing,                           &
                    cdf2d(1,1,1,idxlgsp) )

      PRINT *,' Storing convective rainfall'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, rainc,                           &
                    cdf2d(1,1,1,idxcp) )

      PRINT *,' Storing total accumulated rainfall'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, rain,                            &
                    cdf2d(1,1,1,idxtp) )

     ! delete sfu from netcdf file
     !PRINT *, ' Storing near-sfc u wind'
     !CALL cdf_fill(nx,ny,nxcdf,nycdf, uprt(1,1,2),                     &
     !                    cdf2d(1,1,1,idxsfu) )

      ! add write of u wind into 3d array (tdo)
      PRINT *, ' Storing near-sfc u wind into 3 d array'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, uprt(1,1,2),                     &
                    cdf3d(1,1,ncdflvls,idxuw) )

!      delete sfv from netcdf file
!      PRINT *, ' Storing near-sfc v wind'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, vprt(1,1,2),                     &
!                    cdf2d(1,1,1,idxsfv) )

!     add write of v wind into 3d array (tdo)
      PRINT *, ' Storing near-sfc v wind into 3 d array'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, vprt(1,1,2),                     &
                    cdf3d(1,1,ncdflvls,idxvw) )

!     add write of w wind into 3d array (tdo)
      PRINT *, ' Storing near-sfc w wind into 3 d array'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, wprt(1,1,2),                     &
                    cdf3d(1,1,ncdflvls,idxww) )

!      delete sft from netcdf file
!      PRINT *, ' Storing near-sfc temperature'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, tem2(1,1,2),                     &
!                    cdf2d(1,1,1,idxsft) )

!     add the storage of surface data into the 3d array (tdo)
      PRINT *, ' Storing near-sfc temperature into 3d array'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, tem2(1,1,2),                     &
                      cdf3d(1,1,ncdflvls,idxt) )

      PRINT *, ' Storing near-sfc dew point temperature'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,2),                     &
                    cdf2d(1,1,1,idxdpt) )

      PRINT *,' Storing near-sfc specific humidity'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, qvprt(1,1,2),                    &
                    cdf2d(1,1,1,idxsh) )

      PRINT *, ' Storing near-sfc theta-e'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                     &
                    cdf2d(1,1,1,idxept) )

      PRINT *, ' Storing near-sfc LI'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, li,                              &
                    cdf2d(1,1,1,idxsli) )

      PRINT *, ' Storing near-sfc CAPE'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, pbe,                             &
                    cdf2d(1,1,1,idxcape) )

      PRINT *, ' Storing near-sfc CIN'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, nbe,                             &
                    cdf2d(1,1,1,idxcin) )

!      PRINT *, ' Storing near-sfc LCL'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, lcl,                             &
!                    cdf2d(1,1,1,idxlcl) )

!      PRINT *, ' Storing near-sfc LFC'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, lfc,                             &
!                    cdf2d(1,1,1,idxlfc) )

!      PRINT *, ' Storing near-sfc EL'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, el,                              &
!                    cdf2d(1,1,1,idxel) )

!      PRINT *, ' Storing wet bulb temp diff'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, twdf,                            &
!                    cdf2d(1,1,1,idxtwdf) )

!      PRINT *, ' Storing Cap Strength'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, tcap,                            &
!                    cdf2d(1,1,1,idxtcap) )

!      PRINT *, ' Storing 3-7km Shear'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, shr37,                           &
!                    cdf2d(1,1,1,idxshr37) )

      PRINT *, ' Storing u-storm'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, ustrm,                           &
                    cdf2d(1,1,1,idxustrm) )

      PRINT *, ' Storing v-storm'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, vstrm,                           &
                    cdf2d(1,1,1,idxvstrm) )

!      PRINT *, ' Storing low-level SR flow'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, srlfl,                           &
!                    cdf2d(1,1,1,idxsrlfl) )

!      PRINT *, ' Storing mid-level SR flow'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, srmfl,                           &
!                    cdf2d(1,1,1,idxsrmfl) )

      PRINT *, ' Storing helicity'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, heli,                            &
                    cdf2d(1,1,1,idxheli) )

!      PRINT *, ' Storing Bulk Richardson Number'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, brn,                             &
!                    cdf2d(1,1,1,idxbrn) )

!      PRINT *, ' Storing BRN Shear'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, brnu,                            &
!                    cdf2d(1,1,1,idxbrnu) )

!      PRINT *, ' Storing Boundary Layer Convergence'
!      CALL cdf_fill(nx,ny,nxcdf,nycdf, blcon,                           &
!                    cdf2d(1,1,1,idxblcon) )
!

! J.Case, ENSOC Inc. (9/29/2004) -- Added 2D radar variables to output.

      PRINT *, ' Storing Column Max Reflectivity'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, radar_comp,                      &
                    cdf2d(1,1,1,idxcxr) )

      PRINT *, ' Storing Maximum Radar Echo Tops'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, echo_top,                        &
                    cdf2d(1,1,1,idxmret) )

!-----------------------------------------------------------------------
!
!  Calculate constant pressure level data
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Put temperature into tem2 and -ln(p) into tem3.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=ptprt(i,j,k)*                                   &
                        (pprt(i,j,k)/p0)**rddcp
            tem3(i,j,k)=-ALOG(pprt(i,j,k))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!    Calculate temperature (K) at ARPS grid points
!   and 700mb pressure level
!
!-----------------------------------------------------------------------
!
      rln700=-ALOG(70000.0)
      CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                       &
                      tem2,tem3,rln700,tem1)
!
!-----------------------------------------------------------------------
!
!    Calculate sea level pressure (mb)
!    Reduction method: Benjamin and Miller: 1990, MWR, vol.118, No.10,
!                   Page: 2100-2101
!
!-----------------------------------------------------------------------
!


!      gamma=.0065      ! std lapse rate per meter
!      ex2=5.2558774
!      DO 355 j=1,ny-1
!      DO 355 i=1,nx-1
!        p00 = 0.01*(pprt(i,j,2))
!        tem1(i,j,1)=p00*
!    :            (1.0+gamma*zps(i,j,2)/tem1(i,j,1))**ex2
! 355     CONTINUE

      gamma=.0065      ! std lapse rate per meter
      ex1=0.1903643           ! R*gamma/g
      ex2=5.2558774
      DO j=1,ny-1
        DO i=1,nx-1
          p00 = 0.01*(pprt(i,j,2))
          IF (p00 <= 700.0) THEN
            t00 = tem2(i,j,2)
          ELSE
            t00 = tem1(i,j,1)*(p00/700.0)**ex1
          END IF
          tem1(i,j,1)=p00*(1.0+gamma*zps(i,j,2)/t00)**ex2
        END DO
      END DO
      PRINT *, ' Storing MSL Pressure '
      CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1,                            &
                    cdf2d(1,1,1,idxmslp) )
!
!-----------------------------------------------------------------------
!
!  Store terrain data.
!
!-----------------------------------------------------------------------
!
      PRINT *, ' Storing terrain data.'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, zp(1,1,2),                       &
                    static2d(1,1,idxtop) )
      PRINT *, ' Storing Coriolis data'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, coriolis,                        &
                    static2d(1,1,idxcor) )
      PRINT *, ' Storing spacing data.'
      CALL cdf_fill(nx,ny,nxcdf,nycdf, spacing,                         &
                    static2d(1,1,idxspa) )
!
!-----------------------------------------------------------------------
!
!    Calculate atmospheric state variables.
!
!-----------------------------------------------------------------------
!
      DO klev=1,nprlvl
        rlnpcdf=-ALOG(100.*FLOAT(iprlvl(klev)))
        PRINT *, ' Storing netCDF data at pr= ',iprlvl(klev)
!
!        Store gh
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     zps,tem3,rlnpcdf,tem1)
        CALL extrph(nx,ny,nz,zps,tem2,pprt,                             &
                    iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1,                          &
                      cdf3d(1,1,klev,idxgh) )
!
!        Store temperature (tem2)
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                      tem2,tem3,rlnpcdf,tem1)
        CALL extrpt(nx,ny,nz,tem2,pprt,zps,                             &
                    iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1,                          &
                      cdf3d(1,1,klev,idxt) )
!
!        Store relative humidity
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     rh,tem3,rlnpcdf,tem1)
        CALL extrpq(nx,ny,nz,rh,pprt,iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1,                          &
                      cdf3d(1,1,klev,idxrh) )
!
!        Store u and v wind components
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     uprt,tem3,rlnpcdf,tem1(1,1,1))
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     vprt,tem3,rlnpcdf,tem1(1,1,2))
        CALL extrpuv(nx,ny,nz,uprt,vprt,pprt,zps,                       &
                     iprlvl(klev),tem1(1,1,1),tem1(1,1,2))
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                   &
                      cdf3d(1,1,klev,idxuw) )
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,2),                   &
                      cdf3d(1,1,klev,idxvw) )
!
!        Store vertical velocity
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     wprt,tem3,rlnpcdf,tem1)
        CALL extrpq(nx,ny,nz,wprt,pprt,iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                   &
                      cdf3d(1,1,klev,idxww) )

!
!        Store cloud water
!
        IF (P_QC > 0) THEN
          CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                       qscalar(1,1,1,P_QC),tem3,rlnpcdf,tem1)
          CALL extrpq(nx,ny,nz,qscalar(1,1,1,P_QC),pprt,iprlvl(klev),tem1)
          CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                   &
                        cdf3d(1,1,klev,idxcw) )
        END IF
!
!        Store cloud ice
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     qi(1,1,1),tem3,rlnpcdf,tem1)
        CALL extrpq(nx,ny,nz,qi(1,1,1),pprt,iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                   &
                      cdf3d(1,1,klev,idxci) )
!
!        Store turbulent kinetic energy
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     tke,tem3,rlnpcdf,tem1)
        CALL extrpq(nx,ny,nz,tke,pprt,iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                   &
                      cdf3d(1,1,klev,idxtk) )
!
! J.Case, ENSCO Inc. (9/29/2004)
!        Store 3D radar reflectivity
!
        CALL v2dinta(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,                     &
                     radar_3d,tem3,rlnpcdf,tem1)
        CALL extrpq(nx,ny,nz,radar_3d,pprt,iprlvl(klev),tem1)
        CALL cdf_fill(nx,ny,nxcdf,nycdf, tem1(1,1,1),                   &
                      cdf3d(1,1,klev,idxrr) )
      END DO
    END IF   ! Good return from read

    invlen(1)=ncdflvls
    invlen(2)=1
    vecistr(2)=idxtime
    volistart(4)=idxtime
    DO ivar=1,n3dvar
      WRITE(6,'(a,i4,2x,a)') '  Writing 3d ivar:',ivar,                 &
                              v3dlngnam(ivar)

      IF(idxtime == 1) THEN
        istatus=nf_put_vara_text(ncid,v3dlvlid(ivar),                   &
                                 planeistr,lvllen,v3dlvl)
        IF (istatus /= nf_noerr)                                        &
          CALL handle_err('nf_put_vara_text 3d lvl',istatus)
      END IF

      istatus=nf_put_vara_text(ncid,v3dinvid(ivar),                     &
                               vecistr,invlen,inventory)
      IF (istatus /= nf_noerr)                                          &
        CALL handle_err('nf_put_vara_text 3d inv',istatus)

      istatus=nf_put_vara_real(ncid,v3did(ivar),                        &
                   volistart,vollen,cdf3d(1,1,1,ivar))
      IF (istatus /= nf_noerr)                                          &
          CALL handle_err('nf_put_vara_real 3d var',istatus)
    END DO
!
    invlen(1)=1
    invlen(2)=1
    planeistr(4)=idxtime
    DO ivar=1,n2dvar
      WRITE(6,'(a,i4,2x,a)') '  Writing 2d ivar:',ivar,                 &
                                v2dlngnam(ivar)
      IF(adddefs) THEN
        IF(idxtime == 1) THEN
          IF(ivar == idxmslp) THEN
            v2dlvl='MSL'
          ELSE
            v2dlvl='SFC'
          END IF
          istatus=nf_put_vara_text(ncid,v2dlvlid(ivar),                 &
                                 vecistr,sfclen,v2dlvl)
          IF (istatus /= nf_noerr)                                      &
              CALL handle_err('nf_put_vara_text 2d lvl',istatus)
        END IF
      END IF

      istatus=nf_put_vara_text(ncid,v2dinvid(ivar),                     &
                                 vecistr,invlen,inventory_sfc)
      IF (istatus /= nf_noerr)                                          &
            CALL handle_err('nf_put_vara_text 2d inv',istatus)

      istatus=nf_put_vara_real(ncid,v2did(ivar),                        &
                     planeistr,planelen,cdf2d(1,1,1,ivar))
      IF (istatus /= nf_noerr) CALL handle_err('nf_put_vara_real 2d',istatus)
    END DO

  END DO

  IF(adddefs) THEN
!
!  Write static variables
!
    DO ivar=1,nstvar
      WRITE(6,'(a,i4,2x,a)') '  Writing st ivar:',ivar,                 &
                                vstlngnam(ivar)
      istatus=nf_put_vara_real(ncid,vstid(ivar),                        &
                       planeistr,planelen,static2d(1,1,ivar))
      IF (istatus /= nf_noerr)  CALL handle_err('nf_put_vara_real static',istatus)
    END DO

    ls=LEN(origin)
    CALL strlnth(origin,ls)
    istatus=nf_put_vara_text(ncid,origid,1,ls,origin)
    ls=LEN(model)
    CALL strlnth(model,ls)
    istatus=nf_put_vara_text(ncid,modelid,1,ls,model)

!
!   Write 1-d variables if idxtime is 1 (i.e. the first history file)
!
    IF (idxtime == 1) THEN

!
!     set up the first element of the time arrays
!     if more than 1 more valid time will be written to the file
!       foreach time
!         set valtime, vtmref, and reftime using ncdtinterval
!         write out the inventory arrays with blanks
!         AWIPS preprocessor must handle updating inventory arrays if appending
!
      valtime(idxtime)=FLOAT(itime)
      vtmref(idxtime)=nint(time)
      reftime(idxtime)=FLOAT(itimref)

      IF (ncdntime>idxtime) THEN
        DO ntime=idxtime+1,ncdntime
          valtime(ntime)=valtime(idxtime)+FLOAT((ntime-1)*ncdtinterval)
          vtmref(ntime)=vtmref(idxtime)+(ntime-1)*ncdtinterval
          reftime(ntime)=reftime(idxtime)
!
!         write the 3d var inventory with the fill character
!         start array: vecistr(1)=1 vecistr(2)=ntime
!         count array: invlen(1)=1 (1 character) invlen(2)=1
!
          vecistr(1)=1
          vecistr(2)=ntime
          invlen(1)=1
          invlen(2)=1
          DO ivar=1,n3dvar
            istatus=nf_put_vara_text(ncid,v3dinvid(ivar),                     &
                                  vecistr,invlen,nf_fill_char)
            IF (istatus /= nf_noerr)                                          &
              CALL handle_err('nf_put_vara_text 3d inv',istatus)
          END DO
!
!         write the 2d var inventory with the fill character
!         start array: vecistr(1)=1 vecistr(2)=ntime
!         count array: invlen(1)=1 invlen(2)=1
          DO ivar=1,n2dvar
            istatus=nf_put_vara_text(ncid,v2dinvid(ivar),                     &
                                  vecistr,invlen,nf_fill_char)
            IF (istatus /= nf_noerr)                                          &
              CALL handle_err('nf_put_vara_text 2d inv',istatus)
          END DO
!
        END DO
      END IF


      istatus=nf_put_vara_int(ncid,vtmrefid,                                &
                              1,ncdntime,vtmref)
!      istatus=nf_put_vara_double(ncid,reftimid,                             &
!                               1,ncdntime,reftime)
      istatus=nf_put_vara_double(ncid,reftimid,                             &
                               1,1,reftime)
      istatus=nf_put_vara_double(ncid,valtimid,                             &
                               1,ncdntime,valtime)

!    ELSE
!      WRITE(6,'(a,i4)') '  Writing st ivar:',ivar,                     &
    END IF

  END IF
! end idxtime=1

  ntwrt=ntwrt+nhisfile

!
!  Close file
!
  istatus=nf_close(ncid)
  IF (istatus /= nf_noerr) CALL handle_err('nf_close',istatus)
  STOP

  900 CONTINUE
  WRITE(6,'(a)') ' Error reading input data'
  950 CONTINUE
  WRITE(6,'(a)') ' Error setting up netCDF file'
  STOP
END PROGRAM arps2ncdf

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE EXTRPH                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extrph(nx,ny,nz,zps,t,pr,iprlvl,hgtcdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extrapolate height by using a std atmos lapse rate
!  below the last physical level.  Above the domain,
!  assume a constant temperature above 300 mb, otherwise
!  use the std atmos lapse rate.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster, from arps2gem
!  2/19/99
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: zps(nx,ny,nz)
  REAL,    INTENT(IN)  :: t(nx,ny,nz)
  REAL,    INTENT(IN)  :: pr(nx,ny,nz)
  INTEGER, INTENT(IN)  :: iprlvl
  REAL,    INTENT(OUT) :: hgtcdf(nx,ny)

  INCLUDE 'phycst.inc'

  REAL, PARAMETER :: gamma = 0.0065     ! degrees/m  lapse rate
  REAL, PARAMETER :: rddg  = (rd/g), const = (rd*gamma/g)

  INTEGER :: i,j
  REAL    :: prcdf

  prcdf=100.*FLOAT(iprlvl)
  DO j=1,ny-1
    DO i=1,nx-1
      IF(prcdf < pr(i,j,nz-1)) THEN
        IF(pr(i,j,nz-1) <= 30000.) THEN
          hgtcdf(i,j)=zps(i,j,nz-1) +                                   &
              rddg*t(i,j,nz-1)*ALOG(pr(i,j,nz-1)/prcdf)
        ELSE
          hgtcdf(i,j)=zps(i,j,nz-1) + (t(i,j,nz-1)/gamma)*              &
                      (1.-(prcdf/pr(i,j,nz-1))**const)
        END IF
      ELSE IF(prcdf >= pr(i,j,2)) THEN
        hgtcdf(i,j)=zps(i,j,2) + (t(i,j,2)/gamma)*                      &
               (1.-(prcdf/pr(i,j,2))**const)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE extrph
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE EXTRPT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extrpt(nx,ny,nz,t,pr,zps,iprlvl,tcdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extrapolate temperature by using a std atmos lapse rate
!  below the last physical level.  Above the domain,
!  assume a constant temperature above 300 mb, otherwise
!  use the std atmos lapse rate.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster, from arps2gem
!  2/19/99
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: t(nx,ny,nz)
  REAL,    INTENT(IN)  :: pr(nx,ny,nz)
  REAL,    INTENT(IN)  :: zps(nx,ny,nz)
  INTEGER, INTENT(IN)  :: iprlvl
  REAL,    INTENT(OUT) :: tcdf(nx,ny)

  INCLUDE 'phycst.inc'

  REAL, PARAMETER :: gamma = 0.0065           ! degrees/m  lapse rate
  REAL, PARAMETER :: const = (rd*gamma/g)

  INTEGER :: i,j
  REAL    :: prcdf

  prcdf=100.*FLOAT(iprlvl)
  DO j=1,ny-1
    DO i=1,nx-1
      IF(prcdf <= pr(i,j,nz-1)) THEN
        IF(pr(i,j,nz-1) <= 30000.) THEN
          tcdf(i,j)=t(i,j,nz-1)
        ELSE
          tcdf(i,j)=t(i,j,nz-1)*((prcdf/pr(i,j,nz-1))**const)
        END IF
      ELSE IF(prcdf >= pr(i,j,2)) THEN
        tcdf(i,j)=t(i,j,2)*((prcdf/pr(i,j,2))**const)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE extrpt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE EXTRPQ                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extrpq(nx,ny,nz,q,pr,iprlvl,qcdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign a zero-gradient value to scalars below ground.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster, from arps2gem
!  2/19/99
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: q(nx,ny,nz)
  REAL,    INTENT(IN)  :: pr(nx,ny,nz)
  INTEGER, INTENT(IN)  :: iprlvl
  REAL,    INTENT(OUT) :: qcdf(nx,ny)

  INTEGER :: i,j
  REAL    :: prcdf

  prcdf=100.*FLOAT(iprlvl)
  DO j=1,ny-1
    DO i=1,nx-1
      IF(prcdf < pr(i,j,nz-1)) THEN
        qcdf(i,j)=q(i,j,nz-1)
      ELSE IF(prcdf > pr(i,j,2)) THEN
        qcdf(i,j)=q(i,j,2)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE extrpq
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE EXTRPUV                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extrpuv(nx,ny,nz,us,vs,pr,zps,iprlvl,ucdf,vcdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign a constant value to sub-terrainian winds
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster, from arps2gem
!  2/19/99
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: us(nx,ny,nz)
  REAL,    INTENT(IN)  :: vs(nx,ny,nz)
  REAL,    INTENT(IN)  :: pr(nx,ny,nz)
  REAL,    INTENT(IN)  :: zps(nx,ny,nz)
  INTEGER, INTENT(IN)  :: iprlvl
  REAL,    INTENT(OUT) :: ucdf(nx,ny)
  REAL,    INTENT(OUT) :: vcdf(nx,ny)

  INTEGER :: i,j
  REAL    :: prcdf

  prcdf=100.*FLOAT(iprlvl)
  DO j=1,ny-1
    DO i=1,nx-1
      IF(prcdf < pr(i,j,nz-1)) THEN
        ucdf(i,j)=us(i,j,nz-1)
        vcdf(i,j)=vs(i,j,nz-1)
      ELSE IF(prcdf > pr(i,j,2)) THEN
        ucdf(i,j)=us(i,j,2)
        vcdf(i,j)=vs(i,j,2)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE extrpuv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HANDLE_ERR                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE handle_err(string,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Report netCDF error in plain English.  Exit.
!
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  2/19/99
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER,          INTENT(IN) :: istatus

  CHARACTER (LEN=80) :: nf_strerror
  WRITE(6,'(a,a,/a)') 'netCDF error: ',string, nf_strerror(istatus)

  STOP
END SUBROUTINE handle_err
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE CDF_FILL                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cdf_fill(nx,ny,nxcdf,nycdf,arpsvar,cdfvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fill netCDF output array which is the physical domain of the
!  ARPS scalar grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2/19/99
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny
  INTEGER, INTENT(IN)  :: nxcdf,nycdf
  REAL,    INTENT(IN)  :: arpsvar(nx,ny)
  REAL,    INTENT(OUT) :: cdfvar(nxcdf,nycdf)

  INTEGER :: i,j

  DO j=2,ny-2
    DO i=2,nx-2
      cdfvar(i-1,j-1)=arpsvar(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE cdf_fill
