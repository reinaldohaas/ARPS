!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM EXT2ARPS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM ext2arps
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Converts external forecast files coordinates and variables to those
!  used in ARPS and writes the information in a history dump file.
!
!  The actual reading of the external files is done by subroutine
!  rdextfil.  In many cases, it will only be necessary to make changes
!  in extdims.inc and rdextfil to customize the processing for the
!  user's specific data source.
!
!  This program could become more general in the future, but
!  for now we assume incoming data are supplied (or converted to)
!  pressure(Pa), height(m), temperature(K), specific humidty (kg/kg),
!  and u,v (m/s).  Units should be converted from their original
!  forms to these in rdextfil.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  Sept, 1994  after MAPS2ARPS
!
!  MODIFICATION HISTORY:
!  Oct, 1994 (KB)
!  Changed input times to be more flexible.
!
!  17 Jan, 1995 (KB)
!  Modified processing so that mean sounding would be for constant
!  pressure surfaces rather than constant ARPS levels.
!
!  23 Jan, 1995 (KB)
!  Added temporary array avgzs_ext, replacing the use of tem1 which
!  would not work if number of ext levels is greater than ARPS.
!
!  9 August, 1995 (KB)
!  Modified creation of mean soundings and storage of external
!  pressure arrays to allow for incoming data on other than
!  pressure vertical coordinates.  Added switch for buoyancy
!  balancing of vertical pressure gradient.
!
!  26 April, 1996 (KB)
!  Version 2.0, including:  Corrects divergence error
!  using a function either linear in k-level, physical height or
!  potential temperature.  Also, moisture is interpolated and
!  extrapolated using RHstar [ RHstar=sqrt(1. - RH) ], which should
!  be more linear in height, rather than qv.   Extrapolation above
!  and below available model data are constant-in-height perturbation
!  fields for all variables except assigned a standard lapse rate near
!  the ground and dT/dz=0 at the top.
!
!  30 April 1996 (KB)
!  Final preparation of version 2.0, added O'Brien option and
!  order of interpolation option to input namelist.
!
!  8 November 1996 (KB)
!  Version 2.1
!  Changed the polynomial interpolation scheme for better
!  accuracy and eliminated chance for a discontinuity at external
!  grid cell boundaries.  Linear or quadratic interpolation
!  using local polynomial fit are supported (iorder=1 or 2).
!
!  3 March 1996 (KB)
!  Pressure is now interpolated with polynomials of pressure
!  rather than ln(p).
!
!  10 June 1997 (Fanyou Kong -- CMRP)
!  Include surface 2D data analysis from external data sets, a new
!  subroutine MKARPS2D is then added in file 'extlib23.f'.
!
!  16 July 1997 (Zonghui Huo)
!  Added 'getcoamps.f' to read the COAMPS data on sigma-z levels.
!
!  06/16/1998 (Dennis Moon, SSESCO)
!  Added support for RUC 2 on the AWIPS 211 grid. extdopt=7
!
!  06/16/1998 (Dennis Moon, SSESCO)
!  Added support for NCEP global re-analysis on T62 Gaussian lat/lon
!  grid. This does not work for ARPS grids which cross the 0 degree
!  meridian. extdopt=8
!
!  29 Oct 1998 (R. Carpenter, G. Bassett)
!  Added RUC-2 GEMPAK & ETA GEMPAK grid #104.  Four digit years and
!  extdname now used in constructing filenames of GEMPAK files.
!
!  15 Sep 1999 (K. Brewster)
!  Added the microphysical variables to the dtadump calls so they
!  will be used in those cases when available.  Added check for
!  positive definiteness to microphysical variables.
!
!  2000/03/23 (Donghai Wang)
!  Added NCEP AVN GRIB global data, grid #3.
!
!  2000/07/28 (Ming Xue)
!  Converted to F90 free format and implemented dynamic memory
!  allocation.
!
!  2000/08/14 (Gene Bassett)
!  Added multiple soil types, sfcdata variables and grid staggering
!  for use with arps history format of external model data.  Added options
!  allowing ext2arps to produce results similar to arpsintrp (intropt,
!  ext_lbc, ext_vbc).
!
!  2001/05/31 (Gene Bassett)
!  Added exttrnopt, an option to use input grid terrain for output grid.
!  Also corrected a bug introduced when updating the use of arps as an
!  external model (2000/08) which caused only the first external file
!  to have u&v rotated to arps map projection.
!
!  1 June 2002 Eric Kemp
!  Soil variable update.
!
!  6 June 2002 Eric Kemp
!  Added new soil variable interpolator for non-ARPS external data.
!
!  7 June 2002 Eric Kemp
!  Bug fix to 3d soil variable interpolation procedure.
!
!  15 June 2002 Eric Kemp
!  More changes to soil variable interpolation procedure.
!
!  24 June 2002 Eric Kemp
!  Bug fix with call to intrp_soil.  Added more temporary arrays for
!  use with this subroutine.
!
!  4 December 2002 Kevin W. Thomas
!  If one external file is bad, the entire run shouldn't abort.
!
!  2003-08-12 (Richard Carpenter)
!  Added extdopt=14,15 for GFS and Reanlysis data on grid 2 (2.5 deg lat/lon)
!
!  2003/10/23 Fanyou Kong
!  Modified to utilize 2-m temperature and humidity and 10-m wind
!  directly read from Eta GRIB (212) data
!
!  2004/01/15 Fanyou Kong
!  Modified to utilize 2-m temperature and humidity and 10-m wind
!  directly read from Eta AVN GRIB (#3) data
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input and output grid dimensions
!
!-----------------------------------------------------------------------

  INTEGER :: nx ! Number of grid points in the x-dir of output arps grid
  INTEGER :: ny ! Number of grid points in the y-dir of output arps grid
  INTEGER :: nz ! Number of grid points in the z-dir of output arps grid
  INTEGER :: nzsoil

  INTEGER :: nx_ext, ny_ext, nz_ext ! dimensions of external data grid
  INTEGER :: nzsoil_ext, nzsoil_extin

  INTEGER :: nstyps
!
!-----------------------------------------------------------------------
!
!  New arrays for Eta GRIB # 218 ( 12km Eta model Output)
!
!-----------------------------------------------------------------------

  REAL    :: latsw,lonsw, latne,lonne, latse,lonse, latnw,lonnw
                                            ! ARPS lat/lon at corners
  REAL    :: latpnt,lonpnt, latadj,lonadj

  INTEGER :: number_tile(54)                ! Total 54 tiles at most
  INTEGER :: npr, npc

  INTEGER :: ksw, kse, knw, kne, ktmp1, ktmp2
  INTEGER :: kswi, kswj, ksei, ksej, knwi, knwj, knei, knej

  INTEGER, ALLOCATABLE :: domain_tile(:,:)  ! Save tile order info.

  REAL, PARAMETER :: lat_min(54) = (/ &
  12.19,  14.07,  15.56,  16.64,  17.30,  17.33,  16.70,  15.64,  14.34, &
  19.65,  21.68,  23.29,  24.45,  25.16,  25.19,  24.51,  23.37,  21.97, &
  27.16,  29.32,  31.03,  32.26,  33.01,  33.04,  32.32,  31.12,  29.63, &
  34.58,  36.85,  38.63,  39.92,  40.70,  40.74,  39.99,  38.73,  37.17, &
  41.77,  44.12,  45.96,  47.29,  48.09,  48.13,  47.36,  46.06,  44.46, &
  48.61,  51.01,  52.89,  54.23,  55.05,  55.08,  54.30,  52.99,  51.35 /)

  REAL, PARAMETER :: lat_max(54) = (/  &
  21.55,  23.16,  24.33,  25.04,  25.30,  25.30,  25.07,  24.39,  23.25, &
  29.19,  30.90,  32.14,  32.89,  33.16,  33.16,  32.93,  32.20,  30.99, &
  36.71,  38.51,  39.80,  40.59,  40.87,  40.87,  40.62,  39.87,  38.60, &
  43.99,  45.84,  47.17,  47.99,  48.27,  48.27,  48.02,  47.24,  45.94, &
  50.89,  52.77,  54.12,  54.95,  55.23,  55.23,  54.98,  54.19,  52.87, &
  56.95,  58.84,  60.20,  61.02,  61.31,  61.31,  61.06,  60.27,  58.94 /)

  REAL, PARAMETER :: lon_min(54) = (/  &
    -135.76,-128.00,-120.01,-111.85,-103.56, -95.20, -87.33, -79.52, -71.82, &
    -138.38,-130.17,-121.69,-112.99,-104.14, -95.21, -86.84, -78.53, -70.35, &
    -141.35,-132.65,-123.61,-114.31,-104.82, -95.23, -86.28, -77.41, -68.70, &
    -144.74,-135.48,-125.82,-115.82,-105.60, -95.25, -85.63, -76.13, -66.80, &
    -148.63,-138.76,-128.39,-117.60,-106.51, -95.27, -84.89, -74.64, -64.62, &
    -152.88,-142.37,-131.23,-119.57,-107.53, -95.29, -84.01, -72.90, -62.08 /)

  REAL, PARAMETER :: lon_max(54) = (/  &
    -126.21,-118.66,-110.96,-103.16, -95.30, -86.96, -78.67, -70.49, -63.29, &
    -128.14,-120.15,-111.98,-103.68, -95.32, -86.41, -77.56, -68.85, -61.20, &
    -130.33,-121.84,-113.14,-104.28, -95.34, -85.78, -76.28, -66.97, -58.81, &
    -132.81,-123.77,-114.46,-104.97, -95.37, -85.05, -74.81, -64.80, -56.08, &
    -135.66,-125.99,-115.99,-105.76, -95.40, -84.19, -73.09, -62.28, -52.91, &
    -138.96,-128.58,-117.78,-106.69, -95.43, -83.23, -71.17, -59.48, -49.42 /)

!-----------------------------------------------------------------------
!
!  Space for mean sounding
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: lvlprof = 601
  REAL,    PARAMETER :: depthp  = 3.0E4

!
!-----------------------------------------------------------------------
!
!  lvlprof:
!
!  The ARPS interpolates unevenly-spaced data from the external
!  data set to evenly-spaced data for use in defining the base
!  state atmosphere.  In this process, an intermediate sounding is
!  generated at evenly-spaced altitudes, with the accuracy of the
!  associated interpolation controlled by the parameter lvlprof.
!  The larger lvlprof, the more accurate the interpolation
!  (we recommend using lvlprof=200 for a model run with about 50
!  points in the vertical direction).  Using the intermediate
!  sounding, the ARPS then generates a base state model sounding
!  for the particular vertical grid resolution chosen (i.e.,
!  the number of points in the vertical, nz, and the vertical grid
!  spacing, dz).
!
!  depthp:
!
!  The depth of atmosphere over which the interpolated profiles
!  will be defined.  depthp should be greater than or equal to
!  (nz-3)*dz, i.e., larger than the physical depth of the model
!  domain.  Otherwise, the code will extrapolate for gridpoints
!  outside the domain, leading to possible inconsistencies.  At all
!  costs, any such extrapolation should be avoided.
!
!-----------------------------------------------------------------------
!
  REAL ::  psnd(lvlprof),   zsnd(lvlprof),  tsnd(lvlprof),              &
          ptsnd(lvlprof), rhssnd(lvlprof), qvsnd(lvlprof),              &
         rhosnd(lvlprof),   usnd(lvlprof),  vsnd(lvlprof),              &
         dumsnd(lvlprof),  plsnd(lvlprof)
!
!-----------------------------------------------------------------------
!
!  ARPS grid variables
!
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at  a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    wbar     Base state vertical velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of soil model
!    hterain  Terrain height (m)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zp)/dz for soil model
!    j3soilinv  inverse of j3soil
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!
!    wetcanp  Canopy water amount
!
!-----------------------------------------------------------------------
!
  REAL, allocatable :: x(:)       ! The x-coord. of the physical and
                                  ! computational grid. Defined at u-point.
  REAL, allocatable :: y(:)       ! The y-coord. of the physical and
                                  ! computational grid. Defined at v-point.
  REAL, allocatable :: z(:)       ! The z-coord. of the computational grid.
                                  ! Defined at w-point on the staggered grid.
  REAL, allocatable :: xscl(:)    ! The x-coordinate of scalar points.
  REAL, allocatable :: yscl(:)    ! The y-coordinate of scalar points.
  REAL, allocatable :: zp(:,:,:)  ! The physical height coordinate defined at
                                  ! w-point of the staggered grid.
  REAL, allocatable :: zpsoil(:,:,:)  ! The physical height coordinate defined at
                                      ! w-point in the soil.
  REAL, allocatable :: zs(:,:,:)      ! The physical height of scalar points.
  REAL, allocatable :: j1(:,:,:)      ! Coordinate transformation Jacobian -d(zp)/dx.
  REAL, allocatable :: j2(:,:,:)      ! Coordinate transformation Jacobian -d(zp)/dy.
  REAL, allocatable :: j3(:,:,:)      ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL, allocatable :: j3soil(:,:,:)     ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL, allocatable :: j3soilinv(:,:,:)  ! Coordinate transformation Jacobian  d(zp)/dz.

  REAL, allocatable :: aj3z(:,:,:)  ! Coordinate transformation Jacobian defined
                                    ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL, allocatable :: hterain(:,:) ! The height of the terrain. (m)
  REAL, allocatable :: htrn_tmp(:,:)! The height of the terrain. (m)
  REAL, allocatable :: mapfct(:,:,:)! Map factor at scalar, u, and v points

  REAL, allocatable :: u(:,:,:)     ! Total u-velocity (m/s)
  REAL, allocatable :: v(:,:,:)     ! Total v-velocity (m/s)
  REAL, allocatable :: w(:,:,:)     ! Total w-velocity (m/s)
  REAL, allocatable :: pprt(:,:,:)  ! Perturbation pressure from that
                                    ! of base state atmosphere (Pascal)
  REAL, allocatable :: ptprt(:,:,:) ! Perturbation potential temperature
                                    ! from that of base state atmosphere (K)

  REAL, allocatable :: qv(:,:,:)    ! Water vapor specific humidity (kg/kg)

  REAL, allocatable :: qscalar(:,:,:,:)  ! Microphysics scalar array

  REAL, allocatable :: pbar(:,:,:)  ! Base state pressure (Pascal)
  REAL, allocatable :: ptbar(:,:,:) ! Base state potential temperature (K)
  REAL, allocatable :: qvbar(:,:,:) ! Base state water vapor specific humidity
                                    ! (kg/kg)
  REAL, allocatable :: ubar(:,:,:)  ! Base state u-velocity (m/s)
  REAL, allocatable :: vbar(:,:,:)  ! Base state v-velocity (m/s)
  REAL, allocatable :: wbar(:,:,:)  ! Base state w-velocity (m/s)

  REAL, allocatable :: rhobar(:,:,:)! Base state density (kg/m3).
  REAL, allocatable :: rhostr(:,:,:)! Base state density rhobar times j3.
  REAL, allocatable :: wcont(:,:,:) ! Contravariant vertical velocity (m/s)

  REAL, allocatable :: csndsq(:,:,:)! Speed of sound squared (m**2/s**2)

  REAL, allocatable :: tsoil(:,:,:,:)  ! Soil temperature (K)
  REAL, allocatable :: qsoil(:,:,:,:)  ! Soil moisture (m**3/m**3)
  REAL, allocatable :: wetcanp(:,:,:)  ! Water amount on canopy
  REAL, allocatable :: snowdpth(:,:)   ! Snow depth (m)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:)! Soil type
  REAL,    ALLOCATABLE :: stypfrct(:,:,:)! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp (:,:)   ! Vegetation type
  REAL,    ALLOCATABLE :: lai    (:,:)   ! Leaf Area Index
  REAL,    ALLOCATABLE :: roufns (:,:)   ! Surface roughness
  REAL,    ALLOCATABLE :: veg    (:,:)   ! Vegetation fraction

  REAL, allocatable :: tem1(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem2(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem3(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem4(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem5(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem6(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem7(:,:,:)   ! Temporary work array.
  REAL, allocatable :: tem8(:,:,:)   ! Temporary work array.

  REAL, ALLOCATABLE    :: tem9 (:,:)     ! Temporary work array.

  INTEGER, ALLOCATABLE :: ksoil3d(:,:,:) ! EMK 6 June 2002

  REAL, ALLOCATABLE :: xs2d(:,:)
  REAL, ALLOCATABLE :: ys2d(:,:)
  REAL, ALLOCATABLE :: xu2d(:,:)
  REAL, ALLOCATABLE :: yu2d(:,:)
  REAL, ALLOCATABLE :: xv2d(:,:)
  REAL, ALLOCATABLE :: yv2d(:,:)
  REAL, ALLOCATABLE :: xw(:,:)
  REAL, ALLOCATABLE :: yw(:,:)

  INTEGER, ALLOCATABLE :: iscl(:,:) ! i index of scalar point
                                    ! in external array.
  INTEGER, ALLOCATABLE :: jscl(:,:) ! j index of scalar point
                                    ! in external array.
  INTEGER, ALLOCATABLE :: iu(:,:)   ! i index of u-velocity point
                                    ! in external array.
  INTEGER, ALLOCATABLE :: ju(:,:)   ! j index of u-velocity point
                                    ! in external array.
  INTEGER, ALLOCATABLE :: iv(:,:)   ! i index of v-velocity point
                                    ! in external array.
  INTEGER, ALLOCATABLE :: jv(:,:)   ! j index of v-velocity point
                                    ! in external array.

  INTEGER :: iproj_arps              ! ARPS map projection (tem storage)
  REAL    :: scale_arps              ! ARPS map scale (tem storage)
  REAL    :: trlon_arps              ! ARPS true-longitude (tem storage)
  REAL    :: latnot_arps(2)          ! ARPS true-latitude(s) (tem storage)
  REAL    :: xorig_arps,yorig_arps   ! ARPS map projection origin (tem storage)

  REAL, ALLOCATABLE :: temsoil1(:,:,:)
  REAL, ALLOCATABLE :: temsoil2(:,:,:)
  REAL, ALLOCATABLE :: temsoil3(:,:,:)
  REAL, ALLOCATABLE :: temsoil4(:,:,:)
!
!-----------------------------------------------------------------------
!
!  NAMELIST parameter (In arps.input)
!
!-----------------------------------------------------------------------
!
  INCLUDE 'ext2arps.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  External forecast variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext            ! external data map projection
  REAL    :: scale_ext            ! external data map scale factor
  REAL    :: trlon_ext            ! external data true longitude
  REAL    :: latnot_ext(2)        ! external data true latitude(s)
  REAL    :: x0_ext,y0_ext        ! external data origin

  REAL, allocatable :: x_ext(:)       ! external data x-coordinate
  REAL, allocatable :: y_ext(:)       ! external data y-coordinate
  REAL, allocatable :: xu_ext(:)       ! external data u x-coordinate
  REAL, allocatable :: yu_ext(:)       ! external data u y-coordinate
  REAL, allocatable :: xv_ext(:)       ! external data v x-coordinate
  REAL, allocatable :: yv_ext(:)       ! external data v y-coordinate
  REAL, allocatable :: lat_ext(:,:)    ! external data latidude
  REAL, allocatable :: lon_ext(:,:)    ! external data longitude
  REAL, allocatable :: latu_ext(:,:)   ! external data latidude (x-stag)
  REAL, allocatable :: lonu_ext(:,:)   ! external data longitude (x-stag)
  REAL, allocatable :: latv_ext(:,:)   ! external data latidude (y-stag)
  REAL, allocatable :: lonv_ext(:,:)   ! external data longitude (y-stag)
  REAL, allocatable :: zs_ext(:,:,:)   ! external data physical height (m)

  REAL, allocatable :: p_ext(:,:,:)    ! Pressure (Pascals)
  REAL, allocatable :: hgt_ext(:,:,:)  ! Height (m)
  REAL, allocatable :: zp2_ext(:,:,:)  ! external w physical height (m)
  REAL, allocatable :: zp_ext(:,:,:)   ! Height (m) (on arps grid)
  REAL, allocatable :: zpsoil_ext(:,:,:)   ! Height (m) (on arps grid)

  REAL, allocatable :: t_ext(:,:,:)    ! Temperature (K)
  REAL, allocatable :: u_ext(:,:,:)    ! Eastward wind component
  REAL, allocatable :: v_ext(:,:,:)    ! Northward wind component
  REAL, allocatable :: w_ext(:,:,:)    ! Vertical wind component
  REAL, allocatable :: qv_ext(:,:,:)   ! Specific humidity (kg/kg)

  REAL, allocatable :: qscalar_ext(:,:,:,:) ! external scalar array

  REAL, allocatable :: tsoil_ext  (:,:,:,:)   ! Soil temperature  (K)
  REAL, allocatable :: qsoil_ext  (:,:,:,:)   ! Soil moisture (m3/m3)
  REAL, allocatable :: wetcanp_ext(:,:,:)   ! Canopy water amount
  REAL, allocatable :: snowdpth_ext(:,:)  ! Snow depth (m)

  INTEGER, allocatable :: soiltyp_ext (:,:,:)! Soil type
  REAL,    allocatable :: stypfrct_ext(:,:,:)! Soil type fraction
  INTEGER, allocatable :: vegtyp_ext (:,:)   ! Vegetation type
  REAL, allocatable ::    lai_ext    (:,:)   ! Leaf Area Index
  REAL, allocatable ::    roufns_ext (:,:)   ! Surface roughness
  REAL, allocatable ::    veg_ext    (:,:)   ! Vegetation fraction

  REAL, allocatable :: pbar_ext(:,:,:)  ! Base state pressure (Pascal)
  REAL, allocatable :: ptbar_ext(:,:,:) ! Base state potential temperature (K)
  REAL, allocatable :: qvbar_ext(:,:,:) ! Base state water vapor specific
                                        ! humidity (kg/kg)
  REAL, allocatable :: ubar_ext(:,:,:)  ! Base state u-velocity (m/s)
  REAL, allocatable :: vbar_ext(:,:,:)  ! Base state v-velocity (m/s)
  REAL, allocatable :: wbar_ext(:,:,:)  ! Base state w-velocity (m/s)


  REAL, allocatable :: trn_ext    (:,:)   ! External terrain (m)
  REAL, allocatable :: psfc_ext   (:,:)   ! Surface pressure (Pa)

  REAL, allocatable :: t_2m_ext (:,:)     ! 2-m temperature (K)
  REAL, allocatable :: qv_2m_ext(:,:)     ! 2-m specific humidity (kg/kg)
  REAL, allocatable :: u_10m_ext(:,:)     ! 10-m u (m/s)
  REAL, allocatable :: v_10m_ext(:,:)     ! 10-m v (m/s)
  REAL, allocatable :: rain_ext(:,:)     !
  REAL, allocatable :: t_2m (:,:)   ! in arps domain
  REAL, allocatable :: qv_2m(:,:)   ! in arps domain
  REAL, allocatable :: u_10m(:,:)   ! in arps domain
  REAL, allocatable :: v_10m(:,:)   ! in arps domain
  REAL, allocatable :: rain(:,:)   !

  REAL, allocatable :: dxfld(:)
  REAL, allocatable :: dyfld(:)
  REAL, allocatable :: rdxfld(:)
  REAL, allocatable :: rdyfld(:)
  REAL, allocatable :: dxfldu(:)
  REAL, allocatable :: dyfldu(:)
  REAL, allocatable :: rdxfldu(:)
  REAL, allocatable :: rdyfldu(:)
  REAL, allocatable :: dxfldv(:)
  REAL, allocatable :: dyfldv(:)
  REAL, allocatable :: rdxfldv(:)
  REAL, allocatable :: rdyfldv(:)

  REAL, ALLOCATABLE :: rdzsoilfld(:,:,:) ! EMK 17 June 2002

  REAL, ALLOCATABLE :: tem1_ext(:,:,:), tem2_ext(:,:,:), tem3_ext(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem4_ext(:,:,:), tem5_ext(:,:,:)                     ! Temporary work array

  REAL, allocatable :: xa_ext(:,:), ya_ext(:,:), avgzs_ext(:,:,:)

  REAL, ALLOCATABLE :: temsoil1_ext(:,:,:), temsoil2_ext(:,:,:)
  REAL, ALLOCATABLE :: temsoil3_ext(:,:,:), temsoil4_ext(:,:,:)

  INTEGER :: ilite,iltgci
  INTEGER :: ibeg_offset,iend_offset,jbeg_offset,jend_offset
!
!-----------------------------------------------------------------------
!
!  include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'adjust.inc'
!
!-----------------------------------------------------------------------
!
!  Non-dimensional smoothing coefficient
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: smfct1 = 0.5, smfct2 = -0.5
  REAL, PARAMETER :: rhmax  = 1.0
!
!-----------------------------------------------------------------------
!
!  Latitude and longitude for some diagnostic printing,
!  e.g. to compare to an observed sounding
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: latdiag = 34.5606, londiag = -103.0820
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: basdmpfn
  CHARACTER (LEN=19)  :: extdinit
  CHARACTER (LEN=9)   :: extdfcst
  CHARACTER (LEN=9)   :: julfinit
  CHARACTER (LEN=9)   :: julfname

  INTEGER :: i,j,k,nq,ksmth,istatus
  INTEGER :: iyr,imo,iday,ihr,imin,isec,jldy
  INTEGER :: ifhr,ifmin,ifsec,mfhr
  REAL    :: time_ext
  INTEGER :: myr,initsec,iabssec,jabssec,kftime
  INTEGER :: ifile,iprtopt,iniotfu,lbasdmpf,onvf,grdbas
  INTEGER :: first_time
  INTEGER :: iextmn,iextmx,jextmn,jextmx
  INTEGER :: idiag,jdiag
  INTEGER :: idiag_arps,jdiag_arps
  INTEGER :: hdmpgrdfmt

  REAL :: latnot(2)
  REAL :: amin,amax
  REAL :: qvmin,qvmax,qvval
  REAL :: csconst,pconst
  REAL :: deltaz,tv_ext
  REAL :: pres,temp,qvsat,rh,tvbar,qvprt,qtot
  REAL :: xdiag,ydiag,dd,dmin,latd,lond
  REAL :: ppasc,pmb,tc,tdc,theta,smix,e,bige,alge,dir,spd
  REAL :: dx_ext, dy_ext

  CHARACTER (LEN=256) :: soiloutfl,temchar,ternfn
  CHARACTER (LEN=80)  :: timsnd
  INTEGER :: lfn, tmstrln, lternfn
  INTEGER :: wetcout,snowdout
  INTEGER :: tsoilout,qsoilout,zpsoilout
  INTEGER :: isnow,jsnow,ii,jj
  REAL    :: xumin,xumax,yvmin,yvmax
  REAL    :: xsmin,xsmax,ysmin,ysmax

  INTEGER :: lenbin

  INTEGER             :: strlen, ireturn

  CHARACTER (LEN=3)   :: fmtstr(10) = (/'bin','asc','hdf','pak','xxx',  &
                                        'bn2','net','nc ','xxx','grb' /)

  INTEGER :: is,nstyp_ext

  DOUBLE PRECISION :: ntmergeinv, mfac
  INTEGER :: idist
  LOGICAL :: lon_0_360 = .FALSE.  ! global lat/lon grids in which lon runs 0:360 not -180:180

!  INTEGER :: grbfmt = 0, nofixdim = 0
  INTEGER :: nofixdim = 0
  LOGICAL :: median_warn = .FALSE.
  LOGICAL :: stagger_ext = .FALSE.

  CHARACTER(LEN=15), PARAMETER :: subnames(0:51) = (/     'getarps        ',  &
    'getnmcruc87    ','getnmceta212   ','getlaps        ','getgemruc      ',  &
    'getgemeta      ','getcoamps      ','getnmcruc211   ','getreanalt62   ',  &
    'getgemruc2     ','getgemeta2     ','getnmcrucn236  ','getnmcrucp236  ',  &
    'getncepavn3    ','getncepavn2    ','getncepavn2    ','getnmceta218   ',  &
    'getnarr221     ','getnmcrucn236  ','getnmcrucp236  ','getncepgfs     ',  &
    'getrucn130     ','getecmf128     ','getwrfdata     ','getnmceta132   ',  &
    'xxx            ','xxx            ','xxx            ','xxx            ',  &
    'xxx            ','xxx            ','xxx            ','xxx            ',  &
    'xxx            ','xxx            ','xxx            ','xxx            ',  &
    'xxx            ','xxx            ','xxx            ','xxx            ',  &
    'xxx            ','xxx            ','xxx            ','xxx            ',  &
    'xxx            ','xxx            ','xxx            ','xxx            ',  &
    'xxx            ','get_avn_bin    ','getncepavn3    ' /)

  CHARACTER(LEN=40 ) :: qnames_wrf(20)

  CHARACTER(LEN=256) :: name_nml_dir
  INTEGER :: unum         ! unit number for reading in namelist

!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = 1
  nestgrd = 0

  nstyps = 4

  nstyp_ext  = 1
  nzsoil_ext = 2

  DO k=1,lvlprof
    dumsnd(k) = 0.0
  END DO
!
!-----------------------------------------------------------------------
!
!  Call initpara to read in ARPS parameters of namelists
!
!-----------------------------------------------------------------------
!
  name_nml_dir = ' '
  CALL initpara(nx,ny,nz,nzsoil,nstyps, name_nml_dir)

  CALL initadas (name_nml_dir)  ! The adjust namelist block is used by ext2arps.  Other
                ! adas parameters are not used by this program,
                ! but the name list blocks need to be read in sequence on
                ! Cray machines such as J90.

  IF (myproc == 0) THEN
    myr=MOD(year,100)
    WRITE(julfinit,'(i2.2,i3.3,i2.2,i2.2)') myr,jday,hour,minute
    CALL ctim2abss(year,month,day,hour,minute,second,initsec)
  END IF

  CALL mpupdatei(initsec,1)

!
!-----------------------------------------------------------------------
!
!  Read in additional namelists for external file specifications.
!
!-----------------------------------------------------------------------
!
  extdopt  = 0
  extdfmt  = 10
  dir_extd = './'
  extdname = 'may20'
  nextdfil = 1
  extdtime(1) = '1977-05-20.21:00:00+000:00:00'
  iorder   = 2
  intropt  = 1
  nsmooth  = 1
  ext_lbc  = 1
  ext_vbc  = 1
  exttrnopt = 0
  extntmrg = 1
  extsfcopt = 0
  grdbasopt = 1
  i2dfmt = 0
  outheader = 'data2d'
  iboxs = 0
  iboxe = -1
  jboxs = 0
  jboxe = -1

  IF (myproc == 0) THEN
    CALL getunit( unum )
    OPEN(unum,FILE=TRIM(name_nml_dir),STATUS='OLD',FORM='FORMATTED',IOSTAT=istatus)

    IF (istatus == 0) THEN
      WRITE(6,'(1x,3a,/,1x,a,/)') 'Reading namelist EXTDFILE from file - ',   &
              TRIM(name_nml_dir),' ... ','========================================'
    ELSE
      CALL retunit(unum)
      unum = 5
      WRITE(6,'(2(1x,a,/))') 'Waiting namelist EXTDFILE from standard input ... ', &
                             '========================================'
    END IF

    READ (unum,extdfile)
  END IF

  CALL mpupdatei(extdopt,1)
  CALL mpupdatei(extdfmt,1)
  CALL mpupdatec(dir_extd,256)
  CALL mpupdatec(extdname,256)
  CALL mpupdatei(iorder,1)
  CALL mpupdatei(intropt,1)
  CALL mpupdatei(nsmooth,1)
  CALL mpupdatei(hydradj,1)
  CALL mpupdatei(wndadj,1)
  CALL mpupdatei(obropt,1)
  CALL mpupdater(obrzero,1)
  CALL mpupdatei(ext_lbc,1)
  CALL mpupdatei(ext_vbc,1)
  CALL mpupdatei(exttrnopt,1)
  CALL mpupdatei(extntmrg,1)
  CALL mpupdatei(extsfcopt,1)
  CALL mpupdatei(grdbasopt,1)
  CALL mpupdatei(nextdfil,1)
  CALL mpupdatec(extdtime,nextdfil*29)
  CALL mpupdatei(i2dfmt,1)
  CALL mpupdatec(outheader,256)

  CALL mpupdatei(iboxs,1)
  CALL mpupdatei(iboxe,1)
  CALL mpupdatei(jboxs,1)
  CALL mpupdatei(jboxe,1)

  strlen = LEN_TRIM(dir_extd)
  IF ( strlen == 0 .OR. dir_extd(1:strlen) == ' ' ) THEN
    dir_extd = './'
    strlen = 2
  END IF

  IF( dir_extd(strlen:strlen) /= '/' ) THEN
    strlen = strlen + 1
    dir_extd(strlen:strlen) = '/'
  END IF

  IF (myproc == 0) THEN

    WRITE (6, '(/1x,a)')    '&extdfile'

    WRITE (6, '(3x,a,a,a)')                                             &
                       'dir_extd   = ''', TRIM(dir_extd), ''','

    WRITE (6, '(3x,a,a,a)')                                             &
                       'extdname   = ''', TRIM(extdname), ''','

    WRITE (6, '(3x,a,i4,a)')    'extdopt   = ', extdopt,  ','
    WRITE (6, '(3x,a,i4,a)')    'extdfmt   = ', extdfmt,  ','
    WRITE (6, '(3x,a,i4,a)')    'nextdfil  = ', nextdfil, ','

    DO i=1,nextdfil
      WRITE (6, '(3x,a,i2.2,a,a,a)')                                    &
                       'extdtime(',i,') = ''', TRIM(extdtime(i)), ''','
    END DO

    WRITE (6, '(3x,2(a,i4),a)') 'iboxs     = ', iboxs,  ',   iboxe     = ', iboxe, ','
    WRITE (6, '(3x,2(a,i4),a)') 'jboxs     = ', jboxs,  ',   jboxe     = ', jboxe, ','

    WRITE (6, '(3x,a,i4,a)')    'iorder    = ', iorder,   ','
    WRITE (6, '(3x,a,i4,a)')    'intropt   = ', intropt,  ','
    WRITE (6, '(3x,a,i4,a)')    'nsmooth   = ', nsmooth,  ','
    WRITE (6, '(3x,a,i4,a)')    'hydradj   = ', hydradj,  ','
    WRITE (6, '(3x,a,i4,a)')    'wndadj    = ', wndadj,   ','
    WRITE (6, '(3x,a,i4,a)')    'obropt    = ', obropt,   ','
    WRITE (6, '(3x,a,f16.4,a)') 'obrzero   = ', obrzero,  ','
    WRITE (6, '(3x,a,i4,a)')    'ext_lbc   = ', ext_lbc,  ','
    WRITE (6, '(3x,a,i4,a)')    'ext_vbc   = ', ext_vbc,  ','
    WRITE (6, '(3x,a,i4,a)')    'exttrnopt = ', exttrnopt,','
    WRITE (6, '(3x,a,i4,a)')    'extntmrg  = ', extntmrg, ','
    WRITE (6, '(3x,a,i4,a)')    'extsfcopt = ', extsfcopt, ','
    WRITE (6, '(3x,a,i4,a)')    'grdbasopt = ', grdbasopt, ','
    WRITE (6, '(3x,a,i4,a)')    'i2dfmt    = ', i2dfmt,   ','
    WRITE (6, '(3x,a,a,a)')    'outheader = ''', TRIM(outheader), ''','
    WRITE (6, '(1x,a)')         '&END'

  END IF

  nofixdim = extdopt / 100
  extdopt  = MOD(extdopt,100)

  IF (nofixdim == 1 .OR. extdopt == 11 .OR. extdopt == 12 .OR.          &
      extdopt == 6  .OR. extdopt == 22 .OR. extdopt == 23 ) THEN  ! Initialize extdinit & extdfcst
    READ(extdtime(1),'(a19,1x,a9)') extdinit,extdfcst
    IF(extdfcst == '         ') extdfcst='000:00:00'
  END IF

  IF (myproc == 0 .AND. unum /= 5) THEN
    CLOSE(unum)
    CALL retunit(unum)
  END IF

!-----------------------------------------------------------------------
!
!  ALLOCATE and initalize arrays based on dimension parameters
!  read in from the input file
!
!-----------------------------------------------------------------------

  ALLOCATE(x(nx),stat=istatus)
  ALLOCATE(y(ny),stat=istatus)
  ALLOCATE(z(nz),stat=istatus)

  ALLOCATE(xscl(nx),stat=istatus)
  ALLOCATE(yscl(ny),stat=istatus)
  ALLOCATE(zp(nx,ny,nz),stat=istatus)
  ALLOCATE(zpsoil(nx,ny,nzsoil),stat=istatus)

  ALLOCATE(zs(nx,ny,nz),stat=istatus)
  ALLOCATE(j1(nx,ny,nz),stat=istatus)
  ALLOCATE(j2(nx,ny,nz),stat=istatus)
  ALLOCATE(j3(nx,ny,nz),stat=istatus)
  ALLOCATE(j3soil(nx,ny,nzsoil),stat=istatus)
  ALLOCATE(j3soilinv(nx,ny,nzsoil),stat=istatus)

  ALLOCATE(aj3z(nx,ny,nz),stat=istatus)

  ALLOCATE(hterain(nx,ny),stat=istatus)
  ALLOCATE(htrn_tmp(nx,ny),stat=istatus)
  ALLOCATE(mapfct(nx,ny,8),stat=istatus)

  ALLOCATE(u(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:U")
  ALLOCATE(v(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:V")
  ALLOCATE(w(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:W")
  ALLOCATE(pprt(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:PPRT")
  ALLOCATE(ptprt(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:PTPRT")
  ALLOCATE(qv(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:QV")

!  ALLOCATE(qscalar(nx,ny,nz,nscalar),stat=istatus)
!  CALL check_alloc_status(istatus, "ext2arps:QSCALAR")

  ALLOCATE(pbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:PBAR")
  ALLOCATE(ptbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:PTBAR")
  ALLOCATE(qvbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:QVBAR")
  ALLOCATE(ubar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:UBAR")
  ALLOCATE(vbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:VBAR")
  ALLOCATE(wbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:WBAR")
  ALLOCATE(rhobar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:RHOBAR")
  ALLOCATE(rhostr(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:RHOSTR")
  ALLOCATE(wcont(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:WCONT")
  ALLOCATE(csndsq(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:CSNDSQ")

  ALLOCATE(wetcanp(nx,ny,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:WETCANP")
  ALLOCATE(snowdpth(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:SNOWDPTH")

  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:TSOIL")
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:QSOIL")

  ALLOCATE(soiltyp (nx,ny,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:soiltyp")
  ALLOCATE(stypfrct(nx,ny,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:sypfrct")
  ALLOCATE(vegtyp (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:vegtyp")
  ALLOCATE(lai    (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:lai")
  ALLOCATE(roufns (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:roufns")
  ALLOCATE(veg    (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:veg")

  ALLOCATE(t_2m (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:t_2m")
  ALLOCATE(qv_2m(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:qv_2m")
  ALLOCATE(u_10m(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:u_10m")
  ALLOCATE(v_10m(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:v_10m")
  ALLOCATE(rain(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:rain")

  ALLOCATE(tem1(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem1")
  ALLOCATE(tem2(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem2")
  ALLOCATE(tem3(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem3")
  ALLOCATE(tem4(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem4")
  ALLOCATE(tem5(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem5")
  ALLOCATE(tem6(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem6")
  ALLOCATE(tem7(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem7")
  ALLOCATE(tem8(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem8")

  ALLOCATE(tem9(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem9")

  ALLOCATE(ksoil3d(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "ext2arps:ksoil3d")
  ksoil3d(:,:,:) = 0

  ALLOCATE(xs2d(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:xs2d")
  ALLOCATE(ys2d(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:ys2d")
  ALLOCATE(xu2d(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:xu2d")
  ALLOCATE(yu2d(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:yu2d")
  ALLOCATE(xv2d(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:xv2d")
  ALLOCATE(yv2d(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:yv2d")
  ALLOCATE(xw(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:xw")
  ALLOCATE(yw(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:yw")

  ALLOCATE(iscl(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:iscl")
  ALLOCATE(jscl(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:jscl")
  ALLOCATE(iu(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:iu")
  ALLOCATE(ju(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:ju")
  ALLOCATE(iv(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:iv")
  ALLOCATE(jv(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:jv")


  x=0.0
  y=0.0
  z=0.0

  xscl=0.0
  yscl=0.0
  zp=0.0

  zs=0.0
  j1=0.0
  j2=0.0
  j3=0.0
  aj3z=0.0
  hterain=0.0
  htrn_tmp=0.0
  mapfct=0.0

  u=0.0
  v=0.0
  w=0.0
  pprt=0.0
  ptprt=0.0
  qv=0.0
!  qscalar=0.0

  pbar=0.0
  ptbar=0.0
  qvbar=0.0
  ubar=0.0
  vbar=0.0
  wbar=0.0
  rhobar=0.0
  rhostr=0.0
  wcont=0.0
  csndsq=0.0

  tsoil=0.0
  qsoil=0.0
  wetcanp=0.0
  snowdpth=0.0

  soiltyp =0
  stypfrct=0.0
  stypfrct(:,:,1) =1.0
  vegtyp =0
  lai    =0.0
  roufns =0.0
  veg    =0.0

  t_2m =0.0
  qv_2m=0.0
  u_10m=0.0
  v_10m=0.0
  rain=0.0

  tem1=0.0
  tem2=0.0
  tem3=0.0
  tem4=0.0
  tem5=0.0
  tem6=0.0
  tem7=0.0
  tem8=0.0

  xs2d=0.0
  ys2d=0.0
  xu2d=0.0
  yu2d=0.0
  xv2d=0.0
  yv2d=0.0

  iscl=0
  jscl=0
  iu=0
  ju=0
  iv=0
  jv=0
!
!-----------------------------------------------------------------------
!
!  Set up ARPS map projection
!
!-----------------------------------------------------------------------
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
!
!-----------------------------------------------------------------------
!
!  Set up ARPS grid
!
!-----------------------------------------------------------------------
!
  IF (exttrnopt == 1) THEN ! don't really do grid here, just set map proj
    ternopt = -1
    hterain = 0
  END IF

  CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                         &
              hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,                &
              tem2(1,1,:),tem3(1,1,:),tem1)

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        aj3z(i,j,k)=0.5*(j3(i,j,k)+j3(i,j,k-1))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Get lat/lon at ARPS domain corns. They are used for the following data sources
!
!  extdopt = 16 (ETA grid #218)  - To determine the tiltes to be read
!  extdopt = 8, 13, 14, 15 (Global data source)
!            - To set median_warn if Greenwich Meridian goes through the
!              ARPS domain so that the data arrays should be rearranged.
!
!-----------------------------------------------------------------------

  CALL xytoll(1,1, x(1), y(1),latsw,lonsw)
  CALL xytoll(1,1,x(nx), y(1),latse,lonse)
  CALL xytoll(1,1,x(nx),y(ny),latne,lonne)
  CALL xytoll(1,1, x(1),y(ny),latnw,lonnw)

  CALL mpbcastr(latsw,0)
  CALL mpbcastr(lonsw,0)
  CALL mpbcastr(latne,nproc_x*nproc_y-1)
  CALL mpbcastr(lonne,nproc_x*nproc_y-1)
  CALL mpbcastr(latse,nproc_x-1)
  CALL mpbcastr(lonse,nproc_x-1)
  CALL mpbcastr(latnw,nproc_x*(nproc_y-1))
  CALL mpbcastr(lonnw,nproc_x*(nproc_y-1))

  IF (myproc == 0) WRITE(6,'(/a,4(a,F7.2),a/,13x,4(a,F7.2),a/)')        &
              ' ARPS domain ',                                          &
              'NW: (',latnw,',',lonnw, ') NE: (',latne,',',lonne,')',   &
              'SW: (',latsw,',',lonsw, ') SE: (',latse,',',lonse,').'

  IF ( (lonsw < 0 .AND. lonse > 0) .OR. (lonnw < 0 .AND. lonne > 0) )   &
      median_warn = .TRUE.   ! Greenwich Meridian goes through
!
!-----------------------------------------------------------------------
!
!  Set up the height levels on which to define the base-state.
!
!-----------------------------------------------------------------------
!
  deltaz = depthp/(lvlprof-1)

  IF (myproc == 0) WRITE(6,'(/a,/a,f7.2,a/)')                           &
      ' The base state is formed from a mean sounding created',         &
      ' with a ',deltaz,' meter interval.'

  DO k=1,lvlprof
    zsnd(k) = (k-1)*deltaz
  END DO

!-----------------------------------------------------------------------
!
!  Set external grid size
!
!-----------------------------------------------------------------------
! DTD: Overwrite number of microphysical scalars, flags, names, and descriptions
! Set by initpara above. initpara initializes these parameters based on the parameter
! mphyopt in the arps namelist input.  However, in ext2arps these parameters should
! be set according to what is available in the external data.

  nscalar  = 0
  nscalarq = 0
  P_QC = -1; P_QR = -1; P_QI = -1; P_QS = -1; P_QH = -1; P_QG = -1
  P_NC = -1; P_NR = -1; P_NI = -1; P_NS = -1; P_NG = -1; P_NH = -1
             P_ZR = -1; P_ZI = -1; P_ZS = -1; P_ZG = -1; P_ZH = -1;
  qnames(:)= ' '; qdescp(:)= ' '

  SELECT CASE ( extdopt )

    CASE ( 0 )  ! ARPS grid

    IF (myproc == 0) THEN

      time_ext = FLOAT( (ifhr*3600)+(ifmin*60)+ifsec )
      CALL cvttsnd( time_ext, timsnd, tmstrln )

      CALL get_dims_from_data(extdfmt,                                    &
         trim(dir_extd)//trim(extdname)//'.'//TRIM(fmtstr(extdfmt))//timsnd(1:tmstrln), &
         nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext, ireturn)
      nstyp_ext = min(nstyp_ext,nstyps)

      IF( ireturn /= 0 ) CALL arpsstop(                                 &
        'Problem occurred when trying to get dimensions from data.', 1)
    END IF

    CALL mpupdatei(nx_ext,1)
    CALL mpupdatei(ny_ext,1)
    CALL mpupdatei(nz_ext,1)
    CALL mpupdatei(nzsoil_ext,1)
    CALL mpupdatei(nstyp_ext,1)

    CALL mpupdatei(nscalar, 1)
    CALL mpupdatei(nscalarq,1)
    CALL mpupdatei(P_QC,1);
    CALL mpupdatei(P_QR,1)
    CALL mpupdatei(P_QI,1)
    CALL mpupdatei(P_QS,1)
    CALL mpupdatei(P_QG,1)
    CALL mpupdatei(P_QH,1)
    CALL mpupdatei(P_NC,1)
    CALL mpupdatei(P_NR,1)
    CALL mpupdatei(P_NI,1)
    CALL mpupdatei(P_NS,1)
    CALL mpupdatei(P_NG,1)
    CALL mpupdatei(P_NH,1)
    CALL mpupdatei(P_ZR,1)
    CALL mpupdatei(P_ZI,1)
    CALL mpupdatei(P_ZS,1)
    CALL mpupdatei(P_ZG,1)
    CALL mpupdatei(P_ZH,1)

    CALL mpupdatec(qnames,nscalar*40)
    CALL mpupdatec(qdescp,nscalar*40)

    CASE ( 1 )  ! RUC (Hybrid-B, GRIB #87) from NMC is 81x62x25
      nx_ext=81
      ny_ext=62
      nz_ext=25
      nzsoil_ext = 2       ! FIXME:  For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 2 )  ! ETA (GRIB #212, 40km) from NMC is 185x129x39
      nx_ext = 185
      ny_ext = 129
      nz_ext = 39
      nzsoil_ext = 5                ! EMK 15 June 2002
      nstyp_ext = 1                 ! EMK 20 June 2002

      nscalar  = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 3 ) ! OLAPS for 95 is 91x73x41
      nx_ext=91
      ny_ext=73
      nz_ext=41
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 4 ) ! RUC in GEMPAK format from Rossby (e.g., those of 1995)
      nx_ext=81
      ny_ext=62
      nz_ext=41
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 5 ) ! GEMPAK ETA data

      CALL arpsstop(                                                      &
          'GEMPAK ETA: Need to find out the grid size and set in ext2arps.f90.',1)

!      nx_ext = ??
!      ny_ext = ??
!      nz_ext = ??

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 6 ) ! Special COAMPS files

      !PRINT*,'Please specify the size of COAMPS input data grid '
      !PRINT*,'by editing file ext2arps.f90.'
      !call arpsstop("Processing COAMPS files",1)
      !nx_ext = ??
      !ny_ext = ??
      !nz_ext = ??
      IF (myproc == 0) THEN
        CALL getcoampsdims(dir_extd,extdinit,extdfcst,                  &
                           nx_ext, ny_ext, nz_ext,istatus)
      END IF
      CALL mpupdatei(istatus,1)
      IF (istatus /= 0) THEN
        CALL arpsstop('Error retrieve COAMPS dimensions.',1)
      ELSE
        CALL mpupdatei(nx_ext,1)
        CALL mpupdatei(ny_ext,1)
        CALL mpupdatei(nz_ext,1)
        IF (myproc == 0) THEN
          WRITE(6,'(1x,a,3(I4,a))') 'COAMPS dimension are: nx_ext = ',  &
                    nx_ext,', ny_ext = ',ny_ext,', nz_ext = ',nz_ext,'.'
        END IF
      END IF
      nzsoil_ext = 2

      nscalar = 5
      nscalarq = 5
      P_QC = 1
      P_QR = 2
      P_QI = 3
      P_QS = 4
      P_QH = 5
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
      qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
      qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
      qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
      qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'

    CASE ( 7 ) ! Isobaric RUC2 GRIB file (AWIPS #211) from OSO
      nx_ext=93
      ny_ext=65
      nz_ext=37
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 8 ) ! Global Reanalysis on T62 Gaussian grid
      nx_ext = 192
      ny_ext = 94
      nz_ext = 28
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.
      IF(median_warn) THEN
        lon_0_360 = .FALSE.
        WRITE(6,'(4(a,/))')  &
        'The ARPS domain goes through Greenwich Meridian. The data should', &
        'be rearranged to handle this situation just as it has been done',  &
        'for GFS #3 grid. You can do it yourself or submit your request to',&
        'arpssupport@lists.ou.edu.'
        CALL arpsstop('Domain accross Greenwich Meridian.',1)
      ELSE
        lon_0_360 = .TRUE.
      END IF

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 9 ) ! RUC-2 (GEMPAK) from Doplight is 151x113x41
      nx_ext=151
      ny_ext=113
      nz_ext=41
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 10 )  ! ETA (GEMPAK #104, 80km) from NMC is 147x110x39
      nx_ext=147
      ny_ext=110
      nz_ext=39
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 11 ) ! Native RUC2 GRIB file (Grid #236, Grid #252 or Grid #130) from OSO
      CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,                &
                      extdinit,extdfcst,nx_ext,ny_ext,istatus)
      IF (istatus < 0) THEN
        WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',        &
                      'getgrbdims, ireturn = ',istatus,'.'
        CALL arpsstop('ERROR inside getgrbdims',1)
      END IF
      WRITE(6,'(1x,a,2(a,I5),a/)') 'Domain of native RUC2 data ',       &
                       'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'

      !nx_ext=151
      !ny_ext=113
      nz_ext=50
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      nscalar = 5
      nscalarq = 5
      P_QC = 1
      P_QR = 2
      P_QI = 3
      P_QS = 4
      P_QH = 5
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
      qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
      qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
      qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
      qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'

    CASE ( 12 ) ! Isobaric RUC2 GRIB file (Grid #236, Grid #252 or Grid #130) from OSO
      CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,                &
                      extdinit,extdfcst,nx_ext,ny_ext,istatus)
      IF (istatus < 0) THEN
        WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',        &
                      'getgrbdims, ireturn = ',istatus,'.'
        CALL arpsstop('ERROR inside getgrbdims',1)
      END IF
      WRITE(6,'(1x,a,2(a,I5),a/)') 'Domain of Isobaric RUC2 data ',       &
                       'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'

      !nx_ext=151
      !ny_ext=113
      nz_ext=37
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 13 )  ! AVN (GRIB #3, 1x1 degree) from NCEP is 360x181x26
      IF (nofixdim == 1) THEN  ! get dimension from data file
        CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,       &
                        extdinit,extdfcst,nx_ext,ny_ext,istatus)
        IF (istatus < 0) THEN
          WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',      &
                        'getgrbdims, ireturn = ',istatus,'.'
          CALL arpsstop('ERROR inside getgrbdims',1)
        END IF
        WRITE(6,'(a,2(a,I5),a/)') 'Subdomain of GFS global 1x1 data ',  &
                         'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'
      ELSE
        nx_ext = 360
        ny_ext = 181
      END IF

      nz_ext = 26
      nzsoil_ext = 4    ! 0, 0.1, 0.4, 1.0

      nstyp_ext  = 1

      IF(.NOT. median_warn .OR. nofixdim == 1) lon_0_360 = .TRUE.

      nscalar = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 14 )  ! AVN (GRIB #2, 2.5x2.5 degree) from NCEP is 144x73x26
      nx_ext = 144
      ny_ext = 73
      nz_ext = 26

      IF(.NOT. median_warn) lon_0_360 = .TRUE.

      nscalar = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 15 )  ! NCAR/NCEP Reanalysis-2 (GRIB #2, 2.5x2.5 degree) 144x73x16
      nx_ext = 144
      ny_ext = 73
      nz_ext = 17
      IF(.NOT. median_warn) lon_0_360 = .TRUE.

      nscalar = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 16 )  ! ETA (GRIB #218, 12km)

      IF (nofixdim == 1) THEN  ! get dimension from data file
        CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,              &
                        extdinit,extdfcst,nx_ext,ny_ext,istatus)
        IF (istatus < 0) THEN
          WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',      &
                        'getgrbdims, ireturn = ',istatus,'.'
          CALL arpsstop('ERROR inside getgrbdims',1)
        END IF
        IF (myproc == 0) WRITE(6,'(a,2(a,I5),a/)')                      &
          'Subdomain of NAM #218 12km data ',                           &
                         'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'

        npc = 1
        npr = 1
        ALLOCATE(domain_tile(npr,npc), STAT = istatus)
        CALL check_alloc_status(istatus, "ext2arps:domain_tile")
        domain_tile(1,1) = 1

      ELSE    ! find which tiles should be read

!        IF (latsw < lat_min(1)  .OR. latne > lat_max(50) .OR.           &
!            lonsw < lon_min(46) .OR. lonne > lon_max(54) ) THEN
        IF (latsw < lat_min(1)  .OR. latne > lat_max(50) .OR.           &
            latse < lat_min(1)  .OR. latnw > lat_max(50) .OR.           &
            lonsw < lon_min(46) .OR. lonne > lon_max(54) .OR.           &
            lonse < lon_min(46) .OR. lonnw > lon_max(54) ) THEN

          CALL arpsstop('ARPS domain is larger than ETA grid #218.',1)
        END IF

        latadj = 0.5           ! Values to be adjusted for complete coverage
        lonadj = 1.0           ! of the ARPS domain within the tiles to be read.
               ! Users can increase the values if ARPS grid is not complete
               ! covered within the external grid.
        !
        ! Locate ARPS south west corner
        !
        latpnt = latsw - latadj
        lonpnt = lonsw - lonadj

        k = 0                 ! tile no.
        SW_loop: DO j = 1,6   ! 6 rows of tile
          DO i = 1,9          ! 9 columns of tile
            k = k+1
            IF (latpnt > lat_min(k) .AND. latpnt < lat_max(k) .AND.     &
                lonpnt > lon_min(k) .AND. lonpnt < lon_max(k) ) EXIT SW_loop
          END DO
        END DO SW_loop
        ksw = k

        !
        ! Locate ARPS south east corner
        !
        latpnt = latse - latadj
        lonpnt = lonse + lonadj

        k = 0                 ! tile no.
        SE_loop: DO j = 1,6   ! 6 rows of tile
          DO i = 1,9          ! 9 columns of tile
            k = k+1
            IF (latpnt > lat_min(k) .AND. latpnt < lat_max(k) .AND.     &
                lonpnt > lon_min(k) .AND. lonpnt < lon_max(k) ) EXIT SE_loop
          END DO
        END DO SE_loop
        kse = k

        !
        ! Locate ARPS north east corner
        !
        latpnt = latne + latadj
        lonpnt = lonne + lonadj

        k = 0                 ! tile no.
        NE_loop: DO j = 1,6   ! 6 rows of tile
          DO i = 1,9          ! 9 columns of tile
            k = k+1
            IF (latpnt > lat_min(k) .AND. latpnt < lat_max(k) .AND.     &
                lonpnt > lon_min(k) .AND. lonpnt < lon_max(k) ) EXIT NE_loop
          END DO
        END DO NE_loop
        kne = k

        !
        ! Locate ARPS north west corner
        !
        latpnt = latnw + latadj
        lonpnt = lonnw - lonadj

        k = 0                 ! tile no.
        NW_loop: DO j = 1,6   ! 6 rows of tile
          DO i = 1,9          ! 9 columns of tile
            k = k+1
            IF (latpnt > lat_min(k) .AND. latpnt < lat_max(k) .AND.     &
                lonpnt > lon_min(k) .AND. lonpnt < lon_max(k) ) EXIT NW_loop
          END DO
        END DO NW_loop
        knw = k

!
!       Compute the i, j values of the four corner points.
!

        kswi = MOD(ksw-1,9)+1
        ksei = MOD(kse-1,9)+1
        knwi = MOD(knw-1,9)+1
        knei = MOD(kne-1,9)+1

        kswj = int((ksw+8)/9)
        ksej = int((kse+8)/9)
        knwj = int((knw+8)/9)
        knej = int((kne+8)/9)

!
!       Find the smallest "i" value, using SW and NW points.
!
        ktmp1 = kswi
        IF (knwi < ktmp1) ktmp1 = knwi

!
!       Find the smallest "j" value, using SW and SE points.
!
        ktmp2 = kswj
        IF (ksej < ktmp2) ktmp2 = ksej

        idiag = (ktmp2-1) * 9 + ktmp1

!
!       Find the largest "i" value, using NE and SE points.
!
        ktmp1 = knei
        IF (ksei > ktmp1) ktmp1 = ksei

!
!       Find the largest "j" value, using NE and NW points.
!
        ktmp2 = knej
        IF (knwj > ktmp2) ktmp2 = knwj

        jdiag = (ktmp2-1) * 9 + ktmp1

        IF (idiag > jdiag) THEN
           IF (myproc == 0)                                             &
             WRITE(6,*) 'Found ARPS south west corner in tile ',idiag,  &
                      ' But ARPS north east corner in tile ',jdiag
             CALL arpsstop("ETA218 problem",1)
        END IF

        jextmn = (idiag-1)/9
        jextmx = (jdiag-1)/9
        npc = jextmx - jextmn + 1         ! row number

        iextmn = MOD(idiag-1,9)
!        ! check if lonnw is within the NW tile (F. KONG)
!        k = (npc+jextmn-1)*9 + iextmn + 1
!        IF (lonnw < lon_min(k) .AND. iextmn > 0) iextmn=iextmn-1

        iextmx = MOD(jdiag-1,9)
        !IF(iextmx < 8) iextmx=iextmx+1
        npr = iextmx - iextmn + 1         ! column number
        ALLOCATE(domain_tile(npr,npc), STAT = istatus)
        CALL check_alloc_status(istatus, "ext2arps:domain_tile")
        DO j = 1,npc
          DO i = 1,npr
            domain_tile(i,j) = (j+jextmn-1)*9 + (iextmn+i-1) + 1
          END DO
        END DO

        IF ( MOD(domain_tile(npr,1),9) == 0) THEN
          nx_ext = (npr-1)*69 + 62
        ELSE
          nx_ext = npr*69
        END IF

        IF (domain_tile(1,npc) > 45) THEN
          ny_ext = (npc-1)*72 + 68
        ELSE
          ny_ext = npc*72
        END IF

        IF (myproc == 0) THEN
          WRITE(6,'(/,1x,a,/,1x,2(a,I6),a)')                            &
                'Will read in the following tiles from Grid #218 with', &
                ' (nx_ext x ny_ext) = (',nx_ext,' x ',ny_ext,'):'

          DO j = npc,1,-1
            DO i = 1,npr
              WRITE(6,FMT='(4x,I2.2)',ADVANCE='NO') domain_tile(i,j)
            END DO
            WRITE(6,*)
          END DO
          WRITE(6,*)
        END IF

      END IF    ! nofixdim

      nz_ext = 39
      nzsoil_ext = 5
      nstyp_ext = 1

      nscalar = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 17 )  !TINA, NARR regional reanalysis (GRIB #221, 32km) is 349x277x29
      !
      ! dynamic grid size does not work because of AWIP32.fixed is fixed.
      !
!      IF (nofixdim == 1) THEN  ! get dimension from data file
!        CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,       &
!                        extdinit,extdfcst,nx_ext,ny_ext,istatus)
!        IF (istatus < 0) THEN
!          WRITE(6,'(1x,2a,I2,a)') 'Error return from subroutine ',      &
!                        'getgrbdims, ireturn = ',istatus,'.'
!          CALL arpsstop('ERROR inside getgrbdims',1)
!        END IF
!        WRITE(6,'(a,2(a,I5),a/)') 'Subdomain of NARR Grid #221 32km data ',  &
!                         'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'
!      ELSE
        nx_ext = 349
        ny_ext = 277
!      END IF
      nz_ext = 29
      nzsoil_ext = 5
      nstyp_ext = 1

      nscalar = 2
      nscalarq = 2
      P_QC = 1
      P_QI = 2
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
      qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'

    CASE ( 18 ) ! Native 20-km RUC2 GRIB file (Grid #252)
      nx_ext=301
      ny_ext=225
      nz_ext=50
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      nscalar = 5
      nscalarq = 5
      P_QC = 1
      P_QR = 2
      P_QI = 3
      P_QS = 4
      P_QH = 5
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
      qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
      qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
      qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
      qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'

    CASE ( 19 ) ! Isobaric 20-km RUC2 GRIB file (Grid #252) from OSO
      nx_ext=301
      ny_ext=225
      nz_ext=37
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 20 ) ! 0.5 degree GFS GRIB2 file

      IF (nofixdim == 1) THEN  ! get dimension from data file
        CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,       &
                        extdinit,extdfcst,nx_ext,ny_ext,istatus)
        IF (istatus < 0) THEN
          WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',      &
                        'getgrbdims, ireturn = ',istatus,'.'
          CALL arpsstop('ERROR inside getgrbdims',1)
        END IF
        WRITE(6,'(a,2(a,I5),a/)') 'Subdomain of GFS global 0.5x0.5 data ',  &
                         'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'
      ELSE
        nx_ext = 720
        ny_ext = 361
      END IF

      nz_ext = 26
      nzsoil_ext = 4                 ! 0, 0.1, 0.4, 1

      nstyp_ext  = 1

      IF(.NOT. median_warn) lon_0_360 = .TRUE.

      nscalar = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 21 ) ! the lastest Rapid Refresh and High-resolution rapid refresh (HRRR) data sets.
                ! RUC native grid #130
      nx_ext = 451
      ny_ext = 337
      nz_ext = 50
      nzsoil_ext = 6  ! The dataset provides 5 levels, we add the surface layers as the first level.

      nscalar = 5
      nscalarq = 5
      P_QC = 1
      P_QR = 2
      P_QI = 3
      P_QS = 4
      P_QH = 5
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
      qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
      qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
      qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
      qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    CASE ( 22 ) ! ECMWF data for Shenzhen Project

      !nx_ext = 113
      !ny_ext =  81

      CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,                &
                      extdinit,extdfcst,nx_ext,ny_ext,istatus)
      IF (istatus < 0) THEN
        WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',        &
                      'getgrbdims, ireturn = ',istatus,'.'
        CALL arpsstop('ERROR inside getgrbdims',1)
      END IF
      WRITE(6,'(a,2(a,I5),a/)') 'Subdomain of ECMWF data ',             &
                       'nx_ext = ',nx_ext,' ny_ext = ',ny_ext,'.'

      nz_ext =  20
      nzsoil_ext = 5  ! The dataset provides 4 levels, we add the surface layers as the first level.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above
    CASE ( 23 ) ! WRF ARW history files

      CALL getwrfdims(dir_extd,extdname,extdopt,extdfmt,extdinit,extdfcst, &
                      nx_ext,ny_ext,nz_ext,nzsoil_ext,qnames_wrf,istatus)
      IF (istatus < 0) THEN
        WRITE(6,'(1x,2a,I0,a)') 'Error return from subroutine ',        &
                      'getwrfdims, ireturn = ',istatus,'.'
        CALL arpsstop('ERROR inside getgrbdims',1)
      END IF
      IF (myproc == 0) WRITE(6,'(1x,a,4(a,I5),/,27x,a,I5,a/)')          &
                       'External WRF domain size: ',                    &
                       'nx_ext = ',nx_ext,', ny_ext = ',ny_ext,         &
                     ', nz_ext = ',nz_ext,', nzsoil_ext = ',nzsoil_ext, &
                       'nscalar = ',nscalar,'.'

      nstyp_ext = 1
      nzsoil_extin = nzsoil_ext
      nzsoil_ext   = nzsoil_ext + 1

    CASE ( 24 )  ! NAM (GRIB #132, 16km) from NMC is 697x553x39

      IF (nofixdim == 1) THEN  ! get dimension from data file
        CALL getgrbdims(dir_extd,extdname,extdopt,extdfmt,              &
                        extdinit,extdfcst,nx_ext,ny_ext,istatus)
        IF (istatus < 0) THEN
          WRITE(6,'(1x,2a,I2,a)') 'Error return from subroutine ',      &
                        'getgrbdims, ireturn = ',istatus,'.'
          CALL arpsstop('ERROR inside getgrbdims',1)
        END IF
        WRITE(6,'(1x,a,2(a,I5),a/)') 'Subdomain of NAM Grid #132 16km data ',  &
                         'nx_ext = ',nx_ext,', ny_ext = ',ny_ext,'.'
      ELSE
        nx_ext = 697
        ny_ext = 553
      END IF

      IF (iboxe > iboxs) THEN
        nx_ext = iboxe-iboxs+1
      ELSE
        iboxs = 1
        iboxe = nx_ext
      END IF

      IF (jboxe > jboxs) THEN
        ny_ext = jboxe-jboxs+1
      ELSE
        jboxs = 1
        jboxe = ny_ext
      END IF

      nz_ext = 39
      nzsoil_ext = 5
      nstyp_ext = 1

      nscalar  = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE ( 50 ) ! Binary format lat/lon data
      nx_ext = 26
      ny_ext = 16
      nz_ext = 26
      nzsoil_ext = 2 ! FIXME: For now, just use two levels.

      ! No microphysical fields are read in, so nothing needs to be changed
      ! from above

    CASE ( 51 ) ! GFS (1x1 degree grid, South America)
      nx_ext = 118
      ny_ext = 91
      nz_ext = 26       ! NOTE: RH and CLWMR fields are 21 levels (up to 100 mb)
      nzsoil_ext = 3    ! ARPS: multi-layer OUSoil model

      nstyp_ext  = 1

      nscalar = 1
      nscalarq = 1
      P_QC = 1
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'

    CASE DEFAULT

      CALL arpsstop(                                                     &
              'extdopt option was invalid. Please specify a new option.',1)

  END SELECT

  IF (soilmodel_option == 1) nzsoil_ext = nzsoil
                          ! Using old ARPS Force-Restore Soil Model

  IF (myproc == 0) WRITE(6,'(1x,4(a,I5),/)')                             &
           'nx_ext = ',nx_ext,', ny_ext = ',ny_ext,', nz_ext = ',nz_ext, &
           ', nzsoil_ext = ',nzsoil_ext

!
!-----------------------------------------------------------------------
!
!  ALLOCATE scalar arrays based on dimension set above
!
!-----------------------------------------------------------------------
!
  IF(nscalar > 0) THEN
    ALLOCATE(qscalar(nx,ny,nz,nscalar),stat=istatus)
    CALL check_alloc_status(istatus, "ext2arps:QSCALAR")
    qscalar = 0.0

    ALLOCATE(qscalar_ext(nx_ext,ny_ext,nz_ext,nscalar),stat=istatus)
    CALL check_alloc_status(istatus, "ext2arps:qscalar_ext")
    qscalar_ext = 0.0
  END IF

!
!-----------------------------------------------------------------------
!
!  ALLOCATE and initialize external grid variables
!
!-----------------------------------------------------------------------
!
  ALLOCATE(x_ext(nx_ext),          STAT=istatus)
  ALLOCATE(y_ext(ny_ext),          STAT=istatus)
  ALLOCATE(xu_ext(nx_ext),         STAT=istatus)
  ALLOCATE(yu_ext(ny_ext),         STAT=istatus)
  ALLOCATE(xv_ext(nx_ext),         STAT=istatus)
  ALLOCATE(yv_ext(ny_ext),         STAT=istatus)
  ALLOCATE(lat_ext(nx_ext,ny_ext), STAT=istatus)
  ALLOCATE(lon_ext(nx_ext,ny_ext), STAT=istatus)
  ALLOCATE(latu_ext(nx_ext,ny_ext),STAT=istatus)
  ALLOCATE(lonu_ext(nx_ext,ny_ext),STAT=istatus)
  ALLOCATE(latv_ext(nx_ext,ny_ext),STAT=istatus)
  ALLOCATE(lonv_ext(nx_ext,ny_ext),STAT=istatus)
  ALLOCATE(zs_ext(nx,ny,nz_ext),   STAT=istatus)

  ALLOCATE(p_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:p_ext")
  ALLOCATE(hgt_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:hgt_ext")
  ALLOCATE(zp2_ext(nx,ny,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:zp2_ext")
  ALLOCATE(zp_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:zp_ext")
  ALLOCATE(zpsoil_ext(nx_ext,ny_ext,nzsoil_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:zpsoil_ext")

  ALLOCATE(t_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:t_ext")
  ALLOCATE(u_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:u_ext")
  ALLOCATE(v_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:v_ext")
  ALLOCATE(w_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:w_ext")
  ALLOCATE(qv_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:qv_ext")

  ALLOCATE(tsoil_ext (nx_ext,ny_ext,nzsoil_ext,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tsoil_ext")
  ALLOCATE(qsoil_ext (nx_ext,ny_ext,nzsoil_ext,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:qsoil_ext")
  ALLOCATE(wetcanp_ext(nx_ext,ny_ext,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:wetcanp_ext")
  ALLOCATE(snowdpth_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:snowdpth_ext")

  ALLOCATE(soiltyp_ext (nx_ext,ny_ext,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:soiltyp_ext")
  soiltyp_ext = 0
  ALLOCATE(stypfrct_ext(nx_ext,ny_ext,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:stypfrct_ext")
  ALLOCATE(vegtyp_ext (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:vegtyp_ext")
  ALLOCATE(lai_ext    (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:lai_ext")
  ALLOCATE(roufns_ext (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:roufns_ext")
  ALLOCATE(veg_ext    (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:veg_ext")

  ALLOCATE(trn_ext    (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:UBAR")
  ALLOCATE(psfc_ext   (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:UBAR")

  ALLOCATE(t_2m_ext (nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:t_2m_ext")
  ALLOCATE(qv_2m_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:qv_2m_ext")
  ALLOCATE(u_10m_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:u_10m_ext")
  ALLOCATE(v_10m_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:v_10m_ext")
  ALLOCATE(rain_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:rain_ext")

  ALLOCATE(dxfld(nx_ext),stat=istatus)
  ALLOCATE(dyfld(ny_ext),stat=istatus)
  ALLOCATE(rdxfld(nx_ext),stat=istatus)
  ALLOCATE(rdyfld(ny_ext),stat=istatus)
  ALLOCATE(dxfldu(nx_ext),stat=istatus)
  ALLOCATE(dyfldu(ny_ext),stat=istatus)
  ALLOCATE(rdxfldu(nx_ext),stat=istatus)
  ALLOCATE(rdyfldu(ny_ext),stat=istatus)
  ALLOCATE(dxfldv(nx_ext),stat=istatus)
  ALLOCATE(dyfldv(ny_ext),stat=istatus)
  ALLOCATE(rdxfldv(nx_ext),stat=istatus)
  ALLOCATE(rdyfldv(ny_ext),stat=istatus)

  ALLOCATE(rdzsoilfld(nx_ext,ny_ext,nzsoil_ext),STAT=istatus)
  CALL check_alloc_status(istatus, "ext2arps:rdzsoilfld")

  ALLOCATE(tem1_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem1_ext")
  ALLOCATE(tem2_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem2_ext")
  ALLOCATE(tem3_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem3_ext")
  ALLOCATE(tem4_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem4_ext")
  ALLOCATE(tem5_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:tem5_ext")

  ALLOCATE(xa_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:xa_ext")
  ALLOCATE(ya_ext(nx_ext,ny_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:ya_ext")
  ALLOCATE(avgzs_ext(nx,ny,nz_ext),stat=istatus)
  CALL check_alloc_status(istatus, "ext2arps:avgzs_ext")

  x_ext=0.0
  y_ext=0.0
  xu_ext=0.0
  yu_ext=0.0
  xv_ext=0.0
  yv_ext=0.0
  lat_ext=0.0
  lon_ext=0.0
  latu_ext=0.0
  lonu_ext=0.0
  latv_ext=0.0
  lonv_ext=0.0
  zs_ext=0.0

  p_ext=0.0
  hgt_ext=0.0
  zp2_ext=0.0
  zp_ext=0.0
  t_ext=0.0
  u_ext=0.0
  v_ext=0.0
  w_ext=0.0
  qv_ext=0.0

  tsoil_ext   =-999.0
  qsoil_ext   =-999.0
  wetcanp_ext =-999.0
  snowdpth_ext=-999.0

  soiltyp = -999.0
  stypfrct= -999.0
  vegtyp = -999
  lai    = -999.0
  roufns = -999.0
  veg    = -999.0

  trn_ext    =0.0
  psfc_ext   =0.0

  t_2m_ext =0.0
  qv_2m_ext=0.0
  u_10m_ext=0.0
  v_10m_ext=0.0
  rain_ext=0.0

  dxfld=0.0
  dyfld=0.0
  rdxfld=0.0
  rdyfld=0.0
  dxfldu=0.0
  dyfldu=0.0
  rdxfldu=0.0
  rdyfldu=0.0
  dxfldv=0.0
  dyfldv=0.0
  rdxfldv=0.0
  rdyfldv=0.0
  tem1_ext=0.0
  tem2_ext=0.0
  tem3_ext=0.0
  tem4_ext=0.0
  tem5_ext=0.0

  xa_ext=0.0
  ya_ext=0.0
  avgzs_ext=0.0

  IF (soilmodel_option == 1) THEN
    ALLOCATE(temsoil1(nx,ny,0:nstyps), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil1")
    ALLOCATE(temsoil2(nx,ny,0:nstyps), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil2")
    ALLOCATE(temsoil3(nx,ny,0:nstyps), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil3")
    ALLOCATE(temsoil4(nx,ny,0:nstyps), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil4")

    temsoil1(:,:,:) = 0.
    temsoil2(:,:,:) = 0.
    temsoil3(:,:,:) = 0.
    temsoil4(:,:,:) = 0.

    ALLOCATE(temsoil1_ext(nx_ext,ny_ext,0:nstyp_ext), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil1_ext")
    ALLOCATE(temsoil2_ext(nx_ext,ny_ext,0:nstyp_ext), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil2_ext")
    ALLOCATE(temsoil3_ext(nx_ext,ny_ext,0:nstyp_ext), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil3_ext")
    ALLOCATE(temsoil4_ext(nx_ext,ny_ext,0:nstyp_ext), stat = istatus)
    CALL check_alloc_status(istatus, "ext2arps:temsoil4_ext")

    temsoil1_ext(:,:,:) = 0.
    temsoil2_ext(:,:,:) = 0.
    temsoil3_ext(:,:,:) = 0.
    temsoil4_ext(:,:,:) = 0.
  END IF
!
!-----------------------------------------------------------------------
!
!  Loop through the data times provided via NAMELIST.
!
!-----------------------------------------------------------------------
!
  iniotfu = 21  ! FORTRAN unit number used for data output

  first_time = 1

  DO ifile=1,nextdfil
!
!-----------------------------------------------------------------------
!
!  Time conversions.
!  Formats:  extdtime='1994-05-06.18:00:00+000:00:00'
!            julfname='941261800'
!
!-----------------------------------------------------------------------
!
    READ(extdtime(ifile),'(a19,1x,a9)') extdinit,extdfcst
    IF(extdfcst == '         ') extdfcst='000:00:00'
    READ(extdinit,                                                      &
        '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)',ERR=920,END=920)           &
        iyr,imo,iday,ihr,imin,isec
    CALL julday(iyr,imo,iday,jldy)
    myr=MOD(iyr,100)
    ifhr=0
    ifmin=0
    ifsec=0
    IF (extdfcst(3:3) == ":") THEN
      READ(extdfcst,'(i2,1x,i2,1x,i2)',ERR=4,END=4) ifhr,ifmin,ifsec
    ELSE
      READ(extdfcst,'(i3,1x,i2,1x,i2)',ERR=4,END=4) ifhr,ifmin,ifsec
    END IF
    4  CONTINUE
    mfhr=MOD(ifhr,24)
    jldy = jldy + ifhr/24
    WRITE(julfname,'(i2.2,i3.3,i2.2,i2.2)') myr,jldy,ihr,mfhr
    CALL ctim2abss(iyr,imo,iday,ihr,imin,isec,iabssec)
    jabssec=(ifhr*3600) + (ifmin*60) + ifsec + iabssec
    kftime=jabssec-initsec
    curtim=REAL(kftime)
    IF (myproc == 0) WRITE(6,'(3a,a9,a,/25x,a,a19,a/a,a/a,i16,a,/,i27,a)') &
        ' Calling ',subnames(extdopt),' looking for ', extdfcst,        &
        ' hour forecast ', 'initialized at ',extdinit,' UTC',           &
        ' Julian filename: ',julfname,                                  &
        ' Which is  ',jabssec,' abs seconds or  ',kftime,               &
        ' seconds from the ARPS initial time.'

    ! name_nml_dir reused for output directory till to the end
    CALL get_output_dirname(1,dirname,curtim,1,name_nml_dir,istatus)
    ldirnam = LEN_TRIM(name_nml_dir)

!
!-----------------------------------------------------------------------
!
!  Call getextd to reads and converts data to ARPS units
!
!-----------------------------------------------------------------------
!

    select case ( extdopt )

    case ( 0 )  ! ARPS grid

      IF (.NOT. ALLOCATED(pbar_ext))  ALLOCATE(pbar_ext (nx_ext,ny_ext,nz_ext),stat=istatus)
      IF (.NOT. ALLOCATED(ptbar_ext)) ALLOCATE(ptbar_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
      IF (.NOT. ALLOCATED(qvbar_ext)) ALLOCATE(qvbar_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
      IF (.NOT. ALLOCATED(ubar_ext))  ALLOCATE(ubar_ext (nx_ext,ny_ext,nz_ext),stat=istatus)
      IF (.NOT. ALLOCATED(vbar_ext))  ALLOCATE(vbar_ext (nx_ext,ny_ext,nz_ext),stat=istatus)
      IF (.NOT. ALLOCATED(wbar_ext))  ALLOCATE(wbar_ext (nx_ext,ny_ext,nz_ext),stat=istatus)
      IF (ifile == 1) CALL check_alloc_status(istatus, "ext2arps:wbar_ext")

      extsfcopt = 0    ! near surface fields unavailable

      CALL getarps(nx_ext,ny_ext,nz_ext,nzsoil_ext,                       &
                   dir_extd,extdname,extdopt,extdfmt,                     &
                   extdinit,extdfcst,julfname,nstyps,                     &
                   iproj_ext,scale_ext,                                   &
                   trlon_ext,latnot_ext,x0_ext,y0_ext,                    &
                   lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,   &
                   p_ext,hgt_ext,zp_ext,zpsoil_ext,t_ext,qv_ext,          &
                   u_ext,tem4_ext,v_ext,tem3_ext,w_ext,                   &
                   qscalar_ext,                                           &
                   soiltyp_ext,stypfrct_ext,vegtyp_ext,                   &
                   lai_ext,roufns_ext,veg_ext,                            &
                   tsoil_ext,qsoil_ext,wetcanp_ext,                       &
                   snowdpth_ext,ubar_ext,vbar_ext,wbar_ext,               &
                   ptbar_ext,pbar_ext,qvbar_ext,                          &
                   tem1_ext,tem2_ext,istatus)
      ! tem3_ext holds uatv_ext and tem4_ext holds vatu_ext until
      ! u & v rotated to new map projection below (see calls to uvetomp)
      trn_ext(:,:) = zp_ext(:,:,2)

      stagger_ext = .TRUE.
!
!-----------------------------------------------------------------------
!
!  Get NMC RUC GRIB #87
!
!-----------------------------------------------------------------------
!
    case ( 1 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      ! 05/23/2002
      !
      ! ZUWEN getnmcruc only give out the average of tsoil and qsoil
      ! something should be done to interpolate tsoil and qsoil from
      ! the RUC model to the arps
      !
      CALL getnmcruc87(nx_ext,ny_ext,nz_ext,                              &
                       dir_extd,extdname,extdopt,extdfmt,                 &
                       extdinit,extdfcst,julfname,                        &
                       iproj_ext,scale_ext,                               &
                       trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                       lat_ext,lon_ext,                                   &
                       p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                       tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),             &
                       qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),wetcanp_ext, &
                       trn_ext,psfc_ext,                                  &
                       istatus)
!
!-----------------------------------------------------------------------
!
!  Get NMC ETA GRIB #212
!
!-----------------------------------------------------------------------
!
    case ( 2 )

      ! NOTE:  zpsoil_ext is defined as soil depth!!!
      CALL getnmceta212(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext,        &
                        dir_extd,extdname,extdopt,extdfmt,                &
                        extdinit,extdfcst,julfname,                       &
                        iproj_ext,scale_ext,                              &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,               &
                        lat_ext,lon_ext,zpsoil_ext,                       &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,           &
                        qscalar_ext(:,:,:,P_QC),                          &
                        tsoil_ext,qsoil_ext,wetcanp_ext,                  &
                        snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,        &
                        t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,  &
                        istatus)
      nstyp_ext = 1
      stypfrct_ext(:,:,1) = 1.
      stypfrct_ext(:,:,2:nstyps) = 0.
      zpsoil_ext(:,:,:) = -1.*zpsoil_ext(:,:,:)         ! EMK 15 June 2002

!      DO k = 1,nzsoil_ext
!        CALL a3dmax0(tsoil_ext(1,1,k,0),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,&
!               1,1,1,1,amax,amin)
!        IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')     &
!          'tsoil_ext_min= ', amin,', tsoil_ext_max=',amax
!      END DO
!
!-----------------------------------------------------------------------
!
!  Get LAPS data. Need LAPS library (see makearps for ext2arps.laps)
!
!-----------------------------------------------------------------------
!
    case ( 3 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getlaps(nx_ext,ny_ext,nz_ext,dir_extd,                         &
                   extdinit,extdfcst,julfname,iabssec,                    &
                   iproj_ext,scale_ext,                                   &
                   trlon_ext,latnot_ext,x0_ext,y0_ext,                    &
                   lat_ext,lon_ext,                                       &
                   p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                &
                   istatus)
!
!-----------------------------------------------------------------------
!
!  Get RUC/MAPS data in GEMPAK format with a special setup for GEMPAK
!  path and library (see makearps for ext2arps.gemruc).
!
!-----------------------------------------------------------------------
!
    case ( 4 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getgemruc(nx_ext,ny_ext,nz_ext,dir_extd,extdname,              &
                     extdinit,extdfcst,julfname,iabssec,                  &
                     iproj_ext,scale_ext,                                 &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                  &
                     lat_ext,lon_ext,                                     &
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,              &
                     istatus, tem1_ext)

    case ( 5 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getgemeta(nx_ext,ny_ext,nz_ext,dir_extd,extdname,              &
                     extdinit,extdfcst,julfname,iabssec,                  &
                     iproj_ext,scale_ext,                                 &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                  &
                     lat_ext,lon_ext,                                     &
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,              &
                     istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  Get COAMPS data
!
!-----------------------------------------------------------------------
!
    case ( 6 )

      extsfcopt = 1    ! Do support near surface field interpolation

      CALL getcoamps(nx_ext, ny_ext, nz_ext, nzsoil_ext,                &
                     dir_extd, extdinit,extdfcst,                       &
                     iproj_ext,scale_ext,                               &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                     lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,&
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                     qscalar_ext(:,:,:,P_QC),qscalar_ext(:,:,:,P_QR),   &
                     qscalar_ext(:,:,:,P_QI),qscalar_ext(:,:,:,P_QS),   &
                     qscalar_ext(:,:,:,P_QH),                           &
                     tsoil_ext,qsoil_ext,wetcanp_ext,                   &
                     snowdpth_ext,trn_ext,psfc_ext,                     &
                     t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,   &
                     istatus)
      stagger_ext = .TRUE.
!
!-----------------------------------------------------------------------
!
!  Get NMC RUC2 GRIB #211
!
!-----------------------------------------------------------------------
!
    case ( 7 )

      CALL getnmcruc211(nx_ext,ny_ext,nz_ext,                             &
                        dir_extd,extdname,extdopt,extdfmt,                &
                        extdinit,extdfcst,julfname,                       &
                        iproj_ext,scale_ext,                              &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,               &
                        lat_ext,lon_ext,                                  &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,           &
                        tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),            &
                        qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),wetcanp_ext,&
                        trn_ext,psfc_ext, t_2m_ext, qv_2m_ext,            &
                        u_10m_ext, v_10m_ext, istatus)
!
!-----------------------------------------------------------------------
!
!  Get NCEP global re-analysis on T62 Gaussian lat/lon grid
!
!-----------------------------------------------------------------------
!
    case ( 8 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getreanalt62(nx_ext,ny_ext,nz_ext,                             &
                        dir_extd,extdname,extdopt,extdfmt,                &
                        extdinit,extdfcst,julfname,                       &
                        iproj_ext,scale_ext,                              &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,               &
                        lat_ext,lon_ext,                                  &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,           &
                        tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),            &
                        qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),wetcanp_ext,&
                        trn_ext,psfc_ext, istatus)
!
!-----------------------------------------------------------------------
!
!  Get NMC RUC2 GEM
!
!-----------------------------------------------------------------------
!
    case ( 9 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getgemruc2(nx_ext,ny_ext,nz_ext,dir_extd,extdname,             &
                     extdinit,extdfcst,julfname,iabssec,                  &
                     iproj_ext,scale_ext,                                 &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                  &
                     lat_ext,lon_ext,                                     &
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,              &
                     istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  Get NMC ETA GEMPAK #104
!
!-----------------------------------------------------------------------
!
    case ( 10 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getgemeta2(nx_ext,ny_ext,nz_ext,dir_extd,extdname,             &
                      extdinit,extdfcst,julfname,iabssec,                 &
                      iproj_ext,scale_ext,                                &
                      trlon_ext,latnot_ext,x0_ext,y0_ext,                 &
                      lat_ext,lon_ext,                                    &
                      p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,             &
                      istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  Get NCEP RUC2 Native Coordinate GRIB (#236, or #252)
!
!  11 - RUC2 Grid #236
!  18 - 20-km RUC2 Grid #252
!
!  Note: if the option number for extdopt changed, should check the
!        subroutine getnmcrucn236 inside for "extdopt".
!
!-----------------------------------------------------------------------
!
    case ( 11, 18 )

      !extsfcopt = 0    ! Do not support near surface field interpolation

      CALL getnmcrucn236(nx_ext,ny_ext,nz_ext,                          &
                     dir_extd,extdname,extdopt,extdfmt,                 &
                     extdinit,extdfcst,julfname,                        &
                     iproj_ext,scale_ext,                               &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                     lat_ext,lon_ext,                                   &
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                     qscalar_ext(:,:,:,P_QC),qscalar_ext(:,:,:,P_QR),   &
                     qscalar_ext(:,:,:,P_QI),qscalar_ext(:,:,:,P_QS),   &
                     qscalar_ext(:,:,:,P_QH),                           &
                     tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),             &
                     qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),wetcanp_ext, &
                     trn_ext,psfc_ext,snowdpth_ext,t_2m_ext, qv_2m_ext, &
                     u_10m_ext, v_10m_ext, rain_ext, istatus)
!
!-----------------------------------------------------------------------
!
!  Get NCEP RUC2 Isobaric Coordinate GRIB (#236)
!
!-----------------------------------------------------------------------
!
    case ( 12, 19 )

      CALL getnmcrucp236(nx_ext,ny_ext,nz_ext,                          &
                     dir_extd,extdname,extdopt,extdfmt,                 &
                     extdinit,extdfcst,julfname,                        &
                     iproj_ext,scale_ext,                               &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                     lat_ext,lon_ext,                                   &
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                     tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),             &
                     qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),wetcanp_ext, &
                     trn_ext,psfc_ext, t_2m_ext, qv_2m_ext,             &
                     u_10m_ext, v_10m_ext, snowdpth_ext,                &
                     istatus)
!
!-----------------------------------------------------------------------
!
!  Get NCEP AVN GRIB #3
!
!-----------------------------------------------------------------------
!
    case ( 13 )

      !
      ! NOTE:  zpsoil_ext should be defined as soil depth!!!
      !
      CALL getncepavn3(nx_ext,ny_ext,nz_ext,nzsoil_ext,                 &
                       dir_extd,extdname,extdopt,extdfmt,               &
                       nofixdim,lon_0_360,extdinit,extdfcst,julfname,   &
                       iproj_ext,scale_ext,                             &
                       trlon_ext,latnot_ext,x0_ext,y0_ext,              &
                       lat_ext,lon_ext,zpsoil_ext,                      &
                       p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,          &
                       qscalar_ext(:,:,:,P_QC),                         &
                       tsoil_ext(:,:,:,0), qsoil_ext(:,:,:,0),          &
                       wetcanp_ext(:,:,0),                              &
                       snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext(:,:,1),&
                       t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,          &
                       rain_ext,istatus)

      nstyp_ext           = 1
      stypfrct_ext(:,:,1) = 1.0
      IF(nstyps > 1) THEN
        stypfrct_ext(:,:,2:nstyps) = 0.0
        soiltyp_ext (:,:,2:nstyps) = 0
      END IF

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

!-----------------------------------------------------------------------
!
!  Get NCEP AVN GRIB #2
!
!-----------------------------------------------------------------------
!
    case ( 14 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      ! NOTE:  zpsoil_ext should be defined as soil depth!!!
      CALL getncepavn2(nx_ext,ny_ext,nz_ext,nzsoil_ext,                   &
                       dir_extd,extdname,extdopt,extdfmt,lon_0_360,       &
                       extdinit,extdfcst,julfname,                        &
                       iproj_ext,scale_ext,                               &
                       trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                       lat_ext,lon_ext,zpsoil_ext,                        &
                       p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                       qscalar_ext(:,:,:,P_QC),                           &
                       tsoil_ext(:,:,:,0), qsoil_ext(:,:,:,0),            &
                       wetcanp_ext(:,:,0),                                &
                       snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext(:,:,1),  &
                       istatus)

      nstyp_ext           = 1
      stypfrct_ext(:,:,1) = 1.0
      IF(nstyps > 1) THEN
        stypfrct_ext(:,:,2:nstyps) = 0.0
        soiltyp_ext (:,:,2:nstyps) = 0
      END IF

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

!-----------------------------------------------------------------------
!
!  Get NCAR/NCEP Reanalysis-2 GRIB #2
!
!-----------------------------------------------------------------------
!
    case ( 15 )

      extsfcopt = 0    ! Do not support near surface field interpolation

      ! NOTE:  zpsoil_ext should be defined as soil depth!!!
      CALL getncepavn2(nx_ext,ny_ext,nz_ext,nzsoil_ext,                   &
                       dir_extd,extdname,extdopt,extdfmt,lon_0_360,       &
                       extdinit,extdfcst,julfname,                        &
                       iproj_ext,scale_ext,                               &
                       trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                       lat_ext,lon_ext,zpsoil_ext,                        &
                       p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                       qscalar_ext(:,:,:,P_QC),                           &
                       tsoil_ext(:,:,:,0), qsoil_ext(:,:,:,0),            &
                       wetcanp_ext(:,:,0),                                &
                       snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext(:,:,1),  &
                       istatus)

      nstyp_ext           = 1
      stypfrct_ext(:,:,1) = 1.0
      IF(nstyps > 1) THEN
        stypfrct_ext(:,:,2:nstyps) = 0.0
        soiltyp_ext (:,:,2:nstyps) = 0.0
      END IF

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)
!
!-----------------------------------------------------------------------
!
!  Get NMC ETA GRIB #218 (12km Eta model output)
!
!-----------------------------------------------------------------------
!
    case ( 16 )

      ! NOTE:  zpsoil_ext is defined as soil depth!!!
      CALL getnmceta218(nofixdim,nx_ext,ny_ext,nz_ext,nzsoil_ext,       &
                nstyp_ext,dir_extd,extdname,extdopt,extdfmt,            &
                extdinit,extdfcst,domain_tile,npr,npc,                  &
                iproj_ext,scale_ext,                                    &
                trlon_ext,latnot_ext,x0_ext,y0_ext,                     &
                lat_ext,lon_ext,zpsoil_ext,                             &
                p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                 &
                qscalar_ext(:,:,:,P_QC),                                &
                tsoil_ext(:,:,:,0),qsoil_ext(:,:,:,0),                  &
                wetcanp_ext(:,:,0),                                     &
                snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,              &
                t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext,rain_ext,     &
                istatus)

      nstyp_ext = 1
      stypfrct_ext(:,:,1) = 1.
      stypfrct_ext(:,:,2:nstyps) = 0.
      zpsoil_ext(:,:,:) = -1.*zpsoil_ext(:,:,:)

!      DO k = 1,nzsoil_ext
!        CALL a3dmax0(tsoil_ext(1,1,k,0),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,&
!               1,1,1,1,amax,amin)
!        WRITE(6,'(1x,2(a,e13.6))')     &
!          'tsoil_ext_min= ', amin,', tsoil_ext_max=',amax
!      END DO

!
!-----------------------------------------------------------------------
!TINA
!
!  Get NMC GRIB #221 - NARR (North American Regional Reanalysis) data
!
!-----------------------------------------------------------------------
!
    case ( 17 )

      ! NOTE:  zpsoil_ext is defined as soil depth!!!
      CALL getnarr221(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext,        &
                      dir_extd,extdname,extdopt,extdfmt,                &
                      extdinit,extdfcst,julfname,                       &
                      iproj_ext,scale_ext,                              &
                      trlon_ext,latnot_ext,x0_ext,y0_ext,               &
                      lat_ext,lon_ext,zpsoil_ext,                       &
                      p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,           &
                      qscalar_ext(:,:,:,P_QC),qscalar_ext(:,:,:,P_QI),  &
                      tsoil_ext,qsoil_ext,wetcanp_ext,                  &
                      snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,        &
                      t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,  &
                      istatus)
      nstyp_ext = 1
      stypfrct_ext(:,:,1) = 1.
      stypfrct_ext(:,:,2:nstyps) = 0.
      zpsoil_ext(:,:,:) = -1.*zpsoil_ext(:,:,:)

!      DO k = 1,nzsoil_ext
!        CALL a3dmax0(tsoil_ext(1,1,k,0),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,&
!               1,1,1,1,amax,amin)
!        WRITE(6,'(1x,2(a,e13.6))')     &
!          'tsoil_ext_min= ', amin,', tsoil_ext_max=',amax
!      END DO

!-----------------------------------------------------------------------
!
!  Get GFS 0.5 degree global data in GRIB2
!
!-----------------------------------------------------------------------
!
    case ( 20 )

      CALL getncepgfs(nx_ext,ny_ext,nz_ext,nzsoil_ext,                  &
                       dir_extd,extdname,extdopt,extdfmt,               &
                       nofixdim,lon_0_360,extdinit,extdfcst,julfname,   &
                       iproj_ext,scale_ext,                             &
                       trlon_ext,latnot_ext,x0_ext,y0_ext,              &
                       lat_ext,lon_ext,zpsoil_ext,                      &
                       p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,          &
                       qscalar_ext(:,:,:,P_QC),                         &
                       tsoil_ext(:,:,:,0), qsoil_ext(:,:,:,0),          &
                       wetcanp_ext(:,:,0),                              &
                       snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext(:,:,1),&
                      t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,  &
                      istatus)

      nstyp_ext           = 1
      stypfrct_ext(:,:,1) = 1.0
      IF(nstyps > 1) THEN
        stypfrct_ext(:,:,2:nstyps) = 0.0
        soiltyp_ext (:,:,2:nstyps) = 0
      END IF

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

!-----------------------------------------------------------------------
!
!  the latest Rapid Refresh and High-resolution rapid refresh (HRRR) data sets.
!
!-----------------------------------------------------------------------
!
    CASE ( 21 )

      CALL getnmcrucn130(nx_ext,ny_ext,nz_ext,nzsoil_ext,               &
                     dir_extd,extdname,extdopt,extdfmt,                 &
                     extdinit,extdfcst,julfname,                        &
                     iproj_ext,scale_ext,                               &
                     trlon_ext,latnot_ext,x0_ext,y0_ext,                &
                     lat_ext,lon_ext,zpsoil_ext,                        &
                     p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,            &
                     qscalar_ext(:,:,:,P_QC),qscalar_ext(:,:,:,P_QR),   &
                     qscalar_ext(:,:,:,P_QI),qscalar_ext(:,:,:,P_QS),   &
                     qscalar_ext(:,:,:,P_QH),                           &
                     tsoil_ext(:,:,:,0),qsoil_ext(:,:,:,0),wetcanp_ext, &
                     trn_ext,psfc_ext,snowdpth_ext, soiltyp_ext(:,:,1), &
                     t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,   &
                     istatus)
      nstyp_ext           = 1
      stypfrct_ext(:,:,1) = 1.0
      IF(nstyps > 1) THEN
        stypfrct_ext(:,:,2:nstyps) = 0.0
        soiltyp_ext (:,:,2:nstyps) = 0
      END IF

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

!-----------------------------------------------------------------------
!
!  ECMWF data set defined with GRIB table 128
!
!-----------------------------------------------------------------------
!
    CASE (22)

      extsfcopt = 0    ! Do not support near surface field interpolation for Shenzhen project

      ! NOTE:  zpsoil_ext is defined as soil depth!!!
      CALL getecmf128(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext,        &
                        dir_extd,extdname,extdopt,extdfmt,              &
                        extdinit,extdfcst,julfname,                     &
                        iproj_ext,scale_ext,                            &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,             &
                        lat_ext,lon_ext,zpsoil_ext,                     &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,         &
                       !qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,             &
                        tsoil_ext,qsoil_ext,wetcanp_ext,                &
                        snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,      &
                        t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,&
                        istatus)
      nstyp_ext = 1
      stypfrct_ext(:,:,1) = 1.
      stypfrct_ext(:,:,2:nstyps) = 0.
      zpsoil_ext(:,:,:) = -1.*zpsoil_ext(:,:,:)

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

!-----------------------------------------------------------------------
!
!  WRF ARW history files in netCDF format
!
!-----------------------------------------------------------------------
!
    CASE (23)

      ! NOTE:  zpsoil_ext is defined as soil depth!!!
      CALL getwrfdata(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext,nscalar,&
              dir_extd,extdname,extdopt,extdfmt,                        &
              extdinit,extdfcst,julfname,                               &
              soilmodel_option,nzsoil_extin,                            &
              iproj_ext,scale_ext,                                      &
              trlon_ext,latnot_ext,x0_ext,y0_ext,                       &
              lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,      &
              zpsoil_ext,                                               &
              p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                   &
              qscalar_ext,qnames_wrf,                                   &
              tsoil_ext(:,:,:,0),qsoil_ext(:,:,:,0),wetcanp_ext(:,:,0), &
              snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                &
              t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,          &
              istatus)

      stypfrct_ext(:,:,1)        = 1.
      stypfrct_ext(:,:,2:nstyps) = 0.
      zpsoil_ext(:,:,:) = zpsoil_ext(:,:,:)

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

      stagger_ext = .TRUE.

!-----------------------------------------------------------------------
!
!  NAM grid #132 (16 km resolution)
!
!-----------------------------------------------------------------------

    CASE (24)

      ! NOTE:  zpsoil_ext is defined as soil depth!!!
      CALL getnmceta132(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext,      &
                        dir_extd,extdname,extdopt,extdfmt,              &
                        extdinit,extdfcst,julfname,                     &
                        iboxs,iboxe,jboxs,jboxe,                        &
                        iproj_ext,scale_ext,                            &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,             &
                        lat_ext,lon_ext,zpsoil_ext,                     &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,         &
                        qscalar_ext(:,:,:,P_QC),                        &
                        tsoil_ext,qsoil_ext,wetcanp_ext,                &
                        snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,      &
                        t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,&
                        istatus)
      nstyp_ext = 1
      stypfrct_ext(:,:,1) = 1.
      stypfrct_ext(:,:,2:nstyps) = 0.
      zpsoil_ext(:,:,:) = -1.*zpsoil_ext(:,:,:)         ! EMK 15 June 2002

!-----------------------------------------------------------------------
!
!  Binary lat/lon
!
!-----------------------------------------------------------------------

    case ( 50 ) ! Binary format lat/lon data

      extsfcopt = 0    ! Do not support near surface field interpolation

      CALL get_avn_bin(nx_ext,ny_ext,nz_ext,extdopt,extdfmt,            &
                 dir_extd,extdname,extdinit,extdfcst,julfname,          &
                 iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext,&
                 lat_ext,lon_ext,p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,&
                 tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),                 &
                 qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),wetcanp_ext,     &
                 snowdpth_ext,trn_ext,psfc_ext,                         &
                 istatus)

    case ( 51 )

!-----------------------------------------------------------------------
!
!  Get NCEP GFS GRIB, 1 x 1 deg, South America
!
!-----------------------------------------------------------------------

      !
      ! NOTE:  zpsoil_ext should be defined as soil depth!!!
      !
      CALL getncepavn3(nx_ext,ny_ext,nz_ext,nzsoil_ext,                 &
                       dir_extd,extdname,extdopt,extdfmt,               &
                       nofixdim,lon_0_360,extdinit,extdfcst,julfname,   &
                       iproj_ext,scale_ext,                             &
                       trlon_ext,latnot_ext,x0_ext,y0_ext,              &
                       lat_ext,lon_ext,zpsoil_ext,                      &
                       p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,          &
                       qscalar_ext(:,:,:,P_QC),                         &
                       tsoil_ext(:,:,:,0), qsoil_ext(:,:,:,0),          &
                       wetcanp_ext(:,:,0),                              &
                       snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext(:,:,1),&
                       t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,          &
                       rain_ext,istatus)

      nstyp_ext           = 1
      stypfrct_ext(:,:,1) = 1.0
      IF(nstyps > 1) THEN
        stypfrct_ext(:,:,2:nstyps) = 0.0
        soiltyp_ext (:,:,2:nstyps) = 0.0
      END IF

      tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
      qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
      wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

    CASE DEFAULT

      CALL arpsstop( &
        'extdopt option was invalid. Please specify a new option.', 1 )

    END select

!-----------------------------------------------------------------------
!
!   If "istatus" has been set to -888, we had a bad file, but we want to go
!   ahead and still process more data.
!
!-----------------------------------------------------------------------

    IF(istatus == -888) CYCLE  ! file not foun

    IF(istatus /= 1) GO TO 999

    IF (myproc == 0) PRINT*,' '
    CALL a3dmax0(lat_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'lat_ext_min= ', amin,', lat_ext_max=',amax
    CALL a3dmax0(lon_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'lon_ext_min= ', amin,', lon_ext_max=',amax

    CALL a3dmax0(p_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'p_ext_min  = ', amin,', p_ext_max  =',amax
    CALL a3dmax0(hgt_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'hgt_ext_min= ', amin,', hgt_ext_max=',amax
    CALL a3dmax0(t_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          't_ext_min  = ', amin,', t_ext_max  =',amax
    CALL a3dmax0(u_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'u_ext_min  = ', amin,', u_ext_max  =',amax
    CALL a3dmax0(v_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'v_ext_min  = ', amin,', v_ext_max  =',amax
    CALL a3dmax0(w_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'w_ext_min  = ', amin,', w_ext_max  =',amax
    CALL a3dmax0(qv_ext ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,           &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'qv_ext_min = ', amin,', qv_ext_max =',amax

    DO nq=1,nscalar
      CALL a3dmax0(qscalar_ext(1,1,1,nq),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,   &
                 1,nz_ext,1,nz_ext,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,a,e13.6,2a,e13.6)')                           &
            TRIM(qnames(nq))//'_ext_min = ', amin,', ',TRIM(qnames(nq))//'_ext_max =',amax
    END DO

    !WDT note: FIXME add soiltype, etc,  plus add nstyp

    DO k = 1, nzsoil_ext
      CALL a3dmax0(tsoil_ext(1,1,k,0),                                  &
                   1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')                &
          'tsoil_',k,'_ext_min= ', amin,', tsoil_',k,'_ext_max=',amax

      CALL a3dmax0(qsoil_ext(1,1,k,0),                                  &
                   1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')                &
          'qsoil_',k,'_ext_min= ', amin,', qsoil_',k,'_ext_max=',amax
    END DO

    CALL a3dmax0(wetcanp_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,       &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'wetcanp_ext_min= ', amin,', wetcanp_ext_max=',amax
    CALL a3dmax0(snowdpth_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,      &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'snowd_ext_min= ', amin,', snow_ext_max=',amax

    CALL a3dmax0(trn_ext    ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,       &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'trn_ext_min= ', amin,', trn_ext_max=',amax
    CALL a3dmax0(psfc_ext   ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,       &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'psfc_ext_min= ', amin,', psfc_ext_max=',amax

    CALL a3dmax0(t_2m_ext ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'T_2m_ext_min= ', amin,', T_2m_ext_max=',amax
    CALL a3dmax0(qv_2m_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'qv_2m_ext_min= ', amin,', qv_2m_ext_max=',amax
    CALL a3dmax0(u_10m_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'U_10m_ext_min= ', amin,', U_10m_ext_max=',amax
    CALL a3dmax0(v_10m_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'V_10m_ext_min= ', amin,', V_10m_ext_max=',amax
    CALL a3dmax0(rain_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,          &
               1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
          'RAIN_ext_min= ', amin,', RAIN_ext_max=',amax
    IF (myproc == 0) PRINT*,' '
!
!-----------------------------------------------------------------------
!
!    First time through the time loop, calculate grid
!    transformation info.
!
!-----------------------------------------------------------------------
!
    IF(first_time == 1) THEN

      DO i=1,nx-1
        xscl(i)=0.5*(x(i)+x(i+1))
      END DO
      xscl(nx)=2.*xscl(nx-1) - xscl(nx-2)
      DO j=1,ny-1
        yscl(j)=0.5*(y(j)+y(j+1))
      END DO
      yscl(ny)=2.*yscl(ny-1) - yscl(ny-2)
!
!-----------------------------------------------------------------------
!
!  Find lat,lon locations of ARPS scalar, u and v grids.
!
!-----------------------------------------------------------------------
!
      CALL xytoll(nx,ny,xscl,yscl,tem1,tem2)
      CALL xytoll(nx,ny,   x,yscl,tem3,tem4)
      CALL xytoll(nx,ny,xscl,   y,tem5,tem6)
      CALL getmapr(iproj_arps,scale_arps,latnot_arps,                   &
                   trlon_arps,xorig_arps,yorig_arps)
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of external grid.
!
!-----------------------------------------------------------------------
!
      IF (trlon_ext > 180) trlon_ext = trlon_ext-360
      CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
      CALL setorig(1,x0_ext,y0_ext)
      DO j=1,ny_ext
        CALL lltoxy(1,1,lat_ext(1,j),lon_ext(1,j),x_ext(1),y_ext(j))
      END DO
      DO i=1,nx_ext
        CALL lltoxy(1,1,lat_ext(i,1),lon_ext(i,1),x_ext(i),y_ext(1))
      END DO
      IF ( stagger_ext ) THEN
        DO j=1,ny_ext
          CALL lltoxy(1,1,latu_ext(1,j),lonu_ext(1,j),xu_ext(1),yu_ext(j))
        END DO
        DO i=1,nx_ext
          CALL lltoxy(1,1,latu_ext(i,1),lonu_ext(i,1),xu_ext(i),yu_ext(1))
        END DO
        DO j=1,ny_ext
          CALL lltoxy(1,1,latv_ext(1,j),lonv_ext(1,j),xv_ext(1),yv_ext(j))
        END DO
        DO i=1,nx_ext
          CALL lltoxy(1,1,latv_ext(i,1),lonv_ext(i,1),xv_ext(i),yv_ext(1))
        END DO
      END IF

      IF ( lon_0_360 ) THEN
        DO i=1,nx_ext
          IF(x_ext(i) < 0.0)      x_ext(i) = x_ext(i)+360.
          IF(x_ext(i) < x_ext(1)) x_ext(i) = x_ext(i)+360.
                                   ! If Greenwich Meridian goes between
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS scalar grid in terms of the
!  external grid.
!
!-----------------------------------------------------------------------
!
      CALL lltoxy(nx,ny,tem1,tem2,xs2d,ys2d)

      IF ( lon_0_360 ) THEN
        DO j=1,ny
          DO i=1,nx
            IF( xs2d(i,j) < 0. )       xs2d(i,j)= xs2d(i,j) + 360.
            IF( xs2d(i,j) < xs2d(1,j)) THEN
              IF( nofixdim == 1)  THEN
                xs2d(i,j)= xs2d(i,j) + 360.
                                  ! If Greenwich Meridian goes between
              ELSE
                WRITE(6,'(2(2a,/))')    &
                'Wrong longitudes for ARPS grid points in terms of the ',&
                'external grid.','Interpolation may be wrong later if ', &
                'it is ignored.'
                CALL arpsstop('Wrong longitudes.',1)
              END IF
            END IF
          END DO
        END DO
      END IF

      CALL setijloc(nx,ny,nx_ext,ny_ext,xs2d,ys2d,                      &
                    x_ext,y_ext,iscl,jscl)
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS u grid in terms of the
!  external grid.
!
!-----------------------------------------------------------------------
!
      CALL lltoxy(nx,ny,tem3,tem4,xu2d,yu2d)
      IF( lon_0_360 ) THEN
        DO j=1,ny
          DO i=1,nx
            IF( xu2d(i,j) < 0. )       xu2d(i,j)= xu2d(i,j) + 360.
            IF( xu2d(i,j) < xu2d(1,j)) THEN
              IF (nofixdim == 1) THEN
                xu2d(i,j)= xu2d(i,j) + 360. ! extends beyond 360
                                            ! If Greenwich Meridian goes between
              ELSE
                WRITE(6,'(2(2a,/))')   &
                'Wrong longitudes for ARPS grid points in terms of the ',&
                'external grid.','Interpolation may be wrong later if ', &
                'it is ignored.'
                CALL arpsstop('Wrong longitudes.',1)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF ( stagger_ext ) THEN
        CALL setijloc(nx,ny,nx_ext,ny_ext,xu2d,yu2d,                    &
                      xu_ext,yu_ext,iu,ju)
      ELSE
        CALL setijloc(nx,ny,nx_ext,ny_ext,xu2d,yu2d,                    &
                      x_ext,y_ext,iu,ju)
      END IF
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS v grid in terms of the
!  external grid.
!
!-----------------------------------------------------------------------
!
      CALL lltoxy(nx,ny,tem5,tem6,xv2d,yv2d)
      IF ( lon_0_360 ) THEN
        DO j=1,ny
          DO i=1,nx
            IF( xv2d(i,j) < 0. )       xv2d(i,j)= xv2d(i,j) + 360.
            IF( xv2d(i,j) < xv2d(1,j)) THEN
              IF (nofixdim == 1) THEN
                xv2d(i,j)= xv2d(i,j) + 360.
                                   ! If Greenwich Meridian goes between
              ELSE
                WRITE(6,'(2(2a,/))')      &
                'Wrong longitudes for ARPS grid points in terms of the ',&
                'external grid.','Interpolation may be wrong later if ', &
                'it is ignored.'
                CALL arpsstop('Wrong longitudes.',1)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF ( stagger_ext ) THEN
        CALL setijloc(nx,ny,nx_ext,ny_ext,xv2d,yv2d,                    &
                      xv_ext,yv_ext,iv,jv)
      ELSE
        CALL setijloc(nx,ny,nx_ext,ny_ext,xv2d,yv2d,                    &
                      x_ext,y_ext,iv,jv)
      END IF

      xumin=xu2d(1,1)
      xumax=xu2d(1,1)
      xsmin=xs2d(1,1)
      xsmax=xs2d(1,1)
      DO j=1,ny-1
        xumin=MIN(xumin,xu2d(1 ,j))
        xumax=MAX(xumax,xu2d(nx,j))
        xsmin=MIN(xsmin,xs2d(1 ,j))
        xsmax=MAX(xsmax,xs2d(nx,j))
      END DO

      ysmin=ys2d(1,1)
      ysmax=ys2d(1,1)
      yvmin=yv2d(1,1)
      yvmax=yv2d(1,1)
      DO i=1,nx-1
        yvmin=MIN(yvmin,yv2d(i,1 ))
        yvmax=MAX(yvmax,yv2d(i,ny))
        ysmin=MIN(ysmin,ys2d(i,1 ))
        ysmax=MAX(ysmax,ys2d(i,ny))
      END DO

      IF (lon_0_360 .AND. nofixdim /= 1) THEN
        ! x_ext(1) for global grids (e.g., AVN grid 3) should be 0, not 360
        x_ext = MOD(x_ext, 360.0)
        WHERE (ABS(x_ext) < 1.0E-12) x_ext = 0.0
        xumin = MOD(xumin + 360.0, 360.0)
        xumax = MOD(xumax + 360.0, 360.0)
      END IF

      CALL mpmax0(xumax,xumin)
      CALL mpmax0(yvmax,yvmin)

      IF (myproc == 0) THEN
        WRITE (6,'(4(a,F15.3))') '  xu ARPS:',xumin,' - ',xumax,        &
                                  ' yv ARPS:',yvmin,' - ',yvmax
        WRITE (6,'(4(a,F15.3))') '  xu EXT: ',x_ext(1),' - ',x_ext(nx_ext), &
                                  ' yv EXT: ',y_ext(1),' - ',y_ext(ny_ext)

        IF(xumin < x_ext(1) .OR. xumax > x_ext(nx_ext) .OR.             &
           yvmin < y_ext(1) .OR. yvmax > y_ext(ny_ext)) THEN
          WRITE(6,'(/a/)') 'ext2arps aborting ....'
          CALL arpsstop('ARPS domain extends beyond the available external data.',1)
        END IF

      END IF
!
!-----------------------------------------------------------------------
!
!  Test code, for diagnostic testing.
!  Find x,y of Norman sounding in external grid.
!
!----------------------------------------------------------------------
!
      CALL lltoxy(1,1,latdiag,londiag,xdiag,ydiag)

      IF ( lon_0_360 .AND. xdiag < 0. ) xdiag= xdiag + 360.

      dmin=((xdiag-x_ext(1))*(xdiag-x_ext(1))+                          &
            (ydiag-y_ext(1))*(ydiag-y_ext(1)))

      idiag=1
      jdiag=1

      DO j=1,ny_ext
        DO i=1,nx_ext
          dd=((xdiag-x_ext(i))*(xdiag-x_ext(i))+                        &
              (ydiag-y_ext(j))*(ydiag-y_ext(j)))
          IF(dd < dmin) THEN
            dmin=dd
            idiag=i
            jdiag=j
          END IF
        END DO
      END DO
      IF (myproc == 0 ) THEN
      WRITE(6,'(a,f10.4,f10.4,/a,i5,i5,a,f10.4,f10.4)')                 &
            ' Nearest ext pt to diagnostic lat,lon: ',                  &
            latdiag,londiag,                                            &
            ' Diagnostic i,j: ',                                        &
            idiag,jdiag,', lat,lon = ',                                 &
            lat_ext(idiag,jdiag),lon_ext(idiag,jdiag)
      WRITE(6,'(///a,/2x,a)') ' External sounding at idiag,jdiag',      &
          'k   pres    hgt    temp   theta   dewp     u     v     dir    spd'
      END IF
!
!-----------------------------------------------------------------------
!
!  Convert units of external data and write as a sounding.
!
!  Note:  mpi results likely won't be same as non-mpi, as the mpi run
!         doesn't look at all processors.  This is harmless.
!-----------------------------------------------------------------------
!
      DO k=nz_ext,1,-1
        pmb=.01*p_ext(idiag,jdiag,k)
        tc=t_ext(idiag,jdiag,k)-273.15
        theta=t_ext(idiag,jdiag,k)*(p0/p_ext(idiag,jdiag,k))**rddcp
        IF( qv_ext(idiag,jdiag,k) > 0.) THEN
          smix=qv_ext(idiag,jdiag,k)/(1.-qv_ext(idiag,jdiag,k))
          e=(pmb*smix)/(0.62197 + smix)
          bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
          alge = ALOG(bige/6.112)
          tdc = (alge * 243.5) / (17.67 - alge)
        ELSE
          tdc = tc-30.
        END IF

        CALL uvrotdd(1,1,lon_ext(idiag,jdiag),                          &
                     u_ext(idiag,jdiag,k),v_ext(idiag,jdiag,k),         &
                     dir,spd)

        IF (myproc == 0)                                                &
          WRITE(6,'(i4,f6.0,f9.0,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1)')  &
                      k,pmb,                                            &
                      hgt_ext(idiag,jdiag,k),                           &
                      tc,theta,tdc,                                     &
                      u_ext(idiag,jdiag,k),                             &
                      v_ext(idiag,jdiag,k),                             &
                      dir,spd
      END DO
!
!-----------------------------------------------------------------------
!
!  Restore map projection to ARPS grid.
!
!-----------------------------------------------------------------------
!
      CALL setmapr(iproj_arps,scale_arps,latnot_arps,trlon_arps)

      CALL setorig(1,xorig_arps,yorig_arps)
!
!-----------------------------------------------------------------------
!
!  Find the min and max of the external indices that will be used
!  to be sure interpolation is within external dataset.
!
!-----------------------------------------------------------------------
!
      iextmn=iscl(1,1)
      jextmn=jscl(1,1)
      iextmx=iscl(1,1)
      jextmx=jscl(1,1)
      DO j=1,ny
        DO i=1,nx
          iextmn=MIN(iextmn,iscl(i,j))
          jextmn=MIN(jextmn,jscl(i,j))

          iextmx=MAX(iextmx,iscl(i,j))
          jextmx=MAX(jextmx,jscl(i,j))

          iextmn=MIN(iextmn,iu(i,j))
          jextmn=MIN(jextmn,ju(i,j))

          iextmx=MAX(iextmx,iu(i,j))
          jextmx=MAX(jextmx,ju(i,j))

          iextmn=MIN(iextmn,iv(i,j))
          jextmn=MIN(jextmn,jv(i,j))

          iextmx=MAX(iextmx,iv(i,j))
          jextmx=MAX(jextmx,jv(i,j))
        END DO
      END DO

    END IF   ! first_time

!
!-----------------------------------------------------------------------
!
!  External data comes in oriented to true north.
!  Change that to u,v in the new coordinate system.
!  (tem3_ext & tem4_ext from above hold uatv_ext & vatu_ext for use here)
!
!-----------------------------------------------------------------------
!
    IF ( stagger_ext ) THEN
      DO k=1,nz_ext
        CALL uvetomp(nx_ext,ny_ext,u_ext(1,1,k),tem4_ext(1,1,k),       &
                     lonu_ext,tem1_ext(1,1,k),tem2_ext(1,1,k))
      END DO
      u_ext = tem1_ext
      DO k=1,nz_ext
        CALL uvetomp(nx_ext,ny_ext,tem3_ext(1,1,k),v_ext(1,1,k),       &
                     lonv_ext,tem1_ext(1,1,k),tem2_ext(1,1,k))
      END DO
      v_ext = tem2_ext
    ELSE
      DO k=1,nz_ext
        CALL uvetomp(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),           &
                     lon_ext,tem1_ext(1,1,k),tem2_ext(1,1,k))
      END DO
      u_ext = tem1_ext
      v_ext = tem2_ext
      ! F.KONG add
      CALL uvetomp(nx_ext,ny_ext,u_10m_ext,v_10m_ext,                   &
                     lon_ext,tem1_ext(1,1,1),tem2_ext(1,1,1))
      DO j=1,ny_ext
        DO i=1,nx_ext
          u_10m_ext(i,j) = tem1_ext(i,j,1)
          v_10m_ext(i,j) = tem2_ext(i,j,1)
        END DO
      END DO
      IF (myproc == 0) THEN
        WRITE(*,'(1x,2(a,I4),2(a,F8.2))') 'idiag = ',i,', jdiag = ',j,  &
                      ', u_10m_ext = ',u_10m_ext(idiag,jdiag),          &
                      ', v_10m_ext = ',v_10m_ext(idiag,jdiag)
      END IF
      ! end F.KONG mod
    END IF
!
!-----------------------------------------------------------------------
!
!  Process height data
!
!  Interpolate the external heights horizontally to the
!  ARPS x,y.  This is for scalar pts.
!
!-----------------------------------------------------------------------
!
    CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,                       &
                 x_ext,y_ext,dxfld,dyfld,rdxfld,rdyfld)
    IF ( stagger_ext ) THEN
      CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,                     &
                   xu_ext,yu_ext,dxfldu,dyfldu,rdxfldu,rdyfldu)
      CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,                     &
                   xv_ext,yv_ext,dxfldv,dyfldv,rdxfldv,rdyfldv)
      IF ( extdopt == 0 ) THEN
        DO k=1,nz_ext
          CALL fldint2d(nx,ny,nx_ext,ny_ext,1,nx,1,ny,1,nx_ext,1,ny_ext,&
                      iorder,xs2d,ys2d,zp_ext(1,1,k),                   &
                      x_ext,y_ext,iscl,jscl,                            &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext,tem2_ext,tem3_ext,                       &
                      zp2_ext(1,1,k))
        END DO
      END IF
    END IF
    DO k=1,nz_ext
      CALL fldint2d(nx,ny,nx_ext,ny_ext,1,nx,1,ny,1,nx_ext,1,ny_ext,    &
                  iorder,xs2d,ys2d,hgt_ext(1,1,k),                      &
                  x_ext,y_ext,iscl,jscl,                                &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  tem1_ext,tem2_ext,tem3_ext,                           &
                  zs_ext(1,1,k))
    END DO

!
!-----------------------------------------------------------------------
!
!  Calculate x and y coordinates of external grid in ARPS coordinate
!  system.  Store them in xa_ext and ya_ext.
!
!-----------------------------------------------------------------------
!
    CALL lltoxy(nx_ext,ny_ext,lat_ext,lon_ext,xa_ext,ya_ext)

    IF(first_time == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Interpolate external terrain to arps grid (if exttrnopt=1).
!
!-----------------------------------------------------------------------
!
      IF (exttrnopt > 0) THEN ! set arps terrain to match external terrain

        CALL mkarps2d (nx_ext,ny_ext,nx,ny,                             &
                   iorder,iscl,jscl,x_ext,y_ext,                        &
                   xs2d,ys2d,trn_ext,htrn_tmp,                          &
                   dxfld,dyfld,rdxfld,rdyfld,                           &
                   tem1_ext(1,1,1),tem1_ext(1,1,2),                     &
                   tem1_ext(1,1,3),tem1_ext(1,1,4))

        IF (exttrnopt == 1) THEN
          hterain = htrn_tmp
        ELSE
          ntmergeinv = 1.d0/extntmrg
          DO j = 1,ny
            DO i = 1,nx
              idist = max(0,min(extntmrg,i-2,nx-2-i,j-2,ny-2-j))
              mfac = idist*ntmergeinv
              hterain(i,j) = (1.d0-mfac)*htrn_tmp(i,j)               &
                              + mfac*hterain(i,j)
            END DO
          END DO
        END IF

        ternopt = -1
        CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                 &
                    hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,        &
                    tem2(1,1,:),tem3(1,1,:),tem1)

        DO k=2,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              aj3z(i,j,k)=0.5*(j3(i,j,k)+j3(i,j,k-1))
            END DO
          END DO
        END DO

      END IF

!-----------------------------------------------------------------------
!
!  Write out terrain data
!
!-----------------------------------------------------------------------

      IF (terndmp > 0) THEN

        ternfn = runname(1:lfnkey)//".trndata"
        lternfn = lfnkey + 8

        temchar = ternfn
        ternfn = name_nml_dir(1:ldirnam)//temchar
        lternfn  = lternfn + ldirnam + 1

        CALL fnversn(ternfn, lternfn )

        PRINT *, 'Write terrain data to ',ternfn(1:lternfn)

        CALL edgfill(hterain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)

        CALL writtrn(nx,ny,ternfn(1:lternfn), dx,dy,                    &
                   mapproj,trulat1,trulat2,trulon,sclfct,               &
                   ctrlat,ctrlon,hterain)

        IF (terndmp == 1)                                               &
            CALL trncntl(nx,ny,ternfn(1:lternfn), x,y)

      END IF
!
!-----------------------------------------------------------------------
!
!  Set up z grid at scalar vertical levels.
!
!-----------------------------------------------------------------------
!
      onvf = 0

      CALL avgz(zp, onvf, nx,ny,nz, 1,nx, 1,ny, 1,nz-1, zs)
!
      DO j=1,ny-1
        DO i=1,nx-1
          zs(i,j,nz)=2.*zs(i,j,nz-1) - zs(i,j,nz-2)
        END DO
      END DO

      CALL edgfill(zs,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL edgfill(zp,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

    END IF
!
!-----------------------------------------------------------------------
!
!  Convert 2-m temperature and humidity and 10-m wind to arps grid
!
!-----------------------------------------------------------------------
!
!    IF ( extdopt == 2 .OR. extdopt == 13 ) THEN
!                       ! ETA (GRIB #212, 40km) or AVN (#3)
    IF (extsfcopt /= 0 .OR. i2dfmt > 0) THEN

      CALL mkarps2d (nx_ext,ny_ext,nx,ny,                               &
                     iorder,iscl,jscl,x_ext,y_ext,                      &
                     xs2d,ys2d,t_2m_ext,t_2m,                           &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext(1,1,1),tem1_ext(1,1,2),                   &
                     tem1_ext(1,1,3),tem1_ext(1,1,4))

      CALL mkarps2d (nx_ext,ny_ext,nx,ny,                               &
                     iorder,iscl,jscl,x_ext,y_ext,                      &
                     xs2d,ys2d,qv_2m_ext,qv_2m,                         &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext(1,1,1),tem1_ext(1,1,2),                   &
                     tem1_ext(1,1,3),tem1_ext(1,1,4))

      CALL mkarps2d (nx_ext,ny_ext,nx,ny,                               &
                     iorder,iscl,jscl,x_ext,y_ext,                      &
                     xs2d,ys2d,u_10m_ext,u_10m,                         &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext(1,1,1),tem1_ext(1,1,2),                   &
                     tem1_ext(1,1,3),tem1_ext(1,1,4))

      CALL mkarps2d (nx_ext,ny_ext,nx,ny,                               &
                     iorder,iscl,jscl,x_ext,y_ext,                      &
                     xs2d,ys2d,v_10m_ext,v_10m,                         &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext(1,1,1),tem1_ext(1,1,2),                   &
                     tem1_ext(1,1,3),tem1_ext(1,1,4))

      IF (exttrnopt == 0) &
        CALL mkarps2d (nx_ext,ny_ext,nx,ny,                             &
                   iorder,iscl,jscl,x_ext,y_ext,                        &
                   xs2d,ys2d,trn_ext,htrn_tmp,                          &
                   dxfld,dyfld,rdxfld,rdyfld,                           &
                   tem1_ext(1,1,1),tem1_ext(1,1,2),                     &
                   tem1_ext(1,1,3),tem1_ext(1,1,4))

    END IF

    htrn_tmp = amax1(htrn_tmp, 0.0)

    CALL mkarps2d (nx_ext,ny_ext,nx,ny,                                 &
                   iorder,iscl,jscl,x_ext,y_ext,                        &
                   xs2d,ys2d,rain_ext,rain,                             &
                   dxfld,dyfld,rdxfld,rdyfld,                           &
                   tem1_ext(1,1,1),tem1_ext(1,1,2),                     &
                   tem1_ext(1,1,3),tem1_ext(1,1,4))

    DO j=1,ny
      DO i=1,nx
        rain(i,j) = max(rain(i,j), 0.0)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Data transformations
!  Calculate log of pressure for external grid.
!  Change qv to RHstar   RHStar=sqrt(1.-relative humidity)
!
!-----------------------------------------------------------------------
!
    qvmin=999.
    qvmax=-999.
    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          qvmax=AMAX1(qvmax,qv_ext(i,j,k))
          qvmin=AMIN1(qvmin,qv_ext(i,j,k))
          tem5_ext(i,j,k)=ALOG(p_ext(i,j,k))
          qvsat=f_qvsat( p_ext(i,j,k), t_ext(i,j,k) )
          qv_ext(i,j,k)=SQRT(AMAX1(0.,(rhmax-(qv_ext(i,j,k)/qvsat))))
!        u_ext(i,j,k)=tem1_ext(i,j,k)
!        v_ext(i,j,k)=tem2_ext(i,j,k)
        END DO
      END DO
    END DO
    IF (myproc == 0) THEN
      PRINT *, ' qv_ext min = ',qvmin
      PRINT *, ' qv_ext max = ',qvmax
    END IF
!
!-----------------------------------------------------------------------
!
!  Calculate base-state sounding (vertical profile)
!
!-----------------------------------------------------------------------
!
   CALL extmnsnd(nx,ny,nx_ext,ny_ext,nz_ext,lvlprof,1,                  &
                iextmn,iextmx,jextmn,jextmx,1,nz_ext,                   &
                xscl,yscl,xa_ext,ya_ext,                                &
                hgt_ext,tem5_ext,t_ext,                                 &
                qv_ext,u_ext,v_ext,                                     &
                zsnd,plsnd,psnd,tsnd,ptsnd,rhssnd,qvsnd,                &
                rhosnd,usnd,vsnd)

!
!-----------------------------------------------------------------------
!
!  Process Pressure data
!
!-----------------------------------------------------------------------
!
    iprtopt=1
    CALL mkarpsvlz(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,                 &
                   iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,          &
                   hgt_ext,zs_ext,xs2d,ys2d,zs,p_ext,                     &
                   zsnd,plsnd,pbar,pprt,                                  &
                   dxfld,dyfld,rdxfld,rdyfld,                             &
                   tem1_ext,tem2_ext,tem3_ext,                            &
                   tem4_ext)

    IF (nsmooth > 0 .AND. mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    DO ksmth=1,nsmooth
      DO k=1,nz
        CALL smooth9p(pprt(1,1,k), nx,ny, 1,nx,1,ny,0, tem1)
      END DO

      IF (nsmooth > ksmth .AND. mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc,wbc,0,tem1)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc,sbc,0,tem1)
      END IF

    END DO
!
!-----------------------------------------------------------------------
!
!  Process potential temperature data
!
!-----------------------------------------------------------------------
!
    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          tem5_ext(i,j,k)=t_ext(i,j,k)*((p0/p_ext(i,j,k))**rddcp)
        END DO
      END DO
    END DO

    IF ( extsfcopt /= 0 .AND. t_2m_ext(1,1) > -900. ) THEN
                                             ! F.KONG add t_2m ajustment
      DO j=1,ny_ext
        DO i=1,nx_ext
          tem1_ext(i,j,1)=t_2m_ext(i,j)*((p0/psfc_ext(i,j))**rddcp)
        END DO
      END DO

      iprtopt=1
      CALL mkarpsvar1(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,            &
                     iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,      &
                     hgt_ext,zs_ext,xs2d,ys2d,zs,tem5_ext,              &
                     zsnd,ptsnd,ptbar,ptprt,                            &
                     trn_ext,tem1_ext(1,1,1),2.0,                       &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)
    ELSE

      iprtopt=1
      CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
                     iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,      &
                     hgt_ext,zs_ext,xs2d,ys2d,zs,tem5_ext,              &
                     zsnd,ptsnd,ptbar,ptprt,                            &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)

    END IF

    IF (nsmooth > 0 .AND. mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(ptprt,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(ptprt,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,                         &
                    smfct1,zs,ptprt,tem1,ptprt)
    END DO
!
!-----------------------------------------------------------------------
!
!  Process qv data.
!
!-----------------------------------------------------------------------
!
    IF ( extsfcopt /= 0 .AND.                                           &
         (t_2m_ext(1,1) > -900. .AND. qv_2m_ext(1,1) > -1.0) ) THEN
                                 ! Using 2m surface field interpolation

      DO j=1,ny_ext
        DO i=1,nx_ext
          qvsat=f_qvsat(psfc_ext(i,j)-25.0, t_2m_ext(i,j))
                                 ! -25.0 reflects 2m effect
          tem1_ext(i,j,2)=SQRT(AMAX1(0.,(rhmax-(qv_2m_ext(i,j)/qvsat))))
        END DO
      END DO

      iprtopt=0
      CALL mkarpsvar1(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,            &
                     iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,      &
                     hgt_ext,zs_ext,xs2d,ys2d,zs,qv_ext,                &
                     zsnd,rhssnd,qvbar,qv,                              &
                     trn_ext,tem1_ext(1,1,2),2.0,                       &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)

    ELSE

      iprtopt=0
      CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
                     iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,      &
                     hgt_ext,zs_ext,xs2d,ys2d,zs,qv_ext,                &
                     zsnd,rhssnd,qvbar,qv,                              &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)

    END IF

    IF (nsmooth > 0 .AND. mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(qv,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(qv,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,                         &
                    smfct1,zs,   qv,tem1,   qv)
    END DO
!
!-----------------------------------------------------------------------
!
!  Convert rhstar back to qv for writing, further calculations
!
!-----------------------------------------------------------------------
!
    qvmax=-999.
    qvmin=999.
    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          qvsat=f_qvsat( p_ext(i,j,k), t_ext(i,j,k) )
          rh=AMAX1(0.01,rhmax-(qv_ext(i,j,k)*qv_ext(i,j,k)))
          qv_ext(i,j,k)=rh*qvsat
          qvmax=AMAX1(qvmax,qv_ext(i,j,k))
          qvmin=AMIN1(qvmin,qv_ext(i,j,k))
        END DO
      END DO
    END DO
    IF ( myproc == 0 ) THEN
      PRINT *, ' qv_ext min = ',qvmin
      PRINT *, ' qv_ext max = ',qvmax
    END IF

    qvmax=-999.
    qvmin=999.
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          pres = pbar(i,j,k)+pprt(i,j,k)
          temp = (ptbar(i,j,k)+ptprt(i,j,k))*((pres/p0)**rddcp)
          qvsat=f_qvsat( pres, temp )
          IF(qvsat > 1.) THEN
            PRINT *, ' qvsat trouble at: ',i,j,k
            PRINT *, ' temp,press,qvsat: ',temp,pres,qvsat
          END IF
          rh=AMAX1(0.01,rhmax-(qv(i,j,k)*qv(i,j,k)))
          rh=AMIN1(rh,rhmax)
          qv(i,j,k)=rh*qvsat
          IF(qv(i,j,k) > 0.5) THEN
            PRINT *, ' qvsat trouble at: ',i,j,k
            PRINT *, ' temp,press,qvsat: ',temp,pres,qv(i,j,k)
            PRINT *, ' rh= ',rh
          END IF
          qvmax=AMAX1(qvmax,qv(i,j,k))
          qvmin=AMIN1(qvmin,qv(i,j,k))
          pres = pbar(i,j,k)
          temp = ptbar(i,j,k)*((pres/p0)**rddcp)
          qvsat=f_qvsat( pres, temp )
          rh=AMAX1(0.01,rhmax-(qvbar(i,j,k)*qvbar(i,j,k)))
          rh=AMIN1(rh,rhmax)
          qvbar(i,j,k)=rh*qvsat
        END DO
      END DO
    END DO
    CALL mpmax0(qvmax,qvmin)
    IF ( myproc == 0 ) THEN
      PRINT *, ' qv     min = ',qvmin
      PRINT *, ' qv     max = ',qvmax
    END IF
!
!-----------------------------------------------------------------------
!
!  Process qc data.
!  If first data point is less than or equal to -1, then
!  no qc info is provided...set to zero.
!
!-----------------------------------------------------------------------
!
    DO nq = 1,nscalar

      IF(qscalar_ext(1,1,1,nq) > -1.0) THEN
        iprtopt=0
        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,               &
                       iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,        &
                       hgt_ext,zs_ext,xs2d,ys2d,zs,qscalar_ext(:,:,:,nq),   &
                       zsnd,dumsnd,tem1,qscalar(:,:,:,nq),                  &
                       dxfld,dyfld,rdxfld,rdyfld,                           &
                       tem1_ext,tem2_ext,tem3_ext,                          &
                       tem4_ext)

        IF (nsmooth > 0 .AND. mp_opt > 0) THEN
          CALL acct_interrupt(mp_acct)
          CALL mpsendrecv2dew(qscalar(1,1,1,nq),nx,ny,nz,ebc,wbc,0,tem1)
          CALL mpsendrecv2dns(qscalar(1,1,1,nq),nx,ny,nz,nbc,sbc,0,tem1)
        END IF

        DO ksmth=1,nsmooth
          CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,                         &
                        smfct1,zs,qscalar(1,1,1,nq),tem1,qscalar(1,1,1,nq))
        END DO

        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              qscalar(i,j,k,nq)=MAX(qscalar(i,j,k,nq),0.0)
            END DO
          END DO
        END DO

      END IF
    END DO

!
!-----------------------------------------------------------------------
!
!  Process density which has been stuffed into tem2_ext array.
!  We really only need rhobar, so set iprtopt=1 and pass
!  a temporary array in place of rhoprt.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          tv_ext=t_ext(i,j,k) /                                         &
                 ((1.0+qv_ext(i,j,k))*                                  &
                  (1.0-(qv_ext(i,j,k)/(rddrv+qv_ext(i,j,k)))))
          tem5_ext(i,j,k)=p_ext(i,j,k)/(rd*tv_ext)
        END DO
      END DO
    END DO
!
    iprtopt=1
    CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,               &
                   iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,        &
                   hgt_ext,zs_ext,xs2d,ys2d,zs,tem5_ext,                &
                   zsnd,rhosnd,rhobar,tem1,                             &
                   dxfld,dyfld,rdxfld,rdyfld,                           &
                   tem1_ext,tem2_ext,tem3_ext,                          &
                   tem4_ext)

    IF (nsmooth > 0 .AND. mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(rhobar,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(rhobar,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,                         &
                    smfct1,zs,rhobar,tem1,rhobar)
    END DO
!
!-----------------------------------------------------------------------
!
!  Calculate and store the sound wave speed squared in csndsq.
!
!-----------------------------------------------------------------------
!
    IF(csopt == 1) THEN       ! Original definition of sound speed.
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            csndsq(i,j,k)= cpdcv*pbar(i,j,k)/rhobar(i,j,k)
          END DO
        END DO
      END DO

    ELSE IF(csopt == 2) THEN   ! Original sound speed multiplied
                               ! by a factor
      csconst = csfactr**2*cpdcv
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            csndsq(i,j,k)= csconst * pbar(i,j,k)/rhobar(i,j,k)
          END DO
        END DO
      END DO
    ELSE                      ! Specified constant sound speed.
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            csndsq(i,j,k)= csound * csound
          END DO
        END DO
      END DO
    END IF
!
    IF(hydradj == 1 .OR. hydradj == 2) THEN
      pconst=0.5*g*dz
!
!-----------------------------------------------------------------------
!
!  Create thermal buoyancy at each scalar point,
!  which is stored in tem2
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            qvprt=qv(i,j,k)-qvbar(i,j,k)
!            qtot=qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qh(i,j,k)
            qtot=0.0
            DO nq=1,nscalarq
             qtot=qtot+qscalar(i,j,k,nq)
            END DO
            tem2(i,j,k)=j3(i,j,k)*rhobar(i,j,k)*                          &
                      g*( (ptprt(i,j,k)/ptbar(i,j,k)) +                   &
                          (qvprt/(rddrv+qvbar(i,j,k))) -                  &
                          ((qvprt+qtot)/(1.+qvbar(i,j,k))) )
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Average thermal buoyancy to w points
!  which is stored in tem3
!
!-----------------------------------------------------------------------
!
      CALL avgsw(tem2,nx,ny,nz,1,nx,1,ny, tem3)

      IF(hydradj == 1) THEN

        DO j=1,ny
          DO i=1,nx
            tem1(i,j,1)=pprt(i,j,1)
            DO k=2,nz-2
              tem1(i,j,k)=                                                &
                    ( (1. - (pconst*j3(i,j,k-1)/csndsq(i,j,k-1)) )*       &
                        tem1(i,j,k-1) +                                   &
                        dz*tem3(i,j,k) ) /                                &
                        (1. + (pconst*j3(i,j,k)/csndsq(i,j,k)) )
              IF(j == 26 .AND. MOD(i,10) == 0) THEN
                IF(k == nz-2) PRINT *, '            Point i= ',i, '  j=26'
                PRINT *, ' k,pprt,tem1,diff =',                           &
                    k,pprt(i,j,k),tem1(i,j,k),                            &
                    (tem1(i,j,k)-pprt(i,j,k))
              END IF
              pprt(i,j,k)=tem1(i,j,k)
            END DO
            pprt(i,j,nz-1)=pprt(i,j,nz-2)
          END DO
        END DO

      ELSE IF(hydradj == 2) THEN

        DO j=1,ny
          DO i=1,nx
            tem1(i,j,nz-1)=pprt(i,j,nz-1)
            DO k=nz-2,2,-1
              tem1(i,j,k)=                                                &
                    ( (1.+ (pconst*j3(i,j,k+1)/csndsq(i,j,k+1)) )*        &
                        tem1(i,j,k+1) -                                   &
                        dz*tem3(i,j,k+1) ) /                              &
                        (1.- (pconst*j3(i,j,k  )/csndsq(i,j,k  )) )
              IF(j == 26 .AND. MOD(i,10) == 0) THEN
                IF(k == nz-2) PRINT *, '            Point i= ',i, '  j=26'
                PRINT *, ' k,pprt,tem1,diff =',                           &
                    k,pprt(i,j,k),tem1(i,j,k),                            &
                    (tem1(i,j,k)-pprt(i,j,k))
              END IF
              pprt(i,j,k)=tem1(i,j,k)
            END DO
            pprt(i,j,1)=pprt(i,j,2)
          END DO
        END DO
      END IF
    ELSE IF (hydradj == 3) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate total pressure, store in tem1.
!  Calculate virtual temperature, store in tem2.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k) = pbar(i,j,k)+pprt(i,j,k)
            temp = (ptbar(i,j,k)+ptprt(i,j,k))*                           &
                   ((tem1(i,j,k)/p0)**rddcp)
            tem2(i,j,k) = temp*(1.0+rvdrd*qv(i,j,k))/(1.0+qv(i,j,k))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equation, begining at k=2
!
!-----------------------------------------------------------------------
!
      pconst=-g/rd
      DO k=3,nz-1
        DO j=1,ny
          DO i=1,nx
            tvbar=0.5*(tem2(i,j,k)+tem2(i,j,k-1))
            tem1(i,j,k)=tem1(i,j,k-1)*                                    &
                       EXP(pconst*(zs(i,j,k)-zs(i,j,k-1))/tvbar)
            pprt(i,j,k)=tem1(i,j,k)-pbar(i,j,k)
          END DO
        END DO
      END DO
      DO j=1,ny
        DO i=1,nx
          tvbar=0.5*(tem2(i,j,2)+tem2(i,j,1))
          tem1(i,j,1)=tem1(i,j,2)*                                        &
                     EXP(pconst*(zs(i,j,1)-zs(i,j,2))/tvbar)
          pprt(i,j,1)=tem1(i,j,1)-pbar(i,j,1)
        END DO
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
!    Process wind data
!
!-----------------------------------------------------------------------
!
    CALL avgsu(zs_ext,nx,ny,nz_ext, 1,ny, 1,nz_ext, avgzs_ext, tem3_ext)
    CALL avgsu(zs,    nx,ny,    nz, 1,ny, 1,nz,          tem2, tem3)

    IF (loc_x == 1) THEN   ! West boundary
      DO k=1,nz_ext
        DO j=1,ny
          avgzs_ext(1,j,k) = zs_ext(1,j,k)
        END DO
      END DO
      DO k=1,nz
        DO j=1,ny
          tem2(1,j,k) = zs(1,j,k)
        END DO
      END DO
    END IF

    IF (loc_x == nproc_x) THEN
      DO k=1,nz_ext
        DO j=1,ny
          avgzs_ext(nx,j,k) = zs_ext(nx-1,j,k)
        END DO
      END DO
      DO k=1,nz
        DO j=1,ny
          tem2(nx,j,k) = zs(nx-1,j,k)
        END DO
      END DO
    END IF

    IF ( stagger_ext ) THEN ! u is staggered as in arps

      CALL avgsu(hgt_ext,nx_ext,ny_ext,nz_ext,1,ny_ext,1,nz_ext,tem5_ext,  &
         tem3_ext)
      tem5_ext(1,1:ny_ext,1:nz_ext) = hgt_ext(1,1:ny_ext,1:nz_ext)
      tem5_ext(nx_ext,1:ny_ext,1:nz_ext) = hgt_ext(nx_ext-1,1:ny_ext,1:nz_ext)

      iprtopt=0
      CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,                 &
                     iorder,iprtopt,intropt,iu,ju,xu_ext,yu_ext,            &
                     tem5_ext,avgzs_ext,xu2d,yu2d,tem2,u_ext,               &
                     zsnd,usnd,ubar,u,                                      &
                     dxfldu,dyfldu,rdxfldu,rdyfldu,                         &
                     tem1_ext,tem2_ext,tem3_ext,                            &
                     tem4_ext)

    ELSE IF ( extsfcopt /=0 .AND. (u_10m_ext(1,1) > -900.0) ) THEN
                 ! u is on a scalar point,  Using 10m surface interpolation

      iprtopt=0
      CALL mkarpsvar1(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,                &
                     iorder,iprtopt,intropt,iu,ju,x_ext,y_ext,              &
                     hgt_ext,avgzs_ext,xu2d,yu2d,tem2,u_ext,                &
                     zsnd,usnd,ubar,u,                                      &
                     trn_ext,u_10m_ext,10.0,                                &
                     dxfld,dyfld,rdxfld,rdyfld,                             &
                     tem1_ext,tem2_ext,tem3_ext,                            &
                     tem4_ext)

    ELSE  ! u is on a scalar point

      iprtopt=0

      CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,                 &
                     iorder,iprtopt,intropt,iu,ju,x_ext,y_ext,              &
                     hgt_ext,avgzs_ext,xu2d,yu2d,tem2,u_ext,                &
                     zsnd,usnd,ubar,u,                                      &
                     dxfld,dyfld,rdxfld,rdyfld,                             &
                     tem1_ext,tem2_ext,tem3_ext,                            &
                     tem4_ext)
    END IF

    IF (mp_opt > 0) THEN  ! note the difference of condition with others
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(u,nx,ny,nz,ebc,wbc,1,tem1)
      CALL mpsendrecv2dns(u,nx,ny,nz,nbc,sbc,1,tem1)
    END IF

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,1,                           &
                  smfct1,tem2,   u,tem1,   u)
    END DO

    CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                             &
            'umin final= ', amin,',  umax final=',amax
!
!-----------------------------------------------------------------------
!
!  Process v component
!
!-----------------------------------------------------------------------
!
    CALL avgsv(zs_ext,nx,ny,nz_ext, 1,nx, 1,nz_ext, avgzs_ext, tem3_ext)
    CALL avgsv(zs,nx,ny,nz, 1,nx, 1,nz, tem2, tem3)

    IF (loc_y == 1) THEN  ! south boundary
      DO k=1,nz_ext
        DO i=1,nx
          avgzs_ext(i,1,k)=zs_ext(i,1,k)
        END DO
      END DO
      DO k=1,nz
        DO i=1,nx
          tem2(i,1,k)=zs(i,1,k)
        END DO
      END DO
    END IF

    IF (loc_y == nproc_y) THEN
      DO k=1,nz_ext
        DO i=1,nx
          avgzs_ext(i,ny,k)=zs_ext(i,ny-1,k)
        END DO
      END DO
      DO k=1,nz
        DO i=1,nx
          tem2(i,ny,k)=zs(i,ny-1,k)
        END DO
      END DO
    END IF

    IF ( stagger_ext ) THEN ! v is staggered as in arps

      CALL avgsv(hgt_ext,nx_ext,ny_ext,nz_ext,1,nx_ext,1,nz_ext,tem5_ext,  &
         tem3_ext)
      tem5_ext(1:nx_ext,1,1:nz_ext) = hgt_ext(1:nx_ext,1,1:nz_ext)
      tem5_ext(1:nx_ext,ny_ext,1:nz_ext) = hgt_ext(1:nx_ext,ny_ext-1,1:nz_ext)

      iprtopt=0
      CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
                     iorder,iprtopt,intropt,iv,jv,xv_ext,yv_ext,        &
                     tem5_ext,avgzs_ext,xv2d,yv2d,tem2,v_ext,           &
                     zsnd,vsnd,vbar,v,                                  &
                     dxfldv,dyfldv,rdxfldv,rdyfldv,                     &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)

    ELSE IF ( extsfcopt /=0 .AND. (v_10m_ext(1,1) > -900.0) ) THEN
                 ! v is on a scalar point, Using 10m surface itnerpolation

      iprtopt=0
      CALL mkarpsvar1(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,            &
                     iorder,iprtopt,intropt,iv,jv,x_ext,y_ext,          &
                     hgt_ext,avgzs_ext,xv2d,yv2d,tem2,v_ext,            &
                     zsnd,vsnd,vbar,v,                                  &
                     trn_ext,v_10m_ext,10.0,                            &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)

    ELSE  ! v is on a scalar point

      iprtopt=0
      CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
                     iorder,iprtopt,intropt,iv,jv,x_ext,y_ext,          &
                     hgt_ext,avgzs_ext,xv2d,yv2d,tem2,v_ext,            &
                     zsnd,vsnd,vbar,v,                                  &
                     dxfld,dyfld,rdxfld,rdyfld,                         &
                     tem1_ext,tem2_ext,tem3_ext,                        &
                     tem4_ext)

    END IF

    IF (mp_opt > 0) THEN      ! note the difference of condition with others
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(v,nx,ny,nz,ebc,wbc,2,tem1)
      CALL mpsendrecv2dns(v,nx,ny,nz,nbc,sbc,2,tem1)
    END IF

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,2,                         &
                  smfct1,tem2,   v,tem1,   v)
    END DO

    CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'vmin final= ', amin,',  vmax final=',amax
!
!-----------------------------------------------------------------------
!
!  Process 2D surface fields if required (sfcout = 1)
!
!-----------------------------------------------------------------------
!
    ! determine the number of soil types supplied
    IF (nstyp_ext <= 0) THEN
      DO is=1,nstyps
        IF (tsoil_ext(1,1,1,is) > 0) nstyp_ext = is
      END DO
    END IF

    ! This should never be the case (or else ext sfc arrays allocated too small!):
    nstyp_ext = min(nstyp_ext,nstyps)

    wetcout = 0
    tsoilout = 0
    qsoilout = 0
    snowdout = 0

!EMK NEW 6 June 2002
!TODO:  Make sure other external model types are processed correctly.
    IF (soilmodel_option == 1) THEN ! Will use old force-restore soil model.

      !commented out mkarps2d & smooth9's, add call to intrp_soil

      DO is=0,nstyp_ext
        IF ( tsoil_ext(1,1,1,is) > 0.0 ) THEN
          tsoilout = 1
!          CALL mkarps2d (nx_ext,ny_ext,nx,ny,                         &
!                         iorder,iscl,jscl,x_ext,y_ext,                  &
!                         xs2d,ys2d,tsoil_ext(1,1,1,is),tsoil(1,1,1,is), &
!                         dxfld,dyfld,rdxfld,rdyfld,                     &
!                         tem1_ext(1,1,1),tem1_ext(1,1,2),               &
!                         tem1_ext(1,1,3),tem1_ext(1,1,4))
!          DO ksmth=1,nsmooth
!            DO k=1,nzsoil
!              CALL smooth9p(tsoil(1,1,k,is), nx,ny, 1,nx,1,ny, tem1)
!            END DO
!          END DO
        END IF

        IF ( qsoil_ext(1,1,1,is) >= 0.0 ) THEN
          qsoilout = 1
!          CALL mkarps2d (nx_ext,ny_ext,nx,ny,                          &
!                         iorder,iscl,jscl,x_ext,y_ext,                 &
!                         xs2d,ys2d,qsoil_ext(1,1,1,is),qsoil(1,1,1,is),&
!                         dxfld,dyfld,rdxfld,rdyfld,                    &
!                         tem1_ext(1,1,1),tem1_ext(1,1,2),              &
!                         tem1_ext(1,1,3),tem1_ext(1,1,4))
!          DO ksmth=1,nsmooth
!            DO k=1,nzsoil
!              CALL smooth9p(qsoil(1,1,k,is), nx,ny, 1,nx,1,ny, tem1)
!            END DO
!          END DO
        END IF

        IF ( wetcanp_ext(1,1,is) >= 0.0 ) THEN
          wetcout = 1
!          CALL mkarps2d (nx_ext,ny_ext,nx,ny,                          &
!                         iorder,iscl,jscl,x_ext,y_ext,                 &
!                         xs2d,ys2d,wetcanp_ext(1,1,is),wetcanp(1,1,is),&
!                         dxfld,dyfld,rdxfld,rdyfld,                    &
!                         tem1_ext(1,1,1),tem1_ext(1,1,2),              &
!                         tem1_ext(1,1,3),tem1_ext(1,1,4))
!          DO ksmth=1,nsmooth
!            CALL smooth9p(wetcanp(1,1,is), nx,ny, 1,nx,1,ny, tem1)
!          END DO
        END IF
      END DO ! is

     ! DTD: Hmmm, it seems this should come AFTER intrp_soil!

!      CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp_ext,nstyp,tsoil,  &
!                          qsoil,wetcanp)

      xw = 0.
      yw = 0.
      DO j = 1,ny-1
        DO i = 1,nx-1
          xw(i,j) =  MAX(0.,MIN(1.,(x_ext(iscl(i,j)+1)-xs2d(i,j))       &
                      /(x_ext(iscl(i,j)+1)-x_ext(iscl(i,j)))))
          yw(i,j) =  MAX(0.,MIN(1.,(y_ext(jscl(i,j)+1)-ys2d(i,j))       &
                      /(y_ext(jscl(i,j)+1)-y_ext(jscl(i,j)))))
        END DO
      END DO

      DO is = 0,nstyp_ext

        DO j = 1,ny_ext
          DO i = 1,nx_ext
            temsoil1_ext(i,j,is) = tsoil_ext(i,j,1,is)
            temsoil2_ext(i,j,is) = tsoil_ext(i,j,2,is)
            temsoil3_ext(i,j,is) = qsoil_ext(i,j,1,is)
            temsoil4_ext(i,j,is) = qsoil_ext(i,j,2,is)
          END DO
        END DO
      END DO

      CALL intrp_soil(nx_ext,ny_ext,nx,ny,nstyp_ext,nstyp,xw,yw,iscl,jscl, &
                      temsoil1_ext,temsoil2_ext,                           &
                      temsoil3_ext,temsoil4_ext,wetcanp_ext,               &
                      soiltyp_ext,stypfrct_ext,vegtyp_ext,                 &
                      temsoil1,temsoil2,                                   &
                      temsoil3,temsoil4,wetcanp,                           &
                      soiltyp,stypfrct,vegtyp)

      DO is = 0,nstyps
        DO j = 1,ny-1
          DO i = 1,nx-1
            tsoil(i,j,1,is) = temsoil1(i,j,is)
            tsoil(i,j,2,is) = temsoil2(i,j,is)
            qsoil(i,j,1,is) = temsoil3(i,j,is)
            qsoil(i,j,2,is) = temsoil4(i,j,is)
            IF ( is > 0 ) THEN
              IF ( stypfrct(i,j,is) <= 0. ) THEN
                ! Eddie's change to avoid passing 0 values to sfcphysics in case
                ! that arps is to be run later starting from this data set.
                tsoil(i,j,1,is) = tsoil(i,j,1,0)
                tsoil(i,j,2,is) = tsoil(i,j,2,0)
                qsoil(i,j,1,is) = qsoil(i,j,1,0)
                qsoil(i,j,2,is) = qsoil(i,j,2,0)
              END IF
            END IF
          END DO
        END DO
      END DO

      ! DTD: placed call to fix_soil_nstyp here, after intrp_soil routine
      ! previously it was called before, before tsoil, qsoil arrays are filled
      ! using intrp_soil

      CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp_ext,nstyp,tsoil,  &
                          qsoil,wetcanp)
      ! end DTD


!WDT FIXME other data types ...; make soiltyp consistant with
! soil temps, etc.; move this section before  soil temp/moist section
!  INTEGER, allocatable :: soiltyp (:,:,:)! Soil type
!  REAL, allocatable ::    stypfrct(:,:,:)! Soil type fraction
!  INTEGER, allocatable :: vegtyp (:,:)   ! Vegetation type
!  REAL, allocatable ::    lai    (:,:)   ! Leaf Area Index
!  REAL, allocatable ::    roufns (:,:)   ! Surface roughness
!  REAL, allocatable ::    veg    (:,:)   ! Vegetation fraction
!
!      CALL fix_stypfrct_nstyp(nx,ny,nstyp_ext,nstyp,stypfrct)

!WDT FIXME
! Add capability to write out sfcdata file here? ...
!WDT end


!EMK END 18 June 2002

    ELSE ! New OUSoil Model

      IF (extdopt /= 0) THEN ! Not processing ARPS data

        !First convert zpsoil to soil depth.
        DO k = 1,nzsoil
          DO j = 1,ny-1
            DO i = 1,nx-1
              zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
            END DO
          END DO
          IF (lvldbg > 25) THEN
            WRITE(6,'(1x,a,I2,F12.5)') 'ARPS     soil depth are: k = ',k,zpsoil(2,2,k)
          END IF
        END DO

        DO k = 1,nzsoil_ext-1
          IF (lvldbg > 25) THEN
            WRITE(6,'(1x,a,I2,F12.5)') 'External soil depth are: k = ',k,zpsoil_ext(2,2,k)
          END IF
          DO j =1,ny_ext-1
            DO i =1,nx_ext-1
              rdzsoilfld(i,j,k) = 1.0/( zpsoil_ext(i,j,k+1) - zpsoil_ext(i,j,k) )
                         ! Fixed by WYH, reciprocal is expected in intrpsoil
            END DO
          END DO
        END DO
        IF (lvldbg > 25) THEN
          WRITE(6,'(1x,a,I2,F12.5)') 'External soil depth are: k = ',nzsoil_ext,zpsoil_ext(2,2,nzsoil_ext)
        END IF

        tsoilout = 1
        qsoilout = 1
        wetcout = 1

        CALL intrpsoil3d_avg(nx_ext,ny_ext,nzsoil_ext,vegtyp_ext,       &
                           tsoil_ext,qsoil_ext,wetcanp_ext,             &
                           x_ext,y_ext,zpsoil_ext,                      &
                           rdxfld,rdyfld,rdzsoilfld,                    &
                           nx,ny,nzsoil,nstyps,                         &
                           vegtyp,tsoil,qsoil,wetcanp,                  &
                           xs2d,ys2d,zpsoil,iscl,jscl,ksoil3d)

        !Convert zpsoil back from soil depth to MSL
        DO k = 1,nzsoil
          DO j = 1,ny-1
            DO i = 1,nx-1
              zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
            END DO
          END DO
        END DO

      ELSE ! Processing ARPS data

!EMK 15 June 2002

        ! Convert to soil depth
        DO k = 1,nzsoil_ext
          DO j =1,ny_ext
            DO i =1,nx_ext
              zpsoil_ext(i,j,k) = trn_ext(i,j) - zpsoil_ext(i,j,k)
            END DO
          END DO
        END DO
        DO k = 1,nzsoil
          DO j =1,ny
            DO i =1,nx
              zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
            END DO
          END DO
        END DO

        DO k = 1,nzsoil_ext-1
          DO j =1,ny_ext-1
            DO i =1,nx_ext-1
              rdzsoilfld(i,j,k) = 1.0/ ( zpsoil_ext(i,j,k+1) - zpsoil_ext(i,j,k) )
                         ! Fixed by WYH, reciprocal is expected in intrpsoil
            END DO
          END DO
        END DO

        tsoilout = 1
        qsoilout = 1
        wetcout = 1

        CALL intrpsoil3d_pst(nx_ext,ny_ext,nzsoil_ext,nstyp_ext, &
                           soiltyp_ext,stypfrct_ext,vegtyp_ext,  &
                           tsoil_ext,qsoil_ext,wetcanp_ext,      &
                           x_ext,y_ext,zpsoil_ext,               &
                           rdxfld,rdyfld,rdzsoilfld,             &
                           nx,ny,nzsoil,nstyp,soiltyp,stypfrct,  &
                           vegtyp,tsoil,qsoil,wetcanp,xs2d,ys2d, &
                           zpsoil,iscl,jscl,ksoil3d)

        ! Convert back to MSL m
        DO k = 1,nzsoil_ext
          DO j =1,ny_ext
            DO i =1,nx_ext
              zpsoil_ext(i,j,k) = trn_ext(i,j) - zpsoil_ext(i,j,k)
            END DO
          END DO
        END DO
        DO k = 1,nzsoil
          DO j =1,ny
            DO i =1,nx
              zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
            END DO
          END DO
        END DO

      END IF ! IF (exdtopt /= 0)
    END IF ! Force-restore or OUSoil

!EMK END 18 June 2002

!
! Commented out 6 June 2002 by Eric Kemp.  I'm not sure if
! we need this or not.
!    CALL initsoilext(nx,ny,nzsoil,nzsoil_ext,nstyp,zpsoil,      &
!         zpsoil_ext,tsoil,tsoil_ext,qsoil,qsoil_ext)


!EMK END 6 June 2002

    IF ( snowdpth_ext(1,1) >= 0 ) THEN

      snowdout = 1

!
!     CALL mkarps2d (nx_ext,ny_ext,nx,ny,              &
!                    iorder,iscl,jscl,x_ext,y_ext,     &
!                    xs2d,ys2d,snowdpth_ext,snowdpth,  &
!                    dxfld,dyfld,rdxfld,rdyfld,        &
!                    tem1_ext(1,1,1),tem1_ext(1,1,2),  &
!                    tem1_ext(1,1,3),tem1_ext(1,1,4))

      dx_ext = x_ext(2)-x_ext(1)
      dy_ext = y_ext(2)-y_ext(1)

      DO j=1,ny-1
        DO i=1,nx-1
!          dmin=((xs2d(i,j)-x_ext(1))*(xs2d(i,j)-x_ext(1)) + &
!                (ys2d(i,j)-y_ext(1))*(ys2d(i,j)-y_ext(1)))
!          isnow=1
!          jsnow=1
!          DO jj=jscl(i,j),MIN(jscl(i,j)+1,ny_ext)
!            DO ii=iscl(i,j),MIN(iscl(i,j)+1,nx_ext)
!              dd=((xs2d(i,j)-x_ext(ii))*(xs2d(i,j)-x_ext(ii)) + &
!                  (ys2d(i,j)-y_ext(jj))*(ys2d(i,j)-y_ext(jj)))
!              IF(dd < dmin) THEN
!                dmin=dd
!                isnow=ii
!                jsnow=jj
!              END IF
!            END DO
!          END DO

!          dx_ext = x_ext(iscl(i,j)+1)-x_ext(iscl(i,j))
!          dy_ext = y_ext(jscl(i,j)+1)-y_ext(jscl(i,j))

          isnow = MIN( NINT( (xs2d(i,j)-x_ext(iscl(i,j)))/dx_ext ) + iscl(i,j), nx_ext )
          jsnow = MIN( NINT( (ys2d(i,j)-y_ext(jscl(i,j)))/dy_ext ) + jscl(i,j), ny_ext )
          snowdpth(i,j) = snowdpth_ext(isnow,jsnow)
        END DO
      END DO

    END IF

    IF ( wetcout == 1 .OR. snowdout == 1 .OR.                           &
         tsoilout == 1 .OR. qsoilout == 1 ) THEN

      zpsoilout = 1
      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
      lfn = lfnkey + 9 + tmstrln

      !IF (mp_opt > 0 .AND. joindmp <= 0 ) THEN
      !  CALL gtsplitfn(runname(1:lfnkey)//'.soilvar.'//timsnd(1:tmstrln),&
      !                 1,1,loc_x,loc_y,1,1,                             &
      !                 0,0,0,lvldbg,soiloutfl,ireturn)
      !   lfn = LEN_TRIM(soiloutfl)
      !END IF

      temchar = soiloutfl
      soiloutfl = name_nml_dir(1:ldirnam)//temchar
      lfn  = lfn + ldirnam + 1

      IF (soildmp > 0) THEN

        !CALL fnversn(soiloutfl, lfn)

        !IF (myproc == 0) PRINT *, 'Write soil initial data to ',soiloutfl(1:lfn)

        ! blocking inserted for ordering i/o for message passing
        DO i=0,nprocs-1,dumpstride
          IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN
            IF(mp_opt > 0 .AND. joindmp(FINDX_S) > 0) THEN
              CALL wrtjoinsoil(nx,ny,nzsoil,nstyps, soiloutfl(1:lfn),   &
                   dx,dy,zpsoil,                                        &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   1,1,1,1,1,                                           &
                   tsoil,qsoil,wetcanp,snowdpth,soiltyp)
            ELSE
              CALL wrtsoil(nx,ny,nzsoil,nstyps, soiloutfl(1:lfn),       &
                   dx,dy,zpsoil,                                        &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   1,1,1,1,1,                                           &
                   tsoil,qsoil,wetcanp,snowdpth,soiltyp)
            END IF
          END IF
          IF (mp_opt > 0) CALL mpbarrier
        END DO

        IF (soildmp == 1) CALL soilcntl(nx,ny, nzsoil,zpsoil,soiloutfl, &
                    zpsoilout,tsoilout,qsoilout, wetcout,snowdout,x,y)


      END IF

    END IF
!
!-----------------------------------------------------------------------
!
!  Test code, for diagnostic testing.
!  Find x,y of diagnostic sounding location in ARPS grid.
!
!-----------------------------------------------------------------------
!
    CALL lltoxy(1,1,latdiag,londiag,xdiag,ydiag)

    dmin=((xdiag-xscl(1))*(xdiag-xscl(1))+                                &
          (ydiag-yscl(1))*(ydiag-yscl(1)))
    idiag_arps=1
    jdiag_arps=1

    DO j=2,ny-2
      DO i=2,nx-2
        dd=((xdiag-xscl(i))*(xdiag-xscl(i))+                              &
            (ydiag-yscl(j))*(ydiag-yscl(j)))
        IF(dd < dmin) THEN
          dmin=dd
          idiag_arps=i
          jdiag_arps=j
        END IF
      END DO
    END DO
    CALL xytoll(1,1,xscl(idiag_arps),yscl(jdiag_arps),                    &
                latd,lond)
    IF (myproc == 0) WRITE(6,'(a,f10.4,f10.4,/a,i5,i5,a,f10.4,f10.4)')    &
          ' Nearest ARPS pt to diagnostic lat,lon: ',                     &
          latdiag,londiag,                                                &
          ' Diagnostic i,j: ',                                            &
          idiag_arps,jdiag_arps,' lat,lon= ',latd,lond
    IF (myproc == 0) WRITE(6,'(///a,/2x,a)')                              &
        ' ARPS extracted sounding at idiag,jdiag',                        &
        'k   pres   hgt   temp   theta   dewp     u     v     dir    spd'
!
!-----------------------------------------------------------------------
!
!  Convert units of ARPS data and write as a sounding.
!
!-----------------------------------------------------------------------
!
    DO k=nz-2,1,-1
      ppasc=pbar(idiag_arps,jdiag_arps,k)+pprt(idiag_arps,jdiag_arps,k)
      pmb=.01*(pbar(idiag_arps,jdiag_arps,k)+pprt(idiag_arps,jdiag_arps,k))
      theta=ptbar(idiag_arps,jdiag_arps,k)+ptprt(idiag_arps,jdiag_arps,k)
      tc=(theta*((ppasc/p0)**rddcp))-273.15
      IF( qv(idiag_arps,jdiag_arps,k) > 0.) THEN
        smix=qv(idiag_arps,jdiag_arps,k)/(1.-qv(idiag_arps,jdiag_arps,k))
        e=(pmb*smix)/(0.62197 + smix)
        bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
        alge = ALOG(bige/6.112)
        tdc = (alge * 243.5) / (17.67 - alge)
      ELSE
        tdc = tc-30.
      END IF

      CALL uvrotdd(1,1,londiag,                                           &
                   u(idiag_arps,jdiag_arps,k),                            &
                   v(idiag_arps,jdiag_arps,k),                            &
                   dir,spd)

      IF (myproc == 0)                                                    &
        WRITE(6,'(i4,f6.0,f9.0,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1)')      &
              k,pmb,                                                      &
              zs(idiag_arps,jdiag_arps,k),                                &
              tc,theta,tdc,                                               &
              u(idiag_arps,jdiag_arps,k),                                 &
              v(idiag_arps,jdiag_arps,k),                                 &
              dir,spd
    END DO

!
!-----------------------------------------------------------------------
!
!  Find vertical velocity and make any u-v adjustments
!  according to wndadj option.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rhostr(i,j,k)=ABS(j3(i,j,k))*rhobar(i,j,k)
        END DO
      END DO
    END DO
    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(rhostr,nx,ny,nz,4,4,0,tem1)
      CALL mpsendrecv2dns(rhostr,nx,ny,nz,4,4,0,tem1)
    END IF


    IF(wndadj == 0) THEN

      IF ( extdopt == 0 ) THEN  ! ARPS history data

        iprtopt=0
        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
                       iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,      &
                       zp_ext,zp2_ext,xs2d,ys2d,zp,w_ext,                 &
                       zsnd,dumsnd,wbar,w,                                &
                       dxfld,dyfld,rdxfld,rdyfld,                         &
                       tem1_ext,tem2_ext,tem3_ext,                        &
                       tem4_ext)

        IF (nsmooth > 0 .AND. mp_opt > 0) THEN
          CALL acct_interrupt(mp_acct)
          CALL mpsendrecv2dew(w,nx,ny,nz,ebc,wbc,3,tem1)
          CALL mpsendrecv2dns(w,nx,ny,nz,nbc,sbc,3,tem1)
        END IF

        DO ksmth=1,nsmooth
          CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,3,                     &
                        smfct1,zp,w,tem1,w)
        END DO

      ELSE

        w = 0.

      END IF

    ELSE

      CALL adjuvw( nx,ny,nz,u,v,w,wcont,ptprt,ptbar,                      &
                   zp,j1,j2,j3,aj3z,mapfct,rhostr,tem1,                   &
                   wndadj,obropt,obrzero,0,                               &
                   tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

    END IF
!
!-----------------------------------------------------------------------
!
!  Enforce vertical w boundary conditions
!
!-----------------------------------------------------------------------
!
    IF (ext_lbc == 1 .or. ext_vbc == 1)                                   &
       CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

    IF (ext_vbc == 1) CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,             &
              rhostr,tem1,tem2,tem3,j1,j2,j3)
!
!
!-----------------------------------------------------------------------
!
!  Assign zero-gradient horizontal boundary conditions
!  to the horizontal & vertical winds.
!
!-----------------------------------------------------------------------
!
    IF (ext_lbc == 1) THEN
      CALL lbcw(nx,ny,nz,1.0, w,wcont,tem1,tem2,tem3,tem4,3,3,3,3,      &
                3,3,3,3)
      CALL bcu(nx,ny,nz,1.0, u, tem1,tem2,tem3,tem4, 3,3,3,3,1,1,       &
                3,3,3,3)
      CALL bcv(nx,ny,nz,1.0, v, tem1,tem2,tem3,tem4, 3,3,3,3,1,1,       &
                3,3,3,3)
    END IF
!
!-----------------------------------------------------------------------
!
!  Zero out uninitialized fields
!
!-----------------------------------------------------------------------
!
     tem1=0.
!
!-----------------------------------------------------------------------
!
!  Print out the max/min of output variables.
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) WRITE(6,'(/1x,a/)')                                &
        'Min and max of External data interpolated to the ARPS grid:'

    CALL a3dmax0(x,1,nx,1,nx,1,1,1,1,1,1,1,1, amax,amin)
    IF (myproc == 0) WRITE(6,'(/1x,2(a,e13.6))')                        &
            'xmin    = ', amin,',  xmax    =',amax

    CALL a3dmax0(y,1,ny,1,ny,1,1,1,1,1,1,1,1, amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'ymin    = ', amin,',  ymax    =',amax

    CALL a3dmax0(z,1,nz,1,nz,1,1,1,1,1,1,1,1, amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'zmin    = ', amin,',  zmax    =',amax

    CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,                  &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'zpmin   = ', amin,', zpmax    =',amax

    CALL a3dmax0(rhobar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,            &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'rhobarmin=', amin,', rhobarmax=',amax

    CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,                   &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'umin    = ', amin,',  umax    =',amax

    CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,                   &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'vmin    = ', amin,',  vmax    =',amax

    CALL a3dmax0(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,                   &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'wmin    = ', amin,',  wmax    =',amax

    CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,             &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'ptbarmin= ', amin,',  ptbarmax=',amax

    CALL a3dmax0(ptprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,             &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'ptprtmin= ', amin,',  ptprtmax=',amax

    CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,              &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'pbarmin= ', amin,',  pbarmax =',amax

    CALL a3dmax0(pprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,              &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'pprtmin = ', amin,',  pprtmax =',amax

    CALL a3dmax0(qv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                &
              amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'qvmin   = ', amin,',  qvmax   =',amax

    DO nq=1,nscalar
      CALL a3dmax0(qscalar(1,1,1,nq),1,nx,1,nx,1,ny,1,ny,   &
                 1,nz,1,nz,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,a,e13.6,2a,e13.6)')                           &
            TRIM(qnames(nq))//'_min = ', amin,', ',TRIM(qnames(nq))//'_max =',amax
    END DO

    CALL a3dmax0(rain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
            'RAINmin   = ', amin,',  RAINmax   =',amax

    IF(sfcout == 1) THEN

      IF (soilmodel_option == 1) THEN

      CALL a3dmax0(tsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil,     &
                amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'tsoilmin = ', amin,',  tsoilmax =',amax
      CALL a3dmax0(qsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil,     &
                amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'qsoilmin= ', amin,',  qsoilmax=',amax
      CALL a3dmax0(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'wetcanpmin = ', amin,',  wetcanpmax =',amax

      ELSE IF (soilmodel_option == 2) THEN

      DO k=1,nzsoil
      CALL a3dmax0(tsoil(1,1,k,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,        &
               1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')                &
              'tsoil_',k,'_min= ', amin,', tsoil_',k,'_max=',amax
      CALL a3dmax0(qsoil(1,1,k,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,        &
               1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')                &
              'qsoil_',k,'_min = ', amin,',  qsoil_',k,'_max =',amax
      END DO
      CALL a3dmax0(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'wetcanpmin = ', amin,',  wetcanpmax =',amax
      END IF

    END IF

    IF(i2dfmt > 0) THEN
      CALL a3dmax0(t_2m,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              't_2m_min = ', amin,',  t_2m_max =',amax
      CALL a3dmax0(qv_2m,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'qv_2m_min = ', amin,',  qv_2m_max =',amax
      CALL a3dmax0(u_10m,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'u_10m_min = ', amin,',  u_10m_max =',amax
      CALL a3dmax0(v_10m,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                       &
              'v_10m_min = ', amin,',  v_10m_max =',amax
    END IF
!
!-----------------------------------------------------------------------
!
!  Print the mean sounding that was used in setting the
!  mean ARPS variables.
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) WRITE(6,'(a/2a)')                                  &
        '  Mean sounding for ARPS base state variables',                &
        '  k     p(mb)     z(m)    pt(mb)    T(C)   qv(kg/kg) ',        &
        '  RH %  u(m/s)    v(m/s)'

    DO k=lvlprof,1,-50
      pres = psnd(k)
      temp = ptsnd(k)*((pres/p0)**rddcp)
      rh=AMAX1(0.01,rhmax-(rhssnd(k)*rhssnd(k)))
      qvsat=f_qvsat( pres, temp )
      qvval=rh*qvsat
      IF (myproc == 0)                                                  &
        WRITE(6,'(i4,f9.1,f9.1,f9.1,f9.1,f9.5,f9.1,f9.1,f9.1)')         &
          k,(0.01*psnd(k)),zsnd(k),ptsnd(k),(temp-273.15),              &
          qvval,(100.*rh),usnd(k),vsnd(k)
    END DO
!
!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!  First find a unique name basdmpfn(1:lbasdmpf) for the grid and
!  base state array dump file
!
!-----------------------------------------------------------------------
!
    tem1=0.

    IF ( hdmpfmt == 5 ) THEN
      IF (myproc == 0) WRITE (6,'(a/a)')                                &
          'Program ext2arps does not support Savi3D format.',           &
          'Reset to binary format, hdmpfmt=1.'
      hdmpfmt = 1
    END IF

    IF ( hdmpfmt == 9 ) GO TO 450

    CALL gtbasfn(runname(1:lfnkey),name_nml_dir,ldirnam,hdmpfmt,        &
                 mgrid,nestgrd,basdmpfn,lbasdmpf)

    IF (myproc == 0)                                                    &
      PRINT*,'Writing base state history dump ',basdmpfn(1:lbasdmpf)

    grdbas =1
    mstout =1
!    rainout=0
    prcout =0
    trbout =0

    hdmpgrdfmt = hdmpfmt
    IF ( grdbasopt == 0 .OR. (grdbasopt == 1 .AND. ifile > 1) ) hdmpgrdfmt = 0
                        ! only want the first file to dump the base state?

    IF (nstyp_ext < 1) landout = 0

    !  blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride

      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

        CALL dtadump(nx,ny,nz,nzsoil,nstyp,                             &
               hdmpgrdfmt,iniotfu,basdmpfn(1:lbasdmpf),grdbas,filcmprs, &
               u,v,w,ptprt,pprt,qv,qscalar,                             &
               tem1,tem1,tem1,                                          &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                  &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               rain,tem1,tem1,                                          &
               tem1,tem1,tem1,tem1,tem1,                                &
               tem1,tem1,tem1,tem1,                                     &
               tem2,tem3,tem4)

      END IF

      IF (mp_opt > 0) CALL mpbarrier

    END DO

    450 CONTINUE

    IF (myproc == 0) PRINT*,'curtim = ', curtim

    CALL gtdmpfn(runname(1:lfnkey),name_nml_dir,ldirnam,curtim,hdmpfmt, &
                 mgrid,nestgrd, hdmpfn, ldmpf)


    IF (myproc == 0)                                                    &
      WRITE(6,'(1x,a,a)') 'History data dump in file ',hdmpfn(1:ldmpf)

    grdbas=0
    ! blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride

      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

        CALL dtadump(nx,ny,nz,nzsoil,nstyp,                             &
               hdmpfmt,iniotfu,hdmpfn(1:ldmpf),grdbas,filcmprs,         &
               u,v,w,ptprt,pprt,qv,qscalar,                             &
               tem1,tem1,tem1,                                          &
               ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                  &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               rain,tem1,tem1,                                          &
               tem1,tem1,tem1,tem1,tem1,                                &
               tem1,tem1,tem1,tem1,                                     &
               tem2,tem3,tem4)

      END IF

      IF (mp_opt > 0) CALL mpbarrier

    END DO

    first_time = 0

!-----------------------------------------------------------------------
!
!  Calculate and write out certian 2D fields for post-processing
!  (both in binary or other formats)
!
!-----------------------------------------------------------------------
!
    IF(i2dfmt > 0)  THEN

      u = u - ubar
      v = v - vbar
      qv = qv - qvbar

      tem9 = -999.0

      ilite=0
      iltgci=0
      ibeg_offset=0
      iend_offset=0
      jbeg_offset=0
      jend_offset=0

      lenbin = len_trim(outheader)
      IF(hdmpfn(ldmpf-2:ldmpf-2) == '.') ldmpf = ldmpf - 3

      CALL postcore(nx,ny,nz,nzsoil,nstyp,curtim, 2,  & ! hard-coded with mp_physics=2 temporarily, add fine-control when needed
              hdmpfn(1:ldmpf),ldmpf,outheader(1:lenbin),                &
              lenbin,' ',1,                                             &
              9,1,0,i2dfmt,0,                                           &
              ilite,iltgci,0,2,3,3,0,1,                 & ! hard-coded icrtm,isatid,chbgn,chend=0,2,3,3 for no CRTM call needed with ext2arps
              ibeg_offset,iend_offset,jbeg_offset,jend_offset,          &
              x,y,z,zp,u ,v ,w ,ptprt, pprt ,                           &
              qv, qscalar,                                              &
              ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,             &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,                   &
              tsoil,qsoil,wetcanp,snowdpth,                             &
              rain,tem9,tem9,tem9,tem9,                                 &
              t_2m,tem9,qv_2m,u_10m,v_10m,tem9,                         &
              tem9,tem9,tem9,tem9,tem9,tem9,                            &
              tem9,tem9,tem9,                                           &
              tem9,tem9,tem9,tem9,tem9,                                 &
              tem9,tem9,tem9,tem9,                                      &
              tem1,tem2,tem3,tem4,tem5,tem6)

    END IF

  END DO
!
!-----------------------------------------------------------------------
!
!  Friendly exit message.
!
!-----------------------------------------------------------------------
!
  IF (istatus /= 1) GOTO 999

  IF (myproc == 0) WRITE(6,'(a/a,i4,a)')                                &
          ' ==== Normal successful completion of EXT2ARPS. ====',       &
          '      Processed',nextdfil,' file(s)'
  IF (mp_opt > 0) CALL mpexit(0)
  STOP
!
!-----------------------------------------------------------------------
!
!  Problem doing time conversions.
!
!-----------------------------------------------------------------------
!
  920 CONTINUE
  IF (myproc == 0) WRITE(6,'(/,a,/,a,i4,a,i4/,a,a19)')                  &
          ' Aborting, error in time format for external file',          &
          ' File number:',ifile,' of',nextdfil,                         &
          ' External file time provided:',extdtime(ifile)
  CALL arpsstop(' ',1)
!
!-----------------------------------------------------------------------
!
!  Error status returned from rdextfil
!
!-----------------------------------------------------------------------
!
  999 CONTINUE
  IF (myproc == 0) WRITE(6,'(/,a,i6)')                                  &
          ' Aborting, error reading external file. istatus=',           &
            istatus
  CALL arpsstop(' ',1)

  STOP
END PROGRAM ext2arps
