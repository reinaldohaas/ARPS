PROGRAM arps2eta212
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARSP2ETA212                 ######
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
!  Converts ARPS history files to a ETA GRIB #212 format file.
!
!  Reads in a history file produced by ARPS in any ARPS format.
!  Converts variables and interpolates to a fixed set of pressure
!  levels.
!
!  NOTE:
!    Please make sure the ARPS domain covers ETA Grid #212 domain
!    before running this program.
!
!  INPUT Namelist
!    arps2eta212.input
!
!  Libraries:
!    libarps.a
!    libadas.a
!    HDF 4 library (When ARPS is in HDF 4 format)
!
!  Subroutine CALLs:
!    dtaread
!    mkarps2d
!    v2dint
!
!  Subroutine defined in this file
!    extrph
!    extrpt
!    extrpq
!    extrpuv
!    cal_avor
!    gtsinlat_ext
!    mkipds_212
!    mkigds_212
!    wrtgb        --  pack and write GRIB message for Grid #212
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  04/19/2003
!
!  MODIFICATION HISTORY:
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

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,lenfil
  PARAMETER (nhisfile_max=200)
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
  REAL, ALLOCATABLE :: qi     (:,:,:) ! Cloud ice mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:) ! Microphysics scalar array
  REAL, ALLOCATABLE :: qv     (:,:,:) ! Water vapor specific humidity (kg/kg)

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

  REAL, ALLOCATABLE :: raing3hr(:,:)     ! 3 hour precipitation accumulation
  REAL, ALLOCATABLE :: rainc3hr(:,:)     ! 3 hour precipitation accumulation
  REAL, ALLOCATABLE :: rain3hr(:,:)      ! 3 hour precipitation accumulation
  REAL, ALLOCATABLE :: raingaccum(:,:)   ! precipitation accumulation 3 hr before
  REAL, ALLOCATABLE :: raincaccum(:,:)   ! precipitation accumulation 3 hr before
!
!-----------------------------------------------------------------------
!
!  Misc ARPS variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nchr, nchw, nchoutr
  REAL :: latnot(2), swx,swy,ctrx,ctry
  LOGICAL :: init
  REAL :: rlnp, soildp

  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
  REAL, ALLOCATABLE :: zps_km(:,:,:)

  REAL :: gamma,ex1,ex2,rln700,p00,t00
  REAL, ALLOCATABLE :: p_mb(:,:,:)

  REAL, ALLOCATABLE :: mf2d(:,:)

  REAL, ALLOCATABLE :: avor(:,:,:)

  REAL, ALLOCATABLE :: uear(:,:,:)
  REAL, ALLOCATABLE :: vear(:,:,:)
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

  REAL, ALLOCATABLE :: wrk1(:),wrk2(:),wrk3(:),wrk4(:),wrk5(:),wrk6(:)
  REAL, ALLOCATABLE :: wrk7(:),wrk8(:),wrk9(:),wrk10(:),wrk11(:),wrk12(:)
  REAL, ALLOCATABLE :: tem2d1(:,:)
  REAL, ALLOCATABLE :: tem2d2(:,:)
  REAL, ALLOCATABLE :: tem2d3(:,:)
  REAL, ALLOCATABLE :: tem2d4(:,:)

!
!-----------------------------------------------------------------------
!
!  ETA output variables
!
!-----------------------------------------------------------------------
!
!  Constants

  INTEGER, PARAMETER :: nx_eta = 185, ny_eta = 129, nz_eta = 39
  REAL,    PARAMETER :: dx_eta = 40635.25, dy_eta = 40635.25  ! Meters

  INTEGER :: nprlvl
  INTEGER :: iprlvl(nz_eta)
  INTEGER, PARAMETER :: iprbgn = 1000, iprend = 50, iprinc = -25  ! mb

  INTEGER, PARAMETER :: iproj_eta = 3
  INTEGER, PARAMETER :: iproj_eia = 2
  REAL, PARAMETER    :: scale_eta = 1.0
  REAL, PARAMETER    :: latsw   =   12.190, lonsw   = -133.459,         &
                        latne   =   57.290, lonne   =  -49.385,         &
                        lattru_eta(2) = (/25.00, 25.00/),               &
                        lontru_eta =   -95.00
  REAL :: x0_eta, y0_eta

  REAL, ALLOCATABLE :: x_eta(:), y_eta(:)
  REAL, ALLOCATABLE :: lat_eta(:,:), lon_eta(:,:)
  REAL, ALLOCATABLE :: umap(:,:)
  REAL, ALLOCATABLE :: vmap(:,:)

!-----------------------------------------------------------------------
!
!  Working arrays
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:,:)
  REAL, ALLOCATABLE :: sinlat(:,:)

  REAL, ALLOCATABLE :: tem2d_eta(:,:)    ! store data to be packed and written
  REAL, ALLOCATABLE :: tem2d2_eta(:,:)   ! store data to be packed and written

  REAL, ALLOCATABLE    :: wrkhori(:,:,:) ! store horizontal interpolations
  REAL, ALLOCATABLE    :: wrkhori2(:,:,:)! store horizontal interpolations
  REAL, ALLOCATABLE    :: p_eta(:,:,:)   ! store horizontally interpolated pressure
  REAL, ALLOCATABLE    :: lnp_eta(:,:,:) ! store -ln(p) at ETA grid

  INTEGER, ALLOCATABLE :: iloc(:,:)      ! x-index/y-index location of
  INTEGER, ALLOCATABLE :: jloc(:,:)      ! ETA grid point in the ARPS array

  REAL, ALLOCATABLE    :: x2d(:,:)       ! x/y coordinate of ETA grid point
  REAL, ALLOCATABLE    :: y2d(:,:)       ! in ARPS coordinate

  REAL, ALLOCATABLE :: dxfld(:), dyfld(:)
  REAL, ALLOCATABLE :: rdxfld(:), rdyfld(:)
  REAL, ALLOCATABLE :: slopey(:,:),alphay(:,:),betay(:,:)

!-----------------------------------------------------------------------
!
! Grib message variables
!
!-----------------------------------------------------------------------
  INTEGER :: wdlen             ! Length of machine word
  INTEGER :: ipds(25)          ! PDS integer array
  INTEGER :: igds(25)          ! GDS integer array

!-----------------------------------------------------------------------
!
!  Variables to determine the length of machine word, wdlen
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: ctema,ctemb
  INTEGER :: itema,itemb
  EQUIVALENCE ( itema,ctema )
  EQUIVALENCE ( itemb,ctemb )
  DATA ctema / '12345678' /

  REAL,    ALLOCATABLE :: mtem(:,:)
  INTEGER, ALLOCATABLE :: mitem1(:,:), mitem2(:,:)


!-----------------------------------------------------------------------
!
!  Namelists
!
!-----------------------------------------------------------------------
!
  INTEGER :: iorder  ! order of polynomial for interpolation (1, 2 or 3)
  NAMELIST /interpolation_options/ iorder

  INTEGER :: vvelout
  NAMELIST /output/ dirname,grbpkbit,readyfl, vvelout
                     ! these variables are declared in globcst.inc
!
!-----------------------------------------------------------------------
!
!  External functions
!
!-----------------------------------------------------------------------
!
  REAL :: wmr2td,oe,dpt, compute_density
  EXTERNAL wmr2td,oe,dpt, compute_density
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,nq,klev,ifile,ireturn
  INTEGER :: istatus, fcttm
  REAL :: time
  INTEGER :: p1time,p2time
  CHARACTER(LEN=256) :: filname
  CHARACTER(LEN=256) :: filnamr
  INTEGER :: lenstr

  INTEGER :: varid
  REAL    :: rho

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
  WRITE(6,'(9(/5x,a))')                                                 &
      '###############################################################',&
      '###############################################################',&
      '##                                                           ##',&
      '##                Welcome to ARPS2ETA212                     ##',&
      '## a program that reads in history files generated by ARPS   ##',&
      '## and produces a ETA GRIB #212 format file.                 ##',&
      '##                                                           ##',&
      '###############################################################',&
      '###############################################################'
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lengbf = len_trim(grdbasfn)

  CALL get_dims_from_data(hinfmt,TRIM(hisfile(1)),                   &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(3(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz

  nstyp = nstyps
!
!-----------------------------------------------------------------------
!
!  Get output options and the name of the grid/base data set.
!
!-----------------------------------------------------------------------
!
   grbpkbit = 16
   dirname  = './'
   readyfl  = 0
   vvelout  = 2      ! dump omega by default

   READ(5,interpolation_options)
   READ(5,output)

!-----------------------------------------------------------------------
!
! Allocate ARPS arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(x      (nx))
  ALLOCATE(y      (ny))
  ALLOCATE(z      (nz))
  ALLOCATE(zp     (nx,ny,nz))
  ALLOCATE(zpsoil (nx,ny,nzsoil))
  ALLOCATE(zps_km(nx,ny,nz))
  ALLOCATE(p_mb(nx,ny,nz))

  ALLOCATE(uprt   (nx,ny,nz))
  ALLOCATE(vprt   (nx,ny,nz))
  ALLOCATE(wprt   (nx,ny,nz))
  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))
  ALLOCATE(qv     (nx,ny,nz))
  ALLOCATE(qi     (nx,ny,nz))
  ALLOCATE(qscalar (nx,ny,nz,nscalar))
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

  ALLOCATE(raing3hr(nx,ny))
  ALLOCATE(rainc3hr(nx,ny))
  ALLOCATE(rain3hr(nx,ny))
  ALLOCATE(raingaccum(nx,ny))
  ALLOCATE(raincaccum(nx,ny))

  ALLOCATE(radfrc(nx,ny,nz))
  ALLOCATE(radsw (nx,ny))
  ALLOCATE(rnflx (nx,ny))
  ALLOCATE(radswnet(nx,ny))
  ALLOCATE(radlwin(nx,ny))

  ALLOCATE(usflx (nx,ny))
  ALLOCATE(vsflx (nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))

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
  qi     =0.0
  qscalar =0.0
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

  raing3hr = 0.0
  rainc3hr = 0.0
  rain3hr  = 0.0
  raingaccum = 0.0
  raingaccum = 0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet =0.0
  radlwin =0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  ALLOCATE(xs(nx), STAT = istatus)
  ALLOCATE(ys(ny), STAT = istatus)
  ALLOCATE(zps(nx,ny,nz), STAT = istatus)
  zps    = 0.0
  zps_km = 0.0
  p_mb   = 0.0

  ALLOCATE ( uear  (nx,ny,nz))
  ALLOCATE ( vear  (nx,ny,nz))
  uear = 0.0
  vear = 0.0

  ALLOCATE ( e_mb   (nx,ny,nz))    ! vapor pressure in mb

  ALLOCATE ( mix    (nx,ny,nz))    ! total mixing ratio (kg/kg)
  ALLOCATE ( esat_mb(nx,ny,nz))    ! saturation vapor pressure in mb
  ALLOCATE ( rh     (nx,ny,nz))    ! Relative humidity in %
  ALLOCATE ( t_dew  (nx,ny,nz))    ! dewpoint temp. in degrees K

  ALLOCATE ( avor  (nx,ny,nz))     ! absolute vorticity
  avor = 0.0

  ALLOCATE(lcl(nx,ny))
  ALLOCATE(lfc(nx,ny))
  ALLOCATE(el(nx,ny))
  ALLOCATE(twdf(nx,ny))
  ALLOCATE(li(nx,ny))
  ALLOCATE(pbe(nx,ny))
  ALLOCATE(mbe(nx,ny))
  ALLOCATE(nbe(nx,ny))
  ALLOCATE(tcap(nx,ny))
  lcl =0.0
  lfc =0.0
  el  =0.0
  twdf=0.0
  li  =0.0
  pbe =0.0
  mbe =0.0
  nbe =0.0
  tcap=0.0

  ALLOCATE(wrk1(nz),wrk2(nz),wrk3(nz),wrk4(nz),wrk5(nz),wrk6(nz))
  ALLOCATE(wrk7(nz),wrk8(nz),wrk9(nz),wrk10(nz),wrk11(nz),wrk12(nz))

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

  ALLOCATE(tem2d1(nx,ny))
  ALLOCATE(tem2d2(nx,ny))
  ALLOCATE(tem2d3(nx,ny))
  ALLOCATE(tem2d4(nx,ny))
  tem2d1=0.0
  tem2d2=0.0
  tem2d3=0.0
  tem2d4=0.0

  ALLOCATE( mf2d(nx,ny))
  mf2d=0.0

  ALLOCATE(heli (nx,ny))
  ALLOCATE(brn  (nx,ny))
  ALLOCATE(brnu (nx,ny))
  ALLOCATE(srlfl(nx,ny))
  ALLOCATE(srmfl(nx,ny))
  ALLOCATE(shr37(nx,ny))
  ALLOCATE(ustrm(nx,ny))
  ALLOCATE(vstrm(nx,ny))
  ALLOCATE(blcon(nx,ny))

  heli =0.0
  brn  =0.0
  brnu =0.0
  srlfl=0.0
  srmfl=0.0
  shr37=0.0
  ustrm=0.0
  vstrm=0.0
  blcon=0.0
!---------------------------------------------------------------------------
!
! Allocate eta arrays
!
!---------------------------------------------------------------------------

  ALLOCATE(x_eta  (nx_eta), STAT = istatus)
  ALLOCATE(y_eta  (ny_eta), STAT = istatus)
  ALLOCATE(lat_eta(nx_eta,ny_eta), STAT = istatus)
  ALLOCATE(lon_eta(nx_eta,ny_eta), STAT = istatus)

  ALLOCATE ( umap  (nx_eta,ny_eta))
  ALLOCATE ( vmap  (nx_eta,ny_eta))
  umap = 0.0
  vmap = 0.0

!--------------------------------------------------------------------------
!
! Allocate temporary arrays
!
!--------------------------------------------------------------------------

  ALLOCATE(lnp_eta(nx_eta,ny_eta,nz), STAT = istatus)
  lnp_eta = 0.0

  ALLOCATE(p_eta(nx_eta,ny_eta,nz), STAT = istatus)
  p_eta = 0.0

  ALLOCATE(wrkhori(nx_eta,ny_eta,nz), STAT = istatus)
  wrkhori = 0.0
  ALLOCATE(wrkhori2(nx_eta,ny_eta,nz), STAT = istatus)
  wrkhori2 = 0.0

  ALLOCATE(iloc(nx_eta,ny_eta), STAT = istatus)
  ALLOCATE(jloc(nx_eta,ny_eta), STAT = istatus)
  iloc = 0
  jloc = 0
  ALLOCATE(x2d (nx_eta,ny_eta), STAT = istatus)
  ALLOCATE(y2d (nx_eta,ny_eta), STAT = istatus)
  x2d = 0.0
  y2d = 0.0

  ALLOCATE( dxfld(nx), STAT = istatus)
  ALLOCATE( dyfld(ny), STAT = istatus)
  ALLOCATE(rdxfld(nx), STAT = istatus)
  ALLOCATE(rdyfld(ny), STAT = istatus)
  dxfld = 0.0
  dyfld = 0.0
  rdxfld = 0.0
  rdyfld = 0.0

  ALLOCATE(slopey(nx,ny), STAT = istatus)
  ALLOCATE(alphay(nx,ny), STAT = istatus)
  ALLOCATE( betay(nx,ny), STAT = istatus)
  slopey = 0.0
  alphay = 0.0
  betay  = 0.0

  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))
  ALLOCATE(tem4(nx,ny,nz))
  ALLOCATE(sinlat(nx,ny))
  sinlat = 0.0

  tem1 =0.0
  tem2 =0.0
  tem3 =0.0
  tem4 =0.0

  ALLOCATE(tem2d_eta(nx_eta,ny_eta), STAT = istatus)
  tem2d_eta = 0.0
  ALLOCATE(tem2d2_eta(nx_eta,ny_eta), STAT = istatus)
  tem2d2_eta = 0.0

!-----------------------------------------------------------------------
!
!  Calculate the length of machine word, wdlen, in bytes, which will
!  be used to pack the data. Assume the length has only two possible
!  value: 4 and 8
!
!-----------------------------------------------------------------------
!
    itemb = itema
    IF ( ctema /= ctemb ) THEN
      wdlen = 4
    ELSE
      wdlen = 8
    END IF
    WRITE (6,'(a,i2,a)')  &
        'The length of machine word is ',wdlen,' bytes'

!-----------------------------------------------------------------------
!
! Initilize GRIB message variables
!
!-----------------------------------------------------------------------

  CALL mkigds_212(nx_eta,ny_eta,nprlvl,igds)

!
  ALLOCATE(mtem(nx_eta,ny_eta),   STAT = istatus)
  ALLOCATE(mitem1(nx_eta,ny_eta), STAT = istatus)
  ALLOCATE(mitem2(nx_eta,ny_eta), STAT = istatus)
  mtem = 0.0
  mitem1 = 0
  mitem2 = 0

!-----------------------------------------------------------------------
!
!  Find x,y locations of ETA grid.
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj_eia,scale_eta,lattru_eta,lontru_eta)
  CALL setorig(2,latsw,lonsw)
  CALL lltoxy(1,1,latsw,lonsw,x0_eta,y0_eta)
!  CALL setorig(1,x0_eta,y0_eta)

  DO i=1,nx_eta
    x_eta(i)=x0_eta+(i-1)*dx_eta
  END DO
  DO j=1,ny_eta
    y_eta(j)=y0_eta+(j-1)*dy_eta
  END DO
  CALL xytoll(nx_eta,ny_eta,x_eta,y_eta,lat_eta,lon_eta)

  DO klev = iprbgn,iprend,iprinc
    iprlvl((iprbgn-klev)/ABS(iprinc) + 1) = klev
  END DO
  nprlvl = SIZE(iprlvl)
  IF(nprlvl > nz_eta) THEN
    WRITE(6,*) 'ERROR: nz_eta is too small.'
    STOP
  END IF

!--------------------------------------------------------------------
!
! Begin loop over all of the history files
!
!--------------------------------------------------------------------

  lengbf=LEN_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

  init = .false.
  nchw = 10

  DO ifile = 1,nhisfile
    lenfil = LEN(hisfile(ifile))
    CALL strlnth( hisfile(ifile), lenfil)
    WRITE(6,'(/a,/1x,a)')' The data set name is ',                      &
                                           hisfile(ifile)(1:lenfil)

!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    CALL dtaread(nx,ny,nz,nzsoil,nstyps,hinfmt,nchr,grdbasfn(1:lengbf), &
               lengbf,                                                  &
               hisfile(ifile)(1:lenfil),lenfil,time,                    &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt, qscalar, tke, kmh, kmv,                &
               ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,            &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)

!    IF(nx < nx_eta .OR. ny < ny_eta) THEN
!      WRITE(6,'(a)') 'ARPS domain is smaller than ETA 212 domain.'
!      STOP
!    END IF
!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
    IF( ireturn == 0 ) THEN   ! successful read

    !WRITE(6,*) 'time =', time

    fcttm = INT(time/3600)
    CALL gtlfnkey(runname, lenstr)

    WRITE(filname,'(a,a1,I4.4,I2.2,I2.2,I2.2,a,I2.2,a4)')              &
           runname(1:lenstr),'.',year,month,day,hour,'f',fcttm,'.grb'

    p1time = INT((time-10800)/3600)      ! in minutes
    p2time = INT(time/3600)              ! in minutes

    lenstr = LEN_TRIM(dirname)
    IF(lenstr > 0) THEN
      IF(dirname(lenstr:lenstr) /= '/') THEN
        dirname(lenstr+1:lenstr+1) = '/'
        lenstr = lenstr + 1
      END IF
      filname = dirname(1:lenstr) // TRIM(filname)
      print *, 'Output file name is ',filname
    END IF

    lenstr = LEN_TRIM(filname)

    CALL asnctl ('NEWLOCAL', 1, istatus)
    CALL asnfile(filname, '-F f77 -N ieee', istatus)

    CALL getunit( nchw )

    OPEN(UNIT=nchw,FILE=trim(filname),STATUS='unknown',             &
                   FORM='unformatted',IOSTAT= istatus )

    IF(istatus /= 0) THEN
      WRITE(6,*) 'Error while creating file.',TRIM(filname)
      STOP
    END IF

    IF( .NOT.init ) THEN

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
!  Find x,y locations of ETA grid in terms of the ARPS grid.
!
!-----------------------------------------------------------------------
        latnot(1) = trulat1
        latnot(2) = trulat2
        CALL setmapr(mapproj,1.0/sclfct,latnot,trulon)
        CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
        swx = ctrx - (FLOAT((nx-3))/2.) * dx/sclfct
        swy = ctry - (FLOAT((ny-3))/2.) * dy/sclfct
        CALL setorig(1,swx,swy)

        CALL xytomf(nx,ny,xs,ys,mf2d)

        CALL lltoxy(nx_eta,ny_eta,lat_eta,lon_eta,x2d,y2d)

        CALL setijloc(nx_eta,ny_eta,nx,ny,x2d,y2d,                    &
                    xs,ys,iloc,jloc)


        CALL setdxdy(nx,ny,1,nx,1,ny, xs,ys,dxfld,dyfld,rdxfld,rdyfld)

        init=.true.
     END IF  ! NOT init

! Restore to use ETA map projection

     CALL setmapr(iproj_eia,scale_eta,lattru_eta,lontru_eta)
     CALL setorig(2,latsw,lonsw)
!
!-----------------------------------------------------------------------
!
!    Begin ARPS data conversions
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k)=pprt(i,j,k)+pbar(i,j,k)
            ptprt(i,j,k)=ptprt(i,j,k)+ptbar(i,j,k)
            qvprt(i,j,k)=qvprt(i,j,k)+qvbar(i,j,k)
            tem1(i,j,k)=0.5*(uprt(i,j,k)+ubar(i,j,k)+                 &
                        uprt(i+1,j,k)+ubar(i+1,j,k))
            tem2(i,j,k)=0.5*(vprt(i,j,k)+vbar(i,j,k)+                 &
                        vprt(i,j+1,k)+vbar(i,j+1,k))
            tem3(i,j,k)=0.5*(wprt(i,j,k)+wbar(i,j,k)+                 &
                        wprt(i,j,k+1)+wbar(i,j,k+1))

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
!  Put temperature into tem2
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=ptprt(i,j,k)*(pprt(i,j,k)/p0)**rddcp
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!
!  Processing Pressure data
!
!-----------------------------------------------------------------------

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,pprt(:,:,k),p_eta(:,:,k),              &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
        lnp_eta(:,:,k) = -ALOG(p_eta(:,:,k))
      END DO

!-----------------------------------------------------------------------
!
!  Processing Temperature data
!
!-----------------------------------------------------------------------

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,tem2(:,:,k),wrkhori(:,:,k),            &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
        ! here use mkarps2d in the reverse direction
        ! interpolate from ARPS grid to ETA grid instead of the other
      END DO

      PRINT *, ' Storing Temperature GRIB data.'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels
        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrpt(nx_eta,ny_eta,nz,wrkhori,p_eta,                     &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message
        CALL mkipds_212(11,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)
      END DO

!  Find temperature in 2m
      DO k = 1, nz
        DO j = 1,ny_eta
          DO i = 1, nx_eta
            wrkhori2(i,j,k) = z(k)
          END DO
        END DO
      END DO

      PRINT *, ' Storing Temperature at 2m.'

      CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                   wrkhori,wrkhori2,2,tem2d_eta)
      CALL mkipds_212(11,105,2,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                 ipds,igds, mtem,mitem1,mitem2)

!-----------------------------------------------------------------------
!
!  Process Specific Humidity (51)
!
!-----------------------------------------------------------------------

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,          &
                   xs,ys,x2d,y2d,qvprt(:,:,k),wrkhori(:,:,k),         &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,  &
                   tem4)
      END DO

      PRINT *, ' Storing Specific Humidity GRIB data.'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels

        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,      &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrpq(nx_eta,ny_eta,nz,wrkhori,p_eta,                  &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message

        CALL mkipds_212(51,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,      &
                   ipds,igds, mtem,mitem1,mitem2)
      END DO

!-----------------------------------------------------------------------
!
!  Process U wind (33) and V wind (34)
!
!  Wind needs to be rotated to ETA grid
!
!-----------------------------------------------------------------------

! First rotate the wind to true-earth-relative in ARPS map projection
      CALL setmapr(mapproj,1.0/sclfct,latnot,trulon)
      CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
      swx = ctrx - (FLOAT((nx-3))/2.) * dx/sclfct
      swy = ctry - (FLOAT((ny-3))/2.) * dy/sclfct
      CALL setorig(1,swx,swy)
      CALL xytoll(nx,ny,xs,ys,tem2d1,tem2d2)
      DO k = 1,nz
        CALL uvmptoe(nx,ny,uprt(:,:,k),vprt(:,:,k),                 &
                     tem2d2,uear(:,:,k),vear(:,:,k))
      END DO

! Second rotate the true earth wind to ETA map projection

      CALL setmapr(iproj_eia,scale_eta,lattru_eta,lontru_eta)
      CALL setorig(2,latsw,lonsw)

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1

        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,        &
                   xs,ys,x2d,y2d,uear(:,:,k),wrkhori(:,:,k),        &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,&
                   tem4)

        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,        &
                   xs,ys,x2d,y2d,vear(:,:,k),wrkhori2(:,:,k),       &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,&
                   tem4)
      END DO

!-----------------------------------------------------------------------
!
!  10m u and v
!  Dan recommended to use k = 2 is 10m
!
!-----------------------------------------------------------------------

      CALL uvetomp(nx_eta,ny_eta,wrkhori(:,:,2),wrkhori2(:,:,2),    &
                   lon_eta,umap,vmap)

      PRINT *, ' Storing horizontal wind at 10m.'

      CALL mkipds_212(33,105,10,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,umap,nchw,wdlen,grbpkbit,            &
                   ipds,igds, mtem,mitem1,mitem2)

      CALL mkipds_212(34,105,10,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,vmap,nchw,wdlen,grbpkbit,            &
                   ipds,igds, mtem,mitem1,mitem2)

!------------------------------------------------------------------------

      PRINT *, ' Storing Horizontal wind GRIB data.'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels

        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)
        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori2,lnp_eta,rlnp,tem2d2_eta)

        CALL extrpuv(nx_eta,ny_eta,nz,wrkhori,wrkhori2,p_eta,           &
                     iprlvl(klev),tem2d_eta,tem2d2_eta)

        ! 3. Pack and write GRIB message

        CALL uvetomp(nx_eta,ny_eta,tem2d_eta,tem2d2_eta,lon_eta,umap,vmap)

        CALL mkipds_212(33,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,umap,nchw,wdlen,grbpkbit,              &
                   ipds,igds, mtem,mitem1,mitem2)

        CALL mkipds_212(34,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,vmap,nchw,wdlen,grbpkbit,              &
                   ipds,igds, mtem,mitem1,mitem2)
      END DO

!-----------------------------------------------------------------------
!
!  Process Geopotential Heigth (gpm) (7)
!
!-----------------------------------------------------------------------

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,zps(:,:,k),wrkhori(:,:,k),             &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
      END DO

      PRINT *, ' Storing Geopotential Height GRIB data.'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels

        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrph(nx_eta,ny_eta,nz,wrkhori,wrkhori,p_eta,             &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message

        CALL mkipds_212(07,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,        &
                   ipds,igds, mtem,mitem1,mitem2)
      END DO

!-----------------------------------------------------------------------
!
!  Process W wind (40)   (ETA uses Pa S**-1 -- 39)
!
!-----------------------------------------------------------------------

      IF(vvelout == 2) THEN   ! convert W to Omega
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              rho = compute_density(tem2(i,j,k),pprt(i,j,k))
              wprt(i,j,k)= -rho*g*wprt(i,j,k)
            END DO
          END DO
        END DO
        varid = 39
      ELSE IF(vvelout == 1) THEN
        varid = 40
      ELSE
        varid = 0
      END IF

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      IF(varid /= 0) THEN

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,wprt(:,:,k),wrkhori(:,:,k),            &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
      END DO

      PRINT *, ' Storing W-wind GRIB data.'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels

        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrpq(nx_eta,ny_eta,nz,wrkhori,p_eta,                     &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message

        CALL mkipds_212(varid,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)
      END DO

      END IF   ! varid /= 0


!-----------------------------------------------------------------------
!
!  Process soil temperature (85)
!    level type = 111   Depth below land surface
!    ETA use 112 layer between two depths below land surface
!
!-----------------------------------------------------------------------

      PRINT *, ' Storing soil temperature at layer below land surface.'

      DO k = 1,nzsoil
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,tsoil(:,:,k,0),tem2d_eta,              &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
        soildp = 100*(zpsoil(2,2,1)-zpsoil(2,2,k))

        CALL mkipds_212(85,111,soildp,0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

      END DO

!-----------------------------------------------------------------------
!
!  Process soil moisture content (144)
!    level type = 111   Depth below land surface
!    ETA use 112 layer between two depths below land surface
!
!-----------------------------------------------------------------------

      PRINT *, ' Storing soil moisture GRIB data.'

      DO k = 1,nzsoil
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,qsoil(:,:,k,0),tem2d_eta,              &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
        soildp = 100*(zpsoil(2,2,1)-zpsoil(2,2,k))
                              !   layer below land surface, cm
        CALL mkipds_212(144,111,soildp,0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

      END DO

!-----------------------------------------------------------------------
!
!  Water equiv. of accum snow depth (kg/m**2) -- 65
!
!----------------------------------------------------------------------

      PRINT *, ' Storing accum. snow depth.'

      tem1(:,:,1) = 100* snowdpth(:,:)

      CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                  xs,ys,x2d,y2d,tem1(:,:,1),tem2d_eta,                &
                  dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                  tem4)

      CALL mkipds_212(65,01,00,0,p2time,0,10,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

!----------------------------------------------------------------------
!
!  Convective precip (63)
!
!  rainc unit is mm
!  63 require kg m**-2  = rhowater (1000 kg m**-3) * rainc (mm) * 1.0E-3
!
!----------------------------------------------------------------------
      IF (time >= 10800 .AND. MOD(INT(time),10800) == 0) THEN
                       ! precipitation dump every 3 hours
        raing3hr = raing - raingaccum
        rainc3hr = rainc - raincaccum
        raingaccum = raing
        raingaccum = rainc
        rain3hr = raing3hr + rainc3hr

        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                    xs,ys,x2d,y2d,rain3hr,tem2d_eta,                    &
                    dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                    tem4)
        print *, ' Storing total precipitation rainfall.'
        CALL mkipds_212(61,01,00,4,p1time,p2time,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                    xs,ys,x2d,y2d,raing3hr,tem2d_eta,                   &
                    dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                    tem4)
        print *, ' Storing large scale precipitation.'
        CALL mkipds_212(62,01,00,4,p1time,p2time,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                    xs,ys,x2d,y2d,rainc3hr,tem2d_eta,                   &
                    dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                    tem4)
        print *, ' Storing convective rainfall.'
        CALL mkipds_212(63,01,00,4,p1time,p2time,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

      END IF
!-----------------------------------------------------------------------
!
!  Output near-sfc data.
!  Simularity theory or something could be applied here
!  to make these data be valid at sfc instrument height,
!  for now we output level 2.
!
!-----------------------------------------------------------------------

      CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                  xs,ys,x2d,y2d,pprt(:,:,2),tem2d_eta,                &
                  dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                  tem4)
      print *, ' Storing near-sfc presure.'

      CALL mkipds_212(01,01,0,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                 ipds,igds, mtem,mitem1,mitem2)

!-----------------------------------------------------------------------
!
!  Geopotential Height (gpm)
!
!-----------------------------------------------------------------------

      CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                  xs,ys,x2d,y2d,zps(:,:,2),tem2d_eta,                 &
                  dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                  tem4)

      print *, ' Storing near-sfc Geopotential Height'

      CALL mkipds_212(07,01,00,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,       &
                 ipds,igds, mtem,mitem1,mitem2)

!-----------------------------------------------------------------------
!
! Surface Temperature (11)
!
!-----------------------------------------------------------------------

      CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                  xs,ys,x2d,y2d,tem2(:,:,2),tem2d_eta,                &
                  dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,   &
                  tem4)

      print *, ' Storing near-sfc temperature'

      CALL mkipds_212(11,01,00,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                 ipds,igds, mtem,mitem1,mitem2)

!----------------------------------------------------------------------
!
!  Calculation t_dew, rh etc.
!
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
!
!  Process Relative Humidity (52)
!
!-----------------------------------------------------------------------

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,rh(:,:,k),wrkhori(:,:,k),              &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
      END DO

      PRINT *, ' Storing Relative Humidity GRIB data.'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels

        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrpq(nx_eta,ny_eta,nz,wrkhori,p_eta,                     &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message

        DO j = 1, ny_eta
           DO i = 1, nx_eta
             tem2d_eta(i,j) = MAX( MIN(tem2d_eta(i,j), 100.0), 0.0)
           END DO
        END DO

        CALL mkipds_212(52,100,iprlvl(klev),0,p2time,0,0,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

      END DO

!-----------------------------------------------------------------------
!
!  Process Dew-point temperature (17) at 850mb only
!
!-----------------------------------------------------------------------

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,t_dew(:,:,k),wrkhori(:,:,k),           &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
      END DO

      rlnp = -ALOG(100.*FLOAT(850))

      PRINT *, ' Storing Dew-point temp. at pressure = ',850, 'mb.'

      ! 2. Vertically interpolate to pressure levels

      CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

      CALL extrpt(nx_eta,ny_eta,nz,wrkhori,p_eta,                     &
                  850,tem2d_eta)

      ! 3. Pack and write GRIB message

      CALL mkipds_212(17,100,850,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                 ipds,igds, mtem,mitem1,mitem2)

!----------------------------------------------------------------------
!
!  Find dew-point temperature at 2m
!
!---------------------------------------------------------------------

      DO k = 1, nz
        DO j = 1,ny_eta
          DO i = 1, nx_eta
            wrkhori2(i,j,k) = z(k)
          END DO
        END DO
      END DO

      CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,           &
                   wrkhori,wrkhori2,2,tem2d_eta)

      PRINT *, ' Storing dew-point temperature at 2m'

      CALL mkipds_212(17,105,2,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                 ipds,igds, mtem,mitem1,mitem2)

!-----------------------------------------------------------------------
!
!  Process Cloud ice, (58) kg m**-2
!
!-----------------------------------------------------------------------

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            IF(P_QI > 0) THEN
              tem1(i,j,k) = rhobar(i,j,k) * qi(i,j,k)       ! ?
            ELSE
              tem1(i,j,k) = 0.0
            END IF
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,tem1(:,:,k),wrkhori(:,:,k),            &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
      END DO

      PRINT *, ' Storing Could ice GRIB data'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels

        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrpq(nx_eta,ny_eta,nz,wrkhori,p_eta,                     &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message

        CALL mkipds_212(58,100,iprlvl(klev),0,p2time,0,10,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)
      END DO

!-----------------------------------------------------------------------
!
!  Calculate abosolute vorticity
!
!-----------------------------------------------------------------------

      CALL cal_avor(avor,uprt,vprt,x,y,nx,ny,nz,1,0,omega,              &
                     sinlat,xs,ys,tem2d1)

      avor = avor * 1.0E-5

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      DO k = 1,nz-1
        CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,avor(:,:,k),wrkhori(:,:,k),            &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)
      END DO

      PRINT *, ' Storing absolut vor. GRIB data'

      DO klev=1,nprlvl
        rlnp = -ALOG(100.*FLOAT(iprlvl(klev)))

        ! 2. Vertically interpolate to pressure levels
        CALL v2dinta(nx_eta,ny_eta,nz,1,nx_eta,1,ny_eta,1,nz-1,         &
                      wrkhori,lnp_eta,rlnp,tem2d_eta)

        CALL extrpq(nx_eta,ny_eta,nz,wrkhori,p_eta,                     &
                    iprlvl(klev),tem2d_eta)

        ! 3. Pack and write GRIB message


        CALL mkipds_212(41,100,iprlvl(klev),0,p2time,0,10,ipds)
        CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,         &
                   ipds,igds, mtem,mitem1,mitem2)

        IF(istatus /= 0) THEN
           WRITE(6,*) 'Error when presure =', iprlvl(klev)
           STOP
        END IF

      END DO

!
!-----------------------------------------------------------------------
!
!  Put temperature into tem2
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

!-----------------------------------------------------------------------
!
!    Calculate sea level pressure (mb)
!    Reduction method: Benjamin and Miller: 1990, MWR, vol.118, No.10,
!                   Page: 2100-2101
!
!-----------------------------------------------------------------------
!
      gamma=.0065             ! std lapse rate per meter
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

      ! 1. Horizontally interpolate from ARPS grid to ETA grid

      CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,            &
                   xs,ys,x2d,y2d,tem1(:,:,1),tem2d_eta,               &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,  &
                   tem4)

      ! 2. Pack and write GRIB message

      CALL mkipds_212(02,102,00,0,p2time,0,0,ipds)
      CALL wrtgb(nx_eta,ny_eta,100*tem2d_eta,nchw,wdlen,grbpkbit,         &
                 ipds,igds, mtem,mitem1,mitem2)
!
!
!-----------------------------------------------------------------------
!
!  Calculate stability indices.
!  Use level k=2 as the "surface".
!
!-----------------------------------------------------------------------
!
!      imid=(nx/2)+1
!      jmid=(ny/2)+1
      CALL arps_be(nx,ny,nz,                                            &
           pprt,zps_km,tem2,qvprt,                                      &
           lcl,lfc,el,twdf,li,pbe,mbe,nbe,tcap,                         &
           wrk1,wrk2,wrk3,wrk4,wrk5,wrk6,                               &
           wrk7,wrk8,wrk9,wrk10,wrk11,wrk12,tem2d1)

!      PRINT *, ' Sample stability values: '
!      PRINT *, '   lcl,lfc : ',lcl(imid,jmid),lfc(imid,jmid)
!      PRINT *, '   el, twdf: ',el(imid,jmid),twdf(imid,jmid)
!      PRINT *, '   li, pbe : ',li(imid,jmid),pbe(imid,jmid)
!      PRINT *, '   mbe, nbe: ',mbe(imid,jmid),nbe(imid,jmid)
!      PRINT *, '   tcap    : ',tcap(imid,jmid)
!
      CALL calcshr(nx,ny,nz,x,y,zp,mf2d,                                &
          pprt,tem2,uprt,vprt,pbe,                                      &
          shr37,ustrm,vstrm,srlfl,srmfl,heli,brn,brnu,blcon,            &
          tem2d1,tem2d2,tem2d3,tem2d4,tem4)

!      PRINT *, ' Sample shear values: '
!      PRINT *, ' shr37,ustrm,vstrm: ',                                  &
!          shr37(imid,jmid),ustrm(imid,jmid),vstrm(imid,jmid)
!      PRINT *, ' srlfl,srmfl,heli: ',                                   &
!          srlfl(imid,jmid),srmfl(imid,jmid),heli(imid,jmid)
!      PRINT *, ' brn,brnu,blcon: ',                                     &
!          brn(imid,jmid),brnu(imid,jmid),blcon(imid,jmid)


!  process near-sfc LI
       PRINT *, ' Storing near-sfc LI'
       CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,             &
                   xs,ys,x2d,y2d,li,tem2d_eta,                          &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)

       CALL mkipds_212(24,01,00,0,p2time,0,0,ipds)
       CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,          &
                  ipds,igds, mtem,mitem1,mitem2)

!  process near-sfc CAPE
       PRINT *, ' Storing near-sfc CAPE'
       CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,             &
                   xs,ys,x2d,y2d,pbe,tem2d_eta,                         &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)

       CALL mkipds_212(157,01,00,0,p2time,0,0,ipds)
       CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,          &
                  ipds,igds, mtem,mitem1,mitem2)

!  process near-sfc CIN
       PRINT *, ' Storing near-sfc CIN'

       CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,             &
                   xs,ys,x2d,y2d,nbe,tem2d_eta,                         &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)

       CALL mkipds_212(156,01,00,0,p2time,0,0,ipds)
       CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,          &
                   ipds,igds, mtem,mitem1,mitem2)

!  process near-sfc Helicity
       PRINT *, ' Storing helicity'
       CALL mkarps2d(nx,ny,nx_eta,ny_eta, iorder,iloc,jloc,             &
                   xs,ys,x2d,y2d,heli,tem2d_eta,                        &
                   dxfld, dyfld,rdxfld, rdyfld,slopey,alphay, betay,    &
                   tem4)

       CALL mkipds_212(190,01,00,0,p2time,0,0,ipds)
       CALL wrtgb(nx_eta,ny_eta,tem2d_eta,nchw,wdlen,grbpkbit,          &
                  ipds,igds, mtem,mitem1,mitem2)

!-----------------------------------------------------------------------
       CLOSE(UNIT=nchw)
       CALL retunit(nchw)
!
!-----------------------------------------------------------------------
!
!  Create ready file, indicating history dump writing is complete
!
!-----------------------------------------------------------------------
!
       IF( readyfl == 1 ) THEN
         WRITE (filnamr,'(a)') trim(filname) // "_ready"
         CALL getunit( nchoutr )
         OPEN (UNIT=nchoutr,FILE=trim(filnamr))
         WRITE (nchoutr,'(a)') trim(filname)
         CLOSE (nchoutr)
         CALL retunit (nchoutr )
       END IF

    END IF

  END DO


  STOP
END PROGRAM arps2eta212

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
  INTEGER :: nx,ny,nz
  REAL :: zps(nx,ny,nz)
  REAL :: t(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  INTEGER :: iprlvl
  REAL :: hgtcdf(nx,ny)
!
  INCLUDE 'phycst.inc'
!
  REAL :: gamma,rddg,const
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          rddg  = (rd/g),                                               &
              const = (rd*gamma/g) )
!
  INTEGER :: i,j
  REAL :: prcdf
!
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

SUBROUTINE extrpt(nx,ny,nz,t,pr,iprlvl,tcdf)
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
  INTEGER :: nx,ny,nz
  REAL :: t(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  INTEGER :: iprlvl
  REAL :: tcdf(nx,ny)

  INCLUDE 'phycst.inc'

  REAL, PARAMETER :: gamma = 0.0065         ! degrees/m  lapse rate
  REAL, PARAMETER :: const = (rd*gamma/g)

  INTEGER :: i,j
  REAL :: prcdf

  prcdf = 100.*FLOAT(iprlvl)
  DO j=1,ny-1
    DO i=1,nx-1
      IF(prcdf <= pr(i,j,nz-1)) THEN
        IF(pr(i,j,nz-1) <= 30000.) THEN
          tcdf(i,j)=t(i,j,nz-1)
        ELSE
          tcdf(i,j)=t(i,j,nz-1)*                                        &
                    ((prcdf/pr(i,j,nz-1))**const)
        END IF
      ELSE IF(prcdf >= pr(i,j,2)) THEN
        tcdf(i,j)=t(i,j,2)*                                             &
                    ((prcdf/pr(i,j,2))**const)
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE extrpt
!
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
  INTEGER :: nx,ny,nz
  REAL :: q(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  INTEGER :: iprlvl
  REAL :: qcdf(nx,ny)
!
  INTEGER :: i,j
  REAL :: prcdf
!
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

SUBROUTINE extrpuv(nx,ny,nz,us,vs,pr,iprlvl,ucdf,vcdf)

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
  INTEGER :: nx,ny,nz
  REAL :: us(nx,ny,nz)
  REAL :: vs(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  INTEGER :: iprlvl
  REAL :: ucdf(nx,ny)
  REAL :: vcdf(nx,ny)
!
  INTEGER :: i,j
  REAL :: prcdf
!
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
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_avor                 ######
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
!  Calculate absolute Vort*10^5 (1/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

SUBROUTINE cal_avor(tem9,u,v,x,y,nx,ny,nz,mode,flagsin,omega,           &
           sinlat,tem1,tem2,tem3)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz),sinlat(nx,ny)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)
  REAL :: tem1(nx), tem2(ny), tem3(nx,ny)
  REAL :: omega
  INTEGER :: mode,flagsin

  INTEGER :: i,j,k
  REAL :: tmp1
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        tem9(i,j,k)= 1.0E5*(                                            &
               (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/       &
               (4*(x(i+1)-x(i)))-                                       &
               (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/       &
               (4*(y(j+1)-y(j))) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      tem9(i,j,   1)=tem9(i,j,   2)
      tem9(i,j,nz-1)=tem9(i,j,nz-2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-2
      tem9(   1,j,k)=tem9(   2,j,k)
      tem9(nx-1,j,k)=tem9(nx-2,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      tem9(i,   1,k)=tem9(i,   2,k)
      tem9(i,ny-1,k)=tem9(i,ny-2,k)
    END DO
  END DO

  IF(mode == 1) THEN
    IF( flagsin == 0) THEN
      CALL gtsinlat_ext(nx,ny,x,y,sinlat,tem1,tem2, tem3)
      tmp1 = 2.0* omega
      DO j=1,ny
        DO i=1,nx
          sinlat(i,j) = tmp1 * sinlat(i,j)*1.0E5
        END DO
      END DO
!      flagsin=1
    END IF
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem9(i,j,k) = tem9(i,j,k) + sinlat(i,j)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE cal_avor

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTSINLAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtsinlat_ext(nx,ny, x,y, sinlat, xs, ys, temxy)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the sin of the lattitude of each grid point, to be used
!  in the calculation of latitude-dependent Coriolis parameters.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  3/21/95.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!
!  OUTPUT:
!
!    sinlat   Sin of latitude at each grid point
!
!  WORK ARRAYS:
!
!    xs       x-coordinate at scalar points.
!    ys       y-coordinate at scalar points.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny          ! Dimensions of model grid.
  REAL :: x     (nx)        ! The x-coord. of the u points.
  REAL :: y     (ny)        ! The y-coord. of the v points.

  REAL :: sinlat(nx,ny)     ! Sin of latitude at each grid point

  REAL :: xs    (nx)        ! x-coord. of scalar points.
  REAL :: ys    (ny)        ! y-coord. of scalar points.

  REAL :: temxy (nx,ny)     ! Work array.

  REAL :: r2d
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  DO i=1,nx-1
!    xs(i) = (x(i)+x(i+1))*0.5
!  END DO
!  xs(nx) = -xs(nx-2)+2.0*xs(nx-1)
!
!  DO j=1,ny-1
!    ys(j) = (y(j)+y(j+1))*0.5
!  END DO
!  ys(ny) = -ys(ny-2)+2.0*ys(ny-1)
!
  CALL xytoll(nx,ny,xs,ys,sinlat,temxy)

  r2d = ATAN(1.0)*4.0/180.0

  DO j=1,ny
    DO i=1,nx
      sinlat(i,j) = SIN( r2d * sinlat(i,j) )
    END DO
  END DO

  RETURN
END SUBROUTINE gtsinlat_ext

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%                                                              %%
!%%  Subroutines for grib dumping                                %%
!%%                                                              %%
!%%  The follwing subroutines were all rewritten from those      %%
!%%  in src/arps/gribio3d.f90 specifically for ETA #212 grid     %%
!%%                                                              %%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MKIPDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkipds_212(parno,lvltyp,vlvl,timer,p1time,p2time,scaling, ipds)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make integer product data array ipds
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    parno    Indicator of parameter (see Code table 2)
!    lvltyp   Indicator of type of level (see Code table 3)
!    vlvl     Height, pressure, etc. of levels (see Code table 3)
!    timer    Time range indicator (eithor 10 or 4)
!    p1time
!    p2time   Period of time (number of time units)
!             ARPS forcast time in minutes
!    scaling  Units decimal scale factor
!
!  OUTPUT:
!
!    ipds     IPDS array
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: parno,lvltyp,vlvl,timer,p1time,p2time,scaling

  INTEGER, INTENT(OUT) :: ipds(25)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ipds(1) = 28          ! number of bytes in PDS section
  ipds(2) = 2           ! parameter table version number, unknown
            !%% ipds(19) = 2   !? (19)  - VERSION NR OF PARAMETER TABLE

  ipds(3) = 0           ! ID number of originating center. CAPS?
  ipds(4) = 0           ! ID number of model. ARPS/CAPS?
  ipds(5) = 212         ! grid ID, 255 for self-defined GDS
  ipds(6) = 1           ! GDS flag, 1 = GDS section included
            !%%  ipds(1) =  07  ! (1)   - ID OF CENTER
            !%%  ipds(2) = 255  ! (2)   - GENERATING PROCESS ID NUMBER
            !%%  ipds(3) = 212  ! (3)   - GRID DEFINITION
            !%%  ipds(4) = 128  ! (4)   - GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)
                                !         128   GDS included
                                !          64   BMS included
                                !         192   Both included
  ipds(7) = 0           ! BMS flag, 0 = no BMS section included

!%%%
!%%% After called this subroutine, ipds(8), ipds(9) and ipds(11)
!%%% need to be set for each variables specifically.
!%%%
  ipds(8) = parno       ! indicator of parameter and unit
                        ! (table 2), depends on variables
  ipds(9) = lvltyp      ! indicator of type of level (table 3)
                        ! depends on which variable
                        ! 100 for isobar surface
                        ! 103 for z coordinates in meters
                        ! 111 for soil layers in centimeters
  ipds(10) = 0          ! value 1 of level, N/A
  ipds(11) = vlvl       ! value 2 of level (value of level)
  ! both use two octets: ipds(10) & ipds(11)
           !%% ipds(5) = 011   ! (5)   - INDICATOR OF PARAMETER
           !%% ipds(6) = 100   ! (6)   - TYPE OF LEVEL
           !%% ipds(7) = 1000  ! (7)   - HEIGHT/PRESSURE , ETC OF LEVEL

  IF ( year <= 2000 ) THEN
    ipds(12) = year-1900 ! year of century
    ipds(23) = 20        ! century (20, change to 21 on Jan 1, 2001)
  ELSE
    ipds(12) = year-2000 ! year of century
    ipds(23) = 21        ! century (20, change to 21 on Jan 1, 2001)
  END IF

  ipds(13) = month      ! month of year
  ipds(14) = day        ! day of month
  ipds(15) = hour       ! hour of day
  ipds(16) = minute     ! minute of hour
            !%% ipds(8) = year-2000  ! (8)   - YEAR INCLUDING (CENTURY-1)
            !%% ipds(9) = month      ! (9)   - MONTH OF YEAR
            !%% ipds(10)= day        ! (10)  - DAY OF MONTH
            !%% ipds(11)= hour       ! (11)  - HOUR OF DAY
            !%% ipds(12)= minute     ! (12)  - MINUTE OF HOUR
            !%% ipds(21)= 21         ! (21)  - CENTURY OF REFERENCE TIME OF DATA

  ipds(17) = 1          ! forecast time unit, minutes
            !%% ipds(13) = 254  !  (13)  - INDICATOR OF FORECAST TIME UNIT
                                !       254 seconds
                                !         0 minutes
                                !         1 hours
  ipds(18) = p1time     ! forecast time P1, or current time, curtim
  ipds(19) = p2time     ! forecast time P2, N/A
  ipds(20) = timer      ! time range indicator,
                        ! 10 = use two octets from ipds(18)
                        !  0 = Analysis at initial time
                        !  4 = accumulated from p1 to p2
            !%% ipds(14) = 0    !  (14)  - TIME RANGE 1
            !%% ipds(15) = 0    !  (15)  - TIME RANGE 2
            !%% ipds(16) = 1    !  (16)  - TIME RANGE FLAG
  ipds(21) = 0          ! number include in average, N/A
  ipds(22) = 0          ! number missing from average, N/A
            !%% ipds(17) = 0    !  (17)  - NUMBER INCLUDED IN AVERAGE
            !%% ipds(20) = 0    !  (20)  - NR MISSING FROM AVERAGE/ACCUMULATION
  ipds(24) = 0          ! subcenter ID number
  ipds(25) = scaling    ! scaling power of 10,
                        ! depends on which variable
            !%% ipds(23) = 255  !  (23)  - SUBCENTER NUMBER
            !%% ipds(22) = 2    !? (22)  - UNITS DECIMAL SCALE FACTOR



! parameters below are those appeared in W3LIB

  !%% ipds(18) = 1   !? (18)  - VERSION NR OF GRIB SPECIFICATION
  !%% ipds(24) = 0   !  (24)  - PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
                     !             128 IF FORECAST FIELD ERROR
                     !             64 IF BIAS CORRECTED FCST FIELD
                     !             32 IF SMOOTHED FIELD
                     !      WARNING: CAN BE COMBINATION OF MORE THAN 1
  !%% ipds(25) = 0   !  (25)  - PDS BYTE 30, NOT USED

  RETURN
END SUBROUTINE mkipds_212
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MKIGDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkigds_212(nx,ny,nz,igds)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make integer grid description array IGDS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    igds     IPDS array
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny,nz ! Number of grid points in 3 directions

  INTEGER, INTENT(OUT) :: igds(25)

!-----------------------------------------------------------------------
!
!  ETA #212 Grid parameters
!
!-----------------------------------------------------------------------
  REAL, PARAMETER :: swlat = 12.19, swlon = -133.459
  REAL, PARAMETER :: trulat = 25.0, trulon = -95.0
  REAL, PARAMETER :: dx = 40635.25, dy = 40635.25  ! Meters
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

    igds( 1) = nz       ! number of vertical (z) coordinates
    igds( 2) = 255      ! location (in octets) of vertical
                        ! coordinate list
    igds( 3) = 3        ! data representation type (table 6)
    igds( 4) = nx       ! number of x-coordinates
    igds( 5) = ny       ! number of y-coordinates
    igds( 6) = nint(swlat*1000)  ! latitude  at southwest corner
    igds( 7) = nint(swlon*1000)  ! longitude at southwest corner
    igds( 8) = 8        ! u & v relative to x & y direction
    igds( 9) = nint(trulon*1000) ! true longitude of projection
    igds(10) = NINT(dx)
    igds(11) = NINT(dy)

    igds(12) = 0        ! northern pole on plane

    igds(13) = 64       ! scanning flag, 64 for +i & +j direction
    igds(14) = 0        ! unused
    igds(15) = nint(trulat*1000)  ! first true latitude (millidegree)
    igds(16) = nint(trulat*1000)  ! second true latitude (millidegree)
    igds(17) = -90000   ! latitude  of south pole (millidegree)
    igds(18) = nint(trulon*1000)   ! longitude of south pole (millidegree)

    igds(19:25) = 0       ! unused

  RETURN
END SUBROUTINE mkigds_212
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE WRTGB                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtgb(nx,ny,fvar,nchanl,wdlen,grbpkbit,              &
           ipds,igds,tem1,item1,item2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write a 2D float array, fvar, into file pointer nchanl
!  which points to a GRIB file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    fvar     Floating array to be written into GRIB file
!
!    nchanl   FILE pointer of GRIB stream file which was opened by
!             a C program, GOPEN
!
!    wdlen    Length of machine word (i.e. size of integer) in bytes
!
!    ipds     Integer array to carry the parameters for generating
!             PDS (Section 1).
!    igds     Integer array to carry the parameters for generating
!             GDS (Section 3).
!    ibdshd   BDS header
!
!  WORK ARRAY:
!
!    pds      PDS (GRIB Section 1)
!    gds      GDS (GRIB Section 3)
!
!    bds      BDS (GRIB Section 4)
!    ibds     Identical to BDS
!
!    nbufsz   Size of buffer of a GRIB message array
!    mgrib    Working buffer of GRIB message array
!
!
!    tem1     Temporary work array.
!    item1    Integer temporary work array.
!    item2    Integer temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny             ! Number of grid points in x,y,z dir.

  REAL :: fvar(nx,ny)          ! 2D Floating array

  INTEGER :: nchanl            ! FILE pointer indicates the GRIB file
  INTEGER :: wdlen             ! Length of machine word
  INTEGER :: grbpkbit          ! Number of bits in packing GRIB data

  INTEGER :: ipds(25)          ! ipds
  INTEGER :: igds(25)          ! igds
  INTEGER :: ibdshd(4)         ! BDS header

  REAL :: tem1 (nx,ny)         ! Temporary work array

  INTEGER :: item1(nx,ny)      ! Integer temporary work array
  INTEGER :: item2(nx,ny)      ! Integer temporary work array
!
!-----------------------------------------------------------------------
!
!  GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: msglen    ! Length of GRIB message
  INTEGER :: npts

  CHARACTER (LEN=1) :: pds(28)       ! PDS
  CHARACTER (LEN=1) :: gds(42)       ! GDS

  INTEGER, PARAMETER :: nbufsz = 800000     ! Size of GRIB message array
  INTEGER :: ibds(nbufsz/4)          ! Identical to BDS
  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS
  EQUIVALENCE (bds, ibds)

  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Working buffer

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  INTEGER :: itype             ! Data type: 0 - floating, 1 - integer
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  npts = nx*ny

  itype = 0
  item1 = 0
  ibdshd = 0

  CALL gribenc(itype,wdlen,grbpkbit,npts,fvar,item1,                  &
               ipds,igds,ibdshd,pds,gds,                              &
               nbufsz,bds,ibds,msglen,mgrib,                          &
               tem1,item2)

  WRITE (nchanl) (mgrib(i),i=1,msglen)

  RETURN
END SUBROUTINE wrtgb
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION compute_density            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION compute_density(t_k, p_pa) RESULT(rho)

! Computes density (kg/m{-3}) given temperature (K) and pressure (Pa)
  IMPLICIT NONE
  REAL, INTENT(IN) :: t_k
  REAL, INTENT(IN) :: p_pa
  REAL             :: rho

  REAL, PARAMETER  :: Rd =  287.04  ! J deg{-1} kg{-1}
                                    ! Gas constant for dry air

  rho = p_pa / (Rd * t_k)

  RETURN
END FUNCTION compute_density
