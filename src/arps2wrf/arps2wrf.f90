PROGRAM arps2wrf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARPS2WRF                    ######
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
!  Converts ADAS analysis file and/or EXT2ARPS dumps in ARPS grid to
!  the corresponding WRF initialization file (wrfinput_d01) and
!  the lateral boundary file (wrfbdy_d01)
!
!  This program reads in a history file produced by ADAS and/or EXT2ARPS
!  in any ARPS history format (except HDF 4 format), interpolates
!  variables to WRF grid specified in the namelist file - arps2wrf.input
!  and writes out these variables in NetCDF/PHDF5 format as wrfinput_d01
!  and wrfbdy_d01
!
!  Users can specify two options for the surface characteristics, such as
!  vegetation type, soil type, and greenness fraction etc.
!
!  Option 1 sfcinitopt = "ARPS", this program reads surface variables from
!           the file created by ARPSSFC. the data must be at the same
!           grid as that in the ADAS/EXT2ARPS data file.
!
!         o  Variables read in are: "soiltyp", "vegtyp", "veg";
!         o  Variable, "xice", "xland", will be determined after the above
!            variables are interpolated to WRF grid.
!         o  This option cannot set "SHDMAX","SHDMIN","SNOALB","ALBBCK".
!            Should study whether it has negative effect to WRF run???
!
!  Option 2 sfcinitopt = "WRFSI", this program reads those variables from
!           the static file 'static.wrfsi' created by gridgen_model of
!           WRFSI package. The data must be at the same grid as those
!           specified for WRF model in the input file, arps2wrf.input.
!
!         (1). Users only need to run "window_domain_rt.pl" to generate
!              WRFSI static file. (Actually, only gridgen_model.exe component).
!         (2). One static file works for all runs in the same domain.
!
!
!  INPUT
!    arps2wrf.input     NAMELIST
!    runname.bin000000  ADAS time-dependent variable file
!    runname.bingrdbas  ADAS base state file
!    runname.sfcdata    ARPSSFC outputs
!    OR
!    static.wrfsi       WRFSI static file
!
!    OR
!    outputs from ext2arps for lateral boundary conditions.
!
!  OUTPUT
!    wrfinput_d01
!    wrfbdy_d01
!    wrfnamelist.input
!
!  Libraries:
!    libarps.a
!    libadas.a
!    NetCDF library
!
!  Subroutine calls:
!    dtaread
!    Subroutines defined in file interplib.f90 and wrf_iolib.f90
!
!  Use module
!    wrf_metadata       Defined in module_wrf_metadata.f90
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang, Dan Weber, Keith Brewster
!  07/08/2003
!
!  MODIFICATION HISTORY:
!
!  03/10/2004 Yunheng Wang
!  Upgraded to support WRF version 1.4 which was downloaded from
!  WRF tiger team website on Feb. 10, 2004.
!
!  07/28/2004 Yunheng Wang
!  Upgraded WRF V1.4 support to V2.0 support. So the unofficial release
!  WRF V1.4 would not be supported any more since ARPS5.1.3.
!
!  11/20/2004 Yunheng Wang
!  Reorangized to support PHDF5 format using WRF external IO-PHDF5
!  package. Upgraded to support WRFV2.0.3. Removed supports for
!  previous WRF version.
!
!  01/10/2005 Yunheng Wang
!  Added WRF internal binary support. Please note that WRF binary does
!  not support random access. So the dumping order is important for
!  binary data. The current code support WRFV2.0.3.1. The compatiblity
!  between WRF versions will be bad, but it should not affect other
!  IO format.
!
!-----------------------------------------------------------------------
!
!  ARPS DATA ARRAYS
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
!    qscalar  Hydrometeor scalars
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
  USE wrf_metadata             ! WRF constants and metadata

  IMPLICIT NONE
!----------------------------------------------------------------------
!
! ARPS grid variables
!
!---------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Grid dimensions for ARPS.
  INTEGER :: nzsoil            ! Soil levels
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER, PARAMETER :: nhisfile_max = 100, nmax_domains = 100
  INTEGER, PARAMETER :: max_vertical_levels = 100

  CHARACTER(LEN=256) :: hisfile(nhisfile_max)
  CHARACTER(LEN=256) :: bdybasfn(nhisfile_max)
  INTEGER            :: nhisfile,lengbf,lenfil,lenfil0
  CHARACTER(LEN=256) :: grdbasfn

  CHARACTER(LEN=256), DIMENSION(nmax_domains) :: adashisfn
  CHARACTER(LEN=256), DIMENSION(nmax_domains) :: adasbasfn
  INTEGER,            DIMENSION(nmax_domains) :: hinfmt
  INTEGER                                     :: max_dom
  LOGICAL,            DIMENSION(nmax_domains) :: input_from_file

  INTEGER            :: finexist(nhisfile_max)

  ! hisfile(1) and bdybasfn(1) are for WRF input file
  !      They can be any ARPS history dumps including ADAS output,
  !      ARPS output or EXT2ARPS output
  !
  ! hisfile(2:nhisfil), bdybasfn(2:nhisfile) are for
  !      WRF lateral bounday files. They can be either ARPS history
  !      dumps or EXT2ARPS outputs

  CHARACTER(LEN=80) :: execname
!
!-----------------------------------------------------------------------
!
!  ARPS include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
! Variables for mpi jobs
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: fzone_arps = 3, fzone_wrf = 1
  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER :: nxsm, nysm             ! smaller domain
  INTEGER :: nxlg, nylg             ! global domain
  INTEGER :: nxlg_wrf, nylg_wrf

  INTEGER :: nprocx_in, nprocy_in
  INTEGER :: ii,jj,ia,ja

  REAL, DIMENSION(:),        POINTER :: xsm,ysm
  REAL, DIMENSION(:,:,:),    POINTER :: zpsm,zpsoilsm
  REAL, DIMENSION(:,:,:),    POINTER :: uprtsm, vprtsm, wprtsm,            &
                                        ptprtsm, pprtsm, qvprtsm
  REAL, DIMENSION(:,:,:,:),  POINTER :: qscalarsm
  REAL, DIMENSION(:,:,:),    POINTER :: ubarsm, vbarsm, wbarsm,            &
                                        ptbarsm,pbarsm, qvbarsm,rhobarsm

  INTEGER, DIMENSION(:,:,:), POINTER :: soiltypsm
  INTEGER, DIMENSION(:,:),   POINTER :: vegtypsm
  REAL, DIMENSION(:,:,:),    POINTER :: stypfrctsm(:,:,:)
  REAL, DIMENSION(:,:,:,:),  POINTER :: tsoilsm, qsoilsm
  REAL, DIMENSION(:,:,:),    POINTER :: wetcanpsm
  REAL, DIMENSION(:,:),      POINTER :: vegsm,snowdpthsm
!
!-----------------------------------------------------------------------
!
!  ARPS arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, DIMENSION(:), POINTER :: x  ! The x-coord. of the physical and
                           ! computational grid. Defined at u-point.
  REAL, DIMENSION(:), POINTER :: y  ! The y-coord. of the physical and
                           ! computational grid. Defined at v-point.
  REAL, DIMENSION(:), POINTER :: z  ! The z-coord. of the computational
                           ! grid. Defined at w-point.

  REAL, DIMENSION(:,:,:), POINTER :: zp     ! The height of the terrain.
  REAL, DIMENSION(:,:,:), POINTER :: zpsoil ! The height of the soil model.

  REAL, DIMENSION(:,:,:), POINTER :: uprt   ! Perturbation u-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: vprt   ! Perturbation v-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: wprt   ! Perturbation w-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: ptprt  ! Perturbation potential temperature (K)
  REAL, DIMENSION(:,:,:), POINTER :: pprt   ! Perturbation pressure (Pascal)
  REAL, DIMENSION(:,:,:), POINTER :: qvprt  ! Perturbation water vapor specific
                                            ! humidity (kg/kg)
  REAL, DIMENSION(:,:,:,:), POINTER :: qscalar

  REAL, DIMENSION(:,:,:), POINTER :: qv     ! Water vapor specific humidity (kg/kg)

  REAL, DIMENSION(:,:,:), POINTER :: ubar   ! Base state u-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: vbar   ! Base state v-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: wbar   ! Base state w-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: ptbar  ! Base state potential temperature (K)
  REAL, DIMENSION(:,:,:), POINTER :: pbar   ! Base state pressure (Pascal)
  REAL, DIMENSION(:,:,:), POINTER :: rhobar ! Base state air density (kg/m**3)
  REAL, DIMENSION(:,:,:), POINTER :: qvbar  ! Base state water vapor specific

  INTEGER, DIMENSION(:,:,:), POINTER :: soiltyp  ! Soil type
  REAL,    DIMENSION(:,:,:), POINTER :: stypfrct ! Soil type
  INTEGER, DIMENSION(:,:),   POINTER :: vegtyp   ! Vegetation type
  REAL,    DIMENSION(:,:),   POINTER :: veg      ! Vegetation fraction

  REAL, DIMENSION(:,:,:,:), POINTER :: tsoil    ! Soil Temperature (K)
  REAL, DIMENSION(:,:,:,:), POINTER :: qsoil    ! Soil Moisture
  REAL, DIMENSION(:,:,:),   POINTER :: wetcanp  ! Canopy water amount
  REAL, DIMENSION(:,:),     POINTER :: snowdpth ! Snow depth (m)

  REAL, DIMENSION(:,:),     POINTER :: dum2da
  REAL, DIMENSION(:,:,:),   POINTER :: dum3da

!-----------------------------------------------------------------------
!
!  Namelist definitions for ARPS2WRF.input
!
!     sfcdt              Specify surface characteristics
!     bdyspc             Obtain boundary input files (ARPS format)
!     wrf_grid           Define WRF horizontal and vertical grid
!     interp_options     Choose interpolation scheme
!     wrf_opts           WRF options from namelist.input
!     output             Ouput options
!
!-----------------------------------------------------------------------
!
  INTEGER, DIMENSION(nmax_domains) :: nxnm, nynm,nznm,nzsoilnm,nstypsnm
  INTEGER, DIMENSION(nmax_domains) :: nx_wrfnm, ny_wrfnm
  REAL,    DIMENSION(nmax_domains) :: dx_wrfnm, dy_wrfnm

  INTEGER, DIMENSION(nmax_domains) :: i_parent_start, j_parent_start
  INTEGER, DIMENSION(nmax_domains) :: parent_id

  INTEGER           :: nx_wrf         ! = nx-2 if the same grid as ARPS
  INTEGER           :: ny_wrf         ! = ny-2 if the same grid as ARPS
  INTEGER           :: nz_wrf         ! = nz-2 if the same grid as ARPS
                                      ! All are staggered values
  INTEGER           :: nzsoil_wrf     ! WRF number of soil layers
                                      ! see sf_surface_physics
  REAL              :: lattru_wrf(2)  ! array of true latitude of WRF map projection
  REAL              :: lontru_wrf     ! true longitude of WRF map projection
                                      ! = trulon_wrf

! Namelist variable declaration

  CHARACTER(LEN=7), DIMENSION(nmax_domains)  :: sfcinitopt      ! either "ARPS" or "WRFSI"

  INTEGER           :: create_bdy      ! Create WRF boundary file
  INTEGER           :: create_namelist ! Dump WRF namelist.input
  INTEGER, DIMENSION(nmax_domains) :: use_arps_grid   ! Use ARPS horizontal grid as WRF grid

  INTEGER           :: tintv_bdywrf    ! The desired WRF boundary file interval
  INTEGER           :: tintv_bdyin     ! ARPS boundary file interval (in seconds)
  INTEGER           :: mgrdbas         ! Options for grid base file
                                       ! = 0 share same grid base as initial state file
                                       ! = 1 All ARPS boundary files share one grd base
                                       !     file but it is difference from the inital
                                       !     base file as specified using grdbasfn
                                       ! = 2 Each file has its own grid base file

  CHARACTER(LEN=256):: wrfnamelist     ! file name for WRF namelist.input
  REAL              :: wrfversion

  INTEGER :: mapproj_wrf     ! Type of map projection in WRF model grid
                             ! modproj = 1  Polar Stereographic
                             !              projection.
                             !         = 2  Mercator projection.
                             !         = 3  Lambert projection.

  REAL    :: sclfct_wrf      ! Map scale factor.
                             ! Distance on map, between two latitudes
                             ! trulat1 and trulat2,
                             ! is = (Distance on earth)*sclfct.
                             ! For ARPS model runs,
                             ! generally this is 1.0

  REAL    :: trulat1_wrf, trulat2_wrf, trulon_wrf
                             ! 1st, 2nd real true latitude and true longitude
                             ! of WRF map projection
  REAL    :: ctrlat_wrf      ! Center latitude of WRF model domain (deg. N)
  REAL    :: ctrlon_wrf      ! Center longitude of WRF model domain (deg. E)

  REAL    :: dx_wrf          ! WRF Grid spacing in x-direction
  REAL    :: dy_wrf          ! WRF Grid spacing in y-direction

  REAL, DIMENSION(nmax_domains) :: ptop            ! WRF atmosphere top pressure in Pascal
  REAL    :: zlevels_wrf(max_vertical_levels)
                             ! WRF mass levels from 1.0 at surfact to
                             ! 0.0 at atmosphere top

  INTEGER :: iorder          ! order of polynomial for horizontal
                             ! interpolation (1 or 2)
  INTEGER :: korder          ! vertical interpolation order (1 or 2)

  INTEGER :: dyn_opt         ! WRF dynamics option
                             ! only works for = 2 Eulerian mass coordinate
  INTEGER :: diff_opt        ! WRF diffusion option
  INTEGER :: km_opt          ! WRF eddy coefficient option
  REAL    :: dampcoef
  REAL,   DIMENSION(nmax_domains) :: khdif           ! Horizontal diffusion constant (m^2/s)
  REAL,   DIMENSION(nmax_domains) :: kvdif           ! Vertical diffusion constant (m^2/s)
  INTEGER,DIMENSION(nmax_domains) :: mp_physics      ! WRF microphysics options
                             != 2 Lin et al. scheme
                             !   (QVAPOR,QRAIN,QSNOW,QCLOUD,QICE,QGRAUP)
                             != 3 NCEP 3-class simple ice scheme
                             !   (QVAPOR,QCLOUD,QICE,QRAIN,QSNOW)
                             != 4 NCEP 5-class scheme
                             !   (QVAPOR,QCLOUD,QICE,QRAIN,QSNOW)
                             != 5 Ferrier (new Eta) microphysics
                             !   (QVAPOR,QCLOUD)
  INTEGER, DIMENSION(nmax_domains) :: ra_lw_physics   ! Longwave radiaiton option
  INTEGER, DIMENSION(nmax_domains) :: ra_sw_physics   ! Shortwave radiation option
  INTEGER, DIMENSION(nmax_domains) :: sf_sfclay_physics  ! WRF surface-layer option
  INTEGER, DIMENSION(nmax_domains) :: sf_surface_physics ! WRF land-surface option
                                ! = 0 no land-surface
                                !     (DO NOT use)
                                ! = 1 Thermal diffusion scheme
                                !     (nzsoil_wrf = 5)
                                ! = 2 OSU land-surface model
                                !     (nzsoil_wrf = 4)
                                ! = 3 Do not use
                                !     (nzsoil_wrf = 6)
  INTEGER, DIMENSION(nmax_domains) :: bl_pbl_physics     ! boundary-layer option
  INTEGER, DIMENSION(nmax_domains) :: cu_physics         ! cumulus option
  REAL    :: dt                 ! time-step for advection
  INTEGER :: spec_bdy_width     ! number of rows for specified boundary values nudging
  INTEGER :: nprocx_wrf      ! Number of X direction processors for WRF run
  INTEGER :: nprocy_wrf      ! Number of Y direction processors for WRF run
  INTEGER, DIMENSION(nmax_domains) :: frames_per_outfile
  INTEGER :: restart_interval
  REAL,    DIMENSION(nmax_domains) :: radt,cudt
  INTEGER, DIMENSION(nmax_domains) :: parent_time_step_ratio
  INTEGER :: ifsnow, w_damping
  INTEGER :: io_form, io_form_history, io_form_restart
  INTEGER, DIMENSION(nmax_domains) :: history_interval
  INTEGER :: intval2d
  CHARACTER(LEN=40) :: dir2d
  !LOGICAL :: pd_moist
  INTEGER :: moist_adv_opt
  INTEGER :: qx_zero_out
  CHARACTER(LEN=256) :: staticdir, indir, outdir

  CHARACTER(LEN=4), PARAMETER :: fmtstr(7) = (/ 'bin_','xxxx','xxxx','xxxx','HDF5','xxxx','ncdf'/)
!
!-----------------------------------------------------------------------
!
!  ARPS working arrays
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: xs(:),ys(:)    ! ARPS coordinate at scalar points
  REAL, ALLOCATABLE :: zps(:,:,:)     ! ARPS vertical coordinate at scalar points

  REAL, ALLOCATABLE :: temp(:,:,:)    ! ARPS temperature

  REAL, ALLOCATABLE :: tem1(:,:,:), tem2(:,:,:), tem3(:,:,:)

  REAL :: amax,amin

!-----------------------------------------------------------------------
!
!  WRF grid related variables
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: lat_wrf(:,:,:) ! WRF grid point latitudes
  REAL, ALLOCATABLE :: lon_wrf(:,:,:) ! WRF grid point lontitudes
  REAL, ALLOCATABLE :: msft_wrf(:,:)  ! Map scale factor on mass grid
  REAL, ALLOCATABLE :: msfu_wrf(:,:)  ! Map scale factor on u grid
  REAL, ALLOCATABLE :: msfv_wrf(:,:)  ! Map scale factor on v grid

  REAL, ALLOCATABLE :: zlevels_half(:) ! WRF vertical mass points at half levels
  REAL, ALLOCATABLE :: zs_wrf(:)      ! Depth of center of soil layers
  REAL, ALLOCATABLE :: dzs_wrf(:)     ! Thickness of soil layers

!-----------------------------------------------------------------------
!
!  Horizontally interpolation arrays
!
!-----------------------------------------------------------------------

  REAL,    ALLOCATABLE :: x2d(:,:,:)
  REAL,    ALLOCATABLE :: y2d(:,:,:)  ! WRF coordinate at ARPS grid
  INTEGER, ALLOCATABLE :: iloc(:,:,:) ! i index of WRF point in ARPS array.
  INTEGER, ALLOCATABLE :: jloc(:,:,:) ! j index of WRF point in ARPS array.

  REAL,    ALLOCATABLE :: dxfld(:,:)  ! interpolation arrays calculated
  REAL,    ALLOCATABLE :: rdxfld(:,:) ! in advance to speed up the
  REAL,    ALLOCATABLE :: dyfld(:,:)  ! interpolation
  REAL,    ALLOCATABLE :: rdyfld(:,:)
!
!-----------------------------------------------------------------------
!
!  WRF work arrays
!  (These arrays defined at WRF grid horizontally and
!   vertically at ARPS physical heights)
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: t_tmp(:,:,:)
  REAL, ALLOCATABLE :: p_tmp(:,:,:)
  REAL, ALLOCATABLE :: qv_tmp(:,:,:)
  REAL, ALLOCATABLE :: zps_tmp(:,:,:)

  REAL, ALLOCATABLE :: etap(:,:,:,:)     ! air mass at ARPS physical height vertically

  REAL, ALLOCATABLE :: tsoil_tmp(:,:,:)  ! Defined at ARPS soil layers
  REAL, ALLOCATABLE :: qsoil_tmp(:,:,:)
  REAL, ALLOCATABLE :: zpsoil_tmp(:,:,:)

  REAL, ALLOCATABLE :: work3d(:,:,:)
  REAL, ALLOCATABLE :: work2d(:,:)

!-----------------------------------------------------------------------
!
!  WRF ARRAYS
!  (Those are arrays defined at WRF grid both horizontally
!   and vetically)
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: mu(:,:)            ! WRF air mass at surface
  REAL, ALLOCATABLE :: hterain_wrf(:,:)   ! WRF topograph

  REAL, ALLOCATABLE :: u_wrf(:,:,:)   ! u-velocity (m/s)
  REAL, ALLOCATABLE :: v_wrf(:,:,:)   ! v-velocity (m/s)
  REAL, ALLOCATABLE :: w_wrf(:,:,:)   ! w-velocity (m/s)
  REAL, ALLOCATABLE :: ph_wrf(:,:,:)  ! perturbation geopotential
  REAL, ALLOCATABLE :: phb_wrf(:,:,:) ! base state geopotential
  REAL, ALLOCATABLE :: pot_wrf(:,:,:) ! total potential. temp.
  REAL, ALLOCATABLE :: pt_wrf(:,:,:)  ! perturbation potential. temp.
  REAL, ALLOCATABLE :: pt_init_wrf(:,:,:)
  REAL, ALLOCATABLE :: qv_wrf(:,:,:)  ! Water vapor mixing ratio
  REAL, ALLOCATABLE :: p_wrf(:,:,:)   ! perturbation pressure
  REAL, ALLOCATABLE :: pb_wrf(:,:,:)  ! reference pressure
  REAL, ALLOCATABLE :: qscalar_wrf(:,:,:,:)
  REAL, ALLOCATABLE :: mup_wrf(:,:)   ! perturbation air mass
  REAL, ALLOCATABLE :: mub_wrf(:,:)   ! base air mass

  REAL, ALLOCATABLE :: tem1_wrf(:,:,:) ! WRF temporary arrays
  REAL, ALLOCATABLE :: tem2_wrf(:,:,:)
  REAL, ALLOCATABLE :: tem3_wrf(:,:,:)
  REAL, ALLOCATABLE :: tem4_wrf(:,:,:)

  INTEGER, ALLOCATABLE :: soiltyp_wrf(:,:),vegtyp_wrf(:,:)
  REAL,    ALLOCATABLE :: vegfrct_wrf(:,:)
  REAL,    ALLOCATABLE :: xice(:,:),xland(:,:)  ! flags why they are reals?
  REAL,    ALLOCATABLE :: sst_wrf(:,:)
  REAL,    ALLOCATABLE :: hgt_wrf(:,:),tmn_wrf(:,:)
  REAL,    ALLOCATABLE :: shdmin(:,:),shdmax(:,:),albbck(:,:),snoalb(:,:)

  REAL, ALLOCATABLE :: tslb_wrf(:,:,:)
  REAL, ALLOCATABLE :: smois_wrf(:,:,:)
  REAL, ALLOCATABLE :: snowh_wrf(:,:)
  REAL, ALLOCATABLE :: canwat_wrf(:,:)

!----------------------------------------------------------------------
!
! OUTPUT work arrays
!
!----------------------------------------------------------------------

  TYPE(wrf_global_metadata) :: global_meta

  REAL, ALLOCATABLE :: uatv_wrf(:,:,:)
  REAL, ALLOCATABLE :: vatu_wrf(:,:,:)

  REAL      :: lat_ll(4), lat_lr(4), lat_ul(4), lat_ur(4)
  REAL      :: lon_ll(4), lon_lr(4), lon_ul(4), lon_ur(4)

!-----------------------------------------------------------------------
!
!  Boundary arrays
!
!-----------------------------------------------------------------------

  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: ubdy3dtemp1, vbdy3dtemp1,  &
                                    tbdy3dtemp1, pbdy3dtemp1, qbdy3dtemp1
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET :: mbdy2dtemp1

  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: ubdy3dtemp2, vbdy3dtemp2,  &
                                     tbdy3dtemp2,pbdy3dtemp2, qbdy3dtemp2
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET :: mbdy2dtemp2

  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: bdy3d1temp1,bdy3d1temp2,   &
                                     bdy3d1temp3,bdy3d1temp4,bdy3d1temp5
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET :: bdy2d1temp1

  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: bdy3d2temp1,bdy3d2temp2,   &
                                     bdy3d2temp3,bdy3d2temp4,bdy3d2temp5
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET :: bdy2d2temp1

  REAL, DIMENSION(:,:,:), POINTER     :: ubdy3dtempc,vbdy3dtempc,       &
                                         tbdy3dtempc,pbdy3dtempc,qbdy3dtempc
  REAL, DIMENSION(:,:),   POINTER     :: mbdy2dtempc

  REAL, DIMENSION(:,:,:), POINTER     :: ubdy3dtempcn,vbdy3dtempcn,       &
                                         tbdy3dtempcn,pbdy3dtempcn,qbdy3dtempcn
  REAL, DIMENSION(:,:),   POINTER     :: mbdy2dtempcn

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: bdys,bdyn,bdyw,bdye

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: blgs,blgn,blgw,blge


!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,n,ifile, nq
  INTEGER :: istatus, lenstr, ireturn
  INTEGER :: abstime,  abstime1
  INTEGER :: abstimep, abstimec, abstimen, abstimecn
  REAL    :: time, gmthr
  INTEGER :: domid

  INTEGER :: abstdiff1, abstdiff2, abstdiff
  LOGICAL :: hinterp_needed
  LOGICAL :: fexist
  INTEGER :: count_bdy,count_in_loop

  REAL    :: w1, w2
  REAL    :: latnot(2),ctrx, ctry, scalef
  REAL    :: swx,     swy
  REAL    :: swx_wrf, swy_wrf

  CHARACTER(LEN=256), DIMENSION(nmax_domains) :: sfcdtfns
  CHARACTER(LEN=256), DIMENSION(nmax_domains) :: cdummy
  INTEGER, DIMENSION(nmax_domains) :: wrftrnopt

  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=24)  :: times_str
  CHARACTER(LEN=24)  :: nextbdytimes
  INTEGER            :: nchr, ncid, nchout

  REAL               :: rdummy

!  INTEGER :: i_parent_end, j_parent_end

  INTEGER :: grid_ratio

  INTEGER, DIMENSION(nmax_domains) :: parent_grid_ratio
  REAL,    DIMENSION(nmax_domains) :: xsub0, ysub0

  INTEGER, DIMENSION(nmax_domains) :: year1,month1,day1,hour1,minute1,second1
  INTEGER                          :: year2,month2,day2,hour2,minute2,second2

  !WDT 2004-01-20 GMB: set SST
!  REAL      :: cnt,sst_sum,dist_weight

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL mpinit_proc(0)
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(14(/5x,a),/)')                              &
      '###############################################################',&
      '###############################################################',&
      '####                                                       ####',&
      '####              Welcome to ARPS2WRF                      ####',&
      '####                                                       ####',&
      '####  Program that reads in history files generated by     ####',&
      '####   ADAS/ARPS and produces files for WRF to start:      ####',&
      '####                                                       ####',&
      '####       wrfinput_d01, wrfbdy_d01, namelist.input        ####',&
      '####                                                       ####',&
      '###############################################################',&
      '###############################################################'

!-----------------------------------------------------------------------
!
! First, read in namelist input
!
!-----------------------------------------------------------------------

  dyn_opt = 2

  CALL readnamelist(0,max_dom,input_from_file,                          &
     hinfmt,adasbasfn,adashisfn,bdybasfn,hisfile,nhisfile,finexist,     &
     nxnm,nynm,nznm,nzsoilnm,nstypsnm,nprocx_in,nprocy_in,              &
     ncompressx,ncompressy,                                             &
     use_arps_grid,nx_wrfnm,ny_wrfnm,nz_wrf,zlevels_wrf,ptop,           &
     i_parent_start,j_parent_start,parent_id,                           &
     mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf,                      &
     ctrlat_wrf,ctrlon_wrf,dx_wrfnm,dy_wrfnm,dt,                        &
     base_pres,base_temp,base_lapse,iso_temp,                           &
     sfcinitopt,wrftrnopt,sfcdtfns,cdummy,cdummy,rdummy,rdummy,         &
     create_bdy,mgrdbas,tintv_bdywrf,tintv_bdyin,spec_bdy_width,        &
     diff_opt,km_opt,khdif,kvdif,mp_physics,ra_lw_physics,              &
     ra_sw_physics,sf_sfclay_physics,sf_surface_physics,bl_pbl_physics, &
     cu_physics, nprocx_wrf,nprocy_wrf,frames_per_outfile,              &
     restart_interval,radt,cudt,ifsnow,w_damping,parent_time_step_ratio,&
     iorder,korder,io_form,qx_zero_out,create_namelist,wrfnamelist,     &
     io_form_history, io_form_restart,history_interval,                 &
     indir,outdir,staticdir,intval2d,dir2d,moist_adv_opt,               &
     dampcoef,wrfversion,istatus)

  filename = ' '
  execname = ' '
  IF (mp_opt == 0) THEN
    execname = 'ARPS2WRF (serial)'
  ELSE
    WRITE(execname,'(a,2(I4.4,a))') 'ARPS2WRF (Parallel, ',nproc_x,'x',nproc_y,')'
  END IF

!-----------------------------------------------------------------------
!
! Begin to loop for each domain
!
!-----------------------------------------------------------------------

  DO domid = 1,max_dom

    IF (input_from_file(domid) .OR. domid == 1) THEN

    nx = nxnm(domid)
    ny = nynm(domid)
    nz = nznm(domid)
    nzsoil = nzsoilnm(domid)
    nstyps = nstypsnm(domid)

    nx_wrf = nx_wrfnm(domid)
    ny_wrf = ny_wrfnm(domid)
    dx_wrf = dx_wrfnm(domid)
    dy_wrf = dy_wrfnm(domid)

    IF (domid > 1) THEN
      create_bdy      = 0

      nhisfile        = 1
      hisfile (1) = adashisfn(domid)
      bdybasfn(1) = adasbasfn(domid)
    END IF

    WRITE(sfcdtfl,'(a)') sfcdtfns(domid)

!-----------------------------------------------------------------------
!
!  Setting microphysics variables
!
!-----------------------------------------------------------------------

   CALL set_mp_physics_variables(domid,mp_physics(domid),istatus)
   IF (istatus /= 0) CALL arpsstop('ERROR in set_mp_physics_variables.',1)
!
!-----------------------------------------------------------------------
!
! Allocate ARPS arrays
!
!-----------------------------------------------------------------------

!write(0,*) 'Deallocating'
  IF (ASSOCIATED(x)) THEN
    DEALLOCATE(x,y,z,zp,zpsoil)
    DEALLOCATE(uprt,vprt,wprt,ptprt,pprt,qvprt)
    DEALLOCATE(qv,qscalar )
    DEALLOCATE(ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar)
    DEALLOCATE(soiltyp,stypfrct,vegtyp,veg)
    DEALLOCATE(tsoil,qsoil,wetcanp,snowdpth)
    DEALLOCATE(dum2da, dum3da)
  END IF

!write(0,*) 'reallocateing'
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
  ALLOCATE(qscalar(nx,ny,nz,nscalar))
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
  ALLOCATE(veg    (nx,ny))

  ALLOCATE(tsoil  (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil  (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(dum2da(nx,ny))
  ALLOCATE(dum3da(nx,ny,nz))

  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0
  zpsoil =0.0

  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
  qscalar=0.0
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
  vegtyp  =0.0
  veg     =0.0

  tsoil   =0.0
  qsoil   =0.0
  wetcanp =0.0
  snowdpth=0.0

  dum2da = 0.0

!write(0,*) 'second reallocating'
  IF (ALLOCATED(xs)) DEALLOCATE(xs,ys,zps,temp)

  ALLOCATE(xs(nx),         STAT = istatus)
  ALLOCATE(ys(ny),         STAT = istatus)
  ALLOCATE(zps(nx,ny,nz),  STAT = istatus)
  ALLOCATE(temp(nx,ny,nz), STAT = istatus)

  xs   = 0.0
  ys   = 0.0
  zps  = 0.0
  temp = 0.0

!---------------------------------------------------------------------------
!
! Allocate WRF arrays
!
!---------------------------------------------------------------------------

  IF (ALLOCATED(zlevels_half)) THEN
    DEALLOCATE(zlevels_half)
    DEALLOCATE(mu,hterain_wrf)
    DEALLOCATE(soiltyp_wrf,vegtyp_wrf)
    DEALLOCATE(vegfrct_wrf,xice,xland)
    DEALLOCATE(sst_wrf,hgt_wrf,tmn_wrf,shdmin,shdmax,albbck,snoalb)
    DEALLOCATE(snowh_wrf,canwat_wrf)
  END IF

!write(0,*) 'allocating wrf arrys'
  ALLOCATE(zlevels_half(nz_wrf-1),     STAT= istatus)
  DO k = 1, nz_wrf-1
    zlevels_half(k) = (zlevels_wrf(k)+zlevels_wrf(k+1))*0.5
  END DO

  ALLOCATE(mu         (nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(hterain_wrf(nx_wrf,ny_wrf), STAT = istatus)

  ALLOCATE(soiltyp_wrf(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(vegtyp_wrf (nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(vegfrct_wrf(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(xice       (nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(xland      (nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(sst_wrf(nx_wrf,  ny_wrf  ), STAT = istatus)
  ALLOCATE(hgt_wrf(nx_wrf,  ny_wrf  ), STAT = istatus)
  ALLOCATE(tmn_wrf(nx_wrf,  ny_wrf  ), STAT = istatus)
  ALLOCATE(shdmin (nx_wrf,  ny_wrf  ), STAT = istatus)
  ALLOCATE(shdmax (nx_wrf,  ny_wrf  ), STAT = istatus)
  ALLOCATE(albbck (nx_wrf,  ny_wrf  ), STAT = istatus)
  ALLOCATE(snoalb (nx_wrf,  ny_wrf  ), STAT = istatus)

  ALLOCATE(snowh_wrf (nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(canwat_wrf(nx_wrf,ny_wrf), STAT = istatus)

  IF (ALLOCATED(u_wrf)) THEN
    DEALLOCATE(u_wrf,v_wrf,w_wrf,ph_wrf,phb_wrf,pot_wrf)
    DEALLOCATE(pt_wrf,p_wrf,pb_wrf,qv_wrf,pt_init_wrf)
    DEALLOCATE(qscalar_wrf,mup_wrf,mub_wrf)
    DEALLOCATE(tem1_wrf,tem2_wrf,tem3_wrf,tem4_wrf)
    DEALLOCATE(zs_wrf,dzs_wrf,msft_wrf,msfu_wrf,msfv_wrf)
  END IF

  ALLOCATE(u_wrf  (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(v_wrf  (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(w_wrf  (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(ph_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(phb_wrf(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(pot_wrf(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(pt_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(pt_init_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(p_wrf  (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(pb_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(qv_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(qscalar_wrf (nx_wrf,ny_wrf,nz_wrf,nscalar_wrf), STAT = istatus)
  ALLOCATE(mup_wrf(nx_wrf,ny_wrf),        STAT = istatus)
  ALLOCATE(mub_wrf(nx_wrf,ny_wrf),        STAT = istatus)

  ALLOCATE(tem1_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(tem2_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(tem3_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
  ALLOCATE(tem4_wrf (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)

  ALLOCATE(zs_wrf (9),    STAT = istatus)  ! maximum soil layer is 9
  ALLOCATE(dzs_wrf(9),    STAT = istatus)

  ALLOCATE(msft_wrf(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(msfu_wrf(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(msfv_wrf(nx_wrf,ny_wrf), STAT = istatus)

  IF (ALLOCATED(lat_wrf)) DEALLOCATE(lat_wrf,lon_wrf)
  ALLOCATE(lat_wrf(nx_wrf,ny_wrf,3), STAT = istatus)
  ALLOCATE(lon_wrf(nx_wrf,ny_wrf,3), STAT = istatus)

  mu          = 0.0
  hterain_wrf = 0.0
  soiltyp_wrf = 0
  vegtyp_wrf  = 0
  xice        = 0.0
  xland       = 0.0
  u_wrf       = 0.0
  v_wrf       = 0.0
  ph_wrf      = 0.0
  phb_wrf     = 0.0
  pot_wrf     = 0.0
  pt_wrf      = 0.0
  pt_init_wrf = 0.0
  p_wrf       = 0.0
  pb_wrf      = 0.0
  qv_wrf      = 0.0
  qscalar_wrf = 0.0
  mup_wrf     = 0.0
  mub_wrf     = 0.0
  zs_wrf      = 0.0
  dzs_wrf     = 0.0

  snowh_wrf   = 0.0
  canwat_wrf  = 0.0

!-----------------------------------------------------------------------
!
! Allocate temporary arrays
!
!-----------------------------------------------------------------------

  nxlg_wrf = (nx_wrf-fzone_wrf)*nproc_x+fzone_wrf
  nylg_wrf = (ny_wrf-fzone_wrf)*nproc_y+fzone_wrf

  IF (ALLOCATED(tem1)) THEN
    DEALLOCATE(tem1,tem2,tem3)
    DEALLOCATE(t_tmp,p_tmp,qv_tmp,zps_tmp)
    DEALLOCATE(etap,work2d,work3d)
  END IF

!write(0,*) 'allocating temporary arrys'
  ALLOCATE(tem1(nx,ny,nz), STAT= istatus)
  ALLOCATE(tem2(nx,ny,nz), STAT= istatus)
  ALLOCATE(tem3(nx,ny,nz), STAT= istatus)

  ALLOCATE(t_tmp   (nx_wrf,ny_wrf,nz),   STAT = istatus)
  ALLOCATE(p_tmp   (nx_wrf,ny_wrf,nz),   STAT = istatus)
  ALLOCATE(qv_tmp  (nx_wrf,ny_wrf,nz),   STAT = istatus)
  ALLOCATE(zps_tmp (nx_wrf,ny_wrf,nz),   STAT = istatus)

  ALLOCATE(etap    (nx_wrf,ny_wrf,nz,3), STAT = istatus)

  ALLOCATE(work3d(nx_wrf,ny_wrf,nz), STAT = istatus)
  ALLOCATE(work2d(nx_wrf,ny_wrf),    STAT = istatus)

!write(0,*) 'create_bdy = ',create_bdy
  IF(create_bdy >= 1) THEN

    ALLOCATE(ubdy3dtemp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(vbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(tbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(pbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(qbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(mbdy2dtemp1(nx_wrf,ny_wrf),        STAT = istatus)

    ALLOCATE(ubdy3dtemp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(vbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(tbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(pbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(qbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
    ALLOCATE(mbdy2dtemp2(nx_wrf,ny_wrf),        STAT = istatus)


    ALLOCATE(bdys(nx_wrf,nz_wrf,   spec_bdy_width), STAT = istatus)
    ALLOCATE(bdyn(nx_wrf,nz_wrf,   spec_bdy_width), STAT = istatus)
    ALLOCATE(bdyw(ny_wrf,nz_wrf,   spec_bdy_width), STAT = istatus)
    ALLOCATE(bdye(ny_wrf,nz_wrf,   spec_bdy_width), STAT = istatus)

    ALLOCATE(blgs(nxlg_wrf, nz_wrf, spec_bdy_width), STAT = istatus)
    ALLOCATE(blgn(nxlg_wrf, nz_wrf, spec_bdy_width), STAT = istatus)
    ALLOCATE(blgw(nylg_wrf, nz_wrf, spec_bdy_width), STAT = istatus)
    ALLOCATE(blge(nylg_wrf, nz_wrf, spec_bdy_width), STAT = istatus)

    blgs = 0.0
    blgn = 0.0
    blgw = 0.0
    blge = 0.0

    bdys = 0.0
    bdyn = 0.0
    bdyw = 0.0
    bdye = 0.0

  END IF

  IF (ALLOCATED(x2d)) DEALLOCATE(x2d,y2d,iloc,jloc)
  ALLOCATE(x2d    (nx_wrf,ny_wrf,3), STAT = istatus)
  ALLOCATE(y2d    (nx_wrf,ny_wrf,3), STAT = istatus)
  ALLOCATE(iloc   (nx_wrf,ny_wrf,3), STAT = istatus)
  ALLOCATE(jloc   (nx_wrf,ny_wrf,3), STAT = istatus)

  x2d  = 0.0
  y2d  = 0.0
  iloc = 0
  jloc = 0

!
! mpi variables
!
  nxsm = (nx-fzone_arps)/ncompressx + fzone_arps
  nysm = (ny-fzone_arps)/ncompressy + fzone_arps
  nxlg = (nx-fzone_arps)*nproc_x + fzone_arps
  nylg = (ny-fzone_arps)*nproc_y + fzone_arps

  IF(ncompressx > 1 .OR. ncompressy > 1) THEN

    IF (ASSOCIATED(xsm)) THEN
      DEALLOCATE(xsm,ysm,zpsm,zpsoilsm)
      DEALLOCATE(uprtsm,vprtsm,wprtsm,ptprtsm,pprtsm,qvprtsm)
      DEALLOCATE(qscalarsm)
      DEALLOCATE(ubarsm,vbarsm,wbarsm,ptbarsm,pbarsm,rhobarsm,qvbarsm)
      DEALLOCATE(soiltypsm,stypfrctsm,vegtypsm,vegsm)
      DEALLOCATE(tsoilsm,qsoilsm,wetcanpsm,snowdpthsm)
    END IF

    ALLOCATE(xsm   (nxsm),STAT=istatus)
    ALLOCATE(ysm   (nysm),STAT=istatus)
    ALLOCATE(zpsm    (nxsm,nysm,nz),    STAT=istatus)
    ALLOCATE(zpsoilsm(nxsm,nysm,nzsoil),STAT=istatus)

    ALLOCATE(uprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(vprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(wprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(ptprtsm(nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(pprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(qvprtsm(nxsm,nysm,nz),STAT=istatus)

    ALLOCATE(qscalarsm  (nxsm,nysm,nz,nscalar),STAT=istatus)

    ALLOCATE(ubarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(vbarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(wbarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(ptbarsm (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(pbarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(rhobarsm(nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(qvbarsm (nxsm,nysm,nz),STAT=istatus)

    ALLOCATE(soiltypsm (nxsm,nysm,nstyps),STAT=istatus)
    ALLOCATE(stypfrctsm(nxsm,nysm,nstyps),STAT=istatus)
    ALLOCATE(vegtypsm  (nxsm,nysm),STAT=istatus)
    ALLOCATE(vegsm     (nxsm,nysm),STAT=istatus)

    ALLOCATE(tsoilsm   (nxsm,nysm,nzsoil,0:nstyps),STAT=istatus)
    ALLOCATE(qsoilsm   (nxsm,nysm,nzsoil,0:nstyps),STAT=istatus)
    ALLOCATE(wetcanpsm (nxsm,nysm,0:nstyps),STAT=istatus)
    ALLOCATE(snowdpthsm(nxsm,nysm),STAT=istatus)

    CALL check_alloc_status(istatus, "arpsplt:snowdpth")
    xsm     = 0.0
    ysm     = 0.0
    zpsm    = 0.0
    zpsoilsm= 0.0

    uprtsm  = 0.0
    vprtsm  = 0.0
    wprtsm  = 0.0
    ptprtsm = 0.0
    pprtsm  = 0.0
    qvprtsm = 0.0
    qscalarsm    = 0.0
    ubarsm  = 0.0
    vbarsm  = 0.0
    wbarsm  = 0.0
    ptbarsm = 0.0
    pbarsm  = 0.0
    rhobarsm= 0.0
    qvbarsm = 0.0

    soiltypsm = 0.0
    stypfrctsm= 0.0
    vegtypsm  = 0.0
    vegsm     = 0.0

    tsoilsm   = 0.0
    qsoilsm   = 0.0
    wetcanpsm = 0.0
    snowdpthsm= 0.0

  ELSE
    xsm      => x
    ysm      => y
    zpsm     => zp
    zpsoilsm => zpsoil
    uprtsm   => uprt
    vprtsm   => vprt
    wprtsm   => wprt
    ptprtsm  => ptprt
    pprtsm   => pprt
    qvprtsm  => qvprt
    qscalarsm => qscalar
    ubarsm   => ubar
    vbarsm   => vbar
    wbarsm   => wbar
    ptbarsm  => ptbar
    pbarsm   => pbar
    rhobarsm => rhobar
    qvbarsm  => qvbar

    soiltypsm  => soiltyp
    stypfrctsm => stypfrct
    vegtypsm   => vegtyp
    vegsm      => veg

    tsoilsm    => tsoil
    qsoilsm    => qsoil
    wetcanpsm  => wetcanp
    snowdpthsm => snowdpth

  END IF

!--------------------------------------------------------------------
!
! Begin loop over all of the history files
!
!--------------------------------------------------------------------

  count_bdy = 0


  DO ifile = 1,nhisfile

    lenfil = LEN_TRIM(hisfile(ifile))
    CALL strlnth( hisfile(ifile), lenfil)
    lenfil0 = lenfil

    filename = hisfile(ifile)

    IF (finexist(ifile) == 1) THEN    ! file to be read exist

      DO jj = 1, ncompressy
      DO ii = 1, ncompressx
        IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
            CALL gtsplitfn(hisfile(ifile),ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                           0,0,1,lvldbg,filename,istatus)
            lenfil = LEN_TRIM(filename)

          WRITE(6,'(a,I5,2a/)') 'Processor: ',myproc,                   &
                            ' reading file: ', filename(1:lenfil)
        ELSE IF (myproc == 0) THEN
            WRITE(6,'(/2a/)') ' Reading file: ',filename(1:lenfil)
        END IF

!
!-----------------------------------------------------------------------
!
!  Set the gridread parameter to 0 so that the boundary
!  grid/base file will be read.
!
!-----------------------------------------------------------------------
!
        IF( mgrdbas == 0 .OR. ifile == 1) THEN
          grdbasfn = bdybasfn(1)
          IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
            CALL gtsplitfn(bdybasfn(1),ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                           0,0,1,lvldbg,grdbasfn,istatus)
          END IF
          lengbf   = LEN_TRIM(grdbasfn)

          IF (myproc == 0) WRITE(6,'(1x,a,a)') 'The grid/base name is ', grdbasfn(1:lengbf)
          IF (ifile == 1) CALL setgbrd(0)
        ELSE IF (mgrdbas == 1 .AND. ifile == 2) THEN
          grdbasfn = bdybasfn(2)
          IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
            CALL gtsplitfn(bdybasfn(2),ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                           0,0,1,lvldbg,grdbasfn,istatus)
          END IF
          lengbf   = LEN_TRIM(grdbasfn)

          IF (myproc == 0) WRITE(6,'(1x,a,a)') ' The grid/base name is ', grdbasfn(1:lengbf)
          CALL setgbrd(0)
        ELSE IF (mgrdbas == 2 .AND. ifile >= 2) THEN
          grdbasfn = bdybasfn(ifile)
          IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
            CALL gtsplitfn(bdybasfn(ifile),ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                           0,0,1,lvldbg,grdbasfn,istatus)
          END IF
          lengbf   = LEN_TRIM(grdbasfn)
          IF (myproc == 0) WRITE(6,'(1x,a,a)') '  The grid/base name is ', grdbasfn(1:lengbf)
          CALL setgbrd(0)
        END IF
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    CALL dtaread(nxsm,nysm,nz,nzsoil,nstyps,hinfmt(domid),nchr,         &
           grdbasfn(1:lengbf),lengbf, filename(1:lenfil),lenfil,        &
           time,                                                        &
           xsm,ysm,z,zpsm,zpsoilsm,uprtsm,vprtsm,wprtsm,ptprtsm,pprtsm, &
           qvprtsm, qscalarsm, dum3da,dum3da,dum3da,                    &
           ubarsm, vbarsm, wbarsm, ptbarsm, pbarsm, rhobarsm, qvbarsm,  &
           soiltypsm,stypfrctsm,vegtypsm,dum2da,dum2da,vegsm,           &
           tsoilsm,qsoilsm,wetcanpsm,snowdpthsm,                        &
           dum2da,dum2da,dum3da,                                        &
           dum3da,dum2da,dum2da,dum2da,dum2da,                          &
           dum2da,dum2da,dum2da,dum2da,                                 &
           ireturn, tem1,tem2,tem3)

!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
        IF( ireturn /= 0 )  CALL arpsstop('dtaread errors.',1)

        IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! need join

          DO j = 1, nysm
            ja = (jj-1)*(nysm-3)+j
            DO i = 1, nxsm
              ia = (ii-1)*(nxsm-3)+i
              x(ia) = xsm(i)
              vegtyp(ia,ja) = vegtypsm(i,j)
              veg(ia,ja)    = vegsm(i,j)
              snowdpth(ia,ja) = snowdpthsm(i,j)

              zp(ia,ja,:)      = zpsm(i,j,:)
              zpsoil(ia,ja,:)  = zpsoilsm(i,j,:)
              uprt(ia,ja,:)  = uprtsm(i,j,:)
              vprt(ia,ja,:)  = vprtsm(i,j,:)
              wprt(ia,ja,:)  = wprtsm(i,j,:)
              ptprt(ia,ja,:) = ptprtsm(i,j,:)
              pprt(ia,ja,:)  = pprtsm(i,j,:)
              qvprt(ia,ja,:) = qvprtsm(i,j,:)
              qscalar(ia,ja,:,:)  = qscalarsm(i,j,:,:)
              ubar(ia,ja,:)   = ubarsm(i,j,:)
              vbar(ia,ja,:)   = vbarsm(i,j,:)
              wbar(ia,ja,:)   = wbarsm(i,j,:)
              ptbar(ia,ja,:)  = ptbarsm(i,j,:)
              pbar(ia,ja,:)   = pbarsm(i,j,:)
              rhobar(ia,ja,:) = rhobarsm(i,j,:)
              qvbar(ia,ja,:)  = qvbarsm(i,j,:)
              soiltyp(ia,ja,:)  = soiltypsm(i,j,:)
              stypfrct(ia,ja,:) = stypfrctsm(i,j,:)

              tsoil(ia,ja,:,:) = tsoilsm(i,j,:,:)
              qsoil(ia,ja,:,:) = qsoilsm(i,j,:,:)
              wetcanp(ia,ja,:) = wetcanpsm(i,j,:)
            END DO    !i
            y(ja) = ysm(j)
          END DO      !j

        END IF   ! join needed
      END DO
      END DO

    ELSE         ! file to be read does not exist
      WRITE(6,'(1x,3a)') 'WARNING: File ',hisfile(ifile)(1:lenfil),     &
                         ' is skipped because it does not exist.'
      CYCLE      ! read next file, time interpolation needed
    END IF
    ! Data read finished
!CALL a3dmax0(qsoil(1,1,1,0),1,nx,1,nx,1,ny,1,ny,&
!               1,1,1,1,amax,amin)
!IF (myproc == 0) print *,'After dtaread: qsoil_min= ', amin,', qsoil_max=',amax

    CALL ctim2abss( year,month,day,hour,minute,second, abstime )
    abstime = abstime + INT(time)
    CALL abss2ctim( abstime, year, month, day, hour, minute, second )
    CALL julday( year, month, day, jday)

    abstimep = abstimen
    abstimen = abstime

!=======================================================================
!
! Only do the following for the first time level
!
!=======================================================================

    IF( ifile == 1 ) THEN

      gmthr = hour + minute/60. + second/3600.

      IF(use_arps_grid(domid) == 1) THEN
        mapproj_wrf = mapproj
        sclfct_wrf  = sclfct
        trulat1_wrf = trulat1
        trulat2_wrf = trulat2
        trulon_wrf  = trulon
        ctrlat_wrf  = ctrlat
        ctrlon_wrf  = ctrlon
        dx_wrf      = x(3)-x(2)
        dy_wrf      = y(3)-y(2)
        hinterp_needed = .FALSE.

        dx_wrfnm(domid) = dx_wrf
        dy_wrfnm(domid) = dy_wrf
      ELSE
        trulat1_wrf = lattru_wrf(1)
        trulat2_wrf = lattru_wrf(2)
        trulon_wrf  = lontru_wrf
      END IF

      IF (myproc == 0) THEN
      WRITE(6,'(/a/)')  'The ARPS grid and the WRF grid are:'
      WRITE(6,'(a38)')  '              ARPS grid      WRF grid'
      WRITE(6,'(a38)')  '              ==========    =========='
      WRITE(6,'(a,I10,4x,I10)')     '  nx       = ',nx,nx_wrf
      WRITE(6,'(a,I10,4x,I10)')     '  ny       = ',ny,ny_wrf
      WRITE(6,'(a,I10,4x,I10)')     '  mapproj  = ',mapproj,mapproj_wrf
      WRITE(6,'(a,F10.2,4x,F10.2)') '  sclfct   = ',sclfct,sclfct_wrf
      WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat1  = ',trulat1,trulat1_wrf
      WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat2  = ',trulat2,trulat2_wrf
      WRITE(6,'(a,F10.2,4x,F10.2)') '  trulon   = ',trulon,trulon_wrf
      WRITE(6,'(a,F10.2,4x,F10.2)') '  ctrlat   = ',ctrlat,ctrlat_wrf
      WRITE(6,'(a,F10.2,4x,F10.2)') '  ctrlon   = ',ctrlon,ctrlon_wrf
      WRITE(6,'(a,F10.0,4x,F10.0)') '  dx       = ',dx,dx_wrf
      WRITE(6,'(a,F10.0,4x,F10.0)') '  dy       = ',dy,dy_wrf
      END IF

      IF(mapproj == mapproj_wrf .AND.                                   &
         trulat1 == trulat1_wrf .AND. trulat2 == trulat2_wrf .AND.      &
         trulon  == trulon_wrf  .AND.                                   &
         ctrlat  == ctrlat_wrf  .AND. ctrlon  == ctrlon_wrf  .AND.      &
         nx      == nx_wrf + 2  .AND. ny      == ny_wrf + 2  .AND.      &
         ABS(dx_wrf-dx) < a_small_real_number   .AND. ABS(dy_wrf-dy) < a_small_real_number ) THEN
        ! ARPS and WRF are at the same grid
        ! Horizontal interpolations are not necessary
        hinterp_needed = .FALSE.
        IF(myproc == 0)   &
        WRITE(6,'(/a/)') '  NOTE: No horizontal interpolation will be performed.'
      ELSE
        hinterp_needed = .TRUE.
        WRITE(6,'(/a/)') '  NOTE: Horizontal interpolation will be performed.'
      END IF

      IF (hinterp_needed .AND. mp_opt > 0) THEN
        WRITE(6,'(1x,/2a/)') 'ARPS2WRF_mpi does not support horizontal ', &
                  'interpolation still. Please check and try again.'
        CALL arpsstop('mpi mode still does not implemented.',1)
      END IF

      !
      ! It was put here becasue dx_wrf/dy_wrf is finalized only until now.
      !
      CALL get_check_grid_ratio(domid,nx_wrf,ny_wrf,dx_wrf,dy_wrf,          &
           dx_wrfnm(parent_id(domid)),dy_wrfnm(parent_id(domid)),           &
           a_small_real_number,grid_ratio,istatus)

      IF (istatus /= 0) CALL arpsstop('Unacceptable parent_grid_ratio or nesting grid size.',1)

      parent_grid_ratio(domid) = grid_ratio

      year1(domid)  = year        ! Save the initial fields time
      month1(domid) = month
      day1(domid)   = day
      hour1(domid)  = hour
      minute1(domid)= minute
      second1(domid)= second

      year2   = year              ! Just initialize the last boundary time
      month2  = month
      day2    = day
      hour2   = hour
      minute2 = minute
      second2 = second

      abstime1 = abstime
      abstimec = abstime
!
!-----------------------------------------------------------------------
!
!  Establish coordinate for ARPS scalar fields.
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
        DO j=1,ny
          DO i=1,nx
            zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
          END DO
        END DO
      END DO
      zps(:,:,nz) = 2*zps(:,:,nz-1)-zps(:,:,nz-2)
!
!-----------------------------------------------------------------------
!
!  Establish WRF map projection and find latitude and longitude of WRF grid.
!  and compute the lat/lon at domain corners
!
!  1 -- T point, 2 -- U point, 3 -- V point, 4 -- massless point
!
!-----------------------------------------------------------------------
!
      lattru_wrf(1) = trulat1_wrf
      lattru_wrf(2) = trulat2_wrf
      lontru_wrf    = trulon_wrf

      ! sclfct_wrf changes inside
      CALL build_wrf_grid(domid,parent_id(domid),use_arps_grid(domid),  &
                        dx_wrfnm(parent_id(domid)),dy_wrfnm(parent_id(domid)), &
                        xsub0(parent_id(domid)),ysub0(parent_id(domid)),&
                        i_parent_start(domid),j_parent_start(domid),    &
                        nx_wrf,ny_wrf,nxlg_wrf,nylg_wrf,dx_wrf,dy_wrf,  &
                        mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf,   &
                        ctrlat_wrf,ctrlon_wrf,swx_wrf,swy_wrf,          &
                        xsub0(domid),ysub0(domid),                      &
                        lat_wrf,lon_wrf,lat_ll,lat_ul,lat_ur,lat_lr,    &
                        lon_ll,lon_ul,lon_ur,lon_lr,istatus)

!      i_parent_end = i_parent_start(domid) + (nx_wrf-1)/grid_ratio
!      j_parent_end = j_parent_start(domid) + (ny_wrf-1)/grid_ratio
!
!      IF (i_parent_end > nx_wrfnm(parent_id(domid)) .OR.    &
!          j_parent_end > ny_wrfnm(parent_id(domid)) ) THEN
!        WRITE(6,'(1x,2(a,I2),a)') 'Subdomain ',domid,' exceed its parent domain ',  &
!                                    parent_id(domid),' in size.'
!        WRITE(6,'(1x,4(a,I8),a)') 'i_parent_start-end, j_parent_start-end = ',      &
!                                 i_parent_start(domid),' - ',i_parent_end,', ',      &
!                                 j_parent_start(domid),' - ',j_parent_end,'.'
!        WRITE(6,'(1x,a,I2,a,2(I8,a),/)') 'Parent domain ',parent_id(domid),' size: ', &
!                                nx_wrfnm(parent_id(domid)),', ',            &
!                                ny_wrfnm(parent_id(domid)),'.'
!
!        CALL arpsstop('Subdomain is too large.',1)
!      END IF

!-----------------------------------------------------------------------
!
!  Find x,y locations of WRF grid in terms of the ARPS grid.
!
!-----------------------------------------------------------------------

      IF (hinterp_needed .OR. sfcinitopt(domid) == 'ARPS') THEN

        latnot(1) = trulat1
        latnot(2) = trulat2
        scalef    = 1.0/sclfct

        IF (ALLOCATED(dxfld)) DEALLOCATE(dxfld,dyfld,rdxfld,rdyfld)

        ALLOCATE(dxfld (nx,3), STAT = istatus)
        ALLOCATE(rdxfld(nx,3), STAT = istatus)
        ALLOCATE(dyfld (ny,3), STAT = istatus)
        ALLOCATE(rdyfld(ny,3), STAT = istatus)

        CALL prepinterp(nx,ny,nxlg,nylg,nxlg_wrf,nylg_wrf,dx,dy,x,y,xs,ys,&
                  mapproj,scalef,latnot,trulon,ctrlat,ctrlon,swx,swy,   &
                  lat_wrf,lon_wrf,x2d,y2d,iloc,jloc,                    &
                  dxfld,rdxfld,dyfld,rdyfld,istatus)

      END IF

!-----------------------------------------------------------------------
!
!  Read in surface characteristic data to intialize
!
!-----------------------------------------------------------------------

      CALL get_sfcdt(sfcinitopt(domid),sfcdtfl,ncompressx,ncompressy,nxsm,nysm,&
            nx,ny,nstyps,dx,dy,mapproj,trulat1,trulat2,trulon,sclfct,  &
            ctrlat,ctrlon,soiltyp,stypfrct,vegtyp,veg,tem1,            &
            hinterp_needed,nx_wrf,ny_wrf,dx_wrf,dy_wrf,                &
            mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,jday,       &
            x2d(:,:,1),y2d(:,:,1),xs,ys,iloc(:,:,1),jloc(:,:,1),       &
            hgt_wrf,soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xice,xland,     &
            tmn_wrf,shdmin,shdmax,albbck,snoalb,istatus)

      IF (istatus /= 0)  &
         CALL arpsstop('ERROR: reading surface data. Aborted.',1)

!-----------------------------------------------------------------------
!
!  Set up soil model vertical layers
!
!-----------------------------------------------------------------------

      CALL init_soil_depth(sf_surface_physics,zs_wrf,dzs_wrf,nzsoil_wrf)
!
!-----------------------------------------------------------------------
!
!  Get Map scale factor
!
!-----------------------------------------------------------------------
!
      CALL lattomf(nx_wrf,ny_wrf,lat_wrf(:,:,1),msft_wrf)  !mass points
      CALL lattomf(nx_wrf,ny_wrf,lat_wrf(:,:,2),msfu_wrf)  !U points
      CALL lattomf(nx_wrf,ny_wrf,lat_wrf(:,:,3),msfv_wrf)  !V points

    END IF  ! ifile == 1

!=======================================================================
!
!  All files should do the followings
!
!=======================================================================

    ! Restore to WRF map projection

    CALL setmapr(mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf)
    CALL setorig( 1, swx_wrf, swy_wrf)

!-----------------------------------------------------------------------
!
!    Begin ARPS data conversions
!
! NOTE:
!    From now on, pprt, ptprt, qvprt, uprt,vprt, wprt will hold total
!    values of the variables instead of only the perturbation part.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          pprt(i,j,k) = pprt(i,j,k) + pbar(i,j,k)
          ptprt(i,j,k)= ptprt(i,j,k)+ ptbar(i,j,k)
          qvprt(i,j,k)= qvprt(i,j,k)+ qvbar(i,j,k)
          uprt(i,j,k) = uprt(i,j,k) + ubar(i,j,k)
          vprt(i,j,k) = vprt(i,j,k) + vbar(i,j,k)
          wprt(i,j,k) = wprt(i,j,k) + wbar(i,j,k)
          !qi(i,j,k)   = qi(i,j,k) + qs(i,j,k) + qh(i,j,k)
        END DO
      END DO
    END DO
!
!  Put temperature into temp
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          temp(i,j,k)=ptprt(i,j,k)*(pprt(i,j,k)/p0)**rddcp
        END DO
      END DO
    END DO

!----------------------------------------------------------------------
!
!  Check ARPS vertical domain
!
!----------------------------------------------------------------------

    DO j = 1,ny-1
      DO i = 1,nx-1
        IF(pprt(i,j,nz-1) > ptop(domid) + 1000) THEN
          WRITE(6,'(/2a)') 'It seems that WRF vertical domain is ',     &
                           'outside of the ARPS physical domain.'
          WRITE(6,'(/a,F6.0,a,/a,I3,a,I3,a,F6.0,a,a)')                  &
                     'WRF top pressure is ', ptop(domid),' Pascal.',    &
                     'ARPS top pressure at i = ',i,' j = ',j,           &
                     ' is ',pprt(i,j,nz-1), ' Pascal.'
          CALL arpsstop('Vertical domain unmatch.',1)
        END IF
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Horizontally interpolate ARPS state variables to WRF grid
!
!-----------------------------------------------------------------------
!
    IF (hinterp_needed) THEN
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   temp,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),  &
                   t_tmp,tem1,tem2,tem3)
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   pprt,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),  &
                   p_tmp,tem1,tem2,tem3)
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   zps,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),   &
                   zps_tmp,tem1,tem2,tem3)
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   qvprt,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1), &
                   qv_tmp,tem1,tem2,tem3)
    ELSE
      DO k = 1,nz
        DO j = 1,ny_wrf
          DO i = 1,nx_wrf
            t_tmp  (i,j,k) = temp (i+1,j+1,k)
            p_tmp  (i,j,k) = pprt (i+1,j+1,k)
            zps_tmp(i,j,k) = zps  (i+1,j+1,k)
            qv_tmp (i,j,k) = qvprt(i+1,j+1,k)
          END DO
        END DO
      END DO
    END IF

    IF (wrftrnopt(domid) == 0 .OR. sfcinitopt(domid)(1:4) == 'ARPS') THEN
      IF (hinterp_needed) THEN
        CALL hinterp(nx,ny,1,nx_wrf,ny_wrf,iorder,                      &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   zp(:,:,2),dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1), &
                   hterain_wrf,tem1,tem2,tem3)
      ELSE
        DO j = 1,ny_wrf
          DO i = 1,nx_wrf
            hterain_wrf(i,j) = zp(i+1,j+1,2)   ! terrain height of source data
          END DO
        END DO
      END IF

      WHERE(hterain_wrf < 0) hterain_wrf = 0.0

    ELSE IF (wrftrnopt(domid) == 1) THEN

      hterain_wrf(:,:) = hgt_wrf(:,:)

    ELSE
      WRITE(6,'(a,I2,a)') ' **** WARNING: wrong wrftrnopt = ',wrftrnopt(domid),' ****'
      CALL arpsstop('Wrong wrftrnopt setting.',1)
    END IF
!
!-----------------------------------------------------------------------
!
!  Convert specific humidity to mixing ratio  (kg/kg)
!
!  NOTE: qv_tmp will be mixing ratio instead of specific humidity
!        from now on.
!
!-----------------------------------------------------------------------
!
    DO k = 1,nz
      DO j = 1,ny_wrf
        DO i = 1,nx_wrf
          qv_tmp(i,j,k) = qv_tmp(i,j,k) / (1-qv_tmp(i,j,k))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Compute ETAP in ARPS vertical grid
!
!-----------------------------------------------------------------------
!
    CALL compute_eta_3d(nx_wrf,ny_wrf,nz,p_tmp,t_tmp,qv_tmp,zps_tmp,  &
                        hterain_wrf,ptop(domid),etap(:,:,:,1),mu)  ! T points

    DO k = 1, nz
      DO j = 1, ny_wrf
        DO i = 2, nx_wrf
          etap(i,j,k,2) = 0.5*(etap(i-1,j,k,1)+etap(i,j,k,1))
        END DO
        etap(1,j,k,2) = 2*etap(2,j,k,2) - etap(3,j,k,2)     ! U points
      END DO
    END DO

    DO k = 1, nz
      DO j = 2, ny_wrf
        DO i = 1, nx_wrf
          etap(i,j,k,3) = 0.5*(etap(i,j-1,k,1)+etap(i,j,k,1))
        END DO
      END DO
      etap(:,1,k,3) = 2*etap(:,2,k,3) - etap(:,3,k,3)      ! V points
    END DO

!-----------------------------------------------------------------------
!
!  Process U wind
!
!-----------------------------------------------------------------------

    IF (hinterp_needed) THEN
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,2),jloc(:,:,2),x,ys,x2d(:,:,2),y2d(:,:,2),  &
                   uprt,dxfld(:,2),dyfld(:,2),rdxfld(:,2),rdyfld(:,2),  &
                   work3d,tem1,tem2,tem3)
    ELSE
      DO k = 1,nz
        DO j = 1,ny_wrf
          DO i = 1,nx_wrf
            work3d(i,j,k) = uprt(i+1,j+1,k)
          END DO
        END DO
      END DO
    END IF
    CALL vinterp(nx_wrf,ny_wrf,nz,nz_wrf-1,korder,etap(:,:,:,2),        &
                 zlevels_half, work3d, u_wrf, .TRUE.)

!-----------------------------------------------------------------------
!
!  Process V wind
!
!-----------------------------------------------------------------------

    IF (hinterp_needed) THEN
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,3),jloc(:,:,3),xs,y,x2d(:,:,3),y2d(:,:,3),  &
                   vprt,dxfld(:,3),dyfld(:,3),rdxfld(:,3),rdyfld(:,3),  &
                   work3d,tem1,tem2,tem3)
    ELSE
      DO k = 1,nz
        DO j = 1,ny_wrf
          DO i = 1,nx_wrf
            work3d(i,j,k) = vprt(i+1,j+1,k)
          END DO
        END DO
      END DO
    END IF
    CALL vinterp(nx_wrf,ny_wrf,nz,nz_wrf-1,korder,etap(:,:,:,3),        &
                 zlevels_half,work3d, v_wrf, .TRUE.)

!-------------------------------------------------------------------
!
!   Rotate U/V from ARPS grid to WRF grid
!
!------------------------------------------------------------------

    IF(hinterp_needed) THEN

      IF (ALLOCATED(uatv_wrf)) DEALLOCATE(uatv_wrf,vatu_wrf)
      ALLOCATE(uatv_wrf  (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
      ALLOCATE(vatu_wrf  (nx_wrf,ny_wrf,nz_wrf), STAT = istatus)

      CALL rotate_UV(nx_wrf,ny_wrf,nz_wrf,                              &
                     mapproj,scalef,latnot,trulon,swx,swy,mapproj_wrf,  &
                     sclfct_wrf,lattru_wrf,lontru_wrf,swx_wrf,swy_wrf,  &
                     lon_wrf(:,:,2),lon_wrf(:,:,3),u_wrf,v_wrf,         &
                     uatv_wrf,vatu_wrf,                                 &
                     tem1_wrf,tem2_wrf,tem3_wrf,tem4_wrf,istatus)

    END IF

    w_wrf(:,:,:) = 0.0

!-----------------------------------------------------------------------
!
!  Process Potential temperature
!
!-----------------------------------------------------------------------

    IF (hinterp_needed) THEN
      CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                       &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   ptprt,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1), &
                   work3d,tem1,tem2,tem3)
    ELSE
      DO k = 1,nz
        DO j = 1,ny_wrf
          DO i = 1,nx_wrf
            work3d(i,j,k) = ptprt(i+1,j+1,k)
          END DO
        END DO
      END DO
    END IF
    CALL vinterp(nx_wrf,ny_wrf,nz,nz_wrf-1,korder,etap(:,:,:,1),        &
                 zlevels_half, work3d, pot_wrf,.TRUE.)

!-----------------------------------------------------------------------
!
!  Process Water vapor mixing ratio
!
!-----------------------------------------------------------------------

    CALL vinterp(nx_wrf,ny_wrf,nz,nz_wrf-1,korder,etap(:,:,:,1),        &
                 zlevels_half,qv_tmp,qv_wrf,.TRUE.)
    DO k = 1, nz_wrf
      DO j = 1, ny_wrf
        DO i = 1, nx_wrf
          qv_wrf(i,j,k) = MAX(0.0,qv_wrf(i,j,k))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Compute geopotential, air mass etc.
!
!  NOTE:
!    pot_wrf is potential temperature at WRF grid, total
!    qv_wrf  is WRF water vapor mixing ratio
!    pt_wrf  is the output of the purturbation potential temperature
!            pot_wrf - t0
!
!-----------------------------------------------------------------------

    CALL compute_ph(nx_wrf,ny_wrf,nz_wrf,hterain_wrf,zlevels_wrf,ptop(domid),  &
                    mu,pot_wrf,qv_wrf,mup_wrf,mub_wrf,ph_wrf,phb_wrf,   &
                    pt_wrf,p_wrf,pb_wrf,                                &
                    pt_init_wrf,tem1_wrf,tem2_wrf,tem3_wrf)


!=======================================================================
!
!  Input file only block
!
!=======================================================================

    IF (ifile == 1) THEN  !  compute WRF input file variables only

!-----------------------------------------------------------------------
!
!  Microphysics variables
!
!-----------------------------------------------------------------------

      !  Compute QCLOUD and QRAIN, QSNOW, QICE, QGRAUP
      !
      !  NOTE: They are not all zeros according to real.exe.
      !        We can initialize them here use those values from ARPS.
      !
      DO nq = 1,nscalar_wrf

        IF (arps_Q_PTR(nq) > 0) THEN
          IF (arps_Q_PTR(nq) == 22222) THEN  ! reused temporary array temp here
            temp = qscalar(:,:,:,P_QG) + qscalar(:,:,:,P_QH)
          ELSE
            temp = qscalar(:,:,:,arps_Q_PTR(nq))
          END IF

          IF (hinterp_needed) THEN
            CALL hinterp(nx,ny,nz,nx_wrf,ny_wrf,iorder,                   &
                     iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                     temp,                           &
                     dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),       &
                     work3d,tem1,tem2,tem3)
          ELSE
            DO k = 1,nz
              DO j = 1,ny_wrf
                DO i = 1,nx_wrf
                  work3d(i,j,k) = temp(i+1,j+1,k)
                END DO
              END DO
            END DO
          END IF
          CALL vinterp(nx_wrf,ny_wrf,nz,nz_wrf-1,korder,etap(:,:,:,1),    &
                       zlevels_half, work3d, qscalar_wrf(:,:,:,nq),.TRUE.)
          DO k = 1, nz_wrf
            DO j = 1, ny_wrf
              DO i = 1, nx_wrf
                qscalar_wrf(i,j,k,nq) = MAX(0.0,qscalar_wrf(i,j,k,nq))
              END DO
            END DO
          END DO
        ELSE
          qscalar_wrf(:,:,:,nq) = 0.0
        END IF

      END DO

!----------------------------------------------------------------------
!
!  Compute Soil related variables
!
!-----------------------------------------------------------------------

      !
      !  Interpolate soil variables
      !
      IF (ALLOCATED(tsoil_tmp)) DEALLOCATE(tsoil_tmp,qsoil_tmp,zpsoil_tmp,tslb_wrf,smois_wrf)
      ALLOCATE(tsoil_tmp (nx_wrf,ny_wrf,nzsoil), STAT = istatus)
      ALLOCATE(qsoil_tmp (nx_wrf,ny_wrf,nzsoil), STAT = istatus)
      ALLOCATE(zpsoil_tmp(nx_wrf,ny_wrf,nzsoil), STAT = istatus)

      ALLOCATE(tslb_wrf  (nx_wrf,ny_wrf,nzsoil_wrf), STAT = istatus)
      ALLOCATE(smois_wrf (nx_wrf,ny_wrf,nzsoil_wrf), STAT = istatus)

      tslb_wrf    = 0.0
      smois_wrf   = 0.0

      DO k = 1, nzsoil
        DO j = 1,ny
          DO i = 1,nx
            zpsoil(i,j,k) = zp(i,j,2) - zpsoil(i,j,k)
          END DO
        END DO
      END DO

      IF (hinterp_needed) THEN
        CALL hinterp(nx,ny,nzsoil,nx_wrf,ny_wrf,iorder,                 &
                 iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1),   &
          tsoil(:,:,:,0),dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1), &
                 tsoil_tmp,tem1,tem2,tem3)
        !WDT - FIXME - be careful about horizontal interpolation of soiL
        !temperature and moisture without matching soil types!
        !See how ext2arps handles this.  GMB
        CALL hinterp(nx,ny,nzsoil,nx_wrf,ny_wrf,iorder,                 &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
          qsoil(:,:,:,0),dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1), &
                   qsoil_tmp,tem1,tem2,tem3)
        CALL hinterp(nx,ny,nzsoil,nx_wrf,ny_wrf,iorder,                 &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   zpsoil,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),&
                   zpsoil_tmp,tem1,tem2,tem3)
      ELSE
        DO k = 1,nzsoil
          DO j = 1,ny_wrf
            DO i = 1,nx_wrf
!           tsoil_tmp(i,j,k) = tsoil(i+1,j+1,k,0)
!           qsoil_tmp(i,j,k) = qsoil(i+1,j+1,k,0)

!WDT 2004-01-10 GMB: for water, set WRF deep layer to be ARPS top layer
!so soil depth interpolation is disabled
            IF ( vegtyp(i+1,j+1) == 14 .AND. k > 1) THEN
              tsoil_tmp(i,j,k) = tsoil(i+1,j+1,1,0)
              qsoil_tmp(i,j,k) = qsoil(i+1,j+1,1,0)
            ELSE
              tsoil_tmp(i,j,k) = tsoil(i+1,j+1,k,0)
              qsoil_tmp(i,j,k) = qsoil(i+1,j+1,k,0)
            ENDIF
            zpsoil_tmp(i,j,k) = zpsoil(i+1,j+1,k)
            END DO
          END DO
        END DO
      END IF

CALL a3dmax0(qsoil_tmp(1,1,1),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,&
               1,1,1,1,amax,amin)
IF (myproc == 0) print *,'qsoil_tmp_min= ', amin,', qsoil_tmp_max=',amax
IF (myproc == 0) print *,'nzsoil,zpsoil:',nzsoil,zpsoil_tmp(1,1,:)
IF (myproc == 0) print *,'nzsoil_wrf,zs_wrf:',nzsoil_wrf,zs_wrf(1:nzsoil_wrf)

      CALL vinterp_soil(nx_wrf,ny_wrf,nzsoil,zpsoil_tmp,tsoil_tmp,        &
                        qsoil_tmp,nzsoil_wrf,zs_wrf,tslb_wrf,smois_wrf)

CALL a3dmax0(smois_wrf(1,1,1),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,&
               1,1,1,1,amax,amin)
IF (myproc == 0) print *,'smois_wrf= ', amin,', smois_wrf_max=',amax

      WHERE(smois_wrf <= 0  ) smois_wrf = 1.0   ! KFY 2008/01
      WHERE(smois_wrf < 1e-3) smois_wrf = 1e-3  ! KFY 2008/01
!      WHERE(smois_wrf < 0.0) smois_wrf = 0.0
      WHERE(smois_wrf > 1.0) smois_wrf = 1.0
CALL a3dmax0(smois_wrf(1,1,1),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,&
               1,1,1,1,amax,amin)
IF (myproc == 0) print *,'smois_wrf= ', amin,', smois_wrf_max=',amax

      IF (ANY(tslb_wrf < 170.0 .OR. tslb_wrf > 400) ) THEN
        WRITE(6,'(2(1x,a,/),6(10x,a,/),1x,a)')                          &
        '*********************** IMPORTANT **************************', &
        'WARNING: Soil temperature is not valid.',                      &
        'One possible reason is that flag "sfcout" is not turned on in ',&
        'the ARPS data file. Please check the namelist file for flag ', &
        '"sfcout" and regenerate the ARPS data file. The program can ', &
        'use average annual temperature read in from static file ',     &
        'temporarily for a replacement. But you should modify the ',    &
        'source code and use it at your own risk.',                     &
        '*********************** IMPORTANT **************************'

        CALL arpsstop('tsoil no available.',1)

!        DO k = 1,nz_soil_wrf
!          DO j = 1,ny_wrf
!            DO i = 1,nx_wrf
!              tslb_wrf(i,j,k) = tmn_wrf(i,j)
!            END DO
!          END DO
!        END DO

      END IF
      !
      !  Process Water equivalent of accumulated snow
      !
    IF (hinterp_needed) THEN
      CALL hinterp(nx,ny,1,nx_wrf,ny_wrf,iorder,                        &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   snowdpth,dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),  &
                   snowh_wrf,tem1,tem2,tem3)
    ELSE
      snowh_wrf = snowdpth(2:nx-1,2:ny-1)
    END IF

    WHERE(snowh_wrf < 0.0) snowh_wrf = 0.0

    !
    !  Canopy water amount (meter, in ARPS, kg/m**2 in WRF)
    !
    IF (hinterp_needed) THEN
      CALL hinterp(nx,ny,1,nx_wrf,ny_wrf,iorder,                        &
                   iloc(:,:,1),jloc(:,:,1),xs,ys,x2d(:,:,1),y2d(:,:,1), &
                   wetcanp(:,:,0),                                      &
                   dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1),       &
                   canwat_wrf,tem1,tem2,tem3)
    ELSE
      canwat_wrf(:,:) = wetcanp(2:nx-1,2:ny-1,0)
    END IF
    canwat_wrf(:,:) = canwat_wrf(:,:)*1000.

    !
    !  SST, ARPS does not provide, use tsoil at surface to approximate
    !
    sst_wrf(:,:) = 0.0
    DO j = 1, ny_wrf-1
      DO i = 1, nx_wrf-1
        IF(vegtyp_wrf(i,j) == ISWATER) THEN
          sst_wrf(i,j) = tsoil_tmp(i,j,1)

        !WDT 2004-01-20 GMB: set all SST's
        ! Use average of SST of neighbors +-2 points away
        !
        ! NOTE: Since WRF arrays do not provide overlay region for MPI-mode
        !       in current code. User should expect inconsistent with
        !       no-mpi mode. However, the difference should be small.
!
! No meaning to set SST over land, comment out by Y. Wang. on 02/10/2005.
!
!        ELSE
!          cnt = 0.125
!          sst_sum = 0.125*tsoil_tmp(i,j,nzsoil) ! own deep soil temp counts too
!          DO jj = j-2,j+2
!            DO ii = i-2,i+2
!              IF (jj>0 .AND. jj< ny_wrf-1 .AND. ii>0 .AND. ii<nx_wrf-1  &
!                  .AND. vegtyp_wrf(ii,jj)==ISWATER) THEN
!                dist_weight = 1./((jj-j)**2+(ii-i)**2)
!                cnt = cnt + dist_weight
!                sst_sum = sst_sum + dist_weight*tsoil_tmp(i,j,1)
!              ENDIF
!            END DO
!          END DO
!          IF (cnt > 0) sst_wrf(i,j) = sst_sum/cnt
        ENDIF
      END DO
    END DO

!-----------------------------------------------------------------------
!
! Consistence check
!
!-----------------------------------------------------------------------

!  oops1 = 0
!  oops2 = 0
!  DO j = 1,ny_wrf
!    DO i = 1,nx_wrf
!      IF ( ( (xland(i,j) > 1.5) .AND. (vegtyp_wrf(i,j) /= ISWATER .OR. soiltyp_wrf(i,j) /= 14) ) &
!      .OR. ( (xland(i,j) < 1.5) .AND. (vegtyp_wrf(i,j) == ISWATER .OR. soiltyp_wrf(i,j) == 14) ) &
!           ) THEN
!        IF ( sst(i,j) < 1. ) THEN
!           oops1=oops1+1
!           vegtyp_wrf(i,j)  = 5
!           soiltyp_wrf(i,j) = 8
!           xland(i,j) = 1
!        ELSE IF ( sst_wrf(i,j) > 1. ) THEN
!           oops2=oops2+1
!           vegtyp_wrf(i,j)  = ISWATER
!           soiltyp_wrf(i,j) = 14
!           xland(i,j) = 2
!        ELSE
!           print *,'the landmask and soil/veg cats do not match'
!           print *,'i,j=',i,j
!           print *,'xland=',xland(i,j)
!           print *,'ivgtyp=',vegtyp_wrf(i,j)
!           print *,'isltyp=',soiltyp_wrf(i,j)
!           print *,'iswater=', ISWATER
!           print *,'tslb=',tslb_wrf(i,j,:)
!           print *,'sst=',sst_wrf(i,j)
!        END IF
!      END IF
!    END DO
!  END DO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  OUTPUT of WRF input begin
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    filename = ' '
!-----------------------------------------------------------------------
!
! Open WRF input file
!
!-----------------------------------------------------------------------

    WRITE(filename,'(a,I2.2)') 'wrfinput_d',domid
    lenstr = LEN_TRIM(dirname)
    IF(lenstr > 0) filename = dirname(1:lenstr) // TRIM(filename)
    IF(myproc == 0) PRINT *, ' Output file name is ',TRIM(filename)

    CALL open_output_file(filename,'INPUT',wrfversion,io_form,          &
               nxlg_wrf,nylg_wrf,nz_wrf, nzsoil_wrf,spec_bdy_width,ncid)

!-----------------------------------------------------------------------
!
!  Initialize and write global attributes
!
!-----------------------------------------------------------------------

    WRITE(times_str,'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a)')&
          year,'-',month,'-',day,'_',hour,':',minute,':',second,'.0000'

    ! set damp_opt = 0

    CALL set_global_meta(nxlg_wrf,nylg_wrf,nz_wrf,execname,times_str,   &
                       domid,parent_id(domid),i_parent_start(domid),    &
                       j_parent_start(domid),grid_ratio,                &
                       dx_wrf,dy_wrf,dt,dyn_opt,diff_opt,km_opt,0,      &
                       khdif(domid),kvdif(domid),mp_physics(domid),     &
                       ra_lw_physics(domid), ra_sw_physics(domid),      &
                       sf_sfclay_physics(domid),sf_surface_physics(domid),&
                       bl_pbl_physics(domid),cu_physics(domid),         &
                       ctrlat_wrf,ctrlon_wrf,trulat1_wrf,trulat2_wrf,   &
                       trulon_wrf,                                      &
                       year,jday,hour,minute,second,mapproj_wrf,        &
                       dampcoef,wrfversion,global_meta)

    CALL write_global_attribute(ncid,io_form,wrfversion,global_meta,.FALSE.)

    IF (io_form == IO_NET) THEN  ! this call because we do not use WRF
                                 ! external IO package for NetCDF format
      CALL write_times_str(ncid,io_form,'Times',                        &
                     times_str(1:19),times_str(1:19),ifile-1)
    END IF

!-----------------------------------------------------------------------
!
! Data writing
!
!-----------------------------------------------------------------------

    CALL write_wrf_input(ncid,io_form,times_str(1:19),wrfversion,       &
                     sfcinitopt(domid),sfcdtfns(domid),                 &
                     nx_wrf,ny_wrf,nz_wrf,nxlg_wrf,nylg_wrf,            &
                     nzsoil_wrf,spec_bdy_width,fzone_wrf,dx_wrf,dy_wrf, &
                     zlevels_wrf, zlevels_half,ptop(domid),             &
                     mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,    &
                     ctrlat_wrf,ctrlon_wrf,lat_wrf(:,:,1),lon_wrf(:,:,1),&
                     lat_wrf(:,:,2),lon_wrf(:,:,2),lat_wrf(:,:,3),lon_wrf(:,:,3), &
                     lat_ll,lat_ul,lat_ur,lat_lr,                       &
                     lon_ll,lon_ul,lon_ur,lon_lr,                       &
                     msft_wrf,msfu_wrf,msfv_wrf,zs_wrf,dzs_wrf,         &
                     u_wrf,v_wrf,w_wrf,ph_wrf,phb_wrf,pt_wrf,           &
                     pt_init_wrf,p_wrf,pb_wrf,mup_wrf,mub_wrf,mu,       &
                     qv_wrf,qscalar_wrf,                                &
                     soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xland,          &
                     xice,tmn_wrf,shdmax,shdmin,snoalb,                 &
                     snowh_wrf,canwat_wrf,albbck,sst_wrf,               &
                     hterain_wrf,tsoil_tmp(:,:,1),tslb_wrf,smois_wrf,   &
                     tem1_wrf,tem2_wrf,tem3_wrf,tem4_wrf,istatus)

!-----------------------------------------------------------------------
!
!  Close WRF input file
!
!-----------------------------------------------------------------------

      CALL close_output_file(ncid,io_form)    ! close input file

!write(0,*) 'input file closed'
!      IF(ALLOCATED(tsoil_tmp)) DEALLOCATE(tsoil_tmp,qsoil_tmp,zpsoil_tmp)
!      IF(ALLOCATED(vegtyp_wrf))DEALLOCATE(soiltyp_wrf,vegtyp_wrf,vegfrct_wrf)
!      IF(ALLOCATED(xland))     DEALLOCATE(xland,xice)

    END IF  ! ifile == 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! WRF boundary dumping begin, block for all files
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF(create_bdy >=  1) THEN

      IF (spec_bdy_width > MIN(nx_wrf,ny_wrf)) THEN
        WRITE(6,'(a,/,a,I3,a,2(I4,a))') 'WARNING: ' // &
            'Lateral boundary zone is larger than the local domain.',   &
            '         Spec_bdy_width = ',spec_bdy_width,                &
            ' Local domain size = ',nx_wrf,'X',ny_wrf,'.'
        WRITE(6,'(a,/)') 'WRF boundary file cannot be generated.'
        CALL arpsstop('Boundary zone too large.',1)
      END IF

      IF (ifile == 1) THEN

        CALL couple(nx_wrf,ny_wrf,nz_wrf,'U',mup_wrf,mub_wrf,u_wrf,     &
                    msfu_wrf,ubdy3dtemp1,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'V',mup_wrf,mub_wrf,v_wrf,     &
                    msfv_wrf,vbdy3dtemp1,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'T',mup_wrf,mub_wrf,pt_wrf,    &
                    msft_wrf,tbdy3dtemp1,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'H',mup_wrf,mub_wrf,ph_wrf,    &
                    msft_wrf,pbdy3dtemp1,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'T',mup_wrf,mub_wrf,qv_wrf,    &
                    msft_wrf,qbdy3dtemp1,work2d)
        mbdy2dtemp1 = mup_wrf

      ELSE

        CALL couple(nx_wrf,ny_wrf,nz_wrf,'U',mup_wrf,mub_wrf,u_wrf,     &
                    msfu_wrf,ubdy3dtemp2,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'V',mup_wrf,mub_wrf,v_wrf,     &
                    msfv_wrf,vbdy3dtemp2,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'T',mup_wrf,mub_wrf,pt_wrf,    &
                    msft_wrf,tbdy3dtemp2,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'H',mup_wrf,mub_wrf,ph_wrf,    &
                    msft_wrf,pbdy3dtemp2,work2d)
        CALL couple(nx_wrf,ny_wrf,nz_wrf,'T',mup_wrf,mub_wrf,qv_wrf,    &
                    msft_wrf,qbdy3dtemp2,work2d)
        mbdy2dtemp2 = mup_wrf

        count_in_loop = 0
        DO WHILE (abstimec < abstimen)

          count_bdy     = count_bdy + 1          ! count in the output boundary file
          count_in_loop = count_in_loop + 1      ! count in this WHILE loop

          abstimecn = abstimec + tintv_bdywrf
          CALL abss2ctim(abstimecn, year, month, day, hour, minute, second )
          WRITE(nextbdytimes,'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a)')&
              year,'-',month,'-',day,'_',hour,':',minute,':',second,'.0000'

          IF (abstimecn <= abstimep .OR. abstimecn > abstimen ) THEN
            WRITE(0,'(1x,a,I4,/,8x,2a,/,8x,a,I12,/,8x,a,I12,a)')        &
            'ERROR: data time inconsistent for boundary ',count_bdy,    &
            'Expecting data time at or beyond: ',nextbdytimes,          &
            'with abstime           = ',abstimecn,                      &
            'but found data at time = ',abstime,'.'
            CALL arpsstop('Data time inconsistent.',1)
          END IF

          CALL abss2ctim(abstimec, year, month, day, hour, minute, second )
          CALL julday( year, month, day, jday)

          WRITE(times_str,'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a)')&
              year,'-',month,'-',day,'_',hour,':',minute,':',second,'.0000'

          IF ( count_bdy == 1) THEN

!            IF( ABS(abstime - abstime1 - tintv_bdyin) > a_small_real_number) THEN
!              WRITE(6,'(2a)') 'ERROR: the first boundary file must be ',    &
!                        'data valid at: the initial time + tintv_bdyin.'
!              WRITE(6,'(a,I4.4,5(a,I2.2),a,I6)') 'Expected data time:    ', &
!                   year1,'/',month1,'/',day1,'-',hour1,':',minute1,':',     &
!                   second1,' + ',tintv_bdyin
!              WRITE(6,'(a,I4.4,5(a,I2.2))') 'Found data valid time: ',      &
!                   year,'/',month,'/',day,'-',hour,':',minute,':',second
!              CALL arpsstop('Inconsistent time.',1)
!            END IF

            filename = 'wrfbdy_d01'
            lenstr = LEN_TRIM(dirname)
            IF(lenstr > 0) filename = dirname(1:lenstr) // TRIM(filename)
            IF (myproc == 0) PRINT *, ' Output file name is ',TRIM(filename)
            CALL open_output_file(filename,'BDY',wrfversion,io_form,        &
                     nxlg_wrf,nylg_wrf,nz_wrf,nzsoil_wrf,spec_bdy_width,    &
                     ncid)

            !
            !  Initialize and write global attributes
            !

            CALL set_global_meta(nxlg_wrf,nylg_wrf,nz_wrf,execname,       &
                         times_str,1,1,1,1,1,                             &
                         dx_wrf,dy_wrf,dt,dyn_opt,diff_opt,km_opt,0,      &
                         khdif(1),kvdif(1),                               &
                         mp_physics(1),ra_lw_physics(1), ra_sw_physics(1),&
                         sf_sfclay_physics(1),sf_surface_physics(1),      &
                         bl_pbl_physics(1),cu_physics(1),                 &
                         ctrlat_wrf,ctrlon_wrf,trulat1_wrf,trulat2_wrf,   &
                         trulon_wrf,                                      &
                         year,jday,hour,minute,second,mapproj_wrf,        &
                         dampcoef,wrfversion,global_meta)        ! set damp_opt = 0

            CALL write_global_attribute(ncid,io_form,wrfversion,global_meta,.TRUE.)

          END IF   ! initilize boundary file

          IF (myproc == 0) WRITE(6,'(a,I4.4,5(a,I2.2))')                &
            '=== Begin boundary dumps at: ', year,'/',month,'/',day,    &
                                        '-',hour,':',minute,':',second

          IF (io_form == IO_NET) THEN  ! this call because we do not use WRF
                                       ! external IO package for NetCDF format
            CALL write_times_str(ncid,io_form,'Times',                    &
                               times_str(1:19),times_str(1:19),count_bdy)
          END IF

          CALL write_times_str(ncid,io_form,'THISBDYTIME',                &
                         times_str(1:19),times_str(1:19),count_bdy)
          CALL write_times_str(ncid,io_form,'NEXTBDYTIME',                &
                         times_str(1:19),nextbdytimes(1:19),count_bdy)

!-----------------------------------------------------------------------
!
! Prepare boundary data in two time levels, abstimec & abstimecn
!
!-----------------------------------------------------------------------

          IF (count_in_loop == 1) THEN
!write(0,*) 'Count ',count_bdy,count_in_loop,' tempc => temp1'
            ubdy3dtempc => ubdy3dtemp1
            vbdy3dtempc => vbdy3dtemp1
            tbdy3dtempc => tbdy3dtemp1
            pbdy3dtempc => pbdy3dtemp1
            qbdy3dtempc => qbdy3dtemp1
            mbdy2dtempc => mbdy2dtemp1
          ELSE
!write(0,*) 'Count ',count_bdy,count_in_loop,' tempc => tempcn'
            ubdy3dtempc => ubdy3dtempcn
            vbdy3dtempc => vbdy3dtempcn
            tbdy3dtempc => tbdy3dtempcn
            pbdy3dtempc => pbdy3dtempcn
            qbdy3dtempc => qbdy3dtempcn
            mbdy2dtempc => mbdy2dtempcn
          END IF

          IF (abstimecn == abstimen) THEN
!write(0,*) 'Count ',count_bdy,count_in_loop,' tempcn => temp2'
            ubdy3dtempcn => ubdy3dtemp2
            vbdy3dtempcn => vbdy3dtemp2
            tbdy3dtempcn => tbdy3dtemp2
            pbdy3dtempcn => pbdy3dtemp2
            qbdy3dtempcn => qbdy3dtemp2
            mbdy2dtempcn => mbdy2dtemp2
          ELSE             ! time interpolation - linear
            abstdiff1 = abstimecn - abstimep
            abstdiff2 = abstimen  - abstimecn
            abstdiff  = abstimen  - abstimep
            IF (abstdiff / tintv_bdyin > 2) THEN
              WRITE(6,'(/,1x,a,I8,a,/,9x,a)')                               &
                'WARNING: The time gap is too large (',abstdiff,' seconds)',&
                         'Time interpolation may be not accurate enough.'
            END IF

!-----------------------------------------------------------------------
!
! This block is to make sure that we will allocate at most 2 sets of
! temporary arrays, bdy3d?temp1-5 and bdy2d?temp1 (where ? is 1/2).
!
!  case 1. Only do 1 loop, no temporary arrays should be allocated (previous behavior);
!  case 2. Do 2 loops, 1 set of temporary arrays should be allocated;
!  case 3. Do 3 or more loops, 2 sets of temporary arrays will be allocated.
!
!-----------------------------------------------------------------------

            IF (mod(count_in_loop,2) == 0) THEN
              IF (.NOT. ALLOCATED(bdy3d2temp1))   &
                  ALLOCATE(bdy3d2temp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d2temp2))   &
                  ALLOCATE(bdy3d2temp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d2temp3))   &
                  ALLOCATE(bdy3d2temp3(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d2temp4))   &
                  ALLOCATE(bdy3d2temp4(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d2temp5))   &
                  ALLOCATE(bdy3d2temp5(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy2d2temp1))   &
                  ALLOCATE(bdy2d2temp1(nx_wrf,ny_wrf),        STAT = istatus)

!write(0,*) 'Count ',count_bdy,count_in_loop,' tempcn => bdy3d2'
              ubdy3dtempcn => bdy3d2temp1
              vbdy3dtempcn => bdy3d2temp2
              tbdy3dtempcn => bdy3d2temp3
              pbdy3dtempcn => bdy3d2temp4
              qbdy3dtempcn => bdy3d2temp5
              mbdy2dtempcn => bdy2d2temp1
            ELSE
              IF (.NOT. ALLOCATED(bdy3d1temp1))   &
                  ALLOCATE(bdy3d1temp1(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d1temp2))   &
                  ALLOCATE(bdy3d1temp2(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d1temp3))   &
                  ALLOCATE(bdy3d1temp3(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d1temp4))   &
                  ALLOCATE(bdy3d1temp4(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy3d1temp5))   &
                  ALLOCATE(bdy3d1temp5(nx_wrf,ny_wrf,nz_wrf), STAT = istatus)
              IF (.NOT. ALLOCATED(bdy2d1temp1))   &
                  ALLOCATE(bdy2d1temp1(nx_wrf,ny_wrf),        STAT = istatus)

!write(0,*) 'Count ',count_bdy,count_in_loop,' tempcn => bdy3d1'
              ubdy3dtempcn => bdy3d1temp1
              vbdy3dtempcn => bdy3d1temp2
              tbdy3dtempcn => bdy3d1temp3
              pbdy3dtempcn => bdy3d1temp4
              qbdy3dtempcn => bdy3d1temp5
              mbdy2dtempcn => bdy2d1temp1
            END IF

            w2 = abstdiff1*1.0/abstdiff
            w1 = abstdiff2*1.0/abstdiff
!write(0,*) 'Count ',count_bdy,count_in_loop,' w1,w2 = ',w1,w2
            ubdy3dtempcn(:,:,:) = w1*ubdy3dtemp1(:,:,:) + w2*ubdy3dtemp2(:,:,:)
            vbdy3dtempcn(:,:,:) = w1*vbdy3dtemp1(:,:,:) + w2*vbdy3dtemp2(:,:,:)
            tbdy3dtempcn(:,:,:) = w1*tbdy3dtemp1(:,:,:) + w2*tbdy3dtemp2(:,:,:)
            pbdy3dtempcn(:,:,:) = w1*pbdy3dtemp1(:,:,:) + w2*pbdy3dtemp2(:,:,:)
            qbdy3dtempcn(:,:,:) = w1*qbdy3dtemp1(:,:,:) + w2*qbdy3dtemp2(:,:,:)
            mbdy2dtempcn(:,:)   = w1*mbdy2dtemp1(:,:)   + w2*mbdy2dtemp2(:,:)

          END IF

!-----------------------------------------------------------------------
!
! Boundary data writing
!
!-----------------------------------------------------------------------

          CALL write_wrf_bdy(ncid,io_form,times_str(1:19),wrfversion,count_bdy,      &
                             nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,         &
                             nxlg_wrf,nylg_wrf,fzone_wrf,tintv_bdywrf,    &
                             ubdy3dtempc,ubdy3dtempcn,                    &
                             vbdy3dtempc,vbdy3dtempcn,                    &
                             tbdy3dtempc,tbdy3dtempcn,                    &
                             pbdy3dtempc,pbdy3dtempcn,                    &
                             qbdy3dtempc,qbdy3dtempcn,                    &
                             mbdy2dtempc,mbdy2dtempcn,                    &
                             bdys,bdyn,bdyw,bdye,                         &
                             blgw,blge,blgs,blgn,                         &
                             tem1_wrf,tem2_wrf,istatus)

          abstimec = abstimecn

        END DO      ! while loop

        ubdy3dtemp1 = ubdy3dtemp2
        vbdy3dtemp1 = vbdy3dtemp2
        tbdy3dtemp1 = tbdy3dtemp2
        pbdy3dtemp1 = pbdy3dtemp2
        qbdy3dtemp1 = qbdy3dtemp2
        mbdy2dtemp1 = mbdy2dtemp2

        ! close boundary file, The last ARPS file should not be skipped
        IF(ifile == nhisfile) THEN
          CALL close_output_file(ncid,io_form)

!write(0,*) 'deallocating bdy arrays'
!          IF (ALLOCATED(ubdy3dtemp1)) THEN
!            DEALLOCATE(ubdy3dtemp1,vbdy3dtemp1,tbdy3dtemp1)
!            DEALLOCATE(pbdy3dtemp1,qbdy3dtemp1,mbdy2dtemp1)
!            DEALLOCATE(ubdy3dtemp2,vbdy3dtemp2,tbdy3dtemp2)
!            DEALLOCATE(pbdy3dtemp2,qbdy3dtemp1,mbdy2dtemp2)
!            DEALLOCATE(bdys,bdyn,bdyw,bdye)
!            DEALLOCATE(blgs,blgn,blgw,blge)
!          END IF

        END IF
      END IF    ! ifile > 1

    END IF    ! create_wrf_bdy

    IF(ifile == nhisfile .AND. domid == 1) THEN
      year2  = year
      month2 = month
      day2   = day
      hour2  = hour
      minute2= minute
      second2= second
    END IF

  END DO    !  ifile

  CALL io_shutdown(io_form)

!-----------------------------------------------------------------------
!
! Generate a ready file if needed
!
!-----------------------------------------------------------------------

  IF ( readyfl == 1 .AND. myproc == 0) THEN

    CALL getunit( nchout )
    IF(create_bdy >= 1) THEN
      WRITE(filename,'(a,I2.2,3a)') 'wrfinput_bdy_d',domid,'.',fmtstr(io_form),'_ready'
    ELSE
      WRITE (filename,'(a,I2.2,3a)') 'wrfinput_d',domid,'.',fmtstr(io_form),'_ready'
    END IF
    WRITE(filename,'(a)') TRIM(dirname) // TRIM(filename)

    OPEN (UNIT=nchout,FILE=trim(filename))
    WRITE (nchout,'(a,I2.2)') 'wrfinput_d',domid
    IF(create_bdy >= 1)  WRITE(nchout,'(a)') 'wrfbdy_d01'
    CLOSE (nchout)

    CALL retunit (nchout )
  END IF

  ELSE    ! Prepare for namelist.input only

    nx_wrf = nx_wrfnm(domid)
    ny_wrf = ny_wrfnm(domid)
    dx_wrf = dx_wrfnm(domid)
    dy_wrf = dy_wrfnm(domid)

    CALL get_check_grid_ratio(domid,nx_wrf,ny_wrf,dx_wrf,dy_wrf,          &
         dx_wrfnm(parent_id(domid)),dy_wrfnm(parent_id(domid)),           &
         a_small_real_number,grid_ratio,istatus)

    IF (istatus /= 0) CALL arpsstop('Unacceptable parent_grid_ratio or nesting grid size.',1)

    parent_grid_ratio(domid) = grid_ratio

    year1(domid)   = year1(parent_id(domid))
    month1(domid)  = month1(parent_id(domid))
    day1(domid)    = day1(parent_id(domid))
    hour1(domid)   = hour1(parent_id(domid))
    minute1(domid) = minute1(parent_id(domid))
    second1(domid) = second1(parent_id(domid))

  END IF   ! input_from_file

  END DO   ! domid

!-----------------------------------------------------------------
!
! Write out a namelist file for WRF to run if needed.
!
!-----------------------------------------------------------------

  IF (create_namelist == 1 .AND. myproc == 0) THEN

    WRITE(6,'(/,1x,a)') 'Writing namelist input file for the WRF model ...'

    CALL write_wrf_namelist(nchout,TRIM(dirname)//TRIM(wrfnamelist),io_form, &
                   wrfversion,                       &
                   max_dom,input_from_file,nx_wrfnm,ny_wrfnm,nz_wrf,nzsoil_wrf,&
                   nproc_x,nproc_y,fzone_wrf,dx_wrfnm,dy_wrfnm,              &
                   dt,spec_bdy_width,parent_time_step_ratio,                 &
                   i_parent_start,j_parent_start,parent_id,parent_grid_ratio,&
                   year2,month2,day2,hour2,minute2,second2, tintv_bdywrf,    &
                   year1,month1,day1,hour1,minute1,second1,                  &
                   mp_physics,ra_lw_physics,ra_sw_physics,sf_sfclay_physics, &
                   sf_surface_physics,bl_pbl_physics,cu_physics,             &
                   diff_opt,km_opt,khdif,kvdif,nprocx_wrf,nprocy_wrf,        &
                   frames_per_outfile,restart_interval,radt,cudt,            &
                   ifsnow,w_damping,io_form_history,io_form_restart,         &
                   history_interval,intval2d,dir2d,                          &
                   indir,outdir,staticdir,readyfl,moist_adv_opt,dampcoef,istatus)

  END IF

!-----------------------------------------------------------------------
!
! Deallocate arrays, not necessary. For testing purpose only
!
!-----------------------------------------------------------------------

!write(0,*) 'Deallocating'
!  IF (ASSOCIATED(x)) THEN
!    DEALLOCATE(x,y,z,zp,zpsoil)
!    DEALLOCATE(uprt,vprt,wprt,ptprt,pprt,qvprt)
!    DEALLOCATE(qv,qc,qr,qi,qs,qh)
!    DEALLOCATE(ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar)
!    DEALLOCATE(soiltyp,stypfrct,vegtyp,veg)
!    DEALLOCATE(tsoil,qsoil,wetcanp,snowdpth)
!    DEALLOCATE(dum2da, dum3da)
!  END IF
!
!  IF (ALLOCATED(xs)) DEALLOCATE(xs,ys,zps,temp)
!
!  IF (ALLOCATED(zlevels_half)) THEN
!write(0,*) 'deallocating zlevels_half'
!    DEALLOCATE(zlevels_half)
!write(0,*) 'deallocating mu'
!    DEALLOCATE(mu,hterain_wrf)
!write(0,*) 'deallocating integer'
!    DEALLOCATE(soiltyp_wrf,vegtyp_wrf)
!write(0,*) 'deallocating next'
!    DEALLOCATE(vegfrct_wrf,xice,xland)
!write(0,*) 'deallocating sst_wrf'
!    DEALLOCATE(sst_wrf,hgt_wrf,tmn_wrf,shdmin,shdmax,albbck,snoalb)
!write(0,*) 'deallocating snowh_wrf'
!    DEALLOCATE(snowh_wrf,canwat_wrf)
!  END IF
!
!  IF (ALLOCATED(u_wrf)) THEN
!    DEALLOCATE(u_wrf,v_wrf,w_wrf,ph_wrf,phb_wrf,pot_wrf)
!    DEALLOCATE(pt_wrf,p_wrf,pb_wrf,qv_wrf,qc_wrf)
!    DEALLOCATE(qr_wrf,qi_wrf,qs_wrf,qg_wrf,mup_wrf,mub_wrf)
!    DEALLOCATE(tem1_wrf,tem2_wrf,tem3_wrf,tem4_wrf)
!    DEALLOCATE(zs_wrf,dzs_wrf,msft_wrf,msfu_wrf,msfv_wrf)
!  END IF
!
!  IF (ALLOCATED(lat_wrf)) DEALLOCATE(lat_wrf,lon_wrf)
!
!  IF (ALLOCATED(tem1)) THEN
!    DEALLOCATE(tem1,tem2,tem3)
!    DEALLOCATE(t_tmp,p_tmp,qv_tmp,zps_tmp)
!    DEALLOCATE(etap,work2d,work3d)
!  END IF
!
!  IF (ALLOCATED(x2d)) DEALLOCATE(x2d,y2d,iloc,jloc)
!
!!    IF (ASSOCIATED(xsm)) THEN
!!      DEALLOCATE(xsm,ysm,zpsm,zpsoilsm)
!!      DEALLOCATE(uprtsm,vprtsm,wprtsm,ptprtsm,pprtsm,qvprtsm)
!!      DEALLOCATE(qcsm,qrsm,qism,qssm,qhsm)
!!      DEALLOCATE(ubarsm,vbarsm,wbarsm,ptbarsm,pbarsm,rhobarsm,qvbarsm)
!!      DEALLOCATE(soiltypsm,stypfrctsm,vegtypsm,vegsm)
!!      DEALLOCATE(tsoilsm,qsoilsm,wetcanpsm,snowdpthsm)
!!    END IF
!
!  IF (ALLOCATED(dxfld)) DEALLOCATE(dxfld,dyfld,rdxfld,rdyfld)
!
!  IF (ALLOCATED(uatv_wrf)) DEALLOCATE(uatv_wrf,vatu_wrf)
!
!  IF (ALLOCATED(tsoil_tmp)) DEALLOCATE(tsoil_tmp,qsoil_tmp,zpsoil_tmp,tslb_wrf,smois_wrf)
!

  IF (mp_opt > 0) THEN
    IF (myproc == 0) WRITE(6,'(/,5x,a)') '==== ARPS2WRF_MPI terminated normally ===='
    CALL mpexit(0)
  ELSE
    WRITE(6,'(/,5x,a)') '==== ARPS2WRF terminated normally ===='
    STOP
  END IF
END PROGRAM arps2wrf
