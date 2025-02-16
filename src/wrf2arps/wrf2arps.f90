!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM WRF2ARPS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
PROGRAM wrf2arps
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Converts WRF forecast files coordinates and variables
!  to those used in ARPS and writes the information in a
!  history dump file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang, Dan Weber
!  Sept. 20 2003 Based on EXT2ARPS and ARPS2WRF.
!
!  MODIFICATION HISTORY:
!  03/25/2004 Yunheng Wang
!  Added support for message passing when WRF horizontal grid and ARPS
!  horizontal grid are the same. i.e. use_arps_grid = 1.
!
!  07/25/2006 Fanyou Kong
!  Fix a bug of not writing accumulated precipitation RAINC and RAING
!  due to resetting rainout and prcout during writing base state.
!  Also, max/min of RAINC and RAING are printed out.
!
!  Yunheng Wang (03/26/2007)
!  Removed namelist parameter frames_per_outfile. It is not determined by
!  the program automatically.
!
!  05/02/2007 Fanyou Kong
!  Make change to directly use near surface fields (t2m, q2m, u10m, etc)
!  read from WRF output, rather than to re-calculate.
!
!  3/25/2011 Fanyou Kong
!  Add support to M-Y microphysics scheme (mp_physics=9)
!    - combined QGRAUPEL+QHAIL => gh
!    - extract refl_10cm from wrfout* data
!  Special processing Ferrier scheme (mp_physics=5) graupel from
!    F_RIMEF_PHY
!  4/12/2012 Fanyou Kong
!  Revised to include SE2012 desired fields
!    - extract refl_10cm for most MP schemes (mp_physics:2,4,5,6,8,9,10,14,16)
!    - remove k=12 related variables
!    - add wspd1kmmax,crefd_max,pblh,accgrpl,acchail,up_heli16_max,refdm10c_max
!  3/24/2012 Fanyou Kong
!    - add mp_physics=17 support (NSSL 2-moment)
!  4/20/2013 Fanyou Kong
!    - Modify to SE2013 setting
!
!-----------------------------------------------------------------------
!
  USE module_wrf2arps_post

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Output ARPS grid dimensions
!
!-----------------------------------------------------------------------

  INTEGER :: nx ! Number of grid points in the x-dir of output arps grid
  INTEGER :: ny ! Number of grid points in the y-dir of output arps grid
  INTEGER :: nz ! Number of grid points in the z-dir of output arps grid
  INTEGER :: nzsoil
  INTEGER :: nstyps

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
!  data set to evenly-spaced data for use_ in defining the base
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
  REAL :: psnd(lvlprof),   zsnd(lvlprof),   tsnd(lvlprof),              &
          ptsnd(lvlprof),  rhssnd(lvlprof), qvsnd(lvlprof),             &
          rhosnd(lvlprof), usnd(lvlprof),   vsnd(lvlprof),              &
          dumsnd(lvlprof), plsnd(lvlprof)
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
!    qscalar  Hydrometeor scalars
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
  REAL, ALLOCATABLE :: x(:)       ! The x-coord. of the physical and
                                  ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y(:)       ! The y-coord. of the physical and
                                  ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z(:)       ! The z-coord. of the computational grid.
                                  ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp(:,:,:)  ! The physical height coordinate defined at
                                  ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The physical height coordinate defined at
                                  ! w-point in the soil.
  REAL, ALLOCATABLE :: j1(:,:,:)  ! Coordinate transformation Jacobian -d(zp)/dx.
  REAL, ALLOCATABLE :: j2(:,:,:)  ! Coordinate transformation Jacobian -d(zp)/dy.
  REAL, ALLOCATABLE :: j3(:,:,:)  ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL, ALLOCATABLE :: j3soil(:,:,:)     ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)  ! Coordinate transformation Jacobian  d(zp)/dz.

  REAL, ALLOCATABLE :: aj3z(:,:,:)! Coordinate transformation Jacobian defined
                                  ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL, ALLOCATABLE :: hterain(:,:) ! The height of the terrain. (m)
  REAL, ALLOCATABLE :: mapfct(:,:,:)! Map factor at scalar, u, and v points

  REAL, ALLOCATABLE :: u(:,:,:)     ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v(:,:,:)     ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w(:,:,:)     ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: pprt(:,:,:)  ! Perturbation pressure from that
                                    ! of base state atmosphere (Pascal)
  REAL, ALLOCATABLE :: ptprt(:,:,:) ! Perturbation potential temperature
                                    ! from that of base state atmosphere (K)

  REAL, ALLOCATABLE :: qv(:,:,:)    ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)

  REAL, ALLOCATABLE :: tke(:,:,:)

  REAL, ALLOCATABLE :: pbar(:,:,:)  ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: ptbar(:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: qvbar(:,:,:) ! Base state water vapor specific humidity
                                    ! (kg/kg)
  REAL, ALLOCATABLE :: ubar(:,:,:)  ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar(:,:,:)  ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar(:,:,:)  ! Base state w-velocity (m/s)

  REAL, ALLOCATABLE :: rhobar(:,:,:)! Base state density (kg/m3).
  REAL, ALLOCATABLE :: rhostr(:,:,:)! Base state density rhobar times j3.
  REAL, ALLOCATABLE :: wcont(:,:,:) ! Contravariant vertical velocity (m/s)

  REAL, ALLOCATABLE :: tsoil(:,:,:,:)  ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil(:,:,:,:)  ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:)  ! Water amount on canopy
  REAL, ALLOCATABLE :: snowdpth(:,:)   ! Snow depth (m)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct(:,:,:) ! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type

  REAL,    ALLOCATABLE :: lai   (:,:)     ! Leaf Area Index
  REAL,    ALLOCATABLE :: roufns(:,:)     ! Surface roughness
  REAL,    ALLOCATABLE :: veg   (:,:)     ! Vegetation fraction

!------------------------------------------------------------------
!
! Working array on ARPS grid
!
!------------------------------------------------------------------

  REAL, ALLOCATABLE :: xscl(:)    ! The x-coordinate of scalar points.
  REAL, ALLOCATABLE :: yscl(:)    ! The y-coordinate of scalar points.
  REAL, ALLOCATABLE :: zps(:,:,:) ! The physical height of scalar points.

  REAL, ALLOCATABLE :: csndsq(:,:,:)  ! Speed of sound squared (m**2/s**2)

  REAL, ALLOCATABLE :: xs2d(:,:)  ! X coordinate of ARPS scalar point in external grid
  REAL, ALLOCATABLE :: ys2d(:,:)  ! Y coordinate of ARPS scalar point in external grid
  REAL, ALLOCATABLE :: xu2d(:,:)  ! X coordinate of ARPS U point in external grid
  REAL, ALLOCATABLE :: yu2d(:,:)  ! Y
  REAL, ALLOCATABLE :: xv2d(:,:)  ! X                    V
  REAL, ALLOCATABLE :: yv2d(:,:)  ! Y                    V

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

  REAL, ALLOCATABLE :: tem1(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem2(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem3(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem4(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem5(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem6(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem7(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem8(:,:,:)   ! Temporary work array.

!----------------------------------------------------------------------
!
! Input WRF grid dimensions
!
!----------------------------------------------------------------------

  INTEGER :: nx_ext, ny_ext, nz_ext ! dimensions of external data grid
  INTEGER :: nzsoil_ext, nstyp_ext

  INTEGER :: iproj_ext              ! external data map projection
  REAL    :: scale_ext              ! external data map scale factor
  REAL    :: trlon_ext              ! external data true longitude
  REAL    :: trlat1_ext,trlat2_ext,latnot_ext(2)
                                    ! external data true latitude(s)
  REAL    :: ctrlat_ext, ctrlon_ext

  REAL    :: x0_ext,y0_ext          ! external data origin
  REAL    :: dx_ext,dy_ext,dt_ext

  INTEGER :: sf_surface_physics, sf_sfclay_physics, mp_physics
!
!-----------------------------------------------------------------------
!
!  WRF forecast variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x_ext(:)        ! external data x-coordinate
  REAL, ALLOCATABLE :: y_ext(:)        ! external data y-coordinate
  REAL, ALLOCATABLE :: xu_ext(:)       ! external data u x-coordinate
  REAL, ALLOCATABLE :: yu_ext(:)       ! external data u y-coordinate
  REAL, ALLOCATABLE :: xv_ext(:)       ! external data v x-coordinate
  REAL, ALLOCATABLE :: yv_ext(:)       ! external data v y-coordinate

  REAL, ALLOCATABLE :: lat_ext(:,:)    ! external data latidude
  REAL, ALLOCATABLE :: lon_ext(:,:)    ! external data longitude
  REAL, ALLOCATABLE :: latu_ext(:,:)   ! external data latidude (x-stag)
  REAL, ALLOCATABLE :: lonu_ext(:,:)   ! external data longitude (x-stag)
  REAL, ALLOCATABLE :: latv_ext(:,:)   ! external data latidude (y-stag)
  REAL, ALLOCATABLE :: lonv_ext(:,:)   ! external data longitude (y-stag)

  REAL, ALLOCATABLE :: zp_ext(:,:,:)   ! external data physical height (m)
  REAL, ALLOCATABLE :: hgt_ext(:,:,:)  ! Height (m) of scalar points
  REAL, ALLOCATABLE :: zpsoil_ext(:,:,:)! Height (m) of soil layers

  REAL, ALLOCATABLE :: p_ext(:,:,:)    ! Pressure (Pascals)
  REAL, ALLOCATABLE :: pt_ext(:,:,:)   ! Potential Temperature (K)
  REAL, ALLOCATABLE :: t_ext(:,:,:)    ! Temperature (K)
  REAL, ALLOCATABLE :: u_ext(:,:,:)    ! Eastward wind component
  REAL, ALLOCATABLE :: v_ext(:,:,:)    ! Northward wind component
  REAL, ALLOCATABLE :: vatu_ext(:,:,:) ! NOTE: only used when use_wrf_grid /= 1
  REAL, ALLOCATABLE :: uatv_ext(:,:,:) !
  REAL, ALLOCATABLE :: w_ext(:,:,:)    ! Vertical wind component
  REAL, ALLOCATABLE :: qv_ext(:,:,:)   ! Specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar_ext(:,:,:,:)

  REAL, ALLOCATABLE :: tke_ext(:,:,:)

  REAL, ALLOCATABLE :: tsoil_ext  (:,:,:,:)   ! Soil temperature  (K)
  REAL, ALLOCATABLE :: qsoil_ext  (:,:,:,:)   ! Soil moisture (m3/m3)
  REAL, ALLOCATABLE :: wetcanp_ext(:,:,:)     ! Canopy water amount

  REAL, ALLOCATABLE :: snowdpth_ext(:,:)      ! Snow depth (m)
  REAL, ALLOCATABLE :: trn_ext     (:,:)      ! External terrain (m)

  INTEGER, ALLOCATABLE :: soiltyp_ext (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct_ext(:,:,:) ! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp_ext  (:,:)   ! Vegetation type
  REAL,    ALLOCATABLE :: veg_ext     (:,:)   ! Vegetation fraction

!-------------------------------------------------------------------
!
!  Working arrays
!
!-------------------------------------------------------------------

  REAL, ALLOCATABLE :: zps_tmp   (:,:,:)  ! Height (m) (on arps grid)
  REAL, ALLOCATABLE :: zp_tmp    (:,:,:)  ! external w physical height (m)
  REAL, ALLOCATABLE :: avgzps_tmp(:,:,:)  ! WRF averge zps on ARPS grid
  REAL, ALLOCATABLE :: zpsoil_tmp(:,:,:)  ! Height (m) (on arps grid)

  REAL, ALLOCATABLE :: tsoil_tmp(:,:,:)
  REAL, ALLOCATABLE :: qsoil_tmp(:,:,:)
  REAL, ALLOCATABLE :: htrn_tmp(:,:)      ! The height of the terrain. (m)
                                          ! on arps grid
  REAL, ALLOCATABLE :: var_tmp  (:,:,:)
  REAL, ALLOCATABLE :: tem1_tmp (:,:,:)


  REAL, ALLOCATABLE :: dxfld(:)         ! on WRF grid
  REAL, ALLOCATABLE :: dyfld(:)
  REAL, ALLOCATABLE :: rdxfld(:)
  REAL, ALLOCATABLE :: rdyfld(:)
  REAL, ALLOCATABLE :: dxfldu(:)
  REAL, ALLOCATABLE :: dyfldu(:)
  REAL, ALLOCATABLE :: rdxfldu(:)
  REAL, ALLOCATABLE :: rdyfldu(:)
  REAL, ALLOCATABLE :: dxfldv(:)
  REAL, ALLOCATABLE :: dyfldv(:)
  REAL, ALLOCATABLE :: rdxfldv(:)
  REAL, ALLOCATABLE :: rdyfldv(:)

  REAL, ALLOCATABLE :: tem1_ext(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem2_ext(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem3_ext(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem4_ext(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem5_ext(:,:,:)    ! Temporary work array

  REAL, ALLOCATABLE :: xa_ext(:,:)        ! WRF x coordinate on ARPS grid
  REAL, ALLOCATABLE :: ya_ext(:,:)        ! WRF y coordinate on ARPS grid

  REAL, ALLOCATABLE :: znw(:)
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
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  NAMELIST parameter (In wrf2arps.input)
!
!-----------------------------------------------------------------------
!
  LOGICAL, PARAMETER :: SPRING_EXPERIMENT = .TRUE.

  INTEGER            :: nprocx_in, nprocy_in

  CHARACTER(LEN=256) :: dir_extd, dir_extp            ! directory of external data
  INTEGER            :: io_form
  LOGICAL            :: multifile

  CHARACTER(LEN=19)  :: init_time_str,start_time_str,end_time_str
  CHARACTER(LEN=11)  :: history_interval

  INTEGER :: iorder              ! Order of polynomial used for interpolation
                                 ! = 1  Linear
                                 ! = 2  Quadratic
                                 ! = 3  Cubic
  INTEGER :: intropt             ! Option indicating to interpolate
                                 ! perturbation or total variables:
                                 ! = 1  Interpolate perturbation variables
                                 !      and add to base sounding (default);
                                 ! = 2  Interploate total variables (except
                                 !      pressure).
  INTEGER :: nsmooth             ! Number of 27-pt (3-D 3 pt)
                                 ! smoothing passes after interpolation.
                                 ! 1 or 2 recommended.
  INTEGER :: ext_lbc             ! Option to apply lateral boundary conditions
                                 ! to the winds.
                                 ! = 0  no boundary contitions applied;
                                 ! = 1  apply zero-gradient boundary contitions
                                 !      (default).
  INTEGER :: ext_vbc             ! Option to apply vertical boundary conditions
                                 ! to w.
                                 ! = 0  no boundary contitions applied;
                                 ! = 1  apply boudary contitions specified
  INTEGER :: wrfexttrnopt        ! Terrain option for output grid:
                                 ! = 0  No terrain
                                 ! = 1  interpolate terrain from original grid.
                                 ! = 2  terrain data read in from file terndta (defined later)
                                 ! = 3  use_ terrain read in from file but merged to original
                                 !      grid at the boundaries (see also extntmrg)
  INTEGER :: extntmrg            ! Number of zones to merge for wrfexttrnopt=2

  INTEGER :: use_wrf_grid        ! Flag for ARPS grid
                                 ! = 0  Set ARPS grid using parameters in namelist
                                 ! = 1  Set ARPS grid using parameters in data file
  INTEGER :: dmp_out_joined

  INTEGER :: grid_id
  INTEGER :: icape,iaccu,iascii,i2dfmt,igempak
  INTEGER :: ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis
  INTEGER :: ibeg_offset,iend_offset,jbeg_offset,jend_offset
  CHARACTER (LEN=256) :: outheader,gemoutheader

  INTEGER :: magnitude_processor

  NAMELIST /message_passing/ nproc_x, nproc_y, max_fopen,               &
                             nprocx_in, nprocy_in

  NAMELIST /wrfdfile/ dir_extd,dir_extp,init_time_str,io_form,grid_id,  &
                      start_time_str,history_interval,end_time_str,     &
                      magnitude_processor

  NAMELIST /arpsgrid/   use_wrf_grid,nx,ny,nz,nzsoil,nstyp,dx,dy,       &
                        strhopt,dzmin,zrefsfc,dlayer1,dlayer2,strhtune, &
                        zflat,dz,dzsoil,soilmodel_option,soilstrhopt,   &
                        mapproj,trulat1,trulat2,                        &
                        trulon,sclfct,ctrlat,ctrlon

  NAMELIST /intrp_opts/ iorder, intropt, nsmooth, ext_lbc, ext_vbc,     &
                        bbc, tbc, fftopt,                               &
                        wrfexttrnopt, terndta, ternfmt, extntmrg

  NAMELIST /adjust/   csopt,csfactr,csound,hydradj,wndadj,obropt,obrzero

  NAMELIST /comment_lines/ runname,nocmnt,cmnt

  NAMELIST /output/     dirname,dmp_out_joined,exbcdmp,soildmp,terndmp, &
                        hdmpfmt,hdfcompr,                               &
                        qcexout,qrexout,qiexout,qsexout,qhexout,        &
                        basout,grdout,varout,mstout,iceout,tkeout,      &
                        trbout,rainout,sfcout,landout,prcout,           &
                        radout,flxout,filcmprs,readyfl

  NAMELIST /output_2d/outheader,gemoutheader,icape,iaccu,iascii,i2dfmt, &
          igempak,ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis,&
          ibeg_offset,iend_offset,jbeg_offset,jend_offset

  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
!-----------------------------------------------------------------------
!
!  Non-dimensional smoothing coefficient
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: smfct1 = 0.5
  REAL, PARAMETER :: smfct2 = -0.5
  REAL, PARAMETER :: rhmax  = 1.0
!
!-----------------------------------------------------------------------
!
!  Latitude and longitude for some diagnostic printing,
!  e.g. to compare to an observed sounding
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: latdiag = 34.5606
  REAL, PARAMETER :: londiag = -103.0820
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: basdmpfn
  CHARACTER (LEN=256) :: soiloutfl,temchar,ternfn
  CHARACTER (LEN=80) :: timsnd
  INTEGER            :: lfn, tmstrln, lternfn
  INTEGER            :: iprtopt,lbasdmpf,onvf,grdbas
  INTEGER            :: lenbin,lengem

  INTEGER :: wetcout,snowdout
  INTEGER :: tsoilout,qsoilout,zpsoilout

  INTEGER :: isnow,jsnow,ii,jj
  INTEGER :: is, idummy
  INTEGER :: num_scalar

  REAL    :: amin,amax
  REAL    :: qvmin,qvmax,qvval
  REAL    :: csconst,pconst,tv_ext
  REAL    :: pres,temp,qvsat,rh,tvbar,qvprt,qtot

  INTEGER, PARAMETER :: MAXWRFFIL = 100
  CHARACTER(LEN=256) :: extdname(MAXWRFFIL)
  CHARACTER(LEN=256) :: tmpstr

  INTEGER, ALLOCATABLE :: fHndl(:,:)
  CHARACTER(LEN=19)    :: timestr

  REAL    :: latnot(2)
  REAL    :: deltaz

  INTEGER :: i,j,k,nq,ksmth
  INTEGER :: strlen
  INTEGER :: istatus
  INTEGER :: nextdfil
  INTEGER :: ifile, itime
  INTEGER :: iniotfu

  INTEGER :: iextmn,iextmx,jextmn,jextmx

  INTEGER :: iyr,imo,iday,ihr,imin,isec,jldy
  INTEGER :: myr,initsec,iabssec,kftime
  LOGICAL :: rewindyr

  REAL    :: dmin,dd

  LOGICAL :: fexist
  LOGICAL :: first_time

  CHARACTER(LEN=1)   :: ach
  INTEGER            :: unum

  INTEGER            :: idist
  DOUBLE PRECISION   :: ntmergeinv, mfac

  INTEGER :: abstimes,abstimee,abstimei

  INTEGER :: frames_per_outfile(MAXWRFFIL)
  LOGICAL :: exit_early

  INTEGER :: mstout0,rainout0,prcout0,trbout0

!
!-----------------------------------------------------------------------
!
! roufns & lai convert table (see src/arpssfc/sfc_winter.tbl)
!
!-----------------------------------------------------------------------

  REAL, PARAMETER :: rfnstbl(14) = (/0.002, 0.02, 0.01, 0.03, 0.10,     &
                                     0.20,  0.40, 2.00, 0.005,0.01,     &
                                     0.02,  0.06, 0.04, 0.002 /)

  REAL, PARAMETER :: laitbl(14)  = (/0.5,   0.5,  0.02, 0.02, 0.02,     &
                                     0.05,  0.5,  0.5,  0.0,  0.02,     &
                                     0.0,   0.5,  0.5,  0.0   /)
!#
!# ARPS vegetation type -to- ARPS veg, LAI, z0 conversion table
!#
!#vtyp   veg     LAI     z0      vtyp definition
!#------------------------------------------------------------
!01      0.1     0.5     0.002   Desert
!02      0.1     0.5     0.02    Tundra
!03      0.3     0.02    0.01    Grassland
!04      0.2     0.02    0.03    Grassland with shrub cover
!05      0.3     0.02    0.10    Grassland with tree cover
!06      0.5     0.05    0.20    Deciduous forest
!07      0.80    0.5     0.40    Evergreen forest
!08      0.99    0.5     2.00    Rain forest
!09      0.01    0.0     0.005   Ice
!10      0.3     0.02    0.01    Cultivation
!11      0.0     0.0     0.02    Bog or marsh
!12      0.4     0.5     0.06    Dwarf shrub
!13      0.2     0.5     0.04    Semidesert
!14      0.0     0.0     0.002   Water
!
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat  ! compute saturation specific humidity defined in
                   ! thermolib3d.f90

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
!
!-----------------------------------------------------------------------
!
!  Initialize message passing processors.
!
!-----------------------------------------------------------------------
!
  ! Non-MPI defaults: All others initialized in mpinit_var
  mp_opt = 0
  myproc = 0
  nprocx_in  = 1
  nprocy_in  = 1
  dumpstride = 1
  readstride = 1

  CALL mpinit_proc(0)

  IF(myproc == 0) THEN
    WRITE(6,'(10(/5x,a),/)')                                            &
  '###################################################################',&
  '###################################################################',&
  '####                                                           ####',&
  '####                Welcome to WRF2ARPS                        ####',&
  '####                                                           ####',&
  '####   Program that reads in output from WRF model             ####',&
  '####          and produces ARPS history files.                 ####',&
  '####                                                           ####',&
  '###################################################################',&
  '###################################################################'

    unum = COMMAND_ARGUMENT_COUNT()
    IF (unum > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, tmpstr, strlen, istatus )
      IF ( tmpstr(1:1) == ' ' .OR. istatus /= 0 ) THEN  ! Use standard input to be backward-compatible
        unum = 5
      ELSE
        INQUIRE(FILE=TRIM(tmpstr),EXIST=fexist)
        IF (.NOT. fexist) THEN
          WRITE(6,'(1x,3a)') 'WARNING: namelist file - ',               &
                TRIM(tmpstr),' does not exist. Falling back to standard input.'
          unum = 5
        END IF
      END IF
    ELSE
      unum = 5
    END IF

    IF (unum /= 5) THEN
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(tmpstr),STATUS='OLD',FORM='FORMATTED')
      WRITE(*,'(1x,3a,/,1x,a,/)') 'Reading WRF2ARPS namelist from file - ', &
              TRIM(tmpstr),' ... ','========================================'
    ELSE
      WRITE(*,'(2(1x,a,/))') 'Waiting namelist from standard input ... ', &
                             '========================================'
    END IF

  END IF
  mgrid = 1
  nestgrd = 0

  DO k=1,lvlprof
    dumsnd(k) = 0.0
  END DO
!
!-----------------------------------------------------------------------
!
!  Read in message passing options.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ(unum,message_passing)
    WRITE(6,'(a)')'  Namelist block message_passing sucessfully read.'
  END IF
  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(max_fopen,1)
  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)
!
!-----------------------------------------------------------------------
!
!  Initialize message passing variables.
!
!-----------------------------------------------------------------------
!
  CALL mpinit_var

  ebc = 4      ! initialized for smooth3d, avgsu, avgsv etc.
  wbc = 4
  nbc = 4
  sbc = 4

  IF (loc_x /= 1)       wbc = 0
  IF (loc_x /= nproc_x) ebc = 0
  IF (loc_y /= 1)       sbc = 0
  IF (loc_y /= nproc_y) nbc = 0

!
!-----------------------------------------------------------------------
!
!  Read in namelist &wrfdfile
!
!-----------------------------------------------------------------------
!
  dir_extd = './'
  dir_extp = ''

  init_time_str         = '0000-00-00_00:00:00'
  start_time_str        = '0000-00-00_00:00:00'
  history_interval      = '00_00:00:00'
  end_time_str          = '0000-00-00_00:00:00'

  io_form               = 7
  grid_id               = 1
  magnitude_processor   = 4

  multifile             = .FALSE.          ! not in namelist
  IF (myproc == 0) THEN
    READ(unum,wrfdfile)
    WRITE(6,'(a)') '  Namelist wrfdfile read in successfully.'

    strlen = LEN_TRIM(dir_extd)
    IF(strlen > 0) THEN
      IF(dir_extd(strlen:strlen) /= '/' .AND. dir_extd(strlen:strlen) /= '>') THEN
        dir_extd(strlen+1:strlen+1) = '/'
      END IF
    ELSE
      dir_extd = './'
    END IF

    IF (LEN_TRIM(dir_extp) == 0 .OR. dir_extp(1:1) == ' ') THEN
      dir_extp = dir_extd
    ELSE IF ( dir_extp(1:5) == '<dir>') THEN
    ELSE
      strlen = LEN_TRIM(dir_extp)
      IF(strlen > 0) THEN
        IF(dir_extp(strlen:strlen) /= '/') THEN
          dir_extp(strlen+1:strlen+1) = '/'
        END IF
      ELSE
        dir_extp = './'
      END IF
    END IF

    IF (io_form > 100) THEN
      io_form = MOD(io_form,100)
      IF (mp_opt > 0) multifile = .TRUE.
    END IF

  END IF   ! myproc == 0
  CALL mpupdatec(dir_extd,256)
  CALL mpupdatec(dir_extp,256)
  CALL mpupdatei(io_form,1)
  CALL mpupdatel(multifile,1)
  CALL mpupdatec(init_time_str,19)
  CALL mpupdatec(start_time_str,19)
  CALL mpupdatec(end_time_str,19)
  CALL mpupdatec(history_interval,11)
  CALL mpupdatei(grid_id, 1)
  CALL mpupdatei(magnitude_processor,1)

!
!-----------------------------------------------------------------------
!
!  Read in namelist &arpsgrid
!
!-----------------------------------------------------------------------
!
  nx = 67
  ny = 67
  nz = 35
  nzsoil = 2
  nstyp  = 4
  dx     = 3200
  dy     = 3200
  dz     = 500
  dzsoil = 1

  strhopt  = 0
  dzmin    = 100.000
  zrefsfc  =   0.0
  dlayer1  =   0.0
  dlayer2  =   1.0e5
  strhtune =   1.0
  zflat    =   1.0e5

  soilmodel_option = 1
  soilstrhopt  = 0
  mapproj      = 2
  trulat1      = 30.0
  trulat2      = 60.0
  trulon       = -100.0
  sclfct       = 1.0
  ctrlat       = 35.0
  ctrlon       = -98.0

  IF (myproc == 0) THEN
    READ(unum,arpsgrid)
    WRITE(6,'(a)') '  Namelist arpsgrid read in successfully.'
  END IF
  CALL mpupdatei(use_wrf_grid,1)
  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdater(dx,1)
  CALL mpupdater(dy,1)
  CALL mpupdatei(mapproj,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon,1)
  CALL mpupdater(sclfct,1)
  CALL mpupdater(ctrlat,1)
  CALL mpupdater(ctrlon,1)
  CALL mpupdatei(strhopt,1)
  CALL mpupdater(dzmin,1)
  CALL mpupdater(dlayer1,1)
  CALL mpupdater(dlayer2,1)
  CALL mpupdater(strhtune,1)
  CALL mpupdater(zflat,1)
  CALL mpupdatei(nz,1)
  CALL mpupdater(dz,1)
  CALL mpupdatei(nzsoil,1)
  CALL mpupdater(dzsoil,1)
  CALL mpupdatei(nstyp,1)
  CALL mpupdatei(soilmodel_option,1)
  CALL mpupdatei(soilstrhopt,1)

  nstyps    = nstyp

  IF (mp_opt > 0 .AND. use_wrf_grid /= 1) THEN
    WRITE(6,'(a)') 'At present, wrf2arps_mpi only works when WRF &
                 & horizontal grid and ARPS horizontal grid are  &
                 & the same. Please set use_wrf_grid = 1 in      &
                 & wrf2arps.input. Or use_ no-mpi version of     &
                 & wrf2arps. The program_ stopping ...'
    CALL arpsstop('Unsupported MPI work.',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in namelist &intrp_opts
!
!-----------------------------------------------------------------------
!
  iorder   = 2
  intropt  = 1
  nsmooth  = 1
  ext_lbc  = 1
  ext_vbc  = 1
  bbc      = 1
  tbc      = 1
  fftopt   = 2
  wrfexttrnopt= 2
  extntmrg = 1
  terndta  = 'arpstern.data'
  ternfmt  = 1

  IF (myproc == 0) THEN
    READ(unum,intrp_opts)
    WRITE(6,'(a)') '  Namelist intrp_opts read in successfully.'

    IF(wrfexttrnopt == 2 .OR. wrfexttrnopt == 3) THEN
!      IF(ternfmt /= 1) THEN
!        WRITE(6,'(a)') 'WARNING: Only binary tern file is supported ',    &
!                       'at present. Reset ternfmt = 1.'
!        ternfmt = 1
!      END IF
      INQUIRE(FILE=TRIM(terndta), EXIST = fexist)
      IF(.NOT. fexist) THEN
        WRITE(6,*) 'The tern file ', TRIM(terndta),                       &
                   ' you specified does not exist.'
        CALL arpsstop( 'tern file not found',1)
      END IF
    END IF

  END IF
  CALL mpupdatei(iorder,1)
  CALL mpupdatei(intropt,1)
  CALL mpupdatei(nsmooth,1)
  CALL mpupdatei(ext_lbc,1)
  CALL mpupdatei(ext_vbc,1)
  CALL mpupdatei(bbc,1)
  CALL mpupdatei(tbc,1)
  CALL mpupdatei(fftopt,1)
  CALL mpupdatei(wrfexttrnopt,1)
  CALL mpupdatec(terndta,LEN(terndta))
  CALL mpupdatei(ternfmt,1)

!
!-----------------------------------------------------------------------
!
!  Read in namelist &adjust
!
!-----------------------------------------------------------------------
!
  csopt    = 1
  csfactr  = 0.5
  csound   = 150.0
  hydradj  = 0
  wndadj   = 2
  obropt   = 12
  obrzero  = 12000.

  IF (myproc == 0) THEN
    READ(unum,adjust)
    WRITE(6,'(a)') '  Namelist adjust read in successfully.'
  END IF
  CALL mpupdatei(csopt,1)
  CALL mpupdater(csfactr,1)
  CALL mpupdater(csound,1)
  CALL mpupdatei(hydradj,1)
  CALL mpupdatei(wndadj,1)
  CALL mpupdatei(obropt,1)
  CALL mpupdatei(obrzero,1)

!
!-----------------------------------------------------------------------
!
!  Read in namelist &comment_lines
!
!-----------------------------------------------------------------------
!
  runname = 'WRF2ARPS, Version 2.2'
  nocmnt  = 2
  cmnt(1) = ' '
  cmnt(2) = ' '

  IF (myproc == 0) THEN
    READ(unum,comment_lines)
    WRITE(6,'(a)') '  Namelist comment_lines read in successfully.'
  END IF
  CALL mpupdatec(runname,LEN(runname))
  CALL mpupdatei(nocmnt,1)
  CALL mpupdatec(cmnt(1),LEN(cmnt(1)))
  CALL mpupdatec(cmnt(2),LEN(cmnt(2)))

!-----------------------------------------------------------------------
!
!  Read in namelist &output
!
!-----------------------------------------------------------------------
!
  dirname  = './'
  exbcdmp  = 1
  soildmp  = 0
  terndmp  = 0
  hdmpfmt  = 1
  hdfcompr = 0

  qcexout = 0
  qrexout = 0
  qiexout = 0
  qsexout = 0
  qhexout = 0
  grdout   = 0

  basout   = 0
  varout   = 1
  mstout   = 1
  rainout  = 0
  prcout   = 0
  iceout   = 0
  totout   = 1
  tkeout   = 1
  trbout   = 0
  sfcout   = 0
  snowout  = 0
  landout  = 0
  radout   = 0
  flxout   = 0

  IF (myproc == 0) THEN
    READ(unum,output)
    WRITE(6,'(a)') '  Namelist output read in successfully.'

    tmpstr = dirname
    CALL get_output_dirname(0,tmpstr,0,0,dirname,istatus)
  END IF
  CALL mpupdatec(dirname,LEN(dirname))
  CALL mpupdatei(dmp_out_joined,1)
  CALL mpupdatei(exbcdmp,1)
  CALL mpupdatei(qcexout,1)
  CALL mpupdatei(qrexout,1)
  CALL mpupdatei(qiexout,1)
  CALL mpupdatei(qsexout,1)
  CALL mpupdatei(qhexout,1)
  CALL mpupdatei(soildmp,1)
  CALL mpupdatei(terndmp,1)
  CALL mpupdatei(hdmpfmt,1)
  CALL mpupdatei(hdfcompr,1)
  CALL mpupdatei(basout,1)
  CALL mpupdatei(grdout,1)
  CALL mpupdatei(varout,1)
  CALL mpupdatei(mstout,1)
  CALL mpupdatei(iceout,1)
  CALL mpupdatei(tkeout,1)
  CALL mpupdatei(trbout,1)
  CALL mpupdatei(rainout,1)
  CALL mpupdatei(sfcout,1)
  CALL mpupdatei(landout,1)
  CALL mpupdatei(prcout,1)
  CALL mpupdatei(radout,1)
  CALL mpupdatei(flxout,1)
  CALL mpupdatei(filcmprs,1)
  CALL mpupdatei(readyfl,1)

  DO i = FINDX_NUM,1,-1
    j = dmp_out_joined/10**(i-1)
    joindmp(i) = MOD(j,2)
  END DO

  CALL gtlfnkey(runname, lfnkey)
  CALL strlnth( dirname, ldirnam)

  outheader    = 'output_2d'
  gemoutheader = 'output_2d_gem'
  icape = 1
  iaccu = 1
  iascii  = 0
  i2dfmt = 0
  igempak = 0
  ilite = 0
  iltgci = 0
  icrtm = 0
  icitm = 0
  user_emis = 0
  isatid = 0
  chbgn = 0
  chend = 0
  ibeg_offset = 0
  iend_offset = 0
  jbeg_offset = 0
  jend_offset = 0

  IF(myproc == 0) THEN
    READ(unum,output_2d,END=200)
    WRITE(6,'(a)') '  Namelist output_2d read in successfully.'
    !WRITE(6,output_2d)

    strlen = INDEX(outheader,'/',.TRUE.)
    IF (strlen > 0) THEN
      tmpstr = outheader(1:strlen)
      CALL get_output_dirname(0,tmpstr,0,0,temchar,istatus)
    END IF
  GO TO 20
  200   WRITE(6,'(a)')                                                  &
         'Error reading NAMELIST block output_2d.  ',                   &
         'Default values used.'
  20 CONTINUE
  END IF
  CALL mpupdatec(outheader,256)
  CALL mpupdatec(gemoutheader,256)
  CALL mpupdatei(icape,1)
  CALL mpupdatei(iaccu,1)
  CALL mpupdatei(iascii,1)
  CALL mpupdatei(i2dfmt,1)
  CALL mpupdatei(igempak,1)
  CALL mpupdatei(ilite,1)
  CALL mpupdatei(iltgci,1)
  CALL mpupdatei(icrtm,1)
  CALL mpupdatei(icitm,1)
  CALL mpupdatei(user_emis,1)
  CALL mpupdatei(isatid,1)
  CALL mpupdatei(chbgn,1)
  CALL mpupdatei(chend,1)
  CALL mpupdatei(ibeg_offset,1)
  CALL mpupdatei(iend_offset,1)
  CALL mpupdatei(jbeg_offset,1)
  CALL mpupdatei(jend_offset,1)

  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE( unum )
    CALL retunit( unum )
  END IF

  IF (mp_opt > 0) THEN        ! should moved into mpinit_var later
    dumpstride = max_fopen
    readstride = nprocs
    IF (ANY(joindmp > 0)) dumpstride = nprocs   ! join and dump
    IF (multifile)   readstride = max_fopen
  END IF

!=======================================================================
!
! NAMELIST readings are done
!
!=======================================================================

  rewindyr = .FALSE.

  READ(start_time_str,'(I4.4,5(a,I2.2))')      &
                  year,ach,month,ach,day,ach,hour,ach,minute,ach,second
  IF (year < 1960) THEN   ! maybe ideal case
    myr  =  year
    year =  1960
    rewindyr = .TRUE.
  END IF
  CALL ctim2abss(year,month,day,hour,minute,second,abstimes)

  READ(end_time_str,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)')        &
                  year,ach,month,ach,day,ach,hour,ach,minute,ach,second
  IF (rewindyr)  year = 1960
  CALL ctim2abss(year,month,day,hour,minute,second,abstimee)

  READ(history_interval,'(I2.2,a,I2.2,a,I2.2,a,I2.2)')                  &
                                     day,ach,hour,ach,minute,ach,second
  abstimei = day*24*3600+hour*3600+minute*60+second

  IF (multifile) THEN
    IF (MOD(nprocx_in,nproc_x) /= 0 .OR. MOD(nprocy_in,nproc_y) /= 0) THEN
      WRITE(6,*) 'nprocx_in (nprocy_in) must be dividable by nproc_x (nproc_y).'
      CALL arpsstop('WRONG message passing parameter.',1)
    END IF

    ncompressx = nprocx_in/nproc_x
    ncompressy = nprocy_in/nproc_y

  ELSE     ! non-mpi or mpi with one file

    ncompressx = 1
    ncompressy = 1

  END IF

  ALLOCATE(fHndl(ncompressx,ncompressy), STAT = istatus)

!-----------------------------------------------------------------------
!
! initialize 2D output module
!
!-----------------------------------------------------------------------

  CALL init_data2d(i2dfmt,igempak,iaccu,abstimei,                       &
                   outheader,gemoutheader,icape,iascii,                 &
                   ilite,iltgci,icrtm,isatid,                           &
                   chbgn,chend,icitm,user_emis,                         &
                   ibeg_offset,iend_offset,jbeg_offset,jend_offset,     &
                   SPRING_EXPERIMENT, istatus)

!-----------------------------------------------------------------------
!
! Check the availability of files and get parameter frames_per_outfile
!
!-----------------------------------------------------------------------

  frames_per_outfile(:) = 1
  nextdfil = 0

  CALL check_wrf_files(multifile,MAXWRFFIL,grid_id,io_form,             &
       nprocx_in,ncompressx,ncompressy,magnitude_processor,             &
       abstimes,abstimei,abstimee,rewindyr,myr,                         &
       dir_extd,extdname,nextdfil,frames_per_outfile,istatus)

  IF (istatus /= 0) CALL arpsstop('ERROR in check_wrf_files, See STDOUT for details',1)

  IF (ANY(frames_per_outfile > 1) .AND. readstride < nprocs) THEN
    WRITE(6,'(3(a,/))') 'WARNING: WRF2ARPS does not support multi-frame in ',&
               '         one file for split-form WRF history files.',        &
               '         The program is stopping ... ...'
    CALL arpsstop('frames_per_outfile should be 1.',1)
  END IF

  readstride = readstride/(ncompressx*ncompressy)
  IF (readstride < 1) THEN
    WRITE(6,*) 'ERROR: readstride < 1, please check max_fopen in namelist file.'
    WRITE(6,*) '       Please remember that readstride = max_fopen/(ncompressx*ncompressy)'
    CALL arpsstop('max_fopen too small',1)
  END IF

!-----------------------------------------------------------------------
!
! Get dimension parameters from the first input file
!
!-----------------------------------------------------------------------

!  blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,readstride
    IF(myproc >= i .AND. myproc <= i+readstride-1) THEN

      CALL open_wrf_file(TRIM(extdname(1)),io_form,multifile,.TRUE.,    &
                         ncompressx,ncompressy,magnitude_processor,     &
                         fHndl)

      CALL get_wrf_metadata(fHndl,io_form,multifile,.TRUE.,1,1,         &
                         nx_ext,ny_ext,nz_ext,nzsoil_ext,               &
                         iproj_ext,trlat1_ext,trlat2_ext,trlon_ext,     &
                         ctrlat_ext,ctrlon_ext,dx_ext,dy_ext,dt_ext,    &
                         sf_surface_physics,sf_sfclay_physics,          &
                         mp_physics,num_scalar,istatus)

      CALL close_wrf_file(fHndl,io_form,multifile,.TRUE.,ncompressx,ncompressy)

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  scale_ext = 1.0
  nstyp_ext = 1

  IF(myproc == 0) WRITE(6,'(/a,4(a,I4))') ' WRF  grid dimensions:',     &
                       ' nx_ext = ',nx_ext, ', ny_ext = ',ny_ext,       &
                      ', nz_ext = ',nz_ext, ', nzsoil_ext = ', nzsoil_ext

  IF(use_wrf_grid == 1) THEN
    nx      = nx_ext + 2
    ny      = ny_ext + 2
    mapproj = iproj_ext
    sclfct  = scale_ext
    trulat1 = trlat1_ext
    trulat2 = trlat2_ext
    trulon  = trlon_ext
    ctrlat  = ctrlat_ext
    ctrlon  = ctrlon_ext
    dx      = dx_ext
    dy      = dy_ext

    IF(myproc == 0) WRITE(6,'(a,4(a,I4))') ' ARPS grid dimensions:',    &
                       ' nx     = ',nx, ', ny     = ',ny,               &
                      ', nz     = ',nz, ', nzsoil     = ', nzsoil

  END IF

  ! WRF forecast started time
  READ(init_time_str,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)')       &
                 year,ach,month,ach,day,ach,hour,ach,minute,ach,second
  CALL ctim2abss(year,month,day,hour,minute,second,initsec)

  thisdmp = abstimei
  tstart  = 0.0
  tstop   = abstimee-initsec
  latitud = ctrlat

!-----------------------------------------------------------------------
!
!  Allocate and initalize arrays based on dimension parameters
!  read in from the input file
!
!-----------------------------------------------------------------------

  ! Note that for MP version nx & ny here are global values.  They will
  ! be reassigned to their per-processor value below.

  IF (mp_opt > 0) THEN

    IF (nx /= nproc_x*((nx-3)/nproc_x)+3) THEN
      IF (myproc == 0) THEN
        WRITE (6,*) "WARNING: nx does not fit on ",nproc_x," processors."
        WRITE (6,*) "wrf2arps_mpi cannot handle this case at present.", &
                    'Use no-mpi version. Program stopping ...'
      END IF
      CALL arpsstop('nx does not fit.',1)
    END IF
    IF (ny /= nproc_y*int((ny-3)/nproc_y)+3) THEN
      IF (myproc == 0) THEN
        WRITE (6,*) "WARNING: ny does not fit on ",nproc_y," processors."
        WRITE (6,*) "wrf2arps_mpi cannot handle this case at present.", &
                    'Use no-mpi version Program stopping ...'
      END IF
      CALL arpsstop('ny does not fit.',1)
    END IF

    nx = (nx - 3)/nproc_x + 3             ! fake zone for ARPS is 3
    ny = (ny - 3)/nproc_y + 3

    nx_ext = (nx_ext - 1)/nproc_x + 1     ! fake zone for WRF is 1
    ny_ext = (ny_ext - 1)/nproc_y + 1

  ELSE
    nproc_x = 1
    nproc_y = 1
    nprocs  = 1
    joindmp(:) = 0
  END IF

  ALLOCATE(x(nx),stat=istatus)
  ALLOCATE(y(ny),stat=istatus)
  ALLOCATE(z(nz),stat=istatus)

  ALLOCATE(zp(nx,ny,nz),stat=istatus)
  ALLOCATE(zpsoil(nx,ny,nzsoil),stat=istatus)

  ALLOCATE(j1(nx,ny,nz),stat=istatus)
  ALLOCATE(j2(nx,ny,nz),stat=istatus)
  ALLOCATE(j3(nx,ny,nz),stat=istatus)
  ALLOCATE(j3soil(nx,ny,nzsoil),stat=istatus)
  ALLOCATE(j3soilinv(nx,ny,nzsoil),stat=istatus)
  ALLOCATE(aj3z(nx,ny,nz),stat=istatus)

  ALLOCATE(hterain(nx,ny),stat=istatus)
  ALLOCATE(mapfct(nx,ny,8),stat=istatus)

  ALLOCATE(u(nx,ny,nz),stat=istatus)
  ALLOCATE(v(nx,ny,nz),stat=istatus)
  ALLOCATE(w(nx,ny,nz),stat=istatus)
  ALLOCATE(pprt(nx,ny,nz),stat=istatus)
  ALLOCATE(ptprt(nx,ny,nz),stat=istatus)

  ALLOCATE(qv(nx,ny,nz),stat=istatus)
  ALLOCATE(qscalar(nx,ny,nz,nscalar), STAT = istatus)

  ALLOCATE(tke(nx,ny,nz),stat=istatus)

  CALL allocate_data2d_arps(nx,ny,istatus)

  ALLOCATE(pbar(nx,ny,nz),stat=istatus)
  ALLOCATE(ptbar(nx,ny,nz),stat=istatus)
  ALLOCATE(qvbar(nx,ny,nz),stat=istatus)
  ALLOCATE(ubar(nx,ny,nz),stat=istatus)
  ALLOCATE(vbar(nx,ny,nz),stat=istatus)
  ALLOCATE(wbar(nx,ny,nz),stat=istatus)
  ALLOCATE(rhobar(nx,ny,nz),stat=istatus)
  ALLOCATE(rhostr(nx,ny,nz),stat=istatus)
  ALLOCATE(wcont(nx,ny,nz),stat=istatus)

  ALLOCATE(wetcanp(nx,ny,0:nstyps),stat=istatus)
  ALLOCATE(snowdpth(nx,ny),stat=istatus)

  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps),stat=istatus)
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps),stat=istatus)

  ALLOCATE(soiltyp (nx,ny,nstyps),stat=istatus)
  ALLOCATE(stypfrct(nx,ny,nstyps),stat=istatus)
  ALLOCATE(vegtyp (nx,ny),stat=istatus)
  ALLOCATE(lai    (nx,ny),stat=istatus)
  ALLOCATE(roufns (nx,ny),stat=istatus)
  ALLOCATE(veg    (nx,ny),stat=istatus)

  ALLOCATE(xscl(nx),stat=istatus)
  ALLOCATE(yscl(ny),stat=istatus)
  ALLOCATE(zps(nx,ny,nz),stat=istatus)

  ALLOCATE(csndsq(nx,ny,nz),stat=istatus)

  ALLOCATE(xs2d(nx,ny),stat=istatus)
  ALLOCATE(ys2d(nx,ny),stat=istatus)
  ALLOCATE(xu2d(nx,ny),stat=istatus)
  ALLOCATE(yu2d(nx,ny),stat=istatus)
  ALLOCATE(xv2d(nx,ny),stat=istatus)
  ALLOCATE(yv2d(nx,ny),stat=istatus)

  ALLOCATE(iscl(nx,ny),stat=istatus)
  ALLOCATE(jscl(nx,ny),stat=istatus)
  ALLOCATE(iu(nx,ny),stat=istatus)
  ALLOCATE(ju(nx,ny),stat=istatus)
  ALLOCATE(iv(nx,ny),stat=istatus)
  ALLOCATE(jv(nx,ny),stat=istatus)

  ALLOCATE(tem1(nx,ny,nz),stat=istatus)
  ALLOCATE(tem2(nx,ny,nz),stat=istatus)
  ALLOCATE(tem3(nx,ny,nz),stat=istatus)
  ALLOCATE(tem4(nx,ny,nz),stat=istatus)
  ALLOCATE(tem5(nx,ny,nz),stat=istatus)
  ALLOCATE(tem6(nx,ny,nz),stat=istatus)
  ALLOCATE(tem7(nx,ny,nz),stat=istatus)
  ALLOCATE(tem8(nx,ny,nz),stat=istatus)

  x=0.0
  y=0.0
  z=0.0

  zp=0.0
  zpsoil = 0.0

  j1=0.0
  j2=0.0
  j3=0.0
  j3soil = 0.0
  j3soilinv = 0.0
  aj3z=0.0

  hterain=0.0
  mapfct=0.0

  u=0.0
  v=0.0
  w=0.0
  pprt=0.0
  ptprt=0.0
  qv=0.0
  qscalar = 0.0

  tke = 0.0

  pbar=0.0
  ptbar=0.0
  qvbar=0.0
  ubar=0.0
  vbar=0.0
  wbar=0.0
  rhobar=0.0
  rhostr=0.0
  wcont=0.0

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

  xscl=0.0
  yscl=0.0
  zps=0.0

  tem1=0.0
  tem2=0.0
  tem3=0.0
  tem4=0.0
  tem5=0.0
  tem6=0.0
  tem7=0.0
  tem8=0.0

  csndsq=0.0

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
!  Allocate and initialize external grid variables
!
!-----------------------------------------------------------------------
!
  ALLOCATE(x_ext(nx_ext),stat=istatus)
  ALLOCATE(y_ext(ny_ext),stat=istatus)
  ALLOCATE(xu_ext(nx_ext),stat=istatus)
  ALLOCATE(yu_ext(ny_ext),stat=istatus)
  ALLOCATE(xv_ext(nx_ext),stat=istatus)
  ALLOCATE(yv_ext(ny_ext),stat=istatus)

  ALLOCATE(lat_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(lon_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(latu_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(lonu_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(latv_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(lonv_ext(nx_ext,ny_ext),stat=istatus)

  ALLOCATE(zp_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(hgt_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(zpsoil_ext(nx_ext,ny_ext,nzsoil_ext),stat=istatus)

  ALLOCATE(p_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(t_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(u_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(v_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(vatu_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(uatv_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(w_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(pt_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(qv_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(qscalar_ext(nx_ext,ny_ext,nz_ext,nscalar),stat=istatus)

  ALLOCATE(tke_ext (nx_ext,ny_ext,nz_ext),stat=istatus)

  ALLOCATE(tsoil_ext   (nx_ext,ny_ext,nzsoil_ext,0:nstyp_ext),stat=istatus)
  ALLOCATE(qsoil_ext   (nx_ext,ny_ext,nzsoil_ext,0:nstyp_ext),stat=istatus)
  ALLOCATE(wetcanp_ext (nx_ext,ny_ext,0:nstyp_ext),stat=istatus)
  ALLOCATE(soiltyp_ext (nx_ext,ny_ext,nstyp_ext),stat=istatus)
  ALLOCATE(stypfrct_ext(nx_ext,ny_ext,nstyp_ext),stat=istatus)

  ALLOCATE(vegtyp_ext  (nx_ext,ny_ext),stat=istatus)
  ALLOCATE(veg_ext     (nx_ext,ny_ext),stat=istatus)
  ALLOCATE(snowdpth_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(trn_ext     (nx_ext,ny_ext),stat=istatus)

  CALL allocate_data2d_ext( nx_ext,ny_ext, istatus )

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

  zp_ext=0.0
  hgt_ext=0.0
  zpsoil_ext=0.0
  p_ext=0.0
  t_ext=0.0
  u_ext=0.0
  v_ext=0.0
  vatu_ext=0.0
  uatv_ext=0.0
  w_ext=0.0
  pt_ext=0.0
  qv_ext=0.0
  qscalar_ext=0.0

  tke_ext = -999.0

  trn_ext     =0.0

  tsoil_ext   =-999.0
  qsoil_ext   =-999.0
  wetcanp_ext =-999.0
  snowdpth_ext=-999.0

  soiltyp_ext = -999
  stypfrct_ext= -999.0
  vegtyp_ext  = -999
  veg_ext     = -999.0

!-----------------------------------------------------------------------
!
! Allocate working arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(zps_tmp(nx,ny,nz_ext),        stat=istatus)
  ALLOCATE(zp_tmp(nx,ny,nz_ext),         stat=istatus)
  ALLOCATE(zpsoil_tmp(nx,ny,nzsoil_ext), stat=istatus)
  ALLOCATE(htrn_tmp(nx,ny),              stat=istatus)
  ALLOCATE(avgzps_tmp(nx,ny,nz_ext),     stat=istatus)

  zps_tmp    = 0.0
  zp_tmp     = 0.0
  zpsoil_tmp = 0.0
  htrn_tmp   = 0.0
  avgzps_tmp = 0.0

  ALLOCATE(var_tmp(nx,ny,nz_ext),     stat=istatus)
  ALLOCATE(tem1_tmp(nx,ny,nz_ext),    stat=istatus)

  ALLOCATE(tsoil_tmp(nx,ny,nzsoil_ext),stat=istatus)
  ALLOCATE(qsoil_tmp(nx,ny,nzsoil_ext),stat=istatus)

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

  ALLOCATE(tem1_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(tem2_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(tem3_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(tem4_ext(nx_ext,ny_ext,nz_ext),stat=istatus)
  ALLOCATE(tem5_ext(nx_ext,ny_ext,nz_ext),stat=istatus)

  ALLOCATE(xa_ext(nx_ext,ny_ext),stat=istatus)
  ALLOCATE(ya_ext(nx_ext,ny_ext),stat=istatus)

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

  ALLOCATE (znw(nz_ext), STAT = istatus)

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
  IF (wrfexttrnopt == 1) THEN ! don't really do grid here, just set map proj
    ternopt = -1
    hterain = 0
  ELSE IF(wrfexttrnopt == 2 .OR. wrfexttrnopt == 3) THEN
    ternopt = 2
  ELSE
    ternopt = 0
  END IF

  CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                          &
              hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,                 &
              tem2(1,1,:),tem3(1,1,:),tem1)

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        aj3z(i,j,k)=0.5*(j3(i,j,k)+j3(i,j,k-1))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Set up the height levels on which to define the base-state.
!
!-----------------------------------------------------------------------
!
  deltaz = depthp/(lvlprof-1)

  IF(myproc ==0) WRITE(6,'(/2a,f7.2,a/)')                               &
      ' The base state is formed from a mean sounding created',         &
      ' with a ',deltaz,' meter interval.'

  DO k=1,lvlprof
    zsnd(k) = (k-1)*deltaz
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop through the data times provided via NAMELIST.
!
!-----------------------------------------------------------------------
!
  basdmpfn = ' '     ! initializes some variables
  hdmpfn   = ' '

  iniotfu = 21       ! FORTRAN unit number used for data output

  first_time = .TRUE.
  exit_early = .FALSE.

  DO ifile = 1,nextdfil

    IF (myproc == 0) WRITE(6,'(1x,2a,/)') 'Reading file ',TRIM(extdname(ifile))

    IF (frames_per_outfile(ifile) == 1) THEN  ! Finish read in the following block
      !
      !  blocking inserted for ordering i/o for message passing
      !
      DO i=0,nprocs-1,readstride
        IF(myproc >= i .AND. myproc <= i+readstride-1) THEN

          CALL open_wrf_file(TRIM(extdname(ifile)),io_form,multifile,   &
                     .FALSE.,ncompressx,ncompressy,magnitude_processor, &
                     fHndl)

!-----------------------------------------------------------------------
!
!  Retrieve global attributes from the file
!
!-----------------------------------------------------------------------

          CALL get_wrf_metadata(fHndl,io_form,multifile,.FALSE.,        &
                          ncompressx,ncompressy,                        &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj_ext,trlat1_ext,trlat2_ext,trlon_ext,    &
                          ctrlat_ext,ctrlon_ext,dx_ext,dy_ext,dt_ext,   &
                          sf_surface_physics,sf_sfclay_physics,         &
                          mp_physics,num_scalar,istatus)

          scale_ext = 1.0
          latnot_ext(1) = trlat1_ext
          latnot_ext(2) = trlat2_ext

          IF (mp_opt > 0) THEN
            nx_ext = (nx_ext - 1)/nproc_x + 1
            ny_ext = (ny_ext - 1)/nproc_y + 1
          END IF

          itime = 1

          CALL get_wrf_Times(fHndl,io_form,multifile,ncompressx,ncompressy, &
                             itime,timestr)

          CALL get_wrf_1d(fHndl,io_form,multifile,ncompressx,ncompressy,&
                      timestr,itime,'ZNW','Z',                          &
                      'bottom_top_stag',nz_ext,znw,nz_ext, istatus)
!
!-----------------------------------------------------------------------
!
!  Call getwrfd to reads and converts data to ARPS units
!
!  NOTE: u_ext, v_ext are just values extracted from data files. It may
!        need to be rotated or extend to be MPI valid.
!
!-----------------------------------------------------------------------
!
          CALL getwrfd(fHndl,io_form,multifile,ncompressx,ncompressy,itime, &
                   timestr,nx_ext,ny_ext,nz_ext,nzsoil_ext,             &
                   iproj_ext,scale_ext,trlon_ext,latnot_ext,            &
                   ctrlat_ext,ctrlon_ext,dx_ext,dy_ext,x0_ext,y0_ext,   &
                   sf_surface_physics,sf_sfclay_physics,num_scalar,     &
                   lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext, &
                   zp_ext,hgt_ext,zpsoil_ext, p_ext,t_ext,              &
                   u_ext,v_ext, w_ext,                                  &
                   qv_ext,qscalar_ext,tke_ext,                          &
                   tsoil_ext(:,:,:,0),qsoil_ext(:,:,:,0),               &
                   wetcanp_ext(:,:,0),snowdpth_ext,trn_ext,             &
                   soiltyp_ext(:,:,1),vegtyp_ext,veg_ext,               &
                   tem1_ext,tem2_ext,tem3_ext,tem4_ext,istatus)

          CALL get_wrf_data2d(fHndl,io_form,                            &
                   multifile,ncompressx,ncompressy,itime,timestr,       &
                   tem1_ext,istatus)

          CALL close_wrf_file(fHndl,io_form,multifile,.FALSE.,          &
                                ncompressx,ncompressy)

        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

    ELSE     ! frames_per_outfile > 1, only open the file and read meta

      CALL open_wrf_file(TRIM(extdname(ifile)),io_form,multifile,       &
                    .FALSE.,ncompressx,ncompressy,magnitude_processor,  &
                    fHndl)

!-----------------------------------------------------------------------
!
!  Retrieve global attributes from the file
!
!-----------------------------------------------------------------------

      CALL get_wrf_metadata(fHndl,io_form,multifile,.FALSE.,            &
                          ncompressx,ncompressy,                        &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj_ext,trlat1_ext,trlat2_ext,trlon_ext,    &
                          ctrlat_ext,ctrlon_ext,dx_ext,dy_ext,dt_ext,   &
                          sf_surface_physics,sf_sfclay_physics,         &
                          mp_physics,num_scalar,istatus)

      scale_ext = 1.0
      latnot_ext(1) = trlat1_ext
      latnot_ext(2) = trlat2_ext

      IF (mp_opt > 0) THEN
        nx_ext = (nx_ext - 1)/nproc_x + 1
        ny_ext = (ny_ext - 1)/nproc_y + 1
      END IF

      CALL get_wrf_1d(fHndl,io_form,multifile,ncompressx,ncompressy,    &
                      timestr,itime,'ZNW','Z',                          &
                      'bottom_top_stag',nz_ext,znw,nz_ext, istatus)

    END IF

    DO itime = 1, frames_per_outfile(ifile)

      IF (frames_per_outfile(ifile) > 1) THEN
        ! read data only when frames_per_outfile >1
        ! otherwise, it has been read before.
        CALL get_wrf_Times(fHndl,io_form,multifile,ncompressx,ncompressy,&
                           itime,timestr)

        IF (myproc == 0) WRITE(6,'(1x,4a,/)') 'Processing data at ',    &
                               timestr,' in file ',TRIM(extdname(ifile))
!
!-----------------------------------------------------------------------
!
!  Call getwrfd to reads and converts data to ARPS units
!
!  NOTE: u_ext, v_ext are just values extracted from data files. It may
!        need to be rotated or extend to be MPI valid.
!
!-----------------------------------------------------------------------
!
        CALL getwrfd(fHndl,io_form,multifile,ncompressx,ncompressy,itime,&
                   timestr,nx_ext,ny_ext,nz_ext,nzsoil_ext,             &
                   iproj_ext,scale_ext,trlon_ext,latnot_ext,            &
                   ctrlat_ext,ctrlon_ext,dx_ext,dy_ext,x0_ext,y0_ext,   &
                   sf_surface_physics, sf_sfclay_physics,num_scalar,    &
                   lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext, &
                   zp_ext,hgt_ext,zpsoil_ext, p_ext,t_ext,              &
                   u_ext,v_ext, w_ext,                                  &
                   qv_ext,qscalar_ext,tke_ext,                          &
                   tsoil_ext(:,:,:,0),qsoil_ext(:,:,:,0),               &
                   wetcanp_ext(:,:,0),snowdpth_ext,trn_ext,             &
                   soiltyp_ext(:,:,1),vegtyp_ext,veg_ext,               &
                   tem1_ext,tem2_ext,tem3_ext,tem4_ext,istatus)

        CALL get_wrf_data2d(fHndl,io_form,                              &
                   multifile,ncompressx,ncompressy,itime,timestr,       &
                   tem1_ext,istatus)
      END IF

!-----------------------------------------------------------------------
!
! Post-reading processing, most of this block is inside getwrfd before
! It must be moved here because of ordering I/O for message passing
!
!-----------------------------------------------------------------------

    CALL adj_wrfuv(multifile,use_wrf_grid,nx_ext,ny_ext,nz_ext,         &
                iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext, &
                lonu_ext,lonv_ext,u_ext,vatu_ext,uatv_ext,v_ext,        &
                tem1_ext,tem2_ext,istatus)

    pt_ext(:,:,:) = t_ext(:,:,:)*((p0/p_ext(:,:,:))**rddcp)

    stypfrct_ext(:,:,1) = 1.
    tsoil_ext(:,:,:,1) = tsoil_ext(:,:,:,0)
    qsoil_ext(:,:,:,1) = qsoil_ext(:,:,:,0)
    wetcanp_ext(:,:,1) = wetcanp_ext(:,:,0)

    IF(istatus /= 0) GO TO 999    ! seems meaningless here

!-----------------------------------------------------------------------
!
!  Time conversions.
!  Formats:  timestr='1998-05-25_18:00:00
!
!-----------------------------------------------------------------------
!
      READ(timestr,'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2)')    &
                 iyr,ach,imo,ach,iday,ach,ihr,ach,imin,ach,isec
      CALL julday(iyr,imo,iday,jldy)

      CALL ctim2abss(iyr,imo,iday,ihr,imin,isec,iabssec)
      kftime=iabssec - initsec
      curtim=FLOAT(kftime)

      IF(myproc == 0) WRITE(6,'(/,a,a19,/,a,i12,a,i6,a,/)')             &
          ' Returned from getwrfd with data valid at ',timestr,         &
          ' Which is ',iabssec,' abs seconds or ',kftime,               &
          ' seconds from the WRF/ARPS initial time.'

!-----------------------------------------------------------------------

      CALL print_maxmin_ext(myproc, nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext, &
                         lat_ext,lon_ext,p_ext,hgt_ext,t_ext,           &
                         u_ext,v_ext,w_ext,qv_ext,qscalar_ext,          &
                         tsoil_ext,qsoil_ext,wetcanp_ext,snowdpth_ext,  &
                         tke_ext,trn_ext,veg_ext,istatus)

      CALL print_maxmin_data2d_ext( myproc, istatus )
!
!-----------------------------------------------------------------------
!
!    First time through the time loop, calculate grid
!    transformation info.
!
!-----------------------------------------------------------------------
!
      IF(first_time) THEN
        CALL prepare_interpolation( nx,ny,nz,nx_ext,ny_ext,nz_ext,      &
                   use_wrf_grid, x,y,                                   &
                   iproj_ext,scale_ext,latnot_ext,trlon_ext,x0_ext,y0_ext, &
                   lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext, &
                   xscl,yscl,xa_ext,ya_ext,xs2d,ys2d,iscl,jscl,         &
                   xu2d,yu2d,iu,ju,xv2d,yv2d,iv,jv,                     &
                   iextmn,iextmx,jextmn,jextmx,                         &
                   x_ext,y_ext,xu_ext,yu_ext,xv_ext,yv_ext,             &
                   dxfld,dyfld,rdxfld,rdyfld,                           &
                   dxfldu,dyfldu,rdxfldu,rdyfldu,                       &
                   dxfldv,dyfldv,rdxfldv,rdyfldv,                       &
                   tem1, tem2, tem3, tem4, tem5, tem6,                  &
                   latdiag,londiag,znw,trn_ext,zp_ext,                  &
                   hgt_ext,p_ext,t_ext,pt_ext,qv_ext,u_ext,v_ext,       &
                   istatus )

        IF (istatus /= 0) CALL arpsstop('ERROR: while setting interpolation parameters.',1)
      END IF

      IF (use_wrf_grid /= 1) THEN  ! Otherwise, do not need to rotate U,V
                                   ! because they are on the same grid
!-----------------------------------------------------------------------
!
!  External data comes in oriented to true north.
!  Change that to u,v in the new coordinate system.
!
!-----------------------------------------------------------------------
!
        DO k=1,nz_ext
          CALL uvetomp(nx_ext,ny_ext,u_ext(1,1,k),vatu_ext(1,1,k),      &
                       lonu_ext,tem1_ext(1,1,k),tem2_ext(1,1,k))
        END DO
        u_ext = tem1_ext
        DO k=1,nz_ext
          CALL uvetomp(nx_ext,ny_ext,uatv_ext(1,1,k),v_ext(1,1,k),      &
                       lonv_ext,tem1_ext(1,1,k),tem2_ext(1,1,k))
        END DO
        v_ext = tem2_ext

      END IF
!
!-----------------------------------------------------------------------
!
!  Process height data
!
!  Interpolate the external heights horizontally to the
!  ARPS x,y.  This is for scalar pts and W points.
!
!-----------------------------------------------------------------------
!
      IF(use_wrf_grid /= 1) THEN

        DO k=1,nz_ext         ! W points
          CALL fldint2d(nx,ny,nx_ext,ny_ext,                            &
                        1,nx,1,ny,                                      &
                        1,nx_ext-1,1,ny_ext-1,                          &
                        iorder,xs2d,ys2d,zp_ext(:,:,k),                 &
                        x_ext,y_ext,iscl,jscl,                          &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext,tem2_ext,tem3_ext,                     &
                        zp_tmp(:,:,k))
        END DO

        DO k=1,nz_ext         ! scalar points
          CALL fldint2d(nx,ny,nx_ext,ny_ext,                            &
                        1,nx,1,ny,                                      &
                        1,nx_ext-1,1,ny_ext-1,                          &
                        iorder,xs2d,ys2d,hgt_ext(:,:,k),                &
                        x_ext,y_ext,iscl,jscl,                          &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext,tem2_ext,tem3_ext,                     &
                        zps_tmp(:,:,k))
        END DO

        ! NOTE:  zpsoil_ext is defined as soil depth!!!
        DO k=1,nzsoil_ext
          CALL fldint2d(nx,ny,nx_ext,ny_ext,                            &
                        1,nx,1,ny,                                      &
                        1,nx_ext-1,1,ny_ext-1,                          &
                        iorder,xs2d,ys2d,zpsoil_ext(:,:,k),             &
                        x_ext,y_ext,iscl,jscl,                          &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext,tem2_ext,tem3_ext,                     &
                        zpsoil_tmp(:,:,k))
        END DO

      ELSE

        DO k = 1,nz_ext
          DO j = 2,ny-2
            DO i = 2,nx-2
              zp_tmp (i,j,k) = zp_ext(i-1,j-1,k)
              zps_tmp(i,j,k) = hgt_ext(i-1,j-1,k)
            END DO
          END DO
        END DO

        DO k = 1,nzsoil_ext
          DO j = 2,ny-2
            DO i = 2,nx-2
              zpsoil_tmp(i,j,k) = zpsoil_ext(i-1,j-1,k)
            END DO
          END DO
        END DO

        !
        !  Supply data at the edge points (zero gradient, where missing)
        !
        CALL edgfill(zp_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)
        CALL edgfill(zps_tmp,1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)
        CALL edgfill(zpsoil_tmp,1,nx,2,nx-2,1,ny,2,ny-2,                &
                     1,nzsoil_ext,1,nzsoil_ext)
      END IF

      IF(first_time) THEN
!
!-----------------------------------------------------------------------
!
!  Interpolate external terrain to arps grid (if wrfexttrnopt=1 or 3).
!
!-----------------------------------------------------------------------
!
        IF (wrfexttrnopt == 1 .OR. wrfexttrnopt == 3) THEN
                              ! set arps terrain to match external terrain

          IF(use_wrf_grid /= 1) THEN
            CALL mkarps2d(nx_ext,ny_ext,nx,ny,                          &
                          iorder,iscl,jscl,x_ext,y_ext,                 &
                          xs2d,ys2d,trn_ext,htrn_tmp,                   &
                          dxfld,dyfld,rdxfld,rdyfld,                    &
                          tem1_ext(:,:,1),tem1_ext(:,:,2),              &
                          tem1_ext(:,:,3),tem1_ext(:,:,4))
          ELSE
            DO j = 2,ny-2
              DO i = 2,nx-2
                htrn_tmp(i,j) = trn_ext(i-1,j-1)
              END DO
            END DO
            CALL edgfill(htrn_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
          END IF

          IF (wrfexttrnopt == 1) THEN
            hterain(:,:) = htrn_tmp(:,:)
          ELSE
            ntmergeinv = 1.d0/extntmrg
            DO j = 1,ny
              DO i = 1,nx
                idist = max(0,min(extntmrg,i-2,nx-2-i,j-2,ny-2-j))
                mfac = idist*ntmergeinv
                hterain(i,j) = (1.d0-mfac)*htrn_tmp(i,j)                &
                             + mfac*hterain(i,j)
              END DO
            END DO
          END IF

          IF( mp_opt > 0) THEN
            CALL mpsendrecv1dew(hterain,nx,ny,4,4,0,tem1)
            CALL mpsendrecv1dns(hterain,nx,ny,4,4,0,tem1)
          END IF
          ternopt = -1
          CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                  &
                      hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,         &
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

          IF( dirname /= ' ' ) THEN
            temchar = ternfn
            ternfn = dirname(1:ldirnam)//temchar
            lternfn  = lternfn + ldirnam
          END IF

          CALL fnversn(ternfn, lternfn )

          IF(myproc == 0) PRINT *, 'Write terrain data to ',ternfn(1:lternfn)

          IF(mp_opt >0 .AND. joindmp(FINDX_T) > 0) THEN
            CALL writjointrn(nx,ny,ternfn(1:lternfn), dx,dy,            &
                         mapproj,trulat1,trulat2,trulon,sclfct,         &
                         ctrlat,ctrlon,hterain)
          ELSE
            IF (mp_opt > 0) THEN
              tmpstr = ternfn
              CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                      &
                             0,0,0,lvldbg,ternfn,istatus)
              lternfn = LEN_TRIM(ternfn)
            END IF

            ! blocking inserted for ordering i/o for message passing
            DO i=0,nprocs-1,dumpstride
              IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN
                CALL writtrn(nx,ny,ternfn(1:lternfn), dx,dy,            &
                             mapproj,trulat1,trulat2,trulon,sclfct,     &
                             ctrlat,ctrlon,hterain)
                IF(terndmp == 1) CALL trncntl(nx,ny,ternfn(1:lternfn),x,y)
              END IF
              IF (mp_opt > 0) CALL mpbarrier
            END DO

          END IF

        END IF

!
!-----------------------------------------------------------------------
!
!  Set up z grid at scalar vertical levels.
!
!-----------------------------------------------------------------------
!
        onvf = 0
        CALL avgz(zp, onvf, nx,ny,nz, 1,nx, 1,ny, 1,nz-1, zps)

        DO j=1,ny-1
          DO i=1,nx-1
            zps(i,j,nz)=2.*zps(i,j,nz-1) - zps(i,j,nz-2)
          END DO
        END DO

        CALL edgfill(zps, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
        CALL edgfill(zp,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

      END IF     ! first_time
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
          END DO
        END DO
      END DO
      CALL mpmax0(qvmax,qvmin)

      IF(myproc == 0) WRITE(6,'(1x,a,/,3x,2(a,g15.8))')                 &
        'RHstar diagnostic information:',                               &
        'qv_ext min = ',qvmin, ', qv_ext max = ',qvmax
!
!-----------------------------------------------------------------------
!
!  Calculate base-state sounding (vertical profile)
!
!-----------------------------------------------------------------------
!
      CALL extmnsnd(nx,ny,nx_ext,ny_ext,nz_ext,lvlprof,2,                 &
                    iextmn,iextmx,jextmn,jextmx,1,nz_ext,                 &
                    xscl,yscl,xa_ext,ya_ext,                              &
                    hgt_ext,tem5_ext,t_ext,                               &
                    qv_ext,u_ext,v_ext,                                   &
                    zsnd,plsnd,psnd,tsnd,ptsnd,rhssnd,qvsnd,              &
                    rhosnd,usnd,vsnd)
!
!-----------------------------------------------------------------------
!
!  Process Pressure data
!
!-----------------------------------------------------------------------
!
      iprtopt=1
      IF(use_wrf_grid /= 1) THEN
        CALL mkarpsvlz(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
                       iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,      &
                       hgt_ext,zps_tmp,xs2d,ys2d,zps,p_ext,               &
                       zsnd,plsnd,pbar,pprt,                              &
                       dxfld,dyfld,rdxfld,rdyfld,                         &
                       tem1_ext,tem2_ext,tem3_ext,                        &
                       tem4_ext)
      ELSE
        DO k = 1, nz_ext
          DO j = 2,ny-2
            DO i = 2,nx-2
              var_tmp(i,j,k) = p_ext(i-1,j-1,k)
            END DO
          END DO
        END DO
        CALL edgfill(var_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)

        CALL vintrpvlz_wrf(nx,ny,nz_ext,nz,lvlprof,iprtopt,intropt,     &
                           zps_tmp,zps,var_tmp,zsnd,plsnd,pbar,pprt,    &
                           tem1_tmp)
      END IF

      IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc,wbc,0,tem1)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc,sbc,0,tem1)
      END IF

      DO ksmth=1,nsmooth
        DO k=1,nz
          CALL smooth9p(pprt(:,:,k), nx,ny, 1,nx,1,ny,0, tem1)
        END DO

        IF (nsmooth > ksmth .AND. mp_opt > 0) THEN
                              ! Ensure output is good for followed smoothing
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
      iprtopt=1
      IF(use_wrf_grid /= 1) THEN
        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,           &
                       iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,    &
                       hgt_ext,zps_tmp,xs2d,ys2d,zps,pt_ext,            &
                       zsnd,ptsnd,ptbar,ptprt,                          &
                       dxfld,dyfld,rdxfld,rdyfld,                       &
                       tem1_ext,tem2_ext,tem3_ext,                      &
                       tem4_ext)
      ELSE
        DO k = 1, nz_ext
          DO j = 2,ny-2
            DO i = 2,nx-2
              var_tmp(i,j,k) = pt_ext(i-1,j-1,k)
            END DO
          END DO
        END DO
        CALL edgfill(var_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)

        CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,          &
                           iprtopt,intropt,zps_tmp,zps,var_tmp,         &
                           zsnd,ptsnd,ptbar,ptprt,tem1_tmp)
      END IF

      IF (mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(ptprt,nx,ny,nz,ebc,wbc,0,tem1)
        CALL mpsendrecv2dns(ptprt,nx,ny,nz,nbc,sbc,0,tem1)
      END IF

      DO ksmth=1,nsmooth
        CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,                       &
                      smfct1,zps,ptprt,tem1,ptprt)
      END DO
!
!-----------------------------------------------------------------------
!
!  Process qv data
!
!-----------------------------------------------------------------------
!
      !
      ! QV stores RHstar data now
      !

      iprtopt=0        ! produce total field instead of only perturbation
      IF(use_wrf_grid /= 1) THEN
        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,           &
                       iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,    &
                       hgt_ext,zps_tmp,xs2d,ys2d,zps,qv_ext,            &
                       zsnd,rhssnd,qvbar,qv,                            &
                       dxfld,dyfld,rdxfld,rdyfld,                       &
                       tem1_ext,tem2_ext,tem3_ext,                      &
                       tem4_ext)

      ELSE
        DO k = 1, nz_ext
          DO j = 2,ny-2
            DO i = 2,nx-2
              var_tmp(i,j,k) = qv_ext(i-1,j-1,k)
            END DO
          END DO
        END DO
        CALL edgfill(var_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)

        CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,          &
                           iprtopt,intropt,zps_tmp,zps,var_tmp,         &
                           zsnd,rhssnd,qvbar,qv,tem1_tmp)
      END IF

      IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(qv,nx,ny,nz,ebc,wbc,0,tem1)
        CALL mpsendrecv2dns(qv,nx,ny,nz,nbc,sbc,0,tem1)
      END IF

      DO ksmth=1,nsmooth
        CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,0,smfct1,zps,qv,tem1,qv)
      END DO

     !
     !  Convert rhstar back to qv for writing, further calculations
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
      IF(myproc == 0) THEN
        PRINT *, 'External QV diagnostic information:'
        PRINT *, ' qv_ext min = ',qvmin, ' qv_ext max = ',qvmax
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
      IF(myproc == 0) THEN
        PRINT *, 'ARPS QV diagnostic information:'
        PRINT *, ' qv     min = ',qvmin, ' qv     max = ',qvmax
      END IF
!
!-----------------------------------------------------------------------
!
!  Process hydrometeor data.
!  If first data point is less than or equal to -1, then
!  no qc info is provided...set to zero.
!
!-----------------------------------------------------------------------
!
      DO nq = 1,nscalar
        !IF(qscalar_ext(1,1,1,nq) > -1.0) THEN
        IF(qscalar_ext(1,1,1,nq) > -21.0) THEN  ! take consideration of
                                                !  -20 dBZ min for P_REF
          iprtopt=0
          IF(use_wrf_grid /= 1) THEN
            CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,         &
                           iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,  &
                           hgt_ext,zps_tmp,xs2d,ys2d,zps,qscalar_ext(:,:,:,nq), &
                           zsnd,dumsnd,tem1,qscalar(:,:,:,nq),            &
                           dxfld,dyfld,rdxfld,rdyfld,                     &
                           tem1_ext,tem2_ext,tem3_ext,                    &
                           tem4_ext)
          ELSE
            DO k = 1, nz_ext
              DO j = 2,ny-2
                DO i = 2,nx-2
                  var_tmp(i,j,k) = qscalar_ext(i-1,j-1,k,nq)
                END DO
              END DO
            END DO
            CALL edgfill(var_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)

            CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,        &
                               iprtopt,intropt,zps_tmp,zps,var_tmp,       &
                               zsnd,dumsnd,tem1,qscalar(:,:,:,nq),tem1_tmp)

          END IF

          IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
            CALL acct_interrupt(mp_acct)
            CALL mpsendrecv2dew(qscalar(:,:,:,nq),nx,ny,nz,ebc,wbc,0,tem1)
            CALL mpsendrecv2dns(qscalar(:,:,:,nq),nx,ny,nz,nbc,sbc,0,tem1)
          END IF

          DO ksmth=1,nsmooth
            CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,0,smfct1,zps,qscalar(:,:,:,nq),tem1,qscalar(:,:,:,nq))
          END DO

        END IF
      END DO
!
!-----------------------------------------------------------------------
!
!  Process density which has been stuffed into tem5_ext array.
!  We really only need rhobar, so set iprtopt=1 and pass
!  a temporary array in place of rhoprt.
!
!-----------------------------------------------------------------------
!
       DO k=1,nz_ext
         DO j=1,ny_ext
           DO i=1,nx_ext
             tv_ext=t_ext(i,j,k) /                                      &
                    ((1.0+qv_ext(i,j,k))*                               &
                    (1.0-(qv_ext(i,j,k)/(rddrv+qv_ext(i,j,k)))))
             tem5_ext(i,j,k)=p_ext(i,j,k)/(rd*tv_ext)
           END DO
         END DO
       END DO

       iprtopt=1
       IF(use_wrf_grid /= 1) THEN
         CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,          &
                        iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,   &
                        hgt_ext,zps_tmp,xs2d,ys2d,zps,tem5_ext,         &
                        zsnd,rhosnd,rhobar,tem1,                        &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext,tem2_ext,tem3_ext,                     &
                        tem4_ext)
       ELSE
         DO k = 1, nz_ext
           DO j = 2,ny-2
             DO i = 2,nx-2
               var_tmp(i,j,k) = tem5_ext(i-1,j-1,k)
             END DO
           END DO
         END DO
         CALL edgfill(var_tmp, 1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)
         CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,         &
                            iprtopt,intropt,zps_tmp,zps,var_tmp,        &
                            zsnd,rhosnd,rhobar,tem1,tem1_tmp)
       END IF

       IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
         CALL acct_interrupt(mp_acct)
         CALL mpsendrecv2dew(rhobar,nx,ny,nz,ebc,wbc,0,tem1)
         CALL mpsendrecv2dns(rhobar,nx,ny,nz,nbc,sbc,0,tem1)
       END IF

       DO ksmth=1,nsmooth
         CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,smfct1,zps,         &
                       rhobar,tem1,rhobar)
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
!              qtot=qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qh(i,j,k)
              qtot=0.0
              DO nq=1,nscalarq
               qtot=qtot+qscalar(i,j,k,nq)
              END DO
              tem2(i,j,k)=j3(i,j,k)*rhobar(i,j,k)*                      &
                          g*( (ptprt(i,j,k)/ptbar(i,j,k)) +             &
                              (qvprt/(rddrv+qvbar(i,j,k))) -            &
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

          DO i=1,nx
            DO j=1,ny
              tem1(i,j,1)=pprt(i,j,1)
              DO k=2,nz-2
                tem1(i,j,k)=                                            &
                      ( (1. - (pconst*j3(i,j,k-1)/csndsq(i,j,k-1)) )*   &
                          tem1(i,j,k-1) +                               &
                          dz*tem3(i,j,k) ) /                            &
                          (1. + (pconst*j3(i,j,k)/csndsq(i,j,k)) )
                IF(j == 26 .AND. MOD(i,10) == 0) THEN
                  IF(k == nz-2) PRINT *, '            Point i= ',i, '  j=26'
                  PRINT *, ' k,pprt,tem1,diff =',                       &
                      k,pprt(i,j,k),tem1(i,j,k),                        &
                      (tem1(i,j,k)-pprt(i,j,k))
                END IF
                pprt(i,j,k)=tem1(i,j,k)
              END DO
              pprt(i,j,nz-1)=pprt(i,j,nz-2)
            END DO
          END DO

        ELSE IF(hydradj == 2) THEN

          DO i=1,nx
            DO j=1,ny
              tem1(i,j,nz-1)=pprt(i,j,nz-1)
              DO k=nz-2,2,-1
                tem1(i,j,k)=                                            &
                      ( (1.+ (pconst*j3(i,j,k+1)/csndsq(i,j,k+1)) )*    &
                          tem1(i,j,k+1) -                               &
                          dz*tem3(i,j,k+1) ) /                          &
                          (1.- (pconst*j3(i,j,k  )/csndsq(i,j,k  )) )
                IF(j == 26 .AND. MOD(i,10) == 0) THEN
                  IF(k == nz-2) PRINT *, '            Point i= ',i, '  j=26'
                  PRINT *, ' k,pprt,tem1,diff =',                       &
                      k,pprt(i,j,k),tem1(i,j,k),                        &
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
              temp = (ptbar(i,j,k)+ptprt(i,j,k))*                       &
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
              tem1(i,j,k)=tem1(i,j,k-1)*                                &
                         EXP(pconst*(zps(i,j,k)-zps(i,j,k-1))/tvbar)
              pprt(i,j,k)=tem1(i,j,k)-pbar(i,j,k)
            END DO
          END DO
        END DO
        DO j=1,ny
          DO i=1,nx
            tvbar=0.5*(tem2(i,j,2)+tem2(i,j,1))
            tem1(i,j,1)=tem1(i,j,2)*                                    &
                       EXP(pconst*(zps(i,j,1)-zps(i,j,2))/tvbar)
            pprt(i,j,1)=tem1(i,j,1)-pbar(i,j,1)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!    Process U wind data
!
!-----------------------------------------------------------------------
!
      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(zps_tmp,nx,ny,nz_ext,ebc,wbc,0,tem1_tmp)
        CALL mpsendrecv2dew(zps,    nx,ny,nz,    ebc,wbc,0,tem3)
      END IF
      CALL avgsu(zps_tmp,nx,ny,nz_ext,1,ny,1,nz_ext,avgzps_tmp,tem1_tmp)
      CALL avgsu(zps,    nx,ny,nz,    1,ny,1,nz,    tem2,      tem3)

      IF (loc_x == 1) THEN  ! West boundary of the whole domain
                            ! avgsu calls bcsu and set it differently
                            ! as what we need.
        DO k=1,nz_ext
          DO j=1,ny
            avgzps_tmp(1,j,k) = zps_tmp(1,j,k)
          END DO
        END DO
        DO k=1,nz
          DO j=1,ny
            tem2(1,j,k) = zps(1,j,k)
          END DO
        END DO
      END IF

      IF (loc_x == nproc_x) THEN  ! east boundary of the whold domain
                                  ! avgsu calls bcsu and set it differently
                                  ! as what we need.
        DO k=1,nz_ext
          DO j=1,ny
            avgzps_tmp(nx,j,k)= zps_tmp(nx-1,j,k)
          END DO
        END DO
        DO k=1,nz
          DO j=1,ny
            tem2(nx,j,k)= zps(nx-1,j,k)
          END DO
        END DO
      END IF

      iprtopt=0
      IF(use_wrf_grid /= 1) THEN
        ! u is staggered as in arps
        CALL avgsu(hgt_ext,nx_ext,ny_ext,nz_ext,1,ny_ext,1,nz_ext,      &
                   tem5_ext,tem3_ext)
        tem5_ext(     1,:,:) = hgt_ext(       1,:,:)
        tem5_ext(nx_ext,:,:) = hgt_ext(nx_ext-1,:,:)

        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,           &
                       iorder,iprtopt,intropt,iu,ju,xu_ext,yu_ext,      &
                       tem5_ext,avgzps_tmp,xu2d,yu2d,tem2,u_ext,        &
                       zsnd,usnd,ubar,u,                                &
                       dxfldu,dyfldu,rdxfldu,rdyfldu,                   &
                       tem1_ext,tem2_ext,tem3_ext,                      &
                       tem4_ext)
      ELSE
        DO k = 1, nz_ext
          DO j = 2,ny-2
            DO i = 2,nx-1
              var_tmp(i,j,k) = u_ext(i-1,j-1,k)
            END DO
          END DO
        END DO
        CALL edgfill(var_tmp,1,nx,2,nx-1,1,ny,2,ny-2,1,nz_ext,1,nz_ext)
        CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,          &
                           iprtopt,intropt,avgzps_tmp,tem2,var_tmp,     &
                           zsnd,usnd,ubar,u,tem1_tmp)
      END IF

      IF (mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(u,nx,ny,nz,ebc,wbc,1,tem1)
        CALL mpsendrecv2dns(u,nx,ny,nz,nbc,sbc,1,tem1)
      END IF

      DO ksmth=1,nsmooth
        CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,1,smfct1,tem2,u,tem1,u)
      END DO

      CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'umin final= ', amin,',  umax final=',amax
!
!-----------------------------------------------------------------------
!
!  Process v component
!
!-----------------------------------------------------------------------
!
      IF (mp_opt > 0) THEN
        CALL mpsendrecv2dns(zps_tmp,nx,ny,nz_ext,nbc,sbc,0,tem1_tmp)
        CALL mpsendrecv2dns(zps,    nx,ny,nz,    nbc,sbc,0,tem3)
      END IF
      CALL avgsv(zps_tmp,nx,ny,nz_ext,1,nx,1,nz_ext,avgzps_tmp,tem1_tmp)
      CALL avgsv(zps,nx,ny,nz,        1,nx,1,nz,    tem2,      tem3)

      IF (loc_y == 1) THEN     ! south boundary
        DO k=1,nz_ext
          DO i=1,nx
            avgzps_tmp(i,1,k) = zps_tmp(i,1,k)
          END DO
        END DO
        DO k=1,nz
          DO i=1,nx
            tem2(i,1,k)=zps(i,1,k)
          END DO
        END DO
      END IF

      IF (loc_y == nproc_y) THEN     ! north boundary
        DO k=1,nz_ext
          DO i=1,nx
            avgzps_tmp(i,ny,k) = zps_tmp(i,ny-1,k)
          END DO
        END DO
        DO k=1,nz
          DO i=1,nx
            tem2(i,ny,k) = zps(i,ny-1,k)
          END DO
        END DO
      END IF

      iprtopt=0
      IF(use_wrf_grid /= 1) THEN
        ! v is staggered as in arps
        CALL avgsv(hgt_ext,nx_ext,ny_ext,nz_ext,1,nx_ext,1,nz_ext,      &
                   tem5_ext,tem3_ext)
        tem5_ext(1:nx_ext,     1,1:nz_ext) = hgt_ext(1:nx_ext,       1,1:nz_ext)
        tem5_ext(1:nx_ext,ny_ext,1:nz_ext) = hgt_ext(1:nx_ext,ny_ext-1,1:nz_ext)

        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,           &
                       iorder,iprtopt,intropt,iv,jv,xv_ext,yv_ext,      &
                       tem5_ext,avgzps_tmp,xv2d,yv2d,tem2,v_ext,        &
                       zsnd,vsnd,vbar,v,                                &
                       dxfldv,dyfldv,rdxfldv,rdyfldv,                   &
                       tem1_ext,tem2_ext,tem3_ext,                      &
                       tem4_ext)
      ELSE
        DO k = 1, nz_ext
          DO j = 2,ny-1
            DO i = 2,nx-2
              var_tmp(i,j,k) = v_ext(i-1,j-1,k)
            END DO
          END DO
        END DO
        CALL edgfill(var_tmp,1,nx,2,nx-2,1,ny,2,ny-1,1,nz_ext,1,nz_ext)
        CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,          &
                           iprtopt,intropt,avgzps_tmp,tem2,var_tmp,     &
                           zsnd,vsnd,vbar,v,tem1_tmp)
      END IF

      IF (mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(v,nx,ny,nz,ebc,wbc,2,tem1)
        CALL mpsendrecv2dns(v,nx,ny,nz,nbc,sbc,2,tem1)
      END IF

      DO ksmth=1,nsmooth
        CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,2,smfct1,tem2,v,tem1,v)
      END DO

      CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'vmin final= ', amin,',  vmax final=',amax
!
!-----------------------------------------------------------------------
!
!  Process 4D surface fields if required (sfcout = 1)
!
!-----------------------------------------------------------------------
!
      wetcout  = 1
      tsoilout = 1
      qsoilout = 1

      ! Horizontally interpolations

      IF(use_wrf_grid /= 1) THEN

        DO k=1,nzsoil_ext
          ! Soil temperature
          CALL fldint2d(nx,ny,nx_ext,ny_ext,                              &
                        1,nx,1,ny,                                        &
                        1,nx_ext-1,1,ny_ext-1,                            &
                        iorder,xs2d,ys2d,tsoil_ext(:,:,k,0),              &
                        x_ext,y_ext,iscl,jscl,                            &
                        dxfld,dyfld,rdxfld,rdyfld,                        &
                        tem1_ext,tem2_ext,tem3_ext,                       &
                        tsoil_tmp(:,:,k))
          ! Soil moisture
          CALL fldint2d(nx,ny,nx_ext,ny_ext,                              &
                        1,nx,1,ny,                                        &
                        1,nx_ext-1,1,ny_ext-1,                            &
                        iorder,xs2d,ys2d,qsoil_ext(:,:,k,0),              &
                        x_ext,y_ext,iscl,jscl,                            &
                        dxfld,dyfld,rdxfld,rdyfld,                        &
                        tem1_ext,tem2_ext,tem3_ext,                       &
                        qsoil_tmp(:,:,k))
        END DO
        ! wet canopy
        CALL fldint2d(nx,ny,nx_ext,ny_ext,                              &
                      1,nx,1,ny,                                        &
                      1,nx_ext-1,1,ny_ext-1,                            &
                      iorder,xs2d,ys2d,wetcanp_ext(:,:,0),              &
                      x_ext,y_ext,iscl,jscl,                            &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext,tem2_ext,tem3_ext,                       &
                      wetcanp(:,:,0))

      ELSE
        DO k = 1, nzsoil_ext
          DO j = 2,ny-2
            DO i = 2,nx-2
              tsoil_tmp(i,j,k) = tsoil_ext(i-1,j-1,k,0)
              qsoil_tmp(i,j,k) = qsoil_ext(i-1,j-1,k,0)
            END DO
          END DO
        END DO
        DO j = 2,ny-2
          DO i = 2,nx-2
            wetcanp(i,j,0) = wetcanp_ext(i-1,j-1,0)
          END DO
        END DO
        CALL edgfill(tsoil_tmp,1,nx,2,nx-2,1,ny,2,ny-2,                 &
                     1,nzsoil_ext,1,nzsoil_ext)
        CALL edgfill(qsoil_tmp,1,nx,2,nx-2,1,ny,2,ny-2,                 &
                     1,nzsoil_ext,1,nzsoil_ext)
        CALL edgfill(wetcanp(:,:,0),1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      END IF

      ! Cut-off bad values for qsoil_tmp and wetcanp
      DO j = 1,ny
        DO i = 1,nx
          DO k = 1, nzsoil_ext
          if(qsoil_tmp(i,j,k) < 0.02) qsoil_tmp(i,j,k) = 0.02
          if(qsoil_tmp(i,j,k) > 1.00) qsoil_tmp(i,j,k) = 1.00
          END DO
          if(wetcanp(i,j,0) < 0.0) wetcanp(i,j,0) = 0.0
        END DO
      END DO

      ! First convert zpsoil to soil depth.
      DO k = 1,nzsoil
        DO j = 1,ny
          DO i = 1,nx
            zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
          END DO
        END DO
      END DO

      CALL vintrpsoil_wrf(nx,ny,nzsoil_ext,nzsoil,soilmodel_option,     &
                          zpsoil_tmp,zpsoil,tsoil_tmp,qsoil_tmp,        &
                          tsoil(:,:,:,0),qsoil(:,:,:,0))

      CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp_ext,nstyp,                 &
                          tsoil,qsoil,wetcanp)
      DO is = 1,nstyp
        CALL edgfill(tsoil(:,:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,           &
                     1,nzsoil,1,nzsoil)
        CALL edgfill(qsoil(:,:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,           &
                     1,nzsoil,1,nzsoil)
        CALL edgfill(wetcanp(:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      END DO

      ! Convert zpsoil back from soil depth to MSL
      DO k = 1,nzsoil
        DO j = 1,ny
          DO i = 1,nx
            zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
          END DO
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!  Process 3D scalar fields if required (tkeout = 1)
!
!-----------------------------------------------------------------------
!
      IF (tkeout > 0 .AND. tke_ext(1,1,1) >= 0.0 ) THEN

        ! Horizontally interpolations

        IF(use_wrf_grid /= 1) THEN

          DO k=1,nz_ext
            CALL fldint2d(nx,ny,nx_ext,ny_ext,                          &
                          1,nx,1,ny,                                    &
                          1,nx_ext-1,1,ny_ext-1,                        &
                          iorder,xs2d,ys2d,tke_ext(:,:,k),              &
                          x_ext,y_ext,iscl,jscl,                        &
                          dxfld,dyfld,rdxfld,rdyfld,                    &
                          tem1_ext,tem2_ext,tem3_ext,                   &
                          var_tmp(:,:,k))
          END DO

        ELSE
          DO k = 1, nz_ext
            DO j = 2,ny-2
              DO i = 2,nx-2
                var_tmp(i,j,k) = tke_ext(i-1,j-1,k)
              END DO
            END DO
          END DO
          CALL edgfill(var_tmp,1,nx,2,nx-2,1,ny,2,ny-2,1,nz_ext,1,nz_ext)
        END IF

        iprtopt=0
        CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext-1,lvlprof,          &
                          iprtopt,intropt,zps_tmp,zps,var_tmp,          &
                          zsnd,dumsnd,tem1,tke,tem1_tmp)

        IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
          CALL acct_interrupt(mp_acct)
          CALL mpsendrecv2dew(tke,nx,ny,nz,ebc,wbc,0,tem1)
          CALL mpsendrecv2dns(tke,nx,ny,nz,nbc,sbc,0,tem1)
        END IF

        WHERE(tke < 0.0) tke = 0.0

        DO ksmth=1,nsmooth
          CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,0,smfct1,zps,tke,tem1,tke)
        END DO

      END IF
!-----------------------------------------------------------------------
!
! Processing Snow depth, soil type and vegtation type etc.
!
!-----------------------------------------------------------------------

    ! Snow depth using nearest neighbor

      snowdout = 1

      IF(use_wrf_grid /= 1) THEN
        DO j=1,ny
          DO i=1,nx
            dmin=((xs2d(i,j)-x_ext(1))*(xs2d(i,j)-x_ext(1)) +             &
                  (ys2d(i,j)-y_ext(1))*(ys2d(i,j)-y_ext(1)))
            isnow=1
            jsnow=1
            DO jj=jscl(i,j),MAX(jscl(i,j)+1,ny_ext)
              DO ii=iscl(i,j),MAX(iscl(i,j)+1,nx_ext)
                dd=((xs2d(i,j)-x_ext(ii))*(xs2d(i,j)-x_ext(ii)) +         &
                    (ys2d(i,j)-y_ext(jj))*(ys2d(i,j)-y_ext(jj)))
                IF(dd < dmin) THEN
                  dmin=dd
                  isnow=ii
                  jsnow=jj
                END IF
              END DO
            END DO
            snowdpth(i,j)  = snowdpth_ext(isnow,jsnow)
            soiltyp(i,j,1) = soiltyp_ext(isnow,jsnow,1)
            vegtyp(i,j)    = vegtyp_ext(isnow,jsnow)
            veg(i,j)       = veg_ext(isnow,jsnow)
          END DO
        END DO
      ELSE
        DO j = 2,ny-2
          DO i = 2,nx-2
            snowdpth(i,j) = snowdpth_ext(i-1,j-1)
            soiltyp(i,j,1)= soiltyp_ext(i-1,j-1,1)
            vegtyp(i,j)   = vegtyp_ext(i-1,j-1)
            veg(i,j)      = veg_ext(i-1,j-1)
          END DO
        END DO
        CALL edgfill (snowdpth,1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill (veg,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL iedgfill(soiltyp(:,:,1),1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL iedgfill(vegtyp,  1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      END IF

      DO k = 2,nstyps
        DO j = 1,ny
          DO i = 1,nx
            soiltyp(i,j,k) = -999
          END DO
        END DO
      END DO

      DO j = 1,ny
        DO i = 1,nx
          roufns(i,j) = rfnstbl(vegtyp(i,j))
          lai(i,j)    = laitbl(vegtyp(i,j))
        END DO
      END DO

      DO k = 1, nsmooth
        DO j = 1,ny
          DO i = 1,nx
            roufns(i,j) = LOG( MAX(0.00001,roufns(i,j)) )
          END DO
        END DO
        CALL smooth25p( lai,    nx,ny,1,nx,1,ny, tem2 ) ! not mpi-ied
        CALL smooth25p( roufns, nx,ny,1,nx,1,ny, tem2 ) ! not mpi-ied
        DO j = 1,ny
          DO i = 1,nx
            roufns(i,j) = EXP ( roufns(i,j) )
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!
! Ouput soil variables
!
!-----------------------------------------------------------------------


    IF (wetcout  == 1 .OR. snowdout == 1 .OR.                           &
        tsoilout == 1 .OR. qsoilout == 1 ) THEN

      zpsoilout = 1
      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
      lfn = lfnkey + 9 + tmstrln

      IF( dirname /= ' ' ) THEN
        temchar = soiloutfl
        soiloutfl = dirname(1:ldirnam)//temchar
        lfn  = lfn + ldirnam
      END IF

      IF (soildmp > 0) THEN
        !CALL fnversn(soiloutfl, lfn)

        !IF(myproc == 0) PRINT *, 'Write soil initial data to ',soiloutfl(1:lfn)

        IF(mp_opt >0 .AND. joindmp(FINDX_S) > 0) THEN
          CALL wrtjoinsoil(nx,ny,nzsoil,nstyp, soiloutfl,dx,dy,zpsoil,  &
                 mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,   &
                 zpsoilout, tsoilout,qsoilout, wetcout,snowdout,        &
                 tsoil,qsoil,wetcanp,snowdpth, soiltyp)
        ELSE
          !IF (mp_opt > 0) THEN
          !  tmpstr = soiloutfl
          !  CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                  &
          !                 0,0,0,lvldbg,soiloutfl,istatus)
          !END IF

          ! blocking inserted for ordering i/o for message passing
          DO i=0,nprocs-1,dumpstride
            IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

              CALL wrtsoil(nx,ny,nzsoil,nstyp, soiloutfl,dx,dy,zpsoil,      &
                     mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,   &
                     zpsoilout, tsoilout,qsoilout, wetcout,snowdout,        &
                     tsoil,qsoil,wetcanp,snowdpth, soiltyp)

              IF (soildmp == 1) CALL soilcntl(nx,ny, nzsoil,zpsoil,soiloutfl,&
                        zpsoilout,tsoilout,qsoilout, wetcout,snowdout,x,y)

            END IF
            IF (mp_opt > 0) CALL mpbarrier
          END DO

        END IF
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
    CALL output_diagnose_arps(nx,ny,nz,xscl,yscl,latdiag,londiag,       &
                              zps,pbar,pprt,ptbar,ptprt,qv,u,v,         &
                              istatus)
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

      IF ( use_wrf_grid /= 1 ) THEN

        iprtopt=0
        CALL mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,           &
                       iorder,iprtopt,intropt,iscl,jscl,x_ext,y_ext,    &
                       zp_ext,zp_tmp,xs2d,ys2d,zp,w_ext,                &
                       zsnd,dumsnd,wbar,w,                              &
                       dxfld,dyfld,rdxfld,rdyfld,                       &
                       tem1_ext,tem2_ext,tem3_ext,                      &
                       tem4_ext)
      ELSE

        DO k = 1, nz_ext
          DO j = 2,ny-1
            DO i = 2,nx-2
              var_tmp(i,j,k) = w_ext(i-1,j-1,k)
            END DO
          END DO
        END DO
        CALL edgfill(var_tmp,1,nx,2,nx-2,1,ny,2,ny-1,1,nz_ext,1,nz_ext)
        CALL vintrpvar_wrf(nx,ny,nz_ext,nz,1,nz_ext,lvlprof,            &
                           iprtopt,intropt,zp_tmp,zp,var_tmp,           &
                           zsnd,dumsnd,wbar,w,tem1_tmp)
      END IF

      IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc,wbc,3,tem1)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc,sbc,3,tem1)
      END IF

      DO ksmth=1,nsmooth
        CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,3,smfct1,zp,w,tem1,w)
      END DO

    ELSE

      CALL adjuvw( nx,ny,nz, u,v,w,wcont,ptprt,ptbar,                   &
                   zp,j1,j2,j3,aj3z,mapfct,rhostr,tem1,                 &
                   wndadj,obropt,obrzero,0,                             &
                   tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

    END IF
!
!-----------------------------------------------------------------------
!
!  Enforce vertical w boundary conditions
!
!-----------------------------------------------------------------------
!
    IF (ext_lbc == 1 .or. ext_vbc == 1) THEN
      CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)
    END IF


    IF (ext_vbc == 1) THEN
      IF (mp_opt > 0) THEN
        CALL mpsendrecv2dew(u,nx,ny,nz,4,4,1,tem4)
        CALL mpsendrecv2dns(v,nx,ny,nz,4,4,2,tem4)
        CALL mpsendrecv2dew(j1,nx,ny,nz,4,4,0,tem4)
        CALL mpsendrecv2dns(j2,nx,ny,nz,4,4,0,tem4)
      END IF
      CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,           &
                rhostr,tem1,tem2,tem3,j1,j2,j3)
    END IF
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
!  Processing t2m, th2m, qv2m, u10m, v10m, raddn
!  Processing wspd10max,w_up_max,w_dn_max,refd_max,up_heli_max,grpl_max
!             ltg1_max,ltg2_max,ltg3_max
!
!-----------------------------------------------------------------------

    CALL make_data2d(use_wrf_grid,iorder,iscl,jscl,x_ext,y_ext,         &
                     xs2d,ys2d,dxfld,dyfld,rdxfld,rdyfld,               &
                     tem1_ext,istatus)

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

    CALL print_maxmin_arps(myproc,nx,ny,nz,nzsoil,mp_physics,           &
                           x, y, z, zp, rhobar, u, v, w,                &
                           ptbar, ptprt, pbar, pprt, qv, qscalar,       &
                           tke, tsoil,qsoil,wetcanp, istatus)

    CALL print_maxmin_data2d_arps( myproc,istatus )

!
!-----------------------------------------------------------------------
!
!  Print the mean sounding that was used in setting the
!  mean ARPS variables.
!
!-----------------------------------------------------------------------
!
    IF(myproc == 0)  THEN

      WRITE(6,'(a/,a,a)')                                               &
          '  Mean sounding for ARPS base state variables',              &
          '  k     p(mb)     z(m)    pt(mb)    T(C)   qv(kg/kg) ',      &
          '  RH %  u(m/s)    v(m/s)'

      DO k=lvlprof,1,-50
        pres = psnd(k)
        temp = ptsnd(k)*((pres/p0)**rddcp)
        rh=AMAX1(0.01,rhmax-(rhssnd(k)*rhssnd(k)))
        qvsat=f_qvsat( pres, temp )
        qvval=rh*qvsat
        WRITE(6,'(i4,f9.1,f9.1,f9.1,f9.1,f9.5,f9.1,f9.1,f9.1)')         &
            k,(0.01*psnd(k)),zsnd(k),ptsnd(k),(temp-273.15),            &
            qvval,(100.*rh),usnd(k),vsnd(k)
      END DO

    END IF
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
    CALL gtdmpfn(runname(1:lfnkey),dirname,ldirnam,curtim,hdmpfmt,      &
                 mgrid,nestgrd, hdmpfn, ldmpf)

    IF (hdmpfmt /= 0) THEN

    tem1=0.

!    ldirnam=LEN(dirname)
!    CALL strlnth( dirname, ldirnam)

    IF ( hdmpfmt == 5 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Program wrf2arps does not support Savi3D format.',           &
          'Reset to binary format, hdmpfmt=1.'
      hdmpfmt = 1
    END IF

    IF ( hdmpfmt == 9 ) GO TO 450

    IF(first_time .AND. ABS(curtim) < 1.0E-5) THEN

      CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,           &
                   mgrid,nestgrd,basdmpfn,lbasdmpf)

      IF(myproc == 0)        &
        PRINT*,'Writing base state history dump ',basdmpfn(1:lbasdmpf)

      mstout0 = mstout
      rainout0= rainout
      prcout0 = prcout
      trbout0 = trbout

      grdbas =1
      mstout =1
      rainout=0
      prcout =0
      trbout =0

      ! blocking inserted for ordering i/o for message passing
      DO i=0,nprocs-1,dumpstride
        IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

          CALL dtadump(nx,ny,nz,nzsoil,nstyp,                           &
                   hdmpfmt,iniotfu,basdmpfn(1:lbasdmpf),grdbas,filcmprs,&
                   u,v,w,ptprt,pprt,qv,qscalar,                         &
                   tke,tem1,tem1,                                       &
                   ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,              &
                   x,y,z,zp,zpsoil,                                     &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,tem1,                                    &
                   tem1,tem1,tem1,                                      &
                   tem1,tem1,                                           &
                   tem1,tem1,tem1,tem1,                                 &
                   tem2,tem3,tem4)
        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

      mstout = mstout0
      rainout= rainout0
      prcout = prcout0
      trbout = trbout0
      first_time = .FALSE.

    END IF  ! first_time

    450 CONTINUE

    IF(myproc == 0) THEN
      PRINT*,'curtim = ', curtim
      WRITE(6,'(1x,a,a)') 'History data dump in file ',hdmpfn(1:ldmpf)
    END IF

    grdbas=0

    ! blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN
        CALL dtadump(nx,ny,nz,nzsoil,nstyp,                             &
                 hdmpfmt,iniotfu,hdmpfn(1:ldmpf),grdbas,filcmprs,       &
                 u,v,w,ptprt,pprt,qv,qscalar,                           &
                 tke,tem1,tem1,                                         &
                 ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,tem1,                                      &
                 tem1,tem1,tem1,                                        &
                 tem1,tem1,                                             &
                 tem1,tem1,tem1,tem1,                                   &
                 tem2,tem3,tem4)
      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

    END IF

!-----------------------------------------------------------------------
!
!  Calculate and write out certian 2D fields for post-processing
!  (both in binary and in GEMPAK format)
!
!-----------------------------------------------------------------------
!
      CALL process_rain( curtim, iabssec, multifile, use_wrf_grid,      &
                  grid_id, io_form, nprocx_in, ncompressx, ncompressy,  &
                  magnitude_processor, rewindyr,myr, dir_extp,          &
                  nz_ext,iorder,iscl,jscl,x_ext,y_ext,xs2d,ys2d,        &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  tem1_ext,istatus )


      IF(hdmpfn(ldmpf-2:ldmpf-2) == '.') ldmpf = ldmpf - 3

      u  = u-ubar
      v  = v-vbar
      qv = qv-qvbar

      CALL process_post_arps( curtim,hdmpfn,ldmpf,mp_physics,           &
                        nz,nzsoil,nstyp,nscalar,x,y,z,zp,               &
                        u,v,w ,ptprt, pprt, qv, qscalar,                &
                        ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,   &
                        soiltyp,stypfrct,vegtyp,lai,roufns,veg,         &
                        tsoil,qsoil,wetcanp,snowdpth,                   &
                        tem1,tem2,tem3,tem4,tem5,tem6,istatus )

      IF (iabssec >= abstimee) THEN   ! exit if time exceeds required in namelist
        exit_early = .TRUE.
        EXIT
      END IF

    END DO     ! itime

    IF (frames_per_outfile(ifile) > 1)                                  &
      CALL close_wrf_file(fHndl,io_form,multifile,.FALSE.,              &
                          ncompressx,ncompressy)

    IF (exit_early) EXIT

  END DO  ! ifile

  CALL io_shutdown(io_form)

!
!-----------------------------------------------------------------------
!
!  Friendly exit message.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a/a,i4,a)')                                 &
          ' ==== Normal succesful completion of WRF2ARPS. ====',        &
          '      Processed',nextdfil,' file(s)'

  IF (mp_opt > 0) CALL mpexit(0)

  STOP
!
!-----------------------------------------------------------------------
!
!  Error status returned from getwrfd
!
!-----------------------------------------------------------------------
!
  999 CONTINUE
  IF(myproc == 0) WRITE(6,'(a,i6)')                                     &
          ' Aborting, error reading external file. istatus=',istatus

  IF (mp_opt > 0) CALL mpexit(0)

  STOP
END PROGRAM wrf2arps

!#######################################################################

SUBROUTINE prepare_interpolation( nx,ny,nz,nx_ext,ny_ext,nz_ext,        &
                   use_wrf_grid, x,y,                                   &
                   iproj_ext,scale_ext,latnot_ext,trlon_ext,x0_ext,y0_ext, &
                   lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext, &
                   xscl,yscl,xa_ext,ya_ext,xs2d,ys2d,iscl,jscl,         &
                   xu2d,yu2d,iu,ju,xv2d,yv2d,iv,jv,                     &
                   iextmn,iextmx,jextmn,jextmx,                         &
                   x_ext,y_ext,xu_ext,yu_ext,xv_ext,yv_ext,             &
                   dxfld,dyfld,rdxfld,rdyfld,                           &
                   dxfldu,dyfldu,rdxfldu,rdyfldu,                       &
                   dxfldv,dyfldv,rdxfldv,rdyfldv,                       &
                   tem1, tem2, tem3, tem4, tem5, tem6,                  &
                   latdiag,londiag,znw,trn_ext,zp_ext,                  &
                   hgt_ext,p_ext,t_ext,pt_ext,qv_ext,u_ext,v_ext,       &
                   istatus )

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny, nz, nx_ext, ny_ext, nz_ext
  INTEGER, INTENT(IN)  :: use_wrf_grid
  INTEGER, INTENT(IN)  :: iproj_ext
  REAL,    INTENT(IN)  :: scale_ext,latnot_ext(2),trlon_ext,x0_ext,y0_ext
  REAL,    INTENT(IN)  :: x(nx), y(ny)
  REAL,    INTENT(IN)  :: lat_ext(nx_ext,ny_ext),  lon_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: latu_ext(nx_ext,ny_ext), lonu_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: latv_ext(nx_ext,ny_ext), lonv_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: iextmn,iextmx,jextmn,jextmx

  REAL,    INTENT(OUT) :: xscl(nx),yscl(ny)
  REAL,    INTENT(OUT) :: xa_ext(nx_ext,ny_ext),ya_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: xs2d(nx,ny),ys2d(nx,ny),xu2d(nx,ny),yu2d(nx,ny),xv2d(nx,ny),yv2d(nx,ny)
  INTEGER, INTENT(OUT) :: iscl(nx,ny),jscl(nx,ny),iu(nx,ny),ju(nx,ny),iv(nx,ny),jv(nx,ny)

  REAL,    INTENT(OUT) :: x_ext(nx_ext), y_ext(ny_ext),                 &
                          xu_ext(nx_ext),yu_ext(ny_ext),                &
                          xv_ext(nx_ext),yv_ext(ny_ext)

  REAL,    INTENT(OUT) :: dxfld(nx_ext), dyfld(ny_ext), rdxfld(nx_ext), rdyfld(ny_ext),   &
                          dxfldu(nx_ext),dyfldu(ny_ext),rdxfldu(nx_ext),rdyfldu(ny_ext),  &
                          dxfldv(nx_ext),dyfldv(ny_ext),rdxfldv(nx_ext),rdyfldv(ny_ext)

  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz), tem2(nx,ny,nz), tem3(nx,ny,nz)
  REAL,    INTENT(INOUT) :: tem4(nx,ny,nz), tem5(nx,ny,nz), tem6(nx,ny,nz)

  ! for diagnostic output only
  REAL,    INTENT(IN)  :: latdiag, londiag
  REAL,    INTENT(IN)  :: znw(nz_ext)
  REAL,    INTENT(IN)  :: trn_ext(nx_ext,ny_ext), zp_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: hgt_ext(nx_ext,ny_ext,nz_ext), p_ext(nx_ext,ny_ext,nz_ext), &
                          t_ext(nx_ext,ny_ext,nz_ext),  pt_ext(nx_ext,ny_ext,nz_ext), &
                          u_ext(nx_ext,ny_ext,nz_ext),   v_ext(nx_ext,ny_ext,nz_ext), &
                          qv_ext(nx_ext,ny_ext,nz_ext)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i,j
  REAL    :: xumax, xumin, yvmin, yvmax

  INTEGER :: iproj_arps
  REAL    :: scale_arps,latnot_arps(2),trlon_arps,xorig_arps,yorig_arps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

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
  CALL getmapr(iproj_arps,scale_arps,latnot_arps,                       &
               trlon_arps,xorig_arps,yorig_arps)
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of external grid.
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL setorig(1,x0_ext,y0_ext)
  DO j=1,ny_ext
    CALL lltoxy(1,1,lat_ext(1,j),lon_ext(1,j),   x_ext(1),y_ext(j))
  END DO
  DO i=1,nx_ext
    CALL lltoxy(1,1,lat_ext(i,1),lon_ext(i,1),   x_ext(i),y_ext(1))
  END DO
  DO j=1,ny_ext
    CALL lltoxy(1,1,latu_ext(1,j),lonu_ext(1,j), xu_ext(1),yu_ext(j))
  END DO
  DO i=1,nx_ext
    CALL lltoxy(1,1,latu_ext(i,1),lonu_ext(i,1), xu_ext(i),yu_ext(1))
  END DO
  DO j=1,ny_ext
    CALL lltoxy(1,1,latv_ext(1,j),lonv_ext(1,j), xv_ext(1),yv_ext(j))
  END DO
  DO i=1,nx_ext
    CALL lltoxy(1,1,latv_ext(i,1),lonv_ext(i,1), xv_ext(i),yv_ext(1))
  END DO

!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS scalar grid in terms of the external grid.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(nx,ny,tem1,tem2,xs2d,ys2d)
  CALL setijloc(nx,ny,nx_ext,ny_ext,xs2d,ys2d,x_ext,y_ext,iscl,jscl)
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS u grid in terms of the external grid.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(nx,ny,tem3,tem4,xu2d,yu2d)
  CALL setijloc(nx,ny,nx_ext,ny_ext,xu2d,yu2d,xu_ext,yu_ext,iu,ju)
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS v grid in terms of the external grid.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(nx,ny,tem5,tem6,xv2d,yv2d)
  CALL setijloc(nx,ny,nx_ext,ny_ext,xv2d,yv2d,xv_ext,yv_ext,iv,jv)

  IF(use_wrf_grid /= 1) THEN

    xumin=xu2d(1,1)
    xumax=xu2d(1,1)
    DO j=1,ny-1
      xumin=MIN(xumin,xu2d(1 ,j))
      xumax=MAX(xumax,xu2d(nx,j))
    END DO

    yvmin=yv2d(1,1)
    yvmax=yv2d(1,1)
    DO i=1,nx-1
      yvmin=MIN(yvmin,yv2d(i,1 ))
      yvmax=MAX(yvmax,yv2d(i,ny))
    END DO

    IF(xumin < xu_ext(1) .OR. xumax > xu_ext(nx_ext).OR.                &
       yvmin < yv_ext(1) .OR. yvmax > yv_ext(ny_ext)) THEN
      WRITE(6,'(a/a)')                                                  &
        ' ARPS domain extends outside available external data',         &
        ' domain, ext2arps aborted.'
      WRITE (6,*) "xu min & ext:",xumin,xu_ext(1)
      WRITE (6,*) "xu max & ext:",xumax,xu_ext(nx_ext)
      WRITE (6,*) "yv min & ext:",yvmin,yv_ext(1)
      WRITE (6,*) "yv max & ext:",yvmax,yv_ext(ny_ext)
      istatus = -1
      RETURN
    END IF

    !
    ! prepare parameters for horizontal interpolations
    !
    CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,x_ext,y_ext,           &
                 dxfld,dyfld,rdxfld,rdyfld)
    CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,xu_ext,yu_ext,         &
                 dxfldu,dyfldu,rdxfldu,rdyfldu)
    CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,xv_ext,yv_ext,         &
                 dxfldv,dyfldv,rdxfldv,rdyfldv)

  END IF      ! use_wrf_grid /= 1

!
!-----------------------------------------------------------------------
!
!  Test code, for diagnostic testing.
!  Find x,y of Norman sounding in external grid.
!
!-----------------------------------------------------------------------

  CALL output_diagnose_ext( nx_ext,ny_ext,nz_ext,x_ext,y_ext,latdiag,londiag,   &
                            hgt_ext,p_ext,t_ext,pt_ext,qv_ext,u_ext,v_ext,      &
                            zp_ext, trn_ext, znw, istatus )

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
      iextmx=MAX(iextmx,iscl(i,j))
      jextmn=MIN(jextmn,jscl(i,j))
      jextmx=MAX(jextmx,jscl(i,j))

      iextmn=MIN(iextmn,iu(i,j))
      iextmx=MAX(iextmx,iu(i,j))
      jextmn=MIN(jextmn,ju(i,j))
      jextmx=MAX(jextmx,ju(i,j))

      iextmn=MIN(iextmn,iv(i,j))
      iextmx=MAX(iextmx,iv(i,j))
      jextmn=MIN(jextmn,jv(i,j))
      jextmx=MAX(jextmx,jv(i,j))
    END DO
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

  RETURN
END SUBROUTINE prepare_interpolation

!#######################################################################

SUBROUTINE output_diagnose_arps( nx,ny,nz,xscl,yscl,latdiag,londiag,    &
                                 zps,pbar,pprt,ptbar,ptprt,qv,u,v,      &
                                 istatus)
!
!-----------------------------------------------------------------------
!
!  Test code, for diagnostic testing.
!  Find x,y of diagnostic sounding location in ARPS grid.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny, nz
  REAL,    INTENT(IN)  :: xscl(nx), yscl(ny)
  REAL,    INTENT(IN)  :: latdiag, londiag
  REAL,    INTENT(IN)  :: zps(nx,ny,nz)
  REAL,    INTENT(IN)  :: pbar(nx,ny,nz), ptbar(nx,ny,nz)
  REAL,    INTENT(IN)  :: pprt(nx,ny,nz), ptprt(nx,ny,nz), qv(nx,ny,nz)
  REAL,    INTENT(IN)  :: u(nx,ny,nz), v(nx,ny,nz)

  INTEGER, INTENT(OUT) ::  istatus
!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: idiag, jdiag
  REAL    :: xdiag, ydiag
  REAL    :: latd, lond
  REAL    :: dmin,dd
  REAL    :: ppasc,pmb,tc,tdc,theta,smix,e,bige,alge,dir,spd

  INCLUDE 'phycst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL lltoxy(1,1,latdiag,londiag,xdiag,ydiag)

  IF (xdiag > xscl(1) .AND. xdiag < xscl(nx) .AND.                      &
      ydiag > xscl(1) .AND. ydiag < xscl(ny) ) THEN

    dmin=((xdiag-xscl(1))*(xdiag-xscl(1))+                              &
          (ydiag-yscl(1))*(ydiag-yscl(1)))
    idiag=1
    jdiag=1

    DO j=2,ny-2
      DO i=2,nx-2
        dd=((xdiag-xscl(i))*(xdiag-xscl(i))+                            &
            (ydiag-yscl(j))*(ydiag-yscl(j)))
        IF(dd < dmin) THEN
          dmin=dd
          idiag=i
          jdiag=j
        END IF
      END DO
    END DO
    CALL xytoll(1,1,xscl(idiag),yscl(jdiag), latd,lond)
    WRITE(6,'(a,f10.4,f10.4,/a,i5,i5,a,f10.4,f10.4)')                   &
          ' Nearest ARPS pt to diagnostic lat,lon: ',                   &
          latdiag,londiag,                                              &
          ' Diagnostic i,j: ',                                          &
          idiag,jdiag,' lat,lon= ',latd,lond
    WRITE(6,'(///a,/2x,a)')                                             &
        ' ARPS extracted sounding at idiag,jdiag',                      &
        'k   pres   hgt   temp   theta   dewp     u     v     dir    spd'
!
!-----------------------------------------------------------------------
!
!  Convert units of ARPS data and write as a sounding.
!
!-----------------------------------------------------------------------
!
    DO k=nz-2,1,-1
      ppasc=pbar(idiag,jdiag,k)+pprt(idiag,jdiag,k)
      pmb=.01*(pbar(idiag,jdiag,k)+pprt(idiag,jdiag,k))
      theta=ptbar(idiag,jdiag,k)+ptprt(idiag,jdiag,k)
      tc=(theta*((ppasc/p0)**rddcp))-273.15
      IF( qv(idiag,jdiag,k) > 0.) THEN
        smix=qv(idiag,jdiag,k)/(1.-qv(idiag,jdiag,k))
        e=(pmb*smix)/(0.62197 + smix)
        bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
        alge = ALOG(bige/6.112)
        tdc = (alge * 243.5) / (17.67 - alge)
      ELSE
        tdc = tc-30.
      END IF

      CALL uvrotdd(1,1,londiag,u(idiag,jdiag,k),v(idiag,jdiag,k),dir,spd)

      WRITE(6,'(i4,f6.0,f9.0,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1)')      &
              k,pmb,                                                    &
              zps(idiag,jdiag,k),                                       &
              tc,theta,tdc,                                             &
              u(idiag,jdiag,k),                                         &
              v(idiag,jdiag,k),                                         &
              dir,spd
    END DO

  END IF

  RETURN
END SUBROUTINE output_diagnose_arps

!#######################################################################

SUBROUTINE output_diagnose_ext( nx_ext,ny_ext,nz_ext,x_ext,y_ext,latdiag,londiag,   &
                                hgt_ext,p_ext,t_ext,pt_ext,qv_ext,u_ext,v_ext,      &
                                zp_ext, trn_ext,znw, istatus)
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Test code, for diagnostic testing.
!  Find x,y of Norman sounding in external grid.
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN)  :: nx_ext, ny_ext, nz_ext
  REAL,    INTENT(IN)  :: x_ext(nx_ext), y_ext(ny_ext)
  REAL,    INTENT(IN)  :: latdiag, londiag
  REAL,    INTENT(IN)  :: hgt_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: p_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: t_ext(nx_ext,ny_ext,nz_ext), pt_ext(nx_ext,ny_ext,nz_ext), qv_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: u_ext(nx_ext,ny_ext,nz_ext), v_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: zp_ext(nx_ext,ny_ext,nz_ext), trn_ext(nx_ext,ny_ext), znw(nz_ext)

  INTEGER, INTENT(OUT) ::  istatus
!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: idiag, jdiag
  REAL    :: xdiag, ydiag
  REAL    :: latd, lond
  REAL    :: dmin,dd
  REAL    :: pmb,tc,tdc,theta,smix,e,bige,alge,dir,spd

  INCLUDE 'phycst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL lltoxy(1,1,latdiag,londiag,xdiag,ydiag)

  IF (xdiag > x_ext(1) .AND. xdiag < x_ext(nx_ext) .AND.                &
      ydiag > y_ext(1) .AND. ydiag < y_ext(ny_ext) ) THEN

    dmin=((xdiag-x_ext(1))*(xdiag-x_ext(1))+                            &
          (ydiag-y_ext(1))*(ydiag-y_ext(1)))

    idiag=1
    jdiag=1

    DO j=1,ny_ext-1
      DO i=1,nx_ext-1
        dd=((xdiag-x_ext(i))*(xdiag-x_ext(i))+                          &
            (ydiag-y_ext(j))*(ydiag-y_ext(j)))
        IF(dd < dmin) THEN
          dmin=dd
          idiag=i
          jdiag=j
        END IF
      END DO
    END DO
    CALL xytoll(1,1,x_ext(idiag),y_ext(jdiag), latd,lond)
    WRITE(6,'(a,f10.4,f10.4,/a,i5,i5,a,f10.4,f10.4)')                   &
        ' Nearest ext pt to diagnostic lat,lon: ', latdiag,londiag,     &
        ' Diagnostic i,j: ', idiag,jdiag,                               &
        ' lat,lon= ', latd,lond
    WRITE(6,'(///a,/2x,a)') ' External sounding at idiag,jdiag',        &
      'k   pres    hgt    temp   theta   dewp     u      v     dir    spd'

    !
    !  Convert units of external data and write as a sounding.
    !
    !
    DO k=nz_ext,1,-1
      pmb=.01*p_ext(idiag,jdiag,k)
      tc=t_ext(idiag,jdiag,k)-273.15
      theta=pt_ext(idiag,jdiag,k)
      IF( qv_ext(idiag,jdiag,k) > 0.) THEN
        smix=qv_ext(idiag,jdiag,k)/(1.-qv_ext(idiag,jdiag,k))
        e=(pmb*smix)/(0.62197 + smix)
        bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
        alge = ALOG(bige/6.112)
        tdc = (alge * 243.5) / (17.67 - alge)
      ELSE
        tdc = tc-30.
      END IF

      CALL uvrotdd(1,1,lond,u_ext(idiag,jdiag,k),v_ext(idiag,jdiag,k),  &
                   dir,spd)

      WRITE(6,'(i4,f6.0,f9.0,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1)')      &
                    k,pmb,                                              &
                    hgt_ext(idiag,jdiag,k),                             &
                    tc,theta,tdc,                                       &
                    u_ext(idiag,jdiag,k),                               &
                    v_ext(idiag,jdiag,k),                               &
                    dir,spd
    END DO

    !IF (lvldbg > 0 ) THEN
      WRITE(*,'(1x,a)') '============================'
      WRITE(*,'(1x,a,2(I0,a))') 'WRF vertical levels at (',idiag,',',jdiag,').'
      WRITE(*,'(1x,a)') '  No.      ZNW         Z(m)      DZ(m)  '
      WRITE(*,'(1x,a)') '-----  ----------  ----------  ----------'
      DO k = 1,nz_ext-1
        WRITE(*,'(1x,I5,2x,F10.3,2(2x,F10.2))') k, znw(k),          &
           zp_ext(idiag,jdiag,k)-trn_ext(idiag,jdiag),              &
           zp_ext(idiag,jdiag,k+1)-zp_ext(idiag,jdiag,k)
      END DO
      k = nz_ext
      WRITE(*,'(1x,I5,2x,F10.3,2x,F10.2)') k, znw(k),               &
           zp_ext(idiag,jdiag,k)-trn_ext(idiag,jdiag)
    !END IF

  END IF       ! do diagnose

  RETURN
END SUBROUTINE output_diagnose_ext

!#######################################################################

SUBROUTINE print_maxmin_ext( myproc, nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyp_ext, &
                         lat_ext,lon_ext,p_ext,hgt_ext,t_ext,           &
                         u_ext,v_ext,w_ext,qv_ext,qscalar_ext,          &
                         tsoil_ext,qsoil_ext,wetcanp_ext,snowdpth_ext,  &
                         tke_ext,trn_ext,veg_ext,istatus )
  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER, INTENT(IN)  :: myproc
  INTEGER, INTENT(IN)  :: nx_ext,ny_ext,nz_ext, nzsoil_ext, nstyp_ext
  REAL,    INTENT(IN)  :: lat_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: lon_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: p_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: hgt_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: t_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: u_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: v_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: w_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: qv_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: qscalar_ext(nx_ext,ny_ext,nz_ext,nscalar)
  REAL,    INTENT(IN)  :: tsoil_ext(nx_ext,ny_ext,nzsoil_ext,0:nstyp_ext)
  REAL,    INTENT(IN)  :: qsoil_ext(nx_ext,ny_ext,nzsoil_ext,0:nstyp_ext)
  REAL,    INTENT(IN)  :: wetcanp_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: snowdpth_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: tke_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(IN)  :: trn_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: veg_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  REAL    :: amax, amin
  INTEGER :: k, nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL a3dmax0(lat_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,1,1,1,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'lat_ext_min= ', amin,', lat_ext_max=',amax
  CALL a3dmax0(lon_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,1,1,1,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'lon_ext_min= ', amin,', lon_ext_max=',amax

  CALL a3dmax0(p_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'p_ext_min  = ', amin,', p_ext_max  =',amax
  CALL a3dmax0(hgt_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'hgt_ext_min= ', amin,', hgt_ext_max=',amax
  CALL a3dmax0(t_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               't_ext_min  = ', amin,', t_ext_max  =',amax
  CALL a3dmax0(u_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'u_ext_min  = ', amin,', u_ext_max  =',amax
  CALL a3dmax0(v_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'v_ext_min  = ', amin,', v_ext_max  =',amax
  CALL a3dmax0(w_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'w_ext_min  = ', amin,', w_ext_max  =',amax
  CALL a3dmax0(qv_ext ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'qv_ext_min = ', amin,', qv_ext_max =',amax
  DO nq=1,nscalar
    CALL a3dmax0(qscalar_ext(1,1,1,nq),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,   &
               1,nz_ext,1,nz_ext,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,a,e13.6,2a,e13.6)')                   &
          TRIM(qnames(nq))//'_ext_min = ', amin,', ',TRIM(qnames(nq))//'_ext_max =',amax
  END DO

  DO k = 1, nzsoil_ext
    CALL a3dmax0(tsoil_ext(1,1,k,0),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,&
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')                   &
        'tsoil_',k,'_ext_min= ', amin,', tsoil_',k,'_ext_max=',amax
    CALL a3dmax0(qsoil_ext(1,1,k,0),1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,&
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')                   &
        'qsoil_',k,'_ext_min= ', amin,', qsoil_',k,'_ext_max=',amax
  END DO

  CALL a3dmax0(wetcanp_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
               1,1,1,1,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'wetcanp_ext_min = ', amin,', wetcanp_ext_max = ',amax
  CALL a3dmax0(snowdpth_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,        &
               1,1,1,1,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'snowd_ext_min   = ', amin,', snow_ext_max    = ',amax

  CALL a3dmax0(tke_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,nz_ext,1,nz_ext,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'tke_ext_min = ', amin,', tke_ext_max = ',amax

  CALL a3dmax0(trn_ext    ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
               1,1,1,1,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'trn_ext_min     = ', amin,', trn_ext_max     = ',amax

  CALL a3dmax0(veg_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,             &
               1,1,1,1,amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
               'veg_ext_min     = ', amin,', veg_ext_max     = ',amax

  IF(myproc == 0) PRINT*,' '

  RETURN
END SUBROUTINE print_maxmin_ext

!#######################################################################

SUBROUTINE print_maxmin_arps(myproc,nx,ny,nz,nzsoil,mp_physics,         &
                             x, y, z, zp, rhobar, u, v, w,              &
                             ptbar, ptprt, pbar, pprt, qv, qscalar,     &
                             tke, tsoil,qsoil,wetcanp,istatus)

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER, INTENT(IN)  :: myproc
  INTEGER, INTENT(IN)  :: nx,ny,nz, nzsoil
  INTEGER, INTENT(IN)  :: mp_physics

  REAL,    INTENT(IN) :: x (nx)
  REAL,    INTENT(IN) :: y (ny)
  REAL,    INTENT(IN) :: z (nz)
  REAL,    INTENT(IN) :: zp (nx,ny,nz)
  REAL,    INTENT(IN) :: rhobar (nx,ny,nz)
  REAL,    INTENT(IN) :: u (nx,ny,nz)
  REAL,    INTENT(IN) :: v (nx,ny,nz)
  REAL,    INTENT(IN) :: w (nx,ny,nz)
  REAL,    INTENT(IN) :: ptbar (nx,ny,nz)
  REAL,    INTENT(IN) :: ptprt (nx,ny,nz)
  REAL,    INTENT(IN) :: pbar (nx,ny,nz)
  REAL,    INTENT(IN) :: pprt (nx,ny,nz)
  REAL,    INTENT(IN) :: qv (nx,ny,nz)
  REAL,    INTENT(IN) :: qscalar (nx,ny,nz, nscalar)
  REAL,    INTENT(IN) :: tke (nx,ny,nz)
  REAL,    INTENT(IN) :: tsoil (nx,ny,nzsoil,0:nstyp)
  REAL,    INTENT(IN) :: qsoil (nx,ny,nzsoil,0:nstyp)
  REAL,    INTENT(IN) :: wetcanp (nx,ny,0:nstyp)

  INTEGER, INTENT(OUT) :: istatus
!-----------------------------------------------------------------------

  REAL    :: amin, amax
  INTEGER :: k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF(myproc == 0) WRITE(6,'(/1x,a/)')                                   &
      'Min and max of External data interpolated to the ARPS grid:'

  CALL a3dmax0(x,1,nx,1,nx,1,1,1,1,1,1,1,1, amax,amin)
  IF(myproc == 0) WRITE(6,'(/1x,2(a,e13.6))')                           &
          'xmin    = ', amin,',  xmax    =',amax

  CALL a3dmax0(y,1,ny,1,ny,1,1,1,1,1,1,1,1, amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'ymin    = ', amin,',  ymax    =',amax

  CALL a3dmax0(z,1,nz,1,nz,1,1,1,1,1,1,1,1, amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'zmin    = ', amin,',  zmax    =',amax

  CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,                    &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'zpmin   = ', amin,', zpmax    =',amax

  CALL a3dmax0(rhobar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,              &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'rhobarmin=', amin,', rhobarmax=',amax

  CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,                     &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'umin    = ', amin,',  umax    =',amax

  CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,                     &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'vmin    = ', amin,',  vmax    =',amax

  CALL a3dmax0(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,                     &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'wmin    = ', amin,',  wmax    =',amax

  CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,               &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'ptbarmin= ', amin,',  ptbarmax=',amax

  CALL a3dmax0(ptprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,               &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'ptprtmin= ', amin,',  ptprtmax=',amax

  CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'pbarmin= ', amin,',  pbarmax =',amax

  CALL a3dmax0(pprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'pprtmin = ', amin,',  pprtmax =',amax

  CALL a3dmax0(qv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                  &
            amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'qvmin   = ', amin,',  qvmax   =',amax

  !
  ! Rime - special for mp_physics = 5
  !
  IF( mp_physics == 5 ) THEN
    CALL a3dmax0(qscalar(:,:,:,nscalar-1),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, &
              amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
            'RimeFmin   = ', amin,',  RimeFmax   =',amax
  END IF

  !
  ! Reflectivity special for a few of schemes
  !
  SELECT CASE (mp_physics)

    CASE(2,4,5,6,8,9,10,14,16,17)
      CALL a3dmax0(qscalar(:,:,:,nscalar),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, &
                   amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'reflmin   = ', amin,',  reflmax   =',amax
  END SELECT

  CALL a3dmax0(tke,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                 &
               amax,amin)
  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                            &
          'tkemin   = ', amin,',  tkemax   = ',amax

  IF(sfcout == 1) THEN

    IF (soilmodel_option == 1) THEN

      CALL a3dmax0(tsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil,     &
                amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'tsoilmin = ', amin,',  tsoilmax =',amax
      CALL a3dmax0(qsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil,     &
                amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'qsoilmin= ', amin,',  qsoilmax=',amax
      CALL a3dmax0(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'wetcanpmin = ', amin,',  wetcanpmax =',amax

    ELSE IF (soilmodel_option == 2) THEN

      DO k=1,nzsoil
        CALL a3dmax0(tsoil(1,1,k,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,      &
                 1,amax,amin)
        IF(myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')               &
               'tsoil_',k,'_min= ', amin,', tsoil_',k,'_max=',amax
        CALL a3dmax0(qsoil(1,1,k,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,      &
                 1,amax,amin)
        IF(myproc == 0) WRITE(6,'(1x,2(a,I3.3,a,e13.6))')               &
               'qsoil_',k,'_min = ', amin,',  qsoil_',k,'_max =',amax
      END DO
      CALL a3dmax0(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
              'wetcanpmin = ', amin,',  wetcanpmax =',amax
    END IF

  END IF

  RETURN
END SUBROUTINE print_maxmin_arps
