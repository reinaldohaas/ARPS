  PROGRAM difobs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    PROGRAM DIFOBS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  PURPOSE:
!  Calculate difference of ARPS model grid from observations.
!  Report statistics, including bias and RMS by observation source.
!
!  AUTHOR: Keith Brewster and Nazir Said, CAPS
!  Completed 08/24/03
!
!  MODIFICATIONS
!  10/14/2003 (Keith Brewster)
!  Added parsing of iuse lists to screen data.
!
!  7/27/2007 (Keith Brewster)
!  Modified to be compatible with the MPI version of ADAS subroutines.
!
!  12/6/2012 (Keith Brewster)
!  Updated to add MPI and split file reading capability.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'alloc.inc'
  INCLUDE 'adas.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'

  REAL :: exbcbuf( 1 ) ! EXBC buffer array (unused)

!
!-----------------------------------------------------------------------
!
!  Arrays defining model grid
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)      ! The x-coord. of the physical and
                                      ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)      ! The y-coord. of the physical and
                                      ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)      ! The z-coord. of the computational grid.
                                      ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:)  ! The physical height coordinate defined at
                                      ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! Soil level depth.

  REAL, ALLOCATABLE :: hterain(:,:)   ! The height of the terrain.
  REAL, ALLOCATABLE :: mapfct(:,:,:)  ! Map factors at scalar, u and v points

  REAL, ALLOCATABLE :: j1    (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as - d( zp )/d( x ).
  REAL, ALLOCATABLE :: j2    (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as - d( zp )/d( y ).
  REAL, ALLOCATABLE :: j3    (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ).
  REAL, ALLOCATABLE :: j3soil(:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zpsoil )/d( z ).
  REAL, ALLOCATABLE :: aj3x  (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL, ALLOCATABLE :: aj3y  (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL, ALLOCATABLE :: aj3z  (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL, ALLOCATABLE :: j3inv (:,:,:)  ! Inverse of j3
  REAL, ALLOCATABLE :: j3soilinv (:,:,:)! Inverse of j3soil
!
!-----------------------------------------------------------------------
!
!  ARPS Time-dependent variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: u     (:,:,:)  ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v     (:,:,:)  ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w     (:,:,:)  ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: wcont (:,:,:)  ! Contravariant vertical velocity (m/s)
  REAL, ALLOCATABLE :: ptprt (:,:,:)  ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt  (:,:,:)  ! Perturbation pressure (Pascal)

  REAL, ALLOCATABLE :: qv    (:,:,:)  ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)
  REAL, ALLOCATABLE :: tke   (:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL, ALLOCATABLE :: ubar  (:,:,:)  ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar  (:,:,:)  ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar (:,:,:)  ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar  (:,:,:)  ! Base state pressure (Pascal).
  REAL, ALLOCATABLE :: rhostr(:,:,:)  ! Base state density rhobar times j3.
  REAL, ALLOCATABLE :: qvbar (:,:,:)  ! Base state water vapor specific
                                      ! humidity(kg/kg)

  REAL, ALLOCATABLE :: udteb (:,:)    ! T-tendency of u at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtwb (:,:)    ! T-tendency of u at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtnb (:,:)    ! T-tendency of v at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtsb (:,:)    ! T-tendency of v at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: pdteb (:,:)    ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: pdtwb (:,:)    ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: pdtnb (:,:)    ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: pdtsb (:,:)    ! T-tendency of pprt at s-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: trigs1(:)      ! Array containing pre-computed trig
                                      ! function for use in fft.
  REAL, ALLOCATABLE :: trigs2(:)      ! Array containing pre-computed trig
                                      ! function for use in fft.
  INTEGER, ALLOCATABLE :: ifax1(:)    ! Array containing the factors of nx.
  INTEGER, ALLOCATABLE :: ifax2(:)    ! Array containing the factors of ny.
  REAL, ALLOCATABLE :: vwork1 (:,:)   ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: vwork2 (:,:)   ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: wsave1 (:)     ! Work array for fftopt=2.
  REAL, ALLOCATABLE :: wsave2 (:)     ! Work array for fftopt=2.

  REAL, ALLOCATABLE :: ppi(:,:,:)     ! Exner function.
  REAL, ALLOCATABLE :: csndsq(:,:,:)  ! Speed of sound squared

  REAL, ALLOCATABLE :: radfrc(:,:,:)  ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw(:,:)     ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx(:,:)     ! Net absorbed radiation by the surface
  REAL, ALLOCATABLE :: radswnet(:,:)  ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)   ! Incominging longwave radiation
!
!-----------------------------------------------------------------------
!
!  ARPS Surface variables:
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL, ALLOCATABLE ::    stypfrct(:,:,:) ! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp (:,:)    ! Vegetation type
  REAL, ALLOCATABLE ::    lai    (:,:)    ! Leaf Area Index
  REAL, ALLOCATABLE ::    roufns (:,:)    ! Surface roughness
  REAL, ALLOCATABLE ::    veg    (:,:)    ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil(:,:,:,:) ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil(:,:,:,:) ! Soil moisture
  REAL, ALLOCATABLE :: qvsfc(:,:,:)   ! Effective qv at sfc.

  REAL, ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (:)

  REAL, ALLOCATABLE :: raing(:,:)     ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)     ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:) ! precipitation rate (kg/(m**2*s))
                                      ! prcrate(:,:,:) = total precipitation rate
                                      ! prcrate(:,:,:) = grid scale precip. rate
                                      ! prcrate(:,:,:) = cumulus precip. rate
                                      ! prcrate(:,:,:) = microphysics precip. rate

  REAL, ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
! Statistical counters
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: knt(:)
  INTEGER, ALLOCATABLE :: kntsngt(:)
  INTEGER, ALLOCATABLE :: kntuat(:)
  INTEGER, ALLOCATABLE :: kntrett(:)
  INTEGER, ALLOCATABLE :: kntradt(:)
  INTEGER, ALLOCATABLE :: kntsng(:,:)
  INTEGER, ALLOCATABLE :: kntua(:,:)
  INTEGER, ALLOCATABLE :: kntret(:,:)
  INTEGER, ALLOCATABLE :: kntrad(:,:)

  REAL, ALLOCATABLE :: bias(:)
  REAL, ALLOCATABLE :: biassngt(:)
  REAL, ALLOCATABLE :: biasuat(:)
  REAL, ALLOCATABLE :: biasrett(:)
  REAL, ALLOCATABLE :: biasradt(:)
  REAL, ALLOCATABLE :: biassng(:,:)
  REAL, ALLOCATABLE :: biasua(:,:)
  REAL, ALLOCATABLE :: biasret(:,:)
  REAL, ALLOCATABLE :: biasrad(:,:)

  REAL, ALLOCATABLE :: rms(:)
  REAL, ALLOCATABLE :: rmssngt(:)
  REAL, ALLOCATABLE :: rmsuat(:)
  REAL, ALLOCATABLE :: rmsrett(:)
  REAL, ALLOCATABLE :: rmsradt(:)
  REAL, ALLOCATABLE :: rmssng(:,:)
  REAL, ALLOCATABLE :: rmsua(:,:)
  REAL, ALLOCATABLE :: rmsret(:,:)
  REAL, ALLOCATABLE :: rmsrad(:,:)

  REAL, ALLOCATABLE :: rngsqi(:)
  REAL, ALLOCATABLE :: wgtsum(:,:)
  REAL, ALLOCATABLE :: zsum(:,:)
  REAL, ALLOCATABLE :: sqsum(:,:)
!
!-----------------------------------------------------------------------
!
!  Analysis variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: anx(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Cross-category correlations
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: xcor(:,:)
  INTEGER, ALLOCATABLE :: icatg(:,:)
!
!-----------------------------------------------------------------------
!
!  Additional grid variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zs(:,:,:)
  REAL, ALLOCATABLE :: xslg(:)
  REAL, ALLOCATABLE :: yslg(:)
  REAL, ALLOCATABLE :: latgr(:,:)
  REAL, ALLOCATABLE :: longr(:,:)
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem2  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem3  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem4  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem5  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem6  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem7  (:,:,:)     ! Temporary work array.

  INTEGER, ALLOCATABLE :: item1 (:)

  INTEGER, ALLOCATABLE :: tems1di(:)
  INTEGER, ALLOCATABLE :: temr1di(:)
  REAL, ALLOCATABLE :: tems1dr(:)
  REAL, ALLOCATABLE :: temr1dr(:)

  REAL, ALLOCATABLE :: tems2dr(:,:)
  REAL, ALLOCATABLE :: temr2dr(:,:)

  REAL, ALLOCATABLE :: tems2drua(:,:)
  REAL, ALLOCATABLE :: temr2drua(:,:)

  REAL, ALLOCATABLE :: tems3dr(:,:,:)
  REAL, ALLOCATABLE :: temr3dr(:,:,:)

  REAL, ALLOCATABLE :: tem1soil (:,:,:)  ! Temporary work array.
  REAL, ALLOCATABLE :: tem2soil (:,:,:)  ! Temporary work array.
  REAL, ALLOCATABLE :: tem3soil (:,:,:)  ! Temporary work array.
  REAL, ALLOCATABLE :: tem4soil (:,:,:)  ! Temporary work array.
  REAL, ALLOCATABLE :: tem5soil (:,:,:)  ! Temporary work array.

  REAL, ALLOCATABLE :: tem1d (:)  ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  ARPS dimensions:
!
!  nx, ny, nz: Dimensions of analysis grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
  INTEGER :: nzsoil   ! Number of soil levels

  INTEGER :: nxlg     ! Number of grid points in the x-direction large grid
  INTEGER :: nylg     ! Number of grid points in the y-direction large grid

!
!-----------------------------------------------------------------------
!
!  Soil types.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nstyps            ! Number of soil types
!
!-----------------------------------------------------------------------
!
!  Misc
!
!-----------------------------------------------------------------------
!
  INTEGER :: nt                ! Number of time levels of data
  INTEGER :: ncat              ! Number of correlation categories
!
!-----------------------------------------------------------------------
!
!  Analysis parameters
!
!-----------------------------------------------------------------------
!
  REAL :: rhmin
  INTEGER :: iqspr,klim
  PARAMETER (rhmin=0.05,  & ! rh safety net value to prevent neg qv
             iqspr=3,     & ! Use qobs of pstn to combine x,y,elev
             klim= 1)       ! Min of one other sfc station for QC
!
!-----------------------------------------------------------------------
!
!  Indices of specific observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: ipres=3   ! Pressure
  INTEGER, PARAMETER :: iptmp=4   ! Potential Temperature
  INTEGER, PARAMETER :: iqv=5     ! Specific Humidity
!
!-----------------------------------------------------------------------
!
!  ARPS include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'adjust.inc'

  INTEGER, PARAMETER :: nvar_anx = nvar_anx_adas
  INTEGER, PARAMETER :: nvar_rad = nvar_rad_adas
!
!-----------------------------------------------------------------------
!
!  Surface Station variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: isrcsng(mx_sng)
  INTEGER :: icatsng(mx_sng)
  INTEGER :: musesng(0:nsrc_sng)
  REAL :: latsng(mx_sng,ntime)
  REAL :: lonsng(mx_sng,ntime)
  REAL :: hgtsng(mx_sng,ntime)
  REAL :: xsng(mx_sng)
  REAL :: ysng(mx_sng)
  REAL :: trnsng(mx_sng)
  INTEGER :: timesng(mx_sng,ntime)
  CHARACTER (LEN=5) :: stnsng(mx_sng,ntime)
!
!-----------------------------------------------------------------------
!
!  Surface (single-level) read-in observation variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: wx(mx_sng,ntime)
  CHARACTER (LEN=8) :: csrcsng(mx_sng,ntime)
  CHARACTER (LEN=1) :: store_emv(mx_sng,5,ntime)
  CHARACTER (LEN=4) :: store_amt(mx_sng,5,ntime)
  INTEGER :: kloud(mx_sng,ntime),idp3(mx_sng,ntime)
  REAL :: store_hgt(mx_sng,5,ntime)
  REAL :: obrdsng(mx_sng,nvar_sng,ntime)
  REAL :: obsng(nvar_anx,mx_sng)
  REAL :: odifsng(nvar_anx,mx_sng)
  REAL :: oanxsng(nvar_anx,mx_sng)
  REAL :: thesng(mx_sng)
  INTEGER :: ival(mx_sng)
!
!-----------------------------------------------------------------------
!
!  Upper Air Station variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: isrcua(mx_ua)
  INTEGER :: museua(0:nsrc_ua)
  REAL :: xua(mx_ua)
  REAL :: yua(mx_ua)
  REAL :: trnua(mx_ua)
  REAL :: elevua(mx_ua)
  REAL :: hgtua(nz_ua,mx_ua)
  INTEGER :: nlevsua(mx_ua)
  CHARACTER (LEN=5) :: stnua(mx_ua)
!
!-----------------------------------------------------------------------
!
!  Upper-air observation variables
!
!-----------------------------------------------------------------------
!
  REAL :: obsua(nvar_anx,nz_ua,mx_ua)
  REAL :: odifua(nvar_anx,nz_ua,mx_ua)
  REAL :: oanxua(nvar_anx,nz_ua,mx_ua)
  REAL :: theua(nz_ua,mx_ua)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnrad(mx_rad)
  INTEGER :: isrcrad(0:mx_rad)
  INTEGER :: muserad(0:nsrc_rad)
  REAL :: latrad(mx_rad),lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
  REAL :: refelvmin(mx_rad)
  REAL :: refelvmax(mx_rad)
  REAL :: refrngmin(mx_rad)
  REAL :: refrngmax(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: irad(mx_colrad)
  INTEGER :: nlevrad(mx_colrad)
  REAL :: latradc(mx_colrad)
  REAL :: lonradc(mx_colrad)
  REAL :: xradc(mx_colrad)
  REAL :: yradc(mx_colrad)
  REAL :: trnradc(mx_colrad)
  REAL :: distrad(mx_colrad)
  REAL :: uazmrad(mx_colrad)
  REAL :: vazmrad(mx_colrad)
  REAL :: hgtradc(nz_rdr,mx_colrad)
  REAL :: wradc(nz_rdr)
  REAL :: theradc(nz_rdr,mx_colrad)
  REAL :: obsrad(nvar_radin,nz_rdr,mx_colrad)
  REAL :: odifrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: oanxrad(nvar_rad,nz_rdr,mx_colrad)
  REAL ::    dsdr(nz_rdr,mx_colrad)
  REAL ::    dhdr(nz_rdr,mx_colrad)!
!-----------------------------------------------------------------------
!
!  Retrieval radar variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnret(mx_ret)
  INTEGER :: isrcret(0:mx_ret)
  INTEGER :: museret(0:nsrc_ret)
  REAL :: latret(mx_ret),lonret(mx_ret)
  REAL :: elvret(mx_ret)
!
!-----------------------------------------------------------------------
!
!  Retrieval observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iret(mx_colret)
  INTEGER :: nlevret(mx_colret)
  REAL :: latretc(mx_colret),lonretc(mx_colret)
  REAL :: xretc(mx_colret),yretc(mx_colret)
  REAL :: hgtretc(nz_ret,mx_colret)
  REAL :: theretc(nz_ret,mx_colret)
  REAL :: obsret(nvar_anx,nz_ret,mx_colret)
  REAL :: qrret(nz_ret,mx_colret)
  REAL :: odifret(nvar_anx,nz_ret,mx_colret)
  REAL :: oanxret(nvar_anx,nz_ret,mx_colret)
!
!-----------------------------------------------------------------------
!
!  Namen de variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: var_sfc(nvar_sng)
  DATA var_sfc/ 't     ','td    ','dd    ','ff    ',                    &
                'urot  ','vrot  ','pstn  ','pmsl  ',                    &
                'alt   ','ceil  ','low   ','cover ',                    &
                'solar ','vis'/
  CHARACTER (LEN=6) :: var_ua(nvar_ua)
  DATA var_ua / 'press ',                                               &
                'temp  ','dewpt ','dd    ','ff    '/
  CHARACTER (LEN=6) :: var_anx(nvar_anx)
  DATA var_anx/ 'u     ','v     ',                                      &
                'press ','theta ','qv    '/
  CHARACTER (LEN=6) :: var_rad(3)
  DATA var_rad/ 'radvel','reflec','cmpref'/
!
!-----------------------------------------------------------------------
!
!  Source-dependent parameters
!  Qsrc is the expected square error it is used for
!  setting the QC threshold (qclim) and for relative
!  weighting of observations.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: srcback
  REAL :: hqback(nz_tab)
  INTEGER :: nlvqback
  REAL :: qback(nvar_anx,nz_tab)

  CHARACTER (LEN=80) :: srcsng_full(nsrc_sng)
  REAL :: qsrcsng(nvar_anx,nsrc_sng)

  CHARACTER (LEN=80) :: srcua_full(nsrc_ua)
  INTEGER :: nlvuatab(nsrc_ua)
  REAL :: huaqsrc(nz_tab,nsrc_ua)
  REAL :: qsrcua(nvar_anx,nz_tab,nsrc_ua)

  CHARACTER (LEN=80) :: srcrad_full(nsrc_rad)
  REAL :: qsrcrad(nvar_radin,nsrc_rad)

  CHARACTER (LEN=80) :: srcret_full(nsrc_ret)
  INTEGER :: nlvrttab(nsrc_ret)
  REAL :: hrtqsrc(nz_tab,nsrc_ret)
  REAL :: qsrcret(nvar_anx,nz_tab,nsrc_ret)
!
!-----------------------------------------------------------------------
!
!  Quality Control Parameters and Variables
!
!-----------------------------------------------------------------------
!
  REAL :: rmiss
  PARAMETER (rmiss=-99.)
!
  INTEGER :: qualrdsng(mx_sng,nvar_sng,ntime)
  INTEGER :: qualsng(nvar_anx,mx_sng)
  REAL :: qobsng(nvar_anx,mx_sng)
  REAL :: qcmulsng(nsrc_sng)
  REAL :: qcthrsng(nvar_anx,nsrc_sng)
  REAL :: barqclim(nvar_anx,nsrc_sng)
!
  INTEGER :: qualua(nvar_anx,nz_ua,mx_ua)
  REAL :: qobsua(nvar_anx,nz_ua,mx_ua)
  REAL :: qcmulua(nsrc_ua)
  REAL :: qcthrua(nvar_anx,nsrc_ua)
!
  INTEGER :: qualrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: qobsrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: qcmulrad(nsrc_rad)
  REAL :: qcthrrad(nvar_radin,nsrc_rad)
!
  INTEGER :: qualret(nvar_anx,nz_ret,mx_colret)
  REAL :: qobsret(nvar_anx,nz_ret,mx_colret)
  REAL :: qcmulret(nsrc_ret)
  REAL :: qcthrret(nvar_anx,nsrc_ret)
!
!-----------------------------------------------------------------------
!
!  Climin is the minimum possible value of each observed variable.
!  Climax is the maximum possible value of each observed variable.
!  Used for surface data only.
!
!-----------------------------------------------------------------------
!
  REAL :: climin(nvar_sng),climax(nvar_sng)
  DATA climin /  -50.,    -50.,    0.,    0.,                           &
                -100.,   -100.,  700.,  880.,                           &
                 880.,      0.,    0.,    0.,                           &
                   0.,      0./
  DATA climax /  125.,    125.,  360.,  100.,                           &
                 100.,    100., 1100., 1090.,                           &
                1090.,  30000.,30000.,    1.,                           &
                1500.,    150./
!
!-----------------------------------------------------------------------
!
!  Filenames and such
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: veriffn
  INTEGER :: nchdif,lveriffn,istat

  CHARACTER (LEN=12) :: suffix
  CHARACTER (LEN=256) :: froot
  CHARACTER (LEN=256) :: qclist
  CHARACTER (LEN=256) :: qcsave
  CHARACTER (LEN=256) :: stats
  INTEGER :: istatus,jstatus
  INTEGER :: nobsng,nobsrd(ntime)
  INTEGER :: i4timanx
  INTEGER :: lenfnm,maxsuf,lensuf,dotloc,mxroot,lenroot
  INTEGER :: iqclist,iqcsave,iwstat

  CHARACTER (LEN=256) :: status_file
  INTEGER n_used_sng(nsrc_sng), n_used_ua(nsrc_sng)
  INTEGER jsta, klev, ngoodlev
!
!-----------------------------------------------------------------------
!
!  MPI bookkeeping
!
!-----------------------------------------------------------------------
!
  INTEGER :: indexsng(mx_sng)    ! indexes each data point according to
                                 ! which local proc domain it resides
  LOGICAL :: usesng(mx_sng)

  INTEGER, ALLOCATABLE :: ksng(:) ! Number of obs owned by each processor
  INTEGER :: ksngmax             ! Max "ksng" value.

  INTEGER :: indexua(mx_ua)      ! indexes each data point according to
                                 ! which local proc domain it resides
  LOGICAL :: useua(mx_ua)
  INTEGER, ALLOCATABLE :: kua(:) ! Number of obs owned by each processor
  INTEGER :: kuamax              ! Max "kua" value.


  INTEGER :: indexrad(mx_rad)    ! indexes each data point according to
                                 ! which local proc domain it resides
  INTEGER :: oindexrad(mx_colrad)

  LOGICAL :: userad(mx_colrad)
  INTEGER :: indexret(mx_ret)    ! indexes each data point according to
                                 ! which local proc domain it resides

  INTEGER :: iproc, iprocv               ! Radii of influence "i" direction
  INTEGER :: jproc, jprocv               ! Radii of influence "j" direction

  INTEGER :: istat_radar, ncolrad_mpi

  INTEGER, ALLOCATABLE :: mpi_map(:,:) ! Map of communication entries
  INTEGER :: nmap       ! Number of entries in "mpi_map".

  INTEGER, ALLOCATABLE :: mpi_mapv(:,:) ! Same for radial velocity
  INTEGER :: nmapv               ! Number of entries in "mpi_mapv".

  INTEGER, PARAMETER :: iradvel = 1   ! Reduce memory requirements for velocity.
                                      ! Never change the set value.  This is a
                                      ! code development tool only!

  REAL, ALLOCATABLE :: ref_mos_3d (:,:,:)
!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: kts2ms=0.514444, mb2pa=100.,f2crat=(5./9.)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!dir$ inline always f_qvsat
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8)  :: rdsource
!  CHARACTER (LEN=80) :: basdmpfn
  CHARACTER (LEN=80) :: header
  INTEGER :: i,j,k,ios,ktab,isrc,jsrc,ivar,ifile,ipass
  INTEGER :: grdbas,lbasdmpf
  INTEGER :: nobsua,ncolrad,ncolret
  INTEGER :: totsng,totua,totret,totrad,totsrc
  REAL :: tgrid,qvmin,qvsat,qsrcmax,rngsq,refmax
  REAL :: temp,tvbar,qvprt,qtot,pconst,pr1,pr2,p0inv
!
  REAL :: score
  REAL, PARAMETER :: rhinf = 1.0
  REAL, PARAMETER :: rvinf = 1.0
  REAL :: wgtvar(nvar_anx)
  DATA wgtvar /1.0, 1.0, 0.1, 1.0, 1000.0/
!
  INTEGER :: ixrad(mx_rad),jyrad(mx_rad)

  CHARACTER(LEN=256) :: namelist_filename

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!  Initializaitons
!-----------------------------------------------------------------------

  usesng(:) = .TRUE.
  useua(:)  = .TRUE.
  userad(:) = .TRUE.

  nt = 1    ! Number of time levels of data
  ncat = 4

!-----------------------------------------------------------------------
!  Set grid mode to non-nesting and grid index, mgrid, to 1.
!-----------------------------------------------------------------------

  myproc=0
  nestgrd=0
  mgrid=1

!-----------------------------------------------------------------------
!  read in ARPS namelists
!-----------------------------------------------------------------------

  print *, ' Calling initpara '

  namelist_filename = ' '
  CALL initpara(nx,ny,nz,nzsoil,nstyps, namelist_filename)

  print *, ' Back from initpara '

  IF (mp_opt > 0) THEN
    nxlg = (nx - 3) * nproc_x + 3
    nylg = (ny - 3) * nproc_y + 3
  ELSE
    nxlg = nx
    nylg = ny
  END IF

  IF (lbcopt == 2) THEN
    WRITE (*,*) "INITADAS: resetting lbcopt to 1 ",                     &
                "& lateral bc's to 4"
    lbcopt = 1
    ebc = 4
    wbc = 4
    nbc = 4
    sbc = 4
  END IF

!-----------------------------------------------------------------------
!  Read in adas namelist parameters
!-----------------------------------------------------------------------
!
  CALL initadas
!
!-----------------------------------------------------------------------
!
!  Allocate adas arrays
!
!-----------------------------------------------------------------------
!
  ALLOCATE(x(nx),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:x")
  ALLOCATE(y(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:y")
  ALLOCATE(z(nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:z")

  ALLOCATE(hterain(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:hterain")
  ALLOCATE(mapfct (nx,ny,8),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:mapfct")

  ALLOCATE(zp  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:zp")
  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:zpsoil")
  ALLOCATE(j1  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:j1")
  ALLOCATE(j2  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:j2")
  ALLOCATE(j3  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:j3")
  ALLOCATE(j3soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:j3soil")
  ALLOCATE(aj3x(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:aj3x")
  ALLOCATE(aj3y(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:aj3y")
  ALLOCATE(aj3z(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:aj3z")
  ALLOCATE(j3inv(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:j3inv")
  ALLOCATE(j3soilinv(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:j3soilinv")

  ALLOCATE(u    (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:u")
  ALLOCATE(v    (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:v")
  ALLOCATE(w    (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:w")
  ALLOCATE(wcont(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:wcont")
  ALLOCATE(ptprt(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ptprt")
  ALLOCATE(pprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:pprt")

  ALLOCATE(qv   (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:qv")
  ALLOCATE(qscalar   (nx,ny,nz,nscalar),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:qscalar")
  ALLOCATE(tke  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tke")

  ALLOCATE(ubar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ubar")
  ALLOCATE(vbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vbar")
  ALLOCATE(ptbar(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ptbar")
  ALLOCATE(pbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:pbar")
  ALLOCATE(rhostr(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:rhostr")
  ALLOCATE(qvbar(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:qvbar")

  ALLOCATE(udteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:udteb")
  ALLOCATE(udtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:udtbw")
  ALLOCATE(vdtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vdtnb")
  ALLOCATE(vdtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vdtsb")

  ALLOCATE(pdteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:pdteb")
  ALLOCATE(pdtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:pdtwb")
  ALLOCATE(pdtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:pdtnb")
  ALLOCATE(pdtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:pdtsb")

  ALLOCATE(trigs1(3*(nx-1)/2+1),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:trigs1")
  ALLOCATE(trigs2(3*(ny-1)/2+1),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:trigs2")
  ALLOCATE(ifax1(13),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ifax1")
  ALLOCATE(ifax2(13),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ifax2")

  ALLOCATE(vwork1(nx+1,ny+1),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vwork1")
  ALLOCATE(vwork2(ny,nx+1),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vwork2")

  ALLOCATE(wsave1(3*(ny-1)+15),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:wsave1")
  ALLOCATE(wsave2(3*(nx-1)+15),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:wsave2")


  ALLOCATE(ppi   (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ppi")
  ALLOCATE(csndsq(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:csndsq")

  ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:radfrc")
  ALLOCATE(radsw (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:radsw")
  ALLOCATE(rnflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:rnflx")
  ALLOCATE(radswnet (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:radswnet")
  ALLOCATE(radlwin (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:radlwin")

  ALLOCATE(soiltyp(nx,ny,nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:soiltyp")
  ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:stypfrct")
  ALLOCATE(vegtyp (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vegtyp")
  ALLOCATE(lai    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:lai")
  ALLOCATE(roufns (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:roufns")
  ALLOCATE(veg    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:veg")

  ALLOCATE(qvsfc  (nx,ny,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:qvsfc")
  ALLOCATE(tsoil  (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tsoil")
  ALLOCATE(qsoil  (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:qsoil")
  ALLOCATE(wetcanp(nx,ny,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:wetcanp")
  ALLOCATE(snowdpth(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:snowdpth")

  ALLOCATE(raing  (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:raing")
  ALLOCATE(rainc  (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:rainc")
  ALLOCATE(prcrate(nx,ny,4),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:prcrate")

  ALLOCATE(usflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:usflx")
  ALLOCATE(vsflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:vsflx")
  ALLOCATE(ptsflx(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ptsflx")
  ALLOCATE(qvsflx(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:qvsflx")

  ALLOCATE(knt(nvar_anx),STAT=istatus)
  ALLOCATE(kntsngt(nvar_anx),STAT=istatus)
  ALLOCATE(kntuat(nvar_anx),STAT=istatus)
  ALLOCATE(kntrett(nvar_anx),STAT=istatus)
  ALLOCATE(kntradt(nvar_anx),STAT=istatus)
  ALLOCATE(kntsng(nvar_anx,nsrc_sng),STAT=istatus)
  ALLOCATE(kntua(nvar_anx,nsrc_ua),STAT=istatus)
  ALLOCATE(kntret(nvar_anx,nsrc_ret),STAT=istatus)
  ALLOCATE(kntrad(nvar_anx,nsrc_rad),STAT=istatus)

  ALLOCATE(bias(nvar_anx),STAT=istatus)
  ALLOCATE(biassngt(nvar_anx),STAT=istatus)
  ALLOCATE(biasuat(nvar_anx),STAT=istatus)
  ALLOCATE(biasrett(nvar_anx),STAT=istatus)
  ALLOCATE(biasradt(nvar_anx),STAT=istatus)
  ALLOCATE(biassng(nvar_anx,nsrc_sng),STAT=istatus)
  ALLOCATE(biasua(nvar_anx,nsrc_ua),STAT=istatus)
  ALLOCATE(biasret(nvar_anx,nsrc_ret),STAT=istatus)
  ALLOCATE(biasrad(nvar_anx,nsrc_rad),STAT=istatus)

  ALLOCATE(rms(nvar_anx),STAT=istatus)
  ALLOCATE(rmssngt(nvar_anx),STAT=istatus)
  ALLOCATE(rmsuat(nvar_anx),STAT=istatus)
  ALLOCATE(rmsrett(nvar_anx),STAT=istatus)
  ALLOCATE(rmsradt(nvar_anx),STAT=istatus)
  ALLOCATE(rmssng(nvar_anx,nsrc_sng),STAT=istatus)
  ALLOCATE(rmsua(nvar_anx,nsrc_ua),STAT=istatus)
  ALLOCATE(rmsret(nvar_anx,nsrc_ret),STAT=istatus)
  ALLOCATE(rmsrad(nvar_anx,nsrc_rad),STAT=istatus)

  ALLOCATE(rngsqi(nvar_anx),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:rngsqi")
  ALLOCATE(wgtsum(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:wgtsum")
  ALLOCATE(zsum(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:zsum")
  ALLOCATE(sqsum(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:sqsum")

  ALLOCATE(xs(nx),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:xs")
  ALLOCATE(ys(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ys")

  ALLOCATE(xslg(nxlg),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:xs")
  ALLOCATE(yslg(nylg),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:ys")

  ALLOCATE(zs(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:zs")
  ALLOCATE(latgr(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:latgr")
  ALLOCATE(longr(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:longr")

  ALLOCATE(tem1soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem1soil")
  ALLOCATE(tem2soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem2soil")
  ALLOCATE(tem3soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem3soil")
  ALLOCATE(tem4soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem4soil")
  ALLOCATE(tem5soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem5soil")

  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem1")
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem2")
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem3")
  ALLOCATE(tem4(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem4")
  ALLOCATE(tem5(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem5")
  ALLOCATE(tem6(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem6")
  ALLOCATE(tem7(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem7")
  ALLOCATE(tem1d(nz),STAT=istatus)
  CALL check_alloc_status(istatus, "difobs:tem1d")

  ALLOCATE(anx(nx,ny,nz,nvar_anx),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:anx")

  ALLOCATE(xcor(ncat,ncat),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:xcor")
  ALLOCATE(icatg(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:icatg")

!-----------------------------------------------------------------------
!
! Initialize the allocated arrays to zero.
!
!-----------------------------------------------------------------------

  x=0.0
  y=0.0
  z=0.0

  hterain=0.0
  mapfct =0.0

  zp  =0.0
  j1  =0.0
  j2  =0.0
  j3  =0.0
  aj3x=0.0
  aj3y=0.0
  aj3z=0.0
  j3inv=0.0

  u    =0.0
  v    =0.0
  w    =0.0
  wcont=0.0
  ptprt=0.0
  pprt =0.0

  qv   =0.0
  qscalar   =0.0
  tke  =0.0

  ubar =0.0
  vbar =0.0
  ptbar=0.0
  pbar =0.0
  rhostr=0.0
  qvbar=0.0

  udteb=0.0
  udtwb=0.0
  vdtnb=0.0
  vdtsb=0.0

  pdteb=0.0
  pdtwb=0.0
  pdtnb=0.0
  pdtsb=0.0

  trigs1=0.0
  trigs2=0.0
  ifax1=0.0
  ifax2=0.0

  vwork1=0.0
  vwork2=0.0

  wsave1=0.0
  wsave2=0.0

  ppi   =0.0
  csndsq=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0

  soiltyp=0.0
  stypfrct=0.0
  vegtyp =0.0
  lai    =0.0
  roufns =0.0
  veg    =0.0

  qvsfc  =0.0
  tsoil  =0.0
  qsoil  =0.0
  wetcanp=0.0
  snowdpth=0.0

  ptcumsrc=0.0
  qcumsrc =0.0
  w0avg   =0.0
  raing  =0.0
  rainc  =0.0
  prcrate=0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  knt=0.0
  kntsngt=0.0
  kntuat=0.0
  kntrett=0.0
  kntradt=0.0
  kntsng=0.0
  kntua=0.0
  kntret=0.0
  kntrad=0.0

  bias=0.0
  biassngt=0.0
  biasuat=0.0
  biasrett=0.0
  biasradt=0.0
  biassng=0.0
  biasua=0.0
  biasret=0.0
  biasrad=0.0

  rms=0.0
  rmssngt=0.0
  rmsuat=0.0
  rmsrett=0.0
  rmsradt=0.0
  rmssng=0.0
  rmsua=0.0
  rmsret=0.0
  rmsrad=0.0
  rngsqi=0.0
  sqsum=0.0

  xs=0.0
  ys=0.0
  xslg=0.0
  yslg=0.0
  zs=0.0
  latgr=0.0
  longr=0.0

  tem1=0.0
  tem2=0.0
  tem3=0.0
  tem4=0.0
  tem5=0.0
  tem6=0.0
  tem7=0.0
!
!-----------------------------------------------------------------------
!
!  Set expected _squared_ background error for each variable
!  This will depend on the particular background field used,
!  season, age of forecast, etc.
!
!-----------------------------------------------------------------------
!
  IF( myproc == 0 ) THEN
    PRINT *, 'Reading ', TRIM(backerrfil)
    OPEN(4,FILE=trim(backerrfil),STATUS='old')
    READ(4,'(a80)') srcback
    READ(4,'(a80)') header
    WRITE(6,'(//,a,/a)') 'Background Std Error for',srcback
    WRITE(6,'(1x,a)')                                                  &
      '   k    hgt(m)  u(m/s)  v(m/s) pres(mb) temp(K)  RH(%)'
    DO ktab=1,nz_tab
      READ(4,*,END=2) hqback(ktab),                                    &
                  qback(3,ktab),                                       &
                  qback(4,ktab),                                       &
                  qback(5,ktab),                                       &
                  qback(1,ktab),                                       &
                  qback(2,ktab)
      WRITE(6,'(1x,i4,f10.2,5f8.2)') ktab,hqback(ktab),                &
                                 (qback(ivar,ktab),ivar=1,5)
      qback(3,ktab)=100.*qback(3,ktab)
      qback(5,ktab)=0.01*qback(5,ktab)
    END DO
  2 nlvqback=ktab-1
    CLOSE(4)
!
    isrc=1
    DO jsrc=1,nsrc_sng
      IF(srcsng(jsrc) /= 'NULL') THEN
        PRINT *, 'Reading ', TRIM(sngerrfil(jsrc))
        OPEN(4,ERR=3,FILE=trim(sngerrfil(jsrc)),STATUS='old')
        READ(4,'(a8)',ERR=3,END=3) rdsource
        IF(rdsource /= srcsng(jsrc)) THEN
          WRITE(6,'(a,i4,a,a,a,a,a)')                                  &
              ' Mismatch of source names',jsrc,                        &
              ' read ',rdsource,' expected ',srcsng(jsrc)
          STOP
        END IF
        READ(4,'(a80)',ERR=3,END=3) srcsng_full(isrc)
        READ(4,*,ERR=3,END=3) qcmulsng(isrc)
        READ(4,'(a80)',ERR=3,END=3) header
        READ(4,*,ERR=3,END=3)                                          &
               qsrcsng(3,isrc),                                        &
               qsrcsng(4,isrc),                                        &
               qsrcsng(5,isrc),                                        &
               qsrcsng(1,isrc),                                        &
               qsrcsng(2,isrc)
        WRITE(6,'(//,a,a,/a)') 'Single-level std error for ',          &
                              srcsng(isrc),srcsng_full(isrc)
        WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulsng(isrc)
        WRITE(6,'(1x,a)')                                              &
            '   u (m/s)  v (m/s)  pres(mb) temp(K) RH(%)'
        WRITE(6,'(1x,5f8.2)') (qsrcsng(ivar,isrc),ivar=1,5)
        qsrcsng(3,isrc)=100.*qsrcsng(3,isrc)
        qsrcsng(5,isrc)=0.01*qsrcsng(5,isrc)
        CLOSE(4)
        isrc=isrc+1
      END IF
    END DO
  
  3 CONTINUE
!
    isrc=1
    DO jsrc=1,nsrc_ua
      IF(srcua(jsrc) /= 'NULL') THEN
        PRINT *, 'Reading ', TRIM(uaerrfil(jsrc))
        OPEN(4,ERR=6,FILE=trim(uaerrfil(jsrc)),STATUS='old')
        READ(4,'(a8)',ERR=6,END=6) rdsource
        IF(rdsource /= srcua(jsrc)) THEN
          WRITE(6,'(a,i4,a,a,a,a,a)')                                  &
              ' Mismatch of source names',jsrc,                        &
              ' read ',rdsource,' expected ',srcua(jsrc)
          STOP
        END IF
        READ(4,'(a80)',ERR=6,END=6) srcua_full(isrc)
        READ(4,*,ERR=6,END=6) qcmulua(isrc)
        READ(4,'(a80)',ERR=6,END=6) header
        WRITE(6,'(//,a,a,/a)') 'UA data std error for ',               &
                              srcua(isrc),srcua_full(isrc)
        WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulua(isrc)
        WRITE(6,'(1x,a)')                                              &
            '   k    hgt(m)  u(m/s)  v(m/s) pres(mb) temp(K)  RH(%)'
        DO ktab=1,nz_tab
          READ(4,*,ERR=5,END=5) huaqsrc(ktab,isrc),                    &
                  qsrcua(3,ktab,isrc),                                 &
                  qsrcua(4,ktab,isrc),                                 &
                  qsrcua(5,ktab,isrc),                                 &
                  qsrcua(1,ktab,isrc),                                 &
                  qsrcua(2,ktab,isrc)
          WRITE(6,'(1x,i4,f10.2,5f8.2)') ktab,huaqsrc(ktab,isrc),      &
                 (qsrcua(ivar,ktab,isrc),ivar=1,5)
          qsrcua(3,ktab,isrc)=100.*qsrcua(3,ktab,isrc)
          qsrcua(5,ktab,isrc)=0.01*qsrcua(5,ktab,isrc)
        END DO
  5     nlvuatab(isrc)=ktab-1
        CLOSE(4)
        isrc=isrc+1
      END IF
    END DO
  
  6 CONTINUE
  !
    isrc=1
    DO jsrc=1,nsrc_rad
      IF(srcrad(jsrc) /= 'NULL') THEN
        PRINT *, 'Reading ', TRIM(raderrfil(jsrc))
        OPEN(4,ERR=7,FILE=trim(raderrfil(jsrc)),STATUS='old')
        READ(4,'(a8)',ERR=7,END=7) rdsource
        IF(rdsource /= srcrad(jsrc)) THEN
          WRITE(6,'(a,i4,a,a,a,a,a)')                                     &
              ' Mismatch of source names',jsrc,                           &
              ' read ',rdsource,' expected ',srcrad(jsrc)
          STOP
        END IF
        READ(4,'(a80)',ERR=7,END=7) srcrad_full(isrc)
        READ(4,*,ERR=7,END=7) qcmulrad(isrc)
        READ(4,'(a80)',ERR=7,END=7) header
        READ(4,*,ERR=7,END=7)                                             &
               qsrcrad(1,isrc),                                           &
               qsrcrad(2,isrc),                                           &
               qsrcrad(3,isrc)
        WRITE(6,'(//,a,a,/a)') 'Radar data std error for ',               &
                            srcrad(isrc),srcrad_full(isrc)
        WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulrad(isrc)
        WRITE(6,'(1x,a)')                                                 &
            '    ref(dBz) Vrad(m/s) SpWid(m/s)'
        WRITE(6,'(1x,4f8.2)')                                             &
            (qsrcrad(ivar,isrc),ivar=1,nvar_radin)
        CLOSE(4)
        isrc=isrc+1
      END IF
    END DO
  
  7 CONTINUE
!
    isrc=1
    DO jsrc=1,nsrc_ret
      IF(srcret(jsrc) /= 'NULL') THEN
        PRINT *, 'Reading ', TRIM(reterrfil(jsrc))
        OPEN(4,FILE=trim(reterrfil(jsrc)),ERR=10,STATUS='old')
        READ(4,'(a8)',ERR=10,END=10) rdsource
        IF(rdsource /= srcret(jsrc)) THEN
          WRITE(6,'(a,i4,a,a,a,a,a)')                                     &
              ' Mismatch of source names',jsrc,                           &
              ' read ',rdsource,' expected ',srcret(jsrc)
          STOP
        END IF
        READ(4,'(a80)',ERR=10,END=10) srcret_full(isrc)
        READ(4,*,ERR=10,END=10) qcmulret(isrc)
        READ(4,'(a80)',ERR=10,END=10) header
        WRITE(6,'(//,a,a,/a)') 'Retrieval std error for ',                &
                            srcret(isrc),srcret_full(isrc)
        WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulrad(isrc)
        WRITE(6,'(1x,a)')                                                 &
            '   k    hgt(m)  u(m/s)  v(m/s) pres(mb) temp(K)  RH(%)'
        DO ktab=1,nz_tab
          READ(4,*,END=9,ERR=9) hrtqsrc(ktab,isrc),                       &
                  qsrcret(3,ktab,isrc),                                   &
                  qsrcret(4,ktab,isrc),                                   &
                  qsrcret(5,ktab,isrc),                                   &
                  qsrcret(1,ktab,isrc),                                   &
                  qsrcret(2,ktab,isrc)
          WRITE(6,'(1x,i4,f10.2,5f8.2)') ktab,hrtqsrc(ktab,isrc),         &
                        (qsrcret(ivar,ktab,isrc),ivar=1,5)
          qsrcret(3,ktab,isrc)=100.*qsrcret(3,ktab,isrc)
          qsrcret(5,ktab,isrc)=0.01*qsrcret(5,ktab,isrc)
        END DO
 9      nlvrttab(isrc)=ktab-1
        CLOSE(4)
        isrc=isrc+1
      END IF
    END DO

10  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Change standard error to standard error variance by squaring
!  Calculate quality control thresholds
!
!-----------------------------------------------------------------------
!
    DO ktab=1,nlvqback
      DO ivar=1,nvar_anx
        qback(ivar,ktab)=qback(ivar,ktab)*qback(ivar,ktab)
      END DO
    END DO
!
    DO isrc=1,nsrc_sng
      DO ivar=1,nvar_anx
        qsrcsng(ivar,isrc)=qsrcsng(ivar,isrc)*qsrcsng(ivar,isrc)
        qcthrsng(ivar,isrc)=qcmulsng(isrc)*                               &
                         SQRT(qback(ivar,1)+qsrcsng(ivar,isrc))
        barqclim(ivar,isrc)=qcmulsng(isrc)*SQRT(qsrcsng(ivar,isrc))
      END DO
      qcthrsng(iqv,isrc)=min(qcthrsng(iqv,isrc),0.75)
      barqclim(iqv,isrc)=min(barqclim(iqv,isrc),0.75)
    END DO
    musesng=0
    DO isrc=1,nsrc_sng
      DO ipass=1,npass
        musesng(isrc)=max(musesng(isrc),iusesng(isrc,ipass))
      END DO
    END DO
!
    DO isrc=1,nsrc_ua
      DO ivar=1,nvar_anx
        qsrcmax=0.
        DO ktab=1,nlvuatab(isrc)
          qsrcua(ivar,ktab,isrc) =                                        &
              qsrcua(ivar,ktab,isrc)*qsrcua(ivar,ktab,isrc)
          qsrcmax=AMAX1(qsrcmax,qsrcua(ivar,ktab,isrc))
        END DO
        qcthrua(ivar,isrc)=qcmulua(isrc)*                                 &
                           SQRT(qback(ivar,1)+qsrcmax)
      END DO
    END DO
    museua=0
    DO isrc=1,nsrc_ua
      DO ipass=1,npass
        museua(isrc)=max(museua(isrc),iuseua(isrc,ipass))
      END DO
    END DO
!
    DO isrc=1,nsrc_rad
      DO ivar=1,nvar_radin
        qsrcrad(ivar,isrc) =                                              &
            qsrcrad(ivar,isrc)*qsrcrad(ivar,isrc)
        qcthrrad(ivar,isrc)=qcmulrad(isrc)*                               &
                            SQRT(qback(ivar,1)+qsrcrad(ivar,isrc))
      END DO
    END DO
    muserad=0
    DO ipass=1,npass
      DO isrc=1,nsrc_rad
        muserad(isrc)=max(muserad(isrc),iuserad(isrc,ipass))
      END DO
    END DO
!
    DO isrc=1,nsrc_ret
      DO ivar=1,nvar_anx
        qsrcmax=0.
        DO ktab=1,nlvrttab(isrc)
          qsrcret(ivar,ktab,isrc) =                                       &
              qsrcret(ivar,ktab,isrc)*qsrcret(ivar,ktab,isrc)
          qsrcmax=AMAX1(qsrcmax,qsrcret(ivar,ktab,isrc))
        END DO
        qcthrret(ivar,isrc)=qcmulret(isrc)*                               &
                            SQRT(qback(ivar,1)+qsrcmax)
      END DO
    END DO
    museret=0
    DO ipass=1,npass
      DO isrc=1,nsrc_ret
        museret(isrc)=max(museret(isrc),iuseret(isrc,ipass))
      END DO
    END DO
  END IF   ! processor zero

!
!  Broadcast read variables to all procs.
!  Also broadcast other variables that may have been updated.
!

  CALL mpupdater(hqback, nz_tab)
  CALL mpupdater(qback, (nvar_anx*nz_tab))
  CALL mpupdatei(nlvqback, 1)

  CALL mpupdatec(srcsng,   (8*nsrc_sng))
  CALL mpupdater(qcthrsng, (nvar_anx*nsrc_sng))

  CALL mpupdater(qcthrua, (nvar_anx*nsrc_ua))

  CALL mpupdater(qcthrrad, (nvar_radin*nsrc_rad))

  CALL mpupdater(qcthrret, (nvar_anx*nsrc_ret))

  CALL mpupdater(qsrcrad,nvar_radin*nsrc_rad)

  CALL mpupdatei(iusesng,(nsrc_sng+1)*mx_pass)
  CALL mpupdatei(iuseua,(nsrc_ua+1)*mx_pass)
  CALL mpupdatei(iuserad,(nsrc_rad+1)*mx_pass)
  CALL mpupdatei(iuseret,(nsrc_ret+1)*mx_pass)

  CALL make_mpi_map(mpi_map,nmap,iproc,jproc,nx,ny)
  CALL make_mpi_map(mpi_mapv,nmapv,iprocv,jprocv,nx,ny)

!
!-----------------------------------------------------------------------
!
!  Set cross-correlations between numerical categories.
!  Roughly 1=clear,2=some evidence of outflow,3=precip,4=conv precip
!  This could be read in as a table.
!
!-----------------------------------------------------------------------
!
  DO j=1,ncat
    DO i=1,ncat
      xcor(i,j)=1.0-(IABS(i-j))*0.25
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Initialize grids, forming first guess fields based on
!  the model initial option specified in the ARPS init file.
!
!-----------------------------------------------------------------------
!
  CALL readarpsmp(ncompressx,ncompressy,nproc_node,                 &
              nprocx_lw,nprocy_lw,hinfmt,                               &
              grdbasfn,lengbf,hisfile(nf),lenfil,nx,ny,nz,nzsoil,nstyps,&
              time,x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,qvprt,     &
              qscalar,tke,kmh,kmv,                                      &
              ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                   &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,                   &
              tsoil,qsoil,wetcanp,snowdpth,                             &
              raing,rainc,prcrate,radfrc,radsw,rnflx,radswnet,radlwin,  &
              usflx,vsflx,ptsflx,qvsflx,istatus,tem1,tem2,tem3)

  CALL xytoll(nx,ny,x,y,latgr,longr)
!
!-----------------------------------------------------------------------
!
!  Set location of scalar fields.
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=2.0*xs(nx-1)-xs(nx-2)
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=2.0*ys(ny-1)-ys(ny-2)
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zs(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
  DO i=1,nx-1
    DO j=1,ny-1
      zs(i,j,nz)=2.*zs(i,j,nz-1)-zs(i,j,nz-2)
    END DO
  END DO

!-----------------------------------------------------------------------
!
! In the MPI world, some subroutines will need to know characteristics
! of the large grid.  To simplify later subroutine calls (one call not two
! surrounded by IF/ELSE/ENDIF, everybody that always needs the large
! grid will use the large grid variables.  Large grid info is used by
! early decision making when determining what obs are needed.
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0 ) THEN
    CALL mpimerge1dx(xs,nx,xslg)
    CALL mpimerge1dy(ys,ny,yslg)

    CALL mpupdater(xslg,nxlg)
    CALL mpupdater(yslg,nylg)
!
  ELSE
     xslg = xs
     yslg = ys

  END IF

  CALL ctim2abss( year,month,day,hour,minute,second, i4timanx)
!
!-----------------------------------------------------------------------
!
!  Identify the background field correlation category for
!  each 2-d point.   This is based on precip rate, cumulus
!  parameterization switch, and surface relative humidity.
!
!-----------------------------------------------------------------------
!
  CALL setcat(nx,ny,nz,nscalar,ccatopt,zs,                              &
              ptprt,pprt,qv,qscalar,                                    &
              ptbar,pbar,                                               &
              prcrate,icatg)
!
!-----------------------------------------------------------------------
!
!  Load "background fields" for analysis
!  Note, analysis is done at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        anx(i,j,k,1)=0.5*(u(i,j,k)+u(i+1,j,k))
        anx(i,j,k,2)=0.5*(v(i,j,k)+v(i,j+1,k))
        anx(i,j,k,3)=pbar(i,j,k)+pprt(i,j,k)
        anx(i,j,k,4)=ptbar(i,j,k)+ptprt(i,j,k)
        anx(i,j,k,5)=qv(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Filename creation
!
!-----------------------------------------------------------------------
!
  WRITE(stats,'(a,a)') runname(1:lfnkey),'.lst'
  CALL getunit(iwstat)
  PRINT *, 'Writing ', TRIM(stats)
  OPEN(iwstat,IOSTAT=ios,FILE=trim(stats),STATUS='unknown')
  IF(ios /= 0) iwstat=0

  maxsuf=LEN(suffix)
  mxroot=LEN(froot)
  nobsng=0

  IF (myproc == 0) THEN
    DO ifile=1,nsngfil
      PRINT *, 'Processing file ', ifile, ' ', TRIM(sngfname(ifile))
      lenfnm=LEN(sngfname(ifile))
      CALL strlnth(sngfname(ifile),lenfnm)
      CALL exsufx(sngfname(ifile),lenfnm,suffix,maxsuf,dotloc,lensuf)
      IF(lensuf == 3. .AND. suffix(1:3) == 'lso') THEN

        CALL exfroot(sngfname(ifile),lenfnm,froot,mxroot,lenroot)
        WRITE(qclist,'(a,a)') froot(1:lenroot),'.sqc'
        WRITE(qcsave,'(a,a)') froot(1:lenroot),'.lsq'
!
!-----------------------------------------------------------------------
!
!  Open the files for listing QC info
!  To suppress listing set the unit numbers to zero
!
!-----------------------------------------------------------------------
!
        CALL getunit(iqclist)
        PRINT *, 'Opening qclist: ', TRIM(qclist)
        OPEN(iqclist,IOSTAT=ios,FILE=trim(qclist),STATUS='unknown')
        IF(ios /= 0) iqclist=0
        CALL getunit(iqcsave)
        PRINT *, 'Opening qclist: ', TRIM(qcsave)
        OPEN(iqcsave,IOSTAT=ios,FILE=trim(qcsave),STATUS='unknown')
        IF(ios /= 0) iqcsave=0
!
!-----------------------------------------------------------------------
!
!  Read surface data, QC and convert units.
!
!-----------------------------------------------------------------------
!
        rngsq=sfcqcrng*sfcqcrng
        PRINT *, 'Calling prepsfcobs'
        CALL prepsfcobs(ntime,mx_sng,                                  &
          nvar_sng,nvar_anx,nsrc_sng,nobsng,ipres,iptmp,iqv,           &
          sngfname(ifile),sngtmchk(ifile),blackfil,                    &
          var_sfc,var_anx,srcsng,qsrcsng,                              &
          rmiss,iqspr,sprdist,nobsrd,timesng,                          &
          stnsng,csrcsng,isrcsng,latsng,lonsng,hgtsng,xsng,ysng,       &
          nxlg,nylg,xslg,yslg,                                         &
          wx,kloud,idp3,store_emv,store_hgt,store_amt,                 &
          obrdsng,obsng,qualrdsng,qobsng,qualsng,                      &
          ival,climin,climax,                                          &
          rngsq,klim,wlim,qcthrsng,barqclim,                           &
          knt,bias,rms,sqsum,iqclist,iqcsave,jstatus)

        IF(iqclist /= 0) THEN
          CLOSE(iqclist)
          CALL retunit(iqclist)
        END IF
        IF(iqcsave /= 0) THEN
          CLOSE(iqcsave)
          CALL retunit(iqcsave)
        END IF

      ELSE
!
!-----------------------------------------------------------------------
!
!  Read other single-level data and convert units.
!
!-----------------------------------------------------------------------
!
        PRINT *, 'Calling prepsglobs'
        CALL prepsglobs(mx_sng,ntime,srcsng,nsrc_sng,                  &
              sngfname(ifile),stnsng,latsng,lonsng,xsng,ysng,          &
              hgtsng,obsng,qobsng,qualsng,isrcsng,qsrcsng,             &
              csrcsng,nxlg,nylg,nz,nvar_anx,anx,xslg,yslg,zp,          &
              tem1,tem2,tem3,tem4,tem5,tem6,                           &
              rmiss,nobsng,istatus)

      END IF
    END DO
  END IF ! myproc == 0
!
!-----------------------------------------------------------------------
!
!  Broadcast the number of data points, then the data.  If
!  there are no data values, we don't waste anyone's time.
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0) THEN
    CALL mpupdatei(nobsng,1)
    IF (nobsng > 0) THEN
      CALL mpupdatei(timesng,mx_sng*ntime)
      CALL mpupdatec(stnsng,mx_sng*ntime*5)
      CALL mpupdatec(csrcsng,mx_sng*ntime*8)
      CALL mpupdatei(isrcsng,mx_sng)
      CALL mpupdater(latsng,mx_sng*ntime)
      CALL mpupdater(lonsng,mx_sng*ntime)
      CALL mpupdater(hgtsng,mx_sng*ntime)
      CALL mpupdater(xsng,mx_sng)
      CALL mpupdater(ysng,mx_sng)
      CALL mpupdatec(wx,mx_sng*ntime*8)
      CALL mpupdatei(kloud,mx_sng*ntime)
      CALL mpupdatei(idp3,mx_sng*ntime)
      CALL mpupdatec(store_emv,mx_sng*5*ntime)
      CALL mpupdater(store_hgt,mx_sng*5*ntime)
      CALL mpupdatec(store_amt,mx_sng*5*ntime*4)
      CALL mpupdater(obrdsng,mx_sng*nvar_sng*ntime)
      CALL mpupdater(obsng,nvar_anx*mx_sng)
      CALL mpupdater(qobsng,nvar_anx*mx_sng)
      CALL mpupdatei(qualsng,nvar_anx*mx_sng)
      CALL mpupdater(knt,nvar_anx*nz)
      CALL mpupdater(wgtsum,nvar_anx*nz)
      CALL mpupdater(zsum,nvar_anx*nz)

!
!-----------------------------------------------------------------------
!      We need to know which processor "owns" which obs, so they can be
!      consulted on an "as needed" basis.
!-----------------------------------------------------------------------
!
!
      ALLOCATE(item1(nobsng),STAT=istatus)
      CALL check_alloc_status(istatus, "adas:item1:nobsng")
      ALLOCATE(ksng(nprocs),STAT=istatus)
      CALL check_alloc_status(istatus, "adas:ksng:nprocs")

      CALL mpiprocess(nobsng,indexsng,nprocs,ksng,ksngmax,            &
          isrcsng,item1,nx,ny,xsng,ysng,xs,ys)
      DEALLOCATE (item1)
!
!-----------------------------------------------------------------------
!     Mark obs that we don't "own" so we don't do unnecessary
!     computations.
!-----------------------------------------------------------------------
!
!
      ALLOCATE(item1(0:nprocs-1),STAT=istatus)
      CALL check_alloc_status(istatus, "adas:item1:nprocs")
!
      item1 = 0
!
      DO k=1,nmap
        IF (mpi_map(k,1) .NE. -1 ) item1(mpi_map(k,1)) = 1
      END DO
      item1(myproc) = 1

      DO k=1,nobsng
        IF (indexsng(k) == -1) CYCLE
        IF (isrcsng(k) == -1) CYCLE
        IF (item1(indexsng(k)) == 0) isrcsng(k) = 0
      END DO

      DEALLOCATE (item1)
!
!-----------------------------------------------------------------------
!
!  Although we're similar to PREPSFCGLOBS, and need to do one final task
!  (compute heights), we now use the local grid variables, not the
!  global.
!
!-----------------------------------------------------------------------
!
      CALL prepsglmdcrs_sm(mx_sng,ntime,latsng,lonsng,                 &
          xsng,ysng,hgtsng,csrcsng,                                    &
          nx,ny,nz,nvar_anx,anx,xs,ys,zp,                              &
          tem1,tem2,tem3,tem4,tem5,tem6,                               &
          rmiss,nobsng,indexsng,istatus)

    END IF
!
! Collect and distribute the data.
!

    IF (mp_opt > 0 .AND. nobsng > 0) THEN
      ALLOCATE(tems1dr(ksngmax),STAT=istat)
      CALL check_alloc_status(istat, "mpi_1dr_collect:tmps")
      ALLOCATE(temr1dr(ksngmax),STAT=istat)
      CALL check_alloc_status(istat, "mpi_1dr_collect:tmpr")

      CALL mpi_1dr_collect(hgtsng,nobsng,indexsng,                     &
        nprocs, ksng, ksngmax, mpi_map, nmap, tems1dr, temr1dr)
    END IF

  END IF  ! processor zero
!
!-----------------------------------------------------------------------
!
!  Calculate initial obs differences for single-level data
!
!-----------------------------------------------------------------------
!
! IF (myproc==0)
    PRINT *, 'Calling grdtosng'
    CALL grdtosng(nx,ny,nxlg,nylg,nz,nz_tab,mx_sng,nvar_anx,nobsng,    &
                    xs,ys,xslg,yslg,zp,icatg,anx,qback,hqback,nlvqback,&
                    tem1,tem2,tem3,tem4,tem5,tem6,                     &
                    stnsng,isrcsng,icatsng,hgtsng,xsng,ysng,           &
                    obsng,qobsng,qualsng,                              &
                    odifsng,oanxsng,thesng,trnsng,indexsng)
!
!-----------------------------------------------------------------------
!
!      Some data needs to be collected and distributed.
!
!-----------------------------------------------------------------------
!
!  IF ( mp_opt > 0 .AND. nobsng > 0 ) THEN
!    ALLOCATE(tems2dr(nvar_anx,ksngmax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_2dr_collect:tmps")
!    ALLOCATE(temr2dr(nvar_anx,ksngmax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_2dr_collect:tmpr")
!
!!
!!  1D real temporary arrays already allocated, but 1D ints aren't.
!!
!
!    ALLOCATE(tems1di(ksngmax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_1di_collect:tmps")
!    ALLOCATE(temr1di(ksngmax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_1di_collect:tmpr")
!
!    CALL mpi_2dr_collect( qobsng, nvar_anx, mx_sng, nobsng, indexsng,  &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems2dr, temr2dr)
!    CALL mpi_2dr_collect( odifsng, nvar_anx, mx_sng, nobsng, indexsng,  &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems2dr, temr2dr)
!    CALL mpi_2dr_collect( oanxsng, nvar_anx, mx_sng, nobsng, indexsng,  &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems2dr, temr2dr)
!
!    CALL mpi_1di_collect( icatsng, nobsng, indexsng,                    &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems1di, temr1di)
!
!    CALL mpi_1dr_collect( thesng, nobsng, indexsng,                     &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems1dr, temr1dr)
!    CALL mpi_1dr_collect( trnsng, nobsng, indexsng,                     &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems1dr, temr1dr)
!
!!
!!   2D space can't be deallocated, as the arrays get passed to "anxiter".
!!
!
!    DEALLOCATE(tems1di)
!    DEALLOCATE(temr1di)
!    DEALLOCATE(tems1dr)
!    DEALLOCATE(temr1dr)
!  END IF

!-----------------------------------------------------------------------
!
!  Read upper-air data, QC and convert units
!
!-----------------------------------------------------------------------
!
!  IF (myproc == 0) THEN
    PRINT *, 'Calling prepuaobs'
    CALL prepuaobs(nx,ny,nz,nvar_anx,                                  &
                   nz_ua,mx_ua,nz_tab,nsrc_ua,mx_ua_file,              &
                   anx,xs,ys,zp,tem1,tem2,tem3,tem4,tem5,tem6,tem7,    &
                   nuafil,uafname,srcua,                               &
                   stnua,elevua,xua,yua,hgtua,obsua,                   &
                   qsrcua,huaqsrc,nlvuatab,                            &
                   qobsua,qualua,isrcua,nlevsua,                       &
                   rmiss,nobsua,istatus)

!  END IF
!
!-----------------------------------------------------------------------
!
!  MPI issues.  Broadcast the number of data points, then the data.  If
!  there are no data values, we don't waste anyone's time.
!
!-----------------------------------------------------------------------
!
!  IF (mp_opt > 0) THEN
!    CALL mpupdatei(nobsua,1)
!    IF (nobsua > 0) THEN
!      CALL mpupdatec(stnua,mx_ua*5)
!      CALL mpupdater(elevua,mx_ua)
!      CALL mpupdater(xua,mx_ua)
!      CALL mpupdater(yua,mx_ua)
!      CALL mpupdater(hgtua,nz_ua*mx_ua)
!      CALL mpupdater(obsua,nvar_anx*nz_ua*mx_ua)
!      CALL mpupdater(qobsua,nvar_anx*nz_ua*mx_ua)
!      CALL mpupdatei(nlvuatab,nsrc_ua)
!      CALL mpupdatei(qualua,nvar_anx*nz_ua*mx_ua)
!      CALL mpupdatei(isrcua,mx_ua)
!      CALL mpupdatei(nlevsua,mx_ua)
!
!!
!!-----------------------------------------------------------------------
!!      We need to know which processor "owns" which obs, so they can be
!!      consulted on an "as needed" basis.
!!-----------------------------------------------------------------------
!!
!
!       ALLOCATE(item1(nobsua),STAT=istatus)
!       CALL check_alloc_status(istatus, "adas:item1:nobsua")
!       ALLOCATE(kua(nprocs),STAT=istatus)
!       CALL check_alloc_status(istatus, "adas:kua:nprocs")
!
!       CALL mpiprocess(nobsua,indexua,nprocs,kua,kuamax,               &
!          isrcua,item1,nx,ny,xua,yua,xs,ys)
!
!       DEALLOCATE (item1)
!
!!
!!-----------------------------------------------------------------------
!!      Mark obs that we don't "own" so we don't do unnecessary
!!      computations.
!!-----------------------------------------------------------------------
!!
!
!       ALLOCATE(item1(0:nprocs-1),STAT=istatus)
!       CALL check_alloc_status(istatus, "adas:item1:nprocs")
!
!       item1 = 0
!
!       DO k=1,nmap
!         IF (mpi_map(k,1) .NE. -1 ) item1(mpi_map(k,1)) = 1
!       END DO
!       item1(myproc) = 1
!
!!
!!-----------------------------------------------------------------------
!!      Multi-level obs are never combined, so we don't have to worry.
!!-----------------------------------------------------------------------
!!
!       DO k=1,nobsua
!         IF (indexua(k) == -1) CYCLE
!         IF (item1(indexua(k)) == 0) isrcua(k) = 0
!       END DO
!
!       DEALLOCATE (item1)
!
!    END IF
!  END IF


!
!-----------------------------------------------------------------------
!
!  Calculate initial obs differences for upper-air data
!
!-----------------------------------------------------------------------
!
!  IF (myproc == 0)
   PRINT *, 'Calling grdtoua'
  CALL grdtoua(nx,ny,nxlg,nylg,nz,nz_tab,nz_ua,mx_ua,nvar_anx,nobsua,   &
              xs,ys,xslg,yslg,zp,anx,qback,hqback,nlvqback,             &
              tem1,tem2,tem3,tem4,tem5,tem6,                            &
              stnua,isrcua,elevua,xua,yua,hgtua,                        &
              obsua,qobsua,qualua,nlevsua,                              &
              odifua,oanxua,theua,trnua,indexua)
!
!-----------------------------------------------------------------------
!
!      Some data needs to be collected and distributed.
!
!-----------------------------------------------------------------------
!
!  IF ( mp_opt > 0 .AND. nobsua > 0 ) THEN
!    ALLOCATE(tems3dr(nvar_anx,nz_ua,kuamax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_2dr_collect:tmps")
!    ALLOCATE(temr3dr(nvar_anx,nz_ua,kuamax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_2dr_collect:tmpr")
!
!    ALLOCATE(tems1dr(kuamax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_1dr_collect:tmps")
!    ALLOCATE(temr1dr(kuamax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_1dr_collect:tmpr")
!
!    ALLOCATE(tems2drua(nz_ua,kuamax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_2dr_collect:tmpsua")
!    ALLOCATE(temr2drua(nz_ua,kuamax),STAT=istat)
!    CALL check_alloc_status(istat, "mpi_2dr_collect:tmprua")
!
!    CALL mpi_3dr_collect( qobsua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
!      nprocs, kua, kuamax, mpi_map, nmap, tems3dr, temr3dr)
!    CALL mpi_3dr_collect( qualua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
!      nprocs, kua, kuamax, mpi_map, nmap, tems3dr, temr3dr)
!    CALL mpi_3dr_collect( odifua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
!      nprocs, kua, kuamax, mpi_map, nmap, tems3dr, temr3dr)
!    CALL mpi_3dr_collect( oanxua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
!      nprocs, kua, kuamax, mpi_map, nmap, tems3dr, temr3dr)
!    CALL mpi_1dr_collect( trnua, nobsua, indexua,                           &
!       nprocs, kua, kuamax, mpi_map, nmap, tems1dr, temr1dr)
!    CALL mpi_2dr_collect( theua, nz_ua, mx_ua, nobsua, indexua,  &
!      nprocs, kua, kuamax, mpi_map, nmap, tems2drua, temr2drua)
!
!    DEALLOCATE(tems1dr)
!    DEALLOCATE(temr1dr)
!    DEALLOCATE(tems2drua)
!    DEALLOCATE(temr2drua)
!  END IF
!
!-----------------------------------------------------------------------
!
!  Read radar data, unfold and convert units
!
!-----------------------------------------------------------------------
!
  IF(cloudopt > 0) THEN
    PRINT *, ' allocating full ref_mos_3d',nx,ny,nz
    ALLOCATE(ref_mos_3d(nx,ny,nz),STAT=istatus)
  ELSE
    PRINT *, ' allocating placeholder ref_mos_3d'
    ALLOCATE(ref_mos_3d(1,1,1),STAT=istatus)
  END IF
  CALL check_alloc_status(istatus, "adas:ref_mos_3d")
  ref_mos_3d=-99.

  raduvobs=1
!  IF (myproc==0)
  PRINT *, 'Calling prepradar nradfil=',nradfil

  CALL prepradar(nx,ny,nz,nz_tab,nvar_anx,nvar_radin,                   &
             nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad,mx_pass,         &
             raduvobs,radrhobs,radistride,radkstride,cloudopt,          &
             iuserad,npass,refrh,rhradobs,                              &
             xs,ys,zs,hterain,latgr,longr,anx,qback,hqback,nlvqback,    &
             nradfil,radfname,srcrad,isrcrad,qsrcrad,qcthrrad,          &
             stnrad,latrad,lonrad,elvrad,                               &
             latradc,lonradc,xradc,yradc,irad,nlevrad,                  &
             distrad,uazmrad,vazmrad,hgtradc,theradc,trnradc,dsdr,dhdr, &
             obsrad,oanxrad,odifrad,qobsrad,qualrad,                    &
             ncolrad,ncolrad_mpi,indexrad,oindexrad,ref_mos_3d,         &
             iprocv,jprocv,iradvel,                                     &
             istatus,istat_radar,tem1,tem2,tem3,tem4,tem5,tem6)
  print *, ' After prepradar ncolrad=',ncolrad

  DEALLOCATE(ref_mos_3d)
!
!-----------------------------------------------------------------------
!
!  Read retrieval data.
!
!-----------------------------------------------------------------------
!
!  IF (myproc == 0)
  PRINT *, 'Calling prepretr'
  CALL prepretr(nx,ny,nz,nvar_anx,                                      &
                nz_ret,mx_ret,mx_colret,nz_tab,nsrc_ret,                &
                nretfil,retfname,                                       &
                isrcret,srcret,nlvrttab,qsrcret,hrtqsrc,                &
                stnret,latret,lonret,elvret,                            &
                latretc,lonretc,xretc,yretc,iret,nlevret,               &
                hgtretc,obsret,qrret,qobsret,qualret,                   &
                rmiss,ncolret,tem1,istatus)

  PRINT *, 'Calling grdtoret'
  CALL grdtoret(nx,ny,nz,nz_tab,                                        &
                nz_ret,mx_ret,mx_colret,nvar_anx,ncolret,               &
                xs,ys,zp,anx,qback,hqback,nlvqback,                     &
                tem1,tem2,tem3,tem4,tem5,tem6,                          &
                stnret,iret,xretc,yretc,hgtretc,                        &
                obsret,qobsret,qualret,nlevret,                         &
                odifret,oanxret,theretc)
!
!-----------------------------------------------------------------------
!
!  Quality-control observation differences
!
!-----------------------------------------------------------------------
!
!  IF( myproc == 0)
  PRINT *, 'Calling qcdiff'
  CALL qcdiff(nvar_anx,nvar_rad,nvar_radin,mx_sng,nsrc_sng,             &
              indexsng,indexua,indexrad,indexret,                       &
              usesng,useua,userad,                                      &
              nz_ua,mx_ua,nsrc_ua,                                      &
              nz_rdr,mx_rad,mx_colrad,nsrc_rad,                         &
              nz_ret,mx_ret,mx_colret,nsrc_ret,                         &
              nobsng,nobsua,ncolrad,ncolret,var_anx,                    &
              stnsng,isrcsng,hgtsng,obsng,oanxsng,odifsng,              &
              qcthrsng,qualsng,                                         &
              stnua,isrcua,hgtua,obsua,oanxua,odifua,                   &
              qcthrua,qualua,nlevsua,                                   &
              stnrad,irad,isrcrad,hgtradc,obsrad,odifrad,               &
              qcthrrad,qualrad,nlevrad,                                 &
              stnret,iret,isrcret,hgtretc,                              &
              obsret,oanxret,odifret,                                   &
              qcthrret,qualret,nlevret,                                 &
              bias,rms)

!  IF( myproc == 0)
  PRINT *, 'Back from qcdiff'

!  IF ( mp_opt > 0 .AND. nobsng > 0) THEN
!    CALL mpi_2dr_collect( qualsng, nvar_anx, mx_sng, nobsng, indexsng,  &
!      nprocs, ksng, ksngmax, mpi_map, nmap, tems2dr, temr2dr)
!  END IF
!  IF ( mp_opt > 0 .AND. nobsua > 0 ) THEN
!    CALL mpi_3dr_collect( qualua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
!      nprocs, kua, kuamax, mpi_map, nmap, tems3dr, temr3dr)
!  END IF
!
!------------------------------------------------------------------------
!
! Build the mosaiked radar grid
! Output is stored in tem3.
!
!------------------------------------------------------------------------
!
!  IF( myproc == 0)
  PRINT *, 'Calling refmosaic'
  CALL refmosaic(nradfil,nx,ny,nz,mx_rad,                               &
           xs,ys,zs,radfname,lvldbg,tem4,rhinf,rvinf,                   &
           tem1,tem2,tem3,istatus)
!  IF( myproc == 0)
  PRINT *, 'Back from refmosaic'
!
!------------------------------------------------------------------------
!
! Calculate a composite reflectivity.  Max reflectivity in column.
! Store it in the sub-terrainian grid level (k=1) of tem3
!
!------------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      refmax=tem4(i,j,2)
      DO k=3,nz-1
        refmax=max(refmax,tem4(i,j,k))
      END DO
      tem4(i,j,1)=refmax
    END DO
  END DO
!
!------------------------------------------------------------------------
!
! Calculate the reflectivity from the model hydrometeor variables.
! First get the temperature and density.  Store in tem1 and tem2,
! respectively.
!
! Ferrier reflectivity is returned in tem5.
!
!------------------------------------------------------------------------
!
  p0inv=1./p0
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=(ptbar(i,j,k)+ptprt(i,j,k))*                        &
                    (((pbar(i,j,k)+pprt(i,j,k))*p0inv)**rddcp)
        tem2(i,j,k)=(pbar(i,j,k)+pprt(i,j,k))/(rd*tem1(i,j,k))
      END DO
    END DO
  END DO
  CALL reflec_ferrier(nx,ny,nz, tem2, qscalar, tem1, tem5)
!
!------------------------------------------------------------------------
!
! Calculate a composite reflectivity.  Max reflectivity in column.
! Store it in the sub-terrainian grid level (k=1) of tem4
!
!------------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      refmax=tem5(i,j,2)
      DO k=3,nz-1
        refmax=max(refmax,tem5(i,j,k))
      END DO
      tem5(i,j,1)=refmax
    END DO
  END DO
!
!------------------------------------------------------------------------
!
!   Calculate the statistical data
!
!------------------------------------------------------------------------
!
  PRINT *, 'Calling difstats...'
  CALL difstats(nx,ny,nz,                                               &
           nvar_anx,nvar_radin,nvar_rad,nz_ua,nz_rdr,nz_ret,            &
           mx_sng,mx_ua,mx_rad,mx_colrad,mx_ret,mx_colret,              &
           nsrc_sng,nsrc_ua,nsrc_rad,nsrc_ret,                          &
           xs,ys,zp,w,                                                  &
           xsng,ysng,hgtsng,thesng,                                     &
           obsng,odifsng,qobsng,qualsng,isrcsng,musesng,nobsng,         &
           xua,yua,hgtua,theua,                                         &
           obsua,odifua,qobsua,qualua,isrcua,museua,nlevsua,nobsua,     &
           elvrad,xradc,yradc,                                          &
           distrad,uazmrad,vazmrad,hgtradc,theradc,wradc,               &
           obsrad,odifrad,qobsrad,qualrad,                              &
           irad,isrcrad,muserad,nlevrad,ncolrad,                        &
           xretc,yretc,hgtretc,theretc,                                 &
           obsret,odifret,qobsret,qualret,                              &
           iret,isrcret,museret,nlevret,ncolret,                        &
           srcsng,srcua,srcrad,srcret,                                  &
           tem5,tem4,                                                   &
           knt,bias,rms,                                                &
           kntsngt,biassngt,rmssngt,kntuat,biasuat,rmsuat,              &
           kntrett,biasrett,rmsrett,kntradt,biasradt,rmsradt,           &
           kntsng,biassng,rmssng,kntua,biasua,rmsua,kntret,             &
           biasret,rmsret,kntrad,biasrad,rmsrad,                        &
           oanxsng,oanxua,oanxrad,oanxret,                              &
           tem1d,istatus)

!-------------------------------------------------------------------
!
!   Write stats to a file
!
!-------------------------------------------------------------------

  CALL gtlfnkey(runname, lfnkey)

  veriffn  = runname(1:lfnkey)//'.difobs'
  lveriffn = 7 + lfnkey

  WRITE(6,'(1x,a,a,a/,1x,a)')                                          &
        'Check to see if file ',veriffn(1:lveriffn),' already exists.',&
        'If so, append a version number to the filename.'

  CALL fnversn( veriffn, lveriffn )

  CALL getunit ( nchdif)

  WRITE(6,'(1x,a)')                                                    &
        'Writing output statistics to ',veriffn(1:lveriffn)

  OPEN(nchdif,FORM='formatted',STATUS='new',                           &
                  FILE=veriffn(1:lveriffn),IOSTAT=istat)

  WRITE(nchdif,'(a,/a)') '  Global Statistics',                      &
             '   ivar        knt    bias        rms'
  score=0.
  DO k=1,nvar_anx
    score=score+wgtvar(k)*rms(k)
    WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                          &
    var_anx(k),knt(k),bias(k),rms(k)
  END DO
  WRITE(nchdif,'(a,f10.2)') '   rms score = ',score

  WRITE(nchdif,'(a)')                                                &
  '-----------------------------------------------------'

  totsng=0
  DO k=1,nvar_anx
    totsng=totsng+kntsngt(k)
  END DO
  IF(totsng > 0 ) THEN
    WRITE(nchdif,'(/a,/a)') '  Statistics for all single-level data',&
             '   var     knt    bias        rms'
    score=0.
    DO k=1,nvar_anx
      score=score+wgtvar(k)*rmssngt(k)
      WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                         &
            var_anx(k),kntsngt(k),biassngt(k),rmssngt(k)
    END DO
    WRITE(nchdif,'(a,f10.2)') '   rms score = ',score

    WRITE(nchdif,'(/a)') '  Statistics by Source - Single'
    DO isrc=1,nsrc_sng
      totsrc=0
      DO k=1,nvar_anx
        totsrc=totsrc+kntsng(k,isrc)
      END DO
      IF(totsrc > 0 ) THEN
        WRITE(nchdif,'(/a,a,/a)') '   Source: ',srcsng(isrc),        &
             '   var     knt    bias        rms'
        score=0.
        DO k=1,nvar_anx
          score=score+wgtvar(k)*rmssng(k,isrc)
          WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                     &
          var_anx(k),kntsng(k,isrc),biassng(k,isrc),rmssng(k,isrc)
        END DO
        WRITE(nchdif,'(a,f10.2)') '   rms score = ',score
      ELSE
        WRITE(nchdif,'(/a,a,a)')                                     &
              '   Source: ',srcsng(isrc),' zero data used'
      END IF
    END DO
  ELSE
    WRITE(nchdif,'(/a)') '   Zero Single-Level data'
  END IF

  WRITE(nchdif,'(/a)')                                               &
  ' -----------------------------------------------------'

  totua=0
  DO k=1,nvar_anx
    totua=totua+kntuat(k)
  END DO
  IF(totua > 0 ) THEN
    WRITE(nchdif,'(/a,/a)') '  Statistics for all upper air data',   &
             '   var         knt    bias        rms'
    score=0.
    DO k=1,nvar_anx
      score=score+wgtvar(k)*rmsuat(k)
      WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                        &
      var_anx(k),kntuat(k),biasuat(k),rmsuat(k)
    END DO
    WRITE(nchdif,'(a,f10.2)') '   rms score = ',score

    WRITE(nchdif,'(/a)') 'Statistics by Source - Upper Air'
    DO isrc=1,nsrc_ua
      totsrc=0
      DO k=1,nvar_anx
        totsrc=totsrc+kntua(k,isrc)
      END DO
      IF(totsrc > 0 ) THEN
        WRITE(nchdif,'(/a,a,/a)') '  Source: ',srcua(isrc),          &
             '   var         knt    bias        rms'
        score=0.
        DO k=1,nvar_anx
          score=score+wgtvar(k)*rmsua(k,isrc)
          WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                    &
          var_anx(k),kntua(k,isrc),biasua(k,isrc),rmsua(k,isrc)
        END DO
        WRITE(nchdif,'(a,f10.2)') '   rms score = ',score
      ELSE
        WRITE(nchdif,'(/a,a,a)')                                     &
              '   Source: ',srcret(isrc),' zero data used'
      END IF
    END DO
  ELSE
    WRITE(nchdif,'(/a)') '   Zero Upper Air data'
  END IF

  WRITE(nchdif,'(/a)')                                               &
  ' -----------------------------------------------------'

  totret=0
  DO k=1,nvar_anx
    totret=totret+kntrett(k)
  END DO
  IF(totret > 0 ) THEN
    WRITE(nchdif,'(/a,/a)') '  Statistics for all retrieval data',   &
             '   var         knt    bias        rms'
    score=0.
    DO k=1,nvar_anx
      score=score+wgtvar(k)*rmsrett(k)
      WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                        &
      var_anx(k),kntrett(k),biasrett(k),rmsrett(k)
    END DO
    WRITE(nchdif,'(a,f10.2)') '   rms score = ',score

    WRITE(nchdif,'(/a)') '  Statistics by Source - Retrieval'
    DO isrc=1,nsrc_ret
      totsrc=0
      DO k=1,nvar_anx
        totsrc=totsrc+kntret(k,isrc)
      END DO
      IF(totsrc > 0 ) THEN
        WRITE(nchdif,'(/a,a,/a)') '   Source: ',srcret(isrc),        &
             '   var         knt    bias        rms'
        score=0.
        DO k=1,nvar_anx
          score=score+wgtvar(k)*rmsret(k,isrc)
          WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                    &
          var_anx(k),kntret(k,isrc),biasret(k,isrc),rmsret(k,isrc)
        END DO
        WRITE(nchdif,'(a,f10.2)') '   rms score = ',score
      ELSE
        WRITE(nchdif,'(/a,a,a)')                                     &
              '   Source: ',srcret(isrc),' zero data used'
      END IF
    END DO
  ELSE
    WRITE(nchdif,'(/a)') '   Zero Retrieval data'
  END IF

  WRITE(nchdif,'(/a)')                                               &
  ' -----------------------------------------------------'

  totrad=0
  DO k=1,nvar_anx
    totrad=totrad+kntradt(k)
  END DO
  IF(totrad > 0 ) THEN
    WRITE(nchdif,'(/a,/a)')  '  Radar data variables ',              &
             '   var         knt    bias        rms'
    score=0.
    DO k=1,3
      WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                        &
            var_rad(k),kntradt(k),biasradt(k),rmsradt(k)
    END DO

    WRITE(nchdif,'(/a)') '  Statistics by source - Radar data'
    DO isrc=1,nsrc_rad
      totsrc=0
      DO k=1,nvar_anx
        totsrc=totsrc+kntrad(k,isrc)
      END DO
      IF(totsrc > 0 ) THEN
        WRITE(nchdif,'(/a,a,/a)') '  Source: ',srcrad(isrc),         &
             '   var        knt    bias        rms'
        WRITE(nchdif,'(2x,a6,i10,g12.3,g11.3)')                       &
        var_rad(1),kntrad(1,isrc),biasrad(1,isrc),rmsrad(1,isrc)
      ELSE
        WRITE(nchdif,'(/a,a,a)')                                     &
              '   Source: ',srcrad(isrc),' zero data used'
      END IF
    END DO
  ELSE
    WRITE(nchdif,'(/a)') '   Zero Radar data'
  END IF

  WRITE(nchdif,'(/a)')                                               &
  ' -----------------------------------------------------'

  CLOSE(nchdif)

END PROGRAM difobs
!
SUBROUTINE exsufx(fname,lenfnm,suffix,maxsuf,dotloc,lensuf)
  IMPLICIT NONE
  INTEGER :: lenfnm
  INTEGER :: maxsuf
  CHARACTER (LEN=1) :: fname(lenfnm)
  CHARACTER (LEN=1) :: suffix(maxsuf)
  INTEGER :: lensuf
  INTEGER :: dotloc
!
  INTEGER :: i
!
  DO i=lenfnm,1,-1
    IF(fname(i) == '.') GO TO 200
  END DO
  dotloc=0
  lensuf=0
  DO i=1,maxsuf
    suffix(i)=' '
  END DO
  RETURN
  200 CONTINUE
  dotloc=i
  lensuf=MIN(maxsuf,lenfnm-i)
  DO i=1,lensuf
    suffix(i)=fname(dotloc+i)
  END DO
  RETURN
END SUBROUTINE exsufx
!
SUBROUTINE exfroot(fname,lenfnm,froot,mxroot,lenroot)
  IMPLICIT NONE
  INTEGER :: lenfnm
  INTEGER :: mxroot
  CHARACTER (LEN=1) :: fname(lenfnm)
  CHARACTER (LEN=1) :: froot(mxroot)
  INTEGER :: lenroot
!
  INTEGER :: i
  INTEGER :: slashloc
  INTEGER :: dotloc
!
  DO i=lenfnm,1,-1
    IF(fname(i) == '/') EXIT
  END DO
  101 CONTINUE
  slashloc=i
  DO i=slashloc,lenfnm
    IF(fname(i) == '.') EXIT
  END DO
  201 CONTINUE
  dotloc=i
  lenroot=(dotloc-slashloc)-1
  DO i=1,lenroot
    froot(i)=fname(slashloc+i)
  END DO
  RETURN
END SUBROUTINE exfroot
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SETCAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE setcat(nx,ny,nz,nscalar,ccatopt,zs,                          &
           ptprt,pprt,qv,qscalar,                                       &
           ptbar,pbar,                                                  &
           prcrate,icatg)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign categories to grid locations according to the background
!  state to account for decorrelation across areas with active
!  parameterized convection and those without.
!
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  8/7/98
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nscalar
  INTEGER :: ccatopt
  REAL :: zs    (nx,ny,nz)  ! The physical height coordinate defined at
                            ! w-point of the staggered grid.
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate
  INTEGER :: icatg(nx,ny)
!
  REAL :: prctrw,bltop,rhrkul,thevrkul
  PARAMETER (prctrw=1.0E-04,    & !  (kg/(m**2*s))
         bltop=1000.,   & !  m
         rhrkul=0.8,                                                    &
             thevrkul=9.0)  !  K**2
!
  INTEGER :: i,j
  INTEGER :: knt1,knt2a,knt2b,knt3,knt4
  REAL :: pr,temp,qvsat,rh
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!dir$ inline always f_qvsat

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
!
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  knt1=0
  knt2a=0
  knt2b=0
  knt3=0
  knt4=0
!
!-----------------------------------------------------------------------
!
!  Initialize with default of no precip or precip influence
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      icatg(i,j)=1
    END DO
  END DO
!
  IF(ccatopt == 1) THEN
    DO j=1,ny-1
      DO i=1,nx-1
!
!-----------------------------------------------------------------------
!
!  Is convective precip turned on or heavy rain occuring?
!
!-----------------------------------------------------------------------
!
        IF (prcrate(i,j,3) > 0. .OR. prcrate(i,j,1) > prctrw) THEN
          knt4=knt4+1
          icatg(i,j)=4
!      print *, ' set i,j:',i,j,' cat 4',' prcrate(1):',
!    :             prcrate(i,j,1),' prcrate(3):',prcrate(i,j,3)
!
!-----------------------------------------------------------------------
!
!  Is it currently raining?
!
!-----------------------------------------------------------------------
!
        ELSE IF ( prcrate(i,j,1) > 0. ) THEN
          knt3=knt3+1
          icatg(i,j)=3
!      print *, ' set i,j:',i,j,' cat 3',' prcrate(1):',
!    :             prcrate(i,j,1),' prcrate(3):',prcrate(i,j,3)
!
!-----------------------------------------------------------------------
!
!  Evidence of rain-cooled air?
!  Lapse rate close to a moist adiabatic lapse rate in the lowest 2-km,
!  and a high relative humidity in that layer.  Or a high relative
!  humidity at the first level above ground.
!
!-----------------------------------------------------------------------
!
        ELSE
!
!-----------------------------------------------------------------------
!
!  For the layer 0-2km above first grid level (k=2).
!  Calculate the mean relative humidity.
!  Calculate the standard deviation of saturated equivalent
!  potential temperature.
!
!-----------------------------------------------------------------------
!
!      rhsum=0.
!      thesum=0.
!      thesq=0.
!      knt=0
!      DO 40 k=2,nz-1
!        IF(k.gt.2 .and. (zs(i,j,k)-zs(i,j,2)).gt.bltop) GO TO 41
!        knt=knt+1
!        pr=pprt(i,j,k)+pbar(i,j,k)
!        temp=(ptprt(i,j,k)+ptbar(i,j,k)) * (pr/p0)**rddcp
!        qvsat=f_qvsat( pr,temp )
!        rh=min(1.,(qv(i,j,k)/qvsat))
!        the=F_PT2PTE( pr,(ptprt(i,j,k)+ptbar(i,j,k)),qvsat)
!        rhsum=rhsum+rh
!        thesum=thesum+the
!        thesq=thesq+(the*the)
!  40     CONTINUE
!  41     CONTINUE
!      IF(knt.gt.0) THEN
!        rhmean=rhsum/float(knt)
!        themean=thesum/float(knt)
!      ELSE
!        rhmean=0.
!        themean=0.
!      END IF
!      IF(knt.gt.1) THEN
!        thevar=(thesq-((thesum*thesum)/float(knt)))/float(knt-1)
!      ELSE
!        thevar=9999.
!      END IF
!
!      IF (rhmean.gt.rhrkul .and. thevar.lt.thevrkul) THEN
!        knt2a=knt2a+1
!        icatg(i,j)=2
!        print *, ' set i,j:',i,j,' cat 2',' prcrate(1):',
!    :             prcrate(i,j,1),' prcrate(3):',prcrate(i,j,3)
!      ELSE
!
          pr=pprt(i,j,2)+pbar(i,j,2)
          temp=(ptprt(i,j,2)+ptbar(i,j,2)) * (pr/p0)**rddcp
          qvsat=f_qvsat( pr,temp )
          rh=qv(i,j,2)/qvsat
          IF(rh > 0.8) THEN
            knt2b=knt2b+1
            icatg(i,j)=2
!          print *, ' set i,j:',i,j,' cat 2*',' RH(2) = ',rh
          END IF

!      END IF
!      IF (mod(i,20).eq.0 .and. mod(j,10).eq.0)  THEN
!        print *, i,j,' rhmean=',rhmean,
!    :               ' thevar=',thevar,' icat=',icatg(i,j)
!        print *, ' '
!      END IF

        END IF
      END DO
    END DO
    knt1=((nx-1)*(ny-1))-knt2a-knt2b-knt3-knt4
    WRITE(6,'(a)') ' Precip correlation category assignment counts'
    WRITE(6,'(a,i6)') '     cat  1: ',knt1
    WRITE(6,'(a,i6)') '     cat 2a: ',knt2a
    WRITE(6,'(a,i6)') '     cat 2b: ',knt2b
    WRITE(6,'(a,i6)') '     cat  3: ',knt3
    WRITE(6,'(a,i6)') '     cat  4: ',knt4
  END IF
!
  CALL iedgfill(icatg,1,nx,1,nx-1, 1,ny,1,ny-1,1,1,1,1)
!
  RETURN
END SUBROUTINE setcat
