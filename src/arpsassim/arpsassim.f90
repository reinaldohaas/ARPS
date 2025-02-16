!
!##################################################################
!##################################################################
!######                                                      ######
!######      Advanced Regional Prediction System (ARPS)      ######
!######                     Version 4.4                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM ARPSASSIM
!
!-----------------------------------------------------------------------
!
!  CONTACT:
!
!  For further information, contact:
!
!  Center for Analysis and Prediction of Storms
!  University of Oklahoma
!  Sarkeys Energy Center,
!  100 East Boyd
!  Norman, OK  73019  USA
!  Phone:  (405) 325-0453
!  FAX  :  (405) 325-7614
!  E-mail: kdroege@tornado.caps.ou.edu
!      or: mxue@tornado.caps.ou.edu
!
!  User support: arpsuser@ou.edu
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This program is the driver for the Advanced Regional Prediction
!  System (ARPS) model.
!
!  The ARPS is a fully three-dimensional, compressible,
!  nonhydrostatic model based on the Arakawa C grid.  The model
!  solves three momentum equations, a thermodynamic energy
!  equation, a pressure equation, six moisture equations (water
!  vapor, cloud water, rainwater, cloud ice, hail, and snow),
!  and also provides for the statistical representation of
!  unresolvable processes via a turbulence parameterization scheme
!  (user option).  The model is written in a terrain-following
!  coordinate system, and has provisions for several boundary
!  condition options.
!
!  For a full description of ARPS, please consult the ARPS 4.0
!  User's Guide.
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY:
!
!  The modification history can be found in file HISTORY.
!
!-----------------------------------------------------------------------
!
!  GRID STRUCTURE AND STENCIL NUMBERING CONVENTION
!
!  Each of the three velocity components is located at a physical
!  boundary in the model.  The numbering convention for the
!  staggered grid is shown below for each of the coordinate
!  directions, where S denotes a scalar variable.
!
!  Arakawa C-grid is used by this model.
!
!
!  X-Direction (East.......West):
!
!        West boundary                       East boundary
!              |                                   |
!              <--------   Physical Domain   ------>
!  +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
!  U     S     U     S     U     S     U     S     U     S     U
!
! i=  1     1     2     2     . . .     NX-2  NX-2  NX-1  NX-1   NX
!
!
!  Y-Direction (North......South):
!
!         South boundary                     North boundary
!              |                                   |
!              <-------    Physical Domain  ------->
!  +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
!  V     S     V     S     V     S     V     S     V     S     V
!
! j=  1     1     2     2     . . .     NY-2  NY-2  NY-1  NY-1   NY
!
!
!  Z-Direction:
!
!                      NZ  + W
!                          |
!                          |
!                          |
!                     NZ-1 + S
!                          |
!                          |
!                          |
!                     NZ-1 + W --  Physical Top boundary.
!                          |
!                          .
!                          .
!                          .
!
!                          |
!                       2  + S
!                          |
!                          |
!                          |
!                       2  + W --  Physical Ground Level
!                          |
!                          |
!                          |
!                       1  + S
!                          |
!                          |
!                          |
!                 k=    1  + W
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE                ! Force explicit declarations
!
  INTEGER :: nt                ! The no. of time levels of time-dependent
                               ! arrays.
  INTEGER :: nx, ny, nz,nzsoil ! Dimension in x, y and z directions
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present time.
  INTEGER :: tfuture           ! Index of time level for the future time.

  INTEGER :: nxndg    ! Number of x grid points for nudging (1 or nx)
  INTEGER :: nyndg    ! Number of y grid points for nudging (1 or ny)
  INTEGER :: nzndg    ! Number of z grid points for nudging (1 or nz)

  INTEGER :: nstyps            ! Number of soil type

! PARAMETER (nstyps=4)
  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)
!
!-----------------------------------------------------------------------
!
!  Permanent arrays defined in the model.
!
!  These arrays must be saved (can not be overwritten) between
!  time steps.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Time dependent variables:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: u     (:,:,:,:)  ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v     (:,:,:,:)  ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w     (:,:,:,:)  ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: wcont (:,:,:)    ! Co:ravaria: vertical velocity (m/s)
  REAL, ALLOCATABLE :: ptprt (:,:,:,:)  ! Perturbation pote:ial temperature
                                        ! from that of base state atmosphere (K)
  REAL, ALLOCATABLE :: pprt  (:,:,:,:)  ! Perturbation pressure from that
                                        ! of base state atmosphere (Pascal)

  REAL, ALLOCATABLE :: qv    (:,:,:,:)  ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qc    (:,:,:,:)  ! Cloud water mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qr    (:,:,:,:)  ! Rain water mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qi    (:,:,:,:)  ! Cloud ice mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qs    (:,:,:,:)  ! Snow mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qh    (:,:,:,:)  ! Hail mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: tke   (:,:,:,:)  ! Turbule: Kinetic Energy ((m/s)**2)

  REAL, ALLOCATABLE :: pbldpth(:,:,:)   ! Planetary boundary layer depth (m)

  REAL, ALLOCATABLE :: kmh   (:,:,:)    ! Horizo:al turb. mixing coef. for
                                        ! mome:um. ( m**2/s )
  REAL, ALLOCATABLE :: kmv   (:,:,:)    ! Vertical turb. mixing coef. for
                                        ! mome:um. ( m**2/s )
  REAL, ALLOCATABLE :: rprntl(:,:,:)    ! Reciprocal of Prandtl number
!
!-----------------------------------------------------------------------
!
!  Time tendencies of certain time dependent variables at the lateral
!  boundaries, which are used to set the lateral boundary conditions
!  for options 4 and 5.
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: udteb (:,:)     ! Time tendency of u at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtwb (:,:)     ! Time tendency of u at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtnb (:,:)     ! Time tendency of u at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtsb (:,:)     ! Time tendency of u at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: vdteb (:,:)     ! Time tendency of v at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtwb (:,:)     ! Time tendency of v at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtnb (:,:)     ! Time tendency of v at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtsb (:,:)     ! Time tendency of v at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: wdteb (:,:)     ! Time tendency of w at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: wdtwb (:,:)     ! Time tendency of w at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: wdtnb (:,:)     ! Time tendency of w at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: wdtsb (:,:)     ! Time tendency of w at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: pdteb (:,:)     ! Time tendency of pprt at e-boundary (Pascal/s)
  REAL, ALLOCATABLE :: pdtwb (:,:)     ! Time tendency of pprt at w-boundary (Pascal/s)
  REAL, ALLOCATABLE :: pdtnb (:,:)     ! Time tendency of pprt at n-boundary (Pascal/s)
  REAL, ALLOCATABLE :: pdtsb (:,:)     ! Time tendency of pprt at s-boundary (Pascal/s)

  REAL, ALLOCATABLE :: sdteb (:,:)     ! Time tendency of a scalar at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: sdtwb (:,:)     ! Time tendency of a scalar at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: sdtnb (:,:)     ! Time tendency of a scalar at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: sdtsb (:,:)     ! Time tendency of a scalar at s-boundary (m/s**2)
!
!-----------------------------------------------------------------------
!
!  Base state variables:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: ubar  (:,:,:)     ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar  (:,:,:)     ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar (:,:,:)     ! Base state pote:ial temperature (K)
  REAL, ALLOCATABLE :: pbar  (:,:,:)     ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: ptbari(:,:,:)     ! Inverse Base state pot. temperature (K)
  REAL, ALLOCATABLE :: pbari (:,:,:)     ! Inverse Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhostr(:,:,:)     ! Base state density rhobar times j3.
  REAL, ALLOCATABLE :: rhostri(:,:,:)    ! Inv. base state density rhobar times j3.
  REAL, ALLOCATABLE :: qvbar (:,:,:)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL, ALLOCATABLE :: ppi   (:,:,:)     ! Exner function.
  REAL, ALLOCATABLE :: csndsq(:,:,:)     ! Sound wave speed squared.
!
!-----------------------------------------------------------------------
!
!  Arrays related to model grid definition:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)           ! The x-coord. of the physical and
                                           ! computational grid. Defined at u-poi:.
  REAL, ALLOCATABLE :: y     (:)           ! The y-coord. of the physical and
                                           ! computational grid. Defined at v-poi:.
  REAL, ALLOCATABLE :: z     (:)           ! The z-coord. of the computational grid.
                                           ! Defined at w-poi: on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:)       ! The physical height coordinate defined at
                                           ! w-poi: of the staggered grid.

  REAL, ALLOCATABLE :: zpsoil (:,:,:)      ! The physical height coordinate defined at
                                           ! w-poi: of the soil grid.

  REAL, ALLOCATABLE :: hterain(:,:)        ! The height of the terrain.

  REAL, ALLOCATABLE :: mapfct(:,:,:)       ! Map factors at scalar, u and v poi:s

  REAL, ALLOCATABLE :: j1    (:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as - d( zp )/d( x ).
  REAL, ALLOCATABLE :: j2    (:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as - d( zp )/d( y ).
  REAL, ALLOCATABLE :: j3    (:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as d( zp )/d( z ).

  REAL, ALLOCATABLE :: j3soil(:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as d( zpsoil )/d( z ).

  REAL, ALLOCATABLE :: aj3x  (:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL, ALLOCATABLE :: aj3y  (:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL, ALLOCATABLE :: aj3z  (:,:,:)       ! Coordinate transformation Jacobian defined
                                           ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL, ALLOCATABLE :: j3inv (:,:,:)       ! Inverse of J3.
  REAL, ALLOCATABLE :: j3soilinv (:,:,:)   ! Inverse of J3soil.


  REAL, ALLOCATABLE :: trigs1(:)           ! Array containing pre-computed trig
                                           ! function for fftopt=1.
  REAL, ALLOCATABLE :: trigs2(:)           ! Array containing pre-computed trig
                                           ! function for fftopt=1.
  INTEGER :: ifax1(13)                     ! Array containing the factors of nx for
                                           ! fftopt=1.
  INTEGER :: ifax2(13)                     ! Array containing the factors of ny for
                                           ! fftopt=1.

  REAL, ALLOCATABLE :: vwork1 (:,:)        ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: vwork2 (:,:)        ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: wsave1 (:)          ! Work array for fftopt =2.
  REAL, ALLOCATABLE :: wsave2 (:)          ! Work array for fftopt =2.

  REAL, ALLOCATABLE :: sinlat(:,:)         ! Sin of latitude at each grid point

  REAL, ALLOCATABLE :: ptcumsrc(:,:,:)     ! Source term in pt-equation due
                                           ! to cumulus parameterization
  REAL, ALLOCATABLE :: qcumsrc(:,:,:,:)    ! Source term in water equations due
                                           ! to cumulus parameterization:
                                           ! qcumsrc(1,1,1,1) for qv equation
                                           ! qcumsrc(1,1,1,2) for qc equation
                                           ! qcumsrc(1,1,1,3) for qr equation
                                           ! qcumsrc(1,1,1,4) for qi equation
                                           ! qcumsrc(1,1,1,5) for qs equation
  REAL, ALLOCATABLE :: w0avg(:,:,:)        ! a closing running average vertical
                                           ! velocity in 10min for K-F scheme
!
!-----------------------------------------------------------------------
!
!  Arrays for soil model.
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: soiltyp(:,:,:)   ! Soil type at each point
  REAL, ALLOCATABLE :: stypfrct(:,:,:)     ! Fraction of soil type
  INTEGER, ALLOCATABLE :: vegtyp (:,:)     ! Vegetation type at each point
  REAL, ALLOCATABLE :: lai    (:,:)        ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)        ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)        ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil (:,:,:,:)     ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil (:,:,:,:)     ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:)      ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)       ! Snow depth (m)
  REAL, ALLOCATABLE :: qvsfc  (:,:,:)      ! Effective sfc. qv
  REAL, ALLOCATABLE :: ptsfc  (:,:)        ! Ground surface potential temperature (K)

  REAL, ALLOCATABLE :: raing(:,:)          ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)          ! Cumulus convective rain

  REAL, ALLOCATABLE :: prcrate(:,:,:)      ! precipitation rates (kg/(m**2*s))
                                           ! prcrate(1,1,1) = total precip. rate
                                           ! prcrate(1,1,2) = grid scale precip. rate
                                           ! prcrate(1,1,3) = cumulus precip. rate
                                           ! prcrate(1,1,4) = microphysics precip. rate
  REAL, ALLOCATABLE :: raincv    (:,:)     ! K-F convective rainfall (cm)
  INTEGER, ALLOCATABLE :: nca    (:,:)     ! K-F counter for CAPE release
!
!-----------------------------------------------------------------------
!
!  Arrays for radiation
!
!-----------------------------------------------------------------------
!
  INCLUDE 'radcst.inc'                     ! dimension parameters for radiation
                                           ! working arrays and other radiation
                                           ! parameters

  INTEGER :: rbufsz

  REAL, ALLOCATABLE :: radfrc(:,:,:)       ! Radiation forcing (K/s)

  REAL, ALLOCATABLE :: radsw (:,:)         ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)         ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation

  REAL, ALLOCATABLE :: rad2d (:,:,:)       ! 2-D arrays for radiation (see radcst.inc)
                                           ! Buffur array to carry the variables calculated or used in
                                           ! radiation calculation. The last index defines the variables:
                                           ! 1  = nrsirbm,  Solar IR surface albedo for beam
                                           ! 2  = nrsirdf,  Solar IR surface albedo for diffuse
                                           ! 3  = nrsuvbm,  Solar UV surface albedo for beam
                                           ! 4  = nrsuvdf,  Solar UV surface albedo for diffuse
                                           ! 5  = ncosz,    Cosine of zenith
                                           ! 6  = ncosss,   Cosine of angle between sun light and
                                           !                terrain slope
                                           ! 7  = nfdirir,  all-sky direct downward IR flux
                                           !                (0.7-10 micron) at the surface
                                           ! 8  = nfdifir,  all-sky diffuse downward IR flux
                                           !                at the surface
                                           ! 9  = nfdirpar, all-sky direct downward par flux
                                           !                (0.4-0.7 micron) at the surface
                                           ! 10 = nfdifpar, all-sky diffuse downward par flux
                                           !                at the surface


  REAL, ALLOCATABLE :: radbuf(:)           ! Buffer to carry working arrays for
                                           ! radiation computing (see radcst.inc)
!
!-----------------------------------------------------------------------
!
!  Arrays that carry surface fluxes
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: usflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Buffer array that carrys external boundary 3-D arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: exbcbufsz
  REAL, ALLOCATABLE :: exbcbuf(:) ! EXBC buffer array

  REAL, ALLOCATABLE :: bcrlx(:,:)         ! EXBC relaxation function coefficients
!
!-----------------------------------------------------------------------
!
!  Work arrays that do not carry physical meaning in the code
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: temxy1(:,:)       ! Temporary work array, used as phydro

  REAL, ALLOCATABLE :: tem1  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem2  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem3  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem4  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem5  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem6  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem7  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem8  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem9  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem10 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem11 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem12 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem13 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem14 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem15 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem16 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem17 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem18 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem19 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem20 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem21 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem22 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem23 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem24 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem25 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem26 (:,:,:)     ! Temporary work array.

  REAL, ALLOCATABLE :: tem1_0(:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem2_0(:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem3_0(:,:,:)     ! Temporary work array.

  INTEGER :: frstep                      ! Flag for the initial time step

  INTEGER :: mptr                        ! Grid identifier.

!  integer ireturn                       ! Returned value from EXTBDT

  INTEGER :: ierr, i,j,k

  REAL :: cpu0, cpu1, f_cputime, tt

  REAL :: eps                            ! Small value to compensate the roundoff
                                         ! error in checking of the end of time
                                         ! integration.
  DATA eps /0.01/
!
  REAL :: testim1, testim2
  INTEGER :: skip1,skip2, istatus

  CHARACTER(LEN=256) :: namelist_filename
!
!=======================================================================
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'soilcst.inc'     ! Soil-vegetation parameters
  INCLUDE 'exbc.inc'        ! EXBC parameters
  INCLUDE 'assim.inc'       ! EXBC parameters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!MP insert:      call setupmp ()                  !Message passing code.
!MP insert:      call checkmpmx (nx, ny, nz)      !Message passing code.

  cpu1 = f_cputime()

  WRITE(6,'(/ 16(/5x,a)//)')                                            &
     '###############################################################', &
     '###############################################################', &
     '#####                                                     #####', &
     '#####                      Welcome to                     #####', &
     '#####                                                     #####', &
     '#####   The Advanced Regional Prediction System  (ARPS)   #####', &
     '#####                                                     #####', &
     '#####                     Version 4.0                     #####', &
     '#####                                                     #####', &
     '#####                     Developed by                    #####', &
     '#####     Center for Analysis and Prediction of Storms    #####', &
     '#####                University of Oklahoma               #####', &
     '#####                                                     #####', &
     '###############################################################', &
     '###############################################################'

  mptr = 1    ! Set the grid number to 1 for a single grid run.
  nestgrd = 0 ! Set the grid nesting flag to zero, i.e. no nesting.
!
!-----------------------------------------------------------------------
!
!  First, initialize the model parameters.  Most of them are
!  global parameters passed among subroutines through common blocks.
!
!  These parameters are declared in the include files:
!  globcst.inc, bndry.inc, exbc.inc, radcst.inc etc.
!
!-----------------------------------------------------------------------
!
!  CALL initpara(nx,ny,nz,nstyps)
  namelist_filename = ' '
  CALL initpara(nx,ny,nz,nzsoil,nstyps,namelist_filename)

!-----------------------------------------------------------------------
!
!  ALLOCATE the variable and Initialize to zero.
!
!-----------------------------------------------------------------------

  ALLOCATE(u(nx,ny,nz,nt),STAT=istatus)
  u=0
  ALLOCATE(v(nx,ny,nz,nt),STAT=istatus)
  v=0
  ALLOCATE(w(nx,ny,nz,nt),STAT=istatus)
  w=0
  ALLOCATE(wcont(nx,ny,nz),STAT=istatus)
  wcont=0
  ALLOCATE(ptprt(nx,ny,nz,nt),STAT=istatus)
  ptprt=0
  ALLOCATE(pprt(nx,ny,nz,nt),STAT=istatus)
  pprt=0
  ALLOCATE(qv(nx,ny,nz,nt),STAT=istatus)
  qv=0
  ALLOCATE(qc(nx,ny,nz,nt),STAT=istatus)
  qc=0
  ALLOCATE(qr(nx,ny,nz,nt),STAT=istatus)
  qr=0
  ALLOCATE(qi(nx,ny,nz,nt),STAT=istatus)
  qi=0
  ALLOCATE(qs(nx,ny,nz,nt),STAT=istatus)
  qs=0
  ALLOCATE(qh(nx,ny,nz,nt),STAT=istatus)
  qh=0
  ALLOCATE(tke(nx,ny,nz,nt),STAT=istatus)
  tke=0
  ALLOCATE(pbldpth(nx,ny,nt),STAT=istatus)
  pbldpth=0
  ALLOCATE(kmh(nx,ny,nz),STAT=istatus)
  kmh=0
  ALLOCATE(kmv(nx,ny,nz),STAT=istatus)
  kmv=0
  ALLOCATE(rprntl(nx,ny,nz),STAT=istatus)
  rprntl=0
  ALLOCATE(udteb(ny,nz),STAT=istatus)
  udteb=0
  ALLOCATE(udtwb(ny,nz),STAT=istatus)
  udtwb=0
  ALLOCATE(udtnb(nx,nz),STAT=istatus)
  udtnb=0
  ALLOCATE(udtsb(nx,nz),STAT=istatus)
  udtsb=0
  ALLOCATE(vdteb(ny,nz),STAT=istatus)
  vdteb=0
  ALLOCATE(vdtwb(ny,nz),STAT=istatus)
  vdtwb=0
  ALLOCATE(vdtnb(nx,nz),STAT=istatus)
  vdtnb=0
  ALLOCATE(vdtsb(nx,nz),STAT=istatus)
  vdtsb=0
  ALLOCATE(wdteb(ny,nz),STAT=istatus)
  wdteb=0
  ALLOCATE(wdtwb(ny,nz),STAT=istatus)
  wdtwb=0
  ALLOCATE(wdtnb(nx,nz),STAT=istatus)
  wdtnb=0
  ALLOCATE(wdtsb(nx,nz),STAT=istatus)
  wdtsb=0
  ALLOCATE(pdteb(ny,nz),STAT=istatus)
  pdteb=0
  ALLOCATE(pdtwb(ny,nz),STAT=istatus)
  pdtwb=0
  ALLOCATE(pdtnb(nx,nz),STAT=istatus)
  pdtnb=0
  ALLOCATE(pdtsb(nx,nz),STAT=istatus)
  pdtsb=0
  ALLOCATE(sdteb(ny,nz),STAT=istatus)
  sdteb=0
  ALLOCATE(sdtwb(ny,nz),STAT=istatus)
  sdtwb=0
  ALLOCATE(sdtnb(nx,nz),STAT=istatus)
  sdtnb=0
  ALLOCATE(sdtsb(nx,nz),STAT=istatus)
  sdtsb=0
  ALLOCATE(ubar(nx,ny,nz),STAT=istatus)
  ubar=0
  ALLOCATE(vbar(nx,ny,nz),STAT=istatus)
  vbar=0
  ALLOCATE(ptbar(nx,ny,nz),STAT=istatus)
  ptbar=0
  ALLOCATE(pbar(nx,ny,nz),STAT=istatus)
  pbar=0
  ALLOCATE(ptbari(nx,ny,nz),STAT=istatus)
  ptbari=0
  ALLOCATE(pbari(nx,ny,nz),STAT=istatus)
  pbari=0
  ALLOCATE(rhostr(nx,ny,nz),STAT=istatus)
  rhostr=0
  ALLOCATE(rhostri(nx,ny,nz),STAT=istatus)
  rhostri=0
  ALLOCATE(qvbar(nx,ny,nz),STAT=istatus)
  qvbar=0
  ALLOCATE(ppi(nx,ny,nz),STAT=istatus)
  ppi=0
  ALLOCATE(csndsq(nx,ny,nz),STAT=istatus)
  csndsq=0
  ALLOCATE(x(nx),STAT=istatus)
  x=0
  ALLOCATE(y(ny),STAT=istatus)
  y=0
  ALLOCATE(z(nz),STAT=istatus)
  z=0
  ALLOCATE(zp(nx,ny,nz),STAT=istatus)
  zp=0
  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
  zpsoil=0
  ALLOCATE(hterain(nx,ny),STAT=istatus)
  hterain=0
  ALLOCATE(mapfct(nx,ny,8),STAT=istatus)
  mapfct=0
  ALLOCATE(j1(nx,ny,nz),STAT=istatus)
  j1=0
  ALLOCATE(j2(nx,ny,nz),STAT=istatus)
  j2=0
  ALLOCATE(j3(nx,ny,nz),STAT=istatus)
  j3=0
  ALLOCATE(j3soil(nx,ny,nzsoil),STAT=istatus)
  j3soil=0
  ALLOCATE(aj3x(nx,ny,nz),STAT=istatus)
  aj3x=0
  ALLOCATE(aj3y(nx,ny,nz),STAT=istatus)
  aj3y=0
  ALLOCATE(aj3z(nx,ny,nz),STAT=istatus)
  aj3z=0
  ALLOCATE(j3inv(nx,ny,nz),STAT=istatus)
  j3inv=0
  ALLOCATE(j3soilinv(nx,ny,nzsoil),STAT=istatus)
  j3soilinv=0
  ALLOCATE(trigs1(3*(nx-1)/2+1),STAT=istatus)
  trigs1=0
  ALLOCATE(trigs2(3*(ny-1)/2+1),STAT=istatus)
  trigs2=0
  ALLOCATE(vwork1(nx+1,ny+1),STAT=istatus)
  vwork1=0
  ALLOCATE(vwork2(ny,nx+1),STAT=istatus)
  vwork2=0
  ALLOCATE(wsave1(3*(ny-1)+15),STAT=istatus)
  wsave1=0
  ALLOCATE(wsave2(3*(nx-1)+15),STAT=istatus)
  wsave2=0
  ALLOCATE(sinlat(nx,ny),STAT=istatus)
  sinlat=0
  ALLOCATE(ptcumsrc(nx,ny,nz),STAT=istatus)
  ptcumsrc=0
  ALLOCATE(qcumsrc(nx,ny,nz,5),STAT=istatus)
  qcumsrc=0
  ALLOCATE(w0avg(nx,ny,nz),STAT=istatus)
  w0avg=0
  ALLOCATE(soiltyp(nx,ny,nstyps),STAT=istatus)
  soiltyp=0
  ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
  stypfrct=0
  ALLOCATE(vegtyp(nx,ny),STAT=istatus)
  vegtyp=0
  ALLOCATE(lai(nx,ny),STAT=istatus)
  lai=0
  ALLOCATE(roufns(nx,ny),STAT=istatus)
  roufns=0
  ALLOCATE(veg(nx,ny),STAT=istatus)
  veg=0
  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  tsoil=0
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  qsoil=0
  ALLOCATE(wetcanp(nx,ny,0:nstyps),STAT=istatus)
  wetcanp=0
  ALLOCATE(snowdpth(nx,ny),STAT=istatus)
  snowdpth=0
  ALLOCATE(qvsfc(nx,ny,0:nstyps),STAT=istatus)
  qvsfc=0
  ALLOCATE(ptsfc(nx,ny),STAT=istatus)
  ptsfc=0
  ALLOCATE(raing(nx,ny),STAT=istatus)
  raing=0
  ALLOCATE(rainc(nx,ny),STAT=istatus)
  rainc=0
  ALLOCATE(prcrate(nx,ny,4),STAT=istatus)
  prcrate=0
  ALLOCATE(raincv(nx,ny),STAT=istatus)
  raincv=0
  ALLOCATE(nca(nx,ny),STAT=istatus)
  nca=0
  ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)
  radfrc=0
  ALLOCATE(radsw(nx,ny),STAT=istatus)
  radsw=0
  ALLOCATE(rnflx(nx,ny),STAT=istatus)
  rnflx=0
  ALLOCATE(radswnet (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radswnet")
  radswnet = 0
  ALLOCATE(radlwin (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radlwin")
  radlwin = 0

  ALLOCATE(usflx(nx,ny),STAT=istatus)
  usflx=0
  ALLOCATE(vsflx(nx,ny),STAT=istatus)
  vsflx=0
  ALLOCATE(ptsflx(nx,ny),STAT=istatus)
  ptsflx=0
  ALLOCATE(qvsflx(nx,ny),STAT=istatus)
  qvsflx=0
  ALLOCATE(bcrlx(nx,ny),STAT=istatus)
  bcrlx=0
  ALLOCATE(temxy1(nx,ny),STAT=istatus)
  temxy1=0
  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  tem1=0
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  tem2=0
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  tem3=0
  ALLOCATE(tem4(nx,ny,nz),STAT=istatus)
  tem4=0
  ALLOCATE(tem5(nx,ny,nz),STAT=istatus)
  tem5=0
  ALLOCATE(tem6(nx,ny,nz),STAT=istatus)
  tem6=0
  ALLOCATE(tem7(nx,ny,nz),STAT=istatus)
  tem7=0
  ALLOCATE(tem8(nx,ny,nz),STAT=istatus)
  tem8=0
  ALLOCATE(tem9(nx,ny,nz),STAT=istatus)
  tem9=0
  ALLOCATE(tem10(nx,ny,nz),STAT=istatus)
  tem10=0
  ALLOCATE(tem11(nx,ny,nz),STAT=istatus)
  tem11=0
  ALLOCATE(tem12(nx,ny,nz),STAT=istatus)
  tem12=0
  ALLOCATE(tem13(nx,ny,nz),STAT=istatus)
  tem13=0
  ALLOCATE(tem14(nx,ny,nz),STAT=istatus)
  tem14=0
  ALLOCATE(tem15(nx,ny,nz),STAT=istatus)
  tem15=0
  ALLOCATE(tem16(nx,ny,nz),STAT=istatus)
  tem16=0
  ALLOCATE(tem17(nx,ny,nz),STAT=istatus)
  tem17=0
  ALLOCATE(tem18(nx,ny,nz),STAT=istatus)
  tem18=0
  ALLOCATE(tem19(nx,ny,nz),STAT=istatus)
  tem19=0
  ALLOCATE(tem20(nx,ny,nz),STAT=istatus)
  tem20=0
  ALLOCATE(tem21(nx,ny,nz),STAT=istatus)
  tem21=0
  ALLOCATE(tem22(nx,ny,nz),STAT=istatus)
  tem22=0
  ALLOCATE(tem23(nx,ny,nz),STAT=istatus)
  tem23=0
  ALLOCATE(tem24(nx,ny,nz),STAT=istatus)
  tem24=0
  ALLOCATE(tem25(nx,ny,nz),STAT=istatus)
  tem25=0
  ALLOCATE(tem26(nx,ny,nz),STAT=istatus)
  tem26=0
  ALLOCATE(tem1_0(0:nx,0:ny,0:nz),STAT=istatus)
  tem1_0=0
  ALLOCATE(tem2_0(0:nx,0:ny,0:nz),STAT=istatus)
  tem2_0=0
  ALLOCATE(tem3_0(0:nx,0:ny,0:nz),STAT=istatus)
  tem3_0=0
!-----------------------------------------------------------------------


!  Set radiation and external boundary work array buffer sizes:

  IF (radopt == 2) THEN
    rbufsz = n2d_radiat*nx*ny + n3d_radiat*nx*ny*nz
  ELSE
    rbufsz = 1
  END IF

  ALLOCATE(radbuf(rbufsz),STAT=istatus)
  ALLOCATE(rad2d (nx,ny,nrad2d),STAT=istatus)

  IF (lbcopt == 2) THEN
    exbcbufsz = 22*nx*ny*nz
  ELSE
    exbcbufsz = 1
  END IF

  ALLOCATE(exbcbuf( exbcbufsz ),STAT=istatus)

!
  CALL initassim(nx,ny,nz)

  IF ( lbcopt == 2 ) CALL setexbcptr(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Initialize all arrays and variables, and set initial
!  condition for thermal disturbance if desired.
!
!-----------------------------------------------------------------------
!

  CALL initial(mptr,nx,ny,nz,nzsoil,nstyps, exbcbufsz,                  &
               u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,            &
               udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,         &
               wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,         &
               sdteb,sdtwb,sdtnb,sdtsb,                                 &
               ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,rhostri,        &
               qvbar,ppi,csndsq,                                        &
               x,y,z,zp,zpsoil,hterain,mapfct,                          &
               j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,          &
               trigs1,trigs2,ifax1,ifax2,                               &
               wsave1,wsave2,vwork1,vwork2,                             &
               sinlat,kmh,kmv,rprntl,                                   &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,ptsfc,qvsfc,                &
               ptcumsrc,qcumsrc,w0avg,nca,raincv,                       &
               raing,rainc,prcrate, exbcbuf,bcrlx,                      &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,temxy1,                        &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11,tem12,tem13,                       &
               tem14,tem15,tem16,tem17,tem18,tem19,                     &
               tem20,tem21,tem22,tem23,tem24,tem25,tem26,               &
               tem1_0,tem2_0,tem3_0)

!MP insert:      call mpbarrier ()                !Message passing code.
!
!  open(unit=83,file='voltim.input',form='formatted',status='old')
!  read (83,*) voltim1
!  read (83,*) voltim2
!  read (83,*) voltim3
!  close(unit=83)
!
  testim1 = voltim1+dtbig
  testim2 = voltim2+dtbig
  skip1 = 0
  skip2 = 0
!
  WRITE(6,*) '===> voltim1,2,3 = ',voltim1,voltim2,voltim3
  WRITE(6,*) '===> testim1,2 = ',testim1,testim2
!
!=======================================================================
!
!
!-----------------------------------------------------------------------
!
!  The current model time (initial time):
!
!-----------------------------------------------------------------------
!
  curtim = tstart

  WRITE(6,'(1x,a,f13.3,a)')                                             &
       'The initial model time is at ', curtim,' seconds.'
!
!-----------------------------------------------------------------------
!
!  Output the initial fields
!
!-----------------------------------------------------------------------
!

  CALL initout(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,                   &
               u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,            &
               ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,               &
               x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3inv,           &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,exbcbuf,                             &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11)

!MP insert:      call mpbarrier ()                !Message passing code.

  nstep = nint( curtim/dtbig )
!
!-----------------------------------------------------------------------
!
!  Time integration loop begins  ----------------------------------->
!
!-----------------------------------------------------------------------
!

  900   CONTINUE      ! Entry point of time step integration.

!
!-----------------------------------------------------------------------
!
!  Define the time step counter nstep reference to time zero
!  (curtim=0.0). This means for tstart.ne.0 case, nstep for the
!  first step of integration is not 1.
!
!-----------------------------------------------------------------------
!
  nstep = nstep +  1

  IF( (ABS(curtim-tstart) <= 1.0E-10) .AND. (restrt /= 1) ) THEN

    frstep=1           ! Indicate that this is the initial step of
                       ! model integration. For this step, a forward
                       ! time integration scheme will be used.

  ELSE                 ! For non-first step or restart run

    frstep=0           ! When frstep=0, leapfrog scheme is used.

  END IF

!
!-----------------------------------------------------------------------
!
!  Perform one time step integration for all equations:
!  On exit of this routine, all time dependent fields are
!  advanced by one time step.
!
!-----------------------------------------------------------------------
!
  CALL cordintg(mptr,frstep, nx,ny,nz,nzsoil,     &
                nxndg,nyndg,nzndg,                &
                rbufsz,nstyps,exbcbufsz,          &
                u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,               &
                tke,pbldpth,                                            &
                udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,        &
                wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,        &
                sdteb,sdtwb,sdtnb,sdtsb,                                &
                ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,rhostri,       &
                qvbar,ppi,csndsq,                                       &
                x,y,z,zp,zpsoil, mapfct,                                &
                j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,         &
                trigs1,trigs2,ifax1,ifax2,                              &
                wsave1,wsave2,vwork1,vwork2,                            &
                sinlat, kmh,kmv,rprntl,                                 &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,                 &
                tsoil,qsoil,wetcanp,snowdpth,               &
                ptsfc,qvsfc,                                            &
                ptcumsrc,qcumsrc,raing,rainc,prcrate,                   &
                w0avg,nca,raincv,                                       &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                rad2d,radbuf,exbcbuf,bcrlx,                             &
                usflx,vsflx,ptsflx,qvsflx,temxy1,                       &
                tem1,tem2,tem3,tem4,tem5,tem6,tem7,                     &
                tem8,tem9,tem10,tem11,tem12,tem13,                      &
                tem14,tem15,tem16,tem17,tem18,tem19,                    &
                tem20,tem21,tem22,tem23,tem24,tem25,tem26,              &
                tem1_0,tem2_0,tem3_0)

!
!-----------------------------------------------------------------------
!
!  Update physical time and generate output files
!
!-----------------------------------------------------------------------
!
  curtim = curtim + dtbig
!
!-----------------------------------------------------------------------
!
!  Print timestep marker
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a,i8,a,f10.2,a)')                                           &
      '  End of integration time step ', nstep,                         &
      ', Model time=',curtim,' (S)'

  CALL output(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,                    &
             u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,              &
             udteb,udtwb,vdtnb,vdtsb,pdteb,pdtwb,pdtnb,pdtsb,           &
             ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,                 &
             x,y,z,zp,zpsoil,hterain, mapfct, j1,j2,j3,j3soil,j3inv,    &
             j3soilinv,soiltyp,stypfrct,vegtyp,lai,roufns,veg,          &
             tsoil,qsoil,wetcanp,snowdpth,qvsfc,                        &
             ptcumsrc,qcumsrc,w0avg,nca,raincv,                         &
             raing,rainc,prcrate,exbcbuf,                               &
             radfrc,radsw,rnflx,radswnet,radlwin,                       &
             usflx,vsflx,ptsflx,qvsflx,                                 &
             tem1, tem2,tem3,tem4,tem5, tem6,tem7,                      &
             tem8,tem9,tem10,tem11)

!
!-----------------------------------------------------------------------
!
!  Check the stability of the integration.  If the model solution
!  appears to be exceeding the linear stability condition, then
!  perform a data dump for post-mortem.
!
!-----------------------------------------------------------------------
!
  cpu0 = f_cputime()

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem11(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  CALL chkstab(mptr,nx,ny,nz,nzsoil,nstyps,                             &
               u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,            &
               ubar,vbar,ptbar,pbar,tem11,qvbar,kmh,kmv,                &
               x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3soil,          &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  Calculate and apply the adjustment of domain translation
!  speed using either cell-tracking or optimal mean speed algorithm.
!
!-----------------------------------------------------------------------
!
  CALL grdtran(nx,ny,nz,ubar,vbar,u,v,w,ptprt,pprt,                     &
               qv,qc,qr,qi,qs,qh,qvbar,rhostr,x,y,zp,j3,j3inv,          &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11)

! misc_cpu_time = misc_cpu_time + f_cputime() - cpu0
!
!-----------------------------------------------------------------------
!
!  Time integration loop ends, go back to the beginning of the
!  time step loop if stop time is not reached.
!
!-----------------------------------------------------------------------
!
  IF( curtim >= testim1 .AND. skip1 == 0 ) THEN
    skip1 = 1
    nstep = nint(voltim2/dtbig)
    curtim = FLOAT(nstep) * dtbig
  ELSE IF( curtim >= testim2 .AND. skip2 == 0 ) THEN
    skip2 = 1
    nstep = nint(voltim3/dtbig)
    curtim = FLOAT(nstep) * dtbig
  END IF
!
!=======================================================================
!

  IF( curtim+eps*dtbig < tstop ) GO TO 900

!
!-----------------------------------------------------------------------
!
!  End of entire model time integration. The program stops.
!
!  Close opened files.
!
!-----------------------------------------------------------------------
!

  IF (hdmpfmt == 5) THEN

    CALL mclosescheme (gridid, ierr)
    CALL mclosedataset (dsindex, ierr)

  ELSE IF (hdmpfmt == 9) THEN

    CLOSE (nchdmp)

  END IF

  WRITE(6,'(//1x,a,/1x,a,f13.3,a,/1x,a)')                               &
      'ARPS stopped normally in the main program. ',                    &
      'The ending time was ', curtim,' seconds.',                       &
      'Thanks for using ARPS.'

! total_cpu_time = total_cpu_time + f_cputime() - cpu1

! IF(mphyopt == 1) warmrain_cpu_time = warmrain_cpu_time                &
!                  - qfall_cpu_time
! IF(mphyopt == 2) linice_cpu_time   = linice_cpu_time                  &
!                  - qfall_cpu_time
! IF(mphyopt == 3) nemice_cpu_time   = nemice_cpu_time                  &
!                  - qfall_cpu_time
!
!-----------------------------------------------------------------------
!
!  Sum up the CPU times used by all processors on MPP machines.
!
!-----------------------------------------------------------------------
!
!tmp samex
!cMP insert:      call mptotal(total_cpu_time)    !Message passing code.
!cMP insert:      call mptotal(init_cpu_time)     !Message passing code.
!cMP insert:      call mptotal(output_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(advuvw_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(advs_cpu_time)     !Message passing code.
!cMP insert:      call mptotal(coriol_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(smlstp_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(rad_cpu_time)      !Message passing code.
!cMP insert:      call mptotal(soil_cpu_time)     !Message passing code.
!cMP insert:      call mptotal(sfcphy_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(tmix_cpu_time)     !Message passing code.
!cMP insert:      call mptotal(cmix_cpu_time)     !Message passing code.
!cMP insert:      call mptotal(raydmp_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(tkesrc_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(bc_cpu_time)       !Message passing code.
!cMP insert:      call mptotal(qpfgrd_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(kuocum_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(kfcum_cpu_time)    !Message passing code.
!cMP insert:      call mptotal(warmrain_cpu_time) !Message passing code.
!cMP insert:      call mptotal(linice_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(nemice_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(qfall_cpu_time)   !Message passing code.
!cMP insert:      call mptotal(misc_cpu_time)     !Message passing code.
!cMP insert:      call mptotal(buoy_cpu_time)     !Message passing code.
!
!-----------------------------------------------------------------------
!
!  Print out CPU usage by various processes.
!
!-----------------------------------------------------------------------
!
!MP insert:      if (myproc .eq. 0) then        !Message passing code.

! tt = total_cpu_time * 0.01

! WRITE(6,'(a/a,25(/1x,a,e14.6,a,6x,f6.2,a))')                          &
!     ' Process             CPU time used   Percentage',                &
!     '-----------------------------------------------',                &
!     'Initialization  :',init_cpu_time,    's',init_cpu_time/tt,  '%', &
!     'Data output     :',output_cpu_time,  's',output_cpu_time/tt,'%', &
!     'Wind advection  :',advuvw_cpu_time,  's',advuvw_cpu_time/tt,'%', &
!     'Scalar advection:',advs_cpu_time,    's',advs_cpu_time/tt,  '%', &
!     'Coriolis force  :',coriol_cpu_time,  's',coriol_cpu_time/tt,'%', &
!     'Buoyancy term   :',buoy_cpu_time,    's',buoy_cpu_time/tt,  '%', &
!     'Small time steps:',smlstp_cpu_time,  's',smlstp_cpu_time/tt,'%', &
!     'Radiation       :',rad_cpu_time,     's',rad_cpu_time/tt,   '%', &
!     'Soil model      :',soil_cpu_time,    's',soil_cpu_time/tt,  '%', &
!     'Surface physics :',sfcphy_cpu_time,  's',sfcphy_cpu_time/tt,'%', &
!     'Turbulence      :',tmix_cpu_time,    's',tmix_cpu_time/tt,  '%', &
!     'Comput. mixing  :',cmix_cpu_time,    's',cmix_cpu_time/tt,  '%', &
!     'Rayleigh damping:',raydmp_cpu_time,  's',raydmp_cpu_time/tt,'%', &
!     'TKE src terms   :',tkesrc_cpu_time,  's',tkesrc_cpu_time/tt,'%', &
!     'Bound.conditions:',bc_cpu_time,      's',bc_cpu_time/tt,    '%', &
!     'Gridscale precp.:',qpfgrd_cpu_time,  's',qpfgrd_cpu_time/tt,'%', &
!     'Kuo cumulus     :',kuocum_cpu_time,  's',kuocum_cpu_time/tt,'%', &
!     'Kain-Fritsch    :',kfcum_cpu_time,   's',kfcum_cpu_time/tt, '%', &
!     'Warmrain microph:',warmrain_cpu_time,'s',warmrain_cpu_time/tt,   &
!     '%',                                                              &
!     'Lin ice microph :',linice_cpu_time,  's',linice_cpu_time/tt,'%', &
!     'NEM ice microph :',nemice_cpu_time,  's',nemice_cpu_time/tt,'%', &
!     'Hydrometero fall:',qfall_cpu_time,   's',qfall_cpu_time/tt,'%',  &
!     'Miscellaneous   :',misc_cpu_time,    's',misc_cpu_time/tt,  '%'

! WRITE(6,'(/1x,a,e14.6,a,6x,f6.2,a)')                                  &
!     'Entire model    :',total_cpu_time,   's',total_cpu_time/tt, '%'

! WRITE(6,'(/1x,a,f6.2,a)')                                             &
!     'Total CPU time accounted for=',                                  &
!     (init_cpu_time+output_cpu_time+sfcphy_cpu_time+soil_cpu_time+     &
!     rad_cpu_time+                                                     &
!     tmix_cpu_time+cmix_cpu_time+raydmp_cpu_time+advuvw_cpu_time+      &
!     advs_cpu_time+coriol_cpu_time+tkesrc_cpu_time+bc_cpu_time+        &
!     qpfgrd_cpu_time+kuocum_cpu_time+warmrain_cpu_time+smlstp_cpu_time+ &
!     linice_cpu_time+nemice_cpu_time+misc_cpu_time+buoy_cpu_time+      &
!     qfall_cpu_time)/tt,                                               &
!     '%'

!MP insert:      ENDIF        !Message passing code.

  STOP

END PROGRAM ARPSASSIM
