!
!##################################################################
!##################################################################
!######                                                      ######
!######      Advanced Regional Prediction System (ARPS)      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM arps
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
!  Norman, OK  73019 USA
!  Phone:  (405) 325-0453
!  FAX  :  (405) 325-7614
!  E-mail: kkd@ou.edu or: mxue@ou.edu
!
!  User support: arpssupport@ou.edu
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
!                     NZ-1 + S
!                          |
!                          |
!                     NZ-1 + W --  Physical Top boundary.
!                          |
!                          .
!                          .
!
!                          |
!                       2  + S
!                          |
!                          |
!                       2  + W --  Physical Ground Level
!                          |
!                          |
!                       1  + S
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
  IMPLICIT NONE             ! Force explicit declarations
!

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'bndry.inc'
  INCLUDE 'timelvls.inc'
  INCLUDE 'radcst.inc'      ! Includes radiation grid size information.
  INCLUDE 'mp.inc'          ! Message passing parameters.
  INCLUDE 'alloc.inc'       ! Memory allocation declaration and interfaces
  INCLUDE 'nudging.inc'     ! Memory allocation declaration and interfaces
  INCLUDE 'dfilter.inc'     ! for digital filter part
  !wdt tmp buffer
  !INCLUDE "mpif.h"

!
!-----------------------------------------------------------------------
!
!  Time dependent variables:
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: u     (:,:,:,:)   ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v     (:,:,:,:)   ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w     (:,:,:,:)   ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: wcont (:,:,:)     ! Contravariant vertical velocity (m/s)
  REAL, ALLOCATABLE :: ptprt (:,:,:,:)   ! Perturbation potential temperature
                                         ! from that of base state atmosphere (K)
  REAL, ALLOCATABLE :: pprt  (:,:,:,:)   ! Perturbation pressure from that
                                         ! of base state atmosphere (Pascal)

  REAL, ALLOCATABLE :: qv     (:,:,:,:)  ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:,:)
  REAL, ALLOCATABLE :: tke    (:,:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL, ALLOCATABLE :: pbldpth(:,:,:)    ! Planetary boundary layer depth (m)

  REAL, ALLOCATABLE :: kmh   (:,:,:)     ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv   (:,:,:)     ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: rprntl(:,:,:)     ! Reciprocal of Prandtl number
!
!-----------------------------------------------------------------------
!
!  Time tendencies of certain time dependent variables at the lateral
!  boundaries, which are used to set the lateral boundary conditions
!  for options 4 and 5.
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: udteb (:,:) ! Time tendency of u at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtwb (:,:) ! Time tendency of u at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtnb (:,:) ! Time tendency of u at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtsb (:,:) ! Time tendency of u at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: vdteb (:,:) ! Time tendency of v at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtwb (:,:) ! Time tendency of v at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtnb (:,:) ! Time tendency of v at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtsb (:,:) ! Time tendency of v at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: wdteb (:,:) ! Time tendency of w at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: wdtwb (:,:) ! Time tendency of w at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: wdtnb (:,:) ! Time tendency of w at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: wdtsb (:,:) ! Time tendency of w at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: pdteb (:,:) ! Time tendency of pprt at e-boundary (Pascal/s)
  REAL, ALLOCATABLE :: pdtwb (:,:) ! Time tendency of pprt at w-boundary (Pascal/s)
  REAL, ALLOCATABLE :: pdtnb (:,:) ! Time tendency of pprt at n-boundary (Pascal/s)
  REAL, ALLOCATABLE :: pdtsb (:,:) ! Time tendency of pprt at s-boundary (Pascal/s)

  REAL, ALLOCATABLE :: sdteb (:,:) ! Time tendency of a scalar at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: sdtwb (:,:) ! Time tendency of a scalar at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: sdtnb (:,:) ! Time tendency of a scalar at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: sdtsb (:,:) ! Time tendency of a scalar at s-boundary (m/s**2)
!
!-----------------------------------------------------------------------
!
!  Base state variables:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: ubar  (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar  (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar  (:,:,:)   ! Base state w-velocity (m/s)  - Added with hxopt = 1
  REAL, ALLOCATABLE :: ptbar (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar  (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: ptbari(:,:,:) ! Inverse Base state pot. temperature (K)
  REAL, ALLOCATABLE :: pbari (:,:,:) ! Inverse Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhostr(:,:,:) ! Base state density rhobar times j3.
  REAL, ALLOCATABLE :: rhostri(:,:,:)! Inv. base state density rhobar times j3.
  REAL, ALLOCATABLE :: qvbar (:,:,:) ! Base state water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: ppi   (:,:,:) ! Exner function.
  REAL, ALLOCATABLE :: csndsq(:,:,:) ! Sound wave speed squared.
!
!-----------------------------------------------------------------------
!
!  Arrays related to model grid definition:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)     ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)     ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)     ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:) ! The physical height coordinate defined at
                                     ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil (:,:,:)! The physical height coordinate defined
                                     ! at the center of a soil layer(m).

  REAL, ALLOCATABLE :: hterain(:,:)  ! The height of the terrain.

  REAL, ALLOCATABLE :: mapfct(:,:,:) ! Map factors at scalar, u and v points

  REAL, ALLOCATABLE :: j1    (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as - d( zp )/d( x ).
  REAL, ALLOCATABLE :: j2    (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as - d( zp )/d( y ).
  REAL, ALLOCATABLE :: j3    (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as d( zp )/d( z ).
  REAL, ALLOCATABLE :: aj3x  (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL, ALLOCATABLE :: aj3y  (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL, ALLOCATABLE :: aj3z  (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL, ALLOCATABLE :: j3inv (:,:,:) ! Inverse of J3.

  REAL, ALLOCATABLE :: j3soil (:,:,:)! Coordinate transformation Jacobian
                                     ! defined as d( zpsoil )/d( zsoil ).
  REAL, ALLOCATABLE :: j3soilinv (:,:,:) ! Inverse of J3soil.

  REAL, ALLOCATABLE :: trigs1(:)     ! Array containing pre-computed trig
                                     ! function for fftopt=1.
  REAL, ALLOCATABLE :: trigs2(:)     ! Array containing pre-computed trig
                                     ! function for fftopt=1.
  INTEGER, ALLOCATABLE :: ifax1(:)   ! Array containing the factors of nx for
                                     ! fftopt=1.
  INTEGER, ALLOCATABLE :: ifax2(:)   ! Array containing the factors of ny for
                                     ! fftopt=1.

  REAL, ALLOCATABLE :: vwork1 (:,:)  ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: vwork2 (:,:)  ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: wsave1 (:)    ! Work array for fftopt =2.
  REAL, ALLOCATABLE :: wsave2 (:)    ! Work array for fftopt =2.

  REAL, ALLOCATABLE :: sinlat(:,:)   ! Sin of latitude at each grid point

  REAL, ALLOCATABLE :: ptcumsrc(:,:,:) ! Source term in pt-equation due
                                       ! to cumulus parameterization
  REAL, ALLOCATABLE :: qcumsrc(:,:,:,:)! Source term in water equations due
                                       ! to cumulus parameterization:
                                       ! qcumsrc(1,1,1,1) for qv equation
                                       ! qcumsrc(1,1,1,2) for qc equation
                                       ! qcumsrc(1,1,1,3) for qr equation
                                       ! qcumsrc(1,1,1,4) for qi equation
                                       ! qcumsrc(1,1,1,5) for qs equation
  REAL, ALLOCATABLE :: w0avg(:,:,:)    ! a closing running average vertical
                                       ! velocity in 10min for K-F scheme
!
!-----------------------------------------------------------------------
!
!  Arrays for soil model.
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: soiltyp(:,:,:)  ! Soil type at each point
  REAL,    ALLOCATABLE :: stypfrct(:,:,:) ! Fraction of soil type
  INTEGER, ALLOCATABLE :: vegtyp (:,:)    ! Vegetation type at each point
  REAL,    ALLOCATABLE :: lai    (:,:)    ! Leaf Area Index
  REAL,    ALLOCATABLE :: roufns (:,:)    ! Surface roughness
  REAL,    ALLOCATABLE :: veg    (:,:)    ! Vegetation fraction
  REAL,    ALLOCATABLE :: tsoil  (:,:,:,:)   ! soil temperature (K)
  REAL,    ALLOCATABLE :: qsoil  (:,:,:,:)   ! soil moisture (m**3/m**3)
  REAL,    ALLOCATABLE :: wetcanp(:,:,:)     ! Canopy water amount
  REAL,    ALLOCATABLE :: snowdpth(:,:)      ! Snow depth (m)
  REAL,    ALLOCATABLE :: qvsfc  (:,:,:)     ! Effective sfc. qv
  REAL,    ALLOCATABLE :: ptsfc  (:,:)       ! Ground surface potential temperature (K)

  REAL,    ALLOCATABLE :: raing(:,:)         ! Grid supersaturation rain
  REAL,    ALLOCATABLE :: rainc(:,:)         ! Cumulus convective rain

  REAL,    ALLOCATABLE :: prcrate(:,:,:)     ! precipitation rates (kg/(m**2*s))
                                             ! prcrate(1,1,1) = total precip. rate
                                             ! prcrate(1,1,2) = grid scale precip. rate
                                             ! prcrate(1,1,3) = cumulative precip. rate
                                             ! prcrate(1,1,4) = microphysics precip. rate
  REAL,    ALLOCATABLE :: kfraincv(:,:)      ! K-F convective rainfall (cm)
  INTEGER, ALLOCATABLE :: nca(:,:)           ! K-F counter for CAPE release

  REAL, ALLOCATABLE :: cldefi(:,:)    ! BMJ cloud efficiency
  REAL, ALLOCATABLE :: xland(:,:)     ! BMJ land mask
                                      ! (1.0 = land, 2.0 = sea)
  REAL, ALLOCATABLE :: bmjraincv(:,:) ! BMJ convective rainfall (cm)

!
!-----------------------------------------------------------------------
!
!  Arrays for radiation
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: radfrc(:,:,:) ! Radiation forcing (K/s)

  REAL, ALLOCATABLE :: radsw (:,:)   ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)   ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation


  REAL, ALLOCATABLE :: rad2d (:,:,:) ! 2-D arrays for radiation (see radcst.inc)

  REAL, ALLOCATABLE :: radbuf(:)     ! Buffer to carry working arrays for
                                     ! radiation computing (see radcst.inc)
!
!-----------------------------------------------------------------------
!
!  Arrays that carry surface fluxes
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: usflx (:,:)  ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)  ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)  ! Surface heat flux (K*kg/(m**2*s))
  REAL, ALLOCATABLE :: qvsflx(:,:)  ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Buffer arrays that carry external boundary 3-D arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: exbcbuf( : ) ! EXBC buffer array
  REAL, ALLOCATABLE :: bcrlx(:,:)   ! EXBC relaxation function coefficients
!
!-----------------------------------------------------------------------
!
!  Arrays for analysis increment updating (a type of nudging) output
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: uincr(:,:,:)
  REAL, ALLOCATABLE :: vincr(:,:,:)
  REAL, ALLOCATABLE :: wincr(:,:,:)
  REAL, ALLOCATABLE :: pincr(:,:,:)
  REAL, ALLOCATABLE :: ptincr(:,:,:)
  REAL, ALLOCATABLE :: qvincr(:,:,:)
  REAL, ALLOCATABLE :: qcincr(:,:,:)
  REAL, ALLOCATABLE :: qrincr(:,:,:)
  REAL, ALLOCATABLE :: qiincr(:,:,:)
  REAL, ALLOCATABLE :: qsincr(:,:,:)
  REAL, ALLOCATABLE :: qhincr(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Work arrays that do not carry physical meaning in the code
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: temxy1(:,:)     ! Temporary work array, used as phydro

  REAL, ALLOCATABLE :: tem1soil(:,:,:) ! Temporary soil model work array.
  REAL, ALLOCATABLE :: tem2soil(:,:,:) ! Temporary soil model work array.
  REAL, ALLOCATABLE :: tem3soil(:,:,:) ! Temporary soil model work array.
  REAL, ALLOCATABLE :: tem4soil(:,:,:) ! Temporary soil model work array.
  REAL, ALLOCATABLE :: tsdiffus(:,:,:) ! Temporary soil model work array.

  REAL, ALLOCATABLE :: tem1  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem2  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem3  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem4  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem5  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem6  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem7  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem8  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem9  (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem10 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem11 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem12 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem13 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem14 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem15 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem16 (:,:,:)   ! Temporary work array.

  REAL, ALLOCATABLE :: tem4d (:,:,:,:)

  REAL, ALLOCATABLE :: tem1_0(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem2_0(:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem3_0(:,:,:)   ! Temporary work array.

  !for dfilter parer
  REAL, ALLOCATABLE :: df_store_atm3d(:,:,:,:)    ! Temporary work array.
  REAL, ALLOCATABLE :: df_store_soil(:,:,:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: df_store_soil2d(:,:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: df_store_soil2ds(:,:,:)    ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  ARPS dimensions:
!
!  nx, ny, nz: Dimensions of computatioal grid. When run on MPP
!              with PVM or MPI, they represent of the size of the
!              subdomains. See below.
!
!              Given nx, ny and nz, the physical domain size will be
!              xl=(nx-3)*dx by yl=(ny-3)*dy by zh=(nz-3)*dz. The
!              variables nx, ny, nz, dx, dy and dz are read in from
!              the input file by subroutine INITPARA.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
  INTEGER :: nzsoil   ! Number of grid points in the -z-direction


  INTEGER :: nxndg    ! Number of x grid points for nudging (1 or nx)
  INTEGER :: nyndg    ! Number of y grid points for nudging (1 or ny)
  INTEGER :: nzndg    ! Number of z grid points for nudging (1 or nz)
!
!-----------------------------------------------------------------------
!
!  Soil types.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nstyps   ! Number of soil types
!
!-----------------------------------------------------------------------
!
!  Radiation buffer size.  Set equal to
!     n2d_radiat*nx*ny + n3d_radiat*nx*ny*nz
!  in subroutine inipara if radopt=2, otherwise set to 1.
!
!-----------------------------------------------------------------------
!
  INTEGER :: rbufsz
!
!-----------------------------------------------------------------------
!
!  External boundary buffer size.  Set equal to
!     22*nx*ny*nz
!  in subroutine inipara if lbcopt=2, otherwise set to 1.
!
!-----------------------------------------------------------------------
!
  INTEGER :: exbcbufsz



  INTEGER :: frstep            ! Flag for the initial time step
  INTEGER :: mptr              ! Grid identifier.
  INTEGER :: ierr, i,j,k
  REAL :: eps                  ! Small value to compensate the roundoff
                               ! error in checking of the end of time
                               ! integration.
  REAL(KIND=8) :: dcurtim      ! DTD: Added double-precision current time
                               ! to avoid roundoff errors when very small
                               ! timesteps are used.

  DATA eps /0.01/
  INTEGER istatus              ! , nxy, nxyz ! comment out by WYH

!  REAL :: twall0, t_walltime                ! comment out by WYH

  !wdt tmp buffer
  !INTEGER :: n_buffer
  !REAL, ALLOCATABLE :: mpi_buffer(:)

  CHARACTER(LEN=256) :: namelist_filename
!clt Misc. for dfilter
  integer :: df_dstp, df_ststp
  integer :: idfstp
  real    :: addwght
  real    :: df_curtim0
!
!-----------------------------------------------------------------------
!
!  Work varibles for Output H(X)
!
!-----------------------------------------------------------------------
!
  INTEGER :: tmstrln
  CHARACTER (LEN=10)  :: timsnd
  LOGICAL :: rad_file_state,sfc_file_state,snd_file_state,pro_file_state
  INTEGER :: runmode_obs = 0
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
!  First, initialize the model parameters.  Most of them are
!  global parameters passed among subroutines through common blocks.
!
!  These parameters are declared in the include files:
!  globcst.inc, bndry.inc, exbc.inc, radcst.inc etc.
!
!-----------------------------------------------------------------------
!

  CALL acct_init
  CALL set_acct(init_acct)

  namelist_filename = ' '
  CALL initpara(nx,ny,nz,nzsoil,nstyps,namelist_filename)


  ! Set radiation and external boundary work array buffer sizes:

  IF (radopt == 2) THEN
    rbufsz = n2d_radiat*nx*ny + n3d_radiat*nx*ny*nz
  ELSE
    rbufsz = 1
  END IF

  IF (lbcopt == 2) THEN
    exbcbufsz = (12+2*nscalar)*nx*ny*nz
  ELSE
    exbcbufsz = 1
  END IF

  IF ( nudgopt > 0 ) THEN
    nxndg=nx
    nyndg=ny
    nzndg=nz
  ELSE
    nxndg=1
    nyndg=1
    nzndg=1
  END IF
!
!-----------------------------------------------------------------------
!
!  Print a log file in namelist format.
!
!-----------------------------------------------------------------------
!
  CALL acct_interrupt(output_acct)
  IF (myproc == 0) THEN

    IF (mp_opt == 0) THEN
      CALL prtlog(nx,ny,nz,nzsoil,6)
      CALL prtlog(nx,ny,nz,nzsoil,0)
    ELSE
      CALL prtlog(nproc_x*(nx-3)+3,nproc_y*(ny-3)+3,nz,nzsoil,6)
      CALL prtlog(nproc_x*(nx-3)+3,nproc_y*(ny-3)+3,nz,nzsoil,0)
    ENDIF

  END IF
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
! Allocate arrays etc.
!
!-----------------------------------------------------------------------
!
  mptr = 1    ! Set the grid number to 1 for a single grid run.
  nestgrd = 0 ! Set the grid nesting flag to zero, i.e. no nesting.

!  nxy = nx*ny
!  nxyz = nxy*nz

  IF (myproc == 0) WRITE(6,*) 'Allocating arrays'

  !wdt tmp buffer
  !n_buffer = 100*(max(nx,ny)*nz+MPI_BSEND_OVERHEAD)
  !ALLOCATE(mpi_buffer(n_buffer))
  !CALL MPI_BUFFER_ATTACH(mpi_buffer,n_buffer*4,istatus)

  ALLOCATE(u    (nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:u")
  u = 0
  ALLOCATE(v    (nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:v")
  v = 0
  ALLOCATE(w    (nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:w")
  w = 0
  ALLOCATE(wcont(nx,ny,nz   ),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wcont")
  wcont = 0
  ALLOCATE(ptprt(nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptprt")
  ptprt = 0
  ALLOCATE(pprt (nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pprt")
  pprt = 0
  ALLOCATE(qv   (nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qv")
  qv = 0
  ALLOCATE(qscalar(nx,ny,nz,nt,nscalar), STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qscalar")
  qscalar(:,:,:,:,:) = 0.0
  ALLOCATE(tke  (nx,ny,nz,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tke")
  tke = 0

  ALLOCATE(pbldpth(nx,ny,nt),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pbldpth")
  pbldpth = 0
  ALLOCATE(kmh(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:kmh")
  kmh = 0
  ALLOCATE(kmv(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:kmv")
  kmv = 0
  ALLOCATE(rprntl(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rprntl")
  rprntl = 0

  ALLOCATE(udteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:udteb")
  udteb = 0
  ALLOCATE(udtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:udtwb")
  udtwb = 0
  ALLOCATE(udtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:udtnb")
  udtnb = 0
  ALLOCATE(udtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:udtsb")
  udtsb = 0

  ALLOCATE(vdteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vdteb")
  vdteb = 0
  ALLOCATE(vdtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vdtwb")
  vdtwb = 0
  ALLOCATE(vdtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vdtnb")
  vdtnb = 0
  ALLOCATE(vdtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vdtsb")
  vdtsb = 0

  ALLOCATE(wdteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wdteb")
  wdteb = 0
  ALLOCATE(wdtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wdtwb")
  wdtwb = 0
  ALLOCATE(wdtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wdtnb")
  wdtnb = 0
  ALLOCATE(wdtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wdtsb")
  wdtsb = 0

  ALLOCATE(pdteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pdteb")
  pdteb = 0
  ALLOCATE(pdtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pdteb")
  pdtwb = 0
  ALLOCATE(pdtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pdtnb")
  pdtnb = 0
  ALLOCATE(pdtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pdtsb")
  pdtsb = 0

  ALLOCATE(sdteb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:sdteb")
  sdteb = 0
  ALLOCATE(sdtwb(ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:sdtwb")
  sdtwb = 0
  ALLOCATE(sdtnb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:sdtnb")
  sdtnb = 0
  ALLOCATE(sdtsb(nx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:sdtsb")
  sdtsb = 0

  ALLOCATE(ubar  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ubar")
  ubar = 0.0
  ALLOCATE(vbar  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vbar")
  vbar = 0.0
  ALLOCATE(wbar  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wbar")
  wbar = 0.0

  ALLOCATE(ptbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptbar")
  ptbar = 0
  ALLOCATE(pbar  (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pbar")
  pbar = 0
  ALLOCATE(ptbari(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptbari")
  ptbari = 0
  ALLOCATE(pbari (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pbari")
  pbari = 0
  ALLOCATE(rhostr(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rhostr")
  rhostr = 0
  ALLOCATE(rhostri(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rhostri")
  rhostri = 0
  ALLOCATE(qvbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvbar")
  qvbar = 0

  ALLOCATE(ppi   (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ppi")
  ppi = 0
  ALLOCATE(csndsq(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:csndsq")
  csndsq = 0

  ALLOCATE(x(nx),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:x")
  x = 0
  ALLOCATE(y(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:y")
  y = 0
  ALLOCATE(z(nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:z")
  z = 0

  ALLOCATE(zp(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:zp")
  zp = 0

  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:zpsoil")
  zpsoil = 0


  ALLOCATE(hterain(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:hterain")
  hterain = 0

  ALLOCATE(mapfct(nx,ny,8),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:mapfct")
  mapfct = 0

  ALLOCATE(j1   (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:j1")
  j1 = 0
  ALLOCATE(j2   (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:j2")
  j2 = 0
  ALLOCATE(j3   (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:j3")
  j3 = 0.0

  ALLOCATE(j3soil (nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:j3soil")
  j3soil = 0.0

  ALLOCATE(aj3x (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:aj3x")
  aj3x = 0
  ALLOCATE(aj3y (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:aj3y")
  aj3y = 0
  ALLOCATE(aj3z (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:aj3z")
  aj3z = 0
  ALLOCATE(j3inv(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:j3inv")
  j3inv = 0

  ALLOCATE(j3soilinv(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:j3soilinv")
  j3soilinv = 0

  ALLOCATE(trigs1(3*(nx-1)/2+1),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:trigs1")
  trigs1 = 0
  ALLOCATE(trigs2(3*(ny-1)/2+1),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:trigs2")
  trigs2 = 0


  ALLOCATE(ifax1(13),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ifax1")
  ifax1 = 0
  ALLOCATE(ifax2(13),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ifax2")
  ifax2 = 0

  ALLOCATE(vwork1(nx+1,ny+1),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vwork1")
  vwork1 = 0
  ALLOCATE(vwork2(ny,nx+1),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vwork2")
  vwork2 = 0

  ALLOCATE(wsave1(3*(ny-1)+15),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wsave1")
  wsave1 = 0
  ALLOCATE(wsave2(3*(nx-1)+15),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wsave2")
  wsave2 = 0

  ALLOCATE(sinlat(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:sinlat")
  sinlat = 0

  ALLOCATE(ptcumsrc(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptcumsrc")
  ptcumsrc = 0
  ALLOCATE(qcumsrc(nx,ny,nz,5),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qcumsrc")
  qcumsrc = 0

  ALLOCATE(w0avg(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:w0avg")
  w0avg = 0

  ALLOCATE(soiltyp (nx,ny,nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:soiltyp")
  soiltyp = 0
  ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:stypfrct")
  stypfrct = 0

  ALLOCATE(vegtyp(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vegtyp")
  vegtyp = 0
  ALLOCATE(lai   (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:lai")
  lai = 0
  ALLOCATE(roufns(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:roufns")
  roufns = 0
  ALLOCATE(veg   (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:veg")
  veg = 0

  ALLOCATE(tsoil   (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tsoil")
  tsoil = 0
  ALLOCATE(qsoil   (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qsoil")
  qsoil = 0
  ALLOCATE(wetcanp (nx,ny,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wetcanp")
  wetcanp = 0
  ALLOCATE(qvsfc   (nx,ny,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvsfc")
  qvsfc = 0

  ALLOCATE(ptsfc   (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptsfc")
  ptsfc = 0
  ALLOCATE(snowdpth(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:snowdpth")
  snowdpth = 0

  ALLOCATE(raing  (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:raing")
  raing = 0
  ALLOCATE(rainc  (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rainc")
  rainc = 0
  ALLOCATE(kfraincv (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:kfraincv")
  kfraincv = 0
  ALLOCATE(nca    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:nca")
  nca = 0

  ALLOCATE(prcrate(nx,ny,4),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:prcrate")
  prcrate = 0

!EMK BMJ
  ALLOCATE(cldefi(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:cldefi")
  cldefi(:,:) = 1

  ALLOCATE(xland(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:xland")
  xland(:,:) = 0

  ALLOCATE(bmjraincv (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:bmjraincv")
  bmjraincv = 0

! EMK END

  ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radfrc")
  radfrc = 0
  ALLOCATE(radsw (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radsw")
  radsw = 0
  ALLOCATE(rnflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rnflx")
  rnflx = 0
  ALLOCATE(radswnet (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radswnet")
  radswnet = 0
  ALLOCATE(radlwin (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radlwin")
  radlwin = 0
  ALLOCATE(rad2d (nx,ny,nrad2d),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rad2d")
  rad2d = 0
  ALLOCATE(radbuf(rbufsz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radbuf")
  radbuf = 0

  ALLOCATE(usflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:usflx")
  usflx = 0
  ALLOCATE(vsflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vsflx")
  vsflx = 0
  ALLOCATE(ptsflx(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptsflx")
  ptsflx = 0
  ALLOCATE(qvsflx(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvsflx")
  qvsflx = 0

  ALLOCATE(exbcbuf(exbcbufsz ),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:exbcbuf")
  exbcbuf = 0
  ALLOCATE(bcrlx(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:bcrlx")
  bcrlx = 0

  IF (myproc == 0) WRITE(6,*) '  nxndg,nyndg,nzndg = ',nxndg,nyndg,nzndg
  ALLOCATE(uincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:uincr")
  uincr = 0.
  ALLOCATE(vincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vincr")
  vincr = 0.
  ALLOCATE(wincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wincr")
  wincr = 0.
  ALLOCATE(pincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pincr")
  pincr = 0.
  ALLOCATE(ptincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptincr")
  ptincr = 0.
  ALLOCATE(qvincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvincr")
  qvincr = 0.
  ALLOCATE(qcincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qcincr")
  qcincr = 0.
  ALLOCATE(qrincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qrincr")
  qrincr = 0.
  ALLOCATE(qiincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qiincr")
  qiincr = 0.
  ALLOCATE(qsincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qsincr")
  qsincr = 0.
  ALLOCATE(qhincr(nxndg,nyndg,nzndg),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qhincr")
  qhincr = 0.

  ALLOCATE(temxy1(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:temxy1")
  temxy1 = 0

  ALLOCATE(tem1soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem1soil")
  tem1soil = 0
  ALLOCATE(tem2soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem2soil")
  tem2soil = 0
  ALLOCATE(tem3soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem3soil")
  tem3soil = 0
  ALLOCATE(tem4soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem4soil")
  tem4soil = 0
  ALLOCATE(tsdiffus(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tsdiffus")
  tsdiffus = 0

  ALLOCATE(tem1 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem1")
  tem1 = 0
  ALLOCATE(tem2 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem2")
  tem2 = 0
  ALLOCATE(tem3 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem3")
  tem3 = 0
  ALLOCATE(tem4 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem4")
  tem4 = 0
  ALLOCATE(tem5 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem5")
  tem5 = 0
  ALLOCATE(tem6 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem6")
  tem6 = 0
  ALLOCATE(tem7 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem7")
  tem7 = 0
  ALLOCATE(tem8 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem8")
  tem8 = 0
  ALLOCATE(tem9 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem9")
  tem9 = 0
  ALLOCATE(tem10(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem10")
  tem10 = 0
  ALLOCATE(tem11(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem11")
  tem11 = 0
  ALLOCATE(tem12(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem12")
  tem12 = 0
  ALLOCATE(tem13(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem13")
  tem13 = 0
  ALLOCATE(tem14(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem14")
  tem14 = 0
  ALLOCATE(tem15(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem15")
  tem15 = 0
  ALLOCATE(tem16(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem16")
  tem16 = 0

  ALLOCATE(tem4d(nx,ny,nz,MAX(10,nscalar)),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem4d")
  tem4d = 0.0

  ALLOCATE(tem1_0(0:nx,0:ny,0:nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem1_0")
  tem1_0 = 0
  ALLOCATE(tem2_0(0:nx,0:ny,0:nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem2_0")
  tem2_0 = 0
  ALLOCATE(tem3_0(0:nx,0:ny,0:nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem3_0")
  tem3_0 = 0

  IF( dfilter_opt == 1 ) THEN
    numdfvar_atm3d = numdfvar3d + nscalarq
    ALLOCATE(df_store_atm3d  (nx,ny,nz,numdfvar_atm3d),             STAT=istatus)
    CALL check_alloc_status(istatus, "arps:df_store_atm3d")
    ALLOCATE(df_store_soil   (nx,ny,nzsoil,0:nstyps,numdfvar_soil), STAT=istatus)
    CALL check_alloc_status(istatus, "arps:df_store_soil")
    ALLOCATE(df_store_soil2d (nx,ny,0:nstyps,numdfvar_soil2d),      STAT=istatus)
    CALL check_alloc_status(istatus, "arps:df_store_soil2d")
    ALLOCATE(df_store_soil2ds(nx,ny,numdfvar_soil2ds),              STAT=istatus)
    CALL check_alloc_status(istatus, "arps:df_store_soil2ds")
    df_store_atm3d   = 0.
    df_store_soil    = 0.
    df_store_soil2d  = 0.
    df_store_soil2ds = 0.
  END IF

  IF (myproc == 0) WRITE(6,*) 'Done allocating arrays'
!
!-----------------------------------------------------------------------
!
!  Set up exbcbuf pointers.
!
!-----------------------------------------------------------------------
!
  IF( lbcopt == 2 .AND. mptr == 1 ) CALL setexbcptr(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Initialize all arrays and variables, and set initial
!  condition for thermal disturbance if desired.
!
!-----------------------------------------------------------------------
!
  CALL initial(mptr,nx,ny,nz,nzsoil,nxndg,nyndg,nzndg,nstyps,exbcbufsz, &
               u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                   &
               udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,         &
               wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,         &
               sdteb,sdtwb,sdtnb,sdtsb,                                 &
               ubar,vbar,wbar,ptbar,pbar,ptbari,pbari,rhostr,rhostri,   &
               qvbar,ppi,csndsq,                                        &
               x,y,z,zp,zpsoil,hterain,mapfct,                          &
               j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,          &
               trigs1,trigs2,ifax1,ifax2,                               &
               wsave1,wsave2,vwork1,vwork2,                             &
               sinlat,kmh,kmv,rprntl,                                   &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,ptsfc,qvsfc,                &
               ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                     &
               cldefi,xland,bmjraincv,                                  &
               raing,rainc,prcrate,exbcbuf,bcrlx,radfrc,radsw,          &
               rnflx,radswnet,radlwin,usflx,vsflx,ptsflx,qvsflx,        &
               uincr,vincr,wincr,pincr,ptincr,qvincr,                   &
               qcincr,qrincr,qiincr,qsincr,qhincr,                      &
               tem1soil,tem2soil,tem3soil,tem4soil,tsdiffus,            &
               temxy1,tem1,tem2,tem3,tem4,tem5,tem6,tem7,               &
               tem8,tem9,tem10,tem11,tem12,tem13,tem14,tem15,tem16,     &
               tem4d,tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  The current model time (initial time):
!
!-----------------------------------------------------------------------
!
  dcurtim = tstart
  curtim = REAL(dcurtim)

  IF( myproc == 0 ) WRITE(6,'(1x,a,f13.3,a)')                           &
       'The initial model time is at ', curtim,' seconds.'
!
!-----------------------------------------------------------------------
!
!  Output the initial fields
!
!-----------------------------------------------------------------------
!
  CALL set_acct(output_acct)
  CALL initout(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,                   &
               u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                   &
               ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,               &
               x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,                 &
               j3inv,                                                   &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,exbcbuf,                             &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11)

  IF (mp_opt > 0) THEN
    CALL mpbarrier ()
  END IF

  CALL set_acct(misc_acct)

  df_dstp=nint(df_tinv/dtbig)
  df_ststp=nint(df_tstart/dtbig)
  addwght=0
  idfstp=1
  df_curtim0=curtim
!
!-----------------------------------------------------------------------
!
!  Time integration loop begins  ----------------------------------->
!
!-----------------------------------------------------------------------
!
  899 CONTINUE
  curtim = df_curtim0
  dcurtim = curtim
  frstep=1
  nstep = nint( curtim/dtbig )

  IF (myproc == 0) WRITE(6,'(1x,a,I8,F12.1)') 'now the nstep,curtim are ',nstep,curtim
  IF (idfstp > df_nstps .and. dfilter_opt == 1) THEN
    CALL initout(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,                 &
               u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                   &
               ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,               &
               x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,                 &
               j3inv,                                                   &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,exbcbuf,                             &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11)
     IF (myproc == 0) WRITE(6,*)'after initout for the changed fields, set dfilter_opt=-1'
     dfilter_opt = 0
   END IF

   900   CONTINUE      ! Entry point of time step integration.
   IF( dfilter_opt == 1 .and. nstep >= df_ststp .and.                   &
       mod(nstep-df_ststp,df_dstp) == 0 .and. idfstp <= df_nstps)  THEN

     WRITE(6,*) 'thinkdfilter df_tstart,df_tinv,df_nstps = ',           &
                df_tstart,df_tinv,df_nstps
     WRITE(6,*) 'think, df_ststp,dfdstp,df_nstps = ',                   &
                df_ststp,df_dstp,df_nstps
     WRITE(6,*) 'thinkdfilter, nstep,idfstp = ',                        &
                nstep,idfstp,'  df_wght = ',df_wght(idfstp)

     CALL dfilter(mptr,nx,ny,nz,nzsoil,nstyps,                          &
             u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                     &
             tsoil,qsoil,wetcanp,snowdpth,                              &
             numdfvar_atm3d,numdfvar_soil,                              &
             numdfvar_soil2d,numdfvar_soil2ds,                          &
             idfstp,df_nstps,df_wght(idfstp),addwght,                   &
             df_store_atm3d,df_store_soil,                              &
             df_store_soil2d,df_store_soil2ds)

     idfstp=idfstp+1
   END IF
   IF(idfstp > df_nstps .and. dfilter_opt == 1) THEN
     WRITE(6,*) 'dfilter finished, goto 899'
     GOTO 899
   END IF

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
  CALL cordintg(mptr,frstep, nx,ny,nz,nzsoil,nxndg,nyndg,nzndg,         &
                rbufsz,nstyps,exbcbufsz,                                &
                u,v,w,wcont,ptprt,pprt,qv,qscalar,                      &
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
                tsoil,qsoil,wetcanp,                                    &
                snowdpth,ptsfc,qvsfc,                                   &
                ptcumsrc,qcumsrc,raing,rainc,prcrate,                   &
                w0avg,nca,kfraincv,                                     &
                cldefi,xland,bmjraincv,                                 &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                rad2d,radbuf,exbcbuf,bcrlx,                             &
                usflx,vsflx,ptsflx,qvsflx,                              &
                uincr,vincr,wincr,pincr,ptincr,qvincr,                  &
                qcincr,qrincr,qiincr,qsincr,qhincr,                     &
                tem1soil,tem2soil,tem3soil,tem4soil,tsdiffus,           &
                temxy1,tem1,tem2,tem3,tem4,tem5,tem6,tem7,              &
                tem8,tem9,tem10,tem11,tem12,tem13,                      &
                tem14,tem15,tem16,tem4d,                                &
                tem1_0,tem2_0,tem3_0)

  CALL set_acct(misc_acct)
!
!-----------------------------------------------------------------------
!
!  Update physical time and generate output files
!
!-----------------------------------------------------------------------
!
  dcurtim = dcurtim + dtbig
  curtim = REAL(dcurtim)
!
!-----------------------------------------------------------------------
!
!  Additional module for Asynchronous EnKF(AEnKF)
!  Output the H(x) according to the obtained observations
!-----------------------------------------------------------------------
!


!=======================================================================

!
!-----------------------------------------------------------------------
!
!  Print timestep marker
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)  &
    WRITE(6,'(a,i8,a,f10.2,a,f8.3,a)')                                  &
        '  End of integration time step ', nstep,                       &
        ', Model time=',curtim,' s (', curtim/3600.0, ' h)'

  CALL flush (6)   ! GMB wdt

  CALL set_acct(output_acct)

  IF(dfilter_opt /= 1 .OR. idfstp > df_nstps) THEN
    CALL output(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,                  &
             u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                     &
             udteb,udtwb,vdtnb,vdtsb,pdteb,pdtwb,pdtnb,pdtsb,           &
             ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,                 &
             x,y,z,zp,zpsoil,hterain, mapfct,                           &
             j1,j2,j3,j3soil,j3inv,j3soilinv,                           &
             soiltyp,stypfrct,vegtyp,lai,roufns,veg,                    &
             tsoil,qsoil,wetcanp,snowdpth,qvsfc,                        &
             ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                       &
             cldefi,xland,bmjraincv,                                    &
             raing,rainc,prcrate,exbcbuf,                               &
             radfrc,radsw,rnflx, radswnet,radlwin,                      &
             usflx,vsflx,ptsflx,qvsflx,                                 &
             tem1, tem2,tem3,tem4,tem5, tem6,tem7,                      &
             tem8,tem9,tem10,tem11)
  END IF ! no output during the digital filtering

  CALL set_acct(misc_acct)
!
!-----------------------------------------------------------------------
!
!  Check the stability of the integration.  If the model solution
!  appears to be exceeding the linear stability condition, then
!  perform a data dump for post-mortem.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem11(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  CALL chkstab(mptr,nx,ny,nz,nzsoil,nstyps,                             &
               u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                   &
               ubar,vbar,ptbar,pbar,tem11,qvbar,kmh,kmv,                &
               x,y,z,zp,zpsoil,hterain,mapfct, j1,j2,j3,j3soil,         &
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
               qv,qvbar,qscalar,rhostr,x,y,zp,j3,j3inv,                 &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11)
!
!-----------------------------------------------------------------------
!
!  Time integration loop ends, go back to the beginning of the
!  time step loop if stop time is not reached.
!
!-----------------------------------------------------------------------
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
    CALL mclosedataset(dsindex, ierr)

  ELSE IF (hdmpfmt == 9) THEN

    IF ((mp_opt > 0 .AND. joindmp(FINDX_H) <= 0) .OR. myproc == 0) CLOSE (nchdmp)

  END IF

  CALL acct_finish
  CALL acct_report_arps

  IF (myproc == 0) THEN

    WRITE(6,'(//1x,a,/1x,a,f13.3,a,/1x,a)')                             &
      'ARPS stopped normally in the main program. ',                    &
      'The ending time was ', curtim,' seconds.',                       &
      'Thanks for using ARPS.'

    WRITE (6,*) "Maxumum memory allocation (in words):",                &
                max_memory_use
  END IF

  IF (mp_opt > 0) CALL mpexit(0)

  STOP
END PROGRAM ARPS
