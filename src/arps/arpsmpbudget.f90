PROGRAM arpsmpbudget
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSMPBUDGET               ######
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
!  This program calculates microphysics budget statistics from a series
!  of ARPS history dumps and optional additional microphysics rate
!  fields (i.e. evaporation, melting rates) in individual files.  To
!  produce these fields, an additional switch in the arps input file must
!  be turned on (to be implemented).  The program outputs vertical profiles
!  of horizontally and temporally summed total cooling/heating, horizontally
!  and temporally averaged mixing ratios, number concentrations, and reflectivity
!  (the latter two for the multimoment MY scheme), among other statistics.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!    1/11/2008
!
!  MODIFICATION HISTORY:
!
!    1/11/2008 D. Dawson
!    First split off into new program from old hacked code.
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
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
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!  CALCULATED DATA ARRAYS:
!
!    u        x component of total velocity (m/s)
!    v        y component of total velocity (m/s)
!    w        z component of total velocity (m/s)
!    qv       total water vapor mixing ratio (kg/kg)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Work array for Savi3D dump
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nzsoil   ! Grid dimensions.
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,nf,lenfil
  PARAMETER (nhisfile_max=2000)
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
  INCLUDE 'indtflg.inc'
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Parameters to be passed in BN2DUMP and SVIDUMP.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ist,ind,isk,jst,jnd,jsk,kst,knd,ksk
  COMMON /dfndomn/ ist,ind,isk,jst,jnd,jsk,kst,knd,ksk
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
  REAL, ALLOCATABLE :: zs_agl(:,:,:)  ! The height of scalar points above ground.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height of the terrain.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)
  REAL, ALLOCATABLE :: qsclhavg(:,:)
  REAL, ALLOCATABLE :: qsclhtavg(:,:)

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
                                      ! humidity (kg/kg)

  REAL, ALLOCATABLE :: u      (:,:,:) ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v      (:,:,:) ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w      (:,:,:) ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: ws     (:,:,:) ! Total w-velocity averaged to scalar points (m/s)
  REAL, ALLOCATABLE :: qv     (:,:,:) ! Water vapor specific humidity (kg/kg)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct(:,:,:)    ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: lai    (:,:)   ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)   ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)   ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil (:,:,:,:) ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil (:,:,:,:) ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (m)

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
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: mpcoolterms(:,:,:,:)  ! Hydrometeor process rates contributing to cooling (kg/kg/s)
  REAL, ALLOCATABLE :: mpheatterms(:,:,:,:)  ! Hydrometeor process rates contributing to heating (kg/kg/s)
  REAL, ALLOCATABLE :: mptemprate(:,:,:)  ! Microphysical temperature change rate (K/s)
  REAL*8, ALLOCATABLE :: N0(:,:,:,:)   ! Intercept parameter
  REAL*8, ALLOCATABLE :: N0havg(:,:)
  REAL*8, ALLOCATABLE :: N0htavg(:,:)

  REAL, ALLOCATABLE :: Ntx(:,:,:,:)   ! Number concentration (only used for SM schemes)
                                      ! MM schemes have this as part of qscalar array

  REAL, ALLOCATABLE :: Ntxhavg(:,:)   ! Horizontally-averaged Nt (for SM schemes)
  REAL, ALLOCATABLE :: Ntxhtavg(:,:)  ! Horizontally and time-averaged Nt (for SM schemes)

  REAL*8, ALLOCATABLE :: N0eff(:,:,:,:) ! "Effective" intercept parameter
  REAL*8, ALLOCATABLE :: lamda(:,:,:,:) ! Slope parameter
  REAL*8, ALLOCATABLE :: alpha(:,:,:,:) ! Shape parameter
  REAL*8, ALLOCATABLE :: lamdahavg(:,:)
  REAL*8, ALLOCATABLE :: alphahavg(:,:)
  REAL*8, ALLOCATABLE :: lamdahtavg(:,:)
  REAL*8, ALLOCATABLE :: alphahtavg(:,:)

  REAL, ALLOCATABLE :: Dm(:,:,:,:)    ! Mean-mass diameter
  REAL, ALLOCATABLE :: Dmhavg(:,:)
  REAL, ALLOCATABLE :: Dmhtavg(:,:)

  REAL, ALLOCATABLE :: gridvol(:,:,:)        ! Grid cell volume around scalar point (m3)
  REAL, ALLOCATABLE :: gridmass(:,:,:)       ! Mass of dry air around scalar point (kg)
  REAL, ALLOCATABLE :: temscalar(:,:,:,:)    ! Temporary array to hold past qscalar data
  REAL, ALLOCATABLE :: rh(:,:,:)             ! Relative humidity
  REAL, ALLOCATABLE :: rho(:,:,:)            ! Air density

!
!-----------------------------------------------------------------------
!
!  Temporary working arrays:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1avg(:,:,:)
  REAL, ALLOCATABLE :: tem2avg(:,:,:)
  REAL, ALLOCATABLE :: tem3avg(:,:,:)

!
!-----------------------------------------------------------------------
!
!  Temporary working arrays:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: houtfmt, nchin, nchout

  CHARACTER (LEN=132) :: houtfn
  CHARACTER (LEN=132) :: basoutfn
  INTEGER :: lbasoutf, loutf

  INTEGER :: grdbas
  INTEGER :: i,j,k,ireturn

  REAL :: alatnot(2)
  REAL :: ctrx, ctry, swx, swy

  REAL :: time
  INTEGER :: gboutcnt, vroutcnt
  DATA    gboutcnt, vroutcnt /0,0/

  INTEGER :: nfile
  INTEGER,PARAMETER :: nzmax = 256
  CHARACTER(LEN=80) :: inrunname, outrunname

  INTEGER :: grdbasopt           ! 0 - no  grdbas file
                                 ! 1 - one grdbas file, the first file only

  INTEGER :: nq,ctravgi

  INTEGER :: bibgn,biend,bjbgn,bjend,bkbgn,bkend
  INTEGER :: nbudgetsteps, nb
  REAL :: budgtbgn,budgtend
!  CHARACTER(LEN=80) :: boutfile

  REAL, DIMENSION(nzmax) :: evapqctot, evapqrtot, sublqitot, sublqstot, sublqgtot,        &
          sublqhtot, meltqitot, meltqstot, meltqgtot, meltqhtot,        &
          evaptot, subltot, melttot, evapqccooltot, evapqrcooltot,      &
          sublqicooltot, sublqscooltot, sublqgcooltot, sublqhcooltot,   &
          meltqicooltot, meltqscooltot, meltqgcooltot, meltqhcooltot,   &
          evapcooltot, sublcooltot, meltcooltot, mptchgtot, mpcooltot,  &
          mpheattot, mpQchgtot,mpQchgcalctot

  REAL, DIMENSION(nzmax) :: evapqccumtot, evapqrcumtot, sublqicumtot, sublqscumtot, sublqgcumtot,        &
          sublqhcumtot, meltqicumtot, meltqscumtot, meltqgcumtot, meltqhcumtot

  REAL, DIMENSION(nzmax) :: mptchgcumtot,mpQchgcumtot,mpQchgcalccumtot

  REAL :: evapqccumcoltot, evapqrcumcoltot, sublqicumcoltot, sublqscumcoltot, sublqgcumcoltot,        &
          sublqhcumcoltot, meltqicumcoltot, meltqscumcoltot, meltqgcumcoltot, meltqhcumcoltot

  REAL :: mptchgcumcoltot,mpQchgcumcoltot

  REAL, DIMENSION(nzmax) :: condqctot, condqrtot, nuclqitot, depoqitot, depoqstot,        &
          depoqgtot, depoqhtot, frzqcitot, colqcitot, colqcstot,        &
          colqcgtot, colqchtot, frzqrhtot, colqritot, colqrstot,        &
          colqrgtot, colqrhtot, condtot, depotot, frztot,               &
          condqcheattot, condqrheattot, nuclqiheattot, depoqiheattot,   &
          depoqsheattot, depoqgheattot, depoqhheattot, frzqciheattot,   &
          colqciheattot, colqcsheattot, colqcgheattot, colqchheattot,   &
          frzqrhheattot, colqriheattot, colqrsheattot, colqrgheattot,   &
          colqrhheattot, condheattot, depoheattot, frzheattot

  REAL, DIMENSION(nzmax) :: condqccumtot, condqrcumtot, nuclqicumtot, depoqicumtot, depoqscumtot,        &
          depoqgcumtot, depoqhcumtot, frzqcicumtot, colqcicumtot, colqcscumtot,        &
          colqcgcumtot, colqchcumtot, frzqrhcumtot, colqricumtot, colqrscumtot,        &
          colqrgcumtot, colqrhcumtot

  REAL, DIMENSION(nzmax) :: gridvoltot, gridmasstot

  REAL, DIMENSION(nzmax) :: chngqwtot, chngqwcumtot
  REAL :: chngqwcumcoltot

  REAL :: wthresh, zthresh

  REAL, DIMENSION(nzmax) :: vmasstot,cmasstot,rmasstot,imasstot,smasstot,gmasstot,hmasstot
  REAL, DIMENSION(nzmax) :: qvhavg
  REAL, DIMENSION(nzmax) :: rhhavg

  REAL, DIMENSION(nzmax) :: qvhtavg
  REAL, DIMENSION(nzmax) :: rhhtavg

  REAL, DIMENSION(nzmax,7) :: rhotot

  REAL :: condqccumcoltot, condqrcumcoltot, nuclqicumcoltot, depoqicumcoltot, depoqscumcoltot,        &
          depoqgcumcoltot, depoqhcumcoltot, frzqcicumcoltot, colqcicumcoltot, colqcscumcoltot,        &
          colqcgcumcoltot, colqchcumcoltot, frzqrhcumcoltot, colqricumcoltot, colqrscumcoltot,        &
          colqrgcumcoltot, colqrhcumcoltot, gridmasscumcoltot

  REAL, DIMENSION(nzmax) :: loadqctot,loadqrtot,loadqitot,loadqstot,loadqgtot,loadqhtot
  REAL, DIMENSION(nzmax) :: loadqccumtot,loadqrcumtot,loadqicumtot,loadqscumtot,loadqgcumtot,loadqhcumtot

  REAL, PARAMETER :: LV = 2.501e6     ! Latent heat of vaporization/condensation (J/kg)
  REAL, PARAMETER :: LF = 3.34e5       ! Latent heat of melting/freezing (J/kg)
  REAL, PARAMETER :: LS = LV+LF       ! Latent heat of sublimation/deposition (J/kg)
!  REAL, PARAMETER :: CP = 1005.46     ! Specific heat at constant pressure for dry air (J/kg/K)
!  REAL, PARAMETER :: RD = 287.05       ! Gas constant for dry air (J/kg/K)
!  REAL, PARAMETER :: p0 = 1.0e5       ! Reference surface pressure (Pa)
!  REAL, PARAMETER :: p0inv = 1.0/p0
!  REAL, PARAMETER :: rddcp = RD/CP

  REAL :: kgtokt, JtoGJ

  INTEGER :: normalize_opt, wopt,dobudget,doavg,dotimeseries
  CHARACTER (LEN=256) :: hdmpftrailer

  NAMELIST /budget/ dobudget,doavg,dotimeseries,outrunname,dirname,hdmpftrailer,budgtbgn,budgtend,normalize_opt,   &
                    wopt,wthresh,zthresh,tintv_dmpin,tintv,bibgn,biend,bjbgn,bjend,bkbgn,bkend

  CHARACTER (LEN=132) :: budgetfn
  INTEGER :: lbudgetfn, nchbudget, istat
  INTEGER :: tintv_dmpin,tintv

  INTEGER :: ierr

  CHARACTER (LEN=40 ) :: varunits
  CHARACTER (LEN=40 ) :: varname
  CHARACTER (LEN=6 ) :: varid

  CHARACTER (LEN=40) :: mpcoolunits(10)
  CHARACTER (LEN=40) :: mpcoolname(10)
  CHARACTER (LEN=6) :: mpcoolid(10)

  CHARACTER (LEN=40) :: mpheatunits(17)
  CHARACTER (LEN=40) :: mpheatname(17)
  CHARACTER (LEN=6) :: mpheatid(17)

  INTEGER :: istatus, nqscalar

  INTEGER :: nunit_ts,nunit_ts2, lts_fn_out
  CHARACTER (LEN=256) :: ts_fn_out

  INTEGER :: nunit_heatprof,nunit_coolprof,nunit_qnzavgprof,nunit_mpavgprof,lprof_fn_out
  CHARACTER (LEN=256) :: prof_fn_out

  INTEGER :: gridcount
  INTEGER :: qflag(6),qcount(6)

  INTEGER, ALLOCATABLE :: numtimes1(:,:),numtimes2(:),qtimes(:,:)

  REAL :: pi
  REAL, PARAMETER :: CP = 1005.46     ! Specific heat at constant pressure for dry air (J/kg/K)
  REAL, PARAMETER :: RD = 287.05       ! Gas constant for dry air (J/kg/K)
  REAL, PARAMETER :: p0 = 1.0e5       ! Reference surface pressure (Pa)
  REAL, PARAMETER :: p0inv = 1.0/p0
  REAL, PARAMETER :: rddcp = RD/CP

  REAL :: cx(6)     ! (PI/6)*rho_qx
  REAL :: Ntcfix,N0rfix,N0sfix,N0gfix,N0hfix
  REAL :: alpharfix,alphaifix,alphasfix,alphagfix,alphahfix
  REAL :: rhor,rhoi,rhos,rhog,rhoh
  NAMELIST /microph_param/ mphyopt,dirname,hdmpftrailer,Ntcfix,N0rfix,N0sfix,N0gfix,N0hfix,alpharfix,alphaifix,alphasfix,alphagfix,alphahfix,rhor,rhoi,rhos,rhog,rhoh

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(/9(/5x,a)/)')                                               &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                Welcome to ARPSMPBUDGET                  ###', &
     '###      This program calculates various microphysical      ###', &
     '###      budget statistics from ARPS history data.          ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

  CALL mpinit_var
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

!-----------------------------------------------------------------------
!
! Read in namelists
!
!-----------------------------------------------------------------------

  READ(5,budget,END=100)

  GOTO 10

  WRITE(6,'(1x,/a/)') 'Namelist budget read in successfully.'

  100   WRITE(6,'(a)')                                                  &
          'Error reading NAMELIST file. Program ARPSCVT stopped.'
  STOP

  10    CONTINUE

  READ(5,microph_param,ERR=100)

  WRITE(6,'(a)')'Namelist microph_param was successfully read.'

  write(6,microph_param)

  IF(mphyopt >= 5 .and. mphyopt <= 7 .or. mphyopt >= 12) THEN
     write(6,'(a)')'Microphysics option not supported!'
     STOP
  END IF

  pi = 4.0 * atan(1.0)

  cx(1) = (pi/6.)*rhor
  cx(2) = (pi/6.)*rhor
  cx(3) = 440.0      ! Constant value used in MY scheme for ice (bullet rosettes)
  cx(4) = (pi/6.)*rhos
  IF(mphyopt == 2) THEN
    cx(5)=(pi/6.)*rhoh
  ELSE
    cx(5) = (pi/6.)*rhog
  END IF
  cx(6) = (pi/6.)*rhoh

  IF(hinfmt == 0) GO TO 888  ! Do not convert history files
                             ! go ahead for terrain data etc.

!-----------------------------------------------------------------------
!
! Get dimension parameters and allocate arrays
!
!-----------------------------------------------------------------------

  WRITE(6,'(2a)') '  Reading dimensions from input file - ',TRIM(grdbasfn)

  lenfil = len_trim(hisfile(1))

  CALL get_dims_from_data(hinfmt,hisfile(1)(1:lenfil),                    &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  ! Determine number of hydrometeor species (nqscalar)
  IF(nscalar >= 6) THEN
    nqscalar = 6
  ELSE IF(nscalar == 5) THEN
    nqscalar = 5
  ELSE
    nqscalar = 2
  END IF

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(2x,4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz,', nzsoil=',nzsoil

  ALLOCATE(x      (nx))
  ALLOCATE(y      (ny))
  ALLOCATE(z      (nz))
  ALLOCATE(zp     (nx,ny,nz))
  ALLOCATE(zs_agl (nx,ny,nz))
  ALLOCATE(zpsoil (nx,ny,nzsoil))

  ALLOCATE(uprt   (nx,ny,nz))
  ALLOCATE(vprt   (nx,ny,nz))
  ALLOCATE(wprt   (nx,ny,nz))
  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))

  ALLOCATE(qscalar(nx,ny,nz,nscalar))
  ALLOCATE(qsclhavg(nz,nscalar))
  ALLOCATE(qsclhtavg(nz,nscalar))

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
  ALLOCATE(u      (nx,ny,nz))
  ALLOCATE(v      (nx,ny,nz))
  ALLOCATE(w      (nx,ny,nz))
  ALLOCATE(ws     (nx,ny,nz))
  ALLOCATE(qv     (nx,ny,nz))

  ALLOCATE(soiltyp (nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp (nx,ny))
  ALLOCATE(lai    (nx,ny))
  ALLOCATE(roufns (nx,ny))
  ALLOCATE(veg    (nx,ny))

  ALLOCATE(tsoil (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(raing(nx,ny))
  ALLOCATE(rainc(nx,ny))
  ALLOCATE(prcrate(nx,ny,4))

  ALLOCATE(radfrc(nx,ny,nz))
  ALLOCATE(radsw (nx,ny))
  ALLOCATE(rnflx (nx,ny))
  ALLOCATE(radswnet (nx,ny))
  ALLOCATE(radlwin (nx,ny))


  ALLOCATE(usflx (nx,ny))
  ALLOCATE(vsflx (nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))
  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))

  ALLOCATE(mpcoolterms(nx,ny,nz,10))
  ALLOCATE(mpheatterms(nx,ny,nz,17))
  ALLOCATE(mptemprate(nx,ny,nz))
  ALLOCATE(Ntx(nx,ny,nz,nqscalar))
  ALLOCATE(Ntxhavg(nz,nqscalar))
  ALLOCATE(Ntxhtavg(nz,nqscalar))
  ALLOCATE(N0(nx,ny,nz,nqscalar))
  ALLOCATE(N0havg(nz,nqscalar))
  ALLOCATE(N0htavg(nz,nqscalar))
  ALLOCATE(N0eff(nx,ny,nz,nqscalar))
  ALLOCATE(lamda(nx,ny,nz,nqscalar))
  ALLOCATE(lamdahavg(nz,nqscalar))
  ALLOCATE(lamdahtavg(nz,nqscalar))
  ALLOCATE(alpha(nx,ny,nz,nqscalar))
  ALLOCATE(alphahavg(nz,nqscalar))
  ALLOCATE(alphahtavg(nz,nqscalar))
  ALLOCATE(Dm(nx,ny,nz,nqscalar))
  ALLOCATE(Dmhavg(nz,nqscalar))
  ALLOCATE(Dmhtavg(nz,nqscalar))
  ALLOCATE(gridvol(nx,ny,nz))
  ALLOCATE(gridmass(nx,ny,nz))
  ALLOCATE(temscalar(nx,ny,nz,nqscalar))
  ALLOCATE(rh(nx,ny,nz))
  ALLOCATE(rho(nx,ny,nz))
  ALLOCATE(numtimes1(nz,7))
  ALLOCATE(numtimes2(nz))
  ALLOCATE(qtimes(6,nz))

  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0
  zs_agl =0.0
  zpsoil =0.0

  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
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
  u      =0.0
  v      =0.0
  w      =0.0
  ws     =0.0
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

  raing=0.0
  rainc=0.0
  prcrate=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet = 0.0
  radlwin = 0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0
  tem1=0.0
  tem2=0.0
  tem3=0.0

  mpcoolterms=0.0
  mpheatterms=0.0
  mptemprate=0.0
  Ntx=0.0
  N0=0.0
  N0eff=0.0
  lamda=0.0
  alpha=0.0
  Dm=0.0

  gridvol = 0.0
  gridmass = 0.0
  temscalar = 0.0
  rh = 0.0
  rho = 0.0
  !wdt Copyright (c) 2001 Weather Decision Technologies, Inc.: ngbrz,zbrdmp
  ngbrz = 5
  zbrdmp = 10000.0

!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!
  ldirnam=LEN(dirname)
  CALL strlnth( dirname , ldirnam)

  lengbf=len_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

  DO nfile = 1,nhisfile

    lenfil=len_trim(hisfile(nfile))
    WRITE(6,'(/a,a,a)')                                                 &
        ' Data set ', trim(hisfile(nfile)),' to be converted.'

!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!

    102   CONTINUE

    CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
                 hinfmt, nchin,grdbasfn(1:lengbf),lengbf,               &
                 hisfile(nfile)(1:lenfil),lenfil,time,                  &
                 x,y,z,zp,zpsoil, uprt ,vprt ,wprt ,ptprt, pprt ,       &
                 qvprt, qscalar,tke,kmh,kmv,                            &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil, qsoil, wetcanp,snowdpth,                        &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1, tem2, tem3)

    ! Calculate grid cell volume and grid cell mass (of dry air)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ! Calculate the mass of dry air at each grid point
          ! Temperature is in tem1, density in rho

          tem1(i,j,k)=(ptbar(i,j,k)+ptprt(i,j,k))*                        &
                (((pbar(i,j,k)+pprt(i,j,k))*p0inv)**rddcp)
          rho(i,j,k)=(pbar(i,j,k)+pprt(i,j,k))/(RD*tem1(i,j,k))

          gridvol(i,j,k) = (zp(i,j,k+1)-zp(i,j,k))*(x(i+1)-x(i))*(y(j+1)-y(j))
          gridmass(i,j,k) = rho(i,j,k)*gridvol(i,j,k)

          ! Calculate height of scalar points above ground level and
          ! vertical velocity averaged to scalar points
          zs_agl(i,j,k) = (zp(i,j,k)+zp(i,j,k+1))/2.-zp(i,j,2)
          ws(i,j,k) = (wprt(i,j,k)+wprt(i,j,k+1))/2.
        END DO
      END DO
    END DO

    ! Calculate various size-distribution related parameters

    IF(doavg == 1) THEN
    IF(mphyopt == 2) THEN
      DO nq=1,5
        ! Calculate N0 (in this case constant everywhere)

        IF(nq == 1) THEN
          N0(:,:,:,nq) = 0.d0   ! Not applicable to monodisperse cloud
        ELSE IF(nq == 2) THEN
          N0(:,:,:,nq) = dble(N0rfix)
        ELSE IF(nq == 3) THEN
          N0(:,:,:,nq) = 0.d0   ! Not applicable to monodisperse ice
        ELSE IF(nq == 4) THEN
          N0(:,:,:,nq) = dble(N0sfix)
        ELSE IF(nq == 5) THEN
          N0(:,:,:,nq) = dble(N0hfix)
        END IF

        ! Calculate alpha (in this case zero everywhere)

        alpha(:,:,:,nq) = 0.d0

        ! Calculate Nt from values of q and N0

        IF(nq == 1) THEN
          Ntx(:,:,:,nq) = Ntcfix  ! Fixed number concentration for cloud droplets
        ELSE IF(nq == 3) THEN
          Ntx(:,:,:,nq) = 0.0     ! Fletcher equation for ice number concentration: TODO: need to determine this for LIN scheme
        ELSE
          CALL cal_Nt(nx,ny,nz,rho,qscalar(1,1,1,nq),N0(1,1,1,nq),cx(nq),alpha(1,1,1,nq),Ntx(1,1,1,nq))
        END IF

        ! Calculate mean-mass diameter

        CALL cal_Dm(nx,ny,nz,rho,qscalar(1,1,1,nq),Ntx(1,1,1,nq),cx(nq),Dm(1,1,1,nq))

        ! Calculate lamda

        CALL cal_lamda(nx,ny,nz,rho,qscalar(1,1,1,nq),Ntx(1,1,1,nq),cx(nq),alpha(1,1,1,nq),lamda(1,1,1,nq))

      END DO
    ELSE IF(mphyopt >= 5 .and. mphyopt <= 7) THEN
!       DO nq=1,6
!         IF(nq /= 6) THEN
!         END IF
!       END DO
    ELSE IF(mphyopt == 8) THEN
      DO nq=1,6

        ! Calculate alpha

        IF(nq == 1) THEN
          alpha(:,:,:,nq) = 1.d0  ! Fixed alpha = 1 for cloud droplets
        ELSE IF(nq == 2) THEN
          alpha(:,:,:,nq) = dble(alpharfix)
        ELSE IF(nq == 3) THEN
          alpha(:,:,:,nq) = dble(alphaifix)
        ELSE IF(nq == 4) THEN
          alpha(:,:,:,nq) = dble(alphasfix)
        ELSE IF(nq == 5) THEN
          alpha(:,:,:,nq) = dble(alphagfix)
        ELSE IF(nq == 6) THEN
          alpha(:,:,:,nq) = dble(alphahfix)
        END IF

        ! Calculate N0

        IF(nq == 1) THEN
          ! Nt constant for cloud
          Ntx(:,:,:,nq) = 1.0e8
          CALL cal_N0(nx,ny,nz,rho,qscalar(1,1,1,nq),Ntx(1,1,1,nq),cx(nq),alpha(1,1,1,nq),N0(1,1,1,nq),N0eff(1,1,1,nq))
        ELSE IF(nq == 2) THEN
          N0(:,:,:,nq) = dble(N0rfix)
        ELSE IF(nq == 3) THEN
          N0(:,:,:,nq) = 0.d0 ! Not applicable (Ntx calculated from Cooper eqn)
        ELSE IF(nq == 4) THEN
          N0(:,:,:,nq) = dble(N0sfix)
        ELSE IF(nq == 5) THEN
          N0(:,:,:,nq) = dble(N0gfix)
        ELSE IF(nq == 6) THEN
          N0(:,:,:,nq) = dble(N0hfix)
        END IF

        ! Calculate Nt

        IF(nq == 3) THEN
          DO k=1,nz-1
            DO i=1,nx-1
              DO j=1,ny-1
                Ntx(i,j,k,nq) = 5.*exp(0.304*(273.15-max(233.,tem1(i,j,k))))  ! Cooper eqn for Nt for ice
              END DO
            END DO
          END DO
        ELSE IF (nq /= 1) THEN
          CALL cal_Nt(nx,ny,nz,rho,qscalar(1,1,1,nq),N0(1,1,1,nq),cx(nq),alpha(1,1,1,nq),Ntx(1,1,1,nq))
        END IF

        ! Calculate mean-mass diameter

        CALL cal_Dm(nx,ny,nz,rho,qscalar(1,1,1,nq),Ntx(1,1,1,nq),cx(nq),Dm(1,1,1,nq))

        ! Calculate lamda

        CALL cal_lamda(nx,ny,nz,rho,qscalar(1,1,1,nq),Ntx(1,1,1,nq),cx(nq),alpha(1,1,1,nq),lamda(1,1,1,nq))
      END DO
    ELSE IF(mphyopt == 9 .or. mphyopt == 10) THEN

      DO nq=1,6
        ! Calculate Dm

        CALL cal_Dm(nx,ny,nz,rho,qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),cx(nq),Dm(1,1,1,nq))

        ! Calculate alpha

        IF(mphyopt == 9) THEN  ! Fixed alpha
          IF(nq == 1) THEN
            alpha(:,:,:,nq) = 1.d0  ! Fixed alpha = 1 for cloud droplets
          ELSE IF(nq == 2) THEN
            alpha(:,:,:,nq) = dble(alpharfix)
          ELSE IF(nq == 3) THEN
            alpha(:,:,:,nq) = dble(alphaifix)
          ELSE IF(nq == 4) THEN
            alpha(:,:,:,nq) = dble(alphasfix)
          ELSE IF(nq == 5) THEN
            alpha(:,:,:,nq) = dble(alphagfix)
          ELSE IF(nq == 6) THEN
            alpha(:,:,:,nq) = dble(alphahfix)
          END IF
        ELSE IF(mphyopt == 10) THEN ! Diagnostic alpha
          IF(nq == 1) THEN
            alpha(:,:,:,nq) = 1.d0  ! Fixed alpha = 1 for cloud droplets
          ELSE
            CALL diag_alpha(nx,ny,nz,rho,qnames(nq),Dm(1,1,1,nq),alpha(1,1,1,nq))
          END IF
        END IF

        ! Calculate N0 (and "effective" N0)

        CALL cal_N0(nx,ny,nz,rho,qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),cx(nq),alpha(1,1,1,nq),N0(1,1,1,nq),N0eff(1,1,1,nq))

        ! Calculate lamda

        CALL cal_lamda(nx,ny,nz,rho,qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),cx(nq),alpha(1,1,1,nq),lamda(1,1,1,nq))

      END DO
    ELSE IF(mphyopt == 11) THEN

      DO nq=1,6
        ! Calculate Dm
        CALL cal_Dm(nx,ny,nz,rho,qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),cx(nq),Dm(1,1,1,nq))

        ! Calculate alpha
        IF(nq >= 2) THEN
          CALL solve_alpha(nx,ny,nz,rho,cx(nq),qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),qscalar(1,1,1,nq+11),alpha(1,1,1,nq))
        ELSE
          alpha(:,:,:,nq) = 1.d0  ! Fixed alpha for cloud droplets
        END IF

        ! Calculate N0 (and "effective" N0)

        CALL cal_N0(nx,ny,nz,rho,qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),cx(nq),alpha(1,1,1,nq),N0(1,1,1,nq),N0eff(1,1,1,nq))

        ! Calculate lamda

        CALL cal_lamda(nx,ny,nz,rho,qscalar(1,1,1,nq),qscalar(1,1,1,nq+6),cx(nq),alpha(1,1,1,nq),lamda(1,1,1,nq))

      END DO
    END IF

    ! Calculate relative humidity

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ! Store total pressure in tem3
          tem3(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
        END DO
      END DO
    END DO

    tem2=0.0
    ! qvs is in tem2
    CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, tem3,tem1,tem2)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ! Calculate relative humidity
          rh(i,j,k) = MIN( MAX( (qvbar(i,j,k)+qvprt(i,j,k))/tem2(i,j,k), 0.0), 1.0)
          !print*,rh(i,j,k)
        END DO
      END DO
    END DO
    END IF ! doavg == 1

    inrunname = runname
    curtim = time
    IF (nfile == nhisfile) tstop = time

    IF (nscalarin == 0) THEN ! Old version
      nscalar = 5
      P_QC = 1
      P_QR = 2
      P_QI = 3
      P_QS = 4
      P_QH = 5

!   ELSE IF (p_qiin <= 0) THEN    ! we assume 3 possiblility
!     nscalar = 2
!     P_QI = 0               ! 1. only contains qc,qr
!     P_QS = 0               ! 2. only contains qc,qr, qi, qs, qh
!     P_QH = 0               ! 3. contains more qc,qr, qi, qs, qh, ....
!    ELSE IF (p_qxin <= 0) THEN
!      nscalar = 5

    END IF

    IF( hinfmt == 9 .AND. ireturn == 2 ) THEN
      WRITE(6,'(/1x,a/)') 'The end of GrADS file was reached.'
      CLOSE ( nchin )
      CALL retunit( nchin )
      CYCLE
    END IF

    IF( ireturn /= 0 ) GO TO 997             ! Read was unsuccessful

    dx = x(2) - x(1)
    dy = y(2) - y(1)

    IF ( mapproj /= 0 ) THEN
      alatnot(1) = trulat1
      alatnot(2) = trulat2

      CALL setmapr(mapproj, 1.0, alatnot, trulon)
      CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )

      swx = ctrx - (FLOAT(nx-3)/2.) * dx - x(2)
      swy = ctry - (FLOAT(ny-3)/2.) * dy - y(2)

      CALL setorig( 1, swx, swy)

      CALL setcornerll( nx,ny, x,y )       ! set corner lat/lon
    ELSE
      swlats = ctrlat
      swlons = ctrlon
      nelats = ctrlat
      nelons = ctrlon
      swlatu = ctrlat
      swlonu = ctrlon
      nelatu = ctrlat
      nelonu = ctrlon
      swlatv = ctrlat
      swlonv = ctrlon
      nelatv = ctrlat
      nelonv = ctrlon
    END IF

    ternopt = 0
    DO j=1,ny-1
      DO i=1,nx-1
        IF ( zp(i,j,2) /= zp(1,1,2) ) THEN
          ternopt = 2
          GO TO 85
        END IF
      END DO
    END DO

    85    CONTINUE

    IF( LEN_TRIM(outrunname) > 0 ) THEN
      runname = TRIM(outrunname)
    END IF

    IF (curtim == budgtbgn .and. (dobudget == 1 .or. doavg == 1)) THEN

      ! Initialize all variables to contain
      ! cumulative totals to zero

      IF(doavg == 1) THEN
        qvhtavg=0.0
        qsclhtavg=0.0
        Ntxhtavg=0.0
        Dmhtavg=0.0
        rhhtavg=0.0
        alphahtavg=0.0
        N0htavg=0.0
        lamdahtavg=0.0
        numtimes1=0
        numtimes2=0
      END IF

      IF(dobudget == 1) THEN
        evapqccumtot=0.0
        evapqrcumtot=0.0
        sublqicumtot=0.0
        sublqscumtot=0.0
        sublqgcumtot=0.0
        sublqhcumtot=0.0
        meltqicumtot=0.0
        meltqscumtot=0.0
        meltqgcumtot=0.0
        meltqhcumtot=0.0

        evaptot=0.0
        subltot=0.0
        melttot=0.0
        evapqccooltot=0.0
        evapqrcooltot=0.0
        sublqicooltot=0.0
        sublqscooltot=0.0
        sublqgcooltot=0.0
        sublqhcooltot=0.0
        meltqicooltot=0.0
        meltqscooltot=0.0
        meltqgcooltot=0.0
        meltqhcooltot=0.0
        evapcooltot=0.0
        sublcooltot=0.0
        meltcooltot=0.0
        mpcooltot=0.0
        mpheattot=0.0
        mptchgcumtot=0.0
        mpQchgcumtot=0.0

        condqccumtot=0.0
        condqrcumtot=0.0
        nuclqicumtot=0.0
        depoqicumtot=0.0
        depoqscumtot=0.0
        depoqgcumtot=0.0
        depoqhcumtot=0.0
        frzqcicumtot=0.0
        colqcicumtot=0.0
        colqcscumtot=0.0
        colqcgcumtot=0.0
        colqchcumtot=0.0
        frzqrhcumtot=0.0
        colqricumtot=0.0
        colqrscumtot=0.0
        colqrgcumtot=0.0
        colqrhcumtot=0.0

        condtot=0.0
        depotot=0.0
        frztot=0.0
        condqcheattot=0.0
        condqrheattot=0.0
        nuclqiheattot=0.0
        depoqiheattot=0.0
        depoqsheattot=0.0
        depoqgheattot=0.0
        depoqhheattot=0.0
        frzqciheattot=0.0
        colqciheattot=0.0
        colqcsheattot=0.0
        colqcgheattot=0.0
        colqchheattot=0.0
        frzqrhheattot=0.0
        colqriheattot=0.0
        colqrsheattot=0.0
        colqrgheattot=0.0
        colqrhheattot=0.0
        condheattot=0.0
        depoheattot=0.0
        frzheattot=0.0

        gridvoltot=0.0
        gridmasstot=0.0
        chngqwcumtot=0.0

      END IF

      IF(dotimeseries == 1) THEN
        ! Open files for time-dependent output

        CALL gtlfnkey( runname, lfnkey )

        ts_fn_out  = runname(1:lfnkey)//'.heatts'
        lts_fn_out = 7 + lfnkey

        CALL getunit(nunit_ts)

        lts_fn_out = len_trim(ts_fn_out)
        CALL fnversn( ts_fn_out, lts_fn_out )

        OPEN(UNIT=nunit_ts,FILE=trim(ts_fn_out),STATUS='unknown',   &
            FORM='formatted',IOSTAT= istat )

        write(nunit_ts,'(a,a,a,a)')' t, drymass,dryvol,vmass,cmass,rmass,imass,smass,gmass,hmass,',  &
                       'condqc,condqr,nuclqi,depoqi,depoqs,depoqg,depoqh,',     &
                       'frcqci,frzqrh,colqci,colqcs,colqcg,colqch,',           &
                       'colqri,colqrs,colqrg,colqrh'


        CALL gtlfnkey( runname, lfnkey )

        ts_fn_out  = runname(1:lfnkey)//'.coolts'
        lts_fn_out = 7 + lfnkey

        CALL getunit(nunit_ts2)

        lts_fn_out = len_trim(ts_fn_out)
        CALL fnversn( ts_fn_out, lts_fn_out )

        OPEN(UNIT=nunit_ts2,FILE=trim(ts_fn_out),STATUS='unknown',   &
             FORM='formatted',IOSTAT= istat )

        write(nunit_ts2,'(a,a,a)')' t, drymass,dryvol,vmass,cmass,rmass,imass,smass,gmass,hmass,',  &
                         'evapqc,evapqr,sublqi,sublqs,sublqg,sublqh,',     &
                         'meltqi,meltqs,meltqg,meltqh'

      END IF
      ! Open files for vertical profiles

      IF(dobudget == 1) THEN
        CALL gtlfnkey( runname, lfnkey )

        prof_fn_out  = runname(1:lfnkey)//'.heatprof'
        lprof_fn_out = 7 + lfnkey

        CALL getunit(nunit_heatprof)

        lprof_fn_out = len_trim(prof_fn_out)
        CALL fnversn( prof_fn_out, lprof_fn_out )

        OPEN(UNIT=nunit_heatprof,FILE=trim(prof_fn_out),STATUS='unknown',   &
             FORM='formatted',IOSTAT= istat )

        write(nunit_heatprof,'(a,a,a)')' k,zs_agl,condqc,condqr,nuclqi,depoqi,depoqs,depoqg,depoqh,',     &
                         'frzqci,frzqrh,colqci,colqcs,colqcg,colqch,',           &
                         'colqri,colqrs,colqrg,colqrh'

        CALL gtlfnkey( runname, lfnkey )

        prof_fn_out  = runname(1:lfnkey)//'.coolprof'
        lprof_fn_out = 7 + lfnkey

        CALL getunit(nunit_coolprof)

        lprof_fn_out = len_trim(prof_fn_out)
        CALL fnversn( prof_fn_out, lprof_fn_out )

        OPEN(UNIT=nunit_coolprof,FILE=trim(prof_fn_out),STATUS='unknown',   &
             FORM='formatted',IOSTAT= istat )

        write(nunit_coolprof,'(a,a)')' k,zs_agl,evapqc,evapqr,sublqi,sublqs,sublqg,sublqh,',     &
                           'meltqi,meltqs,meltqg,meltqh'

      END IF

      IF(doavg == 1) THEN

        CALL gtlfnkey( runname, lfnkey )

        prof_fn_out  = runname(1:lfnkey)//'.qnzavgprof'
        lprof_fn_out = 7 + lfnkey

        CALL getunit(nunit_qnzavgprof)

        lprof_fn_out = len_trim(prof_fn_out)
        CALL fnversn( prof_fn_out, lprof_fn_out )

        OPEN(UNIT=nunit_qnzavgprof,FILE=trim(prof_fn_out),STATUS='unknown',   &
             FORM='formatted',IOSTAT= istat )

        write(nunit_qnzavgprof,'(a,a,a)')' k,zs_agl,qvhtavg,qchtavg,qrhtavg,qihtavg,qshtavg,qghtavg,qhhtavg,',     &
                       'nchtavg,nrhtavg,nihtavg,nshtavg,nghtavg,nhhtavg,',           &
                       'zrhtavg,zihtavg,zshtavg,zghtavg,zhhtavg'


        CALL gtlfnkey( runname, lfnkey )

        prof_fn_out  = runname(1:lfnkey)//'.mpavgprof'
        lprof_fn_out = 7 + lfnkey

        CALL getunit(nunit_mpavgprof)

        lprof_fn_out = len_trim(prof_fn_out)
        CALL fnversn( prof_fn_out, lprof_fn_out )

        OPEN(UNIT=nunit_mpavgprof,FILE=trim(prof_fn_out),STATUS='unknown',   &
             FORM='formatted',IOSTAT= istat )

        write(nunit_mpavgprof,'(a,a,a,a,a)')' k,zs_agl,N0chtavg,N0rhtavg,N0ihtavg,N0shtavg,N0ghtavg,N0hhtavg,',     &
                       'N0effchtavg,N0effrhtavg,N0effihtavg,N0effshtavg,N0effghtavg,N0effhhtavg,',           &
                       'lamdachtavg,lamdarhtavg,lamdaihtavg,lamdashtavg,lamdaghtavg,lamdahhtavg,',           &
                       'alphachtavg,alpharhtavg,alphaihtavg,alphashtavg,alphaghtavg,alphahhtavg,',           &
                       'dmchtavg,dmrhtavg,dmihtavg,dmshtavg,dmghtavg,dmhhtavg,rhhtavg'

     END IF

   END IF !curtim == tbgn

   IF(curtim >= budgtbgn .and. curtim <= budgtend) THEN

     ! Read in the additional microphysics variables from individual files

     ! Determine if budgeting will be done for smaller time intervals than the history dump
     ! use this if microphysics rate information was dumped at a smaller time interval than the history
     ! dumps.

      nbudgetsteps = INT(tintv_dmpin/tintv)
      print*,'nbudgetsteps = ', nbudgetsteps

      DO nb=1,nbudgetsteps

        IF(dobudget == 1) THEN
        time=curtim+(nb-1)*tintv
        ! Read in fields of mp rates contributing to cooling

        mpcoolid(1) = 'evapqc'
        mpcoolid(2) = 'evapqr'
        mpcoolid(3) = 'sublqi'
        mpcoolid(4) = 'sublqs'
        mpcoolid(5) = 'sublqg'
        mpcoolid(6) = 'sublqh'
        mpcoolid(7) = 'meltqi'
        mpcoolid(8) = 'meltqs'
        mpcoolid(9) = 'meltqg'
        mpcoolid(10) = 'meltqh'

        DO nq=1,10
          CALL readvar3(nx,ny,nz, mpcoolterms(:,:,:,nq), mpcoolid(nq), mpcoolname(nq),     &
                 mpcoolunits(nq), time, outrunname, dirname, hdmpftrailer,3, 0, istatus)
           IF(istatus /= 0) THEN
             WRITE(6,'(a)')                                                  &
             'Error reading microphysics file. Program ARPSMPBUDGET stopped.'
             STOP
           END IF

        END DO

        ! Read in temperature rate-of-change field
        varid = 'mptrat'
!        varid = 'tdiffm'
        CALL readvar3(nx,ny,nz, mptemprate, varid, varname,     &
               varunits, time, outrunname, dirname, hdmpftrailer,3, 0, istatus)

        IF(istatus /= 0) THEN
          WRITE(6,'(a)')                                                  &
          'Error reading temperature rate file. Program will continue.'
        END IF

        ! Read in fields of mp rates contributing to heating

        mpheatid(1)='condqc'
        mpheatid(2)='condqr'
        mpheatid(3)='nuclqi'
        mpheatid(4)='depoqi'
        mpheatid(5)='depoqs'
        mpheatid(6)='depoqg'
        mpheatid(7)='depoqh'
        mpheatid(8)='frzqci'
        mpheatid(9)='colqci'
        mpheatid(10)='colqcs'
        mpheatid(11)='colqcg'
        mpheatid(12)='colqch'
        mpheatid(13)='frzqrh'
        mpheatid(14)='colqri'
        mpheatid(15)='colqrs'
        mpheatid(16)='colqrg'
        mpheatid(17)='colqrh'

        DO nq=1,17

          CALL readvar3(nx,ny,nz,mpheatterms(:,:,:,nq), mpheatid(nq),mpheatname(nq),mpheatunits(nq),time,outrunname,    &
                   dirname,hdmpftrailer,3,0,istatus)
          IF(istatus /= 0 .and. nq /= 3) THEN
            WRITE(6,'(a)')                                                  &
            'Error reading microphysics file. Program ARPSMPBUDGET stopped.'
            STOP
          END IF
        END DO

!        DO nq=1,nqscalar
!          varid = 'n0__'//qnames(nq)
!          varname = qdescp(nq)
!          varunits = 'm^-4'
!          CALL readvar3(nx,ny,nz, N0x(:,:,:,nq), varid, varname,     &
!                 varunits, time, outrunname, dirname,hdmpftrailer, 3, 0, istatus)
!          IF(istatus /= 0) THEN
!            WRITE(6,'(a)')                                                  &
!            'Error reading N0 file. Program will continue.'
!          END IF
!        END DO

        ! Peform budget calculations

        ! For now, ctravgi = 0 (no feature tracking for now)
        ctravgi = 0

        ! Initialize individual time total terms to zero.

        evapqctot=0.0
        evapqrtot=0.0
        sublqitot=0.0
        sublqstot=0.0
        sublqgtot=0.0
        sublqhtot=0.0
        meltqitot=0.0
        meltqstot=0.0
        meltqgtot=0.0
        meltqhtot=0.0
        condqctot=0.0
        condqrtot=0.0
        nuclqitot=0.0
        depoqitot=0.0
        depoqstot=0.0
        depoqgtot=0.0
        depoqhtot=0.0
        frzqcitot=0.0
        colqcitot=0.0
        colqcstot=0.0
        colqcgtot=0.0
        colqchtot=0.0
        frzqrhtot=0.0
        colqritot=0.0
        colqrstot=0.0
        colqrgtot=0.0
        colqrhtot=0.0
        gridvoltot=0.0
        mptchgtot=0.0
        mpQchgtot=0.0
        chngqwtot=0.0

        gridmasstot=0.0
        vmasstot=0.0
        cmasstot=0.0
        rmasstot=0.0
        imasstot=0.0
        smasstot=0.0
        gmasstot=0.0
        hmasstot=0.0

        END IF ! dobudget == 1

        IF(doavg == 1 .and. nb == 1) THEN
          qvhavg=0.0
          qsclhavg=0.0
          Ntxhavg=0.0
          Dmhavg=0.0
          rhhavg=0.0
          alphahavg=0.0
          N0havg=0.0
          lamdahavg=0.0

          rhotot=0.0

          gridcount=0
          qflag=0
          qcount=0

        END IF

        ! Sum up all contributions from microphysics heating and cooling processes as a function
        ! of height (k-index)

        IF(dobudget == 1 .or. doavg == 1 .or. dotimeseries == 1) THEN
        DO k=bkbgn,bkend
          DO i=bibgn,biend
            DO j=bjbgn,bjend

              ! Restrict the budget and averaging calculations to levels below zthresh and
              ! vertical velocity < wthresh
              IF((wopt == 1 .and. zs_agl(i+ctravgi,j,k) <= zthresh .and. ws(i+ctravgi,j,k) <= wthresh)      &
                    .or. (wopt == 2 .and. zs_agl(i+ctravgi,j,k) >= zthresh .and. ws(i+ctravgi,j,k) >= wthresh)) THEN

                gridvoltot(k) = gridvoltot(k) + gridvol(i+ctravgi,j,k)
                gridmasstot(k) = gridmasstot(k) + gridmass(i+ctravgi,j,k)

                !Calculate layer total water mass change.  This is done by
                !Subtracting the present value of qw from the future value
                !at each point, multiplying by the time interval size and grid air mass
                !and summing.  It is done lagged by one time interval compared to the other
                !calculations to avoid having to read two history files at once.

                IF(dobudget == 1 .or. dotimeseries == 1) THEN
                IF(nfile >= 2) THEN
                  DO nq=1,nqscalar
                    chngqwtot(k) = chngqwtot(k) + (qscalar(i+ctravgi,j,k,nq) - temscalar(i+ctravgi,j,k,nq))      &
                                        *tintv*gridmass(i+ctravgi,j,k)
                  END DO
                END IF

                IF(curtim <= (budgtend-tintv_dmpin)) THEN  ! Instantaneous rates are propagated forward in time
                                                   ! so values at last history dump are not used.

                  ! Total cloud evaporation (kg)
                  evapqctot(k) = evapqctot(k) + mpcoolterms(i+ctravgi,j,k,1)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total rain evaporation (kg)
                  evapqrtot(k) = evapqrtot(k) + mpcoolterms(i+ctravgi,j,k,2)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total ice sublimation (kg)
                  sublqitot(k) = sublqitot(k) + mpcoolterms(i+ctravgi,j,k,3)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total snow sublimation (kg)
                  sublqstot(k) = sublqstot(k) + mpcoolterms(i+ctravgi,j,k,4)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total graupel sublimation (kg)
                  sublqgtot(k) = sublqgtot(k) + mpcoolterms(i+ctravgi,j,k,5)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total hail sublimation (kg)
                  sublqhtot(k) = sublqhtot(k) + mpcoolterms(i+ctravgi,j,k,6)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total ice melting (kg)
                  meltqitot(k) = meltqitot(k) + mpcoolterms(i+ctravgi,j,k,7)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total snow melting (kg)
                  meltqstot(k) = meltqstot(k) + mpcoolterms(i+ctravgi,j,k,8)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total graupel melting (kg)
                  meltqgtot(k) = meltqgtot(k) + mpcoolterms(i+ctravgi,j,k,9)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total hail melting (kg)
                  meltqhtot(k) = meltqhtot(k) + mpcoolterms(i+ctravgi,j,k,10)*tintv*gridmass(i+ctravgi,j,k)

                  ! Total cloud condensation (kg)
                  condqctot(k) = condqctot(k) + mpheatterms(i+ctravgi,j,k,1)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total rain condensation (kg)
                  condqrtot(k) = condqrtot(k) + mpheatterms(i+ctravgi,j,k,2)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total cloud ice nucleation (kg)
                  nuclqitot(k) = nuclqitot(k) + mpheatterms(i+ctravgi,j,k,3)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total ice deposition (kg)
                  depoqitot(k) = depoqitot(k) + mpheatterms(i+ctravgi,j,k,4)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total snow deposition (kg)
                  depoqstot(k) = depoqstot(k) + mpheatterms(i+ctravgi,j,k,5)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total graupel deposition (kg)
                  depoqgtot(k) = depoqgtot(k) + mpheatterms(i+ctravgi,j,k,6)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total hail deposition (kg)
                  depoqhtot(k) = depoqhtot(k) + mpheatterms(i+ctravgi,j,k,7)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total cloud ice freezing (kg)
                  frzqcitot(k) = frzqcitot(k) + mpheatterms(i+ctravgi,j,k,8)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total cloud ice collection of cloud water (kg)
                  colqcitot(k) = colqcitot(k) + mpheatterms(i+ctravgi,j,k,9)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total snow collection of cloud water (kg)
                  colqcstot(k) = colqcstot(k) + mpheatterms(i+ctravgi,j,k,10)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total graupel collection of cloud water (kg)
                  colqcgtot(k) = colqcgtot(k) + mpheatterms(i+ctravgi,j,k,11)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total hail collection of cloud water (kg)
                  colqchtot(k) = colqchtot(k) + mpheatterms(i+ctravgi,j,k,12)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total rain freezing (to hail) (kg)
                  frzqrhtot(k) = frzqrhtot(k) + mpheatterms(i+ctravgi,j,k,13)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total ice collection of rain (kg)
                  colqritot(k) = colqritot(k) + mpheatterms(i+ctravgi,j,k,14)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total snow collection of rain (kg)
                  colqrstot(k) = colqrstot(k) + mpheatterms(i+ctravgi,j,k,15)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total graupel collection of rain (kg)
                  colqrgtot(k) = colqrgtot(k) + mpheatterms(i+ctravgi,j,k,16)*tintv*gridmass(i+ctravgi,j,k)
                  ! Total hail collection of rain (kg)
                  colqrhtot(k) = colqrhtot(k) + mpheatterms(i+ctravgi,j,k,17)*tintv*gridmass(i+ctravgi,j,k)

                  ! Total heating/cooling due to microphysics (J)

                  mpQchgtot(k) = mpQchgtot(k) + mptemprate(i+ctravgi,j,k)*tintv*CP*gridmass(i+ctravgi,j,k)

                END IF ! within budgeting time window

                END IF ! dobudget == 1 or dotimeseries == 1

                IF(doavg == 1 .and. nb == 1) THEN
                ! Also calculate the horizontally-averaged values of mixing ratio, number concentration,
                ! and reflectivity for each hydrometeor species.  These are done only for regions where
                ! the mixing ratio of the given species is > 0.

                ! Note, each value of mixing ratio is multiplied by air density to get water content
                ! and then summed up.

                ! First determine if a point contains nonzero q,n,z of a given
                ! hydrometeor and set flag and increment gridcount accordingly

                IF(mphyopt <= 8) THEN
                  DO nq=1,nqscalar
                    IF(qscalar(i+ctravgi,j,k,nq) > 0.0) THEN
                      qflag(nq) = 1
                      qcount(nq) = qcount(nq) + 1
                    ELSE
                      qflag(nq) = 0
                    END IF
                  END DO
                ELSE IF(mphyopt <= 10) THEN
                  DO nq=1,6
                    IF(qscalar(i+ctravgi,j,k,nq) > 0.0 .and. qscalar(i+ctravgi,j,k,nq+6) > 0.0) THEN
                      qflag(nq) = 1
                      qcount(nq) = qcount(nq) + 1
                    ELSE
                      qflag(nq) = 0
                    END IF
                  END DO
                ELSE IF(mphyopt == 11) THEN
                  IF(qscalar(i+ctravgi,j,k,1) > 0.0 .and. qscalar(i+ctravgi,j,k,7) > 0.0) THEN
                    qflag(1) = 1
                    qcount(1) = qcount(1) + 1
                  ELSE
                    qflag(1) = 0
                  END IF

                  DO nq=2,6
                    IF(qscalar(i+ctravgi,j,k,nq) > 0.0 .and. qscalar(i+ctravgi,j,k,nq+6) > 0.0 .and. qscalar(i+ctravgi,j,k,nq+11) > 0.0) THEN
                      qflag(nq) = 1
                      qcount(nq) = qcount(nq) + 1
                    ELSE
                      qflag(nq) = 0
                    END IF
                  END DO
                END IF

                qvhavg(k) = qvhavg(k) + (qvprt(i+ctravgi,j,k)+qvbar(i+ctravgi,j,k))*rho(i+ctravgi,j,k)
                rhotot(k,7) = rhotot(k,7) + rho(i+ctravgi,j,k)

                DO nq=1,nqscalar
                  IF(qflag(nq) == 1) THEN
                    qsclhavg(k,nq) = qsclhavg(k,nq) + qscalar(i+ctravgi,j,k,nq)*rho(i+ctravgi,j,k)
                    rhotot(k,nq) = rhotot(k,nq) + rho(i+ctravgi,j,k)
                    IF(mphyopt == 9 .or. mphyopt == 10 .or. mphyopt == 11) THEN
                      qsclhavg(k,nq+6) = qsclhavg(k,nq+6) + qscalar(i+ctravgi,j,k,nq+6)
                    ELSE IF(mphyopt <= 8) THEN
                      Ntxhavg(k,nq) = Ntxhavg(k,nq) + Ntx(i+ctravgi,j,k,nq)
                    END IF
                    IF(mphyopt == 11 .and. nq /= P_QC) THEN
                      qsclhavg(k,nq+11) = qsclhavg(k,nq+11) + qscalar(i+ctravgi,j,k,nq+11)
                    END IF
                    N0havg(k,nq) = N0havg(k,nq) + N0(i+ctravgi,j,k,nq)
                    lamdahavg(k,nq) = lamdahavg(k,nq) + lamda(i+ctravgi,j,k,nq)
                    alphahavg(k,nq) = alphahavg(k,nq) + alpha(i+ctravgi,j,k,nq)
                    Dmhavg(k,nq) = Dmhavg(k,nq) + Dm(i+ctravgi,j,k,nq)
                  END IF
                END DO

                ! Relative humidity

                rhhavg(k) = rhhavg(k) + rh(i+ctravgi,j,k)

                ! Sum up the layer-total hydrometeor mass

                vmasstot(k) = vmasstot(k)+(qvprt(i+ctravgi,j,k)+qvbar(i+ctravgi,j,k))*gridmass(i+ctravgi,j,k)
                IF(P_QC > 0) THEN
                  cmasstot(k) = cmasstot(k)+qscalar(i+ctravgi,j,k,P_QC)*gridmass(i+ctravgi,j,k)
                END IF
                IF(P_QR > 0) THEN
                  rmasstot(k) = rmasstot(k)+qscalar(i+ctravgi,j,k,P_QR)*gridmass(i+ctravgi,j,k)
                END IF
                IF(P_QI > 0) THEN
                  imasstot(k) = imasstot(k)+qscalar(i+ctravgi,j,k,P_QI)*gridmass(i+ctravgi,j,k)
                END IF
                IF(P_QS > 0) THEN
                  smasstot(k) = smasstot(k)+qscalar(i+ctravgi,j,k,P_QS)*gridmass(i+ctravgi,j,k)
                END IF
                IF(P_QG > 0) THEN
                  gmasstot(k) = gmasstot(k)+qscalar(i+ctravgi,j,k,P_QG)*gridmass(i+ctravgi,j,k)
                END IF
                IF(P_QH > 0) THEN
                  hmasstot(k) = hmasstot(k)+qscalar(i+ctravgi,j,k,P_QH)*gridmass(i+ctravgi,j,k)
                END IF

                gridcount=gridcount+1
                END IF ! doavg == 1

              END IF ! within spatial criteria

            END DO ! j loop
          END DO ! i loop

          IF(doavg == 1 .and. nb == 1) THEN
          ! Complete horizontal average calculations for hydrometeors

          IF(rhotot(k,7) /= 0.0) THEN
            qvhavg(k) = qvhavg(k)/rhotot(k,7)
            numtimes1(k,7) = numtimes1(k,7) + 1
          END IF

          DO nq=1,nqscalar
            IF(qcount(nq) /= 0) THEN
              qsclhavg(k,nq) = qsclhavg(k,nq)/rhotot(k,nq)
              numtimes1(k,nq) = numtimes1(k,nq) + 1
              IF(mphyopt <= 8) THEN
                Ntxhavg(k,nq) = Ntxhavg(k,nq)/float(qcount(nq))
              ELSE
                qsclhavg(k,nq+6)=qsclhavg(k,nq+6)/float(qcount(nq))
              END IF
              IF(mphyopt == 11 .and. nq /= P_QC) THEN
                qsclhavg(k,nq+11)=qsclhavg(k,nq+11)/float(qcount(nq))
              END IF
              N0havg(k,nq) = N0havg(k,nq)/dble(qcount(nq))
              lamdahavg(k,nq) = lamdahavg(k,nq)/dble(qcount(nq))
              alphahavg(k,nq) = alphahavg(k,nq)/dble(qcount(nq))
              Dmhavg(k,nq) = Dmhavg(k,nq)/float(qcount(nq))
            END IF
          END DO

          IF (gridcount /= 0) THEN
            rhhavg(k)=rhhavg(k)/gridcount
          ELSE
            rhhavg(k)=0
          END IF

          gridcount = 0
          qcount = 0
          numtimes2(k) = numtimes2(k) + 1

          END IF ! doavg == 1
        END DO ! k loop

        END IF ! dobudget == 1 or doavg == 1 or dotimeseries == 1
        ! Output time-dependent totals to a file

        IF(dotimeseries == 1 .and. curtim <= (budgtend-tintv)) THEN

        kgtokt = 1.0e-6

        write(nunit_ts,'f8.1,26f24.9')curtim,SUM(gridmasstot)*kgtokt,SUM(gridvoltot)*1.0e-9,SUM(vmasstot)*kgtokt,SUM(cmasstot)*kgtokt,  &
                                    SUM(rmasstot)*kgtokt,SUM(imasstot)*kgtokt,SUM(smasstot)*kgtokt,  &
                                    SUM(gmasstot)*kgtokt,SUM(hmasstot)*kgtokt,SUM(condqctot)*kgtokt,SUM(condqrtot)*kgtokt,SUM(nuclqitot)*kgtokt,      &
                                    SUM(depoqitot)*kgtokt,SUM(depoqstot)*kgtokt,SUM(depoqgtot)*kgtokt,SUM(depoqhtot)*kgtokt,SUM(frzqcitot)*kgtokt,    &
                                    SUM(frzqrhtot)*kgtokt,SUM(colqcitot)*kgtokt,SUM(colqcstot)*kgtokt,SUM(colqcgtot)*kgtokt,SUM(colqchtot)*kgtokt,    &
                                    SUM(colqritot)*kgtokt,SUM(colqrstot)*kgtokt,SUM(colqrgtot)*kgtokt,SUM(colqrhtot)*kgtokt

        write(nunit_ts2,'f8.1,19f24.9')curtim,SUM(gridmasstot)*kgtokt,SUM(gridvoltot)*1.0e-9,SUM(vmasstot)*kgtokt,SUM(cmasstot)*kgtokt,            &
                                     SUM(rmasstot)*kgtokt,SUM(imasstot)*kgtokt,SUM(smasstot)*kgtokt, &
                                     SUM(gmasstot)*kgtokt,SUM(hmasstot)*kgtokt,SUM(evapqctot)*kgtokt,SUM(evapqrtot)*kgtokt,SUM(sublqitot)*kgtokt,       &
                                     SUM(sublqstot)*kgtokt,SUM(sublqgtot)*kgtokt,SUM(sublqhtot)*kgtokt,SUM(meltqitot)*kgtokt,SUM(meltqstot)*kgtokt,     &
                                     SUM(meltqgtot)*kgtokt,SUM(meltqhtot)*kgtokt

        END IF

        ! Add subtotal values to cumulative total values (for each k index)

        DO k=bkbgn,bkend

          !cmasscumtot(k)=cmasscumtot(k)+cmasstot(k)
          !rmasscumtot(k)=rmasscumtot(k)+rmasstot(k)
          !imasscumtot(k)=imasscumtot(k)+imasstot(k)
          !smasscumtot(k)=smasscumtot(k)+smasstot(k)
          !gmasscumtot(k)=gmasscumtot(k)+gmasstot(k)
          !hmasscumtot(k)=hmasscumtot(k)+hmasstot(k)

          IF(doavg == 1 .and. nb == 1) THEN
          qvhtavg(k)=qvhtavg(k)+qvhavg(k)
          DO nq=1,nqscalar
            qsclhtavg(k,nq) = qsclhtavg(k,nq)+qsclhavg(k,nq)
            IF(mphyopt <= 8) THEN
              Ntxhtavg(k,nq) = Ntxhtavg(k,nq)+Ntxhavg(k,nq)
            ELSE
              qsclhtavg(k,nq+6) = qsclhtavg(k,nq+6)+qsclhavg(k,nq+6)
            END IF
            IF(mphyopt == 11 .and. nq /= P_QC) THEN
              qsclhtavg(k,nq+11) = qsclhtavg(k,nq+11)+qsclhavg(k,nq+11)
            END IF
            N0htavg(k,nq)=N0htavg(k,nq)+N0havg(k,nq)
            lamdahtavg(k,nq)=lamdahtavg(k,nq)+lamdahavg(k,nq)
            alphahtavg(k,nq)=alphahtavg(k,nq)+alphahavg(k,nq)
            Dmhtavg(k,nq)=Dmhtavg(k,nq)+Dmhavg(k,nq)
          END DO

          rhhtavg(k)=rhhtavg(k)+rhhavg(k)
          END IF ! doavg == 1

          IF(dobudget == 1) THEN
          IF(curtim <= (budgtend-tintv)) THEN

            evapqccumtot(k)=evapqccumtot(k)+evapqctot(k)
            evapqrcumtot(k)=evapqrcumtot(k)+evapqrtot(k)
            sublqicumtot(k)=sublqicumtot(k)+sublqitot(k)
            sublqscumtot(k)=sublqscumtot(k)+sublqstot(k)
            sublqgcumtot(k)=sublqgcumtot(k)+sublqgtot(k)
            sublqhcumtot(k)=sublqhcumtot(k)+sublqhtot(k)
            meltqicumtot(k)=meltqicumtot(k)+meltqitot(k)
            meltqscumtot(k)=meltqscumtot(k)+meltqstot(k)
            meltqgcumtot(k)=meltqgcumtot(k)+meltqgtot(k)
            meltqhcumtot(k)=meltqhcumtot(k)+meltqhtot(k)
            condqccumtot(k)=condqccumtot(k)+condqctot(k)
            condqrcumtot(k)=condqrcumtot(k)+condqrtot(k)
            nuclqicumtot(k)=nuclqicumtot(k)+nuclqitot(k)
            depoqicumtot(k)=depoqicumtot(k)+depoqitot(k)
            depoqscumtot(k)=depoqscumtot(k)+depoqstot(k)
            depoqgcumtot(k)=depoqgcumtot(k)+depoqgtot(k)
            depoqhcumtot(k)=depoqhcumtot(k)+depoqhtot(k)
            frzqcicumtot(k)=frzqcicumtot(k)+frzqcitot(k)
            colqcicumtot(k)=colqcicumtot(k)+colqcitot(k)
            colqcscumtot(k)=colqcscumtot(k)+colqcstot(k)
            colqcgcumtot(k)=colqcgcumtot(k)+colqcgtot(k)
            colqchcumtot(k)=colqchcumtot(k)+colqchtot(k)
            frzqrhcumtot(k)=frzqrhcumtot(k)+frzqrhtot(k)
            colqricumtot(k)=colqricumtot(k)+colqritot(k)
            colqrscumtot(k)=colqrscumtot(k)+colqrstot(k)
            colqrgcumtot(k)=colqrgcumtot(k)+colqrgtot(k)
            colqrhcumtot(k)=colqrhcumtot(k)+colqrhtot(k)
            mpQchgcumtot(k)=mpQchgcumtot(k)+mpQchgtot(k)

          END IF

          chngqwcumtot(k)=chngqwcumtot(k)+chngqwtot(k)
          END IF ! dobudget == 1

        END DO ! summing k-loop

      END DO ! budget sub-interval

      ! Store new qscalar data in temscalar array for use in next time interval
      DO nq=1,nqscalar
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              temscalar(i,j,k,nq) = qscalar(i,j,k,nq)
            END DO
          END DO
        END DO
      END DO

   END IF ! budgeting time window

   IF(curtim == budgtend) THEN

      IF(doavg == 1) THEN
      ! Complete time-averaging calculations

      DO k=bkbgn,bkend

          IF(numtimes1(k,7) > 0) THEN
            qvhtavg(k)=qvhtavg(k)/numtimes1(k,7)
          END IF

          DO nq=1,nqscalar
            IF(numtimes1(k,nq) > 0) THEN
              qsclhtavg(k,nq)=qsclhtavg(k,nq)/float(numtimes1(k,nq))
              IF(mphyopt <= 8) THEN
                Ntxhtavg(k,nq)=Ntxhtavg(k,nq)/float(numtimes1(k,nq))
              ELSE
                qsclhtavg(k,nq+6)=qsclhtavg(k,nq+6)/float(numtimes1(k,nq))
              END IF
              IF(mphyopt == 11 .and. nq /= P_QC) THEN
                qsclhtavg(k,nq+11)=qsclhtavg(k,nq+11)/float(numtimes1(k,nq))
              END IF
              N0htavg(k,nq)=N0htavg(k,nq)/dble(numtimes1(k,nq))
              lamdahtavg(k,nq)=lamdahtavg(k,nq)/dble(numtimes1(k,nq))
              alphahtavg(k,nq)=alphahtavg(k,nq)/dble(numtimes1(k,nq))
              Dmhtavg(k,nq)=Dmhtavg(k,nq)/float(numtimes1(k,nq))
            END IF
          END DO
          rhhtavg(k)=rhhtavg(k)/float(numtimes2(k))

      END DO

      END IF ! doavg == 1
      ! Write out vertical profiles of horizontally-summed totals of
      ! the various quantities and time and horizontally-averaged
      ! vertical profiles of mixing ratio, number concentration, and reflectivity

      DO k=bkend,bkbgn,-1

        IF(dobudget == 1) THEN
          write(nunit_heatprof,'i3,f8.1,23f24.9') k,zs_agl(2,2,k),condqccumtot(k),     &
              condqrcumtot(k),nuclqicumtot(k),depoqicumtot(k),depoqscumtot(k),depoqgcumtot(k), &
              depoqhcumtot(k),frzqcicumtot(k),frzqrhcumtot(k),colqcicumtot(k),colqcscumtot(k), &
              colqcgcumtot(k),colqchcumtot(k),colqricumtot(k),colqrscumtot(k),colqrgcumtot(k), &
              colqrhcumtot(k)

          write(nunit_coolprof,'i3,f8.1,16f24.9') k,zs_agl(2,2,k),evapqccumtot(k),     &
              evapqrcumtot(k),sublqicumtot(k),sublqscumtot(k),sublqgcumtot(k),sublqhcumtot(k), &
              meltqicumtot(k),meltqscumtot(k),meltqgcumtot(k),meltqhcumtot(k)
        END IF

        IF(doavg == 1) THEN

          IF(mphyopt == 2) THEN
            write(nunit_qnzavgprof,'i3,f8.1,7f17.9,6f21.4,5f17.9') k,zs_agl(2,2,k),qvhtavg(k),     &
                (qsclhtavg(k,nq),nq=1,4),0.0,qsclhtavg(k,5),(Ntxhtavg(k,nq),nq=1,4),0.0,Ntxhtavg(k,5),  &
                0.0,0.0,0.0,0.0,0.0
          ELSE IF(mphyopt == 8) THEN
            write(nunit_qnzavgprof,'i3,f8.1,7f17.9,6f21.4,5f17.9') k,zs_agl(2,2,k),qvhtavg(k),     &
                (qsclhtavg(k,nq),nq=1,nqscalar),(Ntxhtavg(k,nq),nq=1,nqscalar),  &
                0.0,0.0,0.0,0.0,0.0
          ELSE IF(mphyopt == 9 .or. mphyopt == 10) THEN
            write(nunit_qnzavgprof,'i3,f8.1,7f17.9,6f21.4,5f17.9') k,zs_agl(2,2,k),qvhtavg(k),     &
                (qsclhtavg(k,nq),nq=1,nscalar),  &
                0.0,0.0,0.0,0.0,0.0
          ELSE
            write(nunit_qnzavgprof,'i3,f8.1,7f17.9,6f21.4,5f17.9') k,zs_agl(2,2,k),qvhtavg(k),     &
                (qsclhtavg(k,nq),nq=1,nscalar)
          END IF

          IF(mphyopt == 2) THEN
            write(nunit_mpavgprof,'i3,f8.1,12d48.2,12d17.9,7f17.9') k,zs_agl(2,2,k),(N0htavg(k,nq),nq=1,4),0.0,N0htavg(k,5),  &
                0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,    &
                (lamdahtavg(k,nq),nq=1,4),0.d0, lamdahtavg(k,5),    &
                (alphahtavg(k,nq),nq=1,4),0.d0, alphahtavg(k,5),    &
                (Dmhtavg(k,nq),nq=1,4),0.0,Dmhtavg(k,5),rhhtavg(k)
          ELSE IF(mphyopt >= 8) THEN
            write(nunit_mpavgprof,'i3,f8.1,12d48.2,12d17.9,7f17.9') k,zs_agl(2,2,k),(N0htavg(k,nq),nq=1,nqscalar),  &
                0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,    &
                (lamdahtavg(k,nq),nq=1,nqscalar),    &
                (alphahtavg(k,nq),nq=1,nqscalar),    &
                (Dmhtavg(k,nq),nq=1,nqscalar),rhhtavg(k)
          END IF
        END IF
      END DO

      IF(dobudget == 1) THEN
      DO k=bkbgn,bkend

        ! Total evaporation

        evaptot(k) = evapqccumtot(k)+evapqrcumtot(k)

        ! Total sublimation

        subltot(k) = sublqicumtot(k)+sublqscumtot(k)+sublqgcumtot(k)+sublqhcumtot(k)

        ! Total melting

        melttot(k) = meltqicumtot(k)+meltqscumtot(k)+meltqgcumtot(k)+meltqhcumtot(k)

        ! Total cooling for each process

        evapqccooltot(k) = LV*evapqccumtot(k)
        evapqrcooltot(k) = LV*evapqrcumtot(k)
        sublqicooltot(k) = LS*sublqicumtot(k)
        sublqscooltot(k) = LS*sublqscumtot(k)
        sublqgcooltot(k) = LS*sublqgcumtot(k)
        sublqhcooltot(k) = LS*sublqhcumtot(k)
        meltqicooltot(k) = LF*meltqicumtot(k)
        meltqscooltot(k) = LF*meltqscumtot(k)
        meltqgcooltot(k) = LF*meltqgcumtot(k)
        meltqhcooltot(k) = LF*meltqhcumtot(k)

        ! Total evaporative cooling

        evapcooltot(k) = evapqccooltot(k)+evapqrcooltot(k)

        ! Total sublimation cooling

        sublcooltot(k) = sublqicooltot(k)+sublqscooltot(k)+sublqgcooltot(k)+sublqhcooltot(k)

        ! Total melting cooling

        meltcooltot(k) = meltqicooltot(k)+meltqscooltot(k)+meltqgcooltot(k)+meltqhcooltot(k)

        ! Total cooling due to all the above processes

        mpcooltot(k) = evapcooltot(k)+sublcooltot(k)+meltcooltot(k)

        ! Total condensation

        condtot(k)=condqccumtot(k)+condqrcumtot(k)

        ! Total deposition

        depotot(k)=nuclqicumtot(k)+depoqicumtot(k)+depoqscumtot(k)+depoqgcumtot(k)+depoqhcumtot(k)

        ! Total freezing

        frztot(k)=frzqcicumtot(k)+frzqrhcumtot(k)+colqcicumtot(k)+colqcscumtot(k)+colqcgcumtot(k)+colqchcumtot(k)+colqricumtot(k)+colqrscumtot(k)+colqrgcumtot(k)+colqrhcumtot(k)

        ! Total heating for each process

        condqcheattot(k)=LV*condqccumtot(k)
        condqrheattot(k)=LV*condqrcumtot(k)
        nuclqiheattot(k)=LS*nuclqicumtot(k)
        depoqiheattot(k)=LS*depoqicumtot(k)
        depoqsheattot(k)=LS*depoqscumtot(k)
        depoqgheattot(k)=LS*depoqgcumtot(k)
        depoqhheattot(k)=LS*depoqhcumtot(k)
        frzqciheattot(k)=LF*frzqcicumtot(k)
        frzqrhheattot(k)=LF*frzqrhcumtot(k)
        colqciheattot(k)=LF*colqcicumtot(k)
        colqcsheattot(k)=LF*colqcscumtot(k)
        colqcgheattot(k)=LF*colqcgcumtot(k)
        colqchheattot(k)=LF*colqchcumtot(k)
        colqriheattot(k)=LF*colqricumtot(k)
        colqrsheattot(k)=LF*colqrscumtot(k)
        colqrgheattot(k)=LF*colqrgcumtot(k)
        colqrhheattot(k)=LF*colqrhcumtot(k)

        ! Total condensational heating

        condheattot(k)=condqcheattot(k)+condqrheattot(k)

        ! Total depositional heating

        depoheattot(k)=nuclqiheattot(k)+depoqiheattot(k)+depoqsheattot(k)+depoqgheattot(k)+depoqhheattot(k)

        ! Total freezing heating

        frzheattot(k)=frzqciheattot(k)+frzqrhheattot(k)+colqciheattot(k)+colqcsheattot(k)+colqcgheattot(k)     &
                       +colqchheattot(k)+colqriheattot(k)+colqrsheattot(k)+colqrgheattot(k)+colqrhheattot(k)

        ! Total heating due to microphysics

        mpheattot(k) = condheattot(k)+depoheattot(k)+frzheattot(k)

      END DO

      ! Now write out a file containing the budget information
      ! Units are also converted to kilotons (kt) and Gigajoules (GJ)
      ! for readability

      !If normalization is desired, divide all processes in each layer by their sum total
      !for that type of process (i.e., (heating process)/(sum total heating))

      IF(normalize_opt == 1) THEN
        kgtokt = 1.0
        JtoGJ = 1.0

        DO k=bkbgn,bkend
          evapqccumtot(k)=evapqccumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          evapqrcumtot(k)=evapqrcumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          sublqicumtot(k)=sublqicumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          sublqscumtot(k)=sublqscumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          sublqgcumtot(k)=sublqgcumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          sublqhcumtot(k)=sublqhcumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          meltqicumtot(k)=meltqicumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          meltqscumtot(k)=meltqscumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          meltqgcumtot(k)=meltqgcumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          meltqhcumtot(k)=meltqhcumtot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))

          evaptot(k)=evaptot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          subltot(k)=subltot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))
          melttot(k)=melttot(k)/(SUM(evaptot)+SUM(subltot)+SUM(melttot))

          evapqccooltot(k)=evapqccooltot(k)/SUM(mpcooltot)
          evapqrcooltot(k)=evapqrcooltot(k)/SUM(mpcooltot)
          sublqicooltot(k)=sublqicooltot(k)/SUM(mpcooltot)
          sublqscooltot(k)=sublqscooltot(k)/SUM(mpcooltot)
          sublqgcooltot(k)=sublqgcooltot(k)/SUM(mpcooltot)
          sublqhcooltot(k)=sublqhcooltot(k)/SUM(mpcooltot)
          meltqicooltot(k)=meltqicooltot(k)/SUM(mpcooltot)
          meltqscooltot(k)=meltqscooltot(k)/SUM(mpcooltot)
          meltqgcooltot(k)=meltqgcooltot(k)/SUM(mpcooltot)
          meltqhcooltot(k)=meltqhcooltot(k)/SUM(mpcooltot)

          evapcooltot(k)=evapcooltot(k)/SUM(mpcooltot)
          sublcooltot(k)=sublcooltot(k)/SUM(mpcooltot)
          meltcooltot(k)=meltcooltot(k)/SUM(mpcooltot)

          condqccumtot(k)=condqccumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          condqrcumtot(k)=condqrcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          nuclqicumtot(k)=nuclqicumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          depoqicumtot(k)=depoqicumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          depoqscumtot(k)=depoqscumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          depoqgcumtot(k)=depoqgcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          depoqhcumtot(k)=depoqhcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          frzqcicumtot(k)=frzqcicumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqcicumtot(k)=colqcicumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqcscumtot(k)=colqcscumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqcgcumtot(k)=colqcgcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqchcumtot(k)=colqchcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          frzqrhcumtot(k)=frzqrhcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqricumtot(k)=colqricumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqrscumtot(k)=colqrscumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqrgcumtot(k)=colqrgcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          colqrhcumtot(k)=colqrhcumtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))

          condtot(k)=condtot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          depotot(k)=depotot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))
          frztot(k)=frztot(k)/(SUM(condtot)+SUM(depotot)+SUM(frztot))

          condqcheattot(k)=condqcheattot(k)/SUM(mpheattot)
          condqrheattot(k)=condqrheattot(k)/SUM(mpheattot)
          nuclqiheattot(k)=nuclqiheattot(k)/SUM(mpheattot)
          depoqiheattot(k)=depoqiheattot(k)/SUM(mpheattot)
          depoqsheattot(k)=depoqsheattot(k)/SUM(mpheattot)
          depoqgheattot(k)=depoqgheattot(k)/SUM(mpheattot)
          depoqhheattot(k)=depoqhheattot(k)/SUM(mpheattot)
          frzqciheattot(k)=frzqciheattot(k)/SUM(mpheattot)
          colqciheattot(k)=colqciheattot(k)/SUM(mpheattot)
          colqcsheattot(k)=colqcsheattot(k)/SUM(mpheattot)
          colqcgheattot(k)=colqcgheattot(k)/SUM(mpheattot)
          colqchheattot(k)=colqchheattot(k)/SUM(mpheattot)
          frzqrhheattot(k)=frzqrhheattot(k)/SUM(mpheattot)
          colqriheattot(k)=colqriheattot(k)/SUM(mpheattot)
          colqrsheattot(k)=colqrsheattot(k)/SUM(mpheattot)
          colqrgheattot(k)=colqrgheattot(k)/SUM(mpheattot)
          colqrhheattot(k)=colqrhheattot(k)/SUM(mpheattot)

          condheattot(k)=condheattot(k)/SUM(mpheattot)
          depoheattot(k)=depoheattot(k)/SUM(mpheattot)
          frzheattot(k)=frzheattot(k)/SUM(mpheattot)

        END DO
      ELSE
        kgtokt = 1.0e-6
        JtoGJ = 1.0e-9
      END IF

      CALL gtlfnkey(runname, lfnkey)

      budgetfn  = runname(1:lfnkey)//'.budget'
      lbudgetfn = 7 + lfnkey

      write(6,'(1x,a,a,a/,1x,a)')                                     &
      'Check to see if file ',budgetfn(1:lbudgetfn),' already exists.',   &
      'If so, append a version number to the filename.'

      CALL fnversn( budgetfn, lbudgetfn )
      CALL getunit( nchbudget )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(budgetfn(1:lbudgetfn), '-F f77 -N ieee', ierr)

      OPEN(nchbudget,form='formatted',status='new',                   &
                file=budgetfn(1:lbudgetfn),iostat=istat)

      !IF( istat.ne.0) then
      !  print*,'Error opening file ',budgetfn(1:lbudgetfn),', job aborted.'
      !  STOP
      !ENDIF

      write(nchbudget,'(a)') 'Total subdomain budget for heating/cooling due to microphysics'
      IF(normalize_opt == 1) THEN
        write(nchbudget,'(a)') 'The following are summed values normalized by the total cooling or heating'
        write(nchbudget,'(a)') 'in all scalar grid boxes where vertical'
        write(nchbudget,'(a,f4.2,a)') 'velocity is less than ',wthresh, ' m/s and the height is less than '
        write(nchbudget,'(f8.2,a)') zthresh, ' m'
        write(nchbudget,'(a)') 'The following values therefore are unitless'
      ELSE
        write(nchbudget,'(a)') 'The following are summed values'
        write(nchbudget,'(a)') 'in all scalar grid boxes where vertical'
        write(nchbudget,'(a,f4.2,a)') 'velocity is less than ',wthresh, ' m/s and the height is less than '
        write(nchbudget,'(f8.2,a)') zthresh, ' m'
        write(nchbudget,'(a)') ''
      END IF
      write(nchbudget,'(a,6i3)') 'ibgn,iend,jbgn,jend,kbgn,kend',bibgn,biend,bjbgn,bjend,bkbgn,bkend
      write(nchbudget,*) ''
      write(nchbudget,'(a)') 'Microphysical processes contributing to cooling'
      write(nchbudget,1000) 'Cloud evaporation (kt): ',-SUM(evapqccumtot)*kgtokt
      write(nchbudget,1000) 'Rain evaporation (kt): ',-SUM(evapqrcumtot)*kgtokt
      write(nchbudget,1000) 'Ice sublimation (kt): ',-SUM(sublqicumtot)*kgtokt
      write(nchbudget,1000) 'Snow sublimation (kt): ',-SUM(sublqscumtot)*kgtokt
      write(nchbudget,1000) 'Graupel sublimation (kt): ',-SUM(sublqgcumtot)*kgtokt
      write(nchbudget,1000) 'Hail sublimation (kt): ',-SUM(sublqhcumtot)*kgtokt
      write(nchbudget,1000) 'Ice melting (kt): ',-SUM(meltqicumtot)*kgtokt
      write(nchbudget,1000) 'Snow melting (kt): ',-SUM(meltqscumtot)*kgtokt
      write(nchbudget,1000) 'Graupel melting (kt): ',-SUM(meltqgcumtot)*kgtokt
      write(nchbudget,1000) 'Hail melting (kt): ',-SUM(meltqhcumtot)*kgtokt
      write(nchbudget,*) ''
      write(nchbudget,1000) 'Total evaporation (kt): ',-SUM(evaptot)*kgtokt
      write(nchbudget,1000) 'Total sublimation (kt): ',-SUM(subltot)*kgtokt
      write(nchbudget,1000) 'Total melting (kt): ',-SUM(melttot)*kgtokt
      write(nchbudget,*) ''
      write(nchbudget,'(a)') 'Total cooling due to above processes'
      write(nchbudget,1000) 'Cloud evaporation (GJ): ',-SUM(evapqccooltot)*JtoGJ
      write(nchbudget,1000) 'Rain evaporation (GJ): ',-SUM(evapqrcooltot)*JtoGJ
      write(nchbudget,1000) 'Ice sublimation (GJ): ',-SUM(sublqicooltot)*JtoGJ
      write(nchbudget,1000) 'Snow sublimation (GJ): ',-SUM(sublqscooltot)*JtoGJ
      write(nchbudget,1000) 'Graupel sublimation (GJ): ',-SUM(sublqgcooltot)*JtoGJ
      write(nchbudget,1000) 'Hail sublimation (GJ): ',-SUM(sublqhcooltot)*JtoGJ
      write(nchbudget,1000) 'Ice melting (GJ): ',-SUM(meltqicooltot)*JtoGJ
      write(nchbudget,1000) 'Snow melting (GJ): ',-SUM(meltqscooltot)*JtoGJ
      write(nchbudget,1000) 'Graupel melting (GJ): ',-SUM(meltqgcooltot)*JtoGJ
      write(nchbudget,1000) 'Hail melting (GJ): ',-SUM(meltqhcooltot)*JtoGJ
      write(nchbudget,*) ''
      write(nchbudget,'(a)') 'Microphysical processes contributing to heating'
      write(nchbudget,1000) 'Cloud condensation (kt): ',SUM(condqccumtot)*kgtokt
      write(nchbudget,1000) 'Rain condensation (kt): ',SUM(condqrcumtot)*kgtokt
      write(nchbudget,1000) 'Ice nucleation (kt): ',SUM(nuclqicumtot)*kgtokt
      write(nchbudget,1000) 'Ice deposition (kt): ',SUM(depoqicumtot)*kgtokt
      write(nchbudget,1000) 'Snow deposition (kt): ',SUM(depoqscumtot)*kgtokt
      write(nchbudget,1000) 'Graupel deposition (kt): ',SUM(depoqgcumtot)*kgtokt
      write(nchbudget,1000) 'Hail deposition (kt): ',SUM(depoqhcumtot)*kgtokt
      write(nchbudget,1000) 'Cloud-to-ice freezing (kt): ',SUM(frzqcicumtot)*kgtokt
      write(nchbudget,1000) 'Ice collection of cloud (kt): ',SUM(colqcicumtot)*kgtokt
      write(nchbudget,1000) 'Snow collection of cloud (kt): ',SUM(colqcscumtot)*kgtokt
      write(nchbudget,1000) 'Graupel collection of cloud (kt): ',SUM(colqcgcumtot)*kgtokt
      write(nchbudget,1000) 'Hail collection of cloud (kt): ',SUM(colqchcumtot)*kgtokt
      write(nchbudget,1000) 'Rain-to-hail freezing (kt): ',SUM(frzqrhcumtot)*kgtokt
      write(nchbudget,1000) 'Ice collection of rain (kt): ',SUM(colqricumtot)*kgtokt
      write(nchbudget,1000) 'Snow collection of rain (kt): ',SUM(colqrscumtot)*kgtokt
      write(nchbudget,1000) 'Graupel collection of rain (kt): ',SUM(colqrgcumtot)*kgtokt
      write(nchbudget,1000) 'Hail collection of rain (kt): ',SUM(colqrhcumtot)*kgtokt
      write(nchbudget,*) ''
      write(nchbudget,1000) 'Total condensation (kt): ',SUM(condtot)*kgtokt
      write(nchbudget,1000) 'Total deposition (kt): ',SUM(depotot)*kgtokt
      write(nchbudget,1000) 'Total freezing (kt): ',SUM(frztot)*kgtokt
      write(nchbudget,*) ''
      write(nchbudget,'(a)') 'Total heating due to above processes'
      write(nchbudget,1000) 'Cloud condensation (GJ): ',SUM(condqcheattot)*JtoGJ
      write(nchbudget,1000) 'Rain condensation (GJ): ',SUM(condqrheattot)*JtoGJ
      write(nchbudget,1000) 'Ice nucleation (GJ): ',SUM(nuclqiheattot)*JtoGJ
      write(nchbudget,1000) 'Ice deposition (GJ): ',SUM(depoqiheattot)*JtoGJ
      write(nchbudget,1000) 'Snow deposition (GJ): ',SUM(depoqsheattot)*JtoGJ
      write(nchbudget,1000) 'Graupel deposition (GJ): ',SUM(depoqgheattot)*JtoGJ
      write(nchbudget,1000) 'Hail deposition (GJ): ',SUM(depoqhheattot)*JtoGJ
      write(nchbudget,1000) 'Cloud-to-ice freezing (GJ): ',SUM(frzqciheattot)*JtoGJ
      write(nchbudget,1000) 'Ice collection of cloud (GJ): ',SUM(colqciheattot)*JtoGJ
      write(nchbudget,1000) 'Snow collection of cloud (GJ): ',SUM(colqcsheattot)*JtoGJ
      write(nchbudget,1000) 'Graupel collection of cloud (GJ): ',SUM(colqcgheattot)*JtoGJ
      write(nchbudget,1000) 'Hail collection of cloud (GJ): ',SUM(colqchheattot)*JtoGJ
      write(nchbudget,1000) 'Rain-to-hail freezing (GJ): ',SUM(frzqrhheattot)*JtoGJ
      write(nchbudget,1000) 'Ice collection of rain (GJ): ',SUM(colqriheattot)*JtoGJ
      write(nchbudget,1000) 'Snow collection of rain (GJ): ',SUM(colqrsheattot)*JtoGJ
      write(nchbudget,1000) 'Graupel collection of rain (GJ): ',SUM(colqrgheattot)*JtoGJ
      write(nchbudget,1000) 'Hail collection of rain (GJ): ',SUM(colqrhheattot)*JtoGJ
      write(nchbudget,*) ''
      write(nchbudget,1000) 'Total cooling due to evaporation (GJ): ',-SUM(evapcooltot)*JtoGJ
      write(nchbudget,1000) 'Total cooling due to sublimation (GJ): ',-SUM(sublcooltot)*JtoGJ
      write(nchbudget,1000) 'Total cooling due to melting (GJ): ',-SUM(meltcooltot)*JtoGJ
      write(nchbudget,1000) 'Total cooling due to microphysics (GJ): ',-SUM(mpcooltot)*JtoGJ
      write(nchbudget,1000) 'Total heating due to condensation (GJ): ',SUM(condheattot)*JtoGJ
      write(nchbudget,1000) 'Total heating due to deposition (GJ): ',SUM(depoheattot)*JtoGJ
      write(nchbudget,1000) 'Total heating due to freezing (GJ): ',SUM(frzheattot)*JtoGJ
      write(nchbudget,1000) 'Total heating due to microphysics (GJ): ',SUM(mpheattot)*JtoGJ
      write(nchbudget,1000) 'Total thermal energy change due to MP(GJ, actual): ',SUM(mpQchgcumtot)*JtoGJ
      write(nchbudget,1000) 'Total thermal energy change due to MP(GJ, calc): ',(SUM(mpheattot)-SUM(mpcooltot))*JtoGJ
      write(nchbudget,*) ''
      write(nchbudget,1000) 'Total hydrometeor mass change (kt): ',SUM(chngqwtot)*1.0e-6

      CLOSE(nchbudget)

      END IF ! dobudget == 1
   END IF !curtim == budgtend

  END DO

!-----------------------------------------------------------------------
!
! End of the program
!
!-----------------------------------------------------------------------

  888 CONTINUE

  WRITE(6,'(1x,/,a,/)') ' ==== Program ARPSMPBUDGET terminated normally. ===='

  STOP

  997   CONTINUE
  WRITE(6,'(1x,a,i2,/1x,a)')                                            &
      'Data read was unsuccessful. ireturn =', ireturn,                 &
      'Job stopped.'
  STOP

1000 FORMAT(a50,f24.9)

END PROGRAM arpsmpbudget
