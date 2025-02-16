PROGRAM arpscvt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSCVT                    ######
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
!  Program to read history data dump from ARPS and convert it
!  into another format.
!
!  Parameters grdin,basin,mstin,rainin,icein,trbin are read in from
!  the data file itself, therefore are determined internally.
!  Arrays that are not read in retains their initial zero values.
!  These parameters are passed among subroutines through
!  a common block defined in 'indtflg.inc'.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    8/10/1992
!
!  MODIFICATION HISTORY:
!    Consolidated version for ARPS 3.1 release.
!
!    4/13/93 M. Xue.
!    Modified to conform to the new data dump format.
!
!    10/26/93 Y. Liu
!    Added GrADS format.
!
!    10/26/93 Y. Liu
!    Added surface fields.
!
!    Version 1.1. 11/13/94 M. Xue
!    Changed to namelist input format.
!
!    09/25/95 Yuhe Liu
!    Fixed a bug. Previously parameter nhisfile was used before it
!    was read in from input namelist.
!
!    Version 1.2. 12/07/95 (Yuhe Liu)
!
!    12/07/95 (Yuhe Liu)
!    Updated namelist parameters for landout and rainout.
!
!    12/07/95 (Yuhe Liu)
!    Added conversion of GRIB dump format
!
!    3/13/96 (Ming Xue)
!    Added tke. Changed km to kmh and kmv.
!
!    4/30/1997 (Fanyou Kong -- CMRP)
!    Add Vis5D format output
!
!    12/14/1998 (Donghai Wang)
!    Added the snow cover.
!
!    04/17/2000 (Ming Xue)
!    Added an option that allows one to specify input history data
!    at a constant time interval.
!
!    05/19/2000 (Gene Bassett)
!    Converted to F90, creating allocation and main subroutines.
!
!    05/26/2002 (J. Brotzge)
!    Added/modified tsoil/qsoil for new soil scheme.
!
!    09/20/2003 (Y. Wang)
!    Added the capability to convert terrain data, surface data and
!    boundary data etc.
!
!    07/28/2004 (K. W. Thomas)
!    Added "_ready" file capability.
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
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height of the terrain.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)

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

!-----------------------------------------------------------------------
!
! Newly added pointers
!
!-----------------------------------------------------------------------

  REAL,    DIMENSION(:,:), POINTER :: htptr

  INTEGER, DIMENSION(:,:,:), POINTER :: stypptr
  REAL,    DIMENSION(:,:,:), POINTER :: sfrctptr
  INTEGER, DIMENSION(:,:),   POINTER :: vegtypptr
  REAL,    DIMENSION(:,:),   POINTER :: laiptr
  REAL,    DIMENSION(:,:),   POINTER :: roufnsptr
  REAL,    DIMENSION(:,:),   POINTER :: vegptr
  REAL,    DIMENSION(:,:),   POINTER :: ndviptr

  REAL,    DIMENSION(:,:,:),   POINTER :: zpsoilptr
  REAL,    DIMENSION(:,:,:,:), POINTER :: tsoilptr
  REAL,    DIMENSION(:,:,:,:), POINTER :: qsoilptr
  REAL,    DIMENSION(:,:,:),   POINTER :: wetcanpptr
  REAL,    DIMENSION(:,:),     POINTER :: snowdpthptr
  INTEGER, DIMENSION(:,:,:),   POINTER :: soiltypptr

!-----------------------------------------------------------------------
!
!  INTERFACE
!
!-----------------------------------------------------------------------
  INTERFACE
    SUBROUTINE readcvttrn(ternfile,ternfmt,nx,ny,dx,dy,                 &
                    mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,&
                    hterain)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: ternfile    ! Terrain data file name
      INTEGER,          INTENT(IN)  :: ternfmt
      INTEGER,          INTENT(OUT) :: nx          ! Number of grid points in the x-direction
      INTEGER,          INTENT(OUT) :: ny          ! Number of grid points in the y-direction
      REAL,             INTENT(OUT) :: dx          ! Grid interval in x-direction
      REAL,             INTENT(OUT) :: dy          ! Grid interval in y-direction

      INTEGER,          INTENT(OUT) :: mapproj     ! Map projection
      REAL,             INTENT(OUT) :: trulat1     ! 1st real true latitude of map projection
      REAL,             INTENT(OUT) :: trulat2     ! 2nd real true latitude of map projection
      REAL,             INTENT(OUT) :: trulon      ! Real true longitude of map projection
      REAL,             INTENT(OUT) :: sclfct      ! Map scale factor
      REAL,             INTENT(OUT) :: ctrlat      ! Center latitude of the model domain (deg. N)
      REAL,             INTENT(OUT) :: ctrlon      ! Center longitude of the model domain (deg. E)

      REAL,             POINTER     :: hterain(:,:)! Terrain height.
    END SUBROUTINE readcvttrn

    SUBROUTINE readcvtsfc(sfcfile,sfcfmt,nx,ny,nstyps,dx,dy,            &
              mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,      &
              stypin,vtypin,laiin,roufin,vegin,ndviin,                  &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: sfcfile      ! Name of the surface data file
      INTEGER,          INTENT(IN)  :: sfcfmt
      INTEGER,          INTENT(OUT) :: nx           ! Number of grid points in the x-direction
      INTEGER,          INTENT(OUT) :: ny           ! Number of grid points in the y-direction
      INTEGER,          INTENT(OUT) :: nstyps       ! Max number of soil types in a grid box
      REAL,             INTENT(OUT) :: dx
      REAL,             INTENT(OUT) :: dy

      INTEGER,          INTENT(OUT) :: mapproj       ! Map projection
      REAL,             INTENT(OUT) :: trulat1       ! 1st real true latitude of map projection
      REAL,             INTENT(OUT) :: trulat2       ! 2nd real true latitude of map projection
      REAL,             INTENT(OUT) :: trulon        ! Real true longitude of map projection
      REAL,             INTENT(OUT) :: sclfct        ! Map scale factor
      REAL,             INTENT(OUT) :: ctrlat        ! Center latitude of the model domain (deg. N)
      REAL,             INTENT(OUT) :: ctrlon        ! Center longitude of the model domain (deg. E)

      INTEGER,          INTENT(OUT) :: stypin,vtypin,laiin,roufin,vegin,ndviin

      INTEGER,          POINTER :: soiltyp(:,:,:)   ! Soil type in model domain
      REAL,             POINTER :: stypfrct(:,:,:)  ! Fraction of soil types
      INTEGER,          POINTER :: vegtyp (:,:)     ! Vegetation type in model domain

      REAL,             POINTER :: lai    (:,:)     ! Leaf Area Index in model domain
      REAL,             POINTER :: roufns (:,:)     ! NDVI in model domain
      REAL,             POINTER :: veg    (:,:)     ! Vegetation fraction
      REAL,             POINTER :: ndvi   (:,:)     ! NDVI

   END SUBROUTINE readcvtsfc

   SUBROUTINE readcvtsoil(soilfile,soilfmt,nx,ny,nzsoil,nstyps,dx,dy,    &
           zpsoil,mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,   &
           zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,                     &
           tsoil,qsoil,wetcanp,snowdpth,soiltyp)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: soilfile ! Name of the soil file
    INTEGER,          INTENT(IN) :: soilfmt

    INTEGER, INTENT(OUT) :: nx            ! Number of grid points in the x-direction
    INTEGER, INTENT(OUT) :: ny            ! Number of grid points in the y-direction
    INTEGER, INTENT(OUT) :: nzsoil        ! Number of grid points in the soil.
    INTEGER, INTENT(OUT) :: nstyps        ! Number of soil types for each grid point

    REAL,    INTENT(OUT) :: dx
    REAL,    INTENT(OUT) :: dy
    INTEGER, INTENT(OUT) :: mapproj       ! Map projection
    REAL,    INTENT(OUT) :: trulat1       ! 1st real true latitude of map projection
    REAL,    INTENT(OUT) :: trulat2       ! 2nd real true latitude of map projection
    REAL,    INTENT(OUT) :: trulon        ! Real true longitude of map projection
    REAL,    INTENT(OUT) :: sclfct        ! Map scale factor
    REAL,    INTENT(OUT) :: ctrlat        ! Center latitude of the model domain (deg. N)
    REAL,    INTENT(OUT) :: ctrlon        ! Center longitude of the model domain (deg. E)

    INTEGER, INTENT(OUT) :: zpsoilin
    INTEGER, INTENT(OUT) :: tsoilin
    INTEGER, INTENT(OUT) :: qsoilin
    INTEGER, INTENT(OUT) :: wcanpin
    INTEGER, INTENT(OUT) :: snowdin

    REAL,    POINTER     :: zpsoil (:,:,:)        ! Soil depths (m)
    REAL,    POINTER     :: tsoil  (:,:,:,:)      ! Soil temperature (K)
    REAL,    POINTER     :: qsoil  (:,:,:,:)      ! Soil moisture (m3/m3)
    REAL,    POINTER     :: wetcanp(:,:,:)        ! Canopy water amount
    INTEGER, POINTER     :: soiltyp(:,:,:)        ! Soil type in model domain
    REAL,    POINTER     :: snowdpth(:,:)         ! Snow depth (m)

   END SUBROUTINE readcvtsoil

  END INTERFACE
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

  CHARACTER (LEN=256) :: houtfn
  CHARACTER (LEN=256) :: basoutfn
  INTEGER :: lbasoutf, loutf

  INTEGER :: grdbas
  INTEGER :: i,j,k,n,ireturn

  REAL :: alatnot(2)
  REAL :: ctrx, ctry, swx, swy

  REAL :: time
  INTEGER :: gboutcnt, vroutcnt
  DATA    gboutcnt, vroutcnt /0,0/

  INTEGER :: nfile
  CHARACTER(LEN=80) :: outrunname

  INTEGER :: grdbasopt           ! 0 - no  grdbas file
                                 ! 1 - one grdbas file, the first file only

  NAMELIST /output/ dirname,readyfl,outrunname,hdmpfmt,hdfcompr,        &
            grbpkbit, grdbasopt, exbcdmp,exbchdfcompr,ngbrz,zbrdmp,     &
            grdout,basout,varout,mstout,iceout,tkeout,                  &
            trbout,rainout,sfcout,landout,prcout,radout,flxout,         &
            qcexout,qrexout,qiexout,qsexout,qhexout,                    &
            filcmprs,istager,terndmp,sfcdmp,soildmp

  INTEGER :: terninfmt, sfcinfmt, soilinfmt

  NAMELIST /other_data/ terninfmt,terndta,sfcinfmt,sfcdtfl,             &
                        soilinfmt,soilinfl

  INTEGER            :: stypout,vtypout,laiout,roufout,vegout,ndviout
  CHARACTER(LEN=256) :: sfcfn,ternfn,soiloutfl
  CHARACTER(LEN=3)   :: fmtstr
  LOGICAL            :: fexist

  INTEGER :: length, ierr

  INTEGER :: zpsoilout,tsoilout,qsoilout,wcanpout,snowdout

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
     '###                Welcome to ARPSCVT                       ###', &
     '###      This program converts the history dump data        ###', &
     '###      sets generated by ARPS, between various formats.   ###', &
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
! Read in namelist
!
!-----------------------------------------------------------------------

  terninfmt = 0
  sfcinfmt  = 0
  soilinfmt = 0

  READ(5,other_data,END=100)

  WRITE(6,'(1x,/a/)') 'Namelist Other_data read in successfully.'
!
!-----------------------------------------------------------------------
!
!  Set the control parameters for output:
!
!-----------------------------------------------------------------------
!

  readyfl = 0
  grdbasopt = 1

  READ (5,output,END=100)

  WRITE(6,'(1x,/a/)') 'Namelist output read in successfully.'

  houtfmt = hdmpfmt

  IF( houtfmt == 0 ) THEN
    WRITE(6,'(/1x,a/)') 'Output format is 0, no data is dumped.'
    STOP
!  ELSE IF( houtfmt.eq.hinfmt ) THEN
!    write(6,'(/2(1x,a/))')
!    : 'Specified output data format was the same as the input,',
!    : 'no conversion will be done. Job stopped.'
!    STOP
  END IF

  IF ( houtfmt == 10 .AND. grbpkbit <= 0 ) THEN
    WRITE(6,'(a,a,i2/a)')                                               &
        'The bit width for packing GRIB data was invalid, ',            &
        'The old value was ', grbpkbit, 'Reset it to 16 bits'
    grbpkbit = 16
  END IF

  totout = 1

  GO TO 10

  100   WRITE(6,'(a)')                                                  &
          'Error reading NAMELIST file. Program ARPSCVT stopped.'
  STOP

  10    CONTINUE

  IF(hinfmt == 0) GO TO 555  ! Do not convert history files
                             ! go ahead for terrain data etc.

!-----------------------------------------------------------------------
!
! Get dimension parameters and allocate arrays
!
!-----------------------------------------------------------------------

  WRITE(6,'(2a)') '  Reading dimensions from input file - ',TRIM(grdbasfn)

  lenfil=len_trim(hisfile(1))

  CALL get_dims_from_data(hinfmt,hisfile(1)(1:lenfil),                    &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

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
  ALLOCATE(zpsoil (nx,ny,nzsoil))

  ALLOCATE(uprt   (nx,ny,nz))
  ALLOCATE(vprt   (nx,ny,nz))
  ALLOCATE(wprt   (nx,ny,nz))
  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))

  ALLOCATE(qscalar(nx,ny,nz,nscalar))

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

!    IF(grdbasopt == 2 .AND. nfile > 1) THEN
!      IF(nfile == 2) THEN
!        WRITE(grdbasfn,'(a,a,I2.2)') TRIM(grdbasfn),'.',nfile-1
!        lengbf=len_trim(grdbasfn)
!      ELSE
!        ireturn = INDEX(grdbasfn,'.',.TRUE.)
!        WRITE(grdbasfn,'(a,I2.2)') grdbasfn(1:ireturn),nfile-1
!        lengbf=len_trim(grdbasfn)
!      END IF
!
!      WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)
!
!      CALL setgbrd(0)
!    END IF
!
!-----------------------------------------------------------------------
!
!  Set the parameters that define the subdomain and skip points to
!  be used by svidump.
!
!-----------------------------------------------------------------------
!
!  CALL sdmpskip(2,2,1)
!  CALL sdmpdomn(1,nx-1,1,ny-1,1,nz-1)
!  CALL sdmpdomn(12,55,3,47,1,nz-1)
!
!-----------------------------------------------------------------------
!
!  Set the parameters that define the subdomain and skip points to
!  be used by dn2dump.
!
!-----------------------------------------------------------------------
!
!  CALL bdmpskip(2,2,1)
!  CALL bdmpdomn(1,nx-1,1,ny-1,1,nz-1)
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!

    102   CONTINUE

    DO n = 1,50
      cmnt(n) = ' '
    END DO

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

    curtim = time
    IF (nfile == nhisfile) tstop = time

    IF ( basin == 0  ) basout  = 0
    IF ( grdin == 0  ) grdout  = 0
    IF ( varin == 0  ) varout  = 0
    IF ( mstin == 0  ) mstout  = 0
    IF ( icein == 0  ) iceout  = 0
    IF ( tkein == 0  ) tkeout  = 0
    IF ( trbin == 0  ) trbout  = 0
    IF ( rainin == 0 ) rainout = 0
    IF ( sfcin == 0  ) sfcout  = 0
    IF ( landin == 0 ) landout = 0
    IF ( prcin == 0  ) prcout  = 0
    IF ( radin == 0  ) radout  = 0
    IF ( flxin == 0  ) flxout  = 0

    IF ( sfcout == 1) THEN
      snowout = 1
    END IF

    IF ( exbcdmp == 0 .OR. mstout == 0 ) THEN
      qcexout = 0
      qrexout = 0
      qiexout = 0
      qsexout = 0
      qhexout = 0
    ELSE IF ( iceout == 0 ) THEN
      qiexout = 0
      qsexout = 0
      qhexout = 0
    END IF

    IF ( exbchdfcompr > 4 ) rayklow = -1

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

    WRITE(6,'(/a,a)') " Runname is : ",runname

    CALL gtlfnkey( runname, lfnkey )
!
!-----------------------------------------------------------------------
!
!  Calculate total fields from that for base state and perturbations
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          u(i,j,k) = ubar(i,j,k) + uprt(i,j,k)
          v(i,j,k) = vbar(i,j,k) + vprt(i,j,k)
          w(i,j,k) = wbar(i,j,k) + wprt(i,j,k)
          qv(i,j,k) = qvbar(i,j,k) + qvprt(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  If grid/base state data has been written out once, skip
!  the following writing block. Also no need to write out
!  separate data for Savi3D dump.
!
!  For GrADS dump, too.
!
!-----------------------------------------------------------------------
!
    IF(houtfmt == 5.OR.houtfmt == 9.OR.houtfmt == 11) GO TO 300

    IF (grdbasopt == 0 .OR. (grdbasopt == 1 .AND. gboutcnt > 0) ) GOTO 300
                                         ! If done already, skip this part.

    CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,             &
                 1,0,basoutfn,lbasoutf)

    PRINT*
    PRINT*,'Output grid/base state file is ', basoutfn(1:lbasoutf)

    grdbas = 1      ! Dump out grd and base state arrays only

!
!-----------------------------------------------------------------------
!
!  Do the data dump.
!
!-----------------------------------------------------------------------
!
    CALL dtadump(nx,ny,nz,nzsoil, nstyps,                               &
                 houtfmt,nchout,basoutfn(1:lbasoutf),grdbas,filcmprs,   &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,               &
                 ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil, qsoil, wetcanp,snowdpth,                        &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx, radswnet,radlwin,                  &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2,tem3)

    gboutcnt = gboutcnt + 1

    300   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Then the time dependent fields:
!
!-----------------------------------------------------------------------
!
    IF( houtfmt == 5 .AND. vroutcnt == 1 ) GO TO 130
    IF( houtfmt == 9 .AND. vroutcnt == 1 ) GO TO 130
                             ! use the same file name
!
!-----------------------------------------------------------------------
!
!  Reconstruct the file name using the specified directory name
!
!-----------------------------------------------------------------------
!
    CALL gtdmpfn(runname(1:lfnkey),dirname,                             &
                 ldirnam,curtim,hdmpfmt,                                &
                 1,0, houtfn, loutf)

    CALL fnversn(houtfn, loutf )

    130   CONTINUE

    PRINT*
    WRITE(6,'(1x,a,f8.1,a,a)')                                          &
          'Output file at time ',curtim,' (s) is ', houtfn(1:loutf)

    grdbas = 0      ! Not just dump out grd and base state arrays only

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                &
                 houtfmt,nchout,houtfn(1:loutf),grdbas,filcmprs,        &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,               &
                 ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2,tem3)

    vroutcnt = 1 ! Variables have been written out at least once

    IF ( hinfmt == 9 .AND. ireturn /= 2 ) GO TO 102

  END DO

  IF (houtfmt == 5) THEN

    CALL mclosescheme (gridid, ierr)
    CALL mclosedataset (dsindex, ierr)

  END IF

  IF ( houtfmt == 9 ) THEN
    CLOSE (UNIT=nchout)
  END IF

!-----------------------------------------------------------------------
!
! Processing terrain data
!
!-----------------------------------------------------------------------

  555 CONTINUE

  IF(ldirnam > 0) THEN
    IF(dirname(ldirnam:ldirnam) /= '/') THEN
       dirname(ldirnam+1:ldirnam+1) = '/'
       ldirnam = ldirnam + 1
    END IF
  ELSE
     dirname = './'
  END IF

  IF(terninfmt /= 0) THEN

    INQUIRE(FILE=TRIM(terndta), EXIST = fexist )
    IF(.NOT. fexist) THEN
      WRITE(6,*) 'The terrain file ',TRIM(terndta),                     &
                 ' you specified does not exist. Terrain data skipped.'
      GO TO 666
    END IF

    CALL readcvttrn(TRIM(terndta),terninfmt,nx,ny,dx,dy,                &
                    mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,&
                    htptr)

    IF(terndmp == 1) THEN
      WRITE(fmtstr,'(a3)') 'bin'
    ELSE IF(terndmp == 3) THEN
      WRITE(fmtstr,'(a3)') 'hdf'
    ELSE IF(terndmp == 7) THEN
      WRITE(fmtstr,'(a3)') 'net'
    ELSE
      WRITE(6,*) 'Unsupported terrain dump format. Terrain data skipped.'
      GO TO 666
    END IF
    i = INDEX( terndta,'/',.TRUE.)+1
    j = INDEX( terndta,'.',.TRUE.)
    WRITE(ternfn,'(4a)') TRIM(dirname),terndta(i:j),TRIM(fmtstr),'trndata'

    CALL writtrn(nx,ny,TRIM(ternfn),dx,dy,                              &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   htptr)
  END IF

!-----------------------------------------------------------------------
!
! Processing surface characteristics data
!
!-----------------------------------------------------------------------

  666 CONTINUE

  IF(sfcinfmt /= 0) THEN
    INQUIRE(FILE=TRIM(sfcdtfl), EXIST = fexist )
    IF(.NOT. fexist) THEN
      WRITE(6,*) 'The surface file ',TRIM(sfcdtfl),                     &
                 ' you specified does not exist. Surface data skipped.'
      GO TO 777
    END IF

    CALL readcvtsfc(sfcdtfl,sfcinfmt,nx,ny,nstyps,dx,dy,                &
                    mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,&
                    stypout,vtypout,laiout,roufout,vegout,ndviout,      &
                    stypptr,sfrctptr,vegtypptr,laiptr,roufnsptr,vegptr, &
                    ndviptr)

    IF(sfcdmp == 1) THEN
      WRITE(fmtstr,'(a3)') 'bin'
    ELSE IF(sfcdmp == 3) THEN
      WRITE(fmtstr,'(a3)') 'hdf'
    ELSE IF(sfcdmp == 7) THEN
      WRITE(fmtstr,'(a3)') 'net'
    ELSE
      WRITE(6,*) 'Unsupported surface data dump format. Surface data skipped.'
      GO TO 777
    END IF
    i = INDEX( sfcdtfl,'/',.TRUE.)+1
    j = INDEX( sfcdtfl,'.',.TRUE.)
    WRITE(sfcfn,'(4a)') TRIM(dirname),sfcdtfl(i:j),TRIM(fmtstr),'sfcdata'

    nstyp = nstyps

    CALL wrtsfcdt(nx,ny,nstyps,TRIM(sfcfn),dx,dy,                       &
                  mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,  &
                  stypout,vtypout,laiout,roufout,vegout,ndviout,        &
                  stypptr,sfrctptr,vegtypptr,laiptr,roufnsptr,vegptr,   &
                  ndviptr)

  END IF

!-----------------------------------------------------------------------
!
! Processing soil data file
!
!-----------------------------------------------------------------------

  777 CONTINUE

  IF(soilinfmt /= 0) THEN
    INQUIRE(FILE=TRIM(soilinfl), EXIST = fexist )
    IF(.NOT. fexist) THEN
      WRITE(6,*) 'The surface file ',TRIM(soilinfl),                    &
                 ' you specified does not exist. Soil data skipped.'
      GO TO 888
    END IF

    CALL readcvtsoil(TRIM(soilinfl),soilinfmt,nx,ny,nzsoil,nstyps,dx,dy,&
           zpsoilptr,mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,&
           zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,               &
           tsoilptr,qsoilptr,wetcanpptr,snowdpthptr,soiltypptr)

    i = INDEX( soilinfl,'/',.TRUE.)+1
    j = LEN_TRIM(soilinfl)
    WRITE(houtfn,'(a)') soilinfl(i:j)
    IF(soildmp == 1) THEN
      WRITE(fmtstr,'(a3)') 'bin'
    ELSE IF(soildmp == 3) THEN
      WRITE(fmtstr,'(a3)') 'hdf'
    ELSE IF(soildmp == 7) THEN
      WRITE(fmtstr,'(a3)') 'net'
    ELSE
      WRITE(6,*) 'Unsupported soil data dump format. Soil data skipped.'
      GO TO 888
    END IF

    IF (soilinfmt == 1) THEN
      k = INDEX(houtfn,'bin')
    ELSE IF (soilinfmt == 3) THEN
      k = INDEX(houtfn,'hdf')
    ELSE IF (soilinfmt == 7) THEN
      k = INDEX(houtfn,'net')
    END IF

    IF(k >0 .AND. k < (j-i+1)) THEN
      WRITE(houtfn(k:k+2),'(a)') fmtstr
    ELSE
      k = INDEX(houtfn,'.',.TRUE.)
      WRITE(soiloutfl,'(a)')  houtfn(k+1:LEN_TRIM(houtfn))
      WRITE(houtfn,   '(3a)') houtfn(1:k),fmtstr,TRIM(soiloutfl)
    END IF

    WRITE(soiloutfl,'(2a)') TRIM(dirname),TRIM(houtfn)

    nstyp = nstyps

    CALL wrtsoil(nx,ny,nzsoil,nstyps,TRIM(soiloutfl),dx,dy,zpsoilptr,   &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,               &
           tsoilptr,qsoilptr,wetcanpptr,snowdpthptr,soiltypptr)
  END IF

!-----------------------------------------------------------------------
!
! End of the program
!
!-----------------------------------------------------------------------

  888 CONTINUE

  WRITE(6,'(1x,/,a,/)') ' ==== Program ARPSCVT terminated normally. ===='

  STOP

  997   CONTINUE
  WRITE(6,'(1x,a,i2,/1x,a)')                                            &
      'Data read was unsuccessful. ireturn =', ireturn,                 &
      'Job stopped.'
  STOP
END PROGRAM ARPSCVT
