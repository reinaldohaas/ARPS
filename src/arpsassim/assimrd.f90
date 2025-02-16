!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ASSIMRD                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                     &
           isrc,itag,ireturn,assimtim,                                  &
           ubar,vbar,pbar,ptbar,rhostr,qvbar,rhobar,                    &
           u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,                 &
           tem1,tem5,tem6,tem7,tem8,tem9,tem10,tem11,                   &
           tem12,tem13)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinates the ingestion of history data dump files in various
!  formats and the ingestion of 'pre-processed' Lincoln Lab radar data.
!  Other data types can be used by selecting the 'other' option for
!  the 'dtyp' parameter in ASSIM.INPUT. The velocity fields can be
!  used for an insertion, adjustment and/or recovery (see ASSIMCON.F).
!
!  NOTE:  In order to use the 'other' option, corresponding changes
!         are needed in this routine. In this version, only two type
!         of data are processed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Alan Shapiro and Steven Lazarus
!
!  2/25/1993
!
!  MODIFICATION HISTORY:
!
!  02/22/96 (Limin Zhao)
!  Changes are made to match the new version of assimilation package.
!
!  02/26/96 (Steve Lazarus & Limin Zhao)
!  Added the simple moisture parameterization (Kessler formula) for
!  qr and qv.
!
!  02/27/96 (Limin Zhao)
!  A critical value of qr is used to set flags at the area outside of
!  rainwater regions for OSSE type of experiments.
!
!  02/29/96 (Limin Zhao)
!  Added two temporary arrays to avoid overwritten of pbar and ptbar,
!  which will be used in thermodynamic retrieval.
!
!  05/08/96 (Limin Zhao)
!  A bug was fixed when control input parameters are set
!  varopt=0,insrtopt=0 and recovopt =1.
!  A correct filename is given now.
!
!  05/14/96 (Limin Zhao and Alan Shapiro)
!  Modified the stratage to get qv by using background temperature
!  and pressure information.
!
!  05/15/96 (Limin Zhao)
!  A control parameter is hardwired to control the moisture retrieval
!  switch on or off.
!
!  09/24/96 (Limin Zhao)
!  Fix a bug in qr reading with isrc=4. qr should be read in for
!  thermodynamic retrieval.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels.
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of soil levels (m)
!    isrc     Flag, indicating the calling routine
!    itag     Counter indicating file name and time
!
!  OUTPUT:
!
!    ireturn  Flag, indicating read status of input file
!             = 0 Successful read
!             = 1 Unsuccessful read
!
!    assimtim Times of all input files
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state density (kg/m**3) times j3
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rain water mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  WORK ARRAYS:
!
!    tem1      Temporary work array.
!    tem5      Temporary work array.
!    tem6      Temporary work array.
!    tem7      Temporary work array.
!    tem8      Temporary work array.
!    tem9      Temporary work array.
!    tem10     Temporary work array.
!    tem11     Temporary work array.
!    tem12     Temporary work array.
!    tem13     Temporary work array.
!
!    rhobar   Base state density (kg/m**3)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nt           ! Number of time levels of time-dependent arrays.
  INTEGER :: tim          ! Index of time level.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  INTEGER :: nzsoil

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  INTEGER :: isrc              ! Flag indicating source of calling routine
  INTEGER :: itag              ! Counter indicating file name and time
  INTEGER :: nstyps

  REAL :: assimtim(100)        ! Time of data input files

  REAL :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity (kg/kg)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz,nt)  ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz,nt)  ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz,nt)  ! Hail mixing ratio (kg/kg)

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array
  REAL :: tem13 (nx,ny,nz)     ! Temporary work array

  INTEGER :: ireturn

  REAL :: frsvnth, coef
  REAL :: rad,xs,ys,zs
  REAL :: xmove,ymove
  REAL :: plcl

  INTEGER :: in,jn
  INTEGER :: klcl
  INTEGER :: itimesv

  REAL, ALLOCATABLE :: tem4dsoil(:,:,:,:)
  REAL, ALLOCATABLE :: tem3dsoil(:,:,:)

!
!   real dummy(153,99,43)      !Note: temporary fix for extra arrays
!   integer idummy(153,99,43)
!
!   real dummy(133,133,37)      !Note: temporary fix for extra arrays
!   integer idummy(133,133,37)
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL dtahead
  EXTERNAL dtaread
  EXTERNAL radread
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: is      ! Counter, used for reading input data
  SAVE    is
  DATA    is/1/

  INTEGER :: i,j,k,n ! Loop index
  INTEGER :: nchanl  ! FORTRAN I/O channel number for history data output.

  INTEGER :: lengbf,lendtf,lbasdmpf

  REAL :: xor     ! Coordinate (xor,yor,zor) of the model grid
  REAL :: yor     ! origin relative to the radar.
  REAL :: zor     !

  REAL :: storstop         ! Temporarily stores the model stop time

  CHARACTER (LEN=80  ) :: saverunm  ! Temporarily stores the name of this run
  CHARACTER (LEN=256 ) :: filein  ! Input file name for recovery

  REAL :: temp,pres,qvsat

  INTEGER :: imoist
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'assim.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
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
!  If isrc=1, then read the times (assimtim) from 3 data file headers.
!  This is done in order to determine when either a recovery,
!  adjustment, and/or an insertion should be performed. It is
!  assumed that 3 data files exist at the beginning of the
!  assimilation period. For this reason the counter, itag, is
!  initially set to 3 as well as updated in ASSIMCON.F.
!
!  The grid/base state file (inigbf) and input file names (assimdat)
!  are to be provided by the user in ASSIM.INPUT. The grid/base
!  state file is used when model data is ingested only.
!
!  NOTE:  None of the time-dependent variables are read in during
!         this step.
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) 'code in assimrd: is,isrc=',is,isrc

  IF( isrc == 1) THEN    ! Called by ASSIMCON

    DO i=1,itag

      tim = i

      lendtf=LEN(assimdat(i))
      CALL strlnth( assimdat(i), lendtf)
      WRITE(6,*) 'The filename assimdat(',i,')=',                       &
                         assimdat(i)(1:lendtf)
!c
!clz
!c     The option for dtyp=0 is not well considered for SOP, It will
!c     be modified in future
!c
      IF( dtyp == 0) THEN

        saverunm = runname

        CALL dtahead(nx,ny,nz,                                          &
             inifmt,inigbf,lengbf,assimdat(i),                          &
             lendtf,assimtim(i),x,y,z,zp,tem11,tem12,tem13,             &
             ptprt(1,1,1,tim),pprt(1,1,1,tim),tem7,                     &
             qc(1,1,1,tim),qr(1,1,1,tim),qi(1,1,1,tim),                 &
             qs(1,1,1,tim),qh(1,1,1,tim),tem9,                          &
             ubar,vbar,tem8,ptbar,pbar,rhobar,qvbar,                    &
             ireturn,tem1,tem9,tem10)

        runname = saverunm

        IF(ireturn /= 0) RETURN

      ELSE IF(dtyp == 1) THEN

        CALL radread(assimdat(i),lendtf,                                &
                       isrc,assimtim(i),ireturn,                        &
                       nx,ny,nz,dx,dy,dz,                               &
                       xor,yor,zor,                                     &
                       tem1,tem5,tem6,tem7,tem8,                        &
                       tem9,tem10,tem11,tem12,tem13)

        IF(ireturn /= 0) RETURN

      END IF

    END DO

!
!-----------------------------------------------------------------------
!
!  Read in a single header to recover the time. For recovopt=1
!
!  NOTE:  None of the time-dependent variables are read in during
!         this step.
!
!-----------------------------------------------------------------------
!

  ELSE IF(isrc == 2) THEN    ! Called by ASSIMCON

    tim = tpresent

    lendtf=LEN(assimdat(itag))
    CALL strlnth( assimdat(itag), lendtf)
    WRITE(6,*) 'The filename assimdat(',itag,')=',                      &
                assimdat(itag)(1:lendtf)

    IF( dtyp == 0) THEN

      saverunm = runname

      CALL dtahead(nx,ny,nz,                                            &
           inifmt,inigbf,lengbf,assimdat(itag),                         &
           lendtf,assimtim(itag),x,y,z,zp,tem11,tem12,tem13,            &
           ptprt(1,1,1,tim),pprt(1,1,1,tim),tem7,                       &
           qc(1,1,1,tim),qr(1,1,1,tim),qi(1,1,1,tim),                   &
           qs(1,1,1,tim),qh(1,1,1,tim),tem9,                            &
           ubar,vbar,tem8,ptbar,pbar,rhobar,qvbar,                      &
           ireturn,tem1,tem9,tem10)

      runname = saverunm

      IF( ireturn /= 0) RETURN

    ELSE IF(dtyp == 1) THEN

      CALL radread(assimdat(itag),lendtf,                               &
                     isrc,assimtim(itag),ireturn,                       &
                     nx,ny,nz,dx,dy,dz,                                 &
                     xor,yor,zor,                                       &
                     tem1,tem5,tem6,tem7,tem8,                          &
                     tem9,tem10,tem11,tem12,tem13)

      IF(ireturn /= 0) RETURN

    END IF

!
!-----------------------------------------------------------------------
!
!  Read in a single data file. The radial velocities are used
!  to perform a variational adjustment and/or insertion.
!
!  NOTE: If model data is ingested (dtyp=0) all current model
!       prognostic variables are overwritten with the exeption of
!       u,v,w and qv. uprt,vprt,wprt, and qvprt are returned in
!       tem4,tem5,tem6,and tem7 respectively.  The total fields
!       can then reconstructed from these and their respective base
!       state quantities. qv is reconstructed here.
!
!-----------------------------------------------------------------------
!

  ELSE IF(isrc == 3) THEN    ! Called by ASSIMVEL

    tim = tpresent

    IF( dtyp == 0) THEN

      lengbf=LEN(inigbf)
      CALL strlnth( inigbf, lengbf)
      WRITE(6,'(/a,a)')                                                 &
          'The grid/base state file name is ',inigbf(1:lengbf)

      lendtf=LEN(assimdat(is))
      CALL strlnth( assimdat(is), lendtf)
      WRITE(6,'(/a,a)')'The file name is ',assimdat(is)(1:lendtf)

      saverunm = runname

!
!-----------------------------------------------------------------------
!
!  Set the tem8, tem5 and tem6  arrays to zero.  This is done to avoid
!  overwriting wbar, pbar and ptbar on calls to dtaread.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem8(i,j,k) = 0.0
            tem5(i,j,k) = 0.0
            tem6(i,j,k) = 0.0
          END DO
        END DO
      END DO

      storstop = tstop  ! Temporary. For ingestion of model data
                        ! generated prior to 4.0
      CALL ctim2abss( year, month, day, hour, minute, second, itimesv )

      ALLOCATE(tem4dsoil(nx,ny,nzsoil,0:nstyps))
      ALLOCATE(tem3dsoil(nx,ny,nstyps))

      CALL dtaread(nx,ny,nz,nzsoil,nstyps,                              &
           inifmt,nchanl,inigbf,lengbf,assimdat(is),                    &
           lendtf,assimtim(is),x,y,z,zp,zpsoil,tem11,tem12,tem13,       &
           ptprt(1,1,1,tim),pprt(1,1,1,tim),tem7,                       &
           qc(1,1,1,tim),qr(1,1,1,tim),qi(1,1,1,tim),                   &
           qs(1,1,1,tim),qh(1,1,1,tim),tem9,tem9,tem9,                  &
           ubar,vbar,tem8,ptbar,pbar,rhobar,qvbar,                      &
           tem3dsoil,tem3dsoil,tem1(1,1,1),                             &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                         &
           tem4dsoil,tem4dsoil,tem1(1,1,1),tem1(1,1,1),                 &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                         &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                         &
           tem1(1,1,1),tem1(1,1,1),                                     &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),             &
           ireturn,tem6,tem9,tem10)
      CALL abss2ctim( itimesv, year, month, day, hour, minute, second )

      tstop   = storstop
      runname = saverunm

      IF(ireturn /= 0) RETURN
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            rhostr(i,j,k) = rhobar(i,j,k)*j3(i,j,k)
            qv(i,j,k,tim) = tem7(i,j,k) + qvbar(i,j,k)
          END DO
        END DO
      END DO

      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            tem11(i,j,k)  = tem11(i,j,k) + ubar(i,j,k)
            tem12(i,j,k)  = tem12(i,j,k) + vbar(i,j,k)
            tem13(i,j,k)  = tem13(i,j,k)
          END DO
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!
!  Interpolate input 'model' velocity components to scalar grid point.
!  Where:
!              tem7 = scalar u component
!              tem8 = scalar v component
!              tem9 = scalar w component
!
!-----------------------------------------------------------------------
!
      CALL avgx(tem11, 0, nx,ny,nz,                                     &
                1,nx-1,1,ny-1,1,nz-1,tem7)
      CALL avgy(tem12, 0, nx,ny,nz,                                     &
                1,nx-1,1,ny-1,1,nz-1,tem8)
      CALL avgz(tem13, 0, nx,ny,nz,                                     &
                1,nx-1,1,ny-1,1,nz-1,tem9)
!
!-----------------------------------------------------------------------
!
!   Calculate the radial velocity component of the input model data.
!   This is for OSSE type experiments. The input model data may or
!   may not match that of the current model stream.
!
!   WARNING: The insertion below assumes that, the user has supplied
!        the radar coordinates with respect to the lower left hand
!        corner of the grid (See ASSIM.INPUT). The grid origin is
!        assumed to be at (0,0,0).
!
!                          * (0,0,0) grid origin
!                       .  .
!                    .     . ymove
!                 .        .
!      radar   (- ..........
!                   xmove
!
!
!-----------------------------------------------------------------------
!
      umove=0.0
      vmove=0.0
      xmove= xshift - umove*(curtim-assimtim(1))
      ymove= yshift - vmove*(curtim-assimtim(1))

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1

            xs    = 0.5*(x(i)+x(i+1)) - xmove
            ys    = 0.5*(y(j)+y(j+1)) - ymove
            zs    = 0.5*(z(k)+z(k+1)) - zshift

            rad   = SQRT(xs**2+ys**2+zs**2)

!
!
!       tem1(i,j,k) = (tem7(i,j,k)*xs + tem9(i,j,k)*ys
!  :                +  tem9(i,j,k)*zs)/rad

            tem1(i,j,k) = (tem7(i,j,k)*xs + tem8(i,j,k)*ys              &
                        +  tem9(i,j,k)*zs)/rad


          END DO
        END DO
      END DO
!
      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            tem11(i,j,k)  = tem11(i,j,k) - ubar(i,j,k)
            tem12(i,j,k)  = tem12(i,j,k) - vbar(i,j,k)
            tem13(i,j,k)  = tem13(i,j,k)
          END DO
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!  Set flags outside of rainwater area.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            IF(qr(i,j,k,tim) < 0.0005) THEN
              tem1(i,j,k)  = spval
              tem11(i,j,k) = spval
              tem12(i,j,k) = spval
              tem13(i,j,k) = spval
            END IF
          END DO
        END DO
      END DO
!
!clz    for real data ingestion
!
    ELSE IF(dtyp == 1) THEN
!cc
!cc
!ccNote: no overwritten for model variables.
!cc
!cc
      lendtf=LEN(assimdat(is))
      CALL strlnth( assimdat(is), lendtf)
      WRITE(6,'(/a,a)')'The filename', assimdat(is)(1:lendtf)

      CALL radread(assimdat(is),lendtf,                                 &
                   isrc,assimtim(is),ireturn,                           &
                   nx,ny,nz,dx,dy,dz,                                   &
                   xor,yor,zor,                                         &
                   tem1,tem5,tem6,tem7,tem8,                            &
                   tem9,tem10,tem11,tem12,tem13)


      IF(ireturn /= 0) RETURN
!
!-----------------------------------------------------------------------
!
!  Compute the rainwater from the radar reflectivity factor. See
!  Kessler (1969).  Reflectivity is assumed to be in dBz, i.e.,
!  Z (mm**6/m**3), and dBz = 10log10(Z).
!
!  Note: When reflectivity is too small or too large, qr might be
!        reliable from the qr-Z relation, a upper limit and lower
!        limit are set up here.
!
!
!     Set qr max at 50 dbz for AMS99 testing (SSW 12-7-98)
!
!
!-----------------------------------------------------------------------
!
      imoist=1
      IF(imoist == 0) GO TO 911

      WRITE(6,*) 'This run is with qr and qv'

      coef    = LOG(10.)/10.
      frsvnth = 4./7.

      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            IF (tem8(i,j,k) > 50.0.AND.tem8(i,j,k) /= spval) THEN
!           IF (tem8(i,j,k).gt.75.0.and.tem8(i,j,k).ne.spval) THEN
!           IF (tem8(i,j,k).gt.75.0.and.tem8(i,j,k).ne.spval) THEN
!              tem8(i,j,k) = 75.0
              tem8(i,j,k) = 50.0
            END IF
!          IF (tem8(i,j,k).ne.spval.and.tem8(i,j,k).gt.20.0) THEN
!          IF (tem8(i,j,k).ne.spval.and.tem8(i,j,k).gt.0.0) THEN
!          IF (tem8(i,j,k).ne.spval.and.tem8(i,j,k).gt.5.0) THEN
            IF (tem8(i,j,k) /= spval) THEN
              qr(i,j,k,tim) = (.001/rhobar(i,j,k))*                     &
                     ( (1./17300.)*EXP(coef*tem8(i,j,k)) )**frsvnth
            ELSE
              qr(i,j,k,tim) = 0.0
            END IF
          END DO
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!  For real data, dtyp=1, assume qv=qvs in rainwater regions with
!  updraft. In downdraft area, model background mean is used.
!
!  Calculate the saturation water vapor specific humidity,
!  qvs, from pressure and potential temperature using Teten's formula.
!
!
!     Comment out Limin's old qvsat for AMS99 testing (SSW 12-7-98)
!
!
!-----------------------------------------------------------------------
!
!         CALL satmrpt(nx,ny,nz,pbar,ptbar,tem9)    ! calculate qvs

!    DO 195 k=1,nz-1
!    DO 195 j=1,ny-1
!    DO 195 i=1,nx-1


!      pres = pbar(i,j,k)+pprt(i,j,k,tim)
!      temp = (ptbar(i,j,k)+ptprt(i,j,k,tim))*((pres/p0)**rddcp)
!ccxxx          IF(temp.ge.273.16) THEN
!      qvsat=(380. / pres) *               !Teten's formula
!    :             exp( 17.27 * (temp-273.15)/(temp-35.86) )
!ccxxx          ELSE
!ccxxx          qvsat=(380. / pres) *
!ccxxx     :             exp( 21.875* (temp-273.15)/(temp-7.5) )
!ccxxx          ENDIF
!ccxxx
!      IF(tem8(i,j,k).gt.20.0.and.tem8(i,j,k).ne.spval
!    :         .and.w(i,j,k,tim).ge.0.0) THEN
!         qv(i,j,k,tim) = qvsat
!      ELSE IF(tem8(i,j,k).gt.20.0.and.tem8(i,j,k).
!    :              ne.spval.and.w(i,j,k,tim).lt.0.0) THEN
!         qv(i,j,k,tim) = 0.8*qvsat
!      ENDIF

!195      CONTINUE

      911     CONTINUE

    END IF

    is = is + 1
!
!-----------------------------------------------------------------------
!
!  Read in the winds, which may (or may not) be variationally adjusted,
!  at three time steps, i.e.
!
!                    t1         t2         t3
!                  ---|----------|----------|---
!
!  The time interval between these files corresponds to the frequency
!  of the radar data. The time between each radar observation does not
!  have to be the same (see RETRINT.F).
!  NOTE:
!  The adjusted winds are dumped via the dump3d routine (see ASSIMVEL.F).
!  Hence, for the recovery option (recovopt=1), the adjusted winds and
!  associated model data at that time step are always written off in
!
!  a model format as specified by the user in ARPSINIT.INPUT
!
!-----------------------------------------------------------------------
!

  ELSE IF(isrc == 4) THEN    ! Called by RETPTPR

    tim = 1

    lengbf=LEN(inigbf)
    CALL strlnth( inigbf, lengbf)
    WRITE(6,'(/a,a)')                                                   &
        'The grid/base state file name is ',inigbf(1:lengbf)
!
!
!ccxxx
!ccxxx    You have to hardwire the filenames if you choose to do
!ccxxx    thermodynamic retrieval only
!ccxxx
!ccxxx        adjdat(1) = 'assim4.bin000032'
!ccxxx        adjdat(2) = 'assim4.bin000384'
!ccxxx        adjdat(3) = 'assim4.bin000736'
!
!     adjdat(1) = 'ADA_15.bin000000.9408171830'
!     adjdat(2) = 'ADA_15.bin000000.9408171836'
!     adjdat(3) = 'ADA_15.bin000000.9408171841'
!
!     adjdat(1) = 'WANG15.bin000000.940817183548'
!     adjdat(2) = 'WANG15.bin000000.940817183600'
!     adjdat(3) = 'WANG15.bin000000.940817183612'
!
!    adjdat(1) = 'KI1881.bin000000'
!    adjdat(2) = 'KI1881.bin000351'
!    adjdat(3) = 'KI1881.bin000702'

    DO n = itag-3,itag-1

      WRITE(6,*)'itag,n',itag,n

!      IF (varopt.eq.0.and.insrtopt.eq.0) THEN
!         lendtf=len(assimdat(n))
!         CALL strlnth( assimdat(n), lendtf)
!         filein = assimdat(n)(1:lendtf)
!         write(6,'(/a,a)')
!  :      'The input file name is ',filein
!       ELSE

      lendtf=LEN(adjdat(n))
      CALL strlnth( adjdat(n), lendtf)
      filein = adjdat(n)(1:lendtf)
      WRITE(6,'(/a,a)')                                                 &
          'The input file name is ',filein

!       ENDIF

      saverunm = runname
!
!-----------------------------------------------------------------------
!
!  Set the tem8 array to zero.  This is done to avoid overwriting
!  wbar on calls to dtaread.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem8(i,j,k) = 0.0
            tem7(i,j,k) = 0.0
            tem6(i,j,k) = 0.0
          END DO
        END DO
      END DO

      storstop = tstop  ! Temporary. For ingestion of model data
                        ! generated prior to 4.0

      CALL ctim2abss( year, month, day, hour, minute, second, itimesv )

      ALLOCATE(tem4dsoil(nx,ny,nzsoil,0:nstyps))
      ALLOCATE(tem3dsoil(nx,ny,nstyps))

      CALL dtaread(nx,ny,nz,nzsoil,nstyps,                              &
           hdmpfmt,nchanl,inigbf,lengbf,filein,                         &
           lendtf,assimtim(n),x,y,z,zp,zpsoil,tem11,tem12,tem13,        &
           ptprt(1,1,1,tim),pprt(1,1,1,tim),tem7,                       &
           qc(1,1,1,tim),tem6,qi(1,1,1,tim),                            &
           qs(1,1,1,tim),qh(1,1,1,tim),tem1,tem1,tem1,                  &
           ubar,vbar,tem8,ptbar,pbar,rhobar,qvbar,                      &
           tem3dsoil,tem3dsoil,tem1(1,1,1),                            &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                         &
           tem4dsoil,tem4dsoil,tem1(1,1,1),tem1(1,1,1),                 &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                         &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                         &
           tem1(1,1,1),tem1(1,1,1),                                     &
           tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),             &
           ireturn,tem5,tem9,tem10)

      DEALLOCATE(tem4dsoil)
      DEALLOCATE(tem3dsoil)

      CALL abss2ctim( itimesv, year, month, day, hour, minute, second )
!
!
!       if( n.eq.1 ) assimtim(n) = 0.0  ! 18:30z
!       if( n.eq.2 ) assimtim(n) = 360.0 ! 18:36z
!       if( n.eq.3 ) assimtim(n) = 660.0 ! 18:41z
!
!       if( n.eq.1 ) assimtim(n) = 0.0  ! 18:35:48z
!       if( n.eq.2 ) assimtim(n) = 12.0 ! 18:36:00z
!       if( n.eq.3 ) assimtim(n) = 24.0 ! 18:36:12z
!
!      if( n.eq.1 ) assimtim(n) =   0.0 ! 18:35:36z     for SDVR data
!      if( n.eq.2 ) assimtim(n) = 351.0 ! 18:41:26z
!      if( n.eq.3 ) assimtim(n) = 702.0 ! 18:47:17z
!

      tstop   = storstop
      runname = saverunm

      IF( ireturn /= 0) RETURN
!
!-----------------------------------------------------------------------
!
!  The arrays tem11,tem12 and tem13 contain uprt, vprt, and wprt
!  respectively.  The recovery requires the total velocities,i.e.,
!  u = uprt + ubar.
!
!  Reconstruct qv from qvprt which is returned from the call to
!  DTAREAD in the tem7 array.
!
!  Construct rhostr from rhobar and j3.
!
!-----------------------------------------------------------------------
!
      DO i = 1,nx
        DO j = 1,ny-1
          DO k = 1,nz-1
            u(i,j,k,tim) = tem11(i,j,k) + ubar(i,j,k)
          END DO
        END DO
      END DO
      DO i = 1,nx-1
        DO j = 1,ny
          DO k = 1,nz-1
            v(i,j,k,tim) = tem12(i,j,k) + vbar(i,j,k)
          END DO
        END DO
      END DO
      DO i = 1,nx-1
        DO j = 1,ny-1
          DO k = 1,nz
            w(i,j,k,tim) = tem13(i,j,k)
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qv(i,j,k,tim)     = tem7(i,j,k) + qvbar(i,j,k)
            qr(i,j,k,tim)     = tem6(i,j,k)
          END DO
        END DO
      END DO

      tim = tim + 1

    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rhostr(i,j,k) = rhobar(i,j,k)*j3(i,j,k)
        END DO
      END DO
    END DO

  END IF

  RETURN

END SUBROUTINE assimrd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE PALCL                     ######
!######                                                      ######
!######                  Steven Lazarus                      ######
!######               University of Oklahoma                 ######
!######               School of  Meteorology                 ######
!##################################################################
!##################################################################

SUBROUTINE palcl(nx,ny,nz,in,jn,p,pienv,tinit,qinit,klcl,plcl)

!
!-----------------------------------------------------------------------
!
!  Purpose:
!
!  Subroutine calculates the lifting condensation level
!
!-----------------------------------------------------------------------
!
!  Author: Steve Lazarus
!          6/10/93
!
!  Modification History:
!
!-----------------------------------------------------------------------
!
!  Input :
!
!    nx,ny,nz Model grid dimensions
!    in,jn    A given point in the horizontally homogenous domain
!    p      Sounding pressure (Pa)
!    pienv    Sounding pressure (non-dimensional)
!    tinit    Sounding temperature (K)
!    qinit    Sounding mixing ratio (kg/kg)
!
!    Cp       Specific heat for dry air (kg/(m s**2))
!    Rd       Gas constant for dry air  (kg/(m s**2))
!    p0       Surface pressure (Pa)
!
!  Output:
!
!    klcl     Index just above LCL
!    plcl     Pressure at LCL (Pa)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: pie, qvs
  REAL :: qtest
  REAL :: plcl
  REAL :: plcl1
  REAL :: pilcl
  REAL :: pilcl1
  REAL :: pavg
  REAL :: piavg

  INTEGER :: in,jn
  INTEGER :: nx,ny,nz

  REAL :: p(nx,ny,nz),pienv(nx,ny,nz),tinit(nx,ny,nz),qinit(nx,ny,nz)

  INTEGER :: k,klcl
  INTEGER :: iter
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Calculate non-dimensional pressure, pienv.
!
!-----------------------------------------------------------------------
!
  PRINT *,'Entering subroutine palcl'

  DO k=1,nz
    pienv(in,jn,k) = pie(p(in,jn,k),p0,rd,cp)
  END DO

  k=1
  qtest = 999.

  10    CONTINUE
  IF(qtest > qinit(in,jn,2)) THEN
    k= k + 1
    qtest = qvs( p(in,jn,k),pienv(in,jn,k),tinit(in,jn,2))
    GO TO 10
  END IF

  klcl   = k

  IF(klcl == 1) klcl=2

!
!-----------------------------------------------------------------------
!
!  Interpolate to determine the 'exact' pressure of the LCL
!
!-----------------------------------------------------------------------
!
  plcl   = p(in,jn,klcl)
  plcl1  = p(in,jn,klcl-1)
  pilcl  = pienv(in,jn,klcl)
  pilcl1 = pienv(in,jn,klcl-1)
  qtest  = 999.

  iter = 0
  20    CONTINUE
  IF( ABS(qtest - qinit(in,jn,2)) > .00001) THEN
    iter  = iter + 1
    pavg  = (plcl  + plcl1 )*.5
    piavg = (pilcl + pilcl1)*.5
    qtest = qvs(pavg,piavg,tinit(in,jn,2))
    IF(qtest > qinit(in,jn,1)) THEN
      plcl1  = pavg
      pilcl1 = piavg
    ELSE
      plcl   = pavg
      pilcl  = piavg
    END IF
    GO TO 20
  END IF

  plcl = pavg
  pilcl= piavg

  RETURN
END SUBROUTINE palcl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION PIE                       ######
!######                                                      ######
!######                  Steven Lazarus                      ######
!######               University of Oklahoma                 ######
!######               School of  Meteorology                 ######
!##################################################################
!##################################################################

  FUNCTION pie(pdim,p0,rd,cp)
!
!-----------------------------------------------------------------------
!
!  Purpose:
!
!  Calculate the non-dimensional pressure given the dimensional
!  pressure.
!
!-----------------------------------------------------------------------
!
!  Author: Steve Lazarus
!          6/10/93
!
!  Modification History:
!
!-----------------------------------------------------------------------
!
!  Input :
!
!    pdim     Sounding pressure (Pa)
!    p0       Surface pressure (Pa)
!    Rd       Gas constant for dry air  (kg/(m s**2))
!    Cp       Specific heat for dry air (kg/(m s**2))
!
!  Output:
!
!    pie      Non-dimensional pressure
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------
!

  REAL :: pdim,rd,cp,p0
!
!
  pie  =  (pdim/p0)**(rd/cp)
!
! end of function pie
  END FUNCTION pie
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION QVS                       ######
!######                                                      ######
!######                  Steven Lazarus                      ######
!######               University of Oklahoma                 ######
!######               School of  Meteorology                 ######
!##################################################################
!##################################################################

  FUNCTION qvs(pdim,pnon,t)
!
!-----------------------------------------------------------------------
!
!  Purpose:
!
!  Calculate the saturation vapor pressure
!
!-----------------------------------------------------------------------
!
!  Author: Steve Lazarus
!          6/10/93
!
!  Modification History:
!
!-----------------------------------------------------------------------
!
!  Input :
!
!    pdim     Sounding pressure (Pa)
!    pnon     Sounding pressure (nondimensional)
!    t        Sounding Potential temperature at level k (K)
!
!  Output:
!
!    qvs      Saturation vapor pressure
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------
!

  REAL :: pdim,pnon,t
!
  a = 7.5*LOG(10.)
  b = 3.8/pdim
!
  es  = 610.78*EXP(a*(pnon*t - 273.16)/(pnon*t - 35.86))
  qvs = 0.622*es/(pdim-es)
!
! end of function qvs
  END FUNCTION qvs

