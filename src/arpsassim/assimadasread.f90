!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADASREAD                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE adasread(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                    &
           isrc,ireturn,adastim,                                        &
           ubar,vbar,pbar,ptbar,rhostr,qvbar,                           &
           u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,                 &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           tem9,rhobar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinates the ingestion of ADAS background data files.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Limin Zhao
!
!  2/21/1996
!
!  MODIFICATION HISTORY:
!
!  2/17/97 (L. Zhao)
!  Modified the code to ingest the ARPS4.2.4 history format.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Height of soil levels
!
!  OUTPUT:
!
!    ireturn  Flag, indicating read status of input file
!             = 0 Successful read
!             = 1 Unsuccessful read
!
!    adastim  Times of all input files
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state density (kg/m**3) times j3
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    rhobar   Base state density (kg/m**3)
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
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
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
  REAL :: zpsoil(nx,ny,nzsoil) ! Height of soil level

  INTEGER :: isrc             ! Flag indicating source of calling routine
  INTEGER :: nstyps
  REAL :: adastim              ! Time of data input files

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
!  real kte   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)

!   real kmh   (nx,ny,nz)     !Horizontal turb. mixing coef.for momentum ( m**2/s )
!   real kmv   (nx,ny,nz)     !Vertical turb. mixing coef. for momentum ( m**2/s )
  REAL :: j1(nx,ny,nz)
  REAL :: j2(nx,ny,nz)
  REAL :: j3(nx,ny,nz)


  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array

  INTEGER :: ireturn

  REAL, ALLOCATABLE :: tem4dsoil(:,:,:,:)
  REAL, ALLOCATABLE :: tem3dsoil(:,:,:)

!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL dtahead
  EXTERNAL dtaread
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
  INTEGER :: itimesv

  REAL :: storstop         ! Temporarily stores the model stop time

  CHARACTER (LEN=80  ) :: saverunm  ! Temporarily stores the name of this run
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'assim.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'adas.inc'
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
!  If isrc=1, then read the times (adastim) from ADAS data file
!  headers. This is done in order to determine if the background data
!  is available.
!
!-----------------------------------------------------------------------
!
  WRITE(6,*)'code in adasread'

  tim = tpresent

  IF (isrc == 1) THEN

    lendtf=LEN(adasdat(is))
    CALL strlnth( adasdat(is), lendtf)
    WRITE(6,'(/a,a)')'The file name is ',adasdat(is)(1:lendtf)

    saverunm = runname

    CALL dtahead(nx,ny,nz,                                              &
         inifmt,inigbf,lengbf,adasdat(is),                              &
         lendtf,adastim,x,y,z,zp,tem4,tem5,tem6,                        &
         ptprt(1,1,1,tim),pprt(1,1,1,tim),tem7,                         &
         qc(1,1,1,tim),qr(1,1,1,tim),qi(1,1,1,tim),                     &
         qs(1,1,1,tim),qh(1,1,1,tim),tem9,                              &
         ubar,vbar,tem8,ptbar,pbar,rhobar,qvbar,                        &
         ireturn,tem1,tem2,tem3)

    runname = saverunm

    IF(ireturn /= 0) RETURN
!
!-----------------------------------------------------------------------
!
!  Read in a single data file.
!
!-----------------------------------------------------------------------
!
  ELSE IF(isrc == 3) THEN


    lengbf=LEN(inigbf)
    CALL strlnth( inigbf, lengbf)
    WRITE(6,'(/a,a)')                                                   &
        'The grid/base state file name is ',inigbf(1:lengbf)

    lendtf=LEN(adasdat(is))
    CALL strlnth( adasdat(is), lendtf)
    WRITE(6,'(/a,a)')'The file name is ',adasdat(is)(1:lendtf)

    saverunm = runname

    WRITE(6,*)'ADAS background field',adasdat(is)
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
        END DO
      END DO
    END DO

    storstop = tstop  ! Temporary. For ingestion of model data
                      ! generated prior to 4.0

    CALL ctim2abss( year, month, day, hour, minute, second, itimesv )

    ALLOCATE(tem4dsoil(nx,ny,nzsoil,0:nstyps))
    ALLOCATE(tem3dsoil(nx,ny,nstyps))

    CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
         inifmt,nchanl,inigbf,lengbf,adasdat(is),                       &
         lendtf,adastim,x,y,z,zp,zpsoil,tem4,tem5,tem6,                 &
         ptprt(1,1,1,tim),pprt(1,1,1,tim),tem7,                         &
         qc(1,1,1,tim),qr(1,1,1,tim),qi(1,1,1,tim),                     &
         qs(1,1,1,tim),qh(1,1,1,tim),tem9,tem9,tem9,                    &
         ubar,vbar,tem8,ptbar,pbar,rhobar,qvbar,                        &
         tem3dsoil,tem3dsoil,tem9(1,1,1),                               &
         tem9(1,1,1),tem9(1,1,1),tem9(1,1,1),                           &
         tem4dsoil,tem4dsoil,tem9(1,1,1),tem9(1,1,1),                   &
         tem9(1,1,1),tem9(1,1,1),tem9(1,1,1),                           &
         tem9(1,1,1),tem9(1,1,1),tem9(1,1,1),                           &
         tem9(1,1,1),tem9(1,1,1),                                       &
         tem9(1,1,1),tem9(1,1,1),tem9(1,1,1),tem9(1,1,1),               &
         ireturn,tem1,tem2,tem3)

    DEALLOCATE(tem4dsoil)
    DEALLOCATE(tem3dsoil)

    CALL abss2ctim( itimesv, year, month, day, hour, minute, second )
    tstop   = storstop
    runname = saverunm

    IF(ireturn /= 0) RETURN
!
!-----------------------------------------------------------------------
!
!  The arrays tem2,tem3 and tem4 contain uprt, vprt, wprt
!  respectively. The total fields are stored in u and v.
!
!-----------------------------------------------------------------------
!
    DO i = 1,nx
      DO j = 1,ny-1
        DO k = 1,nz-1
          tem2(i,j,k)  = tem4(i,j,k) + ubar(i,j,k)
          u(i,j,k,tim) = tem2(i,j,k)
        END DO
      END DO
    END DO

    DO i = 1,nx-1
      DO j = 1,ny
        DO k = 1,nz-1
          tem3(i,j,k)  = tem5(i,j,k)+ vbar(i,j,k)
          v(i,j,k,tim) = tem3(i,j,k)
        END DO
      END DO
    END DO

    DO i = 1,nx-1
      DO j = 1,ny-1
        DO k = 1,nz
          tem4(i,j,k)  = tem6(i,j,k)
          w(i,j,k,tim) = tem4(i,j,k)
        END DO
      END DO
    END DO

    DO i = 1,nx-1
      DO j = 1,ny-1
        DO k = 1,nz-1
          rhostr(i,j,k) = rhobar(i,j,k)*j3(i,j,k)
          qv(i,j,k,tim) = qvbar(i,j,k) + tem7(i,j,k)
        END DO
      END DO
    END DO

  END IF

  RETURN


END SUBROUTINE adasread
