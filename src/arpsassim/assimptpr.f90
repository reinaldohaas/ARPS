!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE ASSIMPTPR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimptpr(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                   &
           u,v,w,ptprt,pprt,qv,                                         &
           qc,qr,qi,qs,qh,                                              &
           tke,kmh,kmv,                                                 &
           ubar,vbar,ptbar,pbar,qvbar,rhostr,                           &
           uforce,vforce,wforce,j1,j2,j3,j3soil)
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is the driver for retrieval and adjustment
!  (e.g., through blending of retreived fields) of the pressure
!  and temperature fields.
!
!  If any changes to temperature or pressure are generated
!  they are applied to the input fields ptprt and pprt before
!  returning.
!
!  Note: Any statistics or processing involving both the original
!  and adjusted mass fields must be done before exiting this subroutine.
!
!---------------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2/25/93.
!
!  MODIFICATION HISTORY:
!    04/15/93 (Keith Brewster)
!    Modified for new subroutine RETRPTPR and to use common block
!    recvry (assim.inc).
!
!    02/08/96 (Scott Ellis and Steven Lazarus)
!    Modified to dump out retrieved pressure and temperature following
!    the thermodynamic recovery. (Call to DTADUMP)
!
!    03/22/96 (Limin Zhao)
!    Modified to incoporate it into the revised FDDA sysytem.
!
!    03/22/96 (Limin Zhao)
!    Added the subroutines for reading radar data, for blending
!    retrieval with model background.
!
!    05/06/96 (Limin Zhao)
!    Added modified Keith's output code to convert the retrieved ,
!    wind field, thermodynamic and moisture fields into
!    "pseudo-soundings" for use in ADAS.
!
!    08/22/96 (Limin Zhao)
!    Added swatch 'itfil' to turn on/off hole-filler after retrieval T.
!
!    08/23/96 (Limin Zhao)
!    Added working array 'work*' to temporally fix the problem in
!    calling 'assimout'. These arrays are only used in this subroutine,
!    and will need be changed according to 'dims.inc'.
!
!    09/24/96 (Limin Zhao)
!    Fixed a bug in 'assimout'. The retrieval temperature and pressure
!    should be dumped out, and also should be qvprt not qv.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Height of soil level (m)
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        z component of velocity (m/s)
!
!    ptprt    Model's perturbation potential temperature (K)
!    pprt     Model's perturbation pressure (Pascal)
!    qv       Model's water vapor specific humidity (kg/kg)
!    qc       Model's cloud water mixing ratio (kg/kg)
!    qr       Model's rain water mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum. ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum. ( m**2/s )
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    rhostr   Base state air density (kg/m**3) times j3
!
!    uforce   Acoustically inactive forcing terms in the u-momentum
!             equation (kg/(m*s)**2). uforce = -uadv + umix + ucorio
!    vforce   Acoustically inactive forcing terms in the v-momentum
!             equation (kg/(m*s)**2). vforce = -vadv + vmix + vcorio
!    wforce   Acoustically inactive forcing terms in the w-momentum
!             equation, except for buoyancy (kg/(m*s)**2).
!             wforce = -wadv + wmix + wcorio
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz

!  OUTPUT:
!
!    ptprt    Adjusted perturbation potential temperature (K)
!    pprt     Adjusted perturbation pressure (Pascal)
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
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  IMPLICIT NONE                ! Force explicit declarations
  
  INTEGER :: nx, ny, nz        ! Number of grid points in the x,y,z-direction
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nt                ! Number of time levels of time-dependent arrays.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present time.
  INTEGER :: tfuture           ! Index of time level for the future time.
  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

  INTEGER :: nstyps

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-poi:.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-poi:.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-poi: on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-poi: of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! Soil level height.

  REAL :: u(nx,ny,nz,nt)       ! x compone: of velocity (m/s)
  REAL :: v(nx,ny,nz,nt)       ! y compone: of velocity (m/s)
  REAL :: w(nx,ny,nz,nt)       ! z compone: of velocity (m/s)

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation pote:ial temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz,nt)  ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz,nt)  ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz,nt)  ! Hail mixing ratio (kg/kg)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbule: Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizo:al turb. mixing coef. for
                               ! mome:um. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! mome:um. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state pote:ial temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: qvbar (nx,ny,nz)     ! Base state specific humidity (kg/kg).
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-mome:um equation (kg/(m*s)**2)
                               ! uforce= -uadv + umix + ucorio

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-mome:um equation (kg/(m*s)**2)
                               ! vforce= -vadv + vmix + vcorio

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-mome:um equation (kg/(m*s)**2)
                               ! wforce= -wadv + wmix

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian  d(zpsoil)/d(z)

  REAL, ALLOCATABLE :: tem1  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem2  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem3  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem4  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem5  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem6  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem7  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem8  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem9  (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem10  (:,:,:)    ! Temporary work array.
  REAL, ALLOCATABLE :: tem11  (:,:,:)    ! Temporary work array.

! work arrays used by assimout:
  REAL, ALLOCATABLE :: work1(:,:,:)
  REAL, ALLOCATABLE :: work2(:,:,:)
  REAL, ALLOCATABLE :: work3(:,:,:)
  REAL, ALLOCATABLE :: work4(:,:,:)
  REAL, ALLOCATABLE :: work5(:,:,:)
  REAL, ALLOCATABLE :: work6(:,:,:)
  REAL, ALLOCATABLE :: work7(:,:,:)
  REAL, ALLOCATABLE :: work8(:,:,:)
  REAL, ALLOCATABLE :: work9(:,:,:)
  REAL, ALLOCATABLE :: work10(:,:,:)
  REAL, ALLOCATABLE :: work11(:,:,:)

  REAL, ALLOCATABLE :: tem4dsoil(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,tlevel
  INTEGER :: rebc,rwbc,rnbc,rsbc,rtbc,rbbc
  REAL :: dtbig1,dtsml1

  INTEGER :: kdim
  PARAMETER (kdim=100)
  INTEGER :: knt(kdim)
  REAL :: vmin(kdim),vmax(kdim)
  REAL :: avg(kdim),rms(kdim),stdv(kdim)
  INTEGER :: knttot
  REAL :: mintot,maxtot
  REAL :: avgtot,rmstot,stdvtot
  REAL :: ptol

  INTEGER :: grdbas

  REAL :: assimtim(100)                ! Time of (all) input data files

  INTEGER :: icount
  INTEGER :: istat
  INTEGER :: isrc
  INTEGER :: lenstr

  INTEGER :: tim
  INTEGER :: nchout

  INTEGER :: lendtf,ireturn

  REAL :: xor                ! Coordinate (xor,yor,zor) of the model grid
  REAL :: yor                ! origin relative to the radar.
  REAL :: zor                !

  CHARACTER (LEN=256) :: retfname

  INTEGER :: itfil
  REAL :: pres,temp,qvsat
  INTEGER :: count, istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'      ! Global constants that control model
  INCLUDE 'bndry.inc'
  INCLUDE 'assim.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL retrptpr
  EXTERNAL pois3d
  EXTERNAL radptpr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  If recovopt equals 0, do not perform a retrieval.
!
!-----------------------------------------------------------------------
!
  IF ( recovopt == 0 ) RETURN
!
!-----------------------------------------------------------------------
!
!  Allocate the variables and initialize the them to zero
!
!-----------------------------------------------------------------------


  ALLOCATE( tem1  (nx,ny,nz), STAT=istatus)
  tem1 = 0
  ALLOCATE( tem2  (nx,ny,nz), STAT=istatus)
  tem2 = 0
  ALLOCATE( tem3  (nx,ny,nz), STAT=istatus)
  tem3 = 0
  ALLOCATE( tem4  (nx,ny,nz), STAT=istatus)
  tem4 = 0
  ALLOCATE( tem5  (nx,ny,nz), STAT=istatus)
  tem5 = 0
  ALLOCATE( tem6  (nx,ny,nz), STAT=istatus)
  tem6 = 0
  ALLOCATE( tem7  (nx,ny,nz), STAT=istatus)
  tem7 = 0
  ALLOCATE( tem8  (nx,ny,nz), STAT=istatus)
  tem8 = 0
  ALLOCATE( tem9  (nx,ny,nz), STAT=istatus)
  tem9 = 0
  ALLOCATE( tem10 (nx,ny,nz), STAT=istatus)
  tem10 = 0
  ALLOCATE( tem11 (nx,ny,nz), STAT=istatus)
  tem11 = 0


  ALLOCATE( work1 (nx,ny,nz), STAT=istatus)
  work1 = 0
  ALLOCATE( work2 (nx,ny,nz), STAT=istatus)
  work2 = 0
  ALLOCATE( work3 (nx,ny,nz), STAT=istatus)
  work3 = 0
  ALLOCATE( work4 (nx,ny,nz), STAT=istatus)
  work4 = 0
  ALLOCATE( work5 (nx,ny,nz), STAT=istatus)
  work5 = 0
  ALLOCATE( work6 (nx,ny,nz), STAT=istatus)
  work6 = 0
  ALLOCATE( work7 (nx,ny,nz), STAT=istatus)
  work7 = 0
  ALLOCATE( work8 (nx,ny,nz), STAT=istatus)
  work8 = 0
  ALLOCATE( work9 (nx,ny,nz), STAT=istatus)
  work9 = 0
  ALLOCATE( work10(nx,ny,nz), STAT=istatus)
  work10 = 0
  ALLOCATE( work11(nx,ny,nz), STAT=istatus)
  work11 = 0


!
!-----------------------------------------------------------------------
!
!  Has switch been set for thermodynamic retrieval?
!
!-----------------------------------------------------------------------
!
  IF (irecov == 1) THEN
!
!----------------------------------------------------------------------------
!
!  Fudge-up a first guess.  This is just to compare stats
!  with the final product.
!
!----------------------------------------------------------------------------

!
!----------------------------------------------------------------------------
!
!  For real data runs, read the raw data as flags to do hole-filling.
!
!----------------------------------------------------------------------------
!
    isrc=3     !NOTE: temporary fix for 3 files case.
    ii=3
    lendtf=LEN_trim(assimdat(ii-1))
!
!
!    write(6,'(/a,a)')'The file name is ',assimdat(i)(1:lendtf)
    WRITE(6,*) 'lendtf=',lendtf
    WRITE(6,*) 'The filename assimdat(',ii-1,')=',                      &
                             assimdat(ii-1)(1:lendtf)

    CALL radptpr(assimdat(ii-1),lendtf,isrc,                            &
                 assimtim(ii-1),ireturn,                                &
                 nx,ny,nz,dx,dy,dz,xor,yor,zor,tem10,tem9)
!
!--------------------------------------------------------------------------
!
!  Recover potential temperature perturbation (tem1) and
!  perturbation pressure (tem2) from the observed winds and
!  the model forcing.
!
!-----------------------------------------------------------------------
!
    PRINT *, '  Calling retrptpr'

    CALL retrptpr(nx,ny,nz,x,y,z,zp,                                    &
               u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,                      &
               ubar,vbar,ptbar,pbar,rhostr,qvbar,                       &
               uforce,vforce,wforce,j1,j2,j3,                           &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,                      &
               tem8,tem9,tem10,tem11)

    PRINT *, '  Back from retrptpr'
!
!-----------------------------------------------------------------------
!
!  Dump out the recovered temperature (tem1), pressure (tem2),
!  the velocity, and moisture fields at tpresent to file hdmpfn.
!
!-----------------------------------------------------------------------
!
    tim = tpresent

    lenstr = LEN(assimnm)
    CALL strlnth( assimnm, lenstr )
    CALL gtdmpfn(assimnm(1:lenstr),dirnam,ldirnm,curtim,                &
                 hdmpfmt,mgrid,nestgrd, hdmpfn, ldmpf)

    WRITE(6,'(1x,a,a)')                                                 &
         'Retrieved data dump in file ',hdmpfn(1:ldmpf)

    ALLOCATE( tem4dsoil(nx,ny,nzsoil,0:nstyps), STAT=istatus)
    tem4dsoil = 0

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                &
         hdmpfmt,nchout,hdmpfn(1:ldmpf),grdbas,filcmprs,                &
         u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),tem1,                   &
         tem2,qv(1,1,1,tim),qc(1,1,1,tim),qr(1,1,1,tim),                &
         qi(1,1,1,tim),qs(1,1,1,tim),qh(1,1,1,tim),                     &
         tke,kmh,kmv,                                                   &
         ubar,vbar,tem8,ptbar,pbar,tem9,qvbar,                          &
         x,y,z,zp,zpsoil,                                               &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem4dsoil,tem4dsoil,tem11,tem11,                               &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem11,tem11,                                                   &
         tem11,tem11,tem11,tem11,                                       &
         tem5,tem6,tem7)

    DEALLOCATE(tem4dsoil)

!
!----------------------------------------------------------------------------
!
!  For real data runs, or OSSE runs with simulated filled regions,
!  you might like to fill-in the retrieved temperature (tem1) at
!  non-data points. Here the raw data is read as flags(tem8).
!
!----------------------------------------------------------------------------
!ccxxx
    isrc=3     !NOTE: temporary fix for 3 files case.
    ii=3
    lendtf=LEN(assimdat(ii-1))
    CALL strlnth( assimdat(ii-1), lendtf)
!
!
!    write(6,'(/a,a)')'The file name is ',assimdat(i)(1:lendtf)
    WRITE(6,*) 'lendtf=',lendtf
    WRITE(6,*) 'The filename assimdat(',ii-1,')=',                      &
                             assimdat(ii-1)(1:lendtf)

    CALL radptpr(assimdat(ii-1),lendtf,isrc,                            &
                 assimtim(ii-1),ireturn,                                &
                 nx,ny,nz,dx,dy,dz,xor,yor,zor,tem8,tem9)
!
!-----------------------------------------------------------------------
!
!  assume qv=qvs in rainwater regions with updraft ( > 0.1 m/s).
!  In downdraft area, model background mean is used.
!
!  Calculate the saturation water vapor specific humidity,
!  qvs, from pressure and potential temperature using Teten's formula.
!
!
!     Added 01-2-99 by SSW for AMS99 testing
!
!-----------------------------------------------------------------------
!
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pres = pbar(i,j,k)+tem2(i,j,k)
          temp = (ptbar(i,j,k)+tem1(i,j,k))*((pres/p0)**rddcp)
          qvsat=(380. / pres) * & !Teten's formula
               EXP( 17.27 * (temp-273.15)/(temp-35.86) )
!
!SSW/05-22-99 Only saturate updrafts within rainwater region
          IF(qr(i,j,k,tim) > 0.0001  .AND.w(i,j,k,tim) >= 1.0) THEN
!cc       IF(qr(i,j,k,tim).gt.0.00002 ) THEN
!cc       IF(qr(i,j,k,tim).gt.0.00002
!cc  :         .and.w(i,j,k,tim).ge.0.1) THEN
            qv(i,j,k,tim) = qvsat
            work5(i,j,k) = qvsat
          ELSE
            work5(i,j,k) = -999.
          END IF

        END DO
      END DO
    END DO
!
!
!
!
!-----------------------------------------------------------------------
!
!  Dump out the recovered temperature (tem1), pressure (tem2), wind
!  and moisture fields into columns, which will be read and used in
!  ARPS DAS. You do not need this if you will not do re-analysis.
!
!-----------------------------------------------------------------------
!
!ccxxx  hardwired working arrays. You need change them
!ccxxx  using your (nx,ny,nz)       !!!
!ccxxx
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          work1(i,j,k) = u(i,j,k,tim)
          work2(i,j,k) = v(i,j,k,tim)
          work3(i,j,k) = tem1(i,j,k)
          work4(i,j,k) = tem2(i,j,k)
!cc   work5(i,j,k) = qv(i,j,k,tim) - qvbar(i,j,k)  ! pert qv is dumped
          work5(i,j,k) = work5(i,j,k) - qvbar(i,j,k)  ! pert qv is dumped
          work6(i,j,k) = qr(i,j,k,tim)
          work7(i,j,k) = ubar(i,j,k)
          work8(i,j,k) = vbar(i,j,k)
          work9(i,j,k) = ptbar(i,j,k)
          work10(i,j,k) = pbar(i,j,k)
          work11(i,j,k) = qvbar(i,j,k)
        END DO
      END DO
    END DO
!
!

    WRITE (*,*) "XXX BASSIM nx,ny,nz",nx,ny,nz
    CALL assimout(nx,ny,nz,retfname,                                    &
                  x,y,z,zp,                                             &
                  work1,work2,work3,work4,work5,work6,                  &
                  work7,work8,work9,work10,work11,                      &
                  tem3,tem5,tem6,tem7,tem8)  ! tem8 serves as flag

!
!-----------------------------------------------------------------------
!
!ccxxx    hole-filling the retrieval temperature and pressure
!ccxxx    at data void area
!ccxxx
    itfil=0
    IF(itfil == 0) GO TO 6789

    tim = tpresent
    ptol = 0.001
    icount = 0
    DO i= 1,nx-1
      DO j= 1,ny-1
        DO k= 1,nz-1
          tem3(i,j,k) = 0.0
          tem6(i,j,k) = 0.0
          tem7(i,j,k) = 0.0
          IF(tem8(i,j,k) == spval) THEN  ! Fill P and T outside of rain regions
            tem4(i,j,k) = spval
            tem1(i,j,k) = 0.0
            tem2(i,j,k) = 0.0
            icount = icount + 1
          END IF
        END DO
      END DO
    END DO

!   print *,'On call to POIS3D, there are ',icount,' filled T values'

    CALL pois3d(nx,ny,nz,dx,dy,dz,ptol,2.0,tem3,tem4,tem1,tem6,tem7)
    CALL pois3d(nx,ny,nz,dx,dy,dz,ptol,2.0,tem3,tem4,tem2,tem6,tem7)

    6789  CONTINUE
!
!ccxxx  You need this if you will turn off the blending, but you need
!ccxxx  turn it off if you will need blending
!ccxxx
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem3(i,j,k) = tem1(i,j,k)
          tem4(i,j,k) = tem2(i,j,k)
        END DO
      END DO
    END DO
!
!----------------------------------------------------------------------------
!
!  Blend retrieved temperatures with the original forecast.
!  In general blending will depend on the time of the model relative
!  to the data window.
!  Blended temperature and pressure go into tem3 and tem4 respectively.
!
!----------------------------------------------------------------------------
!ccxxx
!ccxxx     CALL assimblnd(nx,ny,nz,
!ccxxx     :                 ptprt(1,1,1,tpresent),tem1,tem8,
!ccxxx     :                 tem5,tem6,tem7,tem3)

!ccxxx        CALL assimblnd(nx,ny,nz,
!ccxxx     :                 pprt(1,1,1,tpresent),tem2,tem8,
!ccxxx     :                 tem5,tem6,tem7,tem4)
!
!     print *, '  Back from dablend'
!
!----------------------------------------------------------------------------
!
!  Determine boundary points.  As a temporary short-cut we
!  specify zero-gradient bc's in the event that radiation (outflow)
!  bc's have been requested.  In this case the temporary arrays passed
!  into the bc routines are dummy addresses and are not used.
!
!  Note: Otherwise, the time tendency are requested.
!
!----------------------------------------------------------------------------
!
    dtbig1=dtbig
    dtsml1=dtsml

    IF(ebc == 4) THEN
      rebc=3
    ELSE
      rebc=ebc
    END IF

    IF(wbc == 4) THEN
      rwbc=3
    ELSE
      rwbc=wbc
    END IF

    IF(nbc == 4) THEN
      rnbc=3
    ELSE
      rnbc=nbc
    END IF

    IF(sbc == 4) THEN
      rsbc=3
    ELSE
      rsbc=sbc
    END IF

    IF(tbc == 4) THEN
      rtbc=3
    ELSE
      rtbc=tbc
    END IF

    IF(bbc == 4) THEN
      rbbc=3
    ELSE
      rbbc=bbc
    END IF
!
!----------------------------------------------------------------------------
!
!   Zero normal gradient condition are used here.
!   tem5,tem6,tem7,tem8 not be used.
!
!----------------------------------------------------------------------------
!
    CALL bcsclr(nx,ny,nz, dtbig1, tem3, tem3, tem3,                     &
                 tem5, tem6, tem7, tem8,                                &
                 rebc,rwbc,rnbc,rsbc,rtbc,rbbc)

    CALL bcp   (nx,ny,nz, dtsml1, tem4,                                 &
                tem5, tem6, tem7, tem8,                                 &
                rebc,rwbc,rnbc,rsbc,rtbc,rbbc)

    PRINT *, '  Back from bc routines'
!
!----------------------------------------------------------------------------
!
!  Set pprt and ptprt at times past and future equal to their respective
!  distributions at the present time.
!
!----------------------------------------------------------------------------
!
    DO tlevel = 1, nt
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ptprt(i,j,k,tlevel)  = tem3(i,j,k)
            pprt(i,j,k,tlevel)   = tem4(i,j,k)
          END DO
        END DO
      END DO
    END DO

!
!
!-------------------------------------------------------------------------
!
!  Dump out the recovered temperature (tem1), pressure (tem2),
!  the velocity, and moisture fields at tpresent to file hdmpfn.
!
!-----------------------------------------------------------------------
!
    tim = tpresent

    lenstr = LEN(assimnm)
    CALL strlnth( assimnm, lenstr )
    CALL gtdmpfn(assimnm(1:lenstr),dirnam,ldirnm,curtim,                &
                 hdmpfmt,mgrid,nestgrd, hdmpfn, ldmpf)

    WRITE(6,'(1x,a,a)')                                                 &
         'Retrieved data dump in file ',hdmpfn(1:ldmpf)

    CALL dtadump(nx,ny,nz,nstyps,                                       &
         hdmpfmt,nchout,hdmpfn(1:ldmpf),grdbas,filcmprs,                &
         u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),ptprt(1,1,1,tim),       &
         pprt(1,1,1,tim),qv(1,1,1,tim),qc(1,1,1,tim),qr(1,1,1,tim),     &
         qi(1,1,1,tim),qs(1,1,1,tim),qh(1,1,1,tim),                     &
         tke,kmh,kmv,                                                   &
         ubar,vbar,tem8,ptbar,pbar,tem9,qvbar,                          &
         x,y,z,zp,                                                      &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,                                             &
         tem11,tem11,tem11,tem11,                                       &
         tem5,tem6,tem7)
!
!
!   write(6,*)'normally stop after assimptpr'
!   STOP
!
!SSW 5-18-99 comment out to run arps forward after assimilation
    WRITE(6,*) 'normally stop after assimptpr'
    STOP
!
!

  END IF
!-----------------------------------------------------------------------
!
!  deallocate the variables and return
!
!-----------------------------------------------------------------------

  DEALLOCATE( tem1 , STAT=istatus)
  DEALLOCATE( tem2 , STAT=istatus)
  DEALLOCATE( tem3 , STAT=istatus)
  DEALLOCATE( tem4 , STAT=istatus)
  DEALLOCATE( tem5 , STAT=istatus)
  DEALLOCATE( tem6 , STAT=istatus)
  DEALLOCATE( tem7 , STAT=istatus)
  DEALLOCATE( tem8 , STAT=istatus)
  DEALLOCATE( tem9 , STAT=istatus)
  DEALLOCATE( tem10, STAT=istatus)
  DEALLOCATE( tem11, STAT=istatus)


  DEALLOCATE( work1 , STAT=istatus)
  DEALLOCATE( work2 , STAT=istatus)
  DEALLOCATE( work3 , STAT=istatus)
  DEALLOCATE( work4 , STAT=istatus)
  DEALLOCATE( work5 , STAT=istatus)
  DEALLOCATE( work6 , STAT=istatus)
  DEALLOCATE( work7 , STAT=istatus)
  DEALLOCATE( work8 , STAT=istatus)
  DEALLOCATE( work9 , STAT=istatus)
  DEALLOCATE( work10, STAT=istatus)
  DEALLOCATE( work11, STAT=istatus)

  RETURN
END SUBROUTINE assimptpr
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RADPTPR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radptpr(filnam,lendtf,                                       &
           isrc,tim,ireturn,                                            &
           nx,ny,nz,dx,dy,dz,                                           &
           xor,yor,zor,                                                 &
           tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into channel nchanl as unformatted binary data.
!
!  All output data are located at the grid cell centers.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus
!  8/10/95.
!
!  MODIFICATION HISTORY:
!
!  22/02/96 (Limin Zhao)
!  Added ubar and vbar for mean velocity. Changes are also made to
!  incorporate the retrieval data.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    filnam  Name of the input file
!    nx      Number of grid points in the x-direction (east/west)
!    ny      Number of grid points in the y-direction (north/south)
!    nz      Number of grid points in the vertical
!
!    dx      x grid spacing
!    dy      y grid spacing
!    dz      z grid spacing
!
!    xor     the location of radar station in x-direction
!    yor     the location of radar station in y-direction
!    zor     the location of radar station in z-direction
!
!    tim     time of retrieval data
!    isrc    Location of calling routine
!
!  OUTPUT:
!
!    tem1   Observed radial velocity
!    tem6   Retrieved u-component
!    tem7   Retrieved v-component
!    tem8   Retrieved w-component
!    tem5   Observed reflectivity
!    ubar   Retrieved mean velocity of u-component
!    vbar   Retrieved mean velocity of v-component
!
!  WORK ARRAYS:
!
!    tem2      Temporary work array.
!    tem3      Temporary work array.
!    tem4      Temporary work array.

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

  INTEGER :: nx,ny,nz       ! Number of model grid points in 3 directions
  INTEGER :: nxr,nyr,nzr    ! Number of real data in 3 directions

  INTEGER :: lendtf,i,j,k,ireturn

  CHARACTER (LEN=256) :: filnam
  CHARACTER (LEN=10)  :: radarsp
  INTEGER :: iyrrr,imorr,idarr,iharr,imarr,isarr

  REAL :: dxr,dyr,dzr
  REAL :: dx,dy,dz
  REAL :: xor, yor, zor
  REAL :: bad
  REAL :: tim
  REAL :: umean, vmean

  REAL :: tem1  (nx,ny,nz)  ! Temporary work array
  REAL :: tem2  (nx,ny,nz)  ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: tref

  INTEGER :: l
  INTEGER :: nchanl            ! FORTRAN I/O channel number for output
  INTEGER :: istat
  INTEGER :: ierr
  INTEGER :: isrc

  LOGICAL :: fexist
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'assim.inc'
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
!  Assign a proper filename, and check to see if the file exist.
!
!  Note: Cray routines are used to force binary data file to be
!        in the IEEE format or compatable with IBM.
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(filnam(1:lendtf), '-F f77 -N ieee', ierr)

  CALL getunit(nchanl)

  WRITE(6,*) 'nchanl,filnam(1:lendtf): ',nchanl,filnam(1:lendtf)

  OPEN(UNIT=nchanl,FILE=filnam(1:lendtf),                               &
       STATUS='unknown', FORM='unformatted',IOSTAT= istat )

  IF( istat /= 0 ) GO TO 998

  100   CONTINUE

!
!-----------------------------------------------------------------------
!
!    Read data in.
!
!-----------------------------------------------------------------------
!
!
  READ(nchanl,ERR=110,END=120) radarsp,iyrrr,imorr,idarr,               &
                                       iharr,imarr,isarr

  READ(nchanl,ERR=110,END=120) tim

  WRITE(6,'(3x,a,f10.2,a)') 'Data read at ',tim,'seconds'

  IF(isrc == 1.OR.isrc == 2) THEN
    CALL retunit( nchanl )
    CLOSE(nchanl)
    ireturn=0
    RETURN
  END IF

  READ(nchanl,ERR=110,END=120) xor, yor, zor

  READ(nchanl,ERR=110,END=120) nxr,nyr,nzr,dxr,dyr,dzr


!
!-----------------------------------------------------------------------
!
!  Check to make sure the dimensions of the input data match that of
!  the current model stream.
!
!-----------------------------------------------------------------------
!
  IF((nx /= nxr).OR.(ny /= nyr).OR.(nz /= nzr)) THEN

    WRITE(6,'(a,/a,i5,a,i5,a,i5,/a,i5,a,i5,a,i5)')                      &
        ' Array dimension(s) of the input file inconsistent with ',     &
        ' model definitions, dimensions in input data were nx=',nxr,    &
        ', ny=',nyr,', nz=',nzr,' the model definitions were nx=',      &
        nx,' ny= ', ny, ' nz= ',nz
    WRITE(6,'(a)') ' Job stopped in subroutine RADPTPR.'
    STOP

  END IF

  IF(ABS(dx-dxr) > 0.1 .OR. ABS(dy-dyr) > 0.1 .OR. ABS(dz-dzr) > 0.1) THEN

    WRITE(6,'(a,/a,f10.2,a,f10.2,a,f10.2,/a,f10.2,a,f10.2,a,f10.2)')    &
        'Grid interval in the input data inconsisent with',             &
        'model definitions. In the input data dx=',dxr,                 &
        ', dy=',dyr,', dz=',dzr,' the model definitions were dx=',      &
        dx,' dy= ', dy, ' dz= ',dz
    WRITE(6,'(a)') ' Job stopped in subroutine RADPTPR.'
    STOP

  END IF
!
!-----------------------------------------------------------------------
!
!  Continue reading
!
!-----------------------------------------------------------------------
!
  READ(nchanl,ERR=110,END=120) tem2      ! input Vr

  READ(nchanl,ERR=110,END=120) tem2      ! input retrieved u component

  READ(nchanl,ERR=110,END=120) tem2      ! input retrieved v component

!   read(nchanl,err=110,end=120) tem1      ! input retrieved w component
!
!   read(nchanl,err=110,end=120) tem2      ! input Zr

  READ(nchanl,ERR=110,END=120) tem2      ! input retrieved w component

  READ(nchanl,ERR=110,END=120) tem1      ! input Zr

  CALL retunit(nchanl)

  CLOSE(nchanl)

  WRITE (*,*) "XXX RADPTPR have in nx,ny,nz,dx,dy,dz",                  &
              nx,ny,nz,dx,dy,dz
  WRITE (*,*) "XXX RADPTPR read in nx,ny,nz,dx,dy,dz,xor,yor,zor",      &
              nxr,nyr,nzr,dxr,dyr,dzr

!
!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!----------------------------------------------------------------------
!
  930   CONTINUE

  WRITE(6,'(/a/)') 'Reading was successfully in RADPTPR'
  ireturn = 0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!
  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in RADPTPR'
  ireturn=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in RADREAD'
  ireturn=2

  998   CONTINUE

  WRITE(6,'(/1x,a,/1x,a/)')                                             &
       'File '//filnam(1:lendtf)                                        &
       //' not found.',                                                 &
       'Program returned from RADAREAD.'
  ireturn=1

  RETURN
END SUBROUTINE radptpr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DABLEND                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dablend(nx,ny,nz,varorg,varnew,varbln,tem1)
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  A weighted average of the array varorg (.e.g., the orginal
!  model forecast) and the array varnew (e.g., the data or retreived
!  data field) is computed and output in array varbln.
!
!  A future modification of this routine may allow for the
!  relative weights to vary in space.
!
!--------------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: varorg(nx,ny,nz)
  REAL :: varnew(nx,ny,nz)
  REAL :: varbln(nx,ny,nz)
  REAL :: tem1  (nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Include file 'assim.inc' for ARPS
!
!  This file contains the control parameters
!  used by the data assimilation routines.
!
!  These parameters are allocated in named common blocks
!  therefore accessible to subroutines that include this file.
!
!-----------------------------------------------------------------------
!
!                        Copyright (c) 1993
!            Center for Analysis and Prediction of Storms
!######                University of Oklahoma                ######
!
!-----------------------------------------------------------------------
!
!  AUTH0R: Keith Brewster
!  02/15/1993
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  darunnm: a character string matching the "runname" used in
!  the creation of the files to be read in as data.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=128  ) :: darunnm  ! Name of run that produced data file
!
  REAL :: datime
!
  INTEGER :: datafmt ! Parameter specifying format of ingested data
                     ! = 1, unformatted binary data;
                     ! = 2, formatted ascii data;
                     ! = 3, NCSA HDF format data;
                     ! = 4, Compressed binary data;

  CHARACTER (LEN=12) :: dtagbfn  ! Ingest data grid/base file name.
  CHARACTER (LEN=12) :: datafnm  ! Ingest data file name.
  INTEGER :: ldtagbf   ! Length of the history data dump file name
  INTEGER :: ldtafnm   ! Length of the history data dump file name

  COMMON /dtaingst/ darunnm, datime, datafmt,                           &
                    dtagbfn, ldtagbf, datafnm, ldtafnm

  REAL :: dastart    ! start time (sec) of data assim
  REAL :: daend      ! end time (sec) of data assim
  INTEGER :: dtimwgt ! type of time weighting ramp
!                  = 1, weight new =1 (complete replacement).
!                  = 2, linear weight, new=0 at dastart, 1 at daend
!                  = 3, linear weight, new=1 at dastart, 0 at daend
!
  COMMON /dtassim/ dastart,daend,dtimwgt
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: dtwindw,wgtnew,wgtorg
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!    For current version, the fields are totally replaced by new fields.
!
!-----------------------------------------------------------------------
!
  dtwindw=(daend-dastart)
  IF(dtwindw == 0. .OR. dtimwgt == 1) THEN
    wgtnew=1.
  ELSE IF(dtimwgt == 2) THEN
    wgtnew=(curtim-dastart)/dtwindw
  ELSE IF(dtimwgt == 3) THEN
    wgtnew=(daend-curtim)/dtwindw
  END IF
  wgtorg=1.-wgtnew

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        varbln(i,j,k) = wgtorg*varorg(i,j,k) + wgtnew*varnew(i,j,k)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE dablend
!
