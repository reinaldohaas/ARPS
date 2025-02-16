!
!################################################################
!################################################################
!#####                                                      #####
!#####               SUBROUTINE ASSIMFIL                    #####
!#####                                                      #####
!#####               Copyright (c) 1993                     #####
!#####     Center for Analysis and Prediction of Storms     #####
!#####                University of Oklahoma                #####
!#####                                                      #####
!################################################################
!################################################################
!

SUBROUTINE assimfil(nx,ny,nz,x,y,z,zp,ptol,                             &
           u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,                    &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,j1,j2,j3,                  &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           tem9,tem10,tem11,tem12)
!
!---------------------------------------------------------------------
!
!  PURPOSE:
!
!---------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steve Lazarus and Alan Shapiro
!  9/29/93
!
!  MODIFICATION HISTORY:
!
!  01/16/96 (Limin Zhao)
!  Strip the hole-fill code from former 'assimvel.f' and modified
!  it into a independent subroutine.
!
!  03/08/96 (Limin Zhao)
!  Added checks for processing real data.
!
!  03/20/96 (Limin Zhao and Alan Shapiro)
!  Modified the code to use background information as bounday
!  condition.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state air density (kg/m**3) times j3
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!    tem7     Observed velocity u-component (m/s)
!    tem8     Observed velocity v-component (m/s)
!    tem9     Observed velocity w-component (m/s)
!
!---------------------------------------------------------------------
!
!  OUTPUT:
!
!    u        hole-filled x component of velocity (m/s)
!    v        hole-filled y component of velocity (m/s)
!    w        hole-filled Vertical component of velocity
!             in Cartesian coordinates (m/s)
!
!---------------------------------------------------------------------
!
!  WORK ARRAYS:
!
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!
!    (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
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
  IMPLICIT NONE     ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Model physical constants
  INCLUDE 'globcst.inc'     ! Model control constants
  INCLUDE 'assim.inc'       ! Assim/Retr control parameters
  INCLUDE 'bndry.inc'       ! Boundary condition parameters
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
  INTEGER :: nt                ! The no. of t-levels of t-dependent arrays
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present time.
  INTEGER :: tfuture           ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

  INTEGER :: nx, ny, nz        ! Number of grid points in x, y and z directions

  INTEGER :: grdbas

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz,nt)  ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz,nt)  ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz,nt)  ! Hail mixing ratio (kg/kg)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)

  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  REAL :: tem4(nx,ny,nz)
  REAL :: tem5(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)
  REAL :: tem7(nx,ny,nz)
  REAL :: tem8(nx,ny,nz)
  REAL :: tem9(nx,ny,nz)
  REAL :: tem10(nx,ny,nz)
  REAL :: tem11(nx,ny,nz)
  REAL :: tem12(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: icount ! Parameter input flag
  SAVE icount
  DATA icount /1/

  INTEGER :: i, j, k, ir, jr
  INTEGER :: istor,jstor,kstor
  INTEGER :: nchout    ! I/O Channel of unformatted binary restart data
  INTEGER :: it,ix
  INTEGER :: itmax     ! Maximum number of iterations
  INTEGER :: istat
  INTEGER :: isrc
  INTEGER :: tim
  INTEGER :: idummy
  INTEGER :: itag      ! Counter indicating input file date/time
  INTEGER :: count
  INTEGER :: mdpt

  REAL :: t1           ! Time of first data file
  REAL :: t2           ! Time of second data file
  REAL :: t3           ! Time of third data file

  REAL :: assimtim(100)! Time of (all) input data files

  REAL :: xs           ! Scalar coordinate x-component
  REAL :: ys           ! Scalar coordinate y-component
  REAL :: zs           ! Scalar coordinate z-component
  REAL :: xs1          ! Scalar coordinate x-component
  REAL :: ys1          ! Scalar coordinate y-component
  REAL :: zs1          ! Scalar coordinate z-component

  REAL :: xmove        ! (x,y) position of radar relative to
  REAL :: ymove        ! moving grid (grid origin at (0,0,0))

  REAL :: znz2         ! = z(nz-2) - zshift
  REAL :: znz1         ! = z(nz-1) - zshift
  REAL :: z2           ! = z(2)    - zshift

  REAL :: ptol         ! Error tolerance for hole-filling

  REAL :: rad          ! Radius of sphere in Cartesian coord
  REAL :: radh         ! Radius of circle in hz. plane
  REAL :: radinv       ! Inverse radius magnitude

  REAL :: umean,vmean,usum
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!  external pois3d
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Hole-fill u,v,w, and vr outside the rainwater field.  Hole-filling
!  is performed on the perturbation velocity fields only. tem4 serves
!  as a 'template' - delineating the region to be filled.  The
!  Dirichlet boundary conditions are set here, outside of POIS3D.
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) 'code in assimfil; spval=',spval
  WRITE(6,*) 'adas mean wind',ubar(5,5,5),vbar(5,5,5)
  umove=0.0
  vmove=0.0
  xmove= xshift - umove*(curtim-assimtim(1))
  ymove= yshift - vmove*(curtim-assimtim(1))

  WRITE(6,*) 'dx,dy,dz: ',dx,dy,dz
  WRITE(6,*) 'xshift,yshift: ',xshift,yshift

  tim = tpresent

!   DO 128 k=1,nz
!   DO 128 j=1,ny
!   DO 128 i=1,nx
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = 0.5*(u(i,j,k,tim)+u(i+1,j,k,tim))
        tem2(i,j,k) = 0.5*(v(i,j,k,tim)+v(i,j+1,k,tim))
        IF(tem7(i,j,k) /= spval) THEN
          tem10(i,j,k) = tem7(i,j,k) - tem1(i,j,k)
        ELSE
          tem10(i,j,k) = tem7(i,j,k)
        END IF
        IF(tem8(i,j,k) /= spval) THEN
          tem11(i,j,k) = tem8(i,j,k) - tem2(i,j,k)
        ELSE
          tem11(i,j,k) = tem8(i,j,k)
        END IF
        tem12(i,j,k) = tem9(i,j,k)

      END DO
    END DO
  END DO

  count = 0
  DO k=1,nz-1
    DO j=1,ny-1
!   DO 140 i=1,nx
      DO i=1,nx-1    ! at scalar points
        tem3(i,j,k) = 0.0
        tem6(i,j,k) = 0.0
        tem7(i,j,k) = 0.0
        IF(tem10(i,j,k) /= spval) THEN     ! u component
          tem5(i,j,k) = tem10(i,j,k)
          tem4(i,j,k)  = 0.0
!     ELSE IF(i.le.2.or.i.ge.nx-1.or.j.eq.1.or.j.eq.ny-1) THEN
        ELSE IF(i == 1.OR.i == nx-1.OR.j == 1.OR.j == ny-1) THEN
          tem5(i,j,k) = 0.0  !Dirichlet for lateral, gradient for vertical
          tem4(i,j,k) = 0.0
        ELSE                 ! Fill u outside rain regions
          tem5(i,j,k) = 0.0
          tem4(i,j,k) = spval
          count = count + 1
        END IF
      END DO
    END DO
  END DO

  PRINT *,'On call to POIS3D, there are ',count,' filled U values'

  CALL pois3d(nx,ny,nz,dx,dy,dz,ptol,2.0,tem3,tem4,tem5,tem6,tem7)
!   CALL POIS3D(nx,ny,nz,dx,dy,dz,ptol,3.0,tem3,tem4,tem5,tem6,tem7)

  DO k=1,nz-1
    DO j=1,ny-1
!   DO 141 i=1,nx
      DO i=1,nx-1
        tem10(i,j,k) = tem5(i,j,k) + tem1(i,j,k)
      END DO
    END DO
  END DO

  count = 0
  DO k=1,nz-1
!   DO 142 j=1,ny
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k) = 0.0
        tem6(i,j,k) = 0.0
        tem7(i,j,k) = 0.0
        IF(tem11(i,j,k) /= spval) THEN     ! v component
          tem5(i,j,k)  = tem11(i,j,k)
          tem4(i,j,k)  = 0.0
!     ELSE IF(i.eq.1.or.i.eq.nx-1.or.j.le.2.or.j.ge.ny-1) THEN
        ELSE IF(i == 1.OR.i == nx-1.OR.j == 1.OR.j == ny-1) THEN
          tem5(i,j,k)  = 0.0 !Dirichlet for lateral, gradient for vertical
          tem4(i,j,k)  = 0.0
        ELSE                 ! Fill v outside rain regions
          tem5(i,j,k) = 0.0
          tem4(i,j,k) = spval
          count = count + 1
        END IF
      END DO
    END DO
  END DO

  PRINT *,'On call to POIS3D, there are ',count,' filled V values'

  CALL pois3d(nx,ny,nz,dx,dy,dz,ptol,2.0,tem3,tem4,tem5,tem6,tem7)
!   CALL POIS3D(nx,ny,nz,dx,dy,dz,ptol,3.0,tem3,tem4,tem5,tem6,tem7)

  DO k=1,nz-1
!   DO 143 j=1,ny
    DO j=1,ny-1
      DO i=1,nx-1
        tem11(i,j,k) = tem5(i,j,k) + tem2(i,j,k)
      END DO
    END DO
  END DO

  count = 0
!   DO 144 k=1,nz
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k) = 0.0
        tem6(i,j,k) = 0.0
        tem7(i,j,k) = 0.0
        IF(tem12(i,j,k) /= spval) THEN   ! w component
          tem5(i,j,k) = tem12(i,j,k)
          tem4(i,j,k) = 0.0
!ccxxx        ELSE IF(k.le.2.or.k.ge.nz-1) THEN
!ccxxx          tem5(i,j,k) = 0.0
!ccxxx          tem4(i,j,k) = 0.0
!     ELSE IF(k.le.2.or.k.ge.nz-1) THEN
!       tem5(i,j,k) = 0.0
!       tem4(i,j,k) = 0.0
        ELSE                          ! Fill w outside rain regions
          tem5(i,j,k) = 0.0
          tem4(i,j,k) = spval
          count = count + 1
        END IF
      END DO
    END DO
  END DO

  PRINT *,'On call to POIS3D, there are ',count,' filled W values'

!   CALL POIS3D(nx,ny,nz,dx,dy,dz,ptol,2.0,tem3,tem4,tem5,tem6,tem7)
!   CALL POIS3D(nx,ny,nz,dx,dy,dz,ptol,3.0,tem3,tem4,tem5,tem6,tem7)

!   DO 145 k=1,nz
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem12(i,j,k) = tem5(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,k) = tem10(i,j,k)
        tem8(i,j,k) = tem11(i,j,k)
        tem9(i,j,k) = tem12(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE assimfil
