!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RETRPTPR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE retrptpr(nx,ny,nz,x,y,z,zp,                                  &
           u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,                          &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           uforce,vforce,wforce,j1,j2,j3,                               &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,                          &
           pc,tem9,flag,tem11)
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is the driver for a modified version of the Gal-Chen/
!  Hane thermodynamic recovery.  This technique derives the perturbation
!  pressure and perturbation potential temperature fields at time t from
!  distributions of u, v and w at three consecutive times: t-dtbig, t and
!  t+dtbig.  The technique was originally described in:
!
!       Gal-Chen, T., 1978: A method for the initialization of the
!       the anelastic equations: Implications for matching models with
!       observations, Mon. Wea. Rev., 587-606;
!
!       Hane, C. E., and B. C. Scott, 1978: Temperature and pressure
!       perturbations within convective clouds derived from detailed
!       air motion information: Preliminary testing. Mon. Wea. Rev.,
!       654-661.
!
!  The main difference between this recovery and the conventional
!  Gal-Chen/Hane scheme is that we make provision for the presence of
!  the perturbation pressure in the buoyancy term in the vertical
!  equation of motion.
!
!  IMPORTANT NOTICE:  The results of the recovery are stored in the
!  arrays tem1 (ptprt) and tem2 (pprt).  The original pprt and ptprt
!  fields are unchanged.  Results are provided at every computational
!  point but the user is required to ensure that the boundary conditions
!  are consistent with their run.  The ARPS subroutine assimptpr, serving
!  as a driver for this subroutine can take care of that.
!
!---------------------------------------------------------------------------
!
!  AUTHORS: Alan Shapiro and Steve Lazarus
!  2/16/93.
!
!  MODIFICATION HISTORY:
!    04/15/93 (Keith Brewster)
!    Results stored in arrays tem1 (ptprt) and tem2 (pprt).
!
!    03/15/96 (Limin Zhao)
!    Added the code to read radar data to flag retrieval
!    temperature and pressure fields.
!
!    04/23/96 (Limin Zhao and Alan Shapiro)
!    Modified the code to hole-fill force terms in data void area
!    instead of using hole-filled velocity values.
!
!    06/18/96 (Limin Zhao and Alan Shapiro)
!    Added an option for using Dirichlet boundary conditions
!    for pressure solver.
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
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state air density (kg/m**3) times j3
!    qvbar    Base state water vapor specific humidity (kg/kg)
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
!
!  OUTPUT:
!
!    tem1     Recovered perturbation potential temperature (K)
!    tem2     Recovered perturbation pressure (Pascal)
!
!  WORK ARRAYS:
!
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    pc       Temporary work array.
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
  IMPLICIT NONE        ! Force explicit declarations

  INTEGER :: isrc         ! Flag indicating source of calling routine.
  INTEGER :: nt           ! Number of time levels of time-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.
  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y, z directions

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: u(nx,ny,nz,nt)       ! x component of velocity (m/s)
  REAL :: v(nx,ny,nz,nt)       ! y component of velocity (m/s)
  REAL :: w(nx,ny,nz,nt)       ! z component of velocity (m/s)

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
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity (kg/kg)
  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
                               ! uforce= -uadv + umix + ucorio

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
                               ! vforce= -vadv + vmix + vcorio

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
                               ! wforce= -wadv + wmix

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)

  INTEGER :: nn
  PARAMETER (nn=500)        ! Dimension of the a(z), b(z) and f(z)
                            ! variables.  nz should be < nn.

  REAL :: f(nn)                ! Vertical structure function for the
                               ! recovered perturbation pressure,
                               ! pprt = solution of Poisson equation + f(z)

  REAL :: a(nn)                ! Coefficient a(z) in the ordinary differential
                               ! equation for f(z): df/dz + a(z)f + b(z) = 0.

  REAL :: b(nn)                ! Coefficient b(z) in the ordinary differential
                               ! equation for f(z): df/dz + a(z)f + b(z) = 0.

  REAL :: pc(nx,ny,nz)         ! a variable related to perturbation pressure,
                               ! the solution of a Poisson equation.

  REAL :: f0                   ! Constant of integration for f(z): f0=f(dz/2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array.
  REAL :: flag  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: tlevel

  INTEGER :: ic        ! Counter, used for reading input data
  SAVE    ic
  DATA    ic/1/

  INTEGER :: istat     ! Error flag for reading input data

  REAL :: dtinv        ! Reciprocal of dtbig
  REAL :: rdenom       ! Reciprocal number of grid points on a horizontal sfc.
  REAL :: dudt         ! Time rate of change of rhostr*u
  REAL :: dvdt         ! Time rate of change of rhostr*v
  REAL :: dwdt         ! Time rate of change of rhostr*w
  REAL :: dpdz         ! Vertical derivative of pc
  REAL :: wfave        ! Vertical average of wforce - dwdt
  REAL :: bave         ! Vertical average of individual terms comprising b(z)

  REAL :: b1, b2, b3, b4, b5       ! Individual terms comprising b(z).
  REAL :: term1(nn), term2(nn)     ! 1-D temporary arrays for use in RETRFZ

  REAL :: assimtim(100)            ! Time of data input files

  REAL :: ptol
  INTEGER :: icount

  REAL :: pzng,ptzng,tzng,qvszng,rhzng,dqvsdpt,faczng
  REAL :: corr,pterr
  INTEGER :: iter,itermax, k1
  REAL :: ubig,vbig
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'assim.inc'       ! Velocity insertion/thermodynamic
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL uvwrho
  EXTERNAL retrpois
  EXTERNAL retrfz
  EXTERNAL retrchk
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (recovopt == 0)  RETURN

  IF (irecov == 0) RETURN

  IF (nz > nn) THEN
    PRINT *, 'Warning from assimthermo: nz .gt.', nn
    PRINT*, 'please reset nn'
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate rhostr*u, rhostr*v and rhostr*w at two time levels: tpast
!  and tfuture.  Store the results in temporary arrays.
!
!-----------------------------------------------------------------------
!
!
!===> Use SDVR winds local time change
!
  tlevel = tpast

  CALL uvwrho(nx,ny,nz,                                                 &
              u(1,1,1,tlevel),v(1,1,1,tlevel),w(1,1,1,tlevel),          &
              rhostr,tem1,tem2,tem3)

  tlevel = tfuture

  CALL uvwrho(nx,ny,nz,                                                 &
              u(1,1,1,tlevel),v(1,1,1,tlevel),w(1,1,1,tlevel),          &
              rhostr,tem5,tem6,tem7)
!
!===> for turbulence frozen assumption.
!
!   tlevel = tpresent
!
!   CALL uvwrho(nx,ny,nz,
!  :            u(1,1,1,tlevel),v(1,1,1,tlevel),w(1,1,1,tlevel),
!  :            rhostr,tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  Calculate the (local) time derivative of rhostr*u, rhostr*v and
!  rhostr*w with a centered time difference.  Subtract these time
!  derivatives from uforce, vforce and wforce, and store the results
!  in temporary arrays.
!
!-----------------------------------------------------------------------
!

!   ubig =    8.3128        ! for 16:59Z
!   vbig = - 10.9770

  ubig =   7.9635        ! for 18:41Z
  vbig = - 7.2613

  dtinv = 0.5/dtbig

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-1
!
!
        dudt = dtinv*(tem5(i,j,k) - tem1(i,j,k))   ! SDVR wind
!      dudt = 0.0                                 ! stationary
!      dudt = -ubig*(tem1(i+1,j,k)-tem1(i-1,j,k))/(2.0*dx)
!  :          -vbig*(tem1(i,j+1,k)-tem1(i,j-1,k))/(2.0*dy)  ! trb. frzn
        uforce(i,j,k) = uforce(i,j,k) - dudt
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=2,nx-2
!
!
        dvdt = dtinv*(tem6(i,j,k) - tem2(i,j,k))   ! SDVR wind
!      dvdt = 0.0                                 ! stationary
!      dvdt = -ubig*(tem2(i+1,j,k)-tem2(i-1,j,k))/(2.0*dx)
!  :          -vbig*(tem2(i,j+1,k)-tem2(i,j-1,k))/(2.0*dy)  ! trb. frzn
        vforce(i,j,k) = vforce(i,j,k) - dvdt
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
!
!
        dwdt = dtinv*(tem7(i,j,k) - tem3(i,j,k))    ! SDVR wind
!     dwdt = 0.0                                  ! stationary
!      dwdt = -ubig*(tem3(i+1,j,k)-tem3(i-1,j,k))/(2.0*dx)
!  :          -vbig*(tem3(i,j+1,k)-tem3(i,j-1,k))/(2.0*dy)  ! trb. frzn
        wforce(i,j,k) = wforce(i,j,k) - dwdt
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k) = 0.0
        tem2(i,j,k) = 0.0
        tem3(i,j,k) = 0.0
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-1
        tem1(i,j,k) = uforce(i,j,k)
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=2,nx-2
        tem2(i,j,k) = vforce(i,j,k)
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
        tem3(i,j,k) = wforce(i,j,k)
      END DO
    END DO
  END DO
!
!----------------------------------------------------------------------------
!
!    Check if it is needed to use the hole-filled values in data void area.
!    itest=1, use the hole-filled A & B values.
!      =0, use the calculated values from hole-filled wind fields.
!
!    Note: The resultes are good for itest=1 and A & B are are set to zero
!
!---------------------------------------------------------------------------
!
!ccxxx
!
  IF(itest == 0) GO TO 911

  ptol = .001
  icount = 0
  DO k= 1,nz-1
    DO j= 1,ny-1
      DO i= 1,nx
        tem6(i,j,k) = 0.0
        tem7(i,j,k) = 0.0
        tem9(i,j,k) = 0.0
        IF(flag(i,j,k) /= spval) THEN
          tem11(i,j,k) = tem1(i,j,k)
!ccxxx        ELSE IF(i.le.1.or.i.ge.nx-1.or.j.eq.1.or.j.eq.ny-1) THEN
!ccxxx           tem11(i,j,k) =  tem1(i,j,k)
!ccxxx           tem11(i,j,k) =  0.0
!ccxxx           tem5(i,j,k) = 0.0
        ELSE
          tem11(i,j,k) =0.0
!ccxxx        tem11(i,j,k) = spval
          tem5(i,j,k) = spval
          icount = icount + 1
        END IF
      END DO
    END DO
  END DO

  PRINT *,'No call to POIS3D to fill A,  A=0 outside wind area'
!ccxxxprint *,'On call to POIS3D, there are ',icount,' filled A values'

!ccxxxCALL POIS3D(nx,ny,nz,dx,dy,dz,ptol,3.0,tem6,tem5,tem11,tem7,tem9)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        tem1(i,j,k) = tem11(i,j,k)
      END DO
    END DO
  END DO

  icount = 0
  DO k= 1,nz-1
    DO j= 1,ny
      DO i= 1,nx-1
        tem6(i,j,k) = 0.0
        tem7(i,j,k) = 0.0
        IF(flag(i,j,k) /= spval) THEN
          tem11(i,j,k) = tem2(i,j,k)
!ccxxx        ELSE IF(i.le.1.or.i.ge.nx-1.or.j.eq.1.or.j.eq.ny-1) THEN
!ccxxx        tem11(i,j,k) =  tem2(i,j,k)
!ccxxx        tem11(i,j,k) =  0.0
!ccxxx        tem5(i,j,k) = 0.0
        ELSE
          tem11(i,j,k) =0.0
!ccxxx        tem11(i,j,k) = spval
          tem5(i,j,k) = spval
          icount = icount + 1
        END IF

      END DO
    END DO
  END DO

  PRINT *,'No call to POIS3D to fill A,  A=0 outside wind area'
!ccxxxprint *,'On call to POIS3D, there are ',icount,' filled B values'

!ccxxxCALL POIS3D(nx,ny,nz,dx,dy,dz,ptol,3.0,tem6,tem5,tem11,tem7,tem9)

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        tem2(i,j,k) = tem11(i,j,k)
      END DO
    END DO
  END DO

  911   CONTINUE

!
!----------------------------------------------------------------------------
!
!  Compute the right hand side of the Poisson equation and store it in
!  tem4.  Recall that tem1 = uforce - dudt; tem2 = vforce - dvdt. The
!  right hand side is dtem1/dx + dtem2/dy.
!
!----------------------------------------------------------------------------
!

  DO k = 2, nz-2
    DO j = 2, ny-2
      DO i = 2, nx-2
        tem4(i,j,k) = dxinv*(tem1(i+1,j,k) - tem1(i,j,k)) +             &
                      dyinv*(tem2(i,j+1,k) - tem2(i,j,k))

      END DO
    END DO
  END DO
!
!----------------------------------------------------------------------------
!
!  Solve the two-dimensional (x-y) Poisson Equation for the variable pc:
!
!        d2pc/dx2 + d2pc/dy2 = dtem1/dx + dtem2/dy
!
!
!  with Neumann boundary conditions:
!
!        dpc/dx = tem1 on x faces
!        dpc/dy = tem2 on y faces
!
!  Recall that tem1 = uforce - dudt; tem2 = vforce - dvdt
!
!  pc differs from the perturbation pressure, pprt, by an arbitrary
!  function of height; that is,  pprt = pc + f(z)
!
!----------------------------------------------------------------------------
!
  tlevel = tpresent

  DO i = 1, nx
    DO j = 1, ny
!ccxxx        pc(i,j,1) = 0.                ! First guess of pc(i,j,1) is 0
      pc(i,j,1) = pprt(i,j,1,tlevel)
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem6(i,j,k) = pprt(i,j,k,tlevel)
      END DO
    END DO
  END DO

  DO k = 2, nz-2

    DO i = 1, nx
      DO j = 1, ny
!ccxxx          pc(i,j,k) = pc(i,j,k-1)   ! First guess of pc(i,j,k) is pc(i,j,k-1)
        pc(i,j,k) = pprt(i,j,k,tlevel)
      END DO
    END DO

    CALL retrpois(nx,ny,nz,k,                                           &
                     tem1,tem2,pc,                                      &
                     tem4,tem5(1,1,1),tem5(1,1,2),tem5(1,1,3),          &
                     tem5(1,1,4),tem6)
  END DO

!
!----------------------------------------------------------------------------
!
!  Check how well the retrieved pressure gradients fit the individual
!  momentum equations.
!
!----------------------------------------------------------------------------
!
  CALL retrchk(nx,ny,nz,tem1,tem2,pc)
!
!----------------------------------------------------------------------------
!
!  Calculate the coefficients a(z) and b(z) appearing in the ordinary
!  differential equation for f(z):  df/dz + a(z)f + b(z) = 0.
!
!  a = rhostr*g/(pbar*cpdcv)
!
!  b = horizontal average of (dpc/dz + a*pc - rhostr*g*ptpert/ptbar - tem3),
!
!  where tem3 = wforce - dwdt.
!
!  In the evaluation of the analytic solution for f(z), it is convenient
!  to define the a(z) coefficient at scalar points and the b(z) coefficient
!  at w points.
!
!  Use bc_opt to control what kind of boundary will be used.
!  bc_opt = 1, use Neumann boundary.
!         = 0, use Dirichlet boundary.
!
!----------------------------------------------------------------------------
!
!ccxxx
  IF (bc_opt == 0) GO TO 9119

  tlevel = tpresent

  rdenom = 1./((nx-3)*(ny-3))

  DO k = 3, nz-2
    a(k) = 0.
    b1 = 0.
    b2 = 0.
    b3 = 0.
    b4 = 0.
    b5 = 0.
    DO j = 2, ny-2
      DO i = 2, nx-2
        a(k) = a(k) + rhostr(i,j,k)*g/(pbar(i,j,k)*cpdcv)

        b1 = b1 + dzinv*(pc(i,j,k) - pc(i,j,k-1))

        bave = 0.5*(rhostr(i,j,k)  *pc(i,j,k)  /pbar(i,j,k) +           &
                    rhostr(i,j,k-1)*pc(i,j,k-1)/pbar(i,j,k-1))

        b2 = b2 + bave

        bave = 0.5*(rhostr(i,j,k)*ptprt(i,j,k,tlevel)/ptbar(i,j,k) +    &
                rhostr(i,j,k-1)*ptprt(i,j,k-1,tlevel)/ptbar(i,j,k-1))

        b3 = b3 + bave

        b4 = b4 + tem3(i,j,k)

        bave = 0.5*(rhostr(i,j,k)*qc(i,j,k,tlevel) +                    &
                    rhostr(i,j,k-1)*qc(i,j,k-1,tlevel) ) +              &
               0.5*(rhostr(i,j,k)*qr(i,j,k,tlevel) +                    &
                    rhostr(i,j,k-1)*qr(i,j,k-1,tlevel) ) -              &
            0.608*0.5*(rhostr(i,j,k)*(qv(i,j,k,tlevel)-qvbar(i,j,k)) +  &
                    rhostr(i,j,k-1)*(qv(i,j,k-1,tlevel)-qvbar(i,j,k)))

        b5 = b5 + bave

      END DO
    END DO

    a(k) = a(k)*rdenom
    b(k) = (b1 + g*b2/cpdcv - g*b3 - b4 + g*b5)*rdenom

  END DO
!
!----------------------------------------------------------------------------
!
!  Calculate the constant of integration for f(z) by equating the horizontal
!  average of the recovered pprt to the model's pprt at the first scalar
!  level above the lower boundary (the lowest level at which we have pc);
!  that is, at z = dz/2.  The constant of integration is then given by:
!
!  f0 = horizontal average of (model's pprt - pc)
!
!----------------------------------------------------------------------------
!
  f0 = 0.

  DO i = 2, nx-2
    DO j = 2, ny-2
      f0 = f0 + pprt(i,j,2,tlevel) - pc(i,j,2)
    END DO
  END DO

  f0 = f0*rdenom
!
!----------------------------------------------------------------------------
!
!  Evaluate the analytic solution of the ordinary differential equation
!  df/dz + a(z)f + b(z) = 0.  The general solution is (cf Braun, 1975:
!  Differential Equations and Their Applications):
!
!  f(z) = exp(-integral(a(z))) * [f(0) - integral(b*exp(+integral(a)))]
!
!  The lower limit of the integration is the first scalar grid point
!  above the lower physical boundary; that is, at z = dz/2 (k=2).
!
!----------------------------------------------------------------------------
!
  CALL retrfz(nz,nn,f0,a,b,f,term1,term2)
!
!----------------------------------------------------------------------------
!
!  Recover the perturbation pressure from:  pprt = pc + f(z).
!  Results are stored in tem2.
!
!----------------------------------------------------------------------------
!
!
  9119  CONTINUE
!ccxxx
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        IF (bc_opt /= 0) THEN
          tem2(i,j,k) = pc(i,j,k) + f(k)
        ELSE
          tem2(i,j,k) = pc(i,j,k)
        END IF
      END DO
    END DO
  END DO
!
!----------------------------------------------------------------------------
!  Options for estimate the perturbation pressure at the lowest level :
!
!  ig=0
!  To estimate the perturbation pressure tem2 at k = 1 and k = nz-1, use
!  the first three terms of a Taylor series expansion of tem2 about k=3
!  and k=nz-3, respectively.
!
!  ig=1
!  Use the zero gradient boundary conditions.
!
!----------------------------------------------------------------------------
!ccxxx
!
  DO j=2,ny-2
    DO i=2,nx-2
      IF(ig == 0) THEN
        k=1
        tem2(i,j,k) = 3.0*tem2(i,j,k+1)                                 &
                     -3.0*tem2(i,j,k+2) + tem2(i,j,k+3)
        k=nz-1
        tem2(i,j,k) = 3.0*tem2(i,j,k-1)                                 &
                     -3.0*tem2(i,j,k-2) + tem2(i,j,k-3)
      ELSE
        tem2(i,j,1) = tem2(i,j,2)
        tem2(i,j,nz-1) = tem2(i,j,nz-2)
      END IF
    END DO
  END DO

!
!----------------------------------------------------------------------------
!
!  Recover the perturbation potential temperature, tem1, from the vertical
!  equation of motion as a residual.  When approximating dpdz, use
!  centered differences over 2 dz intervals.
!
!  Note:  since wforce (and hence tem3) is only known for i between 2
!  and nx-2, and for j between 2 and ny-2 we can't recover tem1 on
!  i=1 or i=nx-1 or j=1 or j=ny-1.  To get tem1 on these boundaries, we
!  must impose the model's boundary conditions (but we don't do it in this
!  subroutine).
!
!----------------------------------------------------------------------------
!

  itermax = 0
  pterr = 0.01  ! tolerence in pt iteration

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        dpdz = 0.5*dzinv*(tem2(i,j,k+1) - tem2(i,j,k-1))
        wfave = 0.5*(tem3(i,j,k+1) + tem3(i,j,k))
!
!
!   qc(i,j,k,tlevel) = 0.0
!    qr(i,j,k,tlevel) = qr(i,j,k,tlevel)    ! same
!cc    qr(i,j,k,tlevel) = 0.7 * qr(i,j,k,tlevel)    ! Kelvin's suggestion
!    qr(i,j,k,tlevel) = 0.0                 ! exclude qr effect
!   qi(i,j,k,tlevel) = 0.0
!   qs(i,j,k,tlevel) = 0.0
!   qh(i,j,k,tlevel) = 0.0

!
!
! ==> Original form
!
        tem1(i,j,k) = ptbar(i,j,k) *                                    &
                          ((-wfave + dpdz)/(rhostr(i,j,k)*g) +          &
                            tem2(i,j,k)/(cpdcv*pbar(i,j,k))             &
                          + qc(i,j,k,tlevel) + qr(i,j,k,tlevel)         &
                          - 0.608*(qv(i,j,k,tlevel) - qvbar(i,j,k)) )
!
!
!     tem1(i,j,k) = ptbar(i,j,k) *
!  :              (  (-wfave + dpdz)/(rhostr(i,j,k)*g) +
!  :                  tem2(i,j,k)/(cpdcv*pbar(i,j,k))
!  :                + ( qc(i,j,k,tlevel) + qr(i,j,k,tlevel)
!  :                   - (1.0-rddrv)*(qv(i,j,k,tlevel)-qvbar(i,j,k))
!  :                   /(rddrv+qvbar(i,j,k)) ) / (1.0+qvbar(i,j,k))  )
!
!#### Scenario No. 1: use new scheme.
!
! ==> New scheme for properly hanling inclusion of qv in retrieving
!  ptprt. Assuming RH unchanged before & after the retrival.
!
!  Change to use Newton iteration method. 3-28-98
!
!**** Calculate RH using original p (NOT tem2), pt & qv.
!
!  pzng  =  pbar(i,j,k) +  pprt(i,j,k,tlevel)
!  ptzng = ptbar(i,j,k) + ptprt(i,j,k,tlevel)
!  tzng  = ptzng * (1.0e-5*pzng)**rddcp
!  qvszng = (380./pzng)*exp( 17.27*(tzng-273.15)/(tzng-35.86) )
!  rhzng = qv(i,j,k,tlevel) / qvszng
!
!+++++++++ begin change of subcloud RH
!
!   if(k.ge.2 .and. k.le.5) then     ! change only made from k=2 to 5
!      do k1=2,nz-2
!         if(qr(i,j,k1,tlevel) .gt. 0.5e-3) then  ! qr > 0.5 g/kg
!c               rhzng = rhzng + 0.15         ! increase/decrease by 15%
!            rhzng = rhzng - 0.15         ! increase/decrease by 15%
!            if(rhzng.lt.0.) rhzng=0.
!            if(rhzng.gt.1.) rhzng=1.
!         endif
!      enddo
!   endif
!+++++++++ end of change of subcloud RH
!
!**** Calculate qvs and d(qvs)/d(pt) at ptbar & new p.
!
!  iter = 0
!  tem1(i,j,k) = ptprt(i,j,k,tlevel)   ! initial guess of tem1

!666   continue

!  pzng  =  pbar(i,j,k) + tem2(i,j,k)
!  ptzng = ptbar(i,j,k) + tem1(i,j,k)
!  tzng  = ptzng * (1.0e-5*pzng)**rddcp
!  qvszng = (380./pzng)*exp( 17.27*(tzng-273.15)/(tzng-35.86) )
!  dqvsdpt = qvszng * (1.0e-5*pzng)**rddcp
!    :                 * 4097.9983/(tzng-35.86)**2
!  faczng = 1.0 + 0.608*rhzng*dqvsdpt*ptbar(i,j,k)
!
!**** calculate new ptprt.
!
!  corr = tem1(i,j,k) -
!    :       ptbar(i,j,k) * ( (-wfave + dpdz)/(rhostr(i,j,k)*g) +
!    :                      tem2(i,j,k)/(cpdcv*pbar(i,j,k))
!    :                    + qc(i,j,k,tlevel) + qr(i,j,k,tlevel)
!    :                    - 0.608*(rhzng*qvszng - qvbar(i,j,k)) )
!  corr = corr / faczng

!  if( abs(corr) .ge. pterr ) then
!    iter = iter + 1
!    tem1(i,j,k) = tem1(i,j,k) - corr
!    goto 666
!  endif

!  if( iter .gt. itermax ) itermax = iter
!
!**** Update qv using tem1 & tem2, and initializing qc,qi,qs,qh.
!  Note: BC of retrieved pt' & p' will be applied afterwards
!  in assimptpr.f, but no such an application for qv etc.
!  So we just assume here that boundary values of qv, qc, ...
!  are the same as what they were.
!

!  qv(i,j,k,tlevel) = rhzng * qvszng

!  qc(i,j,k,tlevel) = 0.0
!   qr(i,j,k,tlevel) = qr(i,j,k,tlevel)    ! same
!   qr(i,j,k,tlevel) = 0.0
!  qi(i,j,k,tlevel) = 0.0
!  qs(i,j,k,tlevel) = 0.0
!  qh(i,j,k,tlevel) = 0.0
!
!
!#### Scenario No. 2: totally don't consider qv in getting ptprt.
!
!     tem1(i,j,k) = ptbar(i,j,k) *
!  :                   ( (-wfave + dpdz)/(rhostr(i,j,k)*g) +
!  :                      tem2(i,j,k)/(cpdcv*pbar(i,j,k))
!  :                    + qc(i,j,k,tlevel) + qr(i,j,k,tlevel) )
!

      END DO
    END DO
  END DO

!   write(6,*) ' information about ptprt retrieval via iteration:'
!   write(6,*) ' itermax, pterr=',itermax, pterr
!
!
!----------------------------------------------------------------------------
!  Options for estimate the perturbation temperature at the lowest level :
!
!  ig=0
!  To estimate the perturbation temperature tem1 at k = 1 and k = nz-1, use
!  the first three terms of a Taylor series expansion of tem1 about k=3
!  and k=nz-3, respectively.
!
!  ig=1
!  Use the zero gradient boundary conditions.
!
!----------------------------------------------------------------------------
!
!ccxxx
  DO j=2,ny-2
    DO i=2,nx-2
      IF(ig == 0) THEN
        k=1
        tem1(i,j,k) = 3.0*tem1(i,j,k+1)                                 &
                     -3.0*tem1(i,j,k+2) + tem1(i,j,k+3)
        k=nz-1
        tem1(i,j,k) = 3.0*tem1(i,j,k-1)                                 &
                     -3.0*tem1(i,j,k-2) + tem1(i,j,k-3)
      ELSE
        tem1(i,j,1) = tem1(i,j,2)
        tem1(i,j,nz-1) = tem1(i,j,nz-2)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE retrptpr
