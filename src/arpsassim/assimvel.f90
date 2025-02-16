!
!################################################################
!################################################################
!#####                                                      #####
!#####               SUBROUTINE ASSIMVEL                    #####
!#####                                                      #####
!#####               Copyright (c) 1993                     #####
!#####     Center for Analysis and Prediction of Storms     ######
!#####                University of Oklahoma                ######
!#####                                                      #####
!################################################################
!################################################################
!

SUBROUTINE assimvel(nx,ny,nz,x,y,z,zp,                                  &
           u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,                    &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,j1,j2,j3,                  &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,                &
           tem10,tem11,tem12,tem13)
!
!---------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine solves a 3-D elliptic equation with
!  variable coeficients:
!
!  c1*d2q/dx2 + c2*d2q/dy2 + c3*d2q/dz2 - c4*d2q/dxdy - c5*d2q/dxdz
!  - c6*d2q/dydz - c7*dq/dx - c8*dq/dy - c9*dq/dz = 'Known Forcing'
!
!  subject to the Dirichlet boundary condition:
!
!  q = 0 on the x & y boundaries, and
!
!  subject to the z-boundary condition:
!
!  dq/dz = (z/x**2+y**2)*[2r(vro - vrm) + (xdq/dx + ydq/dy)]
!  where:    vro - observed radial velocity
!            vrm - model derived radial velocity
!
!  The solution is obtained by applying a direct SOR technique
!  (see Haltiner and Williams "Numerical Prediction and Dynamical
!  Weather Prediction", Ch. 5 pg 159.)
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
!  Limin Zhao (10/18/95)
!  Modified the code to obtain the 'lamda2' and 'lamda1' simultaneously
!  by solving the two related equations interatively. It helps to
!  eliminate the sensitivity of the solution to the different
!  numerical schemes.
!
!  Alan Shapiro and Limin Zhao (10/19/95)
!  Put 'lamda1' in new staggered points and modified the related code.
!  The above two modifications work together to ensure the divergence
!  of the wind field be small enough after the variational adjustment.
!
!  Limin Zhao and Alan Shapiro(11/10/95)
!  Put a 'patch' on the top boundary to avoid the singularity problem.
!  Inside the 'patch', the zero gradient B.C. is used to replace the
!  non-zero gradient B.C.
!
!  Limin Zhao (11/16/95)
!  Set an option for scaling 'lamda1'. When the coefficient 'scoef'
!  is set to 1, no scaling is performed.
!
!  Limin Zhao (1/11/96)
!  Strip out the hole-filling and time interpolation part, and make
!  the program mainly to the velocity adjustment, either use Vr or the
!  variational adjustment.
!
!  Alan Shapiro and Limin Zhao (03/25/96)
!  The boundary values in the patch are calculated linearly from
!  its latent boundary values.
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
!    tem1     Observed radial velocity
!    tem5     Model or retrieved radial velocity
!
!---------------------------------------------------------------------
!
!  OUTPUT:
!
!    tem7     Solution of the Poisson equation (q above).
!
!---------------------------------------------------------------------
!
!  WORK ARRAYS:
!
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
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

  INTEGER :: in,jn,ir0,jr0,ipe,ips,jps,jpe,ipatch
  REAL :: acof,bcof

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
  REAL :: tem13(nx,ny,nz)

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

  REAL :: alpha        ! Over-relaxation coefficient
  REAL :: tol          ! An error tolerance, proportional to grid spacing
                       ! and max velocity error
  REAL :: veltol       ! Maximum tolerable error in velocity
  REAL :: newit        ! New iterate of q at a point
  REAL :: stor         ! Stores the new value of q temporarily
  REAL :: resid        ! the residual of the Poisson equation
  REAL :: diff         ! Difference between the updated q and old q at a point
  REAL :: test         ! Value of the maximum difference (qnew-qold) in grid
  REAL :: ptol         ! Error tolerance for hole-filling
  REAL :: anorm
  REAL :: anormf
  REAL :: test1
  REAL :: iterr

  REAL :: b1,b2,b3,b4  ! Coefficients of the Poisson equation.
  REAL :: b5,b6,b7,b8
  REAL :: b9
  REAL :: c1,c2,c3,c4  ! Coefficient combinations of b1...b9 in the Poisson
  REAL :: c5,c6,c7,c8  ! equation. The coefficients are a function of grid
  REAL :: c9,c10       ! spacing and location.
  REAL :: tm1,tm2,tm3
  REAL :: tm4,tm5,tm6
  REAL :: tm7,tm8,tm9
  REAL :: scoef
!ccxxx   You might need change veltol sometime
!   parameter(itmax=100000,veltol=0.01,scoef=1.0)  ! original
!   parameter(itmax=900000,veltol=0.01,scoef=1.0)  ! original
  PARAMETER(itmax=900000, veltol=0.5,scoef=1.0)
!   parameter(itmax=2, veltol=0.05,scoef=1.0)  ! only for test
  REAL :: rad          ! Radius of sphere in Cartesian coord
  REAL :: radh         ! Radius of circle in hz. plane
  REAL :: radinv       ! Inverse radius magnitude

  REAL :: rdx2         ! Reciprocal of (dx)**2
  REAL :: rdy2         ! Reciprocal of (dy)**2
  REAL :: rdz2         ! Reciprocal of (dz)**2

  REAL :: term1,term2  ! Analytic terms representing the RHS of the
  REAL :: term3,term4  ! Poisson eqn.
  REAL :: term5

  REAL :: vphi,vtheta  ! Model tangential and polar velocities
  REAL :: rcos         ! Squareroot of x**2+y**2

  INTEGER :: newebc    ! East boundary condition parameter.
  INTEGER :: newwbc    ! West boundary condition parameter.
  INTEGER :: newnbc    ! North boundary condition parameter.
  INTEGER :: newsbc    ! South boundary condition parameter.

  REAL :: vro,vrm   ! Used for diagnostics/testing
  REAL :: vradj
  REAL :: fpat      ! temp variable for forcing "patch"
  REAL :: maxdiv    ! Maximum divergence (diagnostic calculation)
  INTEGER :: imx,jmx   ! Gridpoint location of max divp
  INTEGER :: kmx
  INTEGER :: imx1,jmx1,kmx1
  REAL :: maxw
  REAL :: rrr       ! used for (x**2+y**2+z**2)**0.5
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL aamult
  EXTERNAL assimrd
  EXTERNAL retrint
  EXTERNAL avgx
  EXTERNAL avgy
  EXTERNAL avgz
  EXTERNAL dtadump
  EXTERNAL gtdmpfn
  EXTERNAL bcu
  EXTERNAL bcv
  EXTERNAL rhouvw
  EXTERNAL lbcw
  EXTERNAL vbcw
  EXTERNAL smth
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
!
!  The following is a brief summary of the assimilation flags
!  available within this routine.  For a more detailed review
!  of the options see ASSIMCON.F.
!
!  ivar    Variational adjustment flag:
!
!          = 0, NO variational adjustment this model time step
!          = 1, Perform variational adjustment this model time step
!
!  insrt   Insertion flag:
!
!          = 0, Do NOT insert velocities this time step
!          = 1, Insert velocities this time step
!
!  A matrix describing the assimilation options herein:
!
!       ivar   insrt   Comments
!       --------------------------------------------------------
!        0       0     Nothing - Exit assimilation mode
!        0       1     Direct insertion only
!        1       0     Variational adjustment only
!        1       1     Variational adjustment + direct insertion
!
!  NOTE: If irecov=1, then the ingest and time interpolation of
!       3 input data files occurs herein.
!
!-----------------------------------------------------------------------
!
  WRITE(6,*)'code in assimvel'

  IF ( insrt == 0.AND.ivar == 0.) RETURN
!
!-----------------------------------------------------------------------
!
!  Has switch been set for a direct insertion only? If yes, ingest
!  data.
!
!  NOTE: Insertion of data is done here ONLY if there is no
!        variational adjustment (ivopt=0). If ivopt=1, then the
!        adjusted data is inserted at the end of this routine.
!
!-----------------------------------------------------------------------
!
  tim = tpresent

  umove=0.0
  vmove=0.0
  xmove= xshift - umove*(curtim-assimtim(1))
  ymove= yshift - vmove*(curtim-assimtim(1))
!
!-------------------------------------------------------------------------
!
!  Interpolate input 'model' velocity components to scalar grid point.
!  Where:
!              tem2 = scalar u component
!              tem3 = scalar v component
!              tem4 = scalar w component
!
!-------------------------------------------------------------------------
!
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1,nx
        tem7(i,j,k) = u(i,j,k,tim)
        tem8(i,j,k) = v(i,j,k,tim)
        tem9(i,j,k) = w(i,j,k,tim)
      END DO
    END DO
  END DO

!

  CALL avgx(tem7, 0, nx,ny,nz,                                          &
            1,nx-1,1,ny-1,1,nz-1,tem2)
  CALL avgy(tem8, 0, nx,ny,nz,                                          &
              1,nx-1,1,ny-1,1,nz-1,tem3)
  CALL avgz(tem9, 0, nx,ny,nz,                                          &
              1,nx-1,1,ny-1,1,nz-1,tem4)

!
!----------------------------------------------------------------------
!
!    Calculate the first guess radial velocity component using
!    the hole-filled and blended velocity data.
!
!----------------------------------------------------------------------
!

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1


        xs = 0.5*(x(i)+x(i+1)) - xmove
        ys = 0.5*(y(j)+y(j+1)) - ymove
        zs = 0.5*(z(k)+z(k+1)) - zshift

        rad         = SQRT(xs**2+ys**2+zs**2)
        tem5(i,j,k) = (tem2(i,j,k)*xs+tem3(i,j,k)*ys+                   &
                        tem4(i,j,k)*zs)/rad

      END DO
    END DO
  END DO
!
!----------------------------------------------------------------------
!
!  Save rhobar in tem9
!
!----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem9(i,j,k) = rhostr(i,j,k)/j3(i,j,k)
      END DO
    END DO
  END DO
!
!----------------------------------------------------------------------
!
!    Check whether the variational adjusment is needed.
!
!----------------------------------------------------------------------
!
  IF ( insrt == 1.AND.ivar == 0) THEN

    WRITE(6,'(/5x,a)') 'Performing a direct insertion only'
!
!-----------------------------------------------------------------------
!
!  Calculate the model tangential velocity, vphi, and model polar
!  velocity, vtheta, at a scalar point.  The adjusted winds ua,va,
!  and wa, are constructed (at a scalar point) as follows:
!
!     tem6 = ua = Vr*x/rad + vphi*y/rcos - vtheta*x*z/rad*rcos
!     tem7 = va = Vr*y/rad - vphi*x/rcos - vtheta*y*z/rad*rcos
!     tem8 = wa = Vr*z/rad + vtheta*rcos/rad
!
!  Where,
!     rcos = sqrt(x**2+y**2)
!     rad  = sqrt(x**2+y**2+z**2)
!     tem1 = Vr (model (OSSE) or observed radial velocity)
!
!  For more details on this insertion technique, see the write-up
!  "Direct insertion of single-Doppler radial winds into a numerical
!  model" by Alan Shapiro.
!
!-----------------------------------------------------------------------
!

    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1

          xs    = 0.5*(x(i)+x(i+1)) - xmove
          ys    = 0.5*(y(j)+y(j+1)) - ymove
          zs    = 0.5*(z(k)+z(k+1)) - zshift

          rad   = SQRT(xs**2+ys**2+zs**2)

          rcos  = SQRT(xs**2+ys**2)

          vphi  = (1./rcos)*(tem2(i,j,k)*ys - tem3(i,j,k)*xs)

          vtheta= (1./(rad*rcos))*(-tem2(i,j,k)*xs*zs                   &
                                   -tem3(i,j,k)*ys*zs                   &
                                   +tem4(i,j,k)*(rcos**2))

!-----------------------------------------------------
!        Comment out overwrite of obs radial component
!        01-03-99 SSW for AMS99
!
!        tem6(i,j,k) = tem1(i,j,k)*xs/rad + vphi*ys/rcos
!    :                  - vtheta*xs*zs/(rad*rcos)
!        tem7(i,j,k) = tem1(i,j,k)*ys/rad - vphi*xs/rcos
!    :                  - vtheta*ys*zs/(rad*rcos)
!        tem8(i,j,k) = tem1(i,j,k)*zs/rad + vtheta*rcos/rad
!
          tem6(i,j,k) = tem2(i,j,k)
          tem7(i,j,k) = tem3(i,j,k)
          tem8(i,j,k) = tem4(i,j,k)
!-----------------------------------------------------
!
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!    Reconstruct the model velocity field on the Arakawa C-grid. u,v,
!    and w are stored in tem2,tem3, and tem4 respectively until the
!    boundary conditions are applied at the end of this routine.
!
!-----------------------------------------------------------------------
!

    DO k = 2, nz-1
      DO j = 2, ny-1
        DO i = 2, nx-1
          tem2(i,j,k) = .5*(tem6(i-1,j,k)+tem6(i,j,k))
          tem3(i,j,k) = .5*(tem7(i,j-1,k)+tem7(i,j,k))
          tem4(i,j,k) = .5*(tem8(i,j,k-1)+tem8(i,j,k))
        END DO
      END DO
    END DO

    GO TO 270                          ! Apply boundary conditions
!
!----------------------------------------------------------------------
!
!  Perform a variational adjustment
!
!----------------------------------------------------------------------
!
  ELSE

    PRINT *,'Performing a variational adjustment'

!
!----------------------------------------------------------------------
!
!  Redefine the following arrays such that:
!
!            tem1 = tem1*rhostr   tem4 = tem4*rhostr
!            tem2 = tem2*rhostr   tem5 = tem5*rhostr
!            tem3 = tem3*rhostr
!
!----------------------------------------------------------------------
!
    CALL aamult(tem1,tem9,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1)

    CALL avgsu(tem9,nx,ny,nz, 1,ny-1, 1,nz-1, tem6, tem10)
    CALL aamult(u(1,1,1,tim),tem6,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem2)

    CALL avgsv(tem9,nx,ny,nz, 1,nx-1, 1,nz-1, tem6, tem10)
    CALL aamult(v(1,1,1,tim),tem6,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem3)

    CALL avgsw(tem9,nx,ny,nz, 1,nx-1, 1,ny-1, tem6)
    CALL aamult(w(1,1,1,tim),tem6,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem4)

    CALL aamult(tem5,tem9,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem5)

    CALL avgsw(tem9,nx,ny,nz, 1,nx-1, 1,ny-1, tem12)  ! store rhobar at
                                                      ! w points.

    WRITE(6,*) 'after redefine'
!
!----------------------------------------------------------------------
!
!  Set the optimum over-relaxation coefficient, alpha.
!  WARNING: In order to avoid solution 'blow-up', alpha should be
!           near 1.
!
!----------------------------------------------------------------------
!
!     alpha=0.25    ! original value; under relaxation, slow convergence
!    alpha=1.00    ! under relaxation, slow convergence
    alpha=0.75    ! under relaxation, slow convergence
!     alpha=0.15    ! under relaxation, slow convergence
!     alpha=1.5    ! SOR; optimum -> 2 for large dimensions.

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          tem6(i,j,k)  = 0.0
          tem7(i,j,k)  = 0.0
          tem8(i,j,k)  = 0.0
          tem10(i,j,k) = 0.0
          tem11(i,j,k) = 0.0
        END DO
      END DO
    END DO

    dxinv     = 1./dx
    dyinv     = 1./dy
    dzinv     = 1./dz

    rdx2 = dxinv*dxinv        ! reciprocal of (dx)**2
    rdy2 = dyinv*dyinv        ! reciprocal of (dy)**2
    rdz2 = dzinv*dzinv        ! reciprocal of (dz)**2
!
!----------------------------------------------------------------------
!
!    Start the interation loop.
!
!----------------------------------------------------------------------
!
    WRITE(6,*) 'start iteration'

    DO it = 1,itmax
!
!----------------------------------------------------------------------
!
!  Compute the right hand side of the Poisson equation and store it
!  in tem6.  This will be handled externally once the code is
!  completed and tested.
!
!----------------------------------------------------------------------
!
      DO k = 2, nz-2
        DO j = 2, ny-2
          DO i = 2, nx-2
!
!
!       xs1       = x(i+1) - xmove
!       ys1       = y(j+1) - ymove
!       zs1       = z(k+1) - zshift
!
!       xs       = x(i) - xmove
!       ys       = y(j) - ymove
!       zs       = z(k) - zshift
!
!       term1 = dxinv*(tem2(i+1,j,k)-tem2(i,j,k))
!  :          + dyinv*(tem3(i,j+1,k)-tem3(i,j,k))
!  :          + dzinv*(tem4(i,j,k+1)-tem4(i,j,k))
!
!
!       term2 = dxinv*(xs1*(tem8(i+1,j,k)+tem8(i+1,j+1,k)
!  :          +             tem8(i+1,j,k+1)+tem8(i+1,j+1,k+1))
!  :          -         xs*(tem8(i,j,k)+tem8(i,j+1,k)
!  :          +             tem8(i,j,k+1)+tem8(i,j+1,k+1)))
!  :          + dyinv*(ys1*(tem8(i,j+1,k)+tem8(i+1,j+1,k)
!  :          +             tem8(i,j+1,k+1)+tem8(i+1,j+1,k+1))
!  :          -         ys*(tem8(i,j,k)+tem8(i+1,j,k)
!  :          +             tem8(i,j,k+1)+tem8(i+1,j,k+1)))
!  :          + dzinv*(zs1*(tem8(i,j,k+1)+tem8(i+1,j,k+1)
!  :          +             tem8(i,j+1,k+1)+tem8(i+1,j+1,k+1))
!  :          -         zs*(tem8(i,j,k)+tem8(i+1,j,k)
!  :          +             tem8(i,j+1,k)+tem8(i+1,j+1,k)))
!
!       term2       = term2/scoef
!
!       tem6(i,j,k) = -2.0*term1 + 0.25*term2
!
! **** new formulation
!
            xs       = 0.5 * ( x(i) + x(i+1) ) - xmove
            ys       = 0.5 * ( y(j) + y(j+1) ) - ymove
            zs       = 0.5 * ( z(k) + z(k+1) ) - zshift
            rrr      = SQRT(xs*xs + ys*ys + zs*zs)

            term1 = dxinv*(tem2(i+1,j,k)-tem2(i,j,k))                   &
                  + dyinv*(tem3(i,j+1,k)-tem3(i,j,k))                   &
                  + dzinv*(tem4(i,j,k+1)-tem4(i,j,k))

            term2 = xs*dxinv*( (tem8(i+1,j,k)+tem8(i+1,j+1,k)           &
                  +             tem8(i+1,j,k+1)+tem8(i+1,j+1,k+1))      &
                  -            (tem8(i,j,k)+tem8(i,j+1,k)               &
                  +             tem8(i,j,k+1)+tem8(i,j+1,k+1)) )        &
                           + ys*dyinv*( (tem8(i,j+1,k)+tem8(i+1,j+1,k)  &
                           +             tem8(i,j+1,k+1)+tem8(i+1,j+1,k+1)) &
                           -            (tem8(i,j,k)+tem8(i+1,j,k)      &
                  +             tem8(i,j,k+1)+tem8(i+1,j,k+1)) )        &
                           + zs*dzinv*( (tem8(i,j,k+1)+tem8(i+1,j,k+1)  &
                           +             tem8(i,j+1,k+1)+tem8(i+1,j+1,k+1)) &
                           -            (tem8(i,j,k)+tem8(i+1,j,k)      &
                           +             tem8(i,j+1,k)+tem8(i+1,j+1,k)) )

            term2 = 0.25 * term2 + 2.0 * 0.125 *                        &
                     ( tem8(i,  j,  k  )+tem8(i+1,j,  k  ) +            &
                       tem8(i,  j+1,k  )+tem8(i+1,j+1,k  ) +            &
                       tem8(i,  j,  k+1)+tem8(i+1,j,  k+1) +            &
                       tem8(i,  j+1,k+1)+tem8(i+1,j+1,k+1) )

            term2 = term2/rrr

            tem6(i,j,k) = - 2.0 * term1 + term2

          END DO
        END DO
      END DO
!
!----------------------------------------------------------------------
!
!  Solve the elliptic equation via the successive overrelaxation
!  method.  Check for convergence. The following coefficients
!  are defined as:
!
!  b1 = rdx2
!  b2 = rdy2
!  b3 = rdz2
!
!----------------------------------------------------------------------
!
      b1          = rdx2
      b2          = rdy2
      b3          = rdz2
      c1          = -2.*(b1 + b2 + b3)
      c2          = b1
      c3          = b1
      c4          = b2
      c5          = b2
      c6          = b3
      c7          = b3

      test   = 0.0
      DO k = 2, nz-2
        DO j = 2, ny-2
          DO i = 2, nx-2

            tm1         = c2*tem7(i+1,j,k)
            tm2         = c3*tem7(i-1,j,k)
            tm3         = c4*tem7(i,j+1,k)
            tm4         = c5*tem7(i,j-1,k)
            tm5         = c6*tem7(i,j,k+1)
            tm6         = c7*tem7(i,j,k-1)

            resid  = c1*tem7(i,j,k) + tm1 + tm2 + tm3 + tm4             &
                                    + tm5 + tm6 - tem6(i,j,k)

            newit  = tem7(i,j,k) - alpha*resid/c1   ! SOR method

!       diff   = abs(newit      - tem7(i,j,k))
!       diff = abs(newit - tem7(i,j,k)) / tem9(i,j,k)

            diff = ABS(resid/c1) / tem9(i,j,k)   ! alpha excluded;
                                                 ! normakized residual

            IF(diff > test) THEN
              test = diff
              imx = i
              jmx = j
              kmx = k
            END IF

!
!       tem7(i,j,k) = newit

            tem11(i,j,k) = newit

            tem10(i,j,k) = resid          !outputs for test

          END DO
        END DO
      END DO

!      write(6,*) 'after relaxation'

!ccxxx
!ccxxx Test simutenous relaxation
!ccxxx
      DO k = 2, nz-2
        DO j = 2, ny-2
          DO i = 2, nx-2

            tem7(i,j,k) = tem11(i,j,k)

          END DO
        END DO
      END DO
!
!----------------------------------------------------------------------
!
!  Set the boundary conditions on Xl, Xr, Yl, & Yr
!
!----------------------------------------------------------------------
!
      DO k = 2, nz-2          !XL,XR B.C.
        DO j = 2, ny-2
          tem7(1,j,k)    = - tem7(2,j,k)
          tem7(nx-1,j,k) = - tem7(nx-2,j,k)
        END DO
      END DO

      DO k = 2,nz-2          !YL,YR B.C
        DO i = 2,nx-2
          tem7(i,1,k)    = - tem7(i,2,k)
          tem7(i,ny-1,k) = - tem7(i,ny-2,k)
        END DO
      END DO
!
!----------------------------------------------------------------------
!
!  Set the top and bottom boundary conditions.
!
!----------------------------------------------------------------------
!
      z2   = z(2)    - zshift
      znz1 = z(nz-1) - zshift
      DO i=2,nx-2
        DO j=2,ny-2

!       xs = x(i) - xmove
!       ys = y(j) - ymove
!
!       radh = sqrt(xs*xs + ys*ys)
!
!       tem7(i,j,1)    = tem7(i,j,2)          !'lamda2' grad B.C.
!  : -                   0.25*dz*z2*(tem8(i,j,2)+tem8(i,j+1,2)
!  : +                   tem8(i+1,j,2)+tem8(i+1,j+1,2))/scoef
!
!       tem7(i,j,nz-1) = tem7(i,j,nz-2)
!  : +                   0.25*dz*znz1*(tem8(i,j,nz-1)+tem8(i,j+1,nz-1)
!  : +                   tem8(i+1,j,nz-1)+tem8(i+1,j+1,nz-1))/scoef
!
          xs = 0.5 * ( x(i) + x(i+1) ) - xmove
          ys = 0.5 * ( y(j) + y(j+1) ) - ymove

          rrr = SQRT(xs*xs + ys*ys + z2*z2)

          tem7(i,j,1)    = tem7(i,j,2)          & !'lamda2' grad B.C.
          -                   0.25*dz*z2*(tem8(i,j,2)+tem8(i,j+1,2)     &
              +                   tem8(i+1,j,2)+tem8(i+1,j+1,2)) / rrr  &
              +                   2.0 * dz * tem4(i,j,2)

          rrr = SQRT(xs*xs + ys*ys + znz1*znz1)

          tem7(i,j,nz-1) = tem7(i,j,nz-2)                               &
              +                   0.25*dz*znz1*(tem8(i,j,nz-1)+tem8(i,j+1,nz-1) &
              +                   tem8(i+1,j,nz-1)+tem8(i+1,j+1,nz-1)) / rrr &
              -                   2.0 * dz * tem4(i,j,nz-1)


        END DO
      END DO

!ccxxx
!ccxxx  Patch: you might like to change it.
!ccxxx

      ir0 = xmove/dx + 1
      jr0 = ymove/dy + 1
      ipatch =5
      ips=ir0-ipatch
      ipe=ir0+ipatch
      jps=jr0-ipatch
      jpe=jr0+ipatch

      IF((ir0-ipatch) < 2) ips=2
      IF((jr0-ipatch) < 2) jps=2
      IF((ir0+ipatch) > (nx-2)) ipe=nx-2
      IF((jr0+ipatch) > (ny-2)) jpe=ny-2

      DO in = ips,ipe
        DO jn = jps,jpe


          acof = (jn-jps)/(jpe-jps)
          tem7(in,jn,nz-1) = (tem7(in,jpe,nz-1)-tem7(in,jps,nz-1))*acof &
                            +  tem7(in,jps,nz-1)

          acof = (in-ips)/(ipe-ips)
          tem7(in,jn,nz-1) = 0.5*(tem7(in,jn,nz-1) + ((tem7(ipe,jn,nz-1) &
                         - tem7(ips,jn,nz-1))*acof) + tem7(ips,jn,nz-1))

!    acof = (in-ips)/(ipe-ips)
!    bcof = (jn-jps)/(jpe-jps)
!    tem7(in,jn,nz-1) = (1.0-acof)*(1.0-bcof)*tem7(ips,jps,nz-1)
!  :                  + (1.0-acof)*     bcof *tem7(ips,jpe,nz-1)
!  :                  +      acof *     bcof *tem7(ipe,jpe,nz-1)
!  :                  +      acof *(1.0-bcof)*tem7(ipe,jps,nz-1)
!
!    tem7(in,jn,nz-1) = 0.5 * (tem7(in,jn,nz-1) + tem7(in,jn,nz-2))
!
!    tem7(in,jn,nz-1) = tem7(in,jn,nz-2)     ! zeor gradient

        END DO
      END DO


!   ir0 = xmove/dx + 1
!   jr0 = ymove/dy + 1
!   ipatch = 10             ! within a radius of 20 km
!   ips=ir0-ipatch
!   ipe=ir0+ipatch
!   jps=jr0-ipatch
!   jpe=jr0+ipatch

!   if((ir0-ipatch).lt.2) ips=2
!   if((jr0-ipatch).lt.2) jps=2
!   if((ir0+ipatch).gt.(nx-2)) ipe=nx-2
!   if((jr0+ipatch).gt.(ny-2)) jpe=ny-2
!
!   DO in = ips,ipe
!   DO jn = jps,jpe
!     tem7(in,jn,1) = 0.5 * ( tem7(in,jn,1) +
!    :          0.25*(tem7(in-1,jn,1)+tem7(in+1,jn,1)+
!    :                tem7(in,jn-1,1)+tem7(in,jn+1,1)) )
!     tem7(in,jn,nz-1) = 0.5 * ( tem7(in,jn,nz-1) +
!    :          0.25*(tem7(in-1,jn,nz-1)+tem7(in+1,jn,nz-1)+
!    :                tem7(in,jn-1,nz-1)+tem7(in,jn+1,nz-1)) )
!   END DO
!   END DO

!   DO  k = 2, nz-2
!   DO in = ips,ipe
!   DO jn = jps,jpe
!     tem7(in,jn,k) = 0.5 * ( tem7(in,jn,k) +
!    :               (tem7(in-1,jn,k)+tem7(in+1,jn,k)+
!    :                tem7(in,jn-1,k)+tem7(in,jn+1,k)+
!    :                tem7(in,jn,k-1)+tem7(in,jn,k+1))/6.0 )
!   END DO
!   END DO
!   END DO

      DO i=2,nx-2
        tem7(i,1,1)       = - tem7(i,2,1)
        tem7(i,1,nz-1)    = - tem7(i,2,nz-1)
        tem7(i,ny-1,1)    = - tem7(i,ny-2,1)
        tem7(i,ny-1,nz-1) = - tem7(i,ny-2,nz-1)
      END DO

      DO j=2,ny-2
        tem7(1,j,1)       = - tem7(2,j,1)
        tem7(1,j,nz-1)    = - tem7(2,j,nz-1)
        tem7(nx-1,j,1)    = - tem7(nx-2,j,1)
        tem7(nx-1,j,nz-1) = - tem7(nx-2,j,nz-1)
      END DO
!
!----------------------------------------------------------------------
!
!    Set the corner points.
!
!----------------------------------------------------------------------
!
      DO k = 1,nz-1
        tem7(1,1,k)       = - tem7(2,2,k)
        tem7(1,ny-1,k)    = - tem7(2,ny-2,k)
        tem7(nx-1,1,k)    = - tem7(nx-2,2,k)
        tem7(nx-1,ny-1,k) = - tem7(nx-2,ny-2,k)
      END DO
!
!----------------------------------------------------------------------
!
!   Determine the first Lagrange Multiplier, Lamda1. Lamda1 is
!   assumed valid at integer points(i,j,k). The result is
!   stored in tem8.
!
!----------------------------------------------------------------------
!
      test1=0.0
      DO k = 2, nz-1
        DO j = 2, ny-1
          DO i = 2, nx-1

            xs       = x(i) - xmove
            ys       = y(j) - ymove
            zs       = z(k) - zshift

            rad        = SQRT(xs**2+ys**2+zs**2)
            radinv     = 1./(rad**2)


            IF( rad < 1.0 ) THEN
              rad = 1.0
              radinv     = 1./(rad**2)
            END IF

            term1 = dxinv*xs*((tem7(i,j-1,k-1)-tem7(i-1,j-1,k-1))       &
                +                    (tem7(i,j-1,k)-tem7(i-1,j-1,k))    &
                +                    (tem7(i,j,k-1)-tem7(i-1,j,k-1))    &
                +                    (tem7(i,j,k)-tem7(i-1,j,k)))

            term2 = dyinv*ys*((tem7(i-1,j,k-1)-tem7(i-1,j-1,k-1))       &
                +                    (tem7(i,j,k-1)-tem7(i,j-1,k-1))    &
                +                    (tem7(i-1,j,k)-tem7(i-1,j-1,k))    &
                +                    (tem7(i,j,k)-tem7(i,j-1,k)))

            term3 = dzinv*zs*((tem7(i,j-1,k)-tem7(i,j-1,k-1))           &
                +                    (tem7(i-1,j-1,k)-tem7(i-1,j-1,k-1)) &
                +                    (tem7(i-1,j,k)-tem7(i-1,j,k-1))    &
                +                    (tem7(i,j,k)-tem7(i,j,k-1)))

            term4 = 0.125*((tem5(i,j,k)+tem5(i-1,j,k)                   &
                +                  tem5(i-1,j-1,k)+tem5(i,j-1,k)        &
                +                  tem5(i,j,k-1)+tem5(i-1,j,k-1)        &
                +                  tem5(i-1,j-1,k-1)+tem5(i,j-1,k-1))   &
                -                 (tem1(i,j,k)+tem1(i-1,j,k)            &
                +                  tem1(i-1,j-1,k)+tem1(i,j-1,k)        &
                +                  tem1(i,j,k-1)+tem1(i-1,j,k-1)        &
                +                  tem1(i-1,j-1,k-1)+tem1(i,j-1,k-1)))

!        term5  = 2.0*term4/rad + 0.25*radinv*(term1+term2+term3)

            term5  = 2.0*term4 + 0.25*(term1+term2+term3)/rad

            diff = ABS(term5-tem8(i,j,k)) / tem12(i,j,k)
            IF(diff > test1) THEN
              test1 = diff
              imx1 = i
              jmx1 = j
              kmx1 = k
            END IF

            tem8(i,j,k) = scoef*term5

          END DO
        END DO
      END DO

      DO k = 2, nz-1
        DO j = 2, ny-1
          DO i = 2, nx-1
            xs = x(i) - xmove
            ys = y(j) - ymove
            zs = z(k) - zshift
            rad = SQRT(xs**2+ys**2+zs**2)

            IF( rad < 1.0 ) THEN
              tem8(i,j,k) = 0.0
              radinv = 0.0

              IF(i > 2) THEN
                tem8(i,j,k) = tem8(i,j,k) + tem8(i-1,j,k)
                radinv = radinv + 1.
              END IF

              IF(i < nx-1) THEN
                tem8(i,j,k) = tem8(i,j,k) + tem8(i+1,j,k)
                radinv = radinv + 1.
              END IF

              IF(j > 2) THEN
                tem8(i,j,k) = tem8(i,j,k) +  tem8(i,j-1,k)
                radinv = radinv + 1.
              END IF

              IF(j < ny-1) THEN
                tem8(i,j,k) = tem8(i,j,k) +  tem8(i,j+1,k)
                radinv = radinv + 1.
              END IF

              IF(k > 2) THEN
                tem8(i,j,k) = tem8(i,j,k) +  tem8(i,j,k-1)
                radinv = radinv + 1.
              END IF

              IF(k < nz-1) THEN
                tem8(i,j,k) = tem8(i,j,k) +  tem8(i,j,k+1)
                radinv = radinv + 1.
              END IF

              tem8(i,j,k) = tem8(i,j,k) / radinv

              WRITE(6,*) 'Do patching for lamda1 at i,j,k=',i,j,k
            END IF

          END DO
        END DO
      END DO

!
!----------------------------------------------------------------------
!
!  Test for solution convergence.
!
!----------------------------------------------------------------------
!
!     tol = sqrt(dx**2+dy**2+dz**2)*veltol
      tol = MIN( MIN(dx, dy), dz) * veltol

!     ix = mod(it,1)
      ix = MOD(it,50)
      IF(it == 50) THEN
        WRITE(6,*) ' Max diff at iteration: imx,jmx,kmx,test,test1'
      END IF

      IF(ix == 0) THEN
        WRITE(6,'(i9,2x,2i4,i3,e16.6,2x,2i4,i3,e16.6)')                 &
            it,imx,jmx,kmx,test, imx1,jmx1,kmx1,test1
      END IF

!     IF(test.le.tol.and.it.gt.2) THEN
      IF(test <= tol .AND. test1 <= veltol .AND. it > 2) THEN

        WRITE(6,'(/5x,a,i6,/5x,a,2f12.6)')                              &
            'The SOR technique converged. The # of iterations was ',    &
            it,' The tolerance level was set at: ',veltol,tol
        WRITE(6,'(i9,2x,2i4,i3,e16.6,2x,2i4,i3,e16.6)')                 &
            it,imx,jmx,kmx,test, imx1,jmx1,kmx1,test1

        GO TO 258

      END IF

    END DO

    STOP


    258     CONTINUE

    WRITE(6,*) 'patch: xmove,ymove,ir0,jr0,ipatch=',                    &
                       xmove,ymove,ir0,jr0,ipatch
    WRITE(6,*) 'patch: ips,ipe,jps,jpe=',ips,ipe,jps,jpe

    2201    CONTINUE

!
!----------------------------------------------------------------------
!
!   Determine the adjusted velocity, ua, store results in tem2:
!
!           2*rhobar*(Ua-Um) + lamda1*x - dlamda2/dx = 0
!
!   NOTE:  rhobar (tem9) is assumed horizontally homogenous.
!
!----------------------------------------------------------------------
!
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          tem9(i,j,k)= rhostr(i,j,k)/j3(i,j,k)
        END DO
      END DO
    END DO

    CALL avgsu(tem9,nx,ny,nz, 1,ny-1, 1,nz-1, tem6, tem10)

    DO k = 2, nz-2
      DO j = 2, ny-2
        DO i = 2, nx-1

!       xs       = x(i) - xmove
!
!       tem2(i,j,k)  =  u(i,j,k,tim) + (.5/tem6(i,j,k))
!  :                    *(dxinv*(tem7(i,j,k)-tem7(i-1,j,k))
!  :                 -  0.25*xs*(tem8(i,j,k)+tem8(i,j+1,k)
!  :                 +  tem8(i,j,k+1)+tem8(i,j+1,k+1)))
!
          xs  =         x(i)            - xmove
          ys  = 0.5 * ( y(j) + y(j+1) ) - ymove
          zs  = 0.5 * ( z(k) + z(k+1) ) - zshift
          rrr = SQRT(xs*xs + ys*ys + zs*zs)

          tem2(i,j,k)  =  u(i,j,k,tim) + (.5/tem6(i,j,k))               &
                       * ( dxinv*(tem7(i,j,k)-tem7(i-1,j,k))            &
                         - 0.25*(xs/rrr)*(tem8(i,j,k)+tem8(i,j+1,k)     &
                           +tem8(i,j,k+1)+tem8(i,j+1,k+1)) )

        END DO
      END DO
    END DO

!
!----------------------------------------------------------------------
!
!   Determine the adjusted velocity, va, store results in tem3:
!
!           2*rhobar*(Va-Vm) + lamda1*y - dlamda2/dy = 0
!
!   NOTE:  rhobar (tem9) is assumed horizontally homogenous.
!
!----------------------------------------------------------------------
!
    CALL avgsv(tem9,nx,ny,nz, 1,nx-1, 1,nz-1, tem6, tem10)

    DO k = 2, nz-2
      DO j = 2, ny-1
        DO i = 2, nx-2

!        ys       = y(j) - ymove
!
!        tem3(i,j,k)  =  v(i,j,k,tim) + (.5/tem6(i,j,k))
!  :                     *(dyinv*(tem7(i,j,k)-tem7(i,j-1,k))
!  :                  -  0.25*ys*(tem8(i,j,k)+tem8(i+1,j,k)
!  :                  +  tem8(i,j,k+1)+tem8(i+1,j,k+1)))

          xs  = 0.5 * ( x(i) + x(i+1) ) - xmove
          ys  =         y(j)            - ymove
          zs  = 0.5 * ( z(k) + z(k+1) ) - zshift
          rrr = SQRT(xs*xs + ys*ys + zs*zs)

          tem3(i,j,k)  =  v(i,j,k,tim) + (.5/tem6(i,j,k))               &
                       * ( dyinv*(tem7(i,j,k)-tem7(i,j-1,k))            &
                         - 0.25*(ys/rrr)*(tem8(i,j,k)+tem8(i+1,j,k)     &
                           +tem8(i,j,k+1)+tem8(i+1,j,k+1)) )

        END DO
      END DO
    END DO

!
!----------------------------------------------------------------------
!
!   Determine the adjusted velocity, wa, store results in tem4:
!
!           2*rhobar*(Wa-Wm) + lamda1*z - dlamda2/dz = 0
!
!----------------------------------------------------------------------
!
    CALL avgsw(tem9,nx,ny,nz, 1,nx-1, 1,ny-1, tem6)

    WRITE(6,*) 'Maximum vertical velocity after adjustment'
    DO k = 2, nz-1
      maxw=0
      DO j = 2, ny-2
        DO i = 2, nx-2

!       zs       = z(k) - zshift
!
!        tem4(i,j,k)  =  w(i,j,k,tim) + (0.5/tem6(i,j,k))
!  :                     *(dzinv*(tem7(i,j,k)-tem7(i,j,k-1))
!  :                  -  0.25*zs*(tem8(i,j,k)+tem8(i+1,j,k)
!  :                  +  tem8(i,j+1,k)+tem8(i+1,j+1,k)))

          xs  = 0.5 * ( x(i) + x(i+1) ) - xmove
          ys  = 0.5 * ( y(j) + y(j+1) ) - ymove
          zs  =         z(k)            - zshift
          rrr = SQRT(xs*xs + ys*ys + zs*zs)

          tem4(i,j,k)  =  w(i,j,k,tim) + (0.5/tem6(i,j,k))              &
                       * ( dzinv*(tem7(i,j,k)-tem7(i,j,k-1))            &
                         - 0.25*(zs/rrr)*(tem8(i,j,k)+tem8(i+1,j,k)     &
                           +tem8(i,j+1,k)+tem8(i+1,j+1,k)) )

          IF (ABS(tem4(i,j,k)) > maxw) THEN
            maxw = ABS(tem4(i,j,k))
            imx=i
            jmx=j
          END IF
        END DO
      END DO

      WRITE(6,*) '  imx,jmx,k,maxw: ',imx,jmx,k,maxw

    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Exit adjustment/or direct insertion
!
!-----------------------------------------------------------------------
!

  270   CONTINUE

!
!-----------------------------------------------------------------------
!
!  Compare the model radial velocity and observed radial velocity
!
!-----------------------------------------------------------------------
!
  CALL avgx(tem2, 0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem6)
  CALL avgy(tem3, 0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem7)
  CALL avgz(tem4, 0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem8)

  WRITE(6,*) 'Check vr & print out if abs(vradj-vro) > 1.5 m/s'

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        xs    = 0.5*(x(i)+x(i+1)) - xmove
        ys    = 0.5*(y(j)+y(j+1)) - ymove
        zs    = 0.5*(z(k)+z(k+1)) - zshift

        rad   = SQRT(xs**2+ys**2+zs**2)
        vradj = (tem6(i,j,k)*xs + tem7(i,j,k)*ys                        &
                     + tem8(i,j,k)*zs)/rad
        vro   = tem1(i,j,k)/tem9(i,j,k)
        vrm   = tem5(i,j,k)/tem9(i,j,k)

!      IF(abs(vradj-vro).gt.1.5) THEN
!        write(6,*) '  i,j,k,vradj,vro,vrm: ',i,j,k,vradj,vro,vrm
!      ENDIF

      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  NOTE: This next section is designed to determine whether or not
!        mass conservation has been satisfied.  It is not an essential
!        part of the assimilation package.
!
!       Stagger rhostar(stored in tem9).
!
!-----------------------------------------------------------------------
!
  CALL avgsu(tem9,nx,ny,nz, 1,ny-1, 1,nz-1, tem6, tem10)
  CALL avgsv(tem9,nx,ny,nz, 1,nx-1, 1,nz-1, tem7, tem10)
  CALL avgsw(tem9,nx,ny,nz, 1,nx-1, 1,ny-1, tem8)
!
!-----------------------------------------------------------------------
!
!  Calculate del dot rho V using the adjusted velocities.
!
!-----------------------------------------------------------------------
!

  WRITE(6,*) 'background divergence'

  DO k=2,nz-2
    maxdiv = 0.0
    DO j=4,ny-4
      DO i=4,nx-4
        tem1(i,j,k)=                                                    &
            (u(i+1,j,k,tim)*tem6(i+1,j,k)-u(i,j,k,tim)*tem6(i,j,k))*dxinv &
            +(v(i,j+1,k,tim)*tem7(i,j+1,k)-v(i,j,k,tim)*tem7(i,j,k))*dyinv &
            +(w(i,j,k+1,tim)*tem8(i,j,k+1)-w(i,j,k,tim)*tem8(i,j,k))*dzinv
        IF (ABS(tem1(i,j,k)) > maxdiv) THEN
          maxdiv = ABS(tem1(i,j,k))
          imx = i
          jmx = j
        END IF
      END DO
    END DO
    maxdiv=tem1(imx,jmx,k)
    PRINT *,'  maxdiv,imx,jmx,k = ',maxdiv,imx,jmx,k
  END DO

  WRITE(6,*) 'adjusted divergence'

  DO k=2,nz-2
    maxdiv = 0.0
    DO j=4,ny-4
      DO i=4,nx-4
        tem1(i,j,k)=                                                    &
            (tem2(i+1,j,k)*tem6(i+1,j,k)-tem2(i,j,k)*tem6(i,j,k))*dxinv &
            +(tem3(i,j+1,k)*tem7(i,j+1,k)-tem3(i,j,k)*tem7(i,j,k))*dyinv &
            +(tem4(i,j,k+1)*tem8(i,j,k+1)-tem4(i,j,k)*tem8(i,j,k))*dzinv
        IF (ABS(tem1(i,j,k)) > maxdiv) THEN
          maxdiv = ABS(tem1(i,j,k))
          imx = i
          jmx = j
        END IF
      END DO
    END DO
    maxdiv=tem1(imx,jmx,k)
    PRINT *,'  maxdiv,imx,jmx,k = ',maxdiv,imx,jmx,k
  END DO

  RETURN
END SUBROUTINE assimvel
