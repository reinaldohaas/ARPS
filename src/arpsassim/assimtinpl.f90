!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RETRINT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE retrint(nx,ny,nz,                                            &
           u,v,w,ptprt,pprt,qv,qc,qr,                                   &
           t1,t2,t3,                                                    &
           tem1,tem2,tem3,tem4,tem5,                                    &
           tem6,tem7,tem8,tem9)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine orchestrates the quadratic interpolation of input data
!  at three time levels, t1, t2 and t3 (where t1 < t2 < t3) to three
!  intermediate time levels.  The time intervals between data availability,
!  t2-t1 and t3-t2 are not necessarily equal to each other.  The Lagrange
!  three point formula for unequal abscissas is used for the interpolation.
!  (cf page 27 of Carnahan, B., Luther, H. A. and Wilkes, J. O., 1969:
!  Applied Numerical Methods).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Alan Shapiro and Steve Lazarus
!  1/21/93.
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction
!
!    u        x component of velocity (m/s) at times t1, t2 and t3
!    v        y component of velocity (m/s) at times t1, t2 and t3
!    w        z component of velocity (m/s) at times t1, t2 and t3
!
!    ptprt    Perturbation potential temperature (K) at times t1, t2 and t3
!    pprt     Perturbation pressure (Pascals) at times t1, t2 and t3
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rain water mixing ratio (kg/kg)
!
!    t1       Time of first data file
!    t2       Time of second data file
!    t3       Time of third data file
!
!  OUTPUT:
!
!    u        Time-interpolated x component of velocity (m/s)
!    v        Time-interpolated y component of velocity (m/s)
!    w        Time-interpolated z component of velocity (m/s)
!
!    ptprt    Time-interpolated perturbation potential temperature (K)
!    pprt     Time-interpolated perturbation pressure (Pascals)
!    qv       Time-interpolated water vapor specific humidity (kg/kg)
!    qc       Time-interpolated cloud water mixing ratio (kg/kg)
!    qr       Time-interpolated rain water mixing ratio (kg/kg)
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
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Forces explicit declarations

  INTEGER :: nt                ! Number of time levels of time-dependent
                               ! arrays.
  PARAMETER (nt=3)

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)

  REAL :: t1                   ! Time of first data file
  REAL :: t2                   ! Time of second data file
  REAL :: t3                   ! Time of third data file

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascals)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  REAL :: t            ! Model time level to which the variables are
                       ! interpolated.  t takes on three values in succession.
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Model control parameters.
  INCLUDE 'assim.inc'       ! Velocity insertion/thermodynamic recovery
                            ! control parameters.
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL tinterp
!
!-----------------------------------------------------------------------
!
!  INTERP performs the Lagrange three-point interpolation for
!  unequal abscissas. The 3-point formula is of the form:
!
!       f(t) = coeff1*f(t1) + coeff2*f(t2) + coeff3*f(t3),
!
!  where f(t1), f(t2) and f(t3) are the known values of f(t) at t=t1,
!  t=t2 and t=t3, respectively, and
!
!       coeff1 = (t3-t)*(t2-t)/((t3-t1)*(t2-t1))
!
!       coeff2 = (t-t1)*(t3-t)/((t2-t1)*(t3-t2))
!
!       coeff3 = (t-t2)*(t-t1)/((t3-t2)*(t3-t1))
!
!
!  The input variable to INTERP is 4-D (i.e., time-dependent) whereas
!  the output variable is 3-D.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!

  PRINT *, ' Inside RETRINT. curtim, nstep =', curtim, nstep
  PRINT *, ' t1 =', t1, 't2 =', t2, 't3 =', t3

!
!-----------------------------------------------------------------------
!
!  Print a warning message if curtim + dtbig is greater than t3 or
!  curtim - dtbig is less than t1.  If either of these situations
!  occurs the algorithm extrapolates rather than interpolates.
!
!-----------------------------------------------------------------------
!

  IF (curtim + dtbig > t3) THEN

    PRINT *, 'Warning from RETRINT: curtim + dtbig > t3'

  ELSE IF (curtim - dtbig < t1) THEN

    PRINT *, 'Warning from RETRINT: curtim - dtbig < t1'

  END IF

  IF (t1 > t2 .OR. t1 > t3) THEN

    PRINT *, 'Warning from RETRINT: t1 > t2 or t1 > t3'
    STOP

  ELSE IF (t2 > t3) THEN

    PRINT *, 'Warning from RETRINT: t2 > t3'
    STOP

  END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate the three velocity components to the time curtim - dtbig.
!  Store results in temporary arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim - dtbig

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,u,tem1)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,v,tem2)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,w,tem3)
!
!-------------------------------------------------------------------------
!
!  Interpolate the three velocity components to the time curtim.
!  Store results in temporary arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,u,tem4)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,v,tem5)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,w,tem6)

!
!-------------------------------------------------------------------------
!
!  Interpolate the three velocity components to the time curtim + dtbig.
!  Store results in temporary arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim + dtbig

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,u,tem7)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,v,tem8)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,w,tem9)

!
!-------------------------------------------------------------------------
!
!  Put interpolated velocity fields (stored in temporary arrays) into
!  the velocity arrays.  This overwrites the old velocity arrays.
!
!-------------------------------------------------------------------------
!
  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx
        u(i,j,k,1) = tem1(i,j,k)
        v(i,j,k,1) = tem2(i,j,k)
        w(i,j,k,1) = tem3(i,j,k)

        u(i,j,k,2) = tem4(i,j,k)
        v(i,j,k,2) = tem5(i,j,k)
        w(i,j,k,2) = tem6(i,j,k)

        u(i,j,k,3) = tem7(i,j,k)
        v(i,j,k,3) = tem8(i,j,k)
        w(i,j,k,3) = tem9(i,j,k)
      END DO
    END DO
  END DO
!
!-------------------------------------------------------------------------
!
!  Interpolate the perturbation pressure and perturbation potential
!  temperature to the time curtim - dtbig.  Store results in temporary
!  arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim - dtbig

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,ptprt,tem1)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,pprt,tem2)

!
!-------------------------------------------------------------------------
!
!  Interpolate the perturbation pressure and perturbation potential
!  temperature to the time curtim.  Store results in temporary arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,ptprt,tem3)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,pprt,tem4)

!
!-------------------------------------------------------------------------
!
!  Interpolate the perturbation pressure and perturbation potential
!  temperature to the time curtim + dtbig.  Store results in temporary
!  arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim + dtbig

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,ptprt,tem5)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,pprt,tem6)

!
!-------------------------------------------------------------------------
!
!  Put interpolated perturbation pressure and interpolated perturbation
!  potential temperature (stored in temporary arrays) into the pprt and
!  ptprt arrays.  This overwrites the old pprt and ptprt arrays.
!
!-------------------------------------------------------------------------
!
  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx
        ptprt(i,j,k,1) = tem1(i,j,k)
        pprt(i,j,k,1) = tem2(i,j,k)

        ptprt(i,j,k,2) = tem3(i,j,k)
        pprt(i,j,k,2) = tem4(i,j,k)

        ptprt(i,j,k,3) = tem5(i,j,k)
        pprt(i,j,k,3) = tem6(i,j,k)
      END DO
    END DO
  END DO

!
!-------------------------------------------------------------------------
!
!  Interpolate the water vapor, cloud and rain water mixing ratios, qv,
!  qc, and qr, to the time curtim - dtbig.  Store results in temporary
!  temporary arrays.
!
!-------------------------------------------------------------------------
!
  t = curtim - dtbig

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qv,tem1)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qc,tem2)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qr,tem3)

!
!-------------------------------------------------------------------------
!
!  Interpolate the water vapor, cloud and rain water mixing ratios, qv,
!  qc and qr, to the time curtim.  Store results in temporary arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qv,tem4)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qc,tem5)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qr,tem6)

!
!-------------------------------------------------------------------------
!
!  Interpolate the water vapor, cloud and rain water mixing ratios, qv,
!  qc and qr, to the time curtim + dtbig.  Store results in temporary
!  arrays.
!
!-------------------------------------------------------------------------
!

  t = curtim + dtbig

  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qv,tem7)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qc,tem8)
  CALL tinterp(nx,ny,nz,t,t1,t2,t3,qr,tem9)

!
!-------------------------------------------------------------------------
!
!  Put interpolated water vapor, cloud and rain water mixing ratios
!  (stored in temporary arrays) into the qv, qc and qr arrays.  This
!  overwrites the old qv, qc and qr arrays.
!
!-------------------------------------------------------------------------
!
  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx

        qv(i,j,k,1) = tem1(i,j,k)
        qc(i,j,k,1) = tem2(i,j,k)
        qr(i,j,k,1) = tem3(i,j,k)

        qv(i,j,k,2) = tem4(i,j,k)
        qc(i,j,k,2) = tem5(i,j,k)
        qr(i,j,k,2) = tem6(i,j,k)

        qv(i,j,k,3) = tem7(i,j,k)
        qc(i,j,k,3) = tem8(i,j,k)
        qr(i,j,k,3) = tem9(i,j,k)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE retrint
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTERP                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tinterp(nx,ny,nz,t,t1,t2,t3,var,varint)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine quadratically interpolates the values of input data
!  at three time levels, t1, t2 and t3 (where t1 < t2 < t3) to the time t.
!  intermediate time levels.  The time intervals between data availability,
!  t2-t1 and t3-t2 are not necessarily equal to each other.  The Lagrange
!  three point formula for unequal abscissas is used for the interpolation.
!  (cf page 27 of Carnahan, B., Luther, H. A. and Wilkes, J. O., 1969:
!  Applied Numerical Methods).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Alan Shapiro and Steve Lazarus
!  1/21/93.
!
!  MODIFICATION HISTORY:
!
!  4/19/93 (Alan Shapiro and Steve Lazarus)
!
!    Included special cases where t1 = t2 or t2 = t3.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction
!
!
!    t        Time level for data interpolation
!    t1       Time of first data file
!    t2       Time of second data file
!    t3       Time of third data file
!
!    var      A 4-D (time-dependent )array to be interpolated to time t
!
!  OUTPUT:
!
!    varint   var after being interpolated to time t, a 3-D array
!
!
!  WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Forces explicit declarations

  INTEGER :: nt                ! Number of time levels of time-dependent
                               ! arrays.
  PARAMETER (nt=3)

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z directions

  REAL :: t                    ! Model time level to which the variables are
                               ! interpolated.
  REAL :: t1                   ! Time of first data file
  REAL :: t2                   ! Time of second data file
  REAL :: t3                   ! Time of third data file

  REAL :: var(nx,ny,nz,nt)     ! A variable to be interpolated to time t

  REAL :: varint(nx,ny,nz)     ! var after being interpolated to time t
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k


  REAL :: coeff1       ! Coefficient of function of time at t = t1 appearing
                       ! in the 3-point Lagrange interpolation formula

  REAL :: coeff2       ! Coefficient of function of time at t = t2 appearing
                       ! in the 3-point Lagrange interpolation formula

  REAL :: coeff3       ! Coefficient of function of time at t = t3 appearing
                       ! in the 3-point Lagrange interpolation formula
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
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
!  The Lagrange three-point interpolation formula for unequal abscissas
!  is used for the time interpolation.  This 3-point formula is of the
!  form:
!
!       f(t) = coeff1*f(t1) + coeff2*f(t2) + coeff3*f(t3),
!
!  where f(t1), f(t2) and f(t3) are the known values of f(t) at t=t1,
!  t=t2 and t=t3, respectively, and
!
!       coeff1 = (t3-t)*(t2-t)/((t3-t1)*(t2-t1))
!
!       coeff2 = (t-t1)*(t3-t)/((t2-t1)*(t3-t2))
!
!       coeff3 = (t-t2)*(t-t1)/((t3-t2)*(t3-t1))
!
!
!-----------------------------------------------------------------------
!
  IF (t1 == t2 .AND. t2 == t3) THEN        ! Only 1 time level of
                                           ! information available.
    coeff1 = 1.
    coeff2 = 0.
    coeff3 = 0.

  ELSE IF ((t1 == t2 .AND. t2 /= t3) .OR.   & ! Only 2 time levels of information
        (t1 /= t2 .AND. t2 == t3)) THEN  ! available.  Linearly interpolate
                                         ! between t1 and t3.
    coeff1 = (t3 - t)/(t3 - t1)
    coeff2 = 0.
    coeff3 = (t - t1)/(t3 - t1)

  ELSE IF (t1 /= t2 .AND. t2 /= t3) THEN    ! 3 time levels of information
                                            ! available.  Use 3-point formula.

    coeff1 = (t3 - t)*(t2 - t)/((t3 - t1)*(t2 - t1))
    coeff2 = (t - t1)*(t3 - t)/((t2 - t1)*(t3 - t2))
    coeff3 = (t - t2)*(t - t1)/((t3 - t2)*(t3 - t1))

  END IF

  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx

        varint(i,j,k) = coeff1*var(i,j,k,1) + coeff2*var(i,j,k,2) +     &
                        coeff3*var(i,j,k,3)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE tinterp
