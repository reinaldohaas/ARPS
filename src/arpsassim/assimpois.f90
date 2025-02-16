!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE POIS3D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pois3d(nx,ny,nz,dx,dy,dz,eps,hole_id,                        &
           tem3,tem4,tem5,tem6,tem7)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:This program is designed to fill-in the holes for missing
!          single-Doppler data on a Cartesian grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Code originally authored by Tzvi Gal-Chen with mods by
!          Mei Xu 12-93, and Steven Lazarus 6-94.
!
!
!  MODIFICATION HISTORY:
!
!  06/94 (Steven Lazarus)
!  All logical statements removed.
!
!  03/10/96 (Limin Zhao)
!  Added an option to use a 3-D/2-D hole-filler
!
!-----------------------------------------------------------------------
!
!  INPUT ARRAYS:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    dx       Grid spacing in the x-direction (east/west)
!    dy       Grid spacing in the y-direction (north/south)
!    dz       Grid spacing in the z-direction (vertical)
!
!    tem3     Right hand side of the Poisson Eqn.
!    tem4     Radial velocity data
!
!  OUTPUT ARRAYS:
!
!    tem5     Data grid to be filled
!
!  WORK ARRAYS:
!
!    tem6     Check for convergence, where
!    tem7     (n,1) = Difference in tem1 between two successive iter.
!             (n,2) = tem1 at previous iteration
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  INTEGER :: iprnt,ier

  REAL :: tem3  (nx,ny,nz)     ! RHS of Poisson Equation
  REAL :: tem4  (nx,ny,nz)     ! Radial velocity data
  REAL :: tem5  (nx,ny,nz)     ! Data grid to be filled
  REAL :: tem6  (nx,ny,nz)     ! Work array
  REAL :: tem7  (nx,ny,nz)     ! Work array

  PARAMETER(iprnt=10)       ! Controls frequency of printed output

  REAL :: dx,dy,dz             ! Grid spacing
  REAL :: tetx,tety            ! Coefficients of the 3-D Poisson
  REAL :: tetz,delt2           ! equation

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k             ! Do-loop indices
  INTEGER :: iter,itmx,ki      ! Iteration control parameters
  INTEGER :: istor,jstor,kstor ! Grid location of max error

  REAL :: omeg,eps             ! Relaxation coef.,error tolerance
  REAL :: tol                  ! User-specified tolerance level
  REAL :: denom                ! Coef. ued in 3-D solver
  REAL :: dx2,dy2,dz2          ! Square of the grid spacing

  REAL :: hole_id
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'assim.inc'       ! Assim/Retr control parameters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  PRINT *,'Entering POIS3D '
  PRINT *,'The tolerance level was set at ',eps
  PRINT *,'The option for hole-filler is ',hole_id
!
!-----------------------------------------------------------------------
!
!   Initialize the variables
!
!-----------------------------------------------------------------------
!
  DO i=1,nx
    DO j=1,ny
      DO k=1,nz
        tem6(i,j,k) =0.0
        tem7(i,j,k) =0.0
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!   Define iteration constants:
!
!  itmx   maximum number of iterations
!  omeg   relaxation factor
!  eps    solver convergence criteria
!  denom  2.0*[(dy**2)*(dz**2)+(dx**2)*(dz**2)+(dx**2)*(dy**2)]
!  tetx   (dy**2)*(dz**2)/denom
!  tety   (dx**2)*(dz**2)/denom
!  tetz   (dx**2)*(dy**2)/denom
!  delt2  (dx**2)*(dy**2)*(dz**2)/denom
!
!-----------------------------------------------------------------------
!
  itmx=2000
  omeg=1.90

  dx2=dx*dx
  dy2=dy*dy
  dz2=dz*dz

  denom=2.0*(dy2*dz2+dx2*dz2+dx2*dy2)
  tetx=dy2*dz2/denom
  tety=dx2*dz2/denom
  tetz=dx2*dy2/denom
  delt2=dx2*dy2*dz2/denom

  IF (hole_id == 2) THEN
    denom=2.0*(dy2+dx2)
    tetx=dy2/denom
    tety=dx2/denom
    tetz=0.0
    delt2=dx2*dy2/denom
    WRITE(6,*)'I am doing 2-D hole-filling'
  END IF

  WRITE(6,*) 'dx2,dy2, dz2: ',dx2,dy2, dz2
  WRITE(6,*) 'denom,tetx,tety,tetz,delt2: ',                            &
              denom,tetx,tety,tetz,delt2
!
!-----------------------------------------------------------------------
!
!   Begin iterations
!
!-----------------------------------------------------------------------
!
  iter =0

  50    iter =iter + 1

!
!-----------------------------------------------------------------------
!
!   Store the current field every iteration
!
!-----------------------------------------------------------------------
!
  ki = MOD(iter,1)

  IF(ki == 0) THEN

    DO i=1,nx-1
      DO j=1,ny-1
        DO k=1,nz-1
          tem7(i,j,k) = tem5(i,j,k)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!   Terminate the iterations after itmx times
!
!-----------------------------------------------------------------------
!
  IF(iter > itmx) THEN
    WRITE(6,60) itmx
    60      FORMAT(' after',i4,' iterations tem5 did not converge')
    ier =1
    WRITE(6,*)'Please increase iterations or check the code'
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!   Begin iterations.  Check the template (tem4) to ensure that the
!   poisson solver is applied to "non-data" points only. To preserve
!   the original data, tem5 is filled in the tem4 non-data regions and
!   is equal to tem4 in the data regions.
!
!   Fill the non-data points in the interior area.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        IF(tem4(i,j,k) == spval) THEN
          tem5(i,j,k)=(tetx*(tem5(i+1,j,k)+tem5(i-1,j,k))               &
                     + tety*(tem5(i,j+1,k)+tem5(i,j-1,k))               &
                     + tetz*(tem5(i,j,k+1)+tem5(i,j,k-1))               &
                     - delt2*tem3(i,j,k))*omeg+(1.-omeg)*tem5(i,j,k)
        END IF

      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Handle the perimeter points separately. Here we apply mixed
!  boundary conditions, with zero normal gradient at non-data
!  points (Neumann) and Dirichlet at data points.
!
!  NOTE:  The choice of which boundary condition is determined
!         in the driver ASSIMDRIV.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-2
    DO i=2,nx-2
      IF(tem4(i,1,k) == spval) THEN
        tem5(i,1,k) = tem5(i,2,k)
      END IF
      IF(tem4(i,ny-1,k) == spval) THEN
        tem5(i,ny-1,k)= tem5(i,ny-2,k)
      END IF
    END DO
  END DO

  DO k=2,nz-2
    DO j=2,ny-2
      IF(tem4(1,j,k) == spval) THEN
        tem5(1,j,k) = tem5(2,j,k)
      END IF
      IF(tem4(nx-1,j,k) == spval) THEN
        tem5(nx-1,j,k)= tem5(nx-2,j,k)
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!   Apply the top and bottom boundary conditions.
!
!-----------------------------------------------------------------------
!
  DO i= 2,nx-2
    DO j= 2,ny-2
      IF(tem4(i,j,1) == spval) THEN
        tem5(i,j,1) = tem5(i,j,2)
      END IF
      IF(tem4(i,j,nz-1) == spval) THEN
        tem5(i,j,nz-1) = tem5(i,j,nz-2)
      END IF
    END DO
  END DO
  DO j = 2, ny-2
    IF(tem4(1,j,1) == spval) THEN
      tem5(1,j,1) = tem5(2,j,1)
    END IF
    IF(tem4(1,j,nz-1) == spval) THEN
      tem5(1,j,nz-1) = tem5(2,j,nz-1)
    END IF
    IF(tem4(nx-1,j,nz-1) == spval) THEN
      tem5(nx-1,j,nz-1) = tem5(nx-2,j,nz-1)
    END IF
    IF(tem4(nx-1,j,1) == spval) THEN
      tem5(nx-1,j,1) = tem5(nx-2,j,1)
    END IF
  END DO

  DO i = 2, nx-2
    IF(tem4(i,1,1) == spval) THEN
      tem5(i,1,1) = tem5(i,2,1)
    END IF
    IF(tem4(i,1,nz-1) == spval) THEN
      tem5(i,1,nz-1) = tem5(i,2,nz-1)
    END IF
    IF(tem4(i,ny-1,1) == spval) THEN
      tem5(i,ny-1,1)    = tem5(i,ny-2,1)
    END IF
    IF(tem4(i,ny-1,nz-1) == spval) THEN
      tem5(i,ny-1,nz-1) = tem5(i,ny-2,nz-1)
    END IF
  END DO

!
!-----------------------------------------------------------------------
!
!  Handle the corner points separately
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    IF(tem4(1,1,k) == spval) THEN
      tem5(1,1,k) = tem5(2,2,k)
    END IF
    IF(tem4(1,ny-1,k) == spval) THEN
      tem5(1,ny-1,k) = tem5(2,ny-2,k)
    END IF
    IF(tem4(nx-1,1,k) == spval) THEN
      tem5(nx-1,1,k) = tem5(nx-2,2,k)
    END IF
    IF(tem4(nx-1,ny-1,k) == spval) THEN
      tem5(nx-1,ny-1,k) = tem5(nx-2,ny-2,k)
    END IF
  END DO

!-----------------------------------------------------------------------
!
!    Check convergence every iteration
!
!-----------------------------------------------------------------------
!
  ki=MOD(iter,1)

  IF(ki /= 0) GO TO 50

  DO i=2,nx-2
    DO j=2,ny-2
      DO k=2,nz-2
        IF(tem4(i,j,k) == spval) THEN
          tem6(i,j,k)= ABS(tem5(i,j,k)-tem7(i,j,k))
        END IF
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!   Determine the largest value of tem6(i,j,k) and store its location
!
!-----------------------------------------------------------------------
!
  tol =0.0

  DO i=2,nx-2
    DO j=2,ny-2
      DO k=2,nz-2
        IF(tem4(i,j,k) == spval) THEN
          tol=MAX(tol,tem6(i,j,k))
          IF(tol == tem6(i,j,k)) THEN
            istor = i
            jstor = j
            kstor = k
          END IF
        END IF
      END DO
    END DO
  END DO

  IF((MOD(iter,iprnt)) == 0) THEN
    PRINT 110, iter, tol
    110     FORMAT(3X,'At iteration',i4,' tol is',g16.8)
    PRINT 115, istor,jstor,kstor
    115     FORMAT(3X,'The max tol appears at (i,j,k)=',3I4)
  END IF

!
!-----------------------------------------------------------------------
!
!   If the tolerance is greater than the acceptable error (eps)
!   continue iterations.
!
!-----------------------------------------------------------------------
!
  IF((tol >= eps).AND.(iter < itmx)) GO TO 50
!
!-----------------------------------------------------------------------
!
!    Cease iterations when the tolerance is less than the acceptable
!    error (eps). Solution has converged.
!
!-----------------------------------------------------------------------
!
  IF(tol < eps) THEN
    WRITE(6,120) iter
    120     FORMAT('Solution converged, # of iterations in pois3d was ',i4)
    ier=0
    RETURN
  END IF

!
!-----------------------------------------------------------------------
!
!    Cease iterations if the tolerance is greater than the acceptable
!    error and the maximum number of iterations (itmx) has been exceeded.
!
!-----------------------------------------------------------------------
!
  WRITE(6,130) iter
  130   FORMAT(' After',i4,' iterations tem5(i,j,k) did not converge')
  ier=1

  RETURN
END SUBROUTINE pois3d
