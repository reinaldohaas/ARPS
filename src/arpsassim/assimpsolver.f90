!
!################################################################
!################################################################
!#####                                                      #####
!#####               SUBROUTINE RETRPOIS                    #####
!#####                                                      #####
!#####               Copyright (c) 1993                     #####
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!#####                                                      #####
!################################################################
!################################################################
!

SUBROUTINE retrpois(nx,ny,nz,k,                                         &
           biga,bigb,pc,                                                &
           tem4,pnew,d,alpha,beta,pprt)
!
!---------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine solves a 2-D Poisson equation:
!
!  d2p/dx2  +  d2p/dy2  =  dbiga/dx + dbigb/dy
!
!  subject to the Dirichlet conditions:
!
!  p = pbkg on the i = 2 and i = nx-1 physical boundaries
!
!  p = pbkg on the j = 2 and j = ny-1 physical boundaries
!
!
!  The solution is obtained by successively applying an alternating
!  direction implicit (ADI) scheme to solve tri-diagonal matrix
!  equations in the i and j directions.  This scheme is known as the
!  "sweeping method" or Thomas algorithm.  (see P. J. Roache, 1982,
!  "Computational Fluid Dynamics").
!
!---------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Alan Shapiro and Steve Lazarus
!  2/2/93.
!
!  MODIFICATION HISTORY:
!
!  6/18/96 (Limin Zhao and Alan Shapiro)
!  Added the option for Dirichlet boundary condition, which turns
!  out very good for properply retrieving pressure and potential
!  temperature outside of rainwater area.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    biga     Sum of terms in x-momentum equation (except for dp/dx)
!    bigb     Sum of terms in y-momentum equation (except for dp/dy)
!    pc       First guess for the solution of the Poisson equation.
!
!
!---------------------------------------------------------------------
!
!  OUTPUT:
!
!  pc        Solution of the Poisson equation.  pc differs from the
!            perturbation pressure by an arbitrary function of height.
!
!---------------------------------------------------------------------
!
!  WORK ARRAYS:
!
!    alpha    Coefficient in solution to a tri-diagonal matrix equation.
!    beta     Coefficient in solution to a tri-diagonal matrix equation.
!    tem4     Temporary work array.
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
  INCLUDE 'globcst.inc'     ! Model control constants
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz        ! Number of grid points in x, y and z directions

  REAL :: biga(nx,ny,nz)       ! Sum of terms in x-momentum equation
                               ! (except for dp/dx).

  REAL :: bigb(nx,ny,nz)       ! Sum of terms in y-momentum equation
                               ! (except for dp/dy)

  REAL :: pc(nx,ny,nz)         ! Solution of the Poisson equation.
                               ! pc differs from the perturbation pressure
                               ! by an arbitrary function of height.
  REAL :: tem4(nx,ny,nz)
  REAL :: pprt(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: iter      ! Iteration number
  INTEGER :: itmax     ! Maximum number of iterations
  INTEGER :: icon      ! Used in error tolerance check;
                       ! =0, iterates have converged, RETURN to RETRPTPR;
                       ! =1, iterates have not converged, continue iterating

  INTEGER :: inum      ! Number of grid points on a horizontal plane

  REAL :: rdenom       ! Reciprocal denominator, used to simplify a calculation
  REAL :: prms         ! Root mean square value of pc along an i or j line
  REAL :: tay          ! Over-relaxation coefficient
  REAL :: tol          ! An error tolerance parameter
  REAL :: resid        ! Discrepancy between updated  pc and old pc at a point
  REAL :: compar       ! Compare resid with the quantity compar (=tol*prms)
  REAL :: pave         ! Horizontal average of pc

  PARAMETER(itmax=500, tol=0.00001)

  REAL :: pnew(nx,ny)  ! Updated value of pc

  REAL :: d(nx,ny)     ! Inhomogeneous term in a tri-diagonal matrix equation

  REAL :: a, b, c      ! Coefficients in a tri-diagonal matrix equation

  REAL :: alpha(nx*ny) ! Coefficient in recursive solution to a tri-diagonal
                       ! matrix equation

  REAL :: beta(nx*ny)  ! Coefficient in recursive solution to a tri-diagonal
                       ! matrix equation

  REAL :: rdx2         ! Reciprocal of (dx)**2
  REAL :: rdy2         ! Reciprocal of (dy)**2

  REAL :: bc_opt
  INTEGER :: ips
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!----------------------------------------------------------------------
!
!  Estimate the optimum over-relaxation coefficient, tay, from
!  the formula (cf G.J. Haltiner and R.T. Williams, 1980: "Numerical
!  Prediction and Dynamic Meteorology, second edition":
!
!  tay = 2. - pi*sqrt(2/nx**2 + 2/ny**2)
!
!
!----------------------------------------------------------------------
!

  tay = 2. - 3.14159*SQRT((2./(nx*nx) + 2./(ny*ny)))

  PRINT *, 'estimated optimum over-relaxation coeff. =', tay

  rdx2 = dxinv*dxinv        ! reciprocal of (dx)**2
  rdy2 = dyinv*dyinv        ! reciprocal of (dy)**2

!
!-----------------------------------------------------------------------
!
!  Begin the big do-loop for iterating back and forth between the
!  solutions of the two tri-diagonal matrix equations.
!
!-----------------------------------------------------------------------
!
  DO iter = 1, itmax

    IF (iter == itmax) THEN
      PRINT 2000, k
      STOP 1
    END IF

    2000  FORMAT(/, 1X, 'Warning& ! Poisson solver did not converge.',  &
    &        /, 1X, 'klevel =', i5, /)

    IF (MOD(iter,20) == 0) PRINT 3000, k, iter

    3000  FORMAT(2X, 'klevel =', i5, 5X, 'iter = ', i5)
!
!-----------------------------------------------------------------------
!
!  First sweep in the i direction (put in an outer j do loop since we will
!  sweep in the i-direction for each value of j).  We are thus solving
!  a tri-diagonal matrix equation of the form:
!
!          a*pnew(i+1,j) + b*pnew(i,j) + c*pnew(i-1,j) = d(i,j)
!
!    where a, b, c and d are of the form:
!
!    a = c = 1/(dx*dx),    b = -2*(1/(dx*dx) + 1/(dy*dy)),
!
!    d = tem4 - (pc(i,j+1) + pc(i,j-1))/(dy*dy)
!
!
!-----------------------------------------------------------------------
!
    a = rdx2
    c = a
    b = -2.*(rdx2 + rdy2)

    DO j = 2, ny-2

      DO i = 2, nx-2
        d(i,j) = tem4(i,j,k) - rdy2*(pc(i,j+1,k) + pc(i,j-1,k))
      END DO
!
!-----------------------------------------------------------------------
!
!  The solution of the tri-diagonal matrix equation satisfies the
!  recursion relation:
!
!            pnew(i-1,j) = alpha(i)*pnew(i,j) + beta(i)
!
!  where alpha and beta satisfy the recursion formulae:
!
!            alpha(i) = -a(i)/(b + c*alpha(i-1));
!
!            beta(i) = (d(i-1,j) - c*beta(i-1))/(b + c*alpha(i-1))
!
!  Write the boundary values of the alpha and beta coefficients at
!  i=2 from the Neumann boundary condition on the west boundary.  Then
!  compute the alpha and beta coefficients along i-lines from the
!  recursion formulae.
!
!-----------------------------------------------------------------------
!
      bc_opt = 0

      IF (bc_opt == 1) THEN        !Neumann B.C.
        alpha(2) = 1.
        beta(2) = - biga(2,j,k)*dx
      ELSE
        alpha(2) = 0.0
        beta(2)  = pprt(1,j,k)
      END IF

      DO i = 3, nx-1
        rdenom = 1./(b + c*alpha(i-1))
        alpha(i) = - a*rdenom
        beta(i) = (d(i-1,j) - c*beta(i-1))*rdenom
      END DO

!
!-----------------------------------------------------------------------
!
!  Compute the updated pc (pnew), at nx-1 from the Neumann boundary
!  condition at the east boundary.  Sweep downwards in the decreasing
!  i-direction to get the updated pc (pnew) from the recursion formula.
!
!-----------------------------------------------------------------------
!
      IF (bc_opt == 1) THEN
        pnew(nx-1,j) = (biga(nx-1,j,k)*dx                               &
                     + beta(nx-1))/(1. - alpha(nx-1))
      ELSE
        pnew(nx-1,j) = pprt(nx-1,j,k)
      END IF

      DO i = nx-1, 2, -1
        pnew(i-1,j) = alpha(i)*pnew(i,j) + beta(i)
      END DO
!
!-----------------------------------------------------------------------
!
!  Compute the over-relaxed residuals and add them to pc.
!
!-----------------------------------------------------------------------
!

      DO i = 1, nx-1
        resid = tay*(pnew(i,j) - pc(i,j,k))
        pc(i,j,k) = pc(i,j,k) + resid
      END DO

    END DO

!
!-----------------------------------------------------------------------
!
!  Now sweep in the j direction (put in an outer i-do loop since we will
!  sweep in the j-direction for each value of i).  We are thus solving
!  a tri-diagonal matrix equation of the form:
!
!          a*pnew(i,j+1) + b*pnew(i,j) + c*pnew(i,j-1) = d(i,j)
!
!    where a, b, c and d are of the form:
!
!    a = c = 1/(dy*dy),    b = -2*(1/(dx*dx) + 1/(dy*dy)),
!
!    d = tem4 - (pc(i+1,j) + pc(i-1,j))/(dx*dx)
!
!-----------------------------------------------------------------------
!
    a = rdy2
    c = a
    b = -2.*(rdx2 + rdy2)

    DO i = 2, nx-2

      DO j = 2, ny-2
        d(i,j) = tem4(i,j,k) - rdx2*(pc(i+1,j,k) + pc(i-1,j,k))
      END DO
!
!-----------------------------------------------------------------------
!
!  The solution of the tri-diagonal matrix equation satisfies the
!  recursion relation:
!
!            pnew(i,j-1) = alpha(j)*pnew(i,j) + beta(j)
!
!  where alpha and beta satisfy the recursion formulae:
!
!            alpha(j) = -a(j)/(b + c*alpha(j-1));
!
!            beta(j) = (d(i,j-1) - c*beta(j-1))/(b + c*alpha(j-1))
!
!  Write the boundary values of the alpha and beta coefficients at
!  j=2 from the Neumann boundary condition on the south boundary.  Then
!  compute the alpha and beta coefficients along j-lines from the
!  recursion formulae.
!
!-----------------------------------------------------------------------
!
      IF (bc_opt == 1) THEN
        alpha(2) = 1.
        beta(2) = - bigb(i,2,k)*dy
      ELSE
        alpha(2) = 0.0
        beta(2)  = pprt(i,1,k)
      END IF

      DO j = 3, ny-1
        rdenom = 1./(b + c*alpha(j-1))
        alpha(j) = - a*rdenom
        beta(j) = (d(i,j-1) - c*beta(j-1))*rdenom
      END DO
!
!-----------------------------------------------------------------------
!
!  Compute the updated pc (pnew), at ny-1 from the Neumann boundary
!  condition at the north boundary.  Sweep downwards in the decreasing
!  j-direction to get the updated pc (pnew) from the recursion formula.
!
!-----------------------------------------------------------------------
!
      IF (bc_opt == 1) THEN
        pnew(i,ny-1) = (bigb(i,ny-1,k)*dy                               &
                     + beta(ny-1))/(1. - alpha(ny-1))
      ELSE
        pnew(i,ny-1) = pprt(i,ny-1,k)
      END IF

      DO j = ny-1, 2, -1
        pnew(i,j-1) = alpha(j)*pnew(i,j) + beta(j)
      END DO
!
!-----------------------------------------------------------------------
!
!  Compute the rms value of the new pc along a j-line for the error
!  tolerance comparison.  Also compute the over-relaxed residuals and
!  add them to pc.
!
!-----------------------------------------------------------------------
!
      prms = 0.
      DO j = 1, ny-1
        prms = prms + pnew(i,j)*pnew(i,j)
      END DO
      prms = SQRT(prms/(ny-1))
      compar = tol*prms

      icon = 0
      DO j = 1, ny-1
        resid = tay*(pnew(i,j) - pc(i,j,k))
        icon = icon + IFIX(ABS(resid/compar))
        pc(i,j,k) = pc(i,j,k) + resid
      END DO

    END DO

    IF (icon == 0) RETURN
!
!-----------------------------------------------------------------------
!
!  Since the solution of Poisson's equation with Neumann boundary
!  conditions isn't unique (there's an arbitrary constant in the
!  solution), subtract off the average value of pc to avoid a potential
!  drift that could affect the convergence.
!
!-----------------------------------------------------------------------
!
    IF (bc_opt /= 1) CYCLE
    WRITE(6,*)'subtract mean pressue'
    pave = 0.
    inum = 0
    DO i = 1, nx-1
      DO j = 1, ny-1
        pave = pave + pc(i,j,k)
        inum = inum + 1
      END DO
    END DO
    pave = pave/FLOAT(inum)

    DO i = 1, nx-1
      DO j = 1, ny-1
        pc(i,j,k) = pc(i,j,k) - pave
      END DO
    END DO

  END DO

  RETURN
END SUBROUTINE retrpois
