!
!################################################################
!################################################################
!#####                                                      #####
!#####                 SUBROUTINE RETRFZ                    #####
!#####                                                      #####
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!#####                                                      #####
!################################################################
!################################################################
!

SUBROUTINE retrfz(nz,nn,                                                &
           fzero,a,b,f,term1,term2)
!
!---------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine evaluates a function f(z) that is related to the
!  vertical structure of the perturbation pressure.  f(z) is determined
!  as the analytic solution of the first order linear ordinary differential
!  equation (o.d.e.):
!
!  df/dz + a(z) f + b(z) = 0.
!
!  The general solution is (cf Braun, 1975: Dif. Eqns and Their
!  Applications):
!
!  f(z) = exp(-integral(a(z))) * [f(0) - integral(b*exp(+integral(a)))]
!
!  The lower limit of the integration is the first scalar grid point
!  above the lower physical boundary; that is, at z = dz/2 (k=2).
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Alan Shapiro and Steve Lazarus
!  2/2/93.
!
!
!  MODIFICATION HISTORY:
!
!
!---------------------------------------------------------------------
!
!  INPUT:
!
!    nz       Number of grid points in the z-direction (vertical)
!    nn       Dimension of a, b and f
!
!    a        Coefficient a(z) in the o.d.e. for f(z):
!             df/dz + a(z)f + b(z) = 0.
!             a is defined on the scalar points.
!
!    b        Coefficient b(z) in the o.d.e. for f(z):
!             df/dz + a(z)f + b(z) = 0.
!             b is defined on the w points.
!
!    fzero    Constant of integration: Value of f(z) at the
!             first scalar grid point above the lower physical
!             boundary; that is, fzero = f(dz/2).
!
!---------------------------------------------------------------------
!
!  OUTPUT:
!
!    f(z)          The solution of the o.d.e..  f(z) is valid on scalar
!                  points from the first scalar point above the lower
!                  physical boundary; that is, at z = dz/2 (k=2), to the
!                  first scalar point beneath the upper physical boundary;
!                  that is, at k=nz-2.
!
!
!---------------------------------------------------------------------
!
!  WORK ARRAYS:
!
!   term1(nn)
!   term2(nn)
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
!
!---------------------------------------------------------------------
!
!  Variable Declarations
!
!---------------------------------------------------------------------
!
  IMPLICIT NONE     ! Force explicit declarations

  INTEGER :: nz        ! Number of grid points in the z-direction (vertical)
  INTEGER :: nn        ! Dimension of a, b and f

  REAL :: f(nn)        ! Function related to the vertical structure of the
                       ! perturbation pressure.  Determined as the solution
                       ! of the o.d.e.:   df/dz + a(z)f + b(z) = 0

  REAL :: a(nn)        ! Coefficient multiplying f in the o.d.e. for f

  REAL :: b(nn)        ! Inhomogeneous term in the o.d.e. for f

  REAL :: fzero        ! Constant of integration. fzero = f(dz/2)
!
!---------------------------------------------------------------------
!
!  Misc. local variables:
!
!---------------------------------------------------------------------
!
  INTEGER :: k

  REAL :: gral         ! Integral of a(z)
  REAL :: term1(nz-2)  ! exp(gral)
  REAL :: term2(nz-2)  ! Integral of b(z)*term1
!
!----------------------------------------------------------------------
!
!  Include files:
!
!----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
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
!  Evaluate the analytic solution for f(z).  Use trapezoidal rule
!  to evaluate the integrals. f(z) and a(z) are valid on scalar
!  points. b(z) is valid on w points.
!
!  Note:  the integrals extend from the first scalar grid point
!  above the lower physical boundary (k=2); to the first scalar grid point
!  beneath the upper physical boundary (k=nz-2).
!
!----------------------------------------------------------------------
!
  gral = 0.
  term1(2) = 1.
  DO k = 3, nz-2
    gral = gral + a(k)*dz
    term1(k) = EXP(gral)
  END DO

  term2(2) = 0.
  DO k = 3, nz-2
    term2(k) = term2(k-1) + b(k)*term1(k)*dz
  END DO

  DO k = 2, nz-2
    f(k) = (fzero - term2(k))/term1(k)
  END DO

  RETURN
END SUBROUTINE retrfz
