!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RETRCHK                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE retrchk(nx,ny,nz,biga,bigb,pc)
!
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine checks how well the retrieved pressure gradients
!  fit the individual momentum equations. For a description of this
!  technique see:
!
!       Chapter 13 of Instruments and Techniques for Thunderstorm
!       Observation and Analysis;  Vol. 3 of Thunderstorms:  A Social,
!       Scientific and Technological Documentary, 2nd Edition. 1988.
!       University of Oklahoma Press.
!
!       Gal-Chen, T., and R.A. Kropfli, 1984.  Buoyancy and Pressure
!       Perturbations Derived from Dual-Doppler Radar Observations
!       of the Planetary Boundary Layer: Applications for Matching
!       Models with Observations. JAS Vol. 41, pp. 3007-20.
!
!---------------------------------------------------------------------------
!
!  AUTHORS: Steve Lazarus and Alan Shapiro
!           4/12/93
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
!    nz       Number of grid points in the z-direction (vertical)
!
!    pc       Computed perturbation pressure (up to an arbitrary
!             function of height).
!
!    biga     Estimate of dp/dx
!    bigb     Estimate of dp/dy
!
!  OUTPUT:
!
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
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y, z directions

  REAL :: pc(nx,ny,nz)         ! Computed perturbation pressure (up to
                               ! an arbitrary function of height).

  REAL :: biga  (nx,ny,nz)     ! Estimate of dp/dx
  REAL :: bigb  (nx,ny,nz)     ! Estimate of dp/dy
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: conchk, numer, denom, stuff1, stuff2
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'grid.inc'
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
!  Compute a nondimensional parameter conchk (CONsistency CHeck)
!  quantifying how well the retrieved pressure gradients fit the
!  equations of motion.  conchk is defined by the formula:
!
!  conchk = numer/denom
!
!  where
!
!  numer = volume integral of [(dp/dx - biga)**2 + (dp/dy - bigb)**2]
!  denom = volume integral of (biga**2 + bigb**2)
!
!
!-----------------------------------------------------------------------
!
  numer = 0.
  denom = 0.

  DO k = 2, nz-2
    DO j = 2, ny-2
      DO i = 2, nx-2
        stuff1 = dxinv*(pc(i,j,k) - pc(i-1,j,k)) - biga(i,j,k)
        stuff2 = dyinv*(pc(i,j,k) - pc(i,j-1,k)) - bigb(i,j,k)

        numer = numer + stuff1*stuff1 + stuff2*stuff2
        denom = denom + biga(i,j,k)*biga(i,j,k) +                       &
                        bigb(i,j,k)*bigb(i,j,k)

      END DO
    END DO
  END DO

  IF (denom /= 0.) conchk = numer/denom
  IF (denom == 0.) PRINT *, 'warning! denom = 0, conchk blows up'

  PRINT *, '          '
  PRINT *, ' Conchk = ', conchk

  RETURN
END SUBROUTINE retrchk
