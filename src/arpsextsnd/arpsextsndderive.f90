!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TEMPER                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE temper ( nx,ny,nz, ptbar, ptprt, ppert, pbar, t )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Using a version of Poisson's formula, calculate temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Joe Bradley
!    12/05/91
!
!  MODIFICATIONS:
!    Modified by Ming Xue so that arrays are only defined at
!             one time level.
!    6/09/92  Added full documentation and phycst include file for
!             rddcp=Rd/Cp  (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    theta    Potential temperature (degrees Kelvin)
!    ppert    Perturbation pressure (Pascals)
!    pbar     Base state pressure (Pascals)
!
!  OUTPUT:
!
!    t        Temperature (degrees Kelvin)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  REAL :: ptbar(nx,ny,nz)      ! potential temperature (degrees Kelvin)
  REAL :: ptprt(nx,ny,nz)      ! potential temperature (degrees Kelvin)
  REAL :: ppert(nx,ny,nz)      ! perturbation pressure (Pascals)
  REAL :: pbar (nx,ny,nz)      ! base state pressure (Pascals)

  REAL :: t    (nx,ny,nz)      ! temperature (degrees Kelvin)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
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
!  Calculate the temperature using Poisson's formula.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        t(i,j,k) = ( ptbar(i,j,k) + ptprt(i,j,k) ) *                    &
                   (((ppert(i,j,k) + pbar(i,j,k)) / p0) ** rddcp)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE temper

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_rfc                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate rfc value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!  Ming Xue (10/16/2001)
!  Now passing in precalculated reflectivity field instead of calculating
!  it inside.
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_rfc(nx, ny, nz, ref, refc)

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN ) :: ref (nx,ny,nz) ! Reflectivity
  REAL,    INTENT(OUT) :: refc(nx,ny)    ! Composite reflectivity

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j=1,ny
    DO i=1,nx
      refc(i,j)= ref(i,j,1)
      DO k=2,nz-1
        refc(i,j) = MAX(refc(i,j),ref(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_rfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vic                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate vic
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vic(tem7,qscalar,rhobar,zp,nx,ny,nz,nscalar,nscalarq,tem6)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nscalar,nscalarq
  REAL    :: tem7(nx,ny)
  REAL    :: qscalar(nx,ny,nz,nscalar)
  REAL    :: rhobar(nx,ny,nz), zp(nx,ny,nz)
  REAL    :: tem6(nx,ny,nz)

  INTEGER :: i,j,k,nq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  tem6 = 0.0

  DO j=1,ny
    DO i=1,nx
      tem7(i,j)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
!        tem6(i,j,k) = qc(i,j,k) + qr(i,j,k) + qi(i,j,k) +               &
!                      qs(i,j,k) + qh(i,j,k)
        DO nq=1,nscalarq
          tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,nq)
        END DO
        tem7(i,j)   = tem7(i,j) + tem6(i,j,k)*rhobar(i,j,k)             &
                                  *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_vic

