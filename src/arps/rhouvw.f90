!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RHOUVW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rhouvw(nx,ny,nz,rhostr,rhostru,rhostrv,rhostrw)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate rhostr averaged to u, v, and w points.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  3/8/1993.
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    rhostr   j3 times base state density rhobar(kg/m**3).
!
!  OUTPUT:
!
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
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

  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions

  REAL :: rhostr(nx,ny,nz)     ! j3 times base state density rhobar
                               ! (kg/m**3).

  REAL :: rhostru(nx,ny,nz)    ! Average rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Average rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Average rhostr at w points (kg/m**3).

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL avgsu(rhostr,nx,ny,nz, 1,ny-1, 1,nz-1, rhostru, rhostrw) 
    ! rhostrw used here as a temporary array
  CALL avgsv(rhostr,nx,ny,nz, 1,nx-1, 1,nz-1, rhostrv, rhostrw)
    ! rhostrw used here as a temporary array
  CALL avgsw(rhostr,nx,ny,nz, 1,nx-1, 1,ny-1, rhostrw)

  RETURN
END SUBROUTINE rhouvw
