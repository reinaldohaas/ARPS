!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CTR_TO_VBL                ######
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
!  PURPOSE: Transfer from control varibles to analysis varibles.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jidong GAO, CAPS, June, 2000
!
!-----------------------------------------------------------------------
!
SUBROUTINE ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,             &
               pbkg,pscal,pgrd,tem1,tem2,tem3,tem4)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ipass_filt
  REAL,    INTENT(IN) :: hradius
  INTEGER, INTENT(IN) :: nx, ny, nz

  REAL, INTENT(IN)    :: radius_z(nx,ny,nz)
  REAL, INTENT(IN)    :: pbkg (nx,ny,nz)
  REAL, INTENT(IN)    :: pscal(nx,ny,nz)

  REAL, INTENT(INOUT) ::  pgrd(nx,ny,nz)

  REAL    :: tem1(nx,ny,nz)
  REAL    :: tem2(nx,ny,nz)
  REAL    :: tem3(nx,ny,nz)
  REAL    :: tem4(nx,ny,nz)

  INTEGER :: i, j, k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL recurfilt3d(nx,ny,nz,pgrd,ipass_filt,ipass_filt/2,               &
                   hradius,radius_z,tem1,tem2,tem3,tem4)

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        pgrd (i,j,k) = pgrd(i,j,k) * pscal(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        pgrd (i,j,k) =  pgrd(i,j,k)*pbkg(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE ctr_to_vbl
