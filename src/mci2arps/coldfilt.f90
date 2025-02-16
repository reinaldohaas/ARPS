  SUBROUTINE coldfilt(nx,ny,filtlen,t11mu,t11mu_cold)
!-----------------------------------------------------------------------
!
!  PURPOSE: Calculate cold filtered temperatures
!           Adapted from cldinsert program for use in mci2arps.
!
!  AUTHOR: Keith Brewster, CAPS
!          27 MAR 2007
!
!  ARGUMENTS
!    INPUT
!    nx, ny      : x and y dimensions
!    filtlen     : filter distance (m)
!    dx, dy      : x and y grid spacing
!    t11mu       : 10.7 micron IR temperature (K)
!
!    OUTPUT
!    t11mu_cold  : cold-filtered 10.7 micron IR temperature (K)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  REAL, INTENT(IN) :: filtlen
  REAL, INTENT(IN) :: t11mu(nx,ny)
  REAL, INTENT(OUT) :: t11mu_cold(nx,ny)
!
! Misc local variables
!
  INTEGER :: idelt,jdelt
  INTEGER :: ibgn,iend,jbgn,jend
  INTEGER :: i,j,ii,jj
  REAL :: tmin,tmax
!
! Include files
!
  INCLUDE 'grid.inc'
!
  t11mu_cold = 999.
  idelt = max(1,nint(filtlen/dx))
  jdelt = max(1,nint(filtlen/dy))
  DO j = 1,ny
    DO i = 1,nx
      jbgn = max(1 ,j-jdelt)
      jend = min(ny,j+jdelt)
      ibgn = max(1 ,i-idelt)
      iend = min(nx,i+idelt)

      DO jj = jbgn,jend
        DO ii = ibgn,iend
          IF (t11mu(ii,jj) > 0.)                                     &
              t11mu_cold(i,j) = min(t11mu_cold(i,j),t11mu(ii,jj))
        END DO ! ii
      END DO ! jj

      IF(t11mu_cold(i,j) > 350.) t11mu_cold(i,j) = -999.

    END DO ! i
  END DO ! j

  tmin=999.
  tmax=-999.
  DO j=1,ny
  DO i=1,nx
    IF(t11mu(i,j) > 0.) THEN
      tmin=min(tmin,t11mu(i,j))
      tmax=max(tmax,t11mu(i,j))
    END IF
  END DO 
  END DO
  WRITE(6,'(/a,f9.2,a,f9.2)') ' IR Temps Min: ',tmin,'   Max:',tmax

  tmin=999.
  tmax=-999.
  DO j=1,ny
  DO i=1,nx
    IF(t11mu_cold(i,j) > 0.) THEN
      tmin=min(tmin,t11mu_cold(i,j))
      tmax=max(tmax,t11mu_cold(i,j))
    END IF
  END DO 
  END DO
  WRITE(6,'(a,f9.2,a,f9.2/)') ' IR Cold Filt T Min: ',tmin,'   Max:',tmax

  RETURN
END SUBROUTINE COLDFILT
