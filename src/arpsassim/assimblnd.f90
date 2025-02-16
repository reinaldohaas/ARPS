!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASSIMBLND                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimblnd(nx,ny,nz,                                          &
           u_obs,u_bkg,hole_flag,                                       &
           bkg_err,obs_err,wgt,u_blnd)
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  A optimal weights is determined from the OI method, which uses
!  the error information of background fields and observations.
!
!  Note: The error profiles are subjectively given for the current
!        version. A more sophisticated method should be used in future
!        to objectively determine the weights at each vertical level
!        or even at each grid point in future.
!
!-------------------------------------------------------------------------
!
!  AUTH0R: Limin Zhao
!  02/15/1996
!
!  MODIFICATION HISTORY:
!
!--------------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: i,j,k
  REAL :: u_obs(nx,ny,nz)
  REAL :: u_bkg(nx,ny,nz)
  REAL :: hole_flag(nx,ny,nz)
  REAL :: bkg_err(nx,ny,nz)
  REAL :: obs_err(nx,ny,nz)
  REAL :: wgt(nx,ny,nz)
  REAL :: u_blnd(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Model control constants
  INCLUDE 'assim.inc'       ! Assim/Retr control parameters
  INCLUDE 'adas.inc'
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
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
!  Assign the error for each observations. Seperate the hole-filled
!  value from non hole-filled value.
!
!-----------------------------------------------------------------------
!
  WRITE(6,*)'code in assimwgt',vfill_err,v_err,adas_err

  IF (ip_wgt == 1) THEN

!   DO 100 k=1,nz-1
!   DO 100 j=1,ny-1
!   DO 100 i=1,nx-1
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          IF (hole_flag(i,j,k) == spval) THEN
            obs_err(i,j,k) = vfill_err
          ELSE
            obs_err(i,j,k) = v_err
          END IF

        END DO
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Assign the error for each background values.
!
!  Note:  It is assumed constant in Spring Tests, 1996.
!
!-----------------------------------------------------------------------
!ccxxx
!ccxxx
!ccxxx
!
!   DO 150 k=1,nz-1
!   DO 150 j=1,ny-1
!   DO 150 i=1,nx-1
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx

        bkg_err(i,j,k) = adas_err

      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate the optimal weight from background error and observation
!  errors.
!
!-----------------------------------------------------------------------
!
!   DO 200 k=1,nz-1
!   DO 200 j=1,ny-1
!   DO 200 i=1,nx-1
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        wgt(i,j,k) = bkg_err(i,j,k)*bkg_err(i,j,k)                      &
                     /(bkg_err(i,j,k)*bkg_err(i,j,k)                    &
                      +obs_err(i,j,k)*obs_err(i,j,k))
      END DO
    END DO
  END DO

!   DO 210 k=1,nz-1
!   DO 210 j=1,ny-1
!   DO 210 i=1,nx-1
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        u_blnd(i,j,k) = u_bkg(i,j,k)                                    &
                        + wgt(i,j,k)*(u_obs(i,j,k) - u_bkg(i,j,k))

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE assimblnd


SUBROUTINE smth(a,m,n,tem)

! Two-Dimension smooth use nine-points
  IMPLICIT NONE
  INTEGER :: m,n,i,j
  REAL :: s
  REAL :: a(m,n)
  REAL :: tem(m,n)
  s=0.5
  DO j=2,n-1
    DO i=2,m-1
      tem(i,j)=a(i,j)*(1.-s)**2 +                                       &
          ( a(i+1,j)+a(i-1,j)+a(i,j+1)+a(i,j-1) )*0.5*s*(1.-s)          &
          + ( a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1) )             &
          *0.25*s**2

    END DO
  END DO

  DO j=2,n-1
    tem(1,j) = tem(2,j)
    tem(m,j) = tem(m-1,j)
  END DO

  DO i=2,m-1
    tem(i,1) = tem(i,2)
    tem(i,n) = tem(i,n-1)
  END DO

  tem(1,1) = a(1,1)
  tem(1,n) = a(1,n)
  tem(m,1) = a(m,1)
  tem(m,n) = a(m,n)


  DO j=1,n
    DO i=1,m
      a(i,j) = tem(i,j)
    END DO
  END DO


  RETURN
END SUBROUTINE smth

