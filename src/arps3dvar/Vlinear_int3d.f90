!
!##################################################################
!##################################################################
!######                                                      ######
!######       The Tri_linear Interpolation subroutine        ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.                           ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Linear Interpolation in 3D.                         
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!-----------------------------------------------------------------------
!
SUBROUTINE  linearint_3d(nx,ny,nz,vbl3,pgx,pgy,pgz,ownobs,              &
                         nzk,mxobs, nlev, nobs, ivar,                   &
                         pxx,pyy,pz1,pz2,hgt,ihgt,pval )
!
!ivar, px,py,pz,pval,iflag)
! 
!  nx,ny,nz:    Dimension of 3d field 
!  vbl3:        3d field
!  pgx,pgy,pgz: Coordinate of 3d field 
!  ownobs:      Whether this processor owns each observations
!  ivar :       varible type, indicate position
!  pxx,pyx,pz :   observation point position
!  pval :       return value at above point
! 
  IMPLICIT NONE 

  INTEGER :: nx, ny, nz, ivar
  INTEGER :: nzk, mxobs, nobs
  INTEGER :: ii,kk,ik
  REAL    :: vbl3(nx,ny,nz)
  REAL    :: pgx(nx)
  REAL    :: pgy(ny)
  REAL    :: pgz(nx,ny,nz)
  LOGICAL, INTENT(IN) :: ownobs(mxobs)

  INTEGER :: nlev(mxobs)
!  REAL    :: pxx(mxobs), pyy(mxobs)
  DOUBLE PRECISION    :: pxx(mxobs), pyy(mxobs)
  REAL    :: pz1(nzk,mxobs)
  REAL    :: pz2(nzk,mxobs)
  REAL    :: hgt(nzk,mxobs)
  INTEGER :: ihgt(nzk,mxobs)
  REAL    :: pval(nzk,mxobs)

  REAL :: zv1,zv2,zdz2,zdz1
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO ii = 1,nobs
    IF ( ownobs(ii) ) THEN
      DO kk = 1, nlev(ii)

        ik = ihgt(kk,ii)
        IF( ik>0 ) THEN

          CALL linearint_2df(nx,ny,vbl3(1,1,ik  ),pxx(ii),pyy(ii),zv1)
          CALL linearint_2df(nx,ny,vbl3(1,1,ik+1),pxx(ii),pyy(ii),zv2)

          zdz2 =  (hgt(kk,ii)-pz1(kk,ii))/(pz2(kk,ii)-pz1(kk,ii))
          zdz1 = -(hgt(kk,ii)-pz2(kk,ii))/(pz2(kk,ii)-pz1(kk,ii))

          pval(kk,ii) = zdz1 * zv1 + zdz2 * zv2

        ELSE IF( ik == 0) THEN

          CALL linearint_2df(nx,ny,vbl3(1,1,ik+1),pxx(ii),pyy(ii),zv2)

          IF( (ivar == 1) .OR.(ivar == 2).OR.(ivar == 4) .OR.              &
              (ivar == 5) .OR.(ivar == 6) ) THEN
            pval(kk,ii) = zv2
          ELSE IF(ivar == 3 ) THEN
            pval(kk,ii) = zv2 + 1.145E-4*(pz2(kk,ii)-hgt(kk,ii))
          ELSE
            PRINT*,' stop, the analysis variable is not exist!'
            CALL arpsstop('The analysis variable does not exist.',1)
          END IF
        !ELSE
          !PRINT*,'the observation is out of top of model domain!!!!!!'
        END IF
        !print*,'pval============',pval(kk,ii)
      END DO

    END IF   ! ownobs
  END DO

  RETURN
END SUBROUTINE linearint_3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######       The Tri_linear Interpolation subroutine        ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.                           ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE  adlinearint_3d(nx,ny,nz,vbl3,pgx,pgy,pgz,useobs,            &
                          nzk,mxobs, nlev, nobs, ivar,                  &
                          pxx,pyy,pz1,pz2,hgt,ihgt,pval, tem1 )
!
!ivar, px,py,pz,pval,iflag)
! 
!  nx,ny,nz:    Dimension of 3d field 
!  pgx,pgy,pgz: Coordinate of 3d field 
!  vbl3:        3d field
!  ivar :       varible type, indicate position
!  px,py,pz :   observation point position
!  pval :       return value at above point
! 
  IMPLICIT NONE

  INTEGER :: nx, ny, nz, ivar
  INTEGER :: nzk, mxobs, nobs
  INTEGER :: ii,kk,ik
  REAL    :: vbl3(nx,ny,nz)
  REAL    :: pgx(nx)
  REAL    :: pgy(ny)
  REAL    :: pgz(nx,ny,nz)
  LOGICAL, INTENT(IN) :: useobs(mxobs)   ! Observation will affect this subdomain

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nlev(mxobs)
!  REAL    :: pxx(mxobs), pyy(mxobs)
  DOUBLE PRECISION    :: pxx(mxobs), pyy(mxobs)
  REAL    :: pz1(nzk,mxobs)
  REAL    :: pz2(nzk,mxobs)
  REAL    :: hgt(nzk,mxobs)
  INTEGER :: ihgt(nzk,mxobs)
  REAL    :: pval(nzk,mxobs)

  REAL    :: zv1,zv2,zdz2,zdz1

  REAL    :: tem1(nx,ny,nz)              ! working array

  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO ii = 1,nobs

    IF ( useobs(ii) ) THEN

      DO kk = 1, nlev(ii)

        ik = ihgt(kk,ii)
        IF( ik>0 ) THEN

          zdz2 =  (hgt(kk,ii)-pz1(kk,ii))/(pz2(kk,ii)-pz1(kk,ii))
          zdz1 = -(hgt(kk,ii)-pz2(kk,ii))/(pz2(kk,ii)-pz1(kk,ii))

          zv2 =  zdz2 * pval(kk,ii) 
          zv1 =  zdz1 * pval(kk,ii) 
          pval(kk,ii) = 0.

          CALL alinearint_2df(nx,ny,vbl3(1,1,ik  ),pxx(ii),pyy(ii),zv1)

          CALL alinearint_2df(nx,ny,vbl3(1,1,ik+1),pxx(ii),pyy(ii),zv2)

        ELSE IF( ik==0 ) THEN

          IF( (ivar == 1).OR.(ivar == 2).OR.(ivar == 4) .OR.            &
              (ivar == 5).OR.(ivar == 6) ) THEN

            zv2 =  pval(kk,ii)
            pval(kk,ii) = 0.0

          ELSE IF(ivar == 3 ) THEN
            zv2 =  pval(kk,ii)
            pval(kk,ii)=0.0

          ELSE
            PRINT*,' stop, the analysis variable is not exist!'
            CALL arpsstop('The analysis variable does not exist.',1)
          END IF

          CALL alinearint_2df(nx,ny,vbl3(1,1,ik+1),pxx(ii),pyy(ii),zv2)
        !ELSE
        !  PRINT*,'the observation is out of top of model domain!!!!!!'
        END IF
      END DO

    END IF            ! useobs
  END DO

  IF (mp_opt > 0) THEN    ! All arrays are on mass grid
!    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(vbl3, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(vbl3, nx, ny, nz, nbc, sbc, 0, tem1)
!    CALL acct_stop_inter
  END IF

  RETURN
END SUBROUTINE adlinearint_3d
