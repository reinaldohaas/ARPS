!
!##################################################################
!##################################################################
!######                                                      ######
!######            subroutine  linearint_2d                  ######
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
!  linear interpolation in 2D. 
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!-----------------------------------------------------------------------
!
!
SUBROUTINE  linearint_2df(nx,ny,vbl2,pxx,pyy,pval)
! 
  IMPLICIT NONE

  INTEGER :: nx, ny
  REAL    :: vbl2(nx,ny)
!  REAL    :: pxx, pyy
  DOUBLE PRECISION :: pxx, pyy
  REAL    :: pval
 
  INTEGER :: i, j
  REAL    :: deltadx,deltady,deltadxm,deltadym

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 
  i = FLOOR(pxx)
  j = FLOOR(pyy)
!
!    print*,'ij=',i,j,' pxx=',pxx,pyy
!
  IF((i > 0) .AND. (i < nx-1)  .AND. (j > 0) .AND. (j < ny-1)) THEN

    deltadx = pxx - i
    deltady = pyy - j

    deltadxm= 1. - deltadx
    deltadym= 1. - deltady

    pval =    deltadxm*deltadym * vbl2(i,  j  )                         &
            + deltadx *deltadym * vbl2(i+1,j  )                         &
            + deltadxm*deltady  * vbl2(i,  j+1)                         &
            + deltadx *deltady  * vbl2(i+1,j+1)

  ELSE

    WRITE (6,'(a)')          ' WARNING: '
    WRITE (6,'(2(a,f10.2))') ' pxx = ',pxx,' pyy = ',pyy
    WRITE (6,'(a,/)')        ' no interpolation was performed'

  END IF

  RETURN
END SUBROUTINE  linearint_2df
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               subroutine alinearint_2d               ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.                           ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE alinearint_2df(nx, ny, vbl2, pxx, pyy, pval)
!
  IMPLICIT NONE

  INTEGER :: nx, ny
  REAL    :: vbl2 (nx,ny)
!  REAL    :: pxx, pyy
  DOUBLE PRECISION :: pxx, pyy
  REAL    :: pval

  INTEGER :: i, j
  REAL    :: deltadx,deltady,deltadxm,deltadym

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  i = FLOOR (pxx)
  j = FLOOR (pyy)

  IF ((i > 0) .AND. (i < nx-1)  .AND. (j > 0) .AND. (j < ny-1)) THEN

    deltadx = pxx - i
    deltady = pyy - j

    deltadxm= 1.  - deltadx
    deltadym= 1.  - deltady

    vbl2(i+1,j+1)=vbl2(i+1,j+1) + deltadx*deltady *pval
    vbl2(i  ,j+1)=vbl2(i  ,j+1) + deltadxm*deltady*pval
    vbl2(i+1,j  )=vbl2(i+1,j  ) + deltadx*deltadym*pval
    vbl2(i  ,j  )=vbl2(i  ,j  ) + deltadxm*deltadym*pval
    pval = 0.

  ELSE
    WRITE (0,*) 'WARNING: no interpolation was performed'
  END IF

  RETURN
END SUBROUTINE alinearint_2df
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               subroutine  map_to_mod2                ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.                           ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE map_to_mod2(nx,ny,mxobs,nobs,useobs,pgx,pgy,px,py,pxx,pyy)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny, nobs, mxobs
  LOGICAL, INTENT(IN)  :: useobs(mxobs)
  REAL,    INTENT(IN)  :: pgx(nx), pgy(ny)
  REAL,    INTENT(IN)  :: px(mxobs), py(mxobs)
!  REAL,    INTENT(OUT) :: pxx(mxobs),pyy(mxobs)
  DOUBLE PRECISION,    INTENT(OUT) :: pxx(mxobs),pyy(mxobs)

  INTEGER :: i, j, n
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  pxx = -99999.
  pyy = -99999.

  obs_loop: DO n = 1, nobs

    IF ( useobs(n) ) THEN

      DO  j=1,ny-1
        DO  i=1,nx-1
          IF( (px(n) >= pgx(i)) .AND. (px(n) < pgx(i+1)) .AND.            &
              (py(n) >= pgy(j)) .AND. (py(n) < pgy(j+1)) ) THEN

!             pxx(n) = FLOAT(i)+ ( px(n)-pgx(i) )/( pgx(i+1)-pgx(i) )
!             pyy(n) = FLOAT(j)+ ( py(n)-pgy(j) )/( pgy(j+1)-pgy(j) )
             pxx(n) = i + ( px(n)-pgx(i) )/( pgx(i+1)-pgx(i) )
             pyy(n) = j + ( py(n)-pgy(j) )/( pgy(j+1)-pgy(j) )

             CYCLE obs_loop
          END IF
        END DO
      END DO

    END IF
  END DO obs_loop

  RETURN
END SUBROUTINE map_to_mod2


SUBROUTINE map_to_modz(nzk,mxobs,nlev,nobs,nx,ny,nz,pgz,                &
                       pxx, pyy, hgt, ihgt, pz1, pz2)  

  IMPLICIT NONE 

  INTEGER :: nzk,mxobs,nobs 
  INTEGER :: nx,ny,nz
  INTEGER :: nlev(mxobs)
  REAL    :: pgz(nx,ny,nz)
!  REAL    :: pxx(mxobs), pyy(mxobs)
  DOUBLE PRECISION    :: pxx(mxobs), pyy(mxobs)
  REAL    :: pz1(nzk,mxobs),pz2(nzk,mxobs)
  REAL    :: hgt(nzk,mxobs)
  INTEGER :: ihgt(nzk,mxobs)

  INTEGER :: k,ii,kk
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ihgt = -1
  obs_loop: DO ii = 1,nobs
  obslvl_loop:  DO kk = 1, nlev(ii)

      IF( pxx(ii)<-99990.0 .OR. pyy(ii)<-99990.0) CYCLE obs_loop

      DO k = 1, nz-1
        CALL linearint_2df(nx,ny,pgz(1,1,k  ),                          &
                           pxx(ii),pyy(ii),pz1(kk,ii) )
        CALL linearint_2df(nx,ny,pgz(1,1,k+1),                          &
                           pxx(ii),pyy(ii),pz2(kk,ii) )

        IF(  hgt(kk,ii) <= pz1(kk,ii) ) THEN 
          ihgt(kk,ii) = 0
          CYCLE obslvl_loop 
        ELSE IF( (hgt(kk,ii) > pz1(kk,ii) ) .AND.                       &
                 (hgt(kk,ii)<= pz2(kk,ii))) THEN
          ihgt(kk,ii) = k
          CYCLE obslvl_loop 
        END IF
      END DO

    END DO obslvl_loop
  END DO obs_loop

  RETURN
END SUBROUTINE map_to_modz
