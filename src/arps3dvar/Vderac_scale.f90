!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SCALE_FACTOR              ######
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
!
!    define a dirac function at each grid point and
!    get the response function
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, 2000
!
!
!-----------------------------------------------------------------------
!
SUBROUTINE scale_factor(nx,ny,nz,pscalc,ipass_filt,radius,radius_z,dirac,    &
                        tem1,tem2,tem3,tem4)
!
!
! define a dirac function at each grid point and
! get the response function
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz, ipass_filt
  REAL    :: radius
  INTEGER :: ii,jj,kk,i,j,k
  REAL :: pscalc(nx,ny,nz)
  REAL :: dirac(nx,ny,nz)
  REAL :: tem1 (nx,ny,nz)
  REAL :: tem2 (nx,ny,nz)
  REAL :: tem3 (nx,ny,nz)
  REAL :: tem4 (nx,ny,nz)
  REAL :: radius_z(nx,ny,nz)

  REAL :: const

  INCLUDE 'mp.inc'
  INTEGER :: nxlg, nylg
  INTEGER :: imid, jmid, imidproc, jmidproc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dirac(:,:,:) = 0.

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
  jmid = nylg/2
  imid = nxlg/2
  imidproc = (imid-2) / (nx-3) + 1
  jmidproc = (jmid-2) / (ny-3) + 1

  IF (loc_x == imidproc) THEN
    imid = MOD((imid-2),(nx-3)) + 2   ! Local index for global middle point
  ELSE
    imid = -999
  END IF

  IF (loc_y == jmidproc) THEN
    jmid = MOD((jmid-2),(ny-3)) + 2   ! Local index for global middle point
  ELSE
    jmid = -999
  END IF

  kk = nz/2

  IF (loc_x == imidproc .AND. loc_y == jmidproc) THEN
    dirac(imid,jmid,kk) = 1.
  END IF

!  if (loc_x ==imidproc .and. loc_y==jmidproc) &
!write(0,*) 'imid,jmid =',imid,jmid,imidproc, jmidproc

  CALL  recurfilt3d( nx,ny,nz,dirac,ipass_filt,ipass_filt,              &
                     radius,radius_z,tem1,tem2,tem3,tem4)
!  const = sqrt( 1. / dirac(ii,jj,kk))

  IF (loc_x == imidproc .AND. loc_y == jmidproc) THEN
    const = SQRT( 1. / dirac(imid,jmid,kk))
    PRINT *,'pscalc = ',const
  END IF

  CALL mpbcastr(const,(jmidproc-1)*nproc_x+imidproc-1)

!  write(0,*)'pscalc(ii,jj,kk) =',pscalc(ii,jj,kk),dirac(ii,jj,kk)

  DO kk=1, nz
    DO jj=1, ny
      DO ii=1, nx
        pscalc (ii,jj,kk) = const
      END DO
    END DO
  END DO

!  stop
!
! DO kk=1, nz
!   DO jj=1, ny
!     DO ii=1, nx
!
!         dirac = 0.
!
!     dirac(ii,jj,kk) = 1.
!
!     CALL  recurfilt_2d( nx,ny,dirac,ipass_filt,radius )
!
!     CALL  recurfilt_3d( nx,ny,nz,dirac,ipass_filt,radius,radius_z )
!
!     pscalc (ii,jj,kk) = 1. / dirac(ii,jj,kk)
!
!     END DO
!   END DO
! END DO
!
  RETURN
END SUBROUTINE scale_factor
