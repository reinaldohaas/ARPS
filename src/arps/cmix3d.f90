!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CMIX2UVW                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cmix2uvw(nx,ny,nz,                                           &
           u,v,w, ubar,vbar,rhostr,                                     &
           umix,vmix,wmix,                                              &
           tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the second order computational mixing terms for the momentum
!  equations. Computational mixing is applied to velocity perturbations
!  only. These terms are placed in the arrays umix, vmix and wmix.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/3/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  4/2/92 (M. Xue and H. Jin)
!  Modify to include terrain
!
!  5/19/1998 (M. Xue)
!  Reformulated the computational mixing terms to move rhostr
!  outside the inner-most derivative.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Verticle component of Cartesian velocity at a
!             given time level (m/s)
!
!    umix     Array containing turbulent mixing on u
!    vmix     Array containing turbulent mixing on v
!    wmix     Array containing turbulent mixing on w
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  OUTPUT:
!
!    umix     Turbulent and computational mixing on u.
!    vmix     Turbulent and computational mixing on v.
!    wmix     Turbulent and computational mixing on w.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: umix  (nx,ny,nz)     ! Total mixing in u-eq.
  REAL :: vmix  (nx,ny,nz)     ! Total mixing in v-eq.
  REAL :: wmix  (nx,ny,nz)     ! Total mixing in w-eq.

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: temx,temy,temz
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  IF the coefficients of computational mixing in both horizontal
!  and vertical directions are zero, exit this subroutine.
!
!-----------------------------------------------------------------------
!
  IF( cfcmh2 == 0.0 .AND. cfcmv2 == 0.0 ) RETURN

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        tem3(i,j,k)=u(i,j,k)-ubar(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the second order
!  horizontal computational mixing on u'
!
!-----------------------------------------------------------------------
!
  IF( cfcmh2 /= 0.0) THEN

    temx = cfcmh2/(dx*dx)
    temy = cfcmh2/(4.0*dy*dy)

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(tem3(i+1,j,k)-tem3(i,j,k))*rhostr(i,j,k)*temx
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=2,nx-1
          tem2(i,j,k)=(tem3(i,j,k)-tem3(i,j-1,k))*                      &
              (rhostr(i,j  ,k)+rhostr(i-1,j  ,k)                        &
              +rhostr(i,j-1,k)+rhostr(i-1,j-1,k))*temy
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO i=2,nx-1
        tem2(i, 1,k) = 0.0
        tem2(i,ny,k) = 0.0
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          umix(i,j,k)=umix(i,j,k)+(tem1(i,j,k)-tem1(i-1,j,k)            &
                                  +tem2(i,j+1,k)-tem2(i,j,k))
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the second order
!  vertical computational mixing on u'
!
!-----------------------------------------------------------------------
!
  IF( cfcmv2 /= 0.0) THEN
!
    temz = cfcmv2/(4.0*dz*dz)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem2(i,j,k)=(tem3(i,j,k)-tem3(i,j,k-1))*                      &
              (rhostr(i,j,k  )+rhostr(i-1,j,k  )                        &
              +rhostr(i,j,k-1)+rhostr(i-1,j,k-1))*temz
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          umix(i,j,k)=umix(i,j,k)+(tem2(i,j,k+1)-tem2(i,j,k))
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the second order
!  horizontal computational mixing on v'
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        tem3(i,j,k)=v(i,j,k)-vbar(i,j,k)
      END DO
    END DO
  END DO

  IF( cfcmh2 /= 0.0) THEN

    temx = cfcmh2/(4.0*dx*dx)
    temy = cfcmh2/(dy*dy)

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(tem3(i,j,k)-tem3(i-1,j,k))*                      &
              (rhostr(i,j  ,k)+rhostr(i  ,j-1,k)                        &
              +rhostr(i-1,j,k)+rhostr(i-1,j-1,k))*temx
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        tem1(1, j,k) = 0.0
        tem1(nx,j,k) = 0.0
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(tem3(i,j+1,k)-tem3(i,j,k))*rhostr(i,j,k)*temy
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          vmix(i,j,k)=vmix(i,j,k)+(tem1(i+1,j,k)-tem1(i,j,k)            &
                                  +tem2(i,j,k)-tem2(i,j-1,k))
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the second order
!  vertical computational mixing on v'
!
!-----------------------------------------------------------------------
!
  IF( cfcmv2 /= 0.0) THEN

    temz = cfcmv2/(4.0*dz*dz)

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(tem3(i,j,k)-tem3(i,j,k-1))*                      &
              (rhostr(i,j,k  )+rhostr(i,j-1,k  )                        &
              +rhostr(i,j,k-1)+rhostr(i,j-1,k-1))*temz
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          vmix(i,j,k)=vmix(i,j,k)+(tem2(i,j,k+1)-tem2(i,j,k))
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Second order computational mixing on w
!
!-----------------------------------------------------------------------
!
  IF( cfcmh2 /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the horizontal
!  computational mixing on w.
!
!-----------------------------------------------------------------------
!
    temx = cfcmh2/(4.0*dx*dx)
    temy = cfcmh2/(4.0*dy*dy)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(w(i,j,k)-w(i-1,j,k))                             &
                     *(rhostr(i,j  ,k)+rhostr(i-1,j  ,k)+               &
                       rhostr(i,j,k-1)+rhostr(i-1,j,k-1))*temx
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO j=1,ny-1
        tem1(1 ,j,k)=0.0
        tem1(nx,j,k)=0.0
      END DO
    END DO

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(w(i,j,k)-w(i,j-1,k))                             &
                     *(rhostr(i,j,  k)+rhostr(i,j-1,k)+                 &
                       rhostr(i,j,k-1)+rhostr(i,j-1,k-1))*temy
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO i=1,nx-1
        tem2(i,1 ,k)=0.0
        tem2(i,ny,k)=0.0
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Note that tem1 at i=1 and nx are zero, and tem2 at j=1,ny are zero,
!  equivalent to assuming zero normal gradient of s at the boundaries.
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wmix(i,j,k)=wmix(i,j,k)+(tem1(i+1,j,k)-tem1(i,j,k)+           &
                                   tem2(i,j+1,k)-tem2(i,j,k))
        END DO
      END DO
    END DO

  END IF

  IF( cfcmv2 /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the vertical
!  computational mixing on w.
!
!-----------------------------------------------------------------------
!
    temz = cfcmv2/(dz*dz)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(w(i,j,k+1)-w(i,j,k))*rhostr(i,j,k)*temz
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wmix(i,j,k)=wmix(i,j,k)+(tem1(i,j,k)-tem1(i,j,k-1))
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE cmix2uvw

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CMIX2S                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cmix2s(nx,ny,nz, s ,rhostr, smix, tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generic routine to calculate the second order computational mixing
!  for scalar s. The computational mixing is accumulated into array
!  smix on exit.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/3/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  10/17/93 (M. Xue)
!  Bug fixes fix in line 3 of loop 111 and 112.
!
!  10/17/94 (M. Xue)
!  Subroutines CMIX2PT and CMIX2Q combined into a single generic
!  routine CMIX2S.
!
!  5/19/1998 (M. Xue)
!  Reformulated the computational mixing terms to move rhostr
!  outside the inner-most derivative.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Verticle component of Cartesian velocity at a
!             given time level (m/s)
!
!    s        scalar field to which the compuational mixing is to
!             be applied.
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    smix     Array containing turbulent mixing s on input.
!
!  OUTPUT:
!
!    smix     Turbulent mixing plus computational mixing on s.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s     (nx,ny,nz)     ! Scalar to be mixed.

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: smix  (nx,ny,nz)     ! Total mixing on s

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: cfcmh2s,cfcmv2s, temx,temy,temz
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'phycst.inc'
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
!
!-----------------------------------------------------------------------
!
!  If the coefficients of computational mixing in both the horizontal
!  and vertical are zero, exit this subroutine.
!
!-----------------------------------------------------------------------
!
  cfcmh2s = cfcmh2 * scmixfctr
  cfcmv2s = cfcmv2 * scmixfctr

  IF( cfcmh2s == 0.0 .AND. cfcmv2s == 0.0 ) RETURN
!
!-----------------------------------------------------------------------
!
!  Second order computational mixing on s.
!
!-----------------------------------------------------------------------
!
  IF( cfcmh2s /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the horizontal
!  computational mixing.
!
!-----------------------------------------------------------------------
!
    temx = cfcmh2s/(2.0*dx*dx)
    temy = cfcmh2s/(2.0*dy*dy)


    DO k=1,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(s(i,j,k)-s(i-1,j,k))                             &
                     *(rhostr(i,j,k)+rhostr(i-1,j,k))*temx
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        tem1(1 ,j,k)=0.0
        tem1(nx,j,k)=0.0
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(s(i,j,k)-s(i,j-1,k))                             &
                     *(rhostr(i,j,k)+rhostr(i,j-1,k))*temy
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO i=1,nx-1
        tem2(i,1 ,k)=0.0
        tem2(i,ny,k)=0.0
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Note that tem1 at i=1 and nx are zero, and tem2 at j=1,ny are zero,
!  equivalent to assuming zero normal gradient of s at the boundaries.
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          smix(i,j,k)=smix(i,j,k)+(tem1(i+1,j,k)-tem1(i,j,k)+           &
                                   tem2(i,j+1,k)-tem2(i,j,k))
        END DO
      END DO
    END DO

  END IF

  IF( cfcmv2s /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the vertical
!  computational mixing.
!
!-----------------------------------------------------------------------
!
    temz = cfcmv2s/(2.0*dz*dz)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(s(i,j,k)-s(i,j,k-1))                             &
                     *(rhostr(i,j,k)+rhostr(i,j,k-1))*temz
        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,1 )=0.0
        tem1(i,j,nz)=0.0
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          smix(i,j,k)=smix(i,j,k)+(tem1(i,j,k+1)-tem1(i,j,k))
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE cmix2s
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CMIX4UVW                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cmix4uvw(nx,ny,nz,                                           &
           u,v,w, ubar,vbar,rhostr,                                     &
           umix,vmix,wmix,                                              &
           tem1,tem2,tem3,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the fourth order computational mixing terms for the momentum
!  equations. These terms are placed in the arrays umix, vmix and wmix
!  and are applied to momentum perturbations only.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/3/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  4/2/92 (M. Xue and H. Jin)
!  Modify to include terrain
!
!  5/25/1998 (M. Xue)
!  Reformulated the computational mixing terms to move rhostr
!  outside the inner-most derivative.
!
!  11/9/2001 (M. Xue, D. Weber, and X. Jin)
!  Added monotonic computational mixing and 6-th order mixing option
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Verticle component of Cartesian velocity at a
!             given time level (m/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    umix     Array containing turbulent mixing on u (kg/(m*s)**2)
!    vmix     Array containing turbulent mixing on v (kg/(m*s)**2)
!    wmix     Array containing turbulent mixing on w (kg/(m*s)**2)
!
!  OUTPUT:
!
!    umix     Turbulent and computational mixing on u (kg/(m*s)**2)
!    vmix     Turbulent and computational mixing on v (kg/(m*s)**2)
!    wmix     Turbulent and computational mixing on w (kg/(m*s)**2)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: umix  (nx,ny,nz)     ! Total mixing in u-eq. (kg/(m*s)**2)
  REAL :: vmix  (nx,ny,nz)     ! Total mixing in v-eq. (kg/(m*s)**2)
  REAL :: wmix  (nx,ny,nz)     ! Total mixing in w-eq. (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: temx, temy, temz
  INTEGER :: kount
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------

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
!  If the coefficients of computational mixing in both the horizontal
!  and vertical are zero, exit this subroutine.
!
!-----------------------------------------------------------------------
!
  IF( cfcmh4 == 0.0 .AND. cfcmv4 == 0.0 ) RETURN

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        tem5(i,j,k)=u(i,j,k)-ubar(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate ustr=rhostr*u', vstr=rhostr*v', wstr=rhostr*w'
!
!-----------------------------------------------------------------------
!


  IF( cfcmh4 /= 0.0 ) THEN

    temx = cfcmh4/(    dx**4)
    temy = cfcmh4/(4.0*dy**4)

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(tem5(i+1,j,k)-tem5(i,j,k))*rhostr(i,j,k)*temx
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=2,nx-1
          tem2(i,j,k)=(tem5(i,j,k)-tem5(i,j-1,k))*                      &
              (rhostr(i,j  ,k)+rhostr(i-1,j  ,k)                        &
              +rhostr(i,j-1,k)+rhostr(i-1,j-1,k))*temy
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO i=2,nx-1
        tem2(i, 1,k) = 0.0
        tem2(i,ny,k) = 0.0
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          tem3(i,j,k)=tem1(i,j,k)-tem1(i-1,j,k)
          tem4(i,j,k)=tem2(i,j+1,k)-tem2(i,j,k)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem3,nx,ny,nz,ebc,wbc,1,tem1)
      CALL mpsendrecv2dns(tem4,nx,ny,nz,nbc,sbc,1,tem1)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL budifxx(tem3, nx,ny,nz,1,ny-1,2,nz-2,ebc,wbc)
    CALL budifyy(tem4, nx,ny,nz,2,nx-1,2,nz-2,nbc,sbc)
    CALL acct_stop_inter

    IF(cmix_opt == 0) THEN     !  non-monotonic 4th order comp. mixing

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)-(                                   &
                (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k)) &
                +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN  !  4th order monotonic computational mixing

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=1,nx-1
            tem1(i,j,k)=(tem3(i+1,j,k)-tem3(i,j,k))
            IF(-tem1(i,j,k)*(tem5(i+1,j,k)-tem5(i,j,k)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-1
            tem2(i,j,k)=(tem4(i,j,k)-tem4(i,j-1,k))
            IF(-tem2(i,j,k)*(tem5(i,j,k)-tem5(i,j-1,k)) < 0.0) tem2(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)-                                    &
                (tem1(i,j,k)-tem1(i-1,j,k)+tem2(i,j+1,k)-tem2(i,j,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-1
            tem1(i,j,k)=(tem3(i+1,j,k)-tem3(i,j,k))                     &
                       -(tem3(i,j,k)-tem3(i-1,j,k))
            tem2(i,j,k)=(tem4(i,j+1,k)-tem4(i,j,k))                     &
                       -(tem4(i,j,k)-tem4(i,j-1,k))
          END DO
        END DO
      END DO

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,1,tem3)
        CALL mpsendrecv2dns(tem2,nx,ny,nz,nbc,sbc,1,tem3)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL budifxx(tem1, nx,ny,nz,1,ny-1,2,nz-2,ebc,wbc)
      CALL budifyy(tem2, nx,ny,nz,2,nx-1,2,nz-2,nbc,sbc)
      CALL acct_stop_inter

      kount = 0
      DO k=2,nz-2
        DO j=2,ny-2
          DO i=1,nx-1
            tem3(i,j,k)=(tem1(i+1,j,k)-tem1(i,j,k))
            IF(cmix_opt == 3.AND.                                       &
                  tem3(i,j,k)*(tem5(i+1,j,k)-tem5(i,j,k)) < 0.0)THEN
              tem3(i,j,k)=0.0
              IF(j == 2) kount = kount + 1
            END IF
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-1
            tem4(i,j,k)=(tem2(i,j,k)-tem2(i,j-1,k))
            IF(cmix_opt == 3.AND.                                       &
                tem4(i,j,k)*(tem5(i,j,k)-tem5(i,j-1,k)) < 0.0)          &
                tem4(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)+                                    &
                (tem3(i,j,k)-tem3(i-1,j,k)+tem4(i,j+1,k)-tem4(i,j,k))
          END DO
        END DO
      END DO

    END IF

    !wdt update
    IF (sbc == 4) THEN ! do them for open BC only
      j = 1
      DO k=2,nz-2
        DO i=2,nx-1
          umix(i,j,k)=umix(i,j,k)-(                                     &
               (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))  &
              +(tem4(i,j+1,k)-tem4(i,j,k))- tem4(i,j,k))
        END DO
      END DO
    END IF
    IF (nbc == 4) THEN ! do them for open BC only
      j = ny-1
      DO k=2,nz-2
        DO i=2,nx-1
          umix(i,j,k)=umix(i,j,k)-(                                     &
               (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))  &
              +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)))
        END DO
      END DO
    END IF

  END IF

  IF( cfcmv4 /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the vertical
!  computational mixing on u'
!
!-----------------------------------------------------------------------
!
    temz = cfcmv4/(4.0*dz**4)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(tem5(i,j,k)-tem5(i,j,k-1))*                      &
              (rhostr(i,j,k  )+rhostr(i-1,j,k  )                        &
              +rhostr(i,j,k-1)+rhostr(i-1,j,k-1))*temz
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          tem2(i,j,k)=tem1(i,j,k+1)-tem1(i,j,k)
        END DO
      END DO
    END DO

    CALL budifzz(tem2, nx,ny,nz,2,nx-1,1,ny-1, tbc, bbc)

    IF( cmix_opt == 0) THEN  ! no monotonic application.

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)                                     &
                -((tem2(i,j,k+1)-tem2(i,j,k))-(tem2(i,j,k)-tem2(i,j,k-1)))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN !  fourth order monotonic 

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=2,nx-1
            tem1(i,j,k)=tem2(i,j,k)-tem2(i,j,k-1)
            IF(-tem1(i,j,k)*(tem5(i,j,k)-tem5(i,j,k-1)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)-(tem1(i,j,k+1)-tem1(i,j,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=2,nx-1
            tem3(i,j,k)= (tem2(i,j,k+1)-tem2(i,j,k))                    &
                        -(tem2(i,j,k)-tem2(i,j,k-1))
          END DO
        END DO
      END DO
      CALL budifzz(tem3, nx,ny,nz,2,nx-1,1,ny-1, tbc, bbc)

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=2,nx-1
            tem1(i,j,k)=tem3(i,j,k)-tem3(i,j,k-1)
            IF(cmix_opt == 3.AND.                                       &
                  tem1(i,j,k)*(tem5(i,j,k)-tem5(i,j,k-1)) < 0.0)THEN
              tem1(i,j,k)=0.0
              IF(j == 2)kount = kount+1
            END IF
          END DO
        END DO
      END DO

!      PRINT*,'On-off switch active for u at ',                          &
!              FLOAT(kount)/(2*(nx-2)*(nz-2))                            &
!             ,' percent of times at t=', curtim

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)+(tem1(i,j,k+1)-tem1(i,j,k))
          END DO
        END DO
      END DO

    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  Fourth-order computational mixing on v'.
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        tem5(i,j,k)=v(i,j,k)-vbar(i,j,k)
      END DO
    END DO
  END DO

  IF( cfcmh4 /= 0.0 ) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the horizontal
!  computational mixing on v'
!
!-----------------------------------------------------------------------
!

    temx = cfcmh4/(4.0*dx**4)
    temy = cfcmh4/(dy**4)

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(tem5(i,j,k)-tem5(i-1,j,k))*                      &
              (rhostr(i,j  ,k)+rhostr(i  ,j-1,k)                        &
              +rhostr(i-1,j,k)+rhostr(i-1,j-1,k))*temx
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        tem1(1, j,k) = 0.0
        tem1(nx,j,k) = 0.0
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(tem5(i,j+1,k)-tem5(i,j,k))*rhostr(i,j,k)*temy
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=tem1(i+1,j,k)-tem1(i,j,k)
          tem4(i,j,k)=tem2(i,j,k)-tem2(i,j-1,k)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem3,nx,ny,nz,ebc,wbc,2,tem1)
      CALL mpsendrecv2dns(tem4,nx,ny,nz,nbc,sbc,2,tem1)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bvdifxx(tem3, nx,ny,nz,2,ny-1,2,nz-2,ebc, wbc)
    CALL bvdifyy(tem4, nx,ny,nz,1,nx-1,2,nz-2,nbc, sbc)
    CALL acct_stop_inter

    IF( cmix_opt == 0) THEN  !  4th order non-monotonic computational mixing

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-2
            vmix(i,j,k)=vmix(i,j,k)-(                                   &
                (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k)) &
                +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN  ! 4th order monotonic

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-1
            tem1(i,j,k)=(tem3(i,j,k)-tem3(i-1,j,k))
            IF(-tem1(i,j,k)*(tem5(i,j,k)-tem5(i-1,j,k)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=2,nx-2
            tem2(i,j,k)=(tem4(i,j+1,k)-tem4(i,j,k))
            IF(-tem2(i,j,k)*(tem5(i,j+1,k)-tem5(i,j,k)) < 0.0) tem2(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-2
            vmix(i,j,k)=vmix(i,j,k)-                                    &
                 (tem1(i+1,j,k)-tem1(i,j,k)+tem2(i,j,k)-tem2(i,j-1,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-2
            tem1(i,j,k)=(tem3(i+1,j,k)-tem3(i,j,k))                     &
                       -(tem3(i,j,k)-tem3(i-1,j,k))
            tem2(i,j,k)=(tem4(i,j+1,k)-tem4(i,j,k))                     &
                       -(tem4(i,j,k)-tem4(i,j-1,k))
          END DO
        END DO
      END DO

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,2,tem3)
        CALL mpsendrecv2dns(tem2,nx,ny,nz,nbc,sbc,2,tem3)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bvdifxx(tem1, nx,ny,nz,2,ny-1,2,nz-2,ebc, wbc)
      CALL bvdifyy(tem2, nx,ny,nz,1,nx-1,2,nz-2,nbc, sbc)
      CALL acct_stop_inter

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-1
            tem3(i,j,k)=(tem1(i,j,k)-tem1(i-1,j,k))
            IF(cmix_opt == 3.AND.                                       &
                tem3(i,j,k)*(tem5(i,j,k)-tem5(i-1,j,k)) < 0.0)          &
                tem3(i,j,k)=0.0
          END DO
        END DO
      END DO

!      kount = 0

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=2,nx-2
            tem4(i,j,k)=(tem2(i,j+1,k)-tem2(i,j,k))
            IF(cmix_opt == 3.AND.                                       &
                  tem4(i,j,k)*(tem5(i,j+1,k)-tem5(i,j,k)) < 0.0)THEN
              tem4(i,j,k)=0.0
!           if(j.eq.2) kount = kount + 1
            END IF
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-2
            vmix(i,j,k)=vmix(i,j,k)+                                    &
                 (tem3(i+1,j,k)-tem3(i,j,k)+tem4(i,j,k)-tem4(i,j-1,k))
          END DO
        END DO
      END DO

    END IF

    !wdt update
    IF (ebc == 4) THEN ! Do them for open BC only
      i = nx-1
      DO k=2,nz-2
        DO j=2,ny-1
          vmix(i,j,k)=vmix(i,j,k)-(                                     &
              (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)))
        END DO
      END DO
    END IF
    IF (wbc == 4) THEN ! Do them for open BC only
      i = 1
      DO k=2,nz-2
        DO j=2,ny-1
          vmix(i,j,k)=vmix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k))                 &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)))
        END DO
      END DO
    END IF

  END IF
!
  IF( cfcmv4 /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the vertical
!  computational mixing on v'
!
!-----------------------------------------------------------------------
!
    temz = cfcmv4/(4.0*dz**4)

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(tem5(i,j,k)-tem5(i,j,k-1))*                      &
              (rhostr(i,j,k  )+rhostr(i,j-1,k  )                        &
              +rhostr(i,j,k-1)+rhostr(i,j-1,k-1))*temz
        END DO
      END DO
    END DO
!
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=tem1(i,j,k+1)-tem1(i,j,k)
        END DO
      END DO
    END DO

    CALL bvdifzz(tem2, nx,ny,nz,1,nx-1,2,ny-1,tbc, bbc)

    IF( cmix_opt == 0) THEN  !  no monotonic application

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=1,nx-1
            vmix(i,j,k)=vmix(i,j,k)                                     &
                -((tem2(i,j,k+1)-tem2(i,j,k))-(tem2(i,j,k)-tem2(i,j,k-1)))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN  !  4th order mono...

      DO k=2,nz-1
        DO j=2,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=tem2(i,j,k)-tem2(i,j,k-1)
            IF(-tem1(i,j,k)*(tem5(i,j,k)-tem5(i,j,k-1)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=1,nx-1
            vmix(i,j,k)=vmix(i,j,k)-(tem1(i,j,k+1)-tem1(i,j,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=1,nx-1
            tem3(i,j,k)=(tem2(i,j,k+1)-tem2(i,j,k))                     &
                       -(tem2(i,j,k)-tem2(i,j,k-1))
          END DO
        END DO
      END DO
      CALL bvdifzz(tem3, nx,ny,nz,1,nx-1,2,ny-1,tbc, bbc)

      DO k=2,nz-1
        DO j=2,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=tem3(i,j,k)-tem3(i,j,k-1)
            IF(cmix_opt == 3.AND.                                       &
                  tem1(i,j,k)*(tem5(i,j,k)-tem5(i,j,k-1)) < 0.0)THEN
              tem1(i,j,k)=0.0
!           if(j.eq.2) kount = kount + 1
            END IF
          END DO
        END DO
      END DO

!      print*,'On-off switch active for v at ',
!    :            float(kount)/(2*(nx-1)*(nz-2))
!    :           ,' percent of times at t=', curtim


      DO k=2,nz-2
        DO j=2,ny-1
          DO i=1,nx-1
            vmix(i,j,k)=vmix(i,j,k)+(tem1(i,j,k+1)-tem1(i,j,k))
          END DO
        END DO
      END DO

    END IF

  END IF

  IF( cfcmh4 /= 0.0 ) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the horizontal
!  computational mixing on w
!
!-----------------------------------------------------------------------
!
    temx = cfcmh4/(4.0*dx**4)
    temy = cfcmh4/(4.0*dy**4)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(w(i,j,k)-w(i-1,j,k))                             &
                     *(rhostr(i,j  ,k)+rhostr(i-1,j  ,k)+               &
                       rhostr(i,j,k-1)+rhostr(i-1,j,k-1))*temx
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO j=1,ny-1
        tem1(1 ,j,k)=0.0
        tem1(nx,j,k)=0.0
      END DO
    END DO

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(w(i,j,k)-w(i,j-1,k))                             &
                     *(rhostr(i,j,  k)+rhostr(i,j-1,k)+                 &
                       rhostr(i,j,k-1)+rhostr(i,j-1,k-1))*temy
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO i=1,nx-1
        tem2(i,1 ,k)=0.0
        tem2(i,ny,k)=0.0
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Note that tem1 at i=1 and nx are zero, and tem2 at j=1,ny are zero,
!  equivalent to assuming zero normal gradient of s at the boundaries.
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=tem1(i+1,j,k)-tem1(i,j,k)
          tem4(i,j,k)=tem2(i,j+1,k)-tem2(i,j,k)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem3,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(tem4,nx,ny,nz,nbc,sbc,0,tem1)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bwdifxx(tem3, nx,ny,nz,1,ny-1,2,nz-1,ebc, wbc)
    CALL bwdifyy(tem4, nx,ny,nz,1,nx-1,2,nz-1,nbc, sbc)
    CALL acct_stop_inter

    IF( cmix_opt == 0) THEN  !  no monotonic application

      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            wmix(i,j,k)=wmix(i,j,k)-(                                   &
                (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k)) &
                +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN  !  4th order monotonic comp. mixing

      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-1
            tem1(i,j,k)= tem3(i,j,k)-tem3(i-1,j,k)
            IF(-tem1(i,j,k)*(w(i,j,k)-w(i-1,j,k)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-1
        DO j=2,ny-1
          DO i=2,nx-2
            tem2(i,j,k)= tem4(i,j,k)-tem4(i,j-1,k)
            IF(-tem2(i,j,k)*(w(i,j,k)-w(i,j-1,k)) < 0.0) tem2(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            wmix(i,j,k)=wmix(i,j,k)-                                    &
                (tem1(i+1,j,k)-tem1(i,j,k)+tem2(i,j+1,k)-tem2(i,j,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            tem1(i,j,k)=(tem3(i+1,j,k)-tem3(i,j,k))                     &
                       -(tem3(i,j,k)-tem3(i-1,j,k))
            tem2(i,j,k)=(tem4(i,j+1,k)-tem4(i,j,k))                     &
                       -(tem4(i,j,k)-tem4(i,j-1,k))
          END DO
        END DO
      END DO
      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,0,tem3)
        CALL mpsendrecv2dns(tem2,nx,ny,nz,nbc,sbc,0,tem3)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bwdifxx(tem1, nx,ny,nz,1,ny-1,2,nz-1,ebc, wbc)
      CALL bwdifyy(tem2, nx,ny,nz,1,nx-1,2,nz-1,nbc, sbc)
      CALL acct_stop_inter

      kount =0
      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-1
            tem3(i,j,k)= tem1(i,j,k)-tem1(i-1,j,k)
            IF(cmix_opt == 3.AND. tem3(i,j,k)*(w(i,j,k)-w(i-1,j,k)) < 0.0) THEN
              tem3(i,j,k)=0.0
              IF(j == 2) kount = kount+1
            END IF
          END DO
        END DO
      END DO

      DO k=2,nz-1
        DO j=2,ny-1
          DO i=2,nx-2
            tem4(i,j,k)= tem2(i,j,k)-tem2(i,j-1,k)
            IF(cmix_opt == 3.AND. tem4(i,j,k)*(w(i,j,k)-w(i,j-1,k)) < 0.0) &
                tem4(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            wmix(i,j,k)=wmix(i,j,k)+                                    &
                (tem3(i+1,j,k)-tem3(i,j,k)+tem4(i,j+1,k)-tem4(i,j,k))
          END DO
        END DO
      END DO

    END IF

    !wdt update
    IF (ebc == 4) THEN
      i = nx-1
      DO k=2,nz-1
        DO j=2,ny-2
          wmix(i,j,k)=wmix(i,j,k)-(                                     &
              (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
        END DO
      END DO
    END IF
    IF (wbc == 4) THEN
      i = 1
      DO k=2,nz-1
        DO j=2,ny-2
          wmix(i,j,k)=wmix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-tem3(i,j,k)                   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
        END DO
      END DO
    END IF
    IF (nbc == 4) THEN
      j = ny - 1
      DO k=2,nz-1
        DO i=2,nx-2
          wmix(i,j,k)=wmix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
        END DO
      END DO
    END IF
    IF (sbc == 4) THEN
      j = 1
      DO k=2,nz-1
        DO i=2,nx-2
          wmix(i,j,k)=wmix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)              ) )
        END DO
      END DO
    END IF
    IF (wbc == 4.AND.sbc == 4) THEN
      i = 1
      j = 1
      DO k=2,nz-1
        wmix(i,j,k)=wmix(i,j,k)-(                                       &
            (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k))                   &
            +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)) )
      END DO
    END IF
    IF (ebc == 4.AND.sbc == 4) THEN
      i = nx-1
      j = 1
      DO k=2,nz-1
        wmix(i,j,k)=wmix(i,j,k)-(                                       &
            (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))     &
            +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)) )
      END DO
    END IF

    IF (wbc == 4.AND.nbc == 4) THEN
      i = 1
      j = ny-1
      DO k=2,nz-1
        wmix(i,j,k)=wmix(i,j,k)-(                                       &
            (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k))                   &
            +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
      END DO
    END IF
    IF (ebc == 4.AND.nbc == 4) THEN
      i = nx-1
      j = ny-1
      DO k=2,nz-1
        wmix(i,j,k)=wmix(i,j,k)-(                                       &
            (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))     &
            +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
      END DO
    END IF

  END IF

  IF( cfcmv4 /= 0.0 ) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the vertical
!  computational mixing on w
!
!-----------------------------------------------------------------------
!
    temz = cfcmv4/(dz**4)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(w(i,j,k+1)-w(i,j,k))*rhostr(i,j,k)*temz
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=tem1(i,j,k)-tem1(i,j,k-1)
        END DO
      END DO
    END DO

    CALL bwdifzz(tem2, nx,ny,nz,1,nx-1,1,ny-1,tbc, bbc)

    IF(cmix_opt == 0) THEN  ! no monotonic application

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            wmix(i,j,k)=wmix(i,j,k)-(                                   &
                (tem2(i,j,k+1)-tem2(i,j,k))-(tem2(i,j,k)-tem2(i,j,k-1)))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN  !  4th order monotonic 

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=(tem2(i,j,k+1)-tem2(i,j,k))
            IF(-tem1(i,j,k)*(w(i,j,k+1)-w(i,j,k)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            wmix(i,j,k)=wmix(i,j,k)-(tem1(i,j,k)-tem1(i,j,k-1))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem3(i,j,k)=(tem2(i,j,k+1)-tem2(i,j,k))                     &
                       -(tem2(i,j,k)-tem2(i,j,k-1))
          END DO
        END DO
      END DO
      CALL bwdifzz(tem3, nx,ny,nz,1,nx-1,1,ny-1,tbc, bbc)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=(tem3(i,j,k+1)-tem3(i,j,k))
            IF(cmix_opt == 3.AND. tem1(i,j,k)*(w(i,j,k+1)-w(i,j,k)) < 0.0)THEN
              tem1(i,j,k)=0.0
              IF(j == 2) kount=kount+1
            END IF
          END DO
        END DO
      END DO

!      PRINT*,'On-off switch active for w at ',                          &
!              FLOAT(kount)/(2*(nx-1)*(nz-2))                            &
!             ,' percent of times at t=', curtim


      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            wmix(i,j,k)=wmix(i,j,k)+(tem1(i,j,k)-tem1(i,j,k-1))
          END DO
        END DO
      END DO

    END IF

  END IF

  RETURN
END SUBROUTINE cmix4uvw

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CMIX4S                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cmix4s(nx,ny,nz, s ,rhostr, smix, tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generic routine to calculate the fourth order computational mixing
!  for scalar s. The computational mixing is accumulated into array
!  smix on exit.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/3/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  10/17/94 (M. Xue)
!  Subroutines CMIX4PT and CMIX4Q combined into a single generic
!  routine CMIX2S.
!
!  5/25/1998 (M. Xue)
!  Reformulated the computational mixing terms to move rhostr
!  outside the inner-most derivative.
!
!  11/9/2001 (M. Xue, D. Weber, and X. Jin)
!  Added monotonic computational mixing and 6-th order mixing option
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    s        scalar field to which the compuational mixing is to
!             be applied.
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    smix     Array containing turbulent mixing on s.
!
!  OUTPUT:
!
!    smix     Turbulent mixing plus computational mixing on s.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s     (nx,ny,nz)     ! Scalar to be mixed.

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: smix  (nx,ny,nz)     ! Total mixing on s

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j,k
  REAL :: dxinv4, dyinv4
  REAL :: cfcmh4s,cfcmv4s, temx,temy,temz
  INTEGER :: kount
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
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
!  If the coefficients of computational mixing in both the horizontal
!  and vertical are zero, exit this subroutine.
!
!-----------------------------------------------------------------------
!
  cfcmh4s = cfcmh4 * scmixfctr
  cfcmv4s = cfcmv4 * scmixfctr

  IF( cfcmh4s == 0.0 .AND. cfcmv4s == 0.0 ) RETURN

  dxinv4 = dxinv**4
  dyinv4 = dyinv**4

  IF( cfcmh4s /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the 4th order horizontal
!  computational mixing on s.
!
!-----------------------------------------------------------------------
!
    temx = cfcmh4s/(2.0*dx**4)
    temy = cfcmh4s/(2.0*dy**4)


    DO k=1,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem1(i,j,k)=(s(i,j,k)-s(i-1,j,k))                             &
                     *(rhostr(i,j,k)+rhostr(i-1,j,k))*temx
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        tem1(1 ,j,k)=0.0
        tem1(nx,j,k)=0.0
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(s(i,j,k)-s(i,j-1,k))                             &
                     *(rhostr(i,j,k)+rhostr(i,j-1,k))*temy
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO i=1,nx-1
        tem2(i,1 ,k)=0.0
        tem2(i,ny,k)=0.0
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Note that tem1 at i=1 and nx are zero, and tem2 at j=1,ny are zero,
!  equivalent to assuming zero normal gradient of s at the boundaries.
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=tem1(i+1,j,k)-tem1(i,j,k)
          tem4(i,j,k)=tem2(i,j+1,k)-tem2(i,j,k)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem3,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(tem4,nx,ny,nz,nbc,sbc,0,tem1)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bsdifxx(tem3, nx,ny,nz,1,ny-1,2,nz-2,ebc, wbc)
    CALL bsdifyy(tem4, nx,ny,nz,1,nx-1,2,nz-2,nbc, sbc)
    CALL acct_stop_inter

    IF( cmix_opt == 0) THEN  !  no monotonic application

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
            smix(i,j,k)=smix(i,j,k)-(                                   &
                (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k)) &
                +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN   ! 4th order monotonic.

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-1
            tem1(i,j,k)= tem3(i,j,k)-tem3(i-1,j,k)
            IF(-tem1(i,j,k)*(s(i,j,k)-s(i-1,j,k)) < 0.0) tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-2
            tem2(i,j,k)= tem4(i,j,k)-tem4(i,j-1,k)
            IF(-tem2(i,j,k)*(s(i,j,k)-s(i,j-1,k)) < 0.0) tem2(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
            smix(i,j,k)=smix(i,j,k)-                                    &
                (tem1(i+1,j,k)-tem1(i,j,k)+tem2(i,j+1,k)-tem2(i,j,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
            tem1(i,j,k)=(tem3(i+1,j,k)-tem3(i,j,k))                     &
                       -(tem3(i,j,k)-tem3(i-1,j,k))
            tem2(i,j,k)=(tem4(i,j+1,k)-tem4(i,j,k))                     &
                       -(tem4(i,j,k)-tem4(i,j-1,k))
          END DO
        END DO
      END DO

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,0,tem3)
        CALL mpsendrecv2dns(tem2,nx,ny,nz,nbc,sbc,0,tem3)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bsdifxx(tem1, nx,ny,nz,1,ny-1,2,nz-2,ebc, wbc)
      CALL bsdifyy(tem2, nx,ny,nz,1,nx-1,2,nz-2,nbc, sbc)
      CALL acct_stop_inter

      kount= 0
      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-1
            tem3(i,j,k)= tem1(i,j,k)-tem1(i-1,j,k)
!           IF(scmix_opt == 1.AND.cmix_opt == 3.AND.                    &
            IF(cmix_opt == 3.AND.                    &
                  tem3(i,j,k)*(s(i,j,k)-s(i-1,j,k)) < 0.0)THEN
              tem3(i,j,k)=0.0
              IF(j == 2) kount=kount+1
            END IF
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-1
          DO i=2,nx-2
            tem4(i,j,k)= tem2(i,j,k)-tem2(i,j-1,k)
!           IF(scmix_opt == 1.AND.cmix_opt == 3.AND.                    &
            IF(cmix_opt == 3.AND.                    &
                tem4(i,j,k)*(s(i,j,k)-s(i,j-1,k)) < 0.0)                &
                tem4(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
            smix(i,j,k)=smix(i,j,k)+                                    &
                (tem3(i+1,j,k)-tem3(i,j,k)+tem4(i,j+1,k)-tem4(i,j,k))
          END DO
        END DO
      END DO

    END IF

!
!-----------------------------------------------------------------------
!
!  Calculate the mixing term on the boundary assuming that mixed
!  field has zero gradient outside the boundary.
!  The boundary values of smix are needed only for open boundary option.
!
!-----------------------------------------------------------------------
!
    !wdt update
    IF (wbc == 4) THEN
      i = 1
      DO k=2,nz-2
        DO j=2,ny-2
          smix(i,j,k)=smix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-tem3(i,j,k)                   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
        END DO
      END DO
    END IF
    IF (ebc == 4) THEN
      i = nx-1
      DO k=2,nz-2
        DO j=2,ny-2
          smix(i,j,k)=smix(i,j,k)-(                                     &
              (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
        END DO
      END DO
    END IF
    IF (sbc == 4) THEN
      j = 1
      DO k=2,nz-2
        DO i=2,nx-2
          smix(i,j,k)=smix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)              ) )
        END DO
      END DO
    END IF
    IF (nbc == 4) THEN
      j = ny - 1
      DO k=2,nz-2
        DO i=2,nx-2
          smix(i,j,k)=smix(i,j,k)-(                                     &
              (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))   &
              +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
        END DO
      END DO
    END IF
    IF (wbc == 4.AND.sbc == 4) THEN
      i = 1
      j = 1
      DO k=2,nz-2
        smix(i,j,k)=smix(i,j,k)-(                                       &
            (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k))                   &
            +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)) )
      END DO
    END IF
    IF (ebc == 4.AND.sbc == 4) THEN
      i = nx-1
      j = 1
      DO k=2,nz-2
        smix(i,j,k)=smix(i,j,k)-(                                       &
            (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))     &
            +(tem4(i,j+1,k)-tem4(i,j,k))-(tem4(i,j,k)) )
      END DO
    END IF
    IF (wbc == 4.AND.nbc == 4) THEN
      i = 1
      j = ny-1
      DO k=2,nz-2
        smix(i,j,k)=smix(i,j,k)-(                                       &
            (tem3(i+1,j,k)-tem3(i,j,k))-(tem3(i,j,k))                   &
            +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
      END DO
    END IF
    IF (ebc == 4.AND.nbc == 4) THEN
      i = nx-1
      j = ny-1
      DO k=2,nz-2
        smix(i,j,k)=smix(i,j,k)-(                                       &
            (             -tem3(i,j,k))-(tem3(i,j,k)-tem3(i-1,j,k))     &
            +(             -tem4(i,j,k))-(tem4(i,j,k)-tem4(i,j-1,k)) )
      END DO
    END IF

  END IF

  IF( cfcmv4s /= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  If the coefficient is not zero, calculate the vertical
!  computational mixing on s.
!
!-----------------------------------------------------------------------

    temz = cfcmv4s/(2.0*dz**4)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(s(i,j,k)-s(i,j,k-1))                             &
                     *(rhostr(i,j,k)+rhostr(i,j,k-1))*temz
        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,1 )=0.0
        tem1(i,j,nz)=0.0
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=tem1(i,j,k+1)-tem1(i,j,k)
        END DO
      END DO
    END DO

    CALL bsdifzz(tem2, nx,ny,nz,1,nx-1,1,ny-1,tbc, bbc)

    IF( cmix_opt == 0) THEN  !  no mono applied

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            smix(i,j,k)=smix(i,j,k)                                     &
                -((tem2(i,j,k+1)-tem2(i,j,k))-(tem2(i,j,k)-tem2(i,j,k-1)))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 1) THEN  !  4th order mono

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=tem2(i,j,k)-tem2(i,j,k-1)
            IF(-tem1(i,j,k)*(s(i,j,k)-s(i,j,k-1)) < 0.)tem1(i,j,k)=0.0
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            smix(i,j,k)=smix(i,j,k)-(tem1(i,j,k+1)-tem1(i,j,k))
          END DO
        END DO
      END DO

    ELSE IF( cmix_opt == 2 .OR.cmix_opt == 3 ) THEN  ! 6th order 
                                             ! = 2 6th order no mono...
                                             ! = 3 6th order mono...

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            tem3(i,j,k)=(tem2(i,j,k+1)-tem2(i,j,k))                     &
                       -(tem2(i,j,k)-tem2(i,j,k-1))
          END DO
        END DO
      END DO
      CALL bsdifzz(tem3, nx,ny,nz,1,nx-1,1,ny-1,tbc, bbc)

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=tem3(i,j,k)-tem3(i,j,k-1)
!           IF(scmix_opt == 1.AND.cmix_opt == 3.AND.                    &
            IF(cmix_opt == 3.AND.                    &
                  tem1(i,j,k)*(s(i,j,k)-s(i,j,k-1)) < 0.)THEN
              tem1(i,j,k)=0.0
              IF(j == 2) kount=kount+1
            END IF
          END DO
        END DO
      END DO

!      PRINT*,'On-off switch active for pt at ',                         &
!              FLOAT(kount)/(2*(nx-1)*(nz-2))                            &
!             ,' percent of times at t=', curtim


      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            smix(i,j,k)=smix(i,j,k)+(tem1(i,j,k+1)-tem1(i,j,k))
          END DO
        END DO
      END DO

    END IF

  END IF

  RETURN
END SUBROUTINE cmix4s
