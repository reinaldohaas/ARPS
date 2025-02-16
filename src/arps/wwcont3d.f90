!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WCONTRA                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wcontra(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,                 &
           rhostr,rhostru,rhostrv,rhostrw,wcont,ustr,vstr)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate wcont, the contravariant vertical velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  1/4/1993.
!
!  Modification history:
!  8/29/94 (A. Shapiro)
!  Bug fix. Call to vbcwcont moved outside IF block.
!
!  9/9/94 (M. Xue)
!  Optimized.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/28/98 (D. Weber)
!  Added (mapfct(i,j,7-8) and aj3z to improve efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of Cartesian velocity
!             at all time levels (m/s)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!    rhostr   j3 times base state density rhobar(kg/m**3).
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!
!  OUTPUT:
!
!    wcont    Vertical component of contravariant velocity in
!             computational coordinates (m/s)
!
!  WORK ARRAYS:
!
!    ustr     Work array
!    vstr     Work array
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

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: rhostr(nx,ny,nz)     ! j3 times base state density rhobar
                               ! (kg/m**3).
  REAL :: rhostru(nx,ny,nz)    ! Average rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Average rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Average rhostr at w points (kg/m**3).

  REAL :: wcont (nx,ny,nz)     ! Vertical velocity in computational
                               ! coordinates (m/s)

  REAL :: ustr  (nx,ny,nz)     ! temporary work array
  REAL :: vstr  (nx,ny,nz)     ! temporary work array

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( crdtrns == 0 ) THEN  ! No coord. transformation case.

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          wcont(i,j,k)=w(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF( ternopt == 0) THEN

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          wcont(i,j,k)=w(i,j,k)/aj3z(i,j,k)
        END DO
      END DO
    END DO

  ELSE

    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx
          ustr(i,j,k)=u(i,j,k)*rhostru(i,j,k)
        END DO
      END DO
    END DO

    DO k= 1,nz-1
      DO j= 1,ny
        DO i= 1,nx-1
          vstr(i,j,k)=v(i,j,k)*rhostrv(i,j,k)
        END DO
      END DO
    END DO

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          wcont(i,j,k)= (                                               &
              ((ustr(i  ,j,k)+ustr(i  ,j,k-1))*j1(i  ,j,k)              &
              +(ustr(i+1,j,k)+ustr(i+1,j,k-1))*j1(i+1,j,k)              &
              +(vstr(i  ,j,k)+vstr(i  ,j,k-1))*j2(i  ,j,k)              &
              +(vstr(i,j+1,k)+vstr(i,j+1,k-1))*j2(i,j+1,k))             &
              * mapfct(i,j,8)                                           &
              / rhostrw(i,j,k) + w(i,j,k)                               &
              ) /aj3z(i,j,k)
        END DO
      END DO
    END DO

  END IF

  CALL vbcwcont(nx,ny,nz,wcont)

  RETURN
END SUBROUTINE wcontra
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WCTOW                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wctow(nx,ny,nz,u,v,wcont,mapfct,                             &
           j1,j2,j3,aj3z,rhostr,rhostru,rhostrv,rhostrw,ubc,bbc,        &
           w,                                                           &
           ustr,vstr)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate w from wcont, the contravariant vertical velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/28/1998.
!
!  Modification history:
!
!  9/28/98 (D. Weber)
!  Added aj3z to improve efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    wcont    Vertical component of contravariant velocity in
!             computational coordinates (m/s)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!    rhostr   j3 times base state density rhobar(kg/m**3).
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!    ubc, bbc Flags for upper and bottom boundary conditions
!
!  OUTPUT:
!
!    w        Vertical component of Cartesian velocity (m/s)
!
!  WORK ARRAYS:
!
!    ustr     Work array
!    vstr     Work array
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

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Vertical velocity in computational
                               ! coordinates (m/s)

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: rhostr(nx,ny,nz)     ! j3 times base state density rhobar
                               ! (kg/m**3).
  REAL :: rhostru(nx,ny,nz)    ! Average rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Average rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Average rhostr at w points (kg/m**3).
  INTEGER :: ubc, bbc          ! Flags for upper and bottom boundary conditions

  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: ustr  (nx,ny,nz)     ! temporary work array
  REAL :: vstr  (nx,ny,nz)     ! temporary work array

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( crdtrns == 0 ) THEN  ! No coord. transformation case.

    DO k= 1,nz
      DO j= 1,ny-1
        DO i= 1,nx-1
          w(i,j,k)=wcont(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF( ternopt == 0) THEN

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          w(i,j,k)=wcont(i,j,k)*aj3z(i,j,k)
        END DO
      END DO
    END DO

    CALL vbcw(nx,ny,nz,w,wcont,ubc,bbc,u,v,                             &
              rhostr,rhostru,rhostrv,rhostrw,                           &
              j1,j2,j3)

  ELSE

    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx
          ustr(i,j,k)=u(i,j,k)*rhostru(i,j,k)
        END DO
      END DO
    END DO

    DO k= 1,nz-1
      DO j= 1,ny
        DO i= 1,nx-1
          vstr(i,j,k)=v(i,j,k)*rhostrv(i,j,k)
        END DO
      END DO
    END DO

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          w(i,j,k)= wcont(i,j,k)*aj3z(i,j,k) -                          &
              ((ustr(i  ,j,k)+ustr(i  ,j,k-1))*j1(i  ,j,k)              &
              +(ustr(i+1,j,k)+ustr(i+1,j,k-1))*j1(i+1,j,k)              &
              +(vstr(i  ,j,k)+vstr(i  ,j,k-1))*j2(i  ,j,k)              &
              +(vstr(i,j+1,k)+vstr(i,j+1,k-1))*j2(i,j+1,k))             &
              * mapfct(i,j,8)/ rhostrw(i,j,k)

          ! Only wcont should multiply aj3z. Detected by Tina K. Chow

        END DO
      END DO
    END DO

    CALL vbcw(nx,ny,nz,w,wcont,ubc,bbc,u,v,                             &
              rhostr,rhostru,rhostrv,rhostrw,                           &
              j1,j2,j3)

  END IF

  RETURN
END SUBROUTINE wctow
