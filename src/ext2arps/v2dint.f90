!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE V2DINT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE v2dinta(nx,ny,nz, ibeg,iend,jbeg,jend,kbeg,kend,             &
           s,z,zint,ss1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate 3-D data to a given horizontal level.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  12/18/92.
!
!  MODIFICATION HISTORY:
!
!  2/02/95 K. Brewster
!  Corrected bug and changed so that extapolation is allowed.
!  3/27/95 K. Brewster
!  Added the valid ranges to the parameter list.
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
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibeg,iend   Valid range of i index
!    jbeg,jend   Valid range of j index
!    kbeg,kend   Valid range of k index
!
!    s        3-dimensional array of data to contour
!    s1       2-dimensional array of data to contour
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in physical space (m)
!
!    zint     horizontal level to find
!
!  OUTPUT:
!    ss1      interpolated 3-D data to a given horizontal level
!
!-----------------------------------------------------------------------
!
!  Parameters of output
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend

  REAL :: s(nx,ny,nz)          ! 3-dimensional array of data to contour
  REAL :: z(nx,ny,nz)          ! z coordinate of grid points
                               ! in physical space (m)
  REAL :: zint
  REAL :: ss1(nx,ny)           ! interpolated 3-D data to a
                               ! given horizontal level

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
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
!  Find index for interpolation
!
!-----------------------------------------------------------------------
!
  DO j=jbeg,jend
    DO i=ibeg,iend
      DO k=kbeg+1,kend
        IF(zint <= z(i,j,k)) EXIT
      END DO

      ss1(i,j)=s(i,j,k-1)+(s(i,j,k)-s(i,j,k-1))*                   &
               (zint-z(i,j,k-1))/(z(i,j,k)-z(i,j,k-1))

    ! J.Case, ENSCO Inc. (3/16/2005)
    ! JTM added this check for Linux cluster since very small numbers were
    ! causing WD_WPGD to hang up (6-14-02).

    !  if (abs(ss1(i,j)) .lt. 1.e-10) ss1(i,j)=0.

    END DO
  END DO

  RETURN
END SUBROUTINE v2dinta
