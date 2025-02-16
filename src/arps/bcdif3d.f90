!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BUDIFXX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE budifxx(d2dxu,nx,ny,nz,jbgn,jend,kbgn,kend,ebc,wbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the east and west boundary values of del**2(u-ubar)/delx**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the u-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dxu        A 3-D array whose east and west boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dxu.
!    ebc, wbc     Parameter defining the type of boundary conditions
!                 for array d2dxu.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    jbgn         Starting j index
!    jend         Ending j index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dxu        Array d2dxu with updated east and west boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: jbgn,jend         ! Domain for j computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dxu(nx,ny,nz)      ! Input array
  INTEGER :: ebc, wbc          ! Control parameter of east and west b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, k
!
  INCLUDE 'mp.inc'
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
!  West boundary condition
!
!-----------------------------------------------------------------------
!
  IF( wbc == 1) THEN         ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxu(1,j,k) = -d2dxu(3,j,k)
      END DO
    END DO

  ELSE IF( wbc == 2) THEN     ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxu(1,j,k) = d2dxu(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF (wbc /= 0) THEN    ! Any of the 4 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxu(1,j,k) = d2dxu(2,j,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  East boundary condition
!
!-----------------------------------------------------------------------
!
  IF( ebc == 1) THEN         ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxu(nx,j,k) = -d2dxu(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN      ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxu(nx,j,k) = d2dxu(3,j,k)
      END DO
      END DO
    END IF

  ELSE IF (ebc /= 0) THEN    ! Any of the 4 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxu(nx,j,k) = d2dxu(nx-1,j,k)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE budifxx

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BUDIFYY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE budifyy(d2dyu,nx,ny,nz,ibgn,iend,kbgn,kend,nbc,sbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the north and south boundary values of del**2(u-ubar)/dely**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the u-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dyu        A 3-D array whose north and south boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2ydu.
!    nbc, sbc     Parameter defining the type of boundary conditions
!                 for array d2ydu.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dyu        Array d2ydu with updated north and south boundary
!                 values.
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
  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dyu(nx,ny,nz)      ! Input array
  INTEGER :: nbc, sbc          ! Control parameter for north/south b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  South boundary condition
!
!-----------------------------------------------------------------------
!

  IF( sbc == 2) THEN         ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyu(i,1,k) = d2dyu(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF (sbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyu(i,1,k) = d2dyu(i,2,k)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  North boundary condition
!
!-----------------------------------------------------------------------
!
  IF(nbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyu(i,ny-1,k) = d2dyu(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF (nbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyu(i,ny-1,k) = d2dyu(i,ny-2,k)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE budifyy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BUDIFZZ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE budifzz(d2dzu,nx,ny,nz,ibgn,iend,jbgn,jend,tbc,bbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary values of del**2(u-ubar)/delz**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the u-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dzu        A 3-D array whose top and bottom boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dzu
!    tbc, bbc     Parameter defining the type of boundary conditions
!                 for array d2dzu.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    jbgn         Starting j index
!    jend         Ending j index
!
!  OUTPUT:
!
!    d2dzu        Array d2dzu with updated top and bottom boundary
!                 values.
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
  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: jbgn,jend         ! Domain for j computations
  REAL :: d2dzu(nx,ny,nz)      ! Input aray
  INTEGER :: tbc, bbc          ! Control parameter for top/bottom b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Top boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(tbc == 2) THEN          ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzu(i,j,nz-1) = d2dzu(i,j,2)
      END DO
    END DO

  ELSE                       ! Rigid lid and zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzu(i,j,nz-1) = d2dzu(i,j,nz-2)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Bottom boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(bbc == 2) THEN          ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzu(i,j,1) = d2dzu(i,j,nz-2)
      END DO
    END DO

  ELSE                       ! Solid ground and zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzu(i,j,1) = d2dzu(i,j,2)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE budifzz

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BVDIFXX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bvdifxx(d2dxv,nx,ny,nz,jbgn,jend,kbgn,kend,ebc,wbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the east and west boundary values of del**2(v-vbar)/delx**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the v-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dxv        A 3-D array whose east and west boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dxv
!    ebc, wbc     Parameter defining the type of boundary conditions
!                 for array d2dxv.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    jbgn         Starting j index
!    jend         Ending j index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dxv        Array d2dxv with updated east and west boundary
!                 values.
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
!
  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: jbgn,jend         ! Domain for j computations
  INTEGER :: kbgn,kend         ! Domain for k computations
  REAL :: d2dxv(nx,ny,nz)      ! Input aray
  INTEGER :: ebc, wbc          ! Control parameter for east/west b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  West boundary condition
!
!-----------------------------------------------------------------------
!
  IF(wbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxv(1,j,k) = d2dxv(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF (wbc /= 0) THEN    ! All the other 5 types of conditions

    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxv(1,j,k) = d2dxv(2,j,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  East boundary condition
!
!-----------------------------------------------------------------------
!
  IF(ebc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxv(nx-1,j,k) = d2dxv(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF (ebc /= 0) THEN    ! All the other 5 types of conditions

    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxv(nx-1,j,k) = d2dxv(nx-2,j,k)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bvdifxx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BVDIFYY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bvdifyy(d2dyv,nx,ny,nz,ibgn,iend,kbgn,kend,nbc,sbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the north and south boundary values of del**2(v-vbar)/dely**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the v-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dyv        A 3-D array whose north and south boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dyv
!    nbc, sbc     Parameter defining the type of boundary conditions
!                 for array d2dyv.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dyv        Array d2dyv with updated north and south boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dyv(nx,ny,nz)      ! Input array
  INTEGER :: nbc, sbc          ! Control parameter for north/south b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  South boundary condition
!
!-----------------------------------------------------------------------
!
  IF(sbc == 1) THEN         ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyv(i,1,k) = -d2dyv(i,3,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN      ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyv(i,1,k) = d2dyv(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF (sbc /= 0) THEN    ! Any of the 4 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyv(i,1,k) = d2dyv(i,2,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  North boundary condition
!
!-----------------------------------------------------------------------
!
  IF(nbc == 1) THEN         ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyv(i,ny,k) = -d2dyv(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN      ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyv(i,ny,k) = d2dyv(i,3,k)
      END DO
      END DO
    END IF

  ELSE IF (nbc /= 0) THEN    ! Any of the 4 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyv(i,ny,k) = d2dyv(i,ny-1,k)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bvdifyy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BVDIFZZ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bvdifzz(d2dzv,nx,ny,nz,ibgn,iend,jbgn,jend,tbc,bbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary values of del**2(v-vbar)/delz**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the v-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dzv        A 3-D array whose top and bottom boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dvz.
!    tbc, bbc     Parameter defining the type of boundary conditions
!                 for array d2dvz.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    jbgn         Starting j index
!    jend         Ending j index
!
!  OUTPUT:
!
!    d2dzv        Array d2dzv with updated top and bottom boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: jbgn,jend         ! Domain for j computations
  REAL :: d2dzv(nx,ny,nz)      ! Input array
  INTEGER :: tbc, bbc          ! Control parameter for top/bottom b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Top boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(tbc == 2) THEN          ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzv(i,j,nz-1) = d2dzv(i,j,2)
      END DO
    END DO

  ELSE                       ! Rigid lid and zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzv(i,j,nz-1) = d2dzv(i,j,nz-2)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Bottom boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(bbc == 2) THEN          ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzv(i,j,1) = d2dzv(i,j,nz-2)
      END DO
    END DO

  ELSE                       ! Solid ground and zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzv(i,j,1) = d2dzv(i,j,2)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bvdifzz

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BWDIFXX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bwdifxx(d2dxw,nx,ny,nz,jbgn,jend,kbgn,kend,ebc,wbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the east and west boundary values of del**2(w-wbar)/delx**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the w-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dxw        A 3-D array whose east and west boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dxw.
!    ebc, wbc     Parameter defining the type of boundary conditions
!                 for array d2dxw.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    jbgn         Starting j index
!    jend         Ending j index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dxw        Array d2dxw with updated east and west boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: jbgn,jend         ! Domain for j computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dxw(nx,ny,nz)      ! Input array
  INTEGER :: ebc, wbc          ! Control parameter for east/west b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  West boundary condition
!
!-----------------------------------------------------------------------
!
  IF(wbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxw(1,j,k) = d2dxw(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc /= 0) THEN     ! Any of the 5 remaining boundary
                             ! conditions.
    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxw(1,j,k) = d2dxw(2,j,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  East boundary condition
!
!-----------------------------------------------------------------------
!
  IF(ebc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxw(nx-1,j,k) = d2dxw(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc /= 0) THEN     ! Any of the 5 remaining boundary
                             ! conditions.
    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxw(nx-1,j,k) = d2dxw(nx-2,j,k)
      END DO
    END DO

  END IF
  RETURN
END SUBROUTINE bwdifxx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BWDIFYY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bwdifyy(d2dyw,nx,ny,nz,ibgn,iend,kbgn,kend,nbc,sbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the north and south boundary values of del**2(w-wbar)/dely**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the w-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dyw        A 3-D array whose north and south boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dyw.
!    nbc, sbc     Parameter defining the type of boundary conditions
!                 for array d2dyw.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dyw        Array d2dyw with updated north and south boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dyw(nx,ny,nz)      ! Input array
  INTEGER :: nbc, sbc          ! Control parameter for north/south b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  South boundary condition
!
!-----------------------------------------------------------------------
!
  IF(sbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyw(i,1,k) = d2dyw(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF (sbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions.
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyw(i,1,k) = d2dyw(i,2,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  North boundary condition
!
!-----------------------------------------------------------------------
!
  IF(nbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyw(i,ny-1,k) = d2dyw(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF (nbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions.
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dyw(i,ny-1,k) = d2dyw(i,ny-2,k)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bwdifyy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BWDIFZZ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bwdifzz(d2dzw,nx,ny,nz,ibgn,iend,jbgn,jend,tbc,bbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary values of del**2(w-wbar)/delz**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the w-equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  8/31/93 (M. Xue)
!  Corrections made to the IF test statements.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dzw        A 3-D array whose top and bottom boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dzw.
!    tbc, bbc     Parameter defining the type of boundary conditions
!                 for array d2dzw.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    jbgn         Starting j index
!    jend         Ending j index
!
!  OUTPUT:
!
!    d2dzw        Array d2dzw with updated top and bottom boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: jbgn,jend         ! Domain for j computations

  REAL :: d2dzw(nx,ny,nz)      ! Input array
  INTEGER :: tbc, bbc          ! Control parameter for top/bottom b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Top boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(tbc == 1) THEN        ! Rigid lid at the top (mirror condition)

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzw(i,j,nz) = -d2dzw(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN      ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzw(i,j,nz) = d2dzw(i,j,3)
      END DO
    END DO

  ELSE                       ! Zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzw(i,j,nz) = d2dzw(i,j,nz-2)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Bottom boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(bbc == 1) THEN          ! Solid ground (mirror condition)

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzw(i,j,1) = -d2dzw(i,j,3)
      END DO
    END DO

  ELSE IF(bbc == 2) THEN      ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzw(i,j,1) = d2dzw(i,j,nz-2)
      END DO
    END DO

  ELSE                       ! Zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzw(i,j,1) = d2dzw(i,j,3)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bwdifzz

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BSDIFXX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bsdifxx(d2dxs,nx,ny,nz,jbgn,jend,kbgn,kend,ebc,wbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the east and west boundary values of del**2(scalar)/delx**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the scalar equations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dzs        A 3-D array whose east and west boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dzs.
!    ebc, wbc     Parameter defining the type of boundary conditions
!                 for array d2dzs.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    jbgn         Starting j index
!    jend         Ending j index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dzs        Array d2dzs with updated east and west boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: jbgn,jend         ! Domain for j computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dxs(nx,ny,nz)      ! Input array
  INTEGER :: ebc, wbc          ! Control parameter for east/west b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  West boundary condition
!
!-----------------------------------------------------------------------
!
  IF(wbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxs(1,j,k) = d2dxs(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF (wbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxs(1,j,k) = d2dxs(2,j,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  East boundary condition
!
!-----------------------------------------------------------------------
!
  IF(ebc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxs(nx-1,j,k) = d2dxs(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF (ebc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO j=jbgn,jend
        d2dxs(nx-1,j,k) = d2dxs(nx-2,j,k)
      END DO
    END DO

  END IF
  RETURN
END SUBROUTINE bsdifxx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BSDIFYY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bsdifyy(d2dys,nx,ny,nz,ibgn,iend,kbgn,kend,nbc,sbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the north and south boundary values of del**2(scalar)/dely**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the scalar equations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dys        A 3-D array whose north and south boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dys
!    nbc, sbc     Parameter defining the type of boundary conditions
!                 for array d2dys.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    kbgn         Starting k index
!    kend         Ending k index
!
!  OUTPUT:
!
!    d2dys        Array d2dys with updated north and south boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: d2dys(nx,ny,nz)      ! Input array
  INTEGER :: nbc, sbc          ! Control parameter for north/south b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, k
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  South boundary condition
!
!-----------------------------------------------------------------------
!
  IF(sbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dys(i,1,k) = d2dys(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF (sbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dys(i,1,k) = d2dys(i,2,k)
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  North boundary condition
!
!-----------------------------------------------------------------------
!
  IF(nbc == 2) THEN          ! Periodic boundary condition

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        d2dys(i,ny-1,k) = d2dys(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF (nbc /= 0) THEN    ! Any of the 5 remaining boundary
                             ! conditions
    DO k=kbgn,kend
      DO i=ibgn,iend
        d2dys(i,ny-1,k) = d2dys(i,ny-2,k)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bsdifyy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BSDIFZZ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bsdifzz(d2dzs,nx,ny,nz,ibgn,iend,jbgn,jend,tbc,bbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary values of del**2(scalar)/delz**2
!  which appear in the calculation of the 4th order computational
!  mixing term of the scalar equations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/10/1992.
!
!  MODIFICATION HISTORY:
!
!  6/10/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    d2dzs        A 3-D array whose top and bottom boundary values
!                 are set in this routine.
!    nx, ny, nz   Dimensions of array d2dzs
!    tbc, bbc     Parameter defining the type of boundary conditions
!                 for array d2dzs.
!                  = 1 for rigid wall boundary;
!                  = 2 for periodic boundary;
!                  = 3 for zero gradient;
!                  = 4 for radiation (open) boundary;
!                  = 5 for user specified;
!                  = 6 for nested grid.
!
!    ibgn         Starting i index
!    iend         Ending i index
!    jbgn         Starting j index
!    jend         Ending j index
!
!  OUTPUT:
!
!    d2dzs        Array d2dzu with updated top and bottom boundary
!                 values.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: jbgn,jend         ! Domain for j computations

  REAL :: d2dzs(nx,ny,nz)      ! Input array
  INTEGER :: tbc,bbc           ! Control parameter for top/bottom b.c.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Top boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(tbc == 2) THEN          ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzs(i,j,nz-1) = d2dzs(i,j,2)
      END DO
    END DO

  ELSE                       ! Rigid lid and zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzs(i,j,nz-1) = d2dzs(i,j,nz-2)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Bottom boundary condition:
!
!-----------------------------------------------------------------------
!
  IF(bbc == 2) THEN          ! Periodic top and bottom condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzs(i,j,1) = d2dzs(i,j,nz-2)
      END DO
    END DO

  ELSE                       ! Solid ground and zero normal gradient

    DO j=jbgn,jend
      DO i=ibgn,iend
        d2dzs(i,j,1) = d2dzs(i,j,2)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bsdifzz

