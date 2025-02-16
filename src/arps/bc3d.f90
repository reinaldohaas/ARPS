!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCU                        ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcu(nx,ny,nz,dtsml,                                          &
           u, udteb,udtwb,udtnb,udtsb,                                  &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for the u-velocity component. Please
!  note that the values at the corner points may depend on the order
!  that e-w, n-s and t-b boundary conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/04/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/07/92 (M. Xue)
!  Modified to take in u at time tfuture only.
!
!  6/15/92 (M. Xue and H. Jin)
!  Implemented open boundary condition
!
!  10/6/92 (MX)
!  Assignment of the boundary conditions at corner columns moved to
!  the front of the top/bottom condition assignment.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2/19/95 (K. Brewster)
!  Separated the application of external and user-supplied BC
!  processing from radiation BC.  Added loops 1715 and 1725 to
!  correctly apply mixed rigid wall and radiation to SW and SE
!  corner points.
!
!  9/10/1995 (Y. Richardson)
!  Fixed a bug with with the open boundary at the corner points.
!
!  08/04/2003 (Yunheng Wang)
!  Fixed a bug with the periodic boundary at the corner points
!  mpi runs. Please note that the boundary flags (web, ebc, nbc and sbc)
!  have been set to zero in the internal processors. So the
!  corresponded global values should be used.
!
!  This comments also apply to subroutine bcv, lbcw, bcp and bcsclr.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml    The small time step size (s)
!
!    u        Interior domain values of u-velocity at tfuture (m/s)
!
!    udteb    Time tendency of the u field at the east boundary
!    udtwb    Time tendency of the u field at the west boundary
!    udtnb    Time tendency of the u field at the north boundary
!    udtsb    Time tendency of the u field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    u        The u-velocity over the entire domain at tfuture (m/s)
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

  REAL :: dtsml                ! The small time step size (s)

  REAL :: u  (nx,ny,nz)        ! Total u-velocity at tfuture (m/s)

  REAL :: udteb (ny,nz)        ! Time tendency of u at east boundary
  REAL :: udtwb (ny,nz)        ! Time tendency of u at west boundary
  REAL :: udtnb (nx,nz)        ! Time tendency of u at north boundary
  REAL :: udtsb (nx,nz)        ! Time tendency of u at south boundary
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  INCLUDE 'mp.inc'

  INTEGER :: mptag, ierror
  INTEGER :: source, dest

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny-1
        u(1,j,k)=-u(3,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        u(1,j,k)=u(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny-1
        u(1,j,k)=u(3,j,k)
      END DO
    END DO

  ELSE IF(wbc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-2
        u(1,j,k)=u(1,j,k)+udtwb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(wbc == 5 .OR. wbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny-1
        u(1,j,k)=u(3,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCU', wbc
    CALL arpsstop ("arpstop called from bcu west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny-1
        u(nx,j,k)=-u(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
        DO j=1,ny-1
          u(nx,j,k)=u(3,j,k)
        END DO
      END DO
    END IF

  ELSE IF(ebc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny-1
        u(nx,j,k)=u(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-2
        u(nx,j,k)=u(nx,j,k)+udteb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(ebc == 5 .OR. ebc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny-1
        u(nx,j,k)=u(nx-2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCU', ebc
    CALL arpsstop ("arpstop called from bcu east bc",1)

  END IF

  5002  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx
        u(i,ny-1,k)=u(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx
        u(i,ny-1,k)=u(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx
        u(i,ny-1,k)=u(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO i=2,nx-1
        u(i,ny-1,k)=u(i,ny-1,k)+udtnb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(nbc == 5 .OR. nbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx
        u(i,ny-1,k)=u(i,ny-2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCU', nbc
    CALL arpsstop ("arpstop called from bcu north bc",1)

  END IF

  5003  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx
        u(i,1,k)=u(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx
        u(i,1,k)=u(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx
        u(i,1,k)=u(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-2
      DO i=2,nx-1
        u(i,1,k)=u(i,1,k)+udtsb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(sbc == 5 .OR. sbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx
        u(i,1,k)=u(i,2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCU', sbc
    CALL arpsstop ("arpstop called from bcu south bc",1)

  END IF

  5004  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southwest corner based on the
!  boundary condition types on the south and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = (nproc_y -1)*nproc_x
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          u(1,1,k)=u(1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_x - 1
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(nx-2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          u(1,1,k)=u(nx-2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(sbc+wbc /= 0)) THEN

    DO k=2,nz-2
      u(1,1,k)=u(1,1,k)+udtwb(1,k)*dtsml
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      u(1,1,k)=u(1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.wbc == 1) THEN

    DO k=2,nz-2
      u(1,1,k)=-u(3,1,k)
    END DO

  ELSE IF(sbc == 4.AND.(wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      u(1,1,k)=u(3,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southeast corner based on the
!  boundary condition types on the south and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2 .AND. ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_y*nproc_x  -1
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(nx,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(nx,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          u(nx,1,k)=u(nx,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4 .AND. ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = 0
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(3,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(nx,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          u(nx,1,k)=u(3,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4 .OR. sbc == 0).AND.(ebc == 4 .OR. ebc == 0).AND.(sbc+ebc /= 0)) THEN

    DO k=2,nz-2
      u(nx,1,k)=u(nx,1,k)+udteb(1,k)*dtsml
    END DO

  ELSE IF((sbc == 1 .OR. sbc == 3 .OR. sbc == 5) .AND. ebc == 4) THEN

    DO k=2,nz-2
      u(nx,1,k)=u(nx,2,k)
    END DO

  ELSE IF(sbc == 4 .AND. ebc == 1) THEN

    DO k=2,nz-2
      u(nx,1,k)=-u(nx-2,1,k)
    END DO

  ELSE IF(sbc == 4 .AND. (ebc == 3 .OR. ebc == 5)) THEN

    DO k=2,nz-2
      u(nx,1,k)=u(nx-2,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northwest corner based on the
!  boundary condition types on the north and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = 0
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          u(1,ny-1,k)=u(1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x >1) THEN
      source = nproc_y*nproc_x - 1
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(nx-2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          u(1,ny-1,k)=u(nx-2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(nbc+wbc /= 0)) THEN

    DO k=2,nz-2
      u(1,ny-1,k)=u(1,ny-1,k)+udtwb(ny-1,k)*dtsml
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      u(1,ny-1,k)=u(1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.wbc == 1) THEN

    DO k=2,nz-2
      u(1,ny-1,k)=-u(3,ny-1,k)
    END DO

  ELSE IF(nbc == 4.AND.(wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      u(1,ny-1,k)=u(3,ny-1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northeast corner based on the
!  boundary condition types on the north and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_x - 1
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(nx,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(nx,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          u(nx,ny-1,k)=u(nx,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = (nproc_y - 1)* nproc_x
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(u(3,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(u(nx,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          u(nx,ny-1,k)=u(3,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(nbc+ebc /= 0)) THEN

    DO k=2,nz-2
      u(nx,ny-1,k)=u(nx,ny-1,k)+udteb(ny-1,k)*dtsml
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      u(nx,ny-1,k)=u(nx,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.ebc == 1) THEN

    DO k=2,nz-2
      u(nx,ny-1,k)=-u(nx-2,ny-1,k)
    END DO

  ELSE IF(nbc == 4.AND.(ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      u(nx,ny-1,k)=u(nx-2,ny-1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the top boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny-1
      DO i=1,nx
        u(i,j,nz-1)=u(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx
        u(i,j,nz-1)=u(i,j,2)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN  ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx
        u(i,j,nz-1)=u(i,j,nz-2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCU', tbc
    CALL arpsstop ("arpstop called from bcu top bc",1)

  END IF

  5005  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary conditions
!
!-----------------------------------------------------------------------
!

  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny-1
      DO i=1,nx
        u(i,j,1)=u(i,j,2)
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx
        u(i,j,1)=u(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx
        u(i,j,1)=u(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCU', bbc
    CALL arpsstop ("arpstop called from bcu bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE bcu
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCV                        ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcv(nx,ny,nz,dtsml,                                          &
           v, vdteb,vdtwb,vdtnb,vdtsb,                                  &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for the v-velocity component. Please note
!  that the values at the corner points may depend on the order that e-w,
!  n-s and t-b boundary conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/04/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/07/92 (M. Xue)
!  Modified to take in v at time tfuture only.
!
!  6/15/92 (M. Xue and H. Jin)
!  Implemented open boundary condition
!
!  10/6/92 (MX)
!  Assignment of the boundary conditions at corner columns moved to
!  the front of the top/bottom condition assignment.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2/19/95 (K. Brewster)
!  Separated the application of external and user-supplied BC
!  processing from radiation BC.
!
!  9/10/1995 (Y. Richardson)
!  Fixed a bug with with the open boundary at the corner points.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml    The small time step size (s)
!
!    v        Interior domain values of v-velocity at tfuture (m/s)
!
!    vdteb    Time tendency of the v field at the east boundary
!    vdtwb    Time tendency of the v field at the west boundary
!    vdtnb    Time tendency of the v field at the north boundary
!    vdtsb    Time tendency of the v field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    v        The v-velocity over the entire domain at tfuture (m/s)
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

  REAL :: dtsml                ! The small time step size (s)

  REAL :: v     (nx,ny,nz)     ! Total v-velocity at tfuture (m/s)

  REAL :: vdteb (ny,nz)        ! Time tendency of v at east boundary
  REAL :: vdtwb (ny,nz)        ! Time tendency of v at west boundary
  REAL :: vdtnb (nx,nz)        ! Time tendency of v at north boundary
  REAL :: vdtsb (nx,nz)        ! Time tendency of v at south boundary
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
!
  INCLUDE 'mp.inc'

  INTEGER :: mptag, ierror
  INTEGER :: source, dest

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny
        v(1,j,k)=v(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny
        v(1,j,k)=v(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny
        v(1,j,k)=v(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-1
        v(1,j,k)=v(1,j,k)+vdtwb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(wbc == 5 .OR. wbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny
        v(1,j,k)=v(2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCV', wbc
    CALL arpsstop ("arpstop called from bcv west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny
        v(nx-1,j,k)=v(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny
        v(nx-1,j,k)=v(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny
        v(nx-1,j,k)=v(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-1
        v(nx-1,j,k)=v(nx-1,j,k)+vdteb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(ebc == 5 .OR. ebc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny
        v(nx-1,j,k)=v(nx-2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCV', ebc
    CALL arpsstop ("arpstop called from bcv east bc",1)

  END IF

  5002  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx-1
        v(i,ny,k)=-v(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        v(i,ny,k)=v(i,3,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx-1
        v(i,ny,k)=v(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-2
      DO i=2,nx-2
        v(i,ny,k)=v(i,ny,k)+vdtnb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(nbc == 5 .OR. nbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx-1
        v(i,ny,k)=v(i,ny-2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCV', nbc
    CALL arpsstop ("arpstop called from bcv north bc",1)

  END IF

  5003  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx-1
        v(i,1,k)=-v(i,3,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        v(i,1,k)=v(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx-1
        v(i,1,k)=v(i,3,k)
      END DO
    END DO

  ELSE IF(sbc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO i=2,nx-2
        v(i,1,k)=v(i,1,k)+vdtsb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(sbc == 5 .OR. sbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx-1
        v(i,1,k)=v(i,3,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCV', sbc
    CALL arpsstop ("arpstop called from bcv south bc",1)

  END IF

  5004  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southwest corner based on the
!  boundary condition types on the south and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = (nproc_y -1)*nproc_x
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          v(1,1,k)=v(1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_x - 1
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(nx-2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          v(1,1,k)=v(nx-2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(sbc+wbc /= 0)) THEN

    DO k=2,nz-2
      v(1,1,k)=v(1,1,k)+vdtsb(1,k)*dtsml
    END DO

  ELSE IF(sbc == 1.AND.wbc == 4) THEN

    DO k=2,nz-2
      v(1,1,k)=-v(1,3,k)
    END DO

  ELSE IF((sbc == 3.OR.sbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      v(1,1,k)=v(1,3,k)
    END DO

  ELSE IF(sbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      v(1,1,k)=v(2,1,k)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southeast corner based on the
!  boundary condition types on the south and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_y*nproc_x  -1
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(nx-1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          v(nx-1,1,k)=v(nx-1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = 0
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          v(nx-1,1,k)=v(2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(sbc+ebc /= 0)) THEN

    DO k=2,nz-2
      v(nx-1,1,k)=v(nx-1,1,k)+vdtsb(nx-1,k)*dtsml
    END DO

  ELSE IF(sbc == 1.AND.ebc == 4) THEN

    DO k=2,nz-2
      v(nx-1,1,k)=-v(nx-1,3,k)
    END DO

  ELSE IF((sbc == 3.OR.sbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      v(nx-1,1,k)=v(nx-1,3,k)
    END DO

  ELSE IF(sbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      v(nx-1,1,k)=v(nx-2,1,k)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northwest corner based on the
!  boundary condition types on the north and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = 0
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(1,3,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(1,ny,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          v(1,ny,k)=v(1,3,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_y*nproc_x - 1
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(nx-2,ny,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(1,ny,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          v(1,ny,k)=v(nx-2,ny,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(nbc+wbc /= 0)) THEN

    DO k=2,nz-2
      v(1,ny,k)=v(1,ny,k)+vdtnb(1,k)*dtsml
    END DO

  ELSE IF(nbc == 1.AND.wbc == 4) THEN

    DO k=2,nz-2
      v(1,ny,k)=-v(1,ny-2,k)
    END DO

  ELSE IF((nbc == 3.OR.nbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      v(1,ny,k)=v(1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      v(1,ny,k)=v(2,ny,k)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northeast corner based on the
!  boundary condition types on the north and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_x - 1
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(nx-1,3,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(nx-1,ny,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          v(nx-1,ny,k)=v(nx-1,3,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = (nproc_y - 1)* nproc_x
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(v(2,ny,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(v(nx-1,ny,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          v(nx-1,ny,k)=v(2,ny,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(nbc+ebc /= 0)) THEN

    DO k=2,nz-2
      v(nx-1,ny,k)=v(nx-1,ny,k)+vdtnb(nx-1,k)*dtsml
    END DO

  ELSE IF(nbc == 1.AND.ebc == 4) THEN

    DO k=2,nz-2
      v(nx-1,ny,k)=-v(nx-1,ny-2,k)
    END DO

  ELSE IF((nbc == 3.OR.nbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      v(nx-1,ny,k)=v(nx-1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      v(nx-1,ny,k)=v(nx-2,ny,k)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the top boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny
      DO i=1,nx-1
        v(i,j,nz-1)=v(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny
      DO i=1,nx-1
        v(i,j,nz-1)=v(i,j,2)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN  ! Zero normal gradient condition.

    DO j=1,ny
      DO i=1,nx-1
        v(i,j,nz-1)=v(i,j,nz-2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCV', tbc
    CALL arpsstop ("arpstop called from bcv top bc",1)

  END IF

  5005  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary conditions
!
!-----------------------------------------------------------------------
!

  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny
      DO i=1,nx-1
        v(i,j,1)=v(i,j,2)
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny
      DO i=1,nx-1
        v(i,j,1)=v(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny
      DO i=1,nx-1
        v(i,j,1)=v(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCV', bbc
    CALL arpsstop ("arpstop called from bcv bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE bcv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LBCW                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE lbcw(nx,ny,nz,dtsml,                                         &
           w,wcont,wdteb,wdtwb,wdtnb,wdtsb,                             &
           ebc,wbc,nbc,sbc,                                             &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the lateral boundary conditions for the w-velocity component.
!  Please note that the values at the corner points may depend on
!  the order that e-w and n-s boundary conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/04/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/07/92 (M. Xue)
!  Modified to take in w at time tfuture only.
!
!  6/15/92 (M. Xue and H. Jin)
!  Implemented open boundary condition.
!
!  10/6/92 (MX)
!  Assignment of the boundary conditions at corner columns moved to
!  the front of the top/bottom condition assignment.
!
!  Loop bounds of k for the lateral boundary conditions corrected
!  to 2,nz-1.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2/19/95 (K. Brewster)
!  Separated the application of external and user-supplied BC
!  processing from radiation BC.  Corrected application of
!  rigid or zero-gradient conditions mixed with radiation
!  conditions at SE corner.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml    The small time step size (s)
!
!    w        Interior domain values of w-velocity at tfuture (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    wdteb    Time tendency of the w field at the east boundary
!    wdtwb    Time tendency of the w field at the west boundary
!    wdtnb    Time tendency of the w field at the north boundary
!    wdtsb    Time tendency of the w field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    w        The w-velocity over the entire domain at tfuture (m/s)
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

  REAL :: dtsml                ! The small time step size (s)

  REAL :: w     (nx,ny,nz)     ! Total w-velocity at tfuture (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at east boundary
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at west boundary
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at north boundary
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at south boundary
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!   6 for nested grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  INCLUDE 'mp.inc'

  INTEGER :: mptag, ierror
  INTEGER :: source, dest

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO k=2,nz-1
      DO j=1,ny-1
        w(1,j,k)=w(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-1
      DO j=1,ny-1
        w(1,j,k)=w(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3) THEN        ! Zero normal gradient condition.

    DO k=2,nz-1
      DO j=1,ny-1
        w(1,j,k)=w(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-1
      DO j=2,ny-2
        w(1,j,k)=w(1,j,k)+wdtwb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(wbc == 5 .OR. wbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-1
      DO j=1,ny-1
        w(1,j,k)=w(2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCW', wbc
    CALL arpsstop ("arpstop called from bcw west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-1
      DO j=1,ny-1
        w(nx-1,j,k)=w(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-1
      DO j=1,ny-1
        w(nx-1,j,k)=w(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-1
      DO j=1,ny-1
        w(nx-1,j,k)=w(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-1
      DO j=2,ny-2
        w(nx-1,j,k)=w(nx-1,j,k)+wdteb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(ebc == 5 .OR. ebc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-1
      DO j=1,ny-1
        w(nx-1,j,k)=w(nx-2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCW', ebc
    CALL arpsstop ("arpstop called from bcw east bc",1)

  END IF

  5002  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-1
      DO i=1,nx-1
        w(i,ny-1,k)=w(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-1
      DO i=1,nx-1
        w(i,ny-1,k)=w(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-1
      DO i=1,nx-1
        w(i,ny-1,k)=w(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-1
      DO i=2,nx-2
        w(i,ny-1,k)=w(i,ny-1,k)+wdtnb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(nbc == 5 .OR. nbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-1
      DO i=1,nx-1
        w(i,ny-1,k)=w(i,ny-2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCW', nbc
    CALL arpsstop ("arpstop called from bcw north bc",1)

  END IF

  5003  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-1
      DO i=1,nx-1
        w(i,1,k)=w(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-1
      DO i=1,nx-1
        w(i,1,k)=w(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-1
      DO i=1,nx-1
        w(i,1,k)=w(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-1
      DO i=2,nx-2
        w(i,1,k)=w(i,1,k)+wdtsb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(sbc == 5 .OR.sbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-1
      DO i=1,nx-1
        w(i,1,k)=w(i,2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCW', sbc
    CALL arpsstop ("arpstop called from bcw south bc",1)

  END IF

  5004  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southwest corner based on the
!  boundary condition types on the south and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = (nproc_y -1)*nproc_x
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          w(1,1,k)=w(1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_x - 1
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(nx-2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          w(1,1,k)=w(nx-2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(sbc+wbc /= 0)) THEN

    DO k=2,nz-2
      w(1,1,k)=w(1,1,k)+wdtwb(1,k)*dtsml
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      w(1,1,k)=w(1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      w(1,1,k)=w(2,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southeast corner based on the
!  boundary condition types on the south and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_y*nproc_x  -1
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(nx-1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          w(nx-1,1,k)=w(nx-1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = 0
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          w(nx-1,1,k)=w(2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(sbc+ebc /= 0)) THEN

    DO k=2,nz-2
      w(nx-1,1,k)=w(nx-1,1,k)+wdteb(1,k)*dtsml
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      w(nx-1,1,k)=w(nx-1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      w(nx-1,1,k)=w(nx-2,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northwest corner based on the
!  boundary condition types on the north and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = 0
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          w(1,ny-1,k)=w(1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_y*nproc_x - 1
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(nx-2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          w(1,ny-1,k)=w(nx-2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(nbc+wbc /= 0)) THEN

    DO k=2,nz-2
      w(1,ny-1,k)=w(1,ny-1,k)+wdtwb(ny-1,k)*dtsml
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      w(1,ny-1,k)=w(1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      w(1,ny-1,k)=w(2,ny-1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northeast corner based on the
!  boundary condition types on the north and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_x - 1
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(nx-1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(nx-1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          w(nx-1,ny-1,k)=w(nx-1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = (nproc_y - 1)* nproc_x
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(w(2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(w(nx-1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          w(nx-1,ny-1,k)=w(2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(nbc+ebc /= 0)) THEN

    DO k=2,nz-2
      w(nx-1,ny-1,k)=w(nx-1,ny-1,k)+wdteb(ny-1,k)*dtsml
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      w(nx-1,ny-1,k)=w(nx-1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      w(nx-1,ny-1,k)=w(nx-2,ny-1,k)
    END DO

  END IF

!  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE lbcw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VBCW                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,                           &
           rhostr,rhostru,rhostrv,rhostrw,                              &
           j1,j2,j3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary conditions for w.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/04/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/07/92 (M. Xue)
!  Modified to take in w at time tfuture only.
!
!  6/15/92 (M. Xue and H. Jin)
!  Implemented open boundary condition
!
!  10/6/92 (MX)
!  Assignment of the boundary conditions at corner columns moved to
!  the front of the top/bottom condition assignment.
!
!  Loop bounds of k for the lateral boundary conditions corrected
!  to 2,nz-1.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2000/02/24 (Gene Bassett)
!  Fixed minor bug in tbc=1.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    w        Interior domain values of w-velocity at tfuture (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!    u        u-velocity at tfuture (m/s)
!    v        v-velocity at tfuture (m/s)
!
!    rhostr   j3 times base state density rhobar (kg/m**3).
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!    j1       Coordinate transformation Jacobian defined as
!             - d( zp )/d( x ).
!    j2       Coordinate transformation Jacobian defined as
!             - d( zp )/d( y ).
!    j3       Coordinate transformation Jacobian defined as
!             d( zp )/d( z ).
!
!  OUTPUT:
!
!    w        The w-velocity over the entire domain at tfuture (m/s)
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

  REAL :: w     (nx,ny,nz)     ! Total w-velocity at tfuture (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: u     (nx,ny,nz)     ! u-velocity at tfuture (m/s)
  REAL :: v     (nx,ny,nz)     ! v-velocity at tfuture (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!   6 for nested grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: urho1,urho2,vrho1,vrho2,wrho1,tems
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
!  Set the top boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny-1
      DO i=1,nx-1
        urho1=0.5*(u(i,j,nz-2)*rhostru(i,j,nz-2)                        &
              +u(i,j,nz-3)*rhostru(i,j,nz-3))
        urho2=0.5*(u(i+1,j,nz-2)*rhostru(i+1,j,nz-2)                    &
              +u(i+1,j,nz-3)*rhostru(i+1,j,nz-3))
        tems=0.5*(urho1*j1(i,j,nz-2)+urho2*j1(i+1,j,nz-2))

        vrho1=0.5*(v(i,j,nz-2)*rhostrv(i,j,nz-2)                        &
              +v(i,j,nz-3)*rhostrv(i,j,nz-3))
        vrho2=0.5*(v(i,j+1,nz-2)*rhostrv(i,j+1,nz-2)                    &
              +v(i,j+1,nz-3)*rhostrv(i,j+1,nz-3))
        tems=tems+0.5*(vrho1*j2(i,j,nz-2)+vrho2*j2(i,j+1,nz-2))

        wrho1=0.5*(j3(i,j,nz-2)+j3(i,j,nz-3))*rhostrw(i,j,nz-2)         &
              *wcont(i,j,nz-2)

!-----------------------------------------------------------------------
!
!  Note wrho1 above is calculated for k=nz-2, wrho1(nz)=-wrho1(nz-2)
!  based on rigid lid condition.
!
!-----------------------------------------------------------------------

        w(i,j,nz)= (- wrho1 - tems)/rhostrw(i,j,nz-2)
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        w(i,j,nz-1)= 0.0
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        w(i,j,nz)=w(i,j,3)
      END DO
    END DO

  ELSE IF(tbc == 3 .OR. tbc == 4) THEN ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        w(i,j,nz)=w(i,j,nz-1)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCW', tbc
    CALL arpsstop ("arpstop called from bcw top bc",1)

  END IF

  5005  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary conditions
!
!-----------------------------------------------------------------------
!

  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny-1
      DO i=1,nx-1
        urho1=0.5*(u(i,j,2)*rhostru(i,j,2)+u(i,j,3)*rhostru(i,j,3))
        urho2=0.5*(u(i+1,j,2)*rhostru(i+1,j,2)                          &
              +u(i+1,j,3)*rhostru(i+1,j,3))
        tems=0.5*(urho1*j1(i,j,3)+urho2*j1(i+1,j,3))

        vrho1=0.5*(v(i,j,2)*rhostrv(i,j,2)+v(i,j,3)*rhostrv(i,j,3))
        vrho2=0.5*(v(i,j+1,2)*rhostrv(i,j+1,2)                          &
              +v(i,j+1,3)*rhostrv(i,j+1,3))
        tems=tems+0.5*(vrho1*j2(i,j,3)+vrho2*j2(i,j+1,3))

        wrho1=0.5*(j3(i,j,2)+j3(i,j,3))*rhostrw(i,j,3)*wcont(i,j,3)

!-----------------------------------------------------------------------
!
!  Note wrho1 above is calculated for k=3, wrho1(1)=-wrho1(3)
!  based on non-penetrative ground condition.
!
!-----------------------------------------------------------------------

        w(i,j,1)= ( -wrho1 - tems )/rhostrw(i,j,3)
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        w(i,j,2) =-((u(i  ,j,2)+u(i  ,j,1))*j1(i,j,2)                   &
                   +(u(i+1,j,2)+u(i+1,j,1))*j1(i+1,j,2)                 &
                   +(v(i,j  ,2)+v(i,j  ,1))*j2(i,j,2)                   &
                   +(v(i,j+1,2)+v(i,j+1,1))*j2(i,j+1,2)) *0.25

      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        w(i,j,1)=w(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        w(i,j,1)=w(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCW', bbc
    CALL arpsstop ("arpstop called from bcw bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE vbcw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCP                        ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcp(nx,ny,nz,dtsml,                                          &
           pprt, pdteb,pdtwb,pdtnb,pdtsb,pprtk1,                        &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for the perturbation pressure.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/07/92 (M. Xue)
!  Modified to take in pprt at time tfuture only.
!
!  6/15/92 (M. Xue and H. Jin)
!  Implemented open boundary condition
!
!  10/6/92 (MX)
!  Assignment of the boundary conditions at corner columns moved to
!  the front of the top/bottom condition assignment.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2/19/95 (K. Brewster)
!  Separated the application of external and user-supplied BC
!  processing from radiation BC.
!
!  5/10/95 (M. Xue)
!  Changed lower BC for pprt to pprt(1)=2*pprt(2)=pprt(3).
!
!  11/05/97 (D. Weber)
!  Added pprtk1 array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml    The small time step size (s)
!
!    pprt     Interior domain perturbation pressure at tfuture
!             (Pascal)
!
!    pdteb    Time tendency of the pprt field at the east boundary
!    pdtwb    Time tendency of the pprt field at the west boundary
!    pdtnb    Time tendency of the pprt field at the north boundary
!    pdtsb    Time tendency of the pprt field at the south boundary
!
!    pprtk1   Perturbation pressure at k=1 computed using the
!             perturbation hydrostatic relation from the model
!             w-equation.
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    pprt     Purturbation pressure over the entire domain at tfuture
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

  REAL :: dtsml                ! The small time step size (s)

  REAL :: pprt  (nx,ny,nz)     ! Perturbation presure at tfuture
                               ! (Pascal)

  REAL :: pdteb (ny,nz)        ! Time tendency of pprt field at east
                               ! boundary
  REAL :: pdtwb (ny,nz)        ! Time tendency of pprt field at west
                               ! boundary
  REAL :: pdtnb (nx,nz)        ! Time tendency of pprt field at north
                               ! boundary
  REAL :: pdtsb (nx,nz)        ! Time tendency of pprt field at south
                               ! boundary

  REAL :: pprtk1(nx,ny)        ! Perturbation pressure at k=1 computed
                               ! using the perturbation hydrostatic
                               ! relation from the model w-equation.

!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!   6 for nested grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k

  INCLUDE 'mp.inc'

  INTEGER :: mptag, ierror
  INTEGER :: source, dest

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(1,j,k)=pprt(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        pprt(1,j,k)=pprt(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3) THEN        ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(1,j,k)=pprt(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 4) THEN
                                   ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-2
        pprt(1,j,k)=pprt(1,j,k)+pdtwb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(wbc == 5 .OR. wbc == 6) THEN
                                   ! External or user specified condition.

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(1,j,k)=pprt(2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCP', wbc
    CALL arpsstop ("arpstop called from bcp west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(nx-1,j,k)=pprt(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        pprt(nx-1,j,k)=pprt(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(nx-1,j,k)=pprt(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-2
        pprt(nx-1,j,k)=pprt(nx-1,j,k)+pdteb(j,k)*dtsml
      END DO
    END DO

  ELSE IF(ebc == 5 .OR. ebc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny-1
        pprt(nx-1,j,k)=pprt(nx-2,j,k)
      END DO
    END DO

  ELSE
    WRITE(6,900) 'BCP', ebc
    CALL arpsstop ("arpstop called from bcp east bc",1)

  END IF

  5002  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,ny-1,k)=pprt(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,ny-1,k)=pprt(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,ny-1,k)=pprt(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-2
      DO i=2,nx-2
        pprt(i,ny-1,k)=pprt(i,ny-1,k)+pdtnb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(nbc == 5 .OR. nbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,ny-1,k)=pprt(i,ny-2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCP', nbc
    CALL arpsstop ("arpstop called from bcp north bc",1)

  END IF

  5003  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,1,k)=pprt(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,1,k)=pprt(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,1,k)=pprt(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-2
      DO i=2,nx-2
        pprt(i,1,k)=pprt(i,1,k)+pdtsb(i,k)*dtsml
      END DO
    END DO

  ELSE IF(sbc == 5 .OR.sbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx-1
        pprt(i,1,k)=pprt(i,2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCP', sbc
    CALL arpsstop ("arpstop called from bcp south bc",1)

  END IF

  5004  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southwest corner based on the
!  boundary condition types on the south and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = (nproc_y -1)*nproc_x
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          pprt(1,1,k)=pprt(1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_x - 1
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(nx-2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          pprt(1,1,k)=pprt(nx-2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(sbc+wbc /= 0)) THEN

    DO k=2,nz-2
      pprt(1,1,k)=pprt(1,1,k)+pdtwb(1,k)*dtsml
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      pprt(1,1,k)=pprt(1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      pprt(1,1,k)=pprt(2,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southeast corner based on the
!  boundary condition types on the south and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_y*nproc_x  -1
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(nx-1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          pprt(nx-1,1,k)=pprt(nx-1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = 0
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          pprt(nx-1,1,k)=pprt(2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(sbc+ebc /= 0)) THEN

    DO k=2,nz-2
      pprt(nx-1,1,k)=pprt(nx-1,1,k)+pdteb(1,k)*dtsml
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      pprt(nx-1,1,k)=pprt(nx-1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      pprt(nx-1,1,k)=pprt(nx-2,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northwest corner based on the
!  boundary condition types on the north and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = 0
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          pprt(1,ny-1,k)=pprt(1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_y*nproc_x - 1
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(nx-2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          pprt(1,ny-1,k)=pprt(nx-2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(nbc+wbc /= 0)) THEN

    DO k=2,nz-2
      pprt(1,ny-1,k)=pprt(1,ny-1,k)+pdtwb(ny-1,k)*dtsml
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      pprt(1,ny-1,k)=pprt(1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      pprt(1,ny-1,k)=pprt(2,ny-1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northeast corner based on the
!  boundary condition types on the north and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_x - 1
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(nx-1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(nx-1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          pprt(nx-1,ny-1,k)=pprt(nx-1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = (nproc_y - 1)* nproc_x
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(pprt(2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(pprt(nx-1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          pprt(nx-1,ny-1,k)=pprt(2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(nbc+ebc /= 0)) THEN

    DO k=2,nz-2
      pprt(nx-1,ny-1,k)=pprt(nx-1,ny-1,k)+pdteb(ny-1,k)*dtsml
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      pprt(nx-1,ny-1,k)=pprt(nx-1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      pprt(nx-1,ny-1,k)=pprt(nx-2,ny-1,k)
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Set the top boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,nz-1)=pprt(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,nz-1)=pprt(i,j,2)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,nz-1)=pprt(i,j,nz-2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCP', tbc
    CALL arpsstop ("arpstop called from bcp top bc",1)

  END IF

  5005  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary conditions
!
!-----------------------------------------------------------------------
!

  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,1)=pprtk1(i,j)                ! hydrostatic pprt.
!      pprt(i,j,1)=2*pprt(i,j,2)-pprt(i,j,3)  ! Cst gradient extrapolation
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,1)=pprt(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,1)=pprt(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCP', bbc
    CALL arpsstop ("arpstop called from bcp bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE bcp
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCPT                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcpt(nx,ny,nz,dtbig,                                         &
           ptprt, ptdteb,ptdtwb,ptdtnb,ptdtsb,                          &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for the potential temperature perturbation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/04/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig    The large time step size (s)
!
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!
!    ptdteb   Time tendency of the ptprt field at the east boundary
!    ptdtwb   Time tendency of the ptprt field at the west boundary
!    ptdtnb   Time tendency of the ptprt field at the north boundary
!    ptdtsb   Time tendency of the ptprt field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    ptprt    ptprt over the entire domain at tfuture (K)
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

  REAL :: dtbig                ! The big time step size (s)

  INCLUDE 'timelvls.inc'

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)

  REAL :: ptdteb(ny,nz)        ! Time tendency of ptprt at east
                               ! boundary
  REAL :: ptdtwb(ny,nz)        ! Time tendency of ptprt at west
                               ! boundary
  REAL :: ptdtnb(nx,nz)        ! Time tendency of ptprt at north
                               ! boundary
  REAL :: ptdtsb(nx,nz)        ! Time tendency of ptprt at south
                               ! boundary
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top   boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
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
!  Set the boundary conditions for ptprt.
!
!-----------------------------------------------------------------------
!
  CALL bcsclr(nx,ny,nz,dtbig,                                           &
              ptprt(1,1,1,tpast),ptprt(1,1,1,tpresent),                 &
              ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,         &
              ebc,wbc,nbc,sbc,tbc,bbc,                                  &
              ebc_global,wbc_global,nbc_global,sbc_global)

  RETURN
END SUBROUTINE bcpt

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCQ                        ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcq(nx,ny,nz,dtbig,                                          &
           q, qdteb,qdtwb,qdtnb,qdtsb,                                  &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for all of the water quantities
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  6/04/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig    The large time step size (s)
!
!    q        mixing ratio for one of the water variables at all
!             time levels (kg/kg)
!
!    qdteb    Time tendency of the q field at the east boundary
!    qdtwb    Time tendency of the q field at the west boundary
!    qdtnb    Time tendency of the q field at the north boundary
!    qdtsb    Time tendency of the q field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    q        Array q over the entire domain at tfuture (m/s)
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

  REAL :: dtbig                ! Big time step size (s)

  INCLUDE 'timelvls.inc'

  REAL :: q     (nx,ny,nz,nt)  ! mixing ratio for one of the water/ice
                               ! variables (kg/kg)

  REAL :: qdteb(ny,nz)         ! Time tendency of q at east boundary
  REAL :: qdtwb(ny,nz)         ! Time tendency of q at west boundary
  REAL :: qdtnb(nx,nz)         ! Time tendency of q at north boundary
  REAL :: qdtsb(nx,nz)         ! Time tendency of q at south boundary
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
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
!  Set the boundary conditions for q.
!
!-----------------------------------------------------------------------
!
  CALL bcsclr(nx,ny,nz,dtbig,                                           &
              q(1,1,1,tpast),q(1,1,1,tpresent),q(1,1,1,tfuture),        &
              qdteb,qdtwb,qdtnb,qdtsb,                                  &
              ebc,wbc,nbc,sbc,tbc,bbc,                                  &
              ebc_global,wbc_global,nbc_global,sbc_global)

  RETURN
END SUBROUTINE bcq
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCQV                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcqv(nx,ny,nz,dtbig,                                         &
           qv,qvbar,qdteb,qdtwb,qdtnb,qdtsb,                            &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for qv.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  06/27/95
!  Created subroutine BCQV, based on BCQ. BCQV handles special
!  top and bottom boundary conditions for qv.
!
!
!  MODIFICATION HISTORY:
!
!  07/06/95 (Ming Xue)
!  Changed the top and bottom BC for qvprt back to zero gradient.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig    The large time step size (s)
!
!    qv       Water vapor mixing ratio (kg/kg)
!    qvbar    Base-state water vapor mixing ratio (kg/kg)
!
!    qdteb    Time tendency of the q field at the east boundary
!    qdtwb    Time tendency of the q field at the west boundary
!    qdtnb    Time tendency of the q field at the north boundary
!    qdtsb    Time tendency of the q field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    qv       qv over the entire domain at tfuture (m/s)
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

  REAL :: dtbig                ! Big time step size (s)

  INCLUDE 'timelvls.inc'

  REAL :: qv    (nx,ny,nz,nt)  ! water vapor mixing ratio (kg/kg)
  REAL :: qvbar (nx,ny,nz)     ! base-state water vapro mixing ratio (kg/kg)

  REAL :: qdteb(ny,nz)         ! Time tendency of q at east boundary
  REAL :: qdtwb(ny,nz)         ! Time tendency of q at west boundary
  REAL :: qdtnb(nx,nz)         ! Time tendency of q at north boundary
  REAL :: qdtsb(nx,nz)         ! Time tendency of q at south boundary
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
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
!  Set the boundary conditions for qv.
!
!-----------------------------------------------------------------------
!
  CALL bcsclr(nx,ny,nz,dtbig,                                           &
              qv(1,1,1,tpast),qv(1,1,1,tpresent),qv(1,1,1,tfuture),     &
              qdteb,qdtwb,qdtnb,qdtsb,                                  &
              ebc,wbc,nbc,sbc,tbc,bbc,                                  &
              ebc_global,wbc_global,nbc_global,sbc_global)

  RETURN
END SUBROUTINE bcqv

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BCSCLR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcsclr(nx,ny,nz,dtbig,                                       &
           s1,s2,s3,sdteb,sdtwb,sdtnb,sdtsb,                            &
           ebc,wbc,nbc,sbc,tbc,bbc,                                     &
           ebc_global,wbc_global,nbc_global,sbc_global)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a scalar s. Please note that the
!  values at the corner points may depend on the order that the e-w,
!  n-s and t-b boundary conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/01/92 (M. Xue)
!  Added full documentation.
!
!  10/6/92 (MX)
!  Assignment of the boundary conditions at corner columns moved to
!  the front of the top/bottom condition assignment.
!
!  2/1/93 (Hao Jin)
!  Assignment of the boundary conditions at corner columns modified
!  to take care of all possible combination of options.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2/19/95 (K. Brewster)
!  Separated the application of external and user-supplied BC
!  processing from radiation BC.  Corrected application of
!  rigid or zero-gradient conditions mixed with radiation
!  conditions at SE corner.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig    The large time step size (s)
!
!    s1       A scalar variable at time tpast
!    s2       A scalar variable at time tpresent
!    s3       A scalar variable at time tfuture
!
!    sdteb    Time tendency of the s field at the east boundary
!    sdtwb    Time tendency of the s field at the west boundary
!    sdtnb    Time tendency of the s field at the north boundary
!    sdtsb    Time tendency of the s field at the south boundary
!
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!    ebc_global      Global value of ebc.
!    wbc_global      Global value of wbc.
!    nbc_global      Global value of nbc.
!    sbc_global      Global value of sbc.
!                    They will be the same as ebc,wbc, nbc and sbc
!                    for one processor serial runs.
!
!  OUTPUT:
!
!    s3       Scalar array s3 over the entire domain tfuture
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

  REAL :: dtbig                ! The big time step size (s)

  REAL :: s1    (nx,ny,nz)     ! A scalar variable at time tpast.
  REAL :: s2    (nx,ny,nz)     ! A scalar variable at time tpresent.
  REAL :: s3    (nx,ny,nz)     ! A scalar variable at time tfuture.

  REAL :: sdteb(ny,nz)         ! Time tendency of s field at east
                               ! boundary
  REAL :: sdtwb(ny,nz)         ! Time tendency of s field at west
                               ! boundary
  REAL :: sdtnb(nx,nz)         ! Time tendency of s field at north
                               ! boundary
  REAL :: sdtsb(nx,nz)         ! Time tendency of s field at south
                               ! boundary

!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.

  INTEGER :: ebc_global        ! Global value of ebc
  INTEGER :: wbc_global        ! Global value of wbc
  INTEGER :: nbc_global        ! Global value of nbc.
  INTEGER :: sbc_global        ! Global value of sbc.
                               ! They are idential to ebc,wbc,nbc, & sbc
                               ! for non-mpi runs.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  INCLUDE 'mp.inc'

  INTEGER :: mptag, ierror
  INTEGER :: source, dest
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
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny-1
        s3(1,j,k)=s3(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        s3(1,j,k)=s3(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3) THEN        ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny-1
        s3(1,j,k)=s3(2,j,k)
      END DO
    END DO

  ELSE IF(wbc == 4) THEN
                                   ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-2
        s3(1,j,k)=s1(1,j,k)+sdtwb(j,k)*2.*dtbig
      END DO
    END DO

  ELSE IF(wbc == 5 .OR. wbc == 6) THEN
                                   ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny-1
        s3(1,j,k)=s3(2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSCLR', wbc
    CALL arpsstop ("arpstop called from bcsclr west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO j=1,ny-1
        s3(nx-1,j,k)=s3(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        s3(nx-1,j,k)=s3(2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO j=1,ny-1
        s3(nx-1,j,k)=s3(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 4) THEN
                                    ! Radiation condition.
    DO k=2,nz-2
      DO j=2,ny-2
        s3(nx-1,j,k)=s1(nx-1,j,k)+sdteb(j,k)*2.*dtbig
      END DO
    END DO

  ELSE IF(ebc == 5 .OR. ebc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO j=1,ny-1
        s3(nx-1,j,k)=s3(nx-2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSCLR', ebc
    CALL arpsstop ("arpstop called from bcsclr east bc",1)

  END IF

  5002  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx-1
        s3(i,ny-1,k)=s3(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        s3(i,ny-1,k)=s3(i,2,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx-1
        s3(i,ny-1,k)=s3(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-2
      DO i=2,nx-2
        s3(i,ny-1,k)=s1(i,ny-1,k)+sdtnb(i,k)*2.*dtbig
      END DO
    END DO

  ELSE IF(nbc == 5 .OR. nbc == 6 ) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx-1
        s3(i,ny-1,k)=s3(i,ny-2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSCLR', nbc
    CALL arpsstop ("arpstop called from bcsclr north bc",1)

  END IF

  5003  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=2,nz-2
      DO i=1,nx-1
        s3(i,1,k)=s3(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        s3(i,1,k)=s3(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3) THEN         ! Zero normal gradient condition.

    DO k=2,nz-2
      DO i=1,nx-1
        s3(i,1,k)=s3(i,2,k)
      END DO
    END DO

  ELSE IF(sbc == 4) THEN         ! Radiation condition.

    DO k=2,nz-2
      DO i=2,nx-2
        s3(i,1,k)=s1(i,1,k)+sdtsb(i,k)*2.*dtbig
      END DO
    END DO

  ELSE IF(sbc == 5 .OR. sbc == 6) THEN
                                    ! External or user specified condition.
    DO k=2,nz-2
      DO i=1,nx-1
        s3(i,1,k)=s3(i,2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSCLR', sbc
    CALL arpsstop ("arpstop called from bcsclr south bc",1)

  END IF

  5004  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southwest corner based on the
!  boundary condition types on the south and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = (nproc_y -1)*nproc_x
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          s3(1,1,k)=s3(1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_x - 1
      dest   = 0
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(nx-2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          s3(1,1,k)=s3(nx-2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(sbc+wbc /= 0)) THEN

    DO k=2,nz-2
      s3(1,1,k)=s1(1,1,k)+sdtwb(1,k)*2.*dtbig
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      s3(1,1,k)=s3(1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      s3(1,1,k)=s3(2,1,k)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the southeast corner based on the
!  boundary condition types on the south and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (sbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (sbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(sbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_y*nproc_x  -1
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(nx-1,ny-2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          s3(nx-1,1,k)=s3(nx-1,ny-2,k)
        END DO
      END IF
    END IF

  ELSE IF(sbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = 0
      dest   = nproc_x -1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(2,1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(nx-1,1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          s3(nx-1,1,k)=s3(2,1,k)
        END DO
      END IF
    END IF
  END IF

  IF((sbc == 4.OR.sbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(sbc+ebc /= 0)) THEN

    DO k=2,nz-2
      s3(nx-1,1,k)=s1(nx-1,1,k)+sdteb(1,k)*2.*dtbig
    END DO

  ELSE IF((sbc == 1.OR.sbc == 3.OR.sbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      s3(nx-1,1,k)=s3(nx-1,2,k)
    END DO

  ELSE IF(sbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      s3(nx-1,1,k)=s3(nx-2,1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northwest corner based on the
!  boundary condition types on the north and west boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. wbc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. wbc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.wbc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = 0
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == 0) THEN
        DO k=2,nz-2
          s3(1,ny-1,k)=s3(1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.wbc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = nproc_y*nproc_x - 1
      dest   = (nproc_y-1)*nproc_x
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(nx-2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          s3(1,ny-1,k)=s3(nx-2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(wbc == 4.OR.wbc == 0).AND.(nbc+wbc /= 0)) THEN

    DO k=2,nz-2
      s3(1,ny-1,k)=s1(1,ny-1,k)+sdtwb(ny-1,k)*2.*dtbig
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.wbc == 4) THEN

    DO k=2,nz-2
      s3(1,ny-1,k)=s3(1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(wbc == 1.OR.wbc == 3.OR.wbc == 5)) THEN

    DO k=2,nz-2
      s3(1,ny-1,k)=s3(2,ny-1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions at the northeast corner based on the
!  boundary condition types on the north and east boundaries.
!
!-----------------------------------------------------------------------
!
! Ensure all processors increase gentag
!
  IF ( (nbc_global == 2 .AND. ebc_global == 4) .OR.                   &
       (nbc_global == 4 .AND. ebc_global == 2)  ) THEN
    CALL inctag
    mptag = gentag
  END IF

  IF(nbc_global == 2.AND.ebc == 4) THEN

    IF (mp_opt > 0 .AND. nproc_y > 1) THEN
      source = nproc_x - 1
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(nx-1,2,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(nx-1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_x-1) THEN
        DO k=2,nz-2
          s3(nx-1,ny-1,k)=s3(nx-1,2,k)
        END DO
      END IF
    END IF

  ELSE IF(nbc == 4.AND.ebc_global == 2) THEN

    IF (mp_opt > 0 .AND. nproc_x > 1) THEN
      source = (nproc_y - 1)* nproc_x
      dest   = nproc_y*nproc_x - 1
      DO k = 2, nz-2
        mptag = gentag + k
        IF(myproc == source) CALL mpsendr(s3(2,ny-1,k),1,dest,mptag,ierror)
        IF(myproc == dest)   CALL mprecvr(s3(nx-1,ny-1,k),1,source,mptag,ierror)
      END DO
    ELSE
      IF (myproc == nproc_y-1) THEN
        DO k=2,nz-2
          s3(nx-1,ny-1,k)=s3(2,ny-1,k)
        END DO
      END IF
    END IF
  END IF

  IF((nbc == 4.OR.nbc == 0).AND.(ebc == 4.OR.ebc == 0).AND.(nbc+ebc /= 0)) THEN

    DO k=2,nz-2
      s3(nx-1,ny-1,k)=s1(nx-1,ny-1,k)+sdteb(ny-1,k)*2.*dtbig
    END DO

  ELSE IF((nbc == 1.OR.nbc == 3.OR.nbc == 5).AND.ebc == 4) THEN

    DO k=2,nz-2
      s3(nx-1,ny-1,k)=s3(nx-1,ny-2,k)
    END DO

  ELSE IF(nbc == 4.AND.(ebc == 1.OR.ebc == 3.OR.ebc == 5)) THEN

    DO k=2,nz-2
      s3(nx-1,ny-1,k)=s3(nx-2,ny-1,k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the top boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,nz-1)=s3(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,nz-1)=s3(i,j,2)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN  ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,nz-1)=s3(i,j,nz-2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSCLR', tbc
    CALL arpsstop ("arpstop called from bcsclr top bc",1)

  END IF

  5005  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary conditions
!
!-----------------------------------------------------------------------
!

  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,1)=s3(i,j,2)
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,1)=s3(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,1)=s3(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSCLR', bbc
    CALL arpsstop ("arpstop called from bcsclr bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE bcsclr
!

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCSU                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcsu(nx,ny,nz,jbgn,jend,kbgn,kend,ebc,wbc,su)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the x-boundary conditions for su, an array that has been averaged
!  from scalar points to u-points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  5/22/93 (D Weber & MX)
!  Added index bounds into the argument list.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    jbgn     Starting j index
!    jend     Ending j index
!    kbgn     Starting k index
!    kend     Ending k index
!
!  OUTPUT:
!
!    su       Output array at u point
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: jbgn,jend         ! Domain for j computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: su(nx,ny,nz)         ! Averaged array

  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j,k
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
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO j=jbgn,jend
        su(1,j,k)=su(3,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        su(1,j,k)=su(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3 .OR. wbc == 4 .OR. wbc == 5 .OR. wbc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO j=jbgn,jend
        su(1,j,k)=su(2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSU', wbc
    CALL arpsstop ("arpstop called from bcsu west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO j=jbgn,jend
        su(nx,j,k)=su(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        su(nx,j,k)=su(3,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc == 3 .OR. ebc == 4 .OR. ebc == 5 .OR. ebc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO j=jbgn,jend
        su(nx,j,k)=su(nx-1,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSU', ebc
    CALL arpsstop ("arpstop called from bcsu east bc",1)

  END IF

  5002  CONTINUE

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')

  RETURN
END SUBROUTINE bcsu
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCSV                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcsv(nx,ny,nz,ibgn,iend,kbgn,kend,nbc,sbc,sv)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the y-boundary conditions for sv, an array that has been averaged
!  from scalar points to v-points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  5/22/93 (D Weber & MX)
!  Added index bounds into the argument list.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Starting i index
!    iend     Ending i index
!    kbgn     Starting k index
!    kend     Ending k index
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!
!  OUTPUT:
!
!    sv       Output array at u point
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions

  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: kbgn,kend         ! Domain for k computations

  REAL :: sv(nx,ny,nz)         ! Averaged array

  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
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
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO i=ibgn,iend
        sv(i,ny,k)=sv(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        sv(i,ny,k)=sv(i,3,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3 .OR. nbc == 4 .OR. nbc == 5 .OR. nbc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO i=ibgn,iend
        sv(i,ny,k)=sv(i,ny-1,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSV', nbc
    CALL arpsstop ("arpstop called from bcsv north bc",1)

  END IF

  5003  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO i=ibgn,iend
        sv(i,1,k)=sv(i,3,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        sv(i,1,k)=sv(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3 .OR. sbc == 4 .OR. sbc == 5 .OR. sbc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO i=ibgn,iend
        sv(i,1,k)=sv(i,2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSV', sbc
    CALL arpsstop ("arpstop called from bcsv south bc",1)

  END IF

  5004  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE bcsv

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCSW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcsw(nx,ny,nz,ibgn,iend,jbgn,jend,tbc,bbc,sw)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary condition for sw, an array that has
!  been averaged from scalar points to w-points..
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  5/22/93 (D Weber & MX)
!  Added index bounds into the argument list.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Starting i index
!    iend     Ending i index
!    jbgn     Starting j index
!    jend     Ending j index
!
!    tbc      Parameter defining top    boundary condition type.
!    bbc      Parameter defining bottom boundary condition type.
!
!  OUTPUT:
!
!    sw       Averaged array at w point
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions

  INTEGER :: ibgn,iend         ! Domain for i computations
  INTEGER :: jbgn,jend         ! Domain for j computations

  REAL :: sw(nx,ny,nz)         ! Averaged array

  INTEGER :: tbc               ! Parameter defining top    boundary
                               ! condition type.
  INTEGER :: bbc               ! Parameter defining bottom boundary
                               ! condition type.
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
!  Set the top boundary condition
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        sw(i,j,nz)=sw(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        sw(i,j,nz)=sw(i,j,3)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN ! Zero normal gradient condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        sw(i,j,nz)=sw(i,j,nz-1)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSW', tbc
    CALL arpsstop ("arpstop called from bcsw top bc",1)

  END IF

  5005  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary condition
!
!-----------------------------------------------------------------------
!
  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        sw(i,j,1)=sw(i,j,3)
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        sw(i,j,1)=sw(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        sw(i,j,1)=sw(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BCSW', bbc
    CALL arpsstop ("arpstop called from bcsw bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE bcsw

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE VBCWCONT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vbcwcont(nx,ny,nz,wcont)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the top and bottom boundary conditions for wcont.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    wcont    Contravariant vertical velocity (m/s)
!
!  OUTPUT:
!
!    wcont    Top and bottom values of wcont.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions

  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
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
!  Set the top boundary condition
!
!-----------------------------------------------------------------------
!
  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,nz)=-wcont(i,j,nz-2)
        wcont(i,j,nz-1)=0.0
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,nz)=wcont(i,j,3)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,nz)=wcont(i,j,nz-1)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'VBCWCONT', tbc
    CALL arpsstop ("arpstop called from vbcwcont top bc",1)

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary condition
!
!-----------------------------------------------------------------------
!
  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,1)=-wcont(i,j,3)
        wcont(i,j,2)=0.0
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,1)=wcont(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,1)=wcont(i,j,2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'VBCWCONT', bbc
    CALL arpsstop ("arpstop called from vbcwcont bottom bc",1)

  END IF

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE vbcwcont
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCU2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcu2d(nx,ny,a, ebc,wbc,nbc,sbc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a 2-D variable at a u-velocity
!  location.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2000/02/28 (Gene Bassett)
!  Fixed bugs where scaler locations, not staggered locations, were
!  used.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    s        A 2-D scalar whose boundary values are to be set.
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!
!  OUTPUT:
!
!    s        The boundary values of s.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny             ! Number of grid points in x, y
                               ! directions

  REAL :: a(nx,ny)             ! A scalar variable

  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO j=1,ny-1
      a(1,j)=a(3,j)
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny-1
        a(1,j)=a(nx-2,j)
      END DO
    END IF

  ELSE IF(wbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny-1
      a(1,j)=a(3,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO j=1,ny-1
      a(nx,j)=a(nx-2,j)
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny-1
        a(nx,j)=a(3,j)
      END DO
    END IF

  ELSE IF(ebc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny-1
      a(nx,j)=a(nx-2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx
      a(i,ny-1)=a(i,ny-2)
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx
        a(i,ny-1)=a(i,2)
      END DO
    END IF

  ELSE IF(nbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx
      a(i,ny-1)=a(i,ny-2)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx
      a(i,1)=a(i,2)
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx
        a(i,1)=a(i,ny-2)
      END DO
    END IF

  ELSE IF(sbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx
      a(i,1)=a(i,2)
    END DO

  END IF

  RETURN
END SUBROUTINE bcu2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCV2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcv2d(nx,ny,a, ebc,wbc,nbc,sbc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a 2-D variable at a v-velocity
!  location.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  2000/02/28 (Gene Bassett)
!  Fixed bugs where scaler locations, not staggered locations, were
!  used.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    s        A 2-D scalar whose boundary values are to be set.
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!
!  OUTPUT:
!
!    s        The boundary values of s.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny             ! Number of grid points in x, y
                               ! directions

  REAL :: a(nx,ny)             ! A scalar variable

  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO j=1,ny
      a(1,j)=a(2,j)
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny
        a(1,j)=a(nx-2,j)
      END DO
    END IF

  ELSE IF(wbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny
      a(1,j)=a(2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO j=1,ny
      a(nx-1,j)=a(nx-2,j)
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny
        a(nx-1,j)=a(2,j)
      END DO
    END IF

  ELSE IF(ebc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny
      a(nx-1,j)=a(nx-2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx-1
      a(i,ny)=a(i,ny-2)
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx-1
        a(i,ny)=a(i,3)
      END DO
    END IF

  ELSE IF(nbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx-1
      a(i,ny)=a(i,ny-2)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx-1
      a(i,1)=a(i,3)
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx-1
        a(i,1)=a(i,ny-2)
      END DO
    END IF

  ELSE IF(sbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx-1
      a(i,1)=a(i,3)
    END DO

  END IF

  RETURN
END SUBROUTINE bcv2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCS2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcs2d(nx,ny,s, ebc,wbc,nbc,sbc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a 2-D scalar.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    s        A 2-D scalar whose boundary values are to be set.
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!
!  OUTPUT:
!
!    s        The boundary values of s.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny             ! Number of grid points in x, y
                               ! directions

  REAL :: s     (nx,ny)        ! A scalar variable

  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO j=1,ny-1
      s(1,j)=s(2,j)
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny-1
        s(1,j)=s(nx-2,j)
      END DO
    END IF

  ELSE IF(wbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny-1
      s(1,j)=s(2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO j=1,ny-1
      s(nx-1,j)=s(nx-2,j)
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny-1
        s(nx-1,j)=s(2,j)
      END DO
    END IF

  ELSE IF(ebc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny-1
      s(nx-1,j)=s(nx-2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx-1
      s(i,ny-1)=s(i,ny-2)
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx-1
        s(i,ny-1)=s(i,2)
      END DO
    END IF

  ELSE IF(nbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx-1
      s(i,ny-1)=s(i,ny-2)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------
!
  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx-1
      s(i,1)=s(i,2)
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx-1
        s(i,1)=s(i,ny-2)
      END DO
    END IF

  ELSE IF(sbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx-1
      s(i,1)=s(i,2)
    END DO

  END IF

  RETURN
END SUBROUTINE bcs2d
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BCIS2D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bcis2d(nx,ny,is, ebc,wbc,nbc,sbc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a 2-D integer scalar array.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    is       A 2-D integer scalar whose boundary values are to be set.
!    ebc      Parameter defining east   boundary condition type.
!    wbc      Parameter defining west   boundary condition type.
!    nbc      Parameter defining north  boundary condition type.
!    sbc      Parameter defining south  boundary condition type.
!
!  OUTPUT:
!
!    is       The boundary values of is.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny             ! Number of grid points in x, y and z
                               ! directions

  INTEGER :: is (nx,ny)        ! A scalar variable

  INTEGER :: ebc               ! Parameter defining east   boundary
                               ! condition type.
  INTEGER :: wbc               ! Parameter defining west   boundary
                               ! condition type.
  INTEGER :: nbc               ! Parameter defining north  boundary
                               ! condition type.
  INTEGER :: sbc               ! Parameter defining south  boundary
                               ! condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
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
!  Set the west boundary condition
!
!-----------------------------------------------------------------------
!
  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO j=1,ny-1
      is(1,j)=is(2,j)
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny-1
        is(1,j)=is(nx-2,j)
      END DO
    END IF

  ELSE IF(wbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny-1
      is(1,j)=is(2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the east boundary condition
!
!-----------------------------------------------------------------------
!
  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO j=1,ny-1
      is(nx-1,j)=is(nx-2,j)
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO j=1,ny-1
        is(nx-1,j)=is(2,j)
      END DO
    END IF

  ELSE IF(ebc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO j=1,ny-1
      is(nx-1,j)=is(nx-2,j)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the north boundary condition
!
!-----------------------------------------------------------------------
!
  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx-1
      is(i,ny-1)=is(i,ny-2)
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx-1
        is(i,ny-1)=is(i,2)
      END DO
    END IF

  ELSE IF(nbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx-1
      is(i,ny-1)=is(i,ny-2)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the south boundary condition
!
!-----------------------------------------------------------------------
!
  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO i=1,nx-1
      is(i,1)=is(i,2)
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO i=1,nx-1
        is(i,1)=is(i,ny-2)
      END DO
    END IF

  ELSE IF(sbc /= 0) THEN        ! Zero normal gradient condition or
                                ! Radiation or user specified condition.
    DO i=1,nx-1
      is(i,1)=is(i,2)
    END DO

  END IF

  RETURN
END SUBROUTINE bcis2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BOUNDU                   ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE boundu(s,nx,ny,nz,jbgn,jend,kbgn,kend)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a quantity at the first and last u
!  points in the x-direction. Please note that the values at the corner
!  points may depend on the order that e-w, n-s and t-b boundary
!  conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  5/13/93 (D. Weber)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    jbgn     Index to start the j direction
!    jend     Index to end the j direction
!    kbgn     Index to start the k direction
!    kend     Index to end the k direction
!
!    s        Input array defined from i=2 to i=nx-1
!
!  OUTPUT:
!
!    s       Output array at u point for i=1 to i=nx
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions

  INTEGER :: jbgn,jend         ! Bounds for computations in the y
                               ! direction

  INTEGER :: kbgn,kend         ! Bounds for computations in the z
                               ! direction

  REAL :: s(nx,ny,nz)          ! A scalar variable in the x direction
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
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
!  Set the west boundary condition
!
!-----------------------------------------------------------------------
!
  IF(wbc == 0) GO TO 5001

  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO j=jbgn,jend
        s(1,j,k)=-s(3,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        s(1,j,k)=s(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc == 3 .OR. wbc == 4 .OR. wbc == 5 .OR. wbc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO j=jbgn,jend
        s(1,j,k)=s(2,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BOUNDU', wbc
    CALL arpsstop ("arpstop called from boundu west bc",1)

  END IF

  5001  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Set the east boundary condition
!
!-----------------------------------------------------------------------
!
  IF(ebc == 0) GO TO 5002

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO j=jbgn,jend
        s(nx,j,k)=-s(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO j=jbgn,jend
        s(nx,j,k)=s(3,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc == 3 .OR. ebc == 4 .OR. ebc == 5 .OR. ebc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO j=jbgn,jend
        s(nx,j,k)=s(nx-1,j,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BOUNDU', ebc
    CALL arpsstop ("arpstop called from boundu east bc",1)

  END IF

  5002  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE boundu

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BOUNDV                   ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE boundv(s,nx,ny,nz,ibgn,iend,kbgn,kend)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a quantity at the first and last v
!  points in the y-direction. Please note that the values at the corner
!  points may depend on the order that e-w, n-s and t-b boundary
!  conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  5/13/93 (D. Weber)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Index to start the i direction
!    iend     Index to end the i direction
!    kbgn     Index to start the k direction
!    kend     Index to end the k direction
!
!    s        Input array defined from j=2 to j=ny-1
!
!  OUTPUT:
!
!    s       Output array at v point for j=1 to j=ny
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions

  INTEGER :: ibgn,iend         ! Bounds for computations in the x
                               ! direction

  INTEGER :: kbgn,kend         ! Bounds for computations in the z
                               ! direction

  REAL :: s(nx,ny,nz)          ! A scalar variable in the y direction
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!   4 for open (radiation) boundary condition.
!   5 for user (externally) specified boundary condition.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
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
!  Set the north boundary condition
!
!-----------------------------------------------------------------------
!
  IF(nbc == 0) GO TO 5003

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO i=ibgn,iend
        s(i,ny,k)=-s(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        s(i,ny,k)=s(i,3,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc == 3 .OR. nbc == 4 .OR. nbc == 5 .OR. nbc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO i=ibgn,iend
        s(i,ny,k)=s(i,ny-1,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BOUNDV', nbc
    CALL arpsstop ("arpstop called from boundv north bc",1)

  END IF

  5003  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the south boundary condition
!
!-----------------------------------------------------------------------
!
  IF(sbc == 0) GO TO 5004

  IF(sbc == 1) THEN             ! Rigid wall boundary condition

    DO k=kbgn,kend
      DO i=ibgn,iend
        s(i,1,k)=-s(i,3,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=kbgn,kend
      DO i=ibgn,iend
        s(i,1,k)=s(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc == 3 .OR. sbc == 4 .OR. sbc == 5 .OR. sbc == 6)THEN
                                   ! Zero normal gradient condition or
! Radiation or user specified condition.
    DO k=kbgn,kend
      DO i=ibgn,iend
        s(i,1,k)=s(i,2,k)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BOUNDV', sbc
    CALL arpsstop ("arpstop called from boundv south bc",1)

  END IF

  5004  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE boundv

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE BOUNDW                   ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE boundw(s,nx,ny,nz,ibgn,iend,jbgn,jend)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for a quantity at the first and last w
!  points in the z-direction. Please note that the values at the corner
!  points may depend on the order that e-w, n-s and t-b boundary
!  conditions are applied.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  5/13/93 (D. Weber)
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Index to start the i direction
!    iend     Index to end the i direction
!    jbgn     Index to start the j direction
!    jend     Index to end the j direction
!
!    s        Input array defined from k=2 to k=nz-1
!
!  OUTPUT:
!
!    s       Output array at w point for k=1 to k=nz
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions

  INTEGER :: ibgn,iend         ! Bounds for computations in the x
                               ! direction

  INTEGER :: jbgn,jend         ! Bounds for computations in the y
                               ! direction

  REAL :: s(nx,ny,nz)          ! A scalar variable in the z direction
!
!-----------------------------------------------------------------------
!
!  The following integer parameters define the type of condition
!  at each boundary.
!
!   1 for rigid wall (mirror) type boundary condition.
!   2 for periodic boundary condition.
!   3 for zero normal gradient boundary condition.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
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
!  Set the top boundary condition
!
!-----------------------------------------------------------------------
!
  IF(tbc == 0) GO TO 5005

  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        s(i,j,nz)=-s(i,j,nz-2)
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        s(i,j,nz)=s(i,j,3)
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN  ! Zero normal gradient condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        s(i,j,nz)=s(i,j,nz-2)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BOUNDW', tbc
    CALL arpsstop ("arpstop called from boundw top bc",1)

  END IF

  5005  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set the bottom boundary condition
!
!-----------------------------------------------------------------------
!
  IF(bbc == 0) GO TO 5006

  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=jbgn,jend
      DO i=ibgn,iend
        s(i,j,1)=-s(i,j,3)
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        s(i,j,1)=s(i,j,nz-2)
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=jbgn,jend
      DO i=ibgn,iend
        s(i,j,1)=s(i,j,3)
      END DO
    END DO

  ELSE

    WRITE(6,900) 'BOUNDW', bbc
    CALL arpsstop ("arpstop called from boundw bottom bc",1)

  END IF

  5006  CONTINUE

  RETURN

  900   FORMAT(1X,'Invalid boundary condition option found in ',a,      &
               /1X,'The option was ',i3,' Job stopped.')
END SUBROUTINE boundw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE EXTNDSBC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE extndsbc(a,nx,ny,nz,vartyp,ebc,wbc,nbc,sbc,tbc,bbc)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extend the values of variable a at scalar points to the extra
!  fake zone on the boundary based on the boundary conditions. The
!  variable can be either a scalar or vector.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:Ming Xue
!  10/12/1996
!
!  MODIFICATIONS
!
!  12/17/1998 (Pengfei Zhang)
!  Changed this routine in compatible with vector variable at scalar
!  points such as fluxes on scalar grids.
!
!  2000/02/28 (Gene Bassett)
!  Added message passing markers and moved to bc3d.f
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Variable whose boundary values will be set here
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    vartyp   Type of the variable, 1 for vector and 0 for scalar
!
!  OUTPUT:
!
!    a        Boundary values of variable a.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  REAL :: a(0:nx,0:ny,0:nz)

  INTEGER :: vartyp

  INTEGER :: ebc   ! Parameter defining east   boundary condition type.
  INTEGER :: wbc   ! Parameter defining west   boundary condition type.
  INTEGER :: nbc   ! Parameter defining north  boundary condition type.
  INTEGER :: sbc   ! Parameter defining south  boundary condition type.
  INTEGER :: tbc   ! Parameter defining top    boundary condition type.
  INTEGER :: bbc   ! Parameter defining bottom boundary condition type.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k

  REAL, ALLOCATABLE :: mp_tem(:)  ! Temporary message passing array.

  INCLUDE 'mp.inc'
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (mp_opt > 0) THEN

    ALLOCATE (mp_tem( MAX(nx+1,ny+1)*(nz+1)*2 ), STAT = i)

    CALL acct_interrupt(mp_acct)
    CALL mpsendrecvextew(a,nx,ny,nz,ebc,wbc,mp_tem)
    CALL mpsendrecvextns(a,nx,ny,nz,nbc,sbc,mp_tem)

    DEALLOCATE(mp_tem)

  ENDIF


  IF (wbc == 0) THEN

!    do nothing

  ELSE IF (wbc == 1) THEN

    IF (vartyp == 0 ) THEN
      DO k=0,nz
        DO j=0,ny
          a(0,j,k) = a(3,j,k)
        END DO
      END DO
    ELSE
      DO k=0,nz
        DO j=0,ny
          a(0,j,k) = -a(3,j,k)
        END DO
      END DO
    END IF

  ELSE IF (wbc == 2) THEN

    IF (mp_opt == 0) THEN
      DO k=0,nz
        DO j=0,ny
          a(0,j,k)= a(nx-3,j,k)
        END DO
      END DO
    END IF

  ELSE ! including options 3 and 4

    DO k=0,nz
      DO j=0,ny
        a(0,j,k)=a(1,j,k)
      END DO
    END DO

  END IF

  IF (ebc == 0) THEN

!    do nothing

  ELSE IF (ebc == 1) THEN

    IF (vartyp == 0 ) THEN
      DO k=0,nz
        DO j=0,ny
          a(nx,j,k)=a(nx-3,j,k)
        END DO
      END DO
    ELSE
      DO k=0,nz
        DO j=0,ny
          a(nx,j,k) = -a(nx-3,j,k)
        END DO
      END DO
    END IF

  ELSE IF (ebc == 2 ) THEN

    IF (mp_opt == 0) THEN
      DO k=0,nz
        DO j=0,ny
          a(nx,j,k)=a(3,j,k)
        END DO
      END DO
    END IF

  ELSE ! including options 3 and 4

    DO k=0,nz
      DO j=0,ny
        a(nx,j,k)=a(nx-1,j,k)
      END DO
    END DO

  END IF
!
  IF ( sbc == 0) THEN

!    do nothing

  ELSE IF (sbc == 1) THEN

    IF (vartyp == 0 ) THEN
      DO k=0,nz
        DO i=0,nx
          a(i,0,k)=a(i,3,k)
        END DO
      END DO
    ELSE
      DO k=0,nz
        DO i=0,nx
          a(i,0,k) = -a(i,3,k)
        END DO
      END DO
    END IF

  ELSE IF (sbc == 2 ) THEN

    IF (mp_opt == 0) THEN
      DO k=0,nz
        DO i=0,nx
          a(i,0,k)=a(i,ny-3,k)
        END DO
      END DO
    END IF

  ELSE ! including options 3 and 4

    DO k=0,nz
      DO i=0,nx
        a(i,0,k)=a(i,1,k)
      END DO
    END DO

  END IF

  IF (nbc == 0) THEN

!    do nothing

  ELSE IF (nbc == 1) THEN

    IF (vartyp == 0 ) THEN
      DO k=0,nz
        DO i=0,nx
          a(i,ny,k)=a(i,ny-3,k)
        END DO
      END DO
    ELSE
      DO k=0,nz
        DO i=0,nx
          a(i,ny,k) = -a(i,ny-3,k)
        END DO
      END DO
    END IF

  ELSE IF (nbc == 2 ) THEN

    IF (mp_opt == 0) THEN
      DO k=0,nz
        DO i=0,nx
          a(i,ny,k)=a(i,3,k)
        END DO
      END DO
   END IF

  ELSE ! including options 3 and 4

    DO k=0,nz
      DO i=0,nx
        a(i,ny,k)=a(i,ny-1,k)
      END DO
    END DO

  END IF

  IF (bbc == 0) THEN

!    do nothing

  ELSE IF (bbc == 1) THEN

    IF (vartyp == 0 ) THEN
      DO i=0,nx
        DO j=0,ny
          a(i,j,0)=a(i,j,3)
        END DO
      END DO
    ELSE
      DO i=0,nx
        DO j=0,ny
          a(i,j,0) = -a(i,j,3)
        END DO
      END DO
    END IF

  ELSE IF (bbc == 2) THEN

    DO i=0,nx-1
      DO j=0,ny-1
        a(i,j,0)=a(i,j,nz-3)
      END DO
    END DO

  ELSE ! including option 3

    DO i=0,nx
      DO j=0,ny
        a(i,j,0)=a(i,j,1)
      END DO
    END DO

  END IF

  IF (tbc == 0) THEN

!    do nothing

  ELSE IF (tbc == 1) THEN

    IF (vartyp == 0 ) THEN
      DO i=0,nx
        DO j=0,ny
          a(i,j,nz)=a(i,j,nz-3)
        END DO
      END DO
    ELSE
      DO i=0,nx
        DO j=0,ny
          a(i,j,nz) = -a(i,j,nz-3)
        END DO
      END DO
    END IF

  ELSE IF (tbc == 2) THEN

    DO i=0,nx
      DO j=0,ny
        a(i,j,nz)=a(i,j,3)
      END DO
    END DO

  ELSE ! including options 3 and 4

    DO i=0,nx
      DO j=0,ny
        a(i,j,nz)=a(i,j,nz-1)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE extndsbc

