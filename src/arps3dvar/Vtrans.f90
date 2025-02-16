!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRANS                      ######
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
!  Convert (u,v,pprt,ptprt,qv,w) to x
!
!-----------------------------------------------------------------------
!
!  AUTHOR: JIDONG GAO
!  01/17/00
!
!
!-----------------------------------------------------------------------
!
!

SUBROUTINE trans(numctr,nx,ny,nz,u,v,pprt,ptprt,qv,w,x)

!-----------------------------------------------------------------------
!
!  INPUT:
!
!    numctr     The number of components of the control variables.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x-component of velocity at all time levels (m/s).
!    v        y-component of velocity at all time levels (m/s).
!    pprt     Perturbation pressure at all time levels (Pascal)
!    ptprt    Perturbation potential temperature at all time levels (K)
!    qv       Water vapor specific humidity at all time levels (kg/kg)
!    w        z-component of velocity at all time levels (m/s).
!
!    Output:
!    x        The control variable.
!
!    Temporary Variables:
!    itemp    Temporary variable indicating the position of the index.
!    i,j,k
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: numctr         ! The no. of components of the control var.
  INTEGER :: itemp, i,j,k ! Temporary variables.

  INTEGER :: nx,ny,nz          ! The number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)  ! Total u-velocity (m/s).
  REAL :: v     (nx,ny,nz)  ! Total v-velocity (m/s).
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure from that
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg).
  REAL :: w     (nx,ny,nz)  ! Total w-velocity (m/s).
  REAL :: x     (numctr)       ! Control variable.

  INCLUDE 'mp.inc'

  INTEGER :: ibgn, iend, jbgn, jend
  !INTEGER :: iuend,jvend
!
!    only initial conditions are controled.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ibgn = 1
  iend = nx-1
  jbgn = 1
  jend = ny-1

  !iuend = nx
  !jvend = ny

  IF (loc_x > 1) ibgn = 2
  IF (loc_y > 1) jbgn = 2

  IF (loc_x < nproc_x) THEN
    iend  = nx-2
    !iuend = nx-2
  END IF

  IF (loc_y < nproc_y) THEN
    jend  = ny-2
    !jvend = ny-2
  END IF

!-----------------------------------------------------------------------
!
!  We first place the initial conditions for (u,v,w,ptprt,pprt,qv)
!  in x in the order as they appear.
!
!-----------------------------------------------------------------------
!
  itemp = 0

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        x(itemp) = u(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        x(itemp) = v(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        x(itemp)=pprt (i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        x(itemp)=ptprt(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        x(itemp)=qv   (i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        x(itemp)=w    (i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  End of processing initial conditions.
!
!-----------------------------------------------------------------------
!
!
  RETURN
END SUBROUTINE trans
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADTRANS                    ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE adtrans(numctr,nx,ny,nz,u,v,pprt,ptprt,qv,w,x,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert x to (u,v,w,ptprt,pprt,qv)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: JIDONG GAO
!  01/17/00
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    numctr     The number of components of the control variables.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x-component of velocity at all time levels (m/s).
!    v        y-component of velocity at all time levels (m/s).
!    pprt     Perturbation pressure at all time levels (Pascal)
!    ptprt    Perturbation potential temperature at all time levels (K)
!    qv       Water vapor specific humidity at all time levels (kg/kg)
!    w        z-component of velocity at all time levels (m/s).
!
!    Output:
!    x        The control variable.
!
!    Temporary Variables:
!    itemp    Temporary variable indicating the position of the index.
!    i,j,k
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: numctr         ! The no. of components of the control var.
  INTEGER :: itemp, i,j,k ! Temporary variables.

  INTEGER :: nx,ny,nz          ! The number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)  ! Total u-velocity (m/s).
  REAL :: v     (nx,ny,nz)  ! Total v-velocity (m/s).
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure from that
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg).
  REAL :: w     (nx,ny,nz)  ! Total w-velocity (m/s).
  REAL :: x     (numctr)       ! Control variable.

  REAL :: tem1  (nx,ny,nz)  ! working array

  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'

  INTEGER :: ibgn, iend, jbgn, jend
  !INTEGER :: iuend,jvend

!
!    only initial conditions are controled.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ibgn = 1
  iend = nx-1
  jbgn = 1
  jend = ny-1

  !iuend = nx
  !jvend = ny

  IF (loc_x > 1)       ibgn = 2
  IF (loc_y > 1)       jbgn = 2

  IF (loc_x < nproc_x) THEN
    iend  = nx-2
    !iuend = nx-2
  END IF

  IF (loc_y < nproc_y) THEN
    jend  = ny-2
    !jvend = ny-2
  END IF

!-----------------------------------------------------------------------
!
!  We first place the initial conditions for (u,v,w,ptprt,pprt,qv)
!  in x in the order as they appear.
!
!-----------------------------------------------------------------------
!
  itemp = 0

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        u    (i,j,k) = x(itemp)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        v    (i,j,k) = x(itemp)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        pprt (i,j,k) = x(itemp)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        ptprt(i,j,k) = x(itemp)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        qv   (i,j,k) = x(itemp)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        itemp=itemp+1
        w    (i,j,k) = x(itemp)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  End of processing initial conditions.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(u, nx, ny, nz, ebc, wbc, 0, tem1)  ! 1- nx-1 is good
    CALL mpsendrecv2dns(u, nx, ny, nz, nbc, sbc, 0, tem1)

    CALL mpsendrecv2dew(v, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(v, nx, ny, nz, nbc, sbc, 0, tem1)  ! 1 - ny-1 is good

    CALL mpsendrecv2dew(w, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(w, nx, ny, nz, nbc, sbc, 0, tem1)

    CALL mpsendrecv2dew(pprt, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(pprt, nx, ny, nz, nbc, sbc, 0, tem1)

    CALL mpsendrecv2dew(ptprt, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(ptprt, nx, ny, nz, nbc, sbc, 0, tem1)

    CALL mpsendrecv2dew(qv, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(qv, nx, ny, nz, nbc, sbc, 0, tem1)

    !CALL mpsendrecv2dew(u, nx, ny, nz, ebc, wbc, 1, tem1)  ! make sure nx contains good value
    !CALL mpsendrecv2dns(v, nx, ny, nz, nbc, sbc, 2, tem1)  ! make sure ny contains good value
!    CALL acct_stop_inter
  END IF

  RETURN
END SUBROUTINE adtrans
