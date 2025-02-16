!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECV2DEW             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecv2dew(var,nx,ny,nz,ebc,wbc,stagdim,tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive east/west boundary data between processors to
!  update the fake zones.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/19/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ebc      East boundary condition
!    wbc      West boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!              =3, staggered in the z-direction (e.g. w).
!
!  INPUT & OUTPUT
!
!    var      Variable for which boundaries need updating.
!
!  WORK array
!
!    tem      Work array (with a size at least nyXnzX2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in
                                       ! x, y and z directions
  INTEGER, INTENT(IN) :: ebc,wbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
                               !  =3, staggered in the z-direction (e.g. w).
  REAL(P), INTENT(INOUT) :: var(nx,ny,nz)
  REAL(P), INTENT(INOUT) :: tem(ny,nz,2)   ! Work array.


!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: j, k
  INTEGER :: si,sj,sk

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0
  sj = 0
  sk = 0
  IF (stagdim == 1) si = 1
  IF (stagdim == 2) sj = 1
  IF (stagdim == 3) sk = 1

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    IF(ebc == 2) THEN             ! periodic boundary
      dest = proc(1+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN
      source = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of
  ! the east neighbor
  !
  DO k=1,nz           ! -1+sk  send full rank
    DO j=1,ny         ! -1+sj
      tem(j,k,1) = var(nx-2,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),ny*nz,mpi_p,dest,  mptag+tag_w,       &
                    tem(:,:,2),ny*nz,mpi_p,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 .OR. wbc == 2 ) THEN ! .NOT. (loc_x ==1 .AND. wbc /=2))
    DO k=1,nz         ! -1+sk
      DO j=1,ny       ! -1+sj
        var(1,j,k) = tem(j,k,2)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN             ! periodic boundary
      dest = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    IF(ebc == 2) THEN
      source = proc(1+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of
  ! the west neighbor
  !
  DO k=1,nz             ! -1+sk
    DO j=1,ny           ! -1+sj
      tem(j,k,1) = var(2+si,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),ny*nz,MPI_P,dest,  mptag+tag_e,       &
                    tem(:,:,2),ny*nz,MPI_P,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x .OR. ebc == 2 )  THEN
                              !.NOT. (loc_x == nproc_x .AND. ebc /=2))
    DO k=1,nz         ! -1+sk
      DO j=1,ny       ! -1+sj
        var(nx-1+si,j,k) = tem(j,k,2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecv2dew
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECV2DNS             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecv2dns(var,nx,ny,nz,nbc,sbc,stagdim,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive north/south boundary data between processors to
!  update the fake zones.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/19/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    nbc      North boundary condition
!    sbc      South boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!              =3, staggered in the z-direction (e.g. w).
!
!  INPUT & OUTPUT:
!
!    var      Variable for which boundaries need updating.
!
!    tem      Work array (with a size at least nx X nz X 2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in
                                       ! x, y and z directions
  INTEGER, INTENT(IN) :: nbc,sbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
                               !  =3, staggered in the z-direction (e.g. w).
  REAL(P), INTENT(INOUT) :: var(nx,ny,nz)
  REAL(P), INTENT(INOUT) :: tem(nx,nz,2)   ! Work array.

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i,k
  INTEGER :: si,sj,sk

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0
  sj = 0
  sk = 0
  IF (stagdim == 1) si = 1
  IF (stagdim == 2) sj = 1
  IF (stagdim == 3) sk = 1

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    IF(sbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    IF(nbc == 2) THEN
      source = proc(loc_x)
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of
  ! the south neighbor
  !
  DO k=1,nz         ! -1+sk
    DO i=1,nx       ! -1+si
      tem(i,k,1) = var(i,2+sj,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),nx*nz,MPI_P,dest,  mptag+tag_n,       &
                    tem(:,:,2),nx*nz,MPI_P,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y .OR. nbc == 2 )  THEN
                               ! .NOT. (loc_y == nproc_y .AND. nbc /=2))
    DO k=1,nz      ! -1+sk
      DO i=1,nx    ! -1+si
        var(i,ny-1+sj,k) = tem(i,k,2)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    IF(nbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x)
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    IF(sbc == 2) THEN
      source = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of
  ! the north neighbor
  !
  DO k=1,nz     ! -1+sk
    DO i=1,nx   ! -1+si
      tem(i,k,1) = var(i,ny-2,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),nx*nz,MPI_P,dest,  mptag+tag_s,       &
                    tem(:,:,2),nx*nz,MPI_P,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 .OR. sbc == 2 )  THEN
                                    ! .NOT. (loc_y == 1 .AND. sbc /=2))
    DO k=1,nz     ! -1+sk
      DO i=1,nx   ! -1+si
        var(i,1,k) = tem(i,k,2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecv2dns
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECV1DEW             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecv1dew(var,nx,ny,ebc,wbc,stagdim,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive east/west boundary 1D data between processors to
!  update the fake zones.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/24/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    ebc      East boundary condition
!    wbc      West boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!
!  INPUT & OUTPUT
!
!    var      Variable for which boundaries need updating.
!
!  WORK array
!
!    tem      Work array (with a size at least nyX2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny         ! Number of grid points in
                                       ! x and y directions
  INTEGER, INTENT(IN) :: ebc,wbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
  REAL(P), INTENT(INOUT) :: var(nx,ny)
  REAL(P), INTENT(INOUT) :: tem(ny,2)   ! Work array.


!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: j
  INTEGER :: si,sj

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0
  sj = 0
  IF (stagdim == 1) si = 1
  IF (stagdim == 2) sj = 1

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    IF(ebc == 2) THEN             ! periodic boundary
      dest = proc(1+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN
      source = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of
  ! the east neighbor
  !
  DO j=1,ny-1+sj
    tem(j,1) = var(nx-2,j)
  END DO

  CALL mpi_sendrecv(tem(:,1),ny,MPI_P,dest,  mptag+tag_w,       &
                    tem(:,2),ny,MPI_P,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 .OR. wbc == 2 ) THEN ! .NOT. (loc_x ==1 .AND. wbc /=2))
    DO j=1,ny-1+sj
      var(1,j) = tem(j,2)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN             ! periodic boundary
      dest = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    IF(ebc == 2) THEN
      source = proc(1+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of
  ! the west neighbor
  !
  DO j=1,ny-1+sj
    tem(j,1) = var(2+si,j)
  END DO

  CALL mpi_sendrecv(tem(:,1),ny,MPI_P,dest,  mptag+tag_e,       &
                    tem(:,2),ny,MPI_P,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x .OR. ebc == 2 )  THEN
                              !.NOT. (loc_x == nproc_x .AND. ebc /=2))
    DO j=1,ny-1+sj
      var(nx-1+si,j) = tem(j,2)
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecv1dew
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECV1DNS             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecv1dns(var,nx,ny,nbc,sbc,stagdim,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive north/south boundary data between processors to
!  update the fake zones.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/24/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    nbc      North boundary condition
!    sbc      South boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!              =3, staggered in the z-direction (e.g. w).
!
!  INPUT & OUTPUT:
!
!    var      Variable for which boundaries need updating.
!
!    tem      Work array (with a size at least nx X 2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny         ! Number of grid points in
                                       ! x, and y directions
  INTEGER, INTENT(IN) :: nbc,sbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
  REAL(P), INTENT(INOUT) :: var(nx,ny)
  REAL(P), INTENT(INOUT) :: tem(nx,2)   ! Work array.

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i
  INTEGER :: si,sj

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0
  sj = 0
  IF (stagdim == 1) si = 1
  IF (stagdim == 2) sj = 1

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    IF(sbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    IF(nbc == 2) THEN
      source = proc(loc_x)
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of
  ! the south neighbor
  !
  DO i=1,nx-1+si
    tem(i,1) = var(i,2+sj)
  END DO

  CALL mpi_sendrecv(tem(:,1),nx,MPI_P,dest,  mptag+tag_n,       &
                    tem(:,2),nx,MPI_P,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y .OR. nbc == 2 )  THEN
                               ! .NOT. (loc_y == nproc_y .AND. nbc /=2))
    DO i=1,nx-1+si
      var(i,ny-1+sj) = tem(i,2)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    IF(nbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x)
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    IF(sbc == 2) THEN
      source = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of
  ! the north neighbor
  !
  DO i=1,nx-1+si
    tem(i,1) = var(i,ny-2)
  END DO

  CALL mpi_sendrecv(tem(:,1),nx,MPI_P,dest,  mptag+tag_s,       &
                    tem(:,2),nx,MPI_P,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 .OR. sbc == 2 )  THEN
                                    ! .NOT. (loc_y == 1 .AND. sbc /=2))
    DO i=1,nx-1+si
      var(i,1) = tem(i,2)
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecv1dns
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECV1DIEW            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecv1diew(var,nx,ny,ebc,wbc,stagdim,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive east/west boundary 1D data between processors to
!  update the fake zones.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/24/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    ebc      East boundary condition
!    wbc      West boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!
!  INPUT & OUTPUT
!
!    var      Variable for which boundaries need updating.
!
!  WORK array
!
!    tem      Work array (with a size at least nyX2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny         ! Number of grid points in
                                       ! x and y directions
  INTEGER, INTENT(IN) :: ebc,wbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
  INTEGER, INTENT(INOUT) :: var(nx,ny)
  INTEGER, INTENT(INOUT) :: tem(ny,2)   ! Work array.


!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: j
  INTEGER :: si,sj

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0
  sj = 0
  IF (stagdim == 1) si = 1
  IF (stagdim == 2) sj = 1

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    IF(ebc == 2) THEN             ! periodic boundary
      dest = proc(1+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN
      source = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of
  ! the east neighbor
  !
  DO j=1,ny-1+sj
    tem(j,1) = var(nx-2,j)
  END DO

  CALL mpi_sendrecv(tem(:,1),ny,MPI_INTEGER,dest,  mptag+tag_w,       &
                    tem(:,2),ny,MPI_INTEGER,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 .OR. wbc == 2 ) THEN ! .NOT. (loc_x ==1 .AND. wbc /=2))
    DO j=1,ny-1+sj
      var(1,j) = tem(j,2)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN             ! periodic boundary
      dest = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    IF(ebc == 2) THEN
      source = proc(1+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of
  ! the west neighbor
  !
  DO j=1,ny-1+sj
    tem(j,1) = var(2+si,j)
  END DO

  CALL mpi_sendrecv(tem(:,1),ny,MPI_INTEGER,dest,  mptag+tag_e,       &
                    tem(:,2),ny,MPI_INTEGER,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x .OR. ebc == 2 )  THEN
                              !.NOT. (loc_x == nproc_x .AND. ebc /=2))
    DO j=1,ny-1+sj
      var(nx-1+si,j) = tem(j,2)
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecv1diew
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECV1DINS            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecv1dins(var,nx,ny,nbc,sbc,stagdim,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive north/south boundary data between processors to
!  update the fake zones.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/24/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    nbc      North boundary condition
!    sbc      South boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!              =3, staggered in the z-direction (e.g. w).
!
!  INPUT & OUTPUT:
!
!    var      Variable for which boundaries need updating.
!
!    tem      Work array (with a size at least nx X 2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny         ! Number of grid points in
                                       ! x, and y directions
  INTEGER, INTENT(IN) :: nbc,sbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
  INTEGER, INTENT(INOUT) :: var(nx,ny)
  INTEGER, INTENT(INOUT) :: tem(nx,2)   ! Work array.

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i
  INTEGER :: si,sj

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0
  sj = 0
  IF (stagdim == 1) si = 1
  IF (stagdim == 2) sj = 1

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    IF(sbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    IF(nbc == 2) THEN
      source = proc(loc_x)
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of
  ! the south neighbor
  !
  DO i=1,nx-1+si
    tem(i,1) = var(i,2+sj)
  END DO

  CALL mpi_sendrecv(tem(:,1),nx,MPI_INTEGER,dest,  mptag+tag_n,       &
                    tem(:,2),nx,MPI_INTEGER,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y .OR. nbc == 2 )  THEN
                               ! .NOT. (loc_y == nproc_y .AND. nbc /=2))
    DO i=1,nx-1+si
      var(i,ny-1+sj) = tem(i,2)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    IF(nbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x)
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    IF(sbc == 2) THEN
      source = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of
  ! the north neighbor
  !
  DO i=1,nx-1+si
    tem(i,1) = var(i,ny-2)
  END DO

  CALL mpi_sendrecv(tem(:,1),nx,MPI_INTEGER,dest,  mptag+tag_s,       &
                    tem(:,2),nx,MPI_INTEGER,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 .OR. sbc == 2 )  THEN
                                    ! .NOT. (loc_y == 1 .AND. sbc /=2))
    DO i=1,nx-1+si
      var(i,1) = tem(i,2)
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecv1dins
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSENDRECVEXTEW            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecvextew(var,nx,ny,nz,ebc,wbc,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive east/west boundary data between processors to
!  update the fake zones. This is for the extended array whcih
!  has two instead of one fake zones on each boundary (arrays run
!  from 0:nx,0:ny,0:nz).
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/24/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ebc      East boundary condition
!    wbc      West boundary condition
!
!  INPUT & OUTPUT
!
!    var      Variable for which boundaries need updating.
!
!  WORK array
!
!    tem      Work array (with a size at least (ny+1)X(nz+1)X2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in
                                       ! x, y and z directions
  INTEGER, INTENT(IN) :: ebc,wbc

  REAL(P), INTENT(INOUT) :: var(0:nx,0:ny,0:nz)
  REAL(P), INTENT(INOUT) :: tem(0:ny,0:nz,2)   ! Work array.


!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: j, k
  INTEGER :: si,sj,sk

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary
  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    IF(ebc == 2) THEN             ! periodic boundary
      dest = proc(1+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN
      source = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of
  ! the east neighbor
  !
  DO k=0,nz
    DO j=0,ny
      tem(j,k,1) = var(nx-3,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),(ny+1)*(nz+1),MPI_P,dest,  mptag+tag_w, &
                    tem(:,:,2),(ny+1)*(nz+1),MPI_P,source,mptag+tag_w, &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 .OR. wbc == 2 ) THEN ! .NOT. (loc_x ==1 .AND. wbc /=2))
    DO k=0,nz
      DO j=0,ny
        var(0,j,k) = tem(j,k,2)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN             ! periodic boundary
      dest = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    IF(ebc == 2) THEN
      source = proc(1+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of
  ! the west neighbor
  !
  DO k=0,nz
    DO j=0,ny
      tem(j,k,1) = var(3,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),(ny+1)*(nz+1),MPI_P,dest,  mptag+tag_e, &
                    tem(:,:,2),(ny+1)*(nz+1),MPI_P,source,mptag+tag_e, &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x .OR. ebc == 2 )  THEN
                              !.NOT. (loc_x == nproc_x .AND. ebc /=2))
    DO k=0,nz
      DO j=0,ny
        var(nx,j,k) = tem(j,k,2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecvextew
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MPSENDRECVEXTNS             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsendrecvextns(var,nx,ny,nz,nbc,sbc,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive north/south boundary data between processors to
!  update the fake zones. This is for the extended array whcih
!  has two instead of one fake zones on each boundary (arrays run
!  from 0:nx,0:ny,0:nz).
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  11/24/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    nbc      North boundary condition
!    sbc      South boundary condition
!
!  INPUT & OUTPUT:
!
!    var      Variable for which boundaries need updating.
!
!    tem      Work array (with a size at least (nx+1) X (nz+1) X 2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in
                                       ! x, y and z directions
  INTEGER, INTENT(IN) :: nbc,sbc

  REAL(P), INTENT(INOUT) :: var(0:nx,0:ny,0:nz)
  REAL(P), INTENT(INOUT) :: tem(0:nx,0:nz,2)   ! Work array.

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i,k

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    IF(sbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    IF(nbc == 2) THEN
      source = proc(loc_x)
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of
  ! the south neighbor
  !
  DO k=0,nz
    DO i=0,nx
      tem(i,k,1) = var(i,3,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),(nx+1)*(nz+1),MPI_P,dest,  mptag+tag_n,&
                    tem(:,:,2),(nx+1)*(nz+1),MPI_P,source,mptag+tag_n,&
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y .OR. nbc == 2 )  THEN
                               ! .NOT. (loc_y == nproc_y .AND. nbc /=2))
    DO k=0,nz
      DO i=0,nx
        var(i,ny,k) = tem(i,k,2)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    IF(nbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x)
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    IF(sbc == 2) THEN
      source = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of
  ! the north neighbor
  !
  DO k=0,nz
    DO i=0,nx
      tem(i,k,1) = var(i,ny-3,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),(nx+1)*(nz+1),MPI_P,dest,  mptag+tag_s,&
                    tem(:,:,2),(nx+1)*(nz+1),MPI_P,source,mptag+tag_s,&
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 .OR. sbc == 2 )  THEN
                                    ! .NOT. (loc_y == 1 .AND. sbc /=2))
    DO k=0,nz
      DO i=0,nx
        var(i,0,k) = tem(i,k,2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE mpsendrecvextns
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE MPFILLEXTARRAY3D               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpfillextarray3d(var,nx,ny,nz,iext,jext,outvar,              &
                            ebc,wbc,sbc,nbc,stagdim,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fills a extend array from the given local array. The given array
!  contains valid data in section (1:nx,1:ny,nz) and will be
!  extended by iext/jext points in all four horizontal boundaries.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  04/12/2007
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable contains valid data
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    iext     # of grid points to be extended in both west and east boundary
!    jext     # of grid points to be extended in both south and north boundary
!
!    ebc      East boundary condition
!    wbc      West boundary condition
!    sbc      South boundary condition
!    nbc      North boundary condition
!
!    stagdim  Dimension of grid staggering:
!              =0, no staggering;
!              =1, staggered in the x-direction (e.g. u);
!              =2, staggered in the y-direction (e.g. v);
!
!  INPUT & OUTPUT
!
!    outvar      Variable with extended dimensions to be returned
!    istatus     Return status (0 - success, other - fail)
!
!  WORK array
!
!    tem1      Work array in X direction(with a size at least iext*ny*nz).
!    tem2      Work array in Y direction(with a size at least (nx+2*iext)*jext*nz).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in
                                       ! x, y and z directions

  INTEGER, INTENT(IN) :: iext, jext
  INTEGER, INTENT(IN) :: ebc,wbc, sbc,nbc
  INTEGER, INTENT(IN) :: stagdim       ! Dimension of grid staggering:
                               !  =0, no staggering;
                               !  =1, staggered in the x-direction (e.g. u);
                               !  =2, staggered in the y-direction (e.g. v);
                               !  =3, staggered in the z-direction (e.g. w).

  REAL(P), INTENT(IN)  :: var(nx,ny,nz)
  REAL(P), INTENT(OUT) :: outvar(1-iext:nx+iext,1-jext:ny+jext,nz)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i, j, k
  INTEGER :: si,sj

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary
  INTEGER :: mpi_p

  INTEGER :: ii, jj
  INTEGER :: isize, jsize, mwesize, msnsize

  REAL(P), ALLOCATABLE :: temx1(:, :, :)   ! Work array.
  REAL(P), ALLOCATABLE :: temx2(:, :, :)
  REAL(P), ALLOCATABLE :: temy1(:, :, :)
  REAL(P), ALLOCATABLE :: temy2(:, :, :)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  si = 0                     ! valid range - 1:nx-1, 1:ny-1
  sj = 0
  IF (stagdim == 1) si = 1   ! valid range - 1:nx
  IF (stagdim == 2) sj = 1   ! valid range - 1:ny

!-----------------------------------------------------------------------
!
! Fill internal sections
!
!-----------------------------------------------------------------------

  outvar(:,:,:) = 0.0     ! initialize to zeros

  DO k = 1,nz
    DO j = 1,ny-1+sj
      DO i = 1,nx-1+si
        outvar(i,j,k) = var(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Message passing parts
!
!-----------------------------------------------------------------------

  isize = iext+1
  jsize = jext+1

  IF ( (isize > nx-2) .OR. (jsize > ny-2) ) THEN
    WRITE(6,'(1x,a,/,2(8x,a,2(I4,a),/),8x,a)')                          &
      'ERROR: the message zone is too large in mpfillextarray3d.',      &
      'The local path size is nx = ',nx,', ny = ',ny,'.',               &
      'The required message zone is isize = ',isize,', jsize = ',jsize,'.', &
      'Please note that the current framework can only handle message exchanges between neighbors.'
    CALL arpsstop('ERROR too large message zone in mpfillextarray3d.',1)
  END IF

  ALLOCATE(temx1(isize,ny,nz), STAT = istatus)
  ALLOCATE(temx2(isize,ny,nz), STAT = istatus)

  ALLOCATE(temy1(-iext+1:nx+iext,jsize,nz), STAT = istatus)
  ALLOCATE(temy2(-iext+1:nx+iext,jsize,nz), STAT = istatus)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  mwesize = isize*ny*nz
  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    IF(ebc == 2) THEN             ! periodic boundary, use ebc since wbc is
                                  ! usually set to 0 in MPI mode
      dest = proc(1+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN
      source = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of
  ! the east neighbor
  !
  DO k=1,nz           ! -1+sk  send full rank
    DO j=1,ny         ! -1+sj
      DO i = 1,isize  ! = iext+1   !  1    2    3    4       For example, iext = 3
        ii = nx- (isize+2) + i     ! nx-5 nx-4 nx-3 nx-2     Send buffer
                                   ! -2    -1   0    1       Receive buffer
        temx1(i,j,k) = var(ii,j,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(temx1,mwesize,MPI_P,dest,  mptag+tag_w,          &
                    temx2,mwesize,MPI_P,source,mptag+tag_w,          &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 .OR. wbc == 2 ) THEN ! .NOT. (loc_x ==1 .AND. wbc /=2))
    DO k=1,nz         ! -1+sk
      DO j=1,ny       ! -1+sj
        DO i = 1,isize
          ii = -1*iext + i
          outvar(ii,j,k) = temx2(i,j,k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    IF(wbc == 2) THEN             ! periodic boundary
      dest = proc(nproc_x+nproc_x*(loc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    IF(ebc == 2) THEN
      source = proc(1+nproc_x*(loc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of
  ! the west neighbor
  !
  DO k=1,nz
    DO j=1,ny
      DO i = 1,isize                    !  1   2   3    4     For example, iext = 3
        temx1(i,j,k) = var(1+si+i,j,k)  !  2   3   4    5     Send buffer for no-stag
                                        ! nx-1 nx nx+1 nx+2
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(temx1,mwesize,MPI_P,dest,  mptag+tag_e,          &
                    temx2,mwesize,MPI_P,source,mptag+tag_e,          &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x .OR. ebc == 2 )  THEN
                              !.NOT. (loc_x == nproc_x .AND. ebc /=2))
    DO k=1,nz
      DO j=1,ny
        DO i = 1,isize
          outvar(nx-2+si+i,j,k) = temx2(i,j,k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  msnsize = (nx+2*iext)*jsize*nz
  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    IF(sbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    IF(nbc == 2) THEN
      source = proc(loc_x)
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of
  ! the south neighbor
  !
  DO k=1,nz
    DO j = 1,jsize
      DO i=1-iext,nx+iext
        temy1(i,j,k) = outvar(i,j+1+sj,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(temy1,msnsize,MPI_P,dest,  mptag+tag_n,          &
                    temy2,msnsize,MPI_P,source,mptag+tag_n,          &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y .OR. nbc == 2 )  THEN
                               ! .NOT. (loc_y == nproc_y .AND. nbc /=2))
    DO k=1,nz
      DO j = 1,jsize
        DO i=1-iext,nx+iext
          outvar(i,ny-2+j+sj,k) = temy2(i,j,k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    IF(nbc == 2) THEN             ! periodic boundary
      dest = proc(loc_x)
    ELSE
      dest = MPI_PROC_NULL
    END IF
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    IF(sbc == 2) THEN
      source = proc(loc_x+nproc_x*(nproc_y-1))
    ELSE
      source = MPI_PROC_NULL
    END IF
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of
  ! the north neighbor
  !
  DO k=1,nz
    DO j = 1,jsize
      DO i=1-iext,nx+iext
        jj = ny-(jsize+2) + j
        temy1(i,j,k) = outvar(i,jj,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(temy1,msnsize,MPI_P,dest,  mptag+tag_s,          &
                    temy2,msnsize,MPI_P,source,mptag+tag_s,          &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 .OR. sbc == 2 )  THEN
                                    ! .NOT. (loc_y == 1 .AND. sbc /=2))
    DO k=1,nz
      DO j = 1,jsize
        DO i=1-iext,nx+iext
          jj = -1*jext+j
          outvar(i,jj,k) = temy2(i,j,k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Just before return
!
!-----------------------------------------------------------------------

  DEALLOCATE(temx1,temx2)
  DEALLOCATE(temy1,temy2)

  RETURN
END SUBROUTINE mpfillextarray3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPUPDATER                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpupdater(var,num)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast the value of var from process 0 to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Gene Bassett
!  2000/04/24
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    var      Array to update (INPUT on proc 0, OUTPUT for rest).
!    num      Number of elements in the array.
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER :: num
  REAL(P) :: var(num)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat
  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_bcast(var,num,mpi_p,0,mpi_comm_world,imstat)
  IF (imstat /= 0) THEN
    WRITE (6,'(1x,a,I3.3)') 'MPUPDATER: error on processor - ',myproc
  END IF

  RETURN
END SUBROUTINE mpupdater
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPUPDATEI                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpupdatei(var,num)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast the value of var from process 0 to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Gene Bassett
!  2000/04/24
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    var      Variable to update (INPUT on proc 0, OUTPUT for rest).
!    num      Number of elements in the array.
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: num
  INTEGER :: var(num)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpi_bcast(var,num,mpi_integer,0,mpi_comm_world,imstat)
  IF (imstat /= 0) THEN
    WRITE (6,*) "MPUPDATEI: error on processor",myproc
  END IF

  RETURN
END SUBROUTINE mpupdatei
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPUPDATEL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpupdatel(var,num)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast the value of var from process 0 to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    var      Variable to update (INPUT on proc 0, OUTPUT for rest).
!    num      Number of elements in the array.
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: num
  LOGICAL :: var(num)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpi_bcast(var,num,MPI_LOGICAL,0,mpi_comm_world,imstat)
  IF (imstat /= 0) THEN
    WRITE (6,*) 'MPUPDATEL: error on processor ',myproc
  END IF

  RETURN
END SUBROUTINE mpupdatel
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPUPDATEC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpupdatec(str,lenstr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast the string str from process 0 to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Gene Bassett
!  2000/04/24
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    str      String to update (INPUT on proc 0, OUTPUT for rest).
!    lenstr   Length of str.
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: lenstr
  CHARACTER (LEN=lenstr) :: str

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpi_bcast(str,lenstr,mpi_character,0,mpi_comm_world,imstat)
  IF (imstat /= 0) THEN
    WRITE (6,*) "MPUPDATEC: error on processor",myproc
  END IF

  RETURN
END SUBROUTINE mpupdatec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Miscellaneous MPI subroutines (not in ARPS standard format)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mpexit(errcode)

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INTEGER :: errcode
  INTEGER :: imstat

  IF (errcode == 0) THEN
    CALL mpi_finalize (imstat)
  ELSE
    CALL mpi_abort (mpi_comm_world, errcode, imstat)
  ENDIF

  RETURN
END SUBROUTINE mpexit

!#######################################################################

SUBROUTINE inctag

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  !
  ! MPI standard only requires MPI_TAG_UB be no less than 32767.
  !
  IF (gentag < 100 .OR. gentag > 32700) gentag = 100
  gentag = gentag + 100

  RETURN
END SUBROUTINE inctag

!#######################################################################

SUBROUTINE mpbarrier

  INCLUDE 'mpif.h'
  INTEGER :: imstat

  CALL mpi_barrier (mpi_comm_world, imstat)

  RETURN
END SUBROUTINE mpbarrier

!#######################################################################

SUBROUTINE mptotal(var)

  USE arps_precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(P) :: var, vartm
  INTEGER :: i,j,imstat

  INTEGER :: mpi_p

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_allreduce(var,vartm,1,mpi_p,MPI_SUM,mpi_comm_world,imstat)

  var = vartm

  RETURN
END SUBROUTINE mptotal

SUBROUTINE mptotali(var)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: var, vartm
  INTEGER :: imstat

  CALL mpi_allreduce(var, vartm, 1, MPI_INTEGER, MPI_SUM,               &
                     mpi_comm_world, imstat)

  var = vartm

  RETURN
END SUBROUTINE mptotali

SUBROUTINE mpmax0(amax,amin)

  USE arps_precision

!
!  Modified by Dan Weber, May 4, 1998
!  Replaces code above for use on t3d/t3e system.
!  mpi_allreduce is not working properly...
!
  IMPLICIT NONE

  REAL(P), INTENT(INOUT) :: amax,amin

  INCLUDE 'mpif.h'

  REAL(P) :: maxtm, mintm
  INTEGER :: imstat

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    start of executable code....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL mpi_allreduce (amax, maxtm, 1, MPI_REAL, MPI_MAX,
!    :     MPI_COMM_WORLD, imstat)  ! commented out because the T3E
                                 ! has trouble with mpi_allreduce
  CALL mpi_reduce(amax,maxtm,1,mpi_p,MPI_MAX,0, mpi_comm_world,imstat)
  CALL mpi_bcast(maxtm,1,mpi_p,0,mpi_comm_world,imstat)

  amax = maxtm

!  CALL mpi_allreduce (amin, mintm, 1, MPI_REAL, MPI_MIN,
!    :     MPI_COMM_WORLD, imstat)  ! commented out because the T3E
                                 ! has trouble with mpi_allreduce
  CALL mpi_reduce(amin,mintm,1,mpi_p,MPI_MIN,0,mpi_comm_world,imstat)
  CALL mpi_bcast(mintm,1,mpi_p,0,mpi_comm_world,imstat)

  amin = mintm

  RETURN
END SUBROUTINE mpmax0

!#######################################################################

SUBROUTINE mpmax(amax,amin,nx,ny,nz,imax,jmax,kmax,imin,jmin,kmin)

!
!  Modified by Dan Weber, October 23, 1997
!

  USE arps_precision

  IMPLICIT NONE

  INTEGER :: nx,ny,nz,imax,jmax,kmax,imin,jmin,kmin,itema,itemb
  REAL(P) :: amax,amin
  INTEGER :: imstat

  INTEGER :: mpi_2p

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  INTEGER :: mpi_status(mpi_status_size)
  REAL(P) :: maxarr (2), minarr(2)
  REAL(P) :: maxtm  (2), mintm(2)
  INTEGER :: maxpack(3), maxunpack(3)
  INTEGER :: minpack(3), minunpack(3)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    start of executable code....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_2p = MPI_2REAL
  IF (P == DP) mpi_2p = MPI_2DOUBLE_PRECISION

  CALL inctag
  maxtm(1)  = 0.0
  maxtm(2)  = 0.0
  maxarr(1) = amax
  maxarr(2) = REAL(myproc)

!  CALL mpi_allreduce (maxarr, maxtm, 1, MPI_2REAL, MPI_MAXLOC,
!    :     MPI_COMM_WORLD, imstat)  ! commented out because the T3E
                                   ! has trouble with mpi_allreduce
  CALL mpi_reduce(maxarr,maxtm,1,MPI_2p,MPI_MAXLOC,0,                &
                  MPI_COMM_WORLD,imstat)
  CALL mpi_bcast (maxtm,1,MPI_2p,0,MPI_COMM_WORLD,imstat)

  itema = nint(maxtm(2))
  IF(myproc == itema .AND. itema /= 0)THEN    ! send only if
    maxpack(1) = imax + (nx-3)*(loc_x-1)      ! itema .ne. myproc=0!!!
    maxpack(2) = jmax + (ny-3)*(loc_y-1)
    maxpack(3) = kmax
    CALL mpi_send (maxpack,3,MPI_INTEGER,0,                             &
                   gentag,MPI_COMM_WORLD,imstat)

    !wdt forced buffering
    !CALL mpi_bsend (maxpack,3,mpi_integer,0,                           &
    !               gentag,mpi_comm_world,imstat)

  END IF

  IF(myproc == 0 .AND. myproc /= itema)THEN ! receive only if
                                            ! itema .ne. myproc=0
    CALL mpi_recv (maxunpack,3,mpi_integer,itema,                       &
                   gentag,mpi_comm_world,mpi_status,imstat)
    imax = maxunpack(1)
    jmax = maxunpack(2)
    kmax = maxunpack(3)
    amax = maxtm(1)
  END IF


  mintm(1) = 0.0
  mintm(2) = 0.0
  minarr(1) = amin
  minarr(2) = REAL(myproc)

!  CALL mpi_allreduce (minarr, mintm, 1, MPI_2REAL, MPI_MINLOC,
!    :     MPI_COMM_WORLD, imstat)  ! commented out because the T3E
                                    ! has trouble with mpi_allreduce

  CALL mpi_reduce(minarr,mintm,1,MPI_2p,MPI_MINLOC,0,                &
                  MPI_COMM_WORLD,imstat)
  CALL mpi_bcast (mintm,1,MPI_2p,0,MPI_COMM_WORLD,imstat)

  itemb = nint(mintm(2))
  IF (myproc == itemb .AND. itemb /= 0) THEN    ! send only if
    minpack(1) = imin + (nx-3)*(loc_x-1)        ! itema .ne. myproc=0!!!
    minpack(2) = jmin + (ny-3)*(loc_y-1)
    minpack(3) = kmin
    CALL mpi_send (minpack,3,mpi_integer,0,                             &
                     gentag+1,mpi_comm_world,imstat)

    !wdt forced buffering
    !CALL mpi_bsend (minpack,3,mpi_integer,0,                           &
    !                 gentag+1,mpi_comm_world,imstat)
  END IF

  IF (myproc == 0 .AND. myproc /= itemb) THEN ! receive only if
                                              ! itemb .ne. myproc=0
    CALL mpi_recv (minunpack,3,MPI_INTEGER,itemb,                       &
                   gentag+1,MPI_COMM_WORLD,mpi_status,imstat)

    imin = minunpack(1)
    jmin = minunpack(2)
    kmin = minunpack(3)
    amin = mintm(1)
  END IF

  RETURN
END SUBROUTINE mpmax

SUBROUTINE mpinit_proc(ioption)

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: ioption  ! 0 - Initialize from beginning
                                  ! 1 - Do not call mpi_init
  INTEGER :: imstat

  mp_opt = 1

  IF (ioption == 0) CALL mpi_init( imstat )

  CALL mpi_comm_rank( mpi_comm_world, myproc, imstat )

  RETURN
END SUBROUTINE mpinit_proc

SUBROUTINE mpinit_var

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'

  INTEGER :: i,j,k,l,numg,parent
  INTEGER :: mytid,nprocs0
  INTEGER :: imstat

  CALL mpi_comm_size( mpi_comm_world, nprocs0, imstat )

  nprocs = nproc_x * nproc_y

  IF(nprocs > max_proc) THEN

    WRITE (6,*) "ERROR: number of processors exceeds maximum ",     &
                "specified in mp.inc:"
    WRITE (6,*) "nprocs =",nprocs
    WRITE (6,*) "max_proc (in mp.inc) =",max_proc
    CALL arpsstop ("arpsstop called from mpinit_var mismatch in   &
                   & number of processors-too many",1)

  END IF

!
!  This subroutine defines the proc(nproc_x+nproc_x*(nproc_y-1)) array
!  and the myproc variable for each process.
!
  IF(nprocs /= nprocs0)THEN  ! test to see if the input file
                             ! number of processors = nprocs
                             ! and set on the command line.
    IF(myproc == 0)THEN
      PRINT *,'Number of processors chosen on the command line '
      PRINT *,'is different from that given in arps.input, EXITING'
      PRINT *,'requested: ', nprocs0
      PRINT *,'in arps.input: ', nprocs, ' = ',nproc_x,' * ',nproc_y
    END IF
    CALL arpsstop ("arpsstop called from mpinit_var mismatch in   &
                   & number of processors",1)

  END IF

  l = 0
  DO j = 1, nproc_y
    DO i = 1, nproc_x
      proc(i+nproc_x*(j-1)) = l
      l = l + 1
    END DO
  END DO

  loc_x = MOD(myproc, nproc_x) + 1
  loc_y = myproc / nproc_x + 1

  gentag = 0

  RETURN
END SUBROUTINE mpinit_var
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge(locvar,nx,ny,nz,nt,fzone,globvar,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global data files from a multiprocessor run to be compared
!  with a single processor file.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Dan Weber
!  2001/04/11
!
!  MODIFICATION HISTORY:
!  2004/04/02 (Yunheng Wang)
!  Added parameter fzone and globvar to make the code work more flexible.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Local variable
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    fzone    Number of fake zone, 3 for ARPS, 1 for WRF.
!
!  OUTPUT:
!
!    globvar  Global variable
!    tem1     Work array.
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in x, y and z
                               ! directions
  INTEGER :: nt
  INTEGER :: fzone
  REAL(P) :: locvar(nx,ny,nz,nt)
  REAL(P) :: tem1(nx,ny,nz)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  REAL(P) :: globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone,nz)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k
  INTEGER :: mptag             ! Unique MPI id used for this BC update.
  INTEGER :: ia,ja, ic,jc,itemc,itemb,itema

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  DO k=1,nz        ! each processor stores the locvar into tem1
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k) = locvar(i,j,k,nt)
      END DO
    END DO
  END DO

  DO jc=1,nproc_y
    DO ic=1,nproc_x

!   message passing section...

      itemb = mptag + 100 + ic + jc
      IF(ic /=1 .OR. jc /=1) THEN     !  pass data to processor 0

        IF( myproc .EQ. (ic+(jc-1)*nproc_x-1)) THEN
          itema = 0
          ! print *,'sending data',itema,itemb,myproc
          CALL mpi_send (tem1,nx*ny*nz,MPI_p,itema,                  &
                         itemb,MPI_COMM_WORLD,imstat)
        END IF

        itemc = ic+(jc-1)*nproc_x-1
        IF(myproc == 0)THEN            ! receive data
          ! print *,'receiving data',itemc,itemb,myproc
          CALL mpi_recv (tem1,nx*ny*nz,MPI_p,itemc,                  &
                         itemb,MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

!  storage section

      IF(myproc == 0)THEN  ! store data into globvar
        DO k=1,nz
          DO j=1,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx
              ia = i + (ic-1)*(nx-fzone)
              ! print *,ic,jc,ia,ja,i,j,k
              globvar(ia,ja,k) = tem1(i,j,k)
            END DO
          END DO
        END DO
      END IF

      call mpbarrier

    END DO
  END DO

!  IF(myproc ==0 ) THEN   !  write the file.....
!
!     write(char1(length+1:length+5),'(a5)')  '.form'
!!    itemc = 80
!!    CALL strlnth(char1,itemc)
!!    CALL comlnth(char1,itemc)
!!    print *,'inside mpimerge', length,char1(1:length+5)
!     open(10,file=char1(1:length+5),form= 'formatted',status='unknown')
!     DO k=1,nz
!     DO j=1,(ny-3)*nproc_y+3
!     DO i=1,(nx-3)*nproc_x+3
!     write(10,'(3(i5),2x,g17.11)') i,j,k,globvar(i,j,k)
!     END DO
!     END DO
!     END DO
!     close (10)
!
!     write(char1(length+1:length+7),'(a7)')  '.unform'
!!    itemc = 80
!!    CALL comlnth(char1,itemc)
!!    CALL strlnth(char1,itemc)
!!    print *,'inside mpimerge', itemc,length,char1(1:itemc)
!!    print *,'inside mpimerge', length,char1(1:length+7)
!     open(11,file=char1(1:length+7),form= 'unformatted',status='unknown')
!     write (11) globvar
!     close (11)
!
!  END IF

   RETURN
END SUBROUTINE mpimerge
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE1dx                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge1dx(locvar,nx,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Dimension of the array
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx
  REAL(P), INTENT(IN) :: locvar(nx)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  REAL(P), INTENT(OUT):: globvar((nx-3)*nproc_x+3)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia, ic,jc, i0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i

  INTEGER :: mpi_p

  REAL(P),ALLOCATABLE :: tem(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(nx), STAT = imstat)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:) = locvar(:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(jc == 1) THEN

        IF(ic /= 1) THEN     !  pass data to processor 0

          mytag = mptag + 100 + ic + 1
          IF( myproc == ic-1 )THEN
            CALL mpi_send (tem,nx,MPI_p,master,                     &
                           mytag,MPI_COMM_WORLD,imstat)

            !CALL mpi_bsend (tem,nx*ny*nz,MPI_REAL,master,             &
            !                mytag,MPI_COMM_WORLD,imstat)
                                                      !forced buffering
          END IF

          IF(myproc == 0)THEN          ! receive data
            source = ic - 1
            CALL mpi_recv (tem,nx,MPI_p,source,                    &
                           mytag, MPI_COMM_WORLD,stat,imstat)
          END IF

        END IF

        !  storage section

        IF(myproc == 0)THEN  ! store data into globvar

          IF (ic == 1) THEN
            i0 = 1
          ELSE
            i0 = 2
          END IF

          DO i=i0,nx
            ia = i + (ic-1)*(nx-fzone)
            globvar(ia) = tem(i)
          END DO

        END IF

      END IF    ! jc == 1

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge1dx
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE1dy                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge1dy(locvar,ny,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    ny       Dimension of the array
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ny
  REAL(P), INTENT(IN) :: locvar(ny)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  REAL(P), INTENT(OUT):: globvar((ny-3)*nproc_y+3)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ja, ic,jc, j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: j

  INTEGER :: mpi_p

  REAL(P), ALLOCATABLE :: tem(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(ny), STAT = imstat)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:) = locvar(:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic == 1) THEN

      IF(jc /= 1) THEN     !  pass data to processor 0

        mytag = mptag + 100 + jc + 1
        IF(myproc ==  (jc-1)*nproc_x )THEN
          CALL mpi_send (tem,ny,MPI_p,master,                     &
                         mytag,MPI_COMM_WORLD,imstat)

          !CALL mpi_bsend (tem,nx*ny*nz,MPI_REAL,master,             &
          !                mytag,MPI_COMM_WORLD,imstat)
                                                    !forced buffering
        END IF

        IF(myproc == 0)THEN          ! receive data
          source = (jc-1)*nproc_x
          CALL mpi_recv (tem,ny,MPI_p,source,                    &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        IF (jc == 1) THEN
          j0 = 1
        ELSE
          j0 = 2
        END IF

        DO j=j0,ny
          ja = j + (jc-1)*(ny-fzone)
          globvar(ja) = tem(j)
        END DO

      END IF

      END IF    ! ic == 1

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge1dy

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE2D                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge2d(locvar,nx,ny,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny
                            ! Number of grid points in x, y and z
  REAL(P), INTENT(IN) :: locvar(nx,ny)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'


  REAL(P), INTENT(OUT):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3)
                            ! Output array in global domain.
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, i0,j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j

  INTEGER :: mpi_p

  REAL(P), ALLOCATABLE :: tem(:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(nx,ny), STAT = imstat)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_send (tem,nx*ny,MPI_p,master,               &
                         mytag,MPI_COMM_WORLD,imstat)

          !CALL mpi_bsend (tem,nx*ny,MPI_REAL,master,             &
          !                mytag,MPI_COMM_WORLD,imstat)
                                                    !forced buffering
        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny,MPI_p,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        IF (ic == 1) THEN
          i0 = 1
        ELSE
          i0 = 2
        END IF

        IF (jc == 1) THEN
          j0 = 1
        ELSE
          j0 = 2
        END IF

        DO j=j0,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=i0,nx
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja) = tem(i,j)
            END DO
        END DO

      END IF

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE2di                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge2di(locvar,nx,ny,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny
                            ! Number of grid points in x, y and z
  INTEGER, INTENT(IN) :: locvar(nx,ny)

  INTEGER, INTENT(OUT):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, i0,j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j

  INTEGER,ALLOCATABLE :: tem(:,:)
! INTEGER :: tem(nx,ny)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(nx,ny), STAT = imstat)

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc
        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_send (tem,nx*ny,MPI_INTEGER,master,               &
                         mytag,MPI_COMM_WORLD,imstat)

          !CALL mpi_bsend (tem,nx*ny,MPI_INTEGER,master,             &
          !                mytag,MPI_COMM_WORLD,imstat)
                                                    !forced buffering
        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny,MPI_INTEGER,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        IF (ic == 1) THEN
          i0 = 1
        ELSE
          i0 = 2
        END IF

        IF (jc == 1) THEN
          j0 = 1
        ELSE
          j0 = 2
        END IF

        DO j=j0,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=i0,nx
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja) = tem(i,j)
            END DO
        END DO

      END IF

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE3d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!SUBROUTINE mpimerge3d(locvar,nx,ny,nz,istagger,globvar)
SUBROUTINE mpimerge3d(locvar,nx,ny,nz,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!  Yunheng Wang (05/10/2011)
!  Added staggering option
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: nx,ny,nz
                            ! Number of grid points in x, y and z
  REAL(P), INTENT(IN) :: locvar(nx,ny,nz)
  !INTEGER, INTENT(IN) :: istagger

  REAL(P), INTENT(OUT):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3,nz)
                            ! Output array in global domain.
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, i0,j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k

  INTEGER :: mpi_p  !, is, js

  REAL(P), ALLOCATABLE :: tem(:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(nx,ny,nz), STAT = imstat)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  !is = 0
  !js = 0
  !IF (istagger == 1) is = 1
  !IF (istagger == 2) js = 1

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_send (tem,nx*ny*nz,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,imstat)

          !CALL mpi_bsend (tem,nx*ny*nz,MPI_REAL,master,             &
          !                mytag,MPI_COMM_WORLD,imstat)
                                                    !forced buffering
        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nz,mpi_p,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        IF (ic == 1) THEN
          i0 = 1
        ELSE
          i0 = 2 !+ is
        END IF

        IF (jc == 1) THEN
          j0 = 1
        ELSE
          j0 = 2 !+ js
        END IF

        DO k=1,nz
          DO j=j0,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=i0,nx
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja,k) = tem(i,j,k)
            END DO
          END DO
        END DO

      END IF

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE3di                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge3di(locvar,nx,ny,nz,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz
                            ! Number of grid points in x, y and z
  INTEGER, INTENT(IN) :: locvar(nx,ny,nz)

  INTEGER, INTENT(OUT):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3,nz)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, i0,j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k

  INTEGER,ALLOCATABLE :: tem(:,:,:)
! INTEGER :: tem(nx,ny,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(nx,ny,nz), STAT = imstat)

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc
        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_send (tem,nx*ny*nz,MPI_INTEGER,master,               &
                         mytag,MPI_COMM_WORLD,imstat)

          !CALL mpi_bsend (tem,nx*ny*nz,MPI_INTEGER,master,             &
          !                mytag,MPI_COMM_WORLD,imstat)
                                                    !forced buffering
        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nz,MPI_INTEGER,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        IF (ic == 1) THEN
          i0 = 1
        ELSE
          i0 = 2
        END IF

        IF (jc == 1) THEN
          j0 = 1
        ELSE
          j0 = 2
        END IF

        DO k=1,nz
          DO j=j0,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=i0,nx
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja,k) = tem(i,j,k)
            END DO
          END DO
        END DO

      END IF

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge3di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE4d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpimerge4d(locvar,nx,ny,nzsoil,nstyps,globvar)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nzsoil, nstyps
                            ! Number of grid points in x, y and z
  REAL(P), INTENT(IN) :: locvar(nx,ny,nzsoil,nstyps)

  REAL(P), INTENT(OUT):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3,nzsoil,nstyps)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, i0,j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, n

  INTEGER :: mpi_p

  REAL(P), ALLOCATABLE :: tem(:,:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(nx,ny,nzsoil,nstyps), STAT = imstat)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:,:,:) = locvar(:,:,:,:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc
        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_send (tem,nx*ny*nzsoil*nstyps,mpi_p,master,    &
                         mytag,MPI_COMM_WORLD,imstat)

          !CALL mpi_bsend (tem,nx*ny*nz,MPI_REAL,master,             &
          !                mytag,MPI_COMM_WORLD,imstat)
                                                    !forced buffering
        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nzsoil*nstyps,mpi_p,source,    &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        IF (ic == 1) THEN
          i0 = 1
        ELSE
          i0 = 2
        END IF

        IF (jc == 1) THEN
          j0 = 1
        ELSE
          j0 = 2
        END IF

        DO n=1,nstyps
          DO k=1,nzsoil
            DO j=j0,ny
              ja = j + (jc-1)*(ny-fzone)
              DO i=i0,nx
                ia = i + (ic-1)*(nx-fzone)
                globvar(ia,ja,k,n) = tem(i,j,k,n)
              END DO
            END DO
          END DO
        END DO

      END IF

      CALL mpbarrier

    END DO
  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge4d
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT1DX                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit1dx(globvar,nx,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Dimension of the array in subdomain.
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx
  REAL(P), INTENT(IN) :: globvar((nx-3)*nproc_x+3)

  REAL(P), INTENT(OUT) :: var(nx)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: targetp
  INTEGER :: ia,ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into var

        DO i=1,nx
           ia = i + (ic-1)*(nx-fzone)
           var(i) = globvar(ia)
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /= 1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          targetp = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx,mpi_p,targetp,                 &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:) = globvar(1:nx)

  RETURN
END SUBROUTINE mpisplit1dx
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT1DY                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit1dy(globvar,ny,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: ny
  REAL(P), INTENT(IN)  :: globvar((ny-3)*nproc_y+3)

  REAL(P), INTENT(OUT) :: var(ny)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: targetp
  INTEGER :: ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: j

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into tem

          DO j=1,ny
            ja = j + (jc-1)*(ny-fzone)
            var(j) = globvar(ja)
          END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          targetp = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,ny,mpi_p,targetp,              &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,ny,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:) = globvar(1:ny)

  RETURN
END SUBROUTINE mpisplit1dy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT2d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit2d(globvar,nx,ny,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny   ! Number of grid points in x and y
  REAL(P), INTENT(IN)  :: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3)

  REAL(P), INTENT(OUT) :: var(nx,ny)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: dest
  INTEGER :: ia,ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into var

        DO j=1,ny
          ja = j + (jc-1)*(ny-fzone)
          DO i=1,nx
            ia = i + (ic-1)*(nx-fzone)
            var(i,j) = globvar(ia,ja)
          END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          dest = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny,mpi_p,dest,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx*ny,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:) = globvar(1:nx,1:ny)

  RETURN
END SUBROUTINE mpisplit2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT2DI                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit2di(globvar,nx,ny,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny
                            ! Number of grid points in x and y
  INTEGER, INTENT(IN):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3)

  INTEGER, INTENT(OUT) :: var(nx,ny)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: dest
  INTEGER :: ia,ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into tem

        DO j=1,ny
          ja = j + (jc-1)*(ny-fzone)
          DO i=1,nx
            ia = i + (ic-1)*(nx-fzone)
            var(i,j) = globvar(ia,ja)
          END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          dest = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny,MPI_INTEGER,dest,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx*ny,MPI_INTEGER,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:) = globvar(1:nx,1:ny)

  RETURN
END SUBROUTINE mpisplit2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT3d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit3d(globvar,nx,ny,nz,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       3rd dimension of the subdomain array, possible value are
!             vertical grid points (nz in other subroutines), nzsoil,
!             nstyps+1, or 4 (prcrate) or 1 (for 2D arrays)
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz
                            ! Number of grid points in x, y and z
  REAL(P), INTENT(IN)  :: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3,nz)

  REAL(P), INTENT(OUT) :: var(nx,ny,nz)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: dest
  INTEGER :: ia,ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into tem

        DO k=1,nz
          DO j=1,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx
              ia = i + (ic-1)*(nx-fzone)
              var(i,j,k) = globvar(ia,ja,k)
            END DO
          END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          dest = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny*nz,mpi_p,dest,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx*ny*nz,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:,:) = globvar(1:nx,1:ny,1:nz)

  RETURN
END SUBROUTINE mpisplit3d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT3DI                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit3di(globvar,nx,ny,nz,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       3rd dimension of the subdomain array, possible value are
!             vertical grid points (nz in other subroutines), nzsoil,
!             nstyps+1, or 4 (prcrate) or 1 (for 2D arrays)
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz
                            ! Number of grid points in x, y and z
  INTEGER, INTENT(IN):: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3,nz)

  INTEGER, INTENT(OUT) :: var(nx,ny,nz)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: dest
  INTEGER :: ia,ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into var

        DO k=1,nz
          DO j=1,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx
              ia = i + (ic-1)*(nx-fzone)
              var(i,j,k) = globvar(ia,ja,k)
            END DO
          END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          dest = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny*nz,MPI_INTEGER,dest,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx*ny*nz,MPI_INTEGER,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:,:) = globvar(1:nx,1:ny,1:nz)

  RETURN
END SUBROUTINE mpisplit3di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT4d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit4d(globvar,nx,ny,nzsoil,nstyps,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   3rd dimension of the subdomain array, possible value may be
!             as nzsoil in other subroutines.
!    nstyps   4rd dimentsion of the 4D array, possible value may be
!             nstyps (in other subroutines) + 1.
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nzsoil, nstyps
  REAL(P), INTENT(IN) :: globvar((nx-3)*nproc_x+3,(ny-3)*nproc_y+3,nzsoil,nstyps)

  REAL(P), INTENT(OUT):: var(nx,ny,nzsoil,nstyps)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: dest
  INTEGER :: ia,ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, l

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into tem

        DO l=1,nstyps
        DO k=1,nzsoil
          DO j=1,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx
              ia = i + (ic-1)*(nx-fzone)
              var(i,j,k,l) = globvar(ia,ja,k,l)
            END DO
          END DO
        END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          dest = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny*nzsoil*nstyps,mpi_p,dest,     &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx*ny*nzsoil*nstyps,mpi_p,master,     &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

  ! Finally, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:,:,:) = globvar(1:nx,1:ny,1:nzsoil,1:nstyps)

  RETURN
END SUBROUTINE mpisplit4d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE2dns               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge2dx(locvar,nx,nz,yp,globvar,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run.
!  for North/South boundary 2D array.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2003/02/25
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!  10/29/2009 (Y. Wang)
!  Expanded the capability to join in any row "yp" and also simplified
!  the code
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    nx       Local dimension of the array
!    nz       This dimension will not change.
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: yp
  REAL(P), INTENT(IN) :: locvar(nx,nz)
  INTEGER, INTENT(OUT):: istatus

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  REAL(P), INTENT(OUT):: globvar((nx-3)*nproc_x+3,nz)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  REAL(P), ALLOCATABLE :: tem(:,:)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ia, ic,jc, i0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i,k

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ALLOCATE( tem(nx,nz), STAT = istatus )

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem

  DO ic=1,nproc_x

    jc = yp
    source = (jc-1)*nproc_x + ic - 1

    ! message passing section...

    IF(source /= 0) THEN     !  pass data to processor 0

      mytag = mptag + 100 + ic + jc
      IF( myproc == source )THEN
        CALL mpi_send (tem,nx*nz,mpi_p,master,                          &
                       mytag,MPI_COMM_WORLD,imstat)

      END IF

      IF(myproc == 0)THEN          ! receive data
        CALL mpi_recv (tem,nx*nz,mpi_p,source,                          &
                       mytag, MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    !  storage section

    IF(myproc == 0)THEN  ! store data into globvar

      IF (ic == 1) THEN
        i0 = 1
      ELSE
        i0 = 2
      END IF

      DO k=1,nz
        DO i=i0,nx
          ia = i + (ic-1)*(nx-fzone)
          globvar(ia,k) = tem(i,k)
        END DO
      END DO

    END IF

    CALL mpbarrier

  END DO

  DEALLOCATE(tem)
  RETURN
END SUBROUTINE mpimerge2dx
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIMERGE2dew               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpimerge2dy(locvar,ny,nz,xp,globvar,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run
!  for East/West boundary 2D array.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2003/02/25
!  Based on subroutine mpimerge
!
!  MODIFICATION HISTORY:
!  10/29/2009 (Y. Wang)
!  Expanded the capability to join in any column "xp" and also simplified
!  the code!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    localvar Variable to be written.
!
!    ny       Dimension of the array
!    nz
!
!  OUTPUT:
!
!    globvar  global variable to be output
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: xp
  REAL(P), INTENT(IN) :: locvar(ny,nz)
  INTEGER, INTENT(OUT):: istatus

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  REAL(P), INTENT(OUT):: globvar((ny-3)*nproc_y+3,nz)
                            ! Output array in global domain.

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  REAL(P), ALLOCATABLE :: tem(:,:)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER, PARAMETER :: master = 0
  INTEGER :: source
  INTEGER :: ja, ic,jc, j0, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: j,k

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ALLOCATE(tem(ny,nz), STAT = istatus)

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3 !  arps.

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem

  DO jc=1,nproc_y

    ic = xp
    source = (jc-1)*nproc_x+ic-1     ! This processor to pass data to Processor 0

    ! message passing section...
    IF (source /= 0) THEN     !  pass data to processor 0

      mytag = mptag + 100 + jc + ic
      IF(myproc ==  source )THEN
        CALL mpi_send (tem,ny*nz,mpi_p,master,                  &
                       mytag,MPI_COMM_WORLD,imstat)
      END IF

      IF(myproc == 0)THEN          ! receive data
        CALL mpi_recv (tem,ny*nz,mpi_p,source,                 &
                       mytag, MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    !  storage section

    IF(myproc == 0)THEN  ! store data into globvar

      IF (jc == 1) THEN
        j0 = 1
      ELSE
        j0 = 2
      END IF

      DO k=1,nz
        DO j=j0,ny
          ja = j + (jc-1)*(ny-fzone)
          globvar(ja,k) = tem(j,k)
        END DO
      END DO

    END IF

    CALL mpbarrier

  END DO

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE mpimerge2dy
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT2DNS               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit2dns(globvar,nx,nz,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!  for North/South boundary arrays.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2003/02/26
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    nx       Dimension of the array in subdomain.
!    nz
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: nz
  REAL(P), INTENT(IN):: globvar((nx-3)*nproc_x+3,nz)

  REAL(P), INTENT(OUT) :: var(nx,nz)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: target
  INTEGER :: ia,ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i,k

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into var

        DO k=1,nz
        DO i=1,nx
           ia = i + (ic-1)*(nx-fzone)
           var(i,k) = globvar(ia,k)
        END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /= 1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          target = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*nz,mpi_p,target,                &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,nx*nz,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:) = globvar(1:nx,:)

  RETURN
END SUBROUTINE mpisplit2dns
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPISPLIT2DEW               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpisplit2dew(globvar,ny,nz,var)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the global array and scatter to each processor.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2002/08/20
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    globvar  Global array passed in from processor 0.
!
!    ny       Number of grid points in the y-direction (north/south)
!    nz
!
!  OUTPUT:
!
!    var      Subdomain variable in each process.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  REAL(P), INTENT(IN):: globvar((ny-3)*nproc_y+3,nz)

  REAL(P), INTENT(OUT) :: var(ny,nz)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: target
  INTEGER :: ja, ic,jc, fzone
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: j,k

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!
!  fill the globvar array
!
!-----------------------------------------------------------------------

  fzone = 3        !  arps.

  DO jc=1,nproc_y
    DO ic=1,nproc_x

      !  storage section

      IF(myproc == 0)THEN  ! store data into var

        DO k = 1,nz
          DO j=1,ny
            ja = j + (jc-1)*(ny-fzone)
            var(j,k) = globvar(ja,k)
          END DO
        END DO

      END IF

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  receive data from processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == 0)THEN          ! send data

          target = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,ny*nz,mpi_p,target,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN

          CALL mpi_recv (var,ny*nz,mpi_p,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat)

        END IF

      END IF

      CALL mpbarrier

    END DO
  END DO

  !At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:) = globvar(1:ny,:)

  RETURN
END SUBROUTINE mpisplit2dew
!
!
!######################################################################
!
! Wrap subroutines added for ARPSPLT_mpi
!
!   mpsendr  -- CALL mpi_send, REAL array
!   mprecvr  -- CALL mpi_recv, REAL array
!   mpsendi  -- CALL mpi_send, INTEGER array (Used to be scalar)
!   mprecvi  -- CALL mpi_recv, INTEGER array (Used to be scalar)
!
!#####################################################################

SUBROUTINE mpsendr(a,isize,dest,itag,ierror)

  USE arps_precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: isize, dest, itag, ierror
  REAL(P) :: a(isize)

  INTEGER :: mpi_p

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_send(a,isize,mpi_p,dest,itag,MPI_COMM_WORLD,ierror)

!  WRITE(6,'(1x,a,I3)') 'ERROR: inside mpsendr, ierror = ',ierror

  RETURN
END SUBROUTINE mpsendr

SUBROUTINE mpsendi(m,isize,dest,itag,ierror)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: isize, dest, itag, ierror
  INTEGER :: m(isize)

  CALL mpi_send(m,isize,MPI_INTEGER,dest,itag,MPI_COMM_WORLD,ierror)

  RETURN

END SUBROUTINE mpsendi

SUBROUTINE mprecvr(a,isize,source,itag,ierror)

  USE arps_precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: isize, source, itag, ierror
  REAL(P) :: a(isize)

  INTEGER :: stat(MPI_STATUS_SIZE)

  INTEGER :: mpi_p

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_recv(a,isize,mpi_p,source,itag,MPI_COMM_WORLD,stat,ierror)

!  WRITE(6,'(1x,a,I3)') 'ERROR: inside mprecvr, ierror = ',ierror

  RETURN

END SUBROUTINE mprecvr

SUBROUTINE mprecvi(m,isize,source,itag,ierror)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: isize, source, itag, ierror
  INTEGER :: m(isize)

  INTEGER :: stat(MPI_STATUS_SIZE)

  CALL mpi_recv(m,isize,MPI_INTEGER,source,itag,MPI_COMM_WORLD,stat,ierror)

  RETURN

END SUBROUTINE mprecvi

SUBROUTINE mpmaxi(imax)
!##################################################################
!
!  Find the maximum integer of all processors
!
!##################################################################

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: imax

!-----------------------------------------------------------------

  INTEGER :: imstat

  INCLUDE 'mpif.h'

  INTEGER :: maxtm

!---------------------------------------------------------
!
! Start of executable code....
!
!---------------------------------------------------------

  CALL MPI_REDUCE(imax,maxtm,1,MPI_INTEGER,MPI_MAX,0,     &
                  MPI_COMM_WORLD,imstat)
  CALL MPI_BCAST(maxtm,1,MPI_INTEGER,0,MPI_COMM_WORLD,imstat)

  imax = maxtm

  RETURN
END SUBROUTINE mpmaxi
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPBCASTR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpbcastr(var,source)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast a real value from source processor to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2003/07/31
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    var      Real value to broadcast
!    source   source processor rank
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: source
  REAL(P), INTENT(IN) :: var

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_bcast(var,1,mpi_p,source,mpi_comm_world,imstat)

  IF (imstat /= 0) THEN
    WRITE (6,*) "MPBCASTR: error on processor",myproc
  END IF

  RETURN
END SUBROUTINE mpbcastr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPBCASTI                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpbcasti(var,source)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast an integer value from source processor to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2010/06/16
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    var      Real value to broadcast
!    source   source processor rank
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: source
  INTEGER, INTENT(IN) :: var

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpi_bcast(var,1,MPI_INTEGER,source,mpi_comm_world,imstat)

  IF (imstat /= 0) THEN
    WRITE (6,'(1x,a,I4)') 'MPBCASTI: error on processor: ',myproc
  END IF

  RETURN
END SUBROUTINE mpbcasti
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPBCASTRA                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpbcastra(var,num,source)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Broadcast the value of var from process source to all other processes.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Y. Wang
!  2010/06/17
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT :
!
!    var      Array to update (INPUT on proc 0, OUTPUT for rest).
!    num      Number of elements in the array.
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: num, source
  REAL(P), INTENT(INOUT) :: var(num)

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: imstat
  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_bcast(var,num,mpi_p,source,mpi_comm_world,imstat)
  IF (imstat /= 0) THEN
    WRITE (6,*) 'MPBCASTRA: error on processor: ',source
  END IF

  RETURN
END SUBROUTINE mpbcastra
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSUMI                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsumi(ivar,ndim)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ndim
  INTEGER :: ivar(ndim)
  INTEGER :: ivartm(ndim)
  INTEGER :: imstat

  CALL mpi_allreduce(ivar, ivartm, ndim, MPI_INTEGER, MPI_SUM,            &
                     mpi_comm_world, imstat)

  ivar(1:ndim) = ivartm(1:ndim)

  RETURN
END SUBROUTINE mpsumi
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSUMR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsumr(var,ndim)

  USE arps_precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ndim
  REAL(P) :: var(ndim)
  REAL(P) :: vartm(ndim)
  INTEGER :: imstat
  INTEGER :: mpi_p

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_allreduce(var, vartm, ndim, mpi_p, MPI_SUM,               &
                     mpi_comm_world, imstat)

  var(1:ndim) = vartm(1:ndim)

  RETURN
END SUBROUTINE mpsumr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPSUMDP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpsumdp(var,ndim)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER,          INTENT(IN)    :: ndim
  DOUBLE PRECISION, INTENT(INOUT) :: var(ndim)

!-----------------------------------------------------------------------

  DOUBLE PRECISION    :: vartm(ndim)
  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpi_allreduce(var, vartm, ndim, MPI_DOUBLE_PRECISION, MPI_SUM,   &
                     mpi_comm_world, imstat)

  var(1:ndim) = vartm(1:ndim)

  RETURN
END SUBROUTINE mpsumdp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPMINR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpminr(var)

  USE arps_precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(P)   :: var
  REAL(P)   :: vartm

  INTEGER :: imstat

  INTEGER :: mpi_p

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_allreduce(var, vartm, 1, mpi_p, MPI_MIN,   &
                     mpi_comm_world, imstat)

  var = vartm

  RETURN
END SUBROUTINE mpminr

SUBROUTINE mpmini(var)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: var
  INTEGER :: vartm

  INTEGER :: imstat

  CALL mpi_allreduce(var, vartm, 1, MPI_INTEGER, MPI_MIN,   &
                     mpi_comm_world, imstat)

  var = vartm

  RETURN
END SUBROUTINE mpmini
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE mpmaxr                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpmaxr(var)

  USE arps_precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(P)   :: var
  REAL(P)   :: vartm

  INTEGER :: imstat

  INTEGER :: mpi_p

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_allreduce(var, vartm, 1, mpi_p, MPI_MAX,   &
                     mpi_comm_world, imstat)

  var = vartm

  RETURN
END SUBROUTINE mpmaxr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE mpgatheri                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpgatheri(sendbuf,sendcount,recvbuf,recvcount,istatus)
! A simplified version of mpigatheri

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, INTENT(IN)  :: sendcount, recvcount
  INTEGER, INTENT(IN)  :: sendbuf(sendcount)
  INTEGER, INTENT(OUT) :: recvbuf(recvcount)   ! large array, only significant with root process
  INTEGER, INTENT(OUT) :: istatus

!----------------------------------------------------------------------

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpi_barrier( mpi_comm_world, istatus)

  CALL mpi_gather(sendbuf,sendcount,MPI_INTEGER,                        &
                  recvbuf,sendcount,MPI_INTEGER,                        &
                  0, MPI_COMM_WORLD, istatus)

  RETURN
END SUBROUTINE mpgatheri

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE mpscatteri                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpscatteri(sendbuf,sendcount,recvbuf,recvcount,istatus)

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, INTENT(IN)  :: sendcount, recvcount
  INTEGER, INTENT(IN)  :: sendbuf(sendcount)  ! larg array, only significant with root process
  INTEGER, INTENT(OUT) :: recvbuf(recvcount)
  INTEGER, INTENT(OUT) :: istatus

!----------------------------------------------------------------------

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  CALL mpi_scatter(sendbuf,recvcount,MPI_INTEGER,                       &
                   recvbuf,recvcount,MPI_INTEGER,                       &
                   0, MPI_COMM_WORLD, istatus)

  RETURN
END SUBROUTINE mpscatteri
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIGATHERi                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpigatheri(vari,isize,outvari,osize,                         &
                      lsize,ldisp,ROOT,numprocs,istatus)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! PURPOSE:
!   CALL mpi_gather & mpi_gatherv for INTEGER. First used with 88d2arps
!   from wrtgridtilt.
!
!-----------------------------------------------------------------------
!
! MODIFICATION HISTORY:
!
! 11/29/2010 (Y. Wang)
! Wrapped MPI_GATHER & MPI_GATHERV calls here
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN)  :: isize, osize
  INTEGER, INTENT(IN)  :: vari(isize)
  INTEGER, INTENT(IN)  :: ROOT, numprocs
  INTEGER, INTENT(IN)  :: lsize(numprocs), ldisp(numprocs)
  INTEGER, INTENT(OUT) :: outvari(osize)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (ANY(ldisp > 0) ) THEN
    CALL MPI_Gatherv(vari,isize,MPI_INTEGER,                            &
                     outvari,lsize,ldisp,MPI_INTEGER,                   &
                     ROOT,MPI_COMM_WORLD,istatus)
  ELSE
    CALL MPI_Gather(vari,isize,MPI_INTEGER,outvari,isize,MPI_INTEGER,   &
                    ROOT,MPI_COMM_WORLD,istatus)
  END IF

  RETURN
END SUBROUTINE mpigatheri

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPIGATHERR                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpigatherr(varr,isize,outvarr,osize,                         &
                      lsize,ldisp,ROOT,numprocs,istatus)
  USE arps_precision

  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! PURPOSE:
!   CALL mpi_gather & mpi_gatherv for REAL or DOUBLE PRECISION. First
!   used with 88d2arps from wrtgridtilt.
!
!-----------------------------------------------------------------------
!
! MODIFICATION HISTORY:
!
! 11/29/2010 (Y. Wang)
! Wrapped MPI_GATHER & MPI_GATHERV calls here
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN)  :: isize, osize
  REAL(P), INTENT(IN)  :: varr(isize)
  INTEGER, INTENT(IN)  :: ROOT, numprocs
  INTEGER, INTENT(IN)  :: lsize(numprocs), ldisp(numprocs)
  REAL(P), INTENT(OUT) :: outvarr(osize)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL mpi_barrier( mpi_comm_world, istatus)

  IF (ANY(ldisp > 0) ) THEN
    !write(6,*) 'Process: ',myproc,' sending ',isize
    !if (myproc == ROOT) write(6,*) 'Process: ',myproc,' receiving ',lsize, ldisp
    !call flush(6)
    CALL MPI_Gatherv(varr,isize,MPI_P,                            &
                     outvarr,lsize,ldisp,MPI_P,                   &
                     ROOT,MPI_COMM_WORLD,istatus)
  ELSE
    CALL MPI_Gather(varr,isize,MPI_P,                             &
                    outvarr,isize,MPI_P,                          &
                    ROOT,MPI_COMM_WORLD,istatus)
  END IF

  RETURN
END SUBROUTINE mpigatherr

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE mpmaxr                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpmaxar(var,isize)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: isize
  REAL(P)   :: var(isize)

!---------------------------------------------------------------------

  REAL(P), ALLOCATABLE  :: vartm(:)

  INTEGER :: imstat

  INTEGER :: mpi_p

  INCLUDE 'mpif.h'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(vartm(isize), STAT = imstat)

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL MPI_ALLREDUCE(var, vartm, isize, mpi_p, MPI_MAX,               &
                     MPI_COMM_WORLD, imstat)

  var = vartm

  DEALLOCATE(vartm)

  RETURN
END SUBROUTINE mpmaxar

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE mpminr                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpminar(var,isize)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: isize
  REAL(P)   :: var(isize)

!---------------------------------------------------------------------

  REAL(P), ALLOCATABLE  :: vartm(:)

  INTEGER :: imstat

  INTEGER :: mpi_p

  INCLUDE 'mpif.h'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(vartm(isize), STAT = imstat)

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  CALL MPI_ALLREDUCE(var, vartm, isize, mpi_p, MPI_MIN,               &
                     MPI_COMM_WORLD, imstat)

  var = vartm

  DEALLOCATE(vartm)


  RETURN
END SUBROUTINE mpminar
