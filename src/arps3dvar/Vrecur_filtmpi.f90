!-----------------------------------------------------------------------
!
! PURPOSE:
!    Message passing part for the parallel version of 
!    the 3D recursive filter.
!
! CONTAINED SUBROUTINES:
!
!    SUBROUTINE send_next_bdyxs(tem1,nx,ny,nz,bdyxs)
!    SUBROUTINE receive_bdyxs(tem1,nx,ny,nz,bdyxs)
!    SUBROUTINE send_previous_bdyxe(tem1,nx,ny,nz,bdyxe)
!    SUBROUTINE receive_bdyxe(tem1,nx,ny,nz,bdyxe)
!    SUBROUTINE send_up_bdyys(tem1,nx,ny,nz,bdyys)
!    SUBROUTINE receive_bdyys(tem1,nx,ny,nz,bdyys)
!    SUBROUTINE send_down_bdyye(tem1,nx,ny,nz,bdyye)
!    SUBROUTINE receive_bdyye(tem1,nx,ny,nz,bdyye)
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Yunheng Wang, CAPS, OU, 11/01/2007
!
!-----------------------------------------------------------------------

SUBROUTINE send_next_bdyxs(tem1,nx,ny,nz,bdyxs)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny, nz
  REAL(P), INTENT(IN)  :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT) :: bdyxs(ny,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: dest
  INTEGER :: mptag

  INTEGER :: j,k
  INTEGER :: imstat
  INTEGER :: ireq

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  CALL inctag

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

  mptag = gentag

  dest = (loc_x+1) + (loc_y-1) * nproc_x - 1   ! next processor in the same row
  
  DO k = 1,nz
    DO j = 1,ny
      bdyxs(j,k) = tem1((nx-2),j,k)   ! no stag
    END DO
  END DO

!write(0,'(a,I2,a,I2,a,I6)') 'Sending message to ',dest,' from ',myproc,' with tag=',mptag
  CALL mpi_isend(bdyxs,ny*nz,mpi_P,dest,mptag,MPI_COMM_WORLD,ireq,imstat)

  RETURN
END SUBROUTINE send_next_bdyxs

SUBROUTINE receive_bdyxs(tem1,nx,ny,nz,bdyxs)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL(P), INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT)   :: bdyxs(ny,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: source
  INTEGER :: mptag

  INTEGER :: j,k
  INTEGER :: imstat(MPI_STATUS_SIZE)
  INTEGER :: ierr

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag
  
  source = (loc_x-1) + (loc_y -1) * nproc_x - 1

!write(0,'(a,I2,a,I2,a,I6)') 'Recving message from ',source,' by ',myproc,' with tag=',mptag
  CALL mpi_recv(bdyxs,ny*nz,mpi_p,source,mptag,MPI_COMM_WORLD,imstat,ierr)

  DO k = 1,nz
    DO j = 1,ny
      tem1(1,j,k) = bdyxs(j,k)
    END DO
  END DO
  
  RETURN
END SUBROUTINE receive_bdyxs

! XE
SUBROUTINE send_previous_bdyxe(tem1,nx,ny,nz,bdyxe)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL(P), INTENT(IN) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT) :: bdyxe(ny,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: dest
  INTEGER :: mptag

  INTEGER :: j,k
  INTEGER :: imstat
  INTEGER :: ireq

  INTEGER :: mpi_p

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag

  dest = (loc_x-1) + (loc_y-1) * nproc_x - 1
  
  DO k = 1,nz
    DO j = 1,ny
      bdyxe(j,k) = tem1(2,j,k)   ! no stag
    END DO
  END DO

!write(0,'(a,I2,a,I2,a,I6)') 'Sending message to ',dest,' from ',myproc,' with tag=',mptag
  CALL mpi_isend(bdyxe,ny*nz,mpi_p,dest,mptag,MPI_COMM_WORLD,ireq,imstat)

  RETURN
END SUBROUTINE send_previous_bdyxe

SUBROUTINE receive_bdyxe(tem1,nx,ny,nz,bdyxe)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL(P), INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT)   :: bdyxe(ny,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: source
  INTEGER :: mptag

  INTEGER :: mpi_p

  INTEGER :: j,k
  INTEGER :: imstat(MPI_STATUS_SIZE)
  INTEGER :: ierr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag
  
  source = (loc_x+1) + (loc_y -1) * nproc_x - 1

!write(0,'(a,I2,a,I2,a,I6)') 'Recving message from ',source,' by ',myproc,' with tag=',mptag
  CALL mpi_recv(bdyxe,ny*nz,mpi_p,source,mptag,MPI_COMM_WORLD,imstat,ierr)

  DO k = 1,nz
    DO j = 1,ny
      tem1(nx-1,j,k) = bdyxe(j,k)
    END DO
  END DO
  
  RETURN
END SUBROUTINE receive_bdyxe

SUBROUTINE send_up_bdyys(tem1,nx,ny,nz,bdyys)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL(P), INTENT(IN) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT) :: bdyys(nx,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: dest
  INTEGER :: mptag

  INTEGER :: mpi_p

  INTEGER :: i,k
  INTEGER :: imstat
  INTEGER :: ireq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag

  dest = loc_x + (loc_y) * nproc_x - 1
  
  DO k = 1,nz
    DO i = 1,nx
      bdyys(i,k) = tem1(i,(ny-2),k)   ! no stag
    END DO
  END DO

!write(0,'(a,I2,a,I2,a,I6)') 'Sending message to ',dest,' from ',myproc,' with tag=',mptag
  CALL mpi_isend(bdyys,nx*nz,mpi_p,dest,mptag,MPI_COMM_WORLD,ireq,imstat)

  RETURN
END SUBROUTINE send_up_bdyys

SUBROUTINE receive_bdyys(tem1,nx,ny,nz,bdyys)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL(P), INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT)   :: bdyys(nx,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: source
  INTEGER :: mptag

  INTEGER :: mpi_p

  INTEGER :: i,k
  INTEGER :: imstat(MPI_STATUS_SIZE)
  INTEGER :: ierr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag
  
  source = loc_x + (loc_y -2) * nproc_x - 1

!write(0,'(a,I2,a,I2,a,I6)') 'Recving message from ',source,' by ',myproc,' with tag=',mptag
  CALL mpi_recv(bdyys,nx*nz,mpi_p,source,mptag,MPI_COMM_WORLD,imstat,ierr)

  DO k = 1,nz
    DO i = 1,nx
      tem1(i,1,k) = bdyys(i,k)
    END DO
  END DO
  
  RETURN
END SUBROUTINE receive_bdyys

! yE
SUBROUTINE send_down_bdyye(tem1,nx,ny,nz,bdyye)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL(P), INTENT(IN) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT) :: bdyye(nx,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: dest
  INTEGER :: mptag

  INTEGER :: mpi_p

  INTEGER :: i,k
  INTEGER :: imstat
  INTEGER :: ireq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag

  dest = loc_x + (loc_y-2) * nproc_x - 1
  
  DO k = 1,nz
    DO i = 1,nx
      bdyye(i,k) = tem1(i,2,k)   ! no stag
    END DO
  END DO

!write(0,'(a,I2,a,I2,a,I6)') 'Sending message to ',dest,' from ',myproc,' with tag=',mptag

  CALL mpi_isend(bdyye,nx*nz,mpi_p,dest,mptag,MPI_COMM_WORLD,ireq,imstat)

  RETURN
END SUBROUTINE send_down_bdyye

SUBROUTINE receive_bdyye(tem1,nx,ny,nz,bdyye)

  USE arps_precision

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL(P), INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL(P), INTENT(OUT)   :: bdyye(nx,nz)

  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: source
  INTEGER :: mptag

  INTEGER :: mpi_p

  INTEGER :: i,k
  INTEGER :: imstat(MPI_STATUS_SIZE)
  INTEGER :: ierr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  mpi_p = MPI_REAL
  IF (P == DP) mpi_p = MPI_DOUBLE_PRECISION

!  CALL inctag
  mptag = gentag
  
  source = loc_x + loc_y * nproc_x - 1

!write(0,'(a,I2,a,I2,a,I6)') 'Recving message from ',source,' by ',myproc,' with tag=',mptag
  CALL mpi_recv(bdyye,nx*nz,mpi_p,source,mptag,MPI_COMM_WORLD,imstat,ierr)

  DO k = 1,nz
    DO i = 1,nx
      tem1(i,ny-1,k) = bdyye(i,k)
    END DO
  END DO
  
  RETURN
END SUBROUTINE receive_bdyye
