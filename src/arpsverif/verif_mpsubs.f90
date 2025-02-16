!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE VERIF_COLLECT             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE VERIF_COLLECT(model_data,obsrv_data,nhisfile)
               
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Collect verification computations from other processors and merge it
!  into one data set.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Kevin W. Thomas
!
!  Original Coding: 08/02/05
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'vericst.inc'
  INCLUDE 'mp.inc'
  INCLUDE 'mpif.h'

  INTEGER :: nhisfile

  REAL :: model_data(sfcmax,nhisfile,5)
  REAL :: obsrv_data(sfcmax,nhisfile,5)
  REAL :: tem1(sfcmax,nhisfile,5)               ! work array
  REAL :: tem2(sfcmax,nhisfile,5)               ! work array
  CHARACTER(LEN=4) :: sfcstid_tmp(sfcmax)       ! work array
  INTEGER :: sfcstn_lcl

  INTEGER, allocatable :: kount(:)
  INTEGER :: mpi_status(MPI_STATUS_SIZE)

!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istat
  INTEGER :: i,j,k

!
!  Collect the number of data points and station lists for each processor.
!  The rest of the data will be collected last once the master list is
!  updated to remove stations that were outside of the domain.
!

  IF (myproc == 0) THEN

    ALLOCATE(kount(nprocs), STAT=istat)
    CALL check_alloc_status(istat,"verif_collect:kount")

!
!  Save processor 0 data.
!

    sfcstid_tmp(1:sfcstn) = sfcstid(1:sfcstn)

    kount(1) = sfcstn
    k = kount(1)

    DO i=1,nprocs-1

      CALL mpi_recv(kount(i+1),1,MPI_INTEGER,i,100+i,                       &
        MPI_COMM_WORLD,mpi_status,istat)

      CALL mpi_recv(sfcstid_tmp(k+1),kount(i+1)*4,MPI_CHARACTER,i,200+i,    &
        MPI_COMM_WORLD, mpi_status,istat)

      k = k + kount(i+1)

    END DO

!
!  Compute master list.  Output goes to "sfcstid", with "sfcstn" now
!  representing the global domain instead of the local domain.  The "master"
!  variables are not referenced after this point.
!

    sfcstn = 0

    DO j=1,sfcstn_master
      DO i=1,k
        IF (sfcstid_master(j) == sfcstid_tmp(i)) THEN
          sfcstn = sfcstn + 1
          sfcstid(sfcstn) = sfcstid_master(j)
          exit
        END IF
      END DO
   END DO

  ELSE

    CALL mpi_send(sfcstn,1,MPI_INTEGER,0,100+myproc,MPI_COMM_WORLD,istat)

    CALL mpi_send(sfcstid,sfcstn*4,MPI_CHARACTER,0,200+myproc,           &
      MPI_COMM_WORLD,istat)

  END IF


!
!  Collect the rest of the data.
!

  IF (myproc == 0) THEN

!
!  Copy what we've already computed so the data can be sorted.
!

    tem1(:,:,:) = model_data(:,:,:)
    tem2(:,:,:) = obsrv_data(:,:,:)

!
!  Sort what processor 0 has already saved.
!

     CALL verif_sort(model_data,obsrv_data,tem1,tem2,sfcmax,nhisfile,5,  &
        sfcstid,sfcstn,sfcstid_tmp,kount(1))
!    write(6,*)'user list:  ',(sfcstid_tmp(i),i=1,sfcstn)

    k = kount(1)

    DO i=1,nprocs-1

!
!  Receive the model data.
!

      CALL mpi_recv(tem1,sfcmax*nhisfile*5,MPI_REAL,i,300+i,             &
        MPI_COMM_WORLD,mpi_status,istat)

!
!  Receive the observations.
!

      CALL mpi_recv(tem2,sfcmax*nhisfile*5,MPI_REAL,i,400+i,              &
         MPI_COMM_WORLD,mpi_status,istat)

!
!  Put the data in our operational arrays.
!

!     write(6,*) 'test:  ',i,kount(i)+1,kount(i+1)
      CALL verif_sort(model_data,obsrv_data,tem1,tem2,sfcmax,nhisfile,5,   &
        sfcstid,sfcstn,sfcstid_tmp(k+1),kount(i+1))
!     write(6,*)'user list:  ',(sfcstid_tmp(j),j=1,sfcstn)

      k = k + kount(i+1)

    END DO

  ELSE

    CALL mpi_send(model_data,sfcmax*nhisfile*5,MPI_REAL,0,300+myproc, &
      MPI_COMM_WORLD,istat)

    CALL mpi_send(obsrv_data,sfcmax*nhisfile*5,MPI_REAL,0,400+myproc, &
       MPI_COMM_WORLD,istat)

  END IF

  IF (myproc == 0) THEN
    DEALLOCATE(kount)
  END IF

  RETURN
END SUBROUTINE verif_collect
