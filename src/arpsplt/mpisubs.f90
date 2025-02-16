!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPMAX0i                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mpmax0i(amax,amin)

  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! Get global maximum and minimux for Integer scalars.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(INOUT) :: amax,amin

  INTEGER :: imstat
  INTEGER :: maxtm, mintm

  INCLUDE 'mpif.h'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    start of executable code....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Get maximum from all processors
!
!-----------------------------------------------------------------------

!  CALL mpi_allreduce (amax, maxtm, 1, MPI_REAL, MPI_MAX,
!    :     MPI_COMM_WORLD, imstat)  ! commented out because the T3E
                                 ! has trouble with mpi_allreduce

  CALL mpi_reduce(amax,maxtm,1,MPI_INTEGER,MPI_MAX,0,                   &
                  MPI_COMM_WORLD,imstat)
  CALL mpi_bcast(maxtm,1,MPI_INTEGER,0,MPI_COMM_WORLD,imstat)
  amax = maxtm

!-----------------------------------------------------------------------
!
!  Get minimum from all processors
!
!-----------------------------------------------------------------------

!  CALL mpi_allreduce (amin, mintm, 1, MPI_REAL, MPI_MIN,
!    :     MPI_COMM_WORLD, imstat)  ! commented out because the T3E
                                 ! has trouble with mpi_allreduce

  CALL mpi_reduce(amin,mintm,1,MPI_INTEGER,MPI_MIN,0,                   &
                  MPI_COMM_WORLD,imstat)
  CALL mpi_bcast(mintm,1,MPI_INTEGER,0,MPI_COMM_WORLD,imstat)

  amin = mintm

  RETURN
END SUBROUTINE mpmax0i
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE globalpbar                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE globalpbar(pbarmax,ini,inj,klvl,zpc,nx,ny,nz,zpcmax)

!-----------------------------------------------------------------------
!
! Find global maximum pbarmax and its index, ini, inj
! and extract the zpc value from a 3d array at that location
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL,    INTENT(INOUT)  :: pbarmax
  INTEGER, INTENT(INOUT)  :: ini
  INTEGER, INTENT(INOUT)  :: inj
  INTEGER, INTENT(IN)     :: klvl
  INTEGER, INTENT(IN)     :: nx,ny,nz
  REAL,    INTENT(IN)     :: zpc(nx,ny,nz)
  REAL,    INTENT(OUT)    :: zpcmax

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  REAL    :: maxarr(2), maxtm(2)
  INTEGER :: maxsrc

  INTEGER :: imstat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code below ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  maxtm(1)  = 0.0
  maxtm(2)  = 0.0
  maxarr(1) = pbarmax
  maxarr(2) = FLOAT(myproc)

  ! should call mpi_allreduce, however, since T3E has trouble with this
  ! call, we use two calls below to substitute it.
  CALL mpi_reduce(maxarr,maxtm,1,MPI_2REAL,MPI_MAXLOC,0,                &
                  MPI_COMM_WORLD,imstat)
  CALL mpi_bcast (maxtm,1,MPI_2REAL,0,MPI_COMM_WORLD,imstat)

  pbarmax = maxtm(1)
  maxsrc  = NINT(maxtm(2))

  IF (myproc == maxsrc) THEN  ! only processor maxsrc contains what we want.
    IF (ini /= 0 .AND. inj /= 0) THEN
      zpcmax = zpc(ini,inj,klvl)
    ELSE
      zpcmax = -9999.0   ! missing value, will not be used
    END IF
  END IF
  CALL mpbcasti(ini,maxsrc)
  CALL mpbcasti(inj,maxsrc)
  CALL mpbcastr(zpcmax,maxsrc)

  RETURN
END SUBROUTINE globalpbar
