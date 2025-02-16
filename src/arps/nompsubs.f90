
SUBROUTINE mpinit_proc
  RETURN
END SUBROUTINE mpinit_proc

SUBROUTINE mpinit_var
  IMPLICIT NONE

  INCLUDE 'mp.inc'

  myproc  = 0
  nproc_x = 1
  nproc_y = 1
  nprocs  = 1
  loc_x   = 1
  loc_y   = 1

  max_fopen = 1
  readsplit(:) = 0
  joindmp(:)   = 0
  readstride = 1
  dumpstride = 1

  RETURN
END SUBROUTINE mpinit_var

SUBROUTINE mpsendrecv1dew
  RETURN
END SUBROUTINE mpsendrecv1dew

SUBROUTINE mpsendrecv1dns
  RETURN
END SUBROUTINE mpsendrecv1dns

SUBROUTINE mpsendrecv1diew
  RETURN
END SUBROUTINE mpsendrecv1diew

SUBROUTINE mpsendrecv1dins
  RETURN
END SUBROUTINE mpsendrecv1dins

SUBROUTINE mpsendrecvextew
  RETURN
END SUBROUTINE mpsendrecvextew

SUBROUTINE mpsendrecvextns
  RETURN
END SUBROUTINE mpsendrecvextns

SUBROUTINE mpsendrecv2dew
  RETURN
END SUBROUTINE mpsendrecv2dew

SUBROUTINE mpsendrecv2dns
  RETURN
END SUBROUTINE mpsendrecv2dns

SUBROUTINE mpfillextarray3d
  RETURN
END SUBROUTINE mpfillextarray3d

SUBROUTINE mpupdater
  RETURN
END SUBROUTINE mpupdater

SUBROUTINE mpupdatei
  RETURN
END SUBROUTINE mpupdatei

SUBROUTINE mpupdatec
  RETURN
END SUBROUTINE mpupdatec

SUBROUTINE mpupdatel
  RETURN
END SUBROUTINE mpupdatel

SUBROUTINE mpexit(errcode)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: errcode
  STOP
  RETURN
END SUBROUTINE mpexit

SUBROUTINE inctag
  RETURN
END SUBROUTINE inctag

SUBROUTINE mpbarrier
  RETURN
END SUBROUTINE mpbarrier

SUBROUTINE mptotal
  RETURN
END SUBROUTINE mptotal

SUBROUTINE mpmax0
  RETURN
END SUBROUTINE mpmax0

SUBROUTINE mpmax
  RETURN
END SUBROUTINE mpmax

SUBROUTINE mpimerge
  RETURN
END SUBROUTINE mpimerge

SUBROUTINE mpimerge1dx(locvar,nx,globvar)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx
  REAL,    INTENT(IN) :: locvar(nx)
  REAL,    INTENT(OUT):: globvar(nx)

  globvar(:) = locvar(:)

  RETURN
END SUBROUTINE mpimerge1dx

SUBROUTINE mpimerge1dy(locvar,ny,globvar)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ny
  REAL,    INTENT(IN) :: locvar(ny)
  REAL,    INTENT(OUT):: globvar(ny)

  globvar(:) = locvar(:)

  RETURN
END SUBROUTINE mpimerge1dy

SUBROUTINE mpimerge2d(locvar,nx,ny,globvar)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny
  REAL,    INTENT(IN) :: locvar(nx,ny)
  REAL,    INTENT(OUT):: globvar(nx,ny)

  globvar(:,:) = locvar(:,:)

  RETURN
END SUBROUTINE mpimerge2d

SUBROUTINE mpimerge2di(locvar,nx,ny,globvar)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: locvar(nx,ny)
  INTEGER, INTENT(OUT):: globvar(nx,ny)

  globvar(:,:) = locvar(:,:)

  RETURN
END SUBROUTINE mpimerge2di

SUBROUTINE mpimerge2dx(locvar,nx,ny,yp,globvar,istatus)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny
  REAL,    INTENT(IN) :: locvar(nx,ny)
  REAL,    INTENT(OUT):: globvar(nx,ny)
  INTEGER, INTENT(IN) :: yp
  INTEGER, INTENT(OUT):: istatus

  istatus = 0

  globvar(:,:) = locvar(:,:)

  RETURN
END SUBROUTINE mpimerge2dx

SUBROUTINE mpimerge2dy(locvar,nx,ny,xp,globvar,istatus)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny
  REAL,    INTENT(IN) :: locvar(nx,ny)
  REAL,    INTENT(OUT):: globvar(nx,ny)
  INTEGER, INTENT(IN) :: xp
  INTEGER, INTENT(OUT):: istatus

  istatus = 0

  globvar(:,:) = locvar(:,:)

  RETURN
END SUBROUTINE mpimerge2dy

SUBROUTINE mpimerge3d(locvar,nx,ny,nz,globvar)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL,    INTENT(IN) :: locvar(nx,ny,nz)
  REAL,    INTENT(OUT):: globvar(nx,ny,nz)

  globvar(:,:,:) = locvar(:,:,:)

  RETURN
END SUBROUTINE mpimerge3d

SUBROUTINE mpimerge3di(locvar,nx,ny,nz,globvar)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: locvar(nx,ny,nz)
  INTEGER, INTENT(OUT):: globvar(nx,ny,nz)

  globvar(:,:,:) = locvar(:,:,:)

  RETURN
END SUBROUTINE mpimerge3di

SUBROUTINE mpimerge4d(locvar,nx,ny,nz,nn,globvar)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz,nn
  REAL,    INTENT(IN) :: locvar(nx,ny,nz,nn)
  REAL,    INTENT(OUT):: globvar(nx,ny,nz,nn)

  globvar(:,:,:,:) = locvar(:,:,:,:)

  RETURN
END SUBROUTINE mpimerge4d

SUBROUTINE mpisplit1dx(globvar,nx,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx
  REAL,    INTENT(IN)  :: globvar(nx)
  REAL,    INTENT(OUT) :: var(nx)

  var(:) = globvar(:)

  RETURN
END SUBROUTINE mpisplit1dx

SUBROUTINE mpisplit1dy(globvar,ny,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ny
  REAL,    INTENT(IN)  :: globvar(ny)
  REAL,    INTENT(OUT) :: var(ny)

  var(:) = globvar(:)

  RETURN
END SUBROUTINE mpisplit1dy

SUBROUTINE mpisplit2d(globvar,nx,ny,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny
  REAL,    INTENT(IN)  :: globvar(nx,ny)
  REAL,    INTENT(OUT) :: var(nx,ny)

  var(:,:) = globvar(:,:)

  RETURN
END SUBROUTINE mpisplit2d

SUBROUTINE mpisplit2di(globvar,nx,ny,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny
  INTEGER, INTENT(IN)  :: globvar(nx,ny)
  INTEGER, INTENT(OUT) :: var(nx,ny)

  var(:,:) = globvar(:,:)

  RETURN
END SUBROUTINE mpisplit2di

SUBROUTINE mpisplit3d(globvar,nx,ny,nz,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: globvar(nx,ny,nz)
  REAL,    INTENT(OUT) :: var(nx,ny,nz)

  var(:,:,:) = globvar(:,:,:)

  RETURN
END SUBROUTINE mpisplit3d

SUBROUTINE mpisplit3di(globvar,nx,ny,nz,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny,nz
  INTEGER, INTENT(IN)  :: globvar(nx,ny,nz)
  INTEGER, INTENT(OUT) :: var(nx,ny,nz)

  var(:,:,:) = globvar(:,:,:)

  RETURN
END SUBROUTINE mpisplit3di

SUBROUTINE mpisplit4d(globvar,nx,ny,nz,nn,var)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny,nz,nn
  REAL,    INTENT(IN)  :: globvar(nx,ny,nz,nn)
  REAL,    INTENT(OUT) :: var(nx,ny,nz,nn)

  var(:,:,:,:) = globvar(:,:,:,:)

  RETURN
END SUBROUTINE mpisplit4d

SUBROUTINE mpisplit2dns
  RETURN
END SUBROUTINE mpisplit2dns

SUBROUTINE mpisplit2dew
  RETURN
END SUBROUTINE mpisplit2dew

SUBROUTINE mpsendr
  RETURN
END SUBROUTINE mpsendr

SUBROUTINE mpsendi
  RETURN
END SUBROUTINE mpsendi

SUBROUTINE mprecvr
  RETURN
END SUBROUTINE mprecvr

SUBROUTINE mprecvi
  RETURN
END SUBROUTINE mprecvi

SUBROUTINE mpmaxi
  RETURN
END SUBROUTINE mpmaxi

SUBROUTINE mpmini
  RETURN
END SUBROUTINE mpmini

SUBROUTINE mpbcastr
  RETURN
END SUBROUTINE mpbcastr

SUBROUTINE mpbcastra
  RETURN
END SUBROUTINE mpbcastra

SUBROUTINE mpbcasti
  RETURN
END SUBROUTINE mpbcasti

SUBROUTINE mptotali(var)
  IMPLICIT NONE
  INTEGER :: var
  RETURN
END SUBROUTINE mptotali

SUBROUTINE mpsumi(var,ndim)
  IMPLICIT NONE
  INTEGER :: ndim
  INTEGER :: var(ndim)
  RETURN
END SUBROUTINE mpsumi

SUBROUTINE mpsumr(var,ndim)
  IMPLICIT NONE
  INTEGER :: ndim
  REAL    :: var(ndim)
  RETURN
END SUBROUTINE mpsumr

SUBROUTINE mpsumdp(var,ndim)
  IMPLICIT NONE
  INTEGER :: ndim
  DOUBLE PRECISION    :: var(ndim)
  RETURN
END SUBROUTINE mpsumdp

SUBROUTINE mpminr(var)
  IMPLICIT NONE
  REAL :: var
  RETURN
END SUBROUTINE mpminr

SUBROUTINE mpmaxr(var)
  IMPLICIT NONE
  REAL :: var
  RETURN
END SUBROUTINE mpmaxr

SUBROUTINE mpgatheri(sendbuf,sendcount,recvbuf,recvcount,istatus)

  IMPLICIT NONE

  INTEGER :: sendcount, recvcount
  INTEGER :: sendbuf(sendcount)
  INTEGER :: recvbuf(recvcount)
  INTEGER :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i

  DO i = 1, sendcount
    recvbuf(i) = sendbuf(i)
  END DO

  RETURN
END SUBROUTINE mpgatheri

SUBROUTINE mpscatteri(sendbuf,sendcount,recvbuf,recvcount,istatus)

  IMPLICIT NONE

  INTEGER :: sendcount, recvcount
  INTEGER :: sendbuf(sendcount)
  INTEGER :: recvbuf(recvcount)
  INTEGER :: istatus

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

  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN)  :: isize, osize
  INTEGER, INTENT(IN)  :: vari(isize)
  INTEGER, INTENT(IN)  :: ROOT, numprocs
  INTEGER, INTENT(IN)  :: lsize(numprocs), ldisp(numprocs)
  INTEGER, INTENT(OUT) :: outvari(osize)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  DO i = 1,isize
    outvari(i) = vari(i)
  END DO

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

  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN)  :: isize, osize
  REAL(P), INTENT(IN)  :: varr(isize)
  INTEGER, INTENT(IN)  :: ROOT, numprocs
  INTEGER, INTENT(IN)  :: lsize(numprocs), ldisp(numprocs)
  REAL(P), INTENT(OUT) :: outvarr(osize)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  DO i = 1, isize
    outvarr(i) = varr(i)
  END DO

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE mpminar


!SUBROUTINE print3d(unt,varname,var,nx,ny,nz)
!
!  IMPLICIT NONE
!
!  INTEGER,      INTENT(IN)  :: unt
!  CHARACTER(*), INTENT(IN)  :: varname
!  INTEGER,      INTENT(IN)  :: nx, ny,nz
!  REAL,         INTENT(IN)  :: var(nx,ny,nz)
!
!!-----------------------------------------------------------------------
!!
!! Misc. local variables
!!
!!-----------------------------------------------------------------------
!
!  INCLUDE 'mp.inc'
!
!  INTEGER           :: i,j,k
!  INTEGER           :: nxlg, nylg
!
!  REAL, ALLOCATABLE :: temlg(:,:,:)
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!! Begin of executable code
!!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  nxlg = (nx-3)*nproc_x + 3
!  nylg = (ny-3)*nproc_y + 3
!
!  ALLOCATE(temlg(nxlg,nylg,nz))
!
!  CALL mpimerge3d(var,nx,ny,nz,temlg)
!
!  IF( myproc == 0) THEN
!
!    WRITE(unt,'(3a)') '--- ',varname,' ---'
!
!    DO k = 1,nz
!      DO j = 1,nylg
!        WRITE(UNIT=unt,FMT='(2(a,I3),a)')  ' --- k = ',k,', j = ',j,' --- '
!        WRITE(UNIT=unt,FMT='(5f16.5)') (temlg(i,j,k),i=1,nxlg)
!     END DO
!    END DO
!
!  END IF
!
!  DEALLOCATE(temlg)
!
!  RETURN
!END SUBROUTINE print3d
