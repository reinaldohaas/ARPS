!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_SPLIT2d                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrf_split2d(globvar,nx,ny,fzone,var)

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
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny       ! Number of grid points in x and y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL,    INTENT(IN)  ::                                               &
           globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone)

  REAL,    INTENT(OUT) :: var(nx,ny)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: target
  INTEGER :: ia,ja, ic,jc
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

          target = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny,MPI_REAL,target,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_recv (var,nx*ny,MPI_REAL,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat) 

        END IF

      END IF

      CALL mpbarrier
 
    END DO      
  END DO      

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:) = globvar(1:nx,1:ny)     

  RETURN
END SUBROUTINE wrf_split2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_SPLIT2DI               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrf_split2di(globvar,nx,ny,fzone,var)

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
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny       ! Number of grid points in x and y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  INTEGER, INTENT(IN)  ::                                               &
           globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone)

  INTEGER, INTENT(OUT) :: var(nx,ny)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: target
  INTEGER :: ia,ja, ic,jc
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

          target = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny,MPI_INTEGER,target,               &
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
END SUBROUTINE wrf_split2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_SPLIT3d                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrf_split3d(globvar,nx,ny,nz,fzone,var)

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
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y and z
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL,    INTENT(IN)  ::                                               &
           globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone,nz)

  REAL, INTENT(OUT) :: var(nx,ny,nz)


!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: target
  INTEGER :: ia,ja, ic,jc
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

          target = ic+(jc-1)*nproc_x-1

          CALL mpi_send (var,nx*ny*nz,MPI_REAL,target,               &
                         mytag, MPI_COMM_WORLD,imstat)
        END IF

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_recv (var,nx*ny*nz,MPI_REAL,master,               &
                         mytag,MPI_COMM_WORLD,stat,imstat) 

        END IF

      END IF

      CALL mpbarrier
 
    END DO      
  END DO      

!At the end, make sure processor 0 contains correct varlue (ic=1,jc=1)
  IF(myproc == 0) var(:,:,:) = globvar(1:nx,1:ny,1:nz)     

  RETURN
END SUBROUTINE wrf_split3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPEXT_WRF_U                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpext_wrf_u(var,nx,ny,nz,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    To extend WRF U staggered variables and leave one fake zone for 
!    WRF variables
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  01/31/2005
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

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in 
                                       ! x, y and z directions
  REAL, INTENT(INOUT) :: var(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem(ny,nz,2)   ! Work array.


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

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: j, k

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions, 
!  var_current(nx,:,:) = var_east_neighbor(1,:,:)
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice to west neighbor
  !
  DO k=1,nz
    DO j=1,ny
      tem(j,k,1) = var(1,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),ny*nz,MPI_REAL,dest,  mptag+tag_e,       &
                    tem(:,:,2),ny*nz,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x )  THEN

    DO k=1,nz
      DO j=1,ny
        var(nx,j,k) = tem(j,k,2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE mpext_wrf_u
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MPEXT_WRF_V                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mpext_wrf_v(var,nx,ny,nz,tem)

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
!  02/01/2005
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

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in 
                                       ! x, y and z directions
  REAL, INTENT(INOUT) :: var(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem(nx,nz,2)   ! Work array.

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

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i,k

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of 
  ! the south neighbor
  !
  DO k=1,nz
    DO i=1,nx
      tem(i,k,1) = var(i,1,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem(:,:,1),nx*nz,MPI_REAL,dest,  mptag+tag_n,       &
                    tem(:,:,2),nx*nz,MPI_REAL,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO k=1,nz
      DO i=1,nx
        var(i,ny,k) = tem(i,k,2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE mpext_wrf_v
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXTEND_U                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extend_u(u,nx,ny,nz,uext,tem1,tem2)

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Extended WRF variable staggered at U grid points in south boundary
!    and east boundary to facilitate the computation of its values
!    at V grid points.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx, ny, nz
  REAL,    INTENT(IN)  :: u(nx,ny,nz)
  REAL,    INTENT(OUT) :: uext(1:nx+1,0:ny,1:nz)
  REAL,    INTENT(OUT) :: tem1(nx,nz,2), tem2(ny,nz,2)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i,j,k

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begin of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !
  ! Valid values at U grid points
  !
  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx
        uext(i,j,k) = u(i,j,k)
      END DO
    END DO
  END DO

  !
  ! Top boundary
  !
  DO j = 1,ny-1
    DO i = 1,nx
      uext(i,j,nz) = u(i,j,nz-1)
    END DO
  END DO

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------
  !
  ! Suppose north boundary is always correct, check sub wrf_split3d
  ! Not true for multifile input
  !

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of 
  ! the south neighbor
  !
  DO k=1,nz
    DO i=1,nx
      tem1(i,k,1) = uext(i,1,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem1(:,:,1),nx*nz,MPI_REAL,dest,  mptag+tag_n,       &
                    tem1(:,:,2),nx*nz,MPI_REAL,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO k=1,nz
      DO i=1,nx
        uext(i,ny,k) = tem1(i,k,2)
      END DO
    END DO
  ELSE
    DO k=1,nz
      DO i=1,nx
        uext(i,ny,k) = uext(i,ny-1,k)
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
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF
  
  ! 
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of 
  ! the north neighbor
  !
  DO k=1,nz
    DO i=1,nx
      tem1(i,k,1) = uext(i,ny-1,k)
    END DO
  END DO
 
  CALL mpi_sendrecv(tem1(:,:,1),nx*nz,MPI_REAL,dest,  mptag+tag_s,       &
                    tem1(:,:,2),nx*nz,MPI_REAL,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)
 
  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 )  THEN
    DO k=1,nz
      DO i=1,nx
        uext(i,0,k) = tem1(i,k,2)
      END DO
    END DO

  ELSE

    DO k = 1,nz
      DO i = 1,nx
        uext(i,0,k) = uext(i,1,k)
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
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of 
  ! the west neighbor
  !
  DO k=1,nz
    DO j=1,ny
      tem2(j,k,1) = uext(2,j,k)
    END DO
  END DO
 
  CALL mpi_sendrecv(tem2(:,:,1),ny*nz,MPI_REAL,dest,  mptag+tag_e,       &
                    tem2(:,:,2),ny*nz,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)
  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x )  THEN

    DO k=1,nz
      DO j=1,ny
        uext(nx+1,j,k) = tem2(j,k,2)
      END DO
    END DO
    uext(nx+1,0,:) = tem2(1,:,2)

  ELSE

    DO k=1,nz
      DO j=1,ny
        uext(nx+1,j,k) = uext(nx,j,k)
      END DO
    END DO
    uext(nx+1,0,:) = uext(nx,1,:)

  END IF

  RETURN
END SUBROUTINE extend_u
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXTEND_V                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extend_v(v,nx,ny,nz,vext,tem1,tem2)

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Extended WRF variable staggered at V grid points in north boundary
!    and west boundary to facilitate the computation of its values
!    at U grid points.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx, ny, nz
  REAL,    INTENT(IN)  :: v(nx,ny,nz)
  REAL,    INTENT(OUT) :: vext(0:nx,1:ny+1,1:nz)
  REAL,    INTENT(OUT) :: tem1(nx,nz,2), tem2(ny,nz,2)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: i,j,k

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_n for north boundary
                             ! mptag + tag_s for south boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begin of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !
  ! Valid values at V grid points
  !
  DO k = 1,nz-1
    DO j = 1,ny
      DO i = 1,nx-1
        vext(i,j,k) = v(i,j,k)
      END DO
    END DO
  END DO

  !
  ! Top boundary
  !
  DO j = 1,ny
    DO i = 1,nx-1
      vext(i,j,nz) = v(i,j,nz-1)
    END DO
  END DO

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
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of 
  ! the east neighbor
  !
  DO k=1,nz
    DO j=1,ny
      tem2(j,k,1) = vext(nx-1,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem2(:,:,1),ny*nz,MPI_REAL,dest,  mptag+tag_w,       &
                    tem2(:,:,2),ny*nz,MPI_REAL,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 ) THEN
    DO k=1,nz
      DO j=1,ny
        vext(0,j,k) = tem2(j,k,2)
      END DO
    END DO
  ELSE
    DO k = 1,nz
      DO j = 1,ny
        vext(0,j,k) = vext(1,j,k)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------
  !
  ! Suppose east boundary was corretly set by wrf_split3d
  ! otherwise messsage passing is needed.
  !

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of 
  ! the west neighbor
  !
  DO k=1,nz
    DO j=1,ny
      tem2(j,k,1) = vext(1,j,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem2(:,:,1),ny*nz,MPI_REAL,dest,  mptag+tag_e,       &
                    tem2(:,:,2),ny*nz,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x )  THEN
    DO k=1,nz
      DO j=1,ny
        vext(nx,j,k) = tem2(j,k,2)
      END DO
    END DO
  ELSE
    DO k=1,nz
      DO j=1,ny
        vext(nx,j,k) = vext(nx-1,j,k)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of 
  ! the south neighbor
  !
  DO k=1,nz
    DO i=1,nx
      tem1(i,k,1) = vext(i,2,k)
    END DO
  END DO

  CALL mpi_sendrecv(tem1(:,:,1),nx*nz,MPI_REAL,dest,  mptag+tag_n,      &
                    tem1(:,:,2),nx*nz,MPI_REAL,source,mptag+tag_n,      &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO k=1,nz
      DO i=1,nx
        vext(i,ny+1,k) = tem1(i,k,2)
      END DO
    END DO
    vext(0,ny+1,:) = tem1(1,:,2)
  ELSE
    DO k = 1,nz
      DO i = 1,nx
        vext(i,ny+1,k) = vext(i,ny,k)
      END DO
    END DO
    vext(0,ny+1,:) = vext(1,ny,:)   ! this point is not used
  END IF


  RETURN
END SUBROUTINE extend_v

!SUBROUTINE printwrf3d(unt,varname,var,stagger,nx,ny,nz)
!  
!  IMPLICIT NONE
!
!  INTEGER,      INTENT(IN)  :: unt
!  CHARACTER(*), INTENT(IN)  :: varname
!  CHARACTER(*), INTENT(IN)  :: stagger
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
!!    WRITE(unt,'(3a)') '--- ',varname,' ---'
!!
!!    DO k = 1,nz
!!      DO j = 1,ny
!!        WRITE(UNIT=unt,FMT='(2I3,a)')  k,j,' --- '
!!        WRITE(UNIT=unt,FMT='(5f16.5)') (var(i,j,k),i=1,nx)
!!      END DO
!!    END DO
!
!  nxlg = (nx-1)*nproc_x + 1
!  nylg = (ny-1)*nproc_y + 1
!
!
!  IF (stagger == 'X') THEN
!     ALLOCATE(temlg(nxlg,nylg-1,nz-1))
!     CALL wrf_merge3du(var,nx,ny,nz,1,temlg)
!  ELSE IF (stagger == 'Y') THEN
!     ALLOCATE(temlg(nxlg-1,nylg,nz-1))
!     CALL wrf_merge3dv(var,nx,ny,nz,1,temlg)
!  ELSE IF (stagger == 'Z') THEN
!     ALLOCATE(temlg(nxlg-1,nylg-1,nz))
!     CALL wrf_merge3dw(var,nx,ny,nz,1,temlg)
!  ELSE
!     ALLOCATE(temlg(nxlg-1,nylg-1,nz-1))
!     CALL wrf_merge3dt(var,nx,ny,nz,1,temlg)
!  END IF
!
!  IF( myproc == 0) THEN
!
!    WRITE(unt,'(3a)') '--- ',varname,' ---'
!
!    IF (stagger == 'X') THEN
!      DO k = 1,nz-1
!        DO j = 1,nylg-1
!          WRITE(UNIT=unt,FMT='(2I3,a)')  k,j,' --- '
!          WRITE(UNIT=unt,FMT='(5f16.5)') (temlg(i,j,k),i=1,nxlg)
!       END DO
!      END DO
!    ELSE IF (stagger == 'Y') THEN
!      DO k = 1,nz-1
!        DO j = 1,nylg
!          WRITE(UNIT=unt,FMT='(2I3,a)')  k,j,' --- '
!          WRITE(UNIT=unt,FMT='(5f16.5)') (temlg(i,j,k),i=1,nxlg-1)
!       END DO
!      END DO
!    ELSE IF (stagger == 'Z') THEN
!      DO k = 1,nz
!        DO j = 1,nylg-1
!          WRITE(UNIT=unt,FMT='(2I3,a)')  k,j,' --- '
!          WRITE(UNIT=unt,FMT='(5f16.5)') (temlg(i,j,k),i=1,nxlg-1)
!       END DO
!      END DO
!    ELSE
!      DO k = 1,nz-1
!        DO j = 1,nylg-1
!          WRITE(UNIT=unt,FMT='(2I3,a)')  k,j,' --- '
!          WRITE(UNIT=unt,FMT='(5f16.5)') (temlg(i,j,k),i=1,nxlg-1)
!       END DO
!      END DO
!    END IF
!
!  END IF
!
!  DEALLOCATE(temlg)
!
!  RETURN
!END SUBROUTINE printwrf3d
!!
!!##################################################################
!!##################################################################
!!######                                                      ######
!!######                SUBROUTINE WRF_merge3dt               ######
!!######                                                      ######
!!######                     Developed by                     ######
!!######     Center for Analysis and Prediction of Storms     ######
!!######                University of Oklahoma                ######
!!######                                                      ######
!!##################################################################
!!##################################################################
!!
!SUBROUTINE wrf_merge3dt(locvar,nx,ny,nz,fzone,globvar)
!!
!!-----------------------------------------------------------------------
!!
!!  PURPOSE:
!!
!!  Generate global array from a multiprocessor run for variable at 
!!  mass grid points  
!!
!!-----------------------------------------------------------------------
!!
!!
!!  AUTHOR: Yunheng Wang
!!  2004/09/29
!!
!!  MODIFICATION HISTORY:
!!
!!-----------------------------------------------------------------------
!!
!!  INPUT :
!!
!!    localvar Variable to be written.
!!
!!    nx       Number of grid points in the x-direction (east/west)
!!    ny       Number of grid points in the y-direction (north/south)
!!    nz       Number of grid points in the vertical
!!
!!  OUTPUT:
!!
!!    globvar  global variable to be output
!!
!!-----------------------------------------------------------------------
!!
!!  Variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  IMPLICIT NONE
!
!!-----------------------------------------------------------------------
!!
!!  Include files.
!!
!!-----------------------------------------------------------------------
!
!  INCLUDE 'mpif.h'
!  INCLUDE 'mp.inc'
!
!  INTEGER, INTENT(IN) :: nx,ny,nz          
!                            ! Number of stagger grid points in x, y and z
!  INTEGER, INTENT(IN) :: fzone
!  REAL,    INTENT(IN) :: locvar(nx,ny,nz)
!
!  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
!     globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone-1,nz-1)
!!-----------------------------------------------------------------------
!!
!!  Local variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  INTEGER, PARAMETER :: master = 0
!  REAL    :: tem(nx,ny,nz)
!
!  INTEGER :: mptag        ! Unique MPI id.
!  INTEGER :: mytag
!  INTEGER :: source
!  INTEGER :: ia,ja, ic,jc
!  INTEGER :: stat(mpi_status_size)
!  INTEGER :: imstat
!  INTEGER :: i, j, k
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!  Beginning of executable code...
!!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL inctag
!  mptag = gentag
!
!!-----------------------------------------------------------------------
!!
!!  fill the globvar array
!!
!!-----------------------------------------------------------------------
!
!  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
! 
!  DO jc=1,nproc_y
!    DO ic=1,nproc_x
!
!      ! message passing section...
!
!      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0
!
!        mytag = mptag + 100 + ic + jc
!
!        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
!       
!          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
!                         mytag,MPI_COMM_WORLD,imstat) 
!
!        END IF
!
!        IF(myproc == 0)THEN          ! receive data
!
!          source = ic+(jc-1)*nproc_x-1
!
!          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
!                         mytag, MPI_COMM_WORLD,stat,imstat)
!        END IF
!
!      END IF
!
!      !  storage section
!
!      IF(myproc == 0)THEN  ! store data into globvar
!
!        DO k=1,nz-1
!          DO j=1,ny-1
!            ja = j + (jc-1)*(ny-fzone)
!            DO i=1,nx-1
!              ia = i + (ic-1)*(nx-fzone)
!              globvar(ia,ja,k) = tem(i,j,k)
!            END DO          
!          END DO      
!        END DO      
!
!      END IF
!
!      CALL mpbarrier
!
!    END DO      
!  END DO      
!
!  RETURN
!END SUBROUTINE wrf_merge3dt
!!
!!##################################################################
!!##################################################################
!!######                                                      ######
!!######                SUBROUTINE WRF_merge3du               ######
!!######                                                      ######
!!######                     Developed by                     ######
!!######     Center for Analysis and Prediction of Storms     ######
!!######                University of Oklahoma                ######
!!######                                                      ######
!!##################################################################
!!##################################################################
!!
!SUBROUTINE wrf_merge3du(locvar,nx,ny,nz,fzone,globvar)
!!
!!-----------------------------------------------------------------------
!!
!!  PURPOSE:
!!
!!  Generate global array from a multiprocessor run for variable at 
!!  U grid points  
!!
!!-----------------------------------------------------------------------
!!
!!
!!  AUTHOR: Yunheng Wang
!!  2004/09/29
!!
!!  MODIFICATION HISTORY:
!!
!!-----------------------------------------------------------------------
!!
!!  INPUT :
!!
!!    localvar Variable to be written.
!!
!!    nx       Number of grid points in the x-direction (east/west)
!!    ny       Number of grid points in the y-direction (north/south)
!!    nz       Number of grid points in the vertical
!!
!!  OUTPUT:
!!
!!    globvar  global variable to be output
!!
!!-----------------------------------------------------------------------
!!
!!  Variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  IMPLICIT NONE
!
!!-----------------------------------------------------------------------
!!
!!  Include files.
!!
!!-----------------------------------------------------------------------
!
!  INCLUDE 'mpif.h'
!  INCLUDE 'mp.inc'
!
!  INTEGER, INTENT(IN) :: nx,ny,nz          
!                            ! Number of stagger grid points in x, y and z
!  INTEGER, INTENT(IN) :: fzone
!  REAL,    INTENT(IN) :: locvar(nx,ny,nz)
!
!  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
!     globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone-1,nz-1)
!!-----------------------------------------------------------------------
!!
!!  Local variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  INTEGER, PARAMETER :: master = 0
!  REAL    :: tem(nx,ny,nz)
!
!  INTEGER :: mptag        ! Unique MPI id.
!  INTEGER :: mytag
!  INTEGER :: source
!  INTEGER :: ia,ja, ic,jc
!  INTEGER :: stat(mpi_status_size)
!  INTEGER :: imstat
!  INTEGER :: i, j, k
!  INTEGER :: i0
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!  Beginning of executable code...
!!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL inctag
!  mptag = gentag
!
!!-----------------------------------------------------------------------
!!
!!  fill the globvar array
!!
!!-----------------------------------------------------------------------
!
!  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
! 
!  DO jc=1,nproc_y
!    DO ic=1,nproc_x
!
!      ! message passing section...
!
!      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0
!
!        mytag = mptag + 100 + ic + jc
!
!        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
!       
!          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
!                         mytag,MPI_COMM_WORLD,imstat) 
!
!        END IF
!
!        IF(myproc == 0)THEN          ! receive data
!
!          source = ic+(jc-1)*nproc_x-1
!
!          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
!                         mytag, MPI_COMM_WORLD,stat,imstat)
!        END IF
!
!      END IF
!
!      !  storage section
!
!      IF(myproc == 0)THEN  ! store data into globvar
!
!        i0 = 1
!        IF (ic > 1) i0 = 2
!        DO k=1,nz-1
!          DO j=1,ny-1
!            ja = j + (jc-1)*(ny-fzone)
!            DO i=i0,nx
!              ia = i + (ic-1)*(nx-fzone)
!              globvar(ia,ja,k) = tem(i,j,k)
!            END DO          
!          END DO      
!        END DO      
!
!      END IF
!
!      CALL mpbarrier
!
!    END DO      
!  END DO      
!
!  RETURN
!END SUBROUTINE wrf_merge3du
!!
!!##################################################################
!!##################################################################
!!######                                                      ######
!!######                SUBROUTINE WRF_merge3dv               ######
!!######                                                      ######
!!######                     Developed by                     ######
!!######     Center for Analysis and Prediction of Storms     ######
!!######                University of Oklahoma                ######
!!######                                                      ######
!!##################################################################
!!##################################################################
!!
!SUBROUTINE wrf_merge3dv(locvar,nx,ny,nz,fzone,globvar)
!!
!!-----------------------------------------------------------------------
!!
!!  PURPOSE:
!!
!!  Generate global array from a multiprocessor run for variable at 
!!  mass grid points  
!!
!!-----------------------------------------------------------------------
!!
!!
!!  AUTHOR: Yunheng Wang
!!  2004/09/29
!!
!!  MODIFICATION HISTORY:
!!
!!-----------------------------------------------------------------------
!!
!!  INPUT :
!!
!!    localvar Variable to be written.
!!
!!    nx       Number of grid points in the x-direction (east/west)
!!    ny       Number of grid points in the y-direction (north/south)
!!    nz       Number of grid points in the vertical
!!
!!  OUTPUT:
!!
!!    globvar  global variable to be output
!!
!!-----------------------------------------------------------------------
!!
!!  Variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  IMPLICIT NONE
!
!!-----------------------------------------------------------------------
!!
!!  Include files.
!!
!!-----------------------------------------------------------------------
!
!  INCLUDE 'mpif.h'
!  INCLUDE 'mp.inc'
!
!  INTEGER, INTENT(IN) :: nx,ny,nz          
!                            ! Number of stagger grid points in x, y and z
!  INTEGER, INTENT(IN) :: fzone
!  REAL,    INTENT(IN) :: locvar(nx,ny,nz)
!
!  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
!     globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone,nz-1)
!!-----------------------------------------------------------------------
!!
!!  Local variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  INTEGER, PARAMETER :: master = 0
!  REAL    :: tem(nx,ny,nz)
!
!  INTEGER :: mptag        ! Unique MPI id.
!  INTEGER :: mytag
!  INTEGER :: source
!  INTEGER :: ia,ja, ic,jc
!  INTEGER :: stat(mpi_status_size)
!  INTEGER :: imstat
!  INTEGER :: i, j, k
!  INTEGER :: j0
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!  Beginning of executable code...
!!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL inctag
!  mptag = gentag
!
!!-----------------------------------------------------------------------
!!
!!  fill the globvar array
!!
!!-----------------------------------------------------------------------
!
!  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
! 
!  DO jc=1,nproc_y
!    DO ic=1,nproc_x
!
!      ! message passing section...
!
!      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0
!
!        mytag = mptag + 100 + ic + jc
!
!        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
!       
!          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
!                         mytag,MPI_COMM_WORLD,imstat) 
!
!        END IF
!
!        IF(myproc == 0)THEN          ! receive data
!
!          source = ic+(jc-1)*nproc_x-1
!
!          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
!                         mytag, MPI_COMM_WORLD,stat,imstat)
!        END IF
!
!      END IF
!
!      !  storage section
!
!      IF(myproc == 0)THEN  ! store data into globvar
!
!        j0 = 1
!        IF (jc > 1) j0 = 2
!
!        DO k=1,nz-1
!          DO j=j0,ny
!            ja = j + (jc-1)*(ny-fzone)
!            DO i=1,nx-1
!              ia = i + (ic-1)*(nx-fzone)
!              globvar(ia,ja,k) = tem(i,j,k)
!            END DO          
!          END DO      
!        END DO      
!
!      END IF
!
!      CALL mpbarrier
!
!    END DO      
!  END DO      
!
!  RETURN
!END SUBROUTINE wrf_merge3dv
!!
!!##################################################################
!!##################################################################
!!######                                                      ######
!!######                SUBROUTINE WRF_merge3dw               ######
!!######                                                      ######
!!######                     Developed by                     ######
!!######     Center for Analysis and Prediction of Storms     ######
!!######                University of Oklahoma                ######
!!######                                                      ######
!!##################################################################
!!##################################################################
!!
!SUBROUTINE wrf_merge3dw(locvar,nx,ny,nz,fzone,globvar)
!!
!!-----------------------------------------------------------------------
!!
!!  PURPOSE:
!!
!!  Generate global array from a multiprocessor run for variable at 
!!  W grid points  
!!
!!-----------------------------------------------------------------------
!!
!!
!!  AUTHOR: Yunheng Wang
!!  2004/09/29
!!
!!  MODIFICATION HISTORY:
!!
!!-----------------------------------------------------------------------
!!
!!  INPUT :
!!
!!    localvar Variable to be written.
!!
!!    nx       Number of grid points in the x-direction (east/west)
!!    ny       Number of grid points in the y-direction (north/south)
!!    nz       Number of grid points in the vertical
!!
!!  OUTPUT:
!!
!!    globvar  global variable to be output
!!
!!-----------------------------------------------------------------------
!!
!!  Variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  IMPLICIT NONE
!
!!-----------------------------------------------------------------------
!!
!!  Include files.
!!
!!-----------------------------------------------------------------------
!
!  INCLUDE 'mpif.h'
!  INCLUDE 'mp.inc'
!
!  INTEGER, INTENT(IN) :: nx,ny,nz          
!                            ! Number of stagger grid points in x, y and z
!  INTEGER, INTENT(IN) :: fzone
!  REAL,    INTENT(IN) :: locvar(nx,ny,nz)
!
!  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
!     globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone-1,nz)
!!-----------------------------------------------------------------------
!!
!!  Local variable declarations.
!!
!!-----------------------------------------------------------------------
!
!  INTEGER, PARAMETER :: master = 0
!  REAL    :: tem(nx,ny,nz)
!
!  INTEGER :: mptag        ! Unique MPI id.
!  INTEGER :: mytag
!  INTEGER :: source
!  INTEGER :: ia,ja, ic,jc
!  INTEGER :: stat(mpi_status_size)
!  INTEGER :: imstat
!  INTEGER :: i, j, k
!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!  Beginning of executable code...
!!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL inctag
!  mptag = gentag
!
!!-----------------------------------------------------------------------
!!
!!  fill the globvar array
!!
!!-----------------------------------------------------------------------
!
!  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
! 
!  DO jc=1,nproc_y
!    DO ic=1,nproc_x
!
!      ! message passing section...
!
!      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0
!
!        mytag = mptag + 100 + ic + jc
!
!        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
!       
!          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
!                         mytag,MPI_COMM_WORLD,imstat) 
!
!        END IF
!
!        IF(myproc == 0)THEN          ! receive data
!
!          source = ic+(jc-1)*nproc_x-1
!
!          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
!                         mytag, MPI_COMM_WORLD,stat,imstat)
!        END IF
!
!      END IF
!
!      !  storage section
!
!      IF(myproc == 0)THEN  ! store data into globvar
!
!        DO k=1,nz
!          DO j=1,ny-1
!            ja = j + (jc-1)*(ny-fzone)
!            DO i=1,nx-1
!              ia = i + (ic-1)*(nx-fzone)
!              globvar(ia,ja,k) = tem(i,j,k)
!            END DO          
!          END DO      
!        END DO      
!
!      END IF
!
!      CALL mpbarrier
!
!    END DO      
!  END DO      
!
!  RETURN
!END SUBROUTINE wrf_merge3dw
