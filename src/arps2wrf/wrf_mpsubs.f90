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
!######                SUBROUTINE WRF_merge3dt               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3dt(locvar,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  mass grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/29
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
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz          
                            ! Number of stagger grid points in x, y and z
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny,nz)

  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
     globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone-1,nz-1)
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny,nz)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
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

  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        DO k=1,nz-1
          DO j=1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx-1
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja,k) = tem(i,j,k)
            END DO          
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge3dt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge3du               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3du(locvar,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  U grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/29
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
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz          
                            ! Number of stagger grid points in x, y and z
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny,nz)

  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
     globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone-1,nz-1)
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny,nz)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k
  INTEGER :: i0

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

  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        i0 = 1
        IF (ic > 1) i0 = 2
        DO k=1,nz-1
          DO j=1,ny-1
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

  RETURN
END SUBROUTINE wrf_merge3du
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge3dv               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3dv(locvar,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  mass grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/29
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
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz          
                            ! Number of stagger grid points in x, y and z
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny,nz)

  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
     globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone,nz-1)
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny,nz)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k
  INTEGER :: j0

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

  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        j0 = 1
        IF (jc > 1) j0 = 2

        DO k=1,nz-1
          DO j=j0,ny
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx-1
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja,k) = tem(i,j,k)
            END DO          
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge3dv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge3dw               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3dw(locvar,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  W grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/29
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
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz          
                            ! Number of stagger grid points in x, y and z
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny,nz)

  REAL,    INTENT(OUT)::  & ! Output array in global domain, defined on mass points
     globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone-1,nz)
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny,nz)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
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

  tem(:,:,:) = locvar(:,:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny*nz,MPI_REAL,master,               &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny*nz,MPI_REAL,source,               &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        DO k=1,nz
          DO j=1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            DO i=1,nx-1
              ia = i + (ic-1)*(nx-fzone)
              globvar(ia,ja,k) = tem(i,j,k)
            END DO          
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge3dw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge2dt               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2dt(locvar,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  mass grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/04
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
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny ! Number of stagger grid points in x, y
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny)

  REAL,    INTENT(OUT):: globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone-1)
                    ! Output array in global domain, defined on mass points
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
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

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny,MPI_REAL,master,                  &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny,MPI_REAL,source,                     &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        DO j=1,ny-1
          ja = j + (jc-1)*(ny-fzone)
          DO i=1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globvar(ia,ja) = tem(i,j)
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge2dt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge2du               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2du(locvar,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  U grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/04
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

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: nx,ny ! Number of stagger grid points in x, y
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny)

  REAL,    INTENT(OUT):: globvar((nx-fzone)*nproc_x+fzone,(ny-fzone)*nproc_y+fzone-1)
                    ! Output array in global domain, defined on U points
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
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

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny,MPI_REAL,master,                  &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny,MPI_REAL,source,                     &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        DO j=1,ny-1
          ja = j + (jc-1)*(ny-fzone)
          DO i=1,nx
            ia = i + (ic-1)*(nx-fzone)
            globvar(ia,ja) = tem(i,j)
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge2du
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge2dv               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2dv(locvar,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  V grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/04
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

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'mpif.h'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: nx,ny ! Number of stagger grid points in x, y
  INTEGER, INTENT(IN) :: fzone
  REAL,    INTENT(IN) :: locvar(nx,ny)

  REAL,    INTENT(OUT):: globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone)
                    ! Output array in global domain, defined on mass points
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  REAL    :: tem(nx,ny)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
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

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny,MPI_REAL,master,                  &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny,MPI_REAL,source,                     &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        DO j=1,ny
          ja = j + (jc-1)*(ny-fzone)
          DO i=1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globvar(ia,ja) = tem(i,j)
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge2dv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_merge2di               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2di(locvar,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate global array from a multiprocessor run for variable at 
!  mass grid points  
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/04
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
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny ! Number of stagger grid points in x, y
  INTEGER, INTENT(IN) :: fzone
  INTEGER, INTENT(IN) :: locvar(nx,ny)

  INTEGER, INTENT(OUT):: globvar((nx-fzone)*nproc_x+fzone-1,(ny-fzone)*nproc_y+fzone-1)
                    ! Output array in global domain, defined on mass points
!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0
  INTEGER :: tem(nx,ny)

  INTEGER :: mptag        ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
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

  tem(:,:) = locvar(:,:) ! each processor stores the locvar into tem
 
  DO jc=1,nproc_y
    DO ic=1,nproc_x

      ! message passing section...

      IF(ic /=1 .OR. jc /=1)THEN     !  pass data to processor 0

        mytag = mptag + 100 + ic + jc

        IF(myproc == (ic+(jc-1)*nproc_x-1))THEN
       
          CALL mpi_send (locvar,nx*ny,MPI_INTEGER,master,              &
                         mytag,MPI_COMM_WORLD,imstat) 

        END IF

        IF(myproc == 0)THEN          ! receive data

          source = ic+(jc-1)*nproc_x-1

          CALL mpi_recv (tem,nx*ny,MPI_INTEGER,source,                 &
                         mytag, MPI_COMM_WORLD,stat,imstat)
        END IF

      END IF

      !  storage section

      IF(myproc == 0)THEN  ! store data into globvar

        DO j=1,ny-1
          ja = j + (jc-1)*(ny-fzone)
          DO i=1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globvar(ia,ja) = tem(i,j)
          END DO      
        END DO      

      END IF

      CALL mpbarrier

    END DO      
  END DO      

  RETURN
END SUBROUTINE wrf_merge2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_mergebdyu              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mergebdyu(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdyzone,          &
                         fzone,globwt,globet,globsu,globnu,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Merge lateral boundary arrays to U staggered grids.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/08
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    bdyw     West boundary array without stagger (local)
!    bdye     East boundary array without stagger (local)
!    bdys     South boundary array without stagger (local)
!    bdyn     North boundary array without stagger (local)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are local staggered size, i.e. scalar 
!                   grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globwt   West lateral boundary array staggered at global Mass points
!    globet   East lateral boundary array staggered at global Mass points
!    globsu   South lateral boundary array staggered at global U points
!    globnu   North lateral boundary array staggered at global U points
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: bdyzone
  INTEGER, INTENT(IN) :: fzone

  REAL,   INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL,   INTENT(OUT) :: globwt((ny-fzone)*nproc_y+fzone-1,nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globet((ny-fzone)*nproc_y+fzone-1,nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globsu((nx-fzone)*nproc_x+fzone,  nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globnu((nx-fzone)*nproc_x+fzone,  nz-1,bdyzone)
                    ! Output array in global domain, defined on U points

  REAL,   INTENT(OUT) :: tem1(nx,nz,bdyzone)
  REAL,   INTENT(OUT) :: tem2(ny,nz,bdyzone)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0  ! NetCDF is not parallel I/O, So 
                         ! only one processor can access the opened file
  INTEGER :: mptag       ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, i0
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  West/south boundary, no-stagger, processor 0 do merge
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdyw(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdys(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! West boudnary
  !
  ic = 1
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 100 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO j = 1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            globwt(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! South boudnary
  !
  jc = 1
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 200 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      i0 = 1
      IF (ic > 1) i0 = 2

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO i = i0,nx
            ia = i + (ic-1)*(nx-fzone)
            globsu(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

!-----------------------------------------------------------------------
!
!  East/North boundary, U-stagger
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdye(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdyn(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! East boudnary
  !
  ic = nproc_x
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 300 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO j = 1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            globet(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! North boudnary
  !
  jc = nproc_y
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 400 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      i0 = 1
      IF (ic > 1) i0 = 2

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO i = i0,nx
            ia = i + (ic-1)*(nx-fzone)
            globnu(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

  RETURN
END SUBROUTINE wrf_mergebdyu
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_mergebdyv              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mergebdyv(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdyzone,          &
                         fzone,globwv,globev,globst,globnt,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Merge lateral boundary arrays to V staggered grids.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/08
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    bdyw     West boundary array without stagger (local)
!    bdye     East boundary array without stagger (local)
!    bdys     South boundary array without stagger (local)
!    bdyn     North boundary array without stagger (local)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are local staggered size, i.e. scalar 
!                   grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globwv   West lateral boundary array staggered at global V points
!    globev   East lateral boundary array staggered at global V points
!    globst   South lateral boundary array staggered at global mass points
!    globnt   North lateral boundary array staggered at global mass points
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: bdyzone
  INTEGER, INTENT(IN) :: fzone

  REAL,   INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL,   INTENT(OUT) :: globwv((ny-fzone)*nproc_y+fzone,  nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globev((ny-fzone)*nproc_y+fzone,  nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globst((nx-fzone)*nproc_x+fzone-1,nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globnt((nx-fzone)*nproc_x+fzone-1,nz-1,bdyzone)
                    ! Output array in global domain, defined on U points

  REAL,   INTENT(OUT) :: tem1(nx,nz,bdyzone)
  REAL,   INTENT(OUT) :: tem2(ny,nz,bdyzone)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0  ! NetCDF is not parallel I/O, So 
                         ! only one processor can access the opened file
  INTEGER :: mptag       ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc, j0
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  West/south boundary 
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdyw(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdys(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! West boudnary, V- staggered
  !
  ic = 1
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 100 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      j0 = 1
      IF (jc > 1) j0 = 2
      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO j = j0,ny
            ja = j + (jc-1)*(ny-fzone)
            globwv(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! South boudnary, mass grid
  !
  jc = 1
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 200 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO i = 1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globst(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

!-----------------------------------------------------------------------
!
!  East/North boundary, U-stagger
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdye(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdyn(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! East boudnary
  !
  ic = nproc_x
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 300 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      j0 = 1
      IF (jc > 1) j0 = 2
      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO j = j0,ny
            ja = j + (jc-1)*(ny-fzone)
            globev(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! North boudnary
  !
  jc = nproc_y
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 400 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO i = 1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globnt(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

  RETURN
END SUBROUTINE wrf_mergebdyv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_mergebdyt              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mergebdyt(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdyzone,          &
                         fzone,globwt,globet,globst,globnt,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Merge lateral boundary arrays to mass grids.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/08
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    bdyw     West boundary array without stagger (local)
!    bdye     East boundary array without stagger (local)
!    bdys     South boundary array without stagger (local)
!    bdyn     North boundary array without stagger (local)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are local staggered size, i.e. scalar 
!                   grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globwt   West lateral boundary array staggered at global Mass points
!    globet   East lateral boundary array staggered at global Mass points
!    globst   South lateral boundary array staggered at global Mass points
!    globnt   North lateral boundary array staggered at global Mass points
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files.
!
!-----------------------------------------------------------------------

  INCLUDE   'mpif.h'
  INCLUDE   'mp.inc'

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: bdyzone
  INTEGER, INTENT(IN) :: fzone

  REAL,   INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL,   INTENT(OUT) :: globwt((ny-fzone)*nproc_y+fzone-1,nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globet((ny-fzone)*nproc_y+fzone-1,nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globst((nx-fzone)*nproc_x+fzone-1,nz-1,bdyzone)
  REAL,   INTENT(OUT) :: globnt((nx-fzone)*nproc_x+fzone-1,nz-1,bdyzone)
                    ! Output array in global domain, defined on U points

  REAL,   INTENT(OUT) :: tem1(nx,nz,bdyzone)
  REAL,   INTENT(OUT) :: tem2(ny,nz,bdyzone)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0  ! NetCDF is not parallel I/O, So 
                         ! only one processor can access the opened file
  INTEGER :: mptag       ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  West/south boundary, no-stagger, processor 0 do merge
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdyw(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdys(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! West boudnary
  !
  ic = 1
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 100 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO j = 1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            globwt(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! South boudnary
  !
  jc = 1
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 200 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO i = 1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globst(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

!-----------------------------------------------------------------------
!
!  East/North boundary, U-stagger
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdye(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdyn(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! East boudnary
  !
  ic = nproc_x
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 300 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO j = 1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            globet(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! North boudnary
  !
  jc = nproc_y
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 400 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz-1
          DO i = 1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globnt(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

  RETURN
END SUBROUTINE wrf_mergebdyt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_mergebdyw              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mergebdyw(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdyzone,          &
                         fzone,globww,globew,globsw,globnw,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Merge lateral boundary arrays to W staggered grids.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/08
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    bdyw     West boundary array without stagger (local)
!    bdye     East boundary array without stagger (local)
!    bdys     South boundary array without stagger (local)
!    bdyn     North boundary array without stagger (local)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are local staggered size, i.e. scalar 
!                   grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globww   West lateral boundary array staggered at global W points
!    globew   East lateral boundary array staggered at global W points
!    globsw   South lateral boundary array staggered at global W points
!    globnw   North lateral boundary array staggered at global W points
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: bdyzone
  INTEGER, INTENT(IN) :: fzone

  REAL,   INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL,   INTENT(OUT) :: globww((ny-fzone)*nproc_y+fzone-1,nz,bdyzone)
  REAL,   INTENT(OUT) :: globew((ny-fzone)*nproc_y+fzone-1,nz,bdyzone)
  REAL,   INTENT(OUT) :: globsw((nx-fzone)*nproc_x+fzone-1,nz,bdyzone)
  REAL,   INTENT(OUT) :: globnw((nx-fzone)*nproc_x+fzone-1,nz,bdyzone)
                    ! Output array in global domain, defined on U points

  REAL,   INTENT(OUT) :: tem1(nx,nz,bdyzone)
  REAL,   INTENT(OUT) :: tem2(ny,nz,bdyzone)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0  ! NetCDF is not parallel I/O, So 
                         ! only one processor can access the opened file
  INTEGER :: mptag       ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  West/south boundary, no-stagger, processor 0 do merge
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdyw(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdys(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! West boudnary
  !
  ic = 1
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 100 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz
          DO j = 1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            globww(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! South boudnary
  !
  jc = 1
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 200 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz
          DO i = 1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globsw(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

!-----------------------------------------------------------------------
!
!  East/North boundary, U-stagger
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO j = 1,ny
        tem2(j,k,bdy) = bdye(j,k,bdy)
      END DO

      DO i = 1,nx
        tem1(i,k,bdy) = bdyn(i,k,bdy)
      END DO
    END DO
  END DO

  !
  ! East boudnary
  !
  ic = nproc_x
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 300 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz
          DO j = 1,ny-1
            ja = j + (jc-1)*(ny-fzone)
            globew(ja,k,bdy) = tem2(j,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! North boudnary
  !
  jc = nproc_y
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 400 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*nz*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*nz*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO k = 1,nz
          DO i = 1,nx-1
            ia = i + (ic-1)*(nx-fzone)
            globnw(ia,k,bdy) = tem1(i,k,bdy)
          END DO      
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

  RETURN
END SUBROUTINE wrf_mergebdyw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_mergebdy2d             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mergebdy2d(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdyzone,          &
                         fzone,globw2d,globe2d,globs2d,globn2d,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Merge lateral boundary arrays to scalar staggered grids for 2D arrays.
!
!  NOTE: the input unstaggered arrays have already been packed to be
!        3D arrays to be compatible with WRF version 1.3. This may be
!        changed later if WRFV1.3 support is not needed any more.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/11
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    bdyw     West boundary array without stagger (local)
!    bdye     East boundary array without stagger (local)
!    bdys     South boundary array without stagger (local)
!    bdyn     North boundary array without stagger (local)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are local staggered size, i.e. scalar 
!                   grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globw2d   West lateral boundary array
!    globe2d   East lateral boundary array
!    globs2d   South lateral boundary array
!    globn2d   North lateral boundary array
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
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

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: bdyzone
  INTEGER, INTENT(IN) :: fzone

  REAL,   INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL,   INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL,   INTENT(OUT) :: globw2d((ny-fzone)*nproc_y+fzone-1,bdyzone)
  REAL,   INTENT(OUT) :: globe2d((ny-fzone)*nproc_y+fzone-1,bdyzone)
  REAL,   INTENT(OUT) :: globs2d((nx-fzone)*nproc_x+fzone-1,bdyzone)
  REAL,   INTENT(OUT) :: globn2d((nx-fzone)*nproc_x+fzone-1,bdyzone)
                    ! Output array in global domain, defined on U points

  REAL,   INTENT(OUT) :: tem1(nx,bdyzone)
  REAL,   INTENT(OUT) :: tem2(ny,bdyzone)

!-----------------------------------------------------------------------
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: master = 0  ! NetCDF is not parallel I/O, So 
                         ! only one processor can access the opened file
  INTEGER :: mptag       ! Unique MPI id.
  INTEGER :: mytag
  INTEGER :: source
  INTEGER :: ia,ja, ic,jc
  INTEGER :: stat(mpi_status_size)
  INTEGER :: imstat
  INTEGER :: i, j, k, bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  West/south boundary, no-stagger, processor 0 do merge
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO j = 1,ny
      tem2(j,bdy) = bdyw(j,1,bdy)
    END DO

    DO i = 1,nx
      tem1(i,bdy) = bdys(i,1,bdy)
    END DO
  END DO

  !
  ! West boudnary
  !
  ic = 1
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 100 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO j = 1,ny-1
          ja = j + (jc-1)*(ny-fzone)
          globw2d(ja,bdy) = tem2(j,bdy)
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! South boudnary
  !
  jc = 1
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 200 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO i = 1,nx-1
          ia = i + (ic-1)*(nx-fzone)
          globs2d(ia,bdy) = tem1(i,bdy)
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

!-----------------------------------------------------------------------
!
!  East/North boundary, U-stagger
!
!-----------------------------------------------------------------------

  DO bdy = 1,bdyzone
    DO j = 1,ny
      tem2(j,bdy) = bdye(j,1,bdy)
    END DO

    DO i = 1,nx
      tem1(i,bdy) = bdyn(i,1,bdy)
    END DO
  END DO

  !
  ! East boudnary
  !
  ic = nproc_x
  DO jc=1,nproc_y

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 300 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem2,ny*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem2,ny*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO j = 1,ny-1
          ja = j + (jc-1)*(ny-fzone)
          globe2d(ja,bdy) = tem2(j,bdy)
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! jc
  CALL mpbarrier

  !
  ! North boudnary
  !
  jc = nproc_y
  DO ic=1,nproc_x

    source = ic+(jc-1)*nproc_x-1

    IF (source /= master) THEN        ! message passing and receiving section...

      mytag = mptag + 400 + ic + jc

      IF(myproc == (ic+(jc-1)*nproc_x-1))THEN  ! pass data to processor 0
       
        CALL mpi_send(tem1,nx*bdyzone,MPI_REAL,master,mytag,         &
                      MPI_COMM_WORLD,imstat) 

      END If

      IF (myproc == master) THEN               ! receive data

        CALL mpi_recv(tem1,nx*bdyzone,MPI_REAL,source,mytag,         &
                      MPI_COMM_WORLD,stat,imstat)
      END IF

    END IF

    IF (myproc == master) THEN  ! store data into global arrays

      DO bdy = 1,bdyzone
        DO i = 1,nx-1
          ia = i + (ic-1)*(nx-fzone)
          globn2d(ia,bdy) = tem1(i,bdy)
        END DO      
      END DO

    END IF                      ! End of storing section

  END DO     ! ic
  CALL mpbarrier

  RETURN
END SUBROUTINE wrf_mergebdy2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE WRF_MPSENDRECV1DE            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mpsendrecv1de(var,nx,ny,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive east boundary 1D data between processors to 
!  update the fake zones.
!
!  NOTE: After this updating, the scalar array with have one extra 
!        valid row, i.e. var(nx,:) is valid now for scalar array.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  10/11/2004
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
!
!  INPUT & OUTPUT
!
!    var      Variable for which boundaries need updating.
!
!  WORK array
!
!    tem      Work array (with a size at least ny x 2).
!
!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny         ! Number of grid points in 
                                       ! x and y directions
  REAL, INTENT(INOUT) :: var(nx,ny)
  REAL, INTENT(INOUT) :: tem(ny,2)   ! Work array.


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
!!
!-----------------------------------------------------------------------
!
  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
  INTEGER :: j

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
  ! Pack send buffer, the first valid slice
  !
  DO j=1,ny
    tem(j,1) = var(1,j)
  END DO

  CALL mpi_sendrecv(tem(:,1),ny,MPI_REAL,dest,  mptag+tag_e,       &
                    tem(:,2),ny,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update the last row
  !
  IF ( loc_x /= nproc_x )  THEN
    DO j=1,ny
      var(nx,j) = tem(j,2)
    END DO
  ELSE                   ! just copy slice
    DO j = 1,ny
      var(nx,j) = var(nx-1,j)
    END DO
  END IF

  RETURN
END SUBROUTINE wrf_mpsendrecv1de
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE WRF_MPSENDRECV1DN            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_mpsendrecv1dn(var,nx,ny,tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Send & receive north boundary data between processors to 
!  update the fake zones.
!
!  NOTE: After this updating, the scalar array with have one extra 
!        valid row, i.e. var(:,ny is valid now for scalar array.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  10/11/2004
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
  REAL, INTENT(INOUT) :: var(nx,ny)
  REAL, INTENT(INOUT) :: tem(nx,2)     ! Work array.

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
  INTEGER :: i

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
  DO i=1,nx
    tem(i,1) = var(i,1)
  END DO

  CALL mpi_sendrecv(tem(:,1),nx,MPI_REAL,dest,  mptag+tag_n,       &
                    tem(:,2),nx,MPI_REAL,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO i=1,nx
      var(i,ny) = tem(i,2)
    END DO
  ELSE
    DO i = 1,nx
      var(i,ny) = var(i,ny-1)
    END DO
  END IF

  RETURN
END SUBROUTINE wrf_mpsendrecv1dn
