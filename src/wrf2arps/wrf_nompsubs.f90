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
!  For no-mp mode, just copy the array from globvar to var.
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
  REAL,    INTENT(IN)  :: globvar(nx,ny)

  REAL,    INTENT(OUT) :: var(nx,ny)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  var(:,:) = globvar(1:nx,1:ny)     

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
!  For no-mp mode, just copy the array from globvar to var.
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
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny       ! Number of grid points in x and y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  INTEGER, INTENT(IN)  :: globvar(nx,ny)
  INTEGER, INTENT(OUT) :: var(nx,ny)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  var(:,:) = globvar(1:nx,1:ny)     

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
!  For no-mp mode, just copy the array from globvar to var.
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
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y and z
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL,    INTENT(IN)  :: globvar(nx,ny,nz)

  REAL,    INTENT(OUT) :: var(nx,ny,nz)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  var(:,:,:) = globvar(1:nx,1:ny,1:nz)     

  RETURN
END SUBROUTINE wrf_split3d

SUBROUTINE mpext_wrf_u(var,nx,ny,nz,tem)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in 
                                       ! x, y and z directions
  REAL, INTENT(INOUT) :: var(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem(ny,nz,2)   ! Work array.

  RETURN
END SUBROUTINE mpext_wrf_u

SUBROUTINE mpext_wrf_v(var,nx,ny,nz,tem)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz      ! Number of grid points in 
                                       ! x, y and z directions
  REAL, INTENT(INOUT) :: var(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem(ny,nz,2)   ! Work array.

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
  REAL,    INTENT(INOUT) :: tem1(nx,nz,2), tem2(ny,nz,2)

  INTEGER  :: i,j,k

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

  !
  ! South & north boundary
  !
  DO k = 1,nz
    DO i = 1,nx
      uext(i,0,k) = u(i,1,k)
      uext(i,ny,k) = u(i,ny-1,k)
    END DO
  END DO

  !
  ! East boundary
  !
  DO k = 1,nz
    DO j = 1,ny
      uext(nx+1,j,k) = u(nx,j,k)
    END DO
  END DO

  uext(nx+1,0,:) = u(nx,1,:)

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
  REAL,    INTENT(INOUT) :: tem1(nx,nz,2), tem2(ny,nz,2)

  INTEGER  :: i,j,k

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

  !
  ! West & east boundary
  !
  DO k = 1,nz
    DO j = 1,ny
      vext(0,j,k) = v(1,j,k)
      vext(nx,j,k) = v(nx-1,j,k)
    END DO
  END DO

  !
  ! North boundary
  !
  DO k = 1,nz
    DO i = 1,nx
      vext(i,ny+1,k) = v(i,ny,k)
    END DO
  END DO

  vext(0,ny+1,:) = v(1,ny,:)

  RETURN
END SUBROUTINE extend_v

!subroutine printwrf3d(unt,varname,var,stagger,nx,ny,nz)
!  
!  implicit none
!  integer      :: unt
!  character(*) :: varname
!  character(*) :: stagger
!  integer      :: nx, ny,nz
!  real         :: var(nx,ny,nz)
!
!  integer i,j,k
!
!  write(unt,'(3a)') '--- ',varname,' ---'
!
!  IF (stagger == 'X') THEN
!    DO k = 1,nz-1
!      DO j = 1,ny-1
!        WRITE(UNIT=unt,FMT='(2I3,a)') k,j,' --- '
!        WRITE(UNIT=unt,FMT='(5f16.5)') (var(i,j,k),i=1,nx)
!      END DO
!    END DO
!  ELSE IF (stagger == 'Y') THEN
!    DO k = 1,nz-1
!      DO j = 1,ny
!        WRITE(UNIT=unt,FMT='(2I3,a)') k,j,' --- '
!        WRITE(UNIT=unt,FMT='(5f16.5)') (var(i,j,k),i=1,nx-1)
!      END DO
!    END DO
!  ELSE IF (stagger == 'Z') THEN
!    DO k = 1,nz
!      DO j = 1,ny-1
!        WRITE(UNIT=unt,FMT='(2I3,a)') k,j,' --- '
!        WRITE(UNIT=unt,FMT='(5f16.5)') (var(i,j,k),i=1,nx-1)
!      END DO
!    END DO
!  ELSE
!    DO k = 1,nz-1
!      DO j = 1,ny-1
!        WRITE(UNIT=unt,FMT='(2I3,a)') k,j,' --- '
!        WRITE(UNIT=unt,FMT='(5f16.5)') (var(i,j,k),i=1,nx-1)
!      END DO
!    END DO
!  END IF
!
!  return
!end subroutine printwrf3d
