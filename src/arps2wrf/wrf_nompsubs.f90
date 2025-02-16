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
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE3DT               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3dt(var,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at mass grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/30
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       3rd dimension of the subdomain array, possible value are
!             vertical grid points (nz in other subroutines), nzsoil,
!             nstyps+1, or 4 (prcrate) or 1 (for 2D arrays)
!
!  OUTPUT:
!
!    globvar  array defined at mass grid points returned.
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
  REAL,    INTENT(IN)  :: var(nx,ny,nz)

  REAL,    INTENT(OUT) :: globvar(nx-1,ny-1,nz-1)

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        globvar(i,j,k) = var(i,j,k)     
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge3dt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE3DU               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3du(var,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at U grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/30
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       3rd dimension of the subdomain array, possible value are
!             vertical grid points (nz in other subroutines), nzsoil,
!             nstyps+1, or 4 (prcrate) or 1 (for 2D arrays)
!
!  OUTPUT:
!
!    globvar  array defined at U grid points returned.
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
  REAL,    INTENT(IN)  :: var(nx,ny,nz)

  REAL,    INTENT(OUT) :: globvar(nx,ny-1,nz-1)

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx
        globvar(i,j,k) = var(i,j,k)     
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge3du
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE3DV               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3dv(var,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at V grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/30
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       3rd dimension of the subdomain array, possible value are
!             vertical grid points (nz in other subroutines), nzsoil,
!             nstyps+1, or 4 (prcrate) or 1 (for 2D arrays)
!
!  OUTPUT:
!
!    globvar  array defined at V grid points returned.
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
  REAL,    INTENT(IN)  :: var(nx,ny,nz)

  REAL,    INTENT(OUT) :: globvar(nx-1,ny,nz-1)

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = 1,nz-1
    DO j = 1,ny
      DO i = 1,nx-1
        globvar(i,j,k) = var(i,j,k)     
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge3dv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE3DW               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge3dw(var,nx,ny,nz,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at W grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/09/30
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       3rd dimension of the subdomain array, possible value are
!             vertical grid points (nz in other subroutines), nzsoil,
!             nstyps+1, or 4 (prcrate) or 1 (for 2D arrays)
!
!  OUTPUT:
!
!    globvar  array defined at W grid points returned.
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
  REAL,    INTENT(IN)  :: var(nx,ny,nz)

  REAL,    INTENT(OUT) :: globvar(nx-1,ny-1,nz)

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = 1,nz
    DO j = 1,ny-1
      DO i = 1,nx-1
        globvar(i,j,k) = var(i,j,k)     
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge3dw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE2DT               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2dt(var,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at mass grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/03
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    globvar  array defined at mass grid points returned.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny      ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL,    INTENT(IN)  :: var(nx,ny)

  REAL,    INTENT(OUT) :: globvar(nx-1,ny-1)

  INTEGER :: i,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO j = 1,ny-1
    DO i = 1,nx-1
      globvar(i,j) = var(i,j)     
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge2dt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE2DU               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2du(var,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at U grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/03
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    globvar  array defined at U grid points returned.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny      ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL,    INTENT(IN)  :: var(nx,ny)

  REAL,    INTENT(OUT) :: globvar(nx,ny-1)

  INTEGER :: i,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO j = 1,ny-1
    DO i = 1,nx
      globvar(i,j) = var(i,j)     
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge2du
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE2DV               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2dv(var,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at V grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/03
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    globvar  array defined at V grid points returned.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny      ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL,    INTENT(IN)  :: var(nx,ny)

  REAL,    INTENT(OUT) :: globvar(nx-1,ny)

  INTEGER :: i,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO j = 1,ny
    DO i = 1,nx-1
      globvar(i,j) = var(i,j)     
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge2dv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGE2DI               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrf_merge2di(var,nx,ny,fzone,globvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  For no-mp mode, just copy the array from var to globvar.
!  Globvar is defined at mass grid points.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Yunheng Wang
!  2004/10/03
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    var      Variable defined at no-staggered grid.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    globvar  array defined at mass grid points returned.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny      ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  INTEGER, INTENT(IN)  :: var(nx,ny)

  INTEGER, INTENT(OUT) :: globvar(nx-1,ny-1)

  INTEGER :: i,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO j = 1,ny-1
    DO i = 1,nx-1
      globvar(i,j) = var(i,j)     
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_merge2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGEBDYU              ######
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
!  For no-mp mode, just copy the boundary arrays to U staggered grids.
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
!    bdyw     West boundary array without stagger
!    bdye     East boundary array without stagger
!    bdys     South boundary array without stagger
!    bdyn     North boundary array without stagger
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are staggered size, i.e. scalar grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globwt   West lateral boundary array staggered at Mass points
!    globet   East lateral boundary array staggered at Mass points
!    globsu   South lateral boundary array staggered at U points
!    globnu   North lateral boundary array staggered at U points
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: bdyzone
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL, INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL, INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL, INTENT(OUT) :: globwt(ny-1,nz-1,bdyzone)
  REAL, INTENT(OUT) :: globet(ny-1,nz-1,bdyzone)
  REAL, INTENT(OUT) :: globsu(nx,  nz-1,bdyzone)
  REAL, INTENT(OUT) :: globnu(nx,  nz-1,bdyzone)

  REAL, INTENT(INOUT) :: tem1(ny,nz,bdyzone)      ! work arrays
  REAL, INTENT(INOUT) :: tem2(nx,nz,bdyzone)      ! not used.

  INTEGER :: i,j,k,bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO bdy = 1,bdyzone
    DO k = 1,nz-1
      DO i = 1,nx
         globsu(i,k,bdy) = bdys(i,k,bdy)
         globnu(i,k,bdy) = bdyn(i,k,bdy)
      END DO

      DO j = 1,ny-1
         globwt(j,k,bdy) = bdyw(j,k,bdy)
         globet(j,k,bdy) = bdye(j,k,bdy)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_mergebdyu
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGEBDYV              ######
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
!  For no-mp mode, just copy the boundary arrays to V staggered grids.
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
!    bdyw     West boundary array without stagger
!    bdye     East boundary array without stagger
!    bdys     South boundary array without stagger
!    bdyn     North boundary array without stagger
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are staggered size, i.e. scalar grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globwv   West lateral boundary array staggered at V points
!    globev   East lateral boundary array staggered at V points
!    globst   South lateral boundary array staggered at Mass points
!    globnt   North lateral boundary array staggered at Mass points
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: bdyzone
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL, INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL, INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL, INTENT(OUT) :: globwv(ny,  nz-1,bdyzone)
  REAL, INTENT(OUT) :: globev(ny,  nz-1,bdyzone)
  REAL, INTENT(OUT) :: globst(nx-1,nz-1,bdyzone)
  REAL, INTENT(OUT) :: globnt(nx-1,nz-1,bdyzone)

  REAL, INTENT(INOUT) :: tem1(ny,nz,bdyzone)      ! work arrays
  REAL, INTENT(INOUT) :: tem2(nx,nz,bdyzone)      ! not used.

  INTEGER :: i,j,k,bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO bdy = 1,bdyzone
    DO k = 1,nz-1
      DO i = 1,nx-1
         globst(i,k,bdy) = bdys(i,k,bdy)
         globnt(i,k,bdy) = bdyn(i,k,bdy)
      END DO

      DO j = 1,ny
         globwv(j,k,bdy) = bdyw(j,k,bdy)
         globev(j,k,bdy) = bdye(j,k,bdy)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_mergebdyv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGEBDYW              ######
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
!  For no-mp mode, just copy the boundary arrays to W staggered grids.
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
!    bdyw     West boundary array without stagger
!    bdye     East boundary array without stagger
!    bdys     South boundary array without stagger
!    bdyn     North boundary array without stagger
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are staggered size, i.e. scalar grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globww   West lateral boundary array staggered at W points
!    globew   East lateral boundary array staggered at W points
!    globsw   South lateral boundary array staggered at w points
!    globnw   North lateral boundary array staggered at w points
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: bdyzone
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL, INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL, INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL, INTENT(OUT) :: globww(ny-1,nz,bdyzone)
  REAL, INTENT(OUT) :: globew(ny-1,nz,bdyzone)
  REAL, INTENT(OUT) :: globsw(nx-1,nz,bdyzone)
  REAL, INTENT(OUT) :: globnw(nx-1,nz,bdyzone)

  REAL, INTENT(INOUT) :: tem1(ny,nz,bdyzone)      ! work arrays
  REAL, INTENT(INOUT) :: tem2(nx,nz,bdyzone)      ! not used.

  INTEGER :: i,j,k,bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO bdy = 1,bdyzone
    DO k = 1,nz
      DO i = 1,nx-1
         globsw(i,k,bdy) = bdys(i,k,bdy)
         globnw(i,k,bdy) = bdyn(i,k,bdy)
      END DO

      DO j = 1,ny-1
         globww(j,k,bdy) = bdyw(j,k,bdy)
         globew(j,k,bdy) = bdye(j,k,bdy)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_mergebdyw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGEBDYT              ######
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
!  For no-mp mode, just copy the boundary arrays to Mass staggered grids.
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
!    bdyw     West boundary array without stagger
!    bdye     East boundary array without stagger
!    bdys     South boundary array without stagger
!    bdyn     North boundary array without stagger
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are staggered size, i.e. scalar grid size plus 1.
!
!    bdyzone  Number of lateral boundary zone
!    fzone    Message passing overlay (fake) zone. It should be 1 for WRF
!             and 3 for ARPS
!
!  OUTPUT:
!
!    globwt   West lateral boundary array staggered at scalar points
!    globet   East lateral boundary array staggered at scalar points
!    globst   South lateral boundary array staggered at scalar points
!    globnt   North lateral boundary array staggered at scalar points
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: bdyzone
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL, INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL, INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL, INTENT(OUT) :: globwt(ny-1,nz-1,bdyzone)
  REAL, INTENT(OUT) :: globet(ny-1,nz-1,bdyzone)
  REAL, INTENT(OUT) :: globst(nx-1,nz-1,bdyzone)
  REAL, INTENT(OUT) :: globnt(nx-1,nz-1,bdyzone)

  REAL, INTENT(INOUT) :: tem1(ny,nz,bdyzone)      ! work arrays
  REAL, INTENT(INOUT) :: tem2(nx,nz,bdyzone)      ! not used.

  INTEGER :: i,j,k,bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO bdy = 1,bdyzone
    DO k = 1,nz-1
      DO i = 1,nx-1
         globst(i,k,bdy) = bdys(i,k,bdy)
         globnt(i,k,bdy) = bdyn(i,k,bdy)
      END DO

      DO j = 1,ny-1
         globwt(j,k,bdy) = bdyw(j,k,bdy)
         globet(j,k,bdy) = bdye(j,k,bdy)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_mergebdyt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRF_MERGEBDY2d             ######
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
!  For no-mp mode, just copy the boundary arrays to 2D boundary arrays
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
!    bdyw     West boundary array without stagger
!    bdye     East boundary array without stagger
!    bdys     South boundary array without stagger
!    bdyn     North boundary array without stagger
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (botthom/top)
!             Note: All are staggered size, i.e. scalar grid size plus 1.
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

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nx,ny,nz    ! Number of grid points in x, y
  INTEGER, INTENT(IN)  :: bdyzone
  INTEGER, INTENT(IN)  :: fzone       ! number of fake zone
                                      ! 1 for wrf
                                      ! 3 for arps
  REAL, INTENT(IN)  :: bdyw(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdye(ny,nz,bdyzone)
  REAL, INTENT(IN)  :: bdys(nx,nz,bdyzone)
  REAL, INTENT(IN)  :: bdyn(nx,nz,bdyzone)

  REAL, INTENT(OUT) :: globw2d(ny-1,bdyzone)
  REAL, INTENT(OUT) :: globe2d(ny-1,bdyzone)
  REAL, INTENT(OUT) :: globs2d(nx-1,bdyzone)
  REAL, INTENT(OUT) :: globn2d(nx-1,bdyzone)

  REAL, INTENT(INOUT) :: tem1(ny,nz,bdyzone)      ! work arrays
  REAL, INTENT(INOUT) :: tem2(nx,nz,bdyzone)      ! not used.

  INTEGER :: i,j,k,bdy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO bdy = 1,bdyzone
    DO i = 1,nx-1
       globs2d(i,bdy) = bdys(i,1,bdy)
       globn2d(i,bdy) = bdyn(i,1,bdy)
    END DO

    DO j = 1,ny-1
       globw2d(j,bdy) = bdyw(j,1,bdy)
       globe2d(j,bdy) = bdye(j,1,bdy)
    END DO
  END DO

  RETURN
END SUBROUTINE wrf_mergebdy2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE WRF_MPSENDRECV1DE              ######
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
!  Send & receive east/west boundary 1D data between processors to 
!  update the fake zones.
!  
!  NOTE:
!    For no-mpi mode, just copy the last valid row (nx-1) to fake
!    zone (row nx).
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
  REAL, INTENT(INOUT) :: tem(ny,2)     ! Work array, not used for no-mpi mode
 
  INTEGER :: j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO j = 1,ny
    var(nx,j) = var(nx-1,j)
  END DO

  RETURN
END SUBROUTINE wrf_mpsendrecv1de
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE WRF_MPSENDRECV1DN              ######
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
!  Send & receive north boundary 1D data between processors to 
!  update the fake zones.
!  
!  NOTE:
!    For no-mpi mode, just copy the last valid column (ny-1) to fake
!    zone (column ny).
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
!  INPUT & OUTPUT
!
!    var      Variable for which boundaries need updating.
!
!  WORK array
!
!    tem      Work array (with a size at least nx x 2).
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
  REAL, INTENT(INOUT) :: tem(nx,2)     ! Work array, not used for no-mpi mode
 
  INTEGER :: i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO i = 1,nx
    var(i,ny) = var(i,ny-1)
  END DO

  RETURN
END SUBROUTINE wrf_mpsendrecv1dn
