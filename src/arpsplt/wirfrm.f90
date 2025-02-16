
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WIRFRM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wirfrm(a,nx,ist,iend,ny,jst,jend,nz,kst,kend,                &
                  valiso, slab)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Use NCARgraphyic routinue "isosrf" to generate 3-D wirefrm-plot.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/21/1992
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation
!
!  01/13/2003 (Y. Wang)
!  Added message passing code.
!
!-----------------------------------------------------------------------
!
!  Argument list:
!
!  a        3-D array of data that defines the iso-valued surface.
!  nx       first dimension of a
!  ist      index of first i grid point to be used. Global index!!
!  iend     index of last  i grid point to be used. Global index!!
!
!  ny       second dimension of a
!  jst      index of first j grid point to be used. Global index!!
!  jend     index of last  j grid point to be used. Global index!!
!
!  nz       third dimension of a
!  kst      index of first k grid point to be used.
!  kend     index of last  k grid point to be used.
!
!  valiso   the iso-value used to define the surface.
!  slab     a working space for internal storage. It is larger that
!           (max(nx,ny,nz)+2)*2
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! Number of grid points in 3 directions
  INTEGER :: ist,iend,jst,jend,kst,kend
                             ! indexes of grid points in 3 directions
  INTEGER :: iflag,muvw

  REAL :: a(nx,ny,nz)       ! 3-D array of data that defines
                            ! the iso-valued surface
  REAL :: valiso            ! the iso-value used to define the surface
!
!-----------------------------------------------------------------------
!
!  Working arrays used for the wrieframe plotting.
!
!-----------------------------------------------------------------------
!
  REAL :: eye(3)
  REAL :: slab( nx+1, nx+1)

  INCLUDE 'mp.inc'
 
  INTEGER :: nxlg, nylg
  REAL, ALLOCATABLE :: tema(:,:,:), temlab(:,:,:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  muvw=MAX(nxlg,nylg,nz)+2

!-----------------------------------------------------------------------
!
!  EYE: The position of the eye in 3-D space. A tentative choice
!       can be (5.0*nx, 3.5*ny, 2.0*nz)
!
!-----------------------------------------------------------------------

  eye(1)= 1.0*nxlg
  eye(2)=-1.0*nylg
  eye(3)= 1.5*nz
!
!-----------------------------------------------------------------------
!
!  IFLAG: Define the types of lines to be drawn.
!         =1, lines of constant z only
!         =2, lines of constant y only
!         =3, lines of constant z and y.
!         =4, lines of constant x only
!         =5, lines of constant x and z.
!         =6, lines of constant x and y.
!         =7, or more, lines of constant x, y and z.
!  The sign of IFLAG determines what is inside and what is outside
!  the isosurface.
!  IFLAG >0, values greater than valiso are assumed inside the
!            surface.
!  IFLAG <0, values less than valiso are assumed inside the surface.
!
!-----------------------------------------------------------------------
!
  iflag=4
  IF( valiso < 0.0) iflag = -iflag

  ALLOCATE(tema(nxlg,nylg,nz))
  ALLOCATE(temlab(nxlg,nylg,nz))

  CALL mpimerge3d(a,nx,ny,nz,tema)

  IF(myproc == 0) THEN
    CALL xnwfrm
    CALL xnwpic
  
    CALL gselnt(0)
  
    CALL isosrf(tema(ist,jst,kst),nxlg,iend-ist+1,                    &
                nylg,jend-jst+1, kend-kst+1,                          &
                eye,muvw,temlab,valiso,iflag)
  
    CALL xnwfrm
  END IF   ! myproc == 0

  DEALLOCATE(tema, temlab)

  RETURN
END SUBROUTINE wirfrm
