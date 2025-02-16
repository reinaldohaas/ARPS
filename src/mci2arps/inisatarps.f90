!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INISATARPS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE inisatarps(arunnam,nx_arps,ny_arps,                          &
                      dx_arps,dy_arps,xnw,ynw)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the ARPS model grid control parameters for the
!  satellite remapper.  They are read from the input file and
!  stored in common blocks and/or returned.
!
!  Note:
!  In order to remain compatible with the ARPS parameter
!  initialization, changes in INITPARA which affect the
!  NAMELIST blocks used here will need to be reflected in
!  this subroutine.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  9/20/1997.
!
!  MODIFICATION HISTORY:
!
!  2001/03/28 (Gene Bassett)
!  Corrected error in mci2arps (with variables xnw & ynw) which was
!  causing satellite data to be shifted over one grid box in x- and two 
!  in the y-direction.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    none, read from standard input in NAMELISTS
!
!  OUTPUT:
!
!    Control parameters in the "globcst.inc' inlcude file and...
!    arunnam   character string describing run (to first comma or space)
!    nx_arps   x-dimension of the remapping grid (ARPS nx)
!    nx_arps   y-dimension of the remapping grid (ARPS ny)
!    dx_arps   grid spacing in x-dimension of the remapping grid (ARPS dx)
!    dy_arps   grid spacing in y-dimension of the remapping grid (ARPS dy)
!    xnw       x coordinate of the NW (upper-left) scalar point
!    ynw       y coordinate of the NW (upper-left) scalar point
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER (LEN=1) :: arunnam(80)
  INTEGER :: nx,ny,nz
  INTEGER :: nx_arps,ny_arps
  REAL    :: dx_arps,dy_arps
  REAL    :: xnw,ynw
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
  REAL :: ctrx,ctry
  REAL :: latnot(2)

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Namelist declarations
!
!-----------------------------------------------------------------------
!
  NAMELIST /grid_dims/ nx, ny, nz
  
  NAMELIST /jobname/ runname

  NAMELIST /grid/ dx,dy,dz,strhopt,dzmin,zrefsfc,dlayer1,dlayer2,       &
            strhtune,zflat,ctrlat,ctrlon, crdorgnopt

  NAMELIST /projection/ mapproj, trulat1,trulat2,trulon, sclfct,        &
            mpfctopt,mptrmopt,maptest
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Set defaults
!  Some are not used here, but set for completeness.
!
!-----------------------------------------------------------------------
!
  nx=67
  ny=67
  nz=35

  runname = 'may20'
  
  dx=1000.
  dy=1000.
  dz=500.

  strhopt  = 0
  dzmin    = 500.0
  zrefsfc  =   0.0
  dlayer1  =   0.0
  dlayer2  =   1.0E5
  strhtune =   1.0
  zflat    =   1.0E5

  ctrlat  =  35.0
  ctrlon  = -100.0

  mapproj  = 0
  trulat1  =   30.0
  trulat2  =   60.0
  trulon   = -100.0
  sclfct   =    1.0
  mpfctopt = 1
  mptrmopt = 1
  maptest = 0
!
!-----------------------------------------------------------------------
!
!  Read namelists from input file
!
!-----------------------------------------------------------------------
!
  READ (5,grid_dims,END=100)
  WRITE(6,'(/a,a)')' Namelist block grid_dims successfully read.'
  
  READ (5,jobname,END=100)
  WRITE(6,'(/a,a)')' Namelist block jobname successfully read.'
  
  READ (5,grid,END=100)
  WRITE(6,'(/a,a)')' Namelist block grid successfully read.'
  
  READ (5,projection,END=100)
  WRITE(6,'(/a,a)')' Namelist block projection successfully read.'

  GO TO 102

  100 CONTINUE
  WRITE(6,'(a)')                                                        &
      'Error reading NAMELIST file. Default values used'

  102 CONTINUE

!
!-----------------------------------------------------------------------
!
!  Set-up map projection variables
!
!-----------------------------------------------------------------------
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
!
!-----------------------------------------------------------------------
!
!  Find coordinate of the NW scalar corner, that is the
!  origin in satellite projection space.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
  xnw=ctrx-dx*(1.+0.5*(nx-3))
  ynw=ctry+dy*(2.+0.5*(ny-3))
!
!-----------------------------------------------------------------------
!
!  Set return nx,ny,dx,dy variables
!  The C-calling routine doesn't have access to the common blocks.
!
!-----------------------------------------------------------------------
!
  DO i=1,80
    IF (runname(i:i) == ' ' .OR. runname(i:i) == ',') EXIT
    arunnam(i)=runname(i:i)
  END DO
  21 CONTINUE
!
  nx_arps=nx
  ny_arps=ny
  dx_arps=dx
  dy_arps=dy
!
  RETURN
END SUBROUTINE inisatarps
