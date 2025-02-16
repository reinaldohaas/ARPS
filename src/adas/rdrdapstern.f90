PROGRAM rdrdapstrn
!
!-----------------------------------------------------------------------
!
!  Reads KMA RDAPS terrain file and converts it to an arps terrain
!  file at the same resolution and map projection of the
!  original KMA RDAPS grid.
!
!  Compile:
!
!  f77 -o rdrdapstrn rdrdapstrn.f maproj3d.f genlib3d.f outlib3d.f
!
!  The output can then be used in PROGRAM MRGTRN
!
!  Keith Brewster, July, 1996
!
!  MODIFICATIONS
!
!  06/10/97 Keith Brewster
!  Version for RDAPS based on RUC version.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny
  PARAMETER (nx=120,ny=104)
  INTEGER :: mapproj
  REAL :: trulat1,trulat2
  REAL :: trulon, sclfct
  REAL :: dx,dy
  PARAMETER ( mapproj = 2,      & ! Lambert Conformal Conic
              trulat1 = 30.,                                            &
              trulat2 = 60.,                                            &
              trulon  = 125.,                                           &
              sclfct  = 1.,                                             &
              dx = 40000.,                                              &
              dy = 40000.)

  CHARACTER(LEN=80), PARAMETER :: rdapstrn    = 'rgm_topo40.dat'
  CHARACTER(LEN=80), PARAMETER :: rdapslatlon = 'rgm40_crosspnt_latlon.dat'
  CHARACTER(LEN=80), PARAMETER :: rdapsout    = 'arpsrdapstrn.dat'

  REAL :: hterain(nx,ny)
  REAL :: lat(nx,ny)
  REAL :: lon(nx,ny)

!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=80) :: head
  REAL :: x0,y0,x1,y1,xctr,yctr,dxin,dyin
  REAL :: ctrlat,ctrlon
  REAL :: latnot(2)
  REAL :: gamma,knot,tol,rdnot,wdn,rdummy
  INTEGER :: i,j,ierr
  INTEGER :: analtype,ipass,itertype,idummy
!
!-----------------------------------------------------------------------
!
!  Set terrain analysis variables to dummy numbers
!
!-----------------------------------------------------------------------
!
  gamma=0.
  knot=0.
  tol=0.
  rdnot=0.
  wdn=0.
  rdummy=0.
  analtype=0.
  ipass=0
  itertype=0
  idummy=0
!
!-----------------------------------------------------------------------
!
!  Set map projection according to the parameters set above.
!
!-----------------------------------------------------------------------
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
!
!-----------------------------------------------------------------------
!
!  Read the terrain, lat and lon from the maps terrain file,
!  this was obtained from FSL.
!
!-----------------------------------------------------------------------
!
!  CALL asnctl ('NEWLOCAL', 1, ierr)
!  CALL asnfile(rdapslatlon, '-F f77 -N ieee', ierr)

  OPEN(31,FILE=trim(rdapslatlon),STATUS='old',FORM='formatted')
  READ(31,'(a80)')    head
  READ(31,'(20F8.3)') ((lat(i,j),i=1,nx),j=1,ny)
  READ(31,'(a80)')    head
  READ(31,'(20F8.3)') ((lon(i,j),i=1,nx),j=1,ny)
  CLOSE(31)

!  CALL asnctl ('NEWLOCAL', 1, ierr)
!  CALL asnfile(rdapstrn, '-F f77 -N ieee', ierr)
  OPEN(31,FILE=trim(rdapstrn),STATUS='old',FORM='unformatted')
  READ(31)((hterain(i,j),j=1,ny),i=1,nx)
  CLOSE(31)
!
!-----------------------------------------------------------------------
!
!  Calculate ctrlat and ctrlon
!  This definition of "center" is to be consistent with ARPS
!  where the scalar grid runs from 1 to nx-1.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(1,1,lat(1,1),lon(1,1),x0,y0)
  CALL lltoxy(1,1,lat(2,2),lon(2,2),x1,y1)
  xctr=x0+0.5*FLOAT(nx-3)*dx
  yctr=y0+0.5*FLOAT(ny-3)*dy
  CALL xytoll(1,1,xctr,yctr,ctrlat,ctrlon)
  WRITE(6,'(a,2f11.4)') ' LL RDAPS corner ',lat(1,1),lon(1,1)
  WRITE(6,'(a,2f11.4)') ' UR RDAPS corner ',lat(nx,ny),lon(nx,ny)
  WRITE(6,'(a,2f11.4)') ' RDAPS center    ',ctrlat,ctrlon
!
  dxin=x1-x0
  dyin=y1-y0
  WRITE(6,'(a,2f11.4)') ' dx from latlon grid',dxin
  WRITE(6,'(a,2f11.4)') ' dy from latlon grid',dyin
!
!-----------------------------------------------------------------------
!
!  Write the analyzed terrain data (model grid data) into file
!  arpstern.dat.
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(rdapsout, '-F f77 -N ieee', ierr)

  OPEN(11,FILE=trim(rdapsout),FORM='unformatted',                       &
          STATUS='unknown')

  WRITE(11) nx,ny

  idummy = 0
  rdummy = 0.0

  WRITE(11) analtype,mapproj,itertype,ipass,idummy,                     &
            idummy,idummy,idummy,idummy,idummy,                         &
            idummy,idummy,idummy,idummy,idummy,                         &
            idummy,idummy,idummy,idummy,idummy

  WRITE(11) dx    ,dy    ,ctrlat,ctrlon,knot ,                          &
            gamma ,trulat1,trulat2,trulon,sclfct,                       &
            tol   ,wdn   ,rdnot ,rdummy,rdummy,                         &
            rdummy,rdummy,rdummy,rdummy,rdummy

  WRITE(11) hterain

  CLOSE(11)
  STOP
END PROGRAM rdrdapstrn
