PROGRAM rdructrn
!
!-----------------------------------------------------------------------
!
!  Reads RUC terrain file and converts it to an arps terrain
!  file at the same resolution and map projection of the
!  original MAPS/RUC grid.
!
!  Compile:
!
!  f77 -o rdructern rdructern.f maproj3d.f
!
!  The output can then be used in PROGRAM MRGTRN
!
!  Keith Brewster, July, 1996
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny
  PARAMETER (nx=81,ny=62)
  INTEGER :: mapproj
  REAL :: trulat1,trulat2
  REAL :: trulon, sclfct
  REAL :: dx,dy
  PARAMETER ( mapproj = 1,      & ! Polar Stereographic
              trulat1 = 40.,                                            &
              trulat2 = 40.,                                            &
              trulon  = -105.,                                          &
              sclfct  = 1.,                                             &
              dx = 60000.,                                              &
              dy = 60000.)

  CHARACTER (LEN=80), PARAMETER :: ructern   = '/vortex/ructern/ructern.dat'
  CHARACTER (LEN=80), PARAMETER :: ruclatlon = '/vortex/ructern/ruclatlon.dat'

  REAL :: hterain(nx,ny)
  REAL :: lat(nx,ny)
  REAL :: lon(nx,ny)

!  Misc local variables
!
  REAL :: x0,y0,xctr,yctr
  REAL :: ctrlat,ctrlon
  REAL :: latnot(2)
  REAL :: gamma,knot,tol,rdnot,wdn,rdummy
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
  OPEN(31,FILE=trim(ruclatlon),STATUS='old',FORM='formatted')
  READ(31,805) lat
  READ(31,805) lon
  805 FORMAT(10F8.3)
  CLOSE(31)
  OPEN(31,FILE=trim(ructern),STATUS='old',FORM='formatted')
  READ(31,805) hterain
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
  xctr=x0+0.5*FLOAT(nx-3)*dx
  yctr=y0+0.5*FLOAT(ny-3)*dy
  CALL xytoll(1,1,xctr,yctr,ctrlat,ctrlon)
  WRITE(6,'(a,2f11.4)') ' LL RUC corner ',lat(1,1),lon(1,1)
  WRITE(6,'(a,2f11.4)') ' UR RUC corner ',lat(nx,ny),lon(nx,ny)
  WRITE(6,'(a,2f11.4)') ' RUC center    ',ctrlat,ctrlon
!
!-----------------------------------------------------------------------
!
!  Write the analyzed terrain data (model grid data) into file
!  arpstern.dat.
!
!-----------------------------------------------------------------------
!
  OPEN(11,FILE='arpsructern.dat',FORM='unformatted',                    &
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
END PROGRAM rdructrn
