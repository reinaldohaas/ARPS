subroutine getcoordinate(nx,ny,lat,lon)
!
!
!-----------------------------------------------------------------------
!
!  get coordinate of arps grid
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!
!  First written: 04/06/2006.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'

  INTEGER :: nx, ny

  REAL :: x(nx)         ! The x-coord. of the u grid point (m)
  REAL :: y(ny)         ! The y-coord. of the v grid point (m)

  REAL :: xh(nx,ny)      ! X-coord. of the terrain (scalar) point.
  REAL :: yh(nx,ny)      ! Y-coord. of the terrain (scalar) point.
  REAL :: xh1d(0:nx)
  REAL :: yh1d(0:ny)

  REAL :: lat(0:nx,0:ny)    ! Latitude of each terrain point
  REAL :: lon(0:nx,0:ny)    ! Longitude of each terrain point

!  NAMELIST /grid/ dx,dy,ctrlat,ctrlon
!  NAMELIST /projection/ mapproj,trulat1,trulat2,trulon,sclfct

!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
  REAL :: trulat (2)           ! Latitude at which the map projection
                               ! is true (degees N)

  REAL :: latmax,latmin,lonmax,lonmin ! max/min lat/lon of the projected grid

  INTEGER :: i,j,ih,jh

  INTEGER :: imin,imax,jmin,jmax,kmin,kmax, istatus
   
  REAL :: swx,swy,ctrx,ctry
!
!-----------------------------------------------------------------------
!
!  TEMPORARY ARRAYS:
!
!-----------------------------------------------------------------------
!


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( mapproj == 0 ) THEN
    trulat1 = ctrlat
    trulat2 = ctrlat
    trulon  = ctrlon
  END IF

!-----------------------------------------------------------------------
!

  IF( ctrlat > 90 .OR. ctrlat < -90. )then
    PRINT *,'ctrlat is too large or too small. ',      &
            'must be set between 0. AND +90. degree north.'
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  test to see if trulat is within limits...
!
!-----------------------------------------------------------------------
!
  IF ( trulat1 > 90 .OR. trulat1 < -90 ) then
    PRINT *,'Input parameter trulat1 is too large or too small.',       &
            'must be set between 0. AND +90. degree north.'
    STOP
  END IF

  IF ( trulat2 > 90 .OR. trulat2 < -90  ) then
    PRINT *,'Input parameter trulat2 is too large or too small.',       &
            'must be set between 0. AND +90. degree north.'
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Set up map projection
!
!-----------------------------------------------------------------------
!
  trulat(1) = trulat1
  trulat(2) = trulat2

  CALL setmapr( mapproj, 1.0 , trulat , trulon)
                             ! set up parameters for map projection

  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
  swx = ctrx - (FLOAT(nx-3)/2.) * dx
  swy = ctry - (FLOAT(ny-3)/2.) * dy

  CALL setorig( 1, swx, swy)

!
!-----------------------------------------------------------------------
!
!  Set the horizontal coordinates of the model grid
!
!-----------------------------------------------------------------------
!
  DO i=1,nx
    x(i) = dx * (i-2)
  END DO

  DO j=1,ny
    y(j) = dy * (j-2)
  END DO

  DO j=1,ny
    DO i=1,nx
      xh(i,j)=(i-1.5)*dx
      yh(i,j)=(j-1.5)*dy
    END DO
  END DO

  DO j=0,ny
    DO i=0,nx
      xh1d(i)=(i-1.5)*dx
      yh1d(j)=(j-1.5)*dy
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Final out the lat/lon of each scalar (terrain data) point.
!
!-----------------------------------------------------------------------
!
  CALL xytoll(nx+1,ny+1,xh1d,yh1d,lat,lon)
!
!-----------------------------------------------------------------------
!
!  Final out the max/min lat/lon of each scalar (terrain data) point.
!
!-----------------------------------------------------------------------

   write(*,*) ' END of get coordinate'

  return

  980   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'Error occured when reading namelist input file. Program stopped.'

  STOP

  990   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'End of file reached when reading namelist input file. ',         &
      'Program stopped.'

  STOP
END subroutine getCoordinate

subroutine getVertcoordinate(nx,ny,nz,zp)
!
!
!-----------------------------------------------------------------------
!
!  get coordinate of arps grid
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!
!  First written: 04/06/2006.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'

  INTEGER :: nx, ny, nz

  REAL :: z(nz)           ! z-coord. of the computational grid.
                          ! Defined at w-point on the staggered grid.
  REAL :: zp(nx,ny,nz)    ! Physical height coordinate defined at
                          ! w-point of the staggered grid.

  REAL :: hterain(nx,ny)  ! Terrain height.

  REAL :: zp1d (nz)       ! Temporary array
  REAL :: dzp1d(nz)       ! Temporary array
  REAL :: tem1(nx,ny,nz)  ! Temporary array

!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: zflat1,z1,z2

  REAL :: zpmax
!
!-------------------------------
!
!
!
  DO k=1,nz
    z(k) = zrefsfc + (k-2) * dz
  END DO
!
!
!-----------------------------------------------------------------------
!
!  Specify the terrain
!
!-----------------------------------------------------------------------
!
  IF( ternopt == 0 ) THEN     ! No terrain, the ground is flat

    DO j=1,ny-1
      DO i=1,nx-1
        hterain(i,j) = zrefsfc
      END DO
    END DO

  ELSE IF( ternopt == 1 ) THEN  ! Bell-shaped mountain
    write(*,*) ' Wrong choice of ternopt!'
    stop 123
  ELSE IF( ternopt == 2 ) THEN          ! Read from terrain data base

!
!-----------------------------------------------------------------------
!
!    Read in the terrain data.
!
!-----------------------------------------------------------------------
!
!  blocking inserted for ordering i/o for message passing
!    DO i=0,nprocs-1,readstride
!      IF(myproc >= i.AND.myproc <= i+readstride-1)THEN
!
!        IF (mp_opt > 0 .AND. readsplit > 0) THEN
!
!        CALL readsplittrn( nx,ny,dx,dy, terndta,                        &
!               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
!               hterain )
!!
!
!        ELSE

        CALL readtrn( nx,ny,dx,dy, terndta,                             &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
               hterain )

!        END IF
!
!      END IF
!
!      IF (mp_opt > 0) CALL mpbarrier
!    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set up a stretched vertical grid.
!
!  For strhopt=1, function y = a+b*x**3 is used to specify dz as a
!                              function of k.
!  For strhopt=2, function y = c + a*tanh(b*x) is used to specify dz
!                              as a function of k.
!
!-----------------------------------------------------------------------
!
  IF ( strhopt == 0 ) THEN


    DO k=1,nz
      zp1d(k) = z(k)
    END DO

  ELSE IF ( strhopt == 1 .OR.strhopt == 2 ) THEN

    z1 = zrefsfc + MAX(0.0, MIN(dlayer1, z(nz-2)-zrefsfc ))
    z2 = z1      + MAX(0.0, MIN(dlayer2, z(nz-1)-z1      ))

    IF( dlayer1 >= (nz-3)*dzmin ) THEN
      WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                             &
          'Can not setup a vertical grid with uniform dz=',dzmin,       &
          ' over the depth of ',dlayer1,' please specify a smaller ',   &
          'value of dlayer1. Program stopped INIGRD.'
      CALL arpsstop('arpsstop called from INIGRD with ther vertical grid ' &
            ,1)
    END IF

    CALL strhgrd_local(nz,strhopt,zrefsfc,z1,z2,z(nz-1),                      &
                 dzmin,strhtune, zp1d,dzp1d)

  ELSE


    WRITE(6,'(1x,a,i3,a/)')                                             &
        'Invalid vertical grid stretching option, strhopt was ',strhopt, &
        '. Program stopped in INIGRD.'
      CALL arpsstop('arpsstop called from INIGRD with stretching ' ,1)

  END IF
!
!
!-----------------------------------------------------------------------
!
!  Physical height of computational grid defined as
!
!  Zp=(z-zrefsfc)*(Zm-hterain)/(Zm-zrefsfc)+hterain for z=<Zm.
!  ZP=z for z>Zm
!
!  where Zm the height at which the grid levels becomes flat.
!  Hm < Zm =< Ztop, hm is the height of mountain and Ztop the height
!  of model top.
!
!-----------------------------------------------------------------------
!
  DO k=nz-1,2,-1
    IF(zp1d(k) <= zflat) THEN
      zflat1 = zp1d(k)
      EXIT
    END IF
  END DO
  zflat1=MAX(MIN(z(nz-1),zflat1),zrefsfc)

  DO k=2,nz-1

    IF(zp1d(k) > zflat1) THEN
      DO j=1,ny-1
        DO i=1,nx-1
          zp(i,j,k)=zp1d(k)
        END DO
      END DO
    ELSE
      DO j=1,ny-1
        DO i=1,nx-1
          zp(i,j,k)=(zp1d(k)-zrefsfc)*(zflat1-hterain(i,j))             &
                     /(zflat1-zrefsfc)+hterain(i,j)
        END DO
      END DO
    END IF

  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      zp(i,j,2)=hterain(i,j)
      zp(i,j,1)=2.0*zp(i,j,2)-zp(i,j,3)
      zp(i,j,nz)=2.0*zp(i,j,nz-1)-zp(i,j,nz-2)
    END DO
  END DO

END subroutine getVertCoordinate

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRHGRD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE strhgrd_local(nz,strhopt,z0,z1,z2,ztop,dzmin,strhtune, z,dzk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To construct a vertically stretched grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/17/94.
!
!  MODIFICATION HISTORY:
!
!  05/11/95 (Jinxing Zong and MX)
!
!  A bug fix for the case of nonzero zrefsfc, the reference height of
!  the surface. Results not affected for zrefsfc=0 (default value).
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nz       The vertical dimension of ARPS grid.
!
!
!  OUTPUT:
!
!    z        Array containing the height of veritical grid levels.
!    dzk      Array containing the grid spacing between vertical levels
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nz
  INTEGER :: strhopt
  REAL :: z0
  REAL :: z1
  REAL :: z2
  REAL :: ztop
  REAL :: dzmin
  REAL :: strhtune

  REAL :: z   (nz)
  REAL :: dzk (nz)

  REAL :: rnzh,dzm
  REAL :: a,b,c,hnew,zkmid,dzu
  INTEGER :: nzh,nzl,k
  REAL :: dz

  IF( (z1-z0) == (nz-3)*dzmin.AND.(z1-z0) == (ztop-z0) ) THEN

    dz = (ztop-z0)/(nz-3)
    DO k=1,nz-1
      dzk(k)= dz
    END DO
    DO k=1,nz
      z(k)=z0 + (k-2) * dz
    END DO

    WRITE(6,'(/1x,a,f13.3,/a,f13.3)')                                   &
        'Layer 1 depth was as deep as the entire domain. i',            &
        'A uniform vertical grid is assumed with dz=',dz,               &
        ' over the model depth of ',ztop-z0

    RETURN

  END IF

  IF(z1 < z0) z1 = z0
  IF(z2 > ztop) z2 = ztop

  nzl = (z1-z0)/dzmin

  IF( (z1-z0) >= (nz-3)*dzmin ) THEN
    WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                               &
        'Can not setup a vertical grid with uniform dz=',dzmin,         &
        ' over the depth of ',z1-z0,' please specify a smaller',        &
        ' value of dlayer1 '
      CALL arpsstop('arpsstop called from STRHGRD with stretching ' ,1)
  END IF

  IF( z2 >= ztop ) THEN
    dzm = (ztop-z0-nzl*dzmin)/(nz-3-nzl)
!    print*, nzl*dzmin + (nz-3-nzl)*dzm
    nzh = 0
    dzu = 2*dzm - dzmin
  ELSE
    a = 2*(nz-3-nzl)
    b = 2*z0-ztop-z2-(nz-3-3*nzl)*dzmin
    c = dzmin*(z2-z0-nzl*dzmin)
    dzm = (-b + SQRT(b*b-4*a*c) )/(2*a)

    rnzh = (ztop-z2)/(2*dzm-dzmin)
    nzh = INT(rnzh)

    hnew = nzl*dzmin + nzh*(2*dzm-dzmin) +                              &
          (nz-3-nzl-nzh)*dzm + z0

    IF( nzh /= 0 ) THEN
      dzu = (2*dzm-dzmin) + (ztop-hnew)/nzh
    ELSE IF( nz-3-nzl-nzh /= 0 ) THEN
      dzm = dzm + (ztop-hnew)/(nz-3-nzl-nzh)
      dzu = (2*dzm-dzmin)
    END IF

  END IF

  DO k=1,nzl+1
    dzk(k)=dzmin
  END DO


  IF( strhopt == 1 ) THEN

    a   = (dzm-dzmin)
    DO k=nzl+2,nz-2-nzh
      dzk(k)= dzm+a*                                                    &
          ((2.0*FLOAT(k-nzl-2)/FLOAT(nz-4-nzh-nzl)-1.0) )**3
    END DO

  ELSE

    zkmid=0.5*FLOAT( nz-nzh+nzl)

    IF( nzl+2-zkmid == 0.0 ) THEN
      b = 0.0
    ELSE
      b= strhtune* 2.0/(nzl+2-zkmid)
    END IF

    a=(dzmin-dzm)/TANH( strhtune* 2.0)

    DO k=nzl+2,nz-2-nzh
      dzk(k)=dzm + a*TANH(b*(FLOAT(k)-zkmid))
    END DO

  END IF

  DO k=nz-2-nzh+1, nz-2
    dzk(k)= dzu
  END DO

  dzk(nz-1) = dzk(nz-2)
  dzk(nz  ) = dzk(nz-1)


  z(2) = z0
  DO k=3,nz-1
    z(k) = z(k-1)+dzk(k-1)
  END DO

  z(1) =z(2)-dzk(1)
  z(nz-1)=ztop
  z(nz)=z(nz-1)+dzk(nz-1)

  RETURN
END SUBROUTINE strhgrd_local
!
Subroutine  OPEN_Mosaic(mosaicfile,NCID)
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

   CHARACTER*120 ::  mosaicfile
   integer :: NCID
   integer :: status
!
   STATUS=NF_OPEN(trim(mosaicfile),0,NCID)
!
   IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

end Subroutine OPEN_Mosaic
!
Subroutine  CLOSE_Mosaic(NCID)
! 
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

   integer :: NCID
   integer :: status
! 
   STATUS=NF_CLOSE(NCID)
! 
   IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

end Subroutine CLOSE_Mosaic


Subroutine  GET_Mosaic_sngl_Mosaic(NCID,mscNlon,mscNlat,mscNlev,mscValue)
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!  
!  First written: 04/06/2006.
!
!  IN:
!     mscNlon
!     mscNlan
!     mscNlev
!     NCID
!  out:
!     mscValue
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data
  INTEGER ::   mscNlev   ! number of level of mosaic data

  INTEGER ::  NCID, STATUS, MSID

  INTEGER ::   NDIMS
  PARAMETER (NDIMS=4)                  ! number of dimensions
  INTEGER START(NDIMS), COUNT(NDIMS)

  REAL ::   mscValue(mscNlon,mscNlat,1,1)
  INTEGER :: i,j

  START(1)=1
  START(2)=1
  START(3)=mscNlev
  START(4)=1
  COUNT(1)=mscNlon
  COUNT(2)=mscNlat
  COUNT(3)=1
  COUNT(4)=1

  STATUS = NF_INQ_VARID (NCID, 'mrefl_mosaic', MSID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_VARA_REAL (NCID, MSID, START, COUNT, mscValue)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

end subroutine GET_Mosaic_sngl_Mosaic


Subroutine  Check_DIM_ATT_Mosaic8(NCID,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat,mosaic_opt)
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!  
!  First written: 04/06/2006.
!
!  History:
!
!    05/15/2007 Fanyou Kong
!               Add 2D field capability
!
!  IN:
!     mscNlon
!     mscNlan
!     mscNlev
!     lonMin,latMin,lonMax,latMax,dlon,dlat
!
!  out:
!     NCID
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data
  INTEGER ::   mscNlev   ! number of level of mosaic data
  REAL    ::   lonMin,latMin,lonMax,latMax,dlon,dlat
  INTEGER :: mosaic_opt

  INTEGER ::  NCID, STATUS
  INTEGER ::  LONID, LATID, LEVID
  INTEGER ::  LONLEN, LATLEN, LEVLEN

  REAL ::   lonMinVal
  REAL ::   latMinVal
  REAL ::   lonMaxVal
  REAL ::   latMaxVal
  REAL ::   dlonVal
  REAL ::   dlatVal

  STATUS = NF_INQ_DIMID(NCID, 'Lon', LONID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMID(NCID, 'Lat', LATID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  IF(mosaic_opt == 1) THEN
  STATUS = NF_INQ_DIMID(NCID, 'Ht', LEVID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  END IF

  STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  IF(mosaic_opt == 1) THEN
  STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  ELSE
    LEVLEN=1
  END IF

  if( mscNlon.ne.LONLEN .or. mscNlat.ne.LATLEN .or. mscNlev.ne.LEVLEN) then
    write(*,*) 'Dimension is inconsistent'
    write(*,*) 'INPUT:', mscNlon,mscNlat,mscNlev
    write(*,*) 'Decoding', LONLEN,LATLEN,LEVLEN
    stop 123
  endif

  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'Longitude', lonMinVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'Latitude', latMaxVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'LonGridSpacing', dlonVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'LatGridSpacing', dlatVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  if( abs(lonMin-lonMinVal) > 0.00001 .or. &
      abs(latMax-latMaxVal) > 0.00001 .or. &
      abs(dlon-dlonVal) > 0.00001 .or. &
      abs(dlat-dlatVal) > 0.00001) then
    write(*,*) ' Attributes are inconsistent. Check'
    stop 123
  endif

END SUBROUTINE Check_DIM_ATT_Mosaic8

Subroutine  Check_DIM_ATT_Mosaic(NCID,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat)
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!  
!  First written: 04/06/2006.
!
!  IN:
!     mscNlon
!     mscNlan
!     mscNlev
!     lonMin,latMin,lonMax,latMax,dlon,dlat
!
!  out:
!     NCID
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data
  INTEGER ::   mscNlev   ! number of level of mosaic data
  REAL    ::   lonMin,latMin,lonMax,latMax,dlon,dlat

  INTEGER ::  NCID, STATUS
  INTEGER ::  LONID, LATID, LEVID
  INTEGER ::  LONLEN, LATLEN, LEVLEN

  REAL ::   lonMinVal
  REAL ::   latMinVal
  REAL ::   lonMaxVal
  REAL ::   latMaxVal
  REAL ::   dlonVal
  REAL ::   dlatVal

  STATUS = NF_INQ_DIMID(NCID, 'x', LONID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMID(NCID, 'y', LATID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMID(NCID, 'mrefl_mosaic_levels', LEVID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  if( mscNlon.ne.LONLEN .or. mscNlat.ne.LATLEN .or. mscNlev.ne.LEVLEN) then
    write(*,*) 'Dimension is inconsistent'
    write(*,*) 'INPUT:', mscNlon,mscNlat,mscNlev
    write(*,*) 'Decoding', LONLEN,LATLEN,LEVLEN
    stop 123
  endif

  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'xMin', lonMinVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'yMin', latMinVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'xMax', lonMaxVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'yMax', latMaxVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'dx', dlonVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'dy', dlatVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  if( abs(lonMin-lonMinVal) > 0.00001 .or. &
      abs(lonMax-lonMaxVal) > 0.00001 .or. &
      abs(latMin-latMinVal) > 0.00001 .or. &
      abs(latMax-latMaxVal) > 0.00001 .or. &
      abs(dlon-dlonVal) > 0.00001 .or. &
      abs(dlat-dlatVal) > 0.00001) then
    write(*,*) ' Attributes are inconsistent. Check'
    stop 123
  endif

END SUBROUTINE Check_DIM_ATT_Mosaic

Subroutine  GET_DIM_ATT_Mosaic8(mosaicsngle,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat,mosaic_opt)
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!  
!  First written: 04/06/2006.
!
!   New verison of Mosaic file that has 8 tiles
!
!  History:
!
!    05/15/2007 Fanyou Kong
!               Add 2D field capability
!
!  IN:
!     mosaicPath  : path of mosaic file
!     mosaicsngle : name of mosaic file
!  OUT
!     mscNlon
!     mscNlan
!     mscNlev
!     lonMin,latMin,lonMax,latMax,dlon,dlat
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  CHARACTER*120    mosaicsngle

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data
  INTEGER ::   mscNlev   ! number of level of mosaic data
  REAL    ::   lonMin,latMin,lonMax,latMax,dlon,dlat
  INTEGER ::   mosaic_opt

  INTEGER ::  NCID, STATUS
  INTEGER ::  LONID, LATID, LEVID
  INTEGER ::  LONLEN, LATLEN, LEVLEN

  REAL ::   lonMinVal
  REAL ::   latMinVal
  REAL ::   lonMaxVal
  REAL ::   latMaxVal
  REAL ::   dlonVal
  REAL ::   dlatVal

  STATUS = NF_OPEN(trim(mosaicsngle), 0, NCID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  STATUS = NF_INQ_DIMID(NCID, 'Lon', LONID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMID(NCID, 'Lat', LATID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  IF(mosaic_opt == 1) THEN
  STATUS = NF_INQ_DIMID(NCID, 'Ht', LEVID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  END IF

  STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  mscNlon=LONLEN
  mscNlat=LATLEN

  IF(mosaic_opt == 1) THEN
  STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  mscNlev=LEVLEN
  ELSE
  mscNlev=1
  END IF



  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'Longitude', lonMinVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'Latitude', latMaxVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'LatGridSpacing', dlonVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'LatGridSpacing', dlatVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  lonMin=lonMinVal
  latMax=latMaxVal
  dlon=dlonVal
  dlat=dlatVal
  lonMax=lonMinVal+dlon*mscNlon
  latMin=latMaxVal-dlat*mscNlat

  STATUS = NF_CLOSE(NCID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS) 

END SUBROUTINE GET_DIM_ATT_Mosaic8

Subroutine  GET_DIM_ATT_Mosaic(mosaicsngle,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat)
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!  
!  First written: 04/06/2006.
!
!  IN:
!     mosaicPath  : path of mosaic file
!     mosaicsngle : name of mosaic file
!  OUT
!     mscNlon
!     mscNlan
!     mscNlev
!     lonMin,latMin,lonMax,latMax,dlon,dlat
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  CHARACTER*120    mosaicsngle

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data
  INTEGER ::   mscNlev   ! number of level of mosaic data
  REAL    ::   lonMin,latMin,lonMax,latMax,dlon,dlat

  INTEGER ::  NCID, STATUS
  INTEGER ::  LONID, LATID, LEVID
  INTEGER ::  LONLEN, LATLEN, LEVLEN

  REAL ::   lonMinVal
  REAL ::   latMinVal
  REAL ::   lonMaxVal
  REAL ::   latMaxVal
  REAL ::   dlonVal
  REAL ::   dlatVal

  STATUS = NF_OPEN(trim(mosaicsngle), 0, NCID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  STATUS = NF_INQ_DIMID(NCID, 'x', LONID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMID(NCID, 'y', LATID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMID(NCID, 'mrefl_mosaic_levels', LEVID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  mscNlon=LONLEN
  mscNlat=LATLEN
  mscNlev=LEVLEN


  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'xMin', lonMinVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'yMin', latMinVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'xMax', lonMaxVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'yMax', latMaxVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'dx', dlonVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_ATT_REAL (NCID, NF_GLOBAL, 'dy', dlatVal)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

  lonMin=lonMinVal
  lonMax=lonMaxVal
  latMin=latMinVal
  latMax=latMaxVal
  dlon=dlonVal
  dlat=dlatVal

  STATUS = NF_CLOSE(NCID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS) 

END SUBROUTINE GET_DIM_ATT_Mosaic
!
SUBROUTINE HANDLE_ERR_Mosaic(STATUS)
     INCLUDE 'netcdf.inc'
     INTEGER STATUS
     IF (STATUS .NE. NF_NOERR) THEN
       PRINT *, NF_STRERROR(STATUS)
       STOP 'Stopped'
     ENDIF
END SUBROUTINE HANDLE_ERR_Mosaic
SUBROUTINE vert_interp_ref(nx,ny,nz,Nmsclvl,ref_mos_3d,ref_arps_3d,zp)
!
!  vertical interpolation NSSL mosica reflectiivty into arps grid
!
!  by MIng Hu (01/26/2007)
!
  implicit none

  INTEGER ::  nx,ny,nz
  INTEGER :: Nmsclvl
  real :: zp(nx,ny,nz)                  ! height in w grid
  real :: ref_arps_3d(nx,ny,nz)         ! reflectivity in arps grid
  real :: ref_mos_3d(nx,ny,Nmsclvl)     ! reflectivity in mosaic vertical grid

  real :: msclvlAll(31)
  real :: msclvl21(21),msclvl31(31)
  DATA msclvl21/1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7,  &
            8, 9, 10, 11, 12, 13, 14, 15, 16, 17/
  DATA msclvl31/0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, &
                3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, &
                9, 10, 11, 12, 13, 14, 15, 16, 18/
  REAL :: heightGSI,upref,downref,wght
  INTEGER :: ilvl,numref
  INTEGER :: istat_radar
  real :: zs(nz)                  ! height in maa grid

  INTEGER :: i,j, k2, k

    ref_arps_3d=-9999.0
    numref=0
    if (Nmsclvl == 31 ) then
      DO k=1,Nmsclvl
        msclvlAll(k)=msclvl31(k)*1000.0
      ENDDO
    elseif( Nmsclvl == 21 ) then
      msclvlAll=0
      DO k=1,Nmsclvl
        msclvlAll(k)=msclvl21(k)*1000.0
      ENDDO
    else
      write(*,*) ' Wrong vertical radar mosaic levels'
      stop 123
    endif
    DO j=1,ny
    DO i=1,nx
      DO k2=2,nz-1
        zs(k2)=(zp(i,j,k2)+zp(i,j,k2+1))/2
      enddo
      zs(1)=zp(i,j,1)
      zs(nz)=zp(i,j,nz)
    DO k2=1,nz
      heightGSI=zs(k2)
!      heightGSI=zp(i,j,k2)
      if(heightGSI >= msclvlAll(1) .and. heightGSI < msclvlAll(Nmsclvl) ) then
         do k=1,Nmsclvl-1
          if( heightGSI >=msclvlAll(k) .and. heightGSI < msclvlAll(k+1) ) ilvl=k
         enddo
         upref=ref_mos_3d(i,j,ilvl+1)
         downref=ref_mos_3d(i,j,ilvl)
         if(abs(upref) < 100.0 .and. abs(downref) <100.0 ) then
           wght=(heightGSI-msclvlAll(ilvl))/(msclvlAll(ilvl+1)-msclvlAll(ilvl))
           ref_arps_3d(i,j,k2)=(1-wght)*downref + wght*upref
           numref=numref+1
         else
           ref_arps_3d(i,j,k2)=-9999.0
         endif
      else
        ref_arps_3d(i,j,k2)=-9999.0
      endif
    ENDDO !nz
    ENDDO
    ENDDO
    if(numref>0) istat_radar=1

END SUBROUTINE vert_interp_ref
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE WTRADCOL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wtradcol_mosaic(nx,ny,nz,rfname,            &
           iyr,imon,iday,ihr,imin,isec,         &
           grdlatc,grdlonc,zpsc,gridref) 
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Writes gridded radar data to a file as columns with
!  individual lat,lons.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!
!  01/26/07 ming HU
!   for Mosaic file
!-----------------------------------------------------------------------
!
!  INPUT:
!    dmpfmt    file format (1:binary, 2:hdf)
!    iradfmt   binary format
!    hdf4cmpr  hdf4 compression level
!    rfname    radar file name (character*80)
!    radid     radar id (character*4)
!    latrad    latitude of radar (degrees N)
!    lonrad    longitude of radar (degrees E)
!    elvrad    elevation of radar (m MSL)
!    iyr       year
!    imon      month
!    iday      day
!    ihr       hour
!    imin      min
!    isec      sec
!    vcpnum    VCP (scan type) number
!    isource)  source number
!                1= WSR-88D raw
!                2= WSR-88D NIDS
!
!  OUTPUT:
!    data are written to file
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
!
  INTEGER :: dmpfmt
  INTEGER :: iradfmt
  INTEGER :: hdf4cmpr
  CHARACTER (LEN=*) :: rfname
  CHARACTER (LEN=4) :: radid
  INTEGER :: iyr,imon,iday,ihr,imin,isec
  INTEGER :: isource
!
  REAL :: zs(nz)
  REAL :: zpsc(nx,ny,nz)
  REAL :: gridref(nx,ny,nz)
!
  REAL :: outk(nz)
  REAL :: outhgt(nz)
  REAL :: outref(nz)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Radar output descriptors
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Radar output thresholds
!
!-----------------------------------------------------------------------
!
  REAL :: refmin,refmax
  PARAMETER(refmin=-5.0, refmax=100)
  REAL :: misval
  PARAMETER(misval=-999.0)
!
!-----------------------------------------------------------------------
!
!  Radar output variables
!
!-----------------------------------------------------------------------
!
  REAL :: grdlatc(nx,ny)
  REAL :: grdlonc(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunit,myr,itime
  INTEGER :: i,j,k,klev,kk,kntcol,nn
  INTEGER :: idummy
  INTEGER :: istat,sd_id
  INTEGER :: irngmin,irngmax
  REAL :: gridlat,gridlon,elev,rdummy
!  INTEGER(2), allocatable :: itmp(:,:,:) ! Temporary array
!  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=80)  :: runname
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  print *, ' inside wtradcol_mosaic'
  idummy=-999 
  rdummy=-999.
  dmpfmt=1
  iradfmt=1
  hdf4cmpr=1
  radid='MOSC'
  runname='Mosaic2arps'
!
  IF (dmpfmt > 1 .AND. hdf4cmpr > 3) THEN
     write(*,*) 'HDF is not available!'
     stop 123
  END IF 

  if( iyr < 100 ) then
     myr=1900+iyr
  else
     myr=iyr
  endif
  IF(myr < 1960) myr=myr+100
  CALL ctim2abss(myr,imon,iday,ihr,imin,isec,itime)

  IF(dmpfmt == 1)THEN
  CALL getunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
  filename=TRIM(rfname)//'.MscARPS'
!  write(*,*) filename
  OPEN(iunit,FILE=TRIM(filename),STATUS='UNKNOWN',FORM='UNFORMATTED')
!
!-----------------------------------------------------------------------
!
!  Write radar description variables
!
!-----------------------------------------------------------------------
!
  WRITE(iunit) radid
  WRITE(iunit) idummy,itime,idummy,isource,idummy,                     &
               idummy,idummy,idummy,idummy,idummy
!
!-----------------------------------------------------------------------
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  icluding elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!
!-----------------------------------------------------------------------
!
  idummy=0
  rdummy=0.
  WRITE(iunit) runname
  WRITE(iunit) iradfmt,strhopt,mapproj,idummy,idummy,                 &
               idummy,idummy,idummy,idummy,idummy
  WRITE(iunit) dx,dy,dz,dzmin,ctrlat,                                   &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               rdummy,rdummy,rdummy,rdummy,rdummy
  WRITE(iunit) idummy,idummy
  ELSE  !HDF4 format
    !
   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF
!
!-----------------------------------------------------------------------
!
!  For each horizontal grid point form a column of remapped
!  data containing the non-missing grid points
!
!-----------------------------------------------------------------------
!
  IF(dmpfmt==1)THEN
    kntcol=0
    DO j=1,ny
      DO i=1,nx
        DO k=1,nz
          outk(k)=misval
          outhgt(k)=misval
          outref(k)=misval
        END DO
        klev=0
        DO k=2,nz-1
          zs(k)=(zpsc(i,j,k)+zpsc(i,j,k+1))/2
        enddo
        zs(1)=zpsc(i,j,1)
        zs(nz)=zpsc(i,j,nz)
        DO k=1,nz-1
          IF(gridref(i,j,k)>refmin .AND. gridref(i,j,k)<refmax) THEN
            klev=klev+1
            outk(klev)=FLOAT(k)
            outhgt(klev)=zs(k)
            outref(klev)=gridref(i,j,k)
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  If there are data in this column, write them to the file.
!
!-----------------------------------------------------------------------
!
        IF(klev > 0) THEN
          kntcol=kntcol+1
          elev=zpsc(i,j,2)
          WRITE(iunit) i,j,rdummy,rdummy,                &
                       grdlatc(i,j),grdlonc(i,j),elev,klev
          WRITE(iunit) (outk(k),k=1,klev)
          WRITE(iunit) (outhgt(k),k=1,klev)
          WRITE(iunit) (outref(k),k=1,klev)
       END IF
      END DO
    END DO
  ELSE    !HDF4 format
!  
   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF
!
  IF(dmpfmt==1)THEN
    CLOSE(iunit)
    CALL retunit(iunit)
  ELSE
  !  
   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Report on what data were written
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(//a,i4.4,i2.2,i2.2,a1,i2.2,a1,i2.2)')                       &
                    ' Output statistics for time ',                     &
                      iyr,imon,iday,' ',ihr,':',imin
  WRITE(6,'(a,i6,a,/a,i6,a//)')                                         &
           ' There were ',kntcol,' columns written ',                   &
           ' of a total ',(nx*ny),' possible.'
!
  RETURN
END SUBROUTINE wtradcol_mosaic
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE READADCOL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readadcol_Mosaic(nx,ny,nz,rfname,gridref)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  read in gridded radar data to a file as columns with
!  individual lat,lons.
!
!-----------------------------------------------------------------------
!
!  01/28/06 MIng HU
!  Modified to fit the need to write NSSL mosaic reflectiivty 
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    dmpfmt    file format (1:binary, 2:hdf)
!    iradfmt   binary format
!    hdf4cmpr  hdf4 compression level
!    rfname    radar file name (character*80)
!    iyr       year
!    imon      month
!    iday      day
!    ihr       hour
!    imin      min
!    isec      sec
!    isource)  source number
!                1= WSR-88D raw
!                2= WSR-88D NIDS
!                3= NSSL mosaic reflectivity
!
!  OUTPUT:
!    data are written to file
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'grid.inc'
!
  INTEGER :: nx,ny,nz
!
  INTEGER :: dmpfmt
  INTEGER :: iradfmt
  INTEGER :: hdf4cmpr
  CHARACTER (LEN=*) :: rfname
  INTEGER :: iyr,imon,iday,ihr,imin,isec
  INTEGER :: isource
!
  REAL :: zs(nz)
  REAL :: zpsc(nx,ny,nz)
  REAL :: gridref(nx,ny,nz)
!
  REAL :: readk(nz)
  REAL :: readhgt(nz)
  REAL :: readref(nz)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
!  Radar output descriptors
!
!-----------------------------------------------------------------------
!
!  INTEGER :: mxradvr,nradvr
!  PARAMETER(mxradvr=10,nradvr=6)
!  INTEGER :: iradvr(mxradvr)
!  DATA iradvr /1,2,3,4,5,6,0,0,0,0/
!
!-----------------------------------------------------------------------
!
!  Radar output thresholds
!
!-----------------------------------------------------------------------
!
  REAL :: refmin,refmax,velmin,velmax
  PARAMETER(refmin=-5.0, refmax=100.,                                   &
            velmin=-200.,velmax=200.)
  REAL :: misval
  PARAMETER(misval=-999.0)
!
!-----------------------------------------------------------------------
!
!  Radar output variables
!
!-----------------------------------------------------------------------
!
  REAL :: grdlatc(nx,ny)
  REAL :: grdlonc(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4) :: radid
  CHARACTER (LEN=12) :: runname

  INTEGER :: iunit,myr,itime
  INTEGER :: i,j,k,klev,kk,kntcol,nn
  INTEGER :: idummy
  INTEGER :: istat,sd_id
  INTEGER :: ipt
  REAL :: gridlat,gridlon,elev,rdummy
  INTEGER(2), allocatable :: itmp(:,:,:) ! Temporary array
  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  CHARACTER (LEN=256) :: filename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  print *, ' inside readadcol_mosaic'
  idummy=-999 
  rdummy=-999.
  dmpfmt=1
  iradfmt=1
  hdf4cmpr=1
!
  IF (dmpfmt > 1 .AND. hdf4cmpr > 3) THEN
   write(*,*) 'HDF is not available!'
   stop 123
  END IF 

  IF(dmpfmt == 1)THEN
  CALL getunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
  filename=TRIM(rfname)//'.MscARPS'
  write(*,*) filename
  OPEN(iunit,FILE=TRIM(filename),STATUS='UNKNOWN',FORM='UNFORMATTED')
!
!-----------------------------------------------------------------------
!
!  Write radar description variables
!
!-----------------------------------------------------------------------
!
  read(iunit) radid
  read(iunit) idummy,itime,idummy,isource,idummy,                     &
               idummy,idummy,idummy,idummy,idummy
!
  CALL abss2ctim(itime,myr,imon,iday,ihr,imin,isec)
!  write(*,*) 'The time of the data is: ',myr,imon,iday,ihr,imin,isec
!-----------------------------------------------------------------------
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  icluding elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!
!-----------------------------------------------------------------------
!
  idummy=0
  rdummy=0.
  read(iunit) runname
  read(iunit) iradfmt,strhopt,mapproj,idummy,idummy,                   &
               idummy,idummy,idummy,idummy,idummy
  read(iunit) dx,dy,dz,dzmin,ctrlat,                                   &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               rdummy,rdummy,rdummy,rdummy,rdummy
  read(iunit) idummy,idummy
  ELSE  !HDF4 format
!
   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF
!
!-----------------------------------------------------------------------
!
!  For each horizontal grid point form a column of remapped
!  data containing the non-missing grid points
!
!-----------------------------------------------------------------------
!
  IF(dmpfmt==1)THEN

     DO ipt=1,(nx*ny)

       read(iunit,END=51) i,j,rdummy,rdummy,                    &
                   gridlat,gridlon,elev,klev
       read(iunit,END=51) (readk(k),k=1,klev)
       read(iunit,END=51) (readhgt(k),k=1,klev)
       read(iunit,END=51) (readref(k),k=1,klev)
 
       IF(i <= nx.AND.i >= 1 .AND. j <= ny.AND.j >= 1) THEN
          DO kk=1,klev
             k=nint(readk(kk))
             IF(k <= nz.AND.k >= 1) THEN
                gridref(i,j,k)=readref(kk)
             END IF  ! 1 < k < nz
          END DO  ! kk = 1, klev
       END IF  ! 1 < i < nx  & 1 < j < ny

     END DO  ! ipt = 1, nx*ny
 51  continue
     ipt=ipt-1
     WRITE(6,'(a,i6,a)') ' End of file reached after reading',          &
                       ipt,' columns'


  ELSE    !HDF4 format
!
   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF
!
  RETURN
END SUBROUTINE readadcol_Mosaic

Subroutine  GET_Mosaic_cref_Mosaic(NCID,mscNlon,mscNlat,mscValue)
!
!  Read in 2D field CREF
!
!  Author: Fanyou Kong, CAPS. University of Oklahma.
!  First written: 05/15/2007.
!
!  IN:
!     mscNlon
!     mscNlan
!     NCID
!  out:
!     mscValue
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data

  INTEGER ::  NCID, STATUS, MSID

  INTEGER ::   NDIMS
  PARAMETER (NDIMS=2)                  ! number of dimensions
  INTEGER START(NDIMS), COUNT(NDIMS)

  REAL ::   mscValue(mscNlon,mscNlat)
  INTEGER :: i,j

  START(1)=1
  START(2)=1
  COUNT(1)=mscNlon
  COUNT(2)=mscNlat

  STATUS = NF_INQ_VARID (NCID, 'cref', MSID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_VARA_REAL (NCID, MSID, START, COUNT, mscValue)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

end subroutine GET_Mosaic_cref_Mosaic

Subroutine  GET_Mosaic_hsr_Mosaic(NCID,mscNlon,mscNlat,mscValue)
!
!  Read in 2D field 1h_rad_hsr
!
!  Author: Fanyou Kong, CAPS. University of Oklahma.
!  First written: 05/15/2007.
!
!  IN:
!     mscNlon
!     mscNlan
!     NCID
!  out:
!     mscValue
!
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data

  INTEGER ::  NCID, STATUS, MSID

  INTEGER ::   NDIMS
  PARAMETER (NDIMS=2)                  ! number of dimensions
  INTEGER START(NDIMS), COUNT(NDIMS)

  REAL ::   mscValue(mscNlon,mscNlat)
  INTEGER :: i,j

  START(1)=1
  START(2)=1
  COUNT(1)=mscNlon
  COUNT(2)=mscNlat

  STATUS = NF_INQ_VARID (NCID, 'rad_hsr_1h', MSID)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)
  STATUS = NF_GET_VARA_REAL (NCID, MSID, START, COUNT, mscValue)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR_Mosaic(STATUS)

end subroutine GET_Mosaic_hsr_Mosaic
