PROGRAM arpstrn
!
!-----------------------------------------------------------------------
!
!  Program for creating a terrain file for ARPS from
!  USGS 3 arc second data base (covers continetal US and Alaska)
!  or 5min global data base from NCAR.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue, CAPS. University of Oklahma.
!
!  First written: 11/25/1995.
!  Last Modified: 12/05/1995.
!
!  MODIFICATIONS:
!
!  11/1/1998 (Ming Xue).
!  Significantly rewritten to also handle 5min global data.
!
!  The code for accessing USGS FTP server and reading the DEM data
!  patches were written by Burt Freeman and P. Patnaik of S-Cubed.
!
!  2001/06/18 (Gene Bassett)
!  Cleaned up namelist input of grid and map projections parameters.
!
!  2004/09/29 (M. Xue)
!  For 2nd-order interpolation, using NINT instead of INT to find the
!  nearest grid point to interploation.
!
!  2012/09/11 (Y. Wang)
!  Added support for SRTM 90m (~3 second) data sets, both Version 2.1 and
!  Version 4.
!
!-----------------------------------------------------------------------
!
!
!  INSTRUCTIONS:
!
!  A number functions are controled by Makfile.
!
!  To compile, do 'make arpstrn' in ARPS root directory.
!
!  You can also do the compilation and linking directly in this
!  directory, making use of pre-built arps library:
!
!  zxpostf90 -I../../include -o arpstrn arpstrn.f90 -L../../lib -larps
!
!  Here we assumed libarps.a has been produced by 'makearps arps'.
!  We also assume you have ZXPLOT Version 3.0 or later installed.
!  If not, do
!
!  f90 -I../../include -o arpstrn arpstrn.f nozxplot.f -L../../lib -larps
!
!  then no plot of the analyzed terrain field will be produced.
!
!  To execute, do 'arpstrn < arpstrn.input >! arpstrn.output'
!
!  When you use the 3 arc second data, a temporary directory
!  whose name is specified for dir_trndata in the input file
!  will be used to store the original data. The required
!  data will be downloaded from anonymous ftp server
!  caps.ou.edu unless you already have them in that directory.
!
!  To allow automatic login, you need to add
!
!  machine caps.ou.edu login anonymous password your_address
!
!  in your $HOME/.netrc file. Make .netrc readable to the owner only
!  i.e., do 'chmod 600 .netrc'.
!
!  To produce a tar file of the entire package, do 'make arpstrn.tar'.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'

  REAL, ALLOCATABLE :: x(:)         ! The x-coord. of the u grid point (m)
  REAL, ALLOCATABLE :: y(:)         ! The y-coord. of the v grid point (m)

  REAL, ALLOCATABLE :: xh(:,:)      ! X-coord. of the terrain (scalar) point.
  REAL, ALLOCATABLE :: yh(:,:)      ! Y-coord. of the terrain (scalar) point.
  REAL, ALLOCATABLE :: xh1d(:)
  REAL, ALLOCATABLE :: yh1d(:)

  REAL, ALLOCATABLE :: lath(:,:)    ! Latitude of each terrain point
  REAL, ALLOCATABLE :: lonh(:,:)    ! Longitude of each terrain point

  REAL, ALLOCATABLE :: hterain(:,:) ! The height of the terrain.
  REAL, ALLOCATABLE :: tem  (:,:)
  REAL, ALLOCATABLE :: tem2 (:,:)
  REAL, ALLOCATABLE :: xw   (:)
  REAL, ALLOCATABLE :: yw   (:)

  INTEGER :: nx_trn,ny_trn
  REAL, ALLOCATABLE :: h_trn(:,:)

!-----------------------------------------------------------------------
!
! Namelist definitions
!
!-----------------------------------------------------------------------

  NAMELIST /jobname/ runname

  INTEGER :: nx, ny

  NAMELIST /grid/ nx,ny,dx,dy,ctrlat,ctrlon
  NAMELIST /projection/ mapproj,trulat1,trulat2,trulon,sclfct

  INTEGER :: trnanxopt     ! Terrain analysis scheme option
  INTEGER :: trndataopt    ! Terrain data set option
  INTEGER :: nsmth         ! Number of times 9-point smoother is applied
                           ! to the terrain field
  CHARACTER (LEN=256) :: usgs_dem_index_fn, fn_trndata, dir_trndata
  INTEGER :: lat_sample,lon_sample    ! Data sampling frequency in seconds
  REAL    :: latd_sample,lond_sample  ! Data sampling frequency in degrees

  INTEGER :: trn3src

  NAMELIST /dem_trn/ trnanxopt,trndataopt,nsmth,                        &
            lat_sample,lon_sample,dir_trndata,fn_trndata,               &
            usgs_dem_index_fn,trn3src

  INTEGER :: ldir_trndata,lfn,lenstr

  REAL    :: latgrid,longrid
  INTEGER :: lmapfile,nmapfile
  INTEGER, PARAMETER :: maxmap = 10
  CHARACTER(LEN=256) :: mapfile(maxmap)
  REAL    :: hctrinc

  NAMELIST /map_plot/latgrid,longrid,nmapfile,mapfile,hctrinc
  COMMON /mappar2/ latgrid,longrid

  NAMELIST /trn_output/ dirname,terndmp

  INTEGER :: trndatares    ! Terrain data resolution in seconds
!
!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
  REAL :: trulat (2)           ! Latitude at which the map projection
                               ! is true (degees N)

  REAL :: swx,swy,ctrx,ctry    ! Temporary variables.

  REAL :: latmax,latmin,lonmax,lonmin ! max/min lat/lon of the projected grid
  INTEGER :: swlat_trn,swlon_trn,nelat_trn,nelon_trn
                               ! SW and NE coner lat/lon of the
                               ! the assembled terrain patch.
  INTEGER :: i,j,ih,jh

  REAL :: lath1,lonh1,w1,w2,w3,w4,a,b,c,d,h1,h2,h3
  REAL :: badvalue

  REAL :: hmax,hmin,swlat_offset,swlon_offset
  INTEGER :: inith_set, badvfound
  LOGICAL :: fexist
  INTEGER :: imin,imax,jmin,jmax,kmin,kmax, istatus
!
!-----------------------------------------------------------------------
!
!  TEMPORARY ARRAYS:
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=256) :: outfilename


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL mpinit_var          ! calling of writtrn needs some variables
!
!-----------------------------------------------------------------------
!
!  Set up the default values for all the variables to be read in
!  using the namelist method. In case the user does not specify a
!  particular value, this value will be used.
!
!-----------------------------------------------------------------------
!
  nx = 67
  ny = 67

  runname = 'may20'

  dx = 1000.0
  dy = 1000.0

  ctrlat  =  35.0
  ctrlon  = -100.0

  mapproj  = 0
  trulat1  =   30.0
  trulat2  =   60.0
  trulon   = -100.0
  sclfct   =    1.0

  trnanxopt  = 1
  trndataopt = 4
  fn_trndata = 'tbase_global_5min.data'
  dir_trndata= '/work/official/arpstern.data/'
  usgs_dem_index_fn = '../src/arpstrn/usgs_dem.index'
  nsmth      = 1

  lat_sample = 30
  lon_sample = 30

  trn3src = 0               !TINA set default

  hctrinc = 250.0

  latgrid =  10.
  longrid =  10.
  nmapfile = 1
  mapfile(1) = '/home/mxue/zxplot3/data/us_state.mapdata'
  mapfile(2) = '/home/mxue/zxplot3/data/world_us_country.mapdata'

  terndmp = 1

!
!-----------------------------------------------------------------------
!
!  Read and Print of Namelist variables and parameters.
!
!-----------------------------------------------------------------------
!

  READ (5,jobname,ERR=980,END=990)
  WRITE(6,'(/a,a)') 'The name of this run is: ', runname

  CALL gtlfnkey( runname, lfnkey )

  READ(5,grid,ERR=980,END=990)
  WRITE(6,'(/a,f10.3,a)')'Input dx was ',dx,' meters'
  WRITE(6,'(/a,f10.3,a)')'Input dy was ',dy,' meters'

  allocate(x(nx))
  allocate(y(ny))

  allocate(xh(nx,ny))
  allocate(yh(nx,ny))
  allocate(xh1d(0:nx))
  allocate(yh1d(0:ny))

  allocate(lath(0:nx,0:ny))
  allocate(lonh(0:nx,0:ny))

  allocate(hterain(nx,ny))
  allocate(tem  (nx,ny))
  allocate(tem2 (nx,ny))
  allocate(xw   (8*nx))
  allocate(yw   (8*nx))

  READ(5,projection,ERR=980,END=990)

  WRITE(6,'(/a,i4)') 'Input mapproj was ',mapproj
  WRITE(6,'(/a,f10.3,a)')                                               &
       'Input trulat1 was ',trulat1,' degree North'
  WRITE(6,'(/a,f10.3,a)')                                               &
       'Input trulat2 was ',trulat2,' degree North'
  WRITE(6,'(/a,f10.3)')                                                 &
      'The latitude of the center of the model domain was ',ctrlat
  WRITE(6,'(/a,f10.3)')                                                 &
      'The longitude of the center of the model domain was ',ctrlon
  WRITE(6,'(/a,f10.3,a)')                                               &
       'Input trulon was ',trulon,' degree East'
  WRITE(6,'(/a,e15.5)')                                                 &
       'Input sclfct was ',sclfct

  IF ( mapproj == 0 ) THEN
    trulat1 = ctrlat
    trulat2 = ctrlat
    trulon  = ctrlon
  END IF

  trnanxopt = 1
  trndataopt= 2
  nsmth = 0

  fn_trndata = 'tbase_global_5min.data'
  usgs_dem_index_fn = 'usgs_dem.index'

  READ(5,dem_trn,ERR=980,END=990)

  IF( trndataopt == 1) THEN
    trndatares = 3600
    lon_sample = max(3600,lon_sample)
    lat_sample = max(3600,lat_sample)
  ELSE IF( trndataopt == 2) THEN
    trndatares = 300
    lon_sample = max(300,lon_sample)
    lat_sample = max(300,lat_sample)
  ELSE IF( trndataopt == 3) THEN
    trndatares = 30
    lon_sample = max(30,lon_sample)
    lat_sample = max(30,lat_sample)
  ELSE IF( trndataopt == 4) THEN
    trndatares = 3
    lon_sample = max(3,lon_sample)
    lat_sample = max(3,lat_sample)

    IF (trn3src == 1) THEN      !TINA - read in option to use Europe data
      WRITE(6,'(/a/)') ' Will read in 3s Eropean terrain data set.'
    ELSE IF (trn3src == 2) THEN
      WRITE(6,'(/a/)') ' Will read in 90m SRTM Eurasia terrain data.'
    ELSE IF (trn3src == 3) THEN
      WRITE(6,'(/a/)') ' Will read in 90m SRTM global terrain data.'
    ELSE
      WRITE(6,'(/a/)') ' Will read in 3s US terrain data set.'
    END IF

  END IF

  WRITE(6,'(/a,i5,a)')                                                  &
      'The sample internal in longitude dir. was set to ',lon_sample,   &
      ' seconds.'

  WRITE(6,'(/a,i5,a)')                                                  &
      'The sample internal in latitude dir. was set to ',lat_sample,    &
      ' seconds.'

  IF(MOD(lat_sample,trndatares) /= 0 .OR. MOD(lon_sample,trndatares) /= 0) THEN
    lon_sample = MAX(1,lon_sample/trndatares) * trndatares
    lat_sample = MAX(1,lat_sample/trndatares) * trndatares
    WRITE(6,'(/a,i4,a,/a,i6,a,i6/))')                                   &
        'lat_sample and lon_sample have to be divisible by ',           &
        trndatares,'.',                                                 &
        'They were reset to lon_sample=',lon_sample,                    &
        ' lat_sample=',lat_sample
  END IF

  IF(MOD(3600,lat_sample) /= 0 .OR. MOD(3600,lon_sample) /= 0) THEN
    WRITE(6,'(/a,/a/))')                                                &
        '3600 has to be divisible by both lat_sample and lon_sample.',  &
        'Please reset them in the input file. Job stopped.'
    STOP
  END IF

  latgrid = 0.0
  longrid = 0.0
  READ(5,map_plot)

  dirname = './'
  READ(5,trn_output)

  WRITE(6,'(/a,i4)') 'Terrain dumping option set to ',terndmp

  lenstr = LEN_TRIM(dirname)
  IF(lenstr > 0) THEN
    IF(dirname(lenstr:lenstr) /= '/') THEN
      dirname(lenstr+1:lenstr+1) = '/'
      lenstr = lenstr + 1
    END IF
  ELSE
    dirname = './'
  END IF

  WRITE(6,'(/a/)')'Finished reading namelist input file.'

  ldir_trndata = 256
  CALL strlnth( dir_trndata, ldir_trndata)
  IF( dir_trndata(ldir_trndata:ldir_trndata) /= '/') THEN
    ldir_trndata=ldir_trndata+1
    dir_trndata(ldir_trndata:ldir_trndata)='/'
  END IF
!
!-----------------------------------------------------------------------
!
!  setting ctrlon in terms of degrees east
!
!-----------------------------------------------------------------------
!
!  IF(ctrlon.lt.0)THEN
!    ctrlon=360+ctrlon
!  END IF
!
!-----------------------------------------------------------------------
!
!  test to see if reset value is within longitudinal limits..
!
!-----------------------------------------------------------------------
!
!  IF(ctrlon.lt.0.or.ctrlon.gt.360)THEN
!    print *,'ctrlon is incorrectly set in terrain.input, please
!    :           reset, use degrees east'
!    STOP
!  END IF
!
!-----------------------------------------------------------------------
!
!  test to see if the ctrlat is within limits...
!
!-----------------------------------------------------------------------
!

  IF( ctrlat > 90 .OR. ctrlat < -90. .OR.                               &
      (trndataopt ==4 .AND. ctrlat < 0.0) )THEN
    PRINT *,'ctrlat in terrain.input is too large or too small. ',      &
            'must be set between 0. AND +90. degree north.'
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  setting trulon in terms of degrees east
!
!-----------------------------------------------------------------------
!
!  IF(trulon.lt.0)THEN
!    trulon=360+trulon
!  END IF
!
!-----------------------------------------------------------------------
!
!  test to see if reset value is within longitudinal limits..
!
!-----------------------------------------------------------------------
!
!  IF(trulon.lt.0.or.trulon.gt.360)THEN
!    print *,'trulon is incorrectly set in terrain.input, please
!    &           reset, use degrees east'
!    STOP
!  END IF
!
!-----------------------------------------------------------------------
!
!  test to see if trulat is within limits...
!
!-----------------------------------------------------------------------
!

  IF ( trulat1 > 90 .OR. trulat1 < -90 .OR.                             &
       (trndataopt == 4 .AND. trulat1 < 0.0) )THEN
    PRINT *,'Input parameter trulat1 is too large or too small.',       &
            'must be set between 0. AND +90. degree north.'
    STOP
  END IF

  IF ( trulat2 > 90 .OR. trulat2 < -90  .OR.                            &
       (trndataopt == 4 .AND. trulat2 < 0.0) )THEN
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

  WRITE(*,'(1x,a,2F10.2)') 'dx,dy = ',dx,dy
!
!-----------------------------------------------------------------------
!
!  Final out the lat/lon of each scalar (terrain data) point.
!
!-----------------------------------------------------------------------
!
!  Include an extra zone so that data grid covering
!  (lonmax-lonmin)*(latmax-latmin) is large enough area for 2nd order
!   interpolation
!
  CALL xytoll(nx+1,ny+1,xh1d,yh1d,lath,lonh)
!
!-----------------------------------------------------------------------
!
!  Final out the max/min lat/lon of each scalar (terrain data) point.
!
!-----------------------------------------------------------------------


  CALL a3dmax0(lath,0,nx,0,nx,0,ny,0,ny,1,1,1,1,latmax,latmin)
  CALL a3dmax0(lonh,0,nx,0,nx,0,ny,0,ny,1,1,1,1,lonmax,lonmin)

  IF( trndataopt == 3 ) THEN
    swlat_trn = INT((1000.0+latmin-15.0/3600.0)*0.1)*10-1000
    swlon_trn = INT((1000.0+lonmin-15.0/3600.0)*0.1)*10-1000
    nelat_trn = INT((1000.0+latmax+15.0/3600.0)*0.1)*10-1000
    nelon_trn = INT((1000.0+lonmax+15.0/3600.0)*0.1)*10-1000
!   IF( nelat_trn /= latmax ) nelat_trn = nelat_trn + 10
!   IF( nelon_trn /= lonmax ) nelon_trn = nelon_trn + 10
    IF( nelat_trn <  latmax+15.0/3600.0 ) nelat_trn = nelat_trn + 10
    IF( nelon_trn <  lonmax+15.0/3600.0 ) nelon_trn = nelon_trn + 10
    swlat_offset = 15.0/3600.0
    swlon_offset = 15.0/3600.0
  ELSE
    swlat_trn = INT(1000.0+latmin)-1000
    swlon_trn = INT(1000.0+lonmin)-1000
    nelat_trn = INT(1000.0+latmax)-1000
    nelon_trn = INT(1000.0+lonmax)-1000
    IF( nelat_trn /= latmax ) nelat_trn = nelat_trn + 1
    IF( nelon_trn /= lonmax ) nelon_trn = nelon_trn + 1
    swlat_offset = 0.0
    swlon_offset = 0.0
  ENDIF


  swlat_trn=MAX( -90, swlat_trn)
  nelat_trn=MIN(  90, nelat_trn)
  swlon_trn=MAX(-180, swlon_trn)
  nelon_trn=MIN( 180, nelon_trn)

  WRITE(*,'(1x,a,2(F10.5,a,I4,a))') 'lonmin, lonmax=',lonmin, '(',swlon_trn,'),', lonmax,'(',nelon_trn,')'
  WRITE(*,'(1x,a,2(F10.5,a,I4,a))') 'latmin, latmax=',latmin, '(',swlat_trn,'),', latmax,'(',nelat_trn,')'
  WRITE(*,*)
  WRITE(*,'(1x,2(a,I3,a,I4),a)') 'sw lat/lon = (',swlat_trn,',',swlon_trn, &
                              '); ne lat/lon = (',nelat_trn,',',nelon_trn,')'

  IF( trndataopt == 3 ) THEN
    nx_trn = (3600*(nelon_trn-swlon_trn)-30.0)/lon_sample+1
    ny_trn = (3600*(nelat_trn-swlat_trn)-30.0)/lat_sample+1
  ELSE
    nx_trn = (3600*(nelon_trn-swlon_trn))/lon_sample+1
    ny_trn = (3600*(nelat_trn-swlat_trn))/lat_sample+1
  ENDIF

  Write(6,'(/1x,2(a,i7),a/)')                                           &
                'Allocating array h_trn of size ',nx_trn,' x',ny_trn,'.'

  ALLOCATE(h_trn(nx_trn,ny_trn),STAT = istatus)
  CALL check_alloc_status(istatus, "arpstrn:h_trn")

  badvalue = -9999.0   !  Default value of -9999.0
  DO i = 1,nx_trn
    DO j = 1,ny_trn
      h_trn(i,j)=badvalue
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Fill in array h_trn, by reading in 3 arc second DEM/terrain data.
!
!-----------------------------------------------------------------------

  IF( trndataopt == 1) THEN

!    CALL FILLHTRN1DEG(h_trn,nx_trn,ny_trn, &
!         swlon_trn,nelon_trn,lon_sample,swlat_trn,nelat_trn,lat_sample,&
!         dir_trndata(1:ldir_trndata),badvalue)

  ELSE IF( trndataopt == 2) THEN

    CALL fillhtrn5min(h_trn,nx_trn,ny_trn, &
        swlon_trn,nelon_trn,lon_sample,swlat_trn,nelat_trn,lat_sample,  &
        dir_trndata(1:ldir_trndata),fn_trndata,badvalue)

  ELSE IF( trndataopt == 3) THEN

     PRINT*,'nx_trn,ny_trn = ',nx_trn,ny_trn

     CALL FILLHTRN30SEC(h_trn,nx_trn,ny_trn, &
          swlon_trn,nelon_trn,lon_sample,swlat_trn,nelat_trn,lat_sample,&
          dir_trndata(1:ldir_trndata),badvalue)

  ELSE IF( trndataopt == 4) THEN

    lfn = len_trim(usgs_dem_index_fn)

    CALL fillhtrn3sec(h_trn,nx_trn,ny_trn,trn3src,                      &
        swlon_trn,nelon_trn,lon_sample,swlat_trn,nelat_trn,lat_sample,  &
        dir_trndata(1:ldir_trndata),usgs_dem_index_fn(1:lfn),badvalue)

  ELSE

    WRITE(6,'(1x,a/1x,a,i3,a)') 'trndataopt must be from 1 to 4.',      &
        'The input value was ', trndataopt,', Program stopped.'
    STOP 991

  END IF

  CALL a3dmax(h_trn,1,nx_trn,1,nx_trn,1,ny_trn,1,ny_trn,                &
              1,1,1,1,                                                  &
              hmax,hmin,imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/,2(1x,a,g25.12,2(a,i5),/))')                               &
       'hmin =',hmin,' at i=',imin,', j=',jmin,                         &
       'hmax =',hmax,' at i=',imax,', j=',jmax

  WRITE(6,'(/a)') 'Finished reading terrain data.'
  !print*,' nx_trn, ny_trn =', nx_trn, ny_trn

  inith_set = 0
  badvfound = 0

  DO i=1,nx_trn
    DO j=1,ny_trn
      IF( h_trn(i,j) /= badvalue)THEN
        IF(inith_set == 0 ) THEN
          hmax = h_trn(i,j)
          hmin = h_trn(i,j)
          inith_set = 1
        ELSE
          hmax = MAX(hmax, h_trn(i,j))
          hmin = MIN(hmin, h_trn(i,j))
        END IF
      ELSE
        badvfound = 1
      END IF
    END DO
  END DO

  IF( badvfound == 1) THEN
    WRITE(6,'(/a,/a)')                                                  &
        'Missing value(s) found in sampled terrain array.',             &
        'The missing data need to be filled. Job stopped.'
    STOP
  END IF

  IF(inith_set == 0 ) THEN  ! If still not set
    hmax = 0.0
    hmin = 0.0
  END IF

  WRITE(6,'(/1x,a,2f10.3)')                                             &
      'Min and max in the sampled terrain data array = ',hmin, hmax
!
!-----------------------------------------------------------------------
!
!  Initialize ARPS terrain array hterain from h_trn using bi-linear
!  interpolation.
!
!-----------------------------------------------------------------------
!
  PRINT*, 'Southwest corn: latitude = ',swlat_trn, ', longitude = ',swlon_trn
  PRINT*, 'Sample latitude          = ',lat_sample,', sample longitude = ',lon_sample

  lond_sample = lon_sample/3600.0
  latd_sample = lat_sample/3600.0

  IF( trnanxopt == 1) THEN ! 1st order interpolation

    DO i=1,nx-1
      DO j=1,ny-1

        ih=MIN(nx_trn-1,                                                &
            MAX(1,INT((lonh(i,j)-(swlon_trn+swlon_offset))*3600.0/lon_sample)+1))

        jh=MIN(ny_trn-1,                                                &
            MAX(1,INT((lath(i,j)-(swlat_trn+swlat_offset))*3600.0/lat_sample)+1))

        !print*,'i,j,ih,jh=', i,j,ih,jh,h_trn(ih,jh)

        IF(h_trn(ih  ,jh  ) == badvalue .OR. h_trn(ih+1,jh)  == badvalue .OR. &
              h_trn(ih+1,jh+1) == badvalue .OR.                         &
              h_trn(ih  ,jh+1) == badvalue) THEN
          hterain(i,j)=badvalue
        ELSE
          lonh1=lonh(i,j)-(swlon_trn+swlon_offset)-(ih-1)*lond_sample
          lath1=lath(i,j)-(swlat_trn+swlat_offset)-(jh-1)*latd_sample
          w1=(latd_sample-lath1)*(lond_sample-lonh1)
          w2=(latd_sample-lath1)*lonh1
          w3=lath1*lonh1
          w4=lath1*(lond_sample-lonh1)
          hterain(i,j)=(w1*h_trn(ih  ,jh  )+w2*h_trn(ih+1,jh)           &
                       +w3*h_trn(ih+1,jh+1)+w4*h_trn(ih,jh+1))          &
                       /(latd_sample*lond_sample)
        END IF

      END DO
    END DO

  ELSE ! 2nd order interpolation

    DO i=1,nx-1
      DO j=1,ny-1

        ih=MIN(nx_trn-1,                                                &
            MAX(2,NINT((lonh(i,j)-(swlon_trn+swlon_offset))*3600.0/lon_sample)+1))

        jh=MIN(ny_trn-1,                                                &
            MAX(2,NINT((lath(i,j)-(swlat_trn+swlat_offset))*3600.0/lat_sample)+1))

        IF(h_trn(ih  ,jh  ) == badvalue .OR. h_trn(ih+1,jh)  == badvalue .OR. &
              h_trn(ih+1,jh+1) == badvalue .OR.                         &
              h_trn(ih  ,jh+1) == badvalue) THEN
          hterain(i,j)=badvalue
        ELSE
          IF(h_trn(ih-1,jh+1) == badvalue .OR.                          &
                h_trn(ih-1,jh)  == badvalue .OR.                        &
                h_trn(ih-1,jh-1) == badvalue .OR.                       &
                h_trn(ih+1,jh-1) == badvalue .OR.                       &
                h_trn(ih  ,jh-1) == badvalue) THEN ! Use 1st order

            hterain(i,j)=badvalue
            lonh1=lonh(i,j)-(swlon_trn+swlon_offset)-(ih-1)*lond_sample
            lath1=lath(i,j)-(swlat_trn+swlat_offset)-(jh-1)*latd_sample
            w1=(latd_sample-lath1)*(lond_sample-lonh1)
            w2=(latd_sample-lath1)*lonh1
            w3=lath1*lonh1
            w4=lath1*(lond_sample-lonh1)
            hterain(i,j)=(w1*h_trn(ih  ,jh  )+w2*h_trn(ih+1,jh)         &
                         +w3*h_trn(ih+1,jh+1)+w4*h_trn(ih,jh+1))        &
                         /(latd_sample*lond_sample)

          ELSE ! now real 2nd order interpolation

            d=lonh(i,j)-(swlon_trn+swlon_offset)-(ih-1)*lond_sample
            a=-lond_sample
            b=0.0
            c=+lond_sample

            w1= (d-b)*(d-c)/((a-b)*(a-c))
            w2= (d-a)*(d-c)/((b-a)*(b-c))
            w3= (d-a)*(d-b)/((c-a)*(c-b))
            h1=w1*h_trn(ih-1,jh-1)+w2*h_trn(ih,jh-1)                    &
                +w3*h_trn(ih+1,jh-1)
            h2=w1*h_trn(ih-1,jh  )+w2*h_trn(ih,jh  )                    &
                +w3*h_trn(ih+1,jh  )
            h3=w1*h_trn(ih-1,jh+1)+w2*h_trn(ih,jh+1)                    &
                +w3*h_trn(ih+1,jh+1)

            d=lath(i,j)-(swlat_trn+swlat_offset)-(jh-1)*latd_sample
            a=-latd_sample
            b=0.0
            c=+latd_sample
            w1= (d-b)*(d-c)/((a-b)*(a-c))
            w2= (d-a)*(d-c)/((b-a)*(b-c))
            w3= (d-a)*(d-b)/((c-a)*(c-b))
            hterain(i,j)=w1*h1+w2*h2+w3*h3
          END IF
        END IF
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Fill in Great Lakes.  The original data set corresponds to the
!  bottom of the lakes.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      tem(i,j)  = lath(i,j)
      tem2(i,j) = lonh(i,j)
    END DO
  END DO

  CALL fix_lake_eliv(hterain,tem,tem2,nx,ny)

  CALL a3dmax(hterain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,                  &
              hmax,hmin,imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/,2(1x,a,g25.12,2(a,i5),/))')                               &
       'htrnmin =',hmin,' at i=',imin,', j=',jmin,                      &
       'htrnmax =',hmax,' at i=',imax,', j=',jmax

!-----------------------------------------------------------------------
!
!  Apply smoothing on terrain height
!
!-----------------------------------------------------------------------

  DO i=1,nsmth
    CALL smth(hterain,nx,1,nx-1,ny,1,ny-1,2,tem)
!    call SMTH(hterain,nx,1,nx-1,ny,1,ny-1,1,tem)
  END DO

!
!-----------------------------------------------------------------------
!
!  Write out the terrain file on ARPS grid.
!
!-----------------------------------------------------------------------
!
  hmax=hterain(1,1)
  hmin=hterain(1,1)
  DO i=1,nx-1
    DO j=1,ny-1
      hterain(i,j) = MAX(0.0,hterain(i,j))
      hmax=MAX(hmax,hterain(i,j))
      hmin=MIN(hmin,hterain(i,j))
    END DO
  END DO

  WRITE(6,'(/1x,a,2f10.3)')                                             &
      'Max/min in the analyzed terrain array was ',hmax,hmin

  IF (terndmp /= 0) THEN
    CALL writtrn(nx,ny,dirname(1:lenstr)//runname(1:lfnkey)//'.trndata',&
                 dx,dy,mapproj,trulat1,trulat2,trulon,sclfct,           &
                 ctrlat,ctrlon,hterain)
  END IF

  DO j=1,ny
    DO i=1,nx
      xh(i,j)=xh(i,j)*0.001
      yh(i,j)=yh(i,j)*0.001
      !hterain(i,j)=hterain(i,j)*0.001
    END DO
  END DO

  IF (terndmp == 1) CALL trncntl(nx,ny,runname(1:lfnkey)//'.trndata',x,y)
!
!-----------------------------------------------------------------------
!
!  Generate graphical output.
!
!-----------------------------------------------------------------------
!
  outfilename = dirname(1:lenstr)//runname(1:lfnkey)
  lfn = lenstr + lfnkey
  CALL xdevic_new(1,outfilename,lfn,0)
  IF (lfn < 0) GOTO 970  ! skip plotting code

  CALL xctrbadv( 1 )  ! Turn on missing value checking for contouring
  CALL xvtrbadv( 1 )  ! Turn on missing value checking for vector plotting
  CALL xsetclrs_new(6,0)
  CALL xaxfmt('(I6)')
  CALL xbadval ( badvalue ) ! Set the missing value flag to -9999.0

  xl = (nx-3)*dx * 0.001
  yl = (ny-3)*dy * 0.001
  xorig = 0.0
  yorig = 0.0

  CALL xstpjgrd(mapproj,trulat1,trulat2,trulon,                         &
                ctrlat,ctrlon,xl,yl,xorig,yorig)

  CALL plot_trn(hterain,xh,yh,nx,1,nx-1,ny,1,ny-1,                      &
       -dx*0.001,(nx-2)*dx*0.001,dx*0.001,                              &
       -dy*0.001,(ny-2)*dy*0.001,dy*0.001, hctrinc, tem,xw,yw)

  DO i=1,nmapfile

    lmapfile=LEN(mapfile(i))
    CALL strlnth(mapfile(i), lmapfile)

    INQUIRE(FILE=mapfile(i)(1:lmapfile), EXIST = fexist )
    IF( .NOT.fexist) THEN
      WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//              &
          ' not found. Program stopped in ARPSTRN.'
      WRITE(6,'(1x,a,a)') 'Input was ',mapfile(i)(1:lmapfile)
      STOP 001
    END IF

    CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)

  END DO

! ctrlat = 40.28
! ctrlon = -111.56

! CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
! call xcolor(1)
! call xthick(3)
! call xpenup(ctrx*0.001,ctry*0.001-0.03*(ny-3)*dy*0.001)
! call xpendn(ctrx*0.001,ctry*0.001+0.03*(ny-3)*dy*0.001)
! call xpenup(ctrx*0.001-0.03*(nx-3)*dx*0.001,ctry*0.001)
! call xpendn(ctrx*0.001+0.03*(nx-3)*dx*0.001,ctry*0.001)

! ctrlat = 40.5
! ctrlon = -111.5

! CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
! call xcolor(1)
! call xthick(3)
! call xpenup(ctrx*0.001,ctry*0.001-0.03*(ny-3)*dy*0.001)
! call xpendn(ctrx*0.001,ctry*0.001+0.03*(ny-3)*dy*0.001)
! call xpenup(ctrx*0.001-0.03*(nx-3)*dx*0.001,ctry*0.001)
! call xpendn(ctrx*0.001+0.03*(nx-3)*dx*0.001,ctry*0.001)

!
  CALL xgrend

  970   CONTINUE
  WRITE(6,'(/,1x,a)') '==== ARPSTRN terminated normally ===='
  STOP

  980   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'Error occured when reading namelist input file. Program stopped.'

  STOP

  990   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'End of file reached when reading namelist input file. ',         &
      'Program stopped.'

  STOP
END PROGRAM arpstrn

SUBROUTINE fillhtrn3sec(h_trn,nx_trn,ny_trn,trn3src,                    &
           swlon_trn,nelon_trn,lon_sample,                              &
           swlat_trn,nelat_trn,lat_sample,dir_trndata,                  &
           usgs_dem_index_fn,badvalue)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!          11/25/95
!
!  MODIFICATION HISTORY:
!
!  2000/02/03 Gene Bassett
!  Added subroutine to fix the elevations in the Great Lakes.
!
!  09/09/2012 (Y. Wang)
!  Added SRTM data support and improved readability.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx_trn,ny_trn
  REAL,    INTENT(OUT) :: h_trn(nx_trn,ny_trn)

  INTEGER, INTENT(IN)  :: trn3src            !TINA

  INTEGER, INTENT(IN)  :: lat_sample,lon_sample ! Data sampling frequency in data index

  INTEGER, INTENT(IN)  :: swlat_trn,swlon_trn,nelat_trn,nelon_trn
                    ! SW and NE coner lat/lon of the assembled terrain patch.

  CHARACTER (LEN=*), INTENT(IN) :: dir_trndata,usgs_dem_index_fn
  REAL,              INTENT(IN) :: badvalue

!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES:
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: trn_name

  INTEGER :: nlat,nlon         ! Number of data points in one data file
  INTEGER, ALLOCATABLE :: h_block(:,:) ! Array to store 1x1 degree DEM patch

  INTEGER :: istatus
  INTEGER :: ii,jj,i,j,ilat,ilon,iilat,iilon
  REAL    :: tmax0,tmin0
  LOGICAL :: existfile
  CHARACTER (LEN=3) :: dtflag
  INTEGER :: lat_skip, lon_skip,ldir_trndata

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Make sure directory dir_trndata already exists.
!
!-----------------------------------------------------------------------

  ldir_trndata = len_trim(dir_trndata)

!  INQUIRE(FILE=dir_trndata(1:ldir_trndata),EXIST=existfile)
  CALL inquiredir(dir_trndata(1:ldir_trndata), existfile)
  IF( .NOT. existfile ) THEN
    WRITE(6,'(/1x,3a,3(/1x,a))')                                        &
        'Directory ',dir_trndata(1:ldir_trndata),' does not exist.',    &
        'Create the directory or create a link to the data directory.', &
        'This directory has to be writable, unless all required DEM data', &
        'exists locally.'
    CALL arpsstop('ERROR: data source directory not exists.',1)
  END IF

  IF (trn3src == 3 ) THEN
    nlat = 6001
    nlon = 6001
  ELSE
    nlat = 1201
    nlon = 1201
  END IF

  ALLOCATE (h_block(nlon,nlat), stat = istatus)
  CALL check_alloc_status(istatus, "fillhtrn3sec:h_block")

  DO ilat = swlat_trn, nelat_trn-1
    DO ilon = swlon_trn, nelon_trn-1

!-----------------------------------------------------------------------
!
!  Find the name of DEM file whose SW corner lat/lon = ilat/ilon.
!
!-----------------------------------------------------------------------

      CALL lookup_trn_file(ilat,ilon,usgs_dem_index_fn,trn3src,         &
                      trn_name,dtflag,istatus)
      IF( istatus /= 0) GO TO 310

!-----------------------------------------------------------------------
!
!  Acquire file trn_name from anonymous FTP server
!  edcftp.cr.usgs.gov if it does not already exist locally in dir_trndata
!
!-----------------------------------------------------------------------

      CALL get_trn_file(dir_trndata(1:ldir_trndata),trn3src,            &
                   trn_name,dtflag,istatus)
      IF(istatus /= 0) GO TO 310
!
!-----------------------------------------------------------------------
!
!  Read in 1x1 degree DEM data block into array h_block.
!
!-----------------------------------------------------------------------
!
      CALL read_trn_data(dir_trndata(1:ldir_trndata)//trn_name,trn3src, &
                         h_block,nlat,nlon,badvalue,istatus)
      IF(istatus /= 0) GO TO 310

      tmax0=h_block(1,1)
      tmin0=h_block(1,1)

      DO j=1,nlat
        DO i=1,nlon
          tmax0=MAX(tmax0,FLOAT(h_block(i,j)))
          tmin0=MIN(tmin0,FLOAT(h_block(i,j)))
        END DO
      END DO

      !WRITE(*,'(1x,3a,2F10.3)') 'max/min in data file ',                &
      !          dir_trndata(1:ldir_trndata),' = ',tmax0,tmin0
!
!-----------------------------------------------------------------------
!
!  Transfer selected data points from h_block into h_trn.
!
!-----------------------------------------------------------------------
!
      IF (trn3src == 3) THEN        ! 5 degree tiles
        iilat = MOD(ilat,5)*1200+1
        iilon = MOD(ilon,5)*1200+1
      ELSE                          ! 1 degree tiles
        iilat = 1
        iilon = 1
      END IF

      lat_skip = max(1,lat_sample/3)
      lon_skip = max(1,lon_sample/3)

      !ii = (nlon-1)/lon_skip*(ilon-swlon_trn)      ! process 1 degree at 3 sec resolution
      !jj = (nlat-1)/lat_skip*(ilat-swlat_trn)
      ii = 1200/lon_skip*(ilon-swlon_trn)
      jj = 1200/lat_skip*(ilat-swlat_trn)

      !DO j=1,(nlat-1)/lat_skip+1
      !  DO i=1,(nlon-1)/lon_skip+1
      DO j = 1, 1200/lat_skip+1
        DO i = 1, 1200/lon_skip+1
          h_trn(ii+i,jj+j)=                                             &
                FLOAT(h_block((i-1)*lon_skip+iilon,(j-1)*lat_skip+iilat))
        END DO
      END DO

      CYCLE

      310     CONTINUE

      WRITE(6,'(/1x,a,2i5,a,/1x,a)')                                      &
          'DEM file with S-W lat/lon=',ilat,ilon,' was not found or not', &
          'read correctly. The DEM array is filled with missing values.'

      ii = (nlon-1)/lon_skip*(ilon-swlon_trn)
      jj = (nlat-1)/lat_skip*(ilat-swlat_trn)

      DO j=1,(nlat-1)/lat_skip+1
        DO i=1,(nlon-1)/lon_skip+1
          h_trn(ii+i,jj+j)= badvalue
        END DO
      END DO

    END DO
  END DO

  DEALLOCATE(h_block)

  RETURN
END SUBROUTINE fillhtrn3sec

SUBROUTINE fillhtrn5min(h_trn,nx_trn,ny_trn,                            &
           swlon_trn,nelon_trn,lon_sample,                              &
           swlat_trn,nelat_trn,lat_sample,dir_trndata,fn_trndata,       &
           badvalue)

!-----------------------------------------------------------------------
!
!  read the TerrainBase data
!  This original program modified by LDI
!  This promgram used to make the 5min topo input data using NCAR 5min topo
!  output is the 5 min topo for ARPS terrain program
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: itopo(4320)   ! Array for topo along one latitude circle
  INTEGER :: imin,imax,jmin,jmax,i,j

  INTEGER :: nx_trn,ny_trn
  REAL :: h_trn(nx_trn,ny_trn)

  INTEGER :: swlat_trn,swlon_trn,nelat_trn,nelon_trn
              ! SW and NE coner lat/lon of the assembled terrain patch.

  INTEGER :: lon_sample, lat_sample
  CHARACTER (LEN=256) :: fn_trndata, dir_trndata
  INTEGER :: ldir_trndata,lfn_trndata
  REAL    :: badvalue

  REAL    :: hmax, hmin
  INTEGER :: inith_set, iskip, jskip, ii
  LOGICAL :: existfile

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  imin = swlon_trn*12+1
  imax = nelon_trn*12+1

  jmax = (90-swlat_trn)*12+1
  jmin = (90-nelat_trn)*12+1

  PRINT *, 'imin = ',imin,' imax = ',imax,' jmini = ',jmin,' jmax = ',jmax

  ldir_trndata = 120
  CALL strlnth( dir_trndata, ldir_trndata)
  lfn_trndata = 80
  CALL strlnth( fn_trndata, lfn_trndata)

  INQUIRE(FILE=dir_trndata(1:ldir_trndata)//fn_trndata(1:lfn_trndata),  &
          EXIST=existfile)

  IF ( .not. existfile) THEN

    WRITE(6,'(a,a,/a)')                                                 &
    'File ',dir_trndata(1:ldir_trndata)//fn_trndata(1:lfn_trndata),     &
    ' not found. Program stopped in subroutine FILLHTRN5MIN.'
    STOP

  END IF

  OPEN(10,FILE=dir_trndata(1:ldir_trndata)//                            &
       fn_trndata(1:lfn_trndata),STATUS='old',FORM='formatted')

  iskip = lon_sample/300
  jskip = lat_sample/300

  DO j=1,jmax  ! jmax =< 2160

    READ(10,100)(itopo(i),i=1,4320)

    IF(j >= jmin.AND.                                                   &
          MOD((jmax-jmin+1)-(j-jmin+1),jskip) == 0)THEN  ! Fill in array h_trn

      DO i=imin, imax, iskip
        ii = i
        IF( i < 0) ii=4320+i
        h_trn((i-imin)/iskip+1,                                         &
             ((jmax-jmin+1)-(j-jmin+1))/jskip+1)=itopo(ii)
      END DO

    END IF

  END DO

  100   FORMAT(20I6)

!  CALL XDEVIC

  IF( .false.) THEN  ! Plotting is true

    nx_trn = (12*(nelon_trn-swlon_trn))/iskip+1
    ny_trn = (12*(nelat_trn-swlat_trn))/jskip+1

!  print*,'nx_trn,ny_trn,imax-imin+1,jmax-jmin+1=',
!    :        nx_trn,ny_trn,imax-imin+1,jmax-jmin+1

    inith_set = 0

    DO i=1,nx_trn
      DO j=1,ny_trn
        IF( h_trn(i,j) /= badvalue)THEN
          IF(inith_set == 0 ) THEN
            hmax = h_trn(i,j)
            hmin = h_trn(i,j)
            inith_set = 1
          ELSE
            hmax = MAX(hmax, h_trn(i,j))
            hmin = MIN(hmin, h_trn(i,j))
          END IF
        END IF
      END DO
    END DO

    IF(inith_set == 0 ) THEN  ! If still not set
      hmax = 0.0
      hmin = 0.0
    END IF

    PRINT*,' hmin, hmax =', hmin, hmax

!  hmax = h_trn(1,1)
!  hmin = h_trn(1,1)
!  DO 1203 i=1,nx_trn
!  DO 1203 j=1,ny_trn
!     hmax = max(hmax, h_trn(i,j))
!     hmin = min(hmin, h_trn(i,j))
!    1203  CONTINUE
!  print*,' hmin, hmax =', hmin, hmax

  END IF

  CLOSE(2)

  RETURN
END SUBROUTINE fillhtrn5min

SUBROUTINE lookup_trn_file(swlat,swlon,usgs_dem_index_fn,trn3src,       &
                      trn_name,dtflag,ierr)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find a 1x1 degree DEM file whose SW corner is at (swlat,swlon).
!  dtflag= 'US' or 'AK', indicating if the data file is in Alaska.
!  Note: The longitude W is given as positive numbers in the Index.
!  This should be fine as long as this program is used for
!  accessing US 3" data only. The code should be modified slightly
!  to handle global data.
!
!  INPUT : swlat, swlon,usgs_dem_index_fn
!  OUTPUT: trn_name, dtflag, ierr
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,           INTENT(IN) :: swlat,swlon
  INTEGER,           INTENT(IN) :: trn3src               !TINA

  CHARACTER (LEN=*), INTENT(IN)  :: usgs_dem_index_fn
  CHARACTER (LEN=*), INTENT(OUT) :: trn_name
  CHARACTER (LEN=3), INTENT(OUT) :: dtflag

  INTEGER,           INTENT(OUT)  :: ierr

!-----------------------------------------------------------------------

  CHARACTER (LEN=80) :: name1
  CHARACTER (LEN=1 ) :: ch1, ch2

  INTEGER :: selat,selon,lat,lon,ll
  REAL    :: alat,alon

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( trn3src == 3) THEN  ! SRTM data set
    lat = (swlat/5)*5      ! Find latitude with 5 degree tiles
    lon = (swlon/5)*5

    selat = (60-lat)/5     ! Find tile numbers
    selon = (lon+180+5)/5

    WRITE(trn_name,'(a,I2.2,a,I2.2,a)') 'srtm_',selon,'_',selat,'.asc'

    WRITE(6,'(/,1x,2a,2(a,I4),a)') 'Need 5 degree tile "',TRIM(trn_name), &
                    '" with the matching lat/lon = ',swlat,' /',swlon,'.'

  ELSE IF ( trn3src == 2) THEN  ! SRTM Eurasia data set

    selat = ABS(swlat)
    selon = ABS(swlon)
    IF (swlon > 0) THEN
      ch1 = 'E'
    ELSE
      ch1 = 'W'
    END IF

    IF (swlat > 0) THEN
      ch2 = 'N'
    ELSE
      ch2 = 'S'
    END IF

    WRITE(trn_name,'(a,I2.2,a,I3.3,a)') ch2,selat,ch1,selon,'.hgt'

    WRITE(6,'(/,1x,2a,2(a,I4),a)') 'Need tile "',TRIM(trn_name),        &
                    '" with the matching lat/lon = ',swlat,'/',swlon,'.'

  ELSE

    OPEN (50,FILE=trim(usgs_dem_index_fn),STATUS='old',IOSTAT=ierr)

    IF( ierr /= 0) THEN
      WRITE(6,'(/1x,a,a,a,/1x,a)') 'Opening file ',                     &
           trim(usgs_dem_index_fn),' failed. Program stopped.'
      STOP
    END IF

    READ (50,*)


    IF (trn3src == 1) THEN   ! Europe date sets

      selat = swlat
      selon = swlon

      WRITE(6,'(/1x,2(a,i4))')                                          &
              'Searching for a file with the SW corner latitude=',selat,&
              ' longitude=', selon

    ELSE
      selat = swlat
      selon = -(swlon + 1)

      WRITE(6,'(/1x,2(a,i4))')                                          &
              'Searching for a file with the SE corner latitude=',selat,&
              ' longitude=', selon
    END IF

    100   READ (50,'(a3,a30,2f10.1)',END=300) dtflag,name1, alat, alon

    lat=nint(alat)
    lon=nint(alon)

    IF (lat == selat.AND.lon == selon) THEN
      trn_name=name1
      ll=INDEX(trn_name,'.gz')
      WRITE(6,'(1x,a,a,a)') 'Found the entry ',trn_name(1:ll+2),        &
                            ' with the matching lat/lon.'

      CLOSE (50)
      RETURN
    ELSE
      GO TO 100
    END IF

    300   CONTINUE
    CLOSE (50)

    PRINT*,'No entry in data base for this latitude and longitude.'
    ierr = 1

  END IF

  RETURN
END SUBROUTINE lookup_trn_file

SUBROUTINE get_trn_file(dir_trndata,trn3src,trn_name,dtflag,istatus)

  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN)  :: dir_trndata
  CHARACTER (LEN=*), INTENT(IN)  :: trn_name
  CHARACTER (LEN=3), INTENT(IN)  :: dtflag
  INTEGER,           INTENT(IN)  :: trn3src
  INTEGER,           INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER ::ldir_trndata

  CHARACTER (LEN=256) :: ch1,ch2,ch3
!  CHARACTER (LEN=80)  :: dr
  LOGICAL :: existfile,existfile1
  INTEGER :: l1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 1

  ldir_trndata=LEN_TRIM(dir_trndata)

  IF (trn3src == 3) THEN

    l1 = INDEX(trn_name,'.asc')

    INQUIRE(FILE=dir_trndata(1:ldir_trndata)//trn_name(1:l1+3),EXIST=existfile)

    IF (existfile) THEN

      WRITE(*,'(1x,5a)') 'Found "',trn_name(1:l1+3),'" in ',TRIM(dir_trndata),' locally.'

    ELSE

      WRITE(ch1,'(3a)') dir_trndata(1:ldir_trndata),trn_name(1:l1-1),'.zip'

      INQUIRE(FILE=TRIM(ch1),EXIST=existfile)

      IF (existfile) THEN

        WRITE(*,'(1x,3a)') 'Found compressed file "',TRIM(ch1),'" locally. Unzipping'

        ch2 = 'unzip '//TRIM(ch1)
        CALL unixcmd(ch2)

      ELSE

        WRITE(*,'(1x,3a)') 'File "',TRIM(trn_name),'" not found. Fetching from CGIAR-CSI server ...'

        ch3 = trn_name(1:l1-1)//'.zip'
        WRITE(ch2,'(4a)') 'wget ','http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_v41/SRTM_Data_ArcASCII/', &
                          TRIM(ch3)
        CALL unixcmd(ch2)

        WRITE(*,'(1x,3a)') 'Uncompressing file "',TRIM(ch3),'" ...'
        ch2 = 'unzip '//trim(ch3)
        CALL unixcmd(ch2)

        CALL unixcmd('/bin/rm *.prj readme.txt '//TRIM(ch3))

      END IF

    END IF

  ELSE IF (trn3src == 2) THEN

    INQUIRE(FILE=dir_trndata(1:ldir_trndata)//TRIM(trn_name),EXIST=existfile)

    IF (existfile) THEN

      WRITE(*,'(1x,5a)') 'Found "',TRIM(trn_name),'" in ',TRIM(dir_trndata),' locally.'

    ELSE

      WRITE(ch3,'(2a)') TRIM(trn_name),'.zip'
      WRITE(ch1,'(2a)') dir_trndata(1:ldir_trndata),TRIM(ch3)

      INQUIRE(FILE=TRIM(ch1),EXIST=existfile)

      IF (existfile) THEN

        WRITE(*,'(1x,3a)') 'Found compressed file "',TRIM(ch1),'" locally. Unzipping'

        ch2 = 'unzip '//TRIM(ch1)
        CALL unixcmd(ch2)

      ELSE

        WRITE(*,'(/,1x,3a)') 'File "',TRIM(trn_name),'" not found. Fetching from USGS.gov server ...'

        WRITE(ch2,'(2a)') 'wget http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/', &
                          TRIM(ch3)
        CALL unixcmd(ch2)

        WRITE(*,'(1x,3a)') 'Uncompressing file "',TRIM(ch3),'" ...'
        ch2 = 'unzip '//trim(ch3)
        CALL unixcmd(ch2)

        CALL unixcmd('/bin/rm '//TRIM(ch3))

      END IF

    END IF

  ELSE
    l1=INDEX(trn_name,'.gz')

    !  IF( dtflag == 'AK') THEN
    !    dr='Alaska/'//CHAR(ICHAR(trn_name(1:1))-32)
    !  ELSE
    !    dr=trn_name(1:1)
    !    dr=CHAR(ICHAR(dr(1:1))-32)
    !  END IF
    !
    !ch1='cd '//dr
    ch1='cd '//TRIM(dtflag)                  ! use CAPS local server
    ch2='dir '//trn_name(1:l1+2)
    ch3='get '//trn_name(1:l1+2)

    !
    !  check to see if dem file pre-exists
    !
    INQUIRE(FILE=dir_trndata(1:ldir_trndata)//trn_name(1:l1+2),         &
            EXIST=existfile)

    INQUIRE(FILE=dir_trndata(1:ldir_trndata)//trn_name(1:l1-1),         &
            EXIST=existfile1)

    IF (existfile) THEN

      PRINT *,trn_name(1:l1+2),' exists locally'

    ELSE IF (existfile1 ) THEN

      PRINT *,trn_name(1:l1-1),' exists locally'

    ELSE
      !
      !  write the ftp input file
      !
      OPEN (29,FILE='ftp.in')
      WRITE (29,'(a)') 'binary'
     !WRITE (29,'(a)') 'cd pub/data/DEM/250'
      WRITE (29,'(a)') 'cd ARPS/ARPS.data/demtopo3.data'     ! use CAPS ftp server
      WRITE (29,'(a)') TRIM(ch1)
      WRITE (29,'(a)') TRIM(ch2)
      WRITE (29,'(a)') TRIM(ch3)
      WRITE (29,'(a)') 'quit'
      CLOSE (29)
      !
      !  invoke ftp command
      !
      !ch1='ftp edcftp.cr.usgs.gov< ftp.in'
      ch1='ftp caps.ou.edu< ftp.in'

      WRITE(6,'(1x,a)') '=== Connecting to ftp site to get the file ...'
      WRITE(6,'(5x,2a)') '$>',TRIM(ch1)

      CALL unixcmd(ch1)

      CALL unixcmd('mv '//trn_name(1:l1+2) //' '//dir_trndata(1:ldir_trndata))

      WRITE(6,'(1x,3a)') 'File ',trn_name(1:l1+2),                      &
              ' transferred to the local machine now.'

    END IF

  END IF

  istatus=0

  RETURN
END SUBROUTINE get_trn_file

SUBROUTINE read_trn_data(trn_fullname,trn3src,                          &
                         h_block,nlat,nlon,badvalue,istatus)

  USE arps_precision

  IMPLICIT NONE

  CHARACTER (LEN=256), INTENT(IN)  :: trn_fullname
  INTEGER,             INTENT(IN)  :: nlat,nlon
  INTEGER,             INTENT(OUT) :: h_block(nlon,nlat)

  INTEGER,             INTENT(IN)  :: trn3src
  REAL,                INTENT(IN)  :: badvalue
  INTEGER,             INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: ch1, ch2, ch3, ch4
  CHARACTER (LEN=144) :: head

  INTEGER :: i,j,ii,jj,ijump,m,n,l1,ltd,lnd,idum
  REAL    :: xjumpinv,tmax0,tmin0,tmax,tmin,dum,amx,amn
  REAL    :: dum15(15), dum6(6)
  REAL    :: selnsec,seltsec
  LOGICAL :: existfile

  INTEGER(KIND=INT16), ALLOCATABLE :: inblock(:,:)  ! for SRTM Eurasia data interpolation
  INTEGER :: nxin,nyin, nodata     ! for reading 90m SRTM global data
  REAL    :: xll, yll, cellsize
  INTEGER :: badtotal, badnum      ! for SRTM Eurasia data interpolation
  INTEGER :: outdist, maxdist

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus=1

  IF (trn3src == 3) THEN

    WRITE(*,'(1x,3a)') 'Reading ',TRIM(trn_fullname),' ...'

    OPEN (29,FILE=TRIM(trn_fullname),FORM='FORMATTED',STATUS='OLD',     &
             ACTION='READ')

    READ(29,*) ch1,nxin
    READ(29,*) ch2,nyin
    READ(29,*) ch3,xll
    READ(29,*) ch4,yll
    READ(29,*) ch4,cellsize
    READ(29,*) ch4,nodata

    IF (nxin /= nlon .OR. nyin /= nlat) THEN
      WRITE(*,'(1x,2(a,2I4),/)') 'WRONG size, expected ',nlon,nlat,'; in file ',nxin,nyin
      STOP
    END IF

    DO j = nlat,1,-1
      READ(29,*) (h_block(i,j),i=1,nlon)
    END DO

    CLOSE (29)

    IF( .true. ) THEN
      tmax0=h_block(1,1)
      tmin0=h_block(1,1)
      DO i=1,nlon
        DO j=1,nlat
          tmax0=MAX(tmax0,FLOAT(h_block(i,j)))
          tmin0=MIN(tmin0,FLOAT(h_block(i,j)))
        END DO
      END DO

      PRINT*,'READTRN: max/min in the data read in = ',tmax0,tmin0
    END IF

    WRITE(*,'(1x,a)') 'Terrain data stored into the terrain array.'

  ELSE IF (trn3src == 2) THEN

    ALLOCATE(inblock(nlon,nlat), STAT = istatus)
    CALL check_alloc_status(istatus, "read_trn_data:inblock")

    WRITE(*,'(1x,3a)') 'Reading ',TRIM(trn_fullname),' ...'

    OPEN (29,FILE=TRIM(trn_fullname),FORM='UNFORMATTED',STATUS='OLD',   &
             ACCESS='STREAM',ACTION='READ')

    READ(29) inblock

    CLOSE (29)

    IF( .true. ) THEN
      tmax0=inblock(1,1)
      tmin0=inblock(1,1)
      DO j=1,nlat
        DO i=1,nlon
          tmax0=MAX(tmax0,FLOAT(inblock(i,j)))
          tmin0=MIN(tmin0,FLOAT(inblock(i,j)))
        END DO
      END DO

      PRINT*,'READTRN: max/min in the data read in = ',tmax0,tmin0
    END IF

    DO j=1,nlat
      DO i=1,nlon
        ii = i
        jj = nlat-j+1       ! reverse south and north

        IF ( inblock(i,j) == -32768) THEN  ! Average badvalue
          CALL fill_missing(inblock,nlat,nlon,i,j,h_block(ii,jj),outdist,istatus)
          maxdist = MAX(maxdist,outdist)
        ELSE
          h_block(ii,jj)=inblock(i,j)
        END IF
      END DO
    END DO

    WRITE(*,'(1x,a,/,1x,a,I0,a)') 'Terrain data stored into the terrain array.', &
    'To fill the missing grid, the program has searched maximum ',maxdist,' grid points far away.'

    IF( .TRUE.) THEN
      tmax0=h_block(1,1)
      tmin0=h_block(1,1)
      DO j=1,nlat
        DO i=1,nlon
          tmax0=MAX(tmax0,FLOAT(h_block(i,j)))
          tmin0=MIN(tmin0,FLOAT(h_block(i,j)))
        END DO
      END DO

      PRINT*,'max/min after storage                = ',tmax0,tmin0
    END IF

    DEALLOCATE(inblock)

  ELSE

    l1=INDEX(trn_fullname,'.gz')

    INQUIRE(FILE='tmp1.file', EXIST = existfile)
    IF(existfile)THEN
      ch2='/bin/rm tmp1.file'
      CALL unixcmd(ch2)
    ENDIF

    INQUIRE(FILE=trn_fullname(1:l1-1),EXIST=existfile)

    IF(existfile) THEN  ! unzipped version exists
      ch1='cp '//trn_fullname(1:l1-1)//' tmp1.file'
      CALL unixcmd(ch1)
    ELSE                ! need to unzip
      ch1='gunzip -c '//trn_fullname(1:l1+2)//'> tmp1.file'
      PRINT *,'Unzipping the terrain file..'
      CALL unixcmd(ch1)
      PRINT *,'Terrain file unzipping competed..'
    END IF

    INQUIRE(FILE='tmp2.file', EXIST = existfile)
    IF(existfile)THEN
      ch2='/bin/rm tmp2.file'
      CALL unixcmd(ch2)
    ENDIF

    IF (trn3src == 0) THEN              !TINA added

      ch4='dd if=tmp1.file of=tmp2.file ibs=4096 cbs=1024 conv=unblock'
      PRINT *,'Converting the terrain file..'
      CALL unixcmd(ch4)
      PRINT *,'Conversion completed..'

    ELSE                               !TINA

      ch4='cp tmp1.file tmp2.file'
      !ch4='ls -l; cp tmp1.file tmp2.file'
      CALL unixcmd(ch4)
      PRINT *,'Copying file...'

    ENDIF                              !TINA

    INQUIRE(FILE='tmp1.file', EXIST = existfile)
    IF(existfile)THEN
      ch2='/bin/rm tmp1.file'
      CALL unixcmd(ch2)
    ENDIF
   !
   !  Read in terrain data
   !
    IF (trn3src == 0) THEN              !TINA added

      OPEN (29,FILE='tmp2.file')

      READ (29,11) head,idum,idum,idum,idum,(dum15(i),i=1,15),            &
          idum,idum,idum,(dum6(i),i=1,6),selnsec,seltsec,tmin,tmax,dum,   &
          idum,dum,dum,dum,idum,m
      11 FORMAT(a144,4I6,15E24.0,3I6,11E24.0,i6,3E12.0,2I6)

      PRINT*,'Min and max in the input DEM file header=',tmin,tmax

      ltd=seltsec/3600.
      lnd=selnsec/3600.

      WRITE (6,'(/1x,a,2i5)')                                             &
        'South-east corner latitude and longitude from the file ',ltd,lnd

      WRITE (6,'(/1x,a)')                                                 &
        'Reading terrain data and store into the terrain array.'

      IF( m == 1201) THEN
        ijump=1
      ELSE IF( m == 601) THEN
        ijump=2
      ELSE IF( m == 401) THEN
        ijump=3
      ELSE
        WRITE(6,'(/1x,a,/1x,i6,a)')                                       &
          'The DEM file must have 1201, 601 or 401 data in the lon. dir.' &
          ,m,' data points were found. Program stopped.'
        STOP
      END IF

      DO i=1,1201,ijump
        READ (29,'(4I6,5E24.0,146I6,/(170I6))')                         &
                 idum,idum,n,idum,dum,dum,dum,amx,amn,(h_block(i,j),j=1,n)
      END DO

      IF( n /= 1201) THEN
        WRITE(6,'(/1x,a,/1x,i6,a)')                                       &
          'The DEM file must have 1201 data in the lat. dir.'             &
          ,n,' data points were found. Program stopped.'
        STOP
      END IF

    ELSE                                !TINA - options for Europe files

      OPEN (29,FILE='tmp2.file')
      DO j=1,1201
        DO i=1,1201
          READ (29,*) h_block(i,j)
        END DO
      END DO

      ijump=1
      n=1201

    ENDIF                               !TINA block for Europe data

    IF( .true. ) THEN
      tmax0=h_block(1,1)
      tmin0=h_block(1,1)
      DO i=1,1201,ijump
        DO j=1,n
          tmax0=MAX(tmax0,FLOAT(h_block(i,j)))
          tmin0=MIN(tmin0,FLOAT(h_block(i,j)))
        END DO
      END DO

      PRINT*,'READTRN: max/min in the data read in=',tmax0,tmin0
    END IF
    !
    !  Linearly interpolate terrain data to all 3 second grid points
    !  when the input data resolution is coarser than 3 second.
    !
    xjumpinv=1.0/ijump

    DO ii=1,ijump-1
      DO i=1,1201-ijump,ijump
        DO j=1,n
          h_block(i+ii,j)=h_block(i,j)+(h_block(i+ijump,j)              &
                          -h_block(i,j))*ii*xjumpinv
        END DO
      END DO
    END DO


    IF( .false.) THEN
      tmax0=h_block(1,1)
      tmin0=h_block(1,1)
      DO i=1,1201
        DO j=1,n
          tmax0=MAX(tmax0,FLOAT(h_block(i,j)))
          tmin0=MIN(tmin0,FLOAT(h_block(i,j)))
        END DO
      END DO

      PRINT*,'max/min in the data read in=',tmax0,tmin0
    END IF

    PRINT*,'Terrain data stored into the terrain array.'

    CLOSE (29)
    CALL unixcmd('/bin/rm tmp2.file')

  END IF

  istatus = 0

  RETURN
END SUBROUTINE read_trn_data

SUBROUTINE fillhtrn30sec(h_trn,nx_trn,ny_trn,                           &
           swlon_trn,nelon_trn,lon_sample,                              &
           swlat_trn,nelat_trn,lat_sample,dir_trndata,badvalue)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!          1/31/2002 (based on fillhtrn3sec)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlat,nlon         ! Number of data points in the 3"
                               ! one degree DEM data block.
  INTEGER, ALLOCATABLE :: h_block(:,:) ! Array to store 1x1 degree DEM patch

  INTEGER :: nx_trn,ny_trn
  REAL :: h_trn(nx_trn,ny_trn)

  CHARACTER (LEN=256) :: trn_name

  INTEGER :: lat_sample,lon_sample ! Data sampling frequency in data index

  INTEGER :: swlat_trn,swlon_trn,nelat_trn,nelon_trn
                                       ! SW and NE coner lat/lon of the
! the assembled terrain patch.
  INTEGER :: istatus,istatus1,istatus2,ii,jj,i,j,ilat,ilon
  REAL :: badvalue, tmax0,tmin0
  CHARACTER (LEN=3) :: dtflag
  LOGICAL :: existfile
  CHARACTER (LEN=*) :: dir_trndata
  INTEGER :: lat_skip, lon_skip,ldir_trndata
!

!-----------------------------------------------------------------------
!
!  Make sure directory dir_trndata already exists.
!
!-----------------------------------------------------------------------

  print*,'Inside  fillhtrn30sec '
  ldir_trndata = len_trim(dir_trndata)

! INQUIRE(FILE=dir_trndata(1:ldir_trndata),EXIST=existfile)
! IF( .NOT. existfile ) THEN
!   WRITE(6,'(/1x,3a,3(/1x,a))')                                        &
!       'Directory ',dir_trndata(1:ldir_trndata),' does not exist.',    &
!       'Create the directory or create a link to the data directory.', &
!       'This directory has to be writable, unless all required DEM data', &
!       'exists locally.'
! END IF

  nlat = 1200
  nlon = 1200
  ALLOCATE(h_block(nlat,nlon),STAT = istatus)
  CALL check_alloc_status(istatus, "fillhtrn30sec:h_block")

  WRITE(*,'(1x,a,4I4)') 'swlat_trn,nelat_trn,swlon_trn,nelon_trn = ',swlat_trn,nelat_trn,swlon_trn,nelon_trn

  DO ilat=swlat_trn,nelat_trn-10,10
    DO ilon=swlon_trn,nelon_trn-10,10

     print*,' '
     print*,'SW corner lat, lon of 10x10 degree patch to be read=',ilat, ilon

!-----------------------------------------------------------------------
!  Find the name of the 10x10 degree file whose SW corner lat/lon = ilat/ilon.
!-----------------------------------------------------------------------

    IF(ilon.lt.0.0.and.ilat.ge.0.0)write(trn_name,'(a,i3.3,a,i3.3,a)')  &
       'W',ABS(ilon),'N',ABS(ilat),'.30s_ascdat'
    IF(ilon.lt.0.0.and.ilat.lt.0.0)write(trn_name,'(a,i3.3,a,i3.3,a)')  &
       'W',ABS(ilon),'S',ABS(ilat),'.30s_ascdat'
    IF(ilon.ge.0.0.and.ilat.ge.0.0)write(trn_name,'(a,i3.3,a,i3.3,a)')  &
       'E',ABS(ilon),'N',ABS(ilat),'.30s_ascdat'
    IF(ilon.ge.0.0.and.ilat.lt.0.0)write(trn_name,'(a,i3.3,a,i3.3,a)')  &
       'E',ABS(ilon),'S',ABS(ilat),'.30s_ascdat'

!-----------------------------------------------------------------------
!
!  Acquire file trn_name from CAPS anonymous FTP server
!  ftp://ftp.caps.ou.edu/pub/ARPS/ARPS.data/arpstopo30.data if it does
!  not already exist locally in dir_trndata
!
!-----------------------------------------------------------------------

      CALL get_30s_trn(trim(dir_trndata),trim(trn_name),istatus1)
      IF(istatus1 /= 0) GO TO 310
!
!-----------------------------------------------------------------------
!
!  Read in 10x10 degree DEM data block into array h_block.
!
!-----------------------------------------------------------------------
!
      print*,'To read data from ', dir_trndata(1:ldir_trndata)//trim(trn_name)

      CALL read_30s_trn(dir_trndata(1:ldir_trndata)//trim(trn_name),    &
                        h_block,nlat,nlon,istatus2)

      print*,'Done reading data in ', dir_trndata(1:ldir_trndata)//trim(trn_name)

      DO i=1,nlon
        DO j=1,nlat
          IF(h_block(i,j).eq.-9999)h_block(i,j)=0.0 ! in ocean by definition
                                                    ! of the data set
        ENDDO
      ENDDO

      IF(istatus2 /= 0) GO TO 310
!
!-----------------------------------------------------------------------
!
!  Transfer selected data points from h_block into h_trn.
!
!-----------------------------------------------------------------------
!
      lat_skip = max(1,lat_sample/30)
      lon_skip = max(1,lon_sample/30)

      ii = nlon/lon_skip*((ilon-swlon_trn)/10)
      jj = nlat/lat_skip*((ilat-swlat_trn)/10)

      DO j=1,nlat/lat_skip
        DO i=1,nlon/lon_skip
          h_trn(ii+i,jj+j)=                                             &
                FLOAT(h_block((i-1)*lon_skip+1,(j-1)*lat_skip+1))
        END DO
      END DO

      print*,'data values transferred from h_block to h_trn'

      CYCLE

      310   CONTINUE

      WRITE(6,'(/1x,a,2i5,a,/1x,a)')                                    &
          '30second data file with S-W lat/lon=',ilat,ilon,' was not found or not', &
          'read correctly. The data array is filled with missing values.'

      ii = nlon/lon_skip*((ilon-swlon_trn)/10)
      jj = nlat/lat_skip*((ilat-swlat_trn)/10)

      DO j=1,nlat/lat_skip
        DO i=1,nlon/lon_skip
          h_trn(ii+i,jj+j)= badvalue
        END DO
      END DO

    END DO
  END DO

  DEALLOCATE( h_block )

  print*,'end of  fillhtrn30sec '

  RETURN
END SUBROUTINE fillhtrn30sec

SUBROUTINE get_30s_trn(dir_trndata,trn_name,istatus)

  IMPLICIT NONE

  CHARACTER (LEN=*) :: dir_trndata
  CHARACTER (LEN=*) :: trn_name
  INTEGER :: istatus,ldir_trndata

  CHARACTER (LEN=256) :: ch1
  LOGICAL :: existfile,existfile1
  INTEGER :: nunit

  ch1='get '//trim(trn_name)//'.gz'
!
!  check to see if dem file pre-exists
!
  istatus = 1

  INQUIRE(FILE=trim(dir_trndata)//trim(trn_name),EXIST=existfile)
  INQUIRE(FILE=trim(dir_trndata)//trim(trn_name)//'.gz',EXIST=existfile1)

  IF (existfile) THEN
    PRINT *,trim(trn_name),' exists locally'
  ELSE IF (existfile1 ) THEN
    PRINT *,trim(trn_name)//'.gz',' exists locally'
  ELSE

!    PRINT *,trim(trn_name),' or its compressed version not found in '
!    PRINT *,'in directory ',trim(dir_trndata),  &
!            '. Will try to ftp it from ftp.caps.ou.edu.'
!    Print *,'It requires that you have .netrc setup in your home directory'
!    Print *,'for automatic anonymous ftp access to ftp.caps.ou.edu. '

    PRINT *,trim(trn_name),' or its compressed version not found in '
    PRINT *,'in directory ',trim(dir_trndata),  &
            '. Will try to ftp it from caps.ou.edu.'
    Print *,'It requires that you have .netrc setup in your home directory'
    Print *,'for automatic anonymous ftp access to caps.ou.edu. '

!
!   write the ftp input file
!
    call getunit( nunit)
    OPEN (nunit,FILE='ftp.in')
    WRITE (nunit,1) 'binary'
!    WRITE (nunit,1) 'cd pub/ARPS/ARPS.data/arpstopo30.data'
    WRITE (nunit,1) 'cd ARPS/ARPS.data/arpstopo30.data'
    WRITE (nunit,1) TRIM(ch1)
    WRITE (nunit,1) 'quit'
    CLOSE (nunit)
    call retunit(nunit)
1   FORMAT(a)
!
!  invoke ftp command
!
!    ch1='ftp ftp.caps.ou.edu < ftp.in'
    ch1='ftp caps.ou.edu < ftp.in'
    PRINT *,'Connecting to ftp site to get the file...'

    CALL unixcmd(ch1)
    CALL unixcmd('mv '//trim(trn_name)//'.gz '//trim(dir_trndata))

    PRINT *,'File ',trim(trn_name)//'.gz ',                             &
           ' transferred to the local machine now.'

  END IF

  istatus=0
!
  RETURN
END SUBROUTINE get_30s_trn


SUBROUTINE read_30s_trn(trn_fullname,h_block,nlat,nlon,istatus)

  IMPLICIT NONE

  CHARACTER (LEN=*) :: trn_fullname
  INTEGER :: nlat,nlon
  INTEGER :: h_block(nlat,nlon)
  INTEGER :: istatus

  INTEGER :: i,j
  LOGICAL :: existfile
  character(len=256) :: ch1

  character(len=256) :: filename
  Doubleprecision :: ll_lon,ll_lat,lon_inc,lat_inc
  INTEGER :: imissing

  IF( nlat .ne. 1200 .or. nlon.ne. 1200 ) then
    print*,'Array size for array h_block, nlat and nlon must be '
    print*,'equal to 1200 inside read_30s_trn. '
    print*,'Input nlat and nlon =', nlat, nlon
    istatus = 1
    RETURN
  END IF

  istatus=1

  print*,'inside read_30s_trn, trn_fullname=', trim(trn_fullname)

  INQUIRE(FILE=trim(trn_fullname),EXIST=existfile)

  IF(existfile) THEN  ! unzipped version exists
    ch1='cp '//trim(trn_fullname)//' tmp1.file'
    CALL unixcmd(ch1)
  ELSE  ! need to unzip
    INQUIRE(FILE=trim(trn_fullname)//'.gz',EXIST=existfile)
    IF(existfile) then
      ch1='gunzip -c '//trim(trn_fullname)//'>tmp1.file'
      PRINT *,'Unzipping the terrain file..'
      CALL unixcmd(ch1)
      PRINT *,'Terrain file unzipping competed..'
    ELSE
      print*,'File ',trim(trn_fullname),' or its gzipped version not found,'
      print*,'Terrain data read failed.'
      istatus = 1
      STOP
      RETURN
    ENDIF
  END IF
!
! Read in terrain data
!
  call read_10x10_topo30s('tmp1.file',h_block,ll_lon, ll_lat,  &
       lon_inc, lat_inc,imissing)
  write(6,'(4(f20.14),3i7)') ll_lon, ll_lat, lon_inc, lat_inc, imissing
!
  CALL unixcmd('/bin/rm tmp1.file')

  istatus = 0
!
  RETURN
END SUBROUTINE read_30s_trn

SUBROUTINE read_10x10_topo30s(filename,it1,ll_lon,ll_lat, &
           lon_inc,lat_inc,imissing)
!
! Read in 10x10 degrees 30 second encoded ASCII terrain file filename
! and store the values in it1
!
! Author: Ming Xue (1/31/2002)
!
  IMPLICIT none
  character (len=*) :: filename
  INTEGER :: it1(1200,1200)
  Doubleprecision :: ll_lon,ll_lat,lon_inc,lat_inc
  INTEGER :: imissing

  character(len=1), ALLOCATABLE :: ichar1(:,:)
  character(len=1), ALLOCATABLE :: ichar2(:,:)

  INTEGER :: itopomax, itopomin
  INTEGER :: itopomax1,itopomin1
  INTEGER :: i,j,ii
  real :: topomean1

  ALLOCATE(ichar1(1200,1200),ichar2(1200,1200))

  open(21,file=trim(filename),status='old',form='formatted')
  read(21,'(4(f20.14),3i7)')ll_lon,ll_lat,lon_inc,lat_inc, &
       itopomin,itopomax,imissing

  do j=1,1200
  do i=1,1200,20
    read(21,'(40a1)') (ichar1(i-1+ii,j),ichar2(i-1+ii,j),ii=1,20)
  enddo
  enddo

  close(21)

  do j=1,1200
  do i=1,1200
    it1(i,j) = (ichar(ichar1(i,j))-32)*110+ ichar(ichar2(i,j))-32 - 999
    if( it1(i,j) .eq. -999) it1(i,j) = -9999
  enddo
  enddo

  itopomin1 =  100000
  itopomax1 = -100000
  topomean1 = 0.0
  do j=1,1200
  do i=1,1200
    if( it1(i,j).ne.-9999 .and. it1(i,j).gt.itopomax1) itopomax1=it1(i,j)
    if( it1(i,j).ne.-9999 .and. it1(i,j).lt.itopomin1) itopomin1=it1(i,j)
    if( it1(i,j).ne.-9999) topomean1 = topomean1 + it1(i,j)
  enddo
  enddo
  topomean1 = topomean1/(1200*1200)

  if( itopomin1 .eq. 100000 ) itopomin1 = -9999
  if( itopomax1 .eq.-100000 ) itopomax1 = -9999

  write(6,'(a,2i7,f20.14)')  &
          'Min., Max and mean of the 10x10 degree block are ', &
          itopomin1, itopomax1, topomean1

  IF( itopomin1 .ne. itopomin1 .or.itopomax1.ne.itopomax ) then
    print*,'Warning: min or max in data read in do not agree with '
    print*,'those in the header. Min., Max in thedata header were ',  &
            itopomin, itopomax
  ENDIF

  DEALLOCATE(ichar1,ichar2)

END SUBROUTINE read_10x10_topo30s


SUBROUTINE plot_trn(z,x,y,m,m1,m2,n,n1,n2,x1,x2,xstep,                  &
                    y1,y2,ystep,zinc,tem,xw,yw)

  IMPLICIT NONE

  INTEGER :: m,m1,m2,n,n1,n2
  REAL    :: x1,x2,xstep,y1,y2,ystep
  REAL    :: z(m,n),x(m,n),y(m,n),cl(100)

  REAL    :: tem(m,n),xw(8*m),yw(8*m)
  REAL    :: zinc,zmin, zmax, amax, amin

  CHARACTER (LEN=80) :: ch
  INTEGER :: lch, mode, ncl
  REAL    :: pl, pr, pb, pt, px, py, pxc, pyc, xs, ys

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ncl = 0     ! initilize

!
! Setup device and define plotting space.
!
  pl = 0.15
  pr = 0.85
  pb = 0.15
  pt = 0.85
  px = pr - pl
  py = pt - pb
  pxc = (pr+pl)/2
  pyc = (pb+pt)/2
  xs = x2-x1
  ys = y2-y1

  IF( py/px >= ys/xs ) THEN
    py = ys/xs*px
    CALL xpspac(pl, pr, pyc-py/2, pyc+py/2 )
  ELSE
    px = xs/ys*py
    CALL xpspac(pxc-px/2, pxc+px/2, pb, pt)
  END IF

  CALL xmap(x1,x2,y1,y2)
!
! Call XCONTA and XCLIMT to plot a contour map with default options.
!
  cl(1)=0.0
  cl(2)=zinc
  CALL xnctrs(10,100 )
  mode=1

  PRINT*,'m,m1,m2, n,n1,n2=' , m,m1,m2, n,n1,n2

  CALL a3dmax0(z,1,m,m1,m2,1,n,n1,n2,1,1,1,1, amax,amin)
  WRITE(6,'(1x,2(a,e13.6))')                                            &
      'Min. terrain= ', amin,', Max. terrain=',amax

  CALL xctrclr(46,68)
  CALL xctrlim(0.0, amax)
  CALL xctrlim(0.0, 0.0)

! call xcontcopt(2)

  CALL xcolfil(z(m1,n1),x(m1,n1),y(m1,n1),                              &
       tem,xw,yw,m,m2-m1+1,n2-n1+1, cl, ncl,mode)
  CALL xchsiz( 0.025 * (y2-y1))
  CALL xcpalet(2,-9999.9, -9999.0)

  CALL xclfmt('(I4)')
  CALL xctrclr(1,1)
!  CALL XCONTA(Z(m1,n1),X(m1,n1),Y(m1,n1),tem,M,m2-m1+1,n2-n1+1,
!    :            CL, NCL,MODE)
!

  IF(ncl <= 0 .OR. ncl > 100) ncl = 2  ! Avoid BT/TRAP on sooner
  zmax=cl(ncl)
  zmin=cl(1)
  zinc=cl( MIN(2,ncl) )-cl(1)
  CALL xchsiz(0.03*(y2-y1))
  WRITE(ch,'(''HMIN='',F6.1,'' HMAX='',F6.1,'' INC='',F5.1)')           &
        amin, amax, zinc
  lch = 80
  CALL strlnth(ch,lch)
  CALL xcolor(1)

  CALL xcharc((x1+x2)*0.5,y2+0.01*(y2-y1),ch(1:lch))
!
! Draw a border and plot ticks and labels on the border.
!
  CALL xaxsca(x1,x2,xstep,y1,y2,ystep)

  RETURN
END SUBROUTINE plot_trn

SUBROUTINE smth(a,m,m1,m2,n,n1,n2, ismooth,tem)

! Two-Dimension smooth use nine-points
  IMPLICIT NONE
  INTEGER :: m,n,m1,m2,n1,n2
  INTEGER :: ismooth

  INTEGER :: i,j,im,ip,jm,jp
  REAL :: s
  REAL :: a(m,n)
  REAL :: tem(m,n)

  IF(ismooth == 0) RETURN

  DO j=1,n
    DO i=1,m
      tem(i,j) = a(i,j)
    END DO
  END DO

  IF(ismooth == 1) THEN ! 5 point smoothing

    PRINT*,'Performing 5-point smoothing on input array'

    s=0.5
    DO j=n1,n2
      DO i=m1,m2
        im=MAX(m1,i-1)
        ip=MIN(m2,i+1)
        jm=MAX(n1,j-1)
        jp=MIN(n2,j+1)
        tem(i,j)=a(i,j)*(1.-s) +                                        &
            ( a(ip,j)+a(im,j)+a(i,jp)+a(i,jm) )*0.25*s
      END DO
    END DO

  ELSE IF(ismooth == 2) THEN  ! 9 point smoothing

    PRINT*,'Performing 9-point smoothing on input array'

    s=0.5
    DO j=n1,n2
      DO i=m1,m2
        im=MAX(m1,i-1)
        ip=MIN(m2,i+1)
        jm=MAX(n1,j-1)
        jp=MIN(n2,j+1)
        tem(i,j)=a(i,j)*(1.-s)**2 +                                     &
            ( a(ip,j)+a(im,j)+a(i,jp)+a(i,jm) )*0.5*s*(1.-s)            &
            + ( a(im,jm)+a(im,jp)+a(ip,jm)+a(ip,jp) )                   &
            *0.25*s**2
      END DO
    END DO

  ELSE IF(ismooth == 3) THEN  ! horizontal smoothing only

    PRINT*,'Performing horizontal smoothing on input array'

    s=0.5
    DO j=n1,n2
      DO i=m1,m2
        im=MAX(m1,i-1)
        ip=MIN(m2,i+1)
        tem(i,j)=a(i,j)*(1.-s)+(a(ip,j)+a(im,j))*0.5*s
      END DO
    END DO

  ELSE IF(ismooth == 4) THEN  ! vertical smoothing only

    PRINT*,'Performing vertical smoothing on input array'
    s=0.5
    DO j=n1,n2
      DO i=m1,m2
        jm=MAX(n1,j-1)
        jp=MIN(n2,j+1)
        tem(i,j)=a(i,j)*(1.-s)+(a(i,jp)+a(i,jm))*0.5*s
      END DO
    END DO

  END IF


  DO j=n1,n2
    DO i=m1,m2
      a(i,j) = tem(i,j)
    END DO
  END DO


  RETURN
END SUBROUTINE smth

SUBROUTINE Fill_missing(inblock,nlat,nlon,i,j,outval,outdist,istatus)

!#######################################################################
!
! Fill the missing value at grid (i,j) from surrounding valid observations
!
!-----------------------------------------------------------------------
!
! AUTHOR: Y. Wang (09/11/2012)
!
!-----------------------------------------------------------------------
!
! MODIFICATIONS:
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

  INTEGER,             INTENT(IN)  :: nlat, nlon
  INTEGER(KIND=INT16), INTENT(IN)  :: inblock(nlon,nlat)
  INTEGER,             INTENT(IN)  :: i,j
  INTEGER,             INTENT(OUT) :: outval
  INTEGER,             INTENT(OUT) :: istatus
  INTEGER,             INTENT(OUT) :: outdist
!-----------------------------------------------------------------------

  INTEGER :: badtotal, badnum
  INTEGER :: sqrtdist

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  badtotal = 0
  badnum   = 0

  outdist  = 0

  DO WHILE(badnum <= 0)

    outdist  = outdist + 1

    sqrtdist = 0   ! first check no sqrt distance

    IF (i > outdist .AND. inblock(i-outdist,j) /= -32768) THEN
      badtotal = badtotal + inblock(i-outdist,j)
      badnum = badnum + 1
    END IF

    IF (i <= nlon-outdist .AND. inblock(i+outdist,j) /= -32768) THEN
      badtotal = badtotal + inblock(i+outdist,j)
      badnum = badnum + 1
    END IF

    IF (j > outdist .AND. inblock(i,j-outdist) /= -32768) THEN
      badtotal = badtotal + inblock(i,j-outdist)
      badnum = badnum + 1
    END IF

    IF (j <= nlat-outdist .AND. inblock(i,j+outdist) /= -32768) THEN
      badtotal = badtotal + inblock(i,j+outdist)
      badnum = badnum + 1
    END IF

    IF (badnum <= 0) THEN   ! Another set of 4 points with distance sqrt(2)*dx

      sqrtdist = 1  ! now check

      IF (i > outdist .AND. j > outdist .AND. inblock(i-outdist,j-outdist) /= -32768) THEN
        badtotal = badtotal + inblock(i-outdist,j-outdist)
        badnum = badnum + 1
      END IF

      IF (i <= nlon-outdist .AND. j > outdist .AND. inblock(i+outdist,j-outdist) /= -32768) THEN
        badtotal = badtotal + inblock(i+outdist,j-outdist)
        badnum = badnum + 1
      END IF

      IF (i > outdist .AND. j <= nlat-outdist .AND. inblock(i-outdist,j+outdist) /= -32768) THEN
        badtotal = badtotal + inblock(i-outdist,j+outdist)
        badnum = badnum + 1
      END IF

      IF (i <= nlon-outdist .AND. j <= nlat-outdist .AND. inblock(i+outdist,j+outdist) /= -32768) THEN
        badtotal = badtotal + inblock(i+outdist,j+outdist)
        badnum = badnum + 1
      END IF

    END IF

  END DO

  IF (badnum > 0) THEN ! use nearest neighbor or average
    !IF (sqrtdist > 0) THEN
    !  WRITE(*,'(1x,a,2(I5,a),I3,a,I0,a)') '(',i,',',j,') Find ',badnum, &
    !           ' valid points at distance SQRT(2)*',outdist,'.'
    !ELSE
    !  WRITE(*,'(1x,a,2(I5,a),I3,a,I0,a)') '(',i,',',j,') Find ',badnum, &
    !           ' valid points at distance ',outdist,'.'
    !END IF
    outval = badtotal / badnum
  END IF

  RETURN
END SUBROUTINE Fill_missing
