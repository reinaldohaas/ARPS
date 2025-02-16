!
!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM ARPSSFC20                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM arpssfc20
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Prepare surface data file for ARPS model.
!
!  The following data sets will be generated:
!
!  soiltyp --  Soil type, translated from Wilson, M. F. and
!              Henderson-Sellers soil type data set, array styp.
!              Array styp was transferred from 1-byte data file to
!              4-integer ASCII file without changing the definition.
!              (Original source: 1-degree 180x360 GED in CD-ROM)
!
!  vegtyp  --  Vegetation type, translated from Olson, J. S., World
!              Ecosystems, array vtyp. The array vtyp was
!              transferred from 1-byte data file to 4-integer ASCII
!              file without changing the definition.
!              (Original source: 10-minute 1080x2160 GED in CD-ROM)
!
!  lai     --  Leaf Area Index, calculated from the vegetation type,
!              vegtyp and NDVI which is provided by Gallo, K. P.,
!              Normalized Difference Vegetation Index in GED CD-ROM
!              and transferred from 1-byte data file to 4-integer
!              ASCII file without changing the definition.
!              (Original NDVI source: 10-minute 1080x2160 GED in
!               CD-ROM)
!
!  roufns  --  Surface roughness which will be obtained from a data
!              table vs. vegtyp.
!
!  veg     --  Vegetation fraction which will be obtained from a data
!              table vs. vegtyp or the NESDIS green vegetation data
!              set (Gutman and Ignatov, 1998).
!
!  This program will put the these data sets into the ARPS
!  model grid points. You may chose one of three different map
!  projections.
!
!  The structure of this program is as follows:
!
!  1. Determine the model domain and set up the analysis grid.
!     (subroutine setgrd)
!
!  2. Read in the surface data set.
!     (subroutine gtsfcdt)
!
!  3. Calculate:
!
!     soil type (soiltyp);
!     vegetation type (vegtyp);
!     Leaf Area Index (lai);
!     surface roughness (roufns).
!     vegetation fraction (veg).
!
!     (subroutine gtsoiltyp, gtvegtyp, gtndvi, and gtlai, gtrfns, and
!      gtveg)
!
!  4. Output soiltyp, vegtyp, lai, roufns, and veg in ARPS model
!     domain into the file sfcoutfl by calling WRTSFCDT.
!
!     sfcoutfl is constructed to be runname(1:lfnkey)//".sfcdata
!
!  5. Plot the data for diagnostic purpose.
!     (subroutine plot)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  02/20/94
!
!  MODIFICATIONS:
!
!  01/25/1995 (Y. Liu)
!  Added the scheme that produce an artificial distribution of soil
!  and vegetation parameters in two specified regions in ARPS domain.
!
!  02/07/1995 (Y. Liu)
!  Added new 2-D array, veg, to the data set file.
!
!  05/08/1996 (Yuhe Liu)
!  Added an adjustment to make soil type and vegetation type
!  consistent with each other for ice and water surface.
!
!  02/20/1997 (Gene Bassett)
!  Added look-up option for setting veg, lai, & roufns.  Added nsmthsl,
!  the number of smoothing passes to make, and altered code to average
!  the LOG of roufns.
!
!  04/08/1997 (Leilei Wang and Vince Wong)
!  Version 2.0 for STATSGO 1 km soil texture data.
!
!    REFERENCE:
!    Miller, D.A. and R.A. White. 1996. A soils dataset for
!    environmental modeling applications in the United States.
!    (To be submitted to Earth Interactions).
!
!  08/11/1997 (Leilei Wang and Vince Wong)
!  Introduce 1km North American Land Cover Characteristic Data Set
!  (Vegetation type data set and NDVI data set) into ARPS.
!  Use GED data as background field in area outside US or North
!  American.
!
!  10/10/1997 (Yuhe Liu)
!  Added a call to print integer value for integer arrays at each
!  grid point
!
!  04/15/1998 (Dan Weber)
!  Added subroutine vegfrac to read in the vegetation fraction data
!  set from  NESDIS.  This data set will overwrite any pre-existing
!  value when data are present.  The pre-existing data is generated
!  using the veg type/table conversion.
!
!  07/16/1999 (Gene Bassett)
!  Corrected dimensions passed into get_ged, getstyp2 & getndvi2.
!
!  2000/08/23 (Gene Bassett)
!  Added F90 allocation of arrays.
!
!  2001/2/26 (Ming Xue)
!  Cleaned up the code and reduced temporary storage usage.
!
!  2001/6/20 (Yunheng Wang)
!  Modify main program and subroutine to reduce memory usage.
!
!  2001/07/26 (Gene Bassett)
!  Moved calls to getstyp1, getvtyp1, getndvi1 inside get_statsgo
!  and calls to getstyp2 getvtyp2 getvtyp2 inside get_ged to
!  eliminate input variables colst,rowst,colvt,rowvt,colnv,&rownv.
!
!  14 June 2002 (Eric Kemp)
!  Bug fix with placement of nstyp in namelist.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  Input read in from the input file, namelist &projection
!
!  mapproj  Type of map projection used to setup the analysis grid.
!           mapproj = 1  Polar Stereographic projection
!                   = 2  Mercator projection
!                   = 3  Lambert
!
!  trulat1  1st real true latitude of map projection.
!
!  trulat2  2nd real true latitude of map projection.
!
!  trulon   Real true longitude of map projection.
!
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!           distance on map = (Distance on earth) * sclfct
!           For ARPS model runs, generally this is 1.0.
!
!  Input read in from the input file, namelist &modgrid
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!
!  ctrlat   Latitude of the southwest corner of the analysis grid
!           (deg. N)
!  ctrlon   Longitude of the southwest corner of the analysis grid
!           (deg. E)
!
!  Input read in from the input file, namelist &outflag
!
!  stypout  Flag for output of soil type
!  vtypout  Flag for output of vegetation type
!  laiout   Flag for output of Leaf Area Index
!  rfnsout  Flag for output of surface roughness
!  vegout   Flag for output of vegetation fraction
!  ndviout  Flag for output of NDVI
!
!  OUTPUT:
!
!  Output written to the surface data file sfcdtfl:
!
!  nx       Number of model grid points in the x-direction
!  ny       Number of model grid points in the y-direction
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat   Lat. at the origin of the model grid (deg. N)
!  ctrlon   Lon. at the origin of the model grid (deg. E)
!
!  stypout  Flag for output of soiltyp
!  vtypout  Flag for output of vegtyp
!  laiout   Flag for output of Leaf Area Index
!  rfnsout  Flag for output of surface roughness
!  vegout   Flag for output of vegetation fraction
!  ndviout  Flag for output of vegetation fraction
!
!  soiltyp  Soil type in model domain
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!
!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT   NONE

  INCLUDE   'arpssfc.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil
  INTEGER :: nstyps

  INTEGER, PARAMETER :: scol = 360,  srow=180   ! GED data set
  INTEGER, PARAMETER :: vncol = 2160,vnrow=1080 ! STATSGO north american 1km data set

  CHARACTER (LEN=65 ) :: glab     ! label for plots

  CHARACTER (LEN=256) :: fstypfl  ! File name of soil data set
  CHARACTER (LEN=256) :: bstypfl  ! File name of soil data set

  CHARACTER (LEN=256) :: fvtypfl  ! File name of Vegetation data set
  CHARACTER (LEN=256) :: bvtypfl  ! File name of Vegetation data set

  CHARACTER (LEN=256) :: fndvifl  ! File name of foreground NDVI data
  CHARACTER (LEN=256) :: bndvifl  ! File name of background NDVI data

  CHARACTER (LEN=256) :: lkupfl   ! File name of look-up file
  CHARACTER (LEN=256) :: vfrcdr   ! Directory name for vegetation
                                  ! fraction data file.

  INTEGER :: schmopt   ! Scheme of options

  INTEGER :: sdatopt   ! Soil type data  options
  INTEGER :: vdatopt   ! Vegetation type data  options
  INTEGER :: ndatopt   ! NDVI data  options
  INTEGER :: vfrcopt   ! Vegetation fraction option

  INTEGER :: nsmthsl   ! Number of passes to smooth the data

  INTEGER :: fgbgni    ! First index in x-direction in foreground region
  INTEGER :: fgendi    ! Last  index in x-direction in foreground region
  INTEGER :: fgbgnj    ! First index in y-direction in foreground region
  INTEGER :: fgendj    ! Last  index in y-direction in foreground region

  INTEGER :: fgstyp    ! Soil type for the foreground
  INTEGER :: fgvtyp    ! Vegetation type for the foreground
  REAL :: fglai        ! Leaf Area Index for the foreground
  REAL :: fgrfns       ! Surface roughness for the foreground
  REAL :: fgveg        ! Vegetation fraction for the foreground

  INTEGER :: bgstyp    ! Soil type for the background
  INTEGER :: bgvtyp    ! Vegetation type for the background
  REAL    :: bglai     ! Leaf Area Index for the background
  REAL    :: bgrfns    ! Surface roughness for the background
  REAL    :: bgveg     ! Vegetation fraction for the background

  INTEGER :: drawval   ! Option to draw integer value at each grid point

!  INTEGER, ALLOCATABLE :: ged_styp(:,:)! Background soil types
!  INTEGER, ALLOCATABLE :: ged_vtyp(:,:)! Background vegetation types
!  REAL,    ALLOCATABLE :: ged_ndvi(:,:)! Background NDVI

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type in model domain
  REAL,    ALLOCATABLE :: stypfrct(:,:,:) ! Fraction of soil types
  INTEGER, ALLOCATABLE :: vegtyp (:,:) ! Vegetation type in model domain
  REAL,    ALLOCATABLE :: ndvi   (:,:) ! NDVI in model domain
  REAL,    ALLOCATABLE :: lai    (:,:) ! Leaf Area Index in model domain
  REAL,    ALLOCATABLE :: roufns (:,:) ! Surface roughness in model domain
  REAL,    ALLOCATABLE :: veg    (:,:) ! Vegetation fraction in model domain

  REAL,    ALLOCATABLE :: vegin1 (:,:) ! Green Vegetation fraction data.
                                       ! for use in reading in the NESDIS
                                       ! green vegetation data set. month#1
  REAL,    ALLOCATABLE :: vegin2 (:,:) ! Green Vegetation fraction data.
                                       ! for use in reading in the NESDIS
                                       ! green vegetation data set. month#2

  REAL, ALLOCATABLE :: x    (:)        ! Model grid point values in x-dir.
  REAL, ALLOCATABLE :: y    (:)        ! Model grid point values in y-dir.
  REAL, ALLOCATABLE :: x1   (:,:)      ! Model grid point values in x-dir.
  REAL, ALLOCATABLE :: y1   (:,:)      ! Model grid point values in y-dir.
  REAL, ALLOCATABLE :: xs   (:)        ! X coordinates for scalar fields
  REAL, ALLOCATABLE :: ys   (:)        ! Y coordinates for scalar fields

  REAL, ALLOCATABLE :: glat (:,:)      ! Lat. of model grid
  REAL, ALLOCATABLE :: glon (:,:)      ! Lon. of model grid

  REAL, ALLOCATABLE :: vegtbl(:)       ! veg as a fuction of vtype.
  REAL, ALLOCATABLE :: laitbl(:)       ! lai as a fuction of vtype.
  REAL, ALLOCATABLE :: z0tbl (:)       ! z0 (:) as a fuction of vtype.

  INTEGER, ALLOCATABLE :: item3(:,:)   ! Temporary array
  REAL,    ALLOCATABLE :: rtem2(:,:)   ! Temporary array
  REAL,    ALLOCATABLE :: rtem3(:,:)   ! Temporary array

  INTEGER :: itmp(1)            ! Temporary array pointer
  REAL    :: rtmp(1)            ! Temporary array pointers

  !
  ! Added 30s data support (Y. Wang)
  !
  INTEGER              :: vcol,vrow
                               ! 30s USGS landuse, will be determined below
  !TINA, Andreas: Included 30s option by Queenie
  INTEGER, ALLOCATABLE :: v30sdat(:,:)  ! 30 second data array
  REAL,    ALLOCATABLE :: v30slat(:)    ! Lat of 30s data
  REAL,    ALLOCATABLE :: v30slon(:)    ! Lon of 30s data
  INTEGER, ALLOCATABLE :: tem19(:,:)    ! temporary array
  ! Andreas - end
  INTEGER              :: latbgn, latend, lonbgn, lonend
  REAL                 :: latmax, latmin, lonmax, lonmin

!-----------------------------------------------------------------------
!  Misc. local variables:
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: sfcoutfl,  temchar
  INTEGER :: lfn

  INTEGER :: i,j, is, n
  INTEGER :: lenstr

  REAL :: resl

  REAL :: alatpro(2)           ! Latitude at which the map projection
                               ! is true (degees N)
  REAL :: alonpro              ! Longitude at which the map
                               ! projection is true (degrees E)
  REAL :: sclf                 ! Map scale factor (m)
  REAL :: cint                 ! Contour interval

  REAL    :: tema, temb, temc, temd, t1, t2
  INTEGER :: itema,itemb

  CHARACTER(LEN=256) :: namelist_filename

!-----------------------------------------------------------------------
! Include files:
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'ext2arps.inc'  ! Not actually needed by this program.

!-----------------------------------------------------------------------
! Namelists:
!-----------------------------------------------------------------------

  NAMELIST /soil_veg_data/ schmopt,sdatopt,vdatopt,ndatopt,vfrcopt,     &
                           nsmthsl,fgbgni,fgendi,fgbgnj,                &
                           fgendj,fgstyp,fgvtyp,fglai,fgrfns,fgveg,     &
                           bgstyp,bgvtyp,bglai,bgrfns,bgveg,            &
                           stypout,vtypout,laiout,                      &
                           rfnsout,vegout,ndviout,drawval,              &
                           fstypfl,bstypfl,fvtypfl,bvtypfl,             &
                           fndvifl,bndvifl,lkupfl,vfrcdr

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!  Set default value of nsmthsl:
  nsmthsl = 1
  drawval = 0

  namelist_filename = ' '
  CALL initpara(nx,ny,nz,nzsoil,nstyps,namelist_filename)

  nstyp    = nstyps

!wdt
!  CALL initadas ! The adas parameters are not used by this program,
!                ! but the name list blocks need to be read in sequence on
!                ! Cray machines such as J90.
!
!  READ(5, extdfile) ! These ext2arps parameters are not nssed by this program.
!                    ! They are read for the same reason given above.
!  WRITE(6,'(a)')'Namelist block extdfile sucessfully read.'

  IF (sfcdmp == 0) sfcdmp = 1 ! Force the dumping of surface data
                              ! since sfcdmp defaults to 0 in initpara.

  alatpro(1) = trulat1
  alatpro(2) = trulat2
  alonpro    = trulon

  WRITE (6, '(a)') 'Read the namelist file.'

  READ (5, soil_veg_data)

  IF ( schmopt < 0 .OR. schmopt > 4 ) THEN
    WRITE (6,'(a,i2/a)')                                                &
        'Unrecognized scheme option: ', schmopt,                        &
        'Program stopped'
    STOP
  END IF

  IF ( schmopt > 1 .AND. (sdatopt < 1 .OR. sdatopt > 3) ) THEN
                                        !TINA, Andreas: allow sdatopt 3
    WRITE (6,'(a,i2/a)')                                                &
        'Unrecognized data option: ', sdatopt,'Program stopped'
    STOP
  END IF

  IF ( nstyp <= 0 ) THEN
    nstyp = 1
    WRITE (6,'(a,i2)')                                                  &
        'The nstyp must be greater than 0. Reset it to 1.'
  ELSE IF ( nstyp > nstyps ) THEN
    nstyp = nstyps
    WRITE (6,'(a,a,i2)')                                                &
        'The nstyp must be less than, or equal to nstyps. ',            &
        'Reset it to nstyps = ', nstyps
  END IF

  WRITE (6, '(/1x,a)')      '&soil_veg_data'
  WRITE (6, '(3x,a,i4,a)')    'schmopt = ', schmopt,','
  WRITE (6, '(3x,a,i4,a)')    'sdatopt = ', sdatopt,','
  WRITE (6, '(3x,a,i4,a)')    'vdatopt = ', vdatopt,','
  WRITE (6, '(3x,a,i4,a)')    'ndatopt = ', ndatopt,','
  WRITE (6, '(3x,a,i4,a)')    'vfrcopt = ', vfrcopt,','
  WRITE (6, '(3x,a,i4,a)')    'nstyp   = ', nstyp,  ','
  WRITE (6, '(3x,a,i4,a)')    'nsmthsl = ', nsmthsl,','
  WRITE (6, '(3x,a,i4,a)')    'fgbgni  = ', fgbgni, ','
  WRITE (6, '(3x,a,i4,a)')    'fgendi  = ', fgendi, ','
  WRITE (6, '(3x,a,i4,a)')    'fgbgnj  = ', fgbgni, ','
  WRITE (6, '(3x,a,i4,a)')    'fgendj  = ', fgendi, ','
  WRITE (6, '(3x,a,i4,a)')    'fgstyp  = ', fgstyp, ','
  WRITE (6, '(3x,a,i4,a)')    'fgvtyp  = ', fgvtyp, ','
  WRITE (6, '(3x,a,f10.3,a)') 'fglai   = ', fglai,  ','
  WRITE (6, '(3x,a,f10.3,a)') 'fgrfns  = ', fgrfns, ','
  WRITE (6, '(3x,a,f10.3,a)') 'fgveg   = ', fgveg,  ','
  WRITE (6, '(3x,a,i4,a)')    'bgstyp  = ', bgstyp, ','
  WRITE (6, '(3x,a,i4,a)')    'bgvtyp  = ', bgvtyp, ','
  WRITE (6, '(3x,a,f10.3,a)') 'bglai   = ', bglai,  ','
  WRITE (6, '(3x,a,f10.3,a)') 'bgrfns  = ', bgrfns, ','
  WRITE (6, '(3x,a,f10.3,a)') 'bgveg   = ', bgveg,  ','
  WRITE (6, '(3x,a,i4,a)')    'stypout = ', stypout,','
  WRITE (6, '(3x,a,i4,a)')    'vtypout = ', vtypout,','
  WRITE (6, '(3x,a,i4,a)')    'laiout  = ', laiout, ','
  WRITE (6, '(3x,a,i4,a)')    'rfnsout = ', rfnsout,','
  WRITE (6, '(3x,a,i4,a)')    'vegout  = ', vegout, ','
  WRITE (6, '(3x,a,i4,a)')    'ndviout = ', ndviout,','

  lenstr = LEN( fstypfl )
  CALL strlnth( fstypfl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'fstypfl  = ''', fstypfl(1:lenstr), ''','

  lenstr = LEN( fvtypfl )
  CALL strlnth( fvtypfl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'fvtypfl  = ''', fvtypfl(1:lenstr), ''','

  lenstr = LEN( fndvifl )
  CALL strlnth( fndvifl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'fndvifl  = ''', fndvifl(1:lenstr), ''','

  lenstr = LEN( bstypfl )
  CALL strlnth( bstypfl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'bstypfl  = ''', bstypfl(1:lenstr), ''','

  lenstr = LEN( bvtypfl )
  CALL strlnth( bvtypfl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'bvtypfl  = ''', bvtypfl(1:lenstr), ''','

  lenstr = LEN( bndvifl )
  CALL strlnth( bndvifl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'bndvifl  = ''', bndvifl(1:lenstr), ''','

  lenstr = LEN ( lkupfl )
  CALL strlnth( lkupfl, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'lkupfl  = ''', lkupfl(1:lenstr), ''','

  lenstr = LEN ( vfrcdr )
  CALL strlnth( vfrcdr, lenstr)
  WRITE (6, '(3x,a,a,a)')                                               &
            'vfrcdr  = ''', vfrcdr(1:lenstr), ''','

  WRITE (6, '(3x,a,i4,a)') '&END'

!-----------------------------------------------------------------------
! Allocate arrays.
!-----------------------------------------------------------------------

!  ALLOCATE(ged_styp(nx,ny))   ! Not used, Why allocate them?
!  ALLOCATE(ged_vtyp(nx,ny))
!  ALLOCATE(ged_ndvi(nx,ny))
!
!  ged_styp =0
!  ged_vtyp =0
!  ged_ndvi =0

  ALLOCATE(soiltyp(nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp(nx,ny))
  ALLOCATE(ndvi(nx,ny))

  soiltyp = 0
  stypfrct = 0
  vegtyp = 0
  ndvi = 0

  ALLOCATE(lai(nx,ny))
  ALLOCATE(roufns(nx,ny))
  ALLOCATE(veg(nx,ny))

  lai    = 0
  roufns = 0
  veg    = 0

  ALLOCATE(x(nx))
  ALLOCATE(y(ny))
  ALLOCATE(x1(nx,ny))
  ALLOCATE(y1(nx,ny))
  ALLOCATE(xs(nx))
  ALLOCATE(ys(ny))

  ALLOCATE(glat(nx,ny))
  ALLOCATE(glon(nx,ny))

  ALLOCATE(vegtbl(nvegtyp))
  ALLOCATE(laitbl(nvegtyp))
  ALLOCATE(z0tbl(nvegtyp))

  ALLOCATE(item3(nx,ny))
  ALLOCATE(rtem2(nx,ny))
  ALLOCATE(rtem3(nx,ny))

!-----------------------------------------------------------------------
!  Set the grid for the model and the extended buffer area.
!-----------------------------------------------------------------------

  CALL setgrd( nx,ny, x, y )

  DO i = 1,nx
    xs(i) = x(i) + dx/2  ! for scalar coordinates
  END DO

  DO j = 1,ny
    ys(j) = y(j) + dy/2  ! for scalar coordinates
  END DO

  CALL xytoll( nx,ny, xs,ys, glat,glon )

!-----------------------------------------------------------------------
!
! Specific allocation for 30s soil type or 30s USGS landuse options
!
!-----------------------------------------------------------------------

  IF (sdatopt == 3 .OR. vdatopt == 3) THEN

    IF (glon(nx,1) < glon(1,1) .AND. glon(nx,1) < 0) THEN
                                ! International date line in between
      DO j = 1, ny
        DO i = 1, nx
          glon(i,j) = MOD( glon(i,j) + 360., 360. )
        END DO
      END DO
    END IF

    latmax = MAXVAL(glat)/10
    latmin = MINVAL(glat)/10
    lonmax = MAXVAL(glon)/10
    lonmin = MINVAL(glon)/10

    latbgn = FLOOR   (latmin)*10
    latend = (CEILING(latmax)-1)*10
    lonbgn = FLOOR   (lonmin)*10
    lonend = (CEILING(lonmax)-1)*10

    vrow = (latend-latbgn)*120 + 1200 ! 30s interval, 1minute = 2 points
    vcol = (lonend-lonbgn)*120 + 1200 ! 1degree = 60m x 2 = 120

    ALLOCATE(v30slat(vrow))         ! latitude  is row
    ALLOCATE(v30slon(vcol))         ! longitude is column
    ALLOCATE(v30sdat(vcol,vrow))
    ALLOCATE(tem19  (vcol,vrow))

  END IF

!-----------------------------------------------------------------------
!
!  Change longitude rang from -180 -- 180 to 0 -- 360.
!
!-----------------------------------------------------------------------

  DO j = 1, ny
    DO i = 1, nx
      glon(i,j) = MOD( glon(i,j) + 360., 360. )
    END DO
  END DO

!-----------------------------------------------------------------------
!  Check if the specified region is in the ARPS domain
!-----------------------------------------------------------------------

  IF ( schmopt == 1 .OR. schmopt == 2 ) THEN
    IF ( fgbgni > fgendi .OR. fgbgnj > fgendj ) THEN
      WRITE (6, '(a/a)')                                                &
          'The indeces for the specified region were not consistent.',  &
          'Program ARPSSFC12 stopped.'
      STOP
    END IF

    IF ( fgbgni < 1 .OR. fgendi > nx .OR.  fgbgnj < 1 .OR. fgendj > ny ) THEN

      fgbgni = MAX( fgbgni, 1  )
      fgendi = MIN( fgendi, nx )
      fgbgnj = MAX( fgbgnj, 1  )
      fgendj = MIN( fgendj, ny )

      WRITE (6, '(a/a/a,i4,4x,a,i4,4x,a,i4,4x,a,i4/a)')                 &
          'The specified foreground region was out of the ARPS domain,', &
          'The adjusted region will be:',                               &
          'fgbgni = ',fgbgni, 'fgendi = ',fgendi,                       &
          'fgbgnj = ',fgbgnj, 'fgendj = ',fgendj,                       &
          'Program ARPSSFC12 will continue.'

    END IF
  END IF

!-----------------------------------------------------------------------
! For different schmopt, determine the soil and vegetation type,
! roughness and vegetation fraction.
!-----------------------------------------------------------------------

  IF ( schmopt == 0 .OR. schmopt == 1 ) THEN ! User specified values

    IF (nstyp > 1) WRITE(6,'(/,a,I2,/)')     &
        'WARNING: nstyp was reset to 1 because of schmopt = ',schmopt

    nstyp = 1

    DO j = 1, ny
      DO i = 1, nx
        soiltyp(i,j,1)  = bgstyp
        stypfrct(i,j,1) = 1.0
        vegtyp (i,j) = bgvtyp
        lai    (i,j) = bglai
        roufns (i,j) = bgrfns
        veg    (i,j) = bgveg
      END DO
    END DO

  ELSE IF ( schmopt == 2 .OR. schmopt == 3 .OR. schmopt == 4 ) THEN

!-----------------------------------------------------------------------
! Determine the soil type for ARPS
!-----------------------------------------------------------------------

    IF ( stypout /= 0 )THEN

      IF ( sdatopt == 1 ) THEN

!-----------------------------------------------------------------------
!       Get the STATSGO soil texture data and store in a subarray stypdat
!       just enough to cover the model domain.
!  Translate the soil classes data storged in stypdat into the USDA soil
!  type categories. nstyp soil types with the highest percentages are
!  then brought to the model grid points, with each type being associated
!  with a value of soil type fraction.
!-----------------------------------------------------------------------

        ! set up stypfrct for the areas where STATSGO doesn't cover
        stypfrct(:,:,1) = 1.
        stypfrct(:,:,2:nstyp) = 0.

        CALL get_statsgo(1,fstypfl,bstypfl,nx,ny,nstyp,x1,y1,     &
                        scol,srow,glat,glon,resl,soiltyp,stypfrct)

      ELSE IF ( sdatopt == 2 ) THEN

!-----------------------------------------------------------------------
!  Get the GED soil type data
!-----------------------------------------------------------------------

        IF (nstyp > 1) WRITE(6,'(/,a,/)')     &
          'WARNING: nstyp was reset to 1 because GED soil type data was chosen.'

        nstyp = 1  ! 1 only for now for this data set.

!-----------------------------------------------------------------------
!  Get the original translated surface soil data to array styp.
!  Translate the soil classes data (stypdat) into the USDA soil type
!  category and fill them into the model grid array soiltyp.
!-----------------------------------------------------------------------

        CALL get_ged(1,bstypfl,nx,ny,scol,srow,resl,glon,glat,          &
                     soiltyp, rtmp )

        stypfrct(:,:,1) = 1.0

      ELSE IF ( sdatopt == 3 ) THEN

        IF (nstyp > 1) WRITE(6,'(/,a,/)')     &
          'WARNING: nstyp was reset to 1 because 30s global soil type data was chosen.'

        nstyp = 1

        !TINA, Andreas - included 30s option by Queenie

!-----------------------------------------------------------------------
! read 30sec global soil type file
!-----------------------------------------------------------------------

        CALL get_30s (1,fstypfl,latbgn,latend,lonbgn,lonend,            &
                      vcol,vrow,resl,v30slon,v30slat,v30sdat )
        CALL mapty30s(1,vcol,vrow,nx,ny,resl,v30slon,v30slat,glon,glat, &
                      v30sdat,soiltyp,tem19)
        ! set up stypfrct
        stypfrct(:,:,1) = 1.
        stypfrct(:,:,2:nstyp) = 0.

        DO j=1,ny
          DO i=1,nx
            IF ( soiltyp(i,j,1) < 1 .OR. soiltyp(i,j,1) > 13 ) THEN
              print *,'CCL print soiltyp',i,j,soiltyp(i,j,1)
              STOP
            END IF
          END DO
        END DO

        !TINA, Andreas - end

      END IF

    END IF

!-----------------------------------------------------------------------
!   Process vegetation type data
!-----------------------------------------------------------------------

    IF ( vtypout /= 0 .OR. laiout /= 0 .OR. rfnsout /= 0 ) THEN

      IF ( vdatopt == 1 ) THEN

!-----------------------------------------------------------------------
!  Get the USGS North American 1km vegetation data just covering the
!  model domain to array stypdat.
!  Calculate vegetation type array using statistical method to get
!  vegtyp(nx,ny)
!-----------------------------------------------------------------------

        CALL get_statsgo(2,fvtypfl,bvtypfl,nx,ny,nstyp,x1,y1,           &
                         vncol,vnrow,glat,glon,resl,vegtyp,rtmp)

      ELSE IF ( vdatopt == 2 ) THEN

!-----------------------------------------------------------------------
!  Get the original translated World Ecosystem Classes of vegetation
!  data to array vegtyp.
!  Transfer the World Ecosystem Classes data into the 12 vegetation
!  type categories and fill into model grid array, vegtyp(nx,ny)
!-----------------------------------------------------------------------

        CALL get_ged(2, bvtypfl,nx,ny,vncol,vnrow,resl,glon,glat,       &
                     vegtyp,rtmp)

      ELSE IF ( vdatopt == 3 ) THEN

!-----------------------------------------------------------------------
! read 30sec global vegetation file
!-----------------------------------------------------------------------

        !TINA, Andreas - include 30s option by Queenie

        CALL get_30s (2,fvtypfl,latbgn,latend,lonbgn,lonend,            &
                      vcol,vrow,resl,v30slon,v30slat,v30sdat)
        CALL mapty30s(2,vcol,vrow,nx,ny,resl,v30slon,v30slat,glon,glat, &
                      v30sdat,vegtyp,tem19)

        !TINA, Andreas - end


      END IF

      IF ( schmopt == 4 ) THEN

        CALL gtlkuptbl( lkupfl, vegtbl, laitbl, z0tbl )

      ELSE

        IF ( ndatopt == 1 ) THEN

!-----------------------------------------------------------------------
!  Get the USGS 1km monthly byteNDVI data just covering the
!  model domain to array stypdat.
!  Calculate real NDVI type data using statistical method to get
!  ndvi(nx,ny)
!-----------------------------------------------------------------------

          CALL get_statsgo(3,fndvifl,bndvifl,nx,ny,nstyp,x1,y1,   &
                          vncol,vnrow,glat,glon,resl,itmp,ndvi)

        ELSE IF ( ndatopt == 2 ) THEN

!-----------------------------------------------------------------------
!  Get the original translated byte-NDVI data to array ndvibyt.
!  Calculate the real-NDVI data from byte-NDVI and fill into model
!  grid array, ndvi(nx,ny)
!-----------------------------------------------------------------------

          CALL get_ged(3,bndvifl,nx,ny,vncol,vnrow,resl,glon,glat,  &
                       itmp,ndvi)

        END IF

      END IF

    END IF

    IF ( laiout /= 0 ) THEN

!-----------------------------------------------------------------------
!  Calculate the Leaf Area Index array, lai(nx,ny), from the NDVI
!  data and vegetation type data, in the model grid.
!-----------------------------------------------------------------------

      IF ( schmopt == 4 ) THEN

        CALL gtlaitbl( nx,ny, vegtyp, laitbl, lai )

      ELSE

        CALL getlai( nx,ny, vegtyp, ndvi, lai )

      END IF
    END IF

    IF ( rfnsout /= 0 ) THEN
!
!-----------------------------------------------------------------------
!  Get the surface roufness from the tabledata vs vegtyp.
!-----------------------------------------------------------------------
!
      IF ( schmopt == 4 ) THEN

        CALL gtrfnstbl( nx,ny, vegtyp, z0tbl, roufns )

      ELSE

        CALL getrfns( nx,ny, vegtyp, roufns )

      END IF

    END IF

    IF ( vegout /= 0 ) THEN
!
!-----------------------------------------------------------------------
!  Get the vegetation fraction from the tabledata vs vegtyp.
!
!  Allow the vegetation table/fraction conversion to proceed.
!  It will be overwritten with data from NESDIS if gvfrcopt=1, only
!  in the location where data are present.
!  Reference: Gutman and Ignotov, 1998 in press.
!-----------------------------------------------------------------------
!
      CALL getveg( nx,ny, vegtyp, veg )

      IF ( schmopt == 4 ) THEN

!  use the vegetation type conversion table method.

        CALL gtvegtbl( nx,ny, vegtyp, vegtbl, veg ) ! old code...

      END IF

      IF(vfrcopt == 1)THEN  ! over write the table derived values
                            ! where we have data.....

        ALLOCATE(vegin1(gvnx,gvny))
        ALLOCATE(vegin2(gvnx,gvny))
        vegin1 = 0
        vegin2 = 0

        CALL gvegfrac(nx,ny,vfrcdr,glat,glon,vegin1,vegin2,veg)

        DEALLOCATE(vegin1, vegin2)

      END IF

    END IF

!
!-----------------------------------------------------------------------
!  Smooth lai, log(roufns), and veg.
!-----------------------------------------------------------------------
!
    IF ( rfnsout /= 0 ) THEN

      DO j = 1,ny
        DO i = 1,nx
          rtem3(i,j) = LOG ( MAX(0.00001,roufns(i,j)) )
        END DO
      END DO

    END IF

    DO n = 1,nsmthsl

      IF ( laiout /= 0 ) THEN

        CALL smooth25p( lai, nx,ny,1,nx,1,ny, rtem2 )

      END IF

      IF ( rfnsout /= 0 ) THEN

        CALL smooth25p( rtem3, nx,ny,1,nx,1,ny, rtem2 )

      END IF

      IF ( vegout /= 0 ) THEN

        CALL smooth25p( veg, nx,ny,1,nx,1,ny, rtem2 )

      END IF

      IF ( ndviout /= 0 ) THEN

        CALL smooth25p( ndvi, nx,ny,1,nx,1,ny, rtem2 )

      END IF

    END DO

    IF ( rfnsout /= 0 ) THEN

      DO j = 1,ny
        DO i = 1,nx
          roufns(i,j) = EXP ( rtem3(i,j) )
        END DO
      END DO

    END IF

  ELSE
    WRITE (6, '(a,i4/a)') 'Unexpected schmopt: ',schmopt,               &
                        'Program ARPSSFC12 stopped here'
    STOP
  END IF
!
!-----------------------------------------------------------------------
! Set the data sets for the foreground region.
!-----------------------------------------------------------------------
!
  IF ( schmopt == 1 .OR. schmopt == 2 ) THEN
    DO j = fgbgnj, fgendj
      DO i = fgbgni, fgendi
        soiltyp(i,j,1) = fgstyp
        vegtyp (i,j) = fgvtyp
        lai    (i,j) = fglai
        roufns (i,j) = fgrfns
        veg    (i,j) = fgveg
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
! Eliminate the inconsistence over water between soil type and
! vegetation type
!-----------------------------------------------------------------------
!
  DO is=1,nstyp
    DO j=1,ny
      DO i=1,nx
        IF ( vegtyp(i,j) == 14 ) THEN
          soiltyp(i,j,is) = 13
        END IF
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
! Write out a surface property data file
!-----------------------------------------------------------------------
!
  sfcoutfl = runname(1:lfnkey)//".sfcdata"
  lfn = lfnkey + 7

  IF( dirname /= ' ' ) THEN

    temchar = sfcoutfl
    sfcoutfl = dirname(1:ldirnam)//'/'//temchar
    lfn  = lfn + ldirnam + 1

  END IF

  CALL fnversn(sfcoutfl, lfn)

  PRINT *, 'Write surface property data in ',sfcoutfl(1:lfn)

  CALL wrtsfcdt(nx,ny,nstyp, sfcoutfl,dx,dy,                            &
                mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,    &
                stypout,vtypout,laiout,rfnsout,vegout,ndviout,          &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
!
!-----------------------------------------------------------------------
! Creating the GrADS control file for data display
!-----------------------------------------------------------------------
!
  IF (sfcdmp == 1) CALL sfccntl( nx,ny, sfcoutfl,                       &
                stypout,vtypout,laiout,rfnsout,vegout,ndviout,          &
                x,y, rtem2,rtem3 )

!
!-----------------------------------------------------------------------
! Plotting out the results using NCAR graphics
!-----------------------------------------------------------------------
!
  CALL opngks

  tema = glat(1,ny)
  temb = glon(1,ny)
  temc = glat(nx,1)
  temd = glon(nx,1)

  PRINT *, 'The up-left corner (lat,lon) is ',tema,temb
  PRINT *, 'The down-right corner (lat,lon) is ',temc,temd

  DO is=1,nstyp

    WRITE (glab,'(a,i1,a)') 'Soil Type ',is,' in Model Grid'

    IF ( drawval == 0 ) THEN
      DO j = 1, ny
        DO i = 1, nx
          rtem2(i,j) = FLOAT( soiltyp(i,j,is) )
        END DO
      END DO

      cint = 1
      CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,               &
                 tema,temb,temc,temd,glab,rtem2  )
    ELSE
      CALL plotint( nx,ny,nx, mapproj,alatpro,alonpro,                  &
                    tema,temb,temc,temd,glab,soiltyp(1,1,is) )
    END IF

    cint = 0
    WRITE (glab,'(a,i1,a)')                                             &
          'Fraction of Soil Type ',is,' in Model Grid'
    CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,                 &
               tema,temb,temc,temd,glab,stypfrct(1,1,is)  )

  END DO

  WRITE (glab,'(a)') 'Vegetation Type in Model Grid'
  IF ( drawval == 0 ) THEN
    DO j = 1, ny
      DO i = 1, nx
        rtem2(i,j) = FLOAT( vegtyp(i,j) )
      END DO
    END DO

    cint = 1
    CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,                 &
               tema,temb,temc,temd,glab,rtem2  )
  ELSE
    CALL plotint( nx,ny,nx, mapproj,alatpro,alonpro,                    &
                  tema,temb,temc,temd,glab,vegtyp )
  END IF

  cint = 0
  WRITE (glab,'(a)') 'Leaf Area Index in Model Grid'
  CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,                   &
             tema,temb,temc,temd,glab,lai   )

  WRITE (glab,'(a)') 'Surface Roughness'
  CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,                   &
             tema,temb,temc,temd,glab,roufns )

  WRITE (glab,'(a)') 'Vegetation Fraction'
  CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,                   &
             tema,temb,temc,temd,glab,veg    )

  cint = 0.05

  WRITE (glab,'(a)') 'NDVI'

  if( ndviout /= 0 ) then

  CALL plot( nx,ny,nx, mapproj,alatpro,alonpro, cint,                   &
             tema,temb,temc,temd,glab,ndvi   )

  ENDIF

  CALL clsgks

  STOP
END PROGRAM arpssfc20
