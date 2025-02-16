!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE readcvttrn               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE readcvttrn(ternfile,ternfmt,nx,ny,dx,dy,                     &
                    mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,&
                    hterain)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the terrain data into model array hterain from a specified
!  terrain data file. This subroutine does not do grid checking as
!  readtrn does.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/20/2003 based on readtrn.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    ternfile   Terrain data file name
!    ternfmt    Terrain file format
!
!  OUTPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    dx       Grid interval in x-direction
!    dy       Grid interval in y-direction
!
!    mapproj    Type of map projection used to setup the analysis grid.
!    trulat1    1st real true latitude of map projection.
!    trulat2    2nd real true latitude of map projection.
!    trulon     Real true longitude of map projection.
!    sclfct     Map scale factor. At latitude = trulat1 and trulat2
!
!    ctrlat    Lat. at the origin of the model grid (deg. N)
!    ctrlon    Lon. at the origin of the model grid (deg. E)
!
!    hterain  Terrain height (m)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: ternfile    ! Terrain data file name
  INTEGER,          INTENT(IN)  :: ternfmt
  INTEGER,          INTENT(OUT) :: nx          ! Number of grid points in the x-direction
  INTEGER,          INTENT(OUT) :: ny          ! Number of grid points in the y-direction
  REAL,             INTENT(OUT) :: dx          ! Grid interval in x-direction
  REAL,             INTENT(OUT) :: dy          ! Grid interval in y-direction

  INTEGER,          INTENT(OUT) :: mapproj     ! Map projection
  REAL,             INTENT(OUT) :: trulat1     ! 1st real true latitude of map projection
  REAL,             INTENT(OUT) :: trulat2     ! 2nd real true latitude of map projection
  REAL,             INTENT(OUT) :: trulon      ! Real true longitude of map projection
  REAL,             INTENT(OUT) :: sclfct      ! Map scale factor
  REAL,             INTENT(OUT) :: ctrlat      ! Center latitude of the model domain (deg. N)
  REAL,             INTENT(OUT) :: ctrlon      ! Center longitude of the model domain (deg. E)

  REAL,             POINTER     :: hterain(:,:)! Terrain height.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  INTEGER :: inunit,istat
  INTEGER :: idummy,ierr
  REAL    :: rdummy,amin,amax

  INTEGER :: ireturn

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
         ! unused array in hdf routines since NO COMPRESSION

  INTEGER :: sd_id

  REAL, ALLOCATABLE :: var2d(:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Read in the terrain data.
!
!-----------------------------------------------------------------------
!
  WRITE (6,'(1x,/,3a,/)') "READCVTTRN: reading in external supplied ", &
                         "terrain data from file ",trim(ternfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (ternfmt == 1) THEN

    CALL getunit( inunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(ternfile, '-F f77 -N ieee', ierr)

    OPEN(UNIT=inunit,FILE=trim(ternfile),FORM='unformatted',            &
         STATUS='old',IOSTAT=istat)

    IF( istat /= 0) THEN
      WRITE(6,'(/1x,a,a,/1x,a/)')                                       &
          'Error occured when opening terrain data file ',              &
          ternfile,' Job stopped in READCVTTRN.'
      STOP
    END IF

    READ(inunit,ERR=999) nx,ny

    READ(inunit,ERR=999) idummy,mapproj,idummy,idummy,idummy,           &
               idummy,idummy,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy

    READ(inunit,ERR=999) dx ,dy    ,ctrlat,ctrlon,rdummy,               &
               rdummy,trulat1,trulat2,trulon,sclfct,                    &
               rdummy,rdummy,rdummy,rdummy,rdummy,                      &
               rdummy,rdummy,rdummy,rdummy,rdummy

  ELSE IF(ternfmt == 3) THEN

    CALL hdfopen(trim(ternfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "READTRN: ERROR opening ",                            &
                 trim(ternfile)," for reading."
      GO TO 999
    END IF

    CALL hdfrdi(sd_id,"nx",nx,istat)
    CALL hdfrdi(sd_id,"ny",ny,istat)
    CALL hdfrdr(sd_id,"dx",dx,istat)
    CALL hdfrdr(sd_id,"dy",dy,istat)
    CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
    CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
    CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
    CALL hdfrdr(sd_id,"trulon",trulon,istat)
    CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
    CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
    CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)
  
  ELSE IF (ternfmt == 7) THEN

    CALL netopen(TRIM(ternfile),'R',sd_id)
    CALL net_get_trn(sd_id,nx,ny,dx,dy,mapproj,sclfct,                  &
                  trulat1,trulat2,trulon,ctrlat,ctrlon, istat)

  ELSE 
   
    ! alternate dump format ...
    WRITE(6,*) 'The supported terrain data format are ',                &
               'binary (ternfmt=1) and HDF4 no compressed (ternfmt = 3).' 
    CALL arpsstop('Terrain data format is not supported.',1)

  END IF

  ALLOCATE(hterain(nx,ny), STAT = istat)

  IF (ternfmt == 1) THEN

    READ(inunit,ERR=999) hterain

    CALL retunit( inunit )
    CLOSE (UNIT=inunit)

  ELSE IF (ternfmt == 3) THEN

    CALL hdfrd2d(sd_id,"hterain",nx,ny,hterain,istat,itmp)
    CALL hdfclose(sd_id,istat)

  ELSE IF (ternfmt == 7) THEN

    ALLOCATE(var2d(nx-1,ny-1), STAT = istat)

    CALL netread2d(sd_id,0,0,'HTERAIN',nx-1,ny-1,var2d)
    CALL netclose(sd_id)

    DO j = 1,ny-1
      DO i = 1,nx-1
        hterain(i,j) = var2d(i,j)
      END DO
    END DO
    CALL edgfill(hterain,1,nx, 1, nx-1, 1,ny, 1,ny-1, 1,1,1,1)

    DEALLOCATE(var2d)

  END IF

  WRITE(6,'(1x,a/)') 'Minimum and maximum terrain height:'

  CALL a3dmax0(hterain,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1,amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'htrnmin = ', amin,', htrnmax=',amax

  RETURN

  999   WRITE(6,'(1x,a)')                                               &
        'Error in reading terrain data. Job stopped in READTRN.'
  CALL arpsstop("arpsstop called from READTRN reading file",1)

END SUBROUTINE readcvttrn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE READCVTSFC                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE readcvtsfc(sfcfile,sfcfmt,nx,ny,nstyps,dx,dy,                &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           stypin,vtypin,laiin,roufin,vegin,ndviin,                     &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the surface data sets from file sfcfile.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/20/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  sfcfile
!  sfcfmt
!
!  OUTPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
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
!  soiltyp  Soil type in model domain
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: sfcfile      ! Name of the surface data file
  INTEGER,          INTENT(IN)  :: sfcfmt
  INTEGER,          INTENT(OUT) :: nx           ! Number of grid points in the x-direction
  INTEGER,          INTENT(OUT) :: ny           ! Number of grid points in the y-direction
  INTEGER,          INTENT(OUT) :: nstyps       ! Max number of soil types in a grid box

  REAL,             INTENT(OUT) :: dx
  REAL,             INTENT(OUT) :: dy
  INTEGER,          INTENT(OUT) :: mapproj       ! Map projection
  REAL,             INTENT(OUT) :: trulat1       ! 1st real true latitude of map projection
  REAL,             INTENT(OUT) :: trulat2       ! 2nd real true latitude of map projection
  REAL,             INTENT(OUT) :: trulon        ! Real true longitude of map projection
  REAL,             INTENT(OUT) :: sclfct        ! Map scale factor
  REAL,             INTENT(OUT) :: ctrlat        ! Center latitude of the model domain (deg. N)
  REAL,             INTENT(OUT) :: ctrlon        ! Center longitude of the model domain (deg. E)

  INTEGER,          INTENT(OUT) :: stypin,vtypin,laiin,roufin,vegin,ndviin

  INTEGER,          POINTER :: soiltyp(:,:,:)  ! Soil type in model domain
  REAL,             POINTER :: stypfrct(:,:,:) ! Fraction of soil types
  INTEGER,          POINTER :: vegtyp (:,:)         ! Vegetation type in model domain

  REAL,             POINTER :: lai    (:,:)     ! Leaf Area Index in model domain
  REAL,             POINTER :: roufns (:,:)     ! NDVI in model domain
  REAL,             POINTER :: veg    (:,:)     ! Vegetation fraction
  REAL,             POINTER :: ndvi   (:,:)     ! NDVI
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: idummy
  REAL    :: rdummy

  INTEGER :: i,j,is
  INTEGER :: istat, ierr, ireturn
  INTEGER :: sfcunit

  INTEGER :: stat, sd_id

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL    :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER, ALLOCATABLE :: temi(:,:,:)
  REAL,    ALLOCATABLE :: temr(:,:,:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Open the surface data file. Read the parameters first, check if
!  the data are consistant with the model. If everything is OK, then
!  read the surface data, soiltyp, vegtyp, lai, roufns, veg, and
!  ndvi.
!
!-----------------------------------------------------------------------
!

  WRITE (6,'(1x,/,3a)') "READCVTSFC: reading in external supplied ",      &
                      "surface data from file ",trim(sfcfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (sfcfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    CALL getunit( sfcunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(sfcfile, '-F f77 -N ieee', ierr)

    OPEN(UNIT=sfcunit,FILE=trim(sfcfile),FORM='unformatted',            &
         STATUS='old',IOSTAT=istat)

    IF( istat /= 0 ) THEN

      WRITE(6,'(/1x,a,i2,/1x,a/)')                                      &
          'Error occured when opening the surface data file '           &
          //sfcfile//' using FORTRAN unit ',sfcunit,                    &
          ' Program stopped in READCVTSFC.'
      CALL arpsstop("arpsstop called from READCVTSFC opening file",1)

    END IF

    WRITE(6,'(/1x,a,/1x,a,i2/)')                                        &
        'This run will start from an external supplied surface ',       &
        'data file '//sfcfile//' using FORTRAN unit ',sfcunit

    READ (sfcunit,ERR=998) nx,ny

    READ (sfcunit,ERR=998) mapproj,stypin,vtypin,laiin,roufin,          &
                         vegin,  ndviin,nstyps,idummy,idummy,           &
                         idummy, idummy,idummy,idummy,idummy,           &
                         idummy, idummy,idummy,idummy,idummy

    READ (sfcunit,ERR=998) dx,dy, ctrlat,ctrlon,trulat1,                &
                         trulat2,trulon, sclfct,rdummy,rdummy,          &
                         rdummy,rdummy,rdummy,rdummy,rdummy,            &
                         rdummy,rdummy,rdummy,rdummy,rdummy

  ELSE IF (sfcfmt == 3) THEN 

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(sfcfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "READCVTSFC: ERROR opening ",                          &
                 trim(sfcfile)," for reading."
      GO TO 998
    END IF

    CALL hdfrdi(sd_id,"nstyp",nstyps,istat)

    CALL hdfrdi(sd_id,"nx",nx,istat)
    CALL hdfrdi(sd_id,"ny",ny,istat)
    CALL hdfrdr(sd_id,"dx",dx,istat)
    CALL hdfrdr(sd_id,"dy",dy,istat)
    CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
    CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
    CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
    CALL hdfrdr(sd_id,"trulon",trulon,istat)
    CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
    CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
    CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)

  ELSE IF (sfcfmt == 7) THEN

!-----------------------------------------------------------------------
!
!  NetCDF format.
!
!-----------------------------------------------------------------------

    CALL netopen(TRIM(sfcfile), 'R', sd_id)

    CALL net_get_sfc(sd_id,nx,ny,nstyps,dx,dy,                          &
                mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,    &
                stypin,vtypin,laiin,roufin,vegin,ndviin,istat)

  ELSE
   
    ! alternate dump format ...

    WRITE(6,'(1x,3a)') 'The supported surface data format are ',        &
               'binary (sfcfmt=1), HDF4 no compressed (sfcfmt = 3). ',  &
               'and NetCDF (sfcfmt = 7).' 
    CALL arpsstop('Surface data format is not supported.',1)

  END IF                     !sfcfmt loop

  nstyps = MAX( nstyps, 1 )

!-----------------------------------------------------------------------
!
!  Read in the surface data from the surface data file.
!
!-----------------------------------------------------------------------
  ALLOCATE(soiltyp (nx,ny,nstyps), STAT = istat)
  ALLOCATE(stypfrct(nx,ny,nstyps), STAT = istat)
  ALLOCATE(vegtyp  (nx,ny),        STAT = istat)
  ALLOCATE(lai     (nx,ny),        STAT = istat)
  ALLOCATE(roufns  (nx,ny),        STAT = istat)
  ALLOCATE(veg     (nx,ny),        STAT = istat)
  ALLOCATE(ndvi    (nx,ny),        STAT = istat)

  IF (sfcfmt == 1) THEN    ! Fortran unformatted

    WRITE (6, '(a/a,i2/a,i2/a,i2/a,i2/a,i2/a,i2)')                      &
      ' Surface data flags for: ',                                      &
      '        soiltyp --', stypin,                                     &
      '         vegtyp --', vtypin,                                     &
      '            lai --', laiin,                                      &
      '         roufns --', roufin,                                     &
      '            veg --', vegin,                                      &
      '           ndvi --', ndviin

    WRITE (6, '(a/a,i2)')                                               &
      ' Number of soil types in each grid box:',                        &
      '          nstyp --', nstyps

    IF(stypin == 1) THEN
      IF ( nstyps == 1 ) THEN
        WRITE (6, '(a)') 'Read in the soil type data'
        READ (sfcunit,ERR=998) ((soiltyp(i,j,1),i=1,nx),j=1,ny)
        DO j=1,ny
          DO i=1,nx
            stypfrct(i,j,1) = 1.0
          END DO
        END DO
      ELSE
        DO is=1,nstyps
            WRITE (6, '(a)') 'Read in the soil type data'
            READ (sfcunit,ERR=998) ((soiltyp(i,j,is),i=1,nx),j=1,ny)
            WRITE (6, '(a)') 'Read in the fraction of soil types'
            READ (sfcunit,ERR=998) ((stypfrct(i,j,is),i=1,nx),j=1,ny)
        END DO
      END IF
    END IF

    IF(vtypin == 1) THEN
      WRITE (6, '(a)') 'Read in the vegetation type data'
      READ (sfcunit,ERR=998) vegtyp
    END IF

    IF(laiin == 1) THEN
      WRITE (6, '(a)') 'Read in the Leaf Area Index data'
      READ (sfcunit,ERR=998) lai
    END IF

    IF(roufin == 1) THEN
      WRITE (6, '(a)') 'Read in the surface roughness data'
      READ (sfcunit,ERR=998) roufns
    END IF

    IF(vegin == 1) THEN
      WRITE (6, '(a)') 'Read in the vegetatin fraction data'
      READ (sfcunit,ERR=998) veg
    END IF

    IF (ndviin == 1) THEN
      WRITE (6, '(a)') 'Read in the NDVI data'
      READ (sfcunit,ERR=998) ndvi
    END IF

    CLOSE ( sfcunit )
    CALL retunit ( sfcunit )

  ELSE IF (sfcfmt == 3) THEN         

    CALL hdfrd3di(sd_id,"soiltyp",nx,ny,nstyps,soiltyp,stat)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in soiltyp'
      stypin = 1
    ELSE
      WRITE (6, '(a)') 'Variable soiltyp is not in the data set.'
      stypin = 0
    END IF
    CALL hdfrd3d(sd_id,"stypfrct",nx,ny,nstyps,                    &
                 stypfrct,stat,itmp,atmp1,atmp2)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in stypfrct'
    ELSE
      WRITE (6, '(a)') 'Variable stypfrct is not in the data set.'
      stypfrct(:,:,1) = 1.
      stypfrct(:,:,2:nstyps) = 0.
    END IF

    CALL hdfrd2di(sd_id,"vegtyp",nx,ny,vegtyp,stat)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in vegtyp'
      vtypin = 1
    ELSE
      WRITE (6, '(a)') 'Variable vegtyp is not in the data set.'
      vtypin = 0
    END IF

    CALL hdfrd2d(sd_id,"lai",nx,ny,lai,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in lai'
      laiin = 1
    ELSE
      WRITE (6, '(a)') 'Variable lai is not in the data set.'
      laiin = 0
    END IF

    CALL hdfrd2d(sd_id,"roufns",nx,ny,roufns,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in roufns'
      roufin = 1
    ELSE
      WRITE (6, '(a)') 'Variable roufns is not in the data set.'
      roufin = 0
    END IF

    CALL hdfrd2d(sd_id,"veg",nx,ny,veg,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in veg'
      vegin = 1
    ELSE
      WRITE (6, '(a)') 'Variable veg is not in the data set.'
      vegin = 0
    END IF

    CALL hdfrd2d(sd_id,"ndvi",nx,ny,ndvi,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in ndvi'
      ndviin = 1
    ELSE
      WRITE (6, '(a)') 'Variable ndvi is not in the data set.'
      ndviin = 0
    END IF

    CALL hdfclose(sd_id, stat)

  ELSE IF (sfcfmt == 7) THEN    ! NetCDF format

    WRITE (6, '(a/a,i2/a,i2/a,i2/a,i2/a,i2/a,i2)')                      &
                            ' Surface data flags for: ',                &
                            '        soiltyp --', stypin,               &
                            '         vegtyp --', vtypin,               &
                            '            lai --', laiin,                &
                            '         roufns --', roufin,               &
                            '            veg --', vegin,                &
                            '           ndvi --', ndviin

    WRITE (6, '(a/a,i2)') ' Number of soil types in each grid box:',    &
                          '          nstyp --', nstyps

    ALLOCATE(temr(nx-1,ny-1,nstyps), STAT = stat)
    ALLOCATE(temi(nx-1,ny-1,nstyps), STAT = stat)

    IF(stypin == 1) THEN
      CALL netread3di(sd_id,0,0,'SOILTYP',nx-1,ny-1,nstyps,temi)
      CALL netread3d (sd_id,0,0,'STYPFRCT',nx-1,ny-1,nstyps,temr)
      DO is = 1, nstyps
        DO j = 1, ny-1
          DO i = 1, nx-1
            soiltyp (i,j,is) = temi(i,j,is)
            stypfrct(i,j,is) = temr(i,j,is)
          END DO
        END DO
      END DO
      CALL iedgfill(soiltyp, 1,nx,1,nx-1,1,ny,1,ny-1,1,nstyps,1,nstyps)
      CALL  edgfill(stypfrct,1,nx,1,nx-1,1,ny,1,ny-1,1,nstyps,1,nstyps)
    END IF

    IF(vtypin == 1) THEN
      CALL netread2di(sd_id,0,0,'VEGTYP',nx-1,ny-1,temi)
      DO j = 1, ny-1
        DO i = 1, nx-1
          vegtyp (i,j) = temi(i,j,1)
        END DO
      END DO
      CALL iedgfill(vegtyp, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF(laiin == 1) THEN
      CALL netread2d(sd_id,0,0,'LAI',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          lai (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(lai,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF(roufin == 1) THEN
      CALL netread2d(sd_id,0,0,'ROUFNS',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          roufns (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(roufns,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF(vegin == 1) THEN
      CALL netread2d(sd_id,0,0,'VEG',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          veg (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(veg,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF (ndviin == 1) THEN
      CALL netread2d(sd_id,0,0,'NDVI',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          ndvi (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(ndvi,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    CALL netclose(sd_id)

    DEALLOCATE(temi, temr)

  END IF


  RETURN

  998   WRITE (6,'(/a,i2/a)')                                           &
         'READCVTSFC: Read error in surface data file '                 &
         //sfcfile//' with the I/O unit ',sfcunit,                      &
         'The model will STOP in subroutine READCVTSFC.'

  CALL arpsstop("arpsstop called from READCVTSFC reading sfc file",1)

END SUBROUTINE readcvtsfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE READCVTSOIL               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE readcvtsoil(soilfile,soilfmt,nx,ny,nzsoil,nstyps,dx,dy,      &
           zpsoil,mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,  &
           zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,                    &
           tsoil,qsoil,wetcanp,snowdpth,soiltyp )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the soil variables and ARPS grid parameters from file soilfile.
!  It does the same job as subroutine readsoil in src/arps/iolib3d.f90. 
!  However, it skip the step to check the consistence between the passed
!  in grid parameters with those read in from the file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  Changed from readsoil in src/arps/iolib3d.f90.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!  
!  soilfile
!  soilfmt
!
!  OUTPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
!  nzsoil   Number of model grid points in the soil.  
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
!  tsoil    Soil temperature (K)
!  qsoil    Soil moisture (m3/m3) 
!  wetcanp  Canopy water amount
!  snowdpth Snow depth (m)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: soilfile ! Name of the soil file
  INTEGER,          INTENT(IN) :: soilfmt

  INTEGER, INTENT(OUT) :: nx            ! Number of grid points in the x-direction
  INTEGER, INTENT(OUT) :: ny            ! Number of grid points in the y-direction
  INTEGER, INTENT(OUT) :: nzsoil        ! Number of grid points in the soil.  
  INTEGER, INTENT(OUT) :: nstyps        ! Number of soil types for each grid point

  REAL,    INTENT(OUT) :: dx
  REAL,    INTENT(OUT) :: dy
  INTEGER, INTENT(OUT) :: mapproj       ! Map projection
  REAL,    INTENT(OUT) :: trulat1       ! 1st real true latitude of map projection
  REAL,    INTENT(OUT) :: trulat2       ! 2nd real true latitude of map projection
  REAL,    INTENT(OUT) :: trulon        ! Real true longitude of map projection
  REAL,    INTENT(OUT) :: sclfct        ! Map scale factor
  REAL,    INTENT(OUT) :: ctrlat        ! Center latitude of the model domain (deg. N)
  REAL,    INTENT(OUT) :: ctrlon        ! Center longitude of the model domain (deg. E)

  INTEGER, INTENT(OUT) :: zpsoilin
  INTEGER, INTENT(OUT) :: tsoilin
  INTEGER, INTENT(OUT) :: qsoilin  
  INTEGER, INTENT(OUT) :: wcanpin
  INTEGER, INTENT(OUT) :: snowdin

  REAL,    POINTER     :: zpsoil (:,:,:)   ! Soil depths (m) 
  REAL,    POINTER     :: tsoil  (:,:,:,:) ! Soil temperature (K)
  REAL,    POINTER     :: qsoil  (:,:,:,:) ! Soil moisture (m3/m3) 
  REAL,    POINTER     :: wetcanp(:,:,:)   ! Canopy water amount
  INTEGER, POINTER     :: soiltyp(:,:,:)   ! Soil type in model domain
  REAL,    POINTER     :: snowdpth(:,:)    ! Snow depth (m)

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: flunit
  INTEGER :: idummy
  REAL    :: rdummy

  INTEGER :: i,j,k,is
  INTEGER :: istat, ierr

  INTEGER :: ireturn

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL :: atmp1(1),atmp2(1)
  ! unused arrays because of no-compression

  INTEGER :: stat, sd_id
  !
  ! fmtver??: to label each data a version.
  ! intver??: an integer to allow faster comparison than fmtver??,
  !           which are strings.
  !
  ! Verion 5.00: significant change in soil variables since version 4.10.
  !
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='* 004.10 GrADS Soilvar Data',intver410=410)
  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin

  INTEGER :: tsfcin        ! for backward compatibility Zuwen He, 07/01/02 
  INTEGER :: wsfcin        ! for backward compatibility Zuwen He, 07/01/02 
  INTEGER :: wdpin         ! for backward compatibility Zuwen He, 07/01/02 

  INTEGER :: snowcin
  INTEGER :: stypin

!-----------------------------------------------------------------------
!
! Working arrays
!
!-----------------------------------------------------------------------

  REAL,    ALLOCATABLE :: tem1(:,:,:) ! Temporary array

  REAL,    ALLOCATABLE :: var3d (:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)        
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  WRITE (6,'(1x,/,3a)') "READCVTSOIL: reading in external supplied ",   &
                        "soil data from file ",trim(soilfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (soilfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    CALL getunit( flunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(soilfile, '-F f77 -N ieee', ierr)

    OPEN (UNIT=flunit,FILE=trim(soilfile),FORM='unformatted',           &
          STATUS='old',IOSTAT=istat)

    IF( istat /= 0 ) THEN
      WRITE(6,'(/1x,a,i2,/1x,a/)')                                      &
          'Error occured when opening the soil data file '              &
          //soilfile//' using FORTRAN unit ',flunit,                    &
          ' Program stopped in READCVTSOIL.'
      CALL arpsstop("arpsstop called from READCVTSOIL while opening file",1)
    END IF

    WRITE(6,'(/1x,a,/1x,a,i2/)')                                        &
        'This run will start from an external supplied soil ',          &
        'data file '//soilfile//' using FORTRAN unit ',flunit

    READ (flunit,ERR=997) fmtverin
    
    ! The following code is not a safe practice. 
    !
    ! However, this may be the only way to distinguish versions 
    ! prior to 500. 
    
    IF (fmtverin == fmtver500) THEN
      intver=intver500
    ELSE
      WRITE(6,'(/1x,a/)')   & 
          'WARNING: Incoming data format are older than version 5.00!!! '
    END IF

    WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

    READ (flunit,ERR=998) nx,ny,nzsoil

    GOTO 996 

    997 WRITE(6,'(/1x,a,a/)')                        &
      'Incoming data format: fmtver=fmtver410. Data read-in may be wrong.'

    CLOSE (flunit) 
    OPEN (UNIT=flunit,FILE=trim(soilfile),FORM='unformatted',           &
          STATUS='old',IOSTAT=istat)

    READ (flunit,ERR=998) nx,ny
    nzsoil=2
    intver=intver410   ! there is no fmtverin prior to version 500
    fmtver=fmtver410 
    WRITE(6,'(/1x,a/,a/)')   & 
          'WARNING: Incoming data format are to read as version 4.10'

    996 CONTINUE 

    IF (intver == intver410) THEN 

      READ (flunit,ERR=998) mapproj,tsfcin,tsoilin,wsfcin,wdpin,        &
                        wcanpin,snowcin,snowdin,stypin,zpsoilin,        &
                        idummy, idummy, idummy, idummy,idummy,          &
                        idummy, idummy, idummy, idummy,nstyps

    ELSE IF (intver >= intver500) THEN 

      READ (flunit,ERR=998) mapproj,tsoilin,qsoilin,                    &
                        wcanpin,snowcin,snowdin,stypin,zpsoilin,        &
                        idummy, idummy, idummy, idummy,idummy,          &
                        idummy, idummy, idummy, idummy,nstyps

    END IF 

    READ (flunit,ERR=998) dx,   dy,    ctrlat,ctrlon,trulat1,          &
                         trulat2,trulon,sclfct,rdummy,rdummy,           &
                         rdummy,rdummy,rdummy,rdummy,rdummy,           &
                         rdummy,rdummy,rdummy,rdummy,rdummy

  ELSE IF (soilfmt == 3) THEN     !HDF4 

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(soilfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "READCVTSOIL: ERROR opening ",                           &
                 trim(soilfile)," for reading."
      GO TO 998
    END IF

    CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

!
! The following code is a dangerous practice 
! but it may be the only way to distinguish 
! versions prior to 500. 
!
    IF (fmtverin == fmtver500) THEN
      intver=intver500
    ELSE 
      intver=intver410  ! prior to 500, there is no fmtver variable
      istat=0
      WRITE(6,'(/1x,a/,a/)')   & 
          'WARNING: Incoming data format are older than version 5.00!!! ', & 
          'It is to be read as if it version 4.10!!! '
    END IF 

    CALL hdfrdi(sd_id,"nstyp",nstyps,istat)
    CALL hdfrdi(sd_id,"nx",nx,istat)
    CALL hdfrdi(sd_id,"ny",ny,istat)
    
    IF (intver >= intver500) THEN 
      CALL hdfrdi(sd_id,"nzsoil",nzsoil,istat)
    ELSE 
      nzsoil = 2  ! prior to version 500, it is 2 layer soil
    END IF 

    CALL hdfrdr(sd_id,"dx",dx,istat)
    CALL hdfrdr(sd_id,"dy",dy,istat)
    CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
    CALL hdfrdr(sd_id,"trulat1",trulat1, istat)
    CALL hdfrdr(sd_id,"trulat2",trulat2, istat)
    CALL hdfrdr(sd_id,"trulon", trulon,  istat)
    CALL hdfrdr(sd_id,"sclfct", sclfct, istat)
    CALL hdfrdr(sd_id,"ctrlat", ctrlat, istat)
    CALL hdfrdr(sd_id,"ctrlon", ctrlon, istat)

  ELSE IF (soilfmt == 7) THEN

    CALL netopen(TRIM(soilfile), 'R', sd_id)
    CALL net_get_soil(sd_id,nx,ny,nzsoil,nstyps,dx,dy,                  &
                mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,    &
                zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,stypin,istat)

  ELSE
   
    ! alternate dump format ...
    WRITE(6,'(1x,3a)') 'The supported soil data format are ',           &
               'binary (soilfmt=1), HDF4 no compressed (soilfmt = 3).', &
               'and NetCDF (soilfmt = 7).' 
    CALL arpsstop('Soil data format is not supported.',1)

  END IF

  nstyps = MAX(nstyps, 1)

  ALLOCATE (zpsoil(nx,ny,nzsoil),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "READCVTSOIL: ERROR allocating zpsoil, returning"
    RETURN
  END IF

  ALLOCATE (tsoil(nx,ny,nzsoil,0:nstyps),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "READCVTSOIL: ERROR allocating tsoil, returning"
    RETURN
  END IF
  ALLOCATE (qsoil(nx,ny,nzsoil,0:nstyps),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "READCVTSOIL: ERROR allocating qsoil, returning"
    RETURN
  END IF
  ALLOCATE (wetcanp(nx,ny,0:nstyps),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "READCVTSOIL: ERROR allocating wetcanp, returning"
    RETURN
  END IF

  ALLOCATE (soiltyp(nx,ny,nstyps),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "READCVTSOIL: ERROR allocating soiltyp, returning"
    RETURN
  END IF

  ALLOCATE (snowdpth(nx,ny),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "READCVTSOIL: ERROR allocating snowdpth, returning"
    RETURN
  END IF

  WRITE (6,'(a//a,i1/6(a,f7.2/))')                                     &
      ' The map projection and griding information for the soil data:', &
      ' Projection:                 ', mapproj,                         &
      ' The 1st real true latitude: ', trulat1,                         &
      ' The 2nd real true latitude: ', trulat2,                         &
      ' The real true longitude:    ', trulon,                          &
      ' Map scale factor:           ', sclfct,                          &
      ' Latitude  at the origin:    ', ctrlat,                          &
      ' Longitude at the origin:    ', ctrlon

!
!-----------------------------------------------------------------------
!
!  Read in the soil data from the soil data file.
!
!-----------------------------------------------------------------------
!
  IF (intver == intver410) THEN 
    ALLOCATE(tem1(nx,ny,0:nstyps),stat=istat)  ! for reading old version
                                                ! tsfc,tsoil,wetsfc,wetdp
  END IF 

  IF (soilfmt == 1) THEN    ! Fortran unformatted

    IF (intver == intver410) THEN 

      WRITE (6, '(a/8(a,i2/))')  ' Surface data flags for: ',            &
                                 '        zpsoilin --', zpsoilin,        &
                                 '        tsfcin   --', tsfcin,          &
                                 '        tsoilin  --', tsoilin,         &
                                 '        wsfcin   --', wsfcin,          &
                                 '        wdpin    --', wdpin,           &
                                 '        wcanpin  --', wcanpin,         &
                                 '        snowdin  --', snowdin,         &
                                 '        stypin   --', stypin

    ELSE IF (intver == intver500) THEN 

      WRITE (6, '(a/6(a,i2/))') ' Surface data flags for: ',            &
                                '        zpsoilin --', zpsoilin,        &
                                '        tsoilin  --', tsoilin,         &
                                '        qsoilin  --', qsoilin,         &
                                '        wcanpin  --', wcanpin,         &
                                '        snowdin  --', snowdin,         &
                                '        stypin   --', stypin
    END IF 

    IF (intver == intver410) THEN 
      zpsoil(:,:,1)=0.
      zpsoil(:,:,2)=-1. 
    ELSE IF (intver >= intver500) THEN 
      IF (zpsoilin /= 0) THEN
        WRITE(6,'(a)') 'Read in the soil depth data'
        DO k=1,nzsoil
          READ (flunit,ERR=998) ((zpsoil(i,j,k),i=1,nx),j=1,ny) 
        END DO 
      END IF
    END IF  ! intver

    IF (intver == intver410) THEN 
      IF ( tsfcin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the surface skin temperature data'
        IF ( nstyps == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nx),j=1,ny)
          tsoil(:,:,1,0)=tem1(:,:,0) 
        ELSE
          READ (flunit,ERR=998) tem1
          tsoil(:,:,1,:)=tem1(:,:,:) 
        END IF
      ELSE
        WRITE(6,'(a)') 'Variable tsfc is not in the data set.'
      END IF

      IF ( tsoilin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the deep soil temperature data'
        IF ( nstyps == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nx),j=1,ny)
          tsoil(:,:,2,0)=tem1(:,:,0) 
        ELSE
          READ (flunit,ERR=998) tem1
          tsoil(:,:,2,:)=tem1(:,:,:) 
        END IF
      ELSE
        WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
      END IF

      IF ( wsfcin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the skin soil moisture data'
        IF ( nstyps == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nx),j=1,ny)
          qsoil(:,:,1,0)=tem1(:,:,0) 
        ELSE
          READ (flunit,ERR=998) tem1
          qsoil(:,:,1,:)=tem1(:,:,:) 
        END IF
      ELSE
      WRITE(6,'(a)') 'Variable wetsfc is not in the data set.'
      END IF

      IF ( wdpin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the deep soil moisture data'
        IF ( nstyps == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nx),j=1,ny)
          qsoil(:,:,2,0)=tem1(:,:,0) 
        ELSE
          READ (flunit,ERR=998) tem1
          qsoil(:,:,2,:)=tem1(:,:,:) 
        END IF
      ELSE
        WRITE(6,'(a)') 'Variable wetdp is not in the data set.'
      END IF

      IF ( wcanpin /= 0 ) THEN
        IF ( nstyps == 1 ) THEN
          READ (flunit,ERR=998) ((wetcanp(i,j,0),i=1,nx),j=1,ny)
        ELSE
          READ (flunit,ERR=998) wetcanp
        END IF
      ELSE
        WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
      END IF

      IF ( snowcin /= 0 ) THEN
        WRITE (6, '(a)') 'File contains snowcvr -- discarding'
        READ (flunit,ERR=998)
      END IF
  
      IF ( snowdin /= 0 ) THEN
        WRITE (6, '(a)') 'Read in the snow depth data'
        READ (flunit,ERR=998) snowdpth
      ELSE
        WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
      END IF

      IF ( stypin /= 0 ) THEN
        WRITE (6, '(a)') 'Read soil type of soil data.'
        READ (flunit,ERR=998) soiltyp
      END IF

    ELSE IF (intver >= intver500) THEN 

      IF ( tsoilin /= 0 ) THEN
        DO is=0,nstyps
          WRITE(6,'(a,i4)') 'Read in the soil temperature for soil type ',is
          DO k=1,nzsoil
            READ (flunit,ERR=998) ((tsoil(i,j,k,is),i=1,nx),j=1,ny)
          END DO 
        END DO 
      ELSE
        WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
      END IF

      IF ( qsoilin /= 0 ) THEN
        DO is=0,nstyps
          WRITE(6,'(a,i4)') 'Read in the soil moisture data for soil type ',is
          DO k=1,nzsoil
            READ (flunit,ERR=998) ((qsoil(i,j,k,is),i=1,nx),j=1,ny)
          END DO 
        END DO 
      ELSE
        WRITE(6,'(a)') 'Variable qsoil is not in the data set.'
      END IF

      IF ( wcanpin /= 0 ) THEN
        DO is=0,nstyps
          WRITE (6, '(a,i4)') 'Read in the canopy water amount data for soil type ',is
          READ (flunit,ERR=998) ((wetcanp(i,j,is),i=1,nx),j=1,ny)
        END DO 
      ELSE
        WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
      END IF

      IF ( snowcin /= 0 ) THEN
        WRITE (6, '(a)') 'File contains snowcvr -- discarding'
        READ (flunit,ERR=998)
      END IF
  
      IF ( snowdin /= 0 ) THEN
        WRITE (6, '(a)') 'Read in the snow depth data'
        READ (flunit,ERR=998) ((snowdpth(i,j),i=1,nx),j=1,ny)
      ELSE
        WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
      END IF
  
      IF ( stypin /= 0 ) THEN
        DO is=1,nstyps
          WRITE (6, '(a,i4)') 'Read soil type of soil data for soil type ',is
          READ (flunit,ERR=998) ((soiltyp(i,j,is),i=1,nx),j=1,ny)
        END DO 
      END IF

    END IF 

    CLOSE ( flunit )
    CALL retunit ( flunit )

  ELSE IF (soilfmt == 3) THEN    

    IF (intver <= intver410) THEN 

      WRITE(6,'(a)') 'WARNING: No zpsoil is defined in this version. '
      WRITE(6,'(a)') 'Assume zpsoil(,,1)=0 and zpsoil(,,2)=-1.'
      zpsoil(:,:,1)=0.
      zpsoil(:,:,2)=-1.

      CALL hdfrd3d(sd_id,"tsfc",nx,ny,nstyps+1,tem1,stat,                 &
                   itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the surface skin temperature data'
        tsfcin = 1
      ELSE
        WRITE(6,'(a)') 'Variable tsfc is not in the data set.'
        tsfcin = 0
      END IF
      tsoil(:,:,1,:)=tem1(:,:,:) 
  
      CALL hdfrd3d(sd_id,"tsoil",nx,ny,nstyps+1,tem1,stat,                 &
                   itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the deep soil temperature data'
        tsoilin = 1
      ELSE
        WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
        tsoilin = 0
      END IF
      tsoil(:,:,2,:)=tem1(:,:,:) 

      CALL hdfrd3d(sd_id,"wetsfc",nx,ny,nstyps+1,tem1,stat,                &
                   itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the skin soil moisture data'
        wsfcin = 1
      ELSE
        WRITE(6,'(a)') 'Variable wetsfc is not in the data set.'
        wsfcin = 0
      END IF
      qsoil(:,:,1,:)=tem1(:,:,:) 

      CALL hdfrd3d(sd_id,"wetdp",nx,ny,nstyps+1,tem1,stat,                 &
                   itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the deep soil moisture data'
        wdpin = 1
      ELSE
        WRITE(6,'(a)') 'Variable wetdp is not in the data set.'
        wdpin = 0
      END IF
      qsoil(:,:,2,:)=tem1(:,:,:) 

    ELSE IF (intver >= intver500) THEN 

      CALL hdfrd3d(sd_id,"zpsoil",nx,ny,nzsoil,zpsoil,stat,         &
                 itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the soil layer depth data'
        zpsoilin = 1
      ELSE
        WRITE(6,'(a)') 'Variable zpsoil is not in the data set.'
        zpsoilin = 0
      END IF

      CALL hdfrd4d(sd_id,"tsoil",nx,ny,nzsoil,nstyps+1,tsoil,stat,   &
                   itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the soil temperature data'
        tsoilin = 1
      ELSE
        WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
        tsoilin = 0
      END IF
  
      CALL hdfrd4d(sd_id,"qsoil",nx,ny,nzsoil,nstyps+1,qsoil,stat,   &
                   itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE(6,'(a)') 'Read in the soil moisture data'
        qsoilin = 1
      ELSE
        WRITE(6,'(a)') 'Variable qsoil is not in the data set.'
        qsoilin = 0
      END IF

    END IF 

    CALL hdfrd3d(sd_id,"wetcanp",nx,ny,nstyps+1,wetcanp,stat,         &
                 itmp,atmp1,atmp2)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in the canopy water amount data'
      wcanpin = 1
    ELSE
      WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
      wcanpin = 0
    END IF

    CALL hdfrd2d(sd_id,"snowdpth",nx,ny,snowdpth,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read in the snow depth data'
      snowdin = 1
    ELSE
      WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
      snowdin = 0
    END IF

    CALL hdfrd3di(sd_id,"soiltyp",nx,ny,nstyps,soiltyp,stat)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(a)') 'Read soil type of soil data'
      stypin = 1
    ELSE
      WRITE (6, '(a)') 'Soil type of soil data is not in the data set.'
      stypin = 0
    END IF

    CALL hdfclose(sd_id, stat)

  ELSE IF (soilfmt == 7) THEN

    ALLOCATE(var3d (nx-1,ny-1,MAX(nzsoil,nstyps+1)), STAT = istat)
    ALLOCATE(var3di(nx-1,ny-1,nstyps),               STAT = istat)
    ALLOCATE(var4d (nx-1,ny-1,nzsoil,nstyps+1),      STAT = istat)

    IF (zpsoilin == 1) THEN
      CALL netread3d(sd_id,0,0,"ZPSOIL",nx-1,ny-1,nzsoil,var3d)
      WRITE(6,'(1x,a)') 'Read in the soil layer depth data.'
      DO k = 1, nzsoil
        DO j = 1,ny-1
          DO i = 1,nx-1
            zpsoil(i,j,k) = var3d(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(zpsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil)
    END IF

    IF (tsoilin == 1) THEN
      CALL netread4d(sd_id,0,0,"TSOIL",nx-1,ny-1,nzsoil,nstyps+1,var4d)
      WRITE(6,'(1x,a)') 'Read in the soil temperature data.'
      DO is = 0, nstyps
        DO k = 1, nzsoil
          DO j = 1,ny-1
            DO i = 1,nx-1
              tsoil(i,j,k,is) = var4d(i,j,k,is+1)
            END DO
          END DO
        END DO
        CALL edgfill(tsoil(:,:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil)
      END DO

    END IF

    IF (qsoilin == 1) THEN
      CALL netread4d(sd_id,0,0,"QSOIL",nx-1,ny-1,nzsoil,nstyps+1,var4d)
      WRITE(6,'(1x,a)') 'Read in the soil moisture data.'
      DO is = 0, nstyps
        DO k = 1, nzsoil
          DO j = 1,ny-1
            DO i = 1,nx-1
              qsoil(i,j,k,is) = var4d(i,j,k,is+1)
            END DO
          END DO
        END DO
        CALL edgfill(qsoil(:,:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil)
      END DO
    END IF

    IF (wcanpin == 1) THEN
      CALL netread3d(sd_id,0,0,"WETCANP",nx-1,ny-1,nstyps+1,var3d)
      WRITE (6, '(1x,a)') 'Read in the canopy water amount data.'
      DO is = 0, nstyps
        DO j = 1,ny-1
          DO i = 1,nx-1
            wetcanp(i,j,is) = var3d(i,j,is+1)
          END DO
        END DO
      END DO
      CALL edgfill(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,nstyps+1,1,nstyps+1)
    END IF

    IF (snowdin == 1) THEN
      CALL netread2d(sd_id,0,0,"SNOWDPTH",nx-1,ny-1,var3d)
      WRITE (6, '(1x,a)') 'Read in the snow depth data.'
      DO j = 1,ny-1
        DO i = 1,nx-1
          snowdpth(i,j) = var3d(i,j,1)
        END DO
      END DO
      CALL edgfill(snowdpth,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF (stypin == 1) THEN
      CALL netread3di(sd_id,0,0,"SOILTYP",nx-1,ny-1,nstyps,var3di)
      WRITE (6, '(1x,a)') 'Read in soil type of soil data.'
      DO is = 1, nstyps
        DO j = 1,ny-1
          DO i = 1,nx-1
            soiltyp(i,j,is) = var3di(i,j,is)
          END DO
        END DO
      END DO
      CALL iedgfill(soiltyp,1,nx,1,nx-1,1,ny,1,ny-1,1,nstyps,1,nstyps)
    END IF

    CALL netclose(sd_id)

    DEALLOCATE(var3d,var3di)
    DEALLOCATE(var4d)

  END IF

  IF (intver == intver410) THEN 

    IF (tsfcin /= tsoilin .OR. wsfcin /= wdpin) THEN 
      WRITE (6,'(1x,a,a/,a/,a/)') 'READCVTSOIL: WARNING: ',             &
                'The soilvar data is of version ', fmtver410,           & 
                '. The inconsistency flag between tsfcin and tsoilin, ',& 
                ' or between wsfin and wdpin, may cause some problems. '   
    END IF 
    tsoilin = max(tsfcin,tsoilin) 
    qsoilin = max(wsfcin,wdpin) 

    DEALLOCATE(tem1)
  END IF 

  RETURN

  998   WRITE (6,'(/a,i2/a)') '     Read error in soil data file '      &
                    //soilfile//' with the I/O unit ',flunit,           &
                    'The model will STOP in subroutine READCVTSOIL.'

  CALL arpsstop("arpsstop called from READCVTSOIL reading surface data",1)

  RETURN
END SUBROUTINE readcvtsoil
