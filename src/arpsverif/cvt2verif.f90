!TODO: List filenames in the input file rather than having the program
!      assemble the filenames...EMK
!TODO: nlevel, nvar input rather than hardwired

!########################################################################
!######################################################################## 
!######                                                            ######
!######                    PROGRAM CVT2VERIF                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM cvt2verif

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Read in an ARPS history dump or RUC or Eta GRIB file, and:
!   1.  Calculates user specified variables at user specified
!       levels.
!   2.  Vertically interpolates to user specified levels, if
!       necessary.
!   3.  Outputs data to an HDF file for use by VERIFGRID.
!
!----------------------------------------------------------------------- 
!
! HISTORY:
!   1999/11/01	First written by Eric Kemp
!   1999/12/29	Modified for multidimensional arrays by Richard Carpenter
!   2000/03/31  Major changes by Eric Kemp.  Modified the HDF
!               output to include the four corners of the grid (lat/lon 
!               coordinates) for the NCL map plotting; replaced some
!               modules with FORTRAN 77 include files; passed 'model'
!               and 'grid' data from input file to program.
!
!   2000/07/07  Added checks for non-blank characters in model and grid
!               variables (HDF does not like writing blank character
!               variables).
!   2002/03/29  Jason Levit (Williams SRA)
!               Modifications to the "getarps" call to update
!               for use with ARPS 5.0
!
!-----------------------------------------------------------------------
!
! Use modules
!
!-----------------------------------------------------------------------

  USE extdims2
  USE verif

!-----------------------------------------------------------------------
!
! Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE                    ! Force explicit declarations

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'

!-----------------------------------------------------------------------
!
! 2-D verification variables on forecast grid (x,y)
!
!-----------------------------------------------------------------------
  
  REAL, ALLOCATABLE :: var_p(:,:,:,:,:)	! dims: nx, ny, nZ, nT, nVar
  REAL, ALLOCATABLE :: var_sfc(:,:,:,:)	! dims: nx, ny, nT, nVar

!-----------------------------------------------------------------------
!
! External grid variables
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz,nxinput,nyinput,nzinput
                                            
  INTEGER :: nstyps
      
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tem1, tem2, tem3, tem4

  REAL, ALLOCATABLE :: vtem2d(:,:)

  INTEGER :: iproj_ext      ! external data map projection
  REAL :: scale_ext         ! external data map scale factor
  REAL :: trlon_ext         ! external data true longitude
  REAL :: latnot_ext(2)     ! external data true latitude(s)
  REAL :: x0_ext,y0_ext     ! external data origin
  REAL,ALLOCATABLE :: x_ext(:)     ! external data x-coordinate
  REAL,ALLOCATABLE :: y_ext(:)     ! external data y-coordinate
  REAL,ALLOCATABLE :: lat_ext(:,:)    ! external data latidude  
  REAL,ALLOCATABLE :: lon_ext(:,:)    ! external data longitude
  REAL,ALLOCATABLE :: latu_ext(:,:)
  REAL,ALLOCATABLE :: lonu_ext(:,:)
  REAL,ALLOCATABLE :: latv_ext(:,:)
  REAL,ALLOCATABLE :: lonv_ext(:,:)

!-----------------------------------------------------------------------
!
! External forecast variables
!
!-----------------------------------------------------------------------

  REAL,ALLOCATABLE :: p_ext(:,:,:)     ! Pressure (Pascals)
  REAL,ALLOCATABLE :: hgt_ext(:,:,:)   ! Height (m)
  REAL,ALLOCATABLE :: zp_ext(:,:,:)    ! Height (m) (on arps grid)
  REAL,ALLOCATABLE :: t_ext(:,:,:)     ! Temperature (K)
  REAL,ALLOCATABLE :: u_ext(:,:,:)     ! Eastward wind component
  REAL,ALLOCATABLE :: v_ext(:,:,:)     ! Northward wind component
  REAL,ALLOCATABLE :: w_ext(:,:,:)     ! Vertical velocity component
  REAL,ALLOCATABLE :: uatv_ext(:,:,:)  
  REAL,ALLOCATABLE :: vatu_ext(:,:,:) 
  
  REAL,ALLOCATABLE :: qv_ext(:,:,:)    ! Specific humidity (kg/kg)
  REAL,ALLOCATABLE :: qc_ext(:,:,:)    ! Cloud H2O mixing ratio
                                           !(kg/kg)
  REAL,ALLOCATABLE :: qr_ext(:,:,:)    ! Rain H2O mixing ratio
                                           !(kg/kg)
  REAL,ALLOCATABLE :: qi_ext(:,:,:)    ! Ice H2O mixing ratio
                                           !(kg/kg)
  REAL,ALLOCATABLE :: qs_ext(:,:,:)    ! Snow H2O mixing ratio
                                           !(kg/kg)
  REAL,ALLOCATABLE :: qh_ext(:,:,:)    ! Hail H2O mixing ratio
                                           !(kg/kg)
  INTEGER,ALLOCATABLE :: snowcvr_ext(:,:)      ! Snow cover

  INTEGER, ALLOCATABLE :: soiltyp_ext (:,:,:)    ! Soil type
  REAL, ALLOCATABLE ::    stypfrct_ext(:,:,:)    ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp_ext  (:,:)      ! Vegetation type
  REAL, ALLOCATABLE ::  lai_ext     (:,:)        ! Leaf Area Index
  REAL, ALLOCATABLE ::  roufns_ext  (:,:)        ! Surface roughness
  REAL, ALLOCATABLE ::  veg_ext     (:,:)        ! Vegetation fraction

  REAL, ALLOCATABLE  :: tsoil_ext  (:,:,:)  ! Deep soil temperature (K)
  REAL, ALLOCATABLE  :: qsoil_ext  (:,:,:)  ! soil moisture 
  REAL, ALLOCATABLE  :: wetcanp_ext(:,:,:)  ! Canopy water amount
  REAL, ALLOCATABLE  :: snowdpth_ext(:,:)   ! Snow depth (m)
      
  REAL,ALLOCATABLE :: trn_ext (:,:)      ! External terrain (m)
  REAL,ALLOCATABLE :: psfc_ext (:,:)      ! Surface pressure (Pa)

  REAL*4,ALLOCATABLE :: T_2m_ext (:,:)
  REAL*4,ALLOCATABLE :: RH_2m_ext(:,:)
  REAL*4,ALLOCATABLE :: U_10m_ext(:,:)
  REAL*4,ALLOCATABLE :: V_10m_ext(:,:)
  REAL,ALLOCATABLE :: MSLP_ext(:,:)
  REAL,ALLOCATABLE :: vMSLP_ext(:,:)
  REAL,ALLOCATABLE :: RH_ext(:,:,:)
  REAL,ALLOCATABLE :: vRH_ext(:,:,:)

  REAL,ALLOCATABLE :: valgpzc(:,:,:)
  REAL,ALLOCATABLE :: vt700(:,:)

  INTEGER,ALLOCATABLE :: undrgrnd(:,:,:)

  REAL,ALLOCATABLE :: vlat2d(:,:),vlon2d(:,:)
  REAL,ALLOCATABLE :: vx_ext(:),vy_ext(:)

  CHARACTER (LEN=256) :: grdbasfn_ext
  CHARACTER (LEN=3)  :: fmtn
  INTEGER :: lenrun, ldir
  INTEGER :: itmp,ireturn

!-----------------------------------------------------------------------
!
! External grid map variables
!
!-----------------------------------------------------------------------

  REAL :: scale
  REAL :: trulat(2)
  REAL :: scswlat,scswlon
  REAL :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
          mapselat,mapselon,mapnelat,mapnelon

!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  INTEGER :: vnx,vny,vnz

  REAL,ALLOCATABLE :: vp_ext(:,:,:)     ! Pressure (Pascals)
  REAL,ALLOCATABLE :: vhgt_ext(:,:,:)   ! Height (m)
  REAL,ALLOCATABLE :: vt_ext(:,:,:)     ! Temperature (K)
  REAL,ALLOCATABLE :: vu_ext(:,:,:)     ! Eastward wind component
  REAL,ALLOCATABLE :: vv_ext(:,:,:)     ! Northward wind component
  REAL,ALLOCATABLE :: vqv_ext(:,:,:)    ! Specific humidity (kg/kg)

  CHARACTER*19 :: extdinit
  CHARACTER*9 :: extdfcst
  CHARACTER*9 :: julfname

  INTEGER :: istatus
  INTEGER :: initsec,kftime,jabssec,iabssec
  INTEGER :: ifhr,ifmin,ifsec,mfhr
  INTEGER :: jldy
  INTEGER :: myr

  REAL :: qvsat

  INTEGER :: i,j,k,l
  INTEGER :: kk

  REAL :: zlevel,z01, p00, t00

  REAL,PARAMETER :: missing = -9999. ! "Underground"/missing flag 

  INTEGER :: ifile,ntime, date(6), init_date(6)

  INTEGER,ALLOCATABLE :: timesec(:)

  REAL :: hgtkm

  REAL :: mapswx,mapswy,mapnwx,mapnwy,mapsex,mapsey,mapnex,mapney

!-----------------------------------------------------------------------
!
! Functions
!
!-----------------------------------------------------------------------

  REAL,EXTERNAL :: f_qvsatl
!fpp$ expand (f_qvsatl)
!dir$ inline always f_qvsatl

!-----------------------------------------------------------------------
!
! Some universal constants
!
!-----------------------------------------------------------------------

  REAL, PARAMETER ::   kappa=287.053/CP,      &
                       gamma=6.5,             & ! 6.5 K/km
                       ex1=0.1903643,         & ! R*gamma/g
                       ex2=5.2558774,         & ! g/R/gamma
                       mbtopa=100.

!-----------------------------------------------------------------------
!
! Namelists
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxfile=50
  INTEGER :: extdopt, extdfmt, use_multiple_names, &
      nxinput,nyinput,nzinput, nextdfil
  CHARACTER*8 :: model
  CHARACTER*4 :: grid
!  CHARACTER*29 :: extdtime(maxfile)
!  CHARACTER*256 :: extdname(maxfile), dir_extd(maxfile), veriffile

  CHARACTER (LEN=256) :: dir_extd(maxfile),extdname(maxfile)
  CHARACTER (LEN=256) :: extdtime(maxfile),veriffile

  CHARACTER (LEN=256) :: dir_extd_temp,extdname_temp,extdtime_temp

  NAMELIST /infile/ extdopt, model, grid,  &
                    nextdfil, extdname, nxinput, nyinput, nzinput,  &
                    extdfmt, dir_extd, extdtime, use_multiple_names
  
  NAMELIST /output/ veriffile
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Set dummy global grid information for GETARPS.  This data will not
! affect what is output to the verification HDF file.
!
!-----------------------------------------------------------------------

  sclfct = 1.0
  trulat1 = 60.0
  trulat2 = 30.0
  trulon = -100.
  ctrlat = 37.
  ctrlat = -92.
  mapproj = 2

  PRINT '(10(A/))', &
 ' ',&
 '   ###################################################################',&
 '   #                                                                 #',&
 '   # Welcome to CVT2VERIF, a program that reads in ARPS history      #',&
 '   # dumps and isobaric RUC and Eta GRIB files (#236 and #212,       #',&
 '   # respectively), and outputs verif data in .hdf format.           #',&
 '   #                                                                 #',&
 '   ###################################################################'

!-----------------------------------------------------------------------
!
! Read namelists
!
!-----------------------------------------------------------------------

  model =''
  grid = ''
  
  PRINT *, 'Reading NAMELIST infile'
  READ (*, infile, ERR=8000, END=8001)
  PRINT *, 'Reading NAMELIST output'
  READ (*, output, ERR=8000, END=8001)
  GO TO 8010

  8000  CONTINUE
  PRINT *, 'FATAL: Error reading Namelist. '
  WRITE (*, NML=infile)
  WRITE (*, NML=output)
  STOP

  8001  CONTINUE
  PRINT *, 'FATAL: End of file reading Namelist. '
  WRITE (*, NML=infile)
  WRITE (*, NML=output)
  STOP

  8010  CONTINUE

  IF (use_multiple_names == 0) THEN
    dir_extd(2:) = dir_extd(1)
    extdname(2:) = extdname(1)
  END IF

!EMK 7/7/2000
  IF (LEN_TRIM(model).EQ.0) THEN
    IF (extdopt.EQ.0) THEN
      model = 'ARPS'
    ELSE IF (extdopt.EQ.2) THEN
      model = 'Eta'
    ELSE IF (extdopt.EQ.11.OR.extdopt.EQ.12) THEN
      model = 'RUC2'
    ELSE
      model = '???'
    END IF
  END IF
  IF (LEN_TRIM(grid).EQ.0) THEN
    IF (extdopt.EQ.0) THEN
      grid = '???'
    ELSE IF (extdopt.EQ.2) THEN
      grid = '212'
    ELSE IF (extdopt.EQ.11.OR.extdopt.EQ.12) THEN
      grid = '236'
    ELSE
      grid = '???'
    END IF
  END IF
!-----------------------------------------------------------------------
!
! Set dimensions
!
!-----------------------------------------------------------------------

  IF (extdopt == 0) THEN  ! Set ARPS dimensions
    nx = nxinput
    ny = nyinput
    nz = nzinput
  ELSE IF (extdopt == 1) THEN ! Set NMC RUC GRIB #87 dimensions
    nx = ruc87nx
    ny = ruc87ny
    nz = ruc87nz
    WRITE(6,*)'ERROR:  RUC GRIB #87 is not supported by this program',&      
              'at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 2) THEN ! Set NMC Eta GRIB #212 dimensions
    nx = eta212nx
    ny = eta212ny
    nz = eta212nz
  ELSE IF (extdopt == 3) THEN ! Set LAPS dimensions
    nx = nxinput
    ny = nyinput
    nz = nzinput
    WRITE(6,*)'ERROR:  LAPS is not supported by this program',        &      
              'at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 4) THEN ! Set NMC RUC GEMPAK dimensions
    nx = rucgemnx
    ny = rucgemny
    nz = rucgemnz
    WRITE(6,*)'ERROR:  RUC GEMPAK is not supported by this program',  &      
              'at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 5) THEN ! Set NMC Eta GEMPAK dimensions
    nx = etagemnx
    ny = etagemny
    nz = etagemnz
    WRITE(6,*)'ERROR:  Eta GEMPAK is not supported by this program',  &      
              'at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 6) THEN ! Set COAMPS dimensions
    nx = nxinput
    ny = nyinput
    nz = nzinput
    WRITE(6,*)'ERROR:  COAMPS is not supported by this program',      &      
              'at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 7) THEN ! Set NCEP RUC2 GRIB #211 dimensions
    nx = ruc211nx
    ny = ruc211ny
    nz = ruc211nz
    WRITE(6,*)'ERROR:  RUC GRIB #211 is not supported by this program',&      
              'at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 8) THEN ! Set NCEP global re-analysis on T62
                              ! Gaussian lat/lon grid dimensions
    nx = glreannx
    ny = glreanny
    nz = glreannz
    WRITE(6,*)'ERROR:  NCEP reanalysis is not supported by this ',    &      
              'program at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 9) THEN ! Set NCEP RUC2 GEMPAK dimensions
    nx = ruc2gemnx
    ny = ruc2gemny
    nz = ruc2gemnz
    WRITE(6,*)'ERROR:  RUC2 GEMPAK is not supported by this ',        &      
              'program at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 10) THEN ! Set NCEP ETA GEMPAK #104 dimensions
    nx = etagemnx
    ny = etagemny
    nz = etagemnz
    WRITE(6,*)'ERROR:  Eta GEMPAK is not supported by this ',        &      
              'program at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 11) THEN ! Set NCEP native coordinate RUC2 
                               ! GRIB #236 dimensions
    nx = rucn236nx
    ny = rucn236ny
    nz = rucn236nz
    WRITE(6,*)'ERROR:  RUC Native coordinate GRIB #236 is not supported',&
              'by this program at this time.  Exiting...'
    STOP
  ELSE IF (extdopt == 12) THEN ! Set NCEP isobaric RUC2 
                               ! GRIB #236 dimensions
    nx = rucp236nx
    ny = rucp236ny
    nz = rucp236nz
  ELSE
    WRITE(6,*)'Invalid option for extdopt.  Exiting program...'
    STOP
  END IF

!-----------------------------------------------------------------------
!
! Allocate arrays.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'Allocating arrays...'

  vnx = nx
  vny = ny
  vnz = nz
  nstyps = 4
  nstyp = nstyps

  IF (extdopt == 0) THEN ! ARPS
    vnx = nx-3
    vny = ny-3

!---------------------------------------------------------------------
!
! Here, even though multiple times are being used, we read in the data
! from the first "grdbas" file, since we are using the same domain for
! the verification statistics and the nx,ny,nz and nsytps won't change.
!
!---------------------------------------------------------------------

  dir_extd_temp=dir_extd(1)
  extdname_temp=extdname(1)

  lenrun=LEN(dir_extd_temp)
  ldir=lenrun
  CALL strlnth( dir_extd_temp, ldir )

  IF ( ldir == 0 .OR. dir_extd_temp(1:ldir) == ' ' ) THEN
    dir_extd_temp = '.'
    ldir = 1
  END IF

  IF( dir_extd_temp(ldir:ldir) /= '/' .AND.  ldir < lenrun ) THEN
    ldir = ldir + 1
    dir_extd_temp(ldir:ldir) = '/'
  END IF

  lenrun = LEN( extdname_temp )
  CALL strlnth( extdname_temp, lenrun )

  IF( extdfmt == 1 ) THEN
    fmtn = 'bin'
  ELSE IF ( extdfmt == 2 ) THEN
    fmtn = 'asc'
  ELSE IF ( extdfmt == 3 ) THEN
    fmtn = 'hdf'
  ELSE IF ( extdfmt == 4 ) THEN
    fmtn = 'pak'
  ELSE IF ( extdfmt == 6 ) THEN
    fmtn = 'bn2'
  ELSE IF ( extdfmt == 7 ) THEN
    fmtn = 'net'
  ELSE IF ( extdfmt == 8 ) THEN
    fmtn = 'npk'
  ELSE IF ( extdfmt == 9 ) THEN
    fmtn = 'gad'
  ELSE IF ( extdfmt == 10 ) THEN
    fmtn = 'grb'
  ELSE
    WRITE(6,'(a,a,a)')                                                  &
        'Unknown format, ', extdfmt, '. Program stopped in CVT2VERIF.'
      STOP
    END IF

    grdbasfn_ext = dir_extd_temp(1:ldir)//extdname_temp(1:lenrun)       &
                                        //'.'//fmtn//'grdbas'

    CALL get_dims_from_data(extdfmt,grdbasfn_ext,nx,ny,nz,nstyps,ireturn)

    nstyp = nstyps

    ALLOCATE(soiltyp_ext(nx,ny,nstyps),stypfrct_ext(nx,ny,nstyps), &
             vegtyp_ext(nx,ny), &
             lai_ext(nx,ny),roufns_ext(nx,ny),veg_ext(nx,ny))
  
    ALLOCATE (tem1(nx,ny,nz), tem2(nx,ny,nz), &
              tem3(nx,ny,nz), tem4(nx,ny,nz), &
              vtem2d(vnx,vny), valgpzc(vnx,vny,vnz), vt700(vnx,vny), &
	      STAT=istatus)
      IF (istatus /= 0) CALL alloc_fail (istatus, 'tem')
    
  ENDIF

  WRITE(6,*)'nx = ',nx,' vnx = ',vnx
  WRITE(6,*)'ny = ',ny,' vny = ',vny
  WRITE(6,*)'nz = ',nz,' vnz = ',vnz

  ALLOCATE (undrgrnd(nx,ny,nz), vp_ext(vnx,vny,vnz), vhgt_ext(vnx,vny,vnz), &
            vt_ext(vnx,vny,vnz), vqv_ext(vnx,vny,vnz), vu_ext(vnx,vny,vnz), &
            vv_ext(vnx,vny,vnz), &
            p_ext(nx,ny,nz), hgt_ext(nx,ny,nz), t_ext(nx,ny,nz), &
            qv_ext(nx,ny,nz), u_ext(nx,ny,nz), v_ext(nx,ny,nz), &
            qc_ext(nx,ny,nz), qr_ext(nx,ny,nz), qi_ext(nx,ny,nz), &
            qs_ext(nx,ny,nz), qh_ext(nx,ny,nz), & 
            tsoil_ext(nx,ny,0:nstyps), qsoil_ext(nx,ny,0:nstyps), &
            wetcanp_ext(nx,ny,0:nstyps), snowcvr_ext(nx,ny), trn_ext(nx,ny), &
            psfc_ext(nx,ny), T_2m_ext(nx,ny), RH_2m_ext(nx,ny), &
            U_10m_ext(nx,ny), V_10m_ext(nx,ny), MSLP_ext(nx,ny), &
            vMSLP_ext(vnx,vny), RH_ext(nx,ny,nz), vRH_ext(vnx,vny,vnz), &
            x_ext(nx), y_ext(ny), lat_ext(nx,ny), lon_ext(nx,ny), &
            var_p(vnx,vny,nlevel,nextdfil,nvar_p),  &
            var_sfc(vnx,vny,nextdfil,nvar_sfc),  &
            timesec(nextdfil),vlat2d(vnx,vny),vlon2d(vnx,vny), &
            vx_ext(vnx),vy_ext(vny), &
            latu_ext(nx,ny),lonu_ext(nx,ny), &
            latv_ext(nx,ny),lonv_ext(nx,ny), &
            zp_ext(nx,ny,nz), w_ext(nx,ny,nz), &
            vatu_ext(nx,ny,nz), uatv_ext(nx,ny,nz), &
            STAT=istatus) 
    IF (istatus /= 0) CALL alloc_fail (istatus, '_ext')

  var_p = missing
  var_sfc = missing

  WRITE(6,*)'Finished allocating arrays...'

!-----------------------------------------------------------------------
!
! Loop through the data times provided via NAMELIST
!
!-----------------------------------------------------------------------
     
loop_ifile:  DO ifile = 1,nextdfil
 
    PRINT *, 'Processing ', ifile, ' ', TRIM(extdname(ifile)), ' ', &
             TRIM(extdtime(ifile))

    runname = extdname(ifile)

! FIX-ME
    nocmnt = 0

!-----------------------------------------------------------------------
!
!   Determine the Julian date for the data.
!
!-----------------------------------------------------------------------

    READ(extdtime(ifile),'(a19,1x,a9)') extdinit,extdfcst
    IF (extdfcst.eq.'         ') extdfcst='000:00:00'
    READ(extdinit, '(I4,5(1X,I2))', ERR=920,END=920) date(:)

    year   = date(1)
    month  = date(2)
    day    = date(3)
    hour   = date(4)
    minute = date(5)
    second = date(6)

    IF (ifile == 1) init_date(:) = date(:)
        
    CALL julday(year, month, day, jldy)
    myr=mod(year,100)
    ifhr=0
    ifmin=0
    ifsec=0
    READ(extdfcst, '(i3,2(1x,i2))',ERR=920,END=920) ifhr,ifmin,ifsec
    mfhr=mod(ifhr,24)
    jldy = jldy + ifhr/24
    WRITE(julfname, '(i2.2,i3.3,2i2.2)') myr,jldy,hour,mfhr
    CALL ctim2abss(year,month,day,hour,minute,second,iabssec)
    jabssec=(ifhr*3600) + (ifmin*60) + ifsec + iabssec
 
    IF (ifile == 1) THEN
      initsec = jabssec
    ENDIF
    kftime=jabssec-initsec
    curtim=float(kftime)

    timesec(ifile) = curtim

    PRINT *, 'Calling rdextfil, looking for ', TRIM(extdfcst), &
        ' hour forecast initialized at ',extdinit, ' UTC'
    PRINT *, 'ARPS delta time is ', kftime, '  abs sec=', jabssec

!-----------------------------------------------------------------------
!
!   Get external data.
!
!-----------------------------------------------------------------------

! Williams SRA
! JJL modification 03/29/02
! Change call from custom "getarps2" to standard "getarps".
! Added several temporary arrays to update call for new ARPS 5.0 code.

    IF (extdopt == 0) THEN
      WRITE(6,*)'Calling getarps...'
       CALL getarps(nx,ny,nz,nzsoil,                                     &
            dir_extd(ifile),extdname(ifile),extdopt,extdfmt,             &
            extdinit,extdfcst,julfname,nstyps,                           &
            iproj_ext,scale_ext,                                         &
            trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
            lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,         &
            p_ext,hgt_ext,zp_ext,zpsoil_ext,t_ext,qv_ext,                &
            u_ext,vatu_ext,v_ext,uatv_ext,w_ext,                         &
            qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                          &
            soiltyp_ext,stypfrct_ext,vegtyp_ext,                         &
            lai_ext,roufns_ext,veg_ext,                                  &
            tsoil_ext,qsoil_ext,wetcanp_ext,                             &
            tem1,tem2,istatus)

      DO k = 1,vnz
        DO j = 1,vny
          DO i = 1,vnx
            vp_ext(i,j,k) = p_ext(i+1,j+1,k)
            vhgt_ext(i,j,k) = hgt_ext(i+1,j+1,k)
            vt_ext(i,j,k) = t_ext(i+1,j+1,k)
            vqv_ext(i,j,k) = qv_ext(i+1,j+1,k)
            vu_ext(i,j,k) = u_ext(i+1,j+1,k)
            vv_ext(i,j,k) = v_ext(i+1,j+1,k)
          END DO
        END DO
      END DO

    ELSE IF (extdopt == 2 .OR. extdopt == 12) THEN
     
      WRITE(6,*)'Calling getgrib for external model'
      CALL getgrib(nx,ny,nz,                                            &
                        dir_extd(ifile),extdname(ifile),extdopt,extdfmt,&
                        extdinit,extdfcst,julfname,                     &
                        iproj_ext,scale_ext,                            &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,             &
                        lat_ext,lon_ext,                                &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,         &
                        qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,             &
                        tsoil_ext(1,1,1,0),tsoil_ext(1,1,2,0),          &
                        qsoil_ext(1,1,1,0),qsoil_ext(1,1,2,0),          &
                        wetcanp_ext,                                    &
                        snowcvr_ext,trn_ext,psfc_ext,                   &
                        T_2m_ext,RH_2m_ext,U_10m_ext,                   &
                        V_10m_ext,MSLP_ext,RH_ext,                      &
                        undrgrnd,istatus)

      vp_ext = p_ext
      vhgt_ext = hgt_ext
      vt_ext = t_ext
      vqv_ext = qv_ext
      vu_ext = u_ext
      vv_ext = v_ext
      vRH_ext = RH_ext
      vMSLP_ext = MSLP_ext

    ENDIF

!-----------------------------------------------------------------------
!
!   Save grid information
!
!-----------------------------------------------------------------------
 
    IF (ifile == 1) THEN
      mapproj = iproj_ext
      scale = scale_ext
      trulat(1) = latnot_ext(1)
      trulat(2) = latnot_ext(2)
      trulon = trlon_ext
      IF (extdopt == 0) THEN
        scswlat = lat_ext(2,2)
        scswlon = lon_ext(2,2)
      ELSE
        scswlat = lat_ext(1,1)
        scswlon = lon_ext(1,1)
      END IF
      CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
      CALL setorig(1,x0_ext,y0_ext)
      DO j = 1,ny
        CALL lltoxy(1,1,lat_ext(1,j),lon_ext(1,j),x_ext(1),y_ext(j))
      END DO
      DO i = 1,nx
        CALL lltoxy(1,1,lat_ext(i,1),lon_ext(i,1),x_ext(i),y_ext(1))
      END DO
      dx = x_ext(2) - x_ext(1)
      dy = y_ext(2) - y_ext(1)

      DO j = 1,vny
        DO i = 1,vnx
          IF (extdopt == 0) THEN
            vlat2d(i,j) = lat_ext(i+1,j+1)
            vlon2d(i,j) = lon_ext(i+1,j+1)
          ELSE
            vlat2d(i,j) = lat_ext(i,j)
            vlon2d(i,j) = lon_ext(i,j)
          END IF
        END DO 
      END DO 

      DO j = 1,vny
        CALL lltoxy(1,1,vlat2d(1,j),vlon2d(1,j),vx_ext(1),vy_ext(j))
      END DO
      DO i = 1,vnx
        CALL lltoxy(1,1,vlat2d(i,1),vlon2d(i,1),vx_ext(i),vy_ext(1))
      END DO

      CALL mcorner(vnx,vny,vx_ext,vy_ext,1,vnx,1,vny,&
                   dx,dy,mapswlat,mapswlon, &
                   mapnwlat,mapnwlon,mapselat,mapselon,mapnelat, &
                   mapnelon)

    ENDIF

! FIX-ME
!    IF (istatus /= 0) THEN
!      PRINT *, 'ERROR: Bad status from getarps or getgrib: ', istatus
!      PRINT *, 'Skipping ifile= ', ifile
!      CYCLE loop_ifile
!    END IF

!-----------------------------------------------------------------------
!
!   Calculate RH from qv for ARPS
!
!-----------------------------------------------------------------------

    IF (extdopt == 0) THEN
      DO k = 1,vnz
        DO j = 1,vny
          DO i = 1,vnx
            IF (vqv_ext(i,j,k) /= missing) THEN
              IF (vp_ext(i,j,k) > 0. .AND. vt_ext(i,j,k) > 0.) THEN
                qvsat = f_qvsatl(vp_ext(i,j,k),vt_ext(i,j,k))
                vRH_ext(i,j,k) = vqv_ext(i,j,k)/qvsat * 100.0
              ELSE
                vRH_ext(i,j,k) = missing
              ENDIF
            ENDIF
          END DO
        END DO
      END DO
    END IF

!-----------------------------------------------------------------------
!
!   Calculate MSLP for ARPS history dump
!
!-----------------------------------------------------------------------

    IF (extdopt == 0) THEN
      DO k = 1,vnz
        DO j = 1,vny
          DO i = 1,vnx
            IF (vp_ext(i,j,k) > 0.) THEN
              valgpzc(i,j,k) = - ALOG(vp_ext(i,j,k))
            ELSE
              valgpzc(i,j,k) = missing
            ENDIF
          END DO
        END DO
      END DO

      zlevel = -ALOG(70000.0)
      CALL hintrp(vnx,vny,vnz,vt_ext,valgpzc,zlevel,vt700)

      DO j = 1,vny
        DO i = 1,vnx
          p00 = 0.01*vp_ext(i,j,2)
          IF(p00 <= 700.0) THEN
            t00=vt_ext(i,j,2)
          ELSE
            t00 = vt700(i,j)*(p00/700.0)**ex1
          ENDIF
          hgtkm = vhgt_ext(i,j,2)*0.001
          vMSLP_ext(i,j) = p00*(1.0+gamma*hgtkm/t00)**ex2
        END DO
      END DO
    ENDIF

!-----------------------------------------------------------------------
!
!   Calculate 2-D variables.  (This may be put into a separate 
!   subroutine.  Also, code will need to be made more flexible in the 
!   future and calculate ONLY those variables that the user requests.)
!
!-----------------------------------------------------------------------

    IF (extdopt == 2 .OR. extdopt == 12) THEN

      var_sfc(:,:,ifile,id_p)  = vMSLP_ext(:,:)*100.
      var_sfc(:,:,ifile,id_t)  = t_2m_ext(:,:)
      var_sfc(:,:,ifile,id_rh) = rh_2m_ext(:,:)
      var_sfc(:,:,ifile,id_u)  = u_10m_ext(:,:)
      var_sfc(:,:,ifile,id_v)  = v_10m_ext(:,:)

      DO l=1,nlevel

!-----------------------------------------------------------------------
!
!       Find the external model's index of the current level.
!
!-----------------------------------------------------------------------

	kk = 1
	DO k = 1,vnz
	  IF (vp_ext(1,1,k) == pressure(l)) kk = k
	END DO

	DO j=1,vny
	DO i=1,vnx
	  IF (undrgrnd(i,j,kk) == 0) THEN
	    var_p(i,j,l,ifile,id_h) = vhgt_ext(i,j,kk)
	    var_p(i,j,l,ifile,id_t) = vt_ext(i,j,kk)
	    var_p(i,j,l,ifile,id_rh) = rh_ext(i,j,kk)
	    var_p(i,j,l,ifile,id_u) = vu_ext(i,j,kk)
	    var_p(i,j,l,ifile,id_v) = vv_ext(i,j,kk)
	  END IF
	END DO
	END DO

      END DO

!-----------------------------------------------------------------------

    ELSE IF (extdopt == 0) THEN ! ARPS
      
      var_sfc(:,:,ifile,id_p)  = vMSLP_ext(:,:)*100.
      var_sfc(:,:,ifile,id_t)  = vt_ext(:,:,2)
      var_sfc(:,:,ifile,id_rh) = vrh_ext(:,:,2)
      var_sfc(:,:,ifile,id_u)  = vu_ext(:,:,2)
      var_sfc(:,:,ifile,id_v)  = vv_ext(:,:,2)

      DO l=1,nlevel
	z01 = - ALOG(pressure(l))
	CALL hintrp (vnx,vny,vnz,vhgt_ext,valgpzc,z01, var_p(1,1,l,ifile,id_h))
	CALL hintrp (vnx,vny,vnz,vt_ext,valgpzc,z01, var_p(1,1,l,ifile,id_t))
	CALL hintrp (vnx,vny,vnz,vrh_ext,valgpzc,z01, var_p(1,1,l,ifile,id_rh))
	CALL hintrp (vnx,vny,vnz,vu_ext,valgpzc,z01, var_p(1,1,l,ifile,id_u))
	CALL hintrp (vnx,vny,vnz,vv_ext,valgpzc,z01, var_p(1,1,l,ifile,id_v))
      END DO
      
    ENDIF

  END DO loop_ifile	! ifile

!-----------------------------------------------------------------------
!
! Now write variables to HDF verification file.
!
!-----------------------------------------------------------------------

  ntime = nextdfil

  CALL wrtverif ( vnx,vny, nlevel,ntime, nvar_p,nvar_sfc, missing,	&
                  veriffile,model,grid,init_date, timesec, pressure,	&
                  mapproj,scale, trulat(1),trulat(2),trulon, dx,dy,     &
                  scswlat,scswlon, &
                  mapswlat,mapswlon,mapnwlat,mapnwlon, &
                  mapselat,mapselon,mapnelat,mapnelon, &
	          flag_p,varid_p, varname_p, varunit_p, var_p,	        &
	          flag_sfc,varid_sfc,varname_sfc,varunit_sfc,var_sfc)

  WRITE(6,*) 'Program CVT2VERIF successfully completed.'
  STOP

!-----------------------------------------------------------------------
!
! Problem doing time conversions.
!
!-----------------------------------------------------------------------

  920 CONTINUE
  PRINT *, 'Aborting, error in time format for external file: ', extdtime(ifile)
  STOP

END PROGRAM cvt2verif
