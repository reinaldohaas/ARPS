PROGRAM arpssoil
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSSOIL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Program to use history data dump and surface property data to
!  create the soil variables
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!    05/04/1995
!
!  MODIFICATION HISTORY:
!    11/10/1995 Keith Brewster
!    Updated for arps4.1.0 I/O routines
!
!    03/08/1996 Keith Brewster
!    Updated for change in parameter list for mksoilvar.
!
!    March, 1997 John Mewes
!    Added API processing
!
!    04/07/97 Keith Brewster
!    Updated for changes to arps4.2.4
!    Restructured API calculations for ARPS-coding compliance.
!
!    04/12/97 Keith Brewster
!    Added API calculation using NCEP near real-time precip data.
!    Renamed arpsh2s package to version 2.0.
!
!    04/22/97 Keith Brewster
!    Added history dump option for more efficient viewing of results.
!    New history dump variables in input file.
!
!    12/09/1998 (Donghai Wang)
!    Added the snow cover.
!
!    1999/07/16 (Gene Bassett)
!    Added option to set water temperature over a lat/lon box.  Added
!    terrain parameters to input file so an input history file is not
!    required.  Added fnldate for api soil moisture generation.
!
!    2000/01/03 (Gene Bassett)
!    Renamed mstinit to soilinit2 and added capability to update soil
!    data over water and snow cover with data in a soil data file.
!
!    2001/05/24 (Gene Bassett)
!    Added soil and surface input format variables soilfmt & sfcfmt.
!
!    1 June 2002 Eric Kemp
!    Soil variable updates.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points for soil model
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
!    qscalar
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    soil temperature (K)
!    qsoil    soil moisture
!    wetcanp  Canopy water amount
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain(mm)
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net shortwave radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Work array for Savi3D dump
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz,nzsoil
  INTEGER :: nstyps,nstyp1

  INTEGER :: mxdays,mxstns  ! API limits

  REAL, ALLOCATABLE :: x      (:)       ! The x-coord. of the physical and
                            ! computational grid.
                            ! Defined at u-point.
  REAL, ALLOCATABLE :: y      (:)       ! The y-coord. of the physical and
                            ! computational grid.
                            ! Defined at v-point.
  REAL, ALLOCATABLE :: z      (:)       ! The z-coord. of the computational
                            ! grid. Defined at w-point.

  REAL, ALLOCATABLE :: zp     (:,:,:) ! Physical height of grid at w-points.
  REAL, ALLOCATABLE :: zpsoil (:,:,:) ! Physical height of grid at w-points.
                                      ! for soil model
  REAL, ALLOCATABLE :: zsc    (:,:,:) ! Physical height of grid at scalar-points.
  REAL, ALLOCATABLE :: hterain(:,:)    ! Terrain height.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                            ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)

  REAL, ALLOCATABLE :: tke    (:,:,:) ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh    (:,:,:) ! Horizontal turb. mixing coef. for
                            ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv    (:,:,:) ! Vertical turb. mixing coef. for
                            ! momentum. ( m**2/s )

  REAL, ALLOCATABLE :: ubar   (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:) ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:) ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:) ! Base state water vapor specific
                            ! humidity (kg/kg)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:) ! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp (:,:)         ! Vegetation type (currently only 1
                                    ! vegetation type is used in the arps)
  REAL, ALLOCATABLE :: lai    (:,:)         ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)         ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)         ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:)    ! Deep soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:)    ! soil moisture
  REAL, ALLOCATABLE :: wetcanp(:,:,:)    ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)               ! Snow depth (:)

  REAL, ALLOCATABLE :: raing(:,:)         ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)         ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(:,:,1) = total precip. rate
                               ! prcrate(:,:,2) = grid scale precip. rate
                               ! prcrate(:,:,3) = cumulus precip. rate
                               ! prcrate(:,:,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)      ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)        ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)        ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)      ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)       ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  API generation arrays:
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=7),   ALLOCATABLE :: staid(:)
  CHARACTER(LEN=256), ALLOCATABLE :: precfile(:)
  INTEGER, ALLOCATABLE :: iptstn(:)
  INTEGER, ALLOCATABLE :: jptstn(:)
  INTEGER, ALLOCATABLE :: iptapi0(:)
  INTEGER, ALLOCATABLE :: jptapi0(:)
  REAL, ALLOCATABLE :: api0(:)
  REAL, ALLOCATABLE :: obsprec(:,:)
  REAL, ALLOCATABLE :: difprec(:,:)
  REAL, ALLOCATABLE :: totpreca(:)

  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: prec(:,:)
!
!-----------------------------------------------------------------------
!
!  Temporary working arrays:
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: itema(:,:)
  REAL, ALLOCATABLE :: initapi(:,:,:)
  REAL, ALLOCATABLE :: api1(:,:,:)
  REAL, ALLOCATABLE :: api2(:,:,:)
  REAL, ALLOCATABLE :: kk2dep(:,:)
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
!
!-----------------------------------------------------------------------
!
!  API parameters
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile
  CHARACTER (LEN=256) :: soiloutfl
  CHARACTER (LEN=256) :: apifile
  CHARACTER (LEN=256) :: apiinit
  CHARACTER (LEN=256) :: prcpdir
  CHARACTER (LEN=256) :: prcplst

  ! J.Case, ENSCO Inc. -- 9/16/2004 -- Changed to include hour for API.
  CHARACTER (LEN=13) :: inidate
  CHARACTER (LEN=13) :: fnldate

  INTEGER :: hinfmt,soilinit2,prdata
  REAL :: respapi,k1,k2,range,gamma(2),respprec
  CHARACTER (LEN=80) :: runnamesv
  INTEGER :: nbwater,nbw_max
  PARAMETER (nbw_max=128)
  REAL :: tbwater(nbw_max),blat1(nbw_max),blat2(nbw_max)
  REAL :: blon1(nbw_max),blon2(nbw_max)
  COMMON /bwcomi/ nbwater
  COMMON /bwcomr/ tbwater,blat1,blat2,blon1,blon2
  REAL :: twater
  COMMON /realcom/ twater

  INTEGER :: zpsoilout, tsoilout, qsoilout, wcanpout,snowdout

!  COMMON /intgcom/tsfcout,tsoilout,wsfcout,wdpout,wcanpout,snowdout
  COMMON /intgcom/ zpsoilout,tsoilout,qsoilout,wcanpout,snowdout

!  NAMELIST /arpssoilvars/ hinfmt, grdbasfn, hisfile,                    &
!            soilinit2,prdata,apiinit,apifile,prcpdir,prcplst,           &
!            sfcdat,styp,wetcanp0,snowdpth0,twater,ttprt,tbprt,          &
!            wgrat,w2rat,sfcdtfl,soilinfl,soiloutfl,                     &
!            tsfcout,tsoilout,wsfcout,wdpout,wcanpout,snowdout,          &
!            inidate,fnldate,respapi,k1,k2,range,                        &
!            respprec,gamma,soildmp,nx,ny,mxdays,mxstns,nstyp,  &
!            soilfmt,sfcfmt
  NAMELIST /arpssoilvars/ hinfmt, grdbasfn, hisfile,                    &
            soilinit2,prdata,apiinit,apifile,prcpdir,prcplst,           &
            sfcdat,styp,wetcanp0,snowdpth0,twater,ttprt,tbprt,          &
            wgrat,w2rat,sfcdtfl,soilinfl,soiloutfl,                     &
            zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,              &
            inidate,fnldate,respapi,k1,k2,range,                        &
            respprec,gamma,soildmp,nx,ny,mxdays,mxstns,nstyp,           &
            soilfmt,sfcfmt

  NAMELIST /box_water/ nbwater,tbwater,blat1,blat2,blon1,blon2

  NAMELIST /history_dump/ hdmpopt,hdmpfmt,grbpkbit,thisdmp,             &
            tstrtdmp,numhdmp,hdmptim,hdfcompr

  NAMELIST /output/ runname,dirname,tfmtprt,exbcdmp,extdadmp,           &
            grdout,basout,varout,mstout,rainout,iceout,                 &
            tkeout, trbout,sfcout,landout,totout,trstout,tmaxmin,       &
            tenergy,imgopt, timgdmp,pltopt,tplots,filcmprs

  NAMELIST /grid_info/ dx,dy,dz,ctrlat,ctrlon,                          &
            mapproj,trulat1,trulat2,trulon,sclfct

  COMMON /apicom1/ apiinit,apifile,prcpdir,prcplst,inidate,fnldate
  COMMON /apicom2/ soilinit2,prdata
  COMMON /apicom3/ respapi,k1,k2,range,respprec,gamma
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: basdmpfn
  INTEGER :: nchin
  INTEGER :: grdbas

  INTEGER :: i,j,k
  INTEGER :: length,lengbfn,lbasdmpf,lengfn,lenout
  INTEGER :: ireturn

  LOGICAL :: iexist
  INTEGER :: is , nxy, nxyz, nsize
  REAL :: amax, amin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = 1
  nestgrd = 0
!
  hinfmt=1
  grdbasfn='may07.bingrdbas'
  hisfile='may07.bin000000'
  soilinit2=2
  soilfmt=1
  prdata=2
  apiinit='initapi.lst'
  apifile='precfile.lst'
  prcpdir='/vortex/precip/'
  prcplst='/vortex/precip/splains.stns'
  sfcdat=3
  sfcfmt=1
  styp=3
  wetcanp0=0.
  snowdpth0=0
  twater=21.
  ttprt=1.
  tbprt=0.
  wgrat=0.5
  w2rat=0.5
  nstyp = 1
  sfcdtfl='arpssfc.dat'
  soilinfl='arpssoil.dat'
  soiloutfl='arpssoil.dat'
  soildmp = 1
  zpsoilout=1
  tsoilout=1
  qsoilout=1
  wcanpout=1
  snowdout=1

! J.Case, ENSCO Inc. -- 9/16/2004 -- Changed to include hour in API dates.
  inidate='1977-05-20-00'
  fnldate='1977-05-20-00'

  respapi=1.0
  k1=0.75
  k2=0.995
  range=250.0
  respprec=0.1
  gamma(1)=1.00
  gamma(2)=1.00
  nbwater = 0

  hdmpopt  = 1
  hdmpfmt  = 10
  grbpkbit = 16
  hdfcompr = 0
  thisdmp  = 3600.0
  tstrtdmp = 0.0
  numhdmp  = 1
  DO i=1,numhdmp
    hdmptim(i) = 0.
  END DO

  runname  = 'X'
  dirname  = ' '
  exbcdmp  = 0
  extdadmp = 0
  grdout   = 0
  basout   = 0
  varout   = 1
  mstout   = 0
  rainout  = 0
  prcout   = 0
  iceout   = 0
  totout   = 1
  tkeout   = 0
  trbout   = 0
  sfcout   = 1
  landout  = 1
  radout   = 0
  flxout   = 0
  trstout  = 3600.0

!
!-----------------------------------------------------------------------
!
!  Set the control parameters for input:
!
!-----------------------------------------------------------------------
!
  WRITE (6,'(/a/)')                                                     &
      ' Input control parameters read in are:'

  READ(5,arpssoilvars, END=998)
!
  WRITE (6, '(/1x,a)')      '&arpssoilvars'
  WRITE (6, '(3x,a,i4,a)')    'hinfmt   = ', hinfmt, ','

  lengbfn = LEN( grdbasfn )
  CALL strlnth(  grdbasfn, lengbfn )
  WRITE (6, '(3x,a,a,a)')                                               &
                     'grdbasfn  = ''', grdbasfn(1:lengbfn), ''','

  lengfn = LEN( hisfile )
  CALL strlnth(  hisfile, lengfn )
  WRITE (6, '(3x,a,a,a)')                                               &
                     'hisfile   = ''', hisfile(1:lengfn),  ''','

  WRITE (6, '(3x,a,i4,a)')    'soilinit2 = ', soilinit2,  ','
  WRITE (6, '(3x,a,i4,a)')    'soilfmt =   ', soilfmt,    ','
  WRITE (6, '(3x,a,i4,a)')    'prdata    = ', prdata ,  ','

  length = LEN( apiinit )
  CALL strlnth( apiinit,length )
  WRITE (6, '(3x,a,a,a)')     'apiinit   = ', apiinit(1:length),        &
      ','

  length = LEN( apifile)
  CALL strlnth( apifile,length )
  WRITE (6, '(3x,a,a,a)')     'apifile   = ', apifile(1:length),        &
      ','

  length = LEN( prcpdir )
  CALL strlnth( prcpdir,length )
  WRITE (6, '(3x,a,a,a)')     'prcpdir   = ', prcpdir(1:length),        &
      ','

  WRITE (6, '(3x,a,i4,a)')    'sfcdat    = ', sfcdat,   ','
  WRITE (6, '(3x,a,i4,a)')    'sfcfmt    = ', sfcfmt,   ','
  WRITE (6, '(3x,a,i4,a)')    'styp      = ', styp,     ','
  WRITE (6, '(3x,a,e15.5,a)') 'wetcanp0  = ', wetcanp0, ','
  WRITE (6, '(3x,a,e15.5,a)') 'snowdpth0 = ', snowdpth0,','
  WRITE (6, '(3x,a,e15.5,a)') 'twater    = ', twater,   ','
  WRITE (6, '(3x,a,e15.5,a)') 'ttprt     = ', ttprt,    ','
  WRITE (6, '(3x,a,e15.5,a)') 'tbprt     = ', tbprt,    ','
  WRITE (6, '(3x,a,e15.5,a)') 'wgrat     = ', wgrat,    ','
  WRITE (6, '(3x,a,e15.5,a)') 'w2rat     = ', w2rat,    ','

  nstyp = MAX(1,nstyp)
  nstyps = nstyp
  WRITE (6, '(3x,a,i4,a)') 'nstyp     = ', nstyp,   ','

  length = LEN( sfcdtfl )
  CALL strlnth( sfcdtfl,length )
  WRITE (6, '(3x,a,a,a)')                                               &
                     'sfcdtfl   = ''', sfcdtfl(1:length), ''','

  lenout = LEN( soilinfl )
  CALL strlnth( soilinfl,length )
  WRITE (6, '(3x,a,a,a)')                                               &
                     'soilinfl  = ''', soilinfl(1:length), ''','

  lenout = LEN( soiloutfl )
  CALL strlnth( soiloutfl,lenout )
  WRITE (6, '(3x,a,a,a)')                                               &
                     'soiloutfl = ''', soiloutfl(1:lenout), ''','

  WRITE (6, '(3x,a,i4,a)')    'soildmp   = ', soildmp,  ','

  WRITE (6, '(3x,a,i4,a)')    'tsoilout  = ', tsoilout, ','
  WRITE (6, '(3x,a,i4,a)')    'qsoilout  = ', qsoilout, ','
  WRITE (6, '(3x,a,i4,a)')    'wcanpout  = ', wcanpout, ','

  WRITE (6, '(3x,a,a,a)')     'inidate   = ', inidate,  ','
  WRITE (6, '(3x,a,a,a)')     'fnldate   = ', fnldate,  ','
  WRITE (6, '(3x,a,f5.3,a)')  'respapi   = ', respapi,  ','
  WRITE (6, '(3x,a,f5.3,a)')  'k1        = ', k1,       ','
  WRITE (6, '(3x,a,f5.3,a)')  'k2        = ', k2,       ','
  WRITE (6, '(3x,a,f7.2,a)')  'range (km)= ', range,    ','
  WRITE (6, '(3x,a,f5.3,a)')  'respprec  = ', respprec, ','
  WRITE (6, '(3x,a,f5.3,a)')  'gamma(1)  = ', gamma(1), ','
  WRITE (6, '(3x,a,f5.3,a)')  'gamma(2)  = ', gamma(2), ','
  WRITE (6, '(1x,a)')       '&END'

  READ (5,box_water, END=998)
!
  WRITE (6, '(/1x,a)')      '&box_water'
  WRITE (6, '(3x,a,i4,a)')  'nbwater   = ', nbwater,       ','
  IF (nbwater > nbw_max) THEN
    WRITE (*,*) "ERROR in arpsh2s: nbwater > nbw_max."
    WRITE (*,*) "     nbwater =",nbwater,"   nbw_max =",nbw_max
    STOP
  END IF
  DO i = 1, nbwater
    WRITE (6, '(3x,a,i3,a,f7.3,a)')                                     &
                        'tbwater(',i,')  = ', tbwater(i), ','
    WRITE (6, '(3x,a,i3,a,f8.3,a)')                                     &
                        'blat1(',i,')  = ', blat1(i), ','
    WRITE (6, '(3x,a,i3,a,f8.3,a)')                                     &
                        'blat2(',i,')  = ', blat2(i), ','
    WRITE (6, '(3x,a,i3,a,f8.3,a)')                                     &
                        'blon1(',i,')  = ', blon1(i), ','
    WRITE (6, '(3x,a,i3,a,f8.3,a)')                                     &
                        'blon2(',i,')  = ', blon2(i), ','
  END DO
  WRITE (6, '(1x,a)')       '&END'

  READ (5,history_dump,END=998)


  IF( hdmpfmt < 0 .OR. hdmpfmt > 10) THEN

    WRITE(6,'(/5x,a,/5x,a/)')                                           &
         'Invalid input values for history dump format option,',        &
         'Job stopped in ARPSH2S.'
    STOP

  END IF
!
  WRITE(6,'(5x,a,i6)')                                                  &
       'Number of bits in packing GRIB dump data was ',grbpkbit
!
  WRITE(6,'(5x,a,i6)')                                                  &
       'HDF4 compression option was ',hdfcompr
!
!-----------------------------------------------------------------------
!
!  Input model output parameters:
!
!  First, give the name of the directory into which output files
!  will be written:
!
!-----------------------------------------------------------------------
!
  READ(5,output, END=998)

  ldirnam = 256
  CALL strlnth( dirname, ldirnam)

  IF( ldirnam == 0 ) THEN
    dirname = '.'
    ldirnam=1
  END IF

  IF( dirname(1:ldirnam) /= ' ') THEN

    READ(5,grid_info, END=998)
!
!-----------------------------------------------------------------------
!
!  Check if the specified output directory exists, if not,
!  create it.
!
!-----------------------------------------------------------------------
!
    INQUIRE(FILE=dirname(1:ldirnam),EXIST=iexist)

    IF( .NOT.iexist ) THEN

      WRITE(6,'(/5x,a,2(/5x,a))')                                       &
          'Specified output directory '//dirname(1:ldirnam)//           &
          ' not found.',                                                &
          'It was created by the program.'

      CALL unixcmd( 'mkdir -p '//dirname(1:ldirnam) )

    END IF

    WRITE(6,'(/5x,a)')                                                  &
        'Output files will be in directory '//dirname(1:ldirnam)//'.'

  ELSE

    WRITE(6,'(/5x,a)')                                                  &
        'Output files will be in the current work directory.'

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the control parameters for the output of selected fields.
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(/5x,a,i4/)')                                                &
       'The input grid coordinate dump option was ', grdout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input base state array dump option was ', basout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input mass-velocity array dump option was ', varout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input non-ice water array dump option was ',mstout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input rain array dump option was ', rainout
  rainout = rainout * mstout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input precipitation rates array dump option was ',prcout
  prcout = prcout * mstout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input ice array dump option was ', iceout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input TKE dump option was ', tkeout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input eddy mixing coeff dump option was ', trbout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input surface variable dump option was ', sfcout

  WRITE(6,'(/5x,a,i4/)')                                                &
      'The input surface property array dump option was ', landout

!
!-----------------------------------------------------------------------
!
!  Allocate arrays.
!
!-----------------------------------------------------------------------
!
  IF (hinfmt == 0) THEN
    nz = nstyps+1
  ELSE
    CALL get_dims_from_data(hinfmt,trim(hisfile),                    &
         nx,ny,nz,nzsoil,nstyp1, ireturn)
  ENDIF

  ALLOCATE(x(nx))
  ALLOCATE(y(ny))
  ALLOCATE(z(nz))

  ALLOCATE(zp(nx,ny,nz))
  ALLOCATE(zpsoil(nx,ny,nzsoil))
  ALLOCATE(zsc(nx,ny,nz))
  ALLOCATE(hterain(nx,ny))

  ALLOCATE(uprt(nx,ny,nz))
  ALLOCATE(vprt(nx,ny,nz))
  ALLOCATE(wprt(nx,ny,nz))
  ALLOCATE(ptprt(nx,ny,nz))
  ALLOCATE(pprt(nx,ny,nz))
  ALLOCATE(qvprt(nx,ny,nz))
  ALLOCATE(qscalar(nx,ny,nz,nscalar))

  ALLOCATE(tke(nx,ny,nz))
  ALLOCATE(kmh(nx,ny,nz))
  ALLOCATE(kmv(nx,ny,nz))

  ALLOCATE(ubar(nx,ny,nz))
  ALLOCATE(vbar(nx,ny,nz))
  ALLOCATE(wbar(nx,ny,nz))
  ALLOCATE(ptbar(nx,ny,nz))
  ALLOCATE(pbar(nx,ny,nz))
  ALLOCATE(rhobar(nx,ny,nz))
  ALLOCATE(qvbar(nx,ny,nz))

  ALLOCATE(soiltyp(nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp(nx,ny))
  ALLOCATE(lai(nx,ny))
  ALLOCATE(roufns(nx,ny))
  ALLOCATE(veg(nx,ny))

  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(raing(nx,ny))
  ALLOCATE(rainc(nx,ny))
  ALLOCATE(prcrate(nx,ny,4))

  ALLOCATE(radfrc(nx,ny,nz))
  ALLOCATE(radsw(nx,ny))
  ALLOCATE(rnflx(nx,ny))
  ALLOCATE(radswnet(nx,ny))
  ALLOCATE(radlwin(nx,ny))

  ALLOCATE(usflx(nx,ny))
  ALLOCATE(vsflx(nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))

  ALLOCATE(staid(mxstns))
  ALLOCATE(precfile(mxstns))
  ALLOCATE(iptstn(mxstns))
  ALLOCATE(jptstn(mxstns))
  ALLOCATE(iptapi0(mxstns))
  ALLOCATE(jptapi0(mxstns))
  ALLOCATE(api0(mxstns))
  ALLOCATE(obsprec(mxdays,mxstns))
  ALLOCATE(difprec(mxdays,mxstns))
  ALLOCATE(totpreca(mxstns))

  ALLOCATE(xs(nx))
  ALLOCATE(ys(ny))
  ALLOCATE(prec(nx,ny))

  ALLOCATE(itema(nx,ny))
  ALLOCATE(initapi(nx,ny,nstyps))
  ALLOCATE(api1(nx,ny,nstyps))
  ALLOCATE(api2(nx,ny,nstyps))
  ALLOCATE(kk2dep(nx,ny))
  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))
!
!-----------------------------------------------------------------------
!
!  Fill all arrays with zero, which is the default value of the
!  data arrays.
!
!-----------------------------------------------------------------------
!
  nxy = nx*ny
  nxyz= nx*ny*nz

  CALL flzero(x    ,nx)
  CALL flzero(y    ,ny)
  CALL flzero(z    ,nz)
  CALL flzero(zp   ,nxyz)
  CALL flzero(hterain,nxy)

  CALL flzero(uprt ,nxyz)
  CALL flzero(vprt ,nxyz)
  CALL flzero(wprt ,nxyz)
  CALL flzero(ptprt,nxyz)
  CALL flzero(pprt ,nxyz)
  CALL flzero(qvprt,nxyz)
  CALL flzero(qscalar   ,nscalar*nxyz)

  CALL flzero(ubar ,nxyz)
  CALL flzero(vbar ,nxyz)
  CALL flzero(wbar ,nxyz)
  CALL flzero(ptbar,nxyz)
  CALL flzero(pbar ,nxyz)
  CALL flzero(rhobar,nxyz)
  CALL flzero(qvbar,nxyz)

  CALL flzero(tke  ,nxyz)
  CALL flzero(kmh  ,nxyz)
  CALL flzero(kmv  ,nxyz)

  CALL flzero(tem1 ,nxyz)
  CALL flzero(tem2 ,nxyz)

  nsize=nxy*nzsoil*(1+nstyps)
  CALL flzero(tsoil  ,nsize)
  CALL flzero(qsoil  ,nsize)

  nsize=nxy*(1+nstyps)
  CALL flzero(wetcanp,nsize)

  IF ( hinfmt > 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!

    ldirnam=LEN(dirname)
    CALL strlnth( dirname , ldirnam)
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    CALL gtlfnkey( runname, lfnkey )
    IF ( lfnkey > 0 ) THEN
      runnamesv(1:80) =  runname(1:80)
    END IF

    CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
                 hinfmt, nchin,grdbasfn(1:lengbfn),lengbfn,             &
                 hisfile(1:lengfn),lengfn,curtim,                       &
                 x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,             &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1,tem2, tem3)

    IF ( lfnkey > 0 ) THEN
      runname(1:80) =  runnamesv(1:80)
    END IF

    CALL gtlfnkey( runname, lfnkey )

    length = lfnkey + 9

    IF( ireturn /= 0 .AND. hinfmt /= 9 ) GO TO 997
                                            ! Read was unsuccessful
  END IF
!
!-----------------------------------------------------------------------
!
!  Set up coordinates for mksoilvar.
!
!-----------------------------------------------------------------------
!
  IF ( hinfmt > 0 ) THEN

    dx=x(2)-x(1)
    dy=y(2)-y(1)

    DO j=1,ny-1
      DO i=1,nx-1
        hterain(i,j)=zp(i,j,2)
      END DO
    END DO

  ELSE
!
!-----------------------------------------------------------------------
!
!  Set terrain and grid info necessary for api.
!
!-----------------------------------------------------------------------
!
    DO i=1,nx
      x(i) = dx * (i-2)
    END DO

    DO j=1,ny
      y(j) = dy * (j-2)
    END DO

!    Create fake z,zp,zsc to keep program from bombing.

    DO k=1,nz
      z(k) = dz * (k-2)
      DO j=1,ny
        DO i=1,nx
          zp(i,j,k) = z(k)
        END DO
      END DO
    END DO

  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zsc(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO

  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=xs(nx-1)+dx
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=ys(ny-1)+dy

  DO is = 0,nstyp

    WRITE(6,'(/a,i4/)') ' In ARPSH2S, for istype =', is

    CALL a3dmax0(tsoil(1,1,1,is),                                          &
        1,nx,1,nx-1,1,ny,1,ny-1,1,1 ,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')  &
            'tsoil_sfcmin = ', amin,', tsoil_sfcmax  =',amax

    CALL a3dmax0(tsoil(1,1,2,is),                                         &
        1,nx,1,nx-1,1,ny,1,ny-1,1,1 ,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')  &
            'tsoil_dpmin= ', amin,', tsoil_dpmax =',amax

    CALL a3dmax0(qsoil(1,1,1,is),                                        &
        1,nx,1,nx-1,1,ny,1,ny-1,1,1 ,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')  &
            'qsoil_wetsmin = ', amin,', qsoil_wetsmax  =',amax

    CALL a3dmax0(qsoil(1,1,2,is),                                         &
        1,nx,1,nx-1,1,ny,1,ny-1,1,1 ,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')  &
            'qsoil_wetdmin = ', amin,', qsoil_wetdmax  =',amax

    CALL a3dmax0(wetcanp(1,1,is),                                       &
        1,nx,1,nx-1,1,ny,1,ny-1,1,1 ,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')  &
            'wetcmin = ', amin,', wetcmax  =',amax

  END DO

  IF(nstyp == 1) THEN
    DO k=1, nzsoil
      DO j=1,ny-1
        DO i=1,nx-1
          tsoil(i,j,k,1) = tsoil(i,j,k,0)
          qsoil(i,j,k,1) = qsoil(i,j,k,0)
        END DO
      END DO
    END DO
    DO j=1,ny-1
      DO i=1,nx-1
        wetcanp(i,j,1) = wetcanp(i,j,0)
      END DO
    END DO
  END IF

  CALL mksoilvar(nx,ny,nz,nzsoil,nstyps,                                &
                 mxdays,mxstns,                                         &
                 xs,ys,                                                 &
                 hterain,zsc,                                           &
                 ptbar,ptprt,                                           &
                 pbar,pprt,                                             &
                 zpsoil,tsoil,qsoil,wetcanp,snowdpth,                   &
                 soiltyp,stypfrct,                                      &
                 precfile,staid,iptstn,jptstn,iptapi0,jptapi0,          &
                 api0,obsprec,difprec,totpreca,prec,                    &
                 itema,initapi,api1,api2,kk2dep,                        &
                 tem1,tem2,tem3)

  DO k=1,nzsoil
    DO j=1,ny-1
      DO i=1,nx-1
        tsoil  (i,j,k,0) = 0.0
        qsoil  (i,j,k,0) = 0.0
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      wetcanp(i,j,0) = 0.0
    END DO
  END DO
!
  PRINT*,' nstyp = ', nstyp

  DO is = 1,nstyp
    DO k=1,nzsoil
      DO j=1,ny-1
        DO i=1,nx-1
          tsoil  (i,j,k,0) = tsoil  (i,j,k,0)                            &
                       + tsoil  (i,j,k,is) * stypfrct(i,j,is)
          qsoil  (i,j,k,0) = qsoil  (i,j,k,0)                            &
                       + qsoil  (i,j,k,is) * stypfrct(i,j,is)
        END DO
      END DO
    END DO
  END DO

  DO is = 1,nstyp
    DO j=1,ny-1
      DO i=1,nx-1
        wetcanp(i,j,0) = wetcanp(i,j,0)                                 &
                       + wetcanp(i,j,is) * stypfrct(i,j,is)
      END DO
    END DO
  END DO

  DO is = 0,nstyp

    WRITE(6,'(/a,i4/)') ' In ARPSH2S, for istype =', is

    CALL a3dmax0(tsoil(1,1,1,is),                                          &
                 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') &
            'tsoil_sfcmin = ', amin,', tsoil_sfcmax  =',amax

    CALL a3dmax0(tsoil(1,1,2,is),                                         &
                 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') &
            'tsoil_dpmin= ', amin,', tsoil_dpmax =',amax

    CALL a3dmax0(qsoil(1,1,1,is),                                        &
                 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') &
            'qsoil_wetsmin = ', amin,', qsoil_wetsmax  =',amax

    CALL a3dmax0(qsoil(1,1,1,is),                                         &
                 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') &
            'qsoil_wetdmin = ', amin,', qsoil_wetdmax  =',amax

    CALL a3dmax0(wetcanp(1,1,is),                                       &
                 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') &
            'wetcmin = ', amin,', wetcmax  =',amax

  END DO


!
!-----------------------------------------------------------------------
!
!  Write the soil data.
!
!-----------------------------------------------------------------------
!
  IF (soildmp > 0) CALL wrtsoil(nx,ny,nzsoil,nstyps,soiloutfl,dx,dy,    &
               zpsoil,mapproj,trulat1,trulat2,trulon,sclfct,            &
               ctrlat,ctrlon,zpsoilout,tsoilout,qsoilout,wcanpout,      &
               snowdout,tsoil,qsoil,wetcanp,snowdpth,soiltyp)
!
!-----------------------------------------------------------------------
!
!  Create a GrADS control file for soil data display
!
!-----------------------------------------------------------------------
!
! J.Case, ENSCO Inc. (11/8/2004) -- found typos.  Original subroutine
!  does not pass the zpsoil variable.

  IF (soildmp == 1) CALL soilcntl(nx,ny,nzsoil,zpsoil,soiloutfl,        &
                zpsoilout,tsoilout,qsoilout,wcanpout,                   &
                snowdout,x,y)
!
!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        qvprt(i,j,k)=qvbar(i,j,k)+qvprt(i,j,k)
        uprt(i,j,k)=ubar(i,j,k)+uprt(i,j,k)
        vprt(i,j,k)=vbar(i,j,k)+vprt(i,j,k)
        wprt(i,j,k)=wbar(i,j,k)+wprt(i,j,k)
      END DO
    END DO
  END DO
!
  IF(hdmpopt /= 0) THEN
    IF( hdmpfmt == 5 ) GO TO 700
    IF( hdmpfmt == 9 ) GO TO 700
!
!-----------------------------------------------------------------------
!
!  Data dump of the model time dependent arrays at initial time:
!
!-----------------------------------------------------------------------
!
    CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,             &
               mgrid,nestgrd,basdmpfn,lbasdmpf)

    WRITE(6,'(1x,a,a)')                                                 &
        'Data dump of grid and base state arrays into file ',           &
        basdmpfn(1:lbasdmpf)

    grdbas = 1                ! Dump out grid and base state arrays only

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                &
               hdmpfmt,nchdmp,basdmpfn(1:lbasdmpf),                     &
               grdbas,filcmprs,                                         &
               uprt,vprt,wprt,ptprt,pprt,qvprt,                         &
               qscalar,tke,kmh,kmv,                                     &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                  &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  Data dump of the model time dependent arrays
!  Find a unique name hdmpfn(1:ldmpf) for history dump data set
!
!-----------------------------------------------------------------------
!
    700   CONTINUE
!

    CALL gtdmpfn(runname(1:lfnkey),dirname,                             &
                 ldirnam,curtim,hdmpfmt,                                &
                 mgrid,nestgrd, hdmpfn, ldmpf)

    WRITE(6,'(1x,a,a)') 'History data dump in file ',hdmpfn(1:ldmpf)

    grdbas = 0                ! No base state or grid array is dumped.

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                &
               hdmpfmt,nchdmp,hdmpfn(1:ldmpf),                          &
               grdbas,filcmprs,                                         &
               uprt,vprt,wprt,ptprt,pprt,qvprt,                         &
               qscalar,tke,kmh,kmv,                                          &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                  &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3)

  END IF

  STOP

  997   CONTINUE
  WRITE(6,'(1x,a,i2,/1x,a)')                                            &
      'Data read was unsuccessful. ireturn =', ireturn,                 &
      'Job stopped.'
  STOP

  998   WRITE(6,'(a)')                                                  &
          'Error reading NAMELIST file. Program ARPSH2S stopped.'
  STOP

  801   FORMAT(a)

END PROGRAM arpssoil
