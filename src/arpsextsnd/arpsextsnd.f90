PROGRAM arpsextsnd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM EXTSND                      ######
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
!  Extracts soundings from ARPS model history data.
!  User specifies an ARPS history file and either a lat,lon or
!  an x-y location.
!
!  Creates a file with the sounding data in an ASCII format.
!  The format matches one used by UNIDATA's GEMPAK software for ASCII
!  output of sounding data (GEMPAK program SNLIST).
!
!  Because of that compatibility, the output file can be
!  read by certain sounding plotting and analysis programs on the
!  computers at the University of Oklahoma (OU).
!
!  To create skewt plot of the generated sounding:
!  bin/skewtpost(or skewtncar) -sfc -hodo data_file
!
!  where data_file is the output of this program (a name specified
!  in the input of this program).  An PostScript or NCARgraphics
!  gmeta file is created by skewt.
!
!  It is also possible to use this file as input to the standard
!  skewt plotting programs on the metgem and Rossby computers.
!  See the author for details.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  October, 1992 (K. Brewster)
!  Conversion to ARPS 3.0.
!
!  February, 1993 (K. Brewster)
!  Additional documentation for ARPS 3.1 release
!
!  4/13/93 (M. Xue)
!  Modified to conform to the new data dump format.
!
!  10/19/94 (KB)
!  Modified to conform to yet another data dump format.
!  Corrections made for change in x,y definition in dump file.
!
!  06/19/95 (KB)
!  Modified documentation and output formats to update info.
!
!  03/07/96 (KB)
!  Changed calling arguments in FINDLC, moved subroutines to a library.
!  Corrected a bug in ymap(ny) calculation.
!
!  11/06/96 (KB)
!  Unified interpolation calls with other programs using interpolation.
!  Added capability to process multiple soundings in one run.
!
!  2000/05/19 (Gene Bassett)
!  Converted to F90, creating allocation and arpsextsnd main
!  subroutines.
!
!  2000/07/28 (Ming Xue)
!  Converted to F90 free format. Use ALLOCATABLE instead of
!  POINTER allocation to avoid double memory usage.
!
!  Changed to namelist format input.
!
!  2001/12/05 (Ming Hu)
!    Add an output file which is sounding used in
!    ARPS as base field like May20.snd
!
!  05/28/2002 (J. Brotzge)
!    Added tsoil/qsoil to accomodate new soil scheme.
!
!  1 June 2002 Eric Kemp
!    Soil variable updates.
!
!  2005/03/30 (Kevin W. Thomas)
!    MPI this program so that large domains get use splitfiles and still
!    get soundings.
!
!  04/10/2005 (Keith Brewster)
!    Added vertical velocity variable to GEMPAK sounding file to
!    accomodate NWS NSHARP program.
!
!  07/26/2006 (Yunheng Wang)
!    Added RAOB sounding format and extracting multiple ARPS grids
!    sounding.
!
!  09/10/2008 (Keith Brewster)
!    Added processing for radar profile for use by the radar
!    remapper programs (88d2arps, ncrad2arps, etc).
!
!  08/04/2011 (Keith Brewster)
!    Added option for extracting qv, clouds and TKE for use by WFIP.
!
!  05/08/2012 (Y. Wang)
!  Added capability to read command line for namelist file.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points in the soil (m)
!
!    uprt     Perturbation x component of velocity (m/s)
!    vprt     Perturbation y component of velocity (m/s)
!    wprt     Perturbation z component of velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!    qscalar  Hydrometeor scalars
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
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
!  CALCULATED DATA ARRAYS:
!
!    su       Sounding x component of velocity (m/s)
!    sv       Sounding y component of velocity (m/s)
!    sw       Sounding w component of velocity (m/s)
!    stheta   Sounding potential temperature (K)
!    spres    Sounding pressure (mb)
!    stemp    Sounding temperature (C)
!    sdewp    Sounding dew-point (C)
!    sdrct    Sounding wind direction (degrees)
!    ssped    Sounding wind speed (m/s)
!    somega   Sounding omega vertical velocity (Pa/s)
!    shght    Sounding height (m)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!   Temporary arrays are defined and used differently by each
!   subroutine.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nzsoil
  INTEGER :: nxlg,nylg

  INTEGER :: nstyps
  INTEGER, PARAMETER :: maxsnd = 2000
  REAL,    PARAMETER :: zsndmax = 20.0E03
  INTEGER, PARAMETER :: nzsnd = 201
  INTEGER, PARAMETER :: lunit = 15
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)         ! The x-coord. of the physical and
                                         ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)         ! The y-coord. of the physical and
                                         ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)         ! The z-coord. of the computational grid.
                                         ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:)     ! The physical height coordinate defined at
                                         ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)     ! The physical height coordinate defined at
                                         ! w-point of the soil grid.

  REAL, ALLOCATABLE :: uprt   (:,:,:)    ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:)    ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:)    ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:)    ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:)    ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:)    ! Perturbation water vapor specific
                                         ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar (:,:,:,:)


  REAL, ALLOCATABLE :: tke   (:,:,:)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh   (:,:,:)     ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv   (:,:,:)     ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )

  REAL, ALLOCATABLE :: ubar   (:,:,:)    ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:)    ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:)    ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:)    ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:)    ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:)    ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:)    ! Base state water vapor specific
                                         ! humidity (kg/kg)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct(:,:,:) ! Soil fraction
  INTEGER, ALLOCATABLE :: vegtyp  (:,:)   ! Vegetation type
  REAL, ALLOCATABLE :: lai     (:,:)      ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns  (:,:)      ! Surface roughness
  REAL, ALLOCATABLE :: veg     (:,:)      ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil (:,:,:,:)    ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil (:,:,:,:)    ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:)     ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)      ! Snow depth (m)

  REAL, ALLOCATABLE :: raing(:,:)         ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)         ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                                          ! prcrate(1,1,1) = total precip. rate
                                          ! prcrate(1,1,2) = grid scale precip. rate
                                          ! prcrate(1,1,3) = cumulus precip. rate
                                          ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)      ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)        ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)        ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)      ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)       ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s)
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xs(:)      ! x location of scalar points
  REAL, ALLOCATABLE :: ys(:)      ! y location of scalar points
!
!-----------------------------------------------------------------------
!
!  Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: su(:,:)
  REAL, ALLOCATABLE :: sv(:,:)
  REAL, ALLOCATABLE :: sw(:,:)
  REAL, ALLOCATABLE :: stheta(:,:)
  REAL, ALLOCATABLE :: sqv(:,:)
  REAL, ALLOCATABLE :: sqc(:,:)
  REAL, ALLOCATABLE :: sqi(:,:)
  REAL, ALLOCATABLE :: stke(:,:)
  REAL, ALLOCATABLE :: spres(:,:)
  REAL, ALLOCATABLE :: stemp(:,:)
  REAL, ALLOCATABLE :: sdewp(:,:)
  REAL, ALLOCATABLE :: sdrct(:,:)
  REAL, ALLOCATABLE :: ssped(:,:)
  REAL, ALLOCATABLE :: somega(:,:)
  REAL, ALLOCATABLE :: shght(:,:)
!
!-----------------------------------------------------------------------
!
!  Radar sounding variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: usnd(:)
  REAL, ALLOCATABLE :: vsnd(:)
  REAL, ALLOCATABLE :: zsnd(:)
  REAL, ALLOCATABLE :: rfrctsnd(:)
  REAL, ALLOCATABLE :: ktsnd(:)

  REAL, ALLOCATABLE :: zps  (:,:,:)
  REAL, ALLOCATABLE :: usc  (:,:,:)
  REAL, ALLOCATABLE :: vsc  (:,:,:)
  REAL, ALLOCATABLE :: rfrct(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: dxfld(:)
  REAL, ALLOCATABLE :: dyfld(:)
  REAL, ALLOCATABLE :: rdxfld(:)
  REAL, ALLOCATABLE :: rdyfld(:)
  REAL, ALLOCATABLE :: slopey(:,:,:)
  REAL, ALLOCATABLE :: alphay(:,:,:)
  REAL, ALLOCATABLE :: betay(:,:,:)
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)

  REAL                :: xpt(maxsnd),ypt(maxsnd)
  INTEGER             :: ipt(maxsnd),jpt(maxsnd)
  INTEGER             :: is_good(maxsnd)
!
!-----------------------------------------------------------------------
!
!  Sounding indentification variables
!  These are required for the proper operation of the
!  plotting programs on the metgem computer at the
!  University of Oklahoma.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: stid(maxsnd)
  REAL              :: slat(maxsnd)
  REAL              :: slon(maxsnd)
  REAL              :: selev(maxsnd)
  REAL              :: smelev(maxsnd)
  INTEGER           :: istnm(maxsnd)
!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xmap(:)
  REAL, ALLOCATABLE :: ymap(:)
  REAL, ALLOCATABLE :: latgr(:,:)
  REAL, ALLOCATABLE :: longr(:,:)

  INTEGER :: i4time,iyr,imo,iday,ihr,imin,isec
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  REAL    :: time,xctr,yctr,xll,yll,xsndmap,ysndmap,xsnd,ysnd
  REAL    :: ustorm,vstorm,dzsnd
  REAL    :: xmin, xmax, ymin, ymax
  INTEGER :: i,j,k,locopt,isnd,nsnd
  INTEGER :: ireturn,hinfmt,lengbf,lenfil,nchin
  INTEGER :: scondition
  REAL    :: valuethres
  REAL    :: radradius,maxradius,ftimhr
  INTEGER :: outfmt
  INTEGER :: nhisfile_max,nhisfile,nfile
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: ofilehead
  CHARACTER (LEN=80)  :: sndtime,snddate,sndlocation
  CHARACTER (LEN=256) :: oldfile(maxsnd)
  CHARACTER (LEN=256) :: fname
  INTEGER :: oldfmt

  INTEGER :: readsplit_in

  NAMELIST /message_passing/ nproc_x,nproc_y,readsplit_in

  NAMELIST /input/ ustorm,vstorm,locopt,scondition,valuethres,          &
                   nsnd,slat,slon,selev,xpt,ypt,stid,istnm,oldfmt,oldfile

  NAMELIST /output/ dirname, ofilehead, outfmt, maxradius

  INTEGER :: nlevs, istatus
  INTEGER :: nsndx, nsndy
  INTEGER :: xptint, yptint
  INTEGER :: ijpt, i1pt,i2pt,i3pt
  INTEGER :: ilg, jlg
  INTEGER :: FIRST_SND, LAST_SND
  INTEGER :: countsnd

  CHARACTER(LEN=22), PARAMETER :: sconstr(2) =                          &
         (/ 'Vert. Integ Condensate', 'Composite Reflectivity' /)

  CHARACTER(LEN=256) :: ofile, tmpstr
  LOGICAL :: GEMPAKSND, ARPSSND, RAOBSND, WFIPSND
  REAL    :: zhght

  INTEGER :: unum         ! unit number for reading in namelist
  LOGICAL :: fexist

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL mpinit_proc(0)

  IF(myproc == 0) THEN
    WRITE(6,'(/8(/1x,a)/)')                                             &
  '###################################################################',&
  '###################################################################',&
  '###                                                             ###',&
  '### Welcome to EXTSND, a program that reads in a history file   ###',&
  '### generated by ARPS and extracts a sounding at an x-y point.  ###',&
  '###                                                             ###',&
  '###################################################################',&
  '###################################################################'

    unum = COMMAND_ARGUMENT_COUNT()
    IF (unum > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, tmpstr, lenfil, istatus )
      IF ( tmpstr(1:1) == ' ' .OR. istatus /= 0 ) THEN  ! Use standard input to be backward-compatible
        unum = 5
      ELSE
        INQUIRE(FILE=TRIM(tmpstr),EXIST=fexist)
        IF (.NOT. fexist) THEN
          WRITE(6,'(1x,3a)') 'WARNING: namelist file - ',               &
                TRIM(tmpstr),' does not exist. Falling back to standard input.'
          unum = 5
        ELSE
          unum = 0
        END IF
      END IF
    ELSE
      unum = 5
    END IF

    IF (unum /= 5) THEN
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(tmpstr),STATUS='OLD',FORM='FORMATTED')
      WRITE(*,'(1x,3a,/,1x,a,/)') 'Reading ARPSEXTSND namelist from file - ', &
              TRIM(tmpstr),' ... ','========================================'
    ELSE
      WRITE(*,'(2(1x,a,/))') 'Waiting namelist from standard input ... ', &
                             '========================================'
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in message passing options.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ(unum,message_passing)
    WRITE(6,'(1x,a)') 'Namelist block message_passing sucessfully read.'
  END IF
  readsplit(:) = readsplit_in
  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(readsplit,FINDX_NUM)

  IF (mp_opt == 0 ) THEN
    nproc_x = 1
    nproc_y = 1
    readsplit(:) = 0
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize message passing variables.
!
!-----------------------------------------------------------------------
!
  CALL mpinit_var
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    CALL get_input_file_names(unum,hinfmt,grdbasfn,hisfile,nhisfile)
    lengbf = LEN_trim(grdbasfn)
    WRITE(6,'(1x,a)') 'Successfully read in ARPS history file names.'

  END IF
  CALL mpupdatec(grdbasfn,256)
  CALL mpupdatec(hisfile,256*nhisfile_max)
  CALL mpupdatei(nhisfile,1)
  CALL mpupdatei(hinfmt,1)
  CALL mpupdatei(lengbf,1)

  xpt(:) = -999.0
  ypt(:) = -999.0
  slat(:) = -999.0
  slon(:) = -999.0
  selev(:) = -999.0
  oldfmt = 0
  oldfile(:) = ' '

  IF (myproc == 0 ) THEN
    READ(unum,input)
    WRITE(6,'(1x,a)') 'Successfully read in namelist block input.'

    IF( nsnd > maxsnd )  then
      WRITE(6,'(a,/a,i5)')                                                &
      'The number of sounding locations to be extracted exceeded maximum ',&
      'allowed. nsnd is reset to ', maxsnd
      nsnd = maxsnd
    ENDIF

  ENDIF
  CALL mpupdater(ustorm,1)
  CALL mpupdater(vstorm,1)
  CALL mpupdatei(locopt,1)
  CALL mpupdatei(scondition,1)
  CALL mpupdater(valuethres,1)
  CALL mpupdatei(nsnd,1)
  CALL mpupdater(slat,maxsnd)
  CALL mpupdater(slon,maxsnd)
  CALL mpupdater(selev,maxsnd)
  CALL mpupdater(xpt,maxsnd)
  CALL mpupdater(ypt,maxsnd)
  CALL mpupdatec(stid,8*maxsnd)
  CALL mpupdatei(istnm,maxsnd)
  CALL mpupdatei(oldfmt,1)
  IF (oldfmt == 1) CALL mpupdatec(oldfile,256*maxsnd)

  dirname   = './'
  outfmt    = 0
  ofilehead = 'SNLIST'
  maxradius = 230.0E03
  IF (myproc == 0) THEN
    READ(unum,output)
    WRITE(6,'(1x,a)') 'Successfully read in namelist block output.'

     lenfil = LEN_TRIM(dirname)
     IF(lenfil > 0) THEN
       IF(dirname(lenfil:lenfil) /= '/') THEN
         dirname(lenfil+1:lenfil+1) = '/'
         lenfil = lenfil + 1
       END IF
     ELSE
       dirname = './'
     END IF
  END IF
  radradius=maxradius/1.4142
  IF (oldfmt == 1 .AND. outfmt .NE. 5) outfmt = 0
  CALL mpupdatec(dirname,256)
  CALL mpupdatei(outfmt,1)
  CALL mpupdatec(ofilehead,256)
  CALL mpupdater(radradius,1)
  CALL mpupdater(maxradius,1)

  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE( unum )
    CALL retunit( unum )
  END IF

  GEMPAKSND = .FALSE.
  ARPSSND   = .FALSE.
  RAOBSND   = .FALSE.
  WFIPSND   = .FALSE.
  IF (outfmt == 0) THEN
    GEMPAKSND = .TRUE.
    ARPSSND   = .TRUE.
  ELSE IF (outfmt == 1) THEN
    GEMPAKSND = .TRUE.
  ELSE IF (outfmt == 2) THEN
    ARPSSND   = .TRUE.
  ELSE IF (outfmt == 3) THEN
    RAOBSND   = .TRUE.
  ELSE IF (outfmt == 4) THEN
    ALLOCATE(zsnd(nzsnd))
    ALLOCATE(usnd(nzsnd))
    ALLOCATE(vsnd(nzsnd))
    ALLOCATE(rfrctsnd(nzsnd))
    ALLOCATE(ktsnd(nzsnd))
  ELSE IF (outfmt == 5) THEN
    GEMPAKSND = .TRUE.
    ARPSSND   = .TRUE.
    WFIPSND = .TRUE.
  END IF
!
!-----------------------------------------------------------------------
!
!  Obtain the grid dimensions from input data.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0 .AND. readsplit(FINDX_H) == 0) THEN
    tmpstr = grdbasfn
    CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                          &
                   0,0,1,lvldbg,grdbasfn,istatus)
    lengbf = LEN_TRIM(grdbasfn)

    DO nfile = 1,nhisfile
      tmpstr = hisfile(nfile)
      CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                        &
                     0,0,1,lvldbg,hisfile(nfile),istatus)
    END DO
  ENDIF

  IF (myproc == 0) THEN

    CALL get_dims_from_data(hinfmt,hisfile(1:lenfil),                   &
                            nx,ny,nz,nzsoil,nstyps, ireturn)

    IF (mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN
      !
      ! Fiddle with nx/ny, which apparently are wrong.
      !
      nx = (nx - 3) / nproc_x + 3
      ny = (ny - 3) / nproc_y + 3
    END IF

    IF (nstyps <= 0) nstyps = 1
    nstyp = nstyps ! Copy to global variable

    IF( ireturn /= 0 ) THEN
      PRINT*,'Problem occured when trying to get dimensions from data.'
      PRINT*,'Program stopped.'
      CALL arpsstop('Problem with data.',1)
    END IF
  END IF
  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)
  CALL mpupdatei(nzsoil,1)
  CALL mpupdatei(nstyps,1)
  CALL mpupdatei(nstyp,1)

  CALL mpupdatei(nscalar, 1)
  CALL mpupdatei(nscalarq,1)
  CALL mpupdatei(P_QC,1);
  CALL mpupdatei(P_QR,1)
  CALL mpupdatei(P_QI,1)
  CALL mpupdatei(P_QS,1)
  CALL mpupdatei(P_QG,1)
  CALL mpupdatei(P_QH,1)
  CALL mpupdatei(P_NC,1)
  CALL mpupdatei(P_NR,1)
  CALL mpupdatei(P_NI,1)
  CALL mpupdatei(P_NS,1)
  CALL mpupdatei(P_NG,1)
  CALL mpupdatei(P_NH,1)
  CALL mpupdatei(P_ZR,1)
  CALL mpupdatei(P_ZI,1)
  CALL mpupdatei(P_ZS,1)
  CALL mpupdatei(P_ZG,1)
  CALL mpupdatei(P_ZH,1)

  CALL mpupdatec(qnames,nscalar*40)
  CALL mpupdatec(qdescp,nscalar*40)


  IF(outfmt == 4) THEN
    ALLOCATE(zps(nx,ny,nz))
    ALLOCATE(usc(nx,ny,nz))
    ALLOCATE(vsc(nx,ny,nz))
    ALLOCATE(rfrct(nx,ny,nz))
  END IF

  IF (myproc == 0 )                                                     &
    WRITE(6,'(1x,4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz,', nzsoil=',nzsoil

  IF (locopt == 3) THEN
    xptint = INT(xpt(1))
    yptint = INT(ypt(1))

    nsndx = 0
    DO i = 2,nx-2
      ilg = (nx-3)*(loc_x-1) + i
      IF ( MOD((ilg-2),xptint) == 0 ) THEN
        nsndx = nsndx + 1
        ipt(nsndx) = i
      END IF
    END DO

    nsndy = 0
    DO j = 2,ny-2
      jlg = (ny-3)*(loc_y-1)+j
      IF ( MOD((jlg-2),yptint) == 0 ) THEN
        nsndy = nsndy + 1
        jpt((nsndy-1)*nsndx+1) = j
      END IF
    END DO

    nsnd  = nsndx*nsndy
    IF( nsnd > maxsnd )  then
      WRITE(6,'(/1x,a,/2(a,i5)/)')                                         &
      'The number of sounding locations to be extracted exceeded maximum ',&
      ' allowed. nsnd = ', nsnd,' maxsnd = ',maxsnd
      CALL arpsstop('Please increase maxsnd in src/arpsextsnd/arpsextsnd.f90.',1)
    ENDIF

    jpt(2:nsndx) = jpt(1)
    DO j = 1,nsndy
      ipt((j-1)*nsndx+1) = ipt(1)
    END DO

    DO j = 2,nsndy
      DO i = 2,nsndx
        isnd = (j-1)*nsndx + i
        ipt(isnd) = ipt(i)
        jpt(isnd) = jpt((j-1)*nsndx+1)
      END DO
    END DO

  END IF

  !WRITE(0,'(/,a,I2.2,3(a,I4),6(a,I3),/)')     &
  !  'Rank = ',myproc,': nsnd = ',nsndx,' x ',nsndy,' = ',nsnd,  &
  ! ', (',ipt(1),', ',jpt(1),') -- (',ipt(nsnd),',',jpt(nsnd),   &
  ! '). xptint,yptint = ',xptint,', ',yptint
!
!-----------------------------------------------------------------------
!
! Allocate arrays
!
!-----------------------------------------------------------------------
!
  ALLOCATE(x      (nx),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:x')
  ALLOCATE(y      (ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:y')
  ALLOCATE(z      (nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:z')
  ALLOCATE(zp     (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:zp')
  ALLOCATE(zpsoil  (nx,ny,nzsoil),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:zpsoil')
  ALLOCATE(uprt   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:uprt')
  ALLOCATE(vprt   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:vprt')
  ALLOCATE(wprt   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:wprt')
  ALLOCATE(ptprt  (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ptprt')
  ALLOCATE(pprt   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:pprt')
  ALLOCATE(qvprt  (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:qvprt')
  ALLOCATE(qscalar(nx,ny,nz,nscalar),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:qscalar')
  ALLOCATE(tke    (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:tke')
  ALLOCATE(kmh    (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:kmh')
  ALLOCATE(kmv    (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:kmv')
  ALLOCATE(ubar   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ubar')
  ALLOCATE(vbar   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:vbar')
  ALLOCATE(wbar   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:wbar')
  ALLOCATE(ptbar  (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ptbar')
  ALLOCATE(pbar   (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:pbar')
  ALLOCATE(rhobar (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:rhobar')
  ALLOCATE(qvbar  (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:qvbar')

  ALLOCATE(soiltyp (nx,ny,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:soiltyp')
  ALLOCATE(stypfrct(nx,ny,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:stypfrct')
  ALLOCATE(vegtyp(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:vegtyp')
  ALLOCATE(lai    (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:lai')
  ALLOCATE(roufns (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:roufns')
  ALLOCATE(veg    (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:veg')
  ALLOCATE(tsoil  (nx,ny,nzsoil,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:tsoil')
  ALLOCATE(qsoil  (nx,ny,nzsoil,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:qsoil')
  ALLOCATE(wetcanp(nx,ny,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:wetcanp')
  ALLOCATE(snowdpth(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:snowdpth')
  ALLOCATE(raing(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:raing')
  ALLOCATE(rainc(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:rainc')
  ALLOCATE(prcrate(nx,ny,4),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:prcrate')
  ALLOCATE(radfrc(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:radfrc')
  ALLOCATE(radsw (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:radsw')
  ALLOCATE(rnflx (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:rnflx')
  ALLOCATE(radswnet(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:radswnet')
  ALLOCATE(radlwin(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:radlwin')
  ALLOCATE(usflx (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:usflx')
  ALLOCATE(vsflx (nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:vsflx')
  ALLOCATE(ptsflx(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ptsflx')
  ALLOCATE(qvsflx(nx,ny),stat=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:qvsflx')

  ALLOCATE(xs(nx),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:xs')
  ALLOCATE(ys(ny),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ys')
  ALLOCATE(su(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:su')
  ALLOCATE(sv(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sv')
  ALLOCATE(sw(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sw')
  ALLOCATE(stheta(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:stheta')
  ALLOCATE(sqv(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sqv')
  ALLOCATE(sqc(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sqc')
  ALLOCATE(sqi(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sqi')
  ALLOCATE(stke(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:stke')
  ALLOCATE(spres(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:spres')
  ALLOCATE(stemp(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:stemp')
  ALLOCATE(sdewp(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sdewp')
  ALLOCATE(sdrct(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:sdrct')
  ALLOCATE(ssped(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ssped')
  ALLOCATE(somega(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:somega')
  ALLOCATE(shght(nz,maxsnd),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:shght')
  ALLOCATE(dxfld(nx),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:dxfld')
  ALLOCATE(dyfld(ny),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:dyfld')
  ALLOCATE(rdxfld(nx),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:rdxfld')
  ALLOCATE(rdyfld(ny),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:rdyfld')
  ALLOCATE(slopey(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:slopey')
  ALLOCATE(alphay(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:alphay')
  ALLOCATE(betay(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:betay')
  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:tem1')
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:tem2')
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:tem3')

  ALLOCATE(xmap(nx),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:xmap')
  ALLOCATE(ymap(ny),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:ymap')
  ALLOCATE(latgr(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:latgr')
  ALLOCATE(longr(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, 'arpsextsnd:longr')

  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0
  zpsoil =0.0
  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
  qscalar =0.0
  tke    =0.0
  kmh    =0.0
  kmv    =0.0
  ubar   =0.0
  vbar   =0.0
  wbar   =0.0
  ptbar  =0.0
  pbar   =0.0
  radsw =0.0
  qvbar  =0.0

  soiltyp =0
  stypfrct=0.0
  vegtyp=0.0
  lai    =0.0
  roufns =0.0
  veg    =0.0
  tsoil  =0.0
  qsoil  =0.0
  wetcanp=0.0
  snowdpth=0.0
  raing=0.0
  rainc=0.0
  prcrate=0.0
  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet =0.0
  radlwin =0.0
  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  xs=0.0
  ys=0.0
  su=0.0
  sv=0.0
  sw=0.0
  stheta=0.0
  sqv=0.0
  sqc=0.0
  sqi=0.0
  stke=0.0
  spres=0.0
  stemp=0.0
  sdewp=0.0
  sdrct=0.0
  ssped=0.0
  somega=0.0
  shght=0.0
  dxfld=0.0
  dyfld=0.0
  rdxfld=0.0
  rdyfld=0.0
  slopey=0.0
  alphay=0.0
  betay=0.0
  tem1=0.0
  tem2=0.0
  tem3=0.0
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
  DO nfile = 1, nhisfile

    lenfil=len_trim(hisfile(nfile))
    WRITE(6,'(/a,a,a)')                                                 &
        ' Data set ', trim(hisfile(nfile)),' to be read.'

  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
               hinfmt,nchin,grdbasfn(1:lengbf),lengbf,                  &
               hisfile(nfile)(1:lenfil),lenfil,time,                    &
               x,y,z,zp,zpsoil, uprt ,vprt ,wprt ,ptprt, pprt ,         &
               qvprt, qscalar, tke,kmh,kmv,                  &
               ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,            &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
  IF( ireturn == 0 ) THEN   ! successful read

!
!-----------------------------------------------------------------------
!
!  Calculate scalar locations
!
!-----------------------------------------------------------------------
!
    DO i=1,nx-1
      xs(i)=0.5*(x(i)+x(i+1))
    END DO
    DO j=1,ny-1
      ys(j)=0.5*(y(j)+y(j+1))
    END DO
    dx=x(2)-x(1)
    dy=y(2)-y(1)
!
!-----------------------------------------------------------------------
!
!  Compute sounding initial and forecast time and write them
!  into the file name
!
!-----------------------------------------------------------------------
!
    CALL ctim2abss( year,month,day,hour,minute,second,i4time)
    i4time=i4time + nint(time)
    ftimhr=time/3600.
    CALL abss2ctim( i4time, iyr, imo, iday, ihr, imin, isec )
    WRITE(6,'(a,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)')         &
        '  Time of history data: ',                                     &
        iyr,'/',imo,'/',iday,':',ihr,':',imin,':',isec
    WRITE(ofile,'(a,i4.4,2i2.2,a,3i2.2)')                               &
           TRIM(ofilehead),iyr,imo,iday,'_',ihr,imin,isec
!
!-----------------------------------------------------------------------
!
!  Set up the grid locations in grid space and map space.
!
!-----------------------------------------------------------------------
!
    latnot(1)=trulat1
    latnot(2)=trulat2
    CALL setmapr(mapproj,sclfct,latnot,trulon)
    CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
    PRINT *, ' dx= ',dx,' dy= ',dy
    nxlg = (nx-3)*nproc_x+3
    nylg = (ny-3)*nproc_y+3
    xll=xctr-(0.5*(nxlg-3)*dx)
    yll=yctr-(0.5*(nylg-3)*dy)

    DO i=1,nx-1
      xmap(i)=xll+xs(i)
    END DO
    xmap(nx)=2.*xmap(nx-1)-xmap(nx-2)
    DO j=1,ny-1
      ymap(j)=yll+ys(j)
    END DO
    ymap(ny)=2.*ymap(ny-1)-ymap(ny-2)
    CALL xytoll(nx,ny,xmap,ymap,latgr,longr)

    xmin = xs(2)
    xmax = xs(nx-1)
    ymin = ys(2)
    ymax = ys(ny-1)

    print *, ' xmin,xmax: ',xmin,xmax
    print *, ' ymin,ymax: ',ymin,ymax

    IF(outfmt == 4) THEN

      dzsnd=zsndmax/(nzsnd-1)
      DO k=1,nzsnd
        zsnd(k)=(k-1)*dzsnd
      END DO

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
          END DO
        END DO
      END DO

      DO j=1,ny-1
        DO i=1,nx-1
          zps(i,j,nz)=2.0*zp(i,j,nz-1)-zp(i,j,nz-2)
        END DO
      END DO

!-----------------------------------------------------------------------
!
!  Bring u and v to a common grid, the scalar points.
!
!-----------------------------------------------------------------------

      uprt(:,:,:)=ubar(:,:,:)+uprt(:,:,:)
      vprt(:,:,:)=vbar(:,:,:)+vprt(:,:,:)
      CALL avgx(uprt, 0, nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, usc)
      CALL avgy(vprt, 0, nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, vsc)
!
      qvprt(:,:,:)=qvbar(:,:,:)+qvprt(:,:,:)
      CALL refract(nx,ny,nz,ptprt,pprt,qvprt,ptbar,pbar,rfrct)
!
      DO isnd=1,nsnd
        CALL lltoxy(1,1,slat(isnd),slon(isnd),xsnd,ysnd)
        xsnd=xsnd-xll
        ysnd=ysnd-yll
          print *, ' stid: ',trim(stid(isnd)),                        &
                 ' xsnd=',xsnd,'  ysnd=',ysnd
          print *, ' slat: ',slat(isnd),'  slon: ',slon(isnd)
        IF(xsnd >= xmin .AND. xsnd <= xmax .AND.                      &
           ysnd >= ymin .AND. ysnd <= ymax ) THEN
          print *, ' sounding is in grid'
          CALL extenvprf3(nx,ny,nz,nzsnd,                             &
                    xs,ys,zps,usc,vsc,rfrct,                          &
                    xsnd,ysnd,selev(isnd),radradius,maxradius,        &
                    zsnd,ktsnd,usnd,vsnd,rfrctsnd,istatus)
!
          OPEN(3,FILE=trim(dirname)//trim(ofile)//'_'//trim(stid(isnd))//'.radsnd',STATUS='unknown')
          WRITE(3,'(a,a8,3x,a,i4.4,i2.2,i2.2,3(a1,i2.2))')             &
            ' STID: ',stid(isnd),                                      &
            'TIME: ',iyr,imo,iday,' ',ihr,':',imin,':',isec
          WRITE(3,'(a,f10.4,3x,a,f10.4,3x,a,f7.1/)')                   &
            ' SLAT:',slat(isnd),'SLON:',slon(isnd),                    &
            'SELV:',selev(isnd)
          WRITE(3,'(5x,a)') 'zsnd     usnd     vsnd   rfrctsnd'
!
!-----------------------------------------------------------------------
!
!    Print out of sounding
!
!-----------------------------------------------------------------------
!
          DO k=1,nzsnd
            WRITE(3,'(1x,4f9.2)')                                      &
              zsnd(k),usnd(k),vsnd(k),rfrctsnd(k)
          END DO

          CLOSE(3)

        END IF  ! inside grid

      END DO  ! isnd loop

    ELSE ! outfmt /= 4

    tem3(:,:,:) = 0.0
    IF (locopt == 3 .AND. scondition == 1) THEN
      CALL cal_vic(tem3,qscalar,rhobar,zp,nx,ny,nz,nscalar,nscalarq,tem1)
    ELSE IF (locopt == 3 .AND. scondition == 2) THEN
      CALL temper (nx,ny,nz,ptbar, ptprt, pprt ,pbar,tem2)
      CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, tem2, tem1)
!      CALL reflec(nx,ny,nz, rhobar, qr, qs, qh, tem1)
      CALL cal_rfc(nx, ny, nz, tem1, tem3)
    END IF
!
!   Set the range for checking the soundings.  Make sure that more than one
!   processor doesn't check a point, as each processor writes out its own
!   files!
!
    FIRST_SND = nsnd
    LAST_SND  = 1
    countsnd  = 0

    WRITE(6,*)
    DO isnd=1,nsnd
      IF(locopt == 1) THEN

        IF(slat(isnd) < -90.) EXIT

        CALL lltoxy(1,1,slat(isnd),slon(isnd),xpt(isnd),ypt(isnd))
        xpt(isnd)=xpt(isnd)-xll
        ypt(isnd)=ypt(isnd)-yll

      ELSE IF (locopt == 2) THEN
        IF (nfile == 1) THEN
          xpt(isnd)=xpt(isnd)*1000.
          ypt(isnd)=ypt(isnd)*1000.
        END IF
        xsndmap=xpt(isnd)+xll
        ysndmap=ypt(isnd)+yll
        CALL xytoll(1,1,xsndmap,ysndmap,slat(isnd),slon(isnd))
      ELSE IF (locopt == 3) THEN
        xpt(isnd) = xs(ipt(isnd))
        ypt(isnd) = ys(jpt(isnd))
        xsndmap=xpt(isnd)+xll
        ysndmap=ypt(isnd)+yll
        CALL xytoll(1,1,xsndmap,ysndmap,slat(isnd),slon(isnd))
        istnm(isnd) = jpt(isnd)*10000+ipt(isnd)

        ilg = (nx-3)*(loc_x-1) + ipt(isnd)
        jlg = (ny-3)*(loc_y-1) + jpt(isnd)
        ijpt = ilg-2 + (jlg-2)*(nxlg-3)
        i1pt = mod(ijpt,26)         + ICHAR('A')
        i2pt = mod(ijpt/26,26)      + ICHAR('A')
        i3pt = mod(ijpt/(26*26),26) + ICHAR('A')
        WRITE(stid(isnd),'(3a)') CHAR(i3pt),CHAR(i2pt),CHAR(i1pt)

      ELSE
        WRITE(6,'(/,a,i2,/)') 'ERROR: unknown option locopt = ',locopt
        CALL arpsstop('Please check parameter locopt in the namelist input.',1)
      END IF

!  Check to see if the sounding is out of the grid box.  If so, computed
!  data will be garbage!
!
      is_good(isnd)=1
      IF( xpt(isnd) < xmin .or. xpt(isnd) > xmax .or.                   &
          ypt(isnd) < ymin .or. ypt(isnd) > ymax ) THEN
!       WRITE(6,'(/,2x,a,a6,a,/)') 'Station ',stid(isnd),' is outside of the grid.'
        is_good(isnd)=0
        CYCLE
      END IF

!-----------------------------------------------------------------------
!
!  Find location in ARPS grid.
!
!-----------------------------------------------------------------------
!
      !
      ! It is requried that valuethres must > 0 and
      ! tem3(:,:,:) = 0 when scondition = 0
      !
      IF (locopt == 3) THEN
        xsndmap = tem3(ipt(isnd),jpt(isnd),1)
        IF (xsndmap > valuethres) THEN
          WRITE(6,'(3(a,I6),4a,G10.2)') '- isnd =',isnd,                &
                  '  at (',ipt(isnd),',',jpt(isnd),') was skipped ',    &
                  'because ',sconstr(scondition),' = ',xsndmap
          is_good(isnd) = 0
        END IF

      ELSE

        CALL setijloc(1,1,nx,ny,xpt(isnd),ypt(isnd),                    &
                      xs,ys,ipt(isnd),jpt(isnd))
      END IF

      IF (is_good(isnd) == 1) THEN

        countsnd = countsnd + 1

        WRITE(6,'(a,i6,a,I8.8,a,a,2(a,F7.2),a)')                        &
         '+ isnd =',isnd,', stnm = ',istnm(isnd),', stid = ',           &
         TRIM(stid(isnd)),', xpt = ',(xpt(isnd)*0.001),                 &
         ' km, ypt = ',(ypt(isnd)*0.001),' km'

        IF (locopt == 3) THEN

        WRITE (6,'(16x,2(a,i6),2(a,f7.2))')                             &
            'at (',ipt(isnd),',',jpt(isnd),                             &
            '),          lat = ',slat(isnd),',    lon = ',slon(isnd)

        ELSE

        WRITE (6,'(2x,4(a,f12.2),/2X,a,i6,a,i6,a)')                     &
            '   loc x= ',(0.001*xpt(isnd)),', y= ',(0.001*ypt(isnd)),   &
            ', lat= ',slat(isnd),', lon= ',slon(isnd),                  &
            '       found near point (',ipt(isnd),',',jpt(isnd),')'

        WRITE (6,'(12x,4(a,f12.2))')            &
            '      x= ',(0.001*xs(ipt(isnd))),                          &
            ', y= ',(0.001*ys(jpt(isnd))),                              &
            ', lat= ',latgr(ipt(isnd),jpt(isnd)),                       &
            ', lon= ',longr(ipt(isnd),jpt(isnd))

        END IF

        FIRST_SND = MIN(FIRST_SND,isnd)
        LAST_SND  = MAX(LAST_SND,isnd)
      ELSE

        WRITE(6,'(a,i6,a,I8.8,a,a,2(a,F7.2),a)')                        &
         'XXX Problem isnd =',isnd,', stnm = ',istnm(isnd),', stid = ', &
         TRIM(stid(isnd)),', xpt = ',(xpt(isnd)*0.001),                 &
         ' km, ypt = ',(ypt(isnd)*0.001),' km'
      END IF

    END DO
    WRITE(6,'(/,1x,a,I6,/)') 'Total valid soundings = ',countsnd
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
    IF (locopt == 3) THEN
      CALL coextract(nx,ny,nz,nsnd,zp,ipt,jpt,                          &
                  uprt, vprt, wprt, ptprt, pprt, qvprt,                 &
                  ubar, vbar, wbar, ptbar, pbar, qvbar, is_good,        &
                  su,sv,sw,stheta,spres,shght,sqv,                      &
                  smelev,nlevs)
    ELSE

      CALL setdxdy(nx,ny,                                               &
                   1,nx-1,1,ny-1,                                       &
                   xs,ys,dxfld,dyfld,rdxfld,rdyfld)
      CALL colint(nx,ny,nz,maxsnd,nsnd,                                 &
                  xs,ys,zp,xpt,ypt,ipt,jpt,                             &
                  uprt, vprt, wprt, ptprt, pprt, qvprt,                 &
                  ubar, vbar, wbar, ptbar, pbar, qvbar,                 &
                  tke,qscalar(:,:,:,P_QC),qscalar(:,:,:,P_QI),is_good,  &
                  su,sv,sw,stheta,spres,shght,sqv,sqc,sqi,stke,smelev,  &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  slopey,alphay,betay,                                  &
                  tem1,nlevs)
    END IF
!
!-----------------------------------------------------------------------
!
!  Add back storm motion that ARPS subtracts from the sounding winds.
!  Force mixing ratios to be positive definite.
!
!-----------------------------------------------------------------------
!
    DO isnd=1,nsnd
      IF (is_good(isnd) == 0 ) cycle
      DO k=1,nlevs
        su(k,isnd)=su(k,isnd)+ustorm
        sv(k,isnd)=sv(k,isnd)+vstorm
        sqv(k,isnd)=max(1.0E-07,sqv(k,isnd))
        sqc(k,isnd)=max(0.0,sqc(k,isnd))
        sqi(k,isnd)=max(0.0,sqi(k,isnd))
        stke(k,isnd)=max(0.0,stke(k,isnd))
      END DO
!
!-----------------------------------------------------------------------
!
!  Convert sounding to desired units/quantities
!
!-----------------------------------------------------------------------
!
      CALL cnvsnd(su(1,isnd),sv(1,isnd),sw(1,isnd),stheta(1,isnd),      &
                spres(1,isnd),sqv(1,isnd),slon(isnd),                   &
                sdrct(1,isnd),ssped(1,isnd),somega(1,isnd),             &
                stemp(1,isnd),sdewp(1,isnd),nlevs)

      WRITE(ofile,'(a,i4.4,2i2.2,a,3i2.2)')                               &
           TRIM(ofilehead),iyr,imo,iday,'_',ihr,imin,isec
      tmpstr = ofile
      IF (locopt == 1 .OR. locopt == 2) THEN
                           ! append statid for locopt = 1/2
                           ! because data for each station is in one file
        IF (oldfmt == 1) THEN    ! Use explicit file names
          ofile=oldfile(isnd)
        ELSE
          WRITE(ofile,'(3a)') TRIM(tmpstr),'_',TRIM(stid(isnd))
        END IF
      ELSE IF (locopt == 3 .AND. mp_opt > 0) THEN
                           ! append proc. num. for mp_opt >0
                           ! because we want to distinguished file names
        CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                      &
                       0,0,0,lvldbg,ofile,istatus)
      END IF

      IF (GEMPAKSND) THEN
!
!-----------------------------------------------------------------------
!
!  Output sounding to look like GEMPAK SNLIST.FIL
!  to use the Skew-T and Hodograph programs
!  e.g., those by Bill McCaul and Rich Carpenter and NSHARP
!
!  Example of what the header looks like:
!
!23456789012345678901234567890123456789012345678901234567890
!SNPARM = PRES;HGHT;TMPC;DWPC;DRCT;SPED;OMEG
!
!STID = SEP        STNM =    72260   TIME = 920308/1500
!SLAT =    32.21   SLON =   -98.18   SELEV =   399.0
!
!  PRES     HGHT     TMPC     DWPC     DRCT     SPED     OMEG
!
!-----------------------------------------------------------------------
!
        IF(locopt /= 3 .OR. isnd ==  FIRST_SND) THEN
          IF (oldfmt == 1) THEN
            OPEN(3,file=trim(ofile),STATUS='unknown')
          ELSE
            OPEN(3,FILE=trim(dirname)//trim(ofile)//'.sounding',STATUS='unknown')
          END IF
        END IF

        WRITE(3,'(/a/)') ' SNPARM = PRES;HGHT;TMPC;DWPC;DRCT;SPED;OMEG'
!        iyr=MOD(iyr,100)
!wdt update
        WRITE(3,'(a,a8,3x,a,i8,3x,a,i2.2,i2.2,i2.2,a1,i2.2,i2.2)')        &
            ' STID = ',stid(isnd),                                        &
            'STNM = ',istnm(isnd),                                        &
            'TIME = ',MOD(iyr,100),imo,iday,'/',ihr,imin
        WRITE(3,'(a,f8.2,3x,a,f8.2,3x,a,f7.1/)')                          &
            ' SLAT = ',slat(isnd),'SLON = ',slon(isnd),                   &
            'SELV = ',selev(isnd)
        WRITE(3,'(6x,a)')                                                 &
            'PRES     HGHT     TMPC     DWPC     DRCT     SPED     OMEG'
!
!-----------------------------------------------------------------------
!
!    Print out of sounding
!
!-----------------------------------------------------------------------
!
        DO k=2,nlevs
          WRITE(3,'(1x,7(F9.2))')                                       &
                spres(k,isnd),shght(k,isnd),                            &
                stemp(k,isnd),sdewp(k,isnd),                            &
                sdrct(k,isnd),ssped(k,isnd),somega(k,isnd)
        END DO

        IF (locopt /= 3 .OR. isnd == LAST_SND) THEN
          CLOSE(3)
          WRITE(6,'(1x,3a,/)') 'Sounding file ',trim(ofile)//'.sounding', &
                             ' was produced.'
        END IF

     END IF

!-----------------------------------------------------------------------
!
!  OUTPUT sounding according to follow format which can be used in
!         ARPS as base field like May20.snd
!
!              Sounding File Format for ARPS 3.0
!
!  Record 1: a header line, such as "1-D Sounding Input for ARPS"
!            (skipped)
!  Record 2: miscellaneous description of sounding (skipped)
!  Record 3: time of sounding (character*72)
!  Record 4: date of sounding (character*72)
!  Record 5: location of sounding (character*72)
!  Record 6: three character strings designate the sounding
!            data type, e.g.
!            'pressure' 'potential_temperature' 'relative_humidity'.
!            Only the first character of the strings
!            is decoded, and thus the first character should not be
!            left blank.  Note that either upper or lower case may be
!            used.  A more detailed explanation is provided in the
!            portion of the code where these strings are declared.
! ------------------------------
!  The following three strings designate the type of sounding data.
!
!  if height(1:1)='h' or 'H', the sounding is given on height levels
!  if height(1:1)='p' or 'P', the sounding is given on pressure levels
!
!  if therm(1:1) ='t' or 'T', the sounding is specified in temperature
!  if therm(1:1) ='p' or 'P', the sounding is specified in potential
!                             temperature.
!
!  if humid(1:1) ='s' or 'S', the soundings uses specific humidity
!  if humid(1:1) ='r' or 'R', the sounding uses relative humidity,
!  if humid(1:1) ='d' or 'D', the soundings uses dewpoint temperature,
!
!  if wind(1:1) ='x' or 'Y', the sounding is specified in x and y
!                            component of velocity (u and v)
!  if wind(1:1) ='d' or 'D', the sounding is specified in direction
!                            and speed in m/s.
!  if wind(1:1) ='k' or 'K', the sounding is specified in direction
!                            and speed in knots.
!--------------------------------
!  Record 7: Ground-level height (m) and the correspoding pressure
!            (Pascal) when sounding is specified at height levels.
!            When it is given on pressure levels, this record is
!            not used. The last sounding data is assumed to be at
!            the ground level (z=0).
!  Record 8: number of data levels in sounding (variable lvlsnd )
!  Record 9: a line of data column labels (not read in)
!  Record 10 to the end:
!            sounding data in the order 'z/p, pt/t, qv/rh, u, v'
!
!  For records 10  to the end, there is one line of data
!            corresponding to each sounding level.
!
!  Important convention:
!            Line 10 corresponds to the level k = lvlsnd
!            Line 11 corresponds to the level k = (lvlsnd -1)
!            etc.
!
!  Units of the data:
!            pressure: Pascals, height: meters, temperature, dewpoint,
!            and potential temperature: degrees Kelvin, specific
!            humidity kg/kg, relative humidity: 0 to 1.
!
!  CAUTION: The sounding data have to be listed in the order of
!           decreasing height or increasing pressure.
!
!-----------------------------------------------------------------------
!
      IF (ARPSSND) THEN
!        WRITE(*,*)  year,month,day,hour,minute,second
        WRITE(sndtime,'(i2,a,i2,a,a)') ihr,':',imin,':','00 CST'
        IF(ihr <=9 ) sndtime(1:1)='0'
        IF(imin <=9 ) sndtime(4:4)='0'
        WRITE(snddate,'(i3,i3,a,i5,a,f6.1,a,f6.1,a)')  iday,imo,',',year, &
          ' - Domain speed plused (umove=',ustorm,',vmove=',ustorm,')'
        WRITE(sndlocation,'(a,f8.2,3x,a,f8.2,3x,a,f7.1)')                 &
         'SLAT=',slat(isnd),'SLON=',slon(isnd),'SELV=',smelev(isnd)

        IF (locopt /= 3 .OR. isnd == FIRST_SND)  THEN
          IF (oldfmt == 1) THEN
            OPEN(4,FILE=trim(ofile)//'forARPS',STATUS='unknown')
          ELSE
            OPEN(4,FILE=trim(dirname)//trim(ofile)//'.sound',STATUS='unknown')
          ENDIF
        END IF

        WRITE(4,*) '1-D Sounding Input for ARPS'
        WRITE(4,*) 'supercell storm'
        WRITE(4,'(1x,a)')  sndtime
        WRITE(4,'(a)')  snddate
        WRITE(4,'(1x,a)')  sndlocation
        WRITE(4,*) '''height'' ''potential temperature'' ''specific humidity'' ''uv'' '
        WRITE(4,*)   shght(1,isnd), spres(1,isnd)*100
        WRITE(4,*)   nlevs
        WRITE(4,*)   ' height potential  temperature   SH    U         V'

!
        DO k=nlevs,2,-1
          WRITE(4,'(1x,6(F14.5))')                                        &
                shght(k,isnd),                                            &
                stheta(k,isnd),sqv(k,isnd),                               &
                su(k,isnd),sv(k,isnd)
        END DO
        IF (locopt /= 3 .OR. isnd == LAST_SND) THEN
          CLOSE(4)
          WRITE(6,'(1x,3a,/)')                                            &
            'Sounding file ',trim(ofile)//'.sound',' was produced.'
        END IF
      END IF

!-----------------------------------------------------------------------
      IF (RAOBSND) THEN

        IF(locopt /= 3 .OR. isnd == FIRST_SND) THEN
          fname=TRIM(dirname)//TRIM(ofile)//'.snd'
          OPEN(lunit,FILE=trim(fname),FORM='formatted',STATUS='unknown')
        END IF

        WRITE (lunit,'(i12.8,i12,f11.4,f15.4,f15.0,5x,a3)')              &
           istnm(isnd),nlevs,slat(isnd),slon(isnd),smelev(isnd),stid(isnd)
!
!-----------------------------------------------------------------------
!
!    Print out of sounding
!
!-----------------------------------------------------------------------
!
        DO k=2,nlevs

          zhght = zp(ipt(isnd),jpt(isnd),k)-zp(ipt(isnd),jpt(isnd),2)

          IF (locopt == 3 .AND. scondition == 3 .AND. zhght <= valuethres) CYCLE

          WRITE(lunit,'(1x,F12.3,5(F12.4))')                            &
                        shght(k,isnd),spres(k,isnd),                    &
                        stemp(k,isnd),sdewp(k,isnd),                    &
                        sdrct(k,isnd),ssped(k,isnd)
        END DO

        IF (locopt /= 3 .OR. isnd == LAST_SND) THEN
          CLOSE(lunit)
          WRITE(6,'(1x,3a,/)') 'Sounding file ',trim(fname),            &
                             ' was produced.'
        END IF

      END IF  ! RAOBSND

      IF (WFIPSND) THEN

        IF(locopt /= 3 .OR. isnd == FIRST_SND) THEN
          fname=' '
          fname=trim(dirname)//trim(ofile)//'.wfipsnd'
          OPEN(lunit,FILE=trim(fname),STATUS='unknown')
        END IF

        WRITE (lunit,'(f8.2,2(a,i4.4,2i2.2,a,2i2.2))')                 &
           ftimhr,' h fcst Init:',year,month,day,               &
           '_',hour,minute,'  Valid:',iyr,imo,iday,'_',ihr,imin
        WRITE (lunit,'(1x,i5.5,i8,2f10.4,2f10.1,2x,a4)')               &
           istnm(isnd),nlevs,slat(isnd),slon(isnd),selev(isnd),        &
           smelev(isnd),stid(isnd)
        WRITE (lunit,'(2a)')                                           &
          'Hgt(m ASL) Press(hPa) Temp(C) QV(g/kg) u(m/s)  ',           &
          'v(m/s)   w(m/s)  TKE(m2/s2) qc(g/kg)  qi(g/kg)'
!
!-----------------------------------------------------------------------
!
!    Print out of sounding
!    Recall data are on scalar levels, physical space k=2 to nlevs=nz-2
!
!-----------------------------------------------------------------------
!
        DO k=2,nlevs

          WRITE(lunit,'(1x,f8.1,f9.1,f8.1,f10.4,2f8.2,4f10.4)')        &
                        shght(k,isnd),spres(k,isnd),                   &
                        stemp(k,isnd),(1000.*sqv(k,isnd)),             &
                        su(k,isnd),sv(k,isnd),sw(k,isnd),stke(k,isnd), &
                        (1000.*sqc(k,isnd)),(1000.*sqi(k,isnd))
        END DO

        IF (locopt /= 3 .OR. isnd == LAST_SND) THEN
          CLOSE(lunit)
          WRITE(6,'(1x,3a,/)') 'Sounding file ',TRIM(fname),           &
                             ' was produced.'
        END IF

      END IF

    END DO  ! the number of soundings

  END IF

  ELSE

    WRITE(6,'(1x,2a,/)') 'Error reading data file ',hisfile(nfile)

  END IF   ! success reading

  END DO   ! data file loop

  IF (mp_opt > 0) CALL mpexit(0)

  STOP
END PROGRAM ARPSEXTSND
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE COLINT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE colint(nx,ny,nz,maxsnd,nsnd,                                 &
           xs,ys,zp,xpt,ypt,ipt,jpt,                                    &
           uprt, vprt, wprt, ptprt, pprt, qvprt,                        &
           ubar, vbar, wbar, ptbar, pbar, qvbar,                        &
           tke,qc,qi, is_good,                                          &
           su,sv,sw,stheta,spres,shght,sqv,sqc,sqi,stke,selev,          &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           tem1,nlevs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates ARPS history data in the horizontal to create
!  a column of data located at point xpt, ypt.
!
!  Bilinear interpolation is used.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  October, 1992 (K. Brewster)
!  Conversion to ARPS 3.0.
!
!  October, 1994 (K. Brewster)
!  Conversion to ARPS 4.0.
!
!  04/10/2005 (K. Brewster)
!  Addition of vertical velocity fields.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz Dimensions of ARPS grids.
!
!    xs       x coordinate of scalar points in physical/comp. space (m)
!    ys       y coordinate of scalar points in physical/comp. space (m)
!    zp       z coordinate of scalar grid points in physical space (m)
!
!    xpt      x coordinate of desired sounding (m)
!    ypt      y coordinate of desired sounding (m)
!
!    ipt      i index of grid point just west of xpt,ypt
!    jpt      j index of grid point just south of xpt,ypt
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    vprt     z component of perturbation velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    is_good  Make no computations if this sounding is outside the grid.
!
!  OUTPUT:
!
!    su       Interpolated u wind component.  (m/s)
!    sv       Interpolated v wind component.  (m/s)
!    sw       Interpolated vertical velocity. (m/s)
!    stheta   Interpolated potential temperature (K).
!    spres    Interpolated pressure. (Pascals)
!    shght    Interpolated height (meters)
!    sqv      Interpolated water vapor mixing ratio (kg/kg).
!    sqc      Interpolated cloud water (kg/kg).
!    sqi      Interpolated cloud ice (kg/kg).
!    stke     Interpolated TKE turbulence
!    selev    Interpolated surface elevation (m)
!    nlevs    Number of above-ground sounding levels.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
!  Arguments -- location data
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Dimensions of ARPS grids.
  INTEGER :: maxsnd
  INTEGER :: nsnd
  REAL :: xs(nx)               ! x coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: ys(ny)               ! y coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: zp(nx,ny,nz)         ! z coordinate of grid points in
                               ! physical space (m)
  REAL :: xpt(maxsnd)          ! location to find in x coordinate (m)
  REAL :: ypt(maxsnd)          ! location to find in y coordinate (m)
  INTEGER :: ipt(maxsnd)       ! i index to the west of desired
                               ! location
  INTEGER :: jpt(maxsnd)       ! j index to the south of desired
                               ! location
  INTEGER :: is_good(maxsnd)   ! Zero when sounding is outside the grid
!
!-----------------------------------------------------------------------
!
!  Arguments -- model data
!
!-----------------------------------------------------------------------
!
  REAL :: uprt   (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: vprt   (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: wprt   (nx,ny,nz)    ! Perturbation w-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: qvprt  (nx,ny,nz)    ! Perturbation water vapor specific
                               ! humidity (kg/kg)

  REAL :: ubar   (nx,ny,nz)    ! Base state u-velocity (m/s)
  REAL :: vbar   (nx,ny,nz)    ! Base state v-velocity (m/s)
  REAL :: wbar   (nx,ny,nz)    ! Base state w-velocity (m/s)
  REAL :: ptbar  (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz)    ! Base state pressure (Pascal)
  REAL :: qvbar  (nx,ny,nz)    ! Base state water vapor specific
                               ! humidity (kg/kg)

  REAL :: tke    (nx,ny,nz)
  REAL :: qc     (nx,ny,nz)
  REAL :: qi     (nx,ny,nz)

  REAL :: dxfld(nx)
  REAL :: dyfld(ny)
  REAL :: rdxfld(nx)
  REAL :: rdyfld(ny)
  REAL :: slopey(nx,ny,nz)
  REAL :: alphay(nx,ny,nz)
  REAL :: betay(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Arguments -- Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  REAL :: su(nz,maxsnd)
  REAL :: sv(nz,maxsnd)
  REAL :: sw(nz,maxsnd)
  REAL :: stheta(nz,maxsnd)
  REAL :: sqv(nz,maxsnd)
  REAL :: sqc(nz,maxsnd)
  REAL :: sqi(nz,maxsnd)
  REAL :: stke(nz,maxsnd)
  REAL :: spres(nz,maxsnd)
  REAL :: shght(nz,maxsnd)
  REAL :: selev(maxsnd)
  INTEGER :: nlevs
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------
!
  REAL :: aint2d
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk,isnd
  REAL :: w2,w3
  REAL :: t1,t2,t3,hmid,tmid,qvsat,rh
  INTEGER :: iorder
  PARAMETER (iorder = 2)
  REAL :: pntint2d
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nlevs=nz-2

  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,zp,                                         &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      shght(k,isnd)=pntint2d(nx,ny,                                     &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),zp(1,1,k),                           &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
    shght(1,isnd)=shght(2,isnd)-(shght(3,isnd)-shght(2,isnd))
  END DO
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=ptbar(i,j,k)+ptprt(i,j,k)
      END DO
    END DO
  END DO
  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tem1,                                       &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      stheta(k,isnd)=pntint2d(nx,ny,                                    &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tem1(1,1,k),                         &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=pbar(i,j,k)+pprt(i,j,k)
      END DO
    END DO
  END DO
  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tem1,                                       &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      spres(k,isnd)=pntint2d(nx,ny,                                     &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tem1(1,1,k),                         &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
! Process qv
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=qvbar(i,j,k)+qvprt(i,j,k)
      END DO
    END DO
  END DO
  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tem1,                                       &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      sqv(k,isnd)=pntint2d(nx,ny,                                       &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tem1(1,1,k),                         &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
! Process TKE
!
  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tke,                                        &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      stke(k,isnd)=pntint2d(nx,ny,                                      &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tke(1,1,k),                          &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
! Process qc
!
  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,qc,                                         &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      sqc(k,isnd)=pntint2d(nx,ny,                                       &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),qc(1,1,k),                           &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
! Process qi
!
  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,qi,                                         &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      sqi(k,isnd)=pntint2d(nx,ny,                                       &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),qi(1,1,k),                           &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get height at scalar points, since zp was defined at w points.
!
!-----------------------------------------------------------------------
!
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    selev(isnd)=shght(2,isnd)
    DO k=1,nz-1
      shght(k,isnd)=0.5*(shght(k+1,isnd)+shght(k,isnd))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total u wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(ubar(i,j,k)+ubar(i+1,j,k))+                    &
                    0.5*(uprt(i,j,k)+uprt(i+1,j,k))
      END DO
    END DO
  END DO

  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tem1,                                       &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      su(k,isnd)=pntint2d(nx,ny,                                        &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tem1(1,1,k),                         &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total v wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(vbar(i,j,k)+vbar(i,j+1,k)) +                   &
                    0.5*(vprt(i,j,k)+vprt(i,j+1,k))
      END DO
    END DO
  END DO

  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tem1,                                       &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      sv(k,isnd)=pntint2d(nx,ny,                                        &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tem1(1,1,k),                         &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total w wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(wbar(i,j,k)+wbar(i,j,k+1)) +                   &
                    0.5*(wprt(i,j,k)+wprt(i,j,k+1))
      END DO
    END DO
  END DO

  CALL setdrvy(nx,ny,nz,                                                &
               1,nx-1,1,ny-1,1,nz-1,                                    &
               dyfld,rdyfld,tem1,                                       &
               slopey,alphay,betay)
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      sw(k,isnd)=pntint2d(nx,ny,                                        &
               1,nx-1,1,ny-1,                                           &
               iorder,xs,ys,xpt(isnd),ypt(isnd),                        &
               ipt(isnd),jpt(isnd),tem1(1,1,k),                         &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,k),alphay(1,1,k),betay(1,1,k))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get a value at the surface, by extrapolating from the
!  2nd and third levels.
!
!-----------------------------------------------------------------------
!
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    shght(1,isnd)=selev(isnd)
    w3=(shght(1,isnd)-shght(2,isnd))                                    &
        /(shght(3,isnd)-shght(2,isnd))
    w2=(1.-w3)
    su(1,isnd)=w2*    su(2,isnd) + w3*    su(3,isnd)
    sv(1,isnd)=w2*    sv(2,isnd) + w3*    sv(3,isnd)
    sw(1,isnd)=0.5*sw(2,isnd)
    IF(stheta(3,isnd) > stheta(2,isnd)) THEN
      stheta(1,isnd)=w2*stheta(2,isnd) + w3*stheta(3,isnd)
    ELSE
      stheta(1,isnd)=stheta(2,isnd)
    END IF
!
!-----------------------------------------------------------------------
!
!  Integrate downward to get the pressure at level 1.
!
!-----------------------------------------------------------------------
!
    t3=stheta(3,isnd)*(spres(3,isnd)/100000.)**rddcp
    t2=stheta(2,isnd)*(spres(2,isnd)/100000.)**rddcp
    hmid=0.5*(shght(2,isnd)+shght(1,isnd))
    tmid=t3+((shght(3,isnd)-hmid)/                                      &
            (shght(3,isnd)-shght(2,isnd)))*(t2-t3)
    spres(1,isnd)=spres(2,isnd)*                                        &
           EXP(g*(shght(2,isnd)-shght(1,isnd))/(rd*tmid))
!
!-----------------------------------------------------------------------
!
!  Use constant RH to extrapolate qv to level 1.
!
!-----------------------------------------------------------------------
!
    qvsat = f_qvsat( spres(2,isnd), t2 )     ! saturated S.H.
    rh=AMIN1((sqv(2,isnd)/qvsat),1.)         ! R.H.
    t1=stheta(1,isnd)*(spres(1,isnd)/100000.)**rddcp
    qvsat = f_qvsat( spres(1,isnd), t1 )
    sqv(1,isnd)=rh*qvsat
  END DO
  RETURN
END SUBROUTINE colint
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CNVSND                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cnvsnd(su,sv,sw,stheta,spres,sqv,slon,                       &
           sdrct,ssped,somega,stemp,sdewp,nlevs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Converts units of data extracted from ARPS history data to
!  those required of sounding data. Determines direction and
!  speed from u and v wind components.
!
!  Dew-point formula from Bolton, 1980, MWR pp 1046-1053.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  October, 1992 (K. Brewster)
!  Conversion to ARPS 3.0.
!
!  10/28/1992 (K. Brewster)
!  Special allowance for low qv-to-dew pt
!
!  04/10/2005 (K. Brewster)
!  Addition of vertical velocity processing.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    su       Sounding u wind component.  (m/s)
!    sv       Sounding v wind component.  (m/s)
!    sw       Sounding v wind component.  (m/s)
!    stheta   Sounding potential temperature (K).
!    spres    Sounding pressure. (Pascals)
!    sqv      Sounding specific humidity
!    slon     Sounding longitude
!    nlevs    Number of above-ground sounding levels.
!
!  OUTPUT:
!
!    spres    Sounding pressure. (mb)
!    sdrct    Sounding wind direction (degrees from north)
!    ssped    Sounding wind speed (m/s)
!    somega   Sounding omega vertical velocity (Pa/s)
!    stemp    Sounding temperature (degrees C)
!    sdewp    Sounding dew point temperature (degrees C)
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
!  Input arguments
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nlevs             ! Number of above-ground sounding levels
  REAL :: su    (nlevs)        ! Sounding u wind component (m/s)
  REAL :: sv    (nlevs)        ! Sounding v wind component (m/s)
  REAL :: sw    (nlevs)        ! Sounding w wind component (m/s)
  REAL :: stheta(nlevs)        ! Sounding potential temperature (K)
  REAL :: spres (nlevs)        ! Sounding pressure. (Pascals)
  REAL :: sqv   (nlevs)        ! Sounding specific humidity (g/g)
  REAL :: slon                 ! Sounding longitude (degrees E)
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  REAL :: sdrct(nlevs)         ! Sounding wind direction
                               ! (degrees from north)
  REAL :: ssped(nlevs)         ! Sounding wind speed (m/s)
  REAL :: somega(nlevs)        ! Sounding verticel velocity (Pa/s)
  REAL :: stemp(nlevs)         ! Sounding temperature (degrees C)
  REAL :: sdewp(nlevs)         ! Sounding dew point temperature
                               ! (degrees C)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k
  REAL :: smix,e,bige,alge,rho
!
  DO k=1,nlevs
!
!-----------------------------------------------------------------------
!
!  Convert u,v to direction and speed
!
!-----------------------------------------------------------------------
!
    CALL uvrotdd(1,1,slon,su(k),sv(k),sdrct(k),ssped(k))
!
!-----------------------------------------------------------------------
!
!  Convert pressure from Pascals to mb
!
!-----------------------------------------------------------------------
!
    spres(k)=spres(k)*0.01
!
!-----------------------------------------------------------------------
!
!  Convert theta to temperature in degrees C.
!
!-----------------------------------------------------------------------
!
    stemp(k)=stheta(k)*((spres(k)/1000.)**rddcp)
    stemp(k)=stemp(k)-273.15
!
!-----------------------------------------------------------------------
!
!  Convert w to omega.  For simplicity use hydrostatic approximation.
!  Use w=dz/dt omega=(dz/dt)*(dp/dz) and hydrostatic relation for dp/dz
!  dp/dz=-rho*g and rho=p/(rd*Tv)   Tv=T*(1.0+0.61*qv)
!
!-----------------------------------------------------------------------
!
    rho=(100.0*spres(k))/(rd*(stemp(k)+273.15)*(1.0+0.61*sqv(k)))
    somega(k)=-sw(k)*rho*g
!
!-----------------------------------------------------------------------
!
!  Convert specific humidity to dew-point in degrees C.
!
!-----------------------------------------------------------------------
!
    IF( sqv(k) > 0.) THEN
      smix=sqv(k)/(1.-sqv(k))
      e=(spres(k)*smix)/(0.62197 + smix)
      bige=e/( 1.001 + ( (spres(k) - 100.) / 900.) * 0.0034)
      alge = ALOG(bige/6.112)
      sdewp(k)= (alge * 243.5) / (17.67 - alge)
    ELSE
      sdewp(k)= stemp(k)-30.
    END IF
  END DO
!
  RETURN

END SUBROUTINE cnvsnd


!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE COEXTRACT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE coextract(nx,ny,nz,nsnd,zp,ipt,jpt,                   &
           uprt, vprt, wprt, ptprt, pprt, qvprt,                        &
           ubar, vbar, wbar, ptbar, pbar, qvbar, is_good,               &
           su,sv,sw,stheta,spres,shght,sqv,                             &
           selev,nlevs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extract a column of data from ARPS history data in the horizontal
!  located at arrray index ipt, jpt.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  July 2006.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz Dimensions of ARPS grids.
!
!    ipt      i index of grid point for desired sounding
!    jpt      j index of grid point for desired sounding
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    vprt     z component of perturbation velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    is_good  Make no computations if this sounding is not satified
!             the condition.
!
!  OUTPUT:
!
!    su       Extracted u wind component.  (m/s)
!    sv       Extracted v wind component.  (m/s)
!    sw       Extracted vertical velocity. (m/s)
!    stheta   Extracted potential temperature (K).
!    spres    Extracted pressure. (Pascals)
!    shght    Extracted height (meters)
!    sqv      Extracted water vapor mixing ratio (kg/kg).
!    selev    extracted surface elevation (m)
!    nlevs    Number of above-ground sounding levels.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
!  Arguments -- location data
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Dimensions of ARPS grids.
  INTEGER, INTENT(IN) :: nsnd
  REAL,    INTENT(IN) :: zp(nx,ny,nz)      ! z coordinate of grid points in
                                           ! physical space (m)
  INTEGER, INTENT(IN) :: ipt(nsnd)       ! i index to the west of desired
                                         ! location
  INTEGER, INTENT(IN) :: jpt(nsnd)       ! j index to the south of desired
                                         ! location
  INTEGER, INTENT(IN) :: is_good(nsnd)   ! Zero when sounding is outside the grid
!
!-----------------------------------------------------------------------
!
!  Arguments -- model data
!
!-----------------------------------------------------------------------
!
  REAL :: uprt   (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: vprt   (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: wprt   (nx,ny,nz)    ! Perturbation w-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: qvprt  (nx,ny,nz)    ! Perturbation water vapor specific
                               ! humidity (kg/kg)

  REAL :: ubar   (nx,ny,nz)    ! Base state u-velocity (m/s)
  REAL :: vbar   (nx,ny,nz)    ! Base state v-velocity (m/s)
  REAL :: wbar   (nx,ny,nz)    ! Base state w-velocity (m/s)
  REAL :: ptbar  (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz)    ! Base state pressure (Pascal)
  REAL :: qvbar  (nx,ny,nz)    ! Base state water vapor specific
                               ! humidity (kg/kg)

!-----------------------------------------------------------------------
!
!  Arguments -- Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: su(nz,nsnd)
  REAL,    INTENT(OUT) :: sv(nz,nsnd)
  REAL,    INTENT(OUT) :: sw(nz,nsnd)
  REAL,    INTENT(OUT) :: stheta(nz,nsnd)
  REAL,    INTENT(OUT) :: sqv(nz,nsnd)
  REAL,    INTENT(OUT) :: spres(nz,nsnd)
  REAL,    INTENT(OUT) :: shght(nz,nsnd)
  REAL,    INTENT(OUT) :: selev(nsnd)
  INTEGER, INTENT(OUT) :: nlevs
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,isnd
  REAL    :: w2,w3
  REAL    :: t1,t2,t3,hmid,tmid,qvsat,rh
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nlevs = nz-2

  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      shght (k,isnd)= zp(ipt(isnd),jpt(isnd),k)
      stheta(k,isnd)= ptbar(ipt(isnd),jpt(isnd),k)+ptprt(ipt(isnd),jpt(isnd),k)
      spres (k,isnd)=  pbar(ipt(isnd),jpt(isnd),k)+ pprt(ipt(isnd),jpt(isnd),k)
      sqv   (k,isnd)= qvbar(ipt(isnd),jpt(isnd),k)+qvprt(ipt(isnd),jpt(isnd),k)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get height at scalar points, since zp was defined at w points.
!
!-----------------------------------------------------------------------
!
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    selev(isnd)=shght(2,isnd)
    DO k=2,nz-2
      shght(k,isnd)=0.5*(shght(k+1,isnd)+shght(k,isnd))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total u wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    DO k=2,nz-1
      su(k,isnd) = 0.5*(ubar(ipt(isnd),jpt(isnd),k)+ubar(ipt(isnd)+1,jpt(isnd),k))+ &
                   0.5*(uprt(ipt(isnd),jpt(isnd),k)+uprt(ipt(isnd)+1,jpt(isnd),k))
      sv(k,isnd) = 0.5*(vbar(ipt(isnd),jpt(isnd),k)+vbar(ipt(isnd),jpt(isnd)+1,k))+ &
                   0.5*(vprt(ipt(isnd),jpt(isnd),k)+vprt(ipt(isnd),jpt(isnd)+1,k))
      sw(k,isnd) = 0.5*(wbar(ipt(isnd),jpt(isnd),k)+wbar(ipt(isnd),jpt(isnd),k+1))+ &
                   0.5*(wprt(ipt(isnd),jpt(isnd),k)+wprt(ipt(isnd),jpt(isnd),k+1))
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get a value at the surface, by extrapolating from the
!  2nd and third levels.
!
!-----------------------------------------------------------------------
!
  DO isnd=1,nsnd
    if ( is_good(isnd) .eq. 0 ) cycle
    shght(1,isnd)=selev(isnd)
    w3=(shght(1,isnd)-shght(2,isnd))/(shght(3,isnd)-shght(2,isnd))
    w2=(1.-w3)
    su(1,isnd)=w2*    su(2,isnd) + w3*    su(3,isnd)
    sv(1,isnd)=w2*    sv(2,isnd) + w3*    sv(3,isnd)
    sw(1,isnd)=0.5*sw(2,isnd)
    IF(stheta(3,isnd) > stheta(2,isnd)) THEN
      stheta(1,isnd)=w2*stheta(2,isnd) + w3*stheta(3,isnd)
    ELSE
      stheta(1,isnd)=stheta(2,isnd)
    END IF
!
!-----------------------------------------------------------------------
!
!  Integrate downward to get the pressure at level 1.
!
!-----------------------------------------------------------------------
!
    t3=stheta(3,isnd)*(spres(3,isnd)/100000.)**rddcp
    t2=stheta(2,isnd)*(spres(2,isnd)/100000.)**rddcp
    hmid=0.5*(shght(2,isnd)+shght(1,isnd))
    tmid=t3+((shght(3,isnd)-hmid)/                                      &
            (shght(3,isnd)-shght(2,isnd)))*(t2-t3)
    spres(1,isnd)=spres(2,isnd)*                                        &
           EXP(g*(shght(2,isnd)-shght(1,isnd))/(rd*tmid))
!
!-----------------------------------------------------------------------
!
!  Use constant RH to extrapolate qv to level 1.
!
!-----------------------------------------------------------------------
!
    qvsat = f_qvsat( spres(2,isnd), t2 )     ! saturated S.H.
    rh=AMIN1((sqv(2,isnd)/qvsat),1.)         ! R.H.
    t1=stheta(1,isnd)*(spres(1,isnd)/100000.)**rddcp
    qvsat = f_qvsat( spres(1,isnd), t1 )
    sqv(1,isnd)=rh*qvsat

  END DO

  RETURN
END SUBROUTINE coextract
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE EXTENVPRF3                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extenvprf3(nx,ny,nz,lvlprof,                                 &
                     xs,ys,zps,usc,vsc,rfrct,                           &
                     radarx,radary,radarz,radius,rngmax,                &
                     zsnd,ktsnd,usnd,vsnd,rfrctsnd,ienvstat)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Finds environmental profile around radar location from gridded data.
!  This version of extenvprf assumes the 3d fields are all precalculated
!  to save steps when processing multiple soundings in one program.
!
!  AUTHOR:
!  Keith Brewster, CAPS
!  10-Sept-2008
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: lvlprof
  REAL, INTENT(IN)    :: xs(nx)
  REAL, INTENT(IN)    :: ys(ny)
  REAL, INTENT(IN)    :: zps(nx,ny,nz)
  REAL, INTENT(IN)    :: usc(nx,ny,nz)
  REAL, INTENT(IN)    :: vsc(nx,ny,nz)
  REAL, INTENT(IN)    :: rfrct(nx,ny,nz)
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radarz
  REAL, INTENT(IN)    :: radius
  REAL, INTENT(IN)    :: rngmax
  REAL, INTENT(IN)    :: zsnd(lvlprof)
  REAL, INTENT(OUT)   :: ktsnd(lvlprof)
  REAL, INTENT(OUT)   :: usnd(lvlprof)
  REAL, INTENT(OUT)   :: vsnd(lvlprof)
  REAL, INTENT(OUT)   :: rfrctsnd(lvlprof)
  INTEGER, INTENT(OUT):: ienvstat

  INTEGER, PARAMETER :: nptsmin = 3

  INTEGER :: i,j,k,kbot,ktop,kext
  INTEGER :: iradar,jradar
  INTEGER :: ibeg,iend,jbeg,jend
  INTEGER :: knt,knttot
  REAL :: dx,dy
  REAL :: rngavg,radius2,dist2,distmin,whigh,wlow,accept,rfrgrad

  INCLUDE 'mp.inc'

  INTEGER :: ioutdomain

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ienvstat=0

  dx=xs(2)-xs(1)
  dy=ys(2)-ys(1)

!
!-----------------------------------------------------------------------
!
!  Find the index of scalar grid location that is closest to the radar.
!  Set bounds of searching to be within a box of
!  radarx-radius < xs < radarx+radius
!  radary-radius < ys < radary+radius
!
!-----------------------------------------------------------------------
!
  rngavg=radius
  radius2=rngavg*rngavg
  iradar=1+nint((radarx-xs(1))/dx)
  jradar=1+nint((radary-ys(1))/dy)
  IF(myproc == 0) print *, ' iradar=',iradar,', jradar=',jradar

  ibeg=iradar-(int(radius/dx)+1)
  ibeg=min(max(ibeg,2),nx-2)

  iend=iradar+(int(radius/dx)+1)
  iend=min(max(iend,2),nx-2)

  jbeg=jradar-(int(radius/dy)+1)
  jbeg=min(max(jbeg,2),ny-2)

  jend=jradar+(int(radius/dy)+1)
  jend=min(max(jend,2),ny-2)

  IF(myproc == 0) print *, ' ibeg = ',ibeg,', iend= ',iend
  IF(myproc == 0) print *, ' jbeg = ',jbeg,', jend= ',jend
!
! Special check if radar location is outside subdomain
!
  ioutdomain = 0
  IF( iradar > (nx-2) .OR. iradar < 2 .OR.                              &
      jradar > (ny-2) .OR. jradar < 2 ) THEN
    ioutdomain = 1
  END IF
  CALL mptotali(ioutdomain)

  IF (ioutdomain == nprocs) THEN   ! It is outside of all subdomains
    knt=0
    distmin=1.0E32
    DO j=jbeg,jend
      DO i=ibeg,iend
        dist2 = (xs(i)-radarx)*(xs(i)-radarx)                           &
               +(ys(j)-radary)*(ys(j)-radary)
        IF(dist2 < radius2) knt=knt+1
        distmin=min(distmin,dist2)
      END DO
    END DO
    CALL mpminr(distmin)     ! minimum in global domain
    distmin=sqrt(distmin)
    CALL mptotali(knt)      ! all should count.
    IF( knt < nptsmin ) THEN
      IF( distmin < rngmax ) THEN
        IF (myproc == 0) THEN
        WRITE(6,'(a,f9.0,a)')                                           &
        ' Too few points within radar radius ',(0.001*radius),' km.'
        WRITE(6,'(a,f9.0,a)')                                           &
        ' Expanding averaging radius to rngmax ',(0.001*rngmax),' km.'
        END IF
        rngavg=rngmax
        radius2=rngavg*rngavg

        ibeg=iradar-(int(rngmax/dx)+1)
        ibeg=min(max(ibeg,2),nx-2)

        iend=iradar+(int(rngmax/dx)+1)
        iend=min(max(iend,2),nx-2)

        jbeg=jradar-(int(rngmax/dy)+1)
        jbeg=min(max(jbeg,2),ny-2)

        jend=jradar+(int(rngmax/dy)+1)
        jend=min(max(jend,2),ny-2)

      ELSE

        IF (myproc == 0) THEN
        WRITE(6,'(a,f12.1,a)')                                         &
        ' extenvprf: Minimum distance ',(0.001*distmin),' km.'
        WRITE(6,'(a,f9.0,a)')                                          &
        ' No points in the domain within radius ',(0.001*rngmax),' km.'
        WRITE(6,'(a)') ' Radar range outside of domain.  Error exit'
        END IF

        CALL flush(6)
        ienvstat=-1
        RETURN

      END IF

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize sums to zero
!
!-----------------------------------------------------------------------
!
  DO k=1,lvlprof
    ktsnd(k)=0.
    usnd(k)=0.
    vsnd(k)=0.
    rfrctsnd(k)=0.
  END DO
  knt=0
  knttot=0
  DO j=jbeg,jend
    DO i=ibeg,iend
!
!-----------------------------------------------------------------------
!
!  Is this point within the ARPS domain?
!  Since the ARPS grid is Cartesian, need only compare
!  the external grid coordinates to the x and y limits
!  of the ARPS grid.
!
!-----------------------------------------------------------------------
!
      knttot=knttot+1
      dist2 = (xs(i)-radarx)*(xs(i)-radarx)                             &
             +(ys(j)-radary)*(ys(j)-radary)
      IF(dist2 < radius2 ) THEN
        knt=knt+1
!
!-----------------------------------------------------------------------
!
!  Interpolate external data in vertical onto profile
!  arrays.
!
!-----------------------------------------------------------------------
!
        DO k=1,lvlprof
          IF(zps(i,j,2) <= zsnd(k)) THEN
            DO kext=3,nz-1
              IF(zps(i,j,kext) >= zsnd(k)) EXIT
            END DO
            IF(kext > nz-1) EXIT
            whigh=(zsnd(k)-zps(i,j,kext-1))/                            &
                  (zps(i,j,kext)-zps(i,j,kext-1))
            wlow=1.-whigh
            ktsnd(k)=ktsnd(k)+1.
            usnd(k)=usnd(k)+                                            &
                  whigh*usc(i,j,kext)+wlow*usc(i,j,kext-1)
            vsnd(k)=vsnd(k)+                                            &
                  whigh*vsc(i,j,kext)+wlow*vsc(i,j,kext-1)
            rfrctsnd(k)=rfrctsnd(k)+                                    &
                  whigh*rfrct(i,j,kext)+wlow*rfrct(i,j,kext-1)
          END IF
        END DO
      END IF
    END DO
  END DO

  CALL mptotali(knt)
  CALL mptotali(knttot)
  CALL mpsumr(ktsnd,lvlprof)
  CALL mpsumr(usnd,lvlprof)
  CALL mpsumr(vsnd,lvlprof)
  CALL mpsumr(rfrctsnd,lvlprof)

  IF (myproc == 0) WRITE(6,'(/a,i6,a,/a,i6,a/)')                        &
       '  extenvprf: found ',knt,' points within radius',               &
       '  of ',knttot,' checked.'

  IF(knt < nptsmin) THEN
    IF (myproc == 0) THEN
    WRITE(6,'(a,f9.0,a)')                                               &
      ' Too few domain points within radar radius ',(0.001*rngavg),' km.'
    WRITE(6,'(a)') ' Radar range outside of domain.  Error exit'
    END IF
    CALL flush(6)
    ienvstat=-2
    RETURN
  END IF

  accept=0.3*float(knt)
!
!-----------------------------------------------------------------------
!
!  Find lowest height with data
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a)') ' Finding range of mean profile data ...'
  DO k=1,lvlprof-1
    IF(ktsnd(k) > accept) EXIT
  END DO
  kbot=k

!
! Catastrophic failure.  Abort now.
!

  IF (kbot > lvlprof) THEN
    IF(myproc == 0) WRITE(6,'(a)') ' No acceptable data in EXTENVPRF.  Abort.'
    ienvstat=-1
    RETURN
  END IF

!
!-----------------------------------------------------------------------
!
!  Find highest height with data
!
!-----------------------------------------------------------------------
!
  DO k=lvlprof,2,-1
    IF(myproc == 0) WRITE(6,'(a,f10.2,a,f6.0,a,f10.0)') ' z = ',zsnd(k),&
                  ' knt = ',ktsnd(k),' accept = ',accept
    IF(ktsnd(k) > accept) EXIT
  END DO
  ktop=k
!
  IF(myproc == 0) WRITE(6,'(a,f10.2,a,f10.2,a)')                       &
               ' Height of data for wind profile spans from ',         &
                 zsnd(kbot),' to ',zsnd(ktop),' meters.'

  IF((ktop - kbot) < 2 ) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a,i6,a,i6)')                                            &
         ' Problem creating evironmental sounding, kbot= ',kbot,        &
         ' ktop=',ktop
      WRITE(6,'(a)') ' Radar range may be outside of domain.  Error exit'
    END IF
    ienvstat=-2
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Divide through to find average
!
!-----------------------------------------------------------------------
!
  DO k=kbot,ktop
    usnd(k)=usnd(k)/ktsnd(k)
    vsnd(k)=vsnd(k)/ktsnd(k)
    rfrctsnd(k)=rfrctsnd(k)/ktsnd(k)
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "below-ground"
!  Zero gradiant assumed in usnd, vsnd.
!  Constant gradiant assumed in rfrctsnd.
!
!-----------------------------------------------------------------------
!
  rfrgrad=(rfrctsnd(kbot+1)-rfrctsnd(kbot))/(zsnd(kbot+1)-zsnd(kbot))
  DO k=kbot-1,1,-1
    usnd(k)=usnd(kbot)
    vsnd(k)=vsnd(kbot)
    rfrctsnd(k)=rfrctsnd(kbot)+(zsnd(k)-zsnd(kbot))*rfrgrad
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "above-top"
!  Zero gradiant assumed.
!
!-----------------------------------------------------------------------
!
  rfrgrad=(rfrctsnd(ktop)-rfrctsnd(ktop-1))/(zsnd(ktop)-zsnd(ktop-1))
  DO k=ktop+1,lvlprof
    usnd(k)=usnd(ktop)
    vsnd(k)=vsnd(ktop)
    rfrctsnd(k)=rfrctsnd(ktop)+(zsnd(k)-zsnd(ktop))*rfrgrad
  END DO

  RETURN
END SUBROUTINE extenvprf3
