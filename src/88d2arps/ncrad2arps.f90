!########################################################################
!########################################################################
!#########                                                      #########
!#########                  PROGRAM ncrad2arps                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

PROGRAM ncrad2arps

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar data from NetCDF data files and remaps data on ARPS
! grid to be used for analysis or display.
! CASA radar data, for example, are stored in netCDF files.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (May 2006)
!
! MODIFICATIONS:
!
! 1 Feb 2007 Keith Brewster
! Added processing for gzipped files.
!
! 10 April 2008 Keith Brewster
! Added checks for ground clutter filter state.
!
! 24 May 2008 Keith Brewster
! Added -fntime option to assign output filename based on input filename
! rather than the time stored
!
! 07 July 2009 Keith Brewster
! Added processing for OU-PRIME radar data.
! Added -rtopts to eliminate surprise and premature processing of
! input namelist variables.
!
! 17 July 2009 Keith Brewster
! Added processing for KOUN dual-pol NEXRAD radar data.
!
! 03 Sept 2009 Youngsun Jung
! Added time interpolation and gridtilt option for KOUN polarimetric
! variables -> output in EnKF input format.
! Added processing for Zdr, Kdp, RhoHV.
! (Potentially, includes Zvv and Zdp)
!
! 16 March 2010 Keith Brewster
! Installed input options processing via input file.
! Updated MPI strategy for better load-balancing.
!
! 24 August 2010 Keith Brewster
! Added processing for TDWR WDSS-II files.
!
! 15 December 2010 Keith Brewster
! Added reading of OU-PRIME Tier-2a style data
!
! 22 August 2011 Keith Brewster
! Added vadsrc to the call to wtvadprf
!
! 30 September 2011 Keith Brewster
! Added additional flexibility for variable naming with input variables
! refcdfname,velcdfname, etc.  Restored original functionality of
! multiple list files via tintrpopt input variable.
!
! 07 December 2011 Keith Brewster
! Made interpolation for on-tilt (gridtilt output) consistent with
! the on-grid-level (column output) interpolation.
!
! 14 December 2011 Youngsun Jung
! Add EnKF format support to Foray data
!
! 19 April 2012 Keith Brewster
! Renamed *cdfname variables to *varname to accomodate use of names
! in other data file types
!
! Y. Wang (05/08/2012)
! Replace the non-standard calls of getarg with the Fortran 2003
! standard call of GET_COMMAND_ARGUMENT.
!
! Y. Jung (05/10/2012)
! Temporary fix for Foray data from UMass radar in "RadialIntp_rdr".
!
! 07 Sept 2012 Keith Brewster
! Added remapping for dual-pol variables.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Include files
!
!------------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'remapcst.inc'
  INCLUDE 'netcdf.inc'
  INCLUDE 'mp.inc'

!------------------------------------------------------------------------
!
! Dimensions and parameters
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: maxncradfile=90
  INTEGER, PARAMETER :: maxvar = 10
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: itime1970 = 315619200

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER :: nx, ny, nz, nzsoil, nstyps
  REAL    :: radlat, radlon, radalt
  REAL    :: radarx, radary

!------------------------------------------------------------------------
!
! NetCDF Radar Data Variables
!
!------------------------------------------------------------------------
!
   REAL, ALLOCATABLE :: refl(:,:)
   REAL, ALLOCATABLE :: attref(:,:)
   REAL, ALLOCATABLE :: rhohv(:,:)
   REAL, ALLOCATABLE :: zdr(:,:)
   REAL, ALLOCATABLE :: phidp(:,:)
   REAL, ALLOCATABLE :: snr(:,:)
   REAL, ALLOCATABLE :: radv(:,:)
   INTEGER, ALLOCATABLE :: gateqc(:,:)
!
!------------------------------------------------------------------------
!
! Radar Data Variables for Gridtilt Option
!
!------------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: rdr_p(:,:,:,:)    ! in polar coordinates
  REAL, ALLOCATABLE :: rdr_tmp(:,:,:)
  INTEGER, PARAMETER :: numrad = 360     ! 360 deg.
  INTEGER, PARAMETER :: numrad2 = 362     ! 362 deg.( Should overlap)
  INTEGER, PARAMETER :: numvar = 7       ! number of dualpol var
               ! 1:Vr, 2:Zhh, 3:Zvv, 4:Zdr, 5:Zdp, 6:Kdp, 7:RhoHV
  INTEGER :: nlist,ilist,jlist,tlevel,idp
  INTEGER :: gridtilt
  INTEGER :: rdrknt_tot(numvar)
  INTEGER :: rdrknt(numvar,maxlistf)
  INTEGER :: reftime
  CHARACTER (LEN=1) :: ref_garb
  INTEGER :: ref_yr,ref_mon,ref_day,ref_hr,ref_min,ref_sec
  INTEGER :: ref_t,cur_time
  REAL, PARAMETER :: missing = -999.9
  REAL, PARAMETER :: rhohv_thr = 0.6
  CHARACTER (LEN=3) :: dualvarname(numvar)
  INTEGER :: scantime(maxncradfile,2)

  REAL, ALLOCATABLE :: elvang(:)
  REAL, ALLOCATABLE :: azmvol(:,:)    ! azimuth angles for each radial and elevation
  REAL, ALLOCATABLE :: elvvol(:,:)    ! elevations for for each radial and elevation
  REAL, ALLOCATABLE :: gridtilt_height(:,:,:)
  REAL, ALLOCATABLE :: gridtilt_range(:,:,:)
  REAL, ALLOCATABLE :: gridtilt_slr(:,:)
  REAL, ALLOCATABLE :: gridtilt_azm(:,:)
  REAL, ALLOCATABLE :: gridtilt_rdr(:,:,:,:)
  INTEGER, ALLOCATABLE :: gridtilt_time(:)
  REAL, ALLOCATABLE :: gridtilt_rdr_mean(:)

  CHARACTER (LEN=256) :: outfn
  REAL :: dpcheck
!
!------------------------------------------------------------------------
!
! Variables for subroutine remap3dcol
!
!------------------------------------------------------------------------
!
  INTEGER, PARAMETER :: nt=1
  INTEGER, PARAMETER :: exbcbufsz=1
  INTEGER, PARAMETER :: bkgopt=1
  INTEGER, PARAMETER :: shropt=1
  INTEGER, PARAMETER :: rfropt = 1  ! 4/3rds earth refraction option
  INTEGER, PARAMETER :: nsort = 601
  REAL, PARAMETER :: sortmin = -150.0
  REAL, PARAMETER :: sortmax = 150.0
  INTEGER, PARAMETER :: nzsnd = 201
  REAL,    PARAMETER :: dzsnd=100.
  INTEGER :: ktsnd(nzsnd)
  REAL :: zsnd(nzsnd)             ! hgt levels for refractivity sounding
  REAL :: usnd(nzsnd)
  REAL :: vsnd(nzsnd)
  REAL :: rfrsnd(nzsnd)           ! refractivity sounding
  REAL, PARAMETER :: rfrconst = (-1./(4.*6371.E03))   ! -1/(4e)

  REAL, ALLOCATABLE :: arfdata(:,:) ! attenuated refl (dBZ) data
  REAL, ALLOCATABLE :: refdata(:,:) ! refl (dBZ) data
  REAL, ALLOCATABLE :: veldata(:,:) ! velocity data (m/s)
  REAL, ALLOCATABLE :: spwdata(:,:)
  REAL, ALLOCATABLE :: bkgvel(:,:)
  REAL, ALLOCATABLE :: unfvdata(:,:) ! unfolded velocity data (m/s)
  REAL, ALLOCATABLE :: rhvdata(:,:) ! rho-HV data
  REAL, ALLOCATABLE :: zdrdata(:,:) ! Zdr data
  REAL, ALLOCATABLE :: kdpdata(:,:) ! Phi data
  REAL, ALLOCATABLE :: snstvty(:) ! radar sensitivity
  REAL, ALLOCATABLE :: rtem(:,:)  ! temporary array for radar data processing
  INTEGER, ALLOCATABLE :: time(:) ! time offset from itimfrst
  INTEGER, ALLOCATABLE :: gcfst(:) ! Ground clutter filter state
  REAL, ALLOCATABLE :: azim(:)    ! azimuth angle for each radial(degree)
  REAL, ALLOCATABLE :: gtspc(:)   ! gate spacing for each radial (meter)
  REAL, ALLOCATABLE :: elev(:)    ! elevation angle for each radial (degree)
  REAL, ALLOCATABLE :: vnyq(:)    ! Nyquist velocities (meters/sec)
  INTEGER, ALLOCATABLE :: istrgate(:)
  INTEGER, ALLOCATABLE :: bgate(:)
  INTEGER, ALLOCATABLE :: egate(:)
  INTEGER :: itimfrst             ! time of first radial in volume
  INTEGER :: irfirstg             ! range to first gate (meters)
  INTEGER :: igatesp              ! gate spacing (meters)
  INTEGER :: irefgatsp            ! reflectivity gate spacing (meters)
  INTEGER :: ivelgatsp            ! velocity gate spacing (meters)
  REAL    :: elv                  ! elevation angle for this tilt
  REAL    :: rfirstg              ! range to first gate (m)
  REAL    :: gatesp               ! gate spacing (meters)

  INTEGER, ALLOCATABLE :: kntrgat(:,:)
  INTEGER, ALLOCATABLE :: kntrazm(:)
  INTEGER :: kntrelv

  INTEGER, ALLOCATABLE :: kntvgat(:,:)
  INTEGER, ALLOCATABLE :: kntvazm(:)
  INTEGER :: kntvelv

  INTEGER, ALLOCATABLE :: timevolr(:,:)
  INTEGER, ALLOCATABLE :: timevolv(:,:)
  REAL, ALLOCATABLE :: vnyqvol(:,:)

  REAL, ALLOCATABLE :: rngrvol(:,:)
  REAL, ALLOCATABLE :: azmrvol(:,:)
  REAL, ALLOCATABLE :: elvrvol(:,:)
  REAL, ALLOCATABLE :: elvmnrvol(:)
  REAL, ALLOCATABLE :: refvol(:,:,:)
  REAL, ALLOCATABLE :: rhvvol(:,:,:)
  REAL, ALLOCATABLE :: zdrvol(:,:,:)
  REAL, ALLOCATABLE :: kdpvol(:,:,:)
  REAL, ALLOCATABLE :: rngvvol(:,:)
  REAL, ALLOCATABLE :: azmvvol(:,:)
  REAL, ALLOCATABLE :: elvvvol(:,:)
  REAL, ALLOCATABLE :: elvmnvvol(:)
  REAL, ALLOCATABLE :: velvol(:,:,:)
  INTEGER, ALLOCATABLE :: kntbin(:)

  REAL, ALLOCATABLE :: rxvol(:,:,:)
  REAL, ALLOCATABLE :: ryvol(:,:,:)
  REAL, ALLOCATABLE :: rzvol(:,:,:)

  INTEGER, ALLOCATABLE :: icolp(:)
  INTEGER, ALLOCATABLE :: jcolp(:)
  REAL, ALLOCATABLE :: xcolp(:)
  REAL, ALLOCATABLE :: ycolp(:)
  REAL, ALLOCATABLE :: zcolp(:,:)
  LOGICAL, ALLOCATABLE :: havdat(:,:)

  REAL, ALLOCATABLE :: colvel(:,:)
  REAL, ALLOCATABLE :: colref(:,:)
  REAL, ALLOCATABLE :: colnyq(:,:)
  REAL, ALLOCATABLE :: coltim(:,:)
  REAL, ALLOCATABLE :: colrhv(:,:)
  REAL, ALLOCATABLE :: colzdr(:,:)
  REAL, ALLOCATABLE :: colkdp(:,:)
!
!------------------------------------------------------------------------
!
! Variables for VAD profiles.
!
!------------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: vadhgt(:)
  REAL, ALLOCATABLE :: vaddir(:)
  REAL, ALLOCATABLE :: vadspd(:)
!
!------------------------------------------------------------------------
!
! Variables for background and ARPS grid data
!
!------------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x(:)
  REAL, ALLOCATABLE :: y(:)
  REAL, ALLOCATABLE :: z(:)
  REAL, ALLOCATABLE :: zp(:,:,:)
  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
  REAL, ALLOCATABLE :: hterain(:,:)
  REAL, ALLOCATABLE :: mapfct(:,:,:) ! Map factors at scalar, u and v points
  REAL, ALLOCATABLE :: j1    (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as - d( zp )/d( x ).
  REAL, ALLOCATABLE :: j2    (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as - d( zp )/d( y ).
  REAL, ALLOCATABLE :: j3    (:,:,:) ! Coordinate transformation Jacobian defined
                                     ! as d( zp )/d( z ).
  REAL, ALLOCATABLE :: tem1d1(:)
  REAL, ALLOCATABLE :: tem1d2(:)
  REAL, ALLOCATABLE :: tem2d(:,:)
  REAL, ALLOCATABLE :: tem2dyz(:,:)
  REAL, ALLOCATABLE :: tem2dxz(:,:)
  REAL, ALLOCATABLE :: tem2dns(:,:,:)
  REAL, ALLOCATABLE :: tem3d1(:,:,:)
  REAL, ALLOCATABLE :: tem3d2(:,:,:)
  REAL, ALLOCATABLE :: tem3d3(:,:,:)
  REAL, ALLOCATABLE :: tem3d4(:,:,:)
  INTEGER, ALLOCATABLE :: tem2dint(:,:)
  INTEGER, ALLOCATABLE :: soiltyp(:,:,:)
  REAL, ALLOCATABLE :: stypfrct(:,:,:)
  REAL, ALLOCATABLE :: tem3dsoil(:,:,:)
  REAL, ALLOCATABLE :: tem4dsoilns(:,:,:,:)

  REAL, ALLOCATABLE :: qscalar(:,:,:,:)

  REAL, ALLOCATABLE :: zpsoil(:,:,:)
  REAL, ALLOCATABLE :: j3soil(:,:,:)
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)

  REAL, ALLOCATABLE :: u(:,:,:)
  REAL, ALLOCATABLE :: v(:,:,:)
  REAL, ALLOCATABLE :: qv(:,:,:)
  REAL, ALLOCATABLE :: ubar(:,:,:)
  REAL, ALLOCATABLE :: vbar(:,:,:)
  REAL, ALLOCATABLE :: ptbar(:,:,:)
  REAL, ALLOCATABLE :: ptprt(:,:,:)
  REAL, ALLOCATABLE :: pprt(:,:,:)
  REAL, ALLOCATABLE :: pbar(:,:,:)
  REAL, ALLOCATABLE :: qvbar(:,:,:)
  REAL, ALLOCATABLE :: rhostr(:,:,:)

  REAL, ALLOCATABLE :: trigs1(:)
  REAL, ALLOCATABLE :: trigs2(:)
  INTEGER, ALLOCATABLE :: ifax(:)
  REAL, ALLOCATABLE :: wsave1(:)
  REAL, ALLOCATABLE :: wsave2(:)
  REAL, ALLOCATABLE :: vwork1(:,:)
  REAL, ALLOCATABLE :: vwork2(:,:)

  REAL, ALLOCATABLE :: qcumsrc(:,:,:,:)
  REAL, ALLOCATABLE :: prcrate(:,:,:)
  REAL, ALLOCATABLE :: exbcbuf(:)

  REAL, ALLOCATABLE :: grdvel2(:,:,:)
  REAL, ALLOCATABLE :: grdref2(:,:,:)
  REAL, ALLOCATABLE :: wgtvel(:,:,:)
  REAL, ALLOCATABLE :: wgtref(:,:,:)
  REAL, ALLOCATABLE :: wpotvel(:,:,:)
  REAL, ALLOCATABLE :: wpotref(:,:,:)
  REAL, ALLOCATABLE :: gridazm(:,:)
  REAL, ALLOCATABLE :: gridrng(:,:)

  REAL(KIND=8),    ALLOCATABLE :: tem_double(:)
  INTEGER(KIND=2), ALLOCATABLE :: tem_short2d(:,:)
  REAL,            ALLOCATABLE :: tem_real2d(:,:)

!------------------------------------------------------------------------
!
! Data dimensions, determined from the data files
!
!------------------------------------------------------------------------

  INTEGER :: maxgate,maxazim,maxelev

!------------------------------------------------------------------------
!
! Misc. variables
!
!------------------------------------------------------------------------

  CHARACTER(LEN=256) :: infilname
  CHARACTER(LEN=180) :: indir
  CHARACTER(LEN=256) :: fname
  CHARACTER(LEN=256) :: full_fname
  CHARACTER(LEN=256) :: vad_fname
  CHARACTER(LEN=80) :: runname_in
  CHARACTER(LEN=6)  :: varid
  CHARACTER(LEN=6)  :: refid = 'refl3d'
  CHARACTER(LEN=6)  :: velid = 'radv3d'
  CHARACTER(LEN=6)  :: zdrid = 'zdr_3d'
  CHARACTER(LEN=6)  :: kdpid = 'kdp_3d'
  CHARACTER(LEN=6)  :: rhvid = 'rhv_3d'
  CHARACTER(LEN=20) :: varname
  CHARACTER(LEN=20) :: refname = 'Reflectivity'
  CHARACTER(LEN=20) :: velname = 'RawVelocity'
  CHARACTER(LEN=20) :: zdrname = 'Differential Reflect'
  CHARACTER(LEN=20) :: kdpname = 'Differential Phase'
  CHARACTER(LEN=20) :: rhvname = 'Correlation Coeff'
  CHARACTER(LEN=6)  :: refunits = 'dBZ'
  CHARACTER(LEN=6)  :: velunits = 'm/s'
  CHARACTER(LEN=6)  :: zdrunits = 'dBZ'
  CHARACTER(LEN=6)  :: kdpunits = 'deg/km'
  CHARACTER(LEN=6)  :: rhvunits = ' '

  CHARACTER(LEN=NF_MAX_NAME), ALLOCATABLE :: ncvarname(:)
  CHARACTER(LEN=266) :: syscall
  INTEGER :: fvartype(maxncradfile,maxlistf)
  LOGICAL :: gzipped(maxncradfile,maxlistf)
  REAL :: felv(maxncradfile,maxlistf)
  INTEGER :: fitime(maxncradfile,maxlistf)

  INTEGER :: n, i, j, k, kfn, kkfn, istat, istatus
  INTEGER :: nazim,ngate,nvar,nelev,ngate1
  INTEGER :: ivcp,itimcdf,initime,mintimcdf
  INTEGER :: ivar,kelv,icol
  INTEGER :: len_fname
  INTEGER :: velopt,varfill,iradfmt

  LOGICAL :: velproc
  LOGICAL :: dualpproc,dualpdata
  LOGICAL :: vel2d,vel3d
  LOGICAL :: ref2d,ref3d
  LOGICAL :: zdr2d,zdr3d
  LOGICAL :: kdp2d,kdp3d
  LOGICAL :: rhv2d,rhv3d
  LOGICAL :: unfdiag
  LOGICAL :: qc_on
  LOGICAL :: fntime
  LOGICAL :: traditional
  LOGICAL :: rtopts
  LOGICAL :: vad

  LOGICAL :: rho_match

  INTEGER :: iyear,imonth,iday,ihr,imin,isec
  INTEGER :: kyear,kmonth,kday,khr,kmin,ksec
  INTEGER :: kmyear,kmmon,kmday,kmhr,kmmin,kmsec
  INTEGER :: iiyear,iimon,iiday,iihr,iimin,iisec
  INTEGER :: iangle,itilt,dtime
  INTEGER :: kntt2a,kntt2b
  INTEGER :: isource,klev,iradtype
  INTEGER :: nyqset,timeset,vardump
  INTEGER :: kntgate,badgate
  INTEGER :: nxny,ncoltot,ncolp
  INTEGER :: xscale

  REAL :: rdx,dsort,rmin,rmax
  REAL :: rmisval,rngfval,frtime,envtime
  REAL :: refelvmin,refelvmax
  REAL :: refrngmin,refrngmax
  REAL :: badpct

  INTEGER :: remapopt
  REAL :: radius,vconst,vnyquist

  REAL :: tmin,tmax
  REAL :: xrad,yrad

  INTEGER :: maxinfile,infile,lenfn,iloc
  INTEGER :: iarg,jarg,iargc,narg,ifile,jfile,kfile,lfile
  INTEGER :: nprocx_in_cmd, nprocy_in_cmd

  CHARACTER(LEN=4)   :: tsthead
  CHARACTER(LEN=256) :: charg

!------------------------------------------------------------------------
!
! Some tunable parameters - convert to Namelist in future
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: gcfchek = 0
  INTEGER, PARAMETER :: gcopt = 1    ! ground clutter option

  CHARACTER(LEN=8) :: vadsrc
  CHARACTER(LEN=32) :: radarname
  INTEGER :: nncradfn(maxlistf)
  CHARACTER(LEN=256) :: ncradfn(maxncradfile,maxlistf)
  CHARACTER(LEN=256) :: tempstr
  CHARACTER(LEN=256) :: namelistfile

  INTEGER, PARAMETER :: ROOT = 0
  LOGICAL :: back = .TRUE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Initializations
!
!------------------------------------------------------------------------

  refelvmin = 91.0
  refelvmax = -91.0
  dmpfmt=1
  hdf4cmpr=0
  kntt2a=0
  kntt2b=0
  fillref=0
  velproc = .TRUE.
  dualpdata = .FALSE.
  dualpproc = .FALSE.
  unfdiag = .FALSE.
  qc_on = .TRUE.
  ref2d = .FALSE.
  ref3d = .FALSE.
  vel2d = .FALSE.
  vel3d = .FALSE.
  zdr2d = .FALSE.
  zdr3d = .FALSE.
  kdp2d = .FALSE.
  kdp3d = .FALSE.
  rhv2d = .FALSE.
  rhv3d = .FALSE.
  fntime = .FALSE.
  traditional = .FALSE.
  rtopts = .FALSE.
  vadsrc = 'CASA VAD'
  rdrknt_tot=0
  rdrknt=0
  gridtilt = 0
  reftime = 0
  dualvarname=(/'Vr ','Zhh','Zvv','Zdr','Zdp','Kdp','rhv'/)

  medfilt = 0
  irefgatsp = 100   ! default value, replaced by actual data value
  ivelgatsp = 100
  ivcp=0
  myproc = 0
  mp_opt = 0
  nprocx_in_cmd = 1
  nprocy_in_cmd = 1

  dsort = (sortmax-sortmin)/float(nsort-1)

  DO lfile=1,maxlistf
    DO ifile=1,maxncradfile
      ncradfn(ifile,lfile)='dummy'
    END DO
  END DO

  CALL mpinit_proc(0)

!------------------------------------------------------------------------
!
! First check to see if help command line argument is specified.
!
!------------------------------------------------------------------------

  istatus=0
  IF(myproc == ROOT) THEN
    narg = COMMAND_ARGUMENT_COUNT()
    IF(narg > 0 ) THEN
      CALL GET_COMMAND_ARGUMENT(1, charg, istat, istatus )
      IF(charg(1:5) == '-help' .OR. charg(1:6) == '--help') THEN
        WRITE(stdout,'(a,a,a,a,a,/,a)')                                                 &
        ' Usage: ncrad2arps [RADAR_ID] [file_list_file] [-novel] [-dir dirname]',       &
        ' [-hdf n] [-binary] [-ref2d] [-ref3d] [-reffile] [-vel2d] [-vel3d] [-velfile]',&
        ' [-noqc] [-fillref] [-fntime] [-rtopts] [-gridtilt]',                          &
        ' [-reftime yyyymmdd-hhmm] [-input nmfile]', ' < radremap.input',               &
        ' Command line options override the options in the input file.'
        istatus = -1
      END IF
    END IF
  END IF
  CALL mpupdatei(istatus,1)

  IF(istatus < 0) CALL arpsstop ('End of help.',istatus)
!
!------------------------------------------------------------------------
!
! Initialize ncrad_data NAMELIST variables
!
!------------------------------------------------------------------------
!
  radname = 'DMMY'
  runname = 'NULL'
  radband = 3
  nncradfn = 0
  kfile = 0
  namelistfile = ' '
!
!------------------------------------------------------------------------
!
! Go through the arguments to extract "nprocx_in" and "nprocy_in"
! if they are present.  This must be done before call to initremapopt.
!
!------------------------------------------------------------------------
!
  IF(myproc == ROOT) THEN

    IF(narg > 0) THEN
      iarg=1
      DO jarg=1,narg
        CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
        IF(charg(1:10) == '-nprocx_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            READ(charg,'(i3)') nprocx_in_cmd
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (X) = ',nprocx_in_cmd
        ELSE IF(charg(1:10) == '-nprocy_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            READ(charg,'(i3)') nprocy_in_cmd
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (Y) = ',nprocy_in_cmd
        ELSE IF(charg(1:6) == '-input') THEN    ! namelist file
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            namelistfile = charg
          END IF
          WRITE(stdout,'(1x,2a)') 'Namelist file from command line - ',TRIM(namelistfile)
        END IF
        iarg=iarg+1
        IF(iarg > narg) EXIT
      END DO
    END IF ! narg > 0
  END IF ! myproc == 0
!
!------------------------------------------------------------------------
!
! Read the input file and set logical variables based on input integers.
!
!------------------------------------------------------------------------
!
  CALL initremapopt(nx,ny,nz,nzsoil,nstyps,                            &
                      nprocx_in_cmd,nprocy_in_cmd,0,                   &
                      namelistfile,istatus)

  nlist = nlistfil
  gridtilt = wtgrdtiltopt
  IF(myproc == ROOT) THEN

    ref2d = (ref2dopt > 0)
    ref3d = (ref3dopt > 0)
    vel2d = (vel2dopt > 0)
    vel3d = (vel3dopt > 0)
    zdr2d = (zdr2dopt > 0)
    zdr3d = (zdr3dopt > 0)
    kdp2d = (kdp2dopt > 0)
    kdp3d = (kdp3dopt > 0)
    rhv2d = (rhv2dopt > 0)
    rhv3d = (rhv3dopt > 0)
    velproc = (velprocopt > 0)
    unfdiag = (unfdiagopt > 0)
    qc_on = (qcopt > 0)
    traditional = (tradopt > 0)
    vad = (vadopt > 0)
    fntime = (fntimopt > 0)

    WRITE(6,'(a,i4)') 'Obtained from namelists: dmpfmt=',dmpfmt
    WRITE(6,'(26x,a,i4)') 'hdf4cmpr=',hdf4cmpr

    nncradfn = 0
    DO lfile=1,nlistfil
      kfile = 0
      WRITE(stdout,'(a,a)') 'Reading file lists ',TRIM(listfil(lfile))
      OPEN(31,file=TRIM(listfil(lfile)),status='old',                     &
              form='formatted',iostat=istat)
      IF(istat == 0) THEN
        DO ifile=1,maxncradfile
          READ(31,'(a)',IOSTAT=istat) tempstr
          IF(istat /= 0) EXIT
          kfile=kfile+1
          ncradfn(ifile,lfile) = ADJUSTL(tempstr)
          WRITE(stdout,'(a,i4,a,a)')                                      &
                ' kfile:',kfile,' ncradfn:',TRIM(ncradfn(ifile,lfile))
        END DO
        nncradfn(lfile)=kfile
        WRITE(stdout,'(i4,a)') nncradfn(lfile),' file names read'
      ELSE
        WRITE(stdout,'(a,a,i6)')                                          &
              'Error opening ',TRIM(listfil(lfile)),istat
      END IF
    END DO
!
!------------------------------------------------------------------------
!
! Process command line
! For backward comptability allow input of radar file names via
! files specified on the command line
!
!------------------------------------------------------------------------
!
    WRITE(6,'(a,i3)') 'n arguments: ',narg
    IF(narg > 0) THEN
      nlist = 0
      CALL GET_COMMAND_ARGUMENT(1, charg )
      radname=charg(1:4)
      WRITE(stdout,'(a,a)')  '    radname = ', radname
      kfile=0
      nncradfn=0
      iarg=2
      DO jarg=2,narg
        CALL GET_COMMAND_ARGUMENT(iarg, charg )
        IF(charg(1:6) == '-novel') THEN
          velproc=.FALSE.
        ELSE IF(charg(1:4) == '-hdf') THEN
          dmpfmt=3
          hdf4cmpr=0
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i1)') hdf4cmpr
            hdf4cmpr=min(max(hdf4cmpr,0),7)
          END IF
          WRITE(stdout,'(a,i2)')                                            &
          ' Output in hdf format with compression level: ',hdf4cmpr
        ELSE IF(charg(1:7) == '-binary') THEN
          dmpfmt=1
          hdf4cmpr=0
          WRITE(stdout,'(a)') ' Output in binary format'
        ELSE IF(charg(1:6) == '-ref2d') THEN
          ref2d=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 2d reflectivity file of lowest tilt'
        ELSE IF(charg(1:6) == '-ref3d') THEN
          ref3d=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d reflectivity file for plotting'
        ELSE IF(charg(1:8) == '-reffile') THEN
          ref3d=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d reflectivity file for plotting'
        ELSE IF(charg(1:6) == '-vel2d') THEN
          vel2d=.TRUE.
          velproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 2d velocity file of lowest tilt'
        ELSE IF(charg(1:6) == '-vel3d') THEN
          vel3d=.TRUE.
          velproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d velocity file for plotting'
        ELSE IF(charg(1:8) == '-velfile') THEN
          vel3d=.TRUE.
          velproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d velocity file for plotting'
        ELSE IF(charg(1:6) == '-zdr2d') THEN
          zdr2d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 2d Zdr file of lowest tilt'
        ELSE IF(charg(1:6) == '-zdr3d') THEN
          zdr3d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d Zdr file for plotting'
        ELSE IF(charg(1:6) == '-kdp2d') THEN
          kdp2d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 2d Kdp file of lowest tilt'
        ELSE IF(charg(1:6) == '-kdp3d') THEN
          kdp3d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d Kdp file for plotting'
        ELSE IF(charg(1:6) == '-rhv2d') THEN
          rhv2d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 2d Rho-HV file of lowest tilt'
        ELSE IF(charg(1:6) == '-rhv3d') THEN
          rhv3d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Will produce 3d Rho-HV file for plotting'
        ELSE IF(charg(1:8) == '-unfdiag') THEN
          unfdiag=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Write files for velocity unfolding diagnostics'
        ELSE IF(charg(1:5) == '-noqc') THEN
          qc_on=.FALSE.
          WRITE(stdout,'(a)')                                               &
            ' Will skip qc steps in processing'
        ELSE IF(charg(1:8) == '-fillref') THEN
          fillref=1
          WRITE(stdout,'(a)')                                               &
            ' Will fill-in reflectivity below lowest tilt'
        ELSE IF(charg(1:8) == '-medfilt') THEN
          medfilt=1
          WRITE(stdout,'(a)')                                               &
            ' Apply median filter to az-ran data in despeckle step.'
        ELSE IF(charg(1:10) == '-nprocx_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i3)') nprocx_in
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (X) = ',nprocx_in

        ELSE IF(charg(1:10) == '-nprocy_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i3)') nprocy_in
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (Y) = ',nprocy_in

        ELSE IF(charg(1:6) == '-tradi') THEN
          traditional = .TRUE.
          WRITE(stdout,'(a)') ' Use traditional backgroud processing (initgrdvar).'
        ELSE IF(charg(1:7) == '-fntime') THEN
          fntime=.TRUE.
          WRITE(stdout,'(a)') ' Use time in filename for output filename.'

        ELSE IF(charg(1:7) == '-rtopts') THEN
          rtopts=.TRUE.
          WRITE(stdout,'(a)') ' Using real-time options.'
        ELSE IF(charg(1:9) == '-gridtilt') THEN
          gridtilt=1
          WRITE(stdout,'(a)') ''
          WRITE(stdout,'(a)')                                               &
          ' Output will be in EnKF format (on radar elevation).'

        ELSE IF(charg(1:8) == '-reftime') THEN
          reftime = 1
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i4,2i2,a,2i2)') ref_yr,ref_mon,ref_day,ref_garb,   &
                                     ref_hr,ref_min
          END IF
          WRITE(stdout,'(a)') ''
          WRITE(stdout,'(a,i4,a,4(i2,a))') ' Output is valid at ',          &
              ref_yr,'/',ref_mon,'/',ref_day,' ',ref_hr,':',ref_min,' UTC'

          CALL ctim2abss (ref_yr,ref_mon,ref_day,ref_hr,ref_min,0,ref_t)

        ELSE
          nlist = nlist+1
          kfile = 0
          full_fname=charg
          WRITE(stdout,'(a,a)')  ' Reading file lists ',TRIM(full_fname)
          OPEN(31,file=full_fname,status='old',form='formatted',iostat=istat)
          IF(istat == 0) THEN
            DO ifile=1,maxncradfile
              READ(31,'(a)',IOSTAT=istat) tempstr
              IF(istat /= 0) EXIT
              ncradfn(ifile,nlist) = ADJUSTL(tempstr)
              kfile=kfile+1
              WRITE(stdout,'(a,i4,a,a)')                                      &
                 ' kfile:',kfile,' ncradfn: *',TRIM(ncradfn(ifile,nlist))
            END DO
            nncradfn(nlist)=kfile
            WRITE(stdout,'(i4,a)') nncradfn(nlist),' file names read'
          ELSE
            WRITE(stdout,'(a,a,i6)') 'Error opening ',TRIM(full_fname),istat
          END IF
          nlistfil = nlist
        END IF
        iarg=iarg+1
        IF(iarg > narg) EXIT
      END DO ! process arguments loop

    END IF  ! narg > 0

    IF(tintrpopt > 0) THEN
      IF(nlistfil /= 2 .OR. gridtilt /= 1 .or. reftime /= 1) THEN
        WRITE(stdout,'(/,a)') 'Time interpolation option > 0'
        WRITE(stdout,'(/,a)') 'Use two list files and'
        WRITE(stdout,'(a)') ' use -gridtilt & -reftime options.'
        STOP
     END IF

   END IF

  END IF

  IF(nncradfn(1) == 0) THEN
    nlistfil=1
    kfile=1
    nncradfn(kfile)=1
    ncradfn(kfile,1)=TRIM(radfname)
    WRITE(stdout,'(a,i4,a,a)')                                        &
              ' kfile:',kfile,' ncradfn:',TRIM(ncradfn(kfile,1))
  END IF

  CALL mpupdatei(nlistfil,1)
  CALL mpupdatec(radname,4)
  CALL mpupdatel(velproc,1)
  CALL mpupdatei(dmpfmt,1)
  CALL mpupdatei(hdf4cmpr,1)
  CALL mpupdatei(dualpout,1)
  CALL mpupdatel(ref2d,1)
  CALL mpupdatel(ref3d,1)
  CALL mpupdatel(vel2d,1)
  CALL mpupdatel(vel3d,1)
  CALL mpupdatel(zdr2d,1)
  CALL mpupdatel(zdr3d,1)
  CALL mpupdatel(kdp2d,1)
  CALL mpupdatel(kdp3d,1)
  CALL mpupdatel(rhv2d,1)
  CALL mpupdatel(rhv3d,1)
  CALL mpupdatel(unfdiag,1)
  CALL mpupdatel(qc_on,1)
  CALL mpupdatei(nncradfn,maxlistf)
  CALL mpupdatec(ncradfn,256*maxncradfile*maxlistf)
  CALL mpupdatei(fillref,1)
  CALL mpupdatel(traditional,1)
  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)

  istatus = 0
  klev = 0
  rfirstg = 0.
  gatesp  = 100.
  runname_in = runname

  irfirstg = NINT(rfirstg)
  igatesp  = NINT(gatesp)
  refrngmin = rngmin
  refrngmax = rngmax

  IF(dualpout > 0) dualpproc = .TRUE.

!------------------------------------------------------------------------
!
! Loop through the NetCDF radar files.
!
!------------------------------------------------------------------------
  maxinfile = 0
  DO ilist = 1,nlistfil
    IF (nncradfn(ilist) > 0 .AND. nncradfn(ilist) <= maxncradfile) THEN
      maxinfile = max(maxinfile,nncradfn(ilist))
      IF(myproc == ROOT ) WRITE(stdout,'(a,i4,a,i4)')       &
        '  List number: ',ilist,'  nfiles: ',nncradfn(ilist)
    ELSE IF (nncradfn(ilist) > maxncradfile) THEN
      IF(myproc == ROOT ) WRITE(stdout,'(a,i4,a,i4,/,a,i4,a)')       &
        ' WARNING: nncradfn=', nncradfn(ilist),' > maxncradfile=',maxncradfile,  &
        ' Will only read in ', maxncradfile,' netCDF files.'
      nncradfn(ilist)=maxncradfile
      maxinfile = max(maxinfile,nncradfn(ilist))
    END IF
  END DO

  IF (maxinfile > 0) THEN
!
!------------------------------------------------------------------------
!
!   First check each file for variable stored and dimensions of data
!
!------------------------------------------------------------------------
!
    maxgate=1
    maxazim=1
    maxelev=1
    mintimcdf=-1
    DO ilist = 1,nlistfil
      DO infile = 1,maxinfile
        fvartype(infile,ilist) = 0
        gzipped(infile,ilist)  = .FALSE.
      END DO
    END DO

    IF (myproc == ROOT) THEN
      DO ilist = 1,nlistfil
        DO infile = 1,nncradfn(ilist)
          nazim = 0
          ngate = 0
          nvar  = 0
          infilname = ncradfn(infile,ilist)
          lenfn = LEN_TRIM(infilname)
          IF(fntime) THEN
            iloc=index(TRIM(infilname),'.netcdf',back)
            IF(iloc > 0) THEN
              print *, ' Getting time from string: ',                   &
                         infilname((iloc-15):(iloc-1))
              READ(infilname((iloc-15):(iloc-1)),'(i4,2i2,1x,3i2)')     &
                  kyear,kmonth,kday,khr,kmin,ksec
            ELSE
              kyear=0
              kmonth=0
              kday=0
              khr=0
              kmin=0
              ksec=0
            END IF
          END IF
          IF(infilname((lenfn-2):lenfn) == ".gz") THEN
            gzipped(infile,ilist)=.TRUE.
            WRITE(syscall,'(a,1x,a)') 'gunzip',TRIM(infilname)
            CALL system(syscall)
            infilname((lenfn-2):lenfn)='   '
            ncradfn(infile,ilist)=infilname
            WRITE(6,'(a,a)') 'Uncompressed ',TRIM(infilname)
          END IF
          WRITE(stdout,'(a,a)') 'Checking dims in file:',TRIM(infilname)
          CALL extdirn(infilname,indir,fname)
          print *, ' infile: ',TRIM(infilname)
          print *, '  indir: ',TRIM(indir)
          print *, '  fname: ',TRIM(fname)
          CALL get_ncraddims(fname,indir,iradtype,ngate,nazim,nvar,istatus)
          WRITE(6,'(a,i6,a,i6)') ' istatus: ',istatus,'  nvar in file: ',nvar
          isource=1+(iradtype/10)
          IF(istatus == NF_NOERR .AND. nvar > 0) THEN
            ALLOCATE(ncvarname(nvar), stat = istat)

            IF(iradtype == 21) THEN  ! CASA Tier-2a
              vadsrc = 'CASA VAD'
              kntt2a=kntt2a+1
              ALLOCATE(tem_double(nazim), stat = istat)
              CALL check_alloc_status(istat,'ncrad2arps:tem_double')
              CALL get_ncrad2ainfo(fname,indir,iradtype,nazim,nvar,tem_double,  &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,elv,         &
                     ncvarname,istatus)
              DEALLOCATE(tem_double)
              fvartype(infile,ilist)=100
              rdrknt(1,ilist)=rdrknt(1,ilist)+1
              rdrknt(2,ilist)=rdrknt(2,ilist)+1
              rdrknt(4,ilist)=rdrknt(4,ilist)+1
              rdrknt(6,ilist)=rdrknt(6,ilist)+1
              rdrknt(7,ilist)=rdrknt(7,ilist)+1
            ELSE IF(iradtype == 22) THEN  ! CASA Tier-2b
              vadsrc = 'CASA VAD'
              kntt2b=kntt2b+1
              CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,    &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,    &
                     ivcp,elv,ncvarname,istatus)
              DO k=1,nvar
                print *, '  variable(',k,') = ',TRIM(ncvarname(k))
                IF(TRIM(ncvarname(k)) == 'Reflect' .OR.                &
                   TRIM(ncvarname(k)) == TRIM(refvarname)) THEN
                  fvartype(infile,ilist)=2
                  rdrknt(2,ilist)=rdrknt(2,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'RawVelocity' .OR.       &
                        TRIM(ncvarname(k)) == 'AliasedVelocity' .OR.   &
                        TRIM(ncvarname(k)) == TRIM(velvarname)) THEN
                  fvartype(infile,ilist)=1
                  rdrknt(1,ilist)=rdrknt(1,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'DifferentialReflectivity') THEN
                  fvartype(infile,ilist)=4
                  rdrknt(4,ilist)=rdrknt(4,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'DifferentialPhase') THEN
                  fvartype(infile,ilist)=6
                  rdrknt(6,ilist)=rdrknt(6,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'RhoHV' .OR.             &
                        TRIM(ncvarname(k)) == TRIM(rhvvarname)) THEN
                  fvartype(infile,ilist)=7
                  rdrknt(7,ilist)=rdrknt(7,ilist)+1
                  EXIT
                END IF
              END DO
            ELSE IF(iradtype == 31) THEN  ! OU-PRIME Teir-2a style
              vadsrc = 'OUPRMVAD'
              kntt2a=kntt2a+1
              ALLOCATE(tem_double(nazim), stat= istat)
              CALL check_alloc_status(istat,'ncrad2arps:tem_double')
              CALL get_ncrad2ainfo(fname,indir,iradtype,nazim,nvar,tem_double,  &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,elv,         &
                     ncvarname,istatus)
              DEALLOCATE(tem_double)
              fvartype(infile,ilist)=100
              rdrknt(2,ilist)=rdrknt(2,ilist)+1
              rdrknt(1,ilist)=rdrknt(1,ilist)+1
              rdrknt(4,ilist)=rdrknt(4,ilist)+1
              rdrknt(6,ilist)=rdrknt(6,ilist)+1
              rdrknt(7,ilist)=rdrknt(7,ilist)+1
            ELSE IF(iradtype == 32) THEN  ! OU-PRIME EEC radar Tier-2b style
              vadsrc = 'OUPRMVAD'
              kntt2b=kntt2b+1
              CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,    &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,    &
                     ivcp,elv,ncvarname,istatus)
              DO k=1,nvar
                print *, '  variable(',k,') = ',TRIM(ncvarname(k))
                IF(TRIM(ncvarname(k)) == 'Corrected_Intensity' .OR.    &
                   TRIM(ncvarname(k)) == TRIM(refvarname)) THEN
                  fvartype(infile,ilist)=2
                  rdrknt(2,ilist)=rdrknt(2,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'Radial_Velocity' .OR.   &
                        TRIM(ncvarname(k)) == TRIM(velvarname)) THEN
                  fvartype(infile,ilist)=1
                  rdrknt(1,ilist)=rdrknt(1,ilist)+1
                  EXIT
               ELSE IF(TRIM(ncvarname(k)) == 'DifferentialReflectivity') THEN
                  fvartype(infile,ilist)=4
                  rdrknt(4,ilist)=rdrknt(4,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'DifferentialPhase') THEN
                  fvartype(infile,ilist)=6
                  rdrknt(6,ilist)=rdrknt(6,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'RhoHV' .OR.             &
                        TRIM(ncvarname(k)) == TRIM(rhvvarname)) THEN
                  fvartype(infile,ilist)=7
                  rdrknt(7,ilist)=rdrknt(7,ilist)+1
                  EXIT
                END IF
              END DO
            ELSE IF(iradtype == 42) THEN   ! KOUN Nexrad radar
              vadsrc = 'KOUN VAD'
              kntt2b=kntt2b+1
              CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,    &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,    &
                     ivcp,elv,ncvarname,istatus)
              DO k=1,nvar
                print *, '  variable(',k,') = ',TRIM(ncvarname(k))
                IF(TRIM(ncvarname(k)) == 'Reflectivity' .OR.           &
                   TRIM(ncvarname(k)) == TRIM(refvarname)) THEN
                  fvartype(infile,ilist)=2
                  rdrknt(2,ilist)=rdrknt(2,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'Velocity' .OR.          &
                        TRIM(ncvarname(k)) == TRIM(velvarname)) THEN
                  fvartype(infile,ilist)=1
                  rdrknt(1,ilist)=rdrknt(1,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'RhoHV' .OR.             &
                        TRIM(ncvarname(k)) == TRIM(rhvvarname)) THEN
                  fvartype(infile,ilist)=7
                  rdrknt(7,ilist)=rdrknt(7,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'Zdr' .OR.               &
                        TRIM(ncvarname(k)) == TRIM(zdrvarname)) THEN
                  fvartype(infile,ilist)=4
                  rdrknt(4,ilist)=rdrknt(4,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'Kdp' .OR.               &
                        TRIM(ncvarname(k)) == TRIM(kdpvarname)) THEN
                  fvartype(infile,ilist)=6
                  rdrknt(6,ilist)=rdrknt(6,ilist)+1
                  EXIT
                END IF
              END DO
            ELSE IF(iradtype == 51) THEN  ! Sigmet IRIS FSL Format
              vadsrc='IRIS VAD'
              kntt2a=kntt2a+1
              ALLOCATE(tem_double(nazim))
              CALL get_ncrad2ainfo(fname,indir,iradtype,nazim,nvar,tem_double,  &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,elv,         &
                     ncvarname,istatus)
              DEALLOCATE(tem_double)
              fvartype(infile,ilist)=100
              rdrknt(2,ilist)=rdrknt(2,ilist)+1
              rdrknt(1,ilist)=rdrknt(1,ilist)+1
            ELSE IF(iradtype == 52) THEN   ! CWB or TDWR radar
              vadsrc = 'CWB VAD'
              kntt2b=kntt2b+1
              CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,     &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,     &
                     ivcp,elv,ncvarname,istatus)
              IF(radarname(4:4) == ' ') THEN
                vadsrc = 'TDWR VAD'
                DO i=3,1,-1
                  radarname(i+1:i+1)=radarname(i:i)
                END DO
                radarname(1:1)='T'
              END IF
              DO k=1,nvar
                print *, '  variable(',k,') = ',TRIM(ncvarname(k))
                IF(TRIM(ncvarname(k)) == 'Reflectivity' .OR.           &
                   TRIM(ncvarname(k)) == TRIM(refvarname)) THEN
                  fvartype(infile,ilist)=2
                  rdrknt(2,ilist)=rdrknt(2,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'Velocity' .OR.          &
                        TRIM(ncvarname(k)) == TRIM(velvarname)) THEN
                  fvartype(infile,ilist)=1
                  rdrknt(1,ilist)=rdrknt(1,ilist)+1
                  EXIT
                END IF
              END DO
            ELSE IF (iradtype == 61) THEN ! Foray processed data
              vadsrc = 'FIELDVAD'
              kntt2a = kntt2a+1
              ALLOCATE(tem_double(nazim),stat=istat)
              CALL check_alloc_status(istat, 'ncrad2arps:tem_double')
              CALL get_ncrad2ainfo(fname,indir,iradtype,nazim,nvar,tem_double,   &
                      radarname,radlat,radlon,radalt,itimcdf,frtime,elv,         &
                      ncvarname,istatus)
              DEALLOCATE(tem_double)
              fvartype(infile,ilist) = 100
              rdrknt(2,ilist) = rdrknt(2,ilist)+1
              rdrknt(1,ilist) = rdrknt(1,ilist)+1
              rdrknt(4,ilist) = rdrknt(4,ilist)+1
              rdrknt(6,ilist) = rdrknt(6,ilist)+1
              rdrknt(7,ilist) = rdrknt(7,ilist)+1
            ELSE IF(iradtype == 72) THEN  ! TDWR NIDS converted to netCDF
              vadsrc = 'TDWR VAD'
              kntt2b=kntt2b+1
              CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,    &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,    &
                     ivcp,elv,ncvarname,istatus)
              DO k=1,nvar
                print *, '  variable(',k,') = ',TRIM(ncvarname(k))
                IF(TRIM(ncvarname(k)) == 'BaseReflectivity' .OR.       &
                   TRIM(ncvarname(k)) == TRIM(refvarname)) THEN
                  fvartype(infile,ilist)=2
                  rdrknt(2,ilist)=rdrknt(2,ilist)+1
                  EXIT
                ELSE IF(TRIM(ncvarname(k)) == 'RadialVelocity' .OR.    &
                        TRIM(ncvarname(k)) == TRIM(velvarname)) THEN
                  fvartype(infile,ilist)=1
                  rdrknt(1,ilist)=rdrknt(1,ilist)+1
                  EXIT
                END IF
              END DO
            ELSE
              WRITE(stdout,'(a,a)') 'Unknown radar file type:',iradtype
              fvartype(infile,ilist)=0
!              STOP
              istatus = -1
              EXIT
            END IF

            !
            ! Altitude and elevation corrections for a few old data files
            !
            IF(radarname(1:5) == 'cyril') THEN
              IF(radalt < 100.) radalt=432.9
              IF(elv < 0.1) elv=1.0
            ELSE IF(radarname(1:5) == 'chick') THEN
              IF(radalt < 100.) radalt=354.0
              IF(elv < 0.1) elv=1.0
            ELSE IF(radarname(1:5) == 'lawto') THEN
              IF(radalt < 100.) radalt=377.4
              IF(elv < 0.1) elv=1.0
            ELSE IF(radarname(1:5) == 'rushs') THEN
              IF(radalt < 100.) radalt=414.8
              IF(elv < 0.1) elv=1.0
            END IF

            felv(infile,ilist)=elv
            fitime(infile,ilist)=itimcdf

            scantime(infile,ilist)=itimcdf+itime1970

            maxazim=max(nazim,maxazim)
            maxgate=max(ngate,maxgate)

            IF(mintimcdf < 0) THEN
              mintimcdf=itimcdf
              kmyear=kyear
              kmmon=kmonth
              kmday=kday
              kmhr=khr
              kmmin=kmin
              kmsec=ksec
            ELSE IF( itimcdf < mintimcdf) THEN
              mintimcdf=itimcdf
              kmyear=kyear
              kmmon=kmonth
              kmday=kday
              kmhr=khr
              kmmin=kmin
              kmsec=ksec
            END IF

            DEALLOCATE(ncvarname)

          ELSE  ! error return from get dims

            fvartype(infile,ilist)=0

          END IF

        END DO
      END DO
    END IF

    CALL mpupdatei(istatus,1)
    IF (istatus < 0) CALL arpsstop('  ',1)

    CALL mpupdatec(radarname,32)
    CALL mpupdater(radlat,1)
    CALL mpupdater(radlon,1)
    CALL mpupdater(radalt,1)

    CALL mpupdatei(rdrknt,numvar*maxlistf)

    rdrknt_tot(:) = 0
    IF( tintrpopt > 0) THEN
      DO i=1,numvar
        rdrknt_tot(i)=rdrknt(i,1)
      END DO
    ELSE
      DO lfile=1,nlistfil
        DO i=1,numvar
          rdrknt_tot(i)=rdrknt_tot(i)+rdrknt(i,lfile)
        END DO
      END DO
    END IF

    maxelev=maxval(rdrknt_tot)
    CALL mpupdatei(rdrknt_tot,numvar)

    CALL mpupdatei(maxazim,1)
    CALL mpupdatei(maxgate,1)

    IF(rdrknt_tot(1) == 0) THEN
      WRITE(6,'(a,/a,6x,i6//)')                                        &
       ' No velocity files found so turning off velocity processing.', &
       'velknt=',rdrknt_tot(1)
      velproc = .FALSE.
    END IF

    IF(rdrknt_tot(4) == 0 .AND. rdrknt_tot(6) == 0                     &
                          .AND. rdrknt_tot(7) == 0) THEN
      WRITE(6,'(a//)')                                                 &
       ' No Zdr, Kdr nor RhoHV data found so turning off dual-pol processing.'
      dualpproc = .FALSE.
      dualpout = 0
    END IF

    IF (myproc == ROOT) THEN
      DO idp = 1, numvar
        IF(rdrknt_tot(idp) > 0) THEN
          WRITE(6,'(a,i6,2x,a,a,i6)')  ' Number of',idp,               &
                dualvarname(idp),' files:',rdrknt_tot(idp)
        ENDIF
      ENDDO
      WRITE(6,'(/a,i6)')  ' Max azimuth dimension of NetCDF data:',maxazim
      WRITE(6,'( a,i6)')  ' Max gate dimension of NetCDF data:',maxgate
      WRITE(6,'( a,i6/)') ' Max elevations among data:',maxelev
    END IF

    CALL mpupdatei(kntt2a,1)
    CALL mpupdatei(kntt2b,1)

    IF(kntt2a > 0 .AND. kntt2b > 0) THEN
      IF (myproc == ROOT) WRITE(6,'(a)') ' Mixed Tier 2a and Tier 2b files in list, stopping'
      IF (myproc == ROOT) WRITE(6,'(a,i4,a,i4)')' N Tier 2a=',kntt2a,'  N Tier 2b=',kntt2b
      CALL arpsstop(' ',1)
    END IF

    ALLOCATE(x(nx),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:x')
    ALLOCATE(y(ny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:y')
    ALLOCATE(z(nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:z')
    ALLOCATE(zp(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:zp')

    ALLOCATE(xs(nx),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:xs')
    ALLOCATE(ys(ny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:ys')
    ALLOCATE(zps(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:zps')

    ALLOCATE(tem2d(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:tem2d')

    ALLOCATE(j1(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:j1')
    ALLOCATE(j2(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:j2')
    ALLOCATE(j3(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:j3')
    ALLOCATE(hterain(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:hterain')
    ALLOCATE(mapfct(nx,ny,8),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:mapfct')
    ALLOCATE(tem3d1(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:tem3d1')

    ALLOCATE(tem1d1(nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:tem1d1')
    ALLOCATE(tem1d2(nz),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:tem1d2')
    ALLOCATE(rdr_p(maxgate,numrad2,maxelev,numvar))
    ALLOCATE(rdr_tmp(maxgate,numrad,nlistfil))
    rdr_p = missing
    rdr_tmp = missing
    ALLOCATE(elvang(maxelev))

    IF(qc_on .AND. rdrknt_tot(1) > 0) THEN
!------------------------------------------------------------------------
!
! Fill refractivity sounding with constant lapse rate
! This is only used when rfropt > 1, but for completeness and
! future upgrade (when an actual sounding made from gridded data
! would be used) this is filled-in.
!
!------------------------------------------------------------------------

      DO k=1,nzsnd
        zsnd(k)=(k-1)*dzsnd
        rfrsnd(k)=325.+rfrconst*zsnd(k)
      END DO

!------------------------------------------------------------------------
!
! ARPS grid array allocations and variable initializations.
!
!------------------------------------------------------------------------

      ALLOCATE(u    (nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:u')
      ALLOCATE(v    (nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:v')
      ALLOCATE(ptprt(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:ptprt')
      ALLOCATE(pprt (nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:pprt')
      ALLOCATE(qv   (nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:qv')

      ALLOCATE(tem3d2(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:tem3d2')
      ALLOCATE(tem3d3(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:tem3d3')

      IF ((.NOT. traditional) .AND.                                    &
          (initopt == 3 .AND. (inifmt == 3 .OR. inifmt == 1)) ) THEN

        CALL get_gridxyzzp(nx,ny,nz,inigbf,inifmt,nprocx_in,nprocy_in, &
                           x,y,z,zp,istatus);

        CALL mapinit(nx,ny);

        IF (myproc == ROOT)  print *, ' Setting radar coordinate: '

        CALL radcoord(nx,ny,nz,x,y,z,zp,xs,ys,zps,                     &
                      radlat,radlon,radarx,radary);

        IF (myproc == ROOT) THEN
          print *, ' Radar x: ',(0.001*radarx),' km'
          print *, ' Radar y: ',(0.001*radary),' km'
          CALL flush(6)
        END IF

        CALL readuvt(nx,ny,nz,inifile,inifmt,nprocx_in,nprocy_in,      &
                     envtime,u,v,pprt,ptprt,qv,istatus)

        IF (myproc == ROOT) WRITE(6,'(1x,a,f10.2)')                    &
           'Environmental wind averaging radius: ',envavgr

        CALL extenvprf2(nx,ny,nz,nzsnd,x,y,zp,xs,ys,zps,               &
             u,v,ptprt,pprt,qv,tem3d1,tem3d2,tem3d3,                   &
             radarx,radary,radalt,envavgr,rngmax,                      &
             zsnd,ktsnd,usnd,vsnd,rfrsnd,istatus);

        IF (myproc == ROOT) WRITE(6,'(1x,a)') 'Back from extenvprf'

        CALL mpupdatei(istatus, 1)

        IF(istatus < 0) THEN
          IF (myproc == ROOT) WRITE(6,'(a,i3,/a)')                     &
                        ' Status from extenvprf = ',istatus,' STOPPING'
          CALL arpsstop('arpsstop called from ncrad2arps',1)
        END IF

      ELSE

        ALLOCATE(ubar(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:ubar')
        ALLOCATE(vbar(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:vbar')
        ALLOCATE(ptbar(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:ptbar')
        ALLOCATE(pbar(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:pbar')
        ALLOCATE(qvbar(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:qvbar')
        ALLOCATE(rhostr(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:rhostr')

        ALLOCATE(tem2dyz(ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem2dyz')
        ALLOCATE(tem2dxz(nx,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem2dxz')
        ALLOCATE(tem2dns(nx,ny,0:nstyps),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem2dns')

        ALLOCATE(trigs1(3*(nx-1)/2+1),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:trigs1')
        ALLOCATE(trigs2(3*(ny-1)/2+1),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:trigs2')
        ALLOCATE(ifax(13),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:ifax')

        ALLOCATE(wsave1(3*(ny-1)+15),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:wsave1')
        ALLOCATE(wsave2(3*(nx-1)+15),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:wsave2')
        ALLOCATE(vwork1(nx+1,ny+1),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:vwork1')
        ALLOCATE(vwork2(ny,nx+1),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:vwork2')

        ALLOCATE(qcumsrc(nx,ny,nz,5),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:qcumsrc')
        ALLOCATE(prcrate(nx,ny,4),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:prcrate')
        ALLOCATE(exbcbuf(exbcbufsz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:exbcbufsz')
        ALLOCATE(soiltyp(nx,ny,nstyps),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:soiltyp')
        ALLOCATE(stypfrct(nx,ny,nstyps),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:stypfrct')
        ALLOCATE(tem2dint(nx,ny),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem2dint')

        ALLOCATE(tem3d4(nx,ny,nz),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem3d4')

        ALLOCATE(tem3dsoil(nx,ny,nzsoil),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem3dsoil')
        ALLOCATE(tem4dsoilns(nx,ny,nzsoil,0:nstyps),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:tem4dsoilns')
        ALLOCATE(qscalar(nx,ny,nz,nscalar),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:qscalar')

        CALL set_lbcopt(1)

        CALL initgrdvar(nx,ny,nz,nzsoil,nt,nstyps,exbcbufsz,              &
                x,y,z,zp,tem3dsoil,hterain,mapfct,                        &
                j1,j2,j3,tem3dsoil,tem3d1,tem3d1,tem3d1,tem3d1,tem3dsoil, &
                u,v,tem3d1,tem3d1,ptprt,pprt,                             &
                qv,qscalar,tem3d1,                                        &
                tem2dyz,tem2dyz,tem2dxz,tem2dxz,                          &
                tem2dyz,tem2dyz,tem2dxz,tem2dxz,                          &
                trigs1,trigs2,ifax,ifax,                                  &
                wsave1,wsave2,vwork1,vwork2,                              &
                ubar,vbar,tem3d1,ptbar,pbar,tem3d1,tem3d1,                &
                rhostr,tem3d1,qvbar,tem3d1,tem3d1,                        &
                soiltyp,stypfrct,tem2dint,tem2d,tem2d,tem2d,              &
                tem4dsoilns,tem4dsoilns,tem2dns,tem2d,tem2dns,            &
                tem3d1,qcumsrc,tem3d1,tem2dint,tem2d,                     &
                tem2d,tem2d,tem2d,                                        &
                tem2d,tem2d,prcrate,exbcbuf,                              &
                tem3d1,tem2d,tem2d,tem2d,tem2d,                           &
                tem2d,tem2d,tem2d,tem2d,                                  &
                tem3dsoil,tem3dsoil,tem3dsoil,tem3dsoil,tem3dsoil,        &
                tem3d2,tem3d3,tem3d4,tem3d4,tem3d4,                       &
                tem3d4,tem3d4,tem3d4,tem3d4)

        IF (myproc == ROOT) print *, ' Back from initgrdvar'

        DEALLOCATE(trigs1,trigs2,ifax,wsave1,wsave2)
        DEALLOCATE(vwork1,vwork2,qcumsrc,prcrate)
        DEALLOCATE(soiltyp,stypfrct,exbcbuf)
        DEALLOCATE(tem2dint,tem2dyz,tem2dxz,tem2dns)
        DEALLOCATE(tem4dsoilns,tem3dsoil)
        DEALLOCATE(qscalar)

        IF (myproc == ROOT) print *, ' Setting radar coordinate: '

        CALL radcoord(nx,ny,nz,x,y,z,zp,xs,ys,zps, &
                  radlat,radlon,radarx,radary)

        IF (myproc == ROOT) THEN
          print *, ' Radar x: ',(0.001*radarx),' km'
          print *, ' Radar y: ',(0.001*radary),' km'
        END IF

        IF (myproc == ROOT)  &
          print *, ' Environmental wind averaging radius: ',envavgr


        CALL extenvprf(nx,ny,nz,nzsnd,    &
                   x,y,zp,xs,ys,zps,      &
                   u,v,ptprt,pprt,qv,ptbar,pbar,tem3d1,tem3d2,tem3d3, &
                   radarx,radary,radalt,envavgr,rngmax, &
                   zsnd,ktsnd,usnd,vsnd,rfrsnd,istatus)

        print *, ' Back from extenvprf'

        IF(istatus < 0) THEN
          IF (myproc == ROOT) WRITE(6,'(a,i3,/a)')     &
             ' Status from extenvprf = ',istatus,' STOPPING'
          CALL arpsstop('arpsstop called from ncrad2arps',1)
        END IF

        DEALLOCATE(ubar, vbar, pbar, ptbar, qvbar, rhostr)
        DEALLOCATE(tem3d4)
      END IF

      DEALLOCATE(u, v, pprt, ptprt, qv)
      DEALLOCATE(tem3d2, tem3d3)

    ELSE

      ALLOCATE(zpsoil(nx,ny,nzsoil),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:tem4dsoilns')
      ALLOCATE(j3soil(nx,ny,nzsoil),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:tem4dsoilns')
      ALLOCATE(j3soilinv(nx,ny,nzsoil),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:tem4dsoilns')

      IF (myproc == ROOT) print *, ' Setting ARPS grid coordinates ...'
      CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                      &
                  hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,             &
                  tem1d1,tem1d2,tem3d1)

      IF (myproc == ROOT) print *, ' Setting radar coordinate...'

      CALL radcoord(nx,ny,nz,x,y,z,zp,xs,ys,zps,                        &
                    radlat,radlon,radarx,radary)

      DEALLOCATE(zpsoil, j3soil, j3soilinv)

    END IF

    DEALLOCATE(j1, j2, j3, hterain, mapfct, tem3d1)
    DEALLOCATE(tem1d1, tem1d2)

!-----------------------------------------------------------------------
!
! Allocation of analysis arrays
!
!-----------------------------------------------------------------------

    nxny=nx*ny
    ALLOCATE(icolp(nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:icolp')
    ALLOCATE(jcolp(nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:jcolp')
    ALLOCATE(xcolp(nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:xcolp')
    ALLOCATE(ycolp(nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:ycolp')
    ALLOCATE(zcolp(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:zcolp')
    ALLOCATE(havdat(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:havdat')

    ALLOCATE(kntbin(nsort),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:kntbin')

    ALLOCATE(colref(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:colref')
    ALLOCATE(colvel(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:colvel')
    ALLOCATE(colnyq(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:colnyq')
    ALLOCATE(coltim(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:coltim')
    IF (rhv3d .OR. dualpout ) THEN
      ALLOCATE(colrhv(nz,nxny),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:colrhv')
    END IF
    IF (zdr3d .OR. dualpout ) THEN
      ALLOCATE(colzdr(nz,nxny),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:colzdr')
    END IF
    IF (kdp3d .OR. dualpout ) THEN
      ALLOCATE(colkdp(nz,nxny),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:colkdp')
    END IF

!------------------------------------------------------------------------
!
! Radar data array allocations
!
!------------------------------------------------------------------------

    IF(rdrknt_tot(2) > 0) THEN
      ALLOCATE(refdata(maxgate,maxazim),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:refdata')
      ALLOCATE(arfdata(maxgate,maxazim),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:arfdata')
    END IF
    IF(rdrknt_tot(1) > 0) THEN
      ALLOCATE(veldata(maxgate,maxazim),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:veldata')
      IF(qc_on) THEN
        ALLOCATE(spwdata(maxgate,maxazim),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:spwdata')
        ALLOCATE(unfvdata(maxgate,maxazim),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:unfvdata')
        ALLOCATE(bkgvel(maxgate,maxazim),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:bkgvel')
        ALLOCATE(bgate(maxazim),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:bgate')
        ALLOCATE(egate(maxazim),stat=istat)
        CALL check_alloc_status(istat,'ncrad2arps:egate')
        spwdata=0.
      END IF
    END IF

    IF(rdrknt_tot(4) > 0) THEN
      ALLOCATE(zdrdata(maxgate,maxazim),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:zdrdata')
    END IF
    IF(rdrknt_tot(6) > 0) THEN
      ALLOCATE(kdpdata(maxgate,maxazim),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:kdpdata')
    END IF
    IF(rdrknt_tot(7) > 0) THEN
      ALLOCATE(rhvdata(maxgate,maxazim),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:rhvdata')
    END IF

    ALLOCATE(snstvty(maxgate),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:snstvty')
    ALLOCATE(rtem(maxgate,maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:rtem')
    ALLOCATE(time(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:time')
    ALLOCATE(gcfst(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:gcfst')
    ALLOCATE(azim(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:azim')
    ALLOCATE(gtspc(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:gtspc')
    ALLOCATE(vnyq(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:vnyq')
    ALLOCATE(elev(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:elev')
    ALLOCATE(istrgate(maxazim),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:istrgate')

    ALLOCATE(kntrgat(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:kntrgat')
    ALLOCATE(kntrazm(maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:kntrazm')

    ALLOCATE(kntvgat(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:kntvgat')
    ALLOCATE(kntvazm(maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:kntvazm')

    ALLOCATE(vnyqvol(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:vnyqvol')
    ALLOCATE(timevolr(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:timevolr')
    ALLOCATE(timevolv(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:timevolv')

    ALLOCATE(rngrvol(maxgate,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:rngrvol')
    ALLOCATE(azmrvol(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:azmrvol')
    ALLOCATE(elvrvol(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:elvrvol')
    ALLOCATE(elvmnrvol(maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:elvmnrvol')
    ALLOCATE(refvol(maxgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:refvol')
    IF(dualpproc) THEN
      ALLOCATE(rhvvol(maxgate,maxazim,maxelev),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:rhvvol')
      ALLOCATE(zdrvol(maxgate,maxazim,maxelev),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:zdrvol')
      ALLOCATE(kdpvol(maxgate,maxazim,maxelev),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:kdpvol')
    END IF

    ALLOCATE(rngvvol(maxgate,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:rngvvol')
    ALLOCATE(azmvvol(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:azmvvol')
    ALLOCATE(elvvvol(maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:elvvol')
    ALLOCATE(elvmnvvol(maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:elvmnvvol')
    ALLOCATE(velvol(maxgate,maxazim,maxelev),stat=istat)
    print *, ' velvol status: ',istat
    CALL check_alloc_status(istat,'ncrad2arps:velvol')

    ALLOCATE(rxvol(maxgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:rxvol')
    ALLOCATE(ryvol(maxgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:ryvol')
    ALLOCATE(rzvol(maxgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'ncrad2arps:rzvol')

!------------------------------------------------------------------------
!
!     Initialize 3-d radar data arrays
!
!------------------------------------------------------------------------

    CALL mpupdatei(mintimcdf,1)

    itimcdf=mintimcdf+itime1970
    itimfrst=itimcdf

    CALL abss2ctim(itimcdf, iyear, imonth, iday, ihr, imin, isec)

    IF(tintrpopt > 0 .AND. nlistfil > 1) THEN
      READ(ref_time,'(i4,2i2,a,2i2)') ref_yr,ref_mon,ref_day,ref_garb,  &
                                      ref_hr,ref_min
      CALL ctim2abss (ref_yr,ref_mon,ref_day,ref_hr,ref_min,0,ref_t)
    ELSE
      ref_t = itimcdf
    ENDIF

    runname=runname_in
    CALL gtlfnkey(runname, lfnkey)
    curtim=outtime

    CALL rmpinit(nx,ny,nz,maxgate,maxgate,maxazim,maxelev,dualpproc,    &
            kntrgat,kntrazm,kntrelv,                                    &
            kntvgat,kntvazm,kntvelv,                                    &
            vnyqvol,timevolr,timevolv,                                  &
            rngrvol,azmrvol,elvrvol,elvmnrvol,                          &
            refvol,rhvvol,zdrvol,kdpvol,                                &
            rngvvol,azmvvol,elvvvol,elvmnvvol,velvol,                   &
            colvel,colref,colnyq,coltim)

    IF (fntime) THEN
      iiyear = kmyear
      iimon = kmmon
      iiday = kmday
      iihr = kmhr
      iimin = kmmin
      iisec = kmsec
    ELSE
      iiyear = iyear
      iimon = imonth
      iiday = iday
      iihr = ihr
      iimin = imin
      iisec = isec
    END IF

    IF (myproc == ROOT) THEN
      IF(kntt2a > 0) THEN
        DO ilist = 1,nlistfil
        DO infile = 1,nncradfn(ilist)
          infilname = ncradfn(infile,ilist)
          WRITE(6,'(a,i4,a,i4)') ' infile= ',infile,'  type2a=',fvartype(infile,ilist)
          CALL extdirn(infilname,indir,fname)
          CALL get_ncraddims(fname,indir,iradtype,ngate,nazim,nvar,        &
                           istatus)
          ALLOCATE(tem_double(nazim),stat=istat)
          CALL check_alloc_status(istat,'ncrad2arps:tem_double')
          ALLOCATE(ncvarname(nvar),stat=istat)
          CALL check_alloc_status(istat,'ncrad2arps:ncvarname')
          CALL get_ncrad2ainfo(fname,indir,iradtype,nazim,nvar,tem_double,    &
                   radarname,radlat,radlon,radalt,itimcdf,frtime,elv,         &
                   ncvarname,istatus)
          DEALLOCATE(tem_double,ncvarname)
          ALLOCATE(rhohv(ngate,nazim),stat=istat)
          CALL check_alloc_status(istat,'ncrad2arps:rhohv')
          ALLOCATE(zdr(ngate,nazim),stat=istat)
          CALL check_alloc_status(istat,'ncrad2arps:zdr')
          ALLOCATE(phidp(ngate,nazim),stat=istat)
          CALL check_alloc_status(istat,'ncrad2arps:phidp')
          IF(iradtype == 61) THEN
            ALLOCATE(snr(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:snr')
            ALLOCATE(tem_short2d(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:tem_short2d')
             print *, ' calling rdtiltforay...'
            CALL rdtiltforay(nazim,ngate,fname,indir,                  &
                    refvarname,velvarname,rhvvarname,                  &
                    zdrvarname,kdpvarname,snrvarname,                  &
                    dualpdata,rmisval,rngfval,itimcdf,frtime,          &
                    vnyquist,rfirstg,bmwidth,                          &
                    azim,gtspc,vnyq,refdata,veldata,rhvdata,           &
                    zdrdata,kdpdata,snr,tem_short2d)
             print *, ' back from rdtiltforay...'
            DEALLOCATE(snr,tem_short2d)
          ELSE
            ALLOCATE(attref(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:attref')
            ALLOCATE(refl(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:refl')
            ALLOCATE(radv(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:radv')
            ALLOCATE(gateqc(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:gateqc')
            ALLOCATE(tem_real2d(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:tem_real2d')
            ALLOCATE(tem_double(nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:tem_double')
            ngate1=ngate
            CALL rd2atiltcdf(nazim,ngate,iradtype,fname,indir,         &
                         refvarname,velvarname,rhvvarname,snrvarname,  &
                         zdrvarname,kdpvarname,                        &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyquist,rfirstg,bmwidth,istrgate,            &
                         gcfst,azim,gtspc,attref,refl,radv,rhohv,      &
                         zdr,phidp,gateqc,tem_double,tem_real2d)
            DEALLOCATE(tem_double)
            DEALLOCATE(tem_real2d)
            vnyq=vnyquist
            arfdata=refmiss
            refdata=refmiss
            veldata=velmiss
            rhvdata=rhvmiss
            zdrdata=zdrmiss
            kdpdata=kdpmiss
            kntgate=0
            badgate=0
            DO j=1,nazim
              istrgate(j) = max(istrgate(j),1)   ! to avoid 0 index
              IF(ngate /= ngate1) STOP
              IF (gcfst(j) > gcfchek) THEN
                DO i=istrgate(j),ngate
                  kntgate=kntgate+1
                  IF(btest(gateqc(i,j),0)) THEN
                    IF(abs(refl(i,j)) < 99. ) THEN
                      arfdata(i,j)=attref(i,j)
                      refdata(i,j)=refl(i,j)
                    END IF
                    IF(abs(radv(i,j)) < 199. ) veldata(i,j)=radv(i,j)
                    IF(abs(rhohv(i,j)) < 1.1 ) rhvdata(i,j)=rhohv(i,j)
                    IF(abs(  zdr(i,j)) < 20. ) zdrdata(i,j)=zdr(i,j)
                    IF(abs(phidp(i,j)) < 199.) kdpdata(i,j)=phidp(i,j)
                  ELSE
                    badgate=badgate+1
                  END IF
                END DO
              END IF
            END DO

            IF( kntgate > 0) THEN
              badpct=100.*(float(badgate)/float(kntgate))
            ELSE
              badpct=0.
            END IF
            WRITE(6,'(a,i7,a,i7,a,f8.2,a)')    &
              ' Gate check bad/missing flags ',badgate,' of ',kntgate,' gates, ', &
                badpct,' percent.'

            DEALLOCATE(attref, refl, radv, gateqc)
          END IF

          DEALLOCATE(rhohv)
          IF(ALLOCATED(zdr)) DEALLOCATE(zdr)
          IF(ALLOCATED(phidp)) DEALLOCATE(phidp)
!
!------------------------------------------------------------------------
!
!     Correct azimuth orientation and altitude
!
!------------------------------------------------------------------------
!
            radalt=radalt+altoffset
            IF( northazim /= 0.0) THEN
              DO j = 1, nazim
                azim(j) = azim(j)-northazim
                IF(azim(j) < 0.0) azim(j)=azim(j)+360.0
              END DO
            END IF
!
!------------------------------------------------------------------------
!
!     Update Nyquist Velocity with input values, if needed.
!
!------------------------------------------------------------------------
!
          IF( nfixed_nyqv > 0 ) THEN
            kfn=0
            DO kkfn=nfixed_nyqv,1,-1
              IF( elv <= elevlim_fnyqv(kkfn) ) kfn=kkfn
            END DO
            IF( kfn > 0 ) THEN
              WRITE(6,'(2(a,f6.2),a)') ' Resetting Nyquist at elev: ',  &
                          elv,' to ',fixed_nyqv(kfn),' m/s.'
              vnyquist=fixed_nyqv(kfn)
              vnyq(:) = vnyquist
            END IF
          END IF
!
!------------------------------------------------------------------------
!
!     Process reflectivity data.
!
!------------------------------------------------------------------------

          IF(tintrpopt == 0 .OR. ilist == 1) klev=klev+1
          elvang(klev) = elv

          WRITE(stdout,'(a,i4)') &
               ' Processing base reflectivity data... ', klev

          WRITE(stdout,*)'Transferring ',nazim,' reflectivity radials.'

          dtime=(itimcdf+itime1970)-itimfrst
          print *, '  1 dtime= ',dtime,' itimcdf= ',itimcdf,' itimfrst= ', &
                   itimfrst
          refelvmin=min(elv,refelvmin)
          refelvmax=max(elv,refelvmax)
          DO j = 1, maxazim
            time(j) = dtime
            elev(j) = elv
          END DO

          WRITE(stdout,'(a,i5,a,4f9.1)')                                 &
            '  Ref  j   azim   elev   refmin   refmax'
          DO j = 1, nazim, 10
            rmin=999.
            rmax=-999.
            DO i=1,ngate
              IF(refdata(i,j) > refchek) THEN
                rmin=min(refdata(i,j),rmin)
                rmax=max(refdata(i,j),rmax)
              END IF
            END DO
            WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
          END DO

          gatesp = gtspc(1)
          igatesp = NINT(gatesp)
          irfirstg = NINT(rfirstg)
          refrngmin=max(refrngmin,rfirstg)
          refrngmax=min(refrngmax,(rfirstg+((ngate-1)*gatesp)))

          IF( qc_on) THEN
            IF( radband == 3 .AND. snrthr > -91.) THEN
              write(6, *) ' calling SNR sensitivity screen'
              CALL snrchk(maxgate,maxazim,ngate,nazim,refchek,rfirstg,gatesp, &
                     snrthr,snstvty,refdata,veldata,arfdata)
            END IF

            write(6, *) ' calling RhoHV screen'
            CALL rhohvchk(maxgate,maxgate,maxazim,ngate,ngate,nazim,   &
                   igatesp,igatesp,refchek,velchek,rhohvthr,           &
                   refdata,veldata,rhvdata)

            write(6, *) ' calling despekl '
            CALL despekl(maxgate,maxazim,maxgate,maxazim,refchek,medfilt, &
                 refdata,rtem)
          END IF

          IF(gridtilt == 1) THEN
            DO j = 1, nazim
              IF(azim(j) < 0.) azim(j)=360.+azim(j)
            END DO
            CALL RadialIntp_rdr(ngate,numrad,nazim,azim,refdata,       &
                                  rdr_tmp(:,:,ilist))
          END IF

          rdr_p(:,1:numrad,klev,2) = rdr_tmp(:,:,ilist)

!-----------------------------------------------------------------------
!
!   Process polarimetric variables
!
!-----------------------------------------------------------------------
!
          IF(dualpproc .OR. gridtilt == 1) THEN
! Zdr
            WRITE(stdout,'(a,i5,a)')                              &
              '   Zdr   j   azim   elev  zdr min  zdr max'
            DO j = 1, nazim
              rmin=999.
              rmax=-999.
              DO i=1,ngate
                IF(zdrdata(i,j) > -20.) THEN
                  rmin=min(zdrdata(i,j),rmin)
                  rmax=max(zdrdata(i,j),rmax)
                END IF
              END DO
              IF(mod(j,5) == 0) &
                WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
            ENDDO

!            IF( qc_on) THEN
              DO j=1,nazim
                DO i=1,ngate
                  IF (rhvdata(i,j) < 0.8) THEN
                    zdrdata(i,j)= missing
                  END IF
                END DO
              END DO
!            END IF

            CALL RadialIntp_rdr(ngate,numrad,nazim,azim,zdrdata,  &
                                rdr_tmp(:,:,ilist))

            rdr_p(:,1:numrad,klev,4) = rdr_tmp(:,:,ilist)

! Phidp
            WRITE(stdout,'(a,i5,a)')                              &
              '   Kdp j   azim   elev  phi min  phi max'
            DO j = 1, nazim
              rmin=999.
              rmax=-999.
              DO i=1,ngate
                IF(kdpdata(i,j) .ge. -180. .and. kdpdata(i,j) .le. 180.) THEN
                  rmin=min(kdpdata(i,j),rmin)
                  rmax=max(kdpdata(i,j),rmax)
                END IF
              END DO
              IF(mod(j,5) == 0) &
                WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
            ENDDO

!            IF( qc_on) THEN
              DO j=1,nazim
                DO i=1,ngate
                  IF (rhvdata(i,j) < 0.8) THEN
                    kdpdata(i,j)= missing
                  END IF
                END DO
              END DO
!            END IF

            CALL RadialIntp_rdr(ngate,numrad,nazim,azim,kdpdata,  &
                                rdr_tmp(:,:,ilist))

            rdr_p(:,1:numrad,klev,6) = rdr_tmp(:,:,ilist)

! Rhohv
            WRITE(stdout,'(a,i5,a)')                              &
              '   Rhohv j   azim   elev  rhv min  rhv max'
            DO j = 1, nazim
              rmin=999.
              rmax=-999.
              DO i=1,ngate
                IF(rhvdata(i,j) > -0.1) THEN
                  rmin=min(rhvdata(i,j),rmin)
                  rmax=max(rhvdata(i,j),rmax)
                END IF
              END DO
              IF(mod(j,5) == 0) &
                WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
            ENDDO

            CALL RadialIntp_rdr(ngate,numrad,nazim,azim,rhvdata,  &
                                rdr_tmp(:,:,ilist))

            rdr_p(:,1:numrad,klev,7) = rdr_tmp(:,:,ilist)

            IF(gridtilt < 1) THEN
              write(6, *) ' Calling volbuild reflectivity & dualpol'
              write(6, *) '    elv=',elev(1),' nazim =',nazim
              timeset = 0
              nyqset = 0
              CALL volbuilddlp(maxgate,maxazim,maxelev,ngate,nazim,    &
                  nyqset,timeset,                                      &
                  kntrgat,kntrazm,kntrelv,                             &
                  igatesp,irfirstg,refchek,                            &
                  vnyq,time,                                           &
                  azim,elev,refdata,rhvdata,zdrdata,kdpdata,           &
                  vnyqvol,timevolr,                                    &
                  rngrvol,azmrvol,elvrvol,elvmnrvol,                   &
                  refvol,rhvvol,zdrvol,kdpvol)

            END IF

          ELSE ! no dual-pol processing

            write(6, *) ' Calling volbuild reflectivity'
            write(6, *) '    elv=',elev(1),' nazim =',nazim
            timeset = 0
            nyqset = 0
            CALL volbuild(maxgate,maxazim,maxelev,ngate,nazim,         &
                  nyqset,timeset,                                      &
                  kntrgat,kntrazm,kntrelv,                             &
                  igatesp,irfirstg,refchek,                            &
                  vnyq,time,                                           &
                  azim,elev,refdata,                                   &
                  vnyqvol,timevolr,                                    &
                  rngrvol,azmrvol,elvrvol,elvmnrvol,refvol)

          END IF !  dualpproc or gridtilt = 1

!------------------------------------------------------------------------
!
!       Process velocity data.
!
!------------------------------------------------------------------------

          WRITE(stdout,'(a,i4)') 'VCP number for this file:', ivcp
          WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO')    &
                 '    DATE: ', iyear, '/', imonth, '/', iday
          WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')                  &
                 '    TIME: ', ihr, ':', imin, ':', isec

          WRITE(stdout,'(a,i4)')' Processing base velocity data... ', klev

          WRITE(stdout,'(a,i6,a)')'Transferring',nazim,' velocity radials.'

          WRITE(stdout,'(a,f10.2)') ' Nyquist velocity: ',vnyquist

          dtime=(itimcdf+itime1970)-itimfrst
          print *, ' 2 dtime= ',dtime,' itimcdf= ',itimcdf,' itimfrst= ', &
                   itimfrst
          DO j = 1, maxazim
            time(j) = dtime
            elev(j) = elv
          END DO

          WRITE(stdout,'(a,i5,a,4f9.1)')                              &
            '  Vel  j   azim   elev   vr min   vr max'
          DO j = 1, nazim
            rmin=999.
            rmax=-999.
            DO i=1,ngate
              IF(veldata(i,j) > velchek) THEN
                rmin=min(veldata(i,j),rmin)
                rmax=max(veldata(i,j),rmax)
              END IF
            END DO
            IF(mod(j,5) == 0)                                          &
              WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
          END DO

          gatesp = gtspc(1)
          igatesp = NINT(gatesp)
          irfirstg = NINT(rfirstg)
          ivelgatsp=igatesp
          nyqset=1
          timeset=1

          IF( qc_on ) THEN
            IF( unfdiag ) THEN
              iangle=NINT(100.*elev(1))
              varid='rdvela'
              WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
              CALL wrtvel(iangle,itilt,varid,                          &
                          iiyear,iimon,iday,iihr,iimin,iisec,          &
                          igatesp,irfirstg,vnyquist,                   &
                          radname,radlat,radlon,radalt,                &
                          maxgate,maxazim,ngate,nazim,                 &
                          azim,elev,veldata)
            END IF

            print *, ' calling despekl '
            CALL despekl(maxgate,maxazim,maxgate,maxazim,velchek,medfilt,  &
                 veldata,rtem)

            IF( unfdiag ) THEN
              varid='rdvelb'
              WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
              CALL wrtvel(iangle,itilt,varid,                          &
                          iiyear,iimon,iday,iihr,iimin,iisec,          &
                          igatesp,irfirstg,vnyquist,                   &
                          radname,radlat,radlon,radalt,                &
                          maxgate,maxazim,ngate,nazim,                 &
                          azim,elev,veldata)
            END IF

            print *,' calling unfoldnqc 1'
            CALL unfoldnqc(maxgate,maxazim,maxgate,nazim,              &
                        nzsnd,zsnd,usnd,vsnd,rfrsnd,                   &
                        bkgopt,shropt,rfropt,                          &
                        igatesp,irfirstg,iangle,itilt,                 &
                        iiyear,iimon,iiday,iihr,iimin,iisec,           &
                        radlat,radlon,radalt,                          &
                        veldata,spwdata,elev,azim,vnyq,                &
                        unfvdata,bkgvel,bgate,egate,rtem)

            IF( unfdiag ) THEN
              varid='rdvelc'
              WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
              CALL wrtvel(iangle,itilt,varid,                          &
                          iiyear,iimon,iday,iihr,iimin,iisec,          &
                          igatesp,irfirstg,vnyquist,                   &
                          radname,radlat,radlon,radalt,                &
                          maxgate,maxazim,ngate,nazim,                 &
                          azim,elev,unfvdata)
            END IF

            IF(gridtilt == 1) THEN
              DO j = 1, nazim
                IF(azim(j) < 0.) azim(j)=360.+azim(j)
              ENDDO
              CALL RadialIntp_rdr(ngate,numrad,nazim,azim,unfvdata,    &
                                    rdr_tmp(:,:,ilist))
            ELSE 
              print *, ' Calling volbuild velocity 1'
              print *, '    elv=',elev(1),' nazim =',nazim
              print *, '    velvol(1,1,1) =',velvol(1,1,1)
              CALL volbuild(maxgate,maxazim,maxelev,maxgate,nazim,      &
                     nyqset,timeset,                                    &
                     kntvgat,kntvazm,kntvelv,                           &
                     igatesp,irfirstg,velchek,                          &
                     vnyq,time,                                         &
                     azim,elev,unfvdata,                                &
                     vnyqvol,timevolv,                                  &
                     rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)
             END IF

          ELSE

            IF(gridtilt == 1) THEN
              DO j = 1, nazim
                IF(azim(j) < 0.) azim(j)=360.+azim(j)
              END DO
              CALL RadialIntp_rdr(ngate,numrad,nazim,azim,veldata,    &
                                    rdr_tmp(:,:,ilist))

            ELSE            ! IF(gridtilt == 0) THEN

              print *, ' Calling volbuild velocity 2'
              print *, '    elv=',elev(1),' nazim =',nazim
              print *, '    velvol(1,1,1) =',velvol(1,1,1)
              CALL volbuild(maxgate,maxazim,maxelev,maxgate,nazim,      &
                     nyqset,timeset,                                    &
                     kntvgat,kntvazm,kntvelv,                           &
                     igatesp,irfirstg,velchek,                          &
                     vnyq,time,                                         &
                     azim,elev,veldata,                                 &
                     vnyqvol,timevolv,                                  &
                     rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)
            END IF

          END IF

          rdr_p(:,1:numrad,klev,1) = rdr_tmp(:,:,ilist)
!
        END DO  ! infile
        END DO  ! ilist

      ELSE IF(kntt2b > 0) THEN
!
!------------------------------------------------------------------------
!
!   Process reflectivity files
!
!------------------------------------------------------------------------
!
        DO infile = 1,maxinfile
        DO ilist = 1,nlistfil
          IF(fvartype(infile,ilist) == 2) THEN  ! type is reflectivity
            WRITE(6,'(a,i4,a,i4)') ' infile= ',infile,                 &
            ' type=',fvartype(infile,ilist)
            infilname = ncradfn(infile,ilist)
            WRITE(stdout,'(a,/a)')                                     &
              'Reading reflectivity file:',TRIM(infilname)
            CALL extdirn(infilname,indir,fname)
            CALL get_ncraddims(fname,indir,iradtype,ngate,nazim,nvar,  &
                             istatus)
            ALLOCATE(ncvarname(nvar),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:ncvarname')

            CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,      &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,    &
                     ivcp,elv,ncvarname,istatus)
            DEALLOCATE(ncvarname)

            ALLOCATE(refl(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:refl')
            ALLOCATE(gateqc(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:gateqc')
            CALL rdrftiltcdf(nazim,ngate,fname,indir,iradtype,refvarname, &
                     rmisval,rngfval,itimcdf,frtime,initime,rfirstg,      &
                     bmwidth,azim,gtspc,refl,gateqc)
            refdata=rmisval
            DO j=1,nazim
              DO i=1,ngate
                IF (btest(gateqc(i,j),0) .AND. refl(i,j) /= rmisval) THEN
                  refdata(i,j)=refl(i,j)
                END IF
              END DO
            END DO
            DEALLOCATE(refl, gateqc)

            IF(radarname(1:5) == 'cyril') THEN
              IF(radalt < 0.1) radalt=445.3
              IF(elv < 0.1) elv=1.5
            ELSE IF(radarname(1:5) == 'chick') THEN
              IF(radalt < 0.1) radalt=355.4
              IF(elv < 0.1) elv=1.5
            ELSE IF(radarname(1:5) == 'lawto') THEN
              IF(radalt < 0.1) radalt=295.0
              IF(elv < 0.1) elv=1.5
            ELSE IF(radarname(1:5) == 'rushs') THEN
              IF(radalt < 0.1) radalt=365.8
              IF(elv < 0.1) elv=1.5
            END IF

!------------------------------------------------------------------------
!
!       Now begin processing the data.
!
!------------------------------------------------------------------------

            IF(tintrpopt == 0 .OR. ilist == 1) THEN
              klev = klev + 1
            END IF
            elvang(klev) = elv

            WRITE(stdout,'(a,i4)') 'VCP number for this file:', ivcp
            WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO') &
                   '    DATE: ', iyear, '/', imonth, '/', iday
            WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')  &
                   '    TIME: ', ihr, ':', imin, ':', isec

!------------------------------------------------------------------------
!
!       Process reflectivity data.
!
!------------------------------------------------------------------------

            WRITE(stdout,'(a,i4)') &
                 ' Processing base reflectivity data... ', klev

            WRITE(stdout,*)'Transferring ',nazim,' reflectivity radials.'

            dtime=(itimcdf+itime1970)-itimfrst
            print *, ' 3 dtime= ',dtime,' itimcdf= ',itimcdf,' itimfrst= ', &
                   itimfrst
            refelvmin=min(elv,refelvmin)
            refelvmax=max(elv,refelvmax)
            DO j = 1, maxazim
              time(j) = dtime
              elev(j) = elv
            END DO

            WRITE(stdout,'(a,i5,a,4f9.1)')                                 &
              '  Ref  j   azim   elev   refmin   refmax'
            DO j = 20, nazim, 20
              rmin=999.
              rmax=-999.
              DO i=1,ngate
                IF(refdata(i,j) > refchek) THEN
                  rmin=min(refdata(i,j),rmin)
                  rmax=max(refdata(i,j),rmax)
                END IF
              END DO
              WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
            END DO

            gatesp  = gtspc(1)
            igatesp  = NINT(gatesp)
            rfirstg = NINT(rfirstg)
!------------------------------------------------------------------------
!
!         Call radar volume builder.
!
!------------------------------------------------------------------------
            irefgatsp=igatesp
            write(6, *) ' gatesp,rfirstg,refchek=',gatesp,rfirstg,refchek
            refrngmin=max(refrngmin,rfirstg)
            refrngmax=min(refrngmax,(rfirstg+((ngate-1)*gatesp)))
!
!------------------------------------------------------------------------
!
!       Get matching rho-HV data
!
!------------------------------------------------------------------------
!
            rho_match=.FALSE.
            DO jlist=1,nlistfil
            DO jfile=1,nncradfn(jlist)
              IF(fvartype(jfile,jlist) == 7 .AND.                      &
                 fitime(jfile,jlist) == itimcdf .AND.                  &
                 abs(felv(jfile,jlist)-elv) < 0.01 ) THEN
                rho_match=.TRUE.
                infilname = ncradfn(jfile,ilist)
                WRITE(stdout,'(a,a)') 'Reading Rho-HV file:',TRIM(infilname)
                CALL extdirn(infilname,indir,fname)

                ALLOCATE(rhohv(ngate,nazim),stat=istat)
                CALL check_alloc_status(istat,'ncrad2arps:rhohv')
                CALL rdrhotiltcdf(nazim,ngate,fname,indir,rhvvarname,  &
                           rmisval,rngfval,itimcdf,frtime,initime,     &
                           rfirstg,azim,rhohv)
                print *, ' back from rdrhotiltcdf '
                rhvdata=refmiss
                DO j=1,nazim
                  DO i=4,ngate
                    IF(rhohv(i,j) > -2.0 ) then
                      rhvdata(i,j)=rhohv(i,j)
                    END IF
                  END DO
                END DO
                DEALLOCATE(rhohv)
                EXIT
              END IF
            END DO  ! jfile loop
            END DO  ! jlist loop

            IF(radarname(1:5) == 'cyril') THEN
              IF(radalt < 0.1) radalt=445.3
              IF(elv < 0.1) elv=1.5
            ELSE IF(radarname(1:5) == 'chick') THEN
              IF(radalt < 0.1) radalt=355.4
              IF(elv < 0.1) elv=1.5
            ELSE IF(radarname(1:5) == 'lawto') THEN
              IF(radalt < 0.1) radalt=295.0
              IF(elv < 0.1) elv=1.5
            ELSE IF(radarname(1:5) == 'rushs') THEN
              IF(radalt < 0.1) radalt=365.8
              IF(elv < 0.1) elv=1.5
            END IF

            IF( qc_on) THEN
              IF(rho_match) THEN
                print *, ' calling RhoHV screen'
                CALL rhohvchk(maxgate,maxgate,maxazim,ngate,ngate,nazim,&
                              igatesp,igatesp,refchek,velchek,rhohvthr, &
                              refdata,veldata,rhvdata)
              END IF
              print *, ' calling despekl '
              CALL despekl(maxgate,maxazim,maxgate,maxazim,refchek,medfilt, &
                   refdata,rtem)
            END IF

            nyqset=0
            IF ( velproc ) THEN
              timeset=0
            ELSE
              timeset=1
            END IF
            print *, '    elv=',elev(1),' nazim =',nazim
            vnyq=0.

            IF(gridtilt == 1) THEN
              CALL RadialIntp_rdr(ngate,numrad,nazim,azim,refdata,     &
                                  rdr_tmp(:,:,ilist))

            ELSE            ! IF(gridtilt == 0) THEN
              print *, ' calling volbuild reflectivity'
              CALL volbuild(maxgate,maxazim,maxelev,ngate,nazim,       &
                      nyqset,timeset,                                  &
                      kntrgat,kntrazm,kntrelv,                         &
                      igatesp,irfirstg,refchek,                        &
                      vnyq,time,                                       &
                      azim,elev,refdata,                               &
                      vnyqvol,timevolr,                                &
                      rngrvol,azmrvol,elvrvol,elvmnrvol,refvol)
            END IF
          END IF ! Reflectivity
        END DO  ! list loop

        IF(gridtilt == 1) THEN
          IF(tintrpopt > 0 .AND. nlistfil > 1) THEN
            CALL TimeIntp_rdr(ngate,numrad,scantime(infile,:),ref_t,   &
                                rdr_tmp,rdr_p(:,:,klev,2))
          ELSE
            rdr_p(:,1:numrad,klev,2) = rdr_tmp(:,:,1)
          END IF
        END IF

        END DO  ! file in list loop

!------------------------------------------------------------------------
!
!   Process radial velocity files
!
!------------------------------------------------------------------------
!
        klev = 0
        DO infile = 1,maxinfile
        DO ilist = 1,nlistfil
          IF(fvartype(infile,ilist) == 1) THEN  ! type is velocity
            infilname = ncradfn(infile,ilist)
            WRITE(stdout,'(a,a)') 'Reading velocity file:',TRIM(infilname)
            CALL extdirn(infilname,indir,fname)
            CALL get_ncraddims(fname,indir,iradtype,ngate,nazim,nvar,  &
                             istatus)
            ALLOCATE(ncvarname(nvar), stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:ncvarname')
            CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,      &
                     radarname,radlat,radlon,radalt,itimcdf,frtime,          &
                     ivcp,elv,ncvarname,istatus)
            DEALLOCATE(ncvarname)

            veldata=velmiss
            ALLOCATE(radv(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:radv')
            ALLOCATE(gateqc(ngate,nazim),stat=istat)
            CALL check_alloc_status(istat,'ncrad2arps:gateqc')
            CALL rdvrtiltcdf(nazim,ngate,fname,indir,iradtype,velvarname, &
                           rmisval,rngfval,itimcdf,frtime,initime,        &
                           vnyquist,rfirstg,                              &
                           bmwidth,azim,gtspc,vnyq,radv,gateqc)
            DO j=1,nazim
              DO i=1,ngate
                IF (btest(gateqc(i,j),0) .AND. radv(i,j)/=rmisval) then
                  veldata(i,j)=radv(i,j)
                END IF
              END DO
            END DO
            DEALLOCATE(radv, gateqc)

!------------------------------------------------------------------------
!
!       Correction for early experimental CASA data
!
!------------------------------------------------------------------------

            IF(radarname(1:5) == 'cyril' .OR.  &
               radarname(1:5) == 'chick' .OR.  &
               radarname(1:5) == 'lawto' .OR.  &
               radarname(1:5) == 'rushs' ) THEN
              vconst=15./acos(-1.)
              DO j=1,nazim
                DO i=1,ngate
                  IF(veldata(i,j) > -4.0) THEN
                    veldata(i,j)=veldata(i,j)*vconst
                  END IF
                END DO
              END DO
              vnyquist=15.0
              vnyq(:)=vnyquist
            END IF
!------------------------------------------------------------------------
!
! Update Nyquist Velocity with input values, if needed.
!
!------------------------------------------------------------------------

            IF( nfixed_nyqv > 0 ) THEN
              kfn=0
              DO kkfn=nfixed_nyqv,1,-1
                IF( elv <= elevlim_fnyqv(kkfn) ) kfn=kkfn
              END DO
              IF( kfn > 0 ) THEN
                WRITE(6,'(2(a,f6.2),a)') ' Resetting Nyquist at elev: ',  &
                            elv,' to ',fixed_nyqv(kfn),' m/s.'
                vnyquist=fixed_nyqv(kfn)
                vnyq(:) = vnyquist
              END IF
            END IF

!------------------------------------------------------------------------
!
!         Process velocity data.
!
!------------------------------------------------------------------------

              IF(tintrpopt == 0 .OR. ilist == 1) klev = klev + 1

              WRITE(stdout,'(a,i4)') 'VCP number for this file:', ivcp
              WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO')   &
                     '    DATE: ', iyear, '/', imonth, '/', iday
              WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')                 &
                     '    TIME: ', ihr, ':', imin, ':', isec


              WRITE(stdout,'(a,i4)')' Processing base velocity data... ', klev

              WRITE(stdout,'(a,i6,a)')'Transferring',nazim,' velocity radials.'

              WRITE(stdout,'(a,f10.2)') ' Nyquist velocity: ',vnyquist

            dtime=(itimcdf+itime1970)-itimfrst
            print *, ' 4 dtime= ',dtime,' itimcdf= ',itimcdf,' itimfrst= ', &
                   itimfrst
            DO j = 1, maxazim
              time(j) = dtime
              elev(j) = elv
            END DO

            WRITE(stdout,'(a,i5,a,4f9.1)')                                 &
              '  Vel  j   azim   elev   velmin   velmax'
            DO j = 20, nazim, 20
              rmin=999.
              rmax=-999.
              DO i=1,ngate
                IF(veldata(i,j) > velchek) THEN
                  rmin=min(veldata(i,j),rmin)
                  rmax=max(veldata(i,j),rmax)
                END IF
              END DO
              WRITE(stdout,'(2x,i5,4f9.1)') j,azim(j),elev(j),rmin,rmax
            END DO

            gatesp = gtspc(1)
            igatesp = NINT(gtspc(1))
            irfirstg = NINT(rfirstg)
!------------------------------------------------------------------------
!
!         Call radar volume builder.
!
!------------------------------------------------------------------------
            ivelgatsp=igatesp
            nyqset=1
            timeset=1

              IF( qc_on ) THEN
                IF( unfdiag ) THEN
                  iangle=NINT(100.*elev(1))
                  varid='rdvela'
                  WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
                  CALL wrtvel(iangle,itilt,varid,                         &
                            iiyear,iimon,iday,iihr,iimin,iisec,           &
                            igatesp,irfirstg,vnyquist,                    &
                            radname,radlat,radlon,radalt,                 &
                            maxgate,maxazim,ngate,nazim,                  &
                            azim,elev,veldata)
                END IF

              print *, ' calling despekl '
              CALL despekl(maxgate,maxazim,maxgate,maxazim,velchek,medfilt,  &
                   veldata,rtem)

              IF( unfdiag ) THEN
                varid='rdvelb'
                WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
                CALL wrtvel(iangle,itilt,varid,                        &
                          iiyear,iimon,iday,iihr,iimin,iisec,          &
                          igatesp,irfirstg,vnyquist,                   &
                          radname,radlat,radlon,radalt,                &
                          maxgate,maxazim,ngate,nazim,                 &
                          azim,elev,veldata)
              END IF

              print *,' calling unfoldnqc 2'
              CALL unfoldnqc(maxgate,maxazim,maxgate,nazim,            &
                        nzsnd,zsnd,usnd,vsnd,rfrsnd,                   &
                        bkgopt,shropt,rfropt,                          &
                        igatesp,irfirstg,iangle,itilt,                 &
                        iiyear,iimon,iiday,iihr,iimin,iisec,           &
                        radlat,radlon,radalt,                          &
                        veldata,spwdata,elev,azim,vnyq,                &
                        unfvdata,bkgvel,bgate,egate,rtem)

              IF( unfdiag ) THEN
                varid='rdvelc'
                WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
                CALL wrtvel(iangle,itilt,varid,                        &
                          iiyear,iimon,iday,iihr,iimin,iisec,          &
                          igatesp,irfirstg,vnyquist,                   &
                          radname,radlat,radlon,radalt,                &
                          maxgate,maxazim,ngate,nazim,                 &
                          azim,elev,unfvdata)
              END IF

              IF(gridtilt == 1) THEN
                print *, ' Calling RadialIntrp_rdr velocity 3'
                print *, '    elv=',elev(1),' nazim =',nazim
                CALL RadialIntp_rdr(ngate,numrad,nazim,azim,unfvdata,  &
                                  rdr_tmp(:,:,ilist))
              ELSE 
                print *, ' Calling volbuild velocity 3'
                print *, '    elv=',elev(1),' nazim =',nazim
                CALL volbuild(maxgate,maxazim,maxelev,maxgate,nazim,   &
                       nyqset,timeset,                                 &
                       kntvgat,kntvazm,kntvelv,                        &
                       igatesp,irfirstg,velchek,                       &
                       vnyq,time,                                      &
                       azim,elev,unfvdata,                             &
                       vnyqvol,timevolv,                               &
                       rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)
              END IF

            ELSE    ! no qc

              IF(gridtilt == 1) THEN
                print *, ' Calling RadialIntrp_rdr velocity 4'
                print *, '    elv=',elev(1),' nazim =',nazim
                CALL RadialIntp_rdr(ngate,numrad,nazim,azim,veldata,   &
                                  rdr_tmp(:,:,ilist))
              ELSE 
                print *, ' Calling volbuild velocity 4'
                print *, '    elv=',elev(1),' nazim =',nazim
                CALL volbuild(maxgate,maxazim,maxelev,maxgate,nazim,   &
                       nyqset,timeset,                                 &
                       kntvgat,kntvazm,kntvelv,                        &
                       igatesp,irfirstg,velchek,                       &
                       vnyq,time,                                      &
                       azim,elev,veldata,                              &
                       vnyqvol,timevolv,                               &
                       rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)
              END IF
            END IF  ! qc ?
            
          END IF ! velocity

          END DO ! list loop

          IF(gridtilt == 1) THEN
            IF(tintrpopt > 0 .AND. nlistfil > 1) then
              CALL TimeIntp_rdr(ngate,numrad,scantime(infile,:),ref_t,  &
                                rdr_tmp,rdr_p(:,:,klev,1))
            ELSE
              rdr_p(:,1:numrad,klev,1) = rdr_tmp(:,:,1)
            END IF
          END IF

        END DO ! DO ifile = 1,maxinfile
!
!------------------------------------------------------------------------
!
!   Process polarimetric variable files
!
!------------------------------------------------------------------------
!
        IF(gridtilt == 1) THEN

        DO idp = 4,numvar
          klev = 0
          refdata=missing

          IF(idp == 5) CYCLE

          DO infile = 1,maxinfile
            DO ilist = 1,nlistfil
              print*, infile, ilist, fvartype(infile,ilist), idp
              IF(fvartype(infile,ilist) == idp) THEN  ! type is dualvarname(idp)
                infilname = ncradfn(infile,ilist)
                WRITE(stdout,'(a,a,a,/a)')                                         &
                  'Reading ',dualvarname(idp),' file:',TRIM(infilname)
                CALL extdirn(infilname,indir,fname)
                CALL get_ncraddims(fname,indir,iradtype,ngate,nazim,nvar,      &
                                 istatus)
                ALLOCATE(ncvarname(nvar),stat=istat)
                CALL check_alloc_status(istat,'ncrad2arps:ncvarname')
                CALL get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,          &
                         radarname,radlat,radlon,radalt,itimcdf,frtime,        &
                         ivcp,elv,ncvarname,istatus)
                DEALLOCATE(ncvarname)

                ALLOCATE(refl(ngate,nazim),stat=istat)    ! recycle refl array
                CALL check_alloc_status(istat,'ncrad2arps:refl')
                ALLOCATE(gateqc(ngate,nazim),stat=istat)
                CALL check_alloc_status(istat,'ncrad2arps:gateqc')
                CALL rdpftiltcdf(nazim,ngate,fname,indir,idp,                  &
                               zvvvarname,zdrvarname,zdpvarname,               &
                               rhvvarname,kdpvarname,                          &
                               rmisval,rngfval,itimcdf,frtime,initime,rfirstg, &
                               bmwidth,azim,gtspc,refl,gateqc)
                refdata=rmisval

                SELECT CASE (idp)
                CASE(4:6)
                  dpcheck = 20.
                CASE(7)
                  dpcheck = 1.
                END SELECT

                DO j=1,nazim
                  DO i=1,ngate
                    IF (btest(gateqc(i,j),0) .AND. refl(i,j) /= rmisval .and.  &
                        ABS(refl(i,j)) <= dpcheck) THEN
                      refdata(i,j)=refl(i,j)
                    END IF
                  END DO
                END DO
                DEALLOCATE(refl)
                DEALLOCATE(gateqc)

!------------------------------------------------------------------------
!
!         Now begin processing the data.
!
!------------------------------------------------------------------------

                IF(tintrpopt == 0 .OR. ilist == 1) klev = klev + 1

                WRITE(stdout,'(a,i4)') 'VCP number for this file:', ivcp
                WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO') &
                       '    DATE: ', iyear, '/', imonth, '/', iday
                WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')            &
                       '    TIME: ', ihr, ':', imin, ':', isec

!------------------------------------------------------------------------
!
!         Process polarimetric data.
!
!------------------------------------------------------------------------

                WRITE(stdout,'(a,i4)')                                 &
                     ' Processing base polarimetric data... ', klev

                WRITE(stdout,*)'Transferring ',nazim,' dualpol radials.'

                dtime=(itimcdf+itime1970)-itimfrst
                print *, ' 5 dtime= ',dtime,' itimcdf= ',itimcdf,' itimfrst= ', &
                   itimfrst
                refelvmin=min(elv,refelvmin)
                refelvmax=max(elv,refelvmax)
                DO j = 1, maxazim
                  time(j) = dtime
                  elev(j) = elv
                END DO

                WRITE(stdout,'(a,i5,a,4f9.1)')                         &
                  '  var  j   azim   elev   varmin   varmax'
                DO j = 20, nazim, 20
                  rmin=999.
                  rmax=-999.
                  DO i=1,ngate
                    IF(refdata(i,j) > refchek) THEN
                      rmin=min(refdata(i,j),rmin)
                      rmax=max(refdata(i,j),rmax)
                    END IF
                  END DO
                  WRITE(stdout,'(2x,i5,4f12.4)') j,azim(j),elev(j),rmin,rmax
                END DO

                gatesp  = gtspc(1)
                igatesp  = NINT(gatesp)
                rfirstg = NINT(rfirstg)

!------------------------------------------------------------------------
!
!           Call radar volume builder.
!
!------------------------------------------------------------------------
                rho_match=.FALSE.
                DO jlist=1,nlistfil
                DO jfile=1,nncradfn(jlist)
                  IF(fvartype(jfile,jlist) == 7 .AND.                  &
                     fitime(jfile,jlist) == itimcdf .AND.              &
                     abs(felv(jfile,jlist)-elv) < 0.01 ) THEN
                    rho_match=.TRUE.
                    infilname = ncradfn(jfile,jlist)
                    WRITE(stdout,'(a,a)') 'Reading Rho-HV file:',TRIM(infilname)
                    CALL extdirn(infilname,indir,fname)

                    ALLOCATE(rhohv(ngate,nazim),stat=istat)
                    CALL check_alloc_status(istat,'ncrad2arps:rhohv')
                    CALL rdrhotiltcdf(nazim,ngate,fname,indir,rhvvarname,  &
                               rmisval,rngfval,itimcdf,frtime,initime,     &
                               rfirstg,azim,rhohv)
                    print *, ' back from rdrhotiltcdf '
                    rhvdata=refmiss
                    DO j=1,nazim
                      DO i=4,ngate
                        IF(rhohv(i,j) > -2.0) then
                          rhvdata(i,j)=rhohv(i,j)
                        END IF
                      END DO
                    END DO
                    DEALLOCATE(rhohv)
                    EXIT
                  END IF
                END DO  ! jfile loop
                END DO  ! jlist loop

                irefgatsp=igatesp
                write(6, *) ' gatesp,rfirstg,varchek=',gatesp,rfirstg,refchek

                IF( qc_on) THEN
                  IF(idp /= 7) THEN          ! ONLY WHEN NOT RhoHV
                    DO j=1,nazim
                      DO i=1,ngate
                        IF (rhvdata(i,j) < rhohv_thr) THEN
                          refdata(i,j)= missing
                        END IF
                      END DO
                    END DO
                  ENDIF
                END IF

                print *, ' calling volbuild dualpol'
                print *, '    elv=',elev(1),' nazim =',nazim

                IF(gridtilt == 1) THEN
                  CALL RadialIntp_rdr(ngate,numrad,nazim,azim,refdata,        &
                                      rdr_tmp(:,:,ilist))
                ENDIF
              END IF  ! type is dualpol

            END DO  ! list loop

            IF(gridtilt == 1) THEN
              IF(tintrpopt > 0 .AND. nlistfil > 1) then
                CALL TimeIntp_rdr(ngate,numrad,scantime(infile,:),ref_t,    &
                                    rdr_tmp,rdr_p(:,:,klev,idp))
              ELSE
                rdr_p(:,1:numrad,klev,idp) = rdr_tmp(:,:,1)
              END IF
            END IF

          END DO  ! ifile loop
        END DO  ! idp loop
        END IF !  gridtilt = 1

      END IF  ! kntt2b > 0
!
!   Restore any gzipped files
!
      DO ilist=1,nlistfil
        DO infile=1,nncradfn(ilist)
          IF(gzipped(infile,ilist)) THEN
            WRITE(syscall,'(a,1x,a)') 'gzip',TRIM(ncradfn(infile,ilist))
            CALL system(syscall)
          END IF
        END DO
      END DO

    END IF ! Only ROOT READS FILES

  END IF ! IF (maxinfile > 0)

  runname=runname_in
  CALL gtlfnkey(runname, lfnkey)
  curtim=outtime

  IF(gridtilt /= 1) THEN

!------------------------------------------------------------------------
!
! Process reflectivity to remove AP and ground clutter.
!
!------------------------------------------------------------------------

  CALL mpupdatei(klev,1)

  IF (klev > 0) THEN

    IF (myproc == ROOT) THEN
      IF( qc_on ) THEN

        write(6, *) ' Calling apdetect ',irefgatsp,ivelgatsp
        CALL apdetect(maxgate,maxgate,maxazim,maxelev,                     &
                      kntrgat,kntrazm,kntrelv,                             &
                      kntvgat,kntvazm,kntvelv,                             &
                      refchek,velchek,                                     &
                      irefgatsp,ivelgatsp,                                 &
                      winszrad,winszazim,ivcp,gcopt,gcvrlim,               &
                      rngrvol,azmrvol,elvrvol,                             &
                      rngvvol,azmvvol,elvvvol,                             &
                      refvol,velvol,rtem,                                  &
                      istatus)
      END IF

    END IF

    CALL mpupdater(bmwidth,1)

!------------------------------------------------------------------------
!
!   Call remapping routines
!
!------------------------------------------------------------------------

    nyqset=0
    timeset=1
    IF( velproc ) timeset=0
    vardump=0
    IF(ref3d) vardump=1
    varfill=fillref
    ivar=1

    IF (myproc == ROOT) print *, ' Calling remapvol for reflectivity '

    CALL mpupdatei(kntrgat,maxazim*maxelev)
    CALL mpupdatei(kntrazm,maxelev)
    CALL mpupdatei(kntrelv,1)

    CALL mpupdater(rngrvol,maxgate*maxelev)
    CALL mpupdater(azmrvol,maxazim*maxelev)
    CALL mpupdater(elvrvol,maxazim*maxelev)
    CALL mpupdater(refvol,maxgate*maxazim*maxelev)
    CALL mpupdatei(timevolr,maxazim*maxelev)
    CALL mpupdatei(timevolv,maxazim*maxelev)
    CALL mpupdater(vnyqvol,maxazim*maxelev)

    IF(velproc) THEN
      CALL mpupdatei(kntvgat,maxazim*maxelev)
      CALL mpupdatei(kntvazm,maxelev)
      CALL mpupdatei(kntvelv,1)

      CALL mpupdater(rngvvol,maxgate*maxelev)
      CALL mpupdater(azmvvol,maxazim*maxelev)
      CALL mpupdater(elvvvol,maxazim*maxelev)
      CALL mpupdater(velvol,maxgate*maxazim*maxelev)
    END IF
!
!   VAD Wind Processing
!
    IF( (myproc == ROOT) .AND. vad ) THEN

      ALLOCATE(vadhgt(maxelev),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:vadhgt')
      ALLOCATE(vaddir(maxelev),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:vaddir')
      ALLOCATE(vadspd(maxelev),stat=istat)
      CALL check_alloc_status(istat,'ncrad2arps:vadspd')

      CALL mkvadfnm(dirname,radname,                                 &
                    iiyear,iimon,iiday,iihr,iimin,iisec,             &
                    vad_fname)
      WRITE(stdout,'(a,a)')                                          &
            ' Filename for VAD output: ',TRIM(vad_fname)
      CALL flush(stdout)

      CALL vadvol(maxgate,maxazim,maxelev,                           &
                  radalt,velchek,vadradius,vadwidth,vadminknt,       &
                  kntvgat,kntvazm,kntvelv,                           &
                  rngvvol,azmvvol,elvvvol,velvol,                    &
                  vadhgt,vaddir,vadspd)
      CALL wtvadprf(maxelev,vad_fname,radname,vadsrc,                &
                    radlat,radlon,radalt,vadhgt,vaddir,vadspd);

      DEALLOCATE(vadhgt,vaddir,vadspd)

    END IF
!
!   Set-up arrays needed for remapping
!
    IF (myproc == ROOT) WRITE(stdout,'(a)') ' Calling rmpsetup...'
    CALL rmpsetup(maxgate,maxgate,maxgate,maxazim,maxelev,             &
                    nx,ny,nxny,nz,nzsnd,                               &
                    rfropt,refchek,velchek,bmwidth,velproc,            &
                    kntrgat,kntrazm,kntrelv,                           &
                    kntvgat,kntvazm,kntvelv,                           &
                    radlat,radlon,radarx,radary,radalt,                &
                    dazim,rngmin,rngmax,                               &
                    rngrvol,azmrvol,elvrvol,                           &
                    rngvvol,azmvvol,elvvvol,                           &
                    refvol,velvol,rxvol,ryvol,rzvol,                   &
                    xs,ys,zps,zsnd,rfrsnd,ncolp,ncoltot,               &
                    havdat,icolp,jcolp,xcolp,ycolp,zcolp,istatus)
    IF (myproc == ROOT) &
      WRITE(stdout,'(a,i12)') ' Back from rmpsetup, ncoltot=',ncoltot
    WRITE(stdout,'(a,i6,a,i12)')' myproc: ',myproc,'  ncolp: ',ncolp

    DEALLOCATE(havdat)
!
!   For MPI run, distribute the non-missing columns equally
!   among all processors.
!
    IF (nprocs > 1 ) THEN

      CALL distrallcol(nxny,nz,ncolp,ncoltot,                           &
                       icolp,jcolp,xcolp,ycolp,zcolp,lvldbg,istatus)

    END IF
!
!  Call remap for reflectivity.
!  Note here that the rxvol,ryvol and rzvol arrays have previously been
!  set for the reflectivity gates.
!
    IF (myproc == ROOT) print *, ' Calling remap3dcol for reflectivity '

    CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,           &
                  nzsnd,nsort,ncolp,                                   &
                  varfill,ivar,nyqset,timeset,rfropt,                  &
                  refchek,refmiss,bmwidth,refmedl,refdazl,iordref,     &
                  sortmin,dsort,                                       &
                  kntrgat,kntrazm,kntrelv,                             &
                  radlat,radlon,radarx,radary,radalt,dazim,            &
                  rngmin,rngmax,                                       &
                  rngrvol,azmrvol,elvrvol,                             &
                  refvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,           &
                  zsnd,rfrsnd,kntbin,elvmnrvol,                        &
                  xcolp,ycolp,zcolp,                                   &
                  colref,coltim,colnyq,istatus)

    DO kelv=1,kntrelv
      refelvmin=min(elvmnrvol(kelv),refelvmin)
      refelvmax=max(elvmnrvol(kelv),refelvmax)
    END DO

    IF(myproc == ROOT) print *, ' outtime: ',outtime

    IF( ref3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,                   &
                  colref,icolp,jcolp,                                  &
                  refid,refname,refunits,refmiss,outtime,runname,      &
                  dirname,dmpfmt,hdfcompr,istatus)

    IF( ref2d ) THEN
      vardump = 1
      ivar    = 1
      varid   = 'refl2d'
      varname = 'Low-level reflect'
      IF (myproc == ROOT) print *, ' Calling remap2d for reflectivity'
      CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,          &
                   vardump,ivar,rfropt,varid,varname,refunits,         &
                   dmpfmt,hdf4cmpr,                                    &
                   refchek,refmiss,refmedl,refdazl,iordref,            &
                   sortmin,dsort,                                      &
                   kntrgat,kntrazm,kntrelv,                            &
                   radlat,radlon,radarx,radary,radalt,dazim,           &
                   rngmin,rngmax,rngrvol,azmrvol,elvrvol,              &
                   refvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,        &
                   tem2d,istatus)

    END IF

    IF( dualpproc ) THEN
      nyqset=0
      timeset=0
      varfill=0
      IF( rhv3d .OR. dualpout > 0 ) THEN
        ivar=4
        IF (myproc == ROOT) print *, ' Calling remap3dcol for Rho-HV'
        CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,       &
                  nzsnd,nsort,ncolp,                                   &
                  varfill,ivar,nyqset,timeset,rfropt,                  &
                  rhvchek,rhvmiss,bmwidth,rhvmedl,refdazl,iordref,     &
                  sortmin,dsort,                                       &
                  kntrgat,kntrazm,kntrelv,                             &
                  radlat,radlon,radarx,radary,radalt,dazim,            &
                  rngmin,rngmax,                                       &
                  rngrvol,azmrvol,elvrvol,                             &
                  rhvvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,           &
                  zsnd,rfrsnd,kntbin,elvmnrvol,                        &
                  xcolp,ycolp,zcolp,                                   &
                  colrhv,coltim,colnyq,istatus)

        IF( rhv3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,               &
                      colrhv,icolp,jcolp,                              &
                      rhvid,rhvname,rhvunits,rhvmiss,outtime,runname,  &
                      dirname,dmpfmt,hdfcompr,istatus)
      END IF 

      IF( zdr3d .OR. dualpout > 0 ) THEN
        IF (myproc == ROOT) print *, ' Calling remap3dcol for Zdr'
        ivar=5
        CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,       &
                  nzsnd,nsort,ncolp,                                   &
                  varfill,ivar,nyqset,timeset,rfropt,                  &
                  zdrchek,zdrmiss,bmwidth,zdrmedl,refdazl,iordref,     &
                  sortmin,dsort,                                       &
                  kntrgat,kntrazm,kntrelv,                             &
                  radlat,radlon,radarx,radary,radalt,dazim,            &
                  rngmin,rngmax,                                       &
                  rngrvol,azmrvol,elvrvol,                             &
                  zdrvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,           &
                  zsnd,rfrsnd,kntbin,elvmnrvol,                        &
                  xcolp,ycolp,zcolp,                                   &
                  colzdr,coltim,colnyq,istatus)

        IF( zdr3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,               &
                      colzdr,icolp,jcolp,                              &
                      zdrid,zdrname,zdrunits,zdrmiss,outtime,runname,  &
                      dirname,dmpfmt,hdfcompr,istatus)
      END IF  ! zdr3d

      IF( kdp3d .OR. dualpout > 0 ) THEN
        IF (myproc == ROOT) print *, ' Calling remap3dcol for Kdp'
        ivar=6
        CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,       &
                  nzsnd,nsort,ncolp,                                   &
                  varfill,ivar,nyqset,timeset,rfropt,                  &
                  kdpchek,kdpmiss,bmwidth,kdpmedl,refdazl,iordref,     &
                  sortmin,dsort,                                       &
                  kntrgat,kntrazm,kntrelv,                             &
                  radlat,radlon,radarx,radary,radalt,dazim,            &
                  rngmin,rngmax,                                       &
                  rngrvol,azmrvol,elvrvol,                             &
                  kdpvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,           &
                  zsnd,rfrsnd,kntbin,elvmnrvol,                        &
                  xcolp,ycolp,zcolp,                                   &
                  colkdp,coltim,colnyq,istatus)

        IF( kdp3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,               &
                      colkdp,icolp,jcolp,                              &
                      kdpid,kdpname,kdpunits,kdpmiss,outtime,runname,  &
                      dirname,dmpfmt,hdfcompr,istatus)
      END IF  ! kdp3d

      IF( rhv2d ) THEN
        varid   = 'RHV_2d'
        varname = 'Low-level Rho-HV'
        vardump = 1
        ivar    = 4
        IF (myproc == ROOT) print *, ' Calling remap2d for Rho-HV'
        CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,        &
                     vardump,ivar,rfropt,varid,rhvname,rhvunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     rhvchek,rhvmiss,rhvmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     refvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)

      END IF  ! rhv2d

      IF( zdr2d ) THEN
        varid   = 'Zdr_2d'
        varname = 'Low-level Zdr'
        vardump = 1
        ivar    = 5
        IF (myproc == ROOT) print *, ' Calling remap2d for Zdr'
        CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,        &
                     vardump,ivar,rfropt,varid,zdrname,zdrunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     zdrchek,zdrmiss,zdrmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     zdrvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)

      END IF  ! zdr2d

      IF( kdp2d ) THEN
        varid   = 'Kdp_2d'
        varname = 'Low-level Kdp'
        vardump = 1
        ivar    = 6
        IF (myproc == ROOT) print *, ' Calling remap2d for Kdp'
        CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,        &
                     vardump,ivar,rfropt,varid,kdpname,kdpunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     kdpchek,kdpmiss,zdrmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     refvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)

      END IF  ! kdp2d

    ELSE
      WRITE(6,'(1x,a)') 'No remapping of dual-pol variables.'
    END IF

    IF (myproc == ROOT) WRITE(6,'(1x,a,l9)')                           &
       'Beginning velocity processing, velproc =',velproc
    IF( velproc ) THEN
      nyqset  = 1
      timeset = 1
      vardump = 0
      IF(vel3d) vardump = 1
      varfill = 0
      ivar    = 2

      CALL rgatexyz(maxgate,maxgate,maxazim,maxelev,nzsnd,             &
                    rfropt,lvldbg,                                     &
                    kntvgat,kntvazm,kntvelv,                           &
                    radlat,radlon,radarx,radary,radalt,                &
                    rngvvol,azmvvol,elvvvol,                           &
                    zsnd,rfrsnd,                                       &
                    rxvol,ryvol,rzvol,istatus)

      IF (myproc == ROOT) WRITE(stdout,'(1x,a)')                       &
         'Calling remap3dcol for velocity '
      CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,         &
                  nzsnd,nsort,ncolp,                                   &
                  varfill,ivar,nyqset,timeset,rfropt,                  &
                  velchek,velmiss,bmwidth,velmedl,veldazl,iordvel,     &
                  sortmin,dsort,                                       &
                  kntvgat,kntvazm,kntvelv,                             &
                  radlat,radlon,radarx,radary,radalt,dazim,            &
                  rngmin,rngmax,                                       &
                  rngvvol,azmvvol,elvvvol,                             &
                  velvol,timevolv,vnyqvol,rxvol,ryvol,rzvol,           &
                  zsnd,rfrsnd,kntbin,elvmnvvol,                        &
                  xcolp,ycolp,zcolp,                                   &
                  colvel,coltim,colnyq,istatus)

      IF( vel3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,                 &
                    colvel,icolp,jcolp,                                &
                    velid,velname,velunits,velmiss,outtime,runname,    &
                    dirname,dmpfmt,hdfcompr,istatus)
    END IF

    IF( vel2d ) THEN
      vardump  = 1
      ivar     = 1
      varid    = 'radv2d'
      varname  = 'Low-level Velocity'
      IF (myproc == ROOT) WRITE(stdout,'(1x,a)')                       &
         'Calling remap2d for velocity '
      CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,          &
                   vardump,ivar,rfropt,varid,varname,velunits,         &
                   dmpfmt,hdf4cmpr,                                    &
                   velchek,velmiss,velmedl,veldazl,iordvel,            &
                   sortmin,dsort,                                      &
                   kntvgat,kntvazm,kntvelv,                            &
                   radlat,radlon,radarx,radary,radalt,dazim,           &
                   rngmin,rngmax,rngvvol,azmvvol,elvvvol,              &
                   velvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,        &
                   tem2d,istatus)
    END IF

!------------------------------------------------------------------------
!
!   Create filename and write remapped file.
!
!------------------------------------------------------------------------

    full_fname = ' '
    CALL mkradfnm(dmpfmt,dirname,ldirnam,radname,iiyear,iimon,iiday,   &
                  iihr, iimin, iisec, full_fname, len_fname)

    IF (myproc == ROOT) WRITE(stdout,'(a,a)') &
        ' Filename for this volume: ',TRIM(full_fname)
    iradfmt=1
    print *, ' ncolp: ',ncolp,' ncoltot: ',ncoltot

    tmin=1.0E09
    tmax=-1.0E09
    DO icol=1,ncolp
      DO k=1,nz
        IF( coltim(k,icol) > -900.) THEN
          tmin=min(tmin,coltim(k,icol))
          tmax=max(tmax,coltim(k,icol))
        END IF
      END DO
    END DO
    print *, ' Calling wrtradcol, time_min =',tmin,' time_max =',tmax

    CALL wrtradcol(nx,ny,nxny,nz,ncolp,ncoltot,                        &
                  dmpfmt,iradfmt,hdf4cmpr,dmpzero,dualpout,            &
                  full_fname,radname,radlat,radlon,radalt,             &
                  iiyear,iimon,iiday,iihr,iimin,iisec,ivcp,isource,    &
                  refelvmin,refelvmax,refrngmin,refrngmax,             &
                  icolp,jcolp,xcolp,ycolp,zcolp,                       &
                  colref,colvel,colnyq,coltim,colrhv,colzdr,colkdp,    &
                  istatus)

  END IF

!------------------------------------------------------------------------
!
! The End.
!
!------------------------------------------------------------------------

  IF (myproc == ROOT)  &
    WRITE(stdout,'(/a/)') '  === Normal termination of NCRAD2ARPS === '

    CALL mpexit(0)

    STOP

  ELSE  ! gridtilt

    ALLOCATE(azmvol(numrad2,maxelev))
    ALLOCATE(elvvol(numrad2,maxelev))

    ALLOCATE(gridtilt_height(nx,ny,maxelev))
    ALLOCATE(gridtilt_range(nx,ny,maxelev))
    ALLOCATE(gridtilt_slr(nx,ny))
    ALLOCATE(gridtilt_azm(nx,ny))
    ALLOCATE(gridtilt_rdr(nx,ny,maxelev,numvar))
    ALLOCATE(gridtilt_time(maxelev))
    ALLOCATE(gridtilt_rdr_mean(maxelev))

    gridtilt_rdr = missing

    DO k=1,maxelev
      DO j=1,numrad
        azmvol(j,k) = j*1.0
        elvvol(j,k) = elvang(k)
      ENDDO
      kntrazm(k) = numrad2

    azmvol(numrad2-1,k) = 0.0
    elvvol(numrad2-1,k) = elvvol(numrad,k)
    azmvol(numrad2,k) = 1.0
    elvvol(numrad2,k) = elvvol(numrad,k)
    DO idp=1,numvar
      DO i=1,ngate
        rdr_p(i,numrad+1,k,idp) = rdr_p(i,numrad,k,idp)
        rdr_p(i,numrad+2,k,idp) = rdr_p(i,1,k,idp)
      ENDDO
    ENDDO

  ENDDO

  DO idp = 1, numvar

    nelev = rdrknt_tot(idp)

    IF( nelev > 0 ) THEN
      WRITE(6,'(a,a)') 'Calling remap_arpsppi for ', dualvarname(idp)
      CALL remap_arpsppi(ngate,numrad2,nelev,nx,ny,kntrazm,maxelev,    &
         radlat,radlon,radarx,radary,radalt,dazim,0,rngmax,            &
         gatesp,azmvol,elvvol,rdr_p(:,:,:,idp),ref_t,xs,ys,            &
         gridtilt_height,gridtilt_range,gridtilt_slr,gridtilt_azm,     &
         gridtilt_rdr(:,:,:,idp),gridtilt_time,gridtilt_rdr_mean)
    ELSE
      WRITE(6,'(i4,3a)') nelev,' Levels for ',dualvarname(idp),       &
            ' not calling remap_arpsppi'
    END IF

  END DO

  full_fname = ' '
  ref_sec = 0
  IF(tintrpopt > 0 .AND. nlistfil > 1) THEN
    WRITE(full_fname,'(a,a4,a1,i4.4,2(i2.2),a1,3(i2.2))')              &
          dirname(1:ldirnam),radname,'.',ref_yr,ref_mon,ref_day,'.',   &
          ref_hr,ref_min,ref_sec
    itimfrst = ref_t
  ELSE
    WRITE(full_fname,'(a,a4,a1,i4.4,2(i2.2),a1,3(i2.2))')              &
          dirname(1:ldirnam),radname,'.',iiyear,iimon,iiday,'.',       &
          iihr, iimin, iisec
  ENDIF
  PRINT*, "Now writing in EnKF format"
  PRINT*, trim(full_fname)

  IF(myproc == ROOT .AND. gridtilt == 1) THEN
     CALL wrtgridtilt(full_fname,fntimopt,                             &
            maxelev,nx,ny,nz,kntrelv,radname,                          &
            radlat,radlon,radarx,radary,radalt,dazim,rngmin,rngmax,    &
            gridtilt_height,gridtilt_range,gridtilt_slr,gridtilt_azm,  &
            gridtilt_rdr,gridtilt_time,gridtilt_time,itimfrst,elvang,  &
            gridtilt_height,gridtilt_range,gridtilt_slr,gridtilt_azm,  &
            elvang,xs,ys,zps,ivcp,isource,grdtiltver,numvar)
  END IF

  DEALLOCATE(azmvol,elvvol,gridtilt_height,gridtilt_range,gridtilt_slr,&
             gridtilt_azm,gridtilt_rdr,gridtilt_time,gridtilt_rdr_mean)

!------------------------------------------------------------------------
!
! The End.
!
!------------------------------------------------------------------------

    IF (myproc == ROOT)  &
      WRITE(stdout,'(/a/)') '  === Normal termination of NCRAD2ARPS === '

    CALL mpexit(0)

    STOP

  END IF ! gridtilt

END PROGRAM ncrad2arps

SUBROUTINE snrchk(maxgate,maxazim,ngate,nazim,refchek,rfirstg,gatesp,   &
                   snrthr,snstvty,refl,radv,attref)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: nazim
  REAL, INTENT(IN) :: refchek
  REAL, INTENT(IN) :: rfirstg
  REAL, INTENT(IN) :: gatesp
  REAL, INTENT(IN) :: snrthr
  REAL, INTENT(OUT) :: snstvty(maxgate)
  REAL, INTENT(INOUT) :: refl(maxgate,maxazim)
  REAL, INTENT(INOUT) :: radv(maxgate,maxazim)
  REAL, INTENT(IN) :: attref(maxgate,maxazim)
!
  REAL, PARAMETER :: refmiss=-911.0
  REAL, PARAMETER :: velmiss=-911.0
  INTEGER :: igate,jazim
  INTEGER :: kntvalid,kntflag
  REAL :: snr,pctflg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL xbsnstvty(maxgate,rfirstg,gatesp,snstvty)
  DO igate=1,maxgate,10
    print *, ' igate: ',igate,' range=',(0.001*(rfirstg+(igate-1)*gatesp)), &
             ' sens= ',snstvty(igate)
  END DO

  kntvalid=0
  kntflag=0
  DO jazim=1,nazim
    DO igate=1,ngate
      IF(refl(igate,jazim) > refchek ) THEN
        kntvalid=kntvalid+1
        IF( (attref(igate,jazim)-snstvty(igate)) < snrthr ) THEN
          kntflag=kntflag+1
          refl(igate,jazim) = refmiss
          radv(igate,jazim) = velmiss
        END IF
      END IF
    END DO
  END DO

  pctflg=0.
  IF(kntvalid > 0) THEN
    pctflg=100.*float(kntflag)/float(kntvalid)
  END IF
  WRITE(6,'(a,i9,a,i9,a,f5.1,a)') &
   ' SNR screened ',kntflag,' of ',kntvalid,' = ',pctflg, &
   '% reflectivity gates'

  RETURN
END SUBROUTINE snrchk

SUBROUTINE xbsnstvty(ngate,rfirstg,gatesp,sens)
!
!-----------------------------------------------------------------------
! Calculate X-band radar sensitivity (dBZ) as a function of range.
!
! Based on information provided by Francesc Joyvount
!
! Keith Brewster and Erin Fay, CAPS
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngate
  REAL, INTENT(IN) :: rfirstg
  REAL, INTENT(IN) :: gatesp
  REAL, INTENT (OUT):: sens(ngate)     ! Sensitivity values for all dist from radar.
!
! Parameters for CASA radar
!
!  REAL, PARAMETER :: lmbda = 3.0E-02 ! Wavelength (m)
  REAL, PARAMETER :: freqkhz=11505000.0  ! kHz
  REAL, PARAMETER :: ptxdbm=66.68    ! Peak power (dBm)
!  REAL, PARAMETER :: pt = 12500     ! Peak power (Watts)
  REAL, PARAMETER :: G  = 38.0       ! Antenna gain (dB)
  REAL, PARAMETER :: F  = 5.5        ! Noise figure (dB)
!  REAL, PARAMETER :: tau = 0.67E-06  ! Radar pulse length (s)
  REAL, PARAMETER :: tau = 4.06E-07  ! Radar pulse length (s)

  REAL, PARAMETER :: theta = 1.8     ! Antenna half-power beamwidth (deg)
  REAL, PARAMETER :: B = 2.5         ! Bandwidth (MHz)
  REAL, PARAMETER :: lm = 1.0        ! Receiver mis-match loss(dB)
  REAL, PARAMETER :: Kw2 = 0.91      ! Refractive Index of water squarred
  REAL, PARAMETER :: T0 = 300.       ! Rx temperature (K)
  REAL, PARAMETER :: Ta = 200.       ! Antenna temperature (K)
!
! Physical constants
!
  REAL, PARAMETER :: K = 1.38E-23    ! Boltsmann's Constant (J/K)
  REAL, PARAMETER :: c = 2.99792E08  ! Speed of light (m/s)
  REAL, PARAMETER :: rconst =1.0E18  ! m6/m3 to mm6/m3
!
! Misc internal variables
!
  INTEGER :: igate
  REAL :: ln2,pi,pi3,four3,bwrad,rnoise,BHz,Ni,sconst,rlm,rG
  REAL :: pt,lmbda
  REAL :: range
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  pt=0.001*(10.0**(0.1*ptxdbm))
  print *, ' pt= ',pt
  lmbda=c/(freqkhz*1000.)
  print *, ' lmbda (m)= ',lmbda
  ln2=log(2.0)
  pi=acos(-1.)
  pi3=pi*pi*pi
  four3=4.*4.*4.
  bwrad=theta*pi/180.
  rnoise=10.**(0.1*F)
  rlm=10.**(0.1*lm)
  rG=10.**(0.1*G)
  BHz=1.0E06*b
  Ni=K*(Ta+(rnoise-1.0)*T0)*BHz
  sconst=rconst *(Ni/(pi3*Kw2)) * ((8.*ln2)/(bwrad*bwrad)) *         &
         (2./(c*tau)) * ((lmbda*lmbda*four3*rlm)/(Pt*rG*rG))
  DO igate = 1, ngate
     range=rfirstg+((igate-1)*gatesp)
     sens(igate) = 10.*LOG10(range*range*sconst)
  END DO
END SUBROUTINE xbsnstvty
SUBROUTINE extdirn(filename,dirname,fname)
!
! Extract directory name from complete filename and output
! directory name and filename in two character string variables
!
  IMPLICIT NONE

  CHARACTER (LEN=*) :: filename
  CHARACTER (LEN=*) :: dirname
  CHARACTER (LEN=*) :: fname

  INTEGER :: locslash,lenf
  LOGICAL, PARAMETER :: back = .TRUE.

  lenf=LEN_TRIM(filename)
  dirname='./'
  fname=filename
  locslash=INDEX(filename,'/',back)
  IF(locslash > 0) THEN
    dirname=filename(1:(locslash-1))
    fname=filename((locslash+1):lenf)
  END IF

  RETURN
END SUBROUTINE extdirn

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE RADIALINTP_RDR                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                 Unxversity of Oklahoma               ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE RadialIntp_rdr(ngate,numrad,nazim,azim,vars,rdr_tmp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Linear spatial interpolation algorithm for radar in polar
!  coordinates. This is necessary for time interpolation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    09/04/2009
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ngate, numrad, nazim
  REAL :: azim(nazim)
  REAL :: vars(ngate,nazim)
  REAL :: rdr_tmp(ngate,numrad)

! Local variables
  REAL, PARAMETER :: missing = -999.9
  REAL, PARAMETER :: missing_Foray = -99.0
  INTEGER :: i,j
  INTEGER :: flip,jk,jk1,jk2,jj,jj1,jj2
  REAL :: azim2l,azim2r,azimlr

  REAL, ALLOCATABLE :: temp1(:,:)
  REAL, ALLOCATABLE :: temp2(:)
  INTEGER :: istat, k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE (temp1(ngate,nazim), stat = istat)
  CALL check_alloc_status(istat,'ncrad2arps:temp1')
  ALLOCATE (temp2(nazim), stat = istat)
  CALL check_alloc_status(istat,'ncrad2arps:temp2')
  temp1 = missing_Foray; temp2 = missing_Foray
  rdr_tmp = missing

! revert the radial if they are counterclockwise.
  IF(azim(2) > azim(3) .and. azim(3) > azim(4)) THEN
    DO j= 1,nazim
      k = (nazim+1)-j
      DO i = 1,ngate
        temp1(i,k) = vars(i,j)
      END DO
      temp2(k)   = azim(j)
    END DO
  ELSE
    temp1 = vars
    temp2 = azim
  END IF

  flip = 0
  DO j = 1,nazim
    jk1 = j
    IF(j == nazim) THEN
      jk2 = 1
    ELSE
      jk2 = jk1+1
    ENDIF
    IF( temp2(jk2) < temp2(jk1) ) THEN
      flip = 1
      temp2(jk2) = temp2(jk2) + 360.0
    ENDIF

    jj1 = int(temp2(jk1)+1.0)
    jj2 = int(temp2(jk2))

    IF(flip == 1) THEN
      temp2(jk2) = temp2(jk2) - 360.0
      flip = 0
    ENDIF

    IF(int(temp2(jk1)) > jj1 .or. int(temp2(jk2)) /= jj1) THEN
      GOTO 345
    ENDIF

    DO jj = jj1,jj2
      jk = jj
      IF(jk > 360) jk=jk-360
      azim2l = jj - temp2(jk1)
      azim2r = temp2(jk2) - jj
      azimlr = temp2(jk2) - temp2(jk1)

      DO i = 1,ngate
        IF(temp1(i,jk1) > missing_Foray .and. temp1(i,jk2) > missing_Foray) THEN
          rdr_tmp(i,jk)=(temp1(i,jk1)*azim2r &
                              +temp1(i,jk2)*azim2l)/azimlr
        ELSE
          rdr_tmp(i,jk)=missing
        ENDIF
      ENDDO
    ENDDO
345 CONTINUE
  ENDDO

  DEALLOCATE(temp1,temp2)

END SUBROUTINE RadialIntp_rdr

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE TIMEINTP_RDR                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                 Unxversity of Oklahoma               ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE TimeIntp_rdr(ngate,numrad,scantime,ref_t,rdr_tmp,rdr_p)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Linear time interpolation algorithm for radar in polar coordinates.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    09/04/2009
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ngate, numrad
  INTEGER :: scantime(2),ref_t
  REAL :: rdr_tmp(ngate,numrad,2)
  REAL :: rdr_p(ngate,numrad)

! Local variables
  REAL, PARAMETER :: missing = -999.9
  INTEGER :: i,j,stdout
  INTEGER :: time2p,time2f,timepf
  INTEGER :: iyear,imonth,iday,ihr,imin,isec

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  time2p = ref_t - scantime(1)
  time2f = scantime(2) - ref_t
  timepf = scantime(2) - scantime(1)
  WRITE(stdout,'(3i20)') scantime(1),ref_t,scantime(2)
  CALL abss2ctim(scantime(1), iyear, imonth, iday, ihr, imin, isec)
  WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO') &
         '    DATE: ', iyear, '/', imonth, '/', iday
  WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')  &
         '    TIME: ', ihr, ':', imin, ':', isec
  CALL abss2ctim(ref_t, iyear, imonth, iday, ihr, imin, isec)
  WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO') &
         '    DATE: ', iyear, '/', imonth, '/', iday
  WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')  &
         '    TIME: ', ihr, ':', imin, ':', isec
  CALL abss2ctim(scantime(2), iyear, imonth, iday, ihr, imin, isec)
  WRITE(stdout,'(1x,a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO') &
         '    DATE: ', iyear, '/', imonth, '/', iday
  WRITE(stdout,'(1X,a,i2.2,a1,i2.2,a1,i2.2)')  &
         '    TIME: ', ihr, ':', imin, ':', isec

  DO i = 1,ngate
    Do j = 1,numrad
      IF (rdr_tmp(i,j,1) > missing .and. rdr_tmp(i,j,2) > missing) THEN
        rdr_p(i,j) = (rdr_tmp(i,j,1)*time2f  &
                               +rdr_tmp(i,j,2)*time2p)/timepf
      ELSE
        rdr_p(i,j)=missing
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE TimeIntp_rdr
