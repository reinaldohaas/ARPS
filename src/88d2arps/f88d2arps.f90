!########################################################################
!########################################################################
!#########                                                      #########
!#########                  PROGRAM f88d2arps                   #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

PROGRAM f88d2arps

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar data from NEXRAD 88D data files and remaps data on ARPS
! grid to be used for analysis or display.
!
! Fortran version that replicates the function of the original 88d2arps.c,
! developed in the early-90s with updates through 2009.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (July, 2009)
!
! MODIFICATIONS:
!
! Keith Brewster, February, 2010
! Updated logic to provide balanced load sharing on MPI.
!
! Keith Brewster, July, 2010
! Some clean-up after testing for acceptance in official ARPS release.
!
! Keith Brewster
! Tanslated the following Y. Jung's addition to 88d2arps.c to this version
!   Added wtradtiltcol subroutine to store radar data in the colume on
!   radar tilts. To be used in the EnKF analysis.
!
! Yunheng Wang (04/08/2011)
! Fixed solo support.
!
! Keith Brewster (08/22/2011)
! Added vadsrc to wtwadprf subroutine call.
!
! Y. Wang (05/08/2012)
! Replace the non-standard calls of getarg with the Fortran 2003
! standard call of GET_COMMAND_ARGUMENT.
!
! Keith Brewster (09/07/2012)
! Added processing of dual-pol variables, including new input file
! options to control dual-pol processing and output.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'remapcst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Dimensions and parameters
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxvar = 10

  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: iordunf = 2
  INTEGER, PARAMETER :: itime1970 = 315619200

  CHARACTER(LEN=6)  :: velid = 'radv3d'
  CHARACTER(LEN=6)  :: refid = 'refl3d'
  CHARACTER(LEN=6)  :: rhvid = 'rhv_3d'
  CHARACTER(LEN=6)  :: zdrid = 'zdr_3d'
  CHARACTER(LEN=6)  :: kdpid = 'kdp_3d'
  CHARACTER(LEN=20) :: velname = 'RawVelocity'
  CHARACTER(LEN=20) :: refname = 'Reflectivity'
  CHARACTER(LEN=20) :: rhvname = 'Correlation H-V'
  CHARACTER(LEN=20) :: zdrname = 'Zdr'
  CHARACTER(LEN=20) :: kdpname = 'Kdp'
  CHARACTER(LEN=20) :: velunits = 'm/s'
  CHARACTER(LEN=20) :: refunits = 'dBZ'
  CHARACTER(LEN=20) :: rhvunits = ' '
  CHARACTER(LEN=20) :: zdrunits = 'dBZ'
  CHARACTER(LEN=20) :: kdpunits = 'deg/km'

  CHARACTER(LEN=20) :: pname
  INTEGER           :: plen

!-----------------------------------------------------------------------
!
! Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx, ny, nz, nzsoil, nstyps
  REAL    :: radalt, radlat, radlon
  REAL    :: radarx, radary

!-----------------------------------------------------------------------
!
! Radar Data Variables for one tilt
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: refl(:,:)
  REAL, ALLOCATABLE :: radv(:,:)
  REAL, ALLOCATABLE :: spw(:,:)
!
!-----------------------------------------------------------------------
!
! Variables for VAD profiles.
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: vadhgt(:)
  REAL, ALLOCATABLE :: vaddir(:)
  REAL, ALLOCATABLE :: vadspd(:)
!
!------------------------------------------------------------------------
!
! Variables for remapping routines.
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
  REAL, ALLOCATABLE :: rhvdata(:,:) ! Rho-HV data
  REAL, ALLOCATABLE :: zdrdata(:,:) ! Zdr data
  REAL, ALLOCATABLE :: phidata(:,:) ! Phi data
  REAL, ALLOCATABLE :: kdpdata(:,:) ! Kdp data
  REAL, ALLOCATABLE :: veldata(:,:) ! velocity data (m/s)
  REAL, ALLOCATABLE :: spwdata(:,:)
  REAL, ALLOCATABLE :: bkgvel(:,:)
  REAL, ALLOCATABLE :: unfvdata(:,:) ! unfolded velocity data (m/s)
  REAL, ALLOCATABLE :: snstvty(:) ! radar sensitivity
  REAL, ALLOCATABLE :: rtem(:,:)  ! temporary array for radar data processing
  INTEGER, ALLOCATABLE :: time(:) ! time offset from itimfrst
  INTEGER, ALLOCATABLE :: gcfst(:) ! Ground clutter filter state
  REAL, ALLOCATABLE :: refrange(:)   ! reflectivity range (km)
  REAL, ALLOCATABLE :: azim(:)    ! azimuth angle for each radial(degree)
  REAL, ALLOCATABLE :: elev(:)    ! elevation angle for each radial (degree)
  REAL, ALLOCATABLE :: vnyq(:)    ! Nyquist velocities (meters/sec)
  INTEGER, ALLOCATABLE :: istrgate(:)
  INTEGER, ALLOCATABLE :: bgate(:)
  INTEGER, ALLOCATABLE :: egate(:)

  INTEGER :: itimfrst            ! time of first radial in volume
  INTEGER :: rfrst_ref           ! range to first gate (meters)
  INTEGER :: gtspc_ref           ! gate spacing (meters)
  INTEGER :: rfrst_vel           ! range to first gate (meters)
  INTEGER :: gtspc_vel           ! gate spacing (meters)

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
  REAL, ALLOCATABLE :: elvmean(:)
  INTEGER, ALLOCATABLE :: kntbin(:)
  INTEGER, ALLOCATABLE :: ngateref(:)
  INTEGER, ALLOCATABLE :: ngatevel(:)

  REAL, ALLOCATABLE :: rxvol(:,:,:)
  REAL, ALLOCATABLE :: ryvol(:,:,:)
  REAL, ALLOCATABLE :: rzvol(:,:,:)

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
  INTEGER, ALLOCATABLE :: icolp(:)
  INTEGER, ALLOCATABLE :: jcolp(:)
  REAL, ALLOCATABLE :: xcolp(:)
  REAL, ALLOCATABLE :: ycolp(:)
  REAL, ALLOCATABLE :: zcolp(:,:)
  LOGICAL, ALLOCATABLE :: havdat(:,:)

  REAL, ALLOCATABLE :: colvel(:,:)
  REAL, ALLOCATABLE :: colref(:,:)
  REAL, ALLOCATABLE :: colrhv(:,:)
  REAL, ALLOCATABLE :: colzdr(:,:)
  REAL, ALLOCATABLE :: colkdp(:,:)
  REAL, ALLOCATABLE :: colnyq(:,:)
  REAL, ALLOCATABLE :: coltim(:,:)

  REAL, ALLOCATABLE :: x(:)
  REAL, ALLOCATABLE :: y(:)
  REAL, ALLOCATABLE :: z(:)
  REAL, ALLOCATABLE :: xslg(:)
  REAL, ALLOCATABLE :: yslg(:)
  REAL, ALLOCATABLE :: zpslg(:,:,:)

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

  REAL, ALLOCATABLE :: zpsoil(:,:,:)
  REAL, ALLOCATABLE :: j3soil(:,:,:)
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)

  REAL, ALLOCATABLE :: u(:,:,:)
  REAL, ALLOCATABLE :: v(:,:,:)
  REAL, ALLOCATABLE :: qv(:,:,:)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)
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

  REAL(KIND=8), ALLOCATABLE :: tem_double(:)

  REAL, ALLOCATABLE :: grdtiltdata(:,:,:,:)

  REAL, ALLOCATABLE :: grdtilthighref(:,:,:)
  REAL, ALLOCATABLE :: grdrangeref(:,:,:)
  REAL, ALLOCATABLE :: grdslrref(:,:)
  REAL, ALLOCATABLE :: grdazmref(:,:)
  REAL, ALLOCATABLE :: elevmeanref(:)

  INTEGER, ALLOCATABLE :: tilttimeref(:)
  INTEGER, ALLOCATABLE :: tilttimevel(:)

  REAL, ALLOCATABLE :: grdtilthighvel(:,:,:)
  REAL, ALLOCATABLE :: grdrangevel(:,:,:)
  REAL, ALLOCATABLE :: grdslrvel(:,:)
  REAL, ALLOCATABLE :: grdazmvel(:,:)
  REAL, ALLOCATABLE :: elevmeanvel(:)

!------------------------------------------------------------------------
!
! Radar data dimensions
!
!------------------------------------------------------------------------

  INTEGER :: maxrgate,maxvgate,maxgate,maxazim,maxelev

!------------------------------------------------------------------------
!
! Misc. variables
!
!------------------------------------------------------------------------

  CHARACTER(LEN=256) :: full_fname
  CHARACTER(LEN=256) :: vad_fname
  CHARACTER(LEN=266) :: syscall
  CHARACTER(LEN=256) :: errmessage

  INTEGER :: n,i,j,k,igate,istat,istatus
  INTEGER :: refstat,rhvstat,dlpstat,velstat,eofstat,rdstat
  INTEGER :: nxny,nazim
  INTEGER :: ivcp,itimcdf,initime,mintimcdf,ivar
  INTEGER :: len_fname
  INTEGER :: arbvaropt,velopt,varfill,dualpol
  INTEGER :: iopt
  INTEGER :: nprocx_in_cmd,nprocy_in_cmd

  LOGICAL :: velproc
  LOGICAL :: dualpdata,dualpproc
  LOGICAL :: vel2d,vel3d
  LOGICAL :: ref2d,ref3d
  LOGICAL :: rhv2d,rhv3d
  LOGICAL :: zdr2d,zdr3d
  LOGICAL :: kdp2d,kdp3d
  LOGICAL :: unfdiag
  LOGICAL :: qc_on
  LOGICAL :: fntime
  LOGICAL :: traditional
  LOGICAL :: rtopts
  LOGICAL :: vad
  LOGICAL :: wrttilt,wrtgrdtilt
  LOGICAL :: wrtsolo,wrtsoloqc
  LOGICAL :: remap_on

  LOGICAL :: clrvcp
  LOGICAL :: rho_match

  INTEGER, PARAMETER :: iradfmt = 1
  CHARACTER (LEN=8), PARAMETER :: vadsrc = '88D VAD'

  INTEGER :: iyear,imon,iday,ihr,imin,isec
  INTEGER :: kyear,kmon,kday,khr,kmin,ksec
  INTEGER :: kmyear,kmmon,kmday,kmhr,kmmin,kmsec
  INTEGER :: iiyear,iimon,iiday,iihr,iimin,iisec
  INTEGER :: isource,icount,iradtype
  INTEGER :: iscan,itilt,iangle,kref,kvel,ktime
  INTEGER :: kfn,kkfn
  INTEGER :: nyqset,timeset,vardump,dlpmedfilt
  INTEGER :: ng_ref,ng_vel
  INTEGER :: ncoltot,ncolp
  INTEGER :: xscale

  REAL :: rdx,dsort,dtime,rmin,rmax
  REAL :: rmisval,rngfval,frtime
  REAL :: refelvmin,refelvmax
  REAL :: refrngmin,refrngmax
  REAL :: badpct
  REAL, PARAMETER :: phichek = -400.

  REAL :: radius,vconst,vnyquist

  REAL :: xrad,yrad
  REAL :: tmin,tmax

  INTEGER :: maxinfile,infile,lenfn,iloc,jloc,jazim,kelv
  INTEGER :: iarg,jarg,iargc,narg,icol

  CHARACTER (LEN=256) :: charg
  CHARACTER (LEN=6)   :: varid
  CHARACTER (LEN=20)  :: varname
  CHARACTER (LEN=20)  :: varunits

  CHARACTER (LEN=256) :: tempstr

  INTEGER, PARAMETER :: ROOT = 0
  REAL, PARAMETER :: anrflag = -800.
  LOGICAL :: back = .TRUE.
  REAL    :: rtime

  INTEGER :: numvar

  INTEGER, EXTERNAL :: write_sweep

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
!
  refelvmin = 91.0
  refelvmax = -91.0
  refrngmin = 5.0E03
  refrngmax = 230.0E03
  dmpfmt=1
  hdf4cmpr=0
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
  vad = .FALSE.
  wrttilt = .FALSE.
  wrtgrdtilt = .FALSE.
  wrtsolo = .FALSE.
  wrtsoloqc = .FALSE.
  clrvcp = .FALSE.
  grdtiltver = 2
  numvar = 2
  dlpmedfilt = 0

  dsort = (sortmax-sortmin)/float(nsort-1)

  radname = 'DMMY'
  radband = 1
  myproc = 0
  mp_opt = 0
  nprocx_in_cmd = 1
  nprocy_in_cmd = 1

  CALL mpinit_proc(0)
!
!------------------------------------------------------------------------
!
! First check to see if help command line argument is specified.
!
!------------------------------------------------------------------------
!
  istatus=0
  IF(myproc == ROOT) THEN
    narg = COMMAND_ARGUMENT_COUNT()
    IF(narg > 0 ) THEN
      CALL GET_COMMAND_ARGUMENT(1, charg, istat, istatus )
      IF(charg(1:5) == '-help' .OR. charg(1:6) == '--help') THEN
        CALL GET_COMMAND_ARGUMENT(0, charg, istat, istatus )
        WRITE(stdout,'(a,a,a,a,a,a)') &
        ' Usage: ',TRIM(charg),' [RADAR_ID] [fname] [-novel] [-dir dirname] [-hdf 2]',&
        ' [-binary] [-ref2d] [-ref3d] [-reffile] [-vel2d] [-vel3d] [-velfile]',&
        ' [-noqc] [-fillref] [-fntime] [-rtopts] [-radar98] [-input nmfile]', &
        ' < radremap.input'
        istatus = 1
      END IF
    END IF
  END IF

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) CALL arpsstop (' ',istatus)

!------------------------------------------------------------------------
!
! Go through the arguments to extract "nprocx_in" and "nprocy_in"
! if they are present.  This must be done before call to initremapopt.
!
!------------------------------------------------------------------------
!
  tempstr = ' '
  IF(myproc == ROOT) THEN

    IF(narg > 0) THEN
      iarg=1
      DO jarg=1,narg
        CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
        IF(charg(1:10) == '-nprocx_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            READ(charg,'(i)') nprocx_in_cmd
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (X) = ',nprocx_in_cmd
        ELSE IF(charg(1:10) == '-nprocy_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            READ(charg,'(i)') nprocy_in_cmd
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (Y) = ',nprocy_in_cmd
        ELSE IF(charg(1:6) == '-input') THEN    ! namelist file
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            tempstr = charg
          END IF
          WRITE(stdout,'(1x,2a)') 'Namelist file from command line - ',TRIM(tempstr)
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
  CALL initremapopt(nx,ny,nz,nzsoil,nstyps,                             &
                    nprocx_in_cmd,nprocy_in_cmd,0,                      &
                    tempstr,istatus)

  IF(myproc == ROOT) THEN

    ref2d = (ref2dopt > 0)
    ref3d = (ref3dopt > 0)
    vel2d = (vel2dopt > 0)
    vel3d = (vel3dopt > 0)
    zdr2d = (zdr2dopt > 0)
    IF( zdr2d ) dualpproc = .TRUE.
    zdr3d = (zdr3dopt > 0)
    IF( zdr3d ) dualpproc = .TRUE.
    kdp2d = (kdp2dopt > 0)
    IF( kdp2d ) dualpproc = .TRUE.
    kdp3d = (kdp3dopt > 0)
    IF( kdp3d ) dualpproc = .TRUE.
    rhv2d = (rhv2dopt > 0)
    IF( rhv2d ) dualpproc = .TRUE.
    rhv3d = (rhv3dopt > 0)
    IF( rhv3d ) dualpproc = .TRUE.
    velproc = (velprocopt > 0)
    IF (dualpout > 0) dualpproc = .TRUE.
    unfdiag = (unfdiagopt > 0)
    qc_on = (qcopt > 0)
    traditional = (tradopt > 0)
    vad = (vadopt > 0)
    wrttilt = (wttiltopt > 0)
    wrtgrdtilt = (wtgrdtiltopt > 0)
    jarg = MOD(wtsoloopt,100)
    wrtsolo   = ( jarg == 1 .OR. jarg == 3 )
    wrtsoloqc = ( jarg == 2 .OR. jarg == 3 )
    fntime = (fntimopt > 0)

    WRITE(6,'(a,i4)') ' Obtained from namelists: dmpfmt=',dmpfmt
    WRITE(6,'(26x,a,i4)') 'hdf4cmpr=',hdf4cmpr

!------------------------------------------------------------------------
!
! Process command line
! For backward comptability allow input of radar file names via
! files specified on the command line
!
!------------------------------------------------------------------------

    IF(narg > 0) THEN
      iarg=1
      DO jarg=1,narg
        CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
        IF(iarg == 1 .AND. charg(1:1) /= '-') THEN
          radname=charg(1:4)
          WRITE(stdout,'(a,a)')  '    radname = ', radname
        ELSE IF(charg(1:6) == '-input') THEN    ! namelist file
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            tempstr = charg
          END IF
          WRITE(stdout,'(1x,2a)') 'Namelist file from command line - ',TRIM(tempstr)
        ELSE IF(charg(1:6) == '-novel') THEN
          velproc=.FALSE.
        ELSE IF(charg(1:4) == '-hdf') THEN
          dmpfmt=3
          hdf4cmpr=0
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i1)',iostat=istat) hdf4cmpr
            IF(istat == 0) THEN
              hdf4cmpr=min(max(hdf4cmpr,0),7)
            ELSE
              WRITE(stdout,'(a,/a,a)')                                      &
                 '-hdf flag requires compression level 0-7 to follow',      &
                 ' found: ',TRIM(charg)
              errmessage = ' Command line argument error'
              istatus = -1
            END IF
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
            ' Produce a 2d reflectivity file of lowest tilt'
        ELSE IF(charg(1:6) == '-ref3d') THEN
          ref3d=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d reflectivity file for plotting'
        ELSE IF(charg(1:8) == '-reffile') THEN
          ref3d=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d reflectivity file for plotting'
        ELSE IF(charg(1:6) == '-vel2d') THEN
          vel2d=.TRUE.
          velproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 2d velocity file of lowest tilt'
        ELSE IF(charg(1:6) == '-vel3d') THEN
          vel3d=.TRUE.
          velproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d velocity file for plotting'
        ELSE IF(charg(1:8) == '-velfile') THEN
          vel3d=.TRUE.
          velproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d velocity file for plotting'
        ELSE IF(charg(1:6) == '-rhv2d') THEN
          rhv2d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 2d Rho-HV file of lowest tilt'
        ELSE IF(charg(1:6) == '-rhv3d') THEN
          rhv3d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d Rho-HV file for plotting'
        ELSE IF(charg(1:6) == '-zdr2d') THEN
          zdr2d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 2d Zdr file of lowest tilt'
        ELSE IF(charg(1:6) == '-zdr3d') THEN
          zdr3d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d Zdr file for plotting'
        ELSE IF(charg(1:6) == '-kdp2d') THEN
          kdp2d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 2d Kdp file of lowest tilt'
        ELSE IF(charg(1:6) == '-kdp3d') THEN
          kdp3d=.TRUE.
          dualpproc=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Produce a 3d Kdp file for plotting'
        ELSE IF(charg(1:8) == '-unfdiag') THEN
          unfdiag=.TRUE.
          WRITE(stdout,'(a)')                                               &
            ' Write files for velocity unfolding diagnostics'
        ELSE IF(charg(1:5) == '-noqc') THEN
          qc_on=.FALSE.
          WRITE(stdout,'(a)')                                               &
            ' Skip QC steps in processing'
        ELSE IF(charg(1:8) == '-fillref') THEN
          fillref=1
          WRITE(stdout,'(a)')                                               &
            ' Fill-in reflectivity below lowest tilt'
        ELSE IF(charg(1:8) == '-medfilt') THEN
          medfilt=1
          WRITE(stdout,'(a)')                                               &
            ' Apply median filter to az-ran data in despeckle step.'
        ELSE IF(charg(1:10) == '-nprocx_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i)') nprocx_in
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (X) = ',nprocx_in
        ELSE IF(charg(1:10) == '-nprocy_in') THEN
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i)') nprocy_in
          END IF
          WRITE(stdout,'(a,i4)') ' Number patches of data (Y) = ',nprocy_in

        ELSE IF(charg(1:4) == '-dir') THEN
          WRITE(stdout,'(a)') ' Use specified directory for output.'
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, dirname)
            WRITE(stdout,'(a,a)') '  directory name: ',TRIM(dirname)
          ELSE
            WRITE(stdout,'(a,a,a)') ' Switch ',TRIM(charg), &
                                    ' requires an argument, try again'
            errmessage = ' Command line argument error'
            istatus    = -1
          END IF
        ELSE IF(charg(1:7) == '-fntime') THEN
          fntime=.TRUE.
          WRITE(stdout,'(a)') ' Use time in filename for output filename.'

        ELSE IF(charg(1:7) == '-rtopts') THEN
          rtopts=.TRUE.
          clrvcpopt=0
          medfilt=1
          WRITE(stdout,'(a)') ' Use real-time options.'
          WRITE(stdout,'(a)') ' Do not process any data from clear-air VCPs'
          WRITE(stdout,'(a)') ' Apply median filter.'

        ELSE IF(charg(1:8) == '-noclrv') THEN
          clrvcpopt=1
          WRITE(stdout,'(a)') ' Skipping velocities in clear-air VCP data.'

        ELSE IF(charg(1:8) == '-noclear') THEN
          clrvcpopt=0
          WRITE(stdout,'(a)') ' Skipping all clear-air VCP data.'

        ELSE IF(charg(1:8) == '-rad98') THEN
          rad98opt=1
          WRITE(stdout,'(a)') ' Read using WSR-98D CINRAD radar format.'

        ELSE IF (charg(1:5) == '-solo') THEN
          WRITE(stdout,'(a)') &
            ' Write sweep file in DORADE format for working with SoloII'
          wrtsolo=.TRUE.

        ELSE IF (charg(1:7) == '-soloqc') THEN
          WRITE(stdout,'(a)') &
            ' Write sweep file in DORADE format for working with SoloII after QC'
          wrtsoloqc = .TRUE.

        ELSE IF (charg(1:6) == '-sound') THEN
          WRITE(stdout,'(a)') ' Use single-sounding backgroud file'
          initopt=1
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, sndfile )
            WRITE(stdout,'(a,a)') ' Using sounding file ',TRIM(sndfile)
          ELSE
            WRITE(stdout,'(a,a,a)') ' Switch ',TRIM(charg), &
                                    ' requires an argument, try again'
            errmessage = ' Command line argument error'
            istatus=-1
          END IF

        ELSE IF(charg(1:6) == '-tradi') THEN
          traditional = .TRUE.
          WRITE(stdout,'(a)') ' Force use of traditional backgroud processing (initgrdvar).'

        ELSE IF (charg(1:6) == '-vad') THEN
          WRITE(stdout,'(a)') ' VAD wind profile file generated for ADAS'
          vad =.TRUE.

        ELSE IF (charg(1:8) == '-wrttilt') THEN
          WRITE(stdout,'(a)') ' Write radar scans in az-ran format.'
          wrttilt =.TRUE.

        ELSE IF (charg(1:12) == '-wrtgridtilt') THEN
          WRITE(stdout,'(a)') ' Write radar scans in grid-tilt format.'
          wrtgrdtilt =.TRUE.
          IF( iarg < narg ) THEN
            CALL GET_COMMAND_ARGUMENT((iarg+1), charg )
            READ(charg,'(i)',iostat=rdstat) iopt
            IF(rdstat == 0 .AND. iopt < 2) THEN
              iarg=iarg+1
              grdtiltver=1
            END IF
          END IF
          IF( grdtiltver == 1 ) THEN
            WRITE(stdout,'(a)') ' Old grid-tilt version for the EnKF application'
            WRITE(stdout,'(a)') ' Data are stored in the entire ARPS domain'
          ELSE
            WRITE(stdout,'(a)') ' Data are stored in updated grid tilt-column format'
            WRITE(stdout,'(a)') ' For old version use -wrtgridtilt 1'
          END IF
        ELSE
          radfname=charg
          WRITE(stdout,'(a,a)')  ' Radar data file name: ',TRIM(radfname)

        END IF
        iarg=iarg+1
        IF(iarg > narg) EXIT

      END DO
    END IF ! narg > 0
  END IF ! myproc == 0

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) THEN
    CALL arpsstop(TRIM(errmessage),istatus)
  END IF
!
!------------------------------------------------------------------------
!
! After command line and input variables have been read,
!   communicate updated option parameters with all processors.
!
!------------------------------------------------------------------------
!
  CALL mpupdatec(radname,4)
  CALL mpupdatel(velproc,1)
  CALL mpupdatel(dualpproc,1)
  CALL mpupdatei(dualpout,1)
  CALL mpupdatei(dmpfmt,1)
  CALL mpupdatei(hdf4cmpr,1)
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
  CALL mpupdatel(qc_on,1)
  CALL mpupdatei(fillref,1)
  CALL mpupdatel(fntime,1)
  CALL mpupdatel(traditional,1)
  CALL mpupdatel(wrttilt,1)
  CALL mpupdatel(wrtgrdtilt,1)
  CALL mpupdatei(grdtiltver,1)
  CALL mpupdatel(wrtsolo,1)
  CALL mpupdatel(wrtsoloqc,1)
  CALL mpupdatei(medfilt,1)
  CALL mpupdater(outtime,1)
  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)
  CALL mpupdatei(rad98opt,1)

  CALL mpupdatei(wtsoloopt,1)
  remap_on = .TRUE.
  IF (wtsoloopt > 0 .AND. wtsoloopt < 100) remap_on = .FALSE.

  IF (wrtsoloqc .AND. (.NOT. qc_on) ) THEN
    qc_on = .TRUE.
    WRITE(*,'(/,1x,a,L,a,/,10x,a,/)')                                   &
      'WARNING: wrtsoloqc = ',wrtsoloqc,', but qc_on is not on.',       &
      'So the program has turned on qc_on automatically.'
  END IF
!
!------------------------------------------------------------------------
!
! A few more initializations.
! Most of these variables will be reset with info obtained from the data file
!
!------------------------------------------------------------------------
!
  isource = 1
  istatus = 0
  icount  = 0
  gtspc_ref = 1000
  gtspc_vel = 250
  rfrst_ref = 1000
  rfrst_vel = 1000
  maxrgate = 1
  maxvgate = 1
  maxgate = 1
  maxazim = 1
  maxelev = 1
  itimfrst = 0
  refrngmin = rngmin
  refrngmax = rngmax
  kref=0
  kvel=0

  IF (myproc == ROOT) THEN

    CALL set_radar98_f(rad98opt)

    CALL get_88dinfo(radfname,radname,ivcp,itimfrst,                    &
                     radlat,radlon,radalt,istatus)
    WRITE(6,'(1x,a)')        'Back from get_88dinfo.'
    WRITE(6,'(1x,a,i6)')     'VCP: ',ivcp
    WRITE(6,'(1x,a,3f10.2)') 'radlat,radlon,radalt: ',radlat,radlon,radalt
    WRITE(6,'(1x,a,i6)')     'istatus: ',istatus
    IF( istatus /= 0) THEN
      errmessage = 'Bad status opening radar file'
      GOTO 8888
    END IF

    clrvcp = ( (ivcp == 31) .OR. (ivcp == 32) )
    IF( clrvcp .AND. clrvcpopt < 1 ) THEN
      istatus=-1
      errmessage = 'Clear VCP detected and no clear selected, exiting'
      GOTO 8888
    END IF

    IF(fntime) THEN
      lenfn=LEN_TRIM(radfname)
      iloc=index(radfname(1:lenfn),'/',back)
      jloc=iloc+index(radfname((iloc+1):lenfn),radname)
      IF(jloc > 0) THEN
        IF(radfname((jloc+12):(jloc+12)) == '_' ) THEN
          WRITE(6,'(a,a)') ' Getting time from filename substring: ', &
                         radfname((jloc+4):(jloc+18))
          READ(radfname((jloc+4):(jloc+18)),'(i4,2i2,1x,3i2)') &
                kmyear,kmmon,kmday,kmhr,kmmin,kmsec
          print *, ' Filename time variables:',kmyear,kmmon,kmday,kmhr,kmmin,kmsec
        ELSE
          WRITE(6,'(a,a)') ' Getting time from filename substring: ', &
                         radfname((jloc+4):(jloc+15))
          READ(radfname((jloc+4):(jloc+15)),'(i4,4i2)') &
                kmyear,kmmon,kmday,kmhr,kmmin
          kmsec=0
          print *, ' Filename time variables:',kmyear,kmmon,kmday,kmhr,kmmin,kmsec
        END IF
      ELSE
        kmyear=0
        kmmon=0
        kmday=0
        kmhr=0
        kmmin=0
        kmsec=0
      END IF
    END IF

    IF( clrvcp .AND. clrvcpopt < 2 ) THEN
      velproc = .FALSE.
      WRITE(6,'(1x,a)') 'Velocity processing not active for Clear VCP'
    END IF
    CALL get_88draddims(ivcp,maxrgate,maxvgate,maxazim,maxelev,istatus)
    WRITE(6,'(1x,a)')     'Back from 88draddims...'
    WRITE(6,'(1x,a,2i9)') 'maxrgate,maxvgate: ',maxrgate,maxvgate
    WRITE(6,'(1x,a,2i9)') 'maxazim,maxelev: ',maxazim,maxelev
    WRITE(6,'(1x,a,i6)')  'istatus: ',istatus
    IF (istatus /= 0) errmessage = 'return from get_88draddims'

  END IF

  8888 CONTINUE

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) THEN
    CALL arpsstop(TRIM(errmessage),istatus)
  END IF

  CALL mpupdatei(ivcp,1)
  CALL mpupdatei(itimfrst,1)
  CALL mpupdater(radlat,1)
  CALL mpupdater(radlon,1)
  CALL mpupdater(radalt,1)

  CALL mpupdatei(maxrgate,1)
  CALL mpupdatei(maxvgate,1)
  CALL mpupdatei(maxazim,1)
  CALL mpupdatei(maxelev,1)
  CALL mpupdatel(dualpout,1)
  CALL mpupdatel(clrvcp,1)
  CALL mpupdatel(velproc,1)
  CALL mpupdatel(dualpproc,1)

  IF (fntime) THEN
    CALL mpupdatei(kmyear,1)
    CALL mpupdatei(kmmon, 1)
    CALL mpupdatei(kmday, 1)
    CALL mpupdatei(kmhr,  1)
    CALL mpupdatei(kmmin, 1)
    CALL mpupdatei(kmsec, 1)
  END IF

  maxgate=max(maxrgate,maxvgate)

  CALL abss2ctim(itimfrst, iyear, imon, iday, ihr, imin, isec)
  WRITE(6,'(a,i4.4,5(a,i2.2))') ' itimfrst: ',iyear,'-',imon,'-',iday, &
          '_',ihr,':',imin,':',isec

  ALLOCATE(x(nx),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:x')
  ALLOCATE(y(ny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:y')
  ALLOCATE(z(nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:z')
  ALLOCATE(zp(nx,ny,nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:zp')

  ALLOCATE(xs(nx),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:xs')
  ALLOCATE(ys(ny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:ys')
  ALLOCATE(zps(nx,ny,nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:zps')

  ALLOCATE(tem2d(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:tem2d')

  ALLOCATE(j1(nx,ny,nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:j1')
  ALLOCATE(j2(nx,ny,nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:j2')
  ALLOCATE(j3(nx,ny,nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:j3')
  ALLOCATE(hterain(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:hterain')
  ALLOCATE(mapfct(nx,ny,8),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:mapfct')
  ALLOCATE(tem3d1(nx,ny,nz),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:tem3d1')

  IF( qc_on ) THEN
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
    CALL check_alloc_status(istat,'f88d2arps:u')
    ALLOCATE(v    (nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:v')
    ALLOCATE(ptprt(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:ptprt')
    ALLOCATE(pprt (nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:pprt')
    ALLOCATE(qv   (nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:qv')

    ALLOCATE(tem3d2(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:tem3d2')
    ALLOCATE(tem3d3(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:tem3d3')

    IF ((.NOT. traditional) .AND.             &
        (initopt == 3 .AND. (inifmt == 3 .OR. inifmt == 1)) ) THEN

      CALL get_gridxyzzp(nx,ny,nz,inigbf,inifmt,nprocx_in,nprocy_in,&
                           x,y,z,zp,istatus)

      CALL mapinit(nx,ny)

      IF (myproc == ROOT)  print *, ' Setting radar coordinate: '

      CALL radcoord(nx,ny,nz,x,y,z,zp,xs,ys,zps,                       &
                      radlat,radlon,radarx,radary)

      IF (myproc == ROOT) THEN
        print *, ' Radar x: ',(0.001*radarx),' km'
        print *, ' Radar y: ',(0.001*radary),' km'
        CALL flush(stdout)
      END IF

      CALL readuvt(nx,ny,nz,inifile,inifmt,nprocx_in,nprocy_in,        &
                   rtime,u,v,pprt,ptprt,qv,istatus)

      IF (myproc == ROOT) WRITE(6,'(1x,a,f10.2)')                      &
        'Environmental wind averaging radius: ',envavgr

      CALL extenvprf2(nx,ny,nz,nzsnd,x,y,zp,xs,ys,zps,                 &
             u,v,ptprt,pprt,qv,tem3d1,tem3d2,tem3d3,                   &
             radarx,radary,radalt,envavgr,rngmax,                      &
             zsnd,ktsnd,usnd,vsnd,rfrsnd,istatus);

      IF (myproc == ROOT) print *, ' Back from extenvprf'

      CALL mpupdatei(istatus, 1)

      IF(istatus < 0) THEN
        IF (myproc == ROOT) WRITE(6,'(a,i,/a)')  &
                        ' Status from extenvprf = ',istatus,' STOPPING'
        CALL arpsstop('arpsstop called from 88d2arps',1)
      END IF

    ELSE

      ALLOCATE(ubar(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:ubar')
      ALLOCATE(vbar(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:vbar')
      ALLOCATE(ptbar(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:ptbar')
      ALLOCATE(pbar(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:pbar')
      ALLOCATE(qvbar(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:qvbar')
      ALLOCATE(rhostr(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:rhostr')

      ALLOCATE(tem2dyz(ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem2dyz')
      ALLOCATE(tem2dxz(nx,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem2dxz')
      ALLOCATE(tem2dns(nx,ny,0:nstyps),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem2dns')

      ALLOCATE(trigs1(3*(nx-1)/2+1),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:trigs1')
      ALLOCATE(trigs2(3*(ny-1)/2+1),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:trigs2')
      ALLOCATE(ifax(13),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:ifax')

      ALLOCATE(wsave1(3*(ny-1)+15),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:wsave1')
      ALLOCATE(wsave2(3*(nx-1)+15),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:wsave2')
      ALLOCATE(vwork1(nx+1,ny+1),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:vwork1')
      ALLOCATE(vwork2(ny,nx+1),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:vwork2')

      ALLOCATE(qcumsrc(nx,ny,nz,5),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:qcsumsrc')
      ALLOCATE(prcrate(nx,ny,4),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:prcrate')
      ALLOCATE(exbcbuf(exbcbufsz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:exbcbuf')
      ALLOCATE(soiltyp(nx,ny,nstyps),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:soiltyp')
      ALLOCATE(stypfrct(nx,ny,nstyps),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:stypfrct')
      ALLOCATE(tem2dint(nx,ny),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem2dint')

      ALLOCATE(tem3d4(nx,ny,nz),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem3d4')

      ALLOCATE(tem3dsoil(nx,ny,nzsoil),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem3dsoil')
      ALLOCATE(tem4dsoilns(nx,ny,nzsoil,0:nstyps),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:tem4dsoilns')

      ALLOCATE(qscalar(nx,ny,nz,nscalar),stat=istat)
      CALL check_alloc_status(istat,'f88d2arps:qscalar')

      CALL set_lbcopt(1)

      CALL initgrdvar(nx,ny,nz,nzsoil,nt,nstyps,exbcbufsz,                &
                x,y,z,zp,tem3dsoil,hterain,mapfct,                        &
                j1,j2,j3,tem3dsoil,tem3d1,tem3d1,tem3d1,tem3d1,tem3dsoil, &
                u,v,tem3d1,tem3d1,ptprt,pprt,qv,qscalar,tem3d1,           &
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

      DEALLOCATE(trigs1)
      DEALLOCATE(trigs2)
      DEALLOCATE(ifax)
      DEALLOCATE(wsave1)
      DEALLOCATE(wsave2)
      DEALLOCATE(vwork1)
      DEALLOCATE(vwork2)
      DEALLOCATE(qcumsrc)
      DEALLOCATE(prcrate)
      DEALLOCATE(soiltyp)
      DEALLOCATE(stypfrct)
      DEALLOCATE(exbcbuf)
      DEALLOCATE(tem2dint)
      DEALLOCATE(tem2dyz)
      DEALLOCATE(tem2dxz)
      DEALLOCATE(tem2dns)
      DEALLOCATE(tem4dsoilns)
      DEALLOCATE(tem3dsoil)
      DEALLOCATE(qscalar)

      IF (myproc == ROOT) print *, ' Setting radar coordinate: '

      CALL radcoord(nx,ny,nz,x,y,z,zp,xs,ys,zps,                       &
                    radlat,radlon,radarx,radary)

      IF (myproc == ROOT) THEN
        print *, ' Radar x: ',(0.001*radarx),' km'
        print *, ' Radar y: ',(0.001*radary),' km'
      END IF

      IF (myproc == ROOT)  &
        print *, ' Environmental wind averaging radius: ',envavgr


      CALL extenvprf(nx,ny,nz,nzsnd,x,y,zp,xs,ys,zps,                  &
                 u,v,ptprt,pprt,qv,ptbar,pbar,tem3d1,tem3d2,tem3d3,    &
                 radarx,radary,radalt,envavgr,rngmax,                  &
                 zsnd,ktsnd,usnd,vsnd,rfrsnd,istatus)

      print *, ' Back from extenvprf'

      IF(istatus < 0) THEN
        IF (myproc == ROOT) WRITE(6,'(a,i,/a)')                        &
             ' Status from extenvprf = ',istatus,' STOPPING'
        CALL arpsstop('arpsstop called from 88d2arps',1)
      END IF

      DEALLOCATE(ubar, vbar, pbar, ptbar, qvbar, rhostr)
      DEALLOCATE(tem3d4)
    END IF

    DEALLOCATE(u, v, pprt, ptprt, qv)
    DEALLOCATE(tem3d2, tem3d3)

  ELSE

    ALLOCATE(zpsoil(nx,ny,nzsoil),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:zpsoil')
    ALLOCATE(j3soil(nx,ny,nzsoil),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:j3soiloil')
    ALLOCATE(j3soilinv(nx,ny,nzsoil),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:j3soilinv')

    ALLOCATE(tem1d1(nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:tem1d1')
    ALLOCATE(tem1d2(nz),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:tem1d2')

    IF (myproc == ROOT) print *, ' Setting ARPS grid coordinates ...'
    CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                        &
                hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,               &
                tem1d1,tem1d2,tem3d1)

    IF (myproc == ROOT) print *, ' Setting radar coordinate...'

    CALL radcoord(nx,ny,nz,x,y,z,zp,xs,ys,zps,                          &
                  radlat,radlon,radarx,radary)

    DEALLOCATE(zpsoil, j3soil, j3soilinv)
    DEALLOCATE(tem1d1, tem1d2)

  END IF

  DEALLOCATE(j1, j2, j3)
  DEALLOCATE(hterain)
  DEALLOCATE(mapfct)
  DEALLOCATE(tem3d1)

!-----------------------------------------------------------------------
!
! Allocation of analysis arrays
!
!-----------------------------------------------------------------------

  nxny=nx*ny
  ALLOCATE(icolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:icolp')
  ALLOCATE(jcolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:jcolp')
  ALLOCATE(xcolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:xcolp')
  ALLOCATE(ycolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:ycolp')
  ALLOCATE(zcolp(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:zcolp')
  ALLOCATE(havdat(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:havdat')

  ALLOCATE(kntbin(nsort),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:kntbin')

  ALLOCATE(colref(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:colref')
  IF (dualpproc) THEN
    ALLOCATE(colrhv(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:colrhv')
    ALLOCATE(colzdr(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:colzdr')
    ALLOCATE(colkdp(nz,nxny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:colkdp')
  END IF
  ALLOCATE(colvel(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:colvel')
  ALLOCATE(colnyq(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:colnyq')
  ALLOCATE(coltim(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:coltim')

!------------------------------------------------------------------------
!
! Radar data array allocations
!
!------------------------------------------------------------------------

  ALLOCATE(refdata(maxrgate,maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:refdata')
  ALLOCATE(refrange(maxrgate),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:refrange')
  ALLOCATE(veldata(maxvgate,maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:veldata')
  ALLOCATE(spwdata(maxvgate,maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:spwdata')
  ALLOCATE(rhvdata(maxrgate,maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:rhvdata')
  IF(dualpproc) THEN
    ALLOCATE(zdrdata(maxrgate,maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:zdrdata')
    ALLOCATE(phidata(maxrgate,maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:phidata')
    ALLOCATE(kdpdata(maxrgate,maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:kdpdata')
  END IF

  IF(qc_on) THEN
    ALLOCATE(unfvdata(maxvgate,maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:unfvdata')
    ALLOCATE(bkgvel(maxvgate,maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:bkgvel')
    ALLOCATE(bgate(maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:bgate')
    ALLOCATE(egate(maxazim),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:egate')
  END IF

  ALLOCATE(rtem(maxgate,maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:rtem')
  ALLOCATE(time(maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:time')
  ALLOCATE(gcfst(maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:gcfst')
  ALLOCATE(azim(maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:azim')
  ALLOCATE(vnyq(maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:vnyq')
  ALLOCATE(elev(maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:elev')
  ALLOCATE(istrgate(maxazim),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:istrgate')

  ALLOCATE(kntrgat(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:kntrgat')
  ALLOCATE(kntrazm(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:kntrazm')

  ALLOCATE(kntvgat(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:kntvgat')
  ALLOCATE(kntvazm(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:kntvazm')

  ALLOCATE(vnyqvol(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:vnyqvol')
  ALLOCATE(timevolr(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:timevolr')
  ALLOCATE(timevolv(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:timevolv')

  ALLOCATE(rngrvol(maxrgate,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:rngrvol')
  ALLOCATE(azmrvol(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:azmrvol')
  ALLOCATE(elvrvol(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:elvrvol')
  ALLOCATE(elvmnrvol(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:elvmnrvol')
  ALLOCATE(refvol(maxrgate,maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:refvol')
  IF( dualpproc ) THEN
    ALLOCATE(rhvvol(maxrgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:rhvvol')
    ALLOCATE(zdrvol(maxrgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:zdrvol')
    ALLOCATE(kdpvol(maxrgate,maxazim,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:kdpvol')
  END IF

  ALLOCATE(rngvvol(maxvgate,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:rngvvol')
  ALLOCATE(azmvvol(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:azmvvol')
  ALLOCATE(elvvvol(maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:elvvvol')
  ALLOCATE(elvmnvvol(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:elvmnvvol')
  ALLOCATE(velvol(maxvgate,maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:velvol')
  ALLOCATE(elvmean(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:elvmean')

  ALLOCATE(ngateref(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:ngateref')
  ALLOCATE(ngatevel(maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:ngatevel')

  ALLOCATE(rxvol(maxgate,maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:rxvol')
  ALLOCATE(ryvol(maxgate,maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:ryvol')
  ALLOCATE(rzvol(maxgate,maxazim,maxelev),stat=istat)
  CALL check_alloc_status(istat,'f88d2arps:rzvol')
!
!------------------------------------------------------------------------
!
!  Allocate arrays for grid-tilt writing for EnKF
!
!------------------------------------------------------------------------
!
  IF( wrtgrdtilt ) THEN
    ALLOCATE(grdtilthighref(nx,ny,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdtilthighref')
    ALLOCATE(grdrangeref(nx,ny,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdrangeref')
    ALLOCATE(grdslrref(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdslrref')
    ALLOCATE(grdazmref(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdazmref')
    ALLOCATE(grdtiltdata(nx,ny,maxelev,numvar),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdtiltdata')
    ALLOCATE(elevmeanref(maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:elevmeanref')

    ALLOCATE(tilttimeref(maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:tilttimeref')
    ALLOCATE(tilttimevel(maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:tilttimevel')

    ALLOCATE(grdtilthighvel(nx,ny,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdtilthighvel')
    ALLOCATE(grdrangevel(nx,ny,maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdrangevel')
    ALLOCATE(grdslrvel(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdslrvel')
    ALLOCATE(grdazmvel(nx,ny),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:grdazmvel')
    ALLOCATE(elevmeanvel(maxelev),stat=istat)
    CALL check_alloc_status(istat,'f88d2arps:elevmeanvel')
  END IF
!
!------------------------------------------------------------------------
!
!     Initialize 3-d radar data arrays
!
!------------------------------------------------------------------------
!
  print *, ' calling rmpinit'
  CALL rmpinit(nx,ny,nz,maxrgate,maxvgate,maxazim,maxelev,dualpproc,    &
          kntrgat,kntrazm,kntrelv,                                      &
          kntvgat,kntvazm,kntvelv,                                      &
          vnyqvol,timevolr,timevolv,                                    &
          rngrvol,azmrvol,elvrvol,elvmnrvol,refvol,rhvvol,zdrvol,kdpvol,&
          rngvvol,azmvvol,elvvvol,elvmnvvol,velvol,                     &
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
    iimon = imon
    iiday = iday
    iihr = ihr
    iimin = imin
    iisec = isec
  END IF

  IF (myproc == ROOT) THEN
    DO itilt=1,maxelev
      CALL rdtilt88d(maxrgate,maxvgate,maxazim,dualpproc,radfname,      &
                     refvarname,rhvvarname,zdrvarname,kdpvarname,       &
                     velvarname,spwvarname,itimfrst,                    &
                     ng_ref,rfrst_ref,gtspc_ref,                        &
                     ng_vel,rfrst_vel,gtspc_vel,                        &
                     nazim,azim,elev,time,vnyq,                         &
                     refdata,rhvdata,zdrdata,phidata,veldata,spwdata,   &
                     refstat,rhvstat,dlpstat,velstat,eofstat,istatus)
      IF(istatus /= 0) EXIT

!------------------------------------------------------------------------
!
!     Correct altitude, where necessary.
!     Correct true north, for example from TDWR data given in mag north.
!
!------------------------------------------------------------------------

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
        DO j=1, nazim
          kfn=0
          DO kkfn=nfixed_nyqv,1,-1
            IF( elev(j) <= elevlim_fnyqv(kkfn) ) kfn=kkfn
          END DO
          IF( kfn > 0 ) THEN
            IF(j==1)                                                   &
              WRITE(6,'(2(a,f6.2),a)') ' Resetting Nyquist at elev: ', &
                          elev(j),' to ',fixed_nyqv(kfn),' m/s.'
            vnyq(j)=fixed_nyqv(kfn)
          END IF
        END DO
      END IF

!------------------------------------------------------------------------
!
!     Process reflectivity data.
!
!------------------------------------------------------------------------

      print *, ' Back from rdtilt88d , istatus=',istatus
      print *, '    refstat= ',refstat,'  velstat= ',velstat
      print *, '    rhvstat= ',rhvstat,'  dlpstat= ',dlpstat
      print *, '    nazim: ',nazim,'  elev(1): ',elev(1)
      print *, '    ng_ref: ',ng_ref,'  gtspc_ref:',gtspc_ref
      print *, '    ng_vel: ',ng_vel,'  gtspc_vel:',gtspc_vel
      print *, '    eofstat= ',eofstat

      IF (dlpstat > 0) dualpdata = .TRUE.
      ngateref(itilt)=ng_ref
      ngatevel(itilt)=ng_vel

      IF( refstat == 1 ) THEN
        kref=kref+1
        icount = icount + 1
        WRITE(stdout,'(a,i4)')                                         &
        ' Processing base reflectivity data... ', icount

        WRITE(stdout,*)'Transferring ',nazim,' reflectivity radials.'

        WRITE(stdout,'(a)')                                            &
              '  Ref    j    dtime     azim     elev   refmin   refmax'
        DO j = 1, nazim, 20
          rmin=999.
          rmax=-999.
          DO i=1,ng_ref
            IF(refdata(i,j) > refchek) THEN
              rmin=min(refdata(i,j),rmin)
              rmax=max(refdata(i,j),rmax)
!             IF(refdata(i,j) > 40.0) THEN
!               WRITE(stdout,'(a,2i6,2f10.2)')' BIG-R i,j,azim,refl:',  &
!                     i,j,azim(j),refdata(i,j)
!             END IF
            END IF
          END DO
          WRITE(stdout,'(5x,i5,I9,4f9.1)') j,time(j),azim(j),elev(j),rmin,rmax
        END DO

        IF( qc_on) THEN
          write(6,'(1x,a)') 'Calling anomradial '
          CALL anomradial(maxrgate,maxazim,ng_ref,nazim,               &
                      rfrst_ref,gtspc_ref,refchek,anrflag,azim,refdata)

          IF( rhvstat == 1 ) THEN
            write(6,'(1x,a)') 'Calling RhoHV screen'
            CALL rhohvchk(maxrgate,maxvgate,maxazim,ng_ref,ng_vel,nazim,&
                 gtspc_ref,gtspc_vel,refchek,velchek,rhohvthr,         &
                 refdata,veldata,rhvdata)
          END IF

          write(6,'(1x,a)') 'Calling despekl '
          CALL despekl(maxrgate,maxazim,ng_ref,maxazim,refchek,medfilt,&
                 refdata,rtem)
        END IF
 
        IF( dlpstat == 1 ) THEN
          print *, ' rfrst_ref: ',rfrst_ref,'  gtscpc_ref: ',gtspc_ref
          DO igate=1,ng_ref
            refrange(igate)=0.001*(rfrst_ref+(igate-1)*gtspc_ref)
          END DO

          WRITE(stdout,'(a,i4)')                                       &
          ' Processing base dual-pol variables',icount

          CALL phi2kdp(maxrgate,maxazim,ng_ref,nazim,                  &
                       refrange,refdata,rhvdata,phidata,kdpdata,rtem)

          WRITE(stdout,'(a)')                                          &
              'Rho-HV   j     azim     elev   rhvmin   rhvmax'
          DO j = 1, nazim, 20
            rmin=999.
            rmax=-999.
            DO i=1,ng_ref
              IF(rhvdata(i,j) > rhvchek) THEN
                rmin=min(rhvdata(i,j),rmin)
                rmax=max(rhvdata(i,j),rmax)
              END IF
            END DO
            WRITE(stdout,'(5x,i5,2f9.1,2f9.4)')                        &
                 j,azim(j),elev(j),rmin,rmax
          END DO
          IF( qc_on) THEN
            write(6, *) ' calling despekl '
            CALL despekl(maxrgate,maxazim,ng_ref,maxazim,rhvchek,      &
                 dlpmedfilt,rhvdata,rtem)
          END IF

          WRITE(stdout,'(a)')                                          &
              ' Zdr     j     azim     elev   Zdrmin   Zdrmax'
          DO j = 1, nazim, 20
            rmin=999.
            rmax=-999.
            DO i=1,ng_ref
              IF(zdrdata(i,j) > zdrchek) THEN
                rmin=min(zdrdata(i,j),rmin)
                rmax=max(zdrdata(i,j),rmax)
              END IF
            END DO
            WRITE(stdout,'(5x,i5,2f9.1,2f9.4)')                        &
                 j,azim(j),elev(j),rmin,rmax
          END DO
          IF( qc_on) THEN
            write(6, *) ' calling despekl '
            CALL despekl(maxrgate,maxazim,ng_ref,maxazim,zdrchek,      &
                 dlpmedfilt,zdrdata,rtem)
          END IF

          WRITE(stdout,'(a)')                                          &
              '  Phi   j     azim     elev   phimin   phimax'
          DO j = 1, nazim, 20
            rmin=999.
            rmax=-999.
            DO i=1,ng_ref
              IF(phidata(i,j) > phichek) THEN
                rmin=min(phidata(i,j),rmin)
                rmax=max(phidata(i,j),rmax)
              END IF
            END DO
            WRITE(stdout,'(5x,i5,2f9.1,2f9.4)')                        &
                 j,azim(j),elev(j),rmin,rmax
          END DO

          WRITE(stdout,'(a)')                                          &
              '  Kdp   j     azim     elev   kdpmin   kdpmax'
          DO j = 1, nazim, 20
            rmin=999.
            rmax=-999.
            DO i=1,ng_ref
              IF(kdpdata(i,j) > kdpchek) THEN
                rmin=min(kdpdata(i,j),rmin)
                rmax=max(kdpdata(i,j),rmax)
              END IF
            END DO
            WRITE(stdout,'(5x,i5,2f9.1,2f9.4)')                        &
                 j,azim(j),elev(j),rmin,rmax
          END DO
          IF( qc_on) THEN
            write(6, *) ' calling despekl '
            CALL despekl(maxrgate,maxazim,ng_ref,maxazim,kdpchek, &
                 dlpmedfilt,kdpdata,rtem)
          END IF

          write(6, *) ' calling volbuilddlp for dual-pol data'
          write(6, *) '   elev=',elev(1),' nazim =',nazim,' ng_ref=',ng_ref
          timeset=1
          CALL volbuilddlp(maxrgate,maxazim,maxelev,ng_ref,nazim,      &
                  nyqset,timeset,                                      &
                  kntrgat,kntrazm,kntrelv,                             &
                  gtspc_ref,rfrst_ref,refchek,                         &
                  vnyq,time,                                           &
                  azim,elev,refdata,rhvdata,zdrdata,kdpdata,           &
                  vnyqvol,timevolr,                                    &
                  rngrvol,azmrvol,elvrvol,elvmnrvol,                   &
                  refvol,rhvvol,zdrvol,kdpvol)

        ELSE  ! no dual-pol data processing

          write(6, *) ' calling volbuild reflectivity'
          write(6, *) '   elev=',elev(1),' nazim =',nazim,' ng_ref=',ng_ref
          CALL volbuild(maxrgate,maxazim,maxelev,ng_ref,nazim,         &
                  nyqset,timeset,                                      &
                  kntrgat,kntrazm,kntrelv,                             &
                  gtspc_ref,rfrst_ref,refchek,                         &
                  vnyq,time,                                           &
                  azim,elev,refdata,                                   &
                  vnyqvol,timevolr,                                    &
                  rngrvol,azmrvol,elvrvol,elvmnrvol,refvol)

        END IF ! dlpstat
      END IF   ! refstat

!------------------------------------------------------------------------
!
!       Process velocity data.
!
!------------------------------------------------------------------------

      IF(velstat == 1) THEN
        kvel=kvel+1
        icount = icount + 1

        WRITE(stdout,'(a,i4)')' Processing base velocity data... ', icount

        WRITE(stdout,'(a,i6,a)')'Transferring',nazim,' velocity radials.'

        vnyquist=vnyq(1)
        WRITE(stdout,'(a,f10.2)') ' Nyquist velocity: ',vnyquist

        WRITE(stdout,'(a)')                                              &
            '  Vel     j     dtime    azim     elev     vr min   vr max'
        DO j = 1, nazim, 20
          rmin=999.
          rmax=-999.
          DO i=1,ng_vel
            IF(veldata(i,j) > velchek) THEN
              rmin=min(veldata(i,j),rmin)
              rmax=max(veldata(i,j),rmax)
            END IF
          END DO
          WRITE(stdout,'(2x,i5,I9,4f9.1)') j,time(j),azim(j),elev(j),rmin,rmax
        END DO

!------------------------------------------------------------------------
!
!       Call quality control and radar volume builder.
!
!------------------------------------------------------------------------

        nyqset=1
        timeset=1
        iangle=NINT(100.*elev(1))

        IF( qc_on ) THEN
          IF( unfdiag ) THEN
            varid='rdvela'
            WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
            CALL wrtvel(iangle,itilt,varid,                             &
                        iiyear,iimon,iiday,iihr,iimin,iisec,            &
                        gtspc_vel,rfrst_vel,vnyquist,                   &
                        radname,radlat,radlon,radalt,                   &
                        maxvgate,maxazim,ng_vel,nazim,                  &
                        azim,elev,veldata)
          END IF

          WRITE(6,'(a)') ' Calling despekl '
          CALL despekl(maxvgate,maxazim,ng_vel,maxazim,velchek,medfilt, &
                       veldata,rtem)

          IF( unfdiag ) THEN
            varid='rdvelb'
            WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
            CALL wrtvel(iangle,itilt,varid,                             &
                        iiyear,iimon,iiday,iihr,iimin,iisec,            &
                        gtspc_vel,rfrst_vel,vnyquist,                   &
                        radname,radlat,radlon,radalt,                   &
                        maxvgate,maxazim,ng_vel,nazim,                  &
                        azim,elev,veldata)
          END IF

          WRITE(6,'(a)') ' Calling unfoldnqc'
          CALL unfoldnqc(maxvgate,maxazim,ng_vel,nazim,                 &
                         nzsnd,zsnd,usnd,vsnd,rfrsnd,                   &
                         bkgopt,shropt,rfropt,                          &
                         gtspc_vel,rfrst_vel,iangle,itilt,              &
                         iiyear,iimon,iday,iihr,iimin,iisec,            &
                         radlat,radlon,radalt,                          &
                         veldata,spwdata,elev,azim,vnyq,                &
                         unfvdata,bkgvel,bgate,egate,rtem)

          IF( unfdiag ) THEN
            varid='rdvelu'
            WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
            CALL wrtvel(iangle,itilt,varid,                             &
                        iiyear,iimon,iiday,iihr,iimin,iisec,            &
                        gtspc_vel,rfrst_vel,vnyquist,                   &
                        radname,radlat,radlon,radalt,                   &
                        maxvgate,maxazim,ng_vel,nazim,                  &
                        azim,elev,unfvdata)
          END IF

          print *, ' calling volbuild velocity '
          print *, '    maxvgate=',maxvgate,'  ng_vel=',ng_vel
          print *, '    elv=',elev(1),' nazim =',nazim
          CALL volbuild(maxvgate,maxazim,maxelev,ng_vel,nazim,          &
                   nyqset,timeset,                                      &
                   kntvgat,kntvazm,kntvelv,                             &
                   gtspc_vel,rfrst_vel,velchek,                         &
                   vnyq,time,                                           &
                   azim,elev,unfvdata,                                  &
                   vnyqvol,timevolv,                                    &
                   rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)
        ELSE

          print *, ' calling volbuild velocity '
          print *, '    maxvgate=',maxvgate,'  ng_vel=',ng_vel
          print *, '    elv=',elev(1),' nazim =',nazim
          CALL volbuild(maxvgate,maxazim,maxelev,ng_vel,nazim,          &
                   nyqset,timeset,                                      &
                   kntvgat,kntvazm,kntvelv,                             &
                   gtspc_vel,rfrst_vel,velchek,                         &
                   vnyq,time,                                           &
                   azim,elev,veldata,                                   &
                   vnyqvol,timevolv,                                    &
                   rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)
        END IF

      END IF ! velstat

!------------------------------------------------------------------------
!
!      Write SOLO-II data
!
!------------------------------------------------------------------------
!
      IF( wrtsolo ) THEN
        WRITE(pname,'(a,I0)') 'SOLO, vcp=',ivcp
        plen = LEN_TRIM(pname)

        IF( refstat == 1 ) THEN
          iangle=NINT(100.*elev(1))
          ktime=itimfrst+time(1)
          CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
          CALL init_solo_f(pname,plen,radname,4,                        &
                   radlat,radlon,radalt,ivcp,                           &
                   kyear,kmon,kday,khr,kmin,ksec,                       &
                   iscan,kref,iangle,                                   &
                   maxrgate,maxvgate,maxazim,istatus)
          DO jazim=1,nazim
            ktime=itimfrst+time(jazim)
            CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
            CALL add_radial_to_sweep_f(kref,elev(jazim),azim(jazim),    &
                 1,0,ng_ref,ng_vel,ng_vel,vnyq,                         &
                 kyear,kmon,kday,khr,kmin,ksec,                         &
                 rfrst_ref,gtspc_ref,rfrst_vel,gtspc_vel,               &
                 refdata(1,jazim),veldata(1,jazim),spwdata(1,jazim),    &
                 istatus)
          END DO
          CALL write_sweep_f(istatus)
        END IF ! refstat

        IF( velstat == 1 ) THEN
          iangle=NINT(100.*elev(1))
          ktime=itimfrst+time(1)
          CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
          CALL init_solo_f(pname,plen,radname,4,                        &
                   radlat,radlon,radalt,ivcp,                           &
                   kyear,kmon,kday,khr,kmin,ksec,                       &
                   iscan,kvel,iangle,                                   &
                   maxrgate,maxvgate,maxazim,istatus)
          DO jazim=1,nazim
            ktime=itimfrst+time(jazim)
            CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
            CALL add_radial_to_sweep_f(kvel,elev(jazim),azim(jazim),    &
                 0,1,ng_ref,ng_vel,ng_vel,vnyq,                         &
                 kyear,kmon,kday,khr,kmin,ksec,                         &
                 rfrst_ref,gtspc_ref,rfrst_vel,gtspc_vel,               &
                 refdata(1,jazim),veldata(1,jazim),spwdata(1,jazim),    &
                 istatus)
          END DO
          CALL write_sweep_f(istatus)
        END IF ! refstat
      END IF ! wrtsolo

      IF( wrtsoloqc ) THEN
        WRITE(pname,'(a,I0)') 'SOLOQC, vcp=',ivcp
        plen = LEN_TRIM(pname)

        IF( refstat == 1 ) THEN
          iangle=NINT(100.*elev(1))
          ktime=itimfrst+time(1)
          CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
          CALL init_solo_f(pname,plen,radname,4,                        &
                   radlat,radlon,radalt,ivcp,                           &
                   kyear,kmon,kday,khr,kmin,ksec,                       &
                   iscan,kref,iangle,                                   &
                   maxrgate,maxvgate,maxazim,istatus)
          DO jazim=1,nazim
            ktime=itimfrst+time(jazim)
            CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
            CALL add_radial_to_sweep_f(kref,elev(jazim),azim(jazim),    &
                 1,0,ng_ref,ng_vel,ng_vel,vnyq,                         &
                 kyear,kmon,kday,khr,kmin,ksec,                         &
                 rfrst_ref,gtspc_ref,rfrst_vel,gtspc_vel,               &
                 refdata(1,jazim),veldata(1,jazim),spwdata(1,jazim),    &
                 istatus)
          END DO
          CALL write_sweep_f(istatus)
        END IF ! refstat

        IF( velstat == 1 ) THEN
          iangle=NINT(100.*elev(1))
          ktime=itimfrst+time(1)
          CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
          CALL init_solo_f(pname,plen,radname,4,                        &
                   radlat,radlon,radalt,ivcp,                           &
                   kyear,kmon,kday,khr,kmin,ksec,                       &
                   iscan,kvel,iangle,                                   &
                   maxrgate,maxvgate,maxazim,istatus)
          DO jazim=1,nazim
            ktime=itimfrst+time(jazim)
            CALL abss2ctim(ktime,kyear,kmon,kday,khr,kmin,ksec)
            CALL add_radial_to_sweep_f(kvel,elev(jazim),azim(jazim),    &
                 0,1,ng_ref,ng_vel,ng_vel,vnyq,                         &
                 kyear,kmon,kday,khr,kmin,ksec,                         &
                 rfrst_ref,gtspc_ref,rfrst_vel,gtspc_vel,               &
                 refdata(1,jazim),unfvdata(1,jazim),spwdata(1,jazim),   &
                 istatus)
          END DO
          CALL write_sweep_f(istatus)
        END IF ! refstat
      END IF ! wrtsoloqc

      IF(eofstat == 1) EXIT

    END DO  ! Tilt loop

    WRITE(6,'(a,i4,a,i4)') 'Exited tilt loop, istatus = ',istatus,      &
                           ',  eofstat = ',eofstat

    IF(velproc .OR. vad .AND. qc_on) THEN
      WRITE(6,'(a)') ' Calling quadunf for velocity'
      CALL quadunf(maxvgate,maxazim,maxelev,nzsnd,nsort,                &
                   velchek,velmedl,veldazl,iordunf,rfropt,              &
                   sortmin,dsort,                                       &
                   kntvgat,kntvazm,kntvelv,                             &
                   radarx,radary,radalt,dazim,                          &
                   rngmin,rngmax,zsnd,rfrsnd,                           &
                   rngvvol,azmvvol,elvvvol,elvmnvvol,vnyqvol,           &
                   kntbin,rxvol,ryvol,rzvol,velvol,                     &
                   istatus)
      IF( unfdiag ) THEN
        varid='rdvelz'
        DO kelv=1,kntvelv
          WRITE(6,'(3a,f7.1)')                                          &
           ' Calling wrtvel varid:',varid,' elev: ',elvmnvvol(kelv)
          iangle=100.*elvmnvvol(kelv)
          gtspc_vel=NINT(rngvvol(2,kelv)-rngvvol(1,kelv))
          rfrst_vel=NINT(rngvvol(1,kelv))
          print *, ' rfrst,gtspace: ',rfrst_vel,gtspc_vel
          ng_vel=kntvgat(1,kelv)
          DO jazim=1,kntvazm(kelv)
            ng_vel=max(ng_vel,kntvgat(jazim,kelv))
          END DO
          print *, ' n gates: ',ng_vel
          vnyquist=vnyqvol(1,kelv)
          CALL wrtvel(iangle,kelv,varid,                                &
                    iiyear,iimon,iiday,iihr,iimin,iisec,                &
                    gtspc_vel,rfrst_vel,vnyquist,                       &
                    radname,radlat,radlon,radalt,                       &
                    maxvgate,maxazim,ng_vel,kntvazm(kelv),              &
                    azmvvol(1,kelv),elvvvol(1,kelv),                    &
                    velvol(1,1,kelv))
        END DO
      END IF

    END IF

  END IF ! IF ROOT
!------------------------------------------------------------------------
!
! Create filename for output of remapped radar variables.
!
!------------------------------------------------------------------------

  IF (remap_on) THEN

    CALL mpupdatei(icount,1)

    IF (icount > 0) THEN

      IF (myproc == ROOT) THEN
        IF( qc_on ) THEN

          write(6, *) ' Calling apdetect ',gtspc_ref,gtspc_vel
          CALL apdetect(maxrgate,maxvgate,maxazim,maxelev,              &
                        kntrgat,kntrazm,kntrelv,                        &
                        kntvgat,kntvazm,kntvelv,                        &
                        refchek,velchek,                                &
                        gtspc_ref,gtspc_vel,                            &
                        winszrad,winszazim,ivcp,gclopt,gcvrlim,         &
                        rngrvol,azmrvol,elvrvol,                        &
                        rngvvol,azmvvol,elvvvol,                        &
                        refvol,velvol,rtem,                             &
                        istatus)
        END IF

      END IF

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
      ivar=1
      curtim = outtime

      CALL mpupdatei(kntrgat,maxazim*maxelev)
      CALL mpupdatei(kntrazm,maxelev)
      CALL mpupdatei(kntrelv,1)

      CALL mpupdater(rngrvol,maxrgate*maxelev)
      CALL mpupdater(azmrvol,maxazim*maxelev)
      CALL mpupdater(elvrvol,maxazim*maxelev)
      CALL mpupdater(refvol,maxrgate*maxazim*maxelev)
      CALL mpupdatei(timevolr,maxazim*maxelev)
      CALL mpupdater(vnyqvol,maxazim*maxelev)

      IF(velproc) THEN
        CALL mpupdatei(kntvgat,maxazim*maxelev)
        CALL mpupdatei(kntvazm,maxelev)
        CALL mpupdatei(kntvelv,1)

        CALL mpupdater(rngvvol,maxvgate*maxelev)
        CALL mpupdater(azmvvol,maxazim*maxelev)
        CALL mpupdater(elvvvol,maxazim*maxelev)
        CALL mpupdater(velvol,maxvgate*maxazim*maxelev)
        CALL mpupdatei(timevolv,maxazim*maxelev)
        CALL mpupdater(vnyqvol,maxazim*maxelev)
      END IF

      IF(dualpproc) THEN
        CALL mpupdater(rhvvol,maxrgate*maxazim*maxelev)
        CALL mpupdater(zdrvol,maxrgate*maxazim*maxelev)
        CALL mpupdater(kdpvol,maxrgate*maxazim*maxelev)
      END IF
!
!   VAD Wind Processing
!
      IF( (myproc == ROOT) .AND. vad ) THEN

        ALLOCATE(vadhgt(maxelev),stat=istat)
        CALL check_alloc_status(istat,'f88d2arps:vadhgt')
        ALLOCATE(vaddir(maxelev),stat=istat)
        CALL check_alloc_status(istat,'f88d2arps:vaddir')
        ALLOCATE(vadspd(maxelev),stat=istat)
        CALL check_alloc_status(istat,'f88d2arps:vadspd')

        CALL mkvadfnm(dirname,radname,                                 &
                      iiyear,iimon,iiday,iihr,iimin,iisec,             &
                      vad_fname)
        WRITE(stdout,'(a,a)')                                          &
              ' Filename for VAD output: ',TRIM(vad_fname)
        CALL flush(stdout)

        CALL vadvol(maxvgate,maxazim,maxelev,                          &
                    radalt,velchek,vadradius,vadwidth,vadminknt,       &
                    kntvgat,kntvazm,kntvelv,                           &
                    rngvvol,azmvvol,elvvvol,velvol,                    &
                    vadhgt,vaddir,vadspd)
        CALL wtvadprf(maxelev,vad_fname,radname,vadsrc,                &
                      radlat,radlon,radalt,vadhgt,vaddir,vadspd);

        DEALLOCATE(vadhgt,vaddir,vadspd)

      END IF
!
!     Set-up arrays needed for remapping
!
      IF (myproc == ROOT) WRITE(stdout,'(a)') ' Calling rmpsetup...'
      CALL rmpsetup(maxrgate,maxvgate,maxgate,maxazim,maxelev,         &
                      nx,ny,nxny,nz,nzsnd,                             &
                      rfropt,refchek,velchek,bmwidth,velproc,          &
                      kntrgat,kntrazm,kntrelv,                         &
                      kntvgat,kntvazm,kntvelv,                         &
                      radlat,radlon,radarx,radary,radalt,              &
                      dazim,rngmin,rngmax,                             &
                      rngrvol,azmrvol,elvrvol,                         &
                      rngvvol,azmvvol,elvvvol,                         &
                      refvol,velvol,rxvol,ryvol,rzvol,                 &
                      xs,ys,zps,zsnd,rfrsnd,ncolp,ncoltot,             &
                      havdat,icolp,jcolp,xcolp,ycolp,zcolp,istatus)
      IF (myproc == ROOT) &
        WRITE(stdout,'(a,i12)') ' Back from rmpsetup, ncoltot=',ncoltot
      WRITE(stdout,'(a,i6,a,i12)')' myproc: ',myproc,'  ncolp: ',ncolp

      DEALLOCATE(havdat)

      IF (nprocs > 1 ) THEN

        CALL distrallcol(nxny,nz,ncolp,ncoltot,                        &
                         icolp,jcolp,xcolp,ycolp,zcolp,lvldbg,istatus)

      END IF
!
!  Call remap for reflectivity.
!  Note here that the rxvol,ryvol and rzvol arrays have previously been
!  set for the reflectivity gates.
!
      IF (myproc == ROOT) THEN
        WRITE(stdout,'(a)') ' Calling remap3dcol for reflectivity '
        WRITE(stdout,'(a,i6)') ' Reflectivity tilt count:',kntrelv
      END IF
      varfill=fillref
      CALL remap3dcol(maxrgate,maxgate,maxazim,maxelev,nxny,nz,        &
                    nzsnd,nsort,ncolp,                                 &
                    varfill,ivar,nyqset,timeset,rfropt,                &
                    refchek,refmiss,bmwidth,refmedl,refdazl,iordref,   &
                    sortmin,dsort,                                     &
                    kntrgat,kntrazm,kntrelv,                           &
                    radlat,radlon,radarx,radary,radalt,dazim,          &
                    rngmin,rngmax,                                     &
                    rngrvol,azmrvol,elvrvol,                           &
                    refvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,         &
                    zsnd,rfrsnd,kntbin,elvmean,                        &
                    xcolp,ycolp,zcolp,                                 &
                    colref,coltim,colnyq,istatus)

      DO kelv=1,kntrelv
        refelvmin=min(elvmean(kelv),refelvmin)
        refelvmax=max(elvmean(kelv),refelvmax)
      END DO

      IF(myproc == ROOT) print *, ' outtime: ',outtime

      IF( ref3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,                 &
                    colref,icolp,jcolp,                                &
                    refid,refname,refunits,refmiss,outtime,runname,    &
                    dirname,dmpfmt,hdfcompr,istatus)

      IF( dualpdata .AND. dualpproc ) THEN
        IF( rhv3d .OR. dualpout > 0 ) THEN
          IF (myproc == ROOT)                                          &
            WRITE(stdout,'(/a)') ' Calling remap3dcol for Rho-HV'
          ivar=4
          varfill=0
          CALL remap3dcol(maxrgate,maxgate,maxazim,maxelev,nxny,nz,    &
                    nzsnd,nsort,ncolp,                                 &
                    varfill,ivar,nyqset,timeset,rfropt,                &
                    rhvchek,rhvmiss,bmwidth,rhvmedl,refdazl,iordref,   &
                    sortmin,dsort,                                     &
                    kntrgat,kntrazm,kntrelv,                           &
                    radlat,radlon,radarx,radary,radalt,dazim,          &
                    rngmin,rngmax,                                     &
                    rngrvol,azmrvol,elvrvol,                           &
                    rhvvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,         &
                    zsnd,rfrsnd,kntbin,elvmean,                        &
                    xcolp,ycolp,zcolp,                                 &
                    colrhv,coltim,colnyq,istatus)

          IF ( rhv3d) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,             &
                        colrhv,icolp,jcolp,                            &
                        rhvid,rhvname,rhvunits,rhvmiss,outtime,runname,&
                        dirname,dmpfmt,hdfcompr,istatus)
        END IF
        IF( zdr3d .OR. dualpout > 0 ) THEN
          IF (myproc == ROOT)                                          &
            WRITE(stdout,'(/a)') ' Calling remap3dcol for Zdr'
          ivar=5
          varfill=0
          CALL remap3dcol(maxrgate,maxgate,maxazim,maxelev,nxny,nz,    &
                    nzsnd,nsort,ncolp,                                 &
                    varfill,ivar,nyqset,timeset,rfropt,                &
                    zdrchek,zdrmiss,bmwidth,zdrmedl,refdazl,iordref,   &
                    sortmin,dsort,                                     &
                    kntrgat,kntrazm,kntrelv,                           &
                    radlat,radlon,radarx,radary,radalt,dazim,          &
                    rngmin,rngmax,                                     &
                    rngrvol,azmrvol,elvrvol,                           &
                    zdrvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,         &
                    zsnd,rfrsnd,kntbin,elvmean,                        &
                    xcolp,ycolp,zcolp,                                 &
                    colzdr,coltim,colnyq,istatus)

          IF ( zdr3d) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,             &
                        colzdr,icolp,jcolp,                            &
                        zdrid,zdrname,zdrunits,zdrmiss,outtime,runname,&
                        dirname,dmpfmt,hdfcompr,istatus)
        END IF
        IF( kdp3d .OR. dualpout > 0 ) THEN
          IF (myproc == ROOT)                                          &
            WRITE(stdout,'(/a)') ' Calling remap3dcol for Kdp'
          ivar=7
          varfill=0
          CALL remap3dcol(maxrgate,maxgate,maxazim,maxelev,nxny,nz,    &
                    nzsnd,nsort,ncolp,                                 &
                    varfill,ivar,nyqset,timeset,rfropt,                &
                    kdpchek,kdpmiss,bmwidth,kdpmedl,refdazl,iordref,   &
                    sortmin,dsort,                                     &
                    kntrgat,kntrazm,kntrelv,                           &
                    radlat,radlon,radarx,radary,radalt,dazim,          &
                    rngmin,rngmax,                                     &
                    rngrvol,azmrvol,elvrvol,                           &
                    kdpvol,timevolr,vnyqvol,rxvol,ryvol,rzvol,         &
                    zsnd,rfrsnd,kntbin,elvmean,                        &
                    xcolp,ycolp,zcolp,                                 &
                    colkdp,coltim,colnyq,istatus)

          IF ( kdp3d) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,             &
                        colkdp,icolp,jcolp,                            &
                        kdpid,kdpname,kdpunits,kdpmiss,outtime,runname,&
                        dirname,dmpfmt,hdfcompr,istatus)
        END IF
      END IF  ! dualpproc

      IF( wrtgrdtilt ) THEN
        IF (myproc == ROOT) WRITE(stdout,'(/a)')                       &
                            ' Calling remap2dcts for reflectivity'
        CALL flush(stdout)
        CALL remap2dcts(maxrgate,maxazim,maxelev,nx,ny,nz,nzsnd,       &
                 kntrazm,kntrelv,                                      &
                 radlat,radlon,radarx,radary,radalt,dazim,             &
                 rngmin,rngmax,itimfrst,                               &
                 rngrvol,azmrvol,elvrvol,                              &
                 refvol,timevolr,                                      &
                 xs,ys,zps,zsnd,rfrsnd,                                &
                 grdtilthighref,grdrangeref,grdslrref,grdazmref,       &
                 grdtiltdata(:,:,:,2),tilttimeref,elevmeanref,ngateref,istatus)
        IF (myproc == ROOT) WRITE(stdout,'(a)') ' Back from remap2dts'
        CALL flush(stdout)
      END IF

      IF( ref2d ) THEN
        vardump = 1
        ivar    = 1
        varid   = 'refl2d'
        varname = 'Low-level reflect'
        IF (myproc == ROOT) THEN
          WRITE(stdout,'(/a)') ' Calling remap2d for reflectivity '
        END IF
        CALL remap2d(maxrgate,maxazim,maxelev,nx,ny,nzsnd,nsort,       &
                     vardump,ivar,rfropt,varid,varname,refunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     refchek,refmiss,refmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     refvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)
      END IF

      IF ( zdr2d ) THEN
        IF( dualpdata ) THEN
          vardump = 1
          ivar    = 1
          varid   = 'zdr_2d'
          varname = 'Low-level Zdr'
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(/a)') ' Calling remap2d for Zdr '
          END IF
          CALL remap2d(maxrgate,maxazim,maxelev,nx,ny,nzsnd,nsort,     &
                     vardump,ivar,rfropt,varid,varname,zdrunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     zdrchek,zdrmiss,zdrmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     zdrvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)
        ELSE
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(1x,a)')                                     &
                  'Dual-pol data unavailable, skipping 2D Zdr'
          END IF
        END IF
      END IF

      IF( kdp2d ) THEN
        IF( dualpdata ) THEN
          vardump = 1
          ivar    = 1
          varid   = 'kdp_2d'
          varname = 'Low-level Kdp'
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(/a)') ' Calling remap2d for Kdp '
          END IF
          CALL remap2d(maxrgate,maxazim,maxelev,nx,ny,nzsnd,nsort,     &
                     vardump,ivar,rfropt,varid,varname,refunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     refchek,refmiss,refmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     kdpvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)
        ELSE
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(1x,a)')                                     &
                  'Dual-pol data unavailable, skipping 2D Kdp'
          END IF
        END IF
      END IF

      IF( kdp2d ) THEN
        IF( dualpdata ) THEN
          vardump = 1
          ivar    = 1
          varid   = 'rhv_2d'
          varname = 'Low-level correlation'
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(/a)') ' Calling remap2d for Rho-HV '
          END IF
          CALL remap2d(maxrgate,maxazim,maxelev,nx,ny,nzsnd,nsort,     &
                     vardump,ivar,rfropt,varid,varname,refunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     rhvchek,rhvmiss,rhvmedl,refdazl,iordref,          &
                     sortmin,dsort,                                    &
                     kntrgat,kntrazm,kntrelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngrvol,azmrvol,elvrvol,            &
                     rhvvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)
        ELSE
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(1x,a)')                                     &
                  'Dual-pol data unavailable, skipping 2D Rho-HV'
          END IF
        END IF
      END IF

      IF( velproc ) THEN
        nyqset  = 1
        timeset = 1
        vardump = 0
        IF(vel3d) vardump = 1
        varfill = 0
        ivar    = 2

        CALL rgatexyz(maxvgate,maxgate,maxazim,maxelev,nzsnd,          &
                      rfropt,lvldbg,                                   &
                      kntvgat,kntvazm,kntvelv,                         &
                      radlat,radlon,radarx,radary,radalt,              &
                      rngvvol,azmvvol,elvvvol,                         &
                      zsnd,rfrsnd,                                     &
                      rxvol,ryvol,rzvol,istatus)

        IF (myproc == ROOT) THEN
          WRITE(stdout,'(/a)') ' Calling remap3dcol for velocity '
          WRITE(stdout,'(a,i6)') ' Velocity tilt count:',kntvelv
        END IF
        CALL remap3dcol(maxvgate,maxgate,maxazim,maxelev,nxny,nz,      &
                    nzsnd,nsort,ncolp,                                 &
                    varfill,ivar,nyqset,timeset,rfropt,                &
                    velchek,velmiss,bmwidth,velmedl,veldazl,iordvel,   &
                    sortmin,dsort,                                     &
                    kntvgat,kntvazm,kntvelv,                           &
                    radlat,radlon,radarx,radary,radalt,dazim,          &
                    rngmin,rngmax,                                     &
                    rngvvol,azmvvol,elvvvol,                           &
                    velvol,timevolv,vnyqvol,rxvol,ryvol,rzvol,         &
                    zsnd,rfrsnd,kntbin,elvmean,                        &
                    xcolp,ycolp,zcolp,                                 &
                    colvel,coltim,colnyq,istatus)

        IF( vel3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,               &
                        colvel,icolp,jcolp,                            &
                        velid,velname,velunits,velmiss,outtime,runname,&
                        dirname,dmpfmt,hdfcompr,istatus)

      END IF

      IF( wrtgrdtilt ) THEN

        IF (myproc == ROOT) WRITE(stdout,'(/a)')                       &
                            ' Calling remap2dcts for velocity'
        CALL flush(stdout)
        CALL remap2dcts(maxvgate,maxazim,maxelev,nx,ny,nz,nzsnd,       &
               kntvazm,kntvelv,                                        &
               radlat,radlon,radarx,radary,radalt,dazim,               &
               rngmin,rngmax,itimfrst,                                 &
               rngvvol,azmvvol,elvvvol,                                &
               velvol,timevolv,                                        &
               xs,ys,zps,zsnd,rfrsnd,                                  &
               grdtilthighvel,grdrangevel,grdslrvel,grdazmvel,         &
               grdtiltdata(:,:,:,1),tilttimevel,elevmeanvel,ngatevel,istatus)
        IF (myproc == ROOT) WRITE(stdout,'(a)') ' Back from remap2dts'
        CALL flush(stdout)
      END IF

      IF( vel2d ) THEN
        vardump  = 1
        ivar     = 2
        varid    = 'radv2d'
        varname  = 'Low-level Velocity'
        varunits = 'm/s'
        IF (myproc == ROOT) THEN
          WRITE(stdout,'(a)') ' Calling remap2d for velocity '
        END IF
        CALL remap2d(maxvgate,maxazim,maxelev,nx,ny,nzsnd,nsort,       &
                     vardump,ivar,rfropt,varid,varname,varunits,       &
                     dmpfmt,hdf4cmpr,                                  &
                     velchek,velmiss,velmedl,veldazl,iordvel,          &
                     sortmin,dsort,                                    &
                     kntvgat,kntvazm,kntvelv,                          &
                     radlat,radlon,radarx,radary,radalt,dazim,         &
                     rngmin,rngmax,rngvvol,azmvvol,elvvvol,            &
                     velvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,      &
                     tem2d,istatus)
      END IF

!------------------------------------------------------------------------
!
!   Create radar column data filename and write file.
!
!------------------------------------------------------------------------

      full_fname = ' '
      CALL mkradfnm(dmpfmt,dirname,ldirnam,radname,iiyear,iimon,iiday, &
                    iihr, iimin, iisec, full_fname, len_fname)
      IF (myproc == ROOT) &
        WRITE(stdout,'(1x,2a)') ' Remapped filename for this volume: ',&
          TRIM(full_fname)

      IF(.NOT. wrtgrdtilt) THEN

        WRITE(*,'(1x,a,I3,3(a,I8))') 'myproc:',myproc,                 &
                 ' ncolp:',ncolp,', ncoltot:',ncoltot,', nxny: ',nxny
        dualpol = dualpout
        IF( .NOT. dualpdata) dualpol = 0
        IF( radband == 2 .AND. radname(1:1) == 'T' ) isource = 4    ! 5-cm TDWR
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
        CALL wrtradcol(nx,ny,nxny,nz,ncolp,ncoltot,                    &
                   dmpfmt,iradfmt,hdf4cmpr,dmpzero,dualpol,            &
                   full_fname,radname,radlat,radlon,radalt,            &
                   iiyear,iimon,iiday,iihr,iimin,iisec,ivcp,isource,   &
                   refelvmin,refelvmax,refrngmin,refrngmax,            &
                   icolp,jcolp,xcolp,ycolp,zcolp,                      &
                   colref,colvel,colnyq,coltim,colrhv,colzdr,colkdp,   &
                   istatus)
      END IF

!------------------------------------------------------------------------
!
!   Write tilt data
!
!------------------------------------------------------------------------

      IF( myproc == ROOT .AND. wrttilt ) THEN
         CALL wrttilts(full_fname,maxvgate,maxrgate,maxazim,maxelev,   &
                   rngvvol,azmvvol,elvvvol,velvol,                     &
                   rngrvol,azmrvol,elvrvol,refvol,                     &
                   radalt,radlat,radlon,                               &
                   iiyear,imon,iiday,iihr,iimin,iisec )
      END IF

    END IF

    IF( wrtgrdtilt ) THEN

       WRITE(full_fname,'(a,a4,a1,i4.4,2(i2.2),a1,3(i2.2))')            &
             dirname(1:ldirnam),radname,'.',iiyear,iimon,iiday,'.',     &
             iihr, iimin, iisec

       IF (myproc == ROOT)                                              &
          WRITE(stdout,'(1x,2a)')                                       &
            'Grid-tilt data filename for this volume: ',TRIM(full_fname)

       CALL wrtgridtilt(full_fname,fntimopt,                            &
              maxelev,nx,ny,nz,kntrelv,radname,                         &
              radlat,radlon,radarx,radary,radalt,dazim,rngmin,rngmax,   &
              grdtilthighref,grdrangeref,grdslrref,grdazmref,           &
              grdtiltdata,tilttimeref,tilttimevel,itimfrst,elevmeanref, &
              grdtilthighvel,grdrangevel,grdslrvel,grdazmvel,           &
              elevmeanvel,xs,ys,zps,ivcp,isource,grdtiltver,numvar)

    END IF
  END IF  ! remap_on

!------------------------------------------------------------------------
!
! The End.
!
!------------------------------------------------------------------------

  IF (myproc == ROOT)  &
    WRITE(stdout,'(/a/)') '  === Normal termination of 88D2ARPS === '

  CALL mpexit(0)

  STOP
END PROGRAM f88d2arps
