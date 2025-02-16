!########################################################################
!########################################################################
!#########                                                      #########
!#########                   PROGRAM nids2arps                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

PROGRAM nids2arps

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads WSR 88D radar data from NIDS Level-III raw data files (e.g., provided
! by LDM feed) and remaps the data onto the ARPS Cartesian grid.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! MODIFICATIONS:
!
! Yunheng Wang  (July 2001)
! Converted from nids2arps.c and using NWS public code instead of
! wsibase.c. Based on the similar work on Echo Top and VIL by
! Eric Kemp 2001
!
! Eric Kemp, 9 August 2001
! Added VIL, Echo Top, and DPA product remapping -- these data are
! currently written in ARPSPLT "arbitrary variable" format to allow for
! plotting.  In the future, a better file format will be developed and
! added for used with ADAS.  Also, changed code to use INTEGER*4 and
! similar variable types instead of using the KIND statement.  Finally,
! consolidated parameters into module N2ACONST, and cleaned up code.
!
! Eric Kemp, 17 August 2001
! Modified code to avoid use of INTEGER*4, INTEGER*1, etc.  Based on
! work by Yunheng Wang.
!
! Eric Kemp, 10 September 2001
! Modified driver code to avoid use of unit 11.  This fixes conflict
! with subroutine 'setlookup', which uses unit 11 for reading/writing the
! radar map parameter file.
!
! Eric Kemp, 17 September 2001
! Consolidated code, allowing for all NIDS data files to be listed
! in the nids2arps.input file in the new 'nidsfn' variable array.  Also,
! added an option to adjust the remapped 3D reflectivity field using
! the NIDS Echo Top and VIL products.
!
! Keith Brewster, 8 September 2008
! Replaced the varsort array with kntbin for new method of median
! calculation.
!
! Keith Brewster, 17 March 2010
! Updated code for new input file used for all remapping codes.
! Updated code for new balanced MPI strategy.
!
! Y. Wang (05/08/2012)
! Replace the non-standard calls of getarg with the Fortran 2003
! standard call of GET_COMMAND_ARGUMENT.
!
! Keith Brewster (03/15/2013)
! Updated code for new dual pol variables -- not yet inplemented for NIDS,
! but placeholders needed for certain subroutines shared with other
! radar data remapping programs.
!
!------------------------------------------------------------------------
!
! REMARKS:
!
!------------------------------------------------------------------------
!
  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Include files
!
!------------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'nidscst.inc'
  INCLUDE 'remapcst.inc'
  INCLUDE 'mp.inc'

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER :: nx, ny, nz, nzsoil, itime, nstyps
  REAL    :: radalt, radlat, radlon
  REAL    :: radarx, radary

!------------------------------------------------------------------------
!
! NWS decoder variables
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: itrailerdim = 20000
  INTEGER, PARAMETER :: iheaderdim = 160

  INTEGER :: prodID                         ! product ID code
  INTEGER :: prodtype                       ! product type code

  INTEGER :: ival0, ival1

  INTEGER :: ival(maxlen)
  INTEGER :: iheader(iheaderdim)
  INTEGER :: ifield(maxgate,maxazim)
  INTEGER :: ifield_raster(ndim1,ndim2)
  INTEGER :: ifield_dpa(ndim1,ndim2)

  CHARACTER(LEN=1)  :: itrailer(itrailerdim)

  INTEGER :: msglen,icode,isite,isite_lat,isite_lon,isite_elev,          &
             iyear,imonth,iday,ihr,imin,isec,iprod,ngate,nazim,          &
             maxval2,ihr1,imin1,ihr2,imin2,istm_spd,istm_dir,            &
             len_header,len_trailer
  INTEGER :: iazmuth(maxazim)

  INTEGER :: maxval1, ivcp, ielev
  INTEGER :: iyr1,imonth1,iday1,iyr2,imonth2,iday2

  INTEGER :: icatt(0:15)
  INTEGER :: icats(0:15)
  REAL    :: rcats(0:15)

  INTEGER :: ndims_p,nrowsp,ncolsp

!------------------------------------------------------------------------
!
! Variables for subroutine remap3dcol
!
!------------------------------------------------------------------------
 
  INTEGER, PARAMETER :: rfropt = 1  ! 4/3rds earth refraction option
  INTEGER, PARAMETER :: nsort = 601
  REAL,    PARAMETER :: sortmin=-150.0
  REAL,    PARAMETER :: sortmax= 150.0
  INTEGER, PARAMETER :: nzsnd = 201
  REAL,    PARAMETER :: dzsnd=100.
  REAL :: zsnd(nzsnd)             ! hgt levels for refractivity sounding
  REAL :: rfrsnd(nzsnd)           ! refractivity sounding
  REAL,    PARAMETER :: rfrconst = (-1./(4.*6371.E03))   ! -1/(4e)

  REAL,    ALLOCATABLE :: rdata(:,:) ! refl (dBZ) or velocity data (m/s)
  REAL,    ALLOCATABLE :: rtem(:,:)  ! temporary array for radar data processing
  INTEGER, ALLOCATABLE :: time(:)    ! time offset from itimfrst
  REAL,    ALLOCATABLE :: azim(:)    ! azimuth angle for each radial(degree)
  REAL,    ALLOCATABLE :: elev(:)    ! elevation angle for each radial (degree)
  REAL,    ALLOCATABLE :: vnyq(:)    ! Nyquist velocity for this tilt
  INTEGER :: itimfrst             ! time of first radial in volume
  INTEGER :: rfirstg              ! range to first gate (meters)
  INTEGER :: gatesp               ! gate spacing (meters)
  REAL    :: eleva                ! elevation angle for this tilt

  INTEGER, ALLOCATABLE :: kntrgat(:,:)
  INTEGER, ALLOCATABLE :: kntrazm(:)
  INTEGER :: kntrelv

  INTEGER, ALLOCATABLE :: kntvgat(:,:)
  INTEGER, ALLOCATABLE :: kntvazm(:)
  INTEGER :: kntvelv

  INTEGER, ALLOCATABLE :: timevolr(:,:)
  INTEGER, ALLOCATABLE :: timevolv(:,:)
  REAL,    ALLOCATABLE :: nyqvvol(:,:)

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

  REAL, ALLOCATABLE :: x(:)
  REAL, ALLOCATABLE :: y(:)
  REAL, ALLOCATABLE :: z(:)
  REAL, ALLOCATABLE :: zp(:,:,:)
  REAL, ALLOCATABLE :: zpsoil(:,:,:)
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
  REAL, ALLOCATABLE :: j3soil(:,:,:)
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)
  REAL, ALLOCATABLE :: tem1d1(:)
  REAL, ALLOCATABLE :: tem1d2(:)
!  REAL, ALLOCATABLE :: tem1d3(:)
!  REAL, ALLOCATABLE :: tem1d4(:)
!  REAL, ALLOCATABLE :: tem1d5(:)
!  REAL, ALLOCATABLE :: tem1d6(:)
  REAL, ALLOCATABLE :: tem2d(:,:)
  REAL, ALLOCATABLE :: tem3d(:,:,:)

  REAL, ALLOCATABLE :: grdvel2(:,:,:)
  REAL, ALLOCATABLE :: grdref2(:,:,:)
  REAL, ALLOCATABLE :: wgtvel(:,:,:)
  REAL, ALLOCATABLE :: wgtref(:,:,:)
  REAL, ALLOCATABLE :: wpotvel(:,:,:)
  REAL, ALLOCATABLE :: wpotref(:,:,:)
  REAL, ALLOCATABLE :: gridazm(:,:)
  REAL, ALLOCATABLE :: gridrng(:,:)

!------------------------------------------------------------------------
!
! Arrays for raster (Echo Top and VIL) and DPA data.
!
!------------------------------------------------------------------------

  REAL, ALLOCATABLE :: rasterdata(:,:)
  REAL, ALLOCATABLE :: rasterx(:,:)
  REAL, ALLOCATABLE :: rastery(:,:)
  REAL, ALLOCATABLE :: rasterlat(:,:)
  REAL, ALLOCATABLE :: rasterlon(:,:)
  REAL, ALLOCATABLE :: remapped_raster(:,:)

  REAL, ALLOCATABLE :: dpadata(:,:)
  REAL, ALLOCATABLE :: dpax(:,:)
  REAL, ALLOCATABLE :: dpay(:,:)
  REAL, ALLOCATABLE :: dpalat(:,:)
  REAL, ALLOCATABLE :: dpalon(:,:)
  REAL, ALLOCATABLE :: remapped_dpa(:,:)

  REAL, ALLOCATABLE :: arpslat(:,:)
  REAL, ALLOCATABLE :: arpslon(:,:)
  REAL, ALLOCATABLE :: arpsrasterx(:,:)
  REAL, ALLOCATABLE :: arpsrastery(:,:)
  REAL, ALLOCATABLE :: arpsdpax(:,:)
  REAL, ALLOCATABLE :: arpsdpay(:,:)

!------------------------------------------------------------------------
!
! New variables for adjusting reflectivity (6 September 2001)
!
!------------------------------------------------------------------------

  REAL, ALLOCATABLE :: nidset(:,:)
  REAL, ALLOCATABLE :: nidsvil(:,:)

!------------------------------------------------------------------------
!
! Misc. variables
!
!------------------------------------------------------------------------

  CHARACTER(LEN=1) :: cval(maxlen)

  CHARACTER(LEN=256) :: infile

  INTEGER :: n, i, j, k, kelv, ios, iskip, istat, istatus
  INTEGER :: ilist,iopstat,irdstat
  INTEGER :: irefgatsp,ivelgatsp
  INTEGER :: unitin

!  INTEGER :: lenmpflstr

  CHARACTER(LEN=256) :: full_fname
  INTEGER :: len_dir, len_fname

  LOGICAL :: velproc,dualpproc
  LOGICAL :: fntime
  LOGICAL :: ncdchdr
  LOGICAL :: vel2d,vel3d
  LOGICAL :: ref2d,ref3d
  LOGICAL :: noclrv
  LOGICAL :: qc_on
  LOGICAL :: traditional
  LOGICAL :: vad
  LOGICAL :: wrttilt

  LOGICAL :: i_first_scan

  INTEGER :: iiyear,iimonth,iiday,iihr,iimin,iisec
  INTEGER :: isource, icount, nxny
  INTEGER :: ncolp,ncoltot,iradfmt,dualpol
  INTEGER :: nyqset,timeset,vardump

  INTEGER :: xscale
  REAL    :: rdx
  INTEGER :: dtime
  REAL    :: dsort
  REAL    :: refelvmin,refelvmax
  REAL    :: refrngmin,refrngmax

  INTEGER :: remapopt
  REAL :: radius

  REAL :: xrad,yrad

  REAL :: dBA

  REAL :: dpatrulat(2)

  INTEGER :: maxinfile,infilenum
  LOGICAL :: assigned_vil, assigned_et

  INTEGER :: iarg,jarg,narg,ifile,kfile
  INTEGER :: varfill,ivar
  CHARACTER(LEN=4)   :: tsthead
  CHARACTER(LEN=256) :: charg
  CHARACTER(LEN=256) :: listfile

  CHARACTER (LEN=6)  :: varid
  CHARACTER (LEN=20) :: varname
  CHARACTER (LEN=20) :: varunits

!------------------------------------------------------------------------
!
! Namelist declaration
!
!------------------------------------------------------------------------

  CHARACTER(LEN=4)   :: radar_name
  INTEGER            :: nnidsfn
  CHARACTER(LEN=256) :: nidsfn(maxnidsfile)
  CHARACTER(LEN=256) :: dir_name
  CHARACTER(LEN=256) :: etfn,vilfn,dpafn
  INTEGER :: et_remapopt,vil_remapopt,dpa_remapopt,crf_remapopt
  REAL    :: et_radius,vil_radius,dpa_radius,crf_radius
  INTEGER :: arbvaropt
  INTEGER :: adjreflopt

  NAMELIST /nids_data/ radar_name,nnidsfn,nidsfn, dir_name, &
                       et_remapopt, et_radius,              &
                       vil_remapopt,vil_radius,             &
                       dpa_remapopt,dpa_radius,             &
                       crf_remapopt,crf_radius, arbvaropt,adjreflopt

  INTEGER, PARAMETER :: ROOT = 0

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

  radband = 1
  isource = 2                  ! 1 - WSR-88D raw     2 - WSR-88D NIDS

  refelvmin = 20.0
  refelvmax = -1.0
  refrngmin = rngmin
  refrngmax = rngmax
  dsort = (sortmax-sortmin)/float(nsort-1)
  dmpfmt=1
  hdf4cmpr=0
  unitin = 1
  i_first_scan = .TRUE.
  noclrv = .FALSE.
  velproc = .TRUE.
  dualpproc = .FALSE.
  qc_on = .TRUE.
  ref2d = .FALSE.
  ref3d = .FALSE.
  vel2d = .FALSE.
  vel3d = .FALSE.
  ncdchdr = .FALSE.
  myproc = 0
  mp_opt = 0

!------------------------------------------------------------------------
!
! Initialize nids_data NAMELIST variables
!
!------------------------------------------------------------------------

  radar_name = 'DMMY'
  nnidsfn = 0
  nlistfil = 0
  ncdcfile = 0
  dir_name = './'
  len_dir = LEN(TRIM(dir_name))
  etfn = 'dummy'
  vilfn = 'dummy'
  dpafn = 'dummy'
  et_remapopt = 1
  et_radius = 4000.
  vil_remapopt = 1
  vil_radius = 4000.
  dpa_remapopt = 1
  dpa_radius = 4000.
  crf_remapopt = 1
  crf_radius = 4000.
  arbvaropt = 1
  adjreflopt = 0
  medfilt = 0

  DO ifile=1,maxnidsfile
    nidsfn(ifile)='dummy'
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
        WRITE(stdout,'(a,a,a,a)') &
        ' Usage: nids2arps [RADAR_ID] [fname] [-novel] [-dir dirname] [-hdf 2]',&
        ' [-binary] [-ref2d] [-ref3d] [-reffile] [-vel2d] [-vel3d] [-velfile]',&
        ' [-noqc] [-fillref] [-fntime] [-rtopts] [-input nmfile]', &
        ' < radremap.input'
        istatus = -1
      END IF
    END IF
  END IF
  CALL mpupdatei(istatus,1)
  CALL mpupdatei(narg,1)

  IF(istatus < 0) THEN
    CALL arpsstop ('.',istatus)
  END IF
!------------------------------------------------------------------------
!
! Read the input file and set logical variables based on input integers.
!
!------------------------------------------------------------------------

  infile = ' '
  IF(myproc == ROOT) THEN
    IF(narg > 0) THEN
      iarg=1
      DO jarg=1,narg
        CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
        IF(charg(1:6) == '-input') THEN    ! namelist file
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg, istat, istatus )
            infile = charg
          END IF
          WRITE(stdout,'(1x,2a)') 'Namelist file from command line - ',TRIM(infile)
          EXIT
        END IF
        iarg=iarg+1
        IF(iarg > narg) EXIT
      END DO
    END IF ! narg > 0
  END IF ! myproc == 0

  CALL initremapopt(nx,ny,nz,nzsoil,nstyps,1,1,1,infile,istatus)

  IF(myproc == ROOT) THEN

    ref2d = (ref2dopt > 0)
    ref3d = (ref3dopt > 0)
    vel2d = (vel2dopt > 0)
    vel3d = (vel3dopt > 0)
    velproc = (velprocopt > 0)
    qc_on = (qcopt > 0)
    traditional = (tradopt > 0)
    vad = (vadopt > 0)
    wrttilt = (wttiltopt > 0)
    fntime = (fntimopt > 0)
    ncdchdr = (ncdcfile > 0)

    WRITE(6,'(a,i4)') ' Obtained from namelists: dmpfmt=',dmpfmt
    WRITE(6,'(26x,a,a4)') 'radname= ',radname
    WRITE(6,'(26x,a,i4)') 'hdf4cmpr=',hdf4cmpr
    WRITE(6,'(26x,a,i4)') 'nlistfil=',nlistfil
    radar_name=radname

    IF(nlistfil > 0) THEN
      kfile=0
      DO ilist=1,nlistfil
        WRITE(stdout,'(a,a)')  ' Reading file lists ',listfil(ilist)
        OPEN(31,file=listfil(ilist),iostat=iopstat,status='old',form='formatted')
        IF(iopstat /= 0) THEN
          WRITE(stdout,'(a,a)') ' Error opening list file: ',listfil(ilist)
          CALL arpsstop('Error in list file.',iopstat)
        END IF
        DO ifile=(nnidsfn+1),maxnidsfile
          READ(31,'(a)',iostat=irdstat) nidsfn(ifile)
          IF(irdstat /= 0) EXIT
          kfile=kfile+1
          WRITE(stdout,'(a,i4,a,a)')                                      &
                 ' kfile:',kfile,'   nidsfn:',TRIM(nidsfn(ifile))
        END DO
        CLOSE(31)
        nnidsfn=kfile
      END DO
    END IF
    WRITE(6,'(a,i4)') '   nnidsfn=',nnidsfn

!------------------------------------------------------------------------
!
! Process command line
! For backward comptability allow input of radar file names via
! files specified on the command line
!
!------------------------------------------------------------------------

    istatus = 0
    WRITE(stdout,'(a,i4)') ' Number of command-line arguments: ',narg
    IF(narg > 1 ) THEN
      CALL GET_COMMAND_ARGUMENT(1, charg )
      radar_name=charg(1:4)
      WRITE(stdout,'(a,a)')  '    radar_name = ', radar_name
      kfile=0
      fillref=0
      iarg=2
      DO jarg=2,narg
        CALL GET_COMMAND_ARGUMENT(iarg, charg )
        IF(charg(1:6) == '-novel') THEN
          velproc=.FALSE.
        ELSE IF(charg(1:4) == '-hdf') THEN
          dmpfmt=2
          hdf4cmpr=0
          IF(iarg < narg) THEN
            iarg=iarg+1
            CALL GET_COMMAND_ARGUMENT(iarg, charg )
            READ(charg,'(i1)') hdf4cmpr
            hdf4cmpr=min(max(hdf4cmpr,0),7)
          END IF
          WRITE(stdout,'(a,i2)')                                        &
          ' Output in hdf format with compression level: ',hdf4cmpr
        ELSE IF(charg(1:7) == '-binary') THEN
          dmpfmt=1
          hdf4cmpr=0
          WRITE(stdout,'(a)') ' Output in binary format'
        ELSE IF(charg(1:5) == '-ncdc') THEN
          ncdchdr=.TRUE.
          WRITE(stdout,'(a)') ' Files are NCDC NIDS files'
        ELSE IF(charg(1:6) == '-ref2d') THEN
          ref2d=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will produce 2d reflectivity file of lowest tilt'
        ELSE IF(charg(1:6) == '-ref3d') THEN
          ref3d=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will produce 3d reflectivity file for plotting'
        ELSE IF(charg(1:8) == '-reffile') THEN
          ref3d=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will produce 3d reflectivity file for plotting'
        ELSE IF(charg(1:8) == '-fillref') THEN
          fillref=1
          WRITE(stdout,'(a)')                                           &
            ' Will fill-in reflectivity below lowest tilt'
        ELSE IF(charg(1:6) == '-vel2d') THEN
          vel2d=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will produce 2d velocity file of lowest tilt'
        ELSE IF(charg(1:6) == '-vel3d') THEN
          vel3d=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will produce 3d velocity file for plotting'
        ELSE IF(charg(1:8) == '-velfile') THEN
          vel3d=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will produce 3d velocity file for plotting'
        ELSE IF(charg(1:7) == '-noclrv') THEN
          noclrv=.TRUE.
          WRITE(stdout,'(a)')                                           &
            ' Will not velocities for clear-air VCPs'
        ELSE IF(charg(1:5) == '-noqc') THEN
          qc_on=.FALSE.
          WRITE(stdout,'(a)')                                           &
            ' Will skip qc steps in processing'
        ELSE IF(charg(1:8) == '-medfilt') THEN
          medfilt=1
          WRITE(stdout,'(a)')                                           &
            ' Apply median filter to az-ran data in despeckle step.'
        ELSE
          listfile=charg
          WRITE(stdout,'(a,a)')  ' Reading file lists ',charg
          OPEN(31,file=listfile,status='old',form='formatted')
          DO ifile=(nnidsfn+1),maxnidsfile
            READ(31,'(a)',END=101) nidsfn(ifile)
            kfile=kfile+1
            WRITE(stdout,'(a,i4,a,a)')                                      &
                 ' kfile:',kfile,' nidsfn: ',TRIM(nidsfn(ifile))
          END DO
          101 CONTINUE
          nnidsfn=kfile
        END IF
        iarg=iarg+1
        IF(iarg > narg) EXIT
      END DO
      WRITE(stdout,'(1x,a,i4)') 'Found number of file names - ',nnidsfn
    END IF
  END IF  ! myproc is ROOT

  CALL mpupdatei(istatus,1)
  IF (istatus < 0) CALL arpsstop(' ',1)

  CALL mpupdatec(radar_name,4)
  CALL mpupdatei(nnidsfn,1)
  CALL mpupdatec(nidsfn,256*maxnidsfile)
  CALL mpupdatei(fillref,1)
  CALL mpupdatel(velproc,1)
  CALL mpupdatei(dmpfmt,1)
  CALL mpupdatei(hdf4cmpr,1)
  CALL mpupdatel(ncdchdr,1)
  CALL mpupdatel(ref2d,1)
  CALL mpupdatel(ref3d,1)
  CALL mpupdatel(vel2d,1)
  CALL mpupdatel(vel3d,1)
  CALL mpupdatel(qc_on,1)

!------------------------------------------------------------------------
!
! Array allocations and variable initializations.
!
!------------------------------------------------------------------------

  ALLOCATE(kntrgat(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:kntrgat')
  ALLOCATE(kntrazm(maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:kntrazm')

  ALLOCATE(kntvgat(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:kntvgat')
  ALLOCATE(kntvazm(maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:kntvazm')

  ALLOCATE(nyqvvol(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:nyqvvol')
  ALLOCATE(timevolr(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:timevolr')
  ALLOCATE(timevolv(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:timevolv')

  ALLOCATE(rngrvol(maxgate,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:rngrvol')
  ALLOCATE(azmrvol(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:azmrvol')
  ALLOCATE(elvrvol(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:elvrvol')
  ALLOCATE(elvmnrvol(maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:elvmnrvol')
  ALLOCATE(refvol(maxgate,maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:refvol')

  ALLOCATE(rngvvol(maxgate,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:rngvvol')
  ALLOCATE(azmvvol(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:azmvvol')
  ALLOCATE(elvvvol(maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:elvvvol')
  ALLOCATE(elvmnvvol(maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:elvmnvvol')
  ALLOCATE(velvol(maxgate,maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:velvol')

  ALLOCATE(rxvol(maxgate,maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:rxvol')
  ALLOCATE(ryvol(maxgate,maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:ryvol')
  ALLOCATE(rzvol(maxgate,maxazim,maxelev), STAT=istat)
  CALL check_alloc_status(istat,'nids2arps:rzvol')

!-----------------------------------------------------------------------
!
! Allocation of analysis arrays
!
!-----------------------------------------------------------------------

  nxny=nx*ny
  ALLOCATE(icolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:icolp')
  ALLOCATE(jcolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:jcolp')
  ALLOCATE(xcolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:xcolp')
  ALLOCATE(ycolp(nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:ycolp')
  ALLOCATE(zcolp(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:zcolp')
  ALLOCATE(havdat(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:havdat')

  ALLOCATE(kntbin(nsort),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:kntbin')

  ALLOCATE(colref(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:colref')
  ALLOCATE(colvel(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:colvel')
  ALLOCATE(colnyq(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:colnyq')
  ALLOCATE(coltim(nz,nxny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:coltim')

  ALLOCATE(x(nx))
  CALL check_alloc_status(istat,'nids2arps:x')
  ALLOCATE(y(ny))
  CALL check_alloc_status(istat,'nids2arps:y')
  ALLOCATE(z(nz))
  CALL check_alloc_status(istat,'nids2arps:z')
  ALLOCATE(zp(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:zp')
  ALLOCATE(zpsoil(nx,ny,nzsoil))
  CALL check_alloc_status(istat,'nids2arps:zpsoil')
  ALLOCATE(xs(nx))
  CALL check_alloc_status(istat,'nids2arps:xs')
  ALLOCATE(ys(ny))
  CALL check_alloc_status(istat,'nids2arps:ys')
  ALLOCATE(zps(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:zps')
  ALLOCATE(hterain(nx,ny))
  CALL check_alloc_status(istat,'nids2arps:hterain')
  ALLOCATE(mapfct(nx,ny,8))
  CALL check_alloc_status(istat,'nids2arps:mapfct')
  ALLOCATE(j1(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:j1')
  ALLOCATE(j2(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:j2')
  ALLOCATE(j3(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:j3')
  ALLOCATE(j3soil(nx,ny,nzsoil))
  CALL check_alloc_status(istat,'nids2arps:j3soil')
  ALLOCATE(j3soilinv(nx,ny,nzsoil))
  CALL check_alloc_status(istat,'nids2arps:nzsoil')

  ALLOCATE(tem1d1(nz))
  CALL check_alloc_status(istat,'nids2arps:tem1d1')
  ALLOCATE(tem1d2(nz))
  CALL check_alloc_status(istat,'nids2arps:tem1d2')
  ALLOCATE(tem2d(nx,ny))
  CALL check_alloc_status(istat,'nids2arps:tem2d')
  ALLOCATE(tem3d(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:tem3d')

  ALLOCATE(grdvel2(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:grdvel2')
  ALLOCATE(grdref2(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:grdref2')
  ALLOCATE(wgtvel(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:wgtvel')
  ALLOCATE(wgtref(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:wgtref')
  ALLOCATE(wpotvel(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:wpotvel')
  ALLOCATE(wpotref(nx,ny,nz))
  CALL check_alloc_status(istat,'nids2arps:wpotref')
  ALLOCATE(gridazm(nx,ny))
  CALL check_alloc_status(istat,'nids2arps:gridazm')
  ALLOCATE(gridrng(nx,ny))
  CALL check_alloc_status(istat,'nids2arps:gridrng')

  ALLOCATE(rdata(maxgate,maxazim))
  CALL check_alloc_status(istat,'nids2arps:rdata')
  ALLOCATE(rtem(maxgate,maxazim))
  CALL check_alloc_status(istat,'nids2arps:rtem')
  ALLOCATE(time(maxazim))
  CALL check_alloc_status(istat,'nids2arps:time')
  ALLOCATE(azim(maxazim))
  CALL check_alloc_status(istat,'nids2arps:azim')
  ALLOCATE(vnyq(maxazim))
  CALL check_alloc_status(istat,'nids2arps:vnyq')
  ALLOCATE(elev(maxazim))
  CALL check_alloc_status(istat,'nids2arps:elev')

  IF (adjreflopt == 1) THEN
    ALLOCATE(nidset(nx,ny))
    CALL check_alloc_status(istat,'nids2arps:nidset')
    ALLOCATE(nidsvil(nx,ny))
    CALL check_alloc_status(istat,'nids2arps:nidsvil')
    nidset(:,:) = -9999.
    nidsvil(:,:) = -9999.
  END IF
  assigned_vil = .FALSE.
  assigned_et = .FALSE.
!
! Fill refractivity sounding with constant lapse rate
! This is only used when rfropt > 1, but for completeness and
! future upgrade (when an actual sounding made from gridded data
! would be used) this is filled-in.
!
  DO k=1,nzsnd
    zsnd(k)=(k-1)*dzsnd
    rfrsnd(k)=325.+rfrconst*zsnd(k)
  END DO

  icount = 0
  rfirstg = 1000
  gatesp = 1000

  CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                  &
              hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,         &
              tem1d1,tem1d2,tem3d)

  DEALLOCATE(mapfct, j1, j2, j3)
  DEALLOCATE(j3soil, j3soilinv, zpsoil)

  ALLOCATE(arpslat(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:arpslat')
  ALLOCATE(arpslon(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:arpslon')
  ALLOCATE(arpsrasterx(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:arpsrasterx')
  ALLOCATE(arpsrastery(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:arpsrastery')
  ALLOCATE(arpsdpax(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:arpsdpax')
  ALLOCATE(arpsdpay(nx,ny),stat=istat)
  CALL check_alloc_status(istat,'nids2arps:arpsdpay')

!------------------------------------------------------------------------
!
! Loop through the NIDS files.
!
!------------------------------------------------------------------------

  IF (nnidsfn > 0 .AND. nnidsfn <= maxnidsfile) THEN
    maxinfile = nnidsfn
  ELSE IF (nnidsfn > maxnidsfile) THEN
    WRITE(stdout,'(a,i4,a,i4)') &
      ' WARNING: nnidsfn=', nnidsfn,' > maxnidsfile=',maxnidsfile
    WRITE(stdout,'(a,i4,a)') &
      ' Will only read in ', maxnidsfile,' NIDS files.'
    maxinfile = maxnidsfile
  ELSE
    maxinfile = 0
  END IF

  IF (maxinfile > 0) THEN
    DO infilenum = 1,maxinfile
      infile = nidsfn(infilenum)

      IF (myproc == ROOT) THEN
        nazim = 0
        ielev = 0

        cval(:) = ' ' ! Initialize array
        ival(:) = 0   ! Initialize array

        WRITE(6,'(/,1x,a,I4,a,/,1x,2a)')                                &
                '==== No. ',infilenum,' ==============================',&
                'File name - ',TRIM(infile)
        unitin = 12
        OPEN (UNIT=unitin, FILE=TRIM(infile), STATUS='OLD',             &
              ACCESS='DIRECT',RECL=1, FORM='UNFORMATTED', ERR=9981,     &
              IOSTAT=ios)
!
!------------------------------------------------------------------------
!
!     If this file contains an extra NCDC header the radar name should
!     be in character records 8-11.  Otherwise, those records will contain
!     something else.
!
!------------------------------------------------------------------------
!
        DO i = 1, 4
          n=7+i
          READ (UNIT=unitin, REC=n, IOSTAT=ios) cval(i)
          IF (ios /= 0) EXIT
        END DO
        WRITE(tsthead,'(4a1)') (cval(i),i=1,4)
!       WRITE(6,'(a,a4)') '  Test header string is ',tsthead
        IF( ncdchdr .OR. tsthead == radar_name ) THEN
          WRITE(6,'(a)') ' This is a file from NCDC, skipping 30 bytes'
          iskip=30
        ELSE
          WRITE(6,'(a)')                                                &
           ' File does not contain extra NCDC header, no skip.'
          iskip=0
        END IF
!
!------------------------------------------------------------------------
!
!     Read product into array ival() below.
!
!------------------------------------------------------------------------
!
        DO i = 1, maxlen
          n=iskip+i
          READ (UNIT=unitin, REC=n, IOSTAT=ios) cval(i)
          IF (ios /= 0) EXIT
        END DO
        CLOSE(UNIT=unitin)

        msglen = i - 1
        DO i = 1, msglen
          ival(i) = ICHAR(cval(i))
        END DO

        WRITE (stdout,'(1x,a)')   'Finished reading ...'
        WRITE (stdout,'(1x,a,i5)', ADVANCE='NO') '  System code: ',ios
        WRITE (stdout,'(a,i10)') ',    Message length: ',msglen

        GOTO 9980       ! Success reading file

        9981 CONTINUE   ! Error opening file.

        WRITE(stdout,'(/a,a,/a)')' ERROR:  Cannot open file ',TRIM(infile),    &
                     ' Aborting...'
        istatus = -1

        9980 CONTINUE

      END IF

      CALL mpupdatei(istatus, 1)
      IF (istatus < 0) CALL arpsstop(' ',1)

      CALL mpupdatei(msglen,1)
      CALL mpupdatei(ival,msglen)

!------------------------------------------------------------------------
!
!     Get product ID code.
!
!------------------------------------------------------------------------

      ival0 = ival(1);   ival1 = ival(2)
      IF (ival1 < 0) ival1 = 256 + ival1
      IF (ival0 < 0) ival0 = 256 + ival0
      prodID = ival1 + (256*ival0)
      IF (myproc == ROOT) THEN
         WRITE(stdout,'(1x,a,i5)',ADVANCE='NO') '  Product  ID: ',prodID
      END IF

!------------------------------------------------------------------------
!
!     Get product type code.
!
!------------------------------------------------------------------------

      ival0 = ival(137);  ival1 = ival(138)
      IF (ival1 < 0) ival1 = 256 + ival1
      IF (ival0 < 0) ival0 = 256 + ival0
      prodtype = ival1 + (256*ival0)
      IF (myproc == ROOT) THEN
        WRITE (stdout,FMT='(a,i10)') ',    Type      code: ', prodtype
      END IF

      icode = -1   ! It should be reset to 0 if everything is ok below
!------------------------------------------------------------------------
!
!     Call appropriate decoding subroutine based on product type
!     (radial, raster, or DPA).  Note that the hybrid scan reflectivity
!     routine could be added here in the future.
!
!------------------------------------------------------------------------

      print *, ' prodtype =', prodtype

!     IF (prodtype == ipradial1 .OR. prodtype == ipradial2 .OR. prodtype == ipradial3 .OR. &
!         prodtype == ipradial4 .OR. prodtype == ipradial5 ) THEN

      IF (prodtype == ipradial1) THEN

        IF (myproc == ROOT) WRITE(stdout,'(1x,a)')'About to call get_radial_data ...'

        CALL get_radial_data(iheaderdim,itrailerdim,msglen,ival,        &
                             isite,isite_lat,isite_lon,isite_elev,      &
                             iyear,imonth,iday,ihr,imin,isec,           &
                             iprod,maxval1,ivcp,ifield,maxgate,maxazim, &
                             iazmuth,ielev,ngate,nazim,icatt,icats,     &
                             maxval2,iyr1,imonth1,iday1,ihr1,imin1,     &
                             iyr2,imonth2,iday2,ihr2,imin2,istm_spd,    &
                             istm_dir,iheader,len_header,itrailer,      &
                             len_trailer,icode)

        IF (myproc == ROOT) THEN
          IF (icode == 2) THEN
            WRITE(stdout,'(a)')' ERROR -- Not a radial graphic product.'
          ELSE IF (icode == 3) THEN
            WRITE(stdout,'(a)')' ERROR -- ncols/nrows too small for product.'
          ELSE IF (icode == 4) THEN
            WRITE(stdout,'(a)')' ERROR -- msglen too small to store product.'
          END IF
        END IF

      ELSE IF (prodtype == ipraster1 .OR. prodtype == ipraster2) THEN

        IF (myproc == ROOT) WRITE(stdout,'(a)')'About to call get_raster_data...'

        CALL get_raster_data(iheaderdim,itrailerdim,msglen,ival,        &
                           isite,isite_lat,isite_lon,isite_elev,        &
                           iyear,imonth,iday,ihr,imin,isec,             &
                           iprod,maxval1,ivcp,ifield_raster,ndim1,ndim2,&
                           icatt,icats,ndims_p,iheader,len_header,      &
                           itrailer,len_trailer,icode)

        IF (myproc == ROOT) THEN
          IF (icode == 2) THEN
            WRITE(stdout,'(a)')'ERROR -- Not a raster graphic product.'
          ELSE IF (icode == 3) THEN
            WRITE(stdout,'(a)')'ERROR -- ncols/nrows too small for product.'
          ELSE IF (icode == 4) THEN
            WRITE(stdout,'(a)')'ERROR -- msglen too small to store product.'
          END IF
        END IF

      ELSE IF (prodID == DPAprodID) THEN

        IF (myproc == ROOT) WRITE(stdout,'(a)')'About to call get_dpa_data...'

        CALL get_dpa_data(iheaderdim,itrailerdim,msglen,ival,           &
                          isite,isite_lat,isite_lon,isite_elev,         &
                          iyear,imonth,iday,ihr,imin,isec,              &
                          iprod,maxval1,ivcp,ifield_dpa,ndim1,ndim2,    &
                          ncolsp,nrowsp,iheader,len_header,itrailer,    &
                          len_trailer,icode)

        IF (myproc == ROOT) THEN
          IF (icode == 2) THEN
            WRITE(stdout,'(a)')'ERROR -- Not a DPA product.'
          ELSE IF (icode == 3) THEN
            WRITE(stdout,'(a)')                                         &
                'ERROR -- ncols/nrows too small for product.'
          ELSE IF (icode == 4) THEN
            WRITE(stdout,'(a)')                                         &
                  'ERROR -- msglen too small to store product.'
          END IF
        END IF

      ELSE
        IF (myproc == ROOT) WRITE(stdout,'(a,/a)')                      &
            'ERROR: Current file does not contain NIDS data.',          &
            'Moving to next file...'
      END IF

      IF (icode /= 0) CYCLE

!------------------------------------------------------------------------
!
!     Save radar lat/lon, elevation, and number of radials.
!
!------------------------------------------------------------------------

      radalt = isite_elev*mpfoot
      radlat = isite_lat*0.001
      radlon = isite_lon*0.001

      eleva = ielev*0.1

      IF (myproc == ROOT) THEN
        WRITE(stdout,'(1x,a)') 'Radar information: '
        WRITE(stdout,'(1x,a,a,a,f8.1)', ADVANCE='NO') '  NAME = ', radar_name,       &
                                     ', ALTITUDE = ', radalt
        WRITE(stdout,'(a,f10.4,a,f10.4)') ', LATITUDE = ', radlat,       &
                                          ', LONGITUDE = ', radlon

        WRITE(stdout,'(a,i3,a,i4)')  &
                   ('   Threshold ', i, ' icats = ', icats(i), i=0,15)
      END IF

      DO i=0,15
        rcats(i)=float(icats(i))
      END DO

!------------------------------------------------------------------------
!
!     VCP Checking.
!     As a quick fix for RT-2008 exit if in clear-air mode.
!     Later use commend line option to make a decision here.
!
!------------------------------------------------------------------------

!      IF (myproc == ROOT) WRITE(stdout,'(1x,a,i4)') ' VCP =',ivcp

      CALL mpupdatei(ivcp,1)

      IF( noclrv .AND. (ivcp == 31 .OR. ivcp == 32)) THEN
!        IF (myproc == ROOT) WRITE(stdout,'(1x,a,i4,a)') &
!             ' Clear-air VCP detected:',ivcp,' exiting'
!        CALL mpexit(0)
!        CALL arpsstop(' ',2)
         IF (myproc == ROOT) WRITE(stdout,'(1x,a,i4)') &
              ' Clear-air VCP detected:',ivcp
         velproc=.FALSE.
         CALL mpupdatel(velproc,1)
      END IF
!------------------------------------------------------------------------
!
!     Initialize 3-d radar data arrays
!
!------------------------------------------------------------------------

      CALL ctim2abss(iyear, imonth, iday, ihr, imin, isec, itime)

      IF (i_first_scan) THEN

        CALL radcoord(nx,ny,nz,                                         &
            x,y,z,zp,xs,ys,zps,                                         &
            radlat,radlon,radarx,radary)

        print *, ' calling rmpinit'
        CALL rmpinit(nx,ny,nz,maxgate,maxgate,maxazim,maxelev,dualpproc,&
                kntrgat,kntrazm,kntrelv,                                &
                kntvgat,kntvazm,kntvelv,                                &
                nyqvvol,timevolr,timevolv,                              &
                rngrvol,azmrvol,elvrvol,elvmnrvol,                      &
                refvol,rhvvol,zdrvol,kdpvol,                            &
                rngvvol,azmvvol,elvvvol,elvmnvvol,velvol,               &
                colvel,colref,colnyq,coltim)
        print *, ' back from rmpinit'

        iiyear = iyear
        iimonth = imonth
        iiday = iday
        iihr = ihr
        iimin = imin
        iisec = isec

        itimfrst=itime

        i_first_scan = .FALSE.

      END IF

!------------------------------------------------------------------------
!
!     Now begin processing the data.
!
!------------------------------------------------------------------------

!     IF (prodtype == ipradial1 .OR. prodtype == ipradial2 .OR. prodtype == ipradial3 .OR. &
!         prodtype == ipradial4 .OR. prodtype == ipradial5 ) THEN

      IF (prodtype == ipradial1 ) THEN

        icount = icount + 1

        IF (myproc == ROOT) THEN
          WRITE(stdout,'(1x,a,i4)',ADVANCE='NO') 'VCP number for this file: ', ivcp
          WRITE(stdout,'(a,i4.4,a1,i2.2,a1,i2.2)', ADVANCE='NO')     &
                 ', DATE: ', iyear, '/', imonth, '/', iday
          WRITE(stdout,'(a,i2.2,a1,i2.2,a1,i2.2)')                   &
                 ', TIME: ', ihr, ':', imin, ':', isec
        END IF
!------------------------------------------------------------------------
!
!       Process base reflectivity data.
!
!------------------------------------------------------------------------

        print *, ' prodID = ',prodID
        IF ( (prodID >= BREFprodID1 .AND. prodID <= BREFprodID2) .OR. &
              prodID == TDWRrefID1 .OR.  prodID == TDWRrefID2) THEN

          IF (myproc == ROOT) THEN
            WRITE(stdout,'(a,i4)') &
               ' Processing base reflectivity data... ', icount

            WRITE(stdout,*)'Transferring ',nazim,' reflectivity radials.'
          END IF

          DO j = 1, maxazim
            DO i = 1, maxgate
              IF (ifield(i,j) > 0 .AND. ifield(i,j) < 16) THEN
                rdata(i,j) = rcats(ifield(i,j))
              ELSE
                rdata(i,j) = r_missing
              END IF
            END DO
          END DO

          dtime=itime-itimfrst
          refelvmin=min(eleva,refelvmin)
          refelvmax=max(eleva,refelvmax)
          DO j = 1, maxazim
            time(j) = dtime
            azim(j) = iazmuth(j)*0.10
            elev(j) = eleva
          END DO

          IF (myproc == ROOT) THEN
            DO j = 1, nazim, 20
              WRITE(stdout,'(a,i5,a,f8.1)') '   Ref j =',j,'  azim =',azim(j)
            END DO
          END IF

          IF (iprod < 0 .OR. iprod > 200) THEN
            gatesp  = 0                 ! iprod: product code number
            rfirstg = 0                 ! should be product->code
          ELSE                          ! in struct Proddesc?
            gatesp  = NINT(1000.0*res(iprod))
            rfirstg = NINT(1000.0*res(iprod))
          END IF

!------------------------------------------------------------------------
!
!       Call radar volume builder.
!
!------------------------------------------------------------------------
          irefgatsp=gatesp

          IF( qc_on) THEN
            IF (myproc == ROOT) print *, ' calling despekl '
            CALL despekl(maxgate,maxazim,maxgate,maxazim,refchek,        &
                         medfilt,rdata,rtem)
          END IF

          nyqset=0
          IF ( velproc ) THEN
            timeset=0
          ELSE
            timeset=1
          END IF
          vnyq=0.
          IF (myproc == ROOT) print *, 'Calling volbuild reflectivity - nazim =',nazim
          CALL volbuild(maxgate,maxazim,maxelev,ngate,nazim,             &
                 nyqset,timeset,                                         &
                 kntrgat,kntrazm,kntrelv,                                &
                 gatesp,rfirstg,refchek,                                 &
                 vnyq,time,                                              &
                 azim,elev,rdata,                                        &
                 nyqvvol,timevolr,                                       &
                 rngrvol,azmrvol,elvrvol,elvmnrvol,refvol)

!------------------------------------------------------------------------
!
!       Process velocity data.
!
!------------------------------------------------------------------------

        ELSE IF ( (prodID >= BVELprodID1 .AND. prodID <= BVELprodID2) .OR.   &
                   prodID == TDWRvelID ) THEN

          IF (myproc == ROOT) THEN
            WRITE(stdout,'(a,i4)')' Processing base velocity data ... ', icount
            WRITE(stdout,'(a,i6,a)')' Transferring ',nazim,' velocity radials.'
          END IF

          DO i = 1,14
            rcats(i)=kts2ms*rcats(i)
          END DO
          vnyq=rcats(14)
          IF (myproc == ROOT) WRITE(stdout,'(a,f10.2)') ' Nyquist velocity: ',vnyq(1)

          DO j = 1, maxazim
            DO i = 1, maxgate
              IF (ifield(i,j) > 0 .AND. ifield(i,j) < 15) THEN
                rdata(i,j) = rcats(ifield(i,j))
              ELSE
                rdata(i,j) = r_missing
              END IF
            END DO
          END DO

          dtime=(itime-itimfrst)
          DO j = 1, maxazim
            time(j) = dtime
            azim(j) = 0.1* float(iazmuth(j))
            elev(j) = eleva
          END DO

          IF (myproc == ROOT) THEN
            DO j = 1, nazim, 20
              WRITE(stdout,'(a,i5,a,f8.1)') '   Vel j =',j,'  azim =',azim(j)
            END DO
          END IF

          IF (iprod < 0 .OR. iprod > 200) THEN
            gatesp = 0                    ! iprod: product code number
            rfirstg = 0                   ! should be product->code
          ELSE                            ! in struct Proddesc?
            gatesp = NINT(1000.0*res(iprod))
            rfirstg = NINT(1000.0*res(iprod))
          END IF

!------------------------------------------------------------------------
!
!       Call radar volume builder.
!
!------------------------------------------------------------------------
          ivelgatsp=gatesp

          IF( qc_on ) THEN
            IF (myproc == ROOT) print *, ' calling despekl '
            CALL despekl(maxgate,maxazim,maxgate,maxazim,velchek,       &
                         medfilt,rdata,rtem)
          END IF

          nyqset=1
          timeset=1
          IF (myproc == ROOT) print *, 'Calling volbuild velocity - nazim =',nazim
          CALL volbuild(maxgate,maxazim,maxelev,maxgate,nazim,          &
                 nyqset,timeset,                                        &
                 kntvgat,kntvazm,kntvelv,                               &
                 gatesp,rfirstg,velchek,                                &
                 vnyq,time,                                             &
                 azim,elev,rdata,                                       &
                 nyqvvol,timevolv,                                      &
                 rngvvol,azmvvol,elvvvol,elvmnvvol,velvol)

!------------------------------------------------------------------------
!
!       Other product ID code for radial data.
!
!------------------------------------------------------------------------

        ELSE
          IF (myproc == ROOT) WRITE(stdout,'(a,i6)')                    &
                ' ERROR -- Unknown NIDS product ID code: ', prodID
          CYCLE
        END IF

      ELSE IF (prodtype == ipraster1 .OR. prodtype == ipraster2) THEN

!------------------------------------------------------------------------
!
!       Save resolution of raster grid.
!
!------------------------------------------------------------------------

        IF (myproc == ROOT) print *, ' ival(147), ival(148): ',ival(147),ival(148)
        ival1 = ival(148)
        IF (ival1 < 0) ival1 = 256 + ival1
        ival0 = ival(147)
        IF (ival0 < 0) ival0 = 256 + ival0
        xscale = ival1 + (256*ival0)
        IF (myproc == ROOT) print *, ' ival0, ival1: ',ival0, ival1
        IF( xscale == 2 ) THEN
          rdx = 4000.
        ELSE
          rdx = xscale*1000.
        END IF
        IF (myproc == ROOT) print *, ' xscale, rdx: ',xscale,rdx
        IF (myproc == ROOT) print *, ' ndims_p: ',ndims_p

!------------------------------------------------------------------------
!
!       Allocate raster arrays and save decoded data.
!
!------------------------------------------------------------------------

        ALLOCATE(rasterdata(ndims_p,ndims_p),STAT=istatus)
        CALL check_alloc_status(istatus,'nids2arps:rasterdata')
        ALLOCATE(rasterx   (ndims_p,ndims_p),STAT=istatus)
        CALL check_alloc_status(istatus,'nids2arps:rasterx')
        ALLOCATE(rastery   (ndims_p,ndims_p),STAT=istatus)
        CALL check_alloc_status(istatus,'nids2arps:rastery')
        ALLOCATE(rasterlat (ndims_p,ndims_p),STAT=istatus)
        CALL check_alloc_status(istatus,'nids2arps:rasterlat')
        ALLOCATE(rasterlon (ndims_p,ndims_p),STAT=istatus)
        CALL check_alloc_status(istatus,'nids2arps:rasterlon')
        ALLOCATE(remapped_raster(nx,ny),     STAT=istatus)
        CALL check_alloc_status(istatus,'nids2arps:remapped_raster')

        rasterdata(:,:)      = rasterchek
        remapped_raster(:,:) = rastermiss

        IF (prodID == ETprodID) THEN ! Echo tops

          IF (myproc == ROOT) WRITE(stdout,'(a)') ' Processing NIDS Echo Top data...'
          remapopt = et_remapopt
          radius   = et_radius

          DO j = 1,ndims_p   ! Column
            DO i = 1,ndims_p ! Row
              rasterdata(i,j)=icats(ifield_raster(i,ndims_p-j+1))*    &
                                                          mpfoot*1000.
            END DO           ! Row
          END DO             ! Column

        ELSE IF (prodID == VILprodID) THEN ! VIL

          IF (myproc == ROOT) WRITE(stdout,'(a)') ' Processing NIDS VIL data...'
          remapopt = vil_remapopt
          radius = vil_radius

          DO j = 1,ndims_p   ! Column
            DO i = 1,ndims_p ! Row
              rasterdata(i,j) = icats(ifield_raster(i,ndims_p-j+1))
            END DO           ! Row
          END DO             ! Column

        ELSE IF (prodID == CREFprodID) THEN ! Composite Reflectivity

          IF (myproc == ROOT) WRITE(stdout,'(a)') ' Processing NIDS Composite Reflectivity...'
          remapopt = crf_remapopt
          radius = crf_radius

          DO j = 1,ndims_p   ! Column
            DO i = 1,ndims_p ! Row
              rasterdata(i,j) = icats(ifield_raster(i,ndims_p-j+1))
            END DO           ! Row
          END DO             ! Column

        ELSE
          IF (myproc == ROOT) THEN
            WRITE(stdout,'(a)') ' ERROR -- Unknown raster data type.'
            WRITE(stdout,'(a)') ' Moving to next file.'
          END IF

          CYCLE
        END IF

!------------------------------------------------------------------------
!
!       Calculate Cartesian coordinates of raster data on raster grid
!       relative to the radar, and set rasterdata to missing if it is
!       beyond radar range.
!
!------------------------------------------------------------------------

        CALL setrastergrid(ndims_p,radlat,radlon,rdx,                   &
                           rasterlat,rasterlon,rasterx,rastery,         &
                           rasterdata,rasterchek)
        xrad = 0. ! Radar at center of raster grid
        yrad = 0. ! Radar at center of raster grid

!------------------------------------------------------------------------
!
!       Determine Cartesian coordinates of ARPS scalar points on radar
!       grid relative to radar.
!
!------------------------------------------------------------------------

        CALL xytoll(nx,ny,xs,ys,arpslat,arpslon)
        CALL getarpsrasterxy(nx,ny,arpslat,arpslon,radlat,radlon,       &
                       arpsrasterx,arpsrastery)

!------------------------------------------------------------------------
!
!       Remap the data and write to file.
!
!------------------------------------------------------------------------

        CALL remapraster(nx,ny,arpsrasterx,arpsrastery,dx,rasterchek,   &
                         rastermiss,xrad,yrad,                          &
                         remapopt,radius,ndims_p,ndims_p,rasterdata,    &
                         rasterx,rastery,remapped_raster,istatus)

        IF (prodID == ETprodID) THEN ! Echo tops

          IF (arbvaropt == 1) THEN
            CALL wrtvar2(nx,ny,1,remapped_raster,radar_name//'et',      &
                       'NIDS Echo Tops','m',000000,runname(1:lfnkey),   &
                       TRIM(dir_name),dmpfmt,hdf4cmpr,joindmp(FINDX_A),istatus)
          END IF

          IF (adjreflopt == 1 .AND. .NOT. assigned_et) THEN
            nidset(:,:) = remapped_raster(:,:)
            assigned_et = .TRUE.
          ELSE IF (adjreflopt == 1 .AND. assigned_et) THEN
            IF (myproc == ROOT) WRITE(stdout,'(a)')                     &
              ' WARNING:  NIDS Echo Top product already remapped!'
          END IF

        ELSE IF (prodID == VILprodID) THEN ! VIL

          IF (arbvaropt == 1) THEN
            CALL wrtvar2(nx,ny,1,remapped_raster,radar_name//'vl',      &
                       'NIDS VIL','kg m^-2',000000,runname(1:lfnkey),   &
                       TRIM(dir_name),dmpfmt,hdf4cmpr,joindmp(FINDX_A),istatus)
          END IF

          IF (adjreflopt == 1 .AND. .NOT. assigned_vil) THEN
            nidsvil(:,:) = remapped_raster(:,:)
            assigned_vil = .TRUE.
          ELSE IF (adjreflopt == 1 .AND. assigned_vil) THEN
            IF (myproc == ROOT) WRITE(stdout,'(a)')                     &
              ' WARNING:  NIDS VIL product already remapped!'
          END IF

        ELSE IF (prodID == CREFprodID) THEN ! Composite Reflectivity

          IF (arbvaropt == 1) THEN
            CALL wrtvar2(nx,ny,1,remapped_raster,radar_name//'cref',    &
                       'NIDS CREF','dBZ',000000,runname(1:lfnkey),      &
                       TRIM(dir_name),dmpfmt,hdf4cmpr,joindmp(FINDX_A),istatus)
          END IF

        END IF

        DEALLOCATE(rasterdata,rasterx,rastery,rasterlat,rasterlon,      &
                   remapped_raster,STAT=istatus)

      ELSE IF (prodID == DPAprodID) THEN

!------------------------------------------------------------------------
!
!       Allocate arrays and save decoded data.
!
!------------------------------------------------------------------------

        IF (myproc == ROOT) WRITE(stdout,'(a)') ' Processing NIDS DPA data...'

        ALLOCATE(dpadata(nrowsp,ncolsp),STAT=istatus)
        ALLOCATE(dpax   (nrowsp,ncolsp),STAT=istatus)
        ALLOCATE(dpay   (nrowsp,ncolsp),STAT=istatus)
        ALLOCATE(dpalat (nrowsp,ncolsp),STAT=istatus)
        ALLOCATE(dpalon (nrowsp,ncolsp),STAT=istatus)
        ALLOCATE(remapped_dpa(nx,ny),   STAT=istatus)

        dpadata(:,:)      = dpachek
        remapped_dpa(:,:) = dpamiss

        remapopt = dpa_remapopt
        radius = dpa_radius

        DO j = 1,ncolsp
          DO i = 1,nrowsp

            IF (ifield_dpa(i,ncolsp-j+1) == 0) THEN
              dpadata(i,j) = 0.
            ELSE IF (ifield_dpa(i,ncolsp-j+1) /= 255) THEN
              dBA = -6.125 + (REAL(ifield_dpa(i,ncolsp-j+1))*0.125)
!              dpadata(i,j) = 10.**(dBA*0.1) ! Rainfall in mm.
              dpadata(i,j) = mm2in*10.**(dBA*0.1) ! Rainfall in inches
              dpadata(i,j) = MAX(0.,dpadata(i,j))
            END IF

          END DO ! i loop
        END DO ! j loop

!------------------------------------------------------------------------
!
!       Determine Cartesian coordinates of DPA points on DPA grid
!       relative to the radar.
!
!------------------------------------------------------------------------

        DO j = 1,nrowsp
          DO i = 1,ncolsp
!            dpax(i,j) = DPAdx*(REAL(i) - 66.)
!            dpay(i,j) = DPAdx*(REAL(j) - 66.)
            dpax(i,j) = DPAdx*(REAL(i) - 67.)   ! Better match with NWS
            dpay(i,j) = DPAdx*(REAL(j) - 67.)   ! graphics.
          END DO ! i loop
        END DO ! j loop

        xrad = 0. ! Radar at center of grid
        yrad = 0. ! Radar at center of grid

!------------------------------------------------------------------------
!
!       Determine Cartesian coordinates of ARPS scalar points in DPA map
!       projection relative to radar.
!
!------------------------------------------------------------------------

        CALL xytoll(nx,ny,xs,ys,arpslat,arpslon)

        dpatrulat(1) = DPAtrulat1
        dpatrulat(2) = DPAtrulat1 ! Not used for polar stereographic
        CALL setmapr(DPAiproj,DPAscale,dpatrulat,DPAtrulon)
        CALL setorig(2,radlat,radlon)
        CALL lltoxy(nx,ny,arpslat,arpslon,arpsdpax,arpsdpay)

!------------------------------------------------------------------------
!
!       Remap the DPA data and write to file.
!
!------------------------------------------------------------------------

        CALL remapraster(nx,ny,arpsdpax,arpsdpay,dx,dpachek,            &
                     dpamiss,xrad,yrad,                                 &
                     remapopt,radius,nrowsp,ncolsp,dpadata,             &
                     dpax,dpay,remapped_dpa,istatus)

        IF (arbvaropt == 1) THEN
          CALL wrtvar2(nx,ny,1,remapped_dpa,radar_name//'dp','NIDS DPA',&
                     'mm',000000,runname(1:lfnkey),TRIM(dir_name),      &
                     dmpfmt,hdf4cmpr,joindmp(FINDX_A),istatus)
        END IF

      END IF
    END DO ! DO ifilenum = 1,maxinfile
  END IF ! IF (maxinfile > 0)

  IF (icount > 0) THEN
!------------------------------------------------------------------------
!
!   Call remapping routines
!
!------------------------------------------------------------------------

    IF( qc_on ) THEN
      IF (myproc == ROOT) print *, ' Calling apdetect '
      CALL apdetect(maxgate,maxgate,maxazim,maxelev,                    &
                    kntrgat,kntrazm,kntrelv,                            &
                    kntvgat,kntvazm,kntvelv,                            &
                    refchek,velchek,                                    &
                    irefgatsp,ivelgatsp,                                &
                    winszrad,winszazim,ivcp,gclopt,gcvrlim,             &
                    rngrvol,azmrvol,elvrvol,                            &
                    rngvvol,azmvvol,elvvvol,                            &
                    refvol,velvol,rtem,                                 &
                    istatus)
    END IF

    nyqset=0
    IF( velproc ) THEN
      timeset=0
    ELSE
      timeset=1
    END IF
    IF(ref3d) THEN
      vardump=1
    ELSE
      vardump=0
    END IF
    varfill=fillref
    ivar=1
!
!   VAD Wind Processing
!
    IF( (myproc == ROOT) .AND. vad ) THEN

      WRITE(stdout,'(a)') 'VAD processing not enabled in nids2arps'

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

    CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,            &
                  nzsnd,nsort,ncolp,                                    &
                  varfill,ivar,nyqset,timeset,rfropt,                   &
                  refchek,refmiss,bmwidth,refmedl,refdazl,iordref,      &
                  sortmin,dsort,                                        &
                  kntrgat,kntrazm,kntrelv,                              &
                  radlat,radlon,radarx,radary,radalt,dazim,             &
                  rngmin,rngmax,                                        &
                  rngrvol,azmrvol,elvrvol,                              &
                  refvol,timevolr,nyqvvol,rxvol,ryvol,rzvol,            &
                  zsnd,rfrsnd,kntbin,elvmnrvol,                         &
                  xcolp,ycolp,zcolp,                                    &
                  colref,coltim,colnyq,istatus)

    DO kelv=1,kntrelv
      refelvmin=min(elvmnrvol(kelv),refelvmin)
      refelvmax=max(elvmnrvol(kelv),refelvmax)
    END DO

    IF(myproc == ROOT) print *, ' outtime: ',outtime

    IF( ref3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,                    &
                      colref,icolp,jcolp,                               &
                      refid,refname,refunits,outtime,runname,           &
                      dirname,dmpfmt,hdfcompr,istatus)

    IF( ref2d ) THEN
      vardump=1
      ivar=1
      varid='refl2d'
      varname='Low-level reflect'
      varunits='dBZ'
      curtim=outtime
      IF (myproc == ROOT) print *, ' Calling remap2d for reflectivity'
      CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,           &
                   vardump,ivar,rfropt,varid,varname,varunits,dmpfmt,hdf4cmpr, &
                   refchek,refmiss,refmedl,refdazl,iordref,             &
                   sortmin,dsort,                                       &
                   kntrgat,kntrazm,kntrelv,                             &
                   radlat,radlon,radarx,radary,radalt,dazim,            &
                   rngmin,rngmax,rngrvol,azmrvol,elvrvol,               &
                   refvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,         &
                   tem2d,istatus)
    END IF

    IF( velproc ) THEN
      nyqset=1
      timeset=1
      IF(vel3d) THEN
        vardump=1
      ELSE
        vardump=0
      END IF
      vardump=0
      varfill=0
      ivar=2

      CALL rgatexyz(maxgate,maxgate,maxazim,maxelev,nzsnd,              &
                    rfropt,lvldbg,                                      &
                    kntvgat,kntvazm,kntvelv,                            &
                    radlat,radlon,radarx,radary,radalt,                 &
                    rngvvol,azmvvol,elvvvol,                            &
                    zsnd,rfrsnd,                                        &
                    rxvol,ryvol,rzvol,istatus)

      IF (myproc == ROOT) WRITE(stdout,'(a)') ' Calling remap3dcol for velocity '
      CALL remap3dcol(maxgate,maxgate,maxazim,maxelev,nxny,nz,          &
                  nzsnd,nsort,ncolp,                                    &
                  varfill,ivar,nyqset,timeset,rfropt,                   &
                  velchek,velmiss,bmwidth,velmedl,veldazl,iordvel,      &
                  sortmin,dsort,                                        &
                  kntvgat,kntvazm,kntvelv,                              &
                  radlat,radlon,radarx,radary,radalt,dazim,             &
                  rngmin,rngmax,                                        &
                  rngvvol,azmvvol,elvvvol,                              &
                  velvol,timevolv,nyqvvol,rxvol,ryvol,rzvol,            &
                  zsnd,rfrsnd,kntbin,elvmnvvol,                         &
                  xcolp,ycolp,zcolp,                                    &
                  colvel,coltim,colnyq,istatus)

      IF( vel3d ) CALL wrtcol2grd(nx,ny,nz,nxny,ncolp,                  &
                      colvel,icolp,jcolp,                               &
                      velid,velname,velunits,outtime,runname,           &
                      dirname,dmpfmt,hdfcompr,istatus)

    END IF

    IF( vel2d ) THEN
      vardump=1
      ivar=2
      varid='radv2d'
      varname='Low-level Velocity'
      varunits='m/s'
      curtim=outtime
      IF (myproc == ROOT) print *, ' Calling remap2d for velocity '
      CALL remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,           &
                   vardump,ivar,rfropt,varid,varname,varunits,dmpfmt,hdf4cmpr, &
                   velchek,velmiss,velmedl,veldazl,iordvel,             &
                   sortmin,dsort,                                       &
                   kntvgat,kntvazm,kntvelv,                             &
                   radlat,radlon,radarx,radary,radalt,dazim,            &
                   rngmin,rngmax,rngvvol,azmvvol,elvvvol,               &
                   velvol,rxvol,ryvol,xs,ys,zsnd,rfrsnd,kntbin,        &
                   tem2d,istatus)
    END IF

!------------------------------------------------------------------------
!
! Now use NIDS echo top and VIL fields to adjust remapped reflectivity.
!
!------------------------------------------------------------------------

! Developed for WDT.
!   IF (adjreflopt == 1) THEN

!     IF (myproc == ROOT) WRITE(stdout,'(a)') ' Adjusting reflectivity field...'

!     CALL adjust_refl(nx,ny,nz,colref,zcol,nidset,nidsvil)

!   END IF

!------------------------------------------------------------------------
!
!   Create filename and write file.
!
!   Force seconds to be 0, as mkradfn will round up, and we can't anticipate
!   that ahead of time.
!
!------------------------------------------------------------------------

    CALL mkradfnm(dmpfmt,dir_name,len_dir,radar_name,iiyear,iimonth,iiday, &
                  iihr, iimin, 0,     full_fname, len_fname)
!                 iihr, iimin, iisec, full_fname, len_fname)

    IF (myproc == ROOT) WRITE(stdout,'(/,1x,2a,/)')                     &
                        'Filename for this volume: ',TRIM(full_fname)

    iradfmt=1
    dualpol=0
    CALL wrtradcol(nx,ny,nxny,nz,ncolp,ncoltot,                         &
                  dmpfmt,iradfmt,hdf4cmpr,dmpzero,dualpol,              &
                  full_fname,radname,radlat,radlon,radalt,              &
                  iiyear,iimonth,iiday,iihr,iimin,iisec,ivcp,isource,   &
                  refelvmin,refelvmax,refrngmin,refrngmax,              &
                  icolp,jcolp,xcolp,ycolp,zcolp,                        &
                  colref,colvel,colnyq,coltim,colrhv,colzdr,colkdp,     &
                  istatus)

  END IF

!------------------------------------------------------------------------
!
! The End.
!
!------------------------------------------------------------------------

  IF (myproc == ROOT)  &
    WRITE(stdout,'(/a/)') '  === Normal termination of NIDS2ARPS === '

  CALL mpexit(0)

  STOP
END PROGRAM nids2arps
