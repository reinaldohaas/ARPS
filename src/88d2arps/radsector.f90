  PROGRAM RADSECTOR
!
!
! Reads a NetCDF file containing radar data and writes the file
! containing only those radials in a sector defined by azimbgn and azimend,
! where azimend is clockwise from azimbgn.  Both are in degrees from North.
!
! The transformation is applied to both the Reflectivity, Velocity and Vorticity files.
! They are presumed to be in the same directory.  Default is the present
! working directory (./), or specified following the -indir flag.  The output 
! files have the same filenames but are stored in a new directory, specified by
! the -outdir flag.
!
! Usage
!
! radsector Filename [-indir input_directory] -outdir output_directory 
!                        -azimbgn beginning_azimuth -azimend ending_azimuth
!
! AUTHOR:
! Keith Brewster, CAPS
! 04/08/2005
!
! MODIFICATIONS:
! 04/12/2005 (Keith Brewster, CAPS)
! Added index file writing.
!
! 05/25/2005 (Keith Brewster, CAPS)
! Added processing of unattenuated reflectivity and velocity.
!
  IMPLICIT NONE
!
! Dimensions
! 
  INTEGER :: nazim,nazimin
  INTEGER :: ngate
!
! Radar data variables
!
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=256) :: rfname
  CHARACTER (LEN=256) :: vfname
  CHARACTER (LEN=256) :: uarfname
  CHARACTER (LEN=256) :: uavfname
  CHARACTER (LEN=256) :: hfname
  CHARACTER (LEN=256) :: idxfname
  CHARACTER (LEN=40)  :: varname
  CHARACTER (LEN=80)  :: indir
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=120) :: cmd
  CHARACTER (LEN=256) :: cdfname
  CHARACTER (LEN=4) :: radname
  CHARACTER (LEN=5) :: cmprext

  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate
  LOGICAL :: fexist
!
  REAL,ALLOCATABLE :: azimin(:)
  REAL,ALLOCATABLE :: beamwin(:)
  REAL,ALLOCATABLE :: gtspcin(:)
  REAL,ALLOCATABLE :: vnyqin(:)
  REAL,ALLOCATABLE :: reflin(:,:)
  REAL,ALLOCATABLE :: radvin(:,:)
  REAL,ALLOCATABLE :: vortin(:,:)
  INTEGER,ALLOCATABLE :: refqc(:,:)
  INTEGER,ALLOCATABLE :: vrqc(:,:)
!
  REAL,ALLOCATABLE :: azim(:)
  REAL,ALLOCATABLE :: beamw(:)
  REAL,ALLOCATABLE :: gtspc(:)
  REAL,ALLOCATABLE :: vnyq(:)
  REAL,ALLOCATABLE :: refl(:,:)
  REAL,ALLOCATABLE :: radv(:,:)
  REAL,ALLOCATABLE :: vort(:,:)
!
!-----------------------------------------------------------------------
!
! NetCDF variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ncmode
  INTEGER :: ncid,nazimid,ngateid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: itim1970=315619200

  CHARACTER(LEN=256) :: argstr
  CHARACTER(LEN=256) :: tmpstr
  REAL :: azimbgn,azimend
  REAL :: dazmax,dazmin,dazim
  INTEGER :: iarg,narg,iloc,ltmp,lfname
  INTEGER :: igate,ielv,felv,jazim,kazim,jstart,irot
  INTEGER :: istatus,iostatus
  INTEGER :: iargc,creidx,openidx,ifsecs,itime,idxunit
  LOGICAL :: back,needfname
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
! Variable initializations
!
  openidx=0
  creidx=1
  azimbgn=0.
  azimend=30.
  indir='.'
  outdir='radsect'
  fname='dummy'
  needfname=.true.
  cmprext='.gz'
!
!-----------------------------------------------------------------------
!
!  Get input arguments
!
!-----------------------------------------------------------------------
!
  narg=iargc()
  IF(narg < 7) THEN
    WRITE(6,'(a,a/a)') &
      ' Usage: radsector Filename [-indir input_directory]', &
      ' -outdir output_directory ', &
      ' -azimbgn beginning_azimuth -azimend ending_azimuth'
    STOP
  END IF
  iarg=1
  DO
    CALL getarg(iarg,argstr)
    iarg=iarg+1
    IF(argstr(1:8) == '-azimbgn') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      READ(argstr,*,iostat=iostatus) azimbgn
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a/a)') &
         ' Usage: radsector Filename [-indir input_directory]', &
         ' -outdir output_directory ', &
         ' -azimbgn beginning_azimuth -azimend ending_azimuth'
        STOP
      END IF
    ELSE IF(argstr(1:8) == '-azimend') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      READ(argstr,*,iostat=iostatus) azimend
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a/a)') ' radsector Filename [-indir input_directory]', &
         ' -outdir output_directory ', &
         ' -azimbgn beginning_azimuth -azimend ending_azimuth'
        STOP
      END IF
    ELSE IF (argstr(1:6) == '-indir') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      indir=argstr(1:80)
    ELSE IF (argstr(1:7) == '-outdir') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      outdir=argstr(1:80)
    ELSE IF (argstr(1:1) /= '-' .AND. needfname) THEN
      fname=argstr
      needfname=.false.
    ELSE
      WRITE(6,'(a,a)') &
         ' Unrecognized command line input: ',TRIM(argstr)
      WRITE(6,'(a,a/a)') ' radsector Filename [-indir input_directory]', &
         ' -outdir output_directory ', &
         ' -azimbgn beginning_azimuth -azimend ending_azimuth'
      STOP
    END IF
    IF(iarg > narg) EXIT
  END DO
!
!-----------------------------------------------------------------------
!
!  Strip directory name, if any from filename
!
!-----------------------------------------------------------------------
!
  back=.true.
  iloc=INDEX(fname,'/',back)
  ltmp=LEN(tmpstr)
  lfname=LEN(TRIM(fname))
  IF(iloc > 0) THEN
    indir=fname(1:(iloc-1))
    tmpstr=TRIM(fname)
    fname=tmpstr((iloc+1):ltmp)
  END IF
  IF(fname(1:7) == 'Reflect') THEN
    rfname=fname
    WRITE(vfname,'(a,a)') 'Velocity',fname(13:lfname) 
    WRITE(uavfname,'(a,a)') 'VelocityUnattenuated',fname(13:lfname) 
    WRITE(uarfname,'(a,a)') 'ReflectivityUnattenuated',fname(13:lfname) 
    WRITE(hfname,'(a,a)') 'Vorticity',fname(13:lfname) 
  ELSE IF(fname(1:8) == 'Velocity') THEN
    vfname=fname
    WRITE(rfname,'(a,a)') 'Reflectivity',fname(9:lfname) 
    WRITE(uarfname,'(a,a)') 'ReflectivityUnattenuated',fname(9:lfname) 
    WRITE(uavfname,'(a,a)') 'VelocityUnattenuated',fname(9:lfname) 
    WRITE(hfname,'(a,a)') 'Vorticity',fname(9:lfname) 
  ELSE IF(fname(1:9) == 'Vorticity') THEN
    hfname=fname
    WRITE(rfname,'(a,a)') 'Reflectivity',fname(10:lfname) 
    WRITE(vfname,'(a,a)') 'Velocity',fname(10:lfname) 
    WRITE(uarfname,'(a,a)') 'ReflectivityUnattenuated',fname(10:lfname) 
    WRITE(uavfname,'(a,a)') 'VelocityUnattenuated',fname(10:lfname) 
  END IF
  WRITE(6,'(a,a)') ' Input directory: ',TRIM(indir)
  WRITE(6,'(a,a)') ' Ouput directory: ',TRIM(outdir)
  WRITE(6,'(a,a)') ' Reflectivity Filename: ',TRIM(rfname)
  WRITE(6,'(a,a)') ' Velocity Filename: ',TRIM(vfname)
  WRITE(6,'(a,a)') ' Unattenuated Reflectivity Filename: ',TRIM(uarfname)
  WRITE(6,'(a,a)') ' Unattenuated Velocity Filename: ',TRIM(uavfname)
  WRITE(6,'(a,a)') ' Vorticity Filename: ',TRIM(hfname)
  WRITE(6,'(a,f6.2)') ' Sector beginning azimuth: ',azimbgn
  WRITE(6,'(a,f6.2//)') ' Sector ending azimuth: ',azimend

  IF(TRIM(outdir) /= '.') THEN
    INQUIRE(FILE=TRIM(outdir),EXIST=fexist)
    IF(.NOT.fexist) THEN
      WRITE(6,'(a,a)') ' Creating directory ',TRIM(outdir)
      WRITE(cmd,'(a,a)') 'mkdir ',TRIM(outdir)
      CALL unixcmd(cmd)
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Open and get reflectivity dimensions from netCDF file
!
!-----------------------------------------------------------------------
!
  nazimin=0
  ngate=0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(rfname)
  write(6,'(a,a)') ' Reading data dimensions from: ',TRIM(cdfname)

  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF(istatus == 0) THEN
  istatus=nf_inq_dimid(ncid,'Azimuth',nazimid)
  istatus=nf_inq_dimid(ncid,'Gate',ngateid)
  istatus=nf_inq_dimlen(ncid,nazimid,nazimin)
  istatus=nf_inq_dimlen(ncid,ngateid,ngate)
  istatus=nf_close(ncid)

  write(6,'(a,i5,a,i5)') ' nazimin=',nazimin,'  ngate=',ngate
  IF( nazimin == 0 .OR. ngate == 0) THEN
    WRITE(6,'(a)') ' Error getting data dimensions.'
    STOP
  END IF

  ALLOCATE(azimin(nazimin))
  ALLOCATE(beamwin(nazimin))
  ALLOCATE(gtspcin(nazimin))
  ALLOCATE(reflin(ngate,nazimin))
  ALLOCATE(refqc(ngate,nazimin))
  azimin=-199.
  beamwin=-199.
  gtspcin=-199.
  reflin=-199.
  runname='dummy'

  CALL rdrftiltcdf(nazimin,ngate,rfname,indir,radname,                  &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azimin,beamwin,gtspcin,reflin,refqc)

  ielv=INT(elv)
  felv=NINT(100.*(elv-ielv))
  CALL gtlfnkey(runname,lfnkey)
  itime=itimcdf+itim1970
  CALL abss2ctim(itime,year,month,day,hour,minute,second)

  IF(azimend > azimbgn) THEN
    dazmax=azimend-azimbgn
  ELSE
    dazmax=360.+azimend-azimbgn
  END IF
  WRITE(6,'(a,f7.2)') ' Sector width (degrees): ',dazmax
!
! Count radials in input arrays within specified sector
!
  nazim=0
  DO jazim=1,nazimin
    dazim=azimin(jazim)-azimbgn
    IF(dazim < 0.) dazim=dazim+360.
    IF(dazim <= dazmax) nazim=nazim+1
  END DO
  WRITE(6,'(i6,a,f9.2,a,f9.2)') &
    nazim,' radials found in range ',azimbgn,' to ',azimend
!
! Find rotating direction and starting azimuth
! 
  irot=1
  dazim=azimin(2)-azimin(1)
  IF(dazim < -180.) dazim=dazim+360.
  IF(dazim < 0.) irot=-1
  IF( irot > 0 ) THEN
    WRITE(6,'(a)') ' Clockwise rotation detected'
  ELSE
    WRITE(6,'(a)') ' Counter-clockwise rotation detected'
  END IF

  jstart=1
  dazmin=361.
  IF(irot > 0) THEN
    DO jazim=1,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  ELSE
    DO jazim=1,nazimin
      dazim=azimend-azimin(jazim)
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  END IF
!
  IF(nazim > 0) THEN
!
!-----------------------------------------------------------------------
!
! Open index record file
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      ifsecs=NINT(curtim)
      write(idxfname,'(a,a,a,a,a,a,i6.6,a)') &
      TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
      idxunit=32
      write(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
      OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
      write(idxunit,'(a)') '<?xml version="1.0" encoding="iso-8859-1" ?>'
      write(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',      &
        TRIM(outdir),'">'
      openidx=1
    END IF

    ALLOCATE(azim(nazim))
    ALLOCATE(beamw(nazim))
    ALLOCATE(gtspc(nazim))
    ALLOCATE(refl(ngate,nazim))
    azim=-199.
    beamw=-199.
    gtspc=-199.
    refl=-199.

    kazim=0
    DO jazim=jstart,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        DO igate=1,ngate
          refl(igate,kazim)=reflin(igate,jazim)
        END DO
      END IF
    END DO 

    DO jazim=1,(jstart-1)
      IF(azimin(jazim) > azimend) THEN
        dazim=azimin(jazim)-azimbgn
      ELSE
        dazim=360.+azimin(jazim)-azimbgn
      END IF
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        DO igate=1,ngate
          refl(igate,kazim)=reflin(igate,jazim)
        END DO
      END IF
    END DO 

    varname='Reflectivity'
    CALL wtrftiltcdf(ngate,nazim,nazim,rfname,outdir,varname,          &
                       radname,radlat,radlon,radelv,ivcp,elv,          &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyquist,rfrgate,                               &
                       azim,beamw,gtspc,refl)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      write(idxunit,'(a)') '<item>'
      write(idxunit,'(a,f8.6,a,i10,a)') &
        '<time fractional="',frtime,'">',itimcdf,'</time>'
      write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',   &
         TRIM(rfname),TRIM(cmprext),'</params>'
      write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,a,i2.2,a1,i2.2,a)')        &
        '<selections>',year,month,day,'-',                            &
        hour,minute,second,' Reflectivity',ielv,'.',felv,             &
        '</selections>'
      write(idxunit,'(a)') '</item>'
    END IF
!
    DEALLOCATE(azimin)
    DEALLOCATE(beamwin)
    DEALLOCATE(gtspcin)
    DEALLOCATE(reflin)
    DEALLOCATE(refqc)
    DEALLOCATE(azim)
    DEALLOCATE(beamw)
    DEALLOCATE(gtspc)
    DEALLOCATE(refl)
  END IF
  ELSE
    WRITE(6,'(a,a)') &
      ' Error opening NetCDF data file: ',TRIM(cdfname)
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Open and get unattenuated reflectivity dimensions from netCDF file
!
!-----------------------------------------------------------------------
!
  nazimin=0
  ngate=0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(uarfname)
  write(6,'(a,a)') ' Reading data dimensions from: ',TRIM(cdfname)

  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF(istatus == 0) THEN
  istatus=nf_inq_dimid(ncid,'Azimuth',nazimid)
  istatus=nf_inq_dimid(ncid,'Gate',ngateid)
  istatus=nf_inq_dimlen(ncid,nazimid,nazimin)
  istatus=nf_inq_dimlen(ncid,ngateid,ngate)
  istatus=nf_close(ncid)

  write(6,'(a,i5,a,i5)') ' nazimin=',nazimin,'  ngate=',ngate
  IF( nazimin == 0 .OR. ngate == 0) THEN
    WRITE(6,'(a)') ' Error getting data dimensions.'
    STOP
  END IF

  ALLOCATE(azimin(nazimin))
  ALLOCATE(beamwin(nazimin))
  ALLOCATE(gtspcin(nazimin))
  ALLOCATE(reflin(ngate,nazimin))
  ALLOCATE(refqc(ngate,nazimin))
  azimin=-199.
  beamwin=-199.
  gtspcin=-199.
  reflin=-199.
  runname='dummy'

  CALL rdrftiltcdf(nazimin,ngate,rfname,indir,radname,                  &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azimin,beamwin,gtspcin,reflin,refqc)

  ielv=INT(elv)
  felv=NINT(100.*(elv-ielv))
  CALL gtlfnkey(runname,lfnkey)
  itime=itimcdf+itim1970
  CALL abss2ctim(itime,year,month,day,hour,minute,second)

  IF(azimend > azimbgn) THEN
    dazmax=azimend-azimbgn
  ELSE
    dazmax=360.+azimend-azimbgn
  END IF
  WRITE(6,'(a,f7.2)') ' Sector width (degrees): ',dazmax
!
! Count radials in input arrays within specified sector
!
  nazim=0
  DO jazim=1,nazimin
    dazim=azimin(jazim)-azimbgn
    IF(dazim < 0.) dazim=dazim+360.
    IF(dazim <= dazmax) nazim=nazim+1
  END DO
  WRITE(6,'(i6,a,f9.2,a,f9.2)') &
    nazim,' radials found in range ',azimbgn,' to ',azimend
!
! Find rotating direction and starting azimuth
! 
  irot=1
  dazim=azimin(2)-azimin(1)
  IF(dazim < -180.) dazim=dazim+360.
  IF(dazim < 0.) irot=-1
  IF( irot > 0 ) THEN
    WRITE(6,'(a)') ' Clockwise rotation detected'
  ELSE
    WRITE(6,'(a)') ' Counter-clockwise rotation detected'
  END IF

  jstart=1
  dazmin=361.
  IF(irot > 0) THEN
    DO jazim=1,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  ELSE
    DO jazim=1,nazimin
      dazim=azimend-azimin(jazim)
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  END IF
!
  IF(nazim > 0) THEN
!
!-----------------------------------------------------------------------
!
! Open index record file
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      ifsecs=NINT(curtim)
      write(idxfname,'(a,a,a,a,a,a,i6.6,a)') &
      TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
      idxunit=32
      write(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
      OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
      write(idxunit,'(a)') '<?xml version="1.0" encoding="iso-8859-1" ?>'
      write(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',      &
        TRIM(outdir),'">'
      openidx=1
    END IF

    ALLOCATE(azim(nazim))
    ALLOCATE(beamw(nazim))
    ALLOCATE(gtspc(nazim))
    ALLOCATE(refl(ngate,nazim))
    azim=-199.
    beamw=-199.
    gtspc=-199.
    refl=-199.

    kazim=0
    DO jazim=jstart,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        DO igate=1,ngate
          refl(igate,kazim)=reflin(igate,jazim)
        END DO
      END IF
    END DO 

    DO jazim=1,(jstart-1)
      IF(azimin(jazim) > azimend) THEN
        dazim=azimin(jazim)-azimbgn
      ELSE
        dazim=360.+azimin(jazim)-azimbgn
      END IF
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        DO igate=1,ngate
          refl(igate,kazim)=reflin(igate,jazim)
        END DO
      END IF
    END DO 

    varname='UnattenuatedReflectivity'
    CALL wtrftiltcdf(ngate,nazim,nazim,rfname,outdir,                  &
                       radname,radlat,radlon,radelv,ivcp,elv,          &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyquist,rfrgate,                               &
                       azim,beamw,gtspc,refl)

!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      write(idxunit,'(a)') '<item>'
      write(idxunit,'(a,f8.6,a,i10,a)') &
        '<time fractional="',frtime,'">',itimcdf,'</time>'
      write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',   &
         TRIM(rfname),TRIM(cmprext),'</params>'
      write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,a,i2.2,a1,i2.2,a)')        &
        '<selections>',year,month,day,'-',                            &
        hour,minute,second,' UnattenuatedReflectivity',ielv,'.',felv, &
        '</selections>'
      write(idxunit,'(a)') '</item>'
    END IF
!
    DEALLOCATE(azimin)
    DEALLOCATE(beamwin)
    DEALLOCATE(gtspcin)
    DEALLOCATE(reflin)
    DEALLOCATE(refqc)
    DEALLOCATE(azim)
    DEALLOCATE(beamw)
    DEALLOCATE(gtspc)
    DEALLOCATE(refl)
  END IF
  ELSE
    WRITE(6,'(a,a)') &
      ' Error opening NetCDF data file: ',TRIM(cdfname)
  END IF
!-----------------------------------------------------------------------
!
!  Open and get velocity dimensions from netCDF file
!
!-----------------------------------------------------------------------
!
  nazimin=0
  ngate=0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(vfname)
  write(6,'(a,a)') ' Reading data dimensions from: ',TRIM(cdfname)

  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF(istatus == 0) THEN
  istatus=nf_inq_dimid(ncid,'Azimuth',nazimid)
  istatus=nf_inq_dimid(ncid,'Gate',ngateid)
  istatus=nf_inq_dimlen(ncid,nazimid,nazimin)
  istatus=nf_inq_dimlen(ncid,ngateid,ngate)
  istatus=nf_close(ncid)

  write(6,'(a,i5,a,i5)') ' nazimin=',nazimin,'  ngate=',ngate
  IF( nazimin == 0 .OR. ngate == 0) THEN
    WRITE(6,'(a)') ' Error getting data dimensions.'
    STOP
  END IF

  ALLOCATE(azimin(nazimin))
  ALLOCATE(beamwin(nazimin))
  ALLOCATE(gtspcin(nazimin))
  ALLOCATE(vnyqin(nazimin))
  ALLOCATE(radvin(ngate,nazimin))
  ALLOCATE(vrqc(ngate,nazimin))
  azimin=-199.
  beamwin=-199.
  gtspcin=-199.
  vnyqin=-199.
  radvin=-199.
  runname='dummy'

  CALL rdvrtiltcdf(nazimin,ngate,vfname,indir,radname,                  &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azimin,beamwin,gtspcin,vnyqin,radvin,vrqc)

  ielv=INT(elv)
  felv=NINT(100.*(elv-ielv))
  CALL gtlfnkey(runname,lfnkey)

  IF(azimend > azimbgn) THEN
    dazmax=azimend-azimbgn
  ELSE
    dazmax=360.+azimend-azimbgn
  END IF
  WRITE(6,'(a,f7.2)') ' Sector width (degrees): ',dazmax
!
! Count radials in input arrays within specified sector
!
  nazim=0
  DO jazim=1,nazimin
    dazim=azimin(jazim)-azimbgn
    IF(dazim < 0.) dazim=dazim+360.
    IF(dazim <= dazmax) nazim=nazim+1
  END DO
  WRITE(6,'(i6,a,f9.2,a,f9.2)') &
    nazim,' radials found in range ',azimbgn,' to ',azimend
!
! Find rotating direction and starting azimuth
! 
  irot=1
  dazim=azimin(2)-azimin(1)
  IF(dazim < -180.) dazim=dazim+360.
  IF(dazim < 0.) irot=-1
  IF( irot > 0 ) THEN
    WRITE(6,'(a)') ' Clockwise rotation detected'
  ELSE
    WRITE(6,'(a)') ' Counter-clockwise rotation detected'
  END IF

  jstart=1
  dazmin=361.
  IF(irot > 0) THEN
    DO jazim=1,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  ELSE
    DO jazim=1,nazimin
      dazim=azimend-azimin(jazim)
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  END IF
!
  IF(nazim > 0) THEN
  
    IF( creidx > 0  .and. openidx < 1) THEN
      ifsecs=NINT(curtim)
      write(idxfname,'(a,a,a,a,a,a,i6.6,a)') &
      TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
      idxunit=32
      write(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
      OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
      write(idxunit,'(a)') '<?xml version="1.0" encoding="iso-8859-1" ?>'
      write(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',      &
        TRIM(outdir),'">'
      openidx = 1
    END IF

    ALLOCATE(azim(nazim))
    ALLOCATE(beamw(nazim))
    ALLOCATE(gtspc(nazim))
    ALLOCATE(vnyq(nazim))
    ALLOCATE(radv(ngate,nazim))
    azim=-199.
    beamw=-199.
    gtspc=-199.
    vnyq=-199.
    radv=-199.

    kazim=0
    DO jazim=jstart,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        vnyq(kazim)=vnyqin(jazim)
        DO igate=1,ngate
          radv(igate,kazim)=radvin(igate,jazim)
        END DO
      END IF
    END DO 

    DO jazim=1,(jstart-1)
      IF(azimin(jazim) > azimend) THEN
        dazim=azimin(jazim)-azimbgn
      ELSE
        dazim=360.+azimin(jazim)-azimbgn
      END IF
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        vnyq(kazim)=vnyqin(jazim)
        DO igate=1,ngate
          radv(igate,kazim)=radvin(igate,jazim)
        END DO
      END IF
    END DO 

    varname='Velocity'
    CALL wtvrtiltcdf(ngate,nazim,nazim,vfname,outdir,varname,          &
                     radname,radlat,radlon,radelv,ivcp,elv,            &
                     rmisval,rngfval,itimcdf,frtime,initime,           &
                     vnyquist,rfrgate,                                 &
                     azim,beamw,gtspc,vnyq,radv)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      write(idxunit,'(a)') '<item>'
      write(idxunit,'(a,f8.6,a,i10,a)') &
       '<time fractional="',frtime,'">',itimcdf,'</time>'
      write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',   &
         TRIM(vfname),TRIM(cmprext),'</params>'
      write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,a,i2.2,a1,i2.2,a)')        &
        '<selections>',year,month,day,'-',                            &
        hour,minute,second,' Velocity',ielv,'.',felv,                 &
        '</selections>'

      write(idxunit,'(a)') '</item>'
    END IF

    DEALLOCATE(azimin)
    DEALLOCATE(beamwin)
    DEALLOCATE(gtspcin)
    DEALLOCATE(vnyqin)
    DEALLOCATE(radvin)
    DEALLOCATE(vrqc)
    DEALLOCATE(azim)
    DEALLOCATE(beamw)
    DEALLOCATE(gtspc)
    DEALLOCATE(vnyq)
    DEALLOCATE(radv)
  END IF
  ELSE
    WRITE(6,'(a,a)') &
      ' Error opening NetCDF data file: ',TRIM(cdfname)
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Open and get unattenuated velocity dimensions from netCDF file
!
!-----------------------------------------------------------------------
!
  nazimin=0
  ngate=0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(uavfname)
  write(6,'(a,a)') ' Reading data dimensions from: ',TRIM(cdfname)

  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF(istatus == 0) THEN
  istatus=nf_inq_dimid(ncid,'Azimuth',nazimid)
  istatus=nf_inq_dimid(ncid,'Gate',ngateid)
  istatus=nf_inq_dimlen(ncid,nazimid,nazimin)
  istatus=nf_inq_dimlen(ncid,ngateid,ngate)
  istatus=nf_close(ncid)

  write(6,'(a,i5,a,i5)') ' nazimin=',nazimin,'  ngate=',ngate
  IF( nazimin == 0 .OR. ngate == 0) THEN
    WRITE(6,'(a)') ' Error getting data dimensions.'
    STOP
  END IF

  ALLOCATE(azimin(nazimin))
  ALLOCATE(beamwin(nazimin))
  ALLOCATE(gtspcin(nazimin))
  ALLOCATE(vnyqin(nazimin))
  ALLOCATE(radvin(ngate,nazimin))
  ALLOCATE(vrqc(ngate,nazimin))
  azimin=-199.
  beamwin=-199.
  gtspcin=-199.
  vnyqin=-199.
  radvin=-199.
  runname='dummy'

  CALL rdvrtiltcdf(nazimin,ngate,vfname,indir,radname,                  &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azimin,beamwin,gtspcin,vnyqin,radvin,vrqc)

  ielv=INT(elv)
  felv=NINT(100.*(elv-ielv))
  CALL gtlfnkey(runname,lfnkey)

  IF(azimend > azimbgn) THEN
    dazmax=azimend-azimbgn
  ELSE
    dazmax=360.+azimend-azimbgn
  END IF
  WRITE(6,'(a,f7.2)') ' Sector width (degrees): ',dazmax
!
! Count radials in input arrays within specified sector
!
  nazim=0
  DO jazim=1,nazimin
    dazim=azimin(jazim)-azimbgn
    IF(dazim < 0.) dazim=dazim+360.
    IF(dazim <= dazmax) nazim=nazim+1
  END DO
  WRITE(6,'(i6,a,f9.2,a,f9.2)') &
    nazim,' radials found in range ',azimbgn,' to ',azimend
!
! Find rotating direction and starting azimuth
! 
  irot=1
  dazim=azimin(2)-azimin(1)
  IF(dazim < -180.) dazim=dazim+360.
  IF(dazim < 0.) irot=-1
  IF( irot > 0 ) THEN
    WRITE(6,'(a)') ' Clockwise rotation detected'
  ELSE
    WRITE(6,'(a)') ' Counter-clockwise rotation detected'
  END IF

  jstart=1
  dazmin=361.
  IF(irot > 0) THEN
    DO jazim=1,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  ELSE
    DO jazim=1,nazimin
      dazim=azimend-azimin(jazim)
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  END IF
!
  IF(nazim > 0) THEN
  
    IF( creidx > 0  .and. openidx < 1) THEN
      ifsecs=NINT(curtim)
      write(idxfname,'(a,a,a,a,a,a,i6.6,a)') &
      TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
      idxunit=32
      write(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
      OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
      write(idxunit,'(a)') '<?xml version="1.0" encoding="iso-8859-1" ?>'
      write(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',      &
        TRIM(outdir),'">'
      openidx = 1
    END IF

    ALLOCATE(azim(nazim))
    ALLOCATE(beamw(nazim))
    ALLOCATE(gtspc(nazim))
    ALLOCATE(vnyq(nazim))
    ALLOCATE(radv(ngate,nazim))
    azim=-199.
    beamw=-199.
    gtspc=-199.
    vnyq=-199.
    radv=-199.

    kazim=0
    DO jazim=jstart,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        vnyq(kazim)=vnyqin(jazim)
        DO igate=1,ngate
          radv(igate,kazim)=radvin(igate,jazim)
        END DO
      END IF
    END DO 

    DO jazim=1,(jstart-1)
      IF(azimin(jazim) > azimend) THEN
        dazim=azimin(jazim)-azimbgn
      ELSE
        dazim=360.+azimin(jazim)-azimbgn
      END IF
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        vnyq(kazim)=vnyqin(jazim)
        DO igate=1,ngate
          radv(igate,kazim)=radvin(igate,jazim)
        END DO
      END IF
    END DO 

    varname='UnattenuatedVelocity'
    CALL wtvrtiltcdf(ngate,nazim,nazim,vfname,outdir,varname,          &
                     radname,radlat,radlon,radelv,ivcp,elv,            &
                     rmisval,rngfval,itimcdf,frtime,initime,           &
                     vnyquist,rfrgate,                                 &
                     azim,beamw,gtspc,vnyq,radv)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      write(idxunit,'(a)') '<item>'
      write(idxunit,'(a,f8.6,a,i10,a)') &
       '<time fractional="',frtime,'">',itimcdf,'</time>'
      write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',   &
         TRIM(vfname),TRIM(cmprext),'</params>'
      write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,a,i2.2,a1,i2.2,a)')        &
        '<selections>',year,month,day,'-',                            &
        hour,minute,second,' UnattenuatedVelocity',ielv,'.',felv,     &
        '</selections>'
      write(idxunit,'(a)') '</item>'
    END IF

    DEALLOCATE(azimin)
    DEALLOCATE(beamwin)
    DEALLOCATE(gtspcin)
    DEALLOCATE(vnyqin)
    DEALLOCATE(radvin)
    DEALLOCATE(vrqc)
    DEALLOCATE(azim)
    DEALLOCATE(beamw)
    DEALLOCATE(gtspc)
    DEALLOCATE(vnyq)
    DEALLOCATE(radv)
  END IF
  ELSE 
    WRITE(6,'(a,a)') &
      ' Error opening NetCDF data file: ',TRIM(cdfname)
  END IF
!-----------------------------------------------------------------------
!
!  Open and get vorticity dimensions from netCDF file
!
!-----------------------------------------------------------------------
!
  nazimin=0
  ngate=0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(hfname)
  write(6,'(a,a)') ' Reading data dimensions from: ',TRIM(cdfname)

  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF(istatus == 0) THEN
  istatus=nf_inq_dimid(ncid,'Azimuth',nazimid)
  istatus=nf_inq_dimid(ncid,'Gate',ngateid)
  istatus=nf_inq_dimlen(ncid,nazimid,nazimin)
  istatus=nf_inq_dimlen(ncid,ngateid,ngate)
  istatus=nf_close(ncid)

  write(6,'(a,i5,a,i5)') ' nazimin=',nazimin,'  ngate=',ngate
  IF( nazimin == 0 .OR. ngate == 0) THEN
    WRITE(6,'(a)') ' Error getting data dimensions.'
    STOP
  END IF

  ALLOCATE(azimin(nazimin))
  ALLOCATE(beamwin(nazimin))
  ALLOCATE(gtspcin(nazimin))
  ALLOCATE(vortin(ngate,nazimin))
  azimin=-199.
  beamwin=-199.
  gtspcin=-199.
  vortin=-199.
  runname='dummy'

  CALL rdvvtiltcdf(nazimin,ngate,hfname,indir,radname,                  &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azimin,beamwin,gtspcin,vortin)

  ielv=INT(elv)
  felv=NINT(100.*(elv-ielv))
  CALL gtlfnkey(runname,lfnkey)

  IF(azimend > azimbgn) THEN
    dazmax=azimend-azimbgn
  ELSE
    dazmax=360.+azimend-azimbgn
  END IF
  WRITE(6,'(a,f7.2)') ' Sector width (degrees): ',dazmax
!
! Count radials in input arrays within specified sector
!
  nazim=0
  DO jazim=1,nazimin
    dazim=azimin(jazim)-azimbgn
    IF(dazim < 0.) dazim=dazim+360.
    IF(dazim <= dazmax) nazim=nazim+1
  END DO
  WRITE(6,'(i6,a,f9.2,a,f9.2)') &
    nazim,' radials found in range ',azimbgn,' to ',azimend
!
! Find rotating direction and starting azimuth
! 
  irot=1
  dazim=azimin(2)-azimin(1)
  IF(dazim < -180.) dazim=dazim+360.
  IF(dazim < 0.) irot=-1
  IF( irot > 0 ) THEN
    WRITE(6,'(a)') ' Clockwise rotation detected'
  ELSE
    WRITE(6,'(a)') ' Counter-clockwise rotation detected'
  END IF

  jstart=1
  dazmin=361.
  IF(irot > 0) THEN
    DO jazim=1,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  ELSE
    DO jazim=1,nazimin
      dazim=azimend-azimin(jazim)
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim < dazmin) THEN
        dazmin=dazim
        jstart=jazim
      END IF
    END DO
  END IF
!
  IF(nazim > 0) THEN
!
!-----------------------------------------------------------------------
!
! Open index record file
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0  .and. openidx < 1) THEN
      ifsecs=NINT(curtim)
      write(idxfname,'(a,a,a,a,a,a,i6.6,a)')                            &
      TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
      idxunit=32
      write(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
      OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
      write(idxunit,'(a)') '<?xml version="1.0" encoding="iso-8859-1" ?>'
      write(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',      &
        TRIM(outdir),'">'
      openidx = 1
    END IF

    ALLOCATE(azim(nazim))
    ALLOCATE(beamw(nazim))
    ALLOCATE(gtspc(nazim))
    ALLOCATE(vort(ngate,nazim))
    azim=-199.
    beamw=-199.
    gtspc=-199.
    vort=-199.

    kazim=0
    DO jazim=jstart,nazimin
      dazim=azimin(jazim)-azimbgn
      IF(dazim < 0.) dazim=dazim+360.
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        DO igate=1,ngate
          vort(igate,kazim)=vortin(igate,jazim)
        END DO
      END IF
    END DO 

    DO jazim=1,(jstart-1)
      IF(azimin(jazim) > azimend) THEN
        dazim=azimin(jazim)-azimbgn
      ELSE
        dazim=360.+azimin(jazim)-azimbgn
      END IF
      IF(dazim <= dazmax) THEN
        kazim=kazim+1
        IF(kazim > nazim) EXIT
        azim(kazim)=azimin(jazim)
        beamw(kazim)=beamwin(jazim)
        gtspc(kazim)=gtspcin(jazim)
        DO igate=1,ngate
          vort(igate,kazim)=vortin(igate,jazim)
        END DO
      END IF
    END DO 

    varname='Vorticity'
    CALL wtvvtiltcdf(ngate,nazim,nazim,hfname,outdir,varname,          &
                     radname,radlat,radlon,radelv,ivcp,elv,            &
                     rmisval,rngfval,itimcdf,frtime,initime,           &
                     vnyquist,rfrgate,                                 &
                     azim,beamw,gtspc,vort)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 ) THEN
      write(idxunit,'(a)') '<item>'
      write(idxunit,'(a,f8.6,a,i10,a)') &
        '<time fractional="',frtime,'">',itimcdf,'</time>'
      write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',   &
         TRIM(hfname),TRIM(cmprext),'</params>'
      write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,a,i2.2,a1,i2.2,a)')        &
        '<selections>',year,month,day,'-',                            &
        hour,minute,second,' Vorticity',ielv,'.',felv,                &
        '</selections>'
      write(idxunit,'(a)') '</item>'
    END IF
!
    DEALLOCATE(azimin)
    DEALLOCATE(beamwin)
    DEALLOCATE(gtspcin)
    DEALLOCATE(vortin)
    DEALLOCATE(azim)
    DEALLOCATE(beamw)
    DEALLOCATE(gtspc)
    DEALLOCATE(vort)
  END IF
  ELSE
    WRITE(6,'(a,a)') &
      ' Error opening NetCDF data file: ',TRIM(cdfname)
  END IF
!
! Finish index file and close it out
!
  IF( openidx > 0 ) THEN
    write(idxunit,'(a)') '</codeindex>'
    CLOSE(idxunit)
  END IF

  STOP

END PROGRAM RADSECTOR
