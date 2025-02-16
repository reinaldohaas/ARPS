  PROGRAM RADBIN2CDF
!
! Purpose:
!   Reads a volume of radar data in binary format and converts it
!   to NetCDF formatted files for use in WDSS-II.
!   Useful for dealing with OSSE files written in hdf which cannot
!   be directly written out in NetCDF due to the incompatibility of
!   HDF and NetCDF libraries
!
!   There is no input file, just two command-line arguments.
!
!   Usage:
!   radbin2cdf [-n nfiles -int interval] First_binary_file Output_directory
!              [-azimbgn beginning_azimuth -azimend ending_azimuth]
!
!   Keith Brewster, CAPS
!   15 November 2004
!
!   23 March 2005 (Keith Brewster, CAPS)
!   Added outer time loop for processing time series of volumes.
!
!   12 April 2005 (Keith Brewster, CAPS)
!   Added logic for sectorized output option.
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Radar tilt variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: azim(:)
  REAL, ALLOCATABLE :: beamw(:)
  REAL, ALLOCATABLE :: gtspc(:)
  REAL, ALLOCATABLE :: vnyq(:)
  REAL, ALLOCATABLE :: radv(:,:)
  REAL, ALLOCATABLE :: uaradv(:,:)
  REAL, ALLOCATABLE :: refl(:,:)
  REAL, ALLOCATABLE :: uarefl(:,:)
  REAL, ALLOCATABLE :: vort(:,:)
!
!-----------------------------------------------------------------------
!
!  Radar volume variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: itimvol(:)
  REAL, ALLOCATABLE :: vnyqvol(:)
  REAL, ALLOCATABLE :: rngvol(:,:)
  REAL, ALLOCATABLE :: azmvol(:,:)
  REAL, ALLOCATABLE :: elvvol(:,:)
  REAL, ALLOCATABLE :: refvol(:,:,:)
  REAL, ALLOCATABLE :: uarefvol(:,:,:)
  REAL, ALLOCATABLE :: velvol(:,:,:)
  REAL, ALLOCATABLE :: uavelvol(:,:,:)
  REAL, ALLOCATABLE :: vorvol(:,:,:)
!
!-----------------------------------------------------------------------
!
!  I/O Options
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: argstr
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=256) :: idxfname
  CHARACTER (LEN=256) :: dirstr
  CHARACTER (LEN=256) :: rfnamevol
  CHARACTER (LEN=40 ) :: varname
  CHARACTER (LEN=4  ) :: radname
  INTEGER :: ifmt,creidx,ipktyp,nbits
  INTEGER :: wrtuaref,wrtuavel,wrtvort,idummy
  CHARACTER (LEN=80) :: outdir
  LOGICAL :: gotfiln,sectorize
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: cmprext
  CHARACTER (LEN=120) :: cmd

  REAL :: radelv,radlat,radlon
  REAL :: beamwid
  REAL :: elv,frtime
  REAL :: rmisval,rngfval,rfrgate,vnyqstl
  REAL :: azimbgn,azimend
  REAL :: dazmax,dazmin,dazim

  INTEGER, PARAMETER :: itim1970=315619200
  INTEGER :: iarg,narg,nfiles,ngatevel,ngateref
  INTEGER :: nazim,nazimin,nelev,ngate
  INTEGER :: ifile,itilt,jazim,igate,idxunit,iunit,intvl
  INTEGER :: ielv,felv,vcp,itime,itimcdf,ifsecs,initime
  INTEGER :: iyear,imon,iday,ihour,imin,isec
  INTEGER :: irot,jstart,kelv,kazim
  INTEGER :: iargc,abstsec,lfname,iostatus

  LOGICAL :: fexist,cmprsd
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
! A few initializations
!
!-----------------------------------------------------------------------
!
  outdir='.'
  frtime=0.
  cmprext='.gz'
  creidx=1
  rfnamevol='KTLX.dat'
  sectorize=.false.
  azimbgn=-999.
  azimend=-999.
  wrtuaref=0
  wrtuavel=0
  wrtvort =0
  idummy  =0
!
!-----------------------------------------------------------------------
!
! Open and read the binary file.
!
!-----------------------------------------------------------------------
!
  narg=iargc()
  IF(narg < 1) THEN
    WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]',   &
    ' first_binary_file [output_directory]', &
    ' [-azimbgn beginning_azimuth -azimend ending_azimuth]'
    STOP
  END IF
  nfiles=1
  intvl=-999
  rfnamevol='dummy'
  gotfiln=.false.
  iarg=1
  DO 
    CALL getarg(iarg,argstr)
    iarg=iarg+1
    IF(argstr(1:2) == '-n') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      READ(argstr,*,iostat=iostatus) nfiles
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]',&
        ' first_binary_file [output_directory]', &
        ' [-azimbgn beginning_azimuth -azimend ending_azimuth]'
        STOP
      END IF
    ELSE IF (argstr(1:4) == '-int') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      READ(argstr,*,iostat=iostatus) intvl
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]',&
        ' first_binary_file [output_directory]', &
        ' [-azimbgn beginning_azimuth -azimend ending_azimuth]'
        STOP
      END IF
    ELSE IF(argstr(1:8) == '-azimbgn') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      READ(argstr,*,iostat=iostatus) azimbgn
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]',&
        ' first_binary_file [output_directory]', &
         ' [-azimbgn beginning_azimuth -azimend ending_azimuth]'
        STOP
      END IF
      sectorize=.true.
    ELSE IF(argstr(1:8) == '-azimend') THEN
      IF(iarg > narg) EXIT
      CALL getarg(iarg,argstr)
      iarg=iarg+1
      READ(argstr,*,iostat=iostatus) azimend
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]',&
        ' first_binary_file [output_directory]', &
         ' -azimbgn beginning_azimuth -azimend ending_azimuth'
        STOP
      END IF
      sectorize=.true.
    ELSE IF (gotfiln) THEN
      outdir=argstr
    ELSE
      rfnamevol=argstr
      gotfiln=.true.
    END IF
    IF(iarg > narg) EXIT
  END DO
  CALL extrctdir(rfnamevol,dirstr)
  lfname=LEN_TRIM(rfnamevol)
  IF(rfnamevol((lfname-2):lfname) == '.gz') THEN
    WRITE(rfnamevol((lfname-2):lfname),'(a3)') '   '
  ELSE IF(rfnamevol((lfname-1):lfname) == '.Z') THEN
    WRITE(rfnamevol((lfname-1):lfname),'(a2)') '  '
  END IF
  WRITE(6,'(a,a)') ' Binary radar data: ',TRIM(rfnamevol)
  WRITE(6,'(a,a)') ' Input directory: ',TRIM(dirstr)
  WRITE(6,'(a,a)') ' Output directory: ',TRIM(outdir)
  WRITE(6,'(a,i4,a)') ' Processing',nfiles,' files'
  WRITE(6,'(a,i6,a)') ' at time interval ',intvl,' secs.'
  IF(nfiles > 1 .AND. intvl < 0) THEN
    WRITE(6,'(a)') ' nfiles > 1 and no time interval supplied'
    WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]', &
      ' first_binary_file [output_directory]', &
      ' -azimbgn beginning_azimuth -azimend ending_azimuth'
    STOP
  END IF

  IF(sectorize) THEN
    WRITE(6,'(a,f9.2,/,10x,a,f9.2)') &
       ' Sector beginning azimuth:',azimbgn, &
       ' Ending azimuth:',azimend
    IF(azimbgn < 0. .OR. azimend < 0.) THEN
      WRITE(6,'(a)') &
      ' Sector specified, must specify both azimbgn and azimend'
      WRITE(6,'(a,a/a)') 'Usage: radbin2cdf [-n nfiles -int interval]',&
      ' first_binary_file [output_directory]', &
      ' -azimbgn beginning_azimuth -azimend ending_azimuth'
      STOP
    END IF
    IF(azimend > azimbgn) THEN
      dazmax=azimend-azimbgn
    ELSE
      dazmax=360.+azimend-azimbgn
    END IF
    WRITE(6,'(a,f7.2)') ' Sector width (degrees): ',dazmax
  ELSE 
    WRITE(6,'(a)') ' No sectorization, writing all data.'
  END IF

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
! Open and read the binary file.
!
!-----------------------------------------------------------------------
!
  DO ifile=1,nfiles
    cmprsd=.false.
    IF(ifile > 1) THEN
      abstsec=abstsec+intvl
      CALL abss2ctim(abstsec,year,month,day,hour,minute,second)
      WRITE(rfnamevol,'(4a,i4.4,2(i2.2),a,3(i2.2),a)')                  &
        TRIM(dirstr),'/',radname,'_',year,month,day,'_',                &
        hour,minute,second,'.vol'
    END IF

    INQUIRE(file=TRIM(rfnamevol),EXIST=fexist)
    IF(.NOT.fexist) THEN
      INQUIRE(file=TRIM(rfnamevol)//'.gz',EXIST=fexist)
      IF(fexist) CALL uncmprs(TRIM(rfnamevol)//'.gz')
      cmprsd=fexist
    END IF
    IF(.NOT.fexist) THEN
      INQUIRE(file=TRIM(rfnamevol)//'.Z',EXIST=fexist)
      IF(fexist) CALL uncmprs(TRIM(rfnamevol)//'.Z')
      cmprsd=fexist
    END IF
    IF(.NOT.fexist) THEN
      WRITE(6,'(/1x,a,a,/a)')                                           &
        'File ',TRIM(rfnamevol),' or its compressed version not found.'
      STOP
    END IF

    iunit=31
    WRITE(6,'(a,a)') ' Opening: ',TRIM(rfnamevol)
    OPEN(iunit,FILE=trim(rfnamevol),STATUS='old',FORM='unformatted')
    READ(iunit) runname
    lfnkey=80
    CALL gtlfnkey(runname,lfnkey)
    WRITE(6,'(a,a)') ' Runname:',runname(1:lfnkey)
    READ(iunit) radname
    WRITE(6,'(a,a)') ' Radar name:',radname
    READ(iunit) ngatevel,ngateref,nazimin,nelev
    ngate=ngatevel
    WRITE(6,'(a,i6,a,i6,a,i6)')                                            &
      ' Ngate : ',ngate,' Nazim: ',nazimin,' Nelev: ',nelev
    READ(iunit) year,month,day,hour,minute,second
    WRITE(6,'(2(a,i2.2),a,i4.4,3(a,i2.2))')                                &
      ' Date: ',month,'/',day,'/',year,' Time: ',hour,':',minute,':',second
    READ(iunit) radelv,radlat,radlon
    WRITE(6,'(a,f9.2,a,f9.2,a,f9.2)')                                      &
      ' Latitude:',radlat,'  Longitude:',radlon,' Elevation:',radelv
    READ(iunit) vcp,beamwid,rmisval,rngfval,curtim 
    WRITE(6,'(a,i5,a,f9.2,a,f10.1,a,f10.1)')                               &
      ' VCP:',vcp,' Beamwidth:',beamwid,                                   &       
      ' MissVal: ',rmisval,' RngfVal :',rngfval
!
    CALL ctim2abss(year,month,day,hour,minute,second,abstsec)
    READ(iunit) wrtuaref,wrtuavel,wrtvort,                                 &
                idummy,idummy,idummy,idummy,idummy
!
    WRITE(6,'(a)') ' Allocating volume memory '
    ALLOCATE(itimvol(nelev))
    ALLOCATE(vnyqvol(nelev))
    ALLOCATE(rngvol(ngate,nelev))
    ALLOCATE(azmvol(nazimin,nelev))
    ALLOCATE(elvvol(nazimin,nelev))
    ALLOCATE(refvol(ngate,nazimin,nelev))
    IF(wrtuaref /= 0) ALLOCATE(uarefvol(ngate,nazimin,nelev))
    ALLOCATE(velvol(ngate,nazimin,nelev))
    IF(wrtuavel /= 0) ALLOCATE(uavelvol(ngate,nazimin,nelev))
    IF(wrtvort  /= 0) ALLOCATE(vorvol(ngate,nazimin,nelev))

    WRITE(6,'(a)') ' Reading Time Array'
    READ(iunit) itimvol
    WRITE(6,'(a)') ' Reading Reflectiving Arrays'
    READ(iunit) rngvol
    READ(iunit) azmvol
    READ(iunit) elvvol
    READ(iunit) refvol
    IF( wrtuaref /= 0 ) READ(iunit) uarefvol

    WRITE(6,'(a)') ' Reading Velocity Arrays'
    READ(iunit) vnyqvol
    READ(iunit) rngvol
    READ(iunit) azmvol
    READ(iunit) elvvol
    READ(iunit) velvol
    IF( wrtuavel /= 0 ) READ(iunit) uavelvol
    IF( wrtvort  /= 0 ) READ(iunit) vorvol
!
    WRITE(6,'(a)') ' Reading successfully completed'
    CLOSE(iunit)
!
    CALL ctim2abss(year,month,day,hour,minute,second,itime)
    itimcdf=itime-itim1970
    ifsecs=NINT(curtim)
    initime=itimcdf-ifsecs
!
! Count radials in input arrays within specified sector
! or set equal to value in the entire data array.
!
    IF(sectorize) THEN
      nazim=0
      DO kelv=1,nelev
        kazim=0
        DO jazim=1,nazimin
          IF(azmvol(jazim,kelv) >= 0.) THEN
            dazim=azmvol(jazim,kelv)-azimbgn
            IF(dazim < 0.) dazim=dazim+360.
            IF(dazim <= dazmax) kazim=kazim+1
          END IF
        END DO
        nazim=max(nazim,kazim)
      END DO
      WRITE(6,'(i6,a,f9.2,a,f9.2)') &
         nazim,' radials found in range ',azimbgn,' to ',azimend
    ELSE
      nazim=nazimin
    END IF

    IF(nazim == 0) CYCLE

    ALLOCATE(azim(nazim))
    ALLOCATE(beamw(nazim))
    ALLOCATE(gtspc(nazim))
    ALLOCATE(vnyq(nazim))
    ALLOCATE(refl(ngate,nazim))
    IF(wrtuaref /= 0) ALLOCATE(uarefl(ngate,nazim))
    ALLOCATE(radv(ngate,nazim))
    IF(wrtuavel /= 0) ALLOCATE(uaradv(ngate,nazim))
    IF(wrtvort  /= 0) ALLOCATE(vort(ngate,nazim))
!
!-----------------------------------------------------------------------
!
! Open index record file
!
!-----------------------------------------------------------------------
!
    IF( ifile == 1 .AND. creidx > 0 ) THEN
     WRITE(idxfname,'(a,a,a,a,a,a,i6.6,a)') &
     TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
     CALL getunit(idxunit)
     WRITE(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
     OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
     WRITE(idxunit,'(a)') '<?xml version="1.0" encoding="iso-8859-1" ?>'
     WRITE(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',   &
         TRIM(outdir),'">'
    END IF

    DO itilt=1,nelev

      azim=0.
      beamw=0.
      gtspc=0.
      refl=0.
      radv=0.
      IF(wrtuaref /= 0) uarefl=0.
      IF(wrtuavel /= 0) uaradv=0.
      IF(wrtvort  /= 0) vort=0.
      CALL abss2ctim(itimvol(itilt),iyear,imon,iday,ihour,imin,isec)

      IF(sectorize) THEN
!
! Find rotating direction and starting azimuth
!
        irot=1
        dazim=azmvol(2,itilt)-azmvol(1,itilt)
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
            IF(azmvol(jazim,itilt) >= 0.) THEN
              dazim=azmvol(jazim,itilt)-azimbgn
              IF(dazim < 0.) dazim=dazim+360.
              IF(dazim < dazmin) THEN
                dazmin=dazim
                jstart=jazim
              END IF
            END IF
          END DO
        ELSE
          DO jazim=1,nazimin
            IF(azmvol(jazim,itilt) >= 0.) THEN
              dazim=azimend-azmvol(jazim,itilt)
              IF(dazim < 0.) dazim=dazim+360.
              IF(dazim < dazmin) THEN
                dazmin=dazim
                jstart=jazim
              END IF
            END IF
          END DO
        END IF
  !
        kazim=0
        DO jazim=jstart,nazimin
          IF(azmvol(jazim,itilt) >= 0.) THEN
            dazim=azmvol(jazim,itilt)-azimbgn
            IF(dazim < 0.) dazim=dazim+360.
            IF(dazim <= dazmax) THEN
              kazim=kazim+1
              azim(kazim)=azmvol(jazim,itilt)
              beamw(kazim)=beamwid
              gtspc(kazim)=rngvol(2,itilt)-rngvol(1,itilt)
              vnyq(kazim)=vnyqvol(itilt)
              DO igate=1,ngate
                refl(igate,kazim)=refvol(igate,jazim,itilt)
                radv(igate,kazim)=velvol(igate,jazim,itilt)
              END DO
              IF(wrtuaref /= 0) THEN
                DO igate=1,ngate
                  uarefl(igate,kazim)=uarefvol(igate,jazim,itilt)
                END DO
              END IF
              IF(wrtuavel /= 0) THEN
                DO igate=1,ngate
                  uaradv(igate,kazim)=uavelvol(igate,jazim,itilt)
                END DO
              END IF
              IF(wrtvort  /= 0) THEN
                DO igate=1,ngate
                  vort(igate,kazim)=vorvol(igate,jazim,itilt)
                END DO
              END IF
            END IF
          END IF
        END DO

        DO jazim=1,(jstart-1)
          IF(azmvol(jazim,itilt) >= 0.) THEN
            dazim=azmvol(jazim,itilt)-azimbgn
            IF(dazim < 0.) dazim=dazim+360.
            IF(dazim <= dazmax) THEN
              kazim=kazim+1
              azim(kazim)=azmvol(jazim,itilt)
              beamw(kazim)=beamwid
              gtspc(kazim)=rngvol(2,itilt)-rngvol(1,itilt)
              vnyq(kazim)=vnyqvol(itilt)
              DO igate=1,ngate
                refl(igate,kazim)=refvol(igate,jazim,itilt)
                radv(igate,kazim)=velvol(igate,jazim,itilt)
              END DO
              IF(wrtuaref /= 0) THEN
                DO igate=1,ngate
                  uarefl(igate,kazim)=uarefvol(igate,jazim,itilt)
                END DO
              END IF
              IF(wrtuavel /= 0) THEN
                DO igate=1,ngate
                  uaradv(igate,kazim)=uavelvol(igate,jazim,itilt)
                END DO
              END IF
              IF(wrtvort  /= 0) THEN
                DO igate=1,ngate
                  vort(igate,kazim)=vorvol(igate,jazim,itilt)
                END DO
              END IF
            END IF
          END IF
        END DO
!
      ELSE   ! no sectorizing
!
        kazim=0
        DO jazim=1,nazimin
          azim(jazim)=azmvol(jazim,itilt)
          beamw(jazim)=beamwid
          gtspc(jazim)=rngvol(2,itilt)-rngvol(1,itilt)
          vnyq(jazim)=vnyqvol(itilt)
          DO igate=1,ngate
            refl(igate,jazim)=refvol(igate,jazim,itilt)
            radv(igate,jazim)=velvol(igate,jazim,itilt)
          END DO
          IF(wrtuaref /= 0) THEN
            DO igate=1,ngate
              uarefl(igate,jazim)=uarefvol(igate,jazim,itilt)
            END DO
          END IF
          IF(wrtuavel /= 0) THEN
            DO igate=1,ngate
              uaradv(igate,jazim)=uavelvol(igate,jazim,itilt)
            END DO
          END IF
          IF(wrtvort  /= 0) THEN
            DO igate=1,ngate
              vort(igate,jazim)=vorvol(igate,jazim,itilt)
            END DO
          END IF
        END DO

      END IF
!
      rfrgate=rngvol(1,itilt)

      elv=elvvol(1,itilt)
      ielv=INT(elv)
      felv=NINT(100.*(elv-ielv))
      vnyqstl=vnyqvol(itilt)
  
      WRITE(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')        &
         'Reflectivity_',ielv,'.',felv,'_',iyear,imon,iday,'-',        &
                     ihour,imin,isec,'.netcdf'
      WRITE(6,'(a,a)') ' Reflectivity data file: ',TRIM(fname)
      varname='Reflectivity'
      CALL wtrftiltcdf(ngate,nazim,nazim,fname,outdir,varname,         &
                       radname,radlat,radlon,radelv,vcp,elv,           &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyqstl,rfrgate,                                &
                       azim,beamw,gtspc,refl)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
      IF( creidx > 0 ) THEN
        WRITE(idxunit,'(a)') '<item>'
        WRITE(idxunit,'(a,f8.6,a,i10,a)') &
          '<time fractional="',frtime,'">',itimcdf,'</time>'
        WRITE(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',  &
           TRIM(fname),TRIM(cmprext),'</params>'
        WRITE(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)') &
          '<selections>',iyear,imon,iday,'-',                          &
          ihour,imin,isec,'Reflectivity',ielv,'.',felv,                &
          '</selections>'
        WRITE(idxunit,'(a)') '</item>'
      END IF
!
      IF(wrtuaref > 0) THEN
        WRITE(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
          'ReflectivityUnatten_',ielv,'.',felv,'_',iyear,imon,iday,'-',&
                      ihour,imin,isec,'.netcdf'
        WRITE(6,'(a,a)') ' Unattenuated Reflectivity: ',TRIM(fname)
        varname='UnattenuatedReflectivity'
        CALL wtrftiltcdf(ngate,nazim,nazim,fname,outdir,varname,       &
                         radname,radlat,radlon,radelv,vcp,elv,         &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                              &
                         azim,beamw,gtspc,uarefl)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         WRITE(idxunit,'(a)') '<item>'
         WRITE(idxunit,'(a,f8.6,a,i10,a)')                             &
            '<time fractional="',frtime,'">',itimcdf,'</time>'
         WRITE(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ', &
             TRIM(fname),TRIM(cmprext),'</params>'
         WRITE(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'UnattenuatedReflectivity',ielv,           &
            '.',felv,'</selections>'
         WRITE(idxunit,'(a)') '</item>'
        END IF
      END IF
!
!-----------------------------------------------------------------------
!
!  Write tilt netCDF file - Velocity
!
!-----------------------------------------------------------------------
!
      WRITE(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')        &
         'Velocity_',ielv,'.',felv,'_',iyear,imon,iday,'-',            &
                  ihour,imin,isec,'.netcdf'
      WRITE(6,'(a,a)') ' Velocity data file: ',TRIM(fname)
      varname='Velocity'
      CALL wtvrtiltcdf(ngate,nazim,nazim,fname,outdir,varname,         &
                       radname,radlat,radlon,radelv,vcp,elv,           &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyqstl,rfrgate,                                &
                       azim,beamw,gtspc,vnyq,radv)
!
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
      IF( creidx > 0 ) THEN
       WRITE(idxunit,'(a)') '<item>'
       WRITE(idxunit,'(a,f8.6,a,i10,a)')                               &
         '<time fractional="',frtime,'">',itimcdf,'</time>'
       WRITE(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ',   &
           TRIM(fname),TRIM(cmprext),'</params>'
       WRITE(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')  &
          '<selections>',iyear,imon,iday,'-',                          &
          ihour,imin,isec,'Velocity',ielv,'.',felv,                    &
          '</selections>'
       WRITE(idxunit,'(a)') '</item>'
      END IF

      IF(wrtuavel > 0) THEN
        WRITE(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'VelocityUnatten_',ielv,'.',felv,'_',iyear,imon,iday,'-',     &
                   ihour,imin,isec,'.netcdf'
        WRITE(6,'(a,a)') ' Unattended Velocity: ',TRIM(fname)
        varname='UnattenuatedVelocity'
        CALL wtvrtiltcdf(ngate,nazim,nazim,fname,outdir,varname,       &
                         radname,radlat,radlon,radelv,vcp,elv,         &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                              &
                         azim,beamw,gtspc,vnyq,radv)
        IF( creidx > 0 ) THEN
         WRITE(idxunit,'(a)') '<item>'
         WRITE(idxunit,'(a,f8.6,a,i10,a)')                             &
           '<time fractional="',frtime,'">',itimcdf,'</time>'
         WRITE(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ', &
             TRIM(fname),TRIM(cmprext),'</params>'
         WRITE(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Velocity',ielv,'.',felv,                  &
            '</selections>'
         WRITE(idxunit,'(a)') '</item>'
        END IF
      END IF

      IF(wrtvort > 0) THEN
        WRITE(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
           'Vorticity_',ielv,'.',felv,'_',iyear,imon,iday,'-',         &
                     ihour,imin,isec,'.netcdf'
        WRITE(6,'(a,a)') ' Vertical Vorticity: ',TRIM(fname)
        varname='VerticalVorticity'
        CALL wtvvtiltcdf(ngate,nazim,nazim,fname,outdir,varname,       &
                         radname,radlat,radlon,radelv,vcp,elv,         &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                             &
                         azim,beamw,gtspc,vort)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         WRITE(idxunit,'(a)') '<item>'
         WRITE(idxunit,'(a,f8.6,a,i10,a)')                             &
            '<time fractional="',frtime,'">',itimcdf,'</time>'
         WRITE(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ', &
             TRIM(fname),TRIM(cmprext),'</params>'
         WRITE(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Vorticity',ielv,'.',felv,                 &
            '</selections>'
         WRITE(idxunit,'(a)') '</item>'
        END IF
      END IF
!
!-----------------------------------------------------------------------
!
!
    END DO  ! itilt
!
!   Free memory as next file may be larger or smaller
!
    DEALLOCATE(itimvol)
    DEALLOCATE(vnyqvol)
    DEALLOCATE(rngvol)
    DEALLOCATE(azmvol)
    DEALLOCATE(elvvol)
    DEALLOCATE(refvol)
    IF( wrtuaref > 0) DEALLOCATE(uarefvol)
    DEALLOCATE(velvol)
    IF( wrtuavel > 0) DEALLOCATE(uavelvol)
    IF( wrtvort  > 0) DEALLOCATE(vorvol)

    DEALLOCATE(azim)
    DEALLOCATE(beamw)
    DEALLOCATE(gtspc)
    DEALLOCATE(vnyq)
    DEALLOCATE(refl)
    IF( wrtuaref > 0) DEALLOCATE(uarefl)
    DEALLOCATE(radv)
    IF( wrtuavel > 0) DEALLOCATE(uaradv)
    IF( wrtvort  > 0) DEALLOCATE(vort)

    IF(cmprsd) CALL cmprsgz(rfnamevol)
!
!-----------------------------------------------------------------------
  END DO  ! ifile
!
! Finish index file and close it out
!
  IF( creidx > 0 ) THEN
    WRITE(idxunit,'(a)') '</codeindex>'
    CLOSE(idxunit)
    CALL retunit(idxunit)
  END IF

END PROGRAM RADBIN2CDF

SUBROUTINE extrctdir(filestr,dirstr)
!
! Extract directory name from filename
! Keith Brewster, CAPS
! March 22, 2005
!
  IMPLICIT NONE
  CHARACTER (LEN=*) :: filestr
  CHARACTER (LEN=*) :: dirstr
  INTEGER :: lenfil
  INTEGER :: lendir
  INTEGER :: i,jloc
  lenfil=LEN(filestr)
  lendir=LEN(dirstr)
  jloc=0
  dirstr='.'
  DO i=1,lenfil
    IF(filestr(i:i) == '/') jloc=i
  END DO
  IF(jloc > 1) THEN
    IF((jloc-1) <= lendir) THEN
      dirstr=filestr(1:(jloc-1)) 
    ELSE
      WRITE(6,'(a,i6)') 'Insufficient length for dirstr, need',(jloc-1)
    END IF
  END IF
  RETURN
END SUBROUTINE extrctdir
