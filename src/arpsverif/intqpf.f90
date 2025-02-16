!########################################################################
!########################################################################
!######                                                            ######
!######                     PROGRAM INTQPF                         ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM intqpf

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Reads in precipitation data from ARPS history dumps and GRIB
! files, calculates total precipitation, and interpolates the
! data to a new grid.
!
! Based on program INTGRB (NCEP) and program ARPSENSCV (CAPS).
!
! AUTHOR:  Eric Kemp, February 2000.
!
! MODIFICATION HISTORY:
! Eric Kemp, March 2000.
! Conversion to FORTRAN 90.
!
! Eric Kemp, 31 March 2000.
! Added lat/lon coordinates of four corners for NCL plotting.
!
!-----------------------------------------------------------------------
!
! Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc' ! FORTRAN 77 include file.

!-----------------------------------------------------------------------
!
! Other variables
!
!-----------------------------------------------------------------------

  INTEGER :: intnx,intny
  REAL,ALLOCATABLE :: tpcp2d(:,:)
  INTEGER,PARAMETER :: maxtimes = 12
  REAL :: qswcorn,qswcorw,qnecorn,qnecorw

  INTEGER,PARAMETER :: IMAXIN=1200,JJMAXIN=900,JMAXIN=IMAXIN*JJMAXIN
  INTEGER,PARAMETER :: IMAXOT=200,JJMAXOT=150,JMAXOT=IMAXOT*JJMAXOT
  INTEGER,PARAMETER :: IBUFSIZE = 50000000

  INTEGER :: KPDS(25),KGDS(25)
  INTEGER :: KPDSO(25),KGDSO(25)
  INTEGER :: KPTR(20)
  INTEGER :: JSTAT(25),IPOPT(20)

  LOGICAL*1 :: KBMS(JMAXIN),KBMSO(JMAXOT)

  REAL :: DATA(JMAXIN),DATAO(JMAXOT),RLAT(JMAXOT),RLON(JMAXOT)
  REAL :: tpcp(JMAXIN),tpcpint(JMAXOT)
  REAL :: tpcpold(JMAXIN)

  REAL, ALLOCATABLE, DIMENSION(:,:) :: lat2d,lon2d
  REAL :: scswlat,scswlon
  REAL :: mapswlat,mapswlon,mapnwlat,mapnwlon
  REAL :: mapselat,mapselon,mapnelat,mapnelon

  CHARACTER*1 :: IBUF(IBUFSIZE),MSGA(200+17*JMAXIN/8),GDS(400)
  CHARACTER :: NFILE*132

  INTEGER :: ILEVL,IPDS16,ILVLT,IPARM
  INTEGER :: IYEAR,IMON,IDAY,IHR,IFHR1,IFHR2
  INTEGER :: date(6)
  INTEGER :: IFHRDIFF
  INTEGER :: kk,i,kj,iret,IGRID,IERR,KBYTES,MERR,j,ll
  INTEGER :: ii,jj
  INTEGER :: LUGBIN,LUGBOT,IBI,IJK,NDATA

  REAL :: RMAX,RMIN

  INTEGER :: LENGDS,JERR2,JSTART,I1,IRCBYTE,KRET
  INTEGER :: ko,ibo
  INTEGER :: foundpcp,addpcp
  INTEGER :: iabssec,isec
  INTEGER :: houroday, minohour, dayoweek, dateomon, month, year
  INTEGER :: countoutfile,ifile
  INTEGER :: lenrun,grbflen,len1

  CHARACTER*13 :: gribtime

  INTEGER :: kftime,initsec,jabssec
  INTEGER :: imin,imo,iyr,jldy

  CHARACTER*9 :: julfname

  INTEGER :: ifhr,ifmin,ifsec,mfhr
  INTEGER :: myr

  CHARACTER*9 :: extdfcst
  CHARACTER*19 :: extdinit

  INTEGER :: ilength
  INTEGER :: ijq
  REAL,ALLOCATABLE,DIMENSION(:,:) :: x2n,y2n
  REAL,PARAMETER :: missing = -9999.

  INTEGER :: tt
  INTEGER :: timesec(maxtimes),tmptimesec

  LOGICAL :: firstfile

  INTEGER :: status
  INTEGER :: grbflen
  INTEGER :: oldifhr1

!-----------------------------------------------------------------------
!
! HDF arrays and variables for output
!
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: flag_p = 0,flag_sfc = 1

  INTEGER,PARAMETER :: npreslevel=1,nvar_p=1,nvar_sfc=1
  REAL,ALLOCATABLE :: tprecip4d(:,:,:,:)
  REAL,ALLOCATABLE :: var_p(:,:,:,:,:)

  CHARACTER*4 :: varid_p(nvar_p)
  CHARACTER*16 :: varname_p(nvar_p)
  CHARACTER*3 :: varunit_p(nvar_p)
  REAL,PARAMETER :: pressure(nvar_p) = -9999.

  CHARACTER*4,PARAMETER :: varid_tpcp(nvar_sfc) = 'APCP'
  CHARACTER*32,PARAMETER :: varname_tpcp(nvar_sfc) = 'Acc. Precip.'
  CHARACTER*32,PARAMETER :: varunit_tpcp(nvar_sfc) = 'mm'

  REAL :: intscale,inttrulat1,inttrulat2,inttrulon,intdx,intdy,  &
          intswlat,intswlon
  INTEGER :: intmapproj

!-----------------------------------------------------------------------
!
! Namelists
!
!-----------------------------------------------------------------------

  INTEGER :: infiletype,arpsfmt,nx,ny,nz,numinfiles
  INTEGER :: iimin,iimax,jjmin,jjmax
  INTEGER :: sumpcp,timepcpopt,numtimes
  INTEGER,PARAMETER :: maxfiles = 36
  CHARACTER*29 :: extdtime(maxfiles)
  CHARACTER*60 :: dir_extd
  CHARACTER*132 :: extdname
  CHARACTER*132 :: infilename(maxfiles)
  INTEGER :: starthr(maxfiles)
  CHARACTER*132 :: outfile
  INTEGER :: igrid,interp,nintopt
  INTEGER,PARAMETER :: maxintopt = 20
  INTEGER :: intopt(maxintopt)
  CHARACTER*8 :: model
  CHARACTER*4 :: cgrid

  NAMELIST /fileinfo/ infiletype,arpsfmt,nx,ny,nz, &
                      iimin,iimax,jjmin,jjmax, &
                      sumpcp,timepcpopt,numtimes, &
                      numinfiles,extdtime, starthr,&
                      dir_extd,extdname,infilename,outfile,igrid, &
                      interp,nintopt,intopt,model,cgrid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  oldifhr1 = 0
  kbms = .true.
  tmptimesec = 0
  DO j = 1,maxtimes
    timesec(j) = 0
  END DO

  firstfile = .true.
  foundpcp = 0
  addpcp = 0
  DO j = 1,jmaxin
    tpcp(j) = REAL(0)
    tpcpold(j) = REAL(0)
  ENDDO

  LUGBIN=11
  LUGBOT=51
  ipopt(1)=-1
  ipopt(2)=-1
  DO i=3,20
    ipopt(i)=0
  END DO

  WRITE(6,*)'Welcome to program INTQPF.  This program reads in'
  WRITE(6,*)'gridded precipitation data and interpolates it to'
  WRITE(6,*)'another grid.'
  WRITE(6,*)

!-----------------------------------------------------------------------
!
! Read input file
!
!-----------------------------------------------------------------------

  READ(5,fileinfo,END=100)
  100  CONTINUE

!-----------------------------------------------------------------------
!
! Perform some error checking...
!
!-----------------------------------------------------------------------

  IF (numinfiles > maxfiles) THEN
    WRITE(6,*)'ERROR -- numinfiles > maxfiles'
    WRITE(6,*)'Reset maxfiles to ',numinfiles,' and recompile.'
    STOP
  ENDIF

  IF (nintopt > maxintopt) THEN
    WRITE(6,*)'ERROR -- nintopt > maxintopt'
    WRITE(6,*)'Change nintopt to 0 <= nintopt <= 20'
    WRITE(6,*)'nintopt = ',nintopt,' maxintopt = ',maxintopt
    STOP
  ELSE
    DO ll = 1,nintopt
      ipopt(ll) = intopt(ll)
    ENDDO
  ENDIF

  IF (sumpcp < 1 .OR. sumpcp > 2) THEN
    WRITE(6,*)'INTQPF: ERROR -- Invalid value for sumpcp.'
    WRITE(6,*)'sumpcp = ',sumpcp
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

  IF (timepcpopt < 1 .OR. timepcpopt > 3) THEN
    WRITE(6,*)'INTQPF: ERROR -- Invalid value for timepcpopt.'
    WRITE(6,*)'timepcpopt = ',timepcpopt
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

  IF (numtimes > numinfiles .AND. infiletype == 1 .AND. &
      timepcpopt == 3) THEN
    WRITE(6,*)'INTQPF: ERROR -- numtimes > numinfiles'
    WRITE(6,*)'numtimes = ',numtimes,' numinfiles = ', &
               numinfiles
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

!-----------------------------------------------------------------------
!
! Create GDS information for new grid
!
!-----------------------------------------------------------------------

  CALL MAKGDS(IGRID,KGDSO,GDS,LENGDS,IERR)
  IF (KGDSO(1) == 1) THEN ! Mercator
    intmapproj = 3
    intscale = 1.
    inttrulat1 = REAL(KGDSO(9))*0.001
    inttrulat2 = REAL(KGDSO(9))*0.001
    inttrulon = missing ! True lon is not included in GRIB!!!
    intdx = REAL(KGDSO(12))
    intdy = REAL(KGDSO(13))
    intnx = KGDSO(2)
    intny = KGDSO(3)
  ELSE IF (KGDSO(1) == 3) THEN ! Lambert Conformal
    intmapproj = 2
    intscale = 1.
    inttrulat1 = REAL(KGDSO(12))*0.001
    inttrulat2 = REAL(KGDSO(13))*0.001
    inttrulon = REAL(KGDSO(7))*0.001
    intdx = REAL(KGDSO(8))
    intdy = REAL(KGDSO(9))
    intnx = KGDSO(2)
    intny = KGDSO(3)
  ELSE IF (KGDSO(1) == 5) THEN ! Polar Stereographic
    intmapproj = 1
    intscale = 1.
    IF (KGDSO(10).eq.0) THEN ! Northern Hemisphere
      inttrulat1 = 60. ! 60. degree latitude in hemisphere
      inttrulat2 = 60.
    ELSE ! Southern Hemisphere
      inttrulat1 = -60.
      inttrulat2 = -60.
    ENDIF
    inttrulon = REAL(KGDSO(7))*0.001
    intdx = REAL(KGDSO(8))
    intdy = REAL(KGDSO(9))
    intnx = KGDSO(2)
    intny = KGDSO(3)
  ELSE ! New projection doesn't have a dx!
    WRITE(6,*)'ERROR:  Invalid Map Projection.'
    intmapproj = missing
    intscale = missing
    inttrulat1 = missing
    inttrulat2 = missing
    inttrulon = missing
    intdx = missing
    intdy = missing
    STOP
  ENDIF

!-----------------------------------------------------------------------
!
! Allocate arrays.
!
!-----------------------------------------------------------------------


  ALLOCATE(tpcp2d(intnx,intny),lat2d(intnx,intny), &
           lon2d(intnx,intny),x2n(intnx,intny),y2n(intnx,intny), &
           tprecip4d(intnx,intny,maxtimes,nvar_sfc), &
           var_p(intnx,intny,npreslevel,maxtimes,nvar_p), &
           STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f1')

!-----------------------------------------------------------------------
!
! Loop through the files.
!
!-----------------------------------------------------------------------

  countoutfile = 1
  tt = 1
  DO ifile = 1,numinfiles

    IF (tt > maxtimes) THEN
      WRITE(6,*)'ERROR:  tt > maxtimes'
      WRITE(6,*)'Increase maxtimes and recompile'
      WRITE(6,*)'tt = ',tt,' maxtimes = ',maxtimes
      STOP
    ENDIF

!-----------------------------------------------------------------------
!
!   Determine the Julian date for the data.
!
!-----------------------------------------------------------------------

    READ(extdtime(ifile),'(a19,1x,a9)') extdinit,extdfcst
    IF (extdfcst.eq.'         ') extdfcst='000:00:00'
    READ(extdinit, &
        '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)',err=920,end=920) &
         iyear,imon,iday,ihr,imin,isec

    CALL julday(iyear,imon,iday,jldy)
    myr=mod(iyear,100)
    ifhr1=0
    ifmin=0
    ifsec=0
    READ(extdfcst, &
        '(i3,1x,i2,1x,i2)',err=4,end=4) ifhr1,ifmin,ifsec
    4 CONTINUE
    mfhr=mod(ifhr1,24)
    jldy = jldy + ifhr1/24
    WRITE(julfname,'(i2.2,i3.3,i2.2,i2.2)') myr,jldy,ihr,mfhr
    CALL ctim2abss(iyear,imon,iday,ihr,imin,isec,iabssec)
    jabssec=(ifhr1*3600) + (ifmin*60) + ifsec + iabssec

    IF (ifile == 1) THEN
      initsec = jabssec
    ENDIF
    kftime=jabssec-initsec
    curtim=float(kftime)

    WRITE(6,*)
    WRITE(6,'(a,a9,a,/19x,a,a19,a/a,i16,a,/,i26,a)') &
       ' Calling rdextfil, looking for ', &
         extdfcst,' hour forecast ', &
        'initialized at ',extdinit,' UTC', &
        ' Which is ',jabssec,' abs seconds or ',kftime, &
        ' seconds from the ARPS initial time.'

!-----------------------------------------------------------------------
!
!   Get precipitation data.
!
!-----------------------------------------------------------------------

    IF (infiletype == 1) THEN ! GRIB

      len1=LEN(dir_extd)
      grbflen=len1

      CALL strlnth( dir_extd, grbflen )

      IF( grbflen .EQ. 0 .OR. dir_extd(1:grbflen) .EQ. ' ' ) THEN
        dir_extd = '.'
        grbflen=1
      END IF

      IF( dir_extd(grbflen:grbflen) .NE. '/'.AND. &
          grbflen .LT. len1) THEN
        grbflen=grbflen+1
        dir_extd(grbflen:grbflen)='/'
      END IF

      CALL strlnth( dir_extd, grbflen )
      NFILE = dir_extd(1:grbflen)//'/'//infilename(ifile)

      ifhr2 = ifhr1
      ifhr1 = starthr(ifile)

      IF (ifhr1 /= oldifhr1 .AND. timepcpopt == 1) THEN
        tpcpold = 0  ! Reset previous accumulation to zero.
      END IF
      CALL rdgrbqpf(sumpcp,timepcpopt,numtimes, &
                    tpcp,tpcpold,foundpcp,addpcp,gribtime, &
                    dir_extd,nfile,kgds, &
                    ibi,kbms,iyear,imon,iday,ihr,ifhr1,ifhr2, &
                    jmaxin,ibufsize)
      oldifhr1 = ifhr1

      IF (IFHR2 > IFHR1) THEN
        tmptimesec = (IFHR2*3600) - (IFHR1*3600) + tmptimesec
      ELSE
        tmptimesec = (IFHR1*3600)
      END IF

      timesec(tt) = tmptimesec

    ELSEIF (infiletype == 2) THEN ! ARPS history dump
      IF ((nx*ny).gt.JMAXIN) THEN
        WRITE(6,*)'ERROR:  Dimensions of ARPS history dump '
        WRITE(6,*)'are larger than JMAXIN.  Reset JMAXIN to'
        WRITE(6,*)'be larger than nx*ny and recompile.'
        WRITE(6,*)'JMAXIN = ',JMAXIN,' nx*ny = ',nx*ny
        WRITE(6,*)'Aborting...'
        STOP
      ENDIF
      CALL rdarpsqpf(nx,ny,nz,tpcp,tpcpold,foundpcp,addpcp,nfile, &
                   kgds,ibi,kbms,iyear,imon,iday,ihr,ifhr1, &
                   ifhr2,jmaxin,ibufsize,dir_extd,extdname, &
                   extdfcst,arpsfmt,iimin,iimax,jjmin,jjmax)

      tmptimesec = (IFHR2*3600) - (IFHR1*3600)
      timesec(tt) = tmptimesec
    ELSE
      WRITE(6,*)'ERROR -- Invalid selection for infiletype.'
      WRITE(6,*)'infiletype = ',infiletype
      WRITE(6,*)'Aborting...'
      STOP
    ENDIF

    IF (firstfile) THEN
      date(1) = IYEAR
      date(2) = IMON
      date(3) = IDAY
      date(4) = IHR
      date(5) = 0
      date(6) = 0
      WRITE(6,*)'Date = ',date
      firstfile = .FALSE.
    ENDIF

!-----------------------------------------------------------------------
!
!   Interpolate field if precipitation data is ready.
!
!-----------------------------------------------------------------------

    IF (foundpcp == 1) THEN
      WRITE(6,*)'Interpolating precipitation field...'!EMK
      WRITE(6,*)

      CALL ipolates(interp,ipopt,kgds,kgdso,jmaxin,jmaxot,1,ibi, &
                    kbms,tpcp,ko,rlat,rlon,ibo,kbmso,datao,iret)

      kbms = .true. ! Reset bit map

      foundpcp = 0
      addpcp = 0
      countoutfile = countoutfile + 1

!-----------------------------------------------------------------------
!
!     Dump 1-D interpolated precip array into 2-D array.
!
!-----------------------------------------------------------------------

      IF (intnx /= kgdso(2).OR. intny /= kgdso(3)) THEN
        WRITE(6,*)'INTGRB: ERROR -- Output dimensions from '
        WRITE(6,*)'subroutine IPOLATES does not match those '
        WRITE(6,*)'set in the main program.  Reset intnx and '
        WRITE(6,*)'intny and then recompile.'
        WRITE(6,*)'intnx = ',intnx,' kgdso(2) = ',kgdso(2)
        WRITE(6,*)'intny = ',intny,' kgdso(3) = ',kgdso(3)
        STOP
      ELSE
        kk = 1
        DO jj = 1,intny
          DO ii = 1,intnx
            IF (kbmso(kk)) THEN
              tprecip4d(ii,jj,tt,1) = datao(kk)
              lat2d(ii,jj) = rlat(kk)
              lon2d(ii,jj) = rlon(kk)
            ELSE
              tprecip4d(ii,jj,tt,1) = missing
              lat2d(ii,jj) = rlat(kk)
              lon2d(ii,jj) = rlon(kk)
            ENDIF
            kk = kk + 1
          END DO
        END DO
        tt = tt + 1
        qswcorn = rlat(1)
        qswcorw = rlon(1) - 360.
        qnecorn = rlat(ko)
        qnecorw = rlon(ko) -360.
      ENDIF

      DO kk=1,25
        kpdso(kk)=kpds(kk)
      ENDDO
      kpdso(3)=igrid
      kpdso(4)=128+64*ibo

      IF (timepcpopt == 3 .AND. infiletype == 1) THEN
        tpcpold = 0.
      ENDIF
    ENDIF
  ENDDO

!-----------------------------------------------------------------------
!
! Dump interpolated precip. field to HDF file.
!
!-----------------------------------------------------------------------

  IF (intmapproj /= missing .AND. intscale /= missing .AND. &
      inttrulat1 /= missing .AND. inttrulat2 /= missing .AND. &
      inttrulon /= missing .AND. intdx /= missing .AND. &
      intdy /= missing) THEN

    tt = tt - 1

    scswlat = lat2d(1,1)
    scswlon = lon2d(1,1) - 360.

    CALL qpfmcorn(intmapproj,intscale,inttrulat1,inttrulat2, &
                    inttrulon,intdx,intdy,intnx,intny,&
                    lat2d,lon2d, &
                    1,intnx,1,intny, &
                    mapswlat,mapswlon,mapnwlat,mapnwlon, &
                    mapselat,mapselon,mapnelat,mapnelon)

    CALL wrtverif(intnx,intny,npreslevel,tt,nvar_p,nvar_sfc,missing, &
                  outfile,model,cgrid,date,timesec,pressure, &
                  intmapproj,intscale,inttrulat1,inttrulat2, &
                  inttrulon,intdx,intdy,scswlat,scswlon, &
                  mapswlat,mapswlon,mapnwlat,mapnwlon, &
                  mapselat,mapselon,mapnelat,mapnelon, &
                  flag_p,varid_p, varname_p, varunit_p, var_p, &
                  flag_sfc,varid_tpcp, varname_tpcp, varunit_tpcp, &
                  tprecip4d)
  ELSE
    WRITE(6,*)'ERROR:  Missing some grid information'
    STOP
  ENDIF

!-----------------------------------------------------------------------
!
! The end.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'Program INTQPF completed.'
  STOP

!-----------------------------------------------------------------------
!
! Problem doing time conversions.
!
!-----------------------------------------------------------------------

  920 CONTINUE
  WRITE(6,*) &
        ' Aborting, error in time format for external file', &
        ' External file time provided:',extdtime
  STOP
END PROGRAM intqpf
