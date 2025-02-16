      SUBROUTINE READMESO(maxsta,yyreq,moreq,ddreq,hhreq,mnreq,
     :              yydat,modat,dddat,hhdat,mndat,
     :              cname,nikname,rlat,rlon,elev,
     :              statime,solar,temp,dewpt,relh,ta9m,
     :              wdir,wspd,wgust,avgspd,wdsd,wssd,ws2m,
     :              stnpr,altim,precip,nsta,
     :              stareq,latest,oneonly,verbose,istat)
C
C     Subroutine to read mesonet data
C
C     Input instructions:
C       Set time request variables for time
C
C     All logic for determining "latest time available", etc
C     must be done by the calling program.  This routine will
C     simple report the status of acquiring data for the
C     requested time.  This is done through the yydat,modat,etc
C     time variables (date/time of latest data in file) and
C     the istat return parameter:
C
C     istat=0        successful open and read of data
C     istat=-1       error opening or reading station table file
C     istat=-2       error opening data file
C     istat=-3       error reading data file
C     istat=-4       no valid data in file
C     istat=-5       no match in station table for station
C
C     Keith Brewster
C     OU School of Meteorology
C     August, 1993
C     Update for changes in data filenames   Sept, 1993         ...KB
C     Added "oneonly" logic, to speed up reading for single     ...KB
C       station time series 
C     Added pass-back of RH and station pressure.  Nov, 1993    ...KB
C
C     Arguments
C
      implicit none
C
      integer maxsta
C
      integer yyreq,moreq,ddreq, ! Requested year,month,day    (opt. input)
     :        hhreq,mnreq        ! Requested hour and minute   (opt. input)
      integer yydat,modat,dddat,hhdat,mndat     ! Data time variables (output)
C
      character*5 cname(maxsta)         ! Station ID's
      character*20 nikname(maxsta)      ! Abbreviated City Names
      real rlat(maxsta),     ! Latitude North (degrees)
     :     rlon(maxsta),     ! Longitude West (degrees)
     :     elev(maxsta)      ! Elevation      (meters MSL)
      character*10 statime(maxsta)      ! Time as yrmodahhmn
      real solar(maxsta),     ! Solar radiation (W/m**2)
     :     temp(maxsta),      ! Temperature (C)
     :     dewpt(maxsta),     ! Dew Point (C)
     :     relh(maxsta),      ! Relative Humidity (%)
     :     ta9m(maxsta),      ! 9m air temp (C)
     :     wdir(maxsta),      ! Average vector wind direction (degrees)
     :     wspd(maxsta),      ! Average vector wind speed     (m/s)
     :     wgust(maxsta),     ! Wind gust      (m/s)
     :     avgspd(maxsta),    ! Average wind speed (m/s)
     :     wdsd(maxsta),      ! Wind direction standard deviation
     :     wssd(maxsta),      ! Wind speed standard deviation (m/s)
     :     ws2m(maxsta),      ! Wind speed at 2m (m/s)
     :     stnpr(maxsta),     ! Station pressure (mb)
     :     altim(maxsta),     ! Altimeter Setting (Inches of Hg)
     :     precip(maxsta)     ! Precipitation (mm since 00z)
      integer nsta
      character*4 stareq
      logical latest,oneonly,verbose
      integer istat
C
C     Location table information
C
      integer maxsite
      parameter (maxsite=200)
c      
      character*4 tname(maxsite)
      character*20 longnm(maxsite)
      real tlat(maxsite),tlon(maxsite),telev(maxsite)
      character*40 mesosites
      parameter (mesosites='/home/kbrews/tables/geomeso.tbl')
      character*80 envname
      character*80 dirname
      character*80 defdirn
      character*4 extnam
      parameter (extnam='.d05')
C
C     Raw, read-in observations
C
      integer idattm
      character*60 mesodata
      character*4 rname
      real rtair,rwspd,rwvec
      integer iwdir,irelh
      real rwdsd,rwssd,rwmax,rrain,rpres
      integer israd
      real rta9m,rws2m
      character*10 rtime
      character*1 dummy
C
C     Misc internal variables
C
      integer iline,ista,jsta,ksta,num,ielev,ntable,ldir
      integer maxdat,maxloc,ios
      character*6 instring
      character*120 inline
      logical tablerd
      data tablerd /.true./
C
C     Initialize environment name
C
      envname='MESODATA'
      defdirn='/data/mesonet/d05'
C
C     Initialize return parameters as zero.
C
      istat=0
      yydat=0
      modat=0
      dddat=0
      hhdat=0
      mndat=0
      nsta=0
C
C     Read in names and site information from table
C
C234 1234  123456789012345678901231234567890123456789012345678912345678->
C126 ARD2  Ardmore        Ardmore         3.5 ENE Carter        34.1926
C23456789012345
C -97.0857  266   2  8  1  0  20040223  20991231
C
      IF(tablerd) THEN
        open(3,file=mesosites,err=121,status='old',action='read')
        DO 90 iline=1,1000
          read(3,801,end=101) instring
 801      FORMAT(a6)
          IF(instring.eq.'[STOP]') GO TO 91
  90    CONTINUE
  91    CONTINUE
        DO 100 ksta=1,maxsite
          read(3,802,end=101) 
     :       num,tname(ksta),longnm(ksta),
     :       tlat(ksta),tlon(ksta),ielev
 802      FORMAT(i4,1x,a4,2x,a13,39x,f8.4,f10.4,i5)
          telev(ksta)=float(ielev)
 100    CONTINUE
 101    CONTINUE
        close(3)
        ntable=ksta-1
        tablerd=.false.
 121    CONTINUE
        istat=0
        IF(tablerd .OR. ntable.lt.80) THEN
          istat=-1
          write(6,803) ntable
 803      FORMAT('  Error reading station table in readmeso.',
     :           '  Info for ',i4,' stations read.  istat=-1')
          RETURN
        END IF
      END IF
C
      call getenv(envname,dirname)
      ldir=len(dirname)
      CALL strlnth(dirname,ldir)
      IF(ldir.le.1 .and. dirname(1:1).eq.' ') THEN
        dirname=defdirn
        ldir=len(dirname)
        CALL strlnth(dirname,ldir)
      END IF
      IF(latest) THEN
        write(mesodata,804) dirname(1:ldir),'/',extnam
 804    FORMAT(a,a,'latest',a4)
      ELSE
        write(mesodata,805) dirname(1:ldir),'/',
     :                      yyreq,moreq,ddreq,hhreq,mnreq,extnam
 805    FORMAT(a,a,i4.4,i2.2,i2.2,i2.2,i2.2,a4)
      END IF
      write(6,'(a,a)') ' mesodata= ',mesodata
      open(3,err=990,file=mesodata,status='old')
      read(3,810,err=995,end=995) dummy
      read(3,810,err=995,end=995) dummy
      read(3,810,err=995,end=995) dummy
      read(3,810,err=995,end=995) dummy
 810  FORMAT(a1)
c
c     Data format looks like this
c
cSTID  YYMMDDhhmm  RELH   TAIR   WSPD   WVEC  WDIR   WDSD   WSSD   WMAX    RAIN     PRES  SRAD   TA9M   WS2M
c12341212345678901231231234567123456712345671231231234567123456712345671234567812345678912345612345671234567
cADAX  9309081800    92   19.6    5.3    5.2   333    9.3    0.9    7.1   31.36   981.98    93  -98.0  -98.0
cALTU  9309081800    83   20.6    5.2    5.1    28   11.6    0.7    6.8    0.73   970.01   234   19.0    4.3
cALVA  9309081800    68   19.9    2.7    2.7   335   12.7    0.8    4.1    8.81   966.82   681   18.6    2.1
cANTL  9309081800    85   24.5    1.4    1.4   262   17.5    0.5    2.5   12.10   993.47   353  -98.0  -98.0
cARDM  9309081800    74   23.1    5.5    5.2   358   16.9    1.0    7.4    7.12   984.26   720  -98.0  -98.0
c
      jsta=0
C
      IF(oneonly) THEN
        DO 150 ista=1,maxsta
          read(3,815,err=995,end=151) inline
 815      FORMAT(a120)
          IF(inline(2:5) .EQ. stareq) THEN
            jsta=1
            istat=-99
            read(inline,820,iostat=ios) 
     :        rname,rtime,irelh,
     :        rtair,rwspd,rwvec,iwdir,rwdsd,rwssd,rwmax,
     :        rrain,rpres,israd,rta9m,rws2m
            IF( ios .EQ. 0 ) THEN
               CALL MESTRANS(rname,rtime,irelh,rtair,rta9m,
     :           iwdir,rwvec,rwmax,rwspd,rwdsd,rwssd,rws2m,
     :           rrain,rpres,israd,
     :           maxsite,ntable,tname,longnm,tlat,tlon,telev,
     :           cname(jsta),nikname(jsta),
     :           rlat(jsta),rlon(jsta),elev(jsta),
     :           statime(jsta),solar(jsta),temp(jsta),
     :           dewpt(jsta),relh(jsta),
     :           ta9m(jsta),wdir(jsta),wspd(jsta),wgust(jsta),
     :           avgspd(jsta),wdsd(jsta),wssd(jsta),ws2m(jsta),
     :           stnpr(jsta),altim(jsta),precip(jsta),istat)
            END IF
            IF(istat.lt.0) THEN
              jsta=0
              IF(verbose) write(6,830) rname
 830    FORMAT('  No match in station table found for station: ',a4,/,
     :         '  Data are discarded.')
            END IF
            GO TO 151
          END IF
 150    CONTINUE
        istat=-4
        IF(verbose) write(6,831) stareq
 831    FORMAT('  No data found for station: ',a4,/)
 151    CONTINUE
        istat=0
      ELSE
        DO 200 ista=1,maxsta
          read(3,815,err=995,end=201) inline
          read(inline,820,iostat=ios) 
     :        rname,rtime,irelh,
     :        rtair,rwspd,rwvec,iwdir,rwdsd,rwssd,rwmax,
     :        rrain,rpres,israd,rta9m,rws2m
 820      FORMAT(1x,a4,4x,a10,3x,i4,
     :         f7.1,f7.1,f7.1,3x,i3,
     :         f7.1,f7.1,f7.1,f8.2,f9.2,i6,f7.1,f7.1)
          istat=-99
          IF( ios .EQ. 0 ) THEN 
            jsta=jsta+1
            CALL MESTRANS(rname,rtime,irelh,rtair,rta9m,
     :        iwdir,rwvec,rwmax,rwspd,rwdsd,rwssd,rws2m,
     :        rrain,rpres,israd,
     :        maxsite,ntable,tname,longnm,tlat,tlon,telev,
     :        cname(jsta),nikname(jsta),
     :        rlat(jsta),rlon(jsta),elev(jsta),
     :        statime(jsta),solar(jsta),temp(jsta),
     :        dewpt(jsta),relh(jsta),
     :        ta9m(jsta),wdir(jsta),wspd(jsta),wgust(jsta),
     :        avgspd(jsta),wdsd(jsta),wssd(jsta),ws2m(jsta),
     :        stnpr(jsta),altim(jsta),precip(jsta),istat)
            IF(istat.lt.0) THEN
              jsta=jsta-1
              IF(verbose) write(6,830) rname
            END IF
          END IF
 200    CONTINUE
 201    CONTINUE
      END IF
      close(3)
      nsta=jsta
C
C     Keep track of max time reported among individual stations
C     This will be returned as obs time.
C
      IF(nsta.gt.0) THEN
        maxdat=-999
        maxloc=0
        DO 300 ista=1,nsta
          idattm=0
          read(statime(ista),840,err=295) 
     :       yydat,modat,dddat,hhdat,mndat
 840      FORMAT(i2,i2,i2,i2,i2)
          call icvtmin(yydat,modat,dddat,hhdat,mndat,idattm)
 295      IF(idattm.gt.maxdat) THEN
            maxdat=idattm
            maxloc=ista
          END IF
 300    CONTINUE
        IF(maxloc.gt.0) read(statime(maxloc),840,err=305) 
     :                     yydat,modat,dddat,hhdat,mndat
 305    CONTINUE
      ELSE
        IF(verbose) THEN
          istat=-4
          write(6,845) mesodata
 845      FORMAT('  No valid data found in mesonet data file.',/,
     :           '  Fname= ',a60,/,'  istat=-4')
        END IF
      END IF
      RETURN
 990  CONTINUE
      IF(verbose) THEN
        istat=-2
        write(6,850) mesodata
 850    FORMAT('  Error opening mesonet data file.',/,
     :           '  Fname= ',a60,/,'  istat=-2')
      END IF
      RETURN
 995  CONTINUE
      IF(verbose) THEN
        istat=-3
        write(6,860) mesodata
 860    FORMAT('  Error reading mesonet data file.',/,
     :           '  Fname= ',a60,/,'  istat=-3')
      END IF
      RETURN
      END
c
      FUNCTION sprtoalt(staprs,elev)
c
c     Converts station pressure (mb) to altimeter setting (inches).
c     From formula in Wallace and Hobbs.
c
c     Keith Brewster, OU School of Meteorology
c     August, 1993
c
      real sprtoalt
      real staprs,elev
      real sltemp,g,rd,tlapse,const,mbtoin
      real presmb
      parameter( sltemp=288.,g=9.80616,rd=287.04,tlapse=6.5,
     +           const=(-g/(tlapse*rd)*1000.),
     +           mbtoin=0.02953)
      presmb=staprs*((1.- (0.001*elev*tlapse/sltemp))**const)
      sprtoalt=mbtoin*presmb
      RETURN
      END
c
      FUNCTION VAPPRS(dewpc)
c
c     Vapor pressure in mb from Bolton's approximation.
c     Dew point in Celcius is input.
c     If temperature is input, this returns saturation
c     vapor pressure.
c
c     Keith Brewster, OU School of Meteorology
c     August, 1993
c
      implicit none
      real vapprs
      real dewpc
      vapprs=6.112 * exp( (17.67 * dewpc) / (dewpc + 243.5))
      RETURN
      END
c
      FUNCTION RHTODPC(rh,tempc)
c
c     Converts relative humidity to dew point in Celcius
c
c     Keith Brewster, OU School of Meteorology
c     August, 1993
c
      implicit none
      real rhtodpc
      real rh,tempc
      real vapor,vapprs,const,logvpr
c
      parameter (const=1.810254)     ! alog(6.112)
      vapor=(0.01*rh)*vapprs(tempc)
      logvpr=alog(vapor)
      rhtodpc= 243.5*(const - logvpr)/(logvpr - const - 17.67)
c
      RETURN
      END
c
      SUBROUTINE icvtmin(yydat,modat,dddat,hhdat,mndat,idattm)
c
c     Converts a time given as 5 integers describing year,
c     month,day,hour and minute into an integer describing
c     the number of minutes since Jan 1, 1980.
c
c     Useful for time differencing.
c
c     Keith Brewster, OU School of Meteorology
c     August, 1993
c
      IMPLICIT NONE
      integer yydat,modat,dddat,hhdat,mndat
      integer idattm
      integer dyear,julday
      integer datjul
c
      julday = datjul(yydat,modat,dddat)
c
      IF(yydat.gt.50) THEN
        dyear=yydat-80
        IF(dyear.gt.1000) dyear=dyear-1900
        IF(dyear.gt.100 ) dyear=dyear-100
      ELSE
        dyear=yydat+20
      END IF
c
      idattm=mndat+(60*hhdat)+(1440*julday)+(525600*dyear)
c
      RETURN
      END
c
      FUNCTION DATJUL(yy,mo,dd)
c
c     Given date as yy (year), mo (month), dd (day),
c     find the Julian day
c
c     Keith Brewster, OU School of Meteorology
c     August, 1993
c
      implicit none
      integer datjul
      integer yy,mo,dd
c
c     cudays is the number of cumulative days
c     preceeding the given month (excluding Feb 29)
c
      integer cudays(12)
      data cudays /0,31,59,90,120,151,181,212,243,273,304,334/
      integer knt
C
      knt=cudays(mo)+dd
      IF(mo.gt.2 .and. mod(yy,4).eq.0 .and. yy.ne.0) knt=knt+1
      datjul=knt
      RETURN
      END
C
      SUBROUTINE MESTRANS(
     :        rname,rtime,irelh,rtair,rta9m,
     :        iwdir,rwvec,rwmax,rwspd,rwdsd,rwssd,rws2m,
     :        rrain,rpres,israd,
     :        maxsite,ntable,tname,longnm,tlat,tlon,telev,
     :        cname,nikname,rlat,rlon,elev,
     :        statime,solar,temp,dewpt,relh,ta9m,
     :        wdir,wspd,wgust,avgspd,wdsd,wssd,ws2m,
     :        stnpr,altim,precip,istat)
      implicit none
C
C     Read-in data
C
      character*4 rname
      character*10 rtime
      integer irelh
      real rtair,rta9m
      integer iwdir
      real rwvec,rwmax,rwspd,rwdsd,rwssd,rws2m
      real rrain,rpres
      integer israd
C
C     Station table info
C
      integer maxsite,ntable
      character*4 tname(maxsite)
      character*20 longnm(maxsite)
      real tlat(maxsite),tlon(maxsite),telev(maxsite)
C
      character*5 cname
      character*20 nikname
      real rlat,rlon,elev
      character*10 statime
      real solar,temp,dewpt,ta9m,relh
      real wdir,wspd,wgust,avgspd,wdsd,wssd,ws2m
      real stnpr,altim,precip
      integer istat
C
C     Functions
C
      real sprtoalt,rhtodpc
C
C     Misc. internal variables
C
      integer ksta
C
      istat=0
C
C     Match station location info with the data read-in
C
C     Transfer data into output arrays doing type conversions
C     and a few data transformations along the way.
C
      DO 150 ksta=1,ntable
        IF(tname(ksta).EQ.rname) THEN
          write(cname,'(a)') rname
          nikname=longnm(ksta)
          rlat=tlat(ksta)
          rlon=tlon(ksta)
          elev=telev(ksta)
          statime=rtime
          temp=rtair
          ta9m=rta9m
          relh=float(irelh)
          IF(rtair.gt.-90. .and. irelh.gt.0) THEN
            dewpt=rhtodpc(float(irelh),rtair)
          ELSE
            dewpt=-98.
          END IF
          wdir=float(iwdir)
          wspd=rwvec
          avgspd=rwspd
          wgust=rwmax
          wdsd=rwdsd
          wssd=rwssd
          ws2m=rws2m
          solar=float(israd)
          stnpr=rpres
          IF(rpres.gt.0.) THEN
            altim=sprtoalt(rpres,elev)
          ELSE
            altim=rpres
          END IF
          precip=rrain
          GO TO 200
        END IF
 150  CONTINUE
C
C       At this point, no matching station has been found
C
      istat=-1      
 200  CONTINUE
      RETURN
      END
