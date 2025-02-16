      SUBROUTINE READMDF(maxsta,yyreq,moreq,ddreq,hhreq,mnreq,
     :              yydat,modat,dddat,hhdat,mndat,
     :              cname,nikname,rlat,rlon,elev,
     :              statime,solar,temp,dewpt,relh,ta9m,
     :              wdir,wspd,wgust,avgspd,wdsd,wssd,ws2m,
     :              stnpr,altim,precip,nsta,
     :              stareq,latest,oneonly,verbose,istat)
C
C     Subroutine to read Oklahoma Mesonet data in mdf format.
C     Based on similar code to read the d05 formatted mesonet data.
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
      parameter (extnam='.mdf')
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
      integer idyear,idmon,idday,idhour,idmin
      integer irmin,ihr,imin
      character*10 rtime
      character*1 dummy
C
C     Misc internal variables
C
      integer iline,ista,jsta,ksta,num,ielev,ntable,ldir
      integer maxdat,maxloc
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
C     The old table looked like this:
C
C 12341234512345671234567123451231234567890
C ADAX  001  34.80  96.67  300   ADA            PONTOTOC
C ALTU  002  34.59  99.34  415   ALTUS          JACKSON
C
C     The newfangled table (end of 1994) looks like
C
C                             1         2         3         4
C234 1234  123456789012345678901234567890123456789012345678901231231212->
C  1 ADAX  ADA            ADA            2 NNE    PONTOTOC       344756->
C12341212 12345
C  964009   296   2  8   1   0
C  2 ALTU  ALTUS          ALTUS          3 S      JACKSON        343514
C  992016   417   2  7   1   0
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
     +       num,tname(ksta),longnm(ksta),
     +       tlat(ksta),tlon(ksta),ielev
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
     +                      yyreq,moreq,ddreq,hhreq,mnreq,extnam
 805    FORMAT(a,a,i4.4,i2.2,i2.2,i2.2,i2.2,a4)
      END IF
      write(6,'(a,a)') 'READMDF: mesodata= ',mesodata
      open(3,err=990,file=mesodata,status='old',action='read')
      read(3,810,err=995,end=995) dummy
 810  FORMAT(a1)
      read(3,812,err=995,end=995) idyear,idmon,idday,idhour,idmin
 812  FORMAT(5x,i4,4i3)
      print *, 'READMDF: time is ',idyear,idmon,idday,idhour,idmin
      idyear=mod(idyear,100)
c
c     MDF Data format looks like this
c
c 101 ! (c) 2007 Oklahoma Climatological Survey - all rights reserved
c2345
c 22 2007 02 13 00 00 00
cSTID  STNM  TIME   RELH   TAIR   WSPD   WVEC  WDIR   WDSD   WSSD   WMAX    RAIN     PRES  SRAD   TA9M   WS2M   TS10   TB10   TS05   TB05   TS30    TR05    TR25    TR60    TR75
c12341234561234561231234123456712345671234567123123123456712345671234567123456781234567891234561234567
c   a     x     i  x   i      f      f      f  x  i      f      f      f       f        f     i      f
cADAX     1   115     99    9.6    1.8    1.8   146    9.6    0.4    2.5    0.25   971.26     0    9.5
c234567
c   0.9 -995.0 -995.0 -995.0 -995.0 -995.0 -995.00 -995.00 -995.00 -995.00
c
      jsta=0
C
      IF(oneonly) THEN
        DO 150 ista=1,maxsta
          read(3,815,err=995,end=151) inline
 815      FORMAT(a120)
          IF(inline(2:5).eq.stareq) THEN
            jsta=1
 820      FORMAT(1x,a4,6x,i6,3x,i4,f7.1,f7.1,f7.1,3x,i3,
     +         f7.1,f7.1,f7.1,f8.2,f9.2,i6,f7.1,f7.1)
            read(inline,820,err=995,end=151) 
     +        rname,irmin,irelh,
     +        rtair,rwspd,rwvec,iwdir,rwdsd,rwssd,rwmax,
     +        rrain,rpres,israd,rta9m,rws2m
            ihr=idhour+(irmin/60)
            imin=idmin+mod(irmin,60)
            WRITE(rtime,822) idyear,idmon,idday,ihr,imin
 822        FORMAT(i2.2,4i2.2)
            CALL MESTRANS(
     :      rname,rtime,irelh,rtair,rta9m,iwdir,rwvec,rwmax,
     :      rwspd,rwdsd,rwssd,rws2m,rrain,rpres,israd,
     :      maxsite,ntable,tname,longnm,tlat,tlon,telev,
     :      cname(jsta),nikname(jsta),rlat(jsta),rlon(jsta),elev(jsta),
     :      statime(jsta),solar(jsta),temp(jsta),dewpt(jsta),relh(jsta),
     :      ta9m(jsta),wdir(jsta),wspd(jsta),wgust(jsta),
     :      avgspd(jsta),wdsd(jsta),wssd(jsta),ws2m(jsta),
     :      stnpr(jsta),altim(jsta),precip(jsta),
     :      istat)
            IF(istat.lt.0) THEN
              jsta=0
              IF(verbose) write(6,830) rname
 830        FORMAT('  No match in station table found for station: ',a4,/,
     :             '  Data are discarded.')
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
        read(3,810,err=995,end=995) dummy
        DO 200 ista=1,maxsta
          read(3,820,err=995,end=201) 
     +      rname,irmin,irelh,
     +      rtair,rwspd,rwvec,iwdir,rwdsd,rwssd,rwmax,
     +      rrain,rpres,israd,rta9m,rws2m
          ihr=idhour+(irmin/60)
          imin=idmin+mod(irmin,60)
          WRITE(rtime,822) idyear,idmon,idday,ihr,imin
          jsta=jsta+1
          CALL MESTRANS(
     :      rname,rtime,irelh,rtair,rta9m,iwdir,rwvec,rwmax,
     :      rwspd,rwdsd,rwssd,rws2m,rrain,rpres,israd,
     :      maxsite,ntable,tname,longnm,tlat,tlon,telev,
     :      cname(jsta),nikname(jsta),rlat(jsta),rlon(jsta),elev(jsta),
     :      statime(jsta),solar(jsta),temp(jsta),dewpt(jsta),relh(jsta),
     :      ta9m(jsta),wdir(jsta),wspd(jsta),wgust(jsta),
     :      avgspd(jsta),wdsd(jsta),wssd(jsta),ws2m(jsta),
     :      stnpr(jsta),altim(jsta),precip(jsta),
     :      istat)
          IF(istat.lt.0) THEN
            jsta=jsta-1
            IF(verbose) write(6,830) rname
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
