      SUBROUTINE readwtx(maxsta,yyreq,moreq,ddreq,hhreq,mnreq,
     :           yydat,modat,dddat,hhdat,mndat,
     :           cname,rlat,rlon,elev,
     :           statime,temp,dewpt,
     :           wdir,wspd,gspd,altim,nsta,istat)
C
C  Subroutine to read surface data from the West Texas Mesonet
C
C  Input instructions:
C    Set time request variables for time
C
C  All logic for determining "latest time available", etc
C  must be done by the calling program.  This routine will
C  simple report the status of acquiring data for the
C  requested time.  This is done through the yydat,modat,etc
C  time variables (date/time of latest data in file) and
C  the istat return parameter:
C
C  istat=0        successful open and read of data
C  istat=-1       error opening or reading station table file
C  istat=-2       error opening data file
C  istat=-3       error reading data file
C  istat=-4       no valid data in file
C  istat=-5       no match in station table for station
C
C  Keith Brewster
C  Center for Analysis and Prediction of Storms
C  April, 2002
C
C  Arguments
C
      IMPLICIT NONE
C
      INTEGER maxsta
C
      INTEGER yyreq,moreq,ddreq,hhreq,mnreq  ! Requested year,month,day
      INTEGER yydat,modat,dddat,hhdat,mndat  ! Data time variables (output)
C
      CHARACTER*5 cname(maxsta)    ! Station ID's
      REAL rlat(maxsta)                  ! Latitude North (degrees)
      REAL rlon(maxsta)                  ! Longitude East (degrees)
      REAL elev(maxsta)                  ! Elevation      (meters MSL)
      CHARACTER*10 statime(maxsta) ! Time as yrmodahhmn
      REAL temp(maxsta)                  ! Temperature (C)
      REAL dewpt(maxsta)                 ! Dew Point (C)
      REAL wdir(maxsta)                  ! Average vector wind direction (degrees)
      REAL wspd(maxsta)                  ! Average vector wind speed     (m/s)
      REAL gspd(maxsta)                  ! Gust wind speed     (m/s)
      REAL altim(maxsta)                 ! Altimeter Setting (Inches of Hg)
      INTEGER nsta
      INTEGER istat
C
C Station table variables
C
      INTEGER maxsite
      PARAMETER(maxsite=100)
      CHARACTER*4 tname(maxsite)    ! Table Station ID's
      REAL tlat(maxsite)                 ! Table Latitude North (degrees)
      REAL tlon(maxsite)                 ! Table Longitude East (degrees)
      REAL telev(maxsite)                ! Table Elevation      (meters MSL)
      INTEGER dtime(maxsite)
C
C Misc internal variables
C
      CHARACTER*80 obstr
      CHARACTER*80 blankstr
      CHARACTER*60 wtxdata
      CHARACTER*4 staid
      INTEGER itime,itimreq,idattm,adt
      INTEGER iob,ista,jsta,ksta,maxdat,maxloc
      INTEGER iyr,imon,iday,ihr,imin,isec,myr
      INTEGER ntabsta,lendir,lenobs,istatus
      REAL dd,ffkt,gffkt,tc,tdc,altin
      CHARACTER*50 dirname
      CHARACTER*1 dummy
      LOGICAL lfound
C
C     dirname='/home/kbrews/mesoarch'
      dirname='/home/kbrews/mesonet'
      blankstr=
     :'                                                                '
C
C Initialize return parameters as zero.
C
      istat=0
      yydat=0
      modat=0
      dddat=0
      hhdat=0
      mndat=0
C
      iyr=yyreq
      myr=mod(iyr,100)
      imon=moreq
      iday=0
      ihr=0
      imin=0
      isec=0
C
      jsta=0 
      nsta=0
C
      call mkitime(yyreq,moreq,ddreq,hhreq,mnreq,isec,itimreq)
C
      OPEN(31,ERR=990,FILE='/home/kbrews/tables/wtxsta.table')
      READ(31,'(a1)') dummy
      DO ista=1,maxsite
        READ(31,'(a4,18x,f7.4,f12.4,f10.1)',END=101)
     :         tname(ista),tlat(ista),tlon(ista),telev(ista)
        dtime(ista)=-99
      END DO
 101  CONTINUE
      CLOSE(31)
      ntabsta=ista-1
C
      lendir=index(dirname,' ')-1
      WRITE(wtxdata,'(a,a1,i4.4,i2.2,i2.2,i2.2,a)') 
     :dirname(1:lendir),'/',yyreq,moreq,ddreq,hhreq,'wtx.txt'
      OPEN(31,ERR=990,FILE=wtxdata,STATUS='old')
      READ(31,'(a1)',END=501) dummy
      WRITE(6,'(a,a)')
     :'Staid  Day Hour Min   dd    ffkt   ',
     :'gffkt     tc     tdc     altin'
      DO iob=1,99999
        obstr=blankstr
        READ(31,'(a80)',iostat=istat) obstr
        IF(istat .ne. 0) GO TO 501
        lenobs=index(obstr,'=')
        IF(lenobs .gt. 12 ) THEN
          CALL DECODWTX(obstr,staid,iday,ihr,imin,dd,ffkt,gffkt,
     :                  tc,tdc,altin,istatus)
          WRITE(6,'(a,3i4,3f8.0,2f8.2,f10.2)')
     :          staid,iday,ihr,imin,dd,ffkt,gffkt,tc,tdc,altin
        END IF
        call mkitime(iyr,imon,iday,ihr,imin,isec,itime)
        adt=abs(itime-itimreq)
        IF(adt .lt. 1201 ) THEN
          lfound=.false.
          DO ista=1,ntabsta
            IF(tname(ista) .eq. staid) THEN
              IF(dtime(ista) < 0) THEN
                lfound=.true.
                jsta=jsta+1
                ksta=jsta
                cname(ksta)=staid
                rlat(ksta)=tlat(ista)
                rlon(ksta)=tlon(ista)
                elev(ksta)=telev(ista)
                dtime(ista)=adt
              ELSE IF (adt < dtime(ista) ) THEN
                DO ksta=1,jsta
                  IF(cname(ksta)(1:4) .eq. staid) THEN
                    lfound=.true.
                    dtime(ista)=adt
                    GO TO 107
                  END IF
                END DO
 107            CONTINUE
              END IF
              GO TO 111
            END IF 
          END DO
 111      CONTINUE
          IF(lfound) THEN
            write(statime(ksta),'(5(i2.2))') myr,imon,iday,ihr,imin
            temp(ksta)=tc
            dewpt(ksta)=tdc
            wdir(ksta)=dd
            wspd(ksta)=ffkt
            gspd(ksta)=gffkt
            altim(ksta)=altin
          END IF
        END IF
      END DO
 501  CONTINUE
      CLOSE(31)

      nsta=jsta
C
C  Keep track of max time reported among individual stations
C  This will be returned as obs time.
C
      IF(nsta .gt. 0) THEN
        maxdat=-999
        maxloc=0
        DO ista=1,nsta
          idattm=0
          READ(statime(ista),840,ERR=295)
     :         yydat,modat,dddat,hhdat,mndat
 840      FORMAT(i2,i2,i2,i2,i2)
          CALL icvtmin(yydat,modat,dddat,hhdat,mndat,idattm)
 295      IF(idattm .gt. maxdat) THEN
            maxdat=idattm
            maxloc=ista
          END IF
        END DO
        IF(maxloc .gt. 0) READ(statime(maxloc),840,ERR=305)
     :                      yydat,modat,dddat,hhdat,mndat
 305    CONTINUE
      ELSE
        istat=-4
        WRITE(6,845) wtxdata
 845    FORMAT('  No valid data found in West Texas data file.',/,
     :                    '  Fname= ',a60,/,'  istat=-4')
      END IF
      RETURN
 990  CONTINUE
      istat=-2
      WRITE(6,850) wtxdata
 850  FORMAT('  Error opening West Texas mesonet data file.',/,
     :                  '  Fname= ',a60,/,'  istat=-2')
      RETURN
 995  CONTINUE
      istat=-3
      WRITE(6,860) wtxdata
 860  FORMAT('  Error reading West Texas mesonet data file.',/,
     :                  '  Fname= ',a60,/,'  istat=-3')
      RETURN
      END
C
      SUBROUTINE DECODWTX(obstr,staid,iday,ihr,imin,
     :                  dd,ffkt,gffkt,tc,tdc,altin,istatus)
!
! An elementary METAR decoder for the obs from the West Texas
! Mesonet.  It is elementary because the West Texas Mesonet METARs
! do not include weather and cloud groups, and weather and cloud
! groups complicate METAR parsing.
!
! Keith Brewster
! CAPS
! November, 2001
!
      IMPLICIT NONE
      CHARACTER*(*) obstr
      CHARACTER*4 staid 
      INTEGER iday
      INTEGER ihr 
      INTEGER imin
      REAL dd
      REAL ffkt
      REAL gffkt
      REAL tc
      REAL tdc
      REAL altin
      INTEGER istatus
C
C       Misc. local variables
C
      INTEGER ivar,iloc,kloc,gloc,mloc,istat
      INTEGER lenobs
C
C       Initialize returned variables
C
      staid='XXXX'
      iday=99
      ihr=99
      imin=99
      dd=-999.
      ffkt=-999.
      gffkt=-999.
      tc=-999.
      tdc=-999.
      altin=-999.
      istatus=-1
      iloc=1
      gloc=0
      kloc=0
      mloc=0
C
C       Find length of obs string without trailing blanks
C
      lenobs=INDEX(obstr,'=')
C
C       Get station ID, date and time
C       
      IF(lenobs .gt. 11 .AND. obstr(12:12) .eq. 'Z') THEN
        read(obstr(1:4),'(a4)',iostat=istat) staid
        IF(istat .ne. 0) THEN
          istatus=-2
          RETURN
        END IF
        read(obstr(6:7),'(i2)',iostat=istat) iday
        IF(istat .ne. 0) THEN
          istatus=-2
          RETURN
        END IF
        read(obstr(8:9),'(i2)',iostat=istat) ihr 
        IF(istat .ne. 0) THEN
          istatus=-2
          RETURN
        END IF
        read(obstr(10:11),'(i2)',iostat=istat) imin
        IF(istat .ne. 0) THEN
          istatus=-2
          RETURN
        END IF
        istatus=0
      END IF
!
C       Get wind direction and speed
C       Note here we are using METARS that are always in KTs.
C       Speed can be two or three digits, and is followed by either
C       a "G", to precede gusts or "KT" as unit indicator.
!
      IF(lenobs .gt. 13) THEN
        read(obstr(14:16),'(i3)',iostat=istat) ivar
        dd=float(ivar)
      END IF
      IF(lenobs .gt. 16) THEN
        gloc=INDEX(obstr(17:lenobs),'G')
        kloc=INDEX(obstr(17:lenobs),'KT')
        IF(gloc .eq. 0 .OR. gloc .gt. kloc ) THEN
          IF(kloc .eq. 3 ) THEN
            read(obstr(17:18),'(i2)',iostat=istat) ivar
            ffkt=float(ivar)
          ELSE IF (kloc .eq. 4) THEN
            read(obstr(17:19),'(i3)',iostat=istat) ivar
            ffkt=float(ivar)
          END IF
        ELSE IF ( gloc .gt. 0 ) THEN
          IF(gloc .eq. 3) THEN
            read(obstr(17:18),'(i2)',iostat=istat) ivar
            ffkt=float(ivar)
          ELSE IF (gloc .eq. 4) THEN
            read(obstr(17:19),'(i3)',iostat=istat) ivar
            ffkt=float(ivar)
          END IF
        END IF
      END IF
      IF(lenobs .gt. 18) THEN
        iloc=18+kloc+2
        IF(gloc .gt. 0 .AND. gloc .lt. kloc) THEN
          IF((kloc-gloc) .eq. 3 ) THEN 
            read(obstr((17+gloc):(18+gloc)),'(i2)',iostat=istat) ivar
            gffkt=float(ivar)
          ELSE IF((kloc-gloc) .eq. 4) THEN
            read(obstr((17+gloc):(19+gloc)),'(i3)',iostat=istat) ivar
            gffkt=float(ivar)
          END IF
        END IF
      END IF
!
C       Get temperature
C       Note M is negative
C       Sample temperature/dewpoint group 19/12
C       Since there are no clouds in these METARS the first "/" is for
C       this group.
!
      IF(lenobs .gt. (iloc+1) ) THEN
        mloc=INDEX(obstr(iloc:lenobs),'M')
        kloc=INDEX(obstr(iloc:lenobs),'/')
        IF(mloc .ne. 0 .AND. mloc .lt. kloc ) THEN
          read(obstr((iloc+kloc-3):(iloc+kloc-2)),'(i2)',iostat=istat)
     :         ivar
          tc=-float(ivar)
        ELSE IF(kloc .gt. 0 ) THEN
          read(obstr((iloc+kloc-3):(iloc+kloc-2)),'(i2)',iostat=istat)
     :         ivar
          tc=float(ivar)
        END IF
        iloc=iloc+kloc
      END IF
!
C       Get dewpoint
C       Note M is negative
!      
      IF(lenobs .gt. (iloc+1)) THEN
        mloc=INDEX(obstr(iloc:lenobs),'M')
        kloc=INDEX(obstr(iloc:lenobs),' ')
        IF(mloc .ne. 0 .AND. mloc .lt. kloc ) THEN
          read(obstr((iloc+kloc-3):(iloc+kloc-2)),'(i2)',iostat=istat) 
     :         ivar
          tdc=-float(ivar)
        ELSE IF(kloc .gt. 0 ) THEN
          read(obstr((iloc+kloc-3):(iloc+kloc-2)),'(i2)',iostat=istat) 
     :         ivar
          tdc=float(ivar)
        END IF
        iloc=iloc+kloc
      END IF
!
C       Get altimeter setting
C       Sample: A3035
!
      IF(lenobs .gt. (iloc+1)) THEN
        kloc=INDEX(obstr(iloc:lenobs),'A')
        IF(lenobs .gt. (iloc+kloc+2) ) THEN
          read(obstr((iloc+kloc):(iloc+kloc+3)),'(i4)',iostat=istat) 
     :         ivar
          altin=0.01*float(ivar)
        END IF
        iloc=iloc+kloc+4
      END IF
!
C       Get temperature and dewpoints in tenths of degrees C from remarks group
C       Note when preceded by a "1" the temperature or dewpoint is negative.
C       Sample: RMK AO1 T01380125
!
      IF(lenobs .gt. iloc+12) THEN
        kloc=INDEX(obstr(iloc:lenobs),'RMK')
        iloc=iloc+kloc+2
        kloc=INDEX(obstr(iloc:lenobs),'T')
        iloc=iloc+kloc
        kloc=INDEX(obstr(iloc:lenobs),' ')
        IF(kloc .eq. 9 ) THEN
          read(obstr(iloc:(iloc+3)),'(i4)',iostat=istat) ivar
          IF(ivar .lt. 1000 ) THEN
            tc=0.1*float(ivar)
          ELSE
            tc=-0.1*float(mod(ivar,1000))
          END IF
          read(obstr((iloc+4):(iloc+7)),'(i4)',iostat=istat) ivar
          IF(ivar .lt. 1000 ) THEN
            tdc=0.1*float(ivar)
          ELSE
            tdc=-0.1*float(mod(ivar,1000))
          END IF
        END IF
      END IF
      RETURN
      END
