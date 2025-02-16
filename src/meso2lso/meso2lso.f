      PROGRAM wrtlso
C
C     Program to create a file of surface data for
C     use in the ADAS or ARPS 3DVAR analysis.
C
C     Keith Brewster, CAPS, January, 1994
C
C     MODIFICATIONS
C     Keith Brewster, July, 1998
C     Changed spatial limits to lat,lon limits specified in
C     a parameter statement for greater portability.  Consideration 
C     of a "grid-box" is thus removed.  For use in LAPS all obs 
C     are reported to be "in grid".  There is no affect on ADAS.
C
C     Keith Brewster, March, 2010
C     Replaced some time functions with ARPS lib versions of same.
C     Updated and cleaned-up.
C
      IMPLICIT NONE
C
C     Output data variables
C
      INTEGER maxsta,maxcld,nsrc,np
      PARAMETER (maxsta=5000,maxcld=3,nsrc=4,np=41)
      CHARACTER*5 stn(maxsta)
      REAL slat(maxsta),slon(maxsta),selv(maxsta)
      CHARACTER*8 obstype(maxsta)
      INTEGER obstime(maxsta)
      CHARACTER*8 wx(maxsta)
      REAL t(maxsta),td(maxsta),dd(maxsta),ff(maxsta),
     :     ddg(maxsta),ffg(maxsta),
     :     pstn(maxsta),pmsl(maxsta),alt(maxsta),rad(maxsta)
      REAL ceil(maxsta),lowcld(maxsta),cover(maxsta),vis(maxsta)
      INTEGER idp3(maxsta)
      INTEGER kloud(maxsta)
      character*1 emv(maxcld,maxsta)
      character*4 amt(maxcld,maxsta)
      real cldhgt(maxcld,maxsta)
      CHARACTER*8 tyname(nsrc)
      data tyname /'MESO    ','SA      ','MOBLS   ','MOBLM   '/
      INTEGER nin,nobsb,nobsg,nmeso,nmesord,nsao,nsaord,nsaog
      INTEGER nbuoyrd,nbuoy
C
C     Raw Input data (before screening, unit conversion, if any
C
      CHARACTER*5 stnin(maxsta)
      REAL slatin(maxsta),slonin(maxsta),selvin(maxsta)
      REAL tin(maxsta),tdin(maxsta),ddin(maxsta),ffin(maxsta),
     :     ffgin(maxsta),
     :     pstnin(maxsta),altin(maxsta),radin(maxsta)
      INTEGER mcode(maxsta)
C
C     Extra Mesonet-reading variables
C
      character*20 nikname(maxsta)
      character*10 statime(maxsta)      ! Time as yrmodahhmn
      real ta9m(maxsta),rh(maxsta)
      real avgspd(maxsta),wdsd(maxsta),wssd(maxsta)
      real ws2m(maxsta),altmeso(maxsta),precip(maxsta)
C
      CHARACTER*16 lsoname
C
C  Time stuff
C
      INTEGER imoreq,idareq,iyrreq,ihrreq,iminreq,isec
      INTEGER imorq,idarq,iyrrq,ihrrq,iminrq
      INTEGER imo,ida,iyr,ihr,imin
      CHARACTER*23 ctime
      integer isaotm,idelt,itrunc
      INTEGER itimex
      integer d5min,d15min
      parameter ( d5min  = 300,
     :            d15min = 900 )
c
      integer minmeso,npass,maxtry
      integer mesopass,ntry
c
c     "npass" is the number of mesonet files to use.
c     The second and subsequent files are used to fill-in those
c     missing in the first file.
c     "maxtry" is maximum number of mesonet files to try for
c     fills or in the event that less than 
c     "minmeso" stations are reported at a given time.  
c     "maxtry" must be greater than or equal to npass.
c
      parameter (minmeso=5)  ! minimum number of reporting mesonet stas
      parameter (npass=3)    ! number of mesonet files to use the
      parameter (maxtry=5)   ! total number of tries for mesonet data
      integer ndelta5(maxtry)
      data ndelta5 / -1,-2,0,-3,1 /
C
C  Some constants
C
      integer iflagd
      real flagd,chkmis
      parameter(iflagd=-99,flagd=-99.9,chkmis=-90.)
      real intomb,fttom,mstokts
      parameter(intomb=(1013.25/29.92),mstokts=1.944,
     :          fttom=0.3048)
      CHARACTER*29 area
C
C     Functions
C
      real tctotf
C
C     Misc internal variables
C
      CHARACTER*4 stareq
      CHARACTER*40 arg
      LOGICAL oneonly,verbose,ilatest,needtime,cmdtime
      REAL pres,tem,dew,altim
      REAL pi,d2rad,hmeters
      INTEGER istat,ierr,nwx,narg,ibgn,iend
      INTEGER ii,icld,icptr,itime,ista,jsta,ksta,msta,k
      INTEGER lunit,istatus

      INTEGER iargc
c
      idelt=0
      isec=0
      itime=0
      nsaog=0
      nsaord=0
      cmdtime=.false.
      needtime=.true.
      stareq='NORM'
      oneonly=.false.
      verbose=.true.
      ilatest=.false.
      pi=4.*atan(1.)
      d2rad=pi/180.
C
C     Get the input argument, the delta-time in hours
C
      narg=iargc()
      IF(narg .gt. 1) THEN
        CALL getarg(1,arg)
        IF(arg(1:2) .eq. '-t') THEN
          CALL getarg(2,arg)
          read(arg,'(i4,i2,i2,i2,i2)',iostat=istat)
     :         iyrreq,imoreq,idareq,ihrreq,iminreq
          IF(istat .eq. 0) THEN
            cmdtime=.true.
            needtime=.false.
            WRITE(6,'(a,i4.4,i2.2,i2.2,i2.2,i2.2)')
     :        'Read time from command line: ',
     :        iyrreq,imoreq,idareq,ihrreq,iminreq
            CALL mkitime(iyrreq,imoreq,idareq,ihrreq,iminreq,isec,
     :                     itime)
            CALL ictime(iyrreq,imoreq,idareq,ihrreq,iminreq,ctime,istat)
            print *, ' ictime output ctime: ',ctime
          END IF
        END IF
      ELSE IF(narg .eq. 1) THEN
        CALL getarg(1,arg)
        IF(arg(1:1).ne." ") THEN
          read(arg,*) idelt
          idelt=idelt*3600
        END IF
      END IF
C
C     Fill all data with missing flag
C
      DO 10 ksta=1,maxsta
        t(ksta)=flagd
        td(ksta)=flagd
        dd(ksta)=flagd
        ff(ksta)=flagd
        ddg(ksta)=flagd
        ffg(ksta)=flagd
        pstn(ksta)=flagd
        pmsl(ksta)=flagd
        alt(ksta)=flagd
c
        kloud(ksta)=0
        ceil(ksta)=flagd
        lowcld(ksta)=flagd
        cover(ksta)=flagd
        vis(ksta)=flagd
        rad(ksta)=flagd
        idp3(ksta)=iflagd
        DO 9 ii=1,maxcld
          emv(ii,ksta)=' '
          amt(ii,ksta)='    '
          cldhgt(ii,ksta)=flagd
   9    CONTINUE
  10  CONTINUE
C
C  Get the time of the obs to use
C
      IF(needtime) THEN
        itrunc=3600
        CALL GETIME(ctime,itime,iyrreq,imoreq,idareq,ihrreq,iminreq,
     :            idelt,itrunc)
      END IF
      write(lsoname,830) iyrreq,imoreq,idareq,ihrreq,iminreq
  830 FORMAT(i4.4,i2.2,i2.2,i2.2,i2.2,'.lso')
      print *, ' testing lsoname: ',lsoname
c
c     Open LSO file.
c
      lunit=9
      open(lunit,file=lsoname,status='unknown')
C
C     Get this hour's mesonet data
C
C     Subtract 5 mins to begin at 55 min.
C     Second pass will fill with data at 50 min.
C
      ntry=0
      jsta=0
      nobsg=0
      IF(cmdtime) THEN
        ibgn=1
        iend=1
      ELSE
        ibgn=1
        iend=npass
      END IF
      DO 50 mesopass=ibgn,iend
   20   CONTINUE
        ntry=ntry+1
        IF(cmdtime) THEN
          itimex=itime
        ELSE
          itimex=itime+ndelta5(ntry)*d5min
        END IF
        CALL cvtitoi(itimex,
     :     iyrrq,imorq,idarq,ihrrq,iminrq,isec)
        write(6,835) iyrrq,imorq,idarq,ihrrq,iminrq
  835   FORMAT(/'  Getting mesonet data for ',i4.4,4(i2.2))
C
        CALL READMESO(maxsta,iyrrq,imorq,idarq,ihrrq,iminrq,
     :              iyr,imo,ida,ihr,imin,
     :              stnin,nikname,
     :              slatin,slonin,selvin,
     :              statime,radin,tin,tdin,rh,ta9m,
     :              ddin,ffin,ffgin,avgspd,wdsd,wssd,ws2m,
     :              pstnin,altmeso,precip,nmesord,
     :              stareq,ilatest,oneonly,verbose,istat)
        print *, ' back from readmeso istat: ',istat
        IF(istat < 0 ) THEN
          CALL READMDF(maxsta,iyrrq,imorq,idarq,ihrrq,iminrq,
     :              iyr,imo,ida,ihr,imin,
     :              stnin,nikname,
     :              slatin,slonin,selvin,
     :              statime,radin,tin,tdin,rh,ta9m,
     :              ddin,ffin,ffgin,avgspd,wdsd,wssd,ws2m,
     :              pstnin,altmeso,precip,nmesord,
     :              stareq,ilatest,oneonly,verbose,istat)
        END IF
c
        print *, ' Number of mesonet stations read-in: ',nmesord

        IF(nmesord.ge.minmeso) THEN
c
          DO 45 msta=1,nmesord
            IF(mesopass.gt.1) THEN    ! check for dupes on second pass
              DO 30 ista=1,jsta
                IF(stnin(msta).eq.stn(ista)) GO TO 45
  30          CONTINUE
            END IF
            jsta=jsta+1
            nobsg=nobsg+1
            tem=tctotf( tin(msta))
            dew=tctotf(tdin(msta))
            pres=flagd
            CALL LIMCHK(pres,tem,dew,altmeso(msta),flagd)
c
            obstype(jsta)='MESO'
            read(statime(msta),840) obstime(jsta)
 840        FORMAT(6x,i4)
            stn(jsta)=stnin(msta)
            slat(jsta)=slatin(msta)
            slon(jsta)=slonin(msta)
            selv(jsta)=selvin(msta)
            wx(jsta)='        '
            t(jsta)=tem
            td(jsta)=dew
            IF(ddin(msta).ge. 0.) THEN
              dd(jsta)=ddin(msta)
              IF(ffin(msta).gt. 0.) THEN
                ff(jsta)=ffin(msta)*mstokts
                IF(ffgin(msta).gt. 0.) THEN
                  ddg(jsta)=ddin(msta)
                  ffg(jsta)=ffgin(msta)*mstokts
                END IF
              END IF
            END IF
            IF(pstnin(msta).gt. 0.) pstn(jsta)=pstnin(msta)
            IF( radin(msta).ge. 0.) rad(jsta)=radin(msta)
  45      CONTINUE  ! station loop
        END IF   ! nmesord.ge.minmeso
        print *, ' ntry,n-meso-stas: ',ntry,jsta
        IF(ntry.ge.maxtry) GO TO 51      ! escape pass loop
        IF(nmesord.lt.minmeso) GO TO 20
  50  CONTINUE   ! mesopass loop
  51  CONTINUE
      nmeso=jsta
C
C  Now get West Texas Mesonet data and append it to the mesonet
C  data.
C
       WRITE(6,842) iyrreq,imoreq,idareq,ihrreq,iminreq
 842   FORMAT(/'  Getting West Texas Mesonet data for ',i4.4,4(i2.2))

       CALL readwtx(maxsta,iyrreq,imoreq,idareq,ihrreq,iminreq,
     :          iyr,imo,ida,ihr,imin,
     :          stnin,slatin,slonin,selvin,
     :          statime,tin,tdin,
     :          ddin,ffin,ffgin,altmeso,nmesord,istat)

      IF(nmesord .gt. 0) THEN
        DO 55 msta=1,nmesord
          jsta=jsta+1
          nobsg=nobsg+1
          tem=tctotf( tin(msta))
          dew=tctotf(tdin(msta))
          pres=flagd
          CALL LIMCHK(pres,tem,dew,altmeso(msta),flagd)
c
          obstype(jsta)='WTXMN'
          read(statime(msta),840) obstime(jsta)
          stn(jsta)=stnin(msta)
          slat(jsta)=slatin(msta)
          slon(jsta)=slonin(msta)
          selv(jsta)=selvin(msta)
          wx(jsta)='        '
          t(jsta)=tem
          td(jsta)=dew
          IF(ddin(msta).ge. 0.) THEN
            dd(jsta)=ddin(msta)
            IF(ffin(msta).gt. 0.) THEN
              ff(jsta)=ffin(msta)
            END IF
          END IF
          IF(altmeso(msta).gt. 0.) 
     :         alt(jsta)=intomb*altmeso(msta)
          pstn(jsta)=flagd
          rad(jsta)=flagd
  55    CONTINUE  ! station loop
      END IF   ! nmesord.ge.0
      nmeso=jsta
      nobsg=nmeso
      nobsb=nmeso
      WRITE(6,'(a,i6)') ' After W. Texas, nmeso = ',nmeso
      write(6,*) ' Total observations: ',nobsb
C
C  Write to .LSO file
C
c
c.....  File open...first write the header.
c
      print *, ' ctime = ',ctime
      write(lunit,800) ctime,              ! data time
     &               nmeso,          ! # of mesonet stations
     &               nmesord,        ! total # mesonet stations possible
     &               nsaog,           ! # of saos in the laps grid
     &               nsaord,       ! total # of saos possible in laps grid
     &               nsao,           ! # of saos in the box
     &               nsaord,       ! total # of saos possible in the box
     &               nobsg,           ! # of obs in the laps grid
     &          (nmesord+nsaord),  ! total # of obs psbl in the laps grid
     &               nobsb,           ! # of obs in the box
     &               maxsta        ! total # of obs possible in the box
800   format(1x,a23,1x,10(1x,i4))
c
c.....  Error trapping for too many stations for array size.
c
      IF (nobsb .gt. maxsta) THEN
          print 880, maxsta,nobsb,ctime
880       format(' +++ ERROR in READ_SURFACE_OBS: maxsta = ',i8,/,
     & ' but there are ',i8,' stations in the ',a24,' obs file.',/)
          print *,'    Increase the value of "maxsta" and try again.'
          istatus = -2
          STOP
      END IF
c
c.....  Now write the station data.
c
      DO 100 k=1,nobsb
        write(lunit,801) stn(k),slat(k),slon(k),selv(k),
     &                   obstype(k),obstime(k),wx(k)
801     format(1x,a5,f6.2,1x,f7.2,1x,f5.0,1x,a8,1x,i4.4,1x,a8)
c
        write(lunit,803) t(k),td(k),dd(k),ff(k),ddg(k),ffg(k),pstn(k),
     &                pmsl(k),alt(k)
803     format(4x,2(f6.1,1x),4(f5.0,1x),3(f6.1,1x))
c
        write(lunit,805) kloud(k),ceil(k),lowcld(k),cover(k),
     &                   vis(k),rad(k),idp3(k)
805       format(4x,i2,2(1x,f7.1),1x,f5.1,1x,f7.3,1x,f6.1,1x,i4)
c
c.....  Write the cloud data if we have any.
c
          IF (kloud(k) .gt. 0) THEN
            DO 75 ii=1,kloud(k)
              write(lunit,807) 
     :          emv(ii,k),amt(ii,k),cldhgt(ii,k)
807           format(5x,a1,1x,a4,1x,f7.1)
 75         CONTINUE !ii
          END IF
c
  100 CONTINUE
c
      close(lunit)
c
c
c     Open file to hold the name of the lso file to be
c     used by various automated procedures.
c
      open(11,file='lsoname.last',status='unknown')
      write(11,'(a)') lsoname
      close(11)
c
      STOP
      END
c
      SUBROUTINE LIMCHK(pres,temp,dewpt,alt,flagd)
      IMPLICIT NONE
      REAL pres,temp,dewpt,alt,flagd
      real intomb
      parameter (intomb=(1013.25/29.92))
      real atem
      IF(pres.LT.890. .OR. pres.GT.1070.) pres=flagd
      IF(temp.LT.-40. .OR. temp.GT.130.) temp=flagd
      IF(dewpt.LT.-40..OR. dewpt.GT.temp) dewpt=flagd
      atem=alt*intomb
      IF(atem.LT.890. .OR. atem.GT.1070.) alt=flagd
      RETURN
      END
c
      FUNCTION tctotf(ctemp)
      IMPLICIT NONE
      real ctemp
      real tctotf
      IF(ctemp.gt.-50.) THEN
        tctotf=1.8*ctemp + 32.
      ELSE
        tctotf=ctemp
      END IF
      RETURN
      END
