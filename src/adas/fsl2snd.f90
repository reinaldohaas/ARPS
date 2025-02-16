!-----------------------------------------------------------------------
!
!The attached program will process a file containing one or more soundings from 
!the FSL (GSD) online radiosonde database and translate the soundings into 
!file(s) that can be used in ADAS or 3DVAR.   For files with multiple data 
!times, the soundings will be sorted by data time and output to the appropriate 
!file, each named according to the date and time.  The program can also use a 
!station file to verify the station locations with the latest global station 
!list information from NOAA – use of this station file is optional.   
!
!You can get that file from http://weather.noaa.gov/data/nsd_cccc.gz
! 
!To set-up the program, just compile it with an f90 compiler.  No other 
!libraries are required.
!
!To use the program just start it and it will prompt for the name of the 
!ASCII file saved from the web download of the data from
!
!   http://raob.fsl.noaa.gov 
!
!When using the web interface select “FSL Format (ASCII)” and wind units of 
!“kts” (the default selections).
!
!To use the output in ADAS or ARPS 3DVAR put the name of the output file 
!in the section of the arps.input file that begins with nuafil, using the 
!variable uafname.
!
!-Keith (05.24.2007)
!
!-----------------------------------------------------------------------

  PROGRAM fsl2snd
!
! Keith Brewster
! CAPS, University of Oklahoma
! May, 2007
!
  IMPLICIT NONE
!
! Header data
!
  INTEGER :: lintyp
  INTEGER :: hour,day,year,minute
  CHARACTER(LEN=3) :: month
  INTEGER :: wban,wmo
  REAL :: rlat,rlon,selev
  INTEGER :: ielev,itime
  INTEGER :: hydro,mxwd,tropl,lines,tindex,source
  CHARACTER(LEN=4) :: staid
  INTEGER :: sonde
  CHARACTER(LEN=2) :: wsunits
  CHARACTER(LEN=1) :: ns,ew

  REAL, PARAMETER :: kts2ms=0.514444
  REAL :: dtr,spcvt,dtdz,dtddz,dudz,dvdz
  REAL :: tempk,dewpk
  INTEGER :: iline,nlines,nlevs,i,j,k,imon,lunit,ntime
  INTEGER :: kbot,ktop,kk,kn
  INTEGER :: idata(6)
  LOGICAL :: fndtop,newtime

  INTEGER, PARAMETER :: maxtime=20
  CHARACTER(LEN=126) :: fname
  CHARACTER(LEN=126) :: outsnd
  CHARACTER(LEN=12) :: timestr
  CHARACTER(LEN=12) :: filtime(maxtime)
  INTEGER :: filunit(maxtime)
  CHARACTER(LEN=3) :: chmon(12) =                                              &
    (/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)
!
! Output data
!
  INTEGER, parameter ::  maxlev=400
  REAL :: press(maxlev)
  REAL :: height(maxlev)
  REAL :: temp(maxlev)
  REAL :: dewpt(maxlev)
  REAL :: wdir(maxlev)
  REAL :: wspd(maxlev)
  REAL :: u(maxlev)
  REAL :: v(maxlev)
!
! Station data table variables
!
  INTEGER, PARAMETER :: maxtable=10000
  CHARACTER (LEN=180) :: statbl(maxtable)
  CHARACTER (LEN=4) :: staidtbl(maxtable)
  INTEGER :: wmotbl(maxtable)
  INTEGER, PARAMETER :: maxsemi=100
  INTEGER :: slocs(maxsemi)
  CHARACTER(LEN=1) :: cdir
  REAL :: slat,slon
  INTEGER :: ilat1,ilat2,ilat3,ilon1,ilon2,ilon3
  INTEGER :: ista,jsta,nsta,istat,wmogrp,wmosta,nsemi
  LOGICAL :: EOF
!
  nsta=0
  minute=0
  EOF=.false.
  dtr=atan(1.)/45.
  DO itime=1,maxtime
    filunit(itime)=20+itime
  END DO
!
  print *, ' Opening NOAA station list file'
  OPEN(19,file='nsd_cccc',iostat=istat,status='old')
  IF(istat == 0) THEN
  EOF=.false.
  DO jsta=1,maxtable
    READ(19,'(a180)',iostat=istat) statbl(jsta)
    IF(istat < 0 ) EOF = .true.
    IF(istat /= 0 ) EXIT
    staidtbl(jsta)=statbl(jsta)(1:4)
    READ(statbl(jsta)(6:11),'(i2,1x,i3)',iostat=istat) wmogrp,wmosta
    IF(istat == 0) THEN
      wmotbl(jsta)=(wmogrp*1000)+wmosta
    ELSE
      wmotbl(jsta)=-999
    END IF
!   print *, ' jsta: ',jsta,'  id:',staidtbl(jsta),' wmo:',wmotbl(jsta)
  END DO
  nsta=jsta-1
  IF(EOF) THEN
    WRITE(6,'(a)') ' End of file reached'
  ELSE
    WRITE(6,'(a/a,i6)') &
      ' Error before end of file or ran out of allocated space', &
      ' Space allocated is',maxtable
  END IF
  WRITE(6,'(a,i6,a)') ' Read',nsta,' stations from NOAA file'
  CLOSE(19)
  ELSE
    nsta=0
    WRITE(6,'(a)') ' Error opening NOAA station table file, nsd_cccc'
    WRITE(6,'(a)') ' Obtain file from'
    WRITE(6,'(a)') '   http://weather.noaa.gov/data/nsd_cccc.gz'
    WRITE(6,'(a)') '   and uncompress it.'
    WRITE(6,'(a)') ' Meanwhile, continuing using FSL assigned lat,lons'
  END IF
!
  write(6,'(a)') ' Input data file name'
  read(5,'(a126)') fname
  write(6,'(a,a)') 'Opening ',fname
  open(11,file=fname,status='old')
  DO ista=1,99999
!
! Read header lines
!
  read(11,'(3i7,6x,a3,i8)',iostat=istat) lintyp,hour,day,month,year
  IF(istat /= 0) EXIT
  read(11,'(3i7,f7.2,a1,f6.2,a1,i6,i7)',iostat=istat)  &
           lintyp,wban,wmo,rlat,ns,rlon,ew,ielev,itime
  IF(istat /= 0) EXIT
  read(11,'(7i7)',iostat=istat) lintyp,hydro,mxwd,tropl,lines,tindex,source
  IF(istat /= 0) EXIT
  read(11,'(i7,10x,a4,14x,i7,5x,a2)',iostat=istat) lintyp,staid,sonde,wsunits
  IF(istat /= 0) EXIT

  IF(wsunits.eq.'kt') THEN
    spcvt=kts2ms
  ELSE
    spcvt=1.0
  END IF
  IF(ns .eq. 'S') rlat=-rlat
  IF(ew .eq. 'W') rlon=-rlon
  selev=float(ielev)
  DO imon=1,11
    IF(chmon(imon) == month) EXIT
  END DO
!
!  Open output file
!
  write(timestr,'(i4.4,i2.2,i2.2,i2.2,i2.2)') year,imon,day,hour,minute
  IF(ista == 1) THEN
    write(outsnd,'(a,a,a)') 'raobfsl',timestr,'.snd'
    filtime(1)=timestr
    lunit=filunit(1)
    open(lunit,file=outsnd,status='unknown')
    ntime=1
  ELSE
    newtime=.true.
    DO itime=1,ntime
      IF(filtime(itime) == timestr) THEN
        newtime=.false.
        lunit=filunit(itime)
        EXIT
      END IF
    END DO 
    IF(newtime) THEN
      IF(ntime < maxtime) THEN
        ntime=ntime+1
        lunit=filunit(ntime)
        filtime(ntime)=timestr
        write(outsnd,'(a,a,a)') 'raobfsl',timestr,'.snd'
        open(lunit,file=outsnd,status='unknown')
      ELSE
        write(6,'(a,i4,a,i4)') ' Ntime',(ntime+1),' exceeds maxtime',maxtime
        EXIT
      END IF
    END IF
  END IF
!
! Get more precise station location info from the data table
!
  slat=-999.
  slon=-999.
  print *, ' Searching loc table for ',staid,wmo
  DO jsta=1,nsta
    IF((staidtbl(jsta) == staid) .OR. (wmotbl(jsta) == wmo)) THEN
      print *, ' Found match at ',jsta
      print *, ' Table values: ',staidtbl(jsta),wmotbl(jsta)
      CALL semilocs(maxsemi,statbl(jsta),nsemi,slocs)
      print *, ' nsemi=',nsemi
      ilat3=0
      cdir='N'
      IF((slocs(8)-slocs(7)) > 8) THEN
        READ(statbl(jsta)((slocs(7)+1):(slocs(7)+9)),'(i2,1x,i2,1x,i2,a1)',iostat=istat) &
                 ilat1,ilat2,ilat3,cdir
        print *, 'Lat read2: ',ilat1,ilat2,istat
        IF(istat /= 0) EXIT
      ELSE
        READ(statbl(jsta)((slocs(7)+1):(slocs(7)+6)),'(i2,1x,i2,a1)',iostat=istat) &
                 ilat1,ilat2,cdir
        print *, 'Lat read2: ',ilat1,ilat2,istat
        IF(istat /= 0) EXIT
      END IF
      slat=float(ilat1)+(float(ilat2)/60.)+(float(ilat3)/3600.)
      IF(cdir == 'S') slat=-slat
      ilon3=0
      cdir='W'
      IF((slocs(9)-slocs(8)) > 9) THEN
        READ(statbl(jsta)((slocs(8)+1):(slocs(8)+10)),'(i3,1x,i2,1x,i2,a1)',iostat=istat) &
           ilon1,ilon2,ilon3,cdir
        print *, 'Lon read1: ',ilon1,ilon2,istat
        IF(istat /= 0) EXIT
      ELSE
        READ(statbl(jsta)((slocs(8)+1):(slocs(8)+7)),'(i3,1x,i2,a1)',iostat=istat) &
             ilon1,ilon2,cdir
        print *, 'Lon read2: ',ilon1,ilon2,istat
        IF(istat /= 0) EXIT
      END IF
      slon=float(ilon1)+(float(ilon2)/60.)+(float(ilon3)/3600.)
      IF(cdir == 'W') slon=-slon
      EXIT
    END IF
  END DO
  print *, ' rlat: ',rlat,'  rlon:',rlon
  print *, ' slat: ',slat,'  slon:',slon
  IF(slat > -900.) rlat=slat
  IF(slon > -900.) rlon=slon
!
! Read Sonde data
!
  k=0
  nlines=lines-4
  DO iline=1,nlines
    read(11,'(7i7)',iostat=istat) lintyp,(idata(j),j=1,6)
    IF(istat /= 0) EXIT
!    print *, lintyp,(idata(j),j=1,6)
    IF(idata(2).ge.ielev .AND. idata(2).lt.99990) THEN
      k=k+1
      IF(idata(1).lt.99990) THEN
        press(k)=0.1*idata(1)
      ELSE
        press(k)=99999.
      END IF
      height(k)=float(idata(2))
      IF(idata(3).lt.99990) THEN
        temp(k)=idata(3)*0.1
      ELSE
        temp(k)=99999.
      END IF
      IF(idata(4).lt.9990) THEN
        dewpt(k)=idata(4)*0.1
      ELSE
        dewpt(k)=99999.
      END IF
      IF(idata(5).lt.99990 .and. idata(6).lt.200) THEN
        wdir(k)=float(idata(5))
        wspd(k)=spcvt*idata(6)
        u(k)=-wspd(k)*sin(dtr*wdir(k))
        v(k)=-wspd(k)*cos(dtr*wdir(k))
      ELSE
        u(k)=99999.
        v(k)=99999.
      END IF
    END IF
  END DO
  nlevs=k
!
!     Interpolate missing temperatures
!
  DO k=2,nlevs
    IF(temp(k).gt.50.) THEN
      kbot=k-1
      fndtop=.false.
      DO kn=k+1,nlevs
        IF(temp(kn).le.50.) THEN
          fndtop=.true.
          ktop=kn
          EXIT
        END IF
      END DO
      IF(fndtop) THEN
        dtdz=(temp(ktop)-temp(kbot))/                        &
             (height(ktop)-height(kbot))
        DO kk=k,ktop-1
          temp(kk)=temp(kbot)+dtdz*(height(kk)-height(kbot))
        END DO
      ELSE
        EXIT
      END IF
    END IF
  END DO
!
! Interpolate missing dew points
!
  DO k=2,nlevs
    IF(dewpt(k).gt.40.) THEN
      kbot=k-1
      fndtop=.false.
      DO kn=k+1,nlevs
        IF(dewpt(kn).le.40.) THEN
          fndtop=.true.
          ktop=kn
          EXIT
        END IF
      END DO
      IF(fndtop) THEN
        dtddz=(dewpt(ktop)-dewpt(kbot))/                        &
              (height(ktop)-height(kbot))
        DO kk=k,ktop-1
          dewpt(kk)=dewpt(kbot)+dtddz*(height(kk)-height(kbot))
        END DO
      ELSE
        EXIT
      END IF
    END IF
  END DO
!
!  Interpolate missing winds
!
  DO k=2,nlevs
    IF(u(k).gt.300.) THEN
      kbot=k-1
      fndtop=.false.
      DO kn=k+1,nlevs
        IF(u(kn).le.300.) THEN
          fndtop=.true.
          ktop=kn
          EXIT
        END IF
      END DO
      IF(fndtop) THEN
        dudz=(u(ktop)-u(kbot))/                           &
             (height(ktop)-height(kbot))
        dvdz=(v(ktop)-v(kbot))/                           &
             (height(ktop)-height(kbot))
        DO kk=k,ktop-1
          u(kk)=u(kbot)+dudz*(height(kk)-height(kbot))
          v(kk)=v(kbot)+dvdz*(height(kk)-height(kbot))
          CALL get_ddff(u(kk),v(kk),wdir(kk),wspd(kk))
        END DO
      ELSE
        EXIT
      END IF
    END IF
  END DO
!
! Write sounding
! Note radiosonde dewpts are always wrt liquid.
!
  write(lunit,'(i12,i12,f11.4,f15.4,f15.0,5x,a4)') &
       wmo,nlevs,rlat,rlon,selev,staid
  DO k = 1,nlevs
    write(lunit,'(f10.1,3f10.2,f10.1,f10.2)') &
       height(k),press(k),temp(k),dewpt(k),wdir(k),wspd(k)
  END DO
  END DO  ! big station loop
  DO itime=1,ntime
    CLOSE(filunit(itime))
  END DO
  WRITE(6,'(a,i5,a,i4,a)') &
   ' Processed ',(ista-1),' sounding(s) at ',ntime,' separate time(s)'
  WRITE(6,'(a)') ' fsl2snd successful completion'
  STOP
  END PROGRAM fsl2snd
!
  SUBROUTINE semilocs(maxsemi,string,nsemi,slocs)
! 
! Parse a string and report the location of all semi-colons
! Keith Brewster, CAPS
!
  IMPLICIT NONE
  INTEGER :: maxsemi
  CHARACTER (LEN=*) string
  INTEGER :: nsemi
  INTEGER :: slocs(maxsemi)
                                                                                
  INTEGER :: iloc,i
  INTEGER :: slen
  nsemi=0
  DO iloc=1,maxsemi
    slocs(iloc)=-99
  END DO
  slen=LEN(TRIM(string))
  DO i=1,slen
    IF(string(i:i) == ';') THEN
      nsemi=nsemi+1
      slocs(nsemi)=i
    END IF
  END DO
  END SUBROUTINE semilocs

!
  SUBROUTINE GET_DDFF(U,V,DD,FF)
  IMPLICIT NONE
  REAL :: U,V,DD,FF
  REAL, PARAMETER :: RAD2D=57.29577951
  REAL, PARAMETER :: SPVAL=9999.
  REAL, PARAMETER :: MIS_VAL=99999.0
  IF(U.LT.SPVAL .AND. V.LT.SPVAL) THEN
    FF = SQRT((U*U + V*V))
    IF(FF.NE.0.) THEN
      DD = RAD2D*ATAN2(U,V)
      DD = DD+180.
      IF (DD.GT.360.) DD=DD-360.
    ELSE
      DD=0.
    END IF
  ELSE
    DD = MIS_VAL
    FF = MIS_VAL
  END IF
  RETURN
  END SUBROUTINE GET_DDFF
