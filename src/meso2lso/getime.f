c     program test
c     implicit none
c     CHARACTER*23 ctime
c     INTEGER itime
c     INTEGER iyr,imon,iday,ihour,imin,ibatch
c     call GETIME(ctime,itime,iyr,imon,iday,ihour,imin,ibatch)
c     STOP
c     END

      SUBROUTINE GETIME(ctime,itime,
     :                  iyr,imon,iday,ihour,imin,idelt,itrunc)
      IMPLICIT NONE
C
C     This is the Unix version for SunOS.
C
C     Keith Brewster
C     March, 1995
C
C  Arguments
C
      CHARACTER*23 ctime
      INTEGER itime
      INTEGER iyr,imon,iday,ihour,imin,idelt,itrunc
C
C  Intrinsic Functions
C
      INTEGER time
C
C  Months
C
      CHARACTER*3 cmon(12)
      DATA cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     :           'JUL','AUG','SEP','OCT','NOV','DEC'/
C
C  Misc internal variables
C
      CHARACTER*23 reqtime
      INTEGER isec,ierr,tty
      INTEGER isatty

      tty=0
C     print *, ' isatty(6):',isatty(6)
C
      IF( tty.eq.0 ) THEN
        print *, ' Using current time for batch process'
        reqtime(1:4)='    '    ! what is checked later for using current time
      ELSE
   6    write(6,803)
 803    FORMAT
     +('  Enter UTC time DD-MON-YYYY HH or <return> for latest:')
        read(5,804) reqtime
 804    FORMAT(A23)
      END IF
C
      IF(reqtime(1:4).EQ.'    ') THEN
C
C       Obtain current system time using time function
C       The time returned is in UTC.
C
        itime=time()
c
c  If this is a batch job, subtract an hour for data deliverly delay
c
        itime=itime+idelt
c
c  Truncate time according to truncation time in argument list
c
        IF(itrunc.gt.0) itime=(itime/itrunc)*itrunc

        CALL ITIM2CI(itime,iyr,imon,iday,ihour,imin,isec,ctime)
        print *, 'yr,mon,day,hour,imin = ',iyr,imon,iday,ihour,imin
        print *, ' Requested UTC time: ',ctime
      ELSE
        reqtime(19:23)='00.00'
        CALL str_upcase(reqtime,reqtime,23)
        call citime(reqtime,iyr,imon,iday,ihour,imin,ierr)
c
        print *, 'yr,mon,day,hour,min = ',iyr,imon,iday,ihour,imin
        isec=0
        call mkitime(iyr,imon,iday,ihour,imin,isec,itime)
        print *, ' itime = ',itime
        print *, ' Requested UTC time: ',reqtime
        ctime=reqtime

        CALL ITIM2CI(itime,iyr,imon,iday,ihour,imin,isec,ctime)
        print *, 'yr,mon,day,hour,min = ',iyr,imon,iday,ihour,imin
        print *, ' Requested UTC time: ',ctime
      END IF
      RETURN
      END
c
      SUBROUTINE ITIM2CH(i4time,ctime)
c
c     Convert i4time to character time
c
      integer i4time
      character*23 ctime
c
      integer tarray(9)
      integer iyr,imon,iday,ihour,imin,isec
C
C  Months
C
      CHARACTER*3 cmon(12)
      DATA cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     :           'JUL','AUG','SEP','OCT','NOV','DEC'/
c
      call gmtime(i4time,tarray)
      isec=tarray(1)
      imin=tarray(2)
      ihour=tarray(3)
      iday=tarray(4)
      imon=tarray(5)+1
      iyr=tarray(6)+1900
      IF(iyr.lt.1970) iyr=iyr+100
      write(ctime,810) iday,cmon(imon),iyr,ihour,imin,isec
 810  FORMAT(i2.2,'-',a3,'-',i4.4,' ',i2.2,':',i2.2,':',i2.2,'.00')
      RETURN
      END
c
      SUBROUTINE ITIM2CI(i4time,iyr,imon,iday,ihour,imin,isec,ctime)
c
c     Convert i4time to character time, also return individual integers
c
      integer i4time
      character*23 ctime
      integer iyr,imon,iday,ihour,imin,isec
c
      integer tarray(9)
C
C  Months
C
      CHARACTER*3 cmon(12)
      DATA cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     :           'JUL','AUG','SEP','OCT','NOV','DEC'/
c
      call gmtime(i4time,tarray)
      isec=tarray(1)
      imin=tarray(2)
      ihour=tarray(3)
      iday=tarray(4)
      imon=tarray(5)+1
      iyr=tarray(6)+1900
      IF(iyr.lt.1970) iyr=iyr+100
      write(ctime,810) iday,cmon(imon),iyr,ihour,imin,isec
 810  FORMAT(i2.2,'-',a3,'-',i4.4,' ',i2.2,':',i2.2,':',i2.2,'.00')
      RETURN
      END
c
      SUBROUTINE ITIM2IL(i4time,iyr,imon,iday,ihour,imin,isec)
c
c     Convert i4time time into individual integers for each 
c     time element
c
      integer i4time
      integer iyr,imon,iday,ihour,imin,isec
c
      integer tarray(9)
c
      call gmtime(i4time,tarray)
      isec=tarray(1)
      imin=tarray(2)
      ihour=tarray(3)
      iday=tarray(4)
      imon=tarray(5)+1
      iyr=tarray(6)+1900
      IF(iyr.lt.1970) iyr=iyr+100
      RETURN
      END
