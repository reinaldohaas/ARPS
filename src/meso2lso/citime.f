      SUBROUTINE CITIME(ctime,iyr,imo,ida,ihr,imin,ierr)
C
C  Converts VAX character string time (ctime) into integer
C  form needed for Gempak subroutines
C     K. Brewster July, 1991
C
      IMPLICIT NONE
C
C  Arguments
C
      CHARACTER*23 ctime
      INTEGER iyr,imo,ida,ihr,imin,ierr
C
C  misc internal variables
C
      CHARACTER*3 strmon,monnam(12)
      DATA monnam /'JAN','FEB','MAR','APR','MAY','JUN',
     +             'JUL','AUG','SEP','OCT','NOV','DEC'/
C
      ierr=0
      read(ctime,850) ida,strmon,iyr,ihr,imin
 850  FORMAT(i2,1x,a3,1x,i4,1x,i2,1x,i2)
      DO 40 imo=1,12
        IF(monnam(imo).eq.strmon) GO TO 41
  40  CONTINUE
      print *, 'Error reading month.  Month= ',strmon
      ierr=-1
      imo=1
  41  CONTINUE
C
C     print *, ' Inside CITIME  ctime: ',ctime
C     print *, ' imo,ida,iyr ',imo,ida,iyr
C     print *, ' ihr,imin ',ihr,imin
C
      RETURN
      END
      SUBROUTINE ICTIME(iyr,imo,ida,ihr,imin,ctime,ierr)
C
C  Converts VAX character string time (ctime) into integer
C  form needed for Gempak subroutines
C     K. Brewster July, 1991
C
      IMPLICIT NONE
C
C  Arguments
C
      CHARACTER*23 ctime
      INTEGER iyr,imo,ida,ihr,imin,ierr
C
C  misc internal variables
C
      CHARACTER*3 strmon,monnam(12)
      DATA monnam /'JAN','FEB','MAR','APR','MAY','JUN',
     +             'JUL','AUG','SEP','OCT','NOV','DEC'/
C
      ierr=0
      IF(imo < 13 .and. imo > 0) THEN
        write(ctime,850) 
     +    ida,'-',monnam(imo),'-',iyr,ihr,':',imin,':00.00'
 850    FORMAT(i2.2,a1,a3,a1,i4.4,1x,i2.2,a1,i2.2,a6)
      ELSE
        ierr=-1
        strmon=monnam(1)
        write(ctime,850) 
     +    ida,'-',strmon,'-',iyr,ihr,':',imin,':00.00'
      END IF
C
C     print *, ' Inside CITIME  ctime: ',ctime
C     print *, ' imo,ida,iyr ',imo,ida,iyr
C     print *, ' ihr,imin ',ihr,imin
C
      RETURN
      END
