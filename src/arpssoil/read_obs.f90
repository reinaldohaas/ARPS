SUBROUTINE read_obs(c9time,                                             &
           xtnam,xtlat,xtlon,xtelev,xts,k2_sta)
!
!  Reads-in surface (soil and radiation) data from mesonet
!  files.
!
  IMPLICIT NONE
  INTEGER :: maxsta,jsta,ksta,k_sta,num,nmesotable,nsta,j,nvar_ob
  INTEGER :: ista,k2_sta,i,day,day1
  PARAMETER (maxsta=500,nvar_ob=3)
  CHARACTER (LEN=4) :: snam(maxsta),stnam(maxsta),xtnam(maxsta)
  CHARACTER (LEN=4) :: srtnam(maxsta)
  CHARACTER (LEN=9) :: c9time
  CHARACTER (LEN=100) :: xs
  CHARACTER (LEN=100) :: s
  REAL :: slat(maxsta),slon(maxsta),selev(maxsta)
  REAL :: ts05(maxsta),ts30(maxsta),srad(maxsta)
  REAL :: xtlat(maxsta),xtlon(maxsta),xtelev(maxsta)
  REAL :: xts(maxsta,nvar_ob)
  CHARACTER (LEN=10) :: temp_day1
  CHARACTER (LEN=3) :: temp_day
!
  CHARACTER (LEN=80) :: filenam
  INTEGER :: inelev,inrad,knt
  REAL :: s1,s2,s3,s4,sum

  PRINT*,'c9time:',c9time(1:9)
!  Initialize xts
!
  DO j=1,nvar_ob
    DO i=1,maxsta
      xts(i,j)=-999.
    END DO
  END DO
!
!  Get mesonet stations from table
!
!  open(3,file='/laps/vortex/source/graf/meso.list',status='old')
  OPEN(3,FILE='/arps/work/craft94/olaps94.data/mesonet_station_list'    &
       ,STATUS='old')

  DO jsta=1,maxsta
    READ(3,801,END=6)                                                   &
         snam(jsta),num,slat(jsta),slon(jsta),inelev,xs
    selev(jsta)=FLOAT(inelev)
    801    FORMAT(a4,i5,f7.2,f7.2,i5,a34)
    slon(jsta)=-slon(jsta)
  END DO
  6  CONTINUE
  CLOSE(3)
  nmesotable=jsta-1
!
!  Read in "d15" file containing mesonet
!  "supplemental" data.
!

  filenam=c9time(1:9)//'.d15'

  PRINT*,'filename:',filenam
  OPEN(7,FILE=filenam,STATUS='old')
  READ(7,*)
  DO ksta=1,maxsta
    READ(7,815,END=9)stnam(ksta),s,ts05(ksta),s1,ts30(ksta),s3,s4
    815  FORMAT(1X,a4,a28,5(f8.2))
  END DO
  9  CONTINUE
  CLOSE(7)
  nsta=ksta-1
!
  k_sta=0
  DO ksta=1,nsta
    DO jsta=1,nmesotable
      IF(stnam(ksta)(1:4) == snam(jsta)(1:4))THEN
        k_sta=k_sta+1
        xtnam(k_sta)(1:4)=snam(jsta)(1:4)
        xtlat(k_sta) =slat(jsta)
        xtlon(k_sta) =slon(jsta)
        xtelev(k_sta)=selev(jsta)
        xts(k_sta,1) =ts05(ksta)
        xts(k_sta,2) =ts30(ksta)
!         print*,xtnam(k_sta)(1:4),xts(k_sta,1),xts(k_sta,2)
        EXIT
      END IF
    END DO
  END DO
!
!  Read SRAD from mesonet *d05 file.
!
  filenam= c9time(1:9)//'.d05'

  PRINT*,'SRAD filename:',filenam
  OPEN(7,FILE=filenam,STATUS='old')
  READ(7,*)
  knt=0
  sum=0.
  DO ksta=1,maxsta
    READ(7,816,END=14)srtnam(ksta),s,inrad,xs
!    print*,srtnam(ksta),inrad
    srad(ksta)=FLOAT(inrad)
    IF(inrad > -1) THEN
      knt=knt+1
      sum=sum+srad(ksta)
    END IF
    816    FORMAT(1X,a4,a84,i6,a14)
  END DO
  14     CONTINUE
  CLOSE(7)
  nsta=ksta-1
!
!  Fast QC: Often srad is reported zero due to an error.
!  If the mean srad is greater than 200 and ob is exactly
!  zero, assume its bogus.
!
  IF(knt > 0) THEN
    sum=sum/FLOAT(knt)
    IF(sum >= 200.) THEN
      DO ksta=1,nsta
        IF(srad(ksta) == 0.) srad(ksta)=-999.
      END DO
    END IF
  END IF

  k2_sta=k_sta
  DO ksta=1,nsta
    DO jsta=1,k_sta
      IF(srtnam(ksta)(1:4) == xtnam(jsta)(1:4))THEN
        xts(jsta,3)=srad(ksta)
        GOTO 17
      END IF
    END DO

    DO ista=1,nmesotable
      IF(srtnam(ksta)(1:4) == snam(ista)(1:4)) THEN
        k2_sta=k2_sta+1
        xtnam(k2_sta)(1:4)=snam(ista)(1:4)
        xtlat(k2_sta) =slat(ista)
        xtlon(k2_sta) =slon(ista)
        xtelev(k2_sta)=selev(ista)
        xts(k2_sta,3)=srad(ksta)
        GOTO 17
      END IF
    END DO
17 CONTINUE
  END DO

! 3)write out station by matching its name to the names in the station table
!  do j= 1,k2_sta
!    write(12,61) xtnam(j)(1:4),xtlat(j),xtlon(j),xtelev(j),
!    .(xts(j,i),i=1,3)
!  end do
!61   format(a4,6(f10.4))
!
  RETURN
END SUBROUTINE read_obs
