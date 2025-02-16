!

SUBROUTINE pltmap(maxpts,fname,latmap,lonmap,xmap,ymap)
  IMPLICIT NONE
  INTEGER :: maxpts
  CHARACTER (LEN=*) :: fname
  REAL :: latmap(maxpts),lonmap(maxpts)
  REAL :: xmap(maxpts),ymap(maxpts)
!
!  Misc internal variables
!
  REAL :: r2deg
  PARAMETER(r2deg=180./3.141592654)
!
  LOGICAL :: eof
  INTEGER :: kk,kkold,j,jstr,iblock,npts
!
  OPEN(17,FILE=TRIM(fname),STATUS='old')
!
!  Records look like this.  Lon,Lat in radians.
!    1 -1.84431  0.69810
!
  eof=.false.
  kkold=1
  jstr=1
  DO iblock=1,10000
    DO j=jstr,maxpts
      READ(17,*,END=101) kk,lonmap(j),latmap(j)
      IF(kk /= kkold) GO TO 110
    END DO
    GO TO 110
    101   CONTINUE
    eof=.true.
    110   CONTINUE
!
    npts=j-1
!
!    print *, '  kkold,npts: ',kkold,npts
    DO j=1,npts
      latmap(j)=r2deg*latmap(j)
      lonmap(j)=r2deg*lonmap(j)
    END DO
!    PRINT *, ' Lat(1),Lon(1) = ',latmap(1),lonmap(1)
    CALL lltoxy(npts,1,latmap,lonmap,xmap,ymap)
    CALL curve(xmap,ymap,npts)
    IF(eof) EXIT
    jstr=2
    kkold=kk
    lonmap(1)=lonmap(npts+1)
    latmap(1)=latmap(npts+1)
  END DO
  401 CONTINUE
!
  CLOSE(17)
!
  RETURN
END SUBROUTINE pltmap
