!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE RDPROF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rdprof(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,proffile,         &
           stnua,elevua,xua,yua,hgtua,obsua,                            &
           qualua,isrcua,nlevsua,                                       &
           rmiss,nprev,ntotal,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read ASCII file containing wind profiler data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  July, 1995
!
!  MODIFICATION HISTORY:
!  9/3/1996  Keith Brewster
!  Added full documentation.
!
!  2/16/1998 Keith Brewster
!  Added jsrc, source number to the variable list
!
!  7/14/2011 Keith Brewster
!  Added srcua to the argument list and read of src from data file.
!  isrcua is then determined from the source name read-in.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nvar       Number of variables in analysis (array dimension)
!    nzua       Maximum number of vertical levels (array dimension)
!    mxua       Maximum number of UA stations (array dimension)
!    nsrcua     Number of upper-air sources (array dimension)
!    srcua      Names of upper-air sources
!    proffile   Name of profiler data file to read
!    rmiss      Missing data fill value
!    nprev      Number of stations read previously into UA observation
!               data arrays (array index)
!
!  OUTPUT :
!
!    stnua      station name (character*4)
!    elevua     station elevation (m MSL)
!    xua        station location x-coordinate (m)
!    yua        station location y-coordinate (m)
!    hgtua      height of data (m MSL)
!    obsua      observation data
!    qualua     observation quality indicator
!    isrcua     data source index
!    nlevsua    number of levels of data
!    ntotal     number of stations
!    istatus    status indicator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar,nzua,mxua,nsrcua
  CHARACTER(LEN=8), INTENT(IN) :: srcua(nsrcua)
  CHARACTER(LEN=8), INTENT(IN) :: dfltsrc
  CHARACTER(LEN=256), INTENT(IN) :: proffile
  CHARACTER(LEN=5), INTENT(OUT) :: stnua(mxua)
  REAL, INTENT(OUT) :: elevua(mxua)
  REAL, INTENT(OUT) :: xua(mxua)
  REAL, INTENT(OUT) :: yua(mxua)
  REAL, INTENT(OUT) :: hgtua(nzua,mxua)
  REAL, INTENT(OUT) :: obsua(nvar,nzua,mxua)
  INTEGER, INTENT(OUT) :: qualua(nvar,nzua,mxua)
  INTEGER, INTENT(OUT) :: isrcua(mxua)
  INTEGER, INTENT(OUT) :: nlevsua(mxua)
  REAL, INTENT(IN) :: rmiss
  INTEGER, INTENT(IN) :: nprev
  INTEGER, INTENT(OUT) :: ntotal
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar,jlev,jend,ksta,mxprof,nprof
  INTEGER :: isrc,ilen,opstat,rdstat
  INTEGER :: numsta
  CHARACTER(LEN=128) :: rdstring
  CHARACTER(LEN=8) :: stnsrc
  LOGICAL :: found
  REAL :: hgtdum,rlat,rlon,wdir,speed,ddrot
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0
  OPEN(12,FILE=trim(proffile),iostat=opstat,STATUS='old')
  IF( opstat == 0 ) THEN
!
    mxprof=mxua-nprev
!
!-----------------------------------------------------------------------
!
!  Main data-reading loop
!
!-----------------------------------------------------------------------
!
    DO ista=1,mxprof
      ksta=ista+nprev
!
!-----------------------------------------------------------------------
!
!  Fill arrays with "missing data" indicator
!
!-----------------------------------------------------------------------
!
      nlevsua(ksta)=0
      DO jlev=1,nzua
        DO ivar=1,nvar
          obsua(ivar,jlev,ksta)=rmiss
          qualua(ivar,jlev,ksta)=-99
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!   Read station header
!
!-----------------------------------------------------------------------
!
      rdstring=' '
      READ(12,'(a)',iostat=rdstat) rdstring
      IF( rdstat /= 0) EXIT
      ilen=LEN_TRIM(rdstring)
      IF(ilen < 77 ) THEN    ! source not included, use default
        READ(rdstring,'(i12,i12,f11.0,f15.0,f15.0,5x,a5)',iostat=rdstat) &
          numsta,nlevsua(ksta),rlat,rlon,elevua(ksta),stnua(ksta)
        IF( rdstat /= 0) EXIT
        stnsrc=dfltsrc
      ELSE  ! header line includes source
        READ(rdstring,'(i12,i12,f11.0,f15.0,f15.0,5x,a5,1x,a8)',iostat=rdstat) &
          numsta,nlevsua(ksta),rlat,rlon,elevua(ksta),stnua(ksta),stnsrc
        IF( rdstat /= 0) EXIT
      END IF

      WRITE(6,'(2a,2x,a)') 'Reading profiler, source: ',stnua(ksta),stnsrc
!
!   Identify source number from source named
!
      found=.false.
      DO isrc=1,nsrcua
        IF(srcua(isrc) == stnsrc) THEN
          found=.true.
          WRITE(6,'(3a,i4)') ' found source: ',TRIM(stnsrc),' isrc=',isrc
          isrcua(ksta)=isrc
          EXIT
        END IF
      END DO
      IF(.NOT. found) THEN
        WRITE(6,'(a,a,a)') ' Could not find ',TRIM(stnsrc),               &
               ' among upper air sources in input file'
        CALL arpsstop('rdprof data problem',1)
      END IF

      CALL lltoxy(1,1,rlat,rlon,xua(ksta),yua(ksta))
!
      jend=MIN(nlevsua(ksta),nzua)
      DO jlev=1,jend
        wdir=-99.
        speed=-99.
        READ(12,*,iostat=rdstat) hgtua(jlev,ksta),wdir,speed

        IF(rdstat == 0 .AND. wdir >= 0. .AND. speed >= 0.) THEN
          CALL ddrotuv(1,rlon,wdir,speed,ddrot,                         &
                      obsua(1,jlev,ksta),obsua(2,jlev,ksta))
!         print *, ' direct,speed,u,v=',direct,speed,                   &
!                      obsua(1,jlev,ksta),obsua(2,jlev,ksta)
          qualua(1,jlev,ksta)=100
          qualua(2,jlev,ksta)=100
        END IF
!
      END DO
!
!-----------------------------------------------------------------------
!
!  Read any extra levels after nzua, but discard them.
!
!-----------------------------------------------------------------------
!
      IF(nlevsua(ksta) > nzua) THEN
        WRITE(6,'(//a,a/,a,i4,a)') ' RDPROF: WARNING profiler: ',       &
           stnua(ksta),' truncated to nzua=',nzua,' levels.'
        WRITE(6,'(a,i4,a)') ' Data file has ',nlevsua(ksta),' levels.'
        WRITE(6,'(a/)') ' Increase nz_ua in adas.inc'
        DO jlev=jend+1,nlevsua(ksta)
          READ(12,*,iostat=rdstat) hgtdum,wdir,speed
        END DO
        nlevsua(ksta)=nzua
      END IF

    END DO
!
!-----------------------------------------------------------------------
!
!  End-of-file destination
!
!-----------------------------------------------------------------------
!
    nprof=ista-1
    ntotal=nprev+nprof
    WRITE(6,'(a,i4,a)') ' Read ',nprof,' profiler sites'
    CLOSE(12)
    istatus=0
  ELSE
!
!-----------------------------------------------------------------------
!
!  Error opening file 
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(a)') ' Error opening profiler file: ',proffile
    ntotal=nprev
    istatus=-1
  END IF
  RETURN
END SUBROUTINE rdprof
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE RDRAOB                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdraob(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,raobfile,         &
           stnua,elevua,xua,yua,hgtua,obsua,                            &
           qualua,isrcua,nlevsua,                                       &
           rmiss,nprev,ntotal,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read ASCII file containing rawinsonde observations (RAOBs).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  Jan, 1994
!
!  MODIFICATION HISTORY:
!  9/3/96  Keith Brewster
!  Restored proper calculation of moisture variable.
!  Added full documentation.
!
!  2/16/98 Keith Brewster
!  Added jsrc, source number to the variable list
!
!  8/18/2011 Keith Brewster
!  Added srcua to the argument list and read of src from data file.
!  isrcua is then determined from the source name read-in, or dfltsrc
!  if not specified.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nvar       Number of variables in analysis (array dimension)
!    nzua       Maximum number of vertical levels (array dimension)
!    mxua       Maximum number of UA stations (array dimension)
!    jsrc       Source number of RAOB data set
!    raobfile   Name of profiler data file to read
!    rmiss      Missing data fill value
!    nprev      Number of stations read previously into UA observation
!               data arrays (array index)
!
!  OUTPUT :
!
!    stnua      station name (character*4)
!    elevua     station elevation (m MSL)
!    xua        station location x-coordinate (m)
!    yua        station location y-coordinate (m)
!    hgtua      height of data (m MSL)
!    obsua      observation data
!    qualua     observation quality indicator
!    isrcua     data source index
!    nlevsua    number of levels of data
!    ntotal     number of stations
!    istatus    status indicator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar,nzua,mxua,nsrcua
  CHARACTER (LEN=8), INTENT(IN) :: srcua(nsrcua)
  CHARACTER (LEN=8), INTENT(IN) :: dfltsrc
  CHARACTER (LEN=256), INTENT(IN) :: raobfile

  CHARACTER (LEN=5), INTENT(OUT) :: stnua(mxua)
  REAL, INTENT(OUT) :: elevua(mxua)
  REAL, INTENT(OUT) :: xua(mxua)
  REAL, INTENT(OUT) :: yua(mxua)
  REAL, INTENT(OUT) :: hgtua(nzua,mxua)
  REAL, INTENT(OUT) :: obsua(nvar,nzua,mxua)
  INTEGER, INTENT(OUT) :: qualua(nvar,nzua,mxua)
  INTEGER, INTENT(OUT) :: isrcua(mxua)
  INTEGER, INTENT(OUT) :: nlevsua(mxua)
  REAL, INTENT(IN) :: rmiss
  INTEGER, INTENT(IN) :: nprev
  INTEGER, INTENT(OUT) :: ntotal
  INTEGER, INTENT(OUT) :: istatus
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: mbtopa=100.
  CHARACTER(LEN=128) :: rdstring
  CHARACTER(LEN=8) :: stnsrc
  INTEGER :: ista,isrc,ivar,jlev,jend,ksta,nraob,mxraob,ilen
  INTEGER :: opstat,rdstat
  INTEGER :: numsta
  INTEGER :: nlevlast
  REAL :: rlat,rlon,hgtdum,press,temp,tdew,wdir,speed,ddrot,tdk
  REAL :: u1last,v1last,pr1last,pt1last,rh1last
  LOGICAL :: found
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0
  OPEN(12,FILE=trim(raobfile),iostat=opstat,STATUS='old')
  IF( opstat == 0) THEN
!
    mxraob=mxua-nprev

    nlevlast=0
    u1last=-9999.
    v1last=-9999.
    pr1last=-9999.
    pt1last=-9999.
    rh1last=-9999.
!
!-----------------------------------------------------------------------
!
!  Main data-reading loop
!
!-----------------------------------------------------------------------
!
    ntotal=nprev
    ksta=nprev+1
    DO ista=1, mxraob
!
!-----------------------------------------------------------------------
!
!  Fill arrays with "missing data" indicator
!
!-----------------------------------------------------------------------
!
      nlevsua(ksta)=0
      DO jlev=1,nzua
        DO ivar=1,nvar
          obsua(ivar,jlev,ksta)=rmiss
          qualua(ivar,jlev,ksta)=-99
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!   Read station header
!
!-----------------------------------------------------------------------
!
      rdstring=' '
      READ(12,'(a)',iostat=rdstat) rdstring
      IF(rdstat /= 0) EXIT
      ilen=LEN_TRIM(rdstring)
      IF(ilen < 77 ) THEN    ! source not included, use default
        READ(rdstring,'(i12,i12,f11.4,f15.4,f15.0,5x,a5)',iostat=rdstat) &
          numsta,nlevsua(ksta),rlat,rlon,elevua(ksta),stnua(ksta)
        IF(rdstat /= 0) EXIT
        stnsrc=dfltsrc
      ELSE  ! header line includes source
        READ(rdstring,'(i12,i12,f11.4,f15.4,f15.0,5x,a5,1x,a8)',iostat=rdstat) &
          numsta,nlevsua(ksta),rlat,rlon,elevua(ksta),stnua(ksta),stnsrc
        IF(rdstat /= 0) EXIT
      END IF

      WRITE(6,'(2a,2x,a)') 'Reading sounding, source',TRIM(stnua(ksta)),TRIM(stnsrc)
!
!   Identify source number from source named
!
      found=.false.
      DO isrc=1,nsrcua
        IF(srcua(isrc) == stnsrc) THEN
          found=.true.
          WRITE(6,'(3a,i4)') ' found source: ',TRIM(stnsrc),' isrc=',isrc
          isrcua(ksta)=isrc
          EXIT
        END IF
      END DO
      IF(.NOT. found) THEN
        WRITE(6,'(a,a,a)') ' Could not find ',TRIM(stnsrc),            &
               ' among upper air sources in input file'
        CALL arpsstop('rdprof data problem',1)
      END IF
!
      CALL lltoxy(1,1,rlat,rlon,xua(ksta),yua(ksta))
!
      jend=MIN(nlevsua(ksta),nzua)
      DO jlev=1,jend
        wdir=-999.
        speed=-999.
        press=-999.
        temp=-999.
        READ(12,*,iostat=rdstat)                                       &
            hgtua(jlev,ksta),press,temp,tdew,wdir,speed
!
!-----------------------------------------------------------------------
!
!  observed u,v
!
!-----------------------------------------------------------------------
!
        IF(rdstat == 0) THEN
          IF( wdir >= 0. .AND. speed >= 0.) THEN
            CALL ddrotuv(1,rlon,wdir,speed,ddrot,                       &
                      obsua(1,jlev,ksta),obsua(2,jlev,ksta))
            qualua(1,jlev,ksta)=100
            qualua(2,jlev,ksta)=100
          END IF
!
!-----------------------------------------------------------------------
!
!  observed press and potential temperature
!
!-----------------------------------------------------------------------
!
          IF(press >= 0.) THEN
            obsua(3,jlev,ksta)=press*mbtopa
            qualua(3,jlev,ksta)=100
            IF(temp > -99.) THEN
              obsua(4,jlev,ksta)=                                       &
                  (temp+273.15)*((p0/obsua(3,jlev,ksta))**rddcp)
              qualua(4,jlev,ksta)=100
            END IF
!
!-----------------------------------------------------------------------
!
!  observed specific humidity
!
!-----------------------------------------------------------------------
!
            IF(temp > -99. .AND. tdew > -99.) THEN
              tdk=tdew + 273.15
              obsua(5,jlev,ksta)=                                         &
                   MAX(1.0E-08,f_qvsat(obsua(3,jlev,ksta),tdk))
              qualua(5,jlev,ksta)=100
            END IF
          END IF  ! non-missing pressure
        END IF ! rdstat
      END DO  ! level loop
!
!-----------------------------------------------------------------------
!
!  Read any extra levels after nzua, but discard them.
!
!-----------------------------------------------------------------------
!
      IF(nlevsua(ksta) > nzua) THEN
        WRITE(6,'(//a,a/,a,i4,a)') ' RDRAOB: WARNING rawinsonde: ',       &
             stnua(ksta),' truncated to nzua=',nzua,' levels.'
        WRITE(6,'(a,i4,a)') ' Data file has ',nlevsua(ksta),' levels.'
        WRITE(6,'(a/)') ' Increase nz_ua in adas.inc'
!
        DO jlev=jend+1,nlevsua(ksta)
          READ(12,*,iostat=rdstat)                                        &
            hgtdum,press,temp,tdew,wdir,speed
        END DO
        nlevsua(ksta)=nzua
      END IF
!
!-----------------------------------------------------------------------
!
!  Check for duplicate sounding to to error some .snd files
!  If dupe is found, do not increment ksta or ntotal counter.
!
!-----------------------------------------------------------------------
!
      IF(nlevsua(ksta) == nlevlast .AND. obsua(1,1,ksta) == u1last .AND.&
          obsua(2,1,ksta) == v1last .AND.                               &
          obsua(3,1,ksta) == pr1last .AND.                              &
          obsua(4,1,ksta) == pt1last .AND.                              &
          obsua(5,1,ksta) == rh1last ) THEN

        WRITE(6,'(a,a)')                                                &
          ' Discarding erroneous duplicate data stored for',            &
          stnua(ksta)
      ELSE

        WRITE(6,'(a,i5,2a)') ' Read', nlevsua(ksta),                    &
          ' levels of data from radiosonde station ', TRIM(stnua(ksta))
        nlevlast=nlevsua(ksta)
        u1last=obsua(1,1,ksta)
        v1last=obsua(2,1,ksta)
        pr1last=obsua(3,1,ksta)
        pt1last=obsua(4,1,ksta)
        rh1last=obsua(5,1,ksta)

        ksta=ksta+1
        ntotal=ntotal+1

      END IF

    END DO  ! station loop
!
!-----------------------------------------------------------------------
!
!  End-of-file destination
!
!-----------------------------------------------------------------------
!
    nraob=ntotal-nprev
    CLOSE(12)
    istatus=0
!
!-----------------------------------------------------------------------
!
!  Error opening file destination
!
!-----------------------------------------------------------------------
!
  ELSE 
    WRITE(6,'(a)') ' Error opening raob file: ',raobfile
    ntotal=nprev
    istatus=-1
  END IF
  RETURN
END SUBROUTINE rdraob
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE RDGPS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rdgps(nx,ny,nz,nvar,nzua,mxua,nsrcua,gpsfile,                &
           anx,xs,ys,zp,su,sv,st,spres,shght,sqv,sqvsat,                &
           srcua,stnua,elevua,xua,yua,hgtua,obsua,                      &
           qualua,isrcua,nlevsua,                                       &
           rmiss,nprev,ntotal,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read ASCII file containing GPS Integrated Precipitable Water
!  data and generate a psuedo-sounding of qv from the background 
!  qv field at the data location and the IPW information.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  June, 2006
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nvar       Number of variables in analysis (array dimension)
!    nzua       Maximum number of vertical levels (array dimension)
!    mxua       Maximum number of UA stations (array dimension)
!    jsrc       Source number of RAOB data set
!    raobfile   Name of profiler data file to read
!    rmiss      Missing data fill value
!    nprev      Number of stations read previously into UA observation
!               data arrays (array index)
!
!  OUTPUT :
!
!    stnua      station name (character*4)
!    elevua     station elevation (m MSL)
!    xua        station location x-coordinate (m)
!    yua        station location y-coordinate (m)
!    hgtua      height of data (m MSL)
!    obsua      observation data
!    qualua     observation quality indicator
!    isrcua     data source index
!    nlevsua    number of levels of data
!    ntotal     number of stations
!    istatus    status indicator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nvar,nzua,mxua,nsrcua
!
  CHARACTER (LEN=132) :: gpsfile
  CHARACTER (LEN=8) :: srcua(nsrcua)
  CHARACTER (LEN=5) :: stnua(mxua)
  REAL :: anx(nvar,nx,ny,nz)
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: su(nz)
  REAL :: sv(nz)
  REAL :: st(nz)
  REAL :: spres(nz)
  REAL :: shght(nz)
  REAL :: sqv(nz)
  REAL :: sqvsat(nz)
  REAL :: elevua(mxua)
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: obsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nlevsua(mxua)
  REAL :: rmiss
  INTEGER :: nprev,ntotal
  INTEGER :: istatus
!
  REAL, PARAMETER :: mbtopa=100.
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
! GPS iteration parameters
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: maxiter = 20
  REAL, PARAMETER :: epsipw = 2.0E-03
  REAL, PARAMETER :: delvmax = 200.
  REAL, PARAMETER :: stdlapse = -(6.5/1000.)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=5) :: stngps
  CHARACTER(LEN=8) :: srcgps
  CHARACTER(LEN=12) :: datetime
  INTEGER :: ista,isrc,ivar,ksta,k,kk,mxgps,ngps,nlevs
  INTEGER :: ipt,jpt,ireturn,iostatus,niter
  REAL :: rlat,rlon,elev,sfcpr,sfct,sfcrh,gpsipw
  REAL :: sfcqv,sfcqvsat,sfcqvflg
  REAL :: xgps,ygps,selev,temp
  LOGICAL :: found
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  OPEN(12,FILE=trim(gpsfile),STATUS='old',iostat=iostatus)
  IF(iostatus /= 0) THEN
    WRITE(6,'(/a,a/)') ' Error opening GPS file: ',gpsfile
    ntotal=nprev
    RETURN
  END IF
  READ(12,'(a12)',iostat=iostatus) datetime
  IF(iostatus /= 0) THEN
    WRITE(6,'(/a,a/)') ' Error reading header in GPS file: ',gpsfile
    ntotal=nprev
    CLOSE(12)
    RETURN
  END IF
  WRITE(6,'(//a,a)') ' Reading GPS data for date-time: ',datetime
!
  mxgps=mxua-nprev
!
!-----------------------------------------------------------------------
!
!  Main data-reading loop
!
!-----------------------------------------------------------------------
!
  ksta=nprev
  ntotal=nprev
  DO ista=1, mxgps
    READ(12,'(a5,1x,a8,9x,f10.5,f11.5,f8.2,f8.2,f6.2,f8.2,f7.1)',   &
         iostat=iostatus)                                              &
        stngps,srcgps,rlat,rlon,elev,sfcpr,gpsipw,sfct,sfcrh
    print *, ' '
    WRITE(6,'(/a,a5,1x,a8,1x,f10.5,f11.5,f8.2)')                       &
      ' stngps,srcgps,rlat,rlon,elev: ',stngps,srcgps,rlat,rlon,elev
    WRITE(6,'(a,4f8.2)')                                               &
      ' sfcpr,gpsipw,sfct,sfcrh: ',sfcpr,gpsipw,sfct,sfcrh
    IF(iostatus /= 0) EXIT
    IF(gpsipw < 0. .OR. sfcpr < 0. ) CYCLE

    sfcpr=mbtopa*sfcpr
!
    CALL lltoxy(1,1,rlat,rlon,xgps,ygps)
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid.
!
!-----------------------------------------------------------------------
!
    CALL findlc(nx,ny,xs,ys,xgps,ygps,ipt,jpt,ireturn)
    WRITE(6,'(a,2f12.1,2i8)') ' xgps,ygps,ipt,jpt:',xgps,ygps,ipt,jpt
!
!-----------------------------------------------------------------------
!
!  Interpolate if extrapoltion is not indicated.
!
!-----------------------------------------------------------------------
!
    IF(ireturn >= 0) THEN
      found=.false.
      DO isrc=1,nsrcua
        IF(srcua(isrc) == srcgps) THEN
          found=.true.
          EXIT
        END IF
      END DO
      IF(.NOT. found) THEN
        WRITE(6,'(a,a,a)') ' Could not find GPS source ',srcgps,       &
             ' among upper air sources in input file'
        STOP
      END IF

      write(6,'(a,a)')  stngps,                                        &
                  ' inside ADAS domain, creating pseudo sounding'
 
      CALL colint(nx,ny,nz,nvar,                                       &
             xs,ys,zp,xgps,ygps,ipt,jpt,anx,                           &
             su,sv,st,spres,shght,sqv,selev,                           &
             nlevs)

      WRITE(6,'(a,f6.1,a,f6.1)') &
            ' GPS Elev: ',elev,' Model elev: ',selev

      IF(abs(elev-selev) > delvmax) CYCLE
!
!     Convert theta to temperature and compute qv_sat
!
      DO k=1,nlevs
        temp=st(k)*((spres(k)/p0)**rddcp)
        sqvsat(k)=f_qvsat(spres(k),temp)
        st(k)=temp
      END DO

      IF(sfct > -50.0 ) THEN
        sfct=sfct+273.15
      ELSE
        sfct=st(1)+stdlapse*(elev-shght(1))
      END IF

      sfcqvsat=f_qvsat(sfcpr,sfct)

      IF(sfcrh > 0.) THEN
        sfcqv=(0.01*sfcrh)*sfcqvsat
        sfcqvflg=1
      ELSE
        sfcqv=sqv(1)
        sfcqvflg=0
      END IF

      CALL ipw2snd(nlevs,shght,spres,st,sqv,sqvsat,                    &
                   elev,sfcpr,sfct,sfcqv,sfcqvsat,sfcqvflg,gpsipw,     &
                   maxiter,epsipw,niter,istatus)

      IF(istatus == 0) THEN
        ksta=ksta+1
        ntotal=ntotal+1
        DO k=1,nlevs
          DO ivar=1,nvar
            obsua(ivar,k,ksta)=rmiss
            qualua(ivar,k,ksta)=-13
          END DO
        END DO
        xua(ksta)=xgps
        yua(ksta)=ygps
        elevua(ksta)=elev
        stnua(ksta)=stngps
        isrcua(ksta)=isrc
        hgtua(1,ksta)=elev
        obsua(5,1,ksta)=sfcqv
        qualua(5,1,ksta)=100
        kk=1
        DO k=1,nlevs
          IF(spres(k) < sfcpr .AND. shght(k) > elev) THEN
            kk=kk+1
            hgtua(kk,ksta)=shght(k)
            obsua(5,kk,ksta)=sqv(k)
            qualua(5,kk,ksta)=100
          END IF
        END DO
        nlevsua(ksta)=kk

      END IF  ! success in IPW adjustment

    END IF  ! location in grid

  END DO
!
!-----------------------------------------------------------------------
!
!  End-of-file destination
!
!-----------------------------------------------------------------------
!
  ngps=ntotal-nprev
  WRITE(6,'(/a,i6,a/)')' Read and processed ',ngps,' GPS stations'
  CLOSE(12)
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error opening file destination
!
!-----------------------------------------------------------------------
!
  RETURN
END SUBROUTINE rdgps
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE IPW2SND                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ipw2snd(nzsnd,htsnd,psnd,tsnd,qvsnd,qvsat,                  &
                   elev,sfcpr,sfct,sfcqv,sfcqvsat,sfcqvflg,obsipw,     &
                   maxiter,epsipw,niter,istatus)
!
! Given a background sounding, a measurement of integrated 
! precipitable water (IPW), surface pressure, temperature
! and specific humidity generate a pseudo-sounding of 
! specific humidity.
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN) :: htsnd(nzsnd)
  REAL, INTENT(IN) :: psnd(nzsnd)
  REAL, INTENT(IN) :: tsnd(nzsnd)
  REAL, INTENT(INOUT) :: qvsnd(nzsnd)
  REAL, INTENT(INOUT) :: qvsat(nzsnd)
  REAL, INTENT(IN) :: elev
  REAL, INTENT(IN) :: sfcpr
  REAL, INTENT(IN) :: sfct
  REAL, INTENT(INOUT) :: sfcqv
  REAL, INTENT(IN) :: sfcqvsat
  INTEGER, INTENT(IN) :: sfcqvflg
  REAL, INTENT(IN) :: obsipw
  INTEGER, INTENT(IN) :: maxiter
  REAL, INTENT(IN) :: epsipw
  INTEGER, INTENT(OUT) :: niter 
  INTEGER, INTENT(OUT) :: istatus
!
!  Misc internal variables
!
  REAL, PARAMETER :: g=9.8
  REAL, PARAMETER :: rhmin=0.05
  REAL, PARAMETER :: rhmax=1.0
  REAL :: pwcst,sfcipw,sndipw,qvbar,dipw,pwratio,qvnew
  REAL :: f_qvsat
  INTEGER :: k,kstart,iter
!
! Initializations
!
  istatus=0
  pwcst=0.1/g
  DO k=1,nzsnd-1
    IF(psnd(k) > 0.) THEN
      qvsat(k)=f_qvsat(psnd(k),tsnd(k))
    ELSE
      qvsat(k)=1.0E-06
    END IF
  END DO

  DO k=1,nzsnd-1
    IF(psnd(k) < sfcpr .AND. htsnd(k) > elev) THEN
      kstart=k
      EXIT
    END IF
  END DO
!
! Iterate until integrated precipitable water difference
! falls below the threshold, epsilon ipw
!
  IF(sfcqvflg == 1) THEN     ! observed surface qv
    DO iter=1,maxiter
      qvbar=0.5*(sfcqv+qvsnd(kstart))
      sndipw=pwcst*qvbar*(sfcpr-psnd(kstart))
      sfcipw=0.5*sndipw
      DO k=kstart,nzsnd-1
        IF(psnd(k) > 0. .AND. psnd(k+1) > 0.) THEN
          qvbar=0.5*(qvsnd(k)+qvsnd(k+1))
          sndipw=sndipw+(pwcst*qvbar*(psnd(k)-psnd(k+1)))
        END IF
      END DO
      dipw=abs(sndipw-obsipw)
      pwratio=(obsipw-sfcipw)/(sndipw-sfcipw)
      WRITE(6,'(a,i4,4(a,f9.2))') &
        ' Iter:',iter,' IPW:',obsipw,' Snd IPW:',sndipw, &
        ' Diff:',dipw,' Ratio:',pwratio
      IF( dipw < epsipw) EXIT
      DO k=1,nzsnd
        IF(psnd(k) > 0.) THEN
          qvnew=pwratio*qvsnd(k)
          qvsnd(k)=min(max(qvnew,(rhmin*qvsat(k))),(rhmax*qvsat(k)))
        END IF
      END DO
    END DO
  ELSE ! sfcqv is estimated
    DO iter=1,maxiter
      qvbar=0.5*(sfcqv+qvsnd(kstart))
      sndipw=pwcst*qvbar*(sfcpr-psnd(kstart))
      DO k=kstart,nzsnd-1
        IF(psnd(k) > 0. .AND. psnd(k+1) > 0.) THEN
          qvbar=0.5*(qvsnd(k)+qvsnd(k+1))
          sndipw=sndipw+(pwcst*qvbar*(psnd(k)-psnd(k+1)))
        END IF
      END DO
      dipw=abs(sndipw-obsipw)
      pwratio=obsipw/sndipw
      WRITE(6,'(a,i4,4(a,f9.2))') &
        ' Iter:',iter,' IPW:',obsipw,' Snd IPW:',sndipw, &
        ' Diff:',dipw,' Ratio:',pwratio
      IF( dipw < epsipw) EXIT
      qvnew=pwratio*sfcqv
      sfcqv=min(max(qvnew,(rhmin*sfcqvsat)),(rhmax*sfcqvsat))
      DO k=1,nzsnd
        IF(psnd(k) > 0.) THEN
          qvnew=pwratio*qvsnd(k)
          qvsnd(k)=min(max(qvnew,(rhmin*qvsat(k))),(rhmax*qvsat(k)))
        END IF
      END DO
    END DO
  END IF

  IF(iter > maxiter) THEN
    niter=maxiter
    istatus=-1
  ELSE
    niter=iter
    istatus=0
  END IF
  RETURN
  
END SUBROUTINE IPW2SND
