SUBROUTINE read_surface_obs(infile,blackfile,                           &
           maxsta,istart,atime,n_meso_g,                                &
           n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,  &
           n_obs_pos_g,n_obs_b,n_obs_pos_b,stn,obstype,lat,lon,elev,wx, &
           t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad, &
           idp3,store_emv,store_amt,store_hgt,vis,obstime,istatus)
!
!*******************************************************************************
!
!    Routine to read mesonet and SAO data for LAPS that has been written
!    into the .LSO file by the 'get_surface_obs' routine.
!
!    Changes:
!            P. Stamus       12-30-92  Original version.
!                            01-07-93  Add/change obs counters.
!                            01-08-93  Add read_header entry.
!            K. Brewster     12-20-99  Moved the initialization of
!                                      cloud amts to missing to
!                                      inside of station read loop
!
!    Input/Output:
!
!     Variable        Var type   I/O   Description
!    ----------      ---------- ----- -------------
!     infile            A*80      I    Directory where LSO file is.
!     maxsta             I        I    Max Number of stations allowed.
!     atime             A*24      O    Data time in dd-mmm-yyyy hh:mm
!     n_meso_g           I        O    Number of FSL mesonet stations
!     n_meso_pos         I        O    Total number mesonet stations psbl
!     n_sao_g            I        O    Number of SAOs in the laps grid
!     n_sao_pos_g        I        O    Total num. of SAOs psbl in laps grid
!     n_sao_b            I        O    Number of SAOs in the box
!     n_sao_pos_b        I        O    Total num of SAOs psbl in the box
!     n_obs_g            I        O    Number of obs in the laps grid
!     n_obs_pos_g        I        O    Total num of obs psbl in the laps grid
!     n_obs_b            I        O    Number of obs in the box
!     n_obs_pos_b        I        O    Total num of obs possible in the box
!     stn               A*3 A     O    Station names (array)
!     lat                RA       O    Station latitude (deg)
!     lon                RA       O    Station longitude (deg)
!     elev               RA       O    Station elevation (m)
!     obstype           A*8 A     O    Observation type (SA, SP, ASOS, etc)
!     obstime            IA       O    Time of observation (hhmm)
!     wx                A*8 A     O    Observed weather
!     t                  RA       O    Temperature (F)
!     td                 RA       O    Dewpoint (F)
!     dd                 RA       O    Wind direction (deg)
!     ff                 RA       O    Wind speed (kt)
!     ddg                RA       O    Gust wind direction (deg)
!     ffg                RA       O    Gust wind speed (kt)
!     pstn               RA       O    Station pressure (mb)
!     pmsl               RA       O    MSL pressure (mb)
!     alt                RA       O    Altimeter setting (mb)
!     kloud              IA       O    Number of cloud layers...max of 5.
!     ceil               RA       O    Ceiling height (m)
!     lowcld             RA       O    Height lowest cloud (m)
!     cover              RA       O    Cloud cover (tenths)
!     vis                RA       O    Visibility (miles)
!     rad                RA       O    Solar radiation.
!     idp3               IA       O    3-h coded pressure change (e.g.,608)
!     store_emv         A*1 A     O    Cloud descriptors: ea. layer, ea. stn
!     store_amt         A*4 A     O    Cloud layer coverage.
!     store_hgt          RA       O    Height of each cloud layer.
!     istatus            I        O    Status flag: 1 = normal
!                                                  -1 = file not found
!                                                  -2 = Arrays too small
!
!    User Notes:
!
!    1.  Arrays should be dimensioned 'maxsta' in the calling program,
!        with maxsta *at least* 120 (for CO domain).
!
!    2.  Pressures are stored as reported, except that altimeters are
!        converted to millibars.
!
!    3.  The 'kloud' variable tells whether there are clouds and how
!        many layers if there are:
!            a) kloud = 0    means   No cloud DATA (but NOT "no clouds").
!            b) kloud = 1    means   CLR or 1 cloud layer.  A height is
!                                    given for CLR which is the maximum valid
!                                    height of the observation (automatic
!                                    stations have limited valid heights).
!            c) kloud = 2-5  means   Two to five cloud layers.
!
!    4.  Thin obscured (-X) is a cloud layer and is given a 'badflag'
!        height, since it is not supposed to have a height (you're supposed
!        to be able to see other clouds and/or sky).
!
!*******************************************************************************
!
  IMPLICIT NONE
  INTEGER :: maxsta
  INTEGER :: istart
  INTEGER :: misval
  PARAMETER       (misval=999)
  REAL :: badflag
  PARAMETER       (badflag = -99.9)
!
  REAL :: lat(maxsta),lon(maxsta),elev(maxsta)
  REAL :: t(maxsta),td(maxsta)
  REAL :: dd(maxsta),ff(maxsta)
  REAL :: ddg(maxsta),ffg(maxsta)
  REAL :: pstn(maxsta),pmsl(maxsta),alt(maxsta)
  REAL :: store_hgt(maxsta,5)
  REAL :: ceil(maxsta),lowcld(maxsta),cover(maxsta)
  REAL :: vis(maxsta),rad(maxsta)
!
  INTEGER :: obstime(maxsta),kloud(maxsta),idp3(maxsta)
!
  CHARACTER (LEN=*)  :: infile
  CHARACTER (LEN=*)  :: blackfile
  CHARACTER (LEN=24) :: atime
  CHARACTER (LEN=5)  :: stn(maxsta)
  CHARACTER (LEN=8)  :: obstype(maxsta)
  CHARACTER (LEN=8)  :: wx(maxsta)
  CHARACTER (LEN=1)  :: store_emv(maxsta,5)
  CHARACTER (LEN=4)  :: store_amt(maxsta,5)
!
  INTEGER :: n_meso_g,n_meso_pos,n_sao_g,n_sao_pos_g
  INTEGER :: n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g
  INTEGER :: n_obs_b,n_obs_pos_b
  INTEGER :: istatus,iostatus
!
  INTEGER :: i,k,ios,ii
!
!.....  Start here.  Set the status to nothing, zero out the cloud storage.
!
  istatus = 0
!
!.....  Open the file.  Check for a 'file not found' or other problem.
!
  OPEN(1,IOSTAT=ios,FILE=trim(infile),STATUS='old',                     &
       ACCESS='sequential',FORM='formatted')
  IF(ios /= 0) THEN     ! error during read
    istatus = -1
    WRITE(6,650) infile
    650       FORMAT(' +++ ERROR opening: ',a80,' +++')
    WRITE(6,651) ios
    651       FORMAT('     IOS code = ',i5)
    RETURN
  END IF
!
!.....  File open...first read the header.
!
  READ(1,900) atime,              & ! data time
           n_meso_g,          & ! # of mesonet stations
           n_meso_pos,        & ! total # mesonet stations possible
           n_sao_g,           & ! # of saos in the laps grid
           n_sao_pos_g,       & ! total # of saos possible in laps grid
           n_sao_b,           & ! # of saos in the box
           n_sao_pos_b,       & ! total # of saos possible in the box
           n_obs_g,           & ! # of obs in the laps grid
           n_obs_pos_g,       & ! total # of obs psbl in the laps grid
           n_obs_b,           & ! # of obs in the box
           n_obs_pos_b        ! total # of obs possible in the box
  900     FORMAT(1X,a24,10(1X,i4))
!wdt tmp
!  1900     FORMAT("XXX header:",1X,a24,10(1X,i4))
!  WRITE(*,1900) atime,              & ! data time
!           n_meso_g,          & ! # of mesonet stations
!           n_meso_pos,        & ! total # mesonet stations possible
!           n_sao_g,           & ! # of saos in the laps grid
!           n_sao_pos_g,       & ! total # of saos possible in laps grid
!           n_sao_b,           & ! # of saos in the box
!           n_sao_pos_b,       & ! total # of saos possible in the box
!           n_obs_g,           & ! # of obs in the laps grid
!           n_obs_pos_g,       & ! total # of obs psbl in the laps grid
!           n_obs_b,           & ! # of obs in the box
!           n_obs_pos_b        ! total # of obs possible in the box
!
!.....  Error trapping for too many stations for array size.
!
  IF((n_obs_b+istart-1) > maxsta) THEN
    PRINT 990, maxsta,(n_obs_b+istart-1),atime
    990       FORMAT(' +++ ERROR in READ_SURFACE_OBS: maxsta =',i9,/, &
           ' but there will be',i9,' stations after reading file',a24/)
    PRINT *,'    Increase the value of "maxsta" and try again.'
    istatus = -2
    RETURN
  END IF
!
!.....  Now read the station data.
!
  DO i=1,n_obs_b
    k=i+istart-1
    READ(1,901,iostat=iostatus) stn(k),lat(k),lon(k),elev(k),obstype(k),                &
        obstime(k),wx(k)
    901       FORMAT(1X,a5,f6.2,1X,f7.2,1X,f5.0,1X,a8,1X,i4,1X,a8)
    IF(iostatus /= 0) EXIT
!wdt tmp
!    WRITE(*,1901) stn(k),lat(k),lon(k),elev(k),obstype(k),                &
!        obstime(k),wx(k)
!    1901       FORMAT("XXX stn dta:",1X,a5,f6.2,1X,f7.2,1X,f5.0,1X,a8,1X,i4,1X,a8)
!
    READ(1,903) t(k),td(k),dd(k),ff(k),ddg(k),ffg(k),pstn(k),           &
                pmsl(k),alt(k)
    903       FORMAT(4X,2(f6.1,1X),4(f5.0,1X),3(f6.1,1X))
!wdt tmp
!    WRITE(*,1903) t(k),td(k),dd(k),ff(k),ddg(k),ffg(k),pstn(k),           &
!                pmsl(k),alt(k)
!    1903       FORMAT("XXX         ",4X,2(f6.1,1X),4(f5.0,1X),3(f6.1,1X))
!
    READ(1,905) kloud(k),ceil(k),lowcld(k),cover(k),vis(k),rad(k),      &
        idp3(k)
    905       FORMAT(4X,i2,2(1X,f7.1),1X,f5.1,1X,f7.3,1X,f6.1,1X,i4)
!wdt tmp
!    WRITE(*,1905) kloud(k),ceil(k),lowcld(k),cover(k),vis(k),rad(k),      &
!        idp3(k)
!    1905       FORMAT("XXX cloud   ",4X,i2,2(1X,f7.1),1X,f5.1,1X,f7.3,1X,f6.1,1X,i4)
!
!.....  Read the cloud data if we have any.
!
    DO ii=1,5
      store_emv(k,ii) = ' '
      store_amt(k,ii) = '    '
      store_hgt(k,ii) = badflag
    END DO ! ii
    IF(kloud(k) > 0) THEN
!wdt update
      DO ii=1,kloud(k)
        IF (ii <= 5) THEN
          READ(1,'(5X,a1,1X,a4,1X,f7.1)')                               &
                        store_emv(k,ii),store_amt(k,ii),store_hgt(k,ii)
!wdt tmp
!        WRITE(*,1907) store_emv(k,ii),store_amt(k,ii),store_hgt(k,ii)
!        1907           FORMAT("XXX storm   ",5X,a1,1X,a4,1X,f7.1)
        ELSE
          READ (1,*)
          WRITE (*,*) "SKIPPING extra cloud layers for station ",stn(k)
        END IF
      END DO !ii
    END IF
!
  END DO !k
!
  CLOSE(1)
!
  CALL blklistsfc(blackfile,maxsta,istart,n_obs_b,misval,badflag,       &
      stn,obstype,wx,                                                   &
      t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad,     &
      idp3,store_emv,store_amt,store_hgt,vis,istatus)
!
!..... End of data gathering. Let's go home...
!
  istatus = 1             ! everything's ok...
  PRINT *, ' Normal completion of READ_SURFACE_OBS'
!
  RETURN
!
!
  ENTRY read_surface_header(infile,atime,n_meso_g,n_meso_pos,           &
      n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g,      &
      n_obs_b,n_obs_pos_b,istatus)
!
!.....  Entry to read and pass back just the header info from the lso file.
!
!.....  Open the file.  Check for a 'file not found' or other problem.
!
  istatus = 0
  OPEN(1,IOSTAT=ios,FILE=trim(infile),STATUS='old',                     &
       ACCESS='sequential',FORM='formatted')
  IF(ios /= 0) THEN     ! error during read
    istatus = -1
    WRITE(6,650) infile
    WRITE(6,651) ios
    RETURN
  END IF
!
!.....  File open...first read the header.
!
  READ(1,900) atime,              & ! data time
           n_meso_g,          & ! # of mesonet stations
           n_meso_pos,        & ! total # mesonet stations possible
           n_sao_g,           & ! # of saos in the laps grid
           n_sao_pos_g,       & ! total # of saos possible in laps grid
           n_sao_b,           & ! # of saos in the box
           n_sao_pos_b,       & ! total # of saos possible in the box
           n_obs_g,           & ! # of obs in the laps grid
           n_obs_pos_g,       & ! total # of obs psbl in the laps grid
           n_obs_b,           & ! # of obs in the box
           n_obs_pos_b        ! total # of obs possible in the box
!
!.....  Rewind and close the file so we can call this again in the same program.
!
  CLOSE(1)
  istatus = 1
!
  RETURN
END SUBROUTINE read_surface_obs
!

SUBROUTINE blklistsfc(blkfile,maxsta,istart,nobs,misval,rmisval,        &
           stn,obstype,wx,                                              &
           t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad, &
           idp3,store_emv,store_amt,store_hgt,vis,istatus)
!
!  Reads list of stations and variables from a file
!  to force variables known to be chronically wrong to "missing".
!
!  Keith Brewster, CAPS
!  May, 1995
!
  IMPLICIT NONE
!
!  Arguments
!
  CHARACTER (LEN=*) :: blkfile
  INTEGER           :: maxsta,istart,nobs,misval
  REAL              :: rmisval
  CHARACTER (LEN=5) :: stn(maxsta)
  REAL :: t(maxsta),td(maxsta)
  REAL :: dd(maxsta),ff(maxsta),ddg(maxsta),ffg(maxsta)
  REAL :: pstn(maxsta),pmsl(maxsta),alt(maxsta)
  REAL :: store_hgt(maxsta,5)
  REAL :: ceil(maxsta),lowcld(maxsta),cover(maxsta)
  REAL :: vis(maxsta),rad(maxsta)
!
  INTEGER :: kloud(maxsta),idp3(maxsta)
!
  CHARACTER (LEN=8) :: wx(maxsta)
  CHARACTER (LEN=8) :: obstype(maxsta)
  CHARACTER (LEN=1) :: store_emv(maxsta,5)
  CHARACTER (LEN=4) :: store_amt(maxsta,5)
  INTEGER :: istatus,iostatus
!
!  Some internal parameters
!
  INTEGER :: maxblk
  PARAMETER (maxblk=10)
  CHARACTER (LEN=40) :: BLANK
  PARAMETER(BLANK='                                        ')
!
!  Misc internal variables
!
  CHARACTER (LEN=5) :: blkstn
  INTEGER :: blkvar(maxblk)
  INTEGER :: i,iunit,ista,ivar,nblack
  LOGICAL :: found
!
  CALL getunit(iunit)
  OPEN(iunit,FILE=trim(blkfile),STATUS='old')
  DO 
    READ(iunit,'(a5,i5)',IOSTAT=iostatus) blkstn,nblack
    IF(iostatus < 0 ) EXIT
    IF(iostatus > 0 ) THEN
      WRITE(6,'(a,a)') ' Error reading blacklist file ',          &
       trim(blkfile)
      EXIT
    END IF
    DO i=1,maxblk
      blkvar(i)=-99
    END DO
    READ(iunit,*) (blkvar(i),i=1,nblack)
    found=.false.
    DO i=1,nobs
      ista=i+istart-1
      IF(stn(ista) == blkstn) THEN
        found=.true.
        DO ivar=1,nblack
          SELECT CASE (blkvar(ivar))
          CASE(1)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Temperature',    &
              blkvar(ivar),' at ',blkstn
            t(ista)=rmisval
          CASE(2)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Dew Point',      &
              blkvar(ivar),' at ',blkstn
            td(ista)=rmisval
          CASE(3)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Wind Variables', &
              blkvar(ivar),' at ',blkstn
            dd(ista)=rmisval
            ff(ista)=rmisval
            ddg(ista)=rmisval
            ffg(ista)=rmisval
          CASE(4)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Press Variables',&
              blkvar(ivar),' at ',blkstn
            pstn(ista)=rmisval
            pmsl(ista)=rmisval
            alt(ista)=rmisval
          CASE(5)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Cloud variables',&
              blkvar(ivar),' at ',blkstn
            kloud(ista)=0
            ceil(ista)=rmisval
            lowcld(ista)=rmisval
            cover(ista)=rmisval
          CASE(6)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Visibility',     &
              blkvar(ivar),' at ',blkstn
            vis(ista)=rmisval
          CASE(7)
            WRITE(6,'(a,i4,a,a)') ' Blacklisting Solar Radiation',&
              blkvar(ivar),' at ',blkstn
            rad(ista)=rmisval
          CASE DEFAULT
            WRITE(6,'(a,i4,a,a)') ' Invalid variable number ',    &
              blkvar(ivar),' at ',blkstn
          END SELECT
        END DO
        CYCLE
      END IF
    END DO
    IF(.NOT. found) WRITE(6,'(a,a,a)')                       &
      ' Blacklisted station ',blkstn,' not found in dataset'
  END DO
  CLOSE (iunit)
  CALL retunit(iunit)
  RETURN
END SUBROUTINE blklistsfc
