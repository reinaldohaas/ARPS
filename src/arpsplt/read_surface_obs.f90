!

SUBROUTINE read_surface_obs(infile,maxsta,atime,n_meso_g,               &
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
!
!    Input/Output:
!
!     Variable        Var type   I/O   Description
!    ----------      ---------- ----- -------------
!     infile            A*70      I    Directory where LSO file is.
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
  REAL :: badflag
  PARAMETER (badflag = -99.9)

  REAL :: lat(maxsta),lon(maxsta),elev(maxsta),t(maxsta),td(maxsta),    &
       dd(maxsta),ff(maxsta),ddg(maxsta)
  REAL :: ffg(maxsta),pstn(maxsta),pmsl(maxsta),                        &
       alt(maxsta),store_hgt(maxsta,5)
  REAL :: ceil(maxsta),lowcld(maxsta),cover(maxsta),                    &
       vis(maxsta),rad(maxsta)

  INTEGER :: obstime(maxsta),kloud(maxsta),idp3(maxsta)

  CHARACTER (LEN=*) :: infile
  CHARACTER (LEN=24) :: atime
  CHARACTER (LEN=5) :: stn(maxsta)
  CHARACTER (LEN=8) :: obstype(maxsta)
  CHARACTER (LEN=8) :: wx(maxsta)
  CHARACTER (LEN=1) :: store_emv(maxsta,5)
  CHARACTER (LEN=4) :: store_amt(maxsta,5)
!
!.....  Start here.  Set the status to nothing, zero out the cloud storage.
!
  istatus = 0
  DO j=1,5
    DO i=1,maxsta
      store_emv(i,j) = ' '
      store_amt(i,j) = '    '
      store_hgt(i,j) = badflag
    END DO !i
  END DO !j
!
!.....  Open the file.  Check for a 'file not found' or other problem.
!
  OPEN(1,IOSTAT=ios,FILE=infile,STATUS='old',ACCESS='sequential',       &
       FORM='formatted')
  IF(ios /= 0) THEN     ! error during read
    istatus = -1
    n_obs_b = 0
    WRITE(6,650) infile
    650       FORMAT(' +++ ERROR opening: ',a70,' +++')
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
!
!.....  Error trapping for too many stations for array size.
!
  IF(n_obs_b > maxsta) THEN
    PRINT 990, maxsta,n_obs_b,atime
    990       FORMAT(' +++ ERROR in READ_SURFACE_OBS: maxstns = ',i8,/, &
           ' but there are ',i8,' stations in the ',a24,' obs file.',/)
    PRINT *,'    Increase the value of "maxstns" and try again.'
    n_obs_b = 0
    istatus = -2
    RETURN
  END IF
!
!.....  Now read the station data.
!
  DO k=1,n_obs_b
    READ(1,901) stn(k),lat(k),lon(k),elev(k),obstype(k),                &
        obstime(k),wx(k)
    901       FORMAT(1X,a5,f6.2,1x,f7.2,1X,f5.0,1X,a8,1X,i4,1X,a8)
!
    READ(1,903) t(k),td(k),dd(k),ff(k),ddg(k),ffg(k),pstn(k),           &
                pmsl(k),alt(k)
    903       FORMAT(4X,2(f6.1,1X),4(f5.0,1X),3(f6.1,1X))
!
    READ(1,905) kloud(k),ceil(k),lowcld(k),cover(k),vis(k),rad(k),      &
        idp3(k)
    905       FORMAT(4X,i2,2(1X,f7.1),1X,f5.1,1X,f7.3,1X,f6.1,1X,i4)
!
!.....  Read the cloud data if we have any.
!
    IF(kloud(k) > 0) THEN
      DO ii=1,kloud(k)
        READ(1,907) store_emv(k,ii),store_amt(k,ii),store_hgt(k,ii      &
            )
        907           FORMAT(5X,a1,1X,a4,1X,f7.1)
      END DO !ii
    END IF
!
  END DO !k
!
  REWIND(1)
  CLOSE(1)
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
  OPEN(1,IOSTAT=ios,FILE=infile,STATUS='old',ACCESS='sequential',       &
       FORM='formatted')
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
  REWIND(1)
  CLOSE(1)
  istatus = 1
!
  RETURN
END SUBROUTINE read_surface_obs
