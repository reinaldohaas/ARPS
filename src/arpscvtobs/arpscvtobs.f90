!########################################################################
!########################################################################
!#########                                                      #########
!#########                 PROGRAM arpscvtobs                   #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

PROGRAM arpscvtobs

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Read in LSO file(s), and output observation data in HDF file(s)
! for use by NCL.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, January 2002.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! List include files.
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'adas.inc'
  INCLUDE 'grid.inc'

!-----------------------------------------------------------------------
!
! Declare external functions.
!
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: PT_WSYM ! GEMPAK function

!fpp$ expand (PT_WSYM)
!dir$ inline always PT_WSYM
!*$*  inline routine (PT_WSYM)

!-----------------------------------------------------------------------
!
! NAMELIST variables.
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=256) :: blacklistfile
  INTEGER :: n_lso_files
  CHARACTER(LEN=256) :: lsoinfile(mx_sng_file) ! Set in adas.inc
  CHARACTER(LEN=256) :: lsohdfoutfile(mx_sng_file)
  NAMELIST /lsoobs/ blacklistfile,n_lso_files, lsoinfile,lsohdfoutfile


!-----------------------------------------------------------------------
!
! Declare variables used with LSO file.
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: infile,blackfile
  INTEGER, PARAMETER :: maxsta = mx_sng ! Set in adas.inc
  INTEGER :: istart
  CHARACTER (LEN=24) :: atime
  INTEGER :: n_meso_g,n_meso_pos,n_sao_g,n_sao_pos_g
  INTEGER :: n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g
  INTEGER :: n_obs_b,n_obs_pos_b
  CHARACTER (LEN=5) :: stn(maxsta)
  CHARACTER (LEN=8) :: obstype(maxsta)
  REAL :: lat(maxsta),lon(maxsta),elev(maxsta)
  CHARACTER (LEN=8) :: wx(maxsta)
  REAL :: t(maxsta),td(maxsta)
  REAL :: dd(maxsta),ff(maxsta)
  REAL :: ddg(maxsta),ffg(maxsta)
  REAL :: pstn(maxsta),pmsl(maxsta),alt(maxsta)
  INTEGER :: kloud(maxsta)
  REAL :: ceil(maxsta),lowcld(maxsta),cover(maxsta),rad(maxsta)
  INTEGER :: idp3(maxsta)
  CHARACTER (LEN=1) :: store_emv(maxsta,5)
  CHARACTER (LEN=4) :: store_amt(maxsta,5)
  REAL :: store_hgt(maxsta,5)
  REAL :: vis(maxsta)
  INTEGER :: obstime(maxsta)

  REAL :: u(maxsta),v(maxsta)
  REAL :: ug(maxsta),vg(maxsta)
  REAL :: xsta(maxsta),ysta(maxsta)

  CHARACTER*5, PARAMETER :: elevunits = 'm'
  CHARACTER*5, PARAMETER :: tunits = 'F'
  CHARACTER*5, PARAMETER :: tdunits = 'F'
  CHARACTER*5, PARAMETER :: ddunits = 'Deg'
  CHARACTER*5, PARAMETER :: ffunits = 'kts'
  CHARACTER*5, PARAMETER :: ddgunits = 'Deg'
  CHARACTER*5, PARAMETER :: ffgunits = 'kts'
  CHARACTER*5, PARAMETER :: pmslunits = 'mb'
  CHARACTER*5, PARAMETER :: ceilunits = 'm'
  CHARACTER*5, PARAMETER :: lowcldunits = 'm'
  CHARACTER*5, PARAMETER :: visunits = 'mi'

  CHARACTER(LEN=5) :: tmp_stn
  CHARACTER(LEN=8) :: tmp_obstype
  CHARACTER(LEN=1) :: stn2d(5,maxsta)
  CHARACTER(LEN=1) :: obstype2d(8,maxsta)

!-----------------------------------------------------------------------
!
! Misc. internal variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: nhisfile_max=200

  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)  
  INTEGER :: nhisfile  

  INTEGER :: hinfmt, nchin, istatus
  INTEGER :: lengbf, lenfil, ireturn
  REAL :: time

  CHARACTER(LEN=256) :: filename

  INTEGER :: i,j,k
  INTEGER :: ifile
  REAL :: latnot(2)
  REAL :: xctr,yctr,x0,y0
  INTEGER :: iob
  REAL :: ddrot
  
  INTEGER :: yrob(maxsta),monthob(maxsta),dayob(maxsta),hrob(maxsta),    &
             minob(maxsta)
  INTEGER :: iwx(maxsta)
  CHARACTER*3 :: cmonth
  INTEGER :: iyr,imonth,iday,ihr,imin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Read NAMELIST input file.
!
!-----------------------------------------------------------------------

  READ(5,lsoobs,END=999)
  WRITE(6,'(/a)') 'Namelist block lsoobs successfully read.'

!-----------------------------------------------------------------------
!
! Loop through LSO files.
!
!-----------------------------------------------------------------------

  DO ifile = 1,n_lso_files

!-----------------------------------------------------------------------
!
!   Read in LSO file.
!
!-----------------------------------------------------------------------

    infile = TRIM(lsoinfile(ifile))

    WRITE(6,*)'Opening ',infile

    blackfile = TRIM(blacklistfile)

    istart = 1
    CALL read_surface_obs(infile,blackfile,                              &
           maxsta,istart,atime,n_meso_g,                                 &
           n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,   &
           n_obs_pos_g,n_obs_b,n_obs_pos_b,stn,obstype,lat,lon,elev,wx,  &
           t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad, &
           idp3,store_emv,store_amt,store_hgt,vis,obstime,istatus)

!-----------------------------------------------------------------------
!
!   Break up the character array atime into integer values for year, 
!   month, day, hour, and minute.
!
!-----------------------------------------------------------------------

    cmonth = atime(4:6)

    SELECT CASE (cmonth)
    CASE('JAN')
      imonth = 1
    CASE('FEB')
      imonth = 2
    CASE('MAR')
      imonth = 3
    CASE('APR')
      imonth = 4
    CASE('MAY')
      imonth = 5
    CASE('JUN')
      imonth = 6
    CASE('JUL')
      imonth = 7
    CASE('AUG')
      imonth = 8
    CASE('SEP')
      imonth = 9
    CASE('OCT')
      imonth = 10
    CASE('NOV')
      imonth = 11
    CASE('DEC')
      imonth = 12
    CASE DEFAULT
      WRITE(6,*)'ERROR:  Illegal month.'
      WRITE(6,*)'cmonth = ',cmonth
      WRITE(6,*)'Aborting...'
      STOP
    END SELECT

    READ(atime(8:11),'(i4)') iyr
    READ(atime(1:2),'(i2)') iday
    READ(atime(13:14),'(i2)') ihr
    READ(atime(16:17),'(i2)') imin

!-----------------------------------------------------------------------
!
!   Loop through all observations.
!
!-----------------------------------------------------------------------

    DO iob = 1,n_obs_b

!-----------------------------------------------------------------------
!
!     Convert character "present weather" code to numerical WMO code.
!     This involves the use of a GEMPAK function.
!
!-----------------------------------------------------------------------

      iwx(iob) = PT_WSYM(wx(iob))

!TODO:  Add code to check for "alternative" WMO weather symbols.

!      WRITE(6,*)'iob: ',iob,' wx: ',wx(iob),' iwx: ',iwx(iob)
 
!-----------------------------------------------------------------------
!
!     Assign year, month, day, hour, and minute for each observation.
!
!-----------------------------------------------------------------------

      yrob(iob) = iyr
      monthob(iob) = imonth
      dayob(iob) = iday
      hrob(iob) = ihr
      minob(iob) = imin

!-----------------------------------------------------------------------
!
!     Convert station ID's to 2 dimensional array for NCL.
!
!-----------------------------------------------------------------------

      tmp_stn = stn(iob)
      READ(tmp_stn,('(1A,1A,1A,1A,1A)')) stn2d(1,iob), stn2d(2,iob), &
                                         stn2d(3,iob), stn2d(4,iob), &
                                         stn2d(5,iob)


!-----------------------------------------------------------------------
!
!     Convert station types to 2 dimensional array for NCL.
!
!-----------------------------------------------------------------------

      tmp_obstype = obstype(iob)
      READ(tmp_obstype,('(1A,1A,1A,1A,1A,1A,1A,1A)')) &
        obstype2d(1,iob), obstype2d(2,iob), obstype2d(3,iob), &
        obstype2d(4,iob), obstype2d(5,iob), obstype2d(6,iob), &
        obstype2d(7,iob), obstype2d(8,iob)

!-----------------------------------------------------------------------
!
!     Change missing flags.
!
!-----------------------------------------------------------------------

      IF (pmsl(iob).EQ.-99.9) THEN
        pmsl(iob) = -9999.
      END IF
      IF (t(iob).EQ.-99.9) THEN
        t(iob) = -9999.
      END IF
      IF (td(iob).EQ.-99.9) THEN
        td(iob) = -9999.
      END IF
      IF (dd(iob).EQ.-100.) THEN
        dd(iob) = -9999.
      END IF
      IF (ff(iob).EQ.-100.) THEN
        ff(iob) = -9999.
      END IF
      IF (ddg(iob).EQ.-100.) THEN
        ddg(iob) = -9999.
      END IF
      IF (ffg(iob).EQ.-100.) THEN
        ffg(iob) = -9999.
      END IF

      IF (cover(iob).EQ.-99.9) THEN
        cover(iob) = -9999.
      END IF
      IF (kloud(iob).GT.0) THEN
        DO i = 1,kloud(iob)
          IF (store_amt(iob,i) == ' X  ') THEN
            cover(iob) = -1.0 ! Obscured flag for NCL.
          ELSE IF (store_amt(iob,i) == '-X  ') THEN
            cover(iob) = -1.0 ! Obscured flag for NCL.
          END IF
        END DO
      END IF

      IF (elev(iob).EQ.-99.9) THEN
        elev(iob) = -9999.
      END IF
      IF (ceil(iob).EQ.-99.9) THEN
        ceil(iob) = -9999.
      END IF
      IF (lowcld(iob).EQ.-99.9) THEN
        lowcld(iob) = -9999.
      END IF
      IF (vis(iob).EQ.-99.9) THEN
        vis(iob) = -9999.
      END IF

    END DO ! DO iob = 1,n_obs_b

!-----------------------------------------------------------------------
!
!   Write HDF file with LSO data.
!
!-----------------------------------------------------------------------

    WRITE(6,*)'EMK: lsohdfoutfile(ifile) = ',TRIM(lsohdfoutfile(ifile))

    CALL wrtlsohdf(maxsta,n_obs_b,yrob,monthob,dayob,hrob,minob,&
                   stn2d,obstype2d,lat,lon, &
                   elev,elevunits, &
                   iwx, &
                   t,tunits, &
                   td,tdunits, &
                   dd, ddunits, &
                   ff, ffunits, &
                   ddg, ddgunits, &
                   ffg, ffgunits, &
                   pmsl, pmslunits,  &
                   cover, &
                   ceil, ceilunits, &
                   lowcld, lowcldunits, &
                   vis, visunits, &
                   lsohdfoutfile(ifile))

   END DO ! DO ifile = 1,n_lso_file

!-----------------------------------------------------------------------
!
!  The end.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'Normal termination of arpscvtobs.'
  STOP

!-----------------------------------------------------------------------
!
! The end.
!
!-----------------------------------------------------------------------

  999 CONTINUE
  WRITE(6,*)'Problem reading NAMELIST block!'
  WRITE(6,*)'Aborting...'
  STOP

END PROGRAM arpscvtobs


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
  CHARACTER (LEN=256) :: infile
  CHARACTER (LEN=256) :: blackfile
  CHARACTER (LEN=24)  :: atime
  CHARACTER (LEN=5) :: stn(maxsta)
  CHARACTER (LEN=8) :: obstype(maxsta)
  CHARACTER (LEN=8) :: wx(maxsta)
  CHARACTER (LEN=1) :: store_emv(maxsta,5)
  CHARACTER (LEN=4) :: store_amt(maxsta,5)
!
  INTEGER :: n_meso_g,n_meso_pos,n_sao_g,n_sao_pos_g
  INTEGER :: n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g
  INTEGER :: n_obs_b,n_obs_pos_b
  INTEGER :: istatus
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
!
!.....  Error trapping for too many stations for array size.
!
  IF(n_obs_b > maxsta) THEN
    PRINT 990, maxsta,n_obs_b,atime
    990       FORMAT(' +++ ERROR in READ_SURFACE_OBS: maxstns = ',i8,/, &
           ' but there are ',i8,' stations in the ',a24,' obs file.',/)
    PRINT *,'    Increase the value of "maxstns" and try again.'
    istatus = -2
    RETURN
  END IF
!
!.....  Now read the station data.
!
  DO i=1,n_obs_b
    k=i+istart-1
    READ(1,901) stn(k),lat(k),lon(k),elev(k),obstype(k),                &
        obstime(k),wx(k)
    901       FORMAT(1X,a5,f6.2,1X,f7.2,1X,f5.0,1X,a8,1X,i4,1X,a8)
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
    DO ii=1,5
      store_emv(k,ii) = ' '
      store_amt(k,ii) = '    '
      store_hgt(k,ii) = badflag
    END DO ! ii
    IF(kloud(k) > 0) THEN
      DO ii=1,kloud(k)
        READ(1,907) store_emv(k,ii),store_amt(k,ii),store_hgt(k,ii)
        907           FORMAT(5X,a1,1X,a4,1X,f7.1)
      END DO !ii
    END IF
!
  END DO !k
!
  REWIND(1)
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
  REWIND(1)
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
  CHARACTER (LEN=256) :: blkfile
  INTEGER :: maxsta,istart,nobs,misval
  REAL :: rmisval
  CHARACTER (LEN=5) :: stn(maxsta)
  REAL :: t(maxsta),td(maxsta)
  REAL :: dd(maxsta),ff(maxsta),ddg(maxsta),ffg(maxsta)
  REAL :: pstn(maxsta),pmsl(maxsta),alt(maxsta)
  REAL :: store_hgt(maxsta,5)
  REAL :: ceil(maxsta),lowcld(maxsta),cover(maxsta)
  REAL :: vis(maxsta),rad(maxsta)
!
  INTEGER*4   kloud(maxsta),idp3(maxsta)
!
  CHARACTER (LEN=8) :: wx(maxsta)
  CHARACTER (LEN=8) :: obstype(maxsta)
  CHARACTER (LEN=1) :: store_emv(maxsta,5)
  CHARACTER (LEN=4) :: store_amt(maxsta,5)
  INTEGER :: istatus
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
  INTEGER :: i,item,ista,ivar,nblack
!
  OPEN(12,FILE=trim(blkfile),ERR=500,STATUS='old')
  DO item=1,maxsta
    READ(12,810,ERR=520,END=550) blkstn,nblack
    810   FORMAT(a5,i5)
    DO i=1,maxblk
      blkvar(i)=-99
    END DO
    READ(12,*,ERR=240,END=240) (blkvar(i),i=1,nblack)
    240   CONTINUE
    PRINT *, ' blkstn: ',blkstn
    PRINT *, ' blkvar: ',blkvar
    DO i=1,nobs
      ista=i+istart-1
      IF(stn(ista) == blkstn) THEN
        DO ivar=1,nblack
          IF(blkvar(ivar) < 0) GOTO 400
          GO TO (301,302,303,304,305,306,307) blkvar(ivar)
          301         CONTINUE
          WRITE(6,831) blkvar(ivar),blkstn
          831           FORMAT(' Blacklisting Temp variable(s) ',i4,' at ',a5)
          t(ista)=rmisval
          CYCLE
          302         CONTINUE
          WRITE(6,832) blkvar(ivar),blkstn
          832           FORMAT(' Blacklisting Dew Pt variable(s) ',i4,' at ',a5)
          td(ista)=rmisval
          CYCLE
          303         CONTINUE
          WRITE(6,833) blkvar(ivar),blkstn
          833           FORMAT(' Blacklisting Wind variable(s) ',i4,' at ',a5)
          dd(ista)=rmisval
          ff(ista)=rmisval
          ddg(ista)=rmisval
          ffg(ista)=rmisval
          CYCLE
          304         CONTINUE
          WRITE(6,834) blkvar(ivar),blkstn
          834           FORMAT(' Blacklisting Pres variable(s) ',i4,' at ',a5)
          pstn(ista)=rmisval
          pmsl(ista)=rmisval
          alt(ista)=rmisval
          CYCLE
          305         CONTINUE
          WRITE(6,835) blkvar(ivar),blkstn
          835           FORMAT(' Blacklisting Cloud variable(s) ',i4,' at ',a5)
          kloud(ista)=0
          ceil(ista)=rmisval
          lowcld(ista)=rmisval
          cover(ista)=rmisval
          CYCLE
          306         CONTINUE
          WRITE(6,836) blkvar(ivar),blkstn
          836           FORMAT(' Blacklisting Vis variable(s) ',i4,' at ',a5)
          vis(ista)=rmisval
          CYCLE
          307         CONTINUE
          WRITE(6,837) blkvar(ivar),blkstn
          837           FORMAT(' Blacklisting Solar variable(s) ',i4,' at ',a5)
          rad(ista)=rmisval
          CYCLE
        END DO
        GOTO 400
      END IF
    END DO
    WRITE(6,830) blkstn
    830   FORMAT(' Blacklisted station ',a5,' not found in dataset')
400 CONTINUE
  END DO
  CLOSE (12)
  RETURN
  500 CONTINUE
  CLOSE (12)
  RETURN
  520 CONTINUE
  CLOSE (12)
  RETURN
  550 CONTINUE
  CLOSE (12)
  RETURN
END SUBROUTINE blklistsfc

