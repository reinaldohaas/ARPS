!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE PREPSGLOBS                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE prepsglobs(mxsng,ntime,srcsng,nsrc_sng,sngfile,              &
           stnsng,latsng,lonsng,xsng,ysng,hgtsng,                       &
           obsng,qobsng,qualsng,isrcsng,qsrc,chsrc,                     &
           nx,ny,nz,nvar,anx,xs,ys,zp,                                  &
           shght,su,sv,spres,stheta,sqv,                                &
           rmiss,nobsng,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine reads non-surface single-level data.
! The data are appended to surface data already collected.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Keith Brewster, CAPS
! October, 1997
!
! MODIFICATION HISTORY:
!
! 5/4/98 Developed based on RDRASS to assign height from background
!        Can also be used to read satellite winds (John Manobianco).
!
! 1998/07/13 (Gene Bassett) Implemented general format for sgl obs
! and used this for MDCRS data.
!
! 1999/09/16 (Dingchen Hou) Modified the interpretion of MDCRS data.
!
! 1999/11/19 (Gene Bassett, Keith Brewster)
!   Modified MDCRS data to use pressure only for deriving other
!   observational quantities but not as a pressure observation in
!   ADAS proper.
!
! 11/9/2001 (Keith Brewster)
!   Added FORTRAN-90 features and rearranged the logic to eliminate
!   GO-TOs, added processing for humidity and pressure data,
!   and generally make it more suitable for additional data 
!   sources other than aircraft data.
!
! August 2002 (Eric Kemp)
!   Fixed bug involving check of contents of chsng array.  Further
!   reorganizations/changes/updates.
!
! December 2005 (Kevin W. Thomas)
!   MPI update.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: mxsng       ! Max possible number of
                                     ! single-point observations.
  INTEGER, INTENT(IN) :: nvar        ! Total number of observed
                                     ! variable types.
  INTEGER, INTENT(IN) :: ntime       ! Maximum number of time levels
                                     ! for observations.
  INTEGER, INTENT(IN) :: nsrc_sng    ! Max number of types of single
                                     ! point observations.
  CHARACTER(LEN=*), INTENT(IN) :: srcsng(nsrc_sng) ! Observation types
  CHARACTER(LEN=*), INTENT(IN) :: sngfile ! Name of SGL input file.
  CHARACTER(LEN=5), INTENT(INOUT) :: stnsng(mxsng) ! ID of observations
  REAL, INTENT(INOUT) :: latsng(mxsng) ! Observation latitude (deg N)   
  REAL, INTENT(INOUT) :: lonsng(mxsng) ! Observation longitude (deg E)
  REAL, INTENT(INOUT) :: xsng(mxsng)   ! Observation x grid coordinate (m)
  REAL, INTENT(INOUT) :: ysng(mxsng)   ! Observation y grid coordinate (m)
  REAL, INTENT(INOUT) :: hgtsng(mxsng) ! Observation height (m MSL)
  REAL, INTENT(INOUT) :: obsng(nvar,mxsng) ! Observation data by variable.
  REAL, INTENT(INOUT) :: qobsng(nvar,mxsng) ! Observation standard error.
  INTEGER, INTENT(INOUT) :: qualsng(nvar,mxsng) ! Observation quality
                                                ! indicator.   
  INTEGER, INTENT(INOUT) :: isrcsng(mxsng) ! Source number for each
                                           ! station.
  REAL, INTENT(IN) :: qsrc(nvar,mxsng) ! Standard error of observation
                                       ! for each source.
  CHARACTER(LEN=8), INTENT(INOUT) :: chsrc(mxsng,ntime) ! Source for each
                                                        ! station.
  INTEGER, INTENT(IN) :: nx,ny,nz ! Grid dimensions.
  REAL, INTENT(IN) :: anx(nx,ny,nz,nvar) ! Gridded analysis data.
  REAL, INTENT(IN) :: xs(nx) ! x-coordinates of grid scalar points (m)  
  REAL, INTENT(IN) :: ys(ny) ! y-coordinates of grid scalar points (m)
  REAL, INTENT(IN) :: zp(nx,ny,nz) ! Height of grid w-points (m MSL)
  REAL, INTENT(INOUT) :: shght(nz)  ! Dummy array for height
  REAL, INTENT(INOUT) :: su(nz)     ! Dummy array for u-winds
  REAL, INTENT(INOUT) :: sv(nz)     ! Dummy array for v-winds
  REAL, INTENT(INOUT) :: spres(nz)  ! Dummy array for pressure
  REAL, INTENT(INOUT) :: stheta(nz) ! Dummy array for pot. temperature
  REAL, INTENT(INOUT) :: sqv(nz)    ! Dummy array for specific humidity
  REAL, INTENT(IN) :: rmiss     ! Flag for missing data.
  INTEGER, INTENT(INOUT) :: nobsng ! Total number of single-point
                                   ! observations stored in arrays
                                   ! (nobsng <= maxsng)  
  INTEGER, INTENT(INOUT) :: istatus ! Status flag

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Internal parameters
!
!-----------------------------------------------------------------------

  REAL, PARAMETER :: ktstoms=0.5144
  REAL, PARAMETER :: mbtopa=100.

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: jsrc
  INTEGER :: nprev
  CHARACTER (LEN=132) :: line

  INTEGER :: ivar,ksta,nsingle
  INTEGER :: ipt,jpt,kgrid
  INTEGER :: indom
  INTEGER :: nlevs

  CHARACTER (LEN=8) :: obs_id
  INTEGER :: obs_num
  CHARACTER (LEN=8) :: date_str
  CHARACTER (LEN=6) :: time_str
  INTEGER :: num_lines
  CHARACTER (LEN=8) :: obs_type

  REAL :: obs_lat,obs_lon,hgt_msl,hgt_agl
  REAL :: obs_t,obs_td,obs_wspd,obs_wdir,obs_gust
  REAL :: obs_psta,obs_pmsl,obs_palt,obs_vis
  REAL :: xsta,ysta,tk,tdk,psta,qv,qvsat
  CHARACTER (LEN=16) :: obs_wx

  CHARACTER (LEN=6) :: fmt_in

  INTEGER :: k,istat
  REAL :: weight
  REAL :: ddrot,selev

  REAL :: altpres   ! pressure derived from altimeter reading
  INTEGER :: iunit

!-----------------------------------------------------------------------
!
!  External functions
!
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: alttostpr
  REAL, EXTERNAL :: msltostpr
  REAL, EXTERNAL :: f_qvsat
  REAL, EXTERNAL :: ztopsa

!fpp$ expand (alttostpr)
!!dir$ inline always alttostpr
!*$*  inline routine (alttostpr)

!fpp$ expand (msltostpr)
!!dir$ inline always msltostpr
!*$*  inline routine (msltostpr)

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!fpp$ expand (ztopsa)
!!dir$ inline always ztopsa
!*$*  inline routine (ztopsa)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL getunit(iunit)

  nprev = nobsng
  WRITE(6,'(A,A)') ' Opening: ',sngfile
  WRITE(6,'(A,I6,I6)') ' In PREPSGLOBS, mxsng, nprev = ',              &
                         mxsng,nprev
  OPEN(UNIT=iunit,FILE=trim(sngfile),STATUS='old',iostat=istat)
  IF ( istat /= 0 ) THEN
    WRITE(6,'(A)') ' Error opening file: ',sngfile
    istatus = -1
    RETURN
  END IF

!-----------------------------------------------------------------------
!
! Main data-reading loop
!
!-----------------------------------------------------------------------

  ksta = nobsng

  DO ! Go line by line in file until end of file is reached

!-----------------------------------------------------------------------
!
!   INPUT FORMAT:
!
!    first line:     obs_id, obs_num, date_str, time_str,
!                    obs_lat, obs_lon, hgt_msl, hgt_agl,
!                    obs_type, num_lines
!
!     See MDCRS format below.
!
!         obs_id     station id (8 characters)
!         obs_num    station id number
!         date_str   yyyymmdd
!         time_str   hhmm
!         obs_lat    latidute (degrees)
!         obs_lon    longitude (degrees)
!         hgt_msl    height above mean sea level (m)
!                       value considered missing if <= rmiss
!         hgt_agl    height above ground level (m), optional,
!                       intended for use with surface obs to help with
!                       mis-matches in terrain
!                       value considered missing if <= rmiss
!         obs_type   ADAS observation type (8 characters)
!         num_lines  number of observation lines to follow
!                       the header line
!
!         currently unused in ADAS:
!            obs_id, obs_num, date_str, time_str, hgt_agl
!
!
!    currently supported observation lines formatted as follows:
!
!  FMT1    20.0 -10.0 215  10.0 100.0 1013.2 1013.2 100. +TSRA
!  ______ _____ _____ ___ _____ _____ ______ ______ ____ ________________
!
!    FMT1            obs_t, obs_td, obs_wdir, obs_wspd, obs_gust,
!                    obs_pmsl, obs_palt, obs_vis, obs_wx
!
!         obs_t      temperature (C)
!                       value considered missing if <= rmiss
!         obs_td     dew point (C)
!                       value considered missing if <= rmiss
!         obs_wdir   wind direction (degrees clockwise from north)
!                       value considered missing if < 0
!         obs_wspd   wind speed (m/s)
!                       value considered missing if <= rmiss
!         obs_wdir   wind direction (degrees clockwise from north)
!                       value considered missing if < 0
!         obs_gust   wind gust (m/s)
!                       value considered missing if <= rmiss
!         obs_pmsl   sea level pressure (mb)
!                       value considered missing if <= rmiss
!         obs_palt   altimiter setting (mb)
!                       value considered missing if <= rmiss
!         obs_vis    visibility (km)
!                       value considered missing if < 0
!         obs_wx     current weather (METAR format)
!
!         currently unused in ADAS:
!            obs_gust, obs_palt, obs_vis, obs_wx
!
!-----------------------------------------------------------------------
!
!  Initially set all non-header line values to missing
!
!-----------------------------------------------------------------------

    obs_t = rmiss
    obs_td = rmiss
    obs_wdir = -1
    obs_wspd = rmiss
    obs_wdir = -1
    obs_gust = rmiss
    obs_psta = rmiss
    obs_pmsl = rmiss
    obs_palt = rmiss
    obs_vis = rmiss
    obs_wx = '                '

    altpres = rmiss
    psta = rmiss

!-----------------------------------------------------------------------
!
!   Read in the info for one station
!   Skip comments
!   Comments begin with # or !. Also, skip lines that begin with a space.
!
! Sample format
!XXXXXXXX 00000000 19990601_153600  42.7167  -84.9517   9448 -999.00 MDCRS     1
!-------- -------- -------- ------ -------- --------- ------ ------- -------- --
! A8,   1X, I8,1X,I4,2I2,1X,2I2,2X,1X,F8.4,1X,F9.4, 1X, I6,1X,F7.2,1X,A8,1X,I2
! MDCRS   -39.8 -999.0 267  28.3                            YYZ  ORD  FSL00307
!-----------------------------------------------------------------------

    READ (iunit,'(A)',IOSTAT=istat) line

    IF (istat /= 0) EXIT ! Get out of DO loop

    IF (line(1:1) == '#' .OR. line(1:1) == '!' .OR. line(1:1) == ' ') CYCLE

    READ(line,                                                              &
      '(a8,1x,i8,1x,a8,1x,a6,1x,f8.4,1x,f9.4,1x,f6.0,1x,f7.2,1x,a8,1x,i2)', &
      IOSTAT=istat)  obs_id, obs_num, date_str, time_str,                   &
      obs_lat, obs_lon, hgt_msl, hgt_agl, obs_type, num_lines
    IF (istat /= 0) CYCLE

    READ(iunit,                                                      &
      '(2x,a6,1x,f6.1,1x,f6.1,1x,f3.0,1x,f5.1)',IOSTAT=istat)           &
      fmt_in, obs_t, obs_td, obs_wdir, obs_wspd
    IF (istat /= 0) CYCLE

!    WRITE(6,'(2F9.2, F7.0, 2F7.1, F5.0, F5.1)')                         &
!      obs_lat, obs_lon, hgt_msl, obs_t, obs_td, obs_wdir, obs_wspd

    CALL lltoxy(1,1,obs_lat,obs_lon,xsta,ysta)
    CALL findlc(nx,ny,xs,ys,xsta,ysta,ipt,jpt,indom)

    IF (indom < 0 ) CYCLE

!-----------------------------------------------------------------------
!
!   If MDCRS data are being processed, we need to convert from
!   pressure altitude to pressure and then estimate the height
!   based on the background pressure field.
!
!-----------------------------------------------------------------------

    IF ( TRIM(obs_type) == 'MDCRS' ) THEN

      altpres = mbtopa*ztopsa(hgt_msl)

!-----------------------------------------------------------------------
!
!     In previous versions, all the MDCRS computations were here.  They
!     are now in subroutine PREPSGLMDCRS, later in this file. 
!
!     If we are an MPI program, we can't do the vertical interpolation
!     at this time, as only processor 0 is handling the I/O code.  We
!     will postpone the computation of "hgt_msl".
!
!     For the MPI case, PREPSGLMDCRS is called from PREPSGLMDCRS_SM
!
!     Note that it "appears" that "hgt_msl" is used later, however, it
!     really isn't.  There is a lot of dead code below (assumed for now
!     to be hooks by others for other data types to be added).
!
!-----------------------------------------------------------------------

       IF (mp_opt < 1 )                                                 &
         CALL  prepsglmdcrs(nx,ny,nz,nvar,anx,xs,ys,zp,                 &
             xsta,ysta,ipt,jpt,altpres,rmiss,                           &
             shght,su,sv,spres,stheta,sqv,                              &
             selev,nlevs,hgt_msl)

    END IF ! IF ( obs_type == 'MDCRS   ' ) THEN

!-----------------------------------------------------------------------
!
!   Use this datapoint if we have a valid height for it.
!
!-----------------------------------------------------------------------

    ksta = ksta + 1

    IF (ksta > mxsng) THEN
      WRITE(6,'(a,i6,a)')                                                &
          ' ERROR: number of single level stations,',mxsng,              &
          'exceeded in PREPSGLOBS.  Increase mxsng and retry.'
      istatus = -2
      CALL arpsstop("array size problem",1)
    END IF

    DO ivar=1,nvar
      obsng(ivar,ksta)=rmiss
      qobsng(ivar,ksta)=999.
      qualsng(ivar,ksta)=-99
    END DO

    stnsng(ksta) = obs_id(1:5)
    latsng(ksta) = obs_lat
    lonsng(ksta) = obs_lon
    xsng(ksta) = xsta
    ysng(ksta) = ysta

    chsrc(ksta,1) = obs_type
    DO jsrc=1,nsrc_sng
      IF (chsrc(ksta,1) == srcsng(jsrc)) THEN
        isrcsng(ksta)=jsrc
        EXIT
      END IF
    END DO

    jsrc=isrcsng(ksta)

    if ( jsrc == 0 ) then
        write(6,*)                                                      &
        'Unable to find data source ',chsrc(ksta,1),' in input file.'
        write(6,*)
        return
    endif

!-----------------------------------------------------------------------
!
!   Set obs height
!
!-----------------------------------------------------------------------

    hgtsng(ksta) = hgt_msl

!-----------------------------------------------------------------------
!
!   Station pressure (mb)
!   For now we'll convert msl pressure to station pressure
!   using the current temperature.
!
!-----------------------------------------------------------------------

    IF (obs_psta > rmiss ) THEN
      psta=mbtopa*obs_psta
      obsng(3,ksta)=psta
      qobsng(3,ksta)=qsrc(3,jsrc)
      qualsng(3,ksta)=100
    ELSE IF (obs_palt > rmiss ) THEN
      psta=mbtopa*alttostpr(obs_palt,hgt_msl)
      obsng(3,ksta)=psta
      qobsng(3,ksta)=qsrc(3,jsrc)
      qualsng(3,ksta)=100
    ELSE IF (obs_pmsl > rmiss .AND.                             &
             obs_t > rmiss .AND.                                &
             hgt_msl < 500. ) THEN
      tk=obs_t + 273.15
      psta=mbtopa*msltostpr(obs_pmsl,tk,hgt_msl)
      obsng(3,ksta)=psta
      qobsng(3,ksta)=qsrc(3,jsrc)
      qualsng(3,ksta)=100
    ELSE IF (altpres > rmiss) THEN ! MDCRS
      psta=altpres
    ELSE
      psta=mbtopa*alttostpr(1013.,hgt_msl)
    END IF

!-----------------------------------------------------------------------
!
!   Set winds.
!
!-----------------------------------------------------------------------

    IF (obs_wdir >= 0. .AND. obs_wspd >= 0.) THEN
      CALL ddrotuv(1,lonsng(ksta),obs_wdir,obs_wspd,ddrot,          &
                   obsng(1,ksta),obsng(2,ksta))
      qualsng(1,ksta)=100
      qualsng(2,ksta)=100
      qobsng(1,ksta)=qsrc(1,jsrc)
      qobsng(2,ksta)=qsrc(2,jsrc)
    END IF

!-----------------------------------------------------------------------
!
!   Potential temperature (K)
!
!-----------------------------------------------------------------------

    IF (obs_t > rmiss .AND. psta > rmiss) THEN
      tk=obs_t + 273.15
      obsng(4,ksta)=tk*((p0/psta)**rddcp)
      qualsng(4,ksta)=100
      qobsng(4,ksta)=qsrc(4,jsrc)
    END IF

!-----------------------------------------------------------------------
!
!   Specific humidity (kg/kg)
!
!-----------------------------------------------------------------------

    IF (obs_td > rmiss .AND.  obs_t > rmiss .AND.                    &
                                  psta > rmiss ) THEN
      tk=obs_t + 273.15
      tdk=obs_td + 273.15
      qv = f_qvsat( psta,tdk )
      qvsat = f_qvsat( psta, tk )
      obsng(5,ksta)=max((0.05*qvsat),min(qv,qvsat))
      qobsng(5,ksta)=qsrc(5,jsrc)
      qualsng(5,ksta)=100
    END IF
  END DO

!-----------------------------------------------------------------------
!
! Update observation counts, close up and exit.
!
!-----------------------------------------------------------------------

  nobsng = ksta
  nsingle = nobsng - nprev
  WRITE(6,'(a,i6,i6,i6)') ' PREPSGLOBS: nsingle,nprev,nobsng',          &
                           nsingle,nprev,nobsng
  istatus=0
  CLOSE(iunit)
  CALL retunit(iunit)
  RETURN
END SUBROUTINE prepsglobs

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE PREPSGLMDCRS               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE prepsglmdcrs(nx,ny,nz,nvar,anx,xs,ys,zp,                     &
           xsta,ysta,ipt,jpt,altpres,rmiss,                             &
           shght,su,sv,spres,stheta,sqv,                                &
           selev,nlevs,hgt_msl)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine computes the height of an aircraft (MDCRS) ob using the
! background data field.  It has been split from PREPSGLOBS.  For MPI,
! the computations need to be on each processor.  PREPSGLOBS only runs on
! processor 0.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! November 30, 2005
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nvar        ! Total number of observed
                                     ! variable types.
  INTEGER, INTENT(IN) :: nx,ny,nz ! Grid dimensions.
  REAL, INTENT(IN) :: anx(nx,ny,nz,nvar) ! Gridded analysis data.
  REAL, INTENT(IN) :: xs(nx) ! x-coordinates of grid scalar points (m)  
  REAL, INTENT(IN) :: ys(ny) ! y-coordinates of grid scalar points (m)
  REAL, INTENT(IN) :: zp(nx,ny,nz) ! Height of grid w-points (m MSL)
  REAL, INTENT(IN) :: xsta,ysta     ! Station coordinates (m)
  INTEGER, INTENT(IN) :: ipt,jpt    ! Index in grid
  REAL, INTENT(IN)    :: altpres    ! pressure derived from altimeter reading
  REAL, INTENT(IN)    :: rmiss      ! Flag for missing data.
  REAL, INTENT(INOUT) :: shght(nz)  ! Dummy array for height
  REAL, INTENT(INOUT) :: su(nz)     ! Dummy array for u-winds
  REAL, INTENT(INOUT) :: sv(nz)     ! Dummy array for v-winds
  REAL, INTENT(INOUT) :: spres(nz)  ! Dummy array for pressure
  REAL, INTENT(INOUT) :: stheta(nz) ! Dummy array for pot. temperature
  REAL, INTENT(INOUT) :: sqv(nz)    ! Dummy array for specific humidity
  REAL, INTENT(OUT)   :: selev      ! Not used
  INTEGER, INTENT(OUT):: nlevs      ! Number of levels of data
  REAL, INTENT(OUT)   :: hgt_msl    ! Height (goal of subroutine)


!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: kgrid

  INTEGER :: k
  REAL :: weight

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Interpolate background field (in the horizontal) 
! for the whole vertical column.
!
!-----------------------------------------------------------------------

  CALL colint(nx,ny,nz,nvar,                                        &
          xs,ys,zp,xsta,ysta,ipt,jpt,anx,                           &
          su,sv,stheta,spres,shght,sqv,selev,                       &
          nlevs)

!-----------------------------------------------------------------------
!
! Interpolate hgts in the vertical from background.
! Interpolation is linear in log-p.
!
!-----------------------------------------------------------------------

  kgrid = 0
  hgt_msl = rmiss 
  DO k=3,nlevs-1
    IF (spres(k) < altpres) THEN
      kgrid = k
      EXIT
    END IF
  END DO

  IF (kgrid > 0) THEN
    weight = (LOG(altpres)-LOG(spres(kgrid)))/                    &
             (LOG(spres(kgrid-1))-LOG(spres(kgrid)))
    hgt_msl = (1.-weight)*shght(kgrid)                            &
                + weight*shght(kgrid-1)
  END IF

  RETURN
END SUBROUTINE prepsglmdcrs

!########################################################################
!########################################################################
!#########                                                      #########
!#########             SUBROUTINE PREPSGLOBSMDCRS_SM            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE prepsglmdcrs_sm(mxsng,ntime,latsng,lonsng,                   &
           xsng,ysng,hgtsng,chsrc,                                      &
           nx,ny,nz,nvar,anx,xs,ys,zp,                                  &
           shght,su,sv,spres,stheta,sqv,                                &
           rmiss,nobsng,indexsng,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Prepare to compute the MDCRS heights on the small grid.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! December 1, 2005
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: mxsng       ! Max possible number of
                                     ! single-point observations.
  INTEGER, INTENT(IN) :: nvar        ! Total number of observed
                                     ! variable types.
  INTEGER, INTENT(IN) :: ntime       ! Maximum number of time levels
                                     ! for observations.
                                     ! point observations.
  REAL, INTENT(INOUT) :: latsng(mxsng) ! Observation latitude (deg N)   
  REAL, INTENT(INOUT) :: lonsng(mxsng) ! Observation longitude (deg E)
  REAL, INTENT(INOUT) :: xsng(mxsng)   ! Observation x grid coordinate (m)
  REAL, INTENT(INOUT) :: ysng(mxsng)   ! Observation y grid coordinate (m)
  REAL, INTENT(INOUT) :: hgtsng(mxsng) ! Observation height (m MSL)
  CHARACTER(LEN=8), INTENT(INOUT) :: chsrc(mxsng,ntime) ! Source for each
                                                        ! station.
  INTEGER, INTENT(IN) :: nx,ny,nz ! Grid dimensions.
  REAL, INTENT(IN) :: anx(nx,ny,nz,nvar) ! Gridded analysis data.
  REAL, INTENT(IN) :: xs(nx) ! x-coordinates of grid scalar points (m)  
  REAL, INTENT(IN) :: ys(ny) ! y-coordinates of grid scalar points (m)
  REAL, INTENT(IN) :: zp(nx,ny,nz) ! Height of grid w-points (m MSL)
  REAL, INTENT(INOUT) :: shght(nz)  ! Dummy array for height
  REAL, INTENT(INOUT) :: su(nz)     ! Dummy array for u-winds
  REAL, INTENT(INOUT) :: sv(nz)     ! Dummy array for v-winds
  REAL, INTENT(INOUT) :: spres(nz)  ! Dummy array for pressure
  REAL, INTENT(INOUT) :: stheta(nz) ! Dummy array for pot. temperature
  REAL, INTENT(INOUT) :: sqv(nz)    ! Dummy array for specific humidity
  REAL, INTENT(IN) :: rmiss     ! Flag for missing data.
  INTEGER, INTENT(INOUT) :: nobsng ! Total number of single-point
                                   ! observations stored in arrays
                                   ! (nobsng <= maxsng)  
  INTEGER, INTENT(IN) :: indexsng(nobsng) ! "Ob" owner
  INTEGER, INTENT(INOUT) :: istatus ! Status flag

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Internal parameters
!
!-----------------------------------------------------------------------

  REAL, PARAMETER :: mbtopa=100.

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=132) :: line

  INTEGER :: ksta
  INTEGER :: ipt,jpt
  INTEGER :: indom
  INTEGER :: nlevs

  REAL :: hgt_msl,selev
  REAL :: xsta,ysta,tk,tdk,psta,qv,qvsat

  INTEGER :: k,istat

  REAL :: altpres   ! pressure derived from altimeter reading

!-----------------------------------------------------------------------
!
!  External functions
!
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: ztopsa

!fpp$ expand (ztopsa)
!!dir$ inline always ztopsa
!*$*  inline routine (ztopsa)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ksta = nobsng

  DO ksta=1,nobsng

!
!   Skip if we don't "own" the ob.
!

    IF (indexsng(ksta) .ne. myproc) CYCLE

!
!   Make sure we're only doing MDCRS data.
!

    IF (chsrc(ksta,1) .ne. "MDCRS   ") CYCLE

!
!   xsng/ysng were computed in SNGPROCESS.
!
!   The "indom" value isn't important, as we've already guaranteed we're
!   in the domain for this processor.  We didn't save ipt/jpt, so we need
!   to recompute.
!

    CALL findlc(nx,ny,xs,ys,xsng(ksta),ysng(ksta),ipt,jpt,indom)

    hgt_msl = hgtsng(ksta)
    altpres = mbtopa*ztopsa(hgt_msl)

    CALL  prepsglmdcrs(nx,ny,nz,nvar,anx,xs,ys,zp,                      &
        xsng(ksta),ysng(ksta),ipt,jpt,altpres,rmiss,                    &
        shght,su,sv,spres,stheta,sqv,                                   &
        selev,nlevs,hgt_msl)

    hgtsng(ksta) = hgt_msl

  END DO

  RETURN
END SUBROUTINE prepsglmdcrs_sm
