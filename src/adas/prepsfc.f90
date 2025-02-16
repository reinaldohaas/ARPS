!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PREPSFCOBS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE prepsfcobs  (ntime,mxstn,                                    &
           nvar_ob,nvar_anx,nsrc_sng,nobs,ipres,iptmp,iqv,              &
           sfcfile,tmchkfile,blackfile,                                 &
           var_ob,var_anx,name_src,qsrc,                                &
           rmiss,iqspr,sprdist,nobsrd,obstime,                          &
           stn,chsrc,isrc,latsta,lonsta,elevsta,xsta,ysta,              &
           nxlg,nylg,xslg,yslg,                                         &
           wx,kloud,idp3,store_emv,store_hgt,store_amt,                 &
           obsrd,obs,rely,qobs,qual,ival,climin,climax,                 &
           range,klim,wlim,qcthr,barqclim,                              &
           knt,wgtsum,zsum,sqsum,iqclist,iqcsave,jstatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Accumulate surface data and prepare them for analysis.
!  Includes reading of data, temporal and horizontal
!  consistency checks, and conversion of units.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  March, 1994
!
!  MODIFICATION HISTORY:
!  July, 1995  (Keith Brewster, CAPS)
!  ARPSanx version
!
!  Jan, 1996 (KB)
!  Updated documentation.
!
!  Aug, 1997 (KB)
!  Modified for new handling of observation error, qsrc.
!  QC is now done on analysis variables.
!
!  Nov, 1997 (KB)
!  Added argument rhmax, tmchkfile, for reorganized ADAS.
!
!  Mar, 1998 (KB)
!  Added argument blackfile
!
!  May, 2004 (KB)
!  Added arguments ipres,iptmp,iqv and barqclim to pass-on to barqc.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  ntime      Number of data times to process (normally 2)
!  mxstn      Size of data storage arrays
!  nvar_ob    Number of observed variables
!  nvar_anx   Number of analysis variables
!  nsrc_sng   Number of single-level data sources
!  sfcfile    Name of surface data file
!  tmchkfile  Name of surface data file used for time consistency check
!  blackfile  Name file containing blacklisted surface stations
!  var_ob     Names of observation variables
!  var_anx    Names of analysis variables
!  name_src   Names of data sources
!  qsrc       Standard error of observation for each source
!  rmiss      Missing value flag
!  iqspr      Super-observation flag
!  sprdist    Super-observation distance
!  nobsrd     Number of observations read-in
!  nobs       Number of observations output
!  ipres      Index of pressure variable
!  iptmp      Index of potential temperature variable
!  iqv        Index of specific humidity variable
!  obstime    Observation time
!  stn        Station name
!  chsrc      Source for each station
!  isrc       Source number for each station
!  latsta     Latitude of each station
!  lonsta     Longitude of each station
!  elevsta    Elevation of each station
!  xsta       X-coordinate of each station
!  ysta       Y-coordinate of each station
!  nxlg       Number of points in the "x-direction", global domain
!  nylg       Number of points in the "y-direction", global domain
!  xslg       X-location of scalar points (m), global domain
!  yslg       Y-location of scalar points (m), global domain
!  wx         Weather string
!  kloud      Number of cloud levels
!  idp3
!  store_emv  For each cloud level
!  store_hgt  Height of each cloud level
!  store_amt  Cloud amt for each cloud level
!  obsrd      Observations as read
!  obs        Analysis observations
!  rely       Reliability as assigned by time-consistency check
!  qobs       Observation standard error
!  qual       Observation quality indicator
!  ival
!  climin     Climatic mininum (used to QC observations)
!  climax     Climatic maxinum (used to QC observations)
!  range      Analysis range -- used for Barnes QC check
!  klim       Minimum number of stations in Barnes QC check
!  wlim       Minimum weight in Barnes QC
!  qcthr      Quality control threshold for each variable and source
!  knt        Temporary array used in Barnes QC routine
!  wgtsum     Temporary array used in Barnes QC routine
!  zsumc      Temporary array used in Barnes QC routine
!  sqsum      Temporary array used in Barnes QC routine
!  iqclist    Unit number/Flag for printing list of rejected obs
!  iqcsave    Unit number/Flag for printing info on QC processing
!  jstatus    Return status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: ntime,mxstn
  INTEGER :: nvar_ob,nvar_anx,nsrc_sng
  INTEGER :: nobs
  INTEGER :: ipres,iptmp,iqv
  INTEGER :: nxlg, nylg

  CHARACTER (LEN=*) :: sfcfile
  CHARACTER (LEN=*) :: tmchkfile
  CHARACTER (LEN=*) :: blackfile

  CHARACTER (LEN=6) :: var_ob(nvar_ob),var_anx(nvar_anx)
  CHARACTER (LEN=8) :: name_src(nsrc_sng)
  REAL :: qsrc(nvar_anx,nsrc_sng)
  REAL :: rmiss
  INTEGER :: iqspr
  REAL :: sprdist

  INTEGER :: nobsrd(ntime)
  CHARACTER (LEN=8) :: chsrc(mxstn,ntime)
  INTEGER :: isrc(mxstn)
  REAL :: latsta(mxstn,ntime),lonsta(mxstn,ntime)
  REAL :: elevsta(mxstn,ntime)
  REAL :: xsta(mxstn),ysta(mxstn)
  REAL :: xslg(nxlg),yslg(nylg)
  INTEGER :: obstime(mxstn,ntime)
  CHARACTER (LEN=5) :: stn(mxstn,ntime)
  CHARACTER (LEN=8) :: wx(mxstn,ntime)
  CHARACTER (LEN=1) :: store_emv(mxstn,5,ntime)
  CHARACTER (LEN=4) :: store_amt(mxstn,5,ntime)
  INTEGER :: kloud(mxstn,ntime),idp3(mxstn,ntime)
  REAL :: store_hgt(mxstn,5,ntime)
  REAL :: obsrd(mxstn,nvar_ob,ntime)
  INTEGER :: rely(mxstn,nvar_ob,ntime)
  REAL :: qobs(nvar_anx,mxstn)
  REAL :: obs(nvar_anx,mxstn)
  INTEGER :: qual(nvar_anx,mxstn)
  INTEGER :: ival(mxstn)
  REAL :: climin(nvar_ob),climax(nvar_ob)

  REAL    :: range,wlim
  INTEGER :: klim
  REAL    :: qcthr(nvar_anx,nsrc_sng)
  REAL    :: barqclim(nvar_anx,nsrc_sng)

  INTEGER :: knt(nvar_anx)
  REAL :: wgtsum(nvar_anx)
  REAL :: zsum(nvar_anx)
  REAL :: sqsum(nvar_anx)
!
  INTEGER :: iqclist,iqcsave
  INTEGER :: jstatus
!
!-----------------------------------------------------------------------
!
!  Misc quality parameters
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: sprdzcst=50.
  INTEGER, PARAMETER :: qcflag=-199
  INTEGER, PARAMETER :: qcsuper=200
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!

  CHARACTER (LEN=24) :: atime
  REAL :: ddrot,sprdzlim
  INTEGER :: n_meso_g,                                                  &
      n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,       &
      n_obs_pos_g,n_obs_b,n_obs_pos_b
  INTEGER :: nxt,iob,ivar,j,istatus
  INTEGER :: ipt,jpt,indom
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Note that the variable names are in var_ob and are
!  reproduced here:
!  1. t    2. td   3. dd  4. ff  5. u_ear  6. v_ear  7. pstn
!  8. pmsl 9. alt 10. ceil  11. hlow  12. cover 13. solar
!  14. vis
!
!  5 and 6 are actually read in us wind gust data, but
!  are converted to u,v in the grid orientation by the loop
!  after the read call.
!
!
!-----------------------------------------------------------------------
!
  jstatus = 1

  IF(strhopt > 0 ) THEN
    sprdzlim=min((2.*dzmin),dz,sprdzcst)
  ELSE
    sprdzlim=min(dz,sprdzcst)
  END IF
!
  nxt=nobs+1
  WRITE(6,'(1x,2a)') 'Getting surface data from: ',TRIM(sfcfile)
  CALL read_surface_obs(sfcfile,blackfile,mxstn,nxt,atime,n_meso_g,     &
      n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,       &
      n_obs_pos_g,n_obs_b,n_obs_pos_b,stn(1,1),chsrc(1,1),              &
      latsta(1,1),lonsta(1,1),elevsta(1,1),wx(1,1),                     &
      obsrd(1,1,1),obsrd(1,2,1),obsrd(1,3,1),obsrd(1,4,1),              &
      obsrd(1,5,1),obsrd(1,6,1),obsrd(1,7,1),obsrd(1,8,1),              &
      obsrd(1,9,1),  kloud(1,1),obsrd(1,10,1),obsrd(1,11,1),            &
      obsrd(1,12,1),obsrd(1,13,1),idp3(1,1),store_emv(1,1,1),           &
      store_amt(1,1,1),store_hgt(1,1,1),obsrd(1,14,1),                  &
      obstime(1,1),istatus)

  nobsrd(1)=n_obs_b
!
  IF(istatus /= 1) THEN       ! surface obs not available
    jstatus = 0
    WRITE(6,'(1x,2a)')     &
    '++Bad status return trying to get surface data from: ',TRIM(sfcfile)
    RETURN
  END IF
!
  WRITE(6,'(2a)') ' LSO data vaild time: ',atime
!
!-----------------------------------------------------------------------
!
!  The winds are read in as dd,ff so we need to convert them to
!  u,v and we'll place the u,v in the slot originally held by
!  the wind gust info ddg,ffg.
!
!  Important detail: At his point winds are "rotated" so that
!  u,v is with respect to the map projection coordinates rather
!  than earth coordinates.   This is different than LAPS and OLAPS
!  where earth coordinate u,v are analysed.
!
!-----------------------------------------------------------------------
!
  DO iob=nxt,(nobs+nobsrd(1))
    obsrd(iob,5,1)=rmiss
    obsrd(iob,6,1)=rmiss
    IF(obsrd(iob,3,1) >= 0. .AND. obsrd(iob,4,1) >= 0.) THEN
      CALL ddrotuv(1,lonsta(iob,1),                                     &
                   obsrd(iob,3,1),obsrd(iob,4,1),                       &
                   ddrot,obsrd(iob,5,1),obsrd(iob,6,1))
    END IF
  END DO
!
!
!-----------------------------------------------------------------------
!
!  Find the x,y location of each station.
!  Note the use of lltoxy requires that map projection
!  be initialized with map parameters.  It is assumed here
!  that this has been done already.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(nobsrd(1),1,latsta(nxt,1),lonsta(nxt,1),                  &
              xsta(nxt),ysta(nxt))

!  write(6,821)
! 821 format(' stn  lat  lon   elev   dd    ff    ugr     vgr ')
!  DO 111 i=1,nobsrd(1)
!    k=i+nobs
!    write(6,822) stn(k,1),latsta(k,1),lonsta(k,1),elevsta(k,1),
!    :          obsrd(k,3,1),obsrd(k,4,1),obsrd(k,5,1),obsrd(k,6,1)
! 822   format(a5,f6.1,f6.1,f6.0,f6.1,f6.1,f6.1,f6.1)
! 111 CONTINUE
!
!
!-----------------------------------------------------------------------
!
!  Temporal QC
!  Read in past data.
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(1x,2a)') 'Getting past surface data from: ',TRIM(tmchkfile)
  CALL read_surface_obs(tmchkfile,blackfile,mxstn,1,atime,n_meso_g,     &
      n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,       &
      n_obs_pos_g,n_obs_b,n_obs_pos_b,stn(1,2),chsrc(1,2),              &
      latsta(1,2),lonsta(1,2),elevsta(1,2),wx(1,2),                     &
      obsrd(1,1,2),obsrd(1,2,2),obsrd(1,3,2),obsrd(1,4,2),              &
      obsrd(1,5,2),obsrd(1,6,2),obsrd(1,7,2),obsrd(1,8,2),              &
      obsrd(1,9,2),  kloud(1,2),obsrd(1,10,2),obsrd(1,11,2),            &
      obsrd(1,12,2),obsrd(1,13,2),idp3(1,2),store_emv(1,1,2),           &
      store_amt(1,1,2),store_hgt(1,1,2),obsrd(1,14,2),                  &
      obstime(1,2),istatus)
  nobsrd(2)=n_obs_b
!
  IF(istatus == 1 .AND. nobsrd(2) > 0 ) THEN   ! old surface OK, do QC
!
!
!-----------------------------------------------------------------------
!
!  The winds are read in as dd,ff so we need to convert them to
!  u,v and we'll place the u,v in the slot originally held by
!  the wind gust info ddg,ffg.
!
!-----------------------------------------------------------------------
!
    DO iob=1,nobsrd(2)
      obsrd(iob,5,2)=rmiss
      obsrd(iob,6,2)=rmiss
      IF(obsrd(iob,3,2) >= 0. .AND. obsrd(iob,4,2) >= 0.)               &
          CALL ddrotuv(1,lonsta(iob,2),                                 &
                     obsrd(iob,3,2),obsrd(iob,4,2),                     &
                     ddrot,obsrd(iob,5,2),obsrd(iob,6,2))
    END DO
!
!    DO 211 i=1,nobsrd(2)
!      write(6,822) stn(i,2),latsta(i,2),lonsta(i,2),elevsta(i,2),
!    :          obsrd(i,3,2),obsrd(i,4,2),obsrd(i,5,2),obsrd(i,6,2)
! 211   CONTINUE
!
  ELSE
    WRITE(6,'(1x,2a,/,3x,a)')  &
    '++Bad status return trying to get surface data from: ',            &
    TRIM(tmchkfile),'  Temporal QC skipped.'
  END IF

!
!-----------------------------------------------------------------------
!
!  Call subroutine to match-up obs and compare hourly changes.
!
!-----------------------------------------------------------------------
!
  CALL qcdata(mxstn,nsrc_sng,nvar_ob,ntime,nobs,nobsrd,stn,chsrc,       &
                    latsta,lonsta,xsta,ysta,elevsta,                    &
                    wx,obsrd,kloud,idp3,store_emv,                      &
                    store_amt,store_hgt,obstime,                        &
                    climin,climax,rely,ival,istatus)
!
!-----------------------------------------------------------------------
!
!  Set isrc based on read-in source.
!
!-----------------------------------------------------------------------
!
  CALL setisrc(mxstn,nsrc_sng,nobsrd(1),name_src,                       &
               chsrc(nxt,1),isrc(nxt))
!
!-----------------------------------------------------------------------
!
!  Create observation array containing observations or
!  derived observations that are to be analysed.
!
!-----------------------------------------------------------------------
!
  DO j=1,nsrc_sng
    PRINT *, ' qsrc = ',(qsrc(ivar,j),ivar=1,nvar_anx)
  END DO
!
  CALL cvtobsanx(mxstn,nvar_ob,nvar_anx,nsrc_sng,nobs,nobsrd(1),        &
                 stn,isrc,latsta,lonsta,elevsta,                        &
                 obsrd,kloud,idp3,store_emv,store_amt,                  &
                 store_hgt,rely,qsrc,rmiss,                             &
                 obs,qobs,qual,istatus)
!
!-----------------------------------------------------------------------
!
!  Call subroutine to do horizontal quality control
!
!-----------------------------------------------------------------------
!
  CALL barqc(mxstn,nvar_anx,nsrc_sng,nobs,nobsrd(1),ipres,iptmp,iqv,    &
             obs,obstime,xsta,ysta,isrc,qual,                           &
             knt,wgtsum,zsum,sqsum,                                     &
             range,wlim,klim,barqclim,var_anx,stn,                      &
             iqclist,iqcsave,qcflag)
!
!-----------------------------------------------------------------------
!
!  Now create superobs by combining observations within a small
!  distance from each other.  This helps eliminate some Barnes
!  artifacts due to the fact that Barnes does not account for the
!  relative distribution of stations that go into the analysis
!  at each grid point.
!
!-----------------------------------------------------------------------
!
  nobs=nobs+nobsrd(1)
  CALL suprob(ntime,mxstn,nvar_anx,nobs,                               &
              stn,latsta(1,1),lonsta(1,1),xsta,ysta,elevsta,           &
              obs,isrc,chsrc,qobs,qual,                                &
              iqspr,sprdist,sprdzlim,qcflag,qcsuper)

  RETURN
END SUBROUTINE prepsfcobs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SUPROB                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
 SUBROUTINE suprob(ntime,mxstn,nvar,nobs,                               &
                   stn,latsta,lonsta,xsta,ysta,elevsta,                 &
                   obs,isrc,chsrc,qobs,qual,                            &
                   iqspr,sprdist,sprdzlim,qcflag,qcsuper)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Combine observations that are within a specified
!  distance (sprdist) from each other.  The combined obs are
!  put into a new "station" and the original obs are set to
!  missing.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  March, 1994
!
!  MODIFICATION HISTORY:
!
!  1/21/96 (K. Brewster)
!  Added full documentation.
!
!  11/01/2011 (K. Brewster)
!  Added calculation of lat,lon for super-ob sites.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'

  INTEGER :: ntime
  INTEGER :: mxstn,nvar,nobs
  CHARACTER (LEN=5) :: stn(mxstn)
  REAL :: latsta(mxstn)
  REAL :: lonsta(mxstn)
  REAL :: xsta(mxstn)
  REAL :: ysta(mxstn)
  REAL :: elevsta(mxstn)
  REAL :: obs(nvar,mxstn)
  INTEGER :: isrc(mxstn)
  CHARACTER (LEN=8) :: chsrc(mxstn,ntime)
  REAL :: qobs(nvar,mxstn)
  INTEGER :: qual(nvar,mxstn)
  INTEGER :: iqspr
  REAL :: sprdist
  REAL :: sprdzlim
  INTEGER :: qcflag
  INTEGER :: qcsuper
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ivar,ista,jsta,jmin,ksta,kntsta
  REAL :: super2,delz,dist1,dist2,dist,distm,qinvi,qinvj,qobloc
!
  CALL zerosrc(nvar,mxstn,obs,qual,isrc,nobs)
!
  super2=sprdist*sprdist
  kntsta=nobs
  ista=0
!
!-----------------------------------------------------------------------
!
!  Loop over all stations
!  Since may be adding stations within the loop we can't
!  use a DO loop with fixed limits.
!
!-----------------------------------------------------------------------
!

  DO
    ista=ista+1
    IF(ista > kntsta) EXIT
    IF(isrc(ista) > 0 .AND. chsrc(ista,1) .ne. 'MDCRS   ') THEN
      distm=1.E20
      jmin=0
!
!-----------------------------------------------------------------------
!
!  Find the nearest neighbor to this station.
!  Only consider neighbors that have height difference less than
!  sprdzlim from this station.  This will allow for stations close
!  in space but on a slope, and also avoids combining sfc data with
!  nearby aircraft data.
!
!-----------------------------------------------------------------------
!
      DO jsta=ista+1,kntsta
        IF(isrc(jsta) > 0 .AND. chsrc(jsta,1) .ne. 'MDCRS   ') THEN
          delz=abs(elevsta(ista)-elevsta(jsta))
          IF(delz < sprdzlim) THEN
            dist1=xsta(ista)-xsta(jsta)
            dist2=ysta(ista)-ysta(jsta)
            dist=dist1*dist1 + dist2*dist2
            IF(dist < distm) THEN
              distm=dist
              jmin=jsta
            END IF
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!
!  If the distance is less than the super-ob threshold
!  distance, combine the neighboring stations into one superob.
!
!-----------------------------------------------------------------------
!
      IF(jmin > 0 .AND. distm < super2) THEN
        WRITE(6,'(a,a,a,a)') ' Combining stations: ',stn(ista),         &
               ' and ',stn(jmin)
        kntsta=kntsta+1
        IF(kntsta > mxstn) THEN
          nobs=mxstn
          EXIT
        END IF
        ksta=kntsta
        isrc(ksta)=isrc(ista)
!
!-----------------------------------------------------------------------
!
!       Note that we use "-1" here instead of "0".  Most of the code
!       checks for a positive number.  For cloud soundings, this presents
!       a problem.  Cloud soundings allow data outside of the analysis domain
!       to be used.  This means we need to flag obs to never be used
!       (SUPROB'd out).
!
!-----------------------------------------------------------------------
!
        isrc(ista)=-1
        isrc(jmin)=-1
        write(stn(ksta),'(a3,a2)') stn(ista)(1:3),stn(jmin)(1:2)
!
!-----------------------------------------------------------------------
!
!  Note here that iqspr designates which variable's q
!  is used to weight the locations and elevation.  General
!  this should be a pressure variable, since the pressure
!  reduction will be affected by the resultant elevation.
!
!-----------------------------------------------------------------------
!
        qinvi=1./qobs(iqspr,ista)
        qinvj=1./qobs(iqspr,jmin)
        qobloc=1./(qinvi + qinvj)
        xsta(ksta)=qobloc*                                              &
                ( qinvi*xsta(ista) + qinvj*xsta(jmin) )
        ysta(ksta)=qobloc*                                              &
                ( qinvi*ysta(ista) + qinvj*ysta(jmin) )
        elevsta(ksta)=qobloc*                                           &
             ( qinvi*elevsta(ista) + qinvj*elevsta(jmin) )
        CALL xytoll(1,1,xsta(ksta),ysta(ksta),latsta(ksta),lonsta(ksta))
!
!-----------------------------------------------------------------------
!
!  Combine each variable, weighting by the inverse of qobs
!
!-----------------------------------------------------------------------
!
        DO ivar=1,nvar
          IF(qual(ivar,ista) > 0 .AND. qual(ivar,jmin) > 0 ) THEN
            qinvi=1./qobs(ivar,ista)
            qinvj=1./qobs(ivar,jmin)
            qobs(ivar,ksta)=1./(qinvi + qinvj)
            obs(ivar,ksta)=qobs(ivar,ksta)*                             &
                          (qinvi*obs(ivar,ista) +                       &
                           qinvj*obs(ivar,jmin) )
            qual(ivar,ista)=qcflag
            qual(ivar,jmin)=qcflag
            qual(ivar,ksta)=qcsuper
          ELSE IF(qual(ivar,jmin) > 0) THEN
            obs(ivar,ksta)=obs(ivar,jmin)
            qobs(ivar,ksta)=qobs(ivar,jmin)
            qual(ivar,ksta)=qual(ivar,jmin)
          ELSE
            obs(ivar,ksta)=obs(ivar,ista)
            qobs(ivar,ksta)=qobs(ivar,ista)
            qual(ivar,ksta)=qual(ivar,ista)
          END IF
        END DO
      END IF
    END IF
  END DO
  nobs=kntsta
  RETURN
END SUBROUTINE suprob
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ZEROSRC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE zerosrc(nvar,mxstn,obs,qual,isrc,nsta)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set isrc to zero if there are no valid data at this
!  obs site.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nvar,mxstn
  REAL :: obs(nvar,mxstn)
  INTEGER :: qual(nvar,mxstn)
  INTEGER :: isrc(mxstn)
  INTEGER :: nsta
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar
  INTEGER num_vars
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO ista=1,nsta
    IF (isrc(ista) == -1) CYCLE
    num_vars = 0
    DO ivar=1,nvar
      IF(qual(ivar,ista) > 0) num_vars = num_vars + 1
    END DO
    IF (num_vars == 0) isrc(ista) = 0
  END DO
  RETURN
END SUBROUTINE zerosrc
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CVTOBSANX                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cvtobsanx(mxstn,nvar_ob,nvar_anx,nsrc,nobs,nobsrd,           &
           stn,isrc,latsta,lonsta,elevsta,                              &
           obsrd,kloud,idp3,store_emv,store_amt,                        &
           store_hgt,rely,qsrc,rmiss,                                   &
           obs,qobs,qual,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert observed surface data to variables to be
!  used in the analysis.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!
!  12/01/98 (K. Brewster)
!  Removed rotation of winds to avoid duplication later.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: mxstn,nvar_ob,nvar_anx,nsrc
  INTEGER :: nobs,nobsrd
  CHARACTER (LEN=5) :: stn(mxstn)
  INTEGER :: isrc(mxstn)
  REAL :: latsta(mxstn)
  REAL :: lonsta(mxstn)
  REAL :: elevsta(mxstn)
  REAL :: obsrd(mxstn,nvar_ob)
  INTEGER :: kloud(mxstn),idp3(mxstn)
  CHARACTER (LEN=1) :: store_emv(mxstn)
  CHARACTER (LEN=4) :: store_amt(mxstn)
  REAL :: store_hgt(mxstn,5)
  INTEGER :: rely(mxstn,nvar_ob)
  REAL :: qsrc(nvar_anx,nsrc)
  REAL :: rmiss
!
  REAL :: obs(nvar_anx,mxstn)
  REAL :: qobs(nvar_anx,mxstn)
  INTEGER :: qual(nvar_anx,mxstn)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Parameters
!
!-----------------------------------------------------------------------
!

  REAL :: ktstoms,mbtopa
  PARAMETER (ktstoms=0.5144,mbtopa=100.)
!
!-----------------------------------------------------------------------
!
!  Conversion functions
!
!-----------------------------------------------------------------------
!
  REAL :: ftoc,alttostpr,msltostpr
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iob,kob,ivar,jsrc,nprev
  REAL :: tf,psta,tk,tdk,qv,qvsat
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'

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
  ftoc(tf)=5.*(tf-32.)/9.
  nprev=nobs
!
  DO ivar=1,nvar_anx
    DO iob=nprev+1,(nprev+nobsrd)
      obs(ivar,iob)=rmiss
      qobs(ivar,iob)=999.
      qual(ivar,iob)=-999
    END DO
  END DO
!
  DO iob=1,nobsrd
    kob=nprev+iob
    IF(isrc(kob) > 0) THEN
      jsrc=isrc(kob)
!
!-----------------------------------------------------------------------
!
!  Winds, u,v kts to m/s
!
!-----------------------------------------------------------------------
!
      IF (obsrd(kob,5) > rmiss .AND. obsrd(kob,6) > rmiss .AND.         &
          rely(kob,5) > 0 .AND. rely(kob,6) > 0 ) THEN
        obs(1,kob)=obsrd(kob,5)*ktstoms
        obs(2,kob)=obsrd(kob,6)*ktstoms
        qobs(1,kob)=qsrc(1,jsrc)
        qobs(2,kob)=qsrc(2,jsrc)
        qual(1,kob)=100
        qual(2,kob)=100
      END IF
!
!-----------------------------------------------------------------------
!
!  Station pressure (mb)
!  For now we'll convert msl pressure to station pressure
!  using the current temperature.
!
!-----------------------------------------------------------------------
!
      IF(obsrd(kob,7) > rmiss .AND. rely(kob,7) > 0) THEN
        psta=mbtopa*obsrd(kob,7)
        obs(3,kob)=psta
        qobs(3,kob)=qsrc(3,jsrc)
        qual(3,kob)=100
      ELSE IF(obsrd(kob,9) > rmiss .AND. rely(kob,9) > 0) THEN
        psta=mbtopa*alttostpr(obsrd(kob,9),elevsta(kob))
        obs(3,kob)=psta
        qobs(3,kob)=qsrc(3,jsrc)
        qual(3,kob)=100
      ELSE IF(obsrd(kob,8) > rmiss .AND.                                &
                obsrd(kob,1) > rmiss .AND.                              &
                rely(kob,8) > 0 .AND. rely(kob,1) > 0 .AND.             &
                elevsta(kob) < 500. ) THEN
        tk=ftoc(obsrd(kob,1)) + 273.15
        psta=mbtopa*msltostpr(obsrd(kob,8),tk,elevsta(kob))
!        print *, ' msl,elev,stpr: ',
!    :                 obsrd(kob,8),elevsta(kob),(0.01*psta)
        obs(3,kob)=psta
        qobs(3,kob)=qsrc(3,jsrc)
        qual(3,kob)=100
      ELSE IF(obsrd(kob,8) > rmiss .AND. rely(kob,8) > 0 .AND.          &
                elevsta(kob) == 0. ) THEN
        psta=mbtopa*obsrd(kob,8)
        obs(3,kob)=psta
        qobs(3,kob)=qsrc(3,jsrc)
        qual(3,kob)=100
      ELSE
        psta=mbtopa*alttostpr(1013.,elevsta(kob))
      END IF
!      print *, ' stn,pr,qobs,qual: ',stn(kob),obs(3,kob),
!    :                                 qobs(3,kob),qual(3,kob)
!
!-----------------------------------------------------------------------
!
!  Potential temperature (K)
!
!-----------------------------------------------------------------------
!
      IF(obsrd(kob,1) > rmiss .AND. rely(kob,1) > 0) THEN
        tk=ftoc(obsrd(kob,1)) + 273.15
        obs(4,kob)=tk*((p0/psta)**rddcp)
        qobs(4,kob)=qsrc(4,jsrc)
        qual(4,kob)=100
      END IF
!
!-----------------------------------------------------------------------
!
!  Specific humidity (kg/kg)
!
!-----------------------------------------------------------------------
!
      IF(obsrd(kob,1) > rmiss .AND. obsrd(kob,2) > rmiss .AND.          &
            rely(kob,1) > 0 .AND. rely(kob,2) > 0 ) THEN
        tk=ftoc(obsrd(kob,1)) + 273.15
        tdk=ftoc(obsrd(kob,2)) + 273.15
        qv = f_qvsat( psta,tdk )
        qvsat = f_qvsat( psta, tk )
        obs(5,kob)=MAX(1.0E-08,qv)
        qobs(5,kob)=qsrc(5,jsrc)
        qual(5,kob)=100
      END IF
    END IF
  END DO
  RETURN
END SUBROUTINE cvtobsanx
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SETISRC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setisrc(mxstn,nsrc,nobs,name_src,                            &
           chsrc,isrc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign source number to observations based on read-in
!  source name.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  Aug, 1997
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: mxstn,nsrc
  INTEGER :: nobs
  CHARACTER (LEN=8) :: name_src(nsrc)
  CHARACTER (LEN=8) :: chsrc(mxstn)
  INTEGER :: isrc(mxstn)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iob,jsrc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO iob=1,nobs
    DO jsrc=1,nsrc
      IF(chsrc(iob) == name_src(jsrc)) GO TO 120
    END DO
    WRITE(6,'(a,a)') chsrc(iob),' source not matched '
    isrc(iob)=0
    CYCLE
    120   CONTINUE
    isrc(iob)=jsrc
!    write(6,'(a,a,i4)') chsrc(iob),  &
!                         ' source matched isrc=',isrc(iob)
  END DO
!
  RETURN
END SUBROUTINE setisrc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE QCDATA                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qcdata(mxstn,nsrc,nvar_ob,ntime,nobs,nobsrd,stn,chsrc,       &
           latsta,lonsta,xsta,ysta,elevsta,                             &
           wx,obsrd,kloud,idp3,store_emv,                               &
           store_amt,store_hgt,obstime,                                 &
           climin,climax,rely,ival,istatus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!   Perform initial quality control checks including a gross range
!   check and a temporal consistency check.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Adapted from FSL/LAPS by Keith Brewster, CAPS
!
!  MODIFICATIONS:
!    2004-09-09 Keith Brewster, based on similar mod by Eric Kemp, TASC
!    Changed logic to always perform gross limits check, previously was
!    skipping that if previous hour's data were missing.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: mxstn,nsrc,nvar_ob,ntime
!
!-----------------------------------------------------------------------
!
!  Arrays for the input obs at two different times
!
!-----------------------------------------------------------------------
!
  INTEGER :: nobs,nobsrd(ntime)
  REAL :: latsta(mxstn,ntime),lonsta(mxstn,ntime)
  REAL :: elevsta(mxstn,ntime)
  REAL :: xsta(mxstn),ysta(mxstn)
  INTEGER :: obstime(mxstn,ntime)
  CHARACTER (LEN=5) :: stn(mxstn,ntime)
  CHARACTER (LEN=8) :: chsrc(mxstn,ntime)
  CHARACTER (LEN=8) :: wx(mxstn,ntime)
  CHARACTER (LEN=1) :: store_emv(mxstn,5,ntime)
  CHARACTER (LEN=4) :: store_amt(mxstn,5,ntime)
  INTEGER :: kloud(mxstn,ntime),idp3(mxstn,ntime)
  REAL :: store_hgt(mxstn,5,ntime)
  REAL :: obsrd(mxstn,nvar_ob,ntime)
  REAL :: climin(nvar_ob),climax(nvar_ob)
  INTEGER :: rely(mxstn,nvar_ob,ntime)
  INTEGER :: ival(mxstn)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: missing
  INTEGER :: ktime,jvar,iob,kob
  INTEGER :: imissing,istatus,jstatus
!
  PARAMETER ( missing = -99.9,                                          &
             imissing = -99)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0
  !WRITE(6,'(a)') ' ** begin temporal qc on surface data '
  DO ktime=1,ntime
    DO jvar=1,nvar_ob
      DO iob=1,mxstn
        rely(iob,jvar,ktime) = imissing
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Climatological extreme test
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a,i5,i6)') ' # of data stns available: ',                   &
                       nobsrd(1),nobsrd(2)
!
  DO jvar = 1,nvar_ob
    DO iob = 1,nobsrd(1)

      kob=nobs+iob
      IF(obsrd(kob,jvar,1) >= climin(jvar) .AND.                        &
          obsrd(kob,jvar,1) <= climax(jvar) )                           &
          rely(kob,jvar,1)=10

    END DO
  END DO

  DO jvar = 1,nvar_ob
    DO iob = 1,nobsrd(2)

      IF(obsrd(iob,jvar,2) >= climin(jvar) .AND.                        &
          obsrd(iob,jvar,2) <= climax(jvar) )                           &
          rely(iob,jvar,2)=10

    END DO
  END DO

  IF(nobsrd(2) > 0) &
     CALL sta_match(mxstn,ntime,nobs,nobsrd,stn,chsrc,ival)
!
!-----------------------------------------------------------------------
!
!  Standard deviation check
!
!-----------------------------------------------------------------------
!
  DO jvar=1,nvar_ob
    CALL dev_ck(mxstn,nvar_ob,ntime,jvar,nobs,nobsrd,                   &
                    obsrd,rely,ival)
  END DO
!
!-----------------------------------------------------------------------
!
!  Temporal QC finished
!
!-----------------------------------------------------------------------
!
  jstatus = 1
  WRITE(6,'(a)')' ** temporal QC of surface data complete. ** '
  RETURN
END SUBROUTINE qcdata
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STA_MATCH                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sta_match(mxstn,ntime,nobs,nobsrd,stn,chsrc,ival)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Match-up stations by their name and source name.
!  This provides a mapping of the data from one time to another.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: mxstn,ntime
  INTEGER :: nobs
  INTEGER :: nobsrd(ntime)
  CHARACTER (LEN=5) :: stn(mxstn,ntime)
  CHARACTER (LEN=8) :: chsrc(mxstn,ntime)
  INTEGER :: ival(mxstn)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: m,n,kntno,kob
!
  kntno=0
  DO n=1,nobsrd(1)
    kob=n+nobs
    ival(kob) = -99
    DO m=1,nobsrd(2)
      IF ( stn(m,2) == stn(kob,1) ) THEN
        IF( chsrc(m,2) == chsrc(kob,1) ) THEN
          ival(kob) = m
          GO TO 20
        END IF
      END IF
    END DO
    kntno=kntno+1
!    print *, ' Cannot find match for ',stn(kob,1)
20 CONTINUE
  END DO
  WRITE(6,'(a,i4,a,/a,i4,a)')                                           &
          ' FYI from sta_match: ',kntno,' stations ',                   &
          ' of ',nobsrd(1),' were not matched.'
  RETURN
END SUBROUTINE sta_match
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DEV_CHK                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dev_ck(mxstn,nvar_ob,ntime,ivar,nobs,nobsrd,                 &
           obsrd,rely,ival)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check hourly changes for an individual observed
!  variable.  A mean hourly change is computed for the whole
!  field and the change at an individual station is compared
!  to the rms change for the group.
!
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Original version: FSL LAPS team
!  Restructured to remove local storage arrays
!  Keith Brewster, CAPS, March, 1994.
!
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: mxstn,nvar_ob,ntime
  INTEGER :: ivar
  INTEGER :: nobs,nobsrd(ntime)
  REAL :: obsrd(mxstn,nvar_ob,ntime)
  INTEGER :: rely(mxstn,nvar_ob,ntime)
  INTEGER :: ival(mxstn)
!
  INTEGER :: m,n,numdiff,kob
  REAL :: missing,diff,sumdiff,sumsq,std,avgdiff
  REAL :: diffnor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  missing = -99.
!
!-----------------------------------------------------------------------
!
!  Compute avg change from prev hr (omit mesonet pressures)
!
!-----------------------------------------------------------------------
!
  std=0.
  avgdiff=0.
  sumdiff=0.
  sumsq=0.
  numdiff = 0
  DO n=1,nobsrd(1)
    kob=n+nobs
    m  = ival(kob)
    IF(m > 0) THEN
      IF(obsrd(kob,ivar,1) > missing .AND. obsrd(m,ivar,2) > missing .AND. &
            rely(m,ivar,2) > 0) THEN
        diff = obsrd(kob,ivar,1) - obsrd(m,ivar,2)
        numdiff = numdiff + 1
        sumdiff = sumdiff + diff
        sumsq=sumsq+diff*diff
      END IF
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!    Compute average and standard deviation of hourly change.
!
!-----------------------------------------------------------------------
!
  IF (numdiff < 2) THEN
    WRITE(6,'(a)') ' numdiff < 2, no dev_ck done.'
    RETURN
  ELSE
    avgdiff = sumdiff/FLOAT(numdiff)
    sumsq = sumsq - sumdiff*avgdiff
    IF(sumsq > 0. .AND. numdiff > 1) std=SQRT(sumsq/FLOAT((numdiff-1)))
  END IF
!
  IF (std == 0.) THEN
    WRITE(6,'(a)') ' std dev = 0., dev_ck not completed...'
    RETURN
  ELSE
!
!-----------------------------------------------------------------------
!
!  Add or subtract reliability points based on the
!  hourly difference normalized by the std.
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(a,i4,a,f10.2,a,f10.2)') ' Hourly change: ivar=',ivar,     &
             ' mean= ',avgdiff,' stddev= ',std
    DO n=1,nobsrd(1)
      kob=n+nobs
      m  = ival(kob)
      IF(m > 0) THEN
        IF(obsrd(kob,ivar,1) > missing .AND. obsrd(m,ivar,2) > missing .AND. &
              rely(m,ivar,2) > 0) THEN
          diff=ABS(obsrd(kob,ivar,1)-obsrd(m,ivar,2)-avgdiff)
          diffnor=diff/std
          IF(diffnor > 4.) THEN
            rely(kob,ivar,1)=rely(kob,ivar,1)-25
          ELSE
            rely(kob,ivar,1)=rely(kob,ivar,1)+25
          END IF
        END IF
      END IF
    END DO
    !WRITE(6,'(a,i4)') ' subr dev_ck complete ', ivar
  END IF
  RETURN
END SUBROUTINE dev_ck
