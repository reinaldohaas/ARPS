!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE get_88dinfo                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
SUBROUTINE get_88dinfo(fname,radname,ivcp,itimfrst,                    &
                       radlat,radlon,radelv, istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Initialize radar file and obtain information about an 88d radar site.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2009/10/27
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: fname
  CHARACTER(LEN=4), INTENT(IN) :: radname

  INTEGER, INTENT(OUT) :: ivcp
  INTEGER, INTENT(OUT) :: itimfrst
  REAL, INTENT(OUT) :: radlat
  REAL, INTENT(OUT) :: radlon
  REAL, INTENT(OUT) :: radelv
  INTEGER, INTENT(OUT) :: istatus

  INTEGER :: irc
  INTEGER :: ialt,ilat,ilon
  INTEGER :: iyr,iyear,imon,iday,ihour,imin,isec
!
!-----------------------------------------------------------------------
!
! Radar File Control Common
!
!-----------------------------------------------------------------------
!
  LOGICAL :: read_next
  COMMON /rd88dcntl/ read_next

  istatus=0
  irc=-9
  read_next=.TRUE.

  WRITE(6,'(1x,2a)') 'RadarName = ',TRIM(radname)

  CALL set_radar_name_f(radname,4)

  WRITE(6,'(1x,2a)') 'fname = ',TRIM(fname)
  CALL radar_open_f(fname, irc )

  CALL get_radar_info_f()

  call get_latitude_f(ilat)
  call get_longitude_f(ilon)
  call get_altitude_f(ialt)

  radlat=  0.00001 * ilat
  radlon= -0.00001 * ilon
  radelv= float(ialt)

  WRITE(6,'(a,f10.4,a,f10.4,a,f8.1)') &
    ' Radar Lat = ',radlat,'  Lon = ',radlon,'  Height = ',radelv

  IF(irc == 0) THEN

    CALL read_radial_f(istatus)
    WRITE(6,'(a,i6)') ' Read_radial returned ',istatus

    IF( istatus == 1 ) THEN

      WRITE(6,'(a)') ' Read_radial returned double eof'

    ELSE

      read_next=.FALSE.
      CALL get_vcp_f(ivcp)
      WRITE(6,'(a,i6)') ' Found vcp as: ',ivcp

      CALL get_year_f(iyr)
      CALL get_month_f(imon)
      CALL get_day_f(iday)
      CALL get_hour_f(ihour)
      CALL get_min_f(imin)
      CALL get_sec_f(isec)

      iyear = iyr+1900
      if ( iyear < 1960 ) iyear = iyear + 100
      WRITE(6,'(a,6i5)') ' Initial Time: ',iyear,imon,iday,ihour,imin,isec

      CALL ctim2abss(iyear,imon,iday,ihour,imin,isec,itimfrst)
      WRITE(6,'(a,i16)') '  itimfrst: ',itimfrst

    END IF

  ELSE

    istatus = irc

  END IF

  RETURN
END SUBROUTINE get_88dinfo
!
!########################################################################
!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE get_88draddims              #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE get_88draddims(ivcp,maxrgate,maxvgate,maxazim,maxelev,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set data dimension variables for 88D NEXRAD radars.
!  Keys off of the VCP number.  See SELECT CASE in this subroutine.
!  Will need to be updated for any new VCP's that are used.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2009/07/24
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: ivcp
  INTEGER, INTENT(OUT) :: maxrgate
  INTEGER, INTENT(OUT) :: maxvgate
  INTEGER, INTENT(OUT) :: maxazim
  INTEGER, INTENT(OUT) :: maxelev
  INTEGER, INTENT(OUT) :: istatus

  REAL, PARAMETER :: max_elev = 24

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  maxrgate = 1840
  maxvgate = 920
  maxazim  = 800

  SELECT CASE (ivcp)
    CASE (11)
      maxelev = 16
    CASE (12)
      maxelev = 17
    CASE (21)
      maxelev = 11
    CASE (31)
      maxelev = 7
    CASE (32)
      maxelev = 7
    CASE (121)
      maxelev = 20
    CASE (211)
      maxelev = 16
    CASE (212)
      maxelev = 17
    CASE (221)
      maxelev = 11
    CASE DEFAULT
      maxelev = max_elev
  END SELECT

  istatus = 0

  RETURN
END SUBROUTINE get_88draddims
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rdtilt88d                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rdtilt88d(maxrgate,maxvgate,maxazim,dualpproc,fname,          &
                     refvarname,rhvvarname,zdrvarname,phivarname,        &
                     velvarname,spwvarname,itimfrst,                     &
                     ng_ref,rfrst_ref,gtspc_ref,                         &
                     ng_vel,rfrst_vel,gtspc_vel,                         &
                     nazim,azim,elev,time,vnyq,                          &
                     refl,rhv,zdr,phi,vel,spw,                           &
                     refstat,rhvstat,dlpstat,velstat,eofstat,istatus)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Uses a2io calls to read 1-tilt of radar data from an 88d data file.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster 24 July 2009
!
! MODIFICATIONS:
! 19 April 2012 Keith Brewster
! Added variable names to argument list to allow flexibility in
! reading files with folded and unfolded velocities, for example.
!
! 07 Sept 2012 Keith Brewster
! Added processing of dual-pol variables.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: maxrgate
  INTEGER, INTENT(IN) :: maxvgate
  INTEGER, INTENT(IN) :: maxazim
  LOGICAL, INTENT(IN) :: dualpproc
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=80), INTENT(IN) :: refvarname
  CHARACTER (LEN=80), INTENT(IN) :: rhvvarname
  CHARACTER (LEN=80), INTENT(IN) :: zdrvarname
  CHARACTER (LEN=80), INTENT(IN) :: phivarname
  CHARACTER (LEN=80), INTENT(IN) :: velvarname
  CHARACTER (LEN=80), INTENT(IN) :: spwvarname
  INTEGER, INTENT(IN) :: itimfrst
  INTEGER, INTENT(OUT) :: ng_ref
  INTEGER, INTENT(OUT) :: rfrst_ref
  INTEGER, INTENT(OUT) :: gtspc_ref
  INTEGER, INTENT(OUT) :: ng_vel
  INTEGER, INTENT(OUT) :: rfrst_vel
  INTEGER, INTENT(OUT) :: gtspc_vel
  INTEGER, INTENT(OUT) :: nazim
!
  REAL, INTENT(OUT) :: azim(maxazim)
  REAL, INTENT(OUT) :: elev(maxazim)
  INTEGER, INTENT(OUT) :: time(maxazim)
  REAL, INTENT(OUT) :: vnyq(maxazim)
  REAL, INTENT(OUT) :: refl(maxrgate,maxazim)
  REAL, INTENT(OUT) :: rhv(maxrgate,maxazim)
  REAL, INTENT(OUT) :: zdr(maxrgate,maxazim)
  REAL, INTENT(OUT) :: phi(maxrgate,maxazim)
  REAL, INTENT(OUT) :: vel(maxvgate,maxazim)
  REAL, INTENT(OUT) :: spw(maxvgate,maxazim)
  INTEGER, INTENT(OUT) :: refstat
  INTEGER, INTENT(OUT) :: rhvstat
  INTEGER, INTENT(OUT) :: dlpstat
  INTEGER, INTENT(OUT) :: velstat
  INTEGER, INTENT(OUT) :: eofstat
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! PARAMETERS
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: GOOD_STATUS=1
  INTEGER, PARAMETER :: max_bad_stat=100
!
!-----------------------------------------------------------------------
!
! Radar File Control Common
!
!-----------------------------------------------------------------------
!
  LOGICAL :: read_next
  COMMON /rd88dcntl/ read_next

!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ref_index,rhv_index,zdr_index,phi_index,vel_index,spw_index
  INTEGER :: iazim,jazim
  INTEGER :: iangle,iscan,itilt,azitime
  INTEGER :: past_scan,past_tilt,past_angle
  INTEGER :: knt_bad_stat
  INTEGER :: iyr,iyear,imon,iday,ihour,imin,isec
  INTEGER :: refstatin,rhvstatin,dlpstatin,velstatin,spwstatin
  INTEGER :: spwstat
  INTEGER :: iread
  INTEGER :: vtmp
  REAL :: rmax

  LOGICAL :: restart,initial_ray
  INTEGER :: iostatus

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ng_ref = 0
  rfrst_ref = -999
  gtspc_ref = -999
  ng_vel = 0
  rfrst_vel = -999
  gtspc_vel = -999
  azim = -999.
  elev = -999.
  time = -999
  vnyq = -999.
  refl = -999.
  vel = -999.
  spw = -999.

  ref_index=-999
  rhv_index=-999
  zdr_index=-999
  phi_index=-999
  vel_index=-999
  spw_index=-999

  nazim = 0
  azitime = 0
  istatus = -1
  refstat = -9
  rhvstat = -9
  dlpstat = -9
  velstat = -9
  spwstat = -9
  eofstat = 0
  knt_bad_stat = 0

  WRITE(6,'(1x,2a)') 'Reading a tilt of 88d radar data from: ', TRIM(fname)

  CALL get_field_num_f(TRIM(refvarname),ref_index)
  IF(ref_index < 0) CALL get_field_num_f('DBZ',ref_index)
  WRITE(6,'(a,i4)') ' Retrieved reflectivity index as ',ref_index

  CALL get_field_num_f(TRIM(rhvvarname),rhv_index)
  IF(rhv_index < 0) CALL get_field_num_f('RHO',rhv_index)
  WRITE(6,'(a,i4)') ' Retrieved Rho-HV index as ',rhv_index

  IF(dualpproc) THEN
    CALL get_field_num_f(TRIM(zdrvarname),zdr_index)
    IF(zdr_index < 0) CALL get_field_num_f('ZDR',zdr_index)
    WRITE(6,'(a,i4)') ' Retrieved Zdr index as ',zdr_index

    CALL get_field_num_f(TRIM(phivarname),phi_index)
    IF(phi_index < 0) CALL get_field_num_f('PHI',phi_index)
    WRITE(6,'(a,i4)') ' Retrieved Phi index as ',phi_index
  END IF

  CALL get_field_num_f(TRIM(velvarname),vel_index)
  IF(vel_index < 0) CALL get_field_num_f('VEL',vel_index)
  WRITE(6,'(a,i4)') ' Retrieved velocity index as ',vel_index

  CALL get_field_num_f(TRIM(spwvarname),spw_index)
  IF(spw_index < 0) CALL get_field_num_f('SPW',spw_index)
  WRITE(6,'(a,i4)') ' Retrieved spectrum width index as ',spw_index
!
! Loop for this tilt
!
  jazim=0
  initial_ray=.TRUE.
  DO iread=1,10000

    refstatin = -1
    rhvstatin = -1
    dlpstatin = -1
    velstatin = -1
!
!-----------------------------------------------------------------------
!
!  Read radial of data
!
!-----------------------------------------------------------------------
!
    IF(read_next) THEN

      CALL read_radial_f(iostatus)

!     WRITE(6,'(a,i6)') ' Read_radial returned:',iostatus
      IF( iostatus == 1 ) THEN
        WRITE(6,'(a)') ' Read_radial returned double eof'
        eofstat = 1
        EXIT
      END IF

    ELSE

      read_next = .TRUE.

    END IF

    CALL get_restart_flag_f(restart)

    IF(restart) THEN
      WRITE(6,'(a)') ' Restart detected, returning to initial azim'
      jazim=0
      initial_ray=.TRUE.
    END IF

    CALL get_status_f(ref_index, refstatin)
    CALL get_status_f(vel_index, velstatin)
    CALL get_status_f(spw_index, spwstatin)
    IF( rhv_index > -1 ) CALL get_status_f(rhv_index, rhvstatin)
    IF( zdr_index > -1 ) CALL get_status_f(zdr_index, dlpstatin)

    IF ( refstatin == GOOD_STATUS .OR. velstatin == GOOD_STATUS ) THEN
      jazim=jazim+1
      knt_bad_stat = 0
      CALL get_scan_f(iscan)
      CALL get_tilt_f(itilt)
      CALL get_fixed_angle_f(iangle)

      IF(initial_ray) THEN

        past_scan  = iscan
        past_tilt  = itilt
        past_angle = iangle

        CALL get_number_of_gates_f(ref_index,ng_ref)
        CALL get_first_gate_f(ref_index,rfrst_ref)
        CALL get_gate_spacing_f(ref_index,gtspc_ref)

        CALL get_number_of_gates_f(vel_index,ng_vel)
        CALL get_first_gate_f(vel_index,rfrst_vel)
        CALL get_gate_spacing_f(vel_index,gtspc_vel)

        rmax = 0
        IF(ng_ref > 0) THEN
          ng_ref=min(ng_ref,maxrgate)
          rmax = FLOAT(rfrst_ref+(ng_ref*gtspc_ref))
        END IF
        IF(ng_vel > 0) THEN
          ng_vel=min(ng_vel,maxvgate)
          rmax = max(rmax,FLOAT(rfrst_vel+(ng_vel*gtspc_vel)))
        END IF

        WRITE(6,'(a,i8,a,i8)')  &
          ' rfrst_ref: ',rfrst_ref,'  rfrst_vel: ',rfrst_vel
        WRITE(6,'(a,i8,a,i8)')  &
          ' gtspc_ref: ',gtspc_ref,'  gtspc_vel: ',gtspc_vel
        WRITE(6,'(a,i8,a,i8)')  &
          ' ng_ref: ',ng_ref,'  ng_vel: ',ng_vel
        WRITE(6,'(a,f10.1,a)')  ' rmax: ',(0.001*rmax),' km'

        initial_ray=.FALSE.

      ELSE

        IF( itilt /= past_tilt ) THEN
          WRITE(6,'(a,i8,a,i8)') ' itilt changed from ',past_tilt,' to ',itilt
          read_next=.FALSE.
          jazim=jazim-1
          EXIT
       END IF

      END IF  ! initial ray

      refstat=refstatin
      rhvstat=rhvstatin
      dlpstat=dlpstatin
      velstat=velstatin
      spwstat=spwstatin
!
!     If we will break an array, just recycle the last slot in the array
!     instead of risking a SIGSEGV.
!
!     This scenerio should only be tripped if there is a malfunction with the
!     88D or we've received corrupted data.  The default array size is larger
!     that what should be received.
!

      IF(jazim > maxazim) THEN
        jazim = maxazim
      END IF

      IF(refstat == GOOD_STATUS) THEN
        CALL get_data_field_f(ref_index,refl(1,jazim),maxrgate,iostatus)
        IF( rhvstat == GOOD_STATUS)                                    &
          CALL get_data_field_f(rhv_index,rhv(1,jazim),maxrgate,iostatus)
        IF( dlpstat == GOOD_STATUS) THEN
          CALL get_data_field_f(zdr_index,zdr(1,jazim),maxrgate,iostatus)
          CALL get_data_field_f(phi_index,phi(1,jazim),maxrgate,iostatus)
        END IF
      END IF
      IF( velstat == GOOD_STATUS) THEN
        CALL get_data_field_f(vel_index, vel(1,jazim),maxvgate,iostatus)
        IF( spwstat == GOOD_STATUS ) &
          CALL get_data_field_f(spw_index, spw(1,jazim),maxvgate,iostatus)
        CALL get_nyquist_f(vtmp)
        vnyq(jazim) = 0.01 * vtmp
!       IF(mod(iread,20) == 0) print *, ' nyquist vel =',vnyq(jazim)
      END IF

      elev(jazim)=0.01*iangle

      CALL get_azi_f(iazim)
      azim(jazim)=0.01*iazim
!     IF(mod(iread,20) == 0) print *, ' iread:',iread,'  azim:',azim(jazim)

      CALL get_year_f(iyr)
      CALL get_month_f(imon)
      CALL get_day_f(iday)
      CALL get_hour_f(ihour)
      CALL get_min_f(imin)
      CALL get_sec_f(isec)

      iyear = iyr+1900
      if ( iyear < 1960 ) iyear = iyear + 100

      CALL ctim2abss(iyear,imon,iday,ihour,imin,isec,azitime)
      time(jazim)=azitime-itimfrst

!     IF( mod(jazim,20) == 0 ) THEN
!       print *, ' refstat=',refstat
!       print *, ' sample reflectivities: ',refl(50,jazim),refl(100,jazim)
!       print *, ' velstat=',velstat
!       print *, ' sample velocity: ',vel(50,jazim),vel(100,jazim)
!     END IF

    ELSE  ! increment bad status counter

      knt_bad_stat=knt_bad_stat+1
      IF( knt_bad_stat > max_bad_stat ) THEN
        jazim=jazim-1
        EXIT
      END IF

    END IF  ! ref and vel status check


  END DO ! iread

  nazim=max(jazim,0)
  print *, ' elev= ',elev(nazim),'   nazim = ',nazim
  istatus = 0

  RETURN

END SUBROUTINE rdtilt88d
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  SUBROUTINE phi2kdp                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
SUBROUTINE phi2kdp(maxgate,maxazim,ngate,nazim,                        &
                   refrange,refl,rhv,phi,kdp,rtem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Convert from observed angle offset, phi, to local derivative of phi, Kdp.
!  Phi data are noisy so smoothing is applied and a least squares line is 
!  used to determine slope.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  2012/09/12
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: nazim
  REAL, INTENT(IN) :: refrange(maxgate)        ! km
  REAL, INTENT(IN) :: refl(maxgate,maxazim)    ! Reflectivity dBZ
  REAL, INTENT(IN) :: rhv(maxgate,maxazim)     ! Rho-HV
  REAL, INTENT(INOUT) :: phi(maxgate,maxazim)  ! degrees
  REAL, INTENT(OUT) :: kdp(maxgate,maxazim)    ! degrees per km
  REAL, INTENT(OUT) :: rtem(maxgate,maxazim)   ! temporary array
!
! Parameters
!
  REAL, PARAMETER :: phichek = -400.
  REAL, PARAMETER :: rhvthresh = 0.80
  REAL, PARAMETER :: trwthresh = 40.0     ! dBZ, TRW reflectivity
  REAL, PARAMETER :: rmiss = -999.
  INTEGER, PARAMETER :: ngate_conv = 4
  INTEGER, PARAMETER :: ngate_strat = 12
!
! Misc Local Variables
!
  INTEGER :: igate,jazim,kgate,idiag
  INTEGER :: ibgn,iend,nlsqgate
  INTEGER :: knt,istatus
  REAL :: avgsum,const,slope
  REAL :: avgpast,offset,phiadj
!
! Array initializations
!
  kdp=rmiss
  rtem=rmiss
!
!   Filter-out low correlation coefficient data
!
  DO jazim=1,nazim
    idiag=0
    DO igate=1,ngate
      IF(phi(igate,jazim) > phichek) THEN
        IF(rhv(igate,jazim) < rhvthresh) phi(igate,jazim)=rmiss
      END IF
    END DO
    IF(mod(jazim,30) == 0) THEN
!     idiag=10
      print *, ' pht2kdp jazim =',jazim
      WRITE(6,'(20f6.1)') (phi(igate,jazim),igate=101,120)
    END IF
!
!   Smooth Phi data using local running average and remove speckles
!   Also account for phi growing greater than 360 degrees, which then
!   is aliased to low phi.
!   Length of local averaging is a function of reflectivity
!
    DO igate=1,ngate
      avgpast=0.
      offset=0.
      IF(phi(igate,jazim) > phichek) THEN
        IF(refl(igate,jazim) < trwthresh) THEN
          ibgn=max(1,(igate-ngate_strat))
          iend=min(ngate,(igate+ngate_strat))
        ELSE
          ibgn=max(1,(igate-ngate_conv))
          iend=min(ngate,(igate+ngate_conv))
        END IF
        knt=0
        avgsum=0.
        DO kgate=ibgn,iend
          IF(phi(kgate,jazim) > phichek) THEN
            knt=knt+1
            phiadj=phi(kgate,jazim)+offset
            IF((phiadj - avgpast) < -180.0) phiadj=phiadj + 360.0
            avgsum=avgsum+phiadj
          END IF
        END DO
        IF( knt > 1 ) THEN
          rtem(igate,jazim)=avgsum/float(knt)
          IF((rtem(igate,jazim) - offset) > 360.) offset = offset + 360.
          avgpast=rtem(igate,jazim)
        END IF
      END IF
    END DO
    IF(mod(jazim,30) == 0) THEN
      print *, ' smoothed jazim =',jazim
      WRITE(6,'(20f6.1)') (rtem(igate,jazim),igate=101,120)
    END IF
!
!   Find local least squares fit to a line.
!   Length of local least-squares fit is a function of reflectivity
!
    DO igate=1,ngate
      IF(rtem(igate,jazim) > phichek) THEN
        IF(refl(igate,jazim) < trwthresh) THEN
          ibgn=max(1,(igate-ngate_strat))
          iend=min(ngate,(igate+ngate_strat))
        ELSE
          ibgn=max(1,(igate-ngate_conv))
          iend=min(ngate,(igate+ngate_conv))
        END IF
        nlsqgate=(iend-ibgn)+1
        IF(nlsqgate > 4) THEN
          CALL lsqline(nlsqgate,phichek,idiag,                       &
                       rtem(ibgn,jazim),refrange(ibgn),              &
                       const,slope,istatus) 
          IF(istatus == 1) kdp(igate,jazim)=slope
        END IF
      END IF
    END DO
    IF(mod(jazim,30) == 0) THEN
      print *, ' kdp jazim =',jazim
      WRITE(6,'(20f6.1)') (kdp(igate,jazim),igate=101,120)
    END IF
 
  END DO  ! azim loop
  RETURN
END SUBROUTINE phi2kdp

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  SUBROUTINE lsqline                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
SUBROUTINE lsqline(npts,ychek,idiag,y,x,const,slope,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Least squares solution for y = const + slope*x
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  2012/09/12
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: npts
  REAL, INTENT(IN)     :: ychek
  INTEGER, INTENT(IN)  :: idiag
  REAL, INTENT(IN)     :: y(npts)
  REAL, INTENT(IN)     :: x(npts)
  REAL, INTENT(OUT)    :: const
  REAL, INTENT(OUT)    :: slope
  INTEGER, INTENT(OUT) :: istatus
!
! LSQ matrix and solution vector
!
  REAL, PARAMETER :: eps=1.0E-20

  REAL :: array(2,2)
  REAL :: rhs(2)
  REAL :: sol(2)
  REAL :: work(2,3)
  REAL :: work1d(3)
!
! Misc local yiables
!
  INTEGER :: i,knt
!
! Initialization
!
  istatus=0
  array=0.0
  rhs=0.0
  sol=0.0
  knt=0
!
! Fill LSQ matrix and RHS vector
!
  IF(idiag > 0) print *, ' npts = ',npts
  DO i=1,npts
    IF(y(i) > ychek) THEN
      knt=knt+1
      array(1,1)=array(1,1)+1.0
      array(1,2)=array(1,2)+x(i)
      array(2,1)=array(2,1)+x(i)
      array(2,2)=array(2,2)+x(i)*x(i)
      rhs(1)=rhs(1)+y(i)
      rhs(2)=rhs(2)+y(i)*x(i)
    END IF
  END DO

  CALL gjelim(2,array,rhs,sol,work,work1d,eps,istatus)

  IF( idiag > 0) print *, ' knt= ',knt,'  istatus= ',istatus

  IF( knt > 4 ) THEN

    IF(istatus == 1) THEN
      const=sol(1)
      slope=sol(2)
      IF(idiag > 0) THEN
        print *, ' const = ',const,'  slope = ',slope
      END IF
    END IF
  ELSE
    istatus=-9
  END IF

  RETURN
END SUBROUTINE lsqline
