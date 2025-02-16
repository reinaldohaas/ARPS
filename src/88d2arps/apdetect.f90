!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE APDETECT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE apdetect(nrefgates,nvelgates,maxazim,maxelev,                &
                    kntrgate,kntrazim,kntrelev,                         &
                    kntvgate,kntvazim,kntvelev,                         &
                    refchek,velchek,                                    &
                    irefgatsp,ivelgatsp,                                &
                    winszrad,winszazim,i_vcp,gcopt,gcvrlim,             &
                    rngrvol,azmrvol,elvrvol,                            &
                    rngvvol,azmvvol,elvvvol,                            &
                    refvol,velvol,tmp,                                  &
                    istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Anomalous propagation check and clutter check for radar data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: LeiLei Wang, CAPS
!
!  MODIFICATION HISTORY:
!   Keith Brewster,
!   Modified to remove additional array allocation
!
!   07/08/2003 Keith Brewster
!   Added setting of istatus return values.
!
!   02/26/2009 Keith Brewster
!   Separate ground clutter detection options and code added.
!
!   03/20/2011 Keith Brewster
!   Adjusted thresholding to better identify wind farm clutter.
!   Cleaned-up clutter code addded by Zhao Kun placing it under gcopt=2.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    dazlim   Maximum value of azimuth difference (grid vs data) to accept
!             Generally should be 30 degrees or less for velocity, 360 for refl
!    rngmin   Minimum range (m) of data to use
!            (10 000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!    rngvvol  Range to gate in velocity 3-D volume
!    azmvvol  Azimuth angle in velocity 3-D volume
!    elvvvol  Elevation angle in velocity 3-D volume
!    varvol   Radar data 3-D volume
!
!    xs       x coordinate of scalar grid points in physical/comp. space (m)
!    ys       y coordinate of scalar grid points in physical/comp. space (m)
!    zps      Vertical coordinate of scalar grid points in physical space(m)
!
!  OUTPUT:
!    varvol   Radar data 3-D volume
!    istatus  Status indicator
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: nrefgates
  INTEGER, INTENT(IN) :: nvelgates
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev

  INTEGER, INTENT(IN) :: kntrgate(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntrazim(maxelev)
  INTEGER, INTENT(IN) :: kntrelev

  INTEGER, INTENT(IN) :: kntvgate(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntvazim(maxelev)
  INTEGER, INTENT(IN) :: kntvelev

  INTEGER, INTENT(IN) :: i_vcp
  INTEGER, INTENT(IN) :: gcopt
  REAL, INTENT(IN)    :: gcvrlim

  REAL, INTENT(IN)    :: winszrad
  REAL, INTENT(IN)    :: winszazim
  REAL, INTENT(IN)    :: refchek
  REAL, INTENT(IN)    :: velchek
  INTEGER, INTENT(IN)    :: irefgatsp
  INTEGER, INTENT(IN)    :: ivelgatsp
  REAL, INTENT(IN)    :: rngrvol(nrefgates,maxelev)
  REAL, INTENT(IN)    :: azmrvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvrvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: rngvvol(nvelgates,maxelev)
  REAL, INTENT(IN)    :: azmvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvvol(maxazim,maxelev)
  REAL, INTENT(INOUT) :: refvol(nrefgates,maxazim,maxelev)
  REAL, INTENT(INOUT) :: velvol(nvelgates,maxazim,maxelev)
  REAL, INTENT(INOUT) :: tmp(MAX(nrefgates,nvelgates),maxazim)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!

  INTEGER :: ii,jj,kk,kkv,kv,kr2,i,j,k,knt
  INTEGER :: mvelok,spin,flags
  INTEGER :: igate,ivgate,jazim,kelev,jazmin,jazmax
  INTEGER :: iigate,jjazim,jazmref,jazmvel
  INTEGER :: irngmin,irngmax
  INTEGER :: igatebgn,igateend,jazmbgn,jazmend,indgate,indazm
  INTEGER :: igspratio,iwinszrad,iwinszazm
  INTEGER :: kntchek,kntapref,kntapvel

  REAL :: azmdiff,maxdbz
  REAL :: summdbz,sumvel,sumvel2,sumtdbz,sumtvel
  REAL :: dbzdiff,sign,refprev,prevvel,veldiff
  REAL :: all_counts,spinchange_counts,signcnt
  REAL :: delev,avgelv,avgelvv,avgelvr,avgelvr2,appct
  REAL :: mdbz,tdbz,deltdbz,spinchange,stdvel,mvel,tvel
  REAL :: elvmax
  LOGICAL :: found

  REAL, PARAMETER :: elvmincmp = 1.1
! REAL, PARAMETER :: elvmincmp = 7.1
  REAL, PARAMETER :: dbzthr = 10.0
  REAL, PARAMETER :: apflag = -888.0
  REAL, PARAMETER :: gcflag = -890.0
  REAL, PARAMETER :: spinthr = 15.0
  REAL, PARAMETER :: ddbzthr = -12.0
  REAL, PARAMETER :: ddbzthr2 = -20.0
  REAL, PARAMETER :: mvelthr = 2.3
  REAL, PARAMETER :: apvelthr = 5.0
!
! Parameters applied for gcopt=2
!
  REAL, PARAMETER :: dbzclutterL = -5
  REAL, PARAMETER :: dbzclutterH = 200
  REAL, PARAMETER :: velclutter = 1.0
  REAL :: drnglim
  INTEGER ::  iii,jjj,kkk,ku
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Misc initializations
!
!-----------------------------------------------------------------------
!
  istatus=0

  IF( kntrelev > 1 ) THEN

    igspratio=max(int(float(irefgatsp)/float(2*ivelgatsp)),1)
    iwinszrad=max(int(winszrad/float(2*irefgatsp)),1)
    iwinszazm=max(int(0.5*winszazim),1)
!
!-----------------------------------------------------------------------
!
!  Establish the index of the reflectivity tilt above elvmincmp degrees.
!  For NEXRAD the first two unique tilt angles should be 0.5
!  and 1.5 degrees.
!
!-----------------------------------------------------------------------
!
    found=.false.
    avgelvr2=0.
    DO kk=2,kntrelev
      avgelvr2=0.
      DO jazim=1,kntrazim(kk)
        avgelvr2=avgelvr2+elvrvol(jazim,kk)
      END DO
      avgelvr2=avgelvr2/float(kntrazim(kk))
      print *, ' kk= ',kk,'  avgelvr2=',avgelvr2
      IF(avgelvr2 > elvmincmp) THEN
        found=.true.
        EXIT
      END IF
    END DO
    IF(found) THEN
    kr2=kk
    write(6,'(1x,a,f9.2,a,f9.2,a,i6)') 'Elev 1st ref tilt above ',     &
          elvmincmp,' deg: ',avgelvr2,' kr2=',kr2


!
!-----------------------------------------------------------------------
!
!  For now, AP detection is only done on any tilts below elvmincmp.
!
!-----------------------------------------------------------------------
!
    DO kk = 1, (kr2-1)

      tmp=0.
      kntchek=0
      kntapref=0
      kntapvel=0

!-----------------------------------------------------------------------
!
!  Establish the velocity tilt that most closely matches the
!  this reflectivity elevation angle.
!
!-----------------------------------------------------------------------
!
      avgelvr=0.
      DO jazim=1,kntrazim(kk)
        avgelvr=avgelvr+elvrvol(jazim,kk)
      END DO
      avgelvr=avgelvr/float(kntrazim(kk))

      kv=1
      delev=999.
      avgelvv=-99.
      DO kkv=1,kntvelev
        avgelv=0.
        DO jazim=1,kntvazim(kkv)
          avgelv=avgelv+elvvvol(jazim,kkv)
        END DO
        avgelv=avgelv/float(kntvazim(kkv))
        IF(abs(avgelv-avgelvr) < delev ) THEN
          delev=abs(avgelv-avgelvr)
          avgelvv=avgelv
          kv=kkv
        END IF
      END DO
      write(6,'(a,f9.2,a,i6)') ' Elev of nearest velocity tilt: ',      &
        avgelvv,' kv=',kv

      DO jjazim = 1,kntrazim(kk)
!
!  Find nearest azimuth in reflectivity data above elvmincmp
!
        jazmref=1
        azmdiff = 999.0
        DO jazim = 1,kntrazim(kr2)
          IF(abs(azmrvol(jazim,kr2)-azmrvol(jjazim,kk))<azmdiff)THEN
            azmdiff = abs(azmrvol(jazim,kr2)-azmrvol(jjazim,kk))
            jazmref = jazim
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth in velocity data
!
!-----------------------------------------------------------------------
!
        jazmvel=1
        azmdiff = 999.0
        DO jazim = 1,kntvazim(kv)
          IF(abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk)) < azmdiff        &
             .and. kntvgate(jazim,kv)>0 ) THEN
            azmdiff = abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk))
            jazmvel = jazim
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Calculate mean ref, texture of ref, spin
!
!-----------------------------------------------------------------------
!
        DO iigate= 1,kntrgate(jjazim,kk)
          IF(refvol(iigate,jjazim,kk) > refchek )THEN
            kntchek=kntchek+1
            summdbz=0.
            sumtdbz=0.
            all_counts=0.
            spinchange_counts=0.
            spin = 1

            irngmin = max((iigate - iwinszrad),1)
            irngmax = min((iigate + iwinszrad),kntrgate(jjazim,kk))
            jazmin = jjazim - iwinszazm
            jazmax = jjazim + iwinszazm
            IF(jazmin <= 0) jazmin = jazmin + kntrazim(kk)-2
            IF(jazmax > kntrazim(kk) )jazmax = jazmax-kntrazim(kk)+2
!
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin to jazmax
!
!-----------------------------------------------------------------------
!
            IF(jazmax > jazmin)THEN
              DO jazim = jazmin,jazmax
                DO igate = irngmin,irngmax
                  IF(refvol(igate,jazim,kk)>refchek)THEN
                    summdbz=summdbz+refvol(igate,jazim,kk)
                    sumtdbz=sumtdbz+dbzdiff*dbzdiff
                    all_counts = all_counts + 1
                    dbzdiff = refvol(igate,jazim,kk)-refprev
                    IF(dbzdiff > dbzthr)THEN
                      IF(spin<0)THEN
                        spinchange_counts = SPINchange_counts + 1
                        spin = 1
                      END IF
                    END IF
                    IF(dbzdiff < -dbzthr)THEN
                      IF(spin>0)THEN
                        spinchange_counts = SPINchange_counts + 1
                        spin = -1
                      END IF
                    END IF
                    IF(dbzdiff > 0.)THEN
                      signcnt = signcnt + 1.
                    ELSE
                      signcnt = signcnt - 1.
                    END IF
                    refprev = refvol(igate,jazim,kk)
                  END IF
                END DO ! igate
              END DO ! jazim
!
!-----------------------------------------------------------------------
!
! if jazmax<kazmin
!
!-----------------------------------------------------------------------
!
            ELSE
              DO jazim = jazmin,kntrazim(kk)
                DO igate = irngmin,irngmax
                  IF(refvol(igate,jazim,kk)>refchek)THEN
                    summdbz=summdbz+refvol(igate,jazim,kk)
                    sumtdbz=sumtdbz+dbzdiff*dbzdiff
                    all_counts = all_counts + 1
                    dbzdiff = refvol(igate,jazim,kk)-refprev
                    IF(dbzdiff > dbzthr)THEN
                      IF(spin<0)THEN
                        spinchange_counts = SPINchange_counts + 1
                        spin = 1
                      END IF
                    END IF
                    IF(dbzdiff < -dbzthr)THEN
                      IF(spin>0)THEN
                        spinchange_counts = SPINchange_counts + 1
                        spin = -1
                      END IF
                    END IF
                    IF(dbzdiff > 0.)THEN
                      signcnt = signcnt + 1.
                    ELSE
                      signcnt = signcnt - 1.
                    END IF
                    refprev = refvol(igate,jazim,kk)
                  END IF
                END DO ! igate
              END DO ! jazim
              DO jazim = 1,jazmax
                DO igate = irngmin,irngmax
                  IF(refvol(igate,jazim,kk)>refchek)THEN
                    summdbz=summdbz+refvol(igate,jazim,kk)
                    sumtdbz=sumtdbz+dbzdiff*dbzdiff
                    all_counts = all_counts + 1
                    dbzdiff = refvol(igate,jazim,kk)-refprev
                    IF(dbzdiff > dbzthr)THEN
                      IF(spin<0)THEN
                        spinchange_counts = SPINchange_counts + 1
                        spin = 1
                      END IF
                    END IF
                    IF(dbzdiff < -dbzthr)THEN
                      IF(spin>0)THEN
                        spinchange_counts = SPINchange_counts + 1
                        spin = -1
                      END IF
                    END IF
                    IF(dbzdiff > 0.)THEN
                      signcnt = signcnt + 1.
                    ELSE
                      signcnt = signcnt - 1.
                    END IF
                    refprev = refvol(igate,jazim,kk)
                  END IF
                END DO ! igate
              END DO ! jazim
            END IF

            IF(all_counts > 0.)THEN
              mdbz= summdbz/all_counts
              tdbz= sumtdbz/all_counts
              spinchange= (spinchange_counts/all_counts)*100.
            END IF
!
!  Calculate difference in reflectivity between this gate and the
!  first level above elvmincmp.
!
            IF(iigate <= kntrgate(jazmref,kr2) .and.   &
               refvol(iigate,jazmref,kr2) > refchek ) THEN
               deltdbz=                            &
                 refvol(iigate,jazmref,kr2)-refvol(iigate,jjazim,kk)
            ELSE
              IF(i_vcp==31 .or. i_vcp==32) THEN
                deltdbz= -32.0-refvol(iigate,jjazim,kk)
              ELSE
                deltdbz= -5.0-refvol(iigate,jjazim,kk)
              END IF
            END IF
!
!  If the delta-dBZ check fails when comparing to the nearest gate,
!  account for tilt of echoes by finding max of gates in the neighborhood.
!
            IF (deltdbz < ddbzthr) THEN
              jazmbgn = max((jazmref-1),1)
              jazmend = min((jazmref+1),kntrazim(kr2))
              maxdbz = refvol(iigate,jazmref,kr2)
              DO indazm = jazmbgn,jazmend
                igatebgn = max((iigate-1),1)
                igateend = min((iigate+1),kntrgate(indazm,kr2))
                DO indgate = igatebgn,igateend
                  maxdbz = max(maxdbz,refvol(indgate,indazm,kr2))
                END DO
              END DO
              IF( maxdbz > refchek )                                    &
                deltdbz = maxdbz - refvol(iigate,jjazim,kk)
              END IF
!
!-----------------------------------------------------------------------
!
!  Find std of vel, texture of vel, mean vel
!
!  First, find nearest velocity gate
!
!  Nearest velocity tilt (kv) and velocity azimuth (jazmvel) were
!  found earlier
!
!-----------------------------------------------------------------------
!
              DO igate=2,kntvgate(jazmvel,kv)-1
                IF(rngvvol(igate,kv) >= rngrvol(iigate,kk)) EXIT
              END DO
              IF(abs(rngvvol(igate,kv)-rngrvol(iigate,kk)) <            &
                 abs(rngvvol(igate-1,kv)-rngrvol(iigate,kk)) ) THEN
                ivgate=igate
              ELSE
                ivgate=igate-1
              END IF

              IF(abs(rngvvol(ivgate,kv)-rngrvol(iigate,kk))<            &
                 4*max(irefgatsp,ivelgatsp) .and.                       &
                 abs(azmvvol(jazmvel,kv)-azmrvol(jjazim,kk))<           &
                 float(iwinszazm))THEN
                sumvel2=0.
                sumvel =0.
                sumtvel=0.
                mvelok=0
                irngmin = max((ivgate-igspratio),1)
                irngmax = min((ivgate+igspratio),kntvgate(jazmvel,kv))
                jazmin = jazmvel - iwinszazm
                jazmax = jazmvel + iwinszazm
                IF( jazmin < 1 ) jazmin = jazmin+kntvazim(kv)
                IF( jazmax > kntvazim(kv) ) jazmax = jazmax-kntvazim(kv)
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin to jazmax
!
!-----------------------------------------------------------------------
!
                IF(jazmax > jazmin)THEN
                  DO jazim = jazmin,jazmax
                    prevvel = velvol(irngmin,jazim,kv)
                    DO igate = irngmin+1,irngmax
                      IF(velvol(igate,jazim,kv)>velchek)THEN
                        sumvel=sumvel+velvol(igate,jazim,kv)
                        sumvel2=sumvel2+                                &
                          velvol(igate,jazim,kv)*velvol(igate,jazim,kv)
                        mvelok=mvelok+1
                        veldiff=velvol(igate,jazim,kv) - prevvel
                        sumtvel=sumtvel+veldiff*veldiff
                        prevvel=velvol(igate,jazim,kv)
                      END IF
                  END DO ! igate
                END DO ! jazim
!
!-----------------------------------------------------------------------
!
! if jazmax<kazmin
!
!-----------------------------------------------------------------------
!
              ELSE
                DO jazim = jazmin,kntvazim(kv)
                  prevvel = velvol(irngmin,jazim,kv)
                  DO igate = irngmin+1,irngmax
                    IF(velvol(igate,jazim,kv)>velchek)THEN
                      sumvel=sumvel+velvol(igate,jazim,kv)
                      sumvel2=sumvel2+                                  &
                        velvol(igate,jazim,kv)*velvol(igate,jazim,kv)
                      mvelok=mvelok+1
                      veldiff=velvol(igate,jazim,kv) - prevvel
                      sumtvel=sumtvel+veldiff*veldiff
                      prevvel=velvol(igate,jazim,kv)
                    END IF
                  END DO ! igate
                END DO ! jazim
                DO jazim = 1,jazmax
                  prevvel = velvol(irngmin,jazim,kv)
                  DO igate = irngmin,irngmax
                    IF(velvol(igate,jazim,kv)>velchek)THEN
                      sumvel=sumvel+velvol(igate,jazim,kv)
                      sumvel2=sumvel2+                                   &
                        velvol(igate,jazim,kv)*velvol(igate,jazim,kv)
                      mvelok=mvelok+1
                      veldiff= velvol(igate,jazim,kv) - prevvel
                      sumtvel=sumtvel+veldiff*veldiff
                      prevvel = velvol(igate,jazim,kv)
                    END IF
                  END DO ! igate
                END DO ! jazim
              END IF
              IF(mvelok>1)THEN
                stdvel = (sumvel2-sumvel*sumvel/mvelok)/(mvelok-1)
                IF(stdvel .lt. 0.0) stdvel = 0.0
                stdvel = sqrt(stdvel)
                tvel= sumtvel/mvelok
                mvel= abs(sumvel/mvelok)
              ELSE
                stdvel=999.
                tvel=999.
                mvel=999.
              END IF
            ELSE  ! if ivgate is too far away from iigate
              stdvel=999.
              tvel=999.
              mvel=999.
            END IF

!
!  If calculated values match two or more characteristics,
!  set the AP marker, tmp=1.  The actual resetting of refvol
!  is deferred until all points have been checked on this level,
!  so that the calculation of spin at subsequent gates is not affected.
!

            flags=0
            IF(spinchange >= spinthr) flags=flags+1
            IF(deltdbz <= ddbzthr) flags=flags+1
            IF(deltdbz <= ddbzthr2) flags=flags+1
            IF(abs(mvel) < mvelthr ) flags=flags+1

            IF(flags > 1 ) tmp(iigate,jjazim)=1.0

          END IF  ! a valid reflectivity at iigate,jjazim,kk
        END DO  !iigate
      END DO  !jjazim

!
!   Where the AP marker has been set, set for reflectivity to apflag,
!   and set corresponding velocities to apflag
!

      DO jjazim=1,kntrazim(kk)
!
        jazmvel=1
        azmdiff = 999.0
        DO jazim = 1,kntvazim(kv)
          IF(abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk))<azmdiff)THEN
            azmdiff = abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk))
            jazmvel = jazim
          END IF
        END DO
!
        DO iigate=1,kntrgate(jjazim,kk)
          IF(tmp(iigate,jjazim) > 0. ) THEN
            kntapref=kntapref+1
            refvol(iigate,jjazim,kk) = apflag
            DO igate=2,kntvgate(jazmvel,kv)-1
              IF(rngvvol(igate,kv) >= rngrvol(iigate,kk)) EXIT
            END DO

            IF(abs(rngvvol(igate,kv)-rngrvol(iigate,kk)) <              &
               abs(rngvvol(igate-1,kv)-rngrvol(iigate,kk)) ) THEN
              ivgate=igate
            ELSE
              ivgate=igate-1
            END IF

            IF(abs(rngvvol(ivgate,kv)-rngrvol(iigate,kk))<              &
               4*max(irefgatsp,ivelgatsp) .and.                         &
               abs(azmvvol(jazmvel,kv)-azmrvol(jjazim,kk))<             &
               float(iwinszazm))THEN
            irngmin = ivgate - int(irefgatsp/float(2*ivelgatsp))
            irngmin = max(irngmin,1)
            irngmax = ivgate + int(irefgatsp/float(2*ivelgatsp))
            irngmax = min(irngmax,kntvgate(jazmvel,kv))
            DO igate=irngmin,irngmax
              kntapvel=kntapvel+1
              velvol(igate,jazmvel,kv)=apflag
            END DO
            ENDIF
          END IF
        END DO  ! iigate loop
      END DO  ! jjazim loop
!
!  For gcopt=2
!     Zhao Kun Add Residual AP clutter Remove 2008.12.9
!
      drnglim=0.5*irefgatsp
      IF(gcopt == 2) THEN
        DO ku=1,kntrelev
          DO jj=1,kntrazim(ku)
            avgelvr2= elvrvol(jj,kk)
            DO ii = 1,kntrgate(jj,ku)
!
!             If the reflectivity is suspected to be clutter or
!             AP clutter then contintue check velocity
!
              IF(    refvol(ii,jj,ku)  < dbzclutterL .OR.               &
                 abs(refvol(ii,jj,ku)) > dbzclutterH) THEN
!               print*,'zhaokun',refvol(ii,jj,ku)
                refvol(ii,jj,ku)=apflag
                DO jjj = 1,kntvazim(ku)
                  IF(abs(azmrvol(jj,ku)-azmvvol(jjj,ku)) < 1.0 .AND.    &
                     abs(elvrvol(jj,ku)-elvvvol(jjj,ku)) < 0.3 ) THEN
                    DO iii = 1,kntvgate(jjj,ku)
                      IF(abs(rngrvol(ii,ku)-rngvvol(iii,ku)) <=         &
                                                        drnglim) THEN
                        IF(abs(velvol(iii,jjj,ku))< velclutter)         &
                          velvol(iii,jjj,ku)= apflag
                      END IF
                    END DO
                  END IF
                END DO

              END IF
            END DO
          END DO
        END DO
      END IF   ! gcopt == 2

      IF ( kntchek > 0 ) THEN
        appct=100.*float(kntapref)/float(kntchek)
      ELSE
        appct=0.
      END IF
      write(6,'(a,i6,a,f6.2,/a,i8,/a,i8,f9.2,a,/a,i8)')                 &
       ' AP detect completed for level ',kk,' elev=',avgelvr,           &
       '   Reflectivity gates checked:',kntchek,                        &
       '               AP flagged ref:',kntapref,appct,' percent',      &
       '               AP flagged vel:',kntapvel

!
!   Apply despekl to the edited data
!   This should help catch AP residue.
!   Note: here force no median filter.
!
      CALL despekl(nrefgates,maxazim,nrefgates,maxazim,refchek,         &
                   0,refvol(1,1,kk),tmp)
    END DO  ! kk

! open(21,file='ap.dat',form='unformatted',status='unknown')
! write(21) nrefgates
! write(21) maxazim
! write(21) maxelev
! write(21) kntrelev
! write(21) kntrazim
! write(21) kntrgate
! write(21) rngrvol
! write(21) azmvvol
! write(21) elvrvol
! write(21) refvol
! close(21)
    ELSE
      write(6,'(1x,a,f9.2,a,f9.2,/3x,a//)') 'Highest elev not above ', &
         elvmincmp,' deg: ',avgelvr2,'Skipping AP Detect'
    END IF

  ELSE
    write(6,'(a)') ' Need at least two reflectivity levels for AP detection'
    write(6,'(a)') ' Skipping AP detect'
  END IF

  IF(gcopt > 0) THEN
    DO kk=1,kntvelev
      elvmax=-99.
      DO jazim=1,kntvazim(kk)
        elvmax=max(elvmax,elvvvol(jazim,kk))
      END DO

      tmp = 0.0
      DO jazim=1,kntvazim(kk)
        DO igate=1, kntvgate(jazim,kk)
          tmp(igate,jazim) = velvol(igate,jazim,kk)
        END DO
      END DO

      DO jazim=2,kntvazim(kk)-1
        DO igate=2, kntvgate(jazim,kk)-1
          IF( abs(tmp(igate+1,jazim)) < 1e-1 .AND.        &
              abs(tmp(igate-1,jazim)) < 1e-1 .AND.        &
              abs(tmp(igate,jazim+1)) < 1e-1 .AND.        &
              abs(tmp(igate,jazim-1)) < 1e-1 )    THEN
                velvol(igate,jazim,kk) = gcflag
          END IF
        END DO
      END DO

      IF(elvmax < elvmincmp) THEN
        DO jazim=1,kntvazim(kk)
          IF(elvvvol(jazim,kk) < elvmincmp) THEN
            DO igate=1,kntvgate(jazim,kk)
              IF(abs(velvol(igate,jazim,kk)) < gcvrlim) THEN
                velvol(igate,jazim,kk) = gcflag
              END IF
            END DO
          END IF
        END DO
      END IF
    END DO
  END IF   ! gcchkopt

  RETURN
END SUBROUTINE apdetect
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ANOMRADIAL                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE anomradial(maxgate,maxazim,ngate,nazim,                      &
                      rfrst_ref,gtspc_ref,refchek,anrflag,azim,refl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Detect and flag anomalous radials that are commonly found with
!  solar interference (sunrise, sunset) and beam blockage.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!          December, 2010
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate  Maximum number of gates in a radial
!    maxazim  Maximum number of radials in a tilt
!    anrflag  Anomalous radial flag
!    refl     Reflectivity
!
!  OUTPUT:
!    refl     Reflectivity
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim

  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: nazim

  INTEGER, INTENT(IN) :: rfrst_ref
  INTEGER, INTENT(IN) :: gtspc_ref

  REAL, INTENT(IN)    :: refchek
  REAL, INTENT(IN)    :: anrflag

  REAL, INTENT(IN)    :: azim(maxazim)
  REAL, INTENT(INOUT) :: refl(maxgate,maxazim)
!
!-----------------------------------------------------------------------
!
! Linear least-squares variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: eps=1.0E-20
  REAL, PARAMETER :: rminlsq=50.0E03
  REAL :: reflin(2,2),refrhs(2),sol(2)
  REAL :: work(2,3),work1d(3)
  REAL :: refcst,drefdr,dref,sum2,rms
  INTEGER :: nref,istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: refnull = 0.0

  INTEGER, ALLOCATABLE :: kntref(:)
  REAL, ALLOCATABLE :: refmean(:)
  REAL, ALLOCATABLE :: refstdv(:)
  REAL :: refrng,sumref,sum2ref
  REAL :: sumrefglb,sum2refglb
  REAL :: refmeanglb,refstdvglb,fngate,fngatem1,vldratio,vldratglb
  REAL :: refvar,refnull2,fngatetot,fngatetotm1
  INTEGER :: i,j,jj,k,kk,kntrefglb,istat
  INTEGER :: kntazmflg,ngatetot
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(kntref(maxazim),stat=istat)
  CALL check_alloc_status(istat,'anomradial:kntref')
  ALLOCATE(refmean(maxazim),stat=istat)
  CALL check_alloc_status(istat,'anomradial:refmean')
  ALLOCATE(refstdv(maxazim),stat=istat)
  CALL check_alloc_status(istat,'anomradial:refstdv')
!
! Compute global statitiscs on this sweep plane
!
  print *, ' rfrst_ref:',rfrst_ref
  print *, ' gtspc_ref:',gtspc_ref
  kntrefglb=0
  sumrefglb=0.
  sum2refglb=0.
  refnull2=refnull*refnull
  fngate=float(ngate)
  fngatem1=float(ngate-1)
  ngatetot=ngate*nazim
  fngatetot=float(ngatetot)
  fngatetotm1=float(ngatetot-1)
  DO j=1,nazim
    kntref(j)=0
    sumref=0.
    sum2ref=0.
    DO i=1,ngate
      IF(refl(i,j) > refnull ) THEN
        kntref(j)=kntref(j)+1
        sumref=sumref+refl(i,j)
        sum2ref=sum2ref+(refl(i,j)*refl(i,j))
      ELSE
        sumref=sumref+refnull
        sum2ref=sum2ref+refnull2
      END IF
    END DO
    kntrefglb=kntrefglb+kntref(j)
    sumrefglb=sumrefglb+sumref
    sum2refglb=sum2refglb+sum2ref
    refmean(j)=sumref/fngate
    refvar=(sum2ref-(sumref*sumref/fngate))/fngatem1
    refstdv(j)=SQRT(refvar)
  END DO
!
  kntazmflg=0
  WRITE(6,'(a,i6)') ' ngate this elevation:',ngate
  vldratglb=float(kntrefglb)/float(nazim*ngate)
  WRITE(6,'(a,f7.0,a)') ' Global ref>0 gates:',(vldratglb*100.),' percent.'
  IF( kntrefglb > 100 ) THEN
    refmeanglb=sumrefglb/fngatetot
    refvar=(sum2refglb-(sumrefglb*sumrefglb/fngatetot))/fngatetotm1
    refstdvglb=SQRT(refvar)
    WRITE(6,'(3(a,f5.1))') ' Global refl mean:',refmeanglb,' Std dev:',refstdvglb
    IF( refmeanglb > 0. .OR. vldratglb > 0.10) THEN
!
!  Precip situation look for a mostly non-zero radial with positive linear trend
!
      DO j=1,nazim
!       WRITE(6,'(a,f7.1,a,i6,2(a,f7.1))')                                      &
!             ' azim: ',azim(j),' kntref: ',kntref(j),                          &
!             ' refmean: ',refmean(j),'  refstdv: ',refstdv(j)
        vldratio=float(kntref(j))/fngate
        IF(vldratio > 0.5 ) THEN  ! mostly non-missing radial
!
!  Least-squares linear trend
!  refl = a + b*range
!
          reflin=0.
          refrhs=0.
          DO i=1,ngate
            refrng=float(rfrst_ref+(i-1)*gtspc_ref)
            IF(refrng > rminlsq .AND. refl(i,j) > refchek) THEN
              reflin(1,1)=reflin(1,1)+1.0
              reflin(1,2)=reflin(1,2)+refrng
              reflin(2,1)=reflin(2,1)+refrng
              reflin(2,2)=reflin(2,2)+refrng*refrng
              refrhs(1)=refrhs(1)+refl(i,j)
              refrhs(2)=refrhs(2)+refl(i,j)*refrng
            END IF
          END DO

          CALL gjelim(2,reflin,refrhs,sol,work,work1d,eps,istatus)

          nref=0
          sum2=0.
          refcst=sol(1)
          drefdr=sol(2)
          DO i=1,ngate
            refrng=float(rfrst_ref+(i-1)*gtspc_ref)
            IF(refrng > rminlsq .AND. refl(i,j) > refchek) THEN
              nref=nref+1
              dref=(refcst+drefdr*refrng) - refl(i,j)
              sum2=sum2+(dref*dref)
            END IF
          END DO
          rms=sqrt(sum2/float(nref))
!         WRITE(6,'(a,f7.1,2x,g12.6,2x,f7.2)')                        &
!             ' Anomradial LSQ: refcst,drefdr,rms: ',refcst,drefdr,rms
          IF( drefdr > 2.0E-05 .AND. drefdr < 9.0E-05 .AND. rms < 4.0) THEN

            WRITE(6,'(a,f7.1,a)') ' *Anomalous radial detected (a).  Azim: ',azim(j),' degrees'
            WRITE(6,'(a,i6,a,f7.3)') '   N valid gates: ',kntref(j),'  vldratio: ',vldratio
!
!   Flag reflectivity gates in this radial
!
            kntazmflg=kntazmflg+1
            DO i=1,ngate
              IF(refl(i,j) > refchek) THEN
                refl(i,j)=anrflag
              END IF
            END DO
          END IF
        END IF
      END DO
    ELSE  ! clear situation, look for a high valid ratio
      DO j=1,nazim
!       WRITE(6,'(a,f7.1,a,i6,2(a,f7.1))')                                           &
!             ' azim: ',azim(j),' kntref: ',kntref(j),                          &
!             ' refmean: ',refmean(j),'  refstdv: ',refstdv(j)
        vldratio=float(kntref(j))/fngate
        IF(vldratio > 0.40 .AND. (refmean(j)-refmeanglb) > 5.0 ) THEN
          reflin=0.
          refrhs=0.
          DO i=1,ngate
            refrng=float(rfrst_ref+(i-1)*gtspc_ref)
            IF(refrng > rminlsq .AND. refl(i,j) > refchek) THEN
              reflin(1,1)=reflin(1,1)+1.0
              reflin(1,2)=reflin(1,2)+refrng
              reflin(2,1)=reflin(2,1)+refrng
              reflin(2,2)=reflin(2,2)+refrng*refrng
              refrhs(1)=refrhs(1)+refl(i,j)
              refrhs(2)=refrhs(2)+refl(i,j)*refrng
            END IF
          END DO

          CALL gjelim(2,reflin,refrhs,sol,work,work1d,eps,istatus)

          nref=0
          sum2=0.
          refcst=sol(1)
          drefdr=sol(2)
          DO i=1,ngate
            refrng=float(rfrst_ref+(i-1)*gtspc_ref)
            IF(refrng > rminlsq .AND. refl(i,j) > refchek) THEN
              nref=nref+1
              dref=(refcst+drefdr*refrng) - refl(i,j)
              sum2=sum2+(dref*dref)
            END IF
          END DO
          rms=sqrt(sum2/float(nref))
!         WRITE(6,'(a,f7.1,2x,g12.6,2x,f7.2)')                        &
!             ' Anomradial LSQ: refcst,drefdr,rms: ',refcst,drefdr,rms
          IF( drefdr > 2.0E-05 .AND. drefdr < 9.0E-05 .AND. rms < 4.0) THEN

            WRITE(6,'(a,f7.1,a)') ' Anamalous radial detected (b).  Azim: ',azim(j),' degrees'
            WRITE(6,'(a,i6,a,f7.3)') '   N valid gates: ',kntref(j),'  vldratio: ',vldratio
!           WRITE(6,'(a,f7.2,a,f7.2)') '   Ref mean:    ',refmean(j),'  Global mean:   ',refmeanglb
!           WRITE(6,'(a,f7.2,a,f7.2)') '   Ref std dev: ',refstdv(j),'  Global std dev:',refstdvglb
!
!   Flag reflectivity gates in this radial
!
            kntazmflg=kntazmflg+1
            DO i=1,ngate
              IF(refl(i,j) > refchek) THEN
                refl(i,j)=anrflag
              END IF
            END DO
          END IF
        END IF
      END DO
    END IF
  END IF

  WRITE(6,'(a,i4,a,i4,a)') ' Anomradial: ',kntazmflg, &
        ' radials flagged of ',nazim,' radials this level'
!
  DEALLOCATE(kntref)
  DEALLOCATE(refmean)
  DEALLOCATE(refstdv)
!
  RETURN
!
END SUBROUTINE anomradial
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE DESPEKL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE despekl(maxgate,maxazim,                                     &
                   ngate,nazim,varchek,medfilt,rdata,tmp)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Despeckle radar data.
!  Can be used for reflectivity or velocity.
!
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    maxgate  Maximum number of gates in a radial
!    maxazim  Maximum number of radials in a tilt
!    ngate    Number of gates in radial
!    nazim    Number of radials
!    varchek  Threshold to determine good data vs. flagged
!    medfilt  Median filter option switch
!    rdata    Doppler radial velocity
!
!  OUTPUT:
!    rdata
!
!  WORK ARRAYS:
!
!    tmp
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: nazim

  REAL,    INTENT(IN) :: varchek
  INTEGER, INTENT(IN) :: medfilt

  REAL, INTENT(INOUT) :: rdata(maxgate,maxazim)
  REAL, INTENT(OUT)   :: tmp(maxgate,maxazim)

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  INTEGER :: NN,m,n,loc,is
  INTEGER :: kntgd,kntdsp
  INTEGER :: nprint
!
  REAL :: sortdata(9)
  REAL :: swp
!
  REAL :: sum,pctdsp
  REAL, PARAMETER :: dspmiss = -991.
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
! Create temporary array where tmp=1 if non-missing tmp=0 if missing
!
!-----------------------------------------------------------------------
!
!
  tmp=0.
  DO j=1,nazim
    DO i=1,ngate
      IF ( rdata(i,j) > varchek ) tmp(i,j)=1.
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! If datum has fewer than 3 non-missing neighbors in a 3x3
! square template set it to missing
!
!-----------------------------------------------------------------------
!
  kntgd=0
  kntdsp=0
  DO j=2,nazim-1
    DO i=2,ngate-1
      IF(rdata(i,j) > varchek ) THEN
        kntgd=kntgd+1
        sum=tmp(i-1,j+1)+tmp(i,j+1)+tmp(i+1,j+1)          &
           +tmp(i-1,j  )+           tmp(i+1,j)            &
           +tmp(i-1,j-1)+tmp(i,j-1)+tmp(i+1,j-1)
        IF( sum < 3. ) THEN
          kntdsp=kntdsp+1
          rdata(i,j) = dspmiss
        END IF
      END IF
    END DO
  END DO
  IF(kntgd > 0 ) THEN
    pctdsp=100.*(float(kntdsp)/float(kntgd))
  ELSE
    pctdsp=0.
  END IF

  write(6,'(a,i8,a,i8,a,f6.1,a)') &
    ' Despeckled ',kntdsp,' of ',kntgd,' data =',pctdsp,' percent.'

!----------------------------------------------------------------------------
!
!  Zhao Kun added a median filter to smooth the reflectivity and velocity
!  Recoded by Keith Brewster to make it optional and streamline a bit.
!
!----------------------------------------------------------------------------

  IF(medfilt > 0) THEN
!   nprint=0
    DO j=1,nazim
      DO i=1,ngate
        tmp(i,j)=rdata(i,j)
      END DO
    END DO

    DO j=2,nazim-1
      DO i=2,ngate-1
        NN=0
        DO m=-1,1
          DO n=-1,1
            swp=tmp(i+m,j+n)
            IF(swp > varchek) THEN
              NN=NN+1
              IF(NN == 1) THEN ! first
                sortdata(NN)=swp
              ELSE
                DO loc=1,NN-1
                  IF(swp < sortdata(loc)) THEN
                    DO is=NN,loc+1,-1
                      sortdata(is)=sortdata(is-1)
                    END DO
                    EXIT
                  END IF
                END DO
                sortdata(loc)=swp
              END IF
            END IF
          END DO
        END DO

        IF (NN > 6) THEN
!         IF(nprint < 100) THEN
!           print *, ' NN: ',NN
!           print *, ' sortdata: ',(sortdata(m),m=1,NN)
!           print*,'median',sortdata((NN+1)/2),rdata(i,j),tmp(i,j)
!           nprint=nprint+1
!         END IF
          rdata(i,j)=sortdata((NN+1)/2)

        END IF

      END DO
    END DO
!
!  First radial
!
    j=1
    DO i=2,ngate-1
      NN=0
      DO m=-1,1
        DO n=0,2
          swp=tmp(i+m,j+n)
          IF(swp > varchek) THEN
            NN=NN+1
            IF(NN == 1) THEN ! first
              sortdata(NN)=swp
            ELSE
              DO loc=1,NN-1
                IF(swp < sortdata(loc)) THEN
                  DO is=NN,loc+1,-1
                    sortdata(is)=sortdata(is-1)
                  END DO
                  EXIT
                END IF
              END DO
              sortdata(loc)=swp
            END IF
          END IF
        END DO
      END DO

      IF (NN > 6) THEN

!       DO m=1,NN
!         print*,'NN',m,sortdata(m)
!       ENDDO
!       print*,'mean',sortdata((NN+1)/2),rdata(i,j),tmp(i,j)
        rdata(i,j)=sortdata((NN+1)/2)
!       stop

      END IF

    END DO
!
!  Last radial
!
    j=nazim
    DO i=2,ngate-1
      NN=0
      DO m=-1,1
        DO n=-2,0
          swp=tmp(i+m,j+n)
          IF(swp > varchek) THEN
            NN=NN+1
            IF(NN == 1) THEN ! first
              sortdata(NN)=swp
            ELSE
              DO loc=1,NN-1
                IF(swp < sortdata(loc)) THEN
                  DO is=NN,loc+1,-1
                    sortdata(is)=sortdata(is-1)
                  END DO
                  EXIT
                END IF
              END DO
              sortdata(loc)=swp
            END IF
          END IF
        END DO
      END DO

      IF (NN > 6) THEN

!       DO m=1,NN
!         print*,'NN',m,sortdata(m)
!       ENDDO
!       print*,'mean',sortdata((NN+1)/2),rdata(i,j),tmp(i,j)
        rdata(i,j)=sortdata((NN+1)/2)
!       stop

      END IF

    END DO
  END IF ! median filter

  RETURN
!
END SUBROUTINE despekl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SPWTHR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE spwthr(maxgate,maxazim,                                      &
                  ngate,nazim,rvchek,spwthrrat,spwflag,vnyq,rvel,spw)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Threshold radial velocity data based on spectrum width.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!  Keith Brewster (22-June-2011)
!  Added spwthrrat to the argument list.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate  Maximum number of gates in a radial
!    maxazim  Maximum number of radials in a tilt
!    ngate    Number of gates in radial
!    nazim    Number of radials
!    rvchek   Threshold for good data vs. flagged
!    spwthrrat Spectrum width threshold as a ratio of Nyquist velocity
!    spwflag  Flag for spectrum width threshold fail
!    vnyq     Nyquist velocity
!    rvel     Doppler radial velocity
!    spw      Doppler spectrum width
!
!  OUTPUT:
!    rvel     Doppler radial velocity
!
!
!  WORK ARRAYS:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: nazim

  REAL, INTENT(IN)    :: rvchek
  REAL, INTENT(IN)    :: spwthrrat
  REAL, INTENT(IN)    :: spwflag
  REAL, INTENT(IN)    :: vnyq(maxazim)

  REAL, INTENT(INOUT) :: rvel(maxgate,maxazim)
  REAL, INTENT(IN)    :: spw(maxgate,maxazim)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: spwthresh
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
! If spectrum width is greater than threshold, set to missing
!
!-----------------------------------------------------------------------
!
  DO j=1,nazim
    IF(vnyq(j) > 0.) THEN
      spwthresh=spwthrrat*vnyq(j)
    DO i=1,ngate
      IF(rvel(i,j) > rvchek .AND. spw(i,j) > spwthresh ) rvel(i,j) = spwflag
    END DO
    END IF
  END DO
!
  RETURN
!
END SUBROUTINE spwthr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE RHOHVCHK                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rhohvchk(maxrgate,maxvgate,maxazim,ng_ref,ng_vel,nazim, &
                    irefgatsp,ivelgatsp,refchek,velchek,rhohvthr,  &
                    refl,radv,rhohv)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Threshold radar data based on Rho-HV
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!  3/19/2013 (Keith Brewster)
!  Updated to include different gate spacing and gate counts 
!  for reflectivity and vel.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxrgate   Maximum number of reflectivity gates in a radial
!    maxvgate   Maximum number of velocity gates in a radial
!    maxazim    Maximum number of radials in a tilt
!    ngate      Number of reflectivity gates in radial
!    ngate      Number of reflectivity gates in radial
!    nazim      Number of radials
!    irefgatsp  Gate spacing of reflectivity data
!    ivelgatsp  Gate spacing of velocity data
!    refchek    Threshold for good data vs. flagged reflectivity
!    velchek    Threshold for good data vs. flagged velocity
!    rhohvthr   Rho-HV threshold
!    refl       Reflectivity
!    rvel       Doppler radial velocity
!    rhohv      Rho-HV
!
!  OUTPUT:
!    refl     Reflectivity
!    rvel     Doppler radial velocity
!
!  WORK ARRAYS:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: maxrgate
  INTEGER, INTENT(IN) :: maxvgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: ng_ref
  INTEGER, INTENT(IN) :: ng_vel
  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: irefgatsp
  INTEGER, INTENT(IN) :: ivelgatsp
  REAL, INTENT(IN) :: refchek
  REAL, INTENT(IN) :: velchek
  REAL, INTENT(IN) :: rhohvthr
  REAL, INTENT(INOUT) :: refl(maxrgate,maxazim)
  REAL, INTENT(INOUT) :: radv(maxvgate,maxazim)
  REAL, INTENT(IN) :: rhohv(maxrgate,maxazim)
!
  REAL, PARAMETER :: refmiss=-901.0
  REAL, PARAMETER :: velmiss=-901.0
  INTEGER, PARAMETER :: nhist=20
  INTEGER :: histknt(nhist)
  INTEGER :: igate,iigate,jazim,ihist
  INTEGER :: kntvalid,kntflag
  INTEGER :: gtspratio,ngincr,iibgn,iiend
  REAL :: histcst,histinv,histbin,pctflg,pcthist
!
! tdsrefmin is minimum reflectivity (dBZ) of tornado debris signature (TDS)
! TDS will have low Rho-HV and we don't want to flag it bad.
! Same for some large tumbling hail.
!
  REAL, PARAMETER :: tdsrefmin = 46.0
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  histknt=0
  histcst=1.0/float(nhist)
  histinv=float(nhist)

  kntvalid=0
  kntflag=0
  IF(ng_vel < 1) THEN   ! no velocity data
    DO jazim=1,nazim
      DO igate=1,ng_ref
        IF(refl(igate,jazim) > refchek ) THEN
          kntvalid=kntvalid+1
          IF( rhohv(igate,jazim) < rhohvthr .AND.                      &
              refl(igate,jazim) < tdsrefmin ) THEN
            kntflag=kntflag+1
            refl(igate,jazim) = refmiss
          END IF
          ihist=1+INT(rhohv(igate,jazim)*histinv)
          ihist=min(max(ihist,1),nhist)
          histknt(ihist)=histknt(ihist)+1
        END IF
      END DO
    END DO
  ELSE 
    print *, ' rhohvchk: ng_ref, ng_vel: ',ng_ref,ng_vel
    print *, ' rhohvchk: irefgatsp,ivelgatsp: ',irefgatsp,ivelgatsp
    gtspratio=NINT(FLOAT(irefgatsp)/FLOAT(ivelgatsp))
    WRITE(6,'(1x,a,i4)') 'rhohvchk: gate space ratio =',gtspratio
    IF(gtspratio == 1 .AND. ng_vel >= ng_ref) THEN  ! vel and ref same gates
      DO jazim=1,nazim
        DO igate=1,ng_ref
          IF(refl(igate,jazim) > refchek ) THEN
            kntvalid=kntvalid+1
            IF( rhohv(igate,jazim) < rhohvthr .AND.                    &
                refl(igate,jazim) < tdsrefmin ) THEN
              kntflag=kntflag+1
              refl(igate,jazim) = refmiss
              IF(radv(igate,jazim) > velchek)                          &
                 radv(igate,jazim) = velmiss
            END IF
            ihist=1+INT(rhohv(igate,jazim)*histinv)
            ihist=min(max(ihist,1),nhist)
            histknt(ihist)=histknt(ihist)+1
          END IF
        END DO
      END DO
    ELSE ! vel and ref different gate spacing
      ngincr=gtspratio-1
      DO jazim=1,nazim
        DO igate=1,ng_ref
          IF(refl(igate,jazim) > refchek ) THEN
            kntvalid=kntvalid+1
            IF( rhohv(igate,jazim) < rhohvthr .AND.                    &
                refl(igate,jazim) < tdsrefmin ) THEN
              kntflag=kntflag+1
              refl(igate,jazim) = refmiss
              iibgn=((igate-1)*gtspratio)+1
              iiend=MIN((iibgn+ngincr),ng_vel)
              DO iigate=iibgn,iiend
                IF(radv(iigate,jazim) > velchek)                       &
                   radv(iigate,jazim) = velmiss
              END DO
            END IF
            ihist=1+INT(rhohv(igate,jazim)*histinv)
            ihist=min(max(ihist,1),nhist)
            histknt(ihist)=histknt(ihist)+1
          END IF
        END DO
      END DO
    END IF
  END IF
!
  pctflg=0.
  IF(kntvalid > 0) THEN
    pctflg=100.*float(kntflag)/float(kntvalid)
  END IF
  WRITE(6,'(a,i9,a,i9,a,f5.1,a)') &
   ' Rho-HV screened ',kntflag,' of ',kntvalid,' = ',pctflg, &
   '% reflectivity gates'
!
  WRITE(6,'(/a)') ' Histogram of Rho-HV (count, percent):'
  DO ihist=1,nhist
    histbin=histcst*(ihist-1)
    pcthist=100.*float(histknt(ihist))/float(kntvalid)
    WRITE(6,'(f7.3,i9,f10.1)') histbin,histknt(ihist),pcthist
  END DO
  RETURN
END SUBROUTINE rhohvchk
