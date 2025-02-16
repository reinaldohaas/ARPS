!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UNFOLDNQC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE unfoldnqc(maxgate,maxazim,ngate,nazim,                     &
                  nzsnd,zsnd,usnd,vsnd,rfrsnd,                        &
                  bkgopt,shropt,rfropt,                               &
                  igatesp,irfirst,iangle,itilt,                       &
                  iyear,imon,iday,ihr,imin,isec,                      &
                  radlat,radlon,radalt,                               &
                  rvel,spw,elev,azim,nqvel,                           &
                  unfvel,bkgvel,bgate,egate,tmp1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Unfold and Quality Control Doppler radial velocities.
!
!  AUTHOR:
!  Keith Brewster, Center for Analysis and Prediction of Storms
!
!  MODIFICATION HISTORY:   Keith Brewster
!  Nov. 2000, First Version
!  Oct. 2001, Updated and added spectrum width and std dev thresholding
!  Jan  2002, Modified to check vs background wind first, added switches
!             for comparison to background and for shear checking
!  Jan  2009, Modified to use Nyquist for each radial (Youngsun Jung)
!  Aug  2010  Added additional diagnostic output.
!             Added variables to argument list to pass along to wrtvel2. 
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate  Maximum number of gates in a radial
!    maxazim  Maximum number of radials in a tilt
!    ngate    Number of gates in radial
!    nazim    Number of radials
!    nzsnd  Number of levels in background wind profile
!    zsnd     Height (m MSL) of winds in background wind profile
!    usnd     U wind component in background wind profile
!    vsnd     V wind component in background wind profile
!    rfrsnd   Refractivity from background temp and moisture profile
!    bkgopt   Option to check data vs. a background profile (0:no,1:yes)
!    shropt   Option to check data using local shear check (0:no,1:yes)
!    rfropt
!    igatesp  Velocity gate spacing (m)
!    irfirst  Range to first gate (m)
!    iangle
!    itilt
!    iyear,imon,iday,ihr,imin,isec
!    radlat   Latitude (degrees N) of radar
!    radlon   Longitude (degrees E) of radar
!    radalt   Altitude (m MSL) of radar
!    rvel     Doppler radial velocity
!    spw      Doppler spectrum width
!    elev     Elevation angle
!    azim     Azimuth
!    nqvel    Nyquist velocity for each radial
!
!  OUTPUT:
!    unfvel   Unfolded Doppler radial velocity
!
!  WORK ARRAYS:
!    bkgvel   Radial velocity of background wind
!    bgate    Counter to first valid gate in each ray
!    bgate    Counter to last valid gate in each ray
!    tmp1     Quality/Missing indicator
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
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN) :: zsnd(nzsnd)
  REAL, INTENT(IN) :: usnd(nzsnd)
  REAL, INTENT(IN) :: vsnd(nzsnd)
  REAL, INTENT(IN) :: rfrsnd(nzsnd)

  INTEGER, INTENT(IN) :: bkgopt
  INTEGER, INTENT(IN) :: shropt
  INTEGER, INTENT(IN) :: rfropt
  INTEGER, INTENT(IN) :: igatesp
  INTEGER, INTENT(IN) :: irfirst
  INTEGER, INTENT(IN) :: iangle
  INTEGER, INTENT(IN) :: itilt
  INTEGER, INTENT(IN) :: iyear,imon,iday,ihr,imin,isec

  REAL, INTENT(IN) :: radlat
  REAL, INTENT(IN) :: radlon
  REAL, INTENT(IN) :: radalt
  REAL, INTENT(IN) :: rvel(maxgate,maxazim)
  REAL, INTENT(IN) :: spw(maxgate,maxazim)
  REAL, INTENT(IN) :: elev(maxazim)
  REAL, INTENT(IN) :: azim(maxazim)
  REAL, INTENT(IN) :: nqvel(maxazim)

  REAL, INTENT(OUT) :: unfvel(maxgate,maxazim)
  REAL, INTENT(OUT) :: bkgvel(maxgate,maxazim)

  INTEGER, INTENT(OUT) :: bgate(maxazim)
  INTEGER, INTENT(OUT) :: egate(maxazim)
  REAL, INTENT(OUT) :: tmp1(maxgate,maxazim)
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,iray,jray,kray,eray,ipass
  INTEGER :: js1,js2,jstart
  INTEGER :: igate,kgate,kkgate
  INTEGER :: bkgate,ekgate,lgate,kgate1,kgate2
  INTEGER :: k,kh,kclose,kstart,kend
  INTEGER :: knt,knt1,knt2,kntall,kntgood,kntfold,kntrej
  INTEGER :: istatwrt,rem_points,found,nhlfazim
  INTEGER :: igmin,igmax,irayd
  REAL :: bknt,fknt,tknt,tpoint
  REAL :: twonyq,inv2nyq,thrpri,thrpr2
  REAL :: dbkvel,sum,vavg,rvwgt,tstdev,tstvel,sum1,vavg1,range
  REAL :: uvel,vvel,bkrvel
  REAL :: hlfbeam,elvtop,elvbot,hgtagl,hgttop,hgtbot,hgtmsl
  REAL :: refvel,elevavg,sfcr
  REAL :: dtr,umean,vmean,veloc,xcomp,ycomp,dotp,dotmin
  REAL :: diff,shdiff,shdif2,shdif3,azdiff,azd,azd1,azd2,azd3,azim2
  REAL :: rngfrst,gatespc
  REAL :: hgtmin,hgtmax
  REAL :: vnyq
  LOGICAL :: unfdiag
!
!-----------------------------------------------------------------------
!
! Constants
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: rvchek = -200.
  REAL, PARAMETER :: spwflag = -931.
  REAL, PARAMETER :: beamwid = 1.0
  REAL, PARAMETER :: bkgwgt = 0.33
  REAL, PARAMETER :: eradius = 6371000.
  INTEGER, PARAMETER :: lukbak = 20
  INTEGER, PARAMETER :: lukbak2 = 25
  INTEGER, PARAMETER :: lukfwd2 = 25
  INTEGER, PARAMETER :: lukbakr = 20
  CHARACTER (LEN=6) :: varid
  CHARACTER (LEN=20):: varname
  CHARACTER (LEN=20):: varunits = 'm/s'
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'remapcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  kntall=0
  kntgood=0
  kntfold=0
  kntrej=0
  hlfbeam=0.5*beamwid
  rvwgt=1.0-bkgwgt
  dtr=atan(1.)/45.
  rngfrst=float(irfirst)
  gatespc=float(igatesp)
  unfdiag=(unfdiagopt > 0)
  nhlfazim=(nazim/2)+1
  WRITE(6,'(a,f10.1)') ' unfoldnqc: nqvel(1)=',nqvel(1)
!
!-----------------------------------------------------------------------
!
! Apply Spectrum width thresholding
!
!-----------------------------------------------------------------------
!
  IF(spwthrrat > 0.) THEN
    CALL spwthr(maxgate,maxazim,                                        &
                ngate,nazim,rvchek,spwthrrat,spwflag,nqvel,rvel,spw)
  END IF
!
!-----------------------------------------------------------------------
!
! Set up tmp1 to be a quality array. 1=good 0=bad/missing
!
!-----------------------------------------------------------------------
!
  igmin=ngate/2
  igmax=igmin+1
!
  DO iray=1,maxazim
    bgate(iray)=1
    egate(iray)=1
  END DO
  DO iray=1,maxazim
    DO igate=1,maxgate
      tmp1(igate,iray)=0.
      bkgvel(igate,iray)=0.
    END DO
  END DO
!
  DO iray=1,nazim
    bgate(iray)=0
    DO igate=1,ngate
      unfvel(igate,iray)=rvel(igate,iray)
      IF( rvel(igate,iray) > rvchek) THEN
        tmp1(igate,iray)=1.
        egate(iray)=igate
        IF(bgate(iray) == 0) bgate(iray)=igate
      END IF
    END DO
    bgate(iray)=max(bgate(iray),1)
    igmin=min(igmin,bgate(iray))
    igmax=max(igmax,egate(iray))
  END DO
!
  elevavg=0.
  DO iray=1,nazim
    elevavg=elevavg+elev(iray)
  END DO
  elevavg=elevavg/float(nazim)
  WRITE(6,'(a,f10.2)') ' Average Elevation Angle: ',elevavg
!
!-----------------------------------------------------------------------
!
! Determine mean environmental wind in layer observed by radar.
! Used to determine starting radial for unfolding.
!
!-----------------------------------------------------------------------
!
  range=rngfrst+(igmin-1)*gatespc
  CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,elevavg,range,        &
                hgtmin,sfcr)
  hgtmin=hgtmin+radalt
!
  range=rngfrst+(igmax-1)*gatespc
  CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,elevavg,range,        &
                hgtmax,sfcr)
  hgtmax=hgtmax+radalt
!
  write(6,'(a,f12.2,f12.2)') ' Height range (MSL): ',hgtmin,hgtmax
!
  knt=0
  umean=0.
  vmean=0.
  DO k=1,nzsnd
    IF(zsnd(k) > hgtmax ) EXIT
    IF(zsnd(k) > hgtmin ) THEN
      knt=knt+1
      umean=umean+usnd(k)
      vmean=vmean+vsnd(k)
    END IF
  END DO
  knt=max(knt,1)
  umean=umean/float(knt)
  vmean=vmean/float(knt)
  write(6,'(a,f12.2,a,f12.2)')                                          &
   ' Mean velocity, umean: ',umean,' vmean: ',vmean

  veloc=sqrt(umean*umean+vmean*vmean)
  IF(veloc > 0.0) THEN
    umean=umean/veloc
    vmean=vmean/veloc
  END IF
!
!-----------------------------------------------------------------------
!
! Find observed radial with direction most perpendicular to the mean wind.
! Smallest dot product between mean wind and azimuth vector.
!
!-----------------------------------------------------------------------
!
  dotmin=99999.
  js1=1
  DO iray=1,nazim
    xcomp=sin(dtr*azim(iray))
    ycomp=cos(dtr*azim(iray))
    dotp=abs(umean*xcomp+vmean*ycomp)
    IF( dotp < dotmin) THEN
      dotmin=dotp
      js1=iray
    END IF
  END DO
!
  write(6,'(a,i6,f12.2)')                                               &
    ' Most perpendicular radial: ',js1,azim(js1)
!
! Next find the radial opposite of this (+/-180)
!
  azim2=azim(js1)+180.
  IF(azim2 > 360.) azim2=azim2-360.
  azdiff=180.
  js2=js1
  DO iray=1,nazim
    azd1=abs(azim(iray)-azim2)
    azd2=abs(azim(iray)-azim2+360.)
    azd3=abs(azim(iray)-azim2-360.)
    azd=min(azd1,azd2,azd3)
    IF( azd < azdiff ) THEN
      azdiff=azd
      js2=iray
    END IF
  END DO
!
  write(6,'(a,i6,f12.2)')                                               &
    ' Opposite radial: ',js2,azim(js2)
!
! Find which zone around js1 or js2 has the most valid data.
! This is the best place to start
!
  knt1=0
  knt2=0
  DO iray=1,nazim
    irayd=abs(iray-js1)
    IF ( irayd < 6 ) THEN
      DO igate=bgate(iray),egate(iray)
        IF(tmp1(igate,iray) > 0.) knt1=knt1+1
      END DO
    END IF
    irayd=abs(iray-js2)
    IF ( irayd < 6 ) THEN
      DO igate=bgate(iray),egate(iray)
        IF(tmp1(igate,iray) > 0.) knt2=knt2+1
      END DO
    END IF
  END DO
!
  write(6,'(a,i6,a,i6)')                                                &
    ' Data count js1: ',knt1,'  Data count js2: ',knt2

  IF( knt2 > knt1 ) THEN
    jstart=js2
  ELSE
    jstart=js1
  END IF
!
!-----------------------------------------------------------------------
!
! Determine background radial velocity at each valid data point.
!
!-----------------------------------------------------------------------
!
  DO iray=1,nazim
    DO igate=bgate(iray),egate(iray)
      IF( tmp1(igate,iray) > 0. ) THEN
        range=rngfrst+(igate-1)*gatespc
        CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,                &
                      elev(iray),range,hgtagl,sfcr)
        hgtmsl=hgtagl+radalt
        elvtop=elev(iray)+hlfbeam
        CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,                &
                      elvtop,range,hgtagl,sfcr)
        hgttop=hgtagl+radalt
        hgttop=max(hgttop,(hgtmsl+200.))
        elvbot=elev(iray)-hlfbeam
        CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,                &
                      elvbot,range,hgtagl,sfcr)
        hgtbot=hgtagl+radalt
        hgtbot=min(hgtbot,(hgtmsl-200.))
        knt=0
        uvel=0.
        vvel=0.
        DO k=1,nzsnd
          IF( zsnd(k) > hgtbot .AND. zsnd(k) < hgttop ) THEN
            knt=knt+1
            uvel=uvel+usnd(k)
            vvel=vvel+vsnd(k)
          END IF
          IF (zsnd(k) > hgttop) EXIT
        END DO
        IF(knt > 0) THEN
          uvel=uvel/float(knt)
          vvel=vvel/float(knt)
          CALL uv2vr(range,elev(iray),azim(iray),                       &
                     radlat,radlon,radalt,                              &
                     uvel,vvel,bkgvel(igate,iray))
        ELSE
          diff = 99999.0
          DO k=1,nzsnd
            IF( abs(zsnd(k)-hgtmsl) < diff ) THEN
               diff = abs(zsnd(k)-hgtmsl)
               kclose = k
            END IF
          END DO
          CALL uv2vr(range,elev(iray),azim(iray),                       &
                     radlat,radlon,radalt,                              &
                     usnd(kclose),vsnd(kclose),bkgvel(igate,iray))
        END IF
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! Checking of data vs background radial velocity
!
!-----------------------------------------------------------------------
!
  IF( bkgopt > 0 ) THEN
    write(6,'(a,i4)')                                                   &
      ' Unfolding using background wind profile, bkgopt=',bkgopt
    DO iray=1,nazim
      !=========================================================
      ! Beginning of the block added by Youngsun Jung (1/28/09)
      !=========================================================
      vnyq = nqvel(iray)
      twonyq = 2.0*vnyq
      inv2nyq = 1./twonyq
      thrpri = min(max((0.5*vnyq),10.),vnyq)
      thrpr2 = min(vnyq,(1.5*thrpri))
      !=========================================================
      ! End of the block added by Youngsun Jung
      !=========================================================
      DO igate=bgate(iray),egate(iray)
        IF( tmp1(igate,iray) > 0. ) THEN
          kntall=kntall+1
          tstdev=twonyq*                                                &
                 NINT((bkgvel(igate,iray)-rvel(igate,iray))*inv2nyq)
          tstvel=rvel(igate,iray)+tstdev
          IF(abs(tstdev) > 0.) THEN
            kntfold=kntfold+1
          ELSE
            kntgood = kntgood+1
          END IF
          unfvel(igate,iray)=tstvel
        END IF
      END DO
    END DO
    write(6,'(/a/,a/,3i12)')                                            &
          '       Comparison to background wind profile:',              &
          '       Gates  Good-Gates     Folded',                        &
                  kntall,kntgood,kntfold
  END IF  ! bkgopt > 0

  IF(unfdiag) THEN
    varid='rdvelc'
    WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
    CALL wrtvel(iangle,itilt,varid,                        &
                iyear,imon,iday,ihr,imin,isec,             &
                igatesp,irfirst,vnyq,                      &
                radname,radlat,radlon,radalt,              &
                maxgate,maxazim,ngate,nazim,               &
                azim,elev,unfvel)
  END IF
!
!-----------------------------------------------------------------------
!
!   Apply gate-to-gate shear checking for velocity folding.
!   Follows Eilts and Smith, 1989, JTech, 118.
!   Loop forward from jstart
!
!-----------------------------------------------------------------------
!
  IF (shropt > 0 ) THEN
    DO ipass = 1,2
    write(6,'(a,i4,a,i3)')                                              &
      ' Unfolding using local shear check, shropt=',shropt,             &
      '  ipass=',ipass
    write(6,'(a,i6,f12.2)')                                             &
      ' Shear-based unfolding will start at: ',jstart,azim(jstart)
    kntall = 0
    kntgood = 0
    kntfold = 0
!
!   First subpass, increasing azimuth index
!
    DO i=1,nazim
      iray=(jstart+i)-1
      IF(iray > nazim) iray=((jstart+i)-1)-nazim
      !=========================================================
      ! Beginning of the block added by Youngsun Jung (1/28/09)
      !=========================================================
      vnyq = nqvel(iray)
      twonyq = 2.0*vnyq
      inv2nyq = 1./twonyq
      thrpri = min(max((0.5*vnyq),10.),vnyq)
      thrpr2 = min(vnyq,(1.5*thrpri))
      !=========================================================
      ! End of the block added by Youngsun Jung
      !=========================================================
      DO igate=bgate(iray),egate(iray)
        IF( tmp1(igate,iray) > 0. ) THEN
          dbkvel=unfvel(igate,iray)-bkgvel(igate,iray)
          kntall=kntall+1
          lgate = 0
          IF(igate > bgate(iray)) THEN
            ekgate = max(bgate(iray),igate-lukbak)
            DO kgate = igate-1,ekgate,-1
              IF(tmp1(kgate,iray) >0. ) THEN
                lgate = kgate
                EXIT
              END IF
            END DO
          END IF
          IF(lgate >0) THEN
            refvel = unfvel(lgate,iray)-bkgvel(lgate,iray)
          ELSE
            refvel = 999.
          END IF
          shdiff=abs(dbkvel-refvel)
          IF(shdiff > thrpri) THEN
            sum=0.
            bknt=0.
            IF(igate > bgate(iray)) THEN
              ekgate=min(max((igate-lukbak2),bgate(iray)),(igate-1))
              DO kgate=igate-1,ekgate,-1
                bknt=bknt+tmp1(kgate,iray)
                sum=sum+tmp1(kgate,iray)*                              &
                    (unfvel(kgate,iray)-bkgvel(kgate,iray))
                IF (bknt > 2. ) EXIT
              END DO
            END IF
            fknt=0.
            eray=max((iray-lukbakr),1)
            DO kray=iray-1,eray,-1
              ekgate=min((igate+lukfwd2),egate(kray))
              DO kgate=igate,ekgate
                fknt=fknt+tmp1(kgate,kray)
                sum=sum+tmp1(kgate,kray)*                             &
                    (unfvel(kgate,kray)-bkgvel(kgate,kray))
                IF ( fknt > 4. ) EXIT
              END DO
            END DO
            tknt=bknt+fknt
            IF( tknt > 0. ) THEN
              vavg=(rvwgt*sum/tknt)+bkgvel(igate,iray)
              shdif2=abs(unfvel(igate,iray)-vavg)
              IF( shdif2 < thrpr2 ) THEN
                kntgood=kntgood+1
              ELSE
                tstdev=twonyq*NINT((vavg-unfvel(igate,iray))*inv2nyq)
                tstvel=unfvel(igate,iray)+tstdev
                shdif3=abs(tstvel-vavg)
                IF( shdif3 < thrpr2 ) THEN
                  IF( abs(tstdev) > 0. ) THEN
                    kntfold=kntfold+1
!                   write(6,'(a,f10.1,a,f10.1,a,f10.1)')            &
!                   ' Unfold: Meas=',unfvel(igate,iray),' Avgvel=', &
!                               vavg,'  New=',tstvel
                    unfvel(igate,iray)=tstvel
                  ELSE
                    kntgood=kntgood+1
                  END IF
                END IF
              END IF               ! unfvel-vavg >thresh
            END IF
          ELSE
            kntgood=kntgood+1
          END IF                !(unfvel(igate,iray) - refvel) > thrpri
        END IF                  ! tmp1(igate,iray) > 0.
      END DO
    END DO
    write(6,'(/a/,a/,3i12)')                                              &
          '       After forward pass shear consistency check',            &
          '       Gates  Passed Check   Folded',                          &
                  kntall,kntgood,kntfold

    IF(unfdiag) THEN
      IF(ipass == 1) THEN
        varid='rdveld'
      ELSE
        varid='rdvelf'
      END IF
      WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
      CALL wrtvel(iangle,itilt,varid,                        &
                  iyear,imon,iday,ihr,imin,isec,             &
                  igatesp,irfirst,vnyq,                      &
                  radname,radlat,radlon,radalt,              &
                  maxgate,maxazim,ngate,nazim,               &
                  azim,elev,unfvel)
    END IF
!
!-----------------------------------------------------------------------
!
!   Apply gate-to-gate shear checking for velocity folding.
!   Loop backward from jstart-1
!
!-----------------------------------------------------------------------
!
    kntall = 0
    kntgood = 0
    kntfold = 0
    kntrej = 0
    DO i=1,nazim
      iray=(jstart-i)+1
      IF(iray <= 0) iray=(jstart-i)+nazim
      !=========================================================
      ! Beginning of the block added by Youngsun Jung (1/28/09)
      !=========================================================
      vnyq = nqvel(iray)
      twonyq = 2.0*vnyq
      inv2nyq = 1./twonyq
      thrpri = min(max((0.5*vnyq),10.),vnyq)
      thrpr2 = min(vnyq,(1.5*thrpri))
      !=========================================================
      ! End of the block added by Youngsun Jung
      !=========================================================
      DO igate=bgate(iray),egate(iray)
        IF( tmp1(igate,iray) > 0. ) THEN
          dbkvel=unfvel(igate,iray)-bkgvel(igate,iray)
          kntall=kntall+1
          lgate = 0
          IF(igate > bgate(iray)) THEN
            ekgate = max((igate-lukbak),bgate(iray))
            DO kgate = igate-1,ekgate,-1
              IF(tmp1(kgate,iray) >0. ) THEN
                lgate = kgate
                EXIT
              END IF
            END DO
          END IF
          IF(lgate > 0) THEN
            refvel = unfvel(lgate,iray)-bkgvel(lgate,iray)
          ELSE
            refvel = 999.
          END IF
          shdiff=abs(dbkvel-refvel)
          IF( shdiff > thrpri) THEN
            sum=0.
            fknt=0.
            bknt=0.
            IF(igate > bgate(iray)) THEN
              ekgate=max((igate-lukbak2),bgate(iray))
              DO kgate=igate-1,ekgate,-1
                bknt=bknt+tmp1(kgate,iray)
                sum=sum+tmp1(kgate,iray)*                                &
                    (unfvel(kgate,iray)-bkgvel(kgate,iray))
                IF (bknt > 2. ) EXIT
              END DO
            END IF
            eray=min((iray+lukbakr),nazim)
            DO kray=iray+1,eray
              ekgate=min((igate+lukfwd2),egate(kray))
              DO kgate=igate,ekgate
                fknt=fknt+tmp1(kgate,kray)
                sum=sum+tmp1(kgate,kray)*                                &
                    (unfvel(kgate,kray)-bkgvel(kgate,kray))
                IF ( fknt > 4. ) EXIT
              END DO
            END DO

            tknt=bknt+fknt
            IF(tknt > 0. ) THEN
              vavg=(rvwgt*sum/tknt)+bkgvel(igate,iray)
              shdif2=abs(unfvel(igate,iray)-vavg)
              IF( shdif2 < thrpr2 ) THEN
                kntgood=kntgood+1
              ELSE
                tstdev=twonyq*NINT((vavg-unfvel(igate,iray))*inv2nyq)
                tstvel=unfvel(igate,iray)+tstdev
                shdif3=(tstvel-vavg)
                IF( shdif3 < thrpr2 ) THEN
                  IF( abs(tstdev) > 0. ) THEN
                    kntfold=kntfold+1
!                   write(6,'(a,f10.1,a,f10.1,a,f10.1)')            &
!                   ' Unfold: Meas=',unfvel(igate,iray),' Avgvel=', &
!                               vavg,'  New=',tstvel
                    unfvel(igate,iray)=tstvel
                  ELSE
                    kntgood=kntgood+1
                  END IF
                ELSE
                  kntrej=kntrej+1
!                 write(6,'(a,i5,a,i5,a,f10.1,a,f10.1,a,f10.1)')  &
!                  'igate=',igate,' iray=',iray,                  &
!                  ' Unresolved: Meas=',unfvel(igate,iray),      &
!                  ' Avgvel=',vavg,'  Try=',tstvel
                  unfvel(igate,iray)=unfvel(igate,iray)-2000.
                  tmp1(igate,iray)=0.
                END IF              ! tstvel-vavg >thresh
              END IF               ! unfvel-vavg >thresh
            END IF           ! tknt > 0
          ELSE
            kntgood=kntgood+1
          END IF                !(unfvel(igate,iray) - refvel) > thrpri
        END IF                  ! tmp1(igate,iray) > 0.
      END DO
    END DO

    write(6,'(/a/,a/,4i12)')                                            &
          '       After reverse pass shear consistency check',          &
          '       Gates  Passed Check   Folded     Rejected',           &
                  kntall,kntgood,kntfold,kntrej
    IF(unfdiag) THEN
      IF(ipass == 1) THEN
        varid='rdvele'
      ELSE
        varid='rdvelg'
      END IF
      WRITE(6,'(a,a6)') ' Calling wrtvel varid:',varid
      CALL wrtvel(iangle,itilt,varid,                        &
                  iyear,imon,iday,ihr,imin,isec,             &
                  igatesp,irfirst,vnyq,                      &
                  radname,radlat,radlon,radalt,              &
                  maxgate,maxazim,ngate,nazim,               &
                  azim,elev,unfvel)
    END IF
  END DO !ipass

  END IF ! shropt > 0

  RETURN
END SUBROUTINE unfoldnqc
