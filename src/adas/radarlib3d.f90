!  program test
!  implicit none
!  integer nzsnd
!  parameter (nzsnd=5)
!  real zsnd(nzsnd)
!  real rfrsnd(nzsnd)
!  real lat1,lon1,lat2,lon2,dist,head
!  real rdralt,elev,range,azim
!  real sfcrng,height,dhdr
!  integer k,rfropt
!  rfropt=3
!  DO k=1,nzsnd
!    zsnd(k)=(k-1)*2000.
!    rfrsnd(k)=300.-20.*(k-1)
!  END DO
!  rdralt=100. 
! 10  CONTINUE
!  print *, ' Enter height (m), sfcrange (km): '
!  read(5,*) height,sfcrng
!  IF(height.lt.-10.) STOP
!  sfcrng=sfcrng*1000.
!  CALL beamelv(height,sfcrng,elev,range)
!  print *, ' elv, range = ',elev,(0.001*range)
!  CALL beamelvn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,height,sfcrng,elev,range)
!  print *, ' elv, rangeN= ',elev,(0.001*range)
!  CALL beamhgt(elev,range,height,sfcrng)
!  print *, ' height,sfcrng = ',height,(0.001*sfcrng)
!  CALL beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,elev,range,height,sfcrng)
!  print *, ' height,sfcrngN= ',height,(0.001*sfcrng)
!
!  CALL bmhgtsfr(elev,sfcrng,height)
!  print *, '       height2 = ',height
!  CALL bmhgtsfrn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,elev,sfcrng,height)
!  print *, '       height2N= ',height
!
!   print *, ' Enter lat1,lon1: '
!   read (5,*) lat1,lon1
!   lat1=35.
!   lon1=-100.
!   print *, ' Enter elev,azim,range '
!   read (5,*) elev,azim,range
!   IF(elev.gt.90.) STOP
!   CALL beamhgt(elev,range,height,sfcrng)
!   print *, ' beam height = ',height
!   print *, ' sfc range   = ',sfcrng
!   CALL dhdrange(elev,range,dhdr)
!   print *, ' local elv   = ',locelva
!   CALL gcircle(lat1,lon1,azim,sfcrng,lat2,lon2)
!   print *, ' gate lat,lon = ',lat2,lon2
!   CALL disthead(lat1,lon1,lat2,lon2,head,dist)
!   print *, ' distance, heading: ',dist,head
!
!  GO TO 10
!  END
!

SUBROUTINE beamhgt(elvang,range,height,sfcrng)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the height of the radar beam and the along-
!  ground distance from the radar as a function
!  distance along the radar beam (range) and radar
!  elevation angle (elvang).
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is from Eq. 2.28 of Doviak and Zrnic', Doppler Radar
!  and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
!
!  OUTPUT:
!    height   Height (meters) of beam above ground.
!    sfcrng   Distance (meters) of point along ground from radar.
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
  REAL, INTENT(IN) :: elvang
  REAL, INTENT(IN) :: range
  REAL, INTENT(OUT) :: height
  REAL, INTENT(OUT) :: sfcrng
!
  DOUBLE PRECISION :: eradius,frthrde,eighthre,fthsq,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             eighthre=(8.*eradius/3.),                                  &
             fthsq=(frthrde*frthrde),                                   &
             deg2rad=(3.14592654/180.))
!
  DOUBLE PRECISION :: elvrad,hgtdb,rngdb,drange
!
  elvrad=deg2rad*DBLE(elvang)
  drange=DBLE(range)
  hgtdb = SQRT(drange*drange + fthsq +                                  &
                eighthre*drange*SIN(elvrad)) -                          &
                frthrde
  height=hgtdb
  rngdb = frthrde *                                                     &
           ASIN (drange*COS(elvrad)/(frthrde + hgtdb) )
  sfcrng=rngdb
  RETURN
END SUBROUTINE beamhgt
!

SUBROUTINE bmhgtsfr(elvang,sfcrng,height)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the height of the radar beam as a function of
!  the elevation angle (elvang) and the along-ground distance from
!  the radar (sfcrng).
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is from Eq. 2.28a of Doviak and Zrnic', Doppler Radar
!  and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  01/14/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    sfcrng   Distance (meters) of point along ground from radar.
!
!  OUTPUT:
!    height   Height (meters) of beam above ground.
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
  REAL, INTENT(IN) :: elvang
  REAL, INTENT(IN) :: sfcrng
  REAL, INTENT(OUT) :: height
!
  DOUBLE PRECISION :: eradius,frthrde,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             deg2rad=(3.14592654/180.))
!
  DOUBLE PRECISION :: elvrad,elvkea,srange,hgtdb
!
  elvrad=deg2rad*DBLE(elvang)
  srange=DBLE(sfcrng)
  elvkea=elvrad+(srange/frthrde)
  hgtdb=frthrde*((cos(elvrad)/cos(elvkea))-1.0)
  height=hgtdb
  RETURN
END SUBROUTINE bmhgtsfr
!
SUBROUTINE beamelv(height,sfcrng,elvang,range)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the elevation angle (elvang) and the along
!  ray-path distance (range) of a radar beam
!  crossing through the given height and along-ground
!  distance.
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is dervied from Eq. 2.28 of Doviak and Zrnic',
!  Doppler Radar and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  10/10/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    height   Height (meters) of beam above ground.
!    sfcrng   Distance (meters) of point along ground from radar.
!
!  OUTPUT
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
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
  REAL, INTENT(IN) :: height
  REAL, INTENT(IN) :: sfcrng
  REAL, INTENT(OUT) :: elvang
  REAL, INTENT(OUT) :: range
!
  DOUBLE PRECISION :: eradius,frthrde,rad2deg
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             rad2deg=(180./3.14592654))
!
  DOUBLE PRECISION :: elvrad,hgtdb,rngdb,drange
!
  IF(sfcrng > 0.) THEN

    hgtdb=frthrde+DBLE(height)
    rngdb=DBLE(sfcrng)/frthrde

    elvrad = ATAN((hgtdb*COS(rngdb) - frthrde)/(hgtdb * SIN(rngdb)))
    drange = (hgtdb*SIN(rngdb))/COS(elvrad)
    elvang=rad2deg*elvrad
    range=drange

  ELSE

    elvang=90.
    range=height

  END IF
  RETURN
END SUBROUTINE beamelv
!

SUBROUTINE dhdrange(elvang,range,dhdr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the local change in height of the radar
!  beam with respect to a change in range.  Due to
!  curvature of the beam and the earth's surface this is
!  generally different what would be calculated from the
!  elevation angle measured at the radar.  This derivative
!  is needed for finding 3-d velocities from radial winds
!  and accounting for terminal velocity of precipitation.
!
!  This formulation, consistent with subroutine beamhgt,
!  assumes a 4/3 earth radius beam curvature.  This formula
!  is obtained by differentiating Eq 2.28 of Doviak and
!  Zrnic', Doppler Radar and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
!
!  OUTPUT:
!    dhdr     Change in height per change in range (non-dimensional)
!
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
  REAL, INTENT(IN) :: range
  REAL, INTENT(IN) :: elvang
  REAL, INTENT(OUT) :: dhdr
!
  DOUBLE PRECISION :: eradius,frthrde,eighthre,fthsq,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             eighthre=(8.*eradius/3.),                                  &
             fthsq=(frthrde*frthrde),                                   &
             deg2rad=(3.14592654/180.))
!
  DOUBLE PRECISION :: sinelv,dhdrdb,drange
!
  drange=DBLE(range)
  sinelv=SIN(deg2rad*DBLE(elvang))
  dhdrdb = (drange+frthrde*sinelv)/                                     &
         SQRT(drange*drange + fthsq + eighthre*drange*sinelv)
  dhdr = dhdrdb
!
  RETURN
END SUBROUTINE dhdrange

SUBROUTINE beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,elvang,range,        &
                    height,sfcrng)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the height of the radar beam and the along-
!  ground distance from the radar as a function
!  distance along the radar beam (range) and radar
!  elevation angle (elvang).
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is from Eq. 2.28 of Doviak and Zrnic', Doppler Radar
!  and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
!    rke      Equivalent earth radius (from dn/dz)
!
!  OUTPUT:
!    height   Height (meters) of beam above ground.
!    sfcrng   Distance (meters) of point along ground from radar.
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
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN) :: zsnd(nzsnd)
  REAL, INTENT(IN) :: rfrsnd(nzsnd)
  REAL, INTENT(IN) :: rdralt
  INTEGER, INTENT(IN) :: rfropt
  REAL, INTENT(IN) :: elvang
  REAL, INTENT(IN) :: range
  REAL, INTENT(OUT) :: height
  REAL, INTENT(OUT) :: sfcrng
!
  DOUBLE PRECISION :: eradius,frthrde,eighthre,fthsq,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             eighthre=(8.*eradius/3.),                                  &
             fthsq=(frthrde*frthrde),                                   &
             deg2rad=(3.14592654/180.))
  REAL :: dndhlim
  PARAMETER (dndhlim=-1./(eradius+10000.))
!
! Misc. Local Variables
!
  DOUBLE PRECISION :: rkear,elvrad,hgtdb,rngdb,drange
  REAL :: delz,hgt,rfrrad,rfrhgt,dndh,rke
  INTEGER :: kh
!
  elvrad=deg2rad*DBLE(elvang)
  drange=DBLE(range)
  hgtdb = SQRT(drange*drange + fthsq +                                &
                eighthre*drange*SIN(elvrad)) -                          &
                frthrde
  height=hgtdb
  rngdb = frthrde *                                                   &
          ASIN (drange*COS(elvrad)/(frthrde + hgtdb) )
  sfcrng=rngdb

  IF( rfropt > 1 .AND. height > 0. ) THEN
!
    delz=zsnd(2)-zsnd(1)
    kh=max((int((rdralt-zsnd(1))/delz)+1),1)
    kh=min(kh,(nzsnd-1))
    rfrrad=rfrsnd(kh)+((rdralt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    hgt=height+rdralt
    kh=max((int((hgt-zsnd(1))/delz)+1),1)
    rfrhgt=rfrsnd(kh)+((hgt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    dndh=1.0E-06*(rfrhgt-rfrrad)/height
    dndh=max(dndh,dndhlim)
    rke=1.0/(1.0+eradius*dndh)
!
    rkear=rke*eradius
    hgtdb = SQRT(drange*drange + rkear*rkear + &
                 2.*rkear*drange*SIN(elvrad)) - rkear
    height=hgtdb
    rngdb = rkear * ASIN (drange*COS(elvrad)/(rkear + hgtdb) )
    sfcrng=rngdb
!
  END IF

  IF( rfropt > 2 .AND. height > 0. ) THEN
!
    hgt=height+rdralt
    kh=max((int((hgt-zsnd(1))/delz)+1),1)
    kh=min(kh,(nzsnd-1))
    rfrhgt=rfrsnd(kh)+((hgt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    dndh=1.0E-06*(rfrhgt-rfrrad)/height
    dndh=max(dndh,dndhlim)
    rke=1.0/(1.0+eradius*dndh)
!
    rkear=rke*eradius
    hgtdb = SQRT(drange*drange + rkear*rkear  + &
                 2.*rkear*drange*SIN(elvrad)) - rkear
    height=hgtdb
    rngdb = rkear * ASIN (drange*COS(elvrad)/(rkear + hgtdb) )
    sfcrng=rngdb
!
  END IF
!
  RETURN
END SUBROUTINE beamhgtn
!

SUBROUTINE bmhgtsfrn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,    &
                     elvang,sfcrng,height)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the height of the radar beam as a function of
!  the elevation angle (elvang) and the along-ground distance from
!  the radar (sfcrng).
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is from Eq. 2.28a of Doviak and Zrnic', Doppler Radar
!  and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  01/14/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    sfcrng   Distance (meters) of point along ground from radar.
!
!  OUTPUT:
!    height   Height (meters) of beam above ground.
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
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN) :: zsnd(nzsnd)
  REAL, INTENT(IN) :: rfrsnd(nzsnd)
  REAL, INTENT(IN) :: rdralt
  INTEGER, INTENT(IN) :: rfropt
  REAL, INTENT(IN) :: elvang
  REAL, INTENT(IN) :: sfcrng
  REAL, INTENT(OUT) :: height
!
  DOUBLE PRECISION :: eradius,frthrde,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             deg2rad=(3.14592654/180.))
!
  INTEGER :: kh
  REAL :: delz,hgt,dndh,rfrrad,rfrhgt,rke
  DOUBLE PRECISION :: rkear,elvrad,elvkea,srange,hgtdb
  REAL :: dndhlim
  PARAMETER (dndhlim=-1./(eradius+10000.))
!
  elvrad=deg2rad*DBLE(elvang)
  srange=DBLE(sfcrng)
  elvkea=elvrad+(srange/frthrde)
  hgtdb=frthrde*((cos(elvrad)/cos(elvkea))-1.0)
  height=hgtdb
  IF ( rfropt > 1 .AND. height > 0. ) THEN
!
    delz=zsnd(2)-zsnd(1)
    kh=max((int((rdralt-zsnd(1))/delz)+1),1)
    kh=min(kh,(nzsnd-1))
    rfrrad=rfrsnd(kh)+((rdralt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    hgt=height+rdralt
    kh=max((int((hgt-zsnd(1))/delz)+1),1)
    rfrhgt=rfrsnd(kh)+((hgt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    dndh=1.0E-06*(rfrhgt-rfrrad)/height
    dndh=max(dndh,dndhlim)
    rke=1.0/(1.0+eradius*dndh)
!
    rkear=rke*eradius
    elvkea=elvrad+(srange/rkear)
    hgtdb=rkear*((cos(elvrad)/cos(elvkea))-1.0)
    height=hgtdb
  END IF

  IF ( rfropt > 2 .AND. height > 0. ) THEN
    hgt=height+rdralt
    kh=max((int((hgt-zsnd(1))/delz)+1),1)
    rfrhgt=rfrsnd(kh)+((hgt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    dndh=1.0E-06*(rfrhgt-rfrrad)/height
    dndh=max(dndh,dndhlim)
    rke=1.0/(1.0+eradius*dndh)
!
    rkear=rke*eradius
    elvkea=elvrad+(srange/rkear)
    hgtdb=rkear*((cos(elvrad)/cos(elvkea))-1.0)
    height=hgtdb
  END IF

  RETURN
END SUBROUTINE bmhgtsfrn
!
SUBROUTINE beamelvn(nzsnd,zsnd,rfrsnd,rdralt,rfropt, &
                    height,sfcrng,elvang,range)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the elevation angle (elvang) and the along
!  ray-path distance (range) of a radar beam
!  crossing through the given height and along-ground
!  distance.
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is dervied from Eq. 2.28 of Doviak and Zrnic',
!  Doppler Radar and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  10/10/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    height   Height (meters) of beam above ground.
!    sfcrng   Distance (meters) of point along ground from radar.
!    rke      Equivalent earth radius (from dn/dz)
!
!  OUTPUT
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
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
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)
  REAL, INTENT(IN)    :: rdralt
  INTEGER, INTENT(IN) :: rfropt
  REAL, INTENT(IN)    :: height
  REAL, INTENT(IN)    :: sfcrng
  REAL, INTENT(OUT)   :: elvang
  REAL, INTENT(OUT)   :: range
!
  DOUBLE PRECISION :: eradius,frthrde,rad2deg
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             rad2deg=(180./3.14592654))
!
  DOUBLE PRECISION :: rkear,elvrad,hgtdb,rngdb,drange
!
  REAL :: dndhlim
  PARAMETER (dndhlim=-1./(eradius+10000.))
!
! Misc Local Variables
!
  INTEGER :: kh
  REAL :: delz,rfrrad,rfrhgt,hgt,dndh,rke
!
  IF(sfcrng > 0.) THEN

    IF ( rfropt == 1 ) THEN
      hgtdb=frthrde+DBLE(height)
      rngdb=DBLE(sfcrng)/frthrde

      elvrad = ATAN((hgtdb*COS(rngdb) - frthrde)/(hgtdb * SIN(rngdb)))
      drange = (hgtdb*SIN(rngdb))/COS(elvrad)
      elvang=rad2deg*elvrad
      range=drange
    ELSE
      delz=zsnd(2)-zsnd(1)
      kh=max((int((rdralt-zsnd(1))/delz)+1),1)
      kh=min(kh,(nzsnd-1))
      rfrrad=rfrsnd(kh)+((rdralt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
      hgt=height+rdralt
      kh=max((int((hgt-zsnd(1))/delz)+1),1)
      rfrhgt=rfrsnd(kh)+((hgt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
      dndh=1.0E-06*(rfrhgt-rfrrad)/height
      dndh=max(dndh,dndhlim)
      rke=1.0/(1.0+eradius*dndh)
   
      rkear=rke*eradius
      hgtdb=rkear+DBLE(height)
      rngdb=DBLE(sfcrng)/rkear

      elvrad = ATAN((hgtdb*COS(rngdb) - rkear)/(hgtdb * SIN(rngdb)))
      drange = (hgtdb*SIN(rngdb))/COS(elvrad)
      elvang=rad2deg*elvrad
      range=drange
    END IF

  ELSE

    elvang=90.
    range=height

  END IF
  RETURN
END SUBROUTINE beamelvn
!

SUBROUTINE dhdrangn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,  &
                     elvang,range,rke,dhdr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the local change in height of the radar
!  beam with respect to a change in range.  Due to
!  curvature of the beam and the earth's surface this is
!  generally different what would be calculated from the
!  elevation angle measured at the radar.  This derivative
!  is needed for finding 3-d velocities from radial winds
!  and accounting for terminal velocity of precipitation.
!
!  This formulation, consistent with subroutine beamhgt,
!  assumes a 4/3 earth radius beam curvature.  This formula
!  is obtained by differentiating Eq 2.28 of Doviak and
!  Zrnic', Doppler Radar and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
!    rke      Equivalent earth radius (from dn/dz)
!
!  OUTPUT:
!    dhdr     Change in height per change in range (non-dimensional)
!
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
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)
  REAL, INTENT(IN)    :: rdralt
  INTEGER, INTENT(IN) :: rfropt
  REAL, INTENT(IN)    :: range
  REAL, INTENT(IN)    :: elvang
  REAL, INTENT(OUT)   :: dhdr
!
  DOUBLE PRECISION :: eradius,frthrde,eighthre,fthsq,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             eighthre=(8.*eradius/3.),                                  &
             fthsq=(frthrde*frthrde),                                   &
             deg2rad=(3.14592654/180.))
  REAL :: dndhlim
  PARAMETER (dndhlim=-1./(eradius+10000.))
!
  INTEGER :: kh
  REAL :: delz,rfrrad,hgt,rfrhgt,dndh,rke,height,srange
  DOUBLE PRECISION :: rkear,sinelv,dhdrdb,drange
!
  IF( rfropt == 1 ) THEN
    drange=DBLE(range)
    sinelv=SIN(deg2rad*DBLE(elvang))
    dhdrdb = (drange+frthrde*sinelv)/                                   &
         SQRT(drange*drange + fthsq + eighthre*drange*sinelv)
    dhdr = dhdrdb
  ELSE
    CALL beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,elvang,range,          &
                  height,srange)
    delz=zsnd(2)-zsnd(1)
    kh=max((int((rdralt-zsnd(1))/delz)+1),1)
    kh=min(kh,(nzsnd-1))
    rfrrad=rfrsnd(kh)+((rdralt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    hgt=height+rdralt
    kh=max((int((hgt-zsnd(1))/delz)+1),1)
    rfrhgt=rfrsnd(kh)+((hgt-zsnd(kh))*(rfrsnd(kh+1)-rfrsnd(kh))/delz)
    dndh=1.0E-06*(rfrhgt-rfrrad)/height
    dndh=max(dndh,dndhlim)
    rke=1.0/(1.0+eradius*dndh)
   
    rkear=rke*eradius
    drange=DBLE(range)
    sinelv=SIN(deg2rad*DBLE(elvang))
    dhdrdb = (drange+rkear*sinelv)/                                     &
         SQRT(drange*drange + rkear*rkear + 2.0*rkear*drange*sinelv)
    dhdr = dhdrdb
  END IF
!
  RETURN
END SUBROUTINE dhdrangn

SUBROUTINE disthead(lat1,lon1,lat2,lon2,headng,dist)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Given a pair of locations specified in lat,lon on the earth's
!  surface find the distance between them and the great circle
!  heading from the first point to the second point.  Spherical
!  geometry is used, which is more than adequate for radar
!  and local modelling applications.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    lat1   Latitude (degrees, north positive) of first point
!    lon1   Latitude (degrees, east positive) of first point
!    lat2   Latitude (degrees, north positive) of second point
!    lon2   Latitude (degrees, east positive) of second point
!
!  OUTPUT:
!
!    headng Heading (degrees, north zero) of great circle path
!           at first point.
!    dist   Distance (meters) between two points along great circle
!           great circle arc.
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
!  Arguments
!
  REAL :: lat1,lon1
  REAL :: lat2,lon2
  REAL :: headng
  REAL :: dist
!
!  Parameters
!
  DOUBLE PRECISION :: pi,deg2rad,rad2deg,eradius,one,mone
  PARAMETER (pi=3.141592654,                                            &
             deg2rad=(pi/180.),                                         &
             rad2deg=(180./pi),                                         &
             eradius=6371000.,      & ! Earth radius in meters
         one=1.,                                                        &
             mone=-1.)
!
!  Misc internal variables
!
  DOUBLE PRECISION :: alat1,alat2,dlon,arcdst,cosdst,coshd,denom
!
!  Find arc length using law of cosines
!
!  cos a = cos b cos c + sin b sin c cos A
!  cos (1 to 2) = sin(lat1) * sin (lat2)
!               +(cos(lat1) * sin (lat2)
!                 * cos (lon1 - lon2)
!
  alat1=deg2rad * DBLE(lat1)
  alat2=deg2rad * DBLE(lat2)
  dlon=deg2rad*DBLE(lon2-lon1)
  cosdst = SIN(alat1) * SIN(alat2)  +                                   &
           COS(alat1) * COS(alat2) * COS(dlon)
  arcdst = ACOS(cosdst)
  dist = eradius*arcdst
!
  denom=COS(alat1)*SIN(arcdst)
  headng=0.
  IF(ABS(denom) > 1.e-06) THEN
    coshd=(SIN(alat2) - SIN(alat1)*cosdst) / denom
    coshd=DMAX1(coshd,mone)
    coshd=DMIN1(coshd,one)
    headng=rad2deg*ACOS(coshd)
    IF( SIN(dlon) < 0 ) headng = 360.-headng
  ELSE IF( ABS(COS(alat1)) < 1.e-06 .AND. alat1 > 0.) THEN
    headng=180.
  END IF
!
  RETURN
END SUBROUTINE disthead
!

SUBROUTINE gcircle(lat1,lon1,head,dist,lat2,lon2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Following a great circle path from a point specified as
!  lat,lon on the earth's surface leaving at heading given
!  by head for a distance given by dist, give the location
!  of the end point.  Useful for finding the lat,lon of a
!  radar gate.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    lat1   Latitude (degrees, north positive) of first point
!    lon1   Latitude (degrees, east positive) of first point
!    head   Heading (degrees, north zero) of great circle path
!           at first point.
!    dist   Distance (meters) between two points along great circle
!           great circle arc.
!
!  OUTPUT:
!
!    lat2   Latitude (degrees, north positive) of second point
!    lon2   Latitude (degrees, east positive) of second point
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
!  Arguments
!
  REAL :: lat1,lon1
  REAL :: head
  REAL :: dist
  REAL :: lat2,lon2
!
!  Parameters
!
  DOUBLE PRECISION :: pi,deg2rad,rad2deg,eradius,one,mone
  PARAMETER (pi=3.141592654,                                            &
             deg2rad=(pi/180.),                                         &
             rad2deg=(180./pi),                                         &
             eradius=6371000.,      & ! Earth radius in meters
         one=1.,                                                        &
             mone=-1.)
!
!  Misc internal variables
!
  DOUBLE PRECISION :: alat1,alat2,dlon,arcdst,cosdst,coshd
  DOUBLE PRECISION :: denom,sinlat2,cosdlon
!
  alat1=deg2rad*DBLE(lat1)
  arcdst=DBLE(dist)/eradius
  cosdst=COS(arcdst)
  coshd=COS(deg2rad*DBLE(head))
!
  sinlat2=coshd*COS(alat1)*SIN(arcdst) + SIN(alat1)*cosdst
  sinlat2=DMAX1(sinlat2,mone)
  sinlat2=DMIN1(sinlat2,one)
  alat2=ASIN(sinlat2)
  lat2=rad2deg*alat2
!
  denom=COS(alat1)*COS(alat2)
  IF(denom /= 0.) THEN
    cosdlon=(cosdst - SIN(alat1)*sinlat2)/(COS(alat1)*COS(alat2))
    cosdlon=DMAX1(cosdlon,mone)
    cosdlon=DMIN1(cosdlon,one)
    dlon=rad2deg*ACOS(cosdlon)
    IF(SIN(deg2rad*head) < 0.) dlon=-dlon
    lon2=lon1+dlon
  ELSE
    lon2=lon1
  END IF
  RETURN
END SUBROUTINE gcircle
!

SUBROUTINE rmvterm(nvar_radin,mx_rad,nz_rdr,mx_colrad,                  &
           latrad,lonrad,elvrad,                                        &
           latradc,lonradc,irad,nlevrad,hgtradc,obsrad,                 &
           ncolrad,istatus)
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Removes the precipitation terminal velocity from
!  the observed radial velocity.  From Zeigler, 1978.
!
!  An empirical formula based on the reflectivity is
!  used to estimate the terminal velocity.
!
!  For simplicity, atmospheric density is approxmated
!  by a standard atmosphere, rho(z)=rhonot exp(-z/h0)
!
!-----------------------------------------------------------------------
!
  INTEGER :: mx_rad,nz_rdr,mx_colrad,nvar_radin
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  REAL :: latrad(mx_rad),lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: irad(mx_colrad)
  INTEGER :: nlevrad(mx_colrad)
  REAL :: latradc(mx_colrad),lonradc(mx_colrad)
  REAL :: hgtradc(nz_rdr,mx_colrad)
  REAL :: obsrad(nvar_radin,nz_rdr,mx_colrad)
  INTEGER :: ncolrad
  INTEGER :: istatus
!
!  Parameters
!
  REAL :: deg2rad
  PARAMETER (deg2rad=(3.14592654/180.))
!
  REAL :: zfrez,zice,rho0,h0,denom,dbzmin,dbzmax,vrmax,refmax
  PARAMETER (zfrez=3000.,                                               &
             zice=8000.,                                                &
             rho0=1.2250,                                               &
             h0=7000.,                                                  &
             denom=(1./(zice-zfrez)),                                   &
             dbzmin=10.,                                                &
             dbzmax=100.,                                               &
             vrmax=80.,                                                 &
             refmax=80.)
!
  INTEGER :: icol,klev
  REAL :: azm,rng,refz,rhofact,s1,s2,vt,dz,eleva,bmrng,dhdr
!
  IF(dbzmin > 30.) PRINT *, ' Warning Terminal Velocity Removal Disabled'
!
  DO icol=1,ncolrad
    IF(irad(icol) > 0) THEN
      CALL disthead(latrad(irad(icol)),lonrad(irad(icol)),              &
                  latradc(icol),lonradc(icol),azm,rng)
      DO klev=1,nlevrad(icol)
        IF(obsrad(1,klev,icol) > dbzmin .AND.                           &
              ABS(obsrad(1,klev,icol)) < dbzmax .AND.                   &
              ABS(obsrad(2,klev,icol)) < vrmax) THEN
          refz=10.**(0.1*obsrad(1,klev,icol))
          rhofact=EXP(0.4*hgtradc(klev,icol)/h0)
          IF(hgtradc(klev,icol) < zfrez) THEN
            vt=2.6*(refz**0.107)*rhofact
          ELSE IF(hgtradc(klev,icol) < zice) THEN
            s1=(zice-hgtradc(klev,icol))*denom
            s2=2.*(hgtradc(klev,icol)-zfrez)*denom
            vt=s1*2.6*(refz**0.107)*rhofact + s2
          ELSE
            vt=2.0
          END IF
          dz=hgtradc(klev,icol)-elvrad(irad(icol))
          CALL beamelv(dz,rng,eleva,bmrng)
          CALL dhdrange(eleva,bmrng,dhdr)
          obsrad(2,klev,icol)=obsrad(2,klev,icol) +                     &
                            vt*dhdr
        END IF
      END DO
    END IF
  END DO
  istatus=1
  RETURN
END SUBROUTINE rmvterm

SUBROUTINE raytest(nzsnd,zsnd,rfrsnd,rdralt)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN) :: zsnd(nzsnd)
  REAL, INTENT(IN) :: rfrsnd(nzsnd)
  REAL, INTENT(IN) :: rdralt
!
  REAL :: elev(3)
  REAL :: hgt1(3)
  REAL :: hgt2(3)
  REAL :: hgt3(3)
  REAL :: range,sfcrng
  INTEGER i,k
!
  elev(1)=0.5
  elev(2)=1.5
  elev(3)=2.5
  open(41,file='raypath.txt',status='unknown')
  write(41,'(a)') &
  ' 0.5a    0.5b   0.5c   1.5a   1.5b   1.5c   2.5a   2.5b   2.5c'
  DO i=1,115
    range=i*2000.
    DO k=1,3
      CALL beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,1,                    &
                    elev(k),range,hgt1(k),sfcrng)
      CALL beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,2,                    &
                    elev(k),range,hgt2(k),sfcrng)
      CALL beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,3,                    &
                    elev(k),range,hgt3(k),sfcrng)
    END DO
    WRITE(41,'(9f10.0)') hgt1(1),hgt2(1),hgt3(1),    &
                         hgt1(2),hgt2(2),hgt3(2),    &
                         hgt1(3),hgt2(3),hgt3(3)
  END DO
  close(41)
  RETURN
END SUBROUTINE raytest
