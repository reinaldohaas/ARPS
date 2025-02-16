!
!##################################################################
!##################################################################
!######                                                      ######
!######            ARPS Map Projection Subsystem.            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!   General Information
!
!   This set of subroutines allows for transformation between
!   lat-lon coordinates and any one of three map projections: Polar
!   Stereographic, Lambert Conformal or Mercator.
!
!   In order for the transformation subroutines to work, the
!   map projection must first be set up by calling setmapr.  The
!   user may wish to call setorig immediately after setmapr to
!   established an origin (given a lat-long or x-y in the default
!   system) other than the default origin (e.g., the north pole).
!
!   All lat-lons are in degrees (positive north, negative south,
!   positive east and negative west).  Note carefully the dimensions
!   of x,y -- it differs among the subroutines to conform to ARPS usage.
!   x,y coordinates are meters on earth but may be changed using the scale
!   parameter in setmapr to change to km (scale=0.001) or to a different
!   sphere (e.g., scale=mars_radius/earth_radius).
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SETMAPR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setmapr(iproj,scale,latnot,orient)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set constants for map projections, which are stored in
!  the common block named /projcst/.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/13/93.
!
!  MODIFICATION HISTORY:
!  03/30/1995 (K. Brewster)
!  Corrected error in Lambert Conformal scaling and added code to
!  allow Lambert Tangent projection (lat1=lat2 in Lambert Conformal).
!  Resulted in redefinition of projc1 for option 2.
!
!  2003-12-23 (Richard Carpenter)
!  Fixed a potential source of round-off error in Lambert Tangent
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    iproj        Map projection number
!                 1=North Polar Stereographic   (-1 South Pole)
!                 2=Northern Lambert Conformal  (-2 Southern)
!                 3=Mercator
!                 4=Lat,Lon
!
!    scale        Map scale factor,  at latitude=latnot
!                 Distance on map = (Distance on earth) * scale
!                 For ARPS model runs, generally this is 1.0
!                 For ARPS plotting this will depend on window
!                 size and the area to be plotted.
!
!    latnot(2)    Real "True" latitude(s) of map projection
!                 (degrees, positive north)
!                 Except for iproj=1, only latnot(1) is used
!
!    orient       Longitude line that runs vertically on the map.
!                 (degrees, negative west, positive east)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'          ! Message passing parameters.
  INTEGER :: iproj
  REAL :: scale                       ! map scale factor
  REAL :: latnot(2)                   ! true latitude (degrees N)
  REAL :: orient                      ! orientation longitude (degrees E)

  REAL, PARAMETER :: d2rad   = 3.141592654/180.,                        &
                     eradius = 6371000.   ! mean earth radius in m

  INTEGER :: jproj,jpole
  REAL    :: trulat(2),rota,scmap,xorig,yorig,                          &
             projc1,projc2,projc3,projc4,projc5

  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: denom1,denom2,denom3

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xorig=0.
  yorig=0.
  jproj=IABS(iproj)
  jpole=ISIGN(1,iproj)
!  print *, ' jpole = ',jpole
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF ( jproj == 0 ) THEN
    IF (myproc == 0) WRITE(6,'(a)') '  No map projection will be used.'
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    trulat(1)=latnot(1)
    rota=orient
    scmap=scale
    projc1=scale*eradius
    projc2=(1. + SIN(d2rad*jpole*trulat(1)) )
    projc3=projc1*projc2
    IF(jpole > 0) THEN
      IF (myproc == 0) WRITE(6,'(a/,a)')  &
          '  Map projection set to Polar Stereographic',                &
          '  X origin, Y origin set to 0.,0. at the North Pole.'
    ELSE
      IF (myproc == 0) WRITE(6,'(a/,a)')  &
          '  Map projection set to Polar Stereographic',                &
          '  X origin, Y origin set to 0.,0. at the South Pole.'
    END IF
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
    trulat(1)=latnot(1)
    trulat(2)=latnot(2)
    rota=orient
    scmap=scale
    projc2=COS(d2rad*trulat(1))
    projc3=TAN(d2rad*(45.-0.5*jpole*trulat(1)))
    denom1=COS(d2rad*trulat(2))
    denom2=TAN(d2rad*(45.-0.5*jpole*trulat(2)))
    !IF(denom2 /= 0.) THEN
    IF (ABS(trulat(1)-trulat(2)) > 0.01 .AND. denom2 /= 0.) THEN
      denom3=ALOG( projc3/denom2 )
    ELSE
      denom3=0.
    END IF
    IF(denom1 /= 0 .AND. denom3 /= 0.) THEN
      projc4=ALOG( projc2/denom1 ) / denom3
!      print *, '  The cone constant is : ',projc4
      IF( projc4 < 0.) THEN
        IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2,/a)')  &
            '  Warning in SETMAPR for Lambert Projection',              &
            '  For the true latitudes provided, ',                      &
            trulat(1),' and ',trulat(2),                                &
            '  projection must be from opposite pole...changing pole.'
        jpole=-jpole
        projc3=TAN(d2rad*(45.-0.5*jpole*trulat(1)) )
        denom2=TAN(d2rad*(45.-0.5*jpole*trulat(2)))
        IF(denom2 /= 0.) THEN
          denom3=ALOG( projc3/denom2 )
        ELSE
          denom3=0.
        END IF
        IF(denom1 /= 0 .AND. denom3 /= 0.) THEN
          projc4=ALOG( projc2/denom1 ) / denom3
!          print *, '  The revised cone constant is : ',projc4
        ELSE
          IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2)')  &
              '  Error (1) in SETMAPR for Lambert Projection',          &
              '  Illegal combination of trulats one: ',                 &
              trulat(1),' and two: ',trulat(2)
          CALL arpsstop("arpsstop called from SETMAPR problems with     &
                       &  trulat",1)
        END IF
      END IF
      projc1=scale*eradius/projc4
    ELSE IF(denom3 == 0. .AND. denom2 /= 0.) THEN   ! tangent
!     IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2)')  &
!         '  Using Tangent Lambert Projection',                         &
!         '  Based on input combination of trulats one: ',              &
!         trulat(1),' and two: ',trulat(2)
      projc4=SIN(d2rad*jpole*trulat(1))
!      print *, '  The cone constant is : ',projc4
      IF( projc4 < 0.) THEN
        IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2,/a)')  &
            '  Warning in SETMAPR for Lambert Projection',              &
            '  For the true latitudes provided, ',                      &
            trulat(1),' and ',trulat(2),                                &
            '  projection must be from opposite pole...changing pole.'
        jpole=-jpole
        projc4=SIN(d2rad*jpole*trulat(1))
      END IF
      projc1=scale*eradius/projc4
    ELSE
      IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2)')  &
          '  Error (1) in SETMAPR for Lambert Projection',              &
          '  Illegal combination of trulats one: ',                     &
          trulat(1),' and two: ',trulat(2)
        CALL arpsstop("arpsstop called from SETMAPR problems with     &
                       &  trulat Lambert Projection",1)
    END IF

    IF(jpole > 0) THEN
!     IF (myproc == 0) WRITE(6,'(a/,a)')  &
!         '  Map projection set to Lambert Conformal',                  &
!         '  X origin, Y origin set to 0.,0. at the North Pole.'
    ELSE
      IF (myproc == 0) WRITE(6,'(a/,a)')  &
          '  Map projection set to Lambert Conformal',                  &
          '  X origin, Y origin set to 0.,0. at the South Pole.'
    END IF
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 3 ) THEN
    trulat(1)=latnot(1)
    rota=orient
    scmap=scale
    projc1=scale*eradius
    projc2=COS(d2rad*trulat(1))
    projc3=projc1*projc2
    IF(projc2 <= 0.) THEN
      IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2)')  &
          '  Error (1) in SETMAPR for Mercator Projection',             &
          '  Illegal true latitude provided: ',trulat(1)
      CALL arpsstop("arpsstop called from SETMAPR problems with         &
                       &  trulat Mercator Projection",1)
    END IF
    IF (myproc == 0) WRITE(6,'(a/,a,f6.1/,a)')  &
        '  Map projection set to Mercator',                             &
        '  X origin, Y origin set to 0.,0. at the equator,',rota,       &
        '  Y positive toward the North Pole.'
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 4 ) THEN
    trulat(1)=latnot(1)
    rota=orient
    scmap=scale
    projc1=scale*eradius
    projc2=COS(d2rad*trulat(1))
    IF(projc2 <= 0.) THEN
      IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2)')  &
          '  Error (1) in SETMAPR for Lat,Lon Projection',              &
          '  Illegal true latitude provided: ',trulat(1)
        CALL arpsstop("arpsstop called from SETMAPR problems with       &
                       &  trulat illegal lat-lon Projection",1)
    END IF
    projc3=projc1*projc2/d2rad
    IF (myproc == 0) WRITE(6,'(a/,a,/a)')  &
        '  Map projection set to Lat, Lon',                             &
        '  X origin, Y origin set to 0.,0. at the equator, 0. long',    &
        '  Y positive toward the North Pole.'
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!                For a 512 km box
!                at a latitude of 75 it is off by 0.8 km
!                                 70              0.4 km
!                                 65              0.3 km
!                                 55              0.1 km
!                                 45              0.05 km
!                                 35              0.02 km
!                                  0              0.01 km
!                at the corners of the box.
!
!  For this projection:
!      projc1 = lat0
!      projc2 = COS(d2rad*lat0)
!      projc3 = lon0
!      projc4 = d2rad*scale*eradius ! deg_to_km
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    trulat(1)=latnot(1)
    rota=orient
    scmap=scale
    projc1 = trulat(1)
    projc2 = COS(d2rad*trulat(1))
    projc3 = orient
    projc4 = d2rad*scale*eradius ! deg_to_km
    IF(trulat(1) <= 0.) THEN
      IF (myproc == 0) WRITE(6,'(a/,a,f9.2,a,f9.2)')  &
          '  Error (1) in SETMAPR for Lat,Lon Projection',              &
          '  Illegal true latitude provided: ',trulat(1)
        CALL arpsstop("arpsstop called from SETMAPR problems with       &
                       &  trulat illegal lat-lon Projection-2",1)
    END IF
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') iproj,' projection is not supported'
      CALL arpsstop("arpsstop called from SETMAPR problems with         &
                     &  iproj",1)
  END IF

  RETURN
END SUBROUTINE setmapr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETMAPR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getmapr(iproj,scale,latnot,orient,x0,y0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Get the constants for the current map projection, which are stored
!  in the common block named /projcst/.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  9/17/94.
!
!  MODIFICATION HISTORY:
!  1/17/96  Corrected retrieval of iproj to assign sign from jpole.
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!    iproj        Map projection number
!                 1=North Polar Stereographic   (-1 South Pole)
!                 2=Northern Lambert Conformal  (-2 Southern)
!                 3=Mercator
!                 4=Lat,Lon
!
!    scale        Map scale factor,  at latitude=latnot
!                 Distance on map = (Distance on earth) * scale
!                 For ARPS model runs, generally this is 1.0
!                 For ARPS plotting this will depend on window
!                 size and the area to be plotted.
!
!    latnot(2)    Real "True" latitude(s) of map projection
!                 (degrees, positive north)
!                 Except for iproj=2, only latnot(1) is used
!
!    orient       Longitude line that runs vertically on the map.
!                 (degrees, negative west, positive east)
!
!    x0           x coordinate of origin
!    y0           y coordinate of origin
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: iproj       ! map projection number
  REAL :: scale          ! map scale factor
  REAL :: latnot(2)      ! true latitude (degrees N)
  REAL :: orient         ! orientation longitude (degrees E)
  REAL :: x0             ! x coordinate of origin
  REAL :: y0             ! y coordinate of origin

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iproj=jproj*jpole
  scale=scmap
  latnot(1)=trulat(1)
  latnot(2)=trulat(2)
  orient=rota
  x0=xorig
  y0=yorig
  RETURN
END SUBROUTINE getmapr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SETORIG                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setorig(iopt,x0,y0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the origin for the map projection.
!  This is call after subroutine mapproj if the origin
!  must be moved from the original position, which is the
!  pole for the polar stereographic projection and the
!  Lambert conformal, and the equator for Mercator.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/20/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    iopt        origin setting option
!                1: origin given in corrdinate x,y
!                2: origin given in lat,lon on earth
!
!    x0          first coordinate of origin
!    y0          second coordinate of origin
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.
  INTEGER :: iopt       ! origin setting option
  REAL :: x0            ! first coordinate of origin
  REAL :: y0            ! second coordinate of origin

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: xnew,ynew,rlat,rlon

!-----------------------------------------------------------------------
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
!  iopt=1 origin is given in x,y in absolute coordinates.
!
!-----------------------------------------------------------------------
!
  IF( iopt == 1 ) THEN
    xorig=x0
    yorig=y0
    CALL xytoll(1,1,0.,0.,rlat,rlon)

!    WRITE(6,'(/a,f18.2,f18.2,/a,f16.2,f16.2/)')                        &
!     '  Coordinate origin set to absolute x,y =',xorig,yorig,          &
!     '    Latitude, longitude= ',rlat,rlon
!
!-----------------------------------------------------------------------
!
!  iopt=2 origin is given in lat,lon on earth
!
!-----------------------------------------------------------------------
!
!
  ELSE IF( iopt == 2 ) THEN
    xorig=0.
    yorig=0.
    CALL lltoxy(1,1,x0,y0,xnew,ynew)
    xorig=xnew
    yorig=ynew
!    write(6,'(/a,f16.2,f16.2,/a,f16.2,f16.2/)')                        &
!     '  Coordinate origin set to absolute x,y =',xorig,yorig,          &
!     '    Latitude, longitude= ',x0,y0

  ELSE
    CALL xytoll(1,1,0.,0.,rlat,rlon)
    IF (myproc == 0) WRITE(6,'(/a,i4,a,/a,f16.2,f16.2,/a,f16.2,f16.2)') &
        ' Setorig option ',iopt,' not supported.',                      &
        '    Coordinate origin unchanged at x,y =',xorig,yorig,         &
        '    Latitude, longitude= ',rlat,rlon
  END IF
  RETURN
END SUBROUTINE setorig
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE XYTOLL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE xytoll(idim,jdim,x,y,rlat,rlon)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine latitude and longitude given X,Y coordinates on
!  map projection.  SETMAPR must be called before this routine
!  to set-up the map projection constants.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/13/93.
!
!  MODIFICATION HISTORY:
!  01/17/96  Bug in southern hemisphere for Polar Stereo and
!            Mercator projections fixed.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    idim     Number of points in x direction.
!    jdim     Number of points in y direction.
!
!    x        Vector of x in map coordinates
!    y        Vector of y in map coordinates
!             Units are meters unless the scale parameter is
!             not equal to 1.0
!
!  OUTPUT:
!
!    rlat     Array of latitude.
!             (degrees, negative south, positive north)
!    rlon     Array of longitude.
!             (degrees, negative west, positive east)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.
  INTEGER :: idim,jdim
  REAL :: x(idim),y(jdim),rlat(idim,jdim),rlon(idim,jdim)

  REAL :: d2rad,r2deg,eradius
  PARAMETER (d2rad=3.141592654/180., r2deg=180./3.141592654,  &
             eradius = 6371000. )  ! mean earth radius in m

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: xabs,yabs,yjp
  REAL :: radius,ratio,dlon

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF ( jproj == 0 ) THEN
    ratio=r2deg/eradius
    DO j = 1, jdim
      DO i = 1, idim
        rlat(i,j) = ratio*(y(j)+yorig)
        rlon(i,j) = ratio*(x(i)+xorig)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim
        yabs=y(j)+yorig
        xabs=x(i)+xorig
        radius=SQRT( xabs*xabs + yabs*yabs )/projc3
        rlat(i,j) = jpole*(90. - 2.*r2deg*ATAN(radius))
        rlat(i,j)=AMIN1(rlat(i,j), 90.)
        rlat(i,j)=AMAX1(rlat(i,j),-90.)

        IF((jpole*yabs) > 0.) THEN
          dlon=180. + r2deg*ATAN(-xabs/yabs)
        ELSE IF((jpole*yabs) < 0.) THEN
          dlon=r2deg*ATAN(-xabs/yabs)
        ELSE IF (xabs > 0.) THEN     ! y=0.
          dlon=90.
        ELSE
          dlon=-90.
        END IF
        rlon(i,j)= rota + jpole*dlon
        IF(rlon(i,j) > 180) rlon(i,j)=rlon(i,j)-360.
        IF(rlon(i,j) < -180) rlon(i,j)=rlon(i,j)+360.
        rlon(i,j)=AMIN1(rlon(i,j), 180.)
        rlon(i,j)=AMAX1(rlon(i,j),-180.)
!
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF ( jproj == 2 ) THEN
    DO j=1,jdim
      DO i=1,idim
        yabs=y(j)+yorig
        xabs=x(i)+xorig
        radius=SQRT( xabs*xabs+ yabs*yabs )
        ratio=projc3*((radius/(projc1*projc2))**(1./projc4))
        rlat(i,j)=jpole*(90. -2.*r2deg*(ATAN(ratio)))
        rlat(i,j)=AMIN1(rlat(i,j), 90.)
        rlat(i,j)=AMAX1(rlat(i,j),-90.)

        yjp=jpole*yabs
        IF(yjp > 0.) THEN
          dlon=180. + r2deg*ATAN(-xabs/yabs)/projc4
        ELSE IF(yjp < 0.) THEN
          dlon=r2deg*ATAN(-xabs/yabs)/projc4
        ELSE IF (xabs > 0.) THEN     ! y=0.
          dlon=90./projc4
        ELSE
          dlon=-90./projc4
        END IF
        rlon(i,j)= rota + jpole*dlon
        IF(rlon(i,j) > 180) rlon(i,j)=rlon(i,j)-360.
        IF(rlon(i,j) < -180) rlon(i,j)=rlon(i,j)+360.
        rlon(i,j)=AMIN1(rlon(i,j), 180.)
        rlon(i,j)=AMAX1(rlon(i,j),-180.)

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 3 ) THEN
    DO j=1,jdim
      DO i=1,idim
        yabs=y(j)+yorig
        xabs=x(i)+xorig
        rlat(i,j)=(90. - 2.*r2deg*ATAN(EXP(-yabs/projc3)))
        rlat(i,j)=AMIN1(rlat(i,j), 90.)
        rlat(i,j)=AMAX1(rlat(i,j),-90.)
        dlon=r2deg*(xabs/projc3)
        rlon(i,j)=rota + dlon
        IF(rlon(i,j) > 180) rlon(i,j)=rlon(i,j)-360.
        IF(rlon(i,j) < -180) rlon(i,j)=rlon(i,j)+360.
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 4 ) THEN
    DO j=1,jdim
      DO i=1,idim
        rlon(i,j)=x(i)+xorig
        rlat(i,j)=y(j)+yorig
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 = lat0
!      projc2 = COS(d2rad*lat0)
!      projc3 = lon0
!      projc4 = d2rad*scale*eradius ! deg_to_km
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    DO j=1,jdim
      DO i=1,idim
        yabs = y(j) + yorig
        xabs = x(i) + xorig
        rlat(i,j) = projc1 + yabs/projc4
        rlat(i,j) = AMIN1(rlat(i,j), 90.)
        rlat(i,j) = AMAX1(rlat(i,j),-90.)
        rlon(i,j) = xabs/projc4/(0.5*(projc2+COS(rlat(i,j)*d2rad))) + projc3
        IF(rlon(i,j) > 180) rlon(i,j)=rlon(i,j)-360.
        IF(rlon(i,j) < -180) rlon(i,j)=rlon(i,j)+360.
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
      CALL arpsstop("arpsstop called from XYTOLL problems with            &
                     &  jproj",1)
  END IF

  RETURN
END SUBROUTINE xytoll
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LLTOXY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE lltoxy(idim,jdim,rlat,rlon,xloc,yloc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine x, y coordinates on map projection from the given latitude
!  and longitude. SETMAPR must be called before this routine to set-up
!  the map projection constants.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/11/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    idim     Array dimension in x direction
!    jdim     Array dimension in y direction
!
!    rlat     Real vector of latitude.
!             (degrees, negative south, positive north)
!
!    rlon     Real vector of longitude.
!             (degrees, negative west, positive east)
!
!  OUTPUT:
!
!    xloc     Real vector of x in map coordinates
!    yloc     Real vector of y in map coordinates
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: idim,jdim
  REAL :: rlat(idim,jdim),rlon(idim,jdim)
  REAL :: xloc(idim,jdim),yloc(idim,jdim)

  REAL :: d2rad,eradius
  PARAMETER (d2rad=3.141592654/180.,                                    &
             eradius = 6371000. )  ! mean earth radius in m

  INTEGER :: jproj,jpole
  REAL :: tem, lat
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: radius,denom,dlon,ratio

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF( jproj == 0 ) THEN
    ratio=d2rad*eradius
    DO j = 1, jdim
      DO i = 1, idim
        xloc(i,j) = ratio*rlon(i,j) - xorig
        yloc(i,j) = ratio*rlat(i,j) - yorig
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim
        denom=(1. + SIN(d2rad*jpole*rlat(i,j)))
        IF(denom == 0.) denom=1.0E-10
        radius=jpole*projc3*COS(d2rad*rlat(i,j))/denom
        dlon=jpole*d2rad*(rlon(i,j)-rota)
        xloc(i,j)= radius*SIN(dlon) - xorig
        yloc(i,j)=-radius*COS(dlon) - yorig
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN

    DO j=1,jdim
      DO i=1,idim

        ! Handle opposite pole
        IF (jpole*rlat(i,j) < -89.9) THEN
          lat = -89.9 * jpole
        ELSE
          lat = rlat(i,j)
        END IF

        radius=projc1*projc2                                            &
              *(TAN(d2rad*(45.-0.5*jpole*lat))/projc3)**projc4
        tem = rlon(i,j)-rota
        IF( tem < -180.0) tem = 360.0+tem
        IF( tem > 180.0) tem = tem-360.0
        dlon=projc4*d2rad*tem
        xloc(i,j)=       radius*SIN(dlon) - xorig
        yloc(i,j)=-jpole*radius*COS(dlon) - yorig
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO j=1,jdim
      DO i=1,idim
        dlon=rlon(i,j)-rota
        IF(dlon < -180.) dlon=dlon+360.
        IF(dlon > 180.) dlon=dlon-360.
        xloc(i,j)=projc3*d2rad*dlon - xorig
        denom=TAN(d2rad*(45. - 0.5*rlat(i,j)))
        IF( denom <= 0. ) denom=1.0E-10
        yloc(i,j)=-projc3*ALOG(denom) - yorig
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO j=1,jdim
      DO i=1,idim
        xloc(i,j)=rlon(i,j)-xorig
        !J.Case, ENSCO Inc. (8/24/2004)
        if (xloc(i,j) < -180.) xloc(i,j)=xloc(i,j)+360.
        if (xloc(i,j) >  180.) xloc(i,j)=xloc(i,j)-360.
        yloc(i,j)=rlat(i,j)-yorig
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 = lat0
!      projc2 = COS(d2rad*lat0)
!      projc3 = lon0
!      projc4 = d2rad*scale*eradius ! deg_to_km
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    DO j=1,jdim
      DO i=1,idim
        xloc(i,j) = projc4*(rlon(i,j) - projc3)  &
                  * 0.5*(projc2 + COS(rlat(i,j)*d2rad)) - xorig
        yloc(i,j) = projc4*(rlat(i,j) - projc1) - yorig
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
      CALL arpsstop("arpsstop called from LLTOXY problems with         &
                     &  jproj not supported",1)
  END IF
  RETURN
END SUBROUTINE lltoxy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LATTOMF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE lattomf(idim,jdim,rlat,emfact)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine the map scale factor, emfact, at a given latitude.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/11/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    idim        Array dimension in x direction
!    jdim        Array dimension in y direction
!
!    rlat        Real vector of latitudes.
!                (degrees, negative south, positive north)
!
!  OUTPUT:
!
!    emfact      Vector of map scale factors corresponding to the
!                input latitudes (map scale includes the projection
!                image scale times the overall scale of the map).
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: idim,jdim         ! dimensions of arrays
  REAL :: rlat(idim,jdim)      ! latitude (degrees)
  REAL :: emfact(idim,jdim)    ! local map scale factor

  REAL :: d2rad
  PARAMETER (d2rad=3.141592654/180.)

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: denom

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF( jproj == 0 ) THEN
    DO j=1,jdim
      DO i=1,idim
        emfact(i,j)=1.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim
        denom=(1. + SIN(d2rad*jpole*rlat(i,j)))
        IF(denom == 0.) denom=1.0E-10
        emfact(i,j)=scmap*projc2/denom
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
    DO j=1,jdim
      DO i=1,idim
        denom=COS( d2rad*rlat(i,j) )
        IF(denom < 1.0E-06) THEN
          emfact(i,j)=1.0E+10
        ELSE
          emfact(i,j)=scmap*(projc2/denom)                              &
                   *(TAN(d2rad*(45.-0.5*jpole*rlat(i,j)))               &
                   /projc3)**projc4
        END IF
        emfact(i,j)=AMAX1(emfact(i,j),1.0E-10)
        emfact(i,j)=AMIN1(emfact(i,j),1.0E+10)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO j=1,jdim
      DO i=1,idim
        denom=COS( d2rad*rlat(i,j) )
        IF(denom == 0.) denom=1.0E-10
        emfact(i,j)=projc2/denom
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO j=1,jdim
      DO i=1,idim
        denom=COS( d2rad*rlat(i,j) )
        IF(denom == 0.) denom=1.0E-10
        emfact(i,j)=projc3/denom
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 is 2 pi times scaled earth's radius divided by 360
!         (i.e. degrees to km conversion factor)
!      projc2 is projc1 * trulat(1)  (i.e. -y0)
!      projc3 is cos of trulat(1)
!      projc4 is trulon
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    DO j=1,jdim
      DO i=1,idim
        emfact(i,j) = 1.0 ! WARNING: the actual map factor for this
                          ! projection has not been determined.
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
    CALL arpsstop('arpsstop called from LATTOMF problems with jproj',1)
  END IF
  RETURN
END SUBROUTINE lattomf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE XYTOMF                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE xytomf(idim,jdim,x,y,emfact)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine the map scale factor, emfact, given x,y in the projected
!  space.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/11/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    idim     Array dimension in x direction.
!    jdim     Array dimension in y direction.
!
!    x        x coordinate values (meters if scmap=1.0)
!    y        y coordinate values (meters if scmap=1.0)
!
!  OUTPUT:
!
!    emfact    Vector of map scale factors corresponding to the
!             input x,y's.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: idim,jdim       ! array dimensions
  REAL :: x(idim)            ! x map coordinate
  REAL :: y(jdim)            ! y map coordinate
  REAL :: emfact(idim,jdim)  ! local map scale factor

  REAL :: d2rad,r2deg
  PARAMETER (d2rad=3.141592654/180.,                                    &
             r2deg=180./3.141592654)

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: xabs,yabs,rlat,ratio,radius,denom
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
  IF( jproj == 0 ) THEN
    DO j=1,jdim
      DO i=1,idim
        emfact(i,j)=1.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim
        xabs=x(i)+xorig
        yabs=y(j)+yorig
        radius=SQRT( xabs*xabs + yabs*yabs )/projc3
        rlat = 90. - 2.*r2deg*ATAN(radius)
        rlat=AMIN1(rlat, 90.)
        rlat=AMAX1(rlat,-90.)
        denom=(1. + SIN(d2rad*rlat))
        IF(denom == 0.) denom=1.0E-10
        emfact(i,j)=scmap*projc2/denom
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
    DO j=1,jdim
      DO i=1,idim
        xabs=x(i)+xorig
        yabs=y(j)+yorig
        radius=SQRT( xabs*xabs+ yabs*yabs )
        ratio=projc3*((radius/(projc1*projc2))**(1./projc4))
        rlat=90. -2.*r2deg*(ATAN(ratio))
        rlat=AMIN1(rlat, 90.)
        rlat=AMAX1(rlat,-90.)
        denom=COS( d2rad*rlat )
        IF(denom == 0.) denom=1.0E-10
        emfact(i,j)=scmap*(projc2/denom)                                &
                   *(TAN(d2rad*(45.-0.5*rlat))/projc3)**projc4
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO j=1,jdim
      yabs=y(j)+yorig
      rlat=90. - 2.*r2deg*ATAN(EXP(-yabs/projc3))
      rlat=AMIN1(rlat, 90.)
      rlat=AMAX1(rlat,-90.)
      denom=COS( d2rad*rlat )
      IF(denom == 0.) denom=1.0E-10
      DO i=1,idim
        emfact(i,j)=projc2/denom
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO j=1,jdim
      yabs=y(j)+yorig
      denom=COS( d2rad*yabs )
      IF(denom == 0.) denom=1.0E-10
      DO i=1,idim
        emfact(i,j)=projc3/denom
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 is 2 pi times scaled earth's radius divided by 360
!         (i.e. degrees to km conversion factor)
!      projc2 is projc1 * trulat(1)  (i.e. -y0)
!      projc3 is cos of trulat(1)
!      projc4 is trulon
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    DO j=1,jdim
      DO i=1,idim
        emfact(i,j) = 1.0 ! WARNING: the actual map factor for this
                          ! projection has not been determined.
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
    CALL arpsstop('arpsstop called from XYTOMF problems with jproj',1)
  END IF
  RETURN
END SUBROUTINE xytomf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DDROTUV                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ddrotuv(nsta,stalon,dd,ff,ddrot,umap,vmap)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Rotate wind from earth direction to map orientation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/20/93.
!
!  MODIFICATION HISTORY:
!  03/30/95  (K. Brewster)
!  Removed the map scale factor from the conversion of winds
!  from u,v on the earth to projection u,v.  Affected argument
!  list of ddrotuv.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nsta      array dimension
!
!    stalon    longitude (degrees E)
!
!    dd        wind direction (degrees from north)
!    ff        wind speed
!
!  OUTPUT:
!
!    ddrot     wind direction rotated to map orientation
!
!    umap      u wind component on map (same units as ff)
!    vmap      v wind component on map (same units as ff)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: nsta               ! array dimension
  REAL :: stalon(nsta)          ! longitude (degrees E)
  REAL :: dd(nsta)              ! wind direction
  REAL :: ff(nsta)              ! speed
  REAL :: ddrot(nsta)           ! wind direction rotated to map orientation
  REAL :: umap(nsta)            ! u wind component on map
  REAL :: vmap(nsta)            ! v wind component on map

  REAL :: d2rad,r2deg
  PARAMETER (d2rad=3.141592654/180.,                                    &
             r2deg=180./3.141592654)

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
  REAL :: arg

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection.
!  Just do conversion from ddff to u,v.
!
!-----------------------------------------------------------------------
!
  IF( jproj == 0 ) THEN
    DO i=1,nsta
      ddrot(i)=dd(i)
      arg = (ddrot(i) * d2rad)
      umap(i) = -ff(i) * SIN(arg)
      vmap(i) = -ff(i) * COS(arg)
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO i=1,nsta
      ddrot(i)=dd(i) + rota - stalon(i)
      arg = (ddrot(i) * d2rad)
      umap(i) = -ff(i) * SIN(arg)
      vmap(i) = -ff(i) * COS(arg)
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
    DO i=1,nsta
      ddrot(i)=dd(i) + projc4*(rota - stalon(i))
      arg = (ddrot(i) * d2rad)
      umap(i) = -ff(i) * SIN(arg)
      vmap(i) = -ff(i) * COS(arg)
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO i=1,nsta
      ddrot(i)=dd(i)
      arg = (ddrot(i) * d2rad)
      umap(i) = -ff(i) * SIN(arg)
      vmap(i) = -ff(i) * COS(arg)
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO i=1,nsta
      ddrot(i)=dd(i)
      arg = (ddrot(i) * d2rad)
      umap(i) = -ff(i) * SIN(arg)
      vmap(i) = -ff(i) * COS(arg)
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 is 2 pi times scaled earth's radius divided by 360
!         (i.e. degrees to km conversion factor)
!      projc2 is projc1 * trulat(1)  (i.e. -y0)
!      projc3 is cos of trulat(1)
!      projc4 is trulon
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    DO i=1,nsta
      ! WARNING: this projection is not conformal.  The following is only
      ! approximately correct.
      ddrot(i)=dd(i)
      arg = (ddrot(i) * d2rad)
      umap(i) = -ff(i) * SIN(arg)
      vmap(i) = -ff(i) * COS(arg)
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
    CALL arpsstop('arpsstop called from DDROTUV problems with jproj',1)
  END IF
  RETURN
END SUBROUTINE ddrotuv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UVROTDD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uvrotdd(idim,jdim,elon,umap,vmap,dd,ff)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Convert winds u, v in map coordinates to wind direction and speed
!  in earth coordinates.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  11/20/93.
!
!  MODIFICATION HISTORY:
!  03/30/95  (K. Brewster)
!  Removed the map scale factor from the conversion of winds
!  from u,v on the earth to projection u,v.  Affected argument
!  list of uvrotdd.
!
!  01/22/03  (K. W. Thomas)
!  Negative wind directions could be returned.  Prevent it from happening.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    idim       Array dimension in the x direction
!    jdim       Array dimension in the y direction
!
!    elon       Earth longitude (degrees E)
!
!    umap       u wind component on map
!    vmap       v wind component on map
!
!  OUTPUT:
!    dd         wind direction on earth
!    ff         wind speed on earth
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: idim,jdim       ! array dimensions
  REAL :: elon(idim,jdim)    ! longitude (degrees E)
  REAL :: umap(idim,jdim)    ! u wind component on map
  REAL :: vmap(idim,jdim)    ! v wind component on map

  REAL :: dd(idim,jdim)      ! direction
  REAL :: ff(idim,jdim)      ! wind speed

  REAL :: d2rad,r2deg
  PARAMETER (d2rad=3.141592654/180.,                                    &
             r2deg=180./3.141592654)

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: dlon

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF( jproj == 0 ) THEN
    DO j=1,jdim
      DO i=1,idim
        ff(i,j) = SQRT(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

        IF(vmap(i,j) > 0.) THEN
          dlon=r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(vmap(i,j) < 0.) THEN
          dlon=180. + r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(umap(i,j) >= 0.) THEN
          dlon=90.
        ELSE
          dlon=-90.
        END IF

        dd(i,j)= dlon + 180.
        dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
        if ( dd(i,j) < 0.0 ) dd(i,j) = dd(i,j) + 360.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim

        ff(i,j) = SQRT(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

        IF(vmap(i,j) > 0.) THEN
          dlon=r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(vmap(i,j) < 0.) THEN
          dlon=180. + r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(umap(i,j) >= 0.) THEN
          dlon=90.
        ELSE
          dlon=-90.
        END IF

        dd(i,j)= dlon + 180. + elon(i,j) - rota
        dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
        if ( dd(i,j) < 0.0 ) dd(i,j) = dd(i,j) + 360.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
    DO j=1,jdim
      DO i=1,idim
        ff(i,j) = SQRT(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

        IF(vmap(i,j) > 0.) THEN
          dlon=r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(vmap(i,j) < 0.) THEN
          dlon=180. + r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(umap(i,j) >= 0.) THEN
          dlon=90.
        ELSE
          dlon=-90.
        END IF

        dd(i,j)= dlon + 180. + projc4*(elon(i,j) - rota)
        dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
	if ( dd(i,j) < 0.0 ) dd(i,j) = dd(i,j) + 360.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO j=1,jdim
      DO i=1,idim
        ff(i,j) = SQRT(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

        IF(vmap(i,j) > 0.) THEN
          dlon=r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(vmap(i,j) < 0.) THEN
          dlon=180. + r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(umap(i,j) >= 0.) THEN
          dlon=90.
        ELSE
          dlon=-90.
        END IF

        dd(i,j)= dlon + 180.
        dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
	if ( dd(i,j) < 0.0 ) dd(i,j) = dd(i,j) + 360.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO j=1,jdim
      DO i=1,idim
        ff(i,j) = SQRT(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

        IF(vmap(i,j) > 0.) THEN
          dlon=r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(vmap(i,j) < 0.) THEN
          dlon=180. + r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(umap(i,j) >= 0.) THEN
          dlon=90.
        ELSE
          dlon=-90.
        END IF

        dd(i,j)= dlon + 180.
        dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
        if ( dd(i,j) < 0.0 ) dd(i,j) = dd(i,j) + 360.0
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 is 2 pi times scaled earth's radius divided by 360
!         (i.e. degrees to km conversion factor)
!      projc2 is projc1 * trulat(1)  (i.e. -y0)
!      projc3 is cos of trulat(1)
!      projc4 is trulon
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    ! WARNING: this projection is not conformal.  The following is only
    ! approximately correct.
    DO j=1,jdim
      DO i=1,idim
        ff(i,j) = SQRT(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

        IF(vmap(i,j) > 0.) THEN
          dlon=r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(vmap(i,j) < 0.) THEN
          dlon=180. + r2deg*ATAN(umap(i,j)/vmap(i,j))
        ELSE IF(umap(i,j) >= 0.) THEN
          dlon=90.
        ELSE
          dlon=-90.
        END IF

        dd(i,j)= dlon + 180.
        dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
	if ( dd(i,j) < 0.0 ) dd(i,j) = dd(i,j) + 360.0
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
    CALL arpsstop('arpsstop called from UVROTDD problems with jproj',1)
  END IF
  RETURN
END SUBROUTINE uvrotdd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UVETOMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uvetomp(idim,jdim,uear,vear,lon,umap,vmap)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Transform u, v wind from earth coordinates to map coordinates.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  04/30/94.
!
!  MODIFICATION HISTORY:
!  03/30/95  (K. Brewster)
!  Removed the map scale factor from the conversion of winds
!  from u,v on the earth to projection u,v.  Affected argument
!  list of uvetomp.
!  04/30/96  (KB)
!  Streamlined the computation for iproj=1 and iproj=2.
!  12/11/96  (KB)
!  Corrected a bug in the computation for iproj=1 and iproj=2.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    idim       Array dimension in the x direction
!    jdim       Array dimension in the y direction
!
!    uear       u (eastward) wind component on earth
!    vear       v (northwrd) wind component on earth
!
!    lon        earth longitude
!
!  OUTPUT:
!
!    umap       u wind component on map
!    vmap       v wind component on map
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: idim,jdim        ! array dimensions
  REAL :: uear(idim,jdim)     ! u (eastward) wind component on earth
  REAL :: vear(idim,jdim)     ! v (northward) wind component on earth
  REAL :: lon(idim,jdim)      ! longitude (degrees east)

  REAL :: umap(idim,jdim)     ! u wind component on map
  REAL :: vmap(idim,jdim)     ! v wind component on map

  REAL :: d2rad,r2deg
  PARAMETER (d2rad=3.141592654/180.,                                    &
             r2deg=180./3.141592654)

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: dlon,arg,dxdlon,dydlon,utmp,vtmp
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF( jproj == 0 ) THEN
    DO j=1,jdim
      DO i=1,idim
        umap(i,j) = uear(i,j)
        vmap(i,j) = vear(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim
        dlon=(lon(i,j)-rota)
        arg=d2rad*dlon
        dxdlon=COS(arg)
        dydlon=SIN(arg)*jpole  ! by mxue on 8/28
        utmp=uear(i,j)
        vtmp=vear(i,j)
        umap(i,j)=utmp*dxdlon - vtmp*dydlon
        vmap(i,j)=vtmp*dxdlon + utmp*dydlon
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
    DO j=1,jdim
      DO i=1,idim
        dlon=(lon(i,j)-rota)
        arg=d2rad*projc4*(dlon - 360.*nint(dlon/360.))
        dxdlon=COS(arg)
        dydlon=SIN(arg)*jpole   ! by mxue on 8/28
        utmp=uear(i,j)
        vtmp=vear(i,j)
        umap(i,j)=utmp*dxdlon - vtmp*dydlon
        vmap(i,j)=vtmp*dxdlon + utmp*dydlon
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO j=1,jdim
      DO i=1,idim
        umap(i,j) = uear(i,j)
        vmap(i,j) = vear(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO j=1,jdim
      DO i=1,idim
        umap(i,j) = uear(i,j)
        vmap(i,j) = vear(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 is 2 pi times scaled earth's radius divided by 360
!         (i.e. degrees to km conversion factor)
!      projc2 is projc1 * trulat(1)  (i.e. -y0)
!      projc3 is cos of trulat(1)
!      projc4 is trulon
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    ! WARNING: this projection is not conformal.  The following is only
    ! approximately correct.
    DO j=1,jdim
      DO i=1,idim
        umap(i,j) = uear(i,j)
        vmap(i,j) = vear(i,j)
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
    CALL arpsstop('arpsstop called from UVEROTDD problems with jproj',1)
  END IF
  RETURN
END SUBROUTINE uvetomp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UVMPTOE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uvmptoe(idim,jdim,umap,vmap,lon,uear,vear)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Transform u, v wind from map coordinates to earth coordinates.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  04/30/94.
!
!  MODIFICATION HISTORY:
!  03/30/95  (K. Brewster)
!  Removed the map scale factor from the conversion of winds
!  from u,v on the map to earth u,v.  Affected argument
!  list of uvmptoe.
!  04/30/96  (KB)
!  Streamlined the computation for iproj=1 and iproj=2.
!  12/11/96  (KB)
!  Corrected a bug in the computation for iproj=1 and iproj=2.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    idim           Array dimension in x direction
!    jdim           Array dimension in y direction
!
!    umap           u wind component on map
!    vmap           v wind component on map
!
!    lon            Longitude (degrees E)
!
!  OUTPUT:
!
!    uear           u (eastward) wind component on earth
!    vear           v (northward) wind component on earth
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'      ! Message passing parameters.

  INTEGER :: idim,jdim       ! array dimensions
  REAL :: lon(idim,jdim)     ! longitude (degrees E)
  REAL :: umap(idim,jdim)    ! u wind component on map
  REAL :: vmap(idim,jdim)    ! v wind component on map

  REAL :: uear(idim,jdim)    ! u (eastward) wind component on earth
  REAL :: vear(idim,jdim)    ! v (northward) wind component on earth

  REAL :: d2rad,r2deg
  PARAMETER (d2rad=3.141592654/180.,                                    &
             r2deg=180./3.141592654)

  INTEGER :: jproj,jpole
  REAL :: trulat(2),rota,scmap,xorig,yorig,                             &
       projc1,projc2,projc3,projc4,projc5
  COMMON /projcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,           &
                   projc1,projc2,projc3,projc4,projc5
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: dlon,arg,utmp,vtmp,dxdlon,dydlon
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  No map projection
!
!-----------------------------------------------------------------------
!
  IF( jproj == 0 ) THEN
    DO j=1,jdim
      DO i=1,idim
        uear(i,j) = umap(i,j)
        vear(i,j) = vmap(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Polar Stereographic projection
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is the numerator of emfact, the map image scale factor.
!      projc3 is projc2 times the scaled earth's radius.
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 1 ) THEN
    DO j=1,jdim
      DO i=1,idim
        dlon=(lon(i,j)-rota)
        arg=d2rad*dlon
        dxdlon=COS(arg)
        dydlon=SIN(arg)*jpole   ! by mxue on 8/28/2008
        utmp=umap(i,j)
        vtmp=vmap(i,j)
        uear(i,j)=utmp*dxdlon + vtmp*dydlon
        vear(i,j)=vtmp*dxdlon - utmp*dydlon
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal Conic Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius/n
!      projc2 is cos of trulat(1)
!      projc3 is tan (45. - trulat/2) a const for local map scale
!      projc4 is the cone constant, n
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 2 ) THEN
!        uear(1,1)=-9.20149E+00
!        vear(1,1)=-4.53746E+00
    DO j=1,jdim
      DO i=1,idim
        dlon=(lon(i,j)-rota)
        arg=d2rad*projc4*(dlon - 360.*nint(dlon/360.))
        dxdlon=COS(arg)
        dydlon=SIN(arg)*jpole  ! by mxue on 8/28/2008
!2001-05-16 GMB: Having umap & uear (or vmap & vear) point to
!the same array causes numerical errors when optimizing.
        uear(i,j) = umap(i,j)*dxdlon + vmap(i,j)*dydlon
        vear(i,j) = vmap(i,j)*dxdlon - umap(i,j)*dydlon
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Mercator Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 3) THEN
    DO j=1,jdim
      DO i=1,idim
        uear(i,j) = umap(i,j)
        vear(i,j) = vmap(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Lat, Lon Projection.
!  For this projection:
!      projc1 is the scaled earth's radius, scale times eradius
!      projc2 is cos of trulat(1)
!      projc3 is projc1 times projc2 times 180/pi
!
!-----------------------------------------------------------------------
!
  ELSE IF(jproj == 4) THEN
    DO j=1,jdim
      DO i=1,idim
        uear(i,j) = umap(i,j)
        vear(i,j) = vmap(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  WDT mapproj
!
!  Approximate flat earth projection (using approximate great circle to
!  compute distances).
!
!  For this projection:
!      projc1 is 2 pi times scaled earth's radius divided by 360
!         (i.e. degrees to km conversion factor)
!      projc2 is projc1 * trulat(1)  (i.e. -y0)
!      projc3 is cos of trulat(1)
!      projc4 is trulon
!
!-----------------------------------------------------------------------
!
  ELSE IF( jproj == 5 ) THEN
    ! WARNING: this projection is not conformal.  The following is only
    ! approximately correct.
    DO j=1,jdim
      DO i=1,idim
        uear(i,j) = umap(i,j)
        vear(i,j) = vmap(i,j)
      END DO
    END DO
  ELSE
    IF (myproc == 0) WRITE(6,'(i4,a)') jproj,' projection is not supported'
    CALL arpsstop('arpsstop called from UVMPTOE problems with jproj',1)
  END IF

  RETURN
END SUBROUTINE uvmptoe
