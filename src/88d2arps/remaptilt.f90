SUBROUTINE radcoord(nx,ny,nz,                                          &
                    x,y,z,zp,xs,ys,zps,                                &
                    radar_lat,radar_lon,radarx,radary)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Intialize coordinate fields for the radar remapping routines.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    radar_lat Latitude (degrees) of radar
!    radar_lon Longitude (positive East) of radar
!
!  OUTPUT:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    xs       x coordinate of scalar grid points in physical/comp. space (m)
!    ys       y coordinate of scalar grid points in physical/comp. space (m)
!    zps      Vertical coordinate of scalar grid points in physical space(m)
!
!    radarx   x location of radar in grid
!    radary   y location of radar in grid
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx
  INTEGER, INTENT(IN)  :: ny
  INTEGER, INTENT(IN)  :: nz

  REAL, INTENT(IN)     :: x(nx)
  REAL, INTENT(IN)     :: y(ny)
  REAL, INTENT(IN)     :: z(nz)
  REAL, INTENT(IN)     :: zp(nx,ny,nz)

  REAL, INTENT(OUT)    :: xs(nx)
  REAL, INTENT(OUT)    :: ys(ny)
  REAL, INTENT(OUT)    :: zps(nx,ny,nz)

  REAL, INTENT(IN)     :: radar_lat
  REAL, INTENT(IN)     :: radar_lon

  REAL, INTENT(OUT)    :: radarx
  REAL, INTENT(OUT)    :: radary

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
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
! Set the coordinates of the scalar grid points
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=2.0*xs(nx-1)-xs(nx-2)

  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=2.0*ys(ny-1)-ys(ny-2)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
  DO j=1,ny-1
    DO i=1,nx-1
      zps(i,j,nz)=2.0*zps(i,j,nz-1)-zps(i,j,nz-2)
    END DO
  END DO
!
  CALL lltoxy(1,1,radar_lat,radar_lon,radarx,radary)

!  print *, ' radar_lat,radar_lon: ',radar_lat,radar_lon
!
!  print *, ' grid radarx = ',radarx,' radary = ',radary

  RETURN
END SUBROUTINE radcoord
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE RMPINIT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rmpinit(nx,ny,nz,                                          &
              mxrefgat,mxvelgat,maxazim,maxelev,dualpproc,            &
              kntrgat,kntrazm,kntrelv,kntvgat,kntvazm,kntvelv,        &
              nyqvvol,timevolr,timevolv,                              &
              rngrvol,azmrvol,elvrvol,elvmnrvol,                      &
              refvol,rhvvol,zdrvol,kdpvol,                            &
              rngvvol,azmvvol,elvvvol,elvmnvvol,velvol,               &
              gridvel,gridref,gridnyq,gridtim)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Intialize fields for the radar remapping routines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!
!  MODIFICATION HISTORY:
!
!  05/17/2002 (Keith Brewster)
!  Added initialization of gridvel,gridref,gridnyq,gridtim
!  to missing value.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    mxrefgat Maximum number of reflectivity gates
!    mxvelgat Maximum number of velocity gates
!    maxazim  Maximum number of azimuth angles per tilt
!    maxelev  Maximum number of elevation angles (tilts)
!    dualpproc  Process dual-pol data?
!
!  OUTPUT
!
!    kntrgat  Number of reflectivity gates
!    kntrazm  Number of reflectivity radials
!    kntvelv  Number of reflectivity tilts
!
!    kntvgat  Number of velocity gates
!    kntvazm  Number of velocity radials
!    kntvelv  Number of velocity tilts
!
!    rngrvol  Range to gate in reflectivity 3-D volume
!    azmrvol  Azimuth angle in reflectivity 3-D volume
!    elvrvol  Elevation angle in reflectivity 3-D volume
!    refvol   Reflectivity 3-D volume
!    rhvvol   Rho-HV 3-D volume
!    zdrvol   Zdr 3-D volume
!    kdpvol   Kdp 3-D volume
!
!    rngvvol  Range to gate in velocity 3-D volume
!    azmvvol  Azimuth angle in velocity 3-D volume
!    elvvvol  Elevation angle in velocity 3-D volume
!    velvol   Velocity 3-D volume
!
!    gridvel  radial velocity at grid points
!    gridref  reflectivity at grid points
!    gridnyq  Nyquist velocity at grid points
!    gridtim  time (offset) at grid points
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

  INTEGER, INTENT(IN)  :: nx
  INTEGER, INTENT(IN)  :: ny
  INTEGER, INTENT(IN)  :: nz

  INTEGER, INTENT(IN)  :: mxrefgat
  INTEGER, INTENT(IN)  :: mxvelgat
  INTEGER, INTENT(IN)  :: maxazim
  INTEGER, INTENT(IN)  :: maxelev

  LOGICAL, INTENT(IN)  :: dualpproc

  INTEGER, INTENT(OUT) :: kntrgat(maxazim,maxelev)
  INTEGER, INTENT(OUT) :: kntrazm(maxelev)
  INTEGER, INTENT(OUT) :: kntrelv

  INTEGER, INTENT(OUT) :: kntvgat(maxazim,maxelev)
  INTEGER, INTENT(OUT) :: kntvazm(maxelev)
  INTEGER, INTENT(OUT) :: kntvelv

  REAL,    INTENT(OUT) :: nyqvvol(maxazim,maxelev)
  INTEGER, INTENT(OUT) :: timevolr(maxazim,maxelev)
  INTEGER, INTENT(OUT) :: timevolv(maxazim,maxelev)
  REAL, INTENT(OUT)    :: rngrvol(mxrefgat,maxelev)
  REAL, INTENT(OUT)    :: azmrvol(maxazim,maxelev)
  REAL, INTENT(OUT)    :: elvrvol(maxazim,maxelev)
  REAL, INTENT(OUT)    :: elvmnrvol(maxelev)
  REAL, INTENT(OUT)    :: refvol(mxrefgat,maxazim,maxelev)
  REAL, INTENT(OUT)    :: rhvvol(mxrefgat,maxazim,maxelev)
  REAL, INTENT(OUT)    :: zdrvol(mxrefgat,maxazim,maxelev)
  REAL, INTENT(OUT)    :: kdpvol(mxrefgat,maxazim,maxelev)
  REAL, INTENT(OUT)    :: rngvvol(mxvelgat,maxelev)
  REAL, INTENT(OUT)    :: azmvvol(maxazim,maxelev)
  REAL, INTENT(OUT)    :: elvvvol(maxazim,maxelev)
  REAL, INTENT(OUT)    :: elvmnvvol(maxelev)
  REAL, INTENT(OUT)    :: velvol(mxvelgat,maxazim,maxelev)

  REAL, INTENT(OUT)    :: gridvel(nz,nx,ny)
  REAL, INTENT(OUT)    :: gridref(nz,nx,ny)
  REAL, INTENT(OUT)    :: gridnyq(nz,nx,ny)
  REAL, INTENT(OUT)    :: gridtim(nz,nx,ny)
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
! Initialize counters to zero and fill data arrays with missing value
!
!-----------------------------------------------------------------------
!
  kntrelv = 0
  kntrazm = 0
  kntrgat = 0
  kntvelv = 0
  kntvazm = 0
  kntvgat = 0

  nyqvvol = 0.
  timevolr = 0
  timevolv = 0
  elvrvol = -999.
  elvmnrvol = -999.
  azmrvol = -999.
  rngrvol = -999.
  refvol  = -999.
  IF( dualpproc ) THEN
    rhvvol = -999.
    zdrvol = -999.
    kdpvol = -999.
  END IF
  elvvvol = -999.
  elvmnvvol = -999.
  azmvvol = -999.
  rngvvol = -999.
  velvol  = -999.

  gridvel = -999.
  gridref = -999.
  gridnyq = -999.
  gridtim = -999.

  RETURN
END SUBROUTINE rmpinit
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE VOLBUILD                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE volbuild(maxgate,maxazim,maxelev,ngate,nazim,              &
                    nyqset,timeset,                                   &
                    kntgate,kntazim,kntelev,                          &
                    gatesp,rfirstg,varchek,                           &
                    vnyquist,time,                                    &
                    azim,elev,vartilt,                                &
                    nyqvvol,timevol,                                  &
                    rngvol,azmvol,elvvol,elvmnvol,varvol)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Build a 3-D volume of radar data from 2-D tilts.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!
!  MODIFICATION HISTORY:
!
!  2002/07/25 Yunheng Wang
!  Changed timevol from REAL to INTEGER.
!
!  2003/07/01 Keith Brewster
!  Changed input time variable to INTEGER.
!  Added a few enhancements to documentation.
!
!  2008/03/12 Keith Brewster
!  Made input vnyquist a function of azimuth.
!  Made output nyqvvol a function of azimuth and elev.
!
!  2008/03/12 Keith Brewster
!  Made input vnyquist a function of azimuth.
!  Made output nyqvvol a function of azimuth and elev.
!
!  2010/08/24 Keith Brewster
!  Added elvmnvol variable and sorting of tilts based on mean elevation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    ngate     Number of gates in each radial
!    nazim     Number of radials in each tilt
!
!    nyqset    Flag whether to save and use nyquist data for this tilt
!              0: don't use   1: save and use
!    timeset   Flag whether to save and use time data for this tilt
!              0: don't use   1: save and use
!
!    gatesp    Gate spacing
!    rfirstg   Range (m) to first gate
!    varchek   Threshold value to determine good vs. flagged data
!
!    vnyquist  Nyquist velocity
!    time      Time (offset)
!
!    azim      Azimuth angles
!    elev      Elevation angles
!    vartilt   Radar variable (reflectivity, velocity or spectrum width)
!
!  OUTPUT:
!
!    kntgate   Number of gates in each radial (3-D)
!    kntazim   Number of radials in each tilt (3-D)
!    kntelev   Number of elevation angles (tilts)
!
!    nyqvvol   Nyquist velocity
!    timevol   Time (offset)
!    rngvol    Range to gate
!    azmvol    Azimuth angle
!    elvvol    Elevation angle
!    varvol    Radar variables accumulated in 3-D array
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

  INTEGER, INTENT(IN)  :: maxgate
  INTEGER, INTENT(IN)  :: maxazim
  INTEGER, INTENT(IN)  :: maxelev
  INTEGER, INTENT(IN)  :: ngate
  INTEGER, INTENT(IN)  :: nazim

  INTEGER, INTENT(IN)  :: nyqset
  INTEGER, INTENT(IN)  :: timeset

  INTEGER, INTENT(INOUT) :: kntgate(maxazim,maxelev)
  INTEGER, INTENT(INOUT) :: kntazim(maxelev)
  INTEGER, INTENT(INOUT) :: kntelev

  INTEGER, INTENT(IN)  :: gatesp
  INTEGER, INTENT(IN)  :: rfirstg
  REAL,    INTENT(IN)  :: varchek

  REAL,    INTENT(IN)  :: vnyquist(maxazim)
  INTEGER, INTENT(IN)  :: time(maxazim)
  REAL,    INTENT(IN)  :: azim(maxazim)
  REAL,    INTENT(IN)  :: elev(maxazim)
  REAL,    INTENT(IN)  :: vartilt(maxgate,maxazim)

  REAL,    INTENT(INOUT) :: nyqvvol(maxazim,maxelev)
  INTEGER, INTENT(INOUT) :: timevol(maxazim,maxelev)
  REAL,    INTENT(INOUT) :: rngvol(maxgate,maxelev)
  REAL,    INTENT(INOUT) :: azmvol(maxazim,maxelev)
  REAL,    INTENT(INOUT) :: elvvol(maxazim,maxelev)
  REAL,    INTENT(INOUT) :: elvmnvol(maxelev)
  REAL,    INTENT(INOUT) :: varvol(maxgate,maxazim,maxelev)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: igate,jazim,kelev,k
  REAL :: elvsum,elvmean
  LOGICAL :: duplevel
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( nazim > 0 ) THEN
!
! First determine mean elevation angle of this tilt and proper
! vertical array index, kelev to store this data in the volume array
!
    elvsum=0.
    DO jazim=1,nazim
      elvsum=elvsum+elev(jazim)
    END DO
    elvmean=elvsum/float(nazim)

    IF(kntelev > 0) THEN

      duplevel=.false.
      DO kelev=1,kntelev
        IF(abs(elvmnvol(kelev) - elvmean) < 0.2 ) THEN
          duplevel=.true.
          EXIT
        END IF
      END DO

      IF(duplevel) THEN

        WRITE(6,'(a,i4,2f8.2)')  &
        ' volbuild duplicate elevation angle: ',kelev,elvmean,elvmnvol(kelev)

!
! Reset variables on kelev to missing
!
        kntazim(kelev)=0
        DO jazim = 1,maxazim
          kntgate(jazim,kelev)=0
          nyqvvol(jazim,kelev)=0.
          timevol(jazim,kelev)=0.
          azmvol(jazim,kelev)=-999.
          elvvol(jazim,kelev)=-999.
        END DO
        DO igate=1,maxgate
          rngvol(igate,kelev)=-999.
        END DO
        DO jazim=1,maxazim
          DO igate=1,maxgate
            varvol(igate,jazim,kelev)=-999.
          END DO
        END DO
      ELSE

        DO kelev=1,kntelev
          IF(elvmnvol(kelev) > elvmean) EXIT
        END DO
!
! When the vertical index is not above all current data in volume,
! shuffle all the current data up one index.
!
        IF( kelev <= kntelev ) THEN
          WRITE(6,'(a,f8.2,a,i4,a,i4)')                                &
            ' volbuild new elev: ',elvmean,' at k index ',kelev,' of ', &
              kntelev,' existing levels'
          DO k=kntelev,kelev,-1
            elvmnvol(k+1)=elvmnvol(k)
            kntazim(k+1)=kntazim(k)
            DO jazim=1,maxazim
              kntgate(jazim,k+1)=kntgate(jazim,k)
              nyqvvol(jazim,k+1)=nyqvvol(jazim,k)
              timevol(jazim,k+1)=timevol(jazim,k)
              azmvol(jazim,k+1)=azmvol(jazim,k)
              elvvol(jazim,k+1)=elvvol(jazim,k)
            END DO
            DO igate=1,maxgate
              rngvol(igate,k+1)=rngvol(igate,k)
            END DO
            DO jazim=1,maxazim
              DO igate=1,maxgate
                varvol(igate,jazim,k+1)=varvol(igate,jazim,k)
              END DO
            END DO
          END DO
!
! Reset variables on kelev to missing
!
          kntazim(kelev)=0
          DO jazim = 1,maxazim
            kntgate(jazim,kelev)=0
            nyqvvol(jazim,kelev)=0.
            timevol(jazim,kelev)=0.
            azmvol(jazim,kelev)=-999.
            elvvol(jazim,kelev)=-999.
          END DO
          DO igate=1,maxgate
            rngvol(igate,kelev)=-999.
          END DO
          DO jazim=1,maxazim
            DO igate=1,maxgate
              varvol(igate,jazim,kelev)=-999.
            END DO
          END DO
        END IF  ! need to insert data level in middle

        WRITE(6,'(a,f8.2)') ' volbuild adding elevation angle: ',elvmean
        kntelev = kntelev + 1

      END IF ! not dupe

    ELSE

      WRITE(6,'(a,f8.2)') ' volbuild first elevation angle to store: ',elvmean
      kntelev = 1
      kelev = 1

    END IF
!
! Store incoming data in vertical level kelev
!
    kntazim(kelev) = nazim
    elvmnvol(kelev) = elvmean
    IF( nyqset > 0 ) THEN
      DO jazim = 1, nazim
        nyqvvol(jazim,kelev)=vnyquist(jazim)
      END DO
    END IF

    DO igate = 1, maxgate
      rngvol(igate,kelev)=rfirstg+(igate-1)*gatesp
    END DO
    WRITE(6,'(a,i9)') ' volbuild ngate:',ngate
    WRITE(6,'(a,f12.1,a,f12.1)')' volbuild min range:',rngvol(1,kelev), &
                                '  max range:',rngvol(ngate,kelev)

    IF( timeset > 0 ) THEN
      DO jazim = 1, nazim
        timevol(jazim,kelev)=time(jazim)
        azmvol(jazim,kelev)=azim(jazim)
        elvvol(jazim,kelev)=elev(jazim)
      END DO
    ELSE
      DO jazim = 1, nazim
        azmvol(jazim,kelev)=azim(jazim)
        elvvol(jazim,kelev)=elev(jazim)
      END DO
    END IF

    DO jazim = 1, nazim
      DO igate=1, ngate
        varvol(igate,jazim,kelev)=vartilt(igate,jazim)
        IF (vartilt(igate,jazim) > varchek)                               &
          kntgate(jazim,kelev)=max(kntgate(jazim,kelev),igate)
      END DO
    END DO

    WRITE(6,'(a)')  ' Current Mean Elevation Angles in Volume '
    WRITE(6,'(a)')  '  Index  Elevation   N azimuths'
    DO k = 1, kntelev
      WRITE(6,'(1x,i6,f12.2,i10)') k,elvmnvol(k),kntazim(k)
    END DO

  END IF

  RETURN
END SUBROUTINE volbuild
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE VOLBUILDDLP                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE volbuilddlp(maxgate,maxazim,maxelev,ngate,nazim,           &
                    nyqset,timeset,                                   &
                    kntgate,kntazim,kntelev,                          &
                    gatesp,rfirstg,refchek,                           &
                    vnyquist,time,                                    &
                    azim,elev,reftilt,rhvtilt,zdrtilt,kdptilt,        &
                    nyqvvol,timevol,                                  &
                    rngvol,azmvol,elvvol,elvmnvol,                    &
                    refvol,rhvvol,zdrvol,kdpvol)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Build a 3-D volume of radar data from 2-D tilts, including dual-pol
!  variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  2012/09/11  After volbuild
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    ngate     Number of gates in each radial
!    nazim     Number of radials in each tilt
!
!    nyqset    Flag whether to save and use nyquist data for this tilt
!              0: don't use   1: save and use
!    timeset   Flag whether to save and use time data for this tilt
!              0: don't use   1: save and use
!
!    gatesp    Gate spacing
!    rfirstg   Range (m) to first gate
!    refchek   Threshold value to determine good vs. flagged reflectivity data
!
!    vnyquist  Nyquist velocity
!    time      Time (offset)
!
!    azim      Azimuth angles
!    elev      Elevation angles
!    reftilt   Reflectivity data on a tilt
!    rhvtilt   Rho-HV data on a tilt
!    zdrtilt   Zdr data on a tilt
!    kdptilt   Kdp data on a tilt
!
!  OUTPUT:
!
!    kntgate   Number of gates in each radial (3-D)
!    kntazim   Number of radials in each tilt (3-D)
!    kntelev   Number of elevation angles (tilts)
!
!    nyqvvol   Nyquist velocity
!    timevol   Time (offset)
!    rngvol    Range to gate
!    azmvol    Azimuth angle
!    elvvol    Elevation angle
!    refvol    Reflectivity data accumulated in 3-D array
!    rhvvol    Rho-HV data accumulated in 3-D array
!    zdrvol    Zdr data accumulated in 3-D array
!    kdpvol    Kdp data accumulated in 3-D array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: maxgate
  INTEGER, INTENT(IN)  :: maxazim
  INTEGER, INTENT(IN)  :: maxelev
  INTEGER, INTENT(IN)  :: ngate
  INTEGER, INTENT(IN)  :: nazim

  INTEGER, INTENT(IN)  :: nyqset
  INTEGER, INTENT(IN)  :: timeset

  INTEGER, INTENT(INOUT) :: kntgate(maxazim,maxelev)
  INTEGER, INTENT(INOUT) :: kntazim(maxelev)
  INTEGER, INTENT(INOUT) :: kntelev

  INTEGER, INTENT(IN)  :: gatesp
  INTEGER, INTENT(IN)  :: rfirstg
  REAL,    INTENT(IN)  :: refchek

  REAL,    INTENT(IN)  :: vnyquist(maxazim)
  INTEGER, INTENT(IN)  :: time(maxazim)
  REAL,    INTENT(IN)  :: azim(maxazim)
  REAL,    INTENT(IN)  :: elev(maxazim)
  REAL,    INTENT(IN)  :: reftilt(maxgate,maxazim)
  REAL,    INTENT(IN)  :: rhvtilt(maxgate,maxazim)
  REAL,    INTENT(IN)  :: zdrtilt(maxgate,maxazim)
  REAL,    INTENT(IN)  :: kdptilt(maxgate,maxazim)

  REAL,    INTENT(INOUT) :: nyqvvol(maxazim,maxelev)
  INTEGER, INTENT(INOUT) :: timevol(maxazim,maxelev)
  REAL,    INTENT(INOUT) :: rngvol(maxgate,maxelev)
  REAL,    INTENT(INOUT) :: azmvol(maxazim,maxelev)
  REAL,    INTENT(INOUT) :: elvvol(maxazim,maxelev)
  REAL,    INTENT(INOUT) :: elvmnvol(maxelev)
  REAL,    INTENT(INOUT) :: refvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(INOUT) :: rhvvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(INOUT) :: zdrvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(INOUT) :: kdpvol(maxgate,maxazim,maxelev)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: igate,jazim,kelev,k
  REAL :: elvsum,elvmean
  LOGICAL :: duplevel
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( nazim > 0 ) THEN
!
! First determine mean elevation angle of this tilt and proper
! vertical array index, kelev to store this data in the volume array
!
    elvsum=0.
    DO jazim=1,nazim
      elvsum=elvsum+elev(jazim)
    END DO
    elvmean=elvsum/float(nazim)

    IF(kntelev > 0) THEN

      duplevel=.false.
      DO kelev=1,kntelev
        IF(abs(elvmnvol(kelev) - elvmean) < 0.2 ) THEN
          duplevel=.true.
          EXIT
        END IF
      END DO

      IF(duplevel) THEN

        WRITE(6,'(a,i4,2f8.2)')  &
        ' volbuild duplicate elevation angle: ',kelev,elvmean,elvmnvol(kelev)

!
! Reset variables on kelev to missing
!
        kntazim(kelev)=0
        DO jazim = 1,maxazim
          kntgate(jazim,kelev)=0
          nyqvvol(jazim,kelev)=0.
          timevol(jazim,kelev)=0.
          azmvol(jazim,kelev)=-999.
          elvvol(jazim,kelev)=-999.
        END DO
        DO igate=1,maxgate
          rngvol(igate,kelev)=-999.
        END DO
        DO jazim=1,maxazim
          DO igate=1,maxgate
            refvol(igate,jazim,kelev)=-999.
            rhvvol(igate,jazim,kelev)=-999.
            zdrvol(igate,jazim,kelev)=-999.
            kdpvol(igate,jazim,kelev)=-999.
          END DO
        END DO
      ELSE

        DO kelev=1,kntelev
          IF(elvmnvol(kelev) > elvmean) EXIT
        END DO
!
! When the vertical index is not above all current data in volume,
! shuffle all the current data up one index.
!
        IF( kelev <= kntelev ) THEN
          WRITE(6,'(a,f8.2,a,i4,a,i4)')                                &
            ' volbuild new elev: ',elvmean,' at k index ',kelev,' of ', &
              kntelev,' existing levels'
          DO k=kntelev,kelev,-1
            elvmnvol(k+1)=elvmnvol(k)
            kntazim(k+1)=kntazim(k)
            DO jazim=1,maxazim
              kntgate(jazim,k+1)=kntgate(jazim,k)
              nyqvvol(jazim,k+1)=nyqvvol(jazim,k)
              timevol(jazim,k+1)=timevol(jazim,k)
              azmvol(jazim,k+1)=azmvol(jazim,k)
              elvvol(jazim,k+1)=elvvol(jazim,k)
            END DO
            DO igate=1,maxgate
              rngvol(igate,k+1)=rngvol(igate,k)
            END DO
            DO jazim=1,maxazim
              DO igate=1,maxgate
                refvol(igate,jazim,k+1)=refvol(igate,jazim,k)
                rhvvol(igate,jazim,k+1)=rhvvol(igate,jazim,k)
                zdrvol(igate,jazim,k+1)=zdrvol(igate,jazim,k)
                kdpvol(igate,jazim,k+1)=kdpvol(igate,jazim,k)
              END DO
            END DO
          END DO
!
! Reset variables on kelev to missing
!
          kntazim(kelev)=0
          DO jazim = 1,maxazim
            kntgate(jazim,kelev)=0
            nyqvvol(jazim,kelev)=0.
            timevol(jazim,kelev)=0.
            azmvol(jazim,kelev)=-999.
            elvvol(jazim,kelev)=-999.
          END DO
          DO igate=1,maxgate
            rngvol(igate,kelev)=-999.
          END DO
          DO jazim=1,maxazim
            DO igate=1,maxgate
              refvol(igate,jazim,kelev)=-999.
              rhvvol(igate,jazim,kelev)=-999.
              zdrvol(igate,jazim,kelev)=-999.
              kdpvol(igate,jazim,kelev)=-999.
            END DO
          END DO
        END IF  ! need to insert data level in middle

        WRITE(6,'(a,f8.2)') ' volbuild adding elevation angle: ',elvmean
        kntelev = kntelev + 1

      END IF ! not dupe

    ELSE

      WRITE(6,'(a,f8.2)') ' volbuild first elevation angle to store: ',elvmean
      kntelev = 1
      kelev = 1

    END IF
!
! Store incoming data in vertical level kelev
!
    kntazim(kelev) = nazim
    elvmnvol(kelev) = elvmean
    IF( nyqset > 0 ) THEN
      DO jazim = 1, nazim
        nyqvvol(jazim,kelev)=vnyquist(jazim)
      END DO
    END IF

    DO igate = 1, maxgate
      rngvol(igate,kelev)=rfirstg+(igate-1)*gatesp
    END DO
    WRITE(6,'(a,i9)') ' volbuild ngate:',ngate
    WRITE(6,'(a,f12.1,a,f12.1)')' volbuild min range:',rngvol(1,kelev), &
                                '  max range:',rngvol(ngate,kelev)

    IF( timeset > 0 ) THEN
      DO jazim = 1, nazim
        timevol(jazim,kelev)=time(jazim)
        azmvol(jazim,kelev)=azim(jazim)
        elvvol(jazim,kelev)=elev(jazim)
      END DO
    ELSE
      DO jazim = 1, nazim
        azmvol(jazim,kelev)=azim(jazim)
        elvvol(jazim,kelev)=elev(jazim)
      END DO
    END IF

    DO jazim = 1, nazim
      DO igate=1, ngate
        refvol(igate,jazim,kelev)=reftilt(igate,jazim)
        IF (reftilt(igate,jazim) > refchek)                               &
          kntgate(jazim,kelev)=igate
        rhvvol(igate,jazim,kelev)=rhvtilt(igate,jazim)
        zdrvol(igate,jazim,kelev)=zdrtilt(igate,jazim)
        kdpvol(igate,jazim,kelev)=kdptilt(igate,jazim)
      END DO
    END DO

    WRITE(6,'(a)')  ' Current Mean Elevation Angles in Volume '
    WRITE(6,'(a)')  '  Index  Elevation   N azimuths'
    DO k = 1, kntelev
      WRITE(6,'(1x,i6,f12.2,i10)') k,elvmnvol(k),kntazim(k)
    END DO

  END IF

  RETURN
END SUBROUTINE volbuilddlp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE RGATEXYZ                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rgatexyz(maxvargate,maxgate,maxazim,maxelev,nzsnd,         &
                    rfropt,lvldbg,                                    &
                    kntgat,kntazm,kntelv,                             &
                    radlat,radlon,radarx,radary,radalt,               &
                    rngvol,azmvol,elvvol,                             &
                    zsnd,rfrsnd,                                      &
                    rxvol,ryvol,rzvol,istatus)
!
!-----------------------------------------------------------------------
!
!  Keith Brewster, CAPS
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxvargate
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nzsnd

  INTEGER, INTENT(IN) :: rfropt
  INTEGER, INTENT(IN) :: lvldbg

  INTEGER, INTENT(IN) :: kntgat(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntazm(maxelev)
  INTEGER, INTENT(IN) :: kntelv

  REAL, INTENT(IN)    :: radlat
  REAL, INTENT(IN)    :: radlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radalt

  REAL,    INTENT(IN) :: rngvol(maxvargate,maxelev)
  REAL,    INTENT(IN) :: azmvol(maxazim,maxelev)
  REAL,    INTENT(IN) :: elvvol(maxazim,maxelev)

  REAL,    INTENT(IN) :: zsnd(nzsnd)
  REAL,    INTENT(IN) :: rfrsnd(nzsnd)

  REAL,    INTENT(OUT):: rxvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(OUT):: ryvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(OUT):: rzvol(maxgate,maxazim,maxelev)

  INTEGER, INTENT(OUT):: istatus
!
! Misc local variables
!
  INTEGER :: igate,jazim,kelev
  REAL :: rmapfct,azmrot,xcomp,ycomp,zagl,sfcr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL lattomf(1,1,radlat,rmapfct)
  istatus=0
  DO kelev=1,kntelv
    DO jazim=1,kntazm(kelev)
      CALL ddrotuv(1,radlon,azmvol(jazim,kelev),rmapfct,              &
                     azmrot,xcomp,ycomp)
      DO igate=1,kntgat(jazim,kelev)
        CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,                &
             elvvol(jazim,kelev),rngvol(igate,kelev),zagl,sfcr)
        rxvol(igate,jazim,kelev)=radarx-xcomp*sfcr
        ryvol(igate,jazim,kelev)=radary-ycomp*sfcr
        rzvol(igate,jazim,kelev)=radalt+zagl
      END DO
    END DO
  END DO
END SUBROUTINE rgatexyz
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE RMPSETUP                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rmpsetup(maxrgate,maxvgate,maxgate,maxazim,maxelev,         &
              nx,ny,nxny,nz,nzsnd,                                     &
              rfropt,refchek,velchek,bmwidth,velproc,                  &
              kntrgat,kntrazm,kntrelv,                                 &
              kntvgat,kntvazm,kntvelv,                                 &
              radlat,radlon,radarx,radary,radalt,                      &
              dazim,rngmin,rngmax,                                     &
              rngrvol,azmrvol,elvrvol,                                 &
              rngvvol,azmvvol,elvvvol,                                 &
              refvol,velvol,rxvol,ryvol,rzvol,                         &
              xs,ys,zps,zsnd,rfrsnd,ncolp,ncoltot,                     &
              havdat,icolp,jcolp,xcolp,ycolp,zcolp,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set-up the x,y,z locations of the radar data and determine the
!  distribution of load among MPI processors, in the case of MPI.
!  This is done as a separate step from the remapping to allow for
!  MPI control and data allocation by the driver program before
!  calling remapvol.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  18-Nov-2009
!
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxrgate   Maximum reflectivity gates in a radial
!    maxvgate   Maximum velocity gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxny     nx*ny
!
!    nzsnd     Number of levels in refractivity sounding arrays
!
!    rfropt   Rafractivity option (1: std atmos lapse, 2: avg actual lapse rate)
!
!    varchek  Threshold for checking data, good vs. flagged
!    bmwidth  Radar beamwidth (degrees)
!
!    kntgate   Number of gates in radials
!    kntazim   Number of azimuth radials for each elevation angle
!    kntelev   Number of elevation angles (tilts)
!
!    radlat   Latitude (degrees N) of radar
!    radlon   Longitude (degrees E) of radar
!    radarx   x-location of radar
!    radary   y-location of radar
!    radalt   Elevation (m MSL) of radar
!    dazim    Typical azimuth spacing for this dataset (degrees)
!
!    rngmin   Minimum range (m) of data to use
!             (5000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!    rngvol   Range to gate in velocity 3-D volume
!    azmvol   Azimuth angle in velocity 3-D volume
!    elvvol   Elevation angle in velocity 3-D volume
!    varvol   Radar data 3-D volume
!
!    xs       x coordinate of scalar grid points in physical/comp. space (m)
!    ys       y coordinate of scalar grid points in physical/comp. space (m)
!    zps      Vertical coordinate of scalar grid points in physical space(m)
!
!    zsnd     Heights of levels in refractivity profile
!    rfrsnd   Refractivity profile
!
!  OUTPUT:
!
!    rxvol    x-coordinate at radar data location
!    ryvol    y-coordinate at radar data location
!    rzvol    z-coordinate at radar data location
!
!    havdat   Flag to indicate valid data available at this column
!    icolp    global i index for MPI load sharing
!    jcolp    gloval j index for MPI load sharing
!    xcolp    x grid location for MPI load sharing
!    ycolp    y grid location for MPI load sharing
!    zcolp    grid height at columns
!    ncolp    Number of columns per processor
!    ncoltot  Total number of columns
!
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
  INTEGER, INTENT(IN) :: maxrgate
  INTEGER, INTENT(IN) :: maxvgate
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: nxny
  INTEGER, INTENT(IN) :: nzsnd

  INTEGER, INTENT(IN) :: rfropt

  REAL, INTENT(IN)    :: refchek
  REAL, INTENT(IN)    :: velchek
  REAL, INTENT(IN)    :: bmwidth
  LOGICAL, INTENT(IN) :: velproc

  INTEGER, INTENT(IN) :: kntrgat(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntrazm(maxelev)
  INTEGER, INTENT(IN) :: kntrelv

  INTEGER, INTENT(IN) :: kntvgat(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntvazm(maxelev)
  INTEGER, INTENT(IN) :: kntvelv

  REAL, INTENT(IN)    :: radlat
  REAL, INTENT(IN)    :: radlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radalt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL,    INTENT(IN) :: rngrvol(maxrgate,maxelev)
  REAL,    INTENT(IN) :: azmrvol(maxazim,maxelev)
  REAL,    INTENT(IN) :: elvrvol(maxazim,maxelev)

  REAL,    INTENT(IN) :: rngvvol(maxvgate,maxelev)
  REAL,    INTENT(IN) :: azmvvol(maxazim,maxelev)
  REAL,    INTENT(IN) :: elvvvol(maxazim,maxelev)

  REAL,    INTENT(IN) :: refvol (maxrgate,maxazim,maxelev)
  REAL,    INTENT(IN) :: velvol (maxvgate,maxazim,maxelev)

  REAL,    INTENT(OUT):: rxvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(OUT):: ryvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(OUT):: rzvol(maxgate,maxazim,maxelev)

  REAL,    INTENT(IN) :: xs(nx)
  REAL,    INTENT(IN) :: ys(ny)
  REAL,    INTENT(IN) :: zps(nx,ny,nz)
  REAL,    INTENT(IN) :: zsnd(nzsnd)
  REAL,    INTENT(IN) :: rfrsnd(nzsnd)

  INTEGER, INTENT(OUT) :: ncolp
  INTEGER, INTENT(OUT) :: ncoltot

  LOGICAL, INTENT(OUT) :: havdat(nx,ny)
  INTEGER, INTENT(OUT) :: icolp(nxny)
  INTEGER, INTENT(OUT) :: jcolp(nxny)
  REAL,    INTENT(OUT) :: xcolp(nxny)
  REAL,    INTENT(OUT) :: ycolp(nxny)
  REAL,    INTENT(OUT) :: zcolp(nz,nxny)

  INTEGER, INTENT(OUT) :: istatus
!
! Misc local variables
!
  INTEGER :: k,kelev,jazim,igate
  INTEGER :: i,ii,idat,ibgn,iend,iwdth
  INTEGER :: j,jj,jdat,jbgn,jend,jwdth
  INTEGER :: idxbgn,idxend
  INTEGER :: ioffset,joffset
!
  REAL :: deg2rad,mapfct,azmrot,xcomp,ycomp,zagl,sfcr
  REAL :: bmwmax,sinbmw,sinbmwx,sinbmwy
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=-1
  deg2rad = atan(1.)/45.
  dxinv=1./dx
  dyinv=1./dy
  bmwmax=max(bmwidth,dazim)
  sinbmw=sin(0.5*deg2rad*bmwmax)
  sinbmwx=dxinv*sinbmw
  sinbmwy=dyinv*sinbmw

  havdat=.FALSE.
  icolp=-99
  jcolp=-99
  xcolp=-999.
  ycolp=-999.
  zcolp=-999.
  ncolp=0
  ncoltot=0
!
  CALL lattomf(1,1,radlat,mapfct)

  IF(velproc) THEN
!
! Find x,y,z location of each velocity observation point
!
    CALL rgatexyz(maxvgate,maxgate,maxazim,maxelev,nzsnd,            &
                  rfropt,lvldbg,                                    &
                  kntvgat,kntvazm,kntvelv,                           &
                  radlat,radlon,radarx,radary,radalt,                &
                  rngvvol,azmvvol,elvvvol,                           &
                  zsnd,rfrsnd,                                       &
                  rxvol,ryvol,rzvol,istatus)
!
! Mark the 2d havdat mask array for each datum having valid data.
! Include a halo the width of the beam around each datum
! For simplicity the halo here is equidistant in x,y directions.
!
    DO kelev=1,kntvelv
      DO jazim=1,kntvazm(kelev)
        DO igate=1,kntvgat(jazim,kelev)
          IF(velvol(igate,jazim,kelev) > velchek .AND.                   &
             rngvvol(igate,kelev) >= rngmin .AND.                          &
             rngvvol(igate,kelev) <= rngmax ) THEN
            idat=1+NINT((rxvol(igate,jazim,kelev)-xs(1))*dxinv)
            jdat=1+NINT((ryvol(igate,jazim,kelev)-ys(1))*dyinv)
            IF(idat > 0 .AND. idat <= nx .AND.                           &
               jdat > 0 .AND. jdat <= ny ) THEN
              iwdth=1+INT(sinbmwx*rngvvol(igate,kelev))
              jwdth=1+INT(sinbmwy*rngvvol(igate,kelev))
              ibgn=MAX(1,(idat-iwdth))
              iend=MIN(nx,(idat+iwdth))
              jbgn=MAX(1,(jdat-jwdth))
              jend=MIN(ny,(jdat+jwdth))
              DO jj=jbgn,jend
                DO ii=ibgn,iend
                  havdat(ii,jj)=.TRUE.
                END DO
              END DO
            END IF
          END IF

        END DO
      END DO
    END DO
  END IF ! velproc
!
! Find x,y,z location of each reflectivity observation point
!
  CALL rgatexyz(maxrgate,maxgate,maxazim,maxelev,nzsnd,                &
                    rfropt,lvldbg,                                    &
                    kntrgat,kntrazm,kntrelv,                           &
                    radlat,radlon,radarx,radary,radalt,                &
                    rngrvol,azmrvol,elvrvol,                           &
                    zsnd,rfrsnd,                                       &
                    rxvol,ryvol,rzvol,istatus)
!
! Mark the 2d havdat mask array for each datum having valid data.
! Include a halo the width of the beam around each datum
! For simplicity the halo here is equidistant in x,y directions.
!
  DO kelev=1,kntrelv
    DO jazim=1,kntrazm(kelev)
      DO igate=1,kntrgat(jazim,kelev)
        IF(refvol(igate,jazim,kelev) > refchek .AND.                   &
           rngrvol(igate,kelev) >= rngmin .AND.                          &
           rngrvol(igate,kelev) <= rngmax ) THEN
          idat=1+NINT((rxvol(igate,jazim,kelev)-xs(1))*dxinv)
          jdat=1+NINT((ryvol(igate,jazim,kelev)-ys(1))*dyinv)
          IF(idat > 0 .AND. idat <= nx .AND.                           &
             jdat > 0 .AND. jdat <= ny ) THEN
            iwdth=1+INT(sinbmwx*rngrvol(igate,kelev))
            jwdth=1+INT(sinbmwy*rngrvol(igate,kelev))
            ibgn=MAX(1,(idat-iwdth))
            iend=MIN(nx,(idat+iwdth))
            jbgn=MAX(1,(jdat-jwdth))
            jend=MIN(ny,(jdat+jwdth))
            DO jj=jbgn,jend
              DO ii=ibgn,iend
                havdat(ii,jj)=.TRUE.
              END DO
            END DO
          END IF
        END IF

      END DO
    END DO
  END DO
!
! Set-up MPI Load Sharing arrays
!
  idat=0
  IF(nprocs == 1) THEN
    ioffset=0
    joffset=0
  ELSE
    ioffset=(loc_x-1)*(nx-3)
    joffset=(loc_y-1)*(ny-3)
!   print *, ' myproc: ',myproc,' loc_x: ',loc_x,'  ioffset:',ioffset
!   print *, ' myproc: ',myproc,' loc_y: ',loc_y,'  joffset:',joffset
  END IF

  DO j=2,ny-2
    DO i=2,nx-2
      IF(havdat(i,j)) THEN
        idat=idat+1
        icolp(idat)=i+ioffset
        jcolp(idat)=j+joffset
        xcolp(idat)=xs(i)
        ycolp(idat)=ys(j)
        DO k=1,nz
          zcolp(k,idat)=zps(i,j,k)
        END DO
      END IF
    END DO
  END DO
  ncolp=idat
  WRITE(6,'(a,i6,a,i10,a,i10)')    &
     ' Myproc',myproc,' found ',ncolp,' i,j points with data out of ',nxny

  ncoltot=ncolp
  CALL mptotali(ncoltot)

  istatus=0
  RETURN
  END SUBROUTINE rmpsetup
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE REMAP2D                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,nsort,         &
                   vardump,ivar,rfropt,                                 &
                   varid,varname,varunits,dmpfmt,hdf4cmpr,              &
                   varchek,varmiss,vmedlim,dazlim,iorder,               &
                   sortmin,dsort,                                       &
                   kntgate,kntazim,kntelev,                             &
                   radlat,radlon,radarx,radary,radalt,dazim,            &
                   rngmin,rngmax,                                       &
                   rngvol,azmvol,elvvol,                                &
                   varvol,rxvol,ryvol,                                  &
                   xs,ys,zsnd,rfrsnd,kntbin,colvar,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from the lowest tilt and remap it onto a plane for
!  purposes of display in comparison to the ARPS k=2 reflectivity.
!  Uses a least squares a local quadratic fit to remap the data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!          June, 2002
!
!  MODIFICATION HISTORY:
!  08-Oct-2003  Keith Brewster
!               Modified search logic to search at least jsrchmn (2) radials
!               from estimated closest radial.  Expanded the area included
!               in least sq fit to be 1.5 times width of radial.  Makes for
!               better horizontal continuity with more overlap of data.
!
!  14-Jan-2004  Keith Brewster
!               Modified logic to avoid possibility of overflow of kntbin.
!               kntbin and elvmean arrays are now allocated outside and
!               passed into this routine.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    nzsnd     Number of levels in refractivity sounding arrays
!    nsort   Number of elements in sorting vector for computing median
!
!    vardump  Switch for writing (1) or not writing (0) remapped data to file
!    ivar     Number of variable:
!             1: Reflectivity
!             2: Radial Velocity
!             3: Spectrum Width
!    rfropt   Rafractivity option (1: std atmos lapse, 2: avg actual lapse rate)
!
!    varid    Radar variable ID for diagnostic file writing
!    varname  Name of radar variable for diagnostic file writing
!    varunits Units of radar variable for diagnostic file writing
!
!    varchek  Threshold for checking data, good vs. flagged
!    varmiss  Value to assign to data for missing
!    vmedlim  Threshold limit for median check
!    dazlim   Maximum value of azimuth difference (grid vs data) to accept
!             Generally should be 30 degrees or less for velocity, 360 for refl
!    iorder   Order of polynomial to fit (1: linear, 2: quadratic)
!
!    kntgate   Number of gates in radials
!    kntazim   Number of azimuth radials for each elevation angle
!    kntelev   Number of elevation angles (tilts)
!
!    radlat   Latitude (degrees N) of radar
!    radlon   Longitude (degrees E) of radar
!    radarx   x-location of radar
!    radary   y-location of radar
!    radalt   Elevation (m MSL) of radar
!    dazim    Typical azimuth spacing for this dataset (degrees)
!
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
!
!    zsnd     Heights of levels in refractivity profile
!    rfrsnd   Refractivity profile
!
!  OUTPUT:
!
!    kntbin  Temporary array used for computing the mean
!    elvmean  Temporary array average elevation angle for each tilt
!
!    rxvol    x-coordinate at radar data location
!    ryvol    y-coordinate at radar data location
!
!    colvar  Radar variable remapped to scalar points
!
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
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nzsnd
  INTEGER, INTENT(IN) :: nsort

  INTEGER, INTENT(IN) :: vardump
  INTEGER, INTENT(IN) :: ivar
  INTEGER, INTENT(IN) :: rfropt

  CHARACTER (LEN=6), INTENT(IN) :: varid
  CHARACTER (LEN=20), INTENT(IN) :: varname
  CHARACTER (LEN=20), INTENT(IN) :: varunits
  INTEGER, INTENT(IN) :: dmpfmt
  INTEGER, INTENT(IN) :: hdf4cmpr

  REAL, INTENT(IN)    :: varchek
  REAL, INTENT(IN)    :: varmiss
  REAL, INTENT(IN)    :: vmedlim
  REAL, INTENT(IN)    :: dazlim
  INTEGER, INTENT(IN) :: iorder

  REAL, INTENT(IN)    :: sortmin
  REAL, INTENT(IN)    :: dsort

  INTEGER, INTENT(IN) :: kntgate(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntazim(maxelev)
  INTEGER, INTENT(IN) :: kntelev

  REAL, INTENT(IN)    :: radlat
  REAL, INTENT(IN)    :: radlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radalt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL, INTENT(IN)    :: rngvol(maxgate,maxelev)
  REAL, INTENT(IN)    :: azmvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: varvol(maxgate,maxazim,maxelev)
  REAL, INTENT(OUT)   :: rxvol(maxgate,maxazim,maxelev)
  REAL, INTENT(OUT)   :: ryvol(maxgate,maxazim,maxelev)

  REAL, INTENT(IN)    :: xs(nx)
  REAL, INTENT(IN)    :: ys(ny)
  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)

  INTEGER, INTENT(OUT)   :: kntbin(nsort)
  REAL, INTENT(OUT)   :: colvar(nx,ny)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: n = 6
  REAL :: avar(n,n)
  REAL :: rhsvar(n)
  REAL :: avel(n,n)
  REAL :: rhsvel(n)
  REAL :: sol(n)
  REAL :: work(n,n+1)
  REAL :: work1d(n+1)

  REAL :: array(3,3)
  REAL :: rhsv(3)
  REAL :: solv(3)

  REAL, PARAMETER :: eps = 1.0E-25
  REAL, PARAMETER :: dxlim0=1.1

  INTEGER :: ii,jj,kk,i,j,k,knt,kinbox,jknt,kntall
  INTEGER :: kok,isort,jsort,mid,maxkok
  INTEGER :: kbgn,kend
  INTEGER :: igate,jazim,jazmin,jmirror,jend
  INTEGER :: istatal,istatwrt

! jsrchmn is the minimum number of radials to search to find data
! that are within the distance limits of the grid center.
!
  INTEGER, PARAMETER :: jsrchmn = 2

  REAL :: deg2rad,rad2deg
  REAL :: sortscale
  REAL :: delx,dely,dazimr,daz,azdiff,ff
  REAL :: ddx,ddxy,ddx2,ddy,ddy2,dxthr,dxthr0
  REAL :: azmrot,xcomp,ycomp,mapfct,sfcr,zagl
  REAL :: sum,sum2,sdev,thresh,slrange,elvidat,azimidat,time
  REAL :: varmax,varmin,varavg,varmean,varmed

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  kk=1
  maxkok=0
  deg2rad = atan(1.)/45.
  rad2deg = 1./deg2rad
  dazimr = dxlim0*dazim*deg2rad
  dxthr0=dxlim0*max((xs(3)-xs(2)),(ys(3)-ys(2)))
  sortscale=1./dsort

  CALL lattomf(1,1,radlat,mapfct)

  time=0.
!
  DO j=1,ny
    DO i=1,nx
      colvar(i,j)=varmiss
    END DO
  END DO
!
  DO jazim=1,kntazim(kk)
    CALL ddrotuv(1,radlon,azmvol(jazim,kk),mapfct,                     &
                   azmrot,xcomp,ycomp)
    DO igate=1,kntgate(jazim,kk)
      CALL beamhgtn(nzsnd,zsnd,rfrsnd,radalt,rfropt,                   &
             elvvol(jazim,kk),rngvol(igate,kk),zagl,sfcr)
      rxvol(igate,jazim,kk)=radarx-xcomp*sfcr
      ryvol(igate,jazim,kk)=radary-ycomp*sfcr
    END DO
  END DO
!
  DO j=1,ny-1
    DO i=1,nx-1
      kok=0
      sum=0.
      sum2=0.
      varavg=999999.
      kntbin=0
      delx=xs(i)-radarx
      dely=ys(j)-radary
      slrange=sqrt(delx*delx+dely*dely)
!
      IF( slrange > rngmin .AND. slrange < rngmax ) THEN
        IF(ivar == 1) colvar(i,j)=0.0
        dxthr=max(dxthr0,((slrange+dxthr0)*dazimr))
!
!-----------------------------------------------------------------------
!
! Determine azimuth to this grid cell
!
!-----------------------------------------------------------------------
!
        CALL uvrotdd(1,1,radlon,-delx,-dely,azimidat,ff)

        varmax=-999.
        varmin=999.
        DO jj=1,n
          DO ii=1,n
            avar(ii,jj)=0.
          END DO
        END DO
!
        DO ii=1,n
          rhsvar(ii)=0.
        END DO
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
        azdiff=181.
        jazmin=1
        DO jazim=1,kntazim(kk)
          daz=azmvol(jazim,kk)-azimidat
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          daz=abs(daz)
          IF(daz < azdiff) THEN
            azdiff=daz
            jazmin=jazim
          END IF
        END DO

        jmirror=jazmin+(kntazim(kk)/2)
        IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  First pass, find median, avg, std dev.
!
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
        jend=kntazim(kk)
        IF(jmirror > jazmin) jend=jmirror-1
        jknt=0
        DO jazim=jazmin,jend
          kinbox=0
          jknt=jknt+1
          daz=azmvol(jazim,kk)-azimidat
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          IF(abs(daz) > dazlim) EXIT
          DO igate=1,kntgate(jazim,kk)
            ddx=rxvol(igate,jazim,kk)-xs(i)
            ddy=ryvol(igate,jazim,kk)-ys(j)
!
            IF( rngvol(igate,kk) > rngmin .AND.                    &
                rngvol(igate,kk) < rngmax .AND.                    &
                abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
              kinbox=kinbox+1
              IF(varvol(igate,jazim,kk) > varchek ) THEN
                isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                isort=max(min(isort,nsort),1)
                kntbin(isort)=kntbin(isort)+1
                sum=sum+varvol(igate,jazim,kk)
                sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                kok=kok+1
              END IF   ! data ok
            END IF  ! inside box
          END DO ! igate
          IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
        END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
        IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==kntazim(kk)) THEN
          DO jazim=1,jmirror-1
            kinbox=0
            jknt=jknt+1
            daz=azmvol(jazim,kk)-azimidat
            IF(daz > 180.) daz=daz-360.
            IF(daz < -180.) daz=daz+360.
            IF(abs(daz) > dazlim) EXIT
            DO igate=1,kntgate(jazim,kk)
!
              ddx=rxvol(igate,jazim,kk)-xs(i)
              ddy=ryvol(igate,jazim,kk)-ys(j)

              IF( rngvol(igate,kk) > rngmin .AND.                    &
                  rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                IF(varvol(igate,jazim,kk) > varchek ) THEN
                  isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                  isort=max(min(isort,nsort),1)
                  kntbin(isort)=kntbin(isort)+1
                  sum=sum+varvol(igate,jazim,kk)
                  sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                  kok=kok+1
                END IF
              END IF
            END DO
            IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
        jend= 1
        IF(jmirror < jazmin) jend=jmirror
        jknt=0
        DO jazim=jazmin-1,jend,-1
          kinbox=0
          jknt=jknt+1
          daz=azmvol(jazim,kk)-azimidat
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          IF(abs(daz) > dazlim) EXIT
          DO igate=1,kntgate(jazim,kk)
!
            ddx=rxvol(igate,jazim,kk)-xs(i)
            ddy=ryvol(igate,jazim,kk)-ys(j)

            IF( rngvol(igate,kk) > rngmin .AND.                    &
                rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

              kinbox=kinbox+1

              IF(varvol(igate,jazim,kk) > varchek ) THEN
                isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                isort=max(min(isort,nsort),1)
                kntbin(isort)=kntbin(isort)+1
                sum=sum+varvol(igate,jazim,kk)
                sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                kok=kok+1
              END IF
            END IF
          END DO
          IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
        END DO
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
        IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) THEN
          DO jazim=kntazim(kk),jmirror,-1
            kinbox=0
            jknt=jknt+1
            daz=azmvol(jazim,kk)-azimidat
            IF(daz > 180.) daz=daz-360.
            IF(daz < -180.) daz=daz+360.
            IF(abs(daz) > dazlim) EXIT
            DO igate=1,kntgate(jazim,kk)
!
              ddx=rxvol(igate,jazim,kk)-xs(i)
              ddy=ryvol(igate,jazim,kk)-ys(j)
!
              IF( rngvol(igate,kk) > rngmin .AND.                    &
                  rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                kinbox=kinbox+1
                IF(varvol(igate,jazim,kk) > varchek ) THEN
                  isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                  isort=max(min(isort,nsort),1)
                  kntbin(isort)=kntbin(isort)+1
                  sum=sum+varvol(igate,jazim,kk)
                  sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                  kok=kok+1
                END IF
              END IF
            END DO  ! igate
            IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
          END DO ! jazim
        END IF
!
        IF( kok == 0 ) CYCLE
        varavg=sum/float(kok)
        mid=(kok/2)+1
        kntall=0
        DO isort=1,nsort-1
          kntall=kntall+kntbin(isort)
          IF(kntall >= mid) EXIT
        END DO
        varmed=sortmin+((isort-1)*dsort)
        IF ( kok > 2 ) THEN
          sdev=(sum2-(sum*sum/float(kok)))/float(kok-1)
          sdev=sqrt(max(sdev,0.))
          thresh=max((2.*sdev),vmedlim)
        ELSE
          thresh=vmedlim
        END IF
        maxkok=max(maxkok,kok)
!
!-----------------------------------------------------------------------
!
!  Process data for local quadratic fit
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
        azdiff=181.
        jazmin=1
        DO jazim=1,kntazim(kk)
          daz=azmvol(jazim,kk)-azimidat
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          daz=abs(daz)
          IF(daz < azdiff) THEN
            azdiff=daz
            jazmin=jazim
          END IF
        END DO

        jmirror=jazmin+(kntazim(kk)/2)
        IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
        jend=kntazim(kk)
        IF(jmirror > jazmin) jend=jmirror-1
        jknt=0
        DO jazim=jazmin,jend
          kinbox=0
          jknt=jknt+1
          daz=azmvol(jazim,kk)-azimidat
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          IF(abs(daz) > dazlim) EXIT
          DO igate=1,kntgate(jazim,kk)
!
            ddx=rxvol(igate,jazim,kk)-xs(i)
            ddy=ryvol(igate,jazim,kk)-ys(j)
!
            IF( rngvol(igate,kk) > rngmin .AND.                    &
                rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

              kinbox=kinbox+1
              ddxy=ddx*ddy
              ddx2=ddx*ddx
              ddy2=ddy*ddy

              IF(varvol(igate,jazim,kk) > varchek .AND.  &
                 abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                 varmax=max(varmax,varvol(igate,jazim,kk))
                 varmin=min(varmin,varvol(igate,jazim,kk))
!
                 rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                 rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                 rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                 rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                 rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                 rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
!
                 avar(1,1)=avar(1,1)+1.
                 avar(1,2)=avar(1,2)+ddx
                 avar(1,3)=avar(1,3)+ddy
                 avar(1,4)=avar(1,4)+ddxy
                 avar(1,5)=avar(1,5)+ddx2
                 avar(1,6)=avar(1,6)+ddy2
!
                 avar(2,1)=avar(2,1)+ddx
                 avar(2,2)=avar(2,2)+ddx2
                 avar(2,3)=avar(2,3)+ddx*ddy
                 avar(2,4)=avar(2,4)+ddx*ddxy
                 avar(2,5)=avar(2,5)+ddx*ddx2
                 avar(2,6)=avar(2,6)+ddx*ddy2
!
                 avar(3,1)=avar(3,1)+ddy
                 avar(3,2)=avar(3,2)+ddy*ddx
                 avar(3,3)=avar(3,3)+ddy2
                 avar(3,4)=avar(3,4)+ddy*ddx2
                 avar(3,5)=avar(3,5)+ddy*ddx2
                 avar(3,6)=avar(3,6)+ddy*ddy2
!
                 avar(4,1)=avar(4,1)+ddxy
                 avar(4,2)=avar(4,2)+ddxy*ddx
                 avar(4,3)=avar(4,3)+ddxy*ddy
                 avar(4,4)=avar(4,4)+ddxy*ddxy
                 avar(4,5)=avar(4,5)+ddxy*ddx2
                 avar(4,6)=avar(4,6)+ddxy*ddy2
!
                 avar(5,1)=avar(5,1)+ddx2
                 avar(5,2)=avar(5,2)+ddx2*ddx
                 avar(5,3)=avar(5,3)+ddx2*ddy
                 avar(5,4)=avar(5,4)+ddx2*ddxy
                 avar(5,5)=avar(5,5)+ddx2*ddx2
                 avar(5,6)=avar(5,6)+ddx2*ddy2
!
                 avar(6,1)=avar(6,1)+ddy2
                 avar(6,2)=avar(6,2)+ddy2*ddx
                 avar(6,3)=avar(6,3)+ddy2*ddy
                 avar(6,4)=avar(6,4)+ddy2*ddxy
                 avar(6,5)=avar(6,5)+ddy2*ddx2
                 avar(6,6)=avar(6,6)+ddy2*ddy2
!
               END IF
!
             END IF
           END DO  ! igate
           IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
         END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
         IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==kntazim(kk)) THEN
           DO jazim=1,jmirror-1
             kinbox=0
             jknt=jknt+1
             daz=azmvol(jazim,kk)-azimidat
             IF(daz > 180.) daz=daz-360.
             IF(daz < -180.) daz=daz+360.
             IF(abs(daz) > dazlim) EXIT
             DO igate=1,kntgate(jazim,kk)
!
               ddx=rxvol(igate,jazim,kk)-xs(i)
               ddy=ryvol(igate,jazim,kk)-ys(j)
!
               IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                 kinbox=kinbox+1
                 ddxy=ddx*ddy
                 ddx2=ddx*ddx
                 ddy2=ddy*ddy

                 IF(varvol(igate,jazim,kk) > varchek .AND.             &
                    abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                   varmax=max(varmax,varvol(igate,jazim,kk))
                   varmin=min(varmin,varvol(igate,jazim,kk))
!
                   rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                   rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                   rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                   rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                   rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                   rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
!
                   avar(1,1)=avar(1,1)+1.
                   avar(1,2)=avar(1,2)+ddx
                   avar(1,3)=avar(1,3)+ddy
                   avar(1,4)=avar(1,4)+ddxy
                   avar(1,5)=avar(1,5)+ddx2
                   avar(1,6)=avar(1,6)+ddy2
!
                   avar(2,1)=avar(2,1)+ddx
                   avar(2,2)=avar(2,2)+ddx2
                   avar(2,3)=avar(2,3)+ddx*ddy
                   avar(2,4)=avar(2,4)+ddx*ddxy
                   avar(2,5)=avar(2,5)+ddx*ddx2
                   avar(2,6)=avar(2,6)+ddx*ddy2
!
                   avar(3,1)=avar(3,1)+ddy
                   avar(3,2)=avar(3,2)+ddy*ddx
                   avar(3,3)=avar(3,3)+ddy2
                   avar(3,4)=avar(3,4)+ddy*ddxy
                   avar(3,5)=avar(3,5)+ddy*ddx2
                   avar(3,6)=avar(3,6)+ddy*ddy2
!
                   avar(5,1)=avar(5,1)+ddxy
                   avar(5,2)=avar(5,2)+ddxy*ddx
                   avar(5,3)=avar(5,3)+ddxy*ddy
                   avar(5,4)=avar(5,4)+ddxy*ddxy
                   avar(5,5)=avar(5,5)+ddxy*ddx2
                   avar(5,6)=avar(5,6)+ddxy*ddy2

                   avar(6,1)=avar(6,1)+ddx2
                   avar(6,2)=avar(6,2)+ddx2*ddx
                   avar(6,3)=avar(6,3)+ddx2*ddy
                   avar(6,4)=avar(6,4)+ddx2*ddxy
                   avar(6,5)=avar(6,5)+ddx2*ddx2
                   avar(6,6)=avar(6,6)+ddx2*ddy2
!
                 END IF
!
               END IF
             END DO  ! igate
             IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
           END DO ! jazim
         END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
         jend= 1
         IF(jmirror < jazmin) jend=jmirror
         jknt=0
         DO jazim=jazmin-1,jend,-1
           kinbox=0
           jknt=jknt+1
           daz=azmvol(jazim,kk)-azimidat
           IF(daz > 180.) daz=daz-360.
           IF(daz < -180.) daz=daz+360.
           IF(abs(daz) > dazlim) EXIT
           DO igate=1,kntgate(jazim,kk)
!
             ddx=rxvol(igate,jazim,kk)-xs(i)
             ddy=ryvol(igate,jazim,kk)-ys(j)
!
             IF( rngvol(igate,kk) > rngmin .AND.                    &
                 rngvol(igate,kk) < rngmax .AND.                    &
                 abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

               kinbox=kinbox+1
               ddxy=ddx*ddy
               ddx2=ddx*ddx
               ddy2=ddy*ddy

               IF(varvol(igate,jazim,kk) > varchek .AND.             &
                  abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                 varmax=max(varmax,varvol(igate,jazim,kk))
                 varmin=min(varmin,varvol(igate,jazim,kk))
!
                 rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                 rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                 rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                 rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                 rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                 rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
!
                 avar(1,1)=avar(1,1)+1.
                 avar(1,2)=avar(1,2)+ddx
                 avar(1,3)=avar(1,3)+ddy
                 avar(1,4)=avar(1,4)+ddxy
                 avar(1,5)=avar(1,5)+ddx2
                 avar(1,6)=avar(1,6)+ddy2
!
                 avar(2,1)=avar(2,1)+ddx
                 avar(2,2)=avar(2,2)+ddx2
                 avar(2,3)=avar(2,3)+ddx*ddy
                 avar(2,4)=avar(2,4)+ddx*ddxy
                 avar(2,5)=avar(2,5)+ddx*ddx2
                 avar(2,6)=avar(2,6)+ddx*ddy2
!
                 avar(3,1)=avar(3,1)+ddy
                 avar(3,2)=avar(3,2)+ddy*ddx
                 avar(3,3)=avar(3,3)+ddy2
                 avar(3,4)=avar(3,4)+ddy*ddxy
                 avar(3,5)=avar(3,5)+ddy*ddx2
                 avar(3,6)=avar(3,6)+ddy*ddy2
!
                 avar(4,1)=avar(4,1)+ddxy
                 avar(4,2)=avar(4,2)+ddxy*ddx
                 avar(4,3)=avar(4,3)+ddxy*ddy
                 avar(4,4)=avar(4,4)+ddxy*ddxy
                 avar(4,5)=avar(4,5)+ddxy*ddx2
                 avar(4,6)=avar(4,6)+ddxy*ddy2
!
                 avar(5,1)=avar(5,1)+ddx2
                 avar(5,2)=avar(5,2)+ddx2*ddx
                 avar(5,3)=avar(5,3)+ddx2*ddy
                 avar(5,4)=avar(5,4)+ddx2*ddxy
                 avar(5,5)=avar(5,5)+ddx2*ddx2
                 avar(5,6)=avar(5,6)+ddx2*ddy2
!
                 avar(6,1)=avar(6,1)+ddy2
                 avar(6,2)=avar(6,2)+ddy2*ddx
                 avar(6,3)=avar(6,3)+ddy2*ddy
                 avar(6,4)=avar(6,4)+ddy2*ddxy
                 avar(6,5)=avar(6,5)+ddy2*ddx2
                 avar(6,6)=avar(6,6)+ddy2*ddy2
!
               END IF
!
             END IF
           END DO  ! igate
           IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
         END DO ! jazim
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
         IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) THEN
         DO jazim=kntazim(kk),jmirror,-1
           kinbox=0
           jknt=jknt+1
           daz=azmvol(jazim,kk)-azimidat
           IF(daz > 180.) daz=daz-360.
           IF(daz < -180.) daz=daz+360.
           IF(abs(daz) > dazlim) EXIT
           DO igate=1,kntgate(jazim,kk)
!
             ddx=rxvol(igate,jazim,kk)-xs(i)
             ddy=ryvol(igate,jazim,kk)-ys(j)
!
             IF( rngvol(igate,kk) > rngmin .AND.                    &
                 rngvol(igate,kk) < rngmax .AND.                    &
                 abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

               kinbox=kinbox+1
               ddxy=ddx*ddy
               ddx2=ddx*ddx
               ddy2=ddy*ddy

               IF(varvol(igate,jazim,kk) > varchek .AND.             &
                  abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                  varmax=max(varmax,varvol(igate,jazim,kk))
                  varmin=min(varmin,varvol(igate,jazim,kk))
!
                  rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                  rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                  rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                  rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                  rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                  rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
!
                  avar(1,1)=avar(1,1)+1.
                  avar(1,2)=avar(1,2)+ddx
                  avar(1,3)=avar(1,3)+ddy
                  avar(1,4)=avar(1,4)+ddxy
                  avar(1,5)=avar(1,5)+ddx2
                  avar(1,6)=avar(1,6)+ddy2
!
                  avar(2,1)=avar(2,1)+ddx
                  avar(2,2)=avar(2,2)+ddx2
                  avar(2,3)=avar(2,3)+ddx*ddy
                  avar(2,4)=avar(2,4)+ddx*ddxy
                  avar(2,5)=avar(2,5)+ddx*ddx2
                  avar(2,6)=avar(2,6)+ddx*ddy2
!
                  avar(3,1)=avar(3,1)+ddy
                  avar(3,2)=avar(3,2)+ddy*ddx
                  avar(3,3)=avar(3,3)+ddy2
                  avar(3,4)=avar(3,4)+ddy*ddxy
                  avar(3,5)=avar(3,5)+ddy*ddx2
                  avar(3,6)=avar(3,6)+ddy*ddy2
!
                  avar(5,1)=avar(5,1)+ddx2
                  avar(5,2)=avar(5,2)+ddx2*ddx
                  avar(5,3)=avar(5,3)+ddx2*ddy
                  avar(5,4)=avar(5,4)+ddx2*ddxy
                  avar(5,5)=avar(5,5)+ddx2*ddx2
                  avar(5,6)=avar(5,6)+ddx2*ddy2
!
                  avar(6,1)=avar(6,1)+ddx2
                  avar(6,2)=avar(6,2)+ddx2*ddx
                  avar(6,3)=avar(6,3)+ddx2*ddy
                  avar(6,4)=avar(6,4)+ddx2*ddxy
                  avar(6,5)=avar(6,5)+ddx2*ddx2
                  avar(6,6)=avar(6,6)+ddx2*ddy2
!
                END IF

              END IF
            END DO  ! igate
            IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
          END DO ! jazim
        END IF
!
!-----------------------------------------------------------------------
!
!   Solve for variable at grid point
!
!-----------------------------------------------------------------------
!
        knt=nint(avar(1,1))
        IF ( iorder > 1 .and. knt > 7 ) THEN
          varmean=rhsvar(1)/avar(1,1)
          CALL GJELIM(n,avar,rhsvar,sol,work,work1d,eps,istatus)
          colvar(i,j)=min(varmax,max(varmin,sol(1)))
        ELSE IF ( iorder > 0 .and. knt > 5 ) THEN
          DO jj=1,3
            DO ii=1,3
              array(ii,jj)=avar(ii,jj)
            END DO
          END DO
          DO ii=1,3
            rhsv(ii)=rhsvar(ii)
          END DO
          CALL GJELIM(3,array,rhsv,solv,work,work1d,eps,istatus)
          colvar(i,j)=min(varmax,max(varmin,solv(1)))
        ELSE IF ( knt > 0 ) THEN
          varmean=rhsvar(1)/avar(1,1)
          colvar(i,j)=varmean
        END IF

      END IF

    END DO   ! i loop
  END DO   ! j loop

  IF( vardump == 1 ) THEN
    CALL wrtvar2(nx,ny,1,colvar,varid,varname,varunits,              &
             curtim,runname,dirname,dmpfmt,hdf4cmpr,1,istatwrt)
  END IF

  RETURN
END SUBROUTINE remap2d

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE GJELIM                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE gjelim(n,a,rhs,sol,work,work1d,eps,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Solve an N x N system of equations.       [a][sol]=[rhs]
!    Using Gauss-Jordan with array pivoting.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  09/26/00
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    n        Matrix order
!    a        Array
!    rhs      Right-hand side vector of matrix equation
!    eps      Error checking threshold
!
!  OUTPUT:
!    sol      Solution vector
!    istatus  Status indicator
!
!  WORK SPACE:
!    work     Work Array
!    work1d   Work vector
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
  INTEGER, INTENT(IN)  :: n
  REAL,    INTENT(IN)  :: a(n,n)
  REAL,    INTENT(IN)  :: rhs(n)
  REAL,    INTENT(OUT) :: sol(n)
  REAL,    INTENT(OUT) :: work(n,n+1)
  REAL,    INTENT(OUT) :: work1d(n+1)
  REAL,    INTENT(IN)  :: eps
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  REAL    :: pivot,const
  INTEGER :: np1
  INTEGER :: i,j,k,m
!
!-----------------------------------------------------------------------
!
! Initialize the work array
! First set all elements to zero.
! Fill nxn with elements from input a
! Fill last column with RHS vector.
!
!-----------------------------------------------------------------------
!
  np1=n+1

  DO j=1, np1
    DO i=1, n
      work(i,j)=0.0
    END DO
  END DO

  DO j=1, n
    DO i=1, n
      work(i,j)=a(i,j)
    END DO
  END DO

  DO i=1,n
    work(i,np1)=rhs(i)
  END DO

  DO j=1, n
!
!-----------------------------------------------------------------------
!
! Find largest element in column j
!
!-----------------------------------------------------------------------
!
    m=j
    pivot=ABS(work(m,j))
    DO i=j+1,n
      IF(ABS(work(i,j)) > pivot ) THEN
        m=i
        pivot=ABS(work(m,j))
      END IF
    END DO
!
!-----------------------------------------------------------------------
!
! Error trapping
!
!-----------------------------------------------------------------------
!
    IF( pivot < eps ) THEN
      DO i=1, n
        sol(i)=0.
      END DO
      istatus=-1
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
! Swap rows
!
!-----------------------------------------------------------------------
!
    IF(m /= j) THEN
      DO k=1,np1
        work1d(k)=work(j,k)
      END DO
      DO k=1,np1
        work(j,k)=work(m,k)
        work(m,k)=work1d(k)
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
! Normalize Row
!
!-----------------------------------------------------------------------
!
    const=1./work(j,j)
    DO k=1,np1
      work(j,k)=const*work(j,k)
    END DO
    work(j,j)=1.0
!
!-----------------------------------------------------------------------
!
! Elimination
!
!-----------------------------------------------------------------------
!
    DO i=1,n
      IF ( i /= j ) THEN
        const=work(i,j)
        DO k=1,np1
          work(i,k)=work(i,k)-const*work(j,k)
        END DO
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! Transfer last column to sol vector
!
!-----------------------------------------------------------------------
!
  DO i=1,n
    sol(i)=work(i,n+1)
  END DO
  istatus = 1

  RETURN
END SUBROUTINE gjelim
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE BSCANPRT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE bscanprt(maxgate,maxazim,ngate,nazim,                      &
                      gatesp,rfirstg,varchek,                         &
                      azim,elev,vartilt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Print a B-scan (azim-range sector) diagnostic from a single
!  tilt of radar data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS 2005/05/30
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!
!    ngate     Number of gates in each radial
!    nazim     Number of radials in each tilt
!
!    kntgate   Number of gates in radials
!    kntazim   Number of azimuth radials for each elevation angle
!
!    gatesp    Gate spacing
!    rfirstg   Range (m) to first gate
!    varchek   Threshold value to determine good vs. flagged data
!
!    azim      Azimuth angles
!    elev      Elevation angles
!    vartilt   Radar variable (reflectivity, velocity or spectrum width)
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

  INTEGER, INTENT(IN)  :: maxgate
  INTEGER, INTENT(IN)  :: maxazim
  INTEGER, INTENT(IN)  :: ngate
  INTEGER, INTENT(IN)  :: nazim

  INTEGER, INTENT(IN)  :: gatesp
  INTEGER, INTENT(IN)  :: rfirstg
  REAL, INTENT(IN)     :: varchek

  REAL, INTENT(IN)     :: azim(maxazim)
  REAL, INTENT(IN)     :: elev(maxazim)
  REAL, INTENT(IN)     :: vartilt(maxgate,maxazim)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER,PARAMETER :: nazprt=15
  CHARACTER (LEN=7) :: varprt(nazprt)
  REAL :: azmprt(nazprt)
  REAL :: elvprt(nazprt)
  REAL,PARAMETER :: azmbgn=315.
  REAL,PARAMETER :: rkmbgn=5.
  REAL,PARAMETER :: rkmend=20.
  REAL :: range
  INTEGER :: ngatprt
  INTEGER :: ibgngate,iendgate,j,jj,jp,jbgn,jend,jmax
  INTEGER :: igate,jazim
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ngatprt=INT((1000.*(rkmend-rkmbgn))/gatesp)+1
  ibgngate=((1000.*rkmbgn)-rfirstg)/gatesp
  ibgngate=max(ibgngate,1)
  iendgate=min((ibgngate+ngatprt),ngate)
  DO jazim = 1, nazim
    IF(abs(azim(jazim)-azmbgn) < 1.2 ) THEN
      jbgn=jazim
      jend=min((jazim+nazprt-1),nazim)
      jend=min(jend,nazim)
      jj=0
      jmax=0
      DO j=jbgn,jend
        jp=j
        IF(j > nazim ) jp=j-nazim
        jj=jj+1
        jmax=jj
        azmprt(jj)=azim(jp)
        elvprt(jj)=elev(jp)
      END DO
      WRITE(51,'(//3x,a,15f7.2)') ' Elev: ',(elvprt(j),j=1,jmax)
      WRITE(51,'(3x,a,15f7.1)') ' Azim: ',(azmprt(j),j=1,jmax)

      DO igate = iendgate,ibgngate,-1
        jj=0
        jmax=0
        DO j=jbgn,jend
          jp=j
          IF(j > nazim ) jp=j-nazim
          jj=jj+1
          jmax=jj
          IF(vartilt(igate,jp) > varchek) THEN
            WRITE(varprt(jj),'(f7.1)') vartilt(igate,jp)
          ELSE
            varprt(jj)='     .  '
          END IF
        END DO

        range=0.001*(rfirstg+(igate-1)*gatesp)
        WRITE(51,'(1x,f8.2,a,15a7)') range,':',                       &
                                 (varprt(jj),jj=1,jmax)
      END DO
      WRITE(51,'(3x,a,15f7.1)') ' Azim: ',(azmprt(j),j=1,jmax)
      EXIT
    END IF
  END DO

  call flush(51)

  RETURN
END SUBROUTINE bscanprt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE REMAPCOL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE remap3dcol(maxvargate,maxgate,maxazim,maxelev,nxny,nz,    &
                   nzsnd,nsort,ncolp,                                &
                   varfill,ivar,nyqset,timeset,rfropt,               &
                   varchek,varmiss,bmwidth,vmedlim,dazlim,iorder,    &
                   sortmin,dsort,                                    &
                   kntgate,kntazim,kntelev,                          &
                   radlat,radlon,radarx,radary,radalt,dazim,         &
                   rngmin,rngmax,                                    &
                   rngvol,azmvol,elvvol,                             &
                   varvol,timevol,nyqvvol,rxvol,ryvol,rzvol,         &
                   zsnd,rfrsnd,kntbin,elvmean,                       &
                   xcolp,ycolp,zcolp,                                &
                   colvar,coltim,colnyq,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from a radar volume in polar coordinates and map it to
!  ARPS Cartesian terrain-following grid.  Uses a least squares
!  local quadratic fit to remap the data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!
!  MODIFICATION HISTORY:
!  08-Oct-2003  Keith Brewster
!               Modified search logic to search at least jsrchmn (2) radials
!               from estimated closest radial.  Expanded the area included
!               in least sq fit to be 1.5 times width of radial.  Makes for
!               better horizontal continuity with more overlap of data.
!
!  14-Jan-2004  Keith Brewster
!               Modified logic to avoid possibility of overflow of varsort.
!               varsort and elvmean arrays are now allocated outside and
!               passed into this routine.
!
!  18-Mar-2005  Keith Brewster
!               Added option to fill reflectivity below the lowest beam
!               under certain conditions and if varfill is set to TRUE.
!               Added arguments varfill and ivar
!  01-Apr-2008  Keith Brewster
!               Corrected typo on line 1192 related to update of
!               avar(3,5) sum in quadratic fit.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxvargate   Maximum gates in a radial for this variable
!    maxgate   Maximum gates in a radial for all variables
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxny     nx*ny
!
!    nzsnd     Number of levels in refractivity sounding arrays
!    nsort   Number of elements in sorting vector for computing median
!
!    vardump  Switch for writing (1) or not writing (0) remapped data to file
!    varfill  Switch for filling data below lowest elevation angle (1=yes, 0=no)
!             Used only for reflectivity as indicated by ivar
!    ivar     Number of variable 1=reflectivity
!                                2=velocity
!                                3=spectrum width
!                                4=Rho-hv
!                                5=Zdr
!                                6=Kdp
!    nyqset   Switch for setting (1) or not setting (0) the nyquist variable
!    timeset  Switch for setting (1) or not setting (0) the time variable
!    rfropt   Rafractivity option (1: std atmos lapse, 2: avg actual lapse rate)
!
!    varchek  Threshold for checking data, good vs. flagged
!    varmiss  Value to assign to data for missing
!    vmedlim  Threshold limit for median check
!    bmwidth  Radar beamwidth (degrees)
!    dazlim   Maximum value of azimuth difference (grid vs data) to accept
!             Generally should be 30 degrees or less for velocity, 360 for refl
!    iorder   Order of polynomial to fit (1: linear, 2: quadratic)
!
!    kntgate   Number of gates in radials
!    kntazim   Number of azimuth radials for each elevation angle
!    kntelev   Number of elevation angles (tilts)
!
!    radlat   Latitude (degrees N) of radar
!    radlon   Longitude (degrees E) of radar
!    radarx   x-location of radar
!    radary   y-location of radar
!    radalt   Elevation (m MSL) of radar
!    dazim    Typical azimuth spacing for this dataset (degrees)
!
!    rngmin   Minimum range (m) of data to use
!             (5000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!    rngvol   Range to gate in velocity 3-D volume
!    azmvol   Azimuth angle in velocity 3-D volume
!    elvvol   Elevation angle in velocity 3-D volume
!    varvol   Radar data 3-D volume
!    timevol  Time 3-D volume
!    nyqvvol  Nyquist vel 3-D volume
!
!    xs       x coordinate of scalar grid points in physical/comp. space (m)
!    ys       y coordinate of scalar grid points in physical/comp. space (m)
!    zps      Vertical coordinate of scalar grid points in physical space(m)
!
!    zsnd     Heights of levels in refractivity profile
!    rfrsnd   Refractivity profile
!
!  OUTPUT:
!
!    kntbin   Temporary array used for computing the median
!    elvmean  Temporary array average elevation angle for each tilt
!
!    rxvol    x-coordinate at radar data location
!    ryvol    y-coordinate at radar data location
!    rzvol    z-coordinate at radar data location
!
!    colvar  Radar variable remapped to scalar points
!    coltim  Radar time (offset) on scalar points
!    colnyq  Radar Nyquist velocity remapped to scalar points
!
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
  INTEGER, INTENT(IN) :: maxvargate
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nxny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: nzsnd
  INTEGER, INTENT(IN) :: nsort
  INTEGER, INTENT(IN) :: ncolp

  INTEGER, INTENT(IN) :: varfill
  INTEGER, INTENT(IN) :: ivar
  INTEGER, INTENT(IN) :: nyqset
  INTEGER, INTENT(IN) :: timeset
  INTEGER, INTENT(IN) :: rfropt

  REAL, INTENT(IN)    :: varchek
  REAL, INTENT(IN)    :: varmiss
  REAL, INTENT(IN)    :: bmwidth
  REAL, INTENT(IN)    :: vmedlim
  REAL, INTENT(IN)    :: dazlim
  INTEGER, INTENT(IN) :: iorder

  REAL, INTENT(IN)    :: sortmin
  REAL, INTENT(IN)    :: dsort

  INTEGER, INTENT(IN) :: kntgate(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntazim(maxelev)
  INTEGER, INTENT(IN) :: kntelev

  REAL, INTENT(IN)    :: radlat
  REAL, INTENT(IN)    :: radlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radalt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL,    INTENT(IN) :: nyqvvol(maxazim,maxelev)
  INTEGER, INTENT(IN) :: timevol(maxazim,maxelev)
  REAL,    INTENT(IN) :: rngvol(maxvargate,maxelev)
  REAL,    INTENT(IN) :: azmvol(maxazim,maxelev)
  REAL,    INTENT(IN) :: elvvol(maxazim,maxelev)
  REAL,    INTENT(IN) :: varvol(maxvargate,maxazim,maxelev)
  REAL,    INTENT(IN) :: rxvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(IN) :: ryvol(maxgate,maxazim,maxelev)
  REAL,    INTENT(IN) :: rzvol(maxgate,maxazim,maxelev)

  REAL,    INTENT(IN) :: zsnd(nzsnd)
  REAL,    INTENT(IN) :: rfrsnd(nzsnd)

  INTEGER, INTENT(OUT) :: kntbin(nsort)
  REAL,    INTENT(OUT) :: elvmean(maxelev)

  REAL, INTENT(IN) :: xcolp(nxny)
  REAL, INTENT(IN) :: ycolp(nxny)
  REAL, INTENT(IN) :: zcolp(nz,nxny)

  REAL,    INTENT(OUT) :: colvar(nz,nxny)
  REAL,    INTENT(OUT) :: colnyq(nz,nxny)
  REAL,    INTENT(OUT) :: coltim(nz,nxny)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: n = 7
  REAL :: avar(n,n)
  REAL :: rhsvar(n)
  REAL :: sol(n)
  REAL :: work(n,n+1)
  REAL :: work1d(n+1)

  REAL :: array(4,4)
  REAL :: rhsv(4)
  REAL :: solv(4)

  REAL, PARAMETER :: eps = 1.0E-25
  REAL, PARAMETER :: dxlim0 = 1.1

  INTEGER :: ii,jj,kk,i,j,k,knt,kinbox
  INTEGER :: kok,isort,jsort,mid,maxkok
  INTEGER :: kbgn,kend,ncolpr
  INTEGER :: igate,jazim,kelev,jazmin,jmirror,jjend,jend,jknt
  INTEGER :: istatal,istatwrt
  INTEGER :: idat
!
! jsrchmn is the minimum number of radials to search to find data
! that are within the distance limits of the grid center.
!
  INTEGER, PARAMETER :: jsrchmn = 2
!
! dzmaxfl is the max vertical distance (m, radar beam to ground)
!         to fill reflectivity when the fill woption is set
!         (accounting for rain falling to the ground)
! refminfl is the min reflectivity (dBZ) required to execute the filling
!
  REAL, PARAMETER :: dzmaxfl = 2100.
  REAL, PARAMETER :: refminfl = 20.

  INTEGER :: klow,kntall
  REAL :: deg2rad,rad2deg
  REAL :: elevmin,elevmax
  REAL :: sortscale
  REAL :: delx,dely,delz,dazimr,daz,azdiff,ff
  REAL :: ddx,ddxy,ddx2,ddy,ddy2,ddz,dxthr,dxthr0
  REAL :: mapfct,xcomp,ycomp,azmrot,sfcr,zagl,rmax
  REAL :: sum,sum2,sdev,thresh,slrange,azimidat,time
  REAL :: varmax,varmin,varavg,varmean,varmed
  REAL :: sfcrng,bmtop,hgtlow,reflow,bmhgt,elvidat
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  deg2rad = atan(1.)/45.
  rad2deg = 1./deg2rad
  dazimr=dxlim0*dazim*deg2rad
  dxthr0=dxlim0*max(dx,dy)
  sortscale=1./dsort
  dxinv=1./dx
  dyinv=1./dy
  maxkok=0

  CALL lattomf(1,1,radlat,mapfct)

  time=0.
  elevmax=-999.
  elevmin=999.
  DO kelev=1,kntelev
    sum=0.
    DO jazim=1,kntazim(kelev)
      sum=sum+elvvol(jazim,kelev)
      elevmin=min(elevmin,elvvol(jazim,kelev))
      elevmax=max(elevmax,elvvol(jazim,kelev))
    END DO
    IF(sum > 0.) THEN
      elvmean(kelev)=sum/float(kntazim(kelev))
    ELSE
      elvmean(kelev)=2.0
    END IF
  END DO
!
! Allow variation within half a beamwidth
!
  elevmin=elevmin-(0.5*bmwidth)
  elevmax=elevmax+(0.5*bmwidth)
  bmtop=elevmin+bmwidth
!
! Initializations
!
  colvar=varmiss

  IF (lvldbg > 0) WRITE(6,'(a,i6,a,i12)') ' myproc: ',myproc,' ncolp: ',ncolp
  DO idat=1,ncolp
!   IF(mod(idat,5000) == 0 ) print *, ' myproc: ',myproc,' idat:',idat
    delx=xcolp(idat)-radarx
    dely=ycolp(idat)-radary
    sfcrng=sqrt(delx*delx+dely*dely)
    DO k=2,nz-2
        kok=0
        sum=0.
        sum2=0.
        varavg=999999.
        kntbin=0
        delz=zcolp(k,idat)-radalt
        CALL beamelvn(nzsnd,zsnd,rfrsnd,radalt,rfropt,                  &
                      delz,sfcrng,elvidat,slrange)
!
        IF( delz >= 0. .AND. slrange > rngmin .AND. slrange < rngmax) THEN

          dxthr=max(dxthr0,((slrange+dxthr0)*dazimr))
          IF(elvidat > elevmin .AND. elvidat < elevmax ) THEN
!
!-----------------------------------------------------------------------
!
! Determine azimuth to this grid cell
!
!-----------------------------------------------------------------------
!
            CALL uvrotdd(1,1,radlon,-delx,-dely,azimidat,ff)
!
            DO kelev=2,kntelev-1
              IF(elvmean(kelev) > elvidat) EXIT
            END DO
!
            varmax=-999.
            varmin=999.
            DO jj=1,n
            DO ii=1,n
              avar(ii,jj)=0.
            END DO
            END DO
!
            DO ii=1,n
              rhsvar(ii)=0.
            END DO
!
            kbgn=max((kelev-1),1)
            kend=min(kelev,kntelev)

            DO kk=kbgn,kend
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
              azdiff=181.
              jazmin=1
              DO jazim=1,kntazim(kk)
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                daz=abs(daz)
                IF(daz < azdiff) THEN
                  azdiff=daz
                  jazmin=jazim
                END IF
              END DO

              IF(nyqset  > 0) colnyq(k,idat) = nyqvvol(jazmin,kk)
              IF(timeset > 0) coltim(k,idat) =                         &
                                 float(timevol(jazmin,kk))

              jmirror=jazmin+(kntazim(kk)/2)
              IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  First pass, find median, avg, std dev.
!
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
              jend=kntazim(kk)
              IF(jmirror > jazmin) jend=jmirror-1
              jknt=0
              DO jazim=jazmin,jend
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)
!
                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                    kinbox=kinbox+1
                    IF(varvol(igate,jazim,kk) > varchek ) THEN
                      isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                      isort=max(min(isort,nsort),1)
                      kntbin(isort)=kntbin(isort)+1
                      sum=sum+varvol(igate,jazim,kk)
                      sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                      kok=kok+1
                    END IF   ! data ok
                  END IF  ! inside box
                END DO ! igate
                IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
              IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==kntazim(kk)) THEN
              DO jazim=1,jmirror-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)

                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                    IF(varvol(igate,jazim,kk) > varchek ) THEN
                      isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                      isort=max(min(isort,nsort),1)
                      kntbin(isort)=kntbin(isort)+1
                      sum=sum+varvol(igate,jazim,kk)
                      sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                      kok=kok+1
                    END IF
                  END IF
                END DO
                IF(kinbox == 0 .AND. jknt > jsrchmn ) EXIT
              END DO
              END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
              jend= 1
              IF(jmirror < jazmin) jend=jmirror
              jknt=0
              DO jazim=jazmin-1,jend,-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)

                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1

                    IF(varvol(igate,jazim,kk) > varchek ) THEN
                      isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                      isort=max(min(isort,nsort),1)
                      kntbin(isort)=kntbin(isort)+1
                      sum=sum+varvol(igate,jazim,kk)
                      sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                      kok=kok+1
                    END IF
                  END IF
                END DO
                IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
              IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) THEN
              DO jazim=kntazim(kk),jmirror,-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)
!
                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                    kinbox=kinbox+1
                    IF(varvol(igate,jazim,kk) > varchek ) THEN
                      isort=1+NINT((varvol(igate,jazim,kk)-sortmin)*sortscale)
                      isort=max(min(isort,nsort),1)
                      kntbin(isort)=kntbin(isort)+1
                      sum=sum+varvol(igate,jazim,kk)
                      sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                      kok=kok+1
                    END IF
                  END IF
                END DO  ! igate
                IF(kinbox == 0 .and. jknt > jsrchmn) EXIT
              END DO ! jazim
              END IF
            END DO ! kk

            IF( kok == 0 ) CYCLE
            varavg=sum/float(kok)
            mid=(kok/2)+1
            kntall=0
            DO isort=1,nsort-1
              kntall=kntall+kntbin(isort)
              IF(kntall >= mid) EXIT
            END DO
            varmed=sortmin+((isort-1)*dsort)
            IF ( kok > 2 ) THEN
              sdev=(sum2-(sum*sum/float(kok)))/float(kok-1)
              sdev=sqrt(max(sdev,0.))
              thresh=max((2.*sdev),vmedlim)
            ELSE
              thresh=vmedlim
            END IF
            maxkok=max(maxkok,kok)
!
!-----------------------------------------------------------------------
!
!  Process data for local quadratic fit
!
!-----------------------------------------------------------------------
!
            DO kk=kbgn,kend
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
              azdiff=181.
              jazmin=1
              DO jazim=1,kntazim(kk)
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                daz=abs(daz)
                IF(daz < azdiff) THEN
                  azdiff=daz
                  jazmin=jazim
                END IF
              END DO

              jmirror=jazmin+(kntazim(kk)/2)
              IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
              jend=kntazim(kk)
              IF(jmirror > jazmin) jend=jmirror-1
              jknt=0
              DO jazim=jazmin,jend
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)

                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)
                  ddz=rzvol(igate,jazim,kk)-zcolp(k,idat)

                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1
                    ddxy=ddx*ddy
                    ddx2=ddx*ddx
                    ddy2=ddy*ddy

                    IF(varvol(igate,jazim,kk) > varchek .AND.  &
                       abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN

                      varmax=max(varmax,varvol(igate,jazim,kk))
                      varmin=min(varmin,varvol(igate,jazim,kk))

                      rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                      rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                      rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                      rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddz
                      rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddxy
                      rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddx2
                      rhsvar(7)=rhsvar(7)+varvol(igate,jazim,kk)*ddy2

                      avar(1,1)=avar(1,1)+1.
                      avar(1,2)=avar(1,2)+ddx
                      avar(1,3)=avar(1,3)+ddy
                      avar(1,4)=avar(1,4)+ddz
                      avar(1,5)=avar(1,5)+ddxy
                      avar(1,6)=avar(1,6)+ddx2
                      avar(1,7)=avar(1,7)+ddy2

                      avar(2,1)=avar(2,1)+ddx
                      avar(2,2)=avar(2,2)+ddx2
                      avar(2,3)=avar(2,3)+ddx*ddy
                      avar(2,4)=avar(2,4)+ddx*ddz
                      avar(2,5)=avar(2,5)+ddx*ddxy
                      avar(2,6)=avar(2,6)+ddx*ddx2
                      avar(2,7)=avar(2,7)+ddx*ddy2

                      avar(3,1)=avar(3,1)+ddy
                      avar(3,2)=avar(3,2)+ddy*ddx
                      avar(3,3)=avar(3,3)+ddy2
                      avar(3,4)=avar(3,4)+ddy*ddz
                      avar(3,5)=avar(3,5)+ddy*ddxy
                      avar(3,6)=avar(3,6)+ddy*ddx2
                      avar(3,7)=avar(3,7)+ddy*ddy2

                      avar(4,1)=avar(4,1)+ddz
                      avar(4,2)=avar(4,2)+ddz*ddx
                      avar(4,3)=avar(4,3)+ddz*ddy
                      avar(4,4)=avar(4,4)+ddz*ddz
                      avar(4,5)=avar(4,5)+ddz*ddxy
                      avar(4,6)=avar(4,6)+ddz*ddx2
                      avar(4,7)=avar(4,7)+ddz*ddy2

                      avar(5,1)=avar(5,1)+ddxy
                      avar(5,2)=avar(5,2)+ddxy*ddx
                      avar(5,3)=avar(5,3)+ddxy*ddy
                      avar(5,4)=avar(5,4)+ddxy*ddz
                      avar(5,5)=avar(5,5)+ddxy*ddxy
                      avar(5,6)=avar(5,6)+ddxy*ddx2
                      avar(5,7)=avar(5,7)+ddxy*ddy2

                      avar(6,1)=avar(6,1)+ddx2
                      avar(6,2)=avar(6,2)+ddx2*ddx
                      avar(6,3)=avar(6,3)+ddx2*ddy
                      avar(6,4)=avar(6,4)+ddx2*ddz
                      avar(6,5)=avar(6,5)+ddx2*ddxy
                      avar(6,6)=avar(6,6)+ddx2*ddx2
                      avar(6,7)=avar(6,7)+ddx2*ddy2

                      avar(7,1)=avar(7,1)+ddy2
                      avar(7,2)=avar(7,2)+ddy2*ddx
                      avar(7,3)=avar(7,3)+ddy2*ddy
                      avar(7,4)=avar(7,4)+ddy2*ddz
                      avar(7,5)=avar(7,5)+ddy2*ddxy
                      avar(7,6)=avar(7,6)+ddy2*ddx2
                      avar(7,7)=avar(7,7)+ddy2*ddy2

                    END IF

                  END IF
                END DO  ! igate
                IF( kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
              IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==kntazim(kk)) THEN
              DO jazim=1,jmirror-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)

                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)
                  ddz=rzvol(igate,jazim,kk)-zcolp(k,idat)

                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1
                    ddxy=ddx*ddy
                    ddx2=ddx*ddx
                    ddy2=ddy*ddy

                    IF(varvol(igate,jazim,kk) > varchek .AND.             &
                       abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN

                      varmax=max(varmax,varvol(igate,jazim,kk))
                      varmin=min(varmin,varvol(igate,jazim,kk))

                      rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                      rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                      rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                      rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddz
                      rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddxy
                      rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddx2
                      rhsvar(7)=rhsvar(7)+varvol(igate,jazim,kk)*ddy2

                      avar(1,1)=avar(1,1)+1.
                      avar(1,2)=avar(1,2)+ddx
                      avar(1,3)=avar(1,3)+ddy
                      avar(1,4)=avar(1,4)+ddz
                      avar(1,5)=avar(1,5)+ddxy
                      avar(1,6)=avar(1,6)+ddx2
                      avar(1,7)=avar(1,7)+ddy2

                      avar(2,1)=avar(2,1)+ddx
                      avar(2,2)=avar(2,2)+ddx2
                      avar(2,3)=avar(2,3)+ddx*ddy
                      avar(2,4)=avar(2,4)+ddx*ddz
                      avar(2,5)=avar(2,5)+ddx*ddxy
                      avar(2,6)=avar(2,6)+ddx*ddx2
                      avar(2,7)=avar(2,7)+ddx*ddy2

                      avar(3,1)=avar(3,1)+ddy
                      avar(3,2)=avar(3,2)+ddy*ddx
                      avar(3,3)=avar(3,3)+ddy2
                      avar(3,4)=avar(3,4)+ddy*ddz
                      avar(3,5)=avar(3,5)+ddy*ddxy
                      avar(3,6)=avar(3,6)+ddy*ddx2
                      avar(3,7)=avar(3,7)+ddy*ddy2

                      avar(4,1)=avar(4,1)+ddz
                      avar(4,2)=avar(4,2)+ddz*ddx
                      avar(4,3)=avar(4,3)+ddz*ddy
                      avar(4,4)=avar(4,4)+ddz*ddz
                      avar(4,5)=avar(4,5)+ddz*ddxy
                      avar(4,6)=avar(4,6)+ddz*ddx2
                      avar(4,7)=avar(4,7)+ddz*ddy2

                      avar(5,1)=avar(5,1)+ddxy
                      avar(5,2)=avar(5,2)+ddxy*ddx
                      avar(5,3)=avar(5,3)+ddxy*ddy
                      avar(5,4)=avar(5,4)+ddxy*ddz
                      avar(5,5)=avar(5,5)+ddxy*ddxy
                      avar(5,6)=avar(5,6)+ddxy*ddx2
                      avar(5,7)=avar(5,7)+ddxy*ddy2

                      avar(6,1)=avar(6,1)+ddx2
                      avar(6,2)=avar(6,2)+ddx2*ddx
                      avar(6,3)=avar(6,3)+ddx2*ddy
                      avar(6,4)=avar(6,4)+ddx2*ddz
                      avar(6,5)=avar(6,5)+ddx2*ddxy
                      avar(6,6)=avar(6,6)+ddx2*ddx2
                      avar(6,7)=avar(6,7)+ddx2*ddy2

                      avar(7,1)=avar(7,1)+ddy2
                      avar(7,2)=avar(7,2)+ddy2*ddx
                      avar(7,3)=avar(7,3)+ddy2*ddy
                      avar(7,4)=avar(7,4)+ddy2*ddz
                      avar(7,5)=avar(7,5)+ddy2*ddxy
                      avar(7,6)=avar(7,6)+ddy2*ddx2
                      avar(7,7)=avar(7,7)+ddy2*ddy2

                    END IF

                  END IF
                END DO  ! igate
                IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
              END DO ! jazim
              END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
              jend= 1
              IF(jmirror < jazmin) jend=jmirror
              jknt=0
              DO jazim=jazmin-1,jend,-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)

                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)
                  ddz=rzvol(igate,jazim,kk)-zcolp(k,idat)

                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1
                    ddxy=ddx*ddy
                    ddx2=ddx*ddx
                    ddy2=ddy*ddy

                    IF(varvol(igate,jazim,kk) > varchek .AND.             &
                       abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN

                      varmax=max(varmax,varvol(igate,jazim,kk))
                      varmin=min(varmin,varvol(igate,jazim,kk))

                      rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                      rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                      rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                      rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddz
                      rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddxy
                      rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddx2
                      rhsvar(7)=rhsvar(7)+varvol(igate,jazim,kk)*ddy2

                      avar(1,1)=avar(1,1)+1.
                      avar(1,2)=avar(1,2)+ddx
                      avar(1,3)=avar(1,3)+ddy
                      avar(1,4)=avar(1,4)+ddz
                      avar(1,5)=avar(1,5)+ddxy
                      avar(1,6)=avar(1,6)+ddx2
                      avar(1,7)=avar(1,7)+ddy2

                      avar(2,1)=avar(2,1)+ddx
                      avar(2,2)=avar(2,2)+ddx2
                      avar(2,3)=avar(2,3)+ddx*ddy
                      avar(2,4)=avar(2,4)+ddx*ddz
                      avar(2,5)=avar(2,5)+ddx*ddxy
                      avar(2,6)=avar(2,6)+ddx*ddx2
                      avar(2,7)=avar(2,7)+ddx*ddy2

                      avar(3,1)=avar(3,1)+ddy
                      avar(3,2)=avar(3,2)+ddy*ddx
                      avar(3,3)=avar(3,3)+ddy2
                      avar(3,4)=avar(3,4)+ddy*ddz
                      avar(3,5)=avar(3,5)+ddy*ddxy
                      avar(3,6)=avar(3,6)+ddy*ddx2
                      avar(3,7)=avar(3,7)+ddy*ddy2

                      avar(4,1)=avar(4,1)+ddz
                      avar(4,2)=avar(4,2)+ddz*ddx
                      avar(4,3)=avar(4,3)+ddz*ddy
                      avar(4,4)=avar(4,4)+ddz*ddz
                      avar(4,5)=avar(4,5)+ddz*ddxy
                      avar(4,6)=avar(4,6)+ddz*ddx2
                      avar(4,7)=avar(4,7)+ddz*ddy2

                      avar(5,1)=avar(5,1)+ddxy
                      avar(5,2)=avar(5,2)+ddxy*ddx
                      avar(5,3)=avar(5,3)+ddxy*ddy
                      avar(5,4)=avar(5,4)+ddxy*ddz
                      avar(5,5)=avar(5,5)+ddxy*ddxy
                      avar(5,6)=avar(5,6)+ddxy*ddx2
                      avar(5,7)=avar(5,7)+ddxy*ddy2

                      avar(6,1)=avar(6,1)+ddx2
                      avar(6,2)=avar(6,2)+ddx2*ddx
                      avar(6,3)=avar(6,3)+ddx2*ddy
                      avar(6,4)=avar(6,4)+ddx2*ddz
                      avar(6,5)=avar(6,5)+ddx2*ddxy
                      avar(6,6)=avar(6,6)+ddx2*ddx2
                      avar(6,7)=avar(6,7)+ddx2*ddy2

                      avar(7,1)=avar(7,1)+ddy2
                      avar(7,2)=avar(7,2)+ddy2*ddx
                      avar(7,3)=avar(7,3)+ddy2*ddy
                      avar(7,4)=avar(7,4)+ddy2*ddz
                      avar(7,5)=avar(7,5)+ddy2*ddxy
                      avar(7,6)=avar(7,6)+ddy2*ddx2
                      avar(7,7)=avar(7,7)+ddy2*ddy2

                    END IF

                  END IF
                END DO  ! igate
                IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
              END DO ! jazim
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
              IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) THEN
              DO jazim=kntazim(kk),jmirror,-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimidat
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)

                  ddx=rxvol(igate,jazim,kk)-xcolp(idat)
                  ddy=ryvol(igate,jazim,kk)-ycolp(idat)
                  ddz=rzvol(igate,jazim,kk)-zcolp(k,idat)

                  IF( rngvol(igate,kk) > rngmin .AND.                   &
                      rngvol(igate,kk) < rngmax .AND.                   &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1
                    ddxy=ddx*ddy
                    ddx2=ddx*ddx
                    ddy2=ddy*ddy

                    IF(varvol(igate,jazim,kk) > varchek .AND.           &
                       abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN

                      varmax=max(varmax,varvol(igate,jazim,kk))
                      varmin=min(varmin,varvol(igate,jazim,kk))

                      rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                      rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                      rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                      rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddz
                      rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddxy
                      rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddx2
                      rhsvar(7)=rhsvar(7)+varvol(igate,jazim,kk)*ddy2

                      avar(1,1)=avar(1,1)+1.
                      avar(1,2)=avar(1,2)+ddx
                      avar(1,3)=avar(1,3)+ddy
                      avar(1,4)=avar(1,4)+ddz
                      avar(1,5)=avar(1,5)+ddxy
                      avar(1,6)=avar(1,6)+ddx2
                      avar(1,7)=avar(1,7)+ddy2

                      avar(2,1)=avar(2,1)+ddx
                      avar(2,2)=avar(2,2)+ddx2
                      avar(2,3)=avar(2,3)+ddx*ddy
                      avar(2,4)=avar(2,4)+ddx*ddz
                      avar(2,5)=avar(2,5)+ddx*ddxy
                      avar(2,6)=avar(2,6)+ddx*ddx2
                      avar(2,7)=avar(2,7)+ddx*ddy2

                      avar(3,1)=avar(3,1)+ddy
                      avar(3,2)=avar(3,2)+ddy*ddx
                      avar(3,3)=avar(3,3)+ddy2
                      avar(3,4)=avar(3,4)+ddy*ddz
                      avar(3,5)=avar(3,5)+ddy*ddxy
                      avar(3,6)=avar(3,6)+ddy*ddx2
                      avar(3,7)=avar(3,7)+ddy*ddy2

                      avar(4,1)=avar(4,1)+ddz
                      avar(4,2)=avar(4,2)+ddz*ddx
                      avar(4,3)=avar(4,3)+ddz*ddy
                      avar(4,4)=avar(4,4)+ddz*ddz
                      avar(4,5)=avar(4,5)+ddz*ddxy
                      avar(4,6)=avar(4,6)+ddz*ddx2
                      avar(4,7)=avar(4,7)+ddz*ddy2

                      avar(5,1)=avar(5,1)+ddx2
                      avar(5,2)=avar(5,2)+ddx2*ddx
                      avar(5,3)=avar(5,3)+ddx2*ddy
                      avar(5,4)=avar(5,4)+ddx2*ddz
                      avar(5,5)=avar(5,5)+ddx2*ddxy
                      avar(5,6)=avar(5,6)+ddx2*ddx2
                      avar(5,7)=avar(5,7)+ddx2*ddy2

                      avar(6,1)=avar(6,1)+ddx2
                      avar(6,2)=avar(6,2)+ddx2*ddx
                      avar(6,3)=avar(6,3)+ddx2*ddy
                      avar(6,4)=avar(6,4)+ddx2*ddz
                      avar(6,5)=avar(6,5)+ddx2*ddxy
                      avar(6,6)=avar(6,6)+ddx2*ddx2
                      avar(6,7)=avar(6,7)+ddx2*ddy2

                      avar(7,1)=avar(7,1)+ddy2
                      avar(7,2)=avar(7,2)+ddy2*ddx
                      avar(7,3)=avar(7,3)+ddy2*ddy
                      avar(7,4)=avar(7,4)+ddy2*ddz
                      avar(7,5)=avar(7,5)+ddy2*ddxy
                      avar(7,6)=avar(7,6)+ddy2*ddx2
                      avar(7,7)=avar(7,7)+ddy2*ddy2

                    END IF

                  END IF
                END DO  ! igate
                IF( kinbox == 0 .AND. jknt > jsrchmn ) EXIT
              END DO ! jazim
              END IF
!
            END DO ! kk
!
!-----------------------------------------------------------------------
!
!   Solve for variable at grid point
!
!-----------------------------------------------------------------------
!
            knt=nint(avar(1,1))
            IF ( iorder > 1 .and. knt > 7 ) THEN
              varmean=rhsvar(1)/avar(1,1)
              CALL GJELIM(n,avar,rhsvar,sol,work,work1d,eps,istatus)
              colvar(k,idat)=min(varmax,max(varmin,sol(1)))
!             write(6,'(3i3,i7,4f7.1)') i,j,k,knt,                    &
!                    varmin,varmean,varmax,colvar(k,idat)
            ELSE IF ( iorder > 0 .and. knt > 5 ) THEN
              DO jj=1,4
                DO ii=1,4
                  array(ii,jj)=avar(ii,jj)
                END DO
              END DO
              DO ii=1,4
                rhsv(ii)=rhsvar(ii)
              END DO
              CALL GJELIM(4,array,rhsv,solv,work,work1d,eps,istatus)
              colvar(k,idat)=min(varmax,max(varmin,solv(1)))
!             IF(colvar(k,idat) > 50.) THEN
!               write(6,'(a,3i6,4f7.1)') 'LARGE Refl:',k,idat,knt, &
!                   varmin,varmean,varmax,colvar(k,idat)
!               write(6,'(5x,a,3f10.2)') 'elv,azim,rngkm:', &
!                   elvidat,azimidat,(slrange*0.001)
!             END IF
            ELSE IF ( knt > 0 ) THEN
              varmean=rhsvar(1)/avar(1,1)
              colvar(k,idat)=varmean
!             write(6,'(3i3,a,i6,6f7.1)') i,j,k,'*',knt,              &
!                 varmin,varmean,varmax,colvar(k,idat)
            END IF

          END IF

        END IF
      END DO   ! k loop
  END DO   ! idat (i,j) loop

  IF (lvldbg >0) print *, ' myproc:',myproc,' exited the idat loop'

  IF( varfill == 1 .and. ivar == 1) THEN

    IF(myproc == 0)                                                     &
      WRITE(6,'(/a/a,f9.2,a/a,f10.2,a/a,f10.1,a/)')                     &
        ' Filling reflectivity data below lowest beam...',              &
        '   Top of lowest beam: ',bmtop,' degrees',                     &
        '   Min reflectivity for fill: ',refminfl,' dBZ',               &
        '   Max fill distance: ',dzmaxfl,' m.'
!
! Fill in radar reflectivity below the lowest beam height using
! zero gradient assumption.
!
    DO idat=1,ncolp
!
!   find height of lowest elev + half beamwidth
!
      delx=xcolp(idat)-radarx
      dely=ycolp(idat)-radary
      sfcrng=sqrt(delx*delx+dely*dely)
      CALL bmhgtsfr(bmtop,sfcrng,bmhgt)
      klow=nz
      reflow=0.
      hgtlow=1.0E06
      DO k=nz-1,2,-1
        IF(colvar(k,idat) > varmiss) THEN
          klow=k
          reflow=colvar(k,idat)
          hgtlow=zcolp(k,idat)
        END IF
      END DO
      IF(hgtlow < bmhgt .and. reflow > refminfl) THEN
        DO k=2,klow-1
          IF((hgtlow-zcolp(k,idat)) < dzmaxfl) colvar(k,idat)=reflow
        END DO
      END IF
    END DO

  END IF
!
! Save composite reflectivity in k=1 level
!
  IF(ivar==1) THEN
    DO idat=1,ncolp
      rmax=varmiss
      DO k=2,nz-1
        rmax=max(rmax,colvar(k,idat))
      END DO
      colvar(1,idat)=rmax
    END DO
  END IF

  RETURN
END SUBROUTINE remap3dcol
