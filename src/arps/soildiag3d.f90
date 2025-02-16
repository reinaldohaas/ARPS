!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOILDIAG                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE soildiag(nx,ny,nzsoil,x,y,zpsoil,                            &
                    soiltyp,vegtyp,lai,roufns,veg,hterain,              &
                    tsoil,qsoil,wetcanp, qvsfc,                         &
                    usflx,vsflx,ptsflx,qvsflx,                          &
                    windsp,psfc,rhoa,precip,                            &
                    tair,qvair,                                         &
                    cdha,cdqa,cdma,                                     &
                    radsw, rnflx,                                       &
                    shflx,lhflx,gflx,ct,                                &
                    evaprg,evaprtr,evaprr, qvsat,                       &
                    qvsata,f34,tem1soil)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate and print out diagnostics for the surface processes.
!
!-----------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  08/02/94
!
!  MODIFICATION HISTORY:
!
!  10/31/94 (Y. Liu)
!  Re-wrote the subroputine to make it more general.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D array, veg(nx,ny), to the diagnostic printing list
!
!  03/27/1995 (Yuhe Liu)
!  Changed the solor radiation used in the calculation of surface
!  resistence factor F1 from the one at the top of atmosphere to the
!  one at the surface.
!
!  03/27/1995 (Yuhe Liu)
!  Added the surface resistence into the data dumping.
!
!  05/14/2002 (J. Brotzge)
!  Added new variables/modified call statements to allow for multiple
!  soil schemes
!
!-----------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the soil
!
!    soiltyp  Soil type at the horizontal grid points
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!    hterain  The height of surface terrain
!    zpsoil   Depth of soil (m)
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    windsp   Wind speed just above the surface (m/s)
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    qvsflx   Surface flux of moisture (K*kg/(m**2*s))
!
!    psfc     Surface pressure (Pascal)
!    rhoa     Near sfc air density
!    prcpln   Precipitation path length
!    tair     Air temperature (K) near the surface
!    qvair    S.H. near the surface
!    cdha     Surface drag coefficient for heat
!    cdqa     Surface drag coefficient for moisture
!    cdma     Surface drag coefficient for momentum
!    zenith   Zenith
!    radsw    Solar radiation at the top of atmosphere
!    f34     Input coefficient: f3*f4, output surface resistance
!
!  OUTPUT:
!
!    rnflx    Net radiation flus
!    shflx    Sensible heat flux
!    lhflx    Latent heat flux
!    gflx     Diffusive heat flux from ground surface to deep soil
!    rsw      Net short wave radiation to the surface
!    rlwu     Up-ward long wave radiation flux
!    rlwd     Down-ward long wave radiation flux
!    trwv     Transmisivity due to water vapor
!    trsw     Total transmisivity
!    alfz     Zenith dependent albedo
!    alf      Albedo
!    ct       Thermal capacity
!    f34     Surface resistence
!    qvsat    Surface specific humidity at saturation
!    evaprg   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evaprr   Direct evaporation from the fraction delta
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny
  INTEGER :: nzsoil

  REAL :: x(nx)                 ! X-coordinates
  REAL :: y(ny)                 ! Y-coordinates
  REAL :: zpsoil (nx,ny,nzsoil) ! Depth of soil

  INTEGER :: soiltyp(nx,ny)  ! Soil type at the horizontal grid points
  INTEGER :: vegtyp (nx,ny)  ! Vegetation type at the horizontal grid points
  REAL :: lai    (nx,ny)     ! Leaf Area Index
  REAL :: roufns (nx,ny)     ! Surface roughness
  REAL :: veg    (nx,ny)     ! Vegetation fraction
  REAL :: hterain(nx,ny)     ! The height of surface terrain

  REAL :: tsoil (nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny)  ! Canopy water amount

  REAL :: qvsfc  (nx,ny)  ! Effective S.H. at sfc.

  REAL :: usflx  (nx,ny)  ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx  (nx,ny)  ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx (nx,ny)  ! surface flux of heat (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)  ! surface flux of moisture (kg/(m**2*s))

  REAL :: windsp (nx,ny)  ! Wind speed just above the surface (m/s)
  REAL :: psfc   (nx,ny)  ! Surface pressure (Pascal)
  REAL :: rhoa   (nx,ny)  ! Near sfc air density
  REAL :: precip (nx,ny)  ! Precipitation flux reaching the surface
  REAL :: tair   (nx,ny)  ! Air temperature (K) near the surface
  REAL :: qvair  (nx,ny)  ! S.H. near the surface

  REAL :: cdha   (nx,ny)  ! Surface drag coefficient for heat
  REAL :: cdqa   (nx,ny)  ! Surface drag coefficient for moisture
  REAL :: cdma   (nx,ny)  ! Surface drag coefficient for momentum

  REAL :: radsw  (nx,ny)  ! Solar radiation to the surface

  REAL :: rnflx  (nx,ny)  ! Net radiation flus
  REAL :: shflx  (nx,ny)  ! Sensible heat flux
  REAL :: lhflx  (nx,ny)  ! Latent heat flux
  REAL :: gflx   (nx,ny)  ! Diffusive heat flux from ground surface to
                          ! deep soil
  REAL :: ct     (nx,ny)  ! Thermal capacity

  REAL :: evaprg (nx,ny)  ! Evaporation from groud surface
  REAL :: evaprtr(nx,ny)  ! Transpiration of the remaining part
                          ! (1-delta) of leaves
  REAL :: evaprr (nx,ny)  ! Direct evaporation from the fraction delta
  REAL :: qvsat  (nx,ny)  ! Surface specific humidity at saturation

  REAL :: qvsata (nx,ny)  ! qvsat(tair) (kg/kg)
  REAL :: f34    (nx,ny)  ! f34 and surface resistance
  REAL :: tem1soil (nx,ny,nzsoil) ! Temporary array

!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!    solarc     Solar constant (W/m**2)
!    emissg     Emissivity of the ground
!    emissa     Emissivity of the atmosphere
!    sbcst      Stefen-Boltzmann constant
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  LOGICAL :: dumpsfc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dumpsfc = .false.

  IF ( (curtim > tstart) .AND. (nhisdmp > 0) )THEN
    IF ( hdmpopt == 1 )THEN
      dumpsfc = (MOD(nstep,nhisdmp) == 0)
    ELSE IF ( hdmpopt == 2 )THEN
      dumpsfc = (nstep == hdmpstp(nhisdmp))
    END IF
  ELSE IF ( curtim == tstart ) THEN
    dumpsfc = .true.
  END IF

  IF ( .NOT. dumpsfc ) THEN
    RETURN
  END IF

  IF(myproc ==0) WRITE (6,'(a,i8,a,f10.2,a)')                           &
      ' Dump surface and soil-veg variables at time step, ',nstep,      &
      ', model time=',curtim,' (s)'
!
!-----------------------------------------------------------------------
!
!  Calculate the saturated specific humidity, qvsats.
!
!-----------------------------------------------------------------------
!
  IF(mp_opt >0 .AND. joindmp(FINDX_A) > 0) THEN

  CALL wrtjoinflx(nx,ny,nzsoil,x,y,zpsoil,                       &
              soiltyp,vegtyp,lai,roufns,veg,hterain,             &
              tsoil,qsoil,wetcanp, qvsfc,                        &
              usflx,vsflx,ptsflx,qvsflx,                         &
              windsp,psfc,rhoa,precip,                           &
              tair,qvair,                                        &
              cdha,cdqa,cdma,                                    &
              radsw, rnflx,                                      &
              shflx,lhflx,gflx,ct,                               &
              evaprg,evaprtr,evaprr,qvsat,                       &
              qvsata, f34,tem1soil)
  ELSE

  CALL wrtflx(nx,ny,nzsoil,x,y,zpsoil,                           &
              soiltyp,vegtyp,lai,roufns,veg,hterain,             &
              tsoil,qsoil,wetcanp, qvsfc,                        &
              usflx,vsflx,ptsflx,qvsflx,                         &
              windsp,psfc,rhoa,precip,                           &
              tair,qvair,                                        &
              cdha,cdqa,cdma,                                    &
              radsw, rnflx,                                      &
              shflx,lhflx,gflx,ct,                               &
              evaprg,evaprtr,evaprr,qvsat,                       &
              qvsata, f34,tem1soil)
  END IF

  RETURN
END SUBROUTINE soildiag
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTFLX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtflx(nx,ny,nzsoil,x,y,zpsoil,                       &
           soiltyp,vegtyp,lai,roufns,veg,hterain,                &
           tsoil,qsoil,wetcanp, qvsfc,                           &
           usflx,vsflx,ptsflx,qvsflx,                            &
           windsp,psfc,rhoa,precip,tair,qvair,                   &
           cdh,cdq,cdm,                                          &
           radsw, rnflx,                                         &
           shflx,lhflx,gflx, ct,                                 &
           evaprg,evaprtr,evaprr,qvsat,                          &
           qvsata,f34,tem1soil)
!
!-----------------------------------------------------------------
!
!  PURPOSE:
!
!  Write surface fields in GrADS format for diagnostic purpose.
!
!-----------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  4/15/1994.
!
!  MODIFICATION HISTORY:
!
!  10/30/94 (Y. Liu)
!  using the real names for variables instead of temporary array
!  names.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D array, veg(nx,ny), to the diagnostic printing list
!
!  05/31/2002 (J. Brotzge)
!  Added new soil variables.
!
!------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the soil
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!    hterain  The height of surface terrain
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy moisture
!    qvsfc    Effective specific humidity at sfc.
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    qvsflx   Surface flux of moisture (K*kg/(m**2*s))
!
!    windsp   Wind speed (m/s)
!    rhosfc   Surface air density (kg/m**3)
!    psfc     Surface pressure (Pascal)
!    preci    Precipitation flux reaching the surface
!    cdh      Surface drag coefficient for heat
!    cdq      Surface drag coefficient for moisture
!    cdm      Surface drag coefficient for momentum
!
!    radsw    Incoming solar radiation flux at surface
!    rnflx    Net radiation flux
!    shflx    Sensible heat flux
!    lhflx    Latent heat flux
!    gflx     Diffusive ground heat flux
!    evaprg   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evaprr   Direct evaporation from the fraction delta
!    f34     Surface resistence
!    ct       Thermal capacity
!    qvsat    Surface specific humidity at saturation, qvs(Ts)
!    qvsata   Surface air specific humidity at saturation, qvs(Ta)
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

  INTEGER :: nx,ny          ! The number grid points in 3 directions
  INTEGER :: nzsoil         ! The number grid points in the soil

  REAL :: x(nx)             ! X-coordinates
  REAL :: y(ny)             ! Y-coordinates
  REAL :: zpsoil (nx,ny,nzsoil)

  INTEGER :: soiltyp(nx,ny) ! Soil type at each point
  INTEGER :: vegtyp (nx,ny) ! Vegetation type at each point

  REAL :: lai    (nx,ny)    ! Leaf Area Index
  REAL :: roufns (nx,ny)    ! Surface roughness
  REAL :: veg    (nx,ny)    ! Vegetation fraction
  REAL :: hterain(nx,ny)    ! The height of surface terrain

  REAL :: qvsfc(nx,ny)        ! Effective S.H. at sfc.
  REAL :: tsoil(nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: qsoil(nx,ny,nzsoil) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny)      ! Canopy water amount

  REAL :: usflx  (nx,ny)    ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx  (nx,ny)    ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx (nx,ny)    ! surface flux of heat (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)    ! surface flux of moisture (kg/(m**2*s))

  REAL :: windsp (nx,ny)    ! Wind speed just above the surface (m/s)
  REAL :: psfc   (nx,ny)    ! Surface pressure (Pascal)
  REAL :: rhoa   (nx,ny)    ! Near sfc air density
  REAL :: precip (nx,ny)    ! Precipitation flux reaching the surface
  REAL :: tair   (nx,ny)    ! Air temperature near the surface
  REAL :: qvair  (nx,ny)    ! Specific humidity near the surface

  REAL :: cdh    (nx,ny)    ! Surface drag coefficient for heat
  REAL :: cdq    (nx,ny)    ! Surface drag coefficient for moisture
  REAL :: cdm    (nx,ny)    ! Surface drag coefficient for momentum

  REAL :: radsw  (nx,ny)    ! Incoming solar radiation at surface
  REAL :: rnflx  (nx,ny)    ! Net radiation flus
  REAL :: shflx  (nx,ny)    ! Sensible heat flux
  REAL :: lhflx  (nx,ny)    ! Latent heat flux
  REAL :: gflx   (nx,ny)    ! Diffusive heat flux from ground surface to
                            ! deep soil
  REAL :: ct     (nx,ny)    ! Thermal capacity

  REAL :: evaprg (nx,ny)    ! Evaporation from groud surface
  REAL :: evaprtr(nx,ny)    ! Transpiration of the remaining part
                            ! (1-delta) of leaves
  REAL :: evaprr (nx,ny)    ! Direct evaporation from the fraction delta
  REAL :: qvsat  (nx,ny)    ! qvs(ts)
  REAL :: qvsata (nx,ny)    ! qvs(ta)
  REAL :: f34    (nx,ny)
  REAL :: tem1soil (nx,ny,nzsoil) ! Temporary array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, m

  INTEGER :: tnum,tint
  INTEGER :: year1,month1,day1, hour1,minute1,second1
  INTEGER :: jday1, loopdy
  CHARACTER (LEN=2) :: dtunit

  INTEGER :: mndys(12)                 ! days for each months

  CHARACTER (LEN=3)   :: monnam(12)
  CHARACTER (LEN=256) :: flnctl, flnflx
  INTEGER :: flxunit, flnctlen, flxlen
  LOGICAL :: firstcall
  INTEGER :: ierr

  REAL :: latmin, latmax, lonmin, lonmax, latinc, loninc
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Save and initialize variables.
!
!-----------------------------------------------------------------------
!
  SAVE firstcall, flxunit,flnflx,flxlen
  DATA firstcall /.true./
  DATA mndys     /0,31,59,90,120,151,181,212,243,273,304,334/
  DATA monnam    /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',             &
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL xytoll(nx,ny,x,y,tem1soil(1,1,1),tem1soil(1,1,2))

  CALL a3dmax0(tem1soil(1,1,1),1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,             &
               latmax,latmin)
  CALL a3dmax0(tem1soil(1,1,2),1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,             &
               lonmax,lonmin)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

  IF ( firstcall ) THEN

    IF ( thisdmp <= 0.0 ) THEN
      WRITE (6, '(/a,a)')                                               &
          'Since thisdmp <= 0, only data at the first time step ',      &
          'will be dumped.'
      tnum = 1
      tint = 1
      dtunit = 'MN'
    ELSE IF ( thisdmp < 60.0 ) THEN
      WRITE (6, '(/a/a)')                                               &
          'GrADS reqiures the smallest uint minute for time interval.', &
          'Here we use uint MN to represent the second.'
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp)
      dtunit = 'MN'
    ELSE IF ( thisdmp < 3600.0 ) THEN
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp/60.)
      dtunit = 'MN'
    ELSE IF ( thisdmp < 86400.0 ) THEN
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp/3600.)
      dtunit = 'HR'
    ELSE
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp/86400.)
      dtunit = 'DY'
    END IF
    IF (tnum < 1) tnum = 1

    IF ( initopt /= 2 ) THEN
      second1 = second
      minute1 = minute
      hour1   = hour
      day1    = day
      month1  = month
      year1   = year
    ELSE
      second1 = MOD( second + nint(tstart), 60 )
      minute1 = ( second + nint(tstart) ) / 60
      minute1 = MOD( minute + minute1, 60 )
      hour1   = ( minute + ( second + nint(tstart) ) / 60 ) /60
      hour1   = MOD( hour + hour1, 24 )
      day1    = ( hour + ( minute                                       &
              + ( second + nint(tstart) ) / 60 ) /60 ) / 24
      jday1   = jday + day1

      loopdy  = 0
      IF ( MOD( year, 4 ) == 0 ) loopdy = 1
      year1 = year + jday1 / ( 365 + loopdy )
      jday1 = MOD( jday1, 365 + loopdy )

      month1 = 1

      DO m = 2, 11
        IF ( jday1 > mndys(m) .AND. jday1 <= mndys(m+1) + loopdy ) month1 = m
      END DO
      day1 = jday1 - mndys(month1)

    END IF

    flnctlen = lfnkey + 7
    flnctl(1:flnctlen) = runname(1:lfnkey)//'.sfcctl'
    CALL fnversn( flnctl, flnctlen )

    flnflx(1:ldirnam) = dirname(1:ldirnam)

    flxlen = ldirnam + lfnkey + 8
    flnflx(1:flxlen) = flnflx(1:ldirnam)//'/'//runname(1:lfnkey)        &
                     //'.sfcflx'
    IF (mp_opt > 0) THEN
      CALL gtsplitfn(runname(1:lfnkey)//'.flx',1,1,loc_x,loc_y,1,1,     &
                     0,0,0,lvldbg,flnflx,ierr)
      flxlen = LEN_TRIM(flnflx)
    END IF
    CALL fnversn( flnflx, flxlen )
!
!-----------------------------------------------------------------------
!
!  Open GrADS data control file for surface variables.
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) THEN

      CALL getunit (flxunit)
      OPEN (UNIT = flxunit, FILE = flnctl(1:flnctlen),                  &
            FORM = 'formatted', STATUS = 'new')

      WRITE (6,'(a,a,a)') 'The GrADS control file for surface ',        &
          'fluxes and other fields is ', flnctl(1:flnctlen)

      WRITE (flxunit,'(a,a)')                                           &
          'TITLE   Surface Fluxes, Temperature and Moisture for run ',  &
          runname(1:lfnkey)
      WRITE (flxunit,'(a)')                                             &
          '*'
      WRITE (flxunit,'(a,a)')                                           &
          'DSET    ', flnflx(1:flxlen)
      WRITE (flxunit,'(a)')                                             &
          'OPTIONS sequential big_endian'
!          'OPTIONS sequential cray_32bit_ieee'
      WRITE (flxunit,'(a)')                                             &
          'UNDEF   -9.e+33'
      WRITE (flxunit,'(a,i8,a,2f10.4)')                                 &
          'XDEF    ', nx, '  LINEAR   ', lonmin, loninc
      WRITE (flxunit,'(a,i8,a,2f10.4)')                                 &
          'YDEF    ', ny, '  LINEAR   ', latmin, latinc
      WRITE (flxunit,'(a,i8,a)')                                        &
          'ZDEF    ',nzsoil,'  LEVELS  '
      WRITE (flxunit,'(8f10.2)')                                        &
          ((zpsoil(1,1,k)+zpsoil(1,1,k-1))/2.,k=nzsoil,2,-1),           &
           zpsoil(1,1,1)/1.

      WRITE (flxunit,'(a,i8,a,i2.2,a,i2.2,a,i2.2,a3,i4.4,3X,i2.2,a)')   &
          'TDEF    ', tnum, '  LINEAR   ',                              &
          hour1,':',minute1,'Z',day1,monnam(month1),year1,tint,dtunit
      WRITE (flxunit,'(a)')                                             &
          '*'
      WRITE (flxunit,'(a)')                                             &
          'VARS   35'


      WRITE (flxunit,'(a)')                                             &
          'styp    0   -1,40,4   Soil type (4-byte integer)'
      WRITE (flxunit,'(a,a)')                                           &
          'vtyp    0   -1,40,4   Vegetation type 4-byte integer)'
      WRITE (flxunit,'(a)')                                             &
          'lai     0        99   Leaf Area Index'
      WRITE (flxunit,'(a)')                                             &
          'rfns    0        99   Surface roughness'
      WRITE (flxunit,'(a)')                                             &
          'veg     0        99   Vegetation fraction'
      WRITE (flxunit,'(a)')                                             &
          'trn     0        99   Surface terrain'
      WRITE (flxunit,'(a)')                                             &
          'va      0        99   Surface wind speed (m/s)'
      WRITE (flxunit,'(a)')                                             &
          'ps      0        99   Surface pressure (Pascal)'
      WRITE (flxunit,'(a)')                                             &
          'rhoa    0        99   Surface air density (kg/m**3)'
      WRITE (flxunit,'(a,a)')                                           &
          'rain    0        99   Surface precipitation rate ',          &
                                '(kg/s/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'ta      0        99   Surface air temperature (K)'
      WRITE (flxunit,'(a,a)')                                           &
          'qva     0        99   Surface specific humidity (k',         &
                                'g/kg) '
      WRITE (flxunit,'(a)')                                             &
          'ct      0        99   Surface Heat Capacity'
      WRITE (flxunit,'(a,a)')                                           &
          'qvsat   0        99   Specific humidity at ',                &
                                'ground surface'
      WRITE (flxunit,'(a)')                                             &
          'qvsata  0        99   Surface air specific humidity'
      WRITE (flxunit,'(a)')                                             &
          'f34     0        99   Surface resistence'
      WRITE (flxunit,'(a)')                                             &
          'cdh     0        99   Cdh'
      WRITE (flxunit,'(a)')                                             &
          'cdq     0        99   Cdq'
      WRITE (flxunit,'(a)')                                             &
          'cdm     0        99   Cdm'
      WRITE (flxunit,'(a)')                                             &
          'eg      0        99   Evaporation from ground'
      WRITE (flxunit,'(a,a)')                                           &
          'etr     0        99   Evaporation directly from ',           &
                                'the foliage'
      WRITE (flxunit,'(a,a)')                                           &
          'er      0        99   Transpiration of the part ',           &
                                'of the leaves'
      WRITE (flxunit,'(a,a)')                                           &
          'radsw   0        99   Incoming solar radiation ',            &
                                '(W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'rn      0        99   Net radiation (W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'h       0        99   Sensible heat flux (W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'le      0        99   Latent heat flux (W/m**2)'
      WRITE (flxunit,'(a,a)')                                           &
          'g       0        99   Ground diffusive heat flux'
      WRITE (flxunit,'(a,i2,a)')                                        &
          'tsoil  ',nzsoil,'        99   Soil temperature (K)'
      WRITE (flxunit,'(a,a)')                                           &
          'qvsfc   0        99   Surface water vapor mixing '
      WRITE (flxunit,'(a,i2,a)')                                        &
          'qsoil  ',nzsoil,'        99   Soil moisture (m**3/m**3)'
      WRITE (flxunit,'(a)')                                             &
          'wr      0        99   Canopy moisture'
      WRITE (flxunit,'(a)')                                             &
          'uflx    0        99   U flux'
      WRITE (flxunit,'(a)')                                             &
          'vflx    0        99   V flux'
      WRITE (flxunit,'(a)')                                             &
          'ptflx   0        99   PT flux'
      WRITE (flxunit,'(a)')                                             &
          'qvflx   0        99   QV flux'
      WRITE (flxunit,'(a)')                                             &
          'ENDVARS'

      CLOSE (flxunit)
      CALL retunit (flxunit)

    END IF

!-----------------------------------------------------------------------
!
!  Open GrADS data file for surface variables.
!
!-----------------------------------------------------------------------
!
    CALL getunit (flxunit)

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(flnflx(1:flxlen), '-F f77 -N ieee', ierr)

    OPEN (UNIT = flxunit, FILE = flnflx(1:flxlen),                      &
          FORM = 'unformatted', STATUS = 'new',ACCESS = 'sequential')

    firstcall = .false.

  END IF

  WRITE (flxunit) soiltyp                ! Soil type
  WRITE (flxunit) vegtyp                 ! Veg. type
  WRITE (flxunit) lai                    ! LAI
  WRITE (flxunit) roufns                 ! Roughness
  WRITE (flxunit) veg                    ! Veg
  WRITE (flxunit) hterain                ! Terrain

  CALL edgfill(windsp, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) windsp                 ! Va

  CALL edgfill(psfc,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) psfc                   ! Psfc

  CALL edgfill(rhoa,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) rhoa                   ! Sfc rhoa

  CALL edgfill(precip, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) precip                 ! Precipitation

  CALL edgfill(tair,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) tair                   ! Tair

  CALL edgfill(qvair,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) qvair                  ! Qvair

  CALL edgfill(ct,     1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) ct                     ! Ct

  CALL edgfill(qvsat,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) qvsat                  ! Qvsat

  CALL edgfill(qvsata, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) qvsata                 ! qvsata

  CALL edgfill(f34,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) f34                    ! f34

  CALL edgfill(cdh,    1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) cdh                    ! cdh

  CALL edgfill(cdq,    1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) cdq                    ! cdq

  CALL edgfill(cdm,    1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) cdm                    ! cdm

  CALL edgfill(evaprg, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) evaprg                 ! Eg

  CALL edgfill(evaprtr,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) evaprtr                ! Etr

  CALL edgfill(evaprr, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) evaprr                 ! Er

  CALL edgfill(radsw,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) radsw                  ! Radsw

  CALL edgfill(rnflx,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) rnflx                  ! Net rad. flux

  CALL edgfill(shflx,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) shflx                  ! H flux

  CALL edgfill(lhflx,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) lhflx                  ! LE flux

  CALL edgfill(gflx,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) gflx                   ! G flux

  CALL edgfill(tsoil,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
  DO k=nzsoil,1,-1
    WRITE (flxunit) ((tsoil(i,j,k),i=1,nx),j=1,ny)     ! Soil temp.
  END DO

  CALL edgfill(qvsfc,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) qvsfc                  ! Eff. SH dif.

  CALL edgfill(qsoil,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
  DO k=nzsoil,1,-1
    WRITE (flxunit) ((qsoil(i,j,k),i=1,nx),j=1,ny)     ! Soil moist
  END DO

  CALL edgfill(wetcanp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) wetcanp                ! Canopy moist

  CALL edgfill(usflx,  1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) usflx                  ! u flux

  CALL edgfill(vsflx,  1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
  WRITE (flxunit) vsflx                  ! v flux

  CALL edgfill(ptsflx, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) ptsflx                 ! pt flux

  CALL edgfill(qvsflx, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE (flxunit) qvsflx                 ! qv flux

  RETURN
END SUBROUTINE wrtflx
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTJOINFLX                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtjoinflx(nx,ny,nzsoil,x,y,zpsoil,                          &
           soiltyp,vegtyp,lai,roufns,veg,hterain,                       &
           tsoil,qsoil,wetcanp, qvsfc,                                  &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           windsp,psfc,rhoa,precip,tair,qvair,                          &
           cdh,cdq,cdm,                                                 &
           radsw, rnflx,                                                &
           shflx,lhflx,gflx, ct,                                        &
           evaprg,evaprtr,evaprr,qvsat,                                 &
           qvsata,f34,tem1soil)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write joined surface fields in GrADS format for parallel runs.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  10/18/2002.
!  Based on subroutine wrtflx.
!
!  MODIFICATION HISTORY:
!
!----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the soil
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!    hterain  The height of surface terrain
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy moisture
!    qvsfc    Effective specific humidity at sfc.
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    qvsflx   Surface flux of moisture (K*kg/(m**2*s))
!
!    windsp   Wind speed (m/s)
!    rhosfc   Surface air density (kg/m**3)
!    psfc     Surface pressure (Pascal)
!    preci    Precipitation flux reaching the surface
!    cdh      Surface drag coefficient for heat
!    cdq      Surface drag coefficient for moisture
!    cdm      Surface drag coefficient for momentum
!
!    radsw    Incoming solar radiation flux at surface
!    rnflx    Net radiation flux
!    shflx    Sensible heat flux
!    lhflx    Latent heat flux
!    gflx     Diffusive ground heat flux
!    evaprg   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evaprr   Direct evaporation from the fraction delta
!    f34     Surface resistence
!    ct       Thermal capacity
!    qvsat    Surface specific humidity at saturation, qvs(Ts)
!    qvsata   Surface air specific humidity at saturation, qvs(Ta)
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

  INTEGER :: nx,ny          ! The number grid points in 3 directions
  INTEGER :: nzsoil         ! The number grid points in the soil

  REAL :: x(nx)             ! X-coordinates
  REAL :: y(ny)             ! Y-coordinates
  REAL :: zpsoil (nx,ny,nzsoil)

  INTEGER :: soiltyp(nx,ny) ! Soil type at each point
  INTEGER :: vegtyp (nx,ny) ! Vegetation type at each point

  REAL :: lai    (nx,ny)    ! Leaf Area Index
  REAL :: roufns (nx,ny)    ! Surface roughness
  REAL :: veg    (nx,ny)    ! Vegetation fraction
  REAL :: hterain(nx,ny)    ! The height of surface terrain

  REAL :: qvsfc(nx,ny)        ! Effective S.H. at sfc.
  REAL :: tsoil(nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: qsoil(nx,ny,nzsoil) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny)      ! Canopy water amount

  REAL :: usflx  (nx,ny)    ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx  (nx,ny)    ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx (nx,ny)    ! surface flux of heat (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)    ! surface flux of moisture (kg/(m**2*s))

  REAL :: windsp (nx,ny)    ! Wind speed just above the surface (m/s)
  REAL :: psfc   (nx,ny)    ! Surface pressure (Pascal)
  REAL :: rhoa   (nx,ny)    ! Near sfc air density
  REAL :: precip (nx,ny)    ! Precipitation flux reaching the surface
  REAL :: tair   (nx,ny)    ! Air temperature near the surface
  REAL :: qvair  (nx,ny)    ! Specific humidity near the surface

  REAL :: cdh    (nx,ny)    ! Surface drag coefficient for heat
  REAL :: cdq    (nx,ny)    ! Surface drag coefficient for moisture
  REAL :: cdm    (nx,ny)    ! Surface drag coefficient for momentum

  REAL :: radsw  (nx,ny)    ! Incoming solar radiation at surface
  REAL :: rnflx  (nx,ny)    ! Net radiation flus
  REAL :: shflx  (nx,ny)    ! Sensible heat flux
  REAL :: lhflx  (nx,ny)    ! Latent heat flux
  REAL :: gflx   (nx,ny)    ! Diffusive heat flux from ground surface to
                            ! deep soil
  REAL :: ct     (nx,ny)    ! Thermal capacity

  REAL :: evaprg (nx,ny)    ! Evaporation from groud surface
  REAL :: evaprtr(nx,ny)    ! Transpiration of the remaining part
                            ! (1-delta) of leaves
  REAL :: evaprr (nx,ny)    ! Direct evaporation from the fraction delta
  REAL :: qvsat  (nx,ny)    ! qvs(ts)
  REAL :: qvsata (nx,ny)    ! qvs(ta)
  REAL :: f34    (nx,ny)
  REAL :: tem1soil (nx,ny,nzsoil) ! Temporary array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, m

  INTEGER :: tnum,tint
  INTEGER :: year1,month1,day1, hour1,minute1,second1
  INTEGER :: jday1, loopdy
  CHARACTER (LEN=2) :: dtunit

  INTEGER :: mndys(12)                 ! days for each months

  CHARACTER (LEN=3)   :: monnam(12)
  CHARACTER (LEN=256) :: flnctl, flnflx
  INTEGER :: flxunit, flnctlen, flxlen
  LOGICAL :: firstcall
  INTEGER :: ierr

  REAL :: latmin, latmax, lonmin, lonmax, latinc, loninc

  INTEGER :: nxlg, nylg
  REAL, ALLOCATABLE :: out2d(:,:), out3d(:,:,:)
  INTEGER, ALLOCATABLE :: out2di(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Save and initialize variables.
!
!-----------------------------------------------------------------------
!
  SAVE firstcall, flxunit,flnflx,flxlen
  DATA firstcall/.true./
  DATA mndys/0,31,59,90,120,151,181,212,243,273,304,334/
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

 nxlg = (nx-3)*nproc_x+3
 nylg = (ny-3)*nproc_y+3

 ALLOCATE(out2d(nxlg, nylg))
 ALLOCATE(out2di(nxlg, nylg))
 ALLOCATE(out3d(nxlg, nylg, nzsoil))

 CALL xytoll(nx,ny,x,y,tem1soil(1,1,1),tem1soil(1,1,2))

 CALL a3dmax0(tem1soil(1,1,1),1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,  &
              latmax,latmin)
 CALL a3dmax0(tem1soil(1,1,2),1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,  &
              lonmax,lonmin)

 IF(myproc == 0) THEN

   latinc = (latmax-latmin)/(nylg-1)
   loninc = (lonmax-lonmin)/(nxlg-1)

   IF ( firstcall ) THEN

    IF ( thisdmp <= 0.0 ) THEN
      WRITE (6, '(/a,a)')                                               &
          'Since thisdmp <= 0, only data at the first time step ',      &
          'will be dumped.'
      tnum = 1
      tint = 1
      dtunit = 'MN'
    ELSE IF ( thisdmp < 60.0 ) THEN
      WRITE (6, '(/a/a)')                                               &
          'GrADS reqiures the smallest uint minute for time interval.', &
          'Here we use uint MN to represent the second.'
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp)
      dtunit = 'MN'
    ELSE IF ( thisdmp < 3600.0 ) THEN
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp/60.)
      dtunit = 'MN'
    ELSE IF ( thisdmp < 86400.0 ) THEN
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp/3600.)
      dtunit = 'HR'
    ELSE
      tnum = nint(tstop/thisdmp)
      tint = nint(thisdmp/86400.)
      dtunit = 'DY'
    END IF
    IF (tnum < 1) tnum = 1

    IF ( initopt /= 2 ) THEN
      second1 = second
      minute1 = minute
      hour1   = hour
      day1    = day
      month1  = month
      year1   = year
    ELSE
      second1 = MOD( second + nint(tstart), 60 )
      minute1 = ( second + nint(tstart) ) / 60
      minute1 = MOD( minute + minute1, 60 )
      hour1   = ( minute + ( second + nint(tstart) ) / 60 ) /60
      hour1   = MOD( hour + hour1, 24 )
      day1    = ( hour + ( minute                                       &
            + ( second + nint(tstart) ) / 60 ) /60 ) / 24
      jday1   = jday + day1

      loopdy  = 0
      IF ( MOD( year, 4 ) == 0 ) loopdy = 1
      year1 = year + jday1 / ( 365 + loopdy )
      jday1 = MOD( jday1, 365 + loopdy )

      month1 = 1

      DO m = 2, 11
        IF ( jday1 > mndys(m) .AND. jday1 <= mndys(m+1) + loopdy ) month1 = m
      END DO
      day1 = jday1 - mndys(month1)

    END IF

    flnctlen = lfnkey + 7
    flnctl(1:flnctlen) = runname(1:lfnkey)//'.sfcctl'
    CALL fnversn( flnctl, flnctlen )

    flnflx(1:ldirnam) = dirname(1:ldirnam)

    flxlen = ldirnam + lfnkey + 8
    flnflx(1:flxlen) = flnflx(1:ldirnam)//'/'//runname(1:lfnkey)        &
                     //'.sfcflx'
    CALL fnversn( flnflx, flxlen )
!
!-----------------------------------------------------------------------
!
!  Open GrADS data control file for surface variables.
!
!-----------------------------------------------------------------------
!
      CALL getunit (flxunit)
      OPEN (UNIT = flxunit, FILE = flnctl(1:flnctlen),                  &
            FORM = 'formatted', STATUS = 'new')

      WRITE (6,'(a,a,a)') 'The GrADS control file for surface ',        &
          'fluxes and other fields is ', flnctl(1:flnctlen)

      WRITE (flxunit,'(a,a)')                                           &
          'TITLE   Surface Fluxes, Temperature and Moisture ',          &
          runname(1:lfnkey)
      WRITE (flxunit,'(a)')                                             &
          '*'
      WRITE (flxunit,'(a,a)')                                           &
          'DSET    ', flnflx(1:flxlen)
      WRITE (flxunit,'(a)')                                             &
          '*OPTIONS sequential cray_32bit_ieee'
      WRITE (flxunit,'(a)')                                             &
          'OPTIONS sequential big_endian'
      WRITE (flxunit,'(a)')                                             &
          'UNDEF   -9.e+33'
      WRITE (flxunit,'(a,i8,a,2f10.4)')                                 &
          'XDEF    ', nxlg, '  LINEAR   ', lonmin, loninc
      WRITE (flxunit,'(a,i8,a,2f10.4)')                                 &
          'YDEF    ', nylg, '  LINEAR   ', latmin, latinc
      WRITE (flxunit,'(a,i8,a)')                                        &
          'ZDEF    ',nzsoil,'  LEVELS  '
      WRITE (flxunit,'(8f10.2)')                                        &
          ((zpsoil(1,1,k)+zpsoil(1,1,k-1))/2.,k=nzsoil,2,-1),           &
           zpsoil(1,1,1)

      WRITE (flxunit,'(a,i8,a,i2.2,a,i2.2,a,i2.2,a3,i4.4,3X,i2.2,a)')   &
          'TDEF    ', tnum, '  LINEAR   ',                              &
          hour1,':',minute1,'Z',day1,monnam(month1),year1,              &
          tint,dtunit
      WRITE (flxunit,'(a)')                                             &
          '*'
      WRITE (flxunit,'(a)')                                             &
          'VARS   35'


      WRITE (flxunit,'(a)')                                             &
          'styp    0   -1,40,4   Soil type (4-byte integer)'
      WRITE (flxunit,'(a,a)')                                           &
          'vtyp    0   -1,40,4   Vegetation type 4-byte integer)'
      WRITE (flxunit,'(a)')                                             &
          'lai     0        99   Leaf Area Index'
      WRITE (flxunit,'(a)')                                             &
          'rfns    0        99   Surface roughness'
      WRITE (flxunit,'(a)')                                             &
          'veg     0        99   Vegetation fraction'
      WRITE (flxunit,'(a)')                                             &
          'trn     0        99   Surface terrain'
      WRITE (flxunit,'(a)')                                             &
          'va      0        99   Surface wind speed (m/s)'
      WRITE (flxunit,'(a)')                                             &
          'ps      0        99   Surface pressure (Pascal)'
      WRITE (flxunit,'(a)')                                             &
          'rhoa    0        99   Surface air density (kg/m**3)'
      WRITE (flxunit,'(a,a)')                                           &
          'rain    0        99   Surface precipitation rate ',          &
                                '(kg/s/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'ta      0        99   Surface air temperature (K)'
      WRITE (flxunit,'(a,a)')                                           &
          'qva     0        99   Surface specific humidity (k',         &
                                'g/kg) '
      WRITE (flxunit,'(a)')                                             &
          'ct      0        99   Surface Heat Capacity'
      WRITE (flxunit,'(a,a)')                                           &
          'qvsat   0        99   Specific humidity at ',                &
                                'ground surface'
      WRITE (flxunit,'(a)')                                             &
          'qvsata  0        99   Surface air specific humidity'
      WRITE (flxunit,'(a)')                                             &
          'f34     0        99   Surface resistence'
      WRITE (flxunit,'(a)')                                             &
          'cdh     0        99   Cdh'
      WRITE (flxunit,'(a)')                                             &
          'cdq     0        99   Cdq'
      WRITE (flxunit,'(a)')                                             &
          'cdm     0        99   Cdm'
      WRITE (flxunit,'(a)')                                             &
          'eg      0        99   Evaporation from ground'
      WRITE (flxunit,'(a,a)')                                           &
          'etr     0        99   Evaporation directly from ',           &
                                'the foliage'
      WRITE (flxunit,'(a,a)')                                           &
          'er      0        99   Transpiration of the part ',           &
                                'of the leaves'
      WRITE (flxunit,'(a,a)')                                           &
          'radsw   0        99   Incoming solar radiation ',            &
                                '(W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'rn      0        99   Net radiation (W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'h       0        99   Sensible heat flux (W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'le      0        99   Latent heat flux (W/m**2)'
      WRITE (flxunit,'(a)')                                             &
          'g       0        99   Ground diffusive heat flux'
      WRITE (flxunit,'(a,i2,a)')                                        &
          'tsoil  ',nzsoil,'        99   Soil temperature (K)'
      WRITE (flxunit,'(a,a)')                                           &
          'qvsfc   0        99   Surface water vapor mixing '
      WRITE (flxunit,'(a,i2,a)')                                        &
          'qsoil  ',nzsoil,'        99   Soil moisture (m**3/m**3)'
      WRITE (flxunit,'(a)')                                             &
          'wr      0        99   Canopy moisture'
      WRITE (flxunit,'(a)')                                             &
          'uflx    0        99   U flux'
      WRITE (flxunit,'(a)')                                             &
          'vflx    0        99   V flux'
      WRITE (flxunit,'(a)')                                             &
          'ptflx   0        99   PT flux'
      WRITE (flxunit,'(a)')                                             &
          'qvflx   0        99   QV flux'
      WRITE (flxunit,'(a)')                                             &
          'ENDVARS'

      CLOSE (flxunit)
      CALL retunit (flxunit)

!-----------------------------------------------------------------------
!
!  Open GrADS data file for surface variables.
!
!-----------------------------------------------------------------------
!
      CALL getunit (flxunit)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(flnflx(1:flxlen), '-F f77 -N ieee', ierr)

      OPEN (UNIT = flxunit, FILE = flnflx(1:flxlen),                    &
            FORM = 'unformatted', STATUS = 'new',ACCESS = 'sequential')

    END IF  ! firstcall
  END IF  ! myproc == 0

  firstcall = .false.

  CALL mpimerge2di(soiltyp,nx,ny,out2di)
  IF(myproc == 0) WRITE (flxunit) out2di                ! Soil type
  CALL mpimerge2di(vegtyp,nx,ny,out2di)
  IF(myproc == 0) WRITE (flxunit) out2di                ! Veg. type
  CALL mpimerge2d(lai,nx,ny,out2d)
  IF(myproc == 0) WRITE (flxunit) out2d                 ! LAI
  CALL mpimerge2d(roufns,nx,ny,out2d)
  IF(myproc == 0) WRITE (flxunit) out2d                 ! Roughness
  CALL mpimerge2d(veg,nx,ny,out2d)
  IF(myproc == 0) WRITE (flxunit) out2d                 ! Veg
  CALL mpimerge2d(hterain,nx,ny,out2d)
  IF(myproc == 0) WRITE (flxunit) out2d                 ! Terrain

  CALL mpimerge2d(windsp,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    WRITE (flxunit) out2d                 ! Va
  END IF

  CALL mpimerge2d(psfc,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    WRITE (flxunit) out2d                 ! Psfc
  END IF

  CALL mpimerge2d(rhoa,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! Sfc rhoa
  END IF

  CALL mpimerge2d(precip,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! Precipitation
  END IF

  CALL mpimerge2d(tair,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                   ! Tair
  END IF

  CALL mpimerge2d(qvair,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! Qvair
  END IF

  CALL mpimerge2d(ct,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,     1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                     ! Ct
  END IF

  CALL mpimerge2d(qvsat,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! Qvsat
  END IF

  CALL mpimerge2d(qvsata,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! qvsata
  END IF

  CALL mpimerge2d(f34,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,   1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                    ! f34
  END IF

  CALL mpimerge2d(cdh,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,    1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                    ! cdh
  END IF

  CALL mpimerge2d(cdq,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,    1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                    ! cdq
  END IF

  CALL mpimerge2d(cdm,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,    1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                    ! cdm
  END IF

  CALL mpimerge2d(evaprg,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! Eg
  END IF

  CALL mpimerge2d(evaprtr,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                ! Etr
  END IF

  CALL mpimerge2d(evaprr,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! Er
  END IF

  CALL mpimerge2d(radsw,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! Radsw
  END IF

  CALL mpimerge2d(rnflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! Net rad. flux
  END IF

  CALL mpimerge2d(shflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! H flux
  END IF

  CALL mpimerge2d(lhflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! LE flux
  END IF

  CALL mpimerge2d(gflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,   1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                   ! G flux
  END IF

  CALL mpimerge3d(tsoil,nx,ny,nzsoil,out3d)
  IF(myproc == 0) THEN
    CALL edgfill(out3d,1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,nzsoil,1,nzsoil)
    DO k=nzsoil,1, -1
      WRITE (flxunit) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)        ! Soil temp.
    END DO
  END IF

  CALL mpimerge2d(qvsfc,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! Eff. SH dif.
  END IF

  CALL mpimerge3d(qsoil,nx,ny,nzsoil,out3d)
  IF(myproc == 0) THEN
    CALL edgfill(out3d,1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,nzsoil,1,nzsoil)
    DO k=nzsoil,1,-1
      WRITE (flxunit) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)       ! Soil moist
    END DO
  END IF

  CALL mpimerge2d(wetcanp,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                ! Canopy moist
  END IF

  CALL mpimerge2d(usflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                  ! u flux
  END IF

  CALL mpimerge2d(vsflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d,  1,nxlg,1,nxlg-1, 1,nylg,1,nylg, 1,1,1,1)
    WRITE (flxunit) out2d                  ! v flux
  END IF

  CALL mpimerge2d(ptsflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! pt flux
  END IF

  CALL mpimerge2d(qvsflx,nx,ny,out2d)
  IF(myproc == 0) THEN
    CALL edgfill(out2d, 1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1)
    WRITE (flxunit) out2d                 ! qv flux
  END IF

  DEALLOCATE(out2d)
  DEALLOCATE(out2di)
  DEALLOCATE(out3d)
  RETURN
END SUBROUTINE wrtjoinflx
