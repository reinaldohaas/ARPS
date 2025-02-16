!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE RADIATION                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radiation(nx,ny,nz,rbufsz,                                   &
           ptprt,pprt,qv,qscalar,ptbar,pbar,ppi,rhostr,                 &
           x,y,z,zp,j3inv,                                              &
           soiltyp, tsfc, wetsfc,snowdpth,                              &
           radfrc, radsw, rnflx,radswnet,radlwin,                       &
           rsirbm,rsirdf,rsuvbm,rsuvdf,                                 &
           cosz, cosss,                                                 &
           fdirir,fdifir,fdirpar,fdifpar,                               &
           tem1, tem2, tem3, tem4, tem5,                                &
           tem6, tem7, tem8, tem9, tem10,                               &
           tem11,tem12,tem13,tem14,tem15,tem16,                         &
           radbuf, sh, tem17)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is the main driver of ARPS radiation package. To
!  control the driver, the parameters and variables should be set:
!
!  Input control variables (in namelist &radiation, arps.input)
!
!  radopt    Option to switch on/off radiation physics
!       = 0, No radiation physics;
!       = 1, Simplified surface radiation physics;
!       = 2, Atmospheric radiation transfer parameterization.
!
!  Notes:    When sfcphy is chosen to 3 or 4, radopt=0 will be adjusted
!         to 1 in order to compute the surface energy balance for
!         soil model. However radopt=2 will not be changed.
!
!  radstgr   Option for radiation computing at staggering points
!       = 0, No staggering; Radiation calculation on x-y plane is at
!            all points;
!       = 1, staggering; Radiation calculation on x-y plane is at
!            (even,even) and (odd,odd) points. The values at
!            (even,odd) (odd,even) points are averaged from the
!            surrounding four points. For example for nx=ny=9, the
!            directly calculation are performed at the "x" points,
!            then calculate radiation variables at "o" by averaging
!            from their surrounding "x" points. This scheme can
!            reduce ALMOST HALF of radiation calculation.
!
!
!                      j
!
!                    9 | x o x o x o x o x
!                    8 | o x o x o x o x o
!                    7 | x o x o x o x o x
!                    6 | o x o x o x o x o
!                    5 | x o x o x o x o x
!                    4 | o x o x o x o x o
!                    3 | x o x o x o x o x
!                    2 | o x o x o x o x o
!                    1 | x o x o x o x o x
!                      +-------------------  i
!                        1 2 3 4 5 6 7 8 9
!
!            On boundary, the zero-gradient is assumed.
!
!  rlwopt    Option to choose the longwave schemes.
!       = 0, high = .false. in code, transmission functions are
!            computed using the k-distribution method with linear
!            pressure scaling.  cooling rates are not calculated
!            accurately for pressures less than 20 mb. The
!            computation is faster with high=.false. than with
!            high=.true.
!       = 1, high = .true. in code, transmission functions in the
!            co2, o3 in the co2, o3, and the three water vapor bands
!            with strong absorption are computed using table look-up.
!            cooling rates are computed accurately from the surface
!            up to 0.01 mb.
!
!  dtrad     Time interval (seconds) to update the radiation forcing
!
!  raddiag   Option to dump radiation variables to a file in GrADS
!            format for diagnostic review. The frequency is
!            controled by dtrad. (Effective when radopt=2)
!       = 0, no such dump
!       = 1, dump to a file with a name like 'runname.radout'
!            and its control file has a name like 'runname.radctl'
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/11/1996
!
!  MODIFICATION HISTORY:
!
!  4/13/98 (D.Weber and V. Wong)
!  Added precipitable water term in the radopt=1 air emissivity
!  coefficient.
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  12/8/1998 (Donghai Wang and Vince Wong)
!  Added a new 2-D permanent array, snowcvr(nx,ny), for snow cover.
!  We just used a simple scheme to consider the snow cover process.
!
!  2000/01/10 (Gene Bassett)
!  Snow cover (0 or 1) changed to a fractional value (0 to 1)
!  determined from snow depth.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'radcst.inc'
  INCLUDE 'soilcst.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: rbufsz
!
!-----------------------------------------------------------------------
!
!  Define ARPS variables
!
!-----------------------------------------------------------------------
!
  REAL :: ptprt (nx,ny,nz)
  REAL :: pprt  (nx,ny,nz)
  REAL :: qv    (nx,ny,nz)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)
  REAL :: pbar  (nx,ny,nz)
  REAL :: ppi   (nx,ny,nz)
  REAL :: rhostr(nx,ny,nz)

  REAL :: x      (nx)
  REAL :: y      (ny)
  REAL :: z      (nz)

  REAL :: zp    (nx,ny,nz)  ! The physical height coordinate defined at
                            ! w-point of staggered grid.
  REAL :: j3inv (nx,ny,nz)

  INTEGER :: soiltyp(nx,ny) ! Soil type at each point
  REAL :: tsfc   (nx,ny)
  REAL :: wetsfc (nx,ny)    ! Surface soil moisture in the top 1 cm layer
  REAL :: snowdpth(nx,ny)   ! Snow depth (m)

  REAL :: radfrc(nx,ny,nz)  ! Radiation forcing (K/s)
  REAL :: radsw  (nx,ny)    ! Solar radiation down to the surface
  REAL :: rnflx  (nx,ny)    ! Net radiation flux absorbed by surface
  REAL :: radswnet (nx,ny)  ! Net solar radiation
  REAL :: radlwin  (nx,ny)  ! Incoming longwave radiation

  REAL :: rsirbm(nx,ny)     ! Solar IR surface albedo for beam radiation
  REAL :: rsirdf(nx,ny)     ! Solar IR surface albedo for diffuse radiation
  REAL :: rsuvbm(nx,ny)     ! Solar UV surface albedo for beam radiation
  REAL :: rsuvdf(nx,ny)     ! Solar UV surface albedo for diffuse radiation

  REAL :: cosz  (nx,ny)     ! Cosine of zenith
  REAL :: cosss (nx,ny)     ! Cosine of angle between sun light and
                            ! surface terrain slope

  REAL :: fdirir (nx,ny)    ! all-sky direct downward IR flux
                            ! (0.7-10 micron) at the surface
  REAL :: fdifir (nx,ny)    ! all-sky diffuse downward IR flux
                            ! at the surface
  REAL :: fdirpar(nx,ny)    ! all-sky direct downward par flux
                            ! (0.4-0.7 micron) at the surface
  REAL :: fdifpar(nx,ny)    ! all-sky diffuse downward par flux
                            ! at the surface
  REAL :: radbuf(rbufsz)    ! temporary arrays used for radiation
                            ! transfer computing
!
!-----------------------------------------------------------------------
!
!  Temporary arrays which have the vertical coordinate inversed, that
!  is, k=1 is for top while k=nz is at the surface.
!
!-----------------------------------------------------------------------
!
  REAL :: tem1 (nx,ny,nz)    ! pinv, slpmag, slpdir
  REAL :: tem2 (nx,ny,nz)    ! tinv
  REAL :: tem3 (nx,ny,nz)    ! qvinv
  REAL :: tem4 (nx,ny,nz)    ! o3a
  REAL :: tem5 (nx,ny,nz)    ! ccld
  REAL :: tem6 (nx,ny,nz)    ! flxir
  REAL :: tem7 (nx,ny,nz)    ! flcir
  REAL :: tem8 (nx,ny,nz)    ! flxuv
  REAL :: tem9 (nx,ny,nz)    ! flcuv
  REAL :: tem10(nx,ny,nz)    ! dfdts
  REAL :: tem11(nx,ny,nz)    ! tauir
  REAL :: tem12(nx,ny,nz)    ! taual
  REAL :: tem13(nx,ny,nz)    ! tauswi
  REAL :: tem14(nx,ny,nz)    ! tauswl
  REAL :: tem15(nx,ny,nz)    ! reffi
  REAL :: tem16(nx,ny,nz)    ! reffl

  ! add variable sh(nx,ny) & tem17
  REAL :: sh(nx,ny)          ! Work array for radiation shade
                             ! shadow induced by a topography
                             ! = 1 not shaded, =0 shaded
  REAL :: tem17(nx,ny,nz)    ! Work array for message passing
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k

  REAL :: albedo, albedoz
! augustin
! add variables saltitude,sazimuth
!     they are outputs of ZENANGL and inputs of SHADE 
! 
  REAL :: saltitude,sazimuth ! solar altitude and azimuth

  REAL :: tema,temb
  REAL :: rad2deg

  REAL :: frac_snowcover
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( radopt == 1 .AND. ( sfcphy /= 3 .AND. sfcphy /= 4 ) ) THEN
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate solar zenith angle. Array radsw is used as a2dr2
!  temporarily.
!
!    zp(*,*,2) = surface terrain
!    cosz    = cosine of zenith
!    cosss   = cosine of angle between sun light and terrain slope
!
!    radsw = a2dr2, square ratio of average distance to the time
!                   dependent distance from the earth to the sun
!
!  Array tem1 and tem2 are used as temporary 2-D arrays
!
!    tem1(*,*,1) = rjday
!    tem1(*,*,2) = tloc
!    tem1(*,*,3) = latscl
!    tem1(*,*,4) = lonscl
!
!    tem2(*,*,1) = slpmag
!    tem2(*,*,2) = slpdir
!    tem2(*,*,3) = tem1(*,*)
!    tem2(*,*,4) = tem2(*,*)
!
!-----------------------------------------------------------------------
!augustin
! add saltitude and sazimuth in the outputs of zenangl

  CALL zenangl( nx,ny, x,y, zp(1,1,2), cosz, cosss, radsw,              &
                tem1(1,1,1),tem1(1,1,2),tem1(1,1,3),tem1(1,1,4),        &
                tem2(1,1,1),tem2(1,1,2),tem2(1,1,3),tem2(1,1,4),        &
                saltitude,sazimuth )

!
!-----------------------------------------------------------------------
! augustin added shade
!
! Computes the shade induced by a topography
!
! radshade = 0  no shade computed
! radshade = 1  shade computed
! radshade = 2  shade computed on the whole topography
!               but only sh(i=1,nx,j=ny/2) is used
!               (useful for idealized norht south oriented valleys )
!
!-----------------------------------------------------------------------
!
  sh(:,:) = 1.0   ! default for radshade == 0
                  ! If radshade = 0 the matrix shade is set to 1 
                  ! at every grid point

  IF ( radshade /= 0 )                    &
     CALL shade(nx,ny,x,y,saltitude,sazimuth,zp(:,:,2),sh)

!
!-----------------------------------------------------------------------
!
!  Calculate surface albedo which is dependent on solar zenith angle
!  and soil moisture. Set the albedo for different types of solar
!  flux to be same.
!
!    rsirbm   Solar IR surface albedo for beam radiation
!    rsirdf   Solar IR surface albedo for diffuse radiation
!    rsuvbm   Solar UV surface albedo for beam radiation
!    rsuvdf   Solar UV surface albedo for diffuse radiation
!
!-----------------------------------------------------------------------
!
  rad2deg = 180.0/3.141592654

  DO j=1,ny-1
    DO i=1,nx-1

      albedoz = 0.01 * ( EXP( 0.003286         & ! zenith dependent albedo
          * SQRT( ( ACOS(cosz(i,j))*rad2deg ) ** 3 ) ) - 1.0 )

      IF ( sfcphy == 0 ) THEN             ! soil type not defined
        tema = 0
      ELSE
        tema = wetsfc(i,j)/wsat(soiltyp(i,j))
      END IF

      frac_snowcover = MIN(snowdpth(i,j)/snowdepth_crit, 1.0)

      IF ( tema > 0.5 ) THEN
        albedo = albedoz + (1.-frac_snowcover)*0.14                     &
                         + frac_snowcover*snow_albedo
      ELSE
        albedo = albedoz + (1.-frac_snowcover)*(0.31 - 0.34 * tema)     &
                         + frac_snowcover*snow_albedo
      END IF

      rsirbm(i,j) = albedo
      rsirdf(i,j) = albedo
      rsuvbm(i,j) = albedo
      rsuvdf(i,j) = albedo

    END DO
  END DO

  IF ( radopt == 1 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k) = (ptbar(i,j,k)+ptprt(i,j,k))*ppi(i,j,k)
          tem3(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Calculate solar radiation at the surface using a simplified
!  scheme.
!
!  tem1(*,*,1) = prcpln, precipitation path length
!
!  tem2        = temperature
!  tem3        = total pressure
!
!  fdirir      = fraction of solar radiation reaching the surface
!
!  Note: added precipitable water contribution to the emissa term
!        for radopt=1.
!
!-----------------------------------------------------------------------
!
    CALL solrsfc( nx,ny,nz,                                             &
                  tem2, tem3, qv, cosz, tem1(1,1,1), fdirir )

    DO j=1,ny-1
      DO i=1,nx-1
        tema = 0.5 * ( tem2(i,j,1) + tem2(i,j,2) ) ! surface air temp.

!
!  new code to include the water part of the downward emissa
!  coefficient....
!
        temb = 0.0

        DO k = 2,nz-2
          temb = temb+rhostr(i,j,k)*qv(i,j,k)*(zp(i,j,k+1)-zp(i,j,k))
        END DO

        IF (temb > 0.0) THEN
          temb = 0.17*ALOG10(temb*0.1)
        END IF
        temb = MIN(1.0,temb+emissa)

!  end of new code for emissa..

! augustin
! if the shading is not taken into account this just mutliplies radsw 
! and rnflx by 1
!
! else it sets radsw to zero in the shaded area
! and it sets to zero the direct compounds of the net radiation flux

        radsw(i,j) = solarc * radsw(i,j) * cosss(i,j) * fdirir(i,j)    &
                     * sh(i,j)

        rnflx(i,j) = radsw(i,j) * (1.0-rsirbm(i,j))                     &
                   + temb   * sbcst * tema**4                           &
                   - emissg * sbcst * tsfc(i,j)**4

        radswnet(i,j) = 0.9*radsw(i,j) 
        radlwin(i,j) = temb * sbcst * tema**4 

      END DO
    END DO

  ELSE IF ( radopt == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  Compute the atmospheric radiation forcing
!
!  Temporary arrays are used for
!
!    tem1  = pinv,  Pressure in mb at scalar points
!    tem2  = tinv,  Temperature
!    tem3  = qvinv, Water vapor mixing ratio (g/g)
!    tem4  = o3a,   Ozone (o3) mixing ratio (g/g)
!    tem5  = ccld,  Cloud coverage (fraction)
!    tem6  = flxir, all-sky net downward LW flux
!    tem7  = flcir, clear-sky net downward LW flux
!    tem8  = flxuv, all-sky solar flux (downward minus upward)
!    tem9  = flcuv, clear-sky solar flux (downward minus upward)
!    tem10 = dfdts, Sensitivity of net downward flux to surface temperature
!    tem11 = tauir,  Cloud optical depth for LW IR
!    tem12 = taual,  Aerosol optical thickness
!    tem13 = tauswi, Cloud optical depth for solar IR for ice particles
!    tem14 = tauswl, Cloud optical depth for solar IR for liquid particles
!    tem15 = reffi,  Effective cloud-particle size for ice particles
!    tem16 = reffl,  Effective cloud-particle size for liquid particles
!
!    cosss = cosine of angle between sun light and terrain slope
!
!    radsw = a2dr2, input,  square ratio of average distance to the
!                           time dependent distance from the earth to
!                           the sun
!          = radsw, output,    solar radiation reaching the surface
!    rnflx = st4,   temporary, longwave upward flux at surface
!          = rnflx, output,    net radiation flux absorbed by surface
!
!-----------------------------------------------------------------------
! augustin add sh 
!
    CALL radtrns(nx,ny,nz, rbufsz,                                      &
                 ptprt,pprt,qv,qscalar,                                 &
                 ptbar,pbar,ppi,rhostr, tsfc,                           &
                 x,y,z,zp, j3inv,                                       &
                 radfrc, radsw,rnflx,radswnet,radlwin, cosss,           &
                 rsirbm,rsirdf,rsuvbm,rsuvdf, cosz,                     &
                 fdirir,fdifir,fdirpar,fdifpar, rnflx,                  &
                 tem1, tem2, tem3, tem4, tem5,                          &
                 tem6, tem7, tem8, tem9, tem10,                         &
                 tem11,tem12,tem13,tem14,tem15,tem16,                   &
                 radbuf,sh, tem17)

  END IF

  RETURN

END SUBROUTINE radiation
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLRSFC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solrsfc(nx,ny,nz,                                            &
           temp, pres, qv, cosz, prcpln, trsw )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the atmosphere transmittance for solar radiation
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/25/1996
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
!
!-----------------------------------------------------------------------
!
!  Define ARPS variables
!
!-----------------------------------------------------------------------
!
  REAL :: temp  (nx,ny,nz)     ! temperature
  REAL :: pres  (nx,ny,nz)     ! total pressure
  REAL :: qv    (nx,ny,nz)     ! Mixing ratio

  REAL :: cosz  (nx,ny)

  REAL :: prcpln(nx,ny)        ! Precipitation path length

  REAL :: trsw  (nx,ny)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k

  REAL :: dirf
  REAL :: trrg
  REAL :: trwv
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
!  Calculate the precipitation path length by a vertical integral if
!  moist > 0.
!
!-----------------------------------------------------------------------
!
  DO j = 1, ny-1
    DO i = 1, nx-1
      prcpln(i,j) = 0.0
    END DO
  END DO

  IF ( moist /= 0 ) THEN
    DO j = 1, ny-1
      DO i = 1, nx-1

        prcpln(i,j) = 0.0

        DO k = 2, nz-2
          prcpln(i,j) = prcpln(i,j)                                     &
                      + 0.5*(qv(i,j,k)+qv(i,j,k+1))                     &
                      * 0.5*(pres(i,j,k)+pres(i,j,k+1))/101300.0        &
                      * SQRT(2.*273.16/(temp(i,j,k)+temp(i,j,k+1)))     &
                      * (pres(i,j,k)-pres(i,j,k+1))
        END DO

        prcpln(i,j) = prcpln(i,j) / g
      END DO
    END DO
  END IF

  DO j = 1, ny-1
    DO i = 1, nx-1
!
!-----------------------------------------------------------------------
!
!  Calculate the direction fractor of transmisivity
!
!-----------------------------------------------------------------------
!
      dirf = 35.0 / SQRT( 1224.0 * cosz(i,j)**2 + 1.0 )
!
!-----------------------------------------------------------------------
!
!  Calculate Rayleigh sccattering and absorption transmisivity, trrg.
!
!-----------------------------------------------------------------------
!
      trrg = 1.021 - 0.084 * SQRT( dirf                                 &
                           * ( 949. * 0.5*(pres(i,j,2)+pres(i,j,1))     &
                           * 1.e-8 + .051 ) )
                                ! Rayleigh scatt. and abs. transmission
! AB, Eq. 2, psfc in Pascal
!
!-----------------------------------------------------------------------
!
!  Calculate the transmisivity of water vapor, trwv. The
!  precipitation path length used in the calculation has been
!  calculated and stored in prcpln(i,j).
!
!  For cloud cover, use constant dirfc = 5/3 as the direction
!  fractor. (How to determine if cloudy or not?)
!
!  For clear sky, use dirf as the direction fractor, instead of dirfc
!
!-----------------------------------------------------------------------
!
      trwv = 1.0 - 2.9 * prcpln(i,j) * dirf                             &
           / ( 5.925 * prcpln(i,j) * dirf                               &
           + ( 1. + 141.5 * prcpln(i,j) * dirf ) ** 0.634 )
!
!-----------------------------------------------------------------------
!
!  Calculate the fraction of solar radiation reaching the surface.
!
!-----------------------------------------------------------------------
!
      trsw(i,j) = trrg * trwv

    END DO
  END DO

  RETURN
END SUBROUTINE solrsfc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ZENANGL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE zenangl(nx,ny, x,y, hterain, cosz, cosss, a2dr2,             &
           rjday,tloc, latscl,lonscl, slpmag,slpdir,                    &
           tem1,tem2,saltitude,sazimuth)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate cosine of solar zenith angle.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  11/16/93
!
!  MODIFICATION HISTORY:
!  2/2/99  Vince Wong and Jik Leong
!  This modification calculates the solar declination angle and
!  equation of time using a method found on page C24 of the
!  1996 Astronomical Almanac.
!  The mothod is good to 0.01 degrees in the sky over the
!  period 1950 to 2050.
!
! augustin
!  8/23/01 Augustin Colette EFML/ Stanford University
!
! Computation of the solar altitude and azimuth
! saltitude and sazimuth are outputs of zenangl to be used in shade
! sources:
! http://www.usc.edu/dept/architecture/mbs/tools/vrsolar/Help/ &
! solar_concepts.html
! http://www.uwinnipeg.ca/~blair/physclim/lab2.htm
! http://ra.stsci.edu/cgi-bin/gethelp.cgi?altaz.src
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    x        X coordinates at scalar points
!    y        Y coordinates at scalar points
!    hterain  Surface terrain
!
!  OUTPUT:
!
!    cosz     Cosine of zenith
!    cosss    Cosine of angle between sun light and terrain slope
!    a2dr2    Square ratio of average distance to the time
!             dependent distance from the earth to the sun
!
!augustin
!    sazimuth  solar azimuth
!    saltitude solar altitude
!
!  WORK ARRAY:
!
!    rjday    Julian day at each grid point
!    tloc     Local time at each grid point
!    latscl   Latitudes  at scalar points
!    lonscl   Longitudes at scalar points
!    slpmag   Surface terrain slope magnitude
!    slpdir   Surface terrain slope direction
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny

  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: hterain(nx,ny)

  REAL :: cosz(nx,ny),   cosss(nx,ny), a2dr2(nx,ny)
  REAL :: rjday(nx,ny),  tloc(nx,ny)
  REAL :: latscl(nx,ny), lonscl(nx,ny)
  REAL :: slpmag(nx,ny), slpdir(nx,ny)
!augustin add saltitude and sazimuth in the outputs of zenanlg
  REAL :: saltitude
  REAL :: sazimuth

  REAL :: tem1(nx,ny), tem2(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j

  REAL :: xs, ys

  REAL :: hour0, yrday
  REAL :: deg2rad, pi, pi2

  REAL :: etau, shrangl, sdeclin
  REAL :: azimuth, sinz
  REAL :: dpsi, sinpsi, cospsi

  REAL :: anncyc

  REAL :: hr, days2k, lsun, gsun, obliq, lambda, xsun, ysun
  REAL :: asun, alpha, rad2deg

  LOGICAL :: firstcall        ! First call flag of this subroutine

  SAVE firstcall, hour0, pi, pi2, deg2rad, yrday, rad2deg
  DATA firstcall/.true./

  ! added to calculate saltitude and sazimuth at the middle domain
  INTEGER :: nxmid, nymid, source 
  SAVE nxmid, nymid, source

  REAL :: shrangl_mid
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (firstcall) THEN
    pi  = 3.14159265358979
    pi2 = 2.0 * pi
    deg2rad = pi/180.0
    rad2deg = 1./deg2rad

    hour0 = FLOAT(hour)                                                 &
          + FLOAT(minute)/60.0                                          &
          + FLOAT(second)/3600.0

    IF ( MOD(year, 4) == 0 ) THEN
      yrday = 366.
    ELSE
      yrday = 365.
    END IF

    nxmid = CEILING( 0.5*((nx-3)*nproc_x + 3) )   ! Middle point index at
    nymid = CEILING( 0.5*((ny-3)*nproc_y + 3) )   ! global domain
    source = proc( (nxmid-2)/(nx-3)+1 + ( (nymid-2)/(ny-3) )*nproc_x )
                                ! source processor contain the middle point.
    nxmid = MOD( (nxmid-2), (nx-3) ) + 2    ! local index of central domain
    nymid = MOD( (nymid-2), (ny-3) ) + 2

    firstcall = .false.
  END IF

  CALL sfcslp( nx,ny, hterain, slpmag,slpdir, tem1,tem2 )

  IF ( mapproj == 0 ) THEN
    DO j=1,ny-1
      DO i=1,nx-1
        latscl(i,j) = ctrlat
        lonscl(i,j) = ctrlon
      END DO
    END DO
  ELSE
    DO j=1,ny-1
      ys = 0.5*(y(j)+y(j+1))
      DO i=1,nx-1
        xs = 0.5*(x(i)+x(i+1))
        CALL xytoll(1,1, xs,ys, latscl(i,j), lonscl(i,j))
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the local time at each grid point. The
!  following formula is based on that the input time is the GMT
!  time at the reference grid point of the center.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1

      latscl(i,j) = deg2rad * latscl(i,j)        ! lat: -90 to 90

      tloc(i,j) = hour0 + (curtim-dtbig)/3600.0                         &
                + lonscl(i,j)/15.0

      rjday(i,j) = jday + INT( tloc(i,j)/24.0 )

      tloc(i,j) = MOD( tloc(i,j), 24.0 )

      IF ( tloc(i,j) < 0. ) THEN
        tloc(i,j) = tloc(i,j) + 24.0            ! Local time
        rjday(i,j) = MOD( rjday(i,j)-1, yrday ) ! Julian day at each pts
      END IF

    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      anncyc = pi2 * ( rjday(i,j) - 1.0 ) / yrday

      a2dr2(i,j) = 1.000110                                             &
                 + 0.034221 * COS(anncyc)                               &
                 + 0.00128  * SIN(anncyc)                               &
                 + 0.000719 * COS(2.*anncyc)                            &
                 + 0.000077 * SIN(2.*anncyc)  ! PX, Eq. 17

      hr = hour + minute / 60.

!  days before (-ve) or after (+ve) 1/1/2000
      days2k = 367 * year - 7 * ( year + ( month + 9 ) / 12 ) / 4       &
              + 275 * month / 9 + day - 730531.5 + hr / 24.

      lsun = 280.461 + 0.9856474 * days2k     ! Mean Longitude of the Sun
      950     IF ( lsun < 0 ) THEN
        lsun = lsun + 360.
        GO TO 950
      ELSE IF ( lsun > 360 ) THEN
        lsun = lsun - 360.
        GO TO 950
      END IF

      gsun = 357.528 + 0.9856003 * days2k     ! Mean anomaly of the Sun
      960     IF ( gsun < 0 ) THEN
        gsun = gsun + 360.
        GO TO 960
      ELSE IF ( gsun > 360 ) THEN
        gsun = gsun - 360.
        GO TO 960
      END IF

      lambda = lsun + 1.915 * SIN(gsun*deg2rad)  & ! Ecliptic longitude
          + 0.02 * SIN(2*gsun*deg2rad)
      970     IF ( lambda < 0 ) THEN
        lambda = lambda + 360.
        GO TO 970
      ELSE IF ( lambda > 360 ) THEN
        lambda = lambda - 360.
        GO TO 970
      END IF

      obliq = 23.439 - 0.0000004 * days2k     ! Obliquity of the ecliptic

      xsun = COS(lambda*deg2rad)
      ysun = COS(obliq*deg2rad) * SIN(lambda*deg2rad)
      asun = ATAN(ysun/xsun)*rad2deg
      IF ( xsun < 0. ) THEN
        alpha = asun + 180   ! Right Ascension (RA)
      ELSE IF ( ( ysun < 0. ) .AND. ( xsun > 0. ) ) THEN
        alpha = asun + 360
      ELSE
        alpha = asun
      END IF

      etau = ( lsun - alpha ) * 4. / 60.      ! Equation of time in hour

!    etau = 0.158 * sin( pi*(rjday(i,j)+10.)/91.25 ) ! Equation of time
!    :       + 0.125 * sin( pi*rjday(i,j)/182.5 )       ! Wong, Eq. 8

      shrangl = 15.0 * deg2rad                        & ! Hour angle
                     * ( tloc(i,j) + etau - 12.0)       ! Wong, Eq. 7

!    sdeclin = 23.5 * deg2rad
!    :          * cos( 2.0*pi*(rjday(i,j)-173.)/yrday ) ! Wong, Eq. 6
      sdeclin = ASIN(SIN(obliq*deg2rad)*SIN(lambda*deg2rad))
                                                ! Declination (in radian)

      cosz(i,j) = COS(latscl(i,j)) * COS(sdeclin) * COS(shrangl)        &
                + SIN(latscl(i,j)) * SIN(sdeclin)

!    print *, cos(latscl(i,j)),cos(sdeclin),cos(shrangl)
!    print *, sin(latscl(i,j)),sin(sdeclin)
!    print *,sdeclin,shrangl

      sinz = SIN ( ACOS(cosz(i,j)) )
!
!-----------------------------------------------------------------------
!
!  Consider the effects of the terrain slope on the solar radiation.
!  The slope magnitude and direction has been computed by subroutine
!  SFCSLP and passed in by slpmag and slpdir.
!
!-----------------------------------------------------------------------
!

      sinpsi = COS(sdeclin) * SIN(shrangl) *COS(latscl(i,j))
      cospsi = COSZ(i,j) * SIN( latscl(i,j) ) - SIN( sdeclin )
      azimuth = ATAN2( sinpsi, cospsi)

      dpsi = azimuth - slpdir(i,j)

      cosss(i,j) = COS( slpmag(i,j) ) * cosz(i,j)                       &
                 + SIN( slpmag(i,j) ) * sinz * COS( dpsi )

      cosz (i,j) = MAX( cosz (i,j), 0.0 )
      cosss(i,j) = MAX( cosss(i,j), 0.0 )

! added to calculate saltitude and sazimuth later for shading
      IF( radshade /= 0) THEN

        IF (i == nxmid .AND. j == nymid) THEN
  
! augustin
! computes the solar altitudes and azimuth,
! sazimuth: 0=North; pi/2=East; Pi=South; 3*pi/2=West
! see the reference listed in the header of subroutine ZENANGL for more 
! information

        saltitude = pi/2-ACOS(cosz(i,j))

!  the following code is not computation safe (division by zero) DBW
!       sazimuth  = ACOS ( (SIN(sdeclin) * COS (latscl(i,j))             &
!                   - COS(sdeclin) * SIN (latscl(i,j))                   &
!                      * COS (shrangl))/sinz)

!  the above is replaced with:   DBW

        sinpsi = COS(sdeclin) * SIN(shrangl) * COS(latscl(i,j))
        cospsi = - COSZ(i,j) * SIN( latscl(i,j) ) + SIN( sdeclin )
        sazimuth = ABS( ATAN2(sinpsi, cospsi) )

!  end or replacement code  DBW

        shrangl_mid = shrangl

       END IF

      END IF  !  end of radshade if block...

    END DO
  END DO

  IF( radshade /= 0) THEN

! broadcast these value from source to all other processors

    CALL mpbcastr(saltitude,   source)
    CALL mpbcastr(sazimuth,    source)
    CALL mpbcastr(shrangl_mid, source)

! modification of shrangl

! In the definition of the hour angle, shrangl  should be
! between -pi and 0 before solar noon and between 0 and pi after.
! Before this error, zenangl gave positive values
! between pi and 2*pi of shrangl
! before solar noon for some latitudes and day.
! for instance, +45 deg north March 21 2001
! it causes not problem in the original ARPS since it is just 
! a problem of module but subroutine shade needs to have shrangl 
! between -pi and pi


    IF ((shrangl_mid > pi) .AND. (shrangl_mid < pi2)) THEN
      shrangl_mid = shrangl_mid - pi2
    END IF
    IF ((shrangl_mid > pi2) .AND. (shrangl_mid < pi2+pi)) THEN
      shrangl_mid = shrangl_mid - pi2
    END IF

    IF (shrangl_mid > 0) THEN
      sazimuth = pi2 - sazimuth
    END IF

  END IF  !  end of radshade if block...

  RETURN
END SUBROUTINE zenangl
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SFCSLP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sfcslp( nx,ny, hterain, slpmag,slpdir, tem1,tem2 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This program is to compute the terrain slope vector in terms of
!  magnitude and direction.
!
!  Magnitude slpmag:
!
!                                   1
!    cos(slpmag) = -----------------------------------,   0 <= slpmag <= pi/2
!                  sqrt( 1 + (dz/dx)**2 + (dz/dy)**2 )
!
!  Direction slpdir:
!
!                  dz/dx
!    tan(slpdir) = -----,                              - pi <= slpdir <= pi
!                  dz/dy
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  2/16/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction (east/west)
!  ny       Number of grid points in the y-direction (north/south)
!
!  hterain  The heights of terrain
!
!  OUTPUT:
!
!  slpmag   Magnitude of terrain surface slope vector
!  slpdir   Direction of terrain surface slope vector from north
!           clockwise
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny

  REAL :: hterain(nx,ny)    ! The heights of terrain

  REAL :: slpmag (nx,ny)    ! Terrain surface slope vector magnitude
  REAL :: slpdir (nx,ny)    ! Terrain surface slope vector direction
                            ! from north clockwise

  REAL :: tem1(nx,ny)       ! 2-D temporary array
  REAL :: tem2(nx,ny)       ! 2-D temporary array
!
!-----------------------------------------------------------------------
!
!  Local declarations
!
!-----------------------------------------------------------------------
!
  REAL :: pi
  PARAMETER ( pi = 3.141592654)

  REAL :: twdxinv, twdyinv, dzsds2

  INTEGER :: i, j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( ternopt == 0 ) THEN

    DO j = 1, ny-1
      DO i = 1, nx-1
        slpmag(i,j) = 0.
        slpdir(i,j) = 0.
      END DO
    END DO

  ELSE
!
!-----------------------------------------------------------------------
!
!                                           dz       dz
!  Use the central difference to calculate ---- and ----.
!                                           dx       dy
!
!-----------------------------------------------------------------------
!
    twdxinv = 1.0 / (2.0*dx)

    DO j = 1, ny-1
      DO i = 2, nx-2

        tem1(i,j) = ( hterain(i+1,j) - hterain(i-1,j) ) * twdxinv

      END DO
    END DO

    twdyinv = 1.0 / (2.0*dy)

    DO j = 2, ny-2
      DO i = 1, nx-1

        tem2(i,j) = ( hterain(i,j+1) - hterain(i,j-1) ) * twdyinv

      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv1dew(tem1,nx,ny,ebc,wbc,0,slpmag)  ! use slpmag as temp
      CALL mpsendrecv1dns(tem2,nx,ny,nbc,sbc,0,slpmag)  ! use slpmag as temp
    END IF

    CALL acct_interrupt(bc_acct)

    IF ( wbc == 1 ) THEN

      DO j = 1,ny-1
        tem1(1,j) = - tem1(2,j)
      END DO

    ELSE IF ( wbc == 2 ) THEN

      IF (mp_opt == 0) THEN
        DO j = 1,ny-1
          tem1(1,j) = tem1(nx-2,j)
        END DO
      END IF

    ELSE IF ( wbc /= 0 ) THEN

      DO j = 1, ny-1
        tem1(1,j) = ( hterain(2,j) - hterain(1,j) ) / dx
      END DO

    END IF

    IF ( ebc == 1 ) THEN

      DO j = 1, ny-1
        tem1(nx-1,j) = - tem1(nx-2,j)
      END DO

    ELSE IF ( ebc == 2 ) THEN

      IF (mp_opt == 0) THEN
        DO j = 1,ny-1
          tem1(nx-1,j) = tem1(2,j)
        END DO
      END IF

    ELSE IF ( ebc /= 0 ) THEN

      DO j = 1, ny-1

        tem1(nx-1,j) = ( hterain(nx-1,j) - hterain(nx-2,j) ) / dx

      END DO

    END IF

    IF ( sbc == 1 ) THEN

      DO i = 1, nx-1
        tem2(i,1) = - tem2(i,2)
      END DO

    ELSE IF ( sbc == 2 ) THEN

      IF (mp_opt == 0) THEN
        DO i = 1, nx-1
          tem2(i,1) = tem2(i,ny-2)
        END DO
      END IF

    ELSE IF ( sbc /= 0 ) THEN

      DO i = 1, nx-1
        tem2(i,1) = ( hterain(i,2) - hterain(i,1) ) / dy
      END DO

    END IF

    IF ( nbc == 1 ) THEN

      DO i = 1, nx-1
        tem2(i,ny-1) = - tem2(i,ny-2)
      END DO

    ELSE IF ( nbc == 2 ) THEN

      IF (mp_opt == 0) THEN
        DO i = 1, nx-1
          tem2(i,ny-1) = tem2(i,2)
        END DO
      END IF

    ELSE IF ( nbc /= 0 ) THEN

      DO i = 1, nx-1
        tem2(i,ny-1) = ( hterain(i,ny-1) - hterain(i,ny-2) ) / dy
      END DO

    END IF

    CALL acct_stop_inter

    DO j = 1, ny-1
      DO i = 1, nx-1

        dzsds2 = tem1(i,j)**2 + tem2(i,j)**2

        slpmag(i,j) = ACOS( 1. / SQRT( 1. + dzsds2 ) )

        IF ( dzsds2 == 0. ) THEN
          slpdir(i,j) = 0.
        ELSE
          slpdir(i,j) = ATAN2( tem1(i,j), tem2(i,j) )
        END IF

      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE sfcslp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SHADE                      ######
!######                                                      ######
!######                 Developed by                         ######
!######        Augustin Colette / Robert L. Street           ######
!######             Stanford University                      ######
!######      Environmental Fluid Mechanics Laboratory        ######
!######                   Stanford, CA 94305                 ######
!######      augustin@stanford.edu, street@stanford.edu      ######
!######                                                      ######
!######                                                      ######
!######                   June 2002                          ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE shade(nx,ny,x,y,saltitude,sazimuth,hterain,sh)

!
!-----------------------------------------------------------------------
!
!  PURPOSE: 
!
!  This subroutine computes the shadow induced by a topography
!  by drawing a line from each node and the sun
!  and checking if it hits the topography
!
!----------------------------------------------------------------------
! Modification:
!
!  2003/02/19 (Yunheng Wang)
!  Added code to do message passing for multiple processors.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx                Number of grid points in the x-direction (east/west)
!    ny                Number of grid points in the y-direction (north/south)
!    x(nx)             x-coordinates
!    y(ny)             y-coordinates
!    sazimuth          azimuth (radians)
!    saltitude         solar altitude (radians)
!    hterain(nx,ny)    topography
!
!  OUTPUT:
!    sh(nx,ny)         Matrix shade: 0: shaded, 1:enlightened.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny
!
!-----------------------------------------------------------------------
!
!  Define ARPS variables
!
!-----------------------------------------------------------------------
!
  REAL :: x(nx)            ! x-coordinates
  REAL :: y(ny)            ! y-coordinates
  REAL :: sazimuth         ! azimuth  (radians)
  REAL :: saltitude        ! solar altitude  (radians)
  REAL :: hterain(nx,ny)   ! topography

!-----------------------------------------------------------------------
!  
!  Output
! 
!-----------------------------------------------------------------------
  
  REAL :: sh(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  REAL :: ztest               ! altitude of a test point on a sun ray
  REAL :: htest               ! averaged topography at the test point
  INTEGER :: i,j,l            ! indices
  INTEGER :: II,IF,dI,JI,JF,dJ   ! extrema of the iteration
  REAL :: a,b,C,D,E,F,H       ! definition of the quadrant (NE,NW,SW,SE)
  REAL :: xg,yg               ! local grid cell where the test
                              ! is performed
  REAL :: Sx,Sy               ! coordinates of the point to be tested
  REAL :: pi,pi2,pi32,p2i     ! multiple of pi
  INTEGER :: xh,yh,xhh,yhh

!  CHARACTER (LEN = 132) :: fname

!----------------------------------------------------------------------
!
! MP variables
!
!----------------------------------------------------------------------
  INTEGER :: nxlg, nylg, istat
  REAL, ALLOCATABLE :: xlg(:), ylg(:)
  REAL, ALLOCATABLE :: hterainlg(:,:)
  REAL, ALLOCATABLE :: shlg(:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  ALLOCATE(xlg(nxlg), STAT=istat)
  ALLOCATE(ylg(nylg), STAT=istat)
  ALLOCATE(hterainlg(nxlg,nylg), STAT=istat)
  ALLOCATE(shlg(nxlg,nylg), STAT=istat)

  CALL mpimerge1dx(x,nx,xlg)
  CALL mpimerge1dy(y,ny,ylg)
  CALL mpimerge2d(hterain, nx,ny, hterainlg)

! initialization
! in the subroutine, sh will be set to zero in the shaded areas
! nothing will be done if the point is illuminated

  shlg(:,:)=1
  
  pi=3.14159265358979;
  pi2=pi/2.0; 
  pi32=3.0 * pi/2.0;
  p2i=2.0 * pi;
  
! modification of azimuth
! in the subroutine zenangl azimuth has the usual definition:
! sazimuth: 0:N pi/2:E, Pi:S, 3*pi/2:W
! 
! in the subroutine shade azimuth is 0:E, pi/2:N, pi:W, 3pi/2:S
! (for better comparison with the trigonometric circle)

! Check if azimuth is > 0 and < 2*pi
  IF (sazimuth > p2i) THEN
!  print*,'aziumth zenangl was > 2pi'
    sazimuth =sazimuth-p2i
  END IF
  IF (sazimuth < 0) THEN
!  print*,'aziumth zenangl was <0'
    sazimuth =p2i+sazimuth
  END IF

! second translate sazimuth 
  IF (( sazimuth >= 0 ) .AND. ( sazimuth < pi/2)) THEN !
     sazimuth=pi/2.0-sazimuth
  ELSE
     sazimuth=p2i+pi/2.0-sazimuth
  END IF

!print*,'azshade',sazimuth
!print*,'sdshade',saltitude

!-----------------------------------------------------------
!
! Four differents cases must be considered
! the iteration will vary if the sun is NE,NW,SW,SE because:
! - the shade can not be computed on the most nothern point
!   if the sun is north
! - when going along the line going towards the sun, the direction
!   of 'propagation' will vary
! 
! consequently the iteration extremas and increment vary
! to avoid writting four times the subroutine, local variables
! are used depending on the azimuth of the sun
!
!---------------------------------------------------------

! First quadrant : NE 
  IF  ( (0 <= sazimuth) .AND. (sazimuth < pi2) ) THEN
        a=1
        b=1
        C=1
        D=1
        E=0
        F=0
        H=0
  END IF 
! Second quadrant NW
  IF  ( (pi2 <= sazimuth) .AND. (sazimuth < pi) ) THEN
     a=-1
     b=1
     C=2
     D=0
     E=1
     F=0
     H=0
  END IF  
! Third quadrant SW
  IF  ((pi <= sazimuth) .AND. (sazimuth < pi32) ) THEN
     a=-1
     b=-1
     C=3
     D=0
     E=0
     F=1
     H=0
  END IF
! Fourth quadrant SE
  IF  ((pi32 <= sazimuth) .AND. (sazimuth <p2i)) THEN
     a=1
     b=-1
     C=4
     D=0
     E=0
     F=0
     H=1
  END IF
  
! Iteration extremas and increment

  II=D*2+E*(nxlg-2)+F*(nxlg-2)+H*2        ! initial i index
  IF=D*(nxlg-3)+E*(3)+F*(3)+H*(nxlg-3)    ! final i index
  dI=a*1                              ! i increment
  JI=D*2+E*2+F*(nylg-2)+H*(nylg-2)        ! initial j index
  JF=D*(nylg-3)+E*(nylg-3)+F*3+H*3        ! final j index
  dJ=b*1                              ! j increment
 
!----------------------------------------------------------------------
!
! Exceptions :  sun ray directly N,S,E,W
!
! If the sun is 'exactly' at the N,S,E or W
! the algorithm must be changed 
! e.g. to avoid the evaluation of tan(pi/2)
!
!-----------------------------------------------------------------------

  IF(myproc == 0) THEN

  IF (( ABS(sazimuth-0) < 0.0001) .OR. (ABS(sazimuth-pi2) < 0.0001) .OR. &
      (ABS(sazimuth-pi) < 0.0001) .OR. (ABS(sazimuth-pi32) < 0.0001) ) THEN

    DO i = II,IF,dI     ! we will test each grid point starting from the
      DO j = JI,JF,dJ   ! farthest node form the sun

         xg=xlg(i)      ! x and y coordinates of the test point along the line
         yg=ylg(j)      ! 

         Sx=xlg(i) + a*dx*D + a*dx*F ! x and y coordinate of the first
         Sy=ylg(j) + b*dy*E + b*dy*H ! intersection of the line going toward the sun
                                   ! and the horizontal grid

         DO l=1,(2*nxlg+2*nylg)    ! incrementation along the line
                                   ! (2*nxlg+2*nylg) is choosen so that we are certain
                                   ! to explore the whole domain
            IF ((xg >= xlg(nxlg)) .OR. (yg >= ylg(nylg))                   &
                .OR. (xg <= xlg(1)) .OR. (yg <= ylg(1)) ) EXIT
              ! if we are outside the topography, then :END
            ztest = hterainlg(i,j)                                         &
                  + SQRT((Sx-xlg(i))*(Sx-xlg(i))+(Sy-ylg(j))*(Sy-ylg(j)))  &
                    * TAN(saltitude)           
              ! ztest is the altitude of the point to be tested
              ! as the altitude of the initial point (hterain(i,j)
              ! plus the elevation gain along the line
              ! as the distance between the initial point
              ! and the current testing point times tan(saltitude)

            xh = 2 + NINT((xg + (D+F)*a*dx) / dx) ! find the x and y coordinates
            yh = 2 + NINT((yg + (E+H)*b*dy) / dy) ! of the point to be tested
                                                  ! xg and yg may not be on the grid !
            IF (xh < 1 .OR. xh > nxlg .OR. yh < 1 .OR. yh > nylg) EXIT ! EMK NEW
                                                   
            htest = hterainlg(xh,yh)
            IF (ztest < htest) THEN      ! perform the test
               shlg(i,j) = 0.0
               EXIT
            END IF         
            xg = xg + a*dx*(D+F)         !if not shaded, increment xg and yg
            yg = yg + b*dy*(E+H)         ! by going toward the sun
            Sx = Sx + a*dx*(D+F)
            Sy = Sy + b*dy*(E+H)
         END DO
      END DO
   END DO
    
!-----------------------------------------------------------------------
!
!  General case: NE,NW,SE,SW
!  if the sun is not exactly N,S,E,W
!
!-----------------------------------------------------------------------

  ELSE 
    DO i=II,IF,dI 
      DO j=JI,JF,dJ

         xg = xlg(i)  ! x and y coordinates of the first point along the line
         yg = ylg(j)

         Sx = xlg(i) + a*(b*(yg+b*dy)-b*ylg(j))*((D+F)*TAN(C*pi2-sazimuth)  &
            + (E+H)*TAN(sazimuth-(C-1)*pi2))  
         Sy = ylg(j) + b*(a*(xg+a*dx)-a*xlg(i))*((E+H)*TAN(C*pi2-sazimuth)  &
            + (D+F)*TAN(sazimuth-(C-1)*pi2))
          ! x and y coordinate of the first intersection
          ! between the line and the grid 

         DO l=1,(2*nxlg+2*nylg)  ! incrementation (See above)
           IF ((xg >= xlg(nxlg)) .OR. (yg >= ylg(nylg))                     & 
               .OR. (xg <= xlg(1)) .OR. (yg <= ylg(1)) ) EXIT   
             ! if we are outside the grid: STOP

           IF (( abs(Sy-(yg+b*dy)) < 0.01) .AND.                            &
               ( abs(Sx-(xg+a*dx)) < 0.01) ) THEN
               ! if the next testing point is very close to a grid point
               ! the test wil be performed with this grid point     
              ztest = hterainlg(i,j)                                        &
                    + SQRT((Sx-xlg(i))*(Sx-xlg(i))+(Sy-ylg(j))*(Sy-ylg(j))) &
                      *TAN(saltitude)
               ! altitude of the test point along the line
              xh = 2+NINT((xg+a*dx)/dx) ! x and y coordinates
              yh = 2+NINT((yg+b*dy)/dy) ! of the point where the test will be performed

              IF (xh < 1 .OR. xh > nxlg .OR. yh < 1 .OR. yh > nylg) EXIT ! EMK NEW

              htest = hterainlg(xh,yh)  ! altitude of this point

              IF (ztest < htest) THEN   ! test itself (see above)
                 shlg(i,j) = 0.0
                 EXIT
              END IF

              xg = xg+a*dx   ! incrementation, the point xg,yg has been tested
              yg = yg+b*dy   ! we can continue to the next one

              Sx = xlg(i) + a*(b*(yg+b*dy)-b*ylg(j))*((D+F)*TAN(C*pi2-sazimuth) &
                 + (E+H)*TAN(sazimuth-(C-1.0)*pi2))
              Sy = ylg(j) + b*(a*(xg+a*dx)-a*xlg(i))*((E+H)*TAN(C*pi2-sazimuth) &
                 + (D+F)*TAN(sazimuth-(C-1.0)*pi2))

                ! x and y's of the next point to be tested

           ELSE IF ( abs(Sy-yg) >dy ) THEN
                   ! if the line hits first a horizontal grid line
                   ! the test will be peformed on the horizontal
                   ! grid cell boundary
                      
              ztest = hterainlg(i,j)                                         &
                    + SQRT((Sx-xlg(i))*(Sx-xlg(i))                           &
                            +((yg+b*dy)-ylg(j))*((yg+b*dy)-ylg(j)))          &
                      *TAN(saltitude)
                   ! altitude along the line
              xh = 2+NINT((xg+a*dx)/dx) ! coordinates of the neighboring 
              yh = 2+NINT((yg+b*dy)/dy) ! test points
              IF (xh < 1 .OR. xh > nxlg .OR. yh < 1 .OR. yh > nylg) EXIT ! EMK NEW
              xhh = 2+NINT(xg/dx)
              yhh = 2+NINT((yg+b*dy)/dy)
              IF (xhh < 1 .OR. xhh > nxlg .OR. yhh < 1 .OR. yhh > nylg) EXIT ! EMK NEW
              htest = (a*(Sx-xg)*hterainlg(xh,yh)                            &
                       +a*((xg+a*dx)-Sx)*hterainlg(xhh,yhh)) / dx
                   ! altitude of the topography
                   ! linearly interpolated along the grid cell boundary
              IF (ztest < htest) THEN   ! test itself
                 shlg(i,j) = 0.0
                 EXIT
              END IF

              xg = xg      ! incrementation since we hit first a horizontal 
              yg = yg+b*dy ! grid cell boundary, we increment only y not x

              Sx = xlg(i)                                                   &
                 + a*(b*(yg+b*dy)-b*ylg(j))                                 &
                   *((D+F)*TAN(C*pi2-sazimuth)+(E+H)*TAN(sazimuth-(C-1.0)*pi2))
              Sy = ylg(j)                                                   &
                 + b*(a*(xg+a*dx)-a*xlg(i))                                 &
                   *((E+H)*TAN(C*pi2-sazimuth)+(D+F)*TAN(sazimuth-(C-1.0)*pi2))
                ! x and y of the next testing point

           ELSE IF (  abs(Sy-yg) < dy ) THEN ! see the previous case for comments

              ztest = hterainlg(i,j)                                        &
                    + SQRT(((xg+a*dx)-xlg(i))*((xg+a*dx)-xlg(i))            &
                           +(Sy-ylg(j))*(Sy-ylg(j)))*TAN(saltitude)

              xh = 2+NINT((xg+a*dx)/dx)
              yh = 2+NINT((yg+b*dy)/dy)
              IF (xh < 1 .OR. xh > nxlg .OR. yh < 1 .OR. yh > nylg) EXIT ! EMK NEW
              xhh = 2+NINT((xg+a*dx)/dx)
              yhh = 2+NINT(yg/dy)
              IF (xhh < 1 .OR. xhh > nxlg .OR. yhh < 1 .OR. yhh > nylg) EXIT ! EMK NEW

              htest = (b*(Sy-yg)*hterainlg(xh,yh)+b*((yg+b*dy)-Sy)            &
                                          *hterainlg(xhh,yhh)) / dy

              IF (ztest < htest) THEN
                 shlg(i,j) = 0.0
                 EXIT
              END IF
              yg = yg
              xg = xg+a*dx

              Sx = xlg(i) + a*(b*(yg+b*dy)-b*ylg(j))*((D+F)*TAN(C*pi2-sazimuth) &
                      +(E+H)*TAN(sazimuth-(C-1.0)*pi2))
              Sy = ylg(j) + b*(a*(xg+a*dx)-a*xlg(i))*((E+H)*TAN(C*pi2-sazimuth) &
                                        +(D+F)*TAN(sazimuth-(C-1.0)*pi2))
                
           END IF

        END DO    ! l
      END DO   ! j
    END DO     ! i

  END IF   
!-----------------------------------------------------------------------
!
! Extrapolation in the direction of the sun
! part of the domain where the shade can not be computed
! the shade is set to its value at the neighboring point
!
!-----------------------------------------------------------------------

   IF ( C==1 ) THEN
      DO i=2,nxlg-3
         shlg(i,nylg-2)=shlg(i,nylg-3)
      END DO
      DO j=2,nylg-3
        shlg(nxlg-2,j)=shlg(nxlg-3,j)
      END DO
      shlg(nxlg-2,nylg-2)=shlg(nxlg-3,nylg-3)
   END IF
  
   IF (C==2) THEN
      DO i=3,nxlg-2
         shlg(i,nylg-2)=shlg(i,nylg-3)
      END DO
      DO j=2,nylg-3
         shlg(2,j)=shlg(3,j)
      END DO
      shlg(2,nylg-2)=shlg(3,nylg-3)
   END IF
  
   IF (C==3) THEN
      DO i=3,nxlg-2
         shlg(i,2)=shlg(i,3)
      END DO
      DO j=3,nylg-2 
         shlg(2,j)=shlg(3,j)
      END DO
      shlg(2,2)=shlg(3,3)
   END IF

   IF (C==4) THEN
      DO i=2,nxlg-3
         shlg(i,2)=shlg(i,3)
      END DO
      DO j=3,nylg-2
         shlg(nxlg-2,j)=shlg(nxlg-3,j)
      END DO
      shlg(nxlg-2,2)=shlg(nxlg-3,3)
   END IF

! Extrapolation at the boundaries, outside the actual grid used in arps 
! i,j=1 and i,j=nxlg (or nylg)
! first for i,j=1  and i,j=N-1
! second for i,j = N

  DO i=2,nxlg-2
    shlg(i,1)=shlg(i,2)
    shlg(i,nylg-1)=shlg(i,nylg-2)
  END DO
  
  DO j=2,nylg-2
    shlg(1,j)=shlg(2,j)
    shlg(nxlg-1,j)=shlg(nxlg-2,j)
  END DO
  
  shlg(1,1)=shlg(2,2)
  shlg(1,nylg-1)=shlg(2,nylg-2)
  shlg(nxlg-1,1)=shlg(nxlg-2,2)
  shlg(nxlg-1,nylg-1)=shlg(nxlg-2,nylg-2)
  DO i=1,nxlg-1
    shlg(i,nylg)=shlg(i,nylg-1)
  END DO
  DO j=1,nylg-1
    shlg(nxlg,j)=shlg(nxlg-1,j)
  END DO
  
  shlg(nxlg,nylg)=shlg(nxlg-1,nylg-1)

! If the topography is two-dimensional in a 3D run(e.g. a periodic valley)
! the shlgade is set to be two-dimensionnal by setting it to
! its value at the middle of the domain

  IF (radshade == 2) THEN
    l=NINT(0.5*nylg)
    DO j=1,nylg
      DO i=1,nxlg
        shlg(i,j)=shlg(i,l)
      END DO
    END DO
  END IF
  
  END IF     ! myproc == 0

  CALL mpisplit2d(shlg,nx,ny,sh)

  IF(myproc == 0 ) print*,'The shade was computed.'
  
  DEALLOCATE(xlg,ylg,hterainlg,shlg)

RETURN

END SUBROUTINE shade            
