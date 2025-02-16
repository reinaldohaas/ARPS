!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SOILEBM                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE soilebm(nx,ny,nz,soiltyp,vegtyp,lai,veg,                     &
           tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,                    &
           qvsfc,windsp,psfc,rhoa,precip,                               &
           tair,qvair,cdha,cdqa,radsw,rnflx,shflx,lhflx,gflx,ct,        &
           evaprg,evaprtr,evaprr,qvsat,qvsata,f34)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Predict the soil surface temperature and moisture contents by solving
!  the surface energy and moisture budget equtions:
!
!  1. the ground surface temperature, Ts -- tsfc
!  2. the deep ground temperature,    T2 -- tdeep
!  3. the surface soil moisture,      wg -- wetsfc
!  4. the deep soil moisture,         w2 -- wetdp
!  5. the canopy moisture,            wr -- wetcanp
!
!-----------------------------------------------------------------------
!
!  The equations are listed as follows.
!
!
!       d Ts                        2PI
!      ------ = Ct (Rn - H - LE) - ----- (Ts - T2)
!       d t                         Tau
!
!
!       d T2      1
!      ------ = ----- (Ts - T2)
!       d t      Tau
!
!
!       d Wg      C1                 C2
!      ------ = ------ (Pg - Eg) - ----- (Wg - Wgeq)
!       d t     ROw d1              Tau
!
!
!       d W2      1
!      ------ = ------ (Pg - Eg - Etr)
!       d t     ROw d2
!
!
!       d Wr
!      ------ = Veg P - Er
!       d t
!
!
!  where
!
!      Tau    -- 1 day in seconds = 86400 seconds
!      PI     -- number of PI = 3.141592654
!      ROw    -- Density of liquid water
!      d1     -- Top layer depth of soil column, 0.01 m
!      d2     -- Deep layer depth of soil column, 1 m
!      Veg    -- Vegetation fraction
!      Ct     -- Thermal capacity
!      Rn     -- Radiation flux, rnflx
!      H      -- Sensible heat flux, shflx
!      LE     -- Latent heat flux, lhflx = latent*(Eg + Ev)
!      Eg     -- Evaporation from ground
!      Ev     -- Evapotranspiration from vegetation, Ev = Etr + Er
!      Etr    -- Transpiration of the remaining part of the leaves
!      Er     -- Evaporation directly from the foliage covered by
!                intercepted water
!      P      -- Precipitation rates
!      Pg     -- Precipitation reaching the ground,
!                Pg = (1 - Veg) P
!      Wgeq   -- Surface volumetric moisture
!      C1, C2 -- Coefficients
!
!  For detailed information about the surface energy budget model,
!  see the articles in the reference list.
!
!  The second-order Rouge-Kutta time integration scheme is used,
!  which is described below.
!
!  Assume a equation in the form of
!
!      d X
!     ----- = F(X, t)
!      d t
!
!  In the forward scheme, we have
!
!      X(1) = X(0) + dt * F[X(0), t0]
!
!  We split one time step into two halves, dt2 = dt/2, and use the
!  forward scheme to calculate the first half step X(1/2).
!
!      X(1/2)  = X(0) + dt2 * F[X(0), t0]
!
!  Then we can calculate the Right Hand Side (RHS) of the equation at
!  the half step, F[X(1/2), t(1/2)]. Finally, we calculate the one
!  step prediction, X(t1), by use of the average of F[X(0), t0] and
!  F[X(1/2), t(1/2)].
!
!      X(1) = X(0) + dt * 0.5 * { F[X(0),t0] + F[X(1/2), t(1/2)] }
!  REFERENCES:
!
!  Jacquemin, B. and J. Noilhan, 1990: Sensitivity Study and
!       Validation of a Land Surface Parameterization Using the
!       HAPEX-MOBILHP Data Set, Boundary-Layer Meteorology, 52,
!       93-134, (JN).
!
!  Noilhan, J. and S. Planton, 1989: A Simple Parameterization of
!       Land Surface Processes for Meteorological Model, Mon. Wea.
!       Rev., 117, 536-549, (NP).
!
!  Pleim, J. E. and A. Xiu, 1993: Development and Testing of a
!       Land-Surface and PBL Model with Explicit Soil Moisture
!       Parameterization, Preprints, Conf. Hydroclimat., AMS, 45-51,
!       (PX).
!
!  Bougeault, P., et al., 1991: An Experiment with an Advanced
!       Surface Parameterization in a Mesobeta-Scale Model. Part I:
!       Implementation, Monthly Weather Review, 119, 2358-2373.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  11/16/93
!
!  MODIFICATION HISTORY:
!
!  10/30/94 (Y. Liu)
!  Fixed a bug reported by Richard Carpenter.
!
!  11/01/1994 (Y. Liu)
!  Subroutine name of soil model were changed from SFCEBM to SOILEBM
!  and the arguments passed were changed to 2-D arrays instead of 3-D
!  temporary arrays.
!
!  12/12/1994 (Y. Liu)
!  Fixed a bug for the final calcultaion of qvsfc.
!
!  12/14/1994 (Y. Liu)
!  Fixed a bug in phycst.inc for the water density value which was
!  previously mistakingly set to 1 kg/m**3. The correct value should
!  be 1000 kg/m**3. This bug largely influenced the integration of Wg
!  and W2.
!
!  12/23/1994 (Y. Liu)
!  Added the runoff calculation for Wr.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the soil model and
!  at the same time delete the table data array veg(14).
!
!  03/27/1995 (Yuhe Liu)
!  Changed the solor radiation used in the calculation of surface
!  resistence factor F1 from the one at the top of atmosphere to the
!  one at the surface.
!
!  Changed the formula of calculating the surface resistence factor
!  F3 to F3=1, instead of varying with qvsat(Tair) and qvair.
!
!  12/8/1998 (Donghai Wang and Vince Wong)
!  Added a new 2-D permanent array, snowcvr(nx,ny), for snow cover.
!  We just used a simple scheme to consider the snow cover process.
!
!  2000/01/10 (Gene Bassett)
!  Snow cover (0 or 1) changed to snow depth (snowdpth).  For simplicity
!  a fractional value for snow cover is not used (simply say the grid
!  point is completely covered with snow if snowdpth > snowdepth_crit,
!  otherwise no snow).
!
!  2000/02/04 (Gene Bassett, Yang Kun)
!  Fixed an error in tsoil integration (rhst2) causing tsoil to change
!  by only 50% of what it should.
!
!  2001/12/07 (Diandong Ren, Ming Xue) 
!  Re-structured the code, moved the calculations of the right hand
!  side terms of the soil-vegetation model into subroutine soilebm_frc.
!
!  Soil seasonal temperature trend to be added.
!
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!  2002/06/9 (Ming Xue)
!  Revoked some modifications to the soil moisture related caculations
!  that Diandong Ren put in since IHOP_2 - the mods need more testing.
!
!  2002/12/13 (Jerry Brotzge) 
!  Updated code to match recommendations by Pleim and Xiu (1995) and 
!  Xiu and Pleim (2001 - JAM).  
!
!-----------------------------------------------------------------------
! 
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (sfc/top)
!
!    soiltyp  Soil type at the horizontal grid points
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    windsp   Wind speed just above the surface (m/s)
!    psfc     Surface pressure (Pascal)
!    rhoa     Near sfc air density
!    precip   Precipitation flux reaching the surface (m/s)
!
!    cdha     Array for cdh, surface drag coefficient for heat
!    cdqa     Array for cdq, surface drag coefficient for moisture
!
!    pres     3-dimensional pressure
!    temp     3-dimensional temperature
!    qv       3-dimensional specific humidity
!
!  INPUT/OUTPUT: 
!
!    tsfc     Temperature at ground surface (K)
!    tdeep    Deep soil temperature (K)
!    wetsfc   Surface soil moisture
!    wetdp    Deep soil moisture
!    wetcanp  Canopy water amount
!
!  OUTPUT:
!
!    qvsfc    Effective S. H. at sfc.
!
!  Local automatic work arrays for storing the right hand forcing terms 
!  of the soil model equations and for storing intermediate values of the
!  soil state variables. 
!
!    frc_tsfc   Temporary array
!    frc_tdeep  Temporary array
!    frc_wsfc   Temporary array
!    frc_wdp    Temporary array
!    frc_wcnp   Temporary array
!
!    tsfcn      Temporary array
!    tdeepn     Temporary array
!    wgn        Temporary array
!    w2n        Temporary array
!    wrn        Temporary array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
  INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
  REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
  REAL :: wetsfc (nx,ny)        ! Surface soil moisture
  REAL :: wetdp  (nx,ny)        ! Deep soil moisture
  REAL :: wetcanp(nx,ny)        ! Canopy water amount
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface

  REAL :: snowdpth(nx,ny)       ! Snow depth (m)

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

  REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
  REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
  REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
  REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
  REAL :: ct     (nx,ny)        ! Soil thermal coefficient

  REAL :: evaprg (nx,ny)        ! Evaporation
  REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
  REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
  REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
  REAL :: qvsat  (nx,ny)        !


  REAL :: frc_tsfc(nx,ny)       ! Right hand side forcing for tsfc eq.
  REAL :: frc_tdeep(nx,ny)      ! Right hand side forcing for tdeep eq.
  REAL :: frc_wsfc(nx,ny)       ! Right hand side forcing for wetsfc eq.
  REAL :: frc_wdp(nx,ny)        ! Right hand side forcing for wetsp eq.
  REAL :: frc_wcnp(nx,ny)       ! Right hand side forcing for wetcanp eq.

  REAL :: tsfcn(nx,ny)          ! Temporary array, tsn
  REAL :: tdeepn(nx,ny)         ! Temporary array, t2n
  REAL :: wgn(nx,ny)            ! Temporary array, wetsfcNEW
  REAL :: w2n(nx,ny)            ! Temporary array, w2n
  REAL :: wrn(nx,ny)            ! Temporary array, wrn
  REAL :: relief          ! Difference between seasonal average skin and deep soil 
!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
!  Parameters and variables are defined in globcst.inc:
!
!    dtsfc      Surface model time step
!    nsfcst     # of surface model time steps
!
!    moist      Moist flag
!
!    year       Reference year
!    month      Reference month
!    day        Reference day
!    jday       Reference Julian day
!    hour       Hour of reference time
!    minute     Minute of reference time
!    second     Second of reference time
!
!    latitud    Latitude at the domain center
!    longitud   Longitude at the domain center
!
!    curtim     Current model time
!    dtbig      Length of big time step
!
!    bslope     Slope of the retention curve
!    cgsat      Soil thermal coefficient for bare ground at saturation
!    cgv        Soil thermal coef. for totally shielded ground by veg.
!    pwgeq      Coefficient of Wgeq formula. NP, Tab. 2
!    awgeq      Coefficient of Wgeq formula. NP, Tab. 2
!    c1sat      Value of C1 at saturation. NP, Tab. 2
!    c2ref      Value of C2 for W2 = .5 * Wsat. NP, Tab. 2
!    wsat       Saturated volumetric moisture content. JN, Tab. 1
!    wfc        Field capacity moisture. JN, Tab. 1
!    wwlt       Wilting volumetric moisture content. JN, Tab. 1
!
!  Parameters and variables are defined in phycst.inc:
!
!    solarc     Solar constant (W/m**2)
!    emissg     Emissivity of the ground
!    emissa     Emissivity of the atmosphere
!    sbcst      Stefen-Boltzmann constant
!
!    rhow       Liquid water reference density (kg/m**3)
!    rd         Gas constant for dry air (kg/(m s**2))
!    cp         Gas heat capacity at constant pressure
!    cv         Gas heat capacity at constant volume
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: pi                         ! Pi
  PARAMETER (pi = 3.141592654)

!  REAL :: tau                        ! Seconds of a day = 24. * 3600.
!  PARAMETER (tau = 86400.)

  REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                      ! dtsfc2 = dtsfc/2.

  REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun

  REAL :: cg          ! Soil thermal coefficient for bare ground

  REAL :: rhsts(nx,ny)       ! Right hand side of Eq. for Ts at current time
  REAL :: rhst2(nx,ny)       ! Right hand side of Eq. for T2 at current time
  REAL :: rhswg(nx,ny)       ! Right hand side of Eq. for Wg at current time
  REAL :: rhsw2(nx,ny)       ! Right hand side of Eq. for W2 at current time
  REAL :: rhswr(nx,ny)       ! Right hand side of Eq. for Wr at current time
  REAL :: wrmax(nx,ny)       ! Maximum value for canopy moisture, wetcanp
  REAL :: c1wg        ! Coefficient in the surface moisture Eq. of Wg
  REAL :: c2wg        ! Coefficient in the surface moisture Eq. of Wg

  REAL :: wgeq        ! Surface moisture when gravity balances the capillarity

  REAL :: wr2max      ! Tendency to reach the maximum wrmax
  REAL :: runoff      ! Runoff of the interception reservoir.
  REAL :: pnet        ! Residual of precip. and evap.
  REAL :: vegp        ! Precip. intercepted by vegetation

  INTEGER :: i, j, it

  REAL :: tema

  LOGICAL :: firstcall        ! First call flag of this subroutine

!  SAVE firstcall, log100, dtsfc2, tauinv
  SAVE firstcall, log100, dtsfc2
  DATA firstcall/.true./

  INTEGER :: jday_min         ! offset value from Jan 01. 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (firstcall) THEN
    log100    = ALOG(100.0)
    dtsfc2    = dtsfc/2.0
    firstcall = .false.
  END IF

  IF ( moist /= 0 ) THEN

!  Calculate saturated specific humidity near the surface, qvsata.

    CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tair,qvsata)

    DO j = 1, ny-1
    DO i = 1, nx-1
!      f34(i,j) = MAX( 0., 1.0 - 0.0016 * (298.0-tair(i,j))**2 )

       IF (tair(i,j) >= 302.15) THEN                    !(XP2001) 
         f34(i,j) = 1.0 / (1.0 + EXP( 1.18 * (tair(i,j)-314.00)))
       ELSE
         f34(i,j) = 1.0 / (1.0 + EXP( -0.41* (tair(i,j)-282.05)))
       ENDIF 

    END DO
    END DO
  END IF

  jday_min = 61       ! may change with latitude
  relief=tsoil_offset_amplitude*sin((jday-jday_min)/365.0*2.0*PI)

!
!-----------------------------------------------------------------------
!
!  Start time integration loop
!
!-----------------------------------------------------------------------
!
  DO it = 1, nsfcst               

    CALL soilebm_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,                   &
         tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,                      &
         qvsfc,windsp,psfc,rhoa,precip,tair,qvair,cdha,cdqa,            &
         radsw,rnflx,shflx,lhflx,gflx,ct,evaprg,evaprtr,                &
         evaprr,qvsat,qvsata,f34,                                       &
         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief)  

    DO j = 1, ny-1
    DO i = 1, nx-1
      tsfcn (i,j) = tsfc(i,j) + dtsfc2 * frc_tsfc(i,j)
      tdeepn(i,j) = tdeep(i,j)+ dtsfc2 * frc_tdeep(i,j)
    END DO
    END DO

    IF ( moist /= 0 ) THEN    

      DO j = 1, ny-1
      DO i = 1, nx-1
        wgn(i,j) = wetsfc(i,j) + dtsfc2 * frc_wsfc(i,j)
        wgn(i,j) = MAX(wgn(i,j), 0.0 )
        wgn(i,j) = MIN(wgn(i,j), wsat(soiltyp(i,j)) )

        w2n(i,j) = wetdp(i,j) + dtsfc2 * frc_wdp(i,j)
        w2n(i,j) = MAX( w2n(i,j), 0.0 )
        w2n(i,j) = MIN(w2n(i,j), wsat(soiltyp(i,j)) )

        wrn(i,j) = wetcanp(i,j) + dtsfc2 * frc_wcnp(i,j)
        wrn(i,j) = MAX(wrn(i,j), 0.0 )
        wrn(i,j) = MIN(wrn(i,j), wrmax(i,j) )
      END DO
      END DO

    ELSE

      DO j = 1, ny-1
        DO i = 1,nx-1
          w2n(i,j) = wetdp(i,j)
        END DO
      END DO

    END IF

    DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        rhsts(i,j) = 0.0
        rhst2(i,j) = 0.0
      ELSE
        rhsts(i,j)=frc_tsfc(i,j)
        rhst2(i,j)=frc_tdeep(i,j)
      END IF
    END DO
    END DO

    IF ( moist /= 0 ) THEN    
      DO j = 1, ny-1
      DO i = 1, nx-1
        rhswg(i,j)=frc_wsfc(i,j)
        rhsw2(i,j)=frc_wdp(i,j)
        rhswr(i,j)=frc_wcnp(i,j)
      END DO
      END DO
    END IF

    CALL soilebm_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,                   &
         tsfcn,tdeepn,wgn,w2n,wrn,snowdpth,                             &
         qvsfc,windsp,psfc,rhoa,precip,tair,qvair,cdha,cdqa,            &
         radsw,rnflx,shflx,lhflx,gflx,ct,                               &
         evaprg,evaprtr,evaprr,qvsat,qvsata,f34,                        &
         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief)  

!  Integration for Ts and T2 at one time step.

    DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        rhsts(i,j) = 0.0
        rhst2(i,j) = 0.0
      ELSE
        rhsts(i,j) = 0.5 * (frc_tsfc(i,j) +rhsts(i,j))
        rhst2(i,j) = 0.5 * (frc_tdeep(i,j)+rhst2(i,j)) 
      END IF
        tsfc(i,j)  = tsfc(i,j) + dtsfc * rhsts(i,j)
        tdeep(i,j) = tdeep(i,j)+ dtsfc * rhst2(i,j)
    END DO
    END DO

    IF ( moist /= 0 ) THEN
      DO j = 1, ny-1
      DO i = 1, nx-1
          IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.         &
            snowdpth(i,j) >= snowdepth_crit ) THEN
            rhswg(i,j) = 0.0
            rhsw2(i,j) = 0.0
            rhswr(i,j) = 0.0
          ELSE
            rhswg(i,j) = 0.5 * (frc_wsfc(i,j)+rhswg(i,j))
            rhsw2(i,j) = 0.5 * (frc_wdp (i,j)+rhsw2(i,j))
            rhswr(i,j) = 0.5 * (frc_wcnp(i,j)+rhswr(i,j))
          END IF
          wetsfc(i,j) = wetsfc(i,j) + dtsfc * rhswg(i,j)
          wetsfc(i,j) = MAX( wetsfc(i,j), 0.0 )
          wetsfc(i,j) = MIN( wetsfc(i,j), wsat(soiltyp(i,j)) )

          wetdp(i,j) = wetdp(i,j) + dtsfc * rhsw2(i,j)
          wetdp(i,j) = MAX( wetdp(i,j), 0.0 )
          wetdp(i,j) = MIN( wetdp(i,j), wsat(soiltyp(i,j)) )

          wetcanp(i,j) = wetcanp(i,j) + dtsfc * rhswr(i,j)
          wetcanp(i,j) = MAX( wetcanp(i,j), 0.0 )
          wetcanp(i,j) = MIN( wetcanp(i,j), wrmax(i,j) )
       END DO
      END DO
    END IF

  END DO  ! TIME INTEGRATION

  DO j = 1, ny-1                               !SOIL
  DO i = 1, nx-1
    tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

    IF (snowdpth(i,j) >= snowdepth_crit) THEN
      ct(i,j) = cg_snow
      gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)                 &
                  *snowflxfac/(tau*ct(i,j))
                                        ! Snow cover
    ELSE
      cg = cgsat(soiltyp(i,j))                                        &
          * ( wsat(soiltyp(i,j))/tema )                               &
          **( bslope(soiltyp(i,j))/log100 )

!      ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv                     &
!                           + veg(i,j) * cg )       ! NP, Eq. 8

      ct(i,j) = cg                                  !(PX1995) 

      gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))
                                        ! Ground diffusive heat flux
    END IF

      shflx(i,j) =rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)             &
         * ( tsfc(i,j) - tair(i,j) *(psfc(i,j)/1.0E5)**rddcp )         !RDDRDD


  END DO
  END DO

  CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)             !HYDROLOGY

  CALL evapflx(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,            &
       soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp,wetcanp,               &
       snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat)

  DO j=1, ny-1
  DO i=1, nx-1
    qvsfc(i,j) = lhflx(i,j) + qvair(i,j)
    evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
    evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
    evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
    lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)           !MOISTURE
    lhflx  (i,j) = lhflx(i,j) * lathv                                   !QE 
  END DO
  END DO

  RETURN
END SUBROUTINE soilebm


!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SOILEBM_FRC               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE soilebm_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,           &
     tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,qvsfc,windsp,psfc,  &
     rhoa,precip,tair,qvair, cdha,cdqa,radsw, rnflx,shflx,lhflx,  &
     gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata,f34,              &
     frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief)  

!------------------------------------------------------------------
!                                                                  
!  AUTHOR: Diandong Ren (dd_ren@rossby.metr.ou.edu)                
!          Ming Xue     (mxue@ou.edu)
!  12/08/2001                                                      
!                                                                  
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!------------------------------------------------------------------
!                                                                  
!  PURPOSE:                                                        
!  To calculate the right hand side forcing terms in the 
!  soil-vegetation model.
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz
!  REAL :: eta_fac,wmax,c1max,sigma_sqr

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
  INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
  REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
  REAL :: wetsfc (nx,ny)        ! Surface soil moisture
  REAL :: wetdp  (nx,ny)        ! Deep soil moisture
  REAL :: wetcanp(nx,ny)        ! Canopy water amount
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface

  REAL :: snowdpth(nx,ny)       ! Snow depth (m)

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

  REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
  REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
  REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
  REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
  REAL :: ct     (nx,ny)        ! Soil thermal coefficient

  REAL :: evaprg (nx,ny)        ! Evaporation
  REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
  REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
  REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
  REAL :: qvsat  (nx,ny)        !

  REAL :: frc_tsfc(nx,ny)      ! Right hand side forcing for tsfc eq.
  REAL :: frc_tdeep(nx,ny)       ! Right hand side forcing for tsoil eq.
  REAL :: frc_wsfc(nx,ny)       ! Right hand side forcing for wetsfc eq. 
  REAL :: frc_wdp(nx,ny)        ! Right hand side forcing for wetsp eq. 
  REAL :: frc_wcnp(nx,ny)       ! Right hand side forcing for wetcanp eq.
  REAL :: wrmax(nx, ny)         ! Maximum value for canopy moisture, wetcanp

  REAL :: c1wg        ! Coefficient in the surface moisture Eq. of Wg
  REAL :: c2wg        ! Coefficient in the surface moisture Eq. of Wg

  REAL :: pi          ! Pi

  PARAMETER (pi = 3.141592654)

!  REAL :: tau         ! Seconds of a day = 24. * 3600.
!  PARAMETER (tau = 86400.)

  REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                      ! dtsfc2 = dtsfc/2.
  REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun
  REAL :: cg          ! Soil thermal coefficient for bare ground
  REAL :: wgeq        ! Surface moisture when gravity balances the capillarity
  REAL :: wr2max      ! Tendency to reach the maximum wrmax
  REAL :: runoff      ! Runoff of the interception reservoir.
  REAL :: pnet        ! Residual of precip. and evap.
  REAL :: vegp        ! Precip. intercepted by vegetation

  INTEGER :: i, j

  REAL :: tema
  REAL :: eta, c1max, wmax, sig2 


  REAL :: relief       ! Difference between seasonal average skin and deep soil 
                      ! temperature  (t_skin - t_deeplayer).
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  log100    = ALOG(100.0)
  dtsfc2    = dtsfc/2.0

  DO j = 1, ny-1
    DO i = 1, nx-1
      tema = MAX( wetdp(i,j), wwlt(soiltyp(i,j)) )
  
      IF (snowdpth(i,j) >= snowdepth_crit) THEN !snow cover
        ct(i,j) = cg_snow
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)                     &
                          *snowflxfac/(tau*ct(i,j))
      ELSE
        cg = cgsat(soiltyp(i,j))                                      &
            * ( wsat(soiltyp(i,j))/tema )                             &
            **( bslope(soiltyp(i,j))/log100 )
        ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv + veg(i,j) * cg )     
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))
      END IF
  
      shflx(i,j) = rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)           &
                  * ( tsfc(i,j) - tair(i,j) )  
    END DO
  END DO

  IF ( moist == 0 ) THEN

    DO j = 1, ny-1
      DO i = 1, nx-1
         evaprg(i,j) = 0.0    
         evaprtr(i,j) = 0.0   
         evaprr(i,j) = 0.0  
         lhflx(i,j) = 0.0   
      END DO
    END DO

  ELSE

    CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)

    CALL evapflx(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,         &
                   soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp,wetcanp,    &
                   snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat)
  END IF

  DO j = 1, ny-1
    DO i = 1, nx-1
      evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
      evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
      evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
      lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)
      lhflx(i,j) = lhflx(i,j) * lathv       ! Latent heat flux
    END DO
  END DO

  DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        frc_tsfc (i,j) = 0.
        frc_tdeep(i,j) = 0.
      ELSE
        frc_tsfc(i,j) = ct(i,j)*(rnflx(i,j)-shflx(i,j)-lhflx(i,j)-gflx(i,j))
        IF ( snowdpth(i,j) >= snowdepth_crit ) THEN
          frc_tdeep(i,j)= 1.03* (tsfc(i,j)-tdeep(i,j)-relief)*tauinv*snowflxfac ! Snow cover
        ELSE
          frc_tdeep(i,j)= 1.03*(tsfc(i,j) - tdeep(i,j)-relief) * tauinv      
        END IF
      END IF
    END DO
  END DO

  IF ( moist /= 0 ) THEN   

    DO j = 1, ny-1                                       !HYDROLOGY
      DO i = 1, nx-1
        wrmax(i,j) = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! meter
        IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.           &
             snowdpth(i,j) >= snowdepth_crit ) THEN
          frc_wsfc(i,j) = 0.
          frc_wdp(i,j)  = 0.
          frc_wcnp(i,j) = 0.
        ELSE

          tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

          c1wg = 0.4* c1sat(soiltyp(i,j))                               &
                 * ( wsat(soiltyp(i,j)) / tema )                        &
                 **( bslope(soiltyp(i,j)) / 2.0 + 1.0)

!--------------------------------------------------------------------------
! Replacement Cl to improve dry soils (NM1996 - A.3)   (JAB)
!--------------------------------------------------------------------------

          IF (wetsfc(i,j) < wwlt(soiltyp(i,j)) ) THEN

            eta = (-0.01815 * tsfc(i,j) + 6.41 ) * wwlt(soiltyp(i,j)) + &
                  ( 0.0065 * tsfc(i,j) - 1.4 )
      
            c1max = (1.19 * wwlt(soiltyp(i,j)) - 5.09)*0.01*tsfc(i,j)   &
                   +( 1.464 * wwlt(soiltyp(i,j)) + 17.86 )
      
            wmax = eta * wwlt(soiltyp(i,j))
      
            sig2 = - ( (2 * alog ( 0.01 / c1max )) / ( wmax*wmax) )
      
            c1wg = c1max * EXP ( -0.5* ( (wetsfc(i,j)-wmax)*(wetsfc(i,j)-wmax)*sig2 ) )
          ENDIF

!----------------------------------------------------------------------------

          c2wg = c2ref(soiltyp(i,j)) * wetdp(i,j)                       &
                / ( wsat(soiltyp(i,j)) - wetdp(i,j) + wetsml )
          wgeq = wetdp(i,j) - wsat(soiltyp(i,j))                        &
                  * awgeq(soiltyp(i,j))                                 &
                  * ( wetdp(i,j) / wsat(soiltyp(i,j)) )                 &
                  ** pwgeq(soiltyp(i,j))                                &
                  * ( 1.0 - ( wetdp(i,j) / wsat(soiltyp(i,j)) )         &
                  ** ( 8 * pwgeq(soiltyp(i,j)) ) )

          frc_wsfc(i,j) = c1wg * ( ( 1.0 - veg(i,j) )                   &
                                  * precip(i,j) - evaprg(i,j) )         &
                          /(rhow * d1)                                  &
                         - c2wg * ( wetsfc(i,j) - wgeq ) * tauinv
          frc_wdp(i,j) = ( ( 1.0 - veg(i,j) ) * precip(i,j)             &
                        - evaprg(i,j) - evaprtr(i,j) )                  &
                      / ( rhow * d2 )     
          wr2max = ( wrmax(i,j) - wetcanp(i,j) ) / dtsfc2    
          vegp = veg(i,j) * precip(i,j)
          pnet = vegp - evaprr(i,j)
          tema = pnet - wr2max * rhow
          runoff = MAX( tema, 0.0 )
          vegp = vegp - runoff
          frc_wcnp(i,j) = ( vegp - evaprr(i,j) ) / rhow 
        END IF

      END DO
    END DO           

  END IF

  RETURN
END SUBROUTINE soilebm_frc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE EVAPFLX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE evapflx(nx,ny, radsw, f34, cdqa,                             &
           windsp,psfc,rhoa,qvair,                                      &
           soiltyp,vegtyp,lai,veg,                                      &
           tsfc,wetsfc,wetdp,wetcanp,snowdpth,                          &
           evaprg,evaprtr,evaprr,qvflx,qvsat)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate:
!
!  1. Evaporation from ground surface,
!
!        evaprg = rhoa * cdq * windsp * evaprg'
!
!  where
!
!        evaprg' = (1.-veg) * (rhgs * qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  2. Direct evaporation from the fraction delta of the foliage covered
!     by intercepted water.
!
!        Er = rhoa * cdq * windsp * Er'
!
!        Er' = delta * veg * (qvsats-qvair)
!
!  3. Transpiration of the remaining part (1-delta) of leaves,
!
!        Etr = rhoa * cdq * windsp * Etr'
!
!  where
!
!        Etr' = veg * (1-delta) * Ra/(Ra+Rs) * (qvsats-qvair )
!
!  and Ra is aerodynamic resistance and Rs is the surface resistance
!
!                     1
!        Ra = ----------------
!               cdq * windsp
!
!                   Rsmin
!        Rs = ------------------
!              LAI*F1*F2*F3*F4
!
!
!               f + Rsmin/Rsmax
!        F1 = -------------------
!                   f + 1
!
!                   Rg    2
!        f = 0.55 ----- -----
!                  Rgl   LAI
!
!             -  1,                             Wfc < W2
!             |
!             |    W2 - Wwlt
!        F2 = -  ------------,              Wwlt <= W2 <= Wfc
!             |   Wfc - Wwlt
!             |
!             -  0,                             W2 < Wwlt
!
!
!               1-0.06*(qvsats-qvair),   qvsats-qvair <= 12.5 g/kg
!        F3 = {
!               0.25,                      otherwise
!
!
!        F4 = 1 - 0.0016 * (298-tair)**2
!
!  4. Water vapor flux, qvflx,
!
!        qvflx = rhoa * cdq * windsp * qvflx'
!
!        qvflx' = (Eg' + Etr' + Er') = (qvsfc - qvair)
!
!     where qvsfc is the effective surface specific humidity
!
!        (qvsfc - qvair) = (Eg' + Etr' + Er')
!
!  This subroutine will solve Eg', Etr', Er', and qvflx'
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  4/20/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    radsw    Incoming solar radiation
!
!    f34      f3*f4;
!             f3 = Fractional conductance of atmospheric vapor pressure
!             f4 = Fractional conductance of air temperature
!
!    cdqa     Array for surface drag coefficient for water vapor
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    tsfc     Surface soil temperature (K)
!    wetsfc   Surface soil moisture
!    wetdp    Deep soil moisture
!    wetcanp  Vegetation moisture
!
!    psfc     Surface pressure
!    qvair    Specific humidity near the surface
!    windsp     Wind speed near the surface
!    rhoa     Air density near the surface
!
!  OUTPUT:
!
!    evaprp   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evaprr   Direct evaporation from the fraction delta
!    qvflx    Water vapor flux
!
!  WORK ARRAY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny
!
  REAL :: radsw  (nx,ny)

  REAL :: f34    (nx,ny)

  REAL :: cdqa   (nx,ny)

  INTEGER :: soiltyp(nx,ny)
  INTEGER :: vegtyp (nx,ny)
  REAL :: lai    (nx,ny)
  REAL :: veg    (nx,ny)

  REAL :: tsfc   (nx,ny)
  REAL :: wetsfc (nx,ny)
  REAL :: wetdp  (nx,ny)
  REAL :: wetcanp(nx,ny)
  REAL :: snowdpth(nx,ny)
!
  REAL :: psfc   (nx,ny)
  REAL :: qvair  (nx,ny)
  REAL :: windsp   (nx,ny)
  REAL :: rhoa   (nx,ny)

  REAL :: evaprg (nx,ny)
  REAL :: evaprtr(nx,ny)
  REAL :: evaprr (nx,ny)
  REAL :: qvflx  (nx,ny)
  REAL :: qvsat  (nx,ny)
!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: pi
  PARAMETER ( pi = 3.141592654 )
!
  INTEGER :: i, j

  REAL :: wrmax       ! Maximum value for canopy moisture, wetcanp
  REAL :: rstcoef     ! Coefficient of resistance
  REAL :: delta
  REAL :: rhgs        ! R.H. at ground surface
!
  REAL :: tema
  REAL :: temb

  REAL :: pterm
  REAL :: mterm
  REAL :: waf         ! Available soil moisture fraction
  REAL :: bw          ! Half point of the available soil mstr frctn curve
  REAL :: ps          ! Shelter factor 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j = 1, ny-1
    DO i = 1, nx-1

      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
!
!-----------------------------------------------------------------------
!
!  Over water and ice, qvsfc should be saturated
!
!-----------------------------------------------------------------------
!
        qvflx(i,j) = qvsat(i,j) - qvair(i,j)

        evaprg (i,j) = qvflx(i,j)
        evaprtr(i,j) = 0.0
        evaprr (i,j) = 0.0

      ELSE                                  ! over land

        wrmax = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! in meter
!
!-----------------------------------------------------------------------
!
!  In order to calculate the qv flux at the surface, we need to
!  calculate some parameters like the resistance coefficient,
!
!      rstcoef = rsta / (rsts + rsta)
!
!  where rsta is the aerodynamic resistance and rsts is the surface
!  resistances.
!
!      rsta = 1 / ( cdq * Va )
!
!      rsts = rsmin(vtyp) / ( lai(vtyp) * f1 * f2 * f3 * f4 )
!
!  f3 * f4 is time-independent and has been calculated previously
!  and stored in f34(i,j)
!
!-----------------------------------------------------------------------
!
!  Calculate f1.
!
!        f = 0.55*(radsw/rgl(vtyp))*(2./lai(vtyp))       ! NP, Eq. 34
!        f1 = ( rsmin(vtyp)/rsmax + f ) / ( 1. + f )     ! NP, Eq. 34
!
!  Note: the incoming solar radiation radsw is stored in radsw(i,j).
!
!-----------------------------------------------------------------------
!
        IF ( lai(i,j) == 0. ) THEN        ! No vegetation, I/W
          rstcoef = 1.
        ELSE
!          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )               &
!               * ( 2.0 / lai(i,j) )

          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )               &
               * ( 2.0 )       !(XP2001 - Eq.8)   (JAB) 

          rstcoef = ( rsmin(vegtyp(i,j)) / rsmax + temb )               &
                  / ( 1.0 + temb )
        END IF
!
!-----------------------------------------------------------------------
!
!  Calculate f2 and f1*f2.
!
!-----------------------------------------------------------------------
!
!        pterm = .5 + SIGN( .5, wetdp(i,j) - wfc(soiltyp(i,j)) )
!        mterm = SIGN( .5, wetdp(i,j) - wwlt(soiltyp(i,j)) )              &
!              - SIGN( .5, wetdp(i,j) - wfc(soiltyp(i,j)) )
!
!        rstcoef = rstcoef * ( pterm + mterm                             &
!                * ( wetdp(i,j)-wwlt(soiltyp(i,j)) )                     &
!                / ( wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) ) )

        waf = ( wetdp(i,j)-wwlt(soiltyp(i,j)) )            &
                / ( wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )

        bw = wwlt(soiltyp(i,j)) + (wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )  &
                / 3.0

        rstcoef = rstcoef / ( 1.0 + EXP( -5.0 * (waf - bw) ) )
!             (XP2001 - Eq.9)                                       (JAB)
!
!-----------------------------------------------------------------------
!
!  Calculate lai*f1*f2*f3*f4 where f3*f4 is stored in f34(i,j).
!
!-----------------------------------------------------------------------
!
!        rstcoef = lai(i,j)*rstcoef*f34(i,j)      ! lai*f1*f2*f3*f4

        ps = 0.3 * lai(i,j) + 0.7                 ! XP2001             (JAB)
        rstcoef = lai(i,j)*rstcoef*f34(i,j)/ps    ! lai*f1*f2*f3*f4/ps (JAB)
!
!-----------------------------------------------------------------------
!
!  Calculate the resistance coefficient, rsta/(rsts+rsta)
!
!        rsts = rsmin(vtyp)/(lai(i,j)*f1*f2*f3*f4) ! Sfc. resistance
!        rsta = 1./(cdh*va)                 ! NP, between Eq. 32 & 33
!        rstcoef = rsta/(rsta+rsts)
!                = 1/(1+rsts/rsta)
!
!-----------------------------------------------------------------------
!
        tema = rsmin(vegtyp(i,j)) * cdqa(i,j) * windsp(i,j)

        IF ( ABS(rstcoef) > 1.0E-30 ) THEN
          rstcoef = 1.0 / (1.0 + tema/rstcoef)
        END IF
!
!-----------------------------------------------------------------------
!
!  1. evaprg'
!
!        evaprg' = (1.-veg) * (rhgs*qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  evaprg will be stored for current and future use to
!  calculate the latent heat flux and soil moisture transports.
!
!-----------------------------------------------------------------------
!
        pterm = .5 + SIGN( .5, wetsfc(i,j)-1.1*wfc(soiltyp(i,j)) )

        rhgs = pterm                                                    &
             + (1.-pterm) * ( 0.25 * ( 1.0 - COS( wetsfc(i,j)           &
                              * pi / (1.1*wfc(soiltyp(i,j))))) ** 2 )

!      IF (snowdpth(i,j) .ge. snowdepth_crit) rhgs=1.0 !Snow cover

!        evaprg(i,j) = ( 1.0 - veg(i,j) )                                &
!                    * ( rhgs * qvsat(i,j) - qvair(i,j) )

        evaprg(i,j) = ( 1.0 - veg(i,j) )                                &
                    *  rhgs * ( qvsat(i,j) - qvair(i,j) )  !XP2001 (JAB) 

!
!-----------------------------------------------------------------------
!
!  2. Transpiration of the remaining part (1-delta) of leaves, Etr',
!
!        Etr' = (1-delta) * veg
!             * Ra/(Ra+Rs) * ( qvsats - qvair )
!
!  3. Direct evaporation from the fraction delta, Er'
!
!        Er' = delta * veg * ( qvsats - qvair )
!
!  Er' and Etr' are stored for future use to calculate the latent heat
!  flux and soil moisture transports.
!
!-----------------------------------------------------------------------
!
        IF ( wrmax == 0.0 ) THEN
          delta = 0.0
        ELSE
          delta = ( wetcanp(i,j) / wrmax ) ** 0.66666667
        END IF

        pterm = .5 + SIGN( .5, qvsat(i,j) - qvair(i,j) )
        delta = pterm * delta + ( 1. - pterm )

        tema = veg(i,j) * ( qvsat(i,j) - qvair(i,j) )

        evaprtr(i,j) = ( 1.0 - delta ) * rstcoef * tema
        evaprr (i,j) = delta * tema
!
!-----------------------------------------------------------------------
!
!  4. Water vapor flux, qvflx',
!
!        qvflx' = evaprg' + evaprtr' + evaprr'
!
!  qvflx will be saved for the future use to calculate the latent
!  heat flux.
!
!-----------------------------------------------------------------------
!
        qvflx(i,j) = evaprg(i,j) + evaprtr(i,j) + evaprr(i,j)
      END IF
                                                   ! NP, expl. Eq. 26-27
    END DO
  END DO

  RETURN
END SUBROUTINE evapflx  

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE OUSOIL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ousoil(nx,ny,nz,nzsoil,zpsoil,j3soil,j3soilinv,soiltyp,      &
           vegtyp,lai,veg,tsoil,qsoil,wetcanp,snowdpth,qvsfc,           &
           windsp,psfc,rhoa,precip,tair,qvair,cdma,cdha,cdqa,           &
           radsw, rnflx, radswnet, radlwin, shflx, lhflx, gflx,         &
           ct,evaprg,evaprtr,evapcan,qvsat,qvsata,f34,tem1soil,          &
           tem2soil,tem3soil,tem4soil,tsdiffus,deltem,rrtem,temple)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Predict the soil surface temperature and moisture contents by solving
!  the surface energy and moisture budget equtions:
!
!  1. the soil temperatures,          tsoil(nx,ny,nz)
!  2. the soil moisture,              qsoil(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jerry Brotzge and Dan Weber
!  1/22/02
!
!  Updated code to allow for implicit soil temp/moisture scheme
!  and to allow for multiple options of soil/temp schemes.
!
!  MODIFICATION HISTORY:
!
!  12/15/02 Jerry Brotzge 
!           Added snow physics, based entirely on ETA code. 
!
!-----------------------------------------------------------------------
!
!  where
!
!      Tau    -- 1 day in seconds = 86400 seconds
!      PI     -- number of PI = 3.141592654
!      ROw    -- Density of liquid water
!      Veg    -- Vegetation fraction
!      Ct     -- Thermal capacity
!      Rn     -- Radiation flux, rnflx
!      H      -- Sensible heat flux, shflx
!      LE     -- Latent heat flux, lhflx = latent*(Eg + Ev)
!
!      Eg     -- Evaporation from ground
!      Ev     -- Evapotranspiration from vegetation, Ev = Etr + Er
!      Etr    -- Transpiration of the remaining part of the leaves
!      Er     -- Evaporation directly from the foliage covered by
!                intercepted water
!
!      P      -- Precipitation rates
!      Pg     -- Precipitation reaching the ground,
!                Pg = (1 - Veg) P
!
!      Wgeq   -- Surface volumetric moisture
!      C1, C2 -- Coefficients
!
!
!  For detailed information about the surface energy budget model,
!  see the articles in the reference list.
!
!-----------------------------------------------------------------------
!
!  REFERENCES:
!
!  Bougeault, P., et al., 1991: An Experiment with an Advanced
!       Surface Parameterization in a Mesobeta-Scale Model. Part I:
!       Implementation, Monthly Weather Review, 119, 2358-2373.
!
!  Chen, F., and J. Dudhia, 2001: Coupling an Advanced Land Surface-
!       Hydrology Model with the Penn State-NCAR MM5 Modeling System.
!       Part 1: Model Implementation and Sensitivity, Mon Wea Rev.,
!       129, 569-585.
!
!  Ek, M., and L. Mahrt, 1991: OSU 1-D PBL model user's guide.
!       Version 1.04, 120 pp.  [Available from the Dept of
!       Atmospheric Sciences, Oregon State Univ., Corvallis, OR
!       97331-2209].
!
!  Jacquemin, B. and J. Noilhan, 1990: Sensitivity Study and
!       Validation of a Land Surface Parameterization Using the
!       HAPEX-MOBILHP Data Set, Boundary-Layer Meteorology, 52,
!       93-134, (JN).
!
!  Noilhan, J. and S. Planton, 1989: A Simple Parameterization of
!       Land Surface Processes for Meteorological Model, Mon. Wea.
!       Rev., 117, 536-549, (NP).
!
!  Peters-Lidard, C.D., E. Blackburn, X. Liang, and E.F. Wood,
!       1998: The effect of soil thermal conductivity parameterization
!       on surface energy fluxes and temperatures, J. Atmos. Sci.,
!       55, 1209-1224.
!
!  Pleim, J. E. and A. Xiu, 1993: Development and Testing of a
!       Land-Surface and PBL Model with Explicit Soil Moisture
!       Parameterization, Preprints, Conf. Hydroclimat., AMS, 45-51,
!       (PX).
!
!  Smirnova, T. G., et al., 1997: Performance of Different Soil
!       Model Configurations in Simulating Ground Surface Temperatures
!       and Surface Fluxes, Mon. Wea. Rev., 125, 1870-1884.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (sfc/top)
!    nzsoil   Number of grid points in the z-direction (sfc/soil substrate)
!    zpsoil   The physical depth defined at w-point in soil
!    j3soil   Coordinate transformation Jacobian, d(zpsoil)/d(zsoil)
!    j3soilinv Inverse of j3soil
!
!    soiltyp  Soil type at the horizontal grid points
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    veg      Green vegetation fraction
!
!    windsp   Wind speed just above the surface (m/s)
!    psfc     Surface pressure (Pascal)
!    rhoa     Near sfc air density
!    precip   Precipitation flux reaching the surface [kg m-2 s-1]  
!
!
!    cdma     Array for cdm, surface drag coefficient for momentum
!    cdha     Array for cdh, surface drag coefficient for heat
!    cdqa     Array for cdq, surface drag coefficient for moisture
!
!    pres     3-dimensional pressure
!    temp     3-dimensional temperature
!    qv       3-dimensional specific humidity
!
!  INPUT and OUTPUT:
!
!    tsfc     Temperature at skin surface (K)
!    qsfc     Moisture at skin surface
!    tsoil    Soil temperatures (K)
!    qsoil    Soil moisture
!    wetcanp  Canopy wetness
!
!  OUTPUT:
!
!    qvsfc    Effective S. H. at sfc.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil                 ! Grid points in the soil
  INTEGER :: snowng(nx,ny)          ! Snow flag 
  INTEGER :: frzgra(nx,ny)          ! Freezing rain flag 

  REAL :: zpsoil (nx,ny,nzsoil)      ! The depth of the soil.
  REAL :: j3soil (nx,ny,nzsoil)      ! Coordinate transformation
  REAL :: j3soilinv (nx,ny,nzsoil)   ! Inverse of j3soil

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
  INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: veg    (nx,ny)        ! Green vegetation fraction

  REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface
  REAL :: tsoil  (nx,ny,nzsoil) ! Soil temperature at nzsoil levels
  REAL :: qsoil  (nx,ny,nzsoil) ! Soil moisture at nzsoil levels
  REAL :: xqsoil (nx,ny,nzsoil) ! Volume of nonfrozen soil moisture 
  REAL :: qice   (nx,ny,nzsoil) ! Ice content of each soil layer 
  REAL :: wetcanp(nx,ny)        ! Canopy wetness

  REAL :: snowdpth(nx,ny)       ! Snow depth (m)

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
  REAL :: cdma   (nx,ny)        ! Surface drag coeff. for momentum
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

  REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
  REAL :: radswnet(nx,ny)       ! Net solar radiation, SWin - SWout  
  REAL :: radlwin(nx,ny)        ! Incoming longwave radiation to sfc
  REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
  REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface [W/m2]
  REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
  REAL :: ct     (nx,ny)        ! Soil thermal coefficient

  REAL :: evaprg (nx,ny)        ! Evaporation
  REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
  REAL :: evapcan(nx,ny)        ! Direct evaporation from canopy leaves
  REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
  REAL :: qvsat  (nx,ny)        !

  REAL :: aatem                 ! Temporary array for est'g Pot Evapo.
  REAL :: fftem   (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: rrtem   (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: deltem  (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: rnettot               ! Temporary array for est'g Pot Evapo.
  REAL :: evappot (nx,ny)       ! Potential Evaporation, [kg/m2/s] 
  REAL :: temple  (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: rsttemp (nx,ny)       ! Canopy resistance used for LE,tr
!  REAL :: etp1    (nx,ny)       ! Potential evaporation * 0.001  


!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
!  Parameters and variables are defined in globcst.inc:
!
!    dtsfc      Surface model time step
!    nsfcst     # of surface model time steps
!
!    moist      Moist flag
!
!    year       Reference year
!    month      Reference month
!    day        Reference day
!    jday       Reference Julian day
!    hour       Hour of reference time
!    minute     Minute of reference time
!    second     Second of reference time
!
!    latitud    Latitude at the domain center
!    longitud   Longitude at the domain center
!
!    curtim     Current model time
!    dtbig      Length of big time step
!
!    bslope     Slope of the retention curve
!    cgsat      Soil thermal coefficient for bare ground at saturation
!    cgv        Soil thermal coef. for totally shielded ground by veg.
!    wsat       Saturated volumetric moisture content. JN, Tab. 1
!    wfc        Field capacity moisture. JN, Tab. 1
!    wwlt       Wilting volumetric moisture content. JN, Tab. 1
!
!  Parameters and variables are defined in phycst.inc:
!
!    solarc     Solar constant (W/m**2)
!    emissg     Emissivity of the ground
!    emissa     Emissivity of the atmosphere
!    sbcst      Stefen-Boltzmann constant
!
!    rhow       Liquid water reference density (kg/m**3)
!    rd         Gas constant for dry air (kg/(m s**2))
!    cp         Gas heat capacity at constant pressure
!    cv         Gas heat capacity at constant volume
!    lathv      LH of vaporization at 0 deg C (Lv = 2.50, [m^2/s^2])
!
!    tsoil0     Initial skin temperature
!    qsoil0     Initial skin wetness
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!

  REAL :: pi                         ! Pi
  PARAMETER (pi = 3.141592654)

  REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun

  REAL :: wrmax       ! Maximum value for canopy moisture, wetcanp
  REAL :: wr2max      ! Tendency to reach the maximum wrmax
  REAL :: runoff      ! Runoff of the interception reservoir.
  REAL :: pnet        ! Residual of precip. and evap.
  REAL :: vegp        ! Precip. intercepted by vegetation

  INTEGER :: i, j, k, l 
  INTEGER :: temcount 

  REAL :: temsum 
  REAL :: pf                    ! log(psi)
  REAL :: cisoil(nx,ny,nzsoil)  ! Volumetric heat capacity
  REAL :: ktsoiltop(nx,ny)      ! Thermal conductivity, top soil level
  REAL :: totaldepth            ! Thickness of total soil column 

  REAL :: tsdiffus(nx,ny,nzsoil)    ! Thermal diffusivity

  REAL :: tem1soil (nx,ny,nzsoil)   ! Tridiagonal parameters
  REAL :: tem2soil (nx,ny,nzsoil)   ! Tridiagonal parameters
  REAL :: tem3soil (nx,ny,nzsoil)   ! Tridiagonal parameters
  REAL :: tem4soil (nx,ny,nzsoil)   ! Tridiagonal parameters
  REAL :: plantcoef(nx,ny) 
  REAL :: sink     (nx,ny,nzsoil)   ! Sink/source term

  REAL :: tempcoefa        ! Temporary arrays for tridiag.
  REAL :: tempcoefb        ! Temporary arrays for tridiag.
  REAL :: tempcoefc        ! Temporary arrays for tridiag.
  REAL :: tempcoefx1       ! Temporary arrays for tridiag.
  REAL :: tempcoefx2       ! Temporary arrays for tridiag.

  REAL :: tempbeta         ! Used to estimate skin temperature
  REAL :: potair  (nx,ny)  ! Potential air temperature
!  REAL :: tempcg           ! Used to compute cg

  REAL :: tema, temb, temc
  REAL :: rhswr
  REAL :: desdt
  REAL :: dqsdt
  REAL :: dew              ! Dew, forms when pot E < 0
  REAL :: sn_new           ! New snow depth  
  REAL :: sndens (nx,ny)   ! Snow density 
  REAL :: sncond (nx,ny)   ! Snow thermal conductivity 
  REAL :: sneqv  (nx,ny)   ! Snow water equivalent (m) 
  REAL :: precip1(nx,ny)   ! Precip including snow melt  
  REAL :: precipdrip(nx,ny) ! Precip that infiltrates ground
  REAL :: snofac (nx,ny)   ! Snow factor 
  REAL :: flx2   (nx,ny)   ! Freezing rain flux from LH release 
  REAL :: beta   (nx,ny)   ! Beta for snow physics 
  REAL :: newsnow
  REAL :: newdens  

  REAL :: denom, t12a, t12b, t12, snomeltrt, snomeltamt
  REAL :: temd, sice  

  REAL :: rsnow, albedo, qtotal
  

  LOGICAL :: firstcall     ! First call flag of this subroutine

  SAVE firstcall, log100 
  DATA firstcall/.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (firstcall) THEN
    log100    = ALOG(100.0)
    firstcall = .false.
  END IF

!
!-----------------------------------------------------------------------
!
!  Calculate saturated specific humidity near the surface, qvsata.
!
!-----------------------------------------------------------------------
!

  IF ( moist /= 0 ) THEN
    CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tair,qvsata)
  END IF


!-----------------------------------------------------------------------
!  Converts standard soil temperatures to potential temperatures.  
!  Estimates nonfrozen soil moisture content (vol fraction), xqsoil. 
!-----------------------------------------------------------------------

  DO k=1,nzsoil 
    DO j=1,ny-1
      DO i=1,nx-1

        IF (soiltyp(i,j) /= 12  .AND. soiltyp(i,j) /= 13 ) THEN 
          tsoil(i,j,k) = tsoil(i,j,k) * (100000.0/psfc(i,j))**0.286 
        END IF

        IF (tsoil(i,j,k) < 273.15) THEN
          xqsoil(i,j,k) = 0.0   
        ELSE
          xqsoil(i,j,k) = qsoil(i,j,k) 
        END IF

!  NOTE: Should xqsoil be set to wwlt(soiltyp(i,j)) if frozen???*****

      END DO
    END DO 
  END DO 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                Beginning of Implicit Scheme (JAB/DBW)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( moist /= 0 ) THEN 

      CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsoil(1,1,1),qvsat)

  END IF 


!---------------------------------------------------------------------
!    Calculate snow density and conductivity 
!---------------------------------------------------------------------

  DO j=1,ny-1
    DO i=1,nx-1 

!      veg(i,j) = 0.41    !Test run 

      IF (snowdpth(i,j) == 0.0) THEN 
        sndens(i,j) = 0.0       ! Snow density (unitless)  
        sneqv (i,j) = 0.0       ! Snow water equivalent 
        sncond(i,j) = 1.0       ! Snow thermal conductivity  

      ELSE 
!       sndens(i,j) = sneqv(i,j)/snowdpth(i,j)  
        sndens(i,j) = 0.10       ! Set arbitrarily until sneqv is available 
        sneqv(i,j) = sndens(i,j)*snowdpth(i,j) 

        !Dyachkova Eqtn (1960) Units: Cal/(cm*hr*C) converted to W/(m*C) 
        sncond(i,j) = 0.11631 * (0.328*10**(2.25*sndens(i,j))) 
      END IF 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    Determine if snowing or freezing rain
!----------------------------------------------------------------------

      snowng(i,j) = 0                  ! Snowing flag 
      frzgra(i,j) = 0                  ! Freezing rain flag  

      IF (precip(i,j) > 0.0) THEN 
        IF (tair(i,j) <= 273.15) THEN 
          snowng(i,j) = 1  
        ELSE
          IF (tsoil(i,j,1) <= 273.15) frzgra(i,j) = 1 
        END IF 
      END IF

!----------------------------------------------------------------------
!     IF precipation is frozen, then update new snowfall depth 
!-----------------------------------------------------------------------

      IF ( (snowng(i,j) == 1) .OR. (frzgra(i,j) == 1) ) THEN 
!        sn_new = precip(i,j) * 1000.0 * dtsfc * 0.001
        sn_new = precip(i,j) * dtsfc * 0.001


        sneqv(i,j) = sneqv(i,j) + sn_new 
        precip1(i,j) = 0.0 

!-----------------------------------------------------------------------
!       Calculate new snowfall density 
!-----------------------------------------------------------------------
        IF (tair(i,j) <= 258.15) THEN
          newdens = 0.05
        ELSE
          newdens = 0.05 + 0.0017*(tair(i,j)-273.15+15.0)**1.5 
        END IF 

!       Adjust snow density based on new snowfall
        newsnow = (sn_new*100.0/newdens)  
        sndens(i,j) = (snowdpth(i,j)*100.0*sndens(i,j)+newsnow*newdens) &
            / (snowdpth(i,j)*100.0 + newsnow)  
        snowdpth(i,j) = 0.01 * (snowdpth(i,j)*100.0 + newsnow)
!----------------------------------------------------------------------------
 
        !Dyachkova Eqtn (1960) Units: Cal/(cm*hr*C) converted to W/(m*C)
        ! Snow conductivity 
        sncond(i,j) = 0.11631 * (0.328*10**(2.25*sndens(i,j))) 

      ELSE
!!        precip1(i,j) = precip(i,j)*1000.0    ![kg m-2 s-1] 
!        precip1(i,j) = precip(i,j)           ![m s-1] 

        precip1(i,j) = precip(i,j)/1000.0      ![m s-1] 


      END IF 

!-----------------------------------------------------------------------
!     Amend albedo as a function of snow cover 
!-----------------------------------------------------------------------

      IF (radsw(i,j) > 0.0) THEN 
        albedo = (radswnet(i,j)-radsw(i,j)) / radsw(i,j)  
        IF (albedo > 1.0) albedo = 1.0 
        IF (albedo < 0.0) albedo = 0.0 
      END IF

        IF (sneqv(i,j) == 0.0) THEN 
!          albedo = alb    !**NOTE: No access to ETA ALB, SNOALB tables***  
        ELSE 
          IF (sneqv(i,j) < snup(vegtyp(i,j))) THEN 
             rsnow = sneqv(i,j) / snup(vegtyp(i,j)) 
             snofac(i,j) = 1.0 - (EXP(-2.6*rsnow) - rsnow*EXP(-2.6))
          ELSE
             snofac(i,j) = 1.0 
          END IF 
         
!           albedo = alb + (1.0-veg(i,j))*snofac(i,j)*(snoalb-alb)
!           IF (albedo > snoalb) albedo = snoalb 

        END IF 

    END DO     ! i index 
  END DO       ! j index  

!-----------------------------------------------------------------------
!    Initialize heat capacity (C) and thermal conductivity (K) for each
!          level
!-----------------------------------------------------------------------

!  Compute soil thermal conductivity and heat capacity at each level k

    CALL getcon(nx, ny, nzsoil, soiltyp, vegtyp, tsoil(1,1,1),       &  
            qsoil(1,1,1), ktsoiltop, cisoil, tsdiffus, xqsoil(1,1,1))   

!------------------------------------------------------------------------
!  Add subsurface heat flux reduction effect from the overlying green
!    canopy.  See Peters-Lidard et al., 1997, JGR, Vol 102(D4)
!------------------------------------------------------------------------

    DO j=1,ny-1
      DO i=1,nx-1
!         ktsoiltop(i,j) = ktsoiltop(i,j) * EXP(-2.0*veg(i,j))

!         ktsoiltop(i,j) = ktsoiltop(i,j) * EXP(-0.5*lai(i,j))

      END DO
    END DO

!-----------------------------------------------------------------------
!    Compute G and the potential evaporation via Penman eqtn.
!-----------------------------------------------------------------------

    CALL penman(nx, ny, nzsoil, soiltyp, tsoil(1,1,1), qsoil(1,1,1), &
            ktsoiltop, fftem, j3soilinv, deltem, rrtem, veg,         &
            temple, evappot, qvair, qvsata, tair, cdha, psfc,       &
            potair, precip, rhoa, windsp, gflx, radswnet, radlwin,  &
            sneqv, snowdpth, snofac, sncond, snowng, frzgra, flx2)


!-------------------------------------------------------------
!    Compute the canopy resistance  
!-------------------------------------------------------------

      CALL canres(nx,ny, radsw, f34, cdqa,                          &
                   windsp,psfc,rhoa,qvair,                          &
                   soiltyp,vegtyp,lai,veg,tair,tsoil(1,1,1),        &
                   qsoil,wetcanp,nzsoil,zpsoil,snowdpth,            &
                   qvsata,rsttemp,plantcoef,cdha)


!*********************************************************************
!    BEGIN NOPAC/SNOPAC CODE ---------------------------------------
!*********************************************************************

!---------------------------------------------------------------
!     Add Snow Physics
!---------------------------------------------------------------

      CALL snophy(nx,ny,nzsoil,tsoil,tair,potair,precip,evappot,     &
                 snowdpth,sneqv,sndens,snowng,snofac,fftem,rrtem,    &
                 ktsoiltop,gflx,psfc,precip1,qvair,cdha,rhoa,flx2,   &
                 beta)

!-------------------------------------------------------------
!    Compute the initial LE flux estimate
!-------------------------------------------------------------


      CALL leflx(nx,ny,nzsoil,cdha,cdqa,windsp, rhoa, qvair, veg, wetcanp,   &
           deltem,rrtem,temple, lhflx, evappot,soiltyp,vegtyp,evaprg,        &
           evaprtr,evapcan,qvsfc,qvsat,rsttemp,plantcoef,zpsoil,         &
           qsoil, xqsoil, qice, precip1, precipdrip)

!--------------------------------------------------------------


!-------------------------------------------------------------
!    Computes skin temperature.  
!-------------------------------------------------------------
!
!      CALL tskin(nx,ny, nzsoil, cdha, rhoa, tair, potair,     &
!               rrtem, fftem, ktsoiltop, psfc, evappot, lhflx, tsoil,  &
!               soiltyp,snowdpth)
!
!
!---------------------------------------------------------------------
!  Calculate the canopy wetness, wr
!----------------------------------------------------------------------

  IF ( moist /= 0) THEN

    CALL canwet(nx, ny, soiltyp, veg, lai, snowdpth, evapcan,         &
              precip, wetcanp)

  END IF

!
!----------------------------------------------------------------------
!   Estimate Soil  Moisture
!-----------------------------------------------------------------------
!

  CALL qdiff(nx,ny,nzsoil, j3soilinv, soiltyp, lai, veg,                 &
           tsoil, qsoil, wetcanp, windsp, rhoa, precip,                  &
           qvair, qvsata, evapcan, evaprg, tem1soil, tem2soil, tem3soil, &
           tem4soil, precipdrip)

!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver to solve the diffusion of soil moisture
!  from k= 2, nzsoil-1.
!
!-----------------------------------------------------------------------
!

  CALL tridiag2(nx,ny,nzsoil,1,nx-1,1,ny-1,2,nzsoil-1,tem1soil,tem2soil, &
           tem3soil,tem4soil)

  DO k=2,nzsoil-1     ! load the final result from tridiag2 into qsoil
    DO j=1,ny-1
      DO i=1,nx-1
        IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN 

          qsoil(i,j,k) = tem4soil(i,j,k)

        END IF 
      END DO
    END DO
  END DO              ! qsoil is updated.........


!--------------------------------------------------------------------
!    Compute bottom boundary condition
!---------------------------------------------------------------------

  DO j=1,ny-1
    DO i=1,nx-1
      IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN 

        qsoil(i,j,nzsoil) = qsoil(i,j,nzsoil-1)

      END IF  ! Soiltyp /= 12 or 13 (Ice and water) 
    END DO
  END DO


!-----------------------------------------------------------------------
!    Initialize heat capacity (C) and thermal conductivity (K) for each
!          level
!-----------------------------------------------------------------------

    CALL getcon(nx, ny, nzsoil, soiltyp, vegtyp, tsoil(1,1,1),       &
            qsoil(1,1,1), ktsoiltop, cisoil, tsdiffus, xqsoil(1,1,1))

!------------------------------------------------------------------------
!  Add subsurface heat flux reduction effect from the overlying green
!    canopy.  See Peters-Lidard et al., 1997, JGR, Vol 102(D4)
!------------------------------------------------------------------------

    DO j=1,ny-1
      DO i=1,nx-1
!         ktsoiltop(i,j) = ktsoiltop(i,j) * EXP(-2.0*veg(i,j))

!         ktsoiltop(i,j) = ktsoiltop(i,j) * EXP(-0.5*lai(i,j))

      END DO
    END DO


!--------------------------------------------------------------------
!    Compute skin temperature.
!--------------------------------------------------------------------

      CALL tskin(nx,ny, nzsoil, cdha, rhoa, tair, potair,     &
               rrtem, fftem, ktsoiltop, psfc, evappot, lhflx, tsoil,  &
               soiltyp,snowdpth,veg)


!------------------------------------------------------------------------
!   Estimate Soil Temperature
!-----------------------------------------------------------------------
!
!  C * dT/dt = d/dz (K * dT/dz)
!
!-----------------------------------------------------------------------
!
!  Calculate coefficients of the tridigonal equation for soil temperature
!
!-----------------------------------------------------------------------
!

!  DO k=2,nzsoil-1
!    DO j=1,ny-1
!      DO i=1,nx-1
!
!      tempcoefa = cnbeta * dtsfc * 0.5 * dzsoilinv2 * j3soilinv(i,j,k)
!      tempcoefb = (1.0 - cnbeta) * dtsfc * 0.5 * dzsoilinv2 * j3soilinv(i,j,k)
!
!      tempcoefx1 = (tsdiffus(i,j,k-1) * j3soilinv(i,j,k-1)) +        &
!                   (tsdiffus(i,j,k  ) * j3soilinv(i,j,k  ))
!
!      tempcoefx2 = (tsdiffus(i,j,k  ) * j3soilinv(i,j,k  )) +        &
!                   (tsdiffus(i,j,k+1) * j3soilinv(i,j,k+1))
!
!        tem1soil(i,j,k) = -tempcoefa * tempcoefx1
!        tem2soil(i,j,k) = 1.0 + tempcoefa * (tempcoefx1 + tempcoefx2)
!        tem3soil(i,j,k) = -tempcoefa * tempcoefx2
!
!        tem4soil(i,j,k) = tempcoefb * tempcoefx1 * tsoil(i,j,k-1) +  &
!        (1.0-tempcoefb*(tempcoefx1 + tempcoefx2))*tsoil(i,j,k)       &
!          + tempcoefb * tempcoefx2 *tsoil(i,j,k+1)
!
!
!        denom = (zpsoil(i,j,k-1)-zpsoil(i,j,k)) * cisoil(i,j,k) 
!        qtotal = -1.0 * denom * tem4soil(i,j,k) 
! 
!        sice = qsoil(i,j,k) - xqsoil(i,j,k) 
!
!        IF ( (sice > 0.0).OR.(tsoil(i,j,1) < 273.15) .OR.              &
!             (tsoil(i,j,2)<273.15).OR.(tsoil(i,j,nzsoil)<273.15) ) THEN
!
!---------------------------------------------------------------------
!    Compute sink/source term for freezing and thawing
!---------------------------------------------------------------------
!
          CALL snksrc(nx,ny,nzsoil, tsoil, qsoil, xqsoil, ktsoiltop,   &
                gflx, zpsoil, sink, soiltyp, tem1soil, tem2soil,       &
                tem3soil, tem4soil, tsdiffus, j3soilinv, cisoil) 


!          tem4soil(i,j,k) = (tempcoefb * tempcoefx1 * tsoil(i,j,k-1) + &
!            (1.0-tempcoefb*(tempcoefx1 + tempcoefx2))*tsoil(i,j,k)     &
!            + tempcoefb * tempcoefx2 *tsoil(i,j,k+1) ) -               & 
!            (sink(i,j,k)*dzsoilinv/cisoil(i,j,k)) 
!
!----------------------------------------------------------------------
!        END IF
!
!
!      END DO
!    END DO
!  END DO

! Set the lower temperature boundary conditions
! note lower boundary is zero gradient
! (tsoil(i,j,nzsoil) =  tsoil(i,j,nz-1)

  DO j=1,ny-1     !  must set the boundry condition for tem2soil(i,j,nzsoil-1)
    DO i=1,nx-1
      tem4soil(i,j,2) = tem4soil(i,j,2) - tem1soil(i,j,2)*tsoil(i,j,1) !Top BC
      tem2soil(i,j,nzsoil-1) = tem2soil(i,j,nzsoil-1) + tem3soil(i,j,nzsoil-1)
               !Bottom, zero gradient
    END DO
  END DO          !  coefficients are ready for the crank-nicholson
                  !  solver.

!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver.
!
!-----------------------------------------------------------------------
!

  CALL tridiag2(nx,ny,nzsoil,1,nx-1,1,ny-1,2,nzsoil-1,tem1soil,   &
       tem2soil,tem3soil,tem4soil)

  DO k=2,nzsoil-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN 

          tsoil(i,j,k) = tem4soil(i,j,k)

        END IF 
      END DO
    END DO
  END DO            ! soil temperature is up to date....

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                    End of qsoil and tsoil implicit solving technique
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-------------------------------------------------------------
!    Compute skin temperature.
!-------------------------------------------------------------

      CALL tskin(nx,ny, nzsoil, cdha, rhoa, tair, potair,     &
               rrtem, fftem, ktsoiltop, psfc, evappot, lhflx, tsoil,  &
               soiltyp,snowdpth,veg)


!---------------------------------------------------------------------
!
!  Calculate the canopy wetness, wr
!
!----------------------------------------------------------------------

  IF ( moist /= 0) THEN

   CALL canwet(nx, ny, soiltyp, veg, lai, snowdpth, evapcan,       &
              precip, wetcanp )

  END IF



!----------------------------------------------------------------------
!     Snow compaction - new estimates of snow depth and density 
!----------------------------------------------------------------------

   CALL snocom(nx, ny, nzsoil, snowdpth, sneqv, sndens, sncond, tsoil)


!-----------------------------------------------------------------------
!
!  Calculate the heat capacity cg and ct, sensible heat flux shflx,
!  and ground heat flux gflx
!
!-----------------------------------------------------------------------
!
    DO j = 1, ny-1
      DO i = 1, nx-1

!        tempcg = MAX( qsoil(i,j,1), wwlt(soiltyp(i,j)) )            !(JAB)

        IF (snowdpth(i,j) >= snowdepth_crit) THEN !snow cover
!          ct(i,j) = cg_snow
!
!          gflx(i,j) = 2.0*pi*(tsoil(i,j,1)-tsoil(i,j,2))                   &
!                            *snowflxfac/(tau*ct(i,j))
                                        ! Snow cover
        ELSE

!-------------------------------------------
!  Estimate 2nd G flux for plotting  
!-------------------------------------------

       totaldepth = zrefsfc - zpsoil(i,j,nzsoil) 

       IF (totaldepth > 0.10) THEN 

         temsum = 0.0 
         temcount = 0 
         
         DO k=2,nzsoil
           zrefsoil = zrefsfc - zpsoil(i,j,k) 
           IF (zrefsoil <= 0.10) THEN 
             temsum = temsum + (tsoil(i,j,k-1) - tsoil(i,j,k))   
             temcount = temcount + 1 
           END IF
         END DO 

         IF (temcount > 0) THEN         
!            gflx(i,j) = ktsoiltop(i,j)*j3soilinv(i,j,2)*dzsoilinv*  &
!               (1.0/temcount) * temsum 

            gflx(i,j) = ktsoiltop(i,j)* EXP(-2.0*veg(i,j)) *    & 
               j3soilinv(i,j,2)*dzsoilinv*                      &
               (1.0/temcount) * temsum 

         ELSE IF (temcount == 0) THEN 

            gflx(i,j) = ktsoiltop(i,j)* EXP(-2.0*veg(i,j)) *     & 
               j3soilinv(i,j,2)*dzsoilinv*                       &
               0.5*(tsoil(i,j,1)-tsoil(i,j,2))

!            gflx(i,j) = ktsoiltop(i,j)*j3soilinv(i,j,2)*dzsoilinv*  &
!               0.5*(tsoil(i,j,1)-tsoil(i,j,2))


         END IF 

        ELSE IF (totaldepth <= 0.10) THEN 
!            gflx(i,j) = ktsoiltop(i,j)*j3soilinv(i,j,2)*dzsoilinv*  &
!               0.5*(tsoil(i,j,1)-tsoil(i,j,2))  


            gflx(i,j) = ktsoiltop(i,j)* EXP(-2.0*veg(i,j)) *     & 
               j3soilinv(i,j,2)*dzsoilinv*                       &
               0.5*(tsoil(i,j,1)-tsoil(i,j,2))  


        END IF   
 
        END IF

        shflx(i,j) = rhoa(i,j)*cp*cdha(i,j)*                &
                    ( tsoil(i,j,1) - potair(i,j) )  ! Sensible heat flux
                                                ! NP, Eq. 26

        IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN
          tsfc(i,j) = tsoil(i,j,1)
        END IF 

      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the water vapor flux, qvflx, and hence the latent heat
!  flux, lhflx. They share the same array lhflx(i,j).
!
!  If moisture flag is off, all moisture fields are set to zero.
!
!-----------------------------------------------------------------------
!
!  Calculate saturated specific humidity at the ground surface.
!
!-----------------------------------------------------------------------

  IF (moist /= 0) THEN 

      CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)

  END IF


    DO j = 1, ny-1
      DO i = 1, nx-1

        cdha(i,j) = cdha(i,j) / MAX(windsp(i,j),0.1)

        IF (snowdpth(i,j) > 0.0) THEN
          lhflx(i,j) = evappot(i,j) * beta(i,j) * lathv 
        ELSE 
!          lhflx(i,j) = evappot(i,j) * lathv 
        END IF 

!        IF (lhflx(i,j) > (evappot(i,j)*lathv)) lhflx(i,j)=evappot(i,j)*lathv

        qvsfc(i,j) = lhflx(i,j)/(lathv*rhoa(i,j)*cdqa(i,j)*windsp(i,j)) &
                + qvair(i,j)

        IF (qvsfc(i,j) > qvsat(i,j)) qvsfc(i,j) = qvsat(i,j)


!-----------------------------------------------------------------------
!
!  Converts potential soil temperatures to standard soil temperatures.  
!
!-----------------------------------------------------------------------

  IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN 

        DO k=1,nzsoil 
          tsoil(i,j,k) = tsoil(i,j,k) * ((100000.0 / psfc(i,j)) **(-0.286))
          tsfc(i,j) = tsoil(i,j,1) 
        END DO 

  END IF 

      END DO
    END DO


  RETURN
END SUBROUTINE ousoil


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CANWET                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE canwet(nx, ny, soiltyp, veg, lai, snowdpth, evapcan,   &
              precip, wetcanp)

!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate the canopy wetness.
!            Uses a forward scheme.
!
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  7/23/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    soiltyp  Soil type
!    veg      Vegetation fraction
!    snowdpth Snow depth (m)
!    evapcan  Evaporation from canopy
!    precip   Precipitation rate
!
!  OUTPUT:
!
!    wetcanp  Soil water content on canopy
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny
  INTEGER :: soiltyp(nx,ny)

  REAL :: precip (nx,ny)        ! Precipitation rate at the surface
  REAL :: veg    (nx,ny)        ! Vegetation fraction
  REAL :: lai    (nx,ny)        ! Leaf area index
  REAL :: snowdpth(nx,ny)       ! Snow depth
  REAL :: evapcan(nx,ny)        ! Evaporation from canopy
  REAL :: wetcanp(nx,ny)        ! Soil water content of canopy

!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j

  REAL :: wrmax
  REAL :: wr2max
  REAL :: rhswr
  REAL :: vegp
  REAL :: pnet
  REAL :: tema
  REAL :: runoff

  REAL :: temwra(nx,ny)   ! Temporary array
  REAL :: temwrb(nx,ny)   ! Temporary array


!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

    DO j = 1, ny-1
      DO i = 1, nx-1

         wrmax = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)      ! meter

         IF ( soiltyp(i,j) == 12 .OR.   soiltyp(i,j) == 13 .OR.         &
             snowdpth(i,j) >= snowdepth_crit ) THEN

           temwra(i,j) = 0.0

         ELSE

           wr2max = ( wrmax - wetcanp(i,j) ) / dtsfc
           vegp = veg(i,j) * precip(i,j)
           pnet = vegp - evapcan(i,j)
           tema = pnet - wr2max * rhow
           runoff = MAX( tema, 0.0 )
           vegp = vegp - runoff
           temwra(i,j) = ( vegp - evapcan(i,j) ) / rhow

         END IF

          temwrb(i,j) = wetcanp(i,j) + dtsfc * temwra(i,j)
          temwrb(i,j) = MAX( temwrb(i,j), 0.0 )
          temwrb(i,j) = MIN( temwrb(i,j), wrmax )


         IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.      &
             snowdpth(i,j) >= snowdepth_crit ) THEN
           rhswr = 0.0
         ELSE
           wr2max = ( wrmax - temwrb(i,j) ) / dtsfc
           vegp = veg(i,j) * precip(i,j)
           pnet = vegp - evapcan(i,j)
           tema = pnet - wr2max * rhow
           runoff = MAX( tema, 0.0 )
           vegp = vegp - runoff
           rhswr = 0.5 * (temwra(i,j)+(vegp - evapcan(i,j) ) / rhow)
         END IF

           wetcanp(i,j) = wetcanp(i,j) + dtsfc * rhswr
           wetcanp(i,j) = MAX( wetcanp(i,j), 0.0 )
           wetcanp(i,j) = MIN( wetcanp(i,j), wrmax )

      END DO
    END DO



!    DO j = 1, ny-1
!      DO i = 1, nx-1
!
!        wrmax = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! meter
!        wr2max = (wrmax - wetcanp(i,j) )/ dtsfc
!
!      IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.         &
!             snowdpth(i,j) >= snowdepth_crit ) THEN
!          rhswr = 0.0
!      ELSE
!
!        vegp = veg(i,j) * precip(i,j)
!        pnet = vegp - evapcan(i,j)
!        tema = pnet - wr2max * rhow
!        runoff = MAX( tema, 0.0 )
!        vegp = vegp - runoff
!        rhswr = (vegp - evapcan(i,j) ) / rhow
!      END IF
!
!          wetcanp(i,j) = wetcanp(i,j) + dtsfc * rhswr
!          wetcanp(i,j) = MAX( wetcanp(i,j), 0.0 )
!          wetcanp(i,j) = MIN( wetcanp(i,j), wrmax )
!
!      END DO
!    END DO


  RETURN
END SUBROUTINE canwet


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE PENMAN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################

!##################################################################
!
!

SUBROUTINE penman(nx, ny, nzsoil, soiltyp, tsoil, qsoil,             &
            ktsoiltop, fftem, j3soilinv, deltem, rrtem, veg,         &
            temple, evappot, qvair, qvsata, tair, cdha, psfc,  &
            potair, precip, rhoa, windsp, gflx, radswnet, radlwin,   &
            sneqv, snowdpth, snofac, sncond, snowng, frzgra, flx2)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate the potential evaporation.
!
!
!      Ep = [ (Rn/rho/Cp/Ch) + (Ta-Ts)] * delta + (r + 1)*A
!           ----------------------------------------------- * (rho*Cp*Ch/Lv)
!                      delta + r + 1.0
!
!            delta = (qsfcsat - qairsat) / (Ts - Tair) * Lv / Cp
!
!
!             r = 4.0 * sigma * Tair**4.0 * Rd
!                 -----------------------------
!                         Psfc * Cp * Ch
!
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  7/18/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the z-direction (soil)
!    soiltyp  Soil type
!    tsoil    Soil temperatures (K)
!    qsoil    Soil moisture
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nzsoil

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point

  INTEGER :: snowng (nx,ny)     ! Snowing flag
  INTEGER :: frzgra (nx,ny)     ! Freezing rain flag 

  REAL :: ktsoiltop(nx,ny)        ! Thermal conductivity, top soil level
  REAL :: j3soilinv (nx,ny,nzsoil)   ! Inverse of j3soil

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: potair(nx,ny)         ! Potential temp (K)
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  REAL :: veg    (nx,ny)        ! Vegetation fraction
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface
  REAL :: tsoil  (nx,ny,nzsoil) ! Soil temperature at nzsoil levels
  REAL :: qsoil  (nx,ny,nzsoil) ! Soil moisture at nzsoil levels

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: gflx   (nx,ny)        ! Ground heat flux
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)

  REAL :: fftem   (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: deltem  (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: evappot (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: temple  (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: rrtem   (nx,ny)       ! Temporary array for est'g Pot Evapo.
  REAL :: flx2    (nx,ny)       ! Temp array for freezing rain 

  REAL :: radswnet(nx,ny)       ! Net shortwave radiation
  REAL :: radlwin (nx,ny)       ! Incoming longwave radiation

  REAL :: sneqv   (nx,ny)       ! Snow water equivalent
  REAL :: snowdpth(nx,ny)       ! Snow depth (m) 
  REAL :: sncond  (nx,ny)       ! Snow conductivity 
  REAL :: snofac  (nx,ny)       ! Snow factor  

!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j, k

  REAL :: rnettot               ! Temporary array for est'g Pot Evapo.
  REAL :: aatem                 ! Temporary array for est'g Pot Evapo.
  REAL :: tempbeta              ! Beta coefficient, soil mstr parameter
  REAL :: dew                   ! Dew

  REAL :: expsno, expsoi
  REAL :: depthtotal 
  REAL :: df1p, rch  

  REAL :: temb, temc, desdt, dqsdt

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  temb = 4.0 * sbcst * rd/cp
  temc = lathv * cpinv

  DO j=1,ny-1
    DO i=1,nx-1

     rch = rhoa(i,j) * cp * cdha(i,j) 

     cdha(i,j) = cdha(i,j) * MAX(windsp(i,j),0.1)

!-------------------------------------------
!  Estimate GH flux term
!-------------------------------------------

     IF (sneqv(i,j) == 0.0) THEN               !No snow cover 

        gflx(i,j)=ktsoiltop(i,j)* EXP(-2.0*veg(i,j)) *     & 
           dzsoilinv*j3soilinv(i,j,2)*                     &
          (tsoil(i,j,1)-tsoil(i,j,2))


!        gflx(i,j)=ktsoiltop(i,j)*dzsoilinv*j3soilinv(i,j,2)*  &
!          (tsoil(i,j,1)-tsoil(i,j,2))


     ELSE                                      !Snow cover (plane parallel) 
    
        depthtotal = snowdpth(i,j) + dzsoil 
        expsno = snowdpth(i,j) / depthtotal
        expsoi = dzsoil / depthtotal 

        df1p = expsno * sncond(i,j) + expsoi * ktsoiltop(i,j)
        ktsoiltop(i,j) = df1p * snofac(i,j) + ktsoiltop(i,j)*(1.0-snofac(i,j))

        gflx(i,j) = ktsoiltop(i,j) * j3soilinv(i,j,2) *        &
          (tsoil(i,j,1) - tsoil(i,j,2)) / depthtotal

      END IF   


!----------------------------------------------------------
!   Compute radiation and potential air temperature
!----------------------------------------------------------

    fftem(i,j) = radswnet(i,j) + radlwin(i,j)

    rnettot = fftem(i,j) - (sbcst * tair(i,j)**4.0) - gflx(i,j)

    aatem = temc * (qvsata(i,j) - qvair(i,j))

    potair(i,j) = tair(i,j) * ((100000.0 / psfc(i,j)) ** 0.286)


!---------------------------------------------------------------------
!   Include the latent heat effects of freezing rain converting to ice
!---------------------------------------------------------------------

    flx2(i,j) = 0.0 
    IF (frzgra(i,j) /= 0) THEN 
!     flx2(i,j) = -3.335E5 * precip(i,j) * 1000.0 
      flx2(i,j) = -3.335E5 * precip(i,j) 
      rnettot=rnettot - flx2(i,j)    
    END IF 

!---------------------------------------------------------------------
!   Adjust for latent heat effects caused by falling precipitation
!---------------------------------------------------------------------

    rrtem(i,j) = temb * (tair(i,j)**4.0) / (psfc(i,j)* cdha(i,j))

    IF (snowng(i,j) /= 1) THEN

!      IF (precip(i,j) > 0.0) rrtem(i,j) = rrtem(i,j) +           &  
!          4.218E3 * (precip(i,j)*1000.0) / rch  

      IF (precip(i,j) > 0.0) rrtem(i,j) = rrtem(i,j) +           &  
          4.218E3 * precip(i,j) / rch  

    ELSE 

!      rrtem(i,j) = rrtem(i,j)+2.106E3*(precip(i,j)*1000.0) / rch  

      rrtem(i,j) = rrtem(i,j)+2.106E3*precip(i,j) / rch  


    END IF 


!---------------------------------------------------------
!  Compute Clausius-Clapyron Eqtn.
!---------------------------------------------------------
   dqsdt = 0.0
   desdt = 0.0

   Do k = 7,2,-1
      desdt = alpha(k)*(k-1) + (tair(i,j)-273.16)*desdt
   END DO

    dqsdt = 0.622 * desdt / psfc(i,j)

    deltem(i,j) = temc * dqsdt

!    evappot(i,j)=( ((rnettot/(rhoa(i,j)*cp*cdha(i,j)))+(potair(i,j)-tair(i,j))) &
!      * deltem(i,j) + ((rrtem(i,j) + 1.0) * aatem) + (aatem - deltem(i,j)*   &
!        tair(i,j)) * (((100000.0/psfc(i,j))**0.286)-1.0) )    &
!        / ( (deltem(i,j) + rrtem(i,j) + 1.0)                  &
!        + (((100000.0/psfc(i,j))**0.286)-1.0) )               &
!        * (rhoa(i,j) * cp *cdha(i,j)/lathv)

    evappot(i,j)=( ((rnettot/rch)+(potair(i,j)-tair(i,j)))     &
        * deltem(i,j) + ((rrtem(i,j) + 1.0) * aatem) )         &
        / ( deltem(i,j) + rrtem(i,j) + 1.0)                    &
        * (rch/lathv)
!  (This 2nd version of evappot is smoother temporally...) 



    IF (soiltyp(i,j) /= 13) THEN
      tempbeta = (qsoil(i,j,1) - wwlt(soiltyp(i,j))) /          &
            (wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )
    ELSE IF (soiltyp(i,j) == 13) THEN
      tempbeta = 1.0
    END IF

!---------------Set upper limit on soil wetness
    IF (tempbeta > 1.0) tempbeta = 1.0

!---------------Set lower limit on soil wetness
    IF (tempbeta < 0.0) tempbeta = 0.0

      temple(i,j) = (1.0 - veg(i,j)) * tempbeta * evappot(i,j)

      END DO
   END DO


  RETURN
END SUBROUTINE penman


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GETCON                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getcon(nx,ny,nzsoil, soiltyp, vegtyp, tsoil, qsoil,       &
           ktsoiltop, cisoil, tsdiffus, xqsoil)        

!
!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate the thermal conductivity and diffusivity  

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  7/18/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the z-direction (soil)
!    soiltyp  Soil type 
!    vegtyp   Vegetation type 
!    tsoil    Soil temperatures (K) 
!    qsoil    Soil moisture  
!
!  OUTPUT:
!
!     cisoil        ! Volumetric heat capacity
!     ktsoiltop     ! Thermal conductivity, top soil level
!     tsdiffus      ! Thermal diffusivity
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nzsoil

  INTEGER :: soiltyp(nx,ny)
  INTEGER :: vegtyp(nx,ny)

  REAL :: qsoil (nx,ny,nzsoil)
  REAL :: tsoil (nx,ny,nzsoil)
  REAL :: xqsoil(nx,ny,nzsoil)    ! Volume of nonfrozen soil water 

  REAL :: tsdiffus(nx,ny,nzsoil)  ! Thermal diffusivity
  REAL :: cisoil(nx,ny,nzsoil)    ! Volumetric heat capacity 
  REAL :: ktsoiltop(nx,ny)        ! Thermal conductivity, top soil level 

!
!-----------------------------------------------------------------------
!  Include files:  
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j, k
  INTEGER :: getopt 

  REAL :: psi              ! Matric water potential
  REAL :: psitemp          ! Temporary array
  REAL :: pf               ! log(psi)
  REAL :: ktsoil           ! Thermal conductivity
  REAL :: totaldepth       ! Thickness of total soil column

  REAL :: satker           ! Saturation of soil, qsoil/porosity
  REAL :: kersten          ! Kersten Number
  REAL :: dryden           ! Dry soil density [kg/m3]
  REAL :: liqfrc           ! Nonfrozen liquid fraction
  REAL :: tcs              ! Solids thermal conductivity
  REAL :: tcko             ! Minerals thermal conductivity
  REAL :: tcdry            ! Dry thermal conductivity
  REAL :: tcsat            ! Saturated thermal conductivity
  REAL :: xu               ! Volume of nonfrozen soil water saturation 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!----------------------------------------------------------------------

  DO k=1,nzsoil    !  set initial ground heat flux
    DO j=1,ny-1
      DO i=1,nx-1

      getopt = 2 

!----------------------------------------------------------------------
!   Option #1 for computing soil thermal conductivity
!   Compute soil thermal conductivity using McCumber and Pielke (1981)
!
!   Kt = 420 * EXP[ -(2.7 + Pf)]       Pf <= 5.1
!      = 0.1744                        Pf > 5.1, Pf < 0
!
!   Pf = log[ psi,sat * (Qsat / Q)**b ]
!
!----------------------------------------------------------------------

  IF (getopt == 1) THEN 

        psitemp = 0.0

        IF (qsoil(i,j,k) > 0.0) THEN
           psitemp = wsat(soiltyp(i,j)) / qsoil(i,j,k)
           psi = psisat(soiltyp(i,j))*(psitemp**bslope(soiltyp(i,j)) )

           IF (psi > 0.0) THEN 
             pf = ALOG10(psi)
             IF (pf <= 5.1) ktsoil = 418.46*EXP(-(pf + 2.7))
             IF (pf > 5.1)  ktsoil = 0.172
             IF (pf < 0.0)  ktsoil = 0.172
           ELSE 
             ktsoil = 0.172 
           END IF
         ELSE
             ktsoil = 0.172
         END IF
 
!!!!        IF (ktsoil >= 1.9) ktsoil = 1.90

        IF (k == 1)  ktsoiltop(i,j) = ktsoil

    END IF 

!-------------------------------------------------------------------------


!        psitemp = 0.0
!
!        IF (qsoil(i,j,k) /= 0.0) THEN
!           psitemp = wsat(soiltyp(i,j)) / qsoil(i,j,k)
!        END IF
!
!        psi = psisat(soiltyp(i,j))*(psitemp**bslope(soiltyp(i,j)) )
!
!        pf = ALOG10(psi)
!
!        IF (pf <= 5.1) ktsoil = 418.46*EXP(-(pf + 2.7))
!
!        IF (pf > 5.1)  ktsoil = 0.172
!        IF (pf < 0.0)  ktsoil = 0.172
!     
!!!!!      IF (ktsoil < 0.2) ktsoil = 0.20 
!        IF (ktsoil >= 1.9) ktsoil = 1.90
!
!        IF (k == 1)  ktsoiltop(i,j) = ktsoil
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!   Option #2 for computing soil thermal conductivity
!   Compute soil thermal cond. using Johansen (1975); Peters-Lidard (1998)
!
!   Kt = KE*(Ksat - Kdry) + Kdry
!   Kdry = (0.135*gammadry + 64.7)/(2700 - 0.947*gammadry)
!      gammadry = (1 - porosity) * 2700
!
!   Ksat = Ks**(1-porosity) * Kice**(porosity-xu) * Kh2o**(xu)
!   Ks = Kq**quartz * Ko**(1-quartz)
!
!   KE = 0.7 * log(SR) + 1.0    SR > 0.05  (coarse)
!      = log(SR) + 1.0          SR > 0.10  (fine)
!-------------------------------------------------------------------------

     IF (getopt == 2) THEN 

!     Compute dry density, dryden   
      dryden = (1.0 - porosity(soiltyp(i,j))) * 2700.0   ! [kg/m3]

!     Compute dry thermal conductivity [W/m/K] 
      tcdry = (0.135*dryden + 64.7) / (2700.0 - 0.947*dryden)

      tcko = 3.0                                   ! [ W/m/K]
      IF (quartz(soiltyp(i,j)) > 0.2) tcko = 2.0   ! [ W/m/K]

!     Conductivity of solids (tcs)
      tcs = (7.7**quartz(soiltyp(i,j))) * (tcko**(1.0-quartz(soiltyp(i,j)) ))

!     Estimate unfrozen (liquid) fraction (1.0 is 100% liquid; wwlt 100% frozen
      liqfrc = (xqsoil(i,j,k) + 1.0E-9) / (qsoil(i,j,k) + 1.0E-9) 


!     Estimate unfrozen volume for saturation (xu) 
      xu = liqfrc * porosity(soiltyp(i,j)) 
 

!     Compute saturated thermal conductivity, (tcsat)
      tcsat = (tcs ** (1.0 - porosity(soiltyp(i,j)) )) *              &
             (2.2 ** (porosity(soiltyp(i,j)) - xu) ) *                &
             (0.57** (xu))

!     Compute saturation ratio, (satker) 
      IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN
        satker = qsoil(i,j,k) / porosity(soiltyp(i,j))
      ELSE
        satker = 1.0
      END IF

!!!      IF (satker > 1.0) satker = 1.0
!!!      IF (satker <= 0.1) satker = 0.11

!     Compute Kersten number 
      IF (satker > 0.1) THEN 
        IF ( (xqsoil(i,j,k) + 0.0005) < qsoil(i,j,k)) THEN  !Frozen 
          kersten = satker
        ELSE                                              !Unfrozen 
          kersten = ALOG10(satker) + 1.0    !Kersten Number, fine soils
        END IF 
        IF (soiltyp(i,j) == 12) kersten = satker
      ELSE
        kersten = 0.0 
      END IF 

!     Compute thermal conductivity 
      ktsoil = kersten * (tcsat-tcdry) + tcdry
!!!      IF (ktsoil > 1.9) ktsoil = 1.90
!!!      IF (ktsoil < 0.2) ktsoil = 0.20

      IF (k == 1) ktsoiltop(i,j) = ktsoil

      END IF  

!---------------------------------------------------------------------------
!    Compute heat capacity as a function of soil moisture
!
!    C = Q * C,h20 + (1 - Qsat)*C,soil  +  (Qsat - Q)*C,air
!
!-----------------------------------------------------------------------

      IF (qsoil(i,j,k) > wfc(soiltyp(i,j))) qsoil(i,j,k) = wfc(soiltyp(i,j))

       IF (qsoil(i,j,k) > 1.0) qsoil(i,j,k) = 1.0

!        cisoil(i,j,k) = (qsoil(i,j,k)*cwater) +                     &
!               (1.0 - wsat(soiltyp(i,j))) * csoil +               &
!               (wsat(soiltyp(i,j)) - qsoil(i,j,k)) * cair


        cisoil(i,j,k) = (xqsoil(i,j,k)*cwater) +                    &
               (1.0 - wsat(soiltyp(i,j))) * csoil +               &
               (wsat(soiltyp(i,j)) - qsoil(i,j,k)) * cair         & 
               + ( qsoil(i,j,k) - xqsoil(i,j,k)) * cice 

        tsdiffus(i,j,k) = ktsoil / cisoil(i,j,k)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE getcon


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE TSKIN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE tskin(nx,ny, nzsoil, cdha, rhoa, tair, potair,     &
                rrtem, fftem, ktsoiltop, psfc, evappot, lhflx, tsoil, &
                soiltyp,snowdpth,veg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:  To estimate the top surface (skin) temperature, (Ek and Mahrt, 1991)
!
!       Ts =  (YY + ZZ * Tsoil(x,y,2)) / (ZZ + 1.0)
!
!            FF = (1 - albedo) * SW,in + LW,in
!
!            RCH = rho * Cp * Ch
!
!       YY = Tair + [(FF-sigma*Tair**4)/RCH + (Tp-Tair)- Beta*(Lv*Ep/RCH)
!                   ----------------------------------------------------
!                                          r + 1
!
!       ZZ = Kt / (-0.5*zsoil(1) * RCH * (r + 1) )
!
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  7/19/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nx         Number of grid points in the x-direction (east/west)
!    ny         Number of grid points in the y-direction (north/south)
!    tsoil      Soil temperatures (K)
!    cdha       Drag coefficient for heat 
!    ktsoiltop  Soil conductivity 
!    evappot    Potential evaporation  
!    lhflx      Initial estimates of LE flux 
!    fftem      Radiative forcing (SW,net + LW,in) 
!
!  OUTPUT:
!
!     tsoil(,,1)    ! Potential skin temperature (K)  
!
!-----------------------------------------------------------------------
!
    IMPLICIT NONE

    INTEGER :: nx,ny,nzsoil
    INTEGER :: soiltyp(nx,ny) 

    REAL :: tsoil (nx,ny,nzsoil)
    REAL :: ktsoiltop(nx,ny)        ! Thermal conductivity, top soil level

    REAL :: cdha(nx,ny)             ! Drag coefficient
    REAL :: rhoa(nx,ny)             ! Density of air
    REAL :: tair(nx,ny)             ! Air temperature (K) 
    REAL :: potair(nx,ny)           ! Potential temp of air (K)
    REAL :: rrtem(nx,ny)            ! 
    REAL :: fftem(nx,ny)            ! Radiative forcing, SW net + LW in
    REAL :: psfc(nx,ny)             ! Atmospheric pressure (pa) 
    REAL :: evappot(nx,ny)          ! Potential evaporation 
    REAL :: lhflx(nx,ny)            ! Latent heat flux (W/m2) 
    REAL :: snowdpth(nx,ny)         ! Snow depth (m) 
    REAL :: veg(nx,ny)              ! Vegetation 

!
!-----------------------------------------------------------------------
!  Include files:  
!-----------------------------------------------------------------------
!
    INCLUDE 'grid.inc'
    INCLUDE 'globcst.inc'
    INCLUDE 'phycst.inc'
    INCLUDE 'soilcst.inc'

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
    INTEGER :: i, j

    REAL :: yy
    REAL :: yynum
    REAL :: zz1 
    REAL :: rch 

    REAL :: tempbeta2 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
    DO j=1,ny-1
      DO i=1,nx-1

      rch = rhoa(i,j) * cp * cdha(i,j)

      IF (snowdpth(i,j) == 0.0) THEN 

      IF (lhflx(i,j) > (evappot(i,j)*lathv)) lhflx(i,j)=evappot(i,j)*lathv

      tempbeta2 = 0.0
      IF ((lhflx(i,j) /= 0.0).and.(evappot(i,j) /= 0.0))  THEN
         tempbeta2 = lhflx(i,j) / (evappot(i,j)*lathv)
      END IF

      yynum = fftem(i,j) - sbcst*(tair(i,j)**4.0) 
      yy = tair(i,j) + (yynum/rch + (potair(i,j)-tair(i,j)) - &
           (tempbeta2*lathv*evappot(i,j)/rch)) / (rrtem(i,j)+1.0)

      zz1 = ktsoiltop(i,j)*EXP(-2.0*veg(i,j)) / (dzsoil * &
           rch * (rrtem(i,j)+1.0)) + 1.0  

!     Note that rrtem(i,j) is r; rrtemp is RR in ETA code. 


      IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13 ) THEN 

         tsoil(i,j,1) = (yy + (zz1-1.0) * tsoil(i,j,2)) / zz1 

       END IF
      END IF 
 
      END DO
    END DO

  RETURN
END SUBROUTINE tskin  


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE QDIFF                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE qdiff(nx,ny,nzsoil, j3soilinv, soiltyp, lai, veg,         &
           tsoil, qsoil, wetcanp, windsp, rhoa, precip,              &
           qvair, qvsata, evapcan, evaprg, tem1soil, tem2soil,        &
           tem3soil, tem4soil, precipdrip)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate the soil moisture profile and input matrix
!                for tridiagonalization.

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  3/06/02
!
!  MODIFICATION HISTORY:
!
!  (4/11/02)  Dan Weber
!  Cleaned up code and modified loop calculations to improve
!  optimization.
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the z-direction (soil)
!
!    cdqa     Array for surface drag coefficient for water vapor
!    cdha     Array for surface drag coefficient for heat
!    veg      Vegetation fraction
!
!    qvair    Specific humidity near the surface
!    windsp     Wind speed near the surface
!    rhoa     Air density near the surface
!
!  OUTPUT:
!
!    tem1soil  Tridiagonal matrix
!    tem2soil  Tridiagonal matrix
!    tem3soil  Tridiagonal matrix
!    tem4soil  Tridiagonal matrix
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nzsoil
  INTEGER :: soiltyp(nx,ny)

  REAL :: precip (nx,ny)
  REAL :: precipdrip(nx,ny)
  REAL :: lai    (nx,ny)
  REAL :: veg    (nx,ny)
  REAL :: wetcanp(nx,ny)

  REAL :: windsp (nx,ny)
  REAL :: rhoa   (nx,ny)
  REAL :: evapcan(nx,ny)
  REAL :: evaprg (nx,ny)

  REAL :: qsdiffustop   (nx,ny)
  REAL :: qsdiffusm    (nx,ny)
  REAL :: qsdiffus     (nx,ny)
  REAL :: qsdiffusa    (nx,ny)

  REAL :: qsconducttop (nx,ny)
  REAL :: qsconductm   (nx,ny)
  REAL :: qsconduct    (nx,ny)
  REAL :: qsconducta   (nx,ny)

  REAL :: qvsata       (nx,ny)
  REAL :: qvair        (nx,ny)

  REAL :: j3soilinv(nx,ny,nzsoil)
  REAL :: qsoil (nx,ny,nzsoil)
  REAL :: tsoil (nx,ny,nzsoil)

  REAL :: tem1soil (nx,ny,nzsoil)
  REAL :: tem2soil (nx,ny,nzsoil)
  REAL :: tem3soil (nx,ny,nzsoil)
  REAL :: tem4soil (nx,ny,nzsoil)


!
!-----------------------------------------------------------------------
!  Include file: phycst.inc
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j, k


  REAL :: tempcoefa
  REAL :: tempcoefb
  REAL :: tempcoefc
  REAL :: tempcoefx1
  REAL :: tempcoefx2
  REAL :: tempcoefz1
  REAL :: tempcoefz2

  REAL :: tempel
  REAL :: tempws
  REAL :: tempinf
  REAL :: tema, temb, temc, temd

  REAL :: wrmax
  REAL :: wr2max
  REAL :: vegp
  REAL :: pnet
  REAL :: runoff


!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!----------------------------------------------------------------------
!   Estimate Soil  Moisture
!-----------------------------------------------------------------------
!
! Compute moisture diffusional coefficient (Dn) and hydraulic (Kn)
! conductivity
!
! Note that the constants (Ksat, Psisat, Qsat, b) are all a f(Veg type)
!
!    dQ/dt = d/dz( D * dQ/dz) + dK/dz
!
!          D = - (b * Ksat * Psisat / Q) * (Q / Qsat)**(b+3)
!
!
!          K = Ksat * (Q / Qsat )**(2*b+3)
!
!-----------------------------------------------------------------------


!------------------------------------------------------------------
!     Compute top moisture boundary conditions
!------------------------------------------------------------------

  tema = 0.2 * 1.e-3
  temb = 1.0/dtsfc 
  temc = dtsfc * dzinv

  DO j=1,ny-1
    DO i=1,nx-1

      qsdiffustop(i,j) = - (bslope(soiltyp(i,j)) * kns(soiltyp(i,j)) *    &
            wsat(soiltyp(i,j)) / qsoil(i,j,1)) *                         &
            ( (qsoil(i,j,1)/wsat(soiltyp(i,j))) ** (bslope(soiltyp(i,j)) &
            + 3.0) )

      qsconducttop(i,j) = (kns(soiltyp(i,j)) / wsat(soiltyp(i,j))) *   &
                ( (qsoil(i,j,1)/wsat(soiltyp(i,j)) ) ** (2.0 *        &
                bslope(soiltyp(i,j)) + 2.0) )

      tempel = evaprg(i,j)

      tempws = rhow * ( (qsdiffustop(i,j) * (qsoil(i,j,1) -     &
               qsoil(i,j,2)) * dzsoilinv) + qsconducttop(i,j) )

      wrmax = tema * veg(i,j) * lai(i,j)        ! m
      wr2max = temb * ( wrmax - wetcanp(i,j))   ! m/s

!      vegp = veg(i,j) * precip(i,j)             ! m/s

      vegp = veg(i,j) * (precip(i,j)/1000.0)     ! m/s


      pnet = vegp - evapcan(i,j)                ! m/s
      temd = pnet - wr2max * rhow               ! m/s  
      runoff = MAX( temd, 0.0 )                 ! m/s

!      tempinf = (precip(i,j)-runoff) * rhow     ! kg/m/m/s
      tempinf = (precip(i,j)/1000.0-runoff) * rhow     ! kg/m/m/s


  IF (soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN 

!      qsoil(i,j,1) = qsoil(i,j,1) + temc * (tempws +      &
!               tempel + tempinf)

       qsoil(i,j,1) = qsoil(i,j,1) + (precipdrip(i,j)/rhow) * temc  

      IF (qsoil(i,j,1) > wfc(soiltyp(i,j))) qsoil(i,j,1) = wfc(soiltyp(i,j)) 

      IF (qsoil(i,j,1) > 1.0) qsoil(i,j,1) = 1.0 

  END IF 

      END DO
   END DO

!------------------------------------------------------------------------
!
!   Compute moisture diffusivity and conductivity terms
!     See Smirnova et al.(1997)
!
!     Qdiff = (-b*Knsat*psisat/theta)*(theta/thetasat)**(b+3)
!
!     Qcond = Knsat * (theta/thetasat)**(2b+3)
!
!------------------------------------------------------------------------

  tema = dtsfc*0.5*dzsoilinv2

  DO k=2,nzsoil-1
    tempcoefc = dtsfc * dzsoilinv

    DO j=1,ny-1
      DO i=1,nx-1

       tempcoefa = cnbeta*tema*j3soilinv(i,j,k) 
       tempcoefb = tema*j3soilinv(i,j,k)-tempcoefa 

       qsdiffusm(i,j) = (bslope(soiltyp(i,j)) * kns(soiltyp(i,j)) *       &
           psisat(soiltyp(i,j)) / qsoil(i,j,k-1)) *                       &
          ( (qsoil(i,j,k-1)/wsat(soiltyp(i,j))) ** (bslope(soiltyp(i,j))  &
               + 3.0) )

        qsdiffus (i,j) = (bslope(soiltyp(i,j)) * kns(soiltyp(i,j)) *      &
              psisat(soiltyp(i,j)) / qsoil(i,j,k)) *                      &
             ( (qsoil(i,j,k)/wsat(soiltyp(i,j))) ** (bslope(soiltyp(i,j)) &
               + 3.0) )

        qsdiffusa(i,j) = (bslope(soiltyp(i,j)) * kns(soiltyp(i,j)) *      &
            psisat(soiltyp(i,j)) / qsoil(i,j,k+1)) *                      &
           ( (qsoil(i,j,k+1)/wsat(soiltyp(i,j))) ** (bslope(soiltyp(i,j)) &
               + 3.0) )

         qsconductm(i,j) = (kns(soiltyp(i,j)) / wsat(soiltyp(i,j))) *   &
                ( (qsoil(i,j,k-1)/wsat(soiltyp(i,j)) ) ** (2.0 *        &
                bslope(soiltyp(i,j)) + 2.0) )

         qsconduct (i,j) = (kns(soiltyp(i,j)) / wsat(soiltyp(i,j))) *   &
                ( (qsoil(i,j,k)/wsat(soiltyp(i,j)) ) ** (2.0 *          &
                bslope(soiltyp(i,j)) + 2.0) )

         qsconducta(i,j) = (kns(soiltyp(i,j)) / wsat(soiltyp(i,j))) *   &
                ( (qsoil(i,j,k+1)/wsat(soiltyp(i,j)) ) ** (2.0 *        &
                bslope(soiltyp(i,j)) + 2.0) )

      END DO
    END DO

!------------------------------------------------------------------------
!
!   Compute coefficients for tridiagonalization of moisture variables
!
!------------------------------------------------------------------------


    DO j=1,ny-1
      DO i=1,nx-1

        tempcoefx1 = (qsdiffusm(i,j) * j3soilinv(i,j,k-1)) +        &
                   (qsdiffus(i,j) * j3soilinv(i,j,k  ))

        tempcoefx2 = (qsdiffus(i,j) * j3soilinv(i,j,k  )) +        &
                   (qsdiffusa(i,j) * j3soilinv(i,j,k+1))

        tempcoefz1 = tempcoefc * qsconductm(i,j) * j3soilinv(i,j,k-1)
        tempcoefz2 = tempcoefc * qsconducta(i,j) * j3soilinv(i,j,k+1)

        tem1soil(i,j,k) = - tempcoefa * tempcoefx1 + tempcoefz1
        tem2soil(i,j,k) = 1.0 + tempcoefa * (tempcoefx1 + tempcoefx2)
        tem3soil(i,j,k) = - tempcoefa * tempcoefx2 - tempcoefz2
        tem4soil(i,j,k) = tempcoefb * tempcoefx1 * qsoil(i,j,k-1) +         &
        (1.0-tempcoefb*(tempcoefx1 + tempcoefx2)) * qsoil(i,j,k)     &
         + (tempcoefb * tempcoefx2) * qsoil(i,j,k+1)

      END DO
    END DO

  END DO         ! end of k loop, the moisture coefficients are set

!-----------------------------------------------------------------
! Set upper and lower moisture boundary conditions
! note lower boundary is zero gradient
!   (qsoil(i,j,nzsoil) =  qsoil(i,j,nz-1)
!-----------------------------------------------------------------

  DO j=1,ny-1     !  must set the boundry condition for tem2soil(i,j,nzsoil-1)
    DO i=1,nx-1

      tem4soil(i,j,2) = tem4soil(i,j,2)-tem1soil(i,j,2)*qsoil(i,j,1) !Top BC
      tem2soil(i,j,nzsoil-1) = tem2soil(i,j,nzsoil-1) +     & 
             tem3soil(i,j,nzsoil-1) !Bottom, zero gradient
    END DO
  END DO

  RETURN
END SUBROUTINE qdiff


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CANRES                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE canres(nx,ny, radsw, f34, cdqa,                              &
           windsp,psfc,rhoa,qvair,                                      &
           soiltyp,vegtyp,lai,veg,                                      &
           tair,tsoil,qsoil,wetcanp,nzsoil,zpsoil,snowdpth,             &
           qvsata,rstcoef,plantcoef,cdha)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate canopy resistance.   
!
!
!                     1
!        Ra = ----------------
!               cdq * windsp
!
!                   Rsmin
!        Rs = ------------------
!              LAI*F1*F2*F3*F4
!
!
!               f + Rsmin/Rsmax
!        F1 = -------------------
!                   f + 1
!
!                   Rg    2
!        f = 0.55 ----- -----
!                  Rgl   LAI
!
!             -  1,                             Wfc < W2
!             |
!             |    W2 - Wwlt
!        F2 = -  ------------,              Wwlt <= W2 <= Wfc
!             |   Wfc - Wwlt
!             |
!             -  0,                             W2 < Wwlt
!
!
!               1-0.06*(qvsats-qvair),   qvsats-qvair <= 12.5 g/kg
!        F3 = {
!               0.25,                      otherwise
!
!
!        F4 = 1 - 0.0016 * (298-tair)**2
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  4/20/94
!
!  MODIFICATION HISTORY:
!
!  07/16/02 Jerry Brotzge
!  Modified according to Chen and Dudhia (2001, MWR)  
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the soil
!
!    radsw    Incoming solar radiation
!
!    f34      f3*f4;
!             f3 = Fractional conductance of atmospheric vapor pressure
!             f4 = Fractional conductance of air temperature
!
!    cdqa     Array for surface drag coefficient for water vapor
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    tsoil    Soil temperatures
!    qsoil    Soil moistures
!
!    psfc     Surface pressure
!    qvair    Specific humidity near the surface
!    windsp     Wind speed near the surface
!    rhoa     Air density near the surface
!
!  OUTPUT:
!
!     rstcoef  Canopy resistance (s/m) 
!
!  WORK ARRAY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
 
  INTEGER :: nx,ny
  INTEGER :: nzsoil
 
  REAL :: radsw  (nx,ny)
  REAL :: f34    (nx,ny)
  REAL :: cdqa   (nx,ny)

  INTEGER :: soiltyp(nx,ny)
  INTEGER :: vegtyp (nx,ny)
  REAL :: lai    (nx,ny)
  REAL :: veg    (nx,ny)
  REAL :: tsoil (nx,ny,nzsoil)
  REAL :: qsoil (nx,ny,nzsoil)
  REAL :: wetcanp(nx,ny)
  REAL :: sumgdz(nx,ny)
  REAL :: zpsoil (nx,ny,nzsoil)
  REAL :: tair   (nx,ny)
  REAL :: snowdpth(nx,ny)
  REAL :: psfc   (nx,ny)
  REAL :: qvair  (nx,ny)
  REAL :: windsp   (nx,ny)
  REAL :: rhoa   (nx,ny)
  REAL :: qvsata (nx,ny) 
  REAL :: rstcoef(nx,ny)
  REAL :: plantcoef(nx,ny) 
  REAL :: cdha(nx,ny) 

!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: pi
  PARAMETER ( pi = 3.141592654 )

  INTEGER :: i, j, k

  REAL :: wrmax       ! Maximum value for canopy moisture, wetcanp
 
  REAL :: temb, temc        ! Used to compute f1*f2
  REAL :: f2          ! Temporary variable for f2

  REAL :: waf         ! Available soil moisture fraction
  REAL :: es          ! saturation vapor pressure 
  REAL :: beta2       ! Used to estimate f1*f2
  REAL :: delzneg

  REAL :: desdt, dqsdt, st1temp, rrtemp, deltatem

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j=1,ny-1
    DO i=1,nx-1
      sumgdz(i,j) = 0.0
      rstcoef(i,j) = 0.0
    END DO
  END DO

  DO j = 1, ny-1
    DO i = 1, nx-1


    IF ( soiltyp(i,j) /= 12 .AND. soiltyp(i,j) /= 13) THEN 

        wrmax = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! in meter
!
!-----------------------------------------------------------------------
!
!  In order to calculate the qv flux at the surface, we need to
!  calculate some parameters like the resistance coefficient,
!
!      rstcoef = rsta / (rsts + rsta)
!
!  where rsta is the aerodynamic resistance and rsts is the surface
!  resistances.
!
!      rsta = 1 / ( cdq * Va )
!
!      rsts = rsmin(vtyp) / ( lai(vtyp) * f1 * f2 * f3 * f4 )
!
!  f3 * f4 is time-independent and has been calculated previously
!  and stored in f34(i,j)
!
!-----------------------------------------------------------------------
!
!  Calculate f1.
!
!        f = 0.55*(radsw/rgl(vtyp))*(2./lai(vtyp))       ! NP, Eq. 34
!        f1 = ( rsmin(vtyp)/rsmax + f ) / ( 1. + f )     ! NP, Eq. 34
!
!  Note: the incoming solar radiation radsw is stored in radsw(i,j).
!
!-----------------------------------------------------------------------
!
        IF ( lai(i,j) == 0. ) THEN        ! No vegetation, I/W
          rstcoef(i,j) = 1.
        ELSE

          IF (radsw(i,j) < 0.0 ) radsw(i,j) = 0.0

          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )          &
               * ( 2.0 / lai(i,j) )

          rstcoef(i,j) = ( rsmin(vegtyp(i,j)) / rsmax + temb )     &
                  / ( 1.0 + temb )

          rstcoef(i,j) = MAX( 0.0001, rstcoef(i,j)) 
        END IF


!---------------------------------------------------------------------------
!
!   Calculate F4.  
!   Replaced force-restore wetdp estimate with integrated soil moisture (JAB)
!     estimated as a function of root zone depth 
!   (See Chen and Dudhia MWR 2001)
!
!---------------------------------------------------------------------------

  DO k = 2,nzsoil-1


        IF ((zpsoil(i,j,1)-zpsoil(i,j,nzsoil)) <= rootzone(vegtyp(i,j))) & 
           rootzone(vegtyp(i,j)) = zpsoil(i,j,1) - zpsoil(i,j,nzsoil)  

        IF (zpsoil(i,j,k) >= (zpsoil(i,j,1)-rootzone(vegtyp(i,j))) ) THEN 
          delzneg = (zpsoil(i,j,k-1)-zpsoil(i,j,k))/              &
                  rootzone(vegtyp(i,j))  
        ELSE  
          delzneg = 0.0
        END IF 

        beta2 = (qsoil(i,j,k) - wwlt(soiltyp(i,j))) /           &
            (wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )

        IF (qsoil(i,j,k) > wfc(soiltyp(i,j)) )  beta2 = 1.0
        IF (qsoil(i,j,k) <= wwlt(soiltyp(i,j))) beta2 = 0.0

        sumgdz(i,j) = sumgdz(i,j) + (beta2 * delzneg)

  END DO

        sumgdz(i,j) = MAX( 0.0001, sumgdz(i,j))  

        rstcoef(i,j)  = rstcoef(i,j) * sumgdz(i,j)



!
!-----------------------------------------------------------------------
!
! Calculate F2. 
! 
!-----------------------------------------------------------------------

      f2 = 1.0 + hsf2(vegtyp(i,j))*(qvsata(i,j) - qvair(i,j) )  
      IF (qvair(i,j) > qvsata(i,j) .OR. f2 < 0) f2 = 1.0 
      IF (f2 > 1.0E-30) THEN    
        f2 = 1.0/f2  
        f2 = MAX( 0.01, f2) 
      ELSE IF (f2 <= 1.0E-30 .AND. f2 /= 1.0) THEN 
        f2 = 0.01 
      END IF 

!
!-----------------------------------------------------------------------
!
!  Calculate f3 * f4, stored in f34(i,j), where
!
!    f3 -- Fractional conductance of atmospheric vapor pressure
!      f3 = 1 - 0.06 * ( qvsata(i,j) - qvair ) * 1.e3  ! use kg/kg for qv
!
!    (We use f3 = 1 in the code instead of the above formula.)
!
!    f4 -- Fractional conductance of air temperature
!      f4 = 1 - 0.0016 * ( 298 - tanem ) ** 2
!
!  f3 * f4 will be used to calculate the resistance coefficient.
!
!-----------------------------------------------------------------------
!
!      f3 = 1.0                                 ! f3 set to 1
!      f34(i,j) = f3 * ( 1.0 - 0.0016 * ( 298. - tair(i,j) ) ** 2 )
!
!-----------------------------------------------------------------------
!

        f34(i,j) = MAX( 0.0001, 1.0 - 0.0016 * (298.0-tair(i,j))**2 )
                                                   ! f3 * f4, JN, Eq. 11

!
!-----------------------------------------------------------------------
!
!  Calculate lai*f1*f2*f3*f4 where f3*f4 is stored in f34(i,j).
!
!-----------------------------------------------------------------------
!

        rstcoef(i,j) = lai(i,j)*rstcoef(i,j)*f2*f34(i,j)      ! lai*f1*f2*f3*f4

!
!-----------------------------------------------------------------------
!
!  Calculate the resistance coefficient, rsta/(rsts+rsta)
!
!        rsts = rsmin(vtyp)/(lai(i,j)*f1*f2*f3*f4) ! Sfc. resistance
!        rsta = 1./(cdh*va)                 ! NP, between Eq. 32 & 33
!        rstcoef = rsta/(rsta+rsts)
!                = 1/(1+rsts/rsta)
!
!-----------------------------------------------------------------------
!
         IF ( ABS(rstcoef(i,j)) > 1.0E-30) THEN
            rstcoef(i,j) = rsmin(vegtyp(i,j)) / rstcoef(i,j)
         END IF

!-------------------------------------------------------------------------
!  Calculate PC - Plant coefficient (PC)
!-------------------------------------------------------------------------

         temc = lathv * cpinv 
         desdt = 0.0
         dqsdt = 0.0 
         DO k = 7,2,-1
           desdt = alpha(k)*(k-1) + (tair(i,j)-273.16)*desdt
         END DO

         dqsdt = 0.622 * desdt / psfc(i,j)

!  TEST **********************************************************
!         es = 6.112*exp( 17.67*(tair(i,j)-273.15)/(tair(i,j)-29.65) )
!         desdt = (es*2.5e6) / (461.0 * tair(i,j) * tair(i,j) ) 
!         dqsdt = desdt * (0.622 * psfc(i,j) / ( (psfc(i,j) -         &
!             0.378*es)**2.0) ) 

         deltatem = temc * dqsdt

         st1temp = 4.0 * sbcst * rd * cpinv  

         rrtemp = (st1temp * (tair(i,j)**4.0)) / (psfc(i,j)*cdha(i,j)) &
           + 1.0  

         plantcoef(i,j) = (rrtemp*deltatem) / (rrtemp             & 
            * (1.0 + rstcoef(i,j) * cdha(i,j)) + deltatem) 


      END IF      !soiltyp /= 12 and 13  

    END DO
  END DO

  RETURN
END SUBROUTINE canres  



!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE LEFLX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE leflx(nx,ny,nzsoil,cdha, cdqa, windsp, rhoa, qvair, veg,  &
           wetcanp,deltem, rrtem, temple, lhflx, evappot, soiltyp,   &
           vegtyp, evaprg, evaprtr, evapcan, qvsfc, qvsat, rsttemp,  &
           plantcoef,zpsoil, qsoil, xqsoil, qice, precip1, precipdrip)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate the latent heat flux.
!
!  Calculate:
!
!  1. Evaporation from ground surface,
!
!        evaprg = rhoa * cdq * windsp * evaprg'
!
!  where
!
!        evaprg' = (1.-veg) * (rhgs * qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  2. Direct evaporation from the fraction delta of the foliage covered
!     by intercepted water.
!
!        Er = rhoa * cdq * windsp * Er'
!
!        Er' = delta * veg * (qvsats-qvair)
!
!  3. Transpiration of the remaining part (1-delta) of leaves,
!
!        Etr = rhoa * cdq * windsp * Etr'
!
!  where
!
!        Etr' = veg * (1-delta) * Ra/(Ra+Rs) * (qvsats-qvair )
!
!  and Ra is aerodynamic resistance and Rs is the surface resistance
!
!                     1
!        Ra = ----------------
!               cdq * windsp
!
!                   Rsmin
!        Rs = ------------------
!              LAI*F1*F2*F3*F4
!
!  4. Water vapor flux, lhflx,
!
!        lhflx = rhoa * cdq * windsp * qvflx'
!
!        lhflx' = (Eg' + Etr' + Er') = (qvsfc - qvair)
!
!     where qvsfc is the effective surface specific humidity
!
!        (qvsfc - qvair) = (Eg' + Etr' + Er')
!
!  This subroutine will solve Eg', Etr', Er', and leflx'
!
!-----------------------------------------------------------------------
!
!  AUTHOR: J. Brotzge
!  3/06/02
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    cdqa     Array for surface drag coefficient for water vapor
!    cdha     Array for surface drag coefficient for heat
!    veg      Vegetation fraction
!
!    qvair    Specific humidity near the surface
!    windsp     Wind speed near the surface
!    rhoa     Air density near the surface
!
!  OUTPUT:
!
!    evaprp   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evapcan  Direct evaporation from the fraction delta
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny, nzsoil 
  INTEGER :: soiltyp(nx,ny)
  INTEGER :: vegtyp (nx,ny)

  REAL :: cdqa   (nx,ny)
  REAL :: cdha   (nx,ny)
  REAL :: veg    (nx,ny)
  REAL :: wetcanp(nx,ny)

  REAL :: qvair  (nx,ny)
  REAL :: qvsfc  (nx,ny)
  REAL :: qvsat  (nx,ny) 
  REAL :: windsp (nx,ny)
  REAL :: rhoa   (nx,ny)
  REAL :: rsttemp(nx,ny)

  REAL :: evappot(nx,ny)
  REAL :: evaprg (nx,ny)
  REAL :: evaprtr(nx,ny)
  REAL :: evapcan(nx,ny)
  REAL :: lhflx  (nx,ny)

  REAL :: deltem (nx,ny)
  REAL :: rrtem  (nx,ny)
  REAL :: temple (nx,ny)

  REAL :: beta3  (nzsoil) 
  REAL :: rtdis  (nzsoil) 
  REAL :: qice(nx,ny,nzsoil)
  REAL :: precip1(nx,ny)  
  REAL :: qsoil  (nx,ny,nzsoil)
  REAL :: xqsoil (nx,ny,nzsoil)  
  REAL :: zpsoil (nx,ny,nzsoil) 
  REAL :: plantcoef(nx,ny) 
  REAL :: precipdrip(nx,ny) 

!
!-----------------------------------------------------------------------
!  Include file: phycst.inc
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc' 
  INCLUDE 'soilcst.inc'
 

!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j, k
  INTEGER :: nzsoiltemp
 
  REAL :: tempbc, tempfx, etp1a, sgxtemp 
  REAL :: drip, dew, trhsct, excess
  REAL :: eta, eta1, wetcan2, rtxtemp, denomtemp  
  REAL :: beta2, tempbeta  

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

    DO j = 1, ny-1
      DO i = 1, nx-1

      IF (evappot(i,j) > 0.0) THEN

!--------------------------------------------------------------------
!    Compute bare soil evaporation
!--------------------------------------------------------------------

        IF (veg(i,j) < 1.0) THEN 

          IF (soiltyp(i,j) /= 13) THEN
            tempbeta = (xqsoil(i,j,1) - wwlt(soiltyp(i,j))) /          &
              (wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )
          ELSE IF (soiltyp(i,j) == 13) THEN
            tempbeta = 1.0
          END IF

          tempbeta = min(tempbeta,1.) 

          IF (tempbeta > 0.0) THEN 
            tempfx = tempbeta**2.0 
            tempfx = max( min(tempfx,1.0), 0.0)
          ELSE 
            tempfx = 0.0
          END IF  
       
          evaprg(i,j) = tempfx * (1.0 - veg(i,j) ) * evappot(i,j) 

        ENDIF


!--------------------------------------------------------------------
!    Compute evapotranspiration
!--------------------------------------------------------------------

        sgxtemp = 0.0
        denomtemp = 0.0  
        nzsoiltemp = 0 
        evaprtr(i,j) = 0.0  
        DO k=1,nzsoil
          beta3(k) = 0.0
        END DO


        IF (veg(i,j) > 0.0) THEN
           
          IF (wetcanp(i,j) /= 0.0) THEN
            etp1a = veg(i,j) * plantcoef(i,j) * evappot(i,j) *      &
              (1.0 - (wetcanp(i,j) / 0.5E-3)**0.5)  
          ELSE
            etp1a = veg(i,j) * plantcoef(i,j) * evappot(i,j)
          END IF 

          sgxtemp = 0.0
          nzsoiltemp = 0
 
          DO k=1,nzsoil
       
            beta3(k) = (xqsoil(i,j,k) - wwlt(soiltyp(i,j))) /           &
              (wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )

            beta3(k) = max( min(beta3(k),1.0), 0.0)   !New******* 

            IF (zpsoil(i,j,nzsoil) >= rootzone(vegtyp(i,j))) THEN
              sgxtemp = sgxtemp + beta3(k) 
              nzsoiltemp = nzsoiltemp + 1  
            END IF
          END DO

          sgxtemp = sgxtemp / real(nzsoiltemp) 

          denomtemp = 0.0

          IF (nzsoiltemp > 1) THEN 
            rtdis(1) = 0.0 
            beta3(1) = 0.0 

            DO k=2,nzsoiltemp 

              IF ((zpsoil(i,j,1)-zpsoil(i,j,nzsoil)) <= rootzone(vegtyp(i,j))) &
                rootzone(vegtyp(i,j)) = zpsoil(i,j,1) - zpsoil(i,j,nzsoil)

              IF (zpsoil(i,j,k) >= (zpsoil(i,j,1)-rootzone(vegtyp(i,j))) ) THEN  
                rtdis(k) = (zpsoil(i,j,k-1)-zpsoil(i,j,k))/            &
                  rootzone(vegtyp(i,j)) 

              ELSE IF ( (zpsoil(i,j,k-1) > (zpsoil(i,j,1)-rootzone(vegtyp(i,j)) )) &
                .AND. (zpsoil(i,j,k) < (zpsoil(i,j,1)-rootzone(vegtyp(i,j)) )) ) THEN
                
                rtdis(k) = (zpsoil(i,j,k-1)-rootzone(vegtyp(i,j))) /    &
                  (zpsoil(i,j,k-1) - zpsoil(i,j,k)) 

              ELSE
                rtdis(k) = 0.0 
              END IF 

              rtxtemp = rtdis(k) + beta3(k) - sgxtemp
              beta3(k) = beta3(k) * max( rtxtemp, 0.0) 
              denomtemp = denomtemp + beta3(k) 
            END DO 
            IF (denomtemp <= 0.0) denomtemp = 1.0
            DO k=2,nzsoiltemp
              evaprtr(i,j) = evaprtr(i,j) + etp1a * beta3(k) / denomtemp
            END DO   
          END IF
          IF (nzsoiltemp == 1) evaprtr(i,j) = 0.0  

!--------------------------------------------------------------------------
!    Compute direct Canopy Evaporation
!--------------------------------------------------------------------------

          IF (wetcanp(i,j) > 0.0) THEN
            evapcan(i,j) = veg(i,j) * ( (wetcanp(i,j)/0.5E-3)**0.5) * evappot(i,j)    
          ELSE
            evapcan(i,j) = 0.0
          END IF

          wetcan2 = wetcanp(i,j) / dtsfc 
          evapcan(i,j) = MIN(wetcan2, evapcan(i,j))

        END IF       ! (veg > 0) 

      ELSE
          evaprg(i,j) = 0.0
          evaprtr(i,j)= 0.0
          evapcan(i,j) = 0.0

      END IF         ! (ETP > 0)  


!-----------------------------------------------------------------------
!    Compute Total Evaporation
!-----------------------------------------------------------------------

      eta1 = evaprg(i,j) + evaprtr(i,j) + evapcan(i,j) 

      lhflx(i,j) = evaprg(i,j) + evaprtr(i,j) + evapcan(i,j)
      lhflx(i,j) = lhflx(i,j) * lathv       ! Latent heat flux

      trhsct = (veg(i,j) * precip1(i,j)*1000.0 - evapcan(i,j)) * dtsfc  
      drip = 0.0
      excess = wetcanp(i,j) + trhsct
      IF (excess > 0.5E-3) drip = excess - 0.5E-3 

!--------------------------------------------------------------------------
!    Compute total precip seeping into the ground.
!--------------------------------------------------------------------------

      precipdrip(i,j) = (1.0 - veg(i,j)) * precip1(i,j)*1000.0 + drip/dtsfc

!--------------------------------------------------------------------------
!    Compute fraction of soil moisture that is ice
!--------------------------------------------------------------------------

      DO k = 1,nzsoil
        qice(i,j,k) = qsoil(i,j,k) - xqsoil(i,j,k) 
      END DO   

      END DO
    END DO

  RETURN
END SUBROUTINE leflx   


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SNOPHY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE snophy(nx,ny,nzsoil,tsoil,tair,potair,precip,evappot,     &
                 snowdpth,sneqv,sndens,snowng,snofac,fftem,rrtem,    &
                 ktsoiltop,gflx,psfc,precip1,qvair,cdha,rhoa,flx2,beta)


!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate snow processes at the land surface.
!
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  12/19/02
!
!  Code obtained from the ETA model, courtesy of NCEP. 
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    sneqv    Snow water equivalent (m)
!    sndens   Snow density
!    snowdpth Snow depth (m)
!    tsoil    Soil temperature; tsoil(1) = tskin
!
!  OUTPUT:
!
!    snowdpth Snow depth (m)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nzsoil
  INTEGER :: snowng(nx,ny) 

  REAL :: sneqv    (nx,ny)
  REAL :: sndens   (nx,ny)
  REAL :: snowdpth (nx,ny)
  REAL :: snofac   (nx,ny) 
  REAL :: tsoil    (nx,ny,nzsoil)
  REAL :: ktsoiltop(nx,ny)

  REAL :: fftem    (nx,ny)
  REAL :: rrtem    (nx,ny)
  REAL :: tair     (nx,ny)
  REAL :: potair   (nx,ny)
  REAL :: gflx     (nx,ny)
  REAL :: psfc     (nx,ny) 
  REAL :: evappot  (nx,ny) 
  REAL :: precip   (nx,ny)      ! [m s-1]
  REAL :: precip1  (nx,ny)      ! [m s-1] 
  REAL :: qvair    (nx,ny)
  REAL :: cdha     (nx,ny)
  REAL :: rhoa     (nx,ny) 
  REAL :: etp1     (nx,ny)
  REAL :: etp2     (nx,ny) 
  REAL :: flx2     (nx,ny) 
  REAL :: beta     (nx,ny) 

!
!-----------------------------------------------------------------------
!  Include file: phycst.inc
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc' 
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc' 

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j 

  REAL :: rch, flx1, flx3, depthtotal
  REAL :: denom, t12a, t12b, t12, snomeltrt, snomeltamt
  REAL :: seh, temd, dew 
  REAL :: etp3, qsat, rsnow, albedo



!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

   DO j = 1, ny-1
     DO i = 1, nx-1

         rch = rhoa(i,j) * cp * cdha(i,j)

!-----------------------------------------------------------------
!     Convert ETP from (kg m-2 s-1) to (m/s) and then to (m)
!     This is for snow-free surfaces.  
!-----------------------------------------------------------------

         IF (sneqv(i,j) == 0.0) THEN            ! No snow cover 

           etp1(i,j) = evappot(i,j) * 0.001 
           dew = 0.0 

           IF (etp1(i,j) < 0.0) THEN    !Dew formation 
             dew = -etp1(i,j)                   ! [m/s]  
             etp1(i,j) = 0.0
             precip1(i,j) = precip1(i,j) + dew  ! [m/s]  
           END IF 

!-----------------------------------------------------------------
!     Convert ETP from (kg m-2 s-1) to (m/s) and then to (m)
!     This is then the "Effective snowpack reduction amount", etp2 (m)
!----------------------------------------------------------------

         ELSE                                  ! Snow cover 

           etp2(i,j) = evappot(i,j) * dtsfc * 0.001
           beta(i,j) = 1.0 

           IF (sneqv(i,j) < etp2(i,j)) THEN
             beta(i,j) = sneqv(i,j) / etp2(i,j)
           END IF

!--------------------------------------------------------------
!     Frost 
!--------------------------------------------------------------
           dew = 0.0
           IF (evappot(i,j) < 0.0) THEN
             dew = - evappot(i,j) * 0.001
           END IF
    
         END IF                              ! Snow cover/No snow cover  

         IF (sneqv(i,j) /= 0.0) THEN         ! Snow cover only  

           etp1(i,j) = 0.0  !No ET if snow cover

!-----------------------------------------------------------------
!     Include heat flux from falling precip
!------------------------------------------------------------------
           flx1 = 0.0                        ! [W/m2] 

           IF (snowng(i,j) /= 1) THEN

!             IF (precip(i,j) > 0.0) flx1 = 4.218E3 * precip(i,j) *     &
!               1000.0 * (tsoil(i,j,1) - tair(i,j))  ! [W/m2] 

             IF (precip(i,j) > 0.0) flx1 = 4.218E3 * precip(i,j) *     &
               (tsoil(i,j,1) - tair(i,j))  ! [W/m2] 


           ELSE
!             flx1 = 2.106E3 * precip(i,j)*1000.0 *(tsoil(i,j,1)-tair(i,j)) 

             flx1 = 2.106E3 * precip(i,j) * (tsoil(i,j,1)-tair(i,j)) 



           END IF

           depthtotal = snowdpth(i,j) + dzsoil   ! [m] 

!---------------------------------------------------------------------
!     Calculate an "effective snow-ground sfc temp"
!---------------------------------------------------------------------

           denom = 1.0 + ktsoiltop(i,j) / (depthtotal*(rrtem(i,j)+1.0)*rch)
           t12a = ((fftem(i,j) - flx1 - flx2(i,j) - sbcst * (tair(i,j)**4.0))/ &
                rch + potair(i,j) - tair(i,j) - beta(i,j) * evappot(i,j) *     &
                lathv / rch ) / (rrtem(i,j)+1.0)
           t12b = ktsoiltop(i,j) * tsoil(i,j,2) / (depthtotal * (rrtem(i,j) &
                +1.0) * rch)
           t12 = (tair(i,j) + t12a + t12b) / denom


!----------------------------------------------------------------------------
!      If soil temp below freezing, then no snow melt; set skin temp to
!        "effective temp", and set the effective precip to zero.
!----------------------------------------------------------------------------

           IF (t12 <= 273.15) THEN
             sneqv(i,j) = MAX(0.0, sneqv(i,j) - etp2(i,j))
             snowdpth(i,j) = sneqv(i,j) / sndens(i,j) 
             tsoil(i,j,1) = t12

!      Update soil heat flux using new skin temperature (tsfc)
             gflx(i,j) = ktsoiltop(i,j) * (tsoil(i,j,1) - tsoil(i,j,2)) / &
                  depthtotal
             flx3 = 0.0
             snomeltrt = 0.0
             snomeltamt = 0.0

!-------------------------------------------------------------------------
!      If soil temps above freezing, then snow melts
!-------------------------------------------------------------------------

           ELSE
             tsoil(i,j,1) = 273.15 * snofac(i,j)+t12*(1.0-snofac(i,j))
             qsat = (0.622*6.11E2) / (psfc(i,j)-0.378*6.11E2)
             evappot(i,j) = rch * (qsat - qvair(i,j)) * cpinv
             etp2(i,j) = evappot(i,j) * 0.001 * dtsfc
             beta(i,j) = 1.0


!-------------------------------------------------------------------------
!       Sublimation
!-------------------------------------------------------------------------
             IF (sneqv(i,j) <= etp2(i,j)) THEN  !More evaporation than snow equiv 
               beta(i,j) = sneqv(i,j) / etp2(i,j) 
               sneqv(i,j) = 0.0
               snowdpth(i,j) = 0.0
               snomeltamt = 0.0
               snomeltrt = 0.0
               gflx(i,j) = ktsoiltop(i,j) * (tsoil(i,j,1)-tsoil(i,j,2))/ &
                      depthtotal

             ELSE                           !More snow equiv than evaporation
               sneqv(i,j) = sneqv(i,j) - etp2(i,j) 
               snowdpth(i,j) = sneqv(i,j) / sndens(i,j)
               etp3 = evappot(i,j) * 2.501E6
               gflx(i,j) = ktsoiltop(i,j) * (tsoil(i,j,1)-tsoil(i,j,2))/ &
                    depthtotal
               seh = rch * (tsoil(i,j,1) - potair(i,j))
               flx3 = fftem(i,j)-flx1-flx2(i,j)-sbcst*(tsoil(i,j,1)**4.0) - &
                 gflx(i,j) - seh - etp3
               IF (flx3 <= 0.0) flx3 = 0.0
               snomeltrt = flx3 * 0.001 / 3.335E5

               IF (snofac(i,j) > 0.05) snomeltrt = snomeltrt * snofac(i,j)
               snomeltamt = snomeltrt * dtsfc
             END IF

!----------------------------------------------------------------------
!         If snow melt amount is less than the snow equivalent
!-----------------------------------------------------------------------

             temd = sneqv(i,j) - 1.0E-6

             IF (snomeltamt < temd) THEN
               sneqv(i,j) = sneqv(i,j) - snomeltamt
               snowdpth(i,j) = sneqv(i,j) / sndens(i,j)
             ELSE
               snomeltrt = sneqv(i,j) / dtsfc     ! [m/s] 
               snomeltamt = sneqv(i,j)            ! [m] 
               sneqv(i,j) = 0.0
               snowdpth(i,j) = 0.0
               flx3 = snomeltrt * 1000.0 * 3.335E5
             END IF
!             precip1(i,j) = precip1(i,j)/1000.0 + snomeltamt 

             precip1(i,j) = precip1(i,j) + snomeltrt  ! [m/s]  

           END IF              ! soil temps <> freezing
         END IF                ! presence of snow cover  
       END DO
     END DO

  RETURN
END SUBROUTINE snophy


!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SNOCOM                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE snocom(nx, ny, nzsoil, snowdpth, sneqv, sndens, sncond, tsoil) 

!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate snow compaction and its effects. 
!
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  12/19/02
!
!  Code adapted from the ETA model, courtesy of NCEP,  
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    sneqv    Snow water equivalent (m) 
!    sndens   Snow density 
!    snowdpth Snow depth (m)
!    sncond   Snow conductivity  
!    tsoil    Soil temperature; tsoil(1) = tskin   
!
!  OUTPUT:
!
!    snowdpth Snow depth (m) 
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nzsoil 

  REAL :: sneqv    (nx,ny)
  REAL :: sndens   (nx,ny) 
  REAL :: snowdpth (nx,ny) 
  REAL :: sncond   (nx,ny)
  REAL :: tsoil    (nx,ny,nzsoil) 


!
!-----------------------------------------------------------------------
!  Include file: phycst.inc
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j, l 

  REAL :: tsoilx, tsnowx, snowdpthx, wx, wxx
  REAL :: tavgtemp, btemp, pexp, dw, dex, dsx 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

   DO j = 1, ny-1
     DO i = 1, nx-1

     IF (sneqv(i,j) > 0.0) THEN                  !Snow cover

       tsnowx = tsoil(i,j,1) - 273.15
       tsoilx = tsoil(i,j,2) - 273.15
       snowdpthx = snowdpth(i,j) * 100.0

       tavgtemp = 0.5 * (tsnowx + tsoilx)
       wx = 100.0 * sneqv(i,j)
       wxx = 1.0E-2
       IF (wx > wxx) wxx = wx
       btemp = (dtsfc/3600.0) * 0.01 * EXP(0.08*tavgtemp -        &
         21.0 * sndens(i,j))

       pexp = 0.0
       DO l=4,1,-1
         pexp = (1.0 + pexp)*btemp*wxx/real(l+1)
       END DO
       pexp = pexp + 1.0
       dsx = sndens(i,j) * pexp
       IF (dsx > 0.40) dsx = 0.40
       IF (dsx < 0.05) dsx = 0.05
       sndens(i,j) = dsx

       IF (tsnowx >= 0.0) THEN
         dw = 0.13 * dtsfc/3600.0/24.0
         sndens(i,j) = sndens(i,j) * (1.0 - dw) + dw
         IF (sndens(i,j) > 0.40) sndens(i,j) = 0.40
       END IF

       snowdpthx = sneqv(i,j) * 100.0 / sndens(i,j)
       snowdpth(i,j) = snowdpthx * 0.01

     ELSE                                        !No snow cover
       sneqv(i,j) = 0.0
       snowdpth(i,j) = 0.0
       sndens(i,j) = 0.0
       sncond(i,j) = 1.0
     END IF

     END DO
   END DO

  RETURN
END SUBROUTINE snocom


!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SNKSRC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE snksrc(nx, ny, nzsoil, tsoil, qsoil, xqsoil, ktsoiltop, &
                    gflx, zpsoil, sink, soiltyp, tem1soil, tem2soil, &
                    tem3soil, tem4soil, tsdiffus, j3soilinv, cisoil) 

!-----------------------------------------------------------------------
!
!  PURPOSE:  To calculate the source and sink terms for the thermal diffusion
!              equation due to freezing and thawing of the soil.  
!
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  01/08/03
!
!  Code adapted from the ETA model, courtesy of NCEP,
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    sneqv    Snow water equivalent (m)
!    sndens   Snow density
!    snowdpth Snow depth (m)
!    sncond   Snow conductivity
!    tsoil    Soil temperature; tsoil(1) = tskin
!
!  OUTPUT:
!
!    snowdpth Snow depth (m)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nzsoil
  INTEGER :: soiltyp(nx,ny) 

  REAL :: sneqv    (nx,ny)
  REAL :: sndens   (nx,ny)
  REAL :: snowdpth (nx,ny)
  REAL :: sncond   (nx,ny)
  REAL :: ktsoiltop(nx,ny)
  REAL :: gflx     (nx,ny) 

  REAL :: tsoil    (nx,ny,nzsoil)
  REAL :: qsoil    (nx,ny,nzsoil)
  REAL :: xqsoil   (nx,ny,nzsoil) 
  REAL :: zpsoil   (nx,ny,nzsoil) 
  REAL :: sink     (nx,ny,nzsoil) 

  REAL :: dza      (nzsoil)
  REAL :: dzh      (nzsoil) 
 
  REAL :: tem1soil (nx,ny,nzsoil)
  REAL :: tem2soil (nx,ny,nzsoil)
  REAL :: tem3soil (nx,ny,nzsoil)
  REAL :: tem4soil (nx,ny,nzsoil)
  REAL :: tsdiffus (nx,ny,nzsoil)
  REAL :: j3soilinv(nx,ny,nzsoil) 
  REAL :: cisoil   (nx,ny,nzsoil) 
 
!
!-----------------------------------------------------------------------
!  Include files 
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'soilcst.inc' 
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc' 

!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i, j, k
  INTEGER :: nlog, kcount 

  REAL :: dtsdz, qtotal, xh2o 
  REAL :: tup, tdn, tm  
  REAL :: xup, xdn, xtemp, tavg 
  REAL :: fk, df, denom, frh2o, swl, swlk 
  REAL :: bx, dswl, error, ck, sice    

  REAL :: tempcoefa, tempcoefb
  REAL :: tempcoefx1, tempcoefx2 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  error = 0.005 
  ck = 8.0 

   DO j = 1, ny-1
     DO i = 1, nx-1
       DO k = 2,nzsoil 
         dza(k) = zpsoil(i,j,k-1) - zpsoil(i,j,k) 
         dzh(k) = 0.5 * dza(k) 
       END DO
     END DO
   END DO 

!------------------------------------------------------------------------
!   Estimate Soil Temperature
!-----------------------------------------------------------------------
!
!  C * dT/dt = d/dz (K * dT/dz)
!
!-----------------------------------------------------------------------
!
!  Calculate coefficients of the tridigonal equation for soil temperature
!
!-----------------------------------------------------------------------
!

  DO k=2,nzsoil-1
    DO j=1,ny-1
      DO i=1,nx-1

      tempcoefa = cnbeta * dtsfc * 0.5 * dzsoilinv2 * j3soilinv(i,j,k)
      tempcoefb = (1.0 - cnbeta) * dtsfc * 0.5 * dzsoilinv2 * j3soilinv(i,j,k)

      tempcoefx1 = (tsdiffus(i,j,k-1) * j3soilinv(i,j,k-1)) +        &
                   (tsdiffus(i,j,k  ) * j3soilinv(i,j,k  ))

      tempcoefx2 = (tsdiffus(i,j,k  ) * j3soilinv(i,j,k  )) +        &
                   (tsdiffus(i,j,k+1) * j3soilinv(i,j,k+1))

      tem1soil(i,j,k) = -tempcoefa * tempcoefx1
      tem2soil(i,j,k) = 1.0 + tempcoefa * (tempcoefx1 + tempcoefx2)
      tem3soil(i,j,k) = -tempcoefa * tempcoefx2

      tem4soil(i,j,k) = tempcoefb * tempcoefx1 * tsoil(i,j,k-1) +  &
        (1.0-tempcoefb*(tempcoefx1 + tempcoefx2))*tsoil(i,j,k)       &
        + tempcoefb * tempcoefx2 *tsoil(i,j,k+1)


      denom = (zpsoil(i,j,k-1)-zpsoil(i,j,k)) * cisoil(i,j,k)
      qtotal = -1.0 * denom * tem4soil(i,j,k)

      sice = qsoil(i,j,k) - xqsoil(i,j,k)

      IF ( (sice > 0.0).OR.(tsoil(i,j,1) < 273.15) .OR.              &
           (tsoil(i,j,2)<273.15).OR.(tsoil(i,j,nzsoil)<273.15) ) THEN

!---------------------------------------------------------------------
!    Compute sink/source term for freezing and thawing
!---------------------------------------------------------------------


!----------------------------------------------------------------------
!    Calculate potential reduction of liquid water content
!----------------------------------------------------------------------

         xh2o = xqsoil(i,j,k) + (qtotal * dtsfc)/(3335.0E5*dza(k))  
 
!-----------------------------------------------------------------------
!    Estimate unfrozen water at an avg temperature
!-----------------------------------------------------------------------

         tup = tsoil(i,j,k-1)
         tdn = tsoil(i,j,k-1) 
         tm = 0.5 * (tup + tdn) 

         IF (tup < 273.15) THEN
           IF (tm < 273.15) THEN
             IF (tdn < 273.15) THEN
               tavg = (tup + 2.0 * tm + tdn) * 0.25    

             ELSE       !Tup,Tm < 0; Td > 0 
               xtemp = (273.15 - tm) * dzh(k) / (tdn - tm) 

               tavg = 0.5 * (tup * dzh(k) + tm * (dzh(k) + xtemp)    &
                 + 273.15 * (2.0*dzh(k) - xtemp)) / dza(k)     
             END IF

           ELSE

             IF (tdn < 273.15) THEN  ! Tup<0, Td <0; Tm >0
               xup = (273.15 - tup) * dzh(k) / (tm - tup) 
               xdn = dzh(k) - (273.15-tm) * dzh(k) / (tdn - tm) 
               tavg = 0.5 * (tup*xup + 273.15 *                &
                 (2.0*dza(k) - xup -xdn)) / dza(k) 

             ELSE

               xup = (273.15 - tup) * dzh(k) / (tm - tup) 
               tavg = 0.5 * (tup * xup + 273.15* (2.0*          &
                 dza(k) - xup)) / dza(k) 
             END IF
           END IF
         ELSE                                    
           IF (tm < 273.15) THEN 
             IF (tdn < 273.15) THEN !Tup>0; Tm,Tdn < 0 
               xup = dzh(k) - (273.15 - tup) * dzh(k) /         &
                 (tm - tup)
               tavg = 0.5 * (273.15 * (dza(k)-xup) +             &  
                 tm * (dzh(k) + xup) + tdn*dzh(k)) / dza(k)

             ELSE 
  
               xup  = dzh(k) - (273.15- tup) * dzh(k) / (tm - tup)
               xdn  = (273.15 - tm) * dzh(k) / (tdn - tm)
               tavg = 0.5 * (273.15*(2.*dza(k) - xup -xdn)+  &
                 tm * (xup + xdn)) / dza(k)  
             ENDIF
           ELSE
             IF (tdn < 273.15) THEN

               xdn  = dzh(k) - (273.15 - tm) * dzh(k) / (tdn - tm)
               tavg = (273.15 * (dza(k) - xdn) + 0.5*(273.15 & 
                 + tdn)*xdn) / dza(k) 

             ELSE
                tavg = (tup + 2.0*tm + tdn) * 0.25

             ENDIF
           ENDIF
         ENDIF

!-----------------------------------------------------------------------
!    Compute amount of supercooled liquid water content 
!-----------------------------------------------------------------------

         bx = bslope(soiltyp(i,j)) 
         IF (bslope(soiltyp(i,j)) > 5.5) bx = 5.5 

         nlog = 0
         kcount = 0
  
         IF (tavg > (273.15-1.0E-3)) THEN
           frh2o = qsoil(i,j,k) 
         ELSE 

!-----------------------------------------------------------------------
!    Option #1: Iterated solution, ck /= 0
!             See Koren et al., JGR 1999; Eqtn #17
!-----------------------------------------------------------------------

           IF (ck /= 0.0) THEN 
             swl = qsoil(i,j,k) - xqsoil(i,j,k) 
             IF (swl > (qsoil(i,j,k)-0.02)) swl = qsoil(i,j,k)-0.02 
             IF (swl < 0.0) swl = 0.0

             DO   
             IF (.NOT. ((nlog < 10).AND.(kcount == 0))) EXIT     
               nlog = nlog + 1

            df = ALOG( (psisat(soiltyp(i,j))*9.81/3.335E5) *      & 
                  ( (1.0+ck*swl)**2.0) *                          &
                  (wsat(soiltyp(i,j))/(qsoil(i,j,k)-swl))**bx)    &
                  - ALOG( -(tavg - 273.15)/tavg) 

            denom = 2.0 * ck / (1.0+ck*swl) + bx / (qsoil(i,j,k)-swl)
            swlk = swl - df/denom 
  
            IF (swlk > (qsoil(i,j,k)-0.02)) swlk = qsoil(i,j,k)-0.02
            IF (swlk < 0.0) swlk = 0.0

             dswl = ABS(swlk - swl)
             swl = swlk 
        
!     If more than 10 iterations required, then the explicit soltn used.
 
             IF (dswl <= error) kcount = kcount+1  

             END DO  
        
             frh2o = qsoil(i,j,k) - swl 

           END IF 

!--------------------------------------------------------------------------
!    Option #2: Explicit solution for Flerchinger Eqtn., ck = 0
!              See Koren et al., JGR 1999; Eqtn #17 
!--------------------------------------------------------------------------

          IF (kcount == 0) THEN

            fk = (((3.335E5 / (9.81 * (-psisat(soiltyp(i,j))))) *  &
              ((tavg - 273.15)/tavg))                              &
              **(-1.0/bx)) * wsat(soiltyp(i,j)) 

            IF (fk < 0.02) fk = 0.02
            frh2o = MIN(fk, qsoil(i,j,k)) 

          END IF  
!---------------------------------------------------------------------------

        END IF   !End of computing amount of supercooled liquid water 
!---------------------------------------------------------------------------
 
          IF ( (xh2o < xqsoil(i,j,k)) .AND. (xh2o < frh2o)) THEN 
            IF (frh2o > xqsoil(i,j,k)) THEN 
              xh2o = xqsoil(i,j,k)
            ELSE
              xh2o = frh2o 
            END IF 
          END IF 

          IF ( (xh2o > xqsoil(i,j,k)) .AND. (xh2o > frh2o)) THEN 
            IF (frh2o < xqsoil(i,j,k)) THEN
              xh2o = xqsoil(i,j,k) 
            ELSE 
              xh2o = frh2o
            END IF
          END IF 

          IF (xh2o < 0.0) xh2o = 0.0
          IF (xh2o > qsoil(i,j,k)) xh2o = qsoil(i,j,k) 

          sink(i,j,k) = -1000.0 * 3.335E5 * dza(k) * (xh2o - xqsoil(i,j,k)) / &  
            dtsfc 

          xqsoil(i,j,k) = xh2o  
!------------------------------------------------------------------------------


!        tem4soil(i,j,k) = (tempcoefb * tempcoefx1 * tsoil(i,j,k-1) + &
!            (1.0-tempcoefb*(tempcoefx1 + tempcoefx2))*tsoil(i,j,k)     &
!            + tempcoefb * tempcoefx2 *tsoil(i,j,k+1) ) 

         tem4soil(i,j,k) = (tempcoefb * tempcoefx1 * tsoil(i,j,k-1) + &
            (1.0-tempcoefb*(tempcoefx1 + tempcoefx2))*tsoil(i,j,k)     &
            + tempcoefb * tempcoefx2 *tsoil(i,j,k+1) ) +               &
            (sink(i,j,k)*dzsoilinv/cisoil(i,j,k))

!----------------------------------------------------------------------
        END IF

        END DO
      END DO
    END DO 

  RETURN
END SUBROUTINE snksrc


