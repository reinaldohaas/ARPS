!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    SOILCST.INC                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'soilcst.inc' for ARPS soil model.
!
!  This file defines the parameters and common blocks to be used
!  in the soil-vegetation energy and moisture model
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  09/25/1995
!
!  MODIFICATION HISTORY:
!
!  02/25/02 (J. Brotzge)
!  Added several new variables to specify soil moisure/temp
!    profiles as a function of soil type.
!    Data taken from Chen and Dudhia, MWR, 2001;
!    Cosby et al. (1984), Wat Res Res.
!    and Rawls et al (1983), J. of Hydraulic Eng.
!
!-----------------------------------------------------------------------
!
  REAL :: bslope(13)         ! Slope of the retention curve
  DATA bslope/  4.05,  4.38,  4.90,  5.30,  5.39,  7.12,                &
                7.75,  8.52, 10.40, 10.40, 11.40, 11.55, 0.00/

  REAL :: cgsat(13)          ! Thermal coefficient for bare ground
                             ! at saturation
  DATA cgsat/ 3.222E-6, 3.057E-6, 3.560E-6, 4.418E-6,                   &
              4.111E-6, 3.670E-6, 3.593E-6, 3.995E-6,                   &
              3.058E-6, 3.729E-6, 3.600E-6, 1.00E-30, 8.0116E-6/

  REAL :: psisat(13)         ! Saturated matric potential (m)
  DATA psisat/ .121, .090, .218, .786, .478, .299,                    &
               .356, .630, .153, .490, .405, .036, .000/

  REAL :: kns(13)            ! Saturated hydraulic conductivity (m/s)
  DATA kns/ 1.07E-6, 1.41E-5, 5.23E-6, 2.81E-6, 3.38E-6, 4.45E-6,     &
              2.04E-6, 2.45E-6, 7.22E-6, 1.34E-6, 9.74E-7, 1.34E-6,   &
              0.0/

  INTEGER :: pwgeq(13)       ! Exponent in Wgeq formula
  DATA pwgeq/  4,  4,  4,  6,  6,  6,                                   &
               8, 10,  8, 10, 12,  0,  0/

  REAL :: awgeq(13)          ! Coefficient in Wgeq formula
  DATA awgeq/ .387, .404, .219, .105, .148, .135,                       &
              .127, .084, .139, .075, .083, .000, .000/

  REAL :: c1sat(13)          ! Value of c1 at saturation (when d1=10cm)
  DATA c1sat/ .082, .098, .132, .153, .191, .213,                       &
              .385, .227, .421, .375, .342, .000, .000/

  REAL :: c2ref(13)          ! Value of c2 at wetdp = 0.5*wsat
  DATA c2ref/ 3.9, 3.7, 1.8,  .8,  .8,  .8,                             &
               .4,  .6,  .3,  .3,  .3,  .0,  .0/

  REAL :: wsat(13)           ! Saturated volumetric moisture content.(ETA)
  DATA wsat / .421, .421, .434, .434, .439, .404,                       &
              .464, .465, .406, .406, .468, 1.0E-20, 1.00/

  REAL :: wfc(13)            ! Field capacity moisture. (ETA) 
  DATA wfc  / .283, .283, .312, .312, .329, .315,                       &
             .387, .382, .338, .338, .412, 1.0E-25, 1.00/

  REAL :: wwlt(13)           ! Wilting volumetric moisture content. (ETA)
  DATA wwlt / .029, .029, .047, .047, .066, .069,                       &
              .119, .103, .020, .020, .139, 1.0E-30, 1.00/


!  REAL :: wsat(13)           ! Saturated volumetric moisture content.
!  DATA wsat / .395, .410, .435, .485, .451, .420,                       &
!              .477, .476, .426, .492, .482, 1.0E-20, 1.00/
!
!  REAL :: wfc(13)            ! Field capacity moisture.
!  DATA wfc  / .135, .150, .195, .255, .240, .255,                       &
!             .322, .325, .310, .370, .367, 1.0E-25, 1.00/
!
!  REAL :: wwlt(13)           ! Wilting volumetric moisture content.
!  DATA wwlt / .068, .075, .114, .179, .155, .175,                       &
!              .218, .250, .219, .283, .286, 1.0E-30, 1.00/

  REAL :: rgl(14)            ! Species-dependent threshold value of
                             ! incoming solar radiation
  DATA rgl  / 100.0, 100.0, 100.0,  30.0,  30.0,  30.0,  30.0,          &
               30.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0/

  REAL :: rsmin(14)          ! Minimum of surface resistance
  DATA rsmin/ 999.0, 500.0,  40.0,  80.0,  80.0, 100.0, 150.0,          &
              100.0, 999.0,  40.0,   5.0,  80.0, 300.0, 1.0E-30/

  REAL :: alpha(7)    ! Polynomial fit used for estimating des/dT
  DATA alpha/ 610.7799961, 44.36518521,1.428945805,            &
           2.650648471E-2, 3.031240396E-4,2.034080948E-6,      &
           6.136820929E-9/

  REAL :: quartz(13)         !Quartz content
  DATA quartz/ 0.92, 0.82, 0.60, 0.25, 0.40, 0.60, 0.10,       &
               0.35, 0.52, 0.10, 0.25, 0.00, 0.00/

  REAL :: porosity(13)       ! Effective porosity (cm3/cm3)
  DATA porosity/ 0.417, 0.421, 0.434, 0.486, 0.439, 0.404,     &
                 0.464, 0.465, 0.406, 0.423, 0.468, 0.000,     &
                 0.000/


  REAL :: rootzone(14)       ! Effective root zone depth (m) 
  DATA rootzone/ 1.0, 0.40, 1.0, 1.0, 1.0, 2.0, 2.0,           &
                 2.0, 0.40, 1.0, 1.0, 1.0, 1.0, 2.0/     

  REAL :: snup(14)           ! Maximum snow depth threshold at max albedo
  DATA snup/ 0.040, 0.025, 0.040, 0.040, 0.040, 0.080, 0.080,  &
             0.080, 0.025, 0.040, 0.040, 0.040, 0.025, 0.000/ 

  REAL :: hsf2(14)           ! Variable used for F2 in Rc estimate 
  DATA hsf2/ 0.00, 42.0, 36.35, 42.0, 54.53, 54.53, 41.69,     &
             60.0, 999.0, 36.35, 60.0, 42.0, 200.0, 0.0/   

  REAL :: cgv                ! Thermal coefficient for vegetation
  DATA cgv/1.5E-05/

  REAL :: dirfc              ! Directional factor of Rayleigh
                             ! scattering
                             ! and absorption transmission function,
                             ! above cloud: dirfc = dirf (not for sfc.)
  DATA dirfc/1.6666667/   ! below & within cloud: dirfc = 5/3

  REAL :: d1                 ! Top layer depth of soil column: 10 cm
  DATA d1/ 0.10/          ! Chaning d1 requires changing C1sat. The
                          ! latter is propotional to the inverse of d1

  REAL :: d2                 ! Deep layer depth of soil column: 1 m
  DATA d2/ 1.0/           ! Deep layer depth of soil column: 1 m

  REAL :: rsmax              ! Maximum surface resistance in s/m
  DATA rsmax/5000./       ! Maximum surface resistence

  REAL :: wetsml             ! NP, Eq. 18, small value to limit
  DATA wetsml/1.0E-30/    ! Small value to limit c2wg at saturation
                          ! w2 = wsat

!  DATA roufns/ 0.002, 0.020, 0.010, 0.100, 0.150, 0.750, 1.000,
!    :             2.000, 0.005, 0.010, 0.020, 0.060, 0.040, 0.001/

!  DATA veg  / 0.10, 0.10, 0.60, 0.40, 0.40, 0.90, 0.99,
!    :            0.99, 0.01, 0.30, 0.99, 0.40, 0.20, 0.00/

  REAL :: snow_albedo        ! Albedo of snow covered ground (was .65, now .6)
  REAL :: snowdepth_crit     ! Depth (m) at which snow cover is considered to
                             ! to be 100%.
  REAL :: snowflxfac         ! Factor by which to modify snow fluxes
  REAL :: cg_snow            ! heat capacity of snow
  PARAMETER (snow_albedo=0.60,snowdepth_crit=0.025,                     &
             cg_snow=1.0/60610.0,snowflxfac=1.0)

