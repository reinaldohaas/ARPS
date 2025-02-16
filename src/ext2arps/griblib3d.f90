!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE GRIBTBLS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gribtbls
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Load the tables for encoding/decoding GRIB Edition 2 parameters.
!  These are used to translate cryptic GRIB message numbers into
!  understandable meteorological descriptions of variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Craig A. Mattocks
!  6/1/1995.
!
!  MODIFICATION HISTORY:
!
!  6/1/1995, version 1.0 (C. Mattocks)
!    Wrote original subroutine.
!
!  12/14/1995 (Yuhe Liu)
!  Changed the subroutine name from TABLES to GRIBTBLS.
!
!  Removed arguments and moved the table declarations to include
!  file gribcst.inc.
!
!  Used some original undefined variable parameters to define ARPS
!  variables.
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!  The following arrays are loaded as per the GRIB tables and
!  carried by a COMMON block, grib001, defined in gribcst.inc:
!
!  Center       Names of oringinating meteorological centers
!  SbCenter     Names of oringinating meteorological sub-centers
!  Model        Names of numerical models/generating processes
!  VarName      Names         of  GRIB variables
!  VarUnit     Units         of  GRIB variables
!  VarAbbr      Abbreviations for GRIB variables
!  VarLvl       Levels        of  GRIB variables
!  TimeUnit     Forecast time units
!  TimeRang     Time range indicators
!  Proj         Names of map projections
!  Scan         Directional scanning modes
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: n                        ! Loop index
  CHARACTER (LEN=2) :: endstr               ! End-of-string marker
  DATA endstr/'  '/
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'gribcst.inc'
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
!  Originating meteorological operational center names:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  center (  1) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center (  2) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center (  3) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center (  4) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center (  5) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center (  6) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center (  7) = 'US NWS - National Meteorological Center' // endstr
  center (  8) = 'US NWS - Telecomms Gateway'              // endstr
  center (  9) = 'US NWS - Field Stations'                 // endstr

  DO n = 10, 33
    center (n) = 'UNDEFINED METEOROLOGICAL CENTER'        // endstr
  END DO

  center ( 34) = 'Japanese Meteorological Agency - Tokyo'  // endstr

  DO n = 35, 51
    center (n) = 'UNDEFINED METEOROLOGICAL CENTER'        // endstr
  END DO

  center ( 52) = 'National Hurricane Center/TPC, Miami'    // endstr
  center ( 53) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center ( 54) = 'Canadian Met Service - Montreal'         // endstr
  center ( 55) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center ( 56) = 'UNDEFINED METEOROLOGICAL CENTER'         // endstr
  center ( 57) = 'U.S. Air Force - Global Weather Center'  // endstr
  center ( 58) = 'US Navy - Fleet Numerical Ocean Center'  // endstr
  center ( 59) = 'NOAA Forecast Systems Lab, Boulder CO'   // endstr

  center ( 60) = 'NCAR, Boulder CO'                        // endstr
  DO n = 61, 73
    center (n) = 'UNDEFINED METEOROLOGICAL CENTER'        // endstr
  END DO

  center ( 74) = 'U.K. Met Office - Bracknell'             // endstr

  DO n = 75, 84
    center (n) = 'UNDEFINED METEOROLOGICAL CENTER'        // endstr
  END DO

  center ( 85) = 'French Weather Service - Toulouse'       // endstr

  DO n = 86, 96
    center (n) = 'UNDEFINED METEOROLOGICAL CENTER'        // endstr
  END DO

  center ( 97) = 'European Space Agency (ESA)'             // endstr
  center ( 98) = 'Eur Ctr for Medium-Range Weather Forcsts'// endstr
  center ( 99) = 'DeBilt, Netherlands'                     // endstr
!
!-----------------------------------------------------------------------
!
!  Originating meteorological operational sub-center names:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  sbcenter(  0)= 'RESERVED METEOROLOGICAL SUB-CENTER'      // endstr
  sbcenter(  1)= 'NMC Re-Analysis Project'                 // endstr
  DO n = 2, 149
    sbcenter (n) = 'UNDEFINED METEOROLOGICAL SUB-CENTER'  // endstr
  END DO

  sbcenter(150) = 'Arkansas-Red River RFC, Tulsa OK'       // endstr
  sbcenter(151) = 'Alaska RFC, Anchorage, AK '             // endstr
  sbcenter(152) = 'Colorado Basin RFC, Salt Lake City, UT' // endstr
  sbcenter(153) = 'California-Nevada RFC, Sacramento, CA'  // endstr
  sbcenter(154) = 'Lower Mississippi RFC, Slidel, LA'      // endstr
  sbcenter(155) = 'Middle Atlantic RFC, State College, PA' // endstr
  sbcenter(156) = 'Missouri Basin RFC, Kansas City, MO'    // endstr
  sbcenter(157) = 'North Central RFC, Minneapolis, MN'     // endstr
  sbcenter(158) = 'Northeast RFC, Hartford, CT'            // endstr
  sbcenter(159) = 'Northwest RFC, Portland, OR'            // endstr

  sbcenter(160) = 'Ohio Basin RFC, Cincinnati, OH'         // endstr
  sbcenter(161) = 'Southeast RFC, Atlanta, GA'             // endstr
  sbcenter(162) = 'West Gulf RFC, Fort Worth, TX'          // endstr

  DO n = 163, 169
    sbcenter (n) = 'UNDEFINED METEOROLOGICAL SUB-CENTER'  // endstr
  END DO

  sbcenter(170) = 'Norman, OK WFO'                         // endstr
!
!-----------------------------------------------------------------------
!
!  Names of numerical model or generating process:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  model  (  1) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  (  2) = 'Ultra-Violet Potential Index Model'      // endstr
  model  (  3) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  (  4) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  (  5) = 'IR Satellite-Derived Precip and Temp'    // endstr
  model  (  6) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  (  7) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  (  8) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  (  9) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr

  model  ( 10) = 'Global Wind-Wave Forecast Model'         // endstr
  DO n = 11, 18
    model (n) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  END DO
  model  ( 19) = 'Limited-area Fine Mesh (LFM) analysis'   // endstr

  DO n = 20, 24
    model (n) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  END DO
  model  ( 25) = 'Snow Cover Analysis'                     // endstr
  DO n = 26, 38
    model (n) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  END DO

  model  ( 39) = 'Nested Grid forecast Model (NGM)'        // endstr

  model  ( 40) = 'ECMWF 213-wave triang, 31-level Spectral'// endstr
  model  ( 41) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 42) = 'Global OI Analysis from Aviation run'    // endstr
  model  ( 43) = 'Global OI Analysis from Final run'       // endstr
  model  ( 44) = 'Sea Surface Temperature (SST) Analysis'  // endstr
  model  ( 45) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 46) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 47) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 48) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 49) = 'Ozone Analysis - TIROS Sat Observations' // endstr

  model  ( 50) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 51) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 52) = 'Ozone Analysis - Nimbus 7 Observations'  // endstr
  model  ( 53) = 'LFM Fourth-Order Forecast Model'         // endstr
  DO n = 54, 63
    model (n) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  END DO

  model  ( 64) = 'Regional Optimum Interpolation Analysis' // endstr
  model  ( 65) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 66) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 67) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 68) = '80-wave triang, 18-layer Spectral AVN'   // endstr
  model  ( 69) = '80-wave triang, 18 layer Spectral MRF'   // endstr

  model  ( 70) = 'Quasi-Lagrangian Hurricane Model (QLM)'  // endstr
  model  ( 71) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 72) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  model  ( 73) = 'Fog Forecast model - Ocean Prod. Center' // endstr
  model  ( 74) = 'Gulf of Mexico Wind/Wave'                // endstr
  model  ( 75) = 'Gulf of Alaska Wind/Wave'                // endstr
  model  ( 76) = 'Bias-corrected Medium Range Forecast'    // endstr
  model  ( 77) = '126-wave triang, 28-layer Spectral AVN'  // endstr
  model  ( 78) = '126-wave triang, 28-layer Spectral MRF'  // endstr
  model  ( 79) = 'Backup from the previous run'            // endstr

  model  ( 80) = '62-wave triang, 18-layer Spectral MRF'   // endstr
  model  ( 81) = 'Spectral Statistical Interpolation AVN'  // endstr
  model  ( 82) = 'Spectral Statistical Interpolation FINAL'// endstr
  model  ( 83) = 'ETA Model - 80 km version'               // endstr
  model  ( 84) = 'ETA Model - 40 km version'               // endstr
  model  ( 85) = 'ETA Model - 30 km version'               // endstr
  model  ( 86) = 'RUC/MAPS Hybrid Sigma/Theta Model, FSL'  // endstr
  model  ( 87) = 'CAC Ensemble Forecasts - Spectral Model' // endstr
  model  ( 88) = 'Ocean Wave Model with add Physics (PWAV)'// endstr

  DO n = 89, 149
    model (n) = 'UNDEFINED NUMERICAL MODEL/PROCESS'       // endstr
  END DO

  model  (150) = 'NWS River Forecast System (NWSRFS)'      // endstr
  model  (151) = 'NWS Flash Flood Guidance Systm (NWSFFGS)'// endstr
  model  (152) = 'WSR-88D Stage II Precipitation Analysis' // endstr
  model  (153) = 'WSR-88D Stage III Precipitation Analysis'// endstr
!
!-----------------------------------------------------------------------
!
!  GRIB variable names:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  varname(  0) = 'RESERVED'                                // endstr
  varname(  1) = 'Pressure'                                // endstr
  varname(  2) = 'Pressure reduced to MSL'                 // endstr
  varname(  3) = 'Pressure tendency'                       // endstr
  varname(  4) = 'UNDEFINED VARIABLE'                      // endstr
  varname(  5) = 'UNDEFINED VARIABLE'                      // endstr
  varname(  6) = 'Geopotential'                            // endstr
  varname(  7) = 'Geopotential height'                     // endstr
  varname(  8) = 'Geometric height'                        // endstr
  varname(  9) = 'Standard deviation of height'            // endstr

  varname( 10) = 'UNDEFINED VARIABLE'                      // endstr
  varname( 11) = 'Temperature'                             // endstr
  varname( 12) = 'Virtual temperature'                     // endstr
  varname( 13) = 'Potential temperature'                   // endstr
  varname( 14) = 'Equivalent potential temperature'        // endstr
  varname( 15) = 'Maximum temperature'                     // endstr
  varname( 16) = 'Minimum temperature'                     // endstr
  varname( 17) = 'Dew point temperature'                   // endstr
  varname( 18) = 'Dew point depression'                    // endstr
  varname( 19) = 'Lapse rate'                              // endstr

  varname( 20) = 'Visibility'                              // endstr
  varname( 21) = 'Radar Spectra (1)'                       // endstr
  varname( 22) = 'Radar Spectra (2)'                       // endstr
  varname( 23) = 'Radar Spectra (3)'                       // endstr
  varname( 24) = 'Total Ozone'                             // endstr
  varname( 25) = 'Temperature anomaly'                     // endstr
  varname( 26) = 'Pressure anomaly'                        // endstr
  varname( 27) = 'Geopotential height anomaly'             // endstr
  varname( 28) = 'Wave Spectra (1)'                        // endstr
  varname( 29) = 'Wave Spectra (2)'                        // endstr

  varname( 30) = 'Wave Spectra (3)'                        // endstr
  varname( 31) = 'Wind direction'                          // endstr
  varname( 32) = 'Wind speed'                              // endstr
  varname( 33) = 'u-component of wind'                     // endstr
  varname( 34) = 'v-component of wind'                     // endstr
  varname( 35) = 'Stream function'                         // endstr
  varname( 36) = 'Velocity potential'                      // endstr
  varname( 37) = 'Montgomery stream function'              // endstr
  varname( 38) = 'Sigma coordinate vertical velocity'      // endstr
  varname( 39) = 'Pressure vertical velocity'              // endstr

  varname( 40) = 'Geometric vertical velocity'             // endstr
  varname( 41) = 'Absolute vorticity'                      // endstr
  varname( 42) = 'Absolute divergence'                     // endstr
  varname( 43) = 'Relative vorticity'                      // endstr
  varname( 44) = 'Relative divergence'                     // endstr
  varname( 45) = 'Vertical u-component shear'              // endstr
  varname( 46) = 'Vertical v-component shear'              // endstr
  varname( 47) = 'Direction of current'                    // endstr
  varname( 48) = 'Speed of current'                        // endstr
  varname( 49) = 'u-component of current'                  // endstr

  varname( 50) = 'v-component of current'                  // endstr
  varname( 51) = 'Specific humidity'                       // endstr
  varname( 52) = 'Relative humidity'                       // endstr
  varname( 53) = 'Humidity mixing ratio'                   // endstr
  varname( 54) = 'Precipitable water'                      // endstr
  varname( 55) = 'Vapor pressure'                          // endstr
  varname( 56) = 'Saturation deficit'                      // endstr
  varname( 57) = 'Evaporation'                             // endstr
  varname( 58) = 'Cloud Ice'                               // endstr
  varname( 59) = 'Precipitation rate'                      // endstr

  varname( 60) = 'Thunderstorm probability'                // endstr
  varname( 61) = 'Total precipitation'                     // endstr
  varname( 62) = 'Large scale precipitation (non-conv)'    // endstr
  varname( 63) = 'Convective precipitation'                // endstr
  varname( 64) = 'Snowfall rate water equivalent'          // endstr
  varname( 65) = 'Water equiv of accum snow depth'         // endstr
  varname( 66) = 'Snow depth'                              // endstr
  varname( 67) = 'Mixed layer depth'                       // endstr
  varname( 68) = 'Transient thermocline depth'             // endstr
  varname( 69) = 'Main thermocline depth'                  // endstr

  varname( 70) = 'Main thermocline anomaly'                // endstr
  varname( 71) = 'Total cloud cover'                       // endstr
  varname( 72) = 'Convective cloud cover'                  // endstr
  varname( 73) = 'Low cloud cover'                         // endstr
  varname( 74) = 'Medium cloud cover'                      // endstr
  varname( 75) = 'High cloud cover'                        // endstr
  varname( 76) = 'Cloud water'                             // endstr
  varname( 77) = 'Condensation pressure'                   // endstr
  varname( 78) = 'Convective snow'                         // endstr
  varname( 79) = 'Large scale snow'                        // endstr

  varname( 80) = 'Water Temperature'                       // endstr
  varname( 81) = 'Land-sea mask (land=1;sea=0)'            // endstr
  varname( 82) = 'Deviation of sea level from mean'        // endstr
  varname( 83) = 'Surface roughness'                       // endstr
  varname( 84) = 'Albedo'                                  // endstr
  varname( 85) = 'Soil temperature'                        // endstr
  varname( 86) = 'Soil moisture content'                   // endstr
  varname( 87) = 'Vegetation'                              // endstr
  varname( 88) = 'Salinity'                                // endstr
  varname( 89) = 'Density'                                 // endstr

  varname( 90) = 'Water runoff'                            // endstr
  varname( 91) = 'Ice concentration (ice=1;no ice=0)'      // endstr
  varname( 92) = 'Ice thickness'                           // endstr
  varname( 93) = 'Direction of ice drift'                  // endstr
  varname( 94) = 'Speed of ice drift'                      // endstr
  varname( 95) = 'u-component of ice drift'                // endstr
  varname( 96) = 'v-component of ice drift'                // endstr
  varname( 97) = 'Ice growth rate'                         // endstr
  varname( 98) = 'Ice divergence'                          // endstr
  varname( 99) = 'Snow melt'                               // endstr

  varname(100) = 'Significant height of wind waves + swell'// endstr
  varname(101) = 'Direction of wind waves'                 // endstr
  varname(102) = 'Significant height of wind waves'        // endstr
  varname(103) = 'Mean period of wind waves'               // endstr
  varname(104) = 'Direction of swell waves'                // endstr
  varname(105) = 'Significant height of swell waves'       // endstr
  varname(106) = 'Mean period of swell waves'              // endstr
  varname(107) = 'Primary wave direction'                  // endstr
  varname(108) = 'Primary wave mean period'                // endstr
  varname(109) = 'Secondary wave direction'                // endstr

  varname(110) = 'Secondary wave mean period'              // endstr
  varname(111) = 'Net short-wave radiation (surface)'      // endstr
  varname(112) = 'Net long-wave  radiation (surface)'      // endstr
  varname(113) = 'Net short-wave radiation (top of atm)'   // endstr
  varname(114) = 'Net long-wave  radiation (top of atm)'   // endstr
  varname(115) = 'Long  wave radiation'                    // endstr
  varname(116) = 'Short wave radiation'                    // endstr
  varname(117) = 'Global radiation'                        // endstr
  varname(118) = 'UNDEFINED VARIABLE'                      // endstr
  varname(119) = 'UNDEFINED VARIABLE'                      // endstr

  varname(120) = 'UNDEFINED VARIABLE'                      // endstr
  varname(121) = 'Latent heat net flux'                    // endstr
  varname(122) = 'Sensible heat net flux'                  // endstr
  varname(123) = 'Boundary layer dissipation'              // endstr
  varname(124) = 'Momentum flux, u component'              // endstr
  varname(125) = 'Momentum flux, v component'              // endstr
  varname(126) = 'Wind mixing energy'                      // endstr
  varname(127) = 'Image data'                              // endstr
  varname(128) = 'Mean Sea Level Pressure (STD atm red)'   // endstr
  varname(129) = 'Mean Sea Level Pressure (MAPS atm red)'  // endstr

  varname(130) = 'Mean Sea Level Pressure (ETA model red)' // endstr
  varname(131) = 'Surface lifted index'                    // endstr
  varname(132) = 'Best (4 layer) lifted index'             // endstr
  varname(133) = 'K index'                                 // endstr
  varname(134) = 'Sweat index'                             // endstr
  varname(135) = 'Horizontal moisture divergence'          // endstr
  varname(136) = 'Vertical speed shear'                    // endstr
  varname(137) = '3-hr pressure tendency (STD atm red)'    // endstr
  varname(138) = 'Brunt-Vaisala frequency (squared)'       // endstr
  varname(139) = 'Potential vorticity  (density weighted)' // endstr

  varname(140) = 'Categorical rain (yes=1; no=0)'          // endstr
  varname(141) = 'Categorical freezing rain (yes=1; no=0)' // endstr
  varname(142) = 'Categorical ice pellets (yes=1; no=0)'   // endstr
  varname(143) = 'Categorical snow (yes=1; no=0)'          // endstr
  varname(144) = 'Volumetric soil moisture content'        // endstr
  varname(145) = 'Potential evaporation rate'              // endstr
  varname(146) = 'Cloud work function'                     // endstr
  varname(147) = 'Zonal flux of gravity wave stress'       // endstr
  varname(148) = 'Meridonal flux of gravity wave stress'   // endstr
  varname(149) = 'Potential vorticity'                     // endstr

  varname(150) = 'UV covariance'                           // endstr
  varname(151) = 'UT covariance'                           // endstr
  varname(152) = 'VT covariance'                           // endstr
  varname(153) = 'Mass Flux Divergence'                    // endstr
  varname(154) = 'UNDEFINED VARIABLE'                      // endstr
  varname(155) = 'Ground Heat Flux'                        // endstr
  varname(156) = 'Convective inhibition'                   // endstr
  varname(157) = 'Convective Available Potential Energy'   // endstr
  varname(158) = 'Turbulent Kinetic Energy'                // endstr
  varname(159) = 'Condensation pressure of lifted parcel'  // endstr

  varname(160) = 'Clear Sky Upward Solar Flux'             // endstr
  varname(161) = 'Clear Sky Downward Solar Flux'           // endstr
  varname(162) = 'Clear Sky upward long wave flux'         // endstr
  varname(163) = 'Clear Sky downward long wave flux'       // endstr
  varname(164) = 'Cloud forcing net solar flux'            // endstr
  varname(165) = 'Cloud forcing net long wave flux'        // endstr
  varname(166) = 'Visible beam downward solar flux'        // endstr
  varname(167) = 'Visible diffuse downward solar flux'     // endstr
  varname(168) = 'Near IR beam downward solar flux'        // endstr
  varname(169) = 'Near IR diffuse downward solar flux'     // endstr

  varname(170) = 'u-component of friction velocity'        // endstr
  varname(171) = 'v-component of friction velocity'        // endstr
  varname(172) = 'Momentum flux'                           // endstr
  varname(173) = 'Mass point model surface'                // endstr
  varname(174) = 'Velocity point model surface'            // endstr
  varname(175) = 'Model layer number (from bottom up)'     // endstr
  varname(176) = 'Latitude north (-90 to +90)'             // endstr
  varname(177) = 'Longitude east (0 to 360)'               // endstr
  varname(178) = 'UNDEFINED VARIABLE'                      // endstr
  varname(179) = 'UNDEFINED VARIABLE'                      // endstr

  varname(180) = 'UNDEFINED VARIABLE'                      // endstr
  varname(181) = 'x-gradient of log pressure'              // endstr
  varname(182) = 'y-gradient of log pressure'              // endstr
  varname(183) = 'x-gradient of height'                    // endstr
  varname(184) = 'y-gradient of height'                    // endstr
!
!-----------------------------------------------------------------------
!
!  The following are defined for ARPS model. They were not defined
!  in the original WMO/NMC table 2
!
!-----------------------------------------------------------------------
!
  varname(185) = 'Horizontal turbulence Mixing coefficient'// endstr
  varname(186) = 'Hail mixing ratio'                       // endstr
  varname(187) = 'Soil type'                               // endstr
  varname(188) = 'Vegetation type'                         // endstr
  varname(189) = 'Leaf Area Index'                         // endstr

  varname(190) = 'Volumetric soil moisture'                // endstr
  varname(191) = 'Amount of liquid water on canopy'        // endstr
  varname(192) = 'Vertical turbulence Mixing coefficient'  // endstr
  varname(193) = 'UNDEFINED VARIABLE'                      // endstr
  varname(194) = 'UNDEFINED VARIABLE'                      // endstr
  varname(195) = 'UNDEFINED VARIABLE'                      // endstr
  varname(196) = 'UNDEFINED VARIABLE'                      // endstr
  varname(197) = 'UNDEFINED VARIABLE'                      // endstr
  varname(198) = 'UNDEFINED VARIABLE'                      // endstr
  varname(199) = 'UNDEFINED VARIABLE'                      // endstr

  varname(200) = 'UNDEFINED VARIABLE'                      // endstr
  varname(201) = 'Ice-free water surface'                  // endstr
  varname(202) = 'UNDEFINED VARIABLE'                      // endstr
  varname(203) = 'UNDEFINED VARIABLE'                      // endstr
  varname(204) = 'Downward short wave radiation flux'      // endstr
  varname(205) = 'Downward long wave radiation flux'       // endstr
  varname(206) = 'Ultra-violet index (1 hr @ solar noon)'  // endstr
  varname(207) = 'Moisture availability'                   // endstr
  varname(208) = 'Exchange coefficient'                    // endstr
  varname(209) = 'Number of mixed layers next to surface'  // endstr

  varname(210) = 'UNDEFINED VARIABLE'                      // endstr
  varname(211) = 'Upward short wave radiation flux'        // endstr
  varname(212) = 'Upward long wave radiation flux'         // endstr
  varname(213) = 'Amount of non-convective cloud'          // endstr
  varname(214) = 'Convective Precipitation rate'           // endstr
  varname(215) = 'Temperature tendency by all physics'     // endstr
  varname(216) = 'Temperature tendency by all radiation'   // endstr
  varname(217) = 'Temperature tendency by non-rad physics' // endstr
  varname(218) = 'Precipitation index (0.00-1.00)'         // endstr
  varname(219) = 'Std. dev. of IR T over 1x1 deg area'     // endstr

  varname(220) = 'Natural logarithm of surface pressure'   // endstr
  varname(221) = 'UNDEFINED VARIABLE'                      // endstr
  varname(222) = '5-wave geopotential height'              // endstr
  varname(223) = 'Plant canopy surface water'              // endstr
  varname(224) = 'UNDEFINED VARIABLE'                      // endstr
  varname(225) = 'UNDEFINED VARIABLE'                      // endstr
  varname(226) = 'Blackadar mixing length scale'           // endstr
  varname(227) = 'Asymptotic mixing length scale'          // endstr
  varname(228) = 'Potential evaporation'                   // endstr
  varname(229) = 'Snow phase-change heat flux'             // endstr

  varname(230) = 'UNDEFINED VARIABLE'                      // endstr
  varname(231) = 'Convective cloud mass flux'              // endstr
  varname(232) = 'Downward total radiation flux'           // endstr
  varname(233) = 'Upward total radiation flux'             // endstr
  varname(234) = 'Baseflow-groundwater runoff'             // endstr
  varname(235) = 'Storm surface runoff'                    // endstr
  varname(236) = 'UNDEFINED VARIABLE'                      // endstr
  varname(237) = 'UNDEFINED VARIABLE'                      // endstr
  varname(238) = 'Snow cover'                              // endstr
  varname(239) = 'Snow temperature'                        // endstr

  varname(240) = 'UNDEFINED VARIABLE'                      // endstr
  varname(241) = 'Large scale condensational heating rate' // endstr
  varname(242) = 'Deep convective heating rate'            // endstr
  varname(243) = 'Deep convective moistening rate'         // endstr
  varname(244) = 'Shallow convective heating rate'         // endstr
  varname(245) = 'Shallow convective moistening rate'      // endstr
  varname(246) = 'Vertical diffusion heating rate'         // endstr
  varname(247) = 'Vertical diffusion zonal acceleration'   // endstr
  varname(248) = 'Vertical diffusion merid acceleration'   // endstr
  varname(249) = 'Vertical diffusion moistening rate'      // endstr

  varname(250) = 'Solar radiative heating rate'            // endstr
  varname(251) = 'Long wave radiative heating rate'        // endstr
  varname(252) = 'Drag coefficient'                        // endstr
  varname(253) = 'Friction velocity'                       // endstr
  varname(254) = 'Richardson number'                       // endstr
  varname(255) = 'MISSING VARIABLE'                        // endstr
!
!-----------------------------------------------------------------------
!
!  GRIB variable units:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER           ----+----+----+----+----+----+----+----+
  varunit(  0) = 'RESERVED UNITS'                          // endstr
  varunit(  1) = 'Pa'                                      // endstr
  varunit(  2) = 'Pa'                                      // endstr
  varunit(  3) = 'Pa/s'                                    // endstr
  varunit(  4) = 'UNDEFINED UNITS'                         // endstr
  varunit(  5) = 'UNDEFINED UNITS'                         // endstr
  varunit(  6) = 'm2/s2'                                   // endstr
  varunit(  7) = 'gpm'                                     // endstr
  varunit(  8) = 'm'                                       // endstr
  varunit(  9) = 'm'                                       // endstr

  varunit( 10) = 'UNDEFINED UNITS'                         // endstr
  varunit( 11) = 'K'                                       // endstr
  varunit( 12) = 'K'                                       // endstr
  varunit( 13) = 'K'                                       // endstr
  varunit( 14) = 'K'                                       // endstr
  varunit( 15) = 'K'                                       // endstr
  varunit( 16) = 'K'                                       // endstr
  varunit( 17) = 'K'                                       // endstr
  varunit( 18) = 'K'                                       // endstr
  varunit( 19) = 'K/m'                                     // endstr

  varunit( 20) = 'm'                                       // endstr
  varunit( 21) = '-'                                       // endstr
  varunit( 22) = '-'                                       // endstr
  varunit( 23) = '-'                                       // endstr
  varunit( 24) = 'Dobson units'                            // endstr
  varunit( 25) = 'K'                                       // endstr
  varunit( 26) = 'Pa'                                      // endstr
  varunit( 27) = 'gpm'                                     // endstr
  varunit( 28) = '-'                                       // endstr
  varunit( 29) = '-'                                       // endstr

  varunit( 30) = '-'                                       // endstr
  varunit( 31) = 'deg'                                     // endstr
  varunit( 32) = 'm/s'                                     // endstr
  varunit( 33) = 'm/s'                                     // endstr
  varunit( 34) = 'm/s'                                     // endstr
  varunit( 35) = 'm2/s'                                    // endstr
  varunit( 36) = 'm2/s'                                    // endstr
  varunit( 37) = 'm2/s2'                                   // endstr
  varunit( 38) = '1/s'                                     // endstr
  varunit( 39) = 'Pa/s'                                    // endstr

  varunit( 40) = 'm/s'                                     // endstr
  varunit( 41) = '1/s'                                     // endstr
  varunit( 42) = '1/s'                                     // endstr
  varunit( 43) = '1/s'                                     // endstr
  varunit( 44) = '1/s'                                     // endstr
  varunit( 45) = '1/s'                                     // endstr
  varunit( 46) = '1/s'                                     // endstr
  varunit( 47) = 'deg'                                     // endstr
  varunit( 48) = 'm/s'                                     // endstr
  varunit( 49) = 'm/s'                                     // endstr

  varunit( 50) = 'm/s'                                     // endstr
  varunit( 51) = 'kg/kg'                                   // endstr
  varunit( 52) = '%'                                       // endstr
  varunit( 53) = 'kg/kg'                                   // endstr
  varunit( 54) = 'kg/m2'                                   // endstr
  varunit( 55) = 'Pa'                                      // endstr
  varunit( 56) = 'Pa'                                      // endstr
  varunit( 57) = 'kg/m2'                                   // endstr
  varunit( 58) = 'kg/m2'                                   // endstr
  varunit( 59) = 'kg/m2/s'                                 // endstr

  varunit( 60) = '%'                                       // endstr
  varunit( 61) = 'kg/m2'                                   // endstr
  varunit( 62) = 'kg/m2'                                   // endstr
  varunit( 63) = 'kg/m2'                                   // endstr
  varunit( 64) = 'kg/m2/s'                                 // endstr
  varunit( 65) = 'kg/m2'                                   // endstr
  varunit( 66) = 'm'                                       // endstr
  varunit( 67) = 'm'                                       // endstr
  varunit( 68) = 'm'                                       // endstr
  varunit( 69) = 'm'                                       // endstr

  varunit( 70) = 'm'                                       // endstr
  varunit( 71) = '%'                                       // endstr
  varunit( 72) = '%'                                       // endstr
  varunit( 73) = '%'                                       // endstr
  varunit( 74) = '%'                                       // endstr
  varunit( 75) = '%'                                       // endstr
  varunit( 76) = 'kg/m2'                                   // endstr
  varunit( 77) = 'Pa'                                      // endstr
  varunit( 78) = 'kg/m2'                                   // endstr
  varunit( 79) = 'kg/m2'                                   // endstr

  varunit( 80) = 'K'                                       // endstr
  varunit( 81) = 'fraction'                                // endstr
  varunit( 82) = 'm'                                       // endstr
  varunit( 83) = 'm'                                       // endstr
  varunit( 84) = '%'                                       // endstr
  varunit( 85) = 'K'                                       // endstr
  varunit( 86) = 'kg/m2'                                   // endstr
  varunit( 87) = '%'                                       // endstr
  varunit( 88) = 'kg/kg'                                   // endstr
  varunit( 89) = 'kg/m3'                                   // endstr

  varunit( 90) = 'kg/m2'                                   // endstr
  varunit( 91) = 'fraction'                                // endstr
  varunit( 92) = 'm'                                       // endstr
  varunit( 93) = 'deg true'                                // endstr
  varunit( 94) = 'm/s'                                     // endstr
  varunit( 95) = 'm/s'                                     // endstr
  varunit( 96) = 'm/s'                                     // endstr
  varunit( 97) = 'm/s'                                     // endstr
  varunit( 98) = '1/s'                                     // endstr
  varunit( 99) = 'kg/m2'                                   // endstr

  varunit(100) = 'm'                                       // endstr
  varunit(101) = 'deg true'                                // endstr
  varunit(102) = 'm'                                       // endstr
  varunit(103) = 's'                                       // endstr
  varunit(104) = 'deg true'                                // endstr
  varunit(105) = 'm'                                       // endstr
  varunit(106) = 's'                                       // endstr
  varunit(107) = 'deg true'                                // endstr
  varunit(108) = 's'                                       // endstr
  varunit(109) = 'deg true'                                // endstr

  varunit(110) = 's'                                       // endstr
  varunit(111) = 'W/m2'                                    // endstr
  varunit(112) = 'W/m2'                                    // endstr
  varunit(113) = 'W/m2'                                    // endstr
  varunit(114) = 'W/m2'                                    // endstr
  varunit(115) = 'W/m2'                                    // endstr
  varunit(116) = 'W/m2'                                    // endstr
  varunit(117) = 'W/m2'                                    // endstr
  varunit(118) = 'UNDEFINED UNITS'                         // endstr
  varunit(119) = 'UNDEFINED UNITS'                         // endstr

  varunit(120) = 'UNDEFINED UNITS'                         // endstr
  varunit(121) = 'W/m2'                                    // endstr
  varunit(122) = 'W/m2'                                    // endstr
  varunit(123) = 'W/m2'                                    // endstr
  varunit(124) = 'N/m2'                                    // endstr
  varunit(125) = 'N/m2'                                    // endstr
  varunit(126) = 'J'                                       // endstr
  varunit(127) = '-'                                       // endstr
  varunit(128) = 'Pa'                                      // endstr
  varunit(129) = 'Pa'                                      // endstr

  varunit(130) = 'Pa'                                      // endstr
  varunit(131) = 'K'                                       // endstr
  varunit(132) = 'K'                                       // endstr
  varunit(133) = 'K'                                       // endstr
  varunit(134) = 'K'                                       // endstr
  varunit(135) = 'kg/kg/s'                                 // endstr
  varunit(136) = '1/s'                                     // endstr
  varunit(137) = 'Pa/s'                                    // endstr
  varunit(138) = '1/s2'                                    // endstr
  varunit(139) = '1/s/m'                                   // endstr

  varunit(140) = 'non-dim'                                 // endstr
  varunit(141) = 'non-dim'                                 // endstr
  varunit(142) = 'non-dim'                                 // endstr
  varunit(143) = 'non-dim'                                 // endstr
  varunit(144) = 'fraction'                                // endstr
  varunit(145) = 'W/m2'                                    // endstr
  varunit(146) = 'J/kg'                                    // endstr
  varunit(147) = 'N/m2'                                    // endstr
  varunit(148) = 'N/m2'                                    // endstr
  varunit(149) = 'm2/s/kg'                                 // endstr

  varunit(150) = 'm2/s2'                                   // endstr
  varunit(151) = 'K*m/s'                                   // endstr
  varunit(152) = 'K*m/s'                                   // endstr
  varunit(153) = 'UNDEFINED UNITS'                         // endstr
  varunit(154) = 'UNDEFINED UNITS'                         // endstr
  varunit(155) = 'W/m2'                                    // endstr
  varunit(156) = 'J/kg'                                    // endstr
  varunit(157) = 'J/kg'                                    // endstr
  varunit(158) = 'J/kg'                                    // endstr
  varunit(159) = 'Pa'                                      // endstr

  varunit(160) = 'W/m2'                                    // endstr
  varunit(161) = 'W/m2'                                    // endstr
  varunit(162) = 'W/m2'                                    // endstr
  varunit(163) = 'W/m2'                                    // endstr
  varunit(164) = 'W/m2'                                    // endstr
  varunit(165) = 'W/m2'                                    // endstr
  varunit(166) = 'W/m2'                                    // endstr
  varunit(167) = 'W/m2'                                    // endstr
  varunit(168) = 'W/m2'                                    // endstr
  varunit(169) = 'W/m2'                                    // endstr

  varunit(170) = 'm/s'                                     // endstr
  varunit(171) = 'm/s'                                     // endstr
  varunit(172) = 'N/m2'                                    // endstr
  varunit(173) = 'non-dim'                                 // endstr
  varunit(174) = 'non-dim'                                 // endstr
  varunit(175) = 'non-dim'                                 // endstr
  varunit(176) = 'deg'                                     // endstr
  varunit(177) = 'deg'                                     // endstr
  varunit(178) = 'UNDEFINED UNITS'                         // endstr
  varunit(179) = 'UNDEFINED UNITS'                         // endstr

  varunit(180) = 'UNDEFINED UNITS'                         // endstr
  varunit(181) = '1/m'                                     // endstr
  varunit(182) = '1/m'                                     // endstr
  varunit(183) = 'm/m'                                     // endstr
  varunit(184) = 'm/m'                                     // endstr
  varunit(185) = 'm**2/s'                                  // endstr
  varunit(186) = '(kg/kg)'                                 // endstr
  varunit(187) = ' '                                       // endstr
  varunit(188) = ' '                                       // endstr
  varunit(189) = ' '                                       // endstr

  varunit(190) = 'm**3/m**3'                               // endstr
  varunit(191) = 'm'                                       // endstr
  varunit(192) = 'UNDEFINED UNITS'                         // endstr
  varunit(193) = 'UNDEFINED UNITS'                         // endstr
  varunit(194) = 'UNDEFINED UNITS'                         // endstr
  varunit(195) = 'UNDEFINED UNITS'                         // endstr
  varunit(196) = 'UNDEFINED UNITS'                         // endstr
  varunit(197) = 'UNDEFINED UNITS'                         // endstr
  varunit(198) = 'UNDEFINED UNITS'                         // endstr
  varunit(199) = 'UNDEFINED UNITS'                         // endstr

  varunit(200) = 'UNDEFINED UNITS'                         // endstr
  varunit(201) = '%'                                       // endstr
  varunit(202) = 'UNDEFINED UNITS'                         // endstr
  varunit(203) = 'UNDEFINED UNITS'                         // endstr
  varunit(204) = 'W/m2'                                    // endstr
  varunit(205) = 'W/m2'                                    // endstr
  varunit(206) = 'J/m2'                                    // endstr
  varunit(207) = '%'                                       // endstr
  varunit(208) = '(kg/m3)(m/s)'                            // endstr
  varunit(209) = 'integer'                                 // endstr

  varunit(210) = 'UNDEFINED UNITS'                         // endstr
  varunit(211) = 'W/m2'                                    // endstr
  varunit(212) = 'W/m2'                                    // endstr
  varunit(213) = '%'                                       // endstr
  varunit(214) = 'kg/m2/s'                                 // endstr
  varunit(215) = 'K/s'                                     // endstr
  varunit(216) = 'K/s'                                     // endstr
  varunit(217) = 'K/s'                                     // endstr
  varunit(218) = 'fraction'                                // endstr
  varunit(219) = 'K'                                       // endstr

  varunit(220) = 'Ln(kPa)'                                 // endstr
  varunit(221) = 'UNDEFINED UNITS'                         // endstr
  varunit(222) = 'gpm'                                     // endstr
  varunit(223) = 'kg/m2'                                   // endstr
  varunit(224) = 'UNDEFINED UNITS'                         // endstr
  varunit(225) = 'UNDEFINED UNITS'                         // endstr
  varunit(226) = 'm'                                       // endstr
  varunit(227) = 'm'                                       // endstr
  varunit(228) = 'kg/m2'                                   // endstr
  varunit(229) = 'W/m2'                                    // endstr

  varunit(230) = 'UNDEFINED UNITS'                         // endstr
  varunit(231) = 'Pa/s'                                    // endstr
  varunit(232) = 'W/m2'                                    // endstr
  varunit(233) = 'W/m2'                                    // endstr
  varunit(234) = 'kg/m2'                                   // endstr
  varunit(235) = 'kg/m2'                                   // endstr
  varunit(236) = 'UNDEFINED UNITS'                         // endstr
  varunit(237) = 'UNDEFINED UNITS'                         // endstr
  varunit(238) = 'percent'                                 // endstr
  varunit(239) = 'K'                                       // endstr

  varunit(240) = 'UNDEFINED UNITS'                         // endstr
  varunit(241) = 'K/s'                                     // endstr
  varunit(242) = 'K/s'                                     // endstr
  varunit(243) = 'kg/kg/s'                                 // endstr
  varunit(244) = 'K/s'                                     // endstr
  varunit(245) = 'kg/kg/s'                                 // endstr
  varunit(246) = 'K/s'                                     // endstr
  varunit(247) = 'm/s2'                                    // endstr
  varunit(248) = 'm/s2'                                    // endstr
  varunit(249) = 'kg/kg/s'                                 // endstr

  varunit(250) = 'K/s'                                     // endstr
  varunit(251) = 'K/s'                                     // endstr
  varunit(252) = 'non-dim'                                 // endstr
  varunit(253) = 'm/s'                                     // endstr
  varunit(254) = 'non-dim'                                 // endstr
  varunit(255) = 'MISSING UNITS'                           // endstr
!
!-----------------------------------------------------------------------
!
!  GRIB variable abbreviations:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER           ----+----+----+----+----+----+----+----+
  varabbr(  0) = 'RESERVED ABBREVIATION'                  // endstr
  varabbr(  1) = 'PRES'                                   // endstr
  varabbr(  2) = 'PRMSL'                                  // endstr
  varabbr(  3) = 'PTEND'                                  // endstr
  varabbr(  4) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(  5) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(  6) = 'GP'                                     // endstr
  varabbr(  7) = 'HGT'                                    // endstr
  varabbr(  8) = 'DIST'                                   // endstr
  varabbr(  9) = 'HSTDV'                                  // endstr

  varabbr( 10) = 'HVAR'                                   // endstr
  varabbr( 11) = 'TMP'                                    // endstr
  varabbr( 12) = 'VTMP'                                   // endstr
  varabbr( 13) = 'POT'                                    // endstr
  varabbr( 14) = 'EPOT'                                   // endstr
  varabbr( 15) = 'T MAX'                                  // endstr
  varabbr( 16) = 'T MIN'                                  // endstr
  varabbr( 17) = 'DPT'                                    // endstr
  varabbr( 18) = 'DEPR'                                   // endstr
  varabbr( 19) = 'LAPR'                                   // endstr

  varabbr( 20) = 'VIS'                                    // endstr
  varabbr( 21) = 'RDSP1'                                  // endstr
  varabbr( 22) = 'RDSP2'                                  // endstr
  varabbr( 23) = 'RDSP3'                                  // endstr
  varabbr( 24) = 'TOTO3'                                  // endstr
  varabbr( 25) = 'TMP A'                                  // endstr
  varabbr( 26) = 'PRESA'                                  // endstr
  varabbr( 27) = 'GP A'                                   // endstr
  varabbr( 28) = 'WVSP1'                                  // endstr
  varabbr( 29) = 'WVSP2'                                  // endstr

  varabbr( 30) = 'WVSP3'                                  // endstr
  varabbr( 31) = 'WDIR'                                   // endstr
  varabbr( 32) = 'WIND'                                   // endstr
  varabbr( 33) = 'U GRD'                                  // endstr
  varabbr( 34) = 'V GRD'                                  // endstr
  varabbr( 35) = 'STRM'                                   // endstr
  varabbr( 36) = 'V POT'                                  // endstr
  varabbr( 37) = 'MNTSF'                                  // endstr
  varabbr( 38) = 'SGCVV'                                  // endstr
  varabbr( 39) = 'V VEL'                                  // endstr

  varabbr( 40) = 'DZDT'                                   // endstr
  varabbr( 41) = 'ABS V'                                  // endstr
  varabbr( 42) = 'ABS D'                                  // endstr
  varabbr( 43) = 'REL V'                                  // endstr
  varabbr( 44) = 'REL D'                                  // endstr
  varabbr( 45) = 'VUCSH'                                  // endstr
  varabbr( 46) = 'VVCSH'                                  // endstr
  varabbr( 47) = 'DIR C'                                  // endstr
  varabbr( 48) = 'SP C'                                   // endstr
  varabbr( 49) = 'UOGRD'                                  // endstr

  varabbr( 50) = 'VOGRD'                                  // endstr
  varabbr( 51) = 'SPF H'                                  // endstr
  varabbr( 52) = 'R H'                                    // endstr
  varabbr( 53) = 'MIXR'                                   // endstr
  varabbr( 54) = 'P WAT'                                  // endstr
  varabbr( 55) = 'VAPP'                                   // endstr
  varabbr( 56) = 'SAT D'                                  // endstr
  varabbr( 57) = 'EVP'                                    // endstr
  varabbr( 58) = 'C ICE'                                  // endstr
  varabbr( 59) = 'PRATE'                                  // endstr

  varabbr( 60) = 'TSTM'                                   // endstr
  varabbr( 61) = 'A PCP'                                  // endstr
  varabbr( 62) = 'NCPCP'                                  // endstr
  varabbr( 63) = 'ACPCP'                                  // endstr
  varabbr( 64) = 'SRWEQ'                                  // endstr
  varabbr( 65) = 'WEASD'                                  // endstr
  varabbr( 66) = 'SNO D'                                  // endstr
  varabbr( 67) = 'MIXHT'                                  // endstr
  varabbr( 68) = 'TTHDP'                                  // endstr
  varabbr( 69) = 'MTHD'                                   // endstr

  varabbr( 70) = 'MTH A'                                  // endstr
  varabbr( 71) = 'T CDC'                                  // endstr
  varabbr( 72) = 'CDCON'                                  // endstr
  varabbr( 73) = 'L CDC'                                  // endstr
  varabbr( 74) = 'M CDC'                                  // endstr
  varabbr( 75) = 'H CDC'                                  // endstr
  varabbr( 76) = 'C WAT'                                  // endstr
  varabbr( 77) = 'CNDP'                                   // endstr
  varabbr( 78) = 'SNO C'                                  // endstr
  varabbr( 79) = 'SNO L'                                  // endstr

  varabbr( 80) = 'WTMP'                                   // endstr
  varabbr( 81) = 'LAND'                                   // endstr
  varabbr( 82) = 'DSL M'                                  // endstr
  varabbr( 83) = 'SFC R'                                  // endstr
  varabbr( 84) = 'ALBDO'                                  // endstr
  varabbr( 85) = 'TSOIL'                                  // endstr
  varabbr( 86) = 'SOIL M'                                 // endstr
  varabbr( 87) = 'VEG'                                    // endstr
  varabbr( 88) = 'SALTY'                                  // endstr
  varabbr( 89) = 'DEN'                                    // endstr

  varabbr( 90) = 'WATR'                                   // endstr
  varabbr( 91) = 'ICE C'                                  // endstr
  varabbr( 92) = 'ICETK'                                  // endstr
  varabbr( 93) = 'DICED'                                  // endstr
  varabbr( 94) = 'SICED'                                  // endstr
  varabbr( 95) = 'U ICE'                                  // endstr
  varabbr( 96) = 'V ICE'                                  // endstr
  varabbr( 97) = 'ICE G'                                  // endstr
  varabbr( 98) = 'ICE D'                                  // endstr
  varabbr( 99) = 'SNO M'                                  // endstr

  varabbr(100) = 'HTSGW'                                  // endstr
  varabbr(101) = 'WVDIR'                                  // endstr
  varabbr(102) = 'WVHGT'                                  // endstr
  varabbr(103) = 'WVPER'                                  // endstr
  varabbr(104) = 'SWDIR'                                  // endstr
  varabbr(105) = 'SWELL'                                  // endstr
  varabbr(106) = 'SWPER'                                  // endstr
  varabbr(107) = 'DIRPW'                                  // endstr
  varabbr(108) = 'PERPW'                                  // endstr
  varabbr(109) = 'DIRSW'                                  // endstr

  varabbr(110) = 'PERSW'                                  // endstr
  varabbr(111) = 'NSWRS'                                  // endstr
  varabbr(112) = 'NLWRS'                                  // endstr
  varabbr(113) = 'NSWRT'                                  // endstr
  varabbr(114) = 'NLWRT'                                  // endstr
  varabbr(115) = 'LWAVR'                                  // endstr
  varabbr(116) = 'SWAVR'                                  // endstr
  varabbr(117) = 'G RAD'                                  // endstr
  varabbr(118) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(119) = 'UNDEFINED ABBREVIATION'                 // endstr

  varabbr(120) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(121) = 'LHTFL'                                  // endstr
  varabbr(122) = 'SHTFL'                                  // endstr
  varabbr(123) = 'BLYDP'                                  // endstr
  varabbr(124) = 'U FLX'                                  // endstr
  varabbr(125) = 'V FLX'                                  // endstr
  varabbr(126) = 'WMIXE'                                  // endstr
  varabbr(127) = 'IMG D'                                  // endstr
  varabbr(128) = 'MSLSA'                                  // endstr
  varabbr(129) = 'MSLMA'                                  // endstr

  varabbr(130) = 'MSLET'                                  // endstr
  varabbr(131) = 'LFT X'                                  // endstr
  varabbr(132) = '4LFTX'                                  // endstr
  varabbr(133) = 'K X'                                    // endstr
  varabbr(134) = 'S X'                                    // endstr
  varabbr(135) = 'MCONV'                                  // endstr
  varabbr(136) = 'VW SH'                                  // endstr
  varabbr(137) = 'TSLSA'                                  // endstr
  varabbr(138) = 'BVF 2'                                  // endstr
  varabbr(139) = 'PV MW'                                  // endstr

  varabbr(140) = 'CRAIN'                                  // endstr
  varabbr(141) = 'CFRZRN'                                 // endstr
  varabbr(142) = 'CICEPL'                                 // endstr
  varabbr(143) = 'CSNOW'                                  // endstr
  varabbr(144) = 'SOILW'                                  // endstr
  varabbr(145) = 'PEVPR'                                  // endstr
  varabbr(146) = 'CWORK'                                  // endstr
  varabbr(147) = 'U-GWD'                                  // endstr
  varabbr(148) = 'V-GWD'                                  // endstr
  varabbr(149) = 'PV'                                     // endstr

  varabbr(150) = 'COVMZ'                                  // endstr
  varabbr(151) = 'COVTZ'                                  // endstr
  varabbr(152) = 'COVTM'                                  // endstr
  varabbr(153) = 'MFXDV'                                  // endstr
  varabbr(154) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(155) = 'GFLUX'                                  // endstr
  varabbr(156) = 'CIN'                                    // endstr
  varabbr(157) = 'CAPE'                                   // endstr
  varabbr(158) = 'TKE'                                    // endstr
  varabbr(159) = 'CONDP'                                  // endstr

  varabbr(160) = 'CSUSF'                                  // endstr
  varabbr(161) = 'CSDSF'                                  // endstr
  varabbr(162) = 'CSULF'                                  // endstr
  varabbr(163) = 'CSDLF'                                  // endstr
  varabbr(164) = 'CFNSF'                                  // endstr
  varabbr(165) = 'CFNLF'                                  // endstr
  varabbr(166) = 'VBDSF'                                  // endstr
  varabbr(167) = 'VDDSF'                                  // endstr
  varabbr(168) = 'NBDSF'                                  // endstr
  varabbr(169) = 'NDDSF'                                  // endstr

  varabbr(170) = 'USTR'                                   // endstr
  varabbr(171) = 'VSTR'                                   // endstr
  varabbr(172) = 'M FLX'                                  // endstr
  varabbr(173) = 'LMH'                                    // endstr
  varabbr(174) = 'LMV'                                    // endstr
  varabbr(175) = 'MLYNO'                                  // endstr
  varabbr(176) = 'NLAT'                                   // endstr
  varabbr(177) = 'ELON'                                   // endstr
  varabbr(178) = 'UMAS'                                   // endstr
  varabbr(179) = 'VMAS'                                   // endstr

  varabbr(180) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(181) = 'LPS X'                                  // endstr
  varabbr(182) = 'LPS Y'                                  // endstr
  varabbr(183) = 'HGT X'                                  // endstr
  varabbr(184) = 'HGT Y'                                  // endstr
  varabbr(185) = 'STDZ'                                   // endstr
  varabbr(186) = 'STDU'                                   // endstr
  varabbr(187) = 'STDV'                                   // endstr
  varabbr(188) = 'STDQ'                                   // endstr
  varabbr(189) = 'STDT'                                   // endstr

  varabbr(190) = 'CBUW'                                   // endstr
  varabbr(191) = 'CBVW'                                   // endstr
  varabbr(192) = 'CBUQ'                                   // endstr
  varabbr(193) = 'CBVQ'                                   // endstr
  varabbr(194) = 'CBTW'                                   // endstr
  varabbr(195) = 'CBQW'                                   // endstr
  varabbr(196) = 'CBMZW'                                  // endstr
  varabbr(197) = 'CBTZW'                                  // endstr
  varabbr(198) = 'CBTMW'                                  // endstr
  varabbr(199) = 'STDRH'                                  // endstr

  varabbr(200) = 'SDTZ'                                   // endstr
  varabbr(201) = 'ICWAT'                                  // endstr
  varabbr(202) = 'SDTU'                                   // endstr
  varabbr(203) = 'SDTV'                                   // endstr
  varabbr(204) = 'DSWRF'                                  // endstr
  varabbr(205) = 'DLWRF'                                  // endstr
  varabbr(206) = 'UVI'                                    // endstr
  varabbr(207) = 'MSTAV'                                  // endstr
  varabbr(208) = 'SFEXC'                                  // endstr
  varabbr(209) = 'MIXLY'                                  // endstr

  varabbr(210) = 'SDTT'                                   // endstr
  varabbr(211) = 'USWRF'                                  // endstr
  varabbr(212) = 'ULWRF'                                  // endstr
  varabbr(213) = 'CDLYR'                                  // endstr
  varabbr(214) = 'CPRAT'                                  // endstr
  varabbr(215) = 'TTDIA'                                  // endstr
  varabbr(216) = 'TTRAD'                                  // endstr
  varabbr(217) = 'TTPHY'                                  // endstr
  varabbr(218) = 'PREIX'                                  // endstr
  varabbr(219) = 'TSD1D'                                  // endstr

  varabbr(220) = 'NLGSP'                                  // endstr
  varabbr(221) = 'SDTRH'                                  // endstr
  varabbr(222) = '5WAVH'                                  // endstr
  varabbr(223) = 'C WAT'                                  // endstr
  varabbr(224) = 'PLTRS'                                  // endstr
  varabbr(225) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(226) = 'BMIXL'                                  // endstr
  varabbr(227) = 'AMIXL'                                  // endstr
  varabbr(228) = 'PEVAP'                                  // endstr
  varabbr(229) = 'SNOHF'                                  // endstr

  varabbr(230) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(231) = 'MFLUX'                                  // endstr
  varabbr(232) = 'DTRF'                                   // endstr
  varabbr(233) = 'UTRF'                                   // endstr
  varabbr(234) = 'BGRUN'                                  // endstr
  varabbr(235) = 'SSRUN'                                  // endstr
  varabbr(236) = 'UNDEFINED ABBREVIATION'                 // endstr
  varabbr(237) = 'OZONE'                                  // endstr
  varabbr(238) = 'SNO C'                                  // endstr
  varabbr(239) = 'SNO T'                                  // endstr

  varabbr(240) = 'GLCR'                                  // endstr
  varabbr(241) = 'LRGHR'                                  // endstr
  varabbr(242) = 'CNVHR'                                  // endstr
  varabbr(243) = 'CNVMR'                                  // endstr
  varabbr(244) = 'SHAHR'                                  // endstr
  varabbr(245) = 'SHAMR'                                  // endstr
  varabbr(246) = 'VDFHR'                                  // endstr
  varabbr(247) = 'VDFUA'                                  // endstr
  varabbr(248) = 'VDFVA'                                  // endstr
  varabbr(249) = 'VDFMR'                                  // endstr

  varabbr(250) = 'SWHR'                                   // endstr
  varabbr(251) = 'LWHR'                                   // endstr
  varabbr(252) = 'CD'                                     // endstr
  varabbr(253) = 'FRICV'                                  // endstr
  varabbr(254) = 'RI'                                     // endstr
  varabbr(255) = 'MISSING ABBREVIATION'                   // endstr
!
!-----------------------------------------------------------------------
!
!  GRIB variable levels/layers:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  varlvl (  0) = 'RESERVED LEVEL'                         // endstr
  varlvl (  1) = 'surface of the Earth'                   // endstr
  varlvl (  2) = 'cloud base level'                       // endstr
  varlvl (  3) = 'cloud top level'                        // endstr
  varlvl (  4) = 'deepest 0 deg (C) isotherm level'       // endstr
  varlvl (  5) = 'adiabatic lifted condensation level'    // endstr
  varlvl (  6) = 'maximum wind speed level'               // endstr
  varlvl (  7) = 'tropopause'                             // endstr
  varlvl (  8) = 'nominal top of atmosphere'              // endstr
  varlvl (  9) = 'sea bottom'                             // endstr

  varlvl ( 10) = 'throughout atmospheric column'          // endstr
  varlvl ( 11) = 'UNDEFINED LEVEL'                        // endstr
  varlvl ( 12) = 'low cloud base level'                   // endstr
  varlvl ( 13) = 'low cloud top level'                    // endstr
  varlvl ( 14) = 'within low cloud layer'                 // endstr
  DO n = 15, 21
    varlvl (n) = 'RESERVED LEVEL'                        // endstr
  END DO

  varlvl ( 22) = 'middle cloud base level'                // endstr
  varlvl ( 23) = 'middle cloud top level'                 // endstr
  varlvl ( 24) = 'within middle cloud layer'              // endstr
  DO n = 25, 31
    varlvl (n) = 'RESERVED LEVEL'                        // endstr
  END DO

  varlvl ( 32) = 'high cloud base level'                  // endstr
  varlvl ( 33) = 'high cloud top level'                   // endstr
  varlvl ( 34) = 'within high cloud layer'                // endstr
  DO n = 35, 99
    varlvl (n) = 'RESERVED LEVEL'                        // endstr
  END DO

  varlvl (100) = 'millibars (hPa) isobaric level'         // endstr
  varlvl (101) = 'kPa @ top:bottom of pressure layer'     // endstr
  varlvl (102) = 'mean sea level (MSL)'                   // endstr
  varlvl (103) = 'meters above mean sea level (MSL)'      // endstr
  varlvl (104) = 'meters above MSL @ top:bottom of Z layer'// endstr
  varlvl (105) = 'meters above ground level (AGL)'        // endstr
  varlvl (106) = 'meters AGL @ top:bottom of Z layer'     // endstr
  varlvl (107) = 'sigma level value (1/10000)'            // endstr
  varlvl (108) = 'top:bottom of sigma layer (1/100)'      // endstr
  varlvl (109) = 'hybrid level number'                    // endstr

  varlvl (110) = 'top:bottom of hybrid layer (level #)'   // endstr
  varlvl (111) = 'centimeters below land surface'         // endstr
  varlvl (112) = 'centimeters BGL @ top:bottom of Z layer'// endstr
  varlvl (113) = 'degrees K isentropic level'             // endstr
  varlvl (114) = 'top:bot isentropic layer (475-theta K)' // endstr
  DO n = 115, 120
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (121) = 'hPa @ top:bot of pressure layer (1100-p)'// endstr
  varlvl (122) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (123) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (124) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (125) = 'centimeters above ground level (AGL)'   // endstr
  varlvl (126) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (127) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (128) = 'top:bot sigma layer (1.1-sigma, 1/1000)'// endstr

  DO n = 129, 140
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (141) = 'top:bot pressure layer (kPa, 1100-p hPa)'// endstr

  DO n = 142, 159
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (160) = 'meters below mean sea level (MSL)'      // endstr

  DO n = 161, 199
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (200) = 'throughout atmospheric column'          // endstr
  varlvl (201) = 'throughout water column'                // endstr
  varlvl (202) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (203) = 'UNDEFINED LEVEL'                        // endstr
  varlvl (204) = 'highest tropospheric freezing level'    // endstr

  DO n = 205, 211
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (212) = 'low cloud base level'                   // endstr
  varlvl (213) = 'low cloud top level'                    // endstr
  varlvl (214) = 'within low cloud layer'                 // endstr

  DO n = 215, 221
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (222) = 'middle cloud base level'                // endstr
  varlvl (223) = 'middle cloud top level'                 // endstr
  varlvl (224) = 'within middle cloud layer'              // endstr

  DO n = 225, 231
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO

  varlvl (232) = 'high cloud base level'                  // endstr
  varlvl (233) = 'high cloud top level'                   // endstr
  varlvl (234) = 'within high cloud layer'                // endstr

  DO n = 235, 255
    varlvl (n) = 'UNDEFINED LEVEL'                       // endstr
  END DO
!
!-----------------------------------------------------------------------
!
!  Forecast time units:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  timeunit(  0)= 'minutes'                                // endstr
  timeunit(  1)= 'hours'                                  // endstr
  timeunit(  2)= 'days'                                   // endstr
  timeunit(  3)= 'months'                                 // endstr
  timeunit(  4)= 'years'                                  // endstr
  timeunit(  5)= 'decades'                                // endstr
  timeunit(  6)= 'normals (30 years)'                     // endstr
  timeunit(  7)= 'centuries'                              // endstr

  DO n = 8, 253
    timeunit(n)= 'RESERVED TIME UNIT'                    // endstr
  END DO

  timeunit(254)= 'seconds'                                // endstr
!
!-----------------------------------------------------------------------
!
!  Time range indicators:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  timerang(  0)= 'Product valid @ ref time + P1 (P1 > 0)' // endstr
  timerang(  1)= 'Analysis product valid @ ref time (P1=0)'// endstr
  timerang(  2)= 'Product valid @ ref time + (P1,...,P2)' // endstr
  timerang(  3)= 'Average of ref time + (P1,...,P2)'      // endstr
  timerang(  4)= 'Accumulation of ref time + (P1,...,P2)' // endstr
  timerang(  5)= 'Difference: (reftime+P2) - (reftime+P1)'// endstr

  DO n = 6, 9
    timeunit(n)='RESERVED TIME RANGE INDICATOR'          // endstr
  END DO

  timerang( 10)= 'Product valid @ reference time + P1'    // endstr

  DO n = 11, 50
    timeunit(n)='RESERVED TIME RANGE INDICATOR'          // endstr
  END DO

  timerang( 51)= 'Climatological mean value: R,...,R+P2'  // endstr

  DO n = 52, 112
    timeunit(n)='RESERVED TIME RANGE INDICATOR'          // endstr
  END DO

  timerang(113)= 'Avg of N forecasts: period=P1, intrvl=P2'// endstr
  timerang(114)= 'Sum of N forecasts: period=P1, intrvl=P2'// endstr
  timerang(115)= 'Avg of N forecasts: P1+P2+...+(N-1)P2 @R'// endstr
  timerang(116)= 'Sum of N forecasts: P1+P2+...+(N-1)P2 @R'// endstr
  timerang(117)= 'Avg of N: P1+(P1-P2)+...+(N-1)(P1-P2)'  // endstr
  timerang(118)= 'Covariance of N analyses:R,...,R+(N-1)P2'// endstr
  timerang(119)= 'RESERVED TIME RANGE INDICATOR'          // endstr

  timerang(120)= 'RESERVED TIME RANGE INDICATOR'          // endstr
  timerang(121)= 'RESERVED TIME RANGE INDICATOR'          // endstr
  timerang(122)= 'RESERVED TIME RANGE INDICATOR'          // endstr
  timerang(123)= 'Avg of N analyses: R,R+P2+...,R+(N-1)P2'// endstr
  timerang(124)= 'Sum of N analyses: R,R+P2+...,R+(N-1)P2'// endstr

  DO n = 125, 254
    timeunit(n)='RESERVED TIME RANGE INDICATOR'          // endstr
  END DO
!
!-----------------------------------------------------------------------
!
!  Names of map projections:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  projs  (  0) = 'Latitude/Longitude Grid'                // endstr
  projs  (  1) = 'Mercator Projection Grid'               // endstr
  projs  (  2) = 'Gnomonic Projection Grid'               // endstr
  projs  (  3) = 'Lambert Conformal Projection Grid'      // endstr
  projs  (  4) = 'Gaussian Latitude/Longitude Grid'       // endstr
  projs  (  5) = 'Polar Stereographic Projection Grid'    // endstr
  projs  (  6) = 'RESERVED PROJECTION'                    // endstr
  projs  (  7) = 'RESERVED PROJECTION'                    // endstr
  projs  (  8) = 'RESERVED PROJECTION'                    // endstr
  projs  (  9) = 'RESERVED PROJECTION'                    // endstr

  projs  ( 10) = 'Latitude/Longitude Grid'                // endstr
  projs  ( 11) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 12) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 13) = 'Oblique Lambert Conformal Projection'   // endstr
  projs  ( 14) = 'Gaussian Latitude/Longitude Grid'       // endstr
  projs  ( 15) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 16) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 17) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 18) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 19) = 'RESERVED PROJECTION'                    // endstr

  projs  ( 20) = 'Latitude/Longitude Grid'                // endstr
  projs  ( 21) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 22) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 23) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 24) = 'Gaussian Latitude/Longitude Grid'       // endstr
  projs  ( 25) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 26) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 27) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 28) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 29) = 'RESERVED PROJECTION'                    // endstr

  projs  ( 30) = 'Latitude/Longitude Grid'                // endstr
  projs  ( 31) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 32) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 33) = 'RESERVED PROJECTION'                    // endstr
  projs  ( 34) = 'Gaussian Latitude/Longitude Grid'       // endstr
  DO n = 35, 49
    projs( n) = 'RESERVED PROJECTION'                    // endstr
  END DO

  projs  ( 50) = 'Spherical Harmonic Coefficients'        // endstr
  DO n = 51, 59
    projs( n) = 'RESERVED PROJECTION'                    // endstr
  END DO

  projs  ( 60) = 'Spherical Harmonic Coefficients'        // endstr
  DO n = 61, 69
    projs( n) = 'RESERVED PROJECTION'                    // endstr
  END DO

  projs  ( 70) = 'Spherical Harmonic Coefficients'        // endstr
  DO n = 71, 79
    projs( n) = 'RESERVED PROJECTION'                    // endstr
  END DO

  projs  ( 80) = 'Spherical Harmonic Coefficients'        // endstr
  DO n = 81, 89
    projs( n) = 'RESERVED PROJECTION'                    // endstr
  END DO

  projs  ( 90) = 'Space view perspective/orthographic grid'// endstr
  DO n = 91, 191
    projs( n) = 'RESERVED PROJECTION'                    // endstr
  END DO
  DO n = 192, 200
    projs( n) = 'PROJECTION AVAILABLE - Consult NWS/NMC.'// endstr
  END DO

  projs  (201) = 'Arakawa semi-staggerd E-grid rot lat,lon'// endstr
  projs  (202) = 'Arakawa filled E-grid rotated lat,lon'  // endstr
  DO n = 203, 254
    projs( n) = 'PROJECTION AVAILABLE - Consult NWS/NMC.'// endstr
  END DO
!
!-----------------------------------------------------------------------
!
!  Description of directional scanning modes:
!
!-----------------------------------------------------------------------
!
!  CHARACTER STRING        10        20        30        40
!  RULER:          ----+----+----+----+----+----+----+----+
  scans  (0,1) = 'Grid points scan in +X direction       '// endstr
  scans  (1,1) = 'Grid points scan in -X direction       '// endstr

  scans  (0,2) = 'Grid points scan in -Y direction       '// endstr
  scans  (1,2) = 'Grid points scan in +Y direction       '// endstr

  scans  (0,3) = 'Consecutive points are in i direction  '// endstr
  scans  (1,3) = 'Consecutive points are in j direction  '// endstr

  scans  (0,4) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (1,4) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (0,5) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (1,5) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (0,6) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (1,6) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (0,7) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (1,7) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (0,8) = 'RESERVED SCANNING MODE                 '// endstr
  scans  (1,8) = 'RESERVED SCANNING MODE                 '// endstr

  RETURN
END SUBROUTINE gribtbls
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GRBGRID                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grbgrid(gridtyp, gdsflg, igds,                               &
           gridesc, iproj, gthin,                                       &
           ni,nj,np, nk,zk, npeq,nit,                                   &
           pi,pj,ipole, di,dj,                                          &
           latsw,lonsw, latne,lonne,                                    &
           latrot,lonrot,angrot,                                        &
           latstr,lonstr,facstr,                                        &
           lattru1,lattru2,lontrue,                                     &
           scanmode, iscan,jscan,kscan,                                 &
           ires,iearth,icomp,                                           &
           jpenta,kpenta,mpenta,ispect,icoeff,                          &
           xp,yp, xo,yo,zo)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set grid parameters according to GRIB Table B or as defined in
!  Section 2 - the Grid Description Section (GDS).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Craig A. Mattocks
!  6/1/1995.
!
!  MODIFICATION HISTORY:
!
!  6/1/1995, version 1.0 (C. Mattocks)
!    Wrote original subroutine.
!
!  01/05/1996 (Yuhe Liu)
!    Modified for NMC GRIB standard and ARPS ext2arps
!
!WDT
!  2003-08-12 (Richard Carpenter)
!    Work around bug in GDS(28) of NCAR/NCEP Reanalysis GRIB files.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  gridtyp  Grid identifier
!  gdsflg   Flag for GDS
!  igds     Integer array for Section 2 of GRIB message, GDS
!
!  OUTPUT:
!
!  gridesc  Grid description string
!
!  Gthin    Indicator of whether the grid is "thinned"
!           (quasi-regular)
!
!  iProj    Projection indicator - WMO Section 2 (GDS) octet #7
!
!  Ni       Number of points along x-axis
!  Nj       Number of points along y-axis
!  Np       Total number of horizontal grid points
!
!  Nk       Number of vertical coordinate parameters
!  Zk       Vertical coordinate parameters
!
!  Npeq     Number of latitude circles from pole to equator
!  Nit      Number of points along x-axis for thinned grid
!
!  Pi       x-coordinate of pole point
!  Pj       y-coordinate of pole point
!  iPole    Projection center flag
!
!  Di       x-direction increment or grid length
!  Dj       y-direction increment or grid length
!
!  LatSW    Latitude  of South West corner point
!  LonSW    Longitude of South West corner point
!  LatNE    Latitude  of North East corner point
!  LonNE    Longitude of North East corner point
!
!  LatRot   Latitude  of southern pole of rotation
!  LonRot   Longitude of southern pole of rotation
!  AngRot   Angle of rotation
!
!  LatStr   Latitude  of the pole of stretching
!  LonStr   Longitude of the pole of stretching
!  FacStr   Stretching factor
!
!  LatTru1  Latitude (1st) at which map projection is "true"
!  LatTru2  Latitude (2nd) at which map projection is "true"
!  LonTrue  Longitude      at which map projection is "true"
!
!  ScanMode Scanning indicator
!  iScan    x-direction   scanning indicator
!  jScan    y-direction   scanning indicator
!  kScan    FORTRAN index scanning indicator
!
!  iRes     Resolution direction increments indicator
!  iEarth   Earth shape indicator: spherical or oblate?
!  iComp    (u,v) components decomposition indicator
!
!  jPenta   J-Pentagonal resolution parameter
!  kPenta   K-Pentagonal resolution parameter
!  mPenta   M-Pentagonal resolution parameter
!
!  iSpect   Spectral representation type
!  iCoeff   Spectral coefficient storage mode
!
!  Xp       X coordinate of sub-satellite point
!  Yp       Y coordinate of sub-satellite point
!  Xo       X coordinate of image sector origin
!  Yo       Y coordinate of image sector origin
!  Zo       Camera altitude from center of Earth
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  GRIB section parameters:
!
!-----------------------------------------------------------------------
!
  INTEGER :: gridtyp       ! Grid indentifier
  INTEGER :: gdsflg        ! GDS indicator
  INTEGER :: igds(*)       ! Integer array for GDS

  CHARACTER (LEN=42) :: gridesc  ! Grid description string

  INTEGER :: gthin         ! Indicator of whether the grid is "thinned"

  INTEGER :: iproj         ! Map projection indicator

  INTEGER :: ni        ! Number of points along x-axis
  INTEGER :: nj        ! Number of points along y-axis
  INTEGER :: np        ! Total number of horizontal grid points

  INTEGER :: nk        ! Number of vertical parameters
  REAL :: zk(*)        ! Vertical coordinate parameters

  INTEGER :: npeq      ! Number of lat circles from pole to equator
  INTEGER :: nit(*)    ! Number of x-points for thinned grid

  REAL :: pi           ! x-coordinate of pole point
  REAL :: pj           ! y-coordinate of pole point
  INTEGER :: ipole     ! Projection center flag

  REAL :: di           ! x-direction increment or grid length
  REAL :: dj           ! y-direction increment or grid length

  REAL :: latsw        ! Latitude  of South West corner point
  REAL :: lonsw        ! Longitude of South West corner point
  REAL :: latne        ! Latitude  of North East corner point
  REAL :: lonne        ! Longitude of North East corner point

  REAL :: lattru1      ! Latitude (1st) at which projection is true
  REAL :: lattru2      ! Latitude (2nd) at which projection is true
  REAL :: lontrue      ! Longitude      at which projection is true

  REAL :: latrot       ! Latitude  of southern pole of rotation
  REAL :: lonrot       ! Longitude of southern pole of rotation
  REAL :: angrot       ! Angle of rotation

  REAL :: latstr       ! Latitude  of the pole of stretching
  REAL :: lonstr       ! Longitude of the pole of stretching
  REAL :: facstr       ! Stretching factor

  INTEGER :: scanmode  ! Scanning indicator
  INTEGER :: iscan     ! x-direction   scanning indicator
  INTEGER :: jscan     ! y-direction   scanning indicator
  INTEGER :: kscan     ! FORTRAN index scanning indicator

  INTEGER :: ires      ! Resolution direction increments indicator
  INTEGER :: iearth    ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp     ! (u,v) components decomposition indicator

  INTEGER :: jpenta    ! J-Pentagonal resolution parameter
  INTEGER :: kpenta    ! K-Pentagonal resolution parameter
  INTEGER :: mpenta    ! M-Pentagonal resolution parameter
  INTEGER :: ispect    ! Spectral representation type
  INTEGER :: icoeff    ! Spectral coefficient storage mode

  REAL :: xp           ! X coordinate of sub-satellite point
  REAL :: yp           ! Y coordinate of sub-satellite point
  REAL :: xo           ! X coordinate of image sector origin
  REAL :: yo           ! Y coordinate of image sector origin
  REAL :: zo           ! Camera altitude from center of Earth

  REAL :: missing      ! Flag to denote a missing datum
  DATA missing/999.0/

  CHARACTER (LEN=2) :: endstr
  DATA endstr/'  '/

  INTEGER :: i,j,jj    ! Loop indices
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
!  Flag variables initially as missing.
!
!-----------------------------------------------------------------------
!
  gthin   = nint(missing)

  ni      = nint(missing)
  nj      = nint(missing)
  np      = nint(missing)

  nk      = nint(missing)
  npeq    = nint(missing)

  pi      =      missing
  pj      =      missing
  ipole   = nint(missing)

  di      =      missing
  dj      =      missing

  latsw   =      missing
  lonsw   =      missing
  latne   =      missing
  lonne   =      missing

  lattru1 =      missing
  lattru2 =      missing
  lontrue =      missing

  latrot  =      missing
  lonrot  =      missing
  angrot  =      missing

  latstr  =      missing
  lonstr  =      missing
  facstr  =      missing

  iscan   = nint(missing)
  jscan   = nint(missing)
  kscan   = nint(missing)

  ires    = nint(missing)
  iearth  = nint(missing)
  icomp   = nint(missing)

  jpenta  = nint(missing)
  kpenta  = nint(missing)
  mpenta  = nint(missing)
  ispect  = nint(missing)
  icoeff  = nint(missing)

  xp      =      missing
  yp      =      missing
  xo      =      missing
  yo      =      missing
  zo      =      missing
!
!-----------------------------------------------------------------------
!
!  Table B. Grid Identification - Master List of NMC Storage Grids,
!  designated by Parameter Description Section (PDS) octet #7 in
!  Section 1 of GRIB message.
!
!  NOTE:  (1) Grids 21-26 and 61-64 are International Exchange grids.
!     (2) Grids 37-44 are 1.25 x 1.25 degree (lat,lon) "thinned"
!         grids, covering the globe by octants of 3447 points.
!         The number of points on each latitudinal row is given
!         by the formula:
!            Npoints = ifix ( 2 + (90/1.25) * cos(Lat) )
!         The latitudinal increment is always 1.25 degrees; this
!         results in 73 rows for each octant. In GRIB terms,
!         these are also known as "quasi-regular" grids.
!     (3) Grids 67-71 and 72-74 are OFFICE NOTE 84 grid types,
!         for local use by NMC. They are NOT International.
!     (4) Grids 201-2nn are AWIPS grids.
!     (5) A grid ID of 255 indicates that the grid is not defined
!         and its description is provided in the Grid Description
!         Section (GDS).
!
!-----------------------------------------------------------------------
!
  GO TO ( 1,  2,  3,  4,  5,  6,254,254,254,254,254,254,254,254,254,    &
       254,254,254,254,254, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,     &
       254,254, 33, 34,254,254, 37, 38, 39, 40, 41, 42, 43, 44, 45,     &
       254,254,254,254, 50,254,254,254,254, 55, 56,254,254,254,254,     &
        61, 62, 63, 64,254,254, 67, 68, 69, 70, 71, 72, 73, 74, 75,     &
        76, 77,254,254,254,254,254,254,254, 85, 86, 87,254,254, 90,     &
        91, 92, 93, 94, 95, 96, 97, 98,254,100,101,254,103,104,105,     &
       106,107,254,254,254,254,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,126,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,201,202,203,204,205,206,207,208,209,210,     &
       211,212,213,214,254,254,254,218,254,254,221,254,254,254,254,     &
       254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,     &
       254,254,254,254,254,254,254,254,254,254,254,254,254,254,255)     &
       gridtyp                                           !TINA added 221
!
!-----------------------------------------------------------------------
!
!  Mercator Tropical Strip 5 deg Longitude
!
!-----------------------------------------------------------------------
!

  1     CONTINUE
  gridesc = 'Mercator Tropical Strip 5 deg Longitude'//endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        1
  gthin   =        0
  ni      =       73
  nj      =       23
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =   513669.
  dj      =   513669.
  latsw   =   -48.09
  lonsw   =     0.00
  latne   =    48.09
  lonne   =   360.00
  lattru1 =    22.50
  lattru2 =  missing
  lontrue =     0.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global 2.5 degree Longitude,Latitude
!
!-----------------------------------------------------------------------
!

  2     CONTINUE
  gridesc = 'Global 2.5 degree Longitude,Latitude' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      144
  nj      =       73
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.50
  dj      =     2.50
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =    -2.50
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global 1.0 degree Longitude,Latitude
!
!-----------------------------------------------------------------------
!

  3     CONTINUE
  gridesc = 'Global 1.0 degree Longitude,Latitude' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      360
  nj      =      181
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     1.00
  dj      =     1.00
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =    -1.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global 0.5 degree Longitude,Latitude
!
!-----------------------------------------------------------------------
!

  4     CONTINUE
  gridesc = 'Global 0.5 degree Longitude,Latitude' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      720
  nj      =      361
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     0.50
  dj      =     0.50
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =    -0.50
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  LFM Analysis N Hem Polar Stereo 190.5 km
!
!-----------------------------------------------------------------------
!

  5     CONTINUE
  gridesc = 'LFM Analysis N Hem Polar Stereo 190.5 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       53
  nj      =       57
  np      =  ni * nj
  pi      =       27.
  pj      =       49.
  di      =   190500.
  dj      =   190500.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!   Grid = 'LFM Forecast N Hem Polar Stereo 190.5 km
!
!-----------------------------------------------------------------------
!

  6     CONTINUE
  gridesc = 'LFM Forecast N Hem Polar Stereo 190.5 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       53
  nj      =       45
  np      =  ni * nj
  pi      =       27.
  pj      =       49.
  di      =   190500.
  dj      =   190500.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!   Grid = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  21    CONTINUE
  gridesc = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =       37
  nj      =       37
  np      =     1333
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     2.50
  latsw   =     0.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =   180.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Grid = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  22    CONTINUE
  gridesc = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =       37
  nj      =       37
  np      =     1333
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     2.50
  latsw   =     0.00
  lonsw   =  -180.00
  latne   =    90.00
  lonne   =     0.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Grid = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  23    CONTINUE
  gridesc = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =       38
  nj      =       36
  np      =     1333
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     2.50
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =     0.00
  lonne   =   180.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exchange FOS 5 by 2.5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  24    CONTINUE
  gridesc = 'Int Exchange FOS 5 by 2.5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =       38
  nj      =       36
  np      =     1333
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     2.50
  latsw   =   -90.00
  lonsw   =  -180.00
  latne   =     0.00
  lonne   =     0.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Grid = 'Int Exchange FOS 5 by 5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  25    CONTINUE
  gridesc = 'Int Exchange FOS 5 by 5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =       72
  nj      =       19
  np      =     1297
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     5.00
  latsw   =     0.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =   355.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exchange FOS 5 by 5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  26    CONTINUE
  gridesc = 'Int Exchange FOS 5 by 5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =       73
  nj      =       18
  np      =     1297
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     5.00
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =     0.00
  lonne   =   355.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  North Hemis Polar Stereographic 381 km
!
!-----------------------------------------------------------------------
!

  27    CONTINUE
  gridesc = 'North Hemis Polar Stereographic 381 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       65
  nj      =       65
  np      =  ni * nj
  pi      =       33.
  pj      =       33.
  di      =   381000.
  dj      =   381000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =   -80.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  South Hemis Polar Stereographic 381 km
!
!-----------------------------------------------------------------------
!

  28    CONTINUE
  gridesc = 'South Hemis Polar Stereographic 381 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       65
  nj      =       65
  np      =  ni * nj
  pi      =       33.
  pj      =       33.
  di      =   381000.
  dj      =   381000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =   -60.00
  lattru2 =  missing
  lontrue =   100.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Northern Hemisphere 2.5 degree Lat,Lon
!
!-----------------------------------------------------------------------
!

  29    CONTINUE
  gridesc = 'Northern Hemisphere 2.5 degree Lat,Lon' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      145
  nj      =       37
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.50
  dj      =     2.50
  latsw   =     0.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =   360.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Southern Hemisphere 2.5 degree Lat,Lon
!
!-----------------------------------------------------------------------
!

  30    CONTINUE
  gridesc = 'Southern Hemisphere 2.5 degree Lat,Lon' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      145
  nj      =       37
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.50
  dj      =     2.50
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =     0.00
  lonne   =   360.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Northern Hemisphere 2.0 degree Lat,Lon
!
!-----------------------------------------------------------------------
!

  33    CONTINUE
  gridesc = 'Northern Hemisphere 2.0 degree Lat,Lon' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      181
  nj      =       46
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     2.00
  latsw   =     0.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =   360.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Southern Hemisphere 2.0 degree Lat,Lon
!
!-----------------------------------------------------------------------
!

  34    CONTINUE
  gridesc = 'Southern Hemisphere 2.0 degree Lat,Lon' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      181
  nj      =       46
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     2.00
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =     0.00
  lonne   =   360.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-I 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  37    CONTINUE
  gridesc = 'Int Xch FOS-I 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =     0.00
  lonsw   =   -30.00
  latne   =    90.00
  lonne   =    60.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-J 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  38    CONTINUE
  gridesc = 'Int Xch FOS-J 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =     0.00
  lonsw   =    60.00
  latne   =    90.00
  lonne   =   150.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-K 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  39    CONTINUE
  gridesc = 'Int Xch FOS-K 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =     0.00
  lonsw   =   150.00
  latne   =    90.00
  lonne   =  -120.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-L 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  40    CONTINUE
  gridesc = 'Int Xch FOS-L 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =     0.00
  lonsw   =  -120.00
  latne   =    90.00
  lonne   =   -30.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-M 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  41    CONTINUE
  gridesc = 'Int Xch FOS-M 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =   -90.00
  lonsw   =   -30.00
  latne   =     0.00
  lonne   =    60.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-N 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  42    CONTINUE
  gridesc = 'Int Xch FOS-N 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =   -90.00
  lonsw   =    60.00
  latne   =     0.00
  lonne   =   150.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-O 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  43    CONTINUE
  gridesc = 'Int Xch FOS-O 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =   -90.00
  lonsw   =   150.00
  latne   =     0.00
  lonne   =  -120.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Xch FOS-P 1.25x1.25 deg thin Lon,Lat
!
!-----------------------------------------------------------------------
!

  44    CONTINUE
  gridesc = 'Int Xch FOS-P 1.25x1.25 deg thin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       73
  nj      =       73
  np      =     3447
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =   -90.00
  lonsw   =  -120.00
  latne   =     0.00
  lonne   =   -30.00
  lattru1 =     0.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global 1.25 degree Latitude,Longitude
!
!-----------------------------------------------------------------------
!

  45    CONTINUE
  gridesc = 'Global 1.25 degree Latitude,Longitude' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      288
  nj      =      145
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     1.25
  dj      =     1.25
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =    -1.25
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exch FOS 2.5x1.25 deg USthin Lon,Lat
!
!-----------------------------------------------------------------------
!

  50    CONTINUE
  gridesc = 'Int Exch FOS 2.5x1.25 deg USthin Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        1
  ni      =       36
  nj      =       33
  np      =      964
  pi      =  missing
  pj      =  missing
  di      =     5.00
  dj      =     1.25
  latsw   =    20.00
  lonsw   =  -140.00
  latne   =    60.00
  lonne   =   -52.50
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  North Amer SFC Anal Polar Stereo 254 km
!
!-----------------------------------------------------------------------
!

  55    CONTINUE
  gridesc = 'North Amer SFC Anal Polar Stereo 254 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       87
  nj      =       71
  np      =  ni * nj
  pi      =       44.
  pj      =       38.
  di      =   254000.
  dj      =   254000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  North Amer Sfc Anal Polar Stereo 127 km
!
!-----------------------------------------------------------------------
!

  56    CONTINUE
  gridesc = 'North Amer Sfc Anal Polar Stereo 127 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       87
  nj      =       71
  np      =  ni * nj
  pi      =       40.
  pj      =       73.
  di      =   127000.
  dj      =   127000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exchange FOS 2 by 2 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  61    CONTINUE
  gridesc = 'Int Exchange FOS 2 by 2 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        0
  ni      =       91
  nj      =       46
  np      =     4096
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     2.00
  latsw   =     0.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =   180.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exchange FOS 2 by 2 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  62    CONTINUE
  gridesc = 'Int Exchange FOS 2 by 2 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        0
  ni      =       91
  nj      =       46
  np      =     4096
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     2.00
  latsw   =     0.00
  lonsw   =  -180.00
  latne   =    90.00
  lonne   =     0.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exchange FOS 2 by 2 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  63    CONTINUE
  gridesc = 'Int Exchange FOS 2 by 2 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        0
  ni      =       92
  nj      =       45
  np      =     4096
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     2.00
  latsw   =   -90.00
  lonsw   =     0.00
  latne   =     0.00
  lonne   =   180.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Int Exchange FOS 2 by 2 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  64    CONTINUE
  gridesc = 'Int Exchange FOS 2 by 2 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        4
  gthin   =        0
  ni      =       92
  nj      =       45
  np      =     4096
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     2.00
  latsw   =   -90.00
  lonsw   =  -180.00
  latne   =     0.00
  lonne   =     0.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  NW Atlantic Polar Stereo 23.8125 km
!
!-----------------------------------------------------------------------
!

  67    CONTINUE
  gridesc = 'NW Atlantic Polar Stereo 23.8125 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      117
  nj      =      117
  np      =  ni * nj
  pi      =        9.
  pj      =      317.
  di      =  23812.5
  dj      =  23812.5
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =   -80.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Gulf of Mexico Polar Stereo 23.8125 km
!
!-----------------------------------------------------------------------
!

  68    CONTINUE
  gridesc = 'Gulf of Mexico Polar Stereo 23.8125 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      117
  nj      =      117
  np      =  ni * nj
  pi      =      -35.
  pj      =      361.
  di      =  23812.5
  dj      =  23812.5
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Gulf of Alaska Polar Stereo 23.8125 km
!
!-----------------------------------------------------------------------
!

  69    CONTINUE
  gridesc = 'Gulf of Alaska Polar Stereo 23.8125 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      117
  nj      =      117
  np      =  ni * nj
  pi      =      177.
  pj      =      209.
  di      =  23812.5
  dj      =  23812.5
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Calif Pacific Polar Stereo 23.8125 km
!
!-----------------------------------------------------------------------
!

  70    CONTINUE
  gridesc = 'Calif Pacific Polar Stereo 23.8125 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      117
  nj      =      117
  np      =  ni * nj
  pi      =      169.
  pj      =      285.
  di      =  23812.5
  dj      =  23812.5
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Mexican Pacific Polar Stereo 23.8125 km
!
!-----------------------------------------------------------------------
!

  71    CONTINUE
  gridesc = 'Mexican Pacific Polar Stereo 23.8125 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      117
  nj      =      117
  np      =  ni * nj
  pi      =      137.
  pj      =      377.
  di      =  23812.5
  dj      =  23812.5
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Mercator 2.5 degree Longitude
!
!-----------------------------------------------------------------------
!

  72    CONTINUE
  gridesc = 'Mercator 2.5 degree Longitude' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        1
  gthin   =        0
  ni      =       29
  nj      =       14
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.50
  dj      =     2.50
  latsw   =    46.40
  lonsw   =   170.00
  latne   =    64.40
  lonne   =  -120.00
  lattru1 =    22.50
  lattru2 =  missing
  lontrue =     0.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global Gaussian R40 wave trans Lon,Lat
!
!-----------------------------------------------------------------------
!

  73    CONTINUE
  gridesc = 'Global Gaussian R40 wave trans Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =       50
  gthin   =        0
  ni      =      128
  nj      =      102
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =  missing
  dj      =  missing
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        0
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  North Hemisphere 2 by 1.5 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  74    CONTINUE
  gridesc = 'North Hemisphere 2 by 1.5 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      180
  nj      =       60
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     2.00
  dj      =     1.50
  latsw   =     0.00
  lonsw   =     0.00
  latne   =    90.00
  lonne   =     0.00
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  QLM Hurricane NH Lambert Conformal 40 km
!
!-----------------------------------------------------------------------
!

  75    CONTINUE
  gridesc = 'QLM Hurricane NH Lambert Conformal 40 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =      111
  nj      =      111
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =    40000.
  dj      =    40000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    30.00
  lattru2 =    60.00
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  QLM Hurricane SH Lambert Conformal 40 km
!
!-----------------------------------------------------------------------
!

  76    CONTINUE
  gridesc = 'QLM Hurricane SH Lambert Conformal 40 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =      111
  nj      =      111
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =    40000.
  dj      =    40000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =   -30.00
  lattru2 =   -60.00
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  QLM Hurricane North Hemis Mercator 40 km
!
!-----------------------------------------------------------------------
!

  77    CONTINUE
  gridesc = 'QLM Hurricane North Hemis Mercator 40 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =      111
  nj      =      111
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =    40000.
  dj      =    40000.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    -22.5
  lattru2 =     22.5
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Northern Hemisphere 1.0 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  85    CONTINUE
  gridesc = 'Northern Hemisphere 1.0 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      360
  nj      =       90
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     1.00
  dj      =     1.00
  latsw   =     0.50
  lonsw   =     0.50
  latne   =    89.50
  lonne   =    -0.50
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Southern Hemisphere 1.0 degree Lon,Lat
!
!-----------------------------------------------------------------------
!

  86    CONTINUE
  gridesc = 'Southern Hemisphere 1.0 degree Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        0
  gthin   =        0
  ni      =      360
  nj      =       90
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =     1.00
  dj      =     1.00
  latsw   =   -89.50
  lonsw   =     0.50
  latne   =    -0.50
  lonne   =    -0.50
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  RUC/MAPS NH Polar Stereo 60 km @ 60N
!
!-----------------------------------------------------------------------
!

  87    CONTINUE
  gridesc = 'RUC/MAPS NH Polar Stereo 60 km @ 60N' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       81
  nj      =       62
  np      =  ni * nj
  pi      =    31.91
  pj      =   112.53
  di      =    68153.
  dj      =    68153.
  latsw   =  22.8756
  lonsw   =-120.4911
  latne   =  46.0172
  lonne   = -60.8284
  lattru1 =    40.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta SS 80 km rot 15/26x14/26 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  90    CONTINUE
  gridesc = 'Eta SS 80 km rot 15/26x14/26 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      201
  gthin   =        1
  ni      =       92
  nj      =      141
!        Np      =  Ni * Nj - Nj/2        ! Arakawa E-Grid Staggering
  np      =    12902
  pi      =  missing
  pj      =  missing
  di      =   15./26.
  dj      =   14./26.
  latsw   =    0.182
  lonsw   = -149.887
  latne   =  missing
  lonne   =  missing
  lattru1 =    52.00
  lattru2 =  missing
  lontrue =  -111.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        0
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta F 80 km rot 15/26x14/26 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  91    CONTINUE
  gridesc = 'Eta F 80 km rot 15/26x14/26 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      202
  gthin   =        0
  ni      =      183
  nj      =      141
!        Np      =  Ni * Nj               ! Arakawa E-Grid Filled
  np      =    25803
  pi      =  missing
  pj      =  missing
  di      =   15./26.
  dj      =   14./26.
  latsw   =    0.182
  lonsw   = -149.887
  latne   =  missing
  lonne   =  missing
  lattru1 =    52.00
  lattru2 =  missing
  lontrue =  -111.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta SS 40 km rot 5/18x15/57 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  92    CONTINUE
  gridesc = 'Eta SS 40 km rot 5/18x15/57 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      201
  gthin   =        1
  ni      =      127
  nj      =      191
!        Np      =  Ni * Nj - Nj/2        ! Arakawa E-Grid Staggering
  np      =    24162
  pi      =  missing
  pj      =  missing
  di      =    5./18.
  dj      =   15./57.
  latsw   =    9.678
  lonsw   = -128.826
  latne   =  missing
  lonne   =  missing
  lattru1 =    41.00
  lattru2 =  missing
  lontrue =   -97.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        0
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta F 40 km rot 5/18 x 15/57 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  93    CONTINUE
  gridesc = 'Eta F 40 km rot 5/18 x 15/57 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      202
  gthin   =        0
  ni      =      253
  nj      =      191
!        Np      =  Ni * Nj               ! Arakawa E-Grid Filled
  np      =    48323
  pi      =  missing
  pj      =  missing
  di      =    5./18.
  dj      =   15./57.
  latsw   =    9.678
  lonsw   = -128.826
  latne   =  missing
  lonne   =  missing
  lattru1 =    41.00
  lattru2 =  missing
  lontrue =   -97.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta SS 29 km rot 7/36 x 5/27 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  94    CONTINUE
  gridesc = 'Eta SS 29 km rot 7/36 x 5/27 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      201
  gthin   =        1
  ni      =      181
  nj      =      271
!        Np      =  Ni * Nj - Nj/2        ! Arakawa E-Grid Staggering
  np      =    48916
  pi      =  missing
  pj      =  missing
  di      =    7./36.
  dj      =    5./27.
  latsw   =    9.678
  lonsw   = -128.826
  latne   =  missing
  lonne   =  missing
  lattru1 =    41.00
  lattru2 =  missing
  lontrue =   -97.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        0
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta F 29 km rot 7/36 x 5/27 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  95    CONTINUE
  gridesc = 'Eta F 29 km rot 7/36 x 5/27 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      202
  gthin   =        0
  ni      =      361
  nj      =      271
!        Np      =  Ni * Nj               ! Arakawa E-Grid Filled
  np      =    97831
  pi      =  missing
  pj      =  missing
  di      =    7./36.
  dj      =    5./27.
  latsw   =    9.678
  lonsw   = -128.826
  latne   =  missing
  lonne   =  missing
  lattru1 =    41.00
  lattru2 =  missing
  lontrue =   -97.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta SS 48 km rot 1/3 x 4/13 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  96    CONTINUE
  gridesc = 'Eta SS 48 km rot 1/3 x 4/13 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      201
  gthin   =        1
  ni      =      160
  nj      =      261
!        Np      =  Ni * Nj - Nj/2        ! Arakawa E-Grid Staggering
  np      =    41630
  pi      =  missing
  pj      =  missing
  di      =     1./3.
  dj      =    4./13.
  latsw   =   -3.441
  lonsw   = -148.799
  latne   =  missing
  lonne   =  missing
  lattru1 =    50.00
  lattru2 =  missing
  lontrue =  -110.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        0
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta F 48 km rot 1/3 x 4/13 deg Lon,Lat
!
!-----------------------------------------------------------------------
!

  97    CONTINUE
  gridesc = 'Eta F 48 km rot 1/3 x 4/13 deg Lon,Lat' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =      202
  gthin   =        0
  ni      =      319
  nj      =      261
!        Np      =  Ni * Nj               ! Arakawa E-Grid Filled
  np      =    83259
  pi      =  missing
  pj      =  missing
  di      =     1./3.
  dj      =    4./13.
  latsw   =   -3.441
  lonsw   = -148.799
  latne   =  missing
  lonne   =  missing
  lattru1 =    50.00
  lattru2 =  missing
  lontrue =  -110.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global Gaussian T62 wave 1.875 degree
!
!-----------------------------------------------------------------------
!

  98    CONTINUE
  gridesc = 'Global Gaussian T62 wave 1.875 degree' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =       50
  gthin   =        0
  ni      =      192
  nj      =       94
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =    1.875
  dj      =    1.875
  latsw   =  -88.542
  lonsw   =    0.000
  latne   =   88.542
  lonne   =   -1.875
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  NGM Original-C Polar Stereo 91.452 km
!
!-----------------------------------------------------------------------
!

  100   CONTINUE
  gridesc = 'NGM Original-C Polar Stereo 91.452 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       83
  nj      =       83
  np      =  ni * nj
  pi      =     40.5
  pj      =     88.5
  di      =    91452.
  dj      =    91452.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  NGM Big-C Polar Stereographic 91.452 km
!
!-----------------------------------------------------------------------
!

  101   CONTINUE
  gridesc = 'NGM Big-C Polar Stereographic 91.452 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      113
  nj      =       91
  np      =  ni * nj
  pi      =     58.5
  pj      =     92.5
  di      =    91452.
  dj      =    91452.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  ARL Northrn Hemis Polar Stereo 91.452 km
!
!-----------------------------------------------------------------------
!

  103   CONTINUE
  gridesc = 'ARL Northrn Hemis Polar Stereo 91.452 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       65
  nj      =       56
  np      =  ni * nj
  pi      =     25.5
  pj      =     84.5
  di      =    91452.
  dj      =    91452.
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  NGM Super-C NH Polar Stereo 90.75464 km
!
!-----------------------------------------------------------------------
!

  104   CONTINUE
  gridesc = 'NGM Super-C NH Polar Stereo 90.75464 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      147
  nj      =      110
  np      =  ni * nj
  pi      =     75.5
  pj      =    109.5
  di      = 90754.64
  dj      = 90754.64
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta Super-C NH Polar Stereo 90.75464 km
!
!-----------------------------------------------------------------------
!

  105   CONTINUE
  gridesc = 'Eta Super-C NH Polar Stereo 90.75464 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       83
  nj      =       83
  np      =  ni * nj
  pi      =     40.5
  pj      =     88.5
  di      = 90754.64
  dj      = 90754.64
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Eta Super-C2 NH Polar Stereo 45.37732 km
!
!-----------------------------------------------------------------------
!

  106   CONTINUE
  gridesc = 'Eta Super-C2 NH Polar Stereo 45.37732 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      165
  nj      =      117
  np      =  ni * nj
  pi      =       80.
  pj      =      176.
  di      = 45377.32
  dj      = 45377.32
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  RUC Super-C2 NH Polar Stereo 45.37732 km
!
!-----------------------------------------------------------------------
!

  107   CONTINUE
  gridesc = 'RUC Super-C2 NH Polar Stereo 45.37732 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      120
  nj      =       92
  np      =  ni * nj
  pi      =       46.
  pj      =      167.
  di      = 45377.32
  dj      = 45377.32
  latsw   =  missing
  lonsw   =  missing
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  Global Gaussian T126 wave 0.9375 degree
!
!-----------------------------------------------------------------------
!

  126   CONTINUE
  gridesc = 'Global Gaussian T126 wave 0.9375 degree' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =       50
  gthin   =        0
  ni      =      384
  nj      =      190
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =   0.9375
  dj      =   0.9375
  latsw   = -89.2770
  lonsw   =   0.0000
  latne   =  89.2770
  lonne   =  -0.9375
  lattru1 =  missing
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        0
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-A North Hemis Polar Stereo 381 km
!
!-----------------------------------------------------------------------
!

  201   CONTINUE
  gridesc = 'AWIPS-A North Hemis Polar Stereo 381 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       65
  nj      =       65
  np      =  ni * nj
  pi      =       33.
  pj      =       33.
  di      =   381000.
  dj      =   381000.
  latsw   =  -20.826
  lonsw   = -150.000
  latne   =  missing
  lonne   =  missing
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-I Nat-CONUS Polar Stereo 190.5 km
!
!-----------------------------------------------------------------------
!

  202   CONTINUE
  gridesc = 'AWIPS-I Nat-CONUS Polar Stereo 190.5 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       65
  nj      =       43
  np      =  ni * nj
  pi      =       33.
  pj      =       45.
  di      =   190500.
  dj      =   190500.
  latsw   =    7.838
  lonsw   = -141.028
  latne   =   35.617
  lonne   =  -18.576
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-J Nat-Alaska Polar Stereo 190.5 km
!
!-----------------------------------------------------------------------
!

  203   CONTINUE
  gridesc = 'AWIPS-J Nat-Alaska Polar Stereo 190.5 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       45
  nj      =       39
  np      =  ni * nj
  pi      =       27.
  pj      =       37.
  di      =   190500.
  dj      =   190500.
  latsw   =   19.132
  lonsw   =  174.163
  latne   =   57.634
  lonne   =  -53.660
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -150.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-K Hawaii Mercator 1.531 deg 160 km
!
!-----------------------------------------------------------------------
!

  204   CONTINUE
  gridesc = 'AWIPS-K Hawaii Mercator 1.531 deg 160 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        1
  gthin   =        0
  ni      =       93
  nj      =       68
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =   160000.
  dj      =   160000.
  latsw   =  -25.000
  lonsw   =  110.000
  latne   =   60.644
  lonne   = -109.129
  lattru1 =    20.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-L PuertoRico Polar Stereo 190.5 km
!
!-----------------------------------------------------------------------
!

  205   CONTINUE
  gridesc = 'AWIPS-L PuertoRico Polar Stereo 190.5 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       45
  nj      =       39
  np      =  ni * nj
  pi      =       27.
  pj      =       57.
  di      =   190500.
  dj      =   190500.
  latsw   =    0.616
  lonsw   =  -84.904
  latne   =   45.620
  lonne   =  -15.000
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =   -60.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-M US-MARD Lambert Conf 81.2705 km
!
!-----------------------------------------------------------------------
!

  206   CONTINUE
  gridesc = 'AWIPS-M US-MARD Lambert Conf 81.2705 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =       51
  nj      =       41
  np      =  ni * nj
  pi      =   30.000
  pj      =  169.745
  di      =  81270.5
  dj      =  81270.5
  latsw   =   22.289
  lonsw   = -117.991
  latne   =   51.072
  lonne   =  -73.182
  lattru1 =    25.00
  lattru2 =    25.00
  lontrue =   -95.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-N Reg-Alaska Polar Stereo 95.25 km
!
!-----------------------------------------------------------------------
!

  207   CONTINUE
  gridesc = 'AWIPS-N Reg-Alaska Polar Stereo 95.25 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       49
  nj      =       35
  np      =  ni * nj
  pi      =       25.
  pj      =       51.
  di      =    95250.
  dj      =    95250.
  latsw   =   42.085
  lonsw   = -175.641
  latne   =   63.976
  lonne   =  -93.689
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -150.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-O Hawaii Mercator 0.766 deg 80 km
!
!-----------------------------------------------------------------------
!

  208   CONTINUE
  gridesc = 'AWIPS-O Hawaii Mercator 0.766 deg 80 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        1
  gthin   =        0
  ni      =       29
  nj      =       27
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =    80000.
  dj      =    80000.
  latsw   =    9.343
  lonsw   = -167.315
  latne   =   28.092
  lonne   = -145.878
  lattru1 =    20.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-S RegMARD Lambert Conf 40.63525 km
!
!-----------------------------------------------------------------------
!

  209   CONTINUE
  gridesc = 'AWIPS-S RegMARD Lambert Conf 40.63525 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =      101
  nj      =       81
  np      =  ni * nj
  pi      =    59.00
  pj      =   338.49
  di      = 40635.25
  dj      = 40635.25
  latsw   =   22.289
  lonsw   = -117.991
  latne   =   51.072
  lonne   =  -73.182
  lattru1 =    25.00
  lattru2 =    25.00
  lontrue =   -95.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-P PuertoRico Mercator .766deg 80km
!
!-----------------------------------------------------------------------
!

  210   CONTINUE
  gridesc = 'AWIPS-P PuertoRico Mercator .766deg 80km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        1
  gthin   =        0
  ni      =       25
  nj      =       25
  np      =  ni * nj
  pi      =  missing
  pj      =  missing
  di      =    80000.
  dj      =    80000.
  latsw   =    9.000
  lonsw   =  -77.000
  latne   =   26.422
  lonne   =  -58.625
  lattru1 =    20.00
  lattru2 =  missing
  lontrue =  missing
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        0
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-Q RegCONUS Lambert Conf 81.2705 km
!
!-----------------------------------------------------------------------
!

  211   CONTINUE
  gridesc = 'AWIPS-Q RegCONUS Lambert Conf 81.2705 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =       93
  nj      =       65
  np      =  ni * nj
  pi      =   53.000
  pj      =  178.745
  di      =  81270.5
  dj      =  81270.5
  latsw   =   12.190
  lonsw   = -133.459
  latne   =   57.290
  lonne   =  -49.385
  lattru1 =    25.00
  lattru2 =    25.00
  lontrue =   -95.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-R RegCONUS Lambert Conf 40.63525km
!
!-----------------------------------------------------------------------
!

  212   CONTINUE
  gridesc = 'AWIPS-R RegCONUS Lambert Conf 40.63525km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        3
  gthin   =        0
  ni      =      185
  nj      =      129
  np      =  ni * nj
  pi      =   105.00
  pj      =   356.49
  di      = 40635.25
  dj      = 40635.25
  latsw   =   12.190
  lonsw   = -133.459
  latne   =   57.290
  lonne   =  -49.385
  lattru1 =    25.00
  lattru2 =    25.00
  lontrue =   -95.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-H Nat-CONUS Polar Stereo 95.250 km
!
!-----------------------------------------------------------------------
!

  213   CONTINUE
  gridesc = 'AWIPS-H Nat-CONUS Polar Stereo 95.250 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =      129
  nj      =       85
  np      =  ni * nj
  pi      =       65.
  pj      =       89.
  di      =    95250.
  dj      =    95250.
  latsw   =    7.838
  lonsw   = -141.028
  latne   =   35.617
  lonne   =  -18.577
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -105.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-T RegAlaska Polar Stereo 47.625 km
!
!-----------------------------------------------------------------------
!

  214   CONTINUE
  gridesc = 'AWIPS-T RegAlaska Polar Stereo 47.625 km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information

  iproj   =        5
  gthin   =        0
  ni      =       97
  nj      =       69
  np      =  ni * nj
  pi      =       49.
  pj      =      101.
  di      =    47625.
  dj      =    47625.
  latsw   =   42.085
  lonsw   = -175.641
  latne   =   63.975
  lonne   =  -93.689
  lattru1 =    60.00
  lattru2 =  missing
  lontrue =  -150.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8
  RETURN
!
!-----------------------------------------------------------------------
!
!  AWIPS-R RegCONUS Lambert Conf 40.63525km
!
!-----------------------------------------------------------------------
!

  218   CONTINUE
  gridesc = 'AWIPS-R RegCONUS Lambert Conf 12km' // endstr
  GO TO 256     ! Use GDS to get grid information
  RETURN

! 
!-----------------------------------------------------------------------
! 
!TINA
!
!  AWIPS 221 Regional NOAMHI - High-res North American Master Grid 
!  (Lambert Conformal) 32 km
! 
!-----------------------------------------------------------------------
! 
  
  221   CONTINUE
  gridesc = 'AWIPS 221 NARR grid Lambert Conf 32km' // endstr
  IF ( gdsflg == 1 ) GO TO 256     ! Use GDS to get grid information
  !the stuff below does not seem to get used
  iproj   =        3
  gthin   =        0
  ni      =      349
  nj      =      277
  np      =  ni * nj
  pi      =  174.507
  pj      =  307.764
  di      = 32463.41
  dj      = 32463.41 
  latsw   =    1.000
  lonsw   = -145.500
  latne   =   46.352
  lonne   =   -2.566
  lattru1 =    50.00
  lattru2 =    50.00
  lontrue =  -107.00
  iscan   =        0
  jscan   =        1
  kscan   =        0
  ires    =        1
  iearth  =        0
  icomp   =        8 !don't know what to set this to!!!
  RETURN
!
!-----------------------------------------------------------------------
!
!  UNDEFINED NMC STORAGE GRID - CONSULT GDS
!
!-----------------------------------------------------------------------
!

  254   CONTINUE

  gridesc = 'UNDEFINED NMC STORAGE GRID - CONSULT GDS' // endstr

  GO TO 256
!
!-----------------------------------------------------------------------
!
!  NON-STANDARD GRID - CONSULT GDS
!
!-----------------------------------------------------------------------
!

  255   CONTINUE

  gridesc = 'NON-STANDARD GRID - CONSULT GDS' // endstr
!
!-----------------------------------------------------------------------
!
!  Grid is NOT in the Master List of NMC Storage Grids! Crack GDS.
!
!-----------------------------------------------------------------------
!

  256   CONTINUE                   ! End of computed GOTO gridtyp

!
!-----------------------------------------------------------------------
!
!  Section 2 (GDS) octet mapping depends on type of projection.
!
!  NOTE:  GRIBEX decoder locations offset from WMO/NMC standard !!!
!
!-----------------------------------------------------------------------
!
  iproj = IABS(igds(1))
!
!-----------------------------------------------------------------------
!
!  Latitude/Longitude grid:
!
!-----------------------------------------------------------------------
!
  IF (iproj ==  0 .OR. iproj == 10 .OR. iproj == 20 .OR. iproj == 30) THEN
    gthin = igds(17)

    ni = igds(2)
    nj = igds(3)
    np =  ni * nj

    latsw = .001 * FLOAT(MIN0(igds(4),igds(7)))
    latne = .001 * FLOAT(MAX0(igds(4),igds(7)))

! Commented out by Y. Wang to handle no fixed dimension source from FTP2U
! see arps.input for extdopt > 100
!
!    IF (MOD(igds(5)+360000,360000) <= MOD(igds(8)+360000,360000) ) THEN
      lonsw = .001 * FLOAT(igds(5))
      lonne = .001 * FLOAT(igds(8))
!    ELSE
!      lonsw = .001 * FLOAT(igds(8))
!      lonne = .001 * FLOAT(igds(5))
!    END IF

    ires   = IAND(igds(6),128)/128
    iearth = IAND(igds(6), 64)/64
    icomp  = IAND(igds(6),  8)/8

    IF ( ires == 1 ) THEN
      di = .001 * FLOAT(igds( 9))
      dj = .001 * FLOAT(igds(10))
    ELSE
      di = missing
      dj = missing
    END IF

!wdt RLC
    ! Work around bug in NCAR/NCEP Reanalysis GRIB files. The lowest 3 bits of 
    ! GDS(28) (scanning mode flag) are incorrectly set to 1,1,1, and other bits
    ! are set randomly. Confirmed by Chi-Fan Shih of NCAR, 2003-08-12.
    IF (igds(11) /= 0 .AND. igds(2) == 144 .AND. igds(3) == 73) THEN
      PRINT '(2A,I4)', 'WARNING: grbgrid: Scanning mode [gds(28)] is ', &
                       'non-zero:', igds(11)
      PRINT '(2A)',   '  Setting it to zero to fix bug in ', &
                         'NCAR/NCEP Reanalysis GRIB files'
      igds(11) = 0
    END IF
    scanmode = igds(11)
    iscan    = IAND(scanmode,128)/128
    jscan    = IAND(scanmode, 64)/64
    kscan    = IAND(scanmode, 32)/32
!
!-----------------------------------------------------------------------
!
!  Gaussian grid:
!
!-----------------------------------------------------------------------
!
  ELSE IF (iproj ==  4 .OR. iproj == 14 .OR.                            &
           iproj == 24 .OR. iproj == 34) THEN
!
!-----------------------------------------------------------------------
!
!  Quasi-regular (thinned) grids introduced in Edition 1.
!
!-----------------------------------------------------------------------
!
    nj    = igds( 3)
    gthin = igds(17)
!    IF (Gthin .eq. 0 .or. iEdit .lt. 1) THEN
    IF (gthin == 0 ) THEN
      ni = igds(2)
    ELSE
      DO jj = 1, nj
        nit(jj) = 0
      END DO
      j = 0
      DO jj = 1, nj
        j = j + 1
        IF (j > nj) GO TO 310
        IF (j == nj) THEN
          nit(j) = igds(j+22)
          GO TO 310
        END IF
        i = 0
        300         i = i + 1
        IF (igds(j+22+1) == igds(j+22)) THEN
          j = j + 1
          GO TO 300
        END IF
        nit(j) = igds(j+22)
      END DO
      310       ni = 0
      np = 0
      DO jj = 1, nj
        ni = MAX0(ni,nit(jj))
        np = np + nit(jj)
      END DO
    END IF

    latsw = .001 * FLOAT(MIN0(igds(4),igds(7)))
    latne = .001 * FLOAT(MAX0(igds(4),igds(7)))
    IF (MOD(igds(5)+360000,360000) <= MOD(igds(8)+360000,360000) ) THEN
      lonsw = .001 * FLOAT(igds(5))
      lonne = .001 * FLOAT(igds(8))
    ELSE
      lonsw = .001 * FLOAT(igds(8))
      lonne = .001 * FLOAT(igds(5))
    END IF


    ires   = IAND(igds(6),128)/128
    iearth = IAND(igds(6), 64)/64
    icomp  = IAND(igds(6),  8)/8

    IF ( ires == 1 ) THEN
      di = .001 * FLOAT(igds(9))
      dj = di
      npeq = igds(10)
    ELSE
      di = missing
      dj = missing
      npeq = nint(missing)
    END IF

    scanmode = igds(11)
    iscan    = IAND(scanmode,128)/128
    jscan    = IAND(scanmode, 64)/64
    kscan    = IAND(scanmode, 32)/32
!
!-----------------------------------------------------------------------
!
!  Spherical harmonic coefficients grid:
!
!-----------------------------------------------------------------------
!
  ELSE IF (iproj == 50 .OR. iproj == 60 .OR.                            &
            iproj == 70 .OR. iproj == 80) THEN
    jpenta = igds( 2)
    kpenta = igds( 3)
    mpenta = igds( 4)
    ispect = igds( 5)
    icoeff = igds( 6)
!
!-----------------------------------------------------------------------
!
!  Polar stereographic grid:
!
!-----------------------------------------------------------------------
!
  ELSE IF (iproj == 5) THEN
    gthin = igds(17)

    ni = igds(2)
    nj = igds(3)
    np =  ni * nj

    latsw = .001 * FLOAT(igds(4))
    lonsw = .001 * FLOAT(igds(5))
    latne = missing
    lonne = missing

    lattru2 = missing
    lontrue = .001 * FLOAT(igds(7))

    di = igds(8)
    dj = igds(9)

    ipole  = IAND(igds(10),128)/128
    IF ( ipole == 0 ) THEN
      lattru1 =  60.
    ELSE
      lattru1 = -60.
    END IF

    ires   = IAND(igds(6),128)/128
    iearth = IAND(igds(6), 64)/64
    icomp  = IAND(igds(6),  8)/8

    scanmode = igds(11)
    iscan    = IAND(scanmode,128)/128
    jscan    = IAND(scanmode, 64)/64
    kscan    = IAND(scanmode, 32)/32
!
!-----------------------------------------------------------------------
!
!  Mercator grid:
!
!-----------------------------------------------------------------------
!
  ELSE IF (iproj ==  1) THEN
    gthin = igds(17)

    ni = igds(2)
    nj = igds(3)
    np =  ni * nj

    latsw = .001 * FLOAT(igds(4))
    lonsw = .001 * FLOAT(igds(5))
    latne = .001 * FLOAT(igds(7))
    lonne = .001 * FLOAT(igds(8))

    lattru1 = .001 * FLOAT(igds(9))
    lattru2 = missing
    lontrue = .001 * FLOAT(igds(7))

    ires   = IAND(igds(6),128)/128
    iearth = IAND(igds(6), 64)/64
    icomp  = IAND(igds(6),  8)/8

    IF ( ires == 1 ) THEN
      di = igds(13)
      dj = igds(14)
    ELSE
      di = missing
      dj = missing
    END IF

    scanmode = igds(11)
    iscan    = IAND(scanmode,128)/128
    jscan    = IAND(scanmode, 64)/64
    kscan    = IAND(scanmode, 32)/32
!
!-----------------------------------------------------------------------
!
!  Lambert conformal, secant or tangent, conical or bi-polar
!      (normal or oblique) or
!  Albers equal-area, secant or tangent, conical or bi-polar
!      (normal or oblique) grid:
!
!-----------------------------------------------------------------------
!
  ELSE IF (iproj == 3 .OR. iproj == 13) THEN
    gthin = igds(17)

    ni = igds(2)
    nj = igds(3)
    np =  ni * nj

    latsw = .001 * FLOAT(igds(4))
    lonsw = .001 * FLOAT(igds(5))
    latne = missing
    lonne = missing

    lontrue = .001 * FLOAT(igds(7))

    di = igds(8)
    dj = igds(9)

    ipole  = IAND(igds(10),128)/128

    lattru1 = .001 * FLOAT(igds(12))
    lattru2 = .001 * FLOAT(igds(13))

    pi = .001 * FLOAT(igds(21))         ! Longitude of South Pole
    pj = .001 * FLOAT(igds(20))         ! Latitude  of South Pole

    ires   = IAND(igds(6),128)/128
    iearth = IAND(igds(6), 64)/64
    icomp  = IAND(igds(6),  8)/8

    scanmode = igds(11)
    iscan    = IAND(scanmode,128)/128
    jscan    = IAND(scanmode, 64)/64
    kscan    = IAND(scanmode, 32)/32
!
!-----------------------------------------------------------------------
!
!  Semi-staggered or filled rotated Arakawa E-grid:
!
!-----------------------------------------------------------------------
!
!  ELSEIF (iProj .eq. 201 .or. iProj .eq. 202) THEN
!    Gthin = igds(17)

!    Ni = igds(7)
!    Nj = igds(8)
!    Np = igds(2)

!    LatSW = .001 * Float(igds(4))
!    LonSW = .001 * Float(igds(5))
!    LatNE = Missing
!    LonNE = Missing

!    iRes   = iand(igds(6),128)/128
!    iEarth = iand(igds(6), 64)/64
!    iComp  = iand(igds(6),  8)/8

!    IF ( iRes.eq.1 ) THEN
!      Di = .001 * Float(igds( 9))
!      Dj = .001 * Float(igds(10))
!    ELSE
!      Di = Missing
!      Dj = Missing
!    ENDIF

!    ScanMode = igds(11)
!    iScan    = iand(ScanMode,128)/128
!    jScan    = iand(ScanMode, 64)/64
!    kScan    = iand(ScanMode, 32)/32
!
!-----------------------------------------------------------------------
!
!  Space view perspective or orthographic.
!
!-----------------------------------------------------------------------
!
!  ELSEIF (iProj .eq. 90) THEN
!    Gthin = igds(17)

!    Ni = igds(2)
!    Nj = igds(3)
!    Np =  Ni * Nj

!    LatSW = .001 * Float(igds(4))     ! Lat of sub-satellite point
!    LonSW = .001 * Float(igds(5))     ! Lon of sub-satellite point
!    LatNE = Missing
!    LonNE = Missing
!    LonTrue = .001 * Float(igds(13))

!    Di = igds(7)        ! Diameter of the earth in x direction
!    Dj = igds(8)        ! Diameter of the earth in y direction

!    iRes   = iand(igds(6),128)/128
!    iEarth = iand(igds(6), 64)/64
!    iComp  = iand(igds(6),  8)/8

!    ScanMode = igds(11)
!    iScan    = iand(ScanMode,128)/128
!    jScan    = iand(ScanMode, 64)/64
!    kScan    = iand(ScanMode, 32)/32

!    Xp = Float(igds( 9))
!    Yp = Float(igds(10))
!    Xo = Float(igds(15))
!    Yo = Float(igds(16))
!    Zo = Float(igds(14))
!
!-----------------------------------------------------------------------
!
!  Reserved or undefined projection:
!
!-----------------------------------------------------------------------
!
  ELSE

    PRINT *,'WARNING - Undefined projection!! iproj = ', iproj

    RETURN

  END IF
!
!-----------------------------------------------------------------------
!
!  Extract vertical coordinate parameters:
!
!-----------------------------------------------------------------------
!
!  IF (igds(12) .ne. 0) THEN
!    Nk = igds(12)
!    Do i = 1, Nk
!    Zk(i) = psec2(i+10)
!    End Do
!  ELSE
!    Nk = Nint(Missing)
!  ENDIF
!
!-----------------------------------------------------------------------
!
!  Rotated and stretched grids were introduced in Edition 1.
!
!-----------------------------------------------------------------------
!
!  IF (iEdit .ge. 1) THEN
!
!-----------------------------------------------------------------------
!
!  Check for rotated grid information.
!
!-----------------------------------------------------------------------
!
  IF (iproj == 10 .OR. iproj == 30 .OR. iproj == 14 .OR. iproj == 34 .OR. &
        iproj == 60 .OR. iproj == 80) THEN
    latrot = .001 * FLOAT(igds(13))
    lonrot = .001 * FLOAT(igds(14))
!      AngRot = .001 *       psec2( 1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Check for stretched grid information.
!
!-----------------------------------------------------------------------
!
  IF (iproj == 20 .OR. iproj == 30 .OR. iproj == 24 .OR. iproj == 34 .OR. &
        iproj == 70 .OR. iproj == 80) THEN

    latstr = .001 * FLOAT(igds(15))
    lonstr = .001 * FLOAT(igds(16))
!      FacStr = psec2(2)
  END IF
!  ENDIF

  RETURN
END SUBROUTINE grbgrid
