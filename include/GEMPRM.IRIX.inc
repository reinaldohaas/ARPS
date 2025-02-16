!************************************************************************
!* GEMPRM.PRM
!*
!* This include file contains parameter definitions for the GEMPAK
!* software.
!*
!* IRIX Version
!**
!* Log:     See gemprm.log
!************************************************************************
!
!     Missing data definitions
!
  PARAMETER     ( rmissd = -9999.0 )
!                          Missing data value
  PARAMETER     ( rdiffd =  0.1    )
!                              Missing value fuzziness
  PARAMETER     ( imissd = -9999   )
!                              Missing integer value
  LOGICAL :: ermiss
!                              Declare for stmt func
!
!      Physical and mathematical constants
!
  PARAMETER     ( pi = 3.14159265  )
  PARAMETER     ( halfpi = pi / 2. )
  PARAMETER     ( twopi  = 2. * pi )
  PARAMETER     ( pi4th  = pi / 4. )
!                  PI,...
  PARAMETER     ( dtr = pi / 180.  )
  PARAMETER     ( rtd = 180. / pi  )
!                              Degrees <--> Radians
  PARAMETER     ( radius = 6371200. )
!                      Earth radius
  PARAMETER     ( omega  = 7.2921E-5 )
!                              Earth angular veclocity
  PARAMETER     ( gravty = 9.80616  )
!                              Acceleration of gravity
  PARAMETER     ( rdgas  = 287.04   )
  PARAMETER     ( rkap   = rdgas / gravty )
!                              Gas constant of dry air
  PARAMETER     ( rkappa = 2. / 7. )
  PARAMETER     ( akappa = 7. / 2. )
!                              Poisson constant;inverse
  PARAMETER     ( gamusd = 6.5 )
!                              US std atmos lapse rate
  PARAMETER     ( tmck   = 273.15 )
!                              Centigrade -> Kelvin
!
!     File information parameters
!
  PARAMETER     ( mmkey  =    12   )
!                          Maximum # of keys
  PARAMETER     ( mmhdrs =  10000  )
!                              Maximum # of headers
  PARAMETER     ( mmprt  =    20   )
!                          Maximum # of parts
  PARAMETER     ( mmlist =    20   )
!                          Maximum search list
  PARAMETER     ( mmfree =    62   )
!                              Number of free pairs
  PARAMETER     ( mmfile =     3   )
!                              Maximum # of open files
  PARAMETER     ( mblksz =   128   )
!                      Block size
  PARAMETER     ( mcache =     8   )
!                          # of cached records
  PARAMETER     ( mmparm =    40   )
!                              Maximum # of parameters
  PARAMETER     ( mmfhdr =    10   )
!                              Maximum # of file hdrs
  PARAMETER     ( mmsrch =    30   )
!                              Max # of cond searches
!
  PARAMETER     ( mmfldp = mmfile * mmprt )
!
  PARAMETER     ( mtvax  =     2   )
  PARAMETER     ( mtsun  =     3   )
  PARAMETER     ( mtiris =     4   )
  PARAMETER     ( mtapol =     5   )
  PARAMETER     ( mtibm  =     6   )
  PARAMETER     ( mtigph =     7   )
  PARAMETER     ( mtultx =     8   )
  PARAMETER       ( mthp   =     9   )
  PARAMETER     ( mtalph =    10   )
  PARAMETER     ( mtmach = mtiris  )
!                      Machine type
!                              Raob, VAS data source
  PARAMETER     ( mmrecl =     1   )
!                              Multiplier for RECL in
!                          file create/open
!                              (usually 4 on UNIX sys)
!
  PARAMETER     ( mdreal =     1   )
  PARAMETER     ( mdintg =     2   )
  PARAMETER     ( mdchar =     3   )
  PARAMETER     ( mdrpck =     4   )
  PARAMETER     ( mdgrid =     5   )
!                              Data types in DM lib
  PARAMETER     ( mdgnon =     0   )
  PARAMETER     ( mdggrb =     1   )
  PARAMETER     ( mdgnmc =     2   )
  PARAMETER     ( mdgdif =     3   )
  PARAMETER     ( mdgdec =     4   )
!                          Grid packing types
  PARAMETER     ( mfsf =       1   )
  PARAMETER     ( mfsn =       2   )
  PARAMETER     ( mfgd =       3   )
!                          Data file types
  PARAMETER     ( mfunkn =    99   )
  PARAMETER     ( mfairw  =    1   )
  PARAMETER     ( mfmetr =     2   )
  PARAMETER     ( mfship =     3   )
!                              Unknown, airways, metar,
!                          ship data source
  PARAMETER      ( mfbuoy =     4   )
  PARAMETER     ( mfsynp =     5   )
  PARAMETER     ( mfraob =     4   )
  PARAMETER     ( mfvas  =     5   )
!
  PARAMETER       ( mfgrid =     6   )
!                                              Grid source
  PARAMETER       ( mftext =   100   )
!                                              Text
!
!     Declarations for array sizes in programs
!
  PARAMETER     ( llmxlv =   500   )
!                              Max # levels/station
  PARAMETER     ( llmxtm =   200   )
!                          Max # times/dataset
  PARAMETER     ( llmxgt =  1000   )
!                          Max # grid times
  PARAMETER     ( llmxst =    20   )
!                              Max # stations in list
  PARAMETER     ( llmxdt = mmparm * llmxlv )
!                          Max # data points
  PARAMETER     ( llmxpt = 10000   )
!                              Max # lines, polygons, ...
  PARAMETER     ( llstfl =  9800   )
!                              Max # stations in file
  PARAMETER     ( llmxgd =  97000  )
!                          Max # grid points
  PARAMETER     ( llmdgg = 2304000 )
!                              Max mem for intern grids
  PARAMETER       ( mxloop =    30   )
!                                              Maximum frames in loop
  PARAMETER       ( llnnav =   256   )
!                                              Grid nav block length
  PARAMETER       ( llnanl =   128   )
!                                              Grid anl block length
  PARAMETER     ( llsthl =    20   )
!                          Max header size
  PARAMETER     ( llgdhd =   128   )
!                          Max grid hdr length
  PARAMETER     ( lloagd =   400   )
!                              Max # grids from 1 OA
  PARAMETER       ( llclev =   50    )
!                              Max # of contour lvls
  PARAMETER       ( llaxis =   64    )
!                              Max # of axis labels
  PARAMETER     ( lltmcx =  100    )
!                              Max # of time xsect pts
!
!     Image file types (defined here and in gemprm.h)
!
  PARAMETER     ( ifinvd = -1        )
  PARAMETER     ( ifarea =  1        )
  PARAMETER     ( ifgini =  2         )
  PARAMETER     ( ifnids =  3        )
  PARAMETER     ( ifnowr =  4        )
!                          File type NONE,
!                      AREA (MCIDAS),
!                          GINI (AWIPS), NIDS
!                          radar, NOWRad (WSI)
!
!     Offsets from local to UTC (GMT) time in HHMM (hour/minute) form
!
  PARAMETER     ( joflst =   500   )
!                          Offset for UTC/EST
  PARAMETER     ( jofldt =   400   )
!                          Offset for UTC/EDT
!
!     ASCII character constants
!
  CHARACTER (LEN=1) :: chnull, chctla, chctlc, chtab,  chlf,   chff
  CHARACTER (LEN=1) :: chcr,   chcan,  chesc,  chfs,   chgs,   chrs
  CHARACTER (LEN=1) :: chus,   chspac, chtlda
!
  PARAMETER       ( chnull = CHAR (0) )
!                                              Null
  PARAMETER       ( chctla = CHAR (1) )
!                                              Control A
  PARAMETER     ( chctlc = CHAR (3) )
!                  Control C
  PARAMETER     ( chtab  = CHAR (9) )
!              Tab
  PARAMETER     ( chlf   = CHAR (10) )
!                  Line feed
  PARAMETER     ( chff   = CHAR (12) )
!                  Form feed
  PARAMETER     ( chcr   = CHAR (13) )
!                          Carriage return
  PARAMETER     ( chcan  = CHAR (24) )
!                      Cancel (CAN)
  PARAMETER     ( chesc  = CHAR (27) )
!                  Escape
  PARAMETER     ( chfs   = CHAR (28) )
!              FS
  PARAMETER     ( chgs   = CHAR (29) )
!              GS
  PARAMETER     ( chrs   = CHAR (30) )
!                          Record Separator
  PARAMETER     ( chus   = CHAR (31) )
!              US
  PARAMETER     ( chspac = CHAR (32) )
!                  Space
  PARAMETER     ( chtlda = CHAR (126) )
!                  Tilda
!
!----------------------- COORDINATE SYSTEMS -----------------------------
  CHARACTER (LEN=10) :: syslo
  CHARACTER (LEN=10) :: sysup
  PARAMETER       ( sysup = 'DNVPLWMQIG', syslo = 'dnvplwmqig' )
  CHARACTER (LEN=1) :: carray (10)
  COMMON          / gpsys / carray
!
!-----------------------GPLT SYMBOL SIZE--------------------------------
  PARAMETER     ( radclm = .13 )
!                         Calm wind circle radius
  PARAMETER     ( radsky = 3.5 )
!                         Sky cover symbol radius
  PARAMETER       ( rszptn = 9.0 )
!                                      Pressure tendency symbol size
!
!-----------------------GPLT BUFFER SIZE--------------------------------
  PARAMETER     ( igbszm  = 20 + ( 4 * llmxgd ) )
!                         Min desired GPLT buffer size
!
  PARAMETER     ( igbsiz = (128*(1+INT((igbszm-1)/128))))
!                         Makes buffer size mul. of 128
!
  PARAMETER     ( igtbsz  = 1100   )
!
!-----------------------DEVICE BUFFER SIZE------------------------------
  PARAMETER     ( igdszm = 20 + ( 2 * llmxpt ) )
!                         Min desired DEVICE buffer size
!
  PARAMETER     ( igdsiz = (128*(1+INT((igdszm-1)/128))))
!                         Makes buffer size mul. of 128
!
!-----------------------CYLINDRICAL MAP TRANSFORMATIONS-----------------
  PARAMETER      ( mccyl = 1 )
  PARAMETER      ( mpcequ = 1, mpcmer = 2, mpcmcd = 3 )
  PARAMETER     ( mscequ = 1 )
!
!-----------------------AZIMUTHAL MAP TRANSFORMATIONS-------------------
  PARAMETER      ( mcazm = 2 )
  PARAMETER      ( mpaequ = 1, mpastr = 2, mpaort = 3,                  &
                      mpalam = 4 )
  PARAMETER      ( mpagno = 5 )
  PARAMETER     ( msanor = 1, msasou = 2 )
!
!-----------------------CONICAL MAP TRANSFORMATIONS---------------------
  PARAMETER      ( mccon = 3 )
  PARAMETER     ( mpcnor = 1, mpcsou = 2 )
!
!-----------------------OBLIQUE MERCATOR--------------------------------
  PARAMETER       ( mcmer = 4 )
  PARAMETER       ( mptmer = 1, mputm = 2, mpoblq = 3 )
!
!-----------------------SATELLITE TRANSFORMATIONS-----------------------
  PARAMETER      ( mcgoes = 6 )
  PARAMETER      ( mpmci = 1 )
!
!-----------------------COLOR TABLE INFORMATION-------------------------
  PARAMETER       ( mxclnm = 50 )
!                              Maximum number of color names saved
!
!-----------------------DEVICE CHARACTERISTICS--------------------------
  PARAMETER     ( ndvchr = 20 )
!                    Size of /DEVCHR/ to be passed
