!************************************************************************
!* GEMPRM.PRM
!*
!* This include file contains parameter definitions for the GEMPAK
!* software.
!*
!* OSF1 Version
!**
!* Log:     See gemprm.log
!************************************************************************
!
!     Missing data definitions
!
  REAL, PARAMETER ::  rmissd = -9999.0  
!                          Missing data value
  REAL, PARAMETER ::  rdiffd =  0.1     
!                              Missing value fuzziness
  INTEGER, PARAMETER ::  imissd = -9999    
!                              Missing integer value
  LOGICAL :: ermiss
!                              Declare for stmt func
!
!      Physical and mathematical constants
!
  REAL, PARAMETER ::  pi = 3.14159265   
  REAL, PARAMETER ::  halfpi = pi / 2.  
  REAL, PARAMETER ::  twopi  = 2. * pi  
  REAL, PARAMETER ::  pi4th  = pi / 4.  ! PI,...
  REAL, PARAMETER ::  dtr = pi / 180.   
  REAL, PARAMETER ::  rtd = 180. / pi   ! Degrees <--> Radians
  REAL, PARAMETER ::  radius = 6371200. ! Earth radius
  REAL, PARAMETER ::  omega  = 7.2921E-5! Earth angular veclocity
  REAL, PARAMETER ::  gravty = 9.80616  ! Acceleration of gravity
  REAL, PARAMETER ::  rdgas  = 287.04    
  REAL, PARAMETER ::  rkap   = rdgas / gravty  ! Gas constant of dry air
  REAL, PARAMETER ::  rkappa = 2. / 7.  
  REAL, PARAMETER ::  akappa = 7. / 2.  ! Poisson constant;inverse
  REAL, PARAMETER ::  gamusd = 6.5      ! US std atmos lapse rate
  REAL, PARAMETER ::  tmck   = 273.15   ! Centigrade -> Kelvin
!
!     File information parameters
!
  INTEGER, PARAMETER ::  mmkey  =    12 ! Maximum # of keys
  INTEGER, PARAMETER ::  mmhdrs = 10000 ! Maximum # of headers
  INTEGER, PARAMETER ::  mmprt  =    20 ! Maximum # of parts
  INTEGER, PARAMETER ::  mmlist =    20 ! Maximum search list
  INTEGER, PARAMETER ::  mmfree =    62 ! Number of free pairs
  INTEGER, PARAMETER ::  mmfile =     3 ! Maximum # of open files
  INTEGER, PARAMETER ::  mblksz =   128 ! Block size
  INTEGER, PARAMETER ::  mcache =     8 ! # of cached records
  INTEGER, PARAMETER ::  mmparm =    40 ! Maximum # of parameters
  INTEGER, PARAMETER ::  mmfhdr =    10 ! Maximum # of file hdrs
  INTEGER, PARAMETER ::  mmsrch =    30 ! Max # of cond searches
!
  INTEGER, PARAMETER ::  mmfldp = mmfile * mmprt  
!
  INTEGER, PARAMETER ::  mtvax  =     2    
  INTEGER, PARAMETER ::  mtsun  =     3    
  INTEGER, PARAMETER ::  mtiris =     4    
  INTEGER, PARAMETER ::  mtapol =     5    
  INTEGER, PARAMETER ::  mtibm  =     6    
  INTEGER, PARAMETER ::  mtigph =     7    
  INTEGER, PARAMETER ::  mtultx =     8    
  INTEGER, PARAMETER ::  mthp   =     9    
  INTEGER, PARAMETER ::  mtalph =    10    
  INTEGER, PARAMETER ::  mtmach =  mtalph  
!                      Machine type
!                              Raob, VAS data source
  INTEGER, PARAMETER ::  mmrecl =     1    
!                              Multiplier for RECL in
!                          file create/open
!                              (usually 4 on UNIX sys)
!
  INTEGER, PARAMETER ::  mdreal =     1    
  INTEGER, PARAMETER ::  mdintg =     2    
  INTEGER, PARAMETER ::  mdchar =     3    
  INTEGER, PARAMETER ::  mdrpck =     4    
  INTEGER, PARAMETER ::  mdgrid =     5    
!                              Data types in DM lib
  INTEGER, PARAMETER ::  mdgnon =     0    
  INTEGER, PARAMETER ::  mdggrb =     1    
  INTEGER, PARAMETER ::  mdgnmc =     2    
  INTEGER, PARAMETER ::  mdgdif =     3    
  INTEGER, PARAMETER ::  mdgdec =     4    
!                          Grid packing types
  INTEGER, PARAMETER ::  mfsf =       1    
  INTEGER, PARAMETER ::  mfsn =       2    
  INTEGER, PARAMETER ::  mfgd =       3    
!                          Data file types
  INTEGER, PARAMETER ::  mfunkn =    99    
  INTEGER, PARAMETER ::  mfairw  =    1    
  INTEGER, PARAMETER ::  mfmetr =     2    
  INTEGER, PARAMETER ::  mfship =     3    
!                              Unknown, airways, metar,
!                          ship data source
  INTEGER, PARAMETER ::  mfbuoy =     4    
  INTEGER, PARAMETER ::  mfsynp =     5    
  INTEGER, PARAMETER ::  mfraob =     4    
  INTEGER, PARAMETER ::  mfvas  =     5    
!
  INTEGER, PARAMETER ::  mfgrid =     6    
!                                              Grid source
  INTEGER, PARAMETER ::  mftext =   100    
!                                              Text
!
!     Declarations for array sizes in programs
!
  INTEGER, PARAMETER ::  llmxlv =   500    ! Max # levels/station
  INTEGER, PARAMETER ::  llmxtm =   200    ! Max # times/dataset
  INTEGER, PARAMETER ::  llmxgt =  1000    ! Max # grid times
  INTEGER, PARAMETER ::  llmxst =    20    ! Max # stations in list
  INTEGER, PARAMETER ::  llmxdt = mmparm * llmxlv  ! Max # data points
  INTEGER, PARAMETER ::  llmxpt = 10000    ! Max # lines, polygons, ...
  INTEGER, PARAMETER ::  llstfl =  9800    ! Max # stations in file
  INTEGER, PARAMETER ::  llmxgd =  97000   ! Max # grid points
  INTEGER, PARAMETER ::  llmdgg = 2304000  ! Max mem for intern grids
  INTEGER, PARAMETER ::  mxloop =    30    ! Maximum frames in loop
  INTEGER, PARAMETER ::  llnnav =   256    ! Grid nav block length
  INTEGER, PARAMETER ::  llnanl =   128    ! Grid anl block length
  INTEGER, PARAMETER ::  llsthl =    20    ! Max header size
  INTEGER, PARAMETER ::  llgdhd =   128    ! Max grid hdr length
  INTEGER, PARAMETER ::  lloagd =   400    ! Max # grids from 1 OA
  INTEGER, PARAMETER ::  llclev =   50     ! Max # of contour lvls
  INTEGER, PARAMETER ::  llaxis =   64     ! Max # of axis labels
  INTEGER, PARAMETER ::  lltmcx =  100     ! Max # of time xsect pts
!
!     Image file types (defined here and in gemprm.h)
!
  INTEGER, PARAMETER ::  ifinvd = -1         
  INTEGER, PARAMETER ::  ifarea =  1         
  INTEGER, PARAMETER ::  ifgini =  2          
  INTEGER, PARAMETER ::  ifnids =  3         
  INTEGER, PARAMETER ::  ifnowr =  4         
!                          File type NONE,
!                      AREA (MCIDAS),
!                          GINI (AWIPS), NIDS
!                          radar, NOWRad (WSI)
!
!     Offsets from local to UTC (GMT) time in HHMM (hour/minute) form
!
  INTEGER, PARAMETER ::  joflst =   500     ! Offset for UTC/EST
  INTEGER, PARAMETER ::  jofldt =   400     ! Offset for UTC/EDT
!
! ASCII character constants
!
  CHARACTER (LEN=1), PARAMETER :: chnull = CHAR (0)   ! Null
  CHARACTER (LEN=1), PARAMETER :: chctla = CHAR (1)   ! Control A
  CHARACTER (LEN=1), PARAMETER :: chctlc = CHAR (3)   ! Control C
  CHARACTER (LEN=1), PARAMETER :: chtab  = CHAR (9)   ! Tab
  CHARACTER (LEN=1), PARAMETER :: chlf   = CHAR (10)  ! Line feed
  CHARACTER (LEN=1), PARAMETER :: chff   = CHAR (12)  ! Form feed
  CHARACTER (LEN=1), PARAMETER :: chcr   = CHAR (13)  ! Carriage return
  CHARACTER (LEN=1), PARAMETER :: chcan  = CHAR (24)  ! Cancel (CAN)
  CHARACTER (LEN=1), PARAMETER :: chesc  = CHAR (27)  ! Escape
  CHARACTER (LEN=1), PARAMETER :: chfs   = CHAR (28)  ! FS
  CHARACTER (LEN=1), PARAMETER :: chgs   = CHAR (29)  ! GS
  CHARACTER (LEN=1), PARAMETER :: chrs   = CHAR (30)  ! Record Separator
  CHARACTER (LEN=1), PARAMETER :: chus   = CHAR (31)  ! US
  CHARACTER (LEN=1), PARAMETER :: chspac = CHAR (32)  ! Space
  CHARACTER (LEN=1), PARAMETER :: chtlda = CHAR (126) ! Tilda
!
!----------------------- COORDINATE SYSTEMS -----------------------------

  CHARACTER (LEN=10), PARAMETER ::  sysup = 'DNVPLWMQIG'
  CHARACTER (LEN=10), PARAMETER ::  syslo = 'dnvplwmqig'
  CHARACTER (LEN=1) :: carray (10)

  COMMON          / gpsys / carray
!
!-----------------------GPLT SYMBOL SIZE--------------------------------
  REAL, PARAMETER ::  radclm = .13   ! Calm wind circle radius
  REAL, PARAMETER ::  radsky = 3.5   ! Sky cover symbol radius
  REAL, PARAMETER ::  rszptn = 9.0   ! Pressure tendency symbol size
!
!-----------------------GPLT BUFFER SIZE--------------------------------

  INTEGER, PARAMETER :: igbszm  = 20 + ( 4 * llmxgd )  
!                Min desired GPLT buffer size
  INTEGER, PARAMETER :: igbsiz = (128*(1+INT((igbszm-1)/128))) 
!                Makes buffer size mul. of 128
  INTEGER, PARAMETER :: igtbsz  = 1100
!
!-----------------------DEVICE BUFFER SIZE------------------------------
  INTEGER, PARAMETER :: igdszm = 20 + ( 2 * llmxpt )
!                         Min desired DEVICE buffer size
!
  INTEGER, PARAMETER :: igdsiz = (128*(1+INT((igdszm-1)/128)))
!                         Makes buffer size mul. of 128
!
!-----------------------CYLINDRICAL MAP TRANSFORMATIONS-----------------
  INTEGER, PARAMETER ::  mccyl = 1  
  INTEGER, PARAMETER ::  mpcequ = 1, mpcmer = 2, mpcmcd = 3  
  INTEGER, PARAMETER ::  mscequ = 1  
!
!-----------------------AZIMUTHAL MAP TRANSFORMATIONS-------------------
  INTEGER, PARAMETER ::  mcazm = 2  
  INTEGER, PARAMETER ::  mpaequ = 1, mpastr = 2, mpaort = 3, mpalam = 4  
  INTEGER, PARAMETER ::  mpagno = 5  
  INTEGER, PARAMETER ::  msanor = 1, msasou = 2  
!
!-----------------------CONICAL MAP TRANSFORMATIONS---------------------
  INTEGER, PARAMETER ::  mccon = 3  
  INTEGER, PARAMETER ::  mpcnor = 1, mpcsou = 2  
!
!-----------------------OBLIQUE MERCATOR--------------------------------
  INTEGER, PARAMETER ::  mcmer = 4  
  INTEGER, PARAMETER ::  mptmer = 1, mputm = 2, mpoblq = 3  
!
!-----------------------SATELLITE TRANSFORMATIONS-----------------------
  INTEGER, PARAMETER ::  mcgoes = 6  
  INTEGER, PARAMETER ::  mpmci = 1  
!
!-----------------------COLOR TABLE INFORMATION-------------------------
  INTEGER, PARAMETER ::  mxclnm = 50  
!                              Maximum number of color names saved
!
!-----------------------DEVICE CHARACTERISTICS--------------------------
  INTEGER, PARAMETER ::  ndvchr = 20 
!                    Size of /DEVCHR/ to be passed
