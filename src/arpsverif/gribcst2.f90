!##################################################################
!##################################################################
!######                                                      ######
!######                  MODULE GRIBDIMS2                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

MODULE gribcst2

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This file defines the GRIB dimension parameters.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Craig A. Mattocks
! 6/1/1995.
!
! MODIFICATION HISTORY:
!
! 6/1/1995, version 1.0 (C. Mattocks)
! Wrote original program.
!
! 12/14/1995 (Yuhe Liu)
! Modified to use ARPS GRIB definitions.
!
! November 1999 (Eric Kemp)
! Modified to remove link to extdims.inc and to prevent allocation
! of arrays dimensioned with nx_ext,ny_ext, and nz_ext.  Also
! added RUC2 GRIB file specifications.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid 87) RUC model:
!
!-----------------------------------------------------------------------

  INTEGER :: ruc87grid              ! Grid identification of RUC
  PARAMETER ( ruc87grid = 87 )   ! RUC Polar Stereo

  INTEGER :: ruc87proj             ! Map projection ID for RUC
  PARAMETER ( ruc87proj = 5 )   ! Polar stereographic projection

  INTEGER :: ruc87nlvt2d      ! Number of RUC level types to be processed
  INTEGER :: ruc87nvs2d       ! Number of RUC 2-D varialbes to be read
  PARAMETER ( ruc87nlvt2d = 1, ruc87nvs2d = 3 )

  INTEGER :: ruc87levs2d(ruc87nlvt2d) ! Level or Layer identification of
                                      ! RUC
  INTEGER :: ruc87var_id2d(ruc87nvs2d,ruc87nlvt2d)

  DATA ruc87levs2d  / 001/       ! Ground surface
  DATA ruc87var_id2d/ 11,  8, 223/
!  11 = Temperature (K)
!   8 = Geometric height (m)
! 223 = Plant canopy surface water (kg/m**2)

  INTEGER :: ruc87nvs3d    ! Number of RUC 3-D varialbes to be read
  INTEGER :: ruc87nlvt3d   ! Number of RUC 3-D level types to be read
  PARAMETER ( ruc87nlvt3d = 2, ruc87nvs3d = 6 )

  INTEGER :: ruc87levs3d(ruc87nlvt3d)  ! Number of RUC 3-D varialbes to be
                                       ! read
  INTEGER :: ruc87var_id3d(ruc87nvs3d, ruc87nlvt3d)

  DATA ruc87levs3d/ 109, 111 /
  DATA ruc87var_id3d/   1,  37,  13, 159,  33,  34,                     &
      85,  86, 255, 255, 255, 255/
!   1 = pressure (Pa)
!  37 = montgomery stream function (m**2/s**2)
!  13 = potential temperature (K)
! 159 = condensation pressure of lifted parcel (Pa)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!  85 = soil temperature (K)
!  86 = soil moisture content (kg/m**2)
! 255 = undefined
!
!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid 211) RUC model (awips grid):
!
!-----------------------------------------------------------------------

  INTEGER :: ruc211grid  ! Grid identification of ruc211
  PARAMETER ( ruc211grid = 211 )   ! RUC Lambert Conformal

  INTEGER :: ruc211proj  ! Map projection ID for ruc211
  PARAMETER ( ruc211proj = 3 )   ! for Lambert Conformal

  INTEGER :: ruc211nlvt2d      ! Number of ruc211 level types to be
                               ! processed
  INTEGER :: ruc211nvs2d       ! Number of ruc211 2-D varialbes to be read
  PARAMETER ( ruc211nvs2d = 4, ruc211nlvt2d = 2 )

  INTEGER :: ruc211levs2d(ruc211nlvt2d) ! Level or Layer identification of
                                        ! RUC
  INTEGER :: ruc211var_id2d(ruc211nvs2d,ruc211nlvt2d)

  DATA ruc211levs2d  / 001, 105/       ! Ground surface
  DATA ruc211var_id2d/ 1,  7, 255, 255,                                 &
      11, 52,  33,  34/
! 1 = Pressure (Pa)
! 7 = Geometric height (m)
! 11= T_2m (C)
! 52= RH_2m
! 33= U_10m
! 34= V_10m

  INTEGER :: ruc211nvs3d       ! Number of RUC 3-D varialbes to be read
  INTEGER :: ruc211nlvt3d      ! Number of RUC 3-D level types to be read
  PARAMETER ( ruc211nvs3d = 5, ruc211nlvt3d= 1 )

  INTEGER :: ruc211levs3d(ruc211nlvt3d)       ! Number of RUC 3-D
                                              ! varialbes to be read
  INTEGER :: ruc211var_id3d(ruc211nvs3d, ruc211nlvt3d)

  DATA ruc211levs3d/ 100 /
  DATA ruc211var_id3d/   7,  11,  52,  33,  34/
!   7 = Pressure surface height (m)
!  11 = Temperature (K)
!  52 = Relative humidity (%)
!  33 = u wind (m/s) (grid relative)
!  34 = v wind (m/s) (grid relative)

!-----------------------------------------------------------------------
!
!  Parameters for NMC GRIB (Grid #212) ETA model (AWIP3D):
!
!-----------------------------------------------------------------------

  INTEGER :: eta212grid              ! Grid identification of ETA
  PARAMETER ( eta212grid = 212 )  ! ETA 40 km Lambert Conformal projection

  INTEGER :: eta212proj             ! Map projection ID for RUC
  PARAMETER ( eta212proj = 3 )   ! Lambert Confornal projection

  INTEGER :: eta212nlvt2d  ! Number of ETA 2-D level types to be processed
  INTEGER :: eta212nvs2d   ! Number of ETA 2-D varialbes   to be processed
!  PARAMETER ( eta212nlvt2d = 1, eta212nvs2d = 6 )
  PARAMETER ( eta212nlvt2d = 3, eta212nvs2d = 6 )

  INTEGER :: eta212levs2d(eta212nlvt2d) ! 2-D Level identification of ETA
  INTEGER :: eta212var_id2d(eta212nvs2d,eta212nlvt2d)

!  DATA eta212levs2d  / 001/       ! At ground surface
  DATA eta212levs2d  / 001,102,105/       ! At ground surface
!  DATA eta212var_id2d/   1,   7,  11, 223,  81,  143/
  DATA eta212var_id2d/   1,   7,  11, 223,  81,  143,     &
                         130, 255, 255, 255, 255, 255,       &
                         11, 52, 33, 34, 255, 255 /

                                 !   1 = Surface pressure (pa)
                                 !   7 = Geopotential height (gpm)
                                 !  11 = Surface temperature (k)
                                 ! 223 = Plant canopy surface water (kg/m**2)
                                 !  81 = Land/sea mask
                                 ! 143 = Categorical snow (yes=1;no=0)
                                 ! 130 = Eta reduced MSLP (pa)
                                 !  11 = 2 m temperature (K)
                                 !  52 = 2 m Relative Humidity (%)
                                 !  33 = 10 m u winds (m/s)
                                 !  34 = 10 m v winds (m/s)

  INTEGER :: eta212nlvt3d        ! Number of ETA 3-D level types to be
                                 ! processed
  INTEGER :: eta212nvs3d         ! Number of ETA 3-D varialbes to be
                                 ! processed

  PARAMETER ( eta212nlvt3d = 2, eta212nvs3d = 6 )

  INTEGER :: eta212levs3d(eta212nlvt3d)   ! 3-D Level identification of
                                          ! ETA

  INTEGER :: eta212var_id3d(eta212nvs3d,eta212nlvt3d)

  DATA eta212levs3d / 100, 112 /     ! sigma vertical coordinates
  DATA eta212var_id3d/  11,  51,  33,  34,   7,  39,                    &
      85, 144, 255, 255, 255, 255/
!  11 = Temperature (K)
!  51 = Specific humidity (kg/kg)
!  33 = u wind (m/s)
!  34 = v wind (m/s)
!   7 = Geopotential height (gpm)
!  39 = Pressure vertical velocity (Pa/s)
!  85 = Soil temperature (K)
! 144 = Volumatric soil moisture (m**3/m**3)
! 255 = undefined

!-----------------------------------------------------------------------
!
!  Parameters for Global Re-Analysis on T62 Gaussian grid:
!
!-----------------------------------------------------------------------

  INTEGER :: reanalt62grid  ! Grid identification of Global Reanalysis
  PARAMETER ( reanalt62grid = 255 )   ! Actually it uses grid 98
                                      ! but 255 shows up in PDS

  INTEGER :: reanalt62proj  ! Map projection ID for Re-Analaysis
  PARAMETER ( reanalt62proj = 4 )   ! for Gaussian latlon grid
  INTEGER :: reanalt62nlvt2d      ! Number of RUCawips level types to be
                                  ! processed
  INTEGER :: reanalt62nvs2d       ! Number of RUCawips 2-D varialbes to be
                                  ! read
  PARAMETER ( reanalt62nvs2d = 2, reanalt62nlvt2d = 1 )

  INTEGER :: reanalt62levs2d(reanalt62nlvt2d) ! Level or Layer
                                              ! identification of RUC
  INTEGER :: reanalt62var_id2d(reanalt62nvs2d,reanalt62nlvt2d)

  DATA reanalt62levs2d  / 1 /       ! Ground surface
  DATA reanalt62var_id2d/ 1,  7/
!   1 = Pressure (Pa)
!   7 = Geometric height (m)

  INTEGER :: reanalt62nvs3d       ! Number of 3-D varialbes to be read
  INTEGER :: reanalt62nlvt3d      ! Number of 3-D level types to be read
  PARAMETER ( reanalt62nvs3d = 4, reanalt62nlvt3d = 1 )

  INTEGER :: reanalt62levs3d(reanalt62nlvt3d)   ! Number of 3-D varialbes
                                                ! to be read
  INTEGER :: reanalt62var_id3d(reanalt62nvs3d, reanalt62nlvt3d)

  DATA reanalt62levs3d/ 107 /
  DATA reanalt62var_id3d/   11,  51,  33,  34/
!  11 = Temperature (K)
!  52 = Specific humidity (%)
!  33 = u wind (m/s) (grid relative)
!  34 = v wind (m/s) (grid relative)


!
!-----------------------------------------------------------------------
!
! Parameters for NCEP RUC Native Coordinate GRIB (Grid 236):
!
!-----------------------------------------------------------------------
!
  INTEGER :: rucn236grid              ! Grid identification of RUC
  PARAMETER ( rucn236grid = 236 )   ! RUC Lambert conformal

  INTEGER :: rucn236proj             ! Map projection ID for RUC
  PARAMETER ( rucn236proj = 3 )   ! Lambert conformal projection

  INTEGER :: rucn236nlvt2d      ! Number of RUC level types to be
                                ! processed
  INTEGER :: rucn236nvs2d       ! Number of RUC 2-D varialbes to be read
!  PARAMETER ( rucn236nlvt2d = 1, rucn236nvs2d = 2 )
  PARAMETER ( rucn236nlvt2d = 2, rucn236nvs2d = 6 )

  INTEGER :: rucn236levs2d(rucn236nlvt2d) ! Level or Layer identification
                                          ! of RUC
  INTEGER :: rucn236var_id2d(rucn236nvs2d,rucn236nlvt2d)

!  DATA rucn236levs2d  / 001/       ! Ground surface
  DATA rucn236levs2d  / 001,102/       ! Ground surface
!  DATA rucn236var_id2d/ 223, 66/
  DATA rucn236var_id2d/ 223, 66, 255, 255, 255, 255,                  &
                        129, 255, 255, 255, 255, 255 /
           ! 223 = Plant canopy surface water (kg/m**2)
           !  66 = Snow depth (m)
           ! 129 = MAPS mean sea level pressure (Pa)
           
  INTEGER :: rucn236nvs3d    ! Number of RUC 3-D varialbes to be read
  INTEGER :: rucn236nlvt3d   ! Number of RUC 3-D level types to be read
!  PARAMETER ( rucn236nlvt3d = 2, rucn236nvs3d = 6 )
  PARAMETER ( rucn236nlvt3d = 2, rucn236nvs3d = 11 )

  INTEGER :: rucn236levs3d(rucn236nlvt3d)  ! Number of RUC 3-D variables
                                            ! to be read
  INTEGER :: rucn236var_id3d(rucn236nvs3d, rucn236nlvt3d)

  DATA rucn236levs3d/ 109, 111 / ! Hybrid coordinates
!  DATA rucn236var_id3d/   1,  7,  189, 185,  33,  34,                   &
!                         85,  144, 255, 255, 255, 255/  
  DATA rucn236var_id3d/   1,  7,  189, 185,  33,  34,   &
                        153,170,178,171,179,            &  ! Hybrid levels
                         85,  144, 255, 255, 255, 255,   & 
                        255,255,255,255,255/               ! Soil levels
           !   7 = height (m)
           ! 189 = virtual potential temperature (K)   
           ! 185 = water vapor mixing ratio (kg/kg)
           !  33 = u wind (m/s)
           !  34 = v wind (m/s)  
           !  85 = soil temperature (K)
           ! 144 = soil volumetric moisture content (fraction)
           ! 255 = undefined

!-----------------------------------------------------------------------
!
!  Parameters for NCEP RUC2 Isobaric GRIB (Grid #236)
!
!-----------------------------------------------------------------------

  INTEGER :: rucp236grid  ! Grid identification of ruc211
  PARAMETER ( rucp236grid = 236 )   ! RUC Lambert Conformal

  INTEGER :: rucp236proj  ! Map projection ID for ruc211
  PARAMETER ( rucp236proj = 3 )   ! for Lambert Conformal

  INTEGER rucp236nlvt2d      ! Number of rucp236 level types to be
                             ! processed
  INTEGER rucp236nvs2d       ! Number of rucp236 2-D varialbes to be
                             ! read
!  PARAMETER ( rucp236nvs2d = 4, rucp236nlvt2d = 2 )
  PARAMETER ( rucp236nvs2d = 4, rucp236nlvt2d = 3 )

  INTEGER :: rucp236levs2d(rucp236nlvt2d) ! Level or Layer identification
                                           ! of RUC
  INTEGER :: rucp236var_id2d(rucp236nvs2d,rucp236nlvt2d)

!  DATA rucp236levs2d  / 001, 105/       ! Ground surface
  DATA rucp236levs2d  / 001, 105, 102/       ! Ground surface
!  DATA rucp236var_id2d/ 1,  7, 66, 255,                               &    
!                       11, 52, 33, 34/
  DATA rucp236var_id2d/ 1,  7, 66, 255,                               &    
                       11, 52, 33, 34,                                &
                       129, 255,255,255/
           ! 1 = Pressure (Pa)
           ! 7 = Geometric height (m)
           ! 66= snow depth (m)
           !255= undefined
           ! 11= T_2m (C)
           ! 52= RH_2m
           ! 33= U_10m
           ! 34= V_10m
           ! 129 = MAPS Mean sea level pressure (Pa)

  INTEGER :: rucp236nvs3d       ! Number of RUC 3-D varialbes to be read
  INTEGER :: rucp236nlvt3d      ! Number of RUC 3-D level types to be
                                ! read
  PARAMETER ( rucp236nvs3d = 5, rucp236nlvt3d= 1 )

  INTEGER rucp236levs3d(rucp236nlvt3d)       ! Number of RUC 3-D
                                             ! varialbes to be read
  INTEGER rucp236var_id3d(rucp236nvs3d, rucp236nlvt3d)

  DATA rucp236levs3d/ 100 /
  DATA rucp236var_id3d/   7,  11,  52,  33,  34/
           !   7 = Pressure surface height (m)
           !  11 = Temperature (K)
           !  52 = Relative humidity (%)
           !  33 = u wind (m/s) (grid relative)
           !  34 = v wind (m/s) (grid relative)

!-----------------------------------------------------------------------
!
!  General definitions for GRIB parameters and variables
!
!-----------------------------------------------------------------------

  INTEGER :: nprods     ! Maximum number of GRIB products (2-D arrays)
  PARAMETER ( nprods = 2000 )

  INTEGER :: n2dvs      ! Maximum number of 2-D variables   to be read
  INTEGER :: n2dlvt     ! Maximum number of 2-D level types to be read
!  PARAMETER ( n2dlvt = 2, n2dvs = 6 )
  PARAMETER ( n2dlvt = 3, n2dvs = 6 )

  INTEGER :: n3dvs      ! Maximum number of 3-D variables   to be read
  INTEGER :: n3dlvt     ! Maximum number of 3-D level types to be read
!  PARAMETER ( n3dlvt = 2, n3dvs = 6 )
  PARAMETER ( n3dlvt = 2, n3dvs = 11 )

  INTEGER :: rcbytes(nprods)  ! record length in bytes
  INTEGER :: rcstart(nprods)  ! record starting byte in a GRIB fil

  INTEGER :: var_nr2d(n2dvs,n2dlvt)
! number of record for all 2-D variables
  INTEGER :: var_nr3d(n3dvs,n3dlvt)  ! number of record for each 3-D
                                     ! variable
!  INTEGER :: var_lev3d(nzgrb,n3dvs,n3dlvt)
! Levels (hybrid) for each 3-D variable

!  REAL :: rcdata(nxygrb)      ! temporary data array
!  REAL :: var_grb2d(nxgrb,nygrb,n2dvs,n2dlvt) ! GRIB variables
!  REAL :: var_grb3d(nxgrb,nygrb,nzgrb,n3dvs,n3dlvt) ! GRIB 3-D variables

!-----------------------------------------------------------------------
!
!  The following variables and arrays must be set before calling
!  subroutine RDNMCGRB.
!
!  (See subroutine GETNMCRUC for RUC and GETNMCETA for ETA in file
!  getextd3d.f as examples.)
!
!-----------------------------------------------------------------------

  INTEGER :: gridtyp    ! Grid ID (working variable)
  INTEGER :: mproj_grb  ! Map projection ID

  INTEGER :: n2dvars    ! Number of 2-D variables (working variable)
  INTEGER :: n2dlvtps   ! Number of 2-D level types (working variable)
  INTEGER :: levtyp2d(n2dlvt)     ! 2-D Level IDs (working variable)
  INTEGER :: var_id2d(n2dvs,n2dlvt)
! Working array for 2-D variable IDs

  INTEGER :: n3dvars    ! Number of 3-D variables (working variable)
  INTEGER :: n3dlvtps   ! Number of 3-D level types (working variable)
  INTEGER :: levtyp3d(n3dlvt)   ! 3-D Level IDs (working variable)
  INTEGER :: var_id3d(n3dvs,n3dlvt)
! Working array for 3-D variable IDs


!-----------------------------------------------------------------------
!
!  GRIB section parameters:
!
!-----------------------------------------------------------------------

  INTEGER :: mxnbufsz               ! Size of GRIB buffer
  INTEGER :: nbufsz                 ! Size of GRIB buffer
  PARAMETER ( mxnbufsz = 800000, nbufsz = 4*(mxnbufsz/4) )

  INTEGER :: npdsz, ngdsz           ! sizes of PDS, GDS
  PARAMETER ( npdsz = 100, ngdsz = 200 )

  INTEGER :: ipdsz, igdsz           ! sizes of IPDS, IGDS
  PARAMETER ( ipdsz = 100, igdsz = 200 )

  INTEGER :: mptrs(15)  ! Integer array containing storage for
                        ! following parameters:
                        ! (1)  total length of grib message
                        ! (2)  length of indicator (SECTION  0)
                        ! (3)  length of pds       (SECTION  1)
                        ! (4)  length of gds       (SECTION  2)
                        ! (5)  length of bms       (SECTION  3)
                        ! (6)  length of bds       (SECTION  4)
                        ! (7)  value of current byte
                        ! (8)  bit pointer
                        ! (9)  GRIB start bit nr
                        ! (10) GRIB/GRID element count
                        ! (11) nr unused bits at end of SECTION 3
                        ! (12) bit map flag (copy of BMS octets 5,6)
                        ! (13) nr unused bits at end of SECTION 2
                        ! (14) BDS flags (right adj copyat end of SECTION 2
                        ! (14) BDS flags (right adj copy of octet 4)
                        ! (15) nr unused bits at end of SECTION 4

  INTEGER :: ipds(ipdsz)       ! PDS integer array
  INTEGER :: igds(igdsz)       ! GDS integer array
!  LOGICAL :: ibms(nxygrb)      ! BMS logical array for data bit map

  INTEGER :: ibdshd(4)         ! IBDS header

  CHARACTER (LEN=1) :: pds(npdsz)    ! PDS ( GRIB Section 1)
  CHARACTER (LEN=1) :: gds(npdsz)    ! GDS ( GRIB Section 2)
  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS ( GRIB Section 4)
  INTEGER :: ibds(nbufsz/4)    ! identical to BDS

!  EQUIVALENCE ( bds,ibds )

  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Buffer to carry GRIB messages

!-----------------------------------------------------------------------
!
!  Parameters which define WMO GRIB Edition 2 format:
!
!  WARNING:  Don't muck with these parameters unless the WMO changes
!            the definition of GRIB Edition 2 format!
!
!-----------------------------------------------------------------------

  INTEGER :: nctrs    ! Maximum number of meteorological centers
  INTEGER :: nsbctrs  ! Maximum number of meteorological sub-centers
  INTEGER :: nmods    ! Maximum number of numerical models/processes
  INTEGER :: ntbvars  ! Maximum number of variables in GRIB tables
  INTEGER :: ntimes   ! Maximum number of time units & range indicators
  INTEGER :: nproj    ! Maximum number of map projection names
  INTEGER :: nscan    ! Maximum number of directional scanning modes

  PARAMETER ( nctrs   = 99  )
  PARAMETER ( nsbctrs = 170 )
  PARAMETER ( nmods   = 153 )
  PARAMETER ( ntbvars = 255 )
  PARAMETER ( ntimes  = 254 )
  PARAMETER ( nproj   = 254 )
  PARAMETER ( nscan   = 8   )

!-----------------------------------------------------------------------
!
!  Table 0 of national/international originating centers, designated
!  by Parameter Description Section (PDS) in Section 1 of GRIB
!  message:
!     Center     = PDS octet # 5
!     Sub-center = PDS octet #26
!
!  NOTE:  (1) Center 1 is a sub-center value for Center 7, the US
!             National Meteorological Center.
!         (2) Centers 150 - 170 are sub-center values for Center 9,
!             the US National Weather Service Field Stations.
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: center(0:nsbctrs)       ! Names of met.
                                                ! sub-centers
  CHARACTER (LEN=42) :: sbcenter(0:nsbctrs)     ! Names of met.
                                                ! sub-centers

!-----------------------------------------------------------------------
!
!  Table A of numerical model or generating process, designated by
!  Parameter Description Section (PDS) octet #6 in Section 1 of GRIB
!  message.
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: model(nmods)         ! Names of numerical models

!-----------------------------------------------------------------------
!
!  Table 2 of GRIB variables, designated by Parameter Description
!  Section (PDS) octet #9 in Section 1 of GRIB message.
!
!  NOTE:  Variables 128 - 254 are RESERVED for use by the originating
!         center. NWS/NMC usage is assumed.
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: varname(0:ntbvars)       ! Names of GRIB variables
  CHARACTER (LEN=42) :: varunit(0:ntbvars)       ! Units of GRIB variables
  CHARACTER (LEN=42) :: varabbr(0:ntbvars)       ! Abbrev of GRIB variables

!-----------------------------------------------------------------------
!
!  Table 3 of levels/layers where variables are defined, designated
!  by Parameter Description Section (PDS) octets #10-12 in Section 1
!  of GRIB message.
!
!  NOTE: The numbering allows for additions within this framework:
!
!        Octet #10   Numerical Precision
!        ---------   -------------------
!          100-119   normal
!          120-139   high
!          140-159   mixed
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: varlvl(0:ntbvars)       ! Levels of GRIB variables

!-----------------------------------------------------------------------
!
!  Tables 4, 5 of time units and range indicators, designated by
!  Parameter Description Section (PDS) in Section 1 of GRIB message:
!     Forecast Time Unit   = PDS octet #18
!     Time Range Indicator = PDS octet #21
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: timeunit(0:ntimes)      ! Forecast time units
  CHARACTER (LEN=42) :: timerang(0:ntimes)      ! Time range indicators

!-----------------------------------------------------------------------
!
!  Table 6 of map projection type, designated by Grid Description
!  Section (GDS) in Section 2 of GRIB message:
!     Projection = GDS octet #8
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: projs(0:nproj)          ! Names of map projections

!-----------------------------------------------------------------------
!
!  Table 8 of directional scanning modes, designated by Grid
!  Description Section (GDS) in Section 2 of GRIB message:
!     Scanning Mode = GDS octet #28
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: scans(0:1,nscan)    ! Directional scanning modes

!-----------------------------------------------------------------------
!
!  End of GRIBCST2.F90
!
!-----------------------------------------------------------------------

END MODULE gribcst2

