!########################################################################
!########################################################################
!#########                                                      #########
!#########                     nidscst.inc                      #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

!------------------------------------------------------------------------
! 
! PURPOSE:
! 
! Declare constants for use with NIDS2ARPS and related NIDS decoding
! routines.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 9 August 2001
! Based on constants listed in C version of NIDS2ARPS, and HPUX 
! FORTRAN 77 software developed by the NWS Techniques Development
! Laboratory.
!
! Eric Kemp, 17 August 2001
! Removed INTEGER*4, INTEGER*1, etc, and avoided use of hexadecimal
! parameters (based on work by Yunheng Wang).  Also, placed module in
! separate file.  Finally, added dimension for array irad in subroutine
! get_radial_data.
!
! Keith Brewster, 20 January 2010
! Added Product IDs and resolutions for TDWR data
!
! Keith Brewster, 17 March 2010
! Converted to an include file to be consistent with all other 
! ARPS programs setting parameters and constants.
!
!------------------------------------------------------------------------
  
!------------------------------------------------------------------------
!
! Julian date of 1 Jan 1970.  Used to convert NEXRAD date to yr/mon/day.
! 
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: jul_1_1_70=2440587
  
!------------------------------------------------------------------------
!
! Internal ID for radial products.
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: ipradial1 = 44831 ! Z'AF1F'
  INTEGER, PARAMETER :: ipradial2 =  5503 ! Z'157F'
  INTEGER, PARAMETER :: ipradial3 = 34175 ! Z'857F'
  INTEGER, PARAMETER :: ipradial4 = 50559 ! Z'c57F'
  INTEGER, PARAMETER :: ipradial5 = 12927 ! Z'327F'

!------------------------------------------------------------------------
! 
! Internal ID's for raster products.
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: ipraster1 = 47631 ! Z'BA0F'
  INTEGER, PARAMETER :: ipraster2 = 47623 ! Z'BA07'

!------------------------------------------------------------------------
!
! Product ID codes.
!  
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: BREFprodID1 = 16  ! Base reflectivity ID code
  INTEGER, PARAMETER :: BREFprodID2 = 21  ! 16 <= ID <= 21
  INTEGER, PARAMETER :: TDWRrefID1  = 181 ! TDWR reflectivity
  INTEGER, PARAMETER :: TDWRrefID2  = 186 ! TDWR long-range reflectivity

  INTEGER, PARAMETER :: BVELprodID1 = 22  ! Base velocity ID code
  INTEGER, PARAMETER :: BVELprodID2 = 27  ! 22 <= ID <+ 27
  INTEGER, PARAMETER :: TDWRvelID   = 182 ! TDWR long-range reflectivity

  INTEGER, PARAMETER :: ETprodID = 41    ! Echo Top ID code

  INTEGER, PARAMETER :: VILprodID = 57   ! VIL ID code

  INTEGER, PARAMETER :: DPAprodID = 81   ! DPA ID code

  INTEGER, PARAMETER :: CREFprodID = 38  ! Composite Refl ID code

!------------------------------------------------------------------------
!
! Unit conversions
!  
!------------------------------------------------------------------------

  REAL, PARAMETER :: rad2deg = 180./3.141592654
  REAL, PARAMETER :: deg2rad = 1./rad2deg 
  REAL, PARAMETER :: mm2in = 0.001*39.37
  REAL, PARAMETER :: mpfoot = 0.3048
  REAL, PARAMETER :: kts2ms = 0.514444

!------------------------------------------------------------------------
!
! Map projection information for NIDS DPA data.
!  
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: DPAiproj = 1 ! Polar stereographic
  REAL, PARAMETER :: DPAtrulat1 = 60.
  REAL, PARAMETER :: DPAtrulon = -105.
  REAL, PARAMETER :: DPAscale = 1.
  REAL, PARAMETER :: DPAdx = 4762.5  ! Grid resolution (m) at trulon
 
!------------------------------------------------------------------------
!
! Maximum dimensions
!  
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: maxlen = 100000    ! max. length of NIDS product

  INTEGER, PARAMETER :: maxelev = 4
  INTEGER, PARAMETER :: maxgate = 600
  INTEGER, PARAMETER :: maxazim = 464
  INTEGER, PARAMETER :: maxgate_tilt = maxgate*maxazim

  INTEGER, PARAMETER :: ndim1 = 600  ! Maximum dimensions of the decoded
  INTEGER, PARAMETER :: ndim2 = 600  ! product (rows/cols or
                                     ! aximuths/ranges)

  INTEGER, PARAMETER :: iraddim = 600

!NEW EMK 13 September 2001
  INTEGER, PARAMETER :: maxnidsfile = 11 ! four tilts reflectivity,
                                         ! four tilts velocity,
                                         ! one VIL, one echo top, and
                                         ! one DPA

!------------------------------------------------------------------------
!
! Character constants.
!  
!------------------------------------------------------------------------

  CHARACTER (LEN=6)  :: refid    = 'refl3d'
  CHARACTER (LEN=20) :: refname  = 'Reflectivity'
  CHARACTER (LEN=20) :: refunits = 'dBZ'

  CHARACTER (LEN=6)  :: velid    = 'radv3d'
  CHARACTER (LEN=20) :: velname  = 'Radial Velocity'
  CHARACTER (LEN=20) :: velunits = 'm/s'

!------------------------------------------------------------------------
!
! Miscellaneous parameters.
!  
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: iscale5 = 40962 ! Z'A002'
   
  INTEGER, PARAMETER :: stdin = 5, stdout=6

  REAL, PARAMETER :: misschk = -199.0
  INTEGER, PARAMETER :: r_missing = -999.

  REAL, PARAMETER :: res(0:200) = (/                             &
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &! 0-9  
               0., 0., 0., 0., 0., 0., 1., 2., 4., 1.,           &! 10-19
               2., 4., 0.25, 0.5, 1., 0.25, 0.5, 1., 0., 0.,     &! 20-29
               0., 0., 0., 0., 0., 1., 4., 1., 4., 0.,           &! 30-39
               0., 4., 0., 0., 0., 0., 0., 0., 0., 0.,           &! 40-49
               0., 0., 0., 0., 0., 0., 0., 4., 0., 0.,           &! 50-59
               0., 0., 0., 4., 4., 4., 4., 0., 0., 0.,           &! 60-69
               0., 0., 0., 0., 0., 0., 0., 0., 2., 2.,           &! 70-79
               2., 0., 0., 0., 0., 0., 0., 0., 0., 4.,           &! 80-89
               4., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &! 90-99
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!100-109
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!110-119
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!120-129
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!130-139
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!140-149
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!150-159
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!160-169
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           &!170-179
               0., 0.15, 0.15, 0., 0., 0., 0.3, 0., 0., 0.,      &!180-189
               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)       !190-200

  REAL, PARAMETER :: rasterchek = -9999.
  REAL, PARAMETER :: rastermiss = -9999.

  REAL, PARAMETER :: dpachek = -9999.
  REAL, PARAMETER :: dpamiss = -9999.
