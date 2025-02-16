MODULE module_wrfgrid_constants
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains only constants to be used in program ARPS4WRF
! and NMM2ARPS
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  INTEGER, PARAMETER :: BINARY = 1
  INTEGER, PARAMETER :: NETCDF = 2

  INTEGER, PARAMETER :: MAXDOM = 100    ! Maximum number of nested domains
  INTEGER, PARAMETER :: MAXFLAGS = 100

  INTEGER, PARAMETER :: MAXDIMENSIONS = 3
  INTEGER, PARAMETER :: MAXFIELDS = 256

  integer, parameter :: WRF_REAL    = 104
  integer, parameter :: WRF_INTEGER = 106

  integer, parameter :: M=1, U=2, V=3, HH=4, VV=5
  integer, parameter :: W=6, SS = 7
  integer, parameter :: OUTSIDE_DOMAIN=1E8, NOT_PROCESSED=1E9, INVALID=1E9

  ! Projection codes for proj_info structure:
  INTEGER, PUBLIC, PARAMETER  :: PROJ_LATLON = 0
  INTEGER, PUBLIC, PARAMETER  :: PROJ_LC = 1
  INTEGER, PUBLIC, PARAMETER  :: PROJ_PS = 2
  INTEGER, PUBLIC, PARAMETER  :: PROJ_PS_WGS84 = 102
  INTEGER, PUBLIC, PARAMETER  :: PROJ_MERC = 3
  INTEGER, PUBLIC, PARAMETER  :: PROJ_GAUSS = 4
  INTEGER, PUBLIC, PARAMETER  :: PROJ_CYL = 5
  INTEGER, PUBLIC, PARAMETER  :: PROJ_CASSINI = 6
  INTEGER, PUBLIC, PARAMETER  :: PROJ_ALBERS_NAD83 = 105
  INTEGER, PUBLIC, PARAMETER  :: PROJ_ROTLL = 203

  ! Interpolation method
  INTEGER, PARAMETER :: NEARNEIGHBOR  = 1
  INTEGER, PARAMETER :: NEARNEIGHBOR1 = 11
  INTEGER, PARAMETER :: NEARNEIGHBOR2 = 21
  INTEGER, PARAMETER :: FOUR_PNT      = 2
  INTEGER, PARAMETER :: SIXTEEN_PNT   = 3
  INTEGER, PARAMETER :: BILINEAR      = 4
  INTEGER, PARAMETER :: QUADRATIC     = 5
  INTEGER, PARAMETER :: ASSIGN_DIRECT = 9999
  INTEGER, PARAMETER :: EXTRACT_SUBDOMAIN = 9998

  ! NMM2ARPS specific
  INTEGER, PARAMETER :: lvlprof = 601
  REAL,    PARAMETER :: depthp  = 3.0E4

  CHARACTER(LEN=128) :: dimwename  = 'west_east'
  CHARACTER(LEN=128) :: dimsnname  = 'south_north'
  CHARACTER(LEN=128) :: dimbtname  = 'bottom_top'
  CHARACTER(LEN=128) :: dimbtsname = 'bottom_top_stag'

END MODULE module_wrfgrid_constants
