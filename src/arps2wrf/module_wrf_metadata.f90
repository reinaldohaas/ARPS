MODULE wrf_metadata

!  INTEGER, PARAMETER :: wrfversion = 2   ! Denote the supported WRF version,
                                ! Should be specified in arps2wrf.input
                                ! = 1  WRFV1.3
                                ! = 2  WRFV2.0.3
  REAL, PARAMETER  :: a_small_real_number = 1.0E-5

  REAL, PARAMETER  :: r_d   = 287.0
  REAL, PARAMETER  :: r_v   = 461.6
  REAL, PARAMETER  :: rvovrd= r_v/r_d
  REAL, PARAMETER  :: cp_wrf= 7.*r_d/2.
  REAL, PARAMETER  :: cv_wrf= cp_wrf - r_d
  REAL, PARAMETER  :: cvpm  = -cv_wrf/cp_wrf
  REAL, PARAMETER  :: g_wrf = 9.81

  REAL, PARAMETER  :: t0      = 300.0
  REAL, PARAMETER  :: p1000mb = 1.0E+5


  REAL, PARAMETER  :: omega_ear = 7.292e-5
  REAL(KIND(0D0)), PARAMETER  :: d2rfactor = 3.14159265358979D0/180.

  INTEGER, PARAMETER  :: LanduseCategories = 24
  INTEGER, PARAMETER  :: SoilCategories    = 16
  INTEGER, PARAMETER  :: ISWATER      = 16
  INTEGER, PARAMETER  :: ISICE        = 24
  INTEGER, PARAMETER  :: ISURBAN      = 1
  INTEGER, PARAMETER  :: ISWATER_SOIL = 14

  ! constants for MMINLU = "UMD"
  INTEGER, PARAMETER  :: ISWATER_UMD  = 14
  INTEGER, PARAMETER  :: ISICE_UMD    = 14
  INTEGER, PARAMETER  :: ISURBAN_UMD  = 13

  INTEGER, PARAMETER  :: WRF_REAL                 = 104
  INTEGER, PARAMETER  :: WRF_DOUBLE               = 105
  INTEGER, PARAMETER  :: WRF_INTEGER              = 106
  INTEGER, PARAMETER  :: WRF_LOGICAL              = 107
  INTEGER, PARAMETER  :: WRF_COMPLEX              = 108
  INTEGER, PARAMETER  :: WRF_DOUBLE_COMPLEX       = 109

  ! Added in WRFV2.1 to initialize WRF, Passed in from namelist
  REAL     :: base_temp  = 290.
  REAL     :: base_pres  = 100000.
  REAL     :: base_lapse = 50.
  REAL     :: iso_temp   = 0.

!   This module defines global attributes for the wrfinput_d01 file

  TYPE wrf_global_metadata
    CHARACTER(LEN=80)   :: title
    CHARACTER(LEN=24)   :: start_date
    INTEGER             :: we_dimension    ! staggered value
    INTEGER             :: sn_dimension
    INTEGER             :: bt_dimension    ! non-staggered value

    INTEGER             :: dyn_opt
    INTEGER             :: diff_opt
    INTEGER             :: km_opt
    INTEGER             :: damp_opt
    REAL                :: dampcoef
    REAL                :: khdif
    REAL                :: kvdif
    INTEGER             :: mp_physics
    INTEGER             :: ra_lw_physics
    INTEGER             :: ra_sw_physics
    INTEGER             :: sf_sfclay_physics
    INTEGER             :: sf_surface_physics
    INTEGER             :: bl_pbl_physics
    INTEGER             :: cu_physics

    INTEGER             :: we_p_unstag_s
    INTEGER             :: we_p_unstag_e
    INTEGER             :: we_p_stag_s
    INTEGER             :: we_p_stag_e
    INTEGER             :: sn_p_unstag_s
    INTEGER             :: sn_p_unstag_e
    INTEGER             :: sn_p_stag_s
    INTEGER             :: sn_p_stag_e
    INTEGER             :: bt_p_unstag_s
    INTEGER             :: bt_p_unstag_e
    INTEGER             :: bt_p_stag_s
    INTEGER             :: bt_p_stag_e

    REAL                :: dx
    REAL                :: dy
    REAL                :: dt
    REAL                :: cen_lat
    REAL                :: cen_lon
    REAL                :: tru_lat1
    REAL                :: tru_lat2
    REAL                :: stand_lon
    REAL                :: moad_cen_lat

    REAL                :: gmt
    INTEGER             :: julyr
    INTEGER             :: julday

    INTEGER             :: num_land_cat
    INTEGER             :: islake
    INTEGER             :: iswater
    INTEGER             :: isice
    INTEGER             :: isurban
    INTEGER             :: isoilwater
    INTEGER             :: map_proj
    CHARACTER(LEN=4)    :: mminlu

    INTEGER             :: grid_id         ! Added in WRFV2.1
    INTEGER             :: parent_id
    INTEGER             :: i_parent_start
    INTEGER             :: j_parent_start
    INTEGER             :: parent_grid_ratio

                                          ! New since WRF2.2
    INTEGER    :: surface_input_source    ! 1: static (fractional),
                                          ! 2: time dependent (dominant)
                                          ! 3: hybrid (not yet implemented)
    INTEGER    :: sst_update              ! update sst from wrflowinp file  0=no, 1=yes
    INTEGER    :: grid_fdda
    INTEGER    :: gfdda_interval_m
    INTEGER    :: gfdda_end_h

    INTEGER    :: grid_sfdda              ! Added since WRFV3.1, since we do not use FFDA, they are ignored.
    INTEGER    :: sgfdda_interval_m
    INTEGER    :: sgfdda_end_h

  END TYPE wrf_global_metadata

  ! this type defines the variable attributes

  TYPE wrf_var_metadata
    CHARACTER(LEN=12)   :: name
    INTEGER             :: fieldType
    CHARACTER(LEN=3)    :: memoryOrder
    CHARACTER(LEN=80)   :: description
    CHARACTER(LEN=25)   :: units
    CHARACTER(LEN=1)    :: stagger
    CHARACTER(LEN=80)   :: dimName1
    CHARACTER(LEN=80)   :: dimName2
    CHARACTER(LEN=80)   :: dimName3
    CHARACTER(LEN=80)   :: coordinates
  END TYPE wrf_var_metadata

  !
  ! IO format constants
  !
  INTEGER, PARAMETER :: IO_NET   = 7        ! In ARPS system, including arps2wrf
  INTEGER, PARAMETER :: IO_INT   = 1
  INTEGER, PARAMETER :: IO_PHDF5 = 5

  INTEGER, PARAMETER :: IO_NET_WRF   = 2    ! WRF system
  INTEGER, PARAMETER :: IO_INT_WRF   = 1
  INTEGER, PARAMETER :: IO_PHDF5_WRF = 4

  !
  ! Microphysics variables
  !

  INTEGER :: P_QV_wrf, P_QC_wrf, P_QR_wrf
  INTEGER :: P_QI_wrf, P_QS_wrf, P_QG_wrf
  INTEGER :: P_QT_wrf, P_QNI_wrf
  INTEGER :: num_moist, num_scalar, nscalar_wrf
  INTEGER :: arps_Q_PTR(20)
  CHARACTER(LEN=8)  :: qnames_wrf(20)
  CHARACTER(LEN=40) :: qdescp_wrf(20)

END MODULE wrf_metadata

!
!###################### WRF subroutines ##############################
!
FUNCTION projrot_latlon(iproj,trulat1,trulat2,trulon,ctrlat,ctrlon,     &
                        rlat,rlon,istatus)
!
! Adapted from WRFSI
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: iproj
  REAL,    INTENT(IN)  :: trulat1,trulat2,trulon
  REAL,    INTENT(IN)  :: ctrlat,ctrlon
  REAL,    INTENT(IN)  :: rlat,rlon
  INTEGER, INTENT(OUT) :: istatus

  REAL :: projrot_latlon

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  REAL :: angdif

  REAL :: s,n
  REAL :: colat1,colat2
  REAL :: rn

  REAL, PARAMETER :: d2rfactor = 3.1415926/180.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(ABS(iproj) == 1) THEN ! polar stereographic

    IF(trulat1 == +90.)then
      projrot_latlon = trulon - rlon
    ELSE IF(trulat1 == -90.)then
      projrot_latlon = rlon - trulon
    ELSE  ! abs(trulat1) .ne. 90.
      IF(ctrlat ==  trulat1 .AND. ctrlon == trulon) THEN
                                     ! grid centered on proj pole
        rn = COS( (90.-trulat1) * d2rfactor)
        projrot_latlon = rn * angdif(trulon,rlon)

      ELSE IF(.TRUE.) THEN
        write(6,*)' ERROR in projrot_latlon: '
        write(6,*)' This type of local',                         &
                  ' stereographic projection not yet supported.'
        write(6,*)' Grid should be centered on projection pole.'

        projrot_latlon = 0.

      ELSE ! .false.
        ! Find dx/lat and dy/lat, then determine projrot_laps

      END IF

    END IF ! trulat1

  ELSE IF(ABS(iproj) == 2) THEN ! lambert conformal

    IF(trulat1 >= 0.0) THEN
      s = +1.
    ELSE
      s = -1.
    END IF

    colat1 = 90. - s * trulat1
    colat2 = 90. - s * trulat2

    IF(trulat1 == trulat2) THEN  ! tangent lambert
      n      = COS( d2rfactor * colat1 )
    ELSE                     ! two standard latitudes
      n = ALOG(COS(d2rfactor * trulat1) / COS( d2rfactor*trulat2) )/     &
          ALOG(TAN(d2rfactor * (45.-s*trulat1/2.) ) /                    &
               TAN(d2rfactor * (45.-s*trulat2/2.) ) )
    END IF

    projrot_latlon = n * s * angdif(trulon,rlon)

  ELSE IF(ABS(iproj) == 3) THEN ! mercator
    projrot_latlon = 0.

  ELSE
    WRITE(6,*) 'projrot_latlon: unrecognized projection ',iproj
    STOP
  END IF

  istatus = 1
  RETURN
END FUNCTION projrot_latlon
!
!##################################################################
!
FUNCTION angdif(x,y)
!
! Difference between two angles, result is between -180. and +180.
!
  IMPLICIT NONE

  REAL, INTENT(IN) :: x,y
  REAL             :: angdif

  angdif = MOD(X-Y+540.,360.)-180.

  RETURN
END FUNCTION angdif
!
!##################################################################
!
SUBROUTINE const_module_initialize (p00,t00,a)

!-----------------------------------------------------------------
!
! PURPOSE:
!
!   Initialize p00,t00 & a from ARPS2WRF namelist input.
!
!----------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  REAL, INTENT(OUT)  :: p00
  REAL, INTENT(OUT)  :: t00
  REAL, INTENT(OUT)  :: a

  p00 = base_pres
  t00 = base_temp
  a   = base_lapse

!  p00 = 100000.   ! constant sea level pressure, Pa
!  t00 = 290.      ! constant sea level temperature, K
!  a   = 50.       ! temperature difference from 1000mb to 300mb, K
END SUBROUTINE

!REAL(KIND(0D0)) FUNCTION d2rfactor
!    d2rfactor = 2.d0*ACOS(0.d0)/180.d0
!END FUNCTION d2rfactor
