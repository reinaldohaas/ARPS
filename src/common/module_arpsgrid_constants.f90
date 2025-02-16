MODULE module_arpsgrid_constants
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains only constants to be used in program ARPS4WRF.
! and NMM2ARPS.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  INTEGER, PARAMETER :: ROOT = 0        ! Root processor

  INTEGER, PARAMETER :: STDIN  = 5
  INTEGER, PARAMETER :: STDOUT = 6
  INTEGER, PARAMETER :: STDERR = 0

  INTEGER, PARAMETER :: MAXHIS = 100    ! Maximum number of history files

  INTEGER, PARAMETER :: MAXFILELEN = 256

  REAL,    PARAMETER :: MISSING_DATA = HUGE(1.0)
  REAL,    PARAMETER :: ARPS_MISSING = -999.0

  REAL, PARAMETER :: p0 = 1.0E5    ! Surface reference pressure, is 100000 Pascal.

  REAL, PARAMETER :: rd = 287.0    ! Gas constant for dry air  (m**2/(s**2*K))
  REAL, PARAMETER :: cp = 1004.0   ! Specific heat of dry air at constant pressure
                                   ! (m**2/(s**2*K)).
  REAL, PARAMETER :: cv = 717.0    ! Specific heat of dry air at constant volume
                                   ! (m**2/(s**2*K)).
  REAL, PARAMETER :: rv = 461.0    ! Gas constant for water vapor  (m**2/(s**2*K)).

  REAL, PARAMETER :: rddcp  = rd/cp
  REAL, PARAMETER :: cpdcv  = cp/cv
  REAL, PARAMETER :: rvdrd  = rv/rd
  REAL, PARAMETER :: rddrv  = rd/rv

  REAL, PARAMETER :: g = 9.81
  REAL, PARAMETER :: govrd = g/rd
  REAL, PARAMETER :: rdovg = rd/g

  REAL, PARAMETER :: gamma = 0.0065   ! degrees/m lapse rate
  REAL, PARAMETER :: const = ( rd*gamma/g )

END MODULE module_arpsgrid_constants
