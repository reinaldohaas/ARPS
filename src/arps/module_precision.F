
MODULE arps_precision

  ! Set precision for floating point variables.

  IMPLICIT NONE

  INTEGER, PARAMETER :: INT16 =  SELECTED_INT_KIND(4)

  ! Choose one to run in single or double precision.
#ifdef DOUBLE
  INTEGER, PARAMETER  :: P  = KIND(1.0D0)         ! Double precision
#else
  INTEGER, PARAMETER  :: P  = KIND(1.0)           ! Single precision
#endif

  ! In some places such as I/O, locally override precision P.
  INTEGER, PARAMETER  :: SP = SELECTED_REAL_KIND(6, 30)
  INTEGER, PARAMETER  :: DP = SELECTED_REAL_KIND(11,30)

  REAL(P), PARAMETER :: UNDEF = 9.99000000000E+33  ! Undefined value

CONTAINS

  SUBROUTINE print_precision

    WRITE(6,'(4(1x,a,I2,a,/))')                                         &
              'ARPS default precision is:  P = ',p, ',',                &
              '     Single  precision is: SP = ',sp,',',                &
              '     Double  precision is: DP = ',dp,'.',                &
              '     Integer kind      is: INT16 = ',INT16,'.'
    RETURN
  END SUBROUTINE print_precision

END MODULE arps_precision
