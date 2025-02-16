!
! This file contains needed WRF subroutine for WRF external
! I/O package, such as PHDF5, binary etc.
!
! These subroutines have been rewritten to remove the dependence
! on WRF framework.
!
! Author: Yunheng Wang (01/07/2005).
!
SUBROUTINE wrf_debug(level, str)

  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: str
  INTEGER ,      INTENT(IN) :: level

  IF (level > 0) WRITE(0,*) str

  RETURN
END SUBROUTINE wrf_debug

SUBROUTINE wrf_message(str)

  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: str

  WRITE(6,*) 'WRF MESSAGE: ',str

  RETURN
END SUBROUTINE wrf_message

SUBROUTINE wrf_error_fatal( str )
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: str

  WRITE(0,*) '-------------- FATAL CALLED ---------------' 
  WRITE(0,*) str
  WRITE(0,*) '-------------------------------------------'
  
  CALL arpsstop('ARPSSTOP called.',1)

END SUBROUTINE wrf_error_fatal
!
! This subroutine uses an extended intrinsic procedure
! sizeof(A). Some compilers, such as PGF90 may not provide this
! intrinsic function.
!
SUBROUTINE wrf_sizeof_integer( retval )
  IMPLICIT NONE
  INTEGER retval

#ifndef NOSIZEOF
  retval = sizeof(retval)       ! works on IBM XL Fortran IFORT_SIZEOF
                                ! Intel IFORT Fortran
                                ! PGI Fortran compiler
                                ! OSF1 Fortran compiler
#else
  retval = 4                    ! for PGF90 etc. PGF_SIZEOF
#endif

  RETURN
END SUBROUTINE wrf_sizeof_integer

SUBROUTINE wrf_sizeof_real( retval )
  IMPLICIT NONE
  INTEGER, INTENT(OUT) ::  retval
  REAL   :: tmp 

#ifndef NOSIZEOF
  retval = sizeof(tmp)            ! IFORT_SIZEOF
#else
  retval = 4                     ! for PGF90 etc. PGF_SIZEOF
#endif

  RETURN
END SUBROUTINE wrf_sizeof_real
