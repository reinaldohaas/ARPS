!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETLAPS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getlaps
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy subroutine of GETLAPS (see getlaps.f).
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE (6,'(a/a)')                                                     &
      'Dummy subroutine was called. Please check option of extdopt ',   &
      'that led to this call. Program stopped in nolaps.f.'
  STOP

END SUBROUTINE getlaps
