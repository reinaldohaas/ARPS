!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE V3DDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Coastal Meteorology Research Project  (CMRP)     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE v5ddump
  WRITE (6,'(/a/a,a/)')                                                 &
      'Option -io v5d was not used when doing makearps.',               &
      'No Vis5d output was produced. ',                                 &
      'Re-do makearps with -io v5d option.'
  CALL arpsstop('Vis5d libary was not linked.',1)
  RETURN
END SUBROUTINE v5ddump
