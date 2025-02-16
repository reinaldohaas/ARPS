!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PAKREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pakread
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A dummy routine to be called instead of the truly functional routine
!  pakread when word packing routines are not available on a given system.
!
!-----------------------------------------------------------------------
!
  WRITE (6,'(5(/a))')                                                   &
      ' We have stopped supporting the packed binrary format',          &
      ' starting from ARPS version 5.0.0Beta4, for the lack of use.',   &
      ' If you wish to use this format, you can obtain the file ',      &
      ' pakio3d.f from the arps5.0.0Beta3. Some code update',           &
      ' may be needed, however.'
  CALL arpsstop('arpsstop called from PAKREAD not supported.',1)

END SUBROUTINE pakread
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PAKDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pakdump
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    A dummy routine to be called instead of the functional routine
!    pakdump when packed binary formmat is not avaiable.
!
!-----------------------------------------------------------------------
!
  WRITE (6,'(5(/a))')                                                   &
      ' We have stopped supporting the packed binrary format',          &
      ' starting from ARPS version 5.0.0Beta4, for the lack of use.',   &
      ' If you wish to use this format, you can obtain the file ',      &
      ' pakio3d.f from the arps5.0.0Beta3. Some code update',           &
      ' may be needed, however.'
  CALL arpsstop('arpsstop called from PAKDUMP not supported.',1)

END SUBROUTINE pakdump

SUBROUTINE decdhdr
  WRITE (6,'(/a/a,a/)')                                                 &
      'Option -io pak was not used when doing makearps.',               &
      'No packed binary output was produced. ',                         &
      'Re-do makearps with -io pak option.'
  RETURN
END SUBROUTINE decdhdr

