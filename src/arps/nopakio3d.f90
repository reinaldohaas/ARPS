
SUBROUTINE pakread
  WRITE (6,'(/a/a,a/)')                                                 &
      'Option -io pak was not used when doing makearps.',               &
      'No packed binary output was produced. ',                         &
      'Re-do makearps with -io pak option.'
  RETURN
END SUBROUTINE pakread

SUBROUTINE pakdump
  WRITE (6,'(/a/a,a/)')                                                 &
      'Option -io pak was not used when doing makearps.',               &
      'No packed binary output was produced. ',                         &
      'Re-do makearps with -io pak option.'
  RETURN
END SUBROUTINE pakdump

SUBROUTINE decdhdr
  WRITE (6,'(/a/a,a/)')                                                 &
      'Option -io pak was not used when doing makearps.',               &
      'No packed binary output was produced. ',                         &
      'Re-do makearps with -io pak option.'
  RETURN
END SUBROUTINE decdhdr
