!-----------------------------------------------------------------------
!
!  In the realtime system, ARPS format files arrive after they become
!  available.  To process them as soon as them become available requires
!  that "arpsintrp" wait for files to arrive instead of assuming that they
!  are present.  This routine does all the "magic".
!
!-----------------------------------------------------------------------

SUBROUTINE check_file( s, ntries, sleeptime )
   
  CHARACTER(LEN=*)   :: s
  CHARACTER(LEN=256) :: s2
  LOGICAL            :: iexist
  INTEGER            :: ntries, sleeptime

! Find first white space, if any.

  last = LEN(s)

  DO i=1,last
    IF ( s(i:i) == ' ' ) THEN
       last = i - 1
       EXIT
    ENDIF
  ENDDO

!
! We want the "ready" file.
!

  s2 = s(1:last) // "_ready"

  DO i=1,ntries
    INQUIRE(FILE=s2,EXIST=iexist)
    IF ( iexist ) RETURN
    WRITE(6,*) "Waiting for ",s2
    CALL flush(6)                        ! so we see the message!
    CALL sleep( sleeptime )
  END DO

! It looks like the file isn't going to arrive, so just exit so the rest
! of the programs can run as far as possible.

  WRITE(6,*) 'arpsintrp:  check_file:  time limit exceeded'

!        Even though this is an error, we need to exit 0 due to realtime needs.
!        Cntl_arps considers an exit code of 1 meaning the program failed to
!        produce any useful output.

  CALL EXIT( 0 )
END
