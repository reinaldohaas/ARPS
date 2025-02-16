!
! --------------------------------------------------------------------
!

SUBROUTINE wrlcm (arr,apr,nw,npack)
!
!  this is a stub for non-wordpacking applications.
!  if this is called we have a problem
!
  WRITE(6,'(''  IN PACKING ROUTINE WITH NPACK = '',I4)')                &
       npack
!
  WRITE(6,'('' ERROR EXIT, NO PACKING SHOULD BE ON '')')
  STOP
  RETURN
END SUBROUTINE wrlcm
!
! --------------------------------------------------------------------
!

SUBROUTINE rdlcm (arr,apr,nw,npack)
!
!  this is a stub for non-wordpacking applications.
!  if this is called we have a problem
!
  WRITE(6,'(''  IN UNPACKING ROUTINE WITH NPACK = '',I4)')              &
       npack
!
  WRITE(6,'('' ERROR EXIT, NO UNPACKING SHOULD BE ON '')')
  STOP
  RETURN
END SUBROUTINE rdlcm
