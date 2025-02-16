PROGRAM joinfile
!
!-----------------------------------------------------------------------
!
!  To join together on set of history dumps files produced by the
!  processors of MPP machines with message passing.
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY.
!
!  06/12/1997 (G. Bassett)
!  Created this simpler version from joinfiles which joins only
!  a single set of files.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: length

!  INTEGER nx,ny,nz,nstyps, ireturn
  INTEGER nx,ny,nz,nzsoil,nstyps, ireturn

  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: filename1

  WRITE (*,*) "Enter the filename (base name):"
  READ (*,'(a)') filename
  WRITE (*,*) "Enter nproc_x, nproc_y:"
  READ (*,*) nproc_x, nproc_y

  filename1 = trim(filename)//"_0101"

!  CALL get_dims_from_data(1,filename1,nx,ny,nz,nstyps, ireturn)
  CALL get_dims_from_data(1,filename1,nx,ny,nz,nzsoil,nstyps, ireturn)

  IF (ireturn /= 0) THEN
    WRITE (6,*) 'JOINFILE: Error returned from get_dims_from_data, aborting', &
       ireturn
    STOP 1
  ENDIF

  WRITE (6, *) 'Joining files...'

!  CALL joindumps (filename,nx,ny,nz)
  CALL joindumps (filename,nx,ny,nz,nzsoil)

  WRITE (6, *) 'Done joining files...'

END PROGRAM joinfile
