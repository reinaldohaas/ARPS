PROGRAM splithdf

!
!-----------------------------------------------------------------------
!
!  To split apart a ARPS history or data file for use by the
!  message passing version of the ARPS.
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY.
!
!  2001/04/23 (G. Bassett) Created.
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
  INTEGER nxsm,nysm
  INTEGER nxlg,nylg,nz,nzsoil
  INTEGER nstyps, ireturn
  INTEGER :: nprocx_in, nprocy_in

  CHARACTER (LEN=256) :: filename

  REAL,    ALLOCATABLE :: buf_r(:,:,:)
  REAL,    ALLOCATABLE :: buf_rsoil(:,:,:,:)
  INTEGER, ALLOCATABLE :: buf_i(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE:: buf_i16(:,:,:,:)

  INTEGER sstat
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE (*,*) "Enter the filename:"
  READ  (*,'(a)') filename
  WRITE (*,*) "Enter nproc_x, nproc_y:"
  READ  (*,*) nprocx_in, nprocy_in

  CALL get_dims_from_hdf(trim(filename),nxlg,nylg,nz,nzsoil,nstyps,ireturn)

  IF (ireturn /= 0) THEN
    WRITE (6,*) 'SPLITHDF: WARNING, error returned from get_dims_from_data', &
       ireturn
  ENDIF

  nxsm = (nxlg - 3)/nprocx_in + 3
  nysm = (nylg - 3)/nprocy_in + 3

  ALLOCATE(buf_r    (nxsm,nysm,nz),              stat = sstat)
  ALLOCATE(buf_rsoil(nxsm,nysm,nzsoil,0:nstyps), stat = sstat)
  ALLOCATE(buf_i    (nxsm,nysm,nz),              stat = sstat)
  ALLOCATE(buf_i16  (nxsm,nysm,MAX(nz,nzsoil),0:nstyps), stat = sstat)
  buf_r     = 0.0
  buf_rsoil = 0.0
  buf_i     = 0
  buf_i16   = 0

  WRITE (6, *) 'Splitting file ...'

  call mpinit_proc(0)
  call mpinit_var()

  nproc_x = nprocx_in
  nproc_y = nprocy_in

  CALL split_hdf(filename,nxsm,nysm,nz,nzsoil,nstyps, &
                 buf_r,buf_rsoil,buf_i,buf_i16,sstat)

  WRITE (6, *) 'Done splitting file ...'

  DEALLOCATE(buf_r)
  DEALLOCATE(buf_rsoil)
  DEALLOCATE(buf_i)
  DEALLOCATE(buf_i16)

  CALL arpsstop("Normal Finish", 0)
END PROGRAM splithdf
