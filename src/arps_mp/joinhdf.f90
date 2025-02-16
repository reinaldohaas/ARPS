PROGRAM joinhdf
!
!-----------------------------------------------------------------------
!
!  To join together a set of ARPS history or data files produced by the
!  processors of MPP machines with message passing.
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY.
!
!  2001/04/23 (G. Bassett) Created.
!  wdt Copyright (c) 2001 Weather Decision Technologies, Inc.
!
!  2003/07/17 (Yunheng Wang)
!  Modified to work with the new 4D soil variables which were
!  added since ARPS version IHOP_3.
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

  INTEGER :: nxsm,nysm,nz,nzsoil
  INTEGER :: nxlg,nylg
  INTEGER :: nstyps, ireturn

  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: filename1

  REAL,    ALLOCATABLE :: buf_r(:,:,:), buf_rsm(:,:,:)
  REAL,    ALLOCATABLE :: buf_rsoil4d(:,:,:,:), buf_rsoil4dsm(:,:,:,:)
  REAL,    ALLOCATABLE :: buf_r1(:), buf_r2(:)
  INTEGER, ALLOCATABLE :: buf_i(:,:,:), buf_ism(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE ::                  &
                          buf_i16(:,:,:), buf_i16sm(:,:,:)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE ::                  &
                          buf_i16soil(:,:,:,:), buf_i16soilsm(:,:,:,:)
  INTEGER :: sstat

  LOGICAL :: fexist
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE (*,*) "Enter the filename (base name):"
  READ (*,'(a)') filename
  WRITE (*,*) "Enter nproc_x, nproc_y:"
  READ (*,*) nproc_x, nproc_y

  filename1 = trim(filename)//"_0101"

  INQUIRE(FILE=TRIM(filename1),EXIST=fexist)
  IF (.NOT. fexist) filename1 = TRIM(filename)//"_001001"

  CALL get_dims_from_hdf(filename1,nxsm,nysm,nz,nzsoil,nstyps, ireturn)

  IF (ireturn /= 0) THEN
    WRITE (6,*) 'JOINFILE: WARNING, error returned from get_dims_from_data', &
       ireturn
  ENDIF

  nxlg = (nxsm-3)*nproc_x+3
  nylg = (nysm-3)*nproc_y+3

  ALLOCATE(buf_r  (nxlg,nylg,nz))
  ALLOCATE(buf_rsm(nxsm,nysm,nz))
  ALLOCATE(buf_rsoil4d  (nxlg,nylg,nzsoil,0:nstyps))
  ALLOCATE(buf_rsoil4dsm(nxsm,nysm,nzsoil,0:nstyps))
  ALLOCATE(buf_r1(nxsm+nysm+nz))
  ALLOCATE(buf_r2(nxlg+nylg+nz))
  ALLOCATE(buf_i  (nxlg,nylg,nz))
  ALLOCATE(buf_ism(nxsm,nysm,nz))
  ALLOCATE(buf_i16  (nxlg,nylg,nz))
  ALLOCATE(buf_i16sm(nxsm,nysm,nz))
  ALLOCATE(buf_i16soil  (nxlg,nylg,nzsoil,0:nstyps))
  ALLOCATE(buf_i16soilsm(nxsm,nysm,nzsoil,0:nstyps))

  WRITE (6, *) 'Joining files ...'

  CALL join_hdf (filename,nxsm,nysm,nz,nzsoil,nstyps,nxlg,nylg,       &
                 buf_r,buf_rsm,buf_rsoil4d,buf_rsoil4dsm,             &
                 buf_i,buf_ism,buf_i16,buf_i16sm,                     &
                 buf_i16soil,buf_i16soilsm,                           &
                 buf_r1,buf_r2,sstat)

  WRITE (6, *) 'Done joining files ...'

  DEALLOCATE(buf_r)
  DEALLOCATE(buf_rsm)
  DEALLOCATE(buf_rsoil4d)
  DEALLOCATE(buf_rsoil4dsm)
  DEALLOCATE(buf_r1)
  DEALLOCATE(buf_r2)
  DEALLOCATE(buf_i)
  DEALLOCATE(buf_ism)
  DEALLOCATE(buf_i16)
  DEALLOCATE(buf_i16sm)
  DEALLOCATE(buf_i16soil)
  DEALLOCATE(buf_i16soilsm)

END PROGRAM joinhdf
