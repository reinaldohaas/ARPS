
PROGRAM joinfiles
!
!-----------------------------------------------------------------------
!
!  To join together history dumps files produced by the processors
!  of MPP machines with message passing.
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY.
!
!  11/06/1995 (M. Xue)
!  Set the start time for file joining to zero instead of tstart.
!  tstart may not be at the history dump time for a restart run.
!  The program will skip the times when the corresponding files
!  are not found.
!
!  2003/07/17 (Y. Wang)
!  Modified to work with the new 4D soil variables which were
!  added since ARPS version IHOP_3.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz,nzsoil,nstyps
  INTEGER :: nxsm,nysm

  INTEGER :: i

  CHARACTER (LEN=256) :: tmplnth  ! Temporary array to store namelist logname

  INTEGER :: ndumps, time
  CHARACTER (LEN=80) :: next

  REAL,    ALLOCATABLE :: buf_r(:,:,:), buf_rsm(:,:,:)
  REAL,    ALLOCATABLE :: buf_rsoil4d(:,:,:,:), buf_rsmsoil4d(:,:,:,:)
  REAL,    ALLOCATABLE :: buf_r1(:), buf_r2(:)
  INTEGER, ALLOCATABLE :: buf_i(:,:,:), buf_ism(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE ::                      &
                          buf_i16(:,:,:), buf_i16sm(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE ::                      &
                          buf_i16soil(:,:,:,:), buf_i16soilsm(:,:,:,:)
  INTEGER :: sstat
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'

  INCLUDE 'globcst.inc'  !  Global constants and parameters, most of
                         !  them specify the model run options.

  INCLUDE 'bndry.inc'    !  Control parameters defining the boundary
                         !  condition types.

  INCLUDE 'phycst.inc'   !  Universal physical constants such as gas constants.

  INCLUDE 'exbc.inc'     !  External boundary parameters and variables.

!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=256) :: namelist_filename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Read the input file to check what input options, which need input
!  data files have been set.
!
!-----------------------------------------------------------------------
!
  namelist_filename = ' '
  CALL initpara(nx, ny, nz, nzsoil, nstyps, namelist_filename)
  nproc_x = nproc_x_in
  nproc_y = nproc_y_in

  ! Convert to processor nx & ny (initpara thinks we're in non-MP mode)

  IF (nx /= nproc_x*int((nx-3)/nproc_x)+3) THEN
    nx = nproc_x*int((nx-3)/nproc_x+0.9999999999999) + 3
    IF (myproc == 0) THEN
      WRITE (6,*) "WARNING: adjusting nx to fit on ",nproc_x," processors:"
      WRITE(6,'(5x,a,i5)') "   new nx =",nx
    ENDIF
  ENDIF
  IF (ny /= nproc_y*int((ny-3)/nproc_y)+3) THEN
    ny = nproc_y*int((ny-3)/nproc_y+0.9999999999999) + 3
    IF (myproc == 0) THEN
      WRITE (6,*) "WARNING: adjusting ny to fit on ",nproc_y," processors:"
      WRITE(6,'(5x,a,i5)') "   new ny =",ny
    ENDIF
  ENDIF

  nxsm = (nx - 3)/nproc_x + 3
  nysm = (ny - 3)/nproc_y + 3

!  IF( thisdmp == 0.0 ) THEN
  IF( thisdmp <= 0.1 ) THEN
    PRINT*,                                                             &
        'The history dump option was off. No file joining was done.'
    STOP
  END IF

  ndumps = nint((tstop-tstart)/thisdmp) + 1

  ALLOCATE(buf_r(nx,ny,nz))
  ALLOCATE(buf_rsm(nxsm,nysm,nz))
  ALLOCATE(buf_rsoil4d(nx,ny,nzsoil,0:nstyps))
  ALLOCATE(buf_rsmsoil4d(nxsm,nysm,nzsoil,0:nstyps))
  ALLOCATE(buf_r1(nxsm+nysm+nz))
  ALLOCATE(buf_r2(nx+ny+nz))
  ALLOCATE(buf_i(nx,ny,nz))
  ALLOCATE(buf_ism(nxsm,nysm,nz))
  ALLOCATE(buf_i16(nx,ny,nz))
  ALLOCATE(buf_i16sm(nxsm,nysm,nz))
  ALLOCATE(buf_i16soil(nx,ny,nzsoil,0:nstyps))
  ALLOCATE(buf_i16soilsm(nxsm,nysm,nzsoil,0:nstyps))
!
!-----------------------------------------------------------------------
!
!  Join the base state data dump
!
!-----------------------------------------------------------------------

  IF (hdmpfmt == 1) THEN
    tmplnth = dirname(1:ldirnam)//'/'//runname(1:lfnkey)//'.bingrdbas'
    CALL joindumps (tmplnth,nxsm,nysm,nz,nzsoil)
  ELSE IF (hdmpfmt == 3) THEN
    tmplnth = dirname(1:ldirnam)//'/'//runname(1:lfnkey)//'.hdfgrdbas'
    CALL join_hdf(tmplnth,nxsm,nysm,nz,nzsoil,nstyps,nx,ny,           &
                  buf_r,buf_rsm, buf_rsoil4d,buf_rsmsoil4d,           &
                  buf_i,buf_ism,buf_i16,buf_i16sm,                    &
                  buf_i16soil,buf_i16soilsm,                          &
                  buf_r1,buf_r2,sstat)
  ELSE
    WRITE (6,*) "History dumps not in compatible format for joining."
    STOP
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Join the history dump files
!
!-----------------------------------------------------------------------
!
  time = INT(tstart)

  DO i = 1, ndumps

    WRITE (next, '(i6.6)') time

    IF (hdmpfmt == 1) THEN
      tmplnth = dirname(1:ldirnam)//'/'//                               &
                runname(1:lfnkey)//'.bin'//next
      CALL joindumps (tmplnth,nxsm,nysm,nz,nzsoil)
    ELSE IF (hdmpfmt == 3) THEN
      tmplnth = dirname(1:ldirnam)//'/'//                               &
                runname(1:lfnkey)//'.hdf'//next
      CALL join_hdf(tmplnth,nxsm,nysm,nz,nzsoil,nstyps,nx,ny,           &
                    buf_r,buf_rsm, buf_rsoil4d,buf_rsmsoil4d,           &
                    buf_i,buf_ism,buf_i16,buf_i16sm,                    &
                    buf_i16soil,buf_i16soilsm,                          &
                    buf_r1,buf_r2,sstat)
    ELSE
      WRITE (6,*) "History dumps not in compatible format for joining."
      STOP
    ENDIF

    time = time + INT (thisdmp)

  END DO

  WRITE (6, *) 'Done joining files...'

  DEALLOCATE(buf_r)
  DEALLOCATE(buf_rsm)
  DEALLOCATE(buf_rsoil4d)
  DEALLOCATE(buf_rsmsoil4d)
  DEALLOCATE(buf_r1)
  DEALLOCATE(buf_r2)
  DEALLOCATE(buf_i)
  DEALLOCATE(buf_ism)
  DEALLOCATE(buf_i16)
  DEALLOCATE(buf_i16soil)
  DEALLOCATE(buf_i16soilsm)

END PROGRAM joinfiles
