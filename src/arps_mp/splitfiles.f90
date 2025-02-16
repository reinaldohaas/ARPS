
PROGRAM splitfiles
!
!-----------------------------------------------------------------------
!
!  Variable Declarations. (Local Variables)
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
  INTEGER :: nx,ny,nz,nstyps,nzsoil
  INTEGER :: nztmp
  INTEGER nxsm,nysm

  INTEGER :: exbcvarsz     ! add by wyh

  INTEGER :: i

  INTEGER :: istat       ! Flag set by open statement on the status
                         ! of file opening
                         ! existence

  INTEGER :: numfiles

  INTEGER :: lfname
  CHARACTER (LEN=256) :: filename

  REAL, ALLOCATABLE :: buf_r(:,:,:)
  REAL, ALLOCATABLE :: buf_rsoil4d(:,:,:,:)
  INTEGER, ALLOCATABLE :: buf_i(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE:: buf_i16(:,:,:,:)
  INTEGER sstat

  integer mp_opt_save

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Global constants and parameters, most of them specify the
!  model run options.
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
!  Control parameters defining the boundary condition types.
!
!-----------------------------------------------------------------------

  INCLUDE 'bndry.inc'

!-----------------------------------------------------------------------
!
!  Universal physical constants such as gas constants.
!
!-----------------------------------------------------------------------

  INCLUDE 'phycst.inc'

!-----------------------------------------------------------------------
!
!  External boundary parameters and variables.
!
!-----------------------------------------------------------------------

  INCLUDE 'exbc.inc'

!-----------------------------------------------------------------------
!
!  Nudging
!
!-----------------------------------------------------------------------

  INCLUDE 'nudging.inc'

!-----------------------------------------------------------------------
!
!  Message passing variables.
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!

  CHARACTER(LEN=256) :: namelist_filename

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Read the input file to check what input options, which need input
!  data files have been set.
!
!-----------------------------------------------------------------------

  namelist_filename = ' '
  CALL initpara(nx, ny, nz, nzsoil,nstyps, namelist_filename)
  nproc_x = nproc_x_in
  nproc_y = nproc_y_in

! Undo what initpara() did so we don't do this twice and get the wrong answers.

  if ( mp_opt > 0 ) then
    nx = ( nx - 3 ) * nproc_x + 3;
    ny = ( ny - 3 ) * nproc_y + 3;
  end if

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

  IF ((inifmt == 3) .or. (ternfmt == 3) .or. (exbcfmt == 3) .or.  &
     (soilfmt == 3) .or. (sfcfmt == 3)) THEN  ! HDF format used
    ALLOCATE(buf_r(nxsm,nysm,nz))
    ALLOCATE(buf_rsoil4d(nxsm,nysm,nzsoil,0:nstyps))
    ALLOCATE(buf_i(nxsm,nysm,nz))
!   ALLOCATE(buf_i16(nxsm,nysm,nzsoil,0:nstyps))
!
!   We need to allocated a 4D array so that "tsoil" and "qsoil" will fit.
!   Although "nzsoil" is the third dimension there, it isn't big enough
!   for other variables.  Historically, "nz" has been used, however, if
!   "nzsoil" is larger than "nz", that isn't good either.
!
!   All of this can be controlled.  What can't be controlled is the size
!   used by the "em" files.  It is controlled by "nstyp" in the "EM" input
!   file, not the "nstyp" in the input to "splitfiles".  The "EM" value
!   of "nstyp" has always been 9, so this may never be a problem.
!
    nztmp = nz
    if ( nzsoil > nztmp ) nztmp = nzsoil
    ALLOCATE(buf_i16(nxsm,nysm,nztmp,0:nstyps))
  ENDIF


!-----------------------------------------------------------------------
!
!  Split the initial data files
!
!-----------------------------------------------------------------------

  IF (initopt == 3) THEN
    WRITE (6, *) 'Splitting initial history dump files...'
    IF (inifmt == 1) THEN
      CALL splitdump (inifile,nxsm,nysm,nz,nzsoil)
      CALL splitdump (inigbf,nxsm,nysm,nz,nzsoil)
    ELSE IF (inifmt == 3) THEN
      CALL split_hdf(inifile,nxsm,nysm,nz,nzsoil,nstyps,buf_r, &
                     buf_rsoil4d,buf_i,buf_i16,sstat)
      CALL split_hdf(inigbf,nxsm,nysm,nz,nzsoil,nstyps,buf_r, &
                     buf_rsoil4d,buf_i,buf_i16,sstat)
    ELSE
      WRITE (6, *) 'File not in binary format. Not split'
    END IF

!------------- add by wyh for split restart files ---------------------
  ELSE IF (initopt == 2) THEN
    IF (lbcopt == 2) THEN
      exbcvarsz = 22
    ELSE
      exbcvarsz = 1
    END IF

    WRITE (6, *) 'Splitting restart file...'
    CALL splitrestart (rstinf,nxsm,nysm,nz,nzsoil, nstyps, exbcvarsz)
!----------------------------------------------------------------------

  END IF
!
!-----------------------------------------------------------------------
!
!  Split the terrain data file
!
!-----------------------------------------------------------------------

  IF (ternopt == 2) THEN
    WRITE (6, *) 'Splitting terrain file...'
    IF (ternfmt == 3) THEN
      CALL split_hdf(terndta,nxsm,nysm,1,1,1,buf_r,buf_rsoil4d, &
                     buf_i,buf_i16,sstat)
    ELSE
      CALL splitterrain(terndta,nxsm,nysm)
    ENDIF
  END IF

!-----------------------------------------------------------------------
!
!  Split the surface and soil data files
!
!-----------------------------------------------------------------------

  IF (sfcdat == 2 .or. sfcdat == 3 ) THEN
    WRITE (6, *) 'Splitting surface data file...'
    IF (sfcfmt == 3) THEN
      CALL split_hdf(sfcdtfl,nxsm,nysm,1,1,1,buf_r,buf_rsoil4d,  &
                     buf_i,buf_i16,sstat)
    ELSE
      CALL split_soil(sfcdtfl,nxsm,nysm,nstyps)
    ENDIF
  END IF

  IF ((soilinit == 2 .or. soilinit == 3 ) .AND. (initopt == 3)) THEN
    WRITE (6, *) 'Splitting soil data file...'
    IF (soilfmt == 3) THEN
      CALL split_hdf(soilinfl,nxsm,nysm,1,nzsoil,nstyps,buf_r, &
                     buf_rsoil4d,buf_i,buf_i16,sstat)
    ELSE
      CALL split_soilini(soilinfl,nxsm,nysm,nzsoil,nstyps)
    ENDIF
  END IF

!-----------------------------------------------------------------------
!
!  Split the external boundary data files
!
!-----------------------------------------------------------------------

  IF (lbcopt == 2) THEN
    WRITE (6, *) 'Splitting EXBC files...'
    numfiles = nint((tstop - tstart - dtbig)/tintvebd) + 10
                           ! Go past the end in case an exbc
                           ! file past the end time is needed.
    CALL ctim2abss ( year,month,day,hour,minute,second, abstinit)
    abststop  = abstinit + nint(tstop)

!
! We will need to turn off "mpi_opt" when getbcfn() is called, as that routine
! will want to try to read the split file which doesn't exist.
!
    mp_opt_save = mp_opt
    DO i = 1, numfiles+1
      abstfcst = abstinit + (i-1) * tintvebd + nint(tstart)
      mp_opt = 0
      CALL getbcfn (abstfcst, exbcname, tinitebd, tintvebd,             &
                    filename, lfname, istat)
      mp_opt = mp_opt_save                ! Restore the value

      !wdt update: filename -> filename(1:lfname) in loop below

      IF (istat == 0) THEN
        IF (exbcfmt == 3) THEN
          CALL split_hdf(filename(1:lfname),nxsm,nysm,nz,1,1,           &
                         buf_r,buf_rsoil4d,buf_i,buf_i16,sstat)
        ELSE
          CALL split_exbc(filename(1:lfname)//' ',nxsm,nysm,nz)
        ENDIF
      END IF
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Split the "incr" files.
!
!-----------------------------------------------------------------------

  IF (nudgopt == 1 ) THEN
      WRITE (6,*) 'splitting incr file'

      CALL split_hdf(incrfnam,nxsm,nysm,nz,nzsoil,nstyps,buf_r, &
                     buf_rsoil4d,buf_i,buf_i16,sstat)
  ENDIF

  IF ((inifmt == 3) .or. (ternfmt == 3) .or. (exbcfmt == 3) .or.  &
     (soilfmt == 3) .or. (sfcfmt == 3)) THEN  ! HDF format used
    DEALLOCATE(buf_r)
    DEALLOCATE(buf_rsoil4d)
    DEALLOCATE(buf_i)
    DEALLOCATE(buf_i16)
  ENDIF

  call arpsstop("Normal Finish", 0)

END PROGRAM splitfiles
