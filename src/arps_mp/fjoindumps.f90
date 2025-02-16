
SUBROUTINE joindumps (fileheader,nx,ny,nz,nzsoil)

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER :: nx,ny,nz,nzsoil

  CHARACTER (LEN=*) :: fileheader
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxlg, nylg, nzlg, nzsoillg
  INTEGER :: lenstr
  CHARACTER (LEN=10 ) :: filetail
  CHARACTER (LEN=256) :: filename
  INTEGER :: fi, fj, i, j, k
  INTEGER :: nxin, nyin, nzin, nzsoilin

  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500,fmtver530
  INTEGER            :: intver,intver410,intver500,intver530
  PARAMETER (fmtver410='004.10 Binary Data',intver410=410)
  PARAMETER (fmtver500='005.00 Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 Binary Data',intver530=530)

  CHARACTER (LEN=80) :: runname, cmnt
  CHARACTER (LEN=10) :: tmunit
  CHARACTER (LEN=12) :: label
  INTEGER :: nocmnt
  REAL :: curtim
  INTEGER :: i01, i02, i03, i04, i05, i06, i07, i08, i09, i10
  INTEGER :: i11, i12, i13, i14, i15, i16, i17, i18, i19, i20
  INTEGER :: i21, i22, i23, i24, i25, i26, i27, i28, i29, i30
  INTEGER :: totin, grdbas
  REAL :: r01, r02, r03, r04, r05, r06, r07, r08, r09, r10
  REAL :: r11, r12, r13, r14, r15, r16, r17, r18, r19, r20
  REAL :: r21, r22, r23, r24, r25, r26, r27, r28, r29, r30

  INTEGER :: ierr
  LOGICAL :: fexist

  REAL, ALLOCATABLE :: xlg(:), ylg(:), z(:)
  REAL, ALLOCATABLE :: xsm(:), ysm(:)
  REAL, ALLOCATABLE :: a3dsoillg(:,:,:), a3dsoilsm(:,:,:)
  REAL, ALLOCATABLE :: a3dlg(:,:,:), a3dsm(:,:,:)
  REAL, ALLOCATABLE :: a2dlg(:,:), a2dsm(:,:)
  INTEGER, ALLOCATABLE :: ai2dlg(:,:), ai2dsm(:,:)
  INTEGER, ALLOCATABLE :: i0(:,:), j0(:,:)

  INTEGER, ALLOCATABLE :: iunit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER :: ii,jj,iiend
  INTEGER :: unit0, maxunit
  PARAMETER (unit0=110,maxunit=60)

  INTEGER :: joff, junit0
  PARAMETER(junit0=11)
  CHARACTER (LEN=256) :: outfile
  CHARACTER (LEN=256) :: outfile_old

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = nproc_x*(nx-3)+3
  nylg = nproc_y*(ny-3)+3
  nzlg = nz
  nzsoillg = nzsoil

  ALLOCATE(xlg(nxlg))
  ALLOCATE(ylg(nylg))
  ALLOCATE(z(nzlg))
  ALLOCATE(xsm(nx))
  ALLOCATE(ysm(ny))
  ALLOCATE(a3dsoillg(nxlg,nylg,nzsoillg))
  ALLOCATE(a3dlg(nxlg,nylg,nzlg))
  ALLOCATE(a3dsoilsm(nx,ny,nzsoil))
  ALLOCATE(a3dsm(nx,ny,nz))
  ALLOCATE(a2dlg(nxlg,nylg))
  ALLOCATE(a2dsm(nx,ny))
  ALLOCATE(ai2dlg(nxlg,nylg))
  ALLOCATE(ai2dsm(nx,ny))
  ALLOCATE(i0(nproc_x,nproc_y))
  ALLOCATE(j0(nproc_x,nproc_y))

  ALLOCATE(iunit(nproc_x*nproc_y))
  ALLOCATE(ffi(nproc_x*nproc_y))
  ALLOCATE(ffj(nproc_x*nproc_y))

  joff = 0

  lenstr = 0
  100   lenstr = lenstr + 1
  IF (fileheader(lenstr:lenstr) /= " ") GO TO 100
  lenstr = lenstr - 1

!
!-----------------------------------------------------------------------
!
!  Open the split files.
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)

  DO fj = 1,nproc_y
    DO fi = 1,nproc_x

      IF (fi == 1) THEN
        i0(fi,fj) = 1
      ELSE
        i0(fi,fj) = 2
      END IF

      IF (fj == 1) THEN
        j0(fi,fj) = 1
      ELSE
        j0(fi,fj) = 2
      END IF

      ii = fi+nproc_x*(fj-1)
      ffi(ii) = fi
      ffj(ii) = fj
      iunit(ii) = unit0 + ii

    END DO
  END DO

  a3dsoillg = 0.0
  a3dsoilsm = 0.0

  DO k = 1,nz
    DO j = 1,nylg
      DO i = 1,nxlg
        a3dlg(i,j,k) = 0.0
      END DO
    END DO
  END DO
  DO j = 1,nylg
    DO i = 1,nxlg
      a2dlg(i,j) = 0.0
      ai2dlg(i,j) = 0
    END DO
  END DO

  DO jj = 1,1+(nproc_x*nproc_y-1)/maxunit

    iiend = MIN(jj*maxunit,nproc_x*nproc_y)

    DO ii=1+(jj-1)*maxunit,iiend

!
!-----------------------------------------------------------------------
!
!   For compatibility with the Cray data formats. The processors
!   read their data in COS format.
!
!-----------------------------------------------------------------------
!
      CALL gtsplitfn(fileheader,nproc_x,nproc_y,1,1,ffi(ii),ffj(ii),    &
                     0,0,1,2,filename,ierr)

      INQUIRE (FILE=filename, EXIST=fexist)
      IF ( .NOT. fexist) THEN
        WRITE (6,*) 'Parts of ',fileheader,' were not found'
        WRITE (6,*) 'No file joining is done for this time.'
        WRITE (6,*) 'Program continues.'
        WRITE (6,*)
        RETURN
      END IF

      CALL asnfile(filename, '-F f77 -N ieee', ierr)
      OPEN (UNIT=iunit(ii),FILE=trim(filename),FORM='unformatted')

    END DO

    outfile_old(1:256) = outfile(1:256)
    IF ( iiend == nproc_x*nproc_y ) THEN
      WRITE(outfile, '(a)') fileheader(1:lenstr)
    ELSE
      WRITE(outfile, '(a,a,i3.3)')                                      &
          fileheader(1:lenstr),'_tmp',iiend
    END IF
    CALL asnfile(outfile, '-F f77 -N ieee', ierr)
    OPEN (UNIT=junit0+joff,FILE=outfile,FORM='unformatted')

    IF (joff > 0 ) OPEN (UNIT=junit0+joff-1,FILE=outfile_old,FORM='unformatted')

!
!-----------------------------------------------------------------------
!
!  Read/write header info
!
!-----------------------------------------------------------------------
!

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) fmtver
    END DO
    IF (fmtver == fmtver500) THEN
      intver = intver500
    ELSE IF (fmtver == fmtver530) THEN
      intver = intver530
    ELSE
      WRITE(6,'(1x,3a)') 'ERROR: fmtver = ',fmtver,' is not supported.'
      CALL arpsstop('Unsupported file format.',1)
    END IF

    IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff) fmtver

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) runname
    END DO
    IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff) runname

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) nocmnt
    END DO
    IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff) nocmnt

    IF ( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        DO ii=1+(jj-1)*maxunit,iiend
          READ (iunit(ii)) cmnt
        END DO
        IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff) cmnt
      END DO
    END IF

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) curtim,tmunit
    END DO
    IF (iiend == nproc_x*nproc_y)   WRITE (junit0+joff) curtim,tmunit

!
!-----------------------------------------------------------------------
!
!  Read/write dimensions of data in binary file and check against
!  the dimensions passed to BINREAD
!
!-----------------------------------------------------------------------
!

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) nxin,nyin,nzin,nzsoilin
    END DO

    IF ((nxin /= nx).OR.(nyin /= ny).OR.(nzin /= nz).OR. &
        (nzsoilin /= nzsoil)) THEN
      WRITE (*,*) "ERROR:  missmatch in sizes."
      WRITE (*,*) "nxin,nyin,nzin,nzsoilin",nxin,nyin,nzin,nzsoilin
      WRITE (*,*) "nx,ny,nz,nzsoil",nx,ny,nz,nzsoil
      STOP
    END IF

    IF (iiend == nproc_x*nproc_y) &
                   WRITE (junit0+joff) nxlg,nylg,nzlg,nzsoillg

!
!-----------------------------------------------------------------------
!
!  Read/write flags for different data groups.
!
!-----------------------------------------------------------------------
!

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii))                                                  &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, totin,                                &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20
    END DO

    IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff)                   &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, totin,                                &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20

    IF (i01 == 1 .AND. i02 == 1) THEN    ! suppose grid and base file
      grdbas = 1
    ELSE
      grdbas = 0
    END IF

    IF (grdbas == 0 .AND. intver >= intver530) THEN
      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                  &
                i01, i02, i03, i04, i05, i06,                           &
                i07, i08, i09, i10, i11, i12,                           &
                i13, i14, i15, i16, i17, i18,                           &
                i19, i20, i21, i22, i23, i24,                           &
                i25, i26, i27, i28, i29, i30
      END DO

      IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff)             &
                i01, i02, i03, i04, i05, i06,                           &
                i07, i08, i09, i10, i11, i12,                           &
                i13, i14, i15, i16, i17, i18,                           &
                i19, i20, i21, i22, i23, i24,                           &
                i25, i26, i27, i28, i29, i30
    END IF

    IF (intver >= intver530) THEN
      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                  &
                r01, r02, r03, r04, r05,                                  &
                r06, r07, r08, r09, r10,                                  &
                r11, r12, r13, r14, r15,                                  &
                r16, r17, r18, r19, r20,                                  &
                r21, r22, r23, r24, r25,                                  &
                r26, r27, r28, r29, r30
      END DO

      IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff)             &
                r01, r02, r03, r04, r05,                                  &
                r06, r07, r08, r09, r10,                                  &
                r11, r12, r13, r14, r15,                                  &
                r16, r17, r18, r19, r20,                                  &
                r21, r22, r23, r24, r25,                                  &
                r26, r27, r28, r29, r30
    ELSE
      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                  &
                r01, r02, r03, r04, r05,                                  &
                r06, r07, r08, r09, r10,                                  &
                r11, r12, r13, r14, r15,                                  &
                r16, r17, r18, r19, r20
      END DO

      IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff)             &
                r01, r02, r03, r04, r05,                                  &
                r06, r07, r08, r09, r10,                                  &
                r11, r12, r13, r14, r15,                                  &
                r16, r17, r18, r19, r20
    END IF

    IF (totin == 1) THEN

      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, i10,                                  &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20
      END DO

      IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff)           &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, i10,                                  &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20

      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                &
              r01, r02, r03, r04, r05,                                  &
              r06, r07, r08, r09, r10,                                  &
              r11, r12, r13, r14, r15,                                  &
              r16, r17, r18, r19, r20
      END DO

      IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff)           &
              r01, r02, r03, r04, r05,                                  &
              r06, r07, r08, r09, r10,                                  &
              r11, r12, r13, r14, r15,                                  &
              r16, r17, r18, r19, r20

    END IF

!
!----------------------------------------------------------------------
!
!  For every 1-, 2-, or 3-d set of data in the input file, read in
!  the arrays from each processor's file and write out the
!  combined data.
!
!----------------------------------------------------------------------
!

    400   CONTINUE
    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii),END=310) label
    END DO
    IF (iiend == nproc_x*nproc_y) WRITE (junit0+joff) label

    IF (label(10:10) == "1") THEN
      IF (label(12:12) == "1") THEN           ! 1-d x
!
!----------------------------------------------------------------------
!
!  x.
!
!----------------------------------------------------------------------
!

        IF (joff > 0 ) READ (junit0+joff-1) xlg
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) xsm
          IF (fj == 1) THEN
            DO i=1,nx
              xlg(i+(fi-1)*(nx-3)) = xsm(i)
            END DO
          END IF
        END DO
        WRITE (junit0+joff) xlg

      ELSE IF (label(12:12) == "2") THEN      ! 1-d y
!
!----------------------------------------------------------------------
!
!  y.
!
!----------------------------------------------------------------------
!

        IF (joff > 0 ) READ (junit0+joff-1) ylg
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) ysm
          IF (fi == 1) THEN
            DO j=1,ny
              ylg(j+(fj-1)*(ny-3)) = ysm(j)
            END DO
          END IF
        END DO
        WRITE (junit0+joff) ylg

      ELSE IF (label(12:12) == "3") THEN      ! 1-d z
!
!----------------------------------------------------------------------
!
!  z.
!
!----------------------------------------------------------------------
!
        IF (joff > 0 ) READ (junit0+joff-1) z
        DO ii=1+(jj-1)*maxunit,iiend
          READ (iunit(ii)) z
        END DO
        WRITE (junit0+joff) z

      ELSE
        GO TO 330
      END IF
    ELSE IF (label(10:10) == "2") THEN
      IF (label(9:9) == "r") THEN             ! 2-d real
!
!----------------------------------------------------------------------
!
!  2-d real array.
!
!----------------------------------------------------------------------
!
        IF (joff > 0 ) READ (junit0+joff-1) a2dlg
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) a2dsm
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              a2dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3)) = a2dsm(i,j)
            END DO
          END DO
        END DO
        WRITE (junit0+joff) a2dlg

      ELSE IF (label(9:9) == "i") THEN        ! 2-d integer
!
!----------------------------------------------------------------------
!
!  2-d integer array.
!
!----------------------------------------------------------------------
!
        IF (joff > 0 ) READ (junit0+joff-1) ai2dlg
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) ai2dsm
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              ai2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3)) = ai2dsm(i,j)
            END DO
          END DO
        END DO
        WRITE (junit0+joff) ai2dlg

      ELSE
        GO TO 330
      END IF

    ELSE IF (label(9:10) == "s3") THEN        ! 3-d with nzsoil
!
!----------------------------------------------------------------------
!
!  3-d real soil array.
!
!----------------------------------------------------------------------
!
      IF (joff > 0 ) READ (junit0+joff-1) a3dsoillg
      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        READ (iunit(ii)) a3dsoilsm
        DO k = 1,nzsoil
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              a3dsoillg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k) =         &
                         a3dsoilsm(i,j,k)
            END DO
          END DO
        END DO
      END DO
      WRITE (junit0+joff) a3dsoillg

    ELSE IF (label(10:10) == "3") THEN        ! 3-d
!
!----------------------------------------------------------------------
!
!  3-d real array.
!
!----------------------------------------------------------------------
!
      IF (joff > 0 ) READ (junit0+joff-1) a3dlg
      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        READ (iunit(ii)) a3dsm
        DO k = 1,nz
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k) =              &
                         a3dsm(i,j,k)
            END DO
          END DO
        END DO
      END DO
      WRITE (junit0+joff) a3dlg

    ELSE
      GO TO 330
    END IF

    GO TO 400

!
!-----------------------------------------------------------------------
!
!  Error free finish.  Close files and return.
!
!----------------------------------------------------------------------
!
    310   CONTINUE

    DO ii=1+(jj-1)*maxunit,iiend
      CLOSE (iunit(ii))
    END DO

    CLOSE (junit0+joff)
    IF (joff > 0) THEN
      CLOSE (junit0+joff-1,STATUS='delete')
    END IF

    joff = joff + 1
  END DO

  DEALLOCATE(xlg)
  DEALLOCATE(ylg)
  DEALLOCATE(z)
  DEALLOCATE(xsm)
  DEALLOCATE(ysm)
  DEALLOCATE(a3dsoillg)
  DEALLOCATE(a3dlg)
  DEALLOCATE(a3dsoilsm)
  DEALLOCATE(a3dsm)
  DEALLOCATE(a2dlg)
  DEALLOCATE(a2dsm)
  DEALLOCATE(ai2dlg)
  DEALLOCATE(ai2dsm)
  DEALLOCATE(i0)
  DEALLOCATE(j0)

  DEALLOCATE(iunit)
  DEALLOCATE(ffi)
  DEALLOCATE(ffj)

  RETURN

!
!-----------------------------------------------------------------------
!
!  Error during read.
!
!----------------------------------------------------------------------
!

!  320   CONTINUE
!  WRITE(6,'(/a/)') ' Error reading data in JOINDUMPS'
!  STOP 320

!
!----------------------------------------------------------------------
!
!  Error with the label.
!
!----------------------------------------------------------------------
!
  330   CONTINUE

  WRITE(6,'(a,a)') ' Error with label in JOINDUMPS:',label
  STOP 330

!
!----------------------------------------------------------------------
!
!  Error with write.
!
!----------------------------------------------------------------------
!
!  340   CONTINUE

!  WRITE(6,'(a,a)') ' Error with write in JOINDUMPS.'
!  STOP 340

END SUBROUTINE joindumps

