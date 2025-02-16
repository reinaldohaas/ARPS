!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE splitdump                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE splitdump(fileheader,nx,ny,nz,nzsoil)

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  CHARACTER (LEN=*) :: fileheader

  INTEGER ::  nx,ny,nz,nzsoil

  INTEGER :: nxlg, nylg, nzlg, nzsoillg
!
!-----------------------------------------------------------------------
!
!  Variables to read in data from the data dumps
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=40) :: fmtver
  CHARACTER (LEN=80) :: runname, cmnt
  CHARACTER (LEN=10) :: tmunit
  CHARACTER (LEN=12) :: label
  INTEGER :: nocmnt
  REAL :: curtim
  INTEGER :: i01, i02, i03, i04, i05, i06, i07, i08, i09, i10
  INTEGER :: i11, i12, i13, i14, i15, i16, i17, i18, i19, i20
  INTEGER :: i21, i22, i23, i24, i25, i26, i27, i28, i29, i30
  REAL :: r01, r02, r03, r04, r05, r06, r07, r08, r09, r10
  REAL :: r11, r12, r13, r14, r15, r16, r17, r18, r19, r20
  REAL :: r21, r22, r23, r24, r25, r26, r27, r28, r29, r30

  INTEGER :: grdbas, totin
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lenstr
  CHARACTER (LEN=256) :: filename
  INTEGER :: fi, fj, i, j, k
  INTEGER :: nxin, nyin, nzin,nzsoilin

  REAL, ALLOCATABLE :: xlg(:), ylg(:), z(:)
  REAL, ALLOCATABLE :: xsm(:), ysm(:)
  REAL, ALLOCATABLE :: a3dlg(:,:,:), a3dsm(:,:,:)
  REAL, ALLOCATABLE :: a2dlg(:,:), a2dsm(:,:)
  REAL, ALLOCATABLE :: ai2dlg(:,:), ai2dsm(:,:)
  REAL, ALLOCATABLE :: a3dsoillg(:,:,:), a3dsoilsm(:,:,:)
!  REAL, ALLOCATABLE :: a4dsoillg(:,:,:,:), a4dsoilsm(:,:,:,:)

  INTEGER, ALLOCATABLE :: ounit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER :: ierr
  INTEGER :: ii,jj,iiend
  INTEGER :: unit0, maxunit
  PARAMETER (unit0=110,maxunit=60)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  if ( mp_opt > 0 ) then
	write(6,*) 'splitdump:  not MP ready'
	call arpsstop('splitdump:   not MP ready', 1)
	return
  endif

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
  nzlg = nz
  nzsoillg = nzsoil

  ALLOCATE(xlg(nxlg))
  ALLOCATE(ylg(nylg))
  ALLOCATE(z(nzlg))
  ALLOCATE(xsm(nx))
  ALLOCATE(ysm(ny))
  ALLOCATE(a3dlg(nxlg,nylg,nzlg))
  a3dlg=0.0
  ALLOCATE(a3dsm(nx,ny,nz))
  a3dsm=0.0
  ALLOCATE(a3dsoillg(nxlg,nylg,nzsoillg))
  a3dsoillg=0.0
  ALLOCATE(a3dsoilsm(nx,ny,nzsoil))
  a3dsoilsm=0.0
  ALLOCATE(a2dlg(nxlg,nylg))
  a2dlg=0.0
  ALLOCATE(a2dsm(nx,ny))
  a2dsm=0.0
  ALLOCATE(ai2dlg(nxlg,nylg))
  ai2dlg=0
  ALLOCATE(ai2dsm(nx,ny))
  ai2dsm=0

  ALLOCATE(ounit(nproc_x*nproc_y))
  ALLOCATE(ffi(nproc_x*nproc_y))
  ALLOCATE(ffj(nproc_x*nproc_y))

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

  DO fj=1,nproc_y
    DO fi=1,nproc_x

      ii = fi+nproc_x*(fj-1)
      ffi(ii) = fi
      ffj(ii) = fj
      ounit(ii) = unit0 + ii

    END DO
  END DO

  DO jj = 1,1+(nproc_x*nproc_y-1)/maxunit

    iiend = MIN(jj*maxunit,nproc_x*nproc_y)

    DO ii=1+(jj-1)*maxunit,iiend

      CALL gtsplitfn(fileheader,1,1,ffi(ii),ffj(ii),1,1,                &
                     0,0,0,2,filename,ierr)

!
!-----------------------------------------------------------------------
!
!  For compatibility with the Cray data formats. The processors
!  read their data in COS format.
!
!-----------------------------------------------------------------------
!
      CALL asnfile(filename, '-F f77 -N ieee', ierr)

      OPEN (UNIT=ounit(ii), FILE=filename, FORM='unformatted')

    END DO

    CALL asnfile(fileheader(1:lenstr), '-F f77 -N ieee', ierr)
    OPEN (UNIT=10, FILE=fileheader(1:lenstr), FORM='unformatted')
!
!-----------------------------------------------------------------------
!
!  Read/write header info.
!
!-----------------------------------------------------------------------
!
    READ (10,ERR=320,END=320) fmtver
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340) fmtver
    END DO

    READ (10,ERR=320,END=320) runname
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340) runname
    END DO

    READ (10,ERR=320,END=320) nocmnt
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340) nocmnt
    END DO

    IF ( nocmnt > 0 ) THEN
      DO i = 1,nocmnt
        READ (10,ERR=320,END=320) cmnt
        DO ii=1+(jj-1)*maxunit,iiend
          WRITE (ounit(ii),ERR=340) cmnt
        END DO
      END DO
    END IF

    READ (10,ERR=320,END=320) curtim,tmunit
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340) curtim,tmunit
    END DO

!
!-----------------------------------------------------------------------
!
!  Read/write the dimensions of data in binary file and check against
!  the dimensions passed to BINREAD.
!
!-----------------------------------------------------------------------
!
    READ (10,ERR=320,END=320) nxin,nyin,nzin,nzsoilin

    IF ((nxin /= nxlg).OR.(nyin /= nylg).OR.(nzin /= nzlg).OR. &
        (nzsoilin /= nzsoillg)) THEN
      WRITE (*,*) "ERROR:  mismatch in sizes."
      WRITE (*,*) "nxin,nyin,nzin,nzsoilin: ",nxin,nyin,nzin,nzsoilin
      WRITE (*,*) "nxlg,nylg,nzlg,nzsoillg: ",nxlg,nylg,nzlg,nzsoillg
      call arpsstop("splitdump:  mismatch", 1)
    END IF

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340) nx,ny,nz,nzsoil
    END DO

!
!-----------------------------------------------------------------------
!
!  Read/write in flags for different data groups
!
!-----------------------------------------------------------------------
!
    READ (10,ERR=320,END=320)                                           &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, totin,                                &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20

    grdbas = 0
    IF (i01 == 1 .AND. i02 == 1) grdbas = 1

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340)                                         &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, totin,                                &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20
    END DO

    IF (grdbas /= 1) THEN
      READ (10,ERR=320,END=320)                                         &
                i01, i02, i03, i04, i05, i06,                           &
                i07, i08, i09, i10, i11, i12,                           &
                i13, i14, i15, i16, i17, i18,                           &
                i19, i20, i21, i22, i23, i24,                           &
                i25, i26, i27, i28, i29, i30

      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii),ERR=340)                                       &
                i01, i02, i03, i04, i05, i06,                           &
                i07, i08, i09, i10, i11, i12,                           &
                i13, i14, i15, i16, i17, i18,                           &
                i19, i20, i21, i22, i23, i24,                           &
                i25, i26, i27, i28, i29, i30
      END DO
    END IF

    IF ( grdbas == 1) THEN
      READ (10,ERR=320,END=320)                                         &
                r01, r02, r03, r04, r05,                                &
                r06, r07, r08, r09, r10,                                &
                r11, r12, r13, r14, r15,                                &
                r16, r17, r18, r19, r20

      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii),ERR=340)                                       &
                r01, r02, r03, r04, r05,                                &
                r06, r07, r08, r09, r10,                                &
                r11, r12, r13, r14, r15,                                &
                r16, r17, r18, r19, r20
      END DO
    ELSE
      READ (10,ERR=320,END=320)                                         &
                r01, r02, r03, r04, r05,                                &
                r06, r07, r08, r09, r10,                                &
                r11, r12, r13, r14, r15,                                &
                r16, r17, r18, r19, r20,                                &
                r21, r22, r23, r24, r25,                                &
                r26, r27, r28, r29, r30

      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii),ERR=340)                                       &
                r01, r02, r03, r04, r05,                                &
                r06, r07, r08, r09, r10,                                &
                r11, r12, r13, r14, r15,                                &
                r16, r17, r18, r19, r20,                                &
                r21, r22, r23, r24, r25,                                &
                r26, r27, r28, r29, r30
      END DO
    END IF

    IF (totin == 1) THEN

      READ (10,ERR=320,END=320)                                         &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, i10,                                  &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20

      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii),ERR=340)                                       &
              i01, i02, i03, i04, i05,                                  &
              i06, i07, i08, i09, i10,                                  &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20
      END DO

      READ (10,ERR=320,END=320)                                         &
              r01, r02, r03, r04, r05,                                  &
              r06, r07, r08, r09, r10,                                  &
              r11, r12, r13, r14, r15,                                  &
              r16, r17, r18, r19, r20

      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii),ERR=340)                                       &
              r01, r02, r03, r04, r05,                                  &
              r06, r07, r08, r09, r10,                                  &
              r11, r12, r13, r14, r15,                                  &
              r16, r17, r18, r19, r20
      END DO

    END IF

!
!----------------------------------------------------------------------
!
!  For every 1-, 2-, or 3-d set of data in the input file, read in
!  the array and then write out each processor's section of the data.
!
!----------------------------------------------------------------------
!

    400   CONTINUE
    READ (10,ERR=320,END=310) label
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii),ERR=340) label
    END DO

    IF (label(10:10) == "1") THEN
      IF (label(12:12) == "1") THEN           ! 1-d x
!
!----------------------------------------------------------------------
!
!  x.
!
!----------------------------------------------------------------------
!
        READ (10,ERR=320,END=310) xlg

        DO ii=1+(jj-1)*maxunit,iiend
          DO i = 1,nx
            xsm(i) = xlg(i+(ffi(ii)-1)*(nx-3))
          END DO
          WRITE (ounit(ii),ERR=340) xsm
        END DO

      ELSE IF (label(12:12) == "2") THEN      ! 1-d y
!
!----------------------------------------------------------------------
!
!  y.
!
!----------------------------------------------------------------------
!
        READ (10,ERR=320,END=310) ylg

        DO ii=1+(jj-1)*maxunit,iiend
          DO j = 1,ny
            ysm(j) = ylg(j+(ffj(ii)-1)*(ny-3))
          END DO
          WRITE (ounit(ii),ERR=340) ysm
        END DO

      ELSE IF (label(12:12) == "3") THEN      ! 1-d z
!
!----------------------------------------------------------------------
!
!  z.
!
!----------------------------------------------------------------------
!
        READ (10,ERR=320,END=310) z

        DO ii=1+(jj-1)*maxunit,iiend
          WRITE (ounit(ii),ERR=340) z
        END DO

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
        READ (10,ERR=320,END=310) a2dlg

        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          DO j = 1,ny
            DO i = 1,nx
              a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
            END DO
          END DO
          WRITE (ounit(ii),ERR=340) a2dsm
        END DO

      ELSE IF (label(9:9) == "i") THEN        ! 2-d integer
!
!----------------------------------------------------------------------
!
!  2-d integer array.
!
!----------------------------------------------------------------------
!
        READ (10,ERR=320,END=310) ai2dlg

        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          DO j = 1,ny
            DO i = 1,nx
              ai2dsm(i,j) = ai2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
            END DO
          END DO
          WRITE (ounit(ii),ERR=340) ai2dsm
        END DO

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
      READ(10,ERR=320,END=310) a3dsoillg
      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nzsoil
          DO j = 1,ny
            DO i = 1,nx
              a3dsoilsm(i,j,k) =                                      &
                  a3dsoillg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO

        WRITE (ounit(ii),ERR=340) a3dsoilsm
      END DO

    ELSE IF (label(10:10) == "3") THEN        ! 3-d
!
!----------------------------------------------------------------------
!
!  3-d real array.
!
!----------------------------------------------------------------------
!
      READ (10,ERR=320,END=310) a3dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              a3dsm(i,j,k) =                                            &
                  a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii),ERR=340) a3dsm
      END DO

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
    CLOSE (10)
    DO ii=1+(jj-1)*maxunit,iiend
      CLOSE (ounit(ii))
    END DO

  END DO   ! jj

  RETURN

!
!-----------------------------------------------------------------------
!
!  Error during read.
!
!----------------------------------------------------------------------
!

  320   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in SPLITDUMP'
  call arpsstop("splitdump:  320 continue", 1)

!
!----------------------------------------------------------------------
!
!  Error with the label.
!
!----------------------------------------------------------------------
!
  330   CONTINUE

  WRITE(6,'(a,a)') ' Error with label in SPLITDUMP:',label
  call arpsstop("splitdump:  330 continue", 1)

!
!----------------------------------------------------------------------
!
!  Error with write.
!
!----------------------------------------------------------------------
!
  340   CONTINUE

  WRITE(6,'(a,a)') ' Error with write in SPLITDUMP.'
  call arpsstop("splitdump:  340 continue", 1)

END SUBROUTINE splitdump

