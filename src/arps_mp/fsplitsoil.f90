
SUBROUTINE split_soil (fileheader,nx,ny,nstyps)

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER nx,ny,nstyps

  CHARACTER (LEN=*) :: fileheader

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxlg, nylg

  INTEGER :: lenstr
  CHARACTER (LEN=256) :: filename
  INTEGER :: fi, fj, i, j
  INTEGER :: nxin, nyin

  INTEGER :: mprojin,stypin,vtypin,laiin,roufin,vegin,ndviin,nstypin
  INTEGER :: idummy
  REAL :: dxin,dyin, ctrlonin,ctrlatin,trlat1in
  REAL :: trlat2in,trlonin,sclin
  REAL :: rdummy

  REAL,    ALLOCATABLE :: a2dlg(:,:), a2dsm(:,:)
  INTEGER, ALLOCATABLE :: i2dlg(:,:), i2dsm(:,:)
  INTEGER, ALLOCATABLE :: ounit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER, PARAMETER :: unit0=110, maxunit=60

  INTEGER :: ierr
  INTEGER :: ii,jj,iiend
  INTEGER :: is

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  if ( mp_opt > 0 ) then
	write(6,*) 'splitsoil:  not MP ready'
	call arpsstop('splitsoil:   not MP ready', 1)
	return
  endif

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3

  ALLOCATE(a2dlg(nxlg,nylg))
  ALLOCATE(i2dlg(nxlg,nylg))
  ALLOCATE(i2dsm(nx,ny))
  ALLOCATE(a2dsm(nx,ny))
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
!  Open the files.
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

!
!-----------------------------------------------------------------------
!
!  Since T3D processors only support COS and IEEE double precision
!  format, we have to translate the files into COS format.
!
!-----------------------------------------------------------------------
!

      CALL gtsplitfn(fileheader,1,1,ffi(ii),ffj(ii),1,1,0,0,0,2,filename,ierr)

      CALL asnfile(filename, '-F f77 -N ieee', ierr)

      OPEN (UNIT=ounit(ii), FILE=filename, FORM='unformatted')

    END DO

    CALL asnfile(fileheader(1:lenstr), '-F f77 -N ieee', ierr)
    OPEN (UNIT=10, FILE=fileheader(1:lenstr), FORM='unformatted')

!
!-----------------------------------------------------------------------
!
!  Read/write the dimensions of data in the file and check against
!  the dimensions passed to this subroutine.
!
!-----------------------------------------------------------------------
!
    READ (10) nxin,nyin
    IF ((nxin /= nxlg).OR.(nyin /= nylg)) THEN
      WRITE (*,*) "ERROR:  mismatch in sizes."
      WRITE (*,*) "nxin,nyin",nxin,nyin
      WRITE (*,*) "nxlg,nylg",nxlg,nylg
      call arpsstop("splitsoil:  mismatch", 1)
    END IF

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) nx,ny
    END DO

!
!-----------------------------------------------------------------------
!
!  Read/write header info.
!
!-----------------------------------------------------------------------
!

    READ (10) mprojin,stypin,vtypin,laiin,roufin,                       &
              vegin, ndviin,nstypin,idummy,idummy,                      &
              idummy, idummy,idummy,idummy,idummy,                      &
              idummy, idummy,idummy,idummy,idummy

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) mprojin,stypin,vtypin,laiin,roufin,             &
                             vegin, ndviin,nstypin,idummy,idummy,       &
                             idummy, idummy,idummy,idummy,idummy,       &
                             idummy, idummy,idummy,idummy,idummy
    END DO

    READ (10) dxin,dyin, ctrlonin,ctrlatin,trlat1in,                    &
              trlat2in,trlonin,sclin,rdummy,rdummy,                     &
              rdummy,rdummy,rdummy,rdummy,rdummy,                       &
              rdummy,rdummy,rdummy,rdummy,rdummy
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) dxin,dyin, ctrlonin,ctrlatin,trlat1in,          &
                             trlat2in,trlonin,sclin,rdummy,rdummy,      &
                             rdummy,rdummy,rdummy,rdummy,rdummy,        &
                             rdummy,rdummy,rdummy,rdummy,rdummy
    END DO
!
!-----------------------------------------------------------------------
!
!  Read in the global data, and write out appropriate sections into
!  each processors file.
!
!-----------------------------------------------------------------------
!

    IF ( stypin /= 0 ) THEN
      IF ( nstypin == 1 ) THEN
!
!----------------------------------------------------------------------
!
!  Read/write soiltyp
!
!----------------------------------------------------------------------
!
        READ (10) i2dlg

        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          DO j = 1,ny
            DO i = 1,nx
              i2dsm(i,j) = i2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
            END DO
          END DO
          WRITE (ounit(ii)) i2dsm
        END DO

      ELSE IF ( nstypin > 1 ) THEN

        DO is=1,nstypin

          READ (10) i2dlg

          DO ii=1+(jj-1)*maxunit,iiend
            fi = ffi(ii)
            fj = ffj(ii)
            DO j = 1,ny
              DO i = 1,nx
                i2dsm(i,j) = i2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
              END DO
            END DO
            WRITE (ounit(ii)) i2dsm
          END DO

          READ (10) a2dlg

          DO ii=1+(jj-1)*maxunit,iiend
            fi = ffi(ii)
            fj = ffj(ii)
            DO j = 1,ny
              DO i = 1,nx
                a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
              END DO
            END DO
            WRITE (ounit(ii)) a2dsm
          END DO

        END DO

      END IF
    END IF

    IF ( vtypin /= 0 ) THEN
!
!----------------------------------------------------------------------
!
!  Read/write vegtyp
!
!----------------------------------------------------------------------
!
      READ (10) i2dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1,ny
          DO i = 1,nx
            i2dsm(i,j) = i2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) i2dsm
      END DO

    END IF

    IF ( laiin /= 0 ) THEN
!
!----------------------------------------------------------------------
!
!  Read/write lai
!
!----------------------------------------------------------------------
!
      READ (10) a2dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1,ny
          DO i = 1,nx
            a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dsm
      END DO

    END IF

    IF ( roufin /= 0 ) THEN
!
!----------------------------------------------------------------------
!
!  Read/write roufns
!
!----------------------------------------------------------------------
!
      READ (10) a2dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1,ny
          DO i = 1,nx
            a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dsm
      END DO

    END IF

    IF ( vegin /= 0 ) THEN
!
!----------------------------------------------------------------------
!
!  Read/write veg
!
!----------------------------------------------------------------------
!
      READ (10) a2dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1,ny
          DO i = 1,nx
            a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dsm
      END DO

    END IF

    IF ( ndviin /= 0 ) THEN
!
!----------------------------------------------------------------------
!
!  Read/write ndvi
!
!----------------------------------------------------------------------
!
      READ (10) a2dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1,ny
          DO i = 1,nx
            a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dsm
      END DO

    END IF

    CLOSE (10)
    DO ii=1+(jj-1)*maxunit,iiend
      CLOSE (ounit(ii))
    END DO

  END DO    ! jj

  RETURN

END SUBROUTINE split_soil
