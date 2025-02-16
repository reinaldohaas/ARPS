
SUBROUTINE splitterrain (fileheader,nx,ny)

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER :: nx,ny

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

  REAL :: rdummy
  INTEGER :: idummy
  REAL :: dxin,dyin

  REAL, ALLOCATABLE :: a2dlg(:,:), a2dsm(:,:)
  INTEGER, ALLOCATABLE :: ounit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER :: ierr
  INTEGER :: ii,jj,iiend
  INTEGER :: unit0, maxunit
  PARAMETER (unit0=110,maxunit=60)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3

  ALLOCATE(a2dlg(nxlg,nylg))
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

  if ( mp_opt > 0 ) then
	write(6,*) 'splitterrain:  not MP ready'
	call arpsstop('splitterrain:   not MP ready')
	return
  endif

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
      call arpsstop("splitterrain:  mismatch", 1)
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
    READ (10)    idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) idummy,idummy,idummy,idummy,idummy,             &
                             idummy,idummy,idummy,idummy,idummy,        &
                             idummy,idummy,idummy,idummy,idummy,        &
                             idummy,idummy,idummy,idummy,idummy
    END DO

    READ (10)    dxin  ,dyin  ,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) dxin  ,dyin  ,rdummy,rdummy,rdummy,             &
                             rdummy,rdummy,rdummy,rdummy,rdummy,        &
                             rdummy,rdummy,rdummy,rdummy,rdummy,        &
                             rdummy,rdummy,rdummy,rdummy,rdummy
    END DO

!
!-----------------------------------------------------------------------
!
!   Read in the global data, and write out appropriate sections into
!   each processors file.
!
!-----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!  Read/write hterain.
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

    CLOSE (10)
    DO ii=1+(jj-1)*maxunit,iiend
      CLOSE (ounit(ii))
    END DO

  END DO     ! jj

  RETURN

END SUBROUTINE splitterrain
