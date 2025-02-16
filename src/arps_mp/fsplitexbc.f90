
SUBROUTINE split_exbc(fileheader,nx,ny,nz)

  IMPLICIT NONE

  CHARACTER (LEN=80) :: fileheader

  INCLUDE 'mp.inc'

  INTEGER nx,ny,nz

  INTEGER :: nxlg, nylg, nzlg
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
  INTEGER :: nxin, nyin, nzin

  REAL :: dxin,dyin,dzin,ctrlatin,ctrlonin
  INTEGER :: ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd
  INTEGER :: qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd
  INTEGER :: ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd
  INTEGER :: zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd,idummy
  INTEGER :: nscalarin, old_v

  CHARACTER (LEN=15) :: ctime

  REAL, ALLOCATABLE :: a3dlg(:,:,:), a3dsm(:,:,:)
  INTEGER, ALLOCATABLE :: ounit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER :: ierr
  INTEGER :: nfields, fcnt
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

  if ( mp_opt > 0 ) then
	write(6,*) 'splitexbc:  not MP ready'
	call arpsstop('splitexbc:   not MP ready', 1)
	return
  endif

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
  nzlg = nz

  ALLOCATE(a3dlg(nxlg,nylg,nzlg))
  ALLOCATE(a3dsm(nx,ny,nz))
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
!  Split the original data file into indivdual files for the
!  processors to read.
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
      CALL gtsplitfn(fileheader,1,1,ffi(ii),ffj(ii),1,1,                &
                     0,0,0,2,filename,ierr)

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
    READ (10)                                                           &
          nxin,nyin,nzin,dxin,dyin,dzin,ctrlatin,ctrlonin
    IF ((nxin /= nxlg).OR.(nyin /= nylg).OR.(nzin /= nzlg)) THEN
      WRITE (*,*) "ERROR:  mismatch in sizes."
      WRITE (*,*) "nxin,nyin,nzin: ",nxin,nyin,nzin
      WRITE (*,*) "nxlg,nylg,nzlg: ",nxlg,nylg,nzlg
      call arpsstop("splitexbc:  mismatch", 1)
    END IF

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii))                                                 &
               nx,ny,nz,dxin,dyin,dzin,ctrlatin,ctrlonin
    END DO

!
!-----------------------------------------------------------------------
!
!  Read/write header info.
!
!-----------------------------------------------------------------------
!
    READ (10) ctime
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) ctime
    END DO

    old_v = 0 ! In case that the EXBC files are of an earlier
              ! version that does not contain water and ice variables,
              ! set old_v to 1. Otherwise, set it to 0.

    IF( old_v == 1 ) THEN

      READ (10) ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd
      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii)) ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd
      END DO

      qcbcrd = 0
      qrbcrd = 0
      qibcrd = 0
      qsbcrd = 0
      qhbcrd = 0

    ELSE

      READ (10)   ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                      &
                  qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,                   &
                  qhbcrd,qgbcrd,ncbcrd,nrbcrd,nibcrd,                   &
                  nsbcrd,ngbcrd,nhbcrd,zrbcrd,zibcrd,                   &
                  zsbcrd,zgbcrd,zhbcrd,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,nscalarin

      DO ii=1+(jj-1)*maxunit,iiend
        WRITE (ounit(ii))                                               &
                  ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                      &
                  qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,                   &
                  qhbcrd,qgbcrd,ncbcrd,nrbcrd,nibcrd,                   &
                  nsbcrd,ngbcrd,nhbcrd,zrbcrd,zibcrd,                   &
                  zsbcrd,zgbcrd,zhbcrd,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,nscalarin
      END DO

    END IF

    nfields = 6
    IF (qcbcrd > 0) nfields = nfields + 1
    IF (qrbcrd > 0) nfields = nfields + 1
    IF (qibcrd > 0) nfields = nfields + 1
    IF (qsbcrd > 0) nfields = nfields + 1
    IF (qhbcrd > 0) nfields = nfields + 1
    IF (qgbcrd > 0) nfields = nfields + 1
    IF (ncbcrd > 0) nfields = nfields + 1
    IF (nrbcrd > 0) nfields = nfields + 1
    IF (nibcrd > 0) nfields = nfields + 1
    IF (nsbcrd > 0) nfields = nfields + 1
    IF (nhbcrd > 0) nfields = nfields + 1
    IF (ngbcrd > 0) nfields = nfields + 1
    IF (zrbcrd > 0) nfields = nfields + 1
    IF (zibcrd > 0) nfields = nfields + 1
    IF (zsbcrd > 0) nfields = nfields + 1
    IF (zhbcrd > 0) nfields = nfields + 1
    IF (zgbcrd > 0) nfields = nfields + 1
!
!-----------------------------------------------------------------------
!
!  Read in the global data, and write out appropriate sections into
!  each processors file.
!
!-----------------------------------------------------------------------
!
    DO fcnt = 1,nfields

      READ (10) a3dlg

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
        WRITE (ounit(ii)) a3dsm
      END DO


    END DO

    CLOSE (10)
    DO ii=1+(jj-1)*maxunit,iiend
      CLOSE (ounit(ii))
    END DO

  END DO    ! jj

  RETURN

END SUBROUTINE split_exbc

