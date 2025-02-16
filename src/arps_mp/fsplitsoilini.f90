!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SPLITSOILINI               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE split_soilini (fileheader,nx,ny,nzsoil,nstyps)

!------------------------------------------------------------------------
! HISTORY:
!
!     Yunheng Wang (09/06/2001)
!     Update with the new version ext2arps changed by Gene Bassett.
!     In general, a variable soiltyp was added in soil initial
!     file (xxxxxx.soilvar.000000) by Gene.
!
!     Deallocate all the allocated variables
!
!     Yunheng Wang (07/15/2002)
!     Updated for new ARPS version (IHOP_3)
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER :: nx,ny,nzsoil,nstyps

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
  INTEGER :: nzsoilin

!
! fmtver**: to label each data a version.
! intver**: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='* 004.10 GrADS Soilvar Data',intver410=410)
  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin

  INTEGER :: mprojin,tsoilin,qsoilin,wcanpin,stypin,zpsoilin
  INTEGER :: snowcin,snowdin,nstypin
  INTEGER :: idummy
  REAL :: dxin,dyin,ctrlonin,ctrlatin,trlat1in,trlat2in,trlonin,sclin
  REAL :: rdummy

  REAL, ALLOCATABLE :: a2dlg(:,:), a2dsm(:,:)
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
	write(6,*) 'splitsoilini:  not MP ready'
	call arpsstop('splitsoilini:   not MP ready', 1)
	return
  endif

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

!----------------------------------------------------------------------
!
!  Read/write the data format string which was added since version IHOP_3
!
!---------------------------------------------------------------------

    READ (10,ERR=997) fmtverin
                                !  To distinguish versions prior to 500.

    WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

    IF (fmtverin == fmtver500) THEN
      intver=intver500
    ELSE
      WRITE(6,'(/1x,a/)')   &
          'WARNING: Incoming data format are older than version 5.00!!! '
    END IF

    intver=intver410      ! there is no fmtverin prior to version 500
    fmtver=fmtver410

    GOTO 996

997 WRITE(6,'(/1x,a/,a/)')   &
          'ERROR: Compatibility with previous verion soil data is still', &
          '       not ready for SPLITFILES.'

    CLOSE (10)
    STOP

996 CONTINUE

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) fmtverin
    END DO
!
!-----------------------------------------------------------------------
!
!  Read/write the dimensions of data in the file and check against
!  the dimensions passed to this subroutine.
!
!-----------------------------------------------------------------------
!
    READ (10) nxin,nyin,nzsoilin

    IF ((nxin /= nxlg).OR.(nyin /= nylg).OR.(nzsoilin /= nzsoil)) THEN
      WRITE (*,*) "ERROR:  mismatch in sizes."
      WRITE (*,*) "nxin,nyin,nzsoilin",nxin,nyin,nzsoilin
      WRITE (*,*) "nxlg,nylg,nzsoil",nxlg,nylg,nzsoil
      call arpsstop("splitsoilini:  mismatch",1)
    END IF

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) nx,ny,nzsoil
    END DO

!
!-----------------------------------------------------------------------
!
!  Read/write header info
!
!-----------------------------------------------------------------------

    READ (10)                mprojin,tsoilin,qsoilin,                   &
                             wcanpin,snowcin,snowdin,stypin,zpsoilin,   &
                             idummy,idummy,idummy,idummy,idummy,        &
                             idummy,idummy,idummy,idummy,nstypin
    DO ii=1+(jj-1)*maxunit,iiend

      WRITE (ounit(ii))      mprojin,tsoilin,qsoilin,                    &
                             wcanpin,snowcin,snowdin,stypin,zpsoilin,   &
                             idummy,idummy,idummy,idummy,idummy,        &
                             idummy,idummy,idummy,idummy,nstypin
    END DO

    READ (10)                dxin,dyin, ctrlonin,ctrlatin,trlat1in,     &
                             trlat2in,trlonin,sclin,rdummy,rdummy,      &
                             rdummy,rdummy,rdummy,rdummy,rdummy,        &
                             rdummy,rdummy,rdummy,rdummy,rdummy
    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii))      dxin,dyin, ctrlonin,ctrlatin,trlat1in,     &
                             trlat2in,trlonin,sclin,rdummy,rdummy,      &
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
    nfields = 0

    IF (zpsoilin >0) nfields = nfields + nzsoilin
    IF (tsoilin > 0) nfields = nfields + nzsoilin*nstyps + nzsoil
    IF (qsoilin > 0) nfields = nfields + nzsoilin*nstyps + nzsoil
    IF (wcanpin > 0) nfields = nfields + nstyps + 1
    IF (snowdin > 0) nfields = nfields + 1
    IF (stypin  > 0) nfields = nfields + nstyps


    DO fcnt = 1,nfields

        READ (10) ((a2dlg(i,j),i=1,nxlg),j=1,nylg)

        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          DO j = 1,ny
            DO i = 1,nx
              a2dsm(i,j) = a2dlg(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3))
            END DO
          END DO

          WRITE (ounit(ii)) ((a2dsm(i,j),i=1,nx),j=1,ny)
        END DO

    END DO


  END DO     ! jj

  DEALLOCATE(a2dlg)
  DEALLOCATE(a2dsm)
  DEALLOCATE(ounit)
  DEALLOCATE(ffi)
  DEALLOCATE(ffj)

  RETURN

END SUBROUTINE split_soilini
