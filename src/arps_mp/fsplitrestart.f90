!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SPLITRESTART               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE splitrestart(fileheader,nx,ny,nz,nzsoil,nstyps, exbcvarsz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Split the original restart file into indivdual files for the
!  processors to read.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Yunheng Wang
!         4/16/2001.
!
! MODIFICATION:
!
!   Yunheng Wang (2003/02/27)
!   Changed the input and output format for exbcbuf. It is now outputted
!   as 22 separated varaibles instead of a big buffer. Actually, it should
!   be a bug fix because it did not consider the contents in exbcbuf before.
!
!-----------------------------------------------------------------------
!
! INPUT:
!
!  fileheader       Restart file name
!  nx               Dimension in x-direction (each processor)
!  ny               Dimension in y-direction (each processor)
!  nz               Dimension in z-direction
!  nzsoil           Number of soil levels
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=*) :: fileheader

  INTEGER :: nx,ny,nz,nzsoil

  INTEGER :: nstyps, exbcvarsz

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lenstr
  CHARACTER (LEN=256) :: filename
  INTEGER :: fi, fj, i, j, k, l, var

  INTEGER :: nxlg, nylg, nzlg,nzsoillg     ! Dimension in the whole
                                           ! computation domain
  INTEGER, ALLOCATABLE :: ounit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER :: ierr
  INTEGER :: nfields, fcnt
  INTEGER :: ii,jj,iiend

  INTEGER, PARAMETER :: unit0=110, maxunit=60

!----------------------------------------------------------------------
!
! Variable to read and write out.
!
!----------------------------------------------------------------------

  REAL    :: ctime

  INTEGER :: nxin, nyin, nzin,nzsoilin

  INTEGER :: basrstout, grdrstout, icerstout, sfcrstout, prcrstout,    &
             rcumout, exbcout, mapfout, radrstout, nstyp, kfrsout,     &
             bmjsout, idummy                        ! Control parameters

  INTEGER :: ij(30)

  INTEGER :: abstfcst0, abstfcst,                                 &
             ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                     &
             qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd
  INTEGER :: ncbcrd,nrbcrd,nibcrd,nsbcrd,nhbcrd,ngbcrd
  INTEGER :: zrbcrd,zibcrd,zsbcrd,zhbcrd,zgbcrd
                                           ! exbc.inc
!  REAL :: exbcbuf (exbcbufsz)

  REAL :: dx, dy, dz, umove, vmove, xgrdorg, ygrdorg, trulat1, trulat2,&
          trulon, sclfct, latitud, ctrlat, ctrlon, rdummy
                                           ! Header info.

  REAL, ALLOCATABLE :: x (:), y (:), z(:), zp (:,:,:), zpsoil(:,:,:)

  REAL, ALLOCATABLE :: xsm (:), ysm (:)

  REAL, ALLOCATABLE :: a3dlg (:,:,:), a3dsm (:, :, :)
  REAL, ALLOCATABLE :: a3dsoillg (:,:,:), a3dsoilsm (:, :, :)

  REAL, ALLOCATABLE :: a4dlg (:,:,:,:), a4dsm (:, :, :, :)

  REAL, ALLOCATABLE :: a2dxlg (:,:), a2dxsm (:, :)

  REAL, ALLOCATABLE :: a2dylg (:,:), a2dysm (:, :)

  REAL, ALLOCATABLE :: a2dxylg (:,:), a2dxysm (:,:)

  REAL, ALLOCATABLE :: a2dxy8lg (:,:,:), a2dxy8sm (:,:,:)

  REAL, ALLOCATABLE :: soil1lg (:,:,:), soil1sm (:,:,:)

  REAL, ALLOCATABLE :: soil2lg (:,:,:), soil2sm (:,:,:)

  REAL, ALLOCATABLE :: soil3lg (:,:,:), soil3sm (:,:,:)

  REAL, ALLOCATABLE :: soil4lg (:,:,:,:), soil4sm (:,:,:,:)

  INTEGER, ALLOCATABLE :: a2dilg(:,:), a2dism(:,:)
  INTEGER, ALLOCATABLE :: a3dilg(:,:,:), a3dism(:,:,:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  WRITE(6,*) "======== Beginning of executable code... ========="

  if ( mp_opt > 0 ) then
        write(6,*) 'splitstart:  not MP ready'
        call arpsstop('splitrestart:   not MP ready')
        return
  endif

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
  nzlg = nz
  nzsoillg = nzsoil

  ALLOCATE(ounit(nproc_x*nproc_y))
  ALLOCATE(ffi(nproc_x*nproc_y))
  ALLOCATE(ffj(nproc_x*nproc_y))

  ALLOCATE(x(nxlg))
  ALLOCATE(xsm(nx))
  ALLOCATE(y(nylg))
  ALLOCATE(ysm(ny))
  ALLOCATE(z(nzlg))
  ALLOCATE(zp(nxlg, nylg, nzlg))
  ALLOCATE(zpsoil(nxlg, nylg, nzsoillg))

  ALLOCATE(a3dlg(nxlg, nylg, nzlg))
  ALLOCATE(a3dsm(nx, ny, nz))
  ALLOCATE(a3dsoillg(nxlg, nylg, nzsoillg))
  ALLOCATE(a3dsoilsm(nx, ny, nzsoil))
  ALLOCATE(a4dlg(nxlg, nylg, nzlg, 5))
  ALLOCATE(a4dsm(nx, ny, nz, 5))
  ALLOCATE(a2dxlg(nxlg, nzlg))
  ALLOCATE(a2dxsm(nx, nz))
  ALLOCATE(a2dylg(nylg, nzlg))
  ALLOCATE(a2dysm(ny,   nz  ))

  ALLOCATE(a2dxylg(nxlg, nylg))
  ALLOCATE(a2dxysm(nx,   ny))
  ALLOCATE(a2dxy8lg(nxlg, nylg, 8))
  ALLOCATE(a2dxy8sm(nx,   ny,   8))
  ALLOCATE(soil1lg(nxlg, nylg, nstyps))
  ALLOCATE(soil1sm(nx,   ny,   nstyps))
  ALLOCATE(soil2lg(nxlg, nylg, 0:nstyps))
  ALLOCATE(soil2sm(nx,   ny,   0:nstyps))
  ALLOCATE(soil3lg(nxlg, nylg, 4))
  ALLOCATE(soil3sm(nx,   ny,   4))
  ALLOCATE(soil4lg(nxlg, nylg, nzsoillg, 0:nstyps))
  ALLOCATE(soil4sm(nx,   ny,   nzsoil,   0:nstyps))

  ALLOCATE(a2dilg(nxlg,nylg))
  ALLOCATE(a2dism(nx,  ny  ))
  ALLOCATE(a3dilg(nxlg,nylg,nstyps))
  ALLOCATE(a3dism(nx,  ny,  nstyps))

  lenstr = 0
  100   lenstr = lenstr + 1
  IF (fileheader(lenstr:lenstr) /= " ") GO TO 100
  lenstr = lenstr - 1

!
!-----------------------------------------------------------------------
!
!  Open the splitting file and splitted files.
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
!  ---   Since T3D processors only support COS and IEEE double precision
!        format, we have to translate the files into COS format.
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
!  Read in time and dimensions. And
!  Check against the dimensions passed to this subroutine.
!
!----------------------------------------------------------------------
!
!    WRITE(6,*) "==== Begin to read ==== filename: ", fileheader(1:lenstr)

    READ (10) ctime

    READ (10) nxin, nyin, nzin, nzsoilin

    IF ((nxin /= nxlg).OR.(nyin /= nylg).OR.(nzin /= nzlg).OR. &
        (nzsoilin /= nzsoillg)) THEN
      WRITE (*,*) "ERROR:  mismatch in sizes."
      WRITE (*,*) "nxin,nyin,nzin,nzsoilin: ",nxin,nyin,nzin,nzsoilin
      WRITE (*,*) "nxlg,nylg,nzlg,nzsoillg: ",nxlg,nylg,nzlg,nzsoillg

      CLOSE (10)
      DO ii=1+(jj-1)*maxunit,iiend
        CLOSE (ounit(ii))
      END DO
      GOTO 999
    END IF

!---------------------------------------------------------------------
!
!  Write time and dimensions for each processor.
!
!---------------------------------------------------------------------

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) ctime
      WRITE (ounit(ii)) nx, ny, nz, nzsoilin
    END DO

!
!-----------------------------------------------------------------------
!
!  Read/write header info.
!
!-----------------------------------------------------------------------
!
    READ (10)  basrstout,grdrstout,icerstout,sfcrstout,prcrstout,     &
                 rcumout,  exbcout,  mapfout,radrstout,    nstyp,     &
                 kfrsout,   idummy,  bmjsout,   idummy,   idummy,     &
                  idummy,                                             &
                  idummy,   idummy,   idummy,   idummy,   idummy

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii))                                                 &
               basrstout,grdrstout,icerstout,sfcrstout,prcrstout,       &
                 rcumout,  exbcout,  mapfout,radrstout,    nstyp,       &
                 kfrsout,   idummy,  bmjsout,   idummy,   idummy,       &
                  idummy,                                               &
                  idummy,   idummy,   idummy,   idummy,   idummy
    END DO

    READ (10)  ij(01), ij(02), ij(03), ij(04), ij(05), ij(06),          &
               ij(07), ij(08), ij(09), ij(10), ij(11), ij(12),          &
               ij(13), ij(14), ij(15), ij(16), ij(17), ij(18),          &
               ij(19), ij(20), ij(21), ij(22), ij(23), ij(24),          &
               ij(25), ij(26), ij(27), ij(28), ij(29), ij(30)

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii))                                                 &
               ij(01), ij(02), ij(03), ij(04), ij(05), ij(06),          &
               ij(07), ij(08), ij(09), ij(10), ij(11), ij(12),          &
               ij(13), ij(14), ij(15), ij(16), ij(17), ij(18),          &
               ij(19), ij(20), ij(21), ij(22), ij(23), ij(24),          &
               ij(25), ij(26), ij(27), ij(28), ij(29), ij(30)
    END DO

    READ (10) dx,dy,dz,umove,vmove,                                  &
              xgrdorg,ygrdorg,trulat1,trulat2,trulon,                &
              sclfct,latitud,ctrlat,ctrlon,rdummy,                   &
              rdummy,rdummy,rdummy,rdummy,rdummy

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE (ounit(ii)) dx,dy,dz,umove,vmove,                        &
              xgrdorg,ygrdorg,trulat1,trulat2,trulon,                &
              sclfct,latitud,ctrlat,ctrlon,rdummy,                   &
              rdummy,rdummy,rdummy,rdummy,rdummy
    END DO

    IF (grdrstout == 1) THEN
      READ(10) x
      READ(10) y
      READ(10) z
      READ(10) zp
      READ(10) zpsoil

      DO ii = 1 + (jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO i = 1, nx
           xsm(i) = x(i+(fi-1)*(nx-3))
        END DO
        DO j = 1, ny
          ysm(j) = y(j+(fj-1)*(ny-3))
        END DO
        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx
              a3dsm(i,j,k) = zp (i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        DO k = 1, nzsoil
          DO j = 1, ny
            DO i = 1, nx
              a3dsoilsm(i,j,k) = zpsoil (i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO

        WRITE (ounit(ii)) xsm
        WRITE (ounit(ii)) ysm
        WRITE (ounit(ii)) z
        WRITE (ounit(ii)) a3dsm
        WRITE (ounit(ii)) a3dsoilsm
      END DO

    END IF

!---------------------------------------------------------------------
!
!  3 dimensions variable READ and WRITE
!
!---------------------------------------------------------------------

!    WRITE(6,*) '======== Processing 3 dimensions variables ==============='

    nfields = 14                         ! 7 for tpast, 7 for tpresent
    IF (basrstout == 1) nfields = nfields + 6      ! ubar, var, ptbar
                                                   ! pbar, rhostr, qvbar
!    IF (icerstout /= 0) nfields = nfields + 6      ! 3 past 3 present
                                                   ! for qi, qs, qh
    DO fcnt = 1,30
      IF (ij(fcnt) > 0) nfields = nfields + 2
    END DO

    DO fcnt = 1,nfields

      READ (10) a3dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              a3dsm(i,j,k) =                                         &
                  a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO

        WRITE (ounit(ii)) a3dsm
      END DO

    END DO

!---------------------------------------------------------------------
!
!  2 dimensions variables
!
!---------------------------------------------------------------------

!    WRITE(6,*) '======== Processing 2 dimensions variables ==============='

  nfields = 2
  DO fcnt = 1, 2

    DO nfields =1, 2            ! udteb, uetwb; pdteb, pdtwb
      READ (10) a2dylg
      DO ii = 1+(jj-1)*maxunit, iiend
        fj = ffj(ii)
        DO k = 1, nz
          DO j =  1, ny
            a2dysm(j,k) = a2dylg (j+(fj-1)*(ny-3), k)
          END DO
        END DO
        WRITE (ounit(ii)) a2dysm
      END DO
    END DO

    DO nfields =1, 2           ! vdtnb, vdtsb; pdtnb, pdtsb
      READ (10) a2dxlg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        DO k = 1, nz
          DO i =  1,nx
            a2dxsm(i,k) = a2dxlg (i+(fi-1)*(nx-3), k)
          END DO
        END DO
        WRITE (ounit(ii)) a2dxsm
      END DO
    END DO

  END DO                      ! fcnt

!------------------------------------------------------------------
!
!  sfc/soil variables
!
!------------------------------------------------------------------


!    WRITE(6,*) '======== Processing sfc/soil variables ==============='

  IF (sfcrstout /= 0) THEN

!    DO fcnt = 1, 2

      READ (10) a3dilg                   ! soiltyp
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1, nstyps
          DO j = 1, ny
            DO i = 1, nx
              a3dism(i, j, k) = a3dilg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a3dism
      END DO

      READ (10) soil1lg                    ! stypfrct
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1, nstyps
          DO j = 1, ny
            DO i = 1, nx
              soil1sm(i, j, k) = soil1lg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) soil1sm
      END DO

!    END DO  ! fcnt

    READ (10) a2dilg              ! vegtyp
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dism(i, j) = a2dilg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dism
    END DO

    DO fcnt = 1, 3       ! lai, roufns, veg
      READ (10) a2dxylg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dxysm(i, j) = a2dxylg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dxysm
      END DO
    END DO  ! fcnt


    ! qvsfc
    READ (10) soil2lg
    DO ii = 1+(jj-1)*maxunit, iiend
      fi = ffi(ii)
      fj = ffj(ii)
      DO k = 0, nstyps
        DO j = 1, ny
          DO i = 1, nx
            soil2sm(i, j, k) = soil2lg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
          END DO
        END DO
      END DO
      WRITE (ounit(ii)) soil2sm
    END DO

    ! tsoil, qsoil
    DO fcnt = 1, 2
      READ (10) soil4lg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO l = 0, nstyps
            DO k = 1, nzsoil
            DO j = 1, ny
              DO i = 1, nx
                soil4sm(i, j, k, l) = &
                  soil4lg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k, l)
                END DO
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) soil4sm
      END DO
    END DO  ! fcnt

    ! wetcanp
    READ (10) soil2lg
    DO ii = 1+(jj-1)*maxunit, iiend
      fi = ffi(ii)
      fj = ffj(ii)
      DO k = 0, nstyps
        DO j = 1, ny
          DO i = 1, nx
            soil2sm(i, j, k) = soil2lg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
          END DO
        END DO
      END DO
      WRITE (ounit(ii)) soil2sm
    END DO

    !snowdpth
    READ (10) a2dxylg
    DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dxysm(i, j) = a2dxylg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dxysm
    END DO
  END IF
!---------------------------------------------------------------------
!
!  prcrstout
!
!---------------------------------------------------------------------

!    WRITE(6,*) '======== Processing prcrstout    variables ==============='

  IF (prcrstout /= 0) THEN

    DO fcnt = 1, 2
      READ (10) a2dxylg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dxysm(i, j) = a2dxylg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dxysm
      END DO
    END DO  ! fcnt

    READ (10) soil3lg
    DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1, 4
          DO j = 1, ny
            DO i = 1, nx
              soil3sm(i, j, k) = soil3lg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) soil3sm
    END DO

   END IF
!---------------------------------------------------------------------
!
!  rcumout
!
!---------------------------------------------------------------------


!    WRITE(6,*) '======== Processing rcumout      variables ==============='

  IF (rcumout /= 0) THEN

      READ (10) a3dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              a3dsm(i,j,k) =                                         &
                  a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a3dsm
      END DO

      READ (10) a4dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO l = 1, 5
          DO k = 1,nz
            DO j = 1,ny
              DO i = 1,nx
                a4dsm(i,j,k,l) =                                     &
                  a4dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k, l)
              END DO
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a4dsm
      END DO

  END IF
!------------------------------------------------------------------
!
!  exbcout /= 0
!
!-----------------------------------------------------------------


!    WRITE(6,*) '======== Processing exbcout variables ==============='

  IF (exbcout /= 0) THEN
    READ (10) abstfcst0, abstfcst,                                 &
              ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                     &
              qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd,    &
                   ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd,      &
                   zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd

    DO ii=1+(jj-1)*maxunit,iiend
      WRITE(ounit(ii)) abstfcst0, abstfcst,                        &
              ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                     &
              qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd,    &
                   ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd,      &
                   zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd
    END DO

!    READ (10) exbcbuf

    DO var = 1, exbcvarsz
      READ (10) a3dlg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              a3dsm(i,j,k) =                                        &
                  a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a3dsm
      END DO

    END DO  ! var

  END IF

!------------------------------------------------------------------
!
!  mapfout == 1
!
!-----------------------------------------------------------------


!    WRITE(6,*) '======== Processing mapfout variables ==============='

  IF (mapfout == 1) THEN

      READ (10) a2dxy8lg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1, 8
          DO j = 1, ny
            DO i = 1, nx
               a2dxy8sm(i, j, k) = a2dxy8lg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a2dxy8sm
      END DO

  END IF

!------------------------------------------------------------------
!
!  radrstout == 1
!
!-----------------------------------------------------------------


!    WRITE(6,*) '======== Processing radsrtout variables ==============='

  IF (radrstout ==1 ) THEN

      READ (10) a3dlg                   ! radfrc

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              a3dsm(i,j,k) =                                        &
                  a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a3dsm
      END DO

    DO fcnt = 1, 4                    ! radsw, rnflx, radswnet, radlwin
      READ (10) a2dxylg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dxysm(i, j) = a2dxylg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dxysm
      END DO    ! ii
    END DO    ! fcnt

  END IF

!------------------------------------------------------------------
!
!  kfrsout  /= 0
!
!-----------------------------------------------------------------


!    WRITE(6,*) '======== Processing kfrsout variables ==============='

  IF (kfrsout /= 0 ) THEN

      READ (10) a3dlg                   ! w0avg

      DO ii=1+(jj-1)*maxunit,iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              a3dsm(i,j,k) =                                        &
                  a3dlg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k)
            END DO
          END DO
        END DO
        WRITE (ounit(ii)) a3dsm
      END DO

      READ (10) a2dilg                ! nca
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dism(i, j) = a2dilg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dism
      END DO    ! ii

      READ (10) a2dxylg               ! kfraincv
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dxysm(i, j) = a2dxylg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dxysm
      END DO    ! ii

  END IF

!EMK BMJ
!------------------------------------------------------------------
!
!  bmjsout  /= 0
!
!-----------------------------------------------------------------


!    WRITE(6,*) '======== Processing bmjsout variables ==============='

  IF (bmjsout /= 0 ) THEN

    DO fcnt = 1, 3                    ! cldefi, xland, bmjraincv
      READ (10) a2dxylg
      DO ii = 1+(jj-1)*maxunit, iiend
        fi = ffi(ii)
        fj = ffj(ii)
        DO j = 1, ny
          DO i = 1, nx
             a2dxysm(i, j) = a2dxylg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3))
          END DO
        END DO
        WRITE (ounit(ii)) a2dxysm
      END DO    ! ii
    END DO    ! fcnt

  END IF

!------------------------------------------------------------------
!
!  Close files
!
!-----------------------------------------------------------------

!    WRITE(6,*) '======== Successful !!!!!!!!!!!!!!!!! ==============='

    CLOSE (10)
    DO ii=1+(jj-1)*maxunit,iiend
      CLOSE (ounit(ii))
    END DO

  END DO    ! jj

  999 CONTINUE

!-------------------------------------------------------------------
!
! Deallocate the variables before return
!
!-------------------------------------------------------------------

  DEALLOCATE(ounit)
  DEALLOCATE(ffi)
  DEALLOCATE(ffj)

  DEALLOCATE(x)
  DEALLOCATE(xsm)
  DEALLOCATE(y)
  DEALLOCATE(ysm)
  DEALLOCATE(z)
  DEALLOCATE(zp)
  DEALLOCATE(zpsoil)
  DEALLOCATE(a3dlg)
  DEALLOCATE(a3dsm)
  DEALLOCATE(a3dsoillg)
  DEALLOCATE(a3dsoilsm)
  DEALLOCATE(a4dlg)
  DEALLOCATE(a4dsm)
  DEALLOCATE(a2dxlg)
  DEALLOCATE(a2dxsm)
  DEALLOCATE(a2dylg)
  DEALLOCATE(a2dysm)

  DEALLOCATE(a2dxylg)
  DEALLOCATE(a2dxysm)
  DEALLOCATE(soil1lg)
  DEALLOCATE(soil1sm)
  DEALLOCATE(soil2lg)
  DEALLOCATE(soil2sm)
  DEALLOCATE(soil3lg)
  DEALLOCATE(soil3sm)
  DEALLOCATE(soil4lg)
  DEALLOCATE(soil4sm)

  DEALLOCATE(a2dilg, a2dism, a3dilg, a3dism)

  RETURN

END SUBROUTINE splitrestart

