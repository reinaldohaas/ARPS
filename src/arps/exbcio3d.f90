!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READEXBC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE readexbc(nx,ny,nz,filename,lfname,ctime,               &
                    u,v,w,pt,pr,qv,qscalar, tem1,ierr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in 6 primary fields for use as external boundary
!  conditions.
!
!  In general these data will come from another model.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  May, 1994
!
!  MODIFICATION HISTORY:
!
!  5/26/94 (Yuhe Liu)
!  Merged into the part of ARPS for external boundary conditions.
!
!  8/8/95 (M. Xue)
!  Added water and ice variables to the EXBC files.
!  To read earlier version EXBC files, one has to set old_v=1.
!
! 12/7/10 (B. Putnam)
!  Added code to remove initial hail and add it to the graupel catagory
!  for the option to turn of hail.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction (east/west)
!  ny       Number of grid points in the y-direction (north/south)
!  nz       Number of grid points in the vertical
!
!  dx       Expected x grid spacing
!  dy       Expected y grid spacing
!  dz       Expected z grid spacing
!
!  ctrlat   Expected center latitude
!  ctrlon   Expected center longitude
!
!  filename File name of EXBC boundary data set.
!  lfname   Length of the filename
!
!  OUTPUT:
!
!  ctime    Charater representation of the time of EXBC data
!
!  ubcrd    Flag indicating (1) if the u  field is valid
!  vbcrd    Flag indicating (1) if the v  field is valid
!  wbcrd    Flag indicating (1) if the w  field is valid
!  ptbcrd   Flag indicating (1) if the pt field is valid
!  prbcrd   Flag indicating (1) if the pr field is valid
!  qvbcrd   Flag indicating (1) if the qv field is valid
!  qcbcrd   Flag indicating (1) if the qc field is valid
!  qrbcrd   Flag indicating (1) if the qr field is valid
!  qibcrd   Flag indicating (1) if the qi field is valid
!  qsbcrd   Flag indicating (1) if the qs field is valid
!  qhbcrd   Flag indicating (1) if the qh field is valid
!
!  u        x component of velocity (m/s)
!  v        y component of velocity (m/s)
!  w        Vertical component of Cartesian velocity (m/s)
!  pt       Potential temperature (K)
!  pr       Pressure (Pascal)
!  qv       Water vapor mixing ratio humidity (kg/kg)
!  qc       Cloud water mixing ratio humidity (kg/kg)
!  qr       Rain water mixing ratio humidity (kg/kg)
!  qi       Cloud ice mixing ratio humidity (kg/kg)
!  qs       Snow mixing ratio humidity (kg/kg)
!  qh       Hail water mixing ratio humidity (kg/kg)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'            ! mpi parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.
  CHARACTER (LEN=* ) :: filename
  INTEGER :: lfname
  CHARACTER (LEN=* ) :: ctime

  REAL :: u(nx,ny,nz)       ! u-velocity (m/s)
  REAL :: v(nx,ny,nz)       ! v-velocity (m/s)
  REAL :: w(nx,ny,nz)       ! w-velocity (m/s)
  REAL :: pt(nx,ny,nz)      ! Potential temperature (K)
  REAL :: pr(nx,ny,nz)      ! Pressure (Pascal)
  REAL :: qv(nx,ny,nz)      ! Water vapor mixing ratio (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin,nzin
  REAL    :: ctrlatin,ctrlonin

  INTEGER :: istat, ierr, idummy, old_v
  REAL    :: amax, amin

  INTEGER :: ireturn
  INTEGER :: strhoptin
  REAL    :: dxin,dyin,dzin,dzminin,zrefsfcin,dlayer1in,                &
             dlayer2in,zflatin,strhtunein
  REAL    :: trulat1in,trulat2in,trulonin,sclfctin
  INTEGER :: maprojin

  INTEGER(2), ALLOCATABLE :: itmp(:,:,:)      ! Temporary array
  REAL,       ALLOCATABLE :: hmax(:), hmin(:) ! Temporary array

  INTEGER :: clipxy, clipz

  INTEGER :: i, j, k

  INTEGER :: nq, nqin, nscalarin

  INTEGER :: qcbcrd, qrbcrd, qibcrd, qsbcrd, qhbcrd, qgbcrd,            &
             ncbcrd, nrbcrd, nibcrd, nsbcrd, nhbcrd, ngbcrd,            &
                     zrbcrd, zibcrd, zsbcrd, zhbcrd, zgbcrd,            &
             ccbcrd

  INTEGER, SAVE :: bfid
  INTEGER, SAVE :: itime = 0
  INTEGER       :: istop

  CHARACTER(LEN=40) :: fmtverin

  CHARACTER(LEN=4)  :: upcase

  REAL, ALLOCATABLE :: var3du(:,:,:)
  REAL, ALLOCATABLE :: var3dv(:,:,:)
  REAL, ALLOCATABLE :: var3dw(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (exbcfmt == 3) THEN
    ALLOCATE (itmp(nx,ny,nz),stat=istat)
    CALL check_alloc_status(istat, "readexbc:itmp")

    ALLOCATE (hmax(nz),stat=istat)
    CALL check_alloc_status(istat, "readexbc:hmax")

    ALLOCATE (hmin(nz),stat=istat)
    CALL check_alloc_status(istat, "readexbc:hmin")
  END IF

  IF (myproc == 0) WRITE (6,'(2a)')         &
          'READEXBC: reading in external boundary data from file ',     &
          filename(1:lfname)

  IF (size(qscalarbcrd) < nscalar) THEN
    IF (myproc == 0) WRITE (6,'(a,/,a,I2)')                             &
                'External boundary flag arrays is too small in size.',  &
                'nscalar = ',nscalar
    CALL arpsstop('qscalarbcrd too small',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------
!
  IF (exbcfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    CALL getunit( bfid )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(filename(1:lfname), '-F f77 -N ieee', ierr)

    OPEN(bfid,FILE=filename(1:lfname),STATUS='old',                &
              FORM='unformatted',IOSTAT=istat)

    IF ( istat /= 0 ) THEN
      ierr = 1
      GO TO 900
    END IF

    READ(bfid,ERR=999) nxin,nyin,nzin,dxin,dyin,dzin,                   &
                       ctrlatin,ctrlonin

    READ(bfid,ERR=999) ctime

    old_v = 0 ! In case that the EXBC files are of an earlier
              ! version that does not contain water and ice variables,
              ! set old_v to 1. Otherwise, set it to 0.

    IF( old_v == 1 ) THEN

      READ(bfid,ERR=999) ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd

      qcbcrd=0
      qrbcrd=0
      qibcrd=0
      qsbcrd=0
      qhbcrd=0

    ELSE

      READ(bfid,ERR=999) ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,               &
                      qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,               &
                      qhbcrd,qgbcrd,ncbcrd,nrbcrd,nibcrd,               &
                      nsbcrd,ngbcrd,nhbcrd,zrbcrd,zibcrd,               &
                      zsbcrd,zgbcrd,zhbcrd,ccbcrd,idummy,               &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,nscalarin

    END IF

    qscalarbcrd(:) = 0
    IF (P_QC > 0) qscalarbcrd(P_QC) = qcbcrd
    IF (P_QR > 0) qscalarbcrd(P_QR) = qrbcrd
    IF (P_QI > 0) qscalarbcrd(P_QI) = qibcrd
    IF (P_QS > 0) qscalarbcrd(P_QS) = qsbcrd
    IF (P_QH > 0) qscalarbcrd(P_QH) = qhbcrd
    IF (P_QG > 0) qscalarbcrd(P_QG) = qgbcrd

    IF (P_NC > 0) qscalarbcrd(P_NC) = ncbcrd
    IF (P_NR > 0) qscalarbcrd(P_NR) = nrbcrd
    IF (P_NI > 0) qscalarbcrd(P_NI) = nibcrd
    IF (P_NS > 0) qscalarbcrd(P_NS) = nsbcrd
    IF (P_NH > 0) qscalarbcrd(P_NH) = nhbcrd
    IF (P_NG > 0) qscalarbcrd(P_NG) = ngbcrd

    IF (P_ZR > 0) qscalarbcrd(P_ZR) = zrbcrd
    IF (P_ZI > 0) qscalarbcrd(P_ZI) = zibcrd
    IF (P_ZS > 0) qscalarbcrd(P_ZS) = zsbcrd
    IF (P_ZH > 0) qscalarbcrd(P_ZH) = zhbcrd
    IF (P_ZG > 0) qscalarbcrd(P_ZG) = zgbcrd

    IF (P_CC > 0) qscalarbcrd(P_CC) = ccbcrd
  ELSE IF (exbcfmt == 3) THEN

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(filename(1:lfname)), 1, bfid)
    IF (bfid < 0) THEN
      WRITE (6,*) "READEXBC: ERROR opening ",                           &
                 trim(filename(1:lfname))," for reading."
      ierr = 1
      GO TO 900
    END IF

    CALL hdfrdc(bfid,40,"fmtver",fmtverin,istat)

                  ! Should check fmtverin here

    CALL hdfrdc(bfid,15,"ctime",ctime,istat)

    CALL hdfrdi(bfid,"nx",nxin,istat)
    CALL hdfrdi(bfid,"ny",nyin,istat)
    CALL hdfrdi(bfid,"nz",nzin,istat)
    CALL hdfrdr(bfid,"dx",dxin,istat)
    CALL hdfrdr(bfid,"dy",dyin,istat)
    CALL hdfrdr(bfid,"dz",dzin,istat)
    CALL hdfrdr(bfid,"dzmin",dzminin,istat)
    CALL hdfrdi(bfid,"strhopt",strhoptin,istat)
    CALL hdfrdr(bfid,"zrefsfc",zrefsfcin,istat)
    CALL hdfrdr(bfid,"dlayer1",dlayer1in,istat)
    CALL hdfrdr(bfid,"dlayer2",dlayer2in,istat)
    CALL hdfrdr(bfid,"zflat",zflatin,istat)
    CALL hdfrdr(bfid,"strhtune",strhtunein,istat)
    CALL hdfrdi(bfid,"mapproj",maprojin,istat)
    CALL hdfrdr(bfid,"trulat1",trulat1in,istat)
    CALL hdfrdr(bfid,"trulat2",trulat2in,istat)
    CALL hdfrdr(bfid,"trulon",trulonin,istat)
    CALL hdfrdr(bfid,"sclfct",sclfctin,istat)
    CALL hdfrdr(bfid,"ctrlat",ctrlatin,istat)
    CALL hdfrdr(bfid,"ctrlon",ctrlonin,istat)

    CALL hdfrdi(bfid,"ubcflg",ubcrd,istat)
    CALL hdfrdi(bfid,"vbcflg",vbcrd,istat)
    CALL hdfrdi(bfid,"wbcflg",wbcrd,istat)
    CALL hdfrdi(bfid,"ptbcflg",ptbcrd,istat)
    CALL hdfrdi(bfid,"prbcflg",prbcrd,istat)
    CALL hdfrdi(bfid,"qvbcflg",qvbcrd,istat)

    CALL hdfrdi(bfid,"nscalar",nscalarin,istat)
    DO nq = 1,nscalar
      CALL hdfrdi(bfid,TRIM(qnames(nq))//'bcflg',qscalarbcrd(nq),istat)
    END DO

    CALL hdfrdi(bfid, 'clipxy', clipxy,istat)
    IF (istat == 0 .AND. clipxy < ngbrz) THEN
      WRITE (6,*) "READEXBC: ERROR, clipxy (ngbrz) in exbc file too small"
      ierr = 1
      GO TO 900
    END IF
    CALL hdfrdi(bfid, 'clipz', clipz,istat)
    IF (istat == 0 .AND. clipz > rayklow) THEN
      WRITE (6,*) "READEXBC: ERROR, clipz (rayklow) in exbc file too large"
      WRITE(6,*) ' rayklow = ',rayklow,', clipz = ',clipz
      ierr = 1
      GO TO 900
    END IF

  ELSE IF (exbcfmt == 7 .OR. exbcfmt == 8) THEN

!-----------------------------------------------------------------------
!
!  NetCDF format
!
!-----------------------------------------------------------------------

    IF (exbcfmt == 7) THEN
      itime = 1
    ELSE
      itime = itime + 1
      istop = NINT( (tstop-tstart)/thisdmp ) + 1
    END IF

    IF (itime == 1)      &
      CALL netopen(TRIM(filename(1:lfname)), 'R', bfid)

      CALL net_get_exbc(bfid,nxin,nyin,nzin,itime,dxin,dyin,dzin,       &
                 dzminin,strhoptin,zrefsfcin,dlayer1in,dlayer2in,       &
                 zflatin,strhtunein,maprojin,sclfctin,trulat1in,        &
                 trulat2in,trulonin,ctrlatin,ctrlonin,                  &
                 ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd,                &
                 nscalarin,qscalarbcrd,ctime,istat)

  ELSE

    ! alternate exbc format ...
    WRITE(6,'(1x,3a)') 'The supported exbc data format are ',           &
               'binary (exbcfmt=1), HDF4 (exbcfmt = 3). ',              &
               'and NetCDF (exbcfmt = 7).'
    CALL arpsstop('Exbc data format is not supported.',1)

  END IF

!-----------------------------------------------------------------------
!
!  Check the data file for consistent grid parameters.
!
!-----------------------------------------------------------------------

  IF (exbcfmt == 1) THEN

    IF (myproc == 0) WRITE (6,*)                                       &
        "READEXBC: WARNING, not checking all map projection parameters"

    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
            dx,dy,ctrlat,ctrlon,                                        &
            mapproj,trulat1,trulat2,trulon,sclfct,                      &
            dxin,dyin,ctrlatin,ctrlonin,                                &
            mapproj,trulat1,trulat2,trulon,sclfct,ireturn)

  ELSE

    CALL checkgrid3d(nx,ny,nz,nxin,nyin,nzin,                           &
            dx,dy,dz,dzmin,ctrlat,ctrlon,                               &
            strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,             &
            mapproj,trulat1,trulat2,trulon,sclfct,                      &
            dxin,dyin,dzin,dzminin,ctrlatin,ctrlonin,                   &
            strhoptin,zrefsfcin,dlayer1in,dlayer2in,zflatin,strhtunein, &
            maprojin,trulat1in,trulat2in,trulonin,sclfctin,ireturn)

  END IF

  IF (ireturn /= 0) THEN
    WRITE (6,*) "READEXBC: ERROR, grid parameter mismatch"
    ierr = 1
    GO TO 900
  END IF

!-----------------------------------------------------------------------
!
!  Read in the external boundary file data
!
!-----------------------------------------------------------------------

  IF (exbcfmt == 1) THEN

    READ(bfid,ERR=999) u
    READ(bfid,ERR=999) v
    READ(bfid,ERR=999) w
    READ(bfid,ERR=999) pt
    READ(bfid,ERR=999) pr

    IF(qvbcrd == 1) READ(bfid,ERR=999) qv

    IF (nscalarin > 0 ) THEN
      DO nqin = 1,nscalarin
        READ(bfid,ERR=999) tem1
        DO nq = 1,nscalar
          IF (nqin == qscalarbcrd(nq)) THEN
            qscalar(:,:,:,nq) = tem1(:,:,:)
            EXIT
          END IF
        END DO
      END DO
      ! DTD: Added call to init_MM to make sure that the
      ! number concentration and/or reflectivity fields are
      ! initialized if only mixing ratio fields are read in
      ! and the user is running the multi-moment scheme.
      ! Otherwise, clouds will be "eaten away" as they move
      ! in from the external boundary.
      IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
        WRITE(6,*) 'Diagnosing additional moments at external boundary for multi-moment scheme'
        CALL init_MM_exbc(nx,ny,nz,pt,pr,qscalar)
        ! Also update qscalarbcrd flags to reflect the newly diagnosed additional moments
        ! even though they were not actually read in from the external boundary file
        DO nqin = 1,nscalarin
          DO nq = 1,nscalar
            IF(nqin == qscalarbcrd(nq)) THEN
              qscalarbcrd(nq) = nq
              qscalarbcrd(nq+6) = nq+6
              IF(mphyopt == 11 .and. nq /= P_QC) THEN
                qscalarbcrd(nq+11) = nq+11
              END IF
              EXIT
            END IF
          END DO
        END DO
      END IF
    ELSE              ! earlier version
      IF(qcbcrd == 1) THEN
        READ(bfid,ERR=999) tem1
        IF (P_QC > 0) qscalar(:,:,:,P_QC) = tem1(:,:,:)
      END IF
      IF(qrbcrd == 1) THEN
        READ(bfid,ERR=999) tem1
        IF (P_QR > 0) qscalar(:,:,:,P_QR) = tem1(:,:,:)
      END IF
      IF(qibcrd == 1) THEN
        READ(bfid,ERR=999) tem1
        IF (P_QI > 0) qscalar(:,:,:,P_QI) = tem1(:,:,:)
      END IF
      IF(qsbcrd == 1) THEN
        READ(bfid,ERR=999) tem1
        IF (P_QS > 0) qscalar(:,:,:,P_QS) = tem1(:,:,:)
      END IF
      IF(qhbcrd == 1) THEN
        READ(bfid,ERR=999) tem1
        IF (P_QH > 0) qscalar(:,:,:,P_QH) = tem1(:,:,:)
      END IF
    END IF


  ELSE IF (exbcfmt == 3) THEN

    IF (ubcrd == 1) THEN
      CALL hdfrd3d(bfid,"u",nx,ny,nz,u,istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 999
    END IF
    IF (vbcrd == 1) THEN
      CALL hdfrd3d(bfid,"v",nx,ny,nz,v,istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 999
    END IF
    IF (wbcrd == 1) THEN
      CALL hdfrd3d(bfid,"w",nx,ny,nz,w,istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 999
    END IF
    IF (ptbcrd == 1) THEN
      CALL hdfrd3d(bfid,"pt",nx,ny,nz,pt,istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 999
    END IF
    IF (prbcrd == 1) THEN
      CALL hdfrd3d(bfid,"p",nx,ny,nz,pr,istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 999
    END IF
    IF (qvbcrd == 1) THEN
      CALL hdfrd3d(bfid,"qv",nx,ny,nz,qv,istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 999
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcrd(nq) > 0) THEN
        CALL hdfrd3d(bfid,qnames(nq),nx,ny,nz,qscalar(:,:,:,nq),        &
                     istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
    END DO
    IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
      WRITE(6,*) 'Diagnosing additional moments at external boundary for multi-moment scheme'
      CALL init_MM_exbc(nx,ny,nz,pt,pr,qscalar)
      ! Also update qscalarbcrd flags to reflect the newly diagnosed additional moments
      ! even though they were not actually read in from the external boundary file
      DO nqin = 1,nscalarin
        DO nq = 1,nscalar
          IF(nqin == qscalarbcrd(nq)) THEN
            qscalarbcrd(nq) = nq
            qscalarbcrd(nq+6) = nq+6
            IF(mphyopt == 11 .and. nq /= P_QC) THEN
              qscalarbcrd(nq+11) = nq+11
            END IF
            EXIT
          END IF
        END DO
      END DO
      DO nq = 1,nscalar
        print*,'nq,qscalarbcrd(nq)',nq,qscalarbcrd(nq)
      END DO
    END IF
  ELSE IF (exbcfmt == 7 .OR. exbcfmt == 8) THEN

    ALLOCATE(var3du(nx,  ny-1,nz-1), STAT = istat)
    ALLOCATE(var3dv(nx-1,ny,  nz-1), STAT = istat)
    ALLOCATE(var3dw(nx-1,ny-1,nz  ), STAT = istat)

    IF (ubcrd == 1) THEN
      CALL netread3d(bfid,0,itime,"U",nx,ny-1,nz-1,var3du)
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx
            u(i,j,k) = var3du(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1)
    END IF

    IF (vbcrd == 1) THEN
      CALL netread3d(bfid,0,itime,"V",nx-1,ny,nz-1,var3dv)
      DO k = 1,nz-1
        DO j = 1,ny
          DO i = 1,nx-1
            v(i,j,k) = var3dv(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1)
    END IF

    IF (wbcrd == 1) THEN
      CALL netread3d(bfid,0,itime,"W",nx-1,ny-1,nz,var3dw)
      DO k = 1,nz
        DO j = 1,ny-1
          DO i = 1,nx-1
            w(i,j,k) = var3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz)
    END IF

    IF (ptbcrd == 1) THEN
      CALL netread3d(bfid,0,itime,"PT",nx-1,ny-1,nz-1,var3dw)
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            pt(i,j,k) = var3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(pt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    END IF

    IF (prbcrd == 1) THEN
      CALL netread3d(bfid,0,itime,"P",nx-1,ny-1,nz-1,var3dw)
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            pr(i,j,k) = var3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(pr,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    END IF

    IF (qvbcrd == 1) THEN
      CALL netread3d(bfid,0,itime,"QV",nx-1,ny-1,nz-1,var3dw)
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            qv(i,j,k) = var3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(qv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcrd(nq) > 0 ) THEN
        CALL netread3d(bfid,0,itime,upcase(qnames(nq)),nx-1,ny-1,nz-1,var3dw)
        DO k = 1,nz-1
          DO j = 1,ny-1
            DO i = 1,nx-1
              qscalar(i,j,k,nq) = var3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      END IF
    END DO
    IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
      WRITE(6,*) 'Diagnosing additional moments at external boundary for multi-moment scheme'
      CALL init_MM_exbc(nx,ny,nz,pt,pr,qscalar)
      ! Also update qscalarbcrd flags to reflect the newly diagnosed additional moments
      ! even though they were not actually read in from the external boundary file
      DO nqin = 1,nscalarin
        DO nq = 1,nscalar
          IF(nqin == qscalarbcrd(nq)) THEN
            qscalarbcrd(nq) = nq
            qscalarbcrd(nq+6) = nq+6
            IF(mphyopt == 11 .and. nq /= P_QC) THEN
              qscalarbcrd(nq+11) = nq+11
            END IF
            EXIT
          END IF
        END DO
      END DO
    END IF
  END IF

  IF(myproc == 0)THEN

    write(6,'(/1x,a/)') 'Max. and Min. of EXBC data variables:'

    IF (ubcrd == 1) THEN
    CALL a3dmax0lcl(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'umin = ', amin,',  umax =',amax
    END IF

    IF (vbcrd == 1) THEN
    CALL a3dmax0lcl(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'vmin = ', amin,',  vmax =',amax
    END IF

    IF (wbcrd == 1) THEN
    CALL a3dmax0lcl(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'wmin = ', amin,',  wmax =',amax
    END IF

    IF (ptbcrd == 1) THEN
    CALL a3dmax0lcl(pt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'ptmin= ', amin,',  ptmax=',amax
    END IF

    IF (prbcrd == 1) THEN
    CALL a3dmax0lcl(pr,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'pmin = ', amin,',  pmax =',amax
    END IF

    IF (qvbcrd == 1) THEN
    CALL a3dmax0lcl(qv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'qvmin= ', amin,',  qvmax=',amax
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcrd(nq) > 0) THEN
      CALL a3dmax0lcl(qscalar(:,:,:,nq),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min= ', amin,       &
                          ',  '//TRIM(qnames(nq))//'max=',amax
      END IF
    END DO

  END IF

  ierr = 0

  GO TO 900

  999   CONTINUE

  WRITE (6,*) "READEXBC: ERROR reading data ",                          &
              "from file ",trim(filename(1:lfname))," returning"
  ierr = 2

  900   CONTINUE

  IF (exbcfmt == 1) THEN

    CLOSE (bfid)
    CALL retunit( bfid )

  ELSE IF (exbcfmt == 3) THEN

    CALL hdfclose(bfid,istat)

    DEALLOCATE (itmp,stat=istat)
    DEALLOCATE (hmax,stat=istat)
    DEALLOCATE (hmin,stat=istat)
    ! alternate dump format ...

  ELSE IF (exbcfmt == 7 .OR. exbcfmt ==8 ) THEN

    IF (exbcfmt == 7 .OR. itime >= istop) CALL netclose(bfid)
    DEALLOCATE(var3du,var3dv,var3dw)

  END IF

  RETURN
END SUBROUTINE readexbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READSPLITEXBC              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readsplitexbc(nx,ny,nz,filename,lfname, ctime,               &
                         u,v,w,pt,pr,qv,qscalar, ierr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in external boundary fields, and Split and scatter to MP
!  processes from the root process.
!
!  In general these data will come from another model.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/30/2002
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction (east/west)
!  ny       Number of grid points in the y-direction (north/south)
!  nz       Number of grid points in the vertical
!
!  dx       Expected x grid spacing
!  dy       Expected y grid spacing
!  dz       Expected z grid spacing
!
!  ctrlat   Expected center latitude
!  ctrlon   Expected center longitude
!
!  filename File name of EXBC boundary data set.
!  lfname   Length of the filename
!
!  OUTPUT:
!
!  ctime    Charater representation of the time of EXBC data
!
!  ubcrd    Flag indicating (1) if the u  field is valid
!  vbcrd    Flag indicating (1) if the v  field is valid
!  wbcrd    Flag indicating (1) if the w  field is valid
!  ptbcrd   Flag indicating (1) if the pt field is valid
!  prbcrd   Flag indicating (1) if the pr field is valid
!  qvbcrd   Flag indicating (1) if the qv field is valid
!  qcbcrd   Flag indicating (1) if the qc field is valid
!  qrbcrd   Flag indicating (1) if the qr field is valid
!  qibcrd   Flag indicating (1) if the qi field is valid
!  qsbcrd   Flag indicating (1) if the qs field is valid
!  qhbcrd   Flag indicating (1) if the qh field is valid
!
!  u        x component of velocity (m/s)
!  v        y component of velocity (m/s)
!  w        Vertical component of Cartesian velocity (m/s)
!  pt       Potential temperature (K)
!  pr       Pressure (Pascal)
!  qv       Water vapor mixing ratio humidity (kg/kg)
!  qc       Cloud water mixing ratio humidity (kg/kg)
!  qr       Rain water mixing ratio humidity (kg/kg)
!  qi       Cloud ice mixing ratio humidity (kg/kg)
!  qs       Snow mixing ratio humidity (kg/kg)
!  qh       Hail water mixing ratio humidity (kg/kg)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'            ! mpi parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER             :: nx,ny,nz  ! Number of grid points in x, y, and z dir.
  CHARACTER (LEN=* )  :: filename
  INTEGER             :: lfname
  CHARACTER (LEN=* )  :: ctime

  REAL :: u(nx,ny,nz)       ! u-velocity (m/s)
  REAL :: v(nx,ny,nz)       ! v-velocity (m/s)
  REAL :: w(nx,ny,nz)       ! w-velocity (m/s)
  REAL :: pt(nx,ny,nz)      ! Potential temperature (K)
  REAL :: pr(nx,ny,nz)      ! Pressure (Pascal)
  REAL :: qv(nx,ny,nz)      ! Water vapor mixing ratio (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin,nzin
  REAL    :: ctrlatin,ctrlonin

  INTEGER :: istat, ierr, idummy, old_v
  REAL    :: amax, amin

  INTEGER :: ireturn
  INTEGER :: strhoptin
  REAL    :: dxin,dyin,dzin,dzminin,zrefsfcin,dlayer1in,                &
             dlayer2in,zflatin,strhtunein
  REAL    :: trulat1in,trulat2in,trulonin,sclfctin
  INTEGER :: maprojin

  INTEGER :: nq, nqin, nscalarin
  INTEGER :: qcbcrd, qrbcrd, qibcrd, qsbcrd, qhbcrd, qgbcrd,            &
             ncbcrd, nrbcrd, nibcrd, nsbcrd, nhbcrd, ngbcrd,            &
                     zrbcrd, zibcrd, zsbcrd, zhbcrd, zgbcrd,            &
             ccbcrd


  INTEGER(2), ALLOCATABLE :: itmp(:,:,:)      ! Temporary array
  REAL,       ALLOCATABLE :: hmax(:), hmin(:) ! Temporary array

  INTEGER :: clipxy, clipz

  INTEGER, SAVE :: bfid
  INTEGER, SAVE :: itime
  INTEGER       :: istop

  INTEGER           :: i, j, k
  INTEGER           :: nxlg, nylg, nzlg
  REAL, ALLOCATABLE :: var3d(:,:,:)

  REAL, ALLOCATABLE :: var3du(:,:,:)          ! for NetCDF I/O only
  REAL, ALLOCATABLE :: var3dv(:,:,:)
  REAL, ALLOCATABLE :: var3dw(:,:,:)

!
! Add the following only for HDF format, which is for
! split_hdf to work properly.
!
! fmtver??: to label each data a version.
!
  CHARACTER(LEN=40) :: fmtverin

  CHARACTER(LEN=4)  :: upcase

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3
  nzlg = nz

  IF (myproc == 0) THEN
    IF (exbcfmt == 3) THEN
      ALLOCATE (itmp(nxlg,nylg,nzlg),stat=istat)
      CALL check_alloc_status(istat, "readsplitexbc:itmp")

      ALLOCATE (hmax(nz),stat=istat)
      CALL check_alloc_status(istat, "readsplitexbc:hmax")

      ALLOCATE (hmin(nz),stat=istat)
      CALL check_alloc_status(istat, "readsplitexbc:hmin")

    END IF

    ALLOCATE(var3d(nxlg, nylg, nzlg), stat=istat)
    CALL check_alloc_status(istat, "readsplitexbc:var3d")


    WRITE (6,'(1x,3a)') 'READSPLITEXBC: reading in external boundary data ',    &
              'from file - ',filename(1:lfname)
!
!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------
!
    IF (exbcfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

      CALL getunit( bfid )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(filename(1:lfname), '-F f77 -N ieee', ierr)

      OPEN(bfid,FILE=filename(1:lfname),STATUS='old',                   &
                FORM='unformatted',IOSTAT=istat)

      IF ( istat /= 0 ) THEN
        ierr = 1
        GO TO 900
      END IF

      READ(bfid,ERR=999) nxin,nyin,nzin,dxin,dyin,dzin,                 &
                         ctrlatin,ctrlonin

      READ(bfid,ERR=999) ctime

      old_v = 0 ! In case that the EXBC files are of an earlier
                ! version that does not contain water and ice variables,
                ! set old_v to 1. Otherwise, set it to 0.

      IF( old_v == 1 ) THEN

        READ(bfid,ERR=999) ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd

        qcbcrd=0
        qrbcrd=0
        qibcrd=0
        qsbcrd=0
        qhbcrd=0

      ELSE

        READ(bfid,ERR=999) ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,             &
                        qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,             &
                        qhbcrd,qgbcrd,ncbcrd,nrbcrd,nibcrd,             &
                        nsbcrd,ngbcrd,nhbcrd,zrbcrd,zibcrd,             &
                        zsbcrd,zgbcrd,zhbcrd,ccbcrd,idummy,             &
                        idummy,idummy,idummy,idummy,idummy,             &
                        idummy,idummy,idummy,idummy,idummy,             &
                        idummy,idummy,idummy,idummy,nscalarin

      END IF

      qscalarbcrd(:) = 0
      IF (P_QC > 0) qscalarbcrd(P_QC) = qcbcrd
      IF (P_QR > 0) qscalarbcrd(P_QR) = qrbcrd
      IF (P_QI > 0) qscalarbcrd(P_QI) = qibcrd
      IF (P_QS > 0) qscalarbcrd(P_QS) = qsbcrd
      IF (P_QH > 0) qscalarbcrd(P_QH) = qhbcrd
      IF (P_QG > 0) qscalarbcrd(P_QG) = qgbcrd

      IF (P_NC > 0) qscalarbcrd(P_NC) = ncbcrd
      IF (P_NR > 0) qscalarbcrd(P_NR) = nrbcrd
      IF (P_NI > 0) qscalarbcrd(P_NI) = nibcrd
      IF (P_NS > 0) qscalarbcrd(P_NS) = nsbcrd
      IF (P_NH > 0) qscalarbcrd(P_NH) = nhbcrd
      IF (P_NG > 0) qscalarbcrd(P_NG) = ngbcrd

      IF (P_ZR > 0) qscalarbcrd(P_ZR) = zrbcrd
      IF (P_ZI > 0) qscalarbcrd(P_ZI) = zibcrd
      IF (P_ZS > 0) qscalarbcrd(P_ZS) = zsbcrd
      IF (P_ZH > 0) qscalarbcrd(P_ZH) = zhbcrd
      IF (P_ZG > 0) qscalarbcrd(P_ZG) = zgbcrd

      IF (P_CC > 0) qscalarbcrd(P_CC) = ccbcrd
    ELSE IF (exbcfmt == 3) THEN           !  HDF4 format.

      CALL hdfopen(trim(filename(1:lfname)), 1, bfid)
      IF (bfid < 0) THEN
        WRITE (6,*) "READSPLITEXBC: ERROR opening ",                    &
                    trim(filename(1:lfname))," for reading."
        ierr = 1
        GO TO 900
      END IF

      CALL hdfrdc(bfid,40,"fmtver",fmtverin,istat)
                    ! Should check fmtverin here

      CALL hdfrdc(bfid,15,"ctime",ctime,istat)

      CALL hdfrdi(bfid,"nx",nxin,istat)
      CALL hdfrdi(bfid,"ny",nyin,istat)
      CALL hdfrdi(bfid,"nz",nzin,istat)
      CALL hdfrdr(bfid,"dx",dxin,istat)
      CALL hdfrdr(bfid,"dy",dyin,istat)
      CALL hdfrdr(bfid,"dz",dzin,istat)
      CALL hdfrdr(bfid,"dzmin",dzminin,istat)
      CALL hdfrdi(bfid,"strhopt",strhoptin,istat)
      CALL hdfrdr(bfid,"zrefsfc",zrefsfcin,istat)
      CALL hdfrdr(bfid,"dlayer1",dlayer1in,istat)
      CALL hdfrdr(bfid,"dlayer2",dlayer2in,istat)
      CALL hdfrdr(bfid,"zflat",zflatin,istat)
      CALL hdfrdr(bfid,"strhtune",strhtunein,istat)
      CALL hdfrdi(bfid,"mapproj",maprojin,istat)
      CALL hdfrdr(bfid,"trulat1",trulat1in,istat)
      CALL hdfrdr(bfid,"trulat2",trulat2in,istat)
      CALL hdfrdr(bfid,"trulon",trulonin,istat)
      CALL hdfrdr(bfid,"sclfct",sclfctin,istat)
      CALL hdfrdr(bfid,"ctrlat",ctrlatin,istat)
      CALL hdfrdr(bfid,"ctrlon",ctrlonin,istat)

      CALL hdfrdi(bfid,"ubcflg",ubcrd,istat)
      CALL hdfrdi(bfid,"vbcflg",vbcrd,istat)
      CALL hdfrdi(bfid,"wbcflg",wbcrd,istat)
      CALL hdfrdi(bfid,"ptbcflg",ptbcrd,istat)
      CALL hdfrdi(bfid,"prbcflg",prbcrd,istat)
      CALL hdfrdi(bfid,"qvbcflg",qvbcrd,istat)

!     CALL hdfrdi(bfid,"nscalar",nscalarin,istat) ! not used in HDF format

      nscalarin = 0
      DO nq = 1,nscalar
        CALL hdfrdi(bfid,TRIM(qnames(nq))//'bcflg',qscalarbcrd(nq),istat)
        IF(qscalarbcrd(nq) > 0) THEN
          nscalarin = nscalarin + 1
        END IF
      END DO

      CALL hdfrdi(bfid, 'clipxy', clipxy,istat)
      IF (istat == 0 .AND. clipxy < ngbrz) THEN
        WRITE (6,*) "READSPLITEXBC: ERROR, clipxy (ngbrz) in exbc file too small"
        ierr = 1
        GO TO 900
      END IF
      CALL hdfrdi(bfid, 'clipz', clipz,istat)
      IF (istat == 0 .AND. clipz > rayklow) THEN
        WRITE (6,*) "READSPLITEXBC: ERROR, clipz (rayklow) in exbc file too large"
        ierr = 1
        GO TO 900
      END IF

    ELSE IF (exbcfmt == 7 .OR. exbcfmt == 8) THEN           !  NetCDF format

      IF (exbcfmt == 7) THEN
        itime = 1
      ELSE
        itime = itime + 1
        istop = NINT((tstop-tstart)/thisdmp) + 1
      END IF

      IF (itime == 1) CALL netopen(TRIM(filename(1:lfname)), 'R', bfid)

      CALL net_get_exbc(bfid,nxin,nyin,nzin,itime,dxin,dyin,dzin,       &
                 dzminin,strhoptin,zrefsfcin,dlayer1in,dlayer2in,       &
                 zflatin,strhtunein,maprojin,sclfctin,trulat1in,        &
                 trulat2in,trulonin,ctrlatin,ctrlonin,                  &
                 ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,qvbcrd,                &
                 nscalarin,qscalarbcrd,ctime,istat)

    ELSE

      ! alternate exbc format ...
      WRITE(6,*) 'The supported exbc data format are ',                &
               'binary (exbcfmt=1) and HDF4 no compressed (exbcfmt = 3).'
      CALL arpsstop('Exbc data format is not supported.',1)

    END IF     ! exbcfmt

    nxin = (nxin-3)/nproc_x + 3
    nyin = (nyin-3)/nproc_y + 3

  END IF     ! myproc == 0

  CALL mpupdatei(nxin, 1)
  CALL mpupdatei(nyin, 1)
  CALL mpupdatei(nzin, 1)
  CALL mpupdater(dxin, 1)
  CALL mpupdater(dyin, 1)
  CALL mpupdater(dzin, 1)
  CALL mpupdater(ctrlatin, 1)
  CALL mpupdater(ctrlonin, 1)
  CALL mpupdatec(ctime, 15)

  CALL mpupdatei(ubcrd, 1)
  CALL mpupdatei(vbcrd, 1)
  CALL mpupdatei(wbcrd, 1)
  CALL mpupdatei(ptbcrd, 1)
  CALL mpupdatei(prbcrd, 1)
  CALL mpupdatei(qvbcrd, 1)
  CALL mpupdatei(qscalarbcrd, nscalar)
  CALL mpupdatei(nscalarin, 1) ! DTD: It seems this is also needed

  IF(exbcfmt == 3 .OR. exbcfmt == 7 .OR. exbcfmt == 8) THEN
    CALL mpupdater(dzminin,   1)
    CALL mpupdatei(strhoptin, 1)
    CALL mpupdater(zrefsfcin, 1)
    CALL mpupdater(dlayer1in, 1)
    CALL mpupdater(dlayer2in, 1)
    CALL mpupdater(zflatin,   1)
    CALL mpupdater(strhtunein,1)
    CALL mpupdatei(maprojin,  1)
    CALL mpupdater(trulat1in, 1)
    CALL mpupdater(trulat2in, 1)
    CALL mpupdater(trulonin,  1)
    CALL mpupdater(sclfctin,  1)
  END IF

  IF (exbcfmt == 3) THEN     ! HDF 4 format specific
    CALL mpupdatec(fmtverin, 40)
    CALL mpupdatei(clipxy, 1)
    CALL mpupdatei(clipz,  1)
  END IF

!-----------------------------------------------------------------------
!
!  Check the data file for consistent grid parameters.
!
!-----------------------------------------------------------------------

  IF (exbcfmt == 1) THEN

    IF (myproc == 0)  WRITE (6,*)                                       &
      "READSPLITEXBC: WARNING, not checking all map projection parameters"

    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
            dx,dy,ctrlat,ctrlon,                                        &
            mapproj,trulat1,trulat2,trulon,sclfct,                      &
            dxin,dyin,ctrlatin,ctrlonin,                                &
            mapproj,trulat1,trulat2,trulon,sclfct,ireturn)

  ELSE

    CALL checkgrid3d(nx,ny,nz,nxin,nyin,nzin,                           &
            dx,dy,dz,dzmin,ctrlat,ctrlon,                               &
            strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,             &
            mapproj,trulat1,trulat2,trulon,sclfct,                      &
            dxin,dyin,dzin,dzminin,ctrlatin,ctrlonin,                   &
            strhoptin,zrefsfcin,dlayer1in,dlayer2in,zflatin,strhtunein, &
            maprojin,trulat1in,trulat2in,trulonin,sclfctin,ireturn)

  END IF

  IF (ireturn /= 0) THEN
    WRITE (6,*) "READSPLITEXBC: ERROR, grid parameter mismatch"
    ierr = 1
    GO TO 999
  END IF

!-----------------------------------------------------------------------
!
!  Read in the external boundary file data
!
!-----------------------------------------------------------------------

  IF (exbcfmt == 1) THEN

    IF (myproc == 0) READ(bfid,ERR=999) var3d
    CALL mpisplit3d(var3d, nx,ny,nz,u)
    IF (myproc == 0) READ(bfid,ERR=999) var3d
    CALL mpisplit3d(var3d, nx,ny,nz,v)
    IF (myproc == 0) READ(bfid,ERR=999) var3d
    CALL mpisplit3d(var3d, nx,ny,nz,w)
    IF (myproc == 0) READ(bfid,ERR=999) var3d
    CALL mpisplit3d(var3d, nx,ny,nz,pt)
    IF (myproc == 0) READ(bfid,ERR=999) var3d
    CALL mpisplit3d(var3d, nx,ny,nz,pr)

    IF(qvbcrd == 1) THEN
      IF (myproc == 0) READ(bfid,ERR=999) var3d
      CALL mpisplit3d(var3d, nx,ny,nz,qv)
    END IF

    IF (nscalarin > 0 ) THEN

      DO nqin = 1,nscalarin
        IF (myproc == 0) READ(bfid,ERR=999) var3d
        DO nq = 1,nscalar
          IF (nqin == qscalarbcrd(nq)) THEN
            CALL mpisplit3d(var3d, nx,ny,nz,qscalar(:,:,:,nq))
            EXIT
          END IF
        END DO
      END DO
      IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
        IF(myproc == 0) WRITE(6,*) 'Diagnosing additional moments at external boundary for multi-moment scheme'
        CALL init_MM_exbc(nx,ny,nz,pt,pr,qscalar)
        ! Also update qscalarbcrd flags to reflect the newly diagnosed additional moments
        ! even though they were not actually read in from the external boundary file
        DO nqin = 1,nscalarin
          DO nq = 1,nscalar
            IF(nqin == qscalarbcrd(nq)) THEN
              qscalarbcrd(nq) = nq
              qscalarbcrd(nq+6) = nq+6
              IF(mphyopt == 11) THEN
                qscalarbcrd(nq+11) = nq+11
              END IF
            END IF
          END DO
        END DO
        CALL mpupdatei(qscalarbcrd, nscalar)
      END IF
    ELSE                        ! earlier version?
      IF(qcbcrd == 1) THEN
        IF (myproc == 0) READ(bfid,ERR=999) var3d
        IF (P_QC > 0) CALL mpisplit3d(var3d, nx,ny,nz,qscalar(:,:,:,P_QC))
      END IF
      IF(qrbcrd == 1) THEN
        IF (myproc == 0) READ(bfid,ERR=999) var3d
        IF (P_QR > 0) CALL mpisplit3d(var3d, nx,ny,nz,qscalar(:,:,:,P_QR))
      END IF
      IF(qibcrd == 1) THEN
        IF (myproc == 0) READ(bfid,ERR=999) var3d
        IF (P_QI > 0) CALL mpisplit3d(var3d, nx,ny,nz,qscalar(:,:,:,P_QI))
      END IF
      IF(qsbcrd == 1) THEN
        IF (myproc == 0) READ(bfid,ERR=999) var3d
        IF (P_QS > 0) CALL mpisplit3d(var3d, nx,ny,nz,qscalar(:,:,:,P_QS))
      END IF
      IF(qhbcrd == 1) THEN
        IF (myproc == 0) READ(bfid,ERR=999) var3d
        IF (P_QH > 0) CALL mpisplit3d(var3d, nx,ny,nz,qscalar(:,:,:,P_QH))
      END IF
    END IF

  ELSE IF (exbcfmt == 3) THEN      ! HDF 4 format

    IF (ubcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL hdfrd3d(bfid,"u",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
      CALL mpisplit3d(var3d,nx, ny, nz,u)
    END IF
    IF (vbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL hdfrd3d(bfid,"v",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
      CALL mpisplit3d(var3d,nx, ny, nz,v)
    END IF
    IF (wbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL hdfrd3d(bfid,"w",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
      CALL mpisplit3d(var3d,nx, ny, nz,w)
    END IF
    IF (ptbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL hdfrd3d(bfid,"pt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
      CALL mpisplit3d(var3d,nx, ny, nz,pt)
    END IF
    IF (prbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL hdfrd3d(bfid,"p",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
      CALL mpisplit3d(var3d,nx, ny, nz,pr)
    END IF
    IF (qvbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL hdfrd3d(bfid,"qv",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat > 1) GO TO 999
      END IF
      CALL mpisplit3d(var3d,nx, ny, nz,qv)
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcrd(nq) > 0) THEN
        IF (myproc == 0) THEN
          CALL hdfrd3d(bfid,qnames(nq),nxlg,nylg,nzlg,var3d,            &
                       istat,itmp,hmax,hmin)
          IF (istat > 1) GO TO 999
        END IF
        CALL mpisplit3d(var3d,nx, ny, nz, qscalar(:,:,:,nq))
      END IF
    END DO
    IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
      IF(myproc == 0) THEN
        WRITE(6,*) 'Diagnosing additional moments at external boundary for multi-moment scheme'
      END IF
      CALL init_MM_exbc(nx,ny,nz,pt,pr,qscalar)
      ! Also update qscalarbcrd flags to reflect the newly diagnosed additional moments
      ! even though they were not actually read in from the external boundary file
      DO nqin = 1,nscalarin
        DO nq = 1,nscalar
          IF(nqin == qscalarbcrd(nq)) THEN
            qscalarbcrd(nq) = nq
            qscalarbcrd(nq+6) = nq+6
            IF(mphyopt == 11 .and. nq /= P_QC) THEN
              qscalarbcrd(nq+11) = nq+11
            END IF
            EXIT
          END IF
        END DO
      END DO
      CALL mpupdatei(qscalarbcrd, nscalar)
    END IF
  ELSE IF (exbcfmt == 7 .OR. exbcfmt == 8) THEN      ! NetCDF format

    ALLOCATE(var3du(nxlg,  nylg-1,nzlg-1), STAT = istat)
    ALLOCATE(var3dv(nxlg-1,nylg,  nzlg-1), STAT = istat)
    ALLOCATE(var3dw(nxlg-1,nylg-1,nzlg),   STAT = istat)

    IF (ubcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(bfid,0,itime,"U",nxlg,nylg-1,nz-1,var3du)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg
              var3d(i,j,k) = var3du(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d, nx, ny, nz, u)
    END IF

    IF (vbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(bfid,0,itime,"V",nxlg-1,nylg,nzlg-1,var3dv)
        DO k = 1,nz-1
          DO j = 1,nylg
            DO i = 1,nxlg-1
              var3d(i,j,k) = var3dv(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d, nx, ny, nz, v)
    END IF

    IF (wbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(bfid,0,itime,"W",nxlg-1,nylg-1,nzlg,var3dw)
        DO k = 1,nz
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = var3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz)
      END IF
      CALL mpisplit3d(var3d, nx, ny, nz, w)
    END IF

    IF (ptbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(bfid,0,itime,"PT",nxlg-1,nylg-1,nz-1,var3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = var3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d, nx, ny, nz, pt)
    END IF

    IF (prbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(bfid,0,itime,"P",nxlg-1,nylg-1,nz-1,var3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = var3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d, nx, ny, nz, pr)
    END IF

    IF (qvbcrd == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(bfid,0,itime,"QV",nxlg-1,nylg-1,nz-1,var3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = var3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d, nx, ny, nz, qv)
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcrd(nq) > 0) THEN
        IF (myproc == 0) THEN
          CALL netread3d(bfid,0,itime,upcase(qnames(nq)),nxlg-1,nylg-1,nz-1,var3dw)
          DO k = 1,nz-1
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var3d(i,j,k) = var3dw(i,j,k)
              END DO
            END DO
          END DO
          CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
        END IF
        CALL mpisplit3d(var3d, nx, ny, nz, qscalar(:,:,:,nq))
      END IF
    END DO
    IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
      IF(myproc == 0) THEN
        WRITE(6,*) 'Diagnosing additional moments at external boundary for multi-moment scheme'
      END IF
      CALL init_MM_exbc(nx,ny,nz,pt,pr,qscalar)
      ! Also update qscalarbcrd flags to reflect the newly diagnosed additional moments
      ! even though they were not actually read in from the external boundary file
      DO nqin = 1,nscalarin
        DO nq = 1,nscalar
          IF(nqin == qscalarbcrd(nq)) THEN
            qscalarbcrd(nq) = nq
            qscalarbcrd(nq+6) = nq+6
            IF(mphyopt == 11 .and. nq /= P_QC) THEN
              qscalarbcrd(nq+11) = nq+11
            END IF
            EXIT
          END IF
        END DO
      END DO
      CALL mpupdatei(qscalarbcrd, nscalar)
    END IF
  END IF

  IF(myproc == 0)THEN

    write(6,'(/1x,a/)') 'Max. and Min. of EXBC data variables:'

    IF (ubcrd == 1) THEN
    CALL a3dmax0lcl(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'umin = ', amin,',  umax =',amax
    END IF

    IF (vbcrd == 1) THEN
    CALL a3dmax0lcl(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'vmin = ', amin,',  vmax =',amax
    END IF

    IF (wbcrd == 1) THEN
    CALL a3dmax0lcl(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'wmin = ', amin,',  wmax =',amax
    END IF

    IF (ptbcrd == 1) THEN
    CALL a3dmax0lcl(pt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'ptmin= ', amin,',  ptmax=',amax
    END IF

    IF (prbcrd == 1) THEN
    CALL a3dmax0lcl(pr,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'pmin = ', amin,',  pmax =',amax
    END IF

    IF (qvbcrd == 1) THEN
    CALL a3dmax0lcl(qv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    write(6,'(1x,2(a,e13.6))') 'qvmin= ', amin,',  qvmax=',amax
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcrd(nq) > 0) THEN
      CALL a3dmax0lcl(qscalar(:,:,:,nq),1,nx,1,nx-1,1,ny,1,ny-1,        &
                      1,nz,1,nz-1,amax,amin)
      WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min= ', amin,       &
                          ',  '//TRIM(qnames(nq))//'max=',amax
      END IF
    END DO

  END IF

  ierr = 0

  GO TO 900

  999   CONTINUE

  WRITE (6,*) "READSPLITEXBC: ERROR reading data ",                     &
              "from file ",filename(1:lfname)," returning"
  ierr = 2

  900   CONTINUE

  IF (myproc == 0) THEN

    IF (exbcfmt == 1) THEN

      CLOSE (bfid)
      CALL retunit( bfid )

    ELSE IF (exbcfmt == 3) THEN

      CALL hdfclose(bfid,istat)

      DEALLOCATE (itmp,stat=istat)
      DEALLOCATE (hmax,stat=istat)
      DEALLOCATE (hmin,stat=istat)

    ELSE IF (exbcfmt == 7 .OR. exbcfmt == 8) THEN

      IF (exbcfmt == 7 .OR. itime >= istop) CALL netclose(bfid)

      DEALLOCATE (var3du)
      DEALLOCATE (var3dv)
      DEALLOCATE (var3dw)

    END IF

    DEALLOCATE(var3d, stat=istat)

  END IF  ! myproc == 0
  CALL mpupdatei(istat,1)
  ierr = istat

  CALL mpupdatei(ierr,1)

  RETURN
END SUBROUTINE readsplitexbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETBCFN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getbcfn( abstsec, exbcnam, tinite, tintve, filename, lfname, &
                    istat )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Get the external boundary data file name.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/1994
!
!  MODIFICATION HISTORY:
!
!  2000/03/24 (Gene Basett)
!  Added HDF4 format dumps.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    abstsec  Absolute seconds from 00:00:00, Jan. 1, 1980
!    exbcnam  A prefix of the external boundary condition file.
!
!    tinite   The boundary forecast initial time in yydddhhmm, In
!             general, the boundary files are named in yydddhhmm.
!    tintve   EXBC forecast time interval in seconds
!
!  OUTPUT:
!
!    filename  File name of EXBC boundary data set.
!    lfname    Length of the filename
!
!    istat     Status of finding the file.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations and COMMON blocks.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: abstsec

  CHARACTER (LEN=* ) :: filename
  INTEGER :: lfname

  CHARACTER (LEN=* ) :: exbcnam
  CHARACTER (LEN=* ) :: tinite
  INTEGER :: tintve

  INTEGER :: lenstr,istat
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=15) :: ctime
  CHARACTER (LEN=1)  :: chr

  INTEGER :: abstsec1
  INTEGER :: bcfcst, bcfcstop
  INTEGER :: iyr, imon, idy, ihr, imin, isec

  INTEGER, PARAMETER ::  maxtry = 10

  LOGICAL :: iexist

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ (tinite,'(i4,a,i2,a,i2,a,i2,a,i2,a,i2)')                         &
                         iyr,chr,imon,chr,idy,chr,ihr,chr,imin,chr,isec

  CALL ctim2abss( iyr,imon,idy,ihr,imin,isec, abstsec1 )

  bcfcstop = abststop + maxtry * tintve

  bcfcst = abstsec - MOD( abstsec-abstsec1, tintve )
  IF ( abstsec < abstsec1 ) bcfcst = bcfcst - tintve

  100   CONTINUE

  bcfcst = bcfcst + tintve

  IF ( bcfcst == abstsec .OR. bcfcst < abstsec1 ) GO TO 100

  CALL abss2ctim( bcfcst, iyr, imon, idy, ihr, imin, isec )

  WRITE (ctime,'(i4.4,2i2.2,a,3i2.2)') iyr,imon,idy,'.',ihr,imin,isec

  lenstr = 256
  CALL strlnth( exbcnam, lenstr )
  lfname = lenstr + 16

  filename(1:lfname) = exbcnam(1:lenstr)//'.'//ctime

  IF (mp_opt > 0 .AND. readsplit(FINDX_B) /= 1) THEN
    CALL gtsplitfn(exbcnam(1:lenstr)//'.'//ctime,1,1,loc_x,loc_y,1,1,   &
                   0,0,1,-1,filename,istat)
    lfname = LEN_TRIM(filename)

    INQUIRE(FILE=filename(1:lfname),EXIST=iexist)
  ELSE
    IF (myproc == 0) INQUIRE(FILE=filename(1:lfname),EXIST=iexist)
    CALL mpupdatel(iexist,1)
  END IF

  IF ( iexist ) THEN
    istat = 0
    IF (myproc == 0) WRITE (6,'(a,a)')                                  &
      'External boundary data file has been found: ', filename(1:lfname)
  ELSE IF ( bcfcst <= bcfcstop ) THEN
    IF (myproc == 0) WRITE (6,'(a,a,a)')                                &
      'External BC data file ', filename(1:lfname),                     &
      ' could not be found. Try another time.'
    GO TO 100
  ELSE
    IF (myproc == 0) WRITE (6,'(a,a)')                                  &
        'No external BC data file could not be found within a time range' &
        ,'Job will stop.'
    istat = 1
  END IF

  RETURN
END SUBROUTINE getbcfn
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRITEXBC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE writexbc(nx,ny,nz,filename,lfname,ctime,                     &
                    ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,               &
                    qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,qhbcdmp,    &
                    qgbcdmp,nqbcdmp,zqbcdmp,                            &
                    u,v,w,pt,pr,qv,qscalar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Output 6 primary fields for use as external boundary conditions.
!  In general these data come from another model.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  May, 1994
!
!  MODIFICATION HISTORY:
!
!  5/26/94 (Yuhe Liu)
!  Merged into the part of ARPS for external boundary conditions.
!
!  2000/03/24 (Gene Basett)
!  Added HDF4 format dumps.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction (east/west)
!  ny       Number of grid points in the y-direction (north/south)
!  nz       Number of grid points in the vertical
!
!  filename File name of EXBC boundary data set.
!  lfname   Length of the filename
!
!  ctime    Charater representation of the time of the EXBC data
!
!  dx       Expected x grid spacing
!  dy       Expected y grid spacing
!  dz       Expected z grid spacing
!
!  ctrlat   Expected center latitude
!  ctrlon   Expected center longitude
!
!  u        x component of velocity (m/s)
!  v        y component of velocity (m/s)
!  w        Vertical component of Cartesian velocity (m/s)
!  pt       Potential temperature (K)
!  pr       Pressure (Pascal)
!  qv       Water vapor specific humidity (kg/kg)
!  qc       Cloud water mixing ratio humidity (kg/kg)
!  qr       Rain water mixing ratio humidity (kg/kg)
!  qi       Cloud ice mixing ratio humidity (kg/kg)
!  qs       Snow mixing ratio humidity (kg/kg)
!  qh       Hail water mixing ratio humidity (kg/kg)
!
!  OUTPUT:
!
!  none
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  CHARACTER (LEN=*) :: filename
  CHARACTER (LEN=*) :: ctime
  INTEGER :: lfname
  INTEGER :: ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp
  INTEGER :: qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,qhbcdmp
  INTEGER :: qgbcdmp,nqbcdmp,zqbcdmp
  REAL    :: u(nx,ny,nz)       ! u-velocity (m/s)
  REAL    :: v(nx,ny,nz)       ! v-velocity (m/s)
  REAL    :: w(nx,ny,nz)       ! w-velocity (m/s)
  REAL    :: pt(nx,ny,nz)      ! Potential temperature (K)
  REAL    :: pr(nx,ny,nz)      ! Pressure (Pascal)
  REAL    :: qv(nx,ny,nz)      ! Specific humidity (kg/kg)

  REAL    :: qscalar(nx,ny,nz,nscalar)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istat, ierr, idummy
  INTEGER :: i,j,k
  INTEGER :: nq

  INTEGER :: qscalarbcdmp(nscalar)
  INTEGER :: ncbcdmp,nrbcdmp,nibcdmp,nsbcdmp,ngbcdmp,nhbcdmp,           &
                     zrbcdmp,zibcdmp,zsbcdmp,zgbcdmp,zhbcdmp,           &
             ccbcdmp
  INTEGER :: nscalarout

  INTEGER                 :: exbccompr,exbccomprtmp
  INTEGER(2), ALLOCATABLE :: itmp(:,:,:)       ! Temporary array
  REAL,       ALLOCATABLE :: hmax(:), hmin(:)  ! Temporary array
  REAL,       ALLOCATABLE :: ctmp(:,:,:)       ! Temporary array

  REAL,       ALLOCATABLE :: var3du(:,:,:)
  REAL,       ALLOCATABLE :: var3dv(:,:,:)
  REAL,       ALLOCATABLE :: var3dw(:,:,:)

  INTEGER, SAVE :: itime = 0
  INTEGER, SAVE :: bfid
  INTEGER       :: istop

!
! Add the following only for HDF formate which is for
! split_hdf to work properly.
!
! fmtver??: to label each data a version.
!
  CHARACTER (LEN=40) :: fmtver,fmtverhdf410,fmtverhdf500,fmtverhdf530

  PARAMETER (fmtverhdf410='004.10 HDF4 Coded Data')
  PARAMETER (fmtverhdf500='005.00 HDF4 Coded Data')
  PARAMETER (fmtverhdf530='005.30 HDF4 Coded Data')

  INTEGER :: npxout, npyout
  INTEGER :: nxout,  nyout
  INTEGER :: ipx,    jpy
  INTEGER :: ia,     ja

  INTEGER :: lenbase
  CHARACTER(LEN=256) :: outfilename

!-----------------------------------------------------------------------

  CHARACTER(LEN=4)  :: upcase
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (exbcdmp == 0) RETURN

  IF (exbcdmp == 3) exbccompr = exbchdfcompr

  IF ( (nproc_x_out > nproc_x .OR. nproc_y_out > nproc_y) .AND.         &
       (splitexbc > 0) ) THEN

    npxout = nproc_x_out/nproc_x
    npyout = nproc_y_out/nproc_y

    nxout = (nx-3)/npxout + 3
    nyout = (ny-3)/npyout + 3

  ELSE
    npxout = 1
    npyout = 1
    nxout  = nx
    nyout  = ny

  END IF

  !lenbase = INDEX(filename,'_',.TRUE.)
  !idummy  = INDEX(filename,'.',.TRUE.)
  !IF (lenbase < idummy) lenbase = LEN_TRIM(filename)

  IF (myproc == 0) WRITE (6,'(1x,/,2a,/)')                              &
    'WRITEXBC: Opening external boundary file ',TRIM(filename(1:lfname))

  qscalarbcdmp(:) = 0
  nscalarout = 0

  IF (qcbcdmp > 0 .AND. P_QC > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QC) = P_QC
    qcbcdmp = P_QC
  ELSE
    qcbcdmp = 0
  END IF

  IF (qrbcdmp > 0 .AND. P_QR > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QR) = P_QR
    qrbcdmp = P_QR
  ELSE
    qrbcdmp = 0
  END IF

  IF (qibcdmp > 0 .AND. P_QI > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QI) = P_QI
    qibcdmp = P_QI
  ELSE
    qibcdmp = 0
  END IF

  IF (qsbcdmp > 0 .AND. P_QS > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QS) = P_QS
    qsbcdmp = P_QS
  ELSE
    qsbcdmp = 0
  END IF

  IF (qhbcdmp > 0 .AND. P_QH > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QH) = P_QH
    qhbcdmp = P_QH
  ELSE
    qhbcdmp = 0
  END IF

  IF (qgbcdmp > 0 .AND. P_QG > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QG) = P_QG
    qgbcdmp = P_QG
  ELSE
    qgbcdmp = 0
  END IF

  ncbcdmp=0; nrbcdmp=0; nibcdmp=0; nsbcdmp=0; ngbcdmp=0; nhbcdmp=0;
             zrbcdmp=0; zibcdmp=0; zsbcdmp=0; zgbcdmp=0; zhbcdmp=0;
  IF (nqbcdmp > 0) THEN
    IF (P_NC > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NC) = P_NC
      ncbcdmp = P_NC
    END IF

    IF (P_NR > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NR) = P_NR
      nrbcdmp = P_NR
    END IF

    IF (P_NI > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NI) = P_NI
      nibcdmp = P_NI
    END IF

    IF (P_NS > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NS) = P_NS
      nsbcdmp = P_NS
    END IF

    IF (P_NG > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NG) = P_NG
      ngbcdmp = P_NG
    END IF

    IF (P_NH > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NH) = P_NH
      nhbcdmp = P_NH
    END IF
  END IF

  IF (zqbcdmp > 0) THEN
    IF (P_ZR > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZR) = P_ZR
      zrbcdmp = P_ZR
    END IF

    IF (P_ZI > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZI) = P_ZI
      zibcdmp = P_ZI
    END IF

    IF (P_ZS > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZS) = P_ZS
      zsbcdmp = P_ZS
    END IF

    IF (P_ZG > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZG) = P_ZG
      zgbcdmp = P_ZG
    END IF

    IF (P_ZH > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZH) = P_ZH
      zhbcdmp = P_ZH
    END IF
  END IF

  IF (P_CC > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_CC) = P_CC
    ccbcdmp = P_CC
  ELSE
    ccbcdmp = 0
  END IF

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (exbcdmp == 1) THEN

    CALL getunit( bfid )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(filename(1:lfname), '-F f77 -N ieee', ierr)

    CALL fnversn(filename, lfname)

    OPEN (bfid,FILE=trim(filename(1:lfname)),STATUS='unknown',          &
               FORM='unformatted')
!
!-----------------------------------------------------------------------
!
!  Write grid and time descriptors
!
!-----------------------------------------------------------------------
!
    WRITE (bfid) nx,ny,nz,dx,dy,dz,ctrlat,ctrlon
    WRITE (bfid) ctime
!
!-----------------------------------------------------------------------
!
!  Write integers which indicate whether each of the
!  variables are valid.  These must be properly set
!  by the calling routine.
!
!-----------------------------------------------------------------------
!
    idummy = 0

    WRITE (bfid) ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,                  &
                 qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,               &
                 qhbcdmp,qgbcdmp,ncbcdmp,nrbcdmp,nibcdmp,               &
                 nsbcdmp,ngbcdmp,nhbcdmp,zrbcdmp,zibcdmp,               &
                 zsbcdmp,zgbcdmp,zhbcdmp,ccbcdmp,idummy,                &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,nscalarout
!
!-----------------------------------------------------------------------
!
!  Write each variable in a separate record
!
!-----------------------------------------------------------------------
!
    IF( ubcdmp == 1) WRITE (bfid) u
    IF( vbcdmp == 1) WRITE (bfid) v
    IF( wbcdmp == 1) WRITE (bfid) w
    IF(ptbcdmp == 1) WRITE (bfid) pt
    IF(prbcdmp == 1) WRITE (bfid) pr
    IF(qvbcdmp == 1) WRITE (bfid) qv
    DO nq = 1,nscalar
      IF(qscalarbcdmp(nq) > 0) WRITE (bfid) qscalar(:,:,:,nq)
    END DO

    CLOSE (bfid)
    CALL retunit( bfid )

  ELSE IF (exbcdmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    ALLOCATE (ctmp(nxout,nyout,nz), STAT = istat)
    CALL check_alloc_status(istat, "WRITEXBC:ctmp")

    IF (exbccompr > 3) THEN
      ALLOCATE (itmp(nxout,nyout,nz),stat=istat)
      CALL check_alloc_status(istat, "WRITEXBC:itmp")

      ALLOCATE (hmax(nz),stat=istat)
      CALL check_alloc_status(istat, "WRITEXBC:hmax")

      ALLOCATE (hmin(nz),stat=istat)
      CALL check_alloc_status(istat, "WRITEXBC:hmin")
    END IF

    DO jpy = 1, npyout
      DO ipx = 1, npxout

        ia = (ipx-1)*(nxout-3)
        ja = (jpy-1)*(nyout-3)

        IF (splitexbc > 0 .OR. mp_opt > 0 ) THEN
          CALL gtsplitfn(filename(1:lfname),npxout,npyout,              &
                         loc_x,loc_y,ipx,jpy,                           &
                         0,0,0,lvldbg,outfilename,istat)
        ELSE
          outfilename = filename
          CALL fnversn(outfilename,lfname)
        END IF

        CALL hdfopen(TRIM(outfilename), 2, bfid)
        IF (bfid < 0) THEN
          WRITE (6,*) "WRITEXBC: ERROR creating HDF4 file: ",           &
                      TRIM(outfilename)
          CALL arpsstop('Error on creating HDF4 file.',1)
        END IF

        fmtver = fmtverhdf530  !  for the time being, only version 5.30

        CALL hdfwrtc(bfid, 40, 'fmtver', fmtver, istat)

        CALL hdfwrtc(bfid, 15, 'ctime', ctime, istat)

        CALL hdfwrti(bfid, 'nx', nxout, istat)
        CALL hdfwrti(bfid, 'ny', nyout, istat)
        CALL hdfwrti(bfid, 'nz', nz, istat)
        CALL hdfwrtr(bfid, 'dx', dx, istat)
        CALL hdfwrtr(bfid, 'dy', dy, istat)
        CALL hdfwrtr(bfid, 'dz', dz, istat)
        CALL hdfwrtr(bfid, 'dzmin',   dzmin, istat)
        CALL hdfwrti(bfid, 'strhopt', strhopt, istat)
        CALL hdfwrtr(bfid, 'zrefsfc', zrefsfc, istat)
        CALL hdfwrtr(bfid, 'dlayer1', dlayer1, istat)
        CALL hdfwrtr(bfid, 'dlayer2', dlayer2, istat)
        CALL hdfwrtr(bfid, 'zflat',   zflat, istat)
        CALL hdfwrtr(bfid, 'strhtune', strhtune, istat)
        CALL hdfwrti(bfid, 'mapproj', mapproj, istat)
        CALL hdfwrtr(bfid, 'trulat1', trulat1, istat)
        CALL hdfwrtr(bfid, 'trulat2', trulat2, istat)
        CALL hdfwrtr(bfid, 'trulon', trulon, istat)
        CALL hdfwrtr(bfid, 'sclfct', sclfct, istat)
        CALL hdfwrtr(bfid, 'ctrlat', ctrlat, istat)
        CALL hdfwrtr(bfid, 'ctrlon', ctrlon, istat)

        CALL hdfwrti(bfid,"ubcflg", ubcdmp,istat)
        CALL hdfwrti(bfid,"vbcflg", vbcdmp,istat)
        CALL hdfwrti(bfid,"wbcflg", wbcdmp,istat)
        CALL hdfwrti(bfid,"ptbcflg",ptbcdmp,istat)
        CALL hdfwrti(bfid,"prbcflg",prbcdmp,istat)
        CALL hdfwrti(bfid,"qvbcflg",qvbcdmp,istat)

        CALL hdfwrti(bfid,"nscalar",nscalarout,istat)
        CALL hdfwrti(bfid,"qcbcflg",qcbcdmp,istat)
        CALL hdfwrti(bfid,"qrbcflg",qrbcdmp,istat)
        CALL hdfwrti(bfid,"qibcflg",qibcdmp,istat)
        CALL hdfwrti(bfid,"qsbcflg",qsbcdmp,istat)
        CALL hdfwrti(bfid,"qhbcflg",qhbcdmp,istat)
        CALL hdfwrti(bfid,"qgbcflg",qgbcdmp,istat)

        CALL hdfwrti(bfid,"ncbcflg",ncbcdmp,istat)
        CALL hdfwrti(bfid,"nrbcflg",nrbcdmp,istat)
        CALL hdfwrti(bfid,"nibcflg",nibcdmp,istat)
        CALL hdfwrti(bfid,"nsbcflg",nsbcdmp,istat)
        CALL hdfwrti(bfid,"ngbcflg",ngbcdmp,istat)
        CALL hdfwrti(bfid,"nhbcflg",nhbcdmp,istat)

        CALL hdfwrti(bfid,"zrbcflg",zrbcdmp,istat)
        CALL hdfwrti(bfid,"zibcflg",zibcdmp,istat)
        CALL hdfwrti(bfid,"zsbcflg",zsbcdmp,istat)
        CALL hdfwrti(bfid,"zgbcflg",zgbcdmp,istat)
        CALL hdfwrti(bfid,"zhbcflg",zhbcdmp,istat)

        CALL hdfwrti(bfid,"ccbcflg",ccbcdmp,istat)

        IF (exbccompr > 4) THEN
          CALL hdfwrti(bfid, 'clipxy', ngbrz,   istat)
          CALL hdfwrti(bfid, 'clipz',  rayklow, istat)
        END IF

        IF( ubcdmp == 1) THEN
          DO k=1,nz
            DO j=1,nyout
              DO i=1,nxout
                ctmp(i,j,k) = u(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          IF (exbccompr > 4) THEN
            ! reset values inside region not used by LBC forcing
            DO k=1,rayklow-1
              DO j=1+ngbrz,nyout-1-ngbrz
                DO i=2+ngbrz,nxout-2-ngbrz
                  ctmp(i,j,k) = ctmp(2,2,k)   ! just need some value within max/min
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccompr,           &
                        'u','u-velocity','m/s',itmp,hmax,hmin)
        END IF

        IF( vbcdmp == 1) THEN

          DO k=1,nz
            DO j=1,nyout
              DO i=1,nxout
                ctmp(i,j,k) = v(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          IF (exbccompr > 4) THEN
            ! reset values inside region not used by LBC forcing
            DO k=1,rayklow-1
              DO j=2+ngbrz,nyout-2-ngbrz
                DO i=1+ngbrz,nxout-1-ngbrz
                  ctmp(i,j,k) = ctmp(2,2,k)   ! just need some value with max/min
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccompr,           &
                        'v','v-velocity','m/s',itmp,hmax,hmin)
        END IF

        IF( wbcdmp == 1) THEN
          DO k=1,nz
            DO j=1,nyout
              DO i=1,nxout
                ctmp(i,j,k) = w(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          IF (exbccompr > 4) THEN
            ! reset values inside region not used by LBC forcing
            DO k=1,rayklow-2
              DO j=1+ngbrz,nyout-1-ngbrz
                DO i=1+ngbrz,nxout-1-ngbrz
                  ctmp(i,j,k) = ctmp(2,2,k)   ! just need some value with max/min
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccompr,           &
                        'w','w-velocity','m/s',itmp,hmax,hmin)
        END IF

        IF(ptbcdmp == 1) THEN

          DO k=1,nz
            DO j=1,nyout
              DO i=1,nxout
                ctmp(i,j,k) = pt(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          IF (exbccompr > 4) THEN
            ! reset values inside region not used by LBC forcing
            DO k=1,rayklow-1
              DO j=1+ngbrz,nyout-1-ngbrz
                DO i=1+ngbrz,nxout-1-ngbrz
                  ctmp(i,j,k) = ctmp(2,2,k)   ! just need some value with max/min
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccompr,           &
                        'pt','Potential temperature','K',itmp,hmax,hmin)
        END IF

        IF(prbcdmp == 1) THEN

          DO k=1,nz
            DO j=1,nyout
              DO i=1,nxout
                ctmp(i,j,k) = pr(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          IF (exbccompr > 4) THEN
            ! reset values inside region not used by LBC forcing
            DO k=1,rayklow-1
              DO j=1+ngbrz,nyout-1-ngbrz
                DO i=1+ngbrz,nxout-1-ngbrz
                  ctmp(i,j,k) = ctmp(2,2,k)   ! just need some value with max/min
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccompr,           &
                        'p','Pressure','Pascal',itmp,hmax,hmin)

        END IF

        IF(qvbcdmp == 1) THEN

          DO k=1,nz
            DO j=1,nyout
              DO i=1,nxout
                ctmp(i,j,k) = qv(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          IF (exbccompr > 4) THEN
            ! reset values inside region not used by LBC forcing
            DO k=1,rayklow-1
              DO j=1+ngbrz,nyout-1-ngbrz
                DO i=1+ngbrz,nxout-1-ngbrz
                  ctmp(i,j,k) = 0.0
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccompr,           &
                        'qv','Water vapor specific humidity','kg/kg',   &
                        itmp,hmax,hmin)

        END IF

        DO nq = 1, nscalar

          IF ( qscalarbcdmp(nq) > 0 ) THEN

            exbccomprtmp = exbccompr

            ! DTD: Turn off bit-packing for Z array for mphyopt == 11
            IF (nq >= 13 .AND. exbccompr > 3) THEN
              exbccomprtmp = exbccompr - 3
            END IF

            DO k=1,nz
              DO j=1,nyout
                DO i=1,nxout
                  ctmp(i,j,k) = qscalar(i+ia,j+ja,k,nq)
                END DO
              END DO
            END DO

            IF (exbccompr > 4) THEN
              ! reset values inside region not used by LBC forcing
              DO k=1,rayklow-1
                DO j=1+ngbrz,nyout-1-ngbrz
                  DO i=1+ngbrz,nxout-1-ngbrz
                    ctmp(i,j,k) = 0.0
                  END DO
                END DO
              END DO
            END IF

            CALL hdfwrt3d(ctmp,nxout,nyout,nz,bfid,1,exbccomprtmp,      &
                          TRIM(qnames(nq)),TRIM(qdescp(nq)),'kg/kg',    &
                          itmp,hmax,hmin)
          END IF

        END DO

        CALL hdfclose(bfid,istat)

      END DO
    END DO

    IF (exbccompr > 3) THEN
      DEALLOCATE (itmp,stat=istat)
      DEALLOCATE (hmax,stat=istat)
      DEALLOCATE (hmin,stat=istat)
    END IF

    DEALLOCATE (ctmp,stat=istat)

  ELSE IF (exbcdmp == 7 .OR. exbcdmp == 8) THEN
!
!-----------------------------------------------------------------------
!
!  Write out in NetCDF format.
!
!-----------------------------------------------------------------------

    IF (exbcdmp == 7) THEN
      itime = 1
    ELSE
      itime = itime + 1
      istop = NINT((tstop-tstart)/thisdmp) + 1
    END IF

    ALLOCATE(var3du(nxout,  nyout-1,nz-1), STAT = istat)
    ALLOCATE(var3dv(nxout-1,nyout,  nz-1), STAT = istat)
    ALLOCATE(var3dw(nxout-1,nyout-1,nz  ), STAT = istat)

    DO jpy = 1, npyout
      DO ipx = 1, npxout

        ia = (ipx-1)*(nxout-3)
        ja = (jpy-1)*(nyout-3)

        IF (itime == 1) THEN

!-----------------------------------------------------------------------
!
!  Define ARPS boundary file dimension and variables
!
!-----------------------------------------------------------------------
          IF (splitexbc > 0 .OR. mp_opt > 0 ) THEN
            CALL gtsplitfn(filename(1:lfname),npxout,npyout,            &
                           loc_x,loc_y,ipx,jpy,                         &
                           0,0,0,lvldbg,outfilename,istat)
          ELSE
            outfilename = filename
            CALL fnversn(outfilename, lfname)
          END IF

          CALL netopen(TRIM(outfilename), 'C', bfid)

          CALL net_define_exbc(bfid,nxout,nyout,nz,itime,dx,dy,dz,       &
                   dzmin,strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune, &
                   mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,  &
                   ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,qvbcdmp,         &
                   nscalarout,qscalarbcdmp,ctime,istat)
        END IF

        IF( ubcdmp == 1) THEN
          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout
                var3du(i,j,k) = u(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,'U',var3du,nxout,nyout-1,nz-1)
        END IF

        IF( vbcdmp == 1) THEN
          DO k = 1,nz-1
            DO j = 1,nyout
              DO i = 1,nxout-1
                var3dv(i,j,k) = v(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,'V',var3dv,nxout-1,nyout,nz-1)
        END IF

        IF( wbcdmp == 1) THEN
          DO k = 1,nz
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = w(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,'W',var3dw,nxout-1,nyout-1,nz)
        END IF

        IF(ptbcdmp == 1) THEN
          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = pt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,'PT',var3dw,nxout-1,nyout-1,nz-1)
        END IF

        IF(prbcdmp == 1) THEN
          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = pr(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,'P',var3dw,nxout-1,nyout-1,nz-1)
        END IF

        IF(qvbcdmp == 1) THEN
          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = qv(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,'QV',var3dw,nxout-1,nyout-1,nz-1)
        END IF

        DO nq = 1, nscalar

          IF (qscalarbcdmp(nq) > 0) THEN
            DO k = 1,nz-1
              DO j = 1,nyout-1
                DO i = 1,nxout-1
                  var3dw(i,j,k) = qscalar(i+ia,j+ja,k,nq)
                END DO
              END DO
            END DO
            CALL netwrt3d(bfid,0,itime,upcase(TRIM(qnames(nq))),var3dw,nxout-1,nyout-1,nz-1)
          END IF

        END DO

        IF (exbcdmp == 7 .OR. itime >= istop)                           &
          CALL netclose(bfid)
      END DO
    END DO

    DEALLOCATE(var3du,var3dv,var3dw)

  ELSE

    ! alternate dump format ...
    WRITE(6,'(1x,3a)') 'The supported exbc data dump format are ',      &
               'binary (exbcdmp = 1), HDF4 (exbcdmp = 3) and ',         &
               'NetCDF format (exbcdmp = 7).'

    CALL arpsstop('EXBC data dump format is not supported.',1)

  END IF

  RETURN
END SUBROUTINE writexbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRITEJOINEXBC              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE writejoinexbc(nx,ny,nz,filename,lfname,ctime,                &
                        ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,           &
                        qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,qhbcdmp,&
                        qgbcdmp,nqbcdmp,zqbcdmp,u,v,w,pt,pr,qv,qscalar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Output 6 primary fields for use as external boundary conditions.
!  In general these data come from another model.
!
!  NOTE:
!    Rewrote from Kevin Thomas' joined subroutine writexbc. Maintain
!    two separate subroutines will keep each subroutine short and
!    managable.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (03/31/2005)
!  Based on subroutine writexbc by Keith Browster and the upgraded
!  version of writexbc by Kevin Thomas.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction (east/west)
!  ny       Number of grid points in the y-direction (north/south)
!  nz       Number of grid points in the vertical
!
!  filename File name of EXBC boundary data set.
!  lfname   Length of the filename
!
!  ctime    Charater representation of the time of the EXBC data
!
!  dx       Expected x grid spacing
!  dy       Expected y grid spacing
!  dz       Expected z grid spacing
!
!  ctrlat   Expected center latitude
!  ctrlon   Expected center longitude
!
!  u        x component of velocity (m/s)
!  v        y component of velocity (m/s)
!  w        Vertical component of Cartesian velocity (m/s)
!  pt       Potential temperature (K)
!  pr       Pressure (Pascal)
!  qv       Water vapor specific humidity (kg/kg)
!  qc       Cloud water mixing ratio humidity (kg/kg)
!  qr       Rain water mixing ratio humidity (kg/kg)
!  qi       Cloud ice mixing ratio humidity (kg/kg)
!  qs       Snow mixing ratio humidity (kg/kg)
!  qh       Hail water mixing ratio humidity (kg/kg)
!
!  OUTPUT:
!
!  none
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=256), INTENT(IN) :: filename
  CHARACTER(LEN=15 ), INTENT(IN) :: ctime
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: lfname
  INTEGER, INTENT(IN) :: ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp
  INTEGER, INTENT(INOUT) :: qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,qhbcdmp
  INTEGER, INTENT(INOUT) :: qgbcdmp
  INTEGER, INTENT(IN) :: nqbcdmp,zqbcdmp
  REAL,    INTENT(IN) :: u(nx,ny,nz)       ! u-velocity (m/s)
  REAL,    INTENT(IN) :: v(nx,ny,nz)       ! v-velocity (m/s)
  REAL,    INTENT(IN) :: w(nx,ny,nz)       ! w-velocity (m/s)
  REAL,    INTENT(IN) :: pt(nx,ny,nz)      ! Potential temperature (K)
  REAL,    INTENT(IN) :: pr(nx,ny,nz)      ! Pressure (Pascal)
  REAL,    INTENT(IN) :: qv(nx,ny,nz)      ! Specific humidity (kg/kg)

  REAL,    INTENT(IN) :: qscalar(nx,ny,nz,nscalar)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istat, ierr, idummy
  INTEGER :: i,j,k
  INTEGER :: nq
  INTEGER :: nxlg, nylg
  REAL    :: tmp_val

  INTEGER :: qscalarbcdmp(nscalar)
  INTEGER :: nscalarout
  INTEGER :: ncbcdmp,nrbcdmp,nibcdmp,nsbcdmp,ngbcdmp,nhbcdmp,           &
                     zrbcdmp,zibcdmp,zsbcdmp,zgbcdmp,zhbcdmp,           &
             ccbcdmp

  INTEGER                 :: exbccompr, exbccomprtmp
  INTEGER(2), ALLOCATABLE :: itmp(:,:,:)       ! Temporary array
  REAL,       ALLOCATABLE :: hmax(:), hmin(:)  ! Temporary array

  INTEGER, SAVE :: bfid
  INTEGER, SAVE :: itime
  INTEGER       :: istop

  REAL,       ALLOCATABLE :: out3d(:,:,:)
  REAL,       ALLOCATABLE :: var3du(:,:,:)
  REAL,       ALLOCATABLE :: var3dv(:,:,:)
  REAL,       ALLOCATABLE :: var3dw(:,:,:)
!
! Add the following only for HDF formate which is for
! split_hdf to work properly.
!
! fmtver??: to label each data a version.
!
  CHARACTER (LEN=40) :: fmtver,fmtverhdf410,fmtverhdf500,fmtverhdf530

  PARAMETER (fmtverhdf410='004.10 HDF4 Coded Data')
  PARAMETER (fmtverhdf500='005.00 HDF4 Coded Data')
  PARAMETER (fmtverhdf530='005.30 HDF4 Coded Data')

  CHARACTER(LEN=4) :: upcase
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (exbcdmp == 0) RETURN

  nxlg = (nx-3) * nproc_x + 3
  nylg = (ny-3) * nproc_y + 3

  ALLOCATE(out3d(nxlg,nylg,nz), STAT = istat)
  CALL check_alloc_status(istat,"WRITEXBC:out3d")

  IF (exbcdmp == 3) exbccompr = exbchdfcompr

  IF (exbcdmp == 3 .AND. exbccompr > 3) THEN

    ALLOCATE (itmp(nxlg,nylg,nz),stat=istat)
    CALL check_alloc_status(istat, "WRITEXBC:itmp")

    ALLOCATE (hmax(nz),stat=istat)
    CALL check_alloc_status(istat, "WRITEXBC:hmax")

    ALLOCATE (hmin(nz),stat=istat)
    CALL check_alloc_status(istat, "WRITEXBC:hmin")

  END IF

  CALL fnversn(filename, lfname)

  IF (myproc == 0) WRITE (6,*)                                          &
                     'WRITEJOINEXBC: Opening external boundary file ',  &
                     trim(filename(1:lfname))

  qscalarbcdmp(:) = 0
  nscalarout = 0

  IF (qcbcdmp > 0 .AND. P_QC > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QC) = P_QC
    qcbcdmp = P_QC
  ELSE
    qcbcdmp = 0
  END IF

  IF (qrbcdmp > 0 .AND. P_QR > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QR) = P_QR
    qrbcdmp = P_QR
  ELSE
    qrbcdmp = 0
  END IF

  IF (qibcdmp > 0 .AND. P_QI > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QI) = P_QI
    qibcdmp = P_QI
  ELSE
    qibcdmp = 0
  END IF

  IF (qsbcdmp > 0 .AND. P_QS > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QS) = P_QS
    qsbcdmp = P_QS
  ELSE
    qsbcdmp = 0
  END IF

  IF (qhbcdmp > 0 .AND. P_QH > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QH) = P_QH
    qhbcdmp = P_QH
  ELSE
    qhbcdmp = 0
  END IF

  IF (qgbcdmp > 0 .AND. P_QG > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_QG) = P_QG
    qgbcdmp = P_QG
  ELSE
    qgbcdmp = 0
  END IF

  ncbcdmp=0; nrbcdmp=0; nibcdmp=0; nsbcdmp=0; ngbcdmp=0; nhbcdmp=0;
             zrbcdmp=0; zibcdmp=0; zsbcdmp=0; zgbcdmp=0; zhbcdmp=0;
  IF (nqbcdmp > 0) THEN
    IF (P_NC > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NC) = P_NC
      ncbcdmp = P_NC
    END IF

    IF (P_NR > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NR) = P_NR
      nrbcdmp = P_NR
    END IF

    IF (P_NI > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NI) = P_NI
      nibcdmp = P_NI
    END IF

    IF (P_NS > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NS) = P_NS
      nsbcdmp = P_NS
    END IF

    IF (P_NG > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NG) = P_NG
      ngbcdmp = P_NG
    END IF

    IF (P_NH > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_NH) = P_NH
      nhbcdmp = P_NH
    END IF
  END IF

  IF (zqbcdmp > 0) THEN
    IF (P_ZR > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZR) = P_ZR
      zrbcdmp = P_ZR
    END IF

    IF (P_ZI > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZI) = P_ZI
      zibcdmp = P_ZI
    END IF

    IF (P_ZS > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZS) = P_ZS
      zsbcdmp = P_ZS
    END IF

    IF (P_ZG > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZG) = P_ZG
      zgbcdmp = P_ZG
    END IF

    IF (P_ZH > 0) THEN
      nscalarout = nscalarout + 1
      qscalarbcdmp(P_ZH) = P_ZH
      zhbcdmp = P_ZH
    END IF
  END IF

  IF ( P_CC > 0) THEN
    nscalarout = nscalarout + 1
    qscalarbcdmp(P_CC) = P_CC
    ccbcdmp = P_CC
  ELSE
    ccbcdmp = 0
  END IF

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (exbcdmp == 1) THEN

    IF (myproc == 0) THEN
      CALL getunit( bfid )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(filename(1:lfname), '-F f77 -N ieee', ierr)

      OPEN (bfid,FILE=trim(filename(1:lfname)),STATUS='unknown',        &
                 FORM='unformatted')
!
!-----------------------------------------------------------------------
!
!  Write grid and time descriptors
!
!-----------------------------------------------------------------------
!
      WRITE (bfid) nxlg,nylg,nz,dx,dy,dz,ctrlat,ctrlon
      WRITE (bfid) ctime
!
!-----------------------------------------------------------------------
!
!  Write integers which indicate whether each of the
!  variables are valid.  These must be properly set
!  by the calling routine.
!
!-----------------------------------------------------------------------
!
      idummy = 0

      WRITE (bfid) ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,                &
                   qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,             &
                   qhbcdmp,qgbcdmp,ncbcdmp,nrbcdmp,nibcdmp,             &
                   nsbcdmp,ngbcdmp,nhbcdmp,zrbcdmp,zibcdmp,             &
                   zsbcdmp,zgbcdmp,zhbcdmp,ccbcdmp,idummy,              &
                   idummy,idummy,idummy,idummy,idummy,                  &
                   idummy,idummy,idummy,idummy,idummy,                  &
                   idummy,idummy,idummy,idummy,nscalarout

    END IF
!
!-----------------------------------------------------------------------
!
!  Write each variable in a separate record
!
!-----------------------------------------------------------------------
!
    IF ( ubcdmp == 1 ) THEN
      CALL mpimerge3d(u,nx,ny,nz,out3d)
      IF (myproc == 0) WRITE (bfid) out3d
    END IF

    IF ( vbcdmp == 1 ) THEN
      CALL mpimerge3d(v,nx,ny,nz,out3d)
      IF (myproc == 0) WRITE (bfid) out3d
    END IF

    IF ( wbcdmp == 1 ) THEN
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF ( myproc == 0 ) WRITE (bfid) out3d
    END IF

    IF ( ptbcdmp == 1 ) THEN
      CALL mpimerge3d(pt,nx,ny,nz,out3d)
      IF ( myproc == 0 ) WRITE (bfid) out3d
    END IF

    IF ( prbcdmp == 1 ) THEN
      CALL mpimerge3d(pr,nx,ny,nz,out3d)
      IF ( myproc == 0 ) WRITE (bfid) out3d
    END IF

    IF ( qvbcdmp == 1 ) THEN
      CALL mpimerge3d(qv,nx,ny,nz,out3d)
      IF ( myproc == 0 ) WRITE (bfid) out3d
    END IF

    DO nq = 1, nscalar
      IF ( qscalarbcdmp(nq) > 0 ) THEN
        CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
        IF ( myproc == 0 ) WRITE (bfid) out3d
      END IF
    END DO

    IF(myproc == 0) THEN
      CLOSE (bfid)
      CALL retunit( bfid )
    END IF

  ELSE IF (exbcdmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    IF ( myproc == 0 ) THEN
      CALL hdfopen(trim(filename(1:lfname)), 2, bfid)
      IF (bfid < 0) THEN
        WRITE (6,*) "WRITEXBC: ERROR creating HDF4 file: ",             &
                  trim(filename(1:lfname))
        CALL arpsstop('Error on creating HDF4 file in writejoinexbc.',1)
      END IF

      fmtver = fmtverhdf530  !  for the time being, only version 5.30

      CALL hdfwrtc(bfid, 40, 'fmtver', fmtver, istat)
      CALL hdfwrtc(bfid, 15, 'ctime', ctime, istat)

      CALL hdfwrti(bfid, 'nx', nxlg, istat)
      CALL hdfwrti(bfid, 'ny', nylg, istat)
      CALL hdfwrti(bfid, 'nz', nz,   istat)
      CALL hdfwrtr(bfid, 'dx', dx,   istat)
      CALL hdfwrtr(bfid, 'dy', dy,   istat)
      CALL hdfwrtr(bfid, 'dz', dz,   istat)
      CALL hdfwrtr(bfid, 'dzmin',    dzmin,    istat)
      CALL hdfwrti(bfid, 'strhopt',  strhopt,  istat)
      CALL hdfwrtr(bfid, 'zrefsfc',  zrefsfc,  istat)
      CALL hdfwrtr(bfid, 'dlayer1',  dlayer1,  istat)
      CALL hdfwrtr(bfid, 'dlayer2',  dlayer2,  istat)
      CALL hdfwrtr(bfid, 'zflat',    zflat,    istat)
      CALL hdfwrtr(bfid, 'strhtune', strhtune, istat)
      CALL hdfwrti(bfid, 'mapproj', mapproj, istat)
      CALL hdfwrtr(bfid, 'trulat1', trulat1, istat)
      CALL hdfwrtr(bfid, 'trulat2', trulat2, istat)
      CALL hdfwrtr(bfid, 'trulon',  trulon,  istat)
      CALL hdfwrtr(bfid, 'sclfct',  sclfct,  istat)
      CALL hdfwrtr(bfid, 'ctrlat',  ctrlat,  istat)
      CALL hdfwrtr(bfid, 'ctrlon',  ctrlon,  istat)

      CALL hdfwrti(bfid,"ubcflg", ubcdmp, istat)
      CALL hdfwrti(bfid,"vbcflg", vbcdmp, istat)
      CALL hdfwrti(bfid,"wbcflg", wbcdmp, istat)
      CALL hdfwrti(bfid,"ptbcflg",ptbcdmp,istat)
      CALL hdfwrti(bfid,"prbcflg",prbcdmp,istat)
      CALL hdfwrti(bfid,"qvbcflg",qvbcdmp,istat)

      CALL hdfwrti(bfid,"nscalar",nscalarout,istat)

      CALL hdfwrti(bfid,"qcbcflg",qcbcdmp,istat)
      CALL hdfwrti(bfid,"qrbcflg",qrbcdmp,istat)
      CALL hdfwrti(bfid,"qibcflg",qibcdmp,istat)
      CALL hdfwrti(bfid,"qsbcflg",qsbcdmp,istat)
      CALL hdfwrti(bfid,"qhbcflg",qhbcdmp,istat)
      CALL hdfwrti(bfid,"qgbcflg",qgbcdmp,istat)

      CALL hdfwrti(bfid,"ncbcflg",ncbcdmp,istat)
      CALL hdfwrti(bfid,"nrbcflg",nrbcdmp,istat)
      CALL hdfwrti(bfid,"nibcflg",nibcdmp,istat)
      CALL hdfwrti(bfid,"nsbcflg",nsbcdmp,istat)
      CALL hdfwrti(bfid,"ngbcflg",ngbcdmp,istat)
      CALL hdfwrti(bfid,"nhbcflg",nhbcdmp,istat)

      CALL hdfwrti(bfid,"zrbcflg",zrbcdmp,istat)
      CALL hdfwrti(bfid,"zibcflg",zibcdmp,istat)
      CALL hdfwrti(bfid,"zsbcflg",zsbcdmp,istat)
      CALL hdfwrti(bfid,"zgbcflg",zgbcdmp,istat)
      CALL hdfwrti(bfid,"zhbcflg",zhbcdmp,istat)

      CALL hdfwrti(bfid,"ccbcflg",ccbcdmp,istat)

      IF (exbccompr > 4) THEN
        CALL hdfwrti(bfid, 'clipxy', ngbrz,   istat)
        CALL hdfwrti(bfid, 'clipz',  rayklow, istat)
      END IF

    END IF

    IF( ubcdmp == 1) THEN
      CALL mpimerge3d(u,nx,ny,nz,out3d)
      IF (myproc == 0 ) THEN
        IF (exbccompr > 4) THEN
          !
          ! reset values inside region not used by LBC forcing
          !
          DO k=1,rayklow-1
            tmp_val = out3d(2,2,k)
            DO j = 1+ngbrz, nylg-1-ngbrz
              DO i = 2+ngbrz, nxlg-2-ngbrz
                out3d(i,j,k) = tmp_val   ! just need some value with max/min
              END DO
            END DO
          END DO
        END IF
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccompr,               &
                      'u','u-velocity','m/s',itmp,hmax,hmin)
      END IF ! myproc == 0
    END IF

    IF( vbcdmp == 1) THEN
      CALL mpimerge3d(v,nx,ny,nz,out3d)
      IF (myproc == 0 ) THEN
        IF (exbccompr > 4) THEN
          !
          ! reset values inside region not used by LBC forcing
          !
          DO k=1,rayklow-1
            tmp_val = out3d(2,2,k)
            DO j= 2+ngbrz, nylg-2-ngbrz
              DO i= 1+ngbrz, nxlg-1-ngbrz
                out3d(i,j,k) = tmp_val   ! just need some value with max/min
              END DO
            END DO
          END DO
        END IF
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccompr,              &
                      'v','v-velocity','m/s',itmp,hmax,hmin)
      END IF
    END IF

    IF( wbcdmp == 1) THEN
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF (myproc == 0 ) THEN
        IF (exbccompr > 4) THEN
          !
          ! reset values inside region not used by LBC forcing
          !
          DO k=1,rayklow-2
            tmp_val = out3d(2,2,k)
            DO j = 1+ngbrz, nylg-1-ngbrz
              DO i = 1+ngbrz, nxlg-1-ngbrz
                out3d(i,j,k) = tmp_val   ! just need some value with max/min
              END DO
            END DO
          END DO
        END IF
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccompr,              &
                      'w','w-velocity','m/s',itmp,hmax,hmin)
      END IF
    END IF

    IF(ptbcdmp == 1) THEN
      CALL mpimerge3d(pt,nx,ny,nz,out3d)
      IF (myproc == 0 ) THEN
        IF (exbccompr > 4) THEN
          !
          ! reset values inside region not used by LBC forcing
          !
          DO k=1,rayklow-1
            tmp_val = out3d(2,2,k)
            DO j = 1+ngbrz, nylg-1-ngbrz
              DO i = 1+ngbrz, nxlg-1-ngbrz
                out3d(i,j,k) = tmp_val   ! just need some value with max/min
              END DO
            END DO
          END DO
        END IF
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccompr,              &
                     'pt','Potential temperature','K',itmp,hmax,hmin)
      END IF
    END IF

    IF(prbcdmp == 1) THEN
      CALL mpimerge3d(pr,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
       IF (exbccompr > 4) THEN          !
          ! reset values inside region not used by LBC forcing
          !
          DO k=1,rayklow-1
            tmp_val = out3d(2,2,k)
            DO j = 1+ngbrz, nylg-1-ngbrz
              DO i = 1+ngbrz, nxlg-1-ngbrz
                out3d(i,j,k) = tmp_val   ! just need some value with max/min
              END DO
            END DO
          END DO
        END IF
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccompr,              &
                      'p','Pressure','Pascal',itmp,hmax,hmin)
      END IF
    END IF

    IF(qvbcdmp == 1) THEN
      CALL mpimerge3d(qv,nx,ny,nz,out3d)
      IF (myproc == 0 ) THEN
        IF (exbccompr > 4) THEN
          !
          ! reset values inside region not used by LBC forcing
          !
          DO k=1,rayklow-1
            DO j = 1+ngbrz, nylg-1-ngbrz
              DO i = 1+ngbrz, nxlg-1-ngbrz
                out3d(i,j,k) = 0.0
              END DO
            END DO
          END DO
        END IF
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccompr,              &
                      'qv','Water vapor specific humidity','kg/kg',     &
                      itmp,hmax,hmin)
      END IF
    END IF

    DO nq = 1,nscalar
      IF(qscalarbcdmp(nq) > 0) THEN

        ! DTD: Turn off bit-packing for Z array if hdfcompr > 3 and mphyopt == 11
        IF(nq >= 13 .and. exbccompr > 3) THEN
          exbccomprtmp = exbccompr - 3
        ELSE
          exbccomprtmp = exbccompr
        END IF

        CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
        IF (myproc == 0 ) THEN
          IF (exbccompr > 4) THEN
          !
          ! reset values inside region not used by LBC forcing
          !
            DO k=1,rayklow-1
              DO j=1+ngbrz,nylg-1-ngbrz
                DO i=1+ngbrz,nxlg-1-ngbrz
                  out3d(i,j,k) = 0.0
                END DO
              END DO
            END DO
          END IF

          CALL hdfwrt3d(out3d,nxlg,nylg,nz,bfid,1,exbccomprtmp,         &
                      TRIM(qnames(nq)),TRIM(qdescp(nq)),'kg/kg',        &
                      itmp,hmax,hmin)
        END IF
      END IF
    END DO

    IF(myproc == 0)  CALL hdfclose(bfid,istat)

    IF (exbccompr > 3) THEN
      DEALLOCATE (itmp,stat=istat)
      DEALLOCATE (hmax,stat=istat)
      DEALLOCATE (hmin,stat=istat)
    END IF

  ELSE IF (exbcdmp == 7 .OR. exbcdmp == 8) THEN
!-----------------------------------------------------------------------
!
!  Write out in NetCDF format.
!
!-----------------------------------------------------------------------

    IF (exbcdmp == 7) THEN
      itime = 1
    ELSE
      itime = itime + 1
      istop = NINT((tstop-tstart)/thisdmp) + 1
    END IF

    IF (myproc == 0 .AND. itime == 1) THEN

!-----------------------------------------------------------------------
!
!  Define ARPS boundary file dimension and variables
!
!-----------------------------------------------------------------------

      CALL netopen(TRIM(filename(1:lfname)), 'C', bfid)

      CALL net_define_exbc(bfid,nxlg,nylg,nz,itime,dx,dy,dz,            &
               dzmin,strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,    &
               mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,     &
               ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,qvbcdmp,            &
               nscalarout,qscalarbcdmp,ctime,istat)
    END IF

    ALLOCATE(var3du(nxlg,  nylg-1,nz-1), STAT = istat)
    ALLOCATE(var3dv(nxlg-1,nylg,  nz-1), STAT = istat)
    ALLOCATE(var3dw(nxlg-1,nylg-1,nz  ), STAT = istat)

    IF( ubcdmp == 1) THEN
      CALL mpimerge3d(u,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg
              var3du(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(bfid,0,itime,'U',var3du,nxlg,nylg-1,nz-1)
      END IF
    END IF

    IF( vbcdmp == 1) THEN
      CALL mpimerge3d(v,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k = 1,nz-1
          DO j = 1,nylg
            DO i = 1,nxlg-1
              var3dv(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(bfid,0,itime,'V',var3du,nxlg-1,nylg,nz-1)
      END IF
    END IF

    IF( wbcdmp == 1) THEN
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k = 1,nz
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(bfid,0,itime,'W',var3dw,nxlg-1,nylg-1,nz)
      END IF
    END IF

    IF(ptbcdmp == 1) THEN
      CALL mpimerge3d(pt,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(bfid,0,itime,'PT',var3dw,nxlg-1,nylg-1,nz-1)
      END IF
    END IF

    IF(prbcdmp == 1) THEN
      CALL mpimerge3d(pr,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(bfid,0,itime,'P',var3dw,nxlg-1,nylg-1,nz-1)
      END IF
    END IF

    IF(qvbcdmp == 1) THEN
      CALL mpimerge3d(qv,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(bfid,0,itime,'QV',var3dw,nxlg-1,nylg-1,nz-1)
      END IF
    END IF

    DO nq = 1,nscalar

      IF(qscalarbcdmp(nq) > 0) THEN
        CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
        IF (myproc == 0) THEN
          DO k = 1,nz-1
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var3dw(i,j,k) = out3d(i,j,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(bfid,0,itime,upcase(TRIM(qnames(nq))),var3dw,nxlg-1,nylg-1,nz-1)
        END IF
      END IF

    END DO

    IF (myproc == 0 .AND.(exbcdmp == 7 .OR. itime >= istop) )           &
      CALL netclose(bfid)

    DEALLOCATE(var3du,var3dv,var3dw)

  ELSE

    ! alternate dump format ...
    IF (myproc == 0) WRITE(6,'(1x,3a)')                                 &
               'The supported exbc data dump format are ',              &
               'binary (exbcdmp = 1), HDF4 (exbcdmp = 3) and ',         &
               'NetCDF format (exbcdmp = 7).'

    CALL arpsstop('EXBC data dump format is not supported.',1)

  END IF

  RETURN
END SUBROUTINE writejoinexbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCDUMP                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcdump( nx,ny,nz,nzsoil,nstyps,                            &
           hisfmt,hisfnm,grdbas,flcmprs,                                &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           u0exb,v0exb,w0exb,pt0exb,pr0exb,qv0exb,qscalar0exb,          &
           udtexb,vdtexb,wdtexb,ptdtexb,prdtexb,qvdtexb,qscalardtexb,   &
           tem1,tem2,tem3 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dump EXBC fields interpolated to the model time in history
!  dump format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/27/94
!
!  MODIFICATION HISTORY:
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/13/2002 (J. Brotzge)
!  Added additional arrays for multiple soil schemes.
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  nx       Number of grid points in the x-direction (east/west)
!  ny       Number of grid points in the y-direction (north/south)
!  nz       Number of grid points in the vertical
!  nzsoil   Number of grid points in the soil
!
!  hisfmt
!  hisfnm
!  grdbas
!  flcmprs
!
!  u        x component of velocity (m/s)
!  v        y component of velocity (m/s)
!  w        Vertical component of Cartesian velocity (m/s)
!  ptprt    Perturbation potential temperature (K)
!  pprt     Perturbation pressure (Pascal)
!  qv       Water vapor specific humidity (kg/kg)
!  qc       Cloud water mixing ratio (kg/kg)
!  qr       Rainwater mixing ratio (kg/kg)
!  qi       Cloud ice mixing ratio (kg/kg)
!  qs       Snow mixing ratio (kg/kg)
!  qh       Hail mixing ratio (kg/kg)
!  tke      Turbulent Kinetic Energy ((m/s)**2)
!
!  kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!  kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!  ubar     Base state zonal velocity component (m/s)
!  vbar     Base state meridional velocity component (m/s)
!  wbar     Base state vertical velocity component (m/s)
!  ptbar    Base state potential temperature (K)
!  pbar     Base state pressure (Pascal)
!  rhobar   Base state density (kg/m**3)
!  qvbar    Base state water vapor specific humidity (kg/kg)
!
!  x coordinate of grid points in physical/comp. space (m)
!  y        y coordinate of grid points in physical/comp. space (m)
!  z        z coordinate of grid points in computational space (m)
!  zp       Vertical coordinate of grid points in physical space (m)
!  zpsoil   Vertical coordinate of grid points in the soil (m)
!
!  soiltyp  Soil type
!  vegtyp   Vegetation type
!  lai      Leaf Area Index
!  roufns   Surface roughness
!  veg      Vegetation fraction
!
!  tsoil    Soil temperature (K)
!  qsoil    Soil moisture (m**3/m**3)
!  wetcanp  Canopy water amount
!
!  raing    Grid supersaturation rain
!  rainc    Cumulus convective rain
!  prcrate  Precipitation rates
!
!  radfrc   Radiation forcing (K/s)
!  radsw    Solar radiation reaching the surface
!  rnflx    Net radiation flux absorbed by surface
!  radswnet Net shortwave radiation
!  radlwin  Incoming longwave radiation
!
!  usflx    Surface flux of u-momentum (kg/(m*s**2))
!  vsflx    Surface flux of v-momentum (kg/(m*s**2))
!  ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!  qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!  None
!
!  TEMPORATY WORKING ARRAY
!
!  tem1
!  tem2
!  tem3
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz
  INTEGER :: nzsoil

  INTEGER :: hisfmt
  CHARACTER (LEN=*) :: hisfnm
  INTEGER :: grdbas
  INTEGER :: flcmprs

  REAL :: u     (nx,ny,nz)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)  ! Horizontal turb. mixing coef. for
  REAL :: kmv   (nx,ny,nz)  ! Vertical turb. mixing coef. for
                            ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)  ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)  ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)  ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)  ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)  ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)  ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)  ! Base state water vapor specific humidity
                            ! (kg/kg)

  REAL :: x     (nx)        ! The x-coord. of the physical and
                            ! computational grid. Defined at u-point.
  REAL :: y     (ny)        ! The y-coord. of the physical and
                            ! computational grid. Defined at v-point.
  REAL :: z     (nz)        ! The z-coord. of the computational grid.
                            ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)  ! The physical height coordinate defined at
                            ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                            ! w-point of the soil.

  INTEGER :: nstyps
  INTEGER :: soiltyp (nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)   ! Soil type fratction
  INTEGER :: vegtyp(nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)           ! Leaf Area Index
  REAL :: roufns (nx,ny)           ! Surface roughness
  REAL :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)     ! Canopy water amount
  REAL :: snowdpth(nx,ny)             ! Snow depth (m)

  REAL :: raing(nx,ny)                ! Grid supersaturation rain
  REAL :: rainc(nx,ny)                ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Reflected shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: u0exb (nx,ny,nz)  ! External boundary u-velocity field
  REAL :: v0exb (nx,ny,nz)  ! External boundary v-velocity field
  REAL :: w0exb (nx,ny,nz)  ! External boundary w-velocity field
  REAL :: pt0exb(nx,ny,nz)  ! External boundary pt field
  REAL :: pr0exb(nx,ny,nz)  ! External boundary p field
  REAL :: qv0exb(nx,ny,nz)  ! External boundary qv field
  REAL :: qscalar0exb(nx,ny,nz,nscalar)

  REAL :: udtexb (nx,ny,nz) ! Time tendency of external boundary u
  REAL :: vdtexb (nx,ny,nz) ! Time tendency of external boundary v
  REAL :: wdtexb (nx,ny,nz) ! Time tendency of external boundary w
  REAL :: ptdtexb(nx,ny,nz) ! Time tendency of external boundary pt
  REAL :: prdtexb(nx,ny,nz) ! Time tendency of external boundary p
  REAL :: qvdtexb(nx,ny,nz) ! Time tendency of external boundary qv
  REAL :: qscalardtexb(nx,ny,nz,nscalar)

  REAL :: tem1 (nx,ny,nz)
  REAL :: tem2 (nx,ny,nz)
  REAL :: tem3 (nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, nq

  INTEGER, SAVE :: nchexbc

  REAL    :: tema
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tema = curtim - ( abstfcst0 - abstinit )

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx
        u(i,j,k) = u0exb(i,j,k) + udtexb(i,j,k) * tema
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny
      DO i = 1, nx-1
        v(i,j,k) = v0exb(i,j,k) + vdtexb(i,j,k) * tema
      END DO
    END DO
  END DO

  DO k = 1, nz
    DO j = 1, ny-1
      DO i = 1, nx-1
        w(i,j,k) = w0exb(i,j,k) + wdtexb(i,j,k) * tema
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        ptprt(i,j,k) = pt0exb(i,j,k) + ptdtexb(i,j,k) * tema
        pprt (i,j,k) = pr0exb(i,j,k) + prdtexb(i,j,k) * tema
        qv   (i,j,k) = qv0exb(i,j,k) + qvdtexb(i,j,k) * tema
!
!  Since we do not have enough tem arrays to store qctem,
!  qrtem, qitem, qstem and qhtem, we pass the model arrays
!  into dtadump.
!
!    qctem(i,j,k) = qc0exb(i,j,k) + qcdtexb(i,j,k) * tema
!    qrtem(i,j,k) = qr0exb(i,j,k) + qrdtexb(i,j,k) * tema
!    qitem(i,j,k) = qi0exb(i,j,k) + qidtexb(i,j,k) * tema
!    qstem(i,j,k) = qs0exb(i,j,k) + qsdtexb(i,j,k) * tema
!    qhtem(i,j,k) = qh0exb(i,j,k) + qhdtexb(i,j,k) * tema
      END DO
    END DO
  END DO

!    blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,dumpstride
    IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

      CALL dtadump( nx,ny,nz,nzsoil,nstyps,                             &
                  hisfmt,nchexbc,hisfnm,grdbas,filcmprs,                &
                  u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,              &
                  ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,               &
                  x,y,z,zp,zpsoil,                                      &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  tem1,tem2,tem3 )

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  RETURN
END SUBROUTINE exbcdump

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE INIT_MM_EXBC             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE init_MM_exbc(nx,ny,nz,pt,pr,qscalar)

  USE my3mom_fncs_mod, ONLY: gammaDP

  IMPLICIT NONE

  INCLUDE  'phycst.inc'
  INCLUDE  'globcst.inc'

  INTEGER :: nx,ny,nz
  REAL :: pr(nx,ny,nz)
  REAL :: ppi(nx,ny,nz)     ! Exner function
  REAL :: pt(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: rho(nx,ny,nz)     ! Air density
  REAL :: t(nx,ny,nz)       ! Air temperature


  REAL :: temNtx            ! temporary number concentration for species x
  REAL :: temZx             ! temporary reflectivity for species x
  REAL :: Gr                ! constant in equation for radar reflectivity
  REAL :: Gi
  REAL :: Gs
  REAL :: Gg
  REAL :: Gh

  REAL, PARAMETER :: pi = 3.14159265

  ! Fixed intercept parameters for rain,snow,graupel,hail

  REAL :: N0r
  REAL :: N0s
  REAL :: N0g
  REAL :: N0h

  ! Fixed cloud number concentration

  REAL :: Ntc

  ! Fixed densities for rain,snow,graupel,hail

  REAL :: rhor
  REAL :: rhoi
  REAL :: rhos
  REAL :: rhog
  REAL :: rhoh

  REAL, PARAMETER :: epsQ = 1.0e-14

  INTEGER :: i,j,k,nq

  REAL :: p0inv

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0r = n0rain
  N0s = n0snow
  N0g = n0grpl
  N0h = n0hail

  Ntc = ntcloud
  rhor = 1000.0
  rhoi = rhoice
  rhos = rhosnow
  rhog = rhogrpl
  rhoh = rhohail

  temNtx = 0.0
  temZx = 0.0

  ! Calculate exner function

!  CALL setppi(nx,ny,nz,nt,tim,pprt,pbar,ppi)

  p0inv=1./p0
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ppi(i,j,k)=(pr(i,j,k)*p0inv)**rddcp
      END DO
    END DO
  END DO

  ! Calculate air temperature and density

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        t(i,j,k) = pt(i,j,k)*ppi(i,j,k)
        rho(i,j,k) = pr(i,j,k)/(rd*t(i,j,k))
      END DO
    END DO
  END DO

  ! For 2 or 3 moment scheme, diagnose the additional moments using the initial values of mixing ratios
  ! Assume alpha=0 (the shape parameter) EDIT: 08/08/08, now uses value of alpha specified in namelist

  ! Precalculate constant in radar equation

!  Gx = 20.

  Gr = ((6.+alpharain)*(5+alpharain)*(4+alpharain))/((3.+alpharain)*(2+alpharain)*(1+alpharain))
  Gi = ((6.+alphaice)*(5+alphaice)*(4+alphaice))/((3.+alphaice)*(2+alphaice)*(1+alphaice))
  Gs = ((6.+alphasnow)*(5+alphasnow)*(4+alphasnow))/((3.+alphasnow)*(2+alphasnow)*(1+alphasnow))
  Gg = ((6.+alphagrpl)*(5+alphagrpl)*(4+alphagrpl))/((3.+alphagrpl)*(2+alphagrpl)*(1+alphagrpl))
  Gh = ((6.+alphahail)*(5+alphahail)*(4+alphahail))/((3.+alphahail)*(2+alphahail)*(1+alphahail))

  IF(mphyopt >= 9) THEN

    IF((hail_ON == 0) .and. (graupel_ON == 1)) THEN
      !BJP NOV 2010 Hail Switch (add initial hail to graupel)
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
             qscalar(i,j,k,P_QG) = qscalar(i,j,k,P_QG) + qscalar(i,j,k,p_QH)
             qscalar(i,j,k,P_QH) = 0.0
          END DO
        END DO
      END DO
    ELSE IF((graupel_ON == 0) .and. (hail_ON == 1) ) THEN
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
             qscalar(i,j,k,P_QH) = qscalar(i,j,k,P_QG) + qscalar(i,j,k,p_QH)
             qscalar(i,j,k,P_QG) = 0.0
          END DO
        END DO
      END DO
    END IF

    DO nq = 1,nscalarq
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1

            temNtx = 0.0
            temZx = 0.0
            IF(qscalar(i,j,k,nq) >= epsQ) THEN
            IF(nq == P_QC) THEN   ! Set Nc to fixed cloud number concentration
              temNtx = Ntc
            ELSE IF(nq == P_QR) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alpharain)))*(N0r**(3./(4.+alpharain)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhor* &
                       sngl(gammaDP(4.d0+dble(alpharain)))))**((1.+alpharain)/(4.+alpharain))
              temZx = Gr/(((pi/6.)*rhor)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QI) THEN
              temNtx = 5.*exp(0.304*(273.15-max(233.,t(i,j,k))))     ! Ice nucleation from Cooper's eqn.
              temZx = Gi/((440.0)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QS) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alphasnow)))*(N0s**(3./(4.+alphasnow)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhos* &
                       sngl(gammaDP(4.d0+dble(alphasnow)))))**((1.+alphasnow)/(4.+alphasnow))
              temZx = Gs/(((pi/6.)*rhos)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QG) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alphagrpl)))*(N0g**(3./(4.+alphagrpl)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhog* &
                       sngl(gammaDP(4.d0+dble(alphagrpl)))))**((1.+alphagrpl)/(4.+alphagrpl))
              temZx = Gg/(((pi/6.)*rhog)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QH) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alphahail)))*(N0h**(3./(4.+alphahail)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhoh* &
                       sngl(gammaDP(4.d0+dble(alphahail)))))**((1.+alphahail)/(4.+alphahail))
              temZx = Gh/(((pi/6.)*rhoh)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            END IF
            END IF
            qscalar(i,j,k,nq+nscalarq) = temNtx
            IF(mphyopt == 11) THEN
              qscalar(i,j,k,nq+nscalarq+(nscalarq-1)) = temZx
            END IF
          END DO
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE init_MM_exbc

