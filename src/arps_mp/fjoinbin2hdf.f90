!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE JOINBIN2HDF                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE joinbin2hdf (fileheader,nx,ny,nz,nzsoil,nstyp,grdbas,hdfcompr)

!
!-----------------------------------------------------------------------
!
!  To join together a set of ARPS history or data files produced by the
!  processors of MPP machines with message passing.
!
!  Input data file is in binary format and the output is in HDF4 format
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (05/16/2002)
!  Based on joindumps.f90
!
!  MODIFICATION HISTORY.
!
!  07/18/2002 (Yunheng Wang)
!  Modified to be compatible with the new soil model (since IHOP_3)
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INCLUDE 'mp.inc'
  INCLUDE 'indtflg.inc'
  INCLUDE 'phycst.inc'
!
!--------------------------------------------------------------------
!
! PARAMETERS
!
!--------------------------------------------------------------------

  CHARACTER (LEN=*) :: fileheader

  INTEGER :: nx,ny,nz,nzsoil,nstyp
  INTEGER, INTENT(IN) :: grdbas

  INTEGER :: hdfcompr

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=40), PARAMETER :: fmtverbin410='004.10 Binary Data'
  CHARACTER(LEN=40), PARAMETER :: fmtverbin500='005.00 Binary Data'
  CHARACTER(LEN=40), PARAMETER :: fmtverbin530='005.30 Binary Data'
  CHARACTER(LEN=40), PARAMETER :: fmtverhdf410='004.10 HDF4 Coded Data'
  CHARACTER(LEN=40), PARAMETER :: fmtverhdf500='005.00 HDF4 Coded Data'
  CHARACTER(LEN=40), PARAMETER :: fmtverhdf530='005.30 HDF4 Coded Data'
  INTEGER, PARAMETER :: intver410 = 410
  INTEGER, PARAMETER :: intver500 = 500
  INTEGER, PARAMETER :: intver530 = 530

  INTEGER :: nxlg, nylg, nzlg
  INTEGER :: nxin, nyin, nzin, nzsoilin

  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: outfile
  CHARACTER (LEN=256) :: outfile_old

  INTEGER, PARAMETER :: unit0=110, maxunit=60
  INTEGER, PARAMETER :: junit0=11
  INTEGER :: sd_ido, stg_dim

  INTEGER :: lenstr, joff
  INTEGER :: fi, fj, i, j, k
  INTEGER :: ii,jj,iiend
  INTEGER :: intver
  CHARACTER (LEN=40) :: fmtver
  CHARACTER (LEN=80) :: runname
  CHARACTER (LEN=10) :: tmunit
  CHARACTER (LEN=12) :: label
  CHARACTER (LEN=10) :: varname

  INTEGER :: nocmnt
  CHARACTER (LEN=80), ALLOCATABLE :: cmnt(:)

  REAL :: curtim
  INTEGER :: i01, i02, i03, i04, i05, i06, i07, i08, i09, i10
  INTEGER :: i11, i12, i13, i14, i15, i16, i17, i18, i19, i20
  REAL :: r01, r02, r03, r04, r05, r06, r07, r08, r09, r10
  REAL :: r11, r12, r13, r14, r15, r16, r17, r18, r19, r20
  INTEGER :: idummy

  REAL, ALLOCATABLE :: xlg(:), ylg(:), z(:)
  REAL, ALLOCATABLE :: xsm(:), ysm(:)
  REAL, ALLOCATABLE :: a3dlg(:,:,:), a3dsm(:,:,:)
  REAL, ALLOCATABLE :: a2dlg(:,:), a2dsm(:,:)

  REAL, ALLOCATABLE :: a3dsoillg(:,:,:),a3dsoilsm(:,:,:)
  REAL, ALLOCATABLE :: qsoil(:,:,:,:), tsoil(:,:,:,:), wetcanp(:,:,:)
  REAL, AlLOCATABLE :: stypfrct(:,:,:)
  INTEGER, ALLOCATABLE :: soiltyp(:,:,:)

  INTEGER, ALLOCATABLE :: ai2dlg(:,:), ai2dsm(:,:)

  INTEGER, ALLOCATABLE :: i0(:,:), j0(:,:)

  INTEGER, ALLOCATABLE :: iunit(:)
  INTEGER, ALLOCATABLE :: ffi(:), ffj(:)

  INTEGER :: ierr, istat
  LOGICAL :: fexist
  INTEGER :: landflg, sfcflg
  INTEGER :: is, nstypvar

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
!  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmpsoil(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmpsoil4d(:,:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp2d(:,:)
  REAL, ALLOCATABLE :: hmax(:), hmin(:)
  REAL, ALLOCATABLE :: hmaxsoil(:), hminsoil(:)

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

  ALLOCATE(xlg(nxlg))
  ALLOCATE(ylg(nylg))
  ALLOCATE(z(nzlg))
  ALLOCATE(xsm(nx))
  ALLOCATE(ysm(ny))
  ALLOCATE(a3dsoillg(nxlg,nylg,nzsoil))
  ALLOCATE(a3dsoilsm(nx,ny,nzsoil))
  ALLOCATE(a3dlg(nxlg,nylg,nzlg))
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

  ALLOCATE(itmp(nxlg, nylg, nz))
!  ALLOCATE(itmpsoil(nxlg, nylg, nzsoil))
  ALLOCATE(itmpsoil4d(nxlg, nylg, nzsoil,0:nstyp))
  ALLOCATE(itmp2d(nxlg, nylg))
  ALLOCATE(hmax(nz))
  ALLOCATE(hmin(nz))
  ALLOCATE(hmaxsoil(nzsoil))
  ALLOCATE(hminsoil(nzsoil))

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

  a3dlg      = 0.0
  a3dsoillg  = 0.0
  a2dlg      = 0.0
  ai2dlg     = 0

  DO jj = 1,1+(nproc_x*nproc_y-1)/maxunit

    nstyp = 0
    is = 0
    nstypvar = 0

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
      WRITE(outfile, '(a,a,i3.3)') fileheader(1:lenstr),'_tmp',iiend
    END IF
    CALL asnfile(outfile, '-F f77 -N ieee', ierr)

    IF (iiend == nproc_x*nproc_y) THEN
      i = index(outfile, 'bin', .TRUE.)
      IF (i <=0 .OR. i == 3) THEN
        WRITE(6,*) "WARNING: Are you sure the input data file is in binary format?"
        WRITE(6,*) "Filename ", outfile, " may be consturcted incorrectly."
      ELSE
        WRITE(outfile,'(a)') outfile(1:i-1)//'hdf'//outfile(i+3:LEN_TRIM(outfile))
      END IF

      CALL hdfopen(outfile,2,sd_ido)
      IF (sd_ido < 0) THEN
        WRITE (6,*) "JOINBIN2HDF: ERROR creating HDF4 file: ", outfile
        CALL arpsstop('arpsstop called from JOINBIN2HDF',1)
      END IF
    ELSE
      OPEN (UNIT=junit0+joff,FILE=outfile,FORM='unformatted')
    END IF

    IF (joff > 0 ) &
      OPEN (UNIT=junit0+joff-1,FILE=outfile_old,FORM='unformatted')

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

    IF (fmtver == fmtverbin410) THEN
      intver = intver410
      WRITE(*,'(a,a)') "WARNING: At present,JOINBIN2HDF only support verion 5.00", &
                       " (bin --> HDF)"
      WRITE(*,'(a,a,a)') "         For verion 4.10 (bin --> HDF), please use ", &
                         "JOINBIN2HDF from a previous ARPS version (like ",     &
                         "arps5.0.0IHOP_2b or earlier). "
      CALL arpsstop("Data format incompatible!!!", 1)
    ELSE IF (fmtver == fmtverbin500) THEN
      intver = intver500
    ELSE IF (fmtver == fmtverbin530) THEN
      intver = intver530
    ELSE
      WRITE(*,*) "ERROR: Data format mismatch."
      WRITE(*,*) "       Expected is ", fmtverbin530
      WRITE(*,*) "       Read in from file [", TRIM(filename), "] is ", fmtver
      CALL arpsstop("Data format incompatible!!!", 1)
    END IF

    IF (iiend == nproc_x*nproc_y)  THEN
      CALL hdfwrtc(sd_ido, 40, 'fmtver', fmtverhdf530, istat)
    END IF

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) runname
    END DO
    IF (iiend == nproc_x*nproc_y) &
      CALL hdfwrtc(sd_ido, 80, 'runname', runname, istat)

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) nocmnt
    END DO
    IF (iiend == nproc_x*nproc_y) &
      CALL hdfwrti(sd_ido, 'nocmnt', nocmnt, istat)

    IF (jj == 1) ALLOCATE (cmnt(nocmnt))

    IF ( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        DO ii=1+(jj-1)*maxunit,iiend
          READ (iunit(ii)) cmnt(i)
        END DO
      END DO
      IF (iiend == nproc_x*nproc_y) &
        CALL hdfwrtc(sd_ido, 80*nocmnt, 'cmnt', cmnt, istat)
    END IF


    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) curtim,tmunit
    END DO
    IF (iiend == nproc_x*nproc_y) THEN
      CALL hdfwrtc(sd_ido, 7, 'tmunit', tmunit, istat)
      CALL hdfwrtr(sd_ido, 'time', curtim, istat)
    END IF
!
!-----------------------------------------------------------------------
!
!  Read/write dimensions of data in binary file and check against
!  the dimensions passed to JOINBIN2HDF
!
!-----------------------------------------------------------------------
!

    DO ii=1+(jj-1)*maxunit,iiend
      READ (iunit(ii)) nxin,nyin,nzin, nzsoilin
    END DO
    IF ((nxin /= nx).OR.(nyin /= ny).OR.(nzin /= nz) .OR. (nzsoilin /= nzsoil)) THEN
      WRITE (*,*) "ERROR:  missmatch in sizes."
      WRITE (*,*) "nxin,nyin,nzin, nzsoilin",nxin,nyin,nzin, nzsoilin
      WRITE (*,*) "nx,ny,nz,nzsoil",nx,ny,nz, nzsoil
      STOP
    END IF
    IF (iiend == nproc_x*nproc_y) THEN
      CALL hdfwrti(sd_ido, 'nx', nxlg, istat)
      CALL hdfwrti(sd_ido, 'ny', nylg, istat)
      CALL hdfwrti(sd_ido, 'nz', nzlg, istat)
      CALL hdfwrti(sd_ido, 'nzsoil', nzsoil, istat)
    END IF

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
              i06, i07, i08, i09, i10,                                  &
              i11, i12, i13, i14, i15,                                  &
              i16, i17, i18, i19, i20
    END DO

    landflg = i09
    sfcflg = i07
    mstin = i04
    icein = i05

    IF (iiend == nproc_x*nproc_y) THEN
      CALL hdfwrti(sd_ido, 'grdflg', i01, istat)
      CALL hdfwrti(sd_ido, 'basflg', i02, istat)
      CALL hdfwrti(sd_ido, 'varflg', i03, istat)
      CALL hdfwrti(sd_ido, 'mstflg', i04, istat)
      CALL hdfwrti(sd_ido, 'iceflg', i05, istat)
      CALL hdfwrti(sd_ido, 'trbflg', i06, istat)
      CALL hdfwrti(sd_ido, 'sfcflg', i07, istat)
      CALL hdfwrti(sd_ido, 'rainflg',i08, istat)
      CALL hdfwrti(sd_ido, 'landflg',i09, istat)
      CALL hdfwrti(sd_ido, 'totflg', i10, istat)
      CALL hdfwrti(sd_ido, 'tkeflg', i11, istat)
      CALL hdfwrti(sd_ido, 'mapproj', i14, istat)
      CALL hdfwrti(sd_ido, 'month',  i15, istat)
      CALL hdfwrti(sd_ido, 'day',    i16, istat)
      CALL hdfwrti(sd_ido, 'year',   i17, istat)
      CALL hdfwrti(sd_ido, 'hour',   i18, istat)
      CALL hdfwrti(sd_ido, 'minute', i19, istat)
      CALL hdfwrti(sd_ido, 'second', i20, istat)
    END IF

    IF (grdbas /= 1 ) THEN
      IF (intver >= intver530) THEN   ! newer version 5.3 or later
        DO ii = 1+(jj-1)*maxunit, iiend
          READ(iunit(ii))                                               &
             p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,                 &
             p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,                 &
             p_zrin,p_ziin,p_zsin,p_zsin,p_zhin,idummy,                 &
             idummy,idummy,idummy,idummy,idummy,idummy,                 &
             idummy,idummy,idummy,idummy,idummy,idummy
        END DO
        nscalarin = i12
      ELSE                      ! version earlier than 5.3
        p_qcin=0; p_qrin=0; p_qiin=0; p_qsin=0; p_qgin=0; p_qhin=0
        p_ncin=0; p_nrin=0; p_niin=0; p_nsin=0; p_nhin=0; p_ngin=0
                  p_zrin=0; p_ziin=0; p_zsin=0; p_zhin=0; p_zgin=0

        nscalarin = 0
        IF (mstin == 1) THEN
          p_qcin=1
          p_qrin=2
          nscalarin = nscalarin + 2
          IF (icein == 2) THEN
            p_qiin=3
            p_qsin=4
            p_qhin=5
            nscalarin = nscalarin + 3
          END IF
        END IF
      END IF

      IF (iiend == nproc_x*nproc_y) THEN
        CALL hdfwrti(sd_ido, 'nscalar', nscalarin, istat)
        CALL hdfwrti(sd_ido, 'P_QC', p_qcin, istat)
        CALL hdfwrti(sd_ido, 'P_QR', p_qrin, istat)
        CALL hdfwrti(sd_ido, 'P_QI', p_qiin, istat)
        CALL hdfwrti(sd_ido, 'P_QS', p_qsin, istat)
        CALL hdfwrti(sd_ido, 'P_QG', p_qgin, istat)
        CALL hdfwrti(sd_ido, 'P_QH', p_qhin, istat)
        CALL hdfwrti(sd_ido, 'P_NC', p_ncin, istat)
        CALL hdfwrti(sd_ido, 'P_NR', p_nrin, istat)
        CALL hdfwrti(sd_ido, 'P_NI', p_niin, istat)
        CALL hdfwrti(sd_ido, 'P_NS', p_nsin, istat)
        CALL hdfwrti(sd_ido, 'P_NG', p_ngin, istat)
        CALL hdfwrti(sd_ido, 'P_NH', p_nhin, istat)
        CALL hdfwrti(sd_ido, 'P_ZR', p_zrin, istat)
        CALL hdfwrti(sd_ido, 'P_ZI', p_ziin, istat)
        CALL hdfwrti(sd_ido, 'P_ZS', p_zsin, istat)
        CALL hdfwrti(sd_ido, 'P_ZG', p_zgin, istat)
        CALL hdfwrti(sd_ido, 'P_ZH', p_zhin, istat)
      END IF

    END IF

    IF (intver < intver530) THEN
      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                                 &
                r01, r02, r03, r04, r05,                                 &
                r06, r07, r08, r09, r10,                                 &
                r11, r12, r13, r14, r15,                                 &
                r16, r17, r18, r19, r20
      END DO
    ELSE
      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                                 &
                r01, r02, r03, r04, r05,                                 &
                r06, r07, r08, ntcloud, n0rain,                          &
                n0snow,n0grpl,n0hail,rhoice,rhosnow,                     &
                rhogrpl,rhohail,alpharain,alphaice,alphasnow,            &
                alphagrpl,alphahail,r09,r10,r11,                         &
                r11, r12, r13, r14, r15,                                 &
                r16, r17, r18, r19, r20
      END DO
    END IF

    IF (iiend == nproc_x*nproc_y) THEN
      CALL hdfwrtr(sd_ido, 'umove',   r01, istat)
      CALL hdfwrtr(sd_ido, 'vmove',   r02, istat)
      CALL hdfwrtr(sd_ido, 'xgrdorg', r03, istat)
      CALL hdfwrtr(sd_ido, 'ygrdorg', r04, istat)

      CALL hdfwrtr(sd_ido, 'trulat1', r05, istat)
      CALL hdfwrtr(sd_ido, 'trulat2', r06, istat)
      CALL hdfwrtr(sd_ido, 'trulon',  r07, istat)
      CALL hdfwrtr(sd_ido, 'sclfct',  r08, istat)
      CALL hdfwrtr(sd_ido, 'tstop',   r16, istat)
      CALL hdfwrtr(sd_ido, 'thisdmp', r17, istat)
      CALL hdfwrtr(sd_ido, 'latitud', r18, istat)
      CALL hdfwrtr(sd_ido, 'ctrlat',  r19, istat)
      CALL hdfwrtr(sd_ido, 'ctrlon',  r20, istat)

      CALL hdfwrtr(sd_ido, 'ntcloud', ntcloud, istat)
      CALL hdfwrtr(sd_ido, 'n0rain',  n0rain,  istat)
      CALL hdfwrtr(sd_ido, 'n0snow',  n0snow,  istat)
      CALL hdfwrtr(sd_ido, 'n0grpl',  n0grpl,  istat)
      CALL hdfwrtr(sd_ido, 'n0hail',  n0hail,  istat)
      CALL hdfwrtr(sd_ido, 'rhoice',  rhoice,  istat)
      CALL hdfwrtr(sd_ido, 'rhosnow', rhosnow, istat)
      CALL hdfwrtr(sd_ido, 'rhogrpl', rhogrpl, istat)
      CALL hdfwrtr(sd_ido, 'rhohail', rhohail, istat)
      CALL hdfwrtr(sd_ido, 'alpharain', alpharain, istat)
      CALL hdfwrtr(sd_ido, 'alphaice',  alphaice,  istat)
      CALL hdfwrtr(sd_ido, 'alphasnow', alphasnow, istat)
      CALL hdfwrtr(sd_ido, 'alphagrpl', alphagrpl, istat)
      CALL hdfwrtr(sd_ido, 'alphahail', alphahail, istat)
    END IF

    IF (i10 == 1) THEN             ! totout == 1

      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                               &
              i01, i02, i03, i04, i05,                                 &
              i06, i07, i08, i09, i10,                                 &
              i11, i12, i13, i14, i15,                                 &
              i16, i17, i18, i19, i20
      END DO

      nstyp = i01

      IF (nstyp < 1) nstyp = 1
      IF (jj == 1) THEN           ! 1st output unit, do allocation
          ALLOCATE(soiltyp(nxlg, nylg, nstyp))
          AlLOCATE(stypfrct(nxlg,nylg, nstyp))
          ALLOCATE(tsoil(  nxlg, nylg, nzsoil, 0:nstyp))
          tsoil = 0.0
          ALLOCATE(qsoil(  nxlg, nylg, nzsoil, 0:nstyp))
          qsoil = 0.0
          ALLOCATE(wetcanp(nxlg, nylg, 0:nstyp))
          wetcanp = 0.0
      END IF

      IF (iiend == nproc_x*nproc_y) THEN
        CALL hdfwrti(sd_ido, 'nstyp',  i01, istat)
        CALL hdfwrti(sd_ido, 'prcflg', i02, istat)
        CALL hdfwrti(sd_ido, 'radflg', i03, istat)
        CALL hdfwrti(sd_ido, 'flxflg', i04, istat)
        CALL hdfwrti(sd_ido, 'snowflg',i06, istat)

      END IF

      DO ii=1+(jj-1)*maxunit,iiend
        READ (iunit(ii))                                               &
              r01, r02, r03, r04, r05,                                 &
              r06, r07, r08, r09, r10,                                 &
              r11, r12, r13, r14, r15,                                 &
              r16, r17, r18, r19, r20
      END DO

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

    WRITE(varname,'(a)') label(1:8)
    DO i = 1,8
     IF(varname(i:i)== " ") EXIT
    END DO
    WRITE(varname,'(a)') varname(1:i-1)

    IF (LEN_TRIM(varname) < 1) THEN
      WRITE(6,'(3a)') 'Can not determind variable name from the label',&
               label, 'Program stopped at JOINBIN2HDF.'
      CALL arpsstop('arpsstop called from JOINBIN2HDF',1)
    END IF

    SELECT CASE (label(12:12))
      CASE ('1')
        stg_dim = 1
      CASE ('2')
        stg_dim = 2
      CASE ('3')
        stg_dim = 3
      CASE DEFAULT
        stg_dim = 0
    END SELECT
!
!--------------------------------------------------------------------------
!
!  Please noted that BIN and HDF have differenct label for:
!
!  prcrat1(prcrate1), prcrat2(prcrate2), prcrat3(prcrate3), prcrat4(prcrate4)
!  and stypfrc(stypfrct)
!
!--------------------------------------------------------------------------

    IF (varname(1:6)== "prcrat") THEN
      WRITE(varname, '(a)') "prcrate"//varname(7:7)
    ELSE IF (varname == "stypfrc") THEN
      WRITE(varname, '(a)') "stypfrct"
    END IF

    WRITE(6,*) "JOINBIN2HDF: ", varname, "being joined. Label in: ", label

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
        IF (iiend == nproc_x*nproc_y) THEN
          CALL hdfwrt1d(xlg,nxlg,sd_ido,'x','x coordinate','m')
!          WRITE(6,*) "Writing x"
        ELSE
          WRITE (junit0+joff) xlg
        END IF

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
        IF (iiend == nproc_x*nproc_y) THEN
          CALL hdfwrt1d(ylg,nylg,sd_ido,'y','y coordinate','m')
!          WRITE(6,*) "Writing y"
        ELSE
          WRITE (junit0+joff) ylg
        END IF

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
        IF (iiend == nproc_x*nproc_y) THEN
          CALL hdfwrt1d(z,nzlg,sd_ido,'z','z coordinate','m')
!          WRITE(6,*) "Writing z"
        ELSE
          WRITE (junit0+joff) z
        END IF

      ELSE
        GO TO 330
      END IF
    ELSE IF (label(10:10) == "2") THEN
      IF (label(9:9) == "r") THEN             ! 2-d real

!---------------------------------------------------------------------
!
! Soil variables (stypfrct, wetcanp)
!
!---------------------------------------------------------------------
      SELECT CASE (varname)

      CASE ("wetcanp")
        IF (sfcflg /= 1) THEN
          WRITE (*,*) "JOINBIN2HDF: Soil variable(",varname,           &
             ")  output is mismatch with sfcflg (", sfcflg,")."
          CALL arpsstop ("ARPSSTOP called from JOINBIN2HDF.", 1)
        END IF

        nstypvar = nstypvar+1

        IF (joff > 0 .AND. nstypvar < 4 ) READ (junit0+joff-1) wetcanp
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) a2dsm
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              wetcanp(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), (nstypvar-1)/3) =   &
                         a2dsm(i,j)
            END DO
          END DO
        END DO

        IF ((nstypvar-1)/3+1 > nstyp .OR. nstyp <= 1) THEN
          IF (iiend == nproc_x*nproc_y) THEN
            CALL hdfwrt3d(wetcanp,nxlg,nylg,nstyp+1,sd_ido,0,hdfcompr,       &
                  'wetcanp',' ',' ', itmp,hmax,hmin)
!            WRITE(6,*) "Writing wetcanp"
          ELSE
            WRITE (junit0+joff) wetcanp
          END IF
        END IF

      CASE ("stypfrct")
        IF (landflg /= 1) THEN
          WRITE (*,*) "JOINBIN2HDF: Soil fraction output is mismatch  &
              & with landflg (", landflg,")."
          CALL arpsstop ("ARPSSTOP called from JOINBIN2HDF.", 1)
        END IF

        IF (joff > 0 .AND. is == 1 ) READ (junit0+joff-1) stypfrct
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) a2dsm
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              stypfrct(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), is) = a2dsm(i,j)
            END DO
          END DO
        END DO

        IF (is >= nstyp) THEN
          IF (iiend == nproc_x*nproc_y) THEN
            CALL hdfwrt3d(stypfrct,nxlg,nylg,nstyp,sd_ido,0,hdfcompr,       &
                  'stypfrct',' ',' ', itmp,hmax,hmin)
!            WRITE(6,*) "Writing stypfrct"
          ELSE
            WRITE (junit0+joff) stypfrct
          END IF
        END IF

      CASE DEFAULT
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
        IF (iiend == nproc_x*nproc_y) THEN
           CALL hdfwrt2d(a2dlg,nxlg,nylg,sd_ido,0,hdfcompr,            &
                  varname,' ',' ',itmp2d)
!            WRITE(6,*) "Writing ", varname
        ELSE
          WRITE (junit0+joff) a2dlg
        END IF

      END SELECT

      ELSE IF (label(9:9) == "i") THEN        ! 2-d integer

!---------------------------------------------------------------------
!
! Soiltyp
!
!---------------------------------------------------------------------
      IF (varname == "soiltyp") THEN

        IF (landflg /= 1) THEN
          WRITE (*,*) "JOINBIN2HDF: Soil type output is mismatch  &
              & with landflg (", landflg,")."
          CALL arpsstop ("ARPSSTOP called from JOINBIN2HDF.", 1)
        END IF

        is = is + 1
        IF (joff > 0 .AND. is == 1 ) READ (junit0+joff-1) soiltyp
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) ai2dsm
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              soiltyp(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), is) = ai2dsm(i,j)
            END DO
          END DO
        END DO

!        WRITE(6,*) " istype: ", is, " of ", nstyp

        IF (is >= nstyp) THEN
          IF (iiend == nproc_x*nproc_y) THEN
            CALL hdfwrt3di(soiltyp,nxlg,nylg,nstyp,sd_ido,0,hdfcompr,       &
                  'soiltyp',' ',' ')
!            WRITE(6,*) "Writing soiltyp"
          ELSE
            WRITE (junit0+joff) soiltyp
          END IF
        END IF

     ELSE
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
        IF (iiend == nproc_x*nproc_y) THEN
          CALL hdfwrt2di(ai2dlg,nxlg,nylg,sd_ido,0,0,                            &
                  varname,' ',' ')
!            WRITE(6,*) "Writing ", varname

        ELSE
          WRITE (junit0+joff) ai2dlg
        END IF

      END IF     ! varname = "soiltyp"

      ELSE
        GO TO 330
      END IF
    ELSE IF (label(10:10) == "3") THEN        ! 3-d
!
!---------------------------------------------------------------------
!
! Soil variables (qsoil, tsoil)
!
!---------------------------------------------------------------------
      SELECT CASE (varname)

      CASE ("tsoil")
        IF (sfcflg /= 1) THEN
          WRITE (*,*) "JOINBIN2HDF: Soil variable(",varname,           &
             ")  output is mismatch with sfcflg (", sfcflg,")."
          CALL arpsstop ("ARPSSTOP called from JOINBIN2HDF.", 1)
        END IF

        nstypvar = nstypvar+1

        IF (joff > 0 .AND. nstypvar < 4 ) READ (junit0+joff-1) tsoil
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) a3dsoilsm
          DO k=1,nzsoil
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              tsoil(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3),k,(nstypvar-1)/3) =  &
                         a3dsoilsm(i,j,k)
            END DO
          END DO
          END DO
        END DO

        IF ((nstypvar-1)/3+1 > nstyp .OR. nstyp <= 1) THEN
          IF (iiend == nproc_x*nproc_y) THEN
            CALL hdfwrt4d(tsoil,nxlg,nylg,nzsoil,nstyp+1,sd_ido,0,hdfcompr,   &
                  'tsoil','Soil temperature','K',                        &
                  itmpsoil4d,hmaxsoil,hminsoil)
!            WRITE(6,*) "Writing tsoil"
          ELSE
            WRITE (junit0+joff) tsoil
          END IF
        END IF

      CASE ("qsoil")
        IF (sfcflg /= 1) THEN
          WRITE (*,*) "JOINBIN2HDF: Soil variable(",varname,           &
             ")  output is mismatch with sfcflg (", sfcflg,")."
          CALL arpsstop ("ARPSSTOP called from JOINBIN2HDF.", 1)
        END IF

        nstypvar = nstypvar+1

        IF (joff > 0 .AND. nstypvar < 4 ) READ (junit0+joff-1) qsoil
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) a3dsoilsm
          DO k=1,nzsoil
          DO j=j0(fi,fj),ny
            DO i=i0(fi,fj),nx
              qsoil(i+(fi-1)*(nx-3),j+(fj-1)*(ny-3),k,(nstypvar-1)/3) =  &
                         a3dsoilsm(i,j,k)
            END DO
          END DO
          END DO
        END DO

        IF ((nstypvar-1)/3+1 > nstyp .OR. nstyp <= 1) THEN
          IF (iiend == nproc_x*nproc_y) THEN
            CALL hdfwrt4d(qsoil,nxlg,nylg,nzsoil,nstyp+1,sd_ido,0,hdfcompr,   &
                  'qsoil','Soil moisture','fraction',                         &
                  itmpsoil4d,hmaxsoil,hminsoil)
!            WRITE(6,*) "Writing qsoil"
          ELSE
            WRITE (junit0+joff) qsoil
          END IF
        END IF

      CASE ("zpsoil")

        IF (joff > 0 ) READ (junit0+joff-1) a3dsoillg
        DO ii=1+(jj-1)*maxunit,iiend
          fi = ffi(ii)
          fj = ffj(ii)
          READ (iunit(ii)) a3dsoilsm
          DO k = 1,nzsoil
            DO j=j0(fi,fj),ny
              DO i=i0(fi,fj),nx
                a3dsoillg(i+(fi-1)*(nx-3), j+(fj-1)*(ny-3), k) =            &
                           a3dsoilsm(i,j,k)
              END DO
            END DO
          END DO
        END DO
        IF (iiend == nproc_x*nproc_y) THEN
          CALL hdfwrt3d(a3dsoillg,nxlg,nylg,nzsoil,sd_ido,stg_dim,hdfcompr,       &
                    varname,' ',' ', itmp,hmax,hmin)
!          WRITE(6,*) "Writing ", varname
        ELSE
          WRITE (junit0+joff) a3dsoillg
        END IF

      CASE DEFAULT

!----------------------------------------------------------------------
!
!  3-d real array.
!
!----------------------------------------------------------------------
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
        IF (iiend == nproc_x*nproc_y) THEN
          CALL hdfwrt3d(a3dlg,nxlg,nylg,nz,sd_ido,stg_dim,hdfcompr,       &
                    varname,' ',' ', itmp,hmax,hmin)
!          WRITE(6,*) "Writing ", varname
        ELSE
          WRITE (junit0+joff) a3dlg
        END IF

      END SELECT

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

    IF (iiend == nproc_x*nproc_y) THEN
      CALL hdfclose(sd_ido,istat)
      IF (istat == 0) THEN
        WRITE(*,'(/a,a)') "JOINBIN2HDF: Successfully dump ", trim(outfile)
      ELSE
        WRITE(*,*) "JOINBIN2HDF: ERROR (status=", istat, ") closing ", trim(outfile)
      END IF
    ELSE
      CLOSE (junit0+joff)
    END IF

    IF (joff > 0) THEN
      CLOSE (junit0+joff-1,STATUS='delete')
    END IF

    joff = joff + 1
  END DO  ! jj

!-----------------------------------------------------------------------
!
! DEALLOCATE the arrays
!
!----------------------------------------------------------------------

  DEALLOCATE(xlg, ylg, z, STAT= istat)
  DEALLOCATE(xsm, ysm,    STAT= istat)
  DEALLOCATE(a3dlg, a3dsm, STAT= istat)
  DEALLOCATE(a2dlg, a2dsm, STAT= istat)
  DEALLOCATE(ai2dlg, ai2dsm, STAT= istat)
  DEALLOCATE(a3dsoillg,      STAT= istat)
  DEALLOCATE(a3dsoilsm,      STAT= istat)
  DEALLOCATE(i0, j0,         STAT= istat)

  DEALLOCATE(iunit, ffi, ffj, STAT= istat)

  DEALLOCATE(itmp, itmp2d, hmax, hmin,     STAT= istat)
!  DEALLOCATE(itmpsoil, hmaxsoil, hminsoil, STAT= istat)
  DEALLOCATE(itmpsoil4d, hmaxsoil, hminsoil, STAT= istat)

  DEALLOCATE(cmnt)

  IF (i10 == 1) THEN
    DEALLOCATE(soiltyp, stypfrct,  STAT= istat)
    DEALLOCATE(tsoil, qsoil, wetcanp, STAT= istat)
  END IF

  RETURN

!
!----------------------------------------------------------------------
!
!  Error with the label.
!
!----------------------------------------------------------------------
!
  330   CONTINUE

  WRITE(6,'(a,a)') ' Error with label in JOINBIN2HDF:',label
  STOP 330


END SUBROUTINE joinbin2hdf

