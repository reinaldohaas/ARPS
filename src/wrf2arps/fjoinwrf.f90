!
!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE check_files_dimensions          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE check_files_dimensions(MAXWRFFIL,grid_id,fileconv,           &
                colon,io_form,magnitude_processor,splitfile,            &
                nprocs,nproc_x,nproc_y,abstimes,abstimei,               &
                abstimee,dir_extd,extdname,nextdfil,                    &
                paddingx,paddingy,istride,jstride,                      &
                ids,ide,idss,idse,jds,jde,jdss,jdse,                    &
                kps,kpe,kpss,kpse,dbglvl,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE: Check the existence of WRF files to be read and return the
!          valid file number, file names and the domain grid indices.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (04/26/2007)
!
!  MODIFICATION HISTORY:
!  William.Gustafson@pnl.gov, 16-Feb-2009: Added nocolons option
!
!  04/16/2012 (Y. Wang)
!  Improved attadj option to read extra column or extra row of patches.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,            INTENT(IN)  :: MAXWRFFIL
  INTEGER,            INTENT(IN)  :: grid_id
  CHARACTER(LEN=1),   INTENT(IN)  :: colon
  INTEGER,            INTENT(IN)  :: magnitude_processor
  LOGICAL,            INTENT(IN)  :: splitfile
  INTEGER,            INTENT(IN)  :: io_form
  INTEGER,            INTENT(IN)  :: fileconv
  INTEGER,            INTENT(IN)  :: abstimes, abstimei, abstimee
  INTEGER,            INTENT(IN)  :: nproc_x,nproc_y
  INTEGER,            INTENT(IN)  :: nprocs(nproc_x*nproc_y)
  CHARACTER(LEN=256), INTENT(IN)  :: dir_extd
  CHARACTER(LEN=256), INTENT(OUT) :: extdname(MAXWRFFIL)
  INTEGER,            INTENT(OUT) :: nextdfil
  LOGICAL,            INTENT(IN)  :: paddingx, paddingy
  INTEGER,            INTENT(IN)  :: istride, jstride
  INTEGER,            INTENT(OUT) :: ids, ide, jds, jde
  INTEGER,            INTENT(OUT) :: idss,idse,jdss,jdse
  INTEGER,            INTENT(OUT) :: kps, kpe, kpss, kpse
  INTEGER,            INTENT(IN)  :: dbglvl
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: year, month, day, hour, minute, second
  INTEGER :: ifile, npx, npy, n
  INTEGER :: ips, ipe, jps, jpe, ipss, ipse, jpss, jpse
  INTEGER :: ipssv, ipesv, jpssv, jpesv
  INTEGER :: imod, jmod
  INTEGER :: nx,ny

  CHARACTER(LEN=256) :: tmpstr
  INTEGER :: abss

  LOGICAL :: fexist
  LOGICAL :: dset = .FALSE., in_a_row = .FALSE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begining of executable code ....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ids  = 99999999
  ide  = 0
  idss = 99999999
  idse = 0

  jds  = 99999999
  jde  = 0
  jdss = 99999999
  jdse = 0

  nextdfil    = 0
  extdname(:) = ' '
  istatus     = 0

  IF (dbglvl > 0) WRITE(6,'(1x,a,/,1x,a,/)')                            &
              '============================','WRF files to be read are:'

  ifile = abstimes
  DO WHILE (ifile <= abstimee)

    IF (ifile < 30*24*3600) THEN
      year = 1
      month = 1
      day = ifile/(24*3600)
      abss = MOD(ifile,24*3600)
      hour = abss / 3600
      abss = MOD( abss, 3600 )
      minute = abss / 60
      second = MOD( abss, 60 )
    ELSE
      CALL abss2ctim(ifile,year,month,day,hour,minute,second)
    END IF

    nextdfil = nextdfil + 1

    CALL get_wrf_filename(fileconv,dir_extd,grid_id,                    &
                          year,month,day,hour,minute,second,            &
                          colon,magnitude_processor,1,.false.,          &
                          extdname(nextdfil),istatus)

    ipssv = 0
    ipesv = 0
    jpssv = 0
    jpesv = 0

    n = 0
    DO npy = 1,nproc_y
      in_a_row = .FALSE.
      DO npx = 1,nproc_x
        IF (npx > 1) in_a_row = .TRUE.

        n = n+1

        CALL get_wrf_filename(fileconv,dir_extd,grid_id,                &
                            year,month,day,hour,minute,second,          &
                            colon,magnitude_processor,nprocs(n),splitfile, &
                            tmpstr,istatus)

        INQUIRE(FILE=TRIM(tmpstr), EXIST = fexist )
        IF(.NOT. fexist) THEN
          WRITE(6,'(1x,3a)') 'ERROR: The WRF file ',                    &
                         TRIM(tmpstr),' does not exist.'
          istatus = -1
          RETURN
        ELSE
          CALL get_wrf_patch_indices(TRIM(tmpstr),io_form,              &
                                  ips,ipe,ipss,ipse,jps,jpe,jpss,jpse,  &
                                  kps,kpe,kpss,kpse,nx,ny,istatus)
          IF (istatus /= 0) EXIT

          IF (.NOT. dset) THEN
            IF (npx == 1) THEN
              ids  = ips
              idss = ipss
            END IF

            IF (npx == nproc_x) THEN
              IF (paddingx) THEN
                imod = MOD(ips-ids,istride)
                IF (istride == 1 .OR. imod == 0) THEN
                  ide  = ips
                ELSE
                  ide = istride-imod + ips
                END IF
                idse = ipss-1
              ELSE
                ide  = ipe
                idse = ipse
              END IF
            END IF

            IF (npy == 1) THEN
              jds  = jps
              jdss = jpss
            END IF

            IF (npy == nproc_y) THEN
              IF (paddingy) THEN
                jmod = MOD(jps-jds,jstride)
                IF (jstride == 1 .OR. jmod == 0) THEN
                  jde  = jps
                ELSE
                  jde = jstride -jmod + jps
                END IF
                jdse = jpss-1
              ELSE
                jde  = jpe
                jdse = jpse
              END IF
            END IF

          END IF

          IF ( n > 1) THEN
            IF (in_a_row) THEN
              IF (jps /= jpssv .OR. jpe /= jpesv .OR. ips /= ipesv+1) THEN
                WRITE(6,'(/,1x,a,I4,2a,/,8x,2(a,I2),a,/,8x,a,/,8x,a,/)')  &
                  'ERROR: Patch ',n,' for file ',TRIM(tmpstr),            &
                  'at relative patch (',npx,',',npy,                      &
                  ') is not aligned in a row with its former patch.',     &
                  'Please check parameter nproc_xin. Either it was specified with a wrong number', &
                  'or the program has made a bad guess about it.'
                istatus = -2
                RETURN
              END IF
            ELSE
              IF (jps /= jpesv+1) THEN
                WRITE(6,'(/,1x,a,I4,2a,/,8x,2(a,I2),a,/,8x,a,/,8x,a,/)')  &
                  'ERROR: Patch ',n,' for file ',TRIM(tmpstr),            &
                  'at relative patch (',npx,',',npy,                      &
                  ') is not aligned in column with its former patch.',    &
                  'Please check parameter nproc_xin. Either it was specified with a wrong number', &
                  'or the program has made a bad guess about it.'
                istatus = -3
                RETURN
              END IF
            END IF
          END IF

          ipssv = ips
          ipesv = ipe
          jpssv = jps
          jpesv = jpe

          IF (dbglvl > 0) WRITE(6,'(3x,a,I2.2,a,I4,a,/,5x,a)')          &
             'WRF file ',nextdfil,': patch - ',n,' =', TRIM(tmpstr)
        END IF
      END DO
    END DO

    ifile = ifile + abstimei
    dset = .TRUE.

    WRITE(*,*)
  END DO


!-----------------------------------------------------------------------
!
! Validate nextdfil before return
!
!-----------------------------------------------------------------------

  IF(nextdfil < 1) THEN
    WRITE(6,'(a)') 'No input WRF file was valid. Please check the input file.'
    istatus = -3
    RETURN
  END IF

  IF (ide < ids .OR. jde < jds) THEN
    WRITE(6,'(1x,2(a,I4),/36x,2(a,I4),a)')                              &
    'ERROR: Domain indices are invalid: ids = ',ids,', ide = ',ide,     &
    '; jds = ',jds,', jde = ',jde,'.'
    istatus = -3
    RETURN
  END IF

  RETURN
END SUBROUTINE check_files_dimensions
!
!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE joinwrfncdf                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE  joinwrfncdf_mem(filenames,nfile,attadj,jointime,fileconv,   &
                  magnitude_processor,procs,npatch,splitfile,           &
                  ids,ide,idss,idse,jds,jde,jdss,jdse,istride,jstride,  &
                  idsinout,ideinout,idssinout,idseinout,                &
                  jdsinout,jdeinout,jdssinout,jdseinout,                &
                  kps,kpe,kpss,kpse,nproc_x_out,nproc_y_out,            &
                  outdirname,filetail,nvarout,varlists,debug,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Join WRF files in netCDF patches into one large piece.
!    This version emphasis io efficiency but use more memory.
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang (03/07/2013)
!
! MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: nfile
  LOGICAL, INTENT(IN)            :: attadj
  LOGICAL, INTENT(IN)            :: jointime, splitfile
  INTEGER, INTENT(IN)            :: fileconv, magnitude_processor
  INTEGER, INTENT(IN)            :: npatch
  INTEGER, INTENT(IN)            :: procs(npatch)
  INTEGER, INTENT(IN)            :: ids,ide,idss,idse,jds,jde,jdss,jdse
  INTEGER, INTENT(IN)            :: idsinout, ideinout, jdsinout, jdeinout
  INTEGER, INTENT(IN)            :: idssinout, idseinout, jdssinout, jdseinout
  INTEGER, INTENT(IN)            :: kps,kpe,kpss,kpse
  INTEGER, INTENT(IN)            :: nproc_x_out, nproc_y_out
  INTEGER, INTENT(IN)            :: istride, jstride
  INTEGER, INTENT(INOUT)         :: nvarout
  INTEGER, INTENT(IN)            :: debug
  INTEGER, INTENT(OUT)           :: istatus

  CHARACTER(LEN=*),  INTENT(IN)  :: filenames(nfile)
  CHARACTER(LEN=*),  INTENT(IN)  :: outdirname
  CHARACTER(LEN=5),  INTENT(IN)  :: filetail
  CHARACTER(LEN=40), INTENT(IN)  :: varlists(nvarout)

!
!-----------------------------------------------------------------------
!
! Including files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nf, nvar, n
  INTEGER :: strlen
  INTEGER :: ispatch
  LOGICAL :: patchset

  INTEGER :: idsout,   ideout,  jdsout, jdeout
  INTEGER :: idssout, idseout, jdssout, jdseout

  CHARACTER(LEN=256) :: infilename, outfilename
  INTEGER :: finid  !, foutid
  INTEGER, ALLOCATABLE :: finids(:), foutids(:)

  !LOGICAL :: paddingx, paddingy
  !INTEGER :: nxdim, nydim
  REAL    :: newctrlat, newctrlon, lat1, lat2,lon1,lon2
  INTEGER :: newlatidx, newlonidx, oldlatidx1, oldlonidx1, oldlatidx2, oldlonidx2
  LOGICAL :: latavg, lonavg
  INTEGER :: latset, lonset

  !
  ! Dimension variables
  !
  CHARACTER(LEN=32), PARAMETER :: xdimname  = 'west_east_stag'
  CHARACTER(LEN=32), PARAMETER :: ydimname  = 'south_north_stag'
  CHARACTER(LEN=32), PARAMETER :: xsdimname = 'west_east'
  CHARACTER(LEN=32), PARAMETER :: ysdimname = 'south_north'
  CHARACTER(LEN=32) :: diminnames(NF_MAX_DIMS)
  CHARACTER(LEN=32) :: dimname

  INTEGER :: nxid, nyid, nxlg, nylg, nxsid, nysid, nxslg, nyslg
  INTEGER :: narrsize, narrisizemax, narrasizemax
  INTEGER :: unlimdimid, unlimdimlen, unlimodimlen, odimid
  INTEGER :: ndims, dimid, dimlen

  INTEGER :: dimina(NF_MAX_DIMS)         ! Dimension size in original file
  !INTEGER :: dimouta(NF_MAX_DIMS)        ! Dimension size in joined files

  !
  ! Attribute variables
  !
  CHARACTER(LEN=32), PARAMETER :: attnm_ips  = 'WEST-EAST_PATCH_START_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_ipe  = 'WEST-EAST_PATCH_END_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_ipss = 'WEST-EAST_PATCH_START_UNSTAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_ipse = 'WEST-EAST_PATCH_END_UNSTAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jps  = 'SOUTH-NORTH_PATCH_START_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jpe  = 'SOUTH-NORTH_PATCH_END_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jpss = 'SOUTH-NORTH_PATCH_START_UNSTAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jpse = 'SOUTH-NORTH_PATCH_END_UNSTAG'
  CHARACTER(LEN=32) :: attname
  INTEGER :: ipsid,  ipeid,  jpsid,  jpeid
  INTEGER :: ipssid, ipseid, jpssid, jpseid
  INTEGER, ALLOCATABLE :: ips(:), ipe(:), ipss(:), ipse(:)
  INTEGER, ALLOCATABLE :: jps(:), jpe(:), jpss(:), jpse(:)

  INTEGER, ALLOCATABLE :: ipsout(:), ipeout(:), ipssout(:), ipseout(:)
  INTEGER, ALLOCATABLE :: jpsout(:), jpeout(:), jpssout(:), jpseout(:)

  INTEGER :: attnum, ngatts
  INTEGER :: istart, jstart, iend, jend, imod, jmod

  CHARACTER(LEN=32), PARAMETER :: attnm_anx = 'WEST-EAST_GRID_DIMENSION'
  CHARACTER(LEN=32), PARAMETER :: attnm_any = 'SOUTH-NORTH_GRID_DIMENSION'
  INTEGER :: anxid, anyid

  CHARACTER(LEN=32), PARAMETER :: attnm_adx = 'DX'
  CHARACTER(LEN=32), PARAMETER :: attnm_ady = 'DY'
  INTEGER :: adxid, adyid
  REAL    :: adx, ady

  INTEGER :: npatchout
  !
  ! Dataset varaibles
  !
  INTEGER, PARAMETER :: MAX_RANK = 4    ! Assume the max rank is 5
  CHARACTER(LEN=32) :: varname
  INTEGER :: varid, nvars, ovarid
  INTEGER :: vartype, varndims, varnatts
  INTEGER :: vardimids(MAX_RANK)
  INTEGER :: startidx(MAX_RANK), countidx(MAX_RANK), strideidx(MAX_RANK)
  INTEGER :: outstart(MAX_RANK), outstride(MAX_RANK)
  INTEGER :: vardim, vdimid

  INTEGER :: varidlists(NF_MAX_VARS), varoutidlists(NF_MAX_VARS)

  INTEGER, ALLOCATABLE :: varari(:)
  REAL,    ALLOCATABLE :: vararr(:)
  CHARACTER(LEN=256)   :: tmpstr,fmtstr

  !
  ! Added box and join-in-memeory feature
  !
  INTEGER :: nxout, nyout, nxyout, nzout
  INTEGER :: nxoutlg, nyoutlg, nxyoutlg
  INTEGER :: icout, jcout, kp

  REAL,    ALLOCATABLE, TARGET :: varlg(:)
  INTEGER, ALLOCATABLE, TARGET :: varlgi(:)

  REAL,    POINTER     :: varout(:)
  INTEGER, POINTER     :: varouti(:)

  INTEGER :: i,j,k,t
  INTEGER :: kin, kout
  INTEGER :: outsizes(MAX_RANK), outlgsizes(MAX_RANK)
  INTEGER :: insize1d, insize2d, insize3d, outsize1d, outsize2d, outsize3d

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  varidlists(:) = 0
  nxlg  = ide-ids+1
  nylg  = jde-jds+1
  nxslg = idse-idss+1
  nyslg = jdse-jdss+1

  nxoutlg  = (ideinout-idsinout)/istride+1
  nyoutlg  = (jdeinout-jdsinout)/jstride+1
  nxyoutlg = nxoutlg*nyoutlg
  nzout    = kpe-kps+1

  ALLOCATE(finids(npatch), STAT = istatus)
  ALLOCATE(ips (npatch),   STAT = istatus)
  ALLOCATE(ipe (npatch),   STAT = istatus)
  ALLOCATE(ipss(npatch),   STAT = istatus)
  ALLOCATE(ipse(npatch),   STAT = istatus)
  ALLOCATE(jps (npatch),   STAT = istatus)
  ALLOCATE(jpe (npatch),   STAT = istatus)
  ALLOCATE(jpss(npatch),   STAT = istatus)
  ALLOCATE(jpse(npatch),   STAT = istatus)

  ALLOCATE(varlg (nxyoutlg*nzout), STAT = istatus)
  ALLOCATE(varlgi(nxyoutlg*nzout), STAT = istatus)

  npatchout = nproc_x_out * nproc_y_out

  ALLOCATE(foutids(npatchout),   STAT = istatus)
  ALLOCATE(ipsout (npatchout),   STAT = istatus)
  ALLOCATE(ipeout (npatchout),   STAT = istatus)
  ALLOCATE(ipssout(npatchout),   STAT = istatus)
  ALLOCATE(ipseout(npatchout),   STAT = istatus)
  ALLOCATE(jpsout (npatchout),   STAT = istatus)
  ALLOCATE(jpeout (npatchout),   STAT = istatus)
  ALLOCATE(jpssout(npatchout),   STAT = istatus)
  ALLOCATE(jpseout(npatchout),   STAT = istatus)

  nxout = (nxoutlg-1)/nproc_x_out
  nyout = (nyoutlg-1)/nproc_y_out
  DO jcout = 1, nproc_y_out
    DO icout = 1, nproc_x_out
      kp = (jcout-1)*nproc_x_out+icout
      ipsout (kp) = nxout*(icout-1)+1
      ipeout (kp) = nxout*icout
      ipssout(kp) = ipsout(kp)
      ipseout(kp) = ipeout(kp)

      jpsout (kp) = nyout*(jcout-1)+1
      jpeout (kp) = nyout*jcout
      jpssout(kp) = jpsout(kp)
      jpseout(kp) = jpeout(kp)

      IF (icout == nproc_x_out) ipeout(kp) = ipeout(kp)+1
    END DO
    IF (jcout == nproc_y_out) jpeout(kp) = jpeout(kp)+1
  END DO

  nxout = nxout+1
  nyout = nyout+1
  IF (npatchout > 1) THEN
    nxyout = nxout*nyout
    ALLOCATE(varout (nxyout*nzout), STAT = istatus)
    ALLOCATE(varouti(nxyout*nzout), STAT = istatus)
  ELSE
    varout  => varlg
    varouti => varlgi
  END IF

  startidx(:)  = 1
  narrisizemax = 0
  narrasizemax = 0

  unlimdimlen  = 1

  istatus = 0
  DO nf = 1,nfile

    WRITE(*,'(1x,a,I3,3a,/)') '-------- ',nf,' : ',TRIM(filenames(nf)),' --------'

!-----------------------------------------------------------------------
!
! First, Create the merged file based on attributes from the first patch
!
!-----------------------------------------------------------------------

    IF (.NOT. jointime .OR. nf == 1) THEN

      strlen = LEN_TRIM(filenames(nf))
      n = INDEX(filenames(nf),'/',.TRUE.)

      IF ( .NOT. splitfile ) THEN
        WRITE(infilename, '(a)')       TRIM(filenames(nf))
      ELSE
        WRITE(fmtstr,'(a,2(I1,a))') '(a,a,I',magnitude_processor,'.',magnitude_processor,')'
        WRITE(infilename, FMT=fmtstr) TRIM(filenames(nf)),'_',procs(1)
      END IF

      IF (debug > 0) WRITE(6,'(1x,2a)') 'Opening file - ',TRIM(infilename)
      istatus = nf_open(infilename,NF_NOWRITE,finid)   ! Open file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      WRITE(outfilename,'(3a)') TRIM(outdirname),                       &
                              filenames(nf)(n+1:strlen),filetail
      DO kp = 1, npatchout
        IF (npatchout > 1) THEN
          WRITE(tmpstr,FMT=fmtstr) TRIM(outfilename),'_',kp-1
        ELSE
          tmpstr = outfilename
        END IF

        IF (debug > 0) WRITE(6,'(1x,2a)') 'Creating file - ',TRIM(tmpstr)
        !istatus = nf_create(TRIM(outfilename),NF_CLOBBER,foutid)                     ! CDF 1
        istatus = NF_CREATE(TRIM(tmpstr),IOR(NF_CLOBBER,NF_64BIT_OFFSET),foutids(kp)) ! CDF 2
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO

      !
      ! Set dimensions
      !
      IF (fileconv < 2) THEN
        istatus = nf_inq_dimid(finid,xdimname,nxid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        istatus = nf_inq_dimid(finid,ydimname,nyid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      ELSE               ! NMM core
        nxid = -1
        nyid = -1
      END IF

      istatus = nf_inq_dimid(finid,xsdimname,nxsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_dimid(finid,ysdimname,nysid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_unlimdim(finid,unlimdimid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_ndims(finid,ndims)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Copying dimensions - ',ndims
      DO dimid = 1,ndims
        istatus = nf_inq_dim(finid,dimid,dimname,dimlen)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        diminnames(dimid) = dimname
        dimina(dimid)     = dimlen           ! Save dimension id and len
        !dimouta(dimid) = dimlen             ! Output dimension id and len
        DO kp = 1, npatchout
          IF (dimid == nxid) THEN
            dimlen = ipeout(kp)-ipsout(kp)+1
            !dimouta(dimid) = dimlen
          ELSE IF (dimid == nxsid) THEN
            dimlen = (ipseout(kp)-ipssout(kp))+1
            !dimouta(dimid) = dimlen
          ELSE IF (dimid == nyid) THEN
            dimlen = jpeout(kp)-jpsout(kp)+1
            !dimouta(dimid) = dimlen
          ELSE IF (dimid == nysid) THEN
            dimlen = (jpseout(kp)-jpssout(kp))+1
            !dimouta(dimid) = dimlen
          ELSE IF (dimid == unlimdimid) THEN
            dimlen = NF_UNLIMITED
          END IF

          IF (debug > 0) WRITE(6,'(9x,a,a20,a,I4,a,I3)')                &
              'Dimension name - ',dimname, ', length = ',dimlen,        &
              ' ==> output patch ',kp
          istatus = nf_def_dim(foutids(kp),dimname,dimlen,odimid)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        END DO
      END DO

      !
      ! Set Global attributes
      !
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ips),ipsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ipe),ipeid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ipss),ipssid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ipse),ipseid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jps),jpsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jpe),jpeid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jpss),jpssid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jpse),jpseid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_anx),anxid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_any),anyid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_adx),adxid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ady),adyid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_real(finid,NF_GLOBAL,TRIM(attnm_adx),adx)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_real(finid,NF_GLOBAL,TRIM(attnm_ady),ady)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_natts(finid,ngatts)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (attadj) THEN
        idsout  = 1
        ideout  = nxoutlg
        idssout = 1
        idseout = (idseinout-idssinout)/istride+1

        jdsout  = 1
        jdeout  = nyoutlg
        jdssout = 1
        jdseout = (jdseinout-jdssinout)/jstride+1

        newlatidx = (jdeout-jdsout)/2    ! Added to modify CEN_LAT & CEN_LON
        newlonidx = (ideout-idsout)/2    ! 0 base

        oldlatidx1 = newlatidx*jstride+jdsinout  ! old domain index
        oldlonidx1 = newlonidx*istride+idsinout  ! old domain index

        latavg = .FALSE.
        lonavg = .FALSE.
        latset = 0
        lonset = 0

        IF (MOD((jdeout-jdsout),2) /= 0) latavg = .TRUE.
        IF (MOD((ideout-idsout),2) /= 0) lonavg = .TRUE.

        oldlatidx2 = -99999999
        oldlonidx2 = -99999999
        IF (.NOT. latavg .AND. .NOT. lonavg) THEN
          oldlatidx2 = oldlatidx1-jstride
          oldlonidx2 = oldlonidx1-istride
        END IF

        IF (debug > 0) WRITE(6,'(1x,a,/,1x,a,2L,a,2I4,/,21x,a,2I4,a))') &
          'It is required to adjust global attribute.',                 &
          'latavg,lonavg = ',latavg,lonavg,                             &
          ', lonidx,latidx = ',newlonidx+1,newlatidx+1,                 &
          ', lonidx,latidx = ',oldlonidx1,oldlatidx1,' (in original domain)'

      ELSE
        idsout  = idsinout
        ideout  = ideinout
        idssout = idssinout
        idseout = idseinout

        jdsout  = jdsinout
        jdeout  = jdeinout
        jdssout = jdssinout
        jdseout = jdseinout
      END IF

      IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Copying global attributes - ',ngatts

      DO attnum = 1,ngatts

        istatus = nf_inq_attname(finid,NF_GLOBAL,attnum,attname)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        IF (debug > 0) WRITE(6,'(9x,2a)') 'Attribute name - ',TRIM(attname)

        DO kp = 1, npatchout

          IF (attnum == ipsid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_ips),NF_INT,1,ipsout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == ipeid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_ipe),NF_INT,1,ipeout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == jpsid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_jps),NF_INT,1,jpsout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == jpeid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_jpe),NF_INT,1,jpeout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == ipssid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_ipss),NF_INT,1,ipssout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == ipseid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_ipse),NF_INT,1,ipseout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == jpssid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_jpss),NF_INT,1,jpssout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == jpseid) THEN
            istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_jpse),NF_INT,1,jpseout(kp))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == anxid) THEN      ! Adjust nx
            IF (istride > 1 .OR. attadj) THEN
              istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_anx),NF_INT,1,nxoutlg)
              IF (istatus /= NF_NOERR) CALL handle_err(istatus)
            ELSE
              istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutids(kp),NF_GLOBAL)
              IF (istatus /= NF_NOERR) CALL handle_err(istatus)
            END IF
          ELSE IF (attnum == anyid) THEN      ! Adjust ny
            IF (jstride > 1 .OR. attadj) THEN
              istatus = NF_PUT_ATT_INT(foutids(kp),NF_GLOBAL,TRIM(attnm_any),NF_INT,1,nyoutlg)
              IF (istatus /= NF_NOERR) CALL handle_err(istatus)
            ELSE
              istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutids(kp),NF_GLOBAL)
              IF (istatus /= NF_NOERR) CALL handle_err(istatus)
            END IF
          ELSE IF (attnum == adxid) THEN      ! adjust dx
            istatus = NF_PUT_ATT_REAL(foutids(kp),NF_GLOBAL,TRIM(attnm_adx),NF_REAL,1,adx*istride)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE IF (attnum == adyid) THEN      ! adjust dy
            istatus = NF_PUT_ATT_REAL(foutids(kp),NF_GLOBAL,TRIM(attnm_ady),NF_REAL,1,ady*jstride)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE
            istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutids(kp),NF_GLOBAL)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          END IF

        END DO

      END DO

      !
      ! Define variables
      !
      istatus = nf_inq_nvars(finid,nvars)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (nvarout >= nvars) THEN
        nvarout = nvars
        DO n = 1,nvars
          varidlists(n) = n
        END DO
      ELSE
        nvar = nvarout         ! suppose to process this number
        nvarout = 0            ! actually got
        DO n = 1,nvar
          istatus = nf_inq_varid(finid,TRIM(varlists(n)),ovarid)
          IF (istatus /= NF_NOERR) THEN
            WRITE(6,'(1x,3a)') 'WARNING: Variable ',TRIM(varlists(n)),' not found. Skipped.'
            CYCLE
          END IF
          nvarout = nvarout + 1
          varidlists(nvarout) = ovarid
        END DO
      END IF

      IF (debug > 0) WRITE(6,'(5x,a,I4)') 'Defining variables - ',nvarout

      DO n = 1,nvarout
        varid = varidlists(n)
        istatus = nf_inq_var(finid,varid,varname,vartype,varndims,vardimids,varnatts)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        IF (debug > 0) WRITE(6,'(9x,2a)') 'Variables - ',TRIM(varname)

        DO kp = 1, npatchout

          ! Dimensions should be in the same order
          istatus = nf_def_var(foutids(kp),varname,vartype,varndims,vardimids,ovarid)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          varoutidlists(n) = ovarid

          DO attnum = 1,varnatts          ! Copy variable attributes
            istatus = nf_inq_attname(finid,varid,attnum,attname)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            istatus = nf_copy_att(finid,varid,attname,foutids(kp),ovarid)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          END DO

        END DO


      END DO

      DO kp = 1, npatchout
        istatus = nf_enddef(foutids(kp))
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO

      IF(debug > 0) WRITE(6,'(1x,a)') 'Merged file have been defined.'

      istatus = nf_close(finid)                              ! Close file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    END IF         ! File created.

    IF (.NOT. jointime) unlimdimlen = 1
    unlimodimlen = 0

!-----------------------------------------------------------------------
!
! Open each patch and get their dimensions
!
!-----------------------------------------------------------------------

    DO n = 1,npatch

      IF ( .NOT. splitfile ) THEN
        WRITE(infilename, '(a)')       TRIM(filenames(nf))
      ELSE
        WRITE(fmtstr,'(a,2(I1,a))') '(a,a,I',magnitude_processor,'.',magnitude_processor,')'
        WRITE(infilename, FMT=fmtstr) TRIM(filenames(nf)),'_',procs(n)
      END IF

      IF (debug > 0) WRITE(6,'(1x,2a)') 'Opening file - ',TRIM(infilename)

      istatus = nf_open(TRIM(infilename),NF_NOWRITE,finids(n))   ! Open file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      !
      ! Get patch indice
      !
      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_ips),ips(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_ipe),ipe(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_ipss),ipss(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_ipse),ipse(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_jps),jps(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_jpe),jpe(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_jpss),jpss(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finids(n),NF_GLOBAL,TRIM(attnm_jpse),jpse(n))
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    END DO

!-----------------------------------------------------------------------
!
! loop over each variable
!
!-----------------------------------------------------------------------

    var_loop: DO nvar = 1,nvarout

      varid = varidlists(nvar)

      ispatch       = 0            ! Wheter this variable is in patches
      patchset      = .FALSE.
      outstart(:)   = 1
      outlgsizes(:) = 1

      patch_loop: DO n = 1,npatch

        vardimids(:) = 0
        istatus = nf_inq_var(finids(n),varid,varname,vartype,varndims,vardimids,varnatts)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Get and save dimension size for this patch
!
!-----------------------------------------------------------------------

        IF (ideinout  < ips(n)  .OR. idsinout  > ipe(n)   .OR.          &
            jdeinout  < jps(n)  .OR. jdsinout  > jpe(n) ) THEN
          IF (debug > 1) WRITE(*,'(1x,a,i4,a,2(/,1x,a,4(a,I4),/,8x,4(a,I4)))')                 &
          'Skip patch ',procs(n),' because it is outside of the output box.',                  &
          'INPUT: ','ips = ',ips(n),', ipe = ',ipe(n),', ipss = ',ipss(n),', ipse = ',ipse(n), &
                    'jps = ',jps(n),', jpe = ',jpe(n),', ipss = ',jpss(n),', jpse = ',jpse(n), &
          'OUTBX: ','ids = ',idsinout,', ide = ',ideinout,', idss = ',idssinout,', idse = ',idseinout,     &
                    'jds = ',jdsinout,', jde = ',jdeinout,', idse = ',jdssinout,', jdse = ',jdseinout

          CYCLE
        END IF

        istatus = nf_inq_ndims(finids(n),ndims)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        dimina(:) = 0
        DO dimid = 1,ndims
          istatus = nf_inq_dim(finids(n),dimid,dimname,dimlen)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          diminnames(dimid) = dimname
          dimina(dimid)     = dimlen             ! Save dimension id and len
        END DO

!-----------------------------------------------------------------------
!
! Handle dimensions
!
!-----------------------------------------------------------------------

        startidx(:)  = 1
        countidx(:)  = 1
        strideidx(:) = 1
        narrsize = 1

        istart = 1  ! in global index of the input domain
        iend   = 1
        jstart = 1
        jend   = 1
        DO vardim = 1, varndims
          vdimid = vardimids(vardim)
          countidx(vardim) = dimina (vdimid)
          IF ( vdimid == nxid) THEN
            IF ( ips(n) > idsinout) THEN
              imod = MOD((ips(n)-idsinout),istride)
              IF (imod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = istride - imod + 1
              END IF
              istart = ips(n)+startidx(vardim)-1    ! start index relative to global domain
            ELSE
              startidx(vardim) = idsinout-ips(n)+1
              istart = idsinout
            END IF

            iend   = MIN(ipe(n),ideinout)
            IF (iend < istart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (iend-istart)/istride + 1
            END IF
            strideidx(vardim) = istride

            !outstart(vardim) = (istart-idsinout)/istride + 1
            IF (.NOT. patchset) THEN
              ispatch = ispatch+1
              outlgsizes(vardim) = nxoutlg
            END IF

          ELSE IF ( vdimid == nyid) THEN
            IF (jps(n) > jdsinout) THEN
              jmod = MOD((jps(n)-jdsinout),jstride)
              IF (jmod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = jstride - jmod + 1
              END IF
              jstart = jps(n)+startidx(vardim)-1
            ELSE
              startidx(vardim) = jdsinout-jps(n)+1
              jstart = jdsinout
            END IF

            jend   = MIN(jpe(n),jdeinout)
            IF (jend < jstart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (jend-jstart)/jstride + 1
            END IF
            strideidx(vardim) = jstride

            !outstart(vardim) = (jstart - jdsinout)/jstride + 1
            IF (.NOT. patchset) THEN
              ispatch = ispatch+10
              outlgsizes(vardim) = nyoutlg
            END IF

          ELSE IF ( vdimid == nxsid) THEN

            IF (ipss(n) > idssinout) THEN
              imod = MOD((ipss(n)-idssinout),istride)
              IF (imod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = istride - imod + 1
              END IF
              istart = ipss(n)+startidx(vardim)-1
            ELSE
              startidx(vardim) = idssinout-ipss(n)+1
              istart = idssinout
            END IF

            iend   = MIN(ipse(n),idseinout)
            IF (iend < istart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (iend-istart)/istride + 1
            END IF
            strideidx(vardim) = istride

            !outstart(vardim) = (istart - idss)/istride + 1
            IF (.NOT. patchset) THEN
              ispatch = ispatch+2
              outlgsizes(vardim) = idseout-idssout+1
            END IF
          ELSE IF ( vdimid == nysid) THEN
            IF (jpss(n) > jdssinout) THEN
              jmod = MOD((jpss(n)-jdssinout),jstride)
              IF (jmod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = jstride - jmod + 1
              END IF
              jstart = jpss(n)+startidx(vardim)-1
            ELSE
              startidx(vardim) = jdssinout-jpss(n)+1
              jstart = jdssinout
            END IF

            jend   = MIN(jpse(n),jdseinout)
            IF (jend < jstart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (jend-jstart)/jstride + 1
            END IF
            strideidx(vardim) = jstride

            !outstart(vardim) = (jstart - jdss)/jstride + 1
            IF (.NOT. patchset) THEN
              ispatch = ispatch+20
              outlgsizes(vardim) = jdseout-jdssout+1
            END IF

          ELSE IF (vdimid == unlimdimid) THEN
            outstart(vardim) = unlimdimlen
            IF (unlimodimlen <= 0) THEN
              unlimodimlen = countidx(vardim)
            ELSE
              IF ( unlimodimlen /= countidx(vardim)) THEN
                WRITE(6,'(1x,a,/)') 'ERROR: Inconsisten size for UNLIMITED dimension.'
                istatus = -1
                RETURN
              END IF
            END IF
            outlgsizes(vardim) = unlimodimlen

          ELSE
            outstart(vardim) = 1
            outlgsizes(vardim) = countidx(vardim)
          END IF

          IF (countidx(vardim) <= 0) THEN
            IF (debug > 0) THEN
              !WRITE(6,'(9x,2a)') 'Processing variables - ',TRIM(varname)
              WRITE(6,'(12x,a,i0,3a)')                                  &
                'Patch ',procs(n),' skipped because dimension "',       &
                TRIM(diminnames(vdimid)),'" has zero length.'
            END IF
            CYCLE patch_loop
          END IF

          narrsize = countidx(vardim)*narrsize
        END DO

        !IF ( narrsize <= 0 ) THEN
        !  IF (debug > 1) WRITE(*,'(1x,a,i4,3a)') 'Skip patch ',procs(n),&
        !         ' for variable ',varname,' because it is non-stagered.'
        !  CYCLE
        !END IF

        IF (.NOT. patchset) THEN
          IF (debug > 2) WRITE(6,'(/,5x,2a,2(a,I2))')                   &
            'Processing variables - ',TRIM(varname),' with rank = ',varndims,', ispatch = ',ispatch

          patchset = .TRUE.
        END IF

        IF ( n > 1 .AND. ( ispatch == 0) ) THEN
          IF (debug > 2) WRITE(6,'(9x,a,i3,a)') 'Skip patch - ',procs(n),' and up.'
          EXIT patch_loop
        ELSE
          IF (debug > 2) THEN
            WRITE(6,'(9x,3(a,I4))') 'Dimensions in Patch : ',procs(n),', istart = ',istart, ', jstart = ',jstart

            DO vardim = 1,varndims
              vdimid = vardimids(vardim)
              WRITE(6,'(12x,a12,5(a,I4))') TRIM(diminnames(vdimid)),             &
              ', outstart = ',outstart(vardim),', size = ', countidx(vardim),    &
              ' <-- start = ',startidx(vardim),', stride = ',strideidx(vardim),  &
              ', size = ',outlgsizes(vardim)
            END DO
          END IF
        END IF

        ! do not have to merge, use values from the first file

!        IF (.NOT. ispatch(nvar)) THEN
!
!          IF (debug > 0) WRITE(6,'(9x,2a)') 'Copying variables - ',TRIM(varname)
!
!          istatus = NF_COPY_VAR(finid,varid,foutid)
!          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
!
!        ELSE
          ovarid = varoutidlists(nvar)

          IF (debug == 1) WRITE(6,'(9x,3a,I4)')                         &
             'Fetching variables - ',TRIM(varname),' from patch: ',procs(n)

          outsize1d  = outlgsizes(1)
          outsize2d  = outlgsizes(1)*outlgsizes(2)
          outsize3d  = outlgsizes(1)*outlgsizes(2)*outlgsizes(3)

          insize1d = countidx(1)
          insize2d = countidx(1)*countidx(2)
          insize3d = countidx(1)*countidx(2)*countidx(3)

          SELECT CASE (vartype)

!-----------------------------------------------------------------------
!
! Integers
!
!-----------------------------------------------------------------------

          CASE (NF_INT)

            IF (narrsize > narrisizemax) THEN   ! Allocate input array only when necessary
              IF (ALLOCATED(varari)) DEALLOCATE(varari, STAT = istatus)
              ALLOCATE(varari(narrsize), STAT = istatus)
              narrisizemax = narrsize
            END IF

            istatus = NF_GET_VARS_INT(finids(n),varid,startidx,countidx,strideidx,varari)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            IF (ispatch > 0) THEN
              DO t = 1, countidx(4)
                DO k = 1, countidx(3)
                  DO j = 1, countidx(2)
                    DO i = 1, countidx(1)
                      !IF ( ispatch == 1) THEN   ! join in x direction only
                      !  kout = istart-idsinout+i+(j-1)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !ELSE IF (ispatch == 2) THEN ! join in y direction only
                      !  kout = jstart-jdsinout+i+(j-1)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !ELSE IF (ispatch == 3) THEN ! join in both x direction and y direction
                        kout = istart-idsinout+i+(jstart-jdsinout+j-1)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !ELSE
                      !  WRITE(*,'(1x,a,i0,a,i0)') 'ERROR: ispatch = ',ispatch,', n= ',n
                      !  STOP
                      !END IF
                      kin = i+(j-1)*insize1d+(k-1)*insize2d+(t-1)*insize3d
                      varlgi(kout) = varari(kin)
                    END DO
                  END DO
                END DO
              END DO

            ELSE
              varlgi(1:narrsize) = varari(1:narrsize)
            END IF

            !istatus = nf_put_vara_INT(foutid,ovarid,outstart,countidx,varari)
            !IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Reals
!
!-----------------------------------------------------------------------

          CASE (NF_FLOAT)

            IF (narrsize > narrasizemax) THEN   ! Allocate input array only when necessary
              IF (ALLOCATED(vararr)) DEALLOCATE(vararr, STAT = istatus)
              ALLOCATE(vararr(narrsize), STAT = istatus)
              narrasizemax = narrsize
            END IF

            istatus = NF_GET_VARS_REAL(finids(n),varid,startidx,countidx,strideidx,vararr)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            IF (ispatch > 0) THEN
              DO t = 1, countidx(4)
                DO k = 1, countidx(3)
                  DO j = 1, countidx(2)
                    DO i = 1, countidx(1)
                      !IF ( ispatch == 1) THEN     ! join in x direction only
                      !  kout = istart-idsinout+i+(j-1)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !ELSE IF (ispatch == 2) THEN ! join in y direction only
                      !  kout = jstart-jdsinout+i+(j-1)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !ELSE IF (ispatch == 3) THEN ! join in both x direction and y direction
                        kout = istart-idsinout+i+(jstart-jdsinout+j-1)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !ELSE
                      !  WRITE(*,'(1x,a,i0,a,i0)') 'ERROR: ispatch = ',ispatch,', n= ',n
                      !  STOP
                      !END IF
                      kin = i+(j-1)*insize1d+(k-1)*insize2d+(t-1)*insize3d
                      !WRITE(*,'(1x,i2,a,4I3,a,I8,a,i8)') n,' : ',i,j,k,t,' - ',kin,' -> ',kout
                      varlg(kout) = vararr(kin)
                    END DO
                  END DO
                END DO
              END DO
            ELSE
              varlg(1:narrsize) = vararr(1:narrsize)
            END IF

            !istatus = nf_put_vara_REAL(foutid,ovarid,outstart,countidx,vararr)
            !IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Character string
!
!-----------------------------------------------------------------------

          CASE (NF_CHAR)

            istatus = NF_GET_VARS_TEXT(finids(n),varid,startidx,countidx,strideidx,tmpstr)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            !WRITE(6,'(3x,5a,I3)') 'Get   ',TRIM(varname),' as ',tmpstr(1:countidx(1)),' from patch ',procs(n)

            !istatus = nf_put_vara_TEXT(foutid,ovarid,outstart,countidx,TRIM(tmpstr))
            !IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          CASE DEFAULT
            WRITE(6,'(1x,a,I2)') 'ERROR: unsupported variable type = ',vartype
            istatus = -4
            RETURN
          END SELECT
!        END IF  ! ispatch(nvar)

        IF (nf == 1 .AND. attadj) THEN
          !
          ! Added to find the new CEN_LON & CEN_LAT
          !
          IF ( (oldlatidx1 >= jps(n) .AND. oldlatidx1 <= jpe(n)) .AND.  &
               (oldlonidx1 >= ips(n) .AND. oldlonidx1 <= ipe(n)) ) THEN

            IF ( .NOT. lonavg .AND. latavg ) THEN      ! at U point

              IF (TRIM(varname) == 'XLONG_U' .OR. TRIM(varname) == 'XLAT_U' ) THEN
                newlatidx = (oldlatidx1-jstart)/jstride     ! 0-base
                newlonidx = (oldlonidx1-istart)/istride+1   ! 1-base
                strlen = newlatidx*countidx(1)+newlonidx

                IF (TRIM(varname) == 'XLONG_U' ) THEN
                  newctrlon = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                     '*new cen_lon at U point. lonidx = ',newlonidx,      &
                     ', newctrlon = ',newctrlon,' from patch ', procs(n)
                  lonset = lonset+1
                ELSE IF (TRIM(varname) == 'XLAT_U') THEN
                  newctrlat = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                       '*new cen_lat at U point. latidx = ',newlatidx,    &
                       ', newctrlat = ',newctrlat,' from patch ', procs(n)
                  latset = latset+1
                END IF
              END IF

            END IF

            IF ( lonavg .AND. .NOT. latavg ) THEN      ! at V point
              IF (TRIM(varname) == 'XLONG_V' .OR. TRIM(varname) == 'XLAT_V'  ) THEN
                newlatidx = (oldlatidx1-jstart)/jstride
                newlonidx = (oldlonidx1-istart)/istride+1
                strlen = newlatidx*countidx(1)+newlonidx

                IF (TRIM(varname) == 'XLONG_V' ) THEN
                  newctrlon = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                     '*new cen_lon at V point. lonidx = ',newlonidx,      &
                     ', newctrlon = ',newctrlon,' from patch ', procs(n)

                  lonset = lonset+1
                ELSE IF (TRIM(varname) == 'XLAT_V') THEN
                  newctrlat = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                     '*new cen_lat at V point. latidx = ',newlatidx+1,    &
                     ', newctrlat = ',newctrlat,' from patch ', procs(n)
                  latset = latset+1
                END IF

              END IF
            END IF

            IF ( lonavg .AND. latavg ) THEN            ! at M point
              IF (TRIM(varname) == 'XLONG' .OR. TRIM(varname) == 'XLAT'  ) THEN
                newlatidx = (oldlatidx1-jstart)/jstride
                newlonidx = (oldlonidx1-istart)/istride+1
                strlen = newlatidx*countidx(1)+newlonidx
                IF (TRIM(varname) == 'XLONG' ) THEN
                  newctrlon = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                     '*new cen_lon at M point. lonidx = ',newlonidx,      &
                     ', newctrlon = ',newctrlon,' from patch ',procs(n)
                  lonset = lonset+1
                ELSE IF (TRIM(varname) == 'XLAT' ) THEN
                  newctrlat = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                     '*new cen_lat at M point. latidx = ',newlatidx,      &
                     ', newctrlat = ',newctrlat,' from patch ', procs(n)
                  latset = latset+1
                END IF
              END IF

            END IF

            IF ( .NOT. lonavg .AND. .NOT. latavg ) THEN ! at grid point
              IF (TRIM(varname) == 'XLONG_U' .OR. TRIM(varname) == 'XLAT_V'  ) THEN
                newlatidx = (oldlatidx1-jstart)/jstride
                newlonidx = (oldlonidx1-istart)/istride+1
                strlen = newlatidx*countidx(1)+newlonidx
                IF (TRIM(varname) == 'XLONG_U' ) THEN
                  lon1 = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')        &
                      '*set lon1 = ',lon1,' at = ',newlatidx+1,newlonidx, &
                      ' from patch ', procs(n)
                ELSE IF (TRIM(varname) == 'XLAT_V' ) THEN
                  lat1 = vararr(strlen)
                  IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')        &
                      '*set lat1 = ',lat1,' at = ',newlatidx+1,newlonidx, &
                      ' from patch ', procs(n)
                END IF
              END IF

            END IF

          END IF

          IF ( (oldlatidx1 >= jps(n) .AND. oldlatidx1 <= jpe(n)) .AND.  &
               (oldlonidx2 >= ips(n) .AND. oldlonidx2 <= ipe(n)) ) THEN

            IF ( .NOT. lonavg .AND. .NOT. latavg ) THEN ! at grid point

              IF (TRIM(varname) == 'XLAT_V') THEN
                newlatidx = (oldlatidx1-jstart)/jstride
                newlonidx = (oldlonidx2-istart)/istride+1
                strlen = newlatidx*countidx(1)+newlonidx
                lat2 = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')          &
                    '*set lat2 = ',lat2,' at = ',newlatidx+1,newlonidx,   &
                    ' from patch ', procs(n)
              END IF

            END IF

          END IF

          IF ( (oldlatidx2 >= jps(n) .AND. oldlatidx2 <= jpe(n)) .AND.  &
               (oldlonidx1 >= ips(n) .AND. oldlonidx1 <= ipe(n)) ) THEN

            IF ( .NOT. lonavg .AND. .NOT. latavg ) THEN ! at grid point

              IF (TRIM(varname) == 'XLONG_U') THEN
                newlatidx = (oldlatidx2-jstart)/jstride
                newlonidx = (oldlonidx1-istart)/istride+1
                strlen = newlatidx*countidx(1)+newlonidx
                lon2 = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')          &
                   '*set lon2 = ',lon2,' at = ',newlatidx+1,newlonidx,    &
                   ' from patch ', procs(n)
              END IF

            END IF

          END IF

        END IF

      END DO patch_loop

!-----------------------------------------------------------------------
!
! Write out the joined Dataset
!
!-----------------------------------------------------------------------

      IF (debug > 0) WRITE(6,*)

      insize1d = outlgsizes(1)
      insize2d = outlgsizes(1)*outlgsizes(2)
      insize3d = outlgsizes(1)*outlgsizes(2)*outlgsizes(3)

      DO kp = 1, npatchout
        outsizes(:) = outlgsizes(:)
        istart = 1
        iend   = outlgsizes(1)
        jstart = 1
        jend   = outlgsizes(2)

        IF (MOD(ispatch,10) == 1) THEN
          istart = ipsout(kp)
          iend   = ipeout(kp)
          outsizes(1) = iend-istart+1
        ELSE IF (MOD(ispatch,10) == 2) THEN
          istart = ipssout(kp)
          iend   = ipseout(kp)
          outsizes(1) = iend-istart+1
        END IF

        IF ( ispatch/10 == 1) THEN
          IF (MOD(ispatch,10) /= 0) THEN
            jstart = jpsout(kp)
            jend   = jpeout(kp)
            outsizes(2) = jend-jstart+1
          ELSE
            istart = jpsout(kp)
            iend   = jpeout(kp)
            outsizes(1) = jend-jstart+1
          END IF
        ELSE IF( ispatch/10 == 2) THEN
          IF (MOD(ispatch,10) /= 0) THEN
            jstart = jpssout(kp)
            jend   = jpseout(kp)
            outsizes(2) = jend-jstart+1
          ELSE
            istart = jpssout(kp)
            iend   = jpseout(kp)
            outsizes(1) = jend-jstart+1
          END IF
        END IF

        IF (debug > 0) THEN
          WRITE(6,'(9x,3a,i2.2)') 'Writing variables - ',TRIM(varname),' to output patch : ',kp
          DO vardim = 1,4
            WRITE(6,'(12x,a,i2,2(a,i4))') 'Dimension ',vardim,          &
                                    ' output size = ',outsizes(vardim), &
                                   ', global size = ',outlgsizes(vardim)
          END DO
        END IF

        outsize1d  = outsizes(1)
        outsize2d  = outsizes(1)*outsizes(2)
        outsize3d  = outsizes(1)*outsizes(2)*outsizes(3)

        SELECT CASE (vartype)

          CASE (NF_INT)
            IF (npatchout > 1) THEN
              DO t = 1, outsizes(4)
                DO k = 1, outsizes(3)
                  DO j = jstart, jend
                    DO i = istart, iend
                      kin = i+(j-1)*insize1d+(k-1)*insize2d+(t-1)*insize3d
                      kout = i-istart+1+(j-jstart)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !WRITE(*,'(1x,i2,a,4I3,a,I8,a,i8)') kp,' : ',i,j,k,t,' - ',kin,' -> ',kout
                      varouti(kout) = varlgi(kin)
                    END DO
                  END DO
                END DO
              END DO
            ELSE
              ! do nothing
            END IF
            istatus = nf_put_vara_INT(foutids(kp),ovarid,outstart,outsizes,varouti)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          CASE (NF_FLOAT)

            IF (npatchout > 1) THEN
              DO t = 1, outsizes(4)
                DO k = 1, outsizes(3)
                  DO j = jstart, jend
                    DO i = istart, iend
                      kin = i+(j-1)*insize1d+(k-1)*insize2d+(t-1)*insize3d
                      kout = i-istart+1+(j-jstart)*outsize1d+(k-1)*outsize2d+(t-1)*outsize3d
                      !WRITE(*,'(1x,i2,a,4I3,a,I8,a,i8)') kp,' : ',i,j,k,t,' - ',kin,' -> ',kout
                      varout(kout) = varlg(kin)
                    END DO
                  END DO
                END DO
              END DO
            ELSE
              ! do nothing
            END IF
            istatus = nf_put_vara_REAL(foutids(kp),ovarid,outstart,outsizes,varout)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          CASE (NF_CHAR)

            !WRITE(6,'(3x,5a)') 'Write ',TRIM(varname),' as ',tmpstr(1:outsizes(1))
            istatus = nf_put_vara_TEXT(foutids(kp),ovarid,outstart,outsizes,TRIM(tmpstr))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          CASE DEFAULT
            WRITE(6,'(1x,a,I2)') 'ERROR: unsupported variable type = ',vartype
            istatus = -4
            RETURN
        END SELECT

      END DO ! out patch loop

    END DO var_loop

    unlimdimlen = unlimdimlen + unlimodimlen                 ! Add # of time levels
                                                             ! in output file
!-----------------------------------------------------------------------
!
! Close the input file handlers
!
!-----------------------------------------------------------------------

    DO n = 1, npatch
      istatus = nf_close(finids(n))                              ! Close file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    END DO

!-----------------------------------------------------------------------
!
! Close the output file if applicable
!
!-----------------------------------------------------------------------

    IF (.NOT. jointime .OR. nf == nfile) THEN

      IF (attadj) THEN

        IF (nf == 1) THEN

          IF ( .NOT. lonavg .AND. .NOT. latavg .AND.                    &
               latset == 0 .AND. lonset == 0) THEN ! at grid point
            newctrlat = 0.5*(lat1+lat2)
            newctrlon = 0.5*(lon1+lon2)
            latset = latset+1
            lonset = lonset+1
          END IF

          IF (latset /= 1 .OR. lonset /= 1) THEN
            WRITE(*,'(1x,2(a,i0),8x,a)')                                &
              'ERROR: latset = ',latset,', lonset = ',lonset,           &
              'Program aborting ...     Please report.'
            STOP
          END IF
          !write(*,*) 'lat/lon = ',lat1,lat2,lon1,lon2,latset,lonset

        END IF

        WRITE(6,'(/,3x,2(a,F8.2),/)')                                   &
             'Changing CEN_LAT & CEN_LON to ',newctrlat,', ',newctrlon

        DO kp = 1, npatchout
          istatus = NF_REDEF( foutids(kp) )
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          istatus = NF_PUT_ATT_REAL(foutids(kp),NF_GLOBAL,'CEN_LAT',NF_REAL,1,newctrlat)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          istatus = NF_PUT_ATT_REAL(foutids(kp),NF_GLOBAL,'MOAD_CEN_LAT',NF_REAL,1,newctrlat)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          istatus = NF_PUT_ATT_REAL(foutids(kp),NF_GLOBAL,'CEN_LON',NF_REAL,1,newctrlon)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          istatus = NF_ENDDEF(foutids(kp))
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        END DO
      END IF

      DO kp = 1, npatchout
        istatus = nf_close(foutids(kp))
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO
    END IF

  END DO

  DEALLOCATE(finids, foutids)
  DEALLOCATE(ips,ipe, jps, jpe, ipss, ipse, jpss, jpse)
  DEALLOCATE(ipsout,ipeout, jpsout, jpeout, ipssout, ipseout, jpssout, jpseout)
  IF (npatchout > 1) DEALLOCATE(varout, varouti)
  DEALLOCATE(varlg,  varlgi)
  IF (ALLOCATED(vararr)) DEALLOCATE(vararr)
  IF (ALLOCATED(varari)) DEALLOCATE(varari)

  RETURN
END SUBROUTINE joinwrfncdf_mem
!
!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE joinwrfncdf                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE  joinwrfncdf_io(filenames,nfile,attadj,jointime,fileconv,    &
                  magnitude_processor,procs,npatch,                     &
                  ids,ide,idss,idse,jds,jde,jdss,jdse,istride,jstride,  &
                  idsinout,ideinout,idssinout,idseinout,                &
                  jdsinout,jdeinout,jdssinout,jdseinout,                &
                  kps,kpe,kpss,kpse,                                    &
                  outdirname,filetail,nvarout,varlists,debug,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Join WRF files in netCDF patches into one large piece.
!    This version does more IO operations but use less memory.
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang (04/27/2007)
!
! MODIFICATIONS:
!
! 03/09/2013 (Y. Wang)
! Added extractions of output box.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: nfile
  LOGICAL, INTENT(IN)            :: attadj
  LOGICAL, INTENT(IN)            :: jointime
  INTEGER, INTENT(IN)            :: fileconv, magnitude_processor
  INTEGER, INTENT(IN)            :: npatch
  INTEGER, INTENT(IN)            :: procs(npatch)
  INTEGER, INTENT(IN)            :: ids,ide,idss,idse,jds,jde,jdss,jdse
  INTEGER, INTENT(IN)            :: idsinout, ideinout, jdsinout, jdeinout
  INTEGER, INTENT(IN)            :: idssinout, idseinout, jdssinout, jdseinout
  INTEGER, INTENT(IN)            :: kps,kpe,kpss,kpse
  INTEGER, INTENT(IN)            :: istride, jstride
  INTEGER, INTENT(INOUT)         :: nvarout
  INTEGER, INTENT(IN)            :: debug
  INTEGER, INTENT(OUT)           :: istatus

  CHARACTER(LEN=*),  INTENT(IN)  :: filenames(nfile)
  CHARACTER(LEN=*),  INTENT(IN)  :: outdirname
  CHARACTER(LEN=5),  INTENT(IN)  :: filetail
  CHARACTER(LEN=40), INTENT(IN)  :: varlists(nvarout)

!
!-----------------------------------------------------------------------
!
! Including files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nf, nvar, n
  INTEGER :: strlen
  LOGICAL :: ispatch(NF_MAX_VARS)

  CHARACTER(LEN=256) :: infilename, outfilename
  INTEGER :: finid, foutid

  INTEGER :: idsout, ideout, jdsout, jdeout
  INTEGER :: idssout, idseout, jdssout, jdseout

  !LOGICAL :: paddingx, paddingy
  !INTEGER :: nxdim, nydim
  REAL    :: newctrlat, newctrlon, lat1, lat2,lon1,lon2
  INTEGER :: newlatidx, newlonidx, oldlatidx1, oldlonidx1, oldlatidx2, oldlonidx2
  LOGICAL :: latavg, lonavg
  INTEGER :: latset, lonset

  !
  ! Dimension variables
  !
  CHARACTER(LEN=32), PARAMETER :: xdimname  = 'west_east_stag'
  CHARACTER(LEN=32), PARAMETER :: ydimname  = 'south_north_stag'
  CHARACTER(LEN=32), PARAMETER :: xsdimname = 'west_east'
  CHARACTER(LEN=32), PARAMETER :: ysdimname = 'south_north'
  CHARACTER(LEN=32) :: diminnames(NF_MAX_DIMS)
  CHARACTER(LEN=32) :: dimname

  INTEGER :: nxid, nyid, nxlg, nylg, nxsid, nysid, nxslg, nyslg
  INTEGER :: narrsize, narrisizemax, narrasizemax
  INTEGER :: unlimdimid, unlimdimlen, unlimodimlen, odimid
  INTEGER :: ndims, dimid, dimlen

  INTEGER :: dimina(NF_MAX_DIMS)         ! Dimension size in original file
  !INTEGER :: dimouta(NF_MAX_DIMS)        ! Dimension size in joined files

  !
  ! Attribute variables
  !
  CHARACTER(LEN=32), PARAMETER :: attnm_ips  = 'WEST-EAST_PATCH_START_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_ipe  = 'WEST-EAST_PATCH_END_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_ipss = 'WEST-EAST_PATCH_START_UNSTAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_ipse = 'WEST-EAST_PATCH_END_UNSTAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jps  = 'SOUTH-NORTH_PATCH_START_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jpe  = 'SOUTH-NORTH_PATCH_END_STAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jpss = 'SOUTH-NORTH_PATCH_START_UNSTAG'
  CHARACTER(LEN=32), PARAMETER :: attnm_jpse = 'SOUTH-NORTH_PATCH_END_UNSTAG'
  CHARACTER(LEN=32) :: attname
  INTEGER :: ipsid,  ipeid,  jpsid,  jpeid
  INTEGER :: ipssid, ipseid, jpssid, jpseid
  INTEGER :: ips, ipe, ipss, ipse
  INTEGER :: jps, jpe, jpss, jpse
  INTEGER :: attnum, ngatts
  INTEGER :: istart, jstart, iend, jend, imod, jmod

  CHARACTER(LEN=32), PARAMETER :: attnm_anx = 'WEST-EAST_GRID_DIMENSION'
  CHARACTER(LEN=32), PARAMETER :: attnm_any = 'SOUTH-NORTH_GRID_DIMENSION'
  INTEGER :: anxid, anyid

  CHARACTER(LEN=32), PARAMETER :: attnm_adx = 'DX'
  CHARACTER(LEN=32), PARAMETER :: attnm_ady = 'DY'
  INTEGER :: adxid, adyid
  REAL    :: adx, ady

  !
  ! Dataset varaibles
  !
  INTEGER, PARAMETER :: MAX_RANK = 4    ! Assume the max rank is 5
  CHARACTER(LEN=32) :: varname
  INTEGER :: varid, nvars, ovarid
  INTEGER :: vartype, varndims, varnatts
  INTEGER :: vardimids(MAX_RANK)
  INTEGER :: startidx(MAX_RANK), countidx(MAX_RANK), strideidx(MAX_RANK)
  INTEGER :: outstart(MAX_RANK), outstride(MAX_RANK)
  INTEGER :: vardim, vdimid

  INTEGER :: varidlists(NF_MAX_VARS), varoutidlists(NF_MAX_VARS)

  INTEGER, ALLOCATABLE :: varari(:)
  REAL,    ALLOCATABLE :: vararr(:)
  CHARACTER(LEN=256)   :: tmpstr,fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  varidlists(:) = 0
  nxlg  = ide-ids+1
  nylg  = jde-jds+1
  nxslg = idse-idss+1
  nyslg = jdse-jdss+1

  startidx(:)  = 1
  narrisizemax = 0
  narrasizemax = 0

  unlimdimlen  = 1

  istatus = 0
  DO nf = 1,nfile

!-----------------------------------------------------------------------
!
! First, Create the merged file based on attributes from the first patch
!
!-----------------------------------------------------------------------

    IF (.NOT. jointime .OR. nf == 1) THEN
      strlen = LEN_TRIM(filenames(nf))
      n = INDEX(filenames(nf),'/',.TRUE.)

      WRITE(outfilename,'(3a)') TRIM(outdirname),                         &
                                filenames(nf)(n+1:strlen),filetail

      IF (jointime .AND. npatch == 1) THEN
        WRITE(infilename, '(a)')       TRIM(filenames(nf))
      ELSE
        WRITE(fmtstr,'(a,2(I1,a))') '(a,a,I',magnitude_processor,'.',magnitude_processor,')'
        WRITE(infilename, FMT=fmtstr) TRIM(filenames(nf)),'_',procs(1)
      END IF

      IF (debug > 0) WRITE(6,'(1x,2a)') 'Opening file - ',TRIM(infilename)
      istatus = nf_open(infilename,NF_NOWRITE,finid)   ! Open file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (debug > 0) WRITE(6,'(1x,2a)') 'Creating file - ',TRIM(outfilename)
      !istatus = nf_create(TRIM(outfilename),NF_CLOBBER,foutid)                     ! CDF 1
      istatus = NF_CREATE(TRIM(outfilename),IOR(NF_CLOBBER,NF_64BIT_OFFSET),foutid) ! CDF 2
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      !
      ! Set dimensions
      !
      IF (fileconv < 2) THEN
        istatus = nf_inq_dimid(finid,xdimname,nxid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        istatus = nf_inq_dimid(finid,ydimname,nyid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      ELSE               ! NMM core
        nxid = -1
        nyid = -1
      END IF

      istatus = nf_inq_dimid(finid,xsdimname,nxsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_dimid(finid,ysdimname,nysid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_unlimdim(finid,unlimdimid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_ndims(finid,ndims)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Copying dimensions - ',ndims
      DO dimid = 1,ndims
        istatus = nf_inq_dim(finid,dimid,dimname,dimlen)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        diminnames(dimid) = dimname
        dimina(dimid)  = dimlen             ! Save dimension id and len
        !dimouta(dimid) = dimlen             ! Output dimension id and len
        IF (dimid == nxid) THEN
          dimlen = (ideinout-idsinout)/istride+1
          !dimouta(dimid) = dimlen
        ELSE IF (dimid == nxsid) THEN
          dimlen = (idseinout-idssinout)/istride+1
          !dimouta(dimid) = dimlen
        ELSE IF (dimid == nyid) THEN
          dimlen = (jdeinout-jdsinout)/jstride+1
          !dimouta(dimid) = dimlen
        ELSE IF (dimid == nysid) THEN
          dimlen = (jdseinout-jdssinout)/jstride+1
          !dimouta(dimid) = dimlen
        ELSE IF (dimid == unlimdimid) THEN
          dimlen = NF_UNLIMITED
        END IF

        IF (debug > 0) WRITE(6,'(9x,a,a20,a,I0)') 'Dimension name - ',dimname, ', length = ',dimlen
        istatus = nf_def_dim(foutid,dimname,dimlen,odimid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO

      !
      ! Set Global attributes
      !
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ips),ipsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ipe),ipeid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ipss),ipssid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ipse),ipseid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jps),jpsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jpe),jpeid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jpss),jpssid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_jpse),jpseid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_anx),anxid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_any),anyid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_adx),adxid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_inq_attid(finid,NF_GLOBAL,TRIM(attnm_ady),adyid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_real(finid,NF_GLOBAL,TRIM(attnm_adx),adx)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_real(finid,NF_GLOBAL,TRIM(attnm_ady),ady)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_natts(finid,ngatts)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (attadj) THEN
        idsout  = 1
        ideout  = (ideinout  - idsinout)/istride  + 1
        idssout = 1
        idseout = (idseinout - idssinout)/istride + 1

        jdsout  = 1
        jdeout  = (jdeinout  - jdsinout)/jstride  + 1
        jdssout = 1
        jdseout = (jdseinout - jdssinout)/jstride + 1

        newlatidx = (jdeout-jdsout)/2    ! Added to modify CEN_LAT & CEN_LON
        newlonidx = (ideout-idsout)/2    ! 0 base

        oldlatidx1 = newlatidx*jstride+jdsinout  ! old domain index
        oldlonidx1 = newlonidx*istride+idsinout  ! old domain index

        latavg = .FALSE.
        lonavg = .FALSE.
        latset = 0
        lonset = 0

        IF (MOD((jdeout-jdsout),2) /= 0) latavg = .TRUE.
        IF (MOD((ideout-idsout),2) /= 0) lonavg = .TRUE.

        oldlatidx2 = -99999999
        oldlonidx2 = -99999999
        IF (.NOT. latavg .AND. .NOT. lonavg) THEN
          oldlatidx2 = oldlatidx1-jstride
          oldlonidx2 = oldlonidx1-istride
        END IF

        IF (debug > 0) WRITE(6,'(1x,a,/,1x,a,2L,a,2I4,/,21x,a,2I4,a))') &
          'It is required to adjust global attribute.',                 &
          'latavg,lonavg = ',latavg,lonavg,                             &
          ', lonidx,latidx = ',newlonidx+1,newlatidx+1,                 &
          ', lonidx,latidx = ',oldlonidx1,oldlatidx1,' (in original domain)'

      ELSE
        idsout  = idsinout
        ideout  = ideinout
        idssout = idssinout
        idseout = idseinout

        jdsout  = jdsinout
        jdeout  = jdeinout
        jdssout = jdssinout
        jdseout = jdseinout
      END IF

      IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Copying global attributes - ',ngatts

      DO attnum = 1,ngatts

        istatus = nf_inq_attname(finid,NF_GLOBAL,attnum,attname)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        IF (debug > 0) WRITE(6,'(9x,2a)') 'Attribute name - ',TRIM(attname)

        IF (attnum == ipsid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_ips),NF_INT,1,idsout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == ipeid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_ipe),NF_INT,1,ideout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == jpsid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_jps),NF_INT,1,jdsout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == jpeid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_jpe),NF_INT,1,jdeout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == ipssid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_ipss),NF_INT,1,idssout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == ipseid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_ipse),NF_INT,1,idseout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == jpssid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_jpss),NF_INT,1,jdssout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == jpseid) THEN
          istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_jpse),NF_INT,1,jdseout)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == anxid) THEN      ! Adjust nx
          IF (istride > 1 .OR. attadj) THEN
            istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_anx),NF_INT,1,(ideout-idsout)+1)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE
            istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutid,NF_GLOBAL)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          END IF
        ELSE IF (attnum == anyid) THEN      ! Adjust ny
          IF (jstride > 1 .OR. attadj) THEN
            istatus = NF_PUT_ATT_INT(foutid,NF_GLOBAL,TRIM(attnm_any),NF_INT,1,(jdeout-jdsout)+1)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          ELSE
            istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutid,NF_GLOBAL)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          END IF
        ELSE IF (attnum == adxid) THEN      ! adjust dx
          istatus = NF_PUT_ATT_REAL(foutid,NF_GLOBAL,TRIM(attnm_adx),NF_REAL,1,adx*istride)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE IF (attnum == adyid) THEN      ! adjust dy
          istatus = NF_PUT_ATT_REAL(foutid,NF_GLOBAL,TRIM(attnm_ady),NF_REAL,1,ady*jstride)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        ELSE
          istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutid,NF_GLOBAL)
          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        END IF

      END DO

      !
      ! Define variables
      !
      istatus = nf_inq_nvars(finid,nvars)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (nvarout >= nvars) THEN
        nvarout = nvars
        DO n = 1,nvars
          varidlists(n) = n
        END DO
      ELSE
        nvar = nvarout         ! suppose to process this number
        nvarout = 0            ! actually got
        DO n = 1,nvar
          istatus = nf_inq_varid(finid,TRIM(varlists(n)),ovarid)
          IF (istatus /= NF_NOERR) THEN
            WRITE(6,'(1x,3a)') 'WARNING: Variable ',TRIM(varlists(n)),' not found. Skipped.'
            CYCLE
          END IF
          nvarout = nvarout + 1
          varidlists(nvarout) = ovarid
        END DO
      END IF

      IF (debug > 0) WRITE(6,'(5x,a,I4)') 'Defining variables - ',nvarout

      DO n = 1,nvarout
        varid = varidlists(n)
        istatus = nf_inq_var(finid,varid,varname,vartype,varndims,vardimids,varnatts)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        IF (debug > 0) WRITE(6,'(9x,2a)') 'Variables - ',TRIM(varname)

        ! Dimensions should be in the same order
        istatus = nf_def_var(foutid,varname,vartype,varndims,vardimids,ovarid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        varoutidlists(n) = ovarid

        DO attnum = 1,varnatts          ! Copy variable attributes
         istatus = nf_inq_attname(finid,varid,attnum,attname)
         IF (istatus /= NF_NOERR) CALL handle_err(istatus)

         istatus = nf_copy_att(finid,varid,attname,foutid,ovarid)
         IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        END DO

      END DO

      istatus = nf_enddef(foutid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF(debug > 0) WRITE(6,'(1x,a)') 'Merged file have been defined.'

      istatus = nf_close(finid)                              ! Close file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    END IF         ! File created.

    IF (.NOT. jointime) unlimdimlen = 1
    unlimodimlen = 0

!-----------------------------------------------------------------------
!
! Write each patch to the merged file
!
!-----------------------------------------------------------------------

    ispatch(:) = .FALSE.
    patch_loop: DO n = 1,npatch
      IF (jointime .AND. npatch == 1) THEN
        WRITE(infilename, '(a)')       TRIM(filenames(nf))
      ELSE
        WRITE(fmtstr,'(a,2(I1,a))') '(a,a,I',magnitude_processor,'.',magnitude_processor,')'
        WRITE(infilename, FMT=fmtstr) TRIM(filenames(nf)),'_',procs(n)
      END IF

      IF (debug > 0) WRITE(6,'(1x,2a)') 'Opening file - ',TRIM(infilename)

      istatus = nf_open(TRIM(infilename),NF_NOWRITE,finid)   ! Open file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      !
      ! Get patch indice
      !
      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_ips),ips)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_ipe),ipe)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_ipss),ipss)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_ipse),ipse)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_jps),jps)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_jpe),jpe)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_jpss),jpss)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      istatus = nf_get_att_int(finid,NF_GLOBAL,TRIM(attnm_jpse),jpse)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Get and save dimension size for this patch
!
!-----------------------------------------------------------------------

      istatus = nf_inq_ndims(finid,ndims)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      dimina(:) = 0
      DO dimid = 1,ndims
        istatus = nf_inq_dim(finid,dimid,dimname,dimlen)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        diminnames(dimid) = dimname
        dimina(dimid)     = dimlen             ! Save dimension id and len
      END DO

!-----------------------------------------------------------------------
!
! loop over each variable
!
!-----------------------------------------------------------------------

      var_loop: DO nvar = 1,nvarout

        varid = varidlists(nvar)

        vardimids(:) = 0
        istatus = nf_inq_var(finid,varid,varname,vartype,varndims,vardimids,varnatts)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Handle dimensions
!
!-----------------------------------------------------------------------

        startidx(:)  = 1
        countidx(:)  = 1
        strideidx(:) = 1
        outstart(:)    = 1
        outstride(:)   = 1
        narrsize = 1

        istart = 1  ! in global index of the input domain
        iend   = 1
        jstart = 1
        jend   = 1
        DO vardim = 1, varndims
          vdimid = vardimids(vardim)
          countidx(vardim) = dimina (vdimid)
          IF ( vdimid == nxid) THEN
            IF ( ips > idsinout) THEN
              imod = MOD((ips-idsinout),istride)
              IF (imod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = istride - imod + 1
              END IF
              istart = ips+startidx(vardim)-1    ! start index relative to global domain
            ELSE
              startidx(vardim) = idsinout-ips+1
              istart = idsinout
            END IF

            iend   = MIN(ipe,ideinout)
            IF (iend < istart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (iend-istart)/istride + 1
            END IF
            strideidx(vardim) = istride

            outstart(vardim) = (istart-idsinout)/istride + 1
            ispatch(nvar) = .TRUE.

          ELSE IF ( vdimid == nyid) THEN
            IF (jps > jdsinout) THEN
              jmod = MOD((jps-jdsinout),jstride)
              IF (jmod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = jstride - jmod + 1
              END IF
              jstart = jps+startidx(vardim)-1
            ELSE
              startidx(vardim) = jdsinout-jps+1
              jstart = jdsinout
            END IF

            jend   = MIN(jpe,jdeinout)
            IF (jend < jstart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (jend-jstart)/jstride + 1
            END IF
            strideidx(vardim) = jstride

            outstart(vardim) = (jstart - jdsinout)/jstride + 1
            ispatch(nvar) = .TRUE.

          ELSE IF ( vdimid == nxsid) THEN
            IF (ipss > idssinout) THEN
              imod = MOD((ipss-idssinout),istride)
              IF (imod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = istride - imod + 1
              END IF
              istart = ipss+startidx(vardim)-1
            ELSE
              startidx(vardim) = idssinout-ipss+1
              istart = idssinout
            END IF

            iend   = MIN(ipse,idseinout)
            IF (iend < istart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (iend-istart)/istride + 1
            END IF
            strideidx(vardim) = istride

            outstart(vardim) = (istart - idssinout)/istride + 1
            ispatch(nvar) = .TRUE.
          ELSE IF ( vdimid == nysid) THEN
            IF (jpss > jdssinout) THEN
              jmod = MOD((jpss-jdssinout),jstride)
              IF (jmod == 0) THEN
                startidx(vardim) = 1
              ELSE
                startidx(vardim) = jstride - jmod + 1
              END IF
              jstart = jpss+startidx(vardim)-1
            ELSE
              startidx(vardim) = jdssinout-jpss+1
              jstart = jdssinout
            END IF

            jend   = MIN(jpse,jdseinout)
            IF (jend < jstart) THEN
              countidx(vardim)  = 0
            ELSE
              countidx(vardim)  = (jend-jstart)/jstride + 1
            END IF
            strideidx(vardim) = jstride

            outstart(vardim) = (jstart - jdssinout)/jstride + 1
            ispatch(nvar) = .TRUE.
          ELSE IF (vdimid == unlimdimid) THEN
            outstart(vardim) = unlimdimlen
            IF (unlimodimlen <= 0) THEN
              unlimodimlen = countidx(vardim)
            ELSE
              IF ( unlimodimlen /= countidx(vardim)) THEN
                WRITE(6,'(1x,a,/)') 'ERROR: Inconsisten size for UNLIMITED dimension.'
                istatus = -1
                RETURN
              END IF
            END IF
          ELSE
            outstart(vardim) = 1
          END IF

          IF (countidx(vardim) <= 0) THEN
            IF (debug > 0) THEN
              WRITE(6,'(9x,2a)') 'Processing variables - ',TRIM(varname)
              WRITE(6,'(12x,a,i0,3a)')                                  &
                'Patch ',procs(n),' skipped because dimension "',       &
                TRIM(diminnames(vdimid)),'" has zero length.'
            END IF
            CYCLE var_loop
          END IF

          narrsize = countidx(vardim)*narrsize
        END DO

        IF ( n > 1 .AND. (.NOT. ispatch(nvar)) ) THEN
          IF (debug > 2) WRITE(6,'(9x,3a)') 'Variable ',TRIM(varname),' skipped.'
          CYCLE
        ELSE
          IF (debug > 2) THEN
            WRITE(6,'(9x,3(a,I4))') 'Dimensions in Patch : ',procs(n),', istart = ',istart, ', jstart = ',jstart

            DO vardim = 1,varndims
              vdimid = vardimids(vardim)
              WRITE(6,'(12x,a,4(a,I4))') diminnames(vdimid),                  &
              ', outstart = ',outstart(vardim),', size = ', countidx(vardim), &
              ' <-- start = ',startidx(vardim),', stride = ',strideidx(vardim)
            END DO
          END IF
        END IF

        ! do not have to merge, use values from the first file

!        IF (.NOT. ispatch(nvar)) THEN
!
!          IF (debug > 0) WRITE(6,'(9x,2a)') 'Copying variables - ',TRIM(varname)
!
!          istatus = NF_COPY_VAR(finid,varid,foutid)
!          IF (istatus /= NF_NOERR) CALL handle_err(istatus)
!
!        ELSE
          ovarid = varoutidlists(nvar)

          IF (debug == 1) WRITE(6,'(9x,3a,I4)')                          &
             'Copying variables - ',TRIM(varname),' from patch: ',procs(n)

          SELECT CASE (vartype)

!-----------------------------------------------------------------------
!
! Integers
!
!-----------------------------------------------------------------------

          CASE (NF_INT)

            IF (narrsize > narrisizemax) THEN   ! Allocate input array only when necessary
              IF (ALLOCATED(varari)) DEALLOCATE(varari, STAT = istatus)
              ALLOCATE(varari(narrsize), STAT = istatus)
              narrisizemax = narrsize
            END IF

            istatus = NF_GET_VARS_INT(finid,varid,startidx,countidx,strideidx,varari)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            istatus = nf_put_vara_INT(foutid,ovarid,outstart,countidx,varari)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Reals
!
!-----------------------------------------------------------------------

          CASE (NF_FLOAT)

            IF (narrsize > narrasizemax) THEN   ! Allocate input array only when necessary
              IF (ALLOCATED(vararr)) DEALLOCATE(vararr, STAT = istatus)
              ALLOCATE(vararr(narrsize), STAT = istatus)
              narrasizemax = narrsize
            END IF

            istatus = NF_GET_VARS_REAL(finid,varid,startidx,countidx,strideidx,vararr)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            istatus = nf_put_vara_REAL(foutid,ovarid,outstart,countidx,vararr)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
! Character string
!
!-----------------------------------------------------------------------

          CASE (NF_CHAR)

            istatus = NF_GET_VARS_TEXT(finid,varid,startidx,countidx,strideidx,tmpstr)
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

            istatus = nf_put_vara_TEXT(foutid,ovarid,outstart,countidx,TRIM(tmpstr))
            IF (istatus /= NF_NOERR) CALL handle_err(istatus)

          CASE DEFAULT
            WRITE(6,'(1x,a,I2)') 'ERROR: unsupported variable type = ',vartype
            istatus = -4
            RETURN
          END SELECT
!        END IF  ! ispatch(nvar)

        !
        ! Added to find the new CEN_LON & CEN_LAT
        !
        IF (attadj .AND. (oldlatidx1 >= jps .AND. oldlatidx1 <= jpe) .AND.  &
                         (oldlonidx1 >= ips .AND. oldlonidx1 <= ipe) ) THEN

          IF ( .NOT. lonavg .AND. latavg ) THEN      ! at U point

            IF (TRIM(varname) == 'XLONG_U' .OR. TRIM(varname) == 'XLAT_U' ) THEN
              newlatidx = (oldlatidx1-jstart)/jstride     ! 0-base
              newlonidx = (oldlonidx1-istart)/istride+1   ! 1-base
              strlen = newlatidx*countidx(1)+newlonidx

              IF (TRIM(varname) == 'XLONG_U' ) THEN
                newctrlon = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                   '*new cen_lon at U point. lonidx = ',newlonidx,      &
                   ', newctrlon = ',newctrlon,' from patch ', procs(n)
                lonset = lonset+1
              ELSE IF (TRIM(varname) == 'XLAT_U') THEN
                newctrlat = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                     '*new cen_lat at U point. latidx = ',newlatidx,    &
                     ', newctrlat = ',newctrlat,' from patch ', procs(n)
                latset = latset+1
              END IF
            END IF

          END IF

          IF ( lonavg .AND. .NOT. latavg ) THEN      ! at V point
            IF (TRIM(varname) == 'XLONG_V' .OR. TRIM(varname) == 'XLAT_V'  ) THEN
              newlatidx = (oldlatidx1-jstart)/jstride
              newlonidx = (oldlonidx1-istart)/istride+1
              strlen = newlatidx*countidx(1)+newlonidx

              IF (TRIM(varname) == 'XLONG_V' ) THEN
                newctrlon = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                   '*new cen_lon at V point. lonidx = ',newlonidx,      &
                   ', newctrlon = ',newctrlon,' from patch ', procs(n)

                lonset = lonset+1
              ELSE IF (TRIM(varname) == 'XLAT_V') THEN
                newctrlat = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                   '*new cen_lat at V point. latidx = ',newlatidx+1,    &
                   ', newctrlat = ',newctrlat,' from patch ', procs(n)
                latset = latset+1
              END IF

            END IF
          END IF

          IF ( lonavg .AND. latavg ) THEN            ! at M point
            IF (TRIM(varname) == 'XLONG' .OR. TRIM(varname) == 'XLAT'  ) THEN
              newlatidx = (oldlatidx1-jstart)/jstride
              newlonidx = (oldlonidx1-istart)/istride+1
              strlen = newlatidx*countidx(1)+newlonidx
              IF (TRIM(varname) == 'XLONG' ) THEN
                newctrlon = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                   '*new cen_lon at M point. lonidx = ',newlonidx,      &
                   ', newctrlon = ',newctrlon,' from patch ',procs(n)
                lonset = lonset+1
              ELSE IF (TRIM(varname) == 'XLAT' ) THEN
                newctrlat = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,I4,a,F8.2,a,i3)')         &
                   '*new cen_lat at M point. latidx = ',newlatidx,      &
                   ', newctrlat = ',newctrlat,' from patch ', procs(n)
                latset = latset+1
              END IF
            END IF

          END IF

          IF ( .NOT. lonavg .AND. .NOT. latavg ) THEN ! at grid point
            IF (TRIM(varname) == 'XLONG_U' .OR. TRIM(varname) == 'XLAT_V'  ) THEN
              newlatidx = (oldlatidx1-jstart)/jstride
              newlonidx = (oldlonidx1-istart)/istride+1
              strlen = newlatidx*countidx(1)+newlonidx
              IF (TRIM(varname) == 'XLONG_U' ) THEN
                lon1 = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')        &
                    '*set lon1 = ',lon1,' at = ',newlatidx+1,newlonidx, &
                    ' from patch ', procs(n)
              ELSE IF (TRIM(varname) == 'XLAT_V' ) THEN
                lat1 = vararr(strlen)
                IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')        &
                    '*set lat1 = ',lat1,' at = ',newlatidx+1,newlonidx, &
                    ' from patch ', procs(n)
              END IF
            END IF

          END IF

        END IF

        IF (attadj .AND. (oldlatidx1 >= jps .AND. oldlatidx1 <= jpe) .AND.  &
                         (oldlonidx2 >= ips .AND. oldlonidx2 <= ipe) ) THEN

          IF ( .NOT. lonavg .AND. .NOT. latavg ) THEN ! at grid point

            IF (TRIM(varname) == 'XLAT_V') THEN
              newlatidx = (oldlatidx1-jstart)/jstride
              newlonidx = (oldlonidx2-istart)/istride+1
              strlen = newlatidx*countidx(1)+newlonidx
              lat2 = vararr(strlen)
              IF (debug > 1) WRITE(6,'(9x,a,F8.2,a,2I4)')               &
                  'set lat2 = ',lat2,' at = ',newlatidx+1,newlonidx
            END IF

          END IF

        END IF

        IF (attadj .AND. (oldlatidx2 >= jps .AND. oldlatidx2 <= jpe) .AND.  &
                         (oldlonidx1 >= ips .AND. oldlonidx1 <= ipe) ) THEN

          IF ( .NOT. lonavg .AND. .NOT. latavg ) THEN ! at grid point

            IF (TRIM(varname) == 'XLONG_U') THEN
              newlatidx = (oldlatidx2-jstart)/jstride
              newlonidx = (oldlonidx1-istart)/istride+1
              strlen = newlatidx*countidx(1)+newlonidx
              lon2 = vararr(strlen)
              IF (debug > 1) WRITE(6,'(1x,a,F8.2,a,2I4,a,i3)')          &
                 '*set lon2 = ',lon2,' at = ',newlatidx+1,newlonidx,    &
                 ' from patch ', procs(n)
            END IF

          END IF

        END IF

      END DO var_loop

      istatus = nf_close(finid)                              ! Close file
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    END DO patch_loop

    unlimdimlen = unlimdimlen + unlimodimlen                 ! Add # of time levels
                                                             ! in output file

!-----------------------------------------------------------------------
!
! Close the output file if applicable
!
!-----------------------------------------------------------------------

    IF (.NOT. jointime .OR. nf == nfile) THEN

      IF (attadj) THEN

        IF ( .NOT. lonavg .AND. .NOT. latavg .AND.                      &
             latset == 0 .AND. lonset == 0) THEN ! at grid point
          newctrlat = 0.5*(lat1+lat2)
          newctrlon = 0.5*(lon1+lon2)
          latset = latset+1
          lonset = lonset+1
        END IF

        IF (latset /= 1 .OR. lonset /= 1) THEN
          WRITE(*,'(1x,2(a,i0),8x,a)')                                  &
            'ERROR: latset = ',latset,', lonset = ',lonset,             &
            'Program aborting ...     Please report.'
          STOP
        END IF

        WRITE(6,'(/,1x,2(a,F8.2),/)')                                   &
           'Changing CEN_LAT & CEN_LON to ',newctrlat,', ',newctrlon

        istatus = NF_REDEF( foutid )
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        istatus = NF_PUT_ATT_REAL(foutid,NF_GLOBAL,'CEN_LAT',NF_REAL,1,newctrlat)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        istatus = NF_PUT_ATT_REAL(foutid,NF_GLOBAL,'MOAD_CEN_LAT',NF_REAL,1,newctrlat)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        istatus = NF_PUT_ATT_REAL(foutid,NF_GLOBAL,'CEN_LON',NF_REAL,1,newctrlon)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)

        istatus = NF_ENDDEF(foutid)
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END IF

      istatus = nf_close(foutid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    END IF

  END DO

  RETURN
END SUBROUTINE joinwrfncdf_io
!
!##################################################################
!##################################################################
!######                                                      ######
!######       SUBROUTINE get_wrf_patch_indices               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_patch_indices(filename,io_form,ips,ipe,ipss,ipse,    &
                                 jps,jpe,jpss,jpse,kps,kpe,kpss,kpse,   &
                                 nx,ny,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Get the size of data patch stored in the WRF data file
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (04/26/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(OUT) :: ips, ipe, jps, jpe
  INTEGER,          INTENT(OUT) :: ipss,ipse,jpss,jpse
  INTEGER,          INTENT(OUT) :: kps, kpe, kpss, kpse
  INTEGER,          INTENT(OUT) :: nx,ny
  INTEGER,          INTENT(OUT) :: istatus

  INCLUDE 'netcdf.inc'

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER           :: ncid
  CHARACTER(LEN=80) :: errmsg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
  IF (io_form == 7) THEN
    istatus = NF_OPEN(TRIM(filename),NF_NOWRITE,ncid)
    IF(istatus /= NF_NOERR)  THEN
       print*,'ERROR with file: ',trim(filename)
       GO TO 999
    END IF

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WEST-EAST_PATCH_START_STAG',ips)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WEST-EAST_PATCH_END_STAG',ipe)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WEST-EAST_PATCH_START_UNSTAG',ipss)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WEST-EAST_PATCH_END_UNSTAG',ipse)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SOUTH-NORTH_PATCH_START_STAG',jps)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SOUTH-NORTH_PATCH_END_STAG',jpe)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SOUTH-NORTH_PATCH_START_UNSTAG',jpss)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SOUTH-NORTH_PATCH_END_UNSTAG',jpse)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WEST-EAST_GRID_DIMENSION',nx)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SOUTH-NORTH_GRID_DIMENSION',ny)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'BOTTOM-TOP_PATCH_START_STAG',kps)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'BOTTOM-TOP_PATCH_END_STAG',kpe)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'BOTTOM-TOP_PATCH_START_UNSTAG',kpss)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'BOTTOM-TOP_PATCH_END_UNSTAG',kpse)
    IF(istatus /= NF_NOERR)  GO TO 999

    istatus = NF_CLOSE(ncid)
    IF(istatus /= NF_NOERR)  GO TO 999
  ELSE
    istatus   = -1
    ips = 0
    ipe = 0
    ipse= 0
    jps = 0
    jpe = 0
    jpse= 0
    WRITE(6,'(1x,a,/)')       &
      'WARNING: Only support netCDF file at present for patch indices.'
  END IF

  RETURN

  999 CONTINUE
  errmsg = NF_STRERROR(istatus)
  WRITE(6,'(1x,2a)') 'NetCDF error: ',errmsg
  istatus = -1

  RETURN
END SUBROUTINE get_wrf_patch_indices
!
! get WRF/METGRID data file name
!
SUBROUTINE get_wrf_filename(filename_convention,dirname,grid_id,        &
                            year,month,day,hour,minute,second,          &
                            colon,magnitude_processor,proc,splitfile,   &
                            filename,istatus)

!#######################################################################

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: filename_convention
  INTEGER, INTENT(IN) :: grid_id
  INTEGER, INTENT(IN) :: year,month,day,hour,minute,second
  INTEGER, INTENT(IN) :: magnitude_processor
  INTEGER, INTENT(IN) :: proc
  LOGICAL, INTENT(IN) :: splitfile

  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  CHARACTER(LEN=1), INTENT(IN)  :: colon
  CHARACTER(LEN=*), INTENT(OUT) :: filename
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  CHARACTER(LEN=3 ) :: faffix
  CHARACTER(LEN=80) :: fheader, ftrailer
  CHARACTER(LEN=80) :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  SELECT CASE (filename_convention)
  CASE (0,3)
    WRITE(fheader,'(a)') 'wrfout_d'
    WRITE(faffix, '(a)') '   '
  CASE (1)
    WRITE(fheader,'(a)') 'met_em.d'
    WRITE(faffix, '(a)') '.nc'
  CASE (2)
    WRITE(fheader,'(a)') 'met_nmm.d'
    WRITE(faffix, '(a)') '.nc'
  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: Unsupported file name convention - ',filename_convention
    istatus = -1
    RETURN
  END SELECT

  IF (splitfile) THEN
    WRITE(fmtstr,'(a,2(I1,a))') '(2a,I',magnitude_processor,'.',magnitude_processor,')'
    WRITE(ftrailer,FMT=TRIM(fmtstr)) TRIM(faffix),'_',proc
  ELSE
    WRITE(ftrailer,FMT='(a)')        TRIM(faffix)
  END IF

  WRITE(filename,FMT='(a,a,I2.2,a,I4.4,5(a,I2.2),a)')                   &
      TRIM(dirname),TRIM(fheader),grid_id,'_',                          &
      year,'-',month,'-',day,'_',hour,colon,minute,colon,second,        &
      TRIM(ftrailer)

  RETURN
END SUBROUTINE get_wrf_filename

!#######################################################################
!
! Handle netCDF error
!
!#######################################################################

SUBROUTINE handle_err(istat)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: istat
  INCLUDE 'netcdf.inc'

  IF (istat /= NF_NOERR) THEN
    PRINT *, TRIM(nf_strerror(istat))
    STOP 'NetCDF error!'
  END IF

  RETURN
END SUBROUTINE handle_err
