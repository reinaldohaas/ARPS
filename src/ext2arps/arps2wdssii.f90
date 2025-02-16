PROGRAM arps2wdssii
!#######################################################################
!
! Convert ARPS NetCDF file into netCDF files suitable for viewing on WDSSII
!
! Keith Brewster, CAPS
! April, 2008
!
!-----------------------------------------------------------------------
!
! MODIFICATION HISTORY
!
!   02/25/2009 Y. Wang
!   Added capability for reading split files. Furthermore the valid time
!   and the file name is also corrected.
!
!   Usage: arps2wdssii [filename] [xpatch ypatch]
!
!   For example:
!
!    o arps2wdssii       # No command line argument, the program will read
!                        # from standard input
!
!    o arps2wdssii arps25may1998.net003600      # read joined file as before
!
!    o arps2wdssii arps25may1998.net003600 2 2  # Read patch files
!
!   02/04/2010 K. Brewster
!   Added support for processing ARPS HDF files.
!
!   03/15/2012(Y. Wang)
!   Removed include file "globcst.inc" to make the code simple and
!   self-explanation. Furthermore, the subroutine refl2d was rewritten
!   with double precsion to avoid float underflow issues.
!
!   04/03/2012 (Y. Wang)
!   Changed refl2d to refl3d and fixed index bugs (nxlg vs. nxlg-1).
!   Changed GETARG to GET_COMMAND_ARGUMENT, which is Fortran 2003 standard.
!
!   02/18/2013 (K. Thomas)
!   Add one more argument, which will specify a new model for output files.
!   This requires the xpatch/ypatch to be specified.  This preserves backwards
!   compatibility.
!
!#######################################################################

  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nx_stag,ny_stag,nz_stag
  INTEGER :: nlat,nlon,nlat_w,nlon_w

  INTEGER :: xpatch, ypatch

  REAL, ALLOCATABLE  :: p(:,:,:)
  REAL, ALLOCATABLE  :: pt(:,:,:)
  REAL, ALLOCATABLE  :: u(:,:,:)
  REAL, ALLOCATABLE  :: v(:,:,:)
  REAL, ALLOCATABLE  :: qr(:,:,:)
  REAL, ALLOCATABLE  :: qs(:,:,:)
  REAL, ALLOCATABLE  :: qh(:,:,:)
  REAL, ALLOCATABLE  :: rho(:,:,:)
  REAL, ALLOCATABLE  :: ref(:,:)
  REAL, ALLOCATABLE  :: reflc(:,:)
  REAL, ALLOCATABLE  :: uwind(:,:)
  REAL, ALLOCATABLE  :: vwind(:,:)

  CHARACTER(LEN=256) :: infile
  CHARACTER(LEN=256) :: reflfile
  CHARACTER(LEN=256) :: windfile
  CHARACTER(LEN=30)  :: DataUnits
  CHARACTER(LEN=30)  :: TypeName
  CHARACTER(LEN=30)  :: DataType
  CHARACTER(LEN=30)  :: attributes
  CHARACTER(LEN=30)  :: unitunit
  CHARACTER(LEN=30)  :: bkgdunit
  CHARACTER(LEN=30)  :: windmunit
  CHARACTER(LEN=30)  :: unitvalue
  CHARACTER(LEN=30)  :: bkgdvalue
  CHARACTER(LEN=30)  :: windmvalue
  CHARACTER(LEN=30)  :: prefix

  REAL(KIND=8) :: latitude
  REAL(KIND=8) :: longitude
  REAL(KIND=8) :: height
  INTEGER      :: itime
  REAL(KIND=8) :: fractime
  REAL(KIND=8) :: dlat
  REAL(KIND=8) :: dlon
  REAL(KIND=8) :: dlat_w
  REAL(KIND=8) :: dlon_w
  INTEGER      :: NumValid
  REAL, PARAMETER :: MissingData = -99900.
  REAL, PARAMETER :: RangeFolded = -99901.
!
! Output level
!
  INTEGER, PARAMETER :: kout = 4     ! = 0, for composite reflectivity
                                     ! > 0, for reflectivity at a specific level
  INTEGER, PARAMETER :: kwout = 4    ! Wind output level
!
! Stride for wind output
!
  INTEGER :: istride
  REAL, PARAMETER :: dx_wout=2000.
!
! NetCDF variables
!
  INTEGER :: refdims(2)
  INTEGER :: ncid,latid,lonid,refid,uwndid,vwndid
  INTEGER :: tid, uid,vid,ptid,pid,qrid,qsid,qhid
  INTEGER :: nxid,nyid,nzid
!
! Input file variables
!
  CHARACTER(LEN=19) :: atime
  INTEGER :: iyear,imon,iday,ihour,imin,isec,itim1970
  INTEGER :: mapproj
  REAL    :: trulat(2)
  REAL    :: sclfct,trulon,ctrlat,ctrlon
  REAL    :: dx,dy
!
!-----------------------------------------------------------------------
!
! Misc internal variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,ilon,jlat,istatus
  INTEGER :: imid,jmid
  INTEGER :: hinfmt
  REAL, PARAMETER :: rscale=1000.
  REAL    :: ctrx,ctry,xlen,ylen
  REAL    :: nwx,nwy,nex,ney,swx,swy,sex,sey
  REAL    :: nwlat,nwlon,nelat,nelon
  REAL    :: swlat,swlon,selat,selon
  REAL    :: qrmax,qhmax,qsmax,refmax
  REAL    :: r2,rscale2,tk

  INTEGER :: ia, ja, iloc, jloc
  INTEGER :: nxlg, nylg
  REAL    :: vtime
  REAL    :: vtime2

  CHARACTER(LEN=256) :: inbuff
  REAL, ALLOCATABLE  :: ubuff(:,:,:)
  REAL, ALLOCATABLE  :: vbuff(:,:,:)
  REAL, ALLOCATABLE  :: rbuff(:,:,:)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL, ALLOCATABLE :: hmin(:)
  REAL, ALLOCATABLE :: hmax(:)

  INTEGER :: idxcdf,idxhdf
  INTEGER :: sd_id
  LOGICAL :: cdffile,hdffile
  LOGICAL :: back = .true.
  INTEGER :: ice
  INTEGER :: new_model

!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'phycst.inc'
  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! Initialize Reflectivity Attributes
! Note: many of this will be overwritten later.
!
  reflfile='wdssrefl.cdf'
  DataUnits='dBZ'
  NumValid=-1
  TypeName='MergedReflectivity'
  DataType='LatLonGrid'
  Latitude=35.0
  Longitude=-99.0
  height=0.
  itime=0
  fractime=0.
  attributes =' BackgroundValue Unit'
  bkgdunit='dimensionless'
  bkgdvalue='-99903'
  unitunit='dimensionless'
  unitvalue='dBZ'
  dlat=0.01
  dlon=0.01
  new_model = 0

!-----------------------------------------------------------------------
!
! Get input arguments
!
!-----------------------------------------------------------------------

  !CALL GETARG(1,infile)
  CALL GET_COMMAND_ARGUMENT(1, infile, k, istatus )

  IF (k <= 0) THEN

    WRITE(6,FMT='(1x,a)',ADVANCE='NO') 'ARPS datafile to be read: '
    READ(5,*) infile

    WRITE(6,FMT='(/,1x,a)',ADVANCE='NO') 'Number of patches in X direction: '
    READ(5,*) xpatch

    WRITE(6,FMT='(/,1x,a)',ADVANCE='NO') 'Number of patches in Y direction: '
    READ(5,*) ypatch

  ELSE
    inbuff = ' '
    !CALL GETARG(2,inbuff)
    CALL GET_COMMAND_ARGUMENT(2, inbuff, k, istatus )
    IF (k > 0) THEN
      READ(inbuff,*) xpatch
      inbuff = ' '
      !CALL GETARG(3,inbuff)
      CALL GET_COMMAND_ARGUMENT(3, inbuff, k, istatus )
      !istatus = LEN_TRIM(inbuff)
      IF (k <= 0) THEN
        WRITE(6,'(1x,a,/,/,1x,a,/)')                                    &
           'ERROR: argument "ypatch" is not found.',                    &
           'Usage: arps2wdssii [filename] [xpatch ypatch] [prefix] [time]'
        STOP
      ELSE
        READ(inbuff,*) ypatch
      END IF
      CALL GET_COMMAND_ARGUMENT(4, prefix, k, istatus)
      IF (k > 0) THEN
        new_model = 1
        vtime = 0
        IF (prefix(1:1) == 'f') THEN
          inbuff = ' '
          CALL GET_COMMAND_ARGUMENT(5, inbuff, k, istatus)
          READ(inbuff,*) vtime2
        ENDIF
      END IF
    ELSE
      xpatch = 1
      ypatch = 1
    END IF
  END IF

  hinfmt=1
  idxcdf=INDEX(TRIM(infile),'.net',back)
  cdffile=(idxcdf > 0)
  IF(cdffile) hinfmt = 4
  idxhdf=INDEX(TRIM(infile),'.hdf',back)
  hdffile=(idxhdf > 0)
  IF(hdffile) hinfmt = 3

  IF( hdffile ) THEN
    WRITE(6,'(a,a,a)') ' file ',TRIM(infile),' is HDF format'
  ELSE IF ( cdffile ) THEN
    WRITE(6,'(a,a,a)') ' file ',TRIM(infile),' is NetCDF'
  ELSE
    WRITE(6,'(a,a)') 'File format of ',TRIM(infile),' not recognized.'
    STOP
  END IF

  ice = 0

  IF( cdffile ) THEN
!-----------------------------------------------------------------------
!
! Read ARPS netCDF File
!
!-----------------------------------------------------------------------

    DO jloc = 1, ypatch
      DO iloc = 1, xpatch
        IF (xpatch > 1 .OR. ypatch > 1) THEN
          CALL gtsplitfn(infile,xpatch,ypatch,1,1,iloc,jloc,            &
                         0,0,1,2,inbuff,istatus)
        ELSE
          WRITE(inbuff,'(a)') TRIM(infile)
        END IF

        WRITE(6,'(1x,3a,/)') 'Opening ARPS netCDF file "',TRIM(inbuff),'"'

        istatus=NF_OPEN(TRIM(inbuff),NF_NOWRITE,ncid)
        IF (istatus /= 0) THEN
          inbuff = NF_STRERROR(istatus)
          WRITE(6,'(/1x,2a)') 'NetCDF error: ',inbuff
          STOP
        END IF

        IF (iloc == 1 .AND. jloc == 1) THEN
          !
          ! Get dimensions
          !
          istatus=NF_INQ_DIMID (ncid,'x',nxid)
          istatus=NF_INQ_DIMLEN(ncid,nxid,nx)
          istatus=NF_INQ_DIMID (ncid,'y',nyid)
          istatus=NF_INQ_DIMLEN(ncid,nyid,ny)
          istatus=NF_INQ_DIMID (ncid,'z',nzid)
          istatus=NF_INQ_DIMLEN(ncid,nzid,nz)
          istatus=NF_INQ_DIMID (ncid,'x_stag',nxid)
          istatus=NF_INQ_DIMLEN(ncid,nxid,nx_stag)
          istatus=NF_INQ_DIMID (ncid,'y_stag',nyid)
          istatus=NF_INQ_DIMLEN(ncid,nyid,ny_stag)
          istatus=NF_INQ_DIMID (ncid,'z_stag',nzid)
          istatus=NF_INQ_DIMLEN(ncid,nzid,nz_stag)

          print *, ' nx,ny,nz:               ',nx,ny,nz
          print *, ' nx_stag,ny_stag,nz_stag:',nx_stag,ny_stag,nz_stag

          nxlg = (nx_stag-3)*xpatch + 3
          nylg = (ny_stag-3)*ypatch + 3
          print *, ' nxlg,nylg:              ',nxlg, nylg
          !
          ! Allocate Arrays
          !
          ALLOCATE (rbuff(nx,      ny,      nz),STAT=istatus)
          ALLOCATE (ubuff(nx_stag, ny,      nz),STAT=istatus)
          ALLOCATE (vbuff(nx,      ny_stag, nz),STAT=istatus)

          ALLOCATE (p  (nxlg-1,nylg-1,nz), STAT=istatus)
          p=0.
          ALLOCATE (pt (nxlg-1,nylg-1,nz), STAT=istatus)
          pt=0.
          ALLOCATE (u  (nxlg,  nylg-1,nz), STAT=istatus)
          u=0.
          ALLOCATE (v  (nxlg-1,nylg,  nz), STAT=istatus)
          v=0.
          ALLOCATE (qr (nxlg-1,nylg-1,nz),   STAT=istatus)
          qr=0.
          ALLOCATE (qs (nxlg-1,nylg-1,nz),   STAT=istatus)
          qs=0.
          ALLOCATE (qh (nxlg-1,nylg-1,nz),   STAT=istatus)
          qh=0.

          !
          print *, ' Done allocating'
          !
          ! Get Global Attributes
          !
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',dx)
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',dy)
          print *, ' dx,dy:                          ',dx,dy
          istatus=NF_GET_ATT_TEXT(ncid,NF_GLOBAL,'INITIAL_TIME',atime)
          print *, ' atime:                             ',atime
          READ(atime,'(i4,5(1x,i2))') iyear,imon,iday,ihour,imin,isec
          WRITE(6,'(1x,a,6I6)') ' time ints:                       ',iyear,imon,iday,ihour,imin,isec
          istatus=NF_GET_ATT_INT(ncid,NF_GLOBAL, 'MAPPROJ', mapproj)
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',  sclfct)
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',trulat(1))
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',trulat(2))
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELON', trulon)
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',  ctrlat)
          istatus=NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',  ctrlon)
          print *, ' mapproj =               ',mapproj
          WRITE(6,'(1x,a,4F12.5)') ' sclfct,trulat1,trulat2,trulon:',sclfct,trulat(1),trulat(2),trulon
          WRITE(6,'(1x,a,2F12.5)') ' ctrlat,ctrlon:                ',ctrlat,ctrlon

          istatus = NF_INQ_VARID(ncid,'Time',tid)
          istatus = NF_GET_VARA_REAL(ncid,tid,(/1/),(/1/),vtime)
          WRITE(6,'(1x,a,F12.5)')  ' Data valid time is:           ',vtime
          WRITE(6,*)
        END IF

        !
        ! Read variables from input file
        !
        istatus=NF_INQ_VARID(ncid,'U',uid)
        istatus=NF_GET_VAR_REAL(ncid,uid,ubuff)
        DO k = 1,nz
          DO j = 1,ny
            ja = (jloc-1)*(ny_stag-3)+j
            DO i = 1,nx_stag
              ia = (iloc-1)*(nx_stag-3)+i
              u(ia,ja,k) = ubuff(i,j,k)
            END DO
          END DO
        END DO

        istatus=NF_INQ_VARID(ncid,'V',vid)
        istatus=NF_GET_VAR_REAL(ncid,vid,vbuff)
        DO k = 1,nz
          DO j = 1,ny_stag
            ja = (jloc-1)*(ny_stag-3)+j
            DO i = 1,nx
              ia = (iloc-1)*(nx_stag-3)+i
              v(ia,ja,k) = vbuff(i,j,k)
            END DO
          END DO
        END DO

        istatus=NF_INQ_VARID(ncid,'PT',ptid)
        istatus=NF_GET_VAR_REAL(ncid,ptid,rbuff)
        DO k = 1,nz
          DO j = 1,ny
            ja = (jloc-1)*(ny_stag-3)+j
            DO i = 1,nx
              ia = (iloc-1)*(nx_stag-3)+i
              pt(ia,ja,k) = rbuff(i,j,k)
            END DO
          END DO
        END DO

        istatus=NF_INQ_VARID(ncid,'P',pid)
        istatus=NF_GET_VAR_REAL(ncid,pid,rbuff)
        DO k = 1,nz
          DO j = 1,ny
            ja = (jloc-1)*(ny_stag-3)+j
            DO i = 1,nx
              ia = (iloc-1)*(nx_stag-3)+i
              p(ia,ja,k) = rbuff(i,j,k)
            END DO
          END DO
        END DO

        istatus=NF_INQ_VARID(ncid,'QR',qrid)
        IF(istatus == 0) THEN
          istatus=NF_GET_VAR_REAL(ncid,qrid,rbuff)
          DO k = 1,nz
            DO j = 1,ny
              ja = (jloc-1)*(ny_stag-3)+j
              DO i = 1,nx
                ia = (iloc-1)*(nx_stag-3)+i
                qr(ia,ja,k) = rbuff(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          WRITE(6,'(a)') ' qr missing, skipping'
        END IF


        istatus=NF_INQ_VARID(ncid,'QS',qsid)
        IF(istatus == 0) THEN
          istatus=NF_GET_VAR_REAL(ncid,qsid,rbuff)
          DO k = 1,nz
            DO j = 1,ny
              ja = (jloc-1)*(ny_stag-3)+j
              DO i = 1,nx
                ia = (iloc-1)*(nx_stag-3)+i
                qs(ia,ja,k) = rbuff(i,j,k)
              END DO
            END DO
          END DO

          ice = 1
        ELSE
          WRITE(6,'(a)') ' qs missing, skipping'
        END IF

        istatus=NF_INQ_VARID(ncid,'QH',qhid)
        IF(istatus == 0) THEN
          istatus=NF_GET_VAR_REAL(ncid,qhid,rbuff)
          DO k = 1,nz
            DO j = 1,ny
              ja = (jloc-1)*(ny_stag-3)+j
              DO i = 1,nx
                ia = (iloc-1)*(nx_stag-3)+i
                qh(ia,ja,k) = rbuff(i,j,k)
              END DO
            END DO
          END DO


        ELSE
          WRITE(6,'(a)') ' qh missing, skipping'
        END IF

        !
        ! Close input NetCDF file
        !
        istatus=NF_CLOSE(ncid)
      END DO
    END DO

    DEALLOCATE(ubuff, vbuff, rbuff)

  ELSE IF ( hdffile ) THEN

!-----------------------------------------------------------------------
!
! Read ARPS HDF File
!
!-----------------------------------------------------------------------

    DO jloc = 1, ypatch
      DO iloc = 1, xpatch
        IF (xpatch > 1 .OR. ypatch > 1) THEN
          CALL gtsplitfn(infile,xpatch,ypatch,1,1,iloc,jloc,            &
                         0,0,1,2,inbuff,istatus)
        ELSE
          WRITE(inbuff,'(a)') TRIM(infile)
        END IF

        WRITE(6,'(1x,3a,/)') 'Processing ARPS HDF file "',TRIM(inbuff),'"'

        CALL hdfopen(inbuff,1,sd_id)

        IF (iloc == 1 .AND. jloc == 1) THEN
          !
          ! Get dimensions
          !
          print *, ' hdf open sd_id: ',sd_id
          CALL hdfrdi(sd_id,'nx',nx,istatus)
          CALL hdfrdi(sd_id,'ny',ny,istatus)
          CALL hdfrdi(sd_id,'nz',nz,istatus)
          print *, ' nx,ny,nz:  ',nx,ny,nz

          nxlg = (nx-3)*xpatch + 3
          nylg = (ny-3)*ypatch + 3
          print *, ' nxlg,nylg: ',nxlg, nylg
          !
          ! Allocate Arrays
          !
          ALLOCATE (rbuff(nx,ny,nz),STAT=istatus)
          ALLOCATE (ubuff(nx,ny,nz),STAT=istatus)
          ALLOCATE (vbuff(nx,ny,nz),STAT=istatus)

          ALLOCATE (p  (nxlg-1,nylg-1,nz), STAT=istatus)
          ALLOCATE (pt (nxlg-1,nylg-1,nz), STAT=istatus)
          ALLOCATE (u  (nxlg,  nylg-1,nz), STAT=istatus)
          ALLOCATE (v  (nxlg-1,nylg,  nz), STAT=istatus)
          pt=0.
          p=0.
          u=0.
          v=0.
          ALLOCATE (qr (nxlg-1,nylg-1,nz), STAT=istatus)
          ALLOCATE (qs (nxlg-1,nylg-1,nz), STAT=istatus)
          ALLOCATE (qh (nxlg-1,nylg-1,nz), STAT=istatus)
          qr=0.
          qs=0.
          qh=0.

          ALLOCATE (itmp(nx,ny,nz),stat=istatus)
          ALLOCATE (hmax(nz),stat=istatus)
          ALLOCATE (hmin(nz),stat=istatus)

          print *, ' Done allocating'
          !
          ! Get Static Attributes
          !
          CALL hdfrdr(sd_id,'dx',dx,istatus)
          CALL hdfrdr(sd_id,'dy',dy,istatus)
          print *, ' dx,dy:                          ',dx,dy
          CALL hdfrdi(sd_id,'year',iyear,istatus)
          CALL hdfrdi(sd_id,'month',imon,istatus)
          CALL hdfrdi(sd_id,'day',iday,istatus)
          CALL hdfrdi(sd_id,'hour',ihour,istatus)
          CALL hdfrdi(sd_id,'minute',imin,istatus)
          CALL hdfrdi(sd_id,'second',isec,istatus)
          WRITE(6,'(1x,a,6I6)') ' time ints:                       ',   &
                                iyear,imon,iday,ihour,imin,isec
          CALL hdfrdi(sd_id,'mapproj',mapproj,istatus)
          CALL hdfrdr(sd_id,'sclfct',sclfct,istatus)
          CALL hdfrdr(sd_id,'trulat1',trulat(1),istatus)
          CALL hdfrdr(sd_id,'trulat2',trulat(2),istatus)
          CALL hdfrdr(sd_id,'trulon',trulon,istatus)
          CALL hdfrdr(sd_id,'ctrlat',ctrlat,istatus)
          CALL hdfrdr(sd_id,'ctrlon',ctrlon,istatus)
          print *, ' mapproj =               ',mapproj
          WRITE(6,'(1x,a,4F12.5)') ' sclfct,trulat1,trulat2,trulon:',   &
                                   sclfct,trulat(1),trulat(2),trulon
          WRITE(6,'(1x,a,2F12.5)') ' ctrlat,ctrlon:                ',   &
                                   ctrlat,ctrlon
          !
          ! Get Time Dependent Attributes
          !
          CALL hdfrdr(sd_id,'time',vtime,istatus)
          WRITE(6,'(1x,a,F12.5,/)') ' Data valid time is:           ',  &
                                    vtime
        END IF
        !
        ! Read variables from input file
        !
        CALL hdfrd3d(sd_id,'u',nx,ny,nz,ubuff,istatus,itmp,hmax,hmin)
        DO k = 1,nz
          DO j = 1,ny-1
            ja = (jloc-1)*(ny-3)+j
            DO i = 1,nx
              ia = (iloc-1)*(nx-3)+i
              u(ia,ja,k) = ubuff(i,j,k)
            END DO
          END DO
        END DO

        CALL hdfrd3d(sd_id,'v',nx,ny,nz,vbuff,istatus,itmp,hmax,hmin)
        DO k = 1,nz
          DO j = 1,ny
            ja = (jloc-1)*(ny-3)+j
            DO i = 1,nx-1
              ia = (iloc-1)*(nx-3)+i
              v(ia,ja,k) = vbuff(i,j,k)
            END DO
          END DO
        END DO

        CALL hdfrd3d(sd_id,'pt',nx,ny,nz,rbuff,istatus,itmp,hmax,hmin)
        DO k = 1,nz
          DO j = 1,ny-1
            ja = (jloc-1)*(ny-3)+j
            DO i = 1,nx-1
              ia = (iloc-1)*(nx-3)+i
              pt(ia,ja,k) = rbuff(i,j,k)
            END DO
          END DO
        END DO

        CALL hdfrd3d(sd_id,'p',nx,ny,nz,rbuff,istatus,itmp,hmax,hmin)
        DO k = 1,nz
          DO j = 1,ny-1
            ja = (jloc-1)*(ny-3)+j
            DO i = 1,nx-1
              ia = (iloc-1)*(nx-3)+i
              p(ia,ja,k) = rbuff(i,j,k)
            END DO
          END DO
        END DO

        CALL hdfrd3d(sd_id,'qr',nx,ny,nz,rbuff,istatus,itmp,hmax,hmin)
        IF(istatus == 0) THEN
          DO k = 1,nz
            DO j = 1,ny-1
              ja = (jloc-1)*(ny-3)+j
              DO i = 1,nx-1
                ia = (iloc-1)*(nx-3)+i
                qr(ia,ja,k) = rbuff(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          WRITE(6,'(a)') ' qr missing, skipping'
        END IF

        CALL hdfrd3d(sd_id,'qs',nx,ny,nz,rbuff,istatus,itmp,hmax,hmin)
        IF(istatus == 0) THEN
          DO k = 1,nz
            DO j = 1,ny-1
              ja = (jloc-1)*(ny-3)+j
              DO i = 1,nx-1
                ia = (iloc-1)*(nx-3)+i
                qs(ia,ja,k) = rbuff(i,j,k)
              END DO
            END DO
          END DO
          ice = 1
        ELSE
          WRITE(6,'(a)') ' qs missing, skipping'
        END IF

        CALL hdfrd3d(sd_id,'qh',nx,ny,nz,rbuff,istatus,itmp,hmax,hmin)
        IF(istatus == 0) THEN
          DO k = 1,nz
            DO j = 1,ny-1
              ja = (jloc-1)*(ny-3)+j
              DO i = 1,nx-1
                ia = (iloc-1)*(nx-3)+i
                qh(ia,ja,k) = rbuff(i,j,k)
              END DO
            END DO
          END DO
          ice = 1
        ELSE
          WRITE(6,'(a)') ' qh missing, skipping'
        END IF
        !
        ! Close hdf data file
        !
        CALL hdfclose(sd_id,istatus)
      END DO
    END DO
    DEALLOCATE(ubuff, vbuff, rbuff)
    DEALLOCATE(itmp, hmax, hmin)

  END IF   ! hdf file

!-----------------------------------------------------------------------
!
! File reading complete.
! Generate some input data stats
!
!-----------------------------------------------------------------------

  DO k=1,nz
    qrmax=-999.0
    qsmax=-999.0
    qhmax=-999.0
    DO j=1,nylg-1
      DO i=1,nxlg-1
        qrmax=max(qrmax,qr(i,j,k))
        qsmax=max(qsmax,qs(i,j,k))
        qhmax=max(qhmax,qh(i,j,k))
      END DO
    END DO
    WRITE(6,'(1x,a,I4,3(a,e12.5))') ' k = ',k,': qrmax:',(1000.*qrmax), &
                     ',  qsmax:',(1000.*qsmax),',  qhmax:',(1000.*qhmax)
  END DO

  !
  ! Calculate air density, needed for reflectivity calculation
  !
  ALLOCATE (rho(nxlg-1,nylg-1,nz), STAT=istatus)
  ALLOCATE (ref(nxlg-1,nylg-1), STAT=istatus)
  rho=0.
  ref=0.

  rho=1.0
  DO k = 1,nz
    DO j=1,nylg-1
      DO i=1,nxlg-1
        tk=pt(i,j,k)*((p(i,j,k)/p0)**rddcp)
        rho(i,j,k)=p(i,j,k)/(rd*tk)
      END DO
    END DO
  END DO
  !
  ! Get reflectivity
  !
  CALL refl3d(nxlg-1,nylg-1,nz,rho,ice,qr,qs,qh,kout,ref)

  DEALLOCATE(rho)

  !
  ! Grid set-up and location calculations
  !
  CALL setmapr(mapproj,sclfct,trulat,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)

  xlen=(nxlg-3)*dx
  ylen=(nylg-3)*dy
  nwx=ctrx-(0.5*xlen)
  nwy=ctry+(0.5*ylen)
  CALL xytoll(1,1,nwx,nwy,nwlat,nwlon)
  WRITE(6,'(1x,2(a,F12.5))') ' NW lat: ',nwlat, ', NW lon: ',nwlon
  Latitude=nwlat
  Longitude=nwlon

  swx=ctrx-(0.5*xlen)
  swy=ctry-(0.5*ylen)
  CALL xytoll(1,1,swx,swy,swlat,swlon)
  WRITE(6,'(1x,2(a,F12.5))') ' SW lat: ',swlat, ', SW lon: ',swlon

  nex=ctrx+(0.5*xlen)
  ney=ctry+(0.5*ylen)
  CALL xytoll(1,1,nex,ney,nelat,nelon)
  WRITE(6,'(1x,2(a,F12.5))') ' NE lat: ',nelat, ', NE lon: ',nelon

  sex=ctrx+(0.5*xlen)
  sey=ctry-(0.5*ylen)
  CALL xytoll(1,1,sex,sey,selat,selon)
  WRITE(6,'(1x,2(a,F12.5))') ' SE lat: ',selat, ', SE lon: ',selon

  dlat=(0.5*((nelat-selat)+(nwlat-swlat)))/float(nylg-3)
  dlon=(0.5*((nelon-nwlon)+(selon-swlon)))/float(nxlg-3)
  WRITE(6,'(1x,2(a,F10.5))') ' refl dlat:',dlat, ', refl dlon:',dlon

  istride=MAX(1,NINT((dx_wout/dx)))
  WRITE(6,'(a,i6,a,f7.2,a)') ' Found istride=',istride,' aiming for', &
         (0.001*dx_wout),' km wind output.'
  WRITE(6,'(a,f7.2,a)') &
         ' Actual wind output will be at dx=',(istride*0.001*dx),' km'
  dlat_w=istride*dlat
  dlon_w=istride*dlon
  WRITE(6,'(1x,2(a,F10.5))') ' wind dlat:',dlat_w, ', wind dlon:',dlon_w

  nlon_w=(nxlg-2)/istride
  nlat_w=(nylg-2)/istride
  ALLOCATE (uwind(nlon_w,nlat_w),STAT=istatus)
  ALLOCATE (vwind(nlon_w,nlat_w),STAT=istatus)

  !
  ! Begin Processing for output
  !
  CALL ctim2abss( iyear,imon,iday,ihour,imin,isec, itime )
  !
  ! In the original naming scheme, the time part of filename included the
  ! *valid* time of the data.  We have a request to make that the
  ! initialization of the data, with the forecast length elsewhre in the name.
  ! Since we're doing cycled runs, we have to pass the timeoffset for the
  ! first cycle so the right things are done.
  !
  IF (new_model == 0) THEN
    itime = itime + nint(vtime)
  ELSE
    itime = itime + nint(vtime2)
  END IF
  print *, ' itime 1960 base: ',itime

  CALL ctim2abss(  1970,   1,   1,    0,   0,   0, itim1970 )
  !print *, ' itim1970=',itim1970

  CALL abss2ctim(itime,iyear,imon,iday,ihour,imin,isec )
  itime=itime-itim1970
  print *, ' itime 1970 base: ',itime

  nlon=nxlg-2
  nlat=nylg-2
  ALLOCATE (reflc(nlon,nlat),STAT=istatus)

  imid=(nlat/2)+1
  jmid=(nlon/2)+1
  rscale2=1.0/(rscale*rscale)

  refmax=-999.0
  DO jlat=1,nlat
    DO ilon=1,nlon
      i=ilon+1
      j=nylg-jlat
      reflc(ilon,jlat)=ref(i,j)
      refmax=max(refmax,reflc(ilon,jlat))
    END DO
  END DO
  DEALLOCATE(ref)
  print *, ' refmax= ',refmax,kout

  k = kwout
  DO jlat=1,nlat_w
    DO ilon=1,nlon_w
      i=(ilon-1)*istride+2
      j=(nylg-(jlat-1)*istride)-1
      uwind(ilon,jlat)=0.5*(u(i,j,k)+u(i+1,j,k))
      vwind(ilon,jlat)=0.5*(v(i,j,k)+v(i,j+1,k))
      !j=(jlat-1)*istride+2
      !uwind(ilon,jlat)=u(i,j,kout)
      !vwind(ilon,jlat)=v(i,j,kout)
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Write WDSSII files
!
!-----------------------------------------------------------------------

  !
  ! Open file
  !
  IF (new_model == 0) THEN 
    WRITE(reflfile,'(i4.4,2(i2.2),a,3(i2.2),a)')                        &
              iyear,imon,iday,'-',ihour,imin,isec,'refl.netcdf'
  ELSE
    WRITE(reflfile,'(a,a,i4.4,2(i2.2),a,3(i2.2),a,i3.3,a)')             &
      TRIM(prefix),'_',iyear,imon,iday,'_',ihour,imin,isec,             &
      '_',INT(vtime-vtime2)/60,'refl.netcdf'
    WRITE(6,'(a,i)'),'  new_model, vtime2=',INT(vtime2)
  END IF
  istatus=NF_CREATE(reflfile,NF_CLOBBER,ncid)
  !
  ! Create dimensions
  !
  !  istatus=NF_REDEF(ncid)
  istatus=NF_DEF_DIM(ncid,'Lat',nlat,latid)
  istatus=NF_DEF_DIM(ncid,'Lon',nlon,lonid)
  !
  ! Create variables
  !
  refdims(1)=lonid
  refdims(2)=latid
  istatus=NF_DEF_VAR(ncid,TypeName,NF_FLOAT,2,refdims,refid)
  !
  ! Write variable attributes
  !
  istatus=NF_PUT_ATT_TEXT(ncid,refid,'Units',LEN_TRIM(DataUnits),TRIM(DataUnits))
  istatus=NF_PUT_ATT_INT(ncid,refid,'NumValidRuns',NF_INT,1,NumValid)
  !
  ! Write global attributes
  !
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'TypeName',LEN_TRIM(TypeName),TRIM(TypeName))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'DataType',LEN_TRIM(DataType),TRIM(DataType))
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Latitude',NF_DOUBLE,1,Latitude)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Longitude',NF_DOUBLE,1,Longitude)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Height',NF_DOUBLE,1,height)
  istatus=NF_PUT_ATT_INT(ncid,NF_GLOBAL,'Time',NF_INT,1,itime)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'FractionalTime',NF_DOUBLE,1,fractime)
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'attributes',LEN_TRIM(attributes),TRIM(attributes))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'BackgroundValue-unit',LEN_TRIM(bkgdunit),TRIM(bkgdunit))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'BackgroundValue-value',LEN_TRIM(bkgdvalue),TRIM(bkgdvalue))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'Unit-unit',LEN_TRIM(unitunit),TRIM(unitunit))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'Unit-value',LEN_TRIM(unitvalue),TRIM(unitvalue))
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'LatGridSpacing',NF_DOUBLE,1,dlat)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'LonGridSpacing',NF_DOUBLE,1,dlon)
  istatus=NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'MissingData',NF_FLOAT,1,MissingData)
  istatus=NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'RangeFolded',NF_FLOAT,1,RangeFolded)
  istatus=NF_ENDDEF(ncid)
  !
  ! Write variables
  !
  istatus=NF_PUT_VAR_REAL(ncid,refid,reflc)
  !
  ! Close reflectivity file
  !
  istatus=NF_CLOSE(ncid)

!-----------------------------------------------------------------------
!
! Wind file
!
!-----------------------------------------------------------------------
  !
  ! Assign wind attributes
  !
  windfile='wdsswind.cdf'
  DataUnits='MetersPerSecond'
  TypeName='3DVarWindField'
  DataType='WindField'
  attributes =' Unit meanwind'
  unitunit='dimensionless'
  unitvalue='MetersPerSecond'
  windmunit='dimensionless'
  windmvalue='5.0'
  !
  ! Open wind file
  !
  IF (new_model == 0) THEN 
    WRITE(windfile,'(i4.4,2(i2.2),a,3(i2.2),a)')                        &
              iyear,imon,iday,'-',ihour,imin,isec,'wind.netcdf'
  ELSE
    WRITE(windfile,'(a,a,i4.4,2(i2.2),a,3(i2.2),a,i3.3,a)')             &
      TRIM(prefix),'_',iyear,imon,iday,'_',ihour,imin,isec,             &
      '_',INT(vtime-vtime2)/60,'wind.netcdf'
  END IF
  istatus=NF_CREATE(windfile,NF_CLOBBER,ncid)
!
! Create dimensions
!
!  istatus=NF_REDEF(ncid)
  istatus=NF_DEF_DIM(ncid,'Lat',nlat_w,latid)
  istatus=NF_DEF_DIM(ncid,'Lon',nlon_w,lonid)
!
! Create variables
!
  refdims(1)=lonid
  refdims(2)=latid
  istatus=NF_DEF_VAR(ncid,'uArray',NF_FLOAT,2,refdims,uwndid)
  istatus=NF_DEF_VAR(ncid,'vArray',NF_FLOAT,2,refdims,vwndid)
!
! Write global attributes
!
  istatus=NF_PUT_ATT_TEXT(ncid,uwndid,'Units',LEN_TRIM(DataUnits),TRIM(DataUnits))
  istatus=NF_PUT_ATT_TEXT(ncid,vwndid,'Units',LEN_TRIM(DataUnits),TRIM(DataUnits))
!
! Write global attributes
!
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'TypeName',LEN_TRIM(TypeName),TRIM(TypeName))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'DataType',LEN_TRIM(DataType),TRIM(DataType))
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Latitude',NF_DOUBLE,1,Latitude)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Longitude',NF_DOUBLE,1,Longitude)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Height',NF_DOUBLE,1,height)
  istatus=NF_PUT_ATT_INT(ncid,NF_GLOBAL,'Time',NF_INT,1,itime)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'FractionalTime',NF_DOUBLE,1,fractime)
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'attributes',LEN_TRIM(attributes),TRIM(attributes))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'Unit-unit',LEN_TRIM(unitunit),TRIM(unitunit))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'Unit-value',LEN_TRIM(unitvalue),TRIM(unitvalue))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'meanwind-unit',LEN_TRIM(windmunit),TRIM(windmunit))
  istatus=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'meanwind-value',LEN_TRIM(windmvalue),TRIM(windmvalue))
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'LatGridSpacing',NF_DOUBLE,1,dlat_w)
  istatus=NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'LonGridSpacing',NF_DOUBLE,1,dlon_w)
  istatus=NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'MissingData',NF_FLOAT,1,MissingData)
  istatus=NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'RangeFolded',NF_FLOAT,1,RangeFolded)
  istatus=NF_ENDDEF(ncid)
!
! Write variables
!
  istatus=NF_PUT_VAR_REAL(ncid,uwndid,uwind)
  istatus=NF_PUT_VAR_REAL(ncid,vwndid,vwind)
!
! Close file
!
  istatus=NF_CLOSE(ncid)

!-----------------------------------------------------------------------
!
! End the program
!
!-----------------------------------------------------------------------

  WRITE(6,'(1x,a)') ' ==== Successful Completion of ARPS2WDSSII ===='

  STOP
END PROGRAM arps2wdssii

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE REFL2D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE refl3d(nx,ny, nz, rhobar, ice, qr, qs, qh, kout, reflc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: K. Droegemeier and M.Xue
!  4/19/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/6/95 (J. Zong and M. Xue)
!  Added qs and qh to the argument list of this subroutine to
!  facilitate inclusion of the contributions of qs and qh to reflec-
!  tivity. A relation between radar reflectivity factor and snow
!  content is adopted from Rogers and Yau (1989) and extended to
!  represent the effects of snow and graupel/hail on the
!  reflectivity. globcst.inc is included to pass the value of ice.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    rhobar   Base state density (kg/m**3)
!    qr       Rainwater mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!  OUTPUT:
!
!    reflc    Radar reflectivity factor.
!
!-----------------------------------------------------------------------
!
  USE arps_precision
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER, INTENT(IN) :: ice
  INTEGER, INTENT(IN) :: kout

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)

  REAL :: reflc (nx,ny)     ! Radar reflectivity (dBZ)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER  :: i,j,k
  REAL(DP) :: coef, svnfrth
  REAL(DP) :: arg
  REAL(DP) :: rho, qrain, qice
  REAL(P)  :: refx

  INTEGER :: kbgn, kend
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!-----------------------------------------------------------------------
!
!  IF ( mphyopt >= 5 .AND. mphyopt <= 7 ) THEN
!    svnfrth = 1.24
!    coef    = 20740.0
!  ELSE
    svnfrth = 7.0D0/4.0D0
    coef    = 17300.0D0
!  END IF


  kbgn = 1             ! composite reflectivity if kout == 0
  kend = nz-1

  IF (kout > 0) THEN   ! output reflectivity at a specific level only
    kbgn = kout
    kend = kout
  END IF

  reflc(:,:) = 0.0
  DO k = kbgn,kend
    DO j = 1,ny
      DO i = 1,nx
        rho   = REAL(rhobar(i,j,k),DP)

        qrain = MAX(0.0D0,REAL(qr(i,j,k),DP))
        arg = coef*( rho*1000.0 * qrain )**svnfrth

        IF (ice == 1) THEN

          qice = MAX(0.0D0,REAL(qs(i,j,k)+qh(i,j,k),DP))
          arg = arg + 38000.0*( rho*1000.0 * qice )**2.2

        END IF

        arg = MAX(arg,1.0D0)
        refx = 10.0*ALOG10( REAL(arg,P) )

        reflc(i,j) = MAX(reflc(i,j),refx)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE refl3d
