!#######################################################################
SUBROUTINE open_hdf(filename,nxpatch,nypatch,iorig,jorig,ioffset,joffset, &
                    iamroot,inpatch,rdwrtflag,dbglvl,fHandles,istatus)

!-----------------------------------------------------------------------
!
! Purpose:
!   Open multiple HDF 4 files.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN) :: filename
  INTEGER,            INTENT(IN) :: nxpatch     ! number of x patches to be opened
  INTEGER,            INTENT(IN) :: nypatch     ! number of y patches to be opened
  INTEGER,            INTENT(IN) :: iorig       ! The relative loation of the lower-left patch
  INTEGER,            INTENT(IN) :: jorig       ! within the whole domain
  INTEGER,            INTENT(IN) :: ioffset
  INTEGER,            INTENT(IN) :: joffset
  LOGICAL,            INTENT(IN) :: iamroot
  LOGICAL,            INTENT(IN) :: inpatch     ! whether the files are in patches
                                                ! The subroutine handles the following cases
                                                !
                                                ! No-MPI
                                                !     one joined file
                                                !     split files
                                                !
                                                ! MPI
                                                !     one patch each processor
                                                !     multiple patches each processor
                                                !     One joined file for all processors (x)

  INTEGER,            INTENT(IN) :: rdwrtflag   ! = 1 Read only
                                                ! = 2 write only
  INTEGER,            INTENT(IN) :: dbglvl
  INTEGER,            INTENT(OUT) :: fHandles(nxpatch,nypatch)
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER            :: ixpatch, iypatch
  CHARACTER(LEN=256) :: tmpfname
  INTEGER            :: ioflag, sd_id

!-----------------------------------------------------------------------

   INCLUDE 'hdf.f90'

   INTEGER :: sfstart

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (rdwrtflag == 1) THEN
    ioflag = DFACC_READ
  ELSE IF (rdwrtflag == 2) THEN
    ioflag = DFACC_CREATE
  ELSE
    istatus = -1
    RETURN
  END IF

  DO iypatch = 1, nypatch
    DO ixpatch = 1, nxpatch

      tmpfname = ' '
      IF (.NOT. inpatch .AND. iamroot) THEN
        tmpfname = filename
      ELSE
        CALL gtsplitfn(filename,nxpatch,nypatch,iorig,jorig,ixpatch,iypatch, &
                       ioffset-1,joffset-1,MOD(rdwrtflag,2),dbglvl,tmpfname,istatus)
      END IF

      sd_id = sfstart(TRIM(tmpfname),ioflag)

      IF (sd_id <= 0) THEN
        istatus = -2
        WRITE (6,'(1x,3a,I5)') 'ERROR: open_hdf in opening ',          &
                                trim(tmpfname),' with sd_id = ',sd_id
        RETURN
      ELSE
        IF (dbglvl > 2) WRITE (6,'(1x,3a,I15)') 'Open ',trim(tmpfname),' with id ',sd_id
        fHandles(ixpatch,iypatch) = sd_id
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE open_hdf

!#######################################################################
SUBROUTINE close_hdf(fHandles,nxpatch,nypatch,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Close multiple HDF 4 files.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,            INTENT(IN) :: nxpatch     ! number of x patches to be opened
  INTEGER,            INTENT(IN) :: nypatch     ! number of y patches to be opened
  INTEGER,            INTENT(IN) :: fHandles(nxpatch,nypatch)
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER            :: ixpatch, iypatch

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfend

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  DO iypatch = 1, nypatch
    DO ixpatch = 1, nxpatch

      istatus = sfend(fHandles(ixpatch,iypatch))
    END DO
  END DO

  RETURN
END SUBROUTINE close_hdf

!#######################################################################
SUBROUTINE copy_global_attributes(iamroot,fHandles,nxpatch,nypatch,     &
                   nprocx_in, nprocy_in,nxpnt,nypnt,runnamein,          &
                   fHndlOut,nxpatchout,nypatchout,                      &
                   nx,ny,nz,nzsoil,nstyps,numvar,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Copy HDF 4 file attributes from multiple files to multiple files.
!   It is assumed that all input file attributes are the same except for
!   'nx' & 'ny'. So only the first input file is opened here.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  LOGICAL, INTENT(IN)  :: iamroot
  INTEGER, INTENT(IN)  :: nxpatch     ! number of x patches to be opened
  INTEGER, INTENT(IN)  :: nypatch     ! number of y patches to be opened
  INTEGER, INTENT(IN)  :: fHandles(nxpatch,nypatch)
  INTEGER, INTENT(IN)  :: nprocx_in, nprocy_in
  INTEGER, INTENT(IN)  :: nxpnt,nypnt
  CHARACTER(LEN=80)    :: runnamein
  INTEGER, INTENT(IN)  :: nxpatchout, nypatchout
  INTEGER, INTENT(IN)  :: fHndlOut(nxpatchout,nypatchout)
  INTEGER, INTENT(OUT) :: nx, ny, nz, nzsoil, nstyps
  INTEGER, INTENT(OUT) :: numvar
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ixpatch, iypatch
  INTEGER :: nxlg, nylg
  INTEGER :: nxout, nyout, nxlgout, nylgout

  INTEGER :: sd_id1, nattr

  CHARACTER(LEN=256) :: aname
  INTEGER :: aindex, dtype, nvalues

  CHARACTER(LEN=1), ALLOCATABLE :: chararr(:)
  INTEGER                       :: charlen

  REAL    :: buf_r
  INTEGER :: buf_i
  INTEGER :: i

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sffinfo         ! Get number of data sets and number of attr.
  INTEGER :: sfgainfo        ! Get attribute name and value
  INTEGER :: sfrcatt         ! Get   character attribute
  INTEGER :: sfscatt         ! Write character attribute
  INTEGER :: sfrnatt         ! Get   float     attribute
  INTEGER :: sfsnatt         ! Write float     attribute

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  sd_id1 = fHandles(1,1)

  istatus = sffinfo(sd_id1,numvar,nattr)
  IF (istatus /= 0) RETURN

  charlen = 1024
  ALLOCATE(chararr(charlen), STAT = istatus)

  DO aindex = 0, nattr-1

    istatus = sfgainfo(sd_id1,aindex,aname,dtype,nvalues)

    SELECT CASE (dtype)

    CASE ( dfnt_char8 )

      IF (TRIM(aname) == 'runname') THEN
        nvalues = LEN_TRIM(runnamein)
        DO i = 1, nvalues
          chararr(i) = runnamein(i:i)
        END DO

      ELSE
        IF (nvalues > charlen) THEN
          DEALLOCATE(chararr)
          ALLOCATE(chararr(nvalues), STAT = istatus)
        END IF
        istatus = sfrcatt(sd_id1, aindex, chararr)
      END IF

      DO iypatch = 1, nypatchout
        DO ixpatch = 1, nxpatchout
          istatus = sfscatt(fHndlOut(ixpatch,iypatch),TRIM(aname),dtype,nvalues,chararr)
        END DO
      END DO

    CASE ( dfnt_float32 )
      istatus = sfrnatt(sd_id1, aindex, buf_r)
      DO iypatch = 1, nypatchout
        DO ixpatch = 1, nxpatchout
          istatus = sfsnatt(fHndlOut(ixpatch,iypatch), TRIM(aname), dtype, nvalues, buf_r)
        END DO
      END DO
    CASE ( dfnt_int32 )
      istatus = sfrnatt(sd_id1, aindex, buf_i)

      SELECT CASE (TRIM(aname))
      CASE ( 'nx' )
        nx   = buf_i
        nxlg = (nx-3)*nprocx_in+3
        nxlgout = nxlg - nxpnt
        nxout   = (nxlgout -3) /nxpatchout + 3
        IF ( MOD((nxlgout -3), nxpatchout) /= 0) THEN
          WRITE(6,'(1x,a,/,8x,2(a,I10))')                               &
            'ERROR: The output dimension size and number of patch does not match.', &
            'nxlgout = ', nxlgout, ', nxpatchout = ',nxpatchout
          istatus = -2
          RETURN
        END IF
        buf_i   = nxout
      CASE ( 'ny' )
        ny = buf_i
        nylg = (ny-3)*nprocy_in+3
        nylgout = nylg - nypnt
        nyout   = (nylgout -3) /nypatchout + 3
        IF ( MOD((nylgout -3), nypatchout) /= 0) THEN
          WRITE(6,'(1x,a,/,8x,2(a,I10))')                               &
            'ERROR: The output dimension size and number of patch does not match.', &
            'nylgout = ', nylgout, ', nypatchout = ',nypatchout
          istatus = -2
          RETURN
        END IF
        buf_i   = nyout
      CASE ( 'nz' )
        nz = buf_i
      CASE ( 'nzsoil' )
        nzsoil = buf_i
      CASE ( 'nstyp' )
        nstyps = buf_i
      END SELECT

      DO iypatch = 1, nypatchout
        DO ixpatch = 1, nxpatchout
          istatus = sfsnatt(fHndlOut(ixpatch,iypatch), TRIM(aname), dtype, nvalues, buf_i)
        END DO
      END DO
    CASE DEFAULT
      WRITE(6,'(1x,a,I10,2a)') 'ERROR: unsupported attribute type ',    &
                               dtype, ' with attribute name ',TRIM(aname)
      istatus = -1
      RETURN
    END SELECT

  END DO

  nx = (nx-3)*nxpatch + 3     ! Change from input patch size to work size
  ny = (ny-3)*nypatch + 3

  DEALLOCATE(chararr)

  RETURN
END SUBROUTINE copy_global_attributes

!#######################################################################
SUBROUTINE peek_hdf_dataset(fHandle,ivar,vname,rank,dims,dtype,ndattr,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Get HDF file data set information
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,            INTENT(IN)  :: fHandle
  INTEGER,            INTENT(IN)  :: ivar
  CHARACTER(LEN=256), INTENT(OUT) :: vname
  INTEGER,            INTENT(OUT) :: rank
  INTEGER,            INTENT(OUT) :: dims(6)
  INTEGER,            INTENT(OUT) :: dtype
  INTEGER,            INTENT(OUT) :: ndattr
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: sds_id

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfselect
  INTEGER :: sfginfo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  sds_id = sfselect(fHandle,ivar)

  istatus = sfginfo(sds_id,vname,rank,dims,dtype,ndattr)

  RETURN
END SUBROUTINE peek_hdf_dataset

SUBROUTINE read_hdf_data_attr(sd_id,vname,comp_code,comp_prm,           &
                              comment,lcmnt,units,lunts,stag_dim,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Get HDF file data set information
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,            INTENT(IN)  :: sd_id
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,            INTENT(OUT) :: comp_code, comp_prm
  CHARACTER(LEN=256), INTENT(OUT) :: comment, units
  INTEGER,            INTENT(OUT) :: lcmnt, lunts
  INTEGER,            INTENT(OUT) :: stag_dim
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: sds_id, sds_index
  INTEGER :: aindex, dtype

  CHARACTER(LEN=356) :: aname

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfselect, sfgainfo, sfn2index, sfendacc
  INTEGER :: sffattr, sfrnatt, sfrcatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  sds_index = sfn2index(sd_id,TRIM(vname))
  IF (sds_index == -1) THEN
    WRITE(6,'(1x,3a)') 'ERROR: varaible ',TRIM(vname),' not found in read_hdf_data_attr.'
    istatus = -1
    RETURN
  END IF

  sds_id  = sfselect(sd_id,sds_index)

  aindex  = sffattr(sds_id, 'hdf_comp_code')
  istatus = sfrnatt(sds_id,aindex,comp_code)

  aindex  = sffattr(sds_id, 'hdf_comp_prm')
  istatus = sfrnatt(sds_id,aindex,comp_prm)

  aindex  = sffattr(sds_id, 'comment')
  istatus = sfgainfo(sds_id,aindex,aname,dtype,lcmnt)
  istatus = sfrcatt(sds_id,aindex,comment)

  aindex  = sffattr(sds_id, 'units')
  istatus = sfgainfo(sds_id,aindex,aname,dtype,lunts)
  istatus = sfrcatt(sds_id,aindex,units)

  aindex  = sffattr(sds_id, 'stag_dim')
  istatus = sfrnatt(sds_id,aindex,stag_dim)

  istatus = sfendacc(sds_id)

  RETURN
END SUBROUTINE read_hdf_data_attr

!#######################################################################
SUBROUTINE read_hdf_dataseti(fHandles,nxpatch,nypatch,vname,ivarin,     &
                             dims,nx,ny,nz,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Read integer array (2D or 3D) from mulitple HDF file and join.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,          INTENT(IN)  :: dims(3)
  INTEGER,          INTENT(IN)  :: nx, ny, nz
  INTEGER,          INTENT(OUT) :: ivarin(nx,ny,nz)
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER, ALLOCATABLE :: ivardta(:,:,:)
  INTEGER :: start(3), stride(3)

  INTEGER :: sds_index, sds_id

  INTEGER :: ix, jy, i, j, ia, ja, k

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  start(1)  = 0; start(2)  = 0; start(3)  = 0
  stride(1) = 1; stride(2) = 1; stride(3) = 1

  ALLOCATE(ivardta(dims(1),dims(2),dims(3)), STAT = istatus)

  DO jy = 1, nypatch
    DO ix = 1, nxpatch
      sds_index = sfn2index(fHandles(ix,jy),TRIM(vname))
      IF (sds_index == -1) THEN
        WRITE(6,'(1x,3a)') 'ERROR: varaible ',TRIM(vname),' not found in read_hdf_dataseti.'
        istatus = -1
        RETURN
      END IF

      sds_id    = sfselect(fHandles(ix,jy),sds_index)

      istatus   = sfrdata(sds_id, start, stride, dims, ivardta)
      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: Reading varaible ',TRIM(vname),' in read_hdf_dataseti, istatus = ',istatus
        RETURN
      END IF
      istatus   = sfendacc(sds_id)
      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: in read_hdf_dataseti for varaible ',TRIM(vname),', istatus = ',istatus
        RETURN
      END IF


      DO k = 1, dims(3)
        DO j = 1, dims(2)
          ja = (jy-1)*(dims(2)-3)+j
          DO i = 1, dims(1)
            ia = (ix-1)*(dims(1)-3)+i
            ivarin(ia,ja,k) = ivardta(i,j,k)
          END DO
        END DO
      END DO

    END DO
  END DO

  DEALLOCATE(ivardta)

  RETURN
END SUBROUTINE read_hdf_dataseti

!#######################################################################
SUBROUTINE write_hdf_dataseti(fHandles,nxpatch,nypatch,vname,           &
             comp_code,comp_prm,comment, lcmnt, units, lunts, stag_dim, &
             ivarlg,nxlg, nylg, ivar, nx,ny,nz,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Write integer array (2D or 3D) to mulitple HDF files. The data array
!   is split first before writing.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN):: vname
  INTEGER,          INTENT(IN)  :: comp_code, comp_prm
  CHARACTER(LEN=256), INTENT(IN):: comment, units
  INTEGER,          INTENT(IN)  :: lcmnt,lunts, stag_dim
  INTEGER,          INTENT(IN)  :: nxlg, nylg
  INTEGER,          INTENT(IN)  :: nx, ny, nz
  INTEGER,          INTENT(IN)  :: ivarlg(nxlg,nylg,nz)
  INTEGER,          INTENT(OUT) :: ivar(nx,ny,nz)
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ix, jy
  INTEGER :: i, j, k, ia, ja

  INTEGER :: ndim, dim_id, sds_id
  INTEGER :: start(3), stride(3), dims(3)

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfcreate, sfwdata, sfendacc
  INTEGER :: sfdimid, sfsdmname
  INTEGER :: sfscompress, sfsnatt, sfscatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  DO jy = 1, nypatch
    DO ix = 1, nxpatch

      DO k = 1, nz
        DO j = 1, ny
          ja = (jy-1)*(ny-3)+j
          DO i = 1, nx
            ia = (ix-1)*(nx-3)+i
            ivar(i,j,k) = ivarlg(ia,ja,k)
          END DO
        END DO
      END DO

      dims(1)   = nx;   dims(2) = ny;   dims(3) = nz
      start(:)  = 0
      stride(:) = 1

      ndim = 3
      IF (nz == 1) ndim = 2

      sds_id = sfcreate(fHandles(ix,jy),vname,dfnt_int32,ndim,dims)

      dim_id  = sfdimid(sds_id,0)
      istatus = sfsdmname(dim_id,'x_dim')

      dim_id  = sfdimid(sds_id,1)
      istatus = sfsdmname(dim_id,'y_dim')

      IF (ndim == 3) THEN
        dim_id  = sfdimid(sds_id,3)
        istatus = sfsdmname(dim_id,'nstyps')
      END IF

      istatus = sfscompress(sds_id,comp_code,comp_prm)
      istatus = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1,    comp_code)
      istatus = sfsnatt(sds_id, 'hdf_comp_prm',  dfnt_int32, 1,    comp_prm)

      istatus = sfscatt(sds_id, 'comment',       dfnt_char8, lcmnt, comment)
      istatus = sfscatt(sds_id, 'units',         dfnt_char8, lunts, units)
      istatus = sfsnatt(sds_id, 'stag_dim',      dfnt_int32, 1,     stag_dim)

      istatus = sfwdata(sds_id, start, stride, dims, ivar)
      IF (istatus /= 0) THEN
        WRITE (6,'(1x,3a)') 'ERROR: writing variable <',trim(vname),'> in write_hdf_dataseti.'
      END IF
      istatus = sfendacc(sds_id)

    END DO
  END DO

  RETURN
END SUBROUTINE write_hdf_dataseti

!#######################################################################
SUBROUTINE read_hdf_dataset(fHandles,nxpatch,nypatch,vname,varin,       &
                            dims,nx,ny,nz,packed16,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Read array (2D or 3D) from mulitple HDF file and join.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,          INTENT(IN)  :: dims(3)
  INTEGER,          INTENT(IN)  :: nx, ny, nz
  REAL,             INTENT(OUT) :: varin(nx,ny,nz)
  INTEGER,          INTENT(OUT) :: packed16
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: vardta(:,:,:)

  INTEGER :: start(3), stride(3)

  INTEGER :: sds_index, sds_id

  INTEGER :: aindex, istat1, istat2, istat3

  INTEGER :: ix, jy, i, j, ia, ja, k

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL(4)                            , ALLOCATABLE :: hmax(:), hmin(:)

  REAL    :: scalef

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc
  INTEGER :: sffattr, sfrnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  start(1)  = 0; start(2)  = 0; start(3)  = 0
  stride(1) = 1; stride(2) = 1; stride(3) = 1

  ALLOCATE(vardta(dims(1),dims(2),dims(3)), STAT = istatus)
  vardta(:,:,:) = 0.0

  DO jy = 1, nypatch
    DO ix = 1, nxpatch
      sds_index = sfn2index(fHandles(ix,jy),TRIM(vname))
      IF (sds_index == -1) THEN
        WRITE(6,'(1x,3a)') 'ERROR: varaible ',TRIM(vname),' not found in read_hdf_dataset.'
        istatus = -1
        RETURN
      END IF

      sds_id = sfselect(fHandles(ix,jy),sds_index)

      aindex = sffattr(sds_id,'packed16')
      IF (aindex > 0) THEN
        IF ( .NOT. ALLOCATED(itmp)) THEN
          ALLOCATE(itmp(dims(1),dims(2),nz), STAT = istatus)
          ALLOCATE(hmax(nz),                 STAT = istatus)
          ALLOCATE(hmin(nz),                 STAT = istatus)
        END IF
        istat1 = sfrnatt(sds_id,aindex,packed16)
        aindex = sffattr(sds_id,'max')
        istat2 = sfrnatt(sds_id,aindex,hmax)
        aindex = sffattr(sds_id,'min')
        istat3 = sfrnatt(sds_id,aindex,hmin)
        IF (istat1 == -1 .OR. istat2 == -1 .OR. istat3 == -1) THEN
          WRITE (6,'(1x,3a)') 'ERROR reading max/min for ',trim(vname),' in read_hdf_dataset.'
          istatus = -2
          RETURN
        END IF
      ELSE
        packed16 = 0
      END IF

      IF (packed16 == 0) THEN
        istatus = sfrdata(sds_id, start, stride, dims, vardta)
      ELSE
        istatus = sfrdata(sds_id, start, stride, dims, itmp)
        DO k = 1,dims(3)
          scalef = (hmax(k)-hmin(k))/65534.0
          DO j = 1, dims(2)
            DO i = 1, dims(1)
              vardta(i,j,k) = scalef * (itmp(i,j,k)+32767) + hmin(k)
            END DO
          END DO
        END DO
      END IF

      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: Reading varaible ',TRIM(vname),' in read_hdf_dataset, istatus = ',istatus
        RETURN
      END IF

      istatus   = sfendacc(sds_id)
      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: in read_hdf_dataset for varaible ',TRIM(vname),', istatus = ',istatus
        RETURN
      END IF

      DO k = 1, dims(3)
        DO j = 1, dims(2)
          ja = (jy-1)*(dims(2)-3)+j
          DO i = 1, dims(1)
            ia = (ix-1)*(dims(1)-3)+i
            varin(ia,ja,k) = vardta(i,j,k)
          END DO
        END DO
      END DO

    END DO
  END DO

  DEALLOCATE(vardta)
  IF (packed16 /= 0) DEALLOCATE(hmin, hmax, itmp)

  RETURN
END SUBROUTINE read_hdf_dataset

!#######################################################################
SUBROUTINE write_hdf_dataset(fHandles,nxpatch,nypatch,vname, packed16,  &
             comp_code,comp_prm,comment, lcmnt, units, lunts, stag_dim, &
             varlg,nxlg, nylg, var, nx,ny,nz,nzin,nzsoilin,nstypsin,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Write integer array (2D or 3D) to mulitple HDF files. The data array
!   is split first before writing.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,            INTENT(IN)  :: nxpatch, nypatch
  INTEGER,            INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,            INTENT(IN)  :: packed16
  INTEGER,            INTENT(IN)  :: comp_code, comp_prm
  CHARACTER(LEN=256), INTENT(IN)  :: comment, units
  INTEGER,            INTENT(IN)  :: lcmnt,lunts, stag_dim
  INTEGER,            INTENT(IN)  :: nxlg, nylg
  INTEGER,            INTENT(IN)  :: nx, ny, nz
  REAL,               INTENT(IN)  :: varlg(nxlg,nylg,nz)
  REAL,               INTENT(OUT) :: var(nx,ny,nz)
  INTEGER,            INTENT(IN)  :: nzin, nzsoilin, nstypsin
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ix, jy
  INTEGER :: i, j, k, ia, ja

  INTEGER :: ndim, dim_id, sds_id
  INTEGER :: start(3), stride(3), dims(3)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL(4),                             ALLOCATABLE :: hmax(:), hmin(:)

  REAL    :: scalef

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfcreate, sfwdata, sfendacc
  INTEGER :: sfdimid, sfsdmname
  INTEGER :: sfscompress, sfsnatt, sfscatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (packed16 /= 0) THEN
    ALLOCATE(itmp(nx,ny,nz), STAT = istatus)
    ALLOCATE(hmin(nz),       STAT = istatus)
    ALLOCATE(hmax(nz),       STAT = istatus)
  END IF

  DO jy = 1, nypatch
    DO ix = 1, nxpatch

      DO k = 1, nz
        DO j = 1, ny
          ja = (jy-1)*(ny-3)+j
          DO i = 1, nx
            ia = (ix-1)*(nx-3)+i
            var(i,j,k) = varlg(ia,ja,k)
          END DO
        END DO
      END DO

!--------------------- Define dataset in the file ----------------------

      dims(1)   = nx;   dims(2) = ny;   dims(3) = nz
      start(:)  = 0
      stride(:) = 1

      ndim = 3
      IF (nz == 1) ndim = 2

      IF (packed16 /= 0) THEN
        sds_id = sfcreate(fHandles(ix,jy),vname,dfnt_int16,ndim,dims)
      ELSE
        sds_id = sfcreate(fHandles(ix,jy),vname,dfnt_float32,ndim,dims)
      END IF

!--------------------- Change dimension names ----------------------

      dim_id  = sfdimid(sds_id,0)
      istatus = sfsdmname(dim_id,'x_dim')

      dim_id  = sfdimid(sds_id,1)
      istatus = sfsdmname(dim_id,'y_dim')

      IF (ndim == 3) THEN
        dim_id  = sfdimid(sds_id,2)
        IF (nz == nzin) THEN
          istatus = sfsdmname(dim_id,'z_dim')
        ELSE IF (nz == nzsoilin) THEN
          istatus = sfsdmname(dim_id,'soil_dim')
        ELSE IF (nz == nstypsin) THEN
          istatus = sfsdmname(dim_id,'nstyps')
        ELSE IF (nz == nstypsin+1) THEN
          istatus = sfsdmname(dim_id,'nstyps_dim')
        ELSE
          WRITE(6,'(1x,3a)')                                            &
            'ERROR: unsupported dimension (3rd dimension) for ',        &
            trim(vname),' in write_hdf_dataset.'
          istatus = -3
          RETURN
        END IF
      END IF

!--------------------- Define dataset attributes ----------------------

      istatus = sfscompress(sds_id,comp_code,comp_prm)
      istatus = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1,    comp_code)
      istatus = sfsnatt(sds_id, 'hdf_comp_prm',  dfnt_int32, 1,    comp_prm)

      istatus = sfscatt(sds_id, 'comment',       dfnt_char8, lcmnt, comment)
      istatus = sfscatt(sds_id, 'units',         dfnt_char8, lunts, units)
      istatus = sfsnatt(sds_id, 'stag_dim',      dfnt_int32, 1,     stag_dim)

!--------------------- Packed and write dataset ----------------------

      IF (packed16 /= 0) THEN
        DO k = 1,nz
          hmax(k) = MAXVAL(var(:,:,k))
          hmin(k) = MINVAL(var(:,:,k))

          IF (ABS(hmax(k)) < 1.0E-25) hmax(k)= 0.0  !Added by Xue and Dawson
          IF (ABS(hmin(k)) < 1.0E-25) hmin(k)= 0.0  !Added by "" and ""

          IF (ABS(hmax(k)-hmin(k)) > 1.0E-10) THEN
            scalef = 65534.0 / (hmax(k) - hmin(k))
          ELSE
            scalef = 65534.0
          END IF

          DO j=1,ny
            DO i=1,nx
              itmp(i,j,k) = nint(scalef * (var(i,j,k) - hmin(k))) - 32767
            END DO
          END DO
        END DO
        istatus = sfsnatt(sds_id, 'packed16', dfnt_int32,    1,    1)
        istatus = sfsnatt(sds_id, 'max',      dfnt_float32, nz, hmax)
        istatus = sfsnatt(sds_id, 'min',      dfnt_float32, nz, hmin)

        istatus = sfwdata(sds_id, start, stride, dims, itmp)
      ELSE
        istatus = sfwdata(sds_id, start, stride, dims, var)
      END IF

      IF (istatus /= 0) THEN
        WRITE (6,'(1x,3a)') 'ERROR: writing variable <',trim(vname),'> in write_hdf_dataset.'
      END IF
      istatus = sfendacc(sds_id)

    END DO
  END DO

  RETURN
END SUBROUTINE write_hdf_dataset

!#######################################################################
SUBROUTINE read_hdf_dataset1d(fHandles,nxpatch,nypatch,vname,varin,     &
                              dims,direction,ni,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Read array (2D or 3D) from mulitple HDF file and join.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,          INTENT(IN)  :: dims(3)
  INTEGER,          INTENT(IN)  :: ni
  CHARACTER(LEN=1), INTENT(IN)  :: direction
  REAL,             INTENT(OUT) :: varin(ni)
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: vardta(:)

  INTEGER :: start, stride

  INTEGER :: sds_index, sds_id

  INTEGER :: ix, jy, i, ia


!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc
  INTEGER :: sffattr, sfrnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  start  = 0
  stride = 1

  ALLOCATE(vardta(dims(1)), STAT = istatus)

  DO jy = 1, nypatch
    DO ix = 1, nxpatch
      sds_index = sfn2index(fHandles(ix,jy),TRIM(vname))
      IF (sds_index == -1) THEN
        WRITE(6,'(1x,3a)') 'ERROR: varaible ',TRIM(vname),' not found in read_hdf_dataset1d.'
        istatus = -1
        RETURN
      END IF

      sds_id = sfselect(fHandles(ix,jy),sds_index)

      istatus = sfrdata(sds_id, start, stride, dims, vardta)

      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: Reading varaible ',TRIM(vname),' in read_hdf_dataset, istatus = ',istatus
        RETURN
      END IF

      istatus   = sfendacc(sds_id)
      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: in read_hdf_dataset for varaible ',TRIM(vname),', istatus = ',istatus
        RETURN
      END IF

      IF (direction == 'x' .AND. jy == 1) THEN
        DO i = 1, dims(1)
          ia = (ix-1)*(dims(1)-3)+i
          varin(ia) = vardta(i)
        END DO
      ELSE IF (direction == 'y' .AND. ix == 1) THEN
        DO i = 1, dims(1)
          ia = (jy-1)*(dims(1)-3)+i
          varin(ia) = vardta(i)
        END DO
      ELSE IF (direction == 'z' .AND. ix == 1 .AND. jy == 1) THEN
        DO i = 1, dims(1)
          varin(i) = vardta(i)
        END DO
      END IF

    END DO
  END DO

  DEALLOCATE(vardta)

  RETURN
END SUBROUTINE read_hdf_dataset1d

!#######################################################################
SUBROUTINE write_hdf_dataset1d(fHandles,nxpatch,nypatch,vname,          &
             comment, lcmnt, units, lunts, direction,                   &
             varlg, nilg, var, ni, istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Write integer array (2D or 3D) to mulitple HDF files. The data array
!   is split first before writing.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  CHARACTER(LEN=256), INTENT(IN)  :: comment, units
  CHARACTER(LEN=1), INTENT(IN)  :: direction
  INTEGER,          INTENT(IN)  :: lcmnt,lunts
  INTEGER,          INTENT(IN)  :: nilg
  INTEGER,          INTENT(IN)  :: ni
  REAL,             INTENT(IN)  :: varlg(nilg)
  REAL,             INTENT(OUT) :: var(ni)
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ix, jy
  INTEGER :: i, ia

  INTEGER :: dim_id, sds_id
  INTEGER :: start, stride, dims

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfcreate, sfwdata, sfendacc
  INTEGER :: sfdimid, sfsdmname
  INTEGER :: sfscompress, sfsnatt, sfscatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  DO jy = 1, nypatch
    DO ix = 1, nxpatch

      IF (direction == 'x') THEN
        DO i = 1, ni
          ia = (ix-1)*(ni-3)+i
          var(i) = varlg(ia)
        END DO
      ELSE IF (direction == 'y') THEN
        DO i = 1, ni
          ia = (jy-1)*(ni-3)+i
          var(i) = varlg(ia)
        END DO
      ELSE
        DO i = 1, ni
          var(i) = varlg(i)
        END DO
      END IF

!--------------------- Define dataset in the file ----------------------

      dims   = ni
      start  = 0
      stride = 1

      sds_id = sfcreate(fHandles(ix,jy),vname,dfnt_float32,1,dims)

!--------------------- Change dimension names ----------------------

      dim_id  = sfdimid(sds_id,0)
      istatus = sfsdmname(dim_id,direction//'_dim')

!--------------------- Define dataset attributes ----------------------

      istatus = sfscatt(sds_id, 'comment',       dfnt_char8, lcmnt, comment)
      istatus = sfscatt(sds_id, 'units',         dfnt_char8, lunts, units)

!--------------------- Packed and write dataset ----------------------

      istatus = sfwdata(sds_id, start, stride, dims, var)

      IF (istatus /= 0) THEN
        WRITE (6,'(1x,3a)') 'ERROR: writing variable <',trim(vname),'> in write_hdf_dataset1d.'
      END IF
      istatus = sfendacc(sds_id)

    END DO
  END DO

  RETURN
END SUBROUTINE write_hdf_dataset1d

!#######################################################################
SUBROUTINE read_hdf_dataset4d(fHandles,nxpatch,nypatch,vname,varin,     &
                            dims,nx,ny,nzsoil,nstyp,packed16,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Read 4d array from mulitple HDF file and join.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,          INTENT(IN)  :: dims(4)
  INTEGER,          INTENT(IN)  :: nx, ny, nzsoil,nstyp
  REAL,             INTENT(OUT) :: varin(nx,ny,nzsoil,nstyp)
  INTEGER,          INTENT(OUT) :: packed16
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: vardta(:,:,:,:)

  INTEGER :: start(4), stride(4)

  INTEGER :: sds_index, sds_id

  INTEGER :: aindex, istat1, istat2, istat3

  INTEGER :: ix, jy, i, j, ia, ja, k, n

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:,:)
  REAL(4)                            , ALLOCATABLE :: hmax(:), hmin(:)

  REAL    :: scalef

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc
  INTEGER :: sffattr, sfrnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  start(1)  = 0; start(2)  = 0; start(3)  = 0; start(4)  = 0
  stride(1) = 1; stride(2) = 1; stride(3) = 1; stride(4) = 1

  ALLOCATE(vardta(dims(1),dims(2),dims(3),dims(4)), STAT = istatus)

  DO jy = 1, nypatch
    DO ix = 1, nxpatch
      sds_index = sfn2index(fHandles(ix,jy),TRIM(vname))
      IF (sds_index == -1) THEN
        WRITE(6,'(1x,3a)') 'ERROR: varaible ',TRIM(vname),' not found in read_hdf_dataset4d.'
        istatus = -1
        RETURN
      END IF

      sds_id = sfselect(fHandles(ix,jy),sds_index)

      aindex = sffattr(sds_id,'packed16')
      IF (aindex > 0) THEN
        IF ( .NOT. ALLOCATED(itmp)) THEN
          ALLOCATE(itmp(dims(1),dims(2),dims(3),dims(4)), STAT = istatus)
          ALLOCATE(hmax(dims(3)),                 STAT = istatus)
          ALLOCATE(hmin(dims(3)),                 STAT = istatus)
        END IF
        istat1 = sfrnatt(sds_id,aindex,packed16)
        aindex = sffattr(sds_id,'max')
        istat2 = sfrnatt(sds_id,aindex,hmax)
        aindex = sffattr(sds_id,'min')
        istat3 = sfrnatt(sds_id,aindex,hmin)
        IF (istat1 == -1 .OR. istat2 == -1 .OR. istat3 == -1) THEN
          WRITE (6,'(1x,3a)') 'ERROR reading max/min for ',trim(vname),' in read_hdf_dataset.'
          istatus = -2
          RETURN
        END IF
      ELSE
        packed16 = 0
      END IF

      IF (packed16 == 0) THEN
        istatus = sfrdata(sds_id, start, stride, dims, vardta)
      ELSE
        istatus = sfrdata(sds_id, start, stride, dims, itmp)
        DO n = 1, nstyp
          DO k = 1,dims(3)
            scalef = (hmax(k)-hmin(k))/65534.0
            DO j = 1, dims(2)
              DO i = 1, dims(1)
                vardta(i,j,k,n) = scalef * (itmp(i,j,k,n)+32767) + hmin(k)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: Reading varaible ',TRIM(vname),' in read_hdf_dataset, istatus = ',istatus
        RETURN
      END IF

      istatus   = sfendacc(sds_id)
      IF (istatus /= 0) THEN
        WRITE(6,'(1x,3a,I4)') 'ERROR: in read_hdf_dataset for varaible ',TRIM(vname),', istatus = ',istatus
        RETURN
      END IF

      DO n = 1, dims(4)
        DO k = 1, dims(3)
          DO j = 1, dims(2)
            ja = (jy-1)*(dims(2)-3)+j
            DO i = 1, dims(1)
              ia = (ix-1)*(dims(1)-3)+i
              varin(ia,ja,k,n) = vardta(i,j,k,n)
            END DO
          END DO
        END DO
      END DO

    END DO
  END DO

  DEALLOCATE(vardta)
  IF (packed16 /= 0) DEALLOCATE(hmin, hmax, itmp)

  RETURN
END SUBROUTINE read_hdf_dataset4d

!#######################################################################
SUBROUTINE write_hdf_dataset4d(fHandles,nxpatch,nypatch,vname, packed16,&
             comp_code,comp_prm,comment, lcmnt, units, lunts, stag_dim, &
             varlg,nxlg,nylg, var,nx,ny,nzsoil,nstyp,istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Write integer array (2D or 3D) to mulitple HDF files. The data array
!   is split first before writing.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nxpatch, nypatch
  INTEGER,          INTENT(IN)  :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN)  :: vname
  INTEGER,          INTENT(IN)  :: packed16
  INTEGER,          INTENT(IN)  :: comp_code, comp_prm
  CHARACTER(LEN=256), INTENT(IN)  :: comment, units
  INTEGER,          INTENT(IN)  :: lcmnt,lunts, stag_dim
  INTEGER,          INTENT(IN)  :: nxlg, nylg
  INTEGER,          INTENT(IN)  :: nx, ny, nzsoil, nstyp
  REAL,             INTENT(IN)  :: varlg(nxlg,nylg,nzsoil,nstyp)
  REAL,             INTENT(OUT) :: var(nx,ny,nzsoil,nstyp)

  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ix, jy
  INTEGER :: i, j, k, n, ia, ja

  INTEGER :: ndim, dim_id, sds_id
  INTEGER :: start(4), stride(4), dims(4)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:,:)
  REAL(4),                             ALLOCATABLE :: hmax(:), hmin(:)

  REAL    :: scalef

!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

  INTEGER :: sfcreate, sfwdata, sfendacc
  INTEGER :: sfdimid, sfsdmname
  INTEGER :: sfscompress, sfsnatt, sfscatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (packed16 /= 0) THEN
    ALLOCATE(itmp(nx,ny,nzsoil,nstyp), STAT = istatus)
    ALLOCATE(hmin(nzsoil),       STAT = istatus)
    ALLOCATE(hmax(nzsoil),       STAT = istatus)
  END IF

  DO jy = 1, nypatch
    DO ix = 1, nxpatch

      DO n = 1,nstyp
        DO k = 1, nzsoil
          DO j = 1, ny
            ja = (jy-1)*(ny-3)+j
            DO i = 1, nx
              ia = (ix-1)*(nx-3)+i
              var(i,j,k,n) = varlg(ia,ja,k,n)
            END DO
          END DO
        END DO
      END DO

!--------------------- Define dataset in the file ----------------------

      dims(1)   = nx; dims(2) = ny; dims(3) = nzsoil; dims(4) = nstyp
      start(:)  = 0
      stride(:) = 1

      ndim = 4

      IF (packed16 /= 0) THEN
        sds_id = sfcreate(fHandles(ix,jy),vname,dfnt_int16,ndim,dims)
      ELSE
        sds_id = sfcreate(fHandles(ix,jy),vname,dfnt_float32,ndim,dims)
      END IF

!--------------------- Change dimension names ----------------------

      dim_id  = sfdimid(sds_id,0)
      istatus = sfsdmname(dim_id,'x_dim')

      dim_id  = sfdimid(sds_id,1)
      istatus = sfsdmname(dim_id,'y_dim')

      dim_id  = sfdimid(sds_id,2)
      istatus = sfsdmname(dim_id,'soil_dim')

      dim_id  = sfdimid(sds_id,3)
      istatus = sfsdmname(dim_id,'nstyps_dim')

!--------------------- Define dataset attributes ----------------------

      istatus = sfscompress(sds_id,comp_code,comp_prm)
      istatus = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1,    comp_code)
      istatus = sfsnatt(sds_id, 'hdf_comp_prm',  dfnt_int32, 1,    comp_prm)

      istatus = sfscatt(sds_id, 'comment',       dfnt_char8, lcmnt, comment)
      istatus = sfscatt(sds_id, 'units',         dfnt_char8, lunts, units)
      istatus = sfsnatt(sds_id, 'stag_dim',      dfnt_int32, 1,     stag_dim)

!--------------------- Packed and write dataset ----------------------

      IF (packed16 /= 0) THEN
        DO n = 1,nstyp
          DO k = 1,nzsoil
            hmax(k) = MAXVAL(var(:,:,k,:))
            hmin(k) = MINVAL(var(:,:,k,:))

            IF (ABS(hmax(k)) < 1.0E-25) hmax(k)= 0.0  !Added by Xue and Dawson
            IF (ABS(hmin(k)) < 1.0E-25) hmin(k)= 0.0  !Added by "" and ""

            IF (ABS(hmax(k)-hmin(k)) > 1.0E-10) THEN
              scalef = 65534.0 / (hmax(k) - hmin(k))
            ELSE
              scalef = 65534.0
            END IF

            DO j=1,ny
              DO i=1,nx
                itmp(i,j,k,n) = nint(scalef * (var(i,j,k,n) - hmin(k))) - 32767
              END DO
            END DO
          END DO
        END DO
        istatus = sfsnatt(sds_id, 'packed16', dfnt_int32,    1,    1)
        istatus = sfsnatt(sds_id, 'max',      dfnt_float32, nzsoil, hmax)
        istatus = sfsnatt(sds_id, 'min',      dfnt_float32, nzsoil, hmin)

        istatus = sfwdata(sds_id, start, stride, dims, itmp)
      ELSE
        istatus = sfwdata(sds_id, start, stride, dims, var)
      END IF

      IF (istatus /= 0) THEN
        WRITE (6,'(1x,3a)') 'ERROR: writing variable <',trim(vname),'> in write_hdf_dataset.'
      END IF
      istatus = sfendacc(sds_id)

    END DO
  END DO

  RETURN
END SUBROUTINE write_hdf_dataset4d
