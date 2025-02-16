!
SUBROUTINE  splitncdf(filenames,nfile,stag, xdimname, ydimname,         &
                      nproc_x,nproc_y,outdirname, debug,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Split files in netCDF format into patches. The patched files will
!    contain the same data as original file but in evenly divided
!    subdomain specified by the user.
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang (11/07/2006)
!
! MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! Variable declaration
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)            :: nfile
  CHARACTER(LEN=*), INTENT(IN)   :: filenames(nfile)
  CHARACTER(LEN=*), INTENT(IN)   :: outdirname
  CHARACTER(LEN=*), INTENT(IN)   :: xdimname, ydimname
  LOGICAL,          INTENT(IN)   :: stag               ! Whether the data contains
       ! both staggered and unstaggered dimensions. If it is true, then xdimname/ydimname is
       ! the name of staggered dimension. The unstaggered dimensions are assumed to be
       ! xdimname/ydimname by tripping the trail '_stag'. Such as xdimname = 'x_stag' then
       ! the unstaggered dimension name is 'x'.

  INTEGER, INTENT(IN)            :: nproc_x, nproc_y
  INTEGER, INTENT(IN)            :: debug
  INTEGER, INTENT(OUT)           :: istatus

!-----------------------------------------------------------------------
!
! Dimensions and work arrays
!
!-----------------------------------------------------------------------

  INTEGER :: nxlg, nylg
  INTEGER :: nx,   ny

  INTEGER, ALLOCATABLE :: variin(:), variout(:)
  INTEGER, ALLOCATABLE :: varain(:), varaout(:)

  INTEGER :: finid
  CHARACTER(LEN=256), ALLOCATABLE :: outfilenames(:,:)
  INTEGER,            ALLOCATABLE :: foutids(:,:)
  INTEGER,            ALLOCATABLE :: dimouts0(:,:,:)
!
!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_RANK = 5    ! Assume the max rank is 5

  INTEGER :: nf, dindx
  INTEGER :: iout, iin

  CHARACTER(LEN=256) :: xsdimname, ysdimname
  INTEGER :: nxid, nyid, nxsid, nysid
  INTEGER :: iloc, jloc

  INTEGER :: ndims, unlimdimid, ngatts, nvars
  INTEGER :: dimid, odimid, attnum, varid, ovarid
  INTEGER :: dimlen
  CHARACTER(LEN=256) :: dimname, attname, varname

  INTEGER :: vartype, varndims, varnatts

  INTEGER :: vardim, varainsize, varaoutsize, varaallsize, varballsize

  INTEGER :: vardimids(MAX_RANK), startidx(MAX_RANK), countidx(MAX_RANK), outidx(MAX_RANK)
  INTEGER :: sin1d, sin2d, sin3d, sin4d
  INTEGER :: sout1d, sout2d, sout3d, sout4d
  INTEGER :: nd1, nd2, nd3, nd4, nd5      ! Assume the max rank is 5
!
!-----------------------------------------------------------------------
!
! Including files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER :: dimina(NF_MAX_DIMS)         ! Dimension size in original file
  INTEGER :: dimouta(NF_MAX_DIMS)        ! Dimension size in split files

  CHARACTER(LEN=256) :: diminnames(NF_MAX_DIMS)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxsid = 0
  nysid = 0
  xsdimname = ' '
  ysdimname = ' '

  IF (stag) THEN
    xsdimname = xdimname(1:INDEX(xdimname,'_stag')-1)
    ysdimname = ydimname(1:INDEX(ydimname,'_stag')-1)
  END IF

!-----------------------------------------------------------------------
!
! Check dimensions first
!
!-----------------------------------------------------------------------

  istatus = nf_open(filenames(1),NF_NOWRITE,finid)
  IF (istatus /= NF_NOERR) CALL handle_err(istatus)

  istatus = nf_inq_dimid(finid,xdimname,nxid)
  IF (istatus /= NF_NOERR) CALL handle_err(istatus)

  istatus = nf_inq_dimlen(finid,nxid,nxlg)
  IF (istatus /= NF_NOERR) CALL handle_err(istatus)

  istatus = nf_inq_dimid(finid,ydimname,nyid)
  IF (istatus /= NF_NOERR) CALL handle_err(istatus)

  istatus = nf_inq_dimlen(finid,nyid,nylg)
  IF (istatus /= NF_NOERR) CALL handle_err(istatus)

  WRITE(6,'(1x,2(a,I5))') 'Dimensions in file to be split are: nx = ',nxlg,', ny = ',nylg

  IF (MOD((nxlg-3),nproc_x) /= 0 ) THEN
    WRITE(6,'(1x,a)') 'ERROR: Wrong dimension size'
    WRITE(6,'(1x,2(a,I5),a)') '       Dimension size in X direction (',nxlg,&
                           ') is not divisible by nproc_x (',nproc_x,').'
    istatus = -1
    RETURN
  END IF

  IF (MOD((nylg-3),nproc_y) /= 0 ) THEN
    WRITE(6,'(1x,a)') 'ERROR: Wrong dimension size'
    WRITE(6,'(1x,2(a,I5),a)') '       Dimension size in Y direction (',nylg,&
                           ') is not divisible by nproc_y (',nproc_y,').'
    istatus = -2
    RETURN
  END IF

  istatus = nf_close(finid)
  IF (istatus /= NF_NOERR) CALL handle_err(istatus)

!-----------------------------------------------------------------------
!
!  Loop over filenames
!
!-----------------------------------------------------------------------

  startidx(:) = 1

  ALLOCATE(outfilenames(nproc_x,nproc_y), STAT = istatus)
  ALLOCATE(foutids(nproc_x,nproc_y),      STAT = istatus)
  ALLOCATE(dimouts0(0:NF_MAX_DIMS,nproc_x,nproc_y),   STAT = istatus)

  DO nf = 1,nfile

    IF (debug > 0) WRITE(6,'(1x,2a)') 'Opening file - ',filenames(nf)
    istatus = nf_open(filenames(nf),NF_NOWRITE,finid)   ! Open file
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    dindx = INDEX(filenames(nf),'/',.TRUE.) + 1
    DO jloc = 1,nproc_y                   ! Create patches
      DO iloc = 1,nproc_x
        CALL gtsplitfn(TRIM(outdirname)//TRIM(filenames(nf)(dindx:)),   &
                       nproc_x,nproc_y,1,1,iloc,jloc,                   &
                       0,0,0,2,outfilenames(iloc,jloc),istatus)

        IF (debug > 0) WRITE(6,'(1x,2a)') 'Creating file - ',TRIM(outfilenames(iloc,jloc))
        istatus = nf_create(outfilenames(iloc,jloc),NF_CLOBBER,foutids(iloc,jloc))
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO
    END DO

    !
    ! Set dimensions
    !
    istatus = nf_inq_dimid(finid,xdimname,nxid)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    istatus = nf_inq_dimlen(finid,nxid,nxlg)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    nx = (nxlg-3)/nproc_x + 3

    istatus = nf_inq_dimid(finid,ydimname,nyid)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    istatus = nf_inq_dimlen(finid,nyid,nylg)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    ny = (nylg-3)/nproc_y + 3

    IF (stag) THEN
      istatus = nf_inq_dimid(finid,xsdimname,nxsid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_dimlen(finid,nxsid,dimlen)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (dimlen /= nxlg -1) THEN
        WRITE(6,'(1x,3a,I4,a,I4,a)') 'ERROR: Wrong size with dimension - ', &
                     TRIM(xsdimname),'(',dimlen,'), expected: ',nxlg-1,'.'
        istatus = -5
        RETURN
      END IF

      istatus = nf_inq_dimid(finid,ysdimname,nysid)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      istatus = nf_inq_dimlen(finid,nysid,dimlen)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (dimlen /= nylg -1) THEN
        WRITE(6,'(1x,3a,I4,a,I4,a)') 'ERROR: Wrong size with dimension - ', &
                     TRIM(ysdimname),'(',dimlen,'), expected: ',nylg-1,'.'
        istatus = -5
        RETURN
      END IF

    END IF

    istatus = nf_inq_ndims(finid,ndims)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    istatus = nf_inq_unlimdim(finid,unlimdimid)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)

    IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Copying dimensions - ',ndims
    DO dimid = 1,ndims
      istatus = nf_inq_dim(finid,dimid,dimname,dimlen)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      diminnames(dimid) = dimname
      dimina(dimid)  = dimlen             ! Save dimension id and len
      dimouta(dimid) = dimlen             ! Output dimension id and len
      IF (dimid == nxid) THEN
        dimlen = nx
        dimouta(dimid) = dimlen
      ELSE IF (dimid == nxsid) THEN
        dimlen = nx-1
        dimouta(dimid) = dimlen
      ELSE IF (dimid == nyid) THEN
        dimlen = ny
        dimouta(dimid) = dimlen
      ELSE IF (dimid == nysid) THEN
        dimlen = ny-1
        dimouta(dimid) = dimlen
      ELSE IF (dimid == unlimdimid) THEN
        dimlen = NF_UNLIMITED
      END IF

      IF (debug > 0) WRITE(6,'(9x,2a)') 'Dimension name - ',TRIM(dimname)
      DO jloc = 1,nproc_y                   ! Write patches Dimensions
        DO iloc = 1,nproc_x
           istatus = nf_def_dim(foutids(iloc,jloc),dimname,dimlen,odimid)
           IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        END DO
      END DO
    END DO

    dimouts0(:,:,:) = 0         ! The starting index of each dimensions for each subdomain
    DO jloc = 1,nproc_y
      DO iloc = 1,nproc_x
        dimouts0(nxid,iloc,jloc) = (iloc-1)*(nx-3)
        dimouts0(nyid,iloc,jloc) = (jloc-1)*(ny-3)
      END DO
    END DO
    IF (stag) THEN
       dimouts0(nxsid,:,:) = dimouts0(nxid,:,:)
       dimouts0(nysid,:,:) = dimouts0(nyid,:,:)
    END IF

    !
    ! Set Global attributes
    !
    istatus = nf_inq_natts(finid,ngatts)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Copying global attributes - ',ngatts
    DO attnum = 1,ngatts
      istatus = nf_inq_attname(finid,NF_GLOBAL,attnum,attname)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (debug > 0) WRITE(6,'(9x,2a)') 'Attribute name - ',TRIM(attname)
      DO jloc = 1,nproc_y                   ! Write patches global attributes
        DO iloc = 1,nproc_x
           istatus = nf_copy_att(finid,NF_GLOBAL,attname,foutids(iloc,jloc),NF_GLOBAL)
           IF (istatus /= NF_NOERR) CALL handle_err(istatus)
        END DO
      END DO

    END DO

    !
    ! Define variables
    !
    istatus = nf_inq_nvars(finid,nvars)
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)
    IF (debug > 0) WRITE(6,'(5x,a,I2)') 'Defining variables - ',nvars

    DO varid = 1,nvars
      istatus = nf_inq_var(finid,varid,varname,vartype,varndims,vardimids,varnatts)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (debug > 0) WRITE(6,'(9x,2a)') 'Variables - ',TRIM(varname)
      DO jloc = 1,nproc_y                  ! Define patches variables
        DO iloc = 1,nproc_x
           ! Dimensions should be in the same order
           istatus = nf_def_var(foutids(iloc,jloc),varname,vartype,varndims,vardimids,ovarid)
           IF (istatus /= NF_NOERR) CALL handle_err(istatus)
           DO attnum = 1,varnatts          ! Copy variable attributes
             istatus = nf_inq_attname(finid,varid,attnum,attname)
             IF (istatus /= NF_NOERR) CALL handle_err(istatus)

             istatus = nf_copy_att(finid,varid,attname,foutids(iloc,jloc),ovarid)
             IF (istatus /= NF_NOERR) CALL handle_err(istatus)
           END DO
        END DO
      END DO

    END DO

    DO jloc = 1,nproc_y                   ! End patches DEF mode
      DO iloc = 1,nproc_x
         istatus = nf_enddef(foutids(iloc,jloc))
         IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO
    END DO
    IF(debug > 0) WRITE(6,'(1x,a)') 'All patches have been defined.'

    !
    ! Split and assign variables, LOOP over variables
    !

    DO varid = 1,nvars
      vardimids(:) = 0
      istatus = nf_inq_var(finid,varid,varname,vartype,varndims,vardimids,varnatts)
      IF (istatus /= NF_NOERR) CALL handle_err(istatus)

      IF (debug >0) WRITE(6,'(9x,2a)') 'Processing variables - ',TRIM(varname)

      varainsize  = 1
      varaoutsize = 1
      outidx(:)   = 1            ! Initial all dimension size to be 1
      countidx(:) = 1
      DO vardim = 1, varndims
        countidx(vardim) = dimina(vardimids(vardim))
        outidx(vardim)   = dimouta(vardimids(vardim))
        varainsize       = countidx(vardim)*varainsize
        varaoutsize      = outidx(vardim)*varaoutsize

        IF (debug > 1) WRITE(6,'(13x,a15,2(a,I4),a)')   &
          TRIM(diminnames(vardimids(vardim))),': ',countidx(vardim),'(in) - ',&
                                                   outidx(vardim),  '(out).'
      END DO

      sin1d = countidx(1)          ! size of one column,  INPUT
      sin2d = sin1d*countidx(2)    ! size of one slice (xy)
      sin3d = sin2d*countidx(3)    ! size of one cell  (xyz)
      sin4d = sin3d*countidx(4)

      sout1d = outidx(1)           ! size of one column, OUTPUT
      sout2d = sout1d*outidx(2)
      sout3d = sout2d*outidx(3)
      sout4d = sout3d*outidx(4)

      SELECT CASE (vartype)
      CASE (NF_INT)

        IF (varainsize > varaallsize) THEN   ! Allocate input array only when necessary
          IF (ALLOCATED(variin)) DEALLOCATE(variin, STAT = istatus)
          ALLOCATE(variin(varainsize), STAT = istatus)
          varaallsize = varainsize
        END IF

        IF (varaoutsize > varballsize) THEN  ! Allocate output array only when necessary
          IF (ALLOCATED(variout)) DEALLOCATE(variout, STAT = istatus)
          ALLOCATE(variout(varaoutsize), STAT = istatus)
          varballsize = varaoutsize
        END IF

        istatus = NF_GET_VARA_INT(finid,varid,startidx,countidx,variin)

        DO jloc = 1,nproc_y                  ! Write patches variables
          DO iloc = 1,nproc_x

             IF (debug > 1) THEN
               WRITE(6,FMT='(13x,a,2(I4,a))',ADVANCE='NO')                  &
                           'Writing to processor - (',iloc,',',jloc,') with: '
               DO vardim = 1,varndims
                 WRITE(6,FMT='(I4,a)',ADVANCE='NO') dimouts0(vardimids(vardim),iloc,jloc),' '
               END DO
               WRITE(6,*)
             END IF

             DO nd5 = 1, outidx(5)        ! Assume max rank is 5, IMPORTANT
             DO nd4 = 1, outidx(4)
             DO nd3 = 1, outidx(3)
             DO nd2 = 1, outidx(2)
             DO nd1 = 1, outidx(1)
               iin =  nd1+dimouts0(vardimids(1),iloc,jloc)              &
                   + (nd2+dimouts0(vardimids(2),iloc,jloc)-1)*sin1d     &
                   + (nd3+dimouts0(vardimids(3),iloc,jloc)-1)*sin2d     &
                   + (nd4+dimouts0(vardimids(4),iloc,jloc)-1)*sin3d     &
                   + (nd5+dimouts0(vardimids(5),iloc,jloc)-1)*sin4d

               iout = nd1 + (nd2-1)*sout1d + (nd3-1)*sout2d             &
                          + (nd4-1)*sout3d + (nd5-1)*sout4d

!               IF (debug > 2) WRITE(6,'(13x,a,2I2,3(a,I4),4I4)')        &
!                  'Processor - ',iloc,jloc,': Extracting from ',iin,    &
!                  ' to ',iout,' at ',nd1,nd2,nd3,nd4,nd5

               variout(iout) = variin(iin)
             END DO
             END DO
             END DO
             END DO
             END DO

             istatus = nf_put_vara_INT(foutids(iloc,jloc),varid,startidx,outidx,variout)
             IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          END DO
        END DO

      CASE (NF_FLOAT)
        IF (varainsize > varaallsize) THEN   ! Allocate input array only when necessary
          IF (ALLOCATED(varain)) DEALLOCATE(varain, STAT = istatus)
          ALLOCATE(varain(varainsize), STAT = istatus)
          varaallsize = varainsize
        END IF

        IF (varaoutsize > varballsize) THEN  ! Allocate output array only when necessary
          IF (ALLOCATED(varaout)) DEALLOCATE(varaout, STAT = istatus)
          ALLOCATE(varaout(varaoutsize), STAT = istatus)
          varballsize = varaoutsize
        END IF

        istatus = NF_GET_VARA_REAL(finid,varid,startidx,countidx,varain)

        DO jloc = 1,nproc_y                  ! Write patches variables
          DO iloc = 1,nproc_x

             IF (debug > 1) THEN
               WRITE(6,FMT='(13x,a,2(I4,a))',ADVANCE='NO')                  &
                           'Writing to processor - (',iloc,',',jloc,') with: '
               DO vardim = 1,varndims
                 WRITE(6,FMT='(I4,a)',ADVANCE='NO') dimouts0(vardimids(vardim),iloc,jloc),' '
               END DO
               WRITE(6,*)
             END IF

             DO nd5 = 1, outidx(5)        ! Assume max rank is 5, IMPORTANT
             DO nd4 = 1, outidx(4)
             DO nd3 = 1, outidx(3)
             DO nd2 = 1, outidx(2)
             DO nd1 = 1, outidx(1)
               iin =  nd1+dimouts0(vardimids(1),iloc,jloc)              &
                   + (nd2+dimouts0(vardimids(2),iloc,jloc)-1)*sin1d     &
                   + (nd3+dimouts0(vardimids(3),iloc,jloc)-1)*sin2d     &
                   + (nd4+dimouts0(vardimids(4),iloc,jloc)-1)*sin3d     &
                   + (nd5+dimouts0(vardimids(5),iloc,jloc)-1)*sin4d

               iout = nd1 + (nd2-1)*sout1d + (nd3-1)*sout2d             &
                          + (nd4-1)*sout3d + (nd5-1)*sout4d

!               IF (debug > 2) WRITE(6,'(13x,a,2I2,3(a,I4),4I4)')        &
!                  'Processor - ',iloc,jloc,': Extracting from ',iin,    &
!                  ' to ',iout,' at ',nd1,nd2,nd3,nd4,nd5
!
               varaout(iout) = varain(iin)
             END DO
             END DO
             END DO
             END DO
             END DO

             istatus = nf_put_vara_real(foutids(iloc,jloc),varid,startidx,outidx,varaout)
             IF (istatus /= NF_NOERR) CALL handle_err(istatus)
          END DO
        END DO
      CASE DEFAULT
        WRITE(6,'(1x,a,I2)') 'ERROR: unsupported variable type = ',vartype
        istatus = -4
        RETURN
      END SELECT

    END DO
    !
    ! Close files
    !

    IF (debug > 0) WRITE(6,'(1x,a)') 'Closing all files ...'

    DO jloc = 1,nproc_y                   ! Close patches
      DO iloc = 1,nproc_x
        istatus = nf_close(foutids(iloc,jloc))
        IF (istatus /= NF_NOERR) CALL handle_err(istatus)
      END DO
    END DO

    istatus = nf_close(finid)                              ! Close file
    IF (istatus /= NF_NOERR) CALL handle_err(istatus)

  END DO

  DEALLOCATE(outfilenames, foutids, STAT = istatus)
  IF (ALLOCATED(varain)) DEALLOCATE(varain, varaout, STAT = istatus)
  IF (ALLOCATED(variin)) DEALLOCATE(variin, variout, STAT = istatus)

  RETURN
END SUBROUTINE splitncdf

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
