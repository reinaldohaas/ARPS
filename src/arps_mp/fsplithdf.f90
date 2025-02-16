!
SUBROUTINE  splithdf(filenames,nfile,dimnamein, xdimname, ydimname,     &
                     varidx_dim, nxidx, nyidx, nproc_x,nproc_y,         &
                     outdirname, debug,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Split files in HDF4 format into patches. The patched files will
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

  LOGICAL, INTENT(IN)            :: dimnamein
                   ! .TRUE.  get attributes xdimname/ydimname for global x size and y size
                   ! .FALSE. get x/y size by reading varidx with dimensions nxidx, nyidx
  INTEGER, INTENT(IN)            :: varidx_dim, nxidx, nyidx
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

  CHARACTER(LEN=1), ALLOCATABLE :: attvalstr(:)
  INTEGER,          ALLOCATABLE :: attvali(:)
  REAL,             ALLOCATABLE :: attvalr(:)

  INTEGER, ALLOCATABLE :: variin(:), variout(:)
  INTEGER, ALLOCATABLE :: varain(:), varaout(:)
  INTEGER(KIND=SELECTED_INT_KIND(4)), ALLOCATABLE :: vari16in(:), vari16out(:)

  INTEGER                         :: finid
  CHARACTER(LEN=256), ALLOCATABLE :: outfilenames(:,:)
  INTEGER,            ALLOCATABLE :: foutids(:,:)
  INTEGER,            ALLOCATABLE :: varido (:,:)
  INTEGER,            ALLOCATABLE :: dimouts0(:,:,:)
!
!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_RANK = 5
  INTEGER :: nf, dindx
  INTEGER :: iout, iin

  INTEGER :: nxid, nyid
  INTEGER :: iloc, jloc

  INTEGER :: comp_code, comp_prm(1)

  INTEGER :: ngatts, nvars               ! file info

  CHARACTER(LEN=256) :: dimname, attname, varname

  INTEGER :: attnum, atttype, attnval    ! attribute info

  INTEGER :: varid, varidx               ! variable info
  INTEGER :: vartype, varndims, varnatts

  INTEGER :: vardim, varinsize, varoutsize
  INTEGER :: variallsize, varaallsize, vari16allsize
  INTEGER :: varibllsize, varabllsize, vari16bllsize
  INTEGER :: vardimsize(MAX_RANK), vardimoutsize(MAX_RANK)
  INTEGER :: startidx(MAX_RANK), stride(MAX_RANK)

  INTEGER :: dimid                       ! dimension info

  INTEGER ::      sin1d,  sin2d,  sin3d,  sin4d
  INTEGER ::      sout1d, sout2d, sout3d, sout4d
  INTEGER :: nd1, nd2,    nd3,    nd4,    nd5      ! Assume the max rank is 5

  INTEGER :: allocstrlen, allocrlen, allocilen
  INTEGER :: resetx, resety
!
!-----------------------------------------------------------------------
!
! Including files
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'

!-----------------------------------------------------------------------
!
!  HDF Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfstart, sfend, sffinfo, sfginfo, sfgainfo, sfgcompress    ! file and quiry
  INTEGER :: sffattr, sfrnatt, sfsnatt, sfrcatt, sfscatt                ! attribute
  INTEGER :: sfdimid, sfsdmname                                         ! Dimension
  INTEGER :: sfselect, sfendacc, sfcreate, sfrdata, sfwdata             ! Variable
  INTEGER :: sfscompress

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!
! Check dimensions first
!
!-----------------------------------------------------------------------

  finid = sfstart(filenames(1),DFACC_READ)
  IF (finid <= 0) CALL print_err(finid,'Cannot open data set.')

  IF (dimnamein) THEN         ! Read global attributes for dimension size
    nxid    = sffattr(finid, xdimname)
    istatus = sfrnatt(finid, nxid, nxlg)

    nyid    = sffattr(finid, ydimname)
    istatus = sfrnatt(finid, nyid, nylg)
  ELSE                      ! Read variable dimensions
    varid   = sfselect(finid,varidx_dim)
    istatus = sfginfo (varid,varname,varndims,vardimsize,vartype,varnatts)
    nxlg = vardimsize(nxidx)
    nylg = vardimsize(nyidx)
  END IF

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

  istatus = sfend(finid)

!-----------------------------------------------------------------------
!
!  Loop over filenames
!
!-----------------------------------------------------------------------


  ALLOCATE(outfilenames(nproc_x,nproc_y), STAT = istatus)
  ALLOCATE(foutids(nproc_x,nproc_y),      STAT = istatus)
  ALLOCATE(varido (nproc_x,nproc_y),      STAT = istatus)
  ALLOCATE(dimouts0(MAX_RANK,nproc_x,nproc_y),   STAT = istatus)

  ALLOCATE(attvalr(1),                    STAT = istatus)
  ALLOCATE(attvali(1),                    STAT = istatus)
  ALLOCATE(attvalstr(1024),               STAT = istatus)
  allocrlen   = 1
  allocilen   = 1
  allocstrlen = 1024

  startidx(:) = 0
  stride(:)   = 1

  variallsize = 0
  vari16allsize = 0
  varaallsize = 0

  varibllsize = 0
  vari16bllsize = 0
  varabllsize = 0

  DO nf = 1,nfile

    WRITE(6,'(/,1x,a,I6)') '== Processing file No. ',nf

    IF (debug > 0) WRITE(6,'(/,1x,2a)') 'Opening file  - ',TRIM(filenames(nf))

    finid = sfstart(filenames(nf),DFACC_READ)
    IF (finid <= 0) CALL print_err(finid,'Cannot open data set.')

    dindx = INDEX(filenames(nf),'/',.TRUE.) + 1
    DO jloc = 1,nproc_y                   ! Create patches
      DO iloc = 1,nproc_x
        CALL gtsplitfn(TRIM(outdirname)//TRIM(filenames(nf)(dindx:)),   &
                       nproc_x,nproc_y,1,1,iloc,jloc,                   &
                       0,0,0,2,outfilenames(iloc,jloc),istatus)
        IF (debug > 0) WRITE(6,'(1x,2a)') 'Creating file - ',TRIM(outfilenames(iloc,jloc))
        foutids(iloc,jloc) = sfstart(outfilenames(iloc,jloc),DFACC_CREATE)
        IF (foutids(iloc,jloc) <= 0) CALL print_err(finid,'Cannot open data set.')
      END DO
    END DO

    !
    ! Set dimensions
    !

    nx = (nxlg-3)/nproc_x + 3
    ny = (nylg-3)/nproc_y + 3

    istatus = sffinfo(finid,nvars,ngatts)

!-----------------------------------------------------------------------
!
!  Read/write Global attributes
!
!-----------------------------------------------------------------------

    IF (debug > 1) WRITE(6,'(5x,a,I2)') 'Total global attributs - ',ngatts

    DO attnum = 0,ngatts-1

      istatus = sfgainfo(finid,attnum,attname,atttype,attnval)

      IF (attnval /= 1 .AND. atttype /= DFNT_CHAR8) THEN
        WRITE(6,'(1x,a,/)') 'ERROR: No support for attribute arrays still.'
        istatus = -3
        RETURN
      ELSE
        IF (debug > 2) WRITE(6,'(9x,2a)') 'Attribute name - ',TRIM(attname)
      END IF

      SELECT CASE (atttype)
      CASE (DFNT_CHAR8)         ! Character string
        IF (attnval > allocstrlen) THEN
          DEALLOCATE(attvalstr,          STAT = istatus)
          ALLOCATE  (attvalstr(attnval), STAT = istatus)
          allocstrlen = attnval
          attvalstr = ' '
        END IF
        istatus = sfrcatt(finid, attnum, attvalstr)

        DO jloc = 1,nproc_y                   ! Write CHARACTER global attributes
          DO iloc = 1,nproc_x
            istatus = sfscatt(foutids(iloc,jloc), TRIM(attname), atttype, attnval, attvalstr)
          END DO
        END DO

      CASE (DFNT_FLOAT32)       ! REAL number
        IF (attnval > allocrlen) THEN
          DEALLOCATE(attvalr,          STAT = istatus)
          ALLOCATE  (attvalr(attnval), STAT = istatus)
          allocrlen = attnval
          attvalr   = 0.0
        END IF
        istatus = sfrnatt(finid, attnum, attvalr)

        DO jloc = 1,nproc_y                   ! Write REAL global attributes
          DO iloc = 1,nproc_x
            istatus = sfsnatt(foutids(iloc,jloc), TRIM(attname), atttype, attnval, attvalr)
          END DO
        END DO

      CASE (DFNT_INT32)         ! INTEGER number
        IF (attnval > allocilen) THEN
          DEALLOCATE(attvali,          STAT = istatus)
          ALLOCATE  (attvali(attnval), STAT = istatus)
          allocilen = attnval
          attvali   = 0
        END IF
        IF ( TRIM(attname) == TRIM(xdimname) ) THEN
          attvali(1) = nx
        ELSE IF ( TRIM(attname) == TRIM(ydimname) ) THEN
          attvali(1) = ny
        ELSE
          istatus = sfrnatt(finid, attnum, attvali)
        END IF

        DO jloc = 1,nproc_y                   ! Write integer global attributes
          DO iloc = 1,nproc_x
            istatus = sfsnatt(foutids(iloc,jloc), TRIM(attname), atttype, attnval, attvali)
          END DO
        END DO

      CASE DEFAULT
        WRITE (6,'(/,1x,a,I2,3a,/)') 'ERROR: unsupported attribute type (',  &
                        atttype,') for attribute ', TRIM(attname),'.'
        istatus = -4
        RETURN
      END SELECT

    END DO

!-----------------------------------------------------------------------
!
!  Read/write each variable
!
!-----------------------------------------------------------------------

    IF (debug > 1) WRITE(6,'(5x,a,I2)') 'Total variables - ',nvars

    DO varidx = 0,nvars-1
      vardimsize(:) = 1

      varid = sfselect(finid,varidx)
      istatus = sfginfo (varid,varname,varndims,vardimsize,vartype,varnatts)

      IF (debug >2) WRITE(6,'(9x,2a)') 'Processing variables - ',TRIM(varname)

      istatus = sfgcompress(varid,comp_code,comp_prm)
      IF (istatus /= SUCCEED) THEN
!        WRITE(6,'(1x,3a)') 'ERROR: HDF error from sfgcompress with variable ',TRIM(varname),'.'
!        RETURN
         comp_code   = 0
         comp_prm(1) = 0
      END IF

      IF (debug >3) WRITE(6,'(13x,2(a,I2))') 'comp_code = ',comp_code,', comp_prm = ',comp_prm(1)

!      attnum = sffattr(varid, 'hdf_comp_prm')     ! Whether the data set is compressed?
!      IF (attnum == FAIL) THEN
!        comp_prm(1) = 0
!      ELSE
!        istatus = sfrnatt(varid, attnum, comp_prm)
!      END IF
!      comp_code = 0
!      IF (comp_prm(1) /= 0) THEN
!        attnum  = sffattr(varid, 'hdf_comp_code')
!        istatus = sfrnatt(varid, attnum, comp_code)
!      END IF

      vardimoutsize(:) = vardimsize(:)

      resetx = 0             ! Assume we only split the first two dimensions
      resety = 0
      IF (vardimsize(1) == nxlg) THEN         ! Check the first dimension of nx/ny
        vardimoutsize(1) = nx
        resetx = 1
      ELSE IF (vardimsize(1) == nylg) THEN
        vardimoutsize(1) = ny
        resety = 1
      END IF
                                              ! Check the second dimension as necessary
      IF (varndims > 1 .AND. vardimsize(2) == nylg) THEN
        vardimoutsize(2) = ny
        resety = 2
      END IF

      varinsize  = 1
      varoutsize = 1
      DO vardim = 1, varndims
        varinsize   = vardimsize(vardim)*varinsize
        varoutsize  = vardimoutsize(vardim)*varoutsize

        IF (debug > 4) WRITE(6,'(13x,2(I4,a))')   &
              vardimsize(vardim),'(in) - ', vardimoutsize(vardim),  '(out).'
      END DO

      sin1d = vardimsize(1)          ! size of one column,  INPUT
      sin2d = sin1d*vardimsize(2)    ! size of one slice (xy)
      sin3d = sin2d*vardimsize(3)    ! size of one cell  (xyz)
      sin4d = sin3d*vardimsize(4)

      sout1d = vardimoutsize(1)           ! size of one column, OUTPUT
      sout2d = sout1d*vardimoutsize(2)
      sout3d = sout2d*vardimoutsize(3)
      sout4d = sout3d*vardimoutsize(4)

      dimouts0(:,:,:) = 0         ! The starting index of each dimensions for each subdomain

      DO jloc = 1,nproc_y                   ! Create data set
        DO iloc = 1,nproc_x
          varido(iloc,jloc) = sfcreate(foutids(iloc,jloc), varname,     &
                                       vartype, varndims, vardimoutsize)
          IF (resetx > 0) THEN
            dimid = sfdimid(varido(iloc,jloc),resetx-1)
            istatus = sfsdmname(dimid,xdimname)
            dimouts0(resetx,iloc,jloc) = (iloc-1)*(nx-3)
          END IF
          IF (resety > 0) THEN
            dimid = sfdimid(varido(iloc,jloc),resety-1)
            istatus = sfsdmname(dimid,ydimname)
            dimouts0(resety,iloc,jloc) = (jloc-1)*(ny-3)
          END IF
          istatus = sfscompress(varido(iloc,jloc), comp_code, comp_prm)
        END DO
      END DO

      SELECT CASE (vartype)
      CASE (dfnt_float32)
        IF (varinsize > varaallsize) THEN   ! Allocate input array only when necessary
          IF (ALLOCATED(varain)) DEALLOCATE(varain, STAT = istatus)
          ALLOCATE(varain(varinsize),               STAT = istatus)
          varaallsize = varinsize
        END IF

        IF (varoutsize > varabllsize) THEN  ! Allocate output array only when necessary
          IF (ALLOCATED(varaout)) DEALLOCATE(varaout, STAT = istatus)
          ALLOCATE(varaout(varoutsize), STAT = istatus)
          varabllsize = varoutsize
        END IF

        istatus = sfrdata(varid, startidx, stride, vardimsize, varain)

        DO jloc = 1,nproc_y                  ! Write patches variables
          DO iloc = 1,nproc_x

             IF (debug > 3) THEN
               WRITE(6,FMT='(13x,a,2(I4,a))',ADVANCE='NO')                  &
                           'Writing to processor - (',iloc,',',jloc,') with: '
               DO vardim = 1,varndims
                 WRITE(6,FMT='(I4,a)',ADVANCE='NO') dimouts0(vardim,iloc,jloc),' '
               END DO
               WRITE(6,*)
             END IF

             DO nd5 = 1, vardimoutsize(5)        ! Assume max rank is 5, IMPORTANT
             DO nd4 = 1, vardimoutsize(4)
             DO nd3 = 1, vardimoutsize(3)
             DO nd2 = 1, vardimoutsize(2)
             DO nd1 = 1, vardimoutsize(1)
               iin =  nd1+dimouts0(1,iloc,jloc)              &
                   + (nd2+dimouts0(2,iloc,jloc)-1)*sin1d     &
                   + (nd3+dimouts0(3,iloc,jloc)-1)*sin2d     &
                   + (nd4+dimouts0(4,iloc,jloc)-1)*sin3d     &
                   + (nd5+dimouts0(5,iloc,jloc)-1)*sin4d

               iout = nd1 + (nd2-1)*sout1d + (nd3-1)*sout2d             &
                          + (nd4-1)*sout3d + (nd5-1)*sout4d

               IF (debug > 3) WRITE(6,'(13x,a,2I2,3(a,I4),4I4)')        &
                  'Processor - ',iloc,jloc,': Extracting from ',iin,    &
                  ' to ',iout,' at ',nd1,nd2,nd3,nd4,nd5

               varaout(iout) = varain(iin)
             END DO
             END DO
             END DO
             END DO
             END DO

             istatus = sfwdata(varido(iloc,jloc), startidx, stride, vardimoutsize, varaout)
             IF (istatus /= SUCCEED) CALL print_err(istatus,'ERROR: sfwdata')
          END DO
        END DO

      CASE (dfnt_int32)

        IF (varinsize > variallsize) THEN   ! Allocate input array only when necessary
          IF (ALLOCATED(variin)) DEALLOCATE(variin, STAT = istatus)
          ALLOCATE(variin(varinsize),               STAT = istatus)
          variallsize = varinsize
        END IF

        IF (varoutsize > varibllsize) THEN  ! Allocate output array only when necessary
          IF (ALLOCATED(variout)) DEALLOCATE(variout, STAT = istatus)
          ALLOCATE(variout(varoutsize), STAT = istatus)
          varibllsize = varoutsize
        END IF

        istatus = sfrdata(varid, startidx, stride, vardimsize, variin)

        DO jloc = 1,nproc_y                  ! Write patches variables
          DO iloc = 1,nproc_x

             IF (debug > 3) THEN
               WRITE(6,FMT='(13x,a,2(I4,a))',ADVANCE='NO')                  &
                           'Writing to processor - (',iloc,',',jloc,') with: '
               DO vardim = 1,varndims
                 WRITE(6,FMT='(I4,a)',ADVANCE='NO') dimouts0(vardim,iloc,jloc),' '
               END DO
               WRITE(6,*)
             END IF

             DO nd5 = 1, vardimoutsize(5)        ! Assume max rank is 5, IMPORTANT
             DO nd4 = 1, vardimoutsize(4)
             DO nd3 = 1, vardimoutsize(3)
             DO nd2 = 1, vardimoutsize(2)
             DO nd1 = 1, vardimoutsize(1)
               iin =  nd1+dimouts0(1,iloc,jloc)              &
                   + (nd2+dimouts0(2,iloc,jloc)-1)*sin1d     &
                   + (nd3+dimouts0(3,iloc,jloc)-1)*sin2d     &
                   + (nd4+dimouts0(4,iloc,jloc)-1)*sin3d     &
                   + (nd5+dimouts0(5,iloc,jloc)-1)*sin4d

               iout = nd1 + (nd2-1)*sout1d + (nd3-1)*sout2d             &
                          + (nd4-1)*sout3d + (nd5-1)*sout4d

               IF (debug > 3) WRITE(6,'(13x,a,2I2,3(a,I4),4I4)')        &
                  'Processor - ',iloc,jloc,': Extracting from ',iin,    &
                  ' to ',iout,' at ',nd1,nd2,nd3,nd4,nd5

               variout(iout) = variin(iin)
             END DO
             END DO
             END DO
             END DO
             END DO

             istatus = sfwdata(varido(iloc,jloc), startidx, stride, vardimoutsize, variout)
             IF (istatus /= SUCCEED) CALL print_err(istatus,'ERROR: sfwdata')
          END DO
        END DO

      CASE (dfnt_int16)

        IF (varinsize > vari16allsize) THEN   ! Allocate input array only when necessary
          IF (ALLOCATED(vari16in)) DEALLOCATE(vari16in, STAT = istatus)
          ALLOCATE(vari16in(varinsize),               STAT = istatus)
          vari16allsize = varinsize
        END IF

        IF (varoutsize > vari16bllsize) THEN  ! Allocate output array only when necessary
          IF (ALLOCATED(vari16out)) DEALLOCATE(vari16out, STAT = istatus)
          ALLOCATE(vari16out(varoutsize), STAT = istatus)
          vari16bllsize = varoutsize
        END IF

        istatus = sfrdata(varid, startidx, stride, vardimsize, vari16in)

        DO jloc = 1,nproc_y                  ! Write patches variables
          DO iloc = 1,nproc_x

             IF (debug > 3) THEN
               WRITE(6,FMT='(13x,a,2(I4,a))',ADVANCE='NO')                  &
                           'Writing to processor - (',iloc,',',jloc,') with: '
               DO vardim = 1,varndims
                 WRITE(6,FMT='(I4,a)',ADVANCE='NO') dimouts0(vardim,iloc,jloc),' '
               END DO
               WRITE(6,*)
             END IF

             DO nd5 = 1, vardimoutsize(5)        ! Assume max rank is 5, IMPORTANT
             DO nd4 = 1, vardimoutsize(4)
             DO nd3 = 1, vardimoutsize(3)
             DO nd2 = 1, vardimoutsize(2)
             DO nd1 = 1, vardimoutsize(1)
               iin =  nd1+dimouts0(1,iloc,jloc)              &
                   + (nd2+dimouts0(2,iloc,jloc)-1)*sin1d     &
                   + (nd3+dimouts0(3,iloc,jloc)-1)*sin2d     &
                   + (nd4+dimouts0(4,iloc,jloc)-1)*sin3d     &
                   + (nd5+dimouts0(5,iloc,jloc)-1)*sin4d

               iout = nd1 + (nd2-1)*sout1d + (nd3-1)*sout2d             &
                          + (nd4-1)*sout3d + (nd5-1)*sout4d

               IF (debug > 3) WRITE(6,'(13x,a,2I2,3(a,I4),4I4)')        &
                  'Processor - ',iloc,jloc,': Extracting from ',iin,    &
                  ' to ',iout,' at ',nd1,nd2,nd3,nd4,nd5

               vari16out(iout) = vari16in(iin)
             END DO
             END DO
             END DO
             END DO
             END DO

             istatus = sfwdata(varido(iloc,jloc), startidx, stride, vardimoutsize, vari16out)
             IF (istatus /= SUCCEED) CALL print_err(istatus,'ERROR: sfwdata')
          END DO
        END DO

      CASE DEFAULT
        WRITE(6,'(1x,a,I2,3a,/)') 'ERROR: Unsupported variable type (', &
                            vartype,') for variable ',TRIM(varname),'.'
        istatus = -5
        RETURN
      END SELECT

      !
      ! Data set attributes
      !

      DO attnum = 0,varnatts-1

        istatus = sfgainfo(varid,attnum,attname,atttype,attnval)

        IF (debug > 2) WRITE(6,'(13x,2a)') 'Attribute name - ',TRIM(attname)

        SELECT CASE (atttype)
        CASE (DFNT_CHAR8)         ! Character string
          IF (attnval > allocstrlen) THEN
            DEALLOCATE(attvalstr,          STAT = istatus)
            ALLOCATE  (attvalstr(attnval), STAT = istatus)
            allocstrlen = attnval
            attvalstr = ' '
          END IF
          istatus = sfrcatt(varid, attnum, attvalstr)

          DO jloc = 1,nproc_y                   ! Write CHARACTER var attributes
            DO iloc = 1,nproc_x
              istatus = sfscatt(varido(iloc,jloc), TRIM(attname), atttype, attnval, attvalstr)
            END DO
          END DO

        CASE (DFNT_FLOAT32)       ! REAL number
          IF (attnval > allocrlen) THEN
            DEALLOCATE(attvalr,          STAT = istatus)
            ALLOCATE  (attvalr(attnval), STAT = istatus)
            allocrlen = attnval
            attvalr   = 0.0
          END IF
          istatus = sfrnatt(varid, attnum, attvalr)

          DO jloc = 1,nproc_y                   ! Write REAL global attributes
            DO iloc = 1,nproc_x
              istatus = sfsnatt(varido(iloc,jloc), TRIM(attname), atttype, attnval, attvalr)
            END DO
          END DO

        CASE (DFNT_INT32)         ! INTEGER number
          IF (attnval > allocilen) THEN
            DEALLOCATE(attvali,          STAT = istatus)
            ALLOCATE  (attvali(attnval), STAT = istatus)
            allocilen = attnval
            attvali   = 0
          END IF
          istatus = sfrnatt(varid, attnum, attvali)

          DO jloc = 1,nproc_y                   ! Write integer global attributes
            DO iloc = 1,nproc_x
              istatus = sfsnatt(varido(iloc,jloc), TRIM(attname), atttype, attnval, attvali)
            END DO
          END DO

        CASE DEFAULT
          WRITE (6,'(/,1x,a,I2,5a,/)') 'ERROR: unsupported attribute type (',  &
                                   atttype,') for attribute ', TRIM(attname),  &
                                   ' of variable ',TRIM(varname),'.'
          istatus = -4
          RETURN
        END SELECT
      END DO

      !
      ! Close data sets
      !
      istatus = sfendacc(varid)
      DO jloc = 1,nproc_y
        DO iloc = 1,nproc_x
          istatus = sfendacc(varido(iloc,jloc))
        END DO
      END DO

    END DO

    !
    ! Close files
    !

    IF (debug > 0) WRITE(6,'(1x,a)') 'Closing all files ...'

    DO jloc = 1,nproc_y                   ! Close patches
      DO iloc = 1,nproc_x
        istatus = sfend(foutids(iloc,jloc))
      END DO
    END DO

    istatus = sfend(finid)                              ! Close file

  END DO

  DEALLOCATE(outfilenames, foutids, varido, STAT = istatus)
  DEALLOCATE(attvalstr, attvali, attvalr,   STAT = istatus)
  DEALLOCATE(dimouts0,                      STAT = istatus)
  IF (ALLOCATED(varain))   DEALLOCATE(varain,   varaout,   STAT = istatus)
  IF (ALLOCATED(variin))   DEALLOCATE(variin,   variout,   STAT = istatus)
  IF (ALLOCATED(vari16in)) DEALLOCATE(vari16in, vari16out, STAT = istatus)

  RETURN
END SUBROUTINE splithdf

SUBROUTINE print_err(istat,message)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: istat
  CHARACTER(*), INTENT(IN) :: message

  PRINT *, TRIM(message)
  STOP 'HDF error!'

  RETURN
END SUBROUTINE print_err
