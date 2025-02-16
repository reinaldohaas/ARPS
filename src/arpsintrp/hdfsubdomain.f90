PROGRAM hdfsubdomain
!
!##################################################################
!##################################################################
!######                                                      ######
!######                PROGRAM HDFSUBDOMAIN                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Program to read history data dump from ARPS HDF 4 format and
!  output a subset of it in a smaller domain.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    02/16/2009
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  USE module_hdfsubdomain_namelist
  IMPLICIT NONE

  INTEGER, PARAMETER :: STDIN = 5, STDOUT = 6
  LOGICAL :: IAMROOT
  INTEGER :: UNTIN, UNTOUT

  INTEGER :: iargc, nargs

  CHARACTER(LEN=256) :: NML_FILE

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istatus
  INTEGER :: nf, nvar, ioutvar

  INTEGER, ALLOCATABLE :: fHandles(:,:), fHndlOut(:,:)

  INTEGER :: nprocx_in, nprocy_in

  INTEGER :: nx, ny, nz, nzsoil, nstyps, nxlg, nylg
  INTEGER :: nxout, nyout, nxlgout, nylgout
  INTEGER :: nxsm, nysm, nksize, nssize
  INTEGER :: iskip, jskip
  INTEGER :: numvar

  LOGICAL :: first_time

  CHARACTER(LEN=256) :: vname
  INTEGER :: rank, dims(6), dtype, ndattr

  INTEGER, ALLOCATABLE :: ivardtain(:,:,:), ivardtaout(:,:,:)
  REAL,    ALLOCATABLE :: vardtain(:,:,:),  vardtaout(:,:,:)
  INTEGER, ALLOCATABLE :: ivarlgin(:,:,:), ivarlgout(:,:,:)
  REAL,    ALLOCATABLE :: varlgin(:,:,:),  varlgout(:,:,:)

  REAL,    ALLOCATABLE :: vardtain4d(:,:,:,:),  vardtaout4d(:,:,:,:)
  REAL,    ALLOCATABLE :: varlgin4d(:,:,:,:),  varlgout4d(:,:,:,:)

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'
  INCLUDE 'hdf.f90'

  INTEGER :: check_var_name

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = 0
  UNTIN   = STDIN
  UNTOUT  = STDOUT

  CALL mpinit_proc(0)

  IAMROOT = .FALSE.
  IF (myproc == 0) IAMROOT = .TRUE.

  IF (IAMROOT) WRITE(6,'(/9(/5x,a)/)')                                  &
     '###############################################################', &
     '###############################################################', &
     '####                                                       ####', &
     '####               Welcome to HDFSUBDOMAIN                 ####', &
     '####     This program subtract a subset of the orignal     ####', &
     '####     datasets in the ARPS HDF 4 format.                ####', &
     '####                                                       ####', &
     '###############################################################', &
     '###############################################################'

!-----------------------------------------------------------------------
!
!  Read namelist.
!
!-----------------------------------------------------------------------
!
  IF ( IAMROOT ) THEN
    nargs = iargc()
    IF (nargs == 1) THEN
      CALL getarg(nargs, NML_FILE)
      CALL getunit( UNTIN )
      WRITE(UNTOUT,'(1x,3a)') 'Reading namelist parameters from file "',&
                              TRIM(NML_FILE),'"'
      OPEN(UNIT=UNTIN,FILE=TRIM(NML_FILE),STATUS='OLD',ACTION='READ',   &
           FORM='FORMATTED',IOSTAT=istatus)
    ELSE
      WRITE(UNTOUT,'(1x,a)') 'Reading namelist parameters from STDIN.'
    END IF
  END IF

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) THEN
    IF (IAMROOT) WRITE(UNTOUT,'(1x,3a,I3)') 'IO ERROR with file ',      &
                 TRIM(NML_FILE),', with error code = ',istatus
    CALL arpsstop('IO ERROR with namelist file.',1)
  END IF

  CALL read_namelist_params(iamroot,untin,untout,istatus)

  IF ( IAMROOT .AND. UNTIN /= STDIN ) THEN
    CLOSE(UNTIN)
    CALL retunit( UNTIN )
  END IF

!-----------------------------------------------------------------------
!
! Get dimensions and allocate working arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(fHandles(ncompressx,ncompressy), STAT = istatus)
  ALLOCATE(fHndlOut(nprocx_out,nprocy_out), STAT = istatus)

!-----------------------------------------------------------------------
!
! Process each files
!
!-----------------------------------------------------------------------

  nprocx_in = nproc_x*ncompressx
  nprocy_in = nproc_y*ncompressy

  first_time = .TRUE.

  DO nf = 1, nhisfile

    IF (lvldbg > 0) WRITE(UNTOUT,'(1x,a,I4,2a,/,31x,2(a,2(I3,a)))')     &
        '=== Processing file No.',nf,' - ',TRIM(hisfile(nf)),           &
        'Number of patches (',nprocx_in,',',nprocy_in,') ',             &
        'starting from (',nprocx_lw,',',nprocy_lw,').'

    CALL open_hdf(hisfile(nf),ncompressx,ncompressy,loc_x,loc_y,        &
                  nprocx_lw,nprocy_lw,iamroot,inpatch,1,lvldbg,         &
                  fHandles,istatus)
    CALL mpmini(istatus)    ! The minimum value of all processes
    IF (istatus /= 0) CYCLE

    IF (IAMROOT) THEN
      CALL open_hdf(outfile(nf),nprocx_out,nprocy_out,loc_x,loc_y,      &
                    1,1,iamroot,outpatch,2,lvldbg,fHndlOut,istatus)
      IF (istatus /= 0) CALL arpsstop('Error inside open_hdf.',1)

      IF (lvldbg > 1) WRITE(UNTOUT,'(3x,a)') '--- Coping global attributes'

      CALL copy_global_attributes(IAMROOT,                              &
                        fHandles,ncompressx,ncompressy,                 &
                        nprocx_in,nprocy_in,nxpnt,nypnt,runname,        &
                        fHndlOut,nprocx_out,nprocy_out,                 &
                        nx,ny,nz,nzsoil,nstyps,numvar,istatus)
    END IF
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) CALL arpsstop('Error in copy_global_attributes.',1)

    CALL mpupdatei(nx,1)
    CALL mpupdatei(ny,1)
    CALL mpupdatei(nz,1)
    CALL mpupdatei(nzsoil,1)
    CALL mpupdatei(nstyps,1)
    CALL mpupdatei(numvar,1)

    IF (first_time) THEN
      nxlg    = (nx-3)*nproc_x+3
      nylg    = (ny-3)*nproc_y+3
      nxlgout = nxlg - nxpnt
      nylgout = nylg - nypnt
      nxout   = (nxlgout-3)/nprocx_out + 3
      nyout   = (nylgout-3)/nprocy_out + 3
      iskip   = nxpnt / 2
      jskip   = nypnt / 2

      nxsm = (nx-3)/ncompressx + 3
      nysm = (ny-3)/ncompressy + 3

      ALLOCATE(vardtain  (nx,     ny,     nz),     STAT = istatus)
      ALLOCATE(vardtaout (nxout,  nyout,  nz),     STAT = istatus)

      ALLOCATE(ivardtain (nx,     ny,     nstyps), STAT = istatus)
      ALLOCATE(ivardtaout(nxout,  nyout,  nstyps), STAT = istatus)

      ALLOCATE(vardtain4d  (nx,    ny,    nzsoil, 0:nstyps), STAT = istatus)
      ALLOCATE(vardtaout4d (nxout, nyout, nzsoil, 0:nstyps), STAT = istatus)

      IF (IAMROOT) THEN
        ALLOCATE(varlgin  (nxlg,   nylg,   nz),     STAT = istatus)
        ALLOCATE(varlgout (nxlgout,nylgout,nz),     STAT = istatus)

        ALLOCATE(ivarlgin (nxlg,   nylg,   nstyps), STAT = istatus)
        ALLOCATE(ivarlgout(nxlgout,nylgout,nstyps), STAT = istatus)

        ALLOCATE(varlgin4d  (nxlg,   nylg,   nzsoil, 0:nstyps),     STAT = istatus)
        ALLOCATE(varlgout4d (nxlgout,nylgout,nzsoil, 0:nstyps),     STAT = istatus)
      END IF

      first_time = .FALSE.
    END IF

    DO nvar = 0, numvar-1

      dims(:) = 1
      CALL peek_hdf_dataset(fHandles(1,1),nvar,vname,rank,dims,dtype,ndattr,istatus)

      ioutvar = check_var_name(vname,varnames,nvarout)

      IF (ioutvar > 0) THEN

        SELECT CASE ( dtype )
        CASE ( dfnt_int32 )  ! soiltyp, vegtyp etc.

          IF (lvldbg > 1) WRITE(UNTOUT,'(3x,2a)') '--- Writing integer variable ',TRIM(vname)

          IF (dims(1) == nxsm .AND. dims(2) == nysm) THEN
            nksize = dims(3)

            CALL read_merge_write_int(IAMROOT,fHandles,ncompressx,ncompressy, &
                         vname,rank,dims,dtype,ivardtain,ivarlgin,            &
                         nx,ny,nksize,nxlg,nylg,iskip,jskip,                  &
                         fHndlOut,nprocx_out,nprocy_out,ivardtaout,ivarlgout, &
                         nxout,nyout,nxlgout,nylgout,istatus)

          ELSE
            WRITE(UNTOUT,'(1x,2a,/,8x,a)')                            &
               'ERROR: Found an unsupported integer array in ARPS history file - ',TRIM(vname), &
               'Program stopped.'
            CALL arpsstop('Unsupported integer variable found.',1)
          END IF

        CASE ( dfnt_float32, dfnt_int16 )  ! all other variable, packed or unpacked

          SELECT CASE (rank)
          CASE (1)
            IF (lvldbg > 1) WRITE(UNTOUT,'(3x,2a)') '--- Writing 1D variable ',TRIM(vname)

            IF (TRIM(vname) == 'x') THEN    ! 'X'
              CALL read_merge_write_real1d(IAMROOT,fHandles,ncompressx,ncompressy, &
                           vname,dims,'x',vardtain,varlgin, nx,nxlg, iskip,   &
                           fHndlOut,nprocx_out,nprocy_out,vardtaout,varlgout, &
                           nxout,nxlgout,lvldbg,istatus)

            ELSE IF (TRIM(vname) == 'y') THEN
              CALL read_merge_write_real1d(IAMROOT,fHandles,ncompressx,ncompressy, &
                           vname,dims,'y',vardtain,varlgin, ny,nylg, jskip,   &
                           fHndlOut,nprocx_out,nprocy_out,vardtaout,varlgout, &
                           nyout,nylgout,lvldbg,istatus)
            ELSE IF (TRIM(vname) == 'z') THEN
              CALL copy_hdf_real1d(IAMROOT,fHandles,ncompressx,ncompressy,    &
                           vname,dims,vardtain,nz,                            &
                           fHndlOut,nprocx_out,nprocy_out,vardtaout,istatus)
            ELSE
              WRITE(UNTOUT,'(1x,2a,/,8x,a)')                                  &
                 'ERROR: Found an unsupported 1D array in ARPS history file - ',TRIM(vname), &
                 'Program stopped.'
              CALL arpsstop('Unsupported variable found.',1)
            END IF

          CASE (2, 3)

            IF (lvldbg > 1) WRITE(UNTOUT,'(3x,2a)') '--- Writing 2D/3D variable ',TRIM(vname)

            IF (dims(1) == nxsm .AND. dims(2) == nysm) THEN
              nksize = dims(3)

              CALL read_merge_write_real(IAMROOT,fHandles,ncompressx,ncompressy, &
                           vname,rank,dims,dtype,vardtain,varlgin,            &
                           nx,ny,nksize,nxlg,nylg, iskip,jskip,               &
                           fHndlOut,nprocx_out,nprocy_out,vardtaout,varlgout, &
                           nxout,nyout,nxlgout,nylgout,nz,nzsoil,nstyps,      &
                           lvldbg,istatus)

            ELSE
              WRITE(UNTOUT,'(1x,2a,/,8x,a)')                            &
                 'ERROR: Found an unsupported array in ARPS history file - ',TRIM(vname), &
                 'Program stopped.'
              CALL arpsstop('Unsupported variable found.',1)
            END IF
          CASE (4)

            IF (lvldbg > 1) WRITE(UNTOUT,'(3x,2a)') '--- Writing 4D variable ',TRIM(vname)

            IF (dims(1) == nxsm .AND. dims(2) == nysm) THEN
              nksize = dims(3)
              nssize = dims(4)

              CALL read_merge_write_real4d(IAMROOT,fHandles,ncompressx,ncompressy, &
                           vname,rank,dims,dtype,vardtain4d,varlgin4d,             &
                           nx,ny,nksize,nssize,nxlg,nylg, iskip,jskip,             &
                           fHndlOut,nprocx_out,nprocy_out,vardtaout4d,varlgout4d,  &
                           nxout,nyout,nxlgout,nylgout,istatus)

            ELSE
              WRITE(UNTOUT,'(1x,2a,/,8x,a)')                            &
                 'ERROR: Found an unsupported array in ARPS history file - ',TRIM(vname), &
                 'Program stopped.'
              CALL arpsstop('Unsupported variable found.',1)
            END IF

          CASE DEFAULT
            WRITE(UNTOUT,'(1x,a,I4,/,8x,a)')                            &
               'ERROR: Found strange rank - ',rank,                     &
               'Program stopped.'
            CALL arpsstop('Unsupported rank found.',1)
          END SELECT

        CASE DEFAULT
            WRITE(UNTOUT,'(1x,a,I4,/,8x,a)')                            &
               'ERROR: Found Unsupported data type - ',dtype,           &
               'Program stopped.'
            CALL arpsstop('Unsupported data type found.',1)

        END SELECT

      END IF

    END DO

    IF (lvldbg > 0) WRITE(UNTOUT,'(1x,3a,2(I2,a),/)')                       &
    '~~~ Closing file ',TRIM(outfile(nf)),' (',nprocx_out,',',nprocy_out,').'

    CALL close_hdf(fHandles,ncompressx,ncompressy,istatus)
    IF (IAMROOT) CALL close_hdf(fHndlOut,nprocx_out,nprocy_out,istatus)

  END DO

!-----------------------------------------------------------------------
!
! Ending of the program
!
!-----------------------------------------------------------------------

  DEALLOCATE(fHandles)
  DEALLOCATE(fHndlOut)

  DEALLOCATE( vardtain,   vardtaout )
  DEALLOCATE( ivardtain,  ivardtaout )
  DEALLOCATE( vardtain4d, vardtaout4d )

  IF (IAMROOT) THEN
    DEALLOCATE( varlgin,   varlgout )
    DEALLOCATE( ivarlgin,  ivarlgout )
    DEALLOCATE( varlgin4d, varlgout4d )
  END IF

  IF (IAMROOT) WRITE(UNTOUT,'(/,1x,a,/)') '=== Program terminated normally ==='

  CALL mpexit(0)
END PROGRAM hdfsubdomain

!########################################################################
FUNCTION check_var_name(vname,varnames,nvarout)

  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN) :: vname
  INTEGER,            INTENT(IN) :: nvarout
  CHARACTER(LEN=20),  INTENT(IN) :: varnames(nvarout)
!-----------------------------------------------------------------------

  INTEGER :: check_var_name

  INTEGER :: n

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (nvarout < 1) THEN
    check_var_name = 1
  ELSE
    check_var_name = 0
    DO n = 1, nvarout
      IF (TRIM(vname) == TRIM(varnames(n))) THEN
        check_var_name = n
        EXIT
      END IF
    END DO
  END IF

  RETURN
END FUNCTION

!########################################################################
SUBROUTINE read_merge_write_int(IAMROOT,fHandles,nxpatch,nypatch,       &
                         vname,rank,dims,dtype,ivarin,ivarlgin,         &
                         nx,ny,nz,nxlg,nylg,iskip,jskip,                &
                         fHndlOut,nprocx_out,nprocy_out,ivarout,ivarlgout, &
                         nxout,nyout,nxlgout,nylgout,istatus)
!-----------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: IAMROOT
  INTEGER, INTENT(IN) :: nxpatch, nypatch
  INTEGER, INTENT(IN) :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN) :: vname
  INTEGER, INTENT(IN) :: rank, dims(6), dtype
  INTEGER, INTENT(IN) :: nx,ny,nz, nxlg, nylg
  INTEGER, INTENT(IN) :: iskip, jskip

  INTEGER, INTENT(IN) :: nprocx_out, nprocy_out
  INTEGER, INTENT(IN) :: fHndlOut(nprocx_out,nprocy_out)
  INTEGER, INTENT(IN) :: nxout, nyout, nxlgout, nylgout

  INTEGER, INTENT(INOUT) :: ivarin(nx,ny,nz)
  INTEGER, INTENT(INOUT) :: ivarlgin(nxlg,nylg,nz)
  INTEGER, INTENT(INOUT) :: ivarout(nxout,nyout,nz)
  INTEGER, INTENT(INOUT) :: ivarlgout(nxlgout,nylgout,nz)

  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: comp_code, comp_prm, stag_dim
  INTEGER :: lcmnt, lunts
  CHARACTER(LEN=256) :: comment, units

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL read_hdf_dataseti(fHandles,nxpatch,nypatch,vname,                &
                         ivarin,dims,nx,ny,nz,istatus)

  CALL mpimerge3di(ivarin,nx,ny,nz,ivarlgin)

  IF (IAMROOT) THEN
    DO k = 1, nz
      DO j = 1,nylgout
        DO i = 1,nxlgout
          ivarlgout(i,j,k) = ivarlgin(i+iskip,j+jskip,k)
        END DO
      END DO
    END DO

    CALL read_hdf_data_attr(fHandles(1,1),vname,comp_code,comp_prm,     &
                            comment, lcmnt, units, lunts, stag_dim,     &
                            istatus)

    CALL write_hdf_dataseti(fHndlOut,nprocx_out,nprocy_out,vname,       &
                  comp_code,comp_prm,comment,lcmnt,units,lunts,stag_dim,&
                  ivarlgout,nxlgout,nylgout,                            &
                  ivarout,nxout,nyout,nz,istatus)
  END IF

  RETURN
END SUBROUTINE read_merge_write_int

!########################################################################
SUBROUTINE read_merge_write_real(IAMROOT,fHandles,nxpatch,nypatch,      &
                         vname,rank,dims,dtype,varin,varlgin,            &
                         nx,ny,nz,nxlg,nylg,iskip,jskip,                 &
                         fHndlOut,nprocx_out,nprocy_out,varout,varlgout, &
                         nxout,nyout,nxlgout,nylgout,                    &
                         nzin,nzsoilin,nstypsin,dbglvl,istatus)
!-----------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: IAMROOT
  INTEGER, INTENT(IN) :: nxpatch, nypatch
  INTEGER, INTENT(IN) :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN) :: vname
  INTEGER, INTENT(IN) :: rank, dims(6), dtype
  INTEGER, INTENT(IN) :: nx,ny,nz, nxlg, nylg
  INTEGER, INTENT(IN) :: iskip, jskip

  INTEGER, INTENT(IN) :: nprocx_out, nprocy_out
  INTEGER, INTENT(IN) :: fHndlOut(nprocx_out,nprocy_out)
  INTEGER, INTENT(IN) :: nxout, nyout, nxlgout, nylgout

  REAL,    INTENT(INOUT) :: varin(nx,ny,nz)
  REAL,    INTENT(INOUT) :: varlgin(nxlg,nylg,nz)
  REAL,    INTENT(INOUT) :: varout(nxout,nyout,nz)
  REAL,    INTENT(INOUT) :: varlgout(nxlgout,nylgout,nz)

  INTEGER, INTENT(IN)    :: nzin, nzsoilin, nstypsin
  INTEGER, INTENT(IN)    :: dbglvl
  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: comp_code, comp_prm, stag_dim
  INTEGER :: lcmnt, lunts
  CHARACTER(LEN=256) :: comment, units
  INTEGER :: packed16
  REAL    :: amin, amax

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL read_hdf_dataset(fHandles,nxpatch,nypatch,vname,varin,           &
                        dims,nx,ny,nz,packed16,istatus)

  CALL mpimerge3d(varin,nx,ny,nz,varlgin)

  IF (IAMROOT) THEN

    IF (dbglvl > 2) THEN
      amin = varlgin(1,1,1)
      amax = varlgin(1,1,1)
      DO k = 1, nz
        DO j = 1, nylg
          DO i = 1, nxlg
            amin = MIN(amin,varlgin(i,j,k))
            amax = MAX(amax,varlgin(i,j,k))
          END DO
        END DO
      END DO
      WRITE(6,'(7x,2(a,F15.5))') 'Read in MINVAL = ',amin,', MAXVAL = ',amax
    END IF

    DO k = 1, nz
      DO j = 1,nylgout
        DO i = 1,nxlgout
          varlgout(i,j,k) = varlgin(i+iskip,j+jskip,k)
        END DO
      END DO
    END DO

    IF (dbglvl > 2) THEN
      amin = varlgout(1,1,1)
      amax = varlgout(1,1,1)
      DO k = 1, nz
        DO j = 1, nylgout
          DO i = 1, nxlgout
            amin = MIN(amin,varlgout(i,j,k))
            amax = MAX(amax,varlgout(i,j,k))
          END DO
        END DO
      END DO
      WRITE(6,'(7x,2(a,F15.5))') 'Write   MINVAL = ',amin,', MAXVAL = ',amax
    END IF

    CALL read_hdf_data_attr(fHandles(1,1),vname,comp_code,comp_prm,     &
                            comment, lcmnt, units, lunts, stag_dim,     &
                            istatus)

    CALL write_hdf_dataset(fHndlOut,nprocx_out,nprocy_out,vname,packed16,&
                  comp_code,comp_prm,comment,lcmnt,units,lunts,stag_dim,&
                  varlgout,nxlgout,nylgout,                             &
                  varout,nxout,nyout,nz,nzin,nzsoilin,nstypsin,istatus)
  END IF

  RETURN
END SUBROUTINE read_merge_write_real

!########################################################################
SUBROUTINE read_merge_write_real1d(IAMROOT,fHandles,nxpatch,nypatch,    &
                       vname,dims,direction,varin,varlgin,ni,nilg,iskip,&
                       fHndlOut,nprocx_out,nprocy_out,varout,varlgout,  &
                       niout,nilgout,dbglvl,istatus)
!-----------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: IAMROOT
  INTEGER, INTENT(IN) :: nxpatch, nypatch
  INTEGER, INTENT(IN) :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN) :: vname
  INTEGER, INTENT(IN) :: dims(6)
  CHARACTER(LEN=1), INTENT(IN) :: direction
  INTEGER, INTENT(IN) :: ni, nilg
  INTEGER, INTENT(IN) :: iskip

  INTEGER, INTENT(IN) :: nprocx_out, nprocy_out
  INTEGER, INTENT(IN) :: fHndlOut(nprocx_out,nprocy_out)
  INTEGER, INTENT(IN) :: niout, nilgout

  REAL,    INTENT(INOUT) :: varin(ni)
  REAL,    INTENT(INOUT) :: varlgin(nilg)
  REAL,    INTENT(INOUT) :: varout(niout)
  REAL,    INTENT(INOUT) :: varlgout(nilgout)

  INTEGER, INTENT(IN)    :: dbglvl
  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: comp_code, comp_prm, stag_dim
  INTEGER :: lcmnt, lunts
  CHARACTER(LEN=256) :: comment, units

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL read_hdf_dataset1d(fHandles,nxpatch,nypatch,vname,varin,           &
                          dims,direction,ni,istatus)

  IF (direction == 'x') THEN
    CALL mpimerge1dx(varin,ni,varlgin)
  ELSE IF (direction == 'y') THEN
    CALL mpimerge1dy(varin,ni,varlgin)
  ELSE
    istatus = -1
    RETURN
  END IF

  IF (IAMROOT) THEN
    IF (dbglvl > 2) THEN
      WRITE(6,'(7x,2(a,F15.5))') 'Read in MINVAL = ',MINVAL(varlgin),', MAXVAL = ',MAXVAL(varlgin)
    END IF

    DO i = 1,nilgout
      varlgout(i) = varlgin(i+iskip)
    END DO

    IF (dbglvl > 2) THEN
      WRITE(6,'(7x,2(a,F15.5))') 'Write   MINVAL = ',MINVAL(varlgout),', MAXVAL = ',MAXVAL(varlgout)
    END IF

    CALL read_hdf_data_attr(fHandles(1,1),vname,comp_code,comp_prm,     &
                            comment, lcmnt, units, lunts, stag_dim,     &
                            istatus)

    CALL write_hdf_dataset1d(fHndlOut,nprocx_out,nprocy_out,vname,      &
                  comment,lcmnt,units,lunts, direction,                 &
                  varlgout,nilgout,varout,niout,istatus)
  END IF

  RETURN
END SUBROUTINE read_merge_write_real1d

!########################################################################
SUBROUTINE copy_hdf_real1d(IAMROOT,fHandles,nxpatch,nypatch,            &
                       vname,dims,varin,nk,                             &
                       fHndlOut,nprocx_out,nprocy_out,varout, istatus)
!-----------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: IAMROOT
  INTEGER, INTENT(IN) :: nxpatch, nypatch
  INTEGER, INTENT(IN) :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN) :: vname
  INTEGER, INTENT(IN) :: dims(6)
  INTEGER, INTENT(IN) :: nk

  INTEGER, INTENT(IN) :: nprocx_out, nprocy_out
  INTEGER, INTENT(IN) :: fHndlOut(nprocx_out,nprocy_out)

  REAL,    INTENT(INOUT) :: varin(nk)
  REAL,    INTENT(INOUT) :: varout(nk)

  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------

  INTEGER :: k

  INTEGER :: comp_code, comp_prm, stag_dim
  INTEGER :: lcmnt, lunts
  CHARACTER(LEN=256) :: comment, units

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL read_hdf_dataset1d(fHandles,nxpatch,nypatch,vname,varin,           &
                          dims,'z',nk,istatus)


  IF (IAMROOT) THEN
    DO k = 1, nk
      varout(k) = varin(k)
    END DO

    CALL read_hdf_data_attr(fHandles(1,1),vname,comp_code,comp_prm,     &
                            comment, lcmnt, units, lunts, stag_dim,     &
                            istatus)

    CALL write_hdf_dataset1d(fHndlOut,nprocx_out,nprocy_out,vname,      &
                  comment,lcmnt,units,lunts, 'z',                       &
                  varout,nk,varout,nk,istatus)
  END IF

  RETURN
END SUBROUTINE copy_hdf_real1d

!########################################################################
SUBROUTINE read_merge_write_real4d(IAMROOT,fHandles,nxpatch,nypatch,    &
                         vname,rank,dims,dtype,varin,varlgin,           &
                         nx,ny,nzsoil,nstyp,nxlg,nylg,iskip,jskip,      &
                         fHndlOut,nprocx_out,nprocy_out,varout,varlgout,&
                         nxout,nyout,nxlgout,nylgout,istatus)
!-----------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: IAMROOT
  INTEGER, INTENT(IN) :: nxpatch, nypatch
  INTEGER, INTENT(IN) :: fHandles(nxpatch,nypatch)
  CHARACTER(LEN=256), INTENT(IN) :: vname
  INTEGER, INTENT(IN) :: rank, dims(6), dtype
  INTEGER, INTENT(IN) :: nx,ny,nzsoil,nstyp, nxlg, nylg
  INTEGER, INTENT(IN) :: iskip, jskip

  INTEGER, INTENT(IN) :: nprocx_out, nprocy_out
  INTEGER, INTENT(IN) :: fHndlOut(nprocx_out,nprocy_out)
  INTEGER, INTENT(IN) :: nxout, nyout, nxlgout, nylgout

  REAL,    INTENT(INOUT) :: varin(nx,ny,nzsoil,nstyp)
  REAL,    INTENT(INOUT) :: varlgin(nxlg,nylg,nzsoil,nstyp)
  REAL,    INTENT(INOUT) :: varout(nxout,nyout,nzsoil,nstyp)
  REAL,    INTENT(INOUT) :: varlgout(nxlgout,nylgout,nzsoil,nstyp)

  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k, n

  INTEGER :: comp_code, comp_prm, stag_dim
  INTEGER :: lcmnt, lunts
  CHARACTER(LEN=256) :: comment, units
  INTEGER :: packed16

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL read_hdf_dataset4d(fHandles,nxpatch,nypatch,vname,varin,           &
                          dims,nx,ny,nzsoil,nstyp,packed16,istatus)

  CALL mpimerge4d(varin,nx,ny,nzsoil,nstyp,varlgin)

  IF (IAMROOT) THEN
    DO n = 1, nstyp
      DO k = 1, nzsoil
        DO j = 1,nylgout
          DO i = 1,nxlgout
            varlgout(i,j,k,n) = varlgin(i+iskip,j+jskip,k,n)
          END DO
        END DO
      END DO
    END DO

    CALL read_hdf_data_attr(fHandles(1,1),vname,comp_code,comp_prm,     &
                            comment, lcmnt, units, lunts, stag_dim,     &
                            istatus)

    CALL write_hdf_dataset4d(fHndlOut,nprocx_out,nprocy_out,vname,packed16,&
                  comp_code,comp_prm,comment,lcmnt,units,lunts,stag_dim,&
                  varlgout,nxlgout,nylgout,                             &
                  varout,nxout,nyout,nzsoil,nstyp,istatus)
  END IF

  RETURN
END SUBROUTINE read_merge_write_real4d
