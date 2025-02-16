MODULE module_arbitrary_vario
  IMPLICIT NONE

  PRIVATE

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    To write arbitary array(s) into a file. The file can be in three format
!    based on the variable, file_fmt
!       file_fmt = 1,      Binary file
!       file_fmt = 3,      HDF4 file,
!         hdfcompr = 0-7, Specify HDF4 compression option, see arps.input
!       file_fmt = 7,      NetCDF file
!       others,           Error message
!
!  NOTE:
!    This module includes "mp.inc", So all variables in mp.inc must be
!    set before calling this subroutine. It also needs to link with the
!    following objects:
!
!           alloclib.o
!           mpisubs.o          or       nompsubs.o
!           hdfio3d.o                   nohdfio3d.o
!           netio3d.o                   nonetio3d.o
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (2009/09/20)
!  Based on wrtvar2.
!
!-----------------------------------------------------------------------

  INTEGER :: no_arbvar, current_var
  INTEGER :: file_fmt, hdf_compr
  INTEGER :: mpi_flag              ! > 0 Join and dump or read and split
                                   ! = 0 nompi or patch in and patch out
  CHARACTER(256) :: file_name
  INTEGER :: file_id

  LOGICAL :: file_exist
  LOGICAL :: file_opened

  INTEGER :: read_write_flag      ! 0 - File does not exist, to be writen
                                  ! 1 - file exist, to be read
  REAL    :: time
  CHARACTER(LEN=80)  :: runname
  CHARACTER(LEN=256) :: dirname
  LOGICAL :: filename_input

  LOGICAL :: initialized
  INTEGER :: current_time, no_times

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! EXPORTED INTERFACE
!
!-----------------------------------------------------------------------

  PUBLIC  :: arbvar_init     ! To be called before any other calls
  PUBLIC  :: arbvar_exit     ! Finalize the module
  PUBLIC  :: arbvar_read     ! Read interface
  PUBLIC  :: arbvar_wrt      ! Write interface
  PUBLIC  :: arbvar_readvar2 ! Provide convenient for traditional readvar2 calls
  PUBLIC  :: arbvar_advance_time

  CONTAINS

!#######################################################################

  SUBROUTINE arbvar_init(number_of_arbvar,form_in,hdfcompr_in,mpiflag_in, &
                         ioflag_in,filename_in,                         &
                         dirname_in,runname_in,first_varid,time_in,     &
                         dbglvl,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Initialize the arbitary var IO interface.
!
!    number_of_arbvar (INPUT) > 0  - will write out "number_of_arbvar"
!                                    arbitary arrays
!                             <=0  - To be determined dynamically
!
!  MODIFICATION:
!  2012/05/03 Keith Brewster, CAPS
!  Moved assignment of read_write_flag to outside of IF block.
!
!-----------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: number_of_arbvar
    INTEGER, INTENT(IN)  :: form_in, mpiflag_in,hdfcompr_in
    INTEGER, INTENT(IN)  :: ioflag_in
                                   ! Read write flag,
                                   ! = 1, To read, file must already exist.
                                   ! = 0, To write
    CHARACTER(LEN=*),   INTENT(IN) :: filename_in
    REAL,               INTENT(IN) :: time_in
    CHARACTER(LEN=*),   INTENT(IN) :: dirname_in,runname_in
    CHARACTER(LEN=6),   INTENT(IN) :: first_varid

    INTEGER, INTENT(IN)  :: dbglvl

    INTEGER, INTENT(OUT) :: istatus

!----------------------------------------------------------------------

    INTEGER :: lendir

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

!-------------------------- Setting file format ------------------------

    no_arbvar = 0
    IF (number_of_arbvar > 0) no_arbvar = number_of_arbvar

    IF (form_in == 1 .OR. form_in == 3 .OR. form_in == 7) THEN
      file_fmt = form_in
    ELSE
      WRITE(6,'(1x,a,I2,a,/,a)') 'Unsupported file format - ',form_in,  &
               ' in module module_arbitary_vario.',                     &
               'Please check the subroutine arguments for arbvar_init.'
      istatus = -1
      RETURN
    END IF

    IF (file_fmt == 3 ) THEN
      IF (hdfcompr_in >=0 .AND. hdfcompr_in <= 7) THEN
        hdf_compr = hdfcompr_in
      ELSE
        WRITE(6,'(1x,a,I2,a,/,a)') 'Argument hdfcompr_in = ',hdfcompr_in,&
                 ' in module module_arbitary_vario is not supported.',   &
                 'Please check the subroutine arguments for arbvar_init.'
        istatus = -2
        RETURN
      END IF
    END IF

    mpi_flag = mpiflag_in
    read_write_flag = ioflag_in

!-------------------------- Setting file name ------------------------

    IF ( LEN_TRIM(dirname_in) == 0 .OR. dirname_in == ' ' ) THEN
      dirname = './'
    ELSE
      dirname = TRIM(dirname_in)
    END IF

    lendir = LEN_TRIM(dirname)
    IF (dirname(lendir:lendir) /= '/') dirname(lendir+1:lendir+1) = '/'

    IF ( LEN_TRIM(filename_in) == 0 .OR. filename_in == ' ') THEN
      filename_input = .FALSE.

      runname = runname_in

      time = time_in

      CALL construct_file_name( first_varid, dbglvl, istatus )
    ELSE
      file_name = TRIM(dirname)//TRIM(filename_in)
      filename_input = .TRUE.
      IF (myproc == 0 .OR. mpi_flag == 0) THEN  ! no MPI mode, myproc == 0
        INQUIRE( FILE = TRIM(file_name), EXIST = file_exist )
      END IF
      IF (mpi_flag > 0) CALL mpupdatel(file_exist,1)
    END IF

!-------------------------- Initialing Module ------------------------

    current_var = 0
    file_opened = .FALSE.

    IF (read_write_flag > 0) THEN  ! to be read
      IF ( .NOT. file_exist ) THEN
        WRITE(6,'(1x,3a)') 'ERROR: File to be read - ',TRIM(file_name),' is not found.'
        istatus = -3
        RETURN
      END IF
    ELSE                           ! to be writen
      IF ( file_exist  ) THEN
        WRITE(6,'(1x,3a)') 'WARNING: File to be created - ',TRIM(file_name),' exists. Data will be appended.'
        istatus = 1
        RETURN
      END IF
    END IF

!-------------------------- Initialize Number of Times -----------------

   IF (form_in == 7 .AND. read_write_flag > 0 ) THEN
     CALL get_number_times( no_times, istatus )
   ELSE
     no_times = 0
   END IF
   current_time = 1

!-------------------------- Initialization Flag ------------------------

    initialized = .TRUE.

    IF (dbglvl > 2) WRITE(6,'(1x,3a)') 'MODULE module_arbitary_vario is initialized.'

    RETURN
  END SUBROUTINE arbvar_init

!#######################################################################

  SUBROUTINE arbvar_advance_time(time_in, dbglvl,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Initialize the arbitary var IO interface.
!
!    number_of_arbvar (INPUT) > 0  - will write out "number_of_arbvar"
!                                    arbitary arrays
!                             <=0  - To be determined dynamically
!
!-----------------------------------------------------------------------

    REAL,               INTENT(IN) :: time_in

    INTEGER, INTENT(IN)  :: dbglvl

    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

!-------------------------- Setting file format ------------------------

    IF (read_write_flag > 0) THEN
      current_time = current_time + 1
      time         = time_in
    ELSE
      ! WRITE file do nothing at present
    END IF

    IF (dbglvl > 2) WRITE(6,'(1x,a,I6)') 'Current time advanced to ',   &
                                         current_time

    RETURN
  END SUBROUTINE arbvar_advance_time

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARBVAR_WRT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  SUBROUTINE arbvar_wrt(nx,ny,nz,var,varid,vardesc,varunits,dbglvl,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To write an array 'var' into a file. The file can be in three format
!
!  AUTHOR: Yunheng Wang (2005/06/14)
!  Based on wrtvar1. The interface is the same as wrtvar1 except that
!  two more arguments, foutfmt, hdfcompr were added.
!
!  MODIFICATION HISTORY:
!
!  INPUT:
!
!    nx, ny, nz   Dimensions of input array 'var'.
!    var          Array to be written out.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!    vardesc      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!
!  OUTPUT:
!
!    istatus       Exit status (0=okay, 1=unknown file format, 2=error)
!
!-----------------------------------------------------------------------
!
    INTEGER,          INTENT(IN)  :: nx,ny,nz
    REAL,    TARGET,  INTENT(IN)  :: var(nx,ny,nz)
    CHARACTER(LEN=6), INTENT(IN)  :: varid
    CHARACTER(LEN=*), INTENT(IN)  :: varunits
    CHARACTER(LEN=*), INTENT(IN)  :: vardesc
    INTEGER,          INTENT(IN)  :: dbglvl
    INTEGER,          INTENT(OUT) :: istatus

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
    CHARACTER(LEN=40)  :: varunits1
    CHARACTER(LEN=40)  :: vardesc1

    INTEGER            :: nxlg, nylg
    REAL,     POINTER  :: varout(:,:,:)

    LOGICAL            :: varout_allocated

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
    istatus = 0

    !
    ! Merge variables if needed
    !

    IF (mp_opt > 0 .AND. mpi_flag > 0) THEN
      nxlg = (nx-3)*nproc_x + 3
      nylg = (ny-3)*nproc_y + 3
    ELSE
      nxlg = nx
      nylg = ny
    END IF

    IF (mp_opt > 0 .AND. mpi_flag > 0) THEN
      ALLOCATE(varout(nxlg,nylg,nz), STAT= istatus)
      CALL mpimerge3d(var,nx,ny,nz,varout)
      varout_allocated = .TRUE.
    ELSE
      varout => var
      varout_allocated = .FALSE.
    END IF

    IF ( .NOT. file_opened ) THEN
      IF ( .NOT. filename_input ) CALL construct_file_name( varid, dbglvl, istatus )
      CALL open_file(0,dbglvl,istatus)
      current_var = 0
    END IF

    vardesc1  = vardesc
    varunits1 = varunits

    CALL write_data(nxlg,nylg,nz,varout,varid,vardesc1,varunits1,dbglvl,istatus)

    current_var = current_var + 1

    IF (current_var >= no_arbvar ) CALL close_file(dbglvl,istatus)

    CALL mpupdatei(istatus, 1)

    IF ( varout_allocated ) DEALLOCATE(varout)

    RETURN
  END SUBROUTINE arbvar_wrt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARBVAR_READ                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

  SUBROUTINE arbvar_read(nx,ny,nz,var, varid,vardesc,varunits,dbglvl,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To read in array 'var' from a file.
!
!  AUTHOR:
!    Yunheng Wang (2005/06/14)
!    Based on readvar1.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!  INPUT:
!
!    nx, ny, nz   Dimensions of input array 'var'.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!
!  OUTPUT:
!
!    var          Array to be written out.
!    vardesc      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!    istatus       Exit status (0=okay, 1=warning, 2=error)
!
!-----------------------------------------------------------------------

    INTEGER,          INTENT(IN)  :: nx,ny,nz
    CHARACTER(LEN=6), INTENT(IN)  :: varid
    CHARACTER(LEN=*), INTENT(OUT) :: vardesc, varunits
    REAL,   TARGET,   INTENT(OUT) :: var(nx,ny,nz)

    INTEGER,          INTENT(IN)  :: dbglvl
    INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    CHARACTER(LEN=40)  :: vardesc1, varunits1

    INTEGER            :: nxlg, nylg
    REAL, POINTER      :: varin(:,:,:)
    LOGICAL            :: varin_allocated
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
    istatus = 0

    IF (mp_opt > 0 .AND. mpi_flag > 0) THEN   ! MPI mode, read and split.
      nxlg = (nx-3)*nproc_x + 3
      nylg = (ny-3)*nproc_y + 3
      ALLOCATE(varin(nxlg,nylg,nz), STAT=istatus)
      CALL check_alloc_status(istatus, 'ARBVAR:varin')
      varin_allocated = .TRUE.
    ELSE                    ! no MPI mode or all process read.
      nxlg = nx
      nylg = ny
      varin => var
      varin_allocated = .FALSE.
    END IF

    varunits  = ' '
    varunits1 = ' '
    vardesc1  = ' '

    8000 CONTINUE

    IF ( .NOT. file_opened ) THEN
      IF ( .NOT. filename_input ) CALL construct_file_name ( varid, dbglvl, istatus )
      CALL open_file( 1, dbglvl, istatus )
      IF (istatus /= 0) RETURN
      current_var = 0
    END IF

    CALL read_data(nxlg,nylg,nz,varid,varin,vardesc1,varunits1,dbglvl,istatus)

    CALL mpupdatei(istatus,1)

    IF (istatus == 0) THEN                     ! successful read
      IF (mp_opt > 0 .AND. mpi_flag > 0) THEN   ! split data
        CALL mpisplit3d(varin,nx,ny,nz,var)
      END IF
      vardesc  = vardesc1
      varunits = varunits1
    ELSE IF (istatus == 9999 .AND. INDEX(file_name,varid) <= 0 .AND. (.NOT. filename_input) ) THEN
                                           ! Variable is not in the current file
      CALL close_file( dbglvl, istatus )
      GOTO 8000
    ELSE
      CALL close_file( dbglvl, istatus )
      IF(myproc == 0) WRITE(6,'(1x,4a)') 'ERROR: reading ',varid,' from ',TRIM(file_name)
      RETURN
    END IF

    current_var = current_var + 1

    IF (current_var >= no_arbvar) CALL close_file( dbglvl, istatus )

    IF ( varin_allocated ) DEALLOCATE(varin)

    RETURN
  END SUBROUTINE arbvar_read

!######################################################################

  SUBROUTINE arbvar_exit(dbglvl,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Finalize the arbitary var IO interface.
!
!-----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    istatus = 0

    IF ( file_opened ) CALL close_file(dbglvl,istatus)

    current_var = 0
    no_arbvar   = 0
    no_times     = 0
    current_time = 0
    file_exist  = .FALSE.
    initialized = .FALSE.

    IF (dbglvl > 2) WRITE(6,'(1x,3a)') 'MODULE module_arbitary_vario exited.'

    RETURN
  END SUBROUTINE arbvar_exit

!######################################################################

  SUBROUTINE open_file(ioflag,dbglvl,istatus)
!----------------------------------------------------------------------
!  Open file, once the file is opend, it should also exist.
!  So module variable "file_exist" is set together with "file_opened".
!----------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: ioflag
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (myproc == 0 .OR. mpi_flag == 0) THEN
      IF (file_fmt == 1) THEN
        CALL getunit( file_id )

        CALL asnctl ('NEWLOCAL', 1, istatus)
        CALL asnfile(file_name, '-F f77 -N ieee', istatus)

        IF ( ioflag > 0) THEN
          OPEN(UNIT=file_id,FILE=trim(file_name),STATUS='old',          &
               FORM='unformatted', ERR=9000, IOSTAT=istatus)
        ELSE
          IF (file_exist) THEN
            OPEN(UNIT=file_id,FILE=TRIM(file_name),STATUS='unknown',    &
                 FORM='unformatted', ERR=9000, IOSTAT=istatus)
            !fseek(file_id,0,2)  ! End of file
            DO WHILE (istatus == 0)  ! Go to the END of file
              READ(file_id,END=8000,IOSTAT=istatus)
            END DO
            8000 CONTINUE
            istatus = 0
          ELSE
            OPEN(UNIT=file_id,FILE=TRIM(file_name),STATUS='new',        &
                 FORM='unformatted', ERR=9000, IOSTAT=istatus)
          END IF
        END IF

        9000 CONTINUE

      ELSE IF (file_fmt == 3) THEN

        IF ( ioflag > 0) THEN   ! to be read
          CALL hdfopen(file_name,1,file_id)
        ELSE
          IF (file_exist) THEN  ! to be appended
            CALL hdfopen(file_name,3,file_id)
          ELSE                  ! to be created
            CALL hdfopen(file_name,2,file_id)
          END IF
        END IF

        IF (file_id <= 0) istatus = file_id

      ELSE IF (file_fmt == 7) THEN

        IF ( ioflag > 0) THEN   ! to be read
          CALL netopen(file_name,'R',file_id)
        ELSE
          IF (file_exist) THEN  ! to be appended
            CALL netopen(file_name,'A',file_id)
          ELSE                  ! to be created
            CALL netopen(file_name,'C',file_id)
          END IF
        END IF

        IF (file_id <= 0) istatus = file_id

      END IF

    END IF

    CALL mpupdatei(istatus,1)
    IF (istatus == 0) THEN
      file_opened = .TRUE.
      file_exist  = .TRUE.

      IF (dbglvl > 0 .AND. myproc == 0)                                 &
        WRITE(6,'(1x,3a)') 'File - ',TRIM(file_name),' opened.'
    ELSE
      IF (myproc == 0) WRITE(6,'(1x,3a,I4)') 'ERROR: File ',TRIM(file_name),' opening error - ',istatus
    END IF

    RETURN
  END SUBROUTINE open_file

!######################################################################

  SUBROUTINE close_file(dbglvl,istatus)

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (myproc == 0 .OR. mpi_flag == 0) THEN

      IF (file_fmt == 1) THEN

        CLOSE(UNIT=file_id)
        CALL retunit(file_id)

      ELSE IF (file_fmt == 3) THEN

        CALL hdfclose(file_id,istatus)

      ELSE IF (file_fmt == 7) THEN

        CALL netclose(file_id)

      END IF

    END IF

    CALL mpupdatei(istatus,1)
    IF (istatus == 0) THEN
      file_opened = .FALSE.
      IF (dbglvl > 0 .AND. myproc == 0)                                 &
        WRITE(6,'(1x,3a)') 'Aribtary variable file - ',TRIM(file_name),' closed.'
    ELSE
      IF (myproc == 0) WRITE(6,'(1x,3a,I4)') 'ERROR: File ',TRIM(file_name),' closing error - ',istatus
    END IF

    RETURN
  END SUBROUTINE close_file

!######################################################################

  SUBROUTINE get_number_times( no_times, istatus )
    INTEGER, INTENT(OUT) :: no_times
    INTEGER, INTENT(OUT) :: istatus

    no_times = 0
    IF (file_fmt == 7 .AND. myproc == 0) THEN
      CALL net_get_unlimit_size( file_name, no_times, istatus)
    END IF

    CALL mpupdatei(no_times,1)

    RETURN
  END SUBROUTINE

!######################################################################

  SUBROUTINE write_data(nx,ny,nz,vardata,varname,vardesc,varunits,dbglvl,istatus)

    INTEGER, INTENT(IN)  :: nx, ny, nz
    REAL,    INTENT(IN)  :: vardata(nx,ny,nz)
    CHARACTER(LEN=6),  INTENT(IN) :: varname
    CHARACTER(LEN=40), INTENT(IN) :: vardesc, varunits
    INTEGER, INTENT(IN)  :: dbglvl

    INTEGER, INTENT(OUT) :: istatus

!----------------------------------------------------------------------

    INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
    REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (myproc == 0 .OR. mpi_flag == 0) THEN

      PRINT *, 'ARBVAR: Writing array ',TRIM(varname),' to file ', TRIM(file_name)

      IF (file_fmt == 1) THEN

        WRITE(file_id, ERR=9000, IOSTAT=istatus) nx,ny,nz
        WRITE(file_id, ERR=9000, IOSTAT=istatus) vardata
        WRITE(file_id, ERR=9000, IOSTAT=istatus) vardesc
        WRITE(file_id, ERR=9000, IOSTAT=istatus) varunits

        GO TO 8000

        ! I/O error handling
        ! Note that IOSTAT < 0 should not occur in this subroutine.

        9000 CONTINUE        ! I/O errors

        IF (istatus < 0) THEN
          PRINT *, 'ARBVAR: I/O ERRORS OCCURRED ',                      &
                   '(possible end of record or file): ', istatus
        ELSE IF (istatus > 0) THEN
          PRINT *, 'ARBVAR: UNRECOVERABLE I/O ERRORS OCCURRED: ',istatus
        END IF

        8000 CONTINUE

      ELSE IF (file_fmt == 3) THEN

        IF (hdf_compr > 3) THEN
          ALLOCATE (itmp(nx,ny,nz),stat=istatus)
          CALL check_alloc_status(istatus, 'ARBVAR:itmp')
          ALLOCATE (hmax(nz),stat=istatus)
          CALL check_alloc_status(istatus, 'ARBVAR:hmax')
          ALLOCATE (hmin(nz),stat=istatus)
          CALL check_alloc_status(istatus, 'ARBVAR:hmin')
        END IF

        CALL hdfwrti(file_id, 'nx', nx, istatus)
        CALL hdfwrti(file_id, 'ny', ny, istatus)
        CALL hdfwrti(file_id, 'nz', nz, istatus)
        CALL hdfwrt3d(vardata,nx,ny,nz,file_id,0,hdf_compr,             &
                      varname,vardesc,varunits,itmp,hmax,hmin)

        IF (hdf_compr > 3) THEN
          DEALLOCATE(itmp)
          DEALLOCATE(hmax,hmin)
        END IF

      ELSE IF (file_fmt == 7) THEN

        CALL net_define_onevar(file_id,nx,ny,nz,varname,vardesc,varunits,istatus)
        CALL netwrt3d(file_id,0,0,varname,vardata,nx,ny,nz)

      END IF

    END IF

    RETURN
  END SUBROUTINE write_data

!######################################################################

  SUBROUTINE read_data(nx,ny,nz,varname,vardata,vardesc,varunits,dbglvl,istatus)

    INTEGER,           INTENT(IN)  :: nx,ny,nz
    CHARACTER(LEN=6),  INTENT(IN)  :: varname
    REAL,              INTENT(OUT) :: vardata(nx,ny,nz)
    CHARACTER(LEN=40), INTENT(OUT) :: vardesc, varunits

    INTEGER,           INTENT(IN)  :: dbglvl
    INTEGER,           INTENT(OUT) :: istatus

!----------------------------------------------------------------------

    INTEGER   :: nx_in,ny_in,nz_in
    INTEGER   :: itime

    INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
    REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF (myproc == 0 .OR. mpi_flag == 0) THEN  ! no MPI mode, myproc == 0
                                             ! MPI mode,  mpiflag = 1 --- Only root process reads
                                             !            mpiflag = 0 --- all processes read
      WRITE(6,'(/,1x,4a)') 'ARBVAR: Reading array ',varname,' from file ',trim(file_name)

      IF (file_fmt == 1) THEN

        READ(file_id, ERR=9000, END=9000, IOSTAT=istatus) nx_in,ny_in,nz_in

        IF(nx_in /= nx .OR. ny_in /= ny .OR. nz_in /= nz) THEN
          WRITE(6,'(a,/a,a,/a,3I5,/a,3I5)')                                   &
              'Warning in subroutine arbvar_read: Dimensions of data file ',  &
              file_name,' do not agree with the expected dimensions.',        &
              'nx, ny and nz in the data are ',nx_in,ny_in,nz_in,             &
              'nx, ny and nz expected    are ',nx,ny,nz
          istatus = 2
          GO TO 9000
        END IF

        READ(file_id, ERR=9000,END=9000,IOSTAT=istatus) vardata
        READ(file_id, ERR=9000,END=9000,IOSTAT=istatus) vardesc
        READ(file_id, ERR=9000,END=9000,IOSTAT=istatus) varunits

        GO TO 8000

        ! I/O error handling
        ! Note that warning (istatus=1) (e.g., end of record or file) is implemented as
        ! error (status=2) because of ambiguities in IOSTAT.

        9000 CONTINUE        ! I/O errors

        IF (istatus < 0) THEN
          PRINT *, 'ARBVAR: I/O ERRORS OCCURRED ',                          &
                   '(possible end of record or file): ', istatus
        ELSE IF (istatus > 0) THEN
          PRINT *, 'ARBVAR: UNRECOVERABLE I/O ERRORS OCCURRED: ', istatus
        END IF
        istatus = 9999

        8000 CONTINUE

      ELSE IF (file_fmt == 3) THEN

        ALLOCATE (itmp(nx,ny,nz),stat=istatus)
        CALL check_alloc_status(istatus, "arbvar_read:itmp")
        ALLOCATE (hmax(nz),stat=istatus)
        CALL check_alloc_status(istatus, "arbvar_read:hmax")
        ALLOCATE (hmin(nz),stat=istatus)
        CALL check_alloc_status(istatus, "arbvar_read:hmin")

        CALL hdfrd3d(file_id,varname,nx,ny,nz,vardata,istatus,itmp,hmax,hmin)
        IF (istatus == 1) THEN  ! Variable is not found
          istatus = 9999
          RETURN
        END IF

        CALL get_var_attr_from_hdf(file_id,varname,'comment',vardesc, 40,istatus)
        CALL get_var_attr_from_hdf(file_id,varname,'units',  varunits,40,istatus)

        DEALLOCATE(itmp)
        DEALLOCATE(hmax,hmin)

      ELSE IF (file_fmt == 7) THEN

        CALL net_get_onevar(file_id,nx_in,ny_in,nz_in,varname,vardesc,varunits,istatus)
        IF (istatus == 9999) RETURN

        itime = 0
        IF (no_times > 0) itime = current_time
        CALL netread3d(file_id,0,itime,varname,nx_in,ny_in,nz_in,vardata)

      END IF

    END IF

  !  CALL dbgwrt(nprocs*100+myproc,var,nx,ny,nz,2,nx-2,2,ny-2,1,1,istatus)

    RETURN
  END SUBROUTINE read_data

!#######################################################################

  SUBROUTINE construct_file_name( varid, dbglvl, istatus )

!-------------------------- Setting file name ------------------------
!
! PURPOSE:
!
!    Construct file name and determine whether it exists or not.
!
!    An file named  dirname/runname.xxxtime will be created or be read.
!
!    where xxx is
!            bin        for foutfmt = 1
!            hdf        for foutfmt = 3
!            net        for foutfmt = 7
!
!   The state of whether file exists or not is denoted by module variable:
!            file_exist
!
!-----------------------------------------------------------------------

    CHARACTER(LEN=6), INTENT(IN)  :: varid
    INTEGER,          INTENT(IN)  :: dbglvl
    INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    CHARACTER(LEN=256) :: tmpstr
    CHARACTER(LEN=256) :: timestring
    INTEGER            :: lfnkey
    INTEGER            :: ltimestring, lenfil

    CHARACTER(LEN=3), PARAMETER :: fmtstr(7) =                          &
                            (/'bin','xxx','hdf','xxx','xxx','xxx','net'/)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    istatus = 0

    file_name = ' '

    CALL gtlfnkey( runname, lfnkey )

    WRITE (file_name,'(3A,A3,A6)') TRIM(dirname),runname(1:lfnkey),'.', &
                               fmtstr(file_fmt),varid
    lenfil = len_trim(file_name)

    CALL cvttsnd( time, timestring, ltimestring )

    file_name(lenfil+1:lenfil+ltimestring) = timestring(1:ltimestring)
    lenfil = lenfil + ltimestring

    IF (mp_opt > 0 .AND. mpi_flag <= 0) THEN
      tmpstr = file_name
      CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,0,0,read_write_flag,dbglvl,file_name,istatus)
      lenfil = LEN_TRIM(file_name)
    END IF

    IF (read_write_flag == 0) CALL fnversn(file_name,lenfil)

    INQUIRE( FILE = TRIM(file_name), EXIST = file_exist )

    RETURN
  END SUBROUTINE construct_file_name

!#######################################################################

  SUBROUTINE arbvar_readvar2(nx,ny,nz,vardata,varid,vardesc,varunits,   &
                      timein,runnamein,dirnamein,filenamein,            &
                      finfmt,mpiflag,dbglvl,istatus)
!-----------------------------------------------------------------------
! Subroutine to be compatible with readvar2
!
!-----------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: nx, ny, nz
    REAL,    INTENT(OUT) :: vardata(nx,ny,nz)
    CHARACTER(LEN=6), INTENT(IN)  :: varid
    CHARACTER(LEN=*), INTENT(OUT) :: vardesc, varunits
    CHARACTER(LEN=*), INTENT(IN)  :: runnamein,dirnamein,filenamein
    REAL,    INTENT(IN)  :: timein
    INTEGER, INTENT(IN)  :: finfmt, mpiflag, dbglvl

    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF (initialized .AND. current_time < no_times) THEN
      IF (timein > time) CALL arbvar_advance_time(timein, dbglvl, istatus)
    ELSE
      CALL arbvar_init(1,finfmt,0,mpiflag,1,filenamein,                 &
                       dirnamein,runnamein,varid,timein,dbglvl,istatus)
    END IF

    CALL arbvar_read(nx,ny,nz,vardata,varid,vardesc,varunits,dbglvl,istatus)

    IF (no_times < 1 .OR. ( current_time >= no_times .AND. current_var >= no_arbvar) ) THEN
      ! not multiple time in one file or
      ! All time levels has been readed.
      CALL arbvar_exit( dbglvl,istatus )
    END IF

    RETURN
  END SUBROUTINE arbvar_readvar2


END MODULE module_arbitrary_vario
