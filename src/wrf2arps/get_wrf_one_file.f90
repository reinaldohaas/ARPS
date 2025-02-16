!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE open_wrf_file              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE open_wrf_one_file(filename,io_form,nidout,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Open a WRF file and return NetCDF file handler. The file handler
!    will only be returned to root processor (myproc == 0). It should
!    be used in no-mpi mode or mpi mode reading joined WRF file.
!
!    NOTE: it is required to call close_wrf_file explicitly to close
!          the opened file in your calling program.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(OUT) :: nidout
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  LOGICAL            :: fexists

  CHARACTER(LEN=80)  :: sysdepinfo
  LOGICAL, SAVE      :: initialized = .FALSE.

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  INQUIRE(FILE = TRIM(filename), EXIST = fexists)
  IF ( .NOT. fexists ) THEN
    WRITE(6,'(3a)') 'File not found: ',filename,' in open_wrf_file'
    CALL arpsstop('WRF file not exist.',1)
  ENDIF

  sysdepinfo = 'DATASET=HISTORY'

  IF (io_form == 1) THEN            ! initialize explicitly

    IF (.NOT. initialized) CALL ext_int_ioinit( SysDepInfo, iStatus )

    IF (myproc == 0)  CALL ext_int_open_for_read(filename, 0, 0,        &
                                        SysDepInfo, nidout, iStatus )

  ELSE IF (io_form == 5) THEN       ! initialized inside open_phdf5_for_read

    CALL open_phdf5_for_read( filename, sysdepinfo, initialized, nidout, iStatus)

  ELSE IF (io_form == 7) THEN           ! no initialization needed

    IF (myproc == 0 ) CALL open_ncd_file(filename,nidout)

  ELSE
    WRITE(0,*) 'Unsupported IO format - ',io_form,'.'
    CAlL arpsstop('Unsupported IO format.',1)
  END IF

  initialized = .TRUE.

  RETURN
END SUBROUTINE open_wrf_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE close_wrf_file               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE close_wrf_one_file(nch,io_form,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!     Close the WRF file which is opened using open_wrf_one_file.
!     Only root processor do the job.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nch
  INTEGER, INTENT(IN) :: io_form
  INTEGER, INTENT(OUT):: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

  IF (io_form == 1) THEN
    IF (myproc == 0) CALL ext_int_ioclose(nch,iStatus)
  ELSE IF (io_form == 5) THEN
    CALL close_phdf5_for_read(nch,istatus)
  ELSE IF(io_form == 7) THEN
    IF (myproc == 0) CALL close_ncd_file(nch)
  END IF

  IF (istatus /= 0) THEN
    WRITE(0,'(1x,2a)') 'ERROR: closing file handler ',nch
    CALL arpsstop('Error in close_wrf_one_file.',1)
  END IF

  RETURN
END SUBROUTINE close_wrf_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE get_wrf_Times                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_Times_from_one_file(nfid,io_form,itime,timestr,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read Date String (Times) in the WRF outputs at specified time.
!    And broadcast to all processors if necessary.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nfid     ! File handler
  INTEGER, INTENT(IN)  :: io_form  ! File format
  INTEGER, INTENT(IN)  :: itime    ! Time dimension value
                                   ! this is the unlimited dimension
  CHARACTER(LEN=*), INTENT(OUT) :: timestr
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (io_form == 5) THEN

     CAlL get_phdf5_next_time(nfid,timestr,istatus)

  ELSE       ! need to broadcast

    IF (io_form == 1) THEN
      IF (myproc == 0) CALL ext_int_get_next_time(nfid,timestr,istatus)
    ELSE IF (io_form == 7) THEN
      IF (myproc == 0) CALL get_ncd_next_time(nfid,itime,timestr,istatus)
    END IF
    CALL mpupdatec(timestr,LEN(timestr))

  END IF

  RETURN
END SUBROUTINE get_wrf_Times_from_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_wrf_metadata              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_meta_from_one_file(nid,io_form,                  &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj,trlat1,trlat2,trlon,ctrlat,ctrlon,      &
                          dx,dy,dt,sfcphys,sfclay,mpphys,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve WRF grib information from the NetCDF file which are stored
!   as Global attributes. Only root processor do the job and broadcast
!   to other processors if necessary.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nid
  INTEGER, INTENT(IN)  :: io_form
  INTEGER, INTENT(OUT) :: nx_ext, ny_ext     ! they are whole domain dimensions
  INTEGER, INTENT(OUT) :: nz_ext, nzsoil_ext
  INTEGER, INTENT(OUT) :: iproj              ! WRF map projection
  REAL,    INTENT(OUT) :: trlat1
  REAL,    INTENT(OUT) :: trlat2
  REAL,    INTENT(OUT) :: trlon
  REAL,    INTENT(OUT) :: ctrlat
  REAL,    INTENT(OUT) :: ctrlon
  REAL,    INTENT(OUT) :: dx
  REAL,    INTENT(OUT) :: dy
  REAL,    INTENT(OUT) :: dt
  INTEGER, INTENT(OUT) :: sfcphys, sfclay
  INTEGER, INTENT(OUT) :: mpphys
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variable
!
!-----------------------------------------------------------------------

  INTEGER           :: ips,  ipe,  jps,  jpe
  INTEGER           :: ips_u,ipe_u,jps_u,jpe_u

  CHARACTER(LEN=80) :: cdump
  INTEGER           :: idump
  REAL              :: rdump

  INCLUDE   'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL get_dom_ti_char_one_file(nid,io_form,'TITLE',     cdump,istatus)
  CALL get_dom_ti_char_one_file(nid,io_form,'START_DATE',cdump,istatus)

  CALL get_dom_ti_integer_one_file(nid,io_form,'WEST-EAST_GRID_DIMENSION',  nx_ext,istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SOUTH-NORTH_GRID_DIMENSION',ny_ext,istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'BOTTOM-TOP_GRID_DIMENSION', nz_ext,istatus)

  CALL get_dom_ti_char_one_file   (nid,io_form,'GRIDTYPE',cdump,istatus)
!  CALL get_dom_ti_integer_one_file(nid,io_form,'DYN_OPT', idump,istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'DIFF_OPT',idump,istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'KM_OPT',  idump,istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'DAMP_OPT',idump,istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'KHDIF',   rdump,istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'KVDIF',   rdump,istatus)

  CALL get_dom_ti_integer_one_file(nid,io_form,'MP_PHYSICS',        mpphys, istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'RA_LW_PHYSICS',     idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'RA_SW_PHYSICS',     idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SF_SFCLAY_PHYSICS', sfclay, istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SF_SURFACE_PHYSICS',sfcphys,istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'BL_PBL_PHYSICS',    idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'CU_PHYSICS',        idump,  istatus)

  CALL get_dom_ti_integer_one_file(nid,io_form,'WEST-EAST_PATCH_START_UNSTAG',  ips_u,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'WEST-EAST_PATCH_END_UNSTAG',    ipe_u,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'WEST-EAST_PATCH_START_STAG',    ips,    istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'WEST-EAST_PATCH_END_STAG',      ipe,    istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SOUTH-NORTH_PATCH_START_UNSTAG',jps_u,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SOUTH-NORTH_PATCH_END_UNSTAG',  jpe_u,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SOUTH-NORTH_PATCH_START_STAG',  jps,    istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'SOUTH-NORTH_PATCH_END_STAG',    jpe,    istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'BOTTOM-TOP_PATCH_START_UNSTAG', idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'BOTTOM-TOP_PATCH_END_UNSTAG',   idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'BOTTOM-TOP_PATCH_START_STAG',   idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'BOTTOM-TOP_PATCH_END_STAG',     idump,  istatus)

  CALL get_dom_ti_real_one_file   (nid,io_form,'DX',dx,    istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'DY',dy,    istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'DT',dt,    istatus)

  CALL get_dom_ti_real_one_file   (nid,io_form,'CEN_LAT',      ctrlat, istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'CEN_LON',      ctrlon, istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'TRUELAT1',     trlat1, istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'TRUELAT2',     trlat2, istatus)

  CALL get_dom_ti_real_one_file   (nid,io_form,'MOAD_CEN_LAT', rdump,  istatus)
  CALL get_dom_ti_real_one_file   (nid,io_form,'STAND_LON',    trlon,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'MAP_PROJ',     iproj,  istatus)

  CALL get_dom_ti_char_one_file   (nid,io_form,'MMINLU',       cdump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'ISWATER',      idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'ISICE',        idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'ISURBAN',      idump,  istatus)
  CALL get_dom_ti_integer_one_file(nid,io_form,'ISOILWATER',   idump,  istatus)

!-----------------------------------------------------------------------
!
!  Determine soil layers from surface physics option
!
!-----------------------------------------------------------------------

  IF (sfcphys == 1) THEN
    nzsoil_ext = 5
  ELSE IF (sfcphys == 2) THEN
    nzsoil_ext = 4
  ELSE IF (sfcphys == 3) THEN
    nzsoil_ext = 6
  ELSE IF (sfcphys == 7) THEN
    nzsoil_ext = 2
  ELSE
    WRITE(6,*) '=============================================='
    WRITE(6,*) 'WARNING: unknown sf_surface_physics = ',sfcphys
    WRITE(6,*) '=============================================='
    nzsoil_ext = 5
  END IF

  nzsoil_ext = nzsoil_ext + 1  ! Use surface as an extra soil layer

  RETURN
END SUBROUTINE get_wrf_meta_from_one_file

SUBROUTINE get_dom_ti_integer_one_file(nid,io_form,element, val, ireturn)
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: element
  INTEGER,      INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'mp.inc'

  INTEGER    :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 5) THEN
    CALL get_phdf5_dom_ti_integer(nid,element,val,ireturn)
  ELSE       ! broadcast necessary
    IF (io_form == 7) THEN
      IF (myproc == 0)    &
         CALL get_ncd_dom_ti_integer(nid,element,val,ireturn)
    ELSE IF (io_form == 1) THEN
      IF (myproc == 0)    &
         CALL ext_int_get_dom_ti_integer(nid,element,val,1,outcount,ireturn)
    END IF
    CALL mpupdatei(val,1)
  END IF

  RETURN
END SUBROUTINE get_dom_ti_integer_one_file

SUBROUTINE  get_dom_ti_real_one_file(nid,io_form,element, val, ireturn)
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: element
  REAL,         INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'mp.inc'

  INTEGER   :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 5) THEN
    CALL get_phdf5_dom_ti_real(nid,element,val,ireturn)
  ELSE
    IF (io_form == 7) THEN
      IF (myproc == 0)    &
         CALL get_ncd_dom_ti_real(nid,element,val,ireturn)
    ELSE IF (io_form == 1) THEN
      IF (myproc == 0)    &
         CALL ext_int_get_dom_ti_real(nid,element,val,1, outcount,ireturn)
    END IF
    CALL mpupdater(val,1)
  END IF

  RETURN
END SUBROUTINE get_dom_ti_real_one_file
!
SUBROUTINE  get_dom_ti_char_one_file(nid,io_form,element, val, ireturn)
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: element
  CHARACTER(*), INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'mp.inc'

  INTEGER   :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 5) THEN
    CALL get_phdf5_dom_ti_char(nid,element,val,ireturn)
  ELSE
    IF (io_form == 7) THEN
      IF (myproc == 0)   &
         CALL get_ncd_dom_ti_char(nid,element,val,ireturn)
    ELSE IF (io_form == 1) THEN
      IF (myproc == 0)   &
         CALL ext_int_get_dom_ti_char(nid,element,val,ireturn)
    END IF
    CALL mpupdatec(val,LEN(val))
  END IF

  RETURN
END SUBROUTINE get_dom_ti_char_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_wrf_dummy                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_dummy_from_one_file(nid,io_form,datestr,itime,varname,   &
                varType,memoryorder,stagger,dimname1,dimname2,dimname3, &
                nx,ny,nz,nxd,nyd,nzd,temtd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in an array from the WRF history file. It just for sequential
!    access of WRF binary file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nid
  INTEGER,          INTENT(IN)  :: io_form
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: varType
  CHARACTER(LEN=*), INTENT(IN)  :: MemoryOrder
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  CHARACTER(LEN=*), INTENT(IN)  :: dimname2
  CHARACTER(LEN=*), INTENT(IN)  :: dimname3
  INTEGER,          INTENT(IN)  :: nx              ! local index
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  INTEGER,          INTENT(IN)  :: nxd,nyd,nzd         ! domain index
  REAL,             INTENT(OUT) :: temtd(nxd*nyd*nzd)   ! domain array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: xdim, ydim, zdim
  INTEGER       :: nxlg, nylg
  INTEGER       :: ilocs,iloce,jlocs,jloce

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF ( io_form /= 1 ) RETURN   ! Only works for binary format

  nxlg = (nx-1)*nproc_x
  nylg = (ny-1)*nproc_y

  IF (stagger == 'X') nxlg = nxlg + 1
  IF (stagger == 'Y') nylg = nylg + 1

  IF (MemoryOrder == 'XZY') THEN
    xdim = 1
    ydim = 3
    zdim = 2
  ELSE IF (MemoryOrder == 'XY' .OR. MemoryOrder == 'XYZ') THEN
    xdim = 1
    ydim = 2
    zdim = 3
  ELSE
    xdim = 1
    ydim = 2
    zdim = 3

    nxlg = nxd
    nylg = nyd
  END IF

  DimNames(xdim) = dimname1
  DimNames(ydim) = dimname2
  DimNames(zdim) = dimname3

  DomainStart(:)  = 1
  DomainEnd(xdim) = nxlg
  DomainEnd(ydim) = nylg
  DomainEnd(zdim) = nzd

  PatchStart(:)  = DomainStart(:)
  PatchEnd(:)    = DomainEnd(:)
  MemoryStart(:) = DomainStart(:)
  MemoryEnd(:)   = DomainEnd(:)

  IF ( myproc == 0 )  THEN

    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading dump variable ', varname

    CALL ext_int_read_field(nid, DateStr, VarName, temtd, varType,      &
                            0, 0, 1, MemoryOrder, Stagger, DimNames,    &
                            DomainStart, DomainEnd,          &
                            MemoryStart, MemoryEnd,          & ! Memory
                            PatchStart,  PatchEnd,           & ! Patch
                            iStatus)

    WRITE(6,'(a)') '       ...  DONE.'

  END IF

  RETURN
END SUBROUTINE get_wrf_dummy_from_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_wrf_1d                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_1d_from_one_file(nfid,io_form,                       &
                      datestr,itime,varname,stagger,                    &
                      dimname1,nz,var1d, nzd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 1D array from the WRF history file. Only root processor
!    reads and then broadcast to all other processors if necessary
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nfid
  INTEGER,          INTENT(IN)  :: io_form
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  INTEGER,          INTENT(IN)  :: nz              ! memory index

  REAL,             INTENT(OUT) :: var1d(nz)
  INTEGER,          INTENT(IN)  :: nzd            ! data index
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 1D variable ', varname

  DimNames(1) = dimname1
  DimNames(2:3) = ''

  DomainStart(1:3) = 1
  DomainEnd(2:3)   = 1
  DomainEnd(1)     = nzd

  MemoryStart(1:3) = 1
  MemoryEnd(2:3)   = 1
  MemoryEnd(1)     = nz

  PatchStart(1:3) = 1
  PatchEnd(2:3)   = 1
  PatchEnd(1)     = nzd

  IF (io_form == 5) THEN                  ! PHDF5 format

    CALL get_phdf5_field(nfid, DateStr, VarName, var1d, WRF_REAL,       &
                         1,'Z  ',stagger,DimNames,                      &
                         DomainStart,DomainEnd,MemoryStart,MemoryEnd,   &
                         PatchStart,PatchEnd,                           &
                         iStatus)

  ELSE            ! need to broadcast

    IF (io_form == 7) THEN                       ! NetCDF format
      IF (myproc == 0) &
        CALL get_ncd_1d(nfid,itime,varname,nzd,var1d,istatus)

    ELSE IF (io_form == 1) THEN

      IF (myproc == 0) &
        CALL ext_int_read_field(nfid, DateStr, VarName, var1d, WRF_REAL,&
                            0, 0, 1, 'Z  ', Stagger , DimNames ,        &
                            DomainStart, DomainEnd,          &
                            MemoryStart, MemoryEnd,          & ! Memory
                            PatchStart,  PatchEnd,           & ! Patch
                            iStatus)
    END IF
    CALL mpupdater(var1d,nz)

  END IF

  IF ( myproc == 0 ) THEN
    IF ( istatus == 0 ) THEN
      WRITE(6,'(a)') '         ...  DONE.'
    ELSE
      WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE get_wrf_1d_from_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_wrf_2d                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_2d_from_one_file(nfid,io_form,                       &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2d,                    &
                      nxd,nyd,temdom,nxlg,nylg,temlg,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 2D array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nfid
  INTEGER,          INTENT(IN)  :: io_form
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  INTEGER,          INTENT(IN)  :: fzone
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  CHARACTER(LEN=*), INTENT(IN)  :: dimname2
  INTEGER,          INTENT(IN)  :: nx              ! local index
  INTEGER,          INTENT(IN)  :: ny

  REAL,             INTENT(OUT) :: var2d(nx,ny)
  INTEGER,          INTENT(IN)  :: nxd,nyd         ! data index
  REAL,             INTENT(OUT) :: temdom(nxd,nyd) ! data array
  INTEGER,          INTENT(IN)  :: nxlg,nylg       ! memory index for the whole domain
  REAL,             INTENT(OUT) :: temlg(nxlg,nylg)! memory array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  INCLUDE 'mp.inc'

  INTEGER            :: i, j

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: nxdim, nydim
  INTEGER       :: ilocs,iloce,jlocs,jloce

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 2D variable ', varname

  nxdim = nxlg - 1    ! domain size
  nydim = nylg - 1
  IF (Stagger == 'X') nxdim = nxlg
  IF (Stagger == 'Y') nydim = nylg

  ilocs = (nx-fzone)*(loc_x-1)+fzone
  jlocs = (ny-fzone)*(loc_y-1)+fzone
  iloce = (nx-fzone)*(loc_x)+fzone
  jloce = (ny-fzone)*(loc_y)+fzone

  DimNames(1) = dimname1
  DimNames(2) = dimname2
  DimNames(3) = ''

  DomainStart(1) = 1
  DomainStart(2) = 1
  DomainStart(3) = 1
  DomainEnd(1)   = nxdim
  DomainEnd(2)   = nydim
  DomainEnd(3)   = 1

  IF (io_form == 5) THEN                  ! PHDF5 format

    MemoryStart(1) = ilocs
    MemoryStart(2) = jlocs
    MemoryStart(3) = 1
    MemoryEnd(1)   = iloce
    MemoryEnd(2)   = jloce
    MemoryEnd(3)   = 1

    PatchStart(1) = ilocs
    PatchEnd(1)   = iloce
    IF (stagger /= 'X' .AND. loc_x == nproc_x) THEN
      PatchEnd(1)   =  iloce - fzone
    END IF

    PatchStart(2) = jlocs
    PatchEnd(2)   = jloce
    IF (stagger /= 'Y' .AND. loc_y == nproc_y) THEN
      PatchEnd(2)   =  jloce - fzone
    END IF

    PatchStart(3) = 1
    PatchEnd(3)   = 1

    CALL get_phdf5_field(nfid, DateStr, VarName, var2d, WRF_REAL,       &
                         1,'XY',Stagger,DimNames,                       &
                         DomainStart,DomainEnd,MemoryStart,MemoryEnd,   &
                         PatchStart,PatchEnd,                           &
                         iStatus)

   IF (loc_x == nproc_x .AND. PatchEnd(1) < MemoryEnd(1) )   &
      var2d(nx,:) = var2d(nx-1,:)

   IF (loc_y == nproc_y .AND. PatchEnd(2) < MemoryEnd(2) )   &
      var2d(:,ny) = var2d(:,ny-1)

  ELSE    ! need to split

    IF (io_form == 7) THEN                       ! NetCDF format
      IF (myproc == 0) &
        CALL get_ncd_2d(nfid,itime,varname,nxd,nyd,temdom,istatus)

    ELSE IF (io_form == 1) THEN

      PatchStart(:)  = DomainStart(:)
      PatchEnd(:)    = DomainEnd(:)
      MemoryStart(:) = PatchStart(:)
      MemoryEnd(:)   = PatchEnd(:)

      IF (myproc == 0)   &
        CALL ext_int_read_field(nfid, DateStr, VarName, temdom,         &
                            WRF_REAL,0, 0, 1, 'XY', Stagger, DimNames,  &
                            DomainStart, DomainEnd,                     &
                            MemoryStart, MemoryEnd,          & ! Memory
                            PatchStart,  PatchEnd,           & ! Patch
                            iStatus)
    END IF

    IF (myproc == 0) THEN
      DO j = 1,nyd
        DO i = 1,nxd
            temlg(i,j) = temdom(i,j)
        END DO
      END DO
      CALL edgfill(temlg,1,nxlg,1,nxd,1,nylg,1,nyd,1,1,1,1)
    END IF

    CALL wrf_split2d(temlg,nx,ny,fzone,var2d)

  END IF

  IF ( myproc == 0 ) THEN
    IF (istatus == 0) THEN
       WRITE(6,'(a)') '         ...  DONE.'
    ELSE
       WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE get_wrf_2d_from_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_wrf_2di                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_2di_from_one_file(nfid,io_form,                      &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2di,                   &
                      nxd,nyd,temdm,nxlg,nylg,temlg,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 2D integer array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nfid
  INTEGER,          INTENT(IN)  :: io_form
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  INTEGER,          INTENT(IN)  :: fzone
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  CHARACTER(LEN=*), INTENT(IN)  :: dimname2
  INTEGER,          INTENT(IN)  :: nx              ! local index
  INTEGER,          INTENT(IN)  :: ny

  INTEGER,          INTENT(OUT) :: var2di(nx,ny)
  INTEGER,          INTENT(IN)  :: nxd,nyd          ! domain index
  INTEGER,          INTENT(OUT) :: temdm(nxd,nyd)   ! domain array
  INTEGER,          INTENT(IN)  :: nxlg,nylg          ! memory
  INTEGER,          INTENT(OUT) :: temlg(nxlg,nylg)   ! memory array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  INCLUDE 'mp.inc'

  INTEGER            :: i, j

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: nxdim, nydim
  INTEGER       :: ilocs,iloce,jlocs,jloce

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 2D integer variable ', varname

  nxdim = nxlg - 1    ! domain size
  nydim = nylg - 1
  IF (Stagger == 'X') nxdim = nxlg
  IF (Stagger == 'Y') nydim = nylg

  ilocs = (nx-fzone)*(loc_x-1)+fzone
  jlocs = (ny-fzone)*(loc_y-1)+fzone
  iloce = (nx-fzone)*(loc_x)+fzone
  jloce = (ny-fzone)*(loc_y)+fzone

  DimNames(1) = dimname1
  DimNames(2) = dimname2
  DimNames(3) = ''

  DomainStart(1) = 1
  DomainStart(2) = 1
  DomainStart(3) = 1
  DomainEnd(1)   = nxdim
  DomainEnd(2)   = nydim
  DomainEnd(3)   = 1

  IF (io_form == 5) THEN                  ! PHDF5 format

    MemoryStart(1) = ilocs
    MemoryStart(2) = jlocs
    MemoryStart(3) = 1
    MemoryEnd(1)   = iloce
    MemoryEnd(2)   = jloce
    MemoryEnd(3)   = 1

    PatchStart(1) = ilocs
    PatchEnd(1)   = iloce
    IF (stagger /= 'X' .AND. loc_x == nproc_x) THEN
      PatchEnd(1)   =  iloce - fzone
    END IF

    PatchStart(2) = jlocs
    PatchEnd(2)   = jloce
    IF (stagger /= 'Y' .AND. loc_y == nproc_y) THEN
      PatchEnd(2)   =  jloce - fzone
    END IF

    PatchStart(3) = 1
    PatchEnd(3)   = 1

    CALL get_phdf5_field(nfid, DateStr, VarName, var2di, WRF_INTEGER,   &
                         1,'XY',Stagger,DimNames,                       &
                         DomainStart,DomainEnd,MemoryStart,MemoryEnd,   &
                         PatchStart,PatchEnd,                           &
                         iStatus)

    IF (loc_x == nproc_x .AND. PatchEnd(1) < MemoryEnd(1) )   &
       var2di(nx,:) = var2di(nx-1,:)

    IF (loc_y == nproc_y .AND. PatchEnd(2) < MemoryEnd(2) )   &
       var2di(:,ny) = var2di(:,ny-1)

  ELSE  ! To be split

    IF (io_form == 7) THEN                       ! NetCDF format
      IF (myproc == 0) &
        CALL get_ncd_2di(nfid,itime,varname,nxd,nyd,temdm,istatus)

    ELSE IF (io_form == 1) THEN

      PatchStart(:) = DomainStart(:)
      PatchEnd(:)   = DomainEnd(:)
      MemoryStart(:) = PatchStart(:)
      MemoryEnd(:)   = PatchEnd(:)

      IF (myproc == 0)   &

        CALL ext_int_read_field(nfid, DateStr, VarName, temdm,          &
                            WRF_INTEGER,0,0,1,'XY',Stagger,DimNames ,   &
                            DomainStart, DomainEnd,                     &
                            MemoryStart, MemoryEnd,          & ! Memory
                            PatchStart,  PatchEnd,           & ! Patch
                            iStatus)
    END IF

    IF (myproc == 0) THEN
      DO j = 1,nyd
        DO i = 1,nxd
            temlg(i,j) = temdm(i,j)
        END DO
      END DO
      CALL iedgfill(temlg,1,nxlg,1,nxd,1,nylg,1,nyd,1,1,1,1)

    END IF      ! myproc == 0

    CALL wrf_split2di(temlg,nx,ny,fzone,var2di)

  END IF

  IF ( myproc == 0 ) THEN
    IF (istatus == 0) THEN
      WRITE(6,'(a)') '         ...  DONE.'
    ELSE
      WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE get_wrf_2di_from_one_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_wrf_3d                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_3d_from_one_file(nfid,io_form,                       &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,dimname3,nx,ny,nz,var3d,        &
                      nxd,nyd,nzd,temdm,nxlg,nylg,nzlg,temlg,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file. Only root processor do
!    the reads, and broadcast to all other processors if necessary.
!
!  Note: the subroutine handles the memeory order 'XZY' for 3d variables.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nfid
  INTEGER,          INTENT(IN)  :: io_form
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  INTEGER,          INTENT(IN)  :: fzone
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  CHARACTER(LEN=*), INTENT(IN)  :: dimname2
  CHARACTER(LEN=*), INTENT(IN)  :: dimname3
  INTEGER,          INTENT(IN)  :: nx              ! local index
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  REAL,             INTENT(OUT) :: var3d(nx,ny,nz)
  INTEGER,          INTENT(IN)  :: nxd,nyd,nzd             ! Data index
  REAL,             INTENT(OUT) :: temdm(nxd*nyd*nzd)      ! domain array
  INTEGER,          INTENT(IN)  :: nxlg,nylg,nzlg          ! domain index
  REAL,             INTENT(OUT) :: temlg(nxlg,nylg,nzlg)   ! memory array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  INCLUDE 'mp.inc'

  INTEGER            :: i, j, k
  INTEGER            :: i1,j1,k1

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: nxdim, nydim, nzdim
  INTEGER       :: ilocs,iloce,jlocs,jloce

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
     WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 3D variable ', varname

  nxdim = nxlg - 1    ! domain size
  nydim = nylg - 1
  nzdim = nzlg - 1
  IF (Stagger == 'X') nxdim = nxlg
  IF (Stagger == 'Y') nydim = nylg
  IF (Stagger == 'Z') nzdim = nzlg

  ilocs = (nx-fzone)*(loc_x-1)+fzone
  jlocs = (ny-fzone)*(loc_y-1)+fzone
  iloce = (nx-fzone)*(loc_x)+fzone
  jloce = (ny-fzone)*(loc_y)+fzone

  IF (io_form == 5) THEN

    DimNames(1) = dimname1
    DimNames(2) = dimname2
    DimNames(3) = dimname3

    DomainStart(1) = 1
    DomainStart(2) = 1
    DomainStart(3) = 1
    DomainEnd(1)   = nxdim
    DomainEnd(2)   = nydim
    DomainEnd(3)   = nzdim

    MemoryStart(1) = ilocs
    MemoryStart(2) = jlocs
    MemoryStart(3) = 1
    MemoryEnd(1)   = iloce
    MemoryEnd(2)   = jloce
    MemoryEnd(3)   = nz

    PatchStart(1) = ilocs
    PatchEnd(1)   = iloce
    IF (stagger /= 'X' .AND. loc_x == nproc_x) THEN
      PatchEnd(1)   =  iloce - fzone
    END IF

    PatchStart(2) = jlocs
    PatchEnd(2)   = jloce
    IF (stagger /= 'Y' .AND. loc_y == nproc_y) THEN
      PatchEnd(2)   =  jloce - fzone
    END IF

    PatchStart(3) = 1
    PatchEnd(3)   = nzd

    CALL get_phdf5_field(nfid, DateStr, VarName, var3d, WRF_REAL,       &
                         1,'XYZ',Stagger,DimNames,                      &
                         DomainStart,DomainEnd,MemoryStart,MemoryEnd,   &
                         PatchStart,PatchEnd,                           &
                         iStatus)

    !
    !  Supply data at the edge points (zero gradient, where missing)
    !
    IF (loc_x == nproc_x .AND. PatchEnd(1) < MemoryEnd(1) )   &
       var3d(nx,:,:) = var3d(nx-1,:,:)

    IF (loc_y == nproc_y .AND. PatchEnd(2) < MemoryEnd(2) )   &
       var3d(:,ny,:) = var3d(:,ny-1,:)

    IF (PatchEnd(3) < MemoryEnd(3))  var3d(:,:,nz) = var3d(:,:,nz-1)

  ELSE

    IF (io_form == 7) THEN

      IF (myproc == 0) THEN
        CALL get_ncd_3d(nfid,itime,varname,nxd,nyd,nzd,temdm,istatus)

        DO k = 1,nzd
          k1 = (k-1)*nxd*nyd
          DO j = 1,nyd
            j1 = (j-1)*nxd
            DO i = 1,nxd
              i1 = i + j1 + k1
              temlg(i,j,k) = temdm(i1)
            END DO
          END DO
        END DO
        CALL edgfill(temlg,1,nxlg,1,nxd,1,nylg,1,nyd,1,nzlg,1,nzd)
      END IF

    ELSE IF (io_form == 1) THEN

      DimNames(1) = dimname1
      DimNames(2) = dimname3
      DimNames(3) = dimname2

      DomainStart(1) = 1
      DomainStart(2) = 1
      DomainStart(3) = 1
      DomainEnd(1)   = nxdim
      DomainEnd(2)   = nzdim
      DomainEnd(3)   = nydim

      IF (myproc == 0) THEN
        CALL ext_int_read_field(nfid, DateStr, VarName, temdm, WRF_REAL,&
                              0, 0, 1, 'XZY', Stagger , DimNames ,      &
                              DomainStart, DomainEnd,                   &
                              DomainStart, DomainEnd,          & ! Memory
                              DomainStart, DomainEnd,          & ! Patch
                              iStatus)
        DO k = 1,nzd
          k1 = (k-1)*nxd
          DO j = 1,nyd
            j1 = (j-1)*nxd*nzd
            DO i = 1,nxd
              i1 = i + j1 + k1
              temlg(i,j,k) = temdm(i1)
            END DO
          END DO
        END DO
        CALL edgfill(temlg,1,nxlg,1,nxd,1,nylg,1,nyd,1,nzlg,1,nzd)
      END IF

    END IF
    CALL wrf_split3d(temlg,nx,ny,nz,fzone,var3d)

  END IF

  IF ( myproc == 0 ) THEN
    IF (istatus == 0) THEN
      WRITE(6,'(a)') '         ...  DONE.'
    ELSE
      WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE get_wrf_3d_from_one_file
