!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%                                                          %%%%
!%%%% NOTE: Before calling subroutines in this file, it is     %%%%
!%%%%       assumed that the program is in MPI multifile mode. %%%%
!%%%%       If ncompressx = ncompressy = 1, the subroutines    %%%%
!%%%%       do not need to join smaller tiles, otherwiser,     %%%%
!%%%%       join smaller WRF tiles to get a local ARPS patch.  %%%%
!%%%%                                                          %%%%
!%%%%       All subroutines in this file do not have to        %%%%
!%%%%       support PHDF5 format because PHDF5 file does       %%%%
!%%%%       not in split files.                                %%%%
!%%%%                                                          %%%%
!%%%%       All file IDs are size of (ncompressx,ncompressy)   %%%%
!%%%%       arrays, except the following attribute related     %%%%
!%%%%       subroutines:                                       %%%%
!%%%%          get_wrf_meta_from_multi_files                   %%%%
!%%%%          get_dom_ti_integer                              %%%%
!%%%%          get_dom_ti_real                                 %%%%
!%%%%          get_dom_ti_char                                 %%%%
!%%%%                                                          %%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
SUBROUTINE open_wrf_multi_files(filename,io_form,for_meta_only,         &
                                ncompressx,ncompressy,numdigits,        &
                                nidout,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Open a WRF file and return NetCDF file handler. Each processor
!    may open multifiles depends on ncompressx/ncompressy. However,
!    for_meta_only is .TRUE., only one file will be opened.
!
!    This mode only supports binary and NetCDF format so far.
!
!    NOTE: it is required to call close_wrf_multi_files explicitly to close
!          the opened file in your calling program.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: for_meta_only
  INTEGER,          INTENT(IN)  :: ncompressx, ncompressy,numdigits
  INTEGER,          INTENT(OUT) :: nidout(ncompressx,ncompressy)
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  CHARACTER(LEN=256), ALLOCATABLE :: tmpstr(:,:)
  LOGICAL                         :: fexists

  CHARACTER(LEN=80)  :: sysdepinfo
  LOGICAL, SAVE      :: initialized = .FALSE.

  INCLUDE 'mp.inc'

  INTEGER  :: loc_proc, iloc, jloc
  INTEGER  :: iloc_x, jloc_y, nprocx_in

  CHARACTER(LEN=20)  :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  iloc_x = (loc_x-1)*ncompressx    ! column of processors
  jloc_y = (loc_y-1)*ncompressy    ! rows of processors
  nprocx_in = ncompressx*nproc_x

  ALLOCATE(tmpstr(ncompressx,ncompressy), STAT = istatus)

  WRITE(fmtstr,'(a,I1,a,I1,a)') '(2a,I',numdigits,'.',numdigits,')'

  DO jloc = 1, ncompressy
    DO iloc = 1, ncompressx

      loc_proc = (jloc_y+jloc-1)*nprocx_in + iloc_x+(iloc-1)
      WRITE(tmpstr(iloc,jloc),FMT=fmtstr) filename,'_',loc_proc

      INQUIRE(FILE = TRIM(tmpstr(iloc,jloc)), EXIST = fexists)
      IF ( .NOT. fexists ) THEN
        WRITE(6,'(3a)') 'File not found: ',tmpstr(iloc,jloc),' in open_wrf_file'
        CALL arpsstop('WRF file not exist.',1)
      ENDIF

    END DO
  END DO

  sysdepinfo = 'DATASET=HISTORY'
  nidout(:,:) = -1

  IF (io_form == 7) THEN              ! not initialize needed

    IF ( for_meta_only )  THEN
      CALL open_ncd_file(tmpstr(1,1),nidout(1,1))
      IF (istatus /= 0) THEN
        WRITE(0,'(1x,2a)') 'ERROR: Opening file ',tmpstr(1,1)
        CALL arpsstop('Open WRF file error.',1)
      END IF
    ELSE
      DO jloc = 1,ncompressy
        DO iloc = 1,ncompressx
          CALL open_ncd_file(tmpstr(iloc,jloc),nidout(iloc,jloc))
          IF (istatus /= 0) THEN
            WRITE(0,'(1x,2a)') 'ERROR: Opening file ',tmpstr(iloc,jloc)
            CALL arpsstop('Open WRF file error.',1)
          END IF
        END DO
      END DO
    END IF

  ELSE IF (io_form == 1) THEN         ! initialize explicitly

    IF (.NOT. initialized) CALL ext_int_ioinit( SysDepInfo, iStatus )

    IF ( for_meta_only )  THEN
      CALL ext_int_open_for_read (tmpstr(1,1), 0, 0, SysDepInfo,        &
                                  nidout(1,1), iStatus )
      IF (istatus /= 0) THEN
        WRITE(0,'(1x,2a)') 'ERROR: Opening file ',tmpstr(1,1)
        CALL arpsstop('Open WRF file error.',1)
      END IF

    ELSE

      DO jloc = 1,ncompressy
        DO iloc = 1,ncompressx
          CALL ext_int_open_for_read(tmpstr(iloc,jloc),0,0, SysDepInfo, &
                                     nidout(iloc,jloc),iStatus )
          IF (istatus /= 0) THEN
            WRITE(0,'(1x,2a)') 'ERROR: Opening file ',tmpstr(iloc,jloc)
            CALL arpsstop('Open WRF file error.',1)
          END IF
        END DO
      END DO

    END IF

  ELSE
    WRITE(0,*) 'Unsupported IO format - ',io_form
    CAlL arpsstop('Unsupported IO format.',1)
  END IF

  initialized = .TRUE.

  DEALLOCATE(tmpstr)

  RETURN
END SUBROUTINE open_wrf_multi_files
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
SUBROUTINE close_wrf_multi_files(nch,io_form,for_meta_only,             &
                                 ncompressx,ncompressy,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!     Close the WRF file which is opened using open_wrf_multi_file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: io_form
  LOGICAL, INTENT(IN) :: for_meta_only
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: nch(ncompressx,ncompressy)
  INTEGER, INTENT(OUT):: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
!
  INTEGER :: iloc, jloc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

  IF(io_form == 7) THEN
    IF (for_meta_only) THEN
      CALL close_ncd_file(nch(1,1))
    ELSE
      DO jloc = 1, ncompressy
        DO iloc = 1, ncompressx
          CALL close_ncd_file(nch(iloc,jloc))
        END DO
      END DO
    END IF
  ELSE IF (io_form == 1) THEN
    IF (for_meta_only) THEN
      CALL ext_int_ioclose(nch(1,1),iStatus)
    ELSE
      DO jloc = 1, ncompressy
        DO iloc = 1, ncompressx
          CALL ext_int_ioclose(nch(iloc,jloc),iStatus)
        END DO
      END DO
    END IF
  END IF

  IF (istatus /= 0) THEN
    WRITE(0,'(1x,2a)') 'ERROR: closing file handler ',nch
    CALL arpsstop('Error in close_wrf_multi_files',1)
  END IF

  RETURN
END SUBROUTINE close_wrf_multi_files
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
SUBROUTINE get_wrf_Times_from_multi_files(nfid,io_form,ncompressx,      &
                                      ncompressy,itime,timestr,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read the the Date String in the WRF outputs at specified time.
!
!    Although the first file only will provide enough information, we
!    have to read through all of the files because the sequencial
!    access requirement of WRF binary file. However, only the last
!    read timestr will be returned.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ncompressx, ncompressy
  INTEGER, INTENT(IN)  :: nfid(ncompressx,ncompressy)     ! file handler
  INTEGER, INTENT(IN)  :: io_form                         ! File format
  INTEGER, INTENT(IN)  :: itime         ! Time dimension value
                                        ! this is the unlimited dimension
  CHARACTER(LEN=*), INTENT(OUT) :: timestr
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
  INTEGER :: iloc, jloc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 7) THEN    ! Read the first file only

    CALL get_ncd_next_time(nfid(1,1),itime,timestr,istatus)

  ELSE IF (io_form == 1) THEN  ! Read through all files. Actually, only
                               ! the last file return valid value.
    DO jloc = 1, ncompressy
      DO iloc = 1, ncompressx
        CALL ext_int_get_next_time(nfid(iloc,jloc),timestr,istatus)
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE get_wrf_Times_from_multi_files
!
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
SUBROUTINE get_wrf_meta_from_multi_files(nid,io_form,dim_check,         &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj,trlat1,trlat2,trlon,ctrlat,ctrlon,      &
                          dx,dy,dt,sfcphys,sfclay,mpphys,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve WRF grib information from the NetCDF file which are stored
!   as Global attributes. Only 1 file ID is enough.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nid
  INTEGER, INTENT(IN)  :: io_form
  LOGICAL, INTENT(IN)  :: dim_check
  INTEGER, INTENT(OUT) :: nx_ext, ny_ext     ! they are whole domain dimensions
  INTEGER, INTENT(OUT) :: nz_ext, nzsoil_ext
  INTEGER, INTENT(OUT) :: iproj
  REAL,    INTENT(OUT) :: trlat1
  REAL,    INTENT(OUT) :: trlat2
  REAL,    INTENT(OUT) :: trlon
  REAL,    INTENT(OUT) :: ctrlat
  REAL,    INTENT(OUT) :: ctrlon
  REAL,    INTENT(OUT) :: dx
  REAL,    INTENT(OUT) :: dy
  REAL,    INTENT(OUT) :: dt
  INTEGER, INTENT(OUT) :: sfcphys,sfclay
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

  INTEGER           :: ilocd, jlocd
  INTEGER           :: ilocs, iloce, jlocs, jloce

  INCLUDE   'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL get_dom_ti_char(nid,io_form,'TITLE',     cdump,istatus)
  CALL get_dom_ti_char(nid,io_form,'START_DATE',cdump,istatus)

  CALL get_dom_ti_integer(nid,io_form,'WEST-EAST_GRID_DIMENSION',  nx_ext,istatus)
  CALL get_dom_ti_integer(nid,io_form,'SOUTH-NORTH_GRID_DIMENSION',ny_ext,istatus)
  CALL get_dom_ti_integer(nid,io_form,'BOTTOM-TOP_GRID_DIMENSION', nz_ext,istatus)

  CALL get_dom_ti_char   (nid,io_form,'GRIDTYPE',cdump,istatus)
!  CALL get_dom_ti_integer(nid,io_form,'DYN_OPT', idump,istatus)
  CALL get_dom_ti_integer(nid,io_form,'DIFF_OPT',idump,istatus)
  CALL get_dom_ti_integer(nid,io_form,'KM_OPT',  idump,istatus)
  CALL get_dom_ti_integer(nid,io_form,'DAMP_OPT',idump,istatus)
  CALL get_dom_ti_real   (nid,io_form,'KHDIF',   rdump,istatus)
  CALL get_dom_ti_real   (nid,io_form,'KVDIF',   rdump,istatus)

  CALL get_dom_ti_integer(nid,io_form,'MP_PHYSICS',        mpphys, istatus)
  CALL get_dom_ti_integer(nid,io_form,'RA_LW_PHYSICS',     idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'RA_SW_PHYSICS',     idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'SF_SFCLAY_PHYSICS', sfclay, istatus)
  CALL get_dom_ti_integer(nid,io_form,'SF_SURFACE_PHYSICS',sfcphys,istatus)
  CALL get_dom_ti_integer(nid,io_form,'BL_PBL_PHYSICS',    idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'CU_PHYSICS',        idump,  istatus)

  CALL get_dom_ti_integer(nid,io_form,'WEST-EAST_PATCH_START_UNSTAG',  ips_u,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'WEST-EAST_PATCH_END_UNSTAG',    ipe_u,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'WEST-EAST_PATCH_START_STAG',    ips,    istatus)
  CALL get_dom_ti_integer(nid,io_form,'WEST-EAST_PATCH_END_STAG',      ipe,    istatus)
  CALL get_dom_ti_integer(nid,io_form,'SOUTH-NORTH_PATCH_START_UNSTAG',jps_u,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'SOUTH-NORTH_PATCH_END_UNSTAG',  jpe_u,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'SOUTH-NORTH_PATCH_START_STAG',  jps,    istatus)
  CALL get_dom_ti_integer(nid,io_form,'SOUTH-NORTH_PATCH_END_STAG',    jpe,    istatus)
  CALL get_dom_ti_integer(nid,io_form,'BOTTOM-TOP_PATCH_START_UNSTAG', idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'BOTTOM-TOP_PATCH_END_UNSTAG',   idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'BOTTOM-TOP_PATCH_START_STAG',   idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'BOTTOM-TOP_PATCH_END_STAG',     idump,  istatus)

  CALL get_dom_ti_real   (nid,io_form,'DX',dx,    istatus)
  CALL get_dom_ti_real   (nid,io_form,'DY',dy,    istatus)
  CALL get_dom_ti_real   (nid,io_form,'DT',dt,    istatus)

  CALL get_dom_ti_real   (nid,io_form,'CEN_LAT',      ctrlat, istatus)
  CALL get_dom_ti_real   (nid,io_form,'CEN_LON',      ctrlon, istatus)
  CALL get_dom_ti_real   (nid,io_form,'TRUELAT1',     trlat1, istatus)
  CALL get_dom_ti_real   (nid,io_form,'TRUELAT2',     trlat2, istatus)

  CALL get_dom_ti_real   (nid,io_form,'MOAD_CEN_LAT', rdump,  istatus)
  CALL get_dom_ti_real   (nid,io_form,'STAND_LON',    trlon,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'MAP_PROJ',     iproj,  istatus)

  CALL get_dom_ti_char   (nid,io_form,'MMINLU',       cdump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'ISWATER',      idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'ISICE',        idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'ISURBAN',      idump,  istatus)
  CALL get_dom_ti_integer(nid,io_form,'ISOILWATER',   idump,  istatus)

!-----------------------------------------------------------------------
!
!  Do some dimension checks, We know mp_opt >0 and it is in multifile mode
!
!-----------------------------------------------------------------------

  IF (dim_check) THEN

    ilocd = (nx_ext-1)/nproc_x
    jlocd = (ny_ext-1)/nproc_y

    ilocs = ilocd*(loc_x-1)+1       ! Patch start for both stag and unstag
    jlocs = jlocd*(loc_y-1)+1
!    iloce = ilocd*loc_x            ! Patch end for unstag
!    jloce = jlocd*loc_y            ! not sure because of readjoin == 1

    IF (ilocs /= ips_u .OR. ilocs /= ips) THEN  ! We are only sure about
                                                ! for the patch start index
      WRITE(0,'(1x,2a,I4.4,a,/,3(7x,a,I4,a,/))')                        &
              'ERROR: Patch size in X direction is not correct in ',    &
              'processor ',myproc,'.',                                  &
              ' Expecting patch started from   ',ilocs, ',',            &
              ' Found unstag points started at ',ips_u, ',',            &
              '     and stag points started at ',ips,   '.'
      CALL arpsstop('Wrong size in split files.',1)
    END IF

    IF (jlocs /= jps_u .OR. jlocs /= jps) THEN
      WRITE(0,'(1x,2a,I4.4,a,/,3(7x,a,I4,a,/))')                        &
              'ERROR: Patch size in Y direction is not correct in ',    &
              ' processor ', myproc,'.',                                &
              ' Expecting patch started from   ',jlocs, ',',            &
              ' Found unstag points started at ',jps_u, ',',            &
              '     and stag points started at ',jps,   '.'
      CALL arpsstop('Wrong size in split files.',1)
    END IF

  END IF

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
END SUBROUTINE get_wrf_meta_from_multi_files

SUBROUTINE get_dom_ti_integer(nid,io_form,element, val, ireturn)
!
!-----------------------------------------------------------------------
!
! NOTE: do not have to broadcast because we are sure it is in
!       MPI multifile mode.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: element
  INTEGER,      INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INTEGER    :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 7) THEN
    CALL get_ncd_dom_ti_integer(nid,element,val,ireturn)
  ELSE IF (io_form == 1) THEN
    CALL ext_int_get_dom_ti_integer(nid,element,val,1,outcount,ireturn)
  END IF

  RETURN
END SUBROUTINE get_dom_ti_integer

SUBROUTINE  get_dom_ti_real(nid,io_form,element, val, ireturn)
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: element
  REAL,         INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INTEGER   :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 7) THEN
    CALL get_ncd_dom_ti_real(nid,element,val,ireturn)
  ELSE IF (io_form == 1) THEN
    CALL ext_int_get_dom_ti_real(nid,element,val,1, outcount,ireturn)
  END IF

  RETURN
END SUBROUTINE get_dom_ti_real
!
SUBROUTINE  get_dom_ti_char(nid,io_form,element, val, ireturn)
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: element
  CHARACTER(*), INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INTEGER   :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (io_form == 7) THEN
    CALL get_ncd_dom_ti_char(nid,element,val,ireturn)
  ELSE IF (io_form == 1) THEN
    CALL ext_int_get_dom_ti_char(nid,element,val,ireturn)
  END IF

  RETURN
END SUBROUTINE get_dom_ti_char
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
SUBROUTINE get_wrf_dummy_from_multi_files(nid,io_form,                  &
                ncompressx,ncompressy,datestr,itime,                    &
                varname,varType,memoryorder,stagger,                    &
                dimname1,dimname2,dimname3,nx,ny,nz,nxd,nyd,nzd,temtd,  &
                istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in an array from the WRF history file. It is just for sequential
!    access of WRF binary file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(IN)  :: ncompressx,ncompressy
  INTEGER,          INTENT(IN)  :: nid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: varType
  CHARACTER(LEN=*), INTENT(IN)  :: MemoryOrder
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  CHARACTER(LEN=*), INTENT(IN)  :: dimname2
  CHARACTER(LEN=*), INTENT(IN)  :: dimname3
  INTEGER,          INTENT(IN)  :: nx              ! ARPS patch size
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  INTEGER,          INTENT(IN)  :: nxd,nyd,nzd     ! WRF data size
  REAL,             INTENT(OUT) :: temtd(nxd*nyd*nzd)  ! data in files
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: xdim, ydim, zdim
  INTEGER       :: nxlg, nylg
  INTEGER       :: ilocs,iloce,jlocs,jloce

  LOGICAL       :: patched
  INTEGER       :: iloc, jloc
  INTEGER       :: nxt, nyt     ! WRF title size
                                ! it should be nxd,nyd except at boundary

  INCLUDE   'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF ( io_form /= 1 ) RETURN

  nxlg = (nx-1)*nproc_x        ! domain size
  nylg = (ny-1)*nproc_y
  IF (stagger == 'X') nxlg = nxlg + 1
  IF (stagger == 'Y') nylg = nylg + 1

  IF (MemoryOrder == 'XZY') THEN
    xdim = 1
    ydim = 3
    zdim = 2
    Patched = .TRUE.
  ELSE IF (MemoryOrder == 'XY' .OR. MemoryOrder == 'XYZ') THEN
    xdim = 1
    ydim = 2
    zdim = 3
    Patched = .TRUE.
  ELSE
    xdim = 1
    ydim = 2
    zdim = 3

    nxlg = nxd   ! 1D arrays, used pass-in instead of domain index
    nylg = nyd

    Patched = .False.
  END IF

  DimNames(xdim) = dimname1
  DimNames(ydim) = dimname2
  DimNames(zdim) = dimname3

  DomainStart(:)  = 1
  DomainEnd(xdim) = nxlg
  DomainEnd(ydim) = nylg
  DomainEnd(zdim) = nzd


  IF(myproc == 0)    &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading dump variable ', varname

  PatchStart(:)  = DomainStart(:)
  PatchEnd(:)    = DomainEnd(:)
  MemoryStart(:) = DomainStart(:)
  MemoryEnd(:)   = DomainEnd(:)

  ilocs = (nx-1)*(loc_x-1)+1      ! ARPS patch start
  jlocs = (ny-1)*(loc_y-1)+1

  IF (patched) THEN

  nxt = (nx-1)/ncompressx         ! WRF tile size
  nyt = (ny-1)/ncompressy
!-----------------------------------------------------------------------
! Do some check, should be commented out before releasing.
!-----------------------------------------------------------------------

  IF (nxt /= nxd .AND. nxt /= nxd-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nxd = ',nxd,', nxt = ',nxt
  END IF
  IF (nyt /= nyd .AND. nyt /= nyd-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nyd = ',nyd,', nyt = ',nyt
  END IF

  END IF   ! patched


  DO jloc = 1, ncompressy
    DO iloc = 1, ncompressx

      IF ( patched ) THEN   ! 2D or 3D arrays
        PatchStart(xdim) = ilocs + nxt*(iloc-1)  ! WRF patch
        PatchStart(ydim) = jlocs + nyt*(jloc-1)
        PatchEnd(xdim)   = ilocs + nxt*iloc - 1
        PatchEnd(ydim)   = jlocs + nyt*jloc - 1
        IF (stagger == 'X' .AND. loc_x == nproc_x .AND.                 &
            iloc == ncompressx) THEN       ! last staggered WRF patch
          PatchEnd(xdim)   = PatchEnd(xdim)+1
        END IF
        IF (stagger == 'Y' .AND. loc_y == nproc_y .AND.                 &
            jloc == ncompressy) THEN       ! last staggered WRF patch
          PatchEnd(ydim)   = PatchEnd(ydim)+1
        END IF

        PatchStart(zdim) = 1
        PatchEnd(zdim) = nzd

        MemoryStart(:) = PatchStart(:)
        MemoryEnd(:) = PatchEnd(:)
      END IF

      CALL ext_int_read_field(nid(iloc,jloc), DateStr, VarName, temtd,  &
                     varType, 0, 0, 1, MemoryOrder, Stagger, DimNames,  &
                     DomainStart, DomainEnd,          &
                     MemoryStart, MemoryEnd,          & ! Memory
                     PatchStart,  PatchEnd,           & ! Patch
                     iStatus)
    END DO
  END DO

  IF(myproc == 0) WRITE(6,'(a)') '       ...  DONE.'

  RETURN
END SUBROUTINE get_wrf_dummy_from_multi_files
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
SUBROUTINE get_wrf_1d_from_multi_files(nfid,io_form,ncompressx,ncompressy, &
                      datestr,itime,varname,stagger,                    &
                      dimname1,nz,var1d,nzd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 1D array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(IN)  :: ncompressx, ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
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

  INTEGER       :: iloc, jloc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 1D variable ', varname

  IF (io_form == 7) THEN                       ! NetCDF format

    CALL get_ncd_1d(nfid(1,1),itime,varname,nzd,var1d,istatus)

  ELSE IF (io_form == 1) THEN

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

    DO jloc = 1, ncompressy
      DO iloc = 1, ncompressx

        CALL ext_int_read_field(nfid(iloc,jloc), DateStr, VarName,      &
                     var1d, WRF_REAL,0, 0, 1, 'Z  ', Stagger, DimNames, &
                     DomainStart, DomainEnd,          &
                     MemoryStart, MemoryEnd,          & ! Memory
                     PatchStart,  PatchEnd,           & ! Patch
                     iStatus)
      END DO
    END DO

  ELSE
    WRITE(0,*) 'ERROR: unsupported io_form = ',io_form
  END IF

  IF ( myproc == 0 )  WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE get_wrf_1d_from_multi_files
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
SUBROUTINE get_wrf_2d_from_multi_files(nfid,io_form,ncompressx,ncompressy, &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2d,                    &
                      nxdin,nydin,temtd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 2D array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(IN)  :: ncompressx, ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN)  :: datestr
  INTEGER,          INTENT(IN)  :: itime
  INTEGER,          INTENT(IN)  :: fzone
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: stagger
  CHARACTER(LEN=*), INTENT(IN)  :: dimname1
  CHARACTER(LEN=*), INTENT(IN)  :: dimname2
  INTEGER,          INTENT(IN)  :: nx              ! ARPS patch size
  INTEGER,          INTENT(IN)  :: ny

  REAL,             INTENT(OUT) :: var2d(nx,ny)
  INTEGER,          INTENT(IN)  :: nxdin,nydin     ! Data patch size
  REAL,             INTENT(OUT) :: temtd(nxdin*nydin)  ! data array
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

  INTEGER       :: ilocs,iloce,jlocs,jloce

  INTEGER       :: iloc, jloc     ! Processor index

!-----------------------------------------------------------------------
!
! Dimensions:
!
!   nxdim, nydim  --  Whole domain size, stagger or unstagger
!   nx,    ny     --  ARPS patch stagger size, (IN)
!   nxp,   nyp    --  ARPS patch size contains acutal data from WRF files
!   nxdin, nydin  --  WRF data patch size, stagger (IN)
!   nxt,   nyt    --  WRF data patch size, unstagger
!   nxd,   nyd    --  WRF actual data size in files, stagger or unstagger
!
!-----------------------------------------------------------------------

  INTEGER       :: nxdim, nydim
  INTEGER       :: nxt, nyt       ! unstagger WRF patch size
  INTEGER       :: nxd, nyd       ! WRF data patch size
                                  ! maybe stagger or unstagger
  INTEGER       :: nxp, nyp

  INTEGER       :: ia, ja, kin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 2D variable ', varname

  nxdim = (nx-1)*nproc_x               ! domain size
  nydim = (ny-1)*nproc_y
  IF (Stagger == 'X') nxdim = nxdim + 1
  IF (Stagger == 'Y') nydim = nydim + 1

  nxp   = nx-1                        ! ARPS patch size
  nyp   = ny-1
  IF (Stagger == 'X' .AND. loc_x == nproc_x) nxp = nxp+1
  IF (Stagger == 'Y' .AND. loc_y == nproc_y) nyp = nyp+1

  nxt = (nx-1)/ncompressx              ! WRF patch size
  nyt = (ny-1)/ncompressy

!-----------------------------------------------------------------------
! Do some check, should be commented out before releasing.
!-----------------------------------------------------------------------

  IF (nxt /= nxdin .AND. nxt /= nxdin-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nxdin = ',nxdin,', nxt = ',nxt
  END IF
  IF (nyt /= nydin .AND. nyt /= nydin-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nydin = ',nydin,', nyt = ',nyt
  END IF

  DimNames(1) = dimname1
  DimNames(2) = dimname2
  DimNames(3) = ''

  DomainStart(1) = 1
  DomainStart(2) = 1
  DomainStart(3) = 1
  DomainEnd(1)   = nxdim
  DomainEnd(2)   = nydim
  DomainEnd(3)   = 1

  PatchStart(3) = 1             ! because it is 2D data
  PatchEnd(3)   = 1

  ilocs = (nx-fzone)*(loc_x-1)+fzone   ! ARPS patch starts
  jlocs = (ny-fzone)*(loc_y-1)+fzone

  DO jloc = 1, ncompressy
    DO iloc = 1, ncompressx

      nxd = nxt
      nyd = nyt
      IF (stagger == 'X' .AND. loc_x == nproc_x .AND.                 &
          iloc == ncompressx)  nxd = nxt+1
      IF (stagger == 'Y' .AND. loc_y == nproc_y .AND.                 &
          jloc == ncompressy)  nyd = nyt+1

      PatchStart(1) = ilocs + (iloc-1)*nxt
      PatchStart(2) = jlocs + (jloc-1)*nyt
      PatchEnd(1)   = PatchStart(1) + nxd - 1
      PatchEnd(2)   = PatchStart(2) + nyd - 1

      MemoryStart(:) = PatchStart(:)
      MemoryEnd(:)   = PatchEnd(:)

      IF (io_form == 7) THEN                       ! NetCDF format

        CALL get_ncd_2d(nfid(iloc,jloc),itime,varname,nxd,nyd,temtd,istatus)

      ELSE IF (io_form == 1) THEN

        CALL ext_int_read_field(nfid(iloc,jloc),DateStr,VarName,temtd,  &
                            WRF_REAL,0, 0, 1, 'XY', Stagger, DimNames,  &
                            DomainStart, DomainEnd,          &
                            MemoryStart, MemoryEnd,          & ! Memory
                            PatchStart,  PatchEnd,           & ! Patch
                            iStatus)
      ELSE
        WRITE(0,*) 'ERROR: unsupported io_form = ',io_form
      END IF

      DO j = 1,nyd          ! join WRF patches into ARPS patch
        DO i = 1,nxd
          kin = i + (j-1)*nxd
          ia  = PatchStart(1) - ilocs + i
          ja  = PatchStart(2) - jlocs + j
          var2d(ia,ja) = temtd(kin)
        END DO
      END DO

    END DO
  END DO
  CALL edgfill(var2d,1,nx,1,nxp,1,ny,1,nyp,1,1,1,1)

  IF ( myproc == 0 ) WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE get_wrf_2d_from_multi_files
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
SUBROUTINE get_wrf_2di_from_multi_files(nfid,io_form,ncompressx,ncompressy, &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2di,nxdin,nydin,temtd, &
                      istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 2D integer array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(IN)  :: ncompressx, ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
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
  INTEGER,          INTENT(IN)  :: nxdin,nydin          ! WRF tiles size,
                                        ! large enough for stagger arrays
  INTEGER,          INTENT(OUT) :: temtd(nxdin*nydin)   ! memory array
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

  INTEGER       :: ilocs,iloce,jlocs,jloce

!-----------------------------------------------------------------------
!
! Dimensions:
!
!   nxdim, nydim  --  Whole domain size, stagger or unstagger
!   nx,    ny     --  ARPS patch stagger size, (IN)
!   nxp,   nyp    --  ARPS patch size contains acutal data from WRF files
!   nxdin, nydin  --  WRF data patch size, stagger (IN)
!   nxt,   nyt    --  WRF data patch size, unstagger
!   nxd,   nyd    --  WRF actual data size in files, stagger or unstagger
!
!-----------------------------------------------------------------------

  INTEGER       :: nxdim, nydim
  INTEGER       :: nxt, nyt       ! unstagger WRF patch size
  INTEGER       :: nxd, nyd       ! WRF data patch size
                                  ! maybe stagger or unstagger
  INTEGER       :: nxp, nyp

  INTEGER       :: iloc, jloc
  INTEGER       :: ia, ja, kin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 2D integer variable ', varname

  nxdim = (nx - 1)*nproc_x    ! domain size
  nydim = (ny - 1)*nproc_y
  IF (Stagger == 'X') nxdim = nxdim + 1
  IF (Stagger == 'Y') nydim = nydim + 1

  nxp   = nx-1                        ! ARPS patch size
  nyp   = ny-1
  IF (Stagger == 'X' .AND. loc_x == nproc_x) nxp = nxp+1
  IF (Stagger == 'Y' .AND. loc_y == nproc_y) nyp = nyp+1

  nxt = (nx-1)/ncompressx             ! WRF patch size
  nyt = (ny-1)/ncompressy

!-----------------------------------------------------------------------
! Do some check, should be commented out before releasing.
!-----------------------------------------------------------------------

  IF (nxt /= nxdin .AND. nxt /= nxdin-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nxdin = ',nxdin,', nxt = ',nxt
  END IF
  IF (nyt /= nydin .AND. nyt /= nydin-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nydin = ',nydin,', nyt = ',nyt
  END IF

  DimNames(1) = dimname1
  DimNames(2) = dimname2
  DimNames(3) = ''

  DomainStart(1) = 1
  DomainStart(2) = 1
  DomainStart(3) = 1
  DomainEnd(1)   = nxdim
  DomainEnd(2)   = nydim
  DomainEnd(3)   = 1

  PatchStart(3) = 1             ! because it is 2D data
  PatchEnd(3)   = 1

  ilocs = (nx-fzone)*(loc_x-1)+fzone   ! ARPS patch starts
  jlocs = (ny-fzone)*(loc_y-1)+fzone

  DO jloc = 1, ncompressy
    DO iloc = 1, ncompressx

      nxd = nxt
      nyd = nyt

      IF (stagger == 'X' .AND. loc_x == nproc_x .AND.                 &
          iloc == ncompressx)  nxd = nxt+1
      IF (stagger == 'Y' .AND. loc_y == nproc_y .AND.                 &
          jloc == ncompressy)  nyd = nyt+1

      PatchStart(1) = ilocs + (iloc-1)*nxt
      PatchStart(2) = jlocs + (jloc-1)*nyt
      PatchEnd(1)   = PatchStart(1) + nxd - 1
      PatchEnd(2)   = PatchStart(2) + nyd - 1

      MemoryStart(:) = PatchStart(:)
      MemoryEnd(:)   = PatchEnd(:)

      IF (io_form == 7) THEN                       ! NetCDF format
        CALL get_ncd_2di(nfid(iloc,jloc),itime,varname,nxd,nyd,temtd,istatus)

      ELSE IF (io_form == 1) THEN

        CALL ext_int_read_field(nfid(iloc,jloc), DateStr, VarName,      &
                       temtd,WRF_INTEGER,0, 0, 1,'XY',Stagger,DimNames, &
                       DomainStart, DomainEnd,          &
                       MemoryStart, MemoryEnd,          & ! Memory
                       PatchStart,  PatchEnd,           & ! Patch
                       iStatus)
      ELSE
        WRITE(0,*) 'ERROR: unsupported io_form = ',io_form
      END IF

      DO j = 1,nyd
        DO i = 1,nxd
          kin = i + (j-1)*nxd
          ia  = PatchStart(1) - ilocs + i
          ja  = PatchStart(2) - jlocs + j
          var2di(ia,ja) = temtd(kin)
        END DO
      END DO

    END DO
  END DO
  CALL iedgfill(var2di,1,nx,1,nxp,1,ny,1,nyp,1,1,1,1)

  IF ( myproc == 0 ) WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE get_wrf_2di_from_multi_files
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
SUBROUTINE get_wrf_3d_from_multi_files(nfid,io_form,ncompressx,ncompressy, &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,dimname3,nx,ny,nz,var3d,        &
                      nxdin,nydin,nzd,temtd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(IN)  :: ncompressx,ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
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
  INTEGER,          INTENT(IN)  :: nxdin,nydin,nzd             ! Data index
  REAL,             INTENT(OUT) :: temtd(nxdin*nydin*nzd)      ! domain array
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

  INTEGER       :: ilocs,iloce,jlocs,jloce

!-----------------------------------------------------------------------
!
! Dimensions:
!
!   nxdim, nydim  --  Whole domain size, stagger or unstagger
!   nx,    ny     --  ARPS patch stagger size, (IN)
!   nxp,   nyp    --  ARPS patch size contains acutal data from WRF files
!   nxdin, nydin  --  WRF data patch size, stagger (IN)
!   nxt,   nyt    --  WRF data patch size, unstagger
!   nxd,   nyd    --  WRF actual data size in files, stagger or unstagger
!
!-----------------------------------------------------------------------

  INTEGER       :: nxdim, nydim
  INTEGER       :: nxt, nyt       ! unstagger WRF patch size
  INTEGER       :: nxd, nyd       ! WRF data patch size
                                  ! maybe stagger or unstagger
  INTEGER       :: nxp, nyp

  INTEGER       :: iloc, jloc
  INTEGER       :: ia, ja, kin
  INTEGER       :: ilocs_t, jlocs_t

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( myproc == 0 )   &
     WRITE(6,FMT='(2a)',ADVANCE='NO') '  Reading 3D variable ', varname

  nxdim = (nx - 1)*nproc_x    ! domain size
  nydim = (ny - 1)*nproc_y
  IF (Stagger == 'X') nxdim = nxdim + 1
  IF (Stagger == 'Y') nydim = nydim + 1

  nxp   = nx-1                        ! ARPS patch size
  nyp   = ny-1
  IF (Stagger == 'X' .AND. loc_x == nproc_x) nxp = nxp+1
  IF (Stagger == 'Y' .AND. loc_y == nproc_y) nyp = nyp+1

  nxt = (nx-1)/ncompressx             ! WRF patch size
  nyt = (ny-1)/ncompressy

!-----------------------------------------------------------------------
! Do some check, should be commented out before releasing.
!-----------------------------------------------------------------------

  IF (nxt /= nxdin .AND. nxt /= nxdin-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nxdin = ',nxdin,', nxt = ',nxt
  END IF
  IF (nyt /= nydin .AND. nyt /= nydin-1) THEN
    WRITE(6,*) 'Wrong in WRF patch size, nydin = ',nydin,', nyt = ',nyt
  END IF

  ilocs = (nx-fzone)*(loc_x-1)+fzone
  jlocs = (ny-fzone)*(loc_y-1)+fzone

  DO jloc = 1,ncompressy
    DO iloc = 1,ncompressx

      nxd = nxt
      nyd = nyt
      IF (stagger == 'X' .AND. loc_x == nproc_x .AND.                 &
          iloc == ncompressx)  nxd = nxt+1
      IF (stagger == 'Y' .AND. loc_y == nproc_y .AND.                 &
          jloc == ncompressy)  nyd = nyt+1

      ilocs_t = (iloc-1)*nxt   ! WRF tile start index within ARPS patch
      jlocs_t = (jloc-1)*nyt

      IF (io_form == 7) THEN

        CALL get_ncd_3d(nfid(iloc,jloc),itime,varname,nxd,nyd,nzd,temtd,istatus)

        DO k = 1,nzd
          k1 = (k-1)*nxd*nyd
          DO j = 1,nyd
            j1 = (j-1)*nxd
            DO i = 1,nxd
              kin = i + j1 + k1
              ia  = ilocs_t + i
              ja  = jlocs_t + j
              var3d(ia,ja,k) = temtd(kin)
            END DO
          END DO
        END DO

      ELSE IF (io_form == 1) THEN

        DimNames(1) = dimname1
        DimNames(2) = dimname3
        DimNames(3) = dimname2

        DomainStart(1) = 1
        DomainStart(2) = 1
        DomainStart(3) = 1
        DomainEnd(1)   = nxdim
        DomainEnd(2)   = nzd
        DomainEnd(3)   = nydim

        PatchStart(1) = ilocs + ilocs_t
        PatchEnd(1)   = PatchStart(1) + nxd - 1
        PatchStart(3) = jlocs + jlocs_t
        PatchEnd(3)   = PatchStart(3) + nyd - 1
        PatchStart(2) = 1
        PatchEnd(2)   = nzd

        MemoryStart(:) = PatchStart(:)
        MemoryEnd(:)   = PatchEnd(:)

        CALL ext_int_read_field(nfid(iloc,jloc), DateStr, VarName,      &
                     temtd, WRF_REAL,0, 0, 1, 'XZY', Stagger, DimNames, &
                     DomainStart, DomainEnd,          &
                     MemoryStart, MemoryEnd,          & ! Memory
                     PatchStart,  PatchEnd,           & ! Patch
                     iStatus)

        DO k = 1,nzd
          k1 = (k-1)*nxd
          DO j = 1,nyd
            j1 = (j-1)*nxd*nzd
            DO i = 1,nxd
              kin = i + j1 + k1
              ia  = ilocs_t + i
              ja  = jlocs_t + j
              var3d(ia,ja,k) = temtd(kin)
            END DO
          END DO
        END DO

      ELSE
        WRITE(0,*) 'ERROR: unsupported io_form = ',io_form
      END IF

    END DO     ! iloc
  END DO       ! jloc

  CALL edgfill(var3d,1,nx,1,nxp,1,ny,1,nyp,1,nz,1,nzd)

  IF ( myproc == 0 )    WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE get_wrf_3d_from_multi_files
