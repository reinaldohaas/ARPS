!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE open_output_file         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE open_output_file(filename,filetype,wrfversion,io_form,       &
                    nx,ny,nz,nzsoil,bdywidth,nout)
!------------------------------------------------------------------
!
!  PURPOSE:
!
!------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  CHARACTER(LEN=*), INTENT(IN)  :: filetype
  REAL,             INTENT(IN)  :: wrfversion
  INTEGER,          INTENT(INOUT):: io_form
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  INTEGER,          INTENT(IN)  :: nzsoil
  INTEGER,          INTENT(IN)  :: bdywidth
  INTEGER,          INTENT(OUT) :: nout

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INCLUDE 'mp.inc'

  INTEGER :: ifile
  INTEGER :: istatus

  CHARACTER(LEN=80) :: sysdepinfo
  LOGICAL           :: initialized
  LOGICAL           :: LargeFile
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (filetype == 'INPUT') THEN
    sysdepinfo = 'DATASET=INPUT'
    initialized = .FALSE.
    ifile = 1
  ELSE IF (filetype == 'STATIC') THEN
    sysdepinfo = 'DATASET=STATIC'     ! check WRF?
    initialized = .FALSE.
    ifile = 3
  ELSE
    sysdepinfo = 'DATASET=BOUNDARY'
    initialized = .TRUE.
    ifile = 2
  END IF

  LargeFile = .TRUE.
  IF (io_form > 100) LargeFile = .FALSE.

  io_form = MOD(io_form,10)

  IF (io_form == 5) THEN              ! PHDF5 format
                              ! Initialize inside open_phdf5_for_write
    CALL open_phdf5_for_write(filename,sysdepinfo,nout,initialized,istatus)

  ELSE IF (io_form == 7) THEN         ! NetCDF format
                              ! CAPS own I/O do not need initialization
    IF(myproc == 0) CALL open_ncd_for_write(filename,LargeFile,ifile,   &
                       wrfversion,nx,ny,nz,nzsoil,bdywidth,nout,istatus)

    CALL mpupdatei(nout,1)

  ELSE IF (io_form == 1) THEN         ! WRF internal format, binary

    IF(.NOT. initialized) CALL ext_int_ioinit( SysDepInfo, istatus )

    IF (myproc == 0) CALL ext_int_open_for_write( FileName,0, 0,        &
                                            SysDepInfo, nout , iStatus )

    CALL mpupdatei(nout,1)

  ELSE
    WRITE(0,*) 'ERROR: unsupport I/O format ', io_form
    istatus = -1
    nout    = -1
  END IF

  RETURN
END SUBROUTINE open_output_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE close_output_file             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE close_output_file(nch,io_form)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!   Close the output file.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nch
  INTEGER,          INTENT(IN)  :: io_form

  INCLUDE 'mp.inc'

  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (io_form == 1) THEN       ! Call WRF I/O API directly
    IF(myproc == 0) CALL ext_int_ioclose( nch, iStatus )
  ELSE IF (io_form == 7) THEN  ! use CAPS own IO API
    IF(myproc == 0) CALL close_ncd_for_write(nch,istatus)
  ELSE IF (io_form == 5) THEN  ! Call a wraper
    CALL close_phdf5_for_write(nch,istatus)
  END IF

  RETURN
END SUBROUTINE close_output_file

SUBROUTINE io_shutdown(io_form)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: io_form
  
  INTEGER :: istatus

  IF (io_form == 1) THEN
    CALL ext_int_ioexit(istatus)
  ELSE IF (io_form == 5) THEN
    CALL shutdown_phdf5_io(istatus)
  END IF

  RETURN
END SUBROUTINE io_shutdown
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE set_global_meta                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE set_global_meta(nx,ny,nz,execname,times_str,                 &
             grid_id,parent_id,i_parent_start,j_parent_start,grid_ratio,&
             dx,dy,dt,dyn_opt,diff_opt,km_opt,damp_opt,khdif,kvdif,     &
             mp_physics,ra_lw_physics, ra_sw_physics,                   &
             sf_sfclay_physics,sf_surface_physics,                      &
             bl_pbl_physics,cu_physics,                                 &
             ctrlat,ctrlon,trulat1,trulat2,trulon,                      &
             year,jday,hour,minute,second,mapproj,                      &
             dampcoef,wrfversion,global_meta)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Initialize WRF global meta data and write as global attributes
!    to the NetCDF file.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: dx,dy,dt
  CHARACTER(*), INTENT(IN) :: execname
  CHARACTER(*), INTENT(IN) :: times_str
  INTEGER, INTENT(IN)  :: grid_id, parent_id
  INTEGER, INTENT(IN)  :: i_parent_start, j_parent_start
  INTEGER, INTENT(IN)  :: grid_ratio
  INTEGER, INTENT(IN)  :: dyn_opt, diff_opt,km_opt, damp_opt
  REAL,    INTENT(IN)  :: khdif,kvdif
  INTEGER, INTENT(IN)  :: mp_physics,ra_lw_physics,ra_sw_physics
  INTEGER, INTENT(IN)  :: sf_sfclay_physics,sf_surface_physics
  INTEGER, INTENT(IN)  :: bl_pbl_physics,cu_physics

  REAL,    INTENT(IN)  :: trulat1,trulat2,trulon,ctrlat,ctrlon
  INTEGER, INTENT(IN)  :: year,jday,hour,minute,second
  INTEGER, INTENT(IN)  :: mapproj      ! ARPS flag
!                  = 1, polar projection;
!                  = 2, Lambert projection;
!                  = 3, Mercator projection.

  REAL,    INTENT(IN)  :: dampcoef

  REAL,    INTENT(IN)  :: wrfversion

  TYPE(wrf_global_metadata), INTENT(OUT) :: global_meta

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  INTEGER :: map_proj      ! WRF flag
     ! 1  -- LAMBERT CONFORMAL
     ! 2  -- POLAR STEREOGRAPHIC
     ! 3  -- MERCATOR
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF(ABS(mapproj) == 1) THEN
    map_proj = 2
  ELSE IF(ABS(mapproj) == 2) THEN
    map_proj = 1
  ELSE IF(ABS(mapproj) == 3) THEN
    map_proj = 3
  ELSE
    WRITE(6,*) 'Unknown map projection, ', mapproj
    STOP
  END IF

! set global_meta
  
  WRITE(global_meta%title,'(a,F3.1)') 'Output from '//TRIM(execname)//' for WRFV',wrfversion
  global_meta%bt_dimension   = nz
  global_meta%cen_lon        = ctrlon
  global_meta%start_date     = times_str
  global_meta%we_dimension   = nx
  global_meta%sn_dimension   = ny
  IF (wrfversion < 3.0) global_meta%dyn_opt        = dyn_opt
  global_meta%diff_opt       = diff_opt
  global_meta%km_opt         = km_opt
  global_meta%damp_opt       = damp_opt
  global_meta%dampcoef       = dampcoef
  global_meta%khdif          = khdif
  global_meta%kvdif          = kvdif
  global_meta%mp_physics     = mp_physics
  global_meta%ra_lw_physics  = ra_lw_physics
  global_meta%ra_sw_physics  = ra_sw_physics
  global_meta%sf_sfclay_physics   = sf_sfclay_physics  ! sf_sfclay_physics  for WRF 2.0
  global_meta%sf_surface_physics  = sf_surface_physics ! sf_surface_physics for WRF 2.0
  global_meta%bl_pbl_physics      = bl_pbl_physics
  global_meta%cu_physics      = cu_physics
  global_meta%we_p_unstag_s   = 1
  global_meta%we_p_unstag_e   = nx-1
  global_meta%we_p_stag_s     = 1
  global_meta%we_p_stag_e     = nx
  global_meta%sn_p_unstag_s   = 1
  global_meta%sn_p_unstag_e   = ny-1
  global_meta%sn_p_stag_s     = 1
  global_meta%sn_p_stag_e     = ny
  global_meta%bt_p_unstag_s   = 1
  global_meta%bt_p_unstag_e   = nz-1
  global_meta%bt_p_stag_s     = 1
  global_meta%bt_p_stag_e     = nz

  global_meta%dx           = dx
  global_meta%dy           = dy
  global_meta%dt           = dt 
  global_meta%cen_lat      = ctrlat
  global_meta%stand_lon    = trulon
  global_meta%moad_cen_lat = ctrlat
  global_meta%tru_lat1     = trulat1
  global_meta%tru_lat2     = trulat2

  global_meta%gmt          = hour + minute/60. + second/3600.
  global_meta%julyr        = year
  global_meta%julday       = jday
  
  global_meta%mminlu   = 'USGS'   ! 'UMD' What is it???

  global_meta%surface_input_source = 1  ! Added since WRFV2.2
  global_meta%sst_update           = 0
  global_meta%grid_fdda            = 0
  global_meta%gfdda_interval_m     = 0
  global_meta%gfdda_end_h          = 0
  global_meta%grid_sfdda           = 0  ! Added since WRFV3.1, since we do not use FFDA, they are ignored.
  global_meta%sgfdda_interval_m    = 0
  global_meta%sgfdda_end_h         = 0
  
  global_meta%grid_id        = grid_id        ! Added since WRFV2.1
  global_meta%parent_grid_ratio = grid_ratio
  IF (grid_id == 1) THEN
    global_meta%parent_id      = 0
    global_meta%i_parent_start = 1
    global_meta%j_parent_start = 1
  ELSE
    global_meta%grid_id        = grid_id        ! Added since WRFV2.1
    global_meta%parent_id      = parent_id
    global_meta%i_parent_start = i_parent_start
    global_meta%j_parent_start = j_parent_start
  END IF

  IF (global_meta%mminlu == 'UMD') THEN
    global_meta%iswater     = ISWATER_UMD
    global_meta%isice       = ISICE_UMD
    global_meta%isurban     = ISURBAN_UMD
  ELSE
    global_meta%num_land_cat = 24
    global_meta%islake      = -1
    global_meta%iswater     = ISWATER
    global_meta%isice       = ISICE
    global_meta%isurban     = ISURBAN
  END IF
  global_meta%isoilwater    = ISWATER_SOIL
  global_meta%map_proj      = map_proj
  
  RETURN
END SUBROUTINE set_global_meta
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE write_global_attribute         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_global_attribute(nfid,io_form,wrfversion,global_meta,bdyflg)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Initialize WRF global meta data and write as global attributes
!    to the NetCDF file.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER,                   INTENT(IN) :: nfid
  INTEGER,                   INTENT(IN) :: io_form
  REAL,                      INTENT(IN) :: wrfversion
  TYPE(wrf_global_metadata), INTENT(IN) :: global_meta
  LOGICAL,                   INTENT(IN) :: bdyflg

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER  :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (io_form == 7) CALL enter_ncd_define(nfid,istatus)

  CALL put_dom_ti_char(nfid,io_form,'TITLE',global_meta%title,istatus)
  CALL put_dom_ti_char(nfid,io_form,'START_DATE',                       &
                         global_meta%start_date(1:19),istatus)
  IF (.NOT. bdyflg) THEN
    CALL put_dom_ti_char(nfid,io_form,'SIMULATION_START_DATE', & !WRFV2.1
                         global_meta%start_date(1:19),istatus)
  END IF

  CALL put_dom_ti_integer(nfid,io_form,'WEST-EAST_GRID_DIMENSION',      &
                          global_meta%we_dimension, istatus)
  CALL put_dom_ti_integer(nfid,io_form,'SOUTH-NORTH_GRID_DIMENSION',    &
                          global_meta%sn_dimension, istatus)
  CALL put_dom_ti_integer(nfid,io_form,'BOTTOM-TOP_GRID_DIMENSION',     &
                          global_meta%bt_dimension, istatus)

  CALL put_dom_ti_real   (nfid,io_form, 'DX',global_meta%dx,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'DY',global_meta%dy,istatus)

  CALL put_dom_ti_char(nfid,io_form,'GRIDTYPE','C', istatus)
  
  IF (wrfversion < 3.1) THEN
  CALL put_dom_ti_integer(nfid,io_form,'DYN_OPT',                       &
                          global_meta%dyn_opt,istatus)
  END IF
  CALL put_dom_ti_integer(nfid,io_form,'DIFF_OPT',                      &
                          global_meta%diff_opt,istatus)
  CALL put_dom_ti_integer(nfid,io_form,'KM_OPT',                        &
                          global_meta%km_opt,istatus)
  CALL put_dom_ti_integer(nfid,io_form,'DAMP_OPT',                      &
                          global_meta%damp_opt,istatus)
  IF (wrfversion >= 3.1) THEN
    CALL put_dom_ti_real (nfid,io_form,'DAMPCOEF',                      &
                          global_meta%dampcoef,istatus)
  END IF
  CALL put_dom_ti_real   (nfid,io_form,'KHDIF',                         &
                          global_meta%khdif,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'KVDIF',                        &
                          global_meta%kvdif,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'MP_PHYSICS',                   &
                          global_meta%mp_physics,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'RA_LW_PHYSICS',                &
                          global_meta%ra_lw_physics,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'RA_SW_PHYSICS',                &
                          global_meta%ra_sw_physics,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SF_SFCLAY_PHYSICS',            &
                          global_meta%sf_sfclay_physics,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SF_SURFACE_PHYSICS',           &
                          global_meta%sf_surface_physics,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'BL_PBL_PHYSICS',               &
                          global_meta%bl_PBL_physics,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'CU_PHYSICS',                   &
                          global_meta%cu_physics,istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'SURFACE_INPUT_SOURCE',         &
                          global_meta%surface_input_source,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SST_UPDATE',                   &
                          global_meta%sst_update,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'GRID_FDDA',                    &
                          global_meta%grid_fdda,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'GFDDA_INTERVAL_M',             &
                          global_meta%gfdda_interval_m,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'GFDDA_END_H',                  &
                          global_meta%gfdda_end_h,istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'GRID_SFDDA',                   &
                          global_meta%grid_sfdda,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SGFDDA_INTERVAL_M',            &
                          global_meta%sgfdda_interval_m,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SGFDDA_END_H',                 &
                          global_meta%sgfdda_end_h,istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'WEST-EAST_PATCH_START_UNSTAG',&
                          global_meta%we_p_unstag_s,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'WEST-EAST_PATCH_END_UNSTAG',  &
                          global_meta%we_p_unstag_e,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'WEST-EAST_PATCH_START_STAG',  &
                          global_meta%we_p_stag_s,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'WEST-EAST_PATCH_END_STAG',    &
                          global_meta%we_p_stag_e,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SOUTH-NORTH_PATCH_START_UNSTAG',&
                          global_meta%sn_p_unstag_s,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SOUTH-NORTH_PATCH_END_UNSTAG', &
                          global_meta%sn_p_unstag_e,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SOUTH-NORTH_PATCH_START_STAG', &
                          global_meta%sn_p_stag_s,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'SOUTH-NORTH_PATCH_END_STAG',   &
                          global_meta%sn_p_stag_e,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'BOTTOM-TOP_PATCH_START_UNSTAG',&
                          global_meta%bt_p_unstag_s,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'BOTTOM-TOP_PATCH_END_UNSTAG',  &
                          global_meta%bt_p_unstag_e,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'BOTTOM-TOP_PATCH_START_STAG',  &
                          global_meta%bt_p_stag_s,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'BOTTOM-TOP_PATCH_END_STAG',    &
                          global_meta%bt_p_stag_e,istatus)

  !WRFV2.1
  CALL put_dom_ti_integer(nfid,io_form, 'GRID_ID',                      &
                          global_meta%grid_id,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'PARENT_ID',                    &
                          global_meta%parent_id,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'I_PARENT_START',               &
                          global_meta%i_parent_start,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'J_PARENT_START',               &
                          global_meta%j_parent_start,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'PARENT_GRID_RATIO',            &
                          global_meta%parent_grid_ratio,istatus)

  CALL put_dom_ti_real   (nfid,io_form, 'DT',global_meta%dt,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'CEN_LAT',                      &
                          global_meta%cen_lat,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'CEN_LON',                      &
                          global_meta%cen_lon,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'TRUELAT1',                     &
                          global_meta%tru_lat1,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'TRUELAT2',                     &
                          global_meta%tru_lat2,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'MOAD_CEN_LAT',                 &
                          global_meta%moad_cen_lat,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'STAND_LON',                    &
                          global_meta%stand_lon,istatus)

  IF (.NOT. bdyflg) THEN
    CALL put_dom_ti_real   (nfid,io_form, 'GMT',                        &
                            global_meta%gmt,istatus)
    CALL put_dom_ti_integer(nfid,io_form, 'JULYR',                      &
                            global_meta%julyr,istatus)
    CALL put_dom_ti_integer(nfid,io_form, 'JULDAY',                     &
                            global_meta%julday,istatus)
  END IF

  CALL put_dom_ti_integer(nfid,io_form, 'MAP_PROJ',                     &
                          global_meta%map_proj,istatus)
  CALL put_dom_ti_char(nfid,io_form,'MMINLU',global_meta%mminlu,istatus)

  IF (wrfversion >= 3.1) THEN
    CALL put_dom_ti_integer(nfid,io_form, 'NUM_LAND_CAT',               &
                            global_meta%num_land_cat,istatus)
    CALL put_dom_ti_integer(nfid,io_form, 'ISLAKE',                     &
                            global_meta%islake,istatus)
  END IF
  CALL put_dom_ti_integer(nfid,io_form, 'ISWATER',                      &
                          global_meta%iswater,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'ISICE',                        &
                          global_meta%isice,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'ISURBAN',                      &
                          global_meta%isurban,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'ISOILWATER',                   &
                          global_meta%isoilwater,istatus)
  
  IF (io_form == 7) CALL exit_ncd_define(nfid,istatus)

  RETURN
END SUBROUTINE write_global_attribute

SUBROUTINE put_dom_ti_char(nfid,io_form,attname,attstr,istatus)

  IMPLICIT NONE

  INTEGER,      INTENT(IN) :: nfid
  INTEGER,      INTENT(IN) :: io_form
  CHARACTER(*), INTENT(IN) :: attname
  CHARACTER(*), INTENT(IN) :: attstr
  INTEGER,      INTENT(OUT):: istatus

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code  ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (io_form)
  CASE (1)     ! binary
    IF(myproc == 0) CALL ext_int_put_dom_ti_char(nfid,attname,attstr,   &
                                  iStatus )
  CASE (5)     ! PHDF5
    CALL put_phdf5_dom_ti_char(nfid,attname,attstr,istatus)
  CASE (7)     ! NetCDF
    IF(myproc == 0) CALL put_ncd_dom_ti_char(nfid,attname,attstr,istatus)
  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  RETURN
END SUBROUTINE put_dom_ti_char

SUBROUTINE put_dom_ti_integer(nfid,io_form,attname,attval,istatus)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nfid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: attname
  INTEGER,      INTENT(IN)  :: attval
  INTEGER,      INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code  ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (io_form)
  CASE (1)     ! binary
    IF(myproc == 0) CALL ext_int_put_dom_ti_integer(nfid,attname,attval,&
                                                    1,iStatus )
  CASE (5)     ! PHDF5
    CALL put_phdf5_dom_ti_integer(nfid,attname,attval,istatus)
  CASE (7)     ! NetCDF
    IF(myproc == 0) CALL put_ncd_dom_ti_integer(nfid,attname,attval,istatus)
  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  RETURN
END SUBROUTINE put_dom_ti_integer

SUBROUTINE put_dom_ti_real(nfid,io_form,attname,attval,istatus)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nfid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: attname
  REAL,         INTENT(IN)  :: attval
  INTEGER,      INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code  ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (io_form)
  CASE (1)     ! binary
    IF (myproc == 0) CALL ext_int_put_dom_ti_real(nfid,attname,attval,  &
                                                   1,iStatus )
  CASE (5)     ! PHDF5
    CALL put_phdf5_dom_ti_real(nfid,attname,attval,1,istatus)
  CASE (7)     ! NetCDF
    IF (myproc == 0) CALL put_ncd_dom_ti_real(nfid,attname,attval,      &
                                                   1,istatus)
  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  RETURN
END SUBROUTINE put_dom_ti_real

SUBROUTINE put_dom_ti_varreal(nfid,io_form,attname,attval,attsiz,istatus)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: nfid
  INTEGER,      INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: attname
  INTEGER,      INTENT(IN)  :: attsiz
  REAL,         INTENT(IN)  :: attval(attsiz)
  INTEGER,      INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code  ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (io_form)
  CASE (1)     ! binary
    IF (myproc == 0) CALL ext_int_put_dom_ti_real(nfid,attname,attval,  &
                                                  attsiz,iStatus )
  CASE (5)     ! PHDF5
    CALL put_phdf5_dom_ti_real(nfid,attname,attval,attsiz,istatus)
  CASE (7)     ! NetCDF
    IF (myproc == 0) CALL put_ncd_dom_ti_real(nfid,attname,attval,   &
                                                 attsiz,istatus)
  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  RETURN
END SUBROUTINE put_dom_ti_varreal
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE write_times_str                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_times_str(nfid,io_form,varname,currDate,DateStr,ntime)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,           INTENT(IN)  :: nfid
  INTEGER,           INTENT(IN)  :: io_form
  CHARACTER(*),      INTENT(IN)  :: varname
  CHARACTER(*),      INTENT(IN)  :: CurrDate
  CHARACTER(*),      INTENT(IN)  :: DateStr
  INTEGER,           INTENT(IN)  :: ntime

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus
  CHARACTER(120)  :: tmpname, lowName

  INTEGER, PARAMETER :: upper_to_lower =IACHAR('a')-IACHAR('A')
  CHARACTER(1)       :: c
  INTEGER            :: i,namelen

  INCLUDE 'mp.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  SELECT CASE (io_form)
  CASE (1)     ! binary
     IF (myproc == 0) CALL ext_int_put_dom_td_char(nfid,varname,CurrDate, &
                                                   DateStr,iStatus )
  CASE (5)     ! PHDF5
    CALL put_phdf5_dom_td_char(nfid,varname,CurrDate,DateStr,istatus) 
  CASE (7)     ! NetCDF
    IF ( varname == 'Times' ) THEN
      tmpname = varname
    ELSE
      namelen = LEN_TRIM(varName)
      DO i=1,namelen
        c = varName(i:i)
        IF( 'A' <= c .AND. c <= 'Z') THEN
          lowName(i:i) = ACHAR(IACHAR(c)+upper_to_lower)
        ELSE IF( c == '-' .OR.  c == ':') THEN
          lowName(i:i) = '_'
        ELSE
          lowName(i:i) = c
        END IF
      END DO
      tmpname = 'md___'//lowName(1:namelen)//'e_x_t_d_o_m_a_i_n_m_e_t_a_data_'
    END IF
    IF (myproc == 0) CALL put_ncd_dom_td_char(nfid,tmpname,             &
                                              DateStr,ntime,istatus) 
  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  RETURN
END SUBROUTINE write_times_str
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write_real               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_real(nout,io_form,var_meta,DateStr,var,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D vector to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,                INTENT(IN)  :: nout
  INTEGER,                INTENT(IN)  :: io_form
  TYPE(wrf_var_metadata), INTENT(IN)  :: var_meta
  CHARACTER(*),           INTENT(IN)  :: DateStr
  REAL,                   INTENT(IN)  :: var
  INTEGER,                INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)  &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing REAL scalar ', var_meta%name

  DimNames(1) = var_meta%dimName1
  DimNames(2) = var_meta%dimName2
  DimNames(3) = var_meta%dimName3

  DomainStart(:) = 1
  DomainEnd(1)   = 1
  DomainEnd(2:3) = 1

  MemoryStart(:) = 1
  MemoryEnd(1)   = 1
  MemoryEnd(2:3) = 1

  PatchStart(:) = 1
  PatchEnd(1)   = 1
  PatchEnd(2:3) = 1

  SELECT CASE (io_form)
  CASE (1)     ! binary
  
    IF (myproc == 0) CALL ext_int_write_field( nout, DateStr, Var_meta%Name, &
                 var, var_meta%FieldType, 0, 0, 1,                      &
                 var_meta%MemoryOrder, var_meta%Stagger, DimNames,      &
                 DomainStart, DomainEnd, MemoryStart, MemoryEnd,        &
                 DomainStart, DomainEnd, iStatus )

  CASE (5)     ! PHDF5

    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var,var_meta%fieldType, 1,                   &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)

  CASE (7)     ! NetCDF

    IF (myproc == 0) CALL write_ncd_real(nout,TRIM(var_meta%name),        &
                                         var,istatus)

  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  IF (myproc == 0)  THEN
    IF (istatus == 0) THEN
      WRITE(6,'(a)') '         ...  DONE.'
    ELSE
      WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE write_real
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write_int                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_int(nout,io_form,var_meta,DateStr,var,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D vector to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,                INTENT(IN)  :: nout
  INTEGER,                INTENT(IN)  :: io_form
  TYPE(wrf_var_metadata), INTENT(IN)  :: var_meta
  CHARACTER(*),           INTENT(IN)  :: DateStr
  INTEGER,                INTENT(IN)  :: var
  INTEGER,                INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)  &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing INT  scalar ', var_meta%name

  DimNames(1) = var_meta%dimName1
  DimNames(2) = var_meta%dimName2
  DimNames(3) = var_meta%dimName3

  DomainStart(:) = 1
  DomainEnd(1)   = 1
  DomainEnd(2:3) = 1

  MemoryStart(:) = 1
  MemoryEnd(1)   = 1
  MemoryEnd(2:3) = 1

  PatchStart(:) = 1
  PatchEnd(1)   = 1
  PatchEnd(2:3) = 1

  SELECT CASE (io_form)
  CASE (1)     ! binary
  
    IF (myproc == 0) CALL ext_int_write_field( nout, DateStr, Var_meta%Name, &
                 var, var_meta%FieldType, 0, 0, 1,                      &
                 var_meta%MemoryOrder, var_meta%Stagger, DimNames,      &
                 DomainStart, DomainEnd, MemoryStart, MemoryEnd,        &
                 DomainStart, DomainEnd, iStatus )

  CASE (5)     ! PHDF5

    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var,var_meta%fieldType, 1,                   &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)

  CASE (7)     ! NetCDF

    IF (myproc == 0) CALL write_ncd_int(nout,TRIM(var_meta%name),        &
                                        var,istatus)

  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  IF (myproc == 0)  THEN
    IF (istatus == 0) THEN
      WRITE(6,'(a)') '         ...  DONE.'
    ELSE
      WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE write_int
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write1d                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write1d(nout,io_form,var_meta,DateStr,var1d,nz,fzone)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D vector to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,                INTENT(IN) :: nout
  INTEGER,                INTENT(IN) :: io_form
  TYPE(wrf_var_metadata), INTENT(IN) :: var_meta
  CHARACTER(*),           INTENT(IN) :: DateStr
  INTEGER,                INTENT(IN) :: nz
  INTEGER,                INTENT(IN) :: fzone
  REAL,                   INTENT(IN) :: var1d(nz)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)  &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 1D variable ', var_meta%name

  DimNames(1) = var_meta%dimName1
  DimNames(2) = var_meta%dimName2
  DimNames(3) = var_meta%dimName3

  DomainStart(:) = 1
  DomainEnd(1)   = nz
  DomainEnd(2:3) = 1

  MemoryStart(:) = 1
  MemoryEnd(1)   = nz
  MemoryEnd(2:3) = 1

  PatchStart(:) = 1
  PatchEnd(1)   = nz
  PatchEnd(2:3) = 1

  SELECT CASE (io_form)
  CASE (1)     ! binary
  
    IF (myproc == 0) CALL ext_int_write_field( nout, DateStr, Var_meta%Name, &
                 var1d, var_meta%FieldType, 0, 0, 1,                    &
                 var_meta%MemoryOrder, var_meta%Stagger, DimNames,      &
                 DomainStart, DomainEnd, MemoryStart, MemoryEnd,        &
                 DomainStart, DomainEnd, iStatus )

  CASE (5)     ! PHDF5

    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var1d,var_meta%fieldType, 1,                 &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)

  CASE (7)     ! NetCDF

    IF (myproc == 0) CALL write_ncd_1d(nout,TRIM(var_meta%name),        &
                                       var1d,nz,istatus)

  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  IF (myproc == 0)  THEN
    IF (istatus == 0) THEN
      WRITE(6,'(a)') '         ...  DONE.'
    ELSE
      WRITE(6,'(a)') '         ...  ERROR.'
    END IF
  END IF

  RETURN
END SUBROUTINE write1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write1di                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write1di(nout,io_form,var_meta,DateStr,var1di,nz,fzone)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D integer vector to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,                INTENT(IN) :: nout
  INTEGER,                INTENT(IN) :: io_form
  TYPE(wrf_var_metadata), INTENT(IN) :: var_meta
  CHARACTER(*),           INTENT(IN) :: DateStr
  INTEGER,                INTENT(IN) :: fzone
  INTEGER,                INTENT(IN) :: nz
  INTEGER,                INTENT(IN) :: var1di(nz)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)  &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 1D integer variable ', var_meta%name

  DimNames(1) = var_meta%dimName1
  DimNames(2) = var_meta%dimName2
  DimNames(3) = var_meta%dimName3

  DomainStart(:) = 1
  DomainEnd(1)   = nz
  DomainEnd(2:3) = 1

  MemoryStart(:) = 1
  MemoryEnd(1)   = nz
  MemoryEnd(2:3) = 1

  PatchStart(:) = 1
  PatchEnd(1)   = nz
  PatchEnd(2:3) = 1

  SELECT CASE (io_form)
  CASE (1)     ! binary

    IF (myproc == 0) CALL ext_int_write_field( nout, DateStr, Var_meta%Name, &
                 var1di, var_meta%FieldType, 0, 0, 1,                   &
                 var_meta%MemoryOrder, var_meta%Stagger, DimNames,      &
                 DomainStart, DomainEnd, MemoryStart, MemoryEnd,        &
                 DomainStart, DomainEnd, iStatus )

  CASE (5)     ! PHDF5

    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var1di,var_meta%fieldType, 1,                &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)

  CASE (7)     ! NetCDF

    IF (myproc == 0) CALL write_ncd_1di(nout,TRIM(var_meta%name),       &
                                        var1di,nz,istatus)

  CASE DEFAULT
    WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
    istatus = -1
  END SELECT

  IF (myproc == 0 .AND. istatus == 0) WRITE(6,'(a)') ' ...  DONE.'

  RETURN
END SUBROUTINE write1di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write2d                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write2d(nout,io_form,var_meta,DateStr,var2d,nx,ny,           &
                   fzone,temlg,nxlg,nylg)
!
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 2D array to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER,                INTENT(IN) :: nout
  INTEGER,                INTENT(IN) :: io_form
  TYPE(wrf_var_metadata), INTENT(IN) :: var_meta
  CHARACTER(*),           INTENT(IN) :: DateStr
  INTEGER,                INTENT(IN) :: nx,ny
  INTEGER,                INTENT(IN) :: nxlg,nylg
  INTEGER,                INTENT(IN) :: fzone
  REAL,                   INTENT(IN) :: var2d(nx,ny)
  REAL,                   INTENT(IN) :: temlg(nxlg,nylg)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus
 
  INCLUDE 'mp.inc'
 
  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
 
  INTEGER       :: ilocs,iloce,jlocs,jloce

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)  &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 2D variable ', var_meta%name

  DimNames(1) = var_meta%dimName1
  DimNames(2) = var_meta%dimName2
  DimNames(3) = var_meta%dimName3
 
  DomainStart(:) = 1
  DomainEnd(1)   = nxlg
  DomainEnd(2)   = nylg
  DomainEnd(3)   = 1                  ! 2d special
 
  IF (io_form == 5) THEN              ! No merge

    ilocs = (nx-fzone)*(loc_x-1)+fzone
    jlocs = (ny-fzone)*(loc_y-1)+fzone
    iloce = (nx-fzone)*(loc_x)+fzone
    jloce = (ny-fzone)*(loc_y)+fzone
 
    MemoryStart(1) = ilocs
    MemoryStart(2) = jlocs
    MemoryStart(3) = 1
    MemoryEnd(1)   = iloce
    MemoryEnd(2)   = jloce
    MemoryEnd(3)   = 1                  ! 2d special

    PatchStart(1) = ilocs
    PatchEnd(1)   = iloce - fzone
    IF (var_meta%stagger == 'X') THEN
      IF (loc_x > 1) PatchStart(1) = ilocs + fzone
      PatchEnd(1)   =  iloce
    END IF
 
    PatchStart(2) = jlocs
    PatchEnd(2)   = jloce - fzone
    IF (var_meta%stagger == 'Y') THEN
      IF (loc_y > 1 ) PatchStart(2) = jlocs + fzone
      PatchEnd(2)   =  jloce
    END IF
 
    PatchStart(3) = 1
    PatchEnd(3)   = 1                  ! 2d special
 
    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var2d,var_meta%fieldType, 1,                 &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)


  ELSE                 ! Need merge

    IF (var_meta%stagger == 'X') THEN
      CALL wrf_merge2du(var2d,nx,ny,fzone,temlg)
    ELSE IF (var_meta%stagger == 'Y') THEN
      CALL wrf_merge2dv(var2d,nx,ny,fzone,temlg)
    ELSE
       CALL wrf_merge2dt(var2d,nx,ny,fzone,temlg)
    END IF
 
    IF (io_form == 1) THEN         ! Binary

      IF (myproc == 0) CALL ext_int_write_field( nout, DateStr,         &
                 Var_meta%Name, temlg, var_meta%FieldType, 0, 0, 1,     &
                 var_meta%MemoryOrder, var_meta%Stagger, DimNames,      &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )
                 ! temlg contain domain data, so Domain, Memeory &
                 ! Patch index are the same

    ELSE IF (io_form == 7) THEN    ! NetCDF
 
      IF (myproc == 0) CALL write_ncd_2d(nout,TRIM(var_meta%name),      &
                 temlg,nxlg,nylg,istatus)
 
    ELSE
      WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
      istatus = -1
    END IF

  END IF

  IF (myproc == 0 .AND. istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write2di                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write2di(nout,io_form,var_meta,DateStr,var2di,nx,ny,         &
                   fzone,temlg,nxlg,nylg)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 2D array to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,                INTENT(IN) :: nout
  INTEGER,                INTENT(IN) :: io_form
  TYPE(wrf_var_metadata), INTENT(IN) :: var_meta
  CHARACTER(*),           INTENT(IN) :: DateStr
  INTEGER,                INTENT(IN) :: nx,ny
  INTEGER,                INTENT(IN) :: nxlg,nylg
  INTEGER,                INTENT(IN) :: fzone
  INTEGER,                INTENT(IN) :: var2di(nx,ny)
  INTEGER,                INTENT(IN) :: temlg(nxlg,nylg)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: ilocs,iloce,jlocs,jloce
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)  &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 2D integer variable ', var_meta%name

  DimNames(1) = var_meta%dimName1
  DimNames(2) = var_meta%dimName2
  DimNames(3) = var_meta%dimName3

  DomainStart(:) = 1
  DomainEnd(1)   = nxlg
  DomainEnd(2)   = nylg
  DomainEnd(3)   = 1               ! 2d special

  IF (io_form == 5) THEN           ! No merge

    ilocs = (nx-fzone)*(loc_x-1)+fzone
    jlocs = (ny-fzone)*(loc_y-1)+fzone
    iloce = (nx-fzone)*(loc_x)+fzone
    jloce = (ny-fzone)*(loc_y)+fzone

    MemoryStart(1) = ilocs
    MemoryStart(2) = jlocs
    MemoryStart(3) = 1
    MemoryEnd(1)   = iloce
    MemoryEnd(2)   = jloce
    MemoryEnd(3)   = 1

    PatchStart(1) = ilocs
    PatchEnd(1)   = iloce - fzone
    IF (var_meta%stagger == 'X') THEN
      IF (loc_x > 1) PatchStart(1) = ilocs + fzone
      PatchEnd(1)   =  iloce
    END IF

    PatchStart(2) = jlocs
    PatchEnd(2)   = jloce - fzone
    IF (var_meta%stagger == 'Y') THEN
      IF (loc_y > 1 ) PatchStart(2) = jlocs + fzone
      PatchEnd(2)   =  jloce
    END IF

    PatchStart(3) = 1
    PatchEnd(3)   = 1                  ! 2d special

    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var2di,var_meta%fieldType, 1,                &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)

  ELSE                  ! Need merge

    IF (var_meta%stagger == 'X') THEN
      WRITE(0,*) 'WARNING: Need to implement merger for integer 2D variable.'
    ELSE IF (var_meta%stagger == 'Y') THEN
      WRITE(0,*) 'WARNING: Need to implement merger for integer 2D variable.'
    ELSE
       CALL wrf_merge2di(var2di,nx,ny,fzone,temlg)
    END IF

    IF (io_form == 1) THEN        ! Binary

      IF (myproc == 0) CALL ext_int_write_field( nout, DateStr,         &
                 Var_meta%Name, temlg, var_meta%FieldType, 0, 0, 1,     &
                 var_meta%MemoryOrder, var_meta%Stagger, DimNames,      &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )
                 ! temlg contain domain data, so Domain, Memeory &
                 ! Patch index are the same

    ELSE IF (io_form == 7) THEN   ! NetCDF

      IF (myproc == 0) CALL write_ncd_2di(nout,TRIM(var_meta%name),     &
                                 temlg,nxlg,nylg,istatus)

    ELSE 
      WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
      istatus = -1
    END IF

  END IF

  IF (myproc == 0 .AND. istatus == 0) WRITE(6,'(a)') ' ...  DONE.'

  RETURN
END SUBROUTINE write2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE write3d                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write3d(nout,io_form,var_meta,DateStr,var3d,nx,ny,nz,       &
                   fzone,temlg,temxzylg,nxlg,nylg,nzlg)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 3D array to the output file.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,                INTENT(IN) :: nout   
  INTEGER,                INTENT(IN) :: io_form
  TYPE(wrf_var_metadata), INTENT(IN) :: var_meta
  CHARACTER(*),           INTENT(IN) :: DateStr
  INTEGER,                INTENT(IN) :: nx,ny,nz
  INTEGER,                INTENT(IN) :: nxlg,nylg,nzlg
  INTEGER,                INTENT(IN) :: fzone
  REAL,                   INTENT(IN) :: var3d(nx,ny,nz)

  REAL,                  INTENT(OUT) :: temlg(nxlg,nylg,nzlg)
  REAL,                  INTENT(OUT) :: temxzylg(nxlg,nzlg,nylg)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus

  INCLUDE 'mp.inc'

  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)

  INTEGER       :: ilocs,iloce,jlocs,jloce
  INTEGER       :: i, j, k
  CHARACTER(3)  :: MemOrder
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF(myproc == 0)   &
    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 3D variable ', var_meta%name

  IF (io_form == 5) THEN             ! No merge

    DimNames(1) = var_meta%dimName1
    DimNames(2) = var_meta%dimName2
    DimNames(3) = var_meta%dimName3

    DomainStart(1) = 1
    DomainStart(2) = 1
    DomainStart(3) = 1
    DomainEnd(1)   = nxlg
    DomainEnd(2)   = nylg
    DomainEnd(3)   = nzlg

    ilocs = (nx-fzone)*(loc_x-1)+fzone
    jlocs = (ny-fzone)*(loc_y-1)+fzone
    iloce = (nx-fzone)*(loc_x)+fzone
    jloce = (ny-fzone)*(loc_y)+fzone

    MemoryStart(1) = ilocs
    MemoryStart(2) = jlocs
    MemoryStart(3) = 1
    MemoryEnd(1)   = iloce
    MemoryEnd(2)   = jloce
    MemoryEnd(3)   = nz
   
    PatchStart(1) = ilocs
    PatchEnd(1)   = iloce - fzone
    IF (var_meta%stagger == 'X') THEN
      IF (loc_x > 1) PatchStart(1) = ilocs + fzone
      PatchEnd(1)   =  iloce
    END IF

    PatchStart(2) = jlocs
    PatchEnd(2)   = jloce - fzone
    IF (var_meta%stagger == 'Y') THEN
      IF (loc_y > 1 ) PatchStart(2) = jlocs + fzone
      PatchEnd(2)   =  jloce
    END IF

    PatchStart(3) = 1
    PatchEnd(3)   = nzlg

    CALL write_phdf5_field(nout,DateStr,var_meta%name,var_meta%description, &
                           var_meta%units,var_meta%stagger,             &
                           var3d,var_meta%fieldType, 1,                 &
                           var_meta%memoryOrder,DimNames,               &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, istatus)

  ELSE                 ! Need merge
  
    IF (var_meta%stagger == 'X') THEN
       CALL wrf_merge3du(var3d,nx,ny,nz,fzone,temlg)
    ELSE IF (var_meta%stagger == 'Y') THEN
       CALL wrf_merge3dv(var3d,nx,ny,nz,fzone,temlg)
    ELSE IF (var_meta%stagger == 'Z') THEN
       CALL wrf_merge3dw(var3d,nx,ny,nz,fzone,temlg)
    ELSE
       CALL wrf_merge3dt(var3d,nx,ny,nz,fzone,temlg)
    END IF

    IF (io_form == 1) THEN            ! Binary

      MemOrder = 'XZY'

      DimNames(1) = var_meta%dimName1         ! X
      DimNames(2) = var_meta%dimName3         ! Z
      DimNames(3) = var_meta%dimName2         ! Y

      DomainStart(1) = 1
      DomainStart(2) = 1
      DomainStart(3) = 1
      DomainEnd(1)   = nxlg
      DomainEnd(2)   = nzlg
      DomainEnd(3)   = nylg

      DO k = 1,nzlg
        DO j = 1,nylg
          DO i = 1,nxlg
            temxzylg(i,k,j) = temlg(i,j,k)
          END DO
        END DO
      END DO

      IF (myproc == 0) CALL ext_int_write_field( nout, DateStr,         &
                 Var_meta%Name, temxzylg, var_meta%FieldType, 0, 0, 1,  &
                 MemOrder, var_meta%Stagger, DimNames,                  &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )
                 ! temlg contain domain data, so Domain, Memeory &
                 ! Patch index are the same

    ELSE IF (io_form == 7) THEN       ! NetCDF

      IF (myproc == 0) CALL write_ncd_3d(nout,TRIM(var_meta%name),      &
                 temlg,nxlg,nylg,nzlg,istatus)

    ELSE
      WRITE(0,*) 'ERROR: unsupport IO format: ',io_form
      istatus = -1
    END IF

  END IF

  IF (myproc == 0 .AND. istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write3d
!
!//////////////////////////////////////////////////////////////////
!
! Read WRF static file
!
!//////////////////////////////////////////////////////////////////
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE open_static_file             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE open_static_file(filename,ncid)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!     Open the static file for read
!
!------------------------------------------------------------------
  IMPLICIT NONE

  CHARACTER(LEN=*),   INTENT(IN)  :: filename
  INTEGER,            INTENT(OUT) :: ncid
!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus
  LOGICAL :: static_exists

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  INQUIRE(FILE=filename, EXIST=static_exists)
  IF (static_exists) THEN
    istatus = NF_OPEN(TRIM(filename),NF_NOWRITE,ncid)
    CALL nf_handle_error(istatus,'NF_OPEN')
  ELSE
    PRINT '(A)', 'Static file not found: ', filename
    STOP 'open_wrfsi_static'
  ENDIF

  RETURN
END SUBROUTINE open_static_file

!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE close_static_file            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
! 
SUBROUTINE close_static_file(ncid)
! 
!------------------------------------------------------------------
!
!  PURPOSE:
!  
!     Close the static file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ncid

  INCLUDE 'netcdf.inc'
!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
!
  INTEGER :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = NF_CLOSE(ncid)
  CALL nf_handle_error(istatus,'NF_CLOSE')

  RETURN
END SUBROUTINE close_static_file

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE check_static_grid            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
! 
SUBROUTINE check_static_grid(staticopt,ncid,nxin,nyin,dxin,dyin,        &
                mapprojin,trulat1in,trulat2in,trulonin,istatus)
! 
!------------------------------------------------------------------
!
!  PURPOSE:
!
!     Read the grid information from the static NetCDF file, and
!     Check whether it is at the same grid as the input map projection.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: staticopt           ! 0: WRFSI static file
                                              ! 1: WRF static file
  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: nxin,nyin
  REAL,    INTENT(IN)  :: dxin,dyin
  INTEGER, INTENT(IN)  :: mapprojin
  REAL,    INTENT(IN)  :: trulat1in,trulat2in
  REAL,    INTENT(IN)  :: trulonin

  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: nx, ny, iproj
  INTEGER :: nxlg, nylg
  REAL    :: dx, dy
  REAL    :: trulat1,trulat2,trulon
  REAL    :: lat1,lon1
  INTEGER :: mapproj

  REAL, PARAMETER :: eps = 0.00001

  INTEGER :: vid
  CHARACTER(LEN=132) :: grid_type

  INCLUDE 'netcdf.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (staticopt == 0) THEN
    istatus = NF_INQ_VARID(ncid, 'Nx', vid)
    CALL nf_handle_error(istatus,'get_static_grid, Nx')
    istatus = NF_GET_VAR_INT(ncid,vid,nx)
    CALL nf_handle_error(istatus,'get_static_grid, Nx')

    istatus = NF_INQ_VARID(ncid, 'Ny', vid)
    CALL nf_handle_error(istatus,'get_static_grid, NF_INQ_VARID')
    istatus = NF_GET_VAR_INT(ncid,vid,ny)
    CALL nf_handle_error(istatus,'get_static_grid, NF_GET_VAR_INT')

    istatus = NF_INQ_VARID(ncid, 'Dx', vid)
    CALL nf_handle_error(istatus,'get_static_grid, NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,dx)
    CALL nf_handle_error(istatus,'get_static_grid, NF_GET_VAR_REAL')

    istatus = NF_INQ_VARID(ncid, 'Dy', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,dy)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

    istatus = NF_INQ_VARID(ncid, 'La1', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,lat1)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

    istatus = NF_INQ_VARID(ncid, 'Lo1', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,lon1)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')
    IF(lon1 > 180) lon1 = lon1 - 360

    istatus = NF_INQ_VARID(ncid, 'LoV', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,trulon)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

    IF(trulon > 180) trulon = trulon - 360

    istatus = NF_INQ_VARID(ncid, 'Latin1', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,trulat1)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

    istatus = NF_INQ_VARID(ncid, 'Latin2', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_REAL(ncid,vid,trulat2)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

    istatus = NF_INQ_VARID(ncid, 'grid_type', vid)
    CALL nf_handle_error(istatus,'NF_INQ_VARID')
    istatus = NF_GET_VAR_TEXT(ncid,vid,grid_type)
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

    IF( INDEX(grid_type,'polar') > 0 )             mapproj = 1
    IF( INDEX(grid_type,'lambert conformal') > 0 ) mapproj = 2
    IF( INDEX(grid_type,'mercator') > 0 )          mapproj = 3

  ELSE IF (staticopt == 1) THEN    ! WRF static file

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WEST-EAST_GRID_DIMENSION', nx)
    CALL nf_handle_error(istatus,'get_static_grid, WEST-EAST_GRID_DIMENSION.')

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SOUTH-NORTH_GRID_DIMENSION',ny)
    CALL nf_handle_error(istatus,'get_static_grid, SOUTH-NORTH_GRID_DIMENSION.')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX', dx)
    CALL nf_handle_error(istatus,'get_static_grid, DX.')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY', dy)
    CALL nf_handle_error(istatus,'get_static_grid, DY.')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1', trulat1)
    CALL nf_handle_error(istatus,'get_static_grid, TRUELAT1.')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2', trulat2)
    CALL nf_handle_error(istatus,'get_static_grid, TRUELAT2.')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'STAND_LON', trulon)
    CALL nf_handle_error(istatus,'get_static_grid, STAND_LON.')
    IF(trulon > 180) trulon = trulon - 360

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAP_PROJ', iproj)
    CALL nf_handle_error(istatus,'get_static_grid, MAP_PROJ.')

    IF(iproj == 0) THEN        ! No projection
      mapproj = 0
    ELSE IF(iproj == 1) THEN   ! LAMBERT CONFORMAL
      mapproj = 2
    ELSE IF(iproj == 2) THEN   ! POLAR STEREOGRAPHIC
      mapproj = 1
    ELSE IF(iproj == 3) THEN   ! MERCATOR
      mapproj = 3
    ELSE
      WRITE(6,*) 'Unknown map projection, ', iproj
      istatus = -2
      RETURN
    END IF

  ELSE

    WRITE(6,'(a)')   &
      'Wrong staticopt (0 for WRFSI static file, 1 for WRF static file).'

    istatus = -3
    RETURN
  END IF

! Determine whether we are at the same grid with the static file
  nxlg = (nxin - 1)*nproc_x + 1
  nylg = (nyin - 1)*nproc_y + 1
  istatus = 0

  IF(nxlg /= nx .OR. nylg /= ny .OR.                          &
     ABS(dxin-dy) > eps .OR. ABS(dyin-dy) > eps .OR.              &
     ABS(mapprojin)  /= mapproj  .OR.                             &
     ABS(trulonin-trulon)  > eps .OR.                             &
     ABS(trulat1in-trulat1)> eps .OR.                             &
     ABS(trulat2in-trulat2)> eps ) THEN

    WRITE(6,*) 'Grid info. from static file', '     WRF grid'
    WRITE(6,*) '          =================', '     ========'
    WRITE(6,'(4x,a,I7,  9x,  I7)')   'nx:             ',nx,nxlg
    WRITE(6,'(4x,a,I7,  9x,  I7)')   'ny:             ',ny,nylg
    WRITE(6,'(4x,a,F7.0,9x,  F7.0)') 'dx:             ',dx,dxin
    WRITE(6,'(4x,a,F7.0,9x,  F7.0)') 'dy:             ',dy,dyin
    WRITE(6,'(4x,a,I7,  5x,a,I7)')   'mapproj:        ',mapproj,  &
                                     '    ',mapprojin
    WRITE(6,'(4x,a,F7.2,5x,a,F7.2)') 'trulat1:        ',trulat1,  &
                                     '    ',trulat1in
    WRITE(6,'(4x,a,F7.2,5x,a,F7.2)') 'trulat2:        ',trulat2,  &
                                    '    ',trulat2in
    WRITE(6,'(4x,a,F7.2,5x,a,F7.2)') 'trulon:         ',trulon,   &
                                   '    ',trulonin
    WRITE(6,'(4x,a,F7.2,5x,a,F7.2)') 'SW corner lat.: ',lat1
    WRITE(6,'(4x,a,F7.2,5x,a,F7.2)') 'SW corner lon.: ',lon1

    WRITE(6,*) 'ERROR: Static NetCDF file is not at the same grid',      &
               ' as WRF grid.'

    istatus = -1

  END IF

  RETURN
END SUBROUTINE check_static_grid
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_static_landusef          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_static_landusef(ncid, var3d)
    
! Reads the individual 2D categorical landuse fraction arrays from
! the static file and populates the 3D variable.
  
  USE wrf_metadata
  
  IMPLICIT NONE
  INTEGER, INTENT(IN)         :: ncid
  REAL,    INTENT(OUT)        :: var3d(:,:,:)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  INTEGER            :: cat
  CHARACTER(LEN=3)   :: varname
  INTEGER            :: varid

  INTEGER            :: istatus
  REAL               :: fillValue

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! For now, we only support 16-category WMO/FAO data set, so hard code the
! info

!  PRINT '(A)', 'Attempting to read categorical landuse fractions...'

  DO cat = 1, LanduseCategories
    WRITE(varname,'("u",I2.2)') cat
    istatus = NF_INQ_VARID(ncid, varname, varid)
    istatus = NF_GET_VAR_REAL(ncid,varid,var3d(:,:,cat))
    CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

!    istatus = NF_GET_ATT_REAL(ncid,varid,'_FillVale',fillValue)
!    WHERE(var3d(:,:,cat) >= fillValue) var3d(:,:,cat) = rmissing
  ENDDO
  
  RETURN
END SUBROUTINE get_static_landusef
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_static_soil              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_static_soil(ncid,static_soiltop,static_soilbot)

! Reads the 2-layer soil categorical fractions and populates the arrays
! static_soiltop and static_soilbot
  
  USE wrf_metadata
  
  IMPLICIT NONE
  INTEGER, INTENT(IN)         :: ncid
  INTEGER, INTENT(OUT)        :: static_soiltop(:,:,:)
  INTEGER, INTENT(OUT)        :: static_soilbot(:,:,:)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  INTEGER                :: vid,istatus
  INTEGER                :: cat
  CHARACTER(LEN=3)       :: vname
  REAL                   :: fillValue

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! For now, we only support 16-category WMO/FAO data set, so hard code the
! info

!  PRINT '(A)', 'Attempting to read categorical soil type fractions...'

  DO cat = 1, SoilCategories
    WRITE(vname,'("b",I2.2)') cat
    istatus = NF_INQ_VARID(ncid, vname, vid)
    istatus = NF_GET_VAR_REAL(ncid,vid,static_soilbot(:,:,cat))
    CALL nf_handle_error(istatus,'get_soilbot')

!    istatus = NF_GET_ATT_REAL(ncid,vid,'_FillVale',fillValue)
!    WHERE(static_soilbot(:,:,cat) >= fillValue) static_soilbot(:,:,cat) = rmissing

    WRITE(vname,'("t",I2.2)') cat
    istatus = NF_INQ_VARID(ncid, vname, vid)
    istatus = NF_GET_VAR_REAL(ncid,vid,static_soiltop(:,:,cat))
    CALL nf_handle_error(istatus,'get_soiltop')

!    istatus = NF_GET_ATT_REAL(ncid,vid,'_FillVale',fillValue)
!    WHERE(static_soiltop(:,:,cat) >= fillValue) static_soiltop(:,:,cat) = rmissing
  END DO
  RETURN
END SUBROUTINE get_static_soil
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_static_monthly           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_static_monthly(staticopt, ncid, vartype, valid_day,      &
                              var2d, varin3d, istatus)

! Returns a time-interpolated (valid for time) 2D array for either
! greenness ("g") or albedo ("a"), based on user-supplied type character.
! The monthly values are 
! valid on the 15th day of each month.  This routine only interpolates to the
! nearest day and does not account for leap years, but this should not be any
! big deal.
  
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: staticopt
  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: valid_day    ! julday of the year
  CHARACTER(LEN=1), INTENT(IN)  :: vartype
  REAL,             INTENT(OUT) :: var2d(:,:)
  REAL,             INTENT(IN)  :: varin3d(:,:,:)
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  INTEGER  :: m, d1, d2, m1, m2
  INTEGER  :: nx,ny
  CHARACTER(LEN=3)  :: varname
  REAL, ALLOCATABLE :: data1(:,:), data2(:,:)
  REAL              :: w1, w2
  INTEGER           :: i, j

  INTEGER :: midmonth_day(12)
  ! midmonth_day is the julian day of the year corresponding to the 15th day
  ! of each month for a standard (non-leap) year
  DATA midmonth_day / 15, 43, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349 /

  INTEGER :: dimid

  INCLUDE 'netcdf.inc'

  INTERFACE
    SUBROUTINE get_static_2d(ncid,varname,var2d)
      INTEGER,          INTENT(IN)  :: ncid
      CHARACTER(LEN=3), INTENT(IN)  :: varname
      REAL,             INTENT(OUT) :: var2d(:,:)
    END SUBROUTINE get_static_2d
  END INTERFACE
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  !
  ! This block is actually not necessary
  !
  IF (staticopt == 0) THEN        ! WRFSI data
    ! get the dimensions of data arrays
    istatus = NF_INQ_DIMID(ncid,'x',dimid)
    CALL nf_handle_error(istatus,'NF_INQ_DIMID')
    istatus = NF_INQ_DIMLEN(ncid,dimid,nx)
    CALL nf_handle_error(istatus,'NF_INQ_DIMLEN')
    istatus = NF_INQ_DIMID(ncid,'y',dimid)
    CALL nf_handle_error(istatus,'NF_INQ_DIMID')
    istatus = NF_INQ_DIMLEN(ncid,dimid,ny)
    CALL nf_handle_error(istatus,'NF_INQ_DIMLEN')

  ELSE  if (staticopt == 1) THEN  ! WRF static data
    istatus = NF_INQ_DIMID(ncid,'west_east_stag',dimid)
    CALL nf_handle_error(istatus,'NF_INQ_DIMID')
    istatus = NF_INQ_DIMLEN(ncid,dimid,nx)
    CALL nf_handle_error(istatus,'NF_INQ_DIMLEN')
    istatus = NF_INQ_DIMID(ncid,'south_north_stag',dimid)
    CALL nf_handle_error(istatus,'NF_INQ_DIMID')
    istatus = NF_INQ_DIMLEN(ncid,dimid,ny)
    CALL nf_handle_error(istatus,'NF_INQ_DIMLEN')
  
  ELSE
    WRITE(6,*) 'Unknown staticopt in get_static_monthly.'
    istatus = -1
    RETURN
  END IF

  ! Check data type character to make sure it is either greenness or albedo.
  IF ((vartype /= "a").AND.(vartype /= "g")) THEN
    PRINT *, 'Unknown data type character passed into get_static_monthly:', &
              vartype
    PRINT *,'Current supported values are a (albedo) and g (greenness fraction).'
    STOP 'get_static_monthly'
  ELSE IF (vartype == 'a') THEN
    PRINT *, ' Getting time-interpolated albedo.'
  ELSE IF (vartype == 'g') THEN
    PRINT *, ' Getting time-interpolated greenness.'
  END IF

  !PRINT *, 'Time-interpolating to day # ', valid_day
  ! Find bounding months
  IF ((valid_day < midmonth_day(1)) .OR. (valid_day > midmonth_day(12))) THEN
      ! December and January are bounding months
    d1 = midmonth_day(12)
    d2 = midmonth_day(1)
    m1 = 12
    m2 = 1
  ELSE
    find_bounds: DO m = 1, 11
      d1 = midmonth_day(m)
      d2 = midmonth_day(m+1)
      IF (valid_day == d1) THEN
         d2 = d1
         m1 = m
         m2 = m1
         EXIT find_bounds
      ELSE IF (valid_day == d2) THEN
         d1 = d2
         m1 = m + 1
         m2 = m1
         EXIT find_bounds
      ELSE IF ((valid_day > d1).AND.(valid_day < d2)) THEN
         m1 = m
         m2 = m + 1
         EXIT find_bounds
      ENDIF
    END DO find_bounds
  END IF 
  
  ! If d1 = d2, then we don't need any interpolation, just get that month's
  ! data values
  IF ( d1 == d2) THEN
    IF (staticopt == 0) THEN
      WRITE(varname, '(A1,I2.2)') vartype,m1
      CALL get_static_2d(ncid,varname,var2d)
    ELSE                ! staticopt == 1
      DO j = 1,ny
        DO i = 1,nx 
          var2d(i,j) = varin3d(i,j,m1)
        END DO
      END DO
    END IF
  ELSE

    ALLOCATE(data1 (nx,ny))
    ALLOCATE(data2 (nx,ny))

    ! We need to get the two months of bounding data and time interpolate
    IF (staticopt == 0) THEN  ! WRFSI data
      WRITE(varname, '(A1,I2.2)') vartype,m1
      CALL get_static_2d(ncid,varname,data1)
      WRITE(varname, '(A1,I2.2)') vartype,m2
      CALL get_static_2d(ncid,varname,data2)
    ELSE                      ! WRF static data
      data1(:,:) = varin3d(:,:,m1)
      data2(:,:) = varin3d(:,:,m2)
    END IF

    ! Compute weights
    IF (d2 > d1) THEN
      w1 = ( FLOAT(d2) - FLOAT(valid_day) ) / FLOAT(d2-d1)
    ELSE ! We must be between Dec 15 and Jan 15
      IF (valid_day < midmonth_day(1)) THEN ! We are in January
        w1 = ( FLOAT(d2) - FLOAT(valid_day) ) / 31.
      ELSE ! We are in December
        w1 = ( 366. - FLOAT(valid_day) + FLOAT(midmonth_day(1)) ) / 31.
      END IF

    END IF
    w2 = 1. - w1 
    DO j = 1,ny
      DO i = 1,nx
        var2d(i,j) = w1*data1(i,j) + w2*data2(i,j)
      END DO
    END DO
!    WRITE(6,*) '**** ',w1,w2,data1(1,1),data2(1,1),var2d(1,1)
!    StOP

    DEALLOCATE(data1)
    DEALLOCATE(data2)
  END IF 
    
  istatus = 0
  !WRITE(6,*) 'Returning from get_static_monthly'
    
  RETURN 
END SUBROUTINE get_static_monthly
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_static_2d                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_static_2d(ncid,varname, var2d)

! Subroutine to get the 2D field from the WRFSI static
! file
  
  USE wrf_metadata
  
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: ncid
  CHARACTER(LEN=3), INTENT(IN)  :: varname
  REAL,             INTENT(OUT) :: var2d(:,:)
  
  INCLUDE 'netcdf.inc'
!------------------------------------------------------------------
!  
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus,vid
  REAL    :: fillValue
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = NF_INQ_VARID(ncid, varname, vid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID')

  istatus = NF_GET_VAR_REAL(ncid,vid,var2d)
  CALL nf_handle_error(istatus,'NF_GET_VAR_REAL')

!  istatus = NF_GET_ATT_REAL(ncid,vid,'_FillVale',fillValue)
!  WHERE(var2d >= fillValue) var2d = rmissing

  RETURN
END SUBROUTINE get_static_2d

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE writebdy                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE writebdy(nout,io_form,DateStr,itime,varname,stagger,desc,   &
                    bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,             &
                    fzone,dbdyw,dbdye,dbdys,dbdyn,                     &
                    nxlg,nylg,nzlg,tem1,tem2)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write the 4 lateral boudnary arrays
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,          INTENT(IN) :: nout   
  INTEGER,          INTENT(IN) :: io_form   
  CHARACTER(LEN=*), INTENT(IN) :: DateStr
  INTEGER,          INTENT(IN) :: itime   
  CHARACTER(LEN=*), INTENT(IN) :: varname
  CHARACTER(LEN=*), INTENT(IN) :: stagger
  CHARACTER(LEN=*), INTENT(IN) :: desc
  INTEGER,          INTENT(IN) :: nx,ny,nz,bdywidth
  INTEGER,          INTENT(IN) :: fzone   
  REAL,             INTENT(INOUT) :: bdys(nx,nz,bdywidth)
  REAL,             INTENT(INOUT) :: bdyn(nx,nz,bdywidth)
  REAL,             INTENT(INOUT) :: bdyw(ny,nz,bdywidth)
  REAL,             INTENT(INOUT) :: bdye(ny,nz,bdywidth)
  INTEGER,          INTENT(IN) :: nxlg,nylg,nzlg
  REAL,             INTENT(IN) :: dbdys(nxlg,nzlg,bdywidth)
  REAL,             INTENT(IN) :: dbdyn(nxlg,nzlg,bdywidth)
  REAL,             INTENT(IN) :: dbdyw(nylg,nzlg,bdywidth)
  REAL,             INTENT(IN) :: dbdye(nylg,nzlg,bdywidth)

  REAL,             INTENT(OUT) :: tem1(nx,nz,bdywidth)
  REAL,             INTENT(OUT) :: tem2(ny,nz,bdywidth)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus

  INCLUDE 'mp.inc'
 
  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
  
  INTEGER       :: ilocs,iloce,jlocs,jloce

  LOGICAL       :: IOFlag

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
    
  DimNames(2) = 'bottom_top'
  DimNames(3) = 'bdy_width'

  DomainStart(2) = 1
  MemoryStart(2) = 1
  PatchStart(2)  = 1
  DomainEnd(2)   = nzlg
  MemoryEnd(2)   = nz
  PatchEnd(2)    = DomainEnd(2)

  DomainStart(3) = 1
  MemoryStart(3) = 1
  PatchStart(3)  = 1
  DomainEnd(3)   = bdywidth
  MemoryEnd(3)   = DomainEnd(3)
  PatchEnd(3)    = DomainEnd(3)

  IF (io_form == 5) THEN

    ilocs = (nx-fzone)*(loc_x-1)+fzone
    jlocs = (ny-fzone)*(loc_y-1)+fzone
    iloce = (nx-fzone)*(loc_x)+fzone
    jloce = (ny-fzone)*(loc_y)+fzone

    !
    ! West boudnary
    !
    DomainStart(1) = 1
    DomainEnd(1)   = nylg
    MemoryStart(1) = jlocs
    MemoryEnd(1)   = jloce
    PatchStart(1)  = jlocs
    PatchEnd(1)    = jloce - fzone
    IF (stagger == 'Y') THEN
      IF (loc_y > 1 ) PatchStart(1) = jlocs + fzone
      PatchEnd(1)   =  jloce
      DimNames(1) = 'south_north_stag'
    ELSE
      DimNames(1) = 'south_north'
    END IF

    IOFLAG = .FALSE.
    IF (loc_x == 1) IOFLAG = .TRUE.

!    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing West boundary for ', varname
     
    CALL write_phdf5_bdy(nout,DateStr,varname//'XS',desc,'',stagger,  &
                           bdyw, WRF_REAL, 1,'XSZ',DimNames,            &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, IOFlag, istatus)
!    IF (istatus == 0) THEN
!      WRITE(6,'(a)') '         ...  DONE.'
!    ELSE
!      WRITE(6,'(a,I3,a)') '         ...  ERROR, ',istatus,'.'
!    END IF

    !
    ! East boudnary
    !

    IOFLAG = .FALSE.
    IF (loc_x == nproc_x) IOFLAG = .TRUE.

!    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing East boundary for ', varname

    CALL write_phdf5_bdy(nout,DateStr,varname//'XE',desc,'',stagger,  &
                           bdye, WRF_REAL, 1,                           &
                           'XEZ',DimNames,                              &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, IOFlag, istatus)

!    IF (istatus == 0) THEN
!      WRITE(6,'(a)') '         ...  DONE.'
!    ELSE
!      WRITE(6,'(a,I3,a)') '         ...  ERROR, ',istatus,'.'
!    END IF

    !
    ! South boudnary
    !
    DomainStart(1) = 1
    DomainEnd(1)   = nxlg
    MemoryStart(1) = ilocs
    MemoryEnd(1)   = iloce
    PatchStart(1)  = ilocs
    PatchEnd(1)    = iloce - fzone
    IF (stagger == 'X') THEN
      IF (loc_x > 1 ) PatchStart(1) = ilocs + fzone
      PatchEnd(1)   =  iloce
      DimNames(1)   = 'west_east_stag'
    ELSE
      DimNames(1)   = 'west_east'
    END IF

    IOFLAG = .FALSE.
    IF (loc_y == 1) IOFLAG = .TRUE.

!    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing South boundary for ', varname

    CALL write_phdf5_bdy(nout,DateStr,varname//'YS',desc,'',stagger,  &
                           bdys, WRF_REAL, 1,                           &
                           'YSZ',DimNames,                              &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, IOFlag, istatus)

!    IF (istatus == 0) THEN
!      WRITE(6,'(a)') '         ...  DONE.'
!    ELSE
!      WRITE(6,'(a,I3,a)') '         ...  ERROR, ',istatus,'.'
!    END IF

    !
    ! North boudnary
    !

    IOFLAG = .FALSE.
    IF (loc_y == nproc_y) IOFLAG = .TRUE.

!    WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing North boundary for ', varname

    CALL write_phdf5_bdy(nout,DateStr,varname//'YE',desc,'',stagger,  &
                           bdyn, WRF_REAL, 1,                           &
                           'YEZ',DimNames,                              &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, IOFlag, istatus)

!    IF (istatus == 0) THEN
!      WRITE(6,'(a)') '         ...  DONE.'
!    ELSE
!      WRITE(6,'(a,I3,a)') '         ...  ERROR, ',istatus,'.'
!    END IF

  ELSE 

    IF (stagger == 'X') THEN
      CALL wrf_mergebdyu(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,         &    
               fzone,dbdyw,dbdye,dbdys,dbdyn,tem1,tem2)
    ELSE IF (stagger == 'Y') THEN
      CALL wrf_mergebdyv(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,         &    
               fzone,dbdyw,dbdye,dbdys,dbdyn,tem1,tem2)
    ELSE IF (stagger == 'Z') THEN
      CALL wrf_mergebdyw(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,         &    
               fzone,dbdyw,dbdye,dbdys,dbdyn,tem1,tem2)
    ELSE
      CALL wrf_mergebdyt(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,         &    
               fzone,dbdyw,dbdye,dbdys,dbdyn,tem1,tem2)
    END IF

    IF (io_form == 1) THEN          ! Binary

      IF (myproc == 0) THEN
        !
        ! West boudnary
        !
        DomainStart(1) = 1
        DomainEnd(1)   = nylg
        IF (stagger == 'Y') THEN
          DimNames(1)   = 'south_north_stag'
        ELSE
          DimNames(1)   = 'south_north'
        END IF

        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing West boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'XS',          &
                 dbdyw, WRF_REAL, 0, 0, 1, 'XSZ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

        !
        ! East boudnary
        !
        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing East boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'XE',          &
                 dbdye, WRF_REAL, 0, 0, 1, 'XEZ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

        !
        ! South boudnary
        !
        IF (stagger == 'X') THEN
          DimNames(1)   = 'west_east_stag'
        ELSE
          DimNames(1)   = 'west_east'
        END IF
        DomainStart(1) = 1
        DomainEnd(1)   = nxlg

        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing South boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'YS',          &
                 dbdys, WRF_REAL, 0, 0, 1, 'YSZ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

        !
        ! North boudnary
        !
        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing North boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'YE',          &
                 dbdyn, WRF_REAL, 0, 0, 1, 'YEZ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

      END IF

    ELSE IF (io_form == 7) THEN     ! NetCDF

      IF (myproc == 0)   & 
        CALL write_ncd_bdy(nout,nxlg,nylg,nzlg, bdywidth,itime,varname, &
                    dbdys,dbdyn,dbdyw,dbdye,istatus)
    END IF

  END IF

  RETURN
END SUBROUTINE writebdy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE writebdy2d               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE writebdy2d(nout,io_form,DateStr,itime,varname,stagger,desc, &
                    bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,             &
                    fzone,dbdyw,dbdye,dbdys,dbdyn,                     &
                    nxlg,nylg,tem1,tem2)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write the 4 lateral boudnary arrays
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER,          INTENT(IN) :: nout   
  INTEGER,          INTENT(IN) :: io_form   
  CHARACTER(LEN=*), INTENT(IN) :: DateStr
  INTEGER,          INTENT(IN) :: itime   
  CHARACTER(LEN=*), INTENT(IN) :: varname
  CHARACTER(LEN=*), INTENT(IN) :: stagger
  CHARACTER(LEN=*), INTENT(IN) :: desc
  INTEGER,          INTENT(IN) :: nx,ny,nz,bdywidth
  INTEGER,          INTENT(IN) :: fzone   
  REAL,             INTENT(IN) :: bdys(nx,nz,bdywidth)
  REAL,             INTENT(IN) :: bdyn(nx,nz,bdywidth)
  REAL,             INTENT(IN) :: bdyw(ny,nz,bdywidth)
  REAL,             INTENT(IN) :: bdye(ny,nz,bdywidth)
  INTEGER,          INTENT(IN) :: nxlg,nylg
  REAL,             INTENT(IN) :: dbdys(nxlg,bdywidth)
  REAL,             INTENT(IN) :: dbdyn(nxlg,bdywidth)
  REAL,             INTENT(IN) :: dbdyw(nylg,bdywidth)
  REAL,             INTENT(IN) :: dbdye(nylg,bdywidth)

  REAL,             INTENT(OUT) :: tem1(nx,bdywidth)
  REAL,             INTENT(OUT) :: tem2(ny,bdywidth)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: istatus

  INCLUDE 'mp.inc'
 
  CHARACTER(80) :: DimNames(3)
  INTEGER       :: DomainStart(3), DomainEnd(3)
  INTEGER       :: MemoryStart(3), MemoryEnd(3)
  INTEGER       :: PatchStart(3),  PatchEnd(3)
  
  INTEGER       :: ilocs,iloce,jlocs,jloce

  LOGICAL       :: IOFLAG
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DimNames(2) = 'bdy_width'
  DimNames(3) = ''

  DomainStart(2) = 1
  DomainEnd(2)   = bdywidth
  DomainStart(3) = 1
  DomainEnd(3)   = 1

  IF (io_form == 5) THEN

    MemoryStart(2) = 1
    MemoryEnd(2)   = DomainEnd(2)
    MemoryStart(3) = 1
    MemoryEnd(3)   = 1

    PatchStart(2)  = 1
    PatchEnd(2)    = DomainEnd(2)
    PatchStart(3)  = 1
    PatchEnd(3)    = 1

    ilocs = (nx-fzone)*(loc_x-1)+fzone
    jlocs = (ny-fzone)*(loc_y-1)+fzone
    iloce = (nx-fzone)*(loc_x)+fzone
    jloce = (ny-fzone)*(loc_y)+fzone

    !
    ! West boudnary
    !
    DomainStart(1) = 1
    DomainEnd(1)   = nylg
    MemoryStart(1) = jlocs
    MemoryEnd(1)   = jloce
    PatchStart(1)  = jlocs
    PatchEnd(1)    = jloce - fzone
    IF (stagger == 'Y') THEN
      IF (loc_y > 1 ) PatchStart(1) = jlocs + fzone
      PatchEnd(1)   =  jloce
      DimNames(1) = 'south_north_stag'
    ELSE
      DimNames(1) = 'south_north'
    END IF

    IOFLAG = .FALSE.
    IF (loc_x == 1) IOFLAG = .TRUE.

    tem2(:,:) = bdyw(:,1,:)
    CALL write_phdf5_bdy(nout,DateStr,varname//'XS',desc,'',stagger,  &
                           tem2, WRF_REAL, 1,'XS',DimNames,             &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd,                         &
                           IOFLAG,istatus)

    !
    ! East boudnary
    !
    IOFLAG = .FALSE.
    IF (loc_x == nproc_x) IOFLAG = .TRUE.

    tem2(:,:) = bdye(:,1,:)
    CALL write_phdf5_bdy(nout,DateStr,varname//'XE',desc,'',stagger,  &
                           tem2, WRF_REAL, 1,'XE',DimNames,             &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd,                         &
                           IOFLAG,istatus)

    !
    ! South boudnary
    !
    DomainStart(1) = 1
    DomainEnd(1)   = nxlg
    MemoryStart(1) = ilocs
    MemoryEnd(1)   = iloce
    PatchStart(1)  = ilocs
    PatchEnd(1)    = iloce - fzone
    IF (stagger == 'X') THEN
      IF (loc_y > 1 ) PatchStart(1) = ilocs + fzone
      PatchEnd(1)   =  iloce
      DimNames(1)   = 'west_east_stag'
    ELSE
      DimNames(1)   = 'west_east'
    END IF

    IOFLAG = .FALSE.
    IF (loc_y == 1) IOFLAG = .TRUE.

    tem1(:,:) = bdys(:,1,:)
    CALL write_phdf5_bdy(nout,DateStr,varname//'YS',desc,'',stagger,  &
                           tem1, WRF_REAL, 1,'YS',DimNames,             &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, IOFLAG,istatus)
    !
    ! North boudnary
    !
    IOFLAG = .FALSE.
    IF (loc_y == nproc_y) IOFLAG = .TRUE.

    tem1(:,:) = bdyn(:,1,:)
    CALL write_phdf5_bdy(nout,DateStr,varname//'YE',desc,'',stagger,  &
                           tem1, WRF_REAL, 1,'YE',DimNames,             &
                           DomainStart,DomainEnd,                       &
                           MemoryStart,MemoryEnd,                       &
                           PatchStart,PatchEnd, IOFLAG, istatus)

  ELSE

    IF (stagger == 'X') THEN
      WRITE(0,*) 'WARNING: To be implemented, wrf_mergebdy2du.'
    ELSE IF (stagger == 'Y') THEN
      WRITE(0,*) 'WARNING: To be implemented, wrf_mergebdy2dv.'
    ELSE IF (stagger == 'Z') THEN
      WRITE(0,*) 'WARNING: To be implemented, wrf_mergebdy2dw.'
    ELSE
      CALL wrf_mergebdy2d(bdyw,bdye,bdys,bdyn,nx,ny,nz,bdywidth,        &
               fzone,dbdyw,dbdye,dbdys,dbdyn,tem1,tem2)
    END IF

    IF (io_form == 1) THEN               ! binary

      IF (myproc == 0) THEN
        !
        ! West boudnary
        !
        DomainStart(1) = 1
        DomainEnd(1)   = nylg
        IF (stagger == 'Y') THEN
          DimNames(1)   = 'south_north_stag'
        ELSE
          DimNames(1)   = 'south_north'
        END IF

        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing West boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'XS',          &
                 dbdyw, WRF_REAL, 0, 0, 1, 'XS ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

        !
        ! East boudnary
        !
        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing East boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'XE',          &
                 dbdye, WRF_REAL, 0, 0, 1, 'XE ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

        !
        ! South boudnary
        !
        IF (stagger == 'X') THEN
          DimNames(1)   = 'west_east_stag'
        ELSE
          DimNames(1)   = 'west_east'
        END IF
        DomainStart(1) = 1
        DomainEnd(1)   = nxlg

        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing South boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'YS',          &
                 dbdys, WRF_REAL, 0, 0, 1, 'YS ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

        !
        ! North boudnary
        !
        WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing North boundary for ', varname
     
        CALL ext_int_write_field( nout, DateStr,varname//'YE',          &
                 dbdyn, WRF_REAL, 0, 0, 1, 'YE ', Stagger, DimNames,    &
                 DomainStart, DomainEnd, DomainStart, DomainEnd,        &
                 DomainStart, DomainEnd, iStatus )

        IF (istatus == 0) WRITE(6,'(a)') '         ...  DONE.'

      END IF

    ELSE IF (io_form == 7) THEN          ! NetCDF
      IF (myproc == 0)  &
        CALL write_ncd_bdy2d(nout,nxlg,nylg,bdywidth,itime,varname,     &
                    dbdys,dbdyn,dbdyw,dbdye,istatus)

    END IF

  END IF

  RETURN
END SUBROUTINE writebdy2d

!
!##################################################################
!##################################################################
!######                                                      ######
!######     SUBROUTINE write_new_gravity_drag_static         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE write_new_gravity_drag_static(nout,io_form,DateStr,         &
                  nxsize,nysize,nzsize,sfcinitopt,sfcdtfl,tem1,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Read new WRF static fields from GEOGRID and writ to wrfinput_d01
!    for new gravity drag option with WRFV3.1.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER,          INTENT(IN) :: nout
  INTEGER,          INTENT(IN) :: io_form
  CHARACTER(LEN=*), INTENT(IN) :: DateStr
  INTEGER,          INTENT(IN) :: nxsize, nysize, nzsize
  CHARACTER(LEN=*), INTENT(IN) :: sfcinitopt
  CHARACTER(LEN=*), INTENT(IN) :: sfcdtfl

  REAL,             INTENT(OUT) :: tem1(nxsize,nysize,nzsize)
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INCLUDE 'mp.inc'

  INTEGER :: nscid

  INTEGER :: n
  INTEGER, PARAMETER :: nsize = 10
  CHARACTER(LEN=3)   :: varnames(nsize) = (/'VAR','CON',                &
                                            'OA1','OA2','OA3','OA4',    &
                                            'OL1','OL2','OL3','OL4' /)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (io_form /= 7) THEN
    WRITE(6,'(/,1x,a,/)') 'WARNING: File format other than netCDF to be implemented '//  &
                          'for new fields with gravity drag option since WRFV3.1.'
    RETURN
  END IF

  IF (myproc == 0) THEN   ! Copy field from static file by GEOGRID

    CALL open_static_file(TRIM(sfcdtfl),nscid)

    WRITE(6,'(2x,a)') 'Copying 3D variable LANDUSEF'
    CALL get_ncd_3d(nscid,1,'LANDUSEF',nxsize,nysize,LanduseCategories,tem1,istatus)
    CALL write_ncd_3d(nout, 'LANDUSEF',tem1,nxsize,nysize,LanduseCategories,istatus)

    WRITE(6,'(2x,a)') 'Copying 3D variable SOILCTOP'
    CALL get_ncd_3d(nscid,1,'SOILCTOP',nxsize,nysize,SoilCategories,tem1,istatus)
    CALL write_ncd_3d(nout, 'SOILCTOP',tem1,nxsize,nysize,SoilCategories,istatus)

    WRITE(6,'(2x,a)') 'Copying 3D variable SOILCBOT'
    CALL get_ncd_3d(nscid,1,'SOILCBOT',nxsize,nysize,SoilCategories,tem1,istatus)
    CALL write_ncd_3d(nout, 'SOILCBOT',tem1,nxsize,nysize,SoilCategories,istatus)

    DO n = 1, nsize
      WRITE(6,'(2x,2a)') 'Copying 2D variable ',varnames(n)
      CALL get_ncd_2d(nscid,1,varnames(n),nxsize,nysize,tem1,istatus)
      CALL write_ncd_2d(nout, varnames(n),tem1,nxsize,nysize,istatus)
    END DO

    CALL close_static_file(nscid)
  END IF

  RETURN
END SUBROUTINE write_new_gravity_drag_static
