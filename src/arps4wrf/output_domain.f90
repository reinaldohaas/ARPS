
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: output_init
!
! Purpose: To initialize the output module. Such initialization may include
!   opening an X window, and making initialization calls to an I/O API.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE output_init( datestr, ndomain, dbglvl,                       &
                        IAMIO_NODE, do_tiled_output,io_form_output,     &
                        output_path,output_fname, numdigits,fhandle, istatus)

  USE module_metgrid
  USE module_geogrid
  USE module_commontypes

  IMPLICIT NONE

  ! Arguments
  INTEGER,            INTENT(IN) :: ndomain, dbglvl
  LOGICAL,            INTENT(IN) :: IAMIO_NODE, do_tiled_output
  CHARACTER(LEN=19),  INTENT(IN) :: datestr
  INTEGER,            INTENT(IN) :: io_form_output

  CHARACTER(LEN=MAXFILELEN),  INTENT(IN)  :: output_path
  CHARACTER(LEN=MAXFILELEN),  INTENT(OUT) :: output_fname
  INTEGER,                    INTENT(IN)  :: numdigits
  INTEGER,                    INTENT(OUT) :: fhandle
  INTEGER,                    INTENT(OUT) :: istatus


  INTERFACE
    SUBROUTINE write_field( fhandle, io_form_output, datestr,           &
                            IAMIO_NODE, do_tiled_output, dbglvl,        &
                            field, is_training, istatus )
      USE module_commontypes
      IMPLICIT NONE

      ! Arguments
      INTEGER,           INTENT(IN)  :: fhandle
      INTEGER,           INTENT(IN)  :: io_form_output
      CHARACTER(LEN=19), INTENT(IN)  :: datestr
      LOGICAL,           INTENT(IN)  :: IAMIO_NODE, do_tiled_output  ! false for serial
      INTEGER,           INTENT(IN)  :: dbglvl
      TYPE( TYPE_FIELD), INTENT(IN)  :: field
      LOGICAL,           INTENT(IN)  :: is_training
      INTEGER,           INTENT(OUT) :: istatus
    END SUBROUTINE write_field
  END INTERFACE


  ! Local variables
  INTEGER :: i, comm_1, comm_2
  LOGICAL :: supports_training
  CHARACTER (LEN=128) :: coption

  CHARACTER(LEN=20)   :: fmtstr
  TYPE(TYPE_FIELD) :: field

  INTEGER :: ntotal

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (IAMIO_NODE) THEN
    istatus = 0
    IF (io_form_output == BINARY) CALL ext_int_ioinit('sysdep info', istatus)
    IF (io_form_output == NETCDF) CALL ext_ncd_ioinit('sysdep info', istatus)

!      CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_ioinit')

    ! Find out what this implementation of WRF I/O API supports
    istatus = 0
    IF (io_form_output == BINARY) coption = 'REQUIRE'
    IF (io_form_output == NETCDF) CALL ext_ncd_inquiry('OPEN_COMMIT_WRITE',coption,istatus)
!      CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_inquiry')

    IF (index(coption,'ALLOW') /= 0) THEN
       supports_training = .false.
    ELSE IF (index(coption,'REQUIRE') /= 0) THEN
       supports_training = .true.
    ELSE IF (index(coption,'NO') /= 0) THEN
       supports_training = .false.
    END IF

!    istatus = 0
!    IF (io_form_output == BINARY) coption = 'YES'
!    IF (io_form_output == NETCDF) CALL ext_ncd_inquiry('SUPPORT_3D_FIELDS',coption,istatus)
!!      CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_inquiry')
!
!    IF (index(coption,'YES') /= 0) THEN
!       supports_3d_fields = .true.
!    ELSE IF (index(coption,'NO') /= 0) THEN
!       supports_3d_fields = .false.
!!         CALL mprintf(.true.,ERROR,'WRF I/O API implementation does NOT support 3-d fields.')
!    END IF

    comm_1 = 1
    comm_2 = 1

    output_fname = ' '
    IF (geogrid%grid_type == 'C') THEN
      IF (io_form_output == BINARY) THEN
         output_fname = trim(output_path)//'met_em.d  .'//trim(datestr)//'.int'
      END IF
      IF (io_form_output == NETCDF) THEN
         output_fname = trim(output_path)//'met_em.d  .'//trim(datestr)//'.nc'
      END IF
      i = LEN_TRIM(output_path)
      WRITE(output_fname(i+9:i+10),'(i2.2)') ndomain
    ELSE IF (geogrid%grid_type == 'E') THEN
      IF (io_form_output == BINARY) THEN
         output_fname = trim(output_path)//'met_nmm.d  .'//trim(datestr)//'.int'
      END IF
      IF (io_form_output == NETCDF) THEN
         output_fname = trim(output_path)//'met_nmm.d  .'//trim(datestr)//'.nc'
      END IF
      i = len_trim(output_path)
      write(output_fname(i+10:i+11),'(i2.2)') ndomain
    END IF

    IF (ntprocs > 1 .AND. do_tiled_output) THEN
      i = len_trim(output_fname) + 1
      WRITE(fmtstr,'(a,I1,a,I1,a)') '(a1,i',numdigits,'.',numdigits,')'
      write(output_fname(i:i+numdigits), FMT=fmtstr) '_', my_proc_id
    END IF

    IF (dbglvl > 3) WRITE(6,'(5x,2a)') '^^^ Output file name is ',TRIM(output_fname)
  END IF

  CALL parallel_bcast_LOGICAL(supports_training)

  ! IF the implementation requires or supports open_for_write begin/commit semantics
  IF (supports_training) THEN

    IF (dbglvl > 3) WRITE(6,'(5x,a)') '^^^ Training file for output ... '

    IF ( IAMIO_NODE ) THEN
      istatus = 0
      IF (io_form_output == BINARY) THEN
        CALL ext_int_open_for_write_begin(trim(output_fname),           &
                         comm_1, comm_2, 'sysdep info', fhandle, istatus)
      ELSE IF (io_form_output == NETCDF) THEN
        CALL ext_ncd_open_for_write_begin(trim(output_fname),           &
                         comm_1, comm_2, 'sysdep info', fhandle, istatus)
      END IF
!            CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_open_for_write_begin.')
    END IF

    DO i=1,metgrid%nfields
      field = metgrid%fields(i)

      IF (dbglvl > 3) WRITE(6,'(5x,a,I2.2,2a)')              &
            '^^^ Training for field (',i,') ',TRIM(field%fieldname)

      CALL write_field( fhandle, io_form_output, datestr,    &
                        IAMIO_NODE, do_tiled_output, dbglvl, &
                        field, .TRUE., istatus )
    END do

    ntotal = geogrid%nfields + metgrid%nfields

    DO i=geogrid%nfields,1,-1
      field = geogrid%fields(i)

      IF (dbglvl > 3) WRITE(6,'(5x,a,I2.2,2a)')              &
            '^^^ Training for field (',ntotal+i-1,') ',TRIM(field%fieldname)

      CALL write_field( fhandle, io_form_output, datestr,    &
                        IAMIO_NODE, do_tiled_output, dbglvl, &
                        field, .TRUE., istatus )
    END do

    IF ( IAMIO_NODE ) THEN
      istatus = 0
      IF (io_form_output == BINARY) THEN
        CALL ext_int_open_for_write_commit(fhandle, istatus)
      ELSE IF (io_form_output == NETCDF) THEN
        CALL ext_ncd_open_for_write_commit(fhandle, istatus)
      END IF
!            CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_write_commit')
    END IF

  ELSE ! No training required

    IF ( IAMIO_NODE ) THEN
      IF (dbglvl > 3) WRITE(6,'(5x,a)') '^^^ Openning for output ...'

      istatus = 0
      IF (io_form_output == BINARY) THEN
          CALL ext_int_open_for_write(trim(output_fname),               &
                         comm_1, comm_2, 'sysdep info', fhandle, istatus)
      ELSE IF (io_form_output == NETCDF) THEN
          CALL ext_ncd_open_for_write(trim(output_fname),               &
                         comm_1, comm_2, 'sysdep info', fhandle, istatus)
      END IF
!            CALL mprintf((istatus /= 0),ERROR,'Error in ext_pkg_open_for_write_begin')
    END IF

  END IF

  RETURN
END SUBROUTINE output_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: write_field
!
! Purpose: This routine writes the provided field to any output devices or APIs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_field( fhandle, io_form_output, datestr,               &
                        IAMIO_NODE, do_tiled_output, dbglvl,            &
                        field, is_training, istatus )

  USE module_geogrid
  USE wrf_parallel_module

  IMPLICIT NONE

  ! Arguments
  INTEGER,           INTENT(IN)  :: fhandle
  INTEGER,           INTENT(IN)  :: io_form_output
  CHARACTER(LEN=19), INTENT(IN)  :: datestr
  LOGICAL,           INTENT(IN)  :: IAMIO_NODE, do_tiled_output  ! false for serial
  INTEGER,           INTENT(IN)  :: dbglvl
  TYPE( TYPE_FIELD), INTENT(IN)  :: field
  LOGICAL,           INTENT(IN)  :: is_training
  INTEGER,           INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  CHARACTER(LEN=128) :: cstagger
  INTEGER, DIMENSION(MAXDIMENSIONS) :: sp,ep,sm,em,sd,ed

  ! Local variables
  INTEGER :: i
  INTEGER :: comm_1, comm_2, domain_desc

  REAL, POINTER :: real_dom_array(:,:,:)
  LOGICAL       :: allocated_real_locally

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  allocated_real_locally = .FALSE.
  istatus = 0

  ! IF we are running distributed memory and need to gather all tiles
  ! onto a single processor for output
  IF ( .NOT. is_training .AND. (ntprocs > 1 .AND. .NOT. do_tiled_output) ) THEN

    ! For the gather routines below, the IO_NODE should give the full
    ! domain dimensions, but the  memory and patch dimensions should
    ! indicate what the processor already has in its patch_array.
    ! This is because an array with dimensions of the full domain
    ! will be allocated, and the patch_array will be copied from local
    ! memory into the full domain array in the area specified by the
    ! patch dimensions.
    sd = field%dom_start
    ed = field%dom_end
    sp = field%patch_start
    ep = field%patch_end
    sm = field%mem_start
    em = field%mem_end

    IF (dbglvl > 4) WRITE(6,'(7x,a)') '--- Gathering global array ...'

    ALLOCATE(real_dom_array(sd(1):ed(1),sd(2):ed(2),sd(3):ed(3)), STAT = istatus)
    allocated_real_locally = .true.
    CALL gather_whole_field_r(field%rdata_arr,                          &
                              sm(1), em(1), sm(2), em(2), sm(3), em(3), &
                              sp(1), ep(1), sp(2), ep(2), sp(3), ep(3), &
                              real_dom_array,                           &
                              sd(1), ed(1), sd(2), ed(2), sd(3), ed(3))
  ELSE
      real_dom_array => field%rdata_arr
  END IF

  ! Now a write CALL is only done if we are the IO_NODE
  IF ( IAMIO_NODE ) THEN
    comm_1 = 1
    comm_2 = 1
    domain_desc = 0

    ! Here, the output array has dimensions of the full grid IF it was
    ! gathered together from all processors
    IF ( ntprocs > 1 .AND. .NOT. do_tiled_output) THEN
      sd = field%dom_start
      ed = field%dom_end
      sm = sd
      em = ed
      sp = sd
      ep = ed
    ELSE
      ! IF we are writing one file per processor, THEN each processor
      ! only writes out the part of the domain that it has in memory
! BUG: Shouldn't we set sd/ed to be domain_start/domain_end?
!      Maybe not, since patch is already adjusted for staggering;
!      but maybe so, and also adjust for staggering IF it is alright
!      to pass true domain dimensions to write_field.
      sd = field%patch_start
      ed = field%patch_end
      sp = field%patch_start
      ep = field%patch_end
      sm = field%mem_start
      em = field%mem_end
    END IF

    istatus = 0
    SELECT CASE (field%stagger)
    CASE (M, HH)
      cstagger = 'M'
    CASE (V, VV)
      cstagger = 'V'
    CASE (U)
      cstagger = 'U'
    END SELECT

    IF (dbglvl > 4) WRITE(6,'(7x,2a)')                                  &
               '--- Writing data for ',trim(field%fieldname)

    IF (io_form_output == BINARY) THEN
      CALL ext_int_write_field(fhandle, datestr, trim(field%fieldname), &
                    real_dom_array, WRF_REAL, comm_1, comm_2,           &
                    domain_desc, trim(field%mem_order), trim(cstagger), &
                    field%dimnames, sd, ed, sm, em, sp, ep, istatus)
    ELSE IF (io_form_output == NETCDF) THEN
      CALL ext_ncd_write_field(fhandle, datestr, trim(field%fieldname), &
                    real_dom_array, WRF_REAL, comm_1, comm_2,           &
                    domain_desc, trim(field%mem_order), trim(cstagger), &
                    field%dimnames, sd, ed, sm, em, sp, ep, istatus)
    END IF
!               CALL mprintf((istatus /= 0),ERROR,'Error in ext_pkg_write_field')

    IF (is_training) THEN
      IF (dbglvl > 4) WRITE(6,'(7x,2a)')                            &
          '--- Writing attributes for ',trim(field%fieldname)

      IF (io_form_output == BINARY) THEN
        CALL ext_int_put_var_ti_char(fhandle, 'units',              &
              trim(field%fieldname), trim(field%units), istatus)
        CALL ext_int_put_var_ti_char(fhandle, 'description',        &
              trim(field%fieldname), trim(field%descr), istatus)
        CALL ext_int_put_var_ti_char(fhandle, 'stagger',            &
              trim(field%fieldname), trim(cstagger), istatus)
        CALL ext_int_put_var_ti_integer(fhandle,'sr_x',             &
              trim(field%fieldname),(/1/),1, istatus)
        CALL ext_int_put_var_ti_integer(fhandle,'sr_y',             &
              trim(field%fieldname),(/1/),1, istatus)
      ELSE IF (io_form_output == NETCDF) THEN
        CALL ext_ncd_put_var_ti_char(fhandle, 'units',              &
              trim(field%fieldname), trim(field%units), istatus)
        CALL ext_ncd_put_var_ti_char(fhandle, 'description',        &
              trim(field%fieldname), trim(field%descr), istatus)
        CALL ext_ncd_put_var_ti_char(fhandle, 'stagger',            &
              trim(field%fieldname), trim(cstagger), istatus)
        CALL ext_ncd_put_var_ti_integer(fhandle,'sr_x',             &
              trim(field%fieldname),(/1/),1, istatus)
        CALL ext_ncd_put_var_ti_integer(fhandle,'sr_y',             &
              trim(field%fieldname),(/1/),1, istatus)
      END IF
    END IF

  END IF

  IF (allocated_real_locally) DEALLOCATE(real_dom_array)

  RETURN
END SUBROUTINE write_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: write_global_attrs
!
! Purpose:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write_global_attrs(fhandle,io_form_output,                   &
                              IAMIO_NODE, do_tiled_output, dbglvl,      &
                              title, start_date, nflags, flags, istatus)

  USE module_geogrid
  USE module_metgrid
  USE wrf_parallel_module

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)  :: fhandle, io_form_output
  INTEGER, INTENT(IN)  :: dbglvl
  LOGICAL, INTENT(IN)  :: IAMIO_NODE, do_tiled_output
  CHARACTER(LEN=*),  INTENT(IN)  :: title, start_date
  INTEGER,           INTENT(OUT) :: istatus

  INTEGER,           INTENT(IN) :: nflags
  CHARACTER(LEN=128),INTENT(IN) :: flags(nflags)


  ! Local variables
  INTEGER :: local_we_patch_s, local_we_patch_s_stag, &
             local_we_patch_e, local_we_patch_e_stag, &
             local_sn_patch_s, local_sn_patch_s_stag, &
             local_sn_patch_e, local_sn_patch_e_stag
  REAL, dimension(16) :: local_corner_lats, local_corner_lons
  INTEGER :: i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (dbglvl > 4) WRITE(6,'(7x,a)') '--- Entering write_global_attrs ...'

  local_we_patch_s      = geogrid%we_patch_s
  local_we_patch_s_stag = geogrid%we_patch_stag_s
  local_we_patch_e      = geogrid%we_patch_e
  local_we_patch_e_stag = geogrid%we_patch_stag_e
  local_sn_patch_s      = geogrid%sn_patch_s
  local_sn_patch_s_stag = geogrid%sn_patch_stag_s
  local_sn_patch_e      = geogrid%sn_patch_e
  local_sn_patch_e_stag = geogrid%sn_patch_stag_e
  local_corner_lats = geogrid%corner_lats
  local_corner_lons = geogrid%corner_lons

  IF (ntprocs > 1) THEN

     IF (dbglvl > 4) WRITE(6,'(7x,a)') '--- Collecting attributes ...'

     IF (.not. do_tiled_output) THEN
        CALL parallel_bcast_int(local_we_patch_s,      processors(0, 0))
        CALL parallel_bcast_int(local_we_patch_s_stag, processors(0, 0))
        CALL parallel_bcast_int(local_sn_patch_s,      processors(0, 0))
        CALL parallel_bcast_int(local_sn_patch_s_stag, processors(0, 0))

        CALL parallel_bcast_int(local_we_patch_e,      processors(nxprocs-1, nyprocs-1))
        CALL parallel_bcast_int(local_we_patch_e_stag, processors(nxprocs-1, nyprocs-1))
        CALL parallel_bcast_int(local_sn_patch_e,      processors(nxprocs-1, nyprocs-1))
        CALL parallel_bcast_int(local_sn_patch_e_stag, processors(nxprocs-1, nyprocs-1))
     END IF

     CALL parallel_bcast_real(local_corner_lats(1),  processors(0,         0))
     CALL parallel_bcast_real(local_corner_lats(2),  processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(3),  processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(4),  processors(nxprocs-1, 0))
     CALL parallel_bcast_real(local_corner_lats(5),  processors(0,         0))
     CALL parallel_bcast_real(local_corner_lats(6),  processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(7),  processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(8),  processors(nxprocs-1, 0))
     CALL parallel_bcast_real(local_corner_lats(9),  processors(0,         0))
     CALL parallel_bcast_real(local_corner_lats(10), processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(11), processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(12), processors(nxprocs-1, 0))
     CALL parallel_bcast_real(local_corner_lats(13), processors(0,         0))
     CALL parallel_bcast_real(local_corner_lats(14), processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(15), processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lats(16), processors(nxprocs-1, 0))

     CALL parallel_bcast_real(local_corner_lons(1),  processors(0,         0))
     CALL parallel_bcast_real(local_corner_lons(2),  processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(3),  processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(4),  processors(nxprocs-1, 0))
     CALL parallel_bcast_real(local_corner_lons(5),  processors(0,         0))
     CALL parallel_bcast_real(local_corner_lons(6),  processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(7),  processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(8),  processors(nxprocs-1, 0))
     CALL parallel_bcast_real(local_corner_lons(9),  processors(0,         0))
     CALL parallel_bcast_real(local_corner_lons(10), processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(11), processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(12), processors(nxprocs-1, 0))
     CALL parallel_bcast_real(local_corner_lons(13), processors(0,         0))
     CALL parallel_bcast_real(local_corner_lons(14), processors(0,         nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(15), processors(nxprocs-1, nyprocs-1))
     CALL parallel_bcast_real(local_corner_lons(16), processors(nxprocs-1, 0))
   END IF

   IF ( IAMIO_NODE ) THEN

     IF (dbglvl > 4) WRITE(6,'(7x,a)') '--- Writing attributes ...'

     CALL ext_put_dom_ti_char          ( fhandle, io_form_output,        &
                        'TITLE', title)
     CALL ext_put_dom_ti_char          ( fhandle, io_form_output,        &
                        'SIMULATION_START_DATE', start_date)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'WEST-EAST_GRID_DIMENSION',   geogrid%west_east_dim)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'SOUTH-NORTH_GRID_DIMENSION', geogrid%south_north_dim)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'BOTTOM-TOP_GRID_DIMENSION',  geogrid%bottom_top_dim)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'WEST-EAST_PATCH_START_UNSTAG',   local_we_patch_s)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'WEST-EAST_PATCH_END_UNSTAG',     local_we_patch_e)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'WEST-EAST_PATCH_START_STAG',     local_we_patch_s_stag)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'WEST-EAST_PATCH_END_STAG',       local_we_patch_e_stag)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'SOUTH-NORTH_PATCH_START_UNSTAG', local_sn_patch_s)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'SOUTH-NORTH_PATCH_END_UNSTAG',   local_sn_patch_e)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'SOUTH-NORTH_PATCH_START_STAG',   local_sn_patch_s_stag)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'SOUTH-NORTH_PATCH_END_STAG',     local_sn_patch_e_stag)
     CALL ext_put_dom_ti_char          ( fhandle, io_form_output, &
                        'GRIDTYPE', geogrid%grid_type)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'DX', geogrid%dx)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'DY', geogrid%dy)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'DYN_OPT', geogrid%dyn_opt)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'CEN_LAT', geogrid%cen_lat)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'CEN_LON', geogrid%cen_lon)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'TRUELAT1', geogrid%truelat1)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'TRUELAT2', geogrid%truelat2)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'MOAD_CEN_LAT', geogrid%moad_cen_lat)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'STAND_LON', geogrid%stand_lon)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'POLE_LAT', geogrid%pole_lat)
     CALL ext_put_dom_ti_real_scalar   ( fhandle, io_form_output, &
                        'POLE_LON', geogrid%pole_lon)
     CALL ext_put_dom_ti_real_vector   ( fhandle, io_form_output, &
                        'corner_lats', local_corner_lats, 16)
     CALL ext_put_dom_ti_real_vector   ( fhandle, io_form_output, &
                        'corner_lons', local_corner_lons, 16)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'MAP_PROJ', geogrid%map_proj)
     CALL ext_put_dom_ti_char          ( fhandle, io_form_output, &
                        'MMINLU', trim(geogrid%mminlu))
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'NUM_LAND_CAT', geogrid%num_land_cat)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'ISWATER', geogrid%iswater)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'ISLAKE', geogrid%islake)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'ISICE', geogrid%isice)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'ISURBAN', geogrid%isurban)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output, &
                        'ISOILWATER', geogrid%isoilwater)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'grid_id',    geogrid%grid_id)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'parent_id',  geogrid%parent_id)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'i_parent_start',    geogrid%i_parent_start)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'j_parent_start',    geogrid%j_parent_start)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'i_parent_end',      geogrid%i_parent_end)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'j_parent_end',      geogrid%j_parent_end)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                        'parent_grid_ratio', geogrid%parent_grid_ratio)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                                         'FLAG_METGRID', 1)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                                         'sr_x', 1)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                                         'sr_y', 1)
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,        &
                                         'NUM_METGRID_SOIL_LEVELS', metgrid%num_soil_levels )

     IF (dbglvl > 4) WRITE(6,'(7x,a)') '--- Writing flags ...'
     DO i=1,nflags
       IF (dbglvl > 5) WRITE(6,'(9x,a,I2,2a)') '&&& Writing #',i,' flag - ',TRIM(flags(i))
       CALL ext_put_dom_ti_integer_scalar( fhandle,io_form_output,      &
                                                  trim(flags(i)), 1)
     END DO

     !
     ! Newly added in WRFV3.3.1, we still did not adapt boundary-only processing
     !
     CALL ext_put_dom_ti_integer_scalar( fhandle, io_form_output,       &
                                         'FLAG_EXCLUDED_MIDDLE', 0 )

  END IF

  RETURN
END SUBROUTINE write_global_attrs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: ext_put_dom_ti_integer
!
! Purpose: Write a domain time-independent INTEGER attribute to output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ext_put_dom_ti_integer_scalar(fhandle, io_form_output, var_name, var_value)

  USE module_wrfgrid_constants
  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)  :: fhandle
  INTEGER, INTENT(IN)  :: io_form_output
  INTEGER, INTENT(IN) :: var_value
  CHARACTER (LEN=*), INTENT(IN) :: var_name

  ! Local variables
  INTEGER :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IF (io_form_output == BINARY) THEN
         CALL ext_int_put_dom_ti_integer(fhandle, trim(var_name), &
                                         var_value, &
                                         1, istatus)
      END IF
      IF (io_form_output == NETCDF) THEN
         CALL ext_ncd_put_dom_ti_integer(fhandle, trim(var_name), &
                                         var_value, &
                                         1, istatus)
      END IF

!      CALL mprintf((istatus /= 0),ERROR,'Error in writing domain time-independent attribute')

   END SUBROUTINE ext_put_dom_ti_integer_scalar


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: ext_put_dom_ti_integer
   !
   ! Purpose: Write a domain time-independent INTEGER attribute to output.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE ext_put_dom_ti_integer_vector(fhandle, io_form_output,var_name, var_value, n)

      USE module_wrfgrid_constants

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: fhandle
      INTEGER, INTENT(IN)  :: io_form_output
      INTEGER, INTENT(IN) :: n
      INTEGER, dimension(n), INTENT(IN) :: var_value
      CHARACTER (LEN=*), INTENT(IN) :: var_name

      ! Local variables
      INTEGER :: istatus

      IF (io_form_output == BINARY) THEN
         CALL ext_int_put_dom_ti_integer(fhandle, trim(var_name), &
                                         var_value, &
                                         n, istatus)
      END IF
      IF (io_form_output == NETCDF) THEN
         CALL ext_ncd_put_dom_ti_integer(fhandle, trim(var_name), &
                                         var_value, &
                                         n, istatus)
      END IF

!      CALL mprintf((istatus /= 0),ERROR,'Error in writing domain time-independent attribute')

END SUBROUTINE ext_put_dom_ti_integer_vector


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: ext_put_dom_ti_real
   !
   ! Purpose: Write a domain time-independent REAL attribute to output.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ext_put_dom_ti_real_scalar(fhandle, io_form_output, var_name, var_value)

  USE module_wrfgrid_constants

      IMPLICIT NONE

      ! Arguments
  INTEGER, INTENT(IN)  :: fhandle
  INTEGER, INTENT(IN)  :: io_form_output
      REAL, INTENT(IN) :: var_value
      CHARACTER (LEN=*), INTENT(IN) :: var_name

      ! Local variables
      INTEGER :: istatus

      IF (io_form_output == BINARY) THEN
         CALL ext_int_put_dom_ti_real(fhandle, trim(var_name), &
                                         var_value, &
                                         1, istatus)
      END IF
      IF (io_form_output == NETCDF) THEN
         CALL ext_ncd_put_dom_ti_real(fhandle, trim(var_name), &
                                         var_value, &
                                         1, istatus)
      END IF

!      CALL mprintf((istatus /= 0),ERROR,'Error in writing domain time-independent attribute')

  RETURN
END SUBROUTINE ext_put_dom_ti_real_scalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: ext_put_dom_ti_real
!
! Purpose: Write a domain time-independent REAL attribute to output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ext_put_dom_ti_real_vector(fhandle, io_form_output, var_name, var_value, n)

  USE module_wrfgrid_constants

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: fhandle
  INTEGER, INTENT(IN)  :: io_form_output

  INTEGER, INTENT(IN) :: n
  REAL, dimension(n), INTENT(IN) :: var_value
  CHARACTER (LEN=*),  INTENT(IN) :: var_name

  ! Local variables
  INTEGER :: istatus

  IF (io_form_output == BINARY) THEN
         CALL ext_int_put_dom_ti_real(fhandle, trim(var_name), &
                                         var_value, &
                                         n, istatus)
  END IF
  IF (io_form_output == NETCDF) THEN
     CALL ext_ncd_put_dom_ti_real(fhandle, trim(var_name), &
                                     var_value, &
                                     n, istatus)
  END IF

!      CALL mprintf((istatus /= 0),ERROR,'Error in writing domain time-independent attribute')

  RETURN
END SUBROUTINE ext_put_dom_ti_real_vector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: ext_put_dom_ti_char
!
! Purpose: Write a domain time-independent CHARACTER attribute to output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ext_put_dom_ti_char(fhandle, io_form_output, var_name, var_value)

  USE module_wrfgrid_constants

  IMPLICIT NONE

  INTEGER,           INTENT(IN) :: fhandle, io_form_output
  CHARACTER (LEN=*), INTENT(IN) :: var_name, var_value

  ! Local variables
  INTEGER :: istatus

  IF (io_form_output == BINARY) THEN
    CALL ext_int_put_dom_ti_char(fhandle, trim(var_name),trim(var_value), istatus)
  END IF

  IF (io_form_output == NETCDF) THEN
    CALL ext_ncd_put_dom_ti_char(fhandle, trim(var_name), trim(var_value), istatus)
  END IF

!      CALL mprintf((istatus /= 0),ERROR,'Error in writing domain time-independent attribute')

  RETURN
END SUBROUTINE ext_put_dom_ti_char


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: output_close
!
! Purpose: Finalizes all output. This may include closing windows, calling I/O
!    API termination routines, or closing files.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE output_close(fhandle, io_form_output, IAMIO_NODE, istatus)

  USE module_wrfgrid_constants

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: fhandle
  INTEGER, INTENT(IN)  :: io_form_output
  LOGICAL, INTENT(IN)  :: IAMIO_NODE
  INTEGER, INTENT(OUT) :: istatus

  ! Local variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF ( IAMIO_NODE ) THEN

    IF (io_form_output == BINARY) THEN
      CALL ext_int_ioclose(fhandle, istatus)
      CALL ext_int_ioexit(istatus)
    ELSE IF (io_form_output == NETCDF) THEN
      CALL ext_ncd_ioclose(fhandle, istatus)
      CALL ext_ncd_ioexit(istatus)
    END IF

!         CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_ioclose')
!         CALL mprintf((istatus /= 0), ERROR, 'Error in ext_pkg_ioexit')
  END IF

  RETURN
END SUBROUTINE output_close
