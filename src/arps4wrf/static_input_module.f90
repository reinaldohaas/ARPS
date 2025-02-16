MODULE static_input_module

  USE module_arpsgrid_constants
  USE module_wrfgrid_constants
  USE queue_module
  USE wrf_parallel_module

  TYPE (queue) :: unit_desc

  INTEGER :: handle
  INTEGER :: num_calls

  CHARACTER(LEN=MAXFILELEN) :: filename
  INTEGER :: io_form
  LOGICAL :: IAMIONODE
  LOGICAL :: tiled_input  ! true for serial run and split input

  CHARACTER (LEN=1) :: internal_gridtype

  CONTAINS

  SUBROUTINE input_init(input_fname,io_form_input,IAMIO_NODE,do_tiled_input, &
                        istatus)
!-----------------------------------------------------------------------
!
! Purpose:
!   Initialize WRF io module (i.e. convert a file name to a file handle)
!
!-----------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: input_fname
    INTEGER, INTENT(IN) :: io_form_input
    LOGICAL, INTENT(IN) :: IAMIO_NODE
    LOGICAL, INTENT(IN) :: do_tiled_input

    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
! Local variables

    INTEGER :: lfinput
    INTEGER :: comm_1, comm_2

    istatus = 0

    IF (IAMIO_NODE) THEN

      IF (io_form_input == BINARY) THEN
        CALL ext_int_ioinit('sysdep info', istatus)
      ELSE IF (io_form_input == NETCDF) THEN
        CALL ext_ncd_ioinit('sysdep info', istatus)
      ELSE
        WRITE(STDOUT,'(1x,a,I3)') 'ERROR: Unknown input format = ',       &
                                  io_form_input
        istatus = -1
        CALL arpsstop('Unknow IO format',1)
      END IF

      filename = input_fname
      IF (do_tiled_input) THEN
        lfinput = LEN_TRIM(filename)
        WRITE(filename(lfinput+1:lfinput+5), '(a1,i4.4)') '_', my_proc_id
      END IF

      istatus = 0
      comm_1  = 1
      comm_2  = 1
      IF (io_form_input == BINARY) THEN
        CALL ext_int_open_for_read(trim(filename), comm_1, comm_2,     &
                                   'sysdep info', handle, istatus)
      ELSE IF (io_form_input == NETCDF) THEN
        CALL ext_ncd_open_for_read(trim(filename), comm_1, comm_2,     &
                                   'sysdep info', handle, istatus)
      END IF

      CALL q_init(unit_desc)

    END IF ! ( IO_NODE )

    IAMIONODE   = IAMIO_NODE              ! initialize module variables
    io_form     = io_form_input
    tiled_input = do_tiled_input
    num_calls   = 0

    RETURN
  END SUBROUTINE input_init

  SUBROUTINE read_next_field(start_patch_i, end_patch_i, &
                             start_patch_j, end_patch_j, &
                             start_patch_k, end_patch_k, &
                             cname, cunits, cdesc, memorder, stagger, &
                             dimnames, real_array, istatus)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: start_patch_i, end_patch_i, &
                            start_patch_j, end_patch_j, &
                            start_patch_k, end_patch_k
    CHARACTER (LEN=*), INTENT(OUT) :: cname, memorder, stagger, cunits, cdesc

    REAL, POINTER, DIMENSION(:,:,:) :: real_array
    CHARACTER (LEN=128), DIMENSION(3), INTENT(INOUT) :: dimnames
    INTEGER, INTENT(INOUT) :: istatus

  !-----------------------------------------------------------------------
  ! Local variables
    INTEGER :: ndim, wrftype
    INTEGER :: sm1, em1, sm2, em2, sm3, em3, sp1, ep1, sp2, ep2, sp3, ep3
    INTEGER, DIMENSION(3) :: domain_start, domain_end
    CHARACTER (LEN=20)    :: datestr

    REAL, POINTER, DIMENSION(:,:,:) :: real_domain

    TYPE(Q_DATA) :: qd

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF (IAMIONODE) THEN

       IF (num_calls == 0 ) then
         IF (io_form == BINARY) THEN
           CALL ext_int_get_next_time(handle, datestr, istatus)
         ELSE IF (io_form == NETCDF) THEN
           CALL ext_ncd_get_next_time(handle, datestr, istatus)
         END IF
       END IF
       num_calls = num_calls + 1

       IF (io_form == BINARY) THEN
         CALL ext_int_get_next_var(handle, cname, istatus)
       ELSE IF (io_form == NETCDF) THEN
         CALL ext_ncd_get_next_var(handle, cname, istatus)
       END IF
    END IF

    IF (ntprocs > 1 .AND. .NOT. TILED_INPUT) CALL parallel_bcast_int(istatus)
    IF (istatus /= 0) RETURN

    IF (IAMIONODE) THEN

       IF (io_form == BINARY) THEN
         CALL ext_int_get_var_info(handle, cname, ndim, memorder,         &
                      stagger, domain_start, domain_end, wrftype, istatus)
       ELSE IF (io_form == NETCDF) THEN
         CALL ext_ncd_get_var_info(handle, cname, ndim, memorder,         &
                      stagger, domain_start, domain_end, wrftype, istatus)
       END IF

       IF (istatus /= 0) CALL arpsstop(                                 &
       'ERROR:In read_next_field(), problems with ext_pkg_get_var_info()',1)

       start_patch_i = domain_start(1)
       start_patch_j = domain_start(2)
       end_patch_i   = domain_end(1)
       end_patch_j   = domain_end(2)
       IF (ndim == 3) THEN
         start_patch_k = domain_start(3)
         end_patch_k   = domain_end(3)
       ELSE
         domain_start(3) = 1
         domain_end(3) = 1
         start_patch_k = 1
         end_patch_k = 1
       END IF

       NULLIFY(real_domain)

       ALLOCATE(real_domain(start_patch_i:end_patch_i,                    &
                            start_patch_j:end_patch_j,                    &
                            start_patch_k:end_patch_k), STAT = istatus)

       IF (io_form == BINARY) THEN
         CALL ext_int_read_field(handle, '0000-00-00_00:00:00', cname,    &
                            real_domain, WRF_REAL,                        &
                            1, 1, 0, memorder, stagger,                   &
                            dimnames, domain_start, domain_end,           &
                            domain_start, domain_end,                     &
                            domain_start, domain_end, istatus)
          qd = q_remove(unit_desc)
          cunits  = qd%units
          cdesc   = qd%description
          stagger = qd%stagger
       ELSE IF (io_form == NETCDF) THEN
         CALL ext_ncd_read_field(handle, '0000-00-00_00:00:00', cname,    &
                            real_domain, WRF_REAL,                        &
                            1, 1, 0, memorder, stagger,                   &
                            dimnames, domain_start, domain_end,           &
                            domain_start, domain_end,                     &
                            domain_start, domain_end, istatus)
         cunits = ' '
         cdesc = ' '
         stagger = ' '
         CALL ext_ncd_get_var_ti_char(handle, 'units', cname, cunits, istatus)
         CALL ext_ncd_get_var_ti_char(handle, 'description', cname, cdesc, istatus)
         CALL ext_ncd_get_var_ti_char(handle, 'stagger', cname, stagger, istatus)
       END IF

     END IF ! (IO_NODE )

     IF (ntprocs > 1 .AND. .NOT. TILED_INPUT) THEN
       call parallel_bcast_char(cname, len(cname))
       call parallel_bcast_char(cunits, len(cunits))
       call parallel_bcast_char(cdesc, len(cdesc))
       call parallel_bcast_char(memorder, len(memorder))
       call parallel_bcast_char(stagger, len(stagger))
       call parallel_bcast_char(dimnames(1), 128)
       call parallel_bcast_char(dimnames(2), 128)
       call parallel_bcast_char(dimnames(3), 128)
       call parallel_bcast_int(domain_start(3))
       call parallel_bcast_int(domain_end(3))

       sp1 = my_minx
       ep1 = my_maxx - 1
       sp2 = my_miny
       ep2 = my_maxy - 1
       sp3 = domain_start(3)
       ep3 = domain_end(3)

       IF (internal_gridtype == 'C') then
          IF (my_x /= nxprocs - 1 .OR. stagger == 'U') THEN
              ep1 = ep1 + 1
          END IF
          IF (my_y /= nyprocs - 1 .OR. stagger == 'V') THEN
              ep2 = ep2 + 1
          END IF
       ELSE IF (internal_gridtype == 'E') THEN
          ep1 = ep1 + 1
          ep2 = ep2 + 1
       END IF

       sm1 = sp1
       em1 = ep1
       sm2 = sp2
       em2 = ep2
       sm3 = sp3
       em3 = ep3

       start_patch_i = sp1
       end_patch_i   = ep1
       start_patch_j = sp2
       end_patch_j   = ep2
       start_patch_k = sp3
       end_patch_k   = ep3

       ALLOCATE(real_array(sm1:em1,sm2:em2,sm3:em3))
       IF (.NOT. IAMIONODE) THEN
          ALLOCATE(real_domain(1,1,1))
          domain_start(1) = 1
          domain_start(2) = 1
          domain_start(3) = 1
          domain_end(1) = 1
          domain_end(2) = 1
          domain_end(3) = 1
       END IF
       CALL scatter_whole_field_r(real_array,                             &
                                  sm1, em1, sm2, em2, sm3, em3,           &
                                  sp1, ep1, sp2, ep2, sp3, ep3,           &
                                  real_domain,                            &
                                  domain_start(1), domain_end(1),         &
                                  domain_start(2), domain_end(2),         &
                                  domain_start(3), domain_end(3))
       DEALLOCATE(real_domain)

    ELSE ! serial or MPI tiled input

       real_array => real_domain

    END IF

    RETURN
  END SUBROUTINE read_next_field

  SUBROUTINE read_global_attrs(title, start_date, grid_type, dyn_opt,      &
                 west_east_dim, south_north_dim, bottom_top_dim,           &
                 we_patch_s, we_patch_e, we_patch_s_stag, we_patch_e_stag, &
                 sn_patch_s, sn_patch_e, sn_patch_s_stag, sn_patch_e_stag, &
                 map_proj, mminlu, num_land_cat,                        &
                 is_water, is_ice, is_urban, isoilwater,                &
                 grid_id, parent_id, i_parent_start, j_parent_start,                   &
                 i_parent_end, j_parent_end, dx, dy, cen_lat, moad_cen_lat, cen_lon,   &
                 stand_lon, truelat1, truelat2, pole_lat, pole_lon, parent_grid_ratio, &
                 corner_lats, corner_lons, istatus)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: dyn_opt, west_east_dim, south_north_dim, bottom_top_dim
    INTEGER, INTENT(OUT) :: we_patch_s, we_patch_e, we_patch_s_stag, we_patch_e_stag
    INTEGER, INTENT(OUT) :: sn_patch_s, sn_patch_e, sn_patch_s_stag, sn_patch_e_stag
    INTEGER, INTENT(OUT) :: map_proj
    INTEGER, INTENT(OUT) :: num_land_cat,is_water, is_ice, is_urban, isoilwater
    INTEGER, INTENT(OUT) :: grid_id, parent_id, i_parent_start, j_parent_start, &
                            i_parent_end, j_parent_end, parent_grid_ratio
    REAL,    INTENT(OUT) :: dx, dy
    REAL,    INTENT(OUT) :: cen_lat, moad_cen_lat, cen_lon, stand_lon,    &
                            truelat1, truelat2, pole_lat, pole_lon
    REAL, DIMENSION(16), INTENT(OUT) :: corner_lats, corner_lons
    CHARACTER (LEN=128), INTENT(OUT) :: title, start_date, grid_type, mminlu

    INTEGER, INTENT(OUT) :: istatus

  !-----------------------------------------------------------------------
  ! Local variables
    integer ::  i
    character (len=128) :: cunits, cdesc, cstagger
    TYPE (Q_DATA) :: qd

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    if (IAMIONODE) then

      if (io_form == BINARY) then
        istatus = 0
        do while (istatus == 0)

           cunits = ' '
           cdesc = ' '
           cstagger = ' '
           CALL ext_int_get_var_ti_char(handle, 'units', 'VAR', cunits, istatus)

           if (istatus == 0) then
             CALL ext_int_get_var_ti_char(handle, 'description', 'VAR', cdesc, istatus)

             if (istatus == 0) then
                CALL ext_int_get_var_ti_char(handle, 'stagger', 'VAR', cstagger, istatus)

                qd%units       = cunits
                qd%description = cdesc
                qd%stagger     = cstagger
                call q_insert(unit_desc, qd)

             end if
           end if

        end do
      end if

      CALL ext_get_dom_ti_char          ('TITLE', title)
      CALL ext_get_dom_ti_char          ('SIMULATION_START_DATE', start_date)
      CALL ext_get_dom_ti_integer_scalar('WEST-EAST_GRID_DIMENSION', west_east_dim)
      CALL ext_get_dom_ti_integer_scalar('SOUTH-NORTH_GRID_DIMENSION', south_north_dim)
      CALL ext_get_dom_ti_integer_scalar('BOTTOM-TOP_GRID_DIMENSION', bottom_top_dim)
      CALL ext_get_dom_ti_integer_scalar('WEST-EAST_PATCH_START_UNSTAG', we_patch_s)
      CALL ext_get_dom_ti_integer_scalar('WEST-EAST_PATCH_END_UNSTAG', we_patch_e)
      CALL ext_get_dom_ti_integer_scalar('WEST-EAST_PATCH_START_STAG', we_patch_s_stag)
      CALL ext_get_dom_ti_integer_scalar('WEST-EAST_PATCH_END_STAG', we_patch_e_stag)
      CALL ext_get_dom_ti_integer_scalar('SOUTH-NORTH_PATCH_START_UNSTAG', sn_patch_s)
      CALL ext_get_dom_ti_integer_scalar('SOUTH-NORTH_PATCH_END_UNSTAG', sn_patch_e)
      CALL ext_get_dom_ti_integer_scalar('SOUTH-NORTH_PATCH_START_STAG', sn_patch_s_stag)
      CALL ext_get_dom_ti_integer_scalar('SOUTH-NORTH_PATCH_END_STAG', sn_patch_e_stag)
      CALL ext_get_dom_ti_char          ('GRIDTYPE', grid_type)
      CALL ext_get_dom_ti_real_scalar   ('DX', dx)
      CALL ext_get_dom_ti_real_scalar   ('DY', dy)
      CALL ext_get_dom_ti_integer_scalar('DYN_OPT', dyn_opt)
      CALL ext_get_dom_ti_real_scalar   ('CEN_LAT', cen_lat)
      CALL ext_get_dom_ti_real_scalar   ('CEN_LON', cen_lon)
      CALL ext_get_dom_ti_real_scalar   ('TRUELAT1', truelat1)
      CALL ext_get_dom_ti_real_scalar   ('TRUELAT2', truelat2)
      CALL ext_get_dom_ti_real_scalar   ('MOAD_CEN_LAT', moad_cen_lat)
      CALL ext_get_dom_ti_real_scalar   ('STAND_LON', stand_lon)
      CALL ext_get_dom_ti_real_scalar   ('POLE_LAT', pole_lat)
      CALL ext_get_dom_ti_real_scalar   ('POLE_LON', pole_lon)
      CALL ext_get_dom_ti_real_vector   ('corner_lats', corner_lats, 16)
      CALL ext_get_dom_ti_real_vector   ('corner_lons', corner_lons, 16)
      CALL ext_get_dom_ti_integer_scalar('MAP_PROJ', map_proj)
      CALL ext_get_dom_ti_char          ('MMINLU', mminlu)
      CALL ext_get_dom_ti_integer_scalar('NUM_LAND_CAT', num_land_cat)
      CALL ext_get_dom_ti_integer_scalar('ISWATER', is_water)
      CALL ext_get_dom_ti_integer_scalar('ISICE', is_ice)
      CALL ext_get_dom_ti_integer_scalar('ISURBAN', is_urban)
      CALL ext_get_dom_ti_integer_scalar('ISOILWATER', isoilwater)
      CALL ext_get_dom_ti_integer_scalar('grid_id', grid_id)
      CALL ext_get_dom_ti_integer_scalar('parent_id', parent_id)
      CALL ext_get_dom_ti_integer_scalar('i_parent_start', i_parent_start)
      CALL ext_get_dom_ti_integer_scalar('j_parent_start', j_parent_start)
      CALL ext_get_dom_ti_integer_scalar('i_parent_end', i_parent_end)
      CALL ext_get_dom_ti_integer_scalar('j_parent_end', j_parent_end)
      CALL ext_get_dom_ti_integer_scalar('parent_grid_ratio', parent_grid_ratio)

    END IF

    IF (ntprocs > 1 .AND. .NOT. TILED_INPUT) THEN

       call parallel_bcast_char(title, len(title))
       call parallel_bcast_char(start_date, len(start_date))
       call parallel_bcast_char(grid_type, len(grid_type))
       call parallel_bcast_int(west_east_dim)
       call parallel_bcast_int(south_north_dim)
       call parallel_bcast_int(bottom_top_dim)
       call parallel_bcast_int(we_patch_s)
       call parallel_bcast_int(we_patch_e)
       call parallel_bcast_int(we_patch_s_stag)
       call parallel_bcast_int(we_patch_e_stag)
       call parallel_bcast_int(sn_patch_s)
       call parallel_bcast_int(sn_patch_e)
       call parallel_bcast_int(sn_patch_s_stag)
       call parallel_bcast_int(sn_patch_e_stag)

       call parallel_bcast_real(dx)
       call parallel_bcast_real(dy)
       call parallel_bcast_int(dyn_opt)
       call parallel_bcast_real(cen_lat)
       call parallel_bcast_real(cen_lon)
       call parallel_bcast_real(truelat1)
       call parallel_bcast_real(truelat2)
       call parallel_bcast_real(pole_lat)
       call parallel_bcast_real(pole_lon)
       call parallel_bcast_real(moad_cen_lat)
       call parallel_bcast_real(stand_lon)
       do i=1,16
          call parallel_bcast_real(corner_lats(i))
          call parallel_bcast_real(corner_lons(i))
       end do
       call parallel_bcast_int(map_proj)
       call parallel_bcast_char(mminlu, len(mminlu))
       call parallel_bcast_int(num_land_cat)
       call parallel_bcast_int(is_water)
       call parallel_bcast_int(is_ice)
       call parallel_bcast_int(is_urban)
       call parallel_bcast_int(isoilwater)
       call parallel_bcast_int(grid_id)
       call parallel_bcast_int(parent_id)
       call parallel_bcast_int(i_parent_start)
       call parallel_bcast_int(i_parent_end)
       call parallel_bcast_int(j_parent_start)
       call parallel_bcast_int(j_parent_end)
       call parallel_bcast_int(parent_grid_ratio)

    END IF

    internal_gridtype = grid_type(1:1)

    RETURN
  END SUBROUTINE read_global_attrs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Name: ext_get_dom_ti_integer
  !
  ! Purpose: Read a domain time-independent integer attribute from input.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ext_get_dom_ti_integer_scalar(var_name, var_value)

    IMPLICIT NONE

    ! Arguments
    integer,          intent(out) :: var_value
    character (len=*), intent(in) :: var_name

    ! Local variables
    integer :: istatus, outcount

    if (io_form == BINARY) then
      CALL ext_int_get_dom_ti_integer(handle, trim(var_name), &
                             var_value, 1, outcount, istatus)
    ELSE if (io_form == NETCDF) then
      CALL ext_ncd_get_dom_ti_integer(handle, trim(var_name), &
                             var_value, 1, outcount, istatus)
    END IF

    RETURN
  END SUBROUTINE ext_get_dom_ti_integer_scalar


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Name: ext_get_dom_ti_integer
  !
  ! Purpose: Read a domain time-independent integer attribute from input.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ext_get_dom_ti_integer_vector(var_name, var_value, n)

    IMPLICIT NONE

    integer,               intent(in)  :: n
    integer, dimension(n), intent(out) :: var_value
    character (len=*),     intent(in)  :: var_name

    ! Local variables
    integer :: istatus, outcount

    if (io_form == BINARY) then
      CALL ext_int_get_dom_ti_integer(handle, trim(var_name), &
                                      var_value, n, outcount, istatus)
    ELSE if (io_form == NETCDF) then
      CALL ext_ncd_get_dom_ti_integer(handle, trim(var_name), &
                                      var_value, n, outcount, istatus)
    end if

    RETURN
  END SUBROUTINE ext_get_dom_ti_integer_vector


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Name: ext_get_dom_ti_real
  !
  ! Purpose: Read a domain time-independent real attribute from input.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ext_get_dom_ti_real_scalar(var_name, var_value)

      implicit none

      ! Arguments
      real, intent(out) :: var_value
      character (len=*), intent(in) :: var_name

      ! Local variables
      integer :: istatus, outcount

      if (io_form == BINARY) then
         CALL ext_int_get_dom_ti_real(handle, trim(var_name), &
                                         var_value, &
                                         1, outcount, istatus)
      else if (io_form == NETCDF) then
         CALL ext_ncd_get_dom_ti_real(handle, trim(var_name), &
                                         var_value, &
                                         1, outcount, istatus)
      end if

      if (istatus /= 0) CALL arpsstop(  &
      'ERROR:Error while reading domain time-independent attribute.',1)

   END SUBROUTINE ext_get_dom_ti_real_scalar


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Name: ext_get_dom_ti_real
  !
  ! Purpose: Read a domain time-independent real attribute from input.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ext_get_dom_ti_real_vector(var_name, var_value, n)

     implicit none

     ! Arguments
     integer, intent(in) :: n
     real, dimension(n), intent(out) :: var_value
     character (len=*), intent(in) :: var_name

     ! Local variables
     integer :: istatus, outcount

     if (io_form == BINARY) then
        CALL ext_int_get_dom_ti_real(handle, trim(var_name), &
                                        var_value, &
                                        n, outcount, istatus)
     else if (io_form == NETCDF) then
        CALL ext_ncd_get_dom_ti_real(handle, trim(var_name), &
                                        var_value, &
                                        n, outcount, istatus)
     end if

     if (istatus /= 0) CALL arpsstop(  &
     'ERROR:Error while reading domain time-independent attribute.',1)

  END SUBROUTINE ext_get_dom_ti_real_vector


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Name: ext_get_dom_ti_char
  !
  ! Purpose: Read a domain time-independent character attribute from input.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ext_get_dom_ti_char(var_name, var_value)

     implicit none

     ! Arguments
     character (len=*), intent(in) :: var_name
     character (len=128), intent(out) :: var_value

     ! Local variables
     integer :: istatus

     if (io_form == BINARY) then
        CALL ext_int_get_dom_ti_char(handle, trim(var_name), &
                                        var_value, &
                                        istatus)
     else if (io_form == NETCDF) then
        CALL ext_ncd_get_dom_ti_char(handle, trim(var_name), &
                                        var_value, &
                                        istatus)
     end if

     if (istatus /= 0) CALL arpsstop(  &
     'ERROR:Error while reading domain time-independent attribute.',1)

  END SUBROUTINE ext_get_dom_ti_char


  SUBROUTINE input_close(istatus)

    IMPLICIT NONE

    integer, INTENT(OUT) :: istatus

    istatus = 0
    IF (IAMIONODE) THEN
       IF (io_form == BINARY) THEN
          CALL ext_int_ioclose(handle, istatus)
          CALL ext_int_ioexit(istatus)
       ELSE IF (io_form == NETCDF) THEN
          CALL ext_ncd_ioclose(handle, istatus)
          CALL ext_ncd_ioexit(istatus)
       END IF
    END IF

    CALL q_destroy(unit_desc)

    RETURN
  END SUBROUTINE input_close

END MODULE static_input_module
