MODULE module_metgrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains data for destination grid (WPS metgrid) for
! program "real.exe" and it also handles WRF(WPS) related IO
! (writes to metgrid files).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  USE module_constants
  USE module_geogrid

  TYPE type_metgrid
    INTEGER :: nfields

    TYPE (type_field), POINTER :: fields(:)

    INTEGER  :: pres_no
    INTEGER  :: tt_no, qv_no
    INTEGER  :: landsea_no
    INTEGER  :: seaice_no
    INTEGER  :: num_soil_levels

    INTEGER            :: nflags
    CHARACTER(LEN=128) :: flags(MAXFLAGS)

  END TYPE type_metgrid

  TYPE (type_metgrid) :: metgrid

  LOGICAL :: met_allocated = .FALSE.

  CONTAINS

  SUBROUTINE alloc_metgrid(istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Allocate arrays in "metgrid" and initialize them to some values
!    (missing).
!
!-----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ALLOCATE(metgrid%fields(MAXFIELDS), STAT = istatus)

    metgrid%nfields = 0
    metgrid%nflags  = 0
    metgrid%pres_no    = 0
    metgrid%tt_no      = 0
    metgrid%qv_no      = 0
    metgrid%landsea_no = 0
    metgrid%seaice_no  = 0

    met_allocated = .TRUE.

    RETURN
  END SUBROUTINE alloc_metgrid


  SUBROUTINE metgrid_add_3d_field(dta_arr,fieldname_in,units_in,descr_in,stagger, &
                                  sp1,ep1,sp2,ep2,ps3,istatus)

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: sp1,ep1, sp2, ep2, ps3
    INTEGER, INTENT(OUT) :: istatus
    REAL,     TARGET,   INTENT(IN) :: dta_arr(sp1:ep1,sp2:ep2,1:ps3)
    CHARACTER(LEN=128), INTENT(IN) :: fieldname_in,units_in,descr_in
    CHARACTER(LEN=1),   INTENT(IN) :: stagger

!-----------------------------------------------------------------------

    INTEGER :: i, j, k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
    metgrid%nfields = metgrid%nfields + 1

    metgrid%fields(metgrid%nfields)%fieldname = fieldname_in
    metgrid%fields(metgrid%nfields)%units     = units_in
    metgrid%fields(metgrid%nfields)%descr     = descr_in

    metgrid%fields(metgrid%nfields)%patch_start(1) = sp1;
    metgrid%fields(metgrid%nfields)%patch_end(1)   = ep1
    metgrid%fields(metgrid%nfields)%patch_start(2) = sp2;
    metgrid%fields(metgrid%nfields)%patch_end(2)   = ep2
    metgrid%fields(metgrid%nfields)%mem_start(1)   = sp1;
    metgrid%fields(metgrid%nfields)%mem_end(1)     = ep1
    metgrid%fields(metgrid%nfields)%mem_start(2)   = sp2;
    metgrid%fields(metgrid%nfields)%mem_end(2)     = ep2
    metgrid%fields(metgrid%nfields)%dom_start(3)   = 1;
    metgrid%fields(metgrid%nfields)%dom_end(3)     = ps3
    metgrid%fields(metgrid%nfields)%patch_start(3) = 1;
    metgrid%fields(metgrid%nfields)%patch_end(3)   = ps3
    metgrid%fields(metgrid%nfields)%mem_start(3)   = 1;
    metgrid%fields(metgrid%nfields)%mem_end(3)     = ps3
    IF (ps3 > 1 ) THEN
      metgrid%fields(metgrid%nfields)%mem_order   = 'XYZ'
      metgrid%fields(metgrid%nfields)%dimnames(3) = 'num_metgrid_levels'
    ELSE
      metgrid%fields(metgrid%nfields)%mem_order = 'XY'
      metgrid%fields(metgrid%nfields)%dimnames(3)    = ' '
    END IF

    IF (geogrid%grid_type == 'C') THEN
      IF (stagger == 'M') THEN
        metgrid%fields(metgrid%nfields)%stagger = M
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      ELSE IF (stagger == 'U') THEN
        metgrid%fields(metgrid%nfields)%stagger = U
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e + 1
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east_stag'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      ELSE IF (stagger == 'V') THEN
        metgrid%fields(metgrid%nfields)%stagger = V
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e + 1
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north_stag'
      ELSE IF (stagger == 'W') THEN  ! Actually all other fields have an extra layers (nz-1) from ARPS
        ! Since usually level (nz-1) is assign meaningful values in the ARPS system, for example ext2arps etc.
        metgrid%fields(metgrid%nfields)%stagger = M
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      ELSE
        istatus = -1
      END IF

    ELSE IF (geogrid%grid_type == 'E') THEN
      metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
      metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
      metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
      metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
      metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
      IF (stagger == 'M') THEN
        metgrid%fields(metgrid%nfields)%stagger = HH
      ELSE IF (stagger == 'U') THEN
        metgrid%fields(metgrid%nfields)%stagger = VV
      ELSE IF (stagger == 'V') THEN
        metgrid%fields(metgrid%nfields)%stagger = VV
      ELSE
        istatus = -1
      END IF
    ELSE
      istatus = -2
    END IF

    IF (istatus == 0) THEN
      ALLOCATE(metgrid%fields(metgrid%nfields)%rdata_arr(sp1:ep1,sp2:ep2,1:ps3), &
                                                        STAT = istatus)
      DO k = 1, ps3
        DO j = sp2, ep2
          DO i = sp1, ep1
            metgrid%fields(metgrid%nfields)%rdata_arr(i,j,k) = dta_arr(i,j,k)
          END DO
        END DO
      END DO
!      metgrid%fields(metgrid%nfields)%rdata_arr => dta_arr

      IF (TRIM(fieldname_in) == 'PRES') THEN
         metgrid%pres_no = metgrid%nfields
      ELSE IF (TRIM(fieldname_in) == 'TT') THEN
         metgrid%tt_no = metgrid%nfields
      ELSE IF (TRIM(fieldname_in) == 'QV') THEN
         metgrid%qv_no = metgrid%nfields
      ELSE IF (TRIM(fieldname_in) == 'LANDSEA') THEN
         metgrid%landsea_no = metgrid%nfields
      ELSE IF (TRIM(fieldname_in) == 'SEAICE') THEN
         metgrid%seaice_no = metgrid%nfields
      END IF

    ELSE
      metgrid%nfields = metgrid%nfields - 1
    END IF

    RETURN
  END SUBROUTINE metgrid_add_3d_field

  SUBROUTINE metgrid_add_3dsoil_field(dta_arr,fieldname_in,units_in,    &
                descr_in,stagger, sp1,ep1,sp2,ep2,ps3,istatus)

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: sp1,ep1, sp2, ep2, ps3
    INTEGER, INTENT(OUT) :: istatus
    REAL,     TARGET,   INTENT(IN) :: dta_arr(sp1:ep1,sp2:ep2,1:ps3)
    CHARACTER(LEN=128), INTENT(IN) :: fieldname_in,units_in,descr_in
    CHARACTER(LEN=1),   INTENT(IN) :: stagger

!-----------------------------------------------------------------------

    INTEGER :: i, j, k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
    metgrid%nfields = metgrid%nfields + 1

    metgrid%fields(metgrid%nfields)%units     = units_in
    metgrid%fields(metgrid%nfields)%descr     = descr_in

    metgrid%fields(metgrid%nfields)%patch_start(1) = sp1;
    metgrid%fields(metgrid%nfields)%patch_end(1)   = ep1
    metgrid%fields(metgrid%nfields)%patch_start(2) = sp2;
    metgrid%fields(metgrid%nfields)%patch_end(2)   = ep2
    metgrid%fields(metgrid%nfields)%mem_start(1)   = sp1;
    metgrid%fields(metgrid%nfields)%mem_end(1)     = ep1
    metgrid%fields(metgrid%nfields)%mem_start(2)   = sp2;
    metgrid%fields(metgrid%nfields)%mem_end(2)     = ep2
    metgrid%fields(metgrid%nfields)%dom_start(3)   = 1;
    metgrid%fields(metgrid%nfields)%dom_end(3)     = ps3
    metgrid%fields(metgrid%nfields)%patch_start(3) = 1;
    metgrid%fields(metgrid%nfields)%patch_end(3)   = ps3
    metgrid%fields(metgrid%nfields)%mem_start(3)   = 1;
    metgrid%fields(metgrid%nfields)%mem_end(3)     = ps3

    metgrid%fields(metgrid%nfields)%fieldname = fieldname_in
    IF (ps3 > 1 ) THEN
      metgrid%fields(metgrid%nfields)%mem_order = 'XYZ'
      IF (geogrid%grid_type == 'E') THEN
        metgrid%fields(metgrid%nfields)%fieldname = TRIM(fieldname_in)//'C_WPS'
      END IF
      IF (TRIM(fieldname_in) == 'ST' .OR. TRIM(fieldname_in) == 'SOIL_LAYERS') THEN
        IF (geogrid%grid_type == 'E') THEN
          metgrid%fields(metgrid%nfields)%dimnames(3) = 'num_st_levels'  ! NMM used the old names
        ELSE
          metgrid%fields(metgrid%nfields)%dimnames(3) = 'num_st_layers'  ! ARW changed to this name since WRFV3.2.1
        END IF
      ELSE
        IF (geogrid%grid_type == 'E') THEN
          metgrid%fields(metgrid%nfields)%dimnames(3) = 'num_sm_levels'
        ELSE
          metgrid%fields(metgrid%nfields)%dimnames(3) = 'num_sm_layers'
        END IF
      END IF
    ELSE
      metgrid%fields(metgrid%nfields)%mem_order = 'XY'
      metgrid%fields(metgrid%nfields)%dimnames(3)    = ' '
    END IF

    IF (geogrid%grid_type == 'C') THEN
      IF (stagger == 'M') THEN
        metgrid%fields(metgrid%nfields)%stagger = M
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      ELSE IF (stagger == 'U') THEN
        metgrid%fields(metgrid%nfields)%stagger = U
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e + 1
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east_stag'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      ELSE IF (stagger == 'V') THEN
        metgrid%fields(metgrid%nfields)%stagger = V
        metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
        metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
        metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
        metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e + 1
        metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north_stag'
      ELSE
        istatus = -1
      END IF

    ELSE IF (geogrid%grid_type == 'E') THEN
      metgrid%fields(metgrid%nfields)%dimnames(1)  = 'west_east'
      metgrid%fields(metgrid%nfields)%dom_start(1) = geogrid%we_dom_s;
      metgrid%fields(metgrid%nfields)%dom_end(1)   = geogrid%we_dom_e
      metgrid%fields(metgrid%nfields)%dimnames(2)  = 'south_north'
      metgrid%fields(metgrid%nfields)%dom_start(2) = geogrid%sn_dom_s;
      metgrid%fields(metgrid%nfields)%dom_end(2)   = geogrid%sn_dom_e
      IF (stagger == 'M') THEN
        metgrid%fields(metgrid%nfields)%stagger = HH
      ELSE IF (stagger == 'U') THEN
        metgrid%fields(metgrid%nfields)%stagger = VV
      ELSE IF (stagger == 'V') THEN
        metgrid%fields(metgrid%nfields)%stagger = VV
      ELSE
        istatus = -1
      END IF
    ELSE
      istatus = -2
    END IF

    IF (istatus == 0) THEN
      ALLOCATE(metgrid%fields(metgrid%nfields)%rdata_arr(sp1:ep1,sp2:ep2,1:ps3), &
                                                        STAT = istatus)
      DO k = 1, ps3
        DO j = sp2, ep2
          DO i = sp1, ep1
            metgrid%fields(metgrid%nfields)%rdata_arr(i,j,k) = dta_arr(i,j,k)
          END DO
        END DO
      END DO
!      metgrid%fields(metgrid%nfields)%rdata_arr => dta_arr

    ELSE
      metgrid%nfields = metgrid%nfields - 1
    END IF

    RETURN
  END SUBROUTINE metgrid_add_3dsoil_field

  SUBROUTINE metgrid_replace_terrain_field(dta_arr,fieldname_in,units_in,descr_in,stagger, &
                                  sp1,ep1,sp2,ep2,istatus)

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: sp1,ep1, sp2, ep2
    INTEGER, INTENT(OUT) :: istatus
    REAL,     TARGET,   INTENT(IN) :: dta_arr(sp1:ep1,sp2:ep2)
    CHARACTER(LEN=128), INTENT(IN) :: fieldname_in,units_in,descr_in
    CHARACTER(LEN=1),   INTENT(IN) :: stagger

!-----------------------------------------------------------------------

    INTEGER :: i, j
    INTEGER :: field_no

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    SELECT CASE (stagger)
    CASE ('M')
      field_no = geogrid%hgt_m_no
    CASE ('U')
      field_no = geogrid%hgt_u_no
    CASE ('V')
      field_no = geogrid%hgt_v_no
    CASE DEFAULT
      WRITE(6,'(1x,3a)') 'ERROR: Unknow stagger = ',stagger,' in metgrid_replace_terrain_field.'
      istatus = -1
      RETURN
    END SELECT

    IF (TRIM(fieldname_in) == TRIM(geogrid%fields(field_no)%fieldname).OR. &
        sp1 == geogrid%fields(field_no)%patch_start(1) .OR.             &
        ep1 == geogrid%fields(field_no)%patch_end(1)   .OR.             &
        sp2 == geogrid%fields(field_no)%patch_start(2) .OR.             &
        ep2 == geogrid%fields(field_no)%patch_end(2) ) THEN
      DO j = sp2,ep2
        DO i = sp1, ep1
          geogrid%fields(field_no)%rdata_arr(i,j,1) = dta_arr(i,j)
        END DO
      END DO
    ELSE
      WRITE(6,'(1x,3a,/,2(7x,a,4(I4,a),/))') 'ERROR: The variable ',TRIM(fieldname_in), &
      ' did not exist or it has wrong dimensions in metgrid_replace_terrain_field.', &
      'Expected (',sp1,'-',ep1,', ',sp2,'-',ep2,'),','but found (',     &
      geogrid%fields(field_no)%patch_start(1),'-',                      &
      geogrid%fields(field_no)%patch_end(1),',',                        &
      geogrid%fields(field_no)%patch_start(2),'-',                      &
      geogrid%fields(field_no)%patch_end(2),').'
      istatus = -2
    END IF

    RETURN
  END SUBROUTINE metgrid_replace_terrain_field


  SUBROUTINE metgrid_add_flag(flagname,always,outgrid,istatus)

!-----------------------------------------------------------------------

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: flagname
    LOGICAL,          INTENT(IN)  :: always
    CHARACTER(LEN=1), INTENT(IN)  :: outgrid
    INTEGER,          INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (always .OR. geogrid%grid_type == outgrid) THEN
      metgrid%nflags = metgrid%nflags + 1
      IF (metgrid%nflags > MAXFLAGS) THEN
        WRITE(6,'(1x,a)') 'ERROR: Too many flags for metgrid - while adding '//flagname
        CALL arpsstop('ERROR: too many flags.',1)
      END IF
      metgrid%flags(metgrid%nflags) = flagname
    END IF

    RETURN
  END SUBROUTINE metgrid_add_flag

  SUBROUTINE dealloc_metgrid(istatus)
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Allocate arrays in "metgrid" and initialize them to some values
!    (missing).
!
!-----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

    INTEGER :: n
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    DO n = 1, metgrid%nfields
      DEALLOCATE(metgrid%fields(n)%rdata_arr, STAT = istatus)
    END DO
    DEALLOCATE(metgrid%fields, STAT = istatus)

    metgrid%nfields = 0
    metgrid%nflags  = 0
    metgrid%pres_no    = 0
    metgrid%tt_no      = 0
    metgrid%qv_no      = 0
    metgrid%landsea_no = 0
    metgrid%seaice_no  = 0

    met_allocated = .FALSE.

    RETURN
  END SUBROUTINE dealloc_metgrid

  SUBROUTINE write_metgrid(IAMROOT,out_path,io_form,ndomain,            &
                           title,datestr,dbglvl,numdigits,out_fname,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write "metgrid" to WPS intermediate format for "real.exe".
!    (missing).
!
!-----------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL,                   INTENT(IN)  :: IAMROOT
    INTEGER,                   INTENT(IN)  :: io_form
    INTEGER,                   INTENT(IN)  :: ndomain
    CHARACTER(LEN=MAXFILELEN), INTENT(IN)  :: out_path
    CHARACTER(LEN=19),         INTENT(IN)  :: datestr
    CHARACTER(LEN=*),          INTENT(IN)  :: title
    INTEGER,                   INTENT(IN)  :: dbglvl
    INTEGER,                   INTENT(IN)  :: numdigits
    CHARACTER(LEN=MAXFILELEN), INTENT(OUT) :: out_fname

    INTEGER, INTENT(OUT) :: istatus
!-----------------------------------------------------------------------

    LOGICAL :: IAMIO_NODE
    LOGICAL :: do_tiled_output
    INTEGER :: io_form_local

    INTEGER :: i, ntotal
    INTEGER :: fhandle

    TYPE(type_field), POINTER :: field

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
    IAMIO_NODE      = .FALSE.
    do_tiled_output = .FALSE.

    IF (io_form > 100) THEN
      io_form_local   = MOD(io_form,100)
      do_tiled_output = .TRUE.
      IAMIO_NODE      = .TRUE.
    ELSE
      IF (IAMROOT) IAMIO_NODE = .TRUE.
      io_form_local = io_form
    END IF

    CALL metgrid_add_flag('FLAG_MF_XY',.FALSE.,'C',istatus)

!-----------------------------------------------------------------------
!
!  Initialize output
!
!-----------------------------------------------------------------------

    IF (dbglvl > 2) WRITE(6,'(3x,a)') '*** Initialize output unit ...'

    CALL output_init(   datestr, ndomain,dbglvl,                        &
                        IAMIO_NODE, do_tiled_output,io_form_local,      &
                        out_path, out_fname, numdigits, fhandle, istatus)

    IF (dbglvl > 2) WRITE(6,'(3x,a)') '*** Writing global attributes ...'

    CALL write_global_attrs(fhandle,io_form_local,                      &
                            IAMIO_NODE, do_tiled_output,dbglvl,         &
                            title, datestr, metgrid%nflags, metgrid%flags, &
                            istatus)

!-----------------------------------------------------------------------
!
! Write fields
!
!-----------------------------------------------------------------------

    DO i=1,metgrid%nfields
      field => metgrid%fields(i)

      IF (dbglvl > 2) WRITE(6,'(3x,a,I2.2,2a)')     &
                      '*** Writing field (',i,') ',TRIM(field%fieldname)

      CALL write_field( fhandle, io_form_local, datestr, IAMIO_NODE,    &
                        do_tiled_output, dbglvl, field, .FALSE., istatus)
    END do

    ntotal = geogrid%nfields+metgrid%nfields
    DO i=geogrid%nfields,1,-1
      field => geogrid%fields(i)

      IF (dbglvl > 2) WRITE(6,'(3x,a,I2.2,2a)')     &
                      '*** Writing field (',ntotal-i+1,') ',TRIM(field%fieldname)

      CALL write_field( fhandle, io_form_local, datestr, IAMIO_NODE,    &
                        do_tiled_output, dbglvl, field, .FALSE., istatus)
    END do

!-----------------------------------------------------------------------
!
! Close output
!
!-----------------------------------------------------------------------
    CALL output_close(fhandle,io_form_local,IAMIO_NODE,istatus)

    RETURN
  END SUBROUTINE write_metgrid

  !#####################################################################

  SUBROUTINE metgrid_qv_to_rh( istatus)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------

    REAL, ALLOCATABLE :: qvs(:,:,:)

    INTEGER :: is,ie, js, je, ks, ke
    INTEGER :: nx, ny, nz
    INTEGER :: i,j,k

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    is = metgrid%fields(metgrid%qv_no)%mem_start(1)
    ie = metgrid%fields(metgrid%qv_no)%mem_end(1)

    js = metgrid%fields(metgrid%qv_no)%mem_start(2)
    je = metgrid%fields(metgrid%qv_no)%mem_end(2)

    ks = metgrid%fields(metgrid%qv_no)%mem_start(3)
    ke = metgrid%fields(metgrid%qv_no)%mem_end(3)

    nx = ie-is+1
    ny = je-js+1
    nz = ke-ks+1

    ALLOCATE(qvs(is:ie,js:je,ks:ke), STAT = istatus)

    CALL getqvs(nx,ny,nz,1,nx,1,ny,1,nz,                    & !dimensions
                metgrid%fields(metgrid%pres_no)%rdata_arr,  & ! pressure
                metgrid%fields(metgrid%tt_no)%rdata_arr,    & ! temperature
                qvs)
    !WRITE(6,*) ks,ke, js,je,is,ie,metgrid%qv_no,metgrid%pres_no,metgrid%tt_no
    DO k = ks,ke
      DO j = js, je
        DO i = is, ie
          !WRITE(6,*) i,j,k,qvs(i,j,k)
          metgrid%fields(metgrid%qv_no)%rdata_arr(i,j,k) = 100          &
                       * metgrid%fields(metgrid%qv_no)%rdata_arr(i,j,k) &
                       / qvs(i,j,k)
        END DO
      END DO
    END DO

    metgrid%fields(metgrid%qv_no)%fieldname = 'RH'
    metgrid%fields(metgrid%qv_no)%units     = '%'
    metgrid%fields(metgrid%qv_no)%descr     = 'Relative Humidity'

    DEALLOCATE(qvs)

    ! Remove FLAG_QV
    k = metgrid%nflags + 1
    DO i = 1, metgrid%nflags
      IF ( metgrid%flags(i) == 'FLAG_QV') THEN
        k = i
      END IF

      IF (i > k) THEN
        metgrid%flags(k) = metgrid%flags(i)
        k = k+1
      END IF
    END DO

    RETURN
  END SUBROUTINE metgrid_qv_to_rh

END MODULE module_metgrid
