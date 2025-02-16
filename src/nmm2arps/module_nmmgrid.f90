MODULE module_nmm2arps_nmmgrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains data for WRF NMM grid
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE module_nmm_input
  USE wrf_llxy_module
  USE wrf_parallel_module

  TYPE TYPE_FIELD
    CHARACTER (LEN=128) :: fieldname
    CHARACTER (LEN=1)   :: stagger
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rdata_arr
  END TYPE TYPE_FIELD

  TYPE type_nmmgrid

    INTEGER :: map_proj
    REAL    :: cen_lat,cen_lon
    INTEGER :: iswater, isice, isurban, isoilwater

    INTEGER :: west_east, south_north, bottom_top_stag  ! stag
    INTEGER :: soil_layer_stag

    REAL    :: dx, dy
    INTEGER :: ips, ipe, jps, jpe
    INTEGER :: ims, ime, jms, jme
    INTEGER :: kps, kpe, kpse, kzs, kze, kzse

    INTEGER :: mp_physics, num_moist, num_qq
    REAL    :: xtime

    TYPE (type_field) :: glat, glon
    TYPE (type_field) :: pres, ww
    TYPE (type_field) :: tt
    TYPE (type_field) :: uu, vv
    TYPE (type_field) :: qv   ! Read in as specific humidity, but changed to RHstar after adjustdata
    TYPE (type_field), ALLOCATABLE :: qscalar(:)
!    TYPE (type_field) :: tke

    TYPE (type_field) :: sldpth

    TYPE (type_field) :: fis  ! Read in as geopotential, but changed to height after calling of adjustdata
    TYPE (type_field) :: ths  ! Read in as potentail temp, but changed to temperaure after calling of adjustdata
    TYPE (type_field) :: qvs
    TYPE (type_field) :: psfc ! may not present in the data file

    TYPE (type_field) :: tsoil, qsoil, canpwet
    TYPE (type_field) :: snowh, vegfrc, z0
    TYPE (type_field) :: acprec, cuprec, th2, u10, v10
    TYPE (type_field) :: tshltr,qshltr,pshltr
!    TYPE (type_field) :: wspd10max_ext,w_up_max_ext,w_dn_max_ext,refd_max_ext
!    TYPE (type_field) :: up_heli_max_ext,grpl_max_ext
    TYPE (type_field) :: w_up_max, w_dn_max, wspd10max, u_10m_max, v_10m_max
    TYPE (type_field) :: up_heli_max, refd_max
    TYPE (type_field) :: rlwdn,rswdn
    LOGICAL :: SPRING2011

    INTEGER, ALLOCATABLE :: isltyp(:,:), ivgtyp(:,:)

    ! derived fields
    TYPE (type_field) :: algp
    TYPE (type_field) :: pp, pt, rho
    TYPE (type_field) :: zp, hgt

    INTEGER :: p_qc, p_qr, p_qi, p_qs, p_qg, p_qh
    INTEGER :: p_nc, p_nr, p_ni, p_ns, p_ng, p_nh, p_nn
    INTEGER :: p_zr, p_zi, p_zs, p_zg, p_zh
    INTEGER :: p_rim

  END TYPE type_nmmgrid

  TYPE (type_nmmgrid) :: nmmgrid

  CONTAINS

  SUBROUTINE nmmgrid_init(IAMROOT,filename,io_form,multifile,           &
                          ncompressx,ncompressy,numdigits,lvldbg,istatus)

!#######################################################################

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: IAMROOT
    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: filename
    INTEGER, INTENT(IN)  :: io_form
    LOGICAL, INTENT(IN)  :: multifile
    INTEGER, INTENT(IN)  :: ncompressx, ncompressy, numdigits
    INTEGER, INTENT(IN)  :: lvldbg

    INTEGER, INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
    INTEGER :: nx,ny,nz,nzsoil
    INTEGER :: iproj
    REAL    :: trlat1,trlat2,trlon
    REAL    :: ctrlat,ctrlon
    REAL    :: dx,dy
    INTEGER :: p1s,p1e,p2s,p2e
    INTEGER :: iswater, isice, isurban, isoilwater

    REAL    :: phi, lambda

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
    nmmgrid%SPRING2011 = .FALSE.

    IF ( .NOT. input_initialized) THEN
      CALL nmm_input_init(IAMROOT,io_form,numdigits,                    &
                          multifile,ncompressx,ncompressy,istatus)
    END IF

    CALL nmm_input_open(filename,lvldbg,istatus)

    CALL nmm_input_get_meta(nx,ny,nz,nzsoil,nmmgrid%mp_physics,         &
                            iproj,ctrlat,ctrlon,dx,dy,                  &
                            p1s,p1e,p2s,p2e,                            &
                            iswater, isice, isurban, isoilwater,istatus)

    CALL nmm_input_close(istatus)

!-----------------------------------------------------------------------

    nmmgrid%west_east   = nx-1
    nmmgrid%south_north = ny-1
    nmmgrid%bottom_top_stag  = nz
    nmmgrid%soil_layer_stag  = nzsoil

    CALL parallel_get_tile_dims(nmmgrid%west_east,nmmgrid%south_north)

    IF (multifile) THEN
      nmmgrid%ips  = p1s
      nmmgrid%ipe  = p1e
      nmmgrid%jps  = p2s
      nmmgrid%jpe  = p2e
    ELSE
      nmmgrid%ips  = my_minx
      nmmgrid%ipe  = my_maxx
      nmmgrid%jps  = my_miny
      nmmgrid%jpe  = my_maxy
    END IF
    nmmgrid%kps    = 1               ! vertical dimension
    nmmgrid%kpe    = nz-1
    nmmgrid%kpse   = nz
    nmmgrid%kzs    = 1               ! Vertical soil dimension
    nmmgrid%kze    = nzsoil-1
    nmmgrid%kzse   = nzsoil

    nmmgrid%map_proj = iproj
    nmmgrid%cen_lat  = ctrlat
    nmmgrid%cen_lon  = ctrlon

    nmmgrid%dx        = dx
    nmmgrid%dy        = dy
    nmmgrid%iswater   = iswater
    nmmgrid%isice     = isice
    nmmgrid%isurban   = isurban
    nmmgrid%isoilwater= isoilwater

!------------------------------------------------------------------------

    nmmgrid%p_qc = 0
    nmmgrid%p_qr = 0
    nmmgrid%p_qi = 0
    nmmgrid%p_qs = 0
    nmmgrid%p_qg = 0
    nmmgrid%p_qh = 0
    nmmgrid%num_qq = 0

    nmmgrid%p_nc = 0
    nmmgrid%p_nr = 0
    nmmgrid%p_ni = 0
    nmmgrid%p_ns = 0
    nmmgrid%p_ng = 0
    nmmgrid%p_nh = 0
    nmmgrid%p_nn = 0
    nmmgrid%num_moist = 0

    nmmgrid%p_zr = 0
    nmmgrid%p_zi = 0
    nmmgrid%p_zs = 0
    nmmgrid%p_zg = 0
    nmmgrid%p_zh = 0

    SELECT CASE (nmmgrid%mp_physics)
    CASE (1, 3)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%num_qq = 2
      nmmgrid%num_moist = 2
    CASE (2, 6, 7)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%p_qg = 5
      nmmgrid%num_qq = 5
      nmmgrid%num_moist = 5
    CASE (4, 13, 85)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%num_qq = 4
      nmmgrid%num_moist = 4
    CASE (5)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qs = 3
      nmmgrid%num_qq = 3
      nmmgrid%p_rim  = 4
      nmmgrid%num_moist = 4
    CASE (8)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%p_qg = 5
      nmmgrid%num_qq = 5

      nmmgrid%p_ni = 6
      nmmgrid%p_nr = 7
      nmmgrid%num_moist = 7
    CASE (9)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%p_qg = 5
      nmmgrid%p_qh = 6
      nmmgrid%num_qq = 6

      nmmgrid%p_nc = 7
      nmmgrid%p_nr = 8
      nmmgrid%p_ni = 9
      nmmgrid%p_ns = 10
      nmmgrid%p_ng = 11
      nmmgrid%p_nh = 12
      nmmgrid%num_moist = 12
    CASE (10)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%p_qg = 5
      nmmgrid%num_qq = 5

      nmmgrid%p_nr = 6
      nmmgrid%p_ni = 7
      nmmgrid%p_ns = 8
      nmmgrid%p_ng = 9
      nmmgrid%num_moist = 9
    CASE (14)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%num_qq = 4

      nmmgrid%p_nc = 5
      nmmgrid%p_nr = 6
      nmmgrid%p_nn = 7
      nmmgrid%num_moist = 7
    CASE (16)
      nmmgrid%p_qc = 1
      nmmgrid%p_qr = 2
      nmmgrid%p_qi = 3
      nmmgrid%p_qs = 4
      nmmgrid%p_qg = 5
      nmmgrid%num_qq = 5

      nmmgrid%p_nc = 6
      nmmgrid%p_nr = 7
      nmmgrid%p_nn = 8
      nmmgrid%num_moist = 8
    CASE DEFAULT
      WRITE(6,'(1x,a,I3)') 'ERROR: Not supported microphysics option: ',nmmgrid%mp_physics
      CALL arpsstop('Unknown mp_physics.',1)
    END SELECT

!-----------------------------------------------------------------------
    phi    = dy*real(nmmgrid%south_north-1)/2.  ! Based on the statements in geogrid
    lambda = dx*real(nmmgrid%west_east-1)       ! file gridinfo_module.F.

    CALL set_domain_projection(iproj, MISSING_DATA, MISSING_DATA, MISSING_DATA, &
                           dx,dy, phi, lambda,                          &
                           nmmgrid%west_east, nmmgrid%south_north,      &
                           nmmgrid%west_east/2.,nmmgrid%south_north/2., &
                           ctrlat, ctrlon, MISSING_DATA, MISSING_DATA )

    RETURN
  END SUBROUTINE nmmgrid_init

  SUBROUTINE nmmgrid_alloc(IAMROOT,lvldbg,istatus)

!#######################################################################
!
!  This one must be called after HALO_WIDTH is determined.
!
    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: IAMROOT
    INTEGER, INTENT(IN)  :: lvldbg

    INTEGER, INTENT(OUT) :: istatus
!-----------------------------------------------------------------------

    INTEGER :: lh_mult, rh_mult, bh_mult, th_mult
    INTEGER :: ims, ime, jms, jme

    INTEGER :: nz

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nz = nmmgrid%bottom_top_stag

!------------------------------------------------------------------------

    ! Compute multipliers for halo width; these must be 0/1
    lh_mult = 1
    if (my_x == 0) lh_mult = 0

    rh_mult = 1
    if (my_x == (nxprocs-1)) rh_mult = 0

    bh_mult = 1
    if (my_y == 0) bh_mult = 0

    th_mult = 1
    if (my_y == (nyprocs-1)) th_mult = 0

    nmmgrid%ims = nmmgrid%ips - HALO_WIDTH*lh_mult
    nmmgrid%ime = nmmgrid%ipe + HALO_WIDTH*rh_mult
    nmmgrid%jms = nmmgrid%jps - HALO_WIDTH*bh_mult
    nmmgrid%jme = nmmgrid%jpe + HALO_WIDTH*th_mult

    ims = nmmgrid%ims     ! used local for allocation
    ime = nmmgrid%ime
    jms = nmmgrid%jms
    jme = nmmgrid%jme

    ALLOCATE(nmmgrid%glat%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%glat%fieldname = 'GLAT'
    nmmgrid%glat%stagger   = ' '

    ALLOCATE(nmmgrid%glon%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%glon%fieldname = 'GLON'
    nmmgrid%glon%stagger   = ' '

    ALLOCATE(nmmgrid%pres%rdata_arr(ims:ime,jms:jme,1:nz), STAT = istatus)
    nmmgrid%pres%fieldname = 'PINT'
    nmmgrid%pres%stagger   = 'Z'

    ALLOCATE(nmmgrid%tt%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%tt%fieldname = 'T'
    nmmgrid%tt%stagger   = ' '

    ALLOCATE(nmmgrid%uu%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%uu%fieldname = 'U'
    nmmgrid%uu%stagger   = ' '

    ALLOCATE(nmmgrid%vv%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%vv%fieldname = 'V'
    nmmgrid%vv%stagger   = ' '

    ALLOCATE(nmmgrid%ww%rdata_arr(ims:ime,jms:jme,1:nz), STAT = istatus)
    nmmgrid%ww%fieldname = 'W'
    nmmgrid%ww%stagger   = 'Z'

    ALLOCATE(nmmgrid%qv%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%qv%fieldname = 'Q'
    nmmgrid%qv%stagger   = ' '

    ALLOCATE(nmmgrid%qscalar(nmmgrid%num_moist), STAT = istatus)

    IF (nmmgrid%p_qc > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_qc)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_qc)%fieldname = 'QCLOUD'
      nmmgrid%qscalar(nmmgrid%p_qc)%stagger   = ' '
    END IF

    IF (nmmgrid%p_qr > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_qr)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_qr)%fieldname = 'QRAIN'
      nmmgrid%qscalar(nmmgrid%p_qr)%stagger   = ' '
    END IF

    IF (nmmgrid%p_qi > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_qi)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_qi)%fieldname = 'QICE'
      nmmgrid%qscalar(nmmgrid%p_qi)%stagger   = ' '
    END IF

    IF (nmmgrid%p_qs > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_qs)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_qs)%fieldname = 'QSNOW'
      nmmgrid%qscalar(nmmgrid%p_qs)%stagger   = ' '
    END IF

    IF (nmmgrid%p_qg > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_qg)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_qg)%fieldname = 'QGRAUP'
      nmmgrid%qscalar(nmmgrid%p_qg)%stagger   = ' '
    END IF

    IF (nmmgrid%p_qh > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_qh)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_qh)%fieldname = 'QHAIL'
      nmmgrid%qscalar(nmmgrid%p_qh)%stagger   = ' '
    END IF

    IF (nmmgrid%p_nc > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_nc)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_nc)%fieldname = 'QNCLOUD'
      nmmgrid%qscalar(nmmgrid%p_nc)%stagger   = ' '
    END IF

    IF (nmmgrid%p_nr > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_nr)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_nr)%fieldname = 'QNR'
      nmmgrid%qscalar(nmmgrid%p_nr)%stagger   = ' '
    END IF

    IF (nmmgrid%p_ni > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_ni)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_ni)%fieldname = 'QNI'
      nmmgrid%qscalar(nmmgrid%p_ni)%stagger   = ' '
    END IF

    IF (nmmgrid%p_ns > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_ns)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_ns)%fieldname = 'QNS'
      nmmgrid%qscalar(nmmgrid%p_ns)%stagger   = ' '
    END IF

    IF (nmmgrid%p_ng > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_ng)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_ng)%fieldname = 'QNG'
      nmmgrid%qscalar(nmmgrid%p_ng)%stagger   = ' '
    END IF

    IF (nmmgrid%p_nh > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_nh)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_nh)%fieldname = 'QNHAIL'
      nmmgrid%qscalar(nmmgrid%p_nh)%stagger   = ' '
    END IF

    IF (nmmgrid%p_nn > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_nn)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_nn)%fieldname = 'QNCCN'
      nmmgrid%qscalar(nmmgrid%p_nn)%stagger   = ' '
    END IF

    IF (nmmgrid%p_rim > 0) THEN
      ALLOCATE(nmmgrid%qscalar(nmmgrid%p_rim)%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
      nmmgrid%qscalar(nmmgrid%p_rim)%fieldname = 'F_RIMEF'
      nmmgrid%qscalar(nmmgrid%p_rim)%stagger   = ' '
    END IF

!    ALLOCATE(nmmgrid%tke%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
!    nmmgrid%tke%fieldname = 'TKE_PBL' !????
!    nmmgrid%tke%stagger   = ' '

    ALLOCATE(nmmgrid%sldpth%rdata_arr(1:nz-1,1,1), STAT = istatus)
    nmmgrid%sldpth%fieldname = 'SLDPTH'
    nmmgrid%sldpth%stagger   = ' '

    ALLOCATE(nmmgrid%fis%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%fis%fieldname = 'FIS'
    nmmgrid%fis%stagger   = ' '

    ALLOCATE(nmmgrid%ths%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%ths%fieldname = 'THS'
    nmmgrid%ths%stagger   = ' '

    ALLOCATE(nmmgrid%qvs%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%qvs%fieldname = 'QS'
    nmmgrid%qvs%stagger   = ' '

    ALLOCATE(nmmgrid%psfc%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%psfc%fieldname = 'PSFC'
    nmmgrid%psfc%stagger   = ' '

    ALLOCATE(nmmgrid%canpwet%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%canpwet%fieldname = 'CMC'   ! An candidate is CANWET
    nmmgrid%canpwet%stagger   = ' '

    ALLOCATE(nmmgrid%snowh%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%snowh%fieldname = 'SNOWH'   ! or SNOW?
    nmmgrid%snowh%stagger   = ' '

    ALLOCATE(nmmgrid%isltyp(ims:ime,jms:jme), STAT = istatus)
    ALLOCATE(nmmgrid%ivgtyp(ims:ime,jms:jme), STAT = istatus)

    ALLOCATE(nmmgrid%vegfrc%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%vegfrc%fieldname = 'VEGFRC'
    nmmgrid%vegfrc%stagger   = ' '

    ALLOCATE(nmmgrid%z0%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%z0%fieldname = 'Z0'
    nmmgrid%z0%stagger   = ' '

    ALLOCATE(nmmgrid%acprec%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%acprec%fieldname = 'ACPREC'
    nmmgrid%acprec%stagger   = ' '

    ALLOCATE(nmmgrid%cuprec%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%cuprec%fieldname = 'CUPREC'
    nmmgrid%cuprec%stagger   = ' '

    ALLOCATE(nmmgrid%th2%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%th2%fieldname = 'TH2'
    nmmgrid%th2%stagger   = ' '

    ALLOCATE(nmmgrid%tshltr%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%tshltr%fieldname = 'TSHLTR'
    nmmgrid%tshltr%stagger   = ' '

    ALLOCATE(nmmgrid%pshltr%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%pshltr%fieldname = 'PSHLTR'
    nmmgrid%pshltr%stagger   = ' '

    ALLOCATE(nmmgrid%qshltr%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%qshltr%fieldname = 'QSHLTR'
    nmmgrid%qshltr%stagger   = ' '

    ALLOCATE(nmmgrid%u10%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%u10%fieldname = 'U10'
    nmmgrid%u10%stagger   = ' '

    ALLOCATE(nmmgrid%v10%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
    nmmgrid%v10%fieldname = 'V10'
    nmmgrid%v10%stagger   = ' '

    IF (nmmgrid%SPRING2011) THEN

      ALLOCATE(nmmgrid%wspd10max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
!      nmmgrid%wspd10max%fieldname = 'WSPD10MAX'
      nmmgrid%wspd10max%fieldname = 'MAX10MW'
      nmmgrid%wspd10max%stagger   = ' '

      ALLOCATE(nmmgrid%u_10m_max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%u_10m_max%fieldname = 'MAX10U'
      nmmgrid%u_10m_max%stagger   = ' '

      ALLOCATE(nmmgrid%v_10m_max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%v_10m_max%fieldname = 'MAX10V'
      nmmgrid%v_10m_max%stagger   = ' '

      ALLOCATE(nmmgrid%w_up_max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%w_up_max%fieldname = 'MAXUPDR'
      nmmgrid%w_up_max%stagger   = ' '

      ALLOCATE(nmmgrid%w_dn_max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%w_dn_max%fieldname = 'MAXDNDR'
      nmmgrid%w_dn_max%stagger   = ' '

      ALLOCATE(nmmgrid%refd_max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%refd_max%fieldname = 'MAXDBZ'
      nmmgrid%refd_max%stagger   = ' '

      ALLOCATE(nmmgrid%up_heli_max%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%up_heli_max%fieldname = 'MAXHLCY'
      nmmgrid%up_heli_max%stagger   = ' '

      ALLOCATE(nmmgrid%rlwdn%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%rlwdn%fieldname = 'RLWIN'
      nmmgrid%rlwdn%stagger   = ' '

      ALLOCATE(nmmgrid%rswdn%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      nmmgrid%rswdn%fieldname = 'RSWIN'
      nmmgrid%rswdn%stagger   = ' '

      !ALLOCATE(nmmgrid%grpl_max_ext%rdata_arr(ims:ime,jms:jme,1), STAT = istatus)
      !nmmgrid%grpl_max_ext%fieldname = 'GRPL_MAX'
      !nmmgrid%grpl_max_ext%stagger   = ' '
    END IF

    !
    ! Derived Fields
    !
    ALLOCATE(nmmgrid%pt%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%pt%fieldname = 'PT'
    nmmgrid%pt%stagger   = ' '

    ALLOCATE(nmmgrid%pp%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%pp%fieldname = 'PP'
    nmmgrid%pp%stagger   = ' '

    ALLOCATE(nmmgrid%rho%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%rho%fieldname = 'RHO'
    nmmgrid%rho%stagger   = ' '

    ALLOCATE(nmmgrid%algp%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%algp%fieldname = 'ALGP'
    nmmgrid%algp%stagger   = ' '

    ALLOCATE(nmmgrid%hgt%rdata_arr(ims:ime,jms:jme,1:nz-1), STAT = istatus)
    nmmgrid%hgt%fieldname = 'HGT'
    nmmgrid%hgt%stagger   = ' '

    ALLOCATE(nmmgrid%zp%rdata_arr(ims:ime,jms:jme,1:nz), STAT = istatus)
    nmmgrid%zp%fieldname = 'ZP'
    nmmgrid%zp%stagger   = ' '

    RETURN
  END SUBROUTINE nmmgrid_alloc

  SUBROUTINE nmmgrid_dealloc(istatus)

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    CALL nmm_input_shutdown( istatus )

    DEALLOCATE(nmmgrid%glat%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%glon%rdata_arr, STAT = istatus)

    DEALLOCATE(nmmgrid%pres%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%tt%rdata_arr,   STAT = istatus)
    DEALLOCATE(nmmgrid%uu%rdata_arr,   STAT = istatus)
    DEALLOCATE(nmmgrid%vv%rdata_arr,   STAT = istatus)
    DEALLOCATE(nmmgrid%ww%rdata_arr,   STAT = istatus)
    DEALLOCATE(nmmgrid%qv%rdata_arr,   STAT = istatus)

    DO nq = 1, nmmgrid%num_moist
      DEALLOCATE(nmmgrid%qscalar(nq)%rdata_arr,   STAT = istatus)
    END DO
    DEALLOCATE(nmmgrid%qscalar,   STAT = istatus)

!    DEALLOCATE(nmmgrid%tke%rdata_arr,   STAT = istatus)

    DEALLOCATE(nmmgrid%fis%rdata_arr,  STAT = istatus)
    DEALLOCATE(nmmgrid%ths%rdata_arr,  STAT = istatus)
    DEALLOCATE(nmmgrid%qvs%rdata_arr,  STAT = istatus)
    DEALLOCATE(nmmgrid%psfc%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%sldpth%rdata_arr, STAT = istatus)

    DEALLOCATE(nmmgrid%tsoil%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%qsoil%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%canpwet%rdata_arr, STAT = istatus)

    DEALLOCATE(nmmgrid%snowh%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%isltyp, STAT = istatus)
    DEALLOCATE(nmmgrid%ivgtyp, STAT = istatus)
    DEALLOCATE(nmmgrid%vegfrc%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%z0%rdata_arr, STAT = istatus)

    DEALLOCATE(nmmgrid%acprec%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%cuprec%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%th2%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%tshltr%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%pshltr%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%qshltr%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%u10%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%v10%rdata_arr, STAT = istatus)

    IF (nmmgrid%SPRING2011) THEN
      DEALLOCATE(nmmgrid%wspd10max%rdata_arr, STAT = istatus)

      DEALLOCATE(nmmgrid%u_10m_max%rdata_arr,   STAT = istatus)
      DEALLOCATE(nmmgrid%v_10m_max%rdata_arr,   STAT = istatus)
      DEALLOCATE(nmmgrid%w_up_max%rdata_arr,    STAT = istatus)
      DEALLOCATE(nmmgrid%w_dn_max%rdata_arr,    STAT = istatus)
      DEALLOCATE(nmmgrid%refd_max%rdata_arr,    STAT = istatus)
      DEALLOCATE(nmmgrid%up_heli_max%rdata_arr, STAT = istatus)
      DEALLOCATE(nmmgrid%rlwdn%rdata_arr, STAT = istatus)
      DEALLOCATE(nmmgrid%rswdn%rdata_arr, STAT = istatus)
      !DEALLOCATE(nmmgrid%grpl_max_ext%rdata_arr, STAT = istatus)
    END IF

    DEALLOCATE(nmmgrid%pt%rdata_arr,   STAT = istatus)
    DEALLOCATE(nmmgrid%pp%rdata_arr,   STAT = istatus)
    DEALLOCATE(nmmgrid%rho%rdata_arr,  STAT = istatus)
    DEALLOCATE(nmmgrid%algp%rdata_arr, STAT = istatus)
    DEALLOCATE(nmmgrid%hgt%rdata_arr,  STAT = istatus)
    DEALLOCATE(nmmgrid%zp%rdata_arr,   STAT = istatus)

    RETURN
  END SUBROUTINE nmmgrid_dealloc

  SUBROUTINE nmmgrid_getdata( datestr, irain, initsec, itime, abstimei, &
                              dir_extp, grid_id, dbglvl, istatus)

!#######################################################################

    IMPLICIT NONE

    CHARACTER(LEN=19), INTENT(IN) :: datestr
    LOGICAL, INTENT(IN)  :: irain   ! 1 - piecewise rain
                                    ! 0 - Accumulated rain just as it is in WRF filename
    INTEGER, INTENT(IN)  :: initsec ! initial absolute time
    INTEGER, INTENT(IN)  :: itime   ! current frame number
    INTEGER, INTENT(IN)  :: abstimei
    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: dir_extp
    INTEGER, INTENT(IN)  :: grid_id

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
    CHARACTER(LEN=128), DIMENSION(MAXDIMENSIONS) :: dimnames
    INTEGER,            DIMENSION(MAXDIMENSIONS) :: mStart
    INTEGER,            DIMENSION(MAXDIMENSIONS) :: mEnd
    INTEGER,            DIMENSION(MAXDIMENSIONS) :: pStart
    INTEGER,            DIMENSION(MAXDIMENSIONS) :: pEnd
    INTEGER,            DIMENSION(MAXDIMENSIONS) :: dStart
    INTEGER,            DIMENSION(MAXDIMENSIONS) :: dEnd

    INTEGER :: i,j,k,nq, nzsoil
    REAL    :: temlocal(1,1,1)

    INTEGER :: abstimec, abstimep
    INTEGER :: year, month, day, hour, minute, second
    CHARACTER(LEN=1) :: ach
    CHARACTER(LEN=MAXFILELEN) :: filename

    REAL, ALLOCATABLE :: raing(:,:), rainc(:,:)

    DOUBLE PRECISION  :: pi, r2d

    REAL :: testa, testb, testc, testd

    CHARACTER(LEN=40) :: tmpwrftime

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    dimnames(1) = dimwename
    dimnames(2) = dimsnname
    dimnames(3) = dimbtname

    mStart(1) = nmmgrid%ims; mStart(2) = nmmgrid%jms; mStart(3) = nmmgrid%kps
    mEnd(1)   = nmmgrid%ime; mEnd(2)   = nmmgrid%jme; mEnd(3)   = nmmgrid%kpe
    pStart(1) = nmmgrid%ips; pStart(2) = nmmgrid%jps; pStart(3) = nmmgrid%kps
    pEnd(1)   = nmmgrid%ipe; pEnd(2)   = nmmgrid%jpe; pEnd(3)   = nmmgrid%kpe
    dStart(1) = 1;  dEnd(1) = nmmgrid%west_east;      dStart(3) = nmmgrid%kps
    dStart(2) = 1;  dEnd(2) = nmmgrid%south_north;    dEnd(3)   = nmmgrid%kpe

!--------------------- LAT/LON ----------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%glat%fieldname,            &
                             nmmgrid%glat%stagger,dimnames,             &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%glat%rdata_arr,        &
                             dbglvl,istatus)

    CALL nmm_input_get_3dvar(datestr,nmmgrid%glon%fieldname,            &
                             nmmgrid%glon%stagger,dimnames,             &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%glon%rdata_arr,        &
                             dbglvl,istatus)
    pi = ACOS(-1.0)
    r2d = 180./pi
    DO j = nmmgrid%jms, nmmgrid%jme
      DO i = nmmgrid%ims, nmmgrid%ime
        nmmgrid%glat%rdata_arr(i,j,1) = nmmgrid%glat%rdata_arr(i,j,1)*r2d
        nmmgrid%glon%rdata_arr(i,j,1) = nmmgrid%glon%rdata_arr(i,j,1)*r2d
      END DO
    END DO

!-- TEST --
!    DO j = 1, 3
!      DO i = 1,3
!        CALL wps_lltoxy(nmmgrid%glat%rdata_arr(i,j,1),nmmgrid%glon%rdata_arr(i,j,1),testa,testb,HH)
!        write(0,*) 'HH',i,j, testa, testb
!!        CALL wps_xytoll(testa, testb,testc,testd,HH)
!!        write(0,*) 'HH',nmmgrid%glat%rdata_arr(i,j),nmmgrid%glon%rdata_arr(i,j), testc, testd
!
!        CALL wps_lltoxy(nmmgrid%glat%rdata_arr(i,j,1),nmmgrid%glon%rdata_arr(i,j,1),testa,testb,VV)
!        write(0,*) 'VV',i,j, testa, testb
! !       CALL wps_xytoll(testa, testb,testc,testd,VV)
! !       write(0,*) 'vv',nmmgrid%glat%rdata_arr(i,j),nmmgrid%glon%rdata_arr(i,j), testc, testd
!      end DO
!    end DO
!
!call arpsstop(' ',0)
!--------------------- FIS ---------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%fis%fieldname,            &
                             nmmgrid%fis%stagger,dimnames,             &
                             mStart,mEnd,pStart,pEnd,                  &
                             dStart,dEnd,nmmgrid%fis%rdata_arr,        &
                             dbglvl,istatus)

!--------------------- THS ----------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%ths%fieldname,            &
                             nmmgrid%ths%stagger,dimnames,             &
                             mStart,mEnd,pStart,pEnd,                  &
                             dStart,dEnd,nmmgrid%ths%rdata_arr,        &
                             dbglvl,istatus)

!--------------------- QVS ----------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%qvs%fieldname,            &
                             nmmgrid%qvs%stagger,dimnames,             &
                             mStart,mEnd,pStart,pEnd,                  &
                             dStart,dEnd,nmmgrid%qvs%rdata_arr,        &
                             dbglvl,istatus)

!--------------------- PSFC ----------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%psfc%fieldname,           &
                             nmmgrid%psfc%stagger,dimnames,            &
                             mStart,mEnd,pStart,pEnd,                  &
                             dStart,dEnd,nmmgrid%psfc%rdata_arr,       &
                             dbglvl,istatus)

!--------------------- Pressure ----------------------------------------
    dimnames(3) = dimbtsname
    mEnd(3) = nmmgrid%kpse; pEnd(3) = nmmgrid%kpse; dEnd(3) = nmmgrid%kpse

    CALL nmm_input_get_3dvar(datestr,nmmgrid%pres%fieldname,            &
                             nmmgrid%pres%stagger,dimnames,             &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%pres%rdata_arr,        &
                             dbglvl,istatus)

!--------------------- Temperature -------------------------------------
    dimnames(3) = dimbtname
    mEnd(3) = nmmgrid%kpe; pEnd(3) = nmmgrid%kpe; dEnd(3) = nmmgrid%kpe

    CALL nmm_input_get_3dvar(datestr,nmmgrid%tt%fieldname,              &
                             nmmgrid%tt%stagger,dimnames,               &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%tt%rdata_arr,          &
                             dbglvl,istatus)

!--------------------- U-Wind -------------------------------------
    dimnames(3) = dimbtname
    mEnd(3) = nmmgrid%kpe; pEnd(3) = nmmgrid%kpe; dEnd(3) = nmmgrid%kpe

    CALL nmm_input_get_3dvar(datestr,nmmgrid%uu%fieldname,              &
                             nmmgrid%uu%stagger,dimnames,               &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%uu%rdata_arr,          &
                             dbglvl,istatus)

!--------------------- V-Wind -------------------------------------
    dimnames(3) = dimbtname
    mEnd(3) = nmmgrid%kpe; pEnd(3) = nmmgrid%kpe; dEnd(3) = nmmgrid%kpe

    CALL nmm_input_get_3dvar(datestr,nmmgrid%vv%fieldname,              &
                             nmmgrid%vv%stagger,dimnames,               &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%vv%rdata_arr,          &
                             dbglvl,istatus)

!--------------------- W-Wind -------------------------------------
    dimnames(3) = dimbtsname
    mEnd(3) = nmmgrid%kpse; pEnd(3) = nmmgrid%kpse; dEnd(3) = nmmgrid%kpse

    CALL nmm_input_get_3dvar(datestr,nmmgrid%ww%fieldname,              &
                             nmmgrid%ww%stagger,dimnames,               &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%ww%rdata_arr,          &
                             dbglvl,istatus)

!--------------------- Q ------------------------------------------
    dimnames(3) = dimbtname
    mEnd(3) = nmmgrid%kpe; pEnd(3) = nmmgrid%kpe; dEnd(3) = nmmgrid%kpe

    CALL nmm_input_get_3dvar(datestr,nmmgrid%qv%fieldname,              &
                             nmmgrid%qv%stagger,dimnames,               &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%qv%rdata_arr,          &
                             dbglvl,istatus)

!--------------------- QSCALAR ---------------------------------------
    DO nq = 1, nmmgrid%num_moist
      dimnames(3) = dimbtname
      mEnd(3) = nmmgrid%kpe; pEnd(3) = nmmgrid%kpe; dEnd(3) = nmmgrid%kpe

      CALL nmm_input_get_3dvar(datestr,nmmgrid%qscalar(nq)%fieldname,   &
                             nmmgrid%qscalar(nq)%stagger,dimnames,      &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%qscalar(nq)%rdata_arr, &
                             dbglvl,istatus)
    END DO

!--------------------- TKE ------------------------------------------
!    dimnames(3) = dimbtname
!    mEnd(3) = nmmgrid%kpe; pEnd(3) = nmmgrid%kpe; dEnd(3) = nmmgrid%kpe
!
!    CALL nmm_input_get_3dvar(datestr,nmmgrid%tke%fieldname,             &
!                             nmmgrid%tke%stagger,dimnames,              &
!                             mStart,mEnd,pStart,pEnd,                   &
!                             dStart,dEnd,nmmgrid%tke%rdata_arr,         &
!                             dbglvl,istatus)
!
!--------------------- SNOWH ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%snowh%fieldname,           &
                             nmmgrid%snowh%stagger,dimnames,            &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%snowh%rdata_arr,       &
                             dbglvl,istatus)

!--------------------- ISLTYP ------------------------------------------

    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvari(datestr,'ISLTYP',' ',dimnames,            &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%isltyp,                &
                             dbglvl,istatus)

!--------------------- IVGTYP ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvari(datestr,'IVGTYP',' ',dimnames,            &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%ivgtyp,                &
                             dbglvl,istatus)

!--------------------- VEGFRC ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%vegfrc%fieldname,          &
                             nmmgrid%vegfrc%stagger,dimnames,           &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%vegfrc%rdata_arr,      &
                             dbglvl,istatus)

!--------------------- Z0    ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%z0%fieldname,              &
                             nmmgrid%z0%stagger,dimnames,               &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%z0%rdata_arr,          &
                             dbglvl,istatus)

!--------------------- ACPREC ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    IF (irain) THEN

      ALLOCATE(raing(nmmgrid%ims:nmmgrid%ime,nmmgrid%jms:nmmgrid%jme), STAT = istatus)
      ALLOCATE(rainc(nmmgrid%ims:nmmgrid%ime,nmmgrid%jms:nmmgrid%jme), STAT = istatus)

      IF (itime > 1) THEN    ! This is not the first file processed
        DO j = nmmgrid%jms, nmmgrid%jme
          DO i = nmmgrid%ims, nmmgrid%ime
            raing(i,j) = nmmgrid%acprec%rdata_arr(i,j,1)/1000.0   ! it is already converted to mm
            rainc(i,j) = nmmgrid%cuprec%rdata_arr(i,j,1)/1000.0   ! So should be converted back
          END DO
        END DO
      ELSE
        READ(datestr,'(I4.4,5(a,I2.2))')      &
                     year,ach,month,ach,day,ach,hour,ach,minute,ach,second
        IF (year < 1960)  year =  1960  ! maybe ideal case
        CALL ctim2abss(year,month,day,hour,minute,second,abstimec)

        IF (abstimec > initsec) THEN  ! try to read precipitation at previous time level
          abstimep = abstimec - abstimei
          CALL abss2ctim(abstimep,year,month,day,hour,minute,second)

          WRITE(tmpwrftime,'(a,I2.2,a,I4.4,5(a,I2.2))')                 &
                 'wrfout_d',grid_id,'_',                                &
                 year,'-',month,'-',day,'_',hour,':',minute,':',second

          IF (dir_extp == '<dir>') THEN
            WRITE(filename,'(3a)') TRIM(tmpwrftime),'/',TRIM(tmpwrftime)
          ELSE
            WRITE(filename,'(2a)') TRIM(dir_extp),TRIM(tmpwrftime)
          END IF

          CALL nmm_get_rain(filename,datestr,                           &
                     nmmgrid%acprec%fieldname,nmmgrid%cuprec%fieldname, &
                          nmmgrid%acprec%stagger, dimnames,             &
                          mStart, mEnd, pStart, pEnd, raing, rainc,     &
                          dStart,dEnd, dbglvl,istatus )

          IF (istatus /= 0) THEN
            raing(:,:) = 0.0
            rainc(:,:) = 0.0
          END IF
        ELSE
          raing(:,:) = 0.0
          rainc(:,:) = 0.0
        END IF
      END IF

    END IF

    CALL nmm_input_get_3dvar(datestr,nmmgrid%acprec%fieldname,          &
                             nmmgrid%acprec%stagger,dimnames,           &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%acprec%rdata_arr,      &
                             dbglvl,istatus)

!--------------------- CUPREC ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%cuprec%fieldname,          &
                             nmmgrid%cuprec%stagger,dimnames,           &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%cuprec%rdata_arr,      &
                             dbglvl,istatus)

    IF ( irain ) THEN
      DO j = nmmgrid%jms, nmmgrid%jme
        DO i = nmmgrid%ims, nmmgrid%ime
          nmmgrid%acprec%rdata_arr(i,j,1) = MAX(0.0,                    &
                          nmmgrid%acprec%rdata_arr(i,j,1) - raing(i,j) )
          nmmgrid%cuprec%rdata_arr(i,j,1) = MAX(0.0,                    &
                          nmmgrid%cuprec%rdata_arr(i,j,1) - rainc(i,j) )
        END DO
      END DO
      DEALLOCATE( raing, rainc )
    END IF

!--------------------- TH2 (TSHLTR) ---------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%tshltr%fieldname,             &
                             nmmgrid%th2%stagger,dimnames,              &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%th2%rdata_arr,         &
                             dbglvl,istatus)

!--------------------- PSHLTR ---------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%pshltr%fieldname,          &
                             nmmgrid%pshltr%stagger,dimnames,           &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%pshltr%rdata_arr,      &
                             dbglvl,istatus)

    nmmgrid%tshltr%rdata_arr(:,:,1) = 0.0
    DO j = nmmgrid%jms, nmmgrid%jme
      DO i = nmmgrid%ims, nmmgrid%ime
        nmmgrid%tshltr%rdata_arr(i,j,1) = nmmgrid%th2%rdata_arr(i,j,1)* &
                           (nmmgrid%pshltr%rdata_arr(i,j,1)/p0)**rddcp
      END DO
    END DO

!--------------------- QSHLTR ---------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%qshltr%fieldname,          &
                             nmmgrid%qshltr%stagger,dimnames,           &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%qshltr%rdata_arr,      &
                             dbglvl,istatus)

!--------------------- U10   ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%u10%fieldname,             &
                             nmmgrid%u10%stagger,dimnames,              &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%u10%rdata_arr,         &
                             dbglvl,istatus)

!--------------------- V10   ------------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%v10%fieldname,             &
                             nmmgrid%v10%stagger,dimnames,              &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%v10%rdata_arr,         &
                             dbglvl,istatus)

!--------------------- SPRING 2011 SPECIFIC FIELDS   -------------------
    IF (nmmgrid%SPRING2011) THEN
      dimnames(3) = ' '
      mEnd(3) = nmmgrid%kps; pEnd(3) = nmmgrid%kps; dEnd(3) = nmmgrid%kps

      CALL nmm_input_get_3dvar(datestr,nmmgrid%u_10m_max%fieldname,     &
                               nmmgrid%u_10m_max%stagger,dimnames,      &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%u_10m_max%rdata_arr, &
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%v_10m_max%fieldname,     &
                               nmmgrid%v_10m_max%stagger,dimnames,      &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%v_10m_max%rdata_arr, &
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%w_up_max%fieldname,      &
                               nmmgrid%w_up_max%stagger,dimnames,       &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%w_up_max%rdata_arr,  &
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%w_dn_max%fieldname,      &
                               nmmgrid%w_dn_max%stagger,dimnames,       &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%w_dn_max%rdata_arr,  &
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%refd_max%fieldname,      &
                               nmmgrid%refd_max%stagger,dimnames,       &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%refd_max%rdata_arr,  &
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%up_heli_max%fieldname,   &
                               nmmgrid%up_heli_max%stagger,dimnames,    &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%up_heli_max%rdata_arr,&
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%rlwdn%fieldname,         &
                               nmmgrid%rlwdn%stagger,dimnames,          &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%rlwdn%rdata_arr,     &
                               dbglvl,istatus)

      CALL nmm_input_get_3dvar(datestr,nmmgrid%rswdn%fieldname,         &
                               nmmgrid%rswdn%stagger,dimnames,          &
                               mStart,mEnd,pStart,pEnd,                 &
                               dStart,dEnd,nmmgrid%rswdn%rdata_arr,     &
                               dbglvl,istatus)
      DO j = nmmgrid%jms, nmmgrid%jme
        DO i = nmmgrid%ims, nmmgrid%ime
          nmmgrid%wspd10max%rdata_arr(i,j,1) = SQRT(                    &
                          nmmgrid%u_10m_max%rdata_arr(i,j,1)**2 +       &
                          nmmgrid%v_10m_max%rdata_arr(i,j,1)**2 )
          nmmgrid%rswdn%rdata_arr(i,j,1) =                              &
                          nmmgrid%rswdn%rdata_arr(i,j,1) +              &
                          nmmgrid%rlwdn%rdata_arr(i,j,1)
          ! rswdn now has total downward radiation flux at surface
        END DO
      END DO

    END IF

!--------------------- SLDPTH ------------------------------------------
    dimnames(1) = dimbtname; dimnames(2) = ' '; dimnames(3) = ' '
    mStart(1) = nmmgrid%kps; pStart(1) = nmmgrid%kps; dStart(1) = nmmgrid%kps
    mStart(2) = nmmgrid%jps; pStart(2) = nmmgrid%jps; dStart(2) = nmmgrid%jps
    mStart(3) = nmmgrid%kps; pStart(3) = nmmgrid%kps; dStart(3) = nmmgrid%kps
    mEnd(1)   = nmmgrid%kpe; pEnd(1)   = nmmgrid%kpe; dEnd(1)   = nmmgrid%kpe
    mEnd(2)   = nmmgrid%jps; pEnd(2)   = nmmgrid%jps; dEnd(2)   = nmmgrid%jps
    mEnd(3)   = nmmgrid%kps; pEnd(3)   = nmmgrid%kps; dEnd(3)   = nmmgrid%kps

    CALL nmm_input_get_3dvar(datestr,nmmgrid%sldpth%fieldname,          &
                             nmmgrid%sldpth%stagger,dimnames,           &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%sldpth%rdata_arr,      &
                             dbglvl,istatus)

    nzsoil = 1
    DO k = nmmgrid%kps+1, nmmgrid%kpe
      IF (nmmgrid%sldpth%rdata_arr(k,1,1) > nmmgrid%sldpth%rdata_arr(k-1,1,1)) THEN
        nzsoil = nzsoil + 1
      ELSE
        EXIT
      END IF
    END DO
    nmmgrid%soil_layer_stag = nzsoil
    nmmgrid%kzs  = 1
    nmmgrid%kzse = nzsoil
    nmmgrid%kze  = nzsoil-1

!--------------------- STC & SMC ------------------------------------------
    dimnames(1) = dimwename; dimnames(2) = dimsnname; dimnames(3) = 'soil_layers_stag'
    mStart(1) = nmmgrid%ims; pStart(1) = nmmgrid%ips; dStart(1) = 1
    mStart(2) = nmmgrid%jms; pStart(2) = nmmgrid%jps; dStart(2) = 1
    mStart(3) = nmmgrid%kzs; pStart(3) = nmmgrid%kzs; dStart(3) = nmmgrid%kzs
    mEnd(1)   = nmmgrid%ime; pEnd(1)   = nmmgrid%ipe; dEnd(1)   = nmmgrid%west_east
    mEnd(2)   = nmmgrid%jme; pEnd(2)   = nmmgrid%jpe; dEnd(2)   = nmmgrid%south_north
    mEnd(3)   = nmmgrid%kzse; pEnd(3)  = nmmgrid%kzse; dEnd(3)  = nmmgrid%kzse

    ALLOCATE(nmmgrid%tsoil%rdata_arr(mStart(1):mEnd(1),mStart(2):mEnd(2),1:nzsoil), STAT = istatus)
    nmmgrid%tsoil%fieldname = 'STC'
    nmmgrid%tsoil%stagger   = 'Z'

    ALLOCATE(nmmgrid%qsoil%rdata_arr(mStart(1):mEnd(1),mStart(2):mEnd(2),1:nzsoil), STAT = istatus)
    nmmgrid%qsoil%fieldname = 'SMC'
    nmmgrid%qsoil%stagger   = 'Z'

    CALL nmm_input_get_3dvar(datestr,nmmgrid%tsoil%fieldname,           &
                             nmmgrid%tsoil%stagger,dimnames,            &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%tsoil%rdata_arr,       &
                             dbglvl,istatus)

    CALL nmm_input_get_3dvar(datestr,nmmgrid%qsoil%fieldname,           &
                             nmmgrid%qsoil%stagger,dimnames,            &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%qsoil%rdata_arr,       &
                             dbglvl,istatus)

!--------------------- CANPWET -----------------------------------------
    dimnames(3) = ' '
    mEnd(3) = nmmgrid%kzs; pEnd(3) = nmmgrid%kzs; dEnd(3) = nmmgrid%kzs

    CALL nmm_input_get_3dvar(datestr,nmmgrid%canpwet%fieldname,         &
                             nmmgrid%canpwet%stagger,dimnames,          &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,nmmgrid%canpwet%rdata_arr,     &
                             dbglvl,istatus)

!--------------------- XTIME ------------------------------------------
    dimnames(1) = ' '; dimnames(2) = ' '; dimnames(3) = ' '
    mStart(1) = 1; pStart(1) = 1; dStart(1) = 1
    mStart(2) = 1; pStart(2) = 1; dStart(2) = 1
    mStart(3) = 1; pStart(3) = 1; dStart(3) = 1
    mEnd(1)   = 1; pEnd(1)   = 1; dEnd(1)   = 1
    mEnd(2)   = 1; pEnd(2)   = 1; dEnd(2)   = 1
    mEnd(3)   = 1; pEnd(3)   = 1; dEnd(3)   = 1

    CALL nmm_input_get_3dvar(datestr,'XTIME',                           &
                             ' ',dimnames,                              &
                             mStart,mEnd,pStart,pEnd,                   &
                             dStart,dEnd,temlocal,                      &
                             dbglvl,istatus)

    nmmgrid%xtime = temlocal(1,1,1)

!--------------------- RETURN   ----------------------------------------

    RETURN
  END SUBROUTINE nmmgrid_getdata

  SUBROUTINE nmmgrid_adjustdata(iamroot, dbglvl, istatus)
!#######################################################################

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: iamroot
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: i,j,k
    INTEGER :: nxs, nys, nzs, nxe, nye, nze

    REAL :: tema, temb, temc
    REAL :: qvmin, qvmax, qvsat

!-----------------------------------------------------------------------

    REAL :: f_qvsat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nxs = nmmgrid%ims
    nxe = nmmgrid%ime
    nys = nmmgrid%jms
    nye = nmmgrid%jme
    nzs = nmmgrid%kps
    nze = nmmgrid%kpe

    ! First rotate winds to earth-relative
    IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Rotating WRF/NMM winds ...'

    DO k = nzs, nze
      CALL map_to_met_nmm(nmmgrid%uu%rdata_arr(:,:,k),1.0,              &
                          nmmgrid%vv%rdata_arr(:,:,k),1.0,              &
                          nmmgrid%ims, nmmgrid%jms,                     &
                          nmmgrid%ime, nmmgrid%jme,                     &
                          nmmgrid%glat%rdata_arr(:,:,1),                &
                          nmmgrid%glon%rdata_arr(:,:,1) )
    END DO

    ! Interpolate pressure to (vertical) scalar points (linear with log(p))
    IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Interpolating scalar pressure ...'

    DO k = nzs, nze
      DO j = nys, nye
        DO i = nxs, nxe
          tema = LOG(nmmgrid%pres%rdata_arr(i,j,k))
          temb = LOG(nmmgrid%pres%rdata_arr(i,j,k+1))
          temc = 0.5*(tema + temb)
          nmmgrid%pp%rdata_arr(i,j,k) = EXP(temc)
          nmmgrid%algp%rdata_arr(i,j,k) = temc
        END DO
      END DO
    END DO

    ! Compute potential temperature based on temperature input
    IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Computing potential temperature ...'

    DO k = nzs, nze
      DO j = nys, nye
        DO i = nxs, nxe
          nmmgrid%pt%rdata_arr(i,j,k) = nmmgrid%tt%rdata_arr(i,j,k)*    &
                               ((p0/nmmgrid%pp%rdata_arr(i,j,k))**rddcp)
        END DO
      END DO
    END DO

!--------------------------- Density -----------------------------------
    ! Compute dentsity
    IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Computing density ...'

    DO k = nzs, nze
      DO j = nys, nye
        DO i = nxs, nxe
          temc = nmmgrid%tt%rdata_arr(i,j,k) /                          &
                 ((1.0+nmmgrid%qv%rdata_arr(i,j,k))*                    &
                  (1.0-(nmmgrid%qv%rdata_arr(i,j,k)/(rddrv+nmmgrid%qv%rdata_arr(i,j,k)))))
          nmmgrid%rho%rdata_arr(i,j,k)=nmmgrid%pp%rdata_arr(i,j,k)/(rd*temc)
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Calculate the height by using hytrostatical equation
!
!-----------------------------------------------------------------------
!
    IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Computing geopotential height ...'

    DO j = nys, nye            ! Terrain height
      DO i= nxs, nxe
        nmmgrid%fis%rdata_arr(i,j,1) = nmmgrid%fis%rdata_arr(i,j,1) / g  ! in meter
      END DO
    END DO

    IF ( ALL(nmmgrid%psfc%rdata_arr < 500.0) ) THEN  ! PSFC not present in data file
      DO j = nys, nye            ! Terrain height
        DO i= nxs, nxe
          !nmmgrid%hgt%rdata_arr(i,j,nzs) = nmmgrid%fis%rdata_arr(i,j,1)
          nmmgrid%zp%rdata_arr(i,j,nzs)  = nmmgrid%fis%rdata_arr(i,j,1)

          tema = nmmgrid%tt%rdata_arr(i,j,nzs)    ! temperature at k-1 is just in the middle of W points k-1 & k
          temb = nmmgrid%qv%rdata_arr(i,j,nzs)
          temc = tema * (1.0 + 0.608*temb)         ! Virtual temperature

          nmmgrid%zp%rdata_arr(i,j,nzs+1) = nmmgrid%zp%rdata_arr(i,j,nzs) + rdovg * temc & ! Height (m)
            * ( log(nmmgrid%pres%rdata_arr(i,j,nzs)) - log(nmmgrid%pres%rdata_arr(i,j,nzs+1)) )

          nmmgrid%hgt%rdata_arr(i,j,nzs) = 0.5* ( nmmgrid%zp%rdata_arr(i,j,nzs) + nmmgrid%zp%rdata_arr(i,j,nzs+1) )
        END DO
      END DO
    ELSE
      ! Compute surface temperature from surface potential temperature
      DO j = nys, nye
        DO i= nxs, nxe
          nmmgrid%ths%rdata_arr(i,j,1) = nmmgrid%ths%rdata_arr(i,j,1)/  &
                             ((p0/nmmgrid%psfc%rdata_arr(i,j,1))**rddcp)
        END DO
      END DO

      DO j = nys, nye
        DO i= nxs, nxe
          tema = 0.5 * ( nmmgrid%ths%rdata_arr(i,j,1) + nmmgrid%tt%rdata_arr(i,j,nzs) )
          temb = 0.5 * ( nmmgrid%qvs%rdata_arr(i,j,1) + nmmgrid%qv%rdata_arr(i,j,nzs) )
          temc = tema * (1.0 + 0.608*temb)         ! Virtual temperature
!          tema = 0.5 * ( nmmgrid%psfc%rdata_arr(i,j,1) + nmmgrid%pp%rdata_arr(i,j,nzs) )
!          temb = temc / tema

          nmmgrid%hgt%rdata_arr(i,j,nzs) = nmmgrid%fis%rdata_arr(i,j,1) &
                                         + rdovg * temc *               &
            ( log(nmmgrid%psfc%rdata_arr(i,j,1)) - nmmgrid%algp%rdata_arr(i,j,nzs) )
                                                              ! Height (m)
          !
          ! W points
          !
!          tema = 0.5 * ( nmmgrid%psfc%rdata_arr(i,j,1) + nmmgrid%pres%rdata_arr(i,j,nzs) )
!          temb = temc / tema

          nmmgrid%zp%rdata_arr(i,j,nzs) = nmmgrid%fis%rdata_arr(i,j,1)
                                                              ! Height (m)
        END DO
      END DO
    END IF

    DO k = nzs+1, nze
      DO j = nys, nye
        DO i = nxs, nxe
          !
          ! W points
          !
          tema = nmmgrid%tt%rdata_arr(i,j,k-1)    ! temperature at k-1 is just in the middle of W points k-1 & k
          temb = nmmgrid%qv%rdata_arr(i,j,k-1)
          temc = tema * (1.0 + 0.608*temb)         ! Virtual temperature

          nmmgrid%zp%rdata_arr(i,j,k) = nmmgrid%zp%rdata_arr(i,j,k-1) + rdovg * temc & ! Height (m)
            * ( log(nmmgrid%pres%rdata_arr(i,j,k-1)) - log(nmmgrid%pres%rdata_arr(i,j,k)) )

          tema = 0.5 * ( nmmgrid%tt%rdata_arr(i,j,k-1) + nmmgrid%tt%rdata_arr(i,j,k) )
          temb = 0.5 * ( nmmgrid%qv%rdata_arr(i,j,k-1) + nmmgrid%qv%rdata_arr(i,j,k) )
          temc = tema * (1.0 + 0.608*temb)         ! Virtual temperature
!          tema = 0.5 * ( nmmgrid%pp%rdata_arr(i,j,k-1) + nmmgrid%pp%rdata_arr(i,j,k) )
!          temb = temc / tema

          nmmgrid%hgt%rdata_arr(i,j,k) = nmmgrid%hgt%rdata_arr(i,j,k-1) + rdovg * temc & ! Height (m)
                   * ( nmmgrid%algp%rdata_arr(i,j,k-1) - nmmgrid%algp%rdata_arr(i,j,k) )



        END DO
      END DO
    END DO
    nmmgrid%zp%rdata_arr(:,:,nze+1)   = 2*nmmgrid%zp%rdata_arr(:,:,nze)  &
                                        - nmmgrid%zp%rdata_arr(:,:,nze-1)
!
!-----------------------------------------------------------------------
!
!  Data transformations - Change qv to RHstar
!    RHStar = sqrt(1. - relative humidity)
!
!-----------------------------------------------------------------------
!
    CALL nmmgrid_qv2rhstar( iamroot, 1, dbglvl, istatus)

!-----------------------------------------------------------------------
!
! Convert precipitation from meter to millimeter
!
!-----------------------------------------------------------------------

    DO j = nys, nye
      DO i = nxs, nxe
        nmmgrid%acprec%rdata_arr(i,j,1) = nmmgrid%acprec%rdata_arr(i,j,1)*1000.
        nmmgrid%cuprec%rdata_arr(i,j,1) = nmmgrid%cuprec%rdata_arr(i,j,1)*1000.
      END DO
    END DO

    RETURN
  END SUBROUTINE nmmgrid_adjustdata

  SUBROUTINE nmmgrid_qv2rhstar( iamroot, direction, dbglvl, istatus )

!#######################################################################

    IMPLICIT NONE
    LOGICAL, INTENT(IN)  :: iamroot
    INTEGER, INTENT(IN)  :: direction
    INTEGER, INTENT(IN)  :: dbglvl

    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    REAL, PARAMETER :: rhmax = 1.0

    INTEGER :: i,j,k
    REAL    :: qvmin, qvmax
    REAL    :: rh, qvsat

    REAL    :: f_qvsat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    qvmax=-999.
    qvmin=999.

    IF (direction > 0) THEN  ! convert qv to rhstar

      IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Data transfering from QV to RHstar ...'

!      IF (dbglvl > 2) THEN
        qvmax = MAXVAL(nmmgrid%qv%rdata_arr)
        qvmin = MINVAL(nmmgrid%qv%rdata_arr)
        CALL mpmax0(qvmax,qvmin)

        IF(IAMROOT) WRITE(6,'(3x,a,2(a,f12.8))')                          &
          'QV     diagnostic information: ',                              &
          'qv_ext min = ',qvmin, ', max = ',qvmax
!      END IF

      DO k = nmmgrid%kps, nmmgrid%kpe
        DO j = nmmgrid%jms, nmmgrid%jme
          DO i = nmmgrid%ims, nmmgrid%ime
            qvsat=f_qvsat( nmmgrid%pp%rdata_arr(i,j,k), nmmgrid%tt%rdata_arr(i,j,k) )
            nmmgrid%qv%rdata_arr(i,j,k)=SQRT(AMAX1(0.,(rhmax-(nmmgrid%qv%rdata_arr(i,j,k)/qvsat))))
          END DO
        END DO
      END DO

!      IF (dbglvl > 2) THEN
        qvmax = MAXVAL(nmmgrid%qv%rdata_arr)
        qvmin = MINVAL(nmmgrid%qv%rdata_arr)
        CALL mpmax0(qvmax,qvmin)

        IF(IAMROOT) WRITE(6,'(3x,a,2(a,f12.8))')                          &
          'RHstar diagnostic information: ',                              &
          'qv_ext min = ',qvmin, ', max = ',qvmax
!      END IF
    ELSE                    ! convert rhstart to qv

      IF (dbglvl > 2) WRITE(6,'(3x,a)') '--- Data transfering from RHstar to QV ...'

      DO k = nmmgrid%kps, nmmgrid%kpe
        DO j = nmmgrid%jms, nmmgrid%jme
          DO i = nmmgrid%ims, nmmgrid%ime
            qvsat=f_qvsat( nmmgrid%pp%rdata_arr(i,j,k), nmmgrid%tt%rdata_arr(i,j,k) )
            rh=AMAX1(0.01,rhmax-(nmmgrid%qv%rdata_arr(i,j,k)*nmmgrid%qv%rdata_arr(i,j,k)))
            nmmgrid%qv%rdata_arr(i,j,k)=rh*qvsat
          END DO
        END DO
      END DO

!      IF ( dbglvl > 2 ) THEN
        qvmax = MAXVAL(nmmgrid%qv%rdata_arr)
        qvmin = MINVAL(nmmgrid%qv%rdata_arr)

        CALL mpmax0(qvmax, qvmin)

        IF ( iamroot ) WRITE(6,'(3x,a,2(a,f12.8))')                     &
          'External QV diagnostic information: ',                       &
          'qv_ext min = ',qvmin, ', max = ',qvmax
!      END IF

    END IF

    RETURN
  END SUBROUTINE nmmgrid_qv2rhstar

END MODULE module_nmm2arps_nmmgrid
