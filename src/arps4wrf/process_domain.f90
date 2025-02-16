!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARPS4WRF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE process_domain(OUNT,domid,IAMROOT,nproc_x,nproc_y,           &
                          nprocx_in,nprocy_in, istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Process a WRF domain for all time levels
!
!-----------------------------------------------------------------------

!  USE module_constants
  USE module_arps4wrf_namelist

  USE module_input_arpsgrid
  USE module_metgrid
  USE module_geogrid

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: OUNT, domid
  LOGICAL, INTENT(IN)  :: IAMROOT
  INTEGER, INTENT(IN)  :: nproc_x, nproc_y
  INTEGER, INTENT(IN)  :: nprocx_in, nprocy_in
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: n
  CHARACTER(LEN=MAXFILELEN) :: out_fname
  CHARACTER(LEN=128)        :: title
  CHARACTER(LEN=19)         :: datestr
  LOGICAL :: arw_matched_grid

!----------------------- External Functions ----------------------------

  LOGICAL :: check_grid_match

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  title   = 'OUTPUT FROM ARPS4WRF'
!-----------------------------------------------------------------------
!
! process static fields
!
!-----------------------------------------------------------------------

  CALL alloc_geogrid(istatus)

  IF (lvldbg > 1) WRITE(ount,'(1x,a)') '--- Reading geo file ... '
  CALL read_static_fields(iamroot,geofile(domid),io_form_input,wrf_core,&
                          lvldbg,istatus)

!-----------------------------------------------------------------------
!
! Loop over all ARPS files
!
!-----------------------------------------------------------------------

  arw_matched_grid = .FALSE.

  DO n = 1, nhisfile(domid)
    IF (lvldbg > 0) WRITE(ount,'(1x,2(a,I2.2),a)')   &
                    '=== Processing #',n,' file for domain d',domid,' ... '
!---------------------------------------------------------------------
! Allocate and initialized ARPS grid
!---------------------------------------------------------------------

    CALL alloc_init_arpsgrid(hisfile(n,domid),hinfmt(domid),iamroot,    &
                             extsoil(domid),soilfile(n,domid),          &
                             nprocx_in,nprocy_in,istatus)

    IF (lvldbg > 1) WRITE(ount,'(1x,a)')             &
                   '--- Reading arps file '//TRIM(hisfile(n,domid))//'... '

    CALL read_arpsgrid(grdbasfn(domid),hisfile(n,domid),hinfmt(domid),  &
                       extsoil(domid),soilfile(n,domid),                &
                       iamroot,datestr,istatus)

    IF (lvldbg > 1) WRITE(OUNT,'(1x,4a)')            &
                   '--- Read in ARPS file - ',TRIM(hisfile(n,domid)),' at ',datestr

    geogrid%bottom_top_dim = arpsgrid%nz-2   ! Vertically the same as ARPS
                                             ! 2 -- nz-1 for W
                                             ! 2 -- nz-2 for P, PT etc.

    IF ( intrp_opt < 100 .AND. n == 1) THEN
      arw_matched_grid = check_grid_match(ount,wrf_core,.TRUE.,lvldbg,istatus)
    END IF

    IF (lvldbg > -1 .AND. n == 1) THEN    ! output corresponding WRF eta values for the ARPS grid
      CALL compute_arpsgrid_eta( iamroot,nproc_x*nproc_y,istatus )
    END IF

!---------------------------------------------------------------------
! Compute loop indices
!---------------------------------------------------------------------

    CALL alloc_metgrid(istatus)

    CALL interpolate_metgrid_from_arpsgrid(ount,wrf_core,lvldbg,        &
                     soilopt(domid),intrp_opt,wtrnopt(domid),           &
                     arw_matched_grid,istatus)

    !IF (wrf_core == 'ARW' .AND. metgrid%qv_no > 0) THEN  ! ARW expects RH, So we make the convertion for it
    IF ( metgrid%qv_no > 0) THEN  ! WRF expects RH, So we make the converting
      CALL metgrid_qv_to_rh( istatus )
    END IF

    IF (lvldbg > 1) WRITE(ount,'(1x,a)') '--- Ouput metgrid ...'

    CALL write_metgrid(iamroot,output_path,io_form_output,domid,        &
                       title,datestr,lvldbg,magnitude_processor,        &
                       out_fname,istatus)

    IF (lvldbg > 1) WRITE(OUNT,'(1x,2a)')            &
                   '--- Data was written into file - ',TRIM(out_fname)

  !-----------------------------------------------------------------------
  ! Clean grids for this time level
  !-----------------------------------------------------------------------

    CALL dealloc_metgrid(istatus)

    CALL dealloc_arpsgrid(istatus)

  END DO
!-----------------------------------------------------------------------
!
! Clean before exit this domain
!
!-----------------------------------------------------------------------

  CALL dealloc_geogrid(istatus)

  RETURN
END SUBROUTINE process_domain

SUBROUTINE interpolate_metgrid_from_arpsgrid(OUNT,wrf_core,dbglvl,      &
                        soilopt,intrp_opt,wtrnopt,arw_matched_grid,istatus)
!
!#######################################################################
!
! Actually do the interpolation for each variable desired
!
!#######################################################################
!
  USE module_input_arpsgrid
  USE module_metgrid
  USE module_geogrid

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: OUNT
  CHARACTER(LEN=3), INTENT(IN) :: wrf_core
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(IN)  :: soilopt, intrp_opt, wtrnopt
  LOGICAL, INTENT(IN)  :: arw_matched_grid
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: bdrwidth
  INTEGER :: do_interp
  LOGICAL ::  extract_only

  CHARACTER (LEN=128) :: fieldname, units, descr
  INTEGER :: iproj

  INTEGER :: mu1s, mu1e, mu2s, mu2e, smu1, smu2
  INTEGER :: mm1s, mm1e, mm2s, mm2e, smm1, smm2
  INTEGER :: mv1s, mv1e, mv2s, mv2e, smv1, smv2

  INTEGER :: pm1s,pm1e,pm2s,pm2e
  INTEGER :: pu1s,pu1e,pu2s,pu2e
  INTEGER :: pv1s,pv1e,pv2s,pv2e
  INTEGER :: pbts

  REAL, ALLOCATABLE :: mx(:,:), my(:,:), ux(:,:), uy(:,:), vx(:,:), vy(:,:)

  INTEGER           :: mxs,mxe, mys, mye
  INTEGER           :: i, j, k, n
  INTEGER           :: doms
  REAL              :: domf, amax
  REAL, ALLOCATABLE :: srcx(:), srcy(:), srcxs(:), srcys(:)
  REAL, ALLOCATABLE :: tem1(:,:,:), tem2(:,:,:)
  REAL, ALLOCATABLE :: tema1(:,:)  ! For interpolating categories
  REAL, ALLOCATABLE :: tem2u(:,:,:), tem2v(:,:,:)
  REAL, ALLOCATABLE :: tem3u(:,:,:), tem3v(:,:,:)
  REAL, ALLOCATABLE :: tem4u(:,:), tem4v(:,:)
  REAL, ALLOCATABLE :: uatv(:,:), vatu(:,:),utem(:,:), vtem(:,:)

  REAL, POINTER :: arpstemqg(:,:,:)
  REAL, ALLOCATABLE :: landsea_external(:,:)

!-----------------------------------------------------------------------

  INTERFACE
    SUBROUTINE interp2d(sdata,mxs,mxe,mys,mye,nx,ny,istagger,sx,sy,dx,dy,   &
                    tx,ty,m1s,m1e,m2s,m2e,interp_method,                    &
                    tdata,ibgn,iend,jbgn,jend,istatus,landsea_mask)

      USE module_wrfgrid_constants

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye,nx,ny
      INTEGER, INTENT(IN)  :: istagger
      REAL,    INTENT(IN)  :: sdata(mxs:mxe,mys:mye)
      REAL,    INTENT(IN)  :: sx(mxs:mxe), sy(mys:mye)
      REAL,    INTENT(IN)  :: dx, dy        ! source grid spacing
      INTEGER, INTENT(IN)  :: m1s,m1e,m2s,m2e
      REAL,    INTENT(IN)  :: tx(m1s:m1e,m2s:m2e), ty(m1s:m1e,m2s:m2e)
      INTEGER, INTENT(IN)  :: interp_method
      INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend

      REAL,    INTENT(OUT) :: tdata(ibgn:iend,jbgn:jend)
      INTEGER, INTENT(OUT) :: istatus
      REAL,    INTENT(IN), OPTIONAL :: landsea_mask(ibgn:iend,jbgn:jend)
    END SUBROUTINE interp2d

  END INTERFACE

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  LOGICAL :: check_grid_match

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

!-----------------------------------------------------------------------
!
! All interpolation to be performed at ARPS grid
!
!-----------------------------------------------------------------------
  IF (dbglvl > 1) WRITE(ount,'(1x,a)') '--- Preparing dimenions for interpolation ...'

  mm1s = geogrid%fields(geogrid%xlat_m_no)%mem_start(1)
  mm1e = geogrid%fields(geogrid%xlat_m_no)%mem_end(1)
  mm2s = geogrid%fields(geogrid%xlat_m_no)%mem_start(2)
  mm2e = geogrid%fields(geogrid%xlat_m_no)%mem_end(2)

  smm1 = mm1e - mm1s + 1  ! memory size 1
  smm2 = mm2e - mm2s + 1  ! memory size 2

  ALLOCATE(mx(mm1s:mm1e,mm2s:mm2e), STAT = istatus)
  ALLOCATE(my(mm1s:mm1e,mm2s:mm2e), STAT = istatus)

  CALL lltoxy(smm1,smm2,geogrid%fields(geogrid%xlat_m_no)%rdata_arr,    &
              geogrid%fields(geogrid%xlong_m_no)%rdata_arr,mx,my)

  mv1s = geogrid%fields(geogrid%xlat_v_no)%mem_start(1)
  mv1e = geogrid%fields(geogrid%xlat_v_no)%mem_end(1)
  mv2s = geogrid%fields(geogrid%xlat_v_no)%mem_start(2)
  mv2e = geogrid%fields(geogrid%xlat_v_no)%mem_end(2)

  smv1 = mv1e - mv1s + 1  ! memory size 1
  smv2 = mv2e - mv2s + 1  ! memory size 2

  ALLOCATE(vx(mv1s:mv1e,mv2s:mv2e), STAT = istatus)
  ALLOCATE(vy(mv1s:mv1e,mv2s:mv2e), STAT = istatus)

  CALL lltoxy(smv1,smv2,geogrid%fields(geogrid%xlat_v_no)%rdata_arr,    &
              geogrid%fields(geogrid%xlong_v_no)%rdata_arr,vx,vy)

  mu1s = geogrid%fields(geogrid%xlat_u_no)%mem_start(1)
  mu1e = geogrid%fields(geogrid%xlat_u_no)%mem_end(1)
  mu2s = geogrid%fields(geogrid%xlat_u_no)%mem_start(2)
  mu2e = geogrid%fields(geogrid%xlat_u_no)%mem_end(2)

  smu1 = mu1e - mu1s + 1  ! memory size 1
  smu2 = mu2e - mu2s + 1  ! memory size 2

  ALLOCATE(ux(mu1s:mu1e,mu2s:mu2e), STAT = istatus)
  ALLOCATE(uy(mu1s:mu1e,mu2s:mu2e), STAT = istatus)

  CALL lltoxy(smu1,smu2,geogrid%fields(geogrid%xlat_u_no)%rdata_arr,  &
              geogrid%fields(geogrid%xlong_u_no)%rdata_arr,ux,uy)

  pm1s = geogrid%we_patch_s
  pm1e = geogrid%we_patch_e
  pm2s = geogrid%sn_patch_s
  pm2e = geogrid%sn_patch_e

  pu1s = geogrid%we_patch_stag_s
  pu1e = geogrid%we_patch_stag_e
  pu2s = geogrid%sn_patch_s
  pu2e = geogrid%sn_patch_e

  pv1s = geogrid%we_patch_s
  pv1e = geogrid%we_patch_e
  pv2s = geogrid%sn_patch_stag_s
  pv2e = geogrid%sn_patch_stag_e

!-----------------------------------------------------------------------
!
! Find halo zone and allocate temporary arrays
!
!-----------------------------------------------------------------------

  extract_only = .FALSE.
  IF (arw_matched_grid) THEN
    bdrwidth = 0
    do_interp = ASSIGN_DIRECT
    IF (dbglvl > 0)  WRITE(ount,'(1x,a)')                               &
          '--- ARW and ARPS grid exact match. So no interpolation involves.'
  ELSE
    CALL find_arpsbdrwidth(mm1s,mm1e,mm2s,mm2e, mu1s,mu1e,mu2s,mu2e,      &
                          mv1s,mv1e,mv2s,mv2e, pm1s,pm1e,pm2s,pm2e,       &
                          pu1s,pu1e,pu2s,pu2e, pv1s,pv1e,pv2s,pv2e,       &
                          arpsgrid%nx,arpsgrid%ny,arpsgrid%dx,arpsgrid%dy,&
                          mx,my,ux,uy,vx,vy,                              &
                          arpsgrid%x,arpsgrid%y,arpsgrid%xs,arpsgrid%ys,  &
                          dbglvl,bdrwidth,istatus)

    IF (intrp_opt < 100)          &
      extract_only = check_grid_match(ount,wrf_core,.FALSE.,dbglvl,istatus)

    IF (extract_only) THEN
      do_interp = EXTRACT_SUBDOMAIN
      IF (dbglvl > 0)  WRITE(ount,'(1x,a,I2,/,5x,a)')                   &
            '--- ARW and ARPS grid match, but ARPS grid is larger with fake zone = ',bdrwidth, &
            'No interpolation involves.'
    ELSE
      do_interp = MOD(intrp_opt,100)
      IF (dbglvl > 0) WRITE(ount,'(1x,a,I2)')                           &
           '--- Message passing fake zone = ',bdrwidth
    END IF

    IF (bdrwidth > arpsgrid%nx/2 .OR. bdrwidth > arpsgrid%ny/2) THEN
      WRITE(6,'(1x,a,3(I4,a))') 'bdrwidth = ',bdrwidth,' is too large, nx = ',arpsgrid%nx,', ny = ',arpsgrid%ny,'.'
      CALL arpsstop('Message passing width is taking over half of the local array size.',1)
    END IF
  END IF

  mxs = 1-bdrwidth
  mxe = arpsgrid%nx+bdrwidth
  mys = 1-bdrwidth
  mye = arpsgrid%ny+bdrwidth

  ALLOCATE(tem1 (mxs:mxe, mys:mye, 1:arpsgrid%nz), STAT = istatus)
  ALLOCATE(srcx (mxs:mxe), STAT = istatus)
  ALLOCATE(srcy (mys:mye), STAT = istatus)
  ALLOCATE(srcxs(mxs:mxe), STAT = istatus)
  ALLOCATE(srcys(mys:mye), STAT = istatus)
  ALLOCATE(landsea_external(mxs:mxe,mys:mye),      STAT = istatus)
  landsea_external(:,:) = 0.0

  CALL set_coord_halo(mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,          &
                      arpsgrid%x,arpsgrid%y,arpsgrid%xs,arpsgrid%ys,    &
                      arpsgrid%dx, arpsgrid%dy,                         &
                      srcx,srcy,srcxs,srcys,istatus)

  pbts = geogrid%bottom_top_dim    ! nz - 2

  ALLOCATE( tem2(pm1s:pm1e,pm2s:pm2e,1:pbts), STAT = istatus )  ! For scalars

!------------------ Pressure -------------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Pressure ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%p(:,:,:)

  CALL exchange_halo(arpsgrid%p,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'PRES'
  units     = ' '
  descr     = ' '
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,pbts,istatus)

!------------------ Height      ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Height ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%zps(:,:,:)

  CALL exchange_halo(arpsgrid%zps,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'GHT'
  units     = 'm'
  descr     = 'Height'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,pbts,istatus)

!------------------ Canwat      ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Canwat ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = arpsgrid%wetcanp(:,:,0)

  CALL exchange_halo(arpsgrid%wetcanp,arpsgrid%nx,arpsgrid%ny,1,        &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'CANWAT'
  units     = 'kg m-2'
  descr     = 'Plant Canopy Surface Water'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)

!------------------ SNOW      ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Snow ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = arpsgrid%snowdpth(:,:)

  CALL exchange_halo(arpsgrid%snowdpth,arpsgrid%nx,arpsgrid%ny,1,       &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'SNOWH'
  units     = 'm'
  descr     = 'Physical Snow Depth'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
  CALL metgrid_add_flag("FLAG_SNOWH",.TRUE.,geogrid%grid_type,istatus)

  ! Convert snow depth in meters to water equiv. of accum. snow depth (kg/m**2)
  ! (where 1 meter liquid water is set equivqlent to 10 meters snow).
  !          0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)
  DO j = mys,mye
    DO i = mxs, mxe
      IF (tem1(i,j,1) /= MISSING_DATA) tem1(i,j,1) = 100*tem1(i,j,1)
    END DO
  END DO

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'SNOW'
  units     = 'kg m-2'
  descr     = 'Water equivalent snow depth'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
  CALL metgrid_add_flag("FLAG_SNOW",.FALSE.,'C',istatus)

!------------------ SOILHGT      ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Soilhgt ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = arpsgrid%zp(:,:,2)

  CALL exchange_halo(arpsgrid%zp(:,:,2),arpsgrid%nx,arpsgrid%ny,1,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'SOILHGT'
  units     = 'm'
  descr     = 'Terrain field of source analysis'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
  CALL metgrid_add_flag("FLAG_SOILHGT",.TRUE.,geogrid%grid_type,istatus)

!------------------ LANDSEA     ---:1:----------------------------------

  ALLOCATE(tema1(1:arpsgrid%nx,1:arpsgrid%ny), STAT = istatus)

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Landsea ...'


  IF ( ANY(arpsgrid%soiltyp(:,:,1) > 0) ) THEN
                      ! ARPS data exists, may be decoded from NAM or GFS
                      ! and only denotes water
    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array from ARPS data ...'

    DO j = 1,arpsgrid%ny      ! Find dominant ARPS soil type
      DO i = 1,arpsgrid%nx
        doms = arpsgrid%soiltyp(i,j,1)
        domf = arpsgrid%stypfrct(i,j,1)
        DO n = 2, arpsgrid%nstyps
          IF (arpsgrid%stypfrct(i,j,n) > domf ) THEN
            doms = arpsgrid%soiltyp(i,j,n)
            domf = arpsgrid%stypfrct(i,j,n)
          END IF
        END DO
        IF (doms == 13) THEN      ! water
          tema1(i,j) = 13.0
        ELSE IF (doms == 12) THEN ! ICE
          tema1(i,j) = 12.0
        ELSE
          tema1(i,j) = 0.0
        END IF
      END DO
    END DO

    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = tema1(:,:)

    CALL exchange_halo(tema1,arpsgrid%nx,arpsgrid%ny,1,    &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                  mx,my,mm1s,mm1e,mm2s,mm2e,NEARNEIGHBOR1,                &
                  tem2,pm1s,pm1e,pm2s,pm2e,istatus)

    !
    ! Remember external land sea flag for later tsoil interpolation
    !
    DO j = 1, arpsgrid%ny
      DO i = 1, arpsgrid%nx
        IF (tema1(i,j) > 0) THEN
          tema1(i,j) = 2.                   ! Water
          landsea_external(i,j) = 2.
        ELSE
          tema1(i,j) = 1.                   ! land
          landsea_external(i,j) = 1.
                                            ! =0, no value
        END IF
      END DO
    END DO

    CALL exchange_halo(tema1,arpsgrid%nx,arpsgrid%ny,1,    &
                       landsea_external,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  ELSE                       ! Get from geogrid
    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array from geogrid data ...'
    DO j = pm2s,pm2e
      DO i = pm1s, pm1e
        tem2(i,j,1) = geogrid%fields(geogrid%LANDMASK_NO)%rdata_arr(i,j,1)
      END DO
    END DO
  END IF

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'LANDSEA'
  units     = 'proprtn'
  descr     = 'Land/Sea flag (1=land,0 or 2=sea)'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)

!------------------ SEAICE      ---:2:----------------------------------
! Reuse temporary array tem1, so keep the order with LANDSEA

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Seaice ...'

  IF ( ALL(arpsgrid%soiltyp(:,:,1) > 0) ) THEN  ! ARPS data exists

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation from ARPS data ...'

    CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                  mx,my,mm1s,mm1e,mm2s,mm2e,NEARNEIGHBOR2,                &
                  tem2,pm1s,pm1e,pm2s,pm2e,istatus,                       &
                  metgrid%fields(metgrid%landsea_no)%rdata_arr)
  ELSE
    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array from geogrid data ...'
    DO j = pm2s,pm2e
      DO i = pm1s, pm1e
        IF( ABS(metgrid%fields(metgrid%LANDSEA_NO)%rdata_arr(i,j,1)) < 1.0E-3 .AND. &
            ABS(geogrid%fields(geogrid%LUINDEX_NO)%rdata_arr(i,j,1)                 &
                                                     -geogrid%isice) < 1.0E-3) THEN
          tem2(i,j,1) = 1.0
        ELSE
          tem2(i,j,1) = 0.0
        END IF
      END DO
    END DO
  END IF

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'SEAICE'
  units     = 'proprtn'
  descr     = 'Ice flag'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)

!------------------ Skintemp    ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Skintemp ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  DO j = 1,arpsgrid%ny               ! may need further consideration
    DO i = 1,arpsgrid%nx
      tema1(i,j)  = arpsgrid%tsoil(i,j,1,0)
      tem1(i,j,1) = arpsgrid%tsoil(i,j,1,0)       ! first soil layer
    END DO
  END DO

  CALL exchange_halo(tema1,arpsgrid%nx,arpsgrid%ny,1,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'SKINTEMP'
  units     = 'K'
  descr     = 'Skin temperature (can use for SST also)'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)

!------------------ PSFC      ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- PSFC ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = arpsgrid%psfc(:,:)

  CALL exchange_halo(arpsgrid%psfc,arpsgrid%nx,arpsgrid%ny,1,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                  &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'PSFC'
  units     = 'Pa'
  descr     = 'Surface Pressure'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
  CALL metgrid_add_flag("FLAG_PSFC",.TRUE.,geogrid%grid_type,istatus)

!------------------ PMSL      ----------------------------------------


  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- PMSL ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  CALL arpsgrid_getpsl(tema1,istatus)

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = tema1(:,:)

  CALL exchange_halo(tema1,arpsgrid%nx,arpsgrid%ny,1,                 &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,         &
                0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                  &
                tem2,pm1s,pm1e,pm2s,pm2e,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'PMSL'
  units     = 'Pa'
  descr     = 'Sea-level Pressure'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',           &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
  CALL metgrid_add_flag("FLAG_SLP",.FALSE.,'C',istatus)

  IF (ALLOCATED(tema1)) DEALLOCATE(tema1, STAT = istatus)    ! masks

!------------------ Temperature ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Potential Temperature ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%pt(:,:,:)

  CALL exchange_halo(arpsgrid%pt,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Converting to temperature ...'
  CALL theta_to_t(tem2,pm1s,pm1e,pm2s,pm2e,1,pbts,                      &
                  metgrid%fields(metgrid%pres_no)%rdata_arr,            &
                  istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'TT'
  units     = 'K'
  descr     = 'Temperature '
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,pbts,istatus)

!------------------ QV ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- QV ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%qv(:,:,:)

  CALL exchange_halo(arpsgrid%qv,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

  DO k = 1,pbts
    DO j = pm2s,pm2e
      DO i = pm1s,pm1e
        tem2(i,j,k) = MAX(0.0,tem2(i,j,k))
      END DO
    END DO
  END DO

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'QV'
  units     = 'kg kg-1'
  descr     = 'Water vapor specific humidity'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                            pm1s,pm1e,pm2s,pm2e,pbts,istatus)
  CALL metgrid_add_flag("FLAG_QV",.TRUE.,geogrid%grid_type,istatus)

!------------------ QC ----------------------------------------

  IF (arpsgrid%P_QC > 0) THEN
    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- QC ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%qscalar(:,:,:,arpsgrid%P_QC)

    CALL exchange_halo(arpsgrid%qscalar(:,:,:,arpsgrid%P_QC),arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                  mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                  tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    DO k = 1,pbts
      DO j = pm2s,pm2e
        DO i = pm1s,pm1e
          tem2(i,j,k) = MAX(0.0,tem2(i,j,k))
        END DO
      END DO
    END DO

    fieldname = 'QC'
    units     = 'kg kg-1'
    descr     = 'Cloud water mixing ratio'
    CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                              pm1s,pm1e,pm2s,pm2e,pbts,istatus)
    CALL metgrid_add_flag("FLAG_QC",.TRUE.,geogrid%grid_type,istatus)
  END IF

!------------------ QR ----------------------------------------

  IF (arpsgrid%P_QR > 0) THEN
    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- QR ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%qscalar(:,:,:,arpsgrid%P_QR)

    CALL exchange_halo(arpsgrid%qscalar(:,:,:,arpsgrid%P_QR),arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                  mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                  tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

    DO k = 1,pbts
      DO j = pm2s,pm2e
        DO i = pm1s,pm1e
          tem2(i,j,k) = MAX(0.0,tem2(i,j,k))
        END DO
      END DO
    END DO

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    fieldname = 'QR'
    units     = 'kg kg-1'
    descr     = 'Rain water mixing ratio'
    CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                              pm1s,pm1e,pm2s,pm2e,pbts,istatus)
    CALL metgrid_add_flag("FLAG_QR",.TRUE.,geogrid%grid_type,istatus)

  END IF

!------------------ QI ----------------------------------------

  IF (arpsgrid%P_QI > 0) THEN
    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- QI ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%qscalar(:,:,:,arpsgrid%P_QI)

    CALL exchange_halo(arpsgrid%qscalar(:,:,:,arpsgrid%P_QI),arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                  mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                  tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

    DO k = 1,pbts
      DO j = pm2s,pm2e
        DO i = pm1s,pm1e
          tem2(i,j,k) = MAX(0.0,tem2(i,j,k))
        END DO
      END DO
    END DO

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    fieldname = 'QI'
    units     = 'kg kg-1'
    descr     = 'Ice water mixing ratio'
    CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                              pm1s,pm1e,pm2s,pm2e,pbts,istatus)
    CALL metgrid_add_flag("FLAG_QI",.TRUE.,geogrid%grid_type,istatus)
  END IF

!------------------ QS ----------------------------------------

  IF (arpsgrid%P_QS > 0) THEN
    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- QS ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%qscalar(:,:,:,arpsgrid%P_QS)

    CALL exchange_halo(arpsgrid%qscalar(:,:,:,arpsgrid%P_QS),arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                  mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                  tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

    DO k = 1,pbts
      DO j = pm2s,pm2e
        DO i = pm1s,pm1e
          tem2(i,j,k) = MAX(0.0,tem2(i,j,k))
        END DO
      END DO
    END DO

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    fieldname = 'QS'
    units     = 'kg kg-1'
    descr     = 'Snow water mixing ratio'
    CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                              pm1s,pm1e,pm2s,pm2e,pbts,istatus)
    CALL metgrid_add_flag("FLAG_QS",.TRUE.,geogrid%grid_type,istatus)
  END IF

!------------------ QH ----------------------------------------

  IF (arpsgrid%P_QG > 0 .OR. arpsgrid%P_QH > 0) THEN

    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- QH ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

    tem1 = 0.0
    IF ( arpsgrid%P_QG > 0 .AND. arpsgrid%P_QH > 0 ) THEN
      ALLOCATE(arpstemqg(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)
      arpstemqg(:,:,:) = arpsgrid%qscalar(:,:,:,arpsgrid%P_QG) + arpsgrid%qscalar(:,:,:,arpsgrid%P_QH)
    ELSE IF ( arpsgrid%P_QG > 0) THEN
      arpstemqg => arpsgrid%qscalar(:,:,:,arpsgrid%P_QG)
    ELSE IF ( arpsgrid%P_QH > 0) THEN
      arpstemqg => arpsgrid%qscalar(:,:,:,arpsgrid%P_QH)
    END IF
    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpstemqg

    CALL exchange_halo(arpstemqg,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,     &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF ( arpsgrid%P_QG > 0 .AND. arpsgrid%P_QH > 0 ) THEN
      DEALLOCATE( arpstemqg )
    END IF

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                  arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                  mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                  tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)

    DO k = 1,pbts
      DO j = pm2s,pm2e
        DO i = pm1s,pm1e
          tem2(i,j,k) = MAX(0.0,tem2(i,j,k))
        END DO
      END DO
    END DO

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    fieldname = 'QG'
    units     = 'kg kg-1'
    descr     = 'Hail water mixing ratio'
    CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'M',             &
                              pm1s,pm1e,pm2s,pm2e,pbts,istatus)
    CALL metgrid_add_flag("FLAG_QG",.TRUE.,geogrid%grid_type,istatus)
  END IF

!------------------ Wind U ---------------------------------------------

  ALLOCATE( tem2u(pu1s:pu1e,pu2s:pu2e,1:pbts), STAT = istatus )  ! For U stagger

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- U ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%u(:,:,:)

  CALL exchange_halo(arpsgrid%u,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                     tem1,mxs,mxe,mys,mye,1,istatus)   ! U stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,1,srcx,srcys,arpsgrid%dx,arpsgrid%dy,       &
                ux,uy,mu1s,mu1e,mu2s,mu2e,do_interp,                    &
                tem2u,pu1s,pu1e,pu2s,pu2e,pbts,istatus)

!------------------ Wind V ---------------------------------------------

  ALLOCATE( tem2v(pv1s:pv1e,pv2s:pv2e,1:pbts), STAT = istatus )  ! For V stagger

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- V ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%v(:,:,:)

  CALL exchange_halo(arpsgrid%v,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,    &
                     tem1,mxs,mxe,mys,mye,2,istatus)   ! V stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,2,srcxs,srcy,arpsgrid%dx,arpsgrid%dy,       &
                vx,vy,mv1s,mv1e,mv2s,mv2e,do_interp,                    &
                tem2v,pv1s,pv1e,pv2s,pv2e,pbts,istatus)

!-----------------------------------------------------------------------
!
! Rotate wind to earth-relative first, and then
! Rotate wind to WRF grid relative
!
!-----------------------------------------------------------------------

  IF (.NOT. arw_matched_grid .AND. .NOT. extract_only) THEN
    ALLOCATE(tem3u(mu1s:mu1e,mu2s:mu2e,1:pbts), STAT = istatus)
    ALLOCATE(tem3v(mv1s:mv1e,mv2s:mv2e,1:pbts), STAT = istatus)

    tem3u(pu1s:pu1e,pu2s:pu2e,1:pbts) = tem2u(pu1s:pu1e,pu2s:pu2e,1:pbts)
    CALL exchange_halo_r(tem3u,mu1s,mu1e,mu2s,mu2e,1,pbts,pu1s,pu1e,pu2s,pu2e,1,pbts)

    tem3v(pv1s:pv1e,pv2s:pv2e,1:pbts) = tem2v(pv1s:pv1e,pv2s:pv2e,1:pbts)
    CALL exchange_halo_r(tem3v,mv1s,mv1e,mv2s,mv2e,1,pbts,pv1s,pv1e,pv2s,pv2e,1,pbts)

    ALLOCATE(tem4u(mu1s:mu1e,mu2s:mu2e), STAT = istatus)
    ALLOCATE(tem4v(mv1s:mv1e,mv2s:mv2e), STAT = istatus)

    IF (geogrid%grid_type == 'C') THEN

      IF (dbglvl > 0) WRITE(6,'(1x,a)') 'Rotating earth winds to C grid winds'

      ALLOCATE(vatu(mu1s:mu1e,mu2s:mu2e), STAT = istatus )
      ALLOCATE(uatv(mv1s:mv1e,mv2s:mv2e), STAT = istatus )

      ALLOCATE(utem(mu1s:mu1e,mu2s:mu2e), STAT = istatus )
      ALLOCATE(vtem(mv1s:mv1e,mv2s:mv2e), STAT = istatus )

      DO k=1,pbts

        ! get u at v grid point locations
        DO j=mu2s+1,mu2e
          DO i=mu1s,mu1e-1
            uatv(i,j) = 0.25*(tem3u(i,j-1,k)+tem3u(i+1,j-1,k) +  &
                              tem3u(i,j,k)  +tem3u(i+1,j,k))
          END DO
        END DO
        DO i=mu1s,mu1e-1
          uatv(i,mu2s) = 0.5*(tem3u(i,mu2s,k)+tem3u(i+1,mu2s,k))
        END DO

        IF (mu1e == mv1e) THEN
          DO j=mu2s+1,mu2e
            uatv(mu1e,j) = 0.5*(tem3u(mu1e,j-1,k)+tem3u(mu1e,j,k))
          END DO
          uatv(mu1e,mu2s) = tem3u(mu1e,mu2s,k)
        END IF
        IF (mv2e > mu2e) THEN
          DO i = mv1s,mv1e
            uatv(i,mv2e) = uatv(i,mv2e-1)
          END DO
        END IF

        ! get v at u grid point locations
        DO j=mv2s,mv2e-1
          DO i=mv1s+1,mv1e
            vatu(i,j) = 0.25*(tem3v(i-1,j,k)+tem3v(i,j,k) +  &
                              tem3v(i-1,j+1,k)+tem3v(i,j+1,k))
          END DO
        END DO
        DO j=mv2s,mv2e-1      ! filling for missing index due to "mv1s+1" above
          vatu(mv1s,j) = 0.5*(tem3v(mv1s,j,k)+tem3v(mv1s,j+1,k))
        END DO

        IF (mv2e == mu2e) THEN     ! note mv2e-1 above
          DO i=mv1s+1,mv1e
            vatu(i,mv2e) = 0.5*(tem3v(i-1,mv2e,k)+tem3v(i,mv2e,k))
          END DO
          vatu(mv1s,mv2e) = tem3v(mv1s,mv2e,k)
        END IF
        IF (mu1e > mv1e) THEN     ! consider U stagger
          DO j = mu2s, mu2e
            vatu(mu1e,j) = vatu(mu1e-1,j)
          END DO
        END IF

        ! ARPS subroutine
        CALL uvmptoe(smu1,smu2,tem3u(:,:,k),vatu,                         &
                     geogrid%fields(geogrid%xlong_u_no)%rdata_arr,        &
                     tem4u,utem)

        CALL uvmptoe(smv1,smv2,uatv,tem3v(:,:,k),                         &
                     geogrid%fields(geogrid%xlong_v_no)%rdata_arr,        &
                     vtem,tem4v)

        ! WPS subroutine
        CALL met_to_map(tem4u,1.0,tem4v,1.0,                              &
                        mu1s,mu2s,mu1e,mu2e, mv1s,mv2s,mv1e,mv2e,         &
                        geogrid%fields(geogrid%xlong_u_no)%rdata_arr,     &
                        geogrid%fields(geogrid%xlong_v_no)%rdata_arr,     &
                        geogrid%fields(geogrid%xlat_u_no)%rdata_arr,      &
                        geogrid%fields(geogrid%xlat_v_no)%rdata_arr )

        DO j = pu2s,pu2e
          DO i = pu1s,pu1e
            tem2u(i,j,k) = tem4u(i,j)
          END DO
        END DO

        DO j = pv2s,pv2e
          DO i = pv1s,pv1e
            tem2v(i,j,k) = tem4v(i,j)
          END DO
        END DO

      END DO

      DEALLOCATE(uatv, vatu)
      DEALLOCATE(vtem, utem)

    ELSE

      IF (dbglvl > 0) WRITE(6,'(1x,a)') 'Rotating earth winds to E grid winds'

      DO k = 1,pbts
        CALL uvmptoe(smv1,smv2,tem3u(:,:,k),tem3v(:,:,k),            &
                     geogrid%fields(geogrid%xlong_v_no)%rdata_arr,   &
                     tem4u(:,:),tem4v(:,:))    ! ARPS sub
        CALL met_to_map_nmm(tem4u,1.0,tem4v,1.0,mu1s,mv2s,mu1e,mv2e, &
                     geogrid%fields(geogrid%xlat_v_no)%rdata_arr,    &
                     geogrid%fields(geogrid%xlong_v_no)%rdata_arr )
                                               ! WPS sub
        DO j = pv2s,pv2e
          DO i = pu1s,pu1e
            tem2u(i,j,k) = tem4u(i,j)
            tem2v(i,j,k) = tem4v(i,j)
          END DO
        END DO

      END DO

    END IF

    DEALLOCATE(tem4u, tem4v, STAT = istatus)
    DEALLOCATE(tem3u, tem3v, STAT = istatus)
  END IF

!-----------------------------------------------------------------------
!
! Added U/V winds to the metgrid
!
!-----------------------------------------------------------------------

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output UU array ...'

  fieldname = 'UU'
  units     = 'm s-1'
  descr     = 'U'
  CALL metgrid_add_3d_field(tem2u,fieldname,units,descr,'U',            &
                            pu1s,pu1e,pu2s,pu2e,pbts,istatus)

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output VV array ...'

  fieldname = 'VV'
  units     = 'm s-1'
  descr     = 'V'
  CALL metgrid_add_3d_field(tem2v,fieldname,units,descr,'V',            &
                            pv1s,pv1e,pv2s,pv2e,pbts,istatus)

!------------------ Wind - W ----------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Wind W ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  tem1(1:arpsgrid%nx,1:arpsgrid%ny,1:arpsgrid%nz) = arpsgrid%w(:,:,:)

  CALL exchange_halo(arpsgrid%zp,arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,   &
                     tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

  CALL interp3d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,           &
                arpsgrid%nz,0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,      &
                mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
                tem2,pm1s,pm1e,pm2s,pm2e,pbts,istatus)  ! pbts = nz-2

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

  fieldname = 'WW'
  units     = 'm s-1'
  descr     = 'W'
  CALL metgrid_add_3d_field(tem2,fieldname,units,descr,'W',             &
                            pm1s,pm1e,pm2s,pm2e,pbts,istatus)

  CALL metgrid_add_flag("USE_INPUT_W",.TRUE.,geogrid%grid_type,istatus)

!-----------------------------------------------------------------------
!
! Process land soil model variable (soil temperature, soil moisture and
! SST etc.)
!
!-----------------------------------------------------------------------

  CALL interpolate_soil_from_arps( wrf_core, soilopt, landsea_external, &
          geogrid%fields(geogrid%LANDMASK_NO)%rdata_arr(:,:,1),         &
          mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,arpsgrid%nzsoil,      &
          srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                          &
          mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                          &
          pm1s,pm1e,pm2s,pm2e,                                          &
          tem1,tem2,dbglvl,istatus )

!-----------------------------------------------------------------------
!
! At last, adjust terrain height if required to do so
!
!-----------------------------------------------------------------------

!------------------ HGT_M ------------------------------------------

  IF (wtrnopt == 0) THEN   ! Use ARPS Terrain

    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- HGT_M ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

    tem1(1:arpsgrid%nx,1:arpsgrid%ny,1) = arpsgrid%zp(:,:,2)

    CALL exchange_halo(arpsgrid%zp(:,:,2),arpsgrid%nx,arpsgrid%ny,1,    &
                       tem1,mxs,mxe,mys,mye,0,istatus)   ! no stagger

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,         &
                  0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                &
                  mx,my,mm1s,mm1e,mm2s,mm2e,do_interp,                  &
                  tem2,pm1s,pm1e,pm2s,pm2e,istatus)

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    fieldname = 'HGT_M'
    units     = 'meters MSL'
    descr     = 'Topography height'
    CALL metgrid_replace_terrain_field(tem2,fieldname,units,descr,'M',  &
                              pm1s,pm1e,pm2s,pm2e,istatus)

!------------------ HGT_V ------------------------------------------

    IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- HGT_V ...'

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,         &
                  0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,                &
                  vx,vy,mv1s,mv1e,mv2s,mv2e,FOUR_PNT ,                  &
                  tem2v,pv1s,pv1e,pv2s,pv2e,istatus)

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

    fieldname = 'HGT_V'
    units     = 'meters MSL'
    descr     = 'Topography height'
    CALL metgrid_replace_terrain_field(tem2v,fieldname,units,descr,'V',  &
                              pv1s,pv1e,pv2s,pv2e,istatus)

    IF (wrf_core == 'ARW') THEN

!------------------ HGT_U ------------------------------------------

      IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- HGT_U ...'

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp2d(tem1,mxs,mxe,mys,mye,arpsgrid%nx,arpsgrid%ny,       &
                    0,srcxs,srcys,arpsgrid%dx,arpsgrid%dy,              &
                    ux,uy,mu1s,mu1e,mu2s,mu2e,FOUR_PNT,                 &
                    tem2u,pu1s,pu1e,pu2s,pu2e,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill output grid ...'

      fieldname = 'HGT_U'
      units     = 'meters MSL'
      descr     = 'Topography height'
      CALL metgrid_replace_terrain_field(tem2u,fieldname,units,descr,'U', &
                                pu1s,pu1e,pu2s,pu2e,istatus)

    END IF

  END IF

!-----------------------------------------------------------------------
!
! Ending of the subroutine
!
!-----------------------------------------------------------------------

  DEALLOCATE(tem2, STAT = istatus)    ! Deallocate scalars
  DEALLOCATE(tem2u, tem2v)

  DEALLOCATE( mx,my, ux,uy, vx,vy )
  DEALLOCATE( srcx, srcxs, srcy, srcys, tem1 )
  DEALLOCATE( landsea_external )

  RETURN
END SUBROUTINE interpolate_metgrid_from_arpsgrid

FUNCTION check_grid_match(OUNT,wrf_core,exact,dbglvl,istatus)
!
!#######################################################################
!
! Check whether it is ARW core and the ARPS grid matches that of the
! WRF grid
!
!#######################################################################
!
  USE module_input_arpsgrid
  USE module_metgrid
  USE module_geogrid

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: OUNT
  LOGICAL, INTENT(IN)  :: exact
  INTEGER, INTENT(IN)  :: dbglvl
  CHARACTER(LEN=3), INTENT(IN)  :: wrf_core
  INTEGER, INTENT(OUT) :: istatus

  LOGICAL :: check_grid_match

!-----------------------------------------------------------------------

  REAL, PARAMETER :: smallreal = 1.0E-4

  INTEGER :: match_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  check_grid_match = .FALSE.

  IF (wrf_core == 'ARW') THEN
    match_grid = 0

    IF (dbglvl > 10 .AND. exact) THEN
      WRITE(6,'(/1x,a/)')  'The ARPS grid and the WRF grid are:'
      WRITE(6,'(a38)')  '              ARPS grid      WRF grid'
      WRITE(6,'(a38)')  '              ==========    =========='
      WRITE(6,'(a,I10,4x,I10)')     '  nx       = ',arpsgrid%nxlg,geogrid%west_east_dim
      WRITE(6,'(a,I10,4x,I10)')     '  ny       = ',arpsgrid%nylg,geogrid%south_north_dim
      WRITE(6,'(a,I10,4x,I10)')     '  mapproj  = ',arpsgrid%mapproj,geogrid%map_proj
      WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat1  = ',arpsgrid%trulat1,geogrid%truelat1
      WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat2  = ',arpsgrid%trulat2,geogrid%truelat2
      WRITE(6,'(a,F10.2,4x,F10.2)') '  trulon   = ',arpsgrid%trulon,geogrid%stand_lon
      WRITE(6,'(a,F10.2,4x,F10.2)') '  ctrlat   = ',arpsgrid%ctrlat,geogrid%cen_lat
      WRITE(6,'(a,F10.2,4x,F10.2)') '  ctrlon   = ',arpsgrid%ctrlon,geogrid%cen_lon
      WRITE(6,'(a,F10.0,4x,F10.0)') '  dx       = ',arpsgrid%dx,geogrid%dx
      WRITE(6,'(a,F10.0,4x,F10.0)') '  dy       = ',arpsgrid%dy,geogrid%dy
    END IF

    IF (ABS(arpsgrid%dx - geogrid%dx) < smallreal    .AND.              &
        ABS(arpsgrid%dy - geogrid%dy) < smallreal    .AND.              &
        ABS(arpsgrid%trulat1 - geogrid%truelat1) < smallreal .AND.      &
        ABS(arpsgrid%trulat2 - geogrid%truelat2) < smallreal .AND.      &
        ABS(arpsgrid%trulon  - geogrid%stand_lon) < smallreal .AND.     &
        ABS(arpsgrid%ctrlat  - geogrid%cen_lat) < smallreal .AND.       &
        ABS(arpsgrid%ctrlon  - geogrid%cen_lon) < smallreal ) THEN

      IF ( (arpsgrid%mapproj == 0 .AND. geogrid%map_proj == 0) .OR.     &
           (arpsgrid%mapproj == 1 .AND. geogrid%map_proj == 2) .OR.     &
           (arpsgrid%mapproj == 2 .AND. geogrid%map_proj == 1) .OR.     &
           (arpsgrid%mapproj == 3 .AND. geogrid%map_proj == 3) ) THEN

        IF (exact) THEN   ! Require same size
          IF(arpsgrid%nxlg == geogrid%west_east_dim+2 .AND.                   &
             arpsgrid%nylg == geogrid%south_north_dim+2 .AND.                 &
             arpsgrid%nx-3 == geogrid%we_patch_e-geogrid%we_patch_s+1 .AND.   &
             arpsgrid%ny-3 == geogrid%sn_patch_e-geogrid%sn_patch_s+1 ) THEN
            match_grid = 1
          END IF
        ELSE
          match_grid = 1
        END IF
      END IF
    END IF

    CALL mpmini(match_grid)
    IF (match_grid > 0) THEN
      check_grid_match = .TRUE.
      IF (dbglvl > 0 .AND. exact) THEN
        WRITE(OUNT,'(/,1x,a,/,9x,a,/)')    &
          '*NOTE*: ARW core detected, furthermore, it is found that the horizontal grid', &
          'of ARPS and WRF are the same. So grid_match flag is set to avoid interpolation.'
      END IF
    END IF
  END IF

  RETURN
END FUNCTION check_grid_match

!#######################################################################

SUBROUTINE interpolate_soil_from_arps( wrf_core,soilopt, landsea_arps,  &
          landmask_geogrid,                                             &
          mxs,mxe,mys,mye, nx,ny,nzsoil,                                &
          srcxs,srcys,dx,dy,                                            &
          mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,                    &
          pm1s,pm1e,pm2s,pm2e,                                          &
          tem1,tem2,dbglvl,istatus )

!-----------------------------------------------------------------------
!
! PURPOSE:
!   Interpolate soil variables (soil temperature & soil moisture) from
!   ARPS grid to the METO grid.
!
! NOTE:
!   It was modied here to add mask capability for better handling of
!   these variables.
!
!-----------------------------------------------------------------------
! AUTHOR: Yunheng Wang (10/25/2010)
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
  USE module_input_arpsgrid
  USE module_metgrid
  !USE module_geogrid

  IMPLICIT NONE

  CHARACTER(LEN=3), INTENT(IN) :: wrf_core
  INTEGER, INTENT(IN)  :: soilopt
  INTEGER, INTENT(IN)  :: nx, ny, nzsoil
  INTEGER, INTENT(IN)  :: mxs, mxe, mys, mye
  INTEGER, INTENT(IN)  :: landsea_arps(mxs:mxe,mys:mye)
  REAL,    INTENT(IN)  :: dx, dy
  REAL,    INTENT(IN)  :: srcxs(mxs:mxe), srcys(mys:mye)

  INTEGER, INTENT(IN)  :: mm1s, mm1e, mm2s, mm2e
  INTEGER, INTENT(IN)  :: pm1s, pm1e, pm2s, pm2e
  REAL,    INTENT(IN)  :: mdesx(mm1s:mm1e,mm2s:mm2e),   & ! Meto grid loation on ARPS grid
                          mdesy(mm1s:mm1e,mm2s:mm2e)
  REAL,    INTENT(IN)  :: landmask_geogrid(pm1s:pm1e,pm2s:pm2e)

  INTEGER, INTENT(IN)  :: do_interp           ! interpolation option

  REAL,    INTENT(INOUT) :: tem1(mxs:mxe,mys:mye,1:nzsoil)
  REAL,    INTENT(INOUT) :: tem2(pm1s:pm1e,pm2s:pm2e,1:nzsoil)

  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER, PARAMETER  :: ount = 6

  CHARACTER (LEN=128) :: fieldname, units, descr

  REAL, ALLOCATABLE :: sst(:,:), soil_layers(:,:,:)
  REAL, ALLOCATABLE :: st_arps(:,:,:), sst_arps(:,:)
  INTEGER :: i,j,k,kk

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ALLOCATE(st_arps (1:nx,1:ny,1:nzsoil-1),   STAT = istatus)
  !ALLOCATE(sst_arps(1:nx,1:ny),              STAT = istatus)

  !ALLOCATE(sst     (pm1s:pm1e,pm2s:pm2e),              STAT = istatus)
  IF (wrf_core == 'ARW') THEN
    IF (soilopt == 3) THEN
      ALLOCATE(soil_layers(pm1s:pm1e,pm2s:pm2e,1:nzsoil-1), STAT = istatus)
    ELSE
      ALLOCATE(soil_layers(pm1s:pm1e,pm2s:pm2e,1:nzsoil), STAT = istatus)
    END IF
  END IF

!-----------------------------------------------------------------------
! SOIL Temperature
!-----------------------------------------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Soil temperature ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  IF (soilopt == 3) THEN ! Keep NAM or GFS soil levels

    SELECT CASE (arpsgrid%nzsoil)
    CASE (6)              ! NAM 4 layer model

      !sst_arps(1:nx,1:ny) = arpsgrid%tsoil(:,:,1,0)  ! first level is surface

      DO k = 2,5
        DO j = 1,ny
          DO i = 1, nx
            kk = nzsoil-k
            st_arps(i,j,kk) = arpsgrid%tsoil(i,j,k,0)
            tem1(i,j,kk)    = arpsgrid%tsoil(i,j,k,0)
          END DO
        END DO
      END DO

      CALL exchange_halo(st_arps,nx,ny,nzsoil-2,            &
                         tem1,mxs,mxe,mys,mye,0,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil-2,0,         &
                    srcxs,srcys,dx,dy,                                  &
                    mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,          &
                    tem2,pm1s,pm1e,pm2s,pm2e,nzsoil-2,                  &
                    landsea_arps,landmask_geogrid,285.0,istatus)

      fieldname = 'ST'
      units     = ' '
      descr     = ' '
      CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',     &
                            pm1s,pm1e,pm2s,pm2e,nzsoil-2,istatus)

      IF (wrf_core == 'ARW') THEN
        !
        ! Added soil_layers
        !
        soil_layers(:,:,1) = 200.
        soil_layers(:,:,2) = 100.
        soil_layers(:,:,3) =  40.
        soil_layers(:,:,4) =  10.

        fieldname = 'SOIL_LAYERS'
        units     = ''
        descr     = ''
        CALL metgrid_add_3dsoil_field(soil_layers,fieldname,units,descr,'M', &
                              pm1s,pm1e,pm2s,pm2e,4,istatus)
        CALL metgrid_add_flag("FLAG_SOIL_LAYERS",.TRUE.,geogrid%grid_type,istatus)
      END IF
      metgrid%num_soil_levels = 4

      !
      ! Added separated soil layers
      !
      IF (dbglvl > 2) WRITE(ount,'(3x,a)')     &
        '~~~ Fill output grid from NAM 4 soil layer directly...'

      fieldname = 'ST000010'
      units     = 'K'
      descr     = 'T 0-10 cm below ground layer (Upper) '
      CALL metgrid_add_3dsoil_field(tem2(:,:,4),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_ST000010",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'ST010040'
      units     = 'K'
      descr     = 'T 10-40 cm below ground layer (Upper) '
      CALL metgrid_add_3dsoil_field(tem2(:,:,3),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_ST010040",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'ST040100'
      units     = 'K'
      descr     = 'T 40-100 cm below ground layer (Upper) '
      CALL metgrid_add_3dsoil_field(tem2(:,:,2),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_ST040100",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'ST100200'
      units     = 'K'
      descr     = 'T 100-200 cm below ground layer (Bottom)'
      CALL metgrid_add_3dsoil_field(tem2(:,:,1),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_ST100200",.TRUE.,geogrid%grid_type,istatus)

    CASE (3)              ! GFS 2 layer model

      !sst_arps(1:nx,1:ny) = arpsgrid%tsoil(:,:,1,0)  ! first level is surface

      DO k = 2,3
        DO j = 1,ny
          DO i = 1, nx
            kk = nzsoil-k+1
            st_arps(i,j,kk) = arpsgrid%tsoil(i,j,k,0)
            tem1(i,j,kk)    = arpsgrid%tsoil(i,j,k,0)
          END DO
        END DO
      END DO

      CALL exchange_halo(st_arps,nx,ny,nzsoil-1,                        &
                         tem1,mxs,mxe,mys,mye,0,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil-1,0,         &
                    srcxs,srcys,dx,dy,                                  &
                    mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,          &
                    tem2,pm1s,pm1e,pm2s,pm2e,nzsoil-1,                  &
                    landsea_arps,landmask_geogrid,285.0,istatus)

      fieldname = 'ST'
      units     = ' '
      descr     = ' '
      CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',     &
                            pm1s,pm1e,pm2s,pm2e,nzsoil-1,istatus)

      IF (wrf_core == 'ARW') THEN
        !
        ! Added soil_layers
        !
        soil_layers(:,:,1) = 200.
        soil_layers(:,:,2) =  10.

        fieldname = 'SOIL_LAYERS'
        units     = ''
        descr     = ''
        CALL metgrid_add_3dsoil_field(soil_layers,fieldname,units,descr,'M', &
                              pm1s,pm1e,pm2s,pm2e,2,istatus)
        CALL metgrid_add_flag("FLAG_SOIL_LAYERS",.TRUE.,geogrid%grid_type,istatus)
      END IF
      metgrid%num_soil_levels = 2
      !
      ! Added separated soil layers
      !

      IF (dbglvl > 2) WRITE(ount,'(3x,a)')     &
        '~~~ Fill output grid from GFS 2 soil layer directly...'


      fieldname = 'ST000010'
      units     = 'K'
      descr     = 'T 0-10 cm below ground layer (Upper)'
      CALL metgrid_add_3dsoil_field(tem2(:,:,2),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_ST000010",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'ST010200'
      units     = 'K'
      descr     = 'T 010-200 cm below ground layer (Bottom)'
      CALL metgrid_add_3dsoil_field(tem2(:,:,1),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_ST010200",.TRUE.,geogrid%grid_type,istatus)

    CASE DEFAULT
      CALL arpsstop('ERROR: Unknown number of soil layers.',1)
    END SELECT

!-----------------------------------------------------------------------
!
! ARPS traditional 2-layer model
!
!-----------------------------------------------------------------------

  ELSE IF (soilopt == 1) THEN

    tem1(1:nx,1:ny,1:nzsoil) = arpsgrid%tsoil(:,:,:,0)

    CALL exchange_halo(arpsgrid%tsoil(:,:,:,0),nx,ny,nzsoil,            &
                       tem1,mxs,mxe,mys,mye,0,istatus)

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil,0,             &
                  srcxs,srcys,dx,dy,                                    &
                  mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,            &
                  tem2,pm1s,pm1e,pm2s,pm2e,nzsoil,                      &
                  landsea_arps,landmask_geogrid,285.0,istatus)

    fieldname = 'ST'
    units     = ' '
    descr     = ' '
    CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',       &
                          pm1s,pm1e,pm2s,pm2e,arpsgrid%nzsoil,istatus)

    IF (wrf_core == 'ARW') THEN
      !
      ! Added soil_layers
      !
      soil_layers(:,:,1) =  10.
      soil_layers(:,:,2) = 200.

      fieldname = 'SOIL_LAYERS'
      units     = ''
      descr     = ''
      CALL metgrid_add_3dsoil_field(soil_layers,fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,2,istatus)
      CALL metgrid_add_flag("FLAG_SOIL_LAYERS",.TRUE.,geogrid%grid_type,istatus)
    END IF
    metgrid%num_soil_levels = 2

    !
    ! Fill separate layers
    !
    IF (dbglvl > 2) WRITE(ount,'(3x,a)')     &
      '~~~ Fill output grid from ARPS 2 soil layer directly...'

    fieldname = 'ST000010'
    units     = 'K'
    descr     = 'T 0-10 cm below ground layer (Upper) '
    CALL metgrid_add_3dsoil_field(tem2(:,:,1),fieldname,units,descr,'M', &
                          pm1s,pm1e,pm2s,pm2e,1,istatus)
    CALL metgrid_add_flag("FLAG_ST000010",.TRUE.,geogrid%grid_type,istatus)

    fieldname = 'ST010200'
    units     = 'K'
    descr     = 'T 10-200 cm below ground layer (Bottom) '
    CALL metgrid_add_3dsoil_field(tem2(:,:,2),fieldname,units,descr,'M', &
                          pm1s,pm1e,pm2s,pm2e,1,istatus)
    CALL metgrid_add_flag("FLAG_ST010200",.TRUE.,geogrid%grid_type,istatus)

!-----------------------------------------------------------------------
!
! ECMWF soil layers 4 layers at
!   1. 0 - 7   cm
!   2. 7 - 28  cm
!   3. 28 - 100 cm
!   4. 100 - 255 cm
!
!-----------------------------------------------------------------------
  ELSE IF (soilopt == 4) THEN

    IF( arpsgrid%nzsoil == 6) THEN   ! ECMWF 4 layer model

      !sst_arps(1:nx,1:ny) = arpsgrid%tsoil(:,:,1,0)  ! first level is surface

      DO k = 2,nzsoil
        DO j = 1,ny
          DO i = 1, nx
            kk = nzsoil-k+1
            st_arps(i,j,kk) = arpsgrid%tsoil(i,j,k,0)
            tem1(i,j,kk)    = arpsgrid%tsoil(i,j,k,0)
          END DO
        END DO
      END DO

      CALL exchange_halo(st_arps,nx,ny,nzsoil-1,            &
                         tem1,mxs,mxe,mys,mye,0,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil-1,0,         &
                    srcxs,srcys,dx,dy,                                  &
                    mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,          &
                    tem2,pm1s,pm1e,pm2s,pm2e,nzsoil-1,                  &
                    landsea_arps,landmask_geogrid,285.0,istatus)

      fieldname = 'ST'
      units     = ' '
      descr     = ' '
      CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',     &
                            pm1s,pm1e,pm2s,pm2e,nzsoil-1,istatus)

      IF (wrf_core == 'ARW') THEN
        !
        ! Added soil_layers
        !
        soil_layers(:,:,1) = 300.   ! an extra fake soil layers
        soil_layers(:,:,2) = 255.
        soil_layers(:,:,3) = 100.
        soil_layers(:,:,4) =  28.
        soil_layers(:,:,5) =   7.

        fieldname = 'SOIL_LAYERS'
        units     = ''
        descr     = ''
        CALL metgrid_add_3dsoil_field(soil_layers,fieldname,units,descr,'M', &
                              pm1s,pm1e,pm2s,pm2e,5,istatus)
        CALL metgrid_add_flag("FLAG_SOIL_LAYERS",.TRUE.,geogrid%grid_type,istatus)
      END IF
      metgrid%num_soil_levels = 5

    ELSE
      CALL arpsstop('ERROR: Unknown number of soil layers.',1)
    END IF

!-----------------------------------------------------------------------
!
! OU Soil to be implemented
!
!-----------------------------------------------------------------------

  ELSE IF (soilopt == 2) THEN
    WRITE(6,'(1x,a)') 'Vertical soil interpolation to be implemented.'
    CALL arpsstop('Unsupported soil option.',1)
  ELSE                    ! ARPS soil levels
    WRITE(6,'(1x,a)') 'Unknow soil option.'
    CALL arpsstop('Unsupported soil option.',1)
  END IF

!-----------------------------------------------------------------------
! SOIL moisture
!-----------------------------------------------------------------------

  IF (dbglvl > 0) WRITE(ount,'(1x,a)') '-=- Soil moisture ...'

  IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Fill memory array ...'

  IF (soilopt == 3) THEN ! Keep NAM or GFS soil levels

    SELECT CASE (arpsgrid%nzsoil)
    CASE (6)              ! NAM 4 layer model

      DO k = 2,5
        DO j = 1,ny
          DO i = 1, nx
            kk = nzsoil-k
            st_arps(i,j,kk) = arpsgrid%qsoil(i,j,k,0)
            tem1(i,j,kk)    = arpsgrid%qsoil(i,j,k,0)
          END DO
        END DO
      END DO

      CALL exchange_halo(st_arps,nx,ny,nzsoil-2,                        &
                         tem1,mxs,mxe,mys,mye,0,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil-2,0,         &
                    srcxs,srcys,dx,dy,                                  &
                    mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,          &
                    tem2,pm1s,pm1e,pm2s,pm2e,nzsoil-2,                  &
                    landsea_arps,landmask_geogrid,1.0,istatus)

      fieldname = 'SM'
      units     = ' '
      descr     = ' '
      CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',     &
                            pm1s,pm1e,pm2s,pm2e,nzsoil-2,istatus)

      !
      ! Added separated soil layers
      !
      IF (dbglvl > 2) WRITE(ount,'(3x,a)')     &
        '~~~ Fill output grid from NAM 4 soil layer directly...'

      fieldname = 'SM000010'
      units     = 'kg m-3'
      descr     = 'Soil Moist 0-10 cm below grn layer (Up)'
      CALL metgrid_add_3dsoil_field(tem2(:,:,4),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_SM000010",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'SM010040'
      units     = 'kg m-3'
      descr     = 'Soil Moist 10-40 cm below grn layer'
      CALL metgrid_add_3dsoil_field(tem2(:,:,3),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_SM010040",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'SM040100'
      units     = 'kg m-3'
      descr     = 'Soil Moist 40-100 cm below grn layer'
      CALL metgrid_add_3dsoil_field(tem2(:,:,2),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_SM040100",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'SM100200'
      units     = 'kg m-3'
      descr     = 'Soil Moist 100-200 cm below gr layer'
      CALL metgrid_add_3dsoil_field(tem2(:,:,1),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_SM100200",.TRUE.,geogrid%grid_type,istatus)

    CASE (3)              ! GFS 2 layer model

      DO k = 2,3
        DO j = 1,ny
          DO i = 1, nx
            kk = nzsoil-k+1
            st_arps(i,j,kk) = arpsgrid%qsoil(i,j,k,0)
            tem1(i,j,kk)    = arpsgrid%qsoil(i,j,k,0)
          END DO
        END DO
      END DO

      CALL exchange_halo(st_arps,nx,ny,nzsoil-1,            &
                         tem1,mxs,mxe,mys,mye,0,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil-1,0,         &
                    srcxs,srcys,dx,dy,                                  &
                    mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,          &
                    tem2,pm1s,pm1e,pm2s,pm2e,nzsoil-1,                  &
                    landsea_arps,landmask_geogrid,1.0,istatus)

      fieldname = 'ST'
      units     = ' '
      descr     = ' '
      CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',     &
                            pm1s,pm1e,pm2s,pm2e,nzsoil-1,istatus)

      !
      ! Added separated soil layers
      !

      IF (dbglvl > 2) WRITE(ount,'(3x,a)')     &
        '~~~ Fill output grid from GFS 2 soil layer directly...'

      fieldname = 'SM000010'
      units     = 'kg m-3'
      descr     = 'Soil Moist 0-10 cm below ground layer (Upper) '
      CALL metgrid_add_3dsoil_field(tem2(:,:,2),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_SM000010",.TRUE.,geogrid%grid_type,istatus)

      fieldname = 'SM010200'
      units     = 'kg m-3'
      descr     = 'Soil Moist 100-200 cm below ground layer (Bottom) '
      CALL metgrid_add_3dsoil_field(tem2(:,:,1),fieldname,units,descr,'M', &
                            pm1s,pm1e,pm2s,pm2e,1,istatus)
      CALL metgrid_add_flag("FLAG_SM010200",.TRUE.,geogrid%grid_type,istatus)

    CASE DEFAULT
      CALL arpsstop('ERROR: Unknown number of soil layers.',1)
    END SELECT

  ELSE IF (soilopt == 4) THEN  ! ECMWF 4 layer dataset

    IF ( arpsgrid%nzsoil == 6 ) THEN

      DO k = 2,nzsoil
        DO j = 1,ny
          DO i = 1, nx
            kk = nzsoil-k+1
            st_arps(i,j,kk) = arpsgrid%qsoil(i,j,k,0)
            tem1(i,j,kk)    = arpsgrid%qsoil(i,j,k,0)
          END DO
        END DO
      END DO

      CALL exchange_halo(st_arps,nx,ny,nzsoil-1,                        &
                         tem1,mxs,mxe,mys,mye,0,istatus)

      IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

      CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil-1,0,         &
                    srcxs,srcys,dx,dy,                                  &
                    mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,          &
                    tem2,pm1s,pm1e,pm2s,pm2e,nzsoil-1,                  &
                    landsea_arps,landmask_geogrid,1.0,istatus)

      fieldname = 'SM'
      units     = ' '
      descr     = ' '
      CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',     &
                            pm1s,pm1e,pm2s,pm2e,nzsoil-1,istatus)

    ELSE
      CALL arpsstop('ERROR: Unknown number of soil layers.',1)
    END IF

  !
  !
  ! ARPS traditional 2-layer model
  !
  !
  ELSE IF (soilopt == 1) THEN

    tem1(1:nx,1:ny,1:nzsoil) = arpsgrid%qsoil(:,:,:,0)

    CALL exchange_halo(arpsgrid%qsoil(:,:,:,0),nx,ny,nzsoil,            &
                       tem1,mxs,mxe,mys,mye,0,istatus)

    IF (dbglvl > 2) WRITE(ount,'(3x,a)') '~~~ Doing interpolation ...'

    CALL interp3d_mask(tem1,mxs,mxe,mys,mye,nx,ny,nzsoil,0,             &
                  srcxs,srcys,dx,dy,                                    &
                  mdesx,mdesy,mm1s,mm1e,mm2s,mm2e,do_interp,            &
                  tem2,pm1s,pm1e,pm2s,pm2e,nzsoil,                      &
                  landsea_arps,landmask_geogrid,1.0,istatus)

    fieldname = 'SM'
    units     = ' '
    descr     = ' '
    CALL metgrid_add_3dsoil_field(tem2,fieldname,units,descr,'M',       &
                          pm1s,pm1e,pm2s,pm2e,arpsgrid%nzsoil,istatus)

    !
    ! Fill separate layers
    !

    IF (dbglvl > 2) WRITE(ount,'(3x,a)')     &
      '~~~ Fill output grid from ARPS 2 layer soil model directly...'


    fieldname = 'SM000010'
    units     = 'kg m-3'
    descr     = 'Soil Moist 0-10 cm below grn layer (Upper) '
    CALL metgrid_add_3dsoil_field(tem2(:,:,1),fieldname,units,descr,'M', &
                          pm1s,pm1e,pm2s,pm2e,1,istatus)
    CALL metgrid_add_flag("FLAG_SM000010",.TRUE.,geogrid%grid_type,istatus)

    fieldname = 'SM010200'
    units     = 'kg m-3'
    descr     = 'Soil Moist 10-200 cm below grn layer (Bottom) '
    CALL metgrid_add_3dsoil_field(tem2(:,:,2),fieldname,units,descr,'M', &
                          pm1s,pm1e,pm2s,pm2e,1,istatus)
    CALL metgrid_add_flag("FLAG_SM010200",.TRUE.,geogrid%grid_type,istatus)

!
!
! OU Soil to be implemented
!
!

  ELSE IF (soilopt == 2) THEN
    WRITE(6,'(1x,a)') 'Vertical soil interpolation to be implemented.'
    CALL arpsstop('Unsupported soil option.',1)
  ELSE                    ! ARPS soil levels
    WRITE(6,'(1x,a)') 'Unknow soil option.'
    CALL arpsstop('Unsupported soil option.',1)
  END IF

!-----------------------------------------------------------------------
! END OF SOIL VARIBLES
!-----------------------------------------------------------------------

  DEALLOCATE(st_arps)
  IF (ALLOCATED(soil_layers)) DEALLOCATE(soil_layers)

  RETURN
END SUBROUTINE interpolate_soil_from_arps
