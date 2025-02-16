
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE getwrfd                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getwrfd(fHndl,io_form,multifile,ncmprx,ncmpry,itime,timestr, &
                   nx_ext,ny_ext,nz_ext,nzsoil_ext,                     &
                   iproj_ext,scale_ext,trlon_ext,latnot_ext,            &
                   ctrlat_ext,ctrlon_ext,dx_ext,dy_ext,x0_ext,y0_ext,   &
                   sfcphys,sfclay,num_scalar,                           &
                   lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext, &
                   zp_ext,hgt_ext,zpsoil_ext, p_ext,t_ext,              &
                   u_ext,v_ext,w_ext,                                   &
                   qv_ext,qscalar_ext,tke_ext,                          &
                   tsoil_ext,qsoil_ext,wetcanp_ext,                     &
                   snowdpth_ext,trn_ext,soiltyp_ext,vegtyp_ext,veg_ext, &
                   tem1_ext,tem2d1,tem2d2,tem2di,istatus)

!------------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in WRF data in one time level (itime) from the opened NetCDF
!  file (fHndl) and convert those data to ARPS units if needed.
!
!------------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (12/29/2003)
!
!  MODIFICATION HISTORY:
!
!  2004-02-17 Richard Carpenter
!  Revised vegetation table.
!
!  2004-02-17 Gene Bassett
!  Revised snow and precip parameters.
!
!  03/25/2004 Yunheng Wang
!  Rewrote for message passing mode.
!
!  03/20/2010 Fanyou Kong
!  Add new WRF fields of hourly maximum (from NSSL WRF mod)
!------------------------------------------------------------------------
  IMPLICIT NONE

!------------------------------------------------------------------------
!
!  Include files
!
!------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!------------------------------------------------------------------------
!
!  Variable declaration
!
!------------------------------------------------------------------

  INTEGER,      INTENT(IN)  :: ncmprx, ncmpry
  INTEGER,      INTENT(IN)  :: fHndl(ncmprx,ncmpry)
  INTEGER,      INTENT(IN)  :: io_form      ! 5 - PHDF5
                                            ! 7 - NetCDF
  LOGICAL,      INTENT(IN)  :: multifile
  INTEGER,      INTENT(IN)  :: itime   ! Recorder number
  CHARACTER(*), INTENT(IN)  :: timestr
  INTEGER, INTENT(IN)  :: nx_ext       !
  INTEGER, INTENT(IN)  :: ny_ext       !
  INTEGER, INTENT(IN)  :: nz_ext       !
  INTEGER, INTENT(IN)  :: nzsoil_ext   ! WRF soil layers
                                       ! Note: it is the actual WRF soil layers
                                       !       plus the data at surface.

  INTEGER, INTENT(IN)  :: iproj_ext    ! Map projection option
                                       ! = 1, polar projection;
                                       ! = 2, Lambert projection;
                                       ! = 3, Mercator projection.
  REAL,    INTENT(IN)  :: scale_ext    ! Map scale factor (should be 1.0)
  REAL,    INTENT(IN)  :: trlon_ext    ! True longitude
  REAL,    INTENT(IN)  :: latnot_ext(2)! True latitude
  REAL,    INTENT(IN)  :: ctrlon_ext, ctrlat_ext
  REAL,    INTENT(IN)  :: dx_ext,dy_ext
  INTEGER, INTENT(IN)  :: sfcphys,sfclay        ! WRF surface physics option
  INTEGER, INTENT(IN)  :: num_scalar     ! number of WRF scalars that can be mapped to ARPS variables

  REAL,    INTENT(OUT) :: x0_ext
  REAL,    INTENT(OUT) :: y0_ext
  REAL,    INTENT(OUT) :: lat_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lon_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: latu_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lonu_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: latv_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lonv_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: zp_ext  (nx_ext,ny_ext,nz_ext)
                                       ! physical height (m) at W points
  REAL,    INTENT(OUT) :: hgt_ext (nx_ext,ny_ext,nz_ext)
                                       ! Physical height (m) at scalar points
  REAL,    INTENT(OUT) :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: p_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: t_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: qv_ext  (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: u_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: v_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: w_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: qscalar_ext(nx_ext,ny_ext,nz_ext,nscalar)
  REAL,    INTENT(OUT) :: tke_ext (nx_ext,ny_ext,nz_ext)

  REAL,    INTENT(OUT) :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)

  REAL,    INTENT(OUT) :: trn_ext     (nx_ext,ny_ext)
  INTEGER, INTENT(OUT) :: soiltyp_ext (nx_ext,ny_ext)
  INTEGER, INTENT(OUT) :: vegtyp_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: veg_ext     (nx_ext,ny_ext)

  REAL,    INTENT(INOUT) :: tem1_ext (nx_ext,  ny_ext,  nz_ext)
  REAL,    INTENT(INOUT) :: tem2d1   (nx_ext,  ny_ext)
  REAL,    INTENT(INOUT) :: tem2d2   (nx_ext,  ny_ext)
  INTEGER, INTENT(INOUT) :: tem2di   (nx_ext,  ny_ext)
  INTEGER, INTENT(OUT)   :: istatus

!------------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER  :: nxlg_ext, nylg_ext
  INTEGER  :: nxd, nyd, nzd            ! dimensions in data file
  INTEGER  :: nxd_stag, nyd_stag, nzd_stag

  INTEGER  :: iproj_orig
  REAL     :: scale_orig
  REAL     :: latnot_orig(2)
  REAL     :: trlon_orig
  REAL     :: x0_orig, y0_orig
  REAL     :: xsub0,ysub0

  REAL, ALLOCATABLE :: x_ext(:)
  REAL, ALLOCATABLE :: y_ext(:)
  REAL, ALLOCATABLE :: xsc_ext(:)
  REAL, ALLOCATABLE :: ysc_ext(:)

  REAL, ALLOCATABLE :: temlg1(:,:,:)
  REAL, ALLOCATABLE :: temlg2(:,:,:)

  REAL, PARAMETER  :: epsl = 10E-2
  REAL, PARAMETER  :: grav = 9.81
  REAL, PARAMETER  :: rd   = 287.0
  REAL, PARAMETER  :: cp   = 1004.0
                      ! Specific heat of dry air at constant pressure
                      ! (m**2/(s**2*K)).
  REAL, PARAMETER  :: p0   = 1.0E5
                      ! Surface reference pressure, is 100000 Pascal.
  REAL, PARAMETER  :: rddcp= rd/cp
  REAL, PARAMETER  :: govrd = grav/rd

!
! These two tables are needed to convert WRF soil and vegetation index
! to those used in ARPS, althought they may not correspond each other
! one to one.
!
! WRF uses 16 soil categories and 24-category (USGS) vegetation
! ARPS uses 13-category soil  and 14-category (ND) vegetation
!
! The following two tables were provided by Jerry Brotzge on Oct. 20, 2003
!
  INTEGER, PARAMETER :: soil_table(17) = (/ 1, 2, 3, 4, 4,              &
                                            5, 6, 7, 8, 9,              &
                                           10,11, 6,13, 1,              &
                                            2, 2/)

! WDT RLC 2004-02-12 Changed veg_table(1)/Urban from 7/EvgrnForest to 1/Desert
! WDT RLC 2004-02-12 Changed veg_table(25)/Playa from 10/Cultiv to 1/Desert
! WDT RLC 2004-02-12 Added veg_table(26:27)
! WDT RLC 2004-02-16 Not accepting CAPS changes:
                     ! WYH - (01.16.2004) veg_table(18) = 11 =>  8
                     !                    veg_table(19) = 1  => 13
  !INTEGER, PARAMETER :: veg_table(25)  = (/ 7,10,10,10,10,              &
  !                                          5, 3,12, 4,12,              &
  !                                          6, 6, 7, 7, 6,              &
  !                                         14,11,11, 1, 2,              &
  !                                          2, 2, 2, 9,10/)
  INTEGER, PARAMETER :: veg_table(27)  = (/ 1,10,10,10,10,              &
                                            5, 3,12, 4,12,              &
                                            6, 6, 7, 7, 6,              &
                                           14,11,11, 1, 2,              &
                                            2, 2, 2, 9, 1,              &
                                            1, 1/)
  !
  ! When the variable does not exist in the NetCDF file, this constant
  ! should be return in variable "istatus" of the calling subroutine.
  !
  INTEGER, PARAMETER :: VAR_NOTEXIST = -1
  INTEGER, PARAMETER :: WRF_REAL     = 104
  INTEGER, PARAMETER :: WRF_INTEGER  = 106

  CHARACTER(LEN=40)  :: qnames_wrf(nscalar)

  REAL               :: tvbot, tvtop, tvbar
  REAL               :: xctr,yctr
  INTEGER            :: i,j,k ,nq

  REAL               :: dzs(nzsoil_ext) ! WRF soil depths
  REAL               :: znw(nz_ext)
  LOGICAL            :: warned


  LOGICAL, PARAMETER :: SPRING2010 = .FALSE.
  LOGICAL, PARAMETER :: SPRING2011 = .FALSE.
  LOGICAL, PARAMETER :: SPRING2012 = .FALSE.
  LOGICAL, PARAMETER :: SPRING2013 = .FALSE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  warned = .FALSE.

  ! remember the map projection paramters which should be restored
  ! before return

  CALL getmapr(iproj_orig,scale_orig,latnot_orig,                       &
               trlon_orig,x0_orig,y0_orig)

!-----------------------------------------------------------------------
!
! Allocate working arrays
!
!-----------------------------------------------------------------------

  nxlg_ext = (nx_ext-1)*nproc_x + 1
  nylg_ext = (ny_ext-1)*nproc_y + 1

  nzd      = nz_ext - 1
  nzd_stag = nz_ext

  IF (mp_opt > 0 .AND. multifile) THEN
    nxd = (nx_ext - 1)/ncmprx
    nyd = (ny_ext - 1)/ncmpry
    nxd_stag = nxd
    nyd_stag = nyd
    IF (loc_x == nproc_x) nxd_stag = nxd + 1
    IF (loc_y == nproc_y) nyd_stag = nyd + 1

    ALLOCATE(temlg2(1,1,1),   STAT = istatus)  ! do not use it
  ELSE
    nxd = nxlg_ext - 1
    nyd = nylg_ext - 1
    nxd_stag = nxlg_ext
    nyd_stag = nylg_ext

    ALLOCATE(temlg2(nxlg_ext,nylg_ext,nz_ext),   STAT = istatus)
  END IF

  ALLOCATE(temlg1(nxd_stag, nyd_stag, nzd_stag), STAT = istatus)

!-----------------------------------------------------------------------
!
! Set up WRF map projection
!
!-----------------------------------------------------------------------

  ALLOCATE(x_ext(nx_ext),   STAT = istatus)
  ALLOCATE(y_ext(ny_ext),   STAT = istatus)
  ALLOCATE(xsc_ext(nx_ext), STAT = istatus)
  ALLOCATE(ysc_ext(ny_ext), STAT = istatus)

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,ctrlat_ext,ctrlon_ext,xctr,yctr)
  x0_ext=xctr - 0.5*nproc_x*(nx_ext-1)*dx_ext
  y0_ext=yctr - 0.5*nproc_y*(ny_ext-1)*dy_ext
  CALL setorig(1,x0_ext,y0_ext)

  xsub0 = dx_ext*(nx_ext-1)*(loc_x-1)
  ysub0 = dy_ext*(ny_ext-1)*(loc_y-1)

  DO i=1,nx_ext
    x_ext(i)= xsub0 + (i-1)*dx_ext
  END DO
  DO j=1,ny_ext
    y_ext(j)= ysub0 + (j-1)*dy_ext
  END DO

  DO i=1,nx_ext-1
    xsc_ext(i)=0.5*(x_ext(i)+x_ext(i+1))
  END DO
  xsc_ext(nx_ext)=2.*xsc_ext(nx_ext-1)-xsc_ext(nx_ext-2)
  DO j=1,ny_ext-1
    ysc_ext(j)=0.5*(y_ext(j)+y_ext(j+1))
  END DO
  ysc_ext(ny_ext)=2.*ysc_ext(ny_ext-1)-ysc_ext(ny_ext-2)

  CALL xytoll(nx_ext,ny_ext,xsc_ext,ysc_ext, lat_ext,lon_ext)
  CALL xytoll(nx_ext,ny_ext,  x_ext,ysc_ext,latu_ext,lonu_ext)
  CALL xytoll(nx_ext,ny_ext,xsc_ext,  y_ext,latv_ext,lonv_ext)

  DEALLOCATE(x_ext,y_ext)
  DEALLOCATE(xsc_ext,ysc_ext)

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LU_INDEX',WRF_REAL,'XY ',' ',       &
  !                   'west_east','sout_north',' ',nx_ext,ny_ext,1,      &
  !                   nxd,nyd,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
!  Horizontal wind
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'U','X',                             &
                   'west_east_stag','south_north','bottom_top',         &
                   nx_ext,    ny_ext,    nz_ext,  u_ext,                &
                   nxd_stag,  nyd,       nzd,     temlg1,               &
                   nxlg_ext,  nylg_ext,  nz_ext,  temlg2, istatus)

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'V','Y',                             &
                   'west_east','south_north_stag','bottom_top',         &
                   nx_ext,    ny_ext,    nz_ext,  v_ext,                &
                   nxd,       nyd_stag,  nzd,     temlg1,               &
                   nxlg_ext,  nylg_ext,  nz_ext,  temlg2, istatus)

!-----------------------------------------------------------------------
!
! Vertical velocity
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'W','Z',                             &
                   'west_east','south_north','bottom_top_stag',         &
                   nx_ext,    ny_ext,    nz_ext,  w_ext,                &
                   nxd,       nyd,       nzd_stag,temlg1,               &
                   nxlg_ext,  nylg_ext,  nz_ext,  temlg2, istatus)

!-----------------------------------------------------------------------
!
! Physical height
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'PH','Z',                            &
                   'west_east','south_north','bottom_top_stag',         &
                   nx_ext,    ny_ext,    nz_ext,  zp_ext,               &
                   nxd,       nyd,       nzd_stag,temlg1,               &
                   nxlg_ext,  nylg_ext,  nz_ext,  temlg2,istatus)

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'PHB','Z',                           &
                   'west_east','south_north','bottom_top_stag',         &
                   nx_ext,   ny_ext,    nz_ext,  tem1_ext,              &
                   nxd,      nyd,       nzd_stag,temlg1,                &
                   nxlg_ext, nylg_ext,  nz_ext,  temlg2,istatus)

  zp_ext(:,:,:) = ( zp_ext(:,:,:) + tem1_ext(:,:,:) )/grav

!
!  Move zp_ext field onto the scalar layers.
!
  DO k=1,nz_ext-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        hgt_ext(i,j,k)=0.5*(zp_ext(i,j,k)+zp_ext(i,j,k+1))
      END DO
    END DO
  END DO
  DO j=1,ny_ext
    DO i=1,nx_ext
      hgt_ext(i,j,nz_ext)=(2.*hgt_ext(i,j,nz_ext-1))                    &
                             -hgt_ext(i,j,nz_ext-2)
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Read potential temperature and convert to temperature
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'T','',                              &
                   'west_east','south_north','bottom_top',              &
                   nx_ext,    ny_ext,    nz_ext,  t_ext,                &
                   nxd,       nyd,       nzd,     temlg1,               &
                   nxlg_ext,  nylg_ext,  nz_ext,  temlg2,istatus)

!-----------------------------------------------------------------------
!
! Dummy arrays
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'MU',                                &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'MUB',                               &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'NEST_POS',                          &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
! Pressure
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'P','',                              &
                   'west_east','south_north','bottom_top',              &
                   nx_ext,    ny_ext,    nz_ext,  p_ext,                &
                   nxd,       nyd,       nzd,     temlg1,               &
                   nxlg_ext,  nylg_ext,  nz_ext,  temlg2,istatus)

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                   timestr,itime,1,'PB','',                             &
                   'west_east','south_north','bottom_top',              &
                   nx_ext,   ny_ext,   nz_ext,  tem1_ext,               &
                   nxd,      nyd,      nzd,     temlg1,                 &
                   nxlg_ext, nylg_ext, nz_ext,  temlg2,istatus)

  p_ext(:,:,:) = p_ext(:,:,:) + tem1_ext(:,:,:)

!-----------------------------------------------------------------------
!
! Convert potential temperature to temperature
!
!-----------------------------------------------------------------------

  t_ext(:,:,:) = (t_ext(:,:,:) + 300.) * ((p_ext(:,:,:)/p0)**rddcp)

!-----------------------------------------------------------------------
!
! Dummy arrays
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SR',                                &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'POTEVP',                            &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SNOPCX',                            &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SOILTB',                            &
  !                   WRF_REAL,'XY ',' ','west_east','south_north',' ',  &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'FNM',WRF_REAL,'Z',' ',              &
  !                   'bottom_top',' ',' ',                              &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'FNP',WRF_REAL,'Z',' ',              &
  !                   'bottom_top',' ',' ',                              &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'RDNW',WRF_REAL,'Z',' ',             &
  !                   'bottom_top',' ',' ',                              &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'RDN',                               &
  !                   WRF_REAL,'Z',' ','bottom_top',' ',' ',             &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'DNW',                               &
  !                   WRF_REAL,'Z',' ','bottom_top',' ',' ',             &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'DN',                                &
  !                   WRF_REAL,'Z',' ','bottom_top',' ',' ',             &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'ZNU',                               &
  !                   WRF_REAL,'Z',' ','bottom_top',' ',' ',             &
  !                   nz_ext,1,1,nzd,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'ZNW',                               &
  !                   WRF_REAL,'Z','Z','bottom_top_stag',' ',' ',        &
  !                   nz_ext,1,1,nzd_stag,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'CFN',                               &
  !                   WRF_REAL,'0',' ',' ',' ',' ',                      &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'CFN1',                              &
  !                   WRF_REAL,'0',' ',' ',' ',' ',                      &
  !                   1,1,1,1,1,1,temlg1,istatus)

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'RDX',WRF_REAL,'0',' ',              &
  !                   ' ',' ',' ',                                       &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'RDY',WRF_REAL,'0',' ',              &
  !                   ' ',' ',' ',                                       &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'RESM',WRF_REAL,'0',' ',             &
  !                   ' ',' ',' ',                                       &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'ZETATOP',WRF_REAL, '0',' ',         &
  !                   ' ',' ',' ',                                       &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'CF1',WRF_REAL,'0',' ',              &
  !                   ' ',' ',' ',                                       &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'CF2',WRF_REAL,'0',' ',              &
  !                   ' ',' ',' ',                                       &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,timestr,itime,'CF3',       &
  !                   WRF_REAL,'0',' ',' ',' ',' ',                      &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,timestr,itime,'ITIMESTEP', &
  !                   WRF_INTEGER,'0',' ',' ',' ',' ',                   &
  !                   1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,timestr,itime,'XTIME',     &
  !                   WRF_REAL,'0',' ',' ',' ',' ',                      &
  !                   1,1,1,1,1,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
! Read water vapor mixing ratio and Convert to specific humidity
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,timestr,itime,1,'QVAPOR','', &
                   'west_east','south_north','bottom_top',              &
                   nx_ext,   ny_ext,   nz_ext, qv_ext,                  &
                   nxd,      nyd,      nzd,    temlg1,                  &
                   nxlg_ext, nylg_ext, nz_ext, temlg2, istatus)

  qv_ext(:,:,:) = qv_ext(:,:,:)/(1+qv_ext(:,:,:))

  WHERE(qv_ext < 0.0) qv_ext = 0.0
  WHERE(qv_ext > 1.0) qv_ext = 1.0

!-----------------------------------------------------------------------
!
!  Make top and bottom mass fields via hydrostatic extrapolation.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny_ext
    DO i=1,nx_ext

      t_ext(i,j,1)=2.*t_ext(i,j,2)-t_ext(i,j,3)
      tvbot=t_ext(i,j,1) * ( 1.0 + 0.608*qv_ext(i,j,1) )
      tvtop=t_ext(i,j,2) * ( 1.0 + 0.608*qv_ext(i,j,2) )
      tvbar=0.5*(tvtop+tvbot)
      p_ext(i,j,1)= p_ext(i,j,2)                                        &
                   *EXP(govrd*(hgt_ext(i,j,2)-hgt_ext(i,j,1))/tvbar)

      t_ext(i,j,nz_ext)=2.*t_ext(i,j,nz_ext-1)-t_ext(i,j,nz_ext-2)
      tvbot=t_ext(i,j,nz_ext-1)*(1.0 + 0.608*qv_ext(i,j,nz_ext-1))
      tvtop=t_ext(i,j,nz_ext)  *(1.0 + 0.608*qv_ext(i,j,nz_ext))
      tvbar=0.5*(tvtop+tvbot)
      p_ext(i,j,nz_ext)= p_ext(i,j,nz_ext-1)                            &
                        *EXP(govrd*(hgt_ext(i,j,nz_ext-1)-              &
                             hgt_ext(i,j,nz_ext))/tvbar)
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Hydrometeor vapor mixing ratio
!
!-----------------------------------------------------------------------

  IF (P_QC > 0)    qnames_WRF(P_QC) = 'QCLOUD'
  IF (P_QR > 0)    qnames_WRF(P_QR) = 'QRAIN'
  IF (P_QI > 0)    qnames_WRF(P_QI) = 'QICE'
  IF (P_QS > 0)    qnames_WRF(P_QS) = 'QSNOW'
  IF (P_QG > 0)    qnames_WRF(P_QG) = 'QGRAUP'
  IF (P_QH > 0)    qnames_WRF(P_QH) = 'QHAIL'

  IF (P_NC > 0)    qnames_WRF(P_NC) = 'QNCLOUD'
  IF (P_NR > 0)    qnames_WRF(P_NR) = 'QNRAIN'
  IF (P_NI > 0)    qnames_WRF(P_NI) = 'QNICE'
  IF (P_NS > 0)    qnames_WRF(P_NS) = 'QNSNOW'
  IF (P_NG > 0)    qnames_WRF(P_NG) = 'QNGRAUPEL'
  IF (P_NH > 0)    qnames_WRF(P_NH) = 'QNHAIL'

  DO nq = num_scalar+1, nscalar      ! Extra arrays to be read
    qnames_WRF(nq) = qnames(nq)
  END DO

  DO nq = 1,nscalar

    CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,             &
                     timestr,itime,1,TRIM(qnames_wrf(nq)),'',           &
                     'west_east','south_north','bottom_top',            &
                     nx_ext,   ny_ext,   nz_ext, qscalar_ext(:,:,:,nq), &
                     nxd,      nyd,      nzd,    temlg1,                &
                     nxlg_ext, nylg_ext, nz_ext, temlg2,istatus)

    !WHERE(qscalar_ext < 0.0) qc_ext = 0.0
    !WHERE(qc_ext > 1.0) qc_ext = 1.0

  END DO

!  IF ( P_QT > 0) THEN
!    CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,           &
!                     timestr,itime,'QT',                                &
!                     WRF_REAL,'XYZ',' ','west_east','south_north','bottom_top',  &
!                     nx_ext,ny_ext,nz_ext,nxd,nyd,nzd,temlg1,istatus)
!  END IF
!
!  IF ( P_QNN > 0) THEN
!    CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,           &
!                     timestr,itime,'QNCCN',                             &
!                     WRF_REAL,'XYZ',' ','west_east','south_north','bottom_top',  &
!                     nx_ext,ny_ext,nz_ext,nxd,nyd,nzd,temlg1,istatus)
!  END IF

!-----------------------------------------------------------------------
!
! Read TKE_PBL
!
!-----------------------------------------------------------------------

  CALL  get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,timestr,itime,1,'TKE_PBL','', &
                   'west_east','south_north','bottom_top',              &
                   nx_ext,   ny_ext,   nz_ext, tke_ext,                 &
                   nxd,      nyd,      nzd,    temlg1,                  &
                   nxlg_ext, nylg_ext, nz_ext, temlg2, istatus)

  IF (sfclay == 5) tke_ext(:,:,:) = tke_ext(:,:,:)*0.5
  !
  ! Get Land mask
  !
  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'LANDMASK','',                        &
                  'west_east','south_north',                            &
                  nx_ext,ny_ext,tem1_ext(:,:,1),nxd,nyd,temlg1,         &
                  nxlg_ext,nylg_ext,temlg2,istatus)

!-----------------------------------------------------------------------
!
! Soil temperature
!
! NOTE: the shape of "temlg1" has been changed below.
!
!-----------------------------------------------------------------------

  CALL get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'TSLB','Z',                           &
                  'west_east','south_north','soil_layers_stag',         &
                  nx_ext,   ny_ext,   nzsoil_ext-1,tsoil_ext(:,:,2),    &
                  nxd,      nyd,      nzsoil_ext-1,temlg1,              &
                  nxlg_ext, nylg_ext, nzsoil_ext-1,temlg2, istatus)
!
!-----------------------------------------------------------------------
!
! Soil depth
!
!-----------------------------------------------------------------------

  CALL get_wrf_1d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,'ZS','Z',                               &
                  'soil_layers_stag',nzsoil_ext-1,dzs,nzsoil_ext-1,istatus)

  DO k = 2, nzsoil_ext
    zpsoil_ext(:,:,k) = dzs(k-1)
  END DO
  zpsoil_ext(:,:,1) = 0.0

  !
  ! Do some checking
  !
  IF(sfcphys == 1) THEN
    IF(nzsoil_ext /= 6) THEN
      WRITE(6,'(/a/)')'================= WARNING ======================'
      WRITE(6,'(2(a,I2),a/,a)') 'Number of soil layers (',nzsoil_ext-1, &
                 ') is not the WRF default for sf_surface_physics =',   &
                 sfcphys,'.','Expecting soil_layers_stag = 5.'
      WRITE(6,'(/a/)')'================= WARNING ======================'
    END IF
  ELSE IF(sfcphys == 2) THEN
    IF(nzsoil_ext /= 5) THEN
      WRITE(6,'(/a/)')'================= WARNING ======================'
      WRITE(6,'(2(a,I2),a/,a)') 'Number of soil layers (',nzsoil_ext-1, &
                 ') is not the WRF default for sf_surface_physics =',   &
                 sfcphys,'.',' Expecting soil_layers_stag = 4.'
      WRITE(6,'(/a/)')'================= WARNING ======================'
    END IF
  ELSE IF(sfcphys == 3) THEN
    !
    ! When do interpolation, we may need to discard the first layer
    ! because of duplication
    !
    IF(nzsoil_ext /= 7) THEN
      WRITE(6,'(/a/)')'================= WARNING ======================'
      WRITE(6,'(2(a,I2),a/,a)') 'Number of soil layers (',nzsoil_ext-1, &
                 ') is not the WRF default for sf_surface_physics = ',  &
                 sfcphys,'.','Expecting soil_layers_stag = 6.'
      WRITE(6,'(/a/)')'================= WARNING ======================'
    END IF
  ELSE IF(sfcphys == 7) THEN
    IF(nzsoil_ext /= 3) THEN
      WRITE(6,'(/a/)')'================= WARNING ======================'
      WRITE(6,'(2(a,I2),a/,a)') 'Number of soil layers (',nzsoil_ext-1, &
                 ') is not the WRF default for sf_surface_physics =',   &
                 sfcphys,'.',' Expecting soil_layers_stag = 2.'
      WRITE(6,'(/a/)')'================= WARNING ======================'
    END IF
  ELSE
    WRITE(6,'(/a/)')'================= WARNING ======================'
    WRITE(6,'(a,I2,a/,a)') 'sf_surface_physics has wrong number = ',    &
          sfcphys, '.','Must be 1, 2 or 3. Please check your data file.'
    WRITE(6,'(/a/)')'================= WARNING ======================'
  END IF

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'DZS',WRF_REAL,'Z','Z',              &
  !                   'soil_layers_stag',' ',' ',                        &
  !                   nzsoil_ext-1,1,1,nzsoil_ext-1,1,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
! Soil moisture  (? fraction (m*3/m*3)
!
!-----------------------------------------------------------------------
!
  CALL get_wrf_3d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'SMOIS','Z',                          &
                  'west_east','south_north','soil_layers_stag',         &
                  nx_ext,   ny_ext,  nzsoil_ext-1, qsoil_ext(:,:,2),    &
                  nxd,      nyd,     nzsoil_ext-1, temlg1,              &
                  nxlg_ext, nylg_ext,nzsoil_ext-1, temlg2, istatus)

!-----------------------------------------------------------------------
!
! Dummy arrays
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SH2O',WRF_REAL,'XZY','Z',           &
  !                   'west_east','south_north','soil_layers_stag',      &
  !                   nx_ext,ny_ext,nzsoil_ext-1,                        &
  !                   nxd,nyd,nzsoil_ext-1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'XICE',WRF_REAL,'XY',' ',            &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SFROFF',WRF_REAL,'XY',' ',          &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'UDROFF',WRF_REAL,'XY',' ',          &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
!  Vegatation type
!
!-----------------------------------------------------------------------
!
  CALL get_wrf_2di(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                  timestr,itime,1,'IVGTYP','',                          &
                  'west_east','south_north',nx_ext,ny_ext,tem2di,       &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
                                         ! temlg? used as INTEGER inside

  IF (ANY(tem2di<=0)) THEN
    WRITE(6,'(/a/)')'================= WARNING ======================'
    WRITE(6,'(a)') 'Find 0s for WRF IVGTYP.'
    WRITE(6,'(/a/)')'================= WARNING ======================'
    vegtyp_ext(:,:) = 0
  ELSE
    DO j = 1, ny_ext
      DO i = 1,nx_ext
        vegtyp_ext(i,j) = veg_table(tem2di(i,j))  ! convert to ARPS categories
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Soil type
!
!-----------------------------------------------------------------------
!
  CALL get_wrf_2di(fHndl,io_form,multifile,ncmprx,ncmpry,               &
                  timestr,itime,1,'ISLTYP','',                          &
                  'west_east','south_north',nx_ext,ny_ext,tem2di,       &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
                                        ! temlg? used as INTEGER inside

  IF (ANY(tem2di <= 0)) THEN
    WRITE(6,'(/a/)')'================= WARNING ======================'
    WRITE(6,'(a)') 'Find 0s for WRF ISLTYP.'
    WRITE(6,'(/a/)')'================= WARNING ======================'
    vegtyp_ext(:,:) = 0
  ELSE
    DO j = 1, ny_ext
      DO i = 1,nx_ext
        soiltyp_ext(i,j) = soil_table(tem2di(i,j)) ! convert to ARPS categories
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Vegetation Fraction
!
!-----------------------------------------------------------------------
!
  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'VEGFRA','',                          &
                  'west_east','south_north',nx_ext,ny_ext,veg_ext,      &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'GRDFLX',WRF_REAL,'XY',' ',          &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
!
!-----------------------------------------------------------------------
!
! Accumulated snow depth (meter)
!
!-----------------------------------------------------------------------

  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'SNOW','',                            &
                  'west_east','south_north',nx_ext,ny_ext,snowdpth_ext, &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  ! Convert water equiv. of accum. snow depth (kg/m**2) to meters
  ! (where 1 meter liquid water is set equivqlent to 10 meters snow).
  !        0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)
  !snowdpth_ext(i,j) = 0.01 * snowdpth_ext(i,j)

  DO j = 1,ny_ext
    DO i = 1,nx_ext
      snowdpth_ext(i,j) = 0.01*snowdpth_ext(i,j)
    END DO
  END DO

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SNOWH',WRF_REAL,'XY',' ',           &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'RHOSN',WRF_REAL,'XY',' ',           &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
! Canopy water amount (What should be the units in WRF? "kg m{-2}")
! (mm = kg m{-2}, magnitude 2-3 mm according to Jerry)
!
!-----------------------------------------------------------------------

  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'CANWAT','',                          &
                  'west_east','south_north',nx_ext,ny_ext,wetcanp_ext,  &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

!
! NOTE: tem1_ext, tem2d1, tem2d2 used as temporary arrays below
!       tem1_ext(:,:,1)     - land mask
!       tem2d1     - Surface skin temperature
!       tem2d2     - SST
!

!
! Get SST
!
  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'SST','',                             &
                  'west_east','south_north',nx_ext,ny_ext,tem2d2,       &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

!-----------------------------------------------------------------------
!
! Dummy arrays
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'MAPFAC_M',WRF_REAL,'XY',' ',        &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'MAPFAC_U',WRF_REAL,'XY','X',        &
  !                   'west_east_stag','south_north',' ',                &
  !                   nx_ext,ny_ext,1,nxd_stag,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'MAPFAC_V',WRF_REAL,'XY','Y',        &
  !                   'west_east','south_north_stag',' ',                &
  !                   nx_ext,ny_ext,1,nxd,nyd_stag,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'F',WRF_REAL,'XY',' ',               &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'E',WRF_REAL,'XY',' ',               &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SINALPHA',WRF_REAL,'XY',' ',        &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'COSALPHA',WRF_REAL,'XY',' ',        &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
!  Terrain height
!
!-----------------------------------------------------------------------
!
  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'HGT','',                             &
                  'west_east','south_north',nx_ext,ny_ext,trn_ext,      &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
!
! Get surface skin temperature
!
  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'TSK','',                             &
                  'west_east','south_north',nx_ext,ny_ext,tem2d1,       &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

!
! Adjust tsoil_ext and qsoil_ext
!
! NOTE: tem1_ext, tem2d1, tem2d2 used as temporary arrays below
!       tem1_ext(:,:,1)     - land mask
!       tem2d1     - Surface skin temperature
!       tem2d2     - SST
!
  DO j = 1, ny_ext
    DO i = 1, nx_ext
      IF(tem1_ext(i,j,1) < 0.5) THEN         ! Water
        DO k = 1,nzsoil_ext
          tsoil_ext(i,j,k) = tem2d2(i,j)     ! SST
        END DO
      ELSE                                   ! Land
        tsoil_ext(i,j,1) = tem2d1(i,j)       ! TSK
      END IF
    END DO
  END DO

  DO j = 1, ny_ext
    DO i = 1, nx_ext
      IF(tem1_ext(i,j,1) < 0.5) THEN           ! Water
        DO k = 1,nzsoil_ext
          qsoil_ext(i,j,k) = 1.0
        END DO
      ELSE                                     ! Land
        qsoil_ext(i,j,1) = qsoil_ext(i,j,2)    ! Assumed to be same
                                               ! as first below ground
                                               ! level.
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Dummy arrays
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'P_TOP',WRF_REAL,'0',' ',            &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LL_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UL_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UR_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LR_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LL_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UL_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UR_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LR_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LL_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UL_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UR_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LR_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LL_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UL_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_UR_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LAT_LR_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LL_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UL_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UR_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LR_T',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LL_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UL_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UR_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LR_U',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LL_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UL_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UR_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LR_V',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LL_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UL_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_UR_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LON_LR_D',WRF_REAL,'0',' ',         &
  !                   '','',' ',1,1,1,1,1,1,temlg1,istatus)
!
!-----------------------------------------------------------------------
!
! Dummy arrays
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SNOWNC',WRF_REAL,'XY',' ',          &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'GRAUPELNC',WRF_REAL,'XY',' ',       &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SWDOWN',WRF_REAL,'XY',' ',          &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'GLW',WRF_REAL,'XY',' ',             &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'OLR',WRF_REAL,'XY',' ',             &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)

!-----------------------------------------------------------------------
!
!  Check whether the latitude and longitude are consistent with the data
!
!-----------------------------------------------------------------------

  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'XLAT','',                            &
                  'west_east','south_north',nx_ext,ny_ext,tem2d1,       &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                  timestr,itime,1,'XLONG','',                           &
                  'west_east','south_north',nx_ext,ny_ext,tem2d2,       &
                  nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2, istatus)

  DO j = 1,ny_ext-1
    DO i = 1,nx_ext-1
      IF( (ABS(lat_ext(i,j)-tem2d1(i,j)) > epsl .OR.                    &
           ABS(lon_ext(i,j)-tem2d2(i,j)) > epsl)                        &
           .AND. .NOT. warned ) THEN
        WRITE(6,*)'================= WARNING ======================'
        WRITE(6,'(2A,2(A,I4),A)') 'Find latitude & longitude inconsistency', &
                   ' at grid point:', ' i = ',i,' j = ',j,'.'
        WRITE(6,'(2(A,F9.3))') 'Expecting           latitude = ',      &
                      lat_ext(i,j), ' longitude = ', lon_ext(i,j)
        WRITE(6,'(2(A,F9.3))') 'Found in data file, latitude = ',      &
                      tem2d1(i,j), ' longitude = ', tem2d2(i,j)
        WRITE(6,*)'================= WARNING ======================'
        warned = .TRUE.
        istatus = -666
        RETURN
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Dummy arrays. It is not necessay if frame_per_file = 1, or io_form /= 1
!
!-----------------------------------------------------------------------

  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'XLAT_U',WRF_REAL,'XY','X',          &
  !                   'west_east_stag','south_north',' ',                &
  !                   nx_ext,ny_ext,1,nxd_stag,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'XLONG_U',WRF_REAL,'XY','X',         &
  !                   'west_east_stag','south_north',' ',                &
  !                   nx_ext,ny_ext,1,nxd_stag,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'XLAT_V',WRF_REAL,'XY','Y',          &
  !                   'west_east','south_north_stag',' ',                &
  !                   nx_ext,ny_ext,1,nxd,nyd_stag,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'XLONG_V',WRF_REAL,'XY','Y',         &
  !                   'west_east','south_north_stag',' ',                &
  !                   nx_ext,ny_ext,1,nxd,nyd_stag,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'ALBEDO',WRF_REAL,'XY',' ',          &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'TMN',WRF_REAL,'XY',' ',             &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'XLAND',WRF_REAL,'XY',' ',           &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'UST',WRF_REAL,'XY',' ',             &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
! ! CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
! !                    timestr,itime,'PBLH',WRF_REAL,'XY',' ',            &
! !                    'west_east','south_north',' ',                     &
! !                    nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'HFX',WRF_REAL,'XY',' ',             &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'QFX',WRF_REAL,'XY',' ',             &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'LH',WRF_REAL,'XY',' ',              &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
  !
  !CALL get_wrf_dummy(fHndl,io_form,multifile,ncmprx,ncmpry,             &
  !                   timestr,itime,'SNOWC',WRF_REAL,'XY',' ',           &
  !                   'west_east','south_north',' ',                     &
  !                   nx_ext,ny_ext,1,nxd,nyd,1,temlg1,istatus)
!
!-----------------------------------------------------------------------
!
!  Restore the original map projection before return
!
!-----------------------------------------------------------------------

  CALL setmapr(iproj_orig,scale_orig,latnot_orig,trlon_orig)
  CALL setorig(1,x0_orig,y0_orig)

!-----------------------------------------------------------------------
!
! Deallocate the working arrays
!
!-----------------------------------------------------------------------

  DEALLOCATE(temlg1,temlg2)

  istatus = 0
  RETURN
END SUBROUTINE getwrfd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE adj_wrfuv                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE adj_wrfuv(multifile,use_wrf_grid,nx_ext,ny_ext,nz_ext,       &
                iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext, &
                lonu_ext,lonv_ext,u_ext,vatu_ext,uatv_ext,v_ext,        &
                tem1_ext,tem2_ext,istatus)

!------------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extended WRF horizontal wind arrays for fake zone.
!  Rotate WRF horizontal wind to be earth relative if needed.
!  This is seperated from getwrfd for ordering I/O purpose in MPI mode.
!
!------------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (04/22/2005)
!
!  MODIFICATION HISTORY:
!
!------------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL,      INTENT(IN)  :: multifile
  INTEGER,      INTENT(IN)  :: use_wrf_grid ! 0 - rotate hori. wind
                                            ! 1 - no rotation
  INTEGER,      INTENT(IN)  :: nx_ext       !
  INTEGER,      INTENT(IN)  :: ny_ext       !
  INTEGER,      INTENT(IN)  :: nz_ext       !

  INTEGER, INTENT(IN)  :: iproj_ext    ! Map projection option
                                       ! = 1, polar projection;
                                       ! = 2, Lambert projection;
                                       ! = 3, Mercator projection.
  REAL,    INTENT(IN)  :: scale_ext    ! Map scale factor (should be 1.0)
  REAL,    INTENT(IN)  :: trlon_ext    ! True longitude
  REAL,    INTENT(IN)  :: latnot_ext(2)! True latitude
  REAL,    INTENT(IN)  :: x0_ext
  REAL,    INTENT(IN)  :: y0_ext

  REAL,    INTENT(IN)  :: lonu_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: lonv_ext(nx_ext,ny_ext)

  REAL,    INTENT(INOUT) :: u_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(INOUT) :: v_ext   (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT)   :: vatu_ext(nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT)   :: uatv_ext(nx_ext,ny_ext,nz_ext)

  REAL,    INTENT(OUT) :: tem1_ext (nx_ext,ny_ext,nz_ext)
  REAL,    INTENT(OUT) :: tem2_ext (nx_ext,ny_ext,nz_ext)
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------------
!
!  Include files
!
!------------------------------------------------------------------

  INCLUDE 'mp.inc'

!------------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER  :: iproj_orig
  REAL     :: scale_orig
  REAL     :: latnot_orig(2)
  REAL     :: trlon_orig
  REAL     :: x0_orig, y0_orig
  REAL     :: xsub0,ysub0

  REAL, ALLOCATABLE :: uext(:,:,:)
  REAL, ALLOCATABLE :: vext(:,:,:)

  INTEGER            :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (multifile) THEN  ! extend values to fake zone
    CALL mpext_wrf_u(u_ext,nx_ext,ny_ext,nz_ext,tem1_ext)
    CALL mpext_wrf_v(v_ext,nx_ext,ny_ext,nz_ext,tem1_ext)
  END IF

  IF (use_wrf_grid == 1) THEN   ! no rotation,
                                ! uatv_ext & vatu_ext will not be used
    uatv_ext(:,:,:) = -9999.0
    vatu_ext(:,:,:) = -9999.0

  ELSE                          ! == 0, rotate hori. wind

    !
    ! remember the map projection paramters which should be restored
    ! before return
    !

    CALL getmapr(iproj_orig,scale_orig,latnot_orig,                     &
                 trlon_orig,x0_orig,y0_orig)

    ! Set up WRF map projection
    CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
    CALL setorig(1,x0_ext,y0_ext)

    ALLOCATE(uext( 1:nx_ext+1, 0:ny_ext,  1:nz_ext), STAT = istatus)
    ALLOCATE(vext( 0:nx_ext,   1:ny_ext+1,1:nz_ext), STAT = istatus)

    CALL extend_u(u_ext,nx_ext,ny_ext,nz_ext,uext,tem1_ext,tem2_ext)
    CALL extend_v(v_ext,nx_ext,ny_ext,nz_ext,vext,tem1_ext,tem2_ext)

    !
    ! get u at V grid point locations
    !
    DO k=1,nz_ext
      DO j=1,ny_ext
       DO i=1,nx_ext
          uatv_ext(i,j,k) = 0.25*(  uext(i,j-1,k) + uext(i+1,j-1,k)     &
                                  + uext(i,j,k)   + uext(i+1,j,k))
        END DO
      END DO
    END DO

    !
    ! get V at U grid point locations
    !

    DO k = 1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          vatu_ext(i,j,k) = 0.25*(  vext(i-1,j,k)   + vext(i,j,k)       &
                                  + vext(i-1,j+1,k) + vext(i,j+1,k))
        END DO
      END DO
    END DO

    !
    !  Orient u & v to true north.
    !

    DO k = 1, nz_ext
      CALL uvmptoe(nx_ext,ny_ext,u_ext(:,:,k),vatu_ext(:,:,k),lonu_ext, &
                   tem1_ext(:,:,k),tem2_ext(:,:,k))
      u_ext   (1:nx_ext,1:ny_ext,k) = tem1_ext(1:nx_ext,1:ny_ext,k)
      vatu_ext(1:nx_ext,1:ny_ext,k) = tem2_ext(1:nx_ext,1:ny_ext,k)

      CALL uvmptoe(nx_ext,ny_ext,uatv_ext(:,:,k),v_ext(:,:,k),lonv_ext, &
                   tem1_ext(:,:,k),tem2_ext(:,:,k))
      uatv_ext(1:nx_ext,1:ny_ext,k) = tem1_ext(1:nx_ext,1:ny_ext,k)
      v_ext   (1:nx_ext,1:ny_ext,k) = tem2_ext(1:nx_ext,1:ny_ext,k)

    END DO

!-----------------------------------------------------------------------
!
! Deallocate the working arrays
!
!-----------------------------------------------------------------------

    DEALLOCATE(uext,vext)
!
!-----------------------------------------------------------------------
!
!  Restore the original map projection before return
!
!-----------------------------------------------------------------------

    CALL setmapr(iproj_orig,scale_orig,latnot_orig,trlon_orig)
    CALL setorig(1,x0_orig,y0_orig)

  END IF   ! use_wrf_grid

  istatus = 0
  RETURN
END SUBROUTINE adj_wrfuv
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
SUBROUTINE open_wrf_file(filename,io_form,multifile,for_meta_only,      &
                         ncmprx,ncmpry,numdigits,nidout)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Open a WRF file and return NetCDF file handler.
!    It will call open_wrf_one_file or open_wrf_multi_files depends
!    on the pass-in parameters
!
!    NOTE: it is required to call close_wrf_file explicitly to close
!          the opened file in your calling program.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: multifile
  LOGICAL,          INTENT(IN)  :: for_meta_only
  INTEGER,          INTENT(IN)  :: ncmprx, ncmpry,numdigits
  INTEGER,          INTENT(OUT) :: nidout(ncmprx,ncmpry)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  INTEGER            :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (multifile) THEN
    CALL open_wrf_multi_files(filename,io_form,for_meta_only,ncmprx,    &
                             ncmpry,numdigits,nidout,istatus)
  ELSE
    CALL open_wrf_one_file(filename,io_form,nidout,istatus)
  END IF

  IF (istatus /= 0) THEN
    WRITE(0,'(1x,2a)') 'ERROR: Opening file ',filename
    CALL arpsstop('Open WRF file error.',1)
  END IF

  RETURN
END SUBROUTINE open_wrf_file
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
SUBROUTINE close_wrf_file(nch,io_form,multifile,for_meta_only,          &
                          ncompressx,ncompressy)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!     Close the WRF file which is opened using open_wrf_file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: io_form
  LOGICAL, INTENT(IN) :: multifile
  LOGICAL, INTENT(IN) :: for_meta_only
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: nch(ncompressx,ncompressy)

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
!
  INTEGER :: istatus

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

  IF(multifile) THEN
    CALL close_wrf_multi_files(nch,io_form,for_meta_only,               &
                               ncompressx,ncompressy,istatus)
  ELSE
    CALL close_wrf_one_file(nch,io_form,istatus)
  END IF

  IF (istatus /= 0) THEN
    WRITE(0,'(1x,2a)') 'ERROR: closing file handler ',nch
    CALL mpexit(1)
  END IF

  RETURN
END SUBROUTINE close_wrf_file

SUBROUTINE io_shutdown(io_form)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: io_form

  INTEGER :: istatus

  istatus = 0

  IF (io_form == 5) THEN
    CALL shutdown_phdf5_io(istatus)
  ELSE IF (io_form == 1) THEN
    CALL ext_int_ioexit( iStatus )
  END IF

  RETURN
END SUBROUTINE io_shutdown
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
SUBROUTINE get_wrf_Times(nfid,io_form,multifile,ncompressx,ncompressy,  &
                         itime,timestr)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read the the Date String in the WRF outputs at specified time
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ncompressx, ncompressy
  INTEGER, INTENT(IN)  :: nfid(ncompressx,ncompressy)     ! file handler
  INTEGER, INTENT(IN)  :: io_form  ! File format
  LOGICAL, INTENT(IN)  :: multifile
  INTEGER, INTENT(IN)  :: itime    ! Time dimension value
                                   ! this is the unlimited dimension
  CHARACTER(LEN=*), INTENT(OUT) :: timestr

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
  INTEGER :: istatus

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (multifile) THEN
     CAlL get_wrf_Times_from_multi_files(nfid,io_form,ncompressx,       &
                                 ncompressy,itime,timestr,istatus)
  ELSE
     CALL get_wrf_Times_from_one_file(nfid,io_form,itime,timestr,istatus)
  END IF

!  write(0,*) 'Next_time = ', timestr
  RETURN
END SUBROUTINE get_wrf_Times
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
SUBROUTINE get_wrf_metadata(nid,io_form,multifile,for_meta_only,        &
                          ncompressx,ncompressy,                        &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          mapproj,trlat1,trlat2,trlon,ctrlat,ctrlon,    &
                          dx,dy,dt,sfcphys,sfclay,mpphys,num_scalar,    &
                          istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve WRF grib information from the NetCDF file which are stored
!   as Global attributes.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncompressx,ncompressy
  INTEGER, INTENT(IN)  :: nid(ncompressx,ncompressy)
  INTEGER, INTENT(IN)  :: io_form
  LOGICAL, INTENT(IN)  :: multifile
  LOGICAL, INTENT(IN)  :: for_meta_only
  INTEGER, INTENT(OUT) :: nx_ext, ny_ext     ! they are whole domain dimensions
  INTEGER, INTENT(OUT) :: nz_ext, nzsoil_ext
  INTEGER, INTENT(OUT) :: mapproj
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
  INTEGER, INTENT(OUT) :: num_scalar  ! number of WRF scalars that mapped to ARPS variables
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variable
!
!-----------------------------------------------------------------------

  INTEGER :: iloc, jloc
  INTEGER :: iproj

  INCLUDE 'globcst.inc'

  INTEGER :: P_QNN, P_QNC, P_QNR, P_QNI, P_QNS, P_QNG, P_QNH
  INTEGER :: P_QT, P_RIM, P_REF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (multifile) THEN

    IF (for_meta_only .OR. io_form == 7) THEN    ! read only one file is enough
      CALL  get_wrf_meta_from_multi_files(nid(1,1),io_form,.TRUE.,      &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj,trlat1,trlat2,trlon,ctrlat,ctrlon,      &
                          dx,dy,dt,sfcphys,sfclay,mpphys,istatus)

    ELSE        ! binary file requires sequential access

      DO jloc = 1,ncompressy
        DO iloc = 1,ncompressx
          CALL  get_wrf_meta_from_multi_files(nid(iloc,jloc),           &
                          io_form,.FALSE.,       &  ! do not do checking
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj,trlat1,trlat2,trlon,ctrlat,ctrlon,      &
                          dx,dy,dt,sfcphys,sfclay,mpphys,istatus)
        END DO
      END DO

    END IF

  ELSE
    CALL get_wrf_meta_from_one_file(nid,io_form,                        &
                          nx_ext,ny_ext,nz_ext,nzsoil_ext,              &
                          iproj,trlat1,trlat2,trlon,ctrlat,ctrlon,      &
                          dx,dy,dt,sfcphys,sfclay,mpphys,istatus)
  END IF

!-----------------------------------------------------------------------
!
!  Convert from WRF map projection parameter to ARPS map projection
!
!-----------------------------------------------------------------------

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
    istatus = -555
    CALL arpsstop('WRONG WRF map projection parameter.',1)
  END IF

  IF(trlat1 < 0.0) mapproj = -1*mapproj

!-----------------------------------------------------------------------
!
! Set up moist and scalar variables to be read based on mp_physics
!
!-----------------------------------------------------------------------

  P_QC = 0
  P_QR = 0
  P_QI = 0
  P_QS = 0
  P_QG = 0
  P_QH = 0
  nscalar  = 0
  nscalarq = 0
  qnames(:)= ' '; qdescp(:)= ' '

  P_QT  = 0
  P_QNN = 0
  P_QNC = 0
  P_QNR = 0
  P_QNI = 0
  P_QNS = 0
  P_QNG = 0
  P_QNH = 0
  P_RIM = 0
  P_REF = 0
  num_scalar = 0

  SELECT CASE (mpphys)
  CASE (0)                           ! passiveqv
    nscalar = 0
  CASE (1,3)                         ! kesslerscheme, wsm3scheme
    P_QC      = 1
    P_QR      = 2
    nscalar   = 2
    nscalarq  = 2
  CASE (2,6)                         ! linscheme, wsm6scheme
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    P_QG      = 5
    nscalar    = 5
    nscalarq   = 5

    P_REF = 1

  CASE (4,99)                        ! wsm5scheme,ncepcloud5
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    nscalar    = 4
    nscalarq   = 4

    P_REF = 1

  CASE (5)                           ! etampnew
    P_QC      = 1
    P_QR      = 2
    P_QS      = 3
    nscalar    = 3
    nscalarq   = 3

    P_QT        = 1
    P_RIM       = 2                  ! scpecial for this scheme

    P_REF = 3

  CASE (8)                           ! thompson
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    P_QG      = 5
    nscalar   = 5
    nscalarq  = 5

    P_QNR      = 1
    P_QNI      = 2

    P_REF = 3

  CASE (9,17)                        ! M-Y scheme, NSSL
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    P_QG      = 5
    P_QH      = 6
    nscalar   = 6
    nscalarq  = 6

    P_QNC      = 1
    P_QNR      = 2
    P_QNI      = 3
    P_QNS      = 4
    P_QNG      = 5
    P_QNH      = 6

    P_REF = 7

  CASE (10)                            ! Morrison DB
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    P_QG = 5
    nscalar  = 5
    nscalarq = 5

    P_QNR = 1
    P_QNI = 2
    P_QNS = 3
    P_QNG = 4

    P_REF = 5

  CASE (14)                            ! DM-5
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    nscalar  = 4
    nscalarq = 4

    P_QNC = 1
    P_QNR = 2
    P_QNN = 3               ! Do not know how to handle it in the ARPS grid?

    P_REF = 4

  CASE (16)                            ! DM-6
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    P_QG = 5
    nscalar  = 5
    nscalarq = 5

    P_QNC = 1
    P_QNR = 2
    P_QNN = 3

    P_REF = 4

  CASE (98)                            ! Old Thompson (07)
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    P_QG = 5
    nscalar  = 5
    nscalarq = 5

    P_QNI = 1

  CASE DEFAULT
    istatus = -1
    WRITE(6,'(/,1x,a,I2,a,/)') 'ERROR: Wrong parameter - mp_physics = ',&
                                mpphys,'.'
  END SELECT

!-----------------------------------------------------------------------
!
! Mapping WRF qscalar to ARPS variables
!
!-----------------------------------------------------------------------

  IF (P_QC > 0) THEN
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
  END IF
  IF (P_QR > 0) THEN
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
  END IF
  IF (P_QI > 0) THEN
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
  END IF
  IF (P_QS > 0) THEN
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
  END IF
  IF (P_QG > 0) THEN
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
  END IF
  IF (P_QH > 0) THEN
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
  END IF

  ! Number of concentrations
  IF (P_QNC > 0) THEN
    P_NC = nscalarq + P_QNC; nscalar = nscalar + 1
    qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water concentrations (#/m3)'
  END IF
  IF (P_QNR > 0) THEN
    P_NR = nscalarq + P_QNR; nscalar = nscalar + 1
    qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Rain water concentrations (#/m3)'
  END IF
  IF (P_QNI > 0) THEN
    P_NI = nscalarq + P_QNI; nscalar = nscalar + 1
    qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Cloud ice concentrations (#/m3)'
  END IF
  IF (P_QNS > 0) THEN
    P_NS = nscalarq + P_QNS; nscalar = nscalar + 1
    qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow concentrations (#/m3)'
  END IF
  IF (P_QNG > 0) THEN
    P_NG = nscalarq + P_QNG; nscalar = nscalar + 1
    qnames(P_NG) = 'ng'; qdescp(P_NG) = 'Graupel concentrations (#/m3)'
  END IF
  IF (P_QNH > 0) THEN
    P_NH = nscalarq + P_QNH; nscalar = nscalar + 1
    qnames(P_NH) = 'nh'; qdescp(P_NH) = 'Hail concentrations (#/m3)'
  END IF

!-----------------------------------------------------------------------
!
! WRF extra variables
!
!-----------------------------------------------------------------------

  num_scalar = nscalar

  IF (P_RIM > 0) THEN         ! only mphysics = 5 contain this variable
    nscalar = nscalar + 1
    qnames(nscalar) = 'F_RIMEF_PHY'; qdescp(nscalar) = 'F_RIMEF_PHY'
  END IF

  IF (P_REF > 0) THEN         ! only mphysics = 9 contain this variable
                              ! all schemes have this variable since 2012
    nscalar = nscalar + 1
    qnames(nscalar) = 'REFL_10CM'; qdescp(nscalar) = 'REFL_10CM'
  END IF

  RETURN
END SUBROUTINE get_wrf_metadata

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
SUBROUTINE get_wrf_dummy(nid,io_form,multifile,ncompressx,ncompressy,   &
                datestr,itime,varname,varType,memoryorder,stagger,      &
                dimname1,dimname2,dimname3,                             &
                nx,ny,nz,nxd,nyd,nzd,temtd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in an array from the WRF history file. It just for sequential
!    access of WRF binary file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncompressx,ncompressy
  INTEGER,          INTENT(IN)  :: nid(ncompressx,ncompressy)
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: multifile
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( io_form /= 1 ) RETURN   ! Only needed for binary format

  IF ( multifile ) THEN
    CALL get_wrf_dummy_from_multi_files(nid,io_form,                    &
                ncompressx,ncompressy,datestr,itime,varname,            &
                varType,memoryorder,stagger,dimname1,dimname2,dimname3, &
                nx,ny,nz,nxd,nyd,nzd,temtd,istatus)
  ELSE
    CALL get_wrf_dummy_from_one_file(nid,io_form,datestr,itime,varname, &
                varType,memoryorder,stagger,dimname1,dimname2,dimname3, &
                nx,ny,nz,nxd,nyd,nzd,temtd,istatus)
  END IF

  RETURN
END SUBROUTINE get_wrf_dummy
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
SUBROUTINE get_wrf_1d(nfid,io_form,multifile,ncompressx,ncompressy,     &
                      datestr,itime,varname,stagger,                    &
                      dimname1,nz,var1d,nzd, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 1D array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncompressx,ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: multifile
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( multifile ) THEN
    CALL get_wrf_1d_from_multi_files(nfid,io_form,ncompressx,ncompressy,&
                      datestr,itime,varname,stagger,                    &
                      dimname1,nz,var1d,nzd,istatus)
  ELSE
    CALL  get_wrf_1d_from_one_file  (nfid,io_form,                      &
                      datestr,itime,varname,stagger,                    &
                      dimname1,nz,var1d,nzd,istatus)
  END IF

  RETURN
END SUBROUTINE get_wrf_1d
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
SUBROUTINE get_wrf_2d(nfid,io_form,multifile,ncompressx,ncompressy,     &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2d,                    &
                      nxd,nyd,temtd,nxlg,nylg,temlg, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 2D array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncompressx,ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: multifile
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
  REAL,             INTENT(OUT) :: temtd(nxd,nyd)  ! data array
  INTEGER,          INTENT(IN)  :: nxlg,nylg       ! memory index for the whole domain
  REAL,             INTENT(OUT) :: temlg(nxlg,nylg)! memory array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( multifile ) THEN

    CALL get_wrf_2d_from_multi_files(nfid,io_form,ncompressx,ncompressy,&
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2d,                    &
                      nxd,nyd,temtd,istatus)
  ELSE
    CALL get_wrf_2d_from_one_file(nfid,io_form,                         &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2d,                    &
                      nxd,nyd,temtd,nxlg,nylg,temlg, istatus)
  END IF

  RETURN
END SUBROUTINE get_wrf_2d
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
SUBROUTINE get_wrf_2di(nfid,io_form,multifile,ncmprx,ncmpry,           &
                      datestr,itime,fzone,varname,stagger,             &
                      dimname1,dimname2,nx,ny,var2di,nxd,nyd,temtd,    &
                      nxlg,nylg,temlg,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 2D integer array from the WRF history file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncmprx,ncmpry
  INTEGER,          INTENT(IN)  :: nfid(ncmprx,ncmpry)
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: multifile
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
  INTEGER,          INTENT(OUT) :: temtd(nxd,nyd)   ! domain array
  INTEGER,          INTENT(IN)  :: nxlg,nylg          ! memory
  INTEGER,          INTENT(OUT) :: temlg(nxlg,nylg)   ! memory array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( multifile ) THEN
    CALL get_wrf_2di_from_multi_files(nfid,io_form,ncmprx,ncmpry,       &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2di,                   &
                      nxd,nyd,temtd,istatus)
  ELSE
    CALL get_wrf_2di_from_one_file(nfid,io_form,                        &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,nx,ny,var2di,                   &
                      nxd,nyd,temtd,nxlg,nylg,temlg,istatus)
  END IF

  RETURN
END SUBROUTINE get_wrf_2di
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
SUBROUTINE get_wrf_3d(nfid,io_form,multifile,ncompressx,ncompressy,     &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,dimname3,nx,ny,nz,var3d,        &
                      nxd,nyd,nzd,temtd,nxlg,nylg,nzlg,temlg,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncompressx,ncompressy
  INTEGER,          INTENT(IN)  :: nfid(ncompressx,ncompressy)
  INTEGER,          INTENT(IN)  :: io_form
  LOGICAL,          INTENT(IN)  :: multifile
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
  REAL,             INTENT(OUT) :: temtd(nxd*nyd*nzd)      ! domain array
  INTEGER,          INTENT(IN)  :: nxlg,nylg,nzlg          ! domain index
  REAL,             INTENT(OUT) :: temlg(nxlg,nylg,nzlg)   ! memory array
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( multifile ) THEN
    CALL get_wrf_3d_from_multi_files(nfid,io_form,ncompressx,ncompressy,&
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,dimname3,nx,ny,nz,var3d,        &
                      nxd,nyd,nzd,temtd,istatus)
  ELSE
    CALL get_wrf_3d_from_one_file(nfid,io_form,                         &
                      datestr,itime,fzone,varname,stagger,              &
                      dimname1,dimname2,dimname3,nx,ny,nz,var3d,        &
                      nxd,nyd,nzd,temtd,nxlg,nylg,nzlg,temlg,istatus)
  END IF

  RETURN
END SUBROUTINE get_wrf_3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE check_wrf_files              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE check_wrf_files(multifile,MAXWRFFIL,grid_id,io_form,         &
           nprocx_in,ncompressx,ncompressy,numdigits,                   &
           abstimes,abstimei,abstimee,rewindyr,myr,                     &
           dir_extd,extdname,nextdfil,frames_per_outfile,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE: Check the existence of WRF files to be read and return the
!          valid file number, file names for each processor, number of
!          frames in each WRF file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (03/26/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL, INTENT(IN)    :: multifile
  INTEGER, INTENT(IN)    :: MAXWRFFIL
  INTEGER, INTENT(IN)    :: grid_id
  INTEGER, INTENT(IN)    :: io_form
  INTEGER, INTENT(IN)    :: nprocx_in
  INTEGER, INTENT(IN)    :: ncompressx, ncompressy
  INTEGER, INTENT(IN)    :: numdigits
  INTEGER, INTENT(IN)    :: abstimes, abstimei, abstimee
  LOGICAL, INTENT(IN)    :: rewindyr
  INTEGER, INTENT(IN)    :: myr
  CHARACTER(LEN=256), INTENT(IN)  :: dir_extd
  CHARACTER(LEN=256), INTENT(OUT) :: extdname(MAXWRFFIL)
  INTEGER,            INTENT(OUT) :: nextdfil
  INTEGER,            INTENT(OUT) :: frames_per_outfile(MAXWRFFIL)
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: iloc, jloc
  INTEGER :: iloc_x, jloc_y, loc_proc

  INTEGER :: year, month, day, hour, minute, second
  INTEGER :: ifile

  CHARACTER(LEN=256) :: tmpstr
  CHARACTER(LEN=40)  :: tmpwrftime
  CHARACTER(LEN=20)  :: fmtstr

  LOGICAL :: fexist
  LOGICAL :: First_in_process

  INTEGER :: itmp

!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nextdfil    = 0
  extdname(:) = ' '
  istatus     = 0

  IF (myproc == 0) WRITE(6,'(1x,a,/,1x,a)')                             &
              '============================','WRF files to be read are:'

  IF (multifile) THEN

    iloc_x = (loc_x-1)*ncompressx    ! column of processors
    jloc_y = (loc_y-1)*ncompressy    ! rows of processors

    WRITE(fmtstr,'(a,I1,a,I1,a)') '(2a,I',numdigits,'.',numdigits,')'
    ifile = abstimes
    nextdfil_loop1: DO WHILE(ifile <= abstimee)

      CALL abss2ctim(ifile,year,month,day,hour,minute,second)
      IF (rewindyr) year = myr

      nextdfil = nextdfil + 1
      WRITE(tmpwrftime,'(a,I2.2,a,I4.4,5(a,I2.2))')                     &
             'wrfout_d',grid_id,'_',                                    &
             year,'-',month,'-',day,'_',hour,':',minute,':',second

      IF (dir_extd == '<dir>') THEN
        WRITE(extdname(nextdfil),'(3a)') TRIM(tmpwrftime),'/',TRIM(tmpwrftime)
      ELSE
        WRITE(extdname(nextdfil),'(2a)') TRIM(dir_extd),TRIM(tmpwrftime)
      END IF

      First_in_process = .TRUE.
      DO jloc = 1,ncompressy
        DO iloc = 1, ncompressx

          loc_proc = (jloc_y+jloc-1)*nprocx_in + iloc_x+(iloc-1)
          WRITE(tmpstr,FMT=fmtstr) TRIM(extdname(nextdfil)),'_',loc_proc

          INQUIRE(FILE=TRIM(tmpstr), EXIST = fexist )
          IF(.NOT. fexist) THEN
            WRITE(6,'(1x,3a)') 'WARNING: The WRF file "',TRIM(tmpstr),  &
                           '" does not exist.'
            nextdfil = nextdfil - 1
            ifile = ifile + abstimei   ! try next time level
            CYCLE nextdfil_loop1
          ELSE
            CALL get_wrf_frames_per_outfile(TRIM(tmpstr),io_form,itmp,istatus)
            IF ( First_in_process) THEN
              frames_per_outfile(nextdfil) = itmp
              First_in_process = .FALSE.
              IF(myproc == 0) WRITE(6,'(5x,a,I2.2,3a,I4,a)')            &
                 'WRF file ',nextdfil,': ', TRIM(tmpstr),               &
                 ', frames = ',frames_per_outfile(nextdfil),'.'
            ELSE
              ! Ensure all files read by this processor have same frames
              IF (itmp /= frames_per_outfile(nextdfil) ) THEN
                WRITE(6,'(1x,a,I4,2(a,I2),a,/,10x,a,I4,a,I2,a,I4,a,/)') &
                  'WARNING: unlimited dimension in processor - ',myproc,&
                  ' for file (',iloc,', ',jloc,                         &
                  ') is different with others files.',                  &
                  'We got iframes = ',itmp,                             &
                  ', but other files have frames_per_outfile(',         &
                  nextdfil,') = ',frames_per_outfile(nextdfil),'.'
                istatus = -1
                RETURN
              END IF
            END IF

          END IF
        END DO
      END DO

      !
      ! Ensure all processors have the same frames
      !
      itmp = frames_per_outfile(nextdfil)
      CALL mpmaxi(itmp)
      IF ( itmp /= frames_per_outfile(nextdfil) ) THEN
        WRITE(6,'(1x,a,/,8x,3(a,I4),/,8x,a,I4,a,/)')                    &
          'ERROR: Each processor has got different numbers for frames in WRF files.', &
          'Processor: ',myproc,' found ',frames_per_outfile(nextdfil),  &
          ' frames in WRF file - ',nextdfil,                            &
          'While other processor may found ',itmp,' frames in that file.'
        istatus = -2
        RETURN
      END IF

      ifile = ifile + frames_per_outfile(nextdfil)*abstimei
    END DO nextdfil_loop1

  ELSE     ! non-mpi or mpi with one file

    IF (myproc == 0) THEN

      ifile = abstimes
      nextdfil_loop2: DO WHILE (ifile <= abstimee)

        CALL abss2ctim(ifile,year,month,day,hour,minute,second)
        IF (rewindyr) year = myr

        nextdfil = nextdfil + 1
        WRITE(tmpwrftime,'(a,I2.2,a,I4.4,5(a,I2.2))')                   &
               'wrfout_d',grid_id,'_',                                  &
               year,'-',month,'-',day,'_',hour,':',minute,':',second

        IF (dir_extd == '<dir>') THEN
          WRITE(extdname(nextdfil),'(3a)') TRIM(tmpwrftime),'/',TRIM(tmpwrftime)
        ELSE
          WRITE(extdname(nextdfil),'(2a)') TRIM(dir_extd),TRIM(tmpwrftime)
        END IF


        INQUIRE(FILE=TRIM(extdname(nextdfil)), EXIST = fexist )
        IF(.NOT. fexist) THEN
          WRITE(6,'(1x,3a)') 'WARNING: The WRF file "',                 &
                         TRIM(extdname(nextdfil)),'" does not exist.'
          ifile = ifile + abstimei
          nextdfil = nextdfil - 1
          CYCLE nextdfil_loop2
        ELSE
          CALL get_wrf_frames_per_outfile(TRIM(extdname(nextdfil)),     &
                                          io_form,itmp,istatus)
          frames_per_outfile(nextdfil) = itmp
          WRITE(6,'(5x,a,I2.2,3a,I4,a)')                                &
             'WRF file ',nextdfil,': ', TRIM(extdname(nextdfil)),       &
             ', frames = ',frames_per_outfile(nextdfil),'.'
        END IF

        ifile = ifile + frames_per_outfile(nextdfil)*abstimei
      END DO nextdfil_loop2

    END IF      ! myproc == 0

    CALL mpupdatei(nextdfil,1)
    CALL mpupdatei(frames_per_outfile,MAXWRFFIL)
    DO ifile = 1, nextdfil
      CALL mpupdatec(extdname(ifile),256)
    END DO

  END IF     ! multifile

!-----------------------------------------------------------------------
!
! Validate nextdfil before return
!
!-----------------------------------------------------------------------

  IF(nextdfil < 1) THEN
    WRITE(6,'(a)') 'No input WRF file was valid. Please check the input file.'
    istatus = -3
    RETURN
  END IF

  itmp = nextdfil
  CALL mpmaxi(nextdfil)
  IF (itmp /= nextdfil) THEN
    WRITE(6,'(a)') 'Each processor may process different number of WRF files.'
    WRITE(6,'(2(a,I4),a,/,a,I4,a)') 'Processor: ',myproc,' found ',   &
                    itmp,' WRF files,',                               &
                    'While other processor may found ',nextdfil,' files.'
    istatus = -4
    RETURN
  END IF

  RETURN
END SUBROUTINE check_wrf_files
!
!##################################################################
!##################################################################
!######                                                      ######
!######       SUBROUTINE get_wrf_frames_per_outfile          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_frames_per_outfile(filename,io_form,numframes,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Get the size of unlimited dimension in the WRF data file
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (03/26/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(OUT) :: numframes
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
  IF (io_form == 7) THEN
    CALL get_ncd_frames_per_outfile(TRIM(filename),numframes,istatus)
  ELSE
    istatus   = 1
    numframes = 1
    WRITE(6,'(1x,a,/,10x,a,/)')       &
      'WARNING: Multiple frames in one file is only support on netCDF file at present', &
      'frames_per_outfile is reset to 1. Pleace check for correctness.'
  END IF

  RETURN
END SUBROUTINE get_wrf_frames_per_outfile
