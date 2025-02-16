PROGRAM wrfstatic
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM WRFSTATIC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!
!  PURPOSE:
!
!    Read static file sets and produce static file for ARPS2WRF & WRF.
!
!  NOTE:
!    o It will replace both grid_gen.exe and staticpost in WRFSI.
!    o This program shares the same namelist input file as ARPS2WRF
!      for convenience with an extra namelist variable for static
!      dataset directory in input/arps2wrf.input.
!
!  Required Static datasets (the same as WRFSI requirements):
!
!    1. albedo_ncep
!    2. greenfrac
!    3. islope
!    4. landuse_30s
!    5. maxsnowalb
!    6. soiltemp_1deg
!    7. soiltype_bot_30s
!    8. soiltype_top_30s
!    *. topo_30s (not required for ARPS2WRF, because terrain will
!       be provided by arpstrn/arpstern.)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (09/28/2005)
!  Created initially based on WRFSIV2.1.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

  USE wrf_metadata             ! WRF constants and metadata

!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!----------------------------------------------------------------------
!
! ARPS grid variables
!
!---------------------------------------------------------------------

  INTEGER, PARAMETER :: nmax_domains = 100
  INTEGER, PARAMETER :: nhisfile_max = 100
  INTEGER, PARAMETER :: max_vertical_levels = 100

  CHARACTER(LEN=256) :: grdbasfn
  CHARACTER(LEN=256) :: bdybasfn(nhisfile_max)
  CHARACTER(LEN=256) :: hisfile(nhisfile_max)
  INTEGER            :: nhisfile,lengbf

  CHARACTER(LEN=256), DIMENSION(nmax_domains) :: adashisfn
  CHARACTER(LEN=256), DIMENSION(nmax_domains) :: adasbasfn
  INTEGER,            DIMENSION(nmax_domains) :: hinfmt

  LOGICAL,            DIMENSION(nmax_domains) :: input_from_file
  INTEGER            :: max_dom

  INTEGER            :: finexist(nhisfile_max)
  CHARACTER(LEN=80)  :: execname
!
!-----------------------------------------------------------------------
!
!  ARPS include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
! Variables for split files
!
!-----------------------------------------------------------------------

  INTEGER :: nprocx_in, nprocy_in

!-----------------------------------------------------------------------
!
!  ARPS variables to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER, DIMENSION(nmax_domains) :: nxnm,nynm,nznm
  INTEGER, DIMENSION(nmax_domains) :: nzsoilnm
  INTEGER, DIMENSION(nmax_domains) :: nstypsnm

  INTEGER :: nx,ny,nz          ! Grid dimensions for ARPS.
  INTEGER :: nzsoil            ! Soil levels
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: mapproj
  REAL    :: sclfct
  REAL    :: trulat1, trulat2, trulon
  REAL    :: ctrlat,  ctrlon
  REAL    :: dx, dy

  REAL    :: time

!-----------------------------------------------------------------------
!
!  Namelist definitions for ARPS2WRF.input
!
!     sfcdt              Specify surface characteristics
!     bdyspc             Obtain boundary input files (ARPS format)
!     wrf_grid           Define WRF horizontal and vertical grid
!     interp_options     Choose interpolation scheme
!     wrf_opts           WRF options from namelist.input
!     output             Ouput options
!
!-----------------------------------------------------------------------
!
  INTEGER, DIMENSION(nmax_domains) :: nx_wrfnm
  INTEGER, DIMENSION(nmax_domains) :: ny_wrfnm
  INTEGER, DIMENSION(nmax_domains) :: i_parent_start
  INTEGER, DIMENSION(nmax_domains) :: j_parent_start
  INTEGER, DIMENSION(nmax_domains) :: parent_id
  INTEGER                          :: i_parent_end   ! they are use as temporary variable
  INTEGER                          :: j_parent_end   ! So do not have to be saved
  INTEGER                          :: parent_grid_ratio

  INTEGER           :: nx_wrf         ! = nx-2 if the same grid as ARPS
  INTEGER           :: ny_wrf         ! = ny-2 if the same grid as ARPS
  INTEGER           :: nz_wrf
  INTEGER           :: nzsoil_wrf

  REAL              :: dx_wrfscl      ! grid length times map scale
  REAL              :: dy_wrfscl
  REAL              :: lattru_wrf(2)  ! array of true latitude of WRF map projection
  REAL              :: lontru_wrf     ! true longitude of WRF map projection
                                      ! = trulon_wrf
  REAL              :: sclf_wrf       ! map scale factor (usually = 1.0)

! Namelist variable declaration

  CHARACTER(LEN=5),  DIMENSION(nmax_domains)  :: sfcinitopt      ! either "ARPS" or "WRFSI"
  CHARACTER(LEN=256)                          :: geogdir
  CHARACTER(LEN=19), DIMENSION(nmax_domains)  :: start_date
  CHARACTER(LEN=256),DIMENSION(nmax_domains)  :: sfcdtfn

  REAL              :: wvln, silwt

  INTEGER, DIMENSION(nmax_domains)  :: use_arps_grid   ! Use ARPS horizontal grid as WRF grid

  INTEGER :: mapproj_wrf     ! Type of map projection in WRF model grid
                             ! modproj = 1  Polar Stereographic
                             !              projection.
                             !         = 2  Mercator projection.
                             !         = 3  Lambert projection.

  REAL    :: sclfct_wrf      ! Map scale factor.
                             ! Distance on map, between two latitudes
                             ! trulat1 and trulat2,
                             ! is = (Distance on earth)*sclfct.
                             ! For ARPS model runs,
                             ! generally this is 1.0

  REAL    :: trulat1_wrf, trulat2_wrf, trulon_wrf
                             ! 1st, 2nd real true latitude and true longitude
                             ! of WRF map projection
  REAL    :: ctrlat_wrf      ! Center latitude of WRF model domain (deg. N)
  REAL    :: ctrlon_wrf      ! Center longitude of WRF model domain (deg. E)

  REAL    :: dx_wrf          ! WRF Grid spacing in x-direction
  REAL    :: dy_wrf          ! WRF Grid spacing in y-direction

  REAL, DIMENSION(nmax_domains)  :: dx_wrfnm
  REAL, DIMENSION(nmax_domains)  :: dy_wrfnm

  INTEGER :: spec_bdy_width
  INTEGER, DIMENSION(nmax_domains) :: mp_physics
  INTEGER :: io_form

  INTEGER, DIMENSION(nmax_domains) :: idumar
  REAL,    DIMENSION(nmax_domains) :: rdumar

  INTEGER            :: idummy
  REAL               :: rdummy, rdumary(max_vertical_levels)
  LOGICAL            :: ldummy
  CHARACTER(LEN=256) :: cdummy

  INTEGER :: mapproj_moad
  REAL    :: ctrlat_moad
  REAL    :: trulat1_moad, trulat2_moad, trulon_moad

  REAL    :: wrfversion

!-----------------------------------------------------------------------
!
!  WRF grid related variables
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: x_wrf(:)       ! X coordinate of WRF U points
  REAL, ALLOCATABLE :: y_wrf(:)       ! Y coordinate of WRF V points
  REAL, ALLOCATABLE :: xs_wrf(:)      ! X coordinate of WRF mass points
  REAL, ALLOCATABLE :: ys_wrf(:)      ! Y coordinate of WRF mass points
  REAL, ALLOCATABLE :: lat_wrf(:,:,:) ! WRF grid point latitudes
  REAL, ALLOCATABLE :: lon_wrf(:,:,:) ! WRF grid point lontitudes
  REAL, ALLOCATABLE :: msft_wrf(:,:)  ! Map scale factor on mass grid
  REAL, ALLOCATABLE :: msfu_wrf(:,:)  ! Map scale factor on u grid
  REAL, ALLOCATABLE :: msfv_wrf(:,:)  ! Map scale factor on v grid

  REAL, ALLOCATABLE :: cosalpha(:,:)
  REAL, ALLOCATABLE :: sinalpha(:,:)
  REAL, ALLOCATABLE :: f(:,:)
  REAL, ALLOCATABLE :: e(:,:)

  REAL, ALLOCATABLE :: hgt(:,:)
  REAL, ALLOCATABLE :: topostdv(:,:)
  REAL, ALLOCATABLE :: toposlpx(:,:)
  REAL, ALLOCATABLE :: toposlpy(:,:)
  REAL, ALLOCATABLE :: landuse(:,:)
  REAL, ALLOCATABLE :: landmask(:,:)
  REAL, ALLOCATABLE :: landusef(:,:,:)
  REAL, ALLOCATABLE :: soilctop(:,:,:)
  REAL, ALLOCATABLE :: soilcbot(:,:,:)
  REAL, ALLOCATABLE :: green12m(:,:,:)
  REAL, ALLOCATABLE :: albdo12m(:,:,:)
  REAL, ALLOCATABLE :: shdmax(:,:)
  REAL, ALLOCATABLE :: shdmin(:,:)
  REAL, ALLOCATABLE :: tmn(:,:)
  REAL, ALLOCATABLE :: slopecat(:,:)
  REAL, ALLOCATABLE :: snoalb(:,:)

  REAL, DIMENSION(nmax_domains) :: xsub0, ysub0

  REAL      :: lat_ll(4), lat_lr(4), lat_ul(4), lat_ur(4)
  REAL      :: lon_ll(4), lon_lr(4), lon_ul(4), lon_ur(4)

!-----------------------------------------------------------------------
!
! Temporary working arrays
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: tem2d1(:,:), tem2d2(:,:), tem2d3(:,:)
  REAL, ALLOCATABLE :: tem3d(:,:,:)

!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,n
  INTEGER :: istatus, lenstr, ireturn
  INTEGER :: abstime
  INTEGER :: iter

  INTEGER :: domid

  LOGICAL :: fexist,dem_data

  REAL    :: latnot(2),ctrx, ctry, swx, swy

  CHARACTER(LEN=256) :: filename, tmpstr
  CHARACTER(LEN=24)  :: times_str
  INTEGER            :: nchr, ncid, nchout

  REAL,              PARAMETER :: eps = 1.0E-5
  CHARACTER(LEN=16), PARAMETER :: setnames(9) = (/ 'topo_30s        ',  &
           'landuse_30s     ', 'islope          ', 'greenfrac       ',  &
           'maxsnowalb      ', 'soiltemp_1deg   ', 'soiltype_bot_30s',  &
           'soiltype_top_30s', 'albedo_ncep     ' /)

  CHARACTER(LEN=2),  PARAMETER :: setindx(9) = (/ '/U',                 &
           '/V',               '/I',               '/G',                &
           '/M',               '/T',               '/O',                &
           '/O',               '/A' /)
  INTEGER,           PARAMETER :: TOPO = 1, LUSE    = 2, ISLOPE = 3,    &
                                  GFRAC= 4, SNOWALB = 5, SOILTMP= 6,    &
                                  SOILTB=7, SOILTT  = 8, ALBEDO = 9

  CHARACTER(LEN=256) :: path_to_soiltemp_1deg

  COMMON /lapscmn/ path_to_soiltemp_1deg

  REAL :: projrot_latlon
  REAL :: projrot_deg

  REAL :: bilinear_interp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(10(/5x,a),/)')                              &
      '###############################################################',&
      '###############################################################',&
      '####                                                       ####',&
      '####              Welcome to WRFSTATIC                     ####',&
      '####                                                       ####',&
      '#### Program that reads in static datasets and generates   ####',&
      '#### WRF static file for ARPS2WRF or WRF.                  ####',&
      '####                                                       ####',&
      '###############################################################',&
      '###############################################################'

!-----------------------------------------------------------------------
!
! First, read namelist input
!
!-----------------------------------------------------------------------

  CALL readnamelist(1,max_dom,input_from_file,                      &
     hinfmt,adashisfn,adasbasfn,bdybasfn,hisfile,nhisfile,finexist, &
     nxnm,nynm,idumar,idumar,idumar,nprocx_in,nprocy_in,idummy,idummy,  &
     use_arps_grid,nx_wrfnm,ny_wrfnm,idummy,rdumary,rdumar,             &
     i_parent_start,j_parent_start,parent_id,                       &
     mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf,                  &
     ctrlat_wrf,ctrlon_wrf,dx_wrfnm,dy_wrfnm,rdummy,                &
     rdummy,rdummy,rdummy,rdummy,                                   &
     sfcinitopt,idumar,sfcdtfn,geogdir,start_date,silwt,wvln,       &
     idummy,idummy,idummy,idummy,spec_bdy_width,                    &
     idummy,idummy,rdumar,rdumar,mp_physics,idumar,                 &
     idumar,idumar,idumar,idumar,                                   &
     idumar,idummy,idummy,idumar,                                   &
     idummy,rdumar,rdumar,idummy,idummy,idumar,                     &
     idummy,idummy,io_form,idummy,idummy,cdummy,                    &
     idummy,idummy,idumar,cdummy,cdummy,cdummy,idummy,cdummy,       &
     idummy,rdummy,wrfversion,istatus)

  trulat1_wrf = lattru_wrf(1)
  trulat2_wrf = lattru_wrf(2)
  trulon_wrf  = lontru_wrf

  nz_wrf     = 1
  nzsoil_wrf = 1

!  io_form = 7

  nproc_x = 1
  nproc_y = 1
  loc_x   = 1
  loc_y   = 1
  myproc  = 0
  readsplit(:) = 0
  execname = 'WRFSTATIC'

!--------------------------------------------------------------------
!
! Begin to set WRF horizontal grid
!
!--------------------------------------------------------------------

  DO domid = 1, max_dom

    WRITE(6,'(1x,a,I2,3a)') '========= Begin of domain: ',domid,      &
            ', for file - ',TRIM(sfcdtfn(domid)),' ==========='

  WRITE(grdbasfn,'(a)') adasbasfn(domid)

  IF (use_arps_grid(domid) == 1) THEN      ! reset nx_wrf,ny_wrf

    lengbf = LEN_TRIM(grdbasfn)

    WRITE(6,'(1x,/2a/)') ' Getting dimensions from file: ',grdbasfn(1:lengbf)

    CALL get_arps_dims_atts(nchr,hinfmt(domid),grdbasfn,lengbf,         &
                            nx,ny,nz,nzsoil,nstyps,                     &
                            year,month,day,hour,minute,second,time,     &
                            mapproj, sclfct,trulat1,trulat2,trulon,     &
                            ctrlat,ctrlon,dx,dy,ireturn)

    IF( ireturn /= 0 )  CALL arpsstop('get_arps_dims_atts errors.',1)

    nx_wrf  = nx - 2
    ny_wrf  = ny - 2
    dx_wrf  = dx
    dy_wrf  = dy

    mapproj_wrf = mapproj
    sclfct_wrf  = sclfct
    trulat1_wrf = trulat1
    trulat2_wrf = trulat2
    trulon_wrf  = trulon
    ctrlat_wrf  = ctrlat
    ctrlon_wrf  = ctrlon

    nx_wrfnm(domid) = nx_wrf
    ny_wrfnm(domid) = ny_wrf
    dx_wrfnm(domid) = dx_wrf
    dy_wrfnm(domid) = dy_wrf

    IF (domid == 1) THEN
      mapproj_moad = mapproj
      trulat1_moad = trulat1
      trulat2_moad = trulat2
      trulon_moad  = trulon
      ctrlat_moad  = ctrlat_wrf
    ELSE
      IF ( ABS(mapproj_wrf-mapproj_moad) > a_small_real_number .OR.    &
           ABS(trulat1_wrf-trulat1_moad) > a_small_real_number .OR.    &
           ABS(trulat2_wrf-trulat2_moad) > a_small_real_number .OR.    &
           ABS(trulon_wrf -trulon_moad)  > a_small_real_number ) THEN
        WRITE(6,'(1x,a,I2,a,/)')  &
        'Inconsistent map project for subdomian ',domid,' with the mother domain.'

        WRITE(6,'(a38)') '             Mother grid   subdomain grid'
        WRITE(6,'(a38)') '              ==========    =========='
        WRITE(6,'(a,I10,4x,I10)') '  mapproj  = ',mapproj_moad,mapproj_wrf
        WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat1  = ',trulat1,trulat1_wrf
        WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat2  = ',trulat2,trulat2_wrf
        WRITE(6,'(a,F10.2,4x,F10.2)') '  trulon   = ',trulon, trulon_wrf
        CALL arpsstop('Inconsistent map projection',1)
      END IF
    END IF
    WRITE(6,'(a38)') '              ARPS grid      WRF grid'
    WRITE(6,'(a38)') '              ==========    =========='
    WRITE(6,'(a,I10,4x,I10)') '  nx       = ',nx,nx_wrf
    WRITE(6,'(a,I10,4x,I10)') '  ny       = ',ny,ny_wrf
    WRITE(6,'(a,I10,4x,I10)') '  mapproj  = ',mapproj,mapproj_wrf
    WRITE(6,'(a,F10.2,4x,F10.2)') '  sclfct   = ',sclfct,sclfct_wrf
    WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat1  = ',trulat1,trulat1_wrf
    WRITE(6,'(a,F10.2,4x,F10.2)') '  trulat2  = ',trulat2,trulat2_wrf
    WRITE(6,'(a,F10.2,4x,F10.2)') '  trulon   = ',trulon,trulon_wrf
    WRITE(6,'(a,F10.2,4x,F10.2)') '  ctrlat   = ',ctrlat,ctrlat_wrf
    WRITE(6,'(a,F10.2,4x,F10.2)') '  ctrlon   = ',ctrlon,ctrlon_wrf
    WRITE(6,'(a,F10.0,4x,F10.0)') '  dx       = ',dx,dx_wrf
    WRITE(6,'(a,F10.0,4x,F10.0)') '  dy       = ',dy,dy_wrf

    CALL ctim2abss( year,month,day,hour,minute,second, abstime )
    abstime = abstime + INT(time)
    CALL abss2ctim( abstime, year, month, day, hour, minute, second )

    WRITE(start_date(domid),'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2)')  &
          year,'-',month,'-',day,'_',hour,':',minute,':',second

  ELSE        ! use parameters got from namelist above.

    nx_wrf = nx_wrfnm(domid)
    ny_wrf = ny_wrfnm(domid)
    dx_wrf = dx_wrfnm(domid)
    dy_wrf = dy_wrfnm(domid)

    WRITE(6,'(/a/)')      '  The WRF grid settings are:'
    WRITE(6,'(a)')       '               WRF grid'
    WRITE(6,'(a)')       '              =========='
    WRITE(6,'(a,I10)')   '  nx_wrf   = ',nx_wrf
    WRITE(6,'(a,I10)')   '  ny_wrf   = ',ny_wrf
    WRITE(6,'(a,I10)')   '  mapproj  = ',mapproj_wrf
    WRITE(6,'(a,F10.2)') '  sclfct   = ',sclfct_wrf
    WRITE(6,'(a,F10.2)') '  trulat1  = ',trulat1_wrf
    WRITE(6,'(a,F10.2)') '  trulat2  = ',trulat2_wrf
    WRITE(6,'(a,F10.2)') '  trulon   = ',trulon_wrf
    WRITE(6,'(a,F10.2)') '  ctrlat   = ',ctrlat_wrf
    WRITE(6,'(a,F10.2)') '  ctrlon   = ',ctrlon_wrf
    WRITE(6,'(a,F10.0)') '  dx_wrf   = ',dx_wrf
    WRITE(6,'(a,F10.0)') '  dy_wrf   = ',dy_wrf

  END IF

  CALL get_check_grid_ratio(domid,nx_wrf,ny_wrf,dx_wrf,dy_wrf,          &
       dx_wrfnm(parent_id(domid)),dy_wrfnm(parent_id(domid)),           &
       a_small_real_number,parent_grid_ratio,istatus)

  IF (istatus /= 0) CALL arpsstop('Unacceptable parent_grid_ratio or nesting grid size.',1)

!---------------------------------------------------------------------------
!
! Allocate WRF arrays
!
!---------------------------------------------------------------------------

  IF (ALLOCATED(x_wrf)) DEALLOCATE(x_wrf)
  IF (ALLOCATED(y_wrf)) DEALLOCATE(y_wrf)
  IF (ALLOCATED(xs_wrf)) DEALLOCATE(xs_wrf)
  IF (ALLOCATED(ys_wrf)) DEALLOCATE(ys_wrf)
  ALLOCATE(x_wrf (nx_wrf), STAT = istatus)
  ALLOCATE(y_wrf (ny_wrf), STAT = istatus)
  ALLOCATE(xs_wrf(nx_wrf), STAT = istatus)
  ALLOCATE(ys_wrf(ny_wrf), STAT = istatus)

  IF (ALLOCATED(lat_wrf)) DEALLOCATE(lat_wrf)
  IF (ALLOCATED(lon_wrf)) DEALLOCATE(lon_wrf)
  IF (ALLOCATED(msft_wrf)) DEALLOCATE(msft_wrf)
  IF (ALLOCATED(msfu_wrf)) DEALLOCATE(msfu_wrf)
  IF (ALLOCATED(msfv_wrf)) DEALLOCATE(msfv_wrf)
  ALLOCATE(lat_wrf(nx_wrf,ny_wrf,3), STAT = istatus)  ! mass grid
  ALLOCATE(lon_wrf(nx_wrf,ny_wrf,3), STAT = istatus)

  ALLOCATE(msft_wrf(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(msfu_wrf(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(msfv_wrf(nx_wrf,ny_wrf), STAT = istatus)

  IF (ALLOCATED(f)) DEALLOCATE(f)
  IF (ALLOCATED(e)) DEALLOCATE(e)
  IF (ALLOCATED(cosalpha)) DEALLOCATE(cosalpha)
  IF (ALLOCATED(sinalpha)) DEALLOCATE(sinalpha)
  ALLOCATE(f(nx_wrf,ny_wrf),        STAT = istatus)
  ALLOCATE(e(nx_wrf,ny_wrf),        STAT = istatus)
  ALLOCATE(cosalpha(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(sinalpha(nx_wrf,ny_wrf), STAT = istatus)

  IF (ALLOCATED(hgt))      DEALLOCATE(hgt)
  IF (ALLOCATED(topostdv)) DEALLOCATE(topostdv)
  IF (ALLOCATED(toposlpx)) DEALLOCATE(toposlpx)
  IF (ALLOCATED(toposlpy)) DEALLOCATE(toposlpy)
  IF (ALLOCATED(landuse))  DEALLOCATE(landuse)
  IF (ALLOCATED(landmask)) DEALLOCATE(landmask)
  IF (ALLOCATED(slopecat)) DEALLOCATE(slopecat)
  IF (ALLOCATED(tmn))      DEALLOCATE(tmn)
  IF (ALLOCATED(snoalb)) DEALLOCATE(snoalb)
  ALLOCATE(hgt(nx_wrf,ny_wrf),      STAT = istatus)
  ALLOCATE(topostdv(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(toposlpx(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(toposlpy(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(landuse (nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(landmask(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(slopecat(nx_wrf,ny_wrf), STAT = istatus)
  ALLOCATE(tmn(nx_wrf,ny_wrf),      STAT = istatus)
  ALLOCATE(snoalb(nx_wrf,ny_wrf),   STAT = istatus)
  IF (ALLOCATED(shdmax))   DEALLOCATE(shdmax)
  IF (ALLOCATED(shdmin))   DEALLOCATE(shdmin)
  IF (ALLOCATED(landusef)) DEALLOCATE(landusef)
  IF (ALLOCATED(soilctop)) DEALLOCATE(soilctop)
  IF (ALLOCATED(soilcbot)) DEALLOCATE(soilcbot)
  IF (ALLOCATED(green12m)) DEALLOCATE(green12m)
  IF (ALLOCATED(albdo12m)) DEALLOCATE(albdo12m)
  ALLOCATE(shdmax(nx_wrf,ny_wrf),   STAT = istatus)
  ALLOCATE(shdmin(nx_wrf,ny_wrf),   STAT = istatus)
  ALLOCATE(landusef(nx_wrf,ny_wrf,LanduseCategories),  STAT = istatus)
  ALLOCATE(soilctop(nx_wrf,ny_wrf,SoilCategories),     STAT = istatus)
  ALLOCATE(soilcbot(nx_wrf,ny_wrf,SoilCategories),     STAT = istatus)
  ALLOCATE(green12m(nx_wrf,ny_wrf,12),                 STAT = istatus)
  ALLOCATE(albdo12m(nx_wrf,ny_wrf,12),                 STAT = istatus)

  IF (ALLOCATED(tem2d1)) DEALLOCATE(tem2d1)
  IF (ALLOCATED(tem2d2)) DEALLOCATE(tem2d2)
  IF (ALLOCATED(tem2d3)) DEALLOCATE(tem2d3)
  IF (ALLOCATED(tem3d))  DEALLOCATE(tem3d)
  ALLOCATE(tem2d1(nx_wrf,ny_wrf),                   STAT = istatus)
  ALLOCATE(tem2d2(nx_wrf,ny_wrf),                   STAT = istatus)
  ALLOCATE(tem2d3(nx_wrf,ny_wrf),                   STAT = istatus)
  ALLOCATE(tem3d (nx_wrf,ny_wrf,LanduseCategories), STAT = istatus)

  msft_wrf = 0.0
  msfu_wrf = 0.0
  msfv_wrf = 0.0

!
!-----------------------------------------------------------------------
!
!  Establish WRF map projection
!
!-----------------------------------------------------------------------
!
  lattru_wrf(1) = trulat1_wrf
  lattru_wrf(2) = trulat2_wrf
  lontru_wrf    = trulon_wrf
  sclf_wrf      = 1.0/sclfct_wrf
  dx_wrfscl     = dx_wrf*sclf_wrf
  dy_wrfscl     = dy_wrf*sclf_wrf

  IF (domid == 1) THEN
    CALL setmapr(mapproj_wrf,sclf_wrf,lattru_wrf,lontru_wrf)
    CALL lltoxy( 1,1, ctrlat_wrf,ctrlon_wrf, ctrx, ctry )
    swx = ctrx - 0.5*(nx_wrf-1) * dx_wrfscl
    swy = ctry - 0.5*(ny_wrf-1) * dy_wrfscl
    CALL setorig( 1, swx, swy)

    xsub0(domid) = 0
    ysub0(domid) = 0
  ELSE
    IF (use_arps_grid(domid) == 1) THEN
      CALL lltoxy( 1,1, ctrlat_wrf,ctrlon_wrf, ctrx, ctry )
      swx = ctrx - 0.5*(nx_wrf-1) * dx_wrfscl
      swy = ctry - 0.5*(ny_wrf-1) * dy_wrfscl

      xsub0(domid) = swx
      ysub0(domid) = swy
      ctrx = ( xsub0(domid)-xsub0(parent_id(domid)) ) / dx_wrfnm(parent_id(domid)) + 1
      ctry = ( ysub0(domid)-ysub0(parent_id(domid)) ) / dy_wrfnm(parent_id(domid)) + 1
      i_parent_start(domid) = NINT(ctrx)
      j_parent_start(domid) = NINT(ctry)

      IF ( ABS(ctrx - i_parent_start(domid)) > a_small_real_number*i_parent_start(domid) .OR.      &
           ABS(ctry - j_parent_start(domid)) > a_small_real_number*j_parent_start(domid) ) THEN
        WRITE(6,'(1x,a,I2,2(a,F12.2),/,1x,2(a,F12.2),a,I2,/,2(1x,2(a,F12.2),a,/))')      &
        'Cannot align subdomain ',domid,' with central lat/lon = ',ctrlat_wrf,', ',ctrlon_wrf, &
        'and low-left corner = ',swx,', ',swy,', inside the parent grid ',parent_id(domid),   &
        'with grid spacing = ',dx_wrfnm(parent_id(domid)),', ',dy_wrfnm(parent_id(domid)),'.',&
        'The parent grid is starting from ',xsub0(parent_id(domid)),', ',ysub0(parent_id(domid)),'.'
        CALL arpsstop('Inconsistent nesting grid',1)
      END IF

      WRITE(6,'(14x,2(a,i6),/)') 'Starting index are: ',i_parent_start(domid),  &
                                                   ', ',j_parent_start(domid)
    ELSE
      xsub0(domid) = xsub0(parent_id(domid)) + (i_parent_start(domid)-1)*dx_wrfnm(parent_id(domid))
      ysub0(domid) = ysub0(parent_id(domid)) + (j_parent_start(domid)-1)*dy_wrfnm(parent_id(domid))

      ctrx = xsub0(domid) + 0.5*(nx_wrf-1)*dx_wrfscl
      ctry = ysub0(domid) + 0.5*(ny_wrf-1)*dy_wrfscl
      CALL xytoll( 1,1, ctrx, ctry, ctrlat_wrf, ctrlon_wrf)

      WRITE(6,'(14x,2(a,F12.7),/)') 'domain center lat/lon is: ',ctrlat_wrf, ', ', ctrlon_wrf
    END IF
  END IF

  i_parent_end = i_parent_start(domid) + (nx_wrf-1)/parent_grid_ratio
  j_parent_end = j_parent_start(domid) + (ny_wrf-1)/parent_grid_ratio

  IF (i_parent_end > nx_wrfnm(parent_id(domid)) .OR.    &
      j_parent_end > ny_wrfnm(parent_id(domid)) ) THEN
    WRITE(6,'(1x,2(a,I2),a)') 'Subdomain ',domid,' exceed its parent domain ',  &
                                parent_id(domid),' in size.'
    WRITE(6,'(1x,4(a,I8),a)') 'i_parent_start-end, j_parent_start-end = ',      &
                            i_parent_start(domid),' - ',i_parent_end,', ',      &
                            j_parent_start(domid),' - ',j_parent_end,'.'
    WRITE(6,'(1x,a,I2,a,2(I8,a),/)') 'Parent domain ',parent_id(domid),' size: ', &
                            nx_wrfnm(parent_id(domid)),', ',            &
                            ny_wrfnm(parent_id(domid)),'.'

    CALL arpsstop('Subdomain is too large.',1)
  END IF

  DO i=1,nx_wrf
    x_wrf(i)= sclf_wrf*xsub0(domid) + (i-1)*dx_wrfscl
  END DO
  DO j=1,ny_wrf
    y_wrf(j)= sclf_wrf*ysub0(domid) + (j-1)*dy_wrfscl
  END DO

  DO i=1,nx_wrf-1
    xs_wrf(i)=0.5*(x_wrf(i)+x_wrf(i+1))
  END DO
  xs_wrf(nx_wrf)=2.*xs_wrf(nx_wrf-1)-xs_wrf(nx_wrf-2)

  DO j=1,ny_wrf-1
    ys_wrf(j)=0.5*(y_wrf(j)+y_wrf(j+1))
  END DO
  ys_wrf(ny_wrf)=2.*ys_wrf(ny_wrf-1)-ys_wrf(ny_wrf-2)

!-----------------------------------------------------------------------
!
!  Find latitude and longitude of WRF grid.
!
!-----------------------------------------------------------------------

  CALL xytoll(nx_wrf,ny_wrf,xs_wrf,ys_wrf,                          &
              lat_wrf(:,:,1),lon_wrf(:,:,1))    ! T points
  CALL xytoll(nx_wrf,ny_wrf,x_wrf,ys_wrf,                           &
              lat_wrf(:,:,2),lon_wrf(:,:,2))    ! U points
  CALL xytoll(nx_wrf,ny_wrf,xs_wrf,y_wrf,                           &
              lat_wrf(:,:,3),lon_wrf(:,:,3))    ! V points

!-----------------------------------------------------------------------
!
!  Find Map scale factor
!
!-----------------------------------------------------------------------

  CALL lattomf(nx_wrf,ny_wrf,lat_wrf(:,:,1),msft_wrf)  !mass points
  CALL lattomf(nx_wrf,ny_wrf,lat_wrf(:,:,2),msfu_wrf)  !U points
  CALL lattomf(nx_wrf,ny_wrf,lat_wrf(:,:,3),msfv_wrf)  !V points


!-----------------------------------------------------------------------
!
!  Latitude and Longitude variables at corns
!
!-----------------------------------------------------------------------
  !
  ! Latitude, 1 -- T point, 2 -- U point, 3 -- V point, 4 -- massless point
  !
  lat_ll(1) = lat_wrf(       1,       1,1)
  lat_ul(1) = lat_wrf(       1,ny_wrf-1,1)
  lat_ur(1) = lat_wrf(nx_wrf-1,ny_wrf-1,1)
  lat_lr(1) = lat_wrf(nx_wrf-1,       1,1)

  lon_ll(1) = lon_wrf(       1,       1,1)
  lon_ul(1) = lon_wrf(       1,ny_wrf-1,1)
  lon_ur(1) = lon_wrf(nx_wrf-1,ny_wrf-1,1)
  lon_lr(1) = lon_wrf(nx_wrf-1,       1,1)

  lat_ll(2) = lat_wrf(       1,       1,2)
  lat_ul(2) = lat_wrf(       1,ny_wrf-1,2)
  lat_ur(2) = lat_wrf(nx_wrf,  ny_wrf-1,2)
  lat_lr(2) = lat_wrf(nx_wrf,         1,2)

  lon_ll(2) = lon_wrf(       1,       1,2)
  lon_ul(2) = lon_wrf(       1,ny_wrf-1,2)
  lon_ur(2) = lon_wrf(nx_wrf,  ny_wrf-1,2)
  lon_lr(2) = lon_wrf(nx_wrf,         1,2)

  lat_ll(3) = lat_wrf(       1,       1,3)
  lat_ul(3) = lat_wrf(       1,  ny_wrf,3)
  lat_ur(3) = lat_wrf(nx_wrf-1,  ny_wrf,3)
  lat_lr(3) = lat_wrf(nx_wrf-1,       1,3)

  lon_ll(3) = lon_wrf(       1,       1,3)
  lon_ul(3) = lon_wrf(       1,  ny_wrf,3)
  lon_ur(3) = lon_wrf(nx_wrf-1,  ny_wrf,3)
  lon_lr(3) = lon_wrf(nx_wrf-1,       1,3)

  CALL xytoll(1,1,x_wrf(1),     y_wrf(1),     lat_ll(4),lon_ll(4))
  CALL xytoll(1,1,x_wrf(1),     y_wrf(ny_wrf),lat_ul(4),lon_ul(4))
  CALL xytoll(1,1,x_wrf(nx_wrf),y_wrf(ny_wrf),lat_ur(4),lon_ur(4))
  CALL xytoll(1,1,x_wrf(nx_wrf),y_wrf(1),     lat_lr(4),lon_lr(4))

  DO j=1,ny_wrf-1
    DO i=1,nx_wrf-1
      f(i,j)  =2*omega_ear*SIN(lat_wrf(i,j,1)*d2rfactor)
      e(i,j)  =2*omega_ear*COS(lat_wrf(i,j,1)*d2rfactor)
    END DO
  END DO

  DO j=1,ny_wrf-1
    DO i=1,nx_wrf-1
      projrot_deg = projrot_latlon(mapproj_wrf,trulat1_wrf,           &
                           trulat2_wrf,trulon_wrf,ctrlat_wrf,ctrlon_wrf,&
                           lat_wrf(i,j,1),lon_wrf(i,j,1),istatus)
      IF(istatus == 1)THEN
        sinalpha(i,j) = SIN(d2rfactor*projrot_deg)
        cosalpha(i,j) = COS(d2rfactor*projrot_deg)
      ELSE
        CALL wrf_debug(1,'WARNING: Problem in projrot_latlon.')
      END IF
    END DO
  END DO

  times_str = start_date(domid)

  WRITE(6,'(/2a/)') ' Preparing static data set starting at ',times_str


!---------------------- Terrain USGS 30 sec --------------------------
! type = U
! symbol = TOPO

  lenstr = LEN_TRIM(geogdir)

  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(TOPO)),setindx(TOPO)

  CALL rdgeodat(nx_wrf,ny_wrf,x_wrf,y_wrf,dx_wrf,                       &
                dy_wrf,TRIM(tmpstr),wvln,silwt,1,                       &
                tem2d1,tem3d,tem2d2,tem2d3,dem_data,istatus)

  hgt(:,:)=0.0

  DO j = 2,ny_wrf
    DO i = 2,nx_wrf
      hgt(i-1,j-1) = bilinear_interp(i,j,nx_wrf,ny_wrf,tem2d1)
    END DO
  END DO
  WHERE( (hgt(:,:)< 0.01) .AND. (hgt(:,:)> -0.01) ) hgt(:,:) = 0.0

  topostdv(:,:) = tem3d (:,:,1)  ! Actually, they are over non-stagger points, why?
  toposlpx(:,:) = tem2d2(:,:)
  toposlpy(:,:) = tem2d3(:,:)

! ------------------- 30s Landuse -------------------------
! type = V
! symbol = LUSE

  wvln  = 2.0    ! they will be fixed from now on
  silwt = 0.0

  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(LUSE)),setindx(LUSE)

  CALL rdgeodat(nx_wrf,ny_wrf,xs_wrf,ys_wrf,dx_wrf,                     &
                dy_wrf,TRIM(tmpstr),wvln,silwt,LanduseCategories,       &
                tem2d1,tem3d,tem2d2,tem2d3,dem_data,istatus)

  landuse(:,:)  = tem2d1(:,:)
  landmask(:,:) = 1.
  where(landuse(:,:) == ISWATER) landmask(:,:) = 0.

  DO n = 1,LanduseCategories
    landusef(:,:,n) = tem3d(:,:,n)
  END DO

! -------------------------------------------------------------------
! potential fix of fictitous islands for certain resolution domains.
! Story here is that west coast terrain can be "funky" possibly due
! to steepness and method to compute avg terrain.

  tem2d1(:,:) = 1.- tem3d(:,:,ISWATER)

! Filter land fraction with 2dx,iter

  iter = MIN(10,INT(8000./dx_wrf))
  DO n = 1,iter
    CALL filter_2dx(tem2d1,nx_wrf,ny_wrf,1, 0.5)
    CALL filter_2dx(tem2d1,nx_wrf,ny_wrf,1,-0.5)
  END DO

  WHERE(tem2d1(:,:) <= 0.1 .AND. hgt(:,:) < 5.0) hgt(:,:) = 0.0

! -------------- Top Layer Soil Type -------------------------
! type = O
!

  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(SOILTT)),setindx(SOILTT)

  CALL rdgeodat(nx_wrf,ny_wrf,xs_wrf,ys_wrf,dx_wrf,                     &
                dy_wrf,TRIM(tmpstr),wvln,silwt,SoilCategories,          &
                tem2d1,tem3d,tem2d2,tem2d3,dem_data,ireturn)

! make water points = 0 for adjust_geog. later well put it back to original

  WHERE(tem2d1 == 14) tem2d1 = 0.0

  CALL adjust_geog(nx_wrf-1,ny_wrf-1,1,'soiltype',ireturn,              &
         lat_wrf (1:nx_wrf-1,1:ny_wrf-1,1), hgt(1:nx_wrf-1,1:ny_wrf-1), &
         landmask(1:nx_wrf-1,1:ny_wrf-1),tem2d1(1:nx_wrf-1,1:ny_wrf-1), &
         istatus)

  WHERE(tem2d1 == 0.0) tem2d1 = 14.

  DO n = 1,SoilCategories
    soilctop(:,:,n) = tem3d(:,:,n)
  END DO

! -------------- Bottom Layer Soil Type -------------------------
! type = O
!

  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(SOILTB)),setindx(SOILTB)

  CALL rdgeodat(nx_wrf,ny_wrf,xs_wrf,ys_wrf,dx_wrf,                     &
                dy_wrf,TRIM(tmpstr),wvln,silwt,SoilCategories,          &
                tem2d1,tem3d,tem2d2,tem2d3,dem_data,ireturn)

! make water points = 0 for adjust_geog. later well put it back to original

  WHERE(tem2d1 == 14) tem2d1 = 0.0

  CALL adjust_geog(nx_wrf-1,ny_wrf-1,1,'soiltype',ireturn,              &
         lat_wrf (1:nx_wrf-1,1:ny_wrf-1,1), hgt(1:nx_wrf-1,1:ny_wrf-1), &
         landmask(1:nx_wrf-1,1:ny_wrf-1),tem2d1(1:nx_wrf-1,1:ny_wrf-1), &
         istatus)

  WHERE(tem2d1 == 0.0) tem2d1 = 14.

  DO n = 1,SoilCategories
    soilcbot(:,:,n) = tem3d(:,:,n)
  END DO

! ----------------- Greenness Fraction --------------------
! type = G
!

  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(GFRAC)),setindx(GFRAC)

  tem3d(:,:,:) = 0.0

  CALL proc_geodat(nx_wrf,ny_wrf,12,TRIM(tmpstr),                       &
           lat_wrf(:,:,1),lon_wrf(:,:,1),landmask,tem3d,istatus)

  DO n = 1,12
    green12m(:,:,n) = tem3d(:,:,n)
  END DO

! annual max/min greenness fraction in domain
! --------------------------------------------
  print*,' compute max/min greenness frac at grid points'

  DO j=1,ny_wrf
    DO i=1,nx_wrf
      shdmax(i,j) = MAXVAL(green12m(i,j,:))
      shdmin(i,j) = MINVAL(green12m(i,j,:))
     END DO
  END DO

! --------------- Deep Soil Temperature -------------------
! type = T
!
  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(SOILTMP)),setindx(SOILTMP)

  path_to_soiltemp_1deg = tmpstr

  CALL proc_geodat(nx_wrf,ny_wrf,1,TRIM(tmpstr),                        &
          lat_wrf(:,:,1),lon_wrf(:,:,1),landmask,                       &
          tmn,ireturn)

  CALL adjust_geog(nx_wrf-1,ny_wrf-1,1,'soiltemp',ireturn,              &
          lat_wrf(1:nx_wrf-1,1:ny_wrf-1,1),hgt(1:nx_wrf-1,1:ny_wrf-1),  &
          landmask(1:nx_wrf-1,1:ny_wrf-1),tmn(1:nx_wrf-1,1:ny_wrf-1),   &
          istatus)

! -------------- Terrain slope index categories -------------------------
! type = I
!
  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(ISLOPE)),setindx(ISLOPE)

  CALL rdgeodat(nx_wrf,ny_wrf,xs_wrf,ys_wrf,dx_wrf,                     &
                dy_wrf,TRIM(tmpstr),wvln,silwt,9,                       &
                tem2d1,tem3d,tem2d2,tem2d3,dem_data,ireturn)

  slopecat(:,:) = tem2d1(:,:)

! put the categories back to the original raw data. if it is a land
! point but islope indicates water, force islope = 1.
  WHERE(slopecat .eq. 8) slopecat = 13.0
  WHERE(slopecat .eq. 9) slopecat = 0.0

  CALL adjust_geog(nx_wrf-1,ny_wrf-1,1,'islope',ireturn,               &
          lat_wrf(1:nx_wrf-1,1:ny_wrf-1,1),hgt(1:nx_wrf-1,1:ny_wrf-1), &
          landmask(1:nx_wrf-1,1:ny_wrf-1),slopecat(1:nx_wrf-1,1:ny_wrf-1),  &
          istatus)

! force land points to have the correct (default) islope value.
   WHERE(slopecat(:,:) == 0.0 .AND. landmask(:,:) == 1.0)               &
         slopecat(:,:) = 1.0
! force water points to have the correct category for islope
   WHERE(landmask(:,:) == 0.0) slopecat(:,:)=0.0


! ---------------Monthly Albedo Climatology--------------------------
! type = A
!
  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(ALBEDO)),setindx(ALBEDO)

  CALL proc_geodat(nx_wrf,ny_wrf,12,TRIM(tmpstr),                       &
          lat_wrf(:,:,1),lon_wrf(:,:,1),landmask(:,:),tem3d,istatus)


  DO n=1,12
     albdo12m(:,:,n) = tem3d(:,:,n)
     WHERE(landmask(:,:) == 0.0) albdo12m(:,:,n) = 0.08
  END DO

! ---------------- Max Snow Albedo ------------------
! type = M
!

  WRITE(tmpstr,'(3a)') geogdir(1:lenstr),                               &
                       TRIM(setnames(SNOWALB)),setindx(SNOWALB)

  CALL proc_geodat(nx_wrf,ny_wrf,1,TRIM(tmpstr),                        &
          lat_wrf(:,:,1),lon_wrf(:,:,1),landmask(:,:),                  &
          tem2d1,istatus)

  snoalb(:,:) = tem2d1
!
! force max albedo = 0.08 over water. force max albedo = 0.7 over ice
!
  WHERE(landmask(:,:) == 0.0)  snoalb = 0.08
  WHERE(landuse(:,:)  == ISICE)snoalb = 0.7

! ---------------------------------------------------------------------------------
! Lets compare the grids to landmask to assess their consistency (or lack thereof)
! ---------------------------------------------------------------------------------

!  CALL gridcompare(nx_wrf,ny_wrf,i,xxx,landmask)


!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  OUTPUT of WRF static begin
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  filename = sfcdtfn(domid)
  lenstr   = LEN_TRIM(dirname)
  IF(lenstr > 0) filename = dirname(1:lenstr) // TRIM(filename)
  PRINT *, ' Output file name is ',TRIM(filename)

  CALL open_output_file(filename,'STATIC',wrfversion,io_form,           &
                nx_wrf,ny_wrf,nz_wrf, nzsoil_wrf,spec_bdy_width,ncid)

!-----------------------------------------------------------------------
!
!  Initialize and write global attributes
!
!-----------------------------------------------------------------------

  CALL write_static_attribute(ncid,io_form,start_date(domid),domid,     &
       parent_id(domid),i_parent_start(domid),j_parent_start(domid),    &
       i_parent_end,j_parent_end,parent_grid_ratio,                     &
                         nx_wrf,ny_wrf,dx_wrf,dy_wrf,                   &
                         mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,&
                         ctrlat_moad,ctrlat_wrf,ctrlon_wrf,             &
                         lat_ll,lat_ul,lat_ur,lat_lr,                   &
                         lon_ll,lon_ul,lon_ur,lon_lr, istatus)

  CALL write_times_str(ncid,io_form,'Times',                            &
                       times_str(1:19),times_str(1:19),1)

  CALL write_wrf_static(ncid,io_form,times_str(1:19),                   &
                     nx_wrf,ny_wrf,LanduseCategories,SoilCategories,    &
                     lat_wrf(:,:,1),lon_wrf(:,:,1),f,e,                 &
                     msft_wrf,msfu_wrf,msfv_wrf,cosalpha,sinalpha,      &
                     hgt,topostdv,toposlpx,toposlpy,landmask,landusef,  &
                     soilctop,soilcbot,tmn,slopecat,snoalb,             &
                     shdmax,shdmin,green12m,albdo12m, istatus)

  CALL close_output_file(ncid,io_form)    ! close input file

  CALL io_shutdown(io_form)

!-----------------------------------------------------------------------
!
! Generate a ready file if needed
!
!-----------------------------------------------------------------------

  IF ( readyfl == 1 ) THEN

    CALL getunit( nchout )
    WRITE (filename,'(2a)') trim(sfcdtfn(domid)),"_ready"
    WRITE(filename,'(a)') TRIM(dirname) // TRIM(filename)
    OPEN (UNIT=nchout,FILE=trim(filename))
    WRITE (nchout,'(a)') trim(sfcdtfn(domid))
    CLOSE (nchout)
    CALL retunit (nchout )
  END IF

  END DO   ! max_dom

  CALL arpsstop('     ==== WRFSTATIC terminated normally. ====',0)
END PROGRAM wrfstatic
